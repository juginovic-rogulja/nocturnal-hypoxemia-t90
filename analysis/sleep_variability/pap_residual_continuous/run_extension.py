"""
PAP residual-T90 extension, requested by a colleague. Report-only, no manuscript edits.

Base analysis being extended (frozen in numbers/treatment_v2.json, produced by
numbers/run_treatment_v2.py): 1,894 patients with sleep T90 >10% before PAP and a second
recording on PAP. 1,326 reduced residual sleep T90 to <=10% on PAP, 568 stayed >10%.
Cox per outcome: x=not-corrected + age (linear, as originally coded) + male +
log1p(pre-treatment sleep T90) + site strata, unpenalized, event floor 40.

STAGE 1  positive control: rebuild the cohort through the original loader path and refit
         the shipped binary model. Must match numbers/treatment_v2.json exactly or abort.
STAGE 2  (a) post-PAP sleep T90 continuous, within-site rank-inverse-normal per 1 SD
             (the paper's continuous convention for T90 exposures in the primary and
             wake/sleep analyses), BH-FDR across estimable non-control outcomes.
             Sensitivities: log1p scale, age-spline (cr df=4, first column dropped).
         (b) graded categories <=5% (ref), >5-10%, >10%, plus binary cuts >5 and >10.
         (c) crude incidence rates per category for 2 outcomes via an independent
             code path (no lifelines, cohort re-derived from the parquets from scratch).

All inputs are the frozen parquets. Local compute only. Rerunnable end to end.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import warnings
warnings.filterwarnings("ignore")
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from patsy import dmatrix
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT = Path(paths.T90_ROOT)
OUTDIR = Path(f"{paths.SV_ROOT}/"
              "pap_residual_continuous")
OUTDIR.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, paths.ANALYSIS_DIR)
from disease_definitions import DISEASES, NEGATIVE_CONTROLS  # noqa: E402

S = ["wake", "sleep", "nrem", "n1", "n2", "n3", "rem"]
STAGE_PARQUET = ROOT / "data_frozen_v8_2026-09" / "cpap_t90_by_stage.parquet"
BASE_PARQUET = ROOT / "data_frozen_v8_2026-09" / "t90_final.parquet"
FROZEN_JSON = Path(paths.NUMBERS_DIR) / "treatment_v2.json"


# ----------------------------------------------------------------------------- cohort
def build_sev():
    """Verbatim replication of the cohort construction in numbers/run_treatment_v2.py.
    Returns the 1,894-patient severe-baseline cohort with notcorr, male, lbase."""
    d = pd.read_parquet(STAGE_PARQUET)
    A = d[d.design == "A_split"].copy(); A["src"] = "A"
    B = d[d.design == "B_pairs"].copy(); B["src"] = "B"
    P = d[d.design == "P_placebo"]
    for s in S:
        for a, b in [("pre", "dx"), ("post", "tx")]:
            for suf in ["t90", "min"]:
                c = f"{b}_{s}_{suf}"
                B[f"{a}_{s}_{suf}"] = B[c] if c in B.columns else np.nan
    keep = (["BDSPPatientID", "SexDSC", "src"]
            + [f"{p}_{s}_{x}" for p in ("pre", "post") for s in S for x in ("t90", "min")])
    cc = pd.concat([A[[c for c in keep if c in A.columns]],
                    B[[c for c in keep if c in B.columns]]], ignore_index=True)
    cc = cc.sort_values("src", ascending=False).drop_duplicates("BDSPPatientID", keep="first")

    base = pd.read_parquet(BASE_PARQUET)
    m = cc.merge(base, on="BDSPPatientID", how="inner")
    m = m[(m.fu_valid == 1) & (m.post_sleep_min >= 30) & (m.pre_sleep_min >= 30)].copy()
    # original QC block preserved verbatim. The nadir/mean columns are not in the keep
    # list above, so both guards are inert here exactly as in the original script.
    for _side in ("pre", "post"):
        _n, _mn = f"{_side}_sleep_nadir", f"{_side}_sleep_mean"
        if _n in m.columns and _mn in m.columns:
            m = m[~(m[_n] > m[_mn])]
        _t = f"{_side}_sleep_t90"
        if _t in m.columns and _mn in m.columns:
            m = m[~(m[_mn].between(50, 80) & (m[_t] > 70))]

    m["male"] = (m.SexDSC.astype(str).str.upper().str[0] == "M").astype(int)
    m["lbase"] = np.log1p(m.pre_sleep_t90)

    sev = m[m.pre_sleep_t90 > 10].copy()
    sev["notcorr"] = (sev.post_sleep_t90 > 10).astype(int)
    counts = {"designs": {"split_night": int(len(A)), "pairs": int(len(B)),
                          "placebo": int(len(P))},
              "n": int(len(sev)), "n_corrected": int((1 - sev.notcorr).sum()),
              "n_not": int(sev.notcorr.sum())}
    return sev, counts


def outcome_frame(sev, k, xcols):
    """Row set and covariates for one outcome, matching the original construction."""
    f = sev[(sev[f"{k}_prevalent"] == 0) & sev[f"{k}_years"].notna() & (sev[f"{k}_years"] > 0)]
    ff = pd.DataFrame({"T": f[f"{k}_years"], "E": f[f"{k}_incident"].astype(int),
                       **{c: f[c] for c in xcols},
                       "age": f.AgeAtVisit, "male": f.male, "lbase": f.lbase,
                       "site": f.site_id}).dropna()
    return ff


def fit_cox(ff, drop=(), fitter_kwargs=None):
    """Unpenalized Cox, site strata, matching the original. For new fits, retry with
    smaller Newton steps on failure and keep the maximum-log-likelihood solution."""
    cols = [c for c in ff.columns if c not in drop]
    tried = []
    for step in (None, 0.5, 0.25, 0.1):
        try:
            c = CoxPHFitter(penalizer=0.0)
            kw = dict(fitter_kwargs or {})
            if step is not None:
                kw["fit_options"] = {"step_size": step}
            c.fit(ff[cols], "T", "E", strata=["site"], **kw)
            tried.append((float(c.log_likelihood_), step, c))
        except Exception:
            continue
        if step is None:
            break  # default converged, use it (matches the original code path)
    if not tried:
        return None, None
    ll, step, c = max(tried, key=lambda t: t[0])
    return c, step


KEYS = list(DISEASES.keys()) + ["death"]


def label_of(k):
    return DISEASES[k][0] if k in DISEASES else "Death from any cause"


# ----------------------------------------------------------------------------- stage 1
def positive_control(sev, counts):
    frozen = json.load(open(FROZEN_JSON))
    fo = frozen["corrected_vs_not"]
    rows, mismatches = [], []
    if (counts["n"], counts["n_corrected"], counts["n_not"]) != \
       (fo["n"], fo["n_corrected"], fo["n_not"]):
        mismatches.append(f"cohort counts {counts} vs frozen "
                          f"{fo['n']}/{fo['n_corrected']}/{fo['n_not']}")
    if counts["designs"] != frozen["designs"]:
        mismatches.append(f"design counts {counts['designs']} vs {frozen['designs']}")

    for k in KEYS:
        if f"{k}_incident" not in sev.columns:
            continue
        lab = label_of(k)
        ff = outcome_frame(sev, k, ["notcorr"]).rename(columns={"notcorr": "x"})
        if ff.E.sum() < 40:
            if lab in fo["outcomes"]:
                mismatches.append(f"{lab}: below floor here but present in frozen JSON")
            continue
        try:
            c = CoxPHFitter(penalizer=0.0).fit(ff, "T", "E", strata=["site"])
        except Exception as e:
            mismatches.append(f"{lab}: fit failed ({e})")
            continue
        lo, hi = np.exp(c.confidence_intervals_.loc["x"])
        hr = float(np.exp(c.params_["x"]))
        p = float(c.summary.loc["x", "p"])
        got = dict(events=int(ff.E.sum()), hr=round(hr, 3), lo=round(float(lo), 3),
                   hi=round(float(hi), 3), p=p)
        ref = fo["outcomes"].get(lab)
        ok = (ref is not None and got["events"] == ref["events"]
              and abs(got["hr"] - ref["hr"]) <= 5e-4 and abs(got["lo"] - ref["lo"]) <= 5e-4
              and abs(got["hi"] - ref["hi"]) <= 5e-4
              and abs(got["p"] - ref["p"]) <= 1e-9 * max(ref["p"], 1e-300))
        if not ok:
            mismatches.append(f"{lab}: got {got} vs frozen {ref}")
        rows.append(dict(outcome=lab, match=ok,
                         events_new=got["events"], hr_new=got["hr"],
                         lo_new=got["lo"], hi_new=got["hi"], p_new=got["p"],
                         events_frozen=ref["events"] if ref else np.nan,
                         hr_frozen=ref["hr"] if ref else np.nan,
                         lo_frozen=ref["lo"] if ref else np.nan,
                         hi_frozen=ref["hi"] if ref else np.nan,
                         p_frozen=ref["p"] if ref else np.nan))
    fitted_labels = {r["outcome"] for r in rows}
    for lab in fo["outcomes"]:
        if lab not in fitted_labels:
            mismatches.append(f"{lab}: in frozen JSON but not fitted here")
    pc = pd.DataFrame(rows)
    pc.to_csv(OUTDIR / "positive_control.csv", index=False)
    return pc, mismatches, fo


# ----------------------------------------------------------------------------- stage 2
def rint(x):
    r = stats.rankdata(x)
    return stats.norm.ppf((r - 0.375) / (len(r) + 0.25))


def run_extensions(sev, fo):
    # exposure codings on the full 1,894 before any outcome subsetting, so every model
    # sees the same coding (mirrors how the paper codes exposures on the analysis set)
    sev = sev.copy()
    sev["z_post"] = sev.groupby("site_id").post_sleep_t90.transform(rint)
    sev["l_post"] = np.log1p(sev.post_sleep_t90)
    sev["cat_mid"] = ((sev.post_sleep_t90 > 5) & (sev.post_sleep_t90 <= 10)).astype(int)
    sev["cat_high"] = (sev.post_sleep_t90 > 10).astype(int)
    sev["over5"] = (sev.post_sleep_t90 > 5).astype(int)

    # age spline on the analysis cohort, first column dropped (primary-analysis rule)
    sp = dmatrix("cr(a, df=4) - 1", {"a": sev.AgeAtVisit.values},
                 return_type="dataframe").iloc[:, 1:]
    for i in range(sp.shape[1]):
        sev[f"age_s{i}"] = sp.iloc[:, i].values
    SPL = [f"age_s{i}" for i in range(sp.shape[1])]

    principal = [lab for lab, v in fo["outcomes"].items()
                 if v["p"] < .05 and not v["negative_control"]]

    rows = []
    for k in KEYS:
        if f"{k}_incident" not in sev.columns:
            continue
        lab = label_of(k)
        neg = k in NEGATIVE_CONTROLS
        ff = outcome_frame(sev, k, ["z_post", "l_post", "cat_mid", "cat_high", "over5",
                                    "notcorr"] + SPL)
        if ff.E.sum() < 40:
            continue
        rec = dict(key=k, outcome=lab, negative_control=neg,
                   n=int(len(ff)), events=int(ff.E.sum()),
                   events_le5=int(ff.loc[(ff.cat_mid == 0) & (ff.cat_high == 0), "E"].sum()),
                   events_5to10=int(ff.loc[ff.cat_mid == 1, "E"].sum()),
                   events_gt10=int(ff.loc[ff.cat_high == 1, "E"].sum()),
                   principal=lab in principal)

        # (a) continuous, RIN per 1 SD, linear age exactly as the original model
        c, step = fit_cox(ff, drop=["l_post", "cat_mid", "cat_high", "over5",
                                    "notcorr"] + SPL)
        if c is not None:
            lo, hi = np.exp(c.confidence_intervals_.loc["z_post"])
            rec.update(hr_sd=float(np.exp(c.params_["z_post"])), lo_sd=float(lo),
                       hi_sd=float(hi), p_sd=float(c.summary.loc["z_post", "p"]),
                       step_sd="default" if step is None else step)

        # sensitivity: log1p scale
        c, step = fit_cox(ff, drop=["z_post", "cat_mid", "cat_high", "over5",
                                    "notcorr"] + SPL)
        if c is not None:
            lo, hi = np.exp(c.confidence_intervals_.loc["l_post"])
            rec.update(hr_log1p=float(np.exp(c.params_["l_post"])), lo_log1p=float(lo),
                       hi_log1p=float(hi), p_log1p=float(c.summary.loc["l_post", "p"]))

        # sensitivity: RIN exposure with the age spline instead of linear age
        c, step = fit_cox(ff, drop=["l_post", "cat_mid", "cat_high", "over5",
                                    "notcorr", "age"])
        if c is not None:
            lo, hi = np.exp(c.confidence_intervals_.loc["z_post"])
            rec.update(hr_sd_spline=float(np.exp(c.params_["z_post"])),
                       lo_sd_spline=float(lo), hi_sd_spline=float(hi),
                       p_sd_spline=float(c.summary.loc["z_post", "p"]))

        # (b) graded categories, ref <=5%
        c, step = fit_cox(ff, drop=["z_post", "l_post", "over5", "notcorr"] + SPL)
        if c is not None:
            for term, tag in (("cat_mid", "5to10"), ("cat_high", "gt10")):
                lo, hi = np.exp(c.confidence_intervals_.loc[term])
                rec[f"hr_{tag}"] = float(np.exp(c.params_[term]))
                rec[f"lo_{tag}"] = float(lo)
                rec[f"hi_{tag}"] = float(hi)
                rec[f"p_{tag}"] = float(c.summary.loc[term, "p"])

        # (b) simple binary cuts
        c, step = fit_cox(ff, drop=["z_post", "l_post", "cat_mid", "cat_high",
                                    "notcorr"] + SPL)
        if c is not None:
            lo, hi = np.exp(c.confidence_intervals_.loc["over5"])
            rec.update(hr_bin5=float(np.exp(c.params_["over5"])), lo_bin5=float(lo),
                       hi_bin5=float(hi), p_bin5=float(c.summary.loc["over5", "p"]))
        c, step = fit_cox(ff, drop=["z_post", "l_post", "cat_mid", "cat_high",
                                    "over5"] + SPL)
        if c is not None:
            lo, hi = np.exp(c.confidence_intervals_.loc["notcorr"])
            rec.update(hr_bin10=float(np.exp(c.params_["notcorr"])), lo_bin10=float(lo),
                       hi_bin10=float(hi), p_bin10=float(c.summary.loc["notcorr", "p"]))
        rows.append(rec)

    res = pd.DataFrame(rows)

    # BH-FDR across estimable non-control outcomes (death included), matching the
    # paper's significance denominator for the treatment family (44 non-control
    # outcomes). Negative controls are a falsification check, not discoveries, so
    # they are outside the family and get no q.
    fam = res[~res.negative_control & res.p_sd.notna()].copy()
    q = multipletests(fam.p_sd.values, method="fdr_bh")[1]
    res["q_sd"] = np.nan
    res.loc[fam.index, "q_sd"] = q
    res["sig_fdr"] = res.q_sd < 0.05

    cont_cols = ["key", "outcome", "negative_control", "principal", "n", "events",
                 "hr_sd", "lo_sd", "hi_sd", "p_sd", "q_sd", "sig_fdr", "step_sd",
                 "hr_log1p", "lo_log1p", "hi_log1p", "p_log1p",
                 "hr_sd_spline", "lo_sd_spline", "hi_sd_spline", "p_sd_spline"]
    grad_cols = ["key", "outcome", "negative_control", "principal", "n", "events",
                 "events_le5", "events_5to10", "events_gt10",
                 "hr_5to10", "lo_5to10", "hi_5to10", "p_5to10",
                 "hr_gt10", "lo_gt10", "hi_gt10", "p_gt10",
                 "hr_bin5", "lo_bin5", "hi_bin5", "p_bin5",
                 "hr_bin10", "lo_bin10", "hi_bin10", "p_bin10"]
    res[cont_cols].to_csv(OUTDIR / "results_continuous.csv", index=False)
    res[grad_cols].to_csv(OUTDIR / "results_graded.csv", index=False)
    return res, sev


# --------------------------------------------------------------- stage 2c, independent
def crude_check(outcome_keys=("hf", "cvd")):
    """Crude incidence per post-PAP category, written as an independent code path:
    cohort re-derived from the parquets with different pandas operations, rates by
    plain numpy arithmetic, no lifelines anywhere."""
    d = pd.read_parquet(STAGE_PARQUET)
    picks = {}
    for _, r in d.iterrows():
        if r["design"] == "A_split":
            pre, post, mnp, mnq = (r.get("pre_sleep_t90"), r.get("post_sleep_t90"),
                                   r.get("pre_sleep_min"), r.get("post_sleep_min"))
        elif r["design"] == "B_pairs":
            pre, post, mnp, mnq = (r.get("dx_sleep_t90"), r.get("tx_sleep_t90"),
                                   r.get("dx_sleep_min"), r.get("tx_sleep_min"))
        else:
            continue
        pid = r["BDSPPatientID"]
        # B_pairs beats A_split when a patient has both, as in the original
        if pid in picks and picks[pid]["design"] == "B_pairs" and r["design"] == "A_split":
            continue
        picks[pid] = dict(design=r["design"], pre=pre, post=post, mnp=mnp, mnq=mnq)
    cw = pd.DataFrame.from_dict(picks, orient="index")
    cw.index.name = "BDSPPatientID"
    cw = cw.reset_index()

    base = pd.read_parquet(BASE_PARQUET)
    j = cw.merge(base, on="BDSPPatientID", how="inner")
    j = j[(j.fu_valid == 1) & (j.mnp >= 30) & (j.mnq >= 30) & (j.pre > 10)]
    cat = np.where(j.post > 10, ">10%", np.where(j.post > 5, ">5-10%", "<=5%"))
    j = j.assign(cat=cat)

    out = []
    for k in outcome_keys:
        lab = label_of(k)
        f = j[(j[f"{k}_prevalent"] == 0) & j[f"{k}_years"].notna() & (j[f"{k}_years"] > 0)]
        f = f[np.isfinite(f[f"{k}_incident"].astype(float))]
        for cval in ["<=5%", ">5-10%", ">10%"]:
            g = f[f.cat == cval]
            ev = float(g[f"{k}_incident"].sum())
            py = float(g[f"{k}_years"].sum())
            out.append(dict(outcome=lab, category=cval, n=int(len(g)), events=int(ev),
                            person_years=round(py, 1),
                            rate_per_1000py=round(1000 * ev / py, 2) if py > 0 else np.nan))
    cr = pd.DataFrame(out)
    cr.to_csv(OUTDIR / "crude_incidence.csv", index=False)
    return cr


# ------------------------------------------------------------------------------ report
def fmt_hr(hr, lo, hi):
    return f"{hr:.2f} ({lo:.2f}-{hi:.2f})"


def write_report(counts, pc, mismatches, res, cr, sev):
    L = []
    n_est = int(res.p_sd.notna().sum())
    n_fam = int((~res.negative_control & res.p_sd.notna()).sum())
    sig = res[res.sig_fdr.fillna(False)].sort_values("q_sd")
    site_sizes = sev.groupby("site_id").size().to_dict()
    cat_n = {"<=5%": int(((sev.post_sleep_t90 <= 5)).sum()),
             ">5-10%": int(((sev.post_sleep_t90 > 5) & (sev.post_sleep_t90 <= 10)).sum()),
             ">10%": int((sev.post_sleep_t90 > 10).sum())}

    L.append("PAP RESIDUAL T90 AS A CONTINUOUS AND GRADED EXPOSURE")
    L.append("Extension of the T90 manuscript treatment analysis. Report only, nothing")
    L.append("in the manuscript, supplement or LATEST folders was touched.")
    L.append("")
    L.append("COHORT (identical to the shipped binary analysis)")
    L.append(f"  {counts['n']:,} patients, pre-PAP sleep T90 >10%, second recording on PAP.")
    L.append(f"  {counts['n_corrected']:,} corrected to <=10% on PAP, {counts['n_not']:,} stayed >10%.")
    L.append(f"  Post-PAP categories: <=5% n={cat_n['<=5%']:,}, >5-10% n={cat_n['>5-10%']:,}, "
             f">10% n={cat_n['>10%']:,}.")
    L.append(f"  Sites in the model strata: {site_sizes}.")
    L.append("")

    L.append("A. CONTINUOUS POST-PAP SLEEP T90 (per 1 SD, rank-inverse-normal within site)")
    L.append(f"  Outcomes estimable at the 40-event floor: {n_est} "
             f"({n_fam} non-control incl death, {n_est - n_fam} negative controls).")
    L.append(f"  Significant after BH-FDR (q<0.05, family = {n_fam} non-control outcomes): "
             f"{len(sig)}.")
    L.append("")
    L.append(f"  {'outcome':<28}{'events':>7}   {'HR/SD (95% CI)':<20}{'p':>10}{'q':>10}")
    for _, r in sig.iterrows():
        L.append(f"  {r.outcome:<28}{r.events:>7}   {fmt_hr(r.hr_sd, r.lo_sd, r.hi_sd):<20}"
                 f"{r.p_sd:>10.2e}{r.q_sd:>10.4f}")
    L.append("")
    ncsub = res[res.negative_control & res.p_sd.notna()]
    if len(ncsub):
        L.append(f"  Negative controls (outside the FDR family): HR/SD "
                 + ", ".join(f"{r.outcome} {r.hr_sd:.2f} (p={r.p_sd:.2f})"
                             for _, r in ncsub.iterrows()))
    L.append("")

    L.append("B. GRADED POST-PAP T90: <=5% (reference), >5-10%, >10%")
    L.append("  Principal outcomes = the 16 conditions plus the cardiovascular composite that")
    L.append("  were significant in the shipped binary analysis.")
    L.append("")
    L.append(f"  {'outcome':<28}{'ev 5/10/>10':>13}   {'HR >5-10%':<20}{'HR >10%':<20}")
    pr = res[res.principal].sort_values("p_sd")
    for _, r in pr.iterrows():
        evs = f"{r.events_le5}/{r.events_5to10}/{r.events_gt10}"
        m1 = fmt_hr(r.hr_5to10, r.lo_5to10, r.hi_5to10) if np.isfinite(r.get("hr_5to10", np.nan)) else "-"
        m2 = fmt_hr(r.hr_gt10, r.lo_gt10, r.hi_gt10) if np.isfinite(r.get("hr_gt10", np.nan)) else "-"
        L.append(f"  {r.outcome:<28}{evs:>13}   {m1:<20}{m2:<20}")
    med_mid = float(pr.hr_5to10.median())
    med_high = float(pr.hr_gt10.median())
    n_mid_sig = int((pr.p_5to10 < .05).sum())
    n_high_sig = int((pr.p_gt10 < .05).sum())
    L.append("")
    L.append(f"  Across the {len(pr)} principal outcomes: median HR {med_mid:.2f} for >5-10% "
             f"({n_mid_sig} of {len(pr)} p<0.05) and {med_high:.2f} for >10% "
             f"({n_high_sig} of {len(pr)} p<0.05).")
    L.append("")
    L.append("  Binary cuts on the same model, principal outcomes:")
    L.append(f"  {'outcome':<28}{'>5 vs <=5':<22}{'>10 vs <=10':<22}")
    for _, r in pr.iterrows():
        b5 = fmt_hr(r.hr_bin5, r.lo_bin5, r.hi_bin5) if np.isfinite(r.get("hr_bin5", np.nan)) else "-"
        b10 = fmt_hr(r.hr_bin10, r.lo_bin10, r.hi_bin10) if np.isfinite(r.get("hr_bin10", np.nan)) else "-"
        L.append(f"  {r.outcome:<28}{b5:<22}{b10:<22}")
    L.append("")

    L.append("C. SANITY CHECK, crude incidence via an independent code path (no lifelines,")
    L.append("   cohort re-derived from the parquets from scratch)")
    for lab in cr.outcome.unique():
        sub = cr[cr.outcome == lab]
        parts = ", ".join(f"{r.category} {r.rate_per_1000py}/1000py "
                          f"({r.events} ev, n={r.n})" for _, r in sub.iterrows())
        L.append(f"  {lab}: {parts}")
    L.append("")

    L.append("POSITIVE CONTROL")
    if mismatches:
        L.append("  FAILED, see mismatches below. Extensions were not run.")
        for msg in mismatches:
            L.append(f"    {msg}")
    else:
        L.append(f"  Rebuilt through the original loader path: designs "
                 f"{counts['designs']}, cohort {counts['n']:,} = "
                 f"{counts['n_corrected']:,} + {counts['n_not']:,}.")
        L.append(f"  All {len(pc)} shipped binary outcomes refit and matched "
                 f"numbers/treatment_v2.json exactly (events, HR, CI at 3 dp, p).")
    L.append("")

    L.append("SPECIFICATION, as in the original treatment model unless stated")
    L.append("  Loader: data_frozen_v8_2026-09/cpap_t90_by_stage.parquet (A_split + B_pairs, B wins")
    L.append("  duplicates) merged into data_frozen_v8_2026-09/t90_final.parquet, fu_valid=1,")
    L.append("  pre and post sleep >=30 min, pre-PAP sleep T90 >10%.")
    L.append("  Covariates: log1p(pre-PAP sleep T90) exactly as originally coded, age as a")
    L.append("  LINEAR term because that is the original treatment-model convention (the")
    L.append("  age spline with one dropped column is the primary-analysis convention and")
    L.append("  is included as a sensitivity), male, site as Cox strata. Unpenalized fits,")
    L.append("  event floor 40, lifelines CoxPHFitter.")
    L.append("  Continuous exposure: rank-inverse-normal within site, per 1 SD, matching the")
    L.append("  paper's continuous T90 convention in the primary and wake/sleep analyses.")
    L.append("  log1p(post T90) sensitivity in results_continuous.csv.")
    L.append("  FDR: BH across the estimable non-control outcomes, one family, death")
    L.append("  included. The shipped binary analysis used raw p<0.05 with negative")
    L.append("  controls as the falsification check, so its counts are not FDR-based.")
    L.append("")
    L.append("PATHS")
    L.append(f"  {OUTDIR}/run_extension.py")
    L.append(f"  {OUTDIR}/results_continuous.csv")
    L.append(f"  {OUTDIR}/results_graded.csv")
    L.append(f"  {OUTDIR}/crude_incidence.csv")
    L.append(f"  {OUTDIR}/positive_control.csv")
    L.append(f"  {OUTDIR}/REPORT.txt")

    (OUTDIR / "REPORT.txt").write_text("\n".join(L) + "\n")
    return "\n".join(L)


def main():
    sev, counts = build_sev()
    pc, mismatches, fo = positive_control(sev, counts)
    print(f"positive control: {int(pc.match.sum())}/{len(pc)} outcomes match, "
          f"{len(mismatches)} mismatches")
    if mismatches:
        for msg in mismatches:
            print("  MISMATCH:", msg)
        (OUTDIR / "REPORT.txt").write_text(
            "POSITIVE CONTROL FAILED. Extensions were not run.\n" + "\n".join(mismatches) + "\n")
        sys.exit(1)

    res, sev2 = run_extensions(sev, fo)
    cr = crude_check()
    report = write_report(counts, pc, mismatches, res, cr, sev2)
    print(report)


if __name__ == "__main__":
    main()
