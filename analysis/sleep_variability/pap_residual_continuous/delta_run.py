"""
Follow-up extension: the CORRECTED amount as the exposure.

Same 1,894-patient PAP cohort. Model per outcome: delta = pre-PAP sleep T90 minus
post-PAP sleep T90 as the exposure, ADJUSTED for post-PAP sleep T90 (the residual),
plus age (linear, as the original treatment model), male, site strata. Unpenalized,
40-event floor, same outcome set, BH-FDR across the non-control family.

Codings
  primary      z_delta + z_post, both within-site rank-inverse-normal per 1 SD
  sensitivity  raw per-10-percentage-points for both delta and post
  sensitivity  z_delta with log1p(post) adjustment

Interpretive check, done numerically, not asserted:
  pre = post + delta identically, so {post, delta} spans the same column space as
  {post, pre} in raw coding. The raw delta coefficient given post must EQUAL the raw
  pre coefficient given post (exact reparametrization). Under rank transforms the
  identity is approximate. Both are quantified here: the coefficient identity across
  all outcomes, and the partial correlation of delta with pre given post.

Gate: the shipped binary analysis must reproduce exactly (positive control from
run_extension.py) before anything else runs. Report-only, no manuscript edits.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

OUTDIR = Path(f"{paths.SV_ROOT}/"
              "pap_residual_continuous")
sys.path.insert(0, str(OUTDIR))
from run_extension import (build_sev, positive_control, outcome_frame, fit_cox,  # noqa: E402
                           rint, KEYS, label_of, STAGE_PARQUET, BASE_PARQUET)
from disease_definitions import NEGATIVE_CONTROLS  # noqa: E402


def term(c, name):
    lo, hi = np.exp(c.confidence_intervals_.loc[name])
    return dict(hr=float(np.exp(c.params_[name])), lo=float(lo), hi=float(hi),
                p=float(c.summary.loc[name, "p"]), beta=float(c.params_[name]))


def run_delta(sev, fo):
    sev = sev.copy()
    sev["delta"] = sev.pre_sleep_t90 - sev.post_sleep_t90
    sev["z_delta"] = sev.groupby("site_id").delta.transform(rint)
    sev["z_post"] = sev.groupby("site_id").post_sleep_t90.transform(rint)
    sev["z_pre"] = sev.groupby("site_id").pre_sleep_t90.transform(rint)
    sev["raw10_delta"] = sev.delta / 10.0
    sev["raw10_post"] = sev.post_sleep_t90 / 10.0
    sev["raw10_pre"] = sev.pre_sleep_t90 / 10.0
    sev["l_post"] = np.log1p(sev.post_sleep_t90)

    XC = ["z_delta", "z_post", "z_pre", "raw10_delta", "raw10_post", "raw10_pre",
          "l_post"]
    principal = [lab for lab, v in fo["outcomes"].items()
                 if v["p"] < .05 and not v["negative_control"]]

    rows, ident = [], []
    for k in KEYS:
        if f"{k}_incident" not in sev.columns:
            continue
        lab = label_of(k)
        ff = outcome_frame(sev, k, XC)
        if ff.E.sum() < 40:
            continue
        rec = dict(key=k, outcome=lab, negative_control=k in NEGATIVE_CONTROLS,
                   principal=lab in principal, n=int(len(ff)), events=int(ff.E.sum()))

        # M1 primary: z_delta exposure, z_post adjustment. lbase is excluded because
        # pre-treatment severity enters through the {post, delta} span by construction.
        c, _ = fit_cox(ff, drop=["z_pre", "raw10_delta", "raw10_post", "raw10_pre",
                                 "l_post", "lbase"])
        if c is not None:
            t = term(c, "z_delta")
            rec.update(hr_delta_sd=t["hr"], lo_delta_sd=t["lo"], hi_delta_sd=t["hi"],
                       p_delta_sd=t["p"], beta_delta_sd=t["beta"])
            tp = term(c, "z_post")
            rec.update(hr_post_sd_inM1=tp["hr"], p_post_sd_inM1=tp["p"])

        # M4 comparison: z_pre exposure, z_post adjustment (same rows, same covariates)
        c, _ = fit_cox(ff, drop=["z_delta", "raw10_delta", "raw10_post", "raw10_pre",
                                 "l_post", "lbase"])
        if c is not None:
            t = term(c, "z_pre")
            rec.update(hr_pre_sd=t["hr"], lo_pre_sd=t["lo"], hi_pre_sd=t["hi"],
                       p_pre_sd=t["p"], beta_pre_sd=t["beta"])

        # M2 sensitivity: raw per 10 percentage points, delta + post
        c, _ = fit_cox(ff, drop=["z_delta", "z_post", "z_pre", "raw10_pre",
                                 "l_post", "lbase"])
        b_delta_raw = np.nan
        if c is not None:
            t = term(c, "raw10_delta")
            rec.update(hr_delta_raw10=t["hr"], lo_delta_raw10=t["lo"],
                       hi_delta_raw10=t["hi"], p_delta_raw10=t["p"])
            b_delta_raw = t["beta"]

        # M5 identity partner: raw pre + raw post. Exact reparametrization of M2, so
        # beta(raw10_pre | post) must equal beta(raw10_delta | post).
        c, _ = fit_cox(ff, drop=["z_delta", "z_post", "z_pre", "raw10_delta",
                                 "l_post", "lbase"])
        if c is not None and np.isfinite(b_delta_raw):
            b_pre_raw = term(c, "raw10_pre")["beta"]
            ident.append(dict(outcome=lab, beta_delta_raw=b_delta_raw,
                              beta_pre_raw=b_pre_raw,
                              absdiff=abs(b_delta_raw - b_pre_raw)))

        # M3 sensitivity: z_delta with log1p(post) adjustment
        c, _ = fit_cox(ff, drop=["z_post", "z_pre", "raw10_delta", "raw10_post",
                                 "raw10_pre", "lbase"])
        if c is not None:
            t = term(c, "z_delta")
            rec.update(hr_delta_lpost=t["hr"], lo_delta_lpost=t["lo"],
                       hi_delta_lpost=t["hi"], p_delta_lpost=t["p"])
        rows.append(rec)

    res = pd.DataFrame(rows)
    fam = res[~res.negative_control & res.p_delta_sd.notna()].copy()
    res["q_delta_sd"] = np.nan
    res.loc[fam.index, "q_delta_sd"] = multipletests(fam.p_delta_sd.values,
                                                     method="fdr_bh")[1]
    res["sig_fdr"] = res.q_delta_sd < 0.05
    res.to_csv(OUTDIR / "delta_results.csv", index=False)
    ident = pd.DataFrame(ident)
    ident.to_csv(OUTDIR / "delta_identity_check.csv", index=False)
    return res, ident, sev


def collinearity(sev):
    """How much of delta is pre, once post is held fixed."""
    out = {}
    # raw scale: pre = post + delta identically, so given post the map delta -> pre
    # is deterministic. Verified, not assumed:
    out["raw_identity_max_abs_err"] = float(
        (sev.pre_sleep_t90 - (sev.post_sleep_t90 + sev.delta)).abs().max())

    def partial(a, b, ctrl):
        A = np.column_stack([np.ones(len(ctrl)), ctrl])
        ra = a - A @ np.linalg.lstsq(A, a, rcond=None)[0]
        rb = b - A @ np.linalg.lstsq(A, b, rcond=None)[0]
        return float(np.corrcoef(ra, rb)[0, 1])

    z = sev[["z_delta", "z_pre", "z_post"]].dropna()
    out["rin_corr_delta_pre"] = float(np.corrcoef(z.z_delta, z.z_pre)[0, 1])
    out["rin_partial_delta_pre_given_post"] = partial(
        z.z_delta.values, z.z_pre.values, z.z_post.values)
    out["n_delta_negative"] = int((sev.delta < 0).sum())
    out["delta_median_iqr"] = [float(sev.delta.median()),
                               float(sev.delta.quantile(.25)),
                               float(sev.delta.quantile(.75))]
    return out


def crude_delta(outcome_keys=("hf", "cvd")):
    """Independent code path: cohort re-derived from the parquets from scratch,
    plain numpy rates. Fixed post stratum = post-PAP sleep T90 <=5% (the largest and
    tightest residual band), delta tertiles within it."""
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
        if pid in picks and picks[pid]["design"] == "B_pairs" and r["design"] == "A_split":
            continue
        picks[pid] = dict(design=r["design"], pre=pre, post=post, mnp=mnp, mnq=mnq)
    cw = pd.DataFrame.from_dict(picks, orient="index")
    cw.index.name = "BDSPPatientID"
    cw = cw.reset_index()
    base = pd.read_parquet(BASE_PARQUET)
    j = cw.merge(base, on="BDSPPatientID", how="inner")
    j = j[(j.fu_valid == 1) & (j.mnp >= 30) & (j.mnq >= 30) & (j.pre > 10)]
    j = j[j.post <= 5].copy()
    j["delta"] = j.pre - j.post
    cuts = j.delta.quantile([1 / 3, 2 / 3]).values
    j["tert"] = np.where(j.delta <= cuts[0], "T1 low",
                         np.where(j.delta <= cuts[1], "T2 mid", "T3 high"))
    out = []
    for k in outcome_keys:
        lab = label_of(k)
        f = j[(j[f"{k}_prevalent"] == 0) & j[f"{k}_years"].notna() & (j[f"{k}_years"] > 0)]
        for tv in ["T1 low", "T2 mid", "T3 high"]:
            g = f[f.tert == tv]
            ev = float(g[f"{k}_incident"].sum())
            py = float(g[f"{k}_years"].sum())
            out.append(dict(outcome=lab, stratum="post<=5%", delta_tertile=tv,
                            n=int(len(g)), events=int(ev),
                            median_delta=round(float(g.delta.median()), 1),
                            median_pre=round(float(g.pre.median()), 1),
                            median_post=round(float(g.post.median()), 2),
                            person_years=round(py, 1),
                            rate_per_1000py=round(1000 * ev / py, 2) if py > 0 else np.nan))
    cr = pd.DataFrame(out)
    cr.to_csv(OUTDIR / "delta_crude.csv", index=False)
    return cr


def fmt(hr, lo, hi):
    return f"{hr:.2f} ({lo:.2f}-{hi:.2f})"


def write_report(counts, res, ident, col, cr, pc_n):
    L = []
    fam = res[~res.negative_control & res.p_delta_sd.notna()]
    sig = res[res.sig_fdr.fillna(False)].sort_values("q_delta_sd")
    n_up = int((fam.hr_delta_sd > 1).sum())
    n_dn = int((fam.hr_delta_sd < 1).sum())

    L.append("DELTA (AMOUNT CORRECTED BY PAP) AS THE EXPOSURE, ADJUSTED FOR RESIDUAL T90")
    L.append("Same 1,894-patient cohort, age + sex + site strata as the original")
    L.append("treatment model, post-PAP sleep T90 in the model, delta = pre minus post")
    L.append("as the exposure. Report only, no manuscript edits.")
    L.append("")
    L.append("READ THIS FIRST, THE MODEL IDENTITY")
    L.append("  pre = post + delta holds row by row (max abs error "
             f"{col['raw_identity_max_abs_err']:.1e}), so in raw coding the model")
    L.append("  post + delta is an exact reparametrization of post + pre. Verified on")
    L.append(f"  the fits: across {len(ident)} outcomes the raw delta coefficient given")
    L.append("  post equals the raw pre coefficient given post to max abs diff "
             f"{ident.absdiff.max():.2e}.")
    L.append("  Under the rank transform the identity is near-exact in this cohort:")
    L.append(f"  corr(z_delta, z_pre) = {col['rin_corr_delta_pre']:.3f}, partial corr")
    L.append(f"  given z_post = {col['rin_partial_delta_pre_given_post']:.3f}.")
    L.append("  So \"larger PAP-related improvement given the same residual\" and \"worse")
    L.append("  pretreatment severity given the same residual\" are the SAME statistical")
    L.append("  quantity here. The model can answer whether risk tracks the corrected")
    L.append("  amount once the residual is fixed. It cannot say whether that reflects")
    L.append("  a benefit of correction or a scar of starting severity, because the two")
    L.append("  are one variable in this design.")
    L.append("")

    L.append("A. PRIMARY: delta per 1 SD (within-site rank-inverse-normal), z_post adjusted")
    L.append(f"  Outcomes estimable at the 40-event floor: {int(res.p_delta_sd.notna().sum())} "
             f"({len(fam)} non-control incl death, "
             f"{int((res.negative_control & res.p_delta_sd.notna()).sum())} negative controls).")
    L.append(f"  Direction across the {len(fam)} non-control outcomes: HR>1 (more correction,")
    L.append(f"  more risk, equivalently worse baseline) in {n_up}, HR<1 in {n_dn}.")
    L.append(f"  Significant after BH-FDR (q<0.05): {len(sig)}.")
    if len(sig):
        L.append("")
        L.append(f"  {'outcome':<28}{'events':>7}   {'HR/SD (95% CI)':<20}{'p':>10}{'q':>10}")
        for _, r in sig.iterrows():
            L.append(f"  {r.outcome:<28}{r.events:>7}   "
                     f"{fmt(r.hr_delta_sd, r.lo_delta_sd, r.hi_delta_sd):<20}"
                     f"{r.p_delta_sd:>10.2e}{r.q_delta_sd:>10.4f}")
    L.append("")
    top = fam.sort_values("p_delta_sd").head(8)
    L.append("  Smallest p regardless of FDR, with the reparametrization partner shown")
    L.append("  (z_pre given z_post on the same rows). The two columns being near-equal")
    L.append("  is the identity above, seen in the estimates:")
    L.append(f"  {'outcome':<28}{'delta/SD':<22}{'pre/SD (partner)':<22}")
    for _, r in top.iterrows():
        L.append(f"  {r.outcome:<28}{fmt(r.hr_delta_sd, r.lo_delta_sd, r.hi_delta_sd):<22}"
                 f"{fmt(r.hr_pre_sd, r.lo_pre_sd, r.hi_pre_sd):<22}")
    L.append("")
    ncsub = res[res.negative_control & res.p_delta_sd.notna()]
    if len(ncsub):
        L.append("  Negative controls (outside the FDR family): "
                 + ", ".join(f"{r.outcome} {r.hr_delta_sd:.2f} (p={r.p_delta_sd:.2f})"
                             for _, r in ncsub.iterrows()))
    L.append("")

    L.append("B. SENSITIVITIES")
    both = fam.sort_values("p_delta_sd").head(8)
    L.append("  Raw per 10 percentage points of delta (post also raw per 10 pp), and")
    L.append("  z_delta with log1p(post) adjustment, same rows:")
    L.append(f"  {'outcome':<28}{'per 10pp':<22}{'z_delta | log1p(post)':<22}")
    for _, r in both.iterrows():
        a = fmt(r.hr_delta_raw10, r.lo_delta_raw10, r.hi_delta_raw10) \
            if np.isfinite(r.get("hr_delta_raw10", np.nan)) else "-"
        b = fmt(r.hr_delta_lpost, r.lo_delta_lpost, r.hi_delta_lpost) \
            if np.isfinite(r.get("hr_delta_lpost", np.nan)) else "-"
        L.append(f"  {r.outcome:<28}{a:<22}{b:<22}")
    agree = int(((fam.p_delta_sd < .05) == (fam.p_delta_lpost < .05)).sum())
    L.append(f"  p<0.05 agreement between primary and log1p(post) coding: {agree} of {len(fam)}.")
    L.append("")

    L.append("C. SANITY, INDEPENDENT CODE PATH")
    L.append("  Crude rates in the fixed residual stratum post<=5%, delta tertiles.")
    L.append("  median_pre per tertile shows delta and pre moving together, which is the")
    L.append("  identity made visible:")
    for lab in cr.outcome.unique():
        sub = cr[cr.outcome == lab]
        for _, r in sub.iterrows():
            L.append(f"  {lab:<26}{r.delta_tertile:<9} n={r.n:<5} delta~{r.median_delta:<6} "
                     f"pre~{r.median_pre:<6} post~{r.median_post:<5} "
                     f"{r.rate_per_1000py}/1000py ({r.events} ev)")
    L.append("")

    L.append("POSITIVE CONTROL")
    L.append(f"  Same gated loader as the continuous extension: cohort {counts['n']:,} = "
             f"{counts['n_corrected']:,} + {counts['n_not']:,}, designs "
             f"{counts['designs']}, all {pc_n} shipped binary outcomes matched")
    L.append("  numbers/treatment_v2.json exactly before this analysis ran.")
    L.append("")
    L.append("SPECIFICATION")
    L.append("  Cohort and loader identical to run_extension.py. Covariates: age (linear,")
    L.append("  original treatment convention), male, site strata. No log1p(pre) term,")
    L.append("  because pre is spanned by post + delta. Unpenalized lifelines fits,")
    L.append("  40-event floor, outcome set = the 53 defined conditions + death (49 clear")
    L.append("  the floor here), BH-FDR across the")
    L.append("  estimable non-control outcomes.")
    L.append("")
    L.append("PATHS")
    for f in ["delta_run.py", "delta_results.csv", "delta_identity_check.csv",
              "delta_crude.csv", "delta_REPORT.txt"]:
        L.append(f"  {OUTDIR}/{f}")
    (OUTDIR / "delta_REPORT.txt").write_text("\n".join(L) + "\n")
    return "\n".join(L)


def main():
    sev, counts = build_sev()
    pc, mismatches, fo = positive_control(sev, counts)
    print(f"positive control: {int(pc.match.sum())}/{len(pc)} match, "
          f"{len(mismatches)} mismatches")
    if mismatches:
        for msg in mismatches:
            print("  MISMATCH:", msg)
        sys.exit(1)
    res, ident, sev2 = run_delta(sev, fo)
    col = collinearity(sev2)
    cr = crude_delta()
    print(write_report(counts, res, ident, col, cr, len(pc)))


if __name__ == "__main__":
    main()
