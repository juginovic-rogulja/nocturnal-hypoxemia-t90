"""
Compute the three quantities Alen's 2026-08-08 decisions require.

B5   Refit MrOS on FUVSDT, days from the sleep visit, instead of FUCDTIME, days from the
     parent-cohort baseline. numbers/freeze_05_external.py has been corrected in place, and
     two further defects in that file were corrected with it: the age basis returned all four
     spline columns when the four sum to one, and eight fits carried a ridge penalty that
     lifelines scales by the sample size. A third defect surfaced while doing so, and is fixed
     here as well: the spline basis was named age_s*, which overwrote SHHS's own age_s1, the
     man's age at visit 1, so the clinical risk model was adjusting for a spline basis column
     rather than for age. This step runs the corrected script and reports old against new.

B9.1 "About a third of the treated group would have been classified as restored by measurement
     noise alone." Computed from the 1,596 placebo splits: PAP-off nights cut at the same clock
     time as the real switch-on, so they carry no treatment and measure only how often T90 falls
     across a night on its own. The treated arm's rule is applied unchanged, T90 above 10% in the
     first segment falling to 10% or below in the second.
     Writes numbers/placebo_restoration_rate.json.

B9.2 "Their removal strengthened every association rather than creating one." Every condition is
     fitted twice, once on the 19,383 cohort that still contains the 210 corrupted oximetry
     recordings and once on the 19,173 cohort after the quality rule removes them, and the two
     sets of estimates are compared.
     Writes numbers/qc_rule_effect.json.

Every Cox fit here is by maximum partial likelihood with no penalty. Several Newton step sizes
are scanned and the highest-likelihood solution is kept. The age basis is a natural cubic spline
with 4 degrees of freedom, one column dropped.

Nothing is deleted. The previous freeze_05_external.py and external.json are kept under
numbers/_superseded_2026-08-08_b5/.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import json
import os
import subprocess
import sys
import warnings

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from patsy import dmatrix
from scipy import stats
import sys as _sys_v8; _sys_v8.path.insert(0, paths.ANALYSIS_DIR)  # v8 sweep 2026-09-12
from cohort_spec import COHORT_N  # v8 sweep 2026-09-12

warnings.filterwarnings("ignore")

ROOT = paths.T90_ROOT
NUM = f"{paths.NUMBERS_DIR}"
HERE = f"{paths.SV_ROOT}/final_verification"
sys.path.insert(0, paths.ANALYSIS_DIR)

STEPS = [None, 0.05, 0.1, 0.25, 0.5, 0.95]
PRE_QC_N, POST_QC_N = 19383, COHORT_N   # v8 sweep: the pre-QC size is the historical ranking_v2 cohort, the current one comes from cohort_spec


def rule(s):
    print("\n" + "=" * 78 + f"\n{s}\n" + "=" * 78)


def fit_ml(f, tcol, ecol, strata=None):
    """Unpenalized Cox fit. Scan Newton step sizes, keep the maximum-likelihood solution."""
    best = None
    for st in STEPS:
        kw = {} if st is None else {"fit_options": {"step_size": st}}
        try:
            c = CoxPHFitter(penalizer=0.0).fit(f, tcol, ecol, strata=strata, **kw)
        except Exception:
            continue
        ll = float(c.log_likelihood_)
        if not np.isfinite(ll):
            continue
        if best is None or ll > best[0] + 1e-9:
            best = (ll, c, st)
    if best is None:
        return None, None
    return best[1], best[2]


# ==================================================================== B5
def b5_external():
    rule("B5  MrOS refit on FUVSDT, plus the spline and penalty defects in the same file")

    old_path = f"{NUM}/_superseded_2026-08-08_b5/external.json.pre_b5"
    old = json.load(open(old_path)) if os.path.exists(old_path) else \
        json.load(open(f"{NUM}/external.json"))

    print("running numbers/freeze_05_external.py (corrected) ...")
    r = subprocess.run([sys.executable, "freeze_05_external.py"], cwd=NUM,
                       capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-4000:])
        print(r.stderr[-4000:])
        raise SystemExit("freeze_05_external.py failed")
    new = json.load(open(f"{NUM}/external.json"))

    cmp_rows = []
    for coh in ("mros", "shhs"):
        for lab, v in new[coh]["outcomes"].items():
            ov = old.get(coh, {}).get("outcomes", {}).get(lab, {})
            for est in ("t90", "ahi", "t90_adj_ahi", "ahi_adj_t90"):
                if est not in v or est not in ov:
                    continue
                cmp_rows.append({
                    "cohort": coh, "outcome": lab, "estimate": est,
                    "events": v[est]["events"], "n": v[est]["n"],
                    "old_hr": ov[est]["hr"], "old_lo": ov[est]["lo"], "old_hi": ov[est]["hi"],
                    "old_p": ov[est]["p"],
                    "new_hr": v[est]["hr"], "new_lo": v[est]["lo"], "new_hi": v[est]["hi"],
                    "new_p": v[est]["p"],
                    "delta_hr": round(v[est]["hr"] - ov[est]["hr"], 4),
                    "pct_change_hr": round(100 * (v[est]["hr"] / ov[est]["hr"] - 1), 2),
                    "sig_old": (ov[est]["lo"] > 1) or (ov[est]["hi"] < 1),
                    "sig_new": (v[est]["lo"] > 1) or (v[est]["hi"] < 1)})
    cmp_df = pd.DataFrame(cmp_rows)
    cmp_df.to_csv(f"{NUM}/external_refit_comparison.csv", index=False)

    fu_new = new["mros"]["followup_years_median"]
    fu_old = new["mros"]["followup_years_median_superseded_FUCDTIME"]
    print(f"\nMrOS n={new['mros']['n']:,}  sites={new['mros']['sites']}")
    print(f"median follow-up  FUVSDT {fu_new} y   (superseded FUCDTIME {fu_old} y)")
    print(f"IQR {new['mros']['followup_years_iqr'][0]}-{new['mros']['followup_years_iqr'][1]} y\n")
    print(f"{'MrOS outcome':<26}{'ev':>6}{'old HR (95% CI)':>24}{'new HR (95% CI)':>24}{'shift':>9}")
    print("-" * 90)
    for lab in new["mros"]["outcomes"]:
        row = cmp_df[(cmp_df.cohort == "mros") & (cmp_df.outcome == lab) &
                     (cmp_df.estimate == "t90")]
        if not len(row):
            continue
        r0 = row.iloc[0]
        o = f"{r0.old_hr:.3f} [{r0.old_lo:.3f}-{r0.old_hi:.3f}]"
        n_ = f"{r0.new_hr:.3f} [{r0.new_lo:.3f}-{r0.new_hi:.3f}]"
        print(f"{lab:<26}{r0.events:>6,}{o:>24}{n_:>24}{r0['pct_change_hr']:>8.2f}%")

    flips = cmp_df[cmp_df.sig_old != cmp_df.sig_new]
    print(f"\nestimates whose significance changed: {len(flips)}")
    for _, r0 in flips.iterrows():
        print(f"  {r0.cohort.upper():<5}{r0.outcome:<28}{r0.estimate:<14}"
              f"{r0.old_hr:.3f}[{r0.old_lo:.3f}-{r0.old_hi:.3f}] -> "
              f"{r0.new_hr:.3f}[{r0.new_lo:.3f}-{r0.new_hi:.3f}]")

    sh = cmp_df[cmp_df.cohort == "shhs"]
    print(f"\nSHHS side effect of dropping the penalty and the spline column: "
          f"median |shift| {sh.pct_change_hr.abs().median():.2f}%, "
          f"max {sh.pct_change_hr.abs().max():.2f}%")
    for lab, v in new["shhs"]["outcomes"].items():
        ov = old["shhs"]["outcomes"].get(lab, {})
        if "clinical_model_gain" in v and "clinical_model_gain" in ov:
            print(f"  C-statistic gain, {lab:<26} {ov['clinical_model_gain']:+.4f} -> "
                  f"{v['clinical_model_gain']:+.4f}")
    return new, old, cmp_df


# ==================================================================== B9.1
def b9_placebo():
    rule("B9.1  restoration under no treatment, from the placebo splits")

    d = pd.read_parquet(f"{paths.TABLES_DIR}/cpap_t90_by_stage.parquet")
    base = pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet")

    def quality(g, sides):
        """The cohort file's oximetry rule, applied to each segment of the night."""
        for sd in sides:
            n, mn, t = f"{sd}_sleep_nadir", f"{sd}_sleep_mean", f"{sd}_sleep_t90"
            if n in g.columns and mn in g.columns:
                g = g[~(g[n] > g[mn])]
            if t in g.columns and mn in g.columns:
                g = g[~(g[mn].between(50, 80) & (g[t] > 70))]
        return g

    def ladder(sub, first, second):
        """Report how many nights each filter costs, so the denominator is auditable."""
        steps = [("on file", len(sub))]
        g = sub[(sub[f"{first}_sleep_min"] >= 30) & (sub[f"{second}_sleep_min"] >= 30)]
        steps.append(("at least 30 min of sleep in each segment", len(g)))
        g = g[g[f"{first}_sleep_t90"].notna() & g[f"{second}_sleep_t90"].notna()]
        steps.append(("T90 measurable in each segment", len(g)))
        g = quality(g, [first, second])
        steps.append(("oximetry quality rule on each segment", len(g)))
        g = g.merge(base[["BDSPPatientID", "fu_valid"]], on="BDSPPatientID", how="inner")
        g = g[g.fu_valid == 1]
        steps.append(("usable follow-up, same as the treated arm", len(g)))
        g = g.drop_duplicates("BDSPPatientID", keep="first")
        steps.append(("one night per patient", len(g)))
        return g, [{"step": s, "n": int(n)} for s, n in steps]

    def rate(g, first, second):
        hi = g[g[f"{first}_sleep_t90"] > 10]
        rest = int((hi[f"{second}_sleep_t90"] <= 10).sum())
        ci = stats.binomtest(rest, len(hi)).proportion_ci(0.95) if len(hi) else None
        return {"n_eligible": int(len(g)), "n_above10_at_baseline": int(len(hi)),
                "n_restored": rest,
                "pct_restored": round(100 * rest / len(hi), 1) if len(hi) else None,
                "ci95_pct": [round(100 * ci.low, 1), round(100 * ci.high, 1)] if ci else None,
                "baseline_t90_median": round(float(hi[f"{first}_sleep_t90"].median()), 2)
                if len(hi) else None,
                "second_t90_median": round(float(hi[f"{second}_sleep_t90"].median()), 2)
                if len(hi) else None}, hi

    P = d[d.design == "P_placebo"].copy()

    # the treated arm is rebuilt exactly as numbers/run_treatment_v2.py builds it: the
    # split-night design plus the paired-study design, the latter mapped from dx/tx onto
    # pre/post, one night per patient with the split night preferred
    S = ["wake", "sleep", "nrem", "n1", "n2", "n3", "rem"]
    A = d[d.design == "A_split"].copy()
    A["src"] = "A"
    B = d[d.design == "B_pairs"].copy()
    B["src"] = "B"
    for s in S:
        for a, bb in [("pre", "dx"), ("post", "tx")]:
            for suf in ["t90", "min", "mean", "nadir"]:
                c = f"{bb}_{s}_{suf}"
                B[f"{a}_{s}_{suf}"] = B[c] if c in B.columns else np.nan
    keep = ["BDSPPatientID", "src"] + [f"{p}_{s}_{x}" for p in ("pre", "post")
                                       for s in S for x in ("t90", "min", "mean", "nadir")]
    T = pd.concat([A[[c for c in keep if c in A.columns]],
                   B[[c for c in keep if c in B.columns]]], ignore_index=True)
    T = T.sort_values("src", ascending=False).drop_duplicates("BDSPPatientID", keep="first")

    pg, plad = ladder(P, "h1", "h2")
    tg, tlad = ladder(T, "pre", "post")
    plac, phi = rate(pg, "h1", "h2")
    treat, thi = rate(tg, "pre", "post")

    # reconciliation against the published treated arm. run_treatment_v2.py names its quality
    # rule on columns that its keep-list never carries, so the rule is a no-op there and its
    # denominator is looser than the one used above.
    Tp = T.merge(base[["BDSPPatientID", "fu_valid"]], on="BDSPPatientID", how="inner")
    Tp = Tp[(Tp.fu_valid == 1) & (Tp.post_sleep_min >= 30) & (Tp.pre_sleep_min >= 30)]
    hi_p = Tp[Tp.pre_sleep_t90 > 10]
    pub = {"n_above10_at_baseline": int(len(hi_p)),
           "n_restored": int((hi_p.post_sleep_t90 <= 10).sum()),
           "pct_restored": round(100 * float((hi_p.post_sleep_t90 <= 10).mean()), 1),
           "note": "matches treatment_v2.json corrected_vs_not: n 1894, corrected 1326"}

    # the treated group starts sicker, and how far T90 drifts on its own depends on where it
    # starts, so the placebo rate is also standardised to the treated baseline distribution
    edges = [10, 15, 25, 50, 1e9]
    names = ["10-15%", "15-25%", "25-50%", ">50%"]
    strata, num, den = [], 0.0, 0
    for lo, hi_, nm in zip(edges[:-1], edges[1:], names):
        pm = phi[(phi.h1_sleep_t90 > lo) & (phi.h1_sleep_t90 <= hi_)]
        tm = thi[(thi.pre_sleep_t90 > lo) & (thi.pre_sleep_t90 <= hi_)]
        pr = float((pm.h2_sleep_t90 <= 10).mean()) if len(pm) else None
        tr = float((tm.post_sleep_t90 <= 10).mean()) if len(tm) else None
        strata.append({"baseline_band": nm, "placebo_n": int(len(pm)),
                       "placebo_restored": int((pm.h2_sleep_t90 <= 10).sum()) if len(pm) else 0,
                       "placebo_pct": round(100 * pr, 1) if pr is not None else None,
                       "treated_n": int(len(tm)),
                       "treated_restored": int((tm.post_sleep_t90 <= 10).sum()) if len(tm) else 0,
                       "treated_pct": round(100 * tr, 1) if tr is not None else None})
        if pr is not None and len(tm):
            num += pr * len(tm)
            den += len(tm)
    standardised = round(100 * num / den, 1) if den else None

    clk = {"placebo_cut_min_median": round(float(P.split_min.median()), 1),
           "placebo_cut_min_iqr": [round(float(P.split_min.quantile(.25)), 1),
                                   round(float(P.split_min.quantile(.75)), 1)],
           "treated_switch_on_min_median": round(float(A.cpap_start_min.median()), 1),
           "treated_switch_on_min_iqr": [round(float(A.cpap_start_min.quantile(.25)), 1),
                                         round(float(A.cpap_start_min.quantile(.75)), 1)],
           "placebo_max_fraction_of_night_on_pap": float(P.cpap_frac.max())}

    p_rate = plac["pct_restored"] / 100
    t_rate = treat["pct_restored"] / 100
    out = {
        "question": "what fraction of the treated group would be called restored by measurement "
                    "noise alone, judged from untreated nights split at the same clock time",
        "rule": "sleep T90 above 10% in the first segment, 10% or below in the second",
        "placebo_designs_on_file": int((d.design == "P_placebo").sum()),
        "placebo": plac,
        "placebo_filter_ladder": plad,
        "treated": treat,
        "treated_filter_ladder": tlad,
        "treated_as_published_in_treatment_v2": pub,
        "clock": clk,
        "by_baseline_band": strata,
        "placebo_restoration_rate_pct": plac["pct_restored"],
        "placebo_restoration_rate_standardised_to_treated_baseline_pct": standardised,
        "treated_restoration_rate_pct": treat["pct_restored"],
        "excess_over_noise_pct_points": round(100 * (t_rate - p_rate), 1),
        "share_of_treated_restorations_noise_alone_reproduces_pct":
            round(100 * p_rate / t_rate, 1) if t_rate else None,
        "filters": "sleep at least 30 min in each segment, the cohort oximetry quality rule "
                   "applied to each segment, usable follow-up, one night per patient",
        "caveat": "the untreated nights are a healthier group, so only "
                  f"{plac['n_above10_at_baseline']} of them start above 10% and the rate rests "
                  "on that many observations",
    }
    json.dump(out, open(f"{NUM}/placebo_restoration_rate.json", "w"), indent=2)

    print("placebo filter ladder")
    for s in plad:
        print(f"  {s['step']:<46}{s['n']:>7,}")
    print(f"\n  above 10% in segment one             {plac['n_above10_at_baseline']:>7,}")
    print(f"  falling to 10% or below by segment two {plac['n_restored']:>5,}"
          f"  = {plac['pct_restored']}%  "
          f"[95% CI, {plac['ci95_pct'][0]}-{plac['ci95_pct'][1]}]")
    print(f"\ntreated arm, identical rule            {treat['n_above10_at_baseline']:,} above 10%,"
          f" {treat['n_restored']:,} restored = {treat['pct_restored']}%"
          f" [95% CI, {treat['ci95_pct'][0]}-{treat['ci95_pct'][1]}]")
    print(f"treated arm as published                {pub['n_above10_at_baseline']:,} above 10%,"
          f" {pub['n_restored']:,} restored = {pub['pct_restored']}%")
    print(f"placebo rate standardised to the treated baseline distribution  {standardised}%")
    print(f"excess over noise                      "
          f"{out['excess_over_noise_pct_points']} percentage points")
    print(f"share of treated restorations that noise alone reproduces  "
          f"{out['share_of_treated_restorations_noise_alone_reproduces_pct']}%")
    print(f"\n{'baseline band':<16}{'placebo':>18}{'treated':>18}")
    for s in strata:
        pp = f"{s['placebo_restored']}/{s['placebo_n']} = {s['placebo_pct']}%"
        tt = f"{s['treated_restored']}/{s['treated_n']} = {s['treated_pct']}%"
        print(f"{s['baseline_band']:<16}{pp:>18}{tt:>18}")
    print(f"\nplacebo cut at {clk['placebo_cut_min_median']} min, real switch-on at "
          f"{clk['treated_switch_on_min_median']} min, "
          f"most PAP on any placebo night {100*clk['placebo_max_fraction_of_night_on_pap']:.2f}%")
    print(f"written -> {NUM}/placebo_restoration_rate.json")
    return out


# ==================================================================== B9.2
def b9_qc():
    rule("B9.2  every association before and after the oximetry quality rule")

    from disease_definitions import DISEASES, NEGATIVE_CONTROLS

    b = pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet")
    pre = b[(b.fu_valid == 1) & b.spo2_pct_below_90.notna()].copy()
    post = pre[pre.oximetry_bad == 0].copy()
    assert len(pre) == PRE_QC_N, f"pre-rule cohort is {len(pre):,}, expected {PRE_QC_N:,}"
    assert len(post) == POST_QC_N, f"post-rule cohort is {len(post):,}, expected {POST_QC_N:,}"
    print(f"before the rule {len(pre):,}   after {len(post):,}   removed {len(pre)-len(post)}")

    def prep(g):
        g = g.copy()
        g["male"] = (g.sex.astype(str).str.upper().str[0] == "M").astype(int)
        sp = dmatrix("cr(a, df=4) - 1", {"a": g.AgeAtVisit.values},
                     return_type="dataframe").iloc[:, 1:]
        for i in range(sp.shape[1]):
            g[f"age_s{i}"] = sp.iloc[:, i].values
        adj = [f"age_s{i}" for i in range(sp.shape[1])] + ["male"]

        def rint(x):
            r = stats.rankdata(x)
            return stats.norm.ppf((r - 0.375) / (len(r) + 0.25))

        g["z"] = g.groupby("site_id").spo2_pct_below_90.transform(rint)
        return g, adj

    pre, ADJ = prep(pre)
    post, _ = prep(post)

    def one(g, adj, key):
        yc, ec, pc = f"{key}_years", f"{key}_incident", f"{key}_prevalent"
        if yc not in g.columns:
            return None
        f = g[(g[pc] == 0) & g[yc].notna() & (g[yc] > 0)]
        dd = pd.DataFrame({"T": f[yc], "E": f[ec].astype(int), "site": f.site_id,
                           **{c: f[c] for c in ["z"] + adj}}).dropna()
        if dd.E.sum() < 60:
            return None
        c, st = fit_ml(dd, "T", "E", ["site"])
        if c is None:
            return None
        r = c.summary.loc["z"]
        return {"n": int(len(dd)), "events": int(dd.E.sum()),
                "hr": round(float(r["exp(coef)"]), 4),
                "lo": round(float(r["exp(coef) lower 95%"]), 4),
                "hi": round(float(r["exp(coef) upper 95%"]), 4),
                "p": float(r["p"]), "coef": float(r["coef"]), "step_size": str(st)}

    keys = list(DISEASES.keys()) + ["death"]
    labels = {k: v[0] for k, v in DISEASES.items()}
    labels["death"] = "Death from any cause"

    rows = []
    for i, k in enumerate(keys, 1):
        a = one(pre, ADJ, k)
        c = one(post, ADJ, k)
        if a is None or c is None:
            print(f"  [{i:>2}/{len(keys)}] {labels.get(k, k):<30} not fitted in both cohorts")
            continue
        sig_a = (a["lo"] > 1) or (a["hi"] < 1)
        sig_c = (c["lo"] > 1) or (c["hi"] < 1)
        rows.append({
            "key": k, "condition": labels.get(k, k),
            "negative_control": k in NEGATIVE_CONTROLS,
            "events_pre": a["events"], "n_pre": a["n"],
            "events_post": c["events"], "n_post": c["n"],
            "hr_pre": a["hr"], "lo_pre": a["lo"], "hi_pre": a["hi"], "p_pre": a["p"],
            "hr_post": c["hr"], "lo_post": c["lo"], "hi_post": c["hi"], "p_post": c["p"],
            "coef_pre": a["coef"], "coef_post": c["coef"],
            "delta_hr": round(c["hr"] - a["hr"], 4),
            "pct_change_hr": round(100 * (c["hr"] / a["hr"] - 1), 2),
            "abs_loghr_pre": abs(a["coef"]), "abs_loghr_post": abs(c["coef"]),
            "moved_away_from_1": abs(c["coef"]) > abs(a["coef"]),
            "hr_rose": c["hr"] > a["hr"],
            "sig_pre": bool(sig_a), "sig_post": bool(sig_c),
            "created_by_removal": bool((not sig_a) and sig_c),
            "lost_by_removal": bool(sig_a and (not sig_c)),
            "sign_flip": bool(np.sign(a["coef"]) != np.sign(c["coef"]))})
        print(f"  [{i:>2}/{len(keys)}] {labels.get(k, k):<30} "
              f"{a['hr']:.3f} -> {c['hr']:.3f}  ({rows[-1]['pct_change_hr']:+.2f}%)")

    t = pd.DataFrame(rows)
    t.to_csv(f"{NUM}/qc_rule_effect.csv", index=False)

    def block(sub, name):
        return {"label": name, "n_outcomes": int(len(sub)),
                "n_strengthened_moved_away_from_1": int(sub.moved_away_from_1.sum()),
                "n_weakened_moved_toward_1": int((~sub.moved_away_from_1).sum()),
                "n_hazard_ratio_rose": int(sub.hr_rose.sum()),
                "n_hazard_ratio_fell": int((~sub.hr_rose).sum()),
                "n_created_by_removal": int(sub.created_by_removal.sum()),
                "n_lost_by_removal": int(sub.lost_by_removal.sum()),
                "n_sign_flip": int(sub.sign_flip.sum()),
                "median_pct_change_hr": round(float(sub.pct_change_hr.median()), 2),
                "min_pct_change_hr": round(float(sub.pct_change_hr.min()), 2),
                "max_pct_change_hr": round(float(sub.pct_change_hr.max()), 2)}

    weak = t[~t.moved_away_from_1]
    created = t[t.created_by_removal]
    out = {
        "claim_under_test": "Their removal strengthened every association rather than creating one.",
        "cohort_before_rule": int(len(pre)), "cohort_after_rule": int(len(post)),
        "recordings_removed": int(len(pre) - len(post)),
        "model": "Cox, site-stratified, natural cubic spline on age with one column dropped, "
                 "sex, rank-based inverse-normal T90 within site, no penalty, "
                 "Newton step size scanned",
        "strengthened_definition": "the hazard ratio moved further from 1 in its own direction",
        "all_outcomes": block(t, "all fitted conditions"),
        "excluding_negative_controls": block(t[~t.negative_control],
                                             "excluding the 5 negative controls"),
        "associations_significant_after_the_rule": block(t[t.sig_post],
                                                         "significant after the rule"),
        "outcomes_that_weakened": weak[["condition", "negative_control", "events_post",
                                        "hr_pre", "hr_post", "pct_change_hr",
                                        "sig_pre", "sig_post"]].to_dict("records"),
        "outcomes_created_by_removal": created[["condition", "negative_control", "events_post",
                                                "hr_pre", "lo_pre", "hi_pre",
                                                "hr_post", "lo_post", "hi_post"]
                                               ].to_dict("records"),
        "per_condition": t.to_dict("records"),
    }
    a = out["all_outcomes"]
    out["verdict"] = (
        "literally true" if a["n_weakened_moved_toward_1"] == 0 and a["n_created_by_removal"] == 0
        else f"not literally true: {a['n_strengthened_moved_away_from_1']} of {a['n_outcomes']} "
             f"strengthened, {a['n_weakened_moved_toward_1']} weakened, "
             f"{a['n_created_by_removal']} became significant only after the removal")
    json.dump(out, open(f"{NUM}/qc_rule_effect.json", "w"), indent=2)

    print(f"\n{'':<44}{'all':>7}{'no controls':>14}{'sig after':>12}")
    for lbl, kk in [("outcomes fitted", "n_outcomes"),
                    ("strengthened, moved away from 1", "n_strengthened_moved_away_from_1"),
                    ("weakened, moved toward 1", "n_weakened_moved_toward_1"),
                    ("hazard ratio rose", "n_hazard_ratio_rose"),
                    ("became significant only after removal", "n_created_by_removal"),
                    ("lost significance after removal", "n_lost_by_removal"),
                    ("changed direction", "n_sign_flip")]:
        print(f"{lbl:<44}{a[kk]:>7}{out['excluding_negative_controls'][kk]:>14}"
              f"{out['associations_significant_after_the_rule'][kk]:>12}")
    print(f"\nmedian change in the hazard ratio {a['median_pct_change_hr']:+.2f}%  "
          f"(range {a['min_pct_change_hr']:+.2f}% to {a['max_pct_change_hr']:+.2f}%)")
    if len(weak):
        print(f"\nthe {len(weak)} that did NOT strengthen")
        for _, r0 in weak.sort_values("pct_change_hr").iterrows():
            tag = "  NEG" if r0.negative_control else ""
            print(f"  {r0.condition:<30}{r0.hr_pre:.3f} -> {r0.hr_post:.3f}"
                  f"{r0.pct_change_hr:>9.2f}%  ev={r0.events_post:,}{tag}")
    if len(created):
        print(f"\nthe {len(created)} that became significant only after the removal")
        for _, r0 in created.iterrows():
            print(f"  {r0.condition:<30}{r0.hr_pre:.3f}[{r0.lo_pre:.3f}-{r0.hi_pre:.3f}] -> "
                  f"{r0.hr_post:.3f}[{r0.lo_post:.3f}-{r0.hi_post:.3f}]")
    print(f"\nverdict: {out['verdict']}")
    print(f"written -> {NUM}/qc_rule_effect.json and qc_rule_effect.csv")
    return out


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("all", "b5"):
        b5_external()
    if which in ("all", "b91"):
        b9_placebo()
    if which in ("all", "b92"):
        b9_qc()
