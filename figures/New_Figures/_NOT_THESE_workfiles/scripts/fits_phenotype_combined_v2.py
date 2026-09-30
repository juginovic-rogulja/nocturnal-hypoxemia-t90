"""
The non-responder phenotype, rebuilt on the corrected treatment frame.

WHY THIS EXISTS
---------------
`fits_phenotype_combined.py` built its sample straight off the cpap_stage pilot CSVs: arm A
only, status ok, pre and post sleep T90 present, pre_sleep_t90 > 10, one row per patient, and
a nadir-must-not-exceed-mean rule. That gives 1,996 patients, 614 of whom did not normalize.

The paper's treatment analysis does not use that frame. `numbers/run_treatment_v2.py` builds it
from `data_frozen_v7_2026-09/cpap_t90_by_stage.parquet`, taking the union of the split-night and paired
designs, requiring fu_valid == 1 and at least 30 minutes of scored sleep on both nights, and
then restricting to pre_sleep_t90 > 10. That gives **1,894 patients, 568 of whom did not
normalize and 1,326 of whom did**, which is what `treatment_v2.json`, Figure 11 and Figure 12
all carry.

Figure 13, eFigure 10 and eTable 12 were still on the 1,996. This script puts them on the
1,894. Nothing else about the specification changes.

WHAT IS REFITTED, AND WHY IT HAD TO BE
--------------------------------------
Every quantity the old file supplied is recomputed here rather than copied, because each one
was frame-dependent:

  * the per-condition marginal odds ratio, adjusted for age, sex and log(1 + pretreatment
    sleep T90), with Benjamini-Hochberg correction across the conditions tested. The old run
    read these from `nonresponder_phenotype_audited.csv`, which is on the 1,996.
  * which conditions clear FDR, which is the candidate set for the joint model.
  * the joint model, the independent subset, the areas under the curve, the risk fifths.
  * the absolute failure rate with and without each condition, which is what panel A prints.
  * the cardiopulmonary dose-response and its trend model, which eFigure 10 panel E prints and
    which the old run copied from `nonresponder_doseresponse.json`, also on the 1,996.

The 5 settled negative controls are fitted alongside the conditions and reported, but they are
never eligible for the joint model.

Writes numbers/nonresponder_combined_v2.json. The v1 file is left where it is.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../..")))  # repository root, where paths.py lives
import paths
import glob
import os
import json
import sys
import warnings

import duckdb
import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")
ROOT = paths.FIGURE_ROOT
sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
from disease_definitions import DISEASES, NEGATIVE_CONTROLS   # noqa: E402
from cohort_spec import (CPAP_STAGE_PARQUET, T90_FINAL, NUMBERS_DIR, OUT_DIR,  # noqa: E402  (v8 sweep 2026-09-12)
                         PREVALENT_CARDIOPULMONARY, sidecar)

CP = ["resp_failure", "copd2", "obesity_hypovent", "hf", "pulm_htn", "asthma"]
STAGES = ["wake", "sleep", "nrem", "n1", "n2", "n3", "rem"]

# ------------------------------------------------------------------ the sample
# Reproduced from numbers/run_treatment_v2.py, verbatim in effect, so the two cannot drift.
d = pd.read_parquet(CPAP_STAGE_PARQUET)
A = d[d.design == "A_split"].copy(); A["src"] = "A"
B = d[d.design == "B_pairs"].copy(); B["src"] = "B"
for s in STAGES:
    for a, b in [("pre", "dx"), ("post", "tx")]:
        for suf in ("t90", "min"):
            c = f"{b}_{s}_{suf}"
            B[f"{a}_{s}_{suf}"] = B[c] if c in B.columns else np.nan
keep = (["BDSPPatientID", "SexDSC", "src"]
        + [f"{p}_{s}_{x}" for p in ("pre", "post") for s in STAGES for x in ("t90", "min")])
cc = pd.concat([A[[c for c in keep if c in A.columns]],
                B[[c for c in keep if c in B.columns]]], ignore_index=True)
cc = cc.sort_values("src", ascending=False).drop_duplicates("BDSPPatientID", keep="first")

base = pd.read_parquet(T90_FINAL).drop_duplicates("BDSPPatientID")
m = cc.merge(base, on="BDSPPatientID", how="inner")
m = m[(m.fu_valid == 1) & (m.post_sleep_min >= 30) & (m.pre_sleep_min >= 30)].copy()
m = m[m.pre_sleep_t90 > 10].copy()
m["nonresp"] = (m.post_sleep_t90 > 10).astype(int)

TREAT = json.load(open(f"{NUMBERS_DIR}/treatment_v2.json"))["corrected_vs_not"]
assert len(m) == TREAT["n"], f"{len(m)} against treatment_v2 {TREAT['n']}"   # v8: the reference is the regenerated treatment_v2.json, no typed count
assert int(m.nonresp.sum()) == TREAT["n_not"], int(m.nonresp.sum())
assert int((1 - m.nonresp).sum()) == TREAT["n_corrected"]
print(f"sample {len(m):,}, of whom {int(m.nonresp.sum()):,} "
      f"({100 * m.nonresp.mean():.1f}%) did not normalize, "
      f"{int((1 - m.nonresp).sum()):,} did")

pcols = [f"{k}_prevalent" for k in DISEASES if f"{k}_prevalent" in m.columns]
for c in pcols:
    m[c] = m[c].fillna(0).astype(int)
m["male"] = (m.SexDSC.astype(str).str.upper().str[0] == "M").astype(int)
m["lbase"] = np.log1p(m.pre_sleep_t90)
m["cp_raw"] = sum(m[f"{k}_prevalent"] for k in CP)
m["cp"] = m.cp_raw.clip(upper=3)
# decision 6 (2026-09-12): the prevalent cardiopulmonary covariate, any condition of the D4 set
CPX = [k for k in PREVALENT_CARDIOPULMONARY if f"{k}_prevalent" in m.columns]
m["cp_any"] = m[[f"{k}_prevalent" for k in CPX]].max(axis=1).astype(int)

# ------------------------------------------------------------------ race
PERSON = glob.glob(f"{paths.OMOP_CACHE_DIR}/**/person_merged.parquet",
                   recursive=True)[0]
p = duckdb.sql(f"""select person_id as "BDSPPatientID", race_source_value
                   from '{PERSON}'""").df()
p = p.drop_duplicates("BDSPPatientID")
m = m.merge(p, on="BDSPPatientID", how="left")
raw = m.race_source_value.astype(str).str.strip().str.upper()
BLACK = {"AFRICAN AMERICAN  OR BLACK", "AFRICAN AMERICAN OR BLACK",
         "BLACK/AFRICAN AMERICAN", "BLACK", "AFRICAN AMERICAN"}
WHITE = {"WHITE", "CAUCASIAN OR WHITE", "CAUCASIAN"}
UNK = {"UNKNOWN, UNAVAILABLE OR UNREPORTED", "UNKNOWN/NOT SPECIFIED", "UNKNOWN",
       "UNABLE TO OBTAIN", "DECLINED TO ANSWER", "NONE", "PREFER NOT TO SAY", "NAN", "",
       "DECLINED", "NOT RECORDED", "PATIENT DECLINED", "NOT SPECIFIED", "UNAVAILABLE",
       "DECLINE TO ANSWER"}
m["race_black"] = raw.isin(BLACK).astype(int)
m["race_asian"] = (raw == "ASIAN").astype(int)
m["race_other"] = (~raw.isin(BLACK | WHITE | UNK) & (raw != "ASIAN")).astype(int)
m["race_unknown"] = raw.isin(UNK).astype(int)
RACE = ["race_black", "race_asian", "race_other", "race_unknown"]
print("race:", {"White": int(raw.isin(WHITE).sum()), "Black": int(m.race_black.sum()),
                "Asian": int(m.race_asian.sum()), "Other": int(m.race_other.sum()),
                "Unknown": int(m.race_unknown.sum())})
DEMO = ["AgeAtVisit", "male"] + RACE


def logit(terms, data=None):
    dd = m if data is None else data
    X = sm.add_constant(dd[terms].astype(float))
    f = pd.concat([X, dd.nonresp], axis=1).dropna()
    r = sm.Logit(f.nonresp, f[X.columns]).fit(disp=0)
    ci = r.conf_int()
    return r, {t: {"or": float(np.exp(r.params[t])),  # v8.3 full precision (2026-09-16)
                   "lo": float(np.exp(ci.loc[t, 0])),
                   "hi": float(np.exp(ci.loc[t, 1])),
                   "p": float(r.pvalues[t])} for t in terms}, int(len(f))


def bh(pvals):
    """Benjamini-Hochberg q values, in the order given."""
    n = len(pvals)
    order = np.argsort(pvals)
    q = np.empty(n, float)
    prev = 1.0
    for rank, i in enumerate(reversed(order), start=1):
        val = pvals[i] * n / (n - rank + 1)
        prev = min(prev, val)
        q[i] = min(prev, 1.0)
    return q


# ------------------------------------------------------------------ marginal per-condition fits
# A condition is testable when at least 20 patients carry it and at least 20 do not, and both
# cells contain at least one failure and one success. Below that a logistic fit is separated and
# reports an odds ratio that is an artifact of the separation, not an estimate.
MIN_CELL = 20
label_of = {k: DISEASES[k][0] for k in DISEASES}
rows = []
for k in DISEASES:
    c = f"{k}_prevalent"
    if c not in m.columns:
        continue
    nw, nwo = int((m[c] == 1).sum()), int((m[c] == 0).sum())
    if nw < MIN_CELL or nwo < MIN_CELL:
        continue
    yw = m.loc[m[c] == 1, "nonresp"]
    if yw.sum() == 0 or yw.sum() == len(yw):
        continue
    _f, res, _n = logit([c, "AgeAtVisit", "male", "lbase"])
    v = res[c]
    # decision 6: the same marginal odds ratio adjusted for any OTHER prevalent cardiopulmonary condition
    others = [f"{j}_prevalent" for j in CPX if j != k]
    m["cp_any_x"] = m[others].max(axis=1).astype(int) if others else 0
    _f2, res2, _n2 = logit([c, "AgeAtVisit", "male", "lbase", "cp_any_x"])
    v2 = res2[c]
    rows.append({"key": k, "Condition": label_of[k], "neg": k in NEGATIVE_CONTROLS,
                 "n_with": nw, "n_without": nwo,
                 "n_failed_with": int(yw.sum()),
                 "n_failed_without": int(m.loc[m[c] == 0, "nonresp"].sum()),
                 "fail_pct_with": round(100 * float(yw.mean()), 1),
                 "fail_pct_without": round(100 * float(m.loc[m[c] == 0, "nonresp"].mean()), 1),
                 "marginal_or": v["or"], "marginal_lo": v["lo"], "marginal_hi": v["hi"],
                 "marginal_p": v["p"],
                 "cp_adj_or": v2["or"], "cp_adj_lo": v2["lo"], "cp_adj_hi": v2["hi"], "cp_adj_p": v2["p"]})
marg = pd.DataFrame(rows)
# FDR is applied across the disease outcomes. The negative controls are reported but are not
# part of the family being corrected, so they cannot move anyone else's q value. They get their
# own correction across the 5 of them, so a control still has a q a figure can test rather than
# a blank that every comparison silently answers False.
qs = []
for flag in (False, True):
    sub = marg[marg.neg == flag].copy()
    sub["marginal_q"] = bh(sub.marginal_p.values)
    sub["cp_adj_q"] = bh(sub.cp_adj_p.values)
    qs.append(sub[["key", "marginal_q", "cp_adj_q"]])
marg = marg.merge(pd.concat(qs, ignore_index=True), on="key", how="left")
print(f"conditions testable at the {MIN_CELL}-patient floor: {len(marg)} "
      f"({int(marg.neg.sum())} of them negative controls)")

SIG = [r.key for r in marg.itertuples() if not r.neg and r.marginal_q < 0.05]
label = {k: label_of[k] for k in SIG}
print(f"conditions significant after FDR, adjusted for age, sex and baseline T90: {len(SIG)}")

# ------------------------------------------------------------------ the joint model
JOINT_TERMS = [f"{k}_prevalent" for k in SIG] + DEMO + ["lbase"]
_r, joint, njoint = logit(JOINT_TERMS)
indep = [k for k in SIG if joint[f"{k}_prevalent"]["p"] < 0.05]
print(f"still independent with everything entered together: {len(indep)}")
for k in sorted(indep, key=lambda k: -joint[f"{k}_prevalent"]["or"]):
    v = joint[f"{k}_prevalent"]
    print(f"   {label[k]:<28}{v['or']:.2f} (95% CI, {v['lo']:.2f}-{v['hi']:.2f})")

# ------------------------------------------------------------------ models and areas
MODELS = {
    "Pretreatment oxygen alone": ["lbase"],
    "Demographics alone": DEMO,
    "Pretreatment oxygen + demographics": ["lbase"] + DEMO,
    "Cardiopulmonary count + oxygen + demographics": ["cp", "lbase"] + DEMO,
    "Prevalent cardiopulmonary disease (any, decision 6 set) + oxygen + demographics": ["cp_any", "lbase"] + DEMO,
    "All FDR-significant conditions + oxygen + demographics": JOINT_TERMS,
    "Independent conditions + oxygen + demographics":
        [f"{k}_prevalent" for k in indep] + DEMO + ["lbase"],
}


def areas(terms, seed=20260807):
    X = m[terms].astype(float).values
    y = m.nonresp.values
    pipe = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000, C=1e6))
    pipe.fit(X, y)
    app = roc_auc_score(y, pipe.predict_proba(X)[:, 1])
    oof = np.zeros(len(y))
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    for tr, te in cv.split(X, y):
        pp = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000, C=1e6))
        pp.fit(X[tr], y[tr])
        oof[te] = pp.predict_proba(X[te])[:, 1]
    return float(app), float(roc_auc_score(y, oof)), oof  # v8.3 full precision (2026-09-16)


print(f"\n{'model':<56}{'terms':>6}{'apparent':>10}{'5-fold CV':>11}")
AUC, oof_store = {}, {}
for name, terms in MODELS.items():
    app, cvv, oof = areas(terms)
    AUC[name] = {"terms": len(terms), "auc_apparent": app, "auc_cv": cvv}
    oof_store[name] = oof
    print(f"{name:<56}{len(terms):>6}{app:>10.4f}{cvv:>11.4f}")

BEST = "Independent conditions + oxygen + demographics"
print(f"\ngain from the parsimonious phenotype over pretreatment oxygen alone: "
      f"{AUC[BEST]['auc_cv'] - AUC['Pretreatment oxygen alone']['auc_cv']:+.4f} "
      f"cross-validated")

# ------------------------------------------------------------------ risk fifths
mm = m.copy()
mm["score"] = oof_store[BEST]
q = pd.qcut(mm.score, 5, labels=False)
risk = [{"fifth": int(i + 1), "n": int((q == i).sum()),
         "failed": int(mm.loc[q == i, "nonresp"].sum()),
         "failed_pct": round(100 * float(mm.loc[q == i, "nonresp"].mean()), 1)}
        for i in range(5)]
print("\nrisk fifths from the cross-validated phenotype score")
for r in risk:
    print(f"   fifth {r['fifth']}   n={r['n']:>4}   failed {r['failed']:>3} ({r['failed_pct']}%)")

# ------------------------------------------------------------------ cardiopulmonary dose-response
CPLAB = ["0", "1", "2", "3 or more"]
dose = {"labels": CPLAB,
        "n": [int((m.cp == i).sum()) for i in range(4)],
        "n_failed": [int(m.loc[m.cp == i, "nonresp"].sum()) for i in range(4)],
        "failed_pct": [round(100 * float(m.loc[m.cp == i, "nonresp"].mean()), 1)
                       for i in range(4)]}
assert sum(dose["n"]) == len(m), dose
_f, tr, ntr = logit(["cp", "AgeAtVisit", "male", "lbase"])
trend = {"model": ("logistic, failure ~ cardiopulmonary count capped at 3 + age + sex + "
                   "log(1 + pretreatment sleep T90)"),
         "or": tr["cp"]["or"], "lo": tr["cp"]["lo"], "hi": tr["cp"]["hi"],
         "p": tr["cp"]["p"], "n": ntr}
print(f"\ncardiopulmonary count: n {dose['n']}, failed % {dose['failed_pct']}, "
      f"trend OR {trend['or']} ({trend['lo']}-{trend['hi']})")

# ------------------------------------------------------------------ the panel-A table
rates = {}
for r in marg.itertuples():
    if r.neg:
        continue
    k = r.key
    rates[r.Condition] = {
        "n_with": r.n_with, "n_without": r.n_without,
        "n_failed_with": r.n_failed_with, "n_failed_without": r.n_failed_without,
        "fail_pct_with": r.fail_pct_with, "fail_pct_without": r.fail_pct_without,
        "marginal_or": r.marginal_or, "marginal_lo": r.marginal_lo,
        "marginal_hi": r.marginal_hi, "marginal_p": r.marginal_p,
        "marginal_q": float(r.marginal_q),
        "fdr_significant": bool(r.marginal_q < 0.05),
        "cp_adjusted_or": r.cp_adj_or, "cp_adjusted_lo": r.cp_adj_lo, "cp_adjusted_hi": r.cp_adj_hi,
        "cp_adjusted_p": r.cp_adj_p, "cp_adjusted_q": float(r.cp_adj_q),
        "adjusted_or": joint[f"{k}_prevalent"]["or"] if k in SIG else None,
        "adjusted_lo": joint[f"{k}_prevalent"]["lo"] if k in SIG else None,
        "adjusted_hi": joint[f"{k}_prevalent"]["hi"] if k in SIG else None,
        "adjusted_p": joint[f"{k}_prevalent"]["p"] if k in SIG else None,
        "independent": k in indep}
controls = {r.Condition: {"n_with": r.n_with, "n_without": r.n_without,
                          "fail_pct_with": r.fail_pct_with,
                          "fail_pct_without": r.fail_pct_without,
                          "marginal_or": r.marginal_or, "marginal_lo": r.marginal_lo,
                          "marginal_hi": r.marginal_hi, "marginal_p": r.marginal_p,
                          "marginal_q": float(r.marginal_q)}
            for r in marg.itertuples() if r.neg}

OUT = {
    "_built_by": "New_Figures/_NOT_THESE_workfiles/scripts/fits_phenotype_combined_v2.py",
    "_supersedes": "numbers/nonresponder_combined_v1.json, which was on the 1,996 frame",
    "_sample": ("Union of the split-night and paired treatment designs in "
                f"{os.path.basename(os.path.dirname(CPAP_STAGE_PARQUET))}/cpap_t90_by_stage.parquet, one record per patient, fu_valid == 1, "
                "at least 30 minutes of scored sleep on both nights, pre_sleep_t90 > 10. "
                "Identical to the frame numbers/run_treatment_v2.py uses for "
                "treatment_v2.json corrected_vs_not, so Figure 11, Figure 12, Figure 13, "
                "eFigure 10 and eTable 12 all stand on one sample."),
    "n": int(len(m)), "n_failed": int(m.nonresp.sum()),
    "n_corrected": int((1 - m.nonresp).sum()),
    "failure_pct": round(100 * float(m.nonresp.mean()), 1),
    "failure_rate_by_condition": rates,
    "negative_controls": controls,
    "n_conditions_tested": int((~marg.neg).sum()),
    "n_conditions_fdr_significant": len(SIG),
    "n_conditions_fdr_significant_cp_adjusted": int(((~marg.neg) & (marg.cp_adj_q < 0.05)).sum()),
    "cp_adjustment": f"marginal odds ratios also adjusted for any other prevalent condition among {CPX} (decision 6)",
    "n_conditions_independent_in_joint_model": len(indep),
    "joint_model": {"terms": JOINT_TERMS, "n": njoint,
                    "conditions": {label[k]: joint[f"{k}_prevalent"] for k in SIG},
                    "demographics": {t: joint[t] for t in DEMO + ["lbase"]}},
    "independent_conditions": {label[k]: joint[f"{k}_prevalent"] for k in indep},
    "auc": AUC,
    "risk_fifths_from_crossvalidated_score": risk,
    "dose_response": dose,
    "dose_response_trend": trend,
}
_pj, _pc = f"{OUT_DIR}/nonresponder_combined_v2.json", f"{OUT_DIR}/nonresponder_phenotype_v2.csv"
json.dump(OUT, open(_pj, "w"), indent=1)
marg.to_csv(_pc, index=False)
for _p in (_pj, _pc):
    sidecar(_p, __file__, extra_inputs=[f"{NUMBERS_DIR}/treatment_v2.json"], note="fits_phenotype_combined_v2 with the decision 6 covariate")
print("\nwrote numbers/nonresponder_combined_v2.json and "
      "numbers/nonresponder_phenotype_v2.csv")
