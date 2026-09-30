"""
Freeze the cardiopulmonary dose-response for failure to normalize oxygen.

Why this exists. The four percentages in the manuscript (20.9, 29.7, 38.0, 60.8) had no frozen
source. The consistency audit (audit/STRESS_5_CONSISTENCY.md, finding A16) traced them to a
hardcoded literal in figures/make_new_efigs.py. That script now recomputes them inline at draw
time, which fixes the drift but still leaves the numbers unfrozen and unavailable to anything
else. The per-condition odds ratio existed nowhere at all.

This script rebuilds the same sample used for numbers/nonresponder_phenotype_audited.csv,
reproduces the four percentages exactly, adds the trend odds ratio, and writes
numbers/nonresponder_doseresponse.json.

Run:  python3 figures_alen/scripts/freeze_doseresponse.py
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import glob
import json

import numpy as np
import pandas as pd
import statsmodels.api as sm

ROOT = paths.FIGURE_ROOT
SRC = f"{paths.X4_DIR}/cpap_stage"   # v8.2: re-extracted, collapse-masked shards (X4_groupP)
OUT = f"{paths.NUMBERS_DIR}/nonresponder_doseresponse.json"

CP = ["resp_failure", "copd2", "obesity_hypovent", "hf", "pulm_htn", "asthma"]

# ---- the sample, identical in construction to nonresponder_phenotype_audited.csv -----------
o = pd.concat([pd.read_csv(f, low_memory=False) for f in sorted(glob.glob(f"{SRC}/out*.csv"))],
              ignore_index=True).drop_duplicates(["BIDSFolder", "SessionID"])
o = o[(o.status == "ok") & (o.arm == "A")
      & o.pre_sleep_t90.notna() & o.post_sleep_t90.notna()]
s = o[o.pre_sleep_t90 > 10].drop_duplicates("BDSPPatientID").copy()
# same oximetry quality rule as the cohort file: a nadir above the mean is not a measurement
s = s[~((s.pre_sleep_nadir > s.pre_sleep_mean) | (s.post_sleep_nadir > s.post_sleep_mean))]
s["nonresp"] = (s.post_sleep_t90 > 10).astype(int)

# 2026-08-07: prevalent flags now come from data_frozen/t90_final.parquet, not
# data_frozen/t90_outcomes_rebuilt.parquet. The two agree on 49 of the 51 shared conditions and
# disagree on obesity hypoventilation (ICD-9 278.01/278.03 correction) and falls (E888/E884
# prefix-collision fix), and t90_final is the corrected one in both cases. t90_outcomes_rebuilt
# also has no column for contact dermatitis or hemorrhoids, so two of the five settled negative
# controls could never be fitted from it.
od = pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet").drop_duplicates("BDSPPatientID")
cols = [f"{k}_prevalent" for k in CP]
m = s.merge(od[["BDSPPatientID"] + cols], on="BDSPPatientID", how="left")
m["cp_raw"] = sum(m[c].fillna(0) for c in cols)
m["cp"] = m["cp_raw"].clip(upper=3)
m["male"] = (m.SexDSC.astype(str).str.upper().str[0] == "M").astype(int)
m["lbase"] = np.log1p(m.pre_sleep_t90)

# ---- the four bars -------------------------------------------------------------------------
nn = [int((m.cp == k).sum()) for k in range(4)]
val = [round(float(m.loc[m.cp == k, "nonresp"].mean()) * 100, 1) for k in range(4)]
nf = [int(m.loc[m.cp == k, "nonresp"].sum()) for k in range(4)]


def logit(terms):
    X = sm.add_constant(m[terms].astype(float))
    d = pd.concat([X, m.nonresp], axis=1).dropna()
    r = sm.Logit(d.nonresp, d[X.columns]).fit(disp=0)
    b, ci = r.params[terms[0]], r.conf_int().loc[terms[0]]
    return {"or": round(float(np.exp(b)), 3), "lo": round(float(np.exp(ci[0])), 3),
            "hi": round(float(np.exp(ci[1])), 3), "p": float(r.pvalues[terms[0]]),
            "n": int(len(d))}


primary = logit(["cp", "AgeAtVisit", "male", "lbase"])
R = {
    "_created": "2026-08-07",
    "_why": "The cardiopulmonary dose-response for failure to normalize oxygen had no frozen "
            "source. STRESS_5_CONSISTENCY.md finding A16 traced it to a hardcoded literal in "
            "figures/make_new_efigs.py. This file freezes the recomputed values and adds the "
            "per-condition odds ratio, which existed nowhere.",
    "_built_by": "figures_alen/scripts/freeze_doseresponse.py",
    "_sample": "Same construction as numbers/nonresponder_phenotype_audited.csv: arm A "
               "split-night records with status ok and non-missing pre and post sleep T90, "
               "restricted to pre_sleep_t90 > 10, one record per patient, oximetry quality rule "
               "(nadir must not exceed mean) applied to both nights.",
    "n_total": int(len(m)),
    "n_corrected": int((1 - m.nonresp).sum()),
    "n_not_corrected": int(m.nonresp.sum()),
    "failure_pct_overall": round(float(m.nonresp.mean()) * 100, 1),
    "cardiopulmonary_conditions": CP,
    "dose_response": {"labels": ["0", "1", "2", "3 or more"], "n": nn,
                      "failed_pct": val, "n_failed": nf},
    "count_uncapped_distribution": {str(k): int(v) for k, v
                                    in m.cp_raw.value_counts().sort_index().items()},
    "trend": {
        "primary": {"model": "logistic, failure ~ cardiopulmonary count capped at 3 + age + sex "
                             "+ log(1 + pretreatment sleep T90)", **primary},
        "alternatives": {
            "capped_unadjusted": logit(["cp"]),
            "capped_age_sex": logit(["cp", "AgeAtVisit", "male"]),
            "uncapped_unadjusted": logit(["cp_raw"]),
            "uncapped_age_sex": logit(["cp_raw", "AgeAtVisit", "male"]),
            "uncapped_age_sex_lbase": logit(["cp_raw", "AgeAtVisit", "male", "lbase"]),
        },
        "_note": "Two models round to 1.57. The primary one is reported because it matches both "
                 "the capped count plotted on the x axis and the adjustment set used for the "
                 "per-condition odds ratios in nonresponder_phenotype_audited.csv.",
    },
}
json.dump(R, open(OUT, "w"), indent=2)
print(f"n={R['n_total']}  corrected={R['n_corrected']}  not={R['n_not_corrected']}")
print(f"dose-response {val} on n={nn}")
print(f"OR per condition {primary['or']} ({primary['lo']}-{primary['hi']})")
print(f"wrote {OUT}")
