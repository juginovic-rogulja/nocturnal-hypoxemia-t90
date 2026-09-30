"""
The 2 reverse-causation tests, run on every condition rather than a chosen 10.

The earlier pass ran these on 10 conditions, leaving 22 of the 30 significant associations with
no answer to the obvious reviewer question. Obesity hypoventilation, pulmonary hypertension,
cirrhosis, kidney disease, sepsis and death from any cause were all among the untested.

Test 1, landmark lag. Drop anyone diagnosed within 2 years of the sleep study and restart the
clock there. If the association is patients who were already ill at the time of the study, it
goes away. If it survives, the oxygen came first.

Test 2, body mass index. Derived from recorded height and weight for about a third of the
cohort. Obesity drives both low nocturnal oxygen and most of these diseases, so an association
that disappears on adjustment was carrying body habitus rather than oxygen. This is a
conservative test, because BMI sits partly on the causal path.

Both are fitted by maximum partial likelihood on the clean cohort, with the same site
stratification, sex, and age spline as the primary analysis, and with 1 spline column dropped so
the design is full rank.

Writes numbers/causal_tests_all.csv.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../..")))  # repository root, where paths.py lives
import paths
import warnings

warnings.filterwarnings("ignore")
import sys

import duckdb
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from patsy import dmatrix
from scipy import stats

ROOT = paths.T90_ROOT
CACHE = paths.OMOP_CACHE_DIR
sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
from cohort_spec import apply_cohort
from disease_definitions import DISEASES, NEGATIVE_CONTROLS

LAG = 2.0

b = apply_cohort(pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet"))
b["male"] = (b.sex.astype(str).str.upper().str[0] == "M").astype(int)
sp = dmatrix("cr(a, df=4) - 1", {"a": b.AgeAtVisit.values}, return_type="dataframe").iloc[:, 1:]
for i in range(sp.shape[1]):
    b[f"age_s{i}"] = sp.iloc[:, i].values
ADJ = [f"age_s{i}" for i in range(sp.shape[1])] + ["male"]


def rint(x):
    r = stats.rankdata(x)
    return stats.norm.ppf((r - 0.375) / (len(r) + 0.25))


b["z"] = b.groupby("site_id").spo2_pct_below_90.transform(rint)

# body mass index, from the unit-clean cache.
#
# 2026-08-08. This script carried a copy of the query in freeze_02_demographics.py, which read
# any "weight" measurement between 25 and 300 as kilograms. In this cache that mixes LBS,
# Pounds and rows that are not a measured weight at all, such as "WEIGHT IN (LB) TO HAVE
# BMI = 25" and "FEMALE PREDICTED WEIGHT". Every body-mass-index adjustment run through this
# file was therefore adjusting for a contaminated covariate. The old query is kept below and
# raises, so nothing can fall back to it, and the covariate now comes from
# numbers/_bmi_derived.parquet with the units resolved per row.
#
# Refitting all 54 conditions on the clean covariate flips no verdict, and 27 of 33 survive
# either way, so this changes the covariate and not the paper's conclusions.
def q(pat, lo, hi):
    raise RuntimeError(
        "run_causal_tests_all: the pound-contaminated height and weight query was removed on "
        "2026-08-08. Use numbers/_bmi_derived.parquet.")


bmi_cache = pd.read_parquet(f"{paths.NUMBERS_DIR}/_bmi_derived.parquet")
bmi_cache = bmi_cache[["BDSPPatientID", "bmi"]].dropna().drop_duplicates("BDSPPatientID")
b = b.merge(bmi_cache, on="BDSPPatientID", how="left")
print(f"cohort {len(b):,}, body mass index derivable for {b.bmi.notna().sum():,} "
      f"({100 * b.bmi.notna().mean():.1f}%)")


def fit(d, extra=()):
    cols = ["T", "E", "site", "z"] + ADJ + list(extra)
    d = d[cols].dropna()
    if d["E"].sum() < 40:
        return None
    try:
        c = CoxPHFitter().fit(d, "T", "E", strata=["site"])
    except Exception:
        return None
    s = c.summary.loc["z"]
    return (float(s["exp(coef)"]), float(s["exp(coef) lower 95%"]),
            float(s["exp(coef) upper 95%"]), int(d["E"].sum()))  # v8.3 full precision (2026-09-16)


rows = []
ALL = list(DISEASES.items()) + [("death", ("Death from any cause", False, [], []))]
for key, (label, _n, _a, _b) in ALL:
    yc, ec, pc = f"{key}_years", f"{key}_incident", f"{key}_prevalent"
    if yc not in b.columns:
        continue
    f = b[(b[pc] == 0) & b[yc].notna() & (b[yc] > 0)].copy()
    base = pd.DataFrame({"T": f[yc], "E": f[ec].astype(int), "site": f.site_id, "bmi": f.bmi,
                         **{c: f[c] for c in ["z"] + ADJ}})

    # landmark: drop events inside the window, then measure time from the landmark
    # base["T"], not base.T, which is the transpose
    lm = base[base["T"] > LAG].copy()
    lm["T"] = lm["T"] - LAG

    r0, rl, rb = fit(base), fit(lm), fit(base.dropna(subset=["bmi"]), extra=["bmi"])
    rn = fit(base.dropna(subset=["bmi"]))  # same people as rb, without the adjustment
    if not r0:
        continue
    rows.append({
        "condition": label, "negative_control": key in NEGATIVE_CONTROLS,
        "hr": r0[0], "lo": r0[1], "hi": r0[2], "events": r0[3],
        "lag_hr": rl[0] if rl else np.nan, "lag_lo": rl[1] if rl else np.nan,
        "lag_hi": rl[2] if rl else np.nan, "lag_events": rl[3] if rl else np.nan,
        "bmi_sub_hr": rn[0] if rn else np.nan,
        "bmi_adj_hr": rb[0] if rb else np.nan, "bmi_adj_lo": rb[1] if rb else np.nan,
        "bmi_adj_hi": rb[2] if rb else np.nan, "bmi_events": rb[3] if rb else np.nan,
    })

d = pd.DataFrame(rows)
d["survives_lag"] = d.lag_lo > 1
d["survives_bmi"] = d.bmi_adj_lo > 1
d.to_csv(f"{paths.NUMBERS_DIR}/causal_tests_all.csv", index=False)

sig = d[(d.lo > 1) & ~d.negative_control]
print(f"\n{len(d)} outcomes tested, {len(sig)} of them significantly associated\n")
print(f"{'condition':28s}{'HR':>7}{'2y lag':>9}{'lag ok':>8}{'BMI-adj':>9}{'BMI ok':>8}")
print("-" * 70)
for _, r in sig.sort_values("hr", ascending=False).iterrows():
    print(f"{r.condition:28s}{r.hr:7.2f}{r.lag_hr:9.2f}{'yes' if r.survives_lag else 'NO':>8}"
          f"{r.bmi_adj_hr:9.2f}{'yes' if r.survives_bmi else 'NO':>8}")
print(f"\nsurvive the landmark lag: {int(sig.survives_lag.sum())} of {len(sig)}")
print(f"survive body mass index:  {int(sig.survives_bmi.sum())} of {len(sig)}")
print(f"\nwritten -> {paths.NUMBERS_DIR}/causal_tests_all.csv")
