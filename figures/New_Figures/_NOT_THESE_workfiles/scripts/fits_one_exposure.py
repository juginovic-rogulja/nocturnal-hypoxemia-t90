"""
One exposure per model, for the escalated eFigure 4.

Alen's objection to the old panel was that a bar labelled "severe sleep apnea" is not a bar
about sleep apnea, because severe apnea and low nocturnal oxygen travel together. His words
were "I want only one in one". This script therefore does two things.

  1. Every exposure is fitted in its OWN model. No bar is an estimate adjusted for another
     exposure that the reader cannot see.
  2. The three exposures he named are additionally restricted to participants whose
     oxygenation was normal, defined as time below 90% saturation of 1% of the night or less,
     which is the reference category used everywhere else in the paper. Inside that
     restriction a bar labelled "apnea" can only be carrying apnea.

The whole-cohort severe-apnea fit is kept alongside the restricted one on purpose. The gap
between them is the confounding Alen was pointing at, and it is worth showing rather than
merely correcting.

Specification is the paper's primary one: unpenalized Cox, stratified by site, natural cubic
age spline with one column dropped so the design is full rank, sex, sleep-study date as time
zero, oximetry_bad == 0. Exposures here are binary contrasts rather than per-SD, because the
question is how one clinically recognisable abnormality compares with another.

Writes numbers/one_exposure_v1.json.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../..")))  # repository root, where paths.py lives
import paths
import json
import sys
import warnings

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from patsy import dmatrix

warnings.filterwarnings("ignore")
ROOT = paths.FIGURE_ROOT
sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
from cohort_spec import apply_cohort, COHORT_N
from disease_definitions import DISEASES, NEGATIVE_CONTROLS, CIRCULAR

b = apply_cohort(pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet"))
assert len(b) == COHORT_N
print(f"cohort {len(b):,}")

b["male"] = (b.sex.astype(str).str.upper().str[0] == "M").astype(int)
sp = dmatrix("cr(a, df=4) - 1", {"a": b.AgeAtVisit.values}, return_type="dataframe").iloc[:, 1:]
for i in range(sp.shape[1]):
    b[f"age_s{i}"] = sp.iloc[:, i].values
ADJ = [f"age_s{i}" for i in range(sp.shape[1])] + ["male"]
assert len(ADJ) == 4, "one age-spline column must be dropped, leaving 3 plus sex"

b["t90"] = b.spo2_pct_below_90
NORMAL = b.t90 <= 1                      # the paper's reference category for oxygenation
print(f"normal oxygenation, time below 90% saturation <= 1% of the night: {NORMAL.sum():,} "
      f"({100 * NORMAL.mean():.1f}%)")

b["x_t90"] = (b.t90 > 10).astype(int)
b["x_short"] = (b.TST_min < 300).astype(float).where(b.TST_min.notna())            # v8.1: missing where TST is masked (split nights)
b["x_ahi"] = (b.AHI >= 30).astype(int)
b["x_se"] = (b.sleep_efficiency_pct < 85).astype(float).where(b.sleep_efficiency_pct.notna())   # v8.1
b["x_n3"] = (b.N3_pct < 5).astype(float).where(b.N3_pct.notna())                    # v8.1
for c in ["x_short", "x_se", "x_n3"]:
    b.loc[b[c].isna(), c] = np.nan

# label, column, restriction mask, printed definition
EXPOSURES = [
    ("Oxygen below 90% for more than 10% of the night", "x_t90", None,
     "spo2_pct_below_90 > 10, whole cohort"),
    ("Sleeping fewer than 5 hours", "x_short", None,
     "TST_min < 300, whole cohort"),
    ("Severe sleep apnea, apnea-hypopnea index 30 or more", "x_ahi", None,
     "AHI >= 30, whole cohort, oxygenation uncontrolled"),
    ("Severe sleep apnea, in normal oxygenation only", "x_ahi", "normal",
     "AHI >= 30 within spo2_pct_below_90 <= 1"),
    ("Sleep efficiency below 85%, in normal oxygenation only", "x_se", "normal",
     "sleep_efficiency_pct < 85 within spo2_pct_below_90 <= 1"),
    ("N3 sleep below 5% of the night, in normal oxygenation only", "x_n3", "normal",
     "N3_pct < 5 within spo2_pct_below_90 <= 1"),
    ("Sleeping fewer than 5 hours, in normal oxygenation only", "x_short", "normal",
     "TST_min < 300 within spo2_pct_below_90 <= 1"),
]

ALL = list(DISEASES.items()) + [("death", ("Death from any cause", False, [], []))]


def fit(frame, xcol, key):
    yc, ec, pc = f"{key}_years", f"{key}_incident", f"{key}_prevalent"
    if yc not in frame.columns:
        return None
    f = frame[(frame[pc] == 0) & frame[yc].notna() & (frame[yc] > 0)]
    d = pd.DataFrame({"T": f[yc], "E": f[ec].astype(int), "site": f.site_id,
                      "x": f[xcol], **{c: f[c] for c in ADJ}}).dropna()
    if d.E.sum() < 40 or d.x.nunique() < 2:
        return None
    # a site with no exposed participants contributes nothing and can break the stratified fit
    keep = d.groupby("site").x.transform("nunique") > 1
    d = d[keep]
    if d.E.sum() < 40:
        return None
    try:
        c = CoxPHFitter(penalizer=0.0).fit(d, "T", "E", strata=["site"])
    except Exception:
        return None
    r = c.summary.loc["x"]
    return {"hr": float(r["exp(coef)"]),  # v8.3 full precision (2026-09-16)
            "lo": float(r["exp(coef) lower 95%"]),
            "hi": float(r["exp(coef) upper 95%"]),
            "p": float(r["p"]), "events": int(d.E.sum()), "n": int(len(d)),
            "n_exposed": int(d.x.sum())}


OUT = {"cohort_n": int(len(b)),
       "normal_oxygen_n": int(NORMAL.sum()),
       "normal_oxygen_definition": "spo2_pct_below_90 <= 1% of the recording",
       "model": ("unpenalized Cox, stratified by site, natural cubic age spline with one "
                 "column dropped, sex; one exposure per model"),
       "exposures": {}}

for label, col, restrict, defn in EXPOSURES:
    frame = b[NORMAL] if restrict == "normal" else b
    n_exposed = int(frame[col].fillna(0).sum())
    rows = {}
    for key, (dlab, _neg, _9, _10) in ALL:
        r = fit(frame, col, key)
        if r is None:
            continue
        r["negative_control"] = key in NEGATIVE_CONTROLS
        r["circular"] = key in CIRCULAR
        r["key"] = key
        rows[dlab] = r
    OUT["exposures"][label] = {"column": col, "definition": defn,
                               "restricted_to_normal_oxygen": restrict == "normal",
                               "n_in_analysis_set": int(len(frame)),
                               "n_exposed": n_exposed,
                               "pct_exposed": round(100 * n_exposed / len(frame), 1),
                               "outcomes": rows}
    sig = sum(1 for v in rows.values() if v["lo"] > 1 and not v["circular"])
    print(f"{label:<58} exposed {n_exposed:>6,}  outcomes {len(rows):>3}  "
          f"significant and raised {sig:>3}")

# ---- the mutually adjusted pair, kept only so the report can say it does not matter -------
mut = {}
for key, (dlab, _n, _a, _c) in ALL:
    yc, ec, pc = f"{key}_years", f"{key}_incident", f"{key}_prevalent"
    if yc not in b.columns:
        continue
    f = b[(b[pc] == 0) & b[yc].notna() & (b[yc] > 0)]
    d = pd.DataFrame({"T": f[yc], "E": f[ec].astype(int), "site": f.site_id,
                      "x_t90": f.x_t90, "x_short": f.x_short,
                      **{c: f[c] for c in ADJ}}).dropna()
    if d.E.sum() < 60:
        continue
    try:
        c = CoxPHFitter(penalizer=0.0).fit(d, "T", "E", strata=["site"])
    except Exception:
        continue
    mut[dlab] = {x: float(np.exp(c.params_[x])) for x in ("x_t90", "x_short")}  # v8.3 full precision (2026-09-16)
OUT["mutually_adjusted_t90_and_short_sleep"] = mut

json.dump(OUT, open(f"{paths.NUMBERS_DIR}/one_exposure_v1.json", "w"), indent=1)
print(f"\nwrote numbers/one_exposure_v1.json")

# ---- consistency check against the frozen mutually adjusted values ------------------------
# results_v2.json is the unpenalized freeze and is the one to reproduce. MASTER_v2.json also
# carries a vs_sleep block, but it was written before the ridge penalty was removed and is
# shrunk, which matters because the old eFigure 4 was drawn from MASTER rather than results.
R2 = json.load(open(f"{paths.NUMBERS_DIR}/results_v2.json"))["vs_sleep_duration"]["outcomes"]
d = [(k, R2[k]["t90"]["hr"], mut[k]["x_t90"]) for k in R2 if k in mut]
diff = max(abs(a - c) for _, a, c in d)
print(f"mutually adjusted refit reproduces results_v2 vs_sleep_duration on {len(d)} "
      f"conditions, largest absolute difference {diff:.4f}")
assert diff <= 5e-4, "the refit does not reproduce the frozen unpenalized values within 5e-4, stop"  # v8.3: reference now full precision, compare within 5e-4 (2026-09-16)

M = json.load(open(f"{paths.NUMBERS_DIR}/MASTER_v2.json"))["vs_sleep"]
dm = [(abs(R2[k]["t90"]["hr"] - M[k]["t90"]), k) for k in M if k in R2]
dm.sort(reverse=True)
OUT["master_v2_vs_sleep_is_stale"] = {
    "n_compared": len(dm),
    "median_abs_difference": round(float(np.median([a for a, _ in dm])), 4),
    "largest": {"condition": dm[0][1], "difference": round(dm[0][0], 3),
                "master": M[dm[0][1]]["t90"], "unpenalized": R2[dm[0][1]]["t90"]["hr"]},
    "n_master_shrunk_toward_null": sum(1 for k in M if k in R2
                                       and M[k]["t90"] < R2[k]["t90"]["hr"]),
    "note": ("MASTER_v2.json vs_sleep is the ridge-penalized freeze. The published eFigure 4 "
             "reads it. results_v2.json vs_sleep_duration is the unpenalized one.")}
print(f"MASTER_v2 vs_sleep is stale: median absolute difference {OUT['master_v2_vs_sleep_is_stale']['median_abs_difference']}, "
      f"largest {dm[0][1]} {dm[0][0]:.3f}, shrunk in "
      f"{OUT['master_v2_vs_sleep_is_stale']['n_master_shrunk_toward_null']} of {len(dm)}")
json.dump(OUT, open(f"{paths.NUMBERS_DIR}/one_exposure_v1.json", "w"), indent=1)
own = OUT["exposures"]["Oxygen below 90% for more than 10% of the night"]["outcomes"]
sh = [(k, mut[k]["x_t90"], own[k]["hr"]) for k in mut if k in own]
md = float(np.median([abs(a - c) for _, a, c in sh]))
print(f"one-exposure vs mutually adjusted T90, median absolute difference in the hazard "
      f"ratio across {len(sh)} conditions: {md:.4f}")
