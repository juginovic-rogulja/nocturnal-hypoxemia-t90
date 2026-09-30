"""
Fill the numeric gaps the second round of comments opened up.

Nothing here re-derives an existing number. It only computes quantities the manuscript is now
asked to state and does not currently have: the age range, total sleep time in each cohort, how
strictly the graded pattern holds, and the treatment comparison expressed as hazard ratios.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../..")))  # repository root, where paths.py lives
import paths
import json
import numpy as np
import sys as _s; _s.path.insert(0, paths.ANALYSIS_DIR)
from cohort_spec import COHORT_N
import pandas as pd

ROOT = paths.T90_ROOT
out = {}

b = pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet")
b = b[(b.fu_valid == 1) & b.spo2_pct_below_90.notna() & (b.oximetry_bad == 0)].copy()
print("cohort rows:", len(b), " (must be 19383)")
assert len(b) == COHORT_N, "cohort does not match the frozen count"

# ---- 1. age range, asked for in Participants ----------------------------------------
a = b["AgeAtVisit"].dropna()
out["age"] = {"min": float(a.min()), "max": float(a.max()), "n_ge_90": int((a >= 90).sum())}
print("age:", out["age"])

# ---- 2. total sleep time, asked for in Participants across cohorts -------------------
t = b["TST_min"].dropna()
out["tst_bdsp"] = {
    "n": int(t.notna().sum()),
    "median_h": round(float(t.median()) / 60, 1),
    "q1_h": round(float(t.quantile(.25)) / 60, 1),
    "q3_h": round(float(t.quantile(.75)) / 60, 1),
    "pct_under_5h": round(float((t < 300).mean() * 100), 1),
    "pct_under_6h": round(float((t < 360).mean() * 100), 1),
    "pct_under_7h": round(float((t < 420).mean() * 100), 1),
}
print("TST:", out["tst_bdsp"])

# ---- 3. how strictly the graded pattern holds ----------------------------------------
g = json.load(open(f"{paths.NUMBERS_DIR}/results_v2.json"))["graded"]["outcomes"]
mono, strict, neg_mono = [], [], []
for cond, v in g.items():
    try:
        hrs = [v[k]["hr"] for k in ["1-5%", "5-10%", ">10%"]]
    except KeyError:
        continue
    if hrs[0] <= hrs[1] <= hrs[2] and hrs[2] > 1:
        (neg_mono if v["negative_control"] else mono).append((cond, hrs))
        if v[">10%"]["lo"] > 1 and not v["negative_control"]:
            strict.append((cond, hrs, v[">10%"]["hr"]))
out["graded"] = {
    "n_evaluated": len(g),
    "n_monotone": len(mono),
    "n_monotone_top_significant": len(strict),
    "n_negative_controls_monotone": len(neg_mono),
    "monotone_conditions": sorted(c for c, _ in mono),
    "top_band_hr_range": [round(min(x[2] for x in strict), 2),
                          round(max(x[2] for x in strict), 2)],
}
print(f"\ngraded: {len(mono)} of {len(g)} rise across all 4 categories, "
      f"{len(strict)} with the top category significant, "
      f"{len(neg_mono)} of the negative controls")

# ---- 4. treatment comparison as hazard ratios, asked for beside the percentages -------
c = json.load(open(f"{paths.NUMBERS_DIR}/treatment_v2.json"))["corrected_vs_not"]
sig = {k: v for k, v in c["outcomes"].items()
       if v["p"] < 0.05 and not v.get("negative_control")}
rows = sorted(sig.items(), key=lambda x: -x[1]["hr"])
out["treatment_sig"] = [
    {"condition": k, "hr": round(v["hr"], 2), "lo": round(v["lo"], 2), "hi": round(v["hi"], 2),
     "pct_lower_if_restored": round((1 - 1 / v["hr"]) * 100)}
    for k, v in rows
]
negs = [v["hr"] for v in c["outcomes"].values() if v.get("negative_control")]
out["treatment_negcontrol"] = {
    "n": len(negs), "geomean": round(float(np.exp(np.mean(np.log(negs)))), 2),
    "range": [round(min(negs), 2), round(max(negs), 2)],
}
out["treatment_median_pct"] = round(float(np.median(
    [r["pct_lower_if_restored"] for r in out["treatment_sig"]])))
print(f"\ntreatment: {len(rows)} conditions, median {out['treatment_median_pct']}% lower")
for r in out["treatment_sig"]:
    print(f"  {r['condition']:28s} HR {r['hr']:.2f} [{r['lo']:.2f}-{r['hi']:.2f}]"
          f"   {r['pct_lower_if_restored']}% lower")
print("negative controls:", out["treatment_negcontrol"])

json.dump(out, open(f"{paths.NUMBERS_DIR}/gaps_v2.json", "w"), indent=1)
print("\nwritten -> gaps_v2.json")
