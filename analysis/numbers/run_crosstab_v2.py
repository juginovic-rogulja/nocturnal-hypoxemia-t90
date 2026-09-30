"""
Authoritative AHI-by-oxygen cross-tabulation.

The frozen version in results_v2.json disagreed with a direct recount at the category
boundaries (mild 3572 vs 3565, moderate 4657 vs 4659, severe 8036 vs 8041). This file
recomputes the table from the cohort with explicit boundaries and asserts that the 4 groups
partition the cohort exactly, so Table 2 and Figure 4 read from a single verified source.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../..")))  # repository root, where paths.py lives
import paths
import json
import sys as _s; _s.path.insert(0, paths.ANALYSIS_DIR)
from cohort_spec import COHORT_N
import pandas as pd

ROOT = paths.T90_ROOT
b = pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet")
b = b[(b.fu_valid == 1) & b.spo2_pct_below_90.notna() & (b.oximetry_bad == 0)].copy()
N_COHORT_ALL = len(b); N_AHI_MISSING = int(b.AHI.isna().sum())   # v7: the index is undefined on nights with under 60 min of true sleep
b = b[b.AHI.notna()].copy()

CATS = [
    ("None (AHI <5)", b.AHI < 5),
    ("Mild (AHI 5 to <15)", (b.AHI >= 5) & (b.AHI < 15)),
    ("Moderate (AHI 15 to <30)", (b.AHI >= 15) & (b.AHI < 30)),
    ("Severe (AHI >=30)", b.AHI >= 30),
]
BANDS = [("0-1%", lambda t: t <= 1), ("1-5%", lambda t: (t > 1) & (t <= 5)),
         ("5-10%", lambda t: (t > 5) & (t <= 10)), (">10%", lambda t: t > 10)]

rows = []
for name, m in CATS:
    s = b[m]
    r = {"category": name, "n": int(len(s)), "t90_median": round(float(s.spo2_pct_below_90.median()), 2)}
    for bn, f in BANDS:
        sel = f(s.spo2_pct_below_90)
        r[f"n_{bn}"] = int(sel.sum())
        r[f"pct_{bn}"] = round(float(sel.mean() * 100), 1)
    rows.append(r)

# R30 (Alen 2026-09-08): three-group version for Fig 4c and 4d, no and mild apnea pooled (AHI under 15)
CATS3 = [
    ("No or mild (AHI <15)", b.AHI < 15),
    ("Moderate (AHI 15 to <30)", (b.AHI >= 15) & (b.AHI < 30)),
    ("Severe (AHI >=30)", b.AHI >= 30),
]
rows3 = []
for name, m in CATS3:
    s = b[m]
    r = {"category": name, "n": int(len(s)), "t90_median": round(float(s.spo2_pct_below_90.median()), 2)}
    for bn, f in BANDS:
        sel = f(s.spo2_pct_below_90)
        r[f"n_{bn}"] = int(sel.sum())
        r[f"pct_{bn}"] = round(float(sel.mean() * 100), 1)
    rows3.append(r)
assert sum(r["n"] for r in rows3) == len(b), "three groups do not partition the cohort"
assert rows3[0]["n"] == rows[0]["n"] + rows[1]["n"] and rows3[1]["n"] == rows[2]["n"] and rows3[2]["n"] == rows[3]["n"]
for r in rows3:
    assert sum(r[f"n_{bn}"] for bn, _ in BANDS) == r["n"], f"bands do not partition {r['category']}"

assert sum(r["n"] for r in rows) == len(b) == N_COHORT_ALL - N_AHI_MISSING, "categories do not partition the cohort"
assert N_COHORT_ALL == COHORT_N
for r in rows:
    assert sum(r[f"n_{bn}"] for bn, _ in BANDS) == r["n"], f"bands do not partition {r['category']}"

out = {"n": int(len(b)), "n_cohort": int(N_COHORT_ALL), "n_ahi_undefined": N_AHI_MISSING, "spearman": float(
    b.AHI.corr(b.spo2_pct_below_90, method="spearman")), "rows": rows, "rows3": rows3,  # v8.3 full precision (2026-09-16)
    "rows3_note": "three apnea groups (Alen 2026-09-08): no and mild pooled as AHI under 15, for Fig 4c, 4d and 4b's left half"}
json.dump(out, open(f"{paths.NUMBERS_DIR}/crosstab_v2.json", "w"), indent=1)

print(f"n = {out['n']}   Spearman = {out['spearman']}\n")
hdr = f"{'category':26s} {'n':>6}" + "".join(f"{b:>9}" for b, _ in BANDS)
print(hdr); print("-" * len(hdr))
for r in rows:
    print(f"{r['category']:26s} {r['n']:6d}" +
          "".join(f"{r[f'pct_{bn}']:8.1f}%" for bn, _ in BANDS))
print("\nthree-group version")
for r in rows3:
    print(f"{r['category']:26s} {r['n']:6d}" +
          "".join(f"{r[f'pct_{bn}']:8.1f}%" for bn, _ in BANDS))
print("\nwritten -> crosstab_v2.json (all partition checks passed)")
