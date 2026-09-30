"""
Recompute the cohort block of MASTER_v2.json from the clean cohort.

MASTER_v2.json was assembled by hand and its cohort block still described the pre-QC set of
19 383, so Table 1, the figure legends and the abstract were all quoting a cohort that no
longer exists. Everything here is derived, so the block can be regenerated at any time and
never has to be typed again.

Writes numbers/MASTER_v2.json in place, touching only the "cohort" key.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../..")))  # repository root, where paths.py lives
import paths
import json
import sys

import pandas as pd

ROOT = paths.T90_ROOT
sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
from cohort_spec import apply_cohort

c = apply_cohort(pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet"))

t = c.spo2_pct_below_90
band = pd.cut(t, [-0.001, 1, 5, 10, 101], labels=["0-1%", "1-5%", "5-10%", ">10%"])
counts = band.value_counts().reindex(["0-1%", "1-5%", "5-10%", ">10%"])

male = (c.sex.astype(str).str.upper().str[0] == "M")
# follow-up in this design runs from the study to death or to the censoring date, which is
# what death_years already holds for every participant
fu = c.death_years

block = {
    "n": int(len(c)),
    "age_median": float(c.AgeAtVisit.median()),
    "age_q1": float(c.AgeAtVisit.quantile(0.25)),
    "age_q3": float(c.AgeAtVisit.quantile(0.75)),
    "male_pct": round(100 * float(male.mean()), 1),
    "t90_median": round(float(t.median()), 2),
    "t90_q1": round(float(t.quantile(0.25)), 2),
    "t90_q3": round(float(t.quantile(0.75)), 2),
    "t90_bands_pct": {k: round(100 * v / len(c), 1) for k, v in counts.items()},
    "t90_bands_n": {k: int(v) for k, v in counts.items()},
    "sites": {k: int(v) for k, v in c.site_id.value_counts().sort_index().items()},
    "fu_median": round(float(fu.median()), 2),
    "deaths": int(c.death_incident.sum()) if "death_incident" in c.columns else None,   # v7: Table 1 prints the deaths row
    "death_pct": round(100 * float(c.death_incident.mean()), 1) if "death_incident" in c.columns else None,
}

path = f"{paths.NUMBERS_DIR}/MASTER_v2.json"
m = json.load(open(path))
old = m.get("cohort", {})
m["cohort"] = block
json.dump(m, open(path, "w"), indent=1)

print(f"cohort block refreshed: {old.get('n'):,} -> {block['n']:,}")
for k in ("age_median", "male_pct", "t90_median", "fu_median"):
    print(f"  {k:14s} {old.get(k)} -> {block[k]}")
print(f"  bands {old.get('t90_bands_n')}\n     -> {block['t90_bands_n']}")
assert sum(block["t90_bands_n"].values()) == len(c), "bands must partition the cohort"
assert sum(block["sites"].values()) == len(c), "sites must partition the cohort"
print("\nbands and sites both partition the cohort exactly")
