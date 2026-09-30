"""
Per-1-SD association of nocturnal oxygen with death from any cause in the primary cohort.

bdsp_diseases_v2.csv covers the 51 conditions but not mortality, and Figure 6 compares the
primary cohort against both community cohorts for mortality. Same specification as every other
primary estimate: site-stratified Cox, sex, and a natural cubic spline on age with 4 df, with
the exposure transformed to a rank-based inverse normal score within hospital.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../..")))  # repository root, where paths.py lives
import paths
import warnings; warnings.filterwarnings("ignore")
import json
import numpy as np
import sys as _s; _s.path.insert(0, paths.ANALYSIS_DIR)
from cohort_spec import COHORT_N
import pandas as pd
from patsy import dmatrix
from scipy import stats
from lifelines import CoxPHFitter

ROOT = paths.T90_ROOT
b = pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet")
b = b[(b.fu_valid == 1) & b.spo2_pct_below_90.notna() & (b.oximetry_bad == 0)].copy()
assert len(b) == COHORT_N

b["male"] = (b.sex.astype(str).str.upper().str[0] == "M").astype(int)
sp = dmatrix("cr(a, df=4) - 1", {"a": b.AgeAtVisit.values}, return_type="dataframe")
sp.columns = [f"age_s{i}" for i in range(4)]; sp.index = b.index
for c in sp.columns:
    b[c] = sp[c]

# rank-based inverse normal within hospital, as everywhere else in the paper
def rint(x):
    r = stats.rankdata(x)
    return stats.norm.ppf((r - 0.375) / (len(r) + 0.25))
b["z"] = b.groupby("site_id").spo2_pct_below_90.transform(rint)

f = b[(b.death_prevalent == 0) & b.death_years.notna() & (b.death_years > 0)].copy()
f["E"] = f.death_incident.astype(int)
f["T"] = f.death_years

cols = ["z", "male"] + [f"age_s{i}" for i in range(4)]
d = pd.DataFrame({"T": f["T"], "E": f["E"], "site": f.site_id,
                  **{c: f[c] for c in cols}}).dropna()
m = CoxPHFitter(penalizer=0.01).fit(d, "T", "E", strata=["site"])
s = m.summary.loc["z"]
out = {"events": int(d.E.sum()), "n": int(len(d)),
       "hr": round(float(s["exp(coef)"]), 3),
       "lo": round(float(s["exp(coef) lower 95%"]), 3),
       "hi": round(float(s["exp(coef) upper 95%"]), 3),
       "p": float(s["p"])}
json.dump(out, open(f"{paths.NUMBERS_DIR}/mortality_v2.json", "w"), indent=1)
print(f"Death from any cause, per 1 SD: HR {out['hr']} [{out['lo']}-{out['hi']}]  "
      f"deaths {out['events']:,} of {out['n']:,}  p={out['p']:.2e}")
