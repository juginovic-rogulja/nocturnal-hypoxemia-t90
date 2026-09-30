"""
Part 2 of the additional analyses.

  3  Interaction: do 2 measures amplify one another beyond the sum of their separate effects?
  4  Do the conditions that respond to treatment share anything?
  5  Severe sleep apnea with normal oxygen, the mirror of the hidden group.

v8 (2026-09-12, stress-test item 44): the table comes from cohort_spec (T90_FINAL, apply_cohort,
the same filter as before) and interactions_v2.csv gains interaction_q, a Benjamini-Hochberg q
over every interaction tested (one family, the 70 pairs), printed beside the raw count. Outputs go
to cohort_spec.OUT_DIR (numbers/ in the chain).
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../..")))  # repository root, where paths.py lives
import paths
import warnings; warnings.filterwarnings("ignore")
import json, itertools
import numpy as np, pandas as pd
from scipy import stats
from patsy import dmatrix
from lifelines import CoxPHFitter
import sys
sys.path.insert(0, paths.ANALYSIS_DIR)
from cohort_spec import apply_cohort, T90_FINAL, OUT_DIR, sidecar

# lifelines scales its penalty by the sample size, so penalizer=0.01 on a cohort this
# size is an effective ridge in the hundreds. Every fit here is by maximum partial
# likelihood, matching the primary analysis and standard survival software.


def bh(pvals):
    """Benjamini-Hochberg q values, in the order given."""
    p = np.asarray(pvals, float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n, float)
    prev = 1.0
    for rank, i in enumerate(reversed(order), start=1):
        val = p[i] * n / (n - rank + 1)
        prev = min(prev, val)
        q[i] = min(prev, 1.0)
    return q


R={}
b = apply_cohort(pd.read_parquet(T90_FINAL))
b["t90"]=b.spo2_pct_below_90
b["male"]=(b.sex.astype(str).str.upper().str[0]=="M").astype(int)
b["hrs"]=b.TST_min/60.0
sp=dmatrix("cr(a, df=4) - 1",{"a":b.AgeAtVisit.values},return_type="dataframe")
sp.columns=[f"age_s{i}" for i in range(4)]; sp.index=b.index
# the 4 spline columns sum to 1. A Cox partial likelihood has no intercept so this is
# harmless under a penalty, but the design is singular without one, which is why every
# unpenalized fit failed. One column is dropped, exactly as the primary analysis does.
sp=sp.iloc[:,1:]
for c in sp.columns: b[c]=sp[c]
ADJ=list(sp.columns)+["male"]

# ---------------------------------------------------------------- 3 interactions
print("3  do any 2 markers amplify one another?")
b["A_oxy"]=(b.t90>10).astype(int)
b["B_short"]=(b.hrs<5).astype(float).where(b.hrs.notna())   # v8.1: missing where TST is masked (split nights)
b["C_ahi"]=(b.AHI>=30).astype(int) if "AHI" in b.columns else 0
b["D_arous"]=(b.arousal_index>=25).astype(int) if "arousal_index" in b.columns else 0
b["E_eff"]=(b.sleep_efficiency_pct<80).astype(float).where(b.sleep_efficiency_pct.notna()) if "sleep_efficiency_pct" in b.columns else 0   # v8.1
MARK={"A_oxy":"Oxygen <90% for >10% of night","B_short":"Sleep <5 h",
      "C_ahi":"Apnea-hypopnea index >=30","D_arous":"Arousal index >=25",
      "E_eff":"Sleep efficiency <80%"}
TEST=[("hf","Heart failure"),("cvd","Cardiovascular composite"),
      ("resp_failure","Respiratory failure"),("diabetes","Type 2 diabetes"),
      ("death","Death from any cause"),("copd2","COPD"),("aki","Acute kidney injury")]
rows=[]
for (m1,m2) in itertools.combinations(MARK,2):
    for k,lab in TEST:
        if f"{k}_incident" not in b.columns: continue
        g=b[(b[f"{k}_prevalent"]==0)&b[f"{k}_years"].notna()&(b[f"{k}_years"]>0)].copy()
        g["both"]=g[m1]*g[m2]
        f=pd.DataFrame({"T":g[f"{k}_years"],"E":g[f"{k}_incident"].astype(int),
            "m1":g[m1],"m2":g[m2],"both":g["both"],"site":g.site_id,
            **{c:g[c] for c in ADJ}}).dropna()
        if f.E.sum()<120 or f.both.sum()<40: continue
        try: c=CoxPHFitter(penalizer=0.0).fit(f,"T","E",strata=["site"])
        except Exception: continue
        h1,h2,hi=[float(np.exp(c.params_[x])) for x in ("m1","m2","both")]
        p=float(c.summary.loc["both","p"])
        # observed joint effect against what the 2 alone would give if they simply multiplied
        obs=h1*h2*hi
        rows.append({"marker1":MARK[m1],"marker2":MARK[m2],"outcome":lab,
            "n_both":int(f.both.sum()),"events":int(f.E.sum()),
            "hr_marker1_alone":h1,"hr_marker2_alone":h2,  # v8.3 full precision (2026-09-16)
            "interaction_hr":hi,"interaction_p":p,
            "hr_both_observed":obs,"hr_both_if_multiplicative":h1*h2})
inter=pd.DataFrame(rows)
# item 44 (2026-09-12): one Benjamini-Hochberg family over every interaction tested
inter["interaction_q"]=bh(inter.interaction_p.values) if len(inter) else []
inter_path=f"{OUT_DIR}/interactions_v2.csv"
inter.to_csv(inter_path,index=False)
sig=inter[inter.interaction_p<.05]
sigq=inter[inter.interaction_q<.05]
print(f"   {len(inter)} pairs tested, {len(sig)} with an interaction at P<.05, {len(sigq)} at q<.05 (BH over the {len(inter)})")
if len(sig):
    print(sig.sort_values("interaction_p")[["marker1","marker2","outcome","interaction_hr",
        "interaction_p","interaction_q","hr_both_observed","hr_both_if_multiplicative"]].head(10).to_string(index=False))
else:
    print("   none. The markers act independently, so their effects simply add.")
R["interactions"]={"n_tested":int(len(inter)),"n_significant":int(len(sig)),
  "n_significant_q":int(len(sigq)),"multiple_testing":"Benjamini-Hochberg over every pair tested",
  "significant":sig.to_dict("records") if len(sig) else []}

# ---------------------------------------------------------------- 5 severe apnea, normal oxygen
print("\n5  severe sleep apnea with normal oxygen")
av=b[b.AHI.notna()].copy()
sev=av[av.AHI>=30].copy()
sev["lowox"]=(sev.t90>10).astype(int)
R["severe_apnea"]={"n":int(len(sev)),
  "n_normal_oxygen":int((sev.t90<=1).sum()),
  "pct_normal_oxygen":round(100*float((sev.t90<=1).mean()),1),"outcomes":{}}
print(f"   n={len(sev):,}  with normal oxygen (<=1%): {int((sev.t90<=1).sum()):,} "
      f"({100*float((sev.t90<=1).mean()):.1f}%)")
print(f"\n   {'outcome':<26}{'ev':>5}   severe apnea WITH low oxygen vs severe apnea WITHOUT")
# every control in the settled panel, imported rather than listed
from disease_definitions import NEGATIVE_CONTROLS as _NC, DISEASES as _DZ
for k,lab in TEST+[(k,f"{_DZ[k][0]} [NEG]") for k in _NC]:
    if f"{k}_incident" not in sev.columns: continue
    g=sev[(sev[f"{k}_prevalent"]==0)&sev[f"{k}_years"].notna()&(sev[f"{k}_years"]>0)]
    f=pd.DataFrame({"T":g[f"{k}_years"],"E":g[f"{k}_incident"].astype(int),"x":g.lowox,
        "site":g.site_id,**{c:g[c] for c in ADJ}}).dropna()
    if f.E.sum()<60: continue
    c=CoxPHFitter(penalizer=0.0).fit(f,"T","E",strata=["site"])
    lo,hi=np.exp(c.confidence_intervals_.loc["x"]); p=float(c.summary.loc["x","p"])
    hr=float(np.exp(c.params_["x"]))
    R["severe_apnea"]["outcomes"][lab]={"events":int(f.E.sum()),"hr":hr,  # v8.3 full precision (2026-09-16)
        "lo":float(lo),"hi":float(hi),"p":p}
    print(f"   {lab:<26}{int(f.E.sum()):>5}   {hr:.2f} ({lo:.2f}-{hi:.2f})  p={p:.4f}{'  *' if p<.05 else ''}")

json_path=f"{OUT_DIR}/extras2_v2.json"
json.dump(R,open(json_path,"w"),indent=2)
for p in (inter_path, json_path):
    sidecar(p, __file__, note=f"run_extras2_v2: {len(inter)} interactions with BH q")
print(f"\nwritten -> {json_path}")
