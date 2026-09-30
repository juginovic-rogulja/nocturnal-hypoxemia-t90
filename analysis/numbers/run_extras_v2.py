"""
Additional analyses requested on the second read.

  1  Where do AHI variants rank, with and without a desaturation requirement? If the versions
     that require desaturation rank higher, the index carries information only insofar as it
     tracks oxygen.
  2  Where would short sleep rank if entered as the clinical threshold rather than as a
     continuous measure?
  3  Do any 2 measures amplify one another beyond the sum of their separate effects?
  4  Do the conditions that respond share anything, or is the pattern arbitrary?
  5  What happens to patients with severe sleep apnea but normal oxygen?
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../..")))  # repository root, where paths.py lives
import paths
import warnings; warnings.filterwarnings("ignore")
import json, itertools
import numpy as np, pandas as pd
from scipy import stats
from patsy import dmatrix
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index

# lifelines scales its penalty by the sample size, so penalizer=0.01 on a cohort this
# size is an effective ridge in the hundreds. Every fit here is by maximum partial
# likelihood, matching the primary analysis and standard survival software.

R={}
b=pd.read_parquet(f'{paths.TABLES_DIR}/t90_final.parquet')
b = b[(b.fu_valid == 1) & b.spo2_pct_below_90.notna() & (b.oximetry_bad == 0)].copy()
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
OUT=[k[:-9] for k in b.columns if k.endswith("_incident")]
OUT=[k for k in OUT if int(b[f"{k}_incident"].sum())>=150]

# ---------------------------------------------------------------- 1 AHI variants
print("1  AHI variants: does the index work only through oxygen?")
raw=pd.read_csv(f"{paths.TABLES_DIR}/"
                "master_cohort.csv",low_memory=False)
cand=[c for c in raw.columns if c.lower().startswith(("ahi","rdi","oai","cai"))]
b2=b.merge(raw[["BDSPPatientID"]+cand],on="BDSPPatientID",how="left",suffixes=("","_r"))
def rint(df,col):
    z=pd.Series(np.nan,index=df.index)
    for s,idx in df.groupby("site_id").groups.items():
        v=df.loc[idx,col]; ok=v.notna()
        if ok.sum()<20: continue
        r=stats.rankdata(v[ok],method="average")
        z.loc[v[ok].index]=stats.norm.ppf((r-0.375)/(ok.sum()+0.25))
    return z
def heldout(df,feat):
    cs=[]
    for out in OUT:
        g=df[(df[f"{out}_prevalent"]==0)&(df[f"{out}_years"]>0)&df[f"{out}_years"].notna()]
        for tr,te in (("I0002","I0006"),("I0006","I0002")):
            A=g[g.site_id==tr][ADJ+([feat] if feat else [])+[f"{out}_years",f"{out}_incident"]].dropna()
            B=g[g.site_id==te][ADJ+([feat] if feat else [])+[f"{out}_years",f"{out}_incident"]].dropna()
            if A[f"{out}_incident"].sum()<30 or B[f"{out}_incident"].sum()<30: continue
            try:
                c=CoxPHFitter(penalizer=0.0).fit(A,duration_col=f"{out}_years",event_col=f"{out}_incident")
                cs.append(concordance_index(B[f"{out}_years"],-c.predict_partial_hazard(B),B[f"{out}_incident"]))
            except Exception: pass
    return float(np.mean(cs)) if cs else np.nan
base=heldout(b2,None)
rows=[]
for c in cand:
    if b2[c].notna().mean()<0.6 or b2[c].nunique()<10: continue
    b2[c+"__z"]=rint(b2,c)
    g=heldout(b2,c+"__z")
    if np.isfinite(g): rows.append({"variant":c,"gain":round(g-base,5)})
ahi=pd.DataFrame(rows).sort_values("gain",ascending=False)
R["ahi_variants"]={"baseline_c":round(base,4),"rows":ahi.to_dict("records")}
print(ahi.head(12).to_string(index=False))

# ---------------------------------------------------------------- 2 short sleep thresholds
print("\n2  short sleep as a clinical threshold")
for thr,name in [(5,"lt5h"),(4,"lt4h"),(6,"lt6h"),(7,"lt7h")]:
    b2[name]=(b2.hrs<thr).astype(float).where(b2.hrs.notna())   # v8.1: missing where TST is masked (split nights)
    g=heldout(b2,name)
    R.setdefault("short_sleep",{})[f"<{thr}h"]={"gain":round(g-base,5),
        "pct":round(100*float((b2.hrs.dropna()<thr).mean()),1)}
    print(f"   sleep <{thr}h  gain {g-base:+.5f}   present in {100*float((b2.hrs.dropna()<thr).mean()):.1f}%")
t90g=heldout(b2.assign(t90_z=rint(b2,"t90")),"t90_z")
R["short_sleep"]["t90_continuous"]={"gain":round(t90g-base,5)}
b2["t90_gt10"]=(b2.t90>10).astype(float)
g=heldout(b2,"t90_gt10")
R["short_sleep"]["t90_gt10"]={"gain":round(g-base,5)}
print(f"   oxygen <90% for >10%  gain {g-base:+.5f}")

# the output goes to the declared output folder, never to the working directory, and carries a
# provenance sidecar (sha256 of itself, of this script and of the four tables) as every other
# output of this folder does.
_out = f"{paths.NUMBERS_DIR}/extras_v2.json"
_os.makedirs(paths.NUMBERS_DIR, exist_ok=True)
json.dump(R, open(_out, "w"), indent=2)
_sys.path.insert(0, paths.ANALYSIS_DIR)
import cohort_spec as _cs
print("wrote", _out, "and", _cs.sidecar(_out, __file__, note="the additional analyses of the second read"))
