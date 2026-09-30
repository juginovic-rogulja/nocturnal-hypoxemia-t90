"""
Every analysis the rewritten paper needs, on the rebuilt outcome set.

One file so the whole paper regenerates from one command. Order follows the argument:
what wins, what it is associated with, where it hides, does correcting it help, and does it
hold elsewhere.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../..")))  # repository root, where paths.py lives
import paths
import warnings; warnings.filterwarnings("ignore")
import json, sys
import numpy as np, pandas as pd
from scipy import stats
from patsy import dmatrix
from lifelines import CoxPHFitter, KaplanMeierFitter
sys.path.insert(0, paths.ANALYSIS_DIR)   # the sibling modules of this folder, whatever the working directory is
from disease_definitions import DISEASES, NEGATIVE_CONTROLS

# lifelines scales its penalty by the sample size, so penalizer=0.01 on a cohort this
# size is an effective ridge in the hundreds. Every fit here is by maximum partial
# likelihood, matching the primary analysis and standard survival software.

R={}
b=pd.read_parquet(f'{paths.TABLES_DIR}/t90_final.parquet')
h=pd.read_parquet(f'{paths.TABLES_DIR}/t90_base4.parquet')[['BDSPPatientID','osa_prevalent']]
b=b.merge(h,on='BDSPPatientID',how='left')
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
BANDS=[(-1,1,"0-1%"),(1,5,"1-5%"),(5,10,"5-10%"),(10,1e9,">10%")]
b["band"]=b.t90.apply(lambda v: next(i for i,(lo,hi,_) in enumerate(BANDS) if lo<v<=hi))

def cox(g,xcols,k,strat=True):
    f=g[(g[f"{k}_prevalent"]==0)&g[f"{k}_years"].notna()&(g[f"{k}_years"]>0)]
    cols={"T":f[f"{k}_years"],"E":f[f"{k}_incident"].astype(int)}
    for c in xcols+ADJ: cols[c]=f[c]
    if strat: cols["site"]=f.site_id
    ff=pd.DataFrame(cols).dropna()
    if ff.E.sum()<60: return None
    try: c=CoxPHFitter(penalizer=0.0).fit(ff,"T","E",strata=["site"] if strat else None)
    except Exception: return None
    out={"n":int(len(ff)),"events":int(ff.E.sum())}
    for x in xcols:
        lo,hi=np.exp(c.confidence_intervals_.loc[x])
        out[x]={"hr":float(np.exp(c.params_[x])),"lo":float(lo),  # v8.3 full precision (2026-09-16)
                "hi":float(hi),"p":float(c.summary.loc[x,"p"])}
    return out

# ---------------------------------------------------------------- 1 cohort description
R["cohort"]={"n":int(len(b)),"age_median":float(b.AgeAtVisit.median()),
  "age_q1":float(b.AgeAtVisit.quantile(.25)),"age_q3":float(b.AgeAtVisit.quantile(.75)),
  "male_pct":round(100*float(b.male.mean()),1),
  "t90_median":round(float(b.t90.median()),2),
  "t90_bands_pct":{lab:round(100*float(((b.t90>lo)&(b.t90<=hi)).mean()),1) for lo,hi,lab in BANDS},
  "t90_bands_n":{lab:int(((b.t90>lo)&(b.t90<=hi)).sum()) for lo,hi,lab in BANDS},
  "sites":{str(k):int(v) for k,v in b.site_id.value_counts().items()},
  "fu_median":round(float(b.death_years.median()),2)}

# ---------------------------------------------------------------- 2 AHI vs T90
av=b[b.AHI.notna()]
AB=[(-1,5,"None, AHI <5"),(5,15,"Mild, AHI 5-14.9"),(15,30,"Moderate, AHI 15-29.9"),(30,1e9,"Severe, AHI >=30")]
R["ahi_vs_t90"]={"n":int(len(av)),
  "spearman":float(stats.spearmanr(av.AHI,av.t90).statistic),"rows":[]}  # v8.3 full precision (2026-09-16)
for lo,hi,lab in AB:
    g=av[(av.AHI>lo)&(av.AHI<=hi)] if lo>=0 else av[av.AHI<hi]
    R["ahi_vs_t90"]["rows"].append({"category":lab,"n":int(len(g)),
      **{f"t90_{l}":round(100*float(((g.t90>a)&(g.t90<=c)).mean()),1) for a,c,l in BANDS},
      "t90_median":round(float(g.t90.median()),2)})

# ---------------------------------------------------------------- 3 hidden hypoxemia
na=av[av.AHI<5].copy(); na["hi"]=(na.t90>10).astype(int)
R["hidden"]={"n":int(len(na)),"n_hi":int(na.hi.sum()),
  "pct_hi":round(100*float(na.hi.mean()),1),
  "age_hi":float(na.loc[na.hi==1,"AgeAtVisit"].median()),
  "age_lo":float(na.loc[na.hi==0,"AgeAtVisit"].median()),"outcomes":{}}
for k,(lab,neg,_,_) in DISEASES.items():
    if f"{k}_incident" not in na.columns: continue
    r=cox(na,["hi"],k)
    if r: R["hidden"]["outcomes"][lab]={"events":r["events"],**r["hi"],"negative_control": k in NEGATIVE_CONTROLS}
r=cox(na,["hi"],"death"); R["hidden"]["outcomes"]["Death from any cause"]={"events":r["events"],**r["hi"],"negative_control":False} if r else None

# ---------------------------------------------------------------- 4 T90 vs sleep duration
b["t90_hi"]=(b.t90>10).astype(int); b["slp_lt5"]=(b.hrs<5).astype(float).where(b.hrs.notna())   # v8.1: missing where TST is masked (split nights), so those rows drop out of every sleep-duration model
R["vs_sleep_duration"]={"n":int(len(b)),
  "pct_t90_hi":round(100*float(b.t90_hi.mean()),1),
  "pct_short":round(100*float(b.slp_lt5.mean()),1),"outcomes":{}}
for k,(lab,neg,_,_) in list(DISEASES.items())+[("death",("Death from any cause",False,[],[]))]:
    if f"{k}_incident" not in b.columns: continue
    r=cox(b,["t90_hi","slp_lt5"],k)
    if r: R["vs_sleep_duration"]["outcomes"][lab]={"events":r["events"],
        "t90":r["t90_hi"],"short_sleep":r["slp_lt5"],"negative_control": k in NEGATIVE_CONTROLS}

# ---------------------------------------------------------------- 5 graded response
R["graded"]={"band_n":R["cohort"]["t90_bands_n"],"outcomes":{}}
for k,(lab,neg,_,_) in list(DISEASES.items())+[("death",("Death from any cause",False,[],[]))]:
    if f"{k}_incident" not in b.columns: continue
    g=b[(b[f"{k}_prevalent"]==0)&b[f"{k}_years"].notna()&(b[f"{k}_years"]>0)]
    X=pd.get_dummies(g.band,prefix="bd").astype(int)
    for i in range(4):
        if f"bd_{i}" not in X: X[f"bd_{i}"]=0
    X=X[[f"bd_{i}" for i in range(4)]].drop(columns=["bd_0"])
    f=pd.concat([pd.DataFrame({"T":g[f"{k}_years"].values,"E":g[f"{k}_incident"].astype(int).values,
        "site":g.site_id.values,**{c:g[c].values for c in ADJ}}),X.reset_index(drop=True)],axis=1).dropna()
    if f.E.sum()<100: continue
    try: c=CoxPHFitter(penalizer=0.0).fit(f,"T","E",strata=["site"])
    except Exception: continue
    row={"events":int(f.E.sum()),"negative_control": k in NEGATIVE_CONTROLS,"0-1%":{"hr":1.0}}
    for i in range(1,4):
        nm=f"bd_{i}"
        if nm in c.params_.index:
            lo,hi=np.exp(c.confidence_intervals_.loc[nm])
            row[BANDS[i][2]]={"hr":float(np.exp(c.params_[nm])),"lo":float(lo),  # v8.3 full precision (2026-09-16)
                              "hi":float(hi),"p":float(c.summary.loc[nm,"p"])}
    ft=pd.DataFrame({"T":g[f"{k}_years"].values,"E":g[f"{k}_incident"].astype(int).values,
        "band":g.band.values,"site":g.site_id.values,**{c_:g[c_].values for c_ in ADJ}}).dropna()
    try:
        ct=CoxPHFitter(penalizer=0.0).fit(ft,"T","E",strata=["site"])
        row["trend_p"]=float(ct.summary.loc["band","p"])
        row["trend_hr"]=float(np.exp(ct.params_["band"]))  # v8.3 full precision (2026-09-16)
    except Exception: pass
    R["graded"]["outcomes"][lab]=row

# the output goes to the declared output folder, never to the working directory, and carries a
# provenance sidecar (sha256 of itself, of this script and of the four tables) as every other
# output of this folder does.
_out = f"{paths.NUMBERS_DIR}/results_v2.json"
_os.makedirs(paths.NUMBERS_DIR, exist_ok=True)
json.dump(R, open(_out, "w"), indent=2)
sys.path.insert(0, paths.ANALYSIS_DIR)
import cohort_spec as _cs
print("wrote", _out, "and", _cs.sidecar(_out, __file__, note="every analysis of the rewritten paper on the rebuilt outcome set"))
print("cohort", R["cohort"]["n"])
print("hidden hypoxemia:", R["hidden"]["n_hi"], "of", R["hidden"]["n"], f"({R['hidden']['pct_hi']}%)")
print("graded outcomes:", len(R["graded"]["outcomes"]))
print("vs sleep duration:", len(R["vs_sleep_duration"]["outcomes"]))
