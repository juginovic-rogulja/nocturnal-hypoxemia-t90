"""
The treatment arm, rebuilt across all 52 conditions.

Two questions, in order:
  1  Among patients who all began with substantial oxygen debt, does it matter whether
     treatment corrected it? Everyone here started above 10%, so severity is not the
     difference between the groups. What differs is the outcome of treatment.
  2  For which conditions does correcting it change risk, and by how much? Reported as the
     percentage reduction in hazard, not only as whether it clears significance, because a
     30% reduction that misses significance in 200 events is still worth reporting as such.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../..")))  # repository root, where paths.py lives
import paths
import warnings; warnings.filterwarnings("ignore")
import json, sys
import numpy as np, pandas as pd
from lifelines import CoxPHFitter
sys.path.insert(0, paths.ANALYSIS_DIR)   # the sibling modules of this folder, whatever the working directory is
from disease_definitions import DISEASES, NEGATIVE_CONTROLS

# lifelines scales its penalty by the sample size, so penalizer=0.01 on a cohort this
# size is an effective ridge in the hundreds. Every fit here is by maximum partial
# likelihood, matching the primary analysis and standard survival software.

S=["wake","sleep","nrem","n1","n2","n3","rem"]
d=pd.read_parquet(f"{paths.TABLES_DIR}/cpap_t90_by_stage.parquet")
A=d[d.design=="A_split"].copy(); A["src"]="A"
B=d[d.design=="B_pairs"].copy(); B["src"]="B"
P=d[d.design=="P_placebo"]
for s in S:
    for a,b in [("pre","dx"),("post","tx")]:
        for suf in ["t90","min"]:
            c=f"{b}_{s}_{suf}"; B[f"{a}_{s}_{suf}"]=B[c] if c in B.columns else np.nan
keep=["BDSPPatientID","SexDSC","src"]+[f"{p}_{s}_{x}" for p in ("pre","post") for s in S for x in ("t90","min")]
cc=pd.concat([A[[c for c in keep if c in A.columns]],B[[c for c in keep if c in B.columns]]],ignore_index=True)
cc=cc.sort_values("src",ascending=False).drop_duplicates("BDSPPatientID",keep="first")

base=pd.read_parquet(f'{paths.TABLES_DIR}/t90_final.parquet')
m=cc.merge(base,on="BDSPPatientID",how="inner")
m=m[(m.fu_valid==1)&(m.post_sleep_min>=30)&(m.pre_sleep_min>=30)].copy()
# same oximetry quality rule as the cohort file: a recording cannot spend most of the
# night below 90% saturation while scoring no desaturation events
for _side in ("pre", "post"):
    _n, _mn = f"{_side}_sleep_nadir", f"{_side}_sleep_mean"
    if _n in m.columns and _mn in m.columns:
        m = m[~(m[_n] > m[_mn])]
    _t = f"{_side}_sleep_t90"
    if _t in m.columns and _mn in m.columns:
        m = m[~(m[_mn].between(50, 80) & (m[_t] > 70))]

m["male"]=(m.SexDSC.astype(str).str.upper().str[0]=="M").astype(int)
m["lbase"]=np.log1p(m.pre_sleep_t90)

R={"designs":{"split_night":int(len(A)),"pairs":int(len(B)),"placebo":int(len(P))}}

# ---- how much does treatment move oxygen, by stage
eff={}
for s in S:
    g=m[(m[f"pre_{s}_min"]>=30)&(m[f"post_{s}_min"]>=30)]
    g=g[g[f"pre_{s}_t90"].notna()&g[f"post_{s}_t90"].notna()]
    if len(g)<30: continue
    pre,post=g[f"pre_{s}_t90"],g[f"post_{s}_t90"]
    from scipy import stats as st
    eff[s]={"n":int(len(g)),"pre_median":round(float(pre.median()),2),
            "post_median":round(float(post.median()),2),
            "pct_reduction":round(100*(1-post.median()/pre.median()),1) if pre.median()>0 else None,
            "fell_pct":round(100*float((post<pre).mean()),1),
            "p":float(st.wilcoxon(pre,post).pvalue)}
R["effect_by_stage"]=eff

# ---- band migration
BANDS=[(-1,1,"0-1%"),(1,5,"1-5%"),(5,10,"5-10%"),(10,1e9,">10%")]
bd=lambda v: next(i for i,(lo,hi,_) in enumerate(BANDS) if lo<v<=hi)
g=m[m.pre_sleep_t90.notna()&m.post_sleep_t90.notna()].copy()
g["pb"]=g.pre_sleep_t90.apply(bd); g["qb"]=g.post_sleep_t90.apply(bd)
R["migration"]={"n":int(len(g)),
  "matrix":[[int(((g.pb==i)&(g.qb==j)).sum()) for j in range(4)] for i in range(4)],
  "pre_totals":[int((g.pb==i).sum()) for i in range(4)],
  "post_totals":[int((g.qb==j).sum()) for j in range(4)],
  "pct_moved_down":round(100*float((g.qb<g.pb).mean()),1),
  "labels":[x[2] for x in BANDS]}
worst=g[g.pb==3]
R["migration"]["from_worst"]={"n":int(len(worst)),
  "to":{BANDS[j][2]:int((worst.qb==j).sum()) for j in range(4)}}

# ---- THE KEY TEST: all started >10%, corrected vs not
sev=m[m.pre_sleep_t90>10].copy()
sev["notcorr"]=(sev.post_sleep_t90>10).astype(int)
R["corrected_vs_not"]={"n":int(len(sev)),
  "n_corrected":int((1-sev.notcorr).sum()),"n_not":int(sev.notcorr.sum()),
  "baseline_corrected":round(float(sev.loc[sev.notcorr==0,"pre_sleep_t90"].median()),1),
  "baseline_not":round(float(sev.loc[sev.notcorr==1,"pre_sleep_t90"].median()),1),
  "outcomes":{}}
KEYS=list(DISEASES.keys())+["death"]
for k in KEYS:
    if f"{k}_incident" not in sev.columns: continue
    lab=DISEASES[k][0] if k in DISEASES else "Death from any cause"
    neg=k in NEGATIVE_CONTROLS   # the raw flag still carries thyroid disease, dropped earlier
    f=sev[(sev[f"{k}_prevalent"]==0)&sev[f"{k}_years"].notna()&(sev[f"{k}_years"]>0)]
    ff=pd.DataFrame({"T":f[f"{k}_years"],"E":f[f"{k}_incident"].astype(int),"x":f.notcorr,
        "age":f.AgeAtVisit,"male":f.male,"lbase":f.lbase,"site":f.site_id}).dropna()
    if ff.E.sum()<40: continue
    try: c=CoxPHFitter(penalizer=0.0).fit(ff,"T","E",strata=["site"])
    except Exception: continue
    lo,hi=np.exp(c.confidence_intervals_.loc["x"])
    hr=float(np.exp(c.params_["x"]))
    R["corrected_vs_not"]["outcomes"][lab]={"events":int(ff.E.sum()),
      "hr":hr,"lo":float(lo),"hi":float(hi),  # v8.3 full precision (2026-09-16)
      "p":float(c.summary.loc["x","p"]),
      "pct_lower_if_corrected":100*(1-1/hr) if hr>0 else None,  # v8.3 full precision (2026-09-16)
      "negative_control":neg}

# the output goes to the declared output folder, never to the working directory, and carries a
# provenance sidecar (sha256 of itself, of this script and of the four tables) as every other
# output of this folder does.
_out = f"{paths.NUMBERS_DIR}/treatment_v2.json"
_os.makedirs(paths.NUMBERS_DIR, exist_ok=True)
json.dump(R, open(_out, "w"), indent=2)
sys.path.insert(0, paths.ANALYSIS_DIR)
import cohort_spec as _cs
print("wrote", _out, "and", _cs.sidecar(_out, __file__, note="the treatment arm across all 52 conditions"))
o=R["corrected_vs_not"]
print(f"all started >10%: n={o['n']:,}  corrected {o['n_corrected']:,}  not {o['n_not']:,}")
print(f"baseline T90 {o['baseline_corrected']}% vs {o['baseline_not']}%\n")
rows=sorted(o["outcomes"].items(),key=lambda kv:-kv[1]["hr"])
print(f"{'condition':<28}{'ev':>6}{'HR':>8}{'95% CI':>17}{'lower if corrected':>22}")
for lab,v in rows:
    if v["events"]<40: continue
    star="*" if v["p"]<.05 else " "
    ci=f"({v['lo']:.2f}-{v['hi']:.2f})"
    red=f"{v['pct_lower_if_corrected']}%" if v["hr"]>1 else "-"
    neg="  NEG" if v["negative_control"] else ""
    print(f"  {lab:<26}{v['events']:>6}{v['hr']:>8.2f}{star}{ci:>16}{red:>22}{neg}")
