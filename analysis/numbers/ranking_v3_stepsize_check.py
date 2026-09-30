"""
Did the fits that decide first place actually converge?

lifelines takes a Newton step whose default size lands on a worse optimum on some fits, and a Cox
model that has quietly converged to the wrong place still returns a concordance, so nothing looks
broken. The whole first-versus-second question turns on a difference of a few thousandths, which is
well inside the range a bad optimum can produce.

So the two contenders are refitted at six step sizes on all 48 outcomes and both folds, the
solution with the highest partial log-likelihood is kept, and the resulting gain is compared with
what ranking_v3 recorded. If they agree, the ranking's fits were already at the right optimum.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../..")))  # repository root, where paths.py lives
import paths
import warnings; warnings.filterwarnings("ignore")
import sys, json
import numpy as np, pandas as pd
from scipy import stats
from patsy import dmatrix
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from joblib import Parallel, delayed

ROOT = paths.T90_ROOT
# Work files stay inside the project. A shared temporary directory was cleared by another
# job mid-run once already, which killed 4,036 completed fits.
WORK = f"{paths.NUMBERS_DIR}/.rank_v3_work"
import os as _os; _os.makedirs(WORK, exist_ok=True)
sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
from disease_definitions import DISEASES, RANKING_EXCLUDE
from cohort_spec import apply_cohort

CIRCULAR = {"insomnia", "rls_plmd", "nocturia", "epilepsy"}
MASTER = f"{paths.FEATURE_OUTPUTS_DIR}/master_cohort_v2.csv"
ADJ = [f"age_s{i}" for i in range(4)] + ["male"]
PKL = f"{WORK}/step_data.pkl"
FEATURES = ["spo2_pct_below_90", "spo2_nadir_corrected", "odi3_total", "spo2_pct_below_88"]
STEPS = [None, 0.05, 0.1, 0.25, 0.5, 0.95]

base = apply_cohort(pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet"))
raw = pd.read_csv(MASTER, usecols=["BDSPPatientID"] + FEATURES, low_memory=False)
left = base[["BDSPPatientID", "site_id", "AgeAtVisit", "sex"] +
            [c for c in base.columns if c.endswith(("_incident", "_years", "_prevalent"))]]
d = left.merge(raw, on="BDSPPatientID", how="left")
d["male"] = (d.sex.astype(str).str.upper().str[0] == "M").astype(int)
sp = dmatrix("cr(a, df=4) - 1", {"a": d.AgeAtVisit.values}, return_type="dataframe")
sp.columns = ADJ[:4]; sp.index = d.index
for c in sp.columns:
    d[c] = sp[c]
for f in FEATURES:
    z = pd.Series(np.nan, index=d.index)
    for s, idx in d.groupby("site_id").groups.items():
        v = d.loc[idx, f]; ok = v.notna()
        if ok.sum() < 20: continue
        r = stats.rankdata(v[ok], method="average")
        z.loc[v[ok].index] = stats.norm.ppf((r - 0.375) / (ok.sum() + 0.25))
    d[f + "__z"] = z

OUT = [k for k in DISEASES if k not in CIRCULAR and k not in RANKING_EXCLUDE
       and f"{k}_incident" in d.columns]
OUT = [k for k in OUT if int(d[f"{k}_incident"].sum()) >= 150] + ["death"]
keep = ADJ + [f + "__z" for f in FEATURES] + ["site_id"] + \
       [f"{o}_{s}" for o in OUT for s in ("years", "incident", "prevalent")]
d[[c for c in dict.fromkeys(keep) if c in d.columns]].to_pickle(PKL)
print(f"cohort {len(d):,}  outcomes {len(OUT)}  measures {len(FEATURES)}  steps {STEPS}")

_CACHE = {}


def one(feat, out):
    if "d" not in _CACHE:
        _CACHE["d"] = pd.read_pickle(PKL)
    dd = _CACHE["d"]
    g = dd[(dd[f"{out}_prevalent"] == 0) & (dd[f"{out}_years"] > 0) & dd[f"{out}_years"].notna()]
    cols = ADJ + ([feat + "__z"] if feat else [])
    default, best = [], []
    changed = 0
    for tr, te in (("I0002", "I0006"), ("I0006", "I0002")):
        A = g[g.site_id == tr][cols + [f"{out}_years", f"{out}_incident"]].dropna()
        B = g[g.site_id == te][cols + [f"{out}_years", f"{out}_incident"]].dropna()
        if A[f"{out}_incident"].sum() < 30 or B[f"{out}_incident"].sum() < 30:
            default.append(np.nan); best.append(np.nan); continue
        cands = []
        for st in STEPS:
            try:
                kw = {} if st is None else {"step_size": st}
                c = CoxPHFitter(penalizer=0.01).fit(A, duration_col=f"{out}_years",
                                                    event_col=f"{out}_incident", **kw)
                ci = concordance_index(B[f"{out}_years"], -c.predict_partial_hazard(B),
                                       B[f"{out}_incident"])
                cands.append((float(c.log_likelihood_), ci, st))
            except Exception:
                continue
        if not cands:
            default.append(np.nan); best.append(np.nan); continue
        dflt = [x for x in cands if x[2] is None]
        default.append(dflt[0][1] if dflt else np.nan)
        top = max(cands, key=lambda x: x[0])
        best.append(top[1])
        if dflt and top[0] - dflt[0][0] > 1e-6:
            changed += 1
    return {"feature": feat or "__BASE__", "outcome": out,
            "c_default": float(np.nanmean(default)), "c_best": float(np.nanmean(best)),
            "folds_where_default_was_not_the_max_likelihood": changed}


jobs = [(None, o) for o in OUT] + [(f, o) for f in FEATURES for o in OUT]
print(f"fits: {len(jobs) * 2 * len(STEPS)}")
rr = Parallel(n_jobs=1, verbose=5)(delayed(one)(f, o) for f, o in jobs)
s = pd.DataFrame(rr)
s.to_csv(f"{paths.NUMBERS_DIR}/ranking_v3_stepsize_check.csv", index=False)

res = {"steps_scanned": [str(x) for x in STEPS], "outcomes": len(OUT),
       "folds_total": int(len(s) * 2),
       "folds_where_default_missed_the_max": int(s.folds_where_default_was_not_the_max_likelihood.sum())}
for key, col in (("default", "c_default"), ("best", "c_best")):
    bas = s[s.feature == "__BASE__"].set_index("outcome")[col]
    f2 = s[s.feature != "__BASE__"].copy()
    f2["gain"] = f2.apply(lambda r: r[col] - bas.get(r.outcome, np.nan), axis=1)
    res[key] = {f: float(f2[f2.feature == f]["gain"].mean()) for f in FEATURES}
res["max_abs_dC_shift"] = max(abs(res["best"][f] - res["default"][f]) for f in FEATURES)
json.dump(res, open(f"{paths.NUMBERS_DIR}/ranking_v3_stepsize_check.json", "w"), indent=2)
print(json.dumps(res, indent=2))
