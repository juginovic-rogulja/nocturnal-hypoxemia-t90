"""
What does a measure that means nothing earn, just for being added?

Every measure in the ranking is scored by how much it raises held-out concordance over age and
sex. A variable carrying no information does not score exactly zero. It scores whatever the
procedure hands out for free: the penalized fit shrinks a useless coefficient towards zero but not
to zero, and the held-out concordance moves by a small amount in whichever direction the noise
happened to fall. That amount is the floor the real measures have to clear, and without it a gain
of +0.0002 cannot be told apart from nothing at all.

This runs the identical procedure on 20 variables drawn from a random number generator, on the
same 19,173 patients and the same 48 outcomes, fitted at one hospital and scored at the other.
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
from cohort_spec import apply_cohort, N_RANKED_OUTCOMES, ranked_outcome_keys  # v8 sweep 2026-09-12; v8.1 ranked set

CIRCULAR = {"insomnia", "rls_plmd", "nocturia", "epilepsy"}
ADJ = [f"age_s{i}" for i in range(4)] + ["male"]
PKL = f"{WORK}/noise_data.pkl"
NNOISE = 20
SEED = 20260807

_tbl = pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet")
base = apply_cohort(_tbl)
d = base[["BDSPPatientID", "site_id", "AgeAtVisit", "sex"] +
         [c for c in base.columns if c.endswith(("_incident", "_years", "_prevalent"))]].copy()
d = d.reset_index(drop=True)
d["male"] = (d.sex.astype(str).str.upper().str[0] == "M").astype(int)
sp = dmatrix("cr(a, df=4) - 1", {"a": d.AgeAtVisit.values}, return_type="dataframe")
sp.columns = ADJ[:4]; sp.index = d.index
for c in sp.columns:
    d[c] = sp[c]

# Noise, then the same site-wise rank-normal transform every real measure goes through, so the
# comparison is to the procedure and not to a different scaling.
rng = np.random.default_rng(SEED)
noise = [f"noise_{i:02d}" for i in range(NNOISE)]
for f in noise:
    d[f] = rng.standard_normal(len(d))
    z = pd.Series(np.nan, index=d.index)
    for s, idx in d.groupby("site_id").groups.items():
        v = d.loc[idx, f]; ok = v.notna()
        if ok.sum() < 20: continue
        r = stats.rankdata(v[ok], method="average")
        z.loc[v[ok].index] = stats.norm.ppf((r - 0.375) / (ok.sum() + 0.25))
    d[f + "__z"] = z

OUT = [k for k in DISEASES if k not in CIRCULAR and k not in RANKING_EXCLUDE
       and f"{k}_incident" in d.columns]
OUT = ranked_outcome_keys(_tbl, OUT) + ["death"]      # v8.1: the ranked set is fixed on the all-nights cohort
assert len(OUT) == N_RANKED_OUTCOMES, len(OUT)
keep = ADJ + [f + "__z" for f in noise] + ["site_id"] + \
       [f"{o}_{s}" for o in OUT for s in ("years", "incident", "prevalent")]
d[[c for c in dict.fromkeys(keep) if c in d.columns]].to_pickle(PKL)
print(f"cohort {len(d):,}   outcomes {len(OUT)}   noise variables {NNOISE}")

_CACHE = {}


def one(feat, out):
    if "d" not in _CACHE:
        _CACHE["d"] = pd.read_pickle(PKL)
    dd = _CACHE["d"]
    g = dd[(dd[f"{out}_prevalent"] == 0) & (dd[f"{out}_years"] > 0) & dd[f"{out}_years"].notna()]
    cols = ADJ + ([feat + "__z"] if feat else [])
    res = []
    for tr, te in (("I0002", "I0006"), ("I0006", "I0002")):
        A = g[g.site_id == tr][cols + [f"{out}_years", f"{out}_incident"]].dropna()
        B = g[g.site_id == te][cols + [f"{out}_years", f"{out}_incident"]].dropna()
        if A[f"{out}_incident"].sum() < 30 or B[f"{out}_incident"].sum() < 30:
            res.append(np.nan); continue
        try:
            c = CoxPHFitter(penalizer=0.01).fit(A, duration_col=f"{out}_years",
                                                event_col=f"{out}_incident")
            r = c.predict_partial_hazard(B)
            res.append(concordance_index(B[f"{out}_years"], -r, B[f"{out}_incident"]))
        except Exception:
            res.append(np.nan)
    return {"feature": feat or "__BASE__", "outcome": out, "c": float(np.nanmean(res))}


jobs = [(None, o) for o in OUT] + [(f, o) for f in noise for o in OUT]
print(f"fits: {len(jobs) * 2}")
rr = Parallel(n_jobs=7, verbose=5, backend="loky")(delayed(one)(f, o) for f, o in jobs)
s = pd.DataFrame(rr)
bas = s[s.feature == "__BASE__"].set_index("outcome")["c"]
f2 = s[s.feature != "__BASE__"].copy()
f2["gain"] = f2.apply(lambda r: r.c - bas.get(r.outcome, np.nan), axis=1)
per = (f2.groupby("feature").agg(dC=("gain", "mean"), worst=("gain", "min"),
       npos=("gain", lambda x: int((x > 0).sum()))).reset_index().sort_values("dC", ascending=False))
f2.to_csv(f"{paths.NUMBERS_DIR}/ranking_v3_noisefloor_percondition.csv", index=False)
per.to_csv(f"{paths.NUMBERS_DIR}/ranking_v3_noisefloor.csv", index=False)

g = f2["gain"]
summ = {
    "built": "2026-08-07", "cohort_n": int(len(d)), "outcomes": len(OUT),
    "n_noise_variables": NNOISE, "seed": SEED,
    "baseline_heldout_C": float(bas.mean()),
    "dC_per_noise_variable_mean": float(per.dC.mean()),
    "dC_per_noise_variable_sd": float(per.dC.std(ddof=1)),
    "dC_per_noise_variable_min": float(per.dC.min()),
    "dC_per_noise_variable_max": float(per.dC.max()),
    "dC_95pct_range_of_a_useless_variable": [float(per.dC.quantile(0.025)),
                                             float(per.dC.quantile(0.975))],
    "single_outcome_gain_sd": float(g.std(ddof=1)),
    "npos_mean_of_48": float(per.npos.mean()),
}
json.dump(summ, open(f"{paths.NUMBERS_DIR}/ranking_v3_noisefloor.json", "w"), indent=2)
print(json.dumps(summ, indent=2))
print(per.to_string(index=False))
