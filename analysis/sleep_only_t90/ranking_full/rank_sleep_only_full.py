#!/usr/bin/env python3
"""The 197-measure ranking (build_ranking_v3 procedure, unchanged arithmetic) rerun inside the
FULL 19,173 cohort, with SLEEP-ONLY oximetry measures entered as extra
candidates: sleep_only_T90, sleep_only_spo2_mean, sleep_only_spo2_nadir (and wake_only_T90 for
contrast). Answers: does T90 measured during sleep still rank first among everything?

Outputs (this folder): ranking_sleeponly_full.csv, ranking_sleeponly_full_percondition.csv,
PROVENANCE.json. Work/checkpoints in .work/. Nothing under numbers/ is written."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import warnings; warnings.filterwarnings("ignore")
import os, sys, json, time, glob, shutil
import numpy as np, pandas as pd
from scipy import stats
from patsy import dmatrix
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from joblib import Parallel, delayed

ROOT = paths.T90_ROOT
sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
from disease_definitions import DISEASES, RANKING_EXCLUDE
from cohort_spec import apply_cohort, COHORT_N
CIRCULAR = {"insomnia", "rls_plmd", "nocturia", "epilepsy"}
MASTER = f"{paths.FEATURE_OUTPUTS_DIR}/master_cohort_v2.csv"
OXY = f"{paths.SV_ROOT}/oxygen_profile/oxyprofile_per_patient.csv"
ADJ = [f"age_s{i}" for i in range(4)] + ["male"]
HERE = os.path.dirname(os.path.abspath(__file__)); WORK = f"{HERE}/.work"; os.makedirs(WORK, exist_ok=True)
PKL = f"{WORK}/frame.pkl"; TAG = "sleeponly_full"
def say(*a): print(f"[{TAG}]", *a, flush=True)

# ---------------------------------------------------------------- cohort = ranking cohort with sleep-only oximetry from the 08-21 full-cohort extraction
base = apply_cohort(pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet"))
assert len(base) == COHORT_N == 19173, len(base)
ox = pd.read_csv(OXY, low_memory=False)
SLEEPCOLS = {"sleep_t90": "sleep_only_T90", "sleep_mean": "sleep_only_spo2_mean", "sleep_nadir": "sleep_only_spo2_nadir", "wake_t90": "wake_only_T90"}
ox = ox[["BDSPPatientID", "sleep_valid_min", "wake_valid_min"] + list(SLEEPCOLS)].rename(columns=SLEEPCOLS)
sub = base.merge(ox, on="BDSPPatientID", how="inner")
sub = sub[(sub.sleep_valid_min >= 10.0) & (sub.wake_valid_min >= 10.0)].dropna(subset=["sleep_only_T90", "wake_only_T90"]).reset_index(drop=True)
say(f"cohort with sleep-only oximetry (>= 10 valid min asleep and awake): {len(sub):,} of 19,173; "
    f"sleep-only vs recording T90 Spearman {stats.spearmanr(sub.sleep_only_T90, sub.spo2_pct_below_90)[0]:.3f}")
assert len(sub) >= 18000, len(sub)

head = pd.read_csv(MASTER, nrows=0).columns.tolist()
FEAT = [c for c in head if c not in ("BDSPPatientID",) and not c.endswith(("_first_date", "_date"))]
raw = pd.read_csv(MASTER, usecols=["BDSPPatientID"] + FEAT, low_memory=False)
left = sub[["BDSPPatientID", "site_id", "AgeAtVisit", "sex"] + list(SLEEPCOLS.values()) +
           [c for c in sub.columns if c.endswith(("_incident", "_years", "_prevalent"))]]
raw = raw.drop(columns=[c for c in raw.columns if c in set(left.columns) - {"BDSPPatientID"}])
FEAT = [c for c in FEAT if c in raw.columns]
d = left.merge(raw, on="BDSPPatientID", how="left")
d["male"] = (d.sex.astype(str).str.upper().str[0] == "M").astype(int)

# the published 197 = the measures ranking_v3 kept; keep exactly those (so ranks are comparable),
# plus the four stage-resolved oximetry candidates
RK = pd.read_csv(f"{paths.NUMBERS_DIR}/ranking_v3.csv", comment="#")
pub = [f for f in RK.feature if f in d.columns]
missing = [f for f in RK.feature if f not in d.columns]
assert not missing, missing[:5]
num = pub + list(SLEEPCOLS.values())
say(f"candidate measures: {len(pub)} published + {len(SLEEPCOLS)} stage-resolved = {len(num)}")

sp = dmatrix("cr(a, df=4) - 1", {"a": d.AgeAtVisit.values}, return_type="dataframe")
sp.columns = [f"age_s{i}" for i in range(4)]; sp.index = d.index
for c in sp.columns: d[c] = sp[c]
for f in num:
    z = pd.Series(np.nan, index=d.index)
    for s, idx in d.groupby("site_id").groups.items():
        v = d.loc[idx, f]; ok = v.notna()
        if ok.sum() < 20: continue
        r = stats.rankdata(v[ok], method="average")
        z.loc[v[ok].index] = stats.norm.ppf((r - 0.375) / (ok.sum() + 0.25))
    d[f + "__z"] = z

# the published 48 outcomes (same rule as ranking_v3, evaluated on the FULL cohort so the set
# is identical), each must still be estimable both ways inside the subset (>= 30 events per side)
OUT_FULL = [k for k in DISEASES if k not in CIRCULAR and k not in RANKING_EXCLUDE and f"{k}_incident" in base.columns]
OUT_FULL = [k for k in OUT_FULL if int(base[f"{k}_incident"].sum()) >= 150] + ["death"]
assert len(OUT_FULL) == 48, len(OUT_FULL)
def estimable(k):
    yc, ec = f"{k}_years", f"{k}_incident"
    g = d[(d[f"{k}_prevalent"] == 0) & (d[yc] > 0) & d[yc].notna()]
    return all(int(g[g.site_id == s][ec].sum()) >= 30 for s in ("I0002", "I0006"))
OUT = [k for k in OUT_FULL if estimable(k)]
say(f"outcomes: {len(OUT)} of the published 48 estimable; dropped: {[k for k in OUT_FULL if k not in OUT]}")
keep = ADJ + [f + "__z" for f in num] + ["site_id"] + [f"{o}_{s}" for o in OUT for s in ("years", "incident", "prevalent")]
d[[c for c in dict.fromkeys(keep) if c in d.columns]].to_pickle(PKL)

_CACHE = {}
def run_outcome(out, feats):
    ck = f"{WORK}/{out}.json"
    if os.path.exists(ck):
        try: return json.load(open(ck))
        except Exception: pass
    if "d" not in _CACHE: _CACHE["d"] = pd.read_pickle(PKL)
    dd = _CACHE["d"]; yc, ec = f"{out}_years", f"{out}_incident"
    g = dd[(dd[f"{out}_prevalent"] == 0) & (dd[yc] > 0) & dd[yc].notna()]
    side = {s: g[g.site_id == s] for s in ("I0002", "I0006")}
    rows = []
    for feat in [None] + list(feats):
        cols = ADJ + ([feat + "__z"] if feat else []); res = []
        for tr, te in (("I0002", "I0006"), ("I0006", "I0002")):
            A = side[tr][cols + [yc, ec]].dropna(); B = side[te][cols + [yc, ec]].dropna()
            if A[ec].sum() < 30 or B[ec].sum() < 30: res.append(np.nan); continue
            try:
                c = CoxPHFitter(penalizer=0.01).fit(A, duration_col=yc, event_col=ec)
                res.append(concordance_index(B[yc], -c.predict_partial_hazard(B), B[ec]))
            except Exception: res.append(np.nan)
        rows.append({"feature": feat or "__BASE__", "outcome": out, "c": float(np.nanmean(res))})
    tmp = ck + ".part"; json.dump(rows, open(tmp, "w")); os.replace(tmp, ck); return rows

done = [o for o in OUT if os.path.exists(f"{WORK}/{o}.json")]
if done: say(f"resuming, {len(done)} of {len(OUT)} outcomes checkpointed")
say(f"fits: {(len(num) + 1) * len(OUT) * 2}"); t0 = time.time()
rr = Parallel(n_jobs=7, verbose=5, backend="loky")(delayed(run_outcome)(o, num) for o in OUT)
say(f"fitting took {(time.time() - t0) / 60:.1f} min")
s = pd.DataFrame([r for chunk in rr for r in chunk])
bas = s[s.feature == "__BASE__"].set_index("outcome")["c"]
f2 = s[s.feature != "__BASE__"].copy(); f2["gain"] = f2.apply(lambda r: r.c - bas.get(r.outcome, np.nan), axis=1)
rank = (f2.groupby("feature").agg(mC=("c", "mean"), dC=("gain", "mean"), worst=("gain", "min"),
        npos=("gain", lambda x: int((x > 0).sum())), nout=("gain", "size")).reset_index().sort_values("dC", ascending=False))
rank.insert(0, "rank", range(1, len(rank) + 1))
rank["rank_published_v3"] = rank.feature.map(RK.set_index("feature")["rank"])
rank["dC_published_v3"] = rank.feature.map(RK.set_index("feature")["dC"])
HDR = ["# ranking in the full cohort with sleep-only oximetry (08-21 extraction), built " + time.strftime("%Y-%m-%d"),
       f"# procedure identical to numbers/build_ranking_v3.py; measures: the published {len(pub)} + sleep_only_T90,",
       "#   sleep_only_spo2_mean, sleep_only_spo2_nadir, wake_only_T90 (from the per-stage oximetry extraction)",
       f"# outcomes: {len(OUT)} of the published 48; baseline held-out C {bas.mean():.6f}",
       "# readers must pass comment='#' to pandas.read_csv."]
for path, frame in ((f"{HERE}/ranking_sleeponly_full.csv", rank), (f"{HERE}/ranking_sleeponly_full_percondition.csv", f2)):
    with open(path, "w") as fh: fh.write("\n".join(HDR) + "\n"); frame.to_csv(fh, index=False)
json.dump({"built": time.strftime("%Y-%m-%d %H:%M"), "n": int(len(sub)), "measures": len(num), "outcomes": OUT,
           "baseline_heldout_C": float(bas.mean()), "sites": sub.site_id.value_counts().to_dict()},
          open(f"{HERE}/PROVENANCE.json", "w"), indent=1)
say(f"baseline held-out C {bas.mean():.6f}")
print(rank.head(12).to_string(index=False), flush=True)
print(rank[rank.feature.isin(["spo2_pct_below_90", "sleep_only_T90", "wake_only_T90", "sleep_only_spo2_mean", "sleep_only_spo2_nadir",
                              "spo2_mean", "spo2_nadir_corrected", "TST_min", "ahi_total", "AHI", "arousal_index", "sleep_efficiency_pct"])].to_string(index=False), flush=True)
say("DONE")
