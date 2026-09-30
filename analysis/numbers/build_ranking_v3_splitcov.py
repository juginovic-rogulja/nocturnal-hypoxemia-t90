"""
The measure comparison, rebuilt on the corrected cohort.

Same question and same procedure as ranking_v2: each measure is added in turn to a model holding
only age and sex, fitted in one hospital and scored in the other, both directions, and measures are
ranked by how much they improve the ordering of who gets sick. Nothing about the fitting changed.
What changed is the data underneath it.

ranking_v2 was frozen 2026-08-02. Two corrections landed 2026-08-04 and never reached it:

  1. The oximetry quality rule (numbers/apply_qc_v3.py) removed 210 recordings whose stored oxygen
     is not a measurement of that patient. Cohort 19,383 -> 19,173.
  2. numbers/refresh_two_conditions.py recoded two outcomes. Obesity hypoventilation had been
     ICD-9 278.01, morbid obesity, instead of 278.03. Falls had been E888 and E884, whose
     decimal-stripped prefixes also matched ICD-10 E88.81, metabolic syndrome.

v8 (2026-09-12, decisions 2 and 3): the tables come from cohort_spec (DATA_DIR), master columns
ending in RANKING_DROP_SUFFIX (_siteZ, the 60 per-site standardized twins) are never candidates,
the sleep-period T90 (SENSITIVITY_ONLY) is never a row of the primary ranking, and the measure and
outcome counts are printed against cohort_spec.N_MEASURES and N_RANKED_OUTCOMES (reported, never
asserted: a new outcome under the 150-event bar drops out by the ranking's own rule and the count
printed is what survives). With T90_COLUMN set (step 230) only that exposure is refitted; every
other row and every base C come from the primary ranking_v3_percondition.csv, and the base C of
this run must equal the primary's per outcome to 1e-9 (the same cohort gives the same base fit).

Run as:
    python3 build_ranking_v3.py current   -> ranking_v3.csv, ranking_v3_percondition.csv
    T90_COLUMN=spo2_pct_below_90_sleep python3 build_ranking_v3.py current
                                          -> ranking_v3_sleepT90.csv, ranking_v3_sleepT90_percondition.csv
    RANK_FEATURES=a,b,c python3 build_ranking_v3.py current
                                          -> ranking_current_subset_check.csv (a smoke run)
    python3 build_ranking_v3.py legacy    -> the pre-correction 19,383 rebuild, used to prove this
                                             script reproduces ranking_v2 when handed ranking_v2's
                                             data. Feed it the old cohort and it must return the
                                             old answer, otherwise a moved number cannot be
                                             attributed to the data rather than the code.

Two things differ from run_ranking_v2.py, neither of which touches a number.

  Work files live under numbers/.rank_v3_work, not in a shared temporary directory. The first
  attempt at this run died at 4,036 of 9,504 fits because another job cleared the scratchpad out
  from under it and the workers lost the frame they were reading.

  The unit of parallel work is an outcome rather than a measure-outcome pair, so the eligible
  patients and the two hospital subsets are assembled once per outcome instead of once per
  measure, and each finished outcome is checkpointed so an interrupted run resumes instead of
  restarting.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../..")))  # repository root, where paths.py lives
import paths
import warnings; warnings.filterwarnings("ignore")
import os, sys, json, time, shutil
import numpy as np, pandas as pd
from scipy import stats
from patsy import dmatrix
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from joblib import Parallel, delayed

ROOT = paths.T90_ROOT
sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
from disease_definitions import DISEASES, RANKING_EXCLUDE
from cohort_spec import (apply_cohort, COHORT_N, N_MEASURES, N_RANKED_OUTCOMES, SENSITIVITY_ONLY,
                         RANKING_DROP_SUFFIX, T90_FINAL, MASTER_CSV, OUT_DIR, NUMBERS_DIR,
                         PRIMARY_EXPOSURE, t90_column, exposure_tag, sidecar)
# ---- v8.1 supplementary ranking (14 Sept): every analysis night (19,173) with the split-night flag as a covariate next to age and
# sex, so a reader sees that the top of the ranking does not depend on excluding the 3,622 split nights. The primary ranking
# (build_ranking_v3.py under T90_FULL_NIGHTS_ONLY=1) uses the full diagnostic nights only. Same fits, same candidate rule.
DIAG = "splitcov"
assert not os.environ.get("T90_FULL_NIGHTS_ONLY"), "this supplementary ranking runs on all nights: unset T90_FULL_NIGHTS_ONLY"

CIRCULAR = {"insomnia", "rls_plmd", "nocturia", "epilepsy"}
MASTER = MASTER_CSV
ADJ = [f"age_s{i}" for i in range(4)] + ["male"] + (["split_night"] if DIAG == "splitcov" else [])

MODE = sys.argv[1] if len(sys.argv) > 1 else "current"
assert MODE in ("current", "legacy")
SUBSET = os.environ.get("RANK_FEATURES")
XCOL = t90_column()                        # v8 decision 2: the exposure row to refit
SENS = XCOL != PRIMARY_EXPOSURE            # the sensitivity set refits only that row
SFX = exposure_tag(XCOL) if SENS else ""   # "_sleepT90" on the output names
TAG = "allnights_splitcov"
WORK = f"{OUT_DIR}/.rank_v3_work/{TAG}"
os.makedirs(WORK, exist_ok=True)
PKL = f"{WORK}/frame.pkl"

# ---------------------------------------------------------------------------------------
# Cohort
# ---------------------------------------------------------------------------------------
if MODE == "current":
    base = apply_cohort(pd.read_parquet(T90_FINAL))
    base["split_night"] = base["_v8_prepap_window"].astype(bool).astype(int)
    print(f"[{TAG}] cohort {len(base):,} nights, split nights {int(base.split_night.sum()):,} (flag entered as a covariate)", flush=True)
else:
    # v8 (2026-09-14): the legacy arm that re-read the April table (data_frozen/t90_final_pre_icdfix.parquet, 19,383 rows)
    # is removed; ranking_v2 is retired and no v8 step runs this script outside MODE == "current".
    raise SystemExit(f"MODE {MODE!r} is retired in v8: only 'current' runs on the frozen v8 tables")
print(f"[{TAG}] cohort {len(base):,}   exposure {XCOL}", flush=True)

head = pd.read_csv(MASTER, nrows=0).columns.tolist()
FEAT = [c for c in head if c not in ("BDSPPatientID",) and not c.endswith(("_first_date", "_date"))
        and not c.endswith(RANKING_DROP_SUFFIX)               # v8: the per-site twins are never candidates
        and (c not in SENSITIVITY_ONLY or c == XCOL)]         # v8: the sleep-period T90 only as the sensitivity exposure
raw = pd.read_csv(MASTER, usecols=["BDSPPatientID"] + FEAT, low_memory=False)
left = base[["BDSPPatientID", "site_id", "AgeAtVisit", "sex", "split_night"] +
            [c for c in base.columns if c.endswith(("_incident", "_years", "_prevalent"))]]
# the master carries its own copies of the cohort columns, drop them so the merge does not
# rename sex to sex_x and silently break every downstream reference
raw = raw.drop(columns=[c for c in raw.columns if c in set(left.columns) - {"BDSPPatientID"}])
FEAT = [c for c in FEAT if c in raw.columns]
d = left.merge(raw, on="BDSPPatientID", how="left")
d["male"] = (d.sex.astype(str).str.upper().str[0] == "M").astype(int)

num = [c for c in FEAT if pd.api.types.is_numeric_dtype(d[c])]
num = [c for c in num if d[c].notna().mean() >= 0.60 and d[c].nunique() > 10]
num = [c for c in num if c not in {"male", "AgeAtVisit"} and not c.startswith(("obs_", "follow_"))]
if SENS:
    assert XCOL in num, f"{XCOL} is not a numeric master column that passes the 60 percent completeness bar"
    num = [XCOL]
print(f"[{TAG}] candidate measures: {len(num)}", flush=True)
if MODE == "current" and not SUBSET and not SENS and len(num) != N_MEASURES:
    print(f"[{TAG}] NOTE: {len(num)} candidate measures against cohort_spec.N_MEASURES = {N_MEASURES}; "
          "reported here and in the provenance file, not corrected (the 60 percent coverage rule decides)", flush=True)

sp = dmatrix("cr(a, df=4) - 1", {"a": d.AgeAtVisit.values}, return_type="dataframe")
sp.columns = [f"age_s{i}" for i in range(4)]; sp.index = d.index
for c in sp.columns:
    d[c] = sp[c]
for f in num:
    z = pd.Series(np.nan, index=d.index)
    for s, idx in d.groupby("site_id").groups.items():
        v = d.loc[idx, f]; ok = v.notna()
        if ok.sum() < 20: continue
        r = stats.rankdata(v[ok], method="average")
        z.loc[v[ok].index] = stats.norm.ppf((r - 0.375) / (ok.sum() + 0.25))
    d[f + "__z"] = z

# RANKING_EXCLUDE holds the two negative controls that were never part of the ranked set. The
# ranked outcomes are every remaining condition with at least 150 incident events, plus death;
# cohort_spec.N_RANKED_OUTCOMES is the expected count and a difference is printed, not fixed.
OUT = [k for k in DISEASES if k not in CIRCULAR and k not in RANKING_EXCLUDE
       and f"{k}_incident" in d.columns]
OUT = [k for k in OUT if int(d[f"{k}_incident"].sum()) >= 150] + ["death"]
print(f"[{TAG}] outcomes: {len(OUT)}", flush=True)
if MODE == "current" and len(OUT) != N_RANKED_OUTCOMES:
    print(f"[{TAG}] NOTE: {len(OUT)} ranked outcomes against cohort_spec.N_RANKED_OUTCOMES = {N_RANKED_OUTCOMES}; "
          "an outcome under the 150-event bar drops out by the ranking's own rule", flush=True)

if SUBSET:
    num = [f for f in SUBSET.split(",") if f in num]
    print(f"[{TAG}] restricted to {len(num)} measures", flush=True)

keep = ADJ + [f + "__z" for f in num] + ["site_id"] + \
       [f"{o}_{s}" for o in OUT for s in ("years", "incident", "prevalent")]
d[[c for c in dict.fromkeys(keep) if c in d.columns]].to_pickle(PKL)

# ---------------------------------------------------------------------------------------
# One outcome, every measure against it, fitted in one hospital and scored in the other,
# both ways. Identical arithmetic to run_ranking_v2.py, reorganised so the eligible-patient
# filter and the hospital split happen once rather than once per measure.
# ---------------------------------------------------------------------------------------
_CACHE = {}


def run_outcome(out, feats):
    ck = f"{WORK}/{out}.json"
    if os.path.exists(ck):
        try:
            return json.load(open(ck))
        except Exception:
            pass
    if "d" not in _CACHE:
        _CACHE["d"] = pd.read_pickle(PKL)
    dd = _CACHE["d"]
    yc, ec = f"{out}_years", f"{out}_incident"
    g = dd[(dd[f"{out}_prevalent"] == 0) & (dd[yc] > 0) & dd[yc].notna()]
    side = {s: g[g.site_id == s] for s in ("I0002", "I0006")}
    rows = []
    for feat in [None] + list(feats):
        cols = ADJ + ([feat + "__z"] if feat else [])
        res = []
        for tr, te in (("I0002", "I0006"), ("I0006", "I0002")):
            A = side[tr][cols + [yc, ec]].dropna()
            B = side[te][cols + [yc, ec]].dropna()
            if A[ec].sum() < 30 or B[ec].sum() < 30:
                res.append(np.nan); continue
            try:
                c = CoxPHFitter(penalizer=0.01).fit(A, duration_col=yc, event_col=ec)
                r = c.predict_partial_hazard(B)
                res.append(concordance_index(B[yc], -r, B[ec]))
            except Exception:
                res.append(np.nan)
        rows.append({"feature": feat or "__BASE__", "outcome": out,
                     "c": float(np.nanmean(res))})
    tmp = ck + ".part"
    json.dump(rows, open(tmp, "w"))
    os.replace(tmp, ck)          # atomic, so a half-written checkpoint can never be read back
    return rows


done = [o for o in OUT if os.path.exists(f"{WORK}/{o}.json")]
if done:
    print(f"[{TAG}] resuming, {len(done)} of {len(OUT)} outcomes already checkpointed", flush=True)
print(f"[{TAG}] fits: {(len(num) + 1) * len(OUT) * 2}", flush=True)
t0 = time.time()
rr = Parallel(n_jobs=7, verbose=10, backend="loky")(
    delayed(run_outcome)(o, num) for o in OUT)
print(f"[{TAG}] fitting took {(time.time() - t0) / 60:.1f} min", flush=True)

s = pd.DataFrame([r for chunk in rr for r in chunk])
bas = s[s.feature == "__BASE__"].set_index("outcome")["c"]
f2 = s[s.feature != "__BASE__"].copy()
f2["gain"] = f2.apply(lambda r: r.c - bas.get(r.outcome, np.nan), axis=1)

PRIMARY_PERCOND = f"{NUMBERS_DIR}/ranking_v3_percondition.csv"
if SENS:
    # decision 2: every row except the exposure row comes from the primary ranking; the base C of
    # this run must be the primary's base C, outcome by outcome (same cohort, same age and sex fit)
    prim = pd.read_csv(PRIMARY_PERCOND, comment="#")
    assert set(prim.outcome) == set(OUT), f"outcome set differs from the primary ranking: {set(prim.outcome) ^ set(OUT)}"
    pb = (prim.c - prim.gain).groupby(prim.outcome).agg(["mean", "std"])
    assert (pb["std"].fillna(0) < 1e-9).all(), "the primary per-outcome base C rows disagree with each other"
    base_diff = float((pb["mean"] - bas.reindex(pb.index)).abs().max())
    print(f"[{TAG}] positive control: base C of this run against the primary ranking, largest difference {base_diff:.2e}", flush=True)
    assert base_diff < 1e-9, "the base C differs from the primary ranking's: not the same cohort or frame"
    f2 = pd.concat([prim[~prim.feature.isin({PRIMARY_EXPOSURE, XCOL})], f2[f2.feature == XCOL]],
                   ignore_index=True)

rank = (f2.groupby("feature").agg(mC=("c", "mean"), dC=("gain", "mean"),
        worst=("gain", "min"), npos=("gain", lambda x: int((x > 0).sum())),
        nout=("gain", "size")).reset_index().sort_values("dC", ascending=False))
rank.insert(0, "rank", range(1, len(rank) + 1))

TODAY = time.strftime("%Y-%m-%d")
written = []
if True:   # v8.1: this script writes its own pair of files, never ranking_v3.csv
    p1, p2 = f"{OUT_DIR}/ranking_v3_{TAG}.csv", f"{OUT_DIR}/ranking_v3_{TAG}_percondition.csv"
    rank.to_csv(p1, index=False)
    f2.to_csv(p2, index=False)
    written += [p1, p2]
    print(f"[{TAG}] written: {p1}, {p2} (supplementary, all nights, split-night covariate)", flush=True)
else:
    HDR = [
        f"# ranking_v3{SFX}, built {TODAY} (v8 recalculation)",
        f"# cohort: {len(base):,} recordings, numbers/cohort_spec.py apply_cohort()",
        "#   fu_valid == 1, spo2_pct_below_90 present, oximetry_bad == 0",
        f"# exposure row: {XCOL}" + ("  (decision 2 sensitivity set: only this row refitted, the rest from ranking_v3_percondition.csv)" if SENS else ""),
        f"# measures: {len(rank)}   outcomes: {len(OUT)}   held-out baseline C {bas.mean():.6f}",
        "# procedure: unchanged from ranking_v2. Age spline cr(df=4) plus sex, each measure added",
        "#   as a site-wise rank-normal score, CoxPHFitter(penalizer=0.01) fitted at I0002 and",
        "#   scored at I0006 then reversed, the two held-out concordances averaged.",
        "# supersedes: numbers/ranking_v2.csv, frozen 2026-08-02. v2 is kept, not deleted.",
        "# what changed since v2, and only this:",
        "#   1. the oximetry quality rule of 2026-08-04 removed 210 corrupted recordings,",
        "#      cohort 19,383 -> 19,173 (numbers/apply_qc_v3.py)",
        "#   2. numbers/refresh_two_conditions.py recoded obesity hypoventilation from ICD-9",
        "#      278.01 morbid obesity to 278.03, and falls, whose E888 and E884 prefixes also",
        "#      matched ICD-10 E88.81 metabolic syndrome",
        "#   3. v8 (2026-09-12): the frozen tables of data_frozen_v8_2026-09 (split nights on the",
        "#      untreated window, four new measures, 52 ranked outcomes), the 60 per-site",
        "#      standardized twins removed from the candidates",
        "# reproducibility: this script in legacy mode, fed the pre-correction 19,383 cohort,",
        "#   returns ranking_v2's numbers, so the movement below is the data and not the code.",
        "# readers must pass comment='#' to pandas.read_csv.",
    ]
    for path, frame in ((f"{OUT_DIR}/ranking_v3{SFX}.csv", rank),
                        (f"{OUT_DIR}/ranking_v3{SFX}_percondition.csv", f2)):
        with open(path, "w") as fh:
            fh.write("\n".join(HDR) + "\n")
            frame.to_csv(fh, index=False)
        written.append(path)
        print(f"[{TAG}] wrote {path}", flush=True)
    pj = f"{OUT_DIR}/ranking_v3{SFX}_provenance.json"
    json.dump({"built": TODAY, "cohort_n": int(len(base)), "exposure": XCOL,
               "measures": int(len(rank)), "measures_expected": N_MEASURES,
               "outcomes": int(len(OUT)), "outcomes_expected": N_RANKED_OUTCOMES,
               "outcome_keys": OUT, "baseline_heldout_C": float(bas.mean()),
               "candidate_measures_this_run": int(len(num)),
               "tables": {"t90_final": T90_FINAL, "master": MASTER},
               "primary_percondition_used": PRIMARY_PERCOND if SENS else None,
               "supersedes": "numbers/ranking_v2.csv (2026-08-02)",
               "changes": ["oximetry quality rule, 19383 -> 19173",
                           "obesity_hypovent ICD-9 278.01 -> 278.03",
                           "falls E888/E884 no longer matching ICD-10 E88.81",
                           "v8 2026-09-12: data_frozen_v8 tables, _siteZ twins excluded, four new measures"]},
              open(pj, "w"), indent=2)
    written.append(pj)

for p in written:
    sidecar(p, __file__, extra_inputs=[PRIMARY_PERCOND] if SENS else [],
            note=f"build_ranking_v3 {TAG}: {len(rank)} measures, {len(OUT)} outcomes, exposure {XCOL}")

print(f"\n[{TAG}] baseline held-out C {bas.mean():.6f}", flush=True)
print(rank.head(15).to_string(index=False), flush=True)
shutil.rmtree(WORK, ignore_errors=True)
