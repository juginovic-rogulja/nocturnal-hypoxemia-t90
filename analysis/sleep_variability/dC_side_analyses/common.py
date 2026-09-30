"""
The published ranking pipeline, re-implemented once so every side analysis in this folder
uses identical arithmetic.

Conventions are copied from T90_Manuscript/numbers/build_ranking_v3.py (which is itself
run_ranking_v2.py reorganised, with the fitting untouched). Nothing here may deviate from
that file without saying so out loud, because s00_positive_control.py gates the whole folder
on reproducing ranking_v3.csv's dC for T90 to 1e-9.

What dC is, exactly
  Base model      age as a 4-column natural cubic spline, cr(df=4) with no intercept, plus sex.
  Feature         added as ONE extra column, the within-site rank inverse normal of the measure.
  Fit             lifelines CoxPHFitter(penalizer=0.01), no strata. Site separation is done by
                  the fold, not by a strata term.
  Folds           fit at I0002 and score at I0006, then the reverse. Two held-out concordances,
                  averaged with nanmean. A fold is dropped when either side has fewer than 30
                  events.
  Eligibility     per outcome, prevalent == 0 and follow-up > 0 and follow-up not missing.
  Row set         dropna over exactly the columns the model uses, so the base model and the
                  feature model are NOT on a common row set. That is what the published code
                  does and it is replicated rather than corrected.
  dC              mean over the 48 outcomes of (feature C minus base C), per outcome.

The age-spline rank trap
  The four cr(df=4) columns sum to one and are singular. The published ranking survives that
  only because penalizer=0.01 regularises it. Every fit in this folder therefore keeps the
  published 4-column basis at penalizer 0.01, and s00 additionally refits with one column
  dropped to prove the published basis is not silently fitting nothing. Any new spline basis
  introduced by a side analysis (the T90 spline in s04) has one column dropped, per the rule.
"""
from __future__ import annotations
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths

import os
import sys
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from patsy import dmatrix
from scipy import stats

T90ROOT = paths.T90_ROOT
HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "_work")
os.makedirs(WORK, exist_ok=True)

sys.path.insert(0, paths.ANALYSIS_DIR)
from cohort_spec import apply_cohort, COHORT_N, ranked_outcome_keys, all_nights_cohort, MIN_RANKED_EVENTS   # noqa: E402  (v8.1)
from disease_definitions import (DISEASES, RANKING_EXCLUDE,          # noqa: E402
                                 ORGAN_GROUP, CVD_COMPONENTS)

MASTER = (f"{paths.TABLES_DIR}/"
          "master_cohort.csv")
PARQUET = f"{paths.TABLES_DIR}/t90_final.parquet"
SHARDS = f"{paths.V8_ROOT}/X4_groupP/cpap_stage/out*.csv"   # v8.2 (2026-09-15): re-extracted, collapse-masked shards (X4_groupP)

CIRCULAR = {"insomnia", "rls_plmd", "nocturia", "epilepsy"}
ADJ = [f"age_s{i}" for i in range(4)] + ["male"]
PEN = 0.01                      # the ranking's penalizer
MIN_FOLD_EVENTS = 30            # the ranking's per-fold floor
FOLDS = (("I0002", "I0006"), ("I0006", "I0002"))

# The paper's 4 oxygen bands, from numbers/run_all_v2.py and numbers/run_treatment_v2.py.
# Membership rule is lo < v <= hi, so a T90 of exactly 1.0 sits in the first band.
BANDS = [(-1, 1, "0-1%"), (1, 5, "1-5%"), (5, 10, "5-10%"), (10, 1e9, ">10%")]


# ------------------------------------------------------------------ transforms
def rint(x) -> np.ndarray:
    """Rank inverse normal, Blom. Identical to the ranking's inline version."""
    x = np.asarray(x, float)
    r = stats.rankdata(x, method="average")
    return stats.norm.ppf((r - 0.375) / (len(r) + 0.25))


def site_rint(frame: pd.DataFrame, col: str, sitecol: str = "site_id") -> pd.Series:
    """
    The ranking's within-site rank inverse normal, transcribed from build_ranking_v3.py.
    Missing values stay missing and do not enter the ranking of their own site. Sites with
    fewer than 20 non-missing values are left entirely missing, as in the published code.
    """
    z = pd.Series(np.nan, index=frame.index)
    for _s, idx in frame.groupby(sitecol).groups.items():
        v = frame.loc[idx, col]
        ok = v.notna()
        if ok.sum() < 20:
            continue
        r = stats.rankdata(v[ok], method="average")
        z.loc[v[ok].index] = stats.norm.ppf((r - 0.375) / (ok.sum() + 0.25))
    return z


def site_z(frame: pd.DataFrame, col: str, sitecol: str = "site_id") -> pd.Series:
    """Within-site standardisation, used only where a side analysis needs to hold the site
    alignment of the published transform while changing its shape."""
    g = frame.groupby(sitecol)[col]
    return (frame[col] - g.transform("mean")) / g.transform("std")


def band_index(v) -> np.ndarray:
    v = np.asarray(v, float)
    out = np.full(len(v), np.nan)
    for i, (lo, hi, _lab) in enumerate(BANDS):
        out[(v > lo) & (v <= hi)] = i
    return out


# ------------------------------------------------------------------ outcome sets
def ranking_outcomes(frame: pd.DataFrame) -> list[str]:
    """
    The ranking's 48. Transcribed from build_ranking_v3.py: every condition that is neither
    circular nor held out of the ranking, with at least 150 incident events in the full
    cohort, plus death.
    """
    out = [k for k in DISEASES if k not in CIRCULAR and k not in RANKING_EXCLUDE
           and f"{k}_incident" in frame.columns]
    out = ranked_outcome_keys(pd.read_parquet(PARQUET), out) + ["death"]   # v8.1: counted on the all-nights cohort, not on the frame passed in
    return out


def combo_outcomes(frame: pd.DataFrame) -> list[str]:
    """The 54 of run_combo_v2.py: every _incident column with 150 or more events."""
    out = [k[:-9] for k in frame.columns if k.endswith("_incident")]
    full = all_nights_cohort(pd.read_parquet(PARQUET))          # v8.1: counted on the all-nights cohort
    return [k for k in out if f"{k}_incident" in full.columns and int(full[f"{k}_incident"].sum()) >= MIN_RANKED_EVENTS]


# ------------------------------------------------------------------ the frame
def build_frame(extra_master_cols=(), cache="frame_main.pkl") -> pd.DataFrame:
    """
    The ranking's analysis frame: the 19,173 cohort joined to master_cohort_v2.csv.

    The measure columns come from the master, not from the frozen parquet, because the
    published ranking merges the master and drops only the columns that collide with the
    cohort key block. spo2_pct_below_90 is therefore the master's copy in every ranking
    number, and it is the master's copy here.
    """
    path = os.path.join(WORK, cache)
    want = list(dict.fromkeys(list(extra_master_cols)))
    if os.path.exists(path):
        d = pd.read_pickle(path)
        if all(c in d.columns for c in want) and len(d) == COHORT_N:   # v8.1c: a cached frame from the other cohort size is rebuilt, never reused
            return d

    base = apply_cohort(pd.read_parquet(PARQUET))
    assert len(base) == COHORT_N, len(base)

    head = pd.read_csv(MASTER, nrows=0).columns.tolist()
    missing = [c for c in want if c not in head]
    assert not missing, f"master_cohort_v2.csv has no column for {missing}"

    left = base[["BDSPPatientID", "site_id", "AgeAtVisit", "sex"] +
                [c for c in base.columns
                 if c.endswith(("_incident", "_years", "_prevalent"))]]
    raw = pd.read_csv(MASTER, usecols=["BDSPPatientID"] + want, low_memory=False)
    # same collision guard as the published merge
    raw = raw.drop(columns=[c for c in raw.columns
                            if c in set(left.columns) - {"BDSPPatientID"}])
    d = left.merge(raw, on="BDSPPatientID", how="left")
    d["male"] = (d.sex.astype(str).str.upper().str[0] == "M").astype(int)

    sp = dmatrix("cr(a, df=4) - 1", {"a": d.AgeAtVisit.values}, return_type="dataframe")
    sp.columns = [f"age_s{i}" for i in range(4)]
    sp.index = d.index
    for c in sp.columns:
        d[c] = sp[c]

    d.to_pickle(path)
    return d


# ------------------------------------------------------------------ the fitter
def _one(pkl, cols, out, tr, te, penalizer, adj):
    dd = pd.read_pickle(pkl)
    yc, ec = f"{out}_years", f"{out}_incident"
    g = dd[(dd[f"{out}_prevalent"] == 0) & (dd[yc] > 0) & dd[yc].notna()]
    use = list(adj) + list(cols)
    res = []
    for _tr, _te in ((tr, te),):
        A = g[g.site_id == _tr][use + [yc, ec]].dropna()
        B = g[g.site_id == _te][use + [yc, ec]].dropna()
        if A[ec].sum() < MIN_FOLD_EVENTS or B[ec].sum() < MIN_FOLD_EVENTS:
            res.append(np.nan)
            continue
        try:
            c = CoxPHFitter(penalizer=penalizer).fit(A, duration_col=yc, event_col=ec)
            r = c.predict_partial_hazard(B)
            res.append(concordance_index(B[yc], -r, B[ec]))
        except Exception:
            res.append(np.nan)
    return res[0]


def heldout(specs: dict, outcomes, frame: pd.DataFrame, penalizer: float = PEN,
            adj=ADJ, pkl_name="fit_frame.pkl", n_jobs: int = 7) -> pd.DataFrame:
    """
    Held-out concordance for every (spec, outcome), on the published two-fold site design.

    specs   {name: [column names]}. Must contain "__base__": [] for gains to be computed.
    Returns a long frame with columns spec, outcome, c, and gain where a base exists.
    """
    pkl = os.path.join(WORK, pkl_name)
    keep = list(dict.fromkeys(
        list(adj) + [c for cols in specs.values() for c in cols] + ["site_id"] +
        [f"{o}_{s}" for o in outcomes for s in ("years", "incident", "prevalent")]))
    frame[[c for c in keep if c in frame.columns]].to_pickle(pkl)

    jobs = [(k, o, tr, te) for k in specs for o in outcomes for tr, te in FOLDS]
    res = Parallel(n_jobs=n_jobs, verbose=0, backend="loky")(
        delayed(_one)(pkl, specs[k], o, tr, te, penalizer, adj) for k, o, tr, te in jobs)

    acc = {}
    for (k, o, _tr, _te), v in zip(jobs, res):
        acc.setdefault((k, o), []).append(v)
    rows = [{"spec": k, "outcome": o, "c": float(np.nanmean(v))}
            for (k, o), v in acc.items()]
    df = pd.DataFrame(rows)
    if "__base__" in specs:
        bas = df[df.spec == "__base__"].set_index("outcome")["c"]
        df["gain"] = df.apply(lambda r: r.c - bas.get(r.outcome, np.nan), axis=1)
    return df


def summarise(df: pd.DataFrame) -> pd.DataFrame:
    """spec-level dC, the ranking's own aggregation: the plain mean of the per-outcome gains."""
    f = df[df.spec != "__base__"]
    s = (f.groupby("spec")
           .agg(mC=("c", "mean"), dC=("gain", "mean"), worst=("gain", "min"),
                best=("gain", "max"), npos=("gain", lambda x: int((x > 0).sum())),
                nout=("gain", "size"))
           .reset_index().sort_values("dC", ascending=False))
    return s


# ------------------------------------------------------------------ organ grouping
# The five organ groups the task calls the plausible oxygen pathway, expressed in the
# vocabulary ORGAN_GROUP actually uses. ORGAN_GROUP says "Cardiac" where the task says
# cardiovascular and "Kidney" where it says renal.
OXYGEN_PATHWAY_GROUPS = ["Cardiac", "Respiratory", "Metabolic", "Kidney"]


def organ_of(key: str) -> str:
    if key == "death":
        return "Death"
    return ORGAN_GROUP.get(key, "UNGROUPED")
