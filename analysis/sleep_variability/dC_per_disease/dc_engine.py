"""
Per-disease discrimination gain for T90, with a stratified nonparametric bootstrap.

The arithmetic is the published ranking's, transcribed from
Sleep_Variability_2026-08/dC_side_analyses/common.py (which is itself
T90_Manuscript/numbers/build_ranking_v3.py). This module only reorganises it so that

  1. the two cross-fit directions are reported separately instead of being averaged away, and
  2. the whole fit-and-evaluate pass is cheap enough to repeat several hundred times.

Nothing here may deviate from common.py's arithmetic. p00_positive_control.py gates the folder
on this engine reproducing common.heldout's per-outcome held-out C to 1e-12 and on
common.heldout itself reproducing numbers/ranking_v3_percondition.csv to 1e-9.

Two age bases are carried, deliberately.

  full   the published four cr(df=4) columns, singular, held together by the 0.01 ridge.
         This is the basis that returns the paper's 0.02033875497249791 and it is used ONLY
         for the positive control and for a cross-reference column.
  drop   the same spline with one basis column removed, which is the rank-safe basis and the
         one this analysis reports. Alen's standing rule for every T90 fit.

The bootstrap
  Patients are resampled with replacement WITHIN hospital, so both arms of the cross-fit keep
  their site sizes and the held-out design survives every replicate. The entire procedure is
  then repeated on the resample, not just the final fit: the age spline knots are re-placed on
  the resampled ages and the within-site rank inverse normal of T90 is re-ranked on the
  resampled values, because both are estimated quantities in the published pipeline.
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
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from patsy import dmatrix
from scipy import stats
import sys as _sys_v8; _sys_v8.path.insert(0, paths.ANALYSIS_DIR)  # v8 sweep 2026-09-12
from cohort_spec import COHORT_N, N_RANKED_OUTCOMES  # v8 sweep 2026-09-12

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "_work")
SIDE = f"{paths.SV_CODE_DIR}/dC_side_analyses"
T90ROOT = paths.T90_ROOT
os.makedirs(WORK, exist_ok=True)

FEAT = "spo2_pct_below_90"
PEN = 0.01                       # the ranking's penalizer
MIN_FOLD_EVENTS = 30             # the ranking's per-fold floor
SITES = ("I0002", "I0006")       # 0, 1
FOLDS = ((0, 1), (1, 0))         # fit at I0002 score at I0006, then the reverse
ARRAYS = os.path.join(WORK, "arrays.pkl")


# ------------------------------------------------------------------ the frame, once
def build_arrays(force: bool = False) -> dict:
    """
    Freeze the ranking's analysis frame into plain numpy, one entry per outcome.

    The frame itself is common.build_frame(), unmodified: the 19,173 cohort from
    data_frozen/t90_final.parquet through numbers/cohort_spec.apply_cohort, joined to
    master_cohort_v2.csv for the measure column.
    """
    if os.path.exists(ARRAYS) and not force:
        return pd.read_pickle(ARRAYS)

    sys.path.insert(0, SIDE)
    from common import build_frame, ranking_outcomes          # noqa: E402

    d = build_frame(extra_master_cols=[FEAT])
    outs = ranking_outcomes(d)
    assert len(outs) == N_RANKED_OUTCOMES, f"outcome set is {len(outs)}, expected {N_RANKED_OUTCOMES}"
    assert len(d) == COHORT_N, len(d)

    site = np.full(len(d), -1, np.int8)
    for i, s in enumerate(SITES):
        site[(d.site_id == s).values] = i
    assert (site >= 0).all(), "a cohort row sits outside the two ranking sites"

    A = {
        "outcomes": outs,
        "site": site,
        "age": d.AgeAtVisit.to_numpy(float),
        "male": d.male.to_numpy(float),
        "t90": d[FEAT].to_numpy(float),
        "years": {o: d[f"{o}_years"].to_numpy(float) for o in outs},
        "event": {o: d[f"{o}_incident"].to_numpy(float) for o in outs},
        "prev": {o: d[f"{o}_prevalent"].to_numpy(float) for o in outs},
    }
    pd.to_pickle(A, ARRAYS)
    return A


# ------------------------------------------------------------------ the transforms
def age_basis(age: np.ndarray) -> np.ndarray:
    """The ranking's age spline: cr(df=4) with no intercept, four columns."""
    sp = dmatrix("cr(a, df=4) - 1", {"a": age}, return_type="dataframe")
    return np.asarray(sp, float)


def site_rint(v: np.ndarray, site: np.ndarray) -> np.ndarray:
    """
    common.site_rint in numpy. Blom rank inverse normal inside each hospital, missing values
    left missing, a site with fewer than 20 non-missing values left entirely missing.
    """
    z = np.full(len(v), np.nan)
    for s in (0, 1):
        m = (site == s) & ~np.isnan(v)
        n = int(m.sum())
        if n < 20:
            continue
        r = stats.rankdata(v[m], method="average")
        z[m] = stats.norm.ppf((r - 0.375) / (n + 0.25))
    return z


# ------------------------------------------------------------------ one held-out fit
def _c(Xtr, ttr, etr, Xte, tte, ete, names) -> float:
    """common._one's body: fit at one site with penalizer 0.01, score at the other."""
    tr = pd.DataFrame(Xtr, columns=names)
    tr["_t"], tr["_e"] = ttr, etr
    te = pd.DataFrame(Xte, columns=names)
    try:
        cph = CoxPHFitter(penalizer=PEN).fit(tr, duration_col="_t", event_col="_e")
        r = cph.predict_partial_hazard(te)
        return float(concordance_index(tte, -np.asarray(r, float), ete))
    except Exception:
        return np.nan


# ------------------------------------------------------------------ one whole pass
def one_pass(A: dict, seed=None, basis: str = "drop", lag: float = 0.0) -> np.ndarray:
    """
    One complete fit-and-evaluate pass over all 48 outcomes and both cross-fit directions.

    seed   None runs the observed data. An integer resamples patients with replacement inside
           each hospital before anything else is computed.
    basis  "drop" keeps age_s1..age_s3 (the rank-safe basis this analysis reports).
           "full" keeps all four published columns.
    lag    landmark in years. 0.0 is the paper's ranking. For lag L the rule is the one the
           published lag ladder uses, transcribed from New_Figures/scripts/lagladder.py and
           recorded in numbers/lag_ladder.json: everyone whose event or censoring arrived at
           or before L leaves the risk set, and the clock restarts at L. The event flag is not
           touched, because every event that survives the filter happens after L. Prevalent
           exclusion, the measure, the age spline and the site design are unchanged.

    Returns (48, 4): C(age+sex) dir1, C(age+sex+T90) dir1, C(age+sex) dir2, C(+T90) dir2,
    where dir1 fits at I0002 and scores at I0006. NaN where a fold misses the 30-event floor.
    """
    outs = A["outcomes"]
    site0 = A["site"]

    if seed is None:
        idx = np.arange(len(site0))
    else:
        rng = np.random.default_rng(seed)
        parts = []
        for s in (0, 1):
            pool = np.flatnonzero(site0 == s)
            parts.append(rng.choice(pool, size=len(pool), replace=True))
        idx = np.concatenate(parts)

    site = site0[idx]
    sp = age_basis(A["age"][idx])
    keep = [0, 1, 2, 3] if basis == "full" else [1, 2, 3]
    names = [f"age_s{i}" for i in keep] + ["male"]
    Xb = np.column_stack([sp[:, keep], A["male"][idx]])
    z = site_rint(A["t90"][idx], site)
    Xf = np.column_stack([Xb, z])
    names_f = names + ["t90_z"]

    ok_cov = ~np.isnan(Xb).any(1)                    # common._one's dropna over the used cols
    ok_cov_f = ok_cov & ~np.isnan(z)

    out = np.full((len(outs), 4), np.nan)
    for i, o in enumerate(outs):
        y = A["years"][o][idx]
        e = A["event"][o][idx]
        p = A["prev"][o][idx]
        elig = (p == 0) & (y > 0) & ~np.isnan(y) & ~np.isnan(e)
        if lag > 0:
            elig = elig & (y > lag)
            y = y - lag
        for k, (tr, te) in enumerate(FOLDS):
            for j, (msk, X, nm) in enumerate(((ok_cov, Xb, names), (ok_cov_f, Xf, names_f))):
                a = elig & msk & (site == tr)
                b = elig & msk & (site == te)
                if e[a].sum() < MIN_FOLD_EVENTS or e[b].sum() < MIN_FOLD_EVENTS:
                    continue
                out[i, 2 * k + j] = _c(X[a], y[a], e[a], X[b], y[b], e[b], nm)
    return out


def dc_from_pass(P: np.ndarray) -> np.ndarray:
    """Per-outcome dC: the gain in each direction, then the nanmean of the two, as the paper."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return np.nanmean(np.column_stack([P[:, 1] - P[:, 0], P[:, 3] - P[:, 2]]), axis=1)


# ------------------------------------------------------------------ event counts
def counts(A: dict, lag: float = 0.0) -> pd.DataFrame:
    """Incident events and patients at risk on the eligible set the folds actually see."""
    rows = []
    for o in A["outcomes"]:
        y, e, p, s = A["years"][o], A["event"][o], A["prev"][o], A["site"]
        elig = (p == 0) & (y > 0) & ~np.isnan(y) & ~np.isnan(e)
        if lag > 0:
            elig = elig & (y > lag)
        rows.append({
            "outcome": o,
            "events": int(e[elig].sum()),
            "n_at_risk": int(elig.sum()),
            "events_I0002": int(e[elig & (s == 0)].sum()),
            "events_I0006": int(e[elig & (s == 1)].sum()),
            "n_I0002": int((elig & (s == 0)).sum()),
            "n_I0006": int((elig & (s == 1)).sum()),
        })
    return pd.DataFrame(rows)
