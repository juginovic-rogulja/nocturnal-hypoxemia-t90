"""
The published per-disease dC engine, with one thing switched: the exposure.

Everything else is the paper's. dc_engine.one_pass (Sleep_Variability_2026-08/dC_per_disease/
dc_engine.py, itself a transcription of T90_Manuscript/numbers/build_ranking_v3.py) is imported
unmodified and used for the frame, the age basis, the site rank inverse normal, the held-out
Cox fit, the fold floor and the landmark rule. one_pass_x below is a line-for-line copy of
dc_engine.one_pass with a single branch added where the exposure column is built:

    exposure="rint"   the published column, site-wise rank inverse normal of spo2_pct_below_90
    exposure="gt10"   1.0 when spo2_pct_below_90 > 10, else 0.0, entered raw

identity_guard() asserts that one_pass_x with exposure="rint" returns dc_engine.one_pass's
array bit for bit, on the observed data and on resamples, so the copy cannot drift.

The T90 column is spo2_pct_below_90 in
sleep-outcome-sandbox/feature_pipeline/outputs/master_cohort_v2.csv, in percent of the
recording. The > 10 rule is the paper's own top oxygen band (dC_side_analyses/common.py BANDS,
membership lo < v <= hi, so the top band is v > 10).
"""
from __future__ import annotations
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../..")))  # repository root, where paths.py lives
import paths

import os
import sys
import warnings

warnings.filterwarnings("ignore")

import numpy as np

PUB = f"{paths.SV_CODE_DIR}/dC_per_disease"
sys.path.insert(0, PUB)

import dc_engine as E                                                     # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "_work")
os.makedirs(WORK, exist_ok=True)

GT10_THRESHOLD = 10.0
EXPOSURES = ("rint", "gt10")
COLNAME = {"rint": "t90_z", "gt10": "t90_gt10"}


def exposure_column(t90: np.ndarray, site: np.ndarray, exposure: str) -> np.ndarray:
    if exposure == "rint":
        return E.site_rint(t90, site)
    if exposure == "gt10":
        z = np.where(np.isnan(t90), np.nan, (t90 > GT10_THRESHOLD).astype(float))
        return z
    raise ValueError(exposure)


def one_pass_x(A: dict, seed=None, basis: str = "drop", lag: float = 0.0,
               exposure: str = "rint") -> np.ndarray:
    """dc_engine.one_pass, with the exposure column chosen by `exposure`. Nothing else moves."""
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
    sp = E.age_basis(A["age"][idx])
    keep = [0, 1, 2, 3] if basis == "full" else [1, 2, 3]
    names = [f"age_s{i}" for i in keep] + ["male"]
    Xb = np.column_stack([sp[:, keep], A["male"][idx]])
    z = exposure_column(A["t90"][idx], site, exposure)
    Xf = np.column_stack([Xb, z])
    names_f = names + [COLNAME[exposure]]

    ok_cov = ~np.isnan(Xb).any(1)
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
        for k, (tr, te) in enumerate(E.FOLDS):
            for j, (msk, X, nm) in enumerate(((ok_cov, Xb, names), (ok_cov_f, Xf, names_f))):
                a = elig & msk & (site == tr)
                b = elig & msk & (site == te)
                if e[a].sum() < E.MIN_FOLD_EVENTS or e[b].sum() < E.MIN_FOLD_EVENTS:
                    continue
                out[i, 2 * k + j] = E._c(X[a], y[a], e[a], X[b], y[b], e[b], nm)
    return out


def exposed_counts(A: dict, lag: float = 0.0) -> "object":
    """Published counts plus the number above 10% on the same eligible set the folds see."""
    import pandas as pd
    base = E.counts(A, lag=lag).set_index("outcome")
    gt = A["t90"] > GT10_THRESHOLD
    rows = []
    for o in A["outcomes"]:
        y, e, p, s = A["years"][o], A["event"][o], A["prev"][o], A["site"]
        elig = (p == 0) & (y > 0) & ~np.isnan(y) & ~np.isnan(e)
        if lag > 0:
            elig = elig & (y > lag)
        rows.append({
            "outcome": o,
            "n_gt10": int((elig & gt).sum()),
            "n_le10": int((elig & ~gt).sum()),
            "pct_gt10": float(100 * (elig & gt).sum() / max(elig.sum(), 1)),
            "events_gt10": int(e[elig & gt].sum()),
            "events_le10": int(e[elig & ~gt].sum()),
            "n_gt10_I0002": int((elig & gt & (s == 0)).sum()),
            "n_gt10_I0006": int((elig & gt & (s == 1)).sum()),
            "events_gt10_I0002": int(e[elig & gt & (s == 0)].sum()),
            "events_gt10_I0006": int(e[elig & gt & (s == 1)].sum()),
        })
    return base.join(pd.DataFrame(rows).set_index("outcome"))


def identity_guard(A: dict, seeds=(None, 20260820, 20260999), verbose=True) -> bool:
    """one_pass_x(exposure='rint') must equal dc_engine.one_pass exactly, or the copy drifted."""
    ok = True
    for basis, lag in (("full", 0.0), ("drop", 0.0), ("drop", 2.0)):
        for sd in seeds:
            a = E.one_pass(A, seed=sd, basis=basis, lag=lag)
            b = one_pass_x(A, seed=sd, basis=basis, lag=lag, exposure="rint")
            same = np.array_equal(a, b, equal_nan=True)
            ok = ok and same
            if verbose:
                print(f"    basis={basis:<4} lag={lag:.0f} seed={str(sd):<8} "
                      f"{'identical' if same else 'DIFFERS'}")
    return ok


# =======================================================================================
# 2026-09-04, second specification: the contrast between the two extreme bands.
#
# T90 above 10% against T90 at or below 1%, with the intermediate 1 to 10 group dropped
# from the cohort entirely rather than folded into the reference. The band edges are not
# retyped here. They are imported from the paper's own definition,
# Sleep_Variability_2026-08/dC_side_analyses/common.py BANDS, which is itself
# T90_Manuscript/numbers/run_all_v2.py line 37 and run_treatment_v2.py line 70:
#
#     BANDS = [(-1, 1, "0-1%"), (1, 5, "1-5%"), (5, 10, "5-10%"), (10, 1e9, ">10%")]
#     membership is  lo < v <= hi
#
# so the reference band is  -1 < v <= 1  and the top band is  v > 10. A T90 of exactly
# 1.0 belongs to the REFERENCE band, and a T90 of exactly 10.0 belongs to the 5-10 band
# and is therefore DROPPED. These are the bands of Figure 2b and Supplementary Table 5,
# and the manuscript writes them as "1% or less" and "over 10%".
# =======================================================================================
SIDE = f"{paths.SV_CODE_DIR}/dC_side_analyses"
sys.path.insert(0, SIDE)
from common import BANDS                                                  # noqa: E402

REF_LO, REF_HI, REF_LAB = BANDS[0]
TOP_LO, TOP_HI, TOP_LAB = BANDS[-1]
assert (REF_HI, TOP_LO) == (1, 10), f"band edges moved: {BANDS}"
assert TOP_LO == GT10_THRESHOLD, "the top band edge and the gt10 threshold disagree"


def band_masks(A: dict):
    """Reference and top band membership under the paper's lo < v <= hi rule."""
    t = A["t90"]
    ref = (t > REF_LO) & (t <= REF_HI)
    top = (t > TOP_LO) & (t <= TOP_HI)
    assert not (ref & top).any(), "a recording landed in both bands"
    return ref, top


def restrict_arrays(A: dict, mask: np.ndarray) -> dict:
    """A copy of the frozen analysis frame holding only the rows in `mask`."""
    assert mask.dtype == bool and len(mask) == len(A["site"])
    B = {"outcomes": A["outcomes"]}
    for k in ("site", "age", "male", "t90"):
        B[k] = A[k][mask]
    for k in ("years", "event", "prev"):
        B[k] = {o: A[k][o][mask] for o in A["outcomes"]}
    return B


def extremes_arrays(A: dict) -> dict:
    """The restricted cohort: the two extreme bands only, the middle dropped before any fit."""
    ref, top = band_masks(A)
    B = restrict_arrays(A, ref | top)
    t = B["t90"]
    assert ((t <= REF_HI) | (t > TOP_LO)).all(), "an intermediate recording survived"
    return B


def band_counts(A: dict, lag: float = 0.0) -> "object":
    """Per outcome, on the eligible set the folds see: n and events in each of the two bands."""
    import pandas as pd
    ref, top = band_masks(A)
    rows = []
    for o in A["outcomes"]:
        y, e, p, s = A["years"][o], A["event"][o], A["prev"][o], A["site"]
        elig = (p == 0) & (y > 0) & ~np.isnan(y) & ~np.isnan(e)
        if lag > 0:
            elig = elig & (y > lag)
        rows.append({
            "outcome": o,
            "n_restricted": int((elig & (ref | top)).sum()),
            "events_restricted": int(e[elig & (ref | top)].sum()),
            "n_le1": int((elig & ref).sum()),
            "events_le1": int(e[elig & ref].sum()),
            "n_gt10": int((elig & top).sum()),
            "events_gt10": int(e[elig & top].sum()),
            "n_dropped_middle": int((elig & ~(ref | top)).sum()),
            "events_dropped_middle": int(e[elig & ~(ref | top)].sum()),
            "n_le1_I0002": int((elig & ref & (s == 0)).sum()),
            "n_le1_I0006": int((elig & ref & (s == 1)).sum()),
            "n_gt10_I0002": int((elig & top & (s == 0)).sum()),
            "n_gt10_I0006": int((elig & top & (s == 1)).sum()),
            "events_restricted_I0002": int(e[elig & (ref | top) & (s == 0)].sum()),
            "events_restricted_I0006": int(e[elig & (ref | top) & (s == 1)].sum()),
        })
    return pd.DataFrame(rows).set_index("outcome")
