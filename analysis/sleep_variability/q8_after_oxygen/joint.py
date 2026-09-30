"""
Q8, part 2. Joint categorical analysis: does nocturnal oxygenation still separate patients
INSIDE every level of the apnea-hypopnea index, of sleep duration and of sleep efficiency?

Part 1 answered the prediction question on a continuous scale and found that the three
conventional measures add nothing once T90 is in the model. That is a statement about
discrimination averaged over a population. It does not tell a clinician whether the oxygen
signal survives in the patients they are actually looking at: the person with no apnea, the
person who slept badly, the person who slept eight hours. This part cuts both variables into
the categories a clinic would use, crosses them, and reads the hazard ratio in every cell.

CATEGORIES
    apnea-hypopnea index   <5, 5 to <15, 15 to <30, >=30 events/h
    total sleep time       <5 h, 5 to <7 h, >=7 h
    sleep efficiency       <85%, >=85%
    T90                    <=1%, >1 to 5%, >5 to 10%, >10%

Each conventional measure is crossed with T90, giving 16, 12 and 8 joint categories. The
reference is always T90 <=1% together with the most favourable conventional category: an
apnea-hypopnea index below 5, seven hours or more of sleep, sleep efficiency of 85% or more.

MODEL
    Cox proportional hazards, time zero the sleep study, adjusted for age, sex and hospital.
    Age is a natural cubic spline with 4 degrees of freedom, ONE COLUMN DROPPED, and the fit is
    UNPENALIZED. The four basis columns sum to one. A Cox partial likelihood has no intercept
    so the redundancy is harmless under a ridge, but an unpenalized fit on a rank-deficient
    design can return coefficients of zero WITHOUT FAILING. Gate G4 below demonstrates that
    failure on real data rather than asserting it. Hospital enters as a stratum, so it is
    adjusted for without assuming its hazards are proportional.

WHAT IS REPORTED, in the order asked for
    1  n and events in every cell, computed and printed BEFORE any model is fitted, with every
       cell too sparse to fit named.
    2  hazard ratio for every joint category against the reference cell.
    3  inside every conventional stratum, the hazard ratio for T90 >10% against T90 <=1%.
       This is the key number. It is a contrast of two coefficients of the same saturated
       model, so it does not depend on which cell is the reference. Gate G6 proves that, and
       gate G7 shows it barely moves when each stratum is refitted on its own.
    4  interaction between T90 and each conventional measure, two ways: a single ordinal
       product term (1 df, the higher-powered test) and the omnibus Wald test on all
       (K-1)x(J-1) interaction contrasts. EVERY interaction P VALUE IS REPORTED BESIDE THE
       MINIMUM DETECTABLE INTERACTION EFFECT, the smallest ratio-of-hazard-ratios that this
       outcome had 80% power to find at alpha 0.05, computed from the observed standard error
       of the contrast and not from its P value. A large P beside a minimum detectable effect
       of 2.5 means the study could not see a doubling, not that there is nothing there.
       The within-stratum hazard ratios and their intervals are the primary evidence.
    5  sensitivity: T90 recoded as >10% against <=10%. In the four-band analysis the key
       contrast is carried only by the 13,204 patients in the top and bottom bands even though
       all 19,173 are in the model. The binary recode puts every patient into the contrast.
    6  the five negative controls are fitted throughout and flagged, never dropped, and
       Benjamini-Hochberg is applied within each analysis, measure and term across the panel,
       twice: q_all over all 52 outcomes and q_main over the 47 that are not controls, because
       a negative control is a diagnostic and not a hypothesis.

OUTCOME PANEL
    The 52 clinical conditions of numbers/disease_definitions.py. That is every entry except
    the cardiovascular composite, which is the union of five conditions already in the panel
    and would double-count them. Death and the composite are fitted anyway and carried in the
    output with in_panel_52 = False, so nothing is silently cut, and they are excluded from the
    Benjamini-Hochberg families so the multiplicity correction is over the 52 that were asked
    for. Four conditions are diagnosed by something the sleep study itself measures (insomnia,
    restless legs, nocturia, epilepsy) and obesity hypoventilation is defined partly by
    nocturnal desaturation. All five stay in the panel and carry circular = True.

Writes joint_results.csv (every row of every table above) and joint_summary.json in this
directory.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import json
import os
import sys
import time
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from lifelines import CoxPHFitter
from patsy import dmatrix
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT = paths.T90_ROOT
HERE = f"{paths.SV_ROOT}/q8_after_oxygen"
SCRATCH = f"{paths.T90_ROOT}/scratch/q8_after_oxygen"
sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
from cohort_spec import COHORT_N, NUMBERS_DIR, apply_cohort           # noqa: E402  # v8 F1: NUMBERS_DIR for the panel check
from disease_definitions import (DISEASES, NEGATIVE_CONTROLS,        # noqa: E402
                                 CIRCULAR, ORGAN_GROUP)

os.makedirs(SCRATCH, exist_ok=True)
STORE = f"{SCRATCH}/q8_joint_cohort.pkl"

SITES = ("I0002", "I0006")
Z_A, Z_B = 1.959963985, 0.841621234      # alpha 0.05 two-sided, 80% power
MIN_EVENTS_TOTAL = 20                    # an outcome must clear this to be modelled at all
SPARSE_CELL_EVENTS = 10                  # fewer than this and the cell HR is flagged unstable
THIN_EPP = 10                            # events per parameter below this and the fit is thin
N_JOBS = 7

# ------------------------------------------------------------------------------------------
# Categories. `score` is the distance from the most favourable category, so 0 is always the
# reference level and larger is worse. Only the ordinal interaction term uses it; the omnibus
# test assumes no ordering at all.
# ------------------------------------------------------------------------------------------
T90_BANDS = {
    "primary_4band": [("t1", "T90 <=1%", 0), ("t2", "T90 >1-5%", 1),
                      ("t3", "T90 >5-10%", 2), ("t4", "T90 >10%", 3)],
    "sens_binary":   [("s1", "T90 <=10%", 0), ("s2", "T90 >10%", 1)],
}
T90_HI = {"primary_4band": "t4", "sens_binary": "s2"}   # the exposed level of the key contrast
T90_LO = {"primary_4band": "t1", "sens_binary": "s1"}   # the unexposed level

MEASURES = {
    "ahi": {"col": "AHI", "name": "Apnea-hypopnea index",
            "cats": [("a1", "AHI <5", 0), ("a2", "AHI 5-<15", 1),
                     ("a3", "AHI 15-<30", 2), ("a4", "AHI >=30", 3)],
            "ref": "a1", "ref_reason": "no sleep apnea by the conventional cut"},
    "tst": {"col": "TST_min", "name": "Total sleep time",
            "cats": [("d1", "TST <5 h", 2), ("d2", "TST 5-<7 h", 1),
                     ("d3", "TST >=7 h", 0)],
            "ref": "d3", "ref_reason": "meets the recommended duration"},
    "se": {"col": "sleep_efficiency_pct", "name": "Sleep efficiency",
           "cats": [("e1", "SE <85%", 1), ("e2", "SE >=85%", 0)],
           "ref": "e2", "ref_reason": "normal sleep efficiency"},
}
ANALYSES = ["primary_4band", "sens_binary"]

PANEL_52 = [k for k in DISEASES if k != "cvd"]
# v8 F1 fix (2026-09-12, integrity gate 6): the panel is checked against the primary analysis' own disease panel,
# numbers/bdsp_diseases_v3.csv (run_primary_unpenalized.py, step 111, provenance sidecar, one row per entry of DISEASES
# plus death), never against a typed number. cohort_spec.N_RANKED_OUTCOMES (52) is the RANKING's outcome count, not this
# panel's: the panel is every disease except the composite, 56 in v8 (52 in v7). The name PANEL_52 is kept.
_PRIM_PANEL = set(pd.read_csv(f"{NUMBERS_DIR}/bdsp_diseases_v3.csv").key) - {"cvd", "death"}
assert set(PANEL_52) == _PRIM_PANEL, (f"joint panel ({len(PANEL_52)}) is not the primary table's disease panel "
                                       f"({len(_PRIM_PANEL)}): {sorted(set(PANEL_52) ^ _PRIM_PANEL)}")
EXTRAS = ["cvd", "death"]
LABEL = {k: v[0] for k, v in DISEASES.items()}
LABEL["death"] = "Death from any cause"
# obesity hypoventilation is defined partly by nocturnal desaturation, so its association with
# T90 is mechanical. It stays in the 52 and is flagged, exactly as the rest of the manuscript.
CIRCULAR_FLAG = set(CIRCULAR) | {"obesity_hypovent"}


def log(m=""):
    print(m, flush=True)


def hms(s):
    return f"{int(s // 60)}m{int(s % 60):02d}s"


# ==========================================================================================
# COHORT
# ==========================================================================================
def build():
    b = apply_cohort(pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet"))
    assert len(b) == COHORT_N
    assert set(b.site_id.unique()) == set(SITES), b.site_id.unique()
    log(f"cohort {len(b):,} at hospitals {SITES[0]} (n={int((b.site_id == SITES[0]).sum()):,}) "
        f"and {SITES[1]} (n={int((b.site_id == SITES[1]).sum()):,})")

    initial = b.sex.astype(str).str.upper().str[0]
    b["male"] = np.where(initial == "M", 1.0, np.where(initial == "F", 0.0, np.nan))
    log(f"sex: {int((b.male == 1).sum()):,} male, {int((b.male == 0).sum()):,} female, "
        f"{int(b.male.isna().sum())} neither")

    # age spline knots are placed on the WHOLE cohort once, not inside each outcome, so the
    # adjustment is the same function of age in every model and the 52 are comparable
    sp = dmatrix("cr(a, df=4) - 1", {"a": b.AgeAtVisit.values}, return_type="dataframe")
    for i in range(4):
        b[f"age_all{i}"] = sp.iloc[:, i].values
    for i in range(3):
        b[f"age_s{i}"] = sp.iloc[:, i + 1].values          # column 0 dropped, full rank

    # ---- categories ----------------------------------------------------------------------
    t = b.spo2_pct_below_90
    assert t.notna().all(), "apply_cohort should already require T90"
    b["cat_primary_4band"] = np.select(
        [t <= 1, (t > 1) & (t <= 5), (t > 5) & (t <= 10), t > 10],
        ["t1", "t2", "t3", "t4"], default="NA")
    b["cat_sens_binary"] = np.select([t <= 10, t > 10], ["s1", "s2"], default="NA")
    a, d, e = b.AHI, b.TST_min, b.sleep_efficiency_pct
    b["cat_ahi"] = np.select([a < 5, (a >= 5) & (a < 15), (a >= 15) & (a < 30), a >= 30],
                             ["a1", "a2", "a3", "a4"], default="NA")
    b["cat_tst"] = np.select([d < 300, (d >= 300) & (d < 420), d >= 420],
                             ["d1", "d2", "d3"], default="NA")
    b["cat_se"] = np.select([e < 85, e >= 85], ["e1", "e2"], default="NA")

    # G2. every banding must be exhaustive and mutually exclusive: NA only where the source
    # value is missing, and the category counts must sum back to the number of non-missing.
    gates = {}
    for col, src, cats in (("cat_primary_4band", t, [c for c, _, _ in T90_BANDS["primary_4band"]]),
                           ("cat_sens_binary", t, [c for c, _, _ in T90_BANDS["sens_binary"]]),
                           ("cat_ahi", a, [c for c, _, _ in MEASURES["ahi"]["cats"]]),
                           ("cat_tst", d, [c for c, _, _ in MEASURES["tst"]["cats"]]),
                           ("cat_se", e, [c for c, _, _ in MEASURES["se"]["cats"]])):
        assert ((b[col] == "NA") == src.isna()).all(), f"{col}: banding is not exhaustive"
        assert set(b[col].unique()) <= set(cats) | {"NA"}, f"{col}: stray level"
        assert int(b[col].isin(cats).sum()) == int(src.notna().sum()), f"{col}: counts differ"
        gates[col] = {"n_categorised": int(src.notna().sum()),
                      "n_missing": int(src.isna().sum())}
    log("G2 all five bandings are exhaustive and mutually exclusive. "
        f"missing: TST {gates['cat_tst']['n_missing']}, SE {gates['cat_se']['n_missing']}, "
        f"AHI {gates['cat_ahi']['n_missing']}, T90 {gates['cat_primary_4band']['n_missing']}")
    return b, gates


# ==========================================================================================
# 1. CELL CENSUS, BEFORE ANY MODEL IS FITTED
# ==========================================================================================
def eligible(b, key):
    """Incident frame for one condition: not prevalent at the study, positive follow-up."""
    yc, ec, pc = f"{key}_years", f"{key}_incident", f"{key}_prevalent"
    if not all(c in b.columns for c in (yc, ec, pc)):
        return None
    return b[(b[pc] == 0) & b[yc].notna() & (b[yc] > 0)]


def census(b, outcomes):
    """n and events in every cell of every crossing, for the cohort and for every outcome."""
    crows = []
    for analysis in ANALYSES:
        tcol = f"cat_{analysis}"
        for mk, M in MEASURES.items():
            ccol = f"cat_{mk}"
            for ci, clab, _ in M["cats"]:
                for tj, tlab, _ in T90_BANDS[analysis]:
                    sel = (b[ccol] == ci) & (b[tcol] == tj)
                    crows.append({"row_type": "cell_cohort", "analysis": analysis,
                                  "measure": mk, "measure_name": M["name"],
                                  "conv_cat": ci, "conv_label": clab,
                                  "t90_cat": tj, "t90_label": tlab,
                                  "is_reference": ci == M["ref"] and tj == T90_LO[analysis],
                                  "n": int(sel.sum()),
                                  "pct_of_cohort": round(100 * float(sel.mean()), 2)})

    rows = []
    for key in outcomes:
        f = eligible(b, key)
        if f is None:
            continue
        yc, ec = f"{key}_years", f"{key}_incident"
        meta = {"outcome": key, "label": LABEL[key],
                "organ_group": ORGAN_GROUP.get(key, "Death"),
                "negative_control": key in NEGATIVE_CONTROLS,
                "circular": key in CIRCULAR_FLAG, "in_panel_52": key in PANEL_52}
        for analysis in ANALYSES:
            tcol = f"cat_{analysis}"
            for mk, M in MEASURES.items():
                ccol = f"cat_{mk}"
                grp = f.groupby([ccol, tcol], observed=True).agg(
                    n=(ec, "size"), events=(ec, "sum"), py=(yc, "sum"))
                for ci, clab, _ in M["cats"]:
                    for tj, tlab, _ in T90_BANDS[analysis]:
                        if (ci, tj) in grp.index:
                            g = grp.loc[(ci, tj)]
                            n, ev, py = int(g.n), int(g.events), float(g.py)
                        else:
                            n, ev, py = 0, 0, 0.0
                        rows.append({
                            "row_type": "cell", "analysis": analysis, "measure": mk,
                            "measure_name": M["name"], **meta,
                            "conv_cat": ci, "conv_label": clab, "t90_cat": tj,
                            "t90_label": tlab,
                            "is_reference": ci == M["ref"] and tj == T90_LO[analysis],
                            "n": n, "events": ev, "person_years": round(py, 1),
                            "rate_per_1000py": round(1000 * ev / py, 2) if py > 0 else np.nan,
                            "zero_events": ev == 0, "sparse_cell": ev < SPARSE_CELL_EVENTS})
    return pd.DataFrame(rows), pd.DataFrame(crows)


# ==========================================================================================
# 2-5. MODEL MACHINERY
# ==========================================================================================
AGE = ["age_s0", "age_s1", "age_s2"]
AGE_ALL = ["age_all0", "age_all1", "age_all2", "age_all3"]


def cellterm(ci, tj):
    return f"X_{ci}_{tj}"


def fit_cox(d, terms, age_cols=AGE):
    cov = list(terms) + ["male"] + list(age_cols)
    return CoxPHFitter(penalizer=0.0).fit(d[cov + ["T", "E", "site"]], "T", "E",
                                          strata=["site"])


def est(m, a, ref=None):
    """HR, 95% CI, P and SE for a coefficient, or for the difference of two coefficients."""
    beta, V = m.params_, m.variance_matrix_
    if ref is None:
        delta, var = float(beta[a]), float(V.loc[a, a])
    else:
        delta = float(beta[a] - beta[ref])
        var = float(V.loc[a, a] + V.loc[ref, ref] - 2 * V.loc[a, ref])
    se = float(np.sqrt(max(var, 0.0)))
    if not np.isfinite(se) or se <= 0:
        return dict(hr=np.nan, lo=np.nan, hi=np.nan, p=np.nan, se=np.nan)
    return dict(hr=float(np.exp(delta)), lo=float(np.exp(delta - Z_A * se)),
                hi=float(np.exp(delta + Z_A * se)),
                p=float(2 * stats.norm.sf(abs(delta / se))), se=se)


def mde(se):
    """Smallest hazard ratio (or ratio of hazard ratios) detectable at 80% power, alpha 0.05.

    exp((z_alpha + z_beta) * SE). Computed from the observed standard error of the contrast,
    which is a statement about the precision this design achieved. It is NOT observed power:
    nothing here is a transformation of the P value."""
    return float(np.exp((Z_A + Z_B) * se)) if np.isfinite(se) and se > 0 else np.nan


def wald(m, L, names):
    """Joint Wald chi-square for L @ beta = 0 over the coefficients `names`."""
    if L.shape[0] == 0:
        return np.nan, 0, np.nan
    b = m.params_[names].values
    V = m.variance_matrix_.loc[names, names].values
    Lb = L @ b
    LVL = L @ V @ L.T
    try:
        chi2 = float(Lb @ np.linalg.solve(LVL, Lb))
    except np.linalg.LinAlgError:
        return np.nan, int(L.shape[0]), np.nan
    if not np.isfinite(chi2) or chi2 < 0:
        return np.nan, int(L.shape[0]), np.nan
    return chi2, int(L.shape[0]), float(stats.chi2.sf(chi2, L.shape[0]))


def prepare(b, key, mk, analysis):
    """The modelling frame for one outcome, one crossing. None if it cannot be built."""
    f = eligible(b, key)
    if f is None:
        return None
    yc, ec = f"{key}_years", f"{key}_incident"
    ccol, tcol = f"cat_{mk}", f"cat_{analysis}"
    d = pd.DataFrame({"T": f[yc].values, "E": f[ec].astype(int).values,
                      "site": f.site_id.values, "male": f.male.values,
                      "conv": f[ccol].values, "t90": f[tcol].values})
    for c in AGE + AGE_ALL:
        d[c] = f[c].values
    d = d[(d.conv != "NA") & (d.t90 != "NA")].dropna()
    return d


def run_job(key, mk, analysis, store):
    """Everything asked for, for one outcome crossed one way, from ONE saturated Cox fit."""
    b = pd.read_pickle(store)
    M = MEASURES[mk]
    tcats = T90_BANDS[analysis]
    ref_c, ref_t = M["ref"], T90_LO[analysis]
    hi_t = T90_HI[analysis]
    base = {"analysis": analysis, "measure": mk, "measure_name": M["name"], "outcome": key,
            "label": LABEL[key], "organ_group": ORGAN_GROUP.get(key, "Death"),
            "negative_control": key in NEGATIVE_CONTROLS, "circular": key in CIRCULAR_FLAG,
            "in_panel_52": key in PANEL_52,
            "reference_cell": f"{ref_c} x {ref_t}"}
    d = prepare(b, key, mk, analysis)
    if d is None or len(d) == 0:
        return [{**base, "row_type": "skipped", "note": "no usable rows"}]
    if int(d.E.sum()) < MIN_EVENTS_TOTAL:
        return [{**base, "row_type": "skipped",
                 "note": f"{int(d.E.sum())} events, below {MIN_EVENTS_TOTAL}"}]

    ev = {(ci, tj): int(d[(d.conv == ci) & (d.t90 == tj)].E.sum())
          for ci, _, _ in M["cats"] for tj, _, _ in tcats}
    nn = {(ci, tj): int(((d.conv == ci) & (d.t90 == tj)).sum())
          for ci, _, _ in M["cats"] for tj, _, _ in tcats}
    if ev[(ref_c, ref_t)] == 0:
        return [{**base, "row_type": "skipped", "note": "reference cell has no events"}]

    # a cell with no events has a maximum-likelihood coefficient of minus infinity, so its
    # members carry zero hazard and leave every risk set. Removing their rows reproduces that
    # fit exactly and keeps the rest of the crossing estimable. Every such cell is recorded.
    order = [(ci, tj) for ci, _, _ in M["cats"] for tj, _, _ in tcats]   # deterministic
    dropped = [c for c in order if ev[c] == 0]
    keepcells = [c for c in order if ev[c] > 0]
    sub = d[pd.Series(list(zip(d.conv, d.t90)), index=d.index).isin(set(keepcells))].copy()
    terms = []
    for ci, tj in keepcells:
        if (ci, tj) == (ref_c, ref_t):
            continue
        nm = cellterm(ci, tj)
        sub[nm] = ((sub.conv == ci) & (sub.t90 == tj)).astype(float)
        terms.append(nm)

    n_par = len(terms) + 1 + len(AGE)
    epp = float(sub.E.sum()) / max(n_par, 1)
    try:
        m = fit_cox(sub, terms)
    except Exception as exc:
        return [{**base, "row_type": "skipped",
                 "note": f"saturated fit failed: {type(exc).__name__}: {exc}"}]

    rows = []
    common = {**base, "n_in_fit": int(len(sub)), "events_in_fit": int(sub.E.sum()),
              "n_cells_dropped_zero_events": len(dropped),
              "cells_dropped_zero_events": ";".join(f"{c}x{t}" for c, t in dropped) or "",
              "events_per_parameter": round(epp, 2), "thin_fit": epp < THIN_EPP}

    # ---- 2. hazard ratio for every joint category ----------------------------------------
    for ci, clab, _ in M["cats"]:
        for tj, tlab, _ in tcats:
            r = {**common, "row_type": "joint_hr", "conv_cat": ci, "conv_label": clab,
                 "t90_cat": tj, "t90_label": tlab, "term": f"{ci}x{tj}",
                 "n": nn[(ci, tj)], "events": ev[(ci, tj)],
                 "sparse_cell": ev[(ci, tj)] < SPARSE_CELL_EVENTS,
                 "is_reference": (ci, tj) == (ref_c, ref_t)}
            if (ci, tj) == (ref_c, ref_t):
                r.update(hr=1.0, lo=np.nan, hi=np.nan, p=np.nan, se=np.nan,
                         note="reference")
            elif (ci, tj) in dropped:
                r.update(hr=np.nan, note="no events, dropped from the fit")
            else:
                r.update(est(m, cellterm(ci, tj)))
            rows.append(r)

    # ---- 3. inside every conventional stratum, T90 high against T90 low ------------------
    for ci, clab, _ in M["cats"]:
        r = {**common, "row_type": "within_stratum", "conv_cat": ci, "conv_label": clab,
             "term": f"t90_{hi_t}_vs_{ref_t}", "t90_cat": f"{hi_t}|{ref_t}",
             "t90_label": f"{dict((c, l) for c, l, _ in tcats)[hi_t]} vs "
                          f"{dict((c, l) for c, l, _ in tcats)[ref_t]}",
             "n_lo": nn[(ci, ref_t)], "n_hi": nn[(ci, hi_t)],
             "events_lo": ev[(ci, ref_t)], "events_hi": ev[(ci, hi_t)],
             "n": nn[(ci, ref_t)] + nn[(ci, hi_t)],
             "events": ev[(ci, ref_t)] + ev[(ci, hi_t)],
             "sparse_cell": min(ev[(ci, ref_t)], ev[(ci, hi_t)]) < SPARSE_CELL_EVENTS}
        if (ci, hi_t) in dropped or (ci, ref_t) in dropped:
            r.update(hr=np.nan, note="a cell of this contrast has no events")
        else:
            hi_term = cellterm(ci, hi_t)
            lo_term = None if (ci, ref_t) == (ref_c, ref_t) else cellterm(ci, ref_t)
            r.update(est(m, hi_term, lo_term))
        r["mde_hr"] = mde(r.get("se", np.nan))
        rows.append(r)

    # ---- 4. interaction ------------------------------------------------------------------
    # gamma_ij = beta_ij - beta_i,ref_t - beta_ref_c,j is the departure of cell (i,j) from what
    # the two main effects predict. Every gamma is a contrast of the SAME saturated fit, so the
    # omnibus test below is the standard interaction test without a second model.
    names = terms
    idx = {nm: k for k, nm in enumerate(names)}

    def vec(ci, tj):
        if (ci, tj) == (ref_c, ref_t):
            return np.zeros(len(names))
        v = np.zeros(len(names))
        v[idx[cellterm(ci, tj)]] = 1.0
        return v

    Lrows, Lnames = [], []
    for ci, _, _ in M["cats"]:
        if ci == ref_c:
            continue
        for tj, _, _ in tcats:
            if tj == ref_t:
                continue
            if any(c in dropped for c in ((ci, tj), (ci, ref_t), (ref_c, tj))):
                continue
            Lrows.append(vec(ci, tj) - vec(ci, ref_t) - vec(ref_c, tj))
            Lnames.append(f"{ci}x{tj}")
    L = np.array(Lrows) if Lrows else np.zeros((0, len(names)))
    chi2, dfree, pj = wald(m, L, names)
    n_possible = (len(M["cats"]) - 1) * (len(tcats) - 1)

    # the corner: how much the T90 high-vs-low hazard ratio differs between the worst and the
    # most favourable conventional stratum. Directly comparable with the within-stratum rows.
    worst = M["cats"][int(np.argmax([s for _, _, s in M["cats"]]))][0]
    corner = {"hr": np.nan, "lo": np.nan, "hi": np.nan, "p": np.nan, "se": np.nan}
    if not any(c in dropped for c in ((worst, hi_t), (worst, ref_t), (ref_c, hi_t))):
        g = vec(worst, hi_t) - vec(worst, ref_t) - vec(ref_c, hi_t)
        bb, VV = m.params_[names].values, m.variance_matrix_.loc[names, names].values
        delta = float(g @ bb)
        se = float(np.sqrt(max(g @ VV @ g, 0.0)))
        if se > 0:
            corner = dict(hr=float(np.exp(delta)), lo=float(np.exp(delta - Z_A * se)),
                          hi=float(np.exp(delta + Z_A * se)),
                          p=float(2 * stats.norm.sf(abs(delta / se))), se=se)

    rows.append({**common, "row_type": "interaction", "term": "omnibus_wald",
                 "chi2": chi2, "df": dfree, "df_possible": n_possible, "p": pj,
                 "n_interaction_contrasts_dropped": n_possible - dfree,
                 "contrasts": ";".join(Lnames),
                 "mde_hr": mde(corner["se"]),
                 "note": "omnibus Wald on every interaction contrast of the saturated model. "
                         "mde_hr is the corner contrast below, the interpretable scale on "
                         "which to judge what this P could and could not have detected"})
    rows.append({**common, "row_type": "interaction", "term": "corner_ratio_of_HR",
                 "conv_cat": worst, "conv_label": dict((c, l) for c, l, _ in M["cats"])[worst],
                 **corner, "mde_hr": mde(corner["se"]),
                 "note": "T90 high-vs-low HR in the worst conventional stratum divided by the "
                         "same HR in the most favourable stratum"})

    # ordinal product term: one degree of freedom, the higher-powered interaction test
    cs = dict((c, s) for c, _, s in M["cats"])
    ts = dict((c, s) for c, _, s in tcats)
    dt = d.copy()
    dt["conv_score"] = dt.conv.map(cs).astype(float)
    dt["t90_score"] = dt.t90.map(ts).astype(float)
    dt["conv_x_t90"] = dt.conv_score * dt.t90_score
    tr = {"hr": np.nan, "lo": np.nan, "hi": np.nan, "p": np.nan, "se": np.nan}
    trend_note = ""
    try:
        mt = fit_cox(dt, ["conv_score", "t90_score", "conv_x_t90"])
        tr = est(mt, "conv_x_t90")
        trend_main = {"trend_conv_hr": float(np.exp(mt.params_["conv_score"])),
                      "trend_t90_hr": float(np.exp(mt.params_["t90_score"]))}
    except Exception as exc:
        trend_main = {}
        trend_note = f"trend fit failed: {type(exc).__name__}"
    rows.append({**common, "row_type": "interaction", "term": "ordinal_product",
                 **tr, "mde_hr": mde(tr["se"]), **trend_main,
                 "n": int(len(dt)), "events": int(dt.E.sum()), "df": 1,
                 "note": trend_note or "change in the T90 hazard ratio per one category step "
                                       "away from the most favourable conventional category"})
    return rows


# ==========================================================================================
# GATES
# ==========================================================================================
def gate_rank_deficiency(b):
    """G4. Show, on real data, what the standing rule protects against."""
    d = prepare(b, "death", "ahi", "primary_4band")
    ev = {(ci, tj): int(d[(d.conv == ci) & (d.t90 == tj)].E.sum())
          for ci, _, _ in MEASURES["ahi"]["cats"] for tj, _, _ in T90_BANDS["primary_4band"]}
    terms = []
    for (ci, tj), v in ev.items():
        if v == 0 or (ci, tj) == ("a1", "t1"):
            continue
        nm = cellterm(ci, tj)
        d[nm] = ((d.conv == ci) & (d.t90 == tj)).astype(float)
        terms.append(nm)
    out = {}
    for tag, cols in (("dropped_one_column", AGE), ("all_four_columns", AGE_ALL)):
        try:
            m = fit_cox(d, terms, age_cols=cols)
            out[tag] = {"fitted": True,
                        "max_abs_age_coef": round(float(np.abs(m.params_[cols]).max()), 6),
                        "log_likelihood": round(float(m.log_likelihood_), 4),
                        "male_hr": round(float(np.exp(m.params_["male"])), 4)}
        except Exception as exc:
            out[tag] = {"fitted": False, "error": f"{type(exc).__name__}: {exc}"}
    ok = out["dropped_one_column"]["fitted"] and out["dropped_one_column"]["max_abs_age_coef"] > 0.1
    assert ok, f"G4: the age spline collapsed in the primary specification: {out}"
    return out


def gate_saturated_equivalence(b, key="death", mk="ahi", analysis="primary_4band"):
    """G5. The cell model and the main-effects-plus-products model are the same fit, and the
    omnibus contrast test equals the Wald test on the explicit product terms."""
    M, tcats = MEASURES[mk], T90_BANDS[analysis]
    ref_c, ref_t = M["ref"], T90_LO[analysis]
    d = prepare(b, key, mk, analysis)
    ev = {(ci, tj): int(d[(d.conv == ci) & (d.t90 == tj)].E.sum())
          for ci, _, _ in M["cats"] for tj, _, _ in tcats}
    assert min(ev.values()) > 0, "G5 needs a crossing with no empty cell"

    cell_terms = []
    dc = d.copy()
    for ci, _, _ in M["cats"]:
        for tj, _, _ in tcats:
            if (ci, tj) == (ref_c, ref_t):
                continue
            nm = cellterm(ci, tj)
            dc[nm] = ((dc.conv == ci) & (dc.t90 == tj)).astype(float)
            cell_terms.append(nm)
    mc = fit_cox(dc, cell_terms)

    dp = d.copy()
    main_c = [f"C_{ci}" for ci, _, _ in M["cats"] if ci != ref_c]
    main_t = [f"T_{tj}" for tj, _, _ in tcats if tj != ref_t]
    prod = []
    for ci, _, _ in M["cats"]:
        if ci == ref_c:
            continue
        dp[f"C_{ci}"] = (dp.conv == ci).astype(float)
    for tj, _, _ in tcats:
        if tj == ref_t:
            continue
        dp[f"T_{tj}"] = (dp.t90 == tj).astype(float)
    for ci, _, _ in M["cats"]:
        if ci == ref_c:
            continue
        for tj, _, _ in tcats:
            if tj == ref_t:
                continue
            dp[f"P_{ci}_{tj}"] = dp[f"C_{ci}"] * dp[f"T_{tj}"]
            prod.append(f"P_{ci}_{tj}")
    mp = fit_cox(dp, main_c + main_t + prod)

    dll = abs(float(mc.log_likelihood_) - float(mp.log_likelihood_))
    assert dll < 1e-4, f"G5: the two parameterisations differ, log-lik gap {dll}"

    idx = {nm: k for k, nm in enumerate(cell_terms)}

    def vec(ci, tj):
        v = np.zeros(len(cell_terms))
        if (ci, tj) != (ref_c, ref_t):
            v[idx[cellterm(ci, tj)]] = 1.0
        return v

    L = np.array([vec(ci, tj) - vec(ci, ref_t) - vec(ref_c, tj)
                  for ci, _, _ in M["cats"] if ci != ref_c
                  for tj, _, _ in tcats if tj != ref_t])
    chi2a, dfa, pa = wald(mc, L, cell_terms)
    Lp = np.zeros((len(prod), len(main_c + main_t + prod)))
    allp = main_c + main_t + prod
    for r, nm in enumerate(prod):
        Lp[r, allp.index(nm)] = 1.0
    chi2b, dfb, pb = wald(mp, Lp, allp)
    assert dfa == dfb and abs(chi2a - chi2b) < 1e-5, f"G5: chi2 {chi2a} vs {chi2b}"
    return {"outcome": key, "measure": mk, "log_lik_cells": float(mc.log_likelihood_),
            "log_lik_products": float(mp.log_likelihood_), "log_lik_gap": dll,
            "omnibus_chi2_from_contrasts": round(chi2a, 6),
            "omnibus_chi2_from_products": round(chi2b, 6), "df": int(dfa),
            "p": float(pa)}


def gate_stratum_only_refit(b, keys=("death", "hf", "diabetes", "htn2")):
    """G7. What the shared age and sex effects cost.

    The within-stratum hazard ratio is a contrast of one saturated model, so age and sex are
    fitted once across the whole crossing. A reader may reasonably ask whether that borrowing
    is doing the work. This refits the same contrast inside each stratum on its own, where age
    and sex are free to differ, and reports how far the two disagree. They are different
    estimators and are not expected to be identical, so the assertion is loose enough to catch
    a coding error and not a modelling nuance."""
    rows = []
    for key in keys:
        for mk, M in MEASURES.items():
            d = prepare(b, key, mk, "primary_4band")
            for ci, clab, _ in M["cats"]:
                s = d[(d.conv == ci) & d.t90.isin(["t1", "t4"])].copy()
                s["hi"] = (s.t90 == "t4").astype(float)
                try:
                    m = fit_cox(s, ["hi"])
                    hr_only = float(np.exp(m.params_["hi"]))
                except Exception:
                    hr_only = np.nan
                rows.append({"outcome": key, "measure": mk, "conv_cat": ci,
                             "conv_label": clab, "hr_stratum_only": hr_only})
    return pd.DataFrame(rows)


def gate_reference_invariance(b, key="death", mk="ahi", analysis="primary_4band"):
    """G6. The within-stratum T90 hazard ratios do not depend on which cell is the reference."""
    M, tcats = MEASURES[mk], T90_BANDS[analysis]
    hi_t, lo_t = T90_HI[analysis], T90_LO[analysis]
    d = prepare(b, key, mk, analysis)
    out = {}
    for altref in [(M["ref"], lo_t), ("a4", "t3")]:
        dc = d.copy()
        cell_terms = []
        for ci, _, _ in M["cats"]:
            for tj, _, _ in tcats:
                if (ci, tj) == altref:
                    continue
                nm = cellterm(ci, tj)
                dc[nm] = ((dc.conv == ci) & (dc.t90 == tj)).astype(float)
                cell_terms.append(nm)
        m = fit_cox(dc, cell_terms)
        hrs = {}
        for ci, _, _ in M["cats"]:
            a = None if (ci, hi_t) == altref else cellterm(ci, hi_t)
            r = None if (ci, lo_t) == altref else cellterm(ci, lo_t)
            if a is None:                       # the high cell is the reference
                hrs[ci] = float(np.exp(-m.params_[r]))
            else:
                hrs[ci] = float(np.exp(m.params_[a] - (0.0 if r is None else m.params_[r])))
        out[f"{altref[0]}x{altref[1]}"] = hrs
    a, bb = list(out.values())
    gap = max(abs(a[k] - bb[k]) for k in a)
    assert gap < 1e-6, f"G6: within-stratum HRs moved by {gap} when the reference changed"
    return {"outcome": key, "measure": mk, "max_abs_difference": gap, "hrs": out}


# ==========================================================================================
# MAIN
# ==========================================================================================
def main():
    t0 = time.time()
    meta = {"run_started": time.strftime("%Y-%m-%d %H:%M:%S"), "inputs": {}}
    for p in (f"{paths.TABLES_DIR}/t90_final.parquet",
              f"{paths.NUMBERS_DIR}/disease_definitions.py",
              f"{paths.NUMBERS_DIR}/cohort_spec.py"):
        s = os.stat(p)
        meta["inputs"][p.rsplit("/", 1)[-1]] = {
            "bytes": s.st_size,
            "modified": time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(s.st_mtime))}
        log(f"input {p.rsplit('/', 1)[-1]:<26}{s.st_size:>12,} bytes  modified "
            f"{meta['inputs'][p.rsplit('/', 1)[-1]]['modified']}")

    b, band_gates = build()
    meta["banding"] = band_gates
    log(f"negative controls kept throughout: {NEGATIVE_CONTROLS}")
    log(f"panel {len(PANEL_52)} clinical outcomes, plus {EXTRAS} carried but outside the "
        f"Benjamini-Hochberg families")

    outcomes = PANEL_52 + EXTRAS
    have = [k for k in outcomes if f"{k}_incident" in b.columns]
    assert set(have) == set(outcomes), f"missing from the parquet: {set(outcomes) - set(have)}"
    b = b.reset_index(drop=True)
    # the workers reload this frame once per job, so the diagnosis-date columns are dropped
    keepcols = [c for c in b.columns if not c.endswith("_first_date")]
    b = b[keepcols]
    b.to_pickle(STORE)
    log(f"worker frame {len(b):,} x {len(keepcols)} columns -> {STORE}")

    # ---------------------------------------------------------------------------------- 1
    log("\n" + "=" * 100)
    log("1. CELL CENSUS, BEFORE ANY MODEL IS FITTED")
    log("=" * 100)
    cells, cohort_cells = census(b, outcomes)
    for analysis in ANALYSES:
        for mk, M in MEASURES.items():
            cc = cohort_cells[(cohort_cells.analysis == analysis) &
                              (cohort_cells.measure == mk)]
            tcats = T90_BANDS[analysis]
            log(f"\n{M['name']} x T90, {analysis}: participants per cell")
            head = f"{'':<16}" + "".join(f"{l:>16}" for _, l, _ in tcats) + f"{'total':>10}"
            log(head)
            for ci, clab, _ in M["cats"]:
                vals = [int(cc[(cc.conv_cat == ci) & (cc.t90_cat == tj)].n.iloc[0])
                        for tj, _, _ in tcats]
                mark = "  <- reference row" if ci == M["ref"] else ""
                log(f"{clab:<16}" + "".join(f"{v:>16,}" for v in vals) +
                    f"{sum(vals):>10,}{mark}")
            tot = [int(cc[cc.t90_cat == tj].n.sum()) for tj, _, _ in tcats]
            log(f"{'total':<16}" + "".join(f"{v:>16,}" for v in tot) + f"{sum(tot):>10,}")
    meta["cohort_cells"] = json.loads(cohort_cells.to_json(orient="records"))

    smallest = (cohort_cells[cohort_cells.analysis == "primary_4band"]
                .sort_values("n").head(5))
    log(f"\nsmallest cohort cells: " + ", ".join(
        f"{r.measure_name} {r.conv_label} & {r.t90_label} n={int(r.n):,}"
        for _, r in smallest.iterrows()))

    # every cell too sparse to fit, named
    c52 = cells[(cells.row_type == "cell") & cells.in_panel_52]
    zero = c52[c52.zero_events]
    sparse = c52[(~c52.zero_events) & c52.sparse_cell]
    log(f"\nof {len(c52):,} outcome-by-cell combinations across the {len(PANEL_52)} outcomes, "
        f"{len(zero):,} have NO events and are dropped from their fit, and {len(sparse):,} "
        f"more have 1 to {SPARSE_CELL_EVENTS - 1} events and are flagged unstable but kept.")
    if len(zero):
        log("\nevery cell too sparse to fit (no events), by outcome:")
        for (an, mk), g in zero.groupby(["analysis", "measure"]):
            for out, gg in g.groupby("outcome"):
                log(f"  {an:<14}{MEASURES[mk]['name']:<22}{LABEL[out]:<26}"
                    + ", ".join(f"{r.conv_label} & {r.t90_label}" for _, r in gg.iterrows()))
    zc = (zero.groupby(["measure", "conv_label", "t90_label"]).size()
          .sort_values(ascending=False).head(8))
    log("\nthe cells that empty out most often: " +
        ", ".join(f"{a[1]} & {a[2]} ({v} outcomes)" for a, v in zc.items()))
    meta["sparse"] = {"n_cell_rows_52": int(len(c52)), "n_zero_event": int(len(zero)),
                      "n_sparse_nonzero": int(len(sparse)),
                      "zero_event_cells": [f"{r.analysis}|{r.measure}|{r.outcome}|"
                                           f"{r.conv_cat}x{r.t90_cat}"
                                           for _, r in zero.iterrows()]}

    # -------------------------------------------------------------------------------- gates
    log("\n" + "=" * 100)
    log("GATES")
    log("=" * 100)
    g4 = gate_rank_deficiency(b)
    log(f"G4 rank deficiency, death x apnea-hypopnea index:")
    for tag, v in g4.items():
        log(f"   {tag:<22}{json.dumps(v)}")
    log("   the primary specification drops one age-spline column and its coefficients are "
        "real, which is\n   what the standing rule requires. The four-column fit is shown "
        "beside it so the failure mode is\n   documented on this data rather than asserted.")
    g5 = gate_saturated_equivalence(b)
    log(f"G5 the cell model and the main-effects-plus-products model are one fit: log-lik gap "
        f"{g5['log_lik_gap']:.2e},\n   omnibus chi2 {g5['omnibus_chi2_from_contrasts']} from "
        f"contrasts against {g5['omnibus_chi2_from_products']} from explicit product terms, "
        f"{g5['df']} df")
    g6 = gate_reference_invariance(b)
    log(f"G6 within-stratum T90 hazard ratios move by at most {g6['max_abs_difference']:.2e} "
        f"when the reference cell\n   changes, so the key number of this analysis does not "
        f"depend on the reference at all")
    g7 = gate_stratum_only_refit(b)
    meta["gates"] = {"G4_rank_deficiency": g4, "G5_saturated_equivalence": g5,
                     "G6_reference_invariance": g6}

    # ------------------------------------------------------------------------------ 2 to 5
    jobs = [(k, mk, an) for an in ANALYSES for mk in MEASURES for k in outcomes]
    log(f"\nfitting {len(jobs):,} outcome-by-crossing jobs "
        f"({len(outcomes)} outcomes x {len(MEASURES)} measures x {len(ANALYSES)} analyses), "
        f"two Cox fits each")
    res = Parallel(n_jobs=N_JOBS, verbose=5, backend="loky")(
        delayed(run_job)(k, mk, an, STORE) for k, mk, an in jobs)
    model_rows = [r for rr in res for r in rr]
    R = pd.DataFrame(model_rows)
    skipped = R[R.row_type == "skipped"]
    log(f"\n{len(skipped)} of {len(jobs)} jobs skipped" +
        ("" if not len(skipped) else ":"))
    for _, r in skipped.iterrows():
        log(f"    {r.analysis} {r.measure} {r.label}: {r['note']}")
    meta["skipped"] = json.loads(skipped.to_json(orient="records")) if len(skipped) else []

    # ---------------------------------------------------------------------------------- 6
    # Benjamini-Hochberg across the 52 outcomes, within each analysis, measure, term AND
    # conventional stratum, so that one family is one test asked of every outcome once. The
    # stratum has to be in the key: without it the four within-stratum tests of the apnea
    # index would be pooled into a family of 208 and the correction would depend on how many
    # strata the measure happens to have. q_all uses all 52. q_main drops the five negative
    # controls, because a control is a diagnostic and not a hypothesis. Death and the
    # composite are outside the panel and carry no q.
    R["q_all"] = np.nan
    R["q_main"] = np.nan
    R["_fam"] = R.conv_cat.fillna("") if "conv_cat" in R.columns else ""
    fam = R[R.row_type.isin(["joint_hr", "within_stratum", "interaction"]) &
            R.in_panel_52.fillna(False) & R.p.notna()]
    nfam, ntest = 0, 0
    for _, g in fam.groupby(["analysis", "measure", "term", "_fam"]):
        if len(g) > 1:
            R.loc[g.index, "q_all"] = multipletests(g.p.values, method="fdr_bh")[1]
            nfam += 1
            ntest += len(g)
        gm = g[~g.negative_control.astype(bool)]
        if len(gm) > 1:
            R.loc[gm.index, "q_main"] = multipletests(gm.p.values, method="fdr_bh")[1]
    R = R.drop(columns=["_fam"])
    log(f"Benjamini-Hochberg applied within {nfam} families of (analysis, measure, term, "
        f"stratum), {ntest:,} tests in all, each family running across the {len(PANEL_52)} outcomes")

    # G7 needs the modelled hazard ratios, so it is closed here rather than beside the others
    wref = R[(R.row_type == "within_stratum") & (R.analysis == "primary_4band")][
        ["outcome", "measure", "conv_cat", "hr"]]
    g7 = g7.merge(wref, on=["outcome", "measure", "conv_cat"], how="left")
    g7["pct_diff"] = 100 * (g7.hr_stratum_only / g7.hr - 1)
    worst = float(g7.pct_diff.abs().max())
    assert worst < 40, f"G7: the joint contrast and the stratum-only refit differ by {worst}%"
    log(f"G7 on {len(g7)} probe contrasts the joint-model within-stratum hazard ratio and a "
        f"Cox refitted inside\n   that stratum alone agree to within {worst:.1f}% "
        f"(median {float(g7.pct_diff.abs().median()):.1f}%), so sharing the age and sex "
        f"effects across\n   strata is not what produces the key number")
    meta["gates"]["G7_stratum_only_refit"] = {
        "n_probes": int(len(g7)), "max_abs_pct_diff": round(worst, 2),
        "median_abs_pct_diff": round(float(g7.pct_diff.abs().median()), 2),
        "probes": json.loads(g7.to_json(orient="records"))}

    # ------------------------------------------------------------------------------ summary
    sumrows = []
    for an in ANALYSES:
        for mk, M in MEASURES.items():
            for ci, clab, _ in M["cats"]:
                w = R[(R.row_type == "within_stratum") & (R.analysis == an) &
                      (R.measure == mk) & (R.conv_cat == ci) & R.in_panel_52.fillna(False) &
                      R.hr.notna()]
                real = w[~w.negative_control.astype(bool)]
                neg = w[w.negative_control.astype(bool)]
                sumrows.append({
                    "row_type": "summary_within_stratum", "analysis": an, "measure": mk,
                    "measure_name": M["name"], "conv_cat": ci, "conv_label": clab,
                    "term": f"t90_{T90_HI[an]}_vs_{T90_LO[an]}",
                    "n_outcomes": int(len(w)), "n_outcomes_noncontrol": int(len(real)),
                    "hr": float(real.hr.median()) if len(real) else np.nan,
                    "lo": float(real.hr.quantile(0.25)) if len(real) else np.nan,
                    "hi": float(real.hr.quantile(0.75)) if len(real) else np.nan,
                    "n_hr_above_1": int((real.hr > 1).sum()),
                    "n_ci_excludes_1": int(((real.lo > 1) | (real.hi < 1)).sum()),
                    "n_ci_above_1": int((real.lo > 1).sum()),
                    "n_q_main_lt_005": int((real.q_main < 0.05).sum()),
                    "control_hr_median": float(neg.hr.median()) if len(neg) else np.nan,
                    "control_hr_max": float(neg.hr.max()) if len(neg) else np.nan,
                    "control_n_ci_above_1": int((neg.lo > 1).sum()) if len(neg) else 0,
                    "mde_hr": float(real.mde_hr.median()) if len(real) else np.nan,
                    "note": "median across the non-control outcomes, interquartile range in "
                            "lo and hi"})
            for tm in ("ordinal_product", "omnibus_wald", "corner_ratio_of_HR"):
                it = R[(R.row_type == "interaction") & (R.analysis == an) &
                       (R.measure == mk) & (R.term == tm) &
                       R.in_panel_52.fillna(False) & R.p.notna()]
                real = it[~it.negative_control.astype(bool)]
                sumrows.append({
                    "row_type": "summary_interaction", "analysis": an, "measure": mk,
                    "measure_name": M["name"], "term": tm, "n_outcomes": int(len(it)),
                    "n_outcomes_noncontrol": int(len(real)),
                    "n_p_lt_005": int((it.p < 0.05).sum()),
                    "n_q_main_lt_005": int((it.q_main < 0.05).sum()),
                    "p": float(it.p.median()), "min_p": float(it.p.min()),
                    "mde_hr": float(it.mde_hr.median()) if it.mde_hr.notna().any() else np.nan,
                    "note": "median P and median minimum detectable interaction effect across "
                            "the panel. A nonsignificant P is not evidence of no interaction"})
    S = pd.DataFrame(sumrows)

    out = pd.concat([cohort_cells, cells, R, S], ignore_index=True, sort=False)
    front = ["row_type", "analysis", "measure", "measure_name", "outcome", "label",
             "organ_group", "negative_control", "circular", "in_panel_52", "conv_cat",
             "conv_label", "t90_cat", "t90_label", "term", "is_reference", "reference_cell",
             "n", "events", "person_years", "rate_per_1000py", "n_lo", "n_hi", "events_lo",
             "events_hi", "hr", "lo", "hi", "p", "se", "q_all", "q_main", "mde_hr",
             "chi2", "df", "df_possible", "sparse_cell", "zero_events", "thin_fit",
             "events_per_parameter", "n_in_fit", "events_in_fit",
             "n_cells_dropped_zero_events", "cells_dropped_zero_events", "note"]
    cols = [c for c in front if c in out.columns] + [c for c in out.columns if c not in front]
    out = out[cols]
    out.to_csv(f"{HERE}/joint_results.csv", index=False)
    log(f"\nwrote joint_results.csv, {len(out):,} rows x {len(cols)} columns")

    # ==========================================================================================
    # REPORT
    # ==========================================================================================
    def wtab(an, mk, keys, title):
        M = MEASURES[mk]
        log(f"\n{title}")
        log(f"{'condition':<26}{'events':>7}" +
            "".join(f"{l:>22}" for _, l, _ in M["cats"]) +
            f"{'P trend':>10}{'P omnibus':>11}{'min det/step':>14}{'min det corner':>16}")
        for k in keys:
            w = R[(R.row_type == "within_stratum") & (R.analysis == an) & (R.measure == mk) &
                  (R.outcome == k)]
            if not len(w):
                continue

            def pick(term, col):
                s = R[(R.row_type == "interaction") & (R.analysis == an) &
                      (R.measure == mk) & (R.outcome == k) & (R.term == term)][col]
                return float(s.iloc[0]) if len(s) and np.isfinite(s.iloc[0]) else np.nan

            cells_txt = []
            for ci, _, _ in M["cats"]:
                r = w[w.conv_cat == ci]
                if not len(r) or not np.isfinite(r.hr.iloc[0]):
                    cells_txt.append(f"{'-':>22}")
                else:
                    r = r.iloc[0]
                    cells_txt.append(f"{r.hr:>8.2f} ({r.lo:.2f},{r.hi:.2f})".rjust(22))
            pv = pick("ordinal_product", "p")
            po = pick("omnibus_wald", "p")
            md = pick("ordinal_product", "mde_hr")
            mc = pick("corner_ratio_of_HR", "mde_hr")
            mark = " *" if k in NEGATIVE_CONTROLS else ""
            ev = int(w.events_in_fit.iloc[0])

            def f(v, w_, d=3):
                return f"{v:>{w_}.{d}f}" if np.isfinite(v) else f"{'-':>{w_}}"

            log(f"{(LABEL[k] + mark)[:25]:<26}{ev:>7,}" + "".join(cells_txt) +
                f(pv, 10) + f(po, 11) + f(md, 14, 2) + f(mc, 16, 2))

    log("\n" + "=" * 100)
    log("3. THE KEY NUMBER. Inside every conventional stratum, T90 >10% against T90 <=1%.")
    log("   Hazard ratio (95% CI), adjusted for age, sex and hospital. Then the two interaction "
        "P values and,\n   beside each, what this outcome could actually have detected at 80% "
        "power: min det/step is the\n   smallest per-category change in the T90 hazard ratio "
        "the trend test could find, min det corner the\n   smallest ratio between the T90 "
        "hazard ratio in the worst stratum and in the most favourable one.")
    log("=" * 100)
    big = (R[(R.row_type == "within_stratum") & (R.analysis == "primary_4band") &
             (R.measure == "ahi") & R.in_panel_52.fillna(False)]
           .drop_duplicates("outcome").sort_values("events_in_fit", ascending=False))
    show = [k for k in big.outcome.tolist() if k not in NEGATIVE_CONTROLS][:14]
    show += [k for k in NEGATIVE_CONTROLS]
    for mk in MEASURES:
        wtab("primary_4band", mk, show,
             f"{MEASURES[mk]['name']} strata, 14 largest outcomes then the five negative "
             f"controls (*)")

    log("\n" + "=" * 100)
    log(f"3b. ACROSS THE {len(PANEL_52)}. Median within-stratum hazard ratio and how often it clears 1.")
    log("=" * 100)
    for an in ANALYSES:
        log(f"\n--- {an}")
        log(f"{'stratum':<34}{'median HR':>11}{'IQR':>20}{'CI above 1':>13}"
            f"{'q<0.05':>9}{'controls, median':>18}{'max':>8}")
        for mk, M in MEASURES.items():
            for ci, clab, _ in M["cats"]:
                r = S[(S.row_type == "summary_within_stratum") & (S.analysis == an) &
                      (S.measure == mk) & (S.conv_cat == ci)].iloc[0]
                log(f"{clab:<34}{r.hr:>11.2f}{f'{r.lo:.2f} to {r.hi:.2f}':>20}"
                    f"{f'{int(r.n_ci_above_1)}/{int(r.n_outcomes_noncontrol)}':>13}"
                    f"{int(r.n_q_main_lt_005):>9}{r.control_hr_median:>18.2f}"
                    f"{r.control_hr_max:>8.2f}")

    log("\n" + "=" * 100)
    log("4. INTERACTION. Read the minimum detectable column before reading the P column.")
    log("   A nonsignificant interaction here is NOT evidence that the T90 hazard ratio is the "
        "same in every\n   stratum. It is a statement that a difference smaller than the "
        "minimum detectable effect would\n   have been invisible. The within-stratum intervals "
        "above are the evidence.")
    log("=" * 100)
    for an in ANALYSES:
        log(f"\n--- {an}")
        log(f"{'measure':<24}{'test':<22}{'outcomes':>9}{'P<0.05':>8}{'q<0.05':>8}"
            f"{'median P':>10}{'min P':>9}{'median min detectable':>23}")
        for mk, M in MEASURES.items():
            for tm in ("ordinal_product", "omnibus_wald", "corner_ratio_of_HR"):
                r = S[(S.row_type == "summary_interaction") & (S.analysis == an) &
                      (S.measure == mk) & (S.term == tm)]
                if not len(r):
                    continue
                r = r.iloc[0]
                md = f"{r.mde_hr:.2f}" if np.isfinite(r.mde_hr) else "-"
                log(f"{M['name']:<24}{tm:<22}{int(r.n_outcomes):>9}{int(r.n_p_lt_005):>8}"
                    f"{int(r.n_q_main_lt_005):>8}{r.p:>10.3f}{r.min_p:>9.4f}{md:>23}")
    strong = R[(R.row_type == "interaction") & (R.term == "ordinal_product") &
               R.in_panel_52.fillna(False) & (R.q_main < 0.05)]
    if len(strong):
        log("\ninteractions surviving Benjamini-Hochberg (q_main < 0.05):")
        for _, r in strong.sort_values("p").iterrows():
            log(f"  {r.analysis:<14}{r.measure_name:<24}{r.label:<26}"
                f"HR per step {r.hr:.3f} ({r.lo:.3f},{r.hi:.3f})  P={r.p:.2e}  "
                f"q={r.q_main:.3f}")
    else:
        log("\nno ordinal interaction survives Benjamini-Hochberg anywhere in the panel. "
            "Given the minimum\ndetectable effects above, that is a limit on what this study "
            "can see, not a demonstration of\nuniformity.")

    log("\n" + "=" * 100)
    log("5. SENSITIVITY. T90 >10% against <=10%, so every participant enters the contrast.")
    log("=" * 100)
    n_lo = int((b.cat_primary_4band == "t1").sum())
    n_hi = int((b.cat_primary_4band == "t4").sum())
    log(f"in the four-band analysis all {len(b):,} patients are in the model but only "
        f"{n_lo + n_hi:,} of them\n({100 * (n_lo + n_hi) / len(b):.1f}%) carry the key "
        f"contrast, the {n_lo:,} at T90 <=1% and the {n_hi:,} above 10%.\n"
        f"The binary recode puts all {len(b):,} into it.")
    log(f"\n{'stratum':<34}{'4-band HR':>22}{'binary HR':>22}{'shift':>9}")
    for mk, M in MEASURES.items():
        for ci, clab, _ in M["cats"]:
            p4 = S[(S.row_type == "summary_within_stratum") & (S.analysis == "primary_4band") &
                   (S.measure == mk) & (S.conv_cat == ci)].iloc[0]
            sb = S[(S.row_type == "summary_within_stratum") & (S.analysis == "sens_binary") &
                   (S.measure == mk) & (S.conv_cat == ci)].iloc[0]
            log(f"{clab:<34}{f'{p4.hr:.2f} [{p4.lo:.2f}-{p4.hi:.2f}]':>22}"
                f"{f'{sb.hr:.2f} [{sb.lo:.2f}-{sb.hi:.2f}]':>22}"
                f"{100 * (sb.hr / p4.hr - 1):>8.1f}%")
    log("median hazard ratio across the non-control outcomes, interquartile range in brackets")

    log("\n" + "=" * 100)
    log("2. JOINT CATEGORIES. Every cell against T90 <=1% with the most favourable "
        "conventional category.")
    log("=" * 100)
    for mk, M in MEASURES.items():
        for k in ("death", "cvd", "hf", "diabetes"):
            j = R[(R.row_type == "joint_hr") & (R.analysis == "primary_4band") &
                  (R.measure == mk) & (R.outcome == k)]
            if not len(j):
                continue
            log(f"\n{LABEL[k]}, {M['name']} x T90 "
                f"(reference {dict((c, l) for c, l, _ in M['cats'])[M['ref']]} & T90 <=1%)"
                + ("" if k in PANEL_52 else "   [outside the 52, carried for reference]"))
            log(f"{'':<16}" + "".join(f"{l:>22}" for _, l, _ in
                                      T90_BANDS["primary_4band"]))
            for ci, clab, _ in M["cats"]:
                txt = []
                for tj, _, _ in T90_BANDS["primary_4band"]:
                    r = j[(j.conv_cat == ci) & (j.t90_cat == tj)]
                    if not len(r) or not np.isfinite(r.hr.iloc[0]):
                        txt.append(f"{'-':>22}")
                    elif bool(r.is_reference.iloc[0]):
                        txt.append(f"{'1.00 (reference)':>22}")
                    else:
                        r = r.iloc[0]
                        txt.append(f"{r.hr:>7.2f} ({r.lo:.2f},{r.hi:.2f})".rjust(22))
                log(f"{clab:<16}" + "".join(txt))

    # ---- headline numbers -----------------------------------------------------------------
    log("\n" + "=" * 100)
    log("HEADLINE")
    log("=" * 100)
    head = {}
    for an in ANALYSES:
        for mk, M in MEASURES.items():
            w = S[(S.row_type == "summary_within_stratum") & (S.analysis == an) &
                  (S.measure == mk)]
            worst_row = w.loc[w.hr.idxmin()]
            best_row = w.loc[w.hr.idxmax()]
            head[f"{an}|{mk}"] = {
                "n_strata": int(len(w)),
                "n_strata_median_hr_above_1": int((w.hr > 1).sum()),
                "lowest_stratum": {"stratum": worst_row.conv_label,
                                   "median_hr": round(float(worst_row.hr), 3),
                                   "n_ci_above_1": int(worst_row.n_ci_above_1),
                                   "n_outcomes": int(worst_row.n_outcomes_noncontrol)},
                "highest_stratum": {"stratum": best_row.conv_label,
                                    "median_hr": round(float(best_row.hr), 3),
                                    "n_ci_above_1": int(best_row.n_ci_above_1),
                                    "n_outcomes": int(best_row.n_outcomes_noncontrol)},
            }
            log(f"{an:<14}{M['name']:<24}median T90 hazard ratio runs "
                f"{worst_row.hr:.2f} ({worst_row.conv_label}) to {best_row.hr:.2f} "
                f"({best_row.conv_label}), above 1 in "
                f"{int((w.hr > 1).sum())} of {len(w)} strata")
    meta["headline"] = head
    meta["runtime_s"] = round(time.time() - t0, 1)
    meta["n_rows_csv"] = int(len(out))
    json.dump(meta, open(f"{HERE}/joint_summary.json", "w"), indent=1, default=str)
    log(f"\nwrote joint_summary.json. total {hms(time.time() - t0)}")


if __name__ == "__main__":
    main()
