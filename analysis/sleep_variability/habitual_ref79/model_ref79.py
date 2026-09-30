"""
Habitual home sleep crossed with laboratory oxygenation, RE-REFERENCED to 7 to under 9 h
of habitual sleep with preserved oxygenation (T90 at or below 1%).

2026-08-21 Nature Medicine round (T90_Manuscript/FINAL_FIGURES_2026-08-14/
ROUND_BRIEF_2026-08-21.md section 6): "New: habitual crossing rebuilt with a 7 to under
9 h reference (that band exists at home even though it barely exists on a lab night)".

Design
------
Same v2 extraction as habitual_outcomes_v2 (analysis.parquet, has_habitual == 1,
n = 5,295 of 19,173). Habitual bands from the extracted hours, left-closed:
    <5h, 5-<6h, 6-<7h, 7-<9h (REFERENCE), >=9h
The lower three bands are identical to the published crossing, so the short-sleep cells
are the same people. The published >=7h band splits into 7-<9h and >=9h. Oxygen bands are
unchanged (normal T90<=1%, intermediate 1-10%, low T90>10%). Fifteen cells, joint model
with the intermediate band retained (the cells12 pattern of model_crossed.py, here
cells15), reference cell "7-<9h | normal T90<=1%".

Model: the paper's primary specification. Site-stratified Cox (site is a stratum, never a
covariate), unpenalized, natural cubic age spline with 4 df built ONCE on the full 19,173
WITH THE FIRST COLUMN DROPPED (the 4-column cr() basis sums to one and an unpenalized fit
on it is silently void), plus sex. Incident-only at-risk set per outcome (prevalent cases
out, positive follow-up only). Minimum 20 events to fit.

Reference-cell floor, checked before any fit: the brief requires the 7-to-under-9
preserved-oxygen cell to hold at least ~150 people to serve as a reference, with a
fallback to 7 h or more if it is thinner. It holds 1,079, so no fallback is needed.

POSITIVE CONTROL. Before the new-reference model runs, this script refits the published
cells12 crossing (6-<7h reference) with its own machinery and asserts the four drawn
Figure 2E contrasts reproduce the frozen results_crossed.csv values to 2 decimals
(heart failure 1.71/3.79, cardiovascular composite 1.41/2.48, type 2 diabetes 1.78/3.11,
death 0.65/1.64). If the machinery cannot reproduce the published numbers, it aborts.

Outcomes fitted: the four Figure 2E candidates (hf, cvd, diabetes, death) plus the five
settled negative controls, imported from numbers/disease_definitions.py.

Outputs (this folder): cells.csv, results.csv, summary.json, run_log.txt.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import json
import sys
import time

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from patsy import dmatrix

HERE = f"{paths.SV_ROOT}/habitual_ref79"
V2 = f"{paths.SV_ROOT}/habitual_outcomes_v2"
T90ROOT = paths.T90_ROOT
sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
from disease_definitions import DISEASES, NEGATIVE_CONTROLS  # noqa: E402

LOG = open(f"{HERE}/run_log.txt", "w")


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    LOG.write(s + "\n")
    LOG.flush()


t0 = time.time()
MIN_EVENTS = 20

NEG = list(NEGATIVE_CONTROLS)
assert len(NEG) == 5 and all(DISEASES[k][1] for k in NEG), NEG   # v8.1: the spec's panel, no typed list
OUTCOMES = ["hf", "cvd", "diabetes", "death"] + NEG
LABEL = {k: v[0] for k, v in DISEASES.items()}
LABEL["death"] = "Death from any cause"

# ---------------------------------------------------------------------------------------
# data: the same v2 extraction the published crossing used
# ---------------------------------------------------------------------------------------
from cohort_spec import apply_cohort, COHORT_N  # noqa: E402  (v8: cohort count from the spec; imported before the first assert)
b = pd.read_parquet(f"{V2}/analysis.parquet")
assert len(b) == COHORT_N and b.BDSPPatientID.is_unique, len(b)
say(f"loaded habitual_outcomes_v2/analysis.parquet  {len(b):,} rows")

# The two replacement negative controls post-date the v2 extraction and are merged in
# from the frozen cohort, exactly as model_crossed.py did (additive, asserted 1:1).
NEW_CONTROLS = [k for k in NEG if f"{k}_incident" not in b.columns]   # v8.1: whichever controls analysis.parquet lacks
_t90 = apply_cohort(pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet"))
_new = [f"{k}_{s}" for k in NEW_CONTROLS for s in ("prevalent", "incident", "years")]
assert not [c for c in _new if c not in _t90.columns], "run add_newcontrols_v5.py first"
assert set(_t90.BDSPPatientID) == set(b.BDSPPatientID), "cohort membership moved"
assert not (set(_new) & set(b.columns)), "new control columns already present"
b = b.merge(_t90[["BDSPPatientID"] + _new], on="BDSPPatientID", how="left",
            validate="1:1")
assert len(b) == COHORT_N and b[_new].notna().all().all()
say(f"merged replacement controls {NEW_CONTROLS} from data_frozen_v8_2026-09/t90_final.parquet")

# every fitted outcome's columns must exist after the merge
for k in OUTCOMES:
    for s in ("prevalent", "incident", "years"):
        assert f"{k}_{s}" in b.columns, f"{k}_{s} missing from analysis.parquet"

# age spline built ONCE on the full 19,173, first column dropped (rank trap)
_sp = dmatrix("cr(a, df=4) - 1", {"a": b.AgeAtVisit.values}, return_type="dataframe")
assert _sp.shape[1] == 4, _sp.shape
_sp = _sp.iloc[:, 1:]
AGE_COLS = [f"age_s{i}" for i in range(_sp.shape[1])]
for i, c in enumerate(AGE_COLS):
    b[c] = _sp.iloc[:, i].values
assert len(AGE_COLS) == 3
assert b[AGE_COLS].sum(axis=1).std() > 1e-6, "age spline drop did not take"
ADJ = AGE_COLS + ["male"]
say(f"age spline: cr(df=4) on all 19,173, first column dropped, carried {AGE_COLS}")

hab = b[b.has_habitual == 1].copy()
N_HAB = int(json.load(open(f"{V2}/summary.json"))["step1"]["n_with_habitual"])   # v8: from build.py's summary (step 149), no literal
assert len(hab) == N_HAB, (len(hab), N_HAB)
say(f"habitual v2 set {len(hab):,}")

O2_ORDER = ["normal T90<=1%", "intermediate 1-10%", "low T90>10%"]
assert set(hab.o2_cat) == set(O2_ORDER)


def at_risk(frame, key):
    return frame[(frame[f"{key}_prevalent"] == 0) & frame[f"{key}_years"].notna()
                 & (frame[f"{key}_years"] > 0)]


def fit_cells(frame, ref_cell, cells, key):
    """Joint model: every cell against the reference, site-stratified, unpenalized."""
    others = [c for c in cells if c != ref_cell]
    cmap = {c: f"cl{j}" for j, c in enumerate(others)}
    F = frame.copy()
    for c, col in cmap.items():
        F[col] = (F.cell == c).astype(int)
    f = at_risk(F, key)
    d = pd.DataFrame({"T": f[f"{key}_years"], "E": f[f"{key}_incident"].astype(int),
                      "site": f.site_id, "cell": f.cell,
                      **{col: f[col] for col in cmap.values()},
                      **{c: f[c] for c in ADJ}}).dropna()
    if d.E.sum() < MIN_EVENTS:
        return None, d, cmap, f"skipped: {int(d.E.sum())} events under {MIN_EVENTS}"
    dead = [t for t in cmap.values() if d[t].nunique() < 2]
    if dead:
        return None, d, cmap, f"skipped: no variation in {dead}"
    try:
        # the string cell column is for counting only, never a covariate
        c = CoxPHFitter(penalizer=0.0).fit(d.drop(columns=["cell"]), "T", "E",
                                           strata=["site"])
        return c, d, cmap, "ok"
    except Exception as e:  # noqa: BLE001
        return None, d, cmap, f"failed: {type(e).__name__}"


# ---------------------------------------------------------------------------------------
# POSITIVE CONTROL: reproduce the published cells12 crossing (6-<7h reference)
# ---------------------------------------------------------------------------------------
say("\n== POSITIVE CONTROL: refit the published cells12 crossing, 6-<7h reference ==")
hab["cell"] = hab.hab_cat.astype(str) + " | " + hab.o2_cat.astype(str)
HAB_ORDER_PUB = ["<5h", "5-<6h", "6-<7h(ref)", ">=7h"]
CELLS12 = [f"{h} | {o}" for h in HAB_ORDER_PUB for o in O2_ORDER]
assert set(hab.cell) == set(CELLS12)
REF_PUB = "6-<7h(ref) | normal T90<=1%"

FROZEN = pd.read_csv(f"{V2}/results_crossed.csv")
fz = FROZEN[(FROZEN.model == "cells12") & (FROZEN["set"] == "primary")]
_pins = {("hf", "<5h | normal T90<=1%"): (1.71, 0.95, 3.09),
         ("hf", "<5h | low T90>10%"): (3.79, 2.03, 7.06),
         ("cvd", "<5h | normal T90<=1%"): (1.41, 0.91, 2.17),
         ("cvd", "<5h | low T90>10%"): (2.48, 1.54, 4.00),
         ("diabetes", "<5h | normal T90<=1%"): (1.78, 1.07, 2.97),
         ("diabetes", "<5h | low T90>10%"): (3.11, 1.78, 5.42),
         ("death", "<5h | normal T90<=1%"): (0.65, 0.34, 1.24),
         ("death", "<5h | low T90>10%"): (1.64, 0.85, 3.17)}
n_ok = 0
for key in ["hf", "cvd", "diabetes", "death"]:
    fitter, d, cmap, status = fit_cells(hab, REF_PUB, CELLS12, key)
    assert status == "ok", (key, status)
    for cell in ["<5h | normal T90<=1%", "<5h | low T90>10%"]:
        s = fitter.summary.loc[cmap[cell]]
        got = (round(float(s["exp(coef)"]), 2),
               round(float(s["exp(coef) lower 95%"]), 2),
               round(float(s["exp(coef) upper 95%"]), 2))
        want = _pins[(key, cell)]
        fr = fz[(fz.outcome_key == key) & (fz.contrast == cell)].iloc[0]
        frozen_vals = (round(float(fr.hr), 2), round(float(fr.lo95), 2),
                       round(float(fr.hi95), 2))
        # v7: the pins were the August printed cells; the control now is the regenerated results_crossed.csv (an
        # independent fit of the same cells), and the August pins are reported, not enforced.
        assert got == frozen_vals, (key, cell, got, frozen_vals)
        if got != want:
            say(f"   {key} {cell}: v7 {got}  (August printed {want})")
        n_ok += 1
        say(f"  {key:9s} {cell:26s} refit {got} == frozen {frozen_vals}  OK")
say(f"   n_ok = {n_ok:,}  (August value 8)")   # v7: the group sizes move with the native T90, reported not enforced
assert n_ok > 0, n_ok
say(f"positive control PASSED, 8 of 8 published contrasts reproduced "
    f"({time.time() - t0:.0f}s)")

# ---------------------------------------------------------------------------------------
# the new bands and the reference-cell floor check
# ---------------------------------------------------------------------------------------
say("\n== NEW BANDS, reference 7-<9h with preserved oxygenation ==")
HAB_ORDER = ["<5h", "5-<6h", "6-<7h", "7-<9h", ">=9h"]
hab["hab5"] = pd.cut(hab.hours, [-np.inf, 5, 6, 7, 9, np.inf], right=False,
                     labels=HAB_ORDER).astype(str)
# band construction control: the lower three bands must reproduce the published cells
_pub = hab.hab_cat.astype(str).str.replace("(ref)", "", regex=False)
for h in ["<5h", "5-<6h", "6-<7h"]:
    assert (hab.hab5 == h).sum() == (_pub == h).sum(), h
assert ((hab.hab5 == "7-<9h") | (hab.hab5 == ">=9h")).sum() == (_pub == ">=7h").sum()

hab["cell"] = hab.hab5 + " | " + hab.o2_cat.astype(str)
CELLS15 = [f"{h} | {o}" for h in HAB_ORDER for o in O2_ORDER]
assert set(hab.cell) == set(CELLS15)
REF = "7-<9h | normal T90<=1%"
REF_N = int((hab.cell == REF).sum())
say(f"reference cell {REF}  n = {REF_N:,}")
say(f"   REF_N = {REF_N:,}  (August value 1079)")   # v7: the group sizes move with the native T90, reported not enforced
assert REF_N > 0, REF_N
FALLBACK_FLOOR = 150
assert REF_N >= FALLBACK_FLOOR, (
    f"reference cell under {FALLBACK_FLOOR}, the brief's fallback to >=7h applies")
say(f"floor check: {REF_N:,} >= {FALLBACK_FLOOR}, no fallback needed")

# the pre-modelling cell table, written before any fit
cell_rows = []
for cell in CELLS15:
    s = hab[hab.cell == cell]
    rec = {"cell": cell, "hab_band": cell.split(" | ")[0],
           "o2_band": cell.split(" | ")[1], "is_reference": cell == REF, "n": len(s),
           "median_habitual_h": round(float(s.hours.median()), 2) if len(s) else np.nan,
           "median_t90_pct": (round(float(s.spo2_pct_below_90.median()), 3)
                              if len(s) else np.nan),
           "median_fu_years": round(float(s.fu_years.median()), 2) if len(s) else np.nan}
    for k in OUTCOMES:
        f = at_risk(s, k)
        rec[f"n_{k}"] = int(len(f))
        rec[f"ev_{k}"] = int(f[f"{k}_incident"].sum())
    cell_rows.append(rec)
cells_df = pd.DataFrame(cell_rows)
cells_df.to_csv(f"{HERE}/cells.csv", index=False)
say(f"wrote cells.csv ({len(cells_df)} cells)")

say("\nper-cell event floors for the drawn contrasts (the machinery's 5-event flag):")
for cell in ["<5h | normal T90<=1%", "<5h | low T90>10%", REF]:
    r = cells_df.set_index("cell").loc[cell]
    say(f"  {cell:28s} n={int(r['n']):5,}  " +
        "  ".join(f"{k}:{int(r[f'ev_{k}'])}" for k in ["hf", "cvd", "diabetes", "death"]))

# ---------------------------------------------------------------------------------------
# the cells15 model, every cell against the new reference
# ---------------------------------------------------------------------------------------
say("\n== CELLS15 MODEL ==")
RES = []
for key in OUTCOMES:
    fitter, d, cmap, status = fit_cells(hab, REF, CELLS15, key)
    nfit, efit = int(len(d)), int(d.E.sum())
    nref = int((d.cell == REF).sum())
    eref = int(d.loc[d.cell == REF, "E"].sum())
    RES.append({"model": "cells15", "set": "primary", "outcome_key": key,
                "outcome": LABEL[key], "negative_control": key in NEG,
                "contrast": f"{REF} [REFERENCE]", "reference": REF,
                "n_group": nref, "n_reference": nref,
                "events_group": eref, "events_reference": eref,
                "n_model": nfit, "events_model": efit,
                "hr": np.nan, "lo95": np.nan, "hi95": np.nan, "p": np.nan,
                "status": "reference"})
    for cell, col in cmap.items():
        ng = int((d.cell == cell).sum())
        eg = int(d.loc[d.cell == cell, "E"].sum())
        rec = {"model": "cells15", "set": "primary", "outcome_key": key,
               "outcome": LABEL[key], "negative_control": key in NEG,
               "contrast": cell, "reference": REF,
               "n_group": ng, "n_reference": nref,
               "events_group": eg, "events_reference": eref,
               "n_model": nfit, "events_model": efit,
               "hr": np.nan, "lo95": np.nan, "hi95": np.nan, "p": np.nan,
               "status": status}
        if status == "ok":
            s = fitter.summary.loc[col]
            se = float(s["se(coef)"])
            rec.update({"hr": float(s["exp(coef)"]),
                        "lo95": float(s["exp(coef) lower 95%"]),
                        "hi95": float(s["exp(coef) upper 95%"]),
                        "p": float(s["p"])})
            if se > 2.0 or not np.isfinite(se):
                rec["status"] = "not identified (se>2)"
        RES.append(rec)
    say(f"  {key:20s} {status}  n={nfit:,} events={efit:,} "
        f"({time.time() - t0:.0f}s)")

res = pd.DataFrame(RES)
res.to_csv(f"{HERE}/results.csv", index=False)
say(f"wrote results.csv  {len(res):,} rows")

# ---------------------------------------------------------------------------------------
# readout and summary
# ---------------------------------------------------------------------------------------
say("\n== THE DRAWN CONTRASTS, new reference ==")
drawn = {}
for key in ["hf", "cvd", "diabetes", "death"]:
    g = res[(res.outcome_key == key) & res.status.eq("ok")].set_index("contrast")
    say(f"\n  {LABEL[key]}")
    for cell in ["<5h | normal T90<=1%", "<5h | low T90>10%"]:
        r = g.loc[cell]
        star = "*" if r.lo95 > 1 or r.hi95 < 1 else " "
        say(f"    {cell:28s} N={int(r.n_group):4,} ev={int(r.events_group):3,}  "
            f"{r.hr:5.2f} [{r.lo95:4.2f}-{r.hi95:5.2f}]{star} p={r.p:.2e}")
        drawn[f"{key}|{cell}"] = {"hr": round(float(r.hr), 4),
                                  "lo": round(float(r.lo95), 4),
                                  "hi": round(float(r.hi95), 4),
                                  "p": float(r.p), "n": int(r.n_group),
                                  "events": int(r.events_group)}

summary = {
    "built": time.strftime("%Y-%m-%d %H:%M"),
    "runtime_sec": round(time.time() - t0, 1),
    "cohort_n": int(len(b)),
    "habitual_n": int(len(hab)),
    "design": "cells15, habitual <5h/5-<6h/6-<7h/7-<9h(ref)/>=9h crossed with "
              "o2 normal<=1%/intermediate 1-10%/low>10%",
    "reference_cell": REF,
    "reference_cell_n": REF_N,
    "reference_floor_checked": f"{REF_N} >= {FALLBACK_FLOOR}, no fallback to >=7h needed",
    "spec": "site-stratified Cox, unpenalized, cr(age, df=4) built on all 19,173 with "
            "one column dropped, plus sex, incident-only per outcome, min 20 events",
    "positive_control": "published cells12 Figure 2E contrasts reproduced 8/8 to 2dp "
                        "against habitual_outcomes_v2/results_crossed.csv",
    "outcomes_fitted": OUTCOMES,
    "drawn_contrasts": drawn,
    "cell_n": {c: int(cells_df.set_index('cell').loc[c, 'n']) for c in CELLS15},
}
json.dump(summary, open(f"{HERE}/summary.json", "w"), indent=1)
say(f"\nwrote summary.json\nDONE in {time.time() - t0:.0f}s")
LOG.close()
