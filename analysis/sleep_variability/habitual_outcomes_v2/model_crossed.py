"""
Habitual home sleep crossed with laboratory oxygenation.

The question this answers, in Alen's words: "These people sleep less than five hours and have
normal oxygenation, what is their risk? These people sleep less than five hours and have terrible
oxygenation, how much is their risk? ... Just so we can understand whether in this cohort, where
you have habitual sleep and lab oxygenation, whether oxygen holds up or not."

Design
------
Exposure A, habitual home sleep from the clinical note, four bands, left-closed:
    <5h, 5-<6h, 6-<7h (reference), >=7h
Exposure B, laboratory oxygenation, T90 as percent of the WHOLE recording below 90%:
    normal <=1%, intermediate 1-10%, low >10%

The primary crossed model uses the two extreme oxygen bands only, so eight cells, reference
6-<7h with normal oxygenation. The intermediate band is not thrown away: a twelve-cell version
is fitted alongside and written to the same file under model="cells12".

Every fit is a Cox proportional hazards model, stratified by site (site is a stratum, never a
covariate, because capture rate and measurement quality both track site), unpenalized, adjusted
for a natural cubic age spline and sex.

THE AGE SPLINE HAS FOUR COLUMNS AND ONE IS DROPPED. A four-column cr() basis sums to one and is
rank deficient. An unpenalized Cox fit on a rank-deficient design silently fits nothing and
returns zero fitted outcomes. The basis is built once on the full 19,173 so the knots are
identical in every subset fit and the strata are comparable.

Model families written to results_crossed.csv
---------------------------------------------
cells8            every cell against 6-<7h with normal oxygenation, joint model, 7 contrasts
cells12           the same with the intermediate oxygen band retained, 11 contrasts
o2_within_hab     WITHIN each habitual band, low versus normal oxygenation. The direct test of
                  whether oxygen holds up regardless of how long people sleep at home.
hab_within_o2     WITHIN each oxygen band, each habitual band against 6-<7h. The mirror test.
                  Also run across all oxygen groups pooled, which is the plain habitual sleep
                  and disease question with no oxygen conditioning.
interaction       formal habitual x oxygen interaction, likelihood ratio test. Two forms: the
                  3 df categorical test on the eight cells, and a 1 df continuous test with more
                  power. Power is poor for both, so the within-stratum estimates and their
                  confidence intervals carry the argument, not these p values.
continuous        per-SD hazard ratios for habitual hours, laboratory T90, laboratory total sleep
                  time and AHI, on a rank-based inverse normal transform within site. This is the
                  head-to-head Alen asked for, "the comparison obviously is how total sleep time
                  does as well", plus the published cohort anchor.

Every family is also refitted on the 4,631 patients whose habitual value was recorded BEFORE the
sleep study (set="pre_study"). 46.8% of the primary set had the sleep duration written down after
the outcome clock started, which runs in the direction that creates an association, so a positive
finding that does not reproduce in the pre-study set is not believed.

Minimum detectable hazard ratio
-------------------------------
Reported on every row, two ways, because with 6,074 patients split across eight cells most nulls
here are uninformative rather than negative and must be labelled as such.
    mdhr_analytic   from the event split alone, exp(2.802 / sqrt(D * p * (1-p))), available even
                    for cells too small to fit
    mdhr_empirical  exp(2.802 * observed standard error), the design-exact figure, which accounts
                    for the age and sex adjustment and the site strata
2.802 = 1.960 + 0.842, two-sided alpha 0.05 and 80% power.

Negative controls and the floor
-------------------------------
All five current controls carry: back pain, cataract, glaucoma, contact dermatitis, hemorrhoids.
The panel is imported from numbers/disease_definitions.py rather than retyped.

THE PANEL CHANGED WHILE THIS RAN. The brief named six controls including osteoarthritis, with a
per-SD floor of 1.133. disease_definitions.py was rewritten on 2026-08-07 to five controls with
osteoarthritis removed, on the ground that published work links osteoarthritis to hypoxia so a
causal path from the exposure exists. Fracture is out for the same class of reason. Both remain
ordinary outcomes here. The two replacements were written into data_frozen_v7_2026-09/t90_final.parquet at
16:41 that day, after this project's analysis.parquet was built at 16:20, so they are merged in
at load time. The per-SD floor is now 1.063, set by glaucoma, read from
numbers/negcontrols_v4_panel.csv. The brief's 1.133 is carried in parallel as per_sd_floor_legacy
and reported in brackets, so the change is visible rather than silent.

For the cell contrasts a per-contrast floor is computed from the controls themselves, since a
categorical cell contrast is not on the per-SD scale. Two versions are written. The symmetric
floor, max(HR, 1/HR) over the controls, is conservative but hostage to one unstable estimate.
The directional pair, the largest and smallest control HR, is the one to read: a harmful
association must clear the upper floor and a protective one must fall below the lower floor.

Outputs
-------
results_crossed.csv     every estimate, long form
cells_crossed.csv       the pre-modelling cell table, N and events, written before any fit
floors_crossed.csv      the negative-control floor for every contrast
summary_crossed.json    headline numbers
model_crossed_log.txt   full transcript
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import json
import os
import sys
import time
import warnings

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from patsy import dmatrix
from scipy import stats

warnings.filterwarnings("ignore")

HERE = f"{paths.SV_ROOT}/habitual_outcomes_v2"
T90ROOT = paths.T90_ROOT
sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
from cohort_spec import apply_cohort, COHORT_N  # noqa: E402  (v8: cohort count from the spec)
from disease_definitions import (DISEASES, ORGAN_GROUP, CIRCULAR,  # noqa: E402
                                 NEGATIVE_CONTROLS)

# ---------------------------------------------------------------------------------------
# standing rules. The control panel is IMPORTED, never retyped, so an upstream change
# surfaces here instead of drifting silently.
#
# It changed under this analysis on 2026-08-07. The brief for this run named six controls
# (back pain, cataract, glaucoma, osteoarthritis, contact dermatitis, hemorrhoids) with a
# per-SD floor of 1.133. numbers/disease_definitions.py was rewritten the same afternoon and
# now names five: back pain, cataract, glaucoma, contact dermatitis, hemorrhoids.
#   OSTEOARTHRITIS IS NO LONGER A CONTROL. Published work links it to hypoxia, so a causal
#   path from the exposure exists. It stays an ordinary outcome.
#   FRACTURE IS NO LONGER A CONTROL either, and was already out under the brief.
#   The two replacements were materialised into data_frozen_v7_2026-09/t90_final.parquet at 16:41 on
#   2026-08-07 by numbers/add_newcontrols_v5.py, after this project's analysis.parquet was
#   built at 16:20, so they are merged in below rather than rebuilt.
# The per-SD floor moves with the panel: 1.063, set by glaucoma, read from
# numbers/negcontrols_v4_panel.csv. The brief's 1.133 is carried alongside as the legacy
# figure so nothing is changed silently and both readings are on the face of the file.
# ---------------------------------------------------------------------------------------
NEG_CONTROLS = list(NEGATIVE_CONTROLS)
assert len(NEG_CONTROLS) == 5 and all(DISEASES[k][1] for k in NEG_CONTROLS), NEG_CONTROLS   # v8.1: the spec's panel, no typed list
NOT_A_CONTROL = ["fracture", "osteoarthritis", "thyroid_dis"]
NEW_CONTROLS = ["dermatitis_contact", "haemorrhoids"]
PER_SD_FLOOR_LEGACY = 1.133
_panel = pd.read_csv(f"{paths.NUMBERS_DIR}/negcontrols_v4_panel.csv")
assert set(_panel.key) == set(NEG_CONTROLS), sorted(set(_panel.key) ^ set(NEG_CONTROLS))
PER_SD_FLOOR = float(_panel.hr.max())
PER_SD_FLOOR_SETBY = str(_panel.loc[_panel.hr.idxmax(), "key"])
MIN_EVENTS = 20           # a condition is fitted at this many events or more
Z_POWER = 1.959963985 + 0.8416212336   # 2.8016, two-sided 0.05 at 80% power

HAB_ORDER = ["<5h", "5-<6h", "6-<7h(ref)", ">=7h"]
HAB_REF = "6-<7h(ref)"
O2_ORDER = ["normal T90<=1%", "intermediate 1-10%", "low T90>10%"]
O2_NORMAL, O2_INTER, O2_LOW = O2_ORDER
O2_REF = O2_NORMAL
REF_CELL = f"{HAB_REF} | {O2_REF}"

OUTCOMES = list(DISEASES) + ["death"]
LABEL = {k: v[0] for k, v in DISEASES.items()}
LABEL["death"] = "Death from any cause"
GROUP = dict(ORGAN_GROUP)
GROUP["death"] = "Mortality"

LOG = open(f"{HERE}/model_crossed_log.txt", "w")


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    LOG.write(s + "\n")
    LOG.flush()


# ---------------------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------------------
t_start = time.time()
b = pd.read_parquet(f"{HERE}/analysis.parquet")
say(f"loaded analysis.parquet  {len(b):,} rows x {b.shape[1]} cols")
assert len(b) == COHORT_N, len(b)
assert b.BDSPPatientID.is_unique

# Merge the two replacement negative controls, which post-date this project's build. Additive
# only: nothing already in analysis.parquet is read back or overwritten, and the merge is
# asserted to be a complete one-to-one match on patient id before anything is fitted.
NEW_CONTROLS = [k for k in NEG_CONTROLS if f"{k}_incident" not in b.columns]   # v8.1: whichever controls analysis.parquet lacks
_t90 = apply_cohort(pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet"))
_new = [f"{k}_{s}" for k in NEW_CONTROLS for s in ("prevalent", "incident", "years")]
_missing = [c for c in _new if c not in _t90.columns]
assert not _missing, f"run numbers/add_newcontrols_v5.py first, missing {_missing}"
assert set(_t90.BDSPPatientID) == set(b.BDSPPatientID), "cohort membership moved upstream"
assert not (set(_new) & set(b.columns)), "new control columns already present, refusing to merge"
b = b.merge(_t90[["BDSPPatientID"] + _new], on="BDSPPatientID", how="left", validate="1:1")
assert len(b) == COHORT_N and b[_new].notna().all().all()
say(f"merged replacement controls {NEW_CONTROLS} from data_frozen_v8_2026-09/t90_final.parquet "
    f"-> {b.shape[1]} cols")
for k in NEW_CONTROLS:
    _f = b[(b[f"{k}_prevalent"] == 0) & (b[f"{k}_years"] > 0)]
    say(f"  {k:22s} {int(_f[f'{k}_incident'].sum()):,} incident in the 19,173, "
        f"{int(b.loc[b.has_habitual == 1, f'{k}_incident'].sum()):,} in the habitual 8,711")

# The age basis is built ONCE on all 19,173 so every subset fit uses identical knots.
# cr(a, df=4) gives four columns that sum to one. ONE IS DROPPED. Fitting the four-column
# basis unpenalized silently returns nothing.
_sp = dmatrix("cr(a, df=4) - 1", {"a": b.AgeAtVisit.values}, return_type="dataframe")
assert _sp.shape[1] == 4, _sp.shape
_sp = _sp.iloc[:, 1:]
AGE_COLS = [f"age_s{i}" for i in range(_sp.shape[1])]
for i, c in enumerate(AGE_COLS):
    b[c] = _sp.iloc[:, i].values
assert len(AGE_COLS) == 3, AGE_COLS
ADJ = AGE_COLS + ["male"]
say(f"age spline: 4-column cr(df=4) basis, first column dropped, {len(AGE_COLS)} carried -> {ADJ}")

# guard: the drop must actually have happened, or every fit below is void
_chk = b[AGE_COLS].sum(axis=1)
assert _chk.std() > 1e-6, "age spline columns still sum to a constant, the drop did not take"

hab = b[b.has_habitual == 1].copy()
say(f"habitual set {len(hab):,}   pre-study only {int(hab.pre_study_only.sum()):,}"
    f"   post-study only {int(hab.post_study_only.sum()):,}")
N_HAB = int(json.load(open(f"{HERE}/summary.json"))["step1"]["n_with_habitual"])   # v8: from build.py's summary (step 149), no literal
assert len(hab) == N_HAB, (len(hab), N_HAB)

hab["cell"] = hab.hab_cat.astype(str) + " | " + hab.o2_cat.astype(str)
CELLS12 = [f"{h} | {o}" for h in HAB_ORDER for o in O2_ORDER]
CELLS8 = [f"{h} | {o}" for h in HAB_ORDER for o in (O2_NORMAL, O2_LOW)]
assert set(hab.cell) == set(CELLS12)

X8 = hab[hab.o2_cat.isin([O2_NORMAL, O2_LOW])].copy()
say(f"eight-cell crossed set {len(X8):,}   twelve-cell set {len(hab):,}")

# rank-based inverse normal, within site, matching the published specification
def rint(x):
    r = stats.rankdata(x)
    return stats.norm.ppf((r - 0.375) / (len(r) + 0.25))


CONT = {"habitual_hours": "hours", "lab_t90": "spo2_pct_below_90",
        "lab_tst": "TST_min", "lab_ahi": "AHI"}
for frame in (b, hab, X8):
    for name, src in CONT.items():
        if src == "hours" and "hours" not in frame.columns:
            continue
        v = frame[src]
        if v.isna().all():
            frame[f"z_{name}"] = np.nan
            continue
        frame[f"z_{name}"] = frame.groupby("site_id")[src].transform(
            lambda s: rint(s.fillna(s.median())))
# habitual hours does not exist outside the habitual set
b["z_habitual_hours"] = np.nan
b.loc[hab.index, "z_habitual_hours"] = hab["z_habitual_hours"]


# ---------------------------------------------------------------------------------------
# STEP 1  the cells, N and events, BEFORE any model is fitted
# ---------------------------------------------------------------------------------------
say("\n" + "=" * 100)
say("STEP 1  CELLS BEFORE MODELLING")
say("=" * 100)

HEADLINE = ["cvd", "death", "htn2", "diabetes", "hf", "dementia"]
cell_rows = []
for design, frame, cells in (("cells8", X8, CELLS8), ("cells12", hab, CELLS12)):
    for cell in cells:
        s = frame[frame.cell == cell]
        rec = {"design": design, "cell": cell,
               "hab_band": cell.split(" | ")[0], "o2_band": cell.split(" | ")[1],
               "is_reference": cell == REF_CELL, "n": len(s),
               "median_t90_pct": round(float(s.spo2_pct_below_90.median()), 3),
               "median_lab_tst_min": round(float(s.TST_min.median()), 1),
               "median_habitual_h": round(float(s.hours.median()), 2),
               "median_fu_years": round(float(s.fu_years.median()), 2),
               "median_age": round(float(s.AgeAtVisit.median()), 1),
               "pct_male": round(100 * float(s.male.mean()), 1),
               "pct_pre_study": round(100 * float(s.pre_study_only.mean()), 1),
               "n_site_I0002": int((s.site_id == "I0002").sum()),
               "n_site_I0006": int((s.site_id == "I0006").sum())}
        for k in OUTCOMES:
            f = s[(s[f"{k}_prevalent"] == 0) & s[f"{k}_years"].notna() & (s[f"{k}_years"] > 0)]
            rec[f"ev_{k}"] = int(f[f"{k}_incident"].sum())
            rec[f"n_{k}"] = int(len(f))
        cell_rows.append(rec)
cells_df = pd.DataFrame(cell_rows)
cells_df.to_csv(f"{HERE}/cells_crossed.csv", index=False)

say(f"\nEIGHT CELLS, reference {REF_CELL}")
say(f"{'cell':40s}{'N':>7}{'T90':>8}{'labTST':>8}{'fu':>6}" +
    "".join(f"{k:>10}" for k in HEADLINE))
for cell in CELLS8:
    r = cells_df[(cells_df.design == "cells8") & (cells_df.cell == cell)].iloc[0]
    mark = " <-REF" if r.is_reference else ""
    say(f"{cell:40s}{r.n:7,}{r.median_t90_pct:8.2f}{r.median_lab_tst_min:8.1f}"
        f"{r.median_fu_years:6.2f}" + "".join(f"{int(r['ev_'+k]):10,}" for k in HEADLINE) + mark)
say(f"{'TOTAL':40s}{cells_df[cells_df.design=='cells8'].n.sum():7,}")

say(f"\nTWELVE CELLS (intermediate oxygen retained)")
for cell in CELLS12:
    r = cells_df[(cells_df.design == "cells12") & (cells_df.cell == cell)].iloc[0]
    say(f"{cell:40s}{r.n:7,}{r.median_t90_pct:8.2f}{r.median_lab_tst_min:8.1f}"
        f"{r.median_fu_years:6.2f}" + "".join(f"{int(r['ev_'+k]):10,}" for k in HEADLINE))

# which cells are too small to fit, stated up front
say("\nCELLS TOO SMALL TO FIT")
small = []
for _, r in cells_df[cells_df.design == "cells8"].iterrows():
    thin = [k for k in OUTCOMES if r[f"ev_{k}"] < 5]
    if thin:
        small.append((r.cell, len(thin)))
        say(f"  {r.cell:40s} N={r.n:,}  {len(thin)} of {len(OUTCOMES)} outcomes under 5 events: "
            f"{', '.join(thin[:8])}{' ...' if len(thin) > 8 else ''}")
if not small:
    say("  none: every one of the eight cells carries at least 5 events for all 52 outcomes")
say(f"\nsmallest cell: " +
    cells_df[cells_df.design == "cells8"].sort_values("n").iloc[0].cell +
    f" at N={cells_df[cells_df.design=='cells8'].n.min():,}")
say(f"lab total sleep time by cell spans "
    f"{cells_df[cells_df.design=='cells8'].median_lab_tst_min.min():.1f} to "
    f"{cells_df[cells_df.design=='cells8'].median_lab_tst_min.max():.1f} min, "
    f"with no ordering by habitual band")


# ---------------------------------------------------------------------------------------
# fitting helpers
# ---------------------------------------------------------------------------------------
def at_risk(frame, key):
    """Incident-case frame for one outcome: prevalent cases out, positive follow-up only."""
    return frame[(frame[f"{key}_prevalent"] == 0) & frame[f"{key}_years"].notna()
                 & (frame[f"{key}_years"] > 0)]


def analytic_mdhr(ev_grp, ev_ref, n_grp, n_ref):
    """Minimum detectable HR from the event split alone, before any model is fitted."""
    D = ev_grp + ev_ref
    n = n_grp + n_ref
    if D < 1 or n < 1:
        return np.nan
    p = n_grp / n
    v = D * p * (1 - p)
    if v <= 0:
        return np.nan
    return float(np.exp(Z_POWER / np.sqrt(v)))


def power_label(mdhr, ev_grp):
    if ev_grp is not None and ev_grp < 5:
        return "uninformative (under 5 events)"
    if not np.isfinite(mdhr):
        return "uninformative (not estimable)"
    if mdhr < 1.35:
        return "informative"
    if mdhr < 1.75:
        return "moderate"
    if mdhr < 2.50:
        return "weak"
    return "uninformative"


def fit_cox(d, terms):
    """Unpenalized, site-stratified Cox. Returns (fitter, status)."""
    keep = list(terms) + ADJ + ["T", "E", "site"]
    dd = d[keep].dropna()
    if dd.E.sum() < MIN_EVENTS:
        return None, dd, f"skipped: {int(dd.E.sum())} events, under {MIN_EVENTS}"
    # a term with no variation cannot be fitted
    dead = [t for t in terms if dd[t].nunique() < 2]
    if dead:
        return None, dd, f"skipped: no variation in {dead}"
    try:
        c = CoxPHFitter(penalizer=0.0).fit(dd, "T", "E", strata=["site"])
        return c, dd, "ok"
    except Exception as e:  # noqa: BLE001
        return None, dd, f"failed: {type(e).__name__}"


def row(model, setname, key, contrast, reference, n_grp, n_ref, ev_grp, ev_ref,
        n_fit, ev_fit, fitter=None, term=None, status="ok", extra=None):
    r = {"model": model, "set": setname, "outcome_key": key, "outcome": LABEL[key],
         "organ_group": GROUP[key],
         "negative_control": key in NEG_CONTROLS, "circular": key in CIRCULAR,
         "contrast": contrast, "reference": reference,
         "n_group": n_grp, "n_reference": n_ref,
         "events_group": ev_grp, "events_reference": ev_ref,
         "n_model": n_fit, "events_model": ev_fit,
         "hr": np.nan, "lo95": np.nan, "hi95": np.nan, "se_loghr": np.nan, "p": np.nan,
         "chi2": np.nan, "df": np.nan,
         "mdhr_analytic": np.nan, "mdhr_empirical": np.nan, "power": "", "epv": np.nan,
         "status": status}
    if ev_grp is not None and ev_ref is not None:
        r["mdhr_analytic"] = analytic_mdhr(ev_grp, ev_ref, n_grp, n_ref)
    if fitter is not None and term is not None:
        s = fitter.summary.loc[term]
        se = float(s["se(coef)"])
        r.update({"hr": float(s["exp(coef)"]), "lo95": float(s["exp(coef) lower 95%"]),
                  "hi95": float(s["exp(coef) upper 95%"]), "se_loghr": se, "p": float(s["p"])})
        r["mdhr_empirical"] = float(np.exp(Z_POWER * se))
        # a standard error this large means the contrast is not identified, not that it is null
        if se > 2.0 or not np.isfinite(se):
            r["status"] = "not identified (se>2)"
        n_par = len(fitter.params_)
        r["epv"] = round(ev_fit / n_par, 2) if n_par else np.nan
    r["power"] = power_label(r["mdhr_empirical"] if np.isfinite(r["mdhr_empirical"])
                             else r["mdhr_analytic"], ev_grp)
    if extra:
        r.update(extra)
    return r


# ---------------------------------------------------------------------------------------
# STEP 2  the crossed cell models
# ---------------------------------------------------------------------------------------
say("\n" + "=" * 100)
say("STEP 2-5  MODELS")
say("=" * 100)

RES = []
SETS = {"primary": lambda f: f,
        "pre_study": lambda f: f[f.pre_study_only == 1]}


def run_cells(design, frame_all, cells, setname):
    ref = REF_CELL
    others = [c for c in cells if c != ref]
    cmap = {c: f"cl{j}" for j, c in enumerate(others)}
    F = frame_all.copy()
    for c, col in cmap.items():
        F[col] = (F.cell == c).astype(int)
    for key in OUTCOMES:
        f = at_risk(F, key)
        d = pd.DataFrame({"T": f[f"{key}_years"], "E": f[f"{key}_incident"].astype(int),
                          "site": f.site_id, "cell": f.cell,
                          **{col: f[col] for col in cmap.values()},
                          **{c: f[c] for c in ADJ}})
        nref = int((d.cell == ref).sum())
        eref = int(d.loc[d.cell == ref, "E"].sum())
        fitter, dd, status = fit_cox(d, list(cmap.values()))
        nfit, efit = int(len(dd)), int(dd.E.sum())
        # the reference cell itself, carried so N and events are on the face of the file
        RES.append(row(design, setname, key, f"{ref} [REFERENCE]", ref, nref, nref, eref, eref,
                       nfit, efit, status="reference"))
        for c in others:
            ng = int((d.cell == c).sum())
            eg = int(d.loc[d.cell == c, "E"].sum())
            RES.append(row(design, setname, key, c, ref, ng, nref, eg, eref, nfit, efit,
                           fitter=fitter if status == "ok" else None,
                           term=cmap[c], status=status))


for setname, sel in SETS.items():
    run_cells("cells8", sel(X8), CELLS8, setname)
    say(f"cells8   set={setname:10s} N={len(sel(X8)):,}   done  ({time.time()-t_start:.0f}s)")
run_cells("cells12", hab, CELLS12, "primary")
say(f"cells12  set=primary    N={len(hab):,}   done  ({time.time()-t_start:.0f}s)")


# ---------------------------------------------------------------------------------------
# STEP 3  WITHIN each habitual band, low versus normal oxygenation
# ---------------------------------------------------------------------------------------
def run_o2_within_hab(frame_all, setname):
    for band in HAB_ORDER:
        F = frame_all[(frame_all.hab_cat == band)
                      & frame_all.o2_cat.isin([O2_NORMAL, O2_LOW])].copy()
        F["low_o2"] = (F.o2_cat == O2_LOW).astype(int)
        for key in OUTCOMES:
            f = at_risk(F, key)
            d = pd.DataFrame({"T": f[f"{key}_years"], "E": f[f"{key}_incident"].astype(int),
                              "site": f.site_id, "low_o2": f.low_o2,
                              **{c: f[c] for c in ADJ}})
            ng = int((d.low_o2 == 1).sum())
            nr = int((d.low_o2 == 0).sum())
            eg = int(d.loc[d.low_o2 == 1, "E"].sum())
            er = int(d.loc[d.low_o2 == 0, "E"].sum())
            fitter, dd, status = fit_cox(d, ["low_o2"])
            RES.append(row("o2_within_hab", setname, key,
                           f"low T90>10% vs normal T90<=1%  within {band}",
                           f"normal T90<=1% within {band}", ng, nr, eg, er,
                           int(len(dd)), int(dd.E.sum()),
                           fitter=fitter if status == "ok" else None,
                           term="low_o2", status=status,
                           extra={"hab_band": band, "o2_band": "low vs normal"}))


for setname, sel in SETS.items():
    run_o2_within_hab(sel(hab), setname)
    say(f"o2_within_hab  set={setname:10s} done  ({time.time()-t_start:.0f}s)")


# ---------------------------------------------------------------------------------------
# STEP 4  WITHIN each oxygenation band, across habitual bands
# ---------------------------------------------------------------------------------------
O2_GROUPS = {O2_NORMAL: [O2_NORMAL], O2_INTER: [O2_INTER], O2_LOW: [O2_LOW],
             "all oxygen groups": O2_ORDER}


def run_hab_within_o2(frame_all, setname):
    others = [h for h in HAB_ORDER if h != HAB_REF]
    hmap = {h: f"hb{j}" for j, h in enumerate(others)}
    for oname, olist in O2_GROUPS.items():
        F = frame_all[frame_all.o2_cat.isin(olist)].copy()
        for h, col in hmap.items():
            F[col] = (F.hab_cat == h).astype(int)
        for key in OUTCOMES:
            f = at_risk(F, key)
            d = pd.DataFrame({"T": f[f"{key}_years"], "E": f[f"{key}_incident"].astype(int),
                              "site": f.site_id, "hab": f.hab_cat,
                              **{col: f[col] for col in hmap.values()},
                              **{c: f[c] for c in ADJ}})
            nr = int((d.hab == HAB_REF).sum())
            er = int(d.loc[d.hab == HAB_REF, "E"].sum())
            fitter, dd, status = fit_cox(d, list(hmap.values()))
            nfit, efit = int(len(dd)), int(dd.E.sum())
            RES.append(row("hab_within_o2", setname, key,
                           f"{HAB_REF} [REFERENCE]  within {oname}",
                           f"{HAB_REF} within {oname}", nr, nr, er, er, nfit, efit,
                           status="reference", extra={"hab_band": HAB_REF, "o2_band": oname}))
            for h in others:
                ng = int((d.hab == h).sum())
                eg = int(d.loc[d.hab == h, "E"].sum())
                RES.append(row("hab_within_o2", setname, key,
                               f"{h} vs {HAB_REF}  within {oname}",
                               f"{HAB_REF} within {oname}", ng, nr, eg, er, nfit, efit,
                               fitter=fitter if status == "ok" else None,
                               term=hmap[h], status=status,
                               extra={"hab_band": h, "o2_band": oname}))


for setname, sel in SETS.items():
    run_hab_within_o2(sel(hab), setname)
    say(f"hab_within_o2  set={setname:10s} done  ({time.time()-t_start:.0f}s)")


# ---------------------------------------------------------------------------------------
# STEP 5  formal interaction, likelihood ratio tests
# ---------------------------------------------------------------------------------------
def lrt(full, reduced, df):
    stat = 2 * (full.log_likelihood_ - reduced.log_likelihood_)
    stat = max(stat, 0.0)
    return float(stat), float(stats.chi2.sf(stat, df))


def run_interaction(frame_all, setname):
    others_h = [h for h in HAB_ORDER if h != HAB_REF]
    hmap = {h: f"hb{j}" for j, h in enumerate(others_h)}
    # categorical, 3 df
    F = frame_all[frame_all.o2_cat.isin([O2_NORMAL, O2_LOW])].copy()
    F["low_o2"] = (F.o2_cat == O2_LOW).astype(int)
    for h, col in hmap.items():
        F[col] = (F.hab_cat == h).astype(int)
    for h, col in hmap.items():
        F[f"{col}_x"] = F[col] * F.low_o2
    main = list(hmap.values()) + ["low_o2"]
    inter = [f"{c}_x" for c in hmap.values()]
    for key in OUTCOMES:
        f = at_risk(F, key)
        d = pd.DataFrame({"T": f[f"{key}_years"], "E": f[f"{key}_incident"].astype(int),
                          "site": f.site_id,
                          **{c: f[c] for c in main + inter}, **{c: f[c] for c in ADJ}})
        fr, dd, st_f = fit_cox(d, main + inter)
        rd, _, st_r = fit_cox(d, main)
        rec = row("interaction", setname, key,
                  "habitual band x oxygenation, 3 df categorical LRT", REF_CELL,
                  int(len(dd)), int(len(dd)), int(dd.E.sum()), int(dd.E.sum()),
                  int(len(dd)), int(dd.E.sum()),
                  status=st_f if st_f != "ok" else st_r)
        if st_f == "ok" and st_r == "ok":
            chi2, p = lrt(fr, rd, 3)
            rec.update({"chi2": chi2, "df": 3, "p": p, "status": "ok"})
        rec["mdhr_analytic"] = np.nan
        rec["power"] = "interaction tests are underpowered here, read the strata"
        RES.append(rec)

    # continuous, 1 df, more power. Uses all three oxygen groups because a continuous
    # exposure needs no cut point, so the intermediate band contributes.
    G = frame_all.copy()
    G["hxo"] = G.z_habitual_hours * G.z_lab_t90
    for key in OUTCOMES:
        f = at_risk(G, key)
        d = pd.DataFrame({"T": f[f"{key}_years"], "E": f[f"{key}_incident"].astype(int),
                          "site": f.site_id,
                          "z_habitual_hours": f.z_habitual_hours, "z_lab_t90": f.z_lab_t90,
                          "hxo": f.hxo, **{c: f[c] for c in ADJ}})
        fr, dd, st_f = fit_cox(d, ["z_habitual_hours", "z_lab_t90", "hxo"])
        rd, _, st_r = fit_cox(d, ["z_habitual_hours", "z_lab_t90"])
        rec = row("interaction", setname, key,
                  "habitual hours x T90, 1 df continuous LRT (per SD, RINT within site)",
                  "no interaction", int(len(dd)), int(len(dd)),
                  int(dd.E.sum()), int(dd.E.sum()), int(len(dd)), int(dd.E.sum()),
                  fitter=fr if st_f == "ok" else None, term="hxo",
                  status=st_f if st_f != "ok" else st_r)
        if st_f == "ok" and st_r == "ok":
            chi2, p = lrt(fr, rd, 1)
            rec.update({"chi2": chi2, "df": 1, "p": p, "status": "ok"})
        rec["mdhr_analytic"] = np.nan
        rec["power"] = "interaction tests are underpowered here, read the strata"
        RES.append(rec)


for setname, sel in SETS.items():
    run_interaction(sel(hab), setname)
    say(f"interaction    set={setname:10s} done  ({time.time()-t_start:.0f}s)")


# ---------------------------------------------------------------------------------------
# STEP 5b  the head-to-head continuous comparison
# ---------------------------------------------------------------------------------------
CONT_SETS = {"habitual_8711": hab, "pre_study_4631": hab[hab.pre_study_only == 1],
             "full_cohort_19173": b}
for sname, frame in CONT_SETS.items():
    for key in OUTCOMES:
        f = at_risk(frame, key)
        for ename in CONT:
            col = f"z_{ename}"
            if col not in f.columns or f[col].isna().all():
                continue
            d = pd.DataFrame({"T": f[f"{key}_years"], "E": f[f"{key}_incident"].astype(int),
                              "site": f.site_id, col: f[col], **{c: f[c] for c in ADJ}})
            fitter, dd, status = fit_cox(d, [col])
            RES.append(row("continuous", sname, key, f"{ename} per SD",
                           "per 1 SD, rank-based inverse normal within site",
                           int(len(dd)), int(len(dd)), int(dd.E.sum()), int(dd.E.sum()),
                           int(len(dd)), int(dd.E.sum()),
                           fitter=fitter if status == "ok" else None, term=col, status=status,
                           extra={"per_sd_floor": PER_SD_FLOOR,
                                  "per_sd_floor_legacy": PER_SD_FLOOR_LEGACY}))
    say(f"continuous     set={sname:18s} N={len(frame):,}  done  ({time.time()-t_start:.0f}s)")

res = pd.DataFrame(RES)
say(f"\n{len(res):,} estimate rows")


# ---------------------------------------------------------------------------------------
# STEP 6  negative controls and the floor
# ---------------------------------------------------------------------------------------
say("\n" + "=" * 100)
say("STEP 6  NEGATIVE CONTROLS AND THE FLOOR")
say("=" * 100)
say(f"controls carried, all 5, imported from numbers/disease_definitions.py: "
    f"{', '.join(NEG_CONTROLS)}")
say(f"NOT controls: {', '.join(NOT_A_CONTROL)}")
say(f"  osteoarthritis and fracture were controls in the brief for this run and were removed")
say(f"  upstream on 2026-08-07. Both are still ordinary outcomes here.")
say(f"per-SD floor {PER_SD_FLOOR:.3f}, set by {PER_SD_FLOOR_SETBY}, from "
    f"numbers/negcontrols_v4_panel.csv")
say(f"the brief's legacy per-SD floor was {PER_SD_FLOOR_LEGACY}, carried alongside")


def dev(h):
    """Distance from the null on a symmetric scale: max(HR, 1/HR)."""
    return np.where(h > 0, np.maximum(h, 1.0 / h), np.nan)


res["deviation"] = dev(res.hr.values)
res["ci_excludes_1"] = (res.lo95 > 1) | (res.hi95 < 1)
ok = res.status.eq("ok") & res.hr.notna()

# Two floors are computed for every contrast, because with only four controls the symmetric
# floor is hostage to a single unstable estimate. Within <5h, glaucoma comes back at 0.405 on
# 11 events, which alone sets a symmetric floor of 2.47 and would erase real findings.
#   nc_floor            symmetric, max(HR, 1/HR) over the controls. Conservative, noisy.
#   nc_floor_upper      the largest control HR. The floor a HARMFUL association must clear.
#   nc_floor_lower      the smallest control HR. The floor a PROTECTIVE association must clear.
# The directional pair is the one to read; the symmetric one is carried so the conservative
# reading is on the face of the file too.
floor_rows = []
for (model, setname, contrast), g in res[ok & res.negative_control].groupby(
        ["model", "set", "contrast"]):
    # v8.1f (2026-09-15): a control sets a cell floor only when it meets the same event minimum as the outcomes it
    # guards (MIN_EVENTS); with the split nights' sleep time masked, inguinal hernia had 14 events in one cell and a
    # hazard ratio of 10.1 (1.3 to 77), which set that cell's floor. With no qualifying control the contrast falls
    # back to the paper's per-SD floor (the merge default below). The earlier rule was 10 events, any control.
    usable = g[g.events_group >= MIN_EVENTS]
    if not len(usable):
        continue
    src = usable
    floor_rows.append({"model": model, "set": setname, "contrast": contrast,
                       "n_controls": int(len(src)),
                       "controls_used": ", ".join(sorted(src.outcome_key)),
                       "floor": float(src.deviation.max()),
                       "floor_upper": float(src.hr.max()),
                       "floor_lower": float(src.hr.min()),
                       "floor_median_deviation": float(src.deviation.median()),
                       "n_controls_ci_off_1": int(src.ci_excludes_1.sum()),
                       "driving_control": str(src.loc[src.deviation.idxmax(), "outcome_key"])})
floors = pd.DataFrame(floor_rows)
floors.to_csv(f"{HERE}/floors_crossed.csv", index=False)

fmap = floors.set_index(["model", "set", "contrast"])
idx = list(zip(res.model, res.set, res.contrast))
for src, dst, dflt in (("floor", "nc_floor", PER_SD_FLOOR),
                       ("floor_upper", "nc_floor_upper", PER_SD_FLOOR),
                       ("floor_lower", "nc_floor_lower", 1 / PER_SD_FLOOR)):
    m = fmap[src].to_dict()
    res[dst] = [m.get(i, np.nan) for i in idx]
    res.loc[res.model.eq("continuous"), dst] = res.loc[res.model.eq("continuous"),
                                                       dst].fillna(dflt)
# symmetric, conservative
res["exceeds_nc_floor"] = res.ci_excludes_1 & (res.deviation > res.nc_floor)
# directional, the one to read
res["exceeds_floor_directional"] = res.ci_excludes_1 & (
    ((res.hr > 1) & (res.hr > res.nc_floor_upper)) |
    ((res.hr < 1) & (res.hr < res.nc_floor_lower)))
for c in ("exceeds_nc_floor", "exceeds_floor_directional"):
    res.loc[res.negative_control, c] = pd.NA

say("\ncontrol floors on the eight-cell contrasts, primary set")
say(f"  {'contrast':42s}{'sym':>7}{'upper':>8}{'lower':>8}{'medDev':>8}  driven by")
for _, r in floors[(floors.model == "cells8") & (floors.set == "primary")].iterrows():
    say(f"  {r.contrast:42s}{r.floor:7.2f}{r.floor_upper:8.2f}{r.floor_lower:8.2f}"
        f"{r.floor_median_deviation:8.2f}  {r.driving_control}"
        f"  ({int(r.n_controls_ci_off_1)} of {int(r.n_controls)} controls off 1)")
say("\ncontrol floors on the within-band oxygen contrasts, primary set")
say(f"  {'contrast':52s}{'sym':>7}{'upper':>8}{'lower':>8}  driven by")
for _, r in floors[(floors.model == "o2_within_hab") & (floors.set == "primary")].iterrows():
    say(f"  {r.contrast:52s}{r.floor:7.2f}{r.floor_upper:8.2f}{r.floor_lower:8.2f}"
        f"  {r.driving_control}"
        f"  ({int(r.n_controls_ci_off_1)} of {int(r.n_controls)} controls off 1)")
say("\nA control whose own confidence interval excludes 1 is a warning about that contrast, not")
say("just a floor. Read the driving_control column in floors_crossed.csv before quoting a floor.")

nc_cont = res[ok & res.negative_control & res.model.eq("continuous")
              & res.set.eq("habitual_8711")]
say("\nnegative controls on the per-SD scale in the habitual set")
for e in CONT:
    s = nc_cont[nc_cont.contrast == f"{e} per SD"]
    if len(s):
        say(f"  {e:16s} controls span {s.hr.min():.3f} to {s.hr.max():.3f}"
            f"   (standing floor {PER_SD_FLOOR})")


# ---------------------------------------------------------------------------------------
# write
# ---------------------------------------------------------------------------------------
ORDER = ["model", "set", "outcome_key", "outcome", "organ_group", "negative_control", "circular",
         "hab_band", "o2_band", "contrast", "reference",
         "n_group", "n_reference", "events_group", "events_reference", "n_model", "events_model",
         "hr", "lo95", "hi95", "se_loghr", "p", "chi2", "df",
         "ci_excludes_1", "deviation", "nc_floor", "nc_floor_upper", "nc_floor_lower",
         "exceeds_nc_floor", "exceeds_floor_directional", "per_sd_floor",
         "per_sd_floor_legacy",
         "mdhr_analytic", "mdhr_empirical", "power", "epv", "status"]
for c in ORDER:
    if c not in res.columns:
        res[c] = np.nan
res = res[ORDER]
for c in ["hr", "lo95", "hi95", "se_loghr", "deviation", "nc_floor", "nc_floor_upper",
          "nc_floor_lower", "mdhr_analytic", "mdhr_empirical"]:
    res[c] = res[c].astype(float).round(4)
res["p"] = res.p.astype(float)
res.to_csv(f"{HERE}/results_crossed.csv", index=False)
say(f"\nwrote results_crossed.csv  {len(res):,} rows x {res.shape[1]} cols")


# ---------------------------------------------------------------------------------------
# readout
# ---------------------------------------------------------------------------------------
low_cells = [c for c in CELLS8 if c.endswith(O2_LOW)]
norm_cells = [c for c in CELLS8 if c.endswith(O2_NORMAL) and c != REF_CELL]


def show(model, setname, keys, order=None, width=50):
    sub = res[(res.model == model) & (res.set == setname) & res.outcome_key.isin(keys)]
    if order:
        sub = sub[sub.contrast.isin(order)]
    for key in keys:
        g = sub[sub.outcome_key == key]
        if not len(g):
            continue
        say(f"\n  {LABEL[key]}  ({int(g.events_model.max()):,} events, "
            f"N={int(g.n_model.max()):,})")
        it = g.set_index("contrast").reindex(order).reset_index() if order else g
        for _, r in it.iterrows():
            if pd.isna(r.get("contrast")):
                continue
            if r.status == "reference":
                say(f"    {r.contrast:{width}s}  N={int(r.n_group):5,} ev={int(r.events_group):4,}"
                    f"   1.00 (reference)")
            elif np.isfinite(r.hr):
                star = "*" if r.ci_excludes_1 else " "
                say(f"    {r.contrast:{width}s}  N={int(r.n_group):5,} ev={int(r.events_group):4,}"
                    f"   {r.hr:5.2f} [{r.lo95:4.2f}-{r.hi95:5.2f}]{star} p={r.p:.2e}"
                    f"   MDHR {r.mdhr_empirical:4.2f}  {r.power}")
            else:
                say(f"    {r.contrast:{width}s}  N={int(r.n_group):5,} ev={int(r.events_group):4,}"
                    f"   {r.status}")


say("\n" + "=" * 100)
say("RESULT 1  THE FOUR CELLS ALEN NAMED, AND THE REST, AGAINST 6-<7h WITH NORMAL OXYGEN")
say("=" * 100)
show("cells8", "primary", HEADLINE, order=[f"{REF_CELL} [REFERENCE]"] +
     [c for c in CELLS8 if c != REF_CELL])

say("\n" + "=" * 100)
say("RESULT 2  DOES OXYGEN HOLD UP INSIDE EVERY HABITUAL SLEEP BAND")
say("=" * 100)
show("o2_within_hab", "primary", HEADLINE,
     order=[f"low T90>10% vs normal T90<=1%  within {h}" for h in HAB_ORDER], width=52)

say("\n" + "=" * 100)
say("RESULT 3  THE MIRROR, HABITUAL SLEEP INSIDE EVERY OXYGEN BAND")
say("=" * 100)
for oname in O2_GROUPS:
    say(f"\n--- within {oname} ---")
    show("hab_within_o2", "primary", HEADLINE,
         order=[f"{HAB_REF} [REFERENCE]  within {oname}"] +
               [f"{h} vs {HAB_REF}  within {oname}" for h in HAB_ORDER if h != HAB_REF],
         width=46)

say("\n" + "=" * 100)
say("RESULT 4  INTERACTION TESTS")
say("=" * 100)
it = res[(res.model == "interaction") & (res.set == "primary")]
for tag, dfree in (("3 df categorical", 3), ("1 df continuous", 1)):
    s = it[it.df == dfree]
    say(f"\n{tag}: {int(s.p.notna().sum())} outcomes tested, "
        f"{int((s.p < 0.05).sum())} with p<0.05, "
        f"{int((s.p < 0.05 / max(int(s.p.notna().sum()), 1)).sum())} surviving Bonferroni")
    for _, r in s.sort_values("p").head(6).iterrows():
        say(f"    {r.outcome:28s} chi2={r.chi2:6.2f} df={int(r.df)} p={r.p:.4f}")

say("\n" + "=" * 100)
say("RESULT 5  HEAD TO HEAD, PER SD, IN THE SAME PATIENTS")
say("=" * 100)
cc = res[(res.model == "continuous") & (res.set == "habitual_8711") & res.status.eq("ok")]
for e in CONT:
    s = cc[cc.contrast == f"{e} per SD"]
    sig = s[s.ci_excludes_1 & ~s.circular]
    above = s[(s.deviation > PER_SD_FLOOR) & s.ci_excludes_1 & ~s.circular & ~s.negative_control]
    say(f"\n{e:16s} {len(s)} outcomes, {len(sig)} with a confidence interval excluding 1, "
        f"{len(above)} also above the {PER_SD_FLOOR} floor")
    for _, r in s.reindex(s.deviation.sort_values(ascending=False).index).head(6).iterrows():
        say(f"    {r.outcome:28s} {r.hr:5.3f} [{r.lo95:.3f}-{r.hi95:.3f}] p={r.p:.2e}"
            f"  ev={int(r.events_model):,}")

say("\n" + "=" * 100)
say("RESULT 6  POWER, HOW MANY OF THESE NULLS MEAN ANYTHING")
say("=" * 100)
for model, setname in (("cells8", "primary"), ("o2_within_hab", "primary"),
                       ("cells8", "pre_study"), ("o2_within_hab", "pre_study")):
    s = res[(res.model == model) & (res.set == setname) & res.status.eq("ok")]
    if not len(s):
        continue
    null = s[~s.ci_excludes_1]
    say(f"\n{model} / {setname}: {len(s):,} fitted contrasts, {len(null):,} not significant")
    say(f"  of those nulls, minimum detectable HR is under 1.35 in {int((null.mdhr_empirical<1.35).sum()):,}, "
        f"1.35-1.75 in {int(((null.mdhr_empirical>=1.35)&(null.mdhr_empirical<1.75)).sum()):,}, "
        f"1.75-2.50 in {int(((null.mdhr_empirical>=1.75)&(null.mdhr_empirical<2.50)).sum()):,}, "
        f"2.50 or worse in {int((null.mdhr_empirical>=2.50).sum()):,}")
    say(f"  median MDHR across all fitted contrasts {s.mdhr_empirical.median():.2f}, "
        f"worst {s.mdhr_empirical.max():.2f}")

say("\n" + "=" * 100)
say(f"RESULT 7  THE WHOLE PANEL, ALL {len(OUTCOMES)} OUTCOMES, NOT JUST THE SIX HEADLINES")
say("=" * 100)


def panel(model, setname, keycol, keys, label):
    say(f"\n{label}   (set={setname})")
    say(f"  {'stratum':26s}{'fitted':>8}{'sig':>6}{'+dirFloor':>11}{'+symFloor':>11}"
        f"{'medHR':>8}{'medMDHR':>9}{'notfit':>8}")
    for k in keys:
        g = res[(res.model == model) & (res.set == setname) & (res[keycol] == k)
                & ~res.negative_control & ~res.circular & ~res.contrast.str.contains("REFERENCE")]
        f = g[g.status.eq("ok") & g.hr.notna()]
        if not len(f):
            continue
        say(f"  {k:26s}{len(f):8,}{int(f.ci_excludes_1.sum()):6,}"
            f"{int((f.exceeds_floor_directional == True).sum()):11,}"
            f"{int((f.exceeds_nc_floor == True).sum()):11,}{f.hr.median():8.2f}"
            f"{f.mdhr_empirical.median():9.2f}{len(g)-len(f):8,}")


panel("o2_within_hab", "primary", "hab_band", HAB_ORDER,
      "OXYGEN (low>10% vs normal<=1%) INSIDE EACH HABITUAL BAND")
panel("o2_within_hab", "pre_study", "hab_band", HAB_ORDER,
      "OXYGEN INSIDE EACH HABITUAL BAND, pre-study exposures only")
panel("hab_within_o2", "primary", "o2_band", list(O2_GROUPS),
      "HABITUAL SLEEP (each band vs 6-<7h) INSIDE EACH OXYGEN BAND")
panel("hab_within_o2", "pre_study", "o2_band", list(O2_GROUPS),
      "HABITUAL SLEEP INSIDE EACH OXYGEN BAND, pre-study exposures only")

say(f"\nEIGHT-CELL MODEL, low-oxygen cells against normal-oxygen cells, all {len(OUTCOMES)} outcomes")
for setname in ("primary", "pre_study"):
    s = res[(res.model == "cells8") & (res.set == setname) & res.status.eq("ok")
            & ~res.negative_control & ~res.circular & res.hr.notna()]
    lo = s[s.contrast.isin(low_cells)]
    no = s[s.contrast.isin(norm_cells)]
    say(f"  {setname:10s} low-oxygen cells  {int(lo.ci_excludes_1.sum()):3,} of {len(lo):3,} "
        f"off 1, median HR {lo.hr.median():.2f}, "
        f"{int((lo.exceeds_floor_directional == True).sum()):3,} above the directional floor, "
        f"{int((lo.exceeds_nc_floor == True).sum()):3,} above the symmetric floor")
    say(f"  {setname:10s} normal-oxygen     {int(no.ci_excludes_1.sum()):3,} of {len(no):3,} "
        f"off 1, median HR {no.hr.median():.2f}, "
        f"{int((no.exceeds_floor_directional == True).sum()):3,} above the directional floor, "
        f"{int((no.exceeds_nc_floor == True).sum()):3,} above the symmetric floor")

say(f"\nPER-SD HEAD TO HEAD, outcomes with a CI off 1 that also clear the confounding floor")
say(f"(current floor {PER_SD_FLOOR:.3f} set by {PER_SD_FLOOR_SETBY}; "
    f"the brief's legacy floor {PER_SD_FLOOR_LEGACY} in brackets)")
say(f"  {'exposure':18s}{'habitual_8711':>16}{'pre_study_4631':>17}{'full_19173':>13}")
for e in CONT:
    line = f"  {e:18s}"
    for sname in ("habitual_8711", "pre_study_4631", "full_cohort_19173"):
        s = res[(res.model == "continuous") & (res.set == sname)
                & (res.contrast == f"{e} per SD") & res.status.eq("ok")
                & ~res.negative_control & ~res.circular]
        n = int((s.ci_excludes_1 & (s.deviation > PER_SD_FLOOR)).sum()) if len(s) else 0
        nl = int((s.ci_excludes_1 & (s.deviation > PER_SD_FLOOR_LEGACY)).sum()) if len(s) else 0
        nd = int((s.exceeds_floor_directional == True).sum()) if len(s) else 0
        cell = f"{nd} / {n} [{nl}]"
        line += f"{cell:>16s}" if sname != "full_cohort_19173" else f"{cell:>13s}"
    say(line)
say("read as: clearing this set's own control band / clearing 1.063 [clearing 1.133]")
say("\nthe control band each exposure actually sits in, habitual_8711")
for e in CONT:
    g = res[(res.model == "continuous") & (res.set == "habitual_8711")
            & (res.contrast == f"{e} per SD") & res.negative_control & res.status.eq("ok")]
    if len(g):
        say(f"  {e:16s} all 5 controls between {g.hr.min():.3f} and {g.hr.max():.3f}. "
            f"An association must fall outside that band to be worth anything.")

say("\nwhere the non-identified contrasts sit")
_ni = res[res.status.eq("not identified (se>2)")]
for (m_, s_), g_ in _ni.groupby(["model", "set"]):
    say(f"  {m_} / {s_}: {len(g_)} rows, median events in the group "
        f"{g_.events_group.median():.0f}")
say("\nfit status across every row")
for st, n in res.status.value_counts().items():
    say(f"  {st:40s}{n:6,}")

# per cell median MDHR, the plain answer to "a null in a cell of 40 people means nothing"
say("\nminimum detectable HR by cell, median over the 52 outcomes, primary eight-cell model")
s = res[(res.model == "cells8") & (res.set == "primary") & res.status.eq("ok")]
for c in CELLS8:
    if c == REF_CELL:
        continue
    g = s[s.contrast == c]
    if not len(g):
        continue
    say(f"  {c:40s} N={int(g.n_group.max()):5,}  median MDHR {g.mdhr_empirical.median():5.2f}"
        f"   best {g.mdhr_empirical.min():4.2f}  worst {g.mdhr_empirical.max():5.2f}")

# ---------------------------------------------------------------------------------------
# summary
# ---------------------------------------------------------------------------------------
c8 = res[(res.model == "cells8") & (res.set == "primary") & res.status.eq("ok")
         & ~res.negative_control & ~res.circular]
o2w = res[(res.model == "o2_within_hab") & (res.set == "primary") & res.status.eq("ok")
          & ~res.negative_control & ~res.circular]
habw = res[(res.model == "hab_within_o2") & (res.set == "primary") & res.status.eq("ok")
           & ~res.negative_control & ~res.circular]
o2p = res[(res.model == "o2_within_hab") & (res.set == "pre_study") & res.status.eq("ok")
          & ~res.negative_control & ~res.circular]
habp = res[(res.model == "hab_within_o2") & (res.set == "pre_study") & res.status.eq("ok")
           & ~res.negative_control & ~res.circular]

summary = {
    "built": time.strftime("%Y-%m-%d %H:%M"),
    "runtime_sec": round(time.time() - t_start, 1),
    "cohort_n": int(len(b)),
    "habitual_n": int(len(hab)),
    "crossed_8cell_n": int(len(X8)),
    "pre_study_n": int(hab.pre_study_only.sum()),
    "pre_study_crossed_8cell_n": int(len(X8[X8.pre_study_only == 1])),
    "reference_cell": REF_CELL,
    "reference_cell_n": int((X8.cell == REF_CELL).sum()),
    "outcomes": len(OUTCOMES),
    "outcomes_ge_20_events_in_crossed_set": int(sum(
        at_risk(X8, k)[f"{k}_incident"].sum() >= MIN_EVENTS for k in OUTCOMES)),
    "cells_with_outcomes_under_5_events": {c: n for c, n in small},
    "smallest_cell": {"cell": cells_df[cells_df.design == "cells8"].sort_values("n").iloc[0].cell,
                      "n": int(cells_df[cells_df.design == "cells8"].n.min())},
    "negative_controls_carried": NEG_CONTROLS,
    "negative_controls_merged_in_from_t90_final": NEW_CONTROLS,
    "not_controls": NOT_A_CONTROL,
    "control_panel_changed_under_this_run": (
        "the brief named 6 controls incl. osteoarthritis with a per-SD floor of 1.133. "
        "numbers/disease_definitions.py was rewritten 2026-08-07 to 5 controls without "
        "osteoarthritis, and the two replacements were materialised into t90_final.parquet "
        "at 16:41 the same day. This run uses the current upstream panel."),
    "per_sd_floor": PER_SD_FLOOR,
    "per_sd_floor_set_by": PER_SD_FLOOR_SETBY,
    "per_sd_floor_legacy_from_brief": PER_SD_FLOOR_LEGACY,
    "low_o2_cells_significant": int(c8[c8.contrast.isin(low_cells) & c8.ci_excludes_1].shape[0]),
    "low_o2_cells_fitted": int(c8[c8.contrast.isin(low_cells)].shape[0]),
    "normal_o2_cells_significant": int(c8[c8.contrast.isin(norm_cells)
                                          & c8.ci_excludes_1].shape[0]),
    "normal_o2_cells_fitted": int(c8[c8.contrast.isin(norm_cells)].shape[0]),
    "median_hr_low_o2_cells": float(c8[c8.contrast.isin(low_cells)].hr.median()),
    "median_hr_normal_o2_cells": float(c8[c8.contrast.isin(norm_cells)].hr.median()),
    "low_o2_cells_above_directional_floor": int(
        (c8[c8.contrast.isin(low_cells)].exceeds_floor_directional == True).sum()),
    "normal_o2_cells_above_directional_floor": int(
        (c8[c8.contrast.isin(norm_cells)].exceeds_floor_directional == True).sum()),
    "low_o2_cells_above_symmetric_floor": int(
        (c8[c8.contrast.isin(low_cells)].exceeds_nc_floor == True).sum()),
    "normal_o2_cells_above_symmetric_floor": int(
        (c8[c8.contrast.isin(norm_cells)].exceeds_nc_floor == True).sum()),
    "o2_within_band": {b_: {
        "outcomes_fitted": int(o2w[o2w.hab_band == b_].shape[0]),
        "significant": int(o2w[(o2w.hab_band == b_) & o2w.ci_excludes_1].shape[0]),
        "significant_and_above_directional_floor": int(
            (o2w[o2w.hab_band == b_].exceeds_floor_directional == True).sum()),
        "significant_and_above_symmetric_floor": int(
            (o2w[o2w.hab_band == b_].exceeds_nc_floor == True).sum()),
        "median_hr": float(o2w[o2w.hab_band == b_].hr.median()),
        "median_mdhr": float(o2w[o2w.hab_band == b_].mdhr_empirical.median()),
    } for b_ in HAB_ORDER},
    "o2_within_band_pre_study": {b_: {
        "outcomes_fitted": int(o2p[o2p.hab_band == b_].shape[0]),
        "significant": int(o2p[(o2p.hab_band == b_) & o2p.ci_excludes_1].shape[0]),
        "median_hr": float(o2p[o2p.hab_band == b_].hr.median()),
        "median_mdhr": float(o2p[o2p.hab_band == b_].mdhr_empirical.median()),
    } for b_ in HAB_ORDER},
    "hab_within_o2": {o_: {
        "outcomes_fitted": int(habw[habw.o2_band == o_].shape[0]),
        "significant": int(habw[(habw.o2_band == o_) & habw.ci_excludes_1].shape[0]),
        "significant_and_above_directional_floor": int(
            (habw[habw.o2_band == o_].exceeds_floor_directional == True).sum()),
        "significant_and_above_symmetric_floor": int(
            (habw[habw.o2_band == o_].exceeds_nc_floor == True).sum()),
        "median_hr": float(habw[habw.o2_band == o_].hr.median()),
        "median_mdhr": float(habw[habw.o2_band == o_].mdhr_empirical.median()),
    } for o_ in O2_GROUPS},
    "hab_within_o2_pre_study": {o_: {
        "outcomes_fitted": int(habp[habp.o2_band == o_].shape[0]),
        "significant": int(habp[(habp.o2_band == o_) & habp.ci_excludes_1].shape[0]),
        "median_hr": float(habp[habp.o2_band == o_].hr.median()),
        "median_mdhr": float(habp[habp.o2_band == o_].mdhr_empirical.median()),
    } for o_ in O2_GROUPS},
    "status_counts": {str(k): int(v) for k, v in res.status.value_counts().items()},
    "mdhr_median_by_cell": {c: float(res[(res.model == "cells8") & (res.set == "primary")
                                         & res.status.eq("ok") & (res.contrast == c)]
                                     .mdhr_empirical.median())
                            for c in CELLS8 if c != REF_CELL},
    "interaction": {
        "categorical_3df_p_lt_05": int(((it.df == 3) & (it.p < 0.05)).sum()),
        "categorical_3df_tested": int(((it.df == 3) & it.p.notna()).sum()),
        "continuous_1df_p_lt_05": int(((it.df == 1) & (it.p < 0.05)).sum()),
        "continuous_1df_tested": int(((it.df == 1) & it.p.notna()).sum()),
    },
    "continuous_head_to_head": {e: {
        "n_significant": int(cc[(cc.contrast == f"{e} per SD") & cc.ci_excludes_1
                                & ~cc.circular].shape[0]),
        "n_above_floor": int(cc[(cc.contrast == f"{e} per SD") & cc.ci_excludes_1
                                & ~cc.circular & ~cc.negative_control
                                & (cc.deviation > PER_SD_FLOOR)].shape[0]),
        "n_above_legacy_floor_1133": int(cc[(cc.contrast == f"{e} per SD") & cc.ci_excludes_1
                                            & ~cc.circular & ~cc.negative_control
                                            & (cc.deviation > PER_SD_FLOOR_LEGACY)].shape[0]),
        "n_fitted": int(cc[cc.contrast == f"{e} per SD"].shape[0]),
        "max_deviation": float(cc[cc.contrast == f"{e} per SD"].deviation.max()),
    } for e in CONT},
    "files": ["results_crossed.csv", "cells_crossed.csv", "floors_crossed.csv",
              "summary_crossed.json", "model_crossed_log.txt"],
}
json.dump(summary, open(f"{HERE}/summary_crossed.json", "w"), indent=1, default=str)
say(f"\nwrote summary_crossed.json")
say(f"\nDONE in {time.time()-t_start:.0f}s")
LOG.close()
