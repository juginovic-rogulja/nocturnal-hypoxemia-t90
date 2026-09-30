"""
Habitual home sleep duration against incident disease, and the head-to-head against the laboratory.

The question. Laboratory total sleep time ranks 190th of 197 sleep-study measures for predicting
what happens to a patient afterwards (held-out concordance gain -0.001016 over age and sex), and
inside these same patients the laboratory night and the patient's habitual home night are
unrelated (r = -0.024). The obvious objection to that result is that one night in a laboratory
simply mismeasured the construct. This script tests the objection directly. If habitual home
sleep, a measure of the construct that the laboratory cannot corrupt, also fails to predict
disease, then the sleep-duration result is about sleep duration and not about the measurement.

What is fitted.
  Cox proportional hazards, UNPENALIZED, stratified by hospital, adjusted for a natural cubic
  age spline and sex. The age basis is cr(age, df=4) with ONE COLUMN DROPPED, which leaves the
  3 columns the analysis plan requires: a 4-column basis is rank-deficient and an unpenalized
  solver silently fits nothing on it.

Exposure scales.
  hab_cat    four habitual bands, reference 6 to under 7 hours
  hab_z      habitual hours, rank-inverse-normal within hospital, so one unit is one SD of the
             normal score. This is the scale the published T90 and lab-TST hazard ratios are on,
             so it is the only scale on which the three can be compared.
  hab_hour   habitual hours, raw, per +1 hour
  hab_sd     habitual hours, raw, per +1 SD of the primary set (2.32 h)
  tst_z      laboratory total sleep time, same within-hospital normal score
  tst_hour   laboratory total sleep time, per +1 hour, so it pairs with hab_hour
  t90_z      laboratory T90, same within-hospital normal score

  The normal scores are computed ONCE on the 8,711-patient primary set and carried unchanged into
  every subset, so a tier A hazard ratio is per SD of the primary-set distribution and not per SD
  of tier A. Recomputing them per subset would silently rescale the comparison.

Direction. For hab and tst a hazard ratio above 1 means LONGER sleep carries more risk. For t90 a
hazard ratio above 1 means WORSE oxygenation carries more risk.

The floor. The negative-control panel settled on 2026-08-07 is back pain, cataract, glaucoma,
contact dermatitis and hemorrhoids, and it is read from disease_definitions.NEGATIVE_CONTROLS
rather than typed here. Fracture and osteoarthritis are carried as outcomes and are NOT
controls: a plausible causal path runs from the exposure to each of them. Thyroid disease is
not a control either. All five are defined in numbers/disease_definitions.py, so the panel is
complete. Two floors are applied. The standing per-SD floor is the largest negative-control
hazard ratio in the published T90 panel, 1.063 set by glaucoma, read from
numbers/negcontrols_final.json. The empirical floor is the largest control effect measured in
this exact model, on this exact scale, and it is the one that governs the categorical models where
a per-SD floor has no meaning.

Concordance. The published -0.001016 is a held-out gain: fit in one hospital, score in the other,
average the two directions, subtract the age-and-sex base model. That protocol is reproduced here
for habitual sleep, laboratory sleep and T90 inside the habitual set. It is a PREDICTION exercise,
not an estimation one, so it keeps the published specification exactly, which means the 4-column
spline and the small ridge the published ranking used. Nothing is read off it as a hazard ratio.

Sets. Primary is all 8,711. Sensitivity in order of increasing measurement quality and decreasing
size: pre-study only (4,631, the set in which the exposure cannot have been recorded after the
outcome clock started), tiers A and B (7,795, drops the bed-time-derived tier C), tier A (1,543),
and the structured template field (1,293, I0002 only so it cannot be hospital-stratified).

Writes results_habitual.csv and five companion files.

v8.1g (15 Sept 2026, step 237, owner decision A8 of round 37). The identical-participant concordance block
(Supplementary Table 19 Panel B) belongs to the concordance family, which the v8.1 rule computes on the 15,551
full diagnostic nights. With T90_FULL_NIGHTS_ONLY=1 in the environment (the psv command of step 237) this script
merges the split-night flag from the frozen t90_final.parquet (analysis.parquet does not carry it), drops the
3,622 split-night recordings before anything is fitted, recomputes the primary set, the within-hospital normal
scores and the continuous set on the full nights, and writes every output with the cohort_tag() suffix
(results_habitual_fullnights.csv, results_habitual_concordance_fullnights.csv, results_habitual_summary_fullnights.json,
model_log_fullnights.txt and the companions), so the all-nights files of step 150 are never overwritten.
Without the flag the script runs exactly as before (step 150). Backup model_habitual_PRE_V8_1g.py.
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
from lifelines.utils import concordance_index
from patsy import dmatrix
from scipy import stats
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")

HERE = f"{paths.SV_ROOT}/habitual_outcomes_v2"
T90ROOT = paths.T90_ROOT
sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
from disease_definitions import DISEASES, NEGATIVE_CONTROLS, CIRCULAR, ORGAN_GROUP  # noqa: E402
from cohort_spec import COHORT_N, NUMBERS_DIR, N_MEASURES, FULL_NIGHTS_ONLY, SPLIT_FLAG, T90_FINAL, cohort_tag  # noqa: E402  v8 sweep 2026-09-12; v8.1g: the full-nights flag, the split flag column, the frozen table path and the output suffix (step 237); F1 moved it up from the last import block (NUMBERS_DIR, COHORT_N are used below); re-applied by lane R1 2026-09-13 after the iCloud rollback

# largest negative-control per-SD hazard ratio in the T90 panel, read never retyped
import json as _json  # noqa: E402
with open(f"{paths.NUMBERS_DIR}/"
          "negcontrols_final.json") as _fh:
    _NCF = _json.load(_fh)
STANDING_FLOOR = float(_NCF["confounding_floor_per_sd"])
STANDING_FLOOR_SETBY = _NCF["confounding_floor_driver"]
# v8 F1 fix (2026-09-12, integrity gate 6), RECONSTRUCTED by lane R1 on 2026-09-13 from F1's STATUS after an iCloud
# rollback (F1's applier is evicted, so this is the documented behaviour, not F1's byte-identical text): the floor is
# checked against the file's own control table, never a typed number: it must be the largest per-SD control hazard
# ratio in negcontrols_final.json and the driver that control. Disagreement with numbers/bdsp_diseases_v3.csv
# (step 111) is LOGGED as a FLOOR NOTE, never asserted.
_NC_HR = {c["label"]: float(c["per_sd"]["hr"]) for c in _NCF["controls"].values()}
assert len(_NC_HR) == int(_NCF["n_controls"]) == len(NEGATIVE_CONTROLS), (_NC_HR, _NCF["n_controls"], NEGATIVE_CONTROLS)
assert STANDING_FLOOR == max(_NC_HR.values()) and _NC_HR[STANDING_FLOOR_SETBY] == STANDING_FLOOR, \
    (STANDING_FLOOR, STANDING_FLOOR_SETBY, _NC_HR)
_BD = pd.read_csv(f"{NUMBERS_DIR}/bdsp_diseases_v3.csv")
_BD = _BD[_BD.negative_control.astype(str).str.lower().isin(["1", "1.0", "true"])]
_BD_FLOOR = float(_BD.t90_hr.max()); _BD_SETBY = str(_BD.loc[_BD.t90_hr.idxmax(), "disease"])
print(f"FLOOR NOTE: standing floor {STANDING_FLOOR} set by {STANDING_FLOOR_SETBY} (negcontrols_final.json) vs largest "
      f"control t90_hr {_BD_FLOOR:.3f} set by {_BD_SETBY} (bdsp_diseases_v3.csv): "
      + ("agree" if abs(STANDING_FLOOR - _BD_FLOOR) < 5e-4 else "DIFFER, logged not asserted"))
MIN_EVENTS = 20                 # an outcome is fitted in a set only above this many incident events
_RK_V8 = __import__("pandas").read_csv(f"{NUMBERS_DIR}/ranking_v3.csv", comment="#").set_index("feature")
_RK_ROWS = ["TST_min", "spo2_pct_below_90", "AHI"]   # the three ranking rows this script reads; the sweep's whole-table std() fails on v8 (C_coupling_* rows carry their own baseline): narrowed by F1, reconstructed by lane R1 2026-09-13
assert (_RK_V8.loc[_RK_ROWS, "mC"] - _RK_V8.loc[_RK_ROWS, "dC"]).std() < 1e-9, "ranking rows disagree on the implied baseline"
PUBLISHED_TST_GAIN = float(_RK_V8.loc["TST_min", "dC"])  # v8 sweep: read from the current ranking, was a typed literal
PUBLISHED_T90_GAIN = float(_RK_V8.loc["spo2_pct_below_90", "dC"])  # v8 sweep: read from the current ranking
PUBLISHED_AHI_GAIN = float(_RK_V8.loc["AHI", "dC"])  # v8 sweep: read from the current ranking

TAG = cohort_tag()   # v8.1g: "_fullnights" under T90_FULL_NIGHTS_ONLY (step 237), "" for the all-nights run (step 150)
LOG = open(f"{HERE}/model_log{TAG}.txt", "w")


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    LOG.write(s + "\n")
    LOG.flush()


# ---------------------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------------------
t_start = time.time()
d = pd.read_parquet(f"{HERE}/analysis.parquet")
say(f"loaded {len(d):,} x {d.shape[1]} from analysis.parquet")
N_ALL_NIGHTS = len(d)
N_PRIMARY_ALL = int((d.is_primary_set == 1).sum())
if FULL_NIGHTS_ONLY:   # v8.1g (15 Sept 2026, step 237): the concordance family runs on the full diagnostic nights
    _flag = pd.read_parquet(T90_FINAL, columns=["BDSPPatientID", SPLIT_FLAG])
    assert _flag.BDSPPatientID.is_unique, "t90_final.parquet has duplicate patients"
    d = d.merge(_flag, on="BDSPPatientID", how="left", validate="1:1")
    assert d[SPLIT_FLAG].notna().all(), "split-night flag missing on some analysis rows"
    d = d[~d[SPLIT_FLAG].astype(bool)].drop(columns=[SPLIT_FLAG]).copy()
    say(f"[v8.1g] T90_FULL_NIGHTS_ONLY: {N_ALL_NIGHTS - len(d):,} split-night recordings dropped, "
        f"{len(d):,} full diagnostic nights kept; every output carries the suffix {TAG!r}")
assert len(d) == COHORT_N, len(d)

P = d[d.is_primary_set == 1].copy()
with open(f"{HERE}/summary.json") as _fh:   # step 149 (build.py) output with sidecar: the v2 primary set size, never retyped (v8 F1 fix, reconstructed by lane R1 2026-09-13)
    N_PRIMARY_REF = int(json.load(_fh)["step1"]["n_with_habitual"])
if FULL_NIGHTS_ONLY:   # v8.1g: the build's count holds on the all-nights frame; the full-nights primary set is that count less the split-night habitual patients
    assert N_PRIMARY_ALL == N_PRIMARY_REF, (N_PRIMARY_ALL, N_PRIMARY_REF)
    say(f"[v8.1g] primary set on the full nights {len(P):,} = the build's {N_PRIMARY_REF:,} habitual patients "
        f"minus {N_PRIMARY_REF - len(P):,} whose recording was a split night")
else:
    assert len(P) == N_PRIMARY_REF, (len(P), N_PRIMARY_REF)   # v2 primary set (v1 was 8711)
say(f"primary set {len(P):,}   of the {len(d):,} analysis cohort "
    f"({len(P) / len(d) * 100:.1f}%)")


def rint(x):
    """Rank-inverse-normal transform, the transform the published hazard ratios are on."""
    r = stats.rankdata(x)
    return stats.norm.ppf((r - 0.375) / (len(r) + 0.25))


# normal scores computed once, on the primary set, within hospital
P["hab_z"] = P.groupby("site_id").hours.transform(rint)
P["tst_z"] = P.groupby("site_id").TST_min.transform(
    lambda x: rint(x.fillna(x.median())))
P["t90_z"] = P.groupby("site_id").spo2_pct_below_90.transform(rint)
P["hab_hour"] = P.hours
P["hab_sd"] = (P.hours - P.hours.mean()) / P.hours.std()
P["tst_hour"] = P.TST_min / 60.0
HOURS_SD = float(P.hours.std())
say(f"habitual hours: mean {P.hours.mean():.3f}  SD {HOURS_SD:.3f}  median "
    f"{P.hours.median():.3f}  range {P.hours.min():.2f} to {P.hours.max():.2f}")
say(f"laboratory TST hours: mean {P.tst_hour.mean():.3f}  SD {P.tst_hour.std():.3f}   "
    f"missing {int(P.TST_min.isna().sum())}")
say(f"laboratory T90 %: median {P.spo2_pct_below_90.median():.3f}   "
    f"missing {int(P.spo2_pct_below_90.isna().sum())}")

# the two patients without a laboratory total sleep time are kept in the categorical model, which
# does not use it, and drop out of every continuous model, which is stated rather than hidden
N_NO_TST = int(P.TST_min.isna().sum())

HAB_LEVELS = ["<5h", "5-<6h", "6-<7h(ref)", ">=7h"]
assert set(P.hab_cat.dropna()) == set(HAB_LEVELS), sorted(set(P.hab_cat.dropna()))
for lev, tag in (("<5h", "hab_lt5"), ("5-<6h", "hab_5to6"), (">=7h", "hab_ge7")):
    P[tag] = (P.hab_cat == lev).astype(int)
CAT_TERMS = ["hab_lt5", "hab_5to6", "hab_ge7"]
say("\nhabitual bands, primary set")
for lev in HAB_LEVELS:
    say(f"   {lev:<12} {int((P.hab_cat == lev).sum()):>6,}")

# the twelve cross cells, habitual band by oxygenation band
O2_LEVELS = ["normal T90<=1%", "intermediate 1-10%", "low T90>10%"]
assert set(P.o2_cat) == set(O2_LEVELS), sorted(set(P.o2_cat))
P["cell"] = P.hab_cat.astype(str) + " x " + P.o2_cat.astype(str)
CELL_REF = "6-<7h(ref) x normal T90<=1%"
CELL_LEVELS = [f"{h} x {o}" for h in HAB_LEVELS for o in O2_LEVELS]
CELL_TERMS = []
for c in CELL_LEVELS:
    if c == CELL_REF:
        continue
    tag = "cell_" + (c.replace(" x ", "_X_").replace(" ", "").replace("<", "lt")
                     .replace(">", "gt").replace("=", "e").replace("%", "")
                     .replace("(", "").replace(")", "").replace("-", "_"))
    P[tag] = (P.cell == c).astype(int)
    CELL_TERMS.append((tag, c))

# ---------------------------------------------------------------------------------------
# outcomes
# ---------------------------------------------------------------------------------------
ALL_KEYS = [k for k in DISEASES if f"{k}_years" in P.columns] + ["death"]
LABEL = {k: v[0] for k, v in DISEASES.items()}
LABEL["death"] = "Death from any cause"
ORGAN = dict(ORGAN_GROUP)
ORGAN["death"] = "Composite"
NONCIRC = [k for k in ALL_KEYS if k not in CIRCULAR]
say(f"\noutcomes carried: {len(ALL_KEYS)}   non-circular: {len(NONCIRC)}   "
    f"circular dropped: {sorted(CIRCULAR)}")
say(f"negative controls in the definitions file: {sorted(NEGATIVE_CONTROLS)}")
say(f"the panel of the definitions file is carried in full: {', '.join(NEGATIVE_CONTROLS)} (v8.1 swap of 14 Sept). "
    "Fracture is an outcome, not a control. "
    "Osteoarthritis is an outcome, not a control. Thyroid disease is not a control.")

ADJ = ["age_s0", "age_s1", "age_s2", "male"]


def add_spline(f):
    """3-column natural cubic age basis. cr(df=4) then drop one column, never all four."""
    sp = dmatrix("cr(a, df=4) - 1", {"a": f.AgeAtVisit.values},
                 return_type="dataframe").iloc[:, 1:]
    assert sp.shape[1] == 3, sp.shape
    out = f.copy()
    for i in range(3):
        out[f"age_s{i}"] = sp.iloc[:, i].values
    return out


def fit(frame, key, terms, require=()):
    """One unpenalized site-stratified Cox. Returns a list of per-term records."""
    yc, ec, pc = f"{key}_years", f"{key}_incident", f"{key}_prevalent"
    f = frame[(frame[pc] == 0) & frame[yc].notna() & (frame[yc] > 0)]
    cols = list(terms) + ADJ
    dd = pd.DataFrame({"T": f[yc].values, "E": f[ec].astype(int).values,
                       "site": f.site_id.values,
                       **{c: f[c].values for c in set(cols) | set(require)}}).dropna()
    ev, n = int(dd.E.sum()), len(dd)
    if ev < MIN_EVENTS:
        return [], ev, n, "below the event floor"
    keep = cols + ["T", "E"]
    try:
        if dd.site.nunique() > 1:
            c = CoxPHFitter().fit(dd[keep + ["site"]], "T", "E", strata=["site"])
        else:
            c = CoxPHFitter().fit(dd[keep], "T", "E")
    except Exception as e:
        return [], ev, n, f"fit failed: {type(e).__name__}"
    out = []
    for t in terms:
        r = c.summary.loc[t]
        out.append({"term": t, "hr": float(r["exp(coef)"]),
                    "lo": float(r["exp(coef) lower 95%"]),
                    "hi": float(r["exp(coef) upper 95%"]), "p": float(r["p"])})
    return out, ev, n, "ok"


rows = []


def run(set_name, frame, model, terms, keys=NONCIRC, require=()):
    """Fit `model` across outcomes and append to the master row list."""
    fr = add_spline(frame)
    nfit = 0
    for key in keys:
        res, ev, n, status = fit(fr, key, terms, require=require)
        if not res:
            rows.append({"set": set_name, "n_set": len(frame), "model": model,
                         "outcome_key": key, "outcome": LABEL[key],
                         "organ_group": ORGAN[key],
                         "negative_control": key in NEGATIVE_CONTROLS,
                         "circular": key in CIRCULAR, "n": n, "events": ev,
                         "term": np.nan, "hr": np.nan, "lo": np.nan, "hi": np.nan,
                         "p": np.nan, "status": status})
            continue
        nfit += 1
        for r in res:
            rows.append({"set": set_name, "n_set": len(frame), "model": model,
                         "outcome_key": key, "outcome": LABEL[key],
                         "organ_group": ORGAN[key],
                         "negative_control": key in NEGATIVE_CONTROLS,
                         "circular": key in CIRCULAR, "n": n, "events": ev,
                         **r, "status": "ok"})
    say(f"   {set_name:<14} {model:<22} fitted {nfit:>2} of {len(keys)} outcomes")


# ---------------------------------------------------------------------------------------
# 1 and 2. the primary, categorical and continuous
# ---------------------------------------------------------------------------------------
say("\n" + "=" * 92)
say("PRIMARY  n = 8,711, unpenalized Cox, stratified by hospital, age spline (3 columns) + sex")
say("=" * 92)

run("primary", P, "habitual_category", CAT_TERMS)

# every continuous model runs on the identical patients, so the head-to-head is a like-for-like
CONT = P[P.TST_min.notna() & P.spo2_pct_below_90.notna() & P.hours.notna()].copy()
say(f"\ncontinuous models run on {len(CONT):,} patients, the {len(P):,} primary set less the "
    f"{N_NO_TST} without a laboratory total sleep time, so habitual, laboratory sleep and "
    f"laboratory oxygen are estimated in exactly the same people")
for term in ["hab_z", "hab_hour", "hab_sd", "tst_z", "tst_hour", "t90_z"]:
    run("primary", CONT, f"continuous_{term}", [term])

# the categorical model repeated on the continuous frame, so nothing rests on those two patients
run("primary_contframe", CONT, "habitual_category", CAT_TERMS)

# ---------------------------------------------------------------------------------------
# 7. measurement quality and timing
# ---------------------------------------------------------------------------------------
say("\n" + "=" * 92)
say("SENSITIVITY SETS")
say("=" * 92)
SETS = {
    "pre_study": P[P.pre_study_only == 1],
    "tier_AB": P[P.is_tier_AB == 1],
    "tier_A": P[P.is_tier_A == 1],
    "template": P[P.is_template == 1],
}
for name, sub in SETS.items():
    sites = sub.site_id.nunique()
    say(f"\n{name}: n = {len(sub):,}   hospitals {sites}   "
        f"median habitual hours {sub.hours.median():.2f}")
    if sites == 1:
        say(f"   only one hospital present, so this set is fitted WITHOUT hospital "
            f"stratification: {sub.site_id.iloc[0]}")
    run(name, sub, "habitual_category", CAT_TERMS)
    sc = sub[sub.TST_min.notna()]
    run(name, sc, "continuous_hab_z", ["hab_z"])
    run(name, sc, "continuous_hab_hour", ["hab_hour"])
    run(name, sc, "continuous_tst_z", ["tst_z"])
    run(name, sc, "continuous_t90_z", ["t90_z"])

# ---------------------------------------------------------------------------------------
# the twelve cells Alen asked for, habitual band crossed with laboratory oxygenation
# ---------------------------------------------------------------------------------------
say("\n" + "=" * 92)
say("CROSS CELLS  habitual band x laboratory oxygenation, reference "
    f"{CELL_REF}")
say("=" * 92)
run("primary", CONT, "cross_cells", [t for t, _ in CELL_TERMS])

R = pd.DataFrame(rows)

# ---------------------------------------------------------------------------------------
# 5 and 6. floors and Benjamini-Hochberg
# ---------------------------------------------------------------------------------------
R["effect"] = np.maximum(R.hr, 1 / R.hr)
R["direction"] = np.where(R.hr.isna(), "", np.where(R.hr > 1, "higher risk", "lower risk"))
R["significant"] = (R.lo > 1) | (R.hi < 1)
R.loc[R.hr.isna(), "significant"] = False

# per-SD terms are the only ones the standing floor is defined on
PER_SD_TERMS = {"hab_z", "hab_sd", "tst_z", "t90_z"}
R["floor_standing"] = np.where(R.term.isin(PER_SD_TERMS), STANDING_FLOOR, np.nan)
R["passes_floor_standing"] = np.where(R.term.isin(PER_SD_TERMS),
                                      R.effect > STANDING_FLOOR, np.nan)

# the empirical floor: the largest control effect in this same set, model and term
emp = (R[R.negative_control & R.hr.notna()]
       .groupby(["set", "model", "term"]).effect.max().rename("floor_empirical"))
n_ctrl = (R[R.negative_control & R.hr.notna()]
          .groupby(["set", "model", "term"]).size().rename("n_controls"))
R = R.merge(emp, on=["set", "model", "term"], how="left")
R = R.merge(n_ctrl, on=["set", "model", "term"], how="left")
R["passes_floor_empirical"] = R.effect > R.floor_empirical
R.loc[R.hr.isna(), ["passes_floor_empirical"]] = np.nan

# BH within each set, model and term. q_all uses every outcome fitted, q_main excludes the
# negative controls because a control is a diagnostic and not a hypothesis
R["q_all"] = np.nan
R["q_main"] = np.nan
for (s, m, t), g in R[R.p.notna()].groupby(["set", "model", "term"]):
    R.loc[g.index, "q_all"] = multipletests(g.p.values, method="fdr_bh")[1]
    gm = g[~g.negative_control]
    if len(gm) > 1:
        R.loc[gm.index, "q_main"] = multipletests(gm.p.values, method="fdr_bh")[1]
R["bh_survivor"] = R.q_main < 0.05
R.loc[R.q_main.isna(), "bh_survivor"] = False
R["survives_everything"] = R.significant & R.bh_survivor & (
    R.passes_floor_empirical.fillna(False).astype(bool))

COLS = ["set", "n_set", "model", "term", "outcome_key", "outcome", "organ_group",
        "negative_control", "circular", "n", "events", "hr", "lo", "hi", "p",
        "q_all", "q_main", "effect", "direction", "significant",
        "floor_standing", "passes_floor_standing", "floor_empirical", "n_controls",
        "passes_floor_empirical", "bh_survivor", "survives_everything", "status"]
R[COLS].to_csv(f"{HERE}/results_habitual{TAG}.csv", index=False)
say(f"\nwrote results_habitual.csv   {len(R):,} rows")

# ---------------------------------------------------------------------------------------
# 3. the head-to-head, three exposures, one set of patients
# ---------------------------------------------------------------------------------------
say("\n" + "=" * 92)
say("HEAD TO HEAD  habitual sleep vs laboratory sleep vs laboratory oxygen, "
    f"same {len(CONT):,} patients, per SD")
say("=" * 92)

h2h = {}
for term, tag in (("hab_z", "hab"), ("tst_z", "tst"), ("t90_z", "t90"),
                  ("hab_hour", "habhr"), ("tst_hour", "tsthr")):
    g = R[(R["set"] == "primary") & (R.term == term) & R.hr.notna()].set_index("outcome_key")
    for c in ["hr", "lo", "hi", "p", "q_main", "significant", "passes_floor_empirical"]:
        h2h[f"{tag}_{c}"] = g[c]
H = pd.DataFrame(h2h)
meta = (R[(R["set"] == "primary") & (R.term == "hab_z")]
        .set_index("outcome_key")[["outcome", "organ_group", "negative_control", "n", "events"]])
H = meta.join(H).sort_values("events", ascending=False)
H.to_csv(f"{HERE}/results_habitual_headtohead{TAG}.csv")
say(f"wrote results_habitual_headtohead.csv   {len(H)} outcomes")

fitted = H[H.hab_hr.notna()]
say(f"\n{'condition':<28}{'events':>7}   {'habitual /SD':>20}{'lab TST /SD':>22}"
    f"{'lab T90 /SD':>22}")
say("-" * 102)
for k, r in fitted.iterrows():
    star = lambda s: "*" if s else " "  # noqa: E731
    say(f"{r.outcome[:27]:<28}{int(r.events):>7}   "
        f"{r.hab_hr:>6.3f} [{r.hab_lo:.3f}-{r.hab_hi:.3f}]{star(r.hab_significant)}"
        f"{r.tst_hr:>7.3f} [{r.tst_lo:.3f}-{r.tst_hi:.3f}]{star(r.tst_significant)}"
        f"{r.t90_hr:>7.3f} [{r.t90_lo:.3f}-{r.t90_hi:.3f}]{star(r.t90_significant)}")

# ---------------------------------------------------------------------------------------
# 4. the count that answers the question
# ---------------------------------------------------------------------------------------
say("\n" + "=" * 92)
say("SCOREBOARD  non-circular outcomes fitted in the same patients")
say("=" * 92)
main = fitted[~fitted.negative_control]
score = {}
for tag, nice in (("hab", "habitual home sleep"), ("tst", "laboratory total sleep time"),
                  ("t90", "laboratory T90")):
    eff = np.maximum(main[f"{tag}_hr"], 1 / main[f"{tag}_hr"])
    sig = int(main[f"{tag}_significant"].sum())
    bh = int((main[f"{tag}_q_main"] < 0.05).sum())
    fl = int((main[f"{tag}_significant"] & main[f"{tag}_passes_floor_empirical"]).sum())
    std = int((main[f"{tag}_significant"] & (eff > STANDING_FLOOR)).sum())
    allthree = int((main[f"{tag}_significant"] & (main[f"{tag}_q_main"] < 0.05)
                    & main[f"{tag}_passes_floor_empirical"]).sum())
    score[nice] = {"outcomes_tested": int(len(main)), "significant": sig,
                   "bh_survivors": bh, "significant_and_above_empirical_floor": fl,
                   "significant_and_above_standing_floor_1133": std,
                   "survives_everything": allthree,
                   "median_abs_hr_per_sd": round(float(eff.median()), 4),
                   "largest_abs_hr_per_sd": round(float(eff.max()), 4)}
    say(f"{nice:<30} {sig:>3} of {len(main)} significant, {bh:>3} survive Benjamini-Hochberg, "
        f"{fl:>3} clear the control floor, {std:>3} clear the standing {STANDING_FLOOR} floor, "
        f"{allthree:>3} clear all three. "
        f"median effect {eff.median():.3f} per SD, largest {eff.max():.3f}")

say("\nnegative controls in the same three models, per SD")
ctrl = fitted[fitted.negative_control]
say(f"{'control':<24}{'events':>7}{'habitual':>22}{'lab TST':>22}{'lab T90':>22}")
for k, r in ctrl.iterrows():
    say(f"{r.outcome[:23]:<24}{int(r.events):>7}"
        f"{r.hab_hr:>8.3f} [{r.hab_lo:.3f}-{r.hab_hi:.3f}]"
        f"{r.tst_hr:>8.3f} [{r.tst_lo:.3f}-{r.tst_hi:.3f}]"
        f"{r.t90_hr:>8.3f} [{r.t90_lo:.3f}-{r.t90_hi:.3f}]")
floors = {}
for tag, nice in (("hab", "habitual"), ("tst", "lab TST"), ("t90", "lab T90")):
    f = float(np.maximum(ctrl[f"{tag}_hr"], 1 / ctrl[f"{tag}_hr"]).max())
    floors[nice] = round(f, 4)
    nsig = int(ctrl[f"{tag}_significant"].sum())
    say(f"empirical floor, {nice:<10} {f:.3f}   controls significant: {nsig} of {len(ctrl)}   "
        f"(standing per-SD floor {STANDING_FLOOR})")

say("\nstruck by the floors: significant results that do not clear the control noise")
struck = []
for tag, nice in (("hab", "habitual"), ("tst", "lab TST"), ("t90", "lab T90")):
    eff = np.maximum(main[f"{tag}_hr"], 1 / main[f"{tag}_hr"])
    bad = main[main[f"{tag}_significant"] & ~main[f"{tag}_passes_floor_empirical"].astype(bool)]
    bad_std = main[main[f"{tag}_significant"] & (eff <= STANDING_FLOOR)]
    struck.append({"exposure": nice,
                   "struck_by_empirical_floor": sorted(bad.outcome.tolist()),
                   "struck_by_standing_floor": sorted(bad_std.outcome.tolist())})
    say(f"   {nice:<10} empirical floor {floors[nice]:.3f} strikes "
        f"{len(bad)}: {', '.join(bad.outcome.tolist()) or 'none'}")
    say(f"   {'':<10} standing floor {STANDING_FLOOR} strikes "
        f"{len(bad_std)}: {', '.join(bad_std.outcome.tolist()) or 'none'}")

# where the two sleep measures point in opposite directions in the same patients
opp = main[((main.hab_hr - 1) * (main.tst_hr - 1)) < 0]
say(f"\nhabitual sleep and laboratory sleep point in OPPOSITE directions for "
    f"{len(opp)} of {len(main)} outcomes")
for k, r in opp.sort_values("events", ascending=False).head(12).iterrows():
    say(f"   {r.outcome[:30]:<32}habitual {r.hab_hr:.3f}   laboratory {r.tst_hr:.3f}")
sp_r, sp_p = stats.spearmanr(CONT.hours, CONT.TST_min)
pe_r = float(np.corrcoef(CONT.hours, CONT.TST_min)[0, 1])
say(f"\nwithin these {len(CONT):,} patients the habitual and laboratory nights correlate at "
    f"Spearman {sp_r:+.4f} (p {sp_p:.3g}), Pearson {pe_r:+.4f}. Median laboratory hours by "
    f"habitual band: " + ", ".join(
        f"{lev} {CONT[CONT.hab_cat == lev].tst_hour.median():.2f}" for lev in HAB_LEVELS))
say("the same correlation by measurement quality, which is where the quoted -0.024 comes from")
corrs = {}
for nm, sub in (("primary", CONT), ("tier A", CONT[CONT.is_tier_A == 1]),
                ("template", CONT[CONT.is_template == 1]),
                ("tier B", CONT[CONT.tier == "B"]), ("tier C", CONT[CONT.tier == "C"])):
    r1, p1 = stats.spearmanr(sub.hours, sub.TST_min)
    r2 = float(np.corrcoef(sub.hours, sub.TST_min)[0, 1])
    corrs[nm] = {"n": int(len(sub)), "spearman": round(float(r1), 4),
                 "spearman_p": float(p1), "pearson": round(r2, 4)}
    say(f"   {nm:<10} n {len(sub):>5,}   Spearman {r1:+.4f} (p {p1:.3g})   Pearson {r2:+.4f}")

# per hour, the scale a clinician reads
say(f"\nper +1 hour of sleep, the two measures side by side, {len(CONT):,} patients")
say(f"{'condition':<30}{'events':>7}{'habitual +1h':>24}{'laboratory +1h':>24}")
hh = H[H.habhr_hr.notna()].sort_values("events", ascending=False)
for k, r in hh.iterrows():
    say(f"{r.outcome[:29]:<30}{int(r.events):>7}   {r.habhr_hr:>6.3f} "
        f"[{r.habhr_lo:.3f}-{r.habhr_hi:.3f}]"
        f"{'*' if r.habhr_significant else ' '}   {r.tsthr_hr:>6.3f} "
        f"[{r.tsthr_lo:.3f}-{r.tsthr_hi:.3f}]{'*' if r.tsthr_significant else ' '}")

# ---------------------------------------------------------------------------------------
# 2 continued. the held-out concordance gain, the published ranking's own currency
# ---------------------------------------------------------------------------------------
say("\n" + "=" * 92)
say("HELD-OUT CONCORDANCE  published protocol: fit in one hospital, score in the other")
say("=" * 92)

C = CONT.copy()
sp4 = dmatrix("cr(a, df=4) - 1", {"a": C.AgeAtVisit.values}, return_type="dataframe")
sp4.columns = [f"aa{i}" for i in range(4)]
for c in sp4.columns:
    C[c] = sp4[c].values
BASE4 = [f"aa{i}" for i in range(4)] + ["male"]


def heldout(feat, key):
    """Returns (mean C, number of the two train/test directions that were scorable)."""
    yc, ec, pc = f"{key}_years", f"{key}_incident", f"{key}_prevalent"
    g = C[(C[pc] == 0) & C[yc].notna() & (C[yc] > 0)]
    cols = BASE4 + ([feat] if feat else [])
    res = []
    for tr, te in (("I0002", "I0006"), ("I0006", "I0002")):
        A = g[g.site_id == tr][cols + [yc, ec]].dropna()
        B = g[g.site_id == te][cols + [yc, ec]].dropna()
        if A[ec].sum() < 30 or B[ec].sum() < 30:
            res.append(np.nan)
            continue
        try:
            m = CoxPHFitter(penalizer=0.01).fit(A, duration_col=yc, event_col=ec)
            res.append(concordance_index(B[yc], -m.predict_partial_hazard(B), B[ec]))
        except Exception:
            res.append(np.nan)
    ok = int(sum(not np.isnan(x) for x in res))
    return (float(np.nanmean(res)) if ok else np.nan), ok


conc_rows = []
for key in NONCIRC:
    base, ndir = heldout(None, key)
    if np.isnan(base):
        continue
    rec = {"outcome_key": key, "outcome": LABEL[key],
           "negative_control": key in NEGATIVE_CONTROLS, "base_c": base,
           "directions_scorable": ndir}
    for feat in ("hab_z", "tst_z", "t90_z"):
        c, _ = heldout(feat, key)
        rec[f"{feat}_c"] = c
        rec[f"{feat}_gain"] = c - base
    conc_rows.append(rec)
CC = pd.DataFrame(conc_rows)
CC.to_csv(f"{HERE}/results_habitual_concordance{TAG}.csv", index=False)
N_BOTH = int((CC.directions_scorable == 2).sum())
say(f"outcomes scorable across both hospitals: {len(CC)} of {len(NONCIRC)}, of which "
    f"{N_BOTH} were scorable in BOTH directions and {len(CC) - N_BOTH} in only one. "
    f"The protocol needs 30 incident events in each hospital and only "
    f"{int((CONT.site_id == 'I0006').sum()):,} of the {len(CONT):,} sit at I0006.")
say(f"baseline held-out concordance, age and sex only: {CC.base_c.mean():.4f}")
conc = {}
for feat, nice in (("hab_z", "habitual home sleep"), ("tst_z", "laboratory total sleep time"),
                   ("t90_z", "laboratory T90")):
    g = CC[f"{feat}_gain"]
    conc[nice] = {"mean_heldout_c": round(float(CC[f"{feat}_c"].mean()), 6),
                  "mean_gain": round(float(g.mean()), 6),
                  "worst_gain": round(float(g.min()), 6),
                  "outcomes_improved": int((g > 0).sum()), "outcomes_tested": int(g.notna().sum())}
    say(f"{nice:<30} mean C {CC[f'{feat}_c'].mean():.4f}   gain {g.mean():+.6f}   "
        f"worst {g.min():+.6f}   improved {int((g > 0).sum())} of {int(g.notna().sum())}")
say(f"\npublished all-cohort reference values, {N_MEASURES}-measure ranking:")
say(f"   laboratory total sleep time  {PUBLISHED_TST_GAIN:+.6f}  rank {int(_RK_V8.loc['TST_min', 'rank'])} of {N_MEASURES}")
say(f"   apnea-hypopnea index         {PUBLISHED_AHI_GAIN:+.6f}  rank {int(_RK_V8.loc['AHI', 'rank'])} of {N_MEASURES}")
say(f"   T90                          {PUBLISHED_T90_GAIN:+.6f}  rank {int(_RK_V8.loc['spo2_pct_below_90', 'rank'])} of {N_MEASURES}")

# ---------------------------------------------------------------------------------------
# 1 continued. the categorical table
# ---------------------------------------------------------------------------------------
say("\n" + "=" * 92)
say("HABITUAL BANDS, primary set, reference 6 to under 7 hours")
say("=" * 92)
cat = R[(R["set"] == "primary") & (R.model == "habitual_category") & R.hr.notna()]
wide = cat.pivot_table(index=["outcome_key", "outcome", "events"], columns="term",
                       values=["hr", "lo", "hi"])
say(f"{'condition':<28}{'events':>7}{'<5h':>24}{'5-<6h':>24}{'>=7h':>24}")
say("-" * 108)
for (k, lab, ev), r in wide.iterrows():
    def cell(t):
        return f"{r[('hr', t)]:>6.2f} [{r[('lo', t)]:.2f}-{r[('hi', t)]:.2f}]"
    say(f"{lab[:27]:<28}{int(ev):>7}      {cell('hab_lt5')}      {cell('hab_5to6')}"
        f"      {cell('hab_ge7')}")

# ---------------------------------------------------------------------------------------
# 7 continued. does the answer move with measurement quality
# ---------------------------------------------------------------------------------------
say("\n" + "=" * 92)
say("MEASUREMENT QUALITY  does the answer change as the exposure gets cleaner")
say("=" * 92)
qual = []
for s in ["primary", "pre_study", "tier_AB", "tier_A", "template"]:
    g = R[(R["set"] == s) & (R.term == "hab_z") & R.hr.notna() & ~R.negative_control]
    gc = R[(R["set"] == s) & (R.model == "habitual_category") & R.hr.notna()
           & ~R.negative_control]
    if not len(g):
        continue
    nset = int(R[R["set"] == s].n_set.iloc[0])
    rec = {"set": s, "n_set": nset, "outcomes_fitted": int(len(g)),
           "significant_per_sd": int(g.significant.sum()),
           "median_abs_hr_per_sd": round(float(g.hr.apply(lambda x: max(x, 1 / x)).median()), 4),
           "largest_abs_hr_per_sd": round(float(g.hr.apply(lambda x: max(x, 1 / x)).max()), 4),
           "cat_terms_fitted": int(len(gc)),
           "cat_terms_significant": int(gc.significant.sum())}
    qual.append(rec)
    say(f"{s:<12} n {nset:>6,}   {rec['outcomes_fitted']:>2} outcomes fitted, "
        f"{rec['significant_per_sd']:>2} significant per SD, median effect "
        f"{rec['median_abs_hr_per_sd']:.3f}, largest {rec['largest_abs_hr_per_sd']:.3f}   "
        f"|  categorical: {rec['cat_terms_significant']} of {rec['cat_terms_fitted']} "
        f"band terms significant")
pd.DataFrame(qual).to_csv(f"{HERE}/results_habitual_quality{TAG}.csv", index=False)

# ---------------------------------------------------------------------------------------
# the twelve cells
# ---------------------------------------------------------------------------------------
say("\n" + "=" * 92)
say(f"TWELVE CELLS  reference {CELL_REF}, the four Alen named are marked")
say("=" * 92)
NAMED = {"<5h x normal T90<=1%", "<5h x low T90>10%",
         ">=7h x normal T90<=1%", ">=7h x low T90>10%"}
cells = R[(R.model == "cross_cells") & R.hr.notna()]
lookup = {t: c for t, c in CELL_TERMS}
for key in ["cvd", "death", "htn2", "diabetes", "hf", "dementia"]:
    g = cells[cells.outcome_key == key]
    if not len(g):
        continue
    say(f"\n{LABEL[key]}   {int(g.events.iloc[0]):,} events in {int(g.n.iloc[0]):,} at risk")
    for _, r in g.iterrows():
        c = lookup[r.term]
        mark = " <-- named" if c in NAMED else ""
        say(f"   {c:<36}{r.hr:>7.2f} [{r.lo:.2f}-{r.hi:.2f}]  p {r.p:.3g}{mark}")
cells.to_csv(f"{HERE}/results_habitual_cells{TAG}.csv", index=False)

# ---------------------------------------------------------------------------------------
# body mass index. The four habitual results that survive everything are obesity, hypertension,
# diabetes and reflux, which is the exact set body mass index would produce on its own. It is
# recorded for only a quarter of these patients, so the adjusted and unadjusted models are both
# refitted on that quarter, otherwise the comparison would confound adjustment with sample change.
# ---------------------------------------------------------------------------------------
say("\n" + "=" * 92)
say("BODY MASS INDEX  same patients, with and without the adjustment")
say("=" * 92)
BM = CONT[CONT.bmi.notna()].copy()
BM["bmi_z"] = (BM.bmi - BM.bmi.mean()) / BM.bmi.std()
say(f"body mass index is recorded for {len(BM):,} of the {len(CONT):,} continuous set "
    f"({len(BM) / len(CONT) * 100:.1f}%), median {BM.bmi.median():.1f}")
BMS = add_spline(BM)
bmi_rows = []
for key in ["obesity", "htn2", "diabetes", "gerd", "alopecia", "cvd", "death", "hf"]:   # v8.1: alopecia took back pain's slot
    yc, ec, pc = f"{key}_years", f"{key}_incident", f"{key}_prevalent"
    f = BMS[(BMS[pc] == 0) & BMS[yc].notna() & (BMS[yc] > 0)]
    dd = pd.DataFrame({"T": f[yc].values, "E": f[ec].astype(int).values,
                       "site": f.site_id.values,
                       **{c: f[c].values for c in ["hab_z", "bmi_z"] + ADJ}}).dropna()
    if dd.E.sum() < MIN_EVENTS:
        say(f"   {LABEL[key]:<26} only {int(dd.E.sum())} events with body mass index, skipped")
        continue
    out = {"outcome_key": key, "outcome": LABEL[key], "n": len(dd), "events": int(dd.E.sum())}
    for tag, cols in (("unadj", ["hab_z"] + ADJ), ("adj", ["hab_z", "bmi_z"] + ADJ)):
        c = CoxPHFitter().fit(dd[cols + ["T", "E", "site"]], "T", "E", strata=["site"])
        r = c.summary.loc["hab_z"]
        out[f"{tag}_hr"] = round(float(r["exp(coef)"]), 4)
        out[f"{tag}_lo"] = round(float(r["exp(coef) lower 95%"]), 4)
        out[f"{tag}_hi"] = round(float(r["exp(coef) upper 95%"]), 4)
        out[f"{tag}_p"] = float(r["p"])
    full = R[(R["set"] == "primary") & (R.term == "hab_z") & (R.outcome_key == key)]
    out["full_set_hr"] = round(float(full.hr.iloc[0]), 4)
    bmi_rows.append(out)
    say(f"   {LABEL[key]:<26} ev {out['events']:>4}   full set {out['full_set_hr']:.3f}   "
        f"BMI subset unadjusted {out['unadj_hr']:.3f} "
        f"[{out['unadj_lo']:.3f}-{out['unadj_hi']:.3f}]   "
        f"BMI adjusted {out['adj_hr']:.3f} [{out['adj_lo']:.3f}-{out['adj_hi']:.3f}]")
BMI_TAB = pd.DataFrame(bmi_rows)
BMI_TAB.to_csv(f"{HERE}/results_habitual_bmi{TAG}.csv", index=False)
say("   obesity as an outcome adjusted for body mass index is close to circular and is shown "
    "only for completeness")

# ---------------------------------------------------------------------------------------
# proportional hazards, and an independent refit
# ---------------------------------------------------------------------------------------
say("\n" + "=" * 92)
say("CHECKS  proportional hazards, and every primary hazard ratio refitted in other software")
say("=" * 92)

from lifelines.statistics import proportional_hazard_test  # noqa: E402
from statsmodels.duration.hazard_regression import PHReg    # noqa: E402
import sys as _sys_v8; _sys_v8.path.insert(0, paths.ANALYSIS_DIR)  # v8 sweep 2026-09-12

CS = add_spline(CONT)
ph_rows = []
for term in ("hab_z", "tst_z", "t90_z"):
    viol = []
    for key in NONCIRC:
        yc, ec, pc = f"{key}_years", f"{key}_incident", f"{key}_prevalent"
        f = CS[(CS[pc] == 0) & CS[yc].notna() & (CS[yc] > 0)]
        dd = pd.DataFrame({"T": f[yc].values, "E": f[ec].astype(int).values,
                           "site": f.site_id.values,
                           **{c: f[c].values for c in [term] + ADJ}}).dropna()
        if dd.E.sum() < MIN_EVENTS:
            continue
        try:
            c = CoxPHFitter().fit(dd, "T", "E", strata=["site"])
            t = proportional_hazard_test(c, dd, time_transform="rank")
            p = float(t.summary.loc[term, "p"])
            ph_rows.append({"term": term, "outcome_key": key, "ph_p": p})
            if p < 0.05:
                viol.append((LABEL[key], p))
        except Exception:
            continue
    say(f"{term:<8} proportional hazards violated at p<0.05 in {len(viol)} of "
        f"{len([r for r in ph_rows if r['term'] == term])} outcomes"
        + (": " + ", ".join(f"{n} ({p:.3f})" for n, p in sorted(viol, key=lambda x: x[1])[:8])
           if viol else ""))
PH = pd.DataFrame(ph_rows)
PH.to_csv(f"{HERE}/results_habitual_ph{TAG}.csv", index=False)

# the same models refitted with statsmodels PHReg, a different solver and a different code path
say("\nindependent refit, statsmodels PHReg against lifelines CoxPHFitter, primary set")
ver, worst = [], 0.0
for term in ("hab_z", "hab_hour", "tst_z", "t90_z"):
    for key in ["cvd", "death", "htn2", "diabetes", "hf", "dementia", "obesity", "gerd",
                "alopecia", "stroke_any"]:   # v8.1
        yc, ec, pc = f"{key}_years", f"{key}_incident", f"{key}_prevalent"
        f = CS[(CS[pc] == 0) & CS[yc].notna() & (CS[yc] > 0)]
        dd = pd.DataFrame({"T": f[yc].values, "E": f[ec].astype(int).values,
                           "site": f.site_id.values,
                           **{c: f[c].values for c in [term] + ADJ}}).dropna()
        X = dd[[term] + ADJ].values
        m = PHReg(dd["T"].values, X, status=dd["E"].values,
                  strata=dd["site"].values, ties="breslow").fit()
        sm_hr = float(np.exp(m.params[0]))
        lf = R[(R["set"] == "primary") & (R.term == term) & (R.outcome_key == key)]
        lf_hr = float(lf.hr.iloc[0])
        diff = abs(sm_hr - lf_hr)
        worst = max(worst, diff)
        ver.append({"term": term, "outcome_key": key, "lifelines_hr": round(lf_hr, 6),
                    "statsmodels_hr": round(sm_hr, 6), "abs_diff": round(diff, 6)})
V = pd.DataFrame(ver)
V.to_csv(f"{HERE}/results_habitual_verify{TAG}.csv", index=False)
say(f"{len(V)} hazard ratios refitted in statsmodels. Largest absolute disagreement "
    f"{worst:.6f}. All agree to 3 decimals: {bool(worst < 5e-4)}")
assert worst < 5e-3, f"two solvers disagree by {worst}"

# read-back: the categorical model on the 8,711 against the same model on the 8,709
a = R[(R["set"] == "primary") & (R.model == "habitual_category") & R.hr.notna()].set_index(
    ["outcome_key", "term"]).hr
bb = R[(R["set"] == "primary_contframe") & (R.model == "habitual_category")
       & R.hr.notna()].set_index(["outcome_key", "term"]).hr
dmax = float((a - bb).abs().max())
say(f"dropping the {N_NO_TST} patients without a laboratory total sleep time moves the "
    f"categorical hazard ratios by at most {dmax:.5f}")

# ---------------------------------------------------------------------------------------
# summary
# ---------------------------------------------------------------------------------------
summary = {
    "cohort_n": int(len(d)),
    "cohort_tag": TAG,   # v8.1g
    "full_nights_only": bool(FULL_NIGHTS_ONLY),   # v8.1g
    "cohort_n_all_nights": int(N_ALL_NIGHTS),   # v8.1g
    "primary_n_all_nights": int(N_PRIMARY_ALL),   # v8.1g
    "primary_n": int(len(P)),
    "continuous_n": int(len(CONT)),
    "n_without_lab_tst": N_NO_TST,
    "outcomes_noncircular": len(NONCIRC),
    "outcomes_fitted_primary": int(len(main)) + int(len(ctrl)),
    "min_events": MIN_EVENTS,
    "standing_per_sd_floor": STANDING_FLOOR,
    "empirical_floors_per_sd": floors,
    "negative_controls_carried": sorted(NEGATIVE_CONTROLS),
    "negative_controls_missing_from_definitions": ["contact dermatitis", "hemorrhoids"],
    "scoreboard": score,
    "concordance": conc,
    "published_reference_gains": {"lab_tst": PUBLISHED_TST_GAIN, "ahi": PUBLISHED_AHI_GAIN,
                                  "t90": PUBLISHED_T90_GAIN},
    "measurement_quality": qual,
    "habitual_hours_sd": round(HOURS_SD, 4),
    "struck_by_floor": struck,
    "opposite_direction_hab_vs_lab": {"n": int(len(opp)), "of": int(len(main)),
                                      "outcomes": sorted(opp.outcome.tolist())},
    "hab_vs_lab_tst_correlation": {"spearman": round(float(sp_r), 4),
                                   "pearson": round(pe_r, 4), "n": int(len(CONT))},
    "hab_vs_lab_tst_correlation_by_set": corrs,
    "bmi": {"n_with_bmi": int(len(BM)), "pct": round(len(BM) / len(CONT) * 100, 1),
            "table": bmi_rows},
    "significant_expected_by_chance_at_005": round(len(main) * 0.05, 2),
    "concordance_directions_both": N_BOTH,
    "ph_violations_p05": {t: int(((PH.term == t) & (PH.ph_p < 0.05)).sum())
                          for t in PH.term.unique()},
    "solver_agreement_max_abs_diff": round(worst, 6),
    "categorical_max_shift_dropping_2": round(dmax, 6),
    "runtime_sec": round(time.time() - t_start, 1),
}
json.dump(summary, open(f"{HERE}/results_habitual_summary{TAG}.json", "w"), indent=1)
say(f"\nwrote results_habitual_summary.json   runtime {time.time() - t_start:.0f}s")
LOG.close()
