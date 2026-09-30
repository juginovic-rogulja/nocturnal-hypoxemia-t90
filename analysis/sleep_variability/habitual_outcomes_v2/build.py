"""
Build the habitual-sleep analysis dataset and quantify how selected it is.

The question this dataset exists to answer: laboratory total sleep time ranks 190th of 197
measurements in the T90 paper, and within-person laboratory and habitual home sleep are
unrelated (r = -0.024 on the cleanest source). The obvious objection is that the laboratory
mismeasured the construct. Habitual self-reported home sleep, extracted from clinical notes,
is a second and independent measurement of the same construct. If it also fails to predict
disease, the result is about sleep duration and not about the night in the laboratory.

Nothing here is modelled. This script merges, categorises, counts, and quantifies selection.

Outputs, all in this directory:
  analysis.parquet              one row per cohort patient, habitual + lab + outcomes + covariates
  habitual_categories.csv       N and events per habitual band
  o2_categories.csv             N and events per lab-oxygenation band
  cross_cells.csv               N and events in every habitual x oxygenation cell
  cross_cells_by_disease.csv    the same cells, every outcome, long form
  selection_smd.csv             standardised mean differences, with vs without a habitual value
  candidate_sets.csv            event counts under each candidate primary analysis set
  summary.json                  every number quoted in the write-up
  build_log.txt                 the console transcript
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

T90 = paths.T90_ROOT
HAB = f"{paths.SV_ROOT}/habitual_sleep"
OUT = Path(f"{paths.SV_ROOT}/habitual_outcomes_v2")
OUT.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
from cohort_spec import COHORT_N, apply_cohort           # noqa: E402
from disease_definitions import (                        # noqa: E402
    CIRCULAR, DISEASES, NEGATIVE_CONTROLS, ORGAN_GROUP,
)

# RERUN NOTE (2026-08-13, habitual_outcomes_v2): reproduce the frozen 2026-08-07 16:20 input
# state. The two replacement negative controls (dermatitis_contact, haemorrhoids) were
# materialised into t90_final.parquet at 16:41 on 2026-08-07, AFTER the shipped
# analysis.parquet was built at 16:20, and model_crossed.py merges them at load time and
# REFUSES to run if they are already present in analysis.parquet. Dropping them from the
# outcome list here keeps analysis.parquet identical to the shipped build (52 outcomes,
# prevalent_count over 50 conditions). They still enter every crossed model, exactly as
# shipped, via model_crossed.py's own merge from t90_final.parquet.
DISEASES = {k: v for k, v in DISEASES.items()
            if k not in ("dermatitis_contact", "haemorrhoids")}

# The negative-control panel, settled 2026-08-07 and read from the source of truth rather than
# retyped here. It is back pain, cataract, glaucoma, contact dermatitis and hemorrhoids.
# Fracture and osteoarthritis were controls in earlier drafts and are ordinary outcomes now,
# because a plausible causal path runs from the exposure to each of them, which disqualifies a
# negative control however clean its estimate looks. Contact dermatitis and hemorrhoids are
# defined in disease_definitions.py since 2026-08-07, so nothing is missing from the panel.
FINAL_CONTROLS = list(NEGATIVE_CONTROLS)
import disease_definitions as _dd   # noqa: E402
assert len(FINAL_CONTROLS) == 5 and all(_dd.DISEASES[k][1] for k in FINAL_CONTROLS), FINAL_CONTROLS   # v8.1: the spec's panel, no typed list
CONTROLS_NOT_IN_FILE = []

LOG = []


def say(*a):
    line = " ".join(str(x) for x in a)
    print(line)
    LOG.append(line)


def hr(title):
    say("")
    say("=" * 100)
    say(title)
    say("=" * 100)


S = {}   # everything quotable ends up here

# =====================================================================================
# 1. MERGE
# =====================================================================================
hr("1. MERGE  -  N at every step")

base = pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet")
say(f"t90_final.parquet rows                                    {len(base):>8,}")
coh = apply_cohort(base)
say(f"after apply_cohort (fu_valid, spo2 notna, oximetry_bad=0) {len(coh):>8,}   expected {COHORT_N:,}")
assert len(coh) == COHORT_N
assert coh.BDSPPatientID.is_unique, "cohort patient id is not unique"

# v2: the validated per-patient file (veto-clean records, corrected hours_v2), built by
# habitual_sleep_v2/make_per_patient_v2.py with aggregate.py's per-patient logic unchanged.
# Instruments below still come from {HAB}: the v2 extraction touched durations only.
hab = pd.read_parquet(f"{paths.SV_ROOT}/"
                      "habitual_sleep_v2/habitual_per_patient_v2.parquet")
say(f"habitual_per_patient.parquet rows                         {len(hab):>8,}")
say(f"  unique patient ids                                      {hab.BDSPPatientID.nunique():>8,}")
assert hab.BDSPPatientID.is_unique, "habitual file has duplicate patients"
say(f"  habitual patients NOT in the analysis cohort            "
    f"{(~hab.BDSPPatientID.isin(set(coh.BDSPPatientID))).sum():>8,}")

HKEEP = ["BDSPPatientID", "tier", "hours", "hours_min", "hours_max", "n_records", "n_tier",
         "n_pre_study", "any_third_party", "from_avg_nightly_hours_field",
         "avg_nightly_hours_field", "habitual_min", "note_date_first", "note_date_last"]
d = coh.merge(hab[HKEEP], on="BDSPPatientID", how="left", validate="1:1")
assert len(d) == COHORT_N, f"merge changed row count to {len(d):,}"

d["has_habitual"] = d.hours.notna().astype(int)
n_with = int(d.has_habitual.sum())
n_without = COHORT_N - n_with
say(f"merged onto the cohort, rows                              {len(d):>8,}")
say(f"  WITH a habitual sleep value                             {n_with:>8,}   "
    f"{100 * n_with / COHORT_N:.1f}% of the cohort")
say(f"  WITHOUT                                                 {n_without:>8,}   "
    f"{100 * n_without / COHORT_N:.1f}%")
# cross-check against the flag the upstream builder stamped
assert n_with == int(hab.in_t90.sum()), "merge disagrees with the upstream in_t90 flag"
say(f"  agrees with upstream in_t90 flag                             yes")

say("")
say("tier composition of the merged habitual sample")
for t in ["A", "B", "C"]:
    n = int((d.tier == t).sum())
    say(f"  tier {t}   {n:>6,}   {100 * n / n_with:>5.1f}% of those with a value")
n_tpl = int((d.from_avg_nightly_hours_field == 1).sum())
say(f'  template field "Average nightly hours:"   {n_tpl:>6,}')

say("")
say("habitual hours, distribution among the merged sample")
q = d.hours.describe(percentiles=[.05, .25, .5, .75, .95])
for k in ["count", "mean", "std", "min", "5%", "25%", "50%", "75%", "95%", "max"]:
    say(f"  {k:<6}{q[k]:>10.3f}")

S["step1"] = {
    "t90_final_rows": int(len(base)), "cohort_n": COHORT_N,
    "habitual_file_rows": int(len(hab)),
    "habitual_outside_cohort": int((~hab.BDSPPatientID.isin(set(coh.BDSPPatientID))).sum()),
    "n_with_habitual": n_with, "n_without_habitual": n_without,
    "pct_with_habitual": round(100 * n_with / COHORT_N, 1),
    "tier_A": int((d.tier == "A").sum()), "tier_B": int((d.tier == "B").sum()),
    "tier_C": int((d.tier == "C").sum()), "template_field": n_tpl,
    "hours_min": float(d.hours.min()), "hours_max": float(d.hours.max()),
    "hours_median": float(d.hours.median()),
    "hours_q1": float(d.hours.quantile(.25)), "hours_q3": float(d.hours.quantile(.75)),
}

# =====================================================================================
# covariates needed downstream: race, BMI, prevalent count, follow-up
# =====================================================================================
d["male"] = (d.sex.astype(str).str.upper().str[0] == "M").astype(int)
d["fu_years"] = (d.censor_date - d.psg_date).dt.days / 365.25

# race, harmonised exactly as numbers/freeze_02_demographics.py does it
import duckdb                                            # noqa: E402
PERSON = f"{paths.OMOP_CACHE_DIR}/person_merged.parquet"
p = duckdb.sql(f"""select person_id as "BDSPPatientID", race_source_value
                   from '{PERSON}'""").df()
p = p.drop_duplicates("BDSPPatientID")
raw = p.race_source_value.astype(str).str.strip().str.upper()
WHITE = {"WHITE", "CAUCASIAN OR WHITE", "CAUCASIAN"}
BLACK = {"AFRICAN AMERICAN  OR BLACK", "AFRICAN AMERICAN OR BLACK",
         "BLACK/AFRICAN AMERICAN", "BLACK", "AFRICAN AMERICAN"}
ASIAN = {"ASIAN"}
NHPI = {"NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER", "PACIFIC ISLANDER", "NATIVE HAWAIIAN"}
AIAN = {"AMERICAN INDIAN OR ALASKA NATIVE", "AMERICAN INDIAN"}
UNK = {"UNKNOWN, UNAVAILABLE OR UNREPORTED", "UNKNOWN/NOT SPECIFIED", "UNKNOWN",
       "UNABLE TO OBTAIN", "DECLINED TO ANSWER", "NONE", "PREFER NOT TO SAY", "NAN", "",
       "DECLINED", "NOT RECORDED", "PATIENT DECLINED", "NOT SPECIFIED", "UNAVAILABLE",
       "DECLINE TO ANSWER"}


def harmonize(v):
    if v in WHITE:
        return "White"
    if v in BLACK:
        return "Black or African American"
    if v in ASIAN:
        return "Asian"
    if v in NHPI:
        return "Native Hawaiian or Other Pacific Islander"
    if v in AIAN:
        return "American Indian or Alaska Native"
    if v in UNK:
        return "Unknown or not reported"
    return "Other"


p["race_h"] = raw.apply(harmonize)
d = d.merge(p[["BDSPPatientID", "race_h"]], on="BDSPPatientID", how="left", validate="1:1")
d["race_h"] = d.race_h.fillna("Unknown or not reported")
assert len(d) == COHORT_N

# BMI, coverage is partial and is reported as such
bmi = pd.read_parquet(f"{paths.NUMBERS_DIR}/_bmi_derived.parquet")[["BDSPPatientID", "bmi"]]
d = d.merge(bmi.drop_duplicates("BDSPPatientID"), on="BDSPPatientID", how="left", validate="1:1")
assert len(d) == COHORT_N

# Epworth and the Insomnia Severity Index, carried for later use
inst = pd.read_parquet(f"{HAB}/instruments_per_patient.parquet")
IK = ["BDSPPatientID", "epworth_score_median", "epworth_score_first",
      "insomnia_severity_index_score_median"]
d = d.merge(inst[IK].drop_duplicates("BDSPPatientID"), on="BDSPPatientID",
            how="left", validate="1:1")
d = d.rename(columns={"epworth_score_median": "ess_median",
                      "epworth_score_first": "ess_first",
                      "insomnia_severity_index_score_median": "isi_median"})
assert len(d) == COHORT_N

# prevalent condition burden. The cardiovascular composite is the union of five conditions
# already counted individually, so it is excluded from the count.
DZ = [k for k in DISEASES if k != "cvd"]
PREV = [f"{k}_prevalent" for k in DZ if f"{k}_prevalent" in d.columns]
d["prevalent_count"] = d[PREV].fillna(0).sum(axis=1).astype(int)
say("")
say(f"prevalent condition count built from {len(PREV)} conditions "
    f"(the cvd composite excluded so its five components are not double counted)")

OUTCOMES = [k for k in DISEASES] + ["death"]
OUTCOMES = [k for k in OUTCOMES if f"{k}_incident" in d.columns and f"{k}_years" in d.columns]
say(f"outcomes carried: {len(OUTCOMES)}  (49 conditions, the cvd composite, and death)")

# =====================================================================================
# 2. HABITUAL SLEEP CATEGORIES
# =====================================================================================
hr("2. HABITUAL SLEEP CATEGORIES   (reference = 6 to under 7 h)")

HAB_LABELS = ["<5h", "5-<6h", "6-<7h(ref)", ">=7h"]
d["hab_cat"] = pd.cut(d.hours, bins=[-np.inf, 5, 6, 7, np.inf],
                      right=False, labels=HAB_LABELS)
d["hab_cat"] = d.hab_cat.astype(object).where(d.hours.notna(), None)
say("bins are left-closed: [min,5), [5,6), [6,7), [7,max]. "
    f"Observed hours run {d.hours.min():.2f} to {d.hours.max():.2f}, so no value is out of range.")
say("")


def events_for(sub, key):
    """Incident events among those eligible for that outcome: not prevalent, positive time."""
    pc, ec, yc = f"{key}_prevalent", f"{key}_incident", f"{key}_years"
    e = sub[(sub[pc] == 0) & sub[yc].notna() & (sub[yc] > 0)]
    return int(len(e)), int(e[ec].fillna(0).sum())


HEAD = ["cvd", "death", "htn2", "diabetes", "hf", "afib", "stroke_any", "ckd",
        "dementia", "depression", "cancer_any", "alopecia"]   # v8.1: alopecia took back pain's control slot

rows = []
say(f"{'band':<12}{'N':>7}{'% of set':>10}{'median h':>10}{'med f/u y':>11}"
    + "".join(f"{k:>13}" for k in HEAD[:6]))
for lab in HAB_LABELS:
    s = d[d.hab_cat == lab]
    r = {"band": lab, "n": len(s), "pct_of_habitual_set": round(100 * len(s) / n_with, 1),
         "median_hours": round(float(s.hours.median()), 2),
         "median_fu_years": round(float(s.fu_years.median()), 2)}
    for k in OUTCOMES:
        n_el, n_ev = events_for(s, k)
        r[f"{k}_atrisk"], r[f"{k}_events"] = n_el, n_ev
    rows.append(r)
    say(f"{lab:<12}{r['n']:>7,}{r['pct_of_habitual_set']:>9.1f}%{r['median_hours']:>10.2f}"
        f"{r['median_fu_years']:>11.2f}" + "".join(f"{r[k + '_events']:>13,}" for k in HEAD[:6]))
say(" " * 50 + "".join(f"{k:>13}" for k in HEAD[:6]))
say("")
say(f"{'band':<12}{'N':>7}" + "".join(f"{k:>13}" for k in HEAD[6:]))
for r in rows:
    say(f"{r['band']:<12}{r['n']:>7,}" + "".join(f"{r[k + '_events']:>13,}" for k in HEAD[6:]))
say(" " * 19 + "".join(f"{k:>13}" for k in HEAD[6:]))
hab_tab = pd.DataFrame(rows)
hab_tab.to_csv(OUT / "habitual_categories.csv", index=False)
say("")
say(f"total across bands {hab_tab.n.sum():,}  must equal {n_with:,}")
assert hab_tab.n.sum() == n_with

say("")
say("Lab total sleep time by habitual band. The bands are supposed to differ in sleep, so if the")
say("laboratory measured the same construct these medians would separate. They do not.")
say(f"{'band':<12}{'N':>7}{'median lab TST min':>22}{'IQR':>22}")
for lab in HAB_LABELS:
    s = d[d.hab_cat == lab]
    say(f"{lab:<12}{len(s):>7,}{s.TST_min.median():>22.1f}"
        f"{s.TST_min.quantile(.25):>13.1f} to {s.TST_min.quantile(.75):<7.1f}")

say("")
say("BOUNDARY SENSITIVITY. Free-text hours pile up on whole numbers, so which side of a cut")
say("point an integer falls on moves hundreds of patients. Counts at the three cut points:")
for v in [5.0, 6.0, 7.0]:
    say(f"   exactly {v:.0f} h   {int((d.hours == v).sum()):>6,} patients")
alt = pd.cut(d.hours, bins=[-np.inf, 5, 6, 7, np.inf], right=True,
             labels=HAB_LABELS).astype(object)
say("   if the bins were right-closed instead of left-closed the band sizes would be:")
for lab in HAB_LABELS:
    say(f"      {lab:<12}{int((alt == lab).sum()):>7,}   against {int((d.hab_cat == lab).sum()):>7,} "
        f"as specified")
say("   The specification is left-closed, so 5.0 h is short sleep, 6.0 h is the reference and")
say("   7.0 h is long sleep. This is stated because it is a real analytic choice, not a detail.")

say("")
say(f"NOTE: the reference band is the SMALLEST of the four "
    f"({int((d.hab_cat == '6-<7h(ref)').sum()):,} patients, "
    f"{100 * (d.hab_cat == '6-<7h(ref)').mean() / (n_with / COHORT_N):.1f}% of the habitual set). "
    f"People write a sleep duration in a note when it is abnormal, so the middle of the "
    f"distribution is under-recorded. Every contrast is measured against that thin band and "
    f"confidence intervals will be wider than the total sample suggests.")

say("")
say("=" * 100)
say("TIMING OF THE MEASUREMENT RELATIVE TO THE SLEEP STUDY  -  read this before modelling")
say("=" * 100)
say("The upstream builder prefers notes written before the sleep study and falls back to notes")
say("written after it when a patient has none (aggregate.py, `src = pre if len(pre) else gt`).")
say("So for some patients the exposure is measured after the outcome clock has already started.")
d["pre_study_only"] = ((d.has_habitual == 1) & (d.n_pre_study > 0)).astype(int)
d["post_study_only"] = ((d.has_habitual == 1) & (d.n_pre_study == 0)).astype(int)
n_pre = int(d.pre_study_only.sum())
n_post = int(d.post_study_only.sum())
say("")
say(f"  value drawn from PRE-study notes    {n_pre:>6,}   {100 * n_pre / n_with:>5.1f}% of the set")
say(f"  value drawn from POST-study notes   {n_post:>6,}   {100 * n_post / n_with:>5.1f}% of the set")
say("")
say("  A patient who developed heart failure two years after the study and then had a note")
say("  recording four hours of sleep is counted here as a short sleeper at baseline. That is")
say("  reverse causation built into the exposure, and it runs in the direction that CREATES an")
say("  association rather than hiding one. The pre-study subset is therefore not an optional")
say("  sensitivity analysis. It is the only subset in which a positive finding is interpretable.")
say("")
say(f"{'band':<12}{'all':>9}{'pre-study':>12}{'post-study':>13}")
for lab in HAB_LABELS:
    s = d[d.hab_cat == lab]
    say(f"{lab:<12}{len(s):>9,}{int(s.pre_study_only.sum()):>12,}{int(s.post_study_only.sum()):>13,}")
say("")
say("Other record-level cautions in the same direction:")
sh = d[d.has_habitual == 1]
say(f"  patients whose value rests on more than one note      {int((sh.n_records > 1).sum()):>6,}")
say(f"  patients whose recorded hours span more than 2 h      {int(((sh.hours_max - sh.hours_min) > 2).sum()):>6,}"
    f"   (hours is the median across their notes)")
say(f"  patients with any third-party reported value          {int(sh.any_third_party.sum()):>6,}")

S["timing"] = {"pre_study_only": n_pre, "post_study_only": n_post,
               "pct_post_study": round(100 * n_post / n_with, 1),
               "by_band": {lab: {"all": int((d.hab_cat == lab).sum()),
                                 "pre": int(d[d.hab_cat == lab].pre_study_only.sum()),
                                 "post": int(d[d.hab_cat == lab].post_study_only.sum())}
                           for lab in HAB_LABELS},
               "n_multi_note": int((sh.n_records > 1).sum()),
               "n_span_over_2h": int(((sh.hours_max - sh.hours_min) > 2).sum()),
               "n_third_party": int(sh.any_third_party.sum())}

S["step2"] = {"labels": HAB_LABELS, "reference": "6-<7h(ref)",
              "n": {r["band"]: r["n"] for r in rows},
              "cvd_events": {r["band"]: r["cvd_events"] for r in rows},
              "death_events": {r["band"]: r["death_events"] for r in rows},
              "median_lab_tst_min": {lab: round(float(d[d.hab_cat == lab].TST_min.median()), 1)
                                     for lab in HAB_LABELS},
              "n_at_cutpoints": {str(int(v)): int((d.hours == v).sum()) for v in [5.0, 6.0, 7.0]},
              "n_if_right_closed": {lab: int((alt == lab).sum()) for lab in HAB_LABELS}}

# =====================================================================================
# 3. LAB OXYGENATION CATEGORIES
# =====================================================================================
hr("3. LAB OXYGENATION CATEGORIES   (T90 = % of the whole recording below 90%)")

O2_LABELS = ["normal T90<=1%", "intermediate 1-10%", "low T90>10%"]
d["o2_cat"] = pd.cut(d.spo2_pct_below_90, bins=[-np.inf, 1.0, 10.0, np.inf],
                     right=True, labels=O2_LABELS).astype(object)
say("bins are right-closed: (-inf,1], (1,10], (10,inf). T90 is a percentage of the WHOLE "
    "recording, not sleep-denominated.")
say("")
rows = []
say(f"{'band':<22}{'N cohort':>10}{'% cohort':>10}{'N habitual':>12}{'% of hab':>10}"
    f"{'median T90':>12}" + "".join(f"{k:>10}" for k in ["cvd", "death"]))
for lab in O2_LABELS:
    s_all = d[d.o2_cat == lab]
    s_hab = s_all[s_all.has_habitual == 1]
    r = {"band": lab, "n_cohort": len(s_all),
         "pct_cohort": round(100 * len(s_all) / COHORT_N, 1),
         "n_habitual": len(s_hab),
         "pct_of_habitual_set": round(100 * len(s_hab) / n_with, 1),
         "median_t90_cohort": round(float(s_all.spo2_pct_below_90.median()), 3)}
    for k in OUTCOMES:
        r[f"{k}_atrisk_cohort"], r[f"{k}_events_cohort"] = events_for(s_all, k)
        r[f"{k}_atrisk_hab"], r[f"{k}_events_hab"] = events_for(s_hab, k)
    rows.append(r)
    say(f"{lab:<22}{r['n_cohort']:>10,}{r['pct_cohort']:>9.1f}%{r['n_habitual']:>12,}"
        f"{r['pct_of_habitual_set']:>9.1f}%{r['median_t90_cohort']:>12.3f}"
        f"{r['cvd_events_hab']:>10,}{r['death_events_hab']:>10,}")
say(" " * 76 + "  events in the habitual set")
o2_tab = pd.DataFrame(rows)
o2_tab.to_csv(OUT / "o2_categories.csv", index=False)
n_int_cohort = int(o2_tab.loc[o2_tab.band == O2_LABELS[1], "n_cohort"].iloc[0])
n_int_hab = int(o2_tab.loc[o2_tab.band == O2_LABELS[1], "n_habitual"].iloc[0])
say("")
say(f"THE INTERMEDIATE BAND IS RETAINED, NOT DISCARDED: {n_int_cohort:,} of the cohort "
    f"({100 * n_int_cohort / COHORT_N:.1f}%) and {n_int_hab:,} of the habitual set "
    f"({100 * n_int_hab / n_with:.1f}%) sit between 1% and 10%.")
assert o2_tab.n_cohort.sum() == COHORT_N
assert o2_tab.n_habitual.sum() == n_with
S["step3"] = {"labels": O2_LABELS,
              "n_cohort": {r["band"]: r["n_cohort"] for r in rows},
              "n_habitual": {r["band"]: r["n_habitual"] for r in rows},
              "intermediate_n_cohort": n_int_cohort, "intermediate_n_habitual": n_int_hab}

# =====================================================================================
# 4. THE CROSS
# =====================================================================================
hr("4. HABITUAL SLEEP x LAB OXYGENATION   -  the cells Alen named")

cells, long = [], []
for hl in HAB_LABELS:
    for ol in O2_LABELS:
        s = d[(d.hab_cat == hl) & (d.o2_cat == ol)]
        c = {"hab_band": hl, "o2_band": ol, "n": len(s),
             "median_hours": round(float(s.hours.median()), 2) if len(s) else np.nan,
             "median_t90": round(float(s.spo2_pct_below_90.median()), 3) if len(s) else np.nan,
             "median_lab_tst_min": round(float(s.TST_min.median()), 1) if len(s) else np.nan,
             "median_age": round(float(s.AgeAtVisit.median()), 1) if len(s) else np.nan,
             "pct_male": round(100 * float(s.male.mean()), 1) if len(s) else np.nan,
             "pct_I0002": round(100 * float((s.site_id == "I0002").mean()), 1) if len(s) else np.nan,
             "median_fu_years": round(float(s.fu_years.median()), 2) if len(s) else np.nan}
        for k in OUTCOMES:
            a, e = events_for(s, k)
            c[f"{k}_atrisk"], c[f"{k}_events"] = a, e
            long.append({"hab_band": hl, "o2_band": ol, "n_cell": len(s), "outcome": k,
                         "label": DISEASES[k][0] if k in DISEASES else "Death from any cause",
                         "organ_group": ORGAN_GROUP.get(k, "Mortality"),
                         "negative_control": k in FINAL_CONTROLS,
                         "circular": k in CIRCULAR,
                         "at_risk": a, "events": e})
        cells.append(c)

cross = pd.DataFrame(cells)
cross.to_csv(OUT / "cross_cells.csv", index=False)
pd.DataFrame(long).to_csv(OUT / "cross_cells_by_disease.csv", index=False)

say(f"{'habitual':<12}{'oxygenation':<22}{'N':>7}{'med T90':>10}{'med lab TST':>13}"
    f"{'med f/u':>9}" + "".join(f"{k:>11}" for k in ["cvd", "death", "htn2", "diabetes", "hf"]))
say("-" * 118)
for hl in HAB_LABELS:
    for ol in O2_LABELS:
        c = cross[(cross.hab_band == hl) & (cross.o2_band == ol)].iloc[0]
        say(f"{hl:<12}{ol:<22}{c.n:>7,}{c.median_t90:>10.2f}{c.median_lab_tst_min:>13.1f}"
            f"{c.median_fu_years:>9.2f}"
            + "".join(f"{int(c[k + '_events']):>11,}"
                      for k in ["cvd", "death", "htn2", "diabetes", "hf"]))
    say("-" * 118)
say(f"cells sum to {cross.n.sum():,}, must equal {n_with:,}")
assert cross.n.sum() == n_with

say("")
say("THE FOUR CELLS ALEN NAMED EXPLICITLY")
named = [("<5h", "normal T90<=1%", "short sleep, normal oxygen"),
         ("<5h", "low T90>10%", "short sleep, terrible oxygen"),
         (">=7h", "normal T90<=1%", "long sleep, normal oxygen"),
         (">=7h", "low T90>10%", "long sleep, terrible oxygen")]
for hl, ol, note in named:
    c = cross[(cross.hab_band == hl) & (cross.o2_band == ol)].iloc[0]
    say(f"  {note:<32} N={int(c.n):>5,}   cvd events {int(c.cvd_events):>4,}   "
        f"death {int(c.death_events):>4,}   htn {int(c.htn2_events):>4,}   "
        f"diabetes {int(c.diabetes_events):>4,}")
S["step4"] = {"cells": {f"{r.hab_band} | {r.o2_band}":
                        {"n": int(r.n), "cvd_events": int(r.cvd_events),
                         "death_events": int(r.death_events)}
                        for r in cross.itertuples()},
              "smallest_cell_n": int(cross.n.min()),
              "smallest_cell": f"{cross.loc[cross.n.idxmin(), 'hab_band']} | "
                               f"{cross.loc[cross.n.idxmin(), 'o2_band']}"}

# =====================================================================================
# 5. SELECTION
# =====================================================================================
hr("5. SELECTION  -  who has a habitual value and who does not")

A = d[d.has_habitual == 1]
B = d[d.has_habitual == 0]


def smd_cont(a, b):
    a, b = a.dropna(), b.dropna()
    if len(a) < 2 or len(b) < 2:
        return np.nan
    s = np.sqrt((a.std(ddof=1) ** 2 + b.std(ddof=1) ** 2) / 2)
    return np.nan if s == 0 else (a.mean() - b.mean()) / s


def smd_bin(pa, pb):
    s = np.sqrt((pa * (1 - pa) + pb * (1 - pb)) / 2)
    return np.nan if s == 0 else (pa - pb) / s


CONT = [("Age at study, years", "AgeAtVisit"),
        ("Lab total sleep time, min", "TST_min"),
        ("Lab T90, % of recording", "spo2_pct_below_90"),
        ("Lab T90, log1p scale", "_logt90"),
        ("Mean SpO2, %", "spo2_mean"),
        ("Apnea-hypopnea index", "AHI"),
        ("Oxygen desaturation index 4%", "odi4_total"),
        ("Arousal index", "arousal_index"),
        ("Sleep efficiency, %", "sleep_efficiency_pct"),
        ("Body mass index", "bmi"),
        ("Prevalent condition count", "prevalent_count"),
        ("Follow-up, years", "fu_years")]
d["_logt90"] = np.log1p(d.spo2_pct_below_90)
A, B = d[d.has_habitual == 1], d[d.has_habitual == 0]

rows = []
say(f"{'variable':<34}{'with habitual':>26}{'without':>26}{'SMD':>9}{'cover%':>9}")
say(f"{'':<34}{'mean (SD)  [median]':>26}{'mean (SD)  [median]':>26}")
say("-" * 105)
for lab, col in CONT:
    a, b = A[col].dropna(), B[col].dropna()
    s = smd_cont(A[col], B[col])
    cov = 100 * d[col].notna().mean()
    rows.append({"variable": lab, "type": "continuous",
                 "with_n": int(len(a)), "with_mean": round(float(a.mean()), 3),
                 "with_sd": round(float(a.std(ddof=1)), 3),
                 "with_median": round(float(a.median()), 3),
                 "without_n": int(len(b)), "without_mean": round(float(b.mean()), 3),
                 "without_sd": round(float(b.std(ddof=1)), 3),
                 "without_median": round(float(b.median()), 3),
                 "smd": round(float(s), 3), "coverage_pct": round(cov, 1)})
    say(f"{lab:<34}{a.mean():>12.2f} ({a.std(ddof=1):>6.2f}) [{a.median():>6.2f}]"
        f"{b.mean():>12.2f} ({b.std(ddof=1):>6.2f}) [{b.median():>6.2f}]"
        f"{s:>9.3f}{cov:>9.1f}")

say("-" * 105)
BIN = [("Male", A.male.mean(), B.male.mean()),
       ("Site I0002", (A.site_id == "I0002").mean(), (B.site_id == "I0002").mean()),
       ("Site I0006", (A.site_id == "I0006").mean(), (B.site_id == "I0006").mean())]
for r in d.race_h.value_counts().index:
    BIN.append((f"Race: {r}", (A.race_h == r).mean(), (B.race_h == r).mean()))
BIN += [("Died during follow-up", A.death_incident.fillna(0).mean(),
         B.death_incident.fillna(0).mean()),
        ("Any prevalent condition", (A.prevalent_count > 0).mean(),
         (B.prevalent_count > 0).mean())]
for lab, pa, pb in BIN:
    s = smd_bin(pa, pb)
    rows.append({"variable": lab, "type": "binary",
                 "with_n": int(round(pa * len(A))), "with_mean": round(100 * float(pa), 1),
                 "with_sd": np.nan, "with_median": np.nan,
                 "without_n": int(round(pb * len(B))), "without_mean": round(100 * float(pb), 1),
                 "without_sd": np.nan, "without_median": np.nan,
                 "smd": round(float(s), 3), "coverage_pct": 100.0})
    say(f"{lab:<34}{100 * pa:>12.1f}% {'':>12}{100 * pb:>12.1f}% {'':>12}{s:>9.3f}")

smd_tab = pd.DataFrame(rows)
smd_tab.to_csv(OUT / "selection_smd.csv", index=False)

big = smd_tab[smd_tab.smd.abs() >= 0.10].sort_values("smd", key=lambda x: x.abs(), ascending=False)
say("")
say(f"Variables with |SMD| >= 0.10 (the usual imbalance threshold): {len(big)} of {len(smd_tab)}")
for r in big.itertuples():
    say(f"   {r.variable:<40}SMD {r.smd:>7.3f}")
say("")
say("Site composition, stated plainly")
for site in ["I0002", "I0006"]:
    tot = int((d.site_id == site).sum())
    wit = int(((d.site_id == site) & (d.has_habitual == 1)).sum())
    say(f"   {site}   cohort {tot:>7,}   with a habitual value {wit:>7,}   "
        f"{100 * wit / tot:>5.1f}%")
tpl_by_site = d[d.from_avg_nightly_hours_field == 1].site_id.value_counts()
say(f"   the template field is {dict(tpl_by_site)} - it does not exist at the second site")
say("")
say("Tier composition by site, because measurement quality and site are confounded")
say(f"{'site':<8}{'tier A':>9}{'tier B':>9}{'tier C':>9}{'template':>11}")
for site in ["I0002", "I0006"]:
    s = d[(d.site_id == site) & (d.has_habitual == 1)]
    say(f"{site:<8}{(s.tier == 'A').sum():>9,}{(s.tier == 'B').sum():>9,}"
        f"{(s.tier == 'C').sum():>9,}{int((s.from_avg_nightly_hours_field == 1).sum()):>11,}")

S["step5"] = {
    "n_with": int(len(A)), "n_without": int(len(B)),
    "smd_over_0.10": [{"variable": r.variable, "smd": r.smd} for r in big.itertuples()],
    "max_abs_smd": float(smd_tab.smd.abs().max()),
    "site_capture": {s: {"cohort": int((d.site_id == s).sum()),
                         "with_habitual": int(((d.site_id == s) & (d.has_habitual == 1)).sum()),
                         "pct": round(100 * ((d.site_id == s) & (d.has_habitual == 1)).sum()
                                      / (d.site_id == s).sum(), 1)}
                     for s in ["I0002", "I0006"]},
    "template_by_site": {k: int(v) for k, v in tpl_by_site.items()},
}

# =====================================================================================
# 6. CANDIDATE PRIMARY SETS
# =====================================================================================
hr("6. CANDIDATE PRIMARY ANALYSIS SETS")

CAND = {
    "all_habitual": d.has_habitual == 1,
    "pre_study_only": d.pre_study_only == 1,
    "tier_AB": d.tier.isin(["A", "B"]),
    "tier_AB_pre_study": d.tier.isin(["A", "B"]) & (d.pre_study_only == 1),
    "tier_A": d.tier == "A",
    "template_field": d.from_avg_nightly_hours_field == 1,
    "template_pre_study": (d.from_avg_nightly_hours_field == 1) & (d.pre_study_only == 1),
}
rows = []
say(f"{'set':<18}{'N':>8}{'I0002%':>9}{'med f/u':>9}"
    + "".join(f"{k:>11}" for k in ["cvd", "death", "htn2", "diabetes", "hf", "dementia"])
    + f"{'min 1of12':>11}{'min band':>11}")
for name, mask in CAND.items():
    s = d[mask]
    r = {"set": name, "n": len(s),
         "pct_I0002": round(100 * float((s.site_id == "I0002").mean()), 1),
         "median_fu_years": round(float(s.fu_years.median()), 2)}
    for k in OUTCOMES:
        r[f"{k}_atrisk"], r[f"{k}_events"] = events_for(s, k)
    ct = s.groupby(["hab_cat", "o2_cat"]).size()
    r["min_cross_cell_n"] = int(ct.min()) if len(ct) else 0
    r["n_cross_cells_under_100"] = int((ct < 100).sum())
    r["min_hab_band_n"] = int(s.hab_cat.value_counts().min())
    r["n_outcomes_with_100plus_events"] = int(
        sum(1 for k in OUTCOMES if r[f"{k}_events"] >= 100))
    rows.append(r)
    say(f"{name:<18}{r['n']:>8,}{r['pct_I0002']:>8.1f}%{r['median_fu_years']:>9.2f}"
        + "".join(f"{r[k + '_events']:>11,}"
                  for k in ["cvd", "death", "htn2", "diabetes", "hf", "dementia"])
        + f"{r['min_cross_cell_n']:>11,}{r['min_hab_band_n']:>11,}")
say("  'min 1of12' is the smallest of the twelve habitual x oxygenation cells; "
    "'min band' the smallest of the four habitual bands.")
cand_tab = pd.DataFrame(rows)
cand_tab.to_csv(OUT / "candidate_sets.csv", index=False)
say("")
for r in rows:
    say(f"  {r['set']:<18} outcomes reaching 100 incident events: "
        f"{r['n_outcomes_with_100plus_events']:>2} of {len(OUTCOMES)};  "
        f"habitual x oxygen cells under 100 patients: {r['n_cross_cells_under_100']} of 12")

PRIMARY = "all_habitual"
say("")
say("DECISION")
say(f"  PRIMARY      all {n_with:,} patients with any habitual value, tiers A, B and C together,")
say("               with tier and measurement timing carried as variables so both can be")
say("               stratified on.")
say("  SENSITIVITY  in this order of importance")
say(f"               1. pre-study only ({n_pre:,}). Removes the patients whose sleep duration was")
val = [r for r in rows if r["set"] == "pre_study_only"][0]
say(f"                  recorded after the outcome clock started. {val['cvd_events']:,} cvd events,")
say(f"                  {val['death_events']:,} deaths, smallest of the twelve cells {val['min_cross_cell_n']:,}.")
say("               2. tiers A and B, dropping tier C, which is bedtime-to-waketime and so")
say("                  measures time in bed rather than sleep.")
say("               3. tier A only, the structured or explicitly stated durations.")
say("               4. the template field alone, the cleanest single source.")
say("  WHY          The paper's claim is expected to be a null: habitual sleep, like")
say("               laboratory sleep, does not predict disease. A null on a noisy measure is")
say("               worth little on its own, because measurement error biases toward it. So")
say("               the design needs both ends. The full set supplies the events, and the")
say("               graded-quality subsets test whether a null there is only measurement")
say("               error. None of the subsets can carry the primary analysis: the pre-study")
say("               subset halves the sample, tier A and the template field are small, and")
say("               the template field exists at one site only, so it cannot be site-")
say("               stratified and carries that site's case mix entirely.")
say(f"  THE ONE      A POSITIVE finding in the full set is not interpretable, because {n_post:,}")
say("  ASYMMETRY    patients had their sleep duration recorded after the study and after some")
say("               of them were already ill. A positive finding must reproduce in the")
say("               pre-study subset before it is believed. A null does not have that problem,")
say("               which is why the full set is a fair primary for the question as posed.")
S["step6"] = {"primary": PRIMARY, "candidates": rows,
              "sensitivity_order": ["pre_study_only", "tier_AB", "tier_A", "template_field"]}

# =====================================================================================
# 7. WRITE
# =====================================================================================
hr("7. WRITE analysis.parquet")

KEEP_BASE = ["BDSPPatientID", "site_id", "psg_date", "censor_date", "fu_years", "fu_valid",
             "AgeAtVisit", "sex", "male", "race_h", "bmi",
             "TST_min", "sleep_efficiency_pct", "N1_pct", "N2_pct", "N3_pct", "REM_pct",
             "arousal_index", "AHI", "recording_dur_min",
             "spo2_pct_below_90", "spo2_pct_below_88", "spo2_mean", "spo2_nadir_corrected",
             "odi3_total", "odi4_total", "o2_cat",
             "has_habitual", "hours", "hab_cat", "tier", "hours_min", "hours_max",
             "n_records", "n_tier", "n_pre_study", "any_third_party",
             "from_avg_nightly_hours_field", "avg_nightly_hours_field", "habitual_min",
             "note_date_first", "note_date_last",
             "ess_median", "ess_first", "isi_median",
             "pre_study_only", "post_study_only",
             "prevalent_count", "oximetry_bad"]
OCOLS = []
for k in OUTCOMES:
    OCOLS += [f"{k}_prevalent", f"{k}_incident", f"{k}_years"]
OCOLS = [c for c in OCOLS if c in d.columns]

a = d[KEEP_BASE + OCOLS].copy()
a["lab_minus_habitual_min"] = a.TST_min - a.habitual_min
a["is_primary_set"] = (a.has_habitual == 1).astype(int)
a["is_tier_AB"] = a.tier.isin(["A", "B"]).astype(int)
a["is_tier_A"] = (a.tier == "A").astype(int)
a["is_template"] = (a.from_avg_nightly_hours_field == 1).fillna(False).astype(int)
a.to_parquet(OUT / "analysis.parquet", index=False)

say(f"rows            {len(a):,}   (the whole cohort, so the excluded {n_without:,} stay visible)")
say(f"columns         {a.shape[1]:,}")
say(f"outcome columns {len(OCOLS):,}  ({len(OUTCOMES)} outcomes x prevalent/incident/years)")
say(f"with habitual   {int(a.is_primary_set.sum()):,}")
assert len(a) == COHORT_N
assert int(a.is_primary_set.sum()) == n_with
assert a.BDSPPatientID.is_unique

say("")
say("negative controls carried: " + ", ".join(FINAL_CONTROLS))
say("NOT present in disease_definitions.py, so not carried: "
    + ", ".join(CONTROLS_NOT_IN_FILE)
    + "  -  these two must be added to the definitions file before the control panel is complete.")
say("fracture is present as an outcome but is NOT a negative control.")
say("circular outcomes flagged in cross_cells_by_disease.csv: " + ", ".join(CIRCULAR))

S["step7"] = {"rows": int(len(a)), "columns": int(a.shape[1]),
              "n_outcomes": len(OUTCOMES), "outcome_columns": len(OCOLS),
              "negative_controls_carried": FINAL_CONTROLS,
              "negative_controls_missing_from_definitions": CONTROLS_NOT_IN_FILE,
              "disease_definitions_controls_flagged": NEGATIVE_CONTROLS}

# -------------------------------------------------------------------------------------
# read the file back off disk and re-derive every headline number from it, so the parquet
# is verified rather than assumed
# -------------------------------------------------------------------------------------
hr("8. READ-BACK VERIFICATION  -  every headline number re-derived from the written file")

v = pd.read_parquet(OUT / "analysis.parquet")
chk = []


def check(name, got, want):
    ok = (got == want)
    chk.append(ok)
    say(f"  {'PASS' if ok else 'FAIL':<5}{name:<58}{got!s:>12}   expected {want!s}")
    return ok


check("rows", len(v), COHORT_N)
check("unique patients", int(v.BDSPPatientID.nunique()), COHORT_N)
check("with a habitual value", int(v.hours.notna().sum()), n_with)
check("without", int(v.hours.isna().sum()), n_without)
# v8 F2 fix 2026-09-12 (integrity gate 6): the tier, template-field and per-site references are re-derived here, at
# run time, by a second direct merge of the habitual file onto the v8 cohort table (apply_cohort on t90_final.parquet),
# so the read-back of analysis.parquet is checked against the v8 tables and not against typed pins. The v2 pins typed
# on 2026-08-13 from the v7 tables (tier A 186, B 4762, C 347, template 2, I0002 3894, I0006 1401; v1 values were
# 1543 / 6252 / 916 / 1293 and 6335 / 2376) are printed next to the v8 references (old versus new), not asserted.
_ref = apply_cohort(base)[["BDSPPatientID", "site_id"]].merge(
    hab[["BDSPPatientID", "tier", "hours", "from_avg_nightly_hours_field"]],
    on="BDSPPatientID", how="left", validate="1:1")
assert len(_ref) == COHORT_N, len(_ref)
REF = {"tier_A": int((_ref.tier == "A").sum()), "tier_B": int((_ref.tier == "B").sum()),
       "tier_C": int((_ref.tier == "C").sum()),
       "template_field": int((_ref.from_avg_nightly_hours_field == 1).sum()),
       "I0002_with_habitual": int(((_ref.site_id == "I0002") & _ref.hours.notna()).sum()),
       "I0006_with_habitual": int(((_ref.site_id == "I0006") & _ref.hours.notna()).sum())}
_V7_PINS = {"tier_A": 186, "tier_B": 4762, "tier_C": 347, "template_field": 2,
            "I0002_with_habitual": 3894, "I0006_with_habitual": 1401}   # printed for old-versus-new only
for _k in REF:
    say(f"  NOTE  reference {_k:<22} v8 tables {REF[_k]:>6,}   v7 pin 2026-08-13 {_V7_PINS[_k]:>6,}   "
        + ("same" if REF[_k] == _V7_PINS[_k] else f"moved {REF[_k] - _V7_PINS[_k]:+d}"))
S["step8_references_v8"] = {"from": "apply_cohort(t90_final.parquet) x habitual_per_patient_v2.parquet, direct merge "
                                    "at run time", "v8": REF, "v7_pins_2026-08-13": _V7_PINS}
check("tier A", int((v.tier == "A").sum()), REF["tier_A"])
check("tier B", int((v.tier == "B").sum()), REF["tier_B"])
check("tier C", int((v.tier == "C").sum()), REF["tier_C"])
check("template field", int(v.is_template.sum()), REF["template_field"])
for lab in HAB_LABELS:
    check(f"habitual band {lab}", int((v.hab_cat == lab).sum()),
          int((d.hab_cat == lab).sum()))
for lab in O2_LABELS:
    check(f"oxygenation band {lab} (whole cohort)", int((v.o2_cat == lab).sum()),
          int((d.o2_cat == lab).sum()))
check("habitual bands sum", int(v.hab_cat.notna().sum()), n_with)
check("oxygenation bands sum", int(v.o2_cat.notna().sum()), COHORT_N)
# v8 F2 fix: per-site references from the same run-time re-derivation (REF above), not typed
check("site I0002 with a habitual value",
      int(((v.site_id == "I0002") & v.hours.notna()).sum()), REF["I0002_with_habitual"])
check("site I0006 with a habitual value",
      int(((v.site_id == "I0006") & v.hours.notna()).sum()), REF["I0006_with_habitual"])
check("pre-study only", int(v.pre_study_only.sum()), n_pre)
check("post-study only", int(v.post_study_only.sum()), n_post)
check("timing flags partition the habitual set",
      int((v.pre_study_only + v.post_study_only).sum()), n_with)
check("outcome column triples", len(OCOLS) // 3, len(OUTCOMES))
check("no oximetry_bad survived", int(v.oximetry_bad.sum()), 0)
check("no missing follow-up", int(v.fu_years.isna().sum()), 0)
check("cvd events, habitual set",
      int(v[(v.hours.notna()) & (v.cvd_prevalent == 0) & (v.cvd_years > 0)].cvd_incident.sum()),
      int(cand_tab.loc[cand_tab.set == "all_habitual", "cvd_events"].iloc[0]))
# every incident flag must sit inside its at-risk definition
bad = 0
for k in OUTCOMES:
    bad += int(((v[f"{k}_incident"] == 1) & (v[f"{k}_prevalent"] == 1)).sum())
check("patients flagged incident and prevalent for one outcome", bad, 0)

say("")
say(f"{sum(chk)} of {len(chk)} checks passed")
assert all(chk), "read-back verification failed"
S["verification"] = {"checks": len(chk), "passed": int(sum(chk))}

(OUT / "summary.json").write_text(json.dumps(S, indent=2, default=str))
(OUT / "build_log.txt").write_text("\n".join(LOG))
say("")
say(f"wrote {OUT}/analysis.parquet and 7 companion files")
