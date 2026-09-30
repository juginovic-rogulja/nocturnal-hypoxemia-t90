#!/usr/bin/env python3
"""
Build the analysis dataset for the stage-specific T90 question.

Question. Is sleep T90 stage-specific when it predicts long-term outcomes, or is one stage
carrying it? Wanted in the general cohort and, for the first time, in a disease-free subgroup.

This script only assembles and describes the dataset. Nothing is modelled here.

What it produces
  analysis.parquet                  primary, one row per patient, strict PAP-off definition
  analysis_papoff_portion.parquet   sensitivity, the looser PAP-off-portion definition that the
                                    2026-07-29 deliverable used, kept so the two can be compared
                                    rather than the choice being buried

PAP-off, two definitions, stated because they differ by thousands of patients
  STRICT   cpap_frac == 0 and has_cpap_chan is true. A PAP channel was recorded and carried no
           pressure, so the whole night is verified PAP-off and every per-stage window is drawn
           from the same untreated recording. This is the definition the task specifies and the
           one analysis.parquet uses.
  PORTION  pap_off_min > 60 and tst_papoff_min > 30, which additionally admits studies with no
           PAP channel at all and split nights where only the PAP-off part was scored. This is
           what the earlier T90_by_state_vs_AHI_all_diseases.csv table was built on.

Traps this script respects (they are documented, not hypothetical)
  - dedupe on BIDSFolder + SessionID before anything else, the shards overlap
  - PAP-off first, then first study per patient
  - the cohort filter is fu_valid == 1 and spo2_pct_below_90 not null and oximetry_bad == 0,
    never fu_valid alone
  - every stage is gated at >= 30 minutes of that stage and the surviving N is reported, the
    N differs a lot by stage and that is itself a result
  - NREM and N2 correlate 0.963, so they must never enter one model. Not enforceable here,
    flagged for the modelling step.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths

import glob
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------------------- paths
SHARD_GLOB = (f"{paths.T90_ROOT}/BDSP_New_Ideas_2026-07-26"
              "/pilots/t90_by_stage/out_*.csv").replace(f"{paths.T90_ROOT}/pilots/t90_by_stage", f"{paths.V8_ROOT}/X4_groupP/t90_by_stage")   # v8.2 re-extracted shards
T90_ROOT = paths.T90_ROOT
COHORT_PARQUET = f"{paths.TABLES_DIR}/t90_final.parquet"
OUT_DIR = Path(f"{paths.SV_ROOT}/stage_specific")

sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
from disease_definitions import DISEASES, NEGATIVE_CONTROLS, ORGAN_GROUP, CIRCULAR  # noqa: E402

STAGES = ["wake", "sleep", "nrem", "n1", "n2", "n3", "rem"]
STAGE_GATE_MIN = 30.0

OUT_DIR.mkdir(parents=True, exist_ok=True)


def hdr(t):
    print("\n" + "=" * 86)
    print(t)
    print("=" * 86)


# ================================================================== 1. pool and dedupe shards
hdr("STEP 1  pool the per-stage shards")

files = sorted(glob.glob(SHARD_GLOB))
raw = pd.concat([pd.read_csv(f, low_memory=False).assign(_shard=Path(f).name) for f in files],
                ignore_index=True)
print(f"shards read                                : {len(files)}")
print(f"rows pooled                                : {len(raw):,}")
print(f"unique patients in the pool                : {raw.BDSPPatientID.nunique():,}")

d = raw.drop_duplicates(["BIDSFolder", "SessionID"], keep="first").copy()
print(f"rows after dedupe on BIDSFolder+SessionID  : {len(d):,}  "
      f"(removed {len(raw) - len(d):,})")

n_ok = int((d.status == "ok").sum())
print(f"\nextraction status")
for k, v in d.status.value_counts(dropna=False).items():
    print(f"  {str(k):22s} {v:7,}")

d = d[d.status == "ok"].copy()
print(f"rows with status ok                        : {len(d):,}")

# quality flag from the pilot, carried and reported, not used as a filter here
d["spo2_suspect"] = (d.t90_all > 90) | (d.spo2_valid_frac < 0.5)

# ------------------------------------------------------------------- the two PAP-off definitions
d["papoff_strict"] = (d.cpap_frac.fillna(-1) == 0) & (d.has_cpap_chan == True)   # noqa: E712
d["papoff_portion"] = (d.pap_off_min > 60) & (d.tst_papoff_min > 30)

print(f"\nPAP channel present                        : "
      f"{int((d.has_cpap_chan == True).sum()):,} of {len(d):,}")
print(f"cpap_frac exactly 0                        : {int((d.cpap_frac.fillna(-1) == 0).sum()):,}")
print(f"cpap_frac missing (no PAP channel)         : {int(d.cpap_frac.isna().sum()):,}")
split = (d.cpap_frac.fillna(0) > 0) & (d.cpap_frac.fillna(0) <= 0.5)
print(f"split nights (0 < cpap_frac <= 0.5)        : {int(split.sum()):,}")
print(f"rows PAP-off STRICT                        : {int(d.papoff_strict.sum()):,}")
print(f"rows PAP-off PORTION                       : {int(d.papoff_portion.sum()):,}")


def first_study(frame):
    """First study per patient, after the PAP-off filter, ordered by SessionID."""
    return (frame.sort_values(["BDSPPatientID", "SessionID"])
                 .drop_duplicates("BDSPPatientID", keep="first"))


psg_strict = first_study(d[d.papoff_strict])
psg_portion = first_study(d[d.papoff_portion])
print(f"\nfirst study per patient, STRICT            : {len(psg_strict):,}")
print(f"  of which the kept row is not SessionID 1 : {int((psg_strict.SessionID != 1).sum()):,}")
print(f"first study per patient, PORTION           : {len(psg_portion):,}")
print(f"  of which the kept row is not SessionID 1 : {int((psg_portion.SessionID != 1).sum()):,}")
print("\nsite before the outcome merge, STRICT")
for k, v in psg_strict.SiteID.value_counts().items():
    print(f"  {k:8s} {v:7,}")


# ================================================================ 2. merge to outcomes + cohort
hdr("STEP 2  merge to the frozen outcome table and apply the cohort filter")

b = pd.read_parquet(COHORT_PARQUET)
print(f"t90_final.parquet                          : {b.shape[0]:,} x {b.shape[1]}")
print(f"  fu_valid == 1                            : {int((b.fu_valid == 1).sum()):,}")
print(f"  + spo2_pct_below_90 not null             : "
      f"{int(((b.fu_valid == 1) & b.spo2_pct_below_90.notna()).sum()):,}")
coh = b[(b.fu_valid == 1) & b.spo2_pct_below_90.notna() & (b.oximetry_bad == 0)].copy()
print(f"  + oximetry_bad == 0  (analysis cohort)    : {len(coh):,}   [cohort_spec COHORT_N 19173]")
print("\ncohort by site")
for k, v in coh.site_id.value_counts().items():
    print(f"  {k:8s} {v:7,}")
print("sites with no follow-up in the frozen file, so unreachable here: "
      + ", ".join(sorted(set(b.site_id.unique()) - set(coh.site_id.unique()))))

OUTCOME_KEYS = sorted({c[:-len("_prevalent")] for c in b.columns if c.endswith("_prevalent")})
OUTCOME_COLS = [f"{k}_{s}" for k in OUTCOME_KEYS for s in ("prevalent", "incident", "years")]
OUTCOME_COLS = [c for c in OUTCOME_COLS if c in b.columns]
COV_COLS = ["site_id", "AgeAtVisit", "sex", "AHI", "TST_min", "spo2_pct_below_90", "spo2_mean",
            "spo2_pct_below_88", "odi3_total", "odi4_total", "spo2_nadir_corrected",
            "recording_dur_min", "sleep_efficiency_pct", "N1_pct", "N2_pct", "N3_pct", "REM_pct",
            "arousal_index", "oximetry_bad", "fu_valid"]
print(f"\noutcome keys carried                       : {len(OUTCOME_KEYS)}  "
      f"({len(OUTCOME_COLS)} columns)")

PSG_KEEP = (["BDSPPatientID", "SiteID", "BIDSFolder", "SessionID", "AgeAtVisit", "SexDSC",
             "status", "spo2_valid_frac", "spo2_suspect", "has_cpap_chan", "cpap_frac",
             "pap_off_min", "papoff_strict", "papoff_portion",
             "rec_min", "tst_min", "staged_min", "sleep_onset_min", "tst_papoff_min",
             "t90_all", "all_mean", "all_nadir"]
            + [f"t90_{s}" for s in STAGES]
            + [f"min_{s}" for s in STAGES]
            + [f"mean_{s}" for s in STAGES]
            + [f"nadir_{s}" for s in STAGES]
            + ["t90_wakepre", "min_wakepre", "t90_wakepost", "min_wakepost"])


def merge(psg, label):
    p = psg[[c for c in PSG_KEEP if c in psg.columns]].rename(
        columns={"AgeAtVisit": "AgeAtVisit_psg", "SexDSC": "SexDSC_psg"})
    m = p.merge(coh[["BDSPPatientID"] + COV_COLS + OUTCOME_COLS], on="BDSPPatientID", how="inner")
    print(f"\n{label}")
    print(f"  first PAP-off studies offered to merge   : {len(psg):,}")
    print(f"  matched into the analysis cohort         : {len(m):,}  "
          f"(lost {len(psg) - len(m):,})")
    print(f"  duplicate patients after merge           : {int(m.BDSPPatientID.duplicated().sum())}")
    for k, v in m.site_id.value_counts().items():
        print(f"    {k:8s} {v:7,}")
    return m


m_strict = merge(psg_strict, "STRICT  cpap_frac == 0 and a PAP channel present")
m_portion = merge(psg_portion, "PORTION pap_off_min > 60 and tst_papoff_min > 30")

# reconciliation against the table already on disk
hdr("STEP 2b  reconciliation against the 2026-07-29 deliverable")
prior = (f"{paths.T90_ROOT}/BDSP_New_Ideas_2026-07-26"
         "/deliverable/T90_by_state_vs_AHI_all_diseases.csv")
if Path(prior).exists():
    pr = pd.read_csv(prior)
    print(f"prior table: {pr.shape[0]} diseases x {pr.shape[1]} columns, recorded n = 13,066")
rec = m_portion[~m_portion.spo2_suspect]
print(f"PORTION set with the pilot spo2_suspect rule applied : {len(rec):,}")
b_nooxqc = b[(b.fu_valid == 1) & b.spo2_pct_below_90.notna()]
rec2 = psg_portion[(~psg_portion.spo2_suspect)
                   & psg_portion.BDSPPatientID.isin(b_nooxqc.BDSPPatientID)]
print(f"  same, before the 2026-08-04 oximetry_bad QC        : {len(rec2):,}")
print("The prior 13,066 sits inside this band. It predates the oximetry_bad QC and used the")
print("PAP-off PORTION rule, not the strict one. The strict rule is the material difference.")

# ============================================================== 3. per-stage minutes and gates
hdr("STEP 3  build the analysis frame")


def finish(m):
    out = m.copy()
    out["male"] = (out.sex.astype(str).str.upper().str[0] == "M").astype(int)
    for s in STAGES:
        out[f"gate_{s}"] = out[f"min_{s}"].fillna(0) >= STAGE_GATE_MIN
        out[f"t90_{s}_g30"] = out[f"t90_{s}"].where(out[f"gate_{s}"])
    out["gate_all"] = out["rec_min"].fillna(0) >= STAGE_GATE_MIN
    out["t90_all_g30"] = out["t90_all"].where(out["gate_all"])
    # Time-zero check. The outcome clock in t90_final starts at the study the paper indexed.
    # If the first PAP-off study we selected is a later night, the two are not the same recording
    # and the exposure is measured after time zero. Age is the only shared field that can detect
    # it, so it is used as the check and the flag is carried, not applied.
    out["age_gap_yr"] = out.AgeAtVisit_psg - out.AgeAtVisit
    out["index_study_match"] = out.age_gap_yr.abs().le(1) | out.age_gap_yr.isna()
    return out


m_strict = finish(m_strict)
m_portion = finish(m_portion)

print("Does the study we picked look like the study the paper indexed? Age is the shared field.")
for lab, out in (("STRICT", m_strict), ("PORTION", m_portion)):
    g = out.age_gap_yr
    print(f"  {lab:8s} n {len(out):6,}   same age {int((g == 0).sum()):6,}   "
          f"gap > 1 y {int((g.abs() > 1).sum()):5,}   gap > 2 y {int((g.abs() > 2).sum()):5,}   "
          f"index_study_match {int(out.index_study_match.sum()):6,}")
    print(f"           SessionID kept: "
          + ", ".join(f"{k}:{v:,}" for k, v in sorted(out.SessionID.value_counts().items())))
print("Those with a gap have the exposure measured after their follow-up clock started. Small,")
print("but it is a lead-time problem, so the modelling step should refit on index_study_match.")


# =============================================================== 4. the healthy subgroup
hdr("STEP 4  healthy (disease-free at baseline) subgroup")

# Built from ORGAN_GROUP rather than a hand-typed list, so it tracks disease_definitions.py.
EXCLUDE_GROUPS = ["Cardiac", "Respiratory", "Kidney", "Liver", "Metabolic"]
# stroke sits in Neuro/psych in ORGAN_GROUP but is cardiovascular disease, and cvd is the
# cardiovascular composite, so both are added by hand and that is stated.
EXTRA_VASCULAR = ["stroke_any", "cvd"]
# cancer sits in "Other" in ORGAN_GROUP, so the two cancer keys are named directly.
CANCER_KEYS = ["cancer_any", "prostate_ca"]

EXCL_KEYS = sorted({k for k, g in ORGAN_GROUP.items() if g in EXCLUDE_GROUPS}
                   | set(EXTRA_VASCULAR) | set(CANCER_KEYS))
EXCL_KEYS = [k for k in EXCL_KEYS if f"{k}_prevalent" in m_strict.columns]

print("Definition H1, the one used for the `healthy` column in analysis.parquet.")
print("Exclude a patient if ANY of these is prevalent at the sleep study.")
print("Organ groups taken from ORGAN_GROUP: " + ", ".join(EXCLUDE_GROUPS))
print("Added by name: stroke_any and cvd (vascular, filed under Neuro/psych and Composite),")
print("               cancer_any and prostate_ca (filed under Other).")
print(f"{len(EXCL_KEYS)} conditions in the exclusion set:")
for k in EXCL_KEYS:
    print(f"  {k:20s} {ORGAN_GROUP.get(k, '?'):12s} {DISEASES[k][0]}")


def add_healthy(out, label):
    prev = out[[f"{k}_prevalent" for k in EXCL_KEYS]].fillna(0).astype(int)
    out["n_prevalent_excl"] = prev.sum(axis=1)
    out["healthy"] = out["n_prevalent_excl"] == 0
    # two alternatives, reported only, so the cost of each choice is visible
    soft = [k for k in EXCL_KEYS if k not in ("obesity", "dyslipid", "thyroid_dis")]
    out["healthy_no_riskfactors"] = out[[f"{k}_prevalent" for k in soft]].fillna(0).sum(axis=1) == 0
    hard = [k for k in EXCL_KEYS if k not in ("thyroid_dis",)]
    out["healthy_no_thyroid"] = out[[f"{k}_prevalent" for k in hard]].fillna(0).sum(axis=1) == 0
    print(f"\n{label}: n = {len(out):,}")
    print(f"  H1 healthy, all {len(EXCL_KEYS)} conditions absent          : "
          f"{int(out.healthy.sum()):,}  ({100 * out.healthy.mean():.1f}%)")
    print(f"  H2 as H1 but obesity, dyslipidaemia, thyroid allowed : "
          f"{int(out.healthy_no_riskfactors.sum()):,}  "
          f"({100 * out.healthy_no_riskfactors.mean():.1f}%)")
    print(f"  H3 as H1 but thyroid allowed                         : "
          f"{int(out.healthy_no_thyroid.sum()):,}  "
          f"({100 * out.healthy_no_thyroid.mean():.1f}%)")
    return out


m_strict = add_healthy(m_strict, "STRICT set")
m_portion = add_healthy(m_portion, "PORTION set")

print("\nWhat each exclusion costs in the STRICT set, prevalent count and how many patients it")
print("alone removes from the healthy group:")
prev = m_strict[[f"{k}_prevalent" for k in EXCL_KEYS]].fillna(0).astype(int)
for k in EXCL_KEYS:
    c = int(prev[f"{k}_prevalent"].sum())
    only = int(((prev.sum(axis=1) == 1) & (prev[f"{k}_prevalent"] == 1)).sum())
    print(f"  {k:20s} prevalent {c:6,}  ({100 * c / len(m_strict):5.1f}%)   "
          f"sole reason for exclusion {only:5,}")

hh = m_strict[m_strict.healthy]
print(f"\nHealthy subgroup (STRICT, H1) against the general cohort, n = {len(hh):,} of "
      f"{len(m_strict):,}")
print(f"  {'':26s}{'healthy':>12s}{'general':>12s}")
for lab, col, fmt in (("age, median", "AgeAtVisit", "{:.0f}"), ("male %", "male", "{:.1%}"),
                      ("AHI, median", "AHI", "{:.1f}"), ("TST min, median", "TST_min", "{:.0f}"),
                      ("T90 all, median", "t90_all", "{:.3f}"),
                      ("T90 sleep, median", "t90_sleep", "{:.3f}"),
                      ("follow-up y, median", "death_years", "{:.2f}")):
    f = (lambda x: fmt.format(x.mean() if col == "male" else x.median()))
    print(f"  {lab:26s}{f(hh[col]):>12s}{f(m_strict[col]):>12s}")
print(f"  {'deaths':26s}{int(hh.death_incident.sum()):>12,}"
      f"{int(m_strict.death_incident.sum()):>12,}")
print(f"  follow-up IQR healthy {hh.death_years.quantile(.25):.2f} to "
      f"{hh.death_years.quantile(.75):.2f} y, general "
      f"{m_strict.death_years.quantile(.25):.2f} to {m_strict.death_years.quantile(.75):.2f} y")
print("  The healthy group is younger and followed for less time. Disease-free at baseline in an")
print("  EHR partly means fewer encounters, so shorter observed follow-up is expected and the")
print("  healthy analysis is power-limited by follow-up as much as by N.")

# ==================================================== 5. stage gates and medians, both cohorts
hdr("STEP 5  how many patients survive each stage gate, and the median T90 in each stage")


def gate_table(out, label):
    """n_gate is how many reach 30 min of the stage. n_model is how many of those also have a
    computable T90 there, which is what a model would actually fit on. They are not the same."""
    print(f"\n{label}  (n = {len(out):,})")
    print(f"  {'stage':7s} {'n >=30min':>10s} {'% of set':>9s} {'n modelable':>12s} "
          f"{'median min':>11s} {'median T90':>11s} {'IQR T90':>19s} {'% exactly 0':>12s}")
    rows = []
    for s in ["all"] + STAGES:
        mcol = "rec_min" if s == "all" else f"min_{s}"
        gcol = f"t90_{s}_g30"
        n = int(out[f"gate_{s}"].sum())
        v = out[gcol].dropna()
        q1, q3 = (v.quantile([.25, .75]) if len(v) else (np.nan, np.nan))
        med = v.median() if len(v) else np.nan
        zero = 100 * float((v == 0).mean()) if len(v) else np.nan
        print(f"  {s:7s} {n:10,} {100 * n / len(out):8.1f}% {len(v):12,} "
              f"{out[mcol].median():11.1f} {med:11.3f} {q1:8.3f} to {q3:<8.3f} {zero:11.1f}%")
        rows.append((s, n, len(v), med))
    return rows


gate_table(m_strict, "GENERAL cohort, STRICT PAP-off")
gate_table(m_strict[m_strict.healthy], "HEALTHY subgroup H1, STRICT PAP-off")
gate_table(m_portion, "GENERAL cohort, PORTION PAP-off (sensitivity)")
gate_table(m_portion[m_portion.healthy], "HEALTHY subgroup H1, PORTION PAP-off (sensitivity)")

print("\nZero-minute stages in the STRICT general cohort, the reason the stage Ns diverge")
for s in STAGES:
    z = int((m_strict[f"min_{s}"].fillna(0) == 0).sum())
    print(f"  {s:6s} zero minutes {z:6,}   under 30 min "
          f"{int((~m_strict[f'gate_{s}']).sum()):6,}")

print("\nSpearman correlation of the gated stage exposures, STRICT general cohort.")
cc = m_strict[[f"t90_{s}_g30" for s in ["all"] + STAGES]].corr(method="spearman")
cc.index = ["all"] + STAGES
cc.columns = ["all"] + STAGES
print(cc.round(3).to_string())
print(f"\nNREM vs N2 = {cc.loc['nrem', 'n2']:.3f}. They must never enter one model.")

# ===================================================================== 6. event counts
hdr("STEP 6  event counts available for modelling")


def events(out, label):
    print(f"\n{label} (n = {len(out):,})")
    rows = []
    for k in OUTCOME_KEYS:
        f = out[(out[f"{k}_prevalent"] == 0) & out[f"{k}_years"].notna() & (out[f"{k}_years"] > 0)]
        rows.append((k, DISEASES.get(k, (k,))[0], len(f), int(f[f"{k}_incident"].sum())))
    r = pd.DataFrame(rows, columns=["key", "label", "at_risk", "events"]).sort_values(
        "events", ascending=False)
    print(f"  outcomes with >= 60 events : {int((r.events >= 60).sum())} of {len(r)}")
    print(f"  outcomes with >= 150 events: {int((r.events >= 150).sum())} of {len(r)}")
    print("  top 12 by events: " + ", ".join(f"{a} {b:,}" for a, b in
                                             zip(r.label.head(12), r.events.head(12))))
    ncs = r[r.key.isin(NEGATIVE_CONTROLS)]
    print("  negative controls: " + ", ".join(f"{a} {b:,}" for a, b in
                                              zip(ncs.label, ncs.events)))
    return r


ev_gen = events(m_strict, "GENERAL cohort, STRICT")
ev_hlt = events(m_strict[m_strict.healthy], "HEALTHY subgroup H1, STRICT")
print(f"\nnegative controls in the panel : {NEGATIVE_CONTROLS}")
print(f"circular, excluded from ranking: {CIRCULAR}")

# The REM gate is the binding constraint. Show what it costs per outcome.
print("\nEvents surviving the REM >= 30 min gate, healthy subgroup H1")
hr = m_strict[m_strict.healthy & m_strict.gate_rem]
print(f"  healthy and REM >= 30 min: n = {len(hr):,}")
r = []
for k in OUTCOME_KEYS:
    f = hr[(hr[f"{k}_prevalent"] == 0) & hr[f"{k}_years"].notna() & (hr[f"{k}_years"] > 0)]
    r.append((DISEASES.get(k, (k,))[0], int(f[f"{k}_incident"].sum())))
r = pd.DataFrame(r, columns=["label", "events"]).sort_values("events", ascending=False)
print(f"  outcomes with >= 60 events: {int((r.events >= 60).sum())} of {len(r)}")
print("  top 10: " + ", ".join(f"{a} {b:,}" for a, b in zip(r.label.head(10), r.events.head(10))))

# ===================================================================== 7. write
hdr("STEP 7  write")

for frame, name in ((m_strict, "analysis.parquet"),
                    (m_portion, "analysis_papoff_portion.parquet")):
    p = OUT_DIR / name
    frame.to_parquet(p, index=False)
    print(f"wrote {p}   rows {len(frame):,}  cols {frame.shape[1]}")

assert not m_strict.BDSPPatientID.duplicated().any(), "duplicate patients in analysis.parquet"
assert m_strict.papoff_strict.all(), "a non PAP-off row survived"
assert (m_strict.oximetry_bad == 0).all() and (m_strict.fu_valid == 1).all(), "cohort filter leak"
print("\nassertions passed: one row per patient, PAP-off only, cohort filter intact")

print("\ncolumn groups in analysis.parquet")
print(f"  identifiers       : BDSPPatientID, site_id, SiteID, BIDSFolder, SessionID")
print(f"  exposures         : t90_all + {len(STAGES)} stages, raw and _g30 gated")
print(f"  stage minutes     : min_<stage> for {len(STAGES)} stages, plus rec_min tst_min "
      f"staged_min tst_papoff_min")
print(f"  stage descriptives: mean_<stage>, nadir_<stage>, wakepre and wakepost")
print(f"  covariates        : AgeAtVisit, sex, male, site_id, AHI, TST_min, plus the paper's "
      f"oximetry panel")
print(f"  outcomes          : {len(OUTCOME_KEYS)} keys x prevalent/incident/years "
      f"= {len(OUTCOME_COLS)} columns")
print(f"  subgroup          : healthy (H1), healthy_no_riskfactors (H2), healthy_no_thyroid (H3),"
      f" n_prevalent_excl")

hdr("MODELLING NOTES for the next step, not run here")
print("  1. Drop one age-spline column. cr(a, df=4) then .iloc[:, 1:], and fit penalizer=0.")
print("     An unpenalized rank-deficient Cox silently fits nothing. Pattern is in")
print("     T90_Manuscript/numbers/run_primary_unpenalized.py line 37.")
print("  2. Site-stratified Cox, sex as a covariate, PSG date is time zero.")
print("  3. Rank-inverse-normal the exposure WITHIN site, then report HR per 1 SD.")
print("  4. Never put NREM and N2 in the same model.")
print("  5. The negative-control panel is the settled five of 2026-08-07 (back pain, cataract,")
print("     glaucoma, contact dermatitis, hemorrhoids). Fracture and osteoarthritis are ordinary")
print("     outcomes. Re-read the controls on this dataset before believing any stage ranking.")
