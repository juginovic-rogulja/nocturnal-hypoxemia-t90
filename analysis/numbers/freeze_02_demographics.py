"""
Phase 1.2b. Race, ethnicity and body mass index in the BDSP analysis cohort.

CORRECTION. The first two versions of this file were wrong in 2 ways, both found by audit
before anything was published:

  1. BMI was searched with an anchored exact match on the measurement source value, which
     matched only literal "BMI" strings. The archive records anthropometrics as height_cm,
     daily_weight_kg, "Weight in kg" and similar, so the search returned patients from 2
     sites that happen to use the literal label and manufactured a false site contrast.
  2. Race was read from a derived demographics CSV holding 1265 populated rows rather than
     from the OMOP person table it was derived from, which holds 19 262.

The corrected position is that race is available for essentially the whole cohort, and BMI is
derivable for a minority whose members differ systematically from the rest.

Writes numbers/demographics_bdsp.json.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../..")))  # repository root, where paths.py lives
import paths
import warnings
warnings.filterwarnings("ignore")
import glob
import json
import numpy as np
import pandas as pd
import duckdb

OUT = paths.NUMBERS_DIR
CACHE = paths.OMOP_CACHE_DIR
PERSON = glob.glob(f"{paths.OMOP_CACHE_DIR}/**/person_merged.parquet",
                   recursive=True)[0]

import sys
sys.path.insert(0, paths.ANALYSIS_DIR)
from cohort_spec import apply_cohort

# the cohort filter lives in one place. t90_base4 predates the oximetry QC stamp, so the
# analysis file is the only one that carries it.
b = pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet",
                    columns=["BDSPPatientID", "fu_valid", "oximetry_bad", "spo2_pct_below_90", "site_id",
                             "AgeAtVisit", "sex", "obesity_prevalent"])
coh = apply_cohort(b)
N = len(coh)

# ------------------------------------------------------------------ race and ethnicity
p = duckdb.sql(f"""select person_id as "BDSPPatientID", race_source_value,
                          ethnicity_source_value from '{PERSON}'""").df()
j = coh.merge(p, on="BDSPPatientID", how="left")
raw = j.race_source_value.astype(str).str.strip().str.upper()

WHITE = {"WHITE", "CAUCASIAN OR WHITE", "CAUCASIAN"}
BLACK = {"AFRICAN AMERICAN  OR BLACK", "AFRICAN AMERICAN OR BLACK",
         "BLACK/AFRICAN AMERICAN", "BLACK", "AFRICAN AMERICAN"}
ASIAN = {"ASIAN"}
NHPI = {"NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER", "PACIFIC ISLANDER",
        "NATIVE HAWAIIAN"}
AIAN = {"AMERICAN INDIAN OR ALASKA NATIVE", "AMERICAN INDIAN"}
UNK = {"UNKNOWN, UNAVAILABLE OR UNREPORTED", "UNKNOWN/NOT SPECIFIED", "UNKNOWN",
       "UNABLE TO OBTAIN", "DECLINED TO ANSWER", "NONE", "PREFER NOT TO SAY",
       "NAN", "", "DECLINED", "NOT RECORDED", "PATIENT DECLINED", "NOT SPECIFIED",
       "UNAVAILABLE", "DECLINE TO ANSWER"}


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


j["race_h"] = raw.apply(harmonize)
known = j[j.race_h != "Unknown or not reported"]

R = {"n_cohort": int(N)}
R["race"] = {
    "n_with_any_value": int(j.race_source_value.notna().sum()),
    "pct_with_any_value": round(100 * float(j.race_source_value.notna().mean()), 1),
    "n_known_category": int(len(known)),
    "pct_known_category": round(100 * len(known) / N, 1),
    "counts": {k: int(v) for k, v in j.race_h.value_counts().items()},
    "pct_of_cohort": {k: round(100 * v / N, 1) for k, v in j.race_h.value_counts().items()},
    "pct_of_known": {k: round(100 * v / len(known), 1)
                     for k, v in known.race_h.value_counts().items()},
    "source": "OMOP person table, race_source_value, harmonized to OMB categories",
}

eth = j.ethnicity_source_value.astype(str).str.strip().str.upper()
R["ethnicity"] = {
    "n_with_any_value": int(j.ethnicity_source_value.notna().sum()),
    "top_values_by_site": {
        str(s): {str(k): int(v) for k, v in
                 eth[j.site_id == s].value_counts().head(6).items()}
        for s in sorted(j.site_id.dropna().unique())},
    "note": ("The 2 hospitals record different quantities under ethnicity. One records the "
             "Hispanic or Latino axis and the other records ancestry, so ethnicity is "
             "reported by site and is not pooled."),
}

# ------------------------------------------------------------------ body mass index
#
# 2026-08-08. The query below is kept so the defect is on the record, and it raises so nothing
# can fall back to it.
#
# It took ANY measurement whose source value matched "weight" with a number between 25 and 300
# and treated it as kilograms. In this OMOP cache that set includes "WEIGHT" recorded in LBS
# (99,072 rows), "Weight" in Pounds (41,049), "Weight in lbs" (21,408), and rows that are not a
# measured weight at all, such as "WEIGHT IN (LB) TO HAVE BMI = 25" and "FEMALE PREDICTED
# WEIGHT". A per-patient median was then taken across that mixture, so a patient with both
# pound and kilogram rows got a number that is neither.
#
# The published Table 1 value, 34.4 (28.8-42.3), came from this. On the unit-clean derivation
# the cohort median is 32.91 (27.88-39.43) on the same 33.0% of the cohort. The published 34.4
# is within rounding of the contaminated MEAN, not its median.
#
# pct_ge30_of_measured carried a second, independent error: it was computed with .mean() over
# the whole cohort frame, where a missing body mass index counts as a non-obese patient. It
# read 22.9%, which is impossible for a population whose median is above 30 and is the reason
# the row could be spotted as wrong from the table alone. On the clean values it is 64.3% of
# those measured.
BMI_CACHE = f"{OUT}/_bmi_derived.parquet"


def _contaminated_weight_query():
    raise RuntimeError(
        "freeze_02_demographics: the pound-contaminated height and weight query was removed on "
        "2026-08-08. It read any 'weight' measurement between 25 and 300 as kilograms, which in "
        "this cache mixes LBS, Pounds and non-weight rows such as 'WEIGHT IN (LB) TO HAVE "
        "BMI = 25'. Use numbers/_bmi_derived.parquet, which carries wt_kg, ht_cm and bmi with "
        "the units resolved per row.")


bmi_cache = pd.read_parquet(BMI_CACHE)
assert {"BDSPPatientID", "bmi"} <= set(bmi_cache.columns), bmi_cache.columns
bmi_cache = bmi_cache[["BDSPPatientID", "bmi"]].dropna().drop_duplicates("BDSPPatientID")
jb = coh.merge(bmi_cache, on="BDSPPatientID", how="left")
have = jb.bmi.notna()

R["bmi"] = {
    "n_derivable": int(have.sum()),
    "pct_derivable": round(100 * float(have.mean()), 1),
    "median": round(float(jb.bmi.median()), 2) if have.sum() else None,
    "q1": round(float(jb.bmi.quantile(.25)), 2) if have.sum() else None,
    "q3": round(float(jb.bmi.quantile(.75)), 2) if have.sum() else None,
    # of the MEASURED, which is what the key says. It used to divide by the whole cohort.
    "pct_ge30_of_measured": (round(100 * float((jb.loc[have, "bmi"] >= 30).mean()), 1)
                             if have.sum() else None),
    "mean": round(float(jb.bmi.mean()), 2) if have.sum() else None,
    "source": "numbers/_bmi_derived.parquet, units resolved per row",
    "by_site": {str(s): round(100 * float(jb.loc[jb.site_id == s, "bmi"].notna().mean()), 1)
                for s in sorted(jb.site_id.dropna().unique())},
    "selection_check": {
        "age_with": round(float(jb.loc[have, "AgeAtVisit"].median()), 1),
        "age_without": round(float(jb.loc[~have, "AgeAtVisit"].median()), 1),
        "t90_with": round(float(jb.loc[have, "spo2_pct_below_90"].median()), 2),
        "t90_without": round(float(jb.loc[~have, "spo2_pct_below_90"].median()), 2),
        "obesity_dx_with_pct": round(100 * float(jb.loc[have, "obesity_prevalent"].mean()), 1),
        "obesity_dx_without_pct": round(100 * float(jb.loc[~have, "obesity_prevalent"].mean()), 1),
    },
    "note": ("Body mass index was derived from recorded height and weight rather than from a "
             "stored index, which is absent from this archive. Participants with a derivable "
             "value differed systematically from those without, being older and more "
             "hypoxemic, because much of the recorded weight data originates from inpatient "
             "encounters. Body mass index was therefore not used for adjustment, and a "
             "recorded diagnosis of obesity, available for the whole cohort, is reported "
             "instead."),
}
R["obesity_diagnosis_pct"] = round(100 * float(coh.obesity_prevalent.mean()), 1)

# The superseded body mass index block, kept so the published values stay on the record and the
# size of the defect is visible without re-running anything. Every value in it is COMPUTED, by
# numbers/_reproduce_contaminated_bmi.py, which reruns the removed query once for the record. It
# reproduces the published Table 1 row exactly.
with open(f"{OUT}/_bmi_contaminated_reproduction.json") as f:
    _old = json.load(f)
R["superseded_bmi_pound_contaminated"] = {
    **_old,
    "retired": "2026-08-08, decision B11",
    "appeared_in": "Table 1 as 34.4 (28.8-42.3), in 33.0%",
    "replaced_by": {"median": R["bmi"]["median"], "q1": R["bmi"]["q1"], "q3": R["bmi"]["q3"],
                    "n_derivable": R["bmi"]["n_derivable"],
                    "pct_ge30_of_measured": R["bmi"]["pct_ge30_of_measured"]},
    "note": ("Two separate errors. The derivation read pound-denominated and non-weight rows as "
             "kilograms, which raised the median from 32.91 to 34.37. Separately, the share at "
             "or above 30 was divided by the whole cohort rather than by those measured, which "
             "printed 22.9% for a population whose median is above 30. On the clean values the "
             "share is 64.3% of those measured. Refitting all 54 conditions on the clean "
             "covariate flips no verdict."),
}

with open(f"{OUT}/demographics_bdsp.json", "w") as f:
    json.dump(R, f, indent=2)

print(f"cohort n={N:,}")
print(f"\nRACE: any value {R['race']['pct_with_any_value']}%   "
      f"known category {R['race']['pct_known_category']}%")
for k, v in R["race"]["pct_of_cohort"].items():
    print(f"   {k:<45}{R['race']['counts'][k]:>7,}  {v:>5}% of cohort")
print(f"\nBMI: derivable for {R['bmi']['n_derivable']:,} ({R['bmi']['pct_derivable']}%)"
      f"   median {R['bmi']['median']}   by site {R['bmi']['by_site']}")
sc = R["bmi"]["selection_check"]
print(f"   selection: with BMI age {sc['age_with']} T90 {sc['t90_with']}, "
      f"without age {sc['age_without']} T90 {sc['t90_without']}")
print(f"   obesity diagnosis: with {sc['obesity_dx_with_pct']}% vs without "
      f"{sc['obesity_dx_without_pct']}%")
print(f"\nrecorded obesity diagnosis in whole cohort: {R['obesity_diagnosis_pct']}%")
print(f"\nwritten -> {OUT}/demographics_bdsp.json")
