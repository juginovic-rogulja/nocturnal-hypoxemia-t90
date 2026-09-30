"""
Phase 1.1 and 1.2. Freeze the cohort description numbers for all three cohorts.

Everything the manuscript states about who was studied comes from this file and nowhere else.
Numbers quoted from the working conversation are treated as unverified until reproduced here.

Writes numbers/cohorts.json.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../..")))  # repository root, where paths.py lives
import paths
import warnings
warnings.filterwarnings("ignore")
import json
import os
import numpy as np
import pandas as pd

OUT = paths.NUMBERS_DIR
MR = os.path.expanduser(paths.MROS_DIR)
SH = paths.SHHS_DIR
R = {}


def q(v, dp=1):
    v = pd.Series(v).dropna()
    return {"n": int(v.notna().sum()), "median": round(float(v.median()), dp),
            "q1": round(float(v.quantile(.25)), dp), "q3": round(float(v.quantile(.75)), dp),
            "mean": round(float(v.mean()), dp), "sd": round(float(v.std()), dp)}


# ----------------------------------------------------------------- BDSP
# t90_base4 carries the observation-window dates but predates the oximetry QC stamp, so the
# flag is merged in from the analysis file rather than either frame being abandoned
DF = paths.TABLES_DIR
b = pd.read_parquet(f"{DF}/t90_base4.parquet")
if "oximetry_bad" not in b.columns:
    qc = pd.read_parquet(f"{DF}/t90_final.parquet", columns=["BDSPPatientID", "oximetry_bad"])
    b = b.merge(qc.drop_duplicates("BDSPPatientID"), on="BDSPPatientID", how="left")
    b["oximetry_bad"] = b.oximetry_bad.fillna(0).astype(int)
R["bdsp"] = {"source": f"{DF}/t90_base4.parquet + oximetry QC from t90_final.parquet",
             "n_master": int(len(b)),
             "n_unique_patients": int(b.BDSPPatientID.nunique())}
v = b[(b.fu_valid == 1) & (b.oximetry_bad == 0)]
R["bdsp"]["n_valid_followup"] = int(len(v))
a = v[v.spo2_pct_below_90.notna()]
R["bdsp"]["n_analysis"] = int(len(a))
R["bdsp"]["sites"] = {str(k): int(x) for k, x in a.site_id.value_counts().items()}
R["bdsp"]["site_pct"] = {str(k): round(100 * x / len(a), 1)
                         for k, x in a.site_id.value_counts().items()}
R["bdsp"]["age"] = q(a.AgeAtVisit, 0)
R["bdsp"]["male_pct"] = round(100 * (a.sex.astype(str).str.upper().str[0] == "M").mean(), 1)
R["bdsp"]["male_n"] = int((a.sex.astype(str).str.upper().str[0] == "M").sum())
R["bdsp"]["t90"] = q(a.spo2_pct_below_90, 2)
R["bdsp"]["t90_above_10_n"] = int((a.spo2_pct_below_90 > 10).sum())
R["bdsp"]["t90_above_10_pct"] = round(100 * (a.spo2_pct_below_90 > 10).mean(), 1)
if "AHI" in a.columns:
    R["bdsp"]["ahi"] = q(a.AHI, 1)
    R["bdsp"]["ahi_above_30_pct"] = round(100 * (a.AHI > 30).mean(), 1)
if "BMI" in a.columns:
    R["bdsp"]["bmi"] = q(a.BMI, 1)
if "TST_min" in a.columns:
    R["bdsp"]["tst_hours"] = q(a.TST_min / 60.0, 2)
if "sleep_efficiency_pct" in a.columns:
    R["bdsp"]["sleep_efficiency_pct"] = q(a.sleep_efficiency_pct, 1)   # v7: read by prep_figure1_profile
if "arousal_index" in a.columns:
    R["bdsp"]["arousal_index"] = q(a.arousal_index, 1)
# follow-up: use death_years as the census clock, and the observation end dates
R["bdsp"]["followup_years"] = q(a.death_years, 2)
oe = pd.to_datetime(a.obs_end, errors="coerce")
oe = oe[(oe > pd.Timestamp("1990-01-01")) & (oe < pd.Timestamp("2027-01-01"))]
R["bdsp"]["final_followup_date_max"] = str(oe.max().date())
R["bdsp"]["final_followup_date_median"] = str(oe.median().date())
pd_ = pd.read_parquet(f"{DF}/t90_final.parquet", columns=["BDSPPatientID", "psg_date"])
psd = pd.to_datetime(pd_.psg_date, errors="coerce")
R["bdsp"]["psg_date_range"] = [str(psd.min().date()), str(psd.max().date())]
# race and ethnicity, required by JAMA
for c in ["race", "Race", "ethnicity", "RaceDSC"]:
    if c in a.columns:
        R["bdsp"]["race_counts"] = {str(k): int(x) for k, x in
                                    a[c].value_counts(dropna=False).head(12).items()}
        break

# ----------------------------------------------------------------- SHHS
d = pd.read_csv(f"{SH}/shhs1-dataset-0.21.0.csv", low_memory=False)
cv = pd.read_csv(f"{SH}/shhs-cvd-summary-dataset-0.21.0.csv", low_memory=False)
cv = cv.drop(columns=[c for c in ["gender", "race", "age_s1"] if c in cv.columns])
s = d[["nsrrid", "pctsa90h", "ahi_a0h4a", "ahi_a0h3a", "age_s1", "gender", "bmi_s1",
       "race", "slpprdp"]].merge(cv, on="nsrrid", how="inner")
R["shhs"] = {"source": f"{SH}/shhs1-dataset-0.21.0.csv + shhs-cvd-summary-dataset-0.21.0.csv",
             "n_analysis": int(len(s)),
             "age": q(s.age_s1, 0),
             "male_pct": round(100 * (s.gender == 1).mean(), 1),
             "male_n": int((s.gender == 1).sum()),
             "bmi": q(s.bmi_s1, 1),
             "t90": q(s.pctsa90h, 2),
             "t90_above_10_pct": round(100 * (s.pctsa90h > 10).mean(), 1),
             "t90_zero_pct": round(100 * (s.pctsa90h == 0).mean(), 1),
             "ahi_4pct": q(s.ahi_a0h4a, 1),
             "ahi_3pct_or_arousal": q(s.ahi_a0h3a, 1),
             "ahi_above_30_pct": round(100 * (s.ahi_a0h4a > 30).mean(), 1),
             "tst_hours": q(s.slpprdp / 60.0, 2),
             "followup_years": q(s.censdate / 365.25, 2),
             "deaths": int((s.vital == 0).sum()),
             "race_counts": {str(k): int(x) for k, x in s.race.value_counts().items()},
             "enrollment_window": "1995-1998 (SHHS visit 1, per NSRR documentation)"}

# ----------------------------------------------------------------- MrOS
pos = pd.read_sas(f"{MR}/POSFEB23.SAS7BDAT", encoding="latin-1")
ef = pd.read_sas(f"{MR}/effeb24.sas7bdat", encoding="latin-1")
vs = pd.read_sas(f"{MR}/vsfeb24.sas7bdat", encoding="latin-1")
m = (pos[["ID", "POPCSA90", "PORDI4P", "POSLPRDP"]]
     .merge(ef[["ID", "SITE", "DADEAD", "DACARDIO", "DACANCER", "FUCDTIME", "FUVSDT"]], on="ID")
     .merge(vs[["ID", "VSAGE1", "HWBMI"]], on="ID", how="left"))
for c in m.columns:
    if c not in ("ID", "SITE"):
        m[c] = pd.to_numeric(m[c], errors="coerce")
R["mros"] = {"source": "MrOS POSFEB23 + EFFEB24 + VSFEB24 (NSRR)",
             "n_analysis": int(len(m)),
             "sites": {str(k): int(x) for k, x in m.SITE.astype(str).value_counts().items()},
             "age": q(m.VSAGE1, 0),
             "male_pct": 100.0,
             "bmi": q(m.HWBMI, 1),
             "t90": q(m.POPCSA90, 2),
             "t90_above_10_pct": round(100 * (m.POPCSA90 > 10).mean(), 1),
             "t90_zero_pct": round(100 * (m.POPCSA90 == 0).mean(), 1),
             "ahi_4pct": q(m.PORDI4P, 1),
             "tst_hours": q(m.POSLPRDP / 60.0, 2),
             # FUCDTIME runs from the parent-cohort baseline visit, a median of 3.4 years
             # before the sleep study, so using it would start the clock before the exposure.
             # FUVSDT runs from the sleep visit itself, which is the correct time zero.
             "followup_years": q(m.FUVSDT / 365.25, 2),
             "followup_years_from_parent_visit": q(m.FUCDTIME / 365.25, 2),
             "deaths": int(m.DADEAD.fillna(0).sum()),
             "cvd_deaths": int(m.DACARDIO.fillna(0).sum()),
             "enrollment_window": "2003-2005 (MrOS sleep visit, per NSRR documentation)"}

R["totals"] = {"n_all_three": R["bdsp"]["n_analysis"] + R["shhs"]["n_analysis"]
               + R["mros"]["n_analysis"]}

with open(f"{OUT}/cohorts.json", "w") as f:
    json.dump(R, f, indent=2)

for k in ("bdsp", "shhs", "mros"):
    c = R[k]
    print(f"\n=== {k.upper()}  n={c['n_analysis']:,}")
    print(f"    age {c['age']['median']:.0f} ({c['age']['q1']:.0f}-{c['age']['q3']:.0f})"
          f"   male {c['male_pct']}%   BMI {c.get('bmi',{}).get('median','-')}")
    print(f"    T90 median {c['t90']['median']}%  above 10%: {c['t90_above_10_pct']}%")
    print(f"    follow-up {c['followup_years']['median']} y")
print(f"\nTOTAL across three cohorts: {R['totals']['n_all_three']:,}")
print(f"BDSP final follow-up: {R['bdsp']['final_followup_date_max']}")
print(f"\nwritten -> {OUT}/cohorts.json")
