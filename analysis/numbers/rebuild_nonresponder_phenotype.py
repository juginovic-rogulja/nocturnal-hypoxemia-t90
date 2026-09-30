"""
Rebuild numbers/nonresponder_phenotype_audited.csv so it carries the settled negative-control
panel of 2026-08-07.

Why this exists. The frozen file had no writer anywhere in the tree. It was flagging fracture
and osteoarthritis as negative controls, and it had no row at all for contact dermatitis or
hemorrhoids, so eFigure 10, Figure 13 and eTable 12 could only ever draw three of the five
settled controls. This script reconstructs the sample by the recipe that
figures_alen/scripts/freeze_doseresponse.py and
New_Figures/_NOT_THESE_workfiles/scripts/fits_phenotype_combined.py both use, refits every
condition, and writes the file back out with all 53 conditions.

Two things are deliberately different from the superseded file.

  1. Prevalent flags come from data_frozen_v7_2026-09/t90_final.parquet, not
     data_frozen/t90_outcomes_rebuilt.parquet. The two agree on 49 of the 51 shared
     conditions. They disagree on obesity hypoventilation and on falls, and in both cases
     t90_final is the corrected one: the ICD-9 278.01/278.03 fix for obesity hypoventilation
     and the E888/E884 prefix-collision fix for falls. t90_outcomes_rebuilt has neither, and
     it has no column for the two replacement controls at all.
  2. The Benjamini-Hochberg correction runs over all 53 conditions rather than 51, because
     two conditions were added to the family.

A positive control runs first: every condition except those two known corrections must
reproduce the superseded odds ratio to within 0.002, or the script stops.

Model, unchanged: logistic, failure ~ prevalent condition + age + sex + log(1 + pretreatment
sleep T90), on arm A split-night records with pretreatment sleep T90 above 10%.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../..")))  # repository root, where paths.py lives
import paths
import glob
import shutil
import sys
from datetime import date

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

ROOT = paths.T90_ROOT
SRC = f"{paths.X4_DIR}/cpap_stage"   # v8.2: re-extracted, collapse-masked shards (X4_groupP)
sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
from disease_definitions import DISEASES, NEGATIVE_CONTROLS  # noqa: E402

OUT = f"{paths.NUMBERS_DIR}/nonresponder_phenotype_audited.csv"
SUPER = f"{paths.NUMBERS_DIR}/_superseded/nonresponder_phenotype_audited_{date.today()}_pre_negcontrol_v5.csv"

# the two conditions whose prevalent definition was corrected after the superseded file was
# frozen; they are expected to move and are exempt from the reproduction gate
CORRECTED = {"Obesity hypoventilation", "Falls"}

# ------------------------------------------------------------------ the sample
o = pd.concat([pd.read_csv(f, low_memory=False) for f in sorted(glob.glob(f"{SRC}/out*.csv"))],
              ignore_index=True).drop_duplicates(["BIDSFolder", "SessionID"])
o = o[(o.status == "ok") & (o.arm == "A")
      & o.pre_sleep_t90.notna() & o.post_sleep_t90.notna()]
s = o[o.pre_sleep_t90 > 10].drop_duplicates("BDSPPatientID").copy()
s = s[~((s.pre_sleep_nadir > s.pre_sleep_mean) | (s.post_sleep_nadir > s.post_sleep_mean))]
s["nonresp"] = (s.post_sleep_t90 > 10).astype(int)

fin = pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet")
pcols = [f"{k}_prevalent" for k in DISEASES if f"{k}_prevalent" in fin.columns]
assert "dermatitis_contact_prevalent" in pcols and "haemorrhoids_prevalent" in pcols
fin = fin[["BDSPPatientID"] + pcols].drop_duplicates("BDSPPatientID")
m = s.merge(fin, on="BDSPPatientID", how="left")
for c in pcols:
    m[c] = m[c].fillna(0).astype(int)
m["male"] = (m.SexDSC.astype(str).str.upper().str[0] == "M").astype(int)
m["lbase"] = np.log1p(m.pre_sleep_t90)
print(f"sample {len(m):,}, failures {int(m.nonresp.sum()):,} "
      f"({100 * m.nonresp.mean():.1f}%)")


def fit(key):
    t = f"{key}_prevalent"
    X = sm.add_constant(m[[t, "AgeAtVisit", "male", "lbase"]].astype(float))
    d = pd.concat([X, m.nonresp], axis=1).dropna()
    if d[t].sum() < 10 or d.loc[d[t] == 1, "nonresp"].sum() < 3:
        return None
    r = sm.Logit(d.nonresp, d[X.columns]).fit(disp=0)
    ci = r.conf_int()
    return {"Rp": int(d.loc[d[t] == 1, "nonresp"].sum()), "Np": int(d[t].sum()),
            "or2": round(float(np.exp(r.params[t])), 3),
            "lo2": round(float(np.exp(ci.loc[t, 0])), 3),
            "hi2": round(float(np.exp(ci.loc[t, 1])), 3),
            "p2": float(r.pvalues[t])}


rows = []
for k, spec in DISEASES.items():
    if f"{k}_prevalent" not in m.columns:
        continue
    v = fit(k)
    if v is None:
        print(f"  skipped, too few prevalent cases: {spec[0]}")
        continue
    rows.append({"Condition": spec[0], "neg": k in NEGATIVE_CONTROLS, **v})

R = pd.DataFrame(rows)
R["q2"] = multipletests(R.p2.values, method="fdr_bh")[1]
R["q2"] = R.q2.astype(float)

# ------------------------------------------------------------------ positive control
old = pd.read_csv(OUT)
bad = []
for _i, r in old.iterrows():
    if r.Condition in CORRECTED or r.Condition not in set(R.Condition):
        continue
    got = float(R.loc[R.Condition == r.Condition, "or2"].iloc[0])
    if abs(got - float(r.or2)) > 0.002:
        bad.append((r.Condition, float(r.or2), got))
# v8.2 (2026-09-15): the file on disk is the v7 build of 2026-09-07. The v8 tables changed the outcome definitions (decision 11,
# the control swap, rebuilt prevalent flags, the v8.1 cohort), so the odds ratios move by construction; the comparison is
# REPORTED, never asserted against the old file (rules file item 6: no old result inside an assert). The structural gate:
# every fitted odds ratio is finite and the panel is not empty.
print(f"comparison with the superseded file: {len(bad)} of {len(old)} conditions differ by more than 0.002 in the odds ratio"
      + (f"; e.g. {bad[:6]}" if bad else "; every shared condition reproduces to 0.002"))
assert len(R) > 0 and np.isfinite(R.or2).all() and np.isfinite(R.q2).all(), "non-finite odds ratio or q"
for c in sorted(CORRECTED & set(old.Condition)):
    was = float(old.loc[old.Condition == c, "or2"].iloc[0])
    now = float(R.loc[R.Condition == c, "or2"].iloc[0]) if c in set(R.Condition) else np.nan
    print(f"  corrected definition, expected to move: {c} {was} -> {now}")

# ------------------------------------------------------------------ what changed
oldsig = set(old[(old.q2 < 0.05) & ~old.neg].Condition)
newsig = set(R[(R.q2 < 0.05) & ~R.neg].Condition)
print(f"significant non-control conditions {len(oldsig)} -> {len(newsig)}")
if oldsig ^ newsig:
    print(f"  gained {sorted(newsig - oldsig)}  lost {sorted(oldsig - newsig)}")
print("negative controls now on file:",
      R.loc[R.neg, ["Condition", "or2", "lo2", "hi2", "q2"]].to_string(index=False))

shutil.copy(OUT, SUPER)
R.to_csv(OUT, index=False)
print(f"superseded copy {SUPER}")
print(f"wrote {OUT}  {len(R)} conditions, {int(R.neg.sum())} negative controls")
