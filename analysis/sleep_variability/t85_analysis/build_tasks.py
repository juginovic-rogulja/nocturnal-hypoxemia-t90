"""
Task list for the T85 extraction.

This is hypoxic_burden/build_tasks.py with one change: the local frozen reference now
also carries spo2_pct_below_88, because the validation gate checks BOTH frozen oximetry
columns and not just T90. The task list itself (BDSPPatientID, site, sub, ses) is
byte-identical to the 2026-08-18 run, which is the point: the same 19,173 patients and
the same H5 session for each.

The validation subset is the pilot's 25 I0002 patients (whose frozen answer is already
known from the 2026-08-13 raw pilot) plus 15 I0006 patients stratified by frozen T90,
drawn with the same seed as the 2026-08-18 run so it is the same 40 records. The GATE is
scored on the 25 pilot patients only, per the task spec.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import os
import sys

import numpy as np
import pandas as pd

T90ROOT = paths.T90_ROOT
HERE = os.path.dirname(os.path.abspath(__file__))
MAP = (f"{paths.FEATURE_WORK}/"
       "atlas_arch_pass/cohort_arch_pass.csv")
PILOT = (f"{paths.SV_ROOT}/"
         "raw_pilot_rederive/pilot_manifest.csv")

sys.path.insert(0, paths.ANALYSIS_DIR)
from cohort_spec import apply_cohort, COHORT_N            # noqa: E402

d = apply_cohort(pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet"))
assert len(d) == COHORT_N, len(d)
d["BDSPPatientID"] = d.BDSPPatientID.astype(int)

m = pd.read_csv(MAP, low_memory=False)[["BDSPPatientID", "bids_col", "SessionID"]]
m["BDSPPatientID"] = m.BDSPPatientID.astype(int)
m = m.drop_duplicates("BDSPPatientID")

t = d[["BDSPPatientID", "site_id", "psg_date", "spo2_pct_below_90", "spo2_pct_below_88",
       "TST_min", "recording_dur_min", "spo2_mean", "AHI"]].merge(m, on="BDSPPatientID",
                                                                  how="left")
assert t.bids_col.notna().all(), t.bids_col.isna().sum()
assert (t.bids_col == "sub-" + t.site_id + t.BDSPPatientID.astype(str)).all()

t = t.rename(columns={"site_id": "site", "bids_col": "sub", "SessionID": "ses"})
t["ses"] = t.ses.astype(int)
tasks = t[["BDSPPatientID", "site", "sub", "ses"]]
tasks.to_csv(f"{HERE}/tasks_all.csv", index=False)
print(f"tasks_all.csv  n={len(tasks)}  sites={t.site.value_counts().to_dict()}")

# the 2026-08-18 task list must be reproduced exactly, or this is a different run
prev = f"{os.path.dirname(HERE)}/hypoxic_burden/tasks_all.csv"
if os.path.exists(prev):
    p = pd.read_csv(prev)
    same = p.equals(tasks.reset_index(drop=True))
    print(f"task list identical to the 2026-08-18 hypoxic-burden run: {same}")
    assert same, "task list drifted from the proven run"

# frozen reference values, kept local, never shipped to EC2
t.to_csv(f"{HERE}/frozen_reference.csv", index=False)

# ------------------------------------------------------------------ validation subset
pil = pd.read_csv(PILOT)[["BDSPPatientID"]]
pil["BDSPPatientID"] = pil.BDSPPatientID.astype(int)
pil.to_csv(f"{HERE}/pilot_ids.csv", index=False)
v1 = tasks[tasks.BDSPPatientID.isin(set(pil.BDSPPatientID))]
print(f"pilot patients found in cohort: {len(v1)} of {len(pil)}")
assert len(v1) == 25, f"expected the pilot's 25, got {len(v1)}"

rng = np.random.RandomState(20260818)          # same seed as the proven run
i6 = t[t.site == "I0006"].copy()
q = i6.spo2_pct_below_90
strata = {"near_zero": i6[q <= 0.1], "mid": i6[(q >= 1.0) & (q <= 10.0)],
          "high": i6[q > 10.0]}
picks = []
for name, pool in strata.items():
    picks.append(pool.sample(n=5, random_state=rng))
    print(f"  I0006 {name}: pool {len(pool)}, took 5")
v2 = pd.concat(picks)[["BDSPPatientID", "site", "sub", "ses"]]

val = pd.concat([v1, v2]).drop_duplicates("BDSPPatientID").reset_index(drop=True)
val.to_csv(f"{HERE}/tasks_validation.csv", index=False)
print(f"tasks_validation.csv  n={len(val)}  sites={val.site.value_counts().to_dict()}")
