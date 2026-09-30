"""
Task list for the overnight-oxygen-profile extraction.

This is t_ladder_up/build_tasks.py with the frozen reference widened. The task list itself
(BDSPPatientID, site, sub, ses) must come out byte-identical to the 2026-08-21 threshold-ladder
run, which is the point: the same 19,173 patients from cohort_spec.apply_cohort and the same H5
session for each, so the deliverable lands on the frozen cohort exactly.

The frozen reference gains spo2_mean (the whole-recording mean saturation, which is this run's
primary validation target) and the four frozen stage percentages N1_pct, N2_pct, N3_pct,
REM_pct, which validate the staging path independently of the oximetry path. It stays on this
machine and is never shipped to EC2.

The validation subset is the pilot's 25 I0002 patients (whose frozen answer the 2026-08-13 raw
pilot already established) plus 15 I0006 patients stratified by frozen T90, drawn with the same
seed as the proven runs so it is the same 40 records. The GATE is scored on the 25 pilot
patients only.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import os
import sys

import numpy as np
import pandas as pd

T90ROOT = paths.T90_ROOT
HERE = os.path.dirname(os.path.abspath(__file__))
SIDE = os.path.dirname(HERE)
MAP = (f"{paths.FEATURE_WORK}/"
       "atlas_arch_pass/cohort_arch_pass.csv")
PILOT = f"{SIDE}/raw_pilot_rederive/pilot_manifest.csv"

sys.path.insert(0, paths.ANALYSIS_DIR)
from cohort_spec import apply_cohort, COHORT_N            # noqa: E402

# repointed to the v7 table 2026-09-08 (the 2026-08-21 run read data_frozen/, the April table)
d = apply_cohort(pd.read_parquet(f"{paths.PREV_TABLES_DIR}/t90_final.parquet"))
assert len(d) == COHORT_N, len(d)
d["BDSPPatientID"] = d.BDSPPatientID.astype(int)

m = pd.read_csv(MAP, low_memory=False)[["BDSPPatientID", "bids_col", "SessionID"]]
m["BDSPPatientID"] = m.BDSPPatientID.astype(int)
m = m.drop_duplicates("BDSPPatientID")

KEEP = ["BDSPPatientID", "site_id", "psg_date", "spo2_mean", "spo2_pct_below_90",
        "spo2_pct_below_88", "TST_min", "recording_dur_min", "N1_pct", "N2_pct",
        "N3_pct", "REM_pct", "AHI"]
t = d[KEEP].merge(m, on="BDSPPatientID", how="left")
assert t.bids_col.notna().all(), t.bids_col.isna().sum()
assert (t.bids_col == "sub-" + t.site_id + t.BDSPPatientID.astype(str)).all()

t = t.rename(columns={"site_id": "site", "bids_col": "sub", "SessionID": "ses"})
t["ses"] = t.ses.astype(int)
tasks = t[["BDSPPatientID", "site", "sub", "ses"]]
tasks.to_csv(f"{HERE}/tasks_all.csv", index=False)
print(f"tasks_all.csv  n={len(tasks)}  sites={t.site.value_counts().to_dict()}")

# the proven task list must be reproduced exactly, or this is a different run
prev = f"{SIDE}/t_ladder_up/tasks_all.csv"
assert os.path.exists(prev), prev
p = pd.read_csv(prev)
same = p.equals(tasks.reset_index(drop=True))
print(f"task list identical to the 2026-08-21 threshold-ladder run: {same}")
assert same, "task list drifted from the proven run"

# frozen reference values, kept local, never shipped to EC2
t.to_csv(f"{HERE}/frozen_reference.csv", index=False)
print(f"frozen_reference.csv  n={len(t)}  cols={list(t.columns)}")

# ------------------------------------------------------------------ validation subset
pil = pd.read_csv(PILOT)[["BDSPPatientID"]]
pil["BDSPPatientID"] = pil.BDSPPatientID.astype(int)
pil.to_csv(f"{HERE}/pilot_ids.csv", index=False)
v1 = tasks[tasks.BDSPPatientID.isin(set(pil.BDSPPatientID))]
print(f"pilot patients found in cohort: {len(v1)} of {len(pil)}")
assert len(v1) == 25, f"expected the pilot's 25, got {len(v1)}"

rng = np.random.RandomState(20260818)          # same seed as the proven runs
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

prev_val = f"{SIDE}/t_ladder_up/tasks_validation.csv"
if os.path.exists(prev_val):
    pv = pd.read_csv(prev_val)
    print(f"validation subset identical to the ladder run: {pv.equals(val)}")
    assert pv.equals(val), "validation subset drifted from the proven run"
