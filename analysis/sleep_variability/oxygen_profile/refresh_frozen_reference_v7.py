"""
Repoint of the oxygen-profile frozen reference to the v7 tables (coordinator, 2026-09-08 22:10 ET).

build_tasks.py (2026-08-21) wrote frozen_reference.csv from data_frozen/t90_final.parquet, the
April table. assemble.py merges spo2_pct_below_90 (and the stage percentages) from that file into
oxyprofile_per_patient.csv, and b01_profile.py bands the cohort on that column, so Supp Fig 2 was
banded on the superseded T90. This script rebuilds ONLY frozen_reference.csv from
data_frozen_v7_2026-09/t90_final.parquet with the same columns; the task list (sub, ses) is taken
from the 2026-08-21 file, which stays as frozen_reference_PRE_V7_2026-08-21.csv. The patient set
must be identical (same 19,173 ids, same site per id) or this stops.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import hashlib
import json
import os
import sys
import time

import pandas as pd

T90ROOT = paths.T90_ROOT
V7 = f"{paths.PREV_TABLES_DIR}/t90_final.parquet"
HERE = os.path.dirname(os.path.abspath(__file__))
OLD = os.path.join(HERE, "frozen_reference_PRE_V7_2026-08-21.csv")
OUT = os.path.join(HERE, "frozen_reference.csv")
sys.path.insert(0, paths.ANALYSIS_DIR)
from cohort_spec import apply_cohort, COHORT_N  # noqa: E402

KEEP = ["BDSPPatientID", "site_id", "psg_date", "spo2_mean", "spo2_pct_below_90",
        "spo2_pct_below_88", "TST_min", "recording_dur_min", "N1_pct", "N2_pct",
        "N3_pct", "REM_pct", "AHI"]

d = apply_cohort(pd.read_parquet(V7))
assert len(d) == COHORT_N, len(d)
d["BDSPPatientID"] = d.BDSPPatientID.astype(int)
missing_cols = [c for c in KEEP if c not in d.columns]
assert not missing_cols, missing_cols
t = d[KEEP].rename(columns={"site_id": "site"})

old = pd.read_csv(OLD, low_memory=False)
old["BDSPPatientID"] = old.BDSPPatientID.astype(int)
assert set(old.BDSPPatientID) == set(t.BDSPPatientID), "patient set differs from the 2026-08-21 run"
t = t.merge(old[["BDSPPatientID", "site", "sub", "ses"]].rename(columns={"site": "site_old"}),
            on="BDSPPatientID", how="left")
assert (t.site == t.site_old).all(), int((t.site != t.site_old).sum())
t = t.drop(columns="site_old")
t["ses"] = t.ses.astype(int)
assert list(t.columns) == list(old.columns), (list(t.columns), list(old.columns))
t = t.sort_values("BDSPPatientID").reset_index(drop=True)
old = old.sort_values("BDSPPatientID").reset_index(drop=True)

# what moved between the April table and v7, for the record
rep = {}
for c in ("spo2_pct_below_90", "spo2_pct_below_88", "spo2_mean", "TST_min", "AHI",
          "N1_pct", "N2_pct", "N3_pct", "REM_pct"):
    a, b = old[c], t[c]
    both = a.notna() & b.notna()
    rep[c] = {"changed_gt_1e-6": int(((a - b).abs() > 1e-6)[both].sum()),
              "old_missing": int(a.isna().sum()), "new_missing": int(b.isna().sum())}
bands = lambda s: {"le1": int((s <= 1).sum()), "1to5": int(((s > 1) & (s <= 5)).sum()),
                   "5to10": int(((s > 5) & (s <= 10)).sum()), "gt10": int((s > 10).sum())}
rep["band_n_old"] = bands(old.spo2_pct_below_90)
rep["band_n_v7"] = bands(t.spo2_pct_below_90)
rep["site_counts"] = t.site.value_counts().to_dict()

t.to_csv(OUT, index=False)
h = hashlib.sha256(open(V7, "rb").read()).hexdigest()
prov = {"output": OUT, "written": time.strftime("%Y-%m-%d %H:%M:%S"), "script": __file__,
        "input_table": V7, "input_sha256": h, "input_mtime": time.strftime(
            "%Y-%m-%d %H:%M:%S", time.localtime(os.path.getmtime(V7))),
        "task_columns_from": OLD, "n": int(len(t)), "report": rep}
json.dump(prov, open(OUT + ".provenance.json", "w"), indent=2)
print(json.dumps(rep, indent=1))
print(f"frozen_reference.csv rewritten from v7, n={len(t)}, sha256(v7 table)={h[:16]}")
