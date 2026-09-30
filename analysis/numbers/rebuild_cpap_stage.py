#!/usr/bin/env python3
"""
Rebuild the per-stage oximetry table that the treatment analysis depends on (chain step 100).

v8.2 (2026-09-15): the 72 cpap_stage shards are RE-EXTRACTED from the BDSP source files (group P re-extraction,
V8_RECALC_2026-09-12/X4_groupP/cpap_stage/: raw/ holds the extractor output as landed from EC2, the folder itself
holds the same shards after mask_collapsed.py, which sets the stage-specific values missing on nights whose sleep
sits in one stage code and flags them _v8_stage_collapsed). analyse.py (the pilot's own assembler, copied unchanged)
builds the three designs from those masked shards and writes /tmp/cpap_t90_by_stage.parquet; this script copies the
result into data_frozen_v8_2026-09/ and writes the provenance sidecar. No count is typed: the design sizes are
printed and recorded, and the only assertions are structural.
Writes data_frozen_v8_2026-09/cpap_t90_by_stage.parquet (+ .provenance.json).
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../..")))  # repository root, where paths.py lives
import paths
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
import numpy as np
import pandas as pd
ROOT = Path(paths.T90_ROOT)
sys.path.insert(0, paths.ANALYSIS_DIR)
from cohort_spec import sidecar  # noqa: E402
SHARDS = ROOT / "V8_RECALC_2026-09-12" / "X4_groupP" / "cpap_stage"      # masked shards + analyse.py (v8.2)
TMP = SHARDS / "cpap_t90_by_stage_v8_2.parquet"                          # analyse.py writes here (CPAP_STAGE_OUT), beside the shards, never under /tmp
OUT = ROOT / "data_frozen_v8_2026-09" / "cpap_t90_by_stage.parquet"
shards = sorted(SHARDS.glob("out*_*.csv"))
assert shards, f"no shards under {SHARDS}"
masks = sorted(SHARDS.glob("MASK_*.json"))
assert masks, "mask_collapsed.py reports missing: the shards were not passed through the collapse rule"
if TMP.exists():
    TMP.unlink()
r = subprocess.run([sys.executable, "analyse.py"], cwd=SHARDS, capture_output=True, text=True, env={**os.environ, "CPAP_STAGE_OUT": str(TMP)})
(SHARDS / "analyse_v8_2.log").write_text(r.stdout + "\n--- stderr ---\n" + r.stderr)
if r.returncode or not TMP.exists():
    print(r.stdout[-2000:], r.stderr[-2000:])
    raise SystemExit("analyse.py did not produce the stage table")
d = pd.read_parquet(TMP)
assert {"A_split", "P_placebo", "B_pairs"} <= set(d.design), d.design.unique()
assert d.loc[d.design != "B_pairs", "status"].eq("ok").all(), d.status.value_counts(dropna=False).to_dict()   # pair rows carry dx_/tx_ values, no status (as in v8.1)
assert "_v8_stage_collapsed" in d.columns and d._v8_stage_collapsed.notna().all(), "the shards were not masked (no _v8_stage_collapsed on every design row)"
counts = d.design.value_counts().to_dict()
ncol = int(d._v8_stage_collapsed.astype(int).sum())
# whole-night covariates from the v8 t90_final, as the v8 table builder did for the v8.1 table (E2_tables/scripts/build_v8_tables.py,
# "cpap_t90_by_stage: whole-night AHI, TST_min and spo2_pct_below_90 refreshed from the v8 t90_final"); its per-stage columns are its own.
from cohort_spec import T90_FINAL  # noqa: E402
v8i = pd.read_parquet(T90_FINAL, columns=["BDSPPatientID", "AHI", "TST_min", "spo2_pct_below_90"]).set_index("BDSPPatientID")
idx = d.BDSPPatientID.isin(v8i.index).to_numpy(); ids = d.loc[idx, "BDSPPatientID"].values
nref = {}
for c in ("AHI", "TST_min", "spo2_pct_below_90"):
    if c in d.columns:
        d[c] = d[c].astype(float)
        before = d.loc[idx, c].to_numpy(); after = v8i.loc[ids, c].to_numpy()
        nref[c] = int((~((np.isnan(before) & np.isnan(after)) | (before == after))).sum())
        d.loc[idx, c] = after
d["_v8_covariates_refreshed"] = idx
d.to_parquet(TMP, index=False)
print(f"whole-night covariates refreshed from the v8 t90_final on {int(idx.sum()):,} of {len(d):,} rows; cells changed {nref}")
# the v8.1 table is kept in data_frozen_v8_2026-09/_v8_1_20260915/ (copied by hand before this step), not beside the live table
shutil.copy(TMP, OUT)
print(f"{len(d):,} recordings, {d.BDSPPatientID.nunique():,} patients")
print(f"designs: {counts}   collapsed-staging recordings (stage values missing): {ncol}")
sidecar(str(TMP), __file__, extra_inputs=[str(p) for p in shards[:3]] + [str(m) for m in masks],
        # the frozen folder carries no beside-sidecar: the driver's flat mirror and PROVENANCE_INDEX record the table; this record sits beside the intermediate
        note=f"v8.2 group P re-extraction: designs {counts}, collapsed {ncol}, shards {len(shards)} under {SHARDS}; whole-night AHI/TST_min/spo2_pct_below_90 from the v8 t90_final on {int(idx.sum())} rows (cells changed {nref})")
print(f"written -> {OUT}")
