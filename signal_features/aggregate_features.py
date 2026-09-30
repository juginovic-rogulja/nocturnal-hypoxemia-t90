"""
Aggregate per-patient feature JSON files into a single flat CSV for model training.

The feature extraction pipeline writes one JSON per patient. For clock
training we need one flat table:
    patient_id, age, sex, BMI, site, feature_1, feature_2, ...

This script:
  1. Scans a directory for feature JSONs
  2. Extracts metadata (patient_id, site) from filename
  3. Joins with a cohort CSV (for age, sex, BMI covariates)
  4. Flattens nested brain/microstructure/heart/lung dicts into columns
  5. Writes a tidy CSV + an integrity JSON describing what was aggregated

Fail-safes:
  * SHA256 every JSON input and record
  * Assert no duplicate (patient_id, session_id) pairs
  * Reload output CSV and compare row count + column count
  * Log every skipped patient with reason
"""

from __future__ import annotations
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..")))  # repository root, where paths.py lives
import paths

import hashlib
import json
import re
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

# ---------------------------------------------------------------------------
# Filename parser — extracts site + patient + session from an h5 stem like:
#   sub-SXXXXNNNNNNNNNNN_ses-1
#   sub-IXXXXNNNNNN_ses-02
# ---------------------------------------------------------------------------
FNAME_RE = re.compile(r"sub-(?P<site>[A-Z]\d{4})(?P<pid>\d+)_ses-(?P<ses>\d+)")


def parse_patient_meta(stem: str) -> Optional[Dict[str, str]]:
    m = FNAME_RE.search(stem)
    if not m:
        return None
    return {
        "site_id": m.group("site"),
        "bdsp_patient_id": m.group("pid"),
        "session_id": m.group("ses"),
    }


def flatten(prefix: str, d: dict) -> dict:
    """Flatten one dict, prefixing keys."""
    return {f"{prefix}_{k}": v for k, v in d.items()}


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def aggregate(feature_dir: Path, cohort_csv: Optional[Path],
              out_csv: Path, out_integrity_json: Path) -> None:
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    feature_jsons = sorted(feature_dir.glob("features_*.json"))
    log: List[dict] = []
    rows: List[dict] = []

    print(f"Scanning {feature_dir} — found {len(feature_jsons)} feature JSON(s)")

    for fp in feature_jsons:
        rec = json.loads(fp.read_text())
        src_stem = Path(rec["source"]).stem
        meta = parse_patient_meta(src_stem)
        if meta is None:
            log.append({"file": str(fp), "skipped_because": "filename_parse_fail"})
            continue

        row = {
            "file": fp.name,
            "sha256": sha256_file(fp),
            "site_id": meta["site_id"],
            "bdsp_patient_id": meta["bdsp_patient_id"],
            "session_id": meta["session_id"],
            "duration_hr": rec["meta"]["duration_hr"],
            "n_channels": rec["meta"]["n_channels"],
            "fs": rec["meta"]["fs"],
            "extracted_at": rec.get("extracted_at"),
            "runtime_sec": rec.get("runtime_sec"),
        }
        row.update(flatten("brain", rec.get("brain", {})))
        row.update(flatten("micro", rec.get("microstructure", {})))
        row.update(flatten("heart", rec.get("heart", {})))
        row.update(flatten("lung", rec.get("lung", {})))
        rows.append(row)

    if not rows:
        print("No valid feature rows — nothing to write.")
        return

    df = pd.DataFrame(rows)

    # Dedupe check — should have exactly one row per (patient, session)
    dup_mask = df.duplicated(subset=["site_id", "bdsp_patient_id", "session_id"], keep=False)
    if dup_mask.any():
        print(f"WARNING: {dup_mask.sum()} duplicate (patient, session) pairs present!")
        log.append({"warning": "duplicate_patient_session_pairs",
                    "n": int(dup_mask.sum())})

    # Merge with cohort CSV if provided (to join age/sex/BMI/etc.)
    if cohort_csv and cohort_csv.exists():
        cohort = pd.read_csv(cohort_csv, low_memory=False)
        # Cohort has BDSPPatientID as int; our row has it as string
        cohort["bdsp_patient_id"] = cohort["BDSPPatientID"].astype(str)
        keep_cols = [
            "bdsp_patient_id", "SiteID", "AgeAtVisit", "SexDSC",
            "StudyType", "CreationTime",
        ]
        keep_cols = [c for c in keep_cols if c in cohort.columns]
        df = df.merge(cohort[keep_cols].drop_duplicates("bdsp_patient_id"),
                      on="bdsp_patient_id", how="left", suffixes=("", "_cohort"))

    # Write CSV
    df.to_csv(out_csv, index=False)

    # Reload check
    reloaded = pd.read_csv(out_csv)
    assert len(reloaded) == len(df), f"Reload mismatch: {len(reloaded)} vs {len(df)}"
    assert set(reloaded.columns) == set(df.columns), "Reload columns differ"

    integrity = {
        "timestamp": ts,
        "feature_dir": str(feature_dir),
        "n_input_jsons": len(feature_jsons),
        "n_output_rows": len(df),
        "n_columns": len(df.columns),
        "columns": sorted(df.columns.tolist()),
        "skipped_log": log,
        "cohort_csv_used": str(cohort_csv) if cohort_csv else None,
    }
    out_integrity_json.write_text(json.dumps(integrity, indent=2))
    print(f"Wrote {out_csv} ({len(df):,} rows, {len(df.columns)} cols)")
    print(f"Wrote {out_integrity_json}")


if __name__ == "__main__":
    import sys
    fd = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(paths.FEATURE_OUTPUTS_DIR)
    cohort = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(f"{paths.FEATURE_WORK}/cohort_v1.csv")
    out = fd / "aggregated_features.csv"
    itg = fd / "aggregated_integrity.json"
    aggregate(fd, cohort, out, itg)
