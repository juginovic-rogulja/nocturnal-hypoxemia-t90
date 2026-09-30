#!/usr/bin/env python3
"""
Master merge for ALL extractions, producing aggregated_features_v5.csv.

Inputs (S3):
  - worker[1-4]_final.csv from atlas-resume bucket (already merged via run_merge.py for v6)
  - workerrecovery_*.csv (8500 missing patients) — when done
  - bch_*.csv (12839 pediatric) — when done
  - arch_checkpoint_*.csv + arch final → atlas_architecture.csv (TST/AHI/sleep_eff/spo2_nadir corrected)
  - heedb_final.csv (919 paired ECG features)
  - hedb_final.csv (250-782 paired EEG features)
  - omop_drug_exposure_v6.csv + omop_condition_occurrence_v6.csv

Outputs:
  - aggregated_features_v5.csv (master adult cohort: ~37K patients with full features + arch + OMOP join)
  - aggregated_features_v5_pediatric.csv (BCH ~12K)
  - aggregated_features_v5_heart.csv (HEEDB paired ECG)
  - aggregated_features_v5_eeg.csv (HEDB paired EEG)
  - merge_v5_integrity.json
"""
from __future__ import annotations
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..")))  # repository root, where paths.py lives
import paths
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(paths.FEATURE_WORK)
OUT = ROOT / "outputs"
ATLAS_BUCKET = paths.S3_ATLAS_RESUME_BUCKET
P11_BUCKET = paths.S3_PAPER11_BUCKET

def s3_cp(s3_path: str, local: Path):
    if local.exists() and local.stat().st_size > 0:
        return
    subprocess.run(["aws", "s3", "cp", s3_path, str(local), "--only-show-errors"], check=True)

def main():
    started = time.time()
    report = {"started": time.strftime("%FT%T")}

    cache = ROOT / "atlas_v6_audit/cache_v5"
    cache.mkdir(exist_ok=True, parents=True)

    # ========== 1. Load existing v6 merged cohort ==========
    print("\n=== 1. Atlas v6 base cohort ===")
    v6 = pd.read_csv(ROOT / "atlas_v6_audit/outputs/aggregated_features_v4_z.csv", low_memory=False)
    print(f"  v6 success cohort: {len(v6):,} rows × {v6.shape[1]} cols")
    report["v6_base"] = len(v6)

    # ========== 2. Recovery (when done) ==========
    print("\n=== 2. Atlas Recovery (8500 missing) ===")
    recovery_local = cache / "worker_recovery_final.csv"
    try:
        s3_cp(f"s3://{ATLAS_BUCKET}/outputs/workerrecovery_final.csv", recovery_local)
        rec = pd.read_csv(recovery_local, low_memory=False)
        print(f"  Recovery final: {len(rec):,} rows")
    except Exception as e:
        print(f"  Recovery final not yet — using checkpoint: {e}")
        # Find latest checkpoint
        ls = subprocess.run(["aws", "s3", "ls", f"s3://{ATLAS_BUCKET}/outputs/"], capture_output=True, text=True).stdout
        cp_files = sorted([l.split()[-1] for l in ls.split("\n") if "workerrecovery_checkpoint_" in l and l.endswith(".csv")])
        if cp_files:
            latest = cp_files[-1]
            s3_cp(f"s3://{ATLAS_BUCKET}/outputs/{latest}", cache / latest)
            rec = pd.read_csv(cache / latest, low_memory=False)
            print(f"  Using {latest}: {len(rec):,} rows")
        else:
            rec = pd.DataFrame()
    report["recovery_rows"] = len(rec)

    # ========== 3. BCH pediatric (when done) ==========
    print("\n=== 3. BCH pediatric ===")
    bch_local = cache / "bch_final.csv"
    try:
        s3_cp(f"s3://{ATLAS_BUCKET}/outputs/bch_final.csv", bch_local)
        bch = pd.read_csv(bch_local, low_memory=False)
        print(f"  BCH final: {len(bch):,} rows")
    except Exception:
        ls = subprocess.run(["aws", "s3", "ls", f"s3://{ATLAS_BUCKET}/outputs/"], capture_output=True, text=True).stdout
        cp_files = sorted([l.split()[-1] for l in ls.split("\n") if "bch_checkpoint_" in l and l.endswith(".csv")])
        if cp_files:
            latest = cp_files[-1]
            s3_cp(f"s3://{ATLAS_BUCKET}/outputs/{latest}", cache / latest)
            bch = pd.read_csv(cache / latest, low_memory=False)
            print(f"  Using {latest}: {len(bch):,} rows")
        else:
            bch = pd.DataFrame()
    report["bch_rows"] = len(bch)

    # ========== 4. Atlas arch pass ==========
    print("\n=== 4. Atlas architecture pass (TST/AHI/sleep_eff) ===")
    arch_local = cache / "atlas_architecture.csv"
    try:
        s3_cp(f"s3://{ATLAS_BUCKET}/outputs/atlas_architecture.csv", arch_local)
        arch = pd.read_csv(arch_local, low_memory=False)
        print(f"  arch final: {len(arch):,} rows")
    except Exception:
        ls = subprocess.run(["aws", "s3", "ls", f"s3://{ATLAS_BUCKET}/outputs/"], capture_output=True, text=True).stdout
        cp_files = sorted([l.split()[-1] for l in ls.split("\n") if "arch_checkpoint_" in l and l.endswith(".csv")])
        if cp_files:
            latest = cp_files[-1]
            s3_cp(f"s3://{ATLAS_BUCKET}/outputs/{latest}", cache / latest)
            arch = pd.read_csv(cache / latest, low_memory=False)
            print(f"  Using {latest}: {len(arch):,} rows (partial)")
        else:
            arch = pd.DataFrame()
    report["arch_rows"] = len(arch)

    # ========== 5. HEEDB heart features ==========
    print("\n=== 5. HEEDB ECG (heart features for paired patients) ===")
    heedb_local = cache / "heedb_final.csv"
    s3_cp(f"s3://{P11_BUCKET}/heedb_outputs/heedb_final.csv", heedb_local)
    heedb = pd.read_csv(heedb_local, low_memory=False)
    print(f"  HEEDB: {len(heedb):,} rows × {heedb.shape[1]} cols")
    report["heedb_rows"] = len(heedb)

    # ========== 6. HEDB EEG features ==========
    print("\n=== 6. HEDB EEG (brain features for paired patients) ===")
    hedb_local = cache / "hedb_final.csv"
    try:
        s3_cp(f"s3://{P11_BUCKET}/hedb_outputs/hedb_final.csv", hedb_local)
        hedb = pd.read_csv(hedb_local, low_memory=False)
        # If size is 1 byte, it's the stale empty placeholder — use latest checkpoint instead
        if len(hedb) < 10:
            ls = subprocess.run(["aws", "s3", "ls", f"s3://{P11_BUCKET}/hedb_outputs/"], capture_output=True, text=True).stdout
            cp_files = sorted([l.split()[-1] for l in ls.split("\n") if "hedb_checkpoint_" in l and l.endswith(".csv")])
            if cp_files:
                latest = cp_files[-1]
                s3_cp(f"s3://{P11_BUCKET}/hedb_outputs/{latest}", cache / latest)
                hedb = pd.read_csv(cache / latest, low_memory=False)
                print(f"  Using {latest}: {len(hedb):,} rows (partial pending v3)")
        else:
            print(f"  HEDB final: {len(hedb):,} rows × {hedb.shape[1]} cols")
    except Exception as e:
        print(f"  HEDB skip: {e}")
        hedb = pd.DataFrame()
    report["hedb_rows"] = len(hedb)

    # ========== 7. OMOP drug + conditions summaries ==========
    print("\n=== 7. OMOP drug + condition summaries ===")
    drug_summary_local = cache / "omop_drug_summary.csv"
    cond_summary_local = cache / "omop_cond_summary.csv"
    s3_cp(f"s3://{ATLAS_BUCKET}/outputs/omop_drug_exposure_summary.csv", drug_summary_local)
    s3_cp(f"s3://{ATLAS_BUCKET}/outputs/omop_condition_occurrence_summary.csv", cond_summary_local)
    drug_summary = pd.read_csv(drug_summary_local, low_memory=False)
    cond_summary = pd.read_csv(cond_summary_local, low_memory=False)
    drug_summary = drug_summary.rename(columns={"person_id":"BDSPPatientID","n_records":"n_drug_records","n_unique_concepts":"n_unique_drugs"})
    cond_summary = cond_summary.rename(columns={"person_id":"BDSPPatientID","n_records":"n_condition_records","n_unique_concepts":"n_unique_conditions"})
    print(f"  drug summary: {len(drug_summary):,}, cond summary: {len(cond_summary):,}")
    report["omop_drug_summary_rows"] = len(drug_summary)
    report["omop_cond_summary_rows"] = len(cond_summary)

    # ========== 8. Build adult master cohort ==========
    print("\n=== 8. Building adult master cohort ===")
    pieces = [v6]
    if len(rec) > 0:
        # Apply same column-set filtering as v6 (drop _warn, drop spo2_nadir)
        rec = rec.drop(columns=[c for c in rec.columns if c.endswith("_warn") or c == "spo2_nadir"], errors="ignore")
        pieces.append(rec)
    master = pd.concat(pieces, ignore_index=True).drop_duplicates(subset=["BDSPPatientID"], keep="last")
    print(f"  v6 + recovery: {len(master):,} unique patients")

    # Left-join arch pass
    if len(arch) > 0:
        arch_keep = arch[[c for c in arch.columns if c not in ["site_id"]]]  # site_id already in master
        master = master.merge(arch_keep, on="BDSPPatientID", how="left")
        print(f"  + arch pass: still {len(master):,} (left join)")

    # Left-join OMOP summaries
    master = master.merge(drug_summary, on="BDSPPatientID", how="left")
    master = master.merge(cond_summary, on="BDSPPatientID", how="left")
    master["n_drug_records"] = master["n_drug_records"].fillna(0).astype(int)
    master["n_unique_drugs"] = master["n_unique_drugs"].fillna(0).astype(int)
    master["n_condition_records"] = master["n_condition_records"].fillna(0).astype(int)
    master["n_unique_conditions"] = master["n_unique_conditions"].fillna(0).astype(int)
    print(f"  + OMOP summaries joined")

    # Left-join HEEDB heart features (paired patients only)
    if len(heedb) > 0:
        heedb["BDSPPatientID"] = pd.to_numeric(heedb["BDSPPatientID"], errors="coerce").astype("Int64")
        heedb_keep = heedb.add_prefix("heedb_").rename(columns={"heedb_BDSPPatientID":"BDSPPatientID"})
        master = master.merge(heedb_keep, on="BDSPPatientID", how="left")
        print(f"  + HEEDB heart features (left join, {heedb['BDSPPatientID'].notna().sum():,} paired)")

    # Left-join HEDB EEG features (paired patients only)
    if len(hedb) > 0:
        hedb["BDSPPatientID"] = pd.to_numeric(hedb["BDSPPatientID"], errors="coerce").astype("Int64")
        hedb_keep = hedb.add_prefix("hedb_").rename(columns={"hedb_BDSPPatientID":"BDSPPatientID"})
        master = master.merge(hedb_keep, on="BDSPPatientID", how="left")
        print(f"  + HEDB EEG features (left join, {hedb['BDSPPatientID'].notna().sum():,} paired)")

    print(f"\n  ADULT MASTER: {len(master):,} rows × {master.shape[1]} cols")
    report["adult_master_rows"] = len(master)
    report["adult_master_cols"] = master.shape[1]

    # ========== 9. Pediatric cohort (BCH) ==========
    print("\n=== 9. BCH pediatric cohort ===")
    if len(bch) > 0:
        bch_clean = bch[bch["error"].isna()] if "error" in bch.columns else bch
        bch_clean.to_csv(OUT / "aggregated_features_v5_pediatric.csv", index=False)
        print(f"  Saved {len(bch_clean):,} pediatric rows to aggregated_features_v5_pediatric.csv")
        report["pediatric_rows"] = len(bch_clean)
    else:
        report["pediatric_rows"] = 0

    # ========== 10. Save master + integrity report ==========
    master_path = OUT / "aggregated_features_v5.csv"
    master.to_csv(master_path, index=False)
    print(f"\n=== 10. Saved master ===")
    print(f"  {master_path} ({master_path.stat().st_size / 1e6:.0f} MB)")

    report["finished"] = time.strftime("%FT%T")
    report["total_runtime_sec"] = round(time.time() - started, 1)
    report_path = OUT / "merge_v5_integrity.json"
    report_path.write_text(json.dumps(report, indent=2, default=str))
    print(f"\n  Integrity: {report_path}")

    print("\n" + "=" * 60)
    print("MERGE V5 SUMMARY")
    print("=" * 60)
    print(f"  Adult master:    {report['adult_master_rows']:,} patients × {report['adult_master_cols']} columns")
    print(f"  Pediatric (BCH): {report.get('pediatric_rows', 0):,}")
    print(f"  HEEDB paired:    {report['heedb_rows']:,}")
    print(f"  HEDB paired:     {report['hedb_rows']:,}")
    print(f"  Runtime:         {report['total_runtime_sec']}s")

if __name__ == "__main__":
    main()
