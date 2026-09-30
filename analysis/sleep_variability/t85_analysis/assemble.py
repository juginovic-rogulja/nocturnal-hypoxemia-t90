"""
Pull the fleet's JSONL shards down from S3 and flatten them into one per-patient CSV.
Numbers only: the shards contain no signal data, which is what the licence requires.

Adapted from hypoxic_burden/assemble.py. The fleet-wide validation now checks BOTH
frozen oximetry columns, and the T85 column gets its own coverage and monotonicity
report because it has no frozen twin to check it against.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import glob
import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
MODE = sys.argv[1] if len(sys.argv) > 1 else "fleet"
BUCKET = paths.S3_OUTPUT_BUCKET
PREFIX = f"outputs/t85/{MODE}"
RAW = os.path.join(HERE, f"_raw_{MODE}")
os.makedirs(RAW, exist_ok=True)

subprocess.run(["aws", "s3", "sync", f"s3://{BUCKET}/{PREFIX}/", RAW,
                "--region", "us-east-1", "--exclude", "*", "--include", "*.jsonl"],
               check=True)

rows = []
for f in sorted(glob.glob(os.path.join(RAW, "*.jsonl"))):
    with open(f) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                rows.append(json.loads(line))
            except Exception:
                pass
print(f"raw records {len(rows):,} from {len(glob.glob(os.path.join(RAW, '*.jsonl')))} shards")

d = pd.DataFrame(rows)
d = d.drop(columns=[c for c in ("h5_structure", "trace", "event_dataset_candidates")
                    if c in d.columns])
d["BDSPPatientID"] = d.BDSPPatientID.astype(int)
d["_ok"] = (d.status == "ok").astype(int)
d = (d.sort_values("_ok").drop_duplicates("BDSPPatientID", keep="last")
       .drop(columns="_ok").reset_index(drop=True))
print(f"unique patients {len(d):,}")
print(d.status.value_counts().to_dict())

ref = pd.read_csv(os.path.join(HERE, "frozen_reference.csv"))
ref["BDSPPatientID"] = ref.BDSPPatientID.astype(int)
missing = set(ref.BDSPPatientID) - set(d.BDSPPatientID)
print(f"cohort patients with no record at all: {len(missing):,}")

j = d.merge(ref[["BDSPPatientID", "spo2_pct_below_90", "spo2_pct_below_88",
                 "TST_min", "AHI"]], on="BDSPPatientID", how="left")
ok = j[j.status == "ok"]

print("\nFLEET-WIDE VALIDATION against the two frozen columns")
for prov, froz, lab in (("prov_t90", "spo2_pct_below_90", "T90"),
                        ("prov_t88", "spo2_pct_below_88", "T88")):
    if prov not in j.columns:
        continue
    j[f"{lab.lower()}_diff_pp"] = j[prov] - j[froz]
    v = ok[prov] - ok[froz]
    v = v.dropna()
    print(f"  {lab}: within 0.5 pp on {int((v.abs()<=0.5).sum()):,}/{len(v):,} "
          f"({100*(v.abs()<=0.5).mean():.1f}%), median |diff| {v.abs().median():.2e} pp, "
          f"rank corr {ok[[prov,froz]].corr(method='spearman').iloc[0,1]:.4f}")

print("\nT85 (no frozen twin, so this is coverage and internal consistency)")
t85 = ok.prov_t85.dropna()
print(f"  n with a value {len(t85):,} of {len(ok):,} ok records "
      f"({100*len(t85)/max(len(ok),1):.2f}%)")
print(f"  median {t85.median():.4f}  IQR {t85.quantile(.25):.4f}-{t85.quantile(.75):.4f} "
      f" p90 {t85.quantile(.90):.3f}  p99 {t85.quantile(.99):.2f}  max {t85.max():.2f}")
m = ok.dropna(subset=["prov_t85", "prov_t88", "prov_t90"])
bad = int(((m.prov_t85 > m.prov_t88 + 1e-9) | (m.prov_t88 > m.prov_t90 + 1e-9)).sum())
print(f"  T85 <= T88 <= T90 violated on {bad} of {len(m):,} records")
print(f"  Spearman T85 vs frozen T90 {ok[['prov_t85','spo2_pct_below_90']].corr(method='spearman').iloc[0,1]:.4f}")
print(f"  Spearman T85 vs prov  T90 {ok[['prov_t85','prov_t90']].corr(method='spearman').iloc[0,1]:.4f}")

j.to_csv(os.path.join(HERE, "t85_per_patient.csv"), index=False)
print(f"\nwrote t85_per_patient.csv  {j.shape}")

# ------------------------------------------------------------------ machine summary
cmp_ = ok.dropna(subset=["prov_t90", "prov_t88", "spo2_pct_below_90", "spo2_pct_below_88"])
summ = {
    "cohort_n": int(len(ref)),
    "n_records": int(len(d)),
    "n_ok": int(len(ok)),
    "n_t85": int(ok.prov_t85.notna().sum()),
    "n_missing_entirely": int(len(missing)),
    "status_counts": {str(k): int(v) for k, v in j.status.value_counts().items()},
    "n_valid_cmp": int(len(cmp_)),
    "t90_within_0p5pp": int(((cmp_.prov_t90 - cmp_.spo2_pct_below_90).abs() <= 0.5).sum()),
    "t88_within_0p5pp": int(((cmp_.prov_t88 - cmp_.spo2_pct_below_88).abs() <= 0.5).sum()),
    "t90_within_0p5pp_pct": float(100 * ((cmp_.prov_t90 - cmp_.spo2_pct_below_90).abs() <= 0.5).mean()),
    "t88_within_0p5pp_pct": float(100 * ((cmp_.prov_t88 - cmp_.spo2_pct_below_88).abs() <= 0.5).mean()),
    "t90_median_abs_diff_pp": float((cmp_.prov_t90 - cmp_.spo2_pct_below_90).abs().median()),
    "t88_median_abs_diff_pp": float((cmp_.prov_t88 - cmp_.spo2_pct_below_88).abs().median()),
    "t90_spearman": float(cmp_[["prov_t90", "spo2_pct_below_90"]].corr(method="spearman").iloc[0, 1]),
    "t88_spearman": float(cmp_[["prov_t88", "spo2_pct_below_88"]].corr(method="spearman").iloc[0, 1]),
    "n_monotone_violations": int(bad),
    "t85_median": float(t85.median()),
    "t85_q75": float(t85.quantile(.75)),
    "t85_max": float(t85.max()),
}
json.dump(summ, open(os.path.join(HERE, "extraction_summary.json"), "w"), indent=2)
print("wrote extraction_summary.json")
