"""
Pull the fleet's JSONL shards down from S3 and flatten them into one per-patient
CSV. Numbers only: the shards contain no signal data, which is what the licence
requires. Adapted from t85_analysis/assemble.py for the six-rung ladder + A95.

Fleet-wide validation:
  * prov_t90 / prov_t88 against their frozen twins (expect the documented
    93-96 percent within 0.5 pp, the signal-tail stratum)
  * t80<=t85<=t88<=t90<=t92<=t95 with ZERO violations allowed
  * prov_t85 against the 2026-08-21 T85 run (same code, same data)
  * prov_a95 + prov_t95 = 100 on every record (independent emission)
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
T85CSV = os.path.join(os.path.dirname(HERE), "t85_analysis", "t85_per_patient.csv")
MODE = sys.argv[1] if len(sys.argv) > 1 else "fleet"
BUCKET = paths.S3_OUTPUT_BUCKET
PREFIX = f"outputs/t_ladder_up/{MODE}"
RAW = os.path.join(HERE, f"_raw_{MODE}")
LADDER = ["prov_t80", "prov_t85", "prov_t88", "prov_t90", "prov_t92", "prov_t95"]
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

prev = pd.read_csv(T85CSV, low_memory=False)
prev = prev[prev.status == "ok"][["BDSPPatientID", "prov_t85"]].rename(
    columns={"prov_t85": "t85_prev_run"})
prev["BDSPPatientID"] = prev.BDSPPatientID.astype(int)

j = d.merge(ref[["BDSPPatientID", "spo2_pct_below_90", "spo2_pct_below_88",
                 "TST_min", "AHI"]], on="BDSPPatientID", how="left")
j = j.merge(prev, on="BDSPPatientID", how="left")
ok = j[j.status == "ok"]

print("\nFLEET-WIDE VALIDATION against the two frozen columns")
val = {}
for prov, froz, lab in (("prov_t90", "spo2_pct_below_90", "T90"),
                        ("prov_t88", "spo2_pct_below_88", "T88")):
    j[f"{lab.lower()}_diff_pp"] = j[prov] - j[froz]
    v = (ok[prov] - ok[froz]).dropna()
    within = int((v.abs() <= 0.5).sum())
    rc = float(ok[[prov, froz]].corr(method="spearman").iloc[0, 1])
    val[lab] = {"n": len(v), "within": within, "pct": 100 * within / len(v),
                "med": float(v.abs().median()), "spearman": rc}
    print(f"  {lab}: within 0.5 pp on {within:,}/{len(v):,} ({100*within/len(v):.1f}%), "
          f"median |diff| {v.abs().median():.2e} pp, rank corr {rc:.4f}")

print("\nT85 against the 2026-08-21 T85 run (same code, same data)")
c85 = ok.dropna(subset=["prov_t85", "t85_prev_run"])
d85 = (c85.prov_t85 - c85.t85_prev_run).abs()
print(f"  both present on {len(c85):,}   within 0.01 pp {int((d85<=0.01).sum()):,} "
      f"({100*(d85<=0.01).mean():.2f}%)   max |diff| {d85.max():.3e} pp")

print("\nLADDER coverage, distribution and monotonicity")
for c in LADDER + ["prov_a95"]:
    v = ok[c].dropna()
    print(f"  {c:<9} n {len(v):,}  median {v.median():8.4f}  p90 {v.quantile(.90):8.3f} "
          f" p99 {v.quantile(.99):7.2f}  max {v.max():6.2f}  exact0 {100*(v==0).mean():.2f}%")
m = ok.dropna(subset=LADDER)
bad = 0
worst_pair, worst_v = "", 0.0
for a, b in zip(LADDER[:-1], LADDER[1:]):
    excess = m[a] - m[b]
    bad += int((excess > 1e-9).sum())
    if excess.max() > worst_v:
        worst_v, worst_pair = float(excess.max()), f"{a}>{b}"
print(f"  monotonicity violated on {bad} of {len(m):,} records "
      f"(worst pairwise excess {worst_v:.3e} pp at {worst_pair})")

mi = ok.dropna(subset=["prov_a95", "prov_t95"])
ident = (mi.prov_a95 + mi.prov_t95 - 100.0).abs()
print(f"  a95 complement identity: max |a95+t95-100| = {ident.max():.3e} pp "
      f"on {len(mi):,} records")

print("\n  Spearman correlations with frozen T90")
for c in LADDER:
    rc = float(ok[[c, "spo2_pct_below_90"]].corr(method="spearman").iloc[0, 1])
    print(f"    {c:<9} {rc:.4f}")

j.to_csv(os.path.join(HERE, "ladder_per_patient.csv"), index=False)
print(f"\nwrote ladder_per_patient.csv  {j.shape}")

summ = {
    "cohort_n": int(len(ref)),
    "n_records": int(len(d)),
    "n_ok": int(len(ok)),
    "n_missing_entirely": int(len(missing)),
    "status_counts": {str(k): int(v) for k, v in j.status.value_counts().items()},
    "coverage": {c: int(ok[c].notna().sum()) for c in LADDER + ["prov_a95"]},
    "t90_within_0p5pp": val["T90"]["within"], "t90_within_0p5pp_pct": val["T90"]["pct"],
    "t88_within_0p5pp": val["T88"]["within"], "t88_within_0p5pp_pct": val["T88"]["pct"],
    "t90_median_abs_diff_pp": val["T90"]["med"], "t88_median_abs_diff_pp": val["T88"]["med"],
    "t90_spearman": val["T90"]["spearman"], "t88_spearman": val["T88"]["spearman"],
    "t85_vs_prev_n": int(len(c85)),
    "t85_vs_prev_within_0p01pp": int((d85 <= 0.01).sum()),
    "t85_vs_prev_max_abs_diff_pp": float(d85.max()) if len(d85) else None,
    "n_monotone_violations": int(bad),
    "worst_mono_excess_pp": worst_v,
    "a95_identity_max_abs_dev_pp": float(ident.max()) if len(mi) else None,
    "medians": {c: float(ok[c].median()) for c in LADDER + ["prov_a95"]},
}
json.dump(summ, open(os.path.join(HERE, "extraction_summary.json"), "w"), indent=2)
print("wrote extraction_summary.json")
