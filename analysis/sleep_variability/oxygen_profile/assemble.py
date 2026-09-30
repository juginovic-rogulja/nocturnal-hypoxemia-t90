"""
Pull the fleet's JSONL shards down from S3 and flatten them into one per-patient CSV. Numbers
only: the shards contain no signal data, which is what the licence requires. Adapted from
t_ladder_up/assemble.py.

Fleet-wide validation, the same three axes the pilot gate used, now on all 19,173:
  * rec_mean against frozen spo2_mean, and prov_t90 against frozen spo2_pct_below_90
    (expect the documented 93-96 percent within 0.5, the signal-tail stratum)
  * the decile partition identities, with ZERO tolerance for a violation
  * the stage count identity, and stage minutes against the 2026-08-21 ladder run
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
LADDER = os.path.join(os.path.dirname(HERE), "t_ladder_up", "ladder_per_patient.csv")
MODE = sys.argv[1] if len(sys.argv) > 1 else "fleet"
BUCKET = paths.S3_OUTPUT_BUCKET
PREFIX = f"outputs/oxygen_profile/{MODE}"
RAW = os.path.join(HERE, f"_raw_{MODE}")
STAGES = ("wake", "n1", "n2", "n3", "rem")
DEC = [f"{i + 1:02d}" for i in range(10)]
os.makedirs(RAW, exist_ok=True)

if MODE == "fleet_v8_2":   # v8.2 (2026-09-15): the collapse-masked shards are placed under _raw_fleet_v8_2 by X4_groupP/assemble_groupP.py;
    # the S3 copy (outputs/groupP_v8_2/oxy) is the unmasked fleet output and is never synced here
    assert glob.glob(os.path.join(RAW, "*.jsonl")), f"no masked shards under {RAW}: run assemble_groupP.py oxy first"
else:
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
d = d.drop(columns=[c for c in ("h5_structure", "trace") if c in d.columns])
d["BDSPPatientID"] = d.BDSPPatientID.astype(int)
d["_ok"] = (d.status == "ok").astype(int)
d = (d.sort_values("_ok").drop_duplicates("BDSPPatientID", keep="last")
       .drop(columns="_ok").reset_index(drop=True))
print(f"unique patients {len(d):,}")
print(d.status.value_counts().to_dict())

ref = pd.read_csv(os.path.join(HERE, f"frozen_reference_{MODE}.csv") if os.path.exists(os.path.join(HERE, f"frozen_reference_{MODE}.csv")) else os.path.join(HERE, "frozen_reference.csv"))   # v8.2: the v8 reference for mode fleet_v8_2
ref["BDSPPatientID"] = ref.BDSPPatientID.astype(int)
missing = set(ref.BDSPPatientID) - set(d.BDSPPatientID)
print(f"cohort patients with no record at all: {len(missing):,}")

prev = pd.read_csv(LADDER, low_memory=False)
prev = prev[prev.status == "ok"][["BDSPPatientID"] + [f"{s}_min" for s in STAGES]].rename(
    columns={f"{s}_min": f"{s}_min_prev" for s in STAGES})
prev["BDSPPatientID"] = prev.BDSPPatientID.astype(int)

j = d.merge(ref[["BDSPPatientID", "spo2_mean", "spo2_pct_below_90", "spo2_pct_below_88",
                 "TST_min", "N1_pct", "N2_pct", "N3_pct", "REM_pct", "AHI"]],
            on="BDSPPatientID", how="left")
j = j.merge(prev, on="BDSPPatientID", how="left")
ok = j[j.status == "ok"].copy()

print("\nFLEET-WIDE VALIDATION against the frozen columns")
val = {}
for got, froz, lab, tol in (("rec_mean", "spo2_mean", "mean SpO2", 0.5),
                            ("prov_t90", "spo2_pct_below_90", "T90", 0.5)):
    v = (ok[got] - ok[froz]).dropna()
    within = int((v.abs() <= tol).sum())
    rc = float(ok[[got, froz]].corr(method="spearman").iloc[0, 1])
    val[lab] = {"n": int(len(v)), "within": within, "pct": 100.0 * within / len(v),
                "med": float(v.abs().median()), "spearman": rc}
    print(f"  {lab}: within {tol} on {within:,}/{len(v):,} ({100*within/len(v):.1f}%), "
          f"median |diff| {v.abs().median():.2e}, rank corr {rc:.4f}")

print("\nDECILE PARTITION, zero tolerance")
# a record with no usable oximetry at all carries no whole-recording mean and no profile, so
# it has nothing to check. Those are counted, not silently folded into a max().
prof = ok[ok.rec_n.notna() & ok.n_deciles_with_data.notna()]
print(f"  records with a profile to check: {len(prof):,} of {len(ok):,} ok "
      f"({len(ok) - len(prof):,} carry no usable oximetry)")
cnt = prof[[f"dec_n_{t}" for t in DEC]].fillna(0)
mns = prof[[f"dec_mean_{t}" for t in DEC]]
mins = prof[[f"dec_min_{t}" for t in DEC]]
sum_n = cnt.sum(axis=1)
d_count = float((sum_n - prof.rec_n).abs().max())
wsum = (cnt.values * np.nan_to_num(mns.values)).sum(axis=1)
wmean = np.where(sum_n.values > 0, wsum / np.maximum(sum_n.values, 1), np.nan)
d_wmean = float(np.nanmax(np.abs(wmean - prof.rec_mean.values)))
d_span = float((mins.sum(axis=1) - prof.spo2_dur_min).abs().max())
n_all_ten = int((prof.n_deciles_with_data == 10).sum())
in_range = bool(mns.apply(lambda c: c.between(50, 100) | c.isna()).all().all())
n_dec_any = int(len(prof))
print(f"  counts sum to rec_n:                 max |diff| {d_count:.3e}")
print(f"  weighted decile mean equals rec_mean: max |diff| {d_wmean:.3e}")
print(f"  spans cover the SpO2 duration:        max |diff| {d_span:.3e} min")
print(f"  all ten deciles carry data on {n_all_ten:,} of {n_dec_any:,} records with a profile")
print(f"  every decile mean inside [50,100]: {in_range}")

print("\nSTAGE IDENTITY and reproduction of the ladder run")
st = ok[ok.staged_valid_n.notna()].copy()
for s in STAGES:
    for c in (f"{s}_n", f"{s}_min"):
        st[c] = st[c].fillna(0.0)
d_stage_n = float((st[[f"{s}_n" for s in STAGES]].sum(axis=1) - st.staged_valid_n).abs().max())
d_stage_min = float((st[[f"{s}_min" for s in STAGES]].sum(axis=1) - st.staged_min).abs().max())
d_tst = float((st[["n1_min", "n2_min", "n3_min", "rem_min"]].sum(axis=1) - st.tst_min).abs().max())
print(f"  five stage counts sum to staged_valid_n: max |diff| {d_stage_n:.3e} on {len(st):,}")
print(f"  five stage minutes sum to staged_min:    max |diff| {d_stage_min:.3e} min")
print(f"  four sleep stages sum to tst_min:        max |diff| {d_tst:.3e} min")
n_pairs, n_within, worst = 0, 0, 0.0
for s in STAGES:
    c = ok.dropna(subset=[f"{s}_min", f"{s}_min_prev"])
    dd = (c[f"{s}_min"] - c[f"{s}_min_prev"]).abs()
    n_pairs += len(dd)
    n_within += int((dd <= 0.01).sum())
    worst = max(worst, float(dd.max()) if len(dd) else 0.0)
    print(f"    {s:<5} both present {len(dd):,}   within 0.01 min {int((dd<=0.01).sum()):,}"
          f"   max |diff| {float(dd.max()) if len(dd) else float('nan'):.3e}")
print(f"  ladder reproduction overall: {n_within:,}/{n_pairs:,} within 0.01 min")

print("\n  stage percentages against the frozen columns (reported)")
stage_pct = {}
for s, pc in (("n1", "N1_pct"), ("n2", "N2_pct"), ("n3", "N3_pct"), ("rem", "REM_pct")):
    got = ok[f"{s}_min"] / ok.tst_min * 100.0
    dd = (got - ok[pc]).abs().dropna()
    stage_pct[s] = {"n": int(len(dd)), "within_2pp": int((dd <= 2).sum()),
                    "median_abs_diff_pp": float(dd.median())}
    print(f"    {s:<5} median |diff| {dd.median():7.4f} pp   within 2 pp "
          f"{int((dd<=2).sum()):,}/{len(dd):,}")

print("\nCOVERAGE and headline distributions")
cov = {"decile_profile": n_all_ten,
       "any_stage_mean": int(ok[[f"{s}_mean" for s in STAGES]].notna().any(axis=1).sum())}
for s in STAGES:
    cov[f"{s}_mean"] = int(ok[f"{s}_mean"].notna().sum())
    v = ok[f"{s}_mean"].dropna()
    print(f"  {s:<5} mean SpO2 n {len(v):,}  median {v.median():6.2f}  "
          f"IQR {v.quantile(.25):6.2f} to {v.quantile(.75):6.2f}   "
          f"minutes median {ok[f'{s}_min'].median():6.1f}")
print("  decile medians across the cohort:")
dec_prof = [float(ok[f"dec_mean_{t}"].median()) for t in DEC]
print("   " + "  ".join(f"{v:.3f}" for v in dec_prof))

j.to_csv(os.path.join(HERE, "oxyprofile_per_patient.csv"), index=False)
print(f"\nwrote oxyprofile_per_patient.csv  {j.shape}")

summ = {
    "cohort_n": int(len(ref)),
    "n_records": int(len(d)),
    "n_ok": int(len(ok)),
    "n_missing_entirely": int(len(missing)),
    "status_counts": {str(k): int(v) for k, v in j.status.value_counts().items()},
    "mean_within_0p5": val["mean SpO2"]["within"],
    "mean_within_0p5_pct": val["mean SpO2"]["pct"],
    "mean_median_abs_diff": val["mean SpO2"]["med"],
    "mean_spearman": val["mean SpO2"]["spearman"],
    "t90_within_0p5pp": val["T90"]["within"],
    "t90_within_0p5pp_pct": val["T90"]["pct"],
    "t90_median_abs_diff_pp": val["T90"]["med"],
    "t90_spearman": val["T90"]["spearman"],
    "decile_count_max_abs_dev": d_count,
    "decile_weighted_mean_max_abs_dev": d_wmean,
    "decile_span_max_abs_dev_min": d_span,
    "records_with_all_ten_deciles": n_all_ten,
    "records_with_a_profile": int(len(prof)),
    "decile_means_in_range": in_range,
    "stage_count_max_abs_dev": d_stage_n,
    "stage_minutes_sum_max_abs_dev_min": d_stage_min,
    "tst_sum_max_abs_dev_min": d_tst,
    "stage_minutes_vs_ladder_pairs": int(n_pairs),
    "stage_minutes_vs_ladder_within_0p01min": int(n_within),
    "stage_minutes_vs_ladder_max_abs_dev_min": worst,
    "stage_pct_vs_frozen": stage_pct,
    "coverage": cov,
    "decile_median_profile": dec_prof,
    "stage_mean_spo2_median": {s: float(ok[f"{s}_mean"].median()) for s in STAGES},
}
json.dump(summ, open(os.path.join(HERE, "extraction_summary.json"), "w"), indent=2)
print("wrote extraction_summary.json")
