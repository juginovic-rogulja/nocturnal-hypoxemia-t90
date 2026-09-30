"""
VALIDATION GATE for the overnight-oxygen-profile run. Nothing goes to the fleet until this
passes. Same shape as t_ladder_up/check_validation.py, with the gates that this run's two new
deliverables actually admit.

No frozen column exists for a decile mean or for saturation within a stage, so neither can be
validated against the parquet directly. Three things can be, and together they cover every step
the deliverables sit on.

  THE SIGNAL PATH. The same lines that produce the decile means also produce rec_mean, the
  whole-recording mean of the cleaned signal, whose frozen twin is spo2_mean, and prov_t90,
  which is the frozen T90 rule verbatim and whose twin is spo2_pct_below_90. If the signal is
  being read, mapped and cleaned correctly then these two reproduce, and the decile means are
  the same samples in ten pieces.

  THE PARTITION. The ten spans partition the recording exactly, so the decile counts must sum
  to the whole-recording valid count and the count-weighted mean of the ten means must equal
  the whole-recording mean. Neither is true by construction of the analysis code: the means and
  counts are computed span by span and the whole-recording figures come from a different line.
  A boundary error, a double-counted sample or a dropped span all break these.

  THE STAGING PATH. The five stage counts must sum to staged_valid_n, which the worker emits
  independently. Stage minutes must reproduce the 2026-08-21 threshold-ladder run's own
  n1_min..rem_min exactly, since that is the same code on the same data, and must agree with
  the frozen N1_pct/N2_pct/N3_pct/REM_pct times TST_min.

GATES, scored on the 25 pilot patients (I0002) whose answer the 2026-08-13 raw pilot
established, except gates 3, 4 and 5 which use every ok record.

  1. rec_mean within 0.5 of frozen spo2_mean on 23 or more of 25
  2. prov_t90 within 0.5 pp of frozen spo2_pct_below_90 on 23 or more of 25
  3. decile partition identities hold on every ok record
  4. stage count identity and stage ranges hold on every ok record with staging
  5. stage minutes reproduce the ladder run within 0.01 min on every record with both
"""
import json
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
LADDER = os.path.join(os.path.dirname(HERE), "t_ladder_up", "ladder_per_patient.csv")
J = os.path.join(HERE, "validate_out_0.jsonl")
STAGES = ("wake", "n1", "n2", "n3", "rem")
DEC = [f"{i + 1:02d}" for i in range(10)]

rows = [json.loads(l) for l in open(J) if l.strip()]
d = pd.DataFrame(rows)
print(f"records {len(d)}")
print(d.status.value_counts().to_dict())
print(d.groupby("site").size().to_dict())

ref = pd.read_csv(os.path.join(HERE, "frozen_reference.csv"))
ref["BDSPPatientID"] = ref.BDSPPatientID.astype(int)
d["BDSPPatientID"] = d.BDSPPatientID.astype(int)
pilot = set(pd.read_csv(os.path.join(HERE, "pilot_ids.csv")).BDSPPatientID.astype(int))

prev = pd.read_csv(LADDER, low_memory=False)
prev = prev[prev.status == "ok"][["BDSPPatientID"] + [f"{s}_min" for s in STAGES]].rename(
    columns={f"{s}_min": f"{s}_min_prev" for s in STAGES})
prev["BDSPPatientID"] = prev.BDSPPatientID.astype(int)

j = d.merge(ref[["BDSPPatientID", "site", "spo2_mean", "spo2_pct_below_90", "TST_min",
                 "N1_pct", "N2_pct", "N3_pct", "REM_pct", "AHI"]],
            on="BDSPPatientID", how="left", suffixes=("", "_frozen"))
j = j.merge(prev, on="BDSPPatientID", how="left")
j["dmean"] = j.rec_mean - j.spo2_mean
j["d90"] = j.prov_t90 - j.spo2_pct_below_90
j["is_pilot"] = j.BDSPPatientID.isin(pilot)
ok = j[j.status == "ok"].copy()
p = ok[ok.is_pilot]

print("\n" + "=" * 82)
print(f"GATES 1-2, on the {len(p)} pilot patients")
print("=" * 82)
n_mean = int((p.dmean.abs() <= 0.5).sum())
n90 = int((p.d90.abs() <= 0.5).sum())
print(f"  frozen mean SpO2 reproduced within 0.5 on {n_mean}/{len(p)}   "
      f"median |diff| {p.dmean.abs().median():.3e}   tight (<=0.05) {int((p.dmean.abs()<=0.05).sum())}")
print(f"  frozen T90 reproduced within 0.5 pp on {n90}/{len(p)}   "
      f"median |diff| {p.d90.abs().median():.3e} pp")

print("\nGATE 3, the decile partition, on all ok records")
cnt = ok[[f"dec_n_{t}" for t in DEC]].fillna(0)
mns = ok[[f"dec_mean_{t}" for t in DEC]]
mins = ok[[f"dec_min_{t}" for t in DEC]]
sum_n = cnt.sum(axis=1)
d_count = (sum_n - ok.rec_n).abs()
wsum = (cnt.values * np.nan_to_num(mns.values)).sum(axis=1)
wmean = np.where(sum_n.values > 0, wsum / np.maximum(sum_n.values, 1), np.nan)
d_wmean = np.abs(wmean - ok.rec_mean.values)
d_span = (mins.sum(axis=1) - ok.spo2_dur_min).abs()
all_ten = (ok.n_deciles_with_data == 10)
in_range = mns.apply(lambda c: c.between(50, 100) | c.isna()).all(axis=1)
print(f"  counts sum to the whole-recording valid count: max |diff| {d_count.max():.3e} "
      f"on {len(ok)} records")
print(f"  count-weighted decile mean equals rec_mean:    max |diff| {np.nanmax(d_wmean):.3e}")
print(f"  the ten spans cover the SpO2 duration:         max |diff| {d_span.max():.3e} min")
print(f"  all ten deciles carry data on {int(all_ten.sum())}/{len(ok)} records")
print(f"  every decile mean inside [50,100]: {bool(in_range.all())}")
# A span that is entirely oximeter dropout has no mean, and the worker reports none rather
# than fabricating one (selftest case 3). That is correct behaviour, so a complete profile is
# a coverage figure, not a pass condition. The figure uses only complete profiles and says so.
print(f"  every record carries at least one decile: {bool((ok.n_deciles_with_data >= 1).all())}")
gate3 = bool(d_count.max() <= 1e-9 and np.nanmax(d_wmean) <= 1e-6
             and d_span.max() <= 1e-6 and in_range.all()
             and (ok.n_deciles_with_data >= 1).all()
             and all_ten.mean() >= 0.85)

print("\nGATE 4, the stage identity, on all ok records with staging")
st = ok.dropna(subset=[f"{s}_n" for s in STAGES], how="all").copy()
for s in STAGES:
    for c in (f"{s}_n", f"{s}_min"):
        st[c] = st[c].fillna(0.0)
d_stage_n = (st[[f"{s}_n" for s in STAGES]].sum(axis=1) - st.staged_valid_n).abs()
d_stage_min = (st[[f"{s}_min" for s in STAGES]].sum(axis=1) - st.staged_min).abs()
d_tst = (st[["n1_min", "n2_min", "n3_min", "rem_min"]].sum(axis=1) - st.tst_min).abs()
smeans = st[[f"{s}_mean" for s in STAGES]]
sm_range = smeans.apply(lambda c: c.between(50, 100) | c.isna()).all(axis=1)
print(f"  five stage counts sum to staged_valid_n:  max |diff| {d_stage_n.max():.3e} "
      f"on {len(st)} records")
print(f"  five stage minutes sum to staged_min:     max |diff| {d_stage_min.max():.3e} min")
print(f"  the four sleep stages sum to tst_min:     max |diff| {d_tst.max():.3e} min")
print(f"  every stage mean inside [50,100]: {bool(sm_range.all())}")
for s in STAGES:
    v = st[f"{s}_mean"].dropna()
    print(f"    {s:<5} mean SpO2 median {v.median():6.2f}  n with a value {len(v):3d}   "
          f"minutes median {st[f'{s}_min'].median():6.1f}")
gate4 = bool(d_stage_n.max() <= 1e-9 and d_stage_min.max() <= 1e-6
             and d_tst.max() <= 1e-6 and sm_range.all() and len(st) >= 38)

print("\nGATE 5, stage minutes against the 2026-08-21 ladder run (same code, same data)")
worst, n_pairs, n_within = 0.0, 0, 0
for s in STAGES:
    c = ok.dropna(subset=[f"{s}_min", f"{s}_min_prev"])
    dd = (c[f"{s}_min"] - c[f"{s}_min_prev"]).abs()
    if len(dd):
        worst = max(worst, float(dd.max()))
        n_pairs += len(dd)
        n_within += int((dd <= 0.01).sum())
    print(f"    {s:<5} records with both {len(dd):3d}   within 0.01 min {int((dd<=0.01).sum()):3d}"
          f"   max |diff| {float(dd.max()) if len(dd) else float('nan'):.3e} min")
gate5 = n_pairs >= 5 * 38 and n_within == n_pairs

print("\n  frozen stage percentages (N_pct times TST_min), reported not gated")
for s, pc in (("n1", "N1_pct"), ("n2", "N2_pct"), ("n3", "N3_pct"), ("rem", "REM_pct")):
    got = ok[f"{s}_min"] / ok.tst_min * 100.0
    dd = (got - ok[pc]).abs().dropna()
    print(f"    {s:<5} median |diff| {dd.median():7.4f} pp   within 2 pp "
          f"{int((dd <= 2).sum())}/{len(dd)}")

print("\n  the ten decile means, median across the 40 validation records")
prof = [float(ok[f"dec_mean_{t}"].median()) for t in DEC]
print("   " + "  ".join(f"{v:.2f}" for v in prof))

print("\n  worst 6 mean-SpO2 differences among the pilot")
print(p.reindex(p.dmean.abs().sort_values(ascending=False).index)
      [["BDSPPatientID", "site", "rec_mean", "spo2_mean", "dmean", "prov_t90",
        "spo2_pct_below_90", "d90"]]
      .head(6).to_string(index=False, float_format=lambda x: f"{x:10.4f}"))

gate1 = n_mean >= 23
gate2 = n90 >= 23

print("\n" + "=" * 82)
print(f"GATE 1 (frozen mean SpO2 reproduces, >=23/25): {'PASS' if gate1 else 'FAIL'}  [{n_mean}/{len(p)}]")
print(f"GATE 2 (frozen T90 reproduces, >=23/25):       {'PASS' if gate2 else 'FAIL'}  [{n90}/{len(p)}]")
print(f"GATE 3 (decile partition identities):          {'PASS' if gate3 else 'FAIL'}")
print(f"GATE 4 (stage identities and ranges):          {'PASS' if gate4 else 'FAIL'}")
print(f"GATE 5 (stage minutes reproduce the ladder):   {'PASS' if gate5 else 'FAIL'}  [{n_within}/{n_pairs}]")
ALL = gate1 and gate2 and gate3 and gate4 and gate5
print(f"VERDICT: {'PROCEED TO FLEET' if ALL else 'STOP'}")
print("=" * 82)

json.dump({"n_records": int(len(d)), "n_ok": int(len(ok)), "n_pilot": int(len(p)),
           "mean_within_0p5": n_mean, "t90_within_0p5pp": n90,
           "mean_median_abs_diff": float(p.dmean.abs().median()),
           "t90_median_abs_diff_pp": float(p.d90.abs().median()),
           "decile_count_max_abs_dev": float(d_count.max()),
           "decile_weighted_mean_max_abs_dev": float(np.nanmax(d_wmean)),
           "decile_span_max_abs_dev_min": float(d_span.max()),
           "records_with_all_ten_deciles": int(all_ten.sum()),
           "stage_count_max_abs_dev": float(d_stage_n.max()),
           "stage_minutes_vs_ladder_pairs": int(n_pairs),
           "stage_minutes_vs_ladder_within_0p01min": int(n_within),
           "stage_minutes_vs_ladder_max_abs_dev_min": worst,
           "decile_median_profile": prof,
           "stage_mean_spo2_median": {s: float(st[f"{s}_mean"].median()) for s in STAGES},
           "gate1": bool(gate1), "gate2": bool(gate2), "gate3": bool(gate3),
           "gate4": bool(gate4), "gate5": bool(gate5),
           "verdict": "PROCEED" if ALL else "STOP"},
          open(os.path.join(HERE, "validation_gate.json"), "w"), indent=2)
ok.drop(columns=[c for c in ("h5_structure", "trace") if c in ok.columns]).to_csv(
    os.path.join(HERE, "validation_detail.csv"), index=False)
print("wrote validation_gate.json, validation_detail.csv")
