"""
VALIDATION GATE for the threshold-ladder run. Nothing goes to the fleet until
this passes. Extends t85_analysis/check_validation.py.

Two frozen oximetry columns exist in the paper's data, spo2_pct_below_90 and
spo2_pct_below_88. T80/T92/T95 have no frozen twin, so they cannot be validated
directly. What CAN be validated is the rule that produces them: the same lines
of code, on the same cleaned signal and the same whole-recording denominator,
also emit prov_t90 and prov_t88 (frozen twins) and prov_t85 (computed by the
2026-08-21 T85 run from the same raw data with the same code, so it must agree
to numerical noise).

GATES, scored on the 25 pilot patients (I0002) whose answer the 2026-08-13 raw
pilot already established, except gate 4 which uses all ok records:

  1. prov_t90 within 0.5 pp of frozen spo2_pct_below_90 on 23 or more of 25
  2. prov_t88 within 0.5 pp of frozen spo2_pct_below_88 on 23 or more of 25
  3. t80 <= t85 <= t88 <= t90 <= t92 <= t95 on every ok record, and all six
     present, finite, in [0,100]. A violation means the threshold logic is
     wrong, and no frozen column could ever have caught it.
  4. prov_t85 agrees with the T85 run's prov_t85 within 0.01 pp on every
     validation record that has both (same code, same data, so any real
     difference is a defect).
  5. prov_a95 + prov_t95 = 100 within 1e-6 pp on every ok record (the
     preserved-saturation complement, emitted independently by the worker).

The 15 I0006 records are reported but not gated on gates 1-2. The pilot never
touched that site, so they are new information rather than a known answer.
"""
import json
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
T85CSV = os.path.join(os.path.dirname(HERE), "t85_analysis", "t85_per_patient.csv")
J = os.path.join(HERE, "validate_out_0.jsonl")
LADDER = ["prov_t80", "prov_t85", "prov_t88", "prov_t90", "prov_t92", "prov_t95"]

rows = [json.loads(l) for l in open(J) if l.strip()]
d = pd.DataFrame(rows)
print(f"records {len(d)}")
print(d.status.value_counts().to_dict())
print(d.groupby("site").size().to_dict())

ref = pd.read_csv(os.path.join(HERE, "frozen_reference.csv"))
ref["BDSPPatientID"] = ref.BDSPPatientID.astype(int)
d["BDSPPatientID"] = d.BDSPPatientID.astype(int)
pilot = set(pd.read_csv(os.path.join(HERE, "pilot_ids.csv")).BDSPPatientID.astype(int))

prev = pd.read_csv(T85CSV, low_memory=False)
prev = prev[prev.status == "ok"][["BDSPPatientID", "prov_t85"]].rename(
    columns={"prov_t85": "t85_prev_run"})
prev["BDSPPatientID"] = prev.BDSPPatientID.astype(int)

j = d.merge(ref[["BDSPPatientID", "site", "spo2_pct_below_90", "spo2_pct_below_88",
                 "TST_min", "recording_dur_min", "AHI"]],
            on="BDSPPatientID", how="left", suffixes=("", "_frozen"))
j = j.merge(prev, on="BDSPPatientID", how="left")
j["d90"] = j.prov_t90 - j.spo2_pct_below_90
j["d88"] = j.prov_t88 - j.spo2_pct_below_88
j["d85_prev"] = j.prov_t85 - j.t85_prev_run
j["is_pilot"] = j.BDSPPatientID.isin(pilot)
ok = j[j.status == "ok"].copy()
p = ok[ok.is_pilot]

print("\n" + "=" * 78)
print(f"GATES 1-2, on the {len(p)} pilot patients")
print("=" * 78)
n90 = int((p.d90.abs() <= 0.5).sum())
n88 = int((p.d88.abs() <= 0.5).sum())
print(f"  frozen T90 reproduced within 0.5 pp on {n90}/{len(p)}   "
      f"median |diff| {p.d90.abs().median():.3e} pp")
print(f"  frozen T88 reproduced within 0.5 pp on {n88}/{len(p)}   "
      f"median |diff| {p.d88.abs().median():.3e} pp")

print("\nGATE 3, ladder sanity on all ok records")
present = all(c in ok.columns for c in LADDER)
if present:
    lv = ok[LADDER]
    rng_ok = bool(lv.notna().all().all() and np.isfinite(lv.values).all()
                  and (lv.values >= 0).all() and (lv.values <= 100).all())
    mono = True
    worst_pair, worst_v = "", 0.0
    for a, b in zip(LADDER[:-1], LADDER[1:]):
        viol = (ok[a] - ok[b]).max()
        if viol > worst_v:
            worst_v, worst_pair = viol, f"{a}>{b}"
        mono = mono and bool((ok[a] <= ok[b] + 1e-9).all())
else:
    rng_ok = mono = False
    worst_pair, worst_v = "columns missing", np.nan
print(f"  all six thresholds present, finite, in [0,100] on all {len(ok)} ok records: {rng_ok}")
print(f"  t80<=t85<=t88<=t90<=t92<=t95 on every record: {mono}   "
      f"(worst pairwise excess {worst_v:.3e} pp at {worst_pair})")
for c in LADDER:
    if c in ok.columns:
        print(f"    {c}  median {ok[c].median():9.4f}   max {ok[c].max():7.2f}")

print("\nGATE 4, prov_t85 against the T85 run (same code, same data)")
c85 = ok.dropna(subset=["prov_t85", "t85_prev_run"])
max85 = float(c85.d85_prev.abs().max()) if len(c85) else np.nan
n85_ok = int((c85.d85_prev.abs() <= 0.01).sum())
print(f"  records with both values: {len(c85)}   within 0.01 pp: {n85_ok}/{len(c85)}   "
      f"max |diff| {max85:.3e} pp")

print("\nGATE 5, the a95 complement identity")
if "prov_a95" in ok.columns:
    ident = (ok.prov_a95 + ok.prov_t95 - 100.0).abs()
    max_id = float(ident.max())
    id_ok = bool((ident <= 1e-6).all())
    print(f"  max |a95 + t95 - 100| over {len(ok)} records: {max_id:.3e} pp   "
          f"all within 1e-6: {id_ok}")
else:
    max_id, id_ok = np.nan, False
    print("  prov_a95 missing")

gate1 = n90 >= 23
gate2 = n88 >= 23
gate3 = mono and rng_ok
gate4 = len(c85) >= 38 and n85_ok == len(c85)
gate5 = id_ok

print("\n  per-site reproduction (I0006 is reported, not gated)")
for s, g in ok.groupby("site"):
    print(f"  {s}  n={len(g):2d}  T90 within 0.5pp {int((g.d90.abs()<=0.5).sum()):2d}/{len(g):2d}   "
          f"T88 within 0.5pp {int((g.d88.abs()<=0.5).sum()):2d}/{len(g):2d}   "
          f"T90 rank corr {g[['prov_t90','spo2_pct_below_90']].corr(method='spearman').iloc[0,1]:.4f}")

print("\n  worst 6 T90 differences among the pilot")
print(p.reindex(p.d90.abs().sort_values(ascending=False).index)
      [["BDSPPatientID", "site", "prov_t90", "spo2_pct_below_90", "d90",
        "prov_t88", "d88", "prov_t85", "d85_prev"]]
      .head(6).to_string(index=False, float_format=lambda x: f"{x:10.4f}"))

print("\n  every pilot record, the full ladder, sorted by frozen T90")
print(p.sort_values("spo2_pct_below_90")
      [["BDSPPatientID", "spo2_pct_below_90", "prov_t80", "prov_t85", "prov_t88",
        "prov_t90", "prov_t92", "prov_t95", "prov_a95"]]
      .to_string(index=False, float_format=lambda x: f"{x:9.4f}"))

print("\n  the 15 I0006 records, new ground")
i6 = ok[~ok.is_pilot]
if len(i6):
    print(i6.sort_values("spo2_pct_below_90")
          [["BDSPPatientID", "spo2_pct_below_90", "prov_t80", "prov_t85", "prov_t88",
            "prov_t90", "prov_t92", "prov_t95"]]
          .to_string(index=False, float_format=lambda x: f"{x:9.4f}"))

print("\n" + "=" * 78)
print(f"GATE 1 (frozen T90 reproduces, >=23/25): {'PASS' if gate1 else 'FAIL'}  [{n90}/{len(p)}]")
print(f"GATE 2 (frozen T88 reproduces, >=23/25): {'PASS' if gate2 else 'FAIL'}  [{n88}/{len(p)}]")
print(f"GATE 3 (six-rung ladder sane, monotone): {'PASS' if gate3 else 'FAIL'}")
print(f"GATE 4 (T85 reproduces the T85 run, <=0.01 pp): {'PASS' if gate4 else 'FAIL'}  [{n85_ok}/{len(c85)}]")
print(f"GATE 5 (a95 complement identity):        {'PASS' if gate5 else 'FAIL'}")
ALL = gate1 and gate2 and gate3 and gate4 and gate5
print(f"VERDICT: {'PROCEED TO FLEET' if ALL else 'STOP'}")
print("=" * 78)

json.dump({"n_records": int(len(d)), "n_ok": int(len(ok)), "n_pilot": int(len(p)),
           "t90_within_0p5pp": n90, "t88_within_0p5pp": n88,
           "t90_median_abs_diff_pp": float(p.d90.abs().median()),
           "t88_median_abs_diff_pp": float(p.d88.abs().median()),
           "ladder_monotone": bool(mono), "ladder_in_range": bool(rng_ok),
           "worst_mono_excess_pp": float(worst_v),
           "t85_vs_prev_n": int(len(c85)), "t85_vs_prev_within_0p01pp": n85_ok,
           "t85_vs_prev_max_abs_diff_pp": max85,
           "a95_identity_max_abs_dev_pp": max_id,
           "gate1": bool(gate1), "gate2": bool(gate2), "gate3": bool(gate3),
           "gate4": bool(gate4), "gate5": bool(gate5),
           "verdict": "PROCEED" if ALL else "STOP"},
          open(os.path.join(HERE, "validation_gate.json"), "w"), indent=2)
ok.drop(columns=[c for c in ("h5_structure", "trace", "event_dataset_candidates")
                 if c in ok.columns]).to_csv(
    os.path.join(HERE, "validation_detail.csv"), index=False)
print("wrote validation_gate.json, validation_detail.csv")
