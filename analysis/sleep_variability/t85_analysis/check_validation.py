"""
VALIDATION GATE. Nothing goes to the fleet until this passes.

Two frozen oximetry columns exist in the paper's data, spo2_pct_below_90 and
spo2_pct_below_88. T85 has no frozen twin, which is the entire reason this run exists,
so T85 itself cannot be validated directly. What CAN be validated is the rule that
produces it: the same lines of code, on the same cleaned signal and the same
whole-recording denominator, also emit prov_t90 and prov_t88. If both of those
reproduce their frozen twins, the rule is right and prov_t85 is the same rule at a
different threshold.

GATE, scored on the 25 pilot patients (I0002) whose answer the 2026-08-13 raw pilot
already established:

  1. prov_t90 within 0.5 pp of frozen spo2_pct_below_90 on 23 or more of 25
  2. prov_t88 within 0.5 pp of frozen spo2_pct_below_88 on 23 or more of 25
  3. prov_t85 is present, finite, in [0,100], and never exceeds prov_t88 (monotone by
     construction: fewer samples sit below 85 than below 88). A violation would mean the
     threshold logic is wrong, and no frozen column could ever have caught it.

The 15 I0006 records are reported but not gated. The pilot never touched that site, so
they are new information rather than a known answer.
"""
import json
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
J = os.path.join(HERE, "validate_out_0.jsonl")

rows = [json.loads(l) for l in open(J) if l.strip()]
d = pd.DataFrame(rows)
print(f"records {len(d)}")
print(d.status.value_counts().to_dict())
print(d.groupby("site").size().to_dict())

ref = pd.read_csv(os.path.join(HERE, "frozen_reference.csv"))
ref["BDSPPatientID"] = ref.BDSPPatientID.astype(int)
d["BDSPPatientID"] = d.BDSPPatientID.astype(int)
pilot = set(pd.read_csv(os.path.join(HERE, "pilot_ids.csv")).BDSPPatientID.astype(int))

j = d.merge(ref[["BDSPPatientID", "site", "spo2_pct_below_90", "spo2_pct_below_88",
                 "TST_min", "recording_dur_min", "AHI"]],
            on="BDSPPatientID", how="left", suffixes=("", "_frozen"))
j["d90"] = j.prov_t90 - j.spo2_pct_below_90
j["d88"] = j.prov_t88 - j.spo2_pct_below_88
j["is_pilot"] = j.BDSPPatientID.isin(pilot)
ok = j[j.status == "ok"].copy()
p = ok[ok.is_pilot]

print("\n" + "=" * 78)
print(f"GATE, on the {len(p)} pilot patients")
print("=" * 78)
n90 = int((p.d90.abs() <= 0.5).sum())
n88 = int((p.d88.abs() <= 0.5).sum())
print(f"  frozen T90 reproduced within 0.5 pp on {n90}/{len(p)}   "
      f"median |diff| {p.d90.abs().median():.3e} pp")
print(f"  frozen T88 reproduced within 0.5 pp on {n88}/{len(p)}   "
      f"median |diff| {p.d88.abs().median():.3e} pp")

t85 = ok.prov_t85
mono = bool((ok.prov_t85 <= ok.prov_t88 + 1e-9).all() and
            (ok.prov_t88 <= ok.prov_t90 + 1e-9).all())
rng_ok = bool(t85.notna().all() and np.isfinite(t85).all() and
              (t85 >= 0).all() and (t85 <= 100).all())
print(f"\n  T85 present and finite in [0,100] on all {len(ok)} ok records: {rng_ok}")
print(f"  T85 <= T88 <= T90 on every record: {mono}")
print(f"  T85 median {t85.median():.4f}   IQR {t85.quantile(.25):.4f}-"
      f"{t85.quantile(.75):.4f}   max {t85.max():.2f}")

gate1 = n90 >= 23
gate2 = n88 >= 23
gate3 = mono and rng_ok

print("\n  per-site reproduction (I0006 is reported, not gated)")
for s, g in ok.groupby("site"):
    print(f"  {s}  n={len(g):2d}  T90 within 0.5pp {int((g.d90.abs()<=0.5).sum()):2d}/{len(g):2d}   "
          f"T88 within 0.5pp {int((g.d88.abs()<=0.5).sum()):2d}/{len(g):2d}   "
          f"T90 rank corr {g[['prov_t90','spo2_pct_below_90']].corr(method='spearman').iloc[0,1]:.4f}")

print("\n  worst 6 T90 differences among the pilot")
print(p.reindex(p.d90.abs().sort_values(ascending=False).index)
      [["BDSPPatientID", "site", "prov_t90", "spo2_pct_below_90", "d90",
        "prov_t88", "spo2_pct_below_88", "d88", "prov_t85"]]
      .head(6).to_string(index=False, float_format=lambda x: f"{x:10.4f}"))

print("\n  every pilot record, sorted by frozen T90")
print(p.sort_values("spo2_pct_below_90")
      [["BDSPPatientID", "spo2_pct_below_90", "prov_t90", "spo2_pct_below_88",
        "prov_t88", "prov_t85", "prov_n_valid"]]
      .to_string(index=False, float_format=lambda x: f"{x:11.4f}"))

print("\n  the 15 I0006 records, new ground")
i6 = ok[~ok.is_pilot]
if len(i6):
    print(i6.sort_values("spo2_pct_below_90")
          [["BDSPPatientID", "spo2_pct_below_90", "prov_t90", "spo2_pct_below_88",
            "prov_t88", "prov_t85"]]
          .to_string(index=False, float_format=lambda x: f"{x:11.4f}"))

print("\n" + "=" * 78)
print(f"GATE 1 (frozen T90 reproduces, >=23/25): {'PASS' if gate1 else 'FAIL'}  [{n90}/{len(p)}]")
print(f"GATE 2 (frozen T88 reproduces, >=23/25): {'PASS' if gate2 else 'FAIL'}  [{n88}/{len(p)}]")
print(f"GATE 3 (T85 sane and monotone):          {'PASS' if gate3 else 'FAIL'}")
print(f"VERDICT: {'PROCEED TO FLEET' if (gate1 and gate2 and gate3) else 'STOP'}")
print("=" * 78)

json.dump({"n_records": int(len(d)), "n_ok": int(len(ok)), "n_pilot": int(len(p)),
           "t90_within_0p5pp": n90, "t88_within_0p5pp": n88,
           "t90_median_abs_diff_pp": float(p.d90.abs().median()),
           "t88_median_abs_diff_pp": float(p.d88.abs().median()),
           "t85_monotone": mono, "t85_in_range": rng_ok,
           "gate1": gate1, "gate2": gate2, "gate3": gate3,
           "verdict": "PROCEED" if (gate1 and gate2 and gate3) else "STOP"},
          open(os.path.join(HERE, "validation_gate.json"), "w"), indent=2)
ok.drop(columns=[c for c in ("h5_structure", "trace", "event_dataset_candidates")
                 if c in ok.columns]).to_csv(
    os.path.join(HERE, "validation_detail.csv"), index=False)
