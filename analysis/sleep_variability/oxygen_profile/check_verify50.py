"""
The 50-patient independent verification gate.

worker_verify50.py re-derived, for 50 patients drawn at random from the full 19,173 cohort with
seed 20260821, the ten per-decile mean saturations and the five per-stage mean saturations and
stage minutes, through a second implementation that shares the definitions with the fleet
worker and none of its code. This compares the two value by value.

TOLERANCES
  means      0.05 percentage points, the gate. Anything above it is listed with its cause.
  minutes    0.05 minutes, reported and gated on the same footing.
  counts     the two implementations cut the ten spans by different arithmetic (rint of a
             linspace against integer division), which can move a single sample across a
             boundary, so a decile count may legitimately differ by one sample per boundary.
             Reported, and gated at 2 samples.

Writes verify50.csv (one row per patient and quantity) and verify50_summary.json.

    python3 check_verify50.py
"""
import json
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
TOL_MEAN_PP = 0.05
TOL_MIN = 0.05
TOL_COUNT = 2
DEC = [f"{i + 1:02d}" for i in range(10)]
STAGES = ("wake", "n1", "n2", "n3", "rem")

draw = json.load(open(os.path.join(HERE, "verify50_draw.json")))
rows = [json.loads(l) for l in open(os.path.join(HERE, "verify50_out_0.jsonl")) if l.strip()]
v = pd.DataFrame(rows)
v["BDSPPatientID"] = v.BDSPPatientID.astype(int)
print(f"seed {draw['seed']}   drew {draw['n']}   verify records {len(v)}   "
      f"status {v.status.value_counts().to_dict()}")
assert len(v) == draw["n"], (len(v), draw["n"])
assert set(v.BDSPPatientID) == set(draw["BDSPPatientID"])

f = pd.read_csv(os.path.join(HERE, "oxyprofile_per_patient.csv"), low_memory=False)
f["BDSPPatientID"] = f.BDSPPatientID.astype(int)
j = v.merge(f, on="BDSPPatientID", how="left", suffixes=("_v", "_f"))
assert len(j) == len(v)
print(f"fleet status on these 50: {j.status_f.value_counts().to_dict()}")

records = []
for _, r in j.iterrows():
    pid = int(r.BDSPPatientID)
    records.append({"BDSPPatientID": pid, "quantity": "whole recording mean SpO2",
                    "unit": "percent", "fleet": r.get("rec_mean"),
                    "independent": r.get("v_rec_mean"), "tol": TOL_MEAN_PP})
    for t in DEC:
        records.append({"BDSPPatientID": pid, "quantity": f"decile {int(t)} mean SpO2",
                        "unit": "percent", "fleet": r.get(f"dec_mean_{t}"),
                        "independent": r.get(f"v_dec_mean_{t}"), "tol": TOL_MEAN_PP})
        records.append({"BDSPPatientID": pid, "quantity": f"decile {int(t)} valid samples",
                        "unit": "samples", "fleet": r.get(f"dec_n_{t}"),
                        "independent": r.get(f"v_dec_n_{t}"), "tol": TOL_COUNT})
    for s in STAGES:
        records.append({"BDSPPatientID": pid, "quantity": f"{s} mean SpO2",
                        "unit": "percent", "fleet": r.get(f"{s}_mean"),
                        "independent": r.get(f"v_{s}_mean"), "tol": TOL_MEAN_PP})
        records.append({"BDSPPatientID": pid, "quantity": f"{s} minutes",
                        "unit": "minutes", "fleet": r.get(f"{s}_min"),
                        "independent": r.get(f"v_{s}_min"), "tol": TOL_MIN})

c = pd.DataFrame(records)
c["fleet"] = pd.to_numeric(c.fleet, errors="coerce")
c["independent"] = pd.to_numeric(c.independent, errors="coerce")
c["diff"] = c.independent - c.fleet
c["abs_diff"] = c["diff"].abs()
c["both_present"] = c.fleet.notna() & c.independent.notna()
c["only_one"] = c.fleet.notna() ^ c.independent.notna()
c["over_tol"] = c.both_present & (c.abs_diff > c.tol)
c.to_csv(os.path.join(HERE, "verify50.csv"), index=False)

print("\nMAX ABSOLUTE DIFFERENCE BY QUANTITY GROUP")
groups = {
    "whole recording mean SpO2 (percent)": c.quantity == "whole recording mean SpO2",
    "decile mean SpO2 (percent)": c.quantity.str.contains("decile") & (c.unit == "percent"),
    "decile valid samples (count)": c.unit == "samples",
    "stage mean SpO2 (percent)": c.quantity.str.contains("mean SpO2")
    & ~c.quantity.str.contains("decile|whole recording"),
    "stage minutes (minutes)": c.unit == "minutes",
}
summary_groups = {}
for name, mask in groups.items():
    g = c[mask & c.both_present]
    mx = float(g.abs_diff.max()) if len(g) else float("nan")
    over = int((g.abs_diff > g.tol).sum()) if len(g) else 0
    summary_groups[name] = {"n_compared": int(len(g)), "max_abs_diff": mx,
                            "median_abs_diff": float(g.abs_diff.median()) if len(g) else None,
                            "n_over_tolerance": over}
    print(f"  {name:<38} n {len(g):5,}  max |diff| {mx:.3e}  median {g.abs_diff.median():.3e}"
          f"  over tolerance {over}")

miss = c[c.only_one]
print(f"\nvalues present in one implementation and not the other: {len(miss)}")
if len(miss):
    print(miss[["BDSPPatientID", "quantity", "fleet", "independent"]]
          .to_string(index=False))

bad = c[c.over_tol]
print(f"\nRECORDS OVER TOLERANCE: {bad.BDSPPatientID.nunique()} patients, {len(bad)} values")
if len(bad):
    print(bad.sort_values("abs_diff", ascending=False)
          [["BDSPPatientID", "quantity", "unit", "fleet", "independent", "abs_diff"]]
          .head(40).to_string(index=False))

worst = c[c.both_present].sort_values("abs_diff", ascending=False).head(8)
print("\nlargest eight differences overall")
print(worst[["BDSPPatientID", "quantity", "unit", "fleet", "independent", "abs_diff"]]
      .to_string(index=False, float_format=lambda x: f"{x:.6f}"))

mean_mask = (c.unit == "percent") & c.both_present
minute_mask = (c.unit == "minutes") & c.both_present
count_mask = (c.unit == "samples") & c.both_present
gate_means = bool((c[mean_mask].abs_diff <= TOL_MEAN_PP).all())
gate_minutes = bool((c[minute_mask].abs_diff <= TOL_MIN).all())
gate_counts = bool((c[count_mask].abs_diff <= TOL_COUNT).all())
gate_presence = len(miss) == 0
ALL = gate_means and gate_minutes and gate_counts and gate_presence

print("\n" + "=" * 78)
print(f"GATE A (every mean within {TOL_MEAN_PP} pp):        {'PASS' if gate_means else 'FAIL'}"
      f"  [{int((c[mean_mask].abs_diff <= TOL_MEAN_PP).sum()):,}/{int(mean_mask.sum()):,}]")
print(f"GATE B (every stage minute within {TOL_MIN} min): {'PASS' if gate_minutes else 'FAIL'}"
      f"  [{int((c[minute_mask].abs_diff <= TOL_MIN).sum()):,}/{int(minute_mask.sum()):,}]")
print(f"GATE C (decile counts within {TOL_COUNT} samples):    "
      f"{'PASS' if gate_counts else 'FAIL'}"
      f"  [{int((c[count_mask].abs_diff <= TOL_COUNT).sum()):,}/{int(count_mask.sum()):,}]")
print(f"GATE D (no value present in only one run):    {'PASS' if gate_presence else 'FAIL'}")
print(f"VERDICT: {'INDEPENDENT VERIFICATION PASSES' if ALL else 'STOP AND DIAGNOSE'}")
print("=" * 78)

json.dump({
    "seed": draw["seed"], "n_drawn": draw["n"],
    "n_patient_ids": len(draw["BDSPPatientID"]),
    "sites": draw["sites"],
    "n_verify_records": int(len(v)),
    "n_values_compared": int(c.both_present.sum()),
    "tolerances": {"mean_pp": TOL_MEAN_PP, "minutes": TOL_MIN, "counts_samples": TOL_COUNT},
    "by_group": summary_groups,
    "n_values_over_tolerance": int(len(bad)),
    "patients_over_tolerance": sorted(int(x) for x in bad.BDSPPatientID.unique()),
    "n_values_present_in_only_one": int(len(miss)),
    "max_abs_diff_any_mean_pp": float(c[mean_mask].abs_diff.max()),
    "max_abs_diff_any_minutes": float(c[minute_mask].abs_diff.max()),
    "max_abs_diff_any_count": float(c[count_mask].abs_diff.max()),
    "gate_means": gate_means, "gate_minutes": gate_minutes,
    "gate_counts": gate_counts, "gate_presence": gate_presence,
    "verdict": "PASS" if ALL else "STOP",
}, open(os.path.join(HERE, "verify50_summary.json"), "w"), indent=2)
print("wrote verify50.csv, verify50_summary.json")
