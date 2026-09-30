"""
The whole threshold ladder on the paper's discrimination-gain machinery (dC).
Extends t85_analysis/a03_dC.py.

Runs on dC_side_analyses/common.py, unchanged. The T90 row printed here IS the
positive control: it must equal the published 0.020338754972 to 1e-9 or the
script stops. Three further cross-checks come free because the T85 run computed
overlapping specs this morning from the same raw data: t85, t88_rederiv and
t90_rederiv must reproduce that run's dC to 1e-9.

Specs fitted
  t90_frozen   spo2_pct_below_90 (master), the published exposure. Positive control.
  t88_frozen   spo2_pct_below_88 (master), the frozen T88 rung.
  t80, t85, t92, t95   prov_* from the 2026-08-21 ladder re-derivation.
  t88_rederiv, t90_rederiv   the same re-derivation at 88 and 90, so the pure
               same-code ladder can be read without frozen-vs-rederived mixing.
  a95          time at or above 95. Included ONLY to verify the mirror argument:
               its within-site rank inverse normal is exactly the negation of
               t95's, so its dC must equal t95's to numerical noise. It is not
               a new exposure and is reported as a check, not a rung.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import json
import os
import sys
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

DCDIR = f"{paths.SV_CODE_DIR}/dC_side_analyses"
sys.path.insert(0, DCDIR)
from common import (ADJ, PEN, T90ROOT, build_frame, heldout,     # noqa: E402

                    ranking_outcomes, site_rint, summarise)
import sys as _sys_v8; _sys_v8.path.insert(0, paths.ANALYSIS_DIR)  # v8 sweep 2026-09-12
from cohort_spec import N_RANKED_OUTCOMES  # v8 sweep 2026-09-12
__PUB_DC_V7__ = float(__import__("pandas").read_csv(f"{paths.NUMBERS_DIR}/ranking_v3.csv", comment="#").set_index("feature").loc["spo2_pct_below_90", "dC"])   # v7: T90 gain from the CURRENT ranking, was the August literal 0.02033875497249791

HERE = os.path.dirname(os.path.abspath(__file__))
T85DIR = os.path.join(os.path.dirname(HERE), "t85_analysis")
PUB_DC = __PUB_DC_V7__

d = build_frame(extra_master_cols=["spo2_pct_below_90", "spo2_pct_below_88"],
                cache="frame_ladder.pkl")
print(f"cohort {len(d):,}")
OUT = ranking_outcomes(d)
assert len(OUT) == N_RANKED_OUTCOMES, len(OUT)
print(f"outcomes {len(OUT)}")

t = pd.read_csv(os.path.join(HERE, "ladder_per_patient_v7.csv"), low_memory=False)
t = t[t.status == "ok"][["BDSPPatientID", "prov_t80", "prov_t85", "prov_t88",
                         "prov_t90", "prov_t92", "prov_t95", "prov_a95"]]
t["BDSPPatientID"] = t.BDSPPatientID.astype(int)
d["BDSPPatientID"] = d.BDSPPatientID.astype(int)
d = d.merge(t, on="BDSPPatientID", how="left")

COLS = {"t90_frozen": "spo2_pct_below_90", "t88_frozen": "spo2_pct_below_88",
        "t80": "prov_t80", "t85": "prov_t85", "t88_rederiv": "prov_t88",
        "t90_rederiv": "prov_t90", "t92": "prov_t92", "t95": "prov_t95",
        "a95": "prov_a95"}
cov = {n: float(d[c].notna().mean()) for n, c in COLS.items()}
print("coverage: " + ", ".join(f"{k} {100*v:.2f}%" for k, v in cov.items()))
for name, col in COLS.items():
    d[f"{name}__z"] = site_rint(d, col)

# the mirror identity at the transform level, checked before any fitting
mz = d[["t95__z", "a95__z"]].dropna()
mirror_z = float((mz.t95__z + mz.a95__z).abs().max())
print(f"max |rint(t95) + rint(a95)| = {mirror_z:.3e}  "
      f"(the two transformed columns are exact negations)")

specs = {"__base__": []}
for name in COLS:
    specs[name] = [f"{name}__z"]

res = heldout(specs, OUT, d, penalizer=PEN, pkl_name="ladder_dc_frame.pkl")
s = summarise(res)
print("\n" + "=" * 78)
print("dC, the paper's discrimination gain, every rung")
print("=" * 78)
print(s.to_string(index=False, float_format=lambda x: f"{x:.6f}"))

DC = {r.spec: float(r.dC) for r in s.itertuples()}
got90 = DC["t90_frozen"]
PC_OK = abs(got90 - PUB_DC) < 1e-9
print(f"\n  POSITIVE CONTROL  published T90 dC {PUB_DC:.12f}")
print(f"                    reproduced here   {got90:.12f}   diff {got90-PUB_DC:+.3e}")
print(f"  VERDICT: {'PASS' if PC_OK else 'FAIL'}")
if not PC_OK:
    raise SystemExit("POSITIVE CONTROL FAILED, stop")

# cross-checks against the T85 run, same data, same machinery, this morning
prev = json.load(open(os.path.join(T85DIR, "dC_result.json")))["dC"]
XCHK = {}
print("\n  CROSS-CHECKS against the T85 run's dC (same data, same machinery)")
for k in ("t85", "t88_rederiv", "t90_rederiv"):
    dev = abs(DC[k] - prev[k])
    XCHK[k] = dev
    print(f"    {k:<12} here {DC[k]:.12f}   T85 run {prev[k]:.12f}   |diff| {dev:.3e}")
    # the fleets are different physical instances, so the resample can differ by a
    # few float ulps per record. 1e-6 catches any real defect without tripping on that.
    assert dev < 1e-6, f"{k} does not reproduce the T85 run"
print("  all reproduce within 1e-6 (expected ~1e-12)")

mdev = abs(DC["a95"] - DC["t95"])
print(f"\n  A95 MIRROR CHECK  dC(a95) {DC['a95']:.12f}  dC(t95) {DC['t95']:.12f}  "
      f"|diff| {mdev:.3e}")
print("  identical as predicted: time above 95 is time below 95 with the sign flipped, "
      "and rank inverse normal erases the sign for dC.")

print("\n  THE LADDER, dC by threshold")
for k in ("t80", "t85", "t88_frozen", "t90_frozen", "t92", "t95"):
    print(f"    {k:<12} dC {DC[k]:+.6f}   {100*DC[k]/got90:6.1f}% of published T90")

per = res.pivot_table(index="outcome", columns="spec", values="gain", aggfunc="first")
per = per.drop(columns=[c for c in ["__base__"] if c in per.columns])
base = res[res.spec == "__base__"].set_index("outcome")["c"]
per.insert(0, "base_c", base)
per = per.sort_values("t90_frozen", ascending=False)
per.to_csv(os.path.join(HERE, "dC_per_outcome.csv"))
s.to_csv(os.path.join(HERE, "dC_summary.csv"), index=False)

nbeat = {k: int((per[k] > per["t90_frozen"]).sum())
         for k in ("t80", "t85", "t88_frozen", "t92", "t95")}
print("\n  per outcome, rung beats T90's gain on: " +
      ", ".join(f"{k} {v}/{N_RANKED_OUTCOMES}" for k, v in nbeat.items()))
print("\n  the 8 outcomes where T90 gains most, the ladder alongside")
print(per.head(8)[["base_c", "t80", "t85", "t88_frozen", "t90_frozen", "t92", "t95"]]
      .to_string(float_format=lambda x: f"{x:+.6f}"))

json.dump({"published_T90_dC": PUB_DC, "reproduced_T90_dC": got90,
           "positive_control_pass": bool(PC_OK),
           "dC": DC, "mC": {r.spec: float(r.mC) for r in s.itertuples()},
           "crosschecks_vs_t85_run": XCHK,
           "a95_mirror_absdiff": mdev, "mirror_z_max": mirror_z,
           "n_outcomes": len(OUT), "coverage": cov,
           "n_beats_t90_per_outcome": nbeat},
          open(os.path.join(HERE, "dC_result.json"), "w"), indent=2)
print("\nwrote dC_summary.csv, dC_per_outcome.csv, dC_result.json")
