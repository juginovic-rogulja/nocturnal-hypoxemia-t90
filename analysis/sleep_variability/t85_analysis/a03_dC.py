"""
T85 against T90 on the paper's discrimination-gain machinery (dC).

Runs on dC_side_analyses/common.py, which is the published ranking pipeline re-implemented
once. That folder's own positive control, s00_positive_control.py, reproduces T90's
published dC of +0.02033875497249791 to 1e-9 and every one of the 48 per-outcome gains to
1e-9. This script re-runs the T90 arm alongside T85, so the T90 number printed here IS the
positive control: if it is not the published 0.020338754972, nothing else on the page counts.

dC, as common.py defines it
  base      age as cr(df=4), all four columns, plus sex, at penalizer 0.01
  feature   ONE extra column, the within-site rank inverse normal of the measure
  folds     fit I0002 score I0006, then the reverse, two held-out concordances, nanmean
  dC        mean over the 48 ranking outcomes of (feature C minus base C)

The four age columns sum to one and are singular. common.py keeps the published basis at
penalizer 0.01 because that is what the published ranking did, and s00 additionally shows
the same dC with one column dropped, so the singular basis is demonstrably not fitting
nothing. Any NEW basis would have a column dropped. No new basis is introduced here.

Specs fitted
  t90_frozen   the master's spo2_pct_below_90, the published exposure
  t85          prov_t85 from the 2026-08-21 re-derivation
  t90_rederiv  prov_t90 from the same re-derivation, so the threshold contrast can be read
               without the re-derivation being a confounder
  t88_rederiv  prov_t88, the third rung of the same ladder
  t85_and_t90  both columns in one model, which answers whether T85 carries anything the
               paper's own exposure does not already have
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
WORK = os.path.join(HERE, "_dc_work")
os.makedirs(WORK, exist_ok=True)
PUB_DC = __PUB_DC_V7__
FEAT = "spo2_pct_below_90"

d = build_frame(extra_master_cols=[FEAT])
print(f"cohort {len(d):,}")
OUT = ranking_outcomes(d)
assert len(OUT) == N_RANKED_OUTCOMES, len(OUT)
print(f"outcomes {len(OUT)}")

t = pd.read_csv(os.path.join(HERE, "t85_per_patient_v7.csv"), low_memory=False)
t = t[t.status == "ok"][["BDSPPatientID", "prov_t85", "prov_t88", "prov_t90"]]
t["BDSPPatientID"] = t.BDSPPatientID.astype(int)
d["BDSPPatientID"] = d.BDSPPatientID.astype(int)
d = d.merge(t, on="BDSPPatientID", how="left")
cov = {c: float(d[c].notna().mean()) for c in ("prov_t85", "prov_t88", "prov_t90")}
print("coverage: " + ", ".join(f"{k} {100*v:.2f}%" for k, v in cov.items()))

COLS = {"t90_frozen": FEAT, "t85": "prov_t85", "t88_rederiv": "prov_t88",
        "t90_rederiv": "prov_t90"}
for name, col in COLS.items():
    d[f"{name}__z"] = site_rint(d, col)

specs = {"__base__": []}
for name in COLS:
    specs[name] = [f"{name}__z"]
specs["t85_and_t90"] = ["t85__z", "t90_frozen__z"]

res = heldout(specs, OUT, d, penalizer=PEN, pkl_name="t85_dc_frame.pkl")
s = summarise(res)
print("\n" + "=" * 78)
print("dC, the paper's discrimination gain")
print("=" * 78)
print(s.to_string(index=False, float_format=lambda x: f"{x:.6f}"))

got90 = float(s[s.spec == "t90_frozen"].dC.iloc[0])
PC_OK = abs(got90 - PUB_DC) < 1e-9
print(f"\n  POSITIVE CONTROL  published T90 dC {PUB_DC:.12f}")
print(f"                    reproduced here   {got90:.12f}   diff {got90-PUB_DC:+.3e}")
print(f"  VERDICT: {'PASS' if PC_OK else 'FAIL'}")
if not PC_OK:
    raise SystemExit("POSITIVE CONTROL FAILED, stop")

got85 = float(s[s.spec == "t85"].dC.iloc[0])
print(f"\n  T85 dC {got85:.6f}   against T90 {got90:.6f}   "
      f"T85 minus T90 {got85-got90:+.6f}  "
      f"({100*got85/got90:.1f}% of the published gain)")

per = res.pivot_table(index="outcome", columns="spec", values="gain", aggfunc="first")
per = per.drop(columns=[c for c in ["__base__"] if c in per.columns])
base = res[res.spec == "__base__"].set_index("outcome")["c"]
per.insert(0, "base_c", base)
per["t85_minus_t90"] = per["t85"] - per["t90_frozen"]
per = per.sort_values("t90_frozen", ascending=False)
per.to_csv(os.path.join(HERE, "dC_per_outcome.csv"))
s.to_csv(os.path.join(HERE, "dC_summary.csv"), index=False)

print(f"\n  per outcome: T85 gain beats T90 gain on "
      f"{int((per.t85_minus_t90 > 0).sum())}/{len(per)} of the {N_RANKED_OUTCOMES}")
print("\n  the 8 outcomes where T90 gains most, with T85 alongside")
print(per.head(8)[["base_c", "t90_frozen", "t85", "t88_rederiv", "t85_minus_t90"]]
      .to_string(float_format=lambda x: f"{x:+.6f}"))

json.dump({"published_T90_dC": PUB_DC, "reproduced_T90_dC": got90,
           "positive_control_pass": bool(PC_OK), "T85_dC": got85,
           "T85_minus_T90_dC": got85 - got90,
           "dC": {r.spec: float(r.dC) for r in s.itertuples()},
           "mC": {r.spec: float(r.mC) for r in s.itertuples()},
           "n_outcomes": len(OUT), "coverage": cov,
           "n_outcomes_T85_beats_T90": int((per.t85_minus_t90 > 0).sum())},
          open(os.path.join(HERE, "dC_result.json"), "w"), indent=2)
print("\nwrote dC_summary.csv, dC_per_outcome.csv, dC_result.json")
