"""
EXTRA. The one thing the positive control turned up on its own.

The published age adjustment is four cr(df=4) columns that sum to one. The basis is singular and
survives only because the ranking penalises at 0.01. s00 refitted T90 with one column dropped,
which makes the design full rank, and T90's dC rose from +0.020339 to +0.021903, a 7.7 percent
rise for a change that touches nothing about the exposure.

That is not by itself a stronger result. Making the age adjustment full rank changes the BASE
model, so it could lift every measure equally and leave the ranking untouched. This refits the
nine measures that sit above the noise floor at the top of ranking_v3.csv under both age bases
and asks whether T90's LEAD moves, which is the only thing that would matter.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import json
import os

import numpy as np
import pandas as pd

from common import ADJ, PEN, T90ROOT, WORK, build_frame, heldout, ranking_outcomes, site_rint
import sys as _sys_v8; _sys_v8.path.insert(0, paths.ANALYSIS_DIR)  # v8 sweep 2026-09-12
from cohort_spec import N_RANKED_OUTCOMES  # v8 sweep 2026-09-12

TOP = ["spo2_pct_below_90", "spo2_pct_below_88", "spo2_nadir_corrected", "odi3_total",
       "odi4_total", "nrem_hr_bpm", "wholenight_hr_bpm", "spo2_mean",
       "wholenight_lf_hf_ratio"]
ADJ_FULLRANK = [f"age_s{i}" for i in range(1, 4)] + ["male"]

print("=" * 78)
print("EXTRA. THE AGE BASIS: SINGULAR (PUBLISHED) AGAINST FULL RANK")
print("=" * 78)

d = build_frame(extra_master_cols=TOP, cache="frame_age.pkl")
OUT = ranking_outcomes(d)
assert len(OUT) == N_RANKED_OUTCOMES
specs = {"__base__": []}
for f in TOP:
    d[f + "__z"] = site_rint(d, f)
    specs[f] = [f + "__z"]

out = {}
for tag, adj in (("published_4col", ADJ), ("fullrank_3col", ADJ_FULLRANK)):
    res = heldout(specs, OUT, d, penalizer=PEN, adj=adj, pkl_name=f"age_{tag}.pkl")
    b = float(res[res.spec == "__base__"].c.mean())
    s = (res[res.spec != "__base__"].groupby("spec").agg(dC=("gain", "mean"))
         .reset_index().sort_values("dC", ascending=False))
    s["rank"] = range(1, len(s) + 1)
    out[tag] = {"baseline": b, "table": s}
    print(f"\n{tag}: held-out baseline {b:.6f}")

RK = pd.read_csv(f"{paths.NUMBERS_DIR}/ranking_v3.csv", comment="#")
pub = RK.set_index("feature").dC

a = out["published_4col"]["table"].set_index("spec")
c = out["fullrank_3col"]["table"].set_index("spec")
print(f"\n{'measure':<26}{'published':>12}{'full rank':>12}{'change':>11}{'%':>8}"
      f"{'rank pub':>10}{'rank full':>11}")
print("-" * 90)
for f in a.sort_values("dC", ascending=False).index:
    print(f"{f:<26}{a.loc[f, 'dC']:>+12.6f}{c.loc[f, 'dC']:>+12.6f}"
          f"{c.loc[f, 'dC'] - a.loc[f, 'dC']:>+11.6f}"
          f"{100 * (c.loc[f, 'dC'] / a.loc[f, 'dC'] - 1):>+8.1f}"
          f"{int(a.loc[f, 'rank']):>10}{int(c.loc[f, 'rank']):>11}")

T = "spo2_pct_below_90"
runner_pub = a.drop(index=T).dC.max()
runner_full = c.drop(index=T).dC.max()
lead_pub = a.loc[T, "dC"] - runner_pub
lead_full = c.loc[T, "dC"] - runner_full
print(f"\nbaseline moves {out['published_4col']['baseline']:.6f} -> "
      f"{out['fullrank_3col']['baseline']:.6f} "
      f"({out['fullrank_3col']['baseline'] - out['published_4col']['baseline']:+.6f})")
print(f"T90 dC moves {a.loc[T, 'dC']:+.6f} -> {c.loc[T, 'dC']:+.6f} "
      f"({100 * (c.loc[T, 'dC'] / a.loc[T, 'dC'] - 1):+.1f}%)")
print(f"T90's lead over the runner-up moves {lead_pub:+.6f} -> {lead_full:+.6f} "
      f"({lead_full - lead_pub:+.6f})")
print(f"mean rise across all {len(TOP)} top measures: "
      f"{float((c.dC - a.dC).mean()):+.6f}   "
      f"T90's own rise: {float(c.loc[T, 'dC'] - a.loc[T, 'dC']):+.6f}")
verdict = ("T90 gains no more than its competitors, the rise is the base model, not the exposure"
           if float(c.loc[T, "dC"] - a.loc[T, "dC"]) <= float((c.dC - a.dC).mean()) + 1e-9
           else "T90 gains more than the average competitor")
print(f"VERDICT: {verdict}")

df = a.join(c, lsuffix="_pub", rsuffix="_full").reset_index()
df.insert(0, "analysis", "EXTRA_agebasis")
df.to_csv(os.path.join(WORK, "EXTRA_agebasis.csv"), index=False)
json.dump({"baseline_published": out["published_4col"]["baseline"],
           "baseline_fullrank": out["fullrank_3col"]["baseline"],
           "t90_dC_published": float(a.loc[T, "dC"]),
           "t90_dC_fullrank": float(c.loc[T, "dC"]),
           "lead_published": float(lead_pub), "lead_fullrank": float(lead_full),
           "mean_rise_top_measures": float((c.dC - a.dC).mean()),
           "t90_rise": float(c.loc[T, "dC"] - a.loc[T, "dC"]),
           "verdict": verdict},
          open(os.path.join(WORK, "EXTRA_agebasis.json"), "w"), indent=2)
print("\nwritten -> _work/EXTRA_agebasis.csv, _work/EXTRA_agebasis.json")
