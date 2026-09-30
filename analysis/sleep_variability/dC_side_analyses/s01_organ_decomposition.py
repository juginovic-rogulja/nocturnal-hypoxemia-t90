"""
A. Where T90's discrimination gain actually comes from.

The published +0.020339 is the plain mean of 48 per-outcome gains. This splits that mean by
the organ system each condition belongs to, using ORGAN_GROUP from numbers/disease_definitions.py
verbatim, so the grouping cannot be tuned to the answer. Death is its own group. The
cardiovascular composite is reported on its own line as well, because ORGAN_GROUP files it as
"Composite" rather than inside Cardiac and it would otherwise be invisible.

Input is the per-outcome frame s00 produced, which reproduces ranking_v3_percondition.csv
exactly, so this is a regrouping of the published numbers and not a refit.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import os

import numpy as np
import pandas as pd

from common import CVD_COMPONENTS, DISEASES, WORK, organ_of
import sys as _sys_v8; _sys_v8.path.insert(0, paths.ANALYSIS_DIR)  # v8 sweep 2026-09-12
from cohort_spec import N_RANKED_OUTCOMES  # v8 sweep 2026-09-12

res = pd.read_csv(os.path.join(WORK, "pc_percondition.csv"))
t = res[res.spec == "t90"][["outcome", "c", "gain"]].copy()
assert len(t) == N_RANKED_OUTCOMES, len(t)

t["organ"] = t.outcome.map(organ_of)
t["label"] = t.outcome.map(lambda k: DISEASES[k][0] if k in DISEASES else "Death from any cause")
assert (t.organ != "UNGROUPED").all(), t[t.organ == "UNGROUPED"]

overall = float(t.gain.mean())
print("=" * 78)
print("A. PER-ORGAN DECOMPOSITION OF T90's dC")
print("=" * 78)
print(f"{N_RANKED_OUTCOMES}-outcome mean, the published number: {overall:+.6f}\n")

g = (t.groupby("organ")
       .agg(n=("gain", "size"), dC=("gain", "mean"), median=("gain", "median"),
            worst=("gain", "min"), best=("gain", "max"),
            npos=("gain", lambda x: int((x > 0).sum())))
       .reset_index().sort_values("dC", ascending=False))
g["vs_overall"] = g.dC - overall

print(f"{'organ group':<16}{'n':>3}{'dC':>12}{'median':>11}{'worst':>11}{'best':>11}"
      f"{'n>0':>5}{'vs ' + str(N_RANKED_OUTCOMES) + '-mean':>12}")
print("-" * 82)
for _i, r in g.iterrows():
    print(f"{r.organ:<16}{int(r.n):>3}{r.dC:>+12.6f}{r['median']:>+11.6f}{r.worst:>+11.6f}"
          f"{r.best:>+11.6f}{int(r.npos):>5}{r.vs_overall:>+12.6f}")

top, bot = g.iloc[0], g.iloc[-1]
print(f"\nspread the {N_RANKED_OUTCOMES}-outcome average hides: {top.organ} {top.dC:+.6f} down to "
      f"{bot.organ} {bot.dC:+.6f}, a range of {top.dC - bot.dC:.6f}")
print(f"that is {top.dC / overall:.1f}-fold the published average at the top and "
      f"{bot.dC / overall:+.2f}-fold at the bottom")

# ------------------------------------------------------------ the cardiovascular composite
cvd = t[t.outcome == "cvd"]
print(f"\ncardiovascular composite (union of {', '.join(CVD_COMPONENTS)}): "
      f"dC {float(cvd.gain.iloc[0]):+.6f}, held-out C {float(cvd.c.iloc[0]):.6f}")
comp = t[t.outcome.isin(CVD_COMPONENTS)]
print(f"its five components separately: mean dC {float(comp.gain.mean()):+.6f}, "
      f"range {float(comp.gain.min()):+.6f} to {float(comp.gain.max()):+.6f}")
for _i, r in comp.sort_values("gain", ascending=False).iterrows():
    print(f"   {r.label:<26}{r.gain:+.6f}")
print("the composite scores below the mean of its own parts: pooling heart failure's "
      f"{float(comp[comp.outcome == 'hf'].gain.iloc[0]):+.6f} with stroke's "
      f"{float(comp[comp.outcome == 'stroke_any'].gain.iloc[0]):+.6f} dilutes it.")

print("\nevery outcome, ordered")
for _i, r in t.sort_values("gain", ascending=False).iterrows():
    print(f"   {r.organ:<14}{r.label:<32}{r.gain:+.6f}")

g.insert(0, "analysis", "A_organ")
t.insert(0, "analysis", "A_organ_percondition")
g.to_csv(os.path.join(WORK, "A_organ_summary.csv"), index=False)
t.to_csv(os.path.join(WORK, "A_organ_percondition.csv"), index=False)
print("\nwritten -> _work/A_organ_summary.csv, _work/A_organ_percondition.csv")
