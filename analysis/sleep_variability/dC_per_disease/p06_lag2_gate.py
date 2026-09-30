"""
POSITIVE CONTROL for the 2-year landmark arm. The lag-2 bootstrap does not run until this
passes.

There is no published dC at lag 2 to reproduce, so the gate is on the windowing instead, which
is the only thing the lag arm changes. numbers/lag_ladder.csv is the paper's landmark ladder,
built by New_Figures/scripts/lagladder.py on the same frozen cohort with the same eligibility.
Its n and events at lag_years == 2 are therefore exactly what this folder's landmark must
produce, outcome by outcome, with no tolerance at all.

  GATE A  lag 0. Every outcome's n and events match the ladder's lag-0 row exactly. This
          proves the eligibility rule is the same one before the landmark is applied.
  GATE B  lag 2. Every outcome's n and events match the ladder's lag-2 row exactly.

Then the estimability check the ranking's own floor implies: at lag 2 an outcome needs 30
events on each side of the cross-fit. Outcomes that fall under it are listed, because their
dC will be NaN rather than small, and that has to be visible rather than silent.
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

import dc_engine as E

T90ROOT = paths.T90_ROOT
LAG = 2.0
MIN = E.MIN_FOLD_EVENTS

print("=" * 78)
print("POSITIVE CONTROL, 2-year landmark arm")
print("=" * 78)

A = E.build_arrays()
from cohort_spec import cohort_tag   # v8.1c
lad = pd.read_csv(f"{paths.NUMBERS_DIR}/lag_ladder{cohort_tag()}.csv")
print("ladder convention, numbers/lag_ladder.json: "
      + json.load(open(f"{paths.NUMBERS_DIR}/lag_ladder{cohort_tag()}.json"))["provenance"]["landmark"])

ok = True
for gate, L in (("A", 0.0), ("B", LAG)):
    mine = E.counts(A, lag=L).set_index("outcome")
    ref = lad[lad.lag_years == L].set_index("key")[["n", "events"]]
    missing = [o for o in A["outcomes"] if o not in ref.index]
    assert not missing, f"outcomes absent from the ladder: {missing}"
    bad = []
    for o in A["outcomes"]:
        if int(mine.loc[o, "n_at_risk"]) != int(ref.loc[o, "n"]):
            bad.append(f"{o} n {int(mine.loc[o, 'n_at_risk'])} vs ladder {int(ref.loc[o, 'n'])}")
        if int(mine.loc[o, "events"]) != int(ref.loc[o, "events"]):
            bad.append(f"{o} events {int(mine.loc[o, 'events'])} vs "
                       f"ladder {int(ref.loc[o, 'events'])}")
    print(f"\nGATE {gate}  lag {L:.0f}: {len(A['outcomes'])} outcomes, "
          f"{'every n and events match the published ladder exactly' if not bad else bad[:6]}")
    ok = ok and not bad

c0 = E.counts(A, lag=0.0).set_index("outcome")
c2 = E.counts(A, lag=LAG).set_index("outcome")
lost = 100 * (1 - c2.events / c0.events)
print(f"\nlag 2 removes {lost.min():.0f}% to {lost.max():.0f}% of events, "
      f"median {lost.median():.0f}%")

thin = [(o, int(c2.loc[o, "events_I0002"]), int(c2.loc[o, "events_I0006"]))
        for o in A["outcomes"]
        if min(c2.loc[o, "events_I0002"], c2.loc[o, "events_I0006"]) < MIN]
print(f"outcomes below the {MIN}-event fold floor at lag 2 on the observed data: {len(thin)}")
for o, a, b in thin:
    print(f"    {o:<20} I0002 {a}, I0006 {b}")

near = [(o, int(min(c2.loc[o, 'events_I0002'], c2.loc[o, 'events_I0006'])))
        for o in A["outcomes"]
        if MIN <= min(c2.loc[o, "events_I0002"], c2.loc[o, "events_I0006"]) < 2 * MIN]
print(f"outcomes within a factor of two of the floor, so at risk of losing a fold in some "
      f"bootstrap replicates: {len(near)}")
for o, m in near:
    print(f"    {o:<20} smaller side {m}")

print("\n" + "=" * 78)
print("VERDICT: " + ("PASS" if ok else "FAIL"))
print("=" * 78)

json.dump({
    "landmark_years": LAG,
    "convention": "event or censoring at or before L removes the patient, clock restarts at L",
    "gate_counts_match_published_ladder": bool(ok),
    "pct_events_removed_at_lag2": {"min": float(lost.min()), "median": float(lost.median()),
                                   "max": float(lost.max())},
    "below_fold_floor_at_lag2": [t[0] for t in thin],
    "near_fold_floor_at_lag2": [t[0] for t in near],
}, open(os.path.join(E.WORK, "lag2_gate.json"), "w"), indent=2)

if not ok:
    raise SystemExit("LAG-2 GATE FAILED, stop")
