"""
The observed-data point estimates for both arms, under the names the assembly expects.

standard  reuses p00_positive_control.py's gated passes, which are the ones the positive
          control checked against the published ranking. They are copied, not recomputed, so
          the delivered numbers are literally the gated ones.
lag2      one pass on the 2-year landmark, reported basis only.
"""
import os
import warnings

warnings.filterwarnings("ignore")

import numpy as np

import dc_engine as E

A = E.build_arrays()

for b in ("drop", "full"):
    src = os.path.join(E.WORK, f"point_{b}.npy")
    dst = os.path.join(E.WORK, f"point_standard_{b}.npy")
    assert os.path.exists(src), f"{src} missing, run p00_positive_control.py first"
    np.save(dst, np.load(src))
    print(f"standard {b:<5} copied from the positive control, "
          f"mean dC {np.nanmean(E.dc_from_pass(np.load(dst))):.12f}")

dst = os.path.join(E.WORK, "point_lag2_drop.npy")
if not os.path.exists(dst):
    P = E.one_pass(A, seed=None, basis="drop", lag=2.0)
    np.save(dst, P)
else:
    P = np.load(dst)
dc = E.dc_from_pass(P)
bad = [o for o, v in zip(A["outcomes"], dc) if not np.isfinite(v)]
print(f"lag2     drop  computed, mean dC over the {int(np.isfinite(dc).sum())} estimable "
      f"{np.nanmean(dc):.12f}")
print(f"lag2     not estimable at the 30-event fold floor: {bad if bad else 'none'}")
