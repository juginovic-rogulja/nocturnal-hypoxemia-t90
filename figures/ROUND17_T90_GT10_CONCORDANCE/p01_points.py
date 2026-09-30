"""
Point estimates for the binary exposure, T90 above 10% of the recording, entered raw as 0/1.

Three passes, all through engine_gt10.one_pass_x, which the self-test has proved is the
published engine with only the exposure column changed:

  gt10_drop_lag0   full follow-up, reported basis (one age-spline column dropped)
  gt10_drop_lag2   two-year landmark, reported basis
  gt10_full_lag0   full follow-up, published four-column basis, cross-reference only

Also writes the exposure counts per outcome on the eligible set the folds actually see.
"""
import os
import sys
import time
import warnings

warnings.filterwarnings("ignore")

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import engine_gt10 as X                                                   # noqa: E402
import dc_engine as E                                                     # noqa: E402

A = E.build_arrays()
t0 = time.time()
for tag, basis, lag in (("gt10_drop_lag0", "drop", 0.0),
                        ("gt10_drop_lag2", "drop", 2.0),
                        ("gt10_full_lag0", "full", 0.0)):
    f = f"{X.WORK}/point_{tag}.npy"
    if os.path.exists(f):
        print(f"{tag}: on disk")
        continue
    P = X.one_pass_x(A, seed=None, basis=basis, lag=lag, exposure="gt10")
    np.save(f, P)
    print(f"{tag}: mean dC {np.nanmean(E.dc_from_pass(P)):+.8f}   "
          f"{time.time() - t0:.0f}s", flush=True)

for lag, tag in ((0.0, "lag0"), (2.0, "lag2")):
    X.exposed_counts(A, lag=lag).to_csv(f"{X.WORK}/counts_{tag}.csv")
print("counts written")
