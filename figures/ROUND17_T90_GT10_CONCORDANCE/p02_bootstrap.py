"""
The bootstrap for the binary exposure. Identical design to the published one
(Sleep_Variability_2026-08/dC_per_disease/p01_bootstrap.py): patients resampled with
replacement WITHIN hospital, the whole procedure repeated on each resample, 500 replicates,
percentile interval.

SEED0 is the published SEED0. Replicate r uses seed 20260820 + r in both this run and the
published continuous run, so replicate r of the binary bootstrap and replicate r of
_work/boot_dc_standard_drop.npy are the same resampled patients. That makes the binary-minus-
continuous difference a paired quantity with its own interval.

Arms
  lag0   full follow-up, reported basis
  lag2   two-year landmark, reported basis, the published landmark rule

Usage: python3 p02_bootstrap.py [n_replicates] [n_jobs] [chunk] [arm]
"""
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"

import numpy as np                                                        # noqa: E402
from joblib import Parallel, delayed                                      # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import engine_gt10 as X                                                   # noqa: E402
import dc_engine as E                                                     # noqa: E402

SEED0 = 20260820          # the published seed base, so the resamples are the paper's own
ARMS = {"lag0": {"lag": 0.0, "ck": "_ck_gt10_lag0"},
        "lag2": {"lag": 2.0, "ck": "_ck_gt10_lag2"}}
LOG = os.path.join(X.WORK, "bootstrap.log")


def log(msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    with open(LOG, "a") as fh:
        fh.write(line + "\n")


def run_chunk(start, n, path, arm):
    if os.path.exists(path):
        return -1.0
    lag = ARMS[arm]["lag"]
    A = E.build_arrays()
    t0 = time.time()
    out = np.full((n, len(A["outcomes"]), 4), np.nan)
    for i in range(n):
        out[i] = X.one_pass_x(A, seed=SEED0 + start + i, basis="drop", lag=lag,
                              exposure="gt10")
    np.savez_compressed(path, seeds=np.arange(SEED0 + start, SEED0 + start + n), drop=out)
    return time.time() - t0


if __name__ == "__main__":
    N_REP = int(sys.argv[1]) if len(sys.argv) > 1 else 500
    N_JOBS = int(sys.argv[2]) if len(sys.argv) > 2 else 7
    CHUNK = int(sys.argv[3]) if len(sys.argv) > 3 else 5
    ARM = sys.argv[4] if len(sys.argv) > 4 else "lag0"
    assert ARM in ARMS
    CK = os.path.join(X.WORK, ARMS[ARM]["ck"])
    os.makedirs(CK, exist_ok=True)
    LOG = os.path.join(X.WORK, f"bootstrap_{ARM}.log")

    A = E.build_arrays()
    log(f"gt10 bootstrap, arm {ARM} (landmark {ARMS[ARM]['lag']:.0f} y): {N_REP} replicates, "
        f"{len(A['outcomes'])} outcomes, n_jobs={N_JOBS}, chunk={CHUNK}, "
        f"seeds {SEED0}..{SEED0 + N_REP - 1}")
    log(f"site sizes: I0002 {(A['site'] == 0).sum():,}  I0006 {(A['site'] == 1).sum():,}")

    jobs = [(s, min(CHUNK, N_REP - s), os.path.join(CK, f"boot_{s:05d}.npz"))
            for s in range(0, N_REP, CHUNK)]
    todo = [j for j in jobs if not os.path.exists(j[2])]
    log(f"{len(jobs)} chunks, {len(jobs) - len(todo)} on disk, {len(todo)} to run")

    t0, done = time.time(), 0
    for i in range(0, len(todo), N_JOBS * 2):
        batch = todo[i:i + N_JOBS * 2]
        Parallel(n_jobs=N_JOBS, backend="loky", verbose=0)(
            delayed(run_chunk)(s, n, p, ARM) for s, n, p in batch)
        done += len(batch)
        el = time.time() - t0
        log(f"  {done}/{len(todo)} chunks  elapsed {el / 60:.1f} min  "
            f"eta {(len(todo) - done) * (el / max(done, 1)) / 60:.1f} min")
    log(f"gt10 bootstrap done in {(time.time() - t0) / 60:.1f} min, arm {ARM}")
