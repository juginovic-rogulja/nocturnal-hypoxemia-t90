"""
The bootstrap. Runs only after p00_positive_control.py has passed.

Design
  Resampling unit is the patient. The frozen cohort is one row per patient (19,173 rows,
  19,173 distinct BDSPPatientID), so a row bootstrap is a patient bootstrap.
  Resampling is stratified by hospital: 10,431 draws with replacement from I0002 and 8,742
  from I0006, so every replicate keeps the site sizes and the cross-fit design survives.
  The whole procedure is repeated on each resample, not just the last fit. The age spline
  knots are re-placed on the resampled ages and T90's within-site rank inverse normal is
  re-ranked inside the resample, because both are estimated on the data in the published
  pipeline.

  Both age bases are carried in the same replicate, sharing one resample:
    drop  age_s1..age_s3 plus sex, the rank-safe basis this analysis reports
    full  the published four-column basis, for a cross-reference column whose point estimate
          is the paper's own

Output
  _work/_ck/boot_XXXXX.npz per chunk, so the run is resumable and can be assembled at any
  point. Each holds (chunk, 48, 4) for each basis: C base dir1, C T90 dir1, C base dir2,
  C T90 dir2, with dir1 fitting at I0002 and scoring at I0006.

Arms
  standard  the paper's windowing, both age bases carried in the same replicate
  lag2      a 2-year landmark, the published ladder's rule exactly (see p06_lag2_gate.py),
            reported basis only, because there is no published four-column number at lag 2
            to cross-reference and the second basis would double an already long run

Usage
  python3 p01_bootstrap.py [n_replicates] [n_jobs] [chunk] [arm]
"""
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"

import numpy as np                                                        # noqa: E402
from joblib import Parallel, delayed                                      # noqa: E402

import dc_engine as E                                                     # noqa: E402

SEED0 = 20260820          # replicate r uses seed SEED0 + r, the same draws in both arms

# Arm settings are module-level data, never read from argv at import time. A loky worker
# re-imports this module with its own argv, so parsing argv here would either pick up the
# worker's flags or crash it. The arm is passed into run_chunk instead.
ARMS = {"standard": {"lag": 0.0, "bases": ("drop", "full"), "ck": "_ck"},
        "lag2": {"lag": 2.0, "bases": ("drop",), "ck": "_ck_lag2"}}


def log(msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    with open(LOG, "a") as fh:
        fh.write(line + "\n")


LOG = os.path.join(E.WORK, "bootstrap.log")     # replaced in __main__ by the arm's own log


def run_chunk(start, n, path, arm="standard"):
    """One chunk of replicates. Writes its own checkpoint and returns its wall time."""
    if os.path.exists(path):
        return -1.0
    lag, bases = ARMS[arm]["lag"], ARMS[arm]["bases"]
    A = E.build_arrays()
    t0 = time.time()
    out = {b: np.full((n, len(A["outcomes"]), 4), np.nan) for b in bases}
    for i in range(n):
        seed = SEED0 + start + i
        for b in bases:
            out[b][i] = E.one_pass(A, seed=seed, basis=b, lag=lag)
    np.savez_compressed(path, seeds=np.arange(SEED0 + start, SEED0 + start + n), **out)
    return time.time() - t0


if __name__ == "__main__":
    N_REP = int(sys.argv[1]) if len(sys.argv) > 1 else 500
    N_JOBS = int(sys.argv[2]) if len(sys.argv) > 2 else 7
    CHUNK = int(sys.argv[3]) if len(sys.argv) > 3 else 5
    ARM = sys.argv[4] if len(sys.argv) > 4 else "standard"
    assert ARM in ARMS, f"arm must be one of {sorted(ARMS)}"
    LAG, BASES = ARMS[ARM]["lag"], ARMS[ARM]["bases"]
    CK = os.path.join(E.WORK, ARMS[ARM]["ck"])
    os.makedirs(CK, exist_ok=True)
    LOG = os.path.join(E.WORK, f"bootstrap_{ARM}.log")

    A = E.build_arrays()
    log(f"bootstrap start, arm {ARM} (landmark {LAG:.0f} y, bases {'+'.join(BASES)}): "
        f"{N_REP} replicates, {len(A['outcomes'])} outcomes, n_jobs={N_JOBS}, chunk={CHUNK}, "
        f"seeds {SEED0}..{SEED0 + N_REP - 1}")
    log(f"site sizes: I0002 {(A['site'] == 0).sum():,}  I0006 {(A['site'] == 1).sum():,}")

    jobs = []
    for s in range(0, N_REP, CHUNK):
        n = min(CHUNK, N_REP - s)
        jobs.append((s, n, os.path.join(CK, f"boot_{s:05d}.npz")))
    todo = [j for j in jobs if not os.path.exists(j[2])]
    log(f"{len(jobs)} chunks, {len(jobs) - len(todo)} already on disk, {len(todo)} to run")

    t0 = time.time()
    done = 0
    for i in range(0, len(todo), N_JOBS * 2):
        batch = todo[i:i + N_JOBS * 2]
        Parallel(n_jobs=N_JOBS, backend="loky", verbose=0)(
            delayed(run_chunk)(s, n, p, ARM) for s, n, p in batch)
        done += len(batch)
        el = time.time() - t0
        rate = el / max(done, 1)
        log(f"  {done}/{len(todo)} chunks  elapsed {el / 60:.1f} min  "
            f"eta {(len(todo) - done) * rate / 60:.1f} min")

    log(f"bootstrap done in {(time.time() - t0) / 60:.1f} min, arm {ARM}")
