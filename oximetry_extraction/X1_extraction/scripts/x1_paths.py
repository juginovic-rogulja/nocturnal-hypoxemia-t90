"""
The ONE spec module of lane X1 (v8 extraction). Every X1 script imports its paths from here
(gate 7 of the integrity rules): no other script names a table by literal path.
"""
from __future__ import annotations
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import os, sys
from pathlib import Path

T90 = paths.T90_ROOT
V8 = f"{paths.V8_ROOT}"
X1 = f"{V8}/X1_extraction"
X1_SCRIPTS = paths.X1_SCRIPTS
X1_WORK = f"{X1}/work"
X1_LOGS = f"{X1}/logs"
X1_DATA = f"{X1}/data"
X1_EC2 = paths.X1_EC2_SCRIPTS
REEX_CODE = paths.REEX_CODE               # step 002: the tarball's copy of the v7 extraction code
REEX_SCRIPTS = paths.REEX_SCRIPTS       # extract_night.py + reex_pipeline/
REEX_TARBALL = f"{paths.T90_ROOT}/EncodingB_Audit_2026-09-06/reextract_197/ec2/reex_code_as_uploaded.tar.gz"

V7D = f"{paths.PREV_TABLES_DIR}"       # PREV_DIR: read only, never written
RX = f"{paths.T90_ROOT}/EncodingB_Audit_2026-09-06/reextract_197"
NUMBERS = f"{paths.NUMBERS_DIR}"
V8_SMOKE = f"{V8}/smoke"
V8_LOGS = f"{V8}/logs"
PSV = f"{V8}/rerun_steps_v8.psv"

# Local mirror OUTSIDE iCloud Drive (sha256-verified copies, see logs/MIRROR.md). iCloud evicted reextract_FINAL.parquet and
# manifest_full_19173.csv on 2026-09-12 between 13:30 and 14:20 while this lane was running. A path resolves to the primary
# when it is healthy (blocks > 0) and to the mirror copy otherwise, never to an evicted file.
MIRROR = paths.MIRROR_ROOT


def _healthy(p):
    try:
        st = os.stat(p)
    except FileNotFoundError:
        return False
    return not (os.path.isfile(p) and st.st_size > 0 and getattr(st, "st_blocks", 1) == 0)


def _pick(primary):
    if _healthy(primary):
        return primary
    alt = primary.replace(T90, MIRROR)
    if _healthy(alt):
        sys.stderr.write(f"[x1_paths] {os.path.basename(primary)}: primary evicted, using the mirror copy {alt}\n")
        return alt
    return primary   # assert_not_evicted will stop the caller


CPAP_STAGE_V7 = _pick(f"{V7D}/cpap_t90_by_stage.parquet")
T90_FINAL_V7 = _pick(f"{V7D}/t90_final.parquet")
REEXTRACT_FINAL = _pick(f"{RX}/data/reextract_FINAL.parquet")
MANIFEST_FULL_V7 = _pick(f"{RX}/data/manifest_full_19173.csv")

# lane outputs (X1 layout; the psv's extract_v8/ paths map onto these, see X1_EC2_READY.md)
MANIFEST_SPLIT = f"{X1_WORK}/manifest_split_3622.csv"
MANIFEST_FULL_V8 = f"{X1_WORK}/manifest_full_19173_v8.csv"
MANIFEST_PROBE12 = f"{X1_WORK}/manifest_probe12.csv"
SAMPLE_300 = f"{X1_WORK}/sample_300.csv"
SAMPLE_300_SPLIT = f"{X1_WORK}/sample_300_split.csv"
HB_CALIBRATION = f"{X1_WORK}/hb_window_calibration.json"
HB_DEFINITION = f"{X1_WORK}/hb_definition.json"

S3_BUCKET = paths.S3_OUTPUT_BUCKET
S3_AP = paths.BDSP_S3_ACCESS_POINT
S3_PSG_BASE = f"s3://{S3_AP}/PSG/bids"
S3_CODE_KEY = "code/v8/reex_code_v8.tar.gz"
S3_INPUT_PREFIX = "inputs/v8"
S3_OUTPUT_PREFIX = "outputs/v8_2026-09"
AWS_PROFILE = os.environ.get("X1_AWS_PROFILE", "default")
AWS_REGION = paths.AWS_REGION

SEED = 20260912
COHORT_N_EXPECTED_SPLIT = 3622           # smoke-test count from the plan (a count, not a result)
def _probe_ids(path):
    """The twelve probe recordings of the v8 extraction (T90_PROBE_IDS_CSV, a csv with a BDSPPatientID column, a data file
    outside the repository); [] when the file is absent."""
    if not os.path.isfile(path):
        return []
    import csv
    with open(path) as f:
        return [int(r["BDSPPatientID"]) for r in csv.DictReader(f)]


PROBE12 = _probe_ids(paths.PROBE_IDS_CSV)

PY = paths.PY


def blocks(path) -> int:
    """Allocated blocks (0 with a non-zero size = iCloud-evicted)."""
    st = os.stat(path)
    return int(getattr(st, "st_blocks", -1))


def assert_not_evicted(*paths):
    bad = []
    for p in paths:
        p = str(p)
        if not os.path.exists(p):
            bad.append(f"MISSING {p}"); continue
        if os.path.isfile(p) and os.path.getsize(p) > 0 and blocks(p) == 0:
            bad.append(f"EVICTED {p}")
    if bad:
        raise SystemExit("input check failed:\n  " + "\n  ".join(bad))


def cohort_frame(columns=None):
    """The v7 cohort frame through numbers/cohort_spec.apply_cohort (the one cohort definition)."""
    import pandas as pd
    if NUMBERS not in sys.path:
        sys.path.insert(0, paths.ANALYSIS_DIR)
    from cohort_spec import apply_cohort  # noqa
    base = ["BDSPPatientID", "fu_valid", "spo2_pct_below_90", "oximetry_bad"]
    cols = base if columns is None else base + [c for c in columns if c not in base]
    assert_not_evicted(T90_FINAL_V7)
    t = pd.read_parquet(T90_FINAL_V7, columns=cols)
    return apply_cohort(t)


def ensure_dirs():
    for d in (X1_WORK, X1_LOGS, X1_DATA, X1_EC2, V8_SMOKE, V8_LOGS):
        Path(d).mkdir(parents=True, exist_ok=True)
