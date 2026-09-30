"""
The ONE spec module of lane E2 (Phase E, steps 050 to 052). Every E2 script takes its paths from here (integrity gate 7).

The cohort, the counts (N_MEASURES, N_RANKED_OUTCOMES, COHORT_N), PREV_DIR (the v7 tables, read only) and the real output folder
come from numbers/cohort_spec.py (v8, installed by lane E1). The lane inputs (X1 assembled passes, X2 outcomes) are named here.
iCloud rule: a file with 0 allocated blocks and a non-zero size is evicted; it is never read and never written onto. A primary
input that is evicted resolves to its copy under the local mirror when one exists, otherwise the caller stops.
"""
from __future__ import annotations
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import hashlib, json, os, platform, sys, time

ROOT = paths.T90_ROOT
NUMBERS = f"{paths.NUMBERS_DIR}"
if NUMBERS not in sys.path:
    sys.path.insert(0, paths.ANALYSIS_DIR)
import cohort_spec as CS  # noqa: E402  (v8 module: PREV_DIR, DATA_DIR, COHORT_N, N_MEASURES, apply_cohort, NEW_MEASURES_V8 ...)

V8 = f"{paths.V8_ROOT}"
E2 = f"{V8}/E2_tables"
E2_SCRIPTS, E2_WORK, E2_LOGS = paths.E2_SCRIPTS, f"{E2}/work", f"{E2}/logs"
X1 = f"{V8}/X1_extraction"
X1_SCRIPTS, X1_DATA, X1_WORK = paths.X1_SCRIPTS, f"{X1}/data", f"{X1}/work"
X2 = f"{V8}/X2_outcomes"
X2_SCRIPTS, X2_WORK = paths.X2_SCRIPTS, f"{X2}/work"
AUD = f"{paths.T90_ROOT}/EncodingB_Audit_2026-09-06"
MIRROR = paths.MIRROR_ROOT

PREV_DIR = CS.PREV_DIR                       # data_frozen_v7_2026-09, read only
REAL_OUT = CS._DEFAULT_DATA_DIR              # data_frozen_v8_2026-09, written only by the real build (step 050)
PREV = {"t90_final": f"{PREV_DIR}/t90_final.parquet", "t90_base4": f"{PREV_DIR}/t90_base4.parquet",
        "master": f"{PREV_DIR}/master_cohort.csv", "cpap": f"{PREV_DIR}/cpap_t90_by_stage.parquet",
        "sha256": f"{PREV_DIR}/sha256.txt", "changelog": f"{PREV_DIR}/CHANGELOG.parquet", "provenance": f"{PREV_DIR}/PROVENANCE.md"}
TABLE_FILES = ("t90_final.parquet", "t90_base4.parquet", "master_cohort.csv", "cpap_t90_by_stage.parquet")

# X1 (Phase A and B) outputs. The assembled parquets are the inputs of the real build; the jsonl passes are accepted too
# (the builder then applies X1's assemble_common rules itself). Candidates in order of preference.
SPLIT_ASSEMBLED_CANDIDATES = [f"{X1_DATA}/split3622_assembled.parquet", f"{X1_DATA}/split_prepap_FINAL.parquet"]
LIGHT_ASSEMBLED_CANDIDATES = [f"{X1_DATA}/light19173_assembled.parquet", f"{X1_DATA}/light_FINAL.parquet"]
SPLIT_JSONL = f"{X1_DATA}/split3622_FINAL.jsonl"
LIGHT_JSONL = f"{X1_DATA}/light19173_FINAL.jsonl"
MANIFEST_SPLIT = f"{X1_WORK}/manifest_split_3622.csv"
MANIFEST_FULL = f"{X1_WORK}/manifest_full_19173_v8.csv"
HB_DEFINITION = f"{X1_WORK}/hb_definition.json"
# pilot equivalents (smoke)
PILOT_SPLIT_JSONL = f"{X1_DATA}/pilot300_split_FINAL.jsonl"
PILOT_SPLIT_ASSEMBLED = f"{X1_DATA}/pilot300_split_assembled.parquet"
PILOT_LIGHT_JSONL = f"{X1_DATA}/light_pilot300_frozen.jsonl"
PILOT_LIGHT_ASSEMBLED = f"{X1_DATA}/light_pilot300_assembled.parquet"
SMOKE_SAMPLE_300 = f"{V8}/smoke/sample_300.csv"

# X2 (Phase D) outputs
OUTCOMES_V8 = f"{X2_WORK}/outcomes_v8.parquet"
EVENT_COUNTS_V8 = f"{X2_WORK}/event_counts_v8.csv"
DEFS_V8 = f"{X2_SCRIPTS}/disease_definitions_v8.py"
# the hospital-record cache the first dates are built from (same literal as X2's x2_spec, which no longer imports: it names cohort_spec.V7_DIR)
OMOP_CACHE = (f"{paths.OMOP_CACHE_DIR}/"
              "condition_occurrence_cohort.parquet")
OUTCOME_SUFFIXES = ("_first_date", "_prevalent", "_incident", "_years")
# v8.1 (Alen, 14 Sept 2026): a split night is half real sleep and half sleep on a mask, so on the split nights the builder sets
# these sleep amounts, stage composition, REM latency and event counts to missing (column _v8_sleep_masked). Oxygen measures,
# AHI, arousal index, PLM index, sleep onset latency and every rate (densities, HR/HRV, spectral) stay.
SLEEP_MASK_COLS = ["TST_min", "recording_dur_min", "wake_min", "WASO_min", "N1_min", "N2_min", "N3_min", "REM_min",
                   "N1_pct", "N2_pct", "N3_pct", "REM_pct", "sleep_efficiency_pct", "rem_latency_min",
                   "F_spindle_n_spindles", "C_spindle_n_spindles", "O_spindle_n_spindles", "F_SO_n_SOs", "C_SO_n_SOs", "O_SO_n_SOs",
                   "wholenight_n_rr", "nrem_n_rr", "resp_events_any_label"]
SLEEP_MASK_REASON = "split night: sleep amounts, composition, REM latency and event counts set missing (v8.1, 14 Sept decision)"

# plausibility bounds (v7 step 5). Evicted on 2026-09-12 and absent from the mirror; the builder falls back to the bounds
# recorded in the v7 CHANGELOG reasons and lists the bounded columns it could not check.
# The plausibility bounds table (column, lo, hi). The audit lane's file is read when it exists; otherwise the copy
# shipped with the repository beside this folder, which is the same table. T90_BOUNDS_CSV overrides both.
_SHIPPED_BOUNDS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "plausibility_bounds.csv")
_AUD_BOUNDS = f"{AUD}/plausibility/impossible_strict_counts.csv"
BOUNDS_CSV = os.environ.get("T90_BOUNDS_CSV") or (_AUD_BOUNDS if os.path.exists(_AUD_BOUNDS) else _SHIPPED_BOUNDS)

STAGE_CODES_SLEEP = (1, 2, 3, 4)             # BDSP legend: 1=N3 2=N2 3=N1 4=REM 5=Wake 0=unscored (read from the file by X1)
MIN_TST_FOR_INDEX = 60.0                     # D1 default (v7 rule): AHI and arousal index undefined under 60 min of true sleep
MIN_OXYGEN_VALID_MIN = 5.0                   # D1 default: oxygen measures need 5 valid minutes of the pre-PAP window
COLLAPSE_MIN_TST_WINDOWED = 60.0             # a windowed night with a single non-REM code and under 60 min sleep is short, not collapsed

PY = paths.PY


def blocks(path):
    return int(getattr(os.stat(path), "st_blocks", -1))


def evicted(path):
    return os.path.isfile(path) and os.path.getsize(path) > 0 and blocks(path) == 0


def pick(primary):
    """The primary path when resident, else the mirror copy, else the primary (require_resident stops the caller)."""
    if os.path.exists(primary) and not evicted(primary):
        return primary
    alt = primary.replace(ROOT, MIRROR)
    if os.path.exists(alt) and not evicted(alt):
        sys.stderr.write(f"[e2_spec] {os.path.basename(primary)}: primary evicted or missing, using the mirror copy\n")
        return alt
    return primary


def require_resident(*paths):
    bad = []
    for p in paths:
        if not os.path.exists(p):
            bad.append(f"MISSING {p}")
        elif evicted(p):
            bad.append(f"EVICTED (0 blocks, {os.path.getsize(p):,} bytes) {p}")
    if bad:
        raise SystemExit("input check failed:\n  " + "\n  ".join(bad))


def first_existing(cands):
    for c in cands:
        if os.path.exists(c) and not evicted(c):
            return c
    return None


def sha256(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def write_sidecar(output, inputs, script, step, extra=None, index_path=None):
    """<output>.provenance.json: sha256 of the output, of every input and of the script, time, environment (gate 8)."""
    rec = {"output": output, "output_sha256": sha256(output) if os.path.isfile(output) else None,
           "output_bytes": os.path.getsize(output) if os.path.isfile(output) else None,
           "step": step, "script": script, "script_sha256": sha256(script) if os.path.isfile(script) else None,
           "inputs": {p: (sha256(p) if os.path.isfile(p) else "DIR_OR_MISSING") for p in inputs},
           "time_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "host": platform.node(),
           "env": {"python": sys.version.split()[0], "platform": platform.platform()},
           "spec": {"module": os.path.abspath(__file__), "PREV_DIR": PREV_DIR, "REAL_OUT": REAL_OUT, "COHORT_N": CS.COHORT_N,
                    "N_MEASURES": CS.N_MEASURES, "N_RANKED_OUTCOMES": CS.N_RANKED_OUTCOMES},
           "extra": extra or {}}
    for m in ("numpy", "pandas", "pyarrow", "duckdb", "scipy"):
        try:
            rec["env"][m] = __import__(m).__version__
        except Exception:
            rec["env"][m] = None
    with open(output + ".provenance.json", "w") as f:
        json.dump(rec, f, indent=1, default=str)
    if index_path:
        with open(index_path, "a") as f:
            f.write(f"{rec['time_utc']}\t{step}\t{output}\t{rec['output_sha256']}\t{os.path.basename(script)}\n")
    return rec
