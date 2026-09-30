"""
X2 lane spec module: the ONE place every X2 script takes its paths from (integrity gate 7).

numbers/ and data_frozen_v7_2026-09/ are READ ONLY for this lane. cohort_spec is imported from numbers/ unchanged;
PREV_DIR is its V7_DIR (the plan's cohort_spec.PREV_DIR does not exist yet because numbers/ is not touched here).
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import datetime
import hashlib
import json
import os
import platform
import sys

ROOT = paths.T90_ROOT
NUMBERS = f"{paths.NUMBERS_DIR}"
if NUMBERS not in sys.path:
    sys.path.insert(0, paths.ANALYSIS_DIR)
import cohort_spec  # noqa: E402  (read-only import)

COHORT_N = cohort_spec.COHORT_N
apply_cohort = cohort_spec.apply_cohort
PREV_DIR = getattr(cohort_spec, "V7_DIR", None) or cohort_spec.PREV_DIR   # v7 frozen tables, read only (the v8 spec names it PREV_DIR)
PREV_T90_FINAL = f"{PREV_DIR}/t90_final.parquet"
PREV_SHA256_TXT = f"{PREV_DIR}/sha256.txt"
V7_DEFINITIONS = paths.V7_DISEASE_DEFINITIONS   # v8.1: the live numbers copy is v8 now; the v7 module is E1's backup

# the local hospital-record cache the v7 first-dates were built from (plan, Phase D)
OMOP_DIR = paths.OMOP_CACHE_DIR
CACHE = f"{OMOP_DIR}/condition_occurrence_cohort.parquet"
PSG_META = {"I0002": f"{OMOP_DIR}/psg_metadata/I0002_psg_metadata_2025-09-08.csv",
            "I0006": f"{OMOP_DIR}/psg_metadata/I0006_psg_metadata_2025-07-09.csv"}
ICD10_ERA_START = "2015-10-01"   # US ICD-10-CM transition; used by the audit only, never by the build rule

V8 = f"{paths.V8_ROOT}"
X2 = f"{V8}/X2_outcomes"
SCRIPTS = paths.X2_SCRIPTS
WORK = f"{X2}/work"
LOGS = f"{X2}/logs"
V8_DEFINITIONS = f"{paths.X2_SCRIPTS}/disease_definitions_v8.py"

OUT_PARQUET = f"{WORK}/outcomes_v8.parquet"
OUT_COUNTS = f"{WORK}/event_counts_v8.csv"
OUT_RECONCILE = f"{WORK}/reconcile_v7_lists_vs_v7_table.csv"
OUT_OLD_VS_NEW = f"{WORK}/old_vs_new_v8.csv"
OUT_INVARIANTS = f"{WORK}/invariants_v8.json"
DIFF_MD = f"{WORK}/DISEASE_DEFINITIONS_DIFF.md"
REPORT_MD = f"{LOGS}/OUTCOMES_V8.md"

OUTCOME_SUFFIXES = ("_first_date", "_prevalent", "_incident", "_years")


def blocks(path):
    st = os.stat(path)
    return st.st_blocks, st.st_size


def require_local(path):
    """iCloud eviction guard: 0 allocated blocks with a non-zero size means dataless. Never read or write onto it."""
    b, s = blocks(path)
    if s > 0 and b == 0:
        raise SystemExit(f"EVICTED (0 blocks, {s:,} bytes), not read: {path}")
    return path


def sha256(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for blk in iter(lambda: f.read(chunk), b""):
            h.update(blk)
    return h.hexdigest()


def write_sidecar(output, inputs, script, extra=None):
    """<output>.provenance.json: sha256 of the output, of every input and of the script, time, environment (gate 8)."""
    import duckdb
    import pandas as pd
    side = {
        "output": output, "output_sha256": sha256(output),
        "written": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
        "script": script, "script_sha256": sha256(script),
        "inputs": {k: {"path": v, "sha256": sha256(v), "st_blocks": blocks(v)[0]} for k, v in inputs.items()},
        "environment": {"python": sys.version.split()[0], "pandas": pd.__version__, "duckdb": duckdb.__version__,
                        "host": platform.node(), "platform": platform.platform()},
        "lane": "X2_outcomes (V8 steps 040, 041)",
    }
    if extra:
        side.update(extra)
    with open(output + ".provenance.json", "w") as f:
        json.dump(side, f, indent=1, default=str)
    return side
