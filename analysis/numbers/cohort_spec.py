"""
The one place the analysis cohort, the frozen tables and the v8 counts are defined.

v8 (2026-09-12, lane X3). Step 053 installs this file as numbers/cohort_spec.py (backup
numbers/cohort_spec_PRE_V8.py). Every v8 script imports its table paths, its counts and its
smoke hooks from here and from nowhere else. The inherited scripts are repointed to DATA_DIR
by repoint_v8.py; the new scripts name no path of their own.

What changed against the v7 module
  DATA_DIR             data_frozen_v8_2026-09, the only table folder a chain step may read
  PREV_DIR             data_frozen_v7_2026-09, read ONLY by the v8 builders (steps 041, 050,
                       240) to carry untouched cells forward and to print old-versus-new
  V7_DIR               removed on purpose: an importer that still names it fails at import
  T90_FINAL, T90_BASE4, MASTER_CSV, CPAP_STAGE_PARQUET   the four tables under DATA_DIR
  N_MEASURES 141       137 non-siteZ v7 measures + hypoxic_burden, plm_index, sol_min,
                       rem_latency_min (decision 3)
  N_RANKED_OUTCOMES 52 48 + polycythemia, Parkinson's disease, interstitial lung disease,
                       ventricular arrhythmia or cardiac arrest (decision 11)
  SENSITIVITY_ONLY     spo2_pct_below_90_sleep: never a row of the primary ranking, it is
                       the exposure of the decision 2 sensitivity set (T90_COLUMN)
  RANKING_DROP_SUFFIX  _siteZ: a master column with this suffix is never a candidate
  PREVALENT_CARDIOPULMONARY   the decision 6 covariate set (D4 default, hypertension out)
  smoke hooks          T90_SMOKE_IDS, T90_DATA_DIR, T90_SMOKE_OUT, T90_COLUMN (below)

Smoke hooks (environment; every one prints a [SMOKE] banner once per process)
  T90_SMOKE_IDS=<csv with a BDSPPatientID column>   apply_cohort returns that subset and
                       prints [SMOKE] instead of asserting COHORT_N (plan section 10 item 15)
  T90_DATA_DIR=<dir>   read the tables from another folder (the v7 tables, for a smoke run
                       before the v8 tables exist). SMOKE_DATA becomes True.
  T90_SMOKE_OUT=<dir>  write outputs there instead of numbers/ (OUT_DIR). Only the new v8
                       scripts and the sweep-edited writers honour it.
  T90_COLUMN=<column>  the exposure column for the decision 2 sensitivity set (steps 230,
                       231). Default PRIMARY_EXPOSURE. exposure_tag() gives the file suffix.
The driver's preflight and repoint_v8.py --check REFUSE to run while any T90_SMOKE_* or
T90_DATA_DIR variable is set, so a smoke setting cannot leak into the real chain.

The oximetry quality rule (2026-08-04) is unchanged: 210 recordings whose stored oxygen is
not a measurement of that patient (mean saturation 50 to 80 percent, more than 70 percent of
the recording below 90, ODI exactly zero) are removed. Cohort 19,173.

v8.1 (Alen, 14 Sept 2026): a split-night recording is half real sleep and half sleep on a mask, so the 3,622 split nights
contribute their pre-CPAP oxygen (and event-rate) measures only. Their sleep amounts, stage composition, REM latency and event
counts are set missing in the v8.1 tables (E2 builder, column _v8_sleep_masked), and the ranking of measurements and every
sleep-duration analysis run on the 15,551 full diagnostic nights: the psv command of those steps sets T90_FULL_NIGHTS_ONLY=1,
apply_cohort then drops the flagged rows and COHORT_N is COHORT_N_FULL. The T90 hazard ratios and the PAP analyses keep all
19,173 nights. Straight swap of two negative controls on the same date: see disease_definitions.SWAPPED_CONTROLS_2026_09_14.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../..")))  # repository root, where paths.py lives
import paths
import hashlib
import json
import os
import sys
import time

COHORT_N_ALL = 19173                                  # every analysis night of the v8 tables
COHORT_N_FULL = 15551                                 # v8.1: the full diagnostic nights (the 3,622 split-night recordings excluded)
SPLIT_FLAG = "_v8_prepap_window"                      # E2 provenance column, true on the split-night recordings
FULL_NIGHTS_ONLY = bool(os.environ.get("T90_FULL_NIGHTS_ONLY"))   # v8.1: set in the psv command of the ranking family and of every sleep-duration analysis
COHORT_N = COHORT_N_FULL if FULL_NIGHTS_ONLY else COHORT_N_ALL

ROOT = paths.T90_ROOT
NUMBERS_DIR = f"{paths.NUMBERS_DIR}"
PREV_DIR = f"{paths.PREV_TABLES_DIR}"          # read only by the v8 builders
_DEFAULT_DATA_DIR = f"{paths.TABLES_DIR}"

# ------------------------------------------------------------------ counts and lists (decisions 3, 11, 2, 6)
N_MEASURES = 141
N_RANKED_OUTCOMES = 52
NEW_MEASURES_V8 = ("hypoxic_burden", "plm_index", "sol_min", "rem_latency_min")
PRIMARY_EXPOSURE = "spo2_pct_below_90"
SENSITIVITY_ONLY = {"spo2_pct_below_90_sleep"}
SENSITIVITY_TAG = {"spo2_pct_below_90_sleep": "sleepT90"}
RANKING_DROP_SUFFIX = ("_siteZ",)
# decision 6, D4 default: any prevalent flag among these twelve. Hypertension (htn2) is left
# out because about half the cohort carries it; the alternative (adding it) is listed in the
# plan and is one edit here.
PREVALENT_CARDIOPULMONARY = ["hf", "ihd", "mi", "afib", "stroke_any", "pad", "pulm_htn", "vte",
                             "copd2", "asthma", "resp_failure", "obesity_hypovent"]
# measure_families.csv ships with the repository beside this file, so it is addressed through the CODE
# location (paths.ANALYSIS_DIR) and not through NUMBERS_DIR, which a reviewer repoints to their own
# output folder. Override with T90_MEASURE_FAMILIES_CSV to supply a different family table.
MEASURE_FAMILIES_CSV = os.environ.get("T90_MEASURE_FAMILIES_CSV") or f"{paths.ANALYSIS_DIR}/measure_families.csv"

# ------------------------------------------------------------------ smoke hooks
_BANNERS = set()


def _banner(msg):
    if msg not in _BANNERS:
        _BANNERS.add(msg)
        print(f"[SMOKE] {msg}", file=sys.stderr, flush=True)


SMOKE_ENV = ("T90_SMOKE_IDS", "T90_DATA_DIR", "T90_SMOKE_OUT")


def smoke_env_active():
    """The names of every smoke variable that is set. The real chain refuses to start on any."""
    return [k for k in SMOKE_ENV if os.environ.get(k)]


def _configure(data_dir=None, out_dir=None):
    """Set every derived name from one data folder and one output folder. Called once at import
    from the environment; a test harness may call it again to point at another folder."""
    g = globals()
    g["DATA_DIR"] = data_dir or os.environ.get("T90_DATA_DIR") or _DEFAULT_DATA_DIR
    g["SMOKE_DATA"] = g["DATA_DIR"] != _DEFAULT_DATA_DIR
    g["T90_FINAL"] = f"{DATA_DIR}/t90_final.parquet"
    g["T90_BASE4"] = f"{DATA_DIR}/t90_base4.parquet"
    g["MASTER_CSV"] = f"{DATA_DIR}/master_cohort.csv"
    g["CPAP_STAGE_PARQUET"] = f"{DATA_DIR}/cpap_t90_by_stage.parquet"
    g["DATA_PROVENANCE_MD"] = f"{DATA_DIR}/PROVENANCE.md"
    g["TABLES"] = (T90_FINAL, T90_BASE4, CPAP_STAGE_PARQUET, MASTER_CSV)
    g["OUT_DIR"] = out_dir or os.environ.get("T90_SMOKE_OUT") or NUMBERS_DIR
    g["SMOKE_OUT"] = g["OUT_DIR"] != NUMBERS_DIR
    if SMOKE_DATA:
        _banner(f"DATA_DIR overridden to {DATA_DIR} (not the v8 tables)")
    if SMOKE_OUT:
        _banner(f"OUT_DIR overridden to {OUT_DIR} (nothing is written into numbers/)")


_configure()


def t90_column():
    """The exposure column: the primary T90, or the T90_COLUMN sensitivity exposure."""
    return os.environ.get("T90_COLUMN") or PRIMARY_EXPOSURE


def exposure_tag(col=None):
    """File-name suffix for the exposure: '' for the primary, '_sleepT90' for the sensitivity set."""
    col = col or t90_column()
    if col == PRIMARY_EXPOSURE:
        return ""
    return "_" + SENSITIVITY_TAG.get(col, col)


def out_path(name):
    """Where an output file goes: numbers/ in the chain, the smoke folder in a smoke run."""
    os.makedirs(OUT_DIR, exist_ok=True)
    return os.path.join(OUT_DIR, name)


def smoke_ids():
    p = os.environ.get("T90_SMOKE_IDS")
    if not p:
        return None
    import pandas as pd
    ids = pd.read_csv(p)
    col = "BDSPPatientID" if "BDSPPatientID" in ids.columns else ids.columns[0]
    return set(ids[col].tolist())


def apply_cohort(b):
    """Filter a loaded cohort frame to the analysis population and assert the expected size.
    With T90_SMOKE_IDS set, return that subset instead and print [SMOKE]."""
    if "oximetry_bad" not in b.columns:
        raise KeyError("oximetry_bad missing: the frozen table must carry the oximetry quality flag")
    out = b[(b.fu_valid == 1) & b.spo2_pct_below_90.notna() & (b.oximetry_bad == 0)].copy()
    if FULL_NIGHTS_ONLY:                                  # v8.1 (14 Sept): split nights contribute oxygen only, never a sleep number
        if SPLIT_FLAG not in out.columns:
            raise KeyError(f"{SPLIT_FLAG} missing: the frozen table must carry the split-night flag")
        out = out[~out[SPLIT_FLAG].astype(bool)].copy()
    ids = smoke_ids()
    if ids is not None:
        out = out[out.BDSPPatientID.isin(ids)].copy()
        _banner(f"apply_cohort returns the {len(out):,}-night subset of T90_SMOKE_IDS, COHORT_N not asserted")
        return out
    assert len(out) == COHORT_N, f"cohort is {len(out):,}, expected {COHORT_N:,}"
    return out


# ------------------------------------------------------------------ provenance sidecar
MIN_RANKED_EVENTS = 150                               # the paper's rule: a condition is ranked when it has at least this many incident events in the cohort


def all_nights_cohort(b):
    """The analysis cohort WITHOUT the v8.1 split-night drop (COHORT_N_ALL rows): the base rule only."""
    if "oximetry_bad" not in b.columns:
        raise KeyError("oximetry_bad missing: the frozen table must carry the oximetry quality flag")
    out = b[(b.fu_valid == 1) & b.spo2_pct_below_90.notna() & (b.oximetry_bad == 0)].copy()
    if smoke_ids() is None:
        assert len(out) == COHORT_N_ALL, f"all-nights cohort is {len(out):,}, expected {COHORT_N_ALL:,}"
    return out


def ranked_outcome_keys(tbl, candidates):
    """v8.1 (14 Sept): the ranked outcome SET is a property of the study, fixed on the all-nights cohort (at least
    MIN_RANKED_EVENTS incident events among the COHORT_N_ALL nights, the Methods rule), whatever cohort the ranking's
    concordance is then computed on. Under T90_FULL_NIGHTS_ONLY four conditions (cirrhosis, Parkinson's disease, inguinal
    hernia, polycythemia) have 132 to 148 events on the 15,551 full nights; they stay ranked. tbl is the frozen table as read."""
    full = all_nights_cohort(tbl)
    return [k for k in candidates if f"{k}_incident" in full.columns and int(full[f"{k}_incident"].sum()) >= MIN_RANKED_EVENTS]


def drop_split_nights_if_full(b):
    """v8.1: with T90_FULL_NIGHTS_ONLY set, drop the split-night recordings (SPLIT_FLAG true); otherwise return b unchanged.
    For the few scripts that filter the cohort by hand instead of calling apply_cohort."""
    if not FULL_NIGHTS_ONLY:
        return b
    if SPLIT_FLAG not in b.columns:
        raise KeyError(f"{SPLIT_FLAG} missing: the frozen table must carry the split-night flag")
    return b[~b[SPLIT_FLAG].astype(bool)].copy()


def cohort_tag():
    """v8.1: '_fullnights' when T90_FULL_NIGHTS_ONLY is set (for outputs that must not overwrite the all-nights file), else ''."""
    return "_fullnights" if FULL_NIGHTS_ONLY else ""


def sha256_of(path, bufsize=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(bufsize), b""):
            h.update(chunk)
    return h.hexdigest()


_SHA_CACHE = {}


def describe(path):
    """path, bytes, mtime, sha256 (hashed once per process per (path, size, mtime))."""
    if not os.path.exists(path):
        return {"path": path, "missing": True}
    st = os.stat(path)
    key = (path, st.st_size, st.st_mtime_ns)
    if key not in _SHA_CACHE:
        _SHA_CACHE[key] = sha256_of(path)
    return {"path": path, "bytes": st.st_size,
            "mtime_local": time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(st.st_mtime)),
            "sha256": _SHA_CACHE[key]}


def sidecar(output, script, extra_inputs=(), note="", record=None):
    """Write <output>.provenance.json beside an output: its sha256, the script's sha256, the
    four tables under DATA_DIR and DATA_DIR/PROVENANCE.md, every extra input, the smoke state.
    The driver's stamper writes the same keys plus the step block after the step; the two
    records agree on every shared field."""
    rec = {"output": describe(output),
           "script": describe(os.path.abspath(script)),
           "inputs": {os.path.basename(p): describe(p) for p in list(TABLES) + [DATA_PROVENANCE_MD]},
           "extra_inputs": {os.path.basename(p): describe(p) for p in extra_inputs},
           "spec": {"module": os.path.abspath(__file__), "DATA_DIR": DATA_DIR, "OUT_DIR": OUT_DIR,
                    "COHORT_N": COHORT_N, "N_MEASURES": N_MEASURES, "N_RANKED_OUTCOMES": N_RANKED_OUTCOMES,
                    "T90_COLUMN": t90_column(), "smoke_env": {k: os.environ.get(k) for k in SMOKE_ENV if os.environ.get(k)}},
           "environment": {"python": sys.version.split()[0], "user": os.environ.get("USER", ""),
                           "written": time.strftime("%Y-%m-%d %H:%M:%S")},
           "note": note}
    if record:
        rec["record"] = record
    side = output + ".provenance.json"
    json.dump(rec, open(side, "w"), indent=1, default=str)
    return side
