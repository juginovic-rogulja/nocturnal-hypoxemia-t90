"""
paths.py: the one place every script of this repository takes its locations from.

Every script imports this module (a two-line bootstrap at the top of each file puts the
repository root on sys.path) and names no absolute path of its own. The values below are
read from environment variables, each with a documented default. Nothing here creates a
folder, touches sys.path or reads a file: importing paths has no side effect.

Two kinds of location are distinguished.

CODE locations are fixed relative to this file and are where the scripts of the repository
live. Scripts that import a sibling module by folder (sys.path.insert) use these.

DATA and WORK locations are where inputs are read and outputs are written. They default
to folders that mirror the layout the analysis ran in, so that a script's relative
references (numbers/, data_frozen_v8_2026-09/, the round folders of the figures) keep
working once the inputs are placed there. Override any of them with the environment
variable named beside it.

    T90_ROOT            the data tree: frozen tables, extraction work, caches
    T90_NUMBERS_DIR     the numbers files (analysis outputs); default beside the analysis scripts
    T90_TABLES_DIR      the four frozen v8 tables (t90_final.parquet, t90_base4.parquet,
                        master_cohort.csv, cpap_t90_by_stage.parquet)
    T90_PREV_TABLES_DIR the v7 tables the v8 builders read (old versus new)
    T90_V8_ROOT         the v8 extraction lane work tree (X1 data/work/logs, X2 work, ...)
    T90_SV_ROOT         outputs of the Sleep_Variability analyses; default beside their scripts
    T90_FIGURE_ROOT     the figure tree (round folders with their scripts, base sheets and
                        lane outputs); default the repository's figures folder
    T90_PSG_DIR         local BDSP .h5 recordings (the extraction reads them here or on S3)
    T90_OMOP_CACHE_DIR  the hospital-record cache (condition_occurrence_cohort.parquet, psg_metadata/)
    T90_SHHS_DIR        the NSRR Sleep Heart Health Study files
    T90_MROS_DIR        the NSRR Osteoporotic Fractures in Men Study files
    T90_MIRROR_ROOT     an optional second copy of the data tree (the extraction lane falls
                        back to it when a primary file is missing); empty = no mirror
    T90_FEATURE_WORK    working folder of the April 2026 feature extraction (merge_v5.py)
    T90_V7_DISEASE_DEFINITIONS  the v7 disease_definitions module (the outcome rebuild proves
                        the untouched lists against it); a reviewer supplies it or skips that gate
    T90_HB_DEFINITION   hb_definition.json (the frozen hypoxic-burden windows per site)
    T90_PROBE_IDS_CSV   the twelve probe recordings of the v8 extraction (a csv with a
                        BDSPPatientID column); a data file, not part of the repository
    T90_OLD_SNAPSHOT_SV the v7 snapshot of the Sleep_Variability outputs (old side of the
                        old-versus-new deltas some figure verifiers print)

    BDSP_S3_ACCESS_POINT  the S3 access point ARN of the credentialed BDSP bucket (no default)
    T90_S3_OUTPUT_BUCKET  the S3 bucket the EC2 workers write to (no default)
    T90_S3_ATLAS_RESUME_BUCKET, T90_S3_PAPER11_BUCKET  the two buckets merge_v5.py read
    AWS_REGION            default us-east-1

    T90_FIGURELAB_DIR   FigureLab (bioglyph), the separate drawing tool the schematic bands
                        were drawn with (not part of this repository)
    T90_PRISM_PLOTTER_SRC  the prism_plotter package (a separate plotting package, not part
                        of this repository) imported by some figure builders
    T90_PY, T90_GS, T90_CHROME, T90_NODE   interpreter and renderer executables
    T90_FONT_DIR        folder holding Arial.ttf, Arial Bold.ttf, Arial Italic.ttf
    T90_DOCX_VALIDATE   optional docx validator script used by the figure assembly
    T90_TESSERACT       the tesseract executable (optical character checks of rendered sheets)
"""
import os
import sys

REPO_ROOT = os.path.dirname(os.path.abspath(__file__))

# ------------------------------------------------------------------ code locations (fixed)
SIGNAL_FEATURES_DIR = os.path.join(REPO_ROOT, "signal_features")
OXIMETRY_DIR = os.path.join(REPO_ROOT, "oximetry_extraction")
X1_DIR = os.path.join(OXIMETRY_DIR, "X1_extraction")
X1_SCRIPTS = os.path.join(X1_DIR, "scripts")
X1_EC2_SCRIPTS = os.path.join(X1_DIR, "ec2")
REEX_CODE = os.path.join(X1_DIR, "reex_code")
REEX_SCRIPTS = os.path.join(REEX_CODE, "scripts")
X2_SCRIPTS = os.path.join(OXIMETRY_DIR, "X2_outcomes", "scripts")
E2_SCRIPTS = os.path.join(OXIMETRY_DIR, "E2_tables", "scripts")
E1_SCRIPTS = os.path.join(OXIMETRY_DIR, "E1_live", "scripts")
X4_DIR = os.path.join(OXIMETRY_DIR, "X4_groupP")
CPAP_STAGE_DIR = os.path.join(OXIMETRY_DIR, "cpap_stage")
ANALYSIS_DIR = os.path.join(REPO_ROOT, "analysis", "numbers")        # cohort_spec.py, disease_definitions.py live here
OTHER_WRITERS_DIR = os.path.join(REPO_ROOT, "analysis", "other_writers")
SV_CODE_DIR = os.path.join(REPO_ROOT, "analysis", "sleep_variability")
SLEEP_ONLY_DIR = os.path.join(REPO_ROOT, "analysis", "sleep_only_t90")
FIGURES_DIR = os.path.join(REPO_ROOT, "figures")


def _env(name, default):
    v = os.environ.get(name)
    return os.path.expanduser(v) if v else default


# ------------------------------------------------------------------ data and work locations (environment)
T90_ROOT = _env("T90_ROOT", os.path.expanduser("~/T90_work"))
NUMBERS_DIR = _env("T90_NUMBERS_DIR", ANALYSIS_DIR)
TABLES_DIR = _env("T90_TABLES_DIR", os.path.join(T90_ROOT, "data_frozen_v8_2026-09"))
PREV_TABLES_DIR = _env("T90_PREV_TABLES_DIR", os.path.join(T90_ROOT, "data_frozen_v7_2026-09"))
V8_ROOT = _env("T90_V8_ROOT", os.path.join(T90_ROOT, "V8_RECALC_2026-09-12"))
SV_ROOT = _env("T90_SV_ROOT", SV_CODE_DIR)
FIGURE_ROOT = _env("T90_FIGURE_ROOT", FIGURES_DIR)
PSG_DIR = _env("T90_PSG_DIR", os.path.join(T90_ROOT, "psg"))
OMOP_CACHE_DIR = _env("T90_OMOP_CACHE_DIR", os.path.join(T90_ROOT, "omop_cache"))
SHHS_DIR = _env("T90_SHHS_DIR", os.path.join(T90_ROOT, "nsrr", "shhs"))
MROS_DIR = _env("T90_MROS_DIR", os.path.join(T90_ROOT, "nsrr", "mros"))
MIRROR_ROOT = _env("T90_MIRROR_ROOT", "")
FEATURE_WORK = _env("T90_FEATURE_WORK", os.path.join(T90_ROOT, "feature_pipeline"))
FEATURE_OUTPUTS_DIR = os.path.join(FEATURE_WORK, "outputs")
V7_DISEASE_DEFINITIONS = _env("T90_V7_DISEASE_DEFINITIONS", os.path.join(T90_ROOT, "v7", "disease_definitions_v7.py"))
HB_DEFINITION = _env("T90_HB_DEFINITION", os.path.join(V8_ROOT, "X1_extraction", "work", "hb_definition.json"))
PROBE_IDS_CSV = _env("T90_PROBE_IDS_CSV", os.path.join(V8_ROOT, "X1_extraction", "work", "probe12_ids.csv"))
OLD_SNAPSHOT_SV = _env("T90_OLD_SNAPSHOT_SV", os.path.join(V8_ROOT, "compare", "old_numbers", "Sleep_Variability_2026-08"))   # the v7 snapshot of the Sleep_Variability outputs (old side of the old-versus-new deltas)

# ------------------------------------------------------------------ cloud storage (environment, no defaults for the private names)
BDSP_S3_ACCESS_POINT = os.environ.get("BDSP_S3_ACCESS_POINT", "")     # the access point identifier the platform issues (region, account and name), no default
S3_PSG_BASE = f"s3://{BDSP_S3_ACCESS_POINT}/PSG/bids"
S3_OUTPUT_BUCKET = os.environ.get("T90_S3_OUTPUT_BUCKET", "")
S3_ATLAS_RESUME_BUCKET = os.environ.get("T90_S3_ATLAS_RESUME_BUCKET", "")
S3_PAPER11_BUCKET = os.environ.get("T90_S3_PAPER11_BUCKET", "")
AWS_REGION = os.environ.get("AWS_REGION", "us-east-1")

# ------------------------------------------------------------------ external tools (environment)
FIGURELAB_DIR = _env("T90_FIGURELAB_DIR", os.path.join(T90_ROOT, "FigureLab"))
PRISM_PLOTTER_SRC = _env("T90_PRISM_PLOTTER_SRC", os.path.join(T90_ROOT, "prism_plotter", "src"))
PY = _env("T90_PY", sys.executable or "python3")
GS = _env("T90_GS", "gs")
CHROME = _env("T90_CHROME", "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome")
NODE = _env("T90_NODE", "node")
FONT_DIR = _env("T90_FONT_DIR", "/System/Library/Fonts/Supplemental")
FONT_ARIAL = os.path.join(FONT_DIR, "Arial.ttf")
FONT_ARIAL_BOLD = os.path.join(FONT_DIR, "Arial Bold.ttf")
FONT_ARIAL_ITALIC = os.path.join(FONT_DIR, "Arial Italic.ttf")
DOCX_VALIDATE = _env("T90_DOCX_VALIDATE", "")
TESSERACT = _env("T90_TESSERACT", "tesseract")
HOME = os.path.expanduser("~")

_EXPORTED = ["T90_ROOT", "NUMBERS_DIR", "TABLES_DIR", "PREV_TABLES_DIR", "V8_ROOT", "SV_ROOT", "FIGURE_ROOT", "PSG_DIR",
             "OMOP_CACHE_DIR", "SHHS_DIR", "MROS_DIR", "MIRROR_ROOT", "FEATURE_WORK", "V7_DISEASE_DEFINITIONS",
             "HB_DEFINITION", "PROBE_IDS_CSV", "BDSP_S3_ACCESS_POINT", "S3_OUTPUT_BUCKET", "S3_ATLAS_RESUME_BUCKET",
             "S3_PAPER11_BUCKET", "AWS_REGION", "FIGURELAB_DIR", "PRISM_PLOTTER_SRC", "PY", "GS", "CHROME", "NODE",
             "FONT_DIR", "DOCX_VALIDATE", "TESSERACT", "REPO_ROOT", "SIGNAL_FEATURES_DIR", "OXIMETRY_DIR", "X1_SCRIPTS", "REEX_SCRIPTS",
             "X2_SCRIPTS", "E2_SCRIPTS", "ANALYSIS_DIR", "SV_CODE_DIR", "FIGURES_DIR"]


def as_shell_exports():
    """Lines 'export T90_<NAME>=...' for the shell scripts of the repository (python paths.py)."""
    lines = []
    for k in _EXPORTED:
        name = k if k.startswith(("T90_", "BDSP_", "AWS_")) else f"T90_{k}"
        v = globals()[k]
        lines.append(f"export {name}={_sh_quote(str(v))}")
    return "\n".join(lines)


def _sh_quote(s):
    return "'" + s.replace("'", "'\\''") + "'"


if __name__ == "__main__":
    print(as_shell_exports())
