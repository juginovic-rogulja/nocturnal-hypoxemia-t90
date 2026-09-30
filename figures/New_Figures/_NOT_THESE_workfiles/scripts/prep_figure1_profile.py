"""
Cohort medians for the measures a clinician uses to qualify a sleep study.

Alen's comment on main Figure 1 panel A: the panel currently shows only the prevalence of low
oxygen, and he wants the medians of the measures people actually read off a sleep report, sleep
efficiency, total sleep time and oxygenation, because a count of diagnoses is skewed by the
patients who carry many.

Nothing is modelled here. Every value is a median and an interquartile range of a stored column,
taken on the frozen analysis cohort defined in numbers/cohort_spec.py, which is the same
19,173 recordings every fit in the paper uses (fu_valid == 1, oxygen present, oximetry_bad == 0).

Total sleep time, sleep efficiency, N3 and REM percentages and the arousal index are already
frozen for the whole cohort in numbers/cohorts.json where they exist there. They are recomputed
here from the same source file so that all 10 rows of the panel come from one pass, and the
script asserts agreement with cohorts.json wherever the two overlap, so a drift would fail loudly
rather than be drawn.

Writes numbers/figure1_sleep_profile_v1.json.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../..")))  # repository root, where paths.py lives
import paths
import json
import sys

import numpy as np
import pandas as pd

ROOT = paths.FIGURE_ROOT
sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
from cohort_spec import COHORT_N, apply_cohort

b = apply_cohort(pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet"))
assert len(b) == COHORT_N

# (key, printed label, unit, source column, transform, decimals, family)
SPEC = [
    ("tst_h", "Total sleep time", "hours", "TST_min", lambda s: s / 60.0, 2, "Sleep quantity and structure"),
    ("sleep_efficiency_pct", "Sleep efficiency", "% of time in bed", "sleep_efficiency_pct", None, 1, "Sleep quantity and structure"),
    ("n3_pct", "Deep sleep (N3)", "% of sleep", "N3_pct", None, 1, "Sleep quantity and structure"),
    ("rem_pct", "REM sleep", "% of sleep", "REM_pct", None, 1, "Sleep quantity and structure"),
    ("arousal_index", "Arousal index", "per hour", "arousal_index", None, 1, "Fragmentation and breathing"),
    ("ahi", "Apnea-hypopnea index", "per hour", "AHI", None, 1, "Fragmentation and breathing"),
    ("spo2_mean", "Mean oxygen saturation", "%", "spo2_mean", None, 1, "Oxygenation"),
    ("spo2_nadir", "Lowest oxygen saturation", "%", "spo2_nadir_corrected", None, 1, "Oxygenation"),
    ("t90", "Time with oxygen <90%", "% of night", "spo2_pct_below_90", None, 2, "Oxygenation"),
    ("t88", "Time with oxygen <88%", "% of night", "spo2_pct_below_88", None, 2, "Oxygenation"),
]

rows = {}
for key, label, unit, col, fn, dec, family in SPEC:
    s = b[col].astype(float)
    s = fn(s) if fn else s
    s = s.dropna()
    rows[key] = {
        "label": label, "unit": unit, "family": family, "source_column": col,
        "n": int(len(s)),
        "median": round(float(s.median()), dec),
        "q1": round(float(s.quantile(0.25)), dec),
        "q3": round(float(s.quantile(0.75)), dec),
        "mean": round(float(s.mean()), dec),
        "sd": round(float(s.std()), dec),
        "p05": round(float(s.quantile(0.05)), dec),
        "p95": round(float(s.quantile(0.95)), dec),
    }

# agreement with the frozen file, for the five quantities it already carries
C = json.load(open(f"{paths.NUMBERS_DIR}/cohorts.json"))["bdsp"]
CHECKS = [("tst_h", C["tst_hours"], 0.01), ("sleep_efficiency_pct", C["sleep_efficiency_pct"], 0.05),
          ("arousal_index", C["arousal_index"], 0.05), ("t90", C["t90"], 0.005),
          ("ahi", C["ahi"], 0.05)]
checks = []
for key, ref, tol in CHECKS:
    got, want = rows[key]["median"], ref["median"]
    ok = abs(got - want) <= tol
    checks.append({"quantity": key, "recomputed_median": got,
                   "cohorts_json_median": want, "agrees": bool(ok)})
    assert ok, f"{key}: recomputed median {got} against frozen {want}"

out = {
    "what_this_is": "Cohort medians and interquartile ranges for the routine sleep-study "
                    "measures drawn on the revised main Figure 1 panel A. Descriptive only, "
                    "no model is fitted.",
    "built": "2026-08-07",
    "source_file": "data_frozen_v8_2026-09/t90_final.parquet",
    "cohort_rule": "numbers/cohort_spec.apply_cohort: fu_valid == 1, oxygen present, "
                   "oximetry_bad == 0",
    "cohort_n": int(len(b)),
    "measures": rows,
    "agreement_with_cohorts_json": checks,
}
path = f"{paths.NUMBERS_DIR}/figure1_sleep_profile_v1.json"
json.dump(out, open(path, "w"), indent=1)
print(f"cohort {len(b):,}")
for k, v in rows.items():
    print(f"  {v['label']:<28s} {v['median']:>8} ({v['q1']}-{v['q3']})  n={v['n']:,}  [{v['unit']}]")
print(f"\nwritten {path}")
