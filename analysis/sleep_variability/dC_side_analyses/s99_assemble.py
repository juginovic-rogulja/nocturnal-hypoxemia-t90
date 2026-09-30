"""
Collects every headline number the five side analyses produced into one tidy file.

results.csv is long: one row per number, with the analysis it came from, what it is compared
against, and the script that computed it. Nothing here recomputes anything.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import json
import os

import pandas as pd

from common import WORK
import sys as _sys_v8; _sys_v8.path.insert(0, paths.ANALYSIS_DIR)  # v8 sweep 2026-09-12
from cohort_spec import N_MEASURES, N_RANKED_OUTCOMES  # v8 sweep 2026-09-12

HERE = os.path.dirname(os.path.abspath(__file__))
R = []


def add(analysis, key, value, comparator=None, delta=None, note="", script=""):
    R.append({"analysis": analysis, "key": key, "value": value,
              "comparator": comparator, "delta": delta, "note": note, "script": script})


PUB = 0.02033875497249791

# ---------------------------------------------------------------- positive control
pc = json.load(open(os.path.join(WORK, "positive_control.json")))
add("00_control", "published_T90_dC", pc["published_dC"], None, None,
    f"numbers/ranking_v3.csv rank 1 of {N_MEASURES}", "s00_positive_control.py")
add("00_control", "reproduced_T90_dC", pc["reproduced_dC"], pc["published_dC"],
    pc["diff_dC"], "exact", "s00_positive_control.py")
add("00_control", "reproduced_baseline", pc["reproduced_baseline"], pc["published_baseline"],
    pc["reproduced_baseline"] - pc["published_baseline"], "held-out C, age and sex",
    "s00_positive_control.py")
add("00_control", "max_abs_per_outcome_gain_diff", pc["max_abs_per_outcome_gain_diff"],
    0.0, None, f"all {N_RANKED_OUTCOMES} outcomes against ranking_v3_percondition.csv",
    "s00_positive_control.py")

# ---------------------------------------------------------------- A
A = pd.read_csv(os.path.join(WORK, "A_organ_summary.csv"))
for _i, r in A.iterrows():
    add("A_organ", f"dC_{r.organ}", r.dC, PUB, r.dC - PUB,
        f"{int(r.n)} outcomes, {int(r.npos)} positive", "s01_organ_decomposition.py")
ap = pd.read_csv(os.path.join(WORK, "A_organ_percondition.csv"))
add("A_organ", "dC_cvd_composite", float(ap[ap.outcome == "cvd"].gain.iloc[0]), PUB,
    float(ap[ap.outcome == "cvd"].gain.iloc[0]) - PUB,
    "union of hf, ihd, mi, stroke_any, pad", "s01_organ_decomposition.py")
comp = ap[ap.outcome.isin(["hf", "ihd", "mi", "stroke_any", "pad"])]
add("A_organ", "dC_mean_of_cvd_components", float(comp.gain.mean()),
    float(ap[ap.outcome == "cvd"].gain.iloc[0]),
    float(comp.gain.mean()) - float(ap[ap.outcome == "cvd"].gain.iloc[0]),
    "the composite scores below the mean of its parts", "s01_organ_decomposition.py")
add("A_organ", "spread_max_minus_min_organ", float(A.dC.max() - A.dC.min()), None, None,
    f"{A.iloc[0].organ} down to {A.iloc[-1].organ}", "s01_organ_decomposition.py")

# ---------------------------------------------------------------- B
_B_PATH = os.path.join(WORK, "B_result.json")
if os.path.exists(_B_PATH):
    B = json.load(open(_B_PATH))
    for k, v in B["dC"].items():
        add("B_window", f"dC_{k}", v, B["dC"]["t90_wholerec"], v - B["dC"]["t90_wholerec"],
            f"{B['subcohort_n']} recordings, {B['n_outcomes_common']} outcomes, "
            f"baseline {B['baseline']:.6f}", "s02_sleep_vs_whole.py")
    for k, v in B["joint"].items():
        add("B_window", f"dC_{k}", v, B["dC"]["t90_wholerec"], v - B["dC"]["t90_wholerec"],
            "both windows in one model", "s02_sleep_vs_whole.py")
else:
    # v8 (2026-09-13): s02_sleep_vs_whole is a group-P step (pilot shards) and was not rerun; its 2026-08-18 result was parked
    # as a stale cache. The v8 sleep-period T90 sensitivity set (all 19,173 nights) is produced by steps 230/231 instead.
    add("B_window", "NOT_RUN_IN_V8", float("nan"), float("nan"), float("nan"),
        "group-P step s02_sleep_vs_whole not rerun; see the v8 sleep-period sensitivity set (steps 230/231)", "s02_sleep_vs_whole.py")
    print("B block: s02_sleep_vs_whole not run in v8; recorded as NOT_RUN_IN_V8")

# ---------------------------------------------------------------- C
C = json.load(open(os.path.join(WORK, "C_result.json")))
for k in ("baseline", "t90", "t90_plus_nremhr"):
    add("C_combo", f"eTable11_{k}_frozen", C["etable11_frozen"][k], None, None,
        "54 outcomes, penalizer 0.05, numbers/combination_v2.json", "s03_combinations.py")
    add("C_combo", f"eTable11_{k}_refit", C["etable11_refit"][k], C["etable11_frozen"][k],
        C["etable11_refit"][k] - C["etable11_frozen"][k], "confirmed", "s03_combinations.py")
for k, v in C["dC"].items():
    add("C_combo", f"dC_{k}", v, PUB, v - PUB,
        f"ranking pipeline, {N_RANKED_OUTCOMES} outcomes, penalizer 0.01", "s03_combinations.py")
add("C_combo", "spearman_t90_vs_nadir", C["spearman_t90_nadir_z"], None, None,
    "near-duplicates, on the modelled scale", "s03_combinations.py")
add("C_combo", "spearman_t90_vs_nremhr", C["spearman_t90_nremhr_z"], None, None,
    "close to independent", "s03_combinations.py")

# ---------------------------------------------------------------- D
D = json.load(open(os.path.join(WORK, "D_result.json")))
for k, v in D["dC"].items():
    add("D_transform", f"dC_{k}", v, D["published_rint_dC"], v - D["published_rint_dC"],
        f"{D['ncols'][k]} column(s)", "s04_nonlinear.py")

# ---------------------------------------------------------------- E
E = pd.read_csv(os.path.join(WORK, "E_summary.csv"))
for _i, r in E.iterrows():
    add("E_restricted", f"T90_dC_{r.subset}", r.t90_dC, PUB, r.t90_dC - PUB,
        f"{int(r.n_outcomes)} outcomes, rank {int(r.t90_rank)} of {int(r.n_measures)}, "
        f"lead over #2 {r.lead_over_runner_up:+.6f}", "s05_restricted_ranking.py")

# ---------------------------------------------------------------- extra
X = json.load(open(os.path.join(WORK, "EXTRA_agebasis.json")))
add("EXTRA_agebasis", "T90_dC_published_4col", X["t90_dC_published"], None, None,
    "singular cr(df=4) basis, the published spec", "s06_agebasis.py")
add("EXTRA_agebasis", "T90_dC_fullrank_3col", X["t90_dC_fullrank"], X["t90_dC_published"],
    X["t90_rise"], "one basis column dropped", "s06_agebasis.py")
add("EXTRA_agebasis", "mean_rise_across_top9_measures", X["mean_rise_top_measures"],
    X["t90_rise"], None, "the rise is the base model, not the exposure", "s06_agebasis.py")
add("EXTRA_agebasis", "T90_lead_fullrank", X["lead_fullrank"], X["lead_published"],
    X["lead_fullrank"] - X["lead_published"], "rank order of the top 9 unchanged",
    "s06_agebasis.py")

df = pd.DataFrame(R)
for c in ("value", "comparator", "delta"):
    df[c] = pd.to_numeric(df[c], errors="coerce").round(8)
out = os.path.join(HERE, "results.csv")
df.to_csv(out, index=False)
print(f"wrote {out}  ({len(df)} rows)")
print(df.to_string(index=False, max_colwidth=44))
