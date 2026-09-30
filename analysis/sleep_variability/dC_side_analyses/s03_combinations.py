"""
C. T90 with a second measurement beside it.

Two jobs.

  1. CONFIRM eTable 11. Its source is numbers/combination_v2.json, written by
     numbers/run_combo_v2.py. Its frozen values are 0.01839 for the best single oxygen
     measurement and 0.02428 for that measurement plus non-REM heart rate, over 54 outcomes at
     penalizer 0.05 against a baseline of 0.6368. Those are NOT the ranking's numbers: the
     ranking scores 48 outcomes at penalizer 0.01. Both configurations are refitted here.

  2. RECOMPUTE CLEANLY in the ranking pipeline, and add T90 plus the corrected nadir, the pair
     that approximates how deep the desaturation went against how long it lasted.

Every measure enters as the within-site rank inverse normal, exactly as in the ranking.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import json
import os

import numpy as np
import pandas as pd

from common import (PEN, T90ROOT, WORK, build_frame, combo_outcomes, heldout,
                    ranking_outcomes, site_rint, summarise)
import sys as _sys_v8; _sys_v8.path.insert(0, paths.ANALYSIS_DIR)  # v8 sweep 2026-09-12
from cohort_spec import N_RANKED_OUTCOMES  # v8 sweep 2026-09-12
__PUB_DC_V7__ = float(__import__("pandas").read_csv(f"{paths.NUMBERS_DIR}/ranking_v3.csv", comment="#").set_index("feature").loc["spo2_pct_below_90", "dC"])   # v7: T90 gain from the CURRENT ranking, was the August literal 0.02033875497249791


T90 = "spo2_pct_below_90"
HR = "nrem_hr_bpm"
NADIR = "spo2_nadir_corrected"
COLS = [T90, HR, NADIR]

print("=" * 78)
print("C. T90 PLUS A SECOND MEASUREMENT")
print("=" * 78)

d = build_frame(extra_master_cols=COLS, cache="frame_combo.pkl")
for c in COLS:
    cov = float(d[c].notna().mean())
    print(f"   {c:<24}complete in {100 * cov:.1f}% of the cohort")
    assert cov >= 0.60, f"{c} is below the pipeline's 60% completeness bar"
    d[c + "__z"] = site_rint(d, c)

Z = lambda f: f + "__z"
OUT48 = ranking_outcomes(d)
OUT54 = combo_outcomes(d)
assert len(OUT48) == N_RANKED_OUTCOMES, (len(OUT48), len(OUT54))  # v8 sweep: OUT54 is every outcome over the 150-event bar, printed, not asserted

# ------------------------------------------------------------------ 1. confirm eTable 11
print("\n" + "-" * 78)
print("1. eTable 11 as published: 54 outcomes, penalizer 0.05")
print("-" * 78)
frozen = json.load(open(f"{paths.NUMBERS_DIR}/combination_v2.json"))
specs11 = {"__base__": [], "t90": [Z(T90)], "t90_plus_nremhr": [Z(T90), Z(HR)]}
r11 = heldout(specs11, OUT54, d, penalizer=0.05, pkl_name="C_e11.pkl")
s11 = summarise(r11)
b11 = float(r11[r11.spec == "__base__"].c.mean())
got = {"baseline": b11,
       "t90": float(s11[s11.spec == "t90"].dC.iloc[0]),
       "t90_plus_nremhr": float(s11[s11.spec == "t90_plus_nremhr"].dC.iloc[0])}
exp = {"baseline": frozen["baseline"],
       "t90": frozen["single"][T90],
       "t90_plus_nremhr": frozen["combinations"]["oxygen_plus_" + HR]}
print(f"{'quantity':<22}{'frozen':>10}{'refit here':>13}{'diff':>12}")
for k in ("baseline", "t90", "t90_plus_nremhr"):
    print(f"{k:<22}{exp[k]:>10.5f}{got[k]:>13.5f}{got[k] - exp[k]:>+12.2e}")
ok11 = all(abs(got[k] - exp[k]) < 5e-5 for k in got)
print(f"eTable 11 CONFIRMED: {ok11}   "
      f"(the published 0.018 -> 0.024 is {exp['t90']:.5f} -> {exp['t90_plus_nremhr']:.5f})")

# the same two configurations on the ranking's 48 at 0.05, which run_combo_v2 also froze
sens = frozen["sensitivity_ranking_outcome_set"]
r11b = heldout(specs11, OUT48, d, penalizer=0.05, pkl_name="C_e11b.pkl")
s11b = summarise(r11b)
print(f"\nrun_combo_v2's own {len(OUT54)}-outcome sensitivity block, refit:")
print(f"   best single oxygen   frozen {sens['best_single_oxygen']:.5f}   "
      f"here {float(s11b[s11b.spec == 't90'].dC.iloc[0]):.5f}")
print(f"   oxygen + non-REM HR  frozen {sens['oxygen_plus_' + HR]:.5f}   "
      f"here {float(s11b[s11b.spec == 't90_plus_nremhr'].dC.iloc[0]):.5f}")

# ------------------------------------------------------------------ 2. the ranking pipeline
print("\n" + "-" * 78)
print(f"2. Recomputed in the ranking pipeline: {N_RANKED_OUTCOMES} outcomes, penalizer {PEN}")
print("-" * 78)
specs = {
    "__base__": [],
    "t90_alone": [Z(T90)],
    "nremhr_alone": [Z(HR)],
    "nadir_alone": [Z(NADIR)],
    "t90_plus_nremhr": [Z(T90), Z(HR)],
    "t90_plus_nadir": [Z(T90), Z(NADIR)],
    "t90_plus_nadir_plus_nremhr": [Z(T90), Z(NADIR), Z(HR)],
}
res = heldout(specs, OUT48, d, penalizer=PEN, pkl_name="C_rank.pkl")
s = summarise(res)
base = float(res[res.spec == "__base__"].c.mean())
print(f"held-out baseline {base:.6f}   (ranking_v3.csv: 0.647540)\n")
print(f"{'model':<30}{'dC':>12}{'n>0 of ' + str(N_RANKED_OUTCOMES):>11}{'over T90 alone':>16}")
print("-" * 70)
t90_dc = float(s[s.spec == "t90_alone"].dC.iloc[0])
for _i, r in s.iterrows():
    print(f"{r.spec:<30}{r.dC:>+12.6f}{int(r.npos):>11}{r.dC - t90_dc:>+16.6f}")

pub = __PUB_DC_V7__
assert abs(t90_dc - pub) < 1e-9, (t90_dc, pub)
print(f"\nT90 alone reproduces the published {pub:.6f} exactly, so the rows above are "
      f"comparable to it.")
g = s.set_index("spec").dC
print(f"\nT90 + non-REM heart rate: {g['t90_plus_nremhr']:+.6f}, "
      f"{g['t90_plus_nremhr'] - t90_dc:+.6f} over T90 alone "
      f"({100 * (g['t90_plus_nremhr'] / t90_dc - 1):+.1f}%)")
print(f"T90 + corrected nadir:    {g['t90_plus_nadir']:+.6f}, "
      f"{g['t90_plus_nadir'] - t90_dc:+.6f} over T90 alone "
      f"({100 * (g['t90_plus_nadir'] / t90_dc - 1):+.1f}%)")
print(f"depth and duration are near-duplicates: Spearman between T90 and the nadir on the "
      f"modelled scale {float(d[[Z(T90), Z(NADIR)]].corr(method='spearman').iloc[0, 1]):.3f}")
print(f"T90 and non-REM heart rate: "
      f"{float(d[[Z(T90), Z(HR)]].corr(method='spearman').iloc[0, 1]):.3f}")

s.insert(0, "analysis", "C")
s.to_csv(os.path.join(WORK, "C_summary.csv"), index=False)
res.to_csv(os.path.join(WORK, "C_percondition.csv"), index=False)
json.dump({"etable11_confirmed": bool(ok11), "etable11_frozen": exp, "etable11_refit": got,
           "ranking_pipeline_baseline": base,
           "dC": {k: float(v) for k, v in g.items()},
           "spearman_t90_nadir_z": float(d[[Z(T90), Z(NADIR)]].corr(method="spearman").iloc[0, 1]),
           "spearman_t90_nremhr_z": float(d[[Z(T90), Z(HR)]].corr(method="spearman").iloc[0, 1])},
          open(os.path.join(WORK, "C_result.json"), "w"), indent=2)
print("\nwritten -> _work/C_summary.csv, _work/C_percondition.csv, _work/C_result.json")
