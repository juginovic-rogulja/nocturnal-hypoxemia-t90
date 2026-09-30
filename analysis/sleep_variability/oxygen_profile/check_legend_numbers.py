"""
Assert every number printed in the legend snippet against the extraction outputs.

The sheet's own build already asserts everything it draws. The legend is prose and would
otherwise be the one place a stale or mistyped number could survive, so each claim in it is
restated here as a check against the source file it came from. Run it after any edit to
_workfiles/legends_updates/eFigureNEW_overnight_oxygen_profile.md.

    python3 check_legend_numbers.py
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import json
import os

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
LEG = (f"{paths.FIGURE_ROOT}/FINAL_FIGURES_2026-08-14/"
       "_workfiles/legends_updates/eFigureNEW_overnight_oxygen_profile.md")

S = json.load(open(f"{HERE}/overnight_profile_summary.json"))
X = json.load(open(f"{HERE}/extraction_summary.json"))
V = json.load(open(f"{HERE}/verify50_summary.json"))
d = pd.read_csv(f"{HERE}/oxyprofile_per_patient.csv", low_memory=False)
ok = d[d.status == "ok"]
text = open(LEG).read()

fails = []


def claim(label, printed, computed, dp=None):
    """`printed` is the string in the legend, `computed` is the value from the source."""
    ok_val = (abs(float(printed.replace(",", "")) - float(computed))
              <= (0.5 * 10 ** -dp if dp is not None else 0.5))
    in_text = printed in text
    print(f"  {'PASS' if ok_val and in_text else 'FAIL'}  {label:<52} legend {printed:>12}"
          f"   source {float(computed):14.6f}   in legend text: {in_text}")
    if not (ok_val and in_text):
        fails.append(label)


ovd = S["deciles_overall"]
st = S["stages_overall"]
le1 = [r["median"] for r in S["deciles_by_band"]["le1"]]
gt10 = [r["median"] for r in S["deciles_by_band"]["gt10"]]

print("cohort and coverage")
claim("cohort n", "19,173", S["cohort_n"])
claim("records with usable oximetry", "19,042", X["records_with_a_profile"])
claim("records with all ten spans", "18,149", X["records_with_all_ten_deciles"])
for key, val in (("le1", "9,917"), ("1to5", "4,315"), ("5to10", "1,654"), ("gt10", "3,287")):
    claim(f"band n {key}", val, S["band_n_cohort"][key])

print("\npanel b counts and medians")
for lab, n in (("Wake", "19,041"), ("N1", "18,551"), ("N2", "18,538"),
               ("N3", "15,201"), ("REM", "17,730")):
    claim(f"{lab} n", n, st[lab]["n"])
claim("recordings with no scored N3", "3,864", int((ok.n3_min == 0).sum()))
for lab, med, q1, q3 in (("Wake", "95.70", "94.33", "96.83"), ("N1", "95.22", "93.71", "96.45"),
                         ("N2", "94.92", "93.35", "96.23"), ("N3", "94.76", "93.05", "96.12"),
                         ("REM", "95.01", "93.11", "96.47")):
    claim(f"{lab} median", med, st[lab]["median"], dp=2)
    claim(f"{lab} q1", q1, st[lab]["q1"], dp=2)
    claim(f"{lab} q3", q3, st[lab]["q3"], dp=2)

print("\npanel a shape")
claim("first span median", "95.4", ovd[0]["median"], dp=1)
claim("second span median", "94.8", ovd[1]["median"], dp=1)
claim("last span median", "95.7", ovd[9]["median"], dp=1)
claim("whole-cohort swing", "0.90", S["headline"]["decile_median_range_pp"], dp=2)
claim("le1 excursion", "0.53", max(le1) - min(le1), dp=2)
claim("le1 first", "96.3", le1[0], dp=1)
claim("le1 lowest", "95.9", min(le1), dp=1)
claim("le1 last", "96.5", le1[-1], dp=1)
claim("gt10 excursion", "2.01", max(gt10) - min(gt10), dp=2)
claim("gt10 first", "92.0", gt10[0], dp=1)
claim("gt10 lowest", "90.7", min(gt10), dp=1)
claim("gt10 last", "92.7", gt10[-1], dp=1)
claim("stage spread", "0.94", S["headline"]["stage_median_spread_pp"], dp=2)
claim("le1 wake minus rem", "0.43", S["headline"]["le1_wake_minus_rem_pp"], dp=2)
claim("gt10 wake minus rem", "1.90", S["headline"]["gt10_wake_minus_rem_pp"], dp=2)

print("\nvalidation quoted in the legend")
claim("mean SpO2 within 0.5", "17,847", X["mean_within_0p5"])
claim("mean SpO2 percent", "93.7", X["mean_within_0p5_pct"], dp=1)
claim("mean SpO2 median abs diff", "0.0004", X["mean_median_abs_diff"], dp=4)
claim("mean SpO2 Spearman", "0.980", X["mean_spearman"], dp=3)
claim("T90 within 0.5", "17,889", X["t90_within_0p5pp"])
claim("T90 percent", "93.9", X["t90_within_0p5pp_pct"], dp=1)
claim("stage minutes reproduced", "95,860", X["stage_minutes_vs_ladder_within_0p01min"])
claim("verification seed", "20260821", V["seed"])
claim("verification patients", "50", V["n_drawn"])
claim("verification values", "1,538", V["n_values_compared"])
claim("verification worst span mean", "0.0016",
      V["by_group"]["decile mean SpO2 (percent)"]["max_abs_diff"], dp=4)

print("\nqualitative claims that must also hold")
checks = [
    ("no value over tolerance in the 50-patient check", V["n_values_over_tolerance"] == 0),
    ("decile count identity exact", X["decile_count_max_abs_dev"] == 0.0),
    ("decile weighted-mean identity under 1e-13",
     X["decile_weighted_mean_max_abs_dev"] < 1e-13),
    ("stage minute agreement under 1e-13", X["stage_minutes_vs_ladder_max_abs_dev_min"] < 1e-13),
    ("stage mean agreement under 1e-8",
     V["by_group"]["stage mean SpO2 (percent)"]["max_abs_diff"] < 1e-8),
    ("REM is the lowest stage above 10 percent",
     min(S["stages_by_band"]["gt10"], key=lambda k: S["stages_by_band"]["gt10"][k]["median"])
     == "REM"),
    ("N3 is the lowest stage in the whole cohort", S["headline"]["stage_lowest"] == "N3"),
    ("the lowest whole-cohort span is the second", S["headline"]["decile_of_lowest_median"] == 2),
    ("no em dash in the legend", "—" not in text),
    ("no semicolon in the legend", ";" not in text),
]
for label, cond in checks:
    print(f"  {'PASS' if cond else 'FAIL'}  {label}")
    if not cond:
        fails.append(label)

print("\n" + "=" * 84)
print("VERDICT: " + ("EVERY LEGEND NUMBER MATCHES ITS SOURCE"
                     if not fails else "FAILURES: " + ", ".join(fails)))
print("=" * 84)
raise SystemExit(1 if fails else 0)
