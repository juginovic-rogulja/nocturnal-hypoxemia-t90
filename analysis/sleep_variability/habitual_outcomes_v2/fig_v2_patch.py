#!/usr/bin/env python3
"""
fig_v2_patch.py -- repoint the copied make_fig26_habitual_x_oxygen.py at the v2 outputs and
replace every v1-pinned identity assert and every rendered v1 number with its v2 value.
Each replacement must occur exactly once or the patch refuses to apply. The (old, new) pairs
below are the complete, machine-readable change log for the figure script.

v2 values were computed from habitual_outcomes_v2/{cells_crossed.csv, results_crossed.csv,
summary_crossed.json, attack_power_cells.csv, attack_headtohead.csv} on 2026-08-13.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
F = (f"{paths.SV_ROOT}/"
     "habitual_outcomes_v2/make_fig26_habitual_x_oxygen.py")

PAIRS = [
    # ---------------------------------------------------------------- source and output paths
    ('SRC = Path(f"{paths.SV_ROOT}/"\n'
     '           "habitual_outcomes")',
     'SRC = Path(f"{paths.SV_ROOT}/"\n'
     '           "habitual_outcomes_v2")'),
    ('RAS = WORK', 'RAS = str(SRC)   # v2: rasters filed with the v2 outputs, not the workfiles'),
    ('name = "eFigure16_habitual_sleep_x_oxygen"', 'name = "eFigure16_habitual_sleep_x_oxygen_v2"'),
    ('pdf_path = f"{_supp_dir()}/{name}.pdf"',
     'pdf_path = f"{SRC}/{name}.pdf"   # v2: never into New_Figures, the main session gates that'),
    ('json.dump(DRAWN, open(WORK / "Figure26_drawn_values.json", "w"), indent=1)',
     'json.dump(DRAWN, open(SRC / "Figure26_drawn_values_v2.json", "w"), indent=1)'),
    ('print("wrote", WORK / "Figure26_drawn_values.json")',
     'print("wrote", SRC / "Figure26_drawn_values_v2.json")'),
    # ------------------------------------------------------------------- identity assert pins
    ('assert S["cohort_n"] == 19173 and S["habitual_n"] == 8711',
     'assert S["cohort_n"] == 19173 and S["habitual_n"] == 5295   # v2 (v1: 8711)'),
    ('assert HEAD_TO_HEAD_N == S["habitual_n"] == 8711, HEAD_TO_HEAD_N',
     'assert HEAD_TO_HEAD_N == S["habitual_n"] == 5295, HEAD_TO_HEAD_N   # v2 (v1: 8711)'),
    ('assert S["reference_cell"] == "6-<7h(ref) | normal T90<=1%" and S["reference_cell_n"] == 644',
     'assert S["reference_cell"] == "6-<7h(ref) | normal T90<=1%" and S["reference_cell_n"] == 402'),
    ('assert len(A) == 12 and int(A.n.sum()) == 8711',
     'assert len(A) == 12 and int(A.n.sum()) == 5295'),
    ('assert gridA["<5h | normal T90<=1%"]["n"] == 1292',
     'assert gridA["<5h | normal T90<=1%"]["n"] == 662'),
    ('assert list(gridA["<5h | normal T90<=1%"]["events"].values()) == [114, 46, 119, 99, 69]',
     'assert list(gridA["<5h | normal T90<=1%"]["events"].values()) == [65, 19, 71, 57, 41]'),
    ('assert gridA["<5h | low T90>10%"]["n"] == 497',
     'assert gridA["<5h | low T90>10%"]["n"] == 257'),
    ('assert list(gridA["<5h | low T90>10%"]["events"].values()) == [67, 40, 47, 52, 49]',
     'assert list(gridA["<5h | low T90>10%"]["events"].values()) == [40, 18, 23, 34, 30]'),
    ('assert gridA["5-<6h | normal T90<=1%"]["n"] == 460 and gridA["5-<6h | low T90>10%"]["n"] == 147',
     'assert gridA["5-<6h | normal T90<=1%"]["n"] == 301 and gridA["5-<6h | low T90>10%"]["n"] == 87'),
    ('assert gridA[REF_CELL]["n"] == 644 and gridA[WEAK_CELL]["n"] == 163',
     'assert gridA[REF_CELL]["n"] == 402 and gridA[WEAK_CELL]["n"] == 96'),
    ('assert gridA[">=7h | normal T90<=1%"]["n"] == 2273',
     'assert gridA[">=7h | normal T90<=1%"]["n"] == 1622'),
    ('assert list(gridA[">=7h | normal T90<=1%"]["events"].values()) == [192, 113, 199, 102, 92]',
     'assert list(gridA[">=7h | normal T90<=1%"]["events"].values()) == [121, 86, 145, 68, 61]'),
    ('assert gridA[">=7h | low T90>10%"]["n"] == 598',
     'assert gridA[">=7h | low T90>10%"]["n"] == 329'),
    ('assert list(gridA[">=7h | low T90>10%"]["events"].values()) == [94, 46, 56, 70, 84]',
     'assert list(gridA[">=7h | low T90>10%"]["events"].values()) == [49, 32, 32, 35, 39]'),
    ('assert weak["outcomes_under_5_events"] == 8 and weak["significant_negative_controls"] == 2',
     'assert weak["outcomes_under_5_events"] == 10 and weak["significant_negative_controls"] == 1'),
    ('assert _cvd == [1.25, 2.14, 1.16, 2.30], _cvd',
     'assert _cvd == [1.41, 2.48, 1.11, 2.52], _cvd'),
    ('assert [round(r["hr"], 2) for r in panelB["Heart failure"]] == [1.36, 2.76, 1.02, 3.50]',
     'assert [round(r["hr"], 2) for r in panelB["Heart failure"]] == [1.71, 3.79, 1.13, 3.79]'),
    ('assert [round(r["hr"], 2) for r in panelB["Death from any cause"]] == [0.84, 2.02, 1.22, 1.99]',
     'assert [round(r["hr"], 2) for r in panelB["Death from any cause"]] == [0.65, 1.64, 1.17, 2.24]'),
    ('assert [round(r["hr"], 2) for r in panelB["Hypertension"]] == [1.28, 2.19, 0.93, 1.86]',
     'assert [round(r["hr"], 2) for r in panelB["Hypertension"]] == [1.36, 1.83, 0.89, 1.64]'),
    ('assert [round(r["hr"], 2) for r in panelB["Type 2 diabetes"]] == [1.57, 2.52, 0.87, 2.36]',
     'assert [round(r["hr"], 2) for r in panelB["Type 2 diabetes"]] == [1.78, 3.11, 0.89, 2.33]'),
    ('assert [round(r["hr"], 2) for r in panelB["Dementia"]] == [1.11, 0.89, 1.24, 1.21]',
     'assert [round(r["hr"], 2) for r in panelB["Dementia"]] == [1.13, 0.66, 1.41, 1.70]'),
    ('assert [bandsum[b]["significant"] for b in HAB] == [23, 15, 18, 27]',
     'assert [bandsum[b]["significant"] for b in HAB] == [18, 9, 15, 26]'),
    ('assert [bandsum[b]["fitted"] for b in HAB] == [44, 33, 37, 45]',
     'assert [bandsum[b]["fitted"] for b in HAB] == [41, 26, 32, 44]'),
    ('assert [bandsum[b]["median_hr"] for b in HAB] == [1.67, 1.50, 1.82, 1.70]',
     'assert [bandsum[b]["median_hr"] for b in HAB] == [1.55, 1.71, 2.10, 1.80]'),
    ('assert mirror["normal T90<=1%"] == {"fitted": 132, "significant": 12, "median_hr": 1.11}',
     'assert mirror["normal T90<=1%"] == {"fitted": 132, "significant": 6, "median_hr": 1.08}'),
    ('assert mirror["intermediate 1-10%"]["significant"] == 6',
     'assert mirror["intermediate 1-10%"]["significant"] == 7'),
    ('assert mirror["low T90>10%"]["significant"] == 9',
     'assert mirror["low T90>10%"]["significant"] == 7'),
    ('assert DRAWN["caveats"]["across_all_outcomes"]["low_oxygen_cells_raised"] == 85',
     'assert DRAWN["caveats"]["across_all_outcomes"]["low_oxygen_cells_raised"] == 71'),
    ('assert DRAWN["caveats"]["across_all_outcomes"]["normal_oxygen_median_hr"] == 1.12',
     'assert DRAWN["caveats"]["across_all_outcomes"]["normal_oxygen_median_hr"] == 1.09'),
    ('assert panelD["z_hab"]["significant"] == 12 and panelD["z_hab"]["surviving_bh"] == 4',
     'assert panelD["z_hab"]["significant"] == 10 and panelD["z_hab"]["surviving_bh"] == 5'),
    ('assert panelD["z_tst"]["significant"] == 6 and panelD["z_tst"]["surviving_bh"] == 0',
     'assert panelD["z_tst"]["significant"] == 1 and panelD["z_tst"]["surviving_bh"] == 0'),
    ('assert panelD["z_t90"]["significant"] == 32 and panelD["z_t90"]["surviving_bh"] == 31',
     'assert panelD["z_t90"]["significant"] == 33 and panelD["z_t90"]["surviving_bh"] == 32'),
    ('assert panelD["z_hab"]["control_lo"] == 0.893 and panelD["z_hab"]["control_hi"] == 0.952',
     'assert panelD["z_hab"]["control_lo"] == 0.864 and panelD["z_hab"]["control_hi"] == 0.930'),
    ('assert DRAWN["caveats"]["confounding_band"]["outcomes_outside_band"] == 6',
     'assert DRAWN["caveats"]["confounding_band"]["outcomes_outside_band"] == 4'),
    # ------------------------------------------------------------------------ rendered prose
    ('FT(0.80, 21.52, "Where the 8,711 patients sit when habitual home sleep is crossed with "',
     'FT(0.80, 21.52, "Where the 5,295 patients sit when habitual home sleep is crossed with "'),
    ('"Reference cell 6 to under 7 h with T90 at or below 1%, N = 644, outlined in black. The 6 to "',
     '"Reference cell 6 to under 7 h with T90 at or below 1%, N = 402, outlined in black. The 6 to "'),
    ('axB.set_xlabel("Hazard ratio against 6 to under 7 h with normal oxygen (N = 644)",',
     'axB.set_xlabel("Hazard ratio against 6 to under 7 h with normal oxygen (N = 402)",'),
    ('"(N = 1,292 at under 5 h and 2,273 at 7 h or more)"),',
     '"(N = 662 at under 5 h and 1,622 at 7 h or more)"),'),
    ('"(N = 497 at under 5 h and 598 at 7 h or more)")],',
     '"(N = 257 at under 5 h and 329 at 7 h or more)")],'),
    ('"band   (N = 497, 147, 163 and 598 across the four bands)"),',
     '"band   (N = 257, 87, 96 and 329 across the four bands)"),'),
    ('"(N = 1,292, 460, 644 and 2,273)")],',
     '"(N = 662, 301, 402 and 1,622)")],'),
    ('DRAWN["panelB"] = {"reference": f"{REF_CELL}, N=644",',
     'DRAWN["panelB"] = {"reference": f"{REF_CELL}, N=402",'),
    # -------------------------------------------------------------- module docstring header
    ('home sleep is what the patient says they sleep, pulled out of clinical notes for 8,711 of the',
     'home sleep is what the patient says they sleep, pulled out of clinical notes for 5,295 of the'),
    ('     8,711 patients, with the negative-control band drawn',
     '     5,295 patients, with the negative-control band drawn'),
]

s = open(F).read()
for old, new in PAIRS:
    n = s.count(old)
    assert n == 1, f"pattern occurs {n} times, refusing:\n{old[:100]}"
    s = s.replace(old, new)
open(F, "w").write(s)
print(f"applied {len(PAIRS)} replacements to {F}")
