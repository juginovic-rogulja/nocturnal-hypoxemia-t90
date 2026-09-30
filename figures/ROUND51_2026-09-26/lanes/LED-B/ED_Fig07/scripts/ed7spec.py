#!/usr/bin/env python3
"""Round 49 LED-B copy (+1 pt rule, Fig 4d records repointed to the round-38 Main_Fig4 lane). ED_Fig07 (V14 lane L4, round 37): the ONE spec module every ED_Fig07 script reads its paths and layout constants from.
Copy of ROUND30_2026-09-08/V7_L4_APNEA/ED_Fig07/scripts/ed7spec.py (beside as ed7spec_PRE_V8_1.py), repointed: lane V14_L4_APNEA,
base sheet = the V13 ED_Fig07 (LF4 layout: panel a at y 8, "HR (95% CI)" column title, panel b in Fig 4d's design, page 900.628),
numbers gated by the v8.1 sidecar, the OLD side = the V13 sheet's own text layer (no v7 snapshot). Panel a geometry = the round-30
measured record of the base forest (unchanged on V13, which places the round-30 panel a byte for byte 8 pt down the page)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import datetime, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import l4lib as L

SHEET = "ED_Fig07"
SD = f"{L.LANE}/{SHEET}"; WORK = f"{SD}/work"; VER = f"{SD}/verify"; CROPS = f"{VER}/crops"; SCRIPTS = f"{L.LANE}/scripts"
BASE_PDF = f"{L.BASE}/{SHEET}.pdf"                       # the V13 sheet (design baseline), never edited in place
BASE_TEXT = f"{WORK}/base_text.json"                      # its probed text layer (00_probe_v13.py)
SRC_A = f"{L.NUM}/severe_apnea_negctrl_v1.json"           # panel a, step 126 (v8.1)
SRC_A_XCHECK = f"{L.NUM}/extras2_v2.json"                 # the same fits, cross-check only (step 115)
SRC_B = f"{L.NUM}/alenfig3_cross_v1.json"                 # panel b, key cross9, step 123 (v8.1)
R30_GEOM_A = f"{L.R30}/V7_L4_APNEA/ED_Fig07/work/panel_a_geometry.json"   # the base forest's measured geometry (round 30), still the V13 panel a
LF4_B_DRAWN = f"{paths.FIGURE_ROOT}/ROUND31_2026-09-10/LF4_FIG4_ED7/ED_Fig07/verify/ED_Fig07_b_drawn.json"   # the V13 panel b strings and layout
LF4_LAYOUT = f"{paths.FIGURE_ROOT}/ROUND31_2026-09-10/LF4_FIG4_ED7/ED_Fig07/work/sheet_layout.json"
_R38 = f"{paths.FIGURE_ROOT}/ROUND38_2026-09-16/figures/lanes"   # round 49: the round-38 Main_Fig4 lane records (the V13 geometry of Fig 4d), read-only
FIG4_D_GEOM = f"{_R38}/Main_Fig4/work/panel_d_geometry.json"; FIG4_D_DRAWN = f"{_R38}/Main_Fig4/verify/Main_Fig4_d_drawn.json"
OUT_PDF = f"{SD}/{SHEET}.pdf"; OUT_PNG = f"{SD}/{SHEET}_150dpi.png"
PANEL_A_RAW = f"{WORK}/panel_a_raw.pdf"; PANEL_A = f"{WORK}/panel_a.pdf"; GEOM_A = f"{WORK}/panel_a_geometry.json"
PANEL_B_RAW = f"{WORK}/panel_b_raw.pdf"; PANEL_B = f"{WORK}/panel_b.pdf"; GEOM_B = f"{WORK}/panel_b_geometry.json"
DRAWN_A = f"{VER}/{SHEET}_a_drawn.json"; DRAWN_B = f"{VER}/{SHEET}_b_drawn.json"
LAYOUT = f"{WORK}/sheet_layout.json"; CHECKS = f"{VER}/checks.txt"; CHANGES = f"{VER}/CHANGES_{SHEET}.csv"

# ---------------------------------------------------------------- sheet layout (page points, y down), the V13 (LF4) layout
EDGE = 4.0 / 25.4 * 72.0; LETTER_X = 12.96; LETTER_SIZE = 13.0 + L.PT_PLUS; LETTER_BASELINE = 21.0   # round 49: letters 14 pt
SHIFT_A = 8.0; PANEL_B_OVERLAP = 8.0
HDR_A = "HR (95% CI)"; HDR_A_FS = 9.5 + L.PT_PLUS         # the round-31 column title on panel a, right-aligned on the estimate column (round 49: 10.5 pt)
GROUPS = ["No or mild", "Moderate", "Severe"]; LEVELS = ["<=1%", ">1-10%", ">10%"]; REF = "No or mild|<=1%"
GROUP_LAB = {"No or mild": ["No or mild", "(AHI <15)"], "Moderate": ["Moderate", "(AHI 15 to <30)"], "Severe": ["Severe", "(AHI 30 or more)"]}
LEVEL_LAB = {"<=1%": "≤1%", ">1-10%": ">1–10%", ">10%": ">10%"}
COLS = ["Cardiovascular composite", "Type 2 diabetes", "Acute kidney injury", "Heart failure"]
COL_LAB = {"Cardiovascular composite": ["Cardiovascular", "composite"], "Type 2 diabetes": ["Type 2", "diabetes"], "Acute kidney injury": ["Acute kidney", "injury"], "Heart failure": ["Heart", "failure"]}
LVL_A = ["Sleep apnea", "severity"]; LVL_O = "T90"
REF_KEY = "Reference, no or mild apnea with T90 ≤1%"; STAR_KEY_B = "* q < 0.05, ** q < 0.01, *** q < 0.001"; RAMP_TITLE = "Hazard ratio"; TICKS_B = ["1", "1.5", "2", "3"]
KEY_STARS_A = "* P < 0.05, ** P < 0.01, *** P < 0.001"; BLOCK_HEAD = "Negative controls"; XLAB_A = "Rate of new diagnosis (95% CI)"


def status(pct, text):
    line = f"- {datetime.datetime.now().strftime('%H:%M')} ET, ED_Fig07 {pct} percent: {text}"
    open(f"{L.LANE}/STATUS.md", "a").write(line + "\n"); print(line, flush=True)


def cluster(values, tol):
    out = []
    for v in sorted(values):
        if out and abs(v - out[-1][-1]) <= tol: out[-1].append(v)
        else: out.append([v])
    return [sum(c) / len(c) for c in out]
