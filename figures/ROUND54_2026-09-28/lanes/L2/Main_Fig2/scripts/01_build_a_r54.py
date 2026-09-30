#!$T90_PY
"""ROUND 54 (2026-09-29), lane L2: panel a of Main_Fig2 at the round-54 sizes (row labels, headers, key, ticks, HR and q values 14 pt,
column headers and the caption 15 pt). Copy of ROUND49 .../L2/Main_Fig2/scripts/01_build_a_r49.py (byte copy beside as
01_build_a_r49_R49_ORIGINAL.py): the DATA section (rows, blocks, BH q, the step-168 full-precision checkpoint guard, the printed
strings) is unchanged; the GEOMETRY section is new. The V13 left header band cannot sit beside 14 pt labels in any layout that
fits the sheet (see ../../PLAN.md), so panel a is ONE column with each organ header as a bold row above its block, the two-line
V13 headers set on one line (declared in the census), a 13.9 pt row pitch (65 lines: 58 rows and 7 headers), the x axis with its
V30 range and tick set, the key at the top and the caption below. This script is also the PLANNER: it writes work/layout.json with
the page height, the a/b split and the b/c split that 01_build_b_r54.py, 01_build_c_r54.py and 02_compose_r54.py read.
No data value, colour, marker size, line width or printed string changes. usage: 01_build_a_r54.py"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths

import json, math, os, sys
import numpy as np, pandas as pd
from fig2_common import *
import v14lib

OLDMODE = False   # the v7 positive control is not run in round 54
SRC = OLD if OLDMODE else NUM
TAG = "_OLD" if OLDMODE else ""
BD_PATH, NCP_PATH, COH_PATH = f"{SRC}/bdsp_diseases_v3.csv", f"{SRC}/negcontrols_v4_panel.csv", f"{SRC}/cohorts.json"
for p in (BD_PATH, NCP_PATH, COH_PATH): hydrated(p)
prov = {} if OLDMODE else {os.path.basename(p): sidecar(p) for p in (BD_PATH, NCP_PATH, COH_PATH)}

BD = pd.read_csv(BD_PATH); NCP = pd.read_csv(NCP_PATH); COH = json.load(open(COH_PATH))
N_COHORT = int(COH["bdsp"]["n_analysis"])
n_ctrl = int(BD.negative_control.sum())
assert n_ctrl == len(NCP) == 5, (n_ctrl, len(NCP))
assert set(NCP.condition) == set(BD[BD.negative_control].disease), "control panel is not the control rows of bdsp_diseases_v3"
BDI = BD.set_index("disease")
for r in NCP.itertuples():
    a = BDI.loc[r.condition]
    assert (round(a.t90_hr, 3), round(a.t90_lo, 3), round(a.t90_hi, 3)) == (round(r.hr, 3), round(r.lo, 3), round(r.hi, 3)), (r.condition, "panel csv drifted from bdsp_diseases_v3")

# ------------------------------------------------------------------ the sheet's block design, read from the V13 text layer (not typed)
G = load_geometry(); T = load_text(); A = [r for r in G["records"] if r["region"] == "a"]; SA = spans_in(T, "a")
def w_(b): return b[2] - b[0]
def h_(b): return b[3] - b[1]
HDR_LINES = {"Cardiovascular": ["Cardiovascular"], "Respiratory": ["Respiratory"], "Metabolic, kidney and liver": ["Metabolic,", "kidney and liver"],
             "Infection, blood and other": ["Infection,", "blood and other"], "Neurological and psychiatric": ["Neurological", "and psychiatric"],
             "Musculoskeletal and sensory": ["Musculoskeletal", "and sensory"], "Negative controls": ["Negative", "controls"]}
BLOCK_ORDER = list(HDR_LINES)
hdr_sp = [s for s in SA if s["font"] == "Arial-BoldMT" and abs(s["size"] - 10.0) < 0.01 and abs(s["origin"][0] - 14.17) < 0.5 and s["origin"][1] > 20.0]
hdr_lines_base = {s["text"]: s["origin"] for s in hdr_sp}
assert set(hdr_lines_base) == set(l for v in HDR_LINES.values() for l in v), set(hdr_lines_base) ^ set(l for v in HDR_LINES.values() for l in v)
labels = sorted([s for s in SA if s["font"] == "ArialMT" and abs(s["size"] - 10.0) < 0.01 and abs(s["origin"][0] - 97.15) < 0.5], key=lambda s: s["origin"][1])
base_lab = [s["text"] for s in labels]
# V13 block membership: a label belongs to the block whose first header line is the nearest one above or level with it
first_line_y = {blk: hdr_lines_base[lines[0]][1] for blk, lines in HDR_LINES.items()}
def block_of_base(y):
    cands = [(blk, fy) for blk, fy in first_line_y.items() if fy <= y + 1.0]
    return max(cands, key=lambda t: t[1])[0]
V13_BLOCK = {s["text"]: block_of_base(s["origin"][1]) for s in labels}
V13_SIZES = {blk: sum(1 for b in V13_BLOCK.values() if b == blk) for blk in BLOCK_ORDER}
assert sum(V13_SIZES.values()) == len(labels), V13_SIZES
# organ map for conditions the V13 sheet did not carry
sys.path.insert(0, paths.ANALYSIS_DIR)
import disease_definitions as DD
ORGAN_TO_BLOCK = {"Cardiac": "Cardiovascular", "Composite": "Cardiovascular", "Respiratory": "Respiratory", "Metabolic": "Metabolic, kidney and liver",
                  "Kidney": "Metabolic, kidney and liver", "Liver": "Metabolic, kidney and liver", "Infection": "Infection, blood and other",
                  "Other": "Infection, blood and other", "Mortality": "Infection, blood and other", "Neuro/psych": "Neurological and psychiatric",
                  "Sensory/MSK": "Musculoskeletal and sensory"}
KEY2NAME = {r.key: r.disease for r in BD.itertuples()}
def block_of(disease, key, is_ctrl):
    if is_ctrl: return "Negative controls"
    if disease in V13_BLOCK and V13_BLOCK[disease] != "Negative controls": return V13_BLOCK[disease]
    og = DD.ORGAN_GROUP.get(key); assert og in ORGAN_TO_BLOCK, (disease, key, og, "no organ group for a condition new to the sheet")
    return ORGAN_TO_BLOCK[og]

outc = BD[~BD.negative_control].copy(); outc["q"] = bh(outc.t90_p.values)
ctrl = BD[BD.negative_control].copy(); ctrl["q"] = bh(ctrl.t90_p.values)
Q = {**outc.set_index("disease")["q"].to_dict(), **ctrl.set_index("disease")["q"].to_dict()}
N_OUT = len(outc); N_SIG = int((outc.q < 0.05).sum()); N_RAISED = int(((outc.q < 0.05) & (outc.t90_hr > 1)).sum())
BELOW = sorted(outc[(outc.q < 0.05) & (outc.t90_hr < 1)].disease)
CTRL_SIG = sorted(ctrl[ctrl.q < 0.05].disease)
FLOOR = float(NCP.hr.max()); FLOOR_BY = NCP.loc[NCP.hr.idxmax(), "condition"]
assert abs(FLOOR - float(ctrl.t90_hr.max())) < 5e-4
KEY_FLOOR = f"Level the negative controls reach, {halfup(FLOOR, 2)}"
if CTRL_SIG: print(f"NOTE: negative control(s) FDR-significant within the control family: {CTRL_SIG} (q {[round(Q[c], 4) for c in CTRL_SIG]}). Printed as the data say, flagged for the coordinator.")

blocks_new = {blk: [] for blk in BLOCK_ORDER}
for r in BD.itertuples():
    blocks_new[block_of(r.disease, r.key, bool(r.negative_control))].append(r.disease)
NEW_TO_SHEET = sorted(set(BD.disease) - set(base_lab)); LEFT_SHEET = sorted(set(base_lab) - set(BD.disease))
print(f"conditions new to the sheet: {NEW_TO_SHEET} -> blocks {[block_of(c, BDI.loc[c].key if 'key' in BDI.columns else None, bool(BDI.loc[c].negative_control)) for c in NEW_TO_SHEET]}; left the sheet: {LEFT_SHEET}")
rows = []      # (block, condition, hr, lo, hi, class, q)
NCPI = NCP.set_index("condition")
for blk in BLOCK_ORDER:
    conds = blocks_new[blk]
    if blk == "Negative controls":
        for c in sorted(conds, key=lambda c: -float(NCPI.loc[c].hr)):
            n = NCPI.loc[c]; rows.append((blk, c, float(n.hr), float(n.lo), float(n.hi), "ctrl", float(Q[c])))
    else:
        for c in sorted(conds, key=lambda c: -BDI.loc[c, "t90_hr"]):
            r = BDI.loc[c]; rows.append((blk, c, float(r.t90_hr), float(r.t90_lo), float(r.t90_hi), "sig" if Q[c] < 0.05 else "ns", float(Q[c])))
BLOCK_SIZES = [len(blocks_new[b]) for b in BLOCK_ORDER]
assert len(rows) == len(BD) and all(n > 0 for n in BLOCK_SIZES), (len(rows), BLOCK_SIZES)
# ------------------------------------------------------------------ FIX A2 (2026-09-15): full-precision strings, as Main_Fig3 panel a
# bdsp_diseases_v3.csv (step 111) stores 3 decimals. The same fit at full precision is the step-168 checkpoint (dC_side_s08_hr_matrix,
# dC_side_analyses/_work/hr_matrix_ck/spo2_pct_below_90.json: the 52 ranked outcomes, HR/lo/hi equal to the csv at 3 dp). The paper's
# rule is half up on full precision, and Main_Fig3 panel a already prints it, so every row present in the checkpoint prints (and draws)
# the checkpoint values; rows absent from it (conditions outside the 52 ranked outcomes) keep the csv values, half up, and are listed.
MATRIX_CK = f"{v14lib.SV}/dC_side_analyses/_work/hr_matrix_ck/spo2_pct_below_90.json"
MATRIX_CSV = f"{v14lib.SV}/dC_side_analyses/hr_matrix_141x52.csv"      # the declared output of step 168 (carries the sidecar)
PRINT_3DP = {c: hr_ci_text(hr, lo, hi) for _g, c, hr, lo, hi, _k, _q in rows}   # half up on the 3-dp csv (the pre-fix strings)
CKINFO = None
if not OLDMODE:
    hydrated(MATRIX_CK); SC168 = sidecar(MATRIX_CSV)
    # round 49 (as round 38, Main_Fig3): step 168 is RESUMABLE (one checkpoint per measurement in _work/hr_matrix_ck/), the checkpoint's
    # mtime (2026-09-14 15:37) predates the 2026-09-16 18:02 re-declaration of the matrix built from it, so the round-37 mtime window no
    # longer holds; the content guard below (every checkpoint HR, lo, hi equals the current bdsp_diseases_v3.csv within 1e-5) replaces it
    CK_WINDOW_S = os.path.getmtime(MATRIX_CSV) - os.path.getmtime(MATRIX_CK); print(f"checkpoint to matrix mtime window {CK_WINDOW_S / 3600:.1f} h (reported, content guard applied)")
    CKL = json.load(open(MATRIX_CK)); assert all(r["measurement"] == "spo2_pct_below_90" for r in CKL)
    _ckd = {r["outcome"]: r for r in CKL}; _k2n = {r.key: r.disease for r in BD.itertuples()}
    _guard = max(abs(float(_ckd[k][a]) - float(BDI.loc[_k2n[k], b])) for k in _ckd if k in _k2n for a, b in (("HR", "t90_hr"), ("lo", "t90_lo"), ("hi", "t90_hi")))
    assert _guard < 1e-5, ("checkpoint differs from bdsp_diseases_v3.csv at full precision: STOP", _guard)
    print(f"content guard: checkpoint = bdsp_diseases_v3.csv on the 52 x 3 values, max |diff| {_guard:.2e}")
    CK = {r["outcome"]: r for r in CKL}; assert len(CK) == 52, len(CK)
    NAME2KEY = {r.disease: r.key for r in BD.itertuples()}
    used, absent, drift = [], [], []
    rows2 = []
    for (grp, c, hr, lo, hi, k, q) in rows:
        key = NAME2KEY[c]
        if key in CK:
            rec = CK[key]
            for a, b in (("HR", hr), ("lo", lo), ("hi", hi)):
                if abs(float(rec[a]) - float(b)) > 5e-4: drift.append((c, a, rec[a], b))
            rows2.append((grp, c, float(rec["HR"]), float(rec["lo"]), float(rec["hi"]), k, q)); used.append(c)
        else:
            rows2.append((grp, c, hr, lo, hi, k, q)); absent.append(c)
    assert not drift, ("checkpoint drifted from the csv beyond 5e-4: STOP", drift)
    assert len(used) == 52 and len(absent) == len(rows) - 52, (len(used), len(absent))
    rows = rows2
    CK_ROWS = {c: NAME2KEY[c] for c in used}
    CKINFO = dict(path=MATRIX_CK, sha256=sha256(MATRIX_CK), step168_declared_output=SC168,
                  rule="half up on the step-168 full-precision checkpoint (equal to bdsp_diseases_v3.csv at 3 dp on all 52 x 3 values, max |diff| below 5e-4, asserted)",
                  rows_from_checkpoint=used, rows_not_in_checkpoint_csv_halfup=absent)
else:
    CK_ROWS = {}
PRINT = {c: hr_ci_text(hr, lo, hi) for _g, c, hr, lo, hi, _k, _q in rows}
if CKINFO is not None:
    CKINFO["moved_strings"] = {c: dict(csv_3dp=PRINT_3DP[c], checkpoint=PRINT[c]) for c in PRINT if PRINT[c] != PRINT_3DP[c]}
    print(f"FIX A2: {len(CK_ROWS)} rows printed from the full-precision checkpoint, {len(CKINFO['rows_not_in_checkpoint_csv_halfup'])} rows not in it (csv half up): {CKINFO['rows_not_in_checkpoint_csv_halfup']}; "
          f"{len(CKINFO['moved_strings'])} strings move against half up on the 3-dp csv: {CKINFO['moved_strings']}")
# display labels: the V13 label column is 119.8 pt wide (x 97.15 to the axes box); a v8.1 name wider than that is shortened for display
# (owner question: 'Ventricular arrhythmia or cardiac arrest' is 172.5 pt at 10 pt Arial and cannot fit on one row of 10.44 pt)
LABEL_SHORT = {"Ventricular arrhythmia or cardiac arrest": "Ventricular arrhythmia/arrest"}
def label_of(c): return LABEL_SHORT.get(c, c)
HALVES = sorted(c for _g, c, hr, lo, hi, _k, _q in rows if any(exact_half_2dp(v) for v in (hr, lo, hi)))


ORDER_C_TITLES = list(json.load(open(hydrated(f"{NUM}/fig2c_selection_v8.json")))["selected_labels"])

# ------------------------------------------------------------------ round-54 geometry (measured text widths, no V13 pin)
G = load_geometry(); A = [r for r in G["records"] if r["region"] == "a"]
def w_(b): return b[2] - b[0]
def h_(b): return b[3] - b[1]
def tw(s, size, weight="normal"): return text_width_pt(None, s, size, weight)
S_LAB, S_VAL, S_HEAD, S_CAP, S_TICK, S_KEY = SIZE["label"], SIZE["value"], SIZE["head"], SIZE["title"], SIZE["tick"], SIZE["key"]
HDR_ONE = {blk: " ".join(lines) for blk, lines in HDR_LINES.items()}          # the V13 two-line headers on one line (declared)
M_L = 14.17                                  # left margin: header rows, row labels, key (the V13 band x)
TITLE_Y, LETTER_Y = 18.0, 44.0               # title strip and letters at 18 pt (stamped by the composer)
KEY_Y = [66.0, 84.0]                         # key rows (baselines), 14 pt
HEAD_Y = 103.0                               # 'HR (95% CI)' and 'q' column headers (baseline), 15 pt
Y0 = 117.0                                   # centre of the first line
PITCH, BLOCK_GAP = 13.9, 7.5                 # 14 pt Arial glyph box (cap to descender) is 12.99 pt
TOP_PAD, BOT_PAD = 8.0, 7.0
W_AX = 158.0                                 # x axis width (V30 144.5): the 14 pt tick labels 0.6 and 0.8 then clear by 3.5 pt
REF_CLEAR, MARK_CLEAR = 4.0, 3.0             # text end to the dashed reference line, text end to the row's own interval
HR_GAP, Q_GAP = 8.0, 7.0                     # widest interval end / axis end to the HR column, HR column to the q column
LAB_DY = centre_dy(pt(S_LAB))                # label baseline below the line centre (text centred on the marker)
# x axis: the V30 range and tick set (round 37 extended the V13 axis to hold every v8.1 interval, ticks 0.6 to 4)
_lo_needed = min(lo for _g, c, hr, lo, hi, _k, _q in rows); _hi_needed = max(hi for _g, c, hr, lo, hi, _k, _q in rows)
XLO, XHI = 0.58, max(4.189078807365091, _hi_needed * 1.02); XTICKS = [0.6, 0.8, 1.0, 1.5, 2.0, 3.0, 4.0]
assert XLO < _lo_needed and _hi_needed < XHI, (_lo_needed, _hi_needed)
R49A = json.load(open(os.path.join(WORK, "r49_panel_a_drawn.json")))
assert abs(R49A["axis"]["xlo"] - XLO) < 1e-6 and abs(R49A["axis"]["xhi"] - XHI) < 1e-6 and R49A["axis"]["ticks"] == XTICKS, (R49A["axis"], XLO, XHI)
AXIS_EXTENDED = True
# label column: the widest header line or row label sets the reference line position
widest_label = max((tw(label_of(c), pt(S_LAB)), c) for _g, c, hr, lo, hi, _k, _q in rows)
widest_hdr = max((tw(HDR_ONE[b], pt(S_LAB), "bold"), b) for b in BLOCK_ORDER)
LAB_X = M_L
REF_X = LAB_X + max(widest_label[0], widest_hdr[0]) + REF_CLEAR
b_ = W_AX / (math.log(XHI) - math.log(XLO)); a_ = REF_X - b_ * math.log(1.0)
AX_L, AX_R = a_ + b_ * math.log(XLO), a_ + b_ * math.log(XHI)
def xpos(v): return a_ + b_ * math.log(v)
ticks = [xpos(v) for v in XTICKS]; ref_x = xpos(1.0)
for i in range(len(ticks) - 1):
    gap = (ticks[i + 1] - ticks[i]) - (tw(f"{XTICKS[i]:g}", pt(S_TICK)) + tw(f"{XTICKS[i + 1]:g}", pt(S_TICK))) / 2
    assert gap >= 2.5, ("tick labels too close", XTICKS[i], XTICKS[i + 1], gap)
# every row label clears its own interval and the reference line
for _g, c, hr, lo, hi, _k, _q in rows:
    end = LAB_X + tw(label_of(c), pt(S_LAB))
    assert end <= xpos(lo) - MARK_CLEAR and end <= ref_x - REF_CLEAR, (c, end, xpos(lo), ref_x)
for b in BLOCK_ORDER: assert LAB_X + tw(HDR_ONE[b], pt(S_LAB), "bold") <= ref_x - REF_CLEAR, b
# HR and q columns
HR_W = max(tw(PRINT[c], pt(S_VAL)) for c in PRINT); Q_W = max(tw(q_text(q), pt(S_VAL), "bold" if q < 0.05 else "normal") for _g, c, hr, lo, hi, _k, q in rows)
COL_R = max(AX_R, max(xpos(hi) for _g, c, hr, lo, hi, _k, _q in rows)) + HR_GAP + HR_W
Q_R = COL_R + Q_GAP + Q_W
A_RIGHT = Q_R
SPLIT_X = round(A_RIGHT + 4.0, 2)            # a/b clip split (invisible), 4 pt right of the q column
# lines: header line then rows, block gaps between blocks
LINES = []                                   # (kind, block, row index or None)
for k, blk in enumerate(BLOCK_ORDER):
    LINES.append(("hdr", blk, None))
    for i, r in enumerate(rows):
        if r[0] == blk: LINES.append(("row", blk, i))
assert len(LINES) == len(rows) + len(BLOCK_ORDER) == 65, len(LINES)
LINE_Y = []; y = Y0; prev_blk = None
for kind, blk, i in LINES:
    if kind == "hdr" and prev_blk is not None: y += BLOCK_GAP
    LINE_Y.append(round(y, 3)); y += PITCH; prev_blk = blk
AX_T = Y0 - TOP_PAD; AX_B = LINE_Y[-1] + BOT_PAD
SLOT_Y = {i: LINE_Y[n] for n, (kind, blk, i) in enumerate(LINES) if kind == "row"}
HDR_Y = {blk: LINE_Y[n] for n, (kind, blk, i) in enumerate(LINES) if kind == "hdr"}
TICK_BASE = AX_B + 3.0 + 0.95 * pt(S_TICK)
CAPTION_BASE = TICK_BASE + 1.35 * pt(S_CAP)
A_BOTTOM = CAPTION_BASE + DESC * pt(S_CAP)
# key geometry: the V13 marks (line length, marker size) with the text at the round-54 widths
kmk = [r for r in A if r["op"] == "B" and 0 < w_(r["bbox"]) < 12 and r["bbox"][0] < 180]
kln = [r for r in A if r["op"] == "S" and h_(r["bbox"]) < 0.05 and r["bbox"][0] < 180]
def kline_len(col, dashed_=False): c = [r for r in kln if hexcol(r["stroke"]) == col and bool(r["dash"]) == dashed_]; assert len(c) == 1, (col, c); return w_(c[0]["bbox"])
KL_SIG, KL_NS, KL_CTRL, KL_FLOOR = kline_len(BLUE), kline_len(PALE), kline_len(GREY), kline_len(GREY, True)
KEY_TXT_GAP = 4.4                            # V13: text origin 4.4 pt after the key line end
KEY_COL_GAP = 14.0
KEY1_LX = M_L; KEY1_TX = KEY1_LX + KL_SIG + KEY_TXT_GAP
KEY2_TX = KEY1_TX + max(tw("Significant after FDR", pt(S_KEY)), tw("Not significant", pt(S_KEY))) + KEY_COL_GAP + KL_CTRL + KEY_TXT_GAP
KEY2_LX = KEY2_TX - KEY_TXT_GAP - KL_CTRL
assert KEY2_TX + tw(KEY_FLOOR, pt(S_KEY)) <= A_RIGHT, "key row 2 runs past the q column"
MARK_W = {"sig": 5.29, "ns": 5.29, "ctrl": 4.8}
def kmark_w(col): c = [r for r in kmk if hexcol(r["fill"]) == col]; assert len(c) == 1, (col, c); return w_(c[0]["bbox"])
KEY_MARK_W = {"sig": kmark_w(BLUE), "ns": kmark_w(PALE), "ctrl": kmark_w(GREY)}   # the V13 key markers (6.0, 6.0, 5.4), larger than the row markers
mk = [r for r in A if r["op"] == "B" and 0 < w_(r["bbox"]) < 12 and h_(r["bbox"]) < 12 and r["bbox"][0] > 200]
for cls, col in (("sig", BLUE), ("ns", PALE), ("ctrl", GREY)):
    ws = [w_(r["bbox"]) for r in mk if hexcol(r["fill"]) == col]; assert all(abs(w - MARK_W[cls]) < 0.03 for w in ws), (cls, sorted(set(round(w, 2) for w in ws)))
CI_LW, MEDGE = 1.7, 0.8
# ------------------------------------------------------------------ the rest of the sheet (planner): panel b above panel c on the right
# panel b: key beside the letter (two rows), rows at PITCH_B with the V30 sub-row spacing (6.44 pt, marker sizes unchanged)
B = dict(x0=SPLIT_X, x1=PAGE_W, letter_xy=[SPLIT_X + 4.0, LETTER_Y], key_y=[LETTER_Y, LETTER_Y + 18.0], key_x0=SPLIT_X + 30.0,
         sub_row_spacing=6.44, pitch=27.0, n_slots=18, gutter_pad=4.0, ax_r=PAGE_W - 12.0, ax_t=LETTER_Y + 18.0 + DESC * pt(S_KEY) + 10.0)
B["ax_b"] = B["ax_t"] + B["pitch"] * B["n_slots"]
B["tick_base"] = B["ax_b"] + 3.0 + 0.95 * pt(S_TICK); B["caption_base"] = B["tick_base"] + 1.35 * pt(S_CAP); B["bottom"] = B["caption_base"] + DESC * pt(S_CAP)
BC_Y = round(B["bottom"] + 8.0, 2)
# panel c: 3 columns x 4 rows, boxes BW x H_BOX, titles centred over the box, the key in the free slots of the last row
C = dict(x0=SPLIT_X, x1=PAGE_W, y0=BC_Y, letter_xy=[SPLIT_X + 4.0, BC_Y + 14.0], ncol=3, nrow=4, bw=120.0, h_box=72.0, col_gap=24.5,
         ylabel_x=SPLIT_X + 15.5, title_dy=4.8, xl_dy=16.0, yt_pad=7.3, row_gap=7.0)
C["title_cap"] = CAP * pt(S_HEAD)
WIDEST_C_TITLE = max(tw(t, pt(S_HEAD)) for t in ORDER_C_TITLES)
C["box_x0"] = C["ylabel_x"] + DESC * pt(S_CAP) + 4.0 + (WIDEST_C_TITLE - C["bw"]) / 2   # first column: the widest centred title clears the y label
C["first_title_base"] = C["letter_xy"][1] + 22.0
C["first_box_y0"] = C["first_title_base"] + C["title_dy"]
C["pitch_row"] = C["h_box"] + C["xl_dy"] + 3.0 + C["row_gap"] + C["title_cap"] + C["title_dy"]
C["last_box_y1"] = C["first_box_y0"] + (C["nrow"] - 1) * C["pitch_row"] + C["h_box"]
C["xtitle_base"] = C["last_box_y1"] + C["xl_dy"] + 3.0 + 4.0 + CAP * pt(S_CAP)
C["bottom"] = C["xtitle_base"] + DESC * pt(S_CAP)
PAGE_H_NEW = round(max(A_BOTTOM, C["bottom"]) + 8.0, 2)
assert PAGE_H_NEW <= 1110.0, ("page taller than the 1110 pt aim", PAGE_H_NEW, A_BOTTOM, C["bottom"])
assert C["box_x0"] + C["ncol"] * C["bw"] + (C["ncol"] - 1) * C["col_gap"] + (WIDEST_C_TITLE - C["bw"]) / 2 <= PAGE_W, "panel c runs off the sheet"
json.dump(dict(page_h=PAGE_H_NEW, page_w=PAGE_W, split_x=SPLIT_X, bc_y=BC_Y, title_y=TITLE_Y, letter_y=LETTER_Y, letter_size=pt(SIZE["letter"]),
               a=dict(lab_x=LAB_X, ref_x=ref_x, ax_box=[AX_L, AX_T, AX_R, AX_B], w_ax=W_AX, xlo=XLO, xhi=XHI, ticks=XTICKS, tick_x=ticks, pitch=PITCH, block_gap=BLOCK_GAP,
                      y0=Y0, key_y=KEY_Y, head_y=HEAD_Y, col_r=COL_R, q_r=Q_R, right=A_RIGHT, bottom=A_BOTTOM, n_lines=len(LINES), n_rows=len(rows),
                      block_sizes=dict(zip(BLOCK_ORDER, BLOCK_SIZES)), headers_one_line=HDR_ONE, tick_base=TICK_BASE, caption_base=CAPTION_BASE),
               b=B, c=C, r49_layout=json.load(open(os.path.join(WORK, "r49_layout.json")))), open(os.path.join(WORK, "layout.json"), "w"), indent=1)
print(f"layout: page {PAGE_W} x {PAGE_H_NEW}, split x {SPLIT_X}, b/c split y {BC_Y}; a: ref line {ref_x:.2f}, axes {AX_L:.2f}-{AX_R:.2f} x {AX_T:.2f}-{AX_B:.2f}, HR right {COL_R:.2f}, q right {Q_R:.2f}, bottom {A_BOTTOM:.2f}; "
      f"widest label {widest_label}, widest header {widest_hdr}; b: axes top {B['ax_t']:.1f} bottom {B['ax_b']:.1f} caption {B['caption_base']:.1f}; c: boxes from y {C['first_box_y0']:.1f} to {C['last_box_y1']:.1f}, x title {C['xtitle_base']:.1f}")

# ------------------------------------------------------------------ draw
COL = {"sig": BLUE, "ns": PALE, "ctrl": GREY}; MARKER = {"sig": "o", "ns": "o", "ctrl": "s"}
drawn = dict(panel="a", mode="new", sources=dict(bdsp_diseases_v3=BD_PATH, negcontrols_v4_panel=NCP_PATH, cohorts=COH_PATH), provenance=prov,
             n_cohort=N_COHORT, n_outcomes=N_OUT, n_fdr_sig=N_SIG, n_fdr_raised=N_RAISED, below_one_sig=BELOW, controls_fdr_sig=CTRL_SIG, control_q={c: float(Q[c]) for c in ctrl.disease},
             floor=FLOOR, floor_set_by=FLOOR_BY, key_floor=KEY_FLOOR, stored_precision="3 decimals in bdsp_diseases_v3.csv and negcontrols_v4_panel.csv; the 52 ranked outcomes print half up on the step-168 full-precision checkpoint (FIX A2), the other rows half up on the csv",
             fullprec_checkpoint=CKINFO,
             exact_half_rows=HALVES, new_to_sheet=NEW_TO_SHEET, left_sheet=LEFT_SHEET, blocks={b: blocks_new[b] for b in BLOCK_ORDER}, v13_blocks=V13_BLOCK,
             layout=dict(page_h=PAGE_H_NEW, pitch=PITCH, block_gap=BLOCK_GAP, ax_box=[AX_L, AX_T, AX_R, AX_B], ax_box_v30=R49A["layout"]["ax_box"], lines=LINES, line_y=LINE_Y),
             axis=dict(box=[AX_L, AX_T, AX_R, AX_B], xlo=XLO, xhi=XHI, extended=AXIS_EXTENDED, ticks=XTICKS, tick_x=ticks, ref_x=ref_x, floor_x=xpos(FLOOR), w_ax=W_AX),
             label_short=LABEL_SHORT,
             r54=dict(plus=PLUS, lab_x=LAB_X, q_r=Q_R, col_r=COL_R, right=A_RIGHT, split_x=SPLIT_X, headers_one_line=HDR_ONE, header_form="bold row above each block (the V13 left band does not fit beside 14 pt labels)",
                      widest_label=widest_label, widest_header=widest_hdr, ref_clear=REF_CLEAR, mark_clear=MARK_CLEAR, key=dict(row_y=KEY_Y, key_mark_w=KEY_MARK_W, col1_line_x=KEY1_LX, col1_text_x=KEY1_TX, col2_line_x=KEY2_LX, col2_text_x=KEY2_TX)),
             texts=[], marks=[])
def T_(x, y, s, size, ha="left", weight="normal", **rec):
    size = pt(size)   # round 54: +4.0 pt, drawn and recorded
    ov.text(x, y, s, fontsize=size, color=INK, ha=ha, va="baseline", fontweight=weight)
    drawn["texts"].append(drawn_rec(s, x, y, ha=ha, size=size, weight=weight, panel="a", **rec))
with plt.rc_context(RC):
    fig, ov = page_figure(PAGE_H_NEW)
    ov.add_patch(matplotlib.patches.Rectangle((AX_L, AX_T), AX_R - AX_L, AX_B - AX_T, facecolor="#ffffff", edgecolor="none", zorder=0))
    ov.plot([ref_x, ref_x], [AX_T, AX_B], color=INK, lw=0.9, ls=(0, (4, 3)), zorder=1, solid_capstyle="butt")
    ov.plot([xpos(FLOOR)] * 2, [AX_T, AX_B], color=GREY, lw=0.9, ls=(0, (1.6, 1.8)), zorder=1, solid_capstyle="butt")
    ov.plot([AX_L, AX_R], [AX_B, AX_B], color=INK, lw=0.8, zorder=4, solid_capstyle="projecting")
    for x, v in zip(ticks, XTICKS):
        ov.plot([x, x], [AX_B, AX_B + 3.0], color=INK, lw=0.8, zorder=4, solid_capstyle="butt")
        T_(x, TICK_BASE, f"{v:g}", S_TICK, ha="center", source_key="x tick (V30 axis: extended in round 37 to hold every v8.1 interval)", rule="text")
    for i, (grp, c, hr, lo, hi, k, q) in enumerate(rows):
        y = SLOT_Y[i]; base = y + LAB_DY
        ov.plot([xpos(lo), xpos(hi)], [y, y], color=COL[k], lw=CI_LW, zorder=2)
        ov.plot([xpos(hr)], [y], marker=MARKER[k], ms=MARK_W[k], ls="None", markerfacecolor=COL[k], markeredgecolor="white", markeredgewidth=MEDGE, zorder=3)
        src = NCP_PATH if k == "ctrl" else BD_PATH
        if c in LABEL_SHORT: T_(LAB_X, base, label_of(c), S_LAB, source_file=src, source_key=f"display label of '{c}' (shortened to fit the label column, owner question)", source_value=label_of(c), rule="text")
        else: T_(LAB_X, base, c, S_LAB, source_file=src, source_key=(f"condition={c}:condition" if k == "ctrl" else f"disease={c}:disease"), source_value=c, rule="text")
        if c in CK_ROWS: T_(COL_R, base, PRINT[c], S_VAL, ha="right", source_file=MATRIX_CK, source_key=f"[outcome={CK_ROWS[c]}]", source_value={"HR": hr, "lo": lo, "hi": hi}, rule="hr_ci_HRlohi", csv_3dp=PRINT_3DP[c], csv_source=f"{src}: disease={c}:t90_hr,t90_lo,t90_hi (equal at 3 dp)")
        else: T_(COL_R, base, PRINT[c], S_VAL, ha="right", source_file=src, source_key=(f"condition={c}:hr,lo,hi" if k == "ctrl" else f"disease={c}:t90_hr,t90_lo,t90_hi"), source_value=[hr, lo, hi], rule="hr_ci_2dp")
        T_(Q_R, base, q_text(q), S_VAL, ha="right", weight="bold" if q < 0.05 else "normal", source_file=BD_PATH,
           source_key=f"BH q of t90_p, computed across the {'5 controls' if k == 'ctrl' else str(N_OUT) + ' non-control outcomes'} (disease={c})", source_value=q, rule="q3")
        drawn["marks"].append(dict(condition=c, block=grp, cls=k, hr=hr, lo=lo, hi=hi, q=q, printed=PRINT[c], x_page=round(xpos(hr), 3), y_page=y, x_lo=round(xpos(lo), 3), x_hi=round(xpos(hi), 3)))
    for blk in BLOCK_ORDER:
        T_(LAB_X, HDR_Y[blk] + LAB_DY, HDR_ONE[blk], S_LAB, weight="bold", source_key="organ block header (V13 words, one line, a bold row above the block)", rule="text", v13_lines=HDR_LINES[blk])
    T_(COL_R, HEAD_Y, "HR (95% CI)", S_HEAD, ha="right"); T_(Q_R, HEAD_Y, "q", S_HEAD, ha="right")
    T_((AX_L + AX_R) / 2, CAPTION_BASE, "Hazard ratio per 1 SD of T90 (95% CI)", S_CAP, ha="center")
    ky = [k - centre_dy(pt(S_KEY)) for k in KEY_Y]
    for (lx, ln, y, col, mkr, ms) in ((KEY1_LX, KL_SIG, ky[0], BLUE, "o", KEY_MARK_W["sig"]), (KEY1_LX, KL_NS, ky[1], PALE, "o", KEY_MARK_W["ns"]), (KEY2_LX, KL_CTRL, ky[0], GREY, "s", KEY_MARK_W["ctrl"])):
        ov.plot([lx, lx + ln], [y, y], color=col, lw=CI_LW, zorder=2)
        ov.plot([lx + ln / 2], [y], marker=mkr, ms=ms, ls="None", markerfacecolor=col, markeredgecolor="white", markeredgewidth=MEDGE, zorder=3)
    ov.plot([KEY2_LX, KEY2_LX + KL_FLOOR], [ky[1], ky[1]], color=GREY, lw=0.9, ls=(0, (1.6, 1.8)), zorder=2, solid_capstyle="butt")
    T_(KEY1_TX, KEY_Y[0], "Significant after FDR", S_KEY); T_(KEY1_TX, KEY_Y[1], "Not significant", S_KEY); T_(KEY2_TX, KEY_Y[0], "Negative controls", S_KEY)
    T_(KEY2_TX, KEY_Y[1], KEY_FLOOR, S_KEY, source_file=NCP_PATH, source_key="key sentence: max of hr over the 5 control rows, half up 2dp (not a pointer)", source_value=KEY_FLOOR, rule="text", floor_value=FLOOR)
    out = os.path.join(WORK, "panel_a.pdf"); fig.savefig(out, format="pdf"); plt.close(fig)
json.dump(drawn, open(os.path.join(WORK, "panel_a_drawn.json"), "w"), indent=1)
print(f"panel_a.pdf written: {len(rows)} rows in blocks {dict(zip(BLOCK_ORDER, BLOCK_SIZES))}, FDR-significant {N_SIG} of {N_OUT} ({N_RAISED} above 1, below 1: {BELOW}), "
      f"floor {FLOOR:.3f} ({FLOOR_BY}) -> '{KEY_FLOOR}', controls FDR-significant: {CTRL_SIG or 'none'}; pitch {PITCH}, gap {BLOCK_GAP}, page {PAGE_H_NEW}; exact-half rows (3dp source): {HALVES}")
