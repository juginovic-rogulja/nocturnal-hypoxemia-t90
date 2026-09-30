#!$T90_PY
"""Round 54 (2026-09-29), lane L3, Main_Fig3: the numbers-only layout planner (no drawing). Derives the width budget of the b + c row
(the binding constraint), the two forest columns (label gutter, box, HR and q columns), the vertical budget with the round-49 band as
the stand-in, the whisker lengths of the narrowest intervals per axis choice, and prints everything against the 1260 pt gate and the
1090 pt aim. The builder (01_build_abc_r54.py) imports PLAN from here so the drawn geometry is the planned geometry."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import json, math, os, sys
import fitz
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lane_common as C

ARIAL, ARIALB = C.ARIAL, C.ARIALB
F, FB = fitz.Font(fontfile=ARIAL), fitz.Font(fontfile=ARIALB)
W = 952.73
R49_DRAWN = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/L3/Main_Fig3/work/Main_Fig3_drawn.json"
DJ49 = json.load(open(C.hydrated(R49_DRAWN)))
ROWS = DJ49["panel_a"]["rows"]; assert len(ROWS) == 52
BAND49 = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LSKETCH/work/build/band_3d_r49.pdf"

# ---------------------------------------------------------------- sizes (round-54 targets by class; V30 size in the comment)
SZ = dict(row=14.0, tick=14.0, hrq=14.0, key=14.0, hdr=15.0, cap=15.0,            # panel a: rows/ticks/key 11.71, HR/q strings and headers 11.17, caption 12.78
          b_title=15.0, b_hdr=15.0, b_cell=14.0, b_row=14.0, b_key=14.0,            # panel b: title 13.65, headers/cells/rows/key 12.5
          c_hdr=15.0, c_sub=13.0, c_row=14.0, c_rowhdr=15.0, c_cell=13.0, c_note=13.0,   # panel c: duration headers 10.91, sub-headers/cells 10.5, row labels 10.5, notes 10.09
          letter=18.0, title=18.0)
def w(s, size, bold=False): return (FB if bold else F).text_length(s, fontsize=size)

# ---------------------------------------------------------------- panel a: two columns of 26 rows
N_COL = 26; LEFT, RIGHT = ROWS[:N_COL], ROWS[N_COL:]
LAB_X0 = 8.88                                                     # the V30 label origin (letter a at x 8)
CLEAR_LABEL = 3.0                                                  # labels end at least 3 pt before the forest box (the round-40 rule)
GAP_BOX_HR = 2.4                                                   # HR string starts this far right of the box (V30: 2.35)
GAP_HR_Q = 4.8                                                     # q string starts this far right of the HR string's right edge (V30: 4.77)
GAP_COL = 20.0                                                     # between the left column's q right edge and the right column's label origin
RIGHT_MARGIN = 8.85                                                # V30: the rightmost string ends at 943.88
HR_W = max(w(r["s_hr"], SZ["hrq"]) for r in ROWS); Q_W = max(w(r["s_q"], SZ["hrq"], bold=r["sig"]) for r in ROWS)
lab_w = {r["label"]: w(r["label"], SZ["row"]) for r in ROWS}
L_LAB = max(lab_w[r["label"]] for r in LEFT); R_LAB = max(lab_w[r["label"]] for r in RIGHT)
L_LONGEST = max(LEFT, key=lambda r: lab_w[r["label"]])["label"]; R_LONGEST = max(RIGHT, key=lambda r: lab_w[r["label"]])["label"]
AFTER_BOX = GAP_BOX_HR + HR_W + GAP_HR_Q + Q_W                     # box right edge to the q right edge
fixed = LAB_X0 + (L_LAB + CLEAR_LABEL) + AFTER_BOX + GAP_COL + (R_LAB + CLEAR_LABEL) + AFTER_BOX + RIGHT_MARGIN
# axis ranges per column: the left 26 rows span 0.63 to 3.78 (the V30 range and ticks), the right 26 rows 0.61 to 1.34 (ticks 0.6, 0.8, 1, 1.5)
XLO = 0.58
L_VALS = [v for r in LEFT for v in (r["lo"], r["hi"], r["tst_lo"], r["tst_hi"])]; R_VALS = [v for r in RIGHT for v in (r["lo"], r["hi"], r["tst_lo"], r["tst_hi"])]
XHI_L = 4.15; XHI_R = 1.55
TICKS_L = [0.6, 0.8, 1, 1.5, 2, 3, 4]; TICKS_R = [0.6, 0.8, 1, 1.5]
assert XLO < min(L_VALS) and max(L_VALS) < XHI_L and XLO < min(R_VALS) and max(R_VALS) < XHI_R, (min(L_VALS), max(L_VALS), min(R_VALS), max(R_VALS))
LN_L, LN_R = math.log(XHI_L / XLO), math.log(XHI_R / XLO)
BOX_TOTAL = math.floor((W - fixed) * 10) / 10                      # the two boxes share what the fixed elements leave
BOX_L = math.floor(BOX_TOTAL * LN_L / (LN_L + LN_R) * 10) / 10; BOX_R = round(BOX_TOTAL - BOX_L, 1)   # equal scale: one ln unit is the same length in both columns
L_BOX0 = math.ceil((LAB_X0 + L_LAB + CLEAR_LABEL) * 10) / 10; L_BOX1 = L_BOX0 + BOX_L
L_HR_R = L_BOX1 + GAP_BOX_HR + HR_W; L_Q_R = L_HR_R + GAP_HR_Q + Q_W
R_LAB_X0 = round(L_Q_R + GAP_COL, 2); R_BOX0 = math.ceil((R_LAB_X0 + R_LAB + CLEAR_LABEL) * 10) / 10; R_BOX1 = R_BOX0 + BOX_R
R_HR_R = R_BOX1 + GAP_BOX_HR + HR_W; R_Q_R = R_HR_R + GAP_HR_Q + Q_W
PT_LN_L = BOX_L / LN_L; PT_LN_R = BOX_R / LN_R; PT_LN_SHARED = BOX_TOTAL / 2 / LN_L; PT_LN_V30 = float(DJ49["panel_a"]["geometry"]["XB"])
BOX_W = (BOX_L, BOX_R)
TICK_GAP_L = min((math.log(b / a) * PT_LN_L - (w(f"{a:g}", SZ["tick"]) + w(f"{b:g}", SZ["tick"])) / 2) for a, b in zip(TICKS_L[:-1], TICKS_L[1:]))   # the clear space between neighbouring tick labels
TICK_GAP_R = min((math.log(b / a) * PT_LN_R - (w(f"{a:g}", SZ["tick"]) + w(f"{b:g}", SZ["tick"])) / 2) for a, b in zip(TICKS_R[:-1], TICKS_R[1:]))
assert TICK_GAP_L >= 3.0 and TICK_GAP_R >= 3.0, ("tick labels closer than 3 pt", TICK_GAP_L, TICK_GAP_R)
def narrowest(rows, ptln):
    m = min(rows, key=lambda r: math.log(r["hi"] / r["lo"])); return m["label"], m["s_hr"], round(math.log(m["hi"] / m["lo"]) * ptln, 2)
MARK_D = 2 * math.sqrt((math.sqrt(28.0) * 1.0706) ** 2 / math.pi)  # the T90 circle's diameter (scatter s in pt^2)

# ---------------------------------------------------------------- panels b and c side by side (the binding width constraint)
B_X0 = 8.88                                                        # the letter b at x 8, the block from the label origin
B_LAB_W = max(w(s, SZ["b_row"]) for s in ("<5 h", "6–7 h", ">7 h")); B_LAB_GAP = 8.0; B_GUT = math.ceil(B_LAB_W + 1 + B_LAB_GAP)
B_CELL_TXT = max(w(f"{n:,} ({p})", SZ["b_cell"]) for r_ in DJ49["panel_b"]["rows"].values() for n, p in [(c["n"], c["pct"]) for c in r_["cells"]])
B_HDR_TXT = max(w(s, SZ["b_hdr"]) for s in ("≤1%", ">1–5%", ">5–10%", ">10%"))
B_GAP_H, B_GAP_V = 6.0, 4.6                                        # V30: 7.2 and 4.6
C_LAB_W = max(w(s, SZ["c_row"]) for s in ("COPD", "Type 2 diabetes", "Heart failure", "Death, any cause", "Cardiovascular", "composite", "Dyslipidemia", "Alopecia", "(negative control)"))
C_ROWHDR_W = max(w(s, SZ["c_rowhdr"]) for s in ("Sleep duration", "Oxygenation"))
C_LAB_GAP = 9.0; C_GUT = math.ceil(max(C_LAB_W, C_ROWHDR_W) + C_LAB_GAP)
C_CELL_TXT = max(w(s, SZ["c_cell"]) for t in DJ49["panel_c"]["rows"].values() for c in t["cells"].values() for s in (c["s_hr"], c["s_ci"]))
C_SUB_TXT = max(w(s, SZ["c_sub"]) for s in ("Preserved", "Low", "T90 ≤1%", "T90 >10%"))
C_PAIR_GAP, C_GROUP_GAP = 1.6, 7.0                                 # V30: 1.6 and 10.41
BC_GAP = 12.0
C_CELL_MARGIN = 1.5                                                # per side, the coordinator accepted 0.8 in round 49
C_CELL = math.ceil((C_CELL_TXT + 2 * C_CELL_MARGIN) * 10) / 10
C_W = C_GUT + 3 * (2 * C_CELL + C_PAIR_GAP) + 2 * C_GROUP_GAP
B_CELL = math.floor(((W - RIGHT_MARGIN) - B_X0 - B_GUT - 3 * B_GAP_H - BC_GAP - C_W) / 4 * 10) / 10
B_W = B_GUT + 4 * B_CELL + 3 * B_GAP_H
C_X0 = round(B_X0 + B_W + BC_GAP, 2); C_CELL_X0 = C_X0 + C_GUT
C_RIGHT = C_CELL_X0 + 3 * (2 * C_CELL + C_PAIR_GAP) + 2 * C_GROUP_GAP
B_MARGIN = (B_CELL - B_CELL_TXT) / 2; B_HDR_MARGIN = (B_CELL - B_HDR_TXT) / 2; C_MARGIN = (C_CELL - C_CELL_TXT) / 2; C_SUB_MARGIN = (C_CELL - C_SUB_TXT) / 2
# b key: lead-in + rungs (swatch 12.42, label 17.42 after the swatch, next swatch 20 after the label's end)
B_KEY_LEAD_GAP = 17.8                                              # V30: the first swatch 17.8 pt after the lead-in text
KEY_BANDS = DJ49["panel_b"]["key_bands"]; kx = B_X0 + w("Percent of the group", SZ["b_key"]) + B_KEY_LEAD_GAP
for lab in KEY_BANDS: kx = kx + 17.42 + w(lab, SZ["b_key"]) + 20.0
B_KEY_END = kx - 20.0
# c legend: bar 124.24 x 9.92 at cell_x0 + 25, reference box 14.5 x 10.74, texts 13 pt
C_LEG_BAR0 = C_CELL_X0 + 25.0; C_LEG_BAR1 = C_LEG_BAR0 + 124.24; C_REF_X = C_LEG_BAR0 + 150.0
C_LEG_END = max(C_REF_X + 14.5 + 5.0 + w("Reference, 6–7 h and T90 ≤1%", SZ["c_note"]), C_REF_X + w("* q < 0.05, ** q < 0.01, *** q < 0.001", SZ["c_note"]))

# ---------------------------------------------------------------- vertical budget (y down)
ASC = 0.905                                                        # Arial Bold ascent (the letters), cap height 0.716
TITLE_BASE = 18.0; A_TOP = 30.0; A_BASE = round(A_TOP + ASC * SZ["letter"], 2)
KEY_BASE = 47.5; KEY_Y = KEY_BASE - 4.3; KEY_X0 = 30.0            # the key marks on the x-height centre of 14-pt text, the row starts 12 pt right of the letter a
HDR_BASE = 66.5
ROW0 = 81.5; PITCH = 14.5; LANE = 3.06 * (PITCH / 12.719) if False else 3.06   # the two lanes (T90 above, sleep time below) keep the V30 offset
LAST_C = ROW0 + (N_COL - 1) * PITCH; AX_TOP = ROW0 - 0.75 * PITCH; AX_BOT = LAST_C + 0.75 * PITCH
TICK_LEN = 3.0 * 1.0706; TICK_BASE = round(AX_BOT + TICK_LEN + 13.5, 3); CAP_BASE = round(TICK_BASE + 18.5, 3)
A_BOTTOM = CAP_BASE + 0.212 * SZ["cap"]
B_TOP = round(A_BOTTOM + 10.0, 1); B_BASE = round(B_TOP + ASC * SZ["letter"], 2)   # letters b and c
# b: title beside the letter, header cells, three rows, key
B_TITLE_BASE = B_BASE + 0.7; B_HDR_Y0 = B_TITLE_BASE + 7.0; B_HDR_H = 19.0; B_HDR_BASE = B_HDR_Y0 + B_HDR_H / 2 + 0.358 * SZ["b_hdr"]
B_ROW_Y0 = B_HDR_Y0 + B_HDR_H + 5.0; B_ROW_H = 52.0; B_ROW_PITCH = B_ROW_H + B_GAP_V
B_KEY_BASE = B_ROW_Y0 + 2 * B_ROW_PITCH + B_ROW_H + 24.0; B_KEY_SW_Y = B_KEY_BASE - 8.08   # swatch 10.76 tall, V30: top at baseline - 8.08
B_BOTTOM = B_KEY_BASE + 0.212 * SZ["b_key"]
# c: header rows at the letter's height, seven rows, legend
C_HDR1_Y0 = B_TOP + 13.0; C_HDR1_H = 26.0; C_HDR2_Y0 = C_HDR1_Y0 + C_HDR1_H + 1.5; C_HDR2_H = 34.0   # the header rows start below the letter c so that "Sleep duration" (right-aligned 2 pt right of the letter) clears it vertically
C_ROW_Y0 = C_HDR2_Y0 + C_HDR2_H + 2.5; C_ROW_PITCH = 34.0; C_ROW_H = 32.0; C_LINE_DY = 13.5
C_ROWS_END = C_ROW_Y0 + 6 * C_ROW_PITCH + C_ROW_H
C_LEG_TITLE_BASE = C_ROWS_END + 13.0; C_BAR_Y0 = C_LEG_TITLE_BASE + 4.0; C_BAR_H = 9.914; C_BAR_TICK_BASE = C_BAR_Y0 + C_BAR_H + 11.0
C_REF_BASE = C_LEG_TITLE_BASE; C_REF_Y0 = C_REF_BASE - 10.05; C_AST_BASE = C_REF_BASE + 15.5   # the 14.5 x 10.74 box centred on the x-height centre of its 13-pt text
C_BOTTOM = max(C_BAR_TICK_BASE, C_AST_BASE) + 0.212 * SZ["c_note"]
BC_BOTTOM = max(B_BOTTOM, C_BOTTOM)
D_TOP = round(BC_BOTTOM + 10.0, 1); D_BASE = round(D_TOP + ASC * SZ["letter"], 2); BAND_TOP = round(D_TOP + 22.0, 3)
bd = fitz.open(C.hydrated(BAND49)); BAND49_H = bd[0].rect.height; bd.close()
H_R49 = round(BAND_TOP + BAND49_H, 3)

PLAN = dict(W=W, sizes=SZ, n_col=N_COL, lab_x0=LAB_X0, clear_label=CLEAR_LABEL, gap_box_hr=GAP_BOX_HR, gap_hr_q=GAP_HR_Q, gap_col=GAP_COL, right_margin=RIGHT_MARGIN,
            hr_w=HR_W, q_w=Q_W, box_w=BOX_W, box_total=BOX_TOTAL, tick_gap=dict(left=TICK_GAP_L, right=TICK_GAP_R), left=dict(lab_x0=LAB_X0, lab_w=L_LAB, longest=L_LONGEST, box0=L_BOX0, box1=L_BOX1, hr_r=L_HR_R, q_r=L_Q_R, xlo=XLO, xhi=XHI_L, ticks=TICKS_L),
            right=dict(lab_x0=R_LAB_X0, lab_w=R_LAB, longest=R_LONGEST, box0=R_BOX0, box1=R_BOX1, hr_r=R_HR_R, q_r=R_Q_R, xlo=XLO, xhi=XHI_R, ticks=TICKS_R),
            pt_per_ln=dict(left=PT_LN_L, right=PT_LN_R, v30=PT_LN_V30), mark_diameter=MARK_D,
            y=dict(title_base=TITLE_BASE, a_top=A_TOP, a_base=A_BASE, key_base=KEY_BASE, key_y=KEY_Y, key_x0=KEY_X0, hdr_base=HDR_BASE, row0=ROW0, pitch=PITCH, lane=LANE, ax_top=AX_TOP, ax_bot=AX_BOT, tick_len=TICK_LEN, tick_base=TICK_BASE, cap_base=CAP_BASE, a_bottom=A_BOTTOM,
                   b_top=B_TOP, b_base=B_BASE, b_title_base=B_TITLE_BASE, b_hdr_y0=B_HDR_Y0, b_hdr_h=B_HDR_H, b_hdr_base=B_HDR_BASE, b_row_y0=B_ROW_Y0, b_row_h=B_ROW_H, b_row_pitch=B_ROW_PITCH, b_key_base=B_KEY_BASE, b_key_sw_y=B_KEY_SW_Y, b_bottom=B_BOTTOM,
                   c_hdr1_y0=C_HDR1_Y0, c_hdr1_h=C_HDR1_H, c_hdr2_y0=C_HDR2_Y0, c_hdr2_h=C_HDR2_H, c_row_y0=C_ROW_Y0, c_row_pitch=C_ROW_PITCH, c_row_h=C_ROW_H, c_line_dy=C_LINE_DY, c_rows_end=C_ROWS_END,
                   c_leg_title_base=C_LEG_TITLE_BASE, c_bar_y0=C_BAR_Y0, c_bar_h=C_BAR_H, c_bar_tick_base=C_BAR_TICK_BASE, c_ref_y0=C_REF_Y0, c_ref_base=C_REF_BASE, c_ast_base=C_AST_BASE, c_bottom=C_BOTTOM,
                   bc_bottom=BC_BOTTOM, d_top=D_TOP, d_base=D_BASE, band_top=BAND_TOP, band49_h=BAND49_H, h_with_band49=H_R49),
            b=dict(x0=B_X0, gut=B_GUT, lab_w=B_LAB_W, key_lead_gap=B_KEY_LEAD_GAP, cell=B_CELL, gap_h=B_GAP_H, gap_v=B_GAP_V, w=B_W, cell_txt=B_CELL_TXT, margin=B_MARGIN, hdr_txt=B_HDR_TXT, hdr_margin=B_HDR_MARGIN, key_end=B_KEY_END, right=B_X0 + B_W),
            c=dict(x0=C_X0, gut=C_GUT, lab_w=C_LAB_W, rowhdr_w=C_ROWHDR_W, cell=C_CELL, pair_gap=C_PAIR_GAP, group_gap=C_GROUP_GAP, cell_x0=C_CELL_X0, right=C_RIGHT, w=C_W, cell_txt=C_CELL_TXT, margin=C_MARGIN, sub_txt=C_SUB_TXT, sub_margin=C_SUB_MARGIN,
                   bar0=C_LEG_BAR0, bar1=C_LEG_BAR1, ref_x=C_REF_X, leg_end=C_LEG_END), bc_gap=BC_GAP)

if __name__ == "__main__":
    print(f"sizes: {SZ}")
    print(f"panel a: HR string {HR_W:.1f}, q (bold '<0.001') {Q_W:.1f}, after-box {AFTER_BOX:.1f}; left gutter {L_LAB:.1f} ({L_LONGEST}), right gutter {R_LAB:.1f} ({R_LONGEST})")
    print(f"  fixed width {fixed:.1f} -> boxes {BOX_L} (left) and {BOX_R} (right), equal scale; left box {L_BOX0}-{L_BOX1:.1f}, HR right {L_HR_R:.1f}, q right {L_Q_R:.1f}; right labels from {R_LAB_X0}, box {R_BOX0}-{R_BOX1:.1f}, HR right {R_HR_R:.1f}, q right {R_Q_R:.1f} (sheet {W}, margin {W - R_Q_R:.2f})")
    print(f"  pt per ln unit: left [{XLO}, {XHI_L}] {PT_LN_L:.1f}, right [{XLO}, {XHI_R}] {PT_LN_R:.1f}, V30 {PT_LN_V30:.1f}; two equal boxes over [{XLO}, {XHI_L}] would give {PT_LN_SHARED:.1f} in both; clear space between neighbouring tick labels left {TICK_GAP_L:.1f}, right {TICK_GAP_R:.1f} pt")
    print(f"  marker diameter {MARK_D:.2f} pt; narrowest T90 interval left {narrowest(LEFT, PT_LN_L)}, right per-column {narrowest(RIGHT, PT_LN_R)}, right shared {narrowest(RIGHT, PT_LN_SHARED)}, V30 {narrowest(ROWS, PT_LN_V30)}")
    print(f"  data ranges: left {min(L_VALS):.3f}-{max(L_VALS):.3f}, right {min(R_VALS):.3f}-{max(R_VALS):.3f}")
    print(f"panel b: gutter {B_GUT} (labels {B_LAB_W:.1f}), cell {B_CELL} x {B_ROW_H} (text {B_CELL_TXT:.1f}, margin {B_MARGIN:.2f}; header text {B_HDR_TXT:.1f}, margin {B_HDR_MARGIN:.2f}), width {B_W:.1f}, block {B_X0}-{B_X0 + B_W:.1f}, key ends {B_KEY_END:.1f}")
    print(f"panel c: x0 {C_X0}, gutter {C_GUT} (labels {C_LAB_W:.1f}, row headers {C_ROWHDR_W:.1f}), cell {C_CELL} x {C_ROW_H} (text {C_CELL_TXT:.1f}, margin {C_MARGIN:.2f}; sub-header {C_SUB_TXT:.1f}, margin {C_SUB_MARGIN:.2f}), cells {C_CELL_X0:.1f}-{C_RIGHT:.1f}, legend ends {C_LEG_END:.1f}")
    print(f"vertical: title {TITLE_BASE}; a letter {A_TOP}/{A_BASE}; key {KEY_BASE}; headers {HDR_BASE}; rows {ROW0} + {PITCH} x 25 = {LAST_C}; axis {AX_TOP:.2f}-{AX_BOT:.2f}; ticks {TICK_BASE}; caption {CAP_BASE}; a bottom {A_BOTTOM:.1f}")
    print(f"  b/c letters {B_TOP}/{B_BASE}; b title {B_TITLE_BASE:.1f}, headers {B_HDR_Y0:.1f}+{B_HDR_H}, rows from {B_ROW_Y0:.1f} pitch {B_ROW_PITCH}, key {B_KEY_BASE:.1f}, b bottom {B_BOTTOM:.1f}")
    print(f"  c headers {C_HDR1_Y0:.1f}+{C_HDR1_H}, {C_HDR2_Y0:.1f}+{C_HDR2_H}; rows from {C_ROW_Y0:.1f} pitch {C_ROW_PITCH} (cell {C_ROW_H}, two 13-pt lines {C_LINE_DY} apart: ink {0.716 * 13 + C_LINE_DY + 0.212 * 13:.1f} in {C_ROW_H}), rows end {C_ROWS_END:.1f}; legend title {C_LEG_TITLE_BASE:.1f}, bar {C_BAR_Y0:.1f}, bar ticks {C_BAR_TICK_BASE:.1f}, asterisks {C_AST_BASE:.1f}; c bottom {C_BOTTOM:.1f}")
    print(f"  d letter {D_TOP}/{D_BASE}; band top {BAND_TOP}; band r49 height {BAND49_H:.3f} -> sheet height {H_R49} (aim 1090, gate 1260, V30 1043.821)")
    json.dump(PLAN, open(f"{C.LANE}/Main_Fig3/work/plan_r54.json", "w"), indent=1, default=float)
    print(f"wrote {C.LANE}/Main_Fig3/work/plan_r54.json")
