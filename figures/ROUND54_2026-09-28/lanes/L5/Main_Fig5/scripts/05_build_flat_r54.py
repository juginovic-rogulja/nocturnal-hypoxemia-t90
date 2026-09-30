#!$T90_PY
"""Round 54, lane L5 (2026-09-29): the round-49 flat builder of panels b and c (05_build_flat_r49.py, itself the round-42 LFIG5b builder) rerun
at the ROUND-54 SIZES with a re-solved horizontal layout. Sizes (R54_BRIEF.md, by role through l5a14.ROLE_PT): row labels, key entries,
tick labels and the HR (95% CI), P and q strings 14 pt; the column heads, the axis titles and the 'Negative controls' head 15 pt; the panel
letters 18 pt. Every string of V30 stays (same strings). What the larger text forced (every change declared in the build record and in
checks_flat_build.txt):
  panel b: the label column and the letter b move 9 pt right (B_X 497.45 -> 506.45) so that panel a's arm text keeps 20 pt of air; the
    log axis keeps its value range (0.298 to 10.55) and its tick values but is compressed (b pt per ln unit solved as the largest value at
    which every label clears its own whisker and the reference line by 2 pt, the widest whisker clears the HR column by 4 pt, the HR and
    P columns keep 6 pt between them and the P column ends at x 940, the sheet's right margin 29 pt instead of the V13 35.7); the rows
    keep their V13 y (pitch 24.670), the two lines of the one two-line label sit 8.5 pt above and below the row (11.6 pt apart before);
    the key's two entries 16.5 pt apart (10.95 before) and the dashed reference line starts below the key (its top 92 instead of 50.84)
    because the 14-pt key would cross it; markers at the cap centre of the 14-pt label (5.0 pt above the baseline, 2.44 before);
  panel c: the log axis keeps its range (0.949 to 3.306) and ticks, compressed the same way with the reference line 4 pt right of the
    widest label and the q column ending 23 pt left of panel b's labels; rows at their V13 y (pitch 21.354);
  both: tick labels 3 pt lower (baseline 690.3, 14-pt caps would touch the tick ends at the V13 baseline), the axis titles 19 pt under
    the tick labels (the two lines of panel b's title 17.4 pt apart), everything else at its V13 position.
Output: work/panels_bc_flat.pdf, work/Main_Fig5_bc_drawn.json, work/Main_Fig5_rows.json, work/flat_build_record.json, verify/checks_flat_build.txt.
Round-49 and round-42 docstrings: the V13 column edges were reconstructed at the V13 sizes and the columns grew leftward (round 49, +1 pt);
panel b keeps ONE contrast per condition, 'Still above 10% on PAP' against the reference (restored to 1% or less), one single-line row per
condition (the 18 largest adj_pub_hr, then the three controls), the upper-line series only (adj_B_hr/lo/hi, P = adj_B_p, bold and the
orange dot filled where P < 0.05), heads 'HR (95% CI)' and 'P', a two-entry key on the left over the labels; the 21 rows spread evenly
between the V13 first-row y and the V13 last-row y at ONE pitch, the 'Negative controls' head in its V13 block; panel c unchanged
(results_continuous.csv, step 147); panel a's seven counts re-read from treatment_v2.json here."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from l5a14 import *
from fontTools.ttLib import TTFont as _TTFont

SHEET = "Main_Fig5"; D = f"{paths.FIGURE_ROOT}/ROUND54_2026-09-28/lanes/L5/Main_Fig5"; W = f"{D}/work"; V = f"{D}/verify"; os.makedirs(f"{V}/crops", exist_ok=True)
BASE = f"{V13}/{SHEET}.pdf"; FLAT = f"{W}/panels_bc_flat.pdf"; N_TOP = 18
PT = ROLE_PT; assert Sheet.BUMP == 0.0
ck = Checks(f"Lane L5 round 54 (labels, values, ticks, keys 14 pt; heads, axis titles 15; letters 18), the round-49/42 flat build of panels b and c re-laid, {SHEET} flat build, {time.strftime('%Y-%m-%d %H:%M')}. BASE (design) = V13 {BASE} (sha256 {sha256(BASE)[:16]}). NEW = L4_fig5b_refit_v3.json (212, adjusted, BH q) for panel b: the {N_TOP} largest adj_pub_hr + 3 controls, HR (95% CI) and P (uncorrected, adj_B_p) columns; results_continuous.csv (147) for panel c; treatment_v2.json (117) for the panel a counts.")
SC = 0.95   # the panels b and c were placed on the sheet at 0.95 of the generators' 171 mm canvases (round 12 recompose)

# ---- round-54 layout constants (declared)
DX_B = 9.0            # panel b's label column, key and letter 9 pt right of V13 (497.45 -> 506.45): panel a's 14-pt arm text ends at 485.6
RIGHT_B = 940.0       # panel b's P column right edge (V13 933.20): the sheet's right margin 28.9 pt
GAP_CB = 20.0         # air between panel c's q column and panel b's labels (V13 25.7)
AIR_LAB_WHISK = 5.0   # a row label ends at least this far left of its own whisker start (round 42: 2; the whisker's round cap reaches 0.85 pt further left, so 4.15 pt of visible air)
CLIP_AC = 385.0       # the compose clips panel a above and panel c below this y (V13 393): panel a's lowest text ends at 364.5, the 18-pt letter c reaches up to 393.6
AIR_LAB_REF_B = 2.0   # and of the dashed reference line (panel b, the round-42 rule)
AIR_LAB_REF_C = 6.0   # panel c: the reference line is the axis's left boundary, 6 pt of air
AIR_WHISK_COL = 4.0   # the widest whisker end to the HR column (round 42: 2)
COL_GAP = 6.0         # HR column right edge to the P or q column left edge (V13: 5.3 and 8.2)
OFF_MARK = -0.716 * PT["label"] / 2      # markers at the cap centre of the 14-pt row label: 5.01 pt above the baseline (V13: 2.44 in b, 3.40 in c)
NUDGE_2L = 8.5        # the two lines of the one two-line label at the row y minus and plus 8.5 pt (17 pt apart)
KEY_DY = 16.5         # the key's second entry 16.5 pt under the first (V13 10.95 at 9.5 pt)
TL_AIR = 3.5          # tick label cap top this far under the tick ends
TITLE_GAP = 19.0      # axis title baseline this far under the tick label baseline
TITLE_DY = 12.10 / 10.45 * PT["axis_title"]   # the two lines of panel b's title: the V13 ratio (12.10 pt at 10.45) at 15 pt = 17.37

# ---- 0 gates
sc_t = gate("treatment_v2"); sc_t3 = gate("treatment_v3"); sc_r = gate("fig5b_refit_v3_json"); sc_c = gate("results_continuous")
dump({"treatment_v2": sc_t, "treatment_v3": sc_t3, "fig5b_refit_v3_json": sc_r, "results_continuous": sc_c}, f"{W}/sidecars.json")
ck.log("treatment_v2.json (117), treatment_v3.json (211), L4_fig5b_refit_v3.json (212, adjusted, --bh) and results_continuous.csv (147) pass the v8.1 sidecar gate and the v8.2 mtime gate", all(s["v8_inputs"] for s in (sc_t, sc_t3, sc_r, sc_c)), f"out {sc_t['output_mtime']}, {sc_t3['output_mtime']}, {sc_r['output_mtime']}, {sc_c['output_mtime']}")

# ---- 1 the V13 text layer (cached probe) and the V13 geometry (the round-37 block, unchanged: the design's row positions and axis ranges)
J = load_json(f"{W}/base_text.json"); spans, info = J["spans"], J["info"]; PW, PH = info["rect"][2], info["rect"][3]
ck.log("V13 page 968.94 x 946.23, nested (454 XObjects, no get_drawings)", abs(PW - 968.94) < 0.05 and abs(PH - 946.2253) < 0.05, info)
letters = {s["text"]: s for s in spans if len(s["text"]) == 1 and s["text"] in "abcd" and s["size"] > 12.5}; assert set(letters) == set("abcd")
title = find_span(spans, "Figure 5"); assert title
labB = [s for s in column(spans, 497, 498, ymin=100, ymax=700, size=10.45) if "Bold" not in s["font"]]; headB = find_span(spans, "Negative controls", xmin=490); assert len(labB) == 22 and headB
B_Y = [s["y"] for s in labB]; B_PITCH = float(np.mean(np.diff(B_Y[:17]))); B_X13 = float(np.mean([s["x"] for s in labB])); B_X = B_X13 + DX_B
hrB = column(spans, 830, 836, ymin=100, ymax=700); pB = column(spans, 900, 926, ymin=100, ymax=700); assert len(hrB) == 44 and len(pB) == 44, (len(hrB), len(pB))
S0 = Sheet(10, 10)   # font metrics only
HR_R_B13 = float(np.median([s["x"] + S0.width(s["text"], s["size"]) for s in hrB])); P_R_B13 = float(np.median([s["x"] + S0.width(s["text"], s["size"], "Bold" in s["font"]) for s in pB]))
assert np.std([s["x"] + S0.width(s["text"], s["size"]) for s in hrB]) < 0.25 and np.std([s["x"] + S0.width(s["text"], s["size"], "Bold" in s["font"]) for s in pB]) < 0.2
keyB = sorted([s for s in spans if abs(s["x"] - 694.31) < 0.5], key=lambda s: s["y"]); assert [s["text"] for s in keyB] == ["Still above 10% on PAP", "1 to 10% on PAP", "Reference: 1% or less on PAP"]
colhB = [find_span(spans, "HR (95% CI)", xmin=800, ymax=80), find_span(spans, "P", xmin=900, ymax=80)]; assert all(colhB)
tlB = sorted([s for s in spans if abs(s["y"] - 687.49) < 0.5 and s["x"] > 600], key=lambda s: s["x"]); assert [s["text"] for s in tlB] == ["0.5", "1", "2", "4", "8"]
titleB = sorted([s for s in spans if s["x"] > 600 and 700 < s["y"] < 720], key=lambda s: s["y"]); assert [s["text"] for s in titleB] == ["Hazard ratio for a new diagnosis", "vs 1% or less on PAP (95% CI)"]
SPINE_B = dict(y=673.85, x0=616.92, x1=824.64, w=0.8 * SC); TICKS_B = [646.98, 687.42, 727.74, 768.18, 808.5]; TICK_LEN = 3.0 * SC; AXTOP_B13 = 50.84; REF_X_B13 = 687.35
SPINE_C = dict(y=673.85, x0=155.52, x1=363.36, w=0.8 * SC); TICKS_C = [164.34, 231.78, 279.66, 347.22]; AXTOP_C = 432.4; REF_X_C13 = 164.3
TICKV_B = [0.5, 1, 2, 4, 8]; TICKV_C = [1, 1.5, 2, 3]
aB13, bB13, rB = fit_axis(TICKV_B, TICKS_B); aC13, bC13, rC = fit_axis(TICKV_C, TICKS_C); assert max(np.abs(rB).max(), np.abs(rC).max()) < 0.15, (rB, rC)
LD = load_json(f"{paths.FIGURE_ROOT}/ROUND30_2026-09-08/LD_fig5/work/markers_geometry.json"); DOT_R = float(np.median([o["page_radius"] for o in LD["markers"]]))
ck.log(f"V13 panel b geometry: 22 rows pitch {B_PITCH:.3f} at x {B_X13:.2f}, HR column right edge {HR_R_B13:.2f}, P/q column right edge {P_R_B13:.2f}; spine y {SPINE_B['y']} x {SPINE_B['x0']} to {SPINE_B['x1']}; log axis a {aB13:.2f} b {bB13:.3f} pt per ln unit (residual {np.abs(rB).max():.3f}); dot radius {DOT_R} (round-30 record)", abs(aB13 - REF_X_B13) < 0.3, "")
labC = column(spans, 36, 37, ymin=430, ymax=680, size=9.5); assert len(labC) == 11; C_Y = [s["y"] for s in labC]; C_PITCH = float(np.mean(np.diff(C_Y))); C_X = float(np.mean([s["x"] for s in labC]))
hrC = column(spans, 370, 373, ymin=430, ymax=680); qC = column(spans, 440, 455, ymin=430, ymax=680); assert len(hrC) == 11 and len(qC) == 11
HR_R_C13 = float(np.mean([s["x"] + S0.width(s["text"], s["size"]) for s in hrC])); Q_R_C13 = float(np.mean([s["x"] + S0.width(s["text"], s["size"], True) for s in qC]))
OFF_C = float(np.mean([h["y"] - s["y"] for h, s in zip(hrC, labC)]))
colhC = [find_span(spans, "HR (95% CI)", xmin=380, xmax=390, ymin=415, ymax=430), find_span(spans, "q", xmin=460, ymin=415, ymax=430)]; assert all(colhC)
tlC = sorted([s for s in spans if abs(s["y"] - 686.82) < 0.5 and s["x"] < 400], key=lambda s: s["x"]); assert [s["text"] for s in tlC] == ["1", "1.5", "2", "3"]
titleC = find_span(spans, "Hazard ratio per 1 SD of residual sleep T90 (95% CI)"); assert titleC
ck.log(f"V13 panel c geometry: 11 rows pitch {C_PITCH:.3f} at x {C_X:.2f}, columns right edges {HR_R_C13:.2f} and {Q_R_C13:.2f}, strings at {OFF_C:+.2f}; spine y {SPINE_C['y']}; log axis a {aC13:.2f} b {bC13:.3f} (residual {np.abs(rC).max():.3f})", abs(aC13 - REF_X_C13) < 0.3, "")
C_DOT, C_EDGE, C_CI = 5.29 * SC, 0.8 * SC, 1.7 * SC; B_DOT, B_EDGE_F, B_EDGE_O, B_CI, B_SQ = 2 * DOT_R, 0.8 * SC, 1.1 * SC, 1.8 * SC, 6.0 * SC
KEY_HANDLE = 0.30 * 72 * SC; KEY_GAP = 0.10 * 72 * SC                        # the builder's 0.30 in handle ending 0.10 in before the label
KEY_H0, KEY_H1, KEY_LABEL_X = B_X, B_X + KEY_HANDLE, B_X + KEY_HANDLE + KEY_GAP   # round 42: the key block left-aligned at the labels' x
KEY_Y = [keyB[0]["y"], keyB[0]["y"] + KEY_DY]
HEAD_Y_B = colhB[0]["y"]; HEAD_Y_C = float(np.mean([s["y"] for s in colhC]))

# ---- 2 panel a: the seven counts re-read from the v8.1 file must equal the V13 strings (the panel is rebuilt by 04_build_panel_a_r54.py)
T, rows = treatment_rows(); mig = T["migration"]; fw = mig["from_worst"]; cvn = T["corrected_vs_not"]
RULES.update({f"band_n:{lab}": (lambda l: (lambda v: f"{l}, n = {int_comma(v)}"))(lab) for lab in ("0–1%", "1–5%", "5–10%", ">10%")})
RULES["label_line"] = lambda v: str(v)      # a two-line label: the drawn text is one line of the file label (checked as a substring by the verifier)
A_REGION = dict(xmin=0, xmax=497, ymin=16, ymax=380)
pa = [("pre block", f"n = {fw['n']:,}", "migration/from_worst/n", fw["n"], "n_eq"),
      ("on-PAP band 0-1%", f"0–1%, n = {fw['to']['0-1%']:,}", "migration/from_worst/to/0-1%", fw["to"]["0-1%"], "band_n:0–1%"),
      ("on-PAP band 1-5%", f"1–5%, n = {fw['to']['1-5%']:,}", "migration/from_worst/to/1-5%", fw["to"]["1-5%"], "band_n:1–5%"),
      ("on-PAP band 5-10%", f"5–10%, n = {fw['to']['5-10%']:,}", "migration/from_worst/to/5-10%", fw["to"]["5-10%"], "band_n:5–10%"),
      ("on-PAP band >10%", f">10%, n = {fw['to']['>10%']:,}", "migration/from_worst/to/>10%", fw["to"][">10%"], "band_n:>10%"),
      ("restored arm", f"n = {cvn['n_corrected']:,} ({round(100 * cvn['n_corrected'] / cvn['n'])}%)", "corrected_vs_not/n_corrected, n", [cvn["n_corrected"], cvn["n"]], "count_pct"),
      ("still-above arm", f"n = {cvn['n_not']:,} ({round(100 * cvn['n_not'] / cvn['n'])}%)", "corrected_vs_not/n_not, n", [cvn["n_not"], cvn["n"]], "count_pct")]
a_found = [find_span(spans, exp, **A_REGION) for _, exp, _, _, _ in pa]
ck.log("panel a: all 7 count strings of the V13 region equal treatment_v2.json (1,894; 613; 498; 215; 568; 1,326 (70%); 568 (30%)): the panel is rebuilt from its generator at 14/15 pt with this read-back as its proof", all(a_found), f"{[exp for _, exp, _, _, _ in pa]}")
A_DRAWN = [dict(text=exp, x=g["x"], baseline=g["y"], ha="left", size=g["size"], panel="a", source_file=SPEC["treatment_v2"], source_key=key, source_value=val, rule=rule, note="V13 probe position; the round-54 panel a page carries the string at 14 pt (checked by the verifier on the composed sheet)") for (lab, exp, key, val, rule), g in zip(pa, a_found) if g]
# the round-54 panel a page (flat, small: PyMuPDF allowed) mapped onto the sheet: its arm text must clear panel b's label column
MAP = load_json(f"{W}/panel_a/map_v13.json"); S_A, TX_A, TY_A = MAP["S"], MAP["TX"], MAP["TY"]; PA_PDF = f"{W}/panel_a/Figure5A_r54.pdf"
pa_sp, pa_info = text_spans(PA_PDF); pa_sheet = [dict(text=s["text"], x0=S_A * s["bbox"][0] + TX_A, y0=S_A * s["bbox"][1] + TY_A, x1=S_A * s["bbox"][2] + TX_A, y1=S_A * s["bbox"][3] + TY_A, size=S_A * s["size"]) for s in pa_sp]
PA_RIGHT = max(s["x1"] for s in pa_sheet); PA_BOTTOM = max(s["y1"] for s in pa_sheet); PA_TOP = min(s["y0"] for s in pa_sheet)
ck.log(f"panel a (round-54 page, mapped by the V13 similarity S {S_A:.5f}): {len(pa_sheet)} strings at {sorted(set(round(s['size'], 2) for s in pa_sheet))} pt on the sheet, the rightmost ends at x {PA_RIGHT:.1f} (arm text at X_ARM 76.5), at least 20 pt left of panel b's label column {B_X:.2f}; the lowest ends at y {PA_BOTTOM:.1f}, above the a/c clip boundary {CLIP_AC} minus 4; the highest starts at y {PA_TOP:.1f}", PA_RIGHT <= B_X - 20 and PA_BOTTOM < CLIP_AC - 4 and PA_TOP > 16 + 4 and sorted(set(round(s["size"], 2) for s in pa_sheet)) == [14.0, 15.0])

# ---- 3 panel b rows (decision 1) and panel c rows (unchanged rules)
R = load_json(SPEC["fig5b_refit_v3_json"]); OC = R["outcomes"]; C = R["counts"]
assert "--bh" in R["built"] and "adjust" in R["built"], R["built"]
EV_MIN = min(v["events"] for v in OC.values()); N_FILE = len(OC)
def view(k, v):
    w = dict(v, outcome_label=k)
    for band in ("A", "B"):
        for f in ("hr", "lo", "hi", "p", "q"): w[f"{band}_{f}"] = v[f"adj_{band}_{f}"]
    w["pub_hr"] = v["adj_pub_hr"]; return w
NC_ALL = sorted(((k, view(k, v)) for k, v in OC.items() if not v["negative_control"]), key=lambda kv: -kv[1]["pub_hr"])
SIG = NC_ALL[:N_TOP]; CTRL = sorted(((k, view(k, v)) for k, v in OC.items() if v["negative_control"]), key=lambda kv: -kv[1]["pub_hr"])
assert len(NC_ALL) > N_TOP and SIG[-1][1]["pub_hr"] > NC_ALL[N_TOP][1]["pub_hr"], "no strict cut at the 18th ratio"
assert [k for k, _ in CTRL] == ["Hemorrhoids", "Glaucoma", "Contact dermatitis"], [k for k, _ in CTRL]
TIES_3DP = [(SIG[i][0], SIG[i + 1][0], SIG[i][1]["pub_hr"], SIG[i + 1][1]["pub_hr"]) for i in range(N_TOP - 1) if halfup(SIG[i][1]["pub_hr"], 3) == halfup(SIG[i + 1][1]["pub_hr"], 3)]
ck.log(f"panel b rows (decision 1): the {N_TOP} non-control conditions with the largest adjusted two-group hazard ratio adj_pub_hr of the refit file ({N_FILE} outcomes, {len(NC_ALL)} non-controls, every outcome with at least 40 events: minimum {EV_MIN}), descending at full precision (18th {SIG[-1][0]} {SIG[-1][1]['pub_hr']:.4f} > 19th {NC_ALL[N_TOP][0]} {NC_ALL[N_TOP][1]['pub_hr']:.4f}), then the 3 controls by the same ratio {[k for k, _ in CTRL]}; consecutive ratios equal at 3 dp (order decided at full precision): {TIES_3DP}",
       EV_MIN >= 40 and len(SIG) == N_TOP and len(CTRL) == 3 and C["n_all"] == cvn["n"] and all(g["match_3dp"] for g in R["positive_control"]), f"rows {[k for k, _ in SIG]}")
n_sigB = sum(v["B_p"] < 0.05 for _, v in SIG)
ck.log(f"panel b: the printed column is the uncorrected P of the contrast (adj_B_p of the adjusted refit); P < 0.05 on the 18 rows: {n_sigB}; controls with P < 0.05: {[k for k, v in CTRL if v['B_p'] < 0.05]}", all(0 <= v["B_p"] <= 1 for _, v in SIG + CTRL), "")
cont, Cc = residual_rows()
ck.log(f"panel c (unchanged): {len(Cc)} Benjamini-Hochberg survivors of results_continuous.csv, all q < 0.05 and lower CI above 1; at the V13 pitch the last row sits above the axis", (Cc.q_sd < 0.05).all() and (Cc.lo_sd > 1.0).all() and C_Y[0] + (len(Cc) - 1) * C_PITCH + 6 < SPINE_C["y"], list(Cc.outcome))
xloB, xhiB = exp((SPINE_B["x0"] - aB13) / bB13), exp((SPINE_B["x1"] - aB13) / bB13); xloC, xhiC = exp((SPINE_C["x0"] - aC13) / bC13), exp((SPINE_C["x1"] - aC13) / bC13)
ck.log("every drawn value (the B contrast) lies inside the V13 axis value ranges (kept: no tick change)", all(xloB <= v["B_lo"] and v["B_hi"] <= xhiB for _, v in SIG + CTRL) and all(xloC <= r.lo_sd and r.hi_sd <= xhiC for r in Cc.itertuples()), f"b {xloB:.3f}-{xhiB:.3f}, c {xloC:.3f}-{xhiC:.3f}")
NS, NCn = len(SIG), len(CTRL)
N_PITCH = (NS - 1) + 3 + (NCn - 1)
P_NEW = (B_Y[-1] - B_Y[0]) / N_PITCH
yS = [B_Y[0] + i * P_NEW for i in range(NS)]; y_head = B_Y[0] + (NS - 1 + 1.5) * P_NEW; yC = [B_Y[0] + (NS - 1 + 3 + j) * P_NEW for j in range(NCn)]
assert abs(yS[0] - B_Y[0]) < 1e-9 and abs(yC[-1] - B_Y[-1]) < 1e-9, (yS[0], yC[-1])
ck.log(f"panel b rows: {NS} single-line rows then the head then {NCn} controls at ONE pitch {P_NEW:.3f} pt (V13 23.598) from the V13 first row {B_Y[0]:.2f} to the V13 last row {B_Y[-1]:.2f}, the head at {y_head:.2f}; the spine at the V13 y {SPINE_B['y']} = panel c's; at 14 pt the pitch leaves {P_NEW - 1.117 * PT['label']:.1f} pt between the ascender of one row and the descender of the next", P_NEW > 1.3 * PT["label"] and SPINE_B["y"] == SPINE_C["y"], f"controls {[round(y, 2) for y in yC]}")

# ---- 3b the round-54 horizontal layout: labels at 14 pt, the axes compressed to the largest scale that keeps every clearance
def label_lines(k):
    if " or " in k and S0.width(k, PT["label"]) > 150: a, b = k.split(" or ", 1); return [a, "or " + b]   # the Figure 4b idiom: the one long label on two lines
    return [k]
LAB = {k: label_lines(k) for k, _ in SIG + CTRL}; LAB_END = {k: max(B_X + S0.width(t, PT["label"]) for t in LAB[k]) for k in LAB}
two = [k for k in LAB if len(LAB[k]) == 2]; assert two == ["Ventricular arrhythmia or cardiac arrest"], two
W_HR_B = max(S0.width(hr_ci(v["B_hr"], v["B_lo"], v["B_hi"]), PT["value"]) for _, v in SIG + CTRL); W_P_B = max(S0.width(fmt_p(v["B_p"]), PT["value"], v["B_p"] < 0.05) for _, v in SIG + CTRL)
W_HR_C = max(S0.width(hr_ci(r.hr_sd, r.lo_sd, r.hi_sd), PT["value"]) for r in Cc.itertuples()); W_Q_C = max(S0.width(fmt_p(r.q_sd), PT["value"], True) for r in Cc.itertuples())
def solve_axis(label_ends, los, his, air_ref, right, b13):
    """The largest b (pt per ln unit, at most the V13 value) at which: ref >= every label end + AIR_LAB_WHISK - b ln(lo) (own whisker) and
    ref >= max label end + air_ref (the reference line); the column chain ref + b ln(max hi) + AIR_WHISK_COL + W_HR + COL_GAP + W_P ends at
    'right'. ref is then the smallest value allowed. Returns (b, ref, binding constraint)."""
    def ref_of(b):
        c1 = [(le + AIR_LAB_WHISK - b * log(lo), f"own whisker of '{k}'") for k, le, lo in zip(label_ends, [label_ends[k] for k in label_ends], los)]
        c2 = [(max(label_ends.values()) + air_ref, f"reference line after the widest label '{max(label_ends, key=label_ends.get)}'")]
        return max(c1 + c2, key=lambda t: t[0])
    def right_of(b): return ref_of(b)[0] + b * log(max(his)) + AIR_WHISK_COL + (W_HR_B if right == RIGHT_B else W_HR_C) + COL_GAP + (W_P_B if right == RIGHT_B else W_Q_C)
    if right_of(b13) <= right: b = b13
    else:
        lo_b, hi_b = 1.0, b13
        for _ in range(80):
            mid = (lo_b + hi_b) / 2
            if right_of(mid) <= right: lo_b = mid
            else: hi_b = mid
        b = lo_b
    r, why = ref_of(b); return b, r, why
losB = [v["B_lo"] for k, v in SIG + CTRL]; hisB = [v["B_hi"] for k, v in SIG + CTRL]
bB, REF_X_B, whyB = solve_axis(LAB_END, losB, hisB, AIR_LAB_REF_B, RIGHT_B, bB13); aB = REF_X_B
X0_B, X1_B = aB + bB * log(xloB), aB + bB * log(xhiB); TKX_B = [aB + bB * log(t) for t in TICKV_B]
WH_END_B = max(aB + bB * log(h) for h in hisB); HR_R_B = WH_END_B + AIR_WHISK_COL + W_HR_B; P_R_B = HR_R_B + COL_GAP + W_P_B
ck.log(f"panel b axis (round 54): the value range {xloB:.3f} to {xhiB:.3f} and the ticks {TICKV_B} kept; scale b {bB:.3f} pt per ln unit = {100 * bB / bB13:.1f} percent of the V13 {bB13:.3f} (the largest at which every 14-pt label clears its own whisker by {AIR_LAB_WHISK} pt and the reference line by {AIR_LAB_REF_B} pt, the widest whisker end {WH_END_B:.2f} clears the HR column by {AIR_WHISK_COL} pt, the HR and P columns {COL_GAP} pt apart, the P column ending at {RIGHT_B}); binding: {whyB}; reference line x {REF_X_B:.2f} (V13 {REF_X_B13}); spine {X0_B:.2f} to {X1_B:.2f} (V13 {SPINE_B['x0']} to {SPINE_B['x1']}); HR column right edge {HR_R_B:.2f} (V13 {HR_R_B13:.2f}), P column right edge {P_R_B:.2f} (V13 {P_R_B13:.2f})",
       bB <= bB13 + 1e-9 and abs(P_R_B - RIGHT_B) < 0.01 and bB > 0.75 * bB13, f"widths: HR {W_HR_B:.1f}, P {W_P_B:.1f}")
LAB_END_C = {r.outcome: C_X + S0.width(plab(r.outcome), PT["label"]) for r in Cc.itertuples()}
RIGHT_C = B_X - GAP_CB
bC, REF_X_C, whyC = solve_axis(LAB_END_C, list(Cc.lo_sd), list(Cc.hi_sd), AIR_LAB_REF_C, RIGHT_C, bC13); aC = REF_X_C
X0_C, X1_C = aC + bC * log(xloC), aC + bC * log(xhiC); TKX_C = [aC + bC * log(t) for t in TICKV_C]
WH_END_C = max(aC + bC * log(h) for h in Cc.hi_sd); HR_R_C = WH_END_C + AIR_WHISK_COL + W_HR_C; Q_R_C = HR_R_C + COL_GAP + W_Q_C
ck.log(f"panel c axis (round 54): the value range {xloC:.3f} to {xhiC:.3f} and the ticks {TICKV_C} kept; scale b {bC:.3f} pt per ln unit = {100 * bC / bC13:.1f} percent of the V13 {bC13:.3f} (every 14-pt label clears the reference line by {AIR_LAB_REF_C} pt and its own whisker by {AIR_LAB_WHISK}, the widest whisker end {WH_END_C:.2f} clears the HR column by {AIR_WHISK_COL}, the q column ends at {RIGHT_C:.2f} = panel b's labels {B_X:.2f} minus {GAP_CB}); binding: {whyC}; reference line x {REF_X_C:.2f} (V13 {REF_X_C13}); spine {X0_C:.2f} to {X1_C:.2f} (V13 {SPINE_C['x0']} to {SPINE_C['x1']}); HR column right edge {HR_R_C:.2f} (V13 {HR_R_C13:.2f}), q column right edge {Q_R_C:.2f} (V13 {Q_R_C13:.2f})",
       bC <= bC13 + 1e-9 and abs(Q_R_C - RIGHT_C) < 0.01 and bC > 0.6 * bC13, f"widths: HR {W_HR_C:.1f}, q {W_Q_C:.1f}")
ROW_LEFT = {k: axis_x(aB, bB, v["B_lo"]) for k, v in SIG + CTRL}
ck.log(f"panel b labels at {PT['label']} pt: the widest line ends at x {max(LAB_END.values()):.1f} ('{max(LAB_END, key=LAB_END.get)}'); every label ends at least {AIR_LAB_WHISK} pt left of its own row's whisker start and {AIR_LAB_REF_B} pt left of the reference line {REF_X_B:.2f}; {len(two)} label on two lines {[(k, LAB[k]) for k in two]}", all(LAB_END[k] <= ROW_LEFT[k] - AIR_LAB_WHISK + 1e-6 and LAB_END[k] <= REF_X_B - AIR_LAB_REF_B + 1e-6 for k in LAB), f"tightest: {min(((round(ROW_LEFT[k] - LAB_END[k], 2), k) for k in LAB))}")
AXTOP_B = KEY_Y[1] + 0.212 * PT["key"] + 6.0     # the dashed reference line starts 6 pt under the key's second line (its 14-pt labels cross the line's x)
TL_Y = SPINE_B["y"] + TICK_LEN + TL_AIR + 0.716 * PT["tick"]
TITLE_Y_C = TL_Y + TITLE_GAP; TITLE_Y_B = [TL_Y + TITLE_GAP, TL_Y + TITLE_GAP + TITLE_DY]
CONTENT_BOTTOM = TITLE_Y_B[1] + 0.212 * PT["axis_title"]
ck.log(f"vertical (round 54): key entries at baselines {[round(y, 2) for y in KEY_Y]} ({KEY_DY} pt apart), the dashed reference line of panel b from y {AXTOP_B:.2f} (under the key; V13 {AXTOP_B13}) to the spine; tick labels at baseline {TL_Y:.2f} (V13 687.49 / 686.82: cap top {TL_AIR} pt under the tick ends); axis titles at {TITLE_Y_C:.2f} (c) and {[round(y, 2) for y in TITLE_Y_B]} (b); the flat content ends at y {CONTENT_BOTTOM:.2f} (the compose puts the band 20 pt lower)", AXTOP_B < B_Y[0] - 20 and TL_Y - 0.716 * PT["tick"] - (SPINE_B["y"] + TICK_LEN) >= TL_AIR - 1e-9)

# ---- 4 draw panels b and c on the flat page
S = Sheet(PW, PH)
S.text(B_X, letters["b"]["y"], "b", PT["letter"], bold=True, color="#" + letters["b"]["color"], panel="b", source_key="panel letter (x = the label column, 9 pt right of V13)", role="letter"); S.text(letters["c"]["x"], letters["c"]["y"], "c", PT["letter"], bold=True, color="#" + letters["c"]["color"], panel="c", source_key="panel letter", role="letter")
S.line(REF_X_B, AXTOP_B, REF_X_B, SPINE_B["y"], INK20, 0.9 * SC, cap=0, dashes=f"[{4 * SC:g} {3 * SC:g}] 0")
KEY_DRAWN = []
KEY2 = [dict(keyB[0], y=KEY_Y[0]), dict(keyB[2], y=KEY_Y[1])]                  # two entries: 'Still above 10% on PAP' and 'Reference: 1% or less on PAP'
for s, kind in zip(KEY2, ("B", "ref")):
    yk = s["y"] + OFF_MARK
    if kind == "ref": S.line(KEY_H0, yk, KEY_H1, yk, INK20, 0.9 * SC, cap=0, dashes=f"[{4 * SC:g} {3 * SC:g}] 0")
    else: S.line(KEY_H0, yk, KEY_H1, yk, BLUE, B_CI, cap=1); S.circle((KEY_H0 + KEY_H1) / 2, yk, B_DOT, ORANGE, WHITE, B_EDGE_F)
    S.text(KEY_LABEL_X, s["y"], s["text"], PT["key"], color="#" + s["color"], panel="b", source_key="key label (V13 text; round 42: two entries, the key block left-aligned at the labels' x; round 54: 14 pt, 16.5 pt apart)", role="key")
    KEY_DRAWN.append(dict(text=s["text"], kind=kind, handle_x=[KEY_H0, KEY_H1], label_x=KEY_LABEL_X, baseline=s["y"], marker_y=yk, v13_label_x=s["x"], v13_baseline=keyB[0]["y"] if kind == "B" else keyB[1]["y"]))
S.text(HR_R_B, HEAD_Y_B, colhB[0]["text"], PT["head"], color="#" + colhB[0]["color"], align="right", panel="b", source_key="static V13 string (column head, right-aligned at the HR column edge)", role="head")
S.text(P_R_B, HEAD_Y_B, colhB[1]["text"], PT["head"], color="#" + colhB[1]["color"], align="right", panel="b", source_key="static V13 string (column head 'P', right-aligned at the P column edge)", role="head")
drawnB = []
for (k, v), y, is_ctrl in [(kv, yy, False) for kv, yy in zip(SIG, yS)] + [(kv, yy, True) for kv, yy in zip(CTRL, yC)]:
    series, name = "B", "Still above 10% on PAP"
    hr, lo, hi, q, p = v[f"{series}_hr"], v[f"{series}_lo"], v[f"{series}_hi"], v[f"{series}_q"], v[f"{series}_p"]; sig = p < 0.05
    ym = y + OFF_MARK
    line_c = GREY if is_ctrl else BLUE; mark_c = GREY if is_ctrl else ORANGE
    S.line(axis_x(aB, bB, lo), ym, axis_x(aB, bB, hi), ym, line_c, B_CI, cap=1)
    if is_ctrl: S.square(axis_x(aB, bB, hr), ym, B_SQ, mark_c if sig else WHITE, WHITE if sig else mark_c, B_EDGE_F if sig else B_EDGE_O)
    else: S.circle(axis_x(aB, bB, hr), ym, B_DOT, mark_c if sig else WHITE, WHITE if sig else mark_c, B_EDGE_F if sig else B_EDGE_O)
    S.text(HR_R_B, y, hr_ci(hr, lo, hi), PT["value"], color=INK20, align="right", panel="b", source_file=SPEC["fig5b_refit_v3_json"], source_key=f"outcomes/{k}/adj_{series}_hr,adj_{series}_lo,adj_{series}_hi", source_value=[hr, lo, hi], rule="hr_ci", note=name + ", " + ADJ_NOTE, role="value")
    S.text(P_R_B, y, fmt_p(p), PT["value"], bold=bool(sig), color=INK20, align="right", panel="b", source_file=SPEC["fig5b_refit_v3_json"], source_key=f"outcomes/{k}/adj_{series}_p", source_value=p, rule="fmt_p", note=name + ", uncorrected P of the adjusted contrast, bold where P < 0.05", role="value")
    drawnB.append(dict(outcome=k, label=plab(k), series=series, hr=hr, lo=lo, hi=hi, q=q, p=p, x=axis_x(aB, bB, hr), y=ym, x_lo=axis_x(aB, bB, lo), x_hi=axis_x(aB, bB, hi), filled=sig, control=is_ctrl, s_hr=hr_ci(hr, lo, hi), s_p=fmt_p(p), s_q=fmt_p(q), colour=dict(line=line_c, mark=mark_c)))
    if len(LAB[k]) == 1: S.text(B_X, y, LAB[k][0], PT["label"], color=GREY if is_ctrl else INK20, panel="b", source_file=SPEC["fig5b_refit_v3_json"], source_key=f"outcomes/{k} (row label; decision 1: among the {N_TOP} largest adj_pub_hr, or a negative control)", source_value=k, rule="label", note="", role="label")
    else:
        for t, off in zip(LAB[k], (-NUDGE_2L, NUDGE_2L)): S.text(B_X, y + off, t, PT["label"], color=GREY if is_ctrl else INK20, panel="b", source_file=SPEC["fig5b_refit_v3_json"], source_key=f"outcomes/{k} (row label on two lines, the Figure 4b idiom; decision 1)", source_value=k, rule="label_line", note=f"line of the file label '{k}'", role="label")
S.text(B_X, y_head, "Negative controls", PT["group_head"], bold=True, color="#" + headB["color"], panel="b", source_key="static V13 string (group head, moved with the rows)", role="group_head")
S.line(X0_B, SPINE_B["y"], X1_B, SPINE_B["y"], INK20, SPINE_B["w"], cap=2)
for tx, lab in zip(TKX_B, tlB): S.line(tx, SPINE_B["y"], tx, SPINE_B["y"] + TICK_LEN, INK20, SPINE_B["w"], cap=0); S.text(tx, TL_Y, lab["text"], PT["tick"], color="#" + lab["color"], align="center", panel="b", source_key="axis tick", role="tick")
AXC_B = (X0_B + X1_B) / 2
for s, yy in zip(titleB, TITLE_Y_B): S.text(AXC_B, yy, s["text"], PT["axis_title"], color="#" + s["color"], align="center", panel="b", source_key="static V13 string (axis title, centred on the axis)", role="axis_title")
# panel c
S.line(REF_X_C, AXTOP_C, REF_X_C, SPINE_C["y"], INK, 0.9 * SC, cap=0, dashes=f"[{4 * SC:g} {3 * SC:g}] 0")
S.text(HR_R_C, HEAD_Y_C, colhC[0]["text"], PT["head"], color="#" + colhC[0]["color"], align="right", panel="c", source_key="static V13 string (column head, right-aligned at the HR column edge)", role="head")
S.text(Q_R_C, HEAD_Y_C, colhC[1]["text"], PT["head"], color="#" + colhC[1]["color"], align="right", panel="c", source_key="static V13 string (column head, right-aligned at the q column edge)", role="head")
drawnC = []
for i, r in enumerate(Cc.itertuples()):
    y = C_Y[0] + i * C_PITCH; ym = y + OFF_MARK
    S.line(axis_x(aC, bC, r.lo_sd), ym, axis_x(aC, bC, r.hi_sd), ym, BLUE, C_CI, cap=1); S.circle(axis_x(aC, bC, r.hr_sd), ym, C_DOT, BLUE, WHITE, C_EDGE)
    S.text(C_X, y, plab(r.outcome), PT["label"], panel="c", source_file=SPEC["results_continuous"], source_key=f"key={r.key}:outcome", source_value=r.outcome, rule="label", role="label")
    S.text(HR_R_C, y + OFF_C, hr_ci(r.hr_sd, r.lo_sd, r.hi_sd), PT["value"], align="right", panel="c", source_file=SPEC["results_continuous"], source_key=f"key={r.key}:hr_sd,lo_sd,hi_sd", source_value=[r.hr_sd, r.lo_sd, r.hi_sd], rule="hr_ci", role="value")
    S.text(Q_R_C, y + OFF_C, fmt_p(r.q_sd), PT["value"], bold=True, align="right", panel="c", source_file=SPEC["results_continuous"], source_key=f"key={r.key}:q_sd", source_value=r.q_sd, rule="fmt_p", role="value")
    drawnC.append(dict(key=r.key, outcome=r.outcome, hr=r.hr_sd, lo=r.lo_sd, hi=r.hi_sd, q=r.q_sd, x=axis_x(aC, bC, r.hr_sd), y=ym, x_lo=axis_x(aC, bC, r.lo_sd), x_hi=axis_x(aC, bC, r.hi_sd), s_hr=hr_ci(r.hr_sd, r.lo_sd, r.hi_sd), s_q=fmt_p(r.q_sd)))
S.line(X0_C, SPINE_C["y"], X1_C, SPINE_C["y"], INK, SPINE_C["w"], cap=2)
for tx, lab in zip(TKX_C, tlC): S.line(tx, SPINE_C["y"], tx, SPINE_C["y"] + TICK_LEN, INK, SPINE_C["w"], cap=0); S.text(tx, TL_Y, lab["text"], PT["tick"], align="center", panel="c", source_key="axis tick", role="tick")
AXC_C = (X0_C + X1_C) / 2
S.text(AXC_C, TITLE_Y_C, titleC["text"], PT["axis_title"], align="center", panel="c", source_key="static V13 string (axis title, centred on the axis)", role="axis_title")
S.save(FLAT); dump(S.drawn + A_DRAWN, f"{W}/{SHEET}_bc_drawn.json"); dump({"b": drawnB, "c": drawnC}, f"{W}/{SHEET}_rows.json")

# ---- 4b round 54: sizes by role, clearances, overlaps (every drawn record carries x_left, width, baseline, size, role)
DRW = S.drawn
def rec(text, panel=None): return [d for d in DRW if d["text"] == text and (panel is None or d["panel"] == panel)]
def right(d): return d["x_left"] + d["width"]
def box(d, pad=0.0): return (d["x_left"] - pad, d["baseline"] - 0.716 * d["size"] - pad, right(d) + pad, d["baseline"] + 0.212 * d["size"] + pad)   # cap top to descender bottom
def inter(a, b): return a[0] < b[2] and b[0] < a[2] and a[1] < b[3] and b[1] < a[3]
roles = sorted(set((d["role"], d["size"]) for d in DRW))
ck.log(f"round 54 sizes: every drawn text carries a role and is at the round-54 size of that role {roles}; no text under 12 pt", all(d["role"] and abs(d["size"] - PT[d["role"]]) < 1e-9 for d in DRW) and min(d["size"] for d in DRW) >= 12.0)
keyR = [right(d) for d in DRW if d["role"] == "key"]
ck.log(f"round 54 key: both entries at {PT['key']} pt end at x {[round(v, 2) for v in keyR]}; they cross the reference line's x {REF_X_B:.2f}, so the line starts at y {AXTOP_B:.2f}, {AXTOP_B - max(box(d)[3] for d in DRW if d['role'] == 'key'):.2f} pt under the key's lowest descender; the key ends at least 8 pt before the head 'HR (95% CI)' starts at {rec('HR (95% CI)', 'b')[0]['x_left']:.2f}; the key sits above the first row ({B_Y[0]:.2f})",
       AXTOP_B - max(box(d)[3] for d in DRW if d["role"] == "key") >= 5.9 and max(keyR) + 8 <= rec("HR (95% CI)", "b")[0]["x_left"] and all(k["baseline"] < B_Y[0] - 20 for k in KEY_DRAWN))
labR = {d["text"]: right(d) for d in DRW if d["panel"] == "b" and d["role"] == "label"}
ck.log(f"round 54 panel b labels at {PT['label']} pt (drawn widths): the widest ends at x {max(labR.values()):.2f} ('{max(labR, key=labR.get)}'), every label at least {AIR_LAB_REF_B} pt left of the reference line {REF_X_B:.2f} and {AIR_LAB_WHISK} pt left of its row's whisker start", all(labR[t] <= REF_X_B - AIR_LAB_REF_B + 1e-6 for t in labR) and all(right(d) <= ROW_LEFT[k] - AIR_LAB_WHISK + 1e-6 for k in LAB for d in DRW if d["panel"] == "b" and d["role"] == "label" and d["text"] in LAB[k]), f"{sorted(((round(v, 1), t) for t, v in labR.items()), reverse=True)[:3]}")
hrL = [d["x_left"] for d in DRW if d["panel"] == "b" and d["rule"] == "hr_ci"]; pL = [d["x_left"] for d in DRW if d["panel"] == "b" and d["rule"] == "fmt_p"]
ck.log(f"round 54 panel b columns at {PT['value']} pt: HR strings right-aligned at {HR_R_B:.2f} start at x >= {min(hrL):.2f} (axis right end {X1_B:.2f}, widest whisker end {max(d['x_hi'] for d in drawnB):.2f}, at least {AIR_WHISK_COL} pt of air); P strings right-aligned at {P_R_B:.2f} start at x >= {min(pL):.2f} (at least {COL_GAP} pt after the HR edge); the head 'HR (95% CI)' ends at {right(rec('HR (95% CI)', 'b')[0]):.2f} at least 6 pt before the head 'P' starts at {rec('P', 'b')[0]['x_left']:.2f}; the P column ends at {P_R_B:.2f} = {RIGHT_B}, {PW - P_R_B:.2f} pt inside the page",
       min(hrL) >= max(d["x_hi"] for d in drawnB) + AIR_WHISK_COL - 0.01 and min(pL) >= HR_R_B + COL_GAP - 0.01 and right(rec("HR (95% CI)", "b")[0]) <= rec("P", "b")[0]["x_left"] - 6 and P_R_B < PW - 20)
labCr = {d["text"]: right(d) for d in DRW if d["panel"] == "c" and d["role"] == "label"}
hrLc = [d["x_left"] for d in DRW if d["panel"] == "c" and d["rule"] == "hr_ci"]; qLc = [d["x_left"] for d in DRW if d["panel"] == "c" and d["rule"] == "fmt_p"]
ck.log(f"round 54 panel c at {PT['label']} / {PT['value']} pt: the widest label ends at x {max(labCr.values()):.2f} ('{max(labCr, key=labCr.get)}'), at least {AIR_LAB_REF_C} pt left of the reference line {REF_X_C:.2f} and {AIR_LAB_WHISK} pt left of every row's whisker start ({min(d['x_lo'] for d in drawnC):.2f} at the least); HR strings start at x >= {min(hrLc):.2f} (widest whisker end {max(d['x_hi'] for d in drawnC):.2f}); q strings start at x >= {min(qLc):.2f} (at least {COL_GAP} pt after the HR edge {HR_R_C:.2f}); the head 'HR (95% CI)' ends at {right(rec('HR (95% CI)', 'c')[0]):.2f} before the head 'q' at {rec('q', 'c')[0]['x_left']:.2f}; the q column ends at {Q_R_C:.2f}, {B_X - Q_R_C:.1f} pt before panel b's labels",
       max(labCr.values()) <= REF_X_C - AIR_LAB_REF_C + 1e-6 and all(right(d) <= dd["x_lo"] - AIR_LAB_WHISK + 1e-6 for d in DRW if d["panel"] == "c" and d["role"] == "label" for dd in drawnC if plab(dd["outcome"]) == d["text"]) and min(hrLc) >= max(d["x_hi"] for d in drawnC) + AIR_WHISK_COL - 0.01 and min(qLc) >= HR_R_C + COL_GAP - 0.01 and right(rec("HR (95% CI)", "c")[0]) <= rec("q", "c")[0]["x_left"] - 6 and B_X - Q_R_C >= GAP_CB - 0.01)
# glyph-level air of the two-line label against its own second line and against the neighbouring rows (Arial glyf outlines)
_ARIAL_TT = _TTFont(ARIAL); _UPM = _ARIAL_TT["head"].unitsPerEm; _CMAP = _ARIAL_TT.getBestCmap(); _GLYF = _ARIAL_TT["glyf"]
def glyph_box_em(ch):
    g = _GLYF[_CMAP[ord(ch)]]
    if g.numberOfContours == 0: return (0.0, 0.0, 0.0, 0.0)
    return (g.xMin / _UPM, g.yMin / _UPM, g.xMax / _UPM, g.yMax / _UPM)
def glyph_boxes(d):
    out = []
    for i, ch in enumerate(d["text"]):
        if ch == " ": continue
        x0e, y0e, x1e, y1e = glyph_box_em(ch); x0 = d["x_left"] + S0.width(d["text"][:i], d["size"], d["bold"])
        out.append((ch, x0 + x0e * d["size"], x0 + x1e * d["size"], d["baseline"] - y1e * d["size"], d["baseline"] - y0e * d["size"]))
    return out
def v_air(upper, lower):
    """Smallest vertical air (pt) between any glyph of the upper record and any x-overlapping glyph of the lower record."""
    g1, g2 = glyph_boxes(upper), glyph_boxes(lower)
    return min(((round(b[3] - a[4], 2), a[0], b[0]) for a in g1 for b in g2 if a[1] < b[2] and b[1] < a[2]), default=(99.0, "", ""))
two_l = sorted([d for d in DRW if d["panel"] == "b" and d["rule"] == "label_line"], key=lambda d: d["baseline"]); d1, d2 = two_l
labels_b = sorted([d for d in DRW if d["panel"] == "b" and d["role"] == "label"], key=lambda d: d["baseline"])
i1 = labels_b.index(d1); prev_lab = labels_b[i1 - 1]; next_lab = labels_b[i1 + 2]
airs = {"line 1 over line 2": v_air(d1, d2), "previous row over line 1": v_air(prev_lab, d1), "line 2 over the next row": v_air(d2, next_lab)}
ck.log(f"round 54 two-line label ('{d1['text']}' / '{d2['text']}' at baselines {d1['baseline']:.2f} and {d2['baseline']:.2f}, {d2['baseline'] - d1['baseline']:.2f} pt apart at {d1['size']} pt): at least 1.5 pt of air between the glyph outlines of the two lines and against the neighbouring rows '{prev_lab['text']}' and '{next_lab['text']}': {airs} (pt of air, glyph above, glyph below)", all(v[0] >= 1.5 for v in airs.values()))
# every pair of drawn strings: no overlapping text boxes (cap top to descender bottom, exact widths)
ov = [(a["text"], b["text"]) for i, a in enumerate(DRW) for b in DRW[i + 1:] if inter(box(a), box(b))]
ck.log(f"round 54 no overlapping text boxes among the {len(DRW)} drawn strings of panels b and c (boxes from the drawn widths, cap top to descender bottom): overlapping pairs {ov}", not ov)
# text against marks: no row label, key, head, value or tick label box intersects a whisker or a marker (the rows' whiskers as boxes ym +- half the marker)
marks = [(d["x_lo"] - B_CI / 2, d["y"] - B_CI / 2, d["x_hi"] + B_CI / 2, d["y"] + B_CI / 2) for d in drawnB] + [(d["x"] - B_SQ / 2, d["y"] - B_SQ / 2, d["x"] + B_SQ / 2, d["y"] + B_SQ / 2) for d in drawnB]   # whiskers with their round caps, markers (the square is the larger)
marks += [(d["x_lo"] - C_CI / 2, d["y"] - C_CI / 2, d["x_hi"] + C_CI / 2, d["y"] + C_CI / 2) for d in drawnC] + [(d["x"] - C_DOT / 2, d["y"] - C_DOT / 2, d["x"] + C_DOT / 2, d["y"] + C_DOT / 2) for d in drawnC]
marks += [(KEY_H0 - 1, k["marker_y"] - B_DOT / 2, KEY_H1 + 1, k["marker_y"] + B_DOT / 2) for k in KEY_DRAWN]
refs = [(REF_X_B - 1, AXTOP_B, REF_X_B + 1, SPINE_B["y"]), (REF_X_C - 1, AXTOP_C, REF_X_C + 1, SPINE_C["y"])]
spines = [(X0_B, SPINE_B["y"] - 1, X1_B, SPINE_B["y"] + TICK_LEN + 1), (X0_C, SPINE_C["y"] - 1, X1_C, SPINE_C["y"] + TICK_LEN + 1)]
tm = [(d["text"], m) for d in DRW if d["role"] != "letter" for m in marks + refs + spines if inter(box(d), m)]
ck.log(f"round 54 no text box touches a whisker (with its round caps), a marker, a key handle, a dashed reference line or an axis spine with its ticks (every drawn string against {len(marks)} mark boxes, 2 reference lines, 2 spines); the smallest air between a row label's right end and its whisker cap: {min((round(m[0] - right(d), 2), d["text"]) for d in DRW if d["role"] == "label" for m in marks if m[1] < d["baseline"] < m[3] + 8 and m[0] > right(d))}: hits {tm[:6]}", not tm)
titR = max(right(d) for d in DRW if d["role"] == "axis_title"); titL = {d["text"]: (d["x_left"], right(d)) for d in DRW if d["role"] == "axis_title"}
ck.log(f"round 54 axis titles at {PT['axis_title']} pt centred on their axes ({AXC_C:.2f} and {AXC_B:.2f}): {titL}; panel c's title ends at least 20 pt before panel b's titles start; right ends at most {titR:.2f} (page width {PW})", titR < PW - 4 and max(v[1] for t, v in titL.items() if t.startswith("Hazard ratio per")) + 20 <= min(v[0] for t, v in titL.items() if not t.startswith("Hazard ratio per")))
tlb = [d for d in DRW if d["role"] == "tick"]
ck.log(f"round 54 tick labels at {PT['tick']} pt at baseline {TL_Y:.2f}: cap tops {TL_AIR} pt under the tick ends ({SPINE_B['y'] + TICK_LEN:.2f}); neighbouring tick labels at least 6 pt apart: {[(a['text'], b['text'], round(b['x_left'] - right(a), 1)) for a, b in zip(sorted(tlb, key=lambda d: d['x_left'])[:-1], sorted(tlb, key=lambda d: d['x_left'])[1:]) if abs(a['baseline'] - b['baseline']) < 0.1 and a['panel'] == b['panel']]}", all(b["x_left"] - right(a) >= 6 for a, b in zip(sorted(tlb, key=lambda d: d["x_left"])[:-1], sorted(tlb, key=lambda d: d["x_left"])[1:]) if a["panel"] == b["panel"]))
letc = rec("c", "c")[0]; hc = [d for d in DRW if d["panel"] == "c" and d["role"] == "head"]
ck.log(f"round 54 letters at {PT['letter']} pt: b at ({B_X:.2f}, {letters['b']['y']}) with the key {KEY_Y[0] - 0.716 * PT['key'] - letters['b']['y']:.1f} pt under its baseline; c at (36.0, {letters['c']['y']}) with its cap box top at {box(letc)[1]:.1f}, under the a/c clip boundary {CLIP_AC} (V13 393; panel a's lowest text ends at {PA_BOTTOM:.1f}), and the heads at {HEAD_Y_C:.2f} to its right", box(letc)[1] > CLIP_AC + 4 and PA_BOTTOM < CLIP_AC - 4 and KEY_Y[0] - 0.716 * PT["key"] - letters["b"]["y"] > 6 and all(h["x_left"] > right(letc) + 20 for h in hc))

# ---- 5 verify the flat panels (get_drawings allowed on the flat page): markers and whiskers read back at the file values on the NEW axes
itemsF = drawings(FLAT); MF = markers(itemsF)
dotsB = sorted([m for m in MF if m["shape"] == "circle" and abs(m["w"] - B_DOT) < 0.06 and m["cx"] > 500 and m["cy"] > 100], key=lambda m: (m["cy"], m["cx"])); sqB = sorted([m for m in MF if m["shape"] == "square" and abs(m["w"] - B_SQ) < 0.06 and m["cx"] > 500], key=lambda m: (m["cy"], m["cx"]))
expB = sorted([d for d in drawnB if not d["control"]], key=lambda d: (d["y"], d["x"])); expBc = sorted([d for d in drawnB if d["control"]], key=lambda d: (d["y"], d["x"]))
ab2, bb2, rb2 = fit_axis([d["hr"] for d in expB], [m["cx"] for m in dotsB]) if len(dotsB) == len(expB) else (0, 0, np.array([9]))
ck.log(f"flat panels: {len(dotsB)} orange dots and {len(sqB)} grey squares in b (one per row) read back at ln(HR) on the round-54 axis (fit a {ab2:.2f} b {bb2:.3f} = the layout's {aB:.2f} / {bB:.3f}, residual under 0.05 pt), filled exactly where P < 0.05 ({sum(not m['open'] for m in dotsB)} filled dots), at y = the label baseline {OFF_MARK:+.2f}",
       len(dotsB) == NS and len(sqB) == NCn and float(np.abs(rb2).max()) < 0.05 and abs(ab2 - aB) < 0.05 and abs(bb2 - bB) < 0.05 and all((not m["open"]) == d["filled"] for m, d in zip(dotsB, expB)) and all((not m["open"]) == d["filled"] for m, d in zip(sqB, expBc)) and all(abs(m["cy"] - d["y"]) < 0.05 for m, d in zip(dotsB, expB)), f"res {np.abs(rb2).max():.3f}")
dotsC = sorted([m for m in MF if m["shape"] == "circle" and abs(m["w"] - C_DOT) < 0.06 and m["cx"] < 400], key=lambda m: m["cy"])
ac2, bc2, rc2 = fit_axis([d["hr"] for d in drawnC], [m["cx"] for m in dotsC]) if len(dotsC) == len(drawnC) else (0, 0, np.array([9]))
ck.log(f"flat panels: {len(drawnC)} blue dots in c at ln(hr_sd) on the round-54 axis (fit a {ac2:.2f} b {bc2:.3f} = the layout's {aC:.2f} / {bC:.3f}, residual under 0.05 pt)", len(dotsC) == len(drawnC) and float(np.abs(rc2).max()) < 0.05 and abs(ac2 - aC) < 0.05 and abs(bc2 - bC) < 0.05, f"res {np.abs(rc2).max():.3f}")
HL = hlines(itemsF, min_len=5)
def whisker_ok(d, width, col):
    hit = [l for l in HL if abs(l["y"] - d["y"]) < 0.05 and abs(l["x0"] - d["x_lo"]) < 0.05 and abs(l["x1"] - d["x_hi"]) < 0.05 and abs(l["width"] - width) < 0.02 and l["stroke"] == col]
    return len(hit) == 1
wb = [whisker_ok(d, B_CI, d["colour"]["line"]) for d in drawnB]; wc = [whisker_ok(d, C_CI, BLUE) for d in drawnC]
ck.log(f"flat panels: every whisker read back (get_drawings) from axis_x(lo) to axis_x(hi) of the file values at the row's marker y, width {B_CI:.2f} / {C_CI:.2f} pt, blue for conditions and grey for controls: b {sum(wb)} of {len(wb)}, c {sum(wc)} of {len(wc)}", all(wb) and all(wc))
colB = {(m["colour"], m["open"]) for m in dotsB} | {(m["colour"], m["open"]) for m in sqB}
ck.log(f"flat panels: marker colours as V13 (orange dots, grey squares, filled where P < 0.05, open otherwise; blue dots in c): {sorted(colB, key=str)} and {sorted({m['colour'] for m in dotsC})}", all(c in (ORANGE, GREY) for c, _ in colB) and {m["colour"] for m in dotsC} == {BLUE})
keyd = sorted([m for m in MF if m["shape"] == "circle" and abs(m["w"] - B_DOT) < 0.06 and m["cy"] < 100], key=lambda m: m["cy"])
sp2, _ = text_spans(FLAT); sp2 = collapse(sp2); keyt = sorted([s for s in sp2 if s["text"] in [k["text"] for k in keyB]], key=lambda s: s["y"])
ck.log(f"key measured on the flat page: two entries {[s['text'] for s in keyt]} at x {[round(s['x'], 2) for s in keyt]} (= labels' x {B_X:.2f} + handle {KEY_HANDLE:.2f} + gap {KEY_GAP:.2f} = {KEY_LABEL_X:.2f}), baselines {[s['y'] for s in keyt]}, key dot at x {[round(m['cx'], 2) for m in keyd]} (= handle centre {(KEY_H0 + KEY_H1) / 2:.2f}); no '1 to 10% on PAP' anywhere",
       [s["text"] for s in keyt] == ["Still above 10% on PAP", "Reference: 1% or less on PAP"] and all(abs(s["x"] - KEY_LABEL_X) < 0.3 for s in keyt) and len(keyd) == 1 and abs(keyd[0]["cx"] - (KEY_H0 + KEY_H1) / 2) < 0.3 and [round(s["y"], 2) for s in keyt] == [round(y, 2) for y in KEY_Y] and not any("1 to 10%" in s["text"] for s in sp2), "")
pcol = column(sp2, P_R_B - 60, P_R_B, ymin=100, ymax=700); exp_p = [(d["s_p"], d["filled"]) for d in drawnB]
ck.log(f"P strings on the flat page: {len(pcol)} in the column, texts and bold faces as drawn (bold exactly where P < 0.05: {sum(b for _, b in exp_p)} of {len(exp_p)}); no 'q' string in panel b", [(s["text"], "Bold" in s["font"]) for s in pcol] == exp_p and not [s for s in sp2 if s["text"] == "q" and s["x"] > 497], "")
spine_lines = sorted(set(round(l["y"], 3) for l in hlines(itemsF, min_len=100) if abs(l["width"] - SPINE_B["w"]) < 0.02))
ck.log(f"axis lines on the flat page (get_drawings, width {SPINE_B['w']:.2f} pt, length over 100 pt): y {spine_lines}: panel b's spine and panel c's spine on the same line {SPINE_C['y']}", spine_lines == [SPINE_C["y"]], "")
sizes_flat = sorted(set(round(s["size"], 2) for s in sp2))
ck.log(f"sizes on the flat page (PyMuPDF on the flat page): {sizes_flat} = the round-54 sizes 14, 15, 18", sizes_flat == [14.0, 15.0, 18.0])
dump({"flat": FLAT, "flat_sha256": sha256(FLAT), "round": 54,
      "sizes": PT, "layout_constants": dict(DX_B=DX_B, RIGHT_B=RIGHT_B, GAP_CB=GAP_CB, AIR_LAB_WHISK=AIR_LAB_WHISK, AIR_LAB_REF_B=AIR_LAB_REF_B, AIR_LAB_REF_C=AIR_LAB_REF_C, AIR_WHISK_COL=AIR_WHISK_COL, COL_GAP=COL_GAP, OFF_MARK=OFF_MARK, NUDGE_2L=NUDGE_2L, KEY_DY=KEY_DY, TL_AIR=TL_AIR, TITLE_GAP=TITLE_GAP, TITLE_DY=TITLE_DY),
      "panel_b_axis": {"a": aB, "b": bB, "a_v13": aB13, "b_v13": bB13, "scale_vs_v13": bB / bB13, "spine_x": [X0_B, X1_B], "spine_x_v13": [SPINE_B["x0"], SPINE_B["x1"]], "spine_y": SPINE_B["y"], "ticks": dict(zip([str(t) for t in TICKV_B], TKX_B)), "ref_x": REF_X_B, "ref_x_v13": REF_X_B13, "ref_top": AXTOP_B, "ref_top_v13": AXTOP_B13, "hr_right": HR_R_B, "p_right": P_R_B, "hr_right_v13": HR_R_B13, "p_right_v13": P_R_B13, "title_baselines": TITLE_Y_B, "title_centre_x": AXC_B, "tick_label_baseline": TL_Y, "binding": whyB},
      "panel_c_axis": {"a": aC, "b": bC, "a_v13": aC13, "b_v13": bC13, "scale_vs_v13": bC / bC13, "spine_x": [X0_C, X1_C], "spine_x_v13": [SPINE_C["x0"], SPINE_C["x1"]], "spine_y": SPINE_C["y"], "ticks": dict(zip([str(t) for t in TICKV_C], TKX_C)), "ref_x": REF_X_C, "ref_x_v13": REF_X_C13, "ref_top": AXTOP_C, "hr_right": HR_R_C, "q_right": Q_R_C, "hr_right_v13": HR_R_C13, "q_right_v13": Q_R_C13, "title_baseline": TITLE_Y_C, "title_centre_x": AXC_C, "tick_label_baseline": TL_Y, "binding": whyC},
      "rows_b": drawnB, "rows_c": drawnC, "panel_a": [dict(item=l, text=e) for l, e, _, _, _ in pa], "panel_a_r54_page": {"pdf": PA_PDF, "sha256": sha256(PA_PDF), "text_on_sheet": pa_sheet, "right": PA_RIGHT, "bottom": PA_BOTTOM, "top": PA_TOP}, "sig_b_p": {"B": n_sigB}, "series": "B only (Still above 10% on PAP vs the reference)",
      "layout_b": {"rows_y": yS, "controls_y": yC, "head_y": y_head, "pitch": P_NEW, "v13_pitch": B_PITCH, "n_pitches": N_PITCH, "spine_y": SPINE_B["y"], "label_x": B_X, "label_x_v13": B_X13, "head_baseline": HEAD_Y_B}, "layout_c": {"rows_y": [C_Y[0] + i * C_PITCH for i in range(len(drawnC))], "pitch": C_PITCH, "label_x": C_X, "head_baseline": HEAD_Y_C},
      "selection": {"rule": f"the {N_TOP} non-control outcomes with the largest adj_pub_hr (L4_fig5b_refit_v3.json), descending at full precision, then the negative controls by adj_pub_hr", "n_file": N_FILE, "n_noncontrol": len(NC_ALL), "events_min": EV_MIN, "rows": [k for k, _ in SIG], "controls": [k for k, _ in CTRL], "adj_pub_hr": {k: v["pub_hr"] for k, v in SIG + CTRL}, "cut_19th": [NC_ALL[N_TOP][0], NC_ALL[N_TOP][1]["pub_hr"]], "ties_3dp": TIES_3DP},
      "key": KEY_DRAWN, "key_geometry": {"handle": KEY_HANDLE, "gap": KEY_GAP, "h0": KEY_H0, "h1": KEY_H1, "label_x": KEY_LABEL_X, "v13_label_x": 694.31, "baselines": KEY_Y, "dy": KEY_DY}, "labels": {"lines": LAB, "two_line": two, "widest": [max(LAB_END, key=LAB_END.get), max(LAB_END.values())], "two_line_nudge_pt": NUDGE_2L, "two_line_air": airs},
      "content_bottom": CONTENT_BOTTOM, "clip_ac": CLIP_AC, "provenance": {"treatment_v2": sc_t, "treatment_v3": sc_t3, "fig5b_refit_v3_json": sc_r, "results_continuous": sc_c}, "letters_v13": letters, "title_v13": title, "page_v13": [PW, PH], "n_drawn": len(S.drawn), "built": time.strftime("%Y-%m-%d %H:%M:%S")}, f"{W}/flat_build_record.json")
ok = ck.write(f"{V}/checks_flat_build.txt"); print("FLAT BUILD RESULT ALL PASS" if ok else "FLAT BUILD RESULT FAIL")
