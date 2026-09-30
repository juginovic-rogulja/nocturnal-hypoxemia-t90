#!$T90_PY
"""Round 42, lane LFIG5b, FINAL spec (Alen, 2026-09-19 late, through the coordinator): panel b keeps ONE contrast per condition, 'Still above 10%
on PAP' against the reference (restored to 1% or less); the '1 to 10% on PAP' line (adj_A_*) is dropped everywhere (rows, key, columns). One
single-line row per condition (the 18 largest adj_pub_hr, then the three controls), the upper-line series only (adj_B_hr/lo/hi, P = adj_B_p,
bold and the orange dot filled where P < 0.05), heads 'HR (95% CI)' and 'P', a two-entry key (Still above 10% on PAP; Reference: 1% or less
on PAP) on the left over the labels. The 21 rows spread evenly between the V13 first-row y and the V13 last-row y at ONE pitch, the
'Negative controls' head in its V13 block (1.5 pitches below the last condition, 1.5 above the first control): 17 + 3 + 2 = 22 pitches of
24.670 pt (V13 23.598); the axis at the V13 spine 673.85 on panel c's line. The earlier two-series P build is kept beside as
05_build_flat_r42b_TWO_SERIES_superseded.py. Its docstring follows.
Round 42, lane LFIG5b (the coordinator's two changes to the LFIG5 sheet, 2026-09-19 late): (i) the printed column is the UNCORRECTED P of each
contrast (adj_B_p, adj_A_p; fmt_p; bold and the orange dot filled where P < 0.05; head 'P'), no q in panel b; (ii) panel b's x-axis line, tick
labels and axis title sit at the V13 positions (spine y 673.85 = panel c's axis), the 21 rows at the V13 pitch from the V13 first-row y, the
spare pitch (22 V13 rows minus 21) absorbed as extra gap between the last condition row and the 'Negative controls' head, so the last control
row sits at the V13 last-row position. Everything else as LFIG5 (copy beside as 05_build_flat_r42_LFIG5_REFERENCE.py). LFIG5 docstring:
Round 42 (2026-09-19, lane LFIG5): Main Fig 5 panels b and c drawn on the flat page work/panels_bc_flat.pdf, the round-37/40 builder
(05_build_fig5.py steps 0 to 5, kept beside as 05_build_flat_r40_REFERENCE.py) with the owner's decisions of 2026-09-19:
  (1) panel b rows = the 18 non-control conditions with the LARGEST adjusted two-group hazard ratio (adj_pub_hr of L4_fig5b_refit_v3.json,
      step 212, adjusted for prevalent cardiopulmonary disease), sorted by that ratio descending at full precision, then the three negative
      controls the file carries (Hemorrhoids, Glaucoma, Contact dermatitis) by the same ratio; the file holds only outcomes with at least 40
      events (asserted); no P-value cut;
  (2) the printed columns are HR (95% CI) and q: the Benjamini-Hochberg q of each contrast (adj_B_q for 'Still above 10% on PAP', adj_A_q for
      '1 to 10% on PAP'), q printed with the sheet's fmt_p idiom, bold where q < 0.05; the orange dots filled where q < 0.05, open otherwise
      (the Figures 2a and 5c convention); the column head 'P' becomes 'q';
  (3) the key (Still above 10% on PAP / 1 to 10% on PAP / Reference: 1% or less on PAP) moves to the LEFT of panel b, above the condition
      labels, left-aligned at the labels' x (the handles start at the labels' x, the V13 handle length and gap kept), stacked at the V13 rows.
Panel c unchanged (results_continuous.csv, step 147). Panel a is the V13 clip (its seven counts re-read from treatment_v2.json here).
Geometry = the V13 design (pitch, fonts, markers, axis), 18 + 3 rows at the V13 pitch, the axis and its title moved with the rows.
Output: work/panels_bc_flat.pdf, work/Main_Fig5_bc_drawn.json, work/Main_Fig5_rows.json, work/flat_build_record.json, verify/checks_flat_build.txt."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from l5a14 import *

SHEET = "Main_Fig5"; D = f"{paths.FIGURE_ROOT}/ROUND42_2026-09-19/lanes/LFIG5b/Main_Fig5"; W = f"{D}/work"; V = f"{D}/verify"; os.makedirs(f"{V}/crops", exist_ok=True)
BASE = f"{V13}/{SHEET}.pdf"; FLAT = f"{W}/panels_bc_flat.pdf"; N_TOP = 18
ck = Checks(f"Lane LFIG5b round 42 (one contrast per condition, P column, V13 axis), {SHEET} flat build, {time.strftime('%Y-%m-%d %H:%M')}. BASE (design) = V13 {BASE} (sha256 {sha256(BASE)[:16]}). NEW = L4_fig5b_refit_v3.json (212, adjusted, BH q) for panel b: the {N_TOP} largest adj_pub_hr + 3 controls, HR (95% CI) and P (uncorrected, adj_B_p / adj_A_p) columns; results_continuous.csv (147) for panel c; treatment_v2.json (117) for the panel a counts.")
SC = 0.95   # the panels b and c were placed on the sheet at 0.95 of the generators' 171 mm canvases (round 12 recompose)

# ---- 0 gates
sc_t = gate("treatment_v2"); sc_t3 = gate("treatment_v3"); sc_r = gate("fig5b_refit_v3_json"); sc_c = gate("results_continuous")
dump({"treatment_v2": sc_t, "treatment_v3": sc_t3, "fig5b_refit_v3_json": sc_r, "results_continuous": sc_c}, f"{W}/sidecars.json")
ck.log("treatment_v2.json (117), treatment_v3.json (211), L4_fig5b_refit_v3.json (212, adjusted, --bh) and results_continuous.csv (147) pass the v8.1 sidecar gate and the v8.2 mtime gate", all(s["v8_inputs"] for s in (sc_t, sc_t3, sc_r, sc_c)), f"out {sc_t['output_mtime']}, {sc_t3['output_mtime']}, {sc_r['output_mtime']}, {sc_c['output_mtime']}")

# ---- 1 the V13 text layer (cached probe) and the clip-render geometry (the round-37 block, unchanged)
J = load_json(f"{W}/base_text.json"); spans, info = J["spans"], J["info"]; PW, PH = info["rect"][2], info["rect"][3]
ck.log("V13 page 968.94 x 946.23, nested (454 XObjects, no get_drawings)", abs(PW - 968.94) < 0.05 and abs(PH - 946.2253) < 0.05, info)
letters = {s["text"]: s for s in spans if len(s["text"]) == 1 and s["text"] in "abcd" and s["size"] > 12.5}; assert set(letters) == set("abcd")
title = find_span(spans, "Figure 5"); assert title
labB = [s for s in column(spans, 497, 498, ymin=100, ymax=700, size=10.45) if "Bold" not in s["font"]]; headB = find_span(spans, "Negative controls", xmin=490); assert len(labB) == 22 and headB
B_Y = [s["y"] for s in labB]; B_PITCH = float(np.mean(np.diff(B_Y[:17]))); B_X = float(np.mean([s["x"] for s in labB]))
hrB = column(spans, 830, 836, ymin=100, ymax=700); pB = column(spans, 900, 926, ymin=100, ymax=700); assert len(hrB) == 44 and len(pB) == 44, (len(hrB), len(pB))
OFF_B_UP = float(np.mean([hrB[2 * i]["y"] - s["y"] for i, s in enumerate(labB)])); OFF_B_LO = float(np.mean([hrB[2 * i + 1]["y"] - s["y"] for i, s in enumerate(labB)]))
S0 = Sheet(10, 10)   # font metrics only
HR_R_B = float(np.median([s["x"] + S0.width(s["text"], s["size"]) for s in hrB])); P_R_B = float(np.median([s["x"] + S0.width(s["text"], s["size"], "Bold" in s["font"]) for s in pB]))
assert np.std([s["x"] + S0.width(s["text"], s["size"]) for s in hrB]) < 0.25 and np.std([s["x"] + S0.width(s["text"], s["size"], "Bold" in s["font"]) for s in pB]) < 0.2
keyB = sorted([s for s in spans if abs(s["x"] - 694.31) < 0.5], key=lambda s: s["y"]); assert [s["text"] for s in keyB] == ["Still above 10% on PAP", "1 to 10% on PAP", "Reference: 1% or less on PAP"]
colhB = [find_span(spans, "HR (95% CI)", xmin=800, ymax=80), find_span(spans, "P", xmin=900, ymax=80)]; assert all(colhB)
tlB = sorted([s for s in spans if abs(s["y"] - 687.49) < 0.5 and s["x"] > 600], key=lambda s: s["x"]); assert [s["text"] for s in tlB] == ["0.5", "1", "2", "4", "8"]
titleB = sorted([s for s in spans if s["x"] > 600 and 700 < s["y"] < 720], key=lambda s: s["y"]); assert [s["text"] for s in titleB] == ["Hazard ratio for a new diagnosis", "vs 1% or less on PAP (95% CI)"]
SPINE_B = dict(y=673.85, x0=616.92, x1=824.64, w=0.8 * SC); TICKS_B = [646.98, 687.42, 727.74, 768.18, 808.5]; TICK_LEN = 3.0 * SC; AXTOP_B = 50.84; REF_X_B = 687.35
SPINE_C = dict(y=673.85, x0=155.52, x1=363.36, w=0.8 * SC); TICKS_C = [164.34, 231.78, 279.66, 347.22]; AXTOP_C = 432.4; REF_X_C = 164.3
aB, bB, rB = fit_axis([0.5, 1, 2, 4, 8], TICKS_B); aC, bC, rC = fit_axis([1, 1.5, 2, 3], TICKS_C); assert max(np.abs(rB).max(), np.abs(rC).max()) < 0.15, (rB, rC)
TL_OFF_B = float(np.mean([s["y"] for s in tlB])) - SPINE_B["y"]
LD = load_json(f"{paths.FIGURE_ROOT}/ROUND30_2026-09-08/LD_fig5/work/markers_geometry.json"); DOT_R = float(np.median([o["page_radius"] for o in LD["markers"]]))
ld_x = sorted(o["page_center"][0] for o in LD["markers"])
ck.log(f"V13 panel b geometry: 22 rows pitch {B_PITCH:.3f} at x {B_X:.2f}, strings at {OFF_B_UP:+.2f} and {OFF_B_LO:+.2f}, HR column right edge {HR_R_B:.2f}, P/q column right edge {P_R_B:.2f}; spine y {SPINE_B['y']} x {SPINE_B['x0']} to {SPINE_B['x1']}; log axis a {aB:.2f} b {bB:.3f} (residual {np.abs(rB).max():.3f}); dot radius {DOT_R} (round-30 record)", abs(aB - REF_X_B) < 0.3, f"tick label offset {TL_OFF_B:.2f}")
labC = column(spans, 36, 37, ymin=430, ymax=680, size=9.5); assert len(labC) == 11; C_Y = [s["y"] for s in labC]; C_PITCH = float(np.mean(np.diff(C_Y))); C_X = float(np.mean([s["x"] for s in labC]))
hrC = column(spans, 370, 373, ymin=430, ymax=680); qC = column(spans, 440, 455, ymin=430, ymax=680); assert len(hrC) == 11 and len(qC) == 11
HR_R_C = float(np.mean([s["x"] + S0.width(s["text"], s["size"]) for s in hrC])); Q_R_C = float(np.mean([s["x"] + S0.width(s["text"], s["size"], True) for s in qC]))
OFF_C = float(np.mean([h["y"] - s["y"] for h, s in zip(hrC, labC)]))
colhC = [find_span(spans, "HR (95% CI)", xmin=380, xmax=390, ymin=415, ymax=430), find_span(spans, "q", xmin=460, ymin=415, ymax=430)]; assert all(colhC)
tlC = sorted([s for s in spans if abs(s["y"] - 686.82) < 0.5 and s["x"] < 400], key=lambda s: s["x"]); assert [s["text"] for s in tlC] == ["1", "1.5", "2", "3"]
titleC = find_span(spans, "Hazard ratio per 1 SD of residual sleep T90 (95% CI)"); assert titleC
ck.log(f"V13 panel c geometry: 11 rows pitch {C_PITCH:.3f} at x {C_X:.2f}, columns right edges {HR_R_C:.2f} and {Q_R_C:.2f}, strings at {OFF_C:+.2f}; spine y {SPINE_C['y']}; log axis a {aC:.2f} b {bC:.3f} (residual {np.abs(rC).max():.3f})", abs(aC - REF_X_C) < 0.3, "")
C_DOT, C_EDGE, C_CI = 5.29 * SC, 0.8 * SC, 1.7 * SC; B_DOT, B_EDGE_F, B_EDGE_O, B_CI, B_SQ = 2 * DOT_R, 0.8 * SC, 1.1 * SC, 1.8 * SC, 6.0 * SC
OFF_C_MARK = -3.58 * SC; OFF_B_MARK = -2.44
ck.log("round-30 LD record: the V13 5b dots sit 2.44 pt above their string baselines and the key dots at x 677.2", any(abs(x - 677.21) < 0.05 for x in ld_x) and any(abs(x - 677.54) < 0.05 for x in ld_x), "")
KEY_HANDLE = 0.30 * 72 * SC; KEY_GAP = 0.10 * 72 * SC                        # the builder's 0.30 in handle ending 0.10 in before the label (V13: 666.95 to 687.47, label 694.31)
KEY_H0, KEY_H1, KEY_LABEL_X = B_X, B_X + KEY_HANDLE, B_X + KEY_HANDLE + KEY_GAP   # round 42: the key block left-aligned at the labels' x

# ---- 2 panel a: the seven counts re-read from the v8.1 file must equal the V13 strings (then the region is kept)
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
ck.log("panel a: all 7 count strings of the V13 region equal treatment_v2.json (1,894; 613; 498; 215; 568; 1,326 (70%); 568 (30%)): the region is kept from V13 with this read-back as its proof", all(a_found), f"{[exp for _, exp, _, _, _ in pa]}")
A_DRAWN = [dict(text=exp, x=g["x"], baseline=g["y"], ha="left", size=g["size"], panel="a", source_file=SPEC["treatment_v2"], source_key=key, source_value=val, rule=rule, note="V13 region kept, value re-read from the v8.1 file") for (lab, exp, key, val, rule), g in zip(pa, a_found) if g]

# ---- 3 panel b rows (decision 1) and panel c rows
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
n_sigB = sum(v["B_p"] < 0.05 for _, v in SIG); n_sigA = sum(v["A_p"] < 0.05 for _, v in SIG)
ck.log(f"panel b (LFIG5b change 1): the printed column is the uncorrected P of each contrast (adj_B_p / adj_A_p of the adjusted refit); P < 0.05 on the 18 rows: still above 10% in {n_sigB}, 1 to 10% in {n_sigA}; controls with P < 0.05: {[k for k, v in CTRL if v['B_p'] < 0.05 or v['A_p'] < 0.05]}; for the record, BH q < 0.05 would mark {sum(v['B_q'] < 0.05 for _, v in SIG)} + {sum(v['A_q'] < 0.05 for _, v in SIG)} cells", all(0 <= v[f'{s}_p'] <= 1 for _, v in SIG + CTRL for s in 'AB'), "")
cont, Cc = residual_rows()
ck.log(f"panel c (unchanged): {len(Cc)} Benjamini-Hochberg survivors of results_continuous.csv, all q < 0.05 and lower CI above 1; at the V13 pitch the last row sits above the axis", (Cc.q_sd < 0.05).all() and (Cc.lo_sd > 1.0).all() and C_Y[0] + (len(Cc) - 1) * C_PITCH + 6 < SPINE_C["y"], list(Cc.outcome))
xloB, xhiB = exp((SPINE_B["x0"] - aB) / bB), exp((SPINE_B["x1"] - aB) / bB); xloC, xhiC = exp((SPINE_C["x0"] - aC) / bC), exp((SPINE_C["x1"] - aC) / bC)
ck.log("every drawn value (the B contrast) lies inside the V13 axis ranges (no tick change)", all(xloB <= v["B_lo"] and v["B_hi"] <= xhiB for _, v in SIG + CTRL) and all(xloC <= r.lo_sd and r.hi_sd <= xhiC for r in Cc.itertuples()), f"b {xloB:.3f}-{xhiB:.3f}, c {xloC:.3f}-{xhiC:.3f}")
NS, NCn = len(SIG), len(CTRL)
N_PITCH = (NS - 1) + 3 + (NCn - 1)                                   # 17 intervals between the conditions, the V13 three-pitch head block, 2 between the controls = 22
P_NEW = (B_Y[-1] - B_Y[0]) / N_PITCH                                 # 24.670 pt: the 21 single-line rows spread evenly from the V13 first row to the V13 last row
yS = [B_Y[0] + i * P_NEW for i in range(NS)]; y_head = B_Y[0] + (NS - 1 + 1.5) * P_NEW; yC = [B_Y[0] + (NS - 1 + 3 + j) * P_NEW for j in range(NCn)]
SPINE_B_NEW = SPINE_B["y"]; DELTA_B = 0.0
assert abs(yS[0] - B_Y[0]) < 1e-9 and abs(yC[-1] - B_Y[-1]) < 1e-9, (yS[0], yC[-1])
ck.log(f"panel b layout (final spec): {NS} single-line rows then the head then {NCn} controls at ONE pitch {P_NEW:.3f} pt (V13 23.598) from the V13 first row {B_Y[0]:.2f} to the V13 last row {B_Y[-1]:.2f} (22 pitches: 17 + the V13 three-pitch head block + 2); the head at {y_head:.2f} (1.5 pitches below the last condition {yS[-1]:.2f}, 1.5 above the first control {yC[0]:.2f}); the axis, tick labels and title at the V13 y (spine {SPINE_B_NEW} = panel c's {SPINE_C['y']})", SPINE_B_NEW == SPINE_C["y"] and P_NEW > B_PITCH, f"controls {[round(y, 2) for y in yC]}")
# row labels: one line at the V13 size; a label wider than the label gutter (497.45 to the spine's left end 616.92: 'Ventricular arrhythmia or
# cardiac arrest', 172 pt; its Figure 2a short form 'Ventricular arrhythmia/arrest' still 132 pt) is printed on TWO lines at the row's two value
# baselines, the Figure 4b idiom of this design family ('Ventricular arrhythmia' / 'or cardiac arrest'), the full file name kept
GUTTER = SPINE_B["x0"] - B_X
def label_lines(k):
    if S0.width(k, 10.45) <= GUTTER - 2: return [k]
    if " or " in k: a, b = k.split(" or ", 1); return [a, "or " + b]
    return [plab(k)]
LAB = {k: label_lines(k) for k, _ in SIG + CTRL}; LAB_END = {k: max(B_X + S0.width(t, 10.45) for t in LAB[k]) for k in LAB}
ROW_LEFT = {k: axis_x(aB, bB, min(v["A_lo"], v["B_lo"])) for k, v in SIG + CTRL}
two = [k for k in LAB if len(LAB[k]) == 2]; worst = max(LAB, key=lambda k: LAB_END[k])
ck.log(f"row labels: {len(two)} label(s) on two lines {[(k, LAB[k]) for k in two]}; the widest line ends at x {LAB_END[worst]:.1f} ('{worst}', a V13 label at its V13 width; the spine's left end is {SPINE_B['x0']}, the reference line {REF_X_B}); every label ends at least 2 pt left of its own row's leftmost interval and of the reference line", all(LAB_END[k] < ROW_LEFT[k] - 2 and LAB_END[k] < REF_X_B - 2 for k in LAB), f"row-left minima {min(ROW_LEFT.values()):.1f}; labels past the spine's left end: {[(k, round(LAB_END[k], 1)) for k in LAB if LAB_END[k] > SPINE_B['x0']]}")

# ---- 4 draw panels b and c on the flat page
S = Sheet(PW, PH)
S.text(letters["b"]["x"], letters["b"]["y"], "b", 13.0, bold=True, color="#" + letters["b"]["color"], panel="b", source_key="panel letter"); S.text(letters["c"]["x"], letters["c"]["y"], "c", 13.0, bold=True, color="#" + letters["c"]["color"], panel="c", source_key="panel letter")
S.line(REF_X_B, AXTOP_B, REF_X_B, SPINE_B_NEW, INK20, 0.9 * SC, cap=0, dashes=f"[{4 * SC:g} {3 * SC:g}] 0")
# the key (decision 3): the V13 rows and text, the block left-aligned at the labels' x (handle from B_X, the V13 handle length and gap)
KEY_DRAWN = []
KEY2 = [dict(keyB[0]), dict(keyB[2], y=keyB[1]["y"])]                  # two entries: 'Still above 10% on PAP' and 'Reference: 1% or less on PAP', stacked at the first two V13 key rows
for s, kind in zip(KEY2, ("B", "ref")):
    yk = s["y"] + OFF_B_MARK
    if kind == "ref": S.line(KEY_H0, yk, KEY_H1, yk, INK20, 0.9 * SC, cap=0, dashes=f"[{4 * SC:g} {3 * SC:g}] 0")
    else: S.line(KEY_H0, yk, KEY_H1, yk, BLUE if kind == "B" else BLUE_LT, B_CI, cap=1); S.circle((KEY_H0 + KEY_H1) / 2, yk, B_DOT, ORANGE, WHITE, B_EDGE_F)
    S.text(KEY_LABEL_X, s["y"], s["text"], s["size"], color="#" + s["color"], panel="b", source_key="key label (V13 text; round 42: two entries, the key block moved to the left, left-aligned at the labels' x)")
    KEY_DRAWN.append(dict(text=s["text"], kind=kind, handle_x=[KEY_H0, KEY_H1], label_x=KEY_LABEL_X, baseline=s["y"], marker_y=yk, v13_label_x=s["x"]))
S.text(colhB[0]["x"], colhB[0]["y"], colhB[0]["text"], colhB[0]["size"], color="#" + colhB[0]["color"], panel="b", source_key="static V13 string")
S.text(colhB[1]["x"], colhB[1]["y"], colhB[1]["text"], colhB[1]["size"], color="#" + colhB[1]["color"], panel="b", source_key="static V13 string (column head 'P', LFIG5b change 1)")
drawnB = []
for (k, v), y, is_ctrl in [(kv, yy, False) for kv, yy in zip(SIG, yS)] + [(kv, yy, True) for kv, yy in zip(CTRL, yC)]:
    for series, off, name in (("B", 0.0, "Still above 10% on PAP"),):                     # final spec: the B contrast only, one line per row on the label baseline
        hr, lo, hi, q, p = v[f"{series}_hr"], v[f"{series}_lo"], v[f"{series}_hi"], v[f"{series}_q"], v[f"{series}_p"]; sig = p < 0.05
        ys_ = y + off; ym = ys_ + OFF_B_MARK
        line_c = GREY if is_ctrl else BLUE; mark_c = GREY if is_ctrl else ORANGE
        S.line(axis_x(aB, bB, lo), ym, axis_x(aB, bB, hi), ym, line_c, B_CI, cap=1)
        if is_ctrl: S.square(axis_x(aB, bB, hr), ym, B_SQ, mark_c if sig else WHITE, WHITE if sig else mark_c, B_EDGE_F if sig else B_EDGE_O)
        else: S.circle(axis_x(aB, bB, hr), ym, B_DOT, mark_c if sig else WHITE, WHITE if sig else mark_c, B_EDGE_F if sig else B_EDGE_O)
        S.text(HR_R_B, ys_, hr_ci(hr, lo, hi), 9.02, color=INK20, align="right", panel="b", source_file=SPEC["fig5b_refit_v3_json"], source_key=f"outcomes/{k}/adj_{series}_hr,adj_{series}_lo,adj_{series}_hi", source_value=[hr, lo, hi], rule="hr_ci", note=name + ", " + ADJ_NOTE)
        S.text(P_R_B, ys_, fmt_p(p), 9.02, bold=bool(sig), color=INK20, align="right", panel="b", source_file=SPEC["fig5b_refit_v3_json"], source_key=f"outcomes/{k}/adj_{series}_p", source_value=p, rule="fmt_p", note=name + ", uncorrected P of the adjusted contrast (LFIG5b change 1), bold where P < 0.05")
        drawnB.append(dict(outcome=k, label=plab(k), series=series, hr=hr, lo=lo, hi=hi, q=q, p=p, x=axis_x(aB, bB, hr), y=ym, x_lo=axis_x(aB, bB, lo), x_hi=axis_x(aB, bB, hi), filled=sig, control=is_ctrl, s_hr=hr_ci(hr, lo, hi), s_p=fmt_p(p), s_q=fmt_p(q)))
    if len(LAB[k]) == 1: S.text(B_X, y, LAB[k][0], 10.45, color=GREY if is_ctrl else INK20, panel="b", source_file=SPEC["fig5b_refit_v3_json"], source_key=f"outcomes/{k} (row label; decision 1: among the {N_TOP} largest adj_pub_hr, or a negative control)", source_value=k, rule="label", note="")
    else:
        for t, off in zip(LAB[k], (OFF_B_UP, OFF_B_LO)): S.text(B_X, y + off, t, 10.45, color=GREY if is_ctrl else INK20, panel="b", source_file=SPEC["fig5b_refit_v3_json"], source_key=f"outcomes/{k} (row label on two lines at the row's value baselines, the Figure 4b idiom; decision 1)", source_value=k, rule="label_line", note=f"line of the file label '{k}'")
S.text(headB["x"], y_head, "Negative controls", headB["size"], bold=True, color="#" + headB["color"], panel="b", source_key="static V13 string (moved with the rows)")
S.line(SPINE_B["x0"], SPINE_B_NEW, SPINE_B["x1"], SPINE_B_NEW, INK20, SPINE_B["w"], cap=2)
for tx, lab in zip(TICKS_B, tlB): S.line(tx, SPINE_B_NEW, tx, SPINE_B_NEW + TICK_LEN, INK20, SPINE_B["w"], cap=0); S.text(tx, lab["y"] + DELTA_B, lab["text"], lab["size"], color="#" + lab["color"], align="center", panel="b", source_key="axis tick")
for s in titleB: S.text(s["x"], s["y"] + DELTA_B, s["text"], s["size"], color="#" + s["color"], panel="b", source_key="static V13 string (moved with the axis)")
# panel c (unchanged)
S.line(REF_X_C, AXTOP_C, REF_X_C, SPINE_C["y"], INK, 0.9 * SC, cap=0, dashes=f"[{4 * SC:g} {3 * SC:g}] 0")
for s in colhC: S.text(s["x"], s["y"], s["text"], s["size"], color="#" + s["color"], panel="c", source_key="static V13 string")
drawnC = []
for i, r in enumerate(Cc.itertuples()):
    y = C_Y[0] + i * C_PITCH; ym = y + OFF_C_MARK
    S.line(axis_x(aC, bC, r.lo_sd), ym, axis_x(aC, bC, r.hi_sd), ym, BLUE, C_CI, cap=1); S.circle(axis_x(aC, bC, r.hr_sd), ym, C_DOT, BLUE, WHITE, C_EDGE)
    S.text(C_X, y, plab(r.outcome), 9.5, panel="c", source_file=SPEC["results_continuous"], source_key=f"key={r.key}:outcome", source_value=r.outcome, rule="label")
    S.text(HR_R_C, y + OFF_C, hr_ci(r.hr_sd, r.lo_sd, r.hi_sd), 9.02, align="right", panel="c", source_file=SPEC["results_continuous"], source_key=f"key={r.key}:hr_sd,lo_sd,hi_sd", source_value=[r.hr_sd, r.lo_sd, r.hi_sd], rule="hr_ci")
    S.text(Q_R_C, y + OFF_C, fmt_p(r.q_sd), 9.02, bold=True, align="right", panel="c", source_file=SPEC["results_continuous"], source_key=f"key={r.key}:q_sd", source_value=r.q_sd, rule="fmt_p")
    drawnC.append(dict(key=r.key, outcome=r.outcome, hr=r.hr_sd, lo=r.lo_sd, hi=r.hi_sd, q=r.q_sd, x=axis_x(aC, bC, r.hr_sd), y=ym))
S.line(SPINE_C["x0"], SPINE_C["y"], SPINE_C["x1"], SPINE_C["y"], INK, SPINE_C["w"], cap=2)
for tx, lab in zip(TICKS_C, tlC): S.line(tx, SPINE_C["y"], tx, SPINE_C["y"] + TICK_LEN, INK, SPINE_C["w"], cap=0); S.text(tx, lab["y"], lab["text"], lab["size"], align="center", panel="c", source_key="axis tick")
S.text(titleC["x"], titleC["y"], titleC["text"], titleC["size"], panel="c", source_key="static V13 string")
S.save(FLAT); dump(S.drawn + A_DRAWN, f"{W}/{SHEET}_bc_drawn.json"); dump({"b": drawnB, "c": drawnC}, f"{W}/{SHEET}_rows.json")

# ---- 5 verify the flat panels (get_drawings allowed on the flat page)
itemsF = drawings(FLAT); MF = markers(itemsF)
dotsB = sorted([m for m in MF if m["shape"] == "circle" and abs(m["w"] - B_DOT) < 0.06 and m["cx"] > 500 and m["cy"] > 100], key=lambda m: (m["cy"], m["cx"])); sqB = sorted([m for m in MF if m["shape"] == "square" and abs(m["w"] - B_SQ) < 0.06 and m["cx"] > 500], key=lambda m: (m["cy"], m["cx"]))
expB = sorted([d for d in drawnB if not d["control"]], key=lambda d: (d["y"], d["x"])); expBc = sorted([d for d in drawnB if d["control"]], key=lambda d: (d["y"], d["x"]))
ab2, bb2, rb2 = fit_axis([d["hr"] for d in expB], [m["cx"] for m in dotsB]) if len(dotsB) == len(expB) else (0, 0, np.array([9]))
ck.log(f"flat panels: {len(dotsB)} orange dots and {len(sqB)} grey squares in b (one per row) read back at ln(HR) on the V13 axis (residual under 0.05 pt), filled exactly where P < 0.05 ({sum(not m['open'] for m in dotsB)} filled dots)",
       len(dotsB) == NS and len(sqB) == NCn and float(np.abs(rb2).max()) < 0.05 and abs(ab2 - aB) < 0.05 and all((not m["open"]) == d["filled"] for m, d in zip(dotsB, expB)) and all((not m["open"]) == d["filled"] for m, d in zip(sqB, expBc)), f"res {np.abs(rb2).max():.3f}")
dotsC = sorted([m for m in MF if m["shape"] == "circle" and abs(m["w"] - C_DOT) < 0.06 and m["cx"] < 400], key=lambda m: m["cy"])
ac2, bc2, rc2 = fit_axis([d["hr"] for d in drawnC], [m["cx"] for m in dotsC]) if len(dotsC) == len(drawnC) else (0, 0, np.array([9]))
ck.log(f"flat panels: {len(drawnC)} blue dots in c at ln(hr_sd) on the V13 axis (residual under 0.05 pt)", len(dotsC) == len(drawnC) and float(np.abs(rc2).max()) < 0.05 and abs(ac2 - aC) < 0.05, f"res {np.abs(rc2).max():.3f}")
keyd = sorted([m for m in MF if m["shape"] == "circle" and abs(m["w"] - B_DOT) < 0.06 and m["cy"] < 100], key=lambda m: m["cy"])
sp2, _ = text_spans(FLAT); sp2 = collapse(sp2); keyt = sorted([s for s in sp2 if s["text"] in [k["text"] for k in keyB]], key=lambda s: s["y"])
ck.log(f"key measured on the flat page: two entries {[s['text'] for s in keyt]} at x {[round(s['x'], 2) for s in keyt]} (= labels' x {B_X:.2f} + handle {KEY_HANDLE:.2f} + gap {KEY_GAP:.2f} = {KEY_LABEL_X:.2f}), baselines {[s['y'] for s in keyt]} (the first two V13 key rows), key dot at x {[round(m['cx'], 2) for m in keyd]} (= handle centre {(KEY_H0 + KEY_H1) / 2:.2f}); the handles start at the labels' x {KEY_H0:.2f}; no '1 to 10% on PAP' anywhere",
       [s["text"] for s in keyt] == ["Still above 10% on PAP", "Reference: 1% or less on PAP"] and all(abs(s["x"] - KEY_LABEL_X) < 0.3 for s in keyt) and len(keyd) == 1 and abs(keyd[0]["cx"] - (KEY_H0 + KEY_H1) / 2) < 0.3 and [s["y"] for s in keyt] == [keyB[0]["y"], keyB[1]["y"]] and not any("1 to 10%" in s["text"] for s in sp2), "")
qhead = [s for s in sp2 if s["text"] == "q" and s["x"] > 497]; phead = [s for s in sp2 if s["text"] == "P" and s["y"] < 80 and s["x"] > 900]
ck.log(f"column head 'P' at its V13 origin ({phead[0]['x'] if phead else None}, {phead[0]['y'] if phead else None}); no 'q' string anywhere in panel b (x > 497)", len(phead) == 1 and abs(phead[0]["x"] - colhB[1]["x"]) < 0.3 and not qhead)
pcol = column(sp2, 895, 935, ymin=100, ymax=700); exp_p = [(d["s_p"], d["filled"]) for d in drawnB]
ck.log(f"P strings on the flat page: {len(pcol)} in the column, texts and bold faces as drawn (bold exactly where P < 0.05: {sum(b for _, b in exp_p)} of {len(exp_p)})", [(s["text"], "Bold" in s["font"]) for s in pcol] == exp_p, "")
spine_lines = sorted(set(round(l["y"], 3) for l in hlines(itemsF, min_len=150) if abs(l["width"] - SPINE_B["w"]) < 0.02))
ck.log(f"axis lines on the flat page (get_drawings, width {SPINE_B['w']:.2f} pt, length over 150 pt): y {spine_lines}: panel b's spine and panel c's spine on the same line {SPINE_C['y']}", spine_lines == [SPINE_C["y"]], "")
dump({"flat": FLAT, "flat_sha256": sha256(FLAT), "panel_b_axis": {"a": aB, "b": bB, "spine_y": SPINE_B_NEW, "delta": DELTA_B, "title_baselines": [s["y"] + DELTA_B for s in titleB], "tick_label_baseline": tlB[0]["y"] + DELTA_B},
      "panel_c_axis": {"a": aC, "b": bC, "spine_y": SPINE_C["y"]}, "rows_b": drawnB, "rows_c": drawnC, "panel_a": [dict(item=l, text=e) for l, e, _, _, _ in pa], "sig_b_p": {"B": n_sigB}, "sig_b_q_for_the_record": {"B": sum(v["B_q"] < 0.05 for _, v in SIG)}, "series": "B only (Still above 10% on PAP vs the reference)", "layout_b": {"rows_y": yS, "controls_y": yC, "head_y": y_head, "pitch": P_NEW, "v13_pitch": B_PITCH, "n_pitches": N_PITCH, "spine_y": SPINE_B_NEW},
      "selection": {"rule": f"the {N_TOP} non-control outcomes with the largest adj_pub_hr (L4_fig5b_refit_v3.json), descending at full precision, then the negative controls by adj_pub_hr", "n_file": N_FILE, "n_noncontrol": len(NC_ALL), "events_min": EV_MIN, "rows": [k for k, _ in SIG], "controls": [k for k, _ in CTRL], "adj_pub_hr": {k: v["pub_hr"] for k, v in SIG + CTRL}, "cut_19th": [NC_ALL[N_TOP][0], NC_ALL[N_TOP][1]["pub_hr"]], "ties_3dp": TIES_3DP},
      "key": KEY_DRAWN, "key_geometry": {"handle": KEY_HANDLE, "gap": KEY_GAP, "h0": KEY_H0, "h1": KEY_H1, "label_x": KEY_LABEL_X, "v13_label_x": 694.31}, "labels": {"lines": LAB, "gutter": GUTTER, "two_line": two, "widest": [worst, LAB_END[worst]]},
      "provenance": {"treatment_v2": sc_t, "treatment_v3": sc_t3, "fig5b_refit_v3_json": sc_r, "results_continuous": sc_c}, "letters_v13": letters, "title_v13": title, "page_v13": [PW, PH], "n_drawn": len(S.drawn), "built": time.strftime("%Y-%m-%d %H:%M:%S")}, f"{W}/flat_build_record.json")
ok = ck.write(f"{V}/checks_flat_build.txt"); print("FLAT BUILD RESULT ALL PASS" if ok else "FLAT BUILD RESULT FAIL")
