#!$T90_PY
"""Round 49, lane L5 (2026-09-26): the edits applied to the copied scripts (recorded here, run once). Copies: scripts/l5_lib.py (round-42
lfig5b_lib.py repointed to this lane, title 14 pt), scripts/gs_text.py, Main_Fig5/scripts/l5a14.py (Sheet.text and Sheet.width +1.0 pt,
BUMP; width(bump=False) for the V13 geometry read-back), 05_build_flat_r49.py (round-42 builder: this lane's paths, the V13 column edges
reconstructed at the V13 sizes, extra clearance checks at +1 pt), 04_build_panel_a_r49.py (the panel a flow generator of
FINAL_FIGURES_2026-08-14/_scripts/figure5A_flow_polish.py writing into this lane, text sizes + BUMP/scale)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import re, os
L5 = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/L5"


def patch(path, pairs, must=True):
    s = open(path).read()
    for old, new in pairs:
        if old not in s:
            if must: raise SystemExit(f"pattern not found in {path}: {old[:80]!r}")
            continue
        s = s.replace(old, new)
    open(path, "w").write(s); print("patched", path)


# ---- l5_lib.py
patch(f"{L5}/scripts/l5_lib.py", [
    ('"""Lane LFIG5b (round 42, 2026-09-19, the P-column variant): shared helpers,', '"""Round 49, lane L5 (2026-09-26, every figure text +1 pt): the round-42 LFIG5b library repointed to this lane (title 14 pt).\nLane LFIG5b (round 42, 2026-09-19, the P-column variant): shared helpers,'),
    ('R42 = f"{paths.FIGURE_ROOT}/ROUND42_2026-09-19"; LANE = f"{R42}/lanes/LFIG5b"; WD = f"{R38}/figures/scripts/wd_run.sh"; WD_LOGDIR = f"{LANE}/logs/watchdog"',
     'R42 = f"{paths.FIGURE_ROOT}/ROUND42_2026-09-19"; R44 = f"{paths.FIGURE_ROOT}/ROUND44_2026-09-21"; R49 = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26"; LANE = f"{R49}/lanes/L5"; WD = f"{LANE}/scripts/wd_run.sh"; WD_LOGDIR = f"{LANE}/logs/watchdog"\nLFIG5B = f"{R42}/lanes/LFIG5b"; V26 = f"{paths.FIGURE_ROOT}/ROUND47_2026-09-24/figures/NEW_FINAL_SET_V26"; R44_VECTOR = f"{R44}/lanes/LSKETCH/Main_Fig5/work/Main_Fig5_r44_vector.pdf"'),
    ('SKETCH = f"{R40}/lanes/LS_SKETCH/work/build"; SKETCH_STATUS = f"{R40}/lanes/LS_SKETCH/STATUS.md"', 'SKETCH = f"{R49}/lanes/LSKETCH/work/build"; SKETCH_R44 = f"{R44}/lanes/LSKETCH/work/build"; SKETCH_STATUS = f"{R49}/lanes/LSKETCH/STATUS.md"'),
    ('TITLE_SIZE, TITLE_BASELINE = 13.0, 14.0', 'TITLE_SIZE, TITLE_BASELINE = 14.0, 14.0          # round 49: the sheet title and the panel letters 13 -> 14 pt\nLETTER_SIZE = 14.0'),
    ('lines = ["# LMAIN STATUS (round 40)"', 'lines = ["# L5 STATUS (round 49)"'),
])
patch(f"{L5}/scripts/gs_text.py", [("import lfig5b_lib as L", "import l5_lib as L")])

# ---- l5a14.py: the page helper, +1 pt on every text
patch(f"{L5}/Main_Fig5/scripts/l5a14.py", [
    ('"""Lane V14_L5a_PAP_MAIN (round 37, 2026-09-15): shared library.', '"""Round 49, lane L5 (2026-09-26): the round-42 copy of this library with the page helper Sheet drawing every text ONE POINT larger\n(Sheet.BUMP = 1.0 inside text() and width(); width(..., bump=False) gives the V13 width for the geometry read-back of the V13 sheet).\nLane V14_L5a_PAP_MAIN (round 37, 2026-09-15): shared library.'),
    ('    def width(self, s, size, bold=False): return (self.bold if bold else self.font).text_length(s, fontsize=size)\n',
     '    BUMP = 1.0   # round 49: every text one point larger than the size the builder asks for (placement rules see the larger text)\n\n    def width(self, s, size, bold=False, bump=True): return (self.bold if bold else self.font).text_length(s, fontsize=size + (self.BUMP if bump else 0.0))\n'),
    ('        w = self.width(s, size, bold); x0 = x - w if align == "right" else (x - w / 2 if align == "center" else x)\n',
     '        size_asked = size; size = size + self.BUMP\n        w = self.width(s, size, bold, bump=False); x0 = x - w if align == "right" else (x - w / 2 if align == "center" else x)\n'),
    ('                               source_file=source_file, source_key=source_key, source_value=source_value, rule=rule, note=note, x_left=round(x0, 3), width=round(w, 3)))',
     '                               source_file=source_file, source_key=source_key, source_value=source_value, rule=rule, note=note, x_left=round(x0, 3), width=round(w, 3), size_asked=size_asked))'),
])

# ---- 05_build_flat_r49.py
P5 = f"{L5}/Main_Fig5/scripts/05_build_flat_r49.py"
patch(P5, [
    ('#!$T90_PY\n"""Round 42, lane LFIG5b, FINAL spec', '#!$T90_PY\n"""Round 49, lane L5 (2026-09-26): the round-42 LFIG5b flat builder of panels b and c rerun with every text one point larger (l5a14.Sheet.BUMP\n= 1.0 inside text() and width()); the V13 column edges (right-aligned HR and P/q columns) are reconstructed at the V13 sizes (width(bump=False))\nso the columns keep their right edges and grow leftward; every other constant unchanged; extra clearance checks at +1 pt appended. Round-42 docstring:\nRound 42, lane LFIG5b, FINAL spec'),
    ('SHEET = "Main_Fig5"; D = f"{paths.FIGURE_ROOT}/ROUND42_2026-09-19/lanes/LFIG5b/Main_Fig5"; W = f"{D}/work"',
     'SHEET = "Main_Fig5"; D = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/L5/Main_Fig5"; W = f"{D}/work"'),
    ('ck = Checks(f"Lane LFIG5b round 42 (one contrast per condition, P column, V13 axis), {SHEET} flat build,',
     'ck = Checks(f"Lane L5 round 49 (every text +1 pt: Sheet.BUMP {Sheet.BUMP}), the round-42 LFIG5b flat build (one contrast per condition, P column, V13 axis), {SHEET} flat build,'),
    ('HR_R_B = float(np.median([s["x"] + S0.width(s["text"], s["size"]) for s in hrB])); P_R_B = float(np.median([s["x"] + S0.width(s["text"], s["size"], "Bold" in s["font"]) for s in pB]))',
     'HR_R_B = float(np.median([s["x"] + S0.width(s["text"], s["size"], bump=False) for s in hrB])); P_R_B = float(np.median([s["x"] + S0.width(s["text"], s["size"], "Bold" in s["font"], bump=False) for s in pB]))'),
    ('assert np.std([s["x"] + S0.width(s["text"], s["size"]) for s in hrB]) < 0.25 and np.std([s["x"] + S0.width(s["text"], s["size"], "Bold" in s["font"]) for s in pB]) < 0.2',
     'assert np.std([s["x"] + S0.width(s["text"], s["size"], bump=False) for s in hrB]) < 0.25 and np.std([s["x"] + S0.width(s["text"], s["size"], "Bold" in s["font"], bump=False) for s in pB]) < 0.2'),
    ('HR_R_C = float(np.mean([s["x"] + S0.width(s["text"], s["size"]) for s in hrC])); Q_R_C = float(np.mean([s["x"] + S0.width(s["text"], s["size"], True) for s in qC]))',
     'HR_R_C = float(np.mean([s["x"] + S0.width(s["text"], s["size"], bump=False) for s in hrC])); Q_R_C = float(np.mean([s["x"] + S0.width(s["text"], s["size"], True, bump=False) for s in qC]))'),
])
# the round-49 clearance checks, appended before the record dump (after S.save)
s = open(P5).read()
anchor = 'S.save(FLAT); dump(S.drawn + A_DRAWN, f"{W}/{SHEET}_bc_drawn.json"); dump({"b": drawnB, "c": drawnC}, f"{W}/{SHEET}_rows.json")\n'
assert anchor in s
extra = anchor + '''
# ---- 4b round 49: clearance at +1 pt (every drawn record carries x_left and width at the drawn size)
DRW = S.drawn
def rec(text, panel=None): return [d for d in DRW if d["text"] == text and (panel is None or d["panel"] == panel)]
def right(d): return d["x_left"] + d["width"]
keyR = [right(d) for d in DRW if d["panel"] == "b" and d["text"] in ("Still above 10% on PAP", "Reference: 1% or less on PAP")]
ck.log(f"round 49 key: both key labels end at x {[round(v, 2) for v in keyR]} at {KEY2[0]['size'] + Sheet.BUMP:.2f} pt, at least 2 pt left of the dashed reference line {REF_X_B} and above the first condition row ({B_Y[0]:.2f}, key baselines {[round(k['baseline'], 2) for k in KEY_DRAWN]})", all(v < REF_X_B - 2 for v in keyR) and all(k["baseline"] < B_Y[0] - 20 for k in KEY_DRAWN))
labR = {d["text"]: right(d) for d in DRW if d["panel"] == "b" and abs(d["x"] - B_X) < 0.01 and d["baseline"] > 100 and d["text"] != "Negative controls"}
ck.log(f"round 49 panel b labels at {10.45 + Sheet.BUMP:.2f} pt: the widest ends at x {max(labR.values()):.2f} ('{max(labR, key=labR.get)}'), every label at least 2 pt left of its row's leftmost interval and of the reference line (rechecked on the drawn widths)", all(labR[t] < REF_X_B - 2 for t in labR) and all(right(d) < ROW_LEFT[k] - 2 for k in LAB for d in DRW if d["panel"] == "b" and abs(d["x"] - B_X) < 0.01 and d["text"] in LAB[k]), f"{sorted(((round(v, 1), t) for t, v in labR.items()), reverse=True)[:3]}")
hrL = [d["x_left"] for d in DRW if d["panel"] == "b" and d["rule"] == "hr_ci"]; pL = [d["x_left"] for d in DRW if d["panel"] == "b" and d["rule"] == "fmt_p"]
ck.log(f"round 49 panel b columns at {9.02 + Sheet.BUMP:.2f} pt: HR strings right-aligned at the V13 edge {HR_R_B:.2f} start at x >= {min(hrL):.2f} (axis right end {SPINE_B['x1']}, widest whisker end {max(d['x_hi'] for d in drawnB):.2f}); P strings right-aligned at {P_R_B:.2f} start at x >= {min(pL):.2f} (> the HR edge + 2); the head 'HR (95% CI)' ends at {right(rec('HR (95% CI)', 'b')[0]):.2f} before the head 'P' at {rec('P', 'b')[0]['x_left']:.2f}", min(hrL) > max(d["x_hi"] for d in drawnB) + 2 and min(pL) > HR_R_B + 2 and right(rec("HR (95% CI)", "b")[0]) < rec("P", "b")[0]["x_left"] - 2)
labC = {d["text"]: right(d) for d in DRW if d["panel"] == "c" and d["rule"] == "label"}
hrLc = [d["x_left"] for d in DRW if d["panel"] == "c" and d["rule"] == "hr_ci"]; qLc = [d["x_left"] for d in DRW if d["panel"] == "c" and d["rule"] == "fmt_p"]
ck.log(f"round 49 panel c at {9.5 + Sheet.BUMP:.2f} / {9.02 + Sheet.BUMP:.2f} pt: the widest label ends at x {max(labC.values()):.2f} ('{max(labC, key=labC.get)}'), at least 2 pt left of the reference line {REF_X_C} and of every row's whisker start ({min(d['x'] for d in drawnC):.2f} at the least); HR strings start at x >= {min(hrLc):.2f} (widest whisker end {max(axis_x(aC, bC, r.hi_sd) for r in Cc.itertuples()):.2f}); q strings start at x >= {min(qLc):.2f} (> the HR edge {HR_R_C:.2f} + 2); the head 'HR (95% CI)' ends at {right(rec('HR (95% CI)', 'c')[0]):.2f} before the head 'q' at {rec('q', 'c')[0]['x_left']:.2f}", max(labC.values()) < REF_X_C - 2 and min(hrLc) > max(axis_x(aC, bC, r.hi_sd) for r in Cc.itertuples()) + 2 and min(qLc) > HR_R_C + 2 and right(rec("HR (95% CI)", "c")[0]) < rec("q", "c")[0]["x_left"] - 2)
two_l = [d for d in DRW if d["panel"] == "b" and d["rule"] == "label_line"]
ck.log(f"round 49 two-line label: lines at baselines {[round(d['baseline'], 2) for d in two_l]} ({abs(two_l[1]['baseline'] - two_l[0]['baseline']):.2f} pt apart) at {two_l[0]['size']:.2f} pt: descender of line 1 (0.212 em) clears the ascender of line 2 (0.716 em cap, 0.75 em with 'd')", len(two_l) == 2 and abs(two_l[1]["baseline"] - two_l[0]["baseline"]) > (0.212 + 0.75) * two_l[0]["size"] + 0.3)
titR = max(right(d) for d in DRW if d["text"] in ("Hazard ratio for a new diagnosis", "vs 1% or less on PAP (95% CI)", "Hazard ratio per 1 SD of residual sleep T90 (95% CI)"))
tb_new = max(d["baseline"] for d in DRW if d["text"] == "vs 1% or less on PAP (95% CI)")
ck.log(f"round 49 axis titles: right ends at most {titR:.2f} (page width {PW}); panel b's second title line baseline {tb_new:.2f} (V13 715.93, the band top 747.80 minus 20 keeps clear)", titR < PW - 4 and tb_new + 0.212 * (10.45 + Sheet.BUMP) < 747.8 - 20)
sizes_drawn = sorted(set((d["size_asked"], d["size"]) for d in DRW))
ck.log(f"round 49 sizes: every drawn text is the asked size + {Sheet.BUMP} pt: {sizes_drawn}", all(abs(b - a - Sheet.BUMP) < 1e-9 for a, b in sizes_drawn))
'''
s = s.replace(anchor, extra); open(P5, "w").write(s); print("patched (appended checks)", P5)

# ---- 04_build_panel_a_r49.py: the flow generator writing into this lane, text sizes + BUMP
P4 = f"{L5}/Main_Fig5/scripts/04_build_panel_a_r49.py"
patch(P4, [
    ('"""\nFigure5A, the pre/post FLOW render (round 3, 2026-08-25):', '"""\nRound 49, lane L5 (2026-09-26): the Figure 5a flow generator (FINAL_FIGURES_2026-08-14/_scripts/figure5A_flow_polish.py, the panel the V13\nsheet carries at scale S_V13 = 1.1342) rerun with every text size + BUMP where BUMP = 1/S_V13 pt on this canvas (= +1.0 pt on the sheet):\nusage 04_build_panel_a_r49.py <bump_on_canvas_pt> <tag>. Writes only into this lane (work/panel_a/). Geometry in data units unchanged, so\nthe flows, blocks and brackets (the data marks) do not move. The Aug-25 docstring follows.\nFigure5A, the pre/post FLOW render (round 3, 2026-08-25):'),
    ('OUT = f"{paths.FIGURE_ROOT}/FINAL_FIGURES_2026-08-14"\nPOLISH = f"{OUT}/Main_Figures_Polished/Figure5"\nWORK = f"{OUT}/_workfiles"\nPNGS = f"{WORK}/polished_pngs"\nfor _d in (POLISH, WORK, PNGS):\n    os.makedirs(_d, exist_ok=True)\nNAME = "Figure5A"',
     'BUMP = float(sys.argv[1]) if len(sys.argv) > 1 else 0.0\nTAG = sys.argv[2] if len(sys.argv) > 2 else "bump0"\nOUT = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/L5/Main_Fig5/work/panel_a"\nPOLISH = OUT\nWORK = OUT\nPNGS = OUT\nfor _d in (POLISH, WORK, PNGS):\n    os.makedirs(_d, exist_ok=True)\nNAME = f"Figure5A_{TAG}"'),
    ('TITLE_PT, TICK_PT, ANN_PT, FLOOR_PT = 11.0, 10.0, 9.5, 9.0\n', 'TITLE_PT, TICK_PT, ANN_PT, FLOOR_PT = 11.0 + BUMP, 10.0 + BUMP, 9.5, 9.0     # round 49: + BUMP on this canvas\n'),
    ('json.dump(DRAWN, open(f"{WORK}/{NAME}_flow_drawn_values.json", "w"), indent=1,\n          default=float)', 'DRAWN["bump_on_canvas_pt"] = BUMP; DRAWN["text_pt_on_canvas"] = {"blocks_arms": TICK_PT, "captions": TITLE_PT}\njson.dump(DRAWN, open(f"{WORK}/{NAME}_flow_drawn_values.json", "w"), indent=1,\n          default=float)'),
])
print("done")
