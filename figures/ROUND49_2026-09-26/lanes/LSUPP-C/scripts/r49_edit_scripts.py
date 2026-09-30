#!$T90_PY
"""Round 49, lane LSUPP-C: exact one-shot edits to the copied builders (repoint to this lane, every text size +1.0 pt through a
PT_PLUS constant, geometry pinned where a text measurement fed a layout constant, the round-40 LNOTES removal folded into
gen_supp13). Every old string must occur exactly once. Records scripts/R49_EDITS.md."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import os, sys
L = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LSUPP-C"
R40 = f"{paths.FIGURE_ROOT}/ROUND40_2026-09-18/lanes/LSUPP"
E = []
def edit(path, old, new, why):
    s = open(path).read(); n = s.count(old); assert n == 1, (path, n, old[:80])
    open(path, "w").write(s.replace(old, new)); E.append((os.path.relpath(path, L), old, new, why))

# ---------------- shared style modules
P = f"{L}/L5b_scripts/supp_polish_common.py"
edit(P, f'SCRIPTS = "{R40}/L5b_scripts"', f'SCRIPTS = "{L}/L5b_scripts"', "repoint")
edit(P, f'OUT = "{R40}/L5b_work/polish_out"', f'OUT = "{L}/L5b_work/polish_out"', "repoint")
edit(P, """HEAD_PT = 11.5                    # the headline, bold
TITLE_PT = 11.0                   # axis titles, regular
PTITLE_PT = 10.5                  # panel titles, regular
TICK_PT = 10.0                    # ticks, row labels, legends
ANN_PT = 9.5                      # load-bearing annotations, printed columns
FLOOR_PT = 9.0
PANEL_PT = 13.0                   # bold panel letters""",
"""PT_PLUS = 1.0                     # R49 (Alen, 2026-09-25): every text one point larger, geometry unchanged
HEAD_PT = 11.5 + PT_PLUS          # the headline, bold
TITLE_PT = 11.0 + PT_PLUS         # axis titles, regular
PTITLE_PT = 10.5 + PT_PLUS        # panel titles, regular
TICK_PT = 10.0 + PT_PLUS          # ticks, row labels, legends
ANN_PT = 9.5 + PT_PLUS            # load-bearing annotations, printed columns
FLOOR_PT = 9.0
PANEL_PT = 13.0 + PT_PLUS         # bold panel letters""", "+1 pt")
P = f"{L}/L5b_scripts/splitstyle.py"
edit(P, "TITLE, LABEL, TICK, ANNOT, SMALL, PANEL = 11.5, 10.0, 9.5, 9.5, 9.0, 13.0",
     "PT_PLUS = 1.0   # R49: every text one point larger\nTITLE, LABEL, TICK, ANNOT, SMALL, PANEL = 11.5 + PT_PLUS, 10.0 + PT_PLUS, 9.5 + PT_PLUS, 9.5 + PT_PLUS, 9.0 + PT_PLUS, 13.0 + PT_PLUS", "+1 pt (rc defaults)")
P = f"{L}/L5b_scripts/legend_capture.py"
edit(P, f'OUT = "{R40}/L5b_work/FIGURE_LEGENDS_lane_capture.md"', f'OUT = "{L}/L5b_work/FIGURE_LEGENDS_lane_capture.md"', "repoint")

# ---------------- gen_supp12
P = f"{L}/L5b_scripts/gen_supp12.py"
edit(P, f'SCRIPTS = "{R40}/L5b_scripts"', f'SCRIPTS = "{L}/L5b_scripts"', "repoint")
edit(P, f'FIN = "{R40}/L5b_work/polish_out"', f'FIN = "{L}/L5b_work/polish_out"', "repoint")
edit(P, "LABF, TCKF, ANNF, FLOOR = 11.0, 10.0, 9.5, 9.0", "PT_PLUS = 1.0   # R49: every text one point larger, geometry unchanged\nLABF, TCKF, ANNF, FLOOR = 11.0 + PT_PLUS, 10.0 + PT_PLUS, 9.5 + PT_PLUS, 9.0", "+1 pt")
edit(P, 'ha="left", va="baseline", fontsize=9.2, color=INK))', 'ha="left", va="baseline", fontsize=9.2 + PT_PLUS, color=INK))', "+1 pt (key labels, literal size)")

# ---------------- gen_supp14 (column widths measured at the round-40 size so the axes rect and every mark stay where V26 has them)
P = f"{L}/L5b_scripts/gen_supp14.py"
edit(P, f'sys.path.insert(0, "{R40}/L5b_scripts")', f'sys.path.insert(0, "{L}/L5b_scripts")', "repoint")
edit(P, f'OUT = "{R40}/L5b_work/polish_out"', f'OUT = "{L}/L5b_work/polish_out"', "repoint")
edit(P, "TITLE_PT, TICK_PT, ANN_PT, FLOOR_PT = 11.0, 10.0, 9.5, 9.0", "PT_PLUS = 1.0   # R49: every text one point larger, geometry unchanged\nTITLE_PT, TICK_PT, ANN_PT, FLOOR_PT = 11.0 + PT_PLUS, 10.0 + PT_PLUS, 9.5 + PT_PLUS, 9.0", "+1 pt")
edit(P, """        t = fig.text(0.5, 0.5, s, fontsize=COL_FS,
                     fontweight="bold" if bold else "normal")""",
"""        t = fig.text(0.5, 0.5, s, fontsize=COL_FS - PT_PLUS,   # R49: the column widths that set the axes rect are measured at the round-40 size (geometry pinned, marks unchanged)
                     fontweight="bold" if bold else "normal")""", "geometry pinned")

# ---------------- gen_supp15
P = f"{L}/L5b_scripts/gen_supp15.py"
edit(P, f'sys.path.insert(0, "{R40}/L5b_scripts")', f'sys.path.insert(0, "{L}/L5b_scripts")', "repoint")
edit(P, f'OUT = "{R40}/L5b_work/polish_out"', f'OUT = "{L}/L5b_work/polish_out"', "repoint")
edit(P, "TITLE_PT, TICK_PT, ANN_PT, FLOOR_PT = 11.0, 10.0, 9.5, 9.0", "PT_PLUS = 1.0   # R49: every text one point larger, geometry unchanged\nTITLE_PT, TICK_PT, ANN_PT, FLOOR_PT = 11.0 + PT_PLUS, 10.0 + PT_PLUS, 9.5 + PT_PLUS, 9.0", "+1 pt")

# ---------------- gen_supp13 (sizes come from supp_polish_common; the twin from this lane; the round-40 LNOTES removal folded in)
P = f"{L}/L5b_scripts/gen_supp13.py"
edit(P, 'LANE = f"{paths.FIGURE_ROOT}/ROUND37_2026-09-14/figures/lanes/V14_L5b_PAP_SUPP"', f'LANE = "{L}/L5b_work"   # R49: the twin regenerated in this lane (work/fullprec_nonresponder_v8_1.json)', "repoint (twin)")
edit(P, """    assert GUT + measure(fig, STAR_KEY, ANN_PT) <= W_IN - EDGE, "asterisk key runs off"
    fig.text(GUT / W_IN, (YA - A_XBAND - 0.03) / H_AC, STAR_KEY, fontsize=ANN_PT,
             ha="left", va="top", color=INK_P)""",
"""    # R49: the asterisk key is NOT printed (round-40 lane LNOTES removed it from the sheet, the sentence lives in the legend);
    # its band A_KEY stays so nothing else moves
    pass""", "round-40 LNOTES removal re-applied at the source")

# ---------------- fullprec twin generator (writes into this lane)
P = f"{L}/L5b_scripts/fullprec_nonresponder_v8_1.py"
edit(P, 'LANE = f"{paths.FIGURE_ROOT}/ROUND37_2026-09-14/figures/lanes/V14_L5b_PAP_SUPP"   # V14_L5b 2026-09-15: lane copy of V7_L12 s2, on the v8 tables', f'LANE = "{L}/L5b_work"   # R49: the twin regenerated in this lane, reconciled against the current v8.3 json', "repoint")

# ---------------- Supp_Fig20 builder
P = f"{L}/Supp_Fig20/scripts/01_build_subsheets.py"
edit(P, f'SCRIPTS = "{R40}/L5b_scripts"', f'SCRIPTS = "{L}/L5b_scripts"', "repoint")
edit(P, f'LANE20 = "{R40}/Supp_Fig20"', f'LANE20 = "{L}/Supp_Fig20"', "repoint")
edit(P, "LABF, TCKF, ANNF, FLOOR = 11.0, 10.0, 10.0, 9.0\nPTLAB = 9.5\nSTARF = 10.0",
     "PT_PLUS = 1.0   # R49: every text one point larger, geometry unchanged\nLABF, TCKF, ANNF, FLOOR = 11.0 + PT_PLUS, 10.0 + PT_PLUS, 10.0 + PT_PLUS, 9.0\nPTLAB = 9.5 + PT_PLUS\nSTARF = 10.0 + PT_PLUS", "+1 pt")
edit(P, "    TITLE_BAND = (ANNF * 1.25 + 6.0) / 72.0", "    TITLE_BAND = ((ANNF - PT_PLUS) * 1.25 + 6.0) / 72.0   # R49: the band that sets panel d's axes and the sub-sheet height stays at the round-40 value (geometry pinned)", "geometry pinned")

with open(f"{L}/scripts/R49_EDITS.md", "w") as f:
    f.write("# R49 edits to the copied builders (exact one-shot replacements, each asserted to occur once)\n\n")
    for p, o, n, w in E:
        f.write(f"- {p} [{w}]: `{o[:110]}` -> `{n[:110]}`\n")
print(f"{len(E)} edits applied")
import subprocess
r = subprocess.run(["grep", "-rn", "ROUND40_2026-09-18/lanes/LSUPP\\|ROUND37_2026-09-14/figures/lanes/V14_L5b", f"{L}/L5b_scripts", f"{L}/Supp_Fig20/scripts/01_build_subsheets.py"], capture_output=True, text=True)
live = [l for l in r.stdout.splitlines() if "#" not in l.split(":", 2)[-1][:5] and "_source" not in l]
print("grep proof, live lines still naming the round-40 or round-37 lane:", len(live)); print("\n".join(live[:20]))
