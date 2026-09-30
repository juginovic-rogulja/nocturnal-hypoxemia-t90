#!$T90_PY
"""Round 54, lane L5 (2026-09-29): the edits applied to the copied round-49 scripts (recorded here, run once). Copies (from
ROUND49_2026-09-26/lanes/L5, nothing there edited): scripts/l5_lib.py (repointed to this lane, V30 and the round-54 sketch folder, title and
letters 18 pt, the gs_flatten recipe of ROUND52 LSUPP-D lsd_common.py), scripts/gs_text.py (docstring), Main_Fig5/scripts/l5a14.py (the page
helper Sheet: BUMP 0.0, a ROLE_PT table with the round-54 sizes, the role recorded in every drawn record), 04_build_panel_a_r54.py (the flow
generator takes the two canvas sizes as arguments and writes into this lane), 04b_panel_a_map_r54.py and 06b_render_r54.py (paths).
The flat builder, the compose, the census and the verifier are rewritten for the round-54 layout (05_build_flat_r54.py, 06_compose_r54.py,
07_census_r54.py, 07b_verify_r54.py) and are not patched here."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import os
L5 = f"{paths.FIGURE_ROOT}/ROUND54_2026-09-28/lanes/L5"


def patch(path, pairs, must=True):
    s = open(path).read()
    for old, new in pairs:
        if old not in s:
            if must: raise SystemExit(f"pattern not found in {path}: {old[:90]!r}")
            continue
        assert s.count(old) == 1, (path, old[:60], s.count(old))
        s = s.replace(old, new)
    open(path, "w").write(s); print("patched", path)


# ---- l5_lib.py
patch(f"{L5}/scripts/l5_lib.py", [
    ('"""Round 49, lane L5 (2026-09-26, every figure text +1 pt): the round-42 LFIG5b library repointed to this lane (title 14 pt).',
     '"""Round 54, lane L5 (2026-09-29, the main figures at 14/15/18 pt): the round-49 L5 library repointed to this lane (title and letters 18 pt,\nV30 as the base set, the round-54 sketch folder, the gs_flatten recipe). Round-49 docstring:\nRound 49, lane L5 (2026-09-26, every figure text +1 pt): the round-42 LFIG5b library repointed to this lane (title 14 pt).'),
    ('R42 = f"{paths.FIGURE_ROOT}/ROUND42_2026-09-19"; R44 = f"{paths.FIGURE_ROOT}/ROUND44_2026-09-21"; R49 = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26"; LANE = f"{R49}/lanes/L5"; WD = f"{LANE}/scripts/wd_run.sh"; WD_LOGDIR = f"{LANE}/logs/watchdog"',
     'R42 = f"{paths.FIGURE_ROOT}/ROUND42_2026-09-19"; R44 = f"{paths.FIGURE_ROOT}/ROUND44_2026-09-21"; R49 = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26"; R52 = f"{paths.FIGURE_ROOT}/ROUND52_2026-09-26"; R54 = f"{paths.FIGURE_ROOT}/ROUND54_2026-09-28"\nLANE = f"{R54}/lanes/L5"; WD = f"{LANE}/scripts/wd_run.sh"; WD_LOGDIR = f"{LANE}/logs/watchdog"\nL5_R49 = f"{R49}/lanes/L5"; R49_VECTOR = f"{L5_R49}/Main_Fig5/work/Main_Fig5_r49_vector.pdf"; V30 = f"{R52}/figures/NEW_FINAL_SET_V30"'),
    ('SKETCH = f"{R49}/lanes/LSKETCH/work/build"; SKETCH_R44 = f"{R44}/lanes/LSKETCH/work/build"; SKETCH_STATUS = f"{R49}/lanes/LSKETCH/STATUS.md"',
     'SKETCH = f"{R54}/lanes/LSKETCH/work/build"; SKETCH_R49 = f"{R49}/lanes/LSKETCH/work/build"; SKETCH_R44 = f"{R44}/lanes/LSKETCH/work/build"\nBANDS_READY = f"{R54}/lanes/LSKETCH/BANDS_READY.txt"; SKETCH_STATUS = f"{R54}/lanes/LSKETCH/STATUS.md"'),
    ('TITLE_SIZE, TITLE_BASELINE = 14.0, 14.0          # round 49: the sheet title and the panel letters 13 -> 14 pt\nLETTER_SIZE = 14.0',
     'TITLE_SIZE, TITLE_BASELINE = 18.0, 18.0          # round 54: the sheet title and the panel letters 14 -> 18 pt (the title baseline 14 -> 18 keeps the cap top 5 pt inside the page)\nLETTER_SIZE = 18.0\nHEIGHT_AIM, HEIGHT_HARD = 1110.0, 1260.0          # round 54: the composed sheet at most 1110 pt tall (width-bound in the manuscript), 1260 the hard gate'),
    ('lines = ["# L5 STATUS (round 49)"', 'lines = ["# L5 STATUS (round 54)"'),
    ('def gs_txt(pdf, txt, name, max_s=600):',
     'def gs_flatten(src, out, name, max_s=1800):\n    """Re-distil the composed (nested) sheet into a single-page pdfwrite file for the coordinator\'s text census: the ROUND52 LSUPP-D\n    lsd_common.gs_flatten recipe (pdfwrite 1.7, /prepress, no downsampling, colours unchanged, fonts embedded), here under the watchdog."""\n    hydrated(src); os.makedirs(os.path.dirname(out), exist_ok=True)\n    if os.path.exists(out): os.remove(out)\n    rc, st = wd(name, [GS, "-dSAFER", "-dBATCH", "-dNOPAUSE", "-dQUIET", "-sDEVICE=pdfwrite", "-dCompatibilityLevel=1.7", "-dPDFSETTINGS=/prepress",\n                       "-dDownsampleColorImages=false", "-dDownsampleGrayImages=false", "-dDownsampleMonoImages=false", "-dAutoFilterColorImages=false",\n                       "-dAutoFilterGrayImages=false", "-dColorImageFilter=/FlateEncode", "-dGrayImageFilter=/FlateEncode", "-dColorConversionStrategy=/LeaveColorUnchanged",\n                       "-dEmbedAllFonts=true", "-dSubsetFonts=true", "-dPreserveAnnots=false", f"-sOutputFile={out}", src], max_s=max_s)\n    assert rc == 0 and os.path.exists(out) and os.path.getsize(out) > 0 and "error" not in wd_out(name).lower(), ("gs pdfwrite failed", name, st, wd_out(name)[-400:])\n    return st\n\n\ndef gs_txt(pdf, txt, name, max_s=600):'),
])
patch(f"{L5}/scripts/gs_text.py", [('"""Text layer of a composed (non-nested) sheet', '"""Round 54 copy (lane L5). Text layer of a composed (non-nested) sheet')])

# ---- l5a14.py: the page helper at the round-54 sizes (a role table instead of a uniform bump)
patch(f"{L5}/Main_Fig5/scripts/l5a14.py", [
    ('"""Round 49, lane L5 (2026-09-26): the round-42 copy of this library with the page helper Sheet drawing every text ONE POINT larger',
     '"""Round 54, lane L5 (2026-09-29): the round-49 copy of this library with the page helper Sheet at BUMP 0.0 and a ROLE_PT table of the\nround-54 sizes (labels, values, ticks, keys 14; column heads, axis titles, group heads 15; panel letters 18); the builder asks for the\nround-54 size by role and every drawn record carries its role (the verifier checks size against role). A uniform +4 would not do here:\nthe V13 sizes of this sheet are the generators\' sizes scaled by 0.95 (9.02, 9.5, 10.45), so +4 gives 13.02, 13.5, 14.45 instead of 14, 14, 15.\nRound-49 docstring: the round-42 copy of this library with the page helper Sheet drawing every text ONE POINT larger'),
    ('    BUMP = 1.0   # round 49: every text one point larger than the size the builder asks for (placement rules see the larger text)\n',
     '    BUMP = 0.0   # round 54: the builder asks for the final round-54 size by role (ROLE_PT), nothing is added here\n'),
    ('    def text(self, x, y, s, size, bold=False, color=INK, align="left", panel="", source_file="static:V13", source_key="", source_value="", rule="text", note=""):\n        size_asked = size; size = size + self.BUMP\n',
     '    def text(self, x, y, s, size, bold=False, color=INK, align="left", panel="", source_file="static:V13", source_key="", source_value="", rule="text", note="", role=""):\n        size_asked = size; size = size + self.BUMP\n        if role: assert abs(size - ROLE_PT[role]) < 1e-9, ("size does not match the role", s, size, role, ROLE_PT[role])\n'),
    ('rule=rule, note=note, x_left=round(x0, 3), width=round(w, 3), size_asked=size_asked))',
     'rule=rule, note=note, x_left=round(x0, 3), width=round(w, 3), size_asked=size_asked, role=role))'),
    ('INK = "#1a1d21"; INK20 = "#1a1d20";',
     '# round 54 sizes (R54_BRIEF.md): tick labels, row labels, keys, cell values 14; axis titles, column heads, group heads 15; panel letters and the title 18\nROLE_PT = {"label": 14.0, "value": 14.0, "tick": 14.0, "key": 14.0, "head": 15.0, "axis_title": 15.0, "group_head": 15.0, "letter": 18.0, "title": 18.0}\nINK = "#1a1d21"; INK20 = "#1a1d20";'),
])

# ---- 04_build_panel_a_r54.py: the flow generator with the two canvas sizes as arguments
PA = f"{L5}/Main_Fig5/scripts/04_build_panel_a_r54.py"
patch(PA, [
    ('"""\nRound 49, lane L5 (2026-09-26): the Figure 5a flow generator (FINAL_FIGURES_2026-08-14/_scripts/figure5A_flow_polish.py, the panel the V13\nsheet carries at scale S_V13 = 1.1342) rerun with every text size + BUMP where BUMP = 1/S_V13 pt on this canvas (= +1.0 pt on the sheet):\nusage 04_build_panel_a_r49.py <bump_on_canvas_pt> <tag>. Writes only into this lane (work/panel_a/). Geometry in data units unchanged, so\nthe flows, blocks and brackets (the data marks) do not move. The Aug-25 docstring follows.',
     '"""\nRound 54, lane L5 (2026-09-29): the Figure 5a flow generator (FINAL_FIGURES_2026-08-14/_scripts/figure5A_flow_polish.py, the panel the V13\nsheet carries at scale S_V13 = 1.13435) rerun with the two text sizes given on the command line IN CANVAS POINTS (the round-54 sheet sizes\n14 and 15 pt divided by S_V13: block and arm strings 12.342, captions 13.223), the arm text column X_ARM optional as a third number:\nusage 04_build_panel_a_r54.py <tick_pt_canvas> <title_pt_canvas> <tag> [x_arm]. With 10 11 bump0 it reproduces the Aug-25 panel (the\nfill remap apart). Writes only into this lane (work/panel_a/). Geometry in data units unchanged, so the flows, blocks and brackets (the\ndata marks) do not move. The round-49 and Aug-25 docstrings follow.\nRound 49, lane L5 (2026-09-26): the Figure 5a flow generator rerun with every text size + BUMP where BUMP = 1/S_V13 pt on this canvas.'),
    ('BUMP = float(sys.argv[1]) if len(sys.argv) > 1 else 0.0\nTAG = sys.argv[2] if len(sys.argv) > 2 else "bump0"\nOUT = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/L5/Main_Fig5/work/panel_a"',
     'TICK_ARG = float(sys.argv[1]) if len(sys.argv) > 1 else 10.0\nTITLE_ARG = float(sys.argv[2]) if len(sys.argv) > 2 else 11.0\nTAG = sys.argv[3] if len(sys.argv) > 3 else "bump0"\nX_ARM_ARG = float(sys.argv[4]) if len(sys.argv) > 4 else 77.5\nOUT = f"{paths.FIGURE_ROOT}/ROUND54_2026-09-28/lanes/L5/Main_Fig5/work/panel_a"'),
    ('TITLE_PT, TICK_PT, ANN_PT, FLOOR_PT = 11.0 + BUMP, 10.0 + BUMP, 9.5, 9.0     # round 49: + BUMP on this canvas',
     'TITLE_PT, TICK_PT, ANN_PT, FLOOR_PT = TITLE_ARG, TICK_ARG, 9.5, 9.0     # round 54: the canvas sizes from the command line (sheet size / S_V13)'),
    ('X_ARM = 77.5                                  # arm text, left-aligned',
     'X_ARM = X_ARM_ARG                             # arm text, left-aligned (Aug-25: 77.5; round 54 may move the text left to clear panel b, text only)'),
    ('DRAWN["bump_on_canvas_pt"] = BUMP; DRAWN["text_pt_on_canvas"] = {"blocks_arms": TICK_PT, "captions": TITLE_PT}',
     'DRAWN["bump_on_canvas_pt"] = TICK_PT - 10.0; DRAWN["text_pt_on_canvas"] = {"blocks_arms": TICK_PT, "captions": TITLE_PT}; DRAWN["x_arm"] = X_ARM'),
])
# ---- 04b map, 06b render: paths
patch(f"{L5}/Main_Fig5/scripts/04b_panel_a_map_r54.py", [
    ('"""Round 49, lane L5: (1) prove', '"""Round 54 copy (lane L5, paths only). Round 49, lane L5: (1) prove'),
    ('L5 = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/L5"', 'L5 = f"{paths.FIGURE_ROOT}/ROUND54_2026-09-28/lanes/L5"'),
])
patch(f"{L5}/Main_Fig5/scripts/06b_render_r54.py", [
    ('"""Round 49, lane L5: Ghostscript renders', '"""Round 54 copy (lane L5, paths only). Round 49, lane L5: Ghostscript renders'),
    ('sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/L5/scripts")', 'sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND54_2026-09-28/lanes/L5/scripts")'),
])
print("done")
