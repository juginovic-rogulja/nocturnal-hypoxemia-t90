#!$T90_PY
"""Round 54 compose of Main_Fig6: the round-54 triangle sheet (work/build/Main_Fig6_tri_r54.pdf, notes 15, captions 14, triangle labels 14)
through the round-44b compose function with the title 'Figure 6' at 18 pt (round 49: 14) in a 20-pt title strip (was 16 pt) with the
baseline at 18 pt (was 14), so the taller title keeps about 5 pt of clearance from the page's top edge. Figure 6 has no panel letters.
The page width is unchanged (968.94 pt), the height is the build's 283.464 + 20 = 303.464 pt (V30: 299.464), declared. Deliverable
ROUND54/lanes/LSKETCH/Main_Fig6/ (Main_Fig6.pdf, Main_Fig6_150dpi.png by Ghostscript, verify/checks.txt continued by verify_fig6_r54.py)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import inspect, sys
sys.dont_write_bytecode = True
T90 = paths.FIGURE_ROOT; R54 = f"{paths.FIGURE_ROOT}/ROUND54_2026-09-28"; V30 = f"{paths.FIGURE_ROOT}/ROUND52_2026-09-26/figures/NEW_FINAL_SET_V30"
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND44_2026-09-21/lanes/LSKETCH/scripts"); import compose_r44b as CR
CR.LANE = f"{R54}/lanes/LSKETCH"; CR.BUILD = f"{R54}/lanes/LSKETCH/work/build"; CR.L.TITLE_SIZE = 18.0
TOP_PAD, TITLE_BASELINE = 20.0, 18.0
src = inspect.getsource(CR.fig6)
for old, new in [("TOP_PAD = 16.0; TITLE_XY = (36.0, L.TITLE_BASELINE)", f"TOP_PAD = {TOP_PAD}; TITLE_XY = (36.0, {TITLE_BASELINE})"),
                 ('f"page {W} x {H} = build + 16 pt strip, title \'Figure 6\'"', 'f"page {W} x {H} = build + 20 pt strip (was 16), title \'Figure 6\' at 18 pt on the 18-pt baseline (was 14 on 14)"'),
                 ('if SHEET == "Main_Fig6": same_size_check(SHEET, OUT, checks)', 'if SHEET == "Main_Fig6": w30, h30 = pagesize(f"{V30}/{SHEET}.pdf"); checks.append(("PASS " if abs(W - w30) < 0.05 else "FAIL ") + f"page width {W:.2f} = V30\'s {w30:.2f}, height {H:.3f} declared (V30 {h30:.3f}, the 4-pt taller title strip)")')]:
    assert src.count(old) == 1, old
    src = src.replace(old, new)
CR.V30 = V30; exec(src, CR.__dict__)   # fig6 with the round-54 strip, in the module's namespace
ok = CR.fig6(f"{CR.BUILD}/Main_Fig6_tri_r54.pdf", "Main_Fig6")
sys.exit(0 if ok else 1)
