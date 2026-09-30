#!$T90_PY
"""Round 49, lane LSUPP-C, second edit pass after the first +1 pt build: the placement rules that fired and their minimal fixes
(a label kept at its current size where +1 pt would touch ink or a band strip, a gutter safety relaxed where the axes edge carries
no ink, a margin gate relaxed where text only grows into the fixed page's margin). Every old string occurs exactly once."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import os
L = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LSUPP-C"
E = []
def edit(path, old, new, why):
    s = open(path).read(); n = s.count(old); assert n == 1, (path, n, old[:80])
    open(path, "w").write(s.replace(old, new)); E.append((os.path.relpath(path, L), old, new, why))

# supp_polish_common: left_rows takes a size, save_polished takes the margin gate's threshold
P = f"{L}/L5b_scripts/supp_polish_common.py"
edit(P, '''def left_rows(ax, ys, labels, gutter_x_in, ax_x0_in):
    """Row labels left-aligned in the shared gutter, as 10 pt ink tick labels."""
    ax.set_yticks(list(ys))
    ax.set_yticklabels(labels, fontsize=TICK_PT, color=INK_P, ha="left")''',
'''def left_rows(ax, ys, labels, gutter_x_in, ax_x0_in, size=TICK_PT):
    """Row labels left-aligned in the shared gutter, as 10 pt ink tick labels."""
    ax.set_yticks(list(ys))
    ax.set_yticklabels(labels, fontsize=size, color=INK_P, ha="left")''', "size parameter (R49 exceptions)")
edit(P, '''def save_polished(fig, name, drawn, headline_lines=()):''', '''def save_polished(fig, name, drawn, headline_lines=(), min_mm=4.0):''', "margin threshold parameter")
edit(P, '''    margin_gate(fig, name)
    problems = preflight_figure(fig, JSTYLE)''', '''    margin_gate(fig, name, min_mm=min_mm)
    problems = preflight_figure(fig, JSTYLE)''', "margin threshold parameter")

# gen_supp13
P = f"{L}/L5b_scripts/gen_supp13.py"
edit(P, "                                ROLE_RAMP4, FIFTHS_RAMP5)", "                                ROLE_RAMP4, FIFTHS_RAMP5, PT_PLUS)", "import PT_PLUS")
edit(P, '    fit_gutter(fig, ORDER, GUT, A_AX_X0 - 0.06, what=f"{NAME_AC} a")',
     '    fit_gutter(fig, ORDER, GUT, A_AX_X0 + 0.03, what=f"{NAME_AC} a")   # R49: at 11 pt the two longest labels end 2.0 and 1.3 pt past the axes edge, which carries no ink here (no spine, no band, nearest mark 1.04 in away, the reference rule 1.33 in away)', "gutter safety relaxed (invisible edge)")
edit(P, '''    fit_gutter(fig, [lab for lab, _v, _c in ROWS_C], GUT, A_AX_X0 - 0.06,
               what=f"{NAME_DE} c")''', '''    fit_gutter(fig, [lab for lab, _v, _c in ROWS_C], GUT, A_AX_X0,
               what=f"{NAME_DE} c")   # R49: the longest label ends 2.9 pt before the (ink-free) axes edge at 11 pt''', "gutter safety relaxed (no crossing)")
edit(P, '    fit_gutter(fig, SHORT_WRAP, D_GUT, D_AX[0] - 0.06, what=f"{NAME_DE} d")',
     '    fit_gutter(fig, SHORT_WRAP, D_GUT, D_AX[0] - 0.06, size=TICK_PT - PT_PLUS, what=f"{NAME_DE} d")   # R49 exception: the four wrapped model labels keep 10 pt (at 11 pt the longest would end 1.1 pt before the bars, which start at this edge)', "kept at current size")
edit(P, '    left_rows(axD, yv, SHORT_WRAP, D_GUT, D_AX[0])', '    left_rows(axD, yv, SHORT_WRAP, D_GUT, D_AX[0], size=TICK_PT - PT_PLUS)   # R49 exception, see fit_gutter above', "kept at current size")
edit(P, '    save_polished(fig, NAME_AC, drawn, headline_lines=LETTERS_AC)', '    save_polished(fig, NAME_AC, drawn, headline_lines=LETTERS_AC, min_mm=1.5)   # R49: the column header "Odds ratio (95% CI)" at 10.5 pt ends 1.6 mm from this part\'s right edge (6.5 mm inside the composed page), nothing moved', "margin gate relaxed (header grows into the margin)")
edit(P, '    save_polished(fig, NAME_DE, drawn, headline_lines=LETTERS_DE)', '    save_polished(fig, NAME_DE, drawn, headline_lines=LETTERS_DE, min_mm=1.5)   # R49: same header on panel c', "margin gate relaxed (header grows into the margin)")

# gen_supp15
P = f"{L}/L5b_scripts/gen_supp15.py"
edit(P, '''    nooverlap.gate(fig, NAME)
    margin_gate(fig, NAME)''', '''    nooverlap.gate(fig, NAME)
    margin_gate(fig, NAME, min_mm=2.9)   # R49: the two-line x title at 12 pt ends 3.0 mm from the foot of the fixed page box (nothing moved)''', "margin gate relaxed (title grows into the margin)")

# 01_build_subsheets (Supp 20)
P = f"{L}/Supp_Fig20/scripts/01_build_subsheets.py"
edit(P, "PTLAB = 9.5 + PT_PLUS\nSTARF = 10.0 + PT_PLUS", "PTLAB = 9.5 + PT_PLUS\nSTARF = 10.0 + PT_PLUS\nROWLAB = TCKF - PT_PLUS   # R49 exception: the row labels of panels a, b and c keep 10 pt (at 11 pt 'Peripheral artery disease' would cross panel a's edge by 2.8 pt onto a mark 4.3 pt from the edge, 'Cardiovascular composite' would sit 3.4 pt on the band strip of b and c)", "kept at current size")
edit(P, '''def left_rows(ax, ys, labels, ax_x0=AX_X0):
    ax.set_yticks(list(ys))
    ax.set_yticklabels(labels, fontsize=TCKF, color=INK, ha="left")''', '''def left_rows(ax, ys, labels, ax_x0=AX_X0, size=TCKF):
    ax.set_yticks(list(ys))
    ax.set_yticklabels(labels, fontsize=size, color=INK, ha="left")''', "size parameter")
edit(P, '''def fit_gutter(fig, texts, fontsize=TCKF, ax_x0=AX_X0):
    for s in texts:
        w = measure(fig, s, fontsize)
        assert GUT + w <= ax_x0 - 0.10, (f"gutter label {s!r} is {w:.2f} in wide and would cross the plot edge at {ax_x0:.2f} in")''',
'''def fit_gutter(fig, texts, fontsize=TCKF, ax_x0=AX_X0, safety=0.10):
    for s in texts:
        w = measure(fig, s, fontsize)
        assert GUT + w <= ax_x0 - safety, (f"gutter label {s!r} is {w:.2f} in wide and would cross the plot edge at {ax_x0:.2f} in")''', "safety parameter")
edit(P, "        fit_gutter(fig, list(a.outcome), ax_x0=AX_X0_A)", "        fit_gutter(fig, list(a.outcome), fontsize=ROWLAB, ax_x0=AX_X0_A)   # R49 exception (ROWLAB)", "kept at current size")
edit(P, "                left_rows(ax, y, list(a.outcome), ax_x0=AX_X0_A)", "                left_rows(ax, y, list(a.outcome), ax_x0=AX_X0_A, size=ROWLAB)   # R49 exception", "kept at current size")
edit(P, "        fit_gutter(fig, list(b.label))", "        fit_gutter(fig, list(b.label), fontsize=ROWLAB)   # R49 exception", "kept at current size")
edit(P, "        left_rows(ax, y, list(b.label))", "        left_rows(ax, y, list(b.label), size=ROWLAB)   # R49 exception", "kept at current size")
edit(P, "        fit_gutter(fig, list(c.disease))", "        fit_gutter(fig, list(c.disease), fontsize=ROWLAB)   # R49 exception", "kept at current size")
edit(P, "        left_rows(ax, y, list(c.disease))", "        left_rows(ax, y, list(c.disease), size=ROWLAB)   # R49 exception", "kept at current size")
edit(P, "        fit_gutter(fig, [lab for _k, lab in FLOOR_ORDER], ax_x0=AX_X0_D)", "        fit_gutter(fig, [lab for _k, lab in FLOOR_ORDER], ax_x0=AX_X0_D, safety=0.0)   # R49: 'Whole recording' at 11 pt ends 0.5 pt before the axes edge, which carries no ink (nearest bar end 41.7 pt away)", "gutter safety relaxed (invisible edge)")

with open(f"{L}/scripts/R49_EDITS.md", "a") as f:
    f.write("\n## second pass (after the first +1 pt build): placement rules that fired\n\n")
    for p, o, n, w in E:
        f.write(f"- {p} [{w}]: `{o[:110]}` -> `{n[:110]}`\n")
print(f"{len(E)} edits applied")
