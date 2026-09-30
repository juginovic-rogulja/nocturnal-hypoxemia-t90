"""
Design-quality layer for the polished supplementary sheets
(eFigure05_polish.py, eFigure07_polish.py, eFigure10_polish.py, eFigure11_polish.py).

The content, data wiring, printed values and text budget of the settled redesigns
(make_eFigure05/07/10/11) are carried over verbatim. This layer changes only how the
settled marks are drawn, to the one design system the polished mains already hold:

  canvas   171 mm wide exactly (savefig.bbox standard, no tight crop), height free for
           supplementary sheets, no ink within 4 mm of any edge (margin_gate).
  type     Arial via splitstyle rc(): axis titles 11 pt regular, panel titles 10.5 pt
           regular, ticks, row labels and legends 10 pt, load-bearing annotations
           9.5 pt, floor 9 pt. Bold only for the headline and the panel letters.
  colour   ink #1a1d21 for text, spines and reference rules. Data colours stay in the
           house roles: #1f4257 primary, #4d7387 secondary, #8a9099 comparison,
           #ccd1d6 de-emphasis, bronze only where the redesign already used it.
           Ordered ramps run the declared house ramps, solid fills, never alpha.
  axes     left+bottom spines 0.8 pt, 3 pt outward ticks, 4 to 6 clean ticks,
           explicit log ticks with NullLocator minors, units in every axis title.
  marks    CI lines 1.7 pt round caps, 6 pt markers with a 0.8 pt white edge,
           HR = 1 rules 0.9 pt dashed ink behind the data, frameless 10 pt legends.

splitstyle.save() is never called and legend_capture is disarmed exactly as in the
existing polish scripts, so nothing writes into the frozen New_Figures tree. Output
goes only to Supplementary_Figures_Polished (PDF) plus the 300 dpi PNG in
_workfiles/polished_pngs and a drawn-values JSON in _workfiles.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import atexit
import json
import os
import sys

sys.dont_write_bytecode = True

SCRIPTS = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LSUPP-C/L5b_scripts"   # V14_L5b 2026-09-15: the lane copies of splitstyle, legend_capture, nooverlap
sys.path.insert(0, SCRIPTS)

from splitstyle import (rc, axis_covers, sep_ok, BLUE, GREY, GREY_PALE, BLUE_MID,  # noqa: E402,F401
                        BRONZE, plt)
import legend_capture as _lc                                                       # noqa: E402
import nooverlap                                                                   # noqa: E402
from prism_plotter.journal import preflight_figure, journal_style                  # noqa: E402
from matplotlib.ticker import NullLocator, NullFormatter, FuncFormatter            # noqa: E402
import matplotlib.text as _mtext                                                   # noqa: E402

# restore the text methods legend_capture patched and disarm its atexit flush, which
# would otherwise append into the frozen New_Figures/FIGURE_LEGENDS.md
from matplotlib.axes import Axes as _Axes                                          # noqa: E402
from matplotlib.figure import Figure as _Figure                                    # noqa: E402
_Figure.text = _lc._orig_text
_Axes.text = _lc._orig_ax_text
_Axes.annotate = _lc._orig_ax_annotate
atexit.unregister(_lc.flush)

ROOT = paths.FIGURE_ROOT
NUM = f"{paths.NUMBERS_DIR}"
OUT = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LSUPP-C/L5b_work/polish_out"   # V14_L5b: nothing is written under FINAL_FIGURES
POLISH = f"{OUT}/pdf"
WORK = f"{OUT}/work"
PNGS = f"{OUT}/png"
for _d in (POLISH, WORK, PNGS):
    os.makedirs(_d, exist_ok=True)

# ---------------------------------------------------------------- the design system
MM = 1 / 25.4
W_IN = 171.0 * MM                 # every polished sheet is exactly 171 mm wide
EDGE = 0.18                       # nothing inside 0.18 in (4.6 mm) of any canvas edge

INK_P = "#1a1d21"                 # the design system's ink, replaces splitstyle's #15181c
GRID = "#eef0f1"
BAND_FAINT = "#eef0f1"            # large shaded regions (negative-control block)
BAND_PALE = "#dfe3e6"             # narrow shaded regions that must read under a scatter
PRIMARY, SECOND, COMPARE, DEEMPH = "#0288d1", "#3f9fd8", "#8a9099", "#ccd1d6"   # V14_L5b: the V13 palette (the sheets were recoloured to the main-figure blues after the August build)
ROLE_RAMP4 = [DEEMPH, COMPARE, SECOND, PRIMARY]          # declared house 4-step ramp
FIFTHS_RAMP5 = ["#e2e7ea", "#b6c1c8", "#8a9aa5", "#4d6577", PRIMARY]   # declared 5-step

PT_PLUS = 1.0                     # R49 (Alen, 2026-09-25): every text one point larger, geometry unchanged
HEAD_PT = 11.5 + PT_PLUS          # the headline, bold
TITLE_PT = 11.0 + PT_PLUS         # axis titles, regular
PTITLE_PT = 10.5 + PT_PLUS        # panel titles, regular
TICK_PT = 10.0 + PT_PLUS          # ticks, row labels, legends
ANN_PT = 9.5 + PT_PLUS            # load-bearing annotations, printed columns
FLOOR_PT = 9.0
PANEL_PT = 13.0 + PT_PLUS         # bold panel letters

CI_LW = 1.7                       # CI lines, round caps
MEDGE = 0.8                       # white marker edge
S_CIRCLE = 28.0                   # scatter area of a 6 pt circle
S_SQUARE = 23.0                   # an equal-weight square
MS_CIRCLE = 6.0                   # legend handle diameters, points
MS_SQUARE = 5.4

JSTYLE = journal_style("jama", columns=2)


def rc_polish():
    r = dict(rc())
    r.update({
        "savefig.bbox": "standard",          # fixed canvas, 171 mm stays 171 mm
        "font.size": TICK_PT, "axes.labelsize": TITLE_PT, "axes.titlesize": PTITLE_PT,
        "xtick.labelsize": TICK_PT, "ytick.labelsize": TICK_PT,
        "legend.fontsize": TICK_PT,
        "axes.linewidth": 0.8, "axes.edgecolor": INK_P,
        "text.color": INK_P, "axes.labelcolor": INK_P,
        "xtick.color": INK_P, "ytick.color": INK_P,
        "xtick.major.width": 0.8, "ytick.major.width": 0.8,
        "xtick.major.size": 3.0, "ytick.major.size": 3.0,
    })
    return r


# ---------------------------------------------------------------- drawing helpers
def ref_rule(ax, x=1.0):
    """The reference rule: 0.9 pt dashed ink, behind the data."""
    ax.axvline(x, color=INK_P, lw=0.9, ls=(0, (4, 3)), zorder=1)


def bare_y(ax):
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)


def logx_clean(ax, ticks):
    """Log x axis with explicit majors and no minors, the package convention."""
    ax.set_xscale("log")
    ax.set_xticks(ticks)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _p: f"{v:g}"))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_minor_formatter(NullFormatter())


def stars(p):
    import numpy as np
    if p is None or not np.isfinite(p):
        return ""
    return "***" if p < .001 else ("**" if p < .01 else ("*" if p < .05 else ""))


def measure(fig, s, size, weight="normal"):
    """Rendered width of a string in inches, on this figure's renderer."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    t = fig.text(0.5, 0.5, s, fontsize=size, fontweight=weight)
    w = t.get_window_extent(renderer=r).width / fig.dpi
    t.remove()
    return w


def fit_gutter(fig, texts, x0_in, limit_in, size=TICK_PT, weight="normal", what=""):
    """Assert every gutter string ends before the plot edge it must not cross."""
    for s in texts:
        w = measure(fig, s, size, weight)
        assert x0_in + w <= limit_in + 1e-6, (
            f"{what}: gutter label {s!r} is {w:.2f} in wide and would cross the plot "
            f"edge at {limit_in:.2f} in (starts {x0_in:.2f})")


def headline(fig, lines, h_in):
    """The bold headline, top left, EDGE margins, one string per line."""
    for i, s in enumerate(lines):
        fig.text(EDGE / W_IN, (h_in - EDGE - i * HEAD_PT * 1.30 / 72.0) / h_in, s,
                 fontsize=HEAD_PT, fontweight="bold", ha="left", va="top", color=INK_P)


def panel_letter(fig, letter, x_in, y_in, h_in):
    """Bold panel letter at an absolute position, va top."""
    fig.text(x_in / W_IN, y_in / h_in, letter, fontsize=PANEL_PT, fontweight="bold",
             ha="left", va="top", color=INK_P)


def head_block(fig, letter, title, x_in, y_in, h_in, size=TITLE_PT, lift_lines=None):
    """Bold panel letter above a regular-weight title, anchored at the panel gutter
    (make_efig08_09_polish.head_block pattern). y_in is the title BOTTOM in inches.
    lift_lines pins the letter height across panels whose titles differ in line count."""
    lines = title.count("\n") + 1 if title else 0
    n = lines if lift_lines is None else lift_lines
    lift = 0.08 + n * size * 1.30 / 72.0
    fig.text(x_in / W_IN, (y_in + lift) / h_in, letter, fontsize=PANEL_PT,
             fontweight="bold", ha="left", va="bottom", color=INK_P)
    if title:
        fig.text(x_in / W_IN, y_in / h_in, title, fontsize=size, ha="left", va="bottom",
                 color=INK_P, linespacing=1.25)


def left_rows(ax, ys, labels, gutter_x_in, ax_x0_in, size=TICK_PT):
    """Row labels left-aligned in the shared gutter, as 10 pt ink tick labels."""
    ax.set_yticks(list(ys))
    ax.set_yticklabels(labels, fontsize=size, color=INK_P, ha="left")
    ax.tick_params(axis="y", pad=(ax_x0_in - gutter_x_in) * 72.0, length=0)


# ---------------------------------------------------------------- gates and writer
def _prose_gate(fig, name, headline_lines):
    """No em-dash, semicolon or multiplication sign, en-dash only inside a numeric
    range, nothing under the 9 pt floor, bold only for the headline and panel letters."""
    allowed_bold = set(headline_lines) | set("ABCDE")
    for t in fig.findobj(_mtext.Text):
        s = t.get_text()
        if not s or not s.strip() or not t.get_visible():
            continue
        for bad, nm in (("—", "em dash"), (";", "semicolon"), ("×", "multiplication sign")):
            assert bad not in s, f"{name}: banned {nm} in on-sheet text: {s!r}"
        for i, ch in enumerate(s):
            if ch == "–":
                assert 0 < i < len(s) - 1 and s[i - 1].isdigit() and s[i + 1].isdigit(), (
                    f"{name}: en dash outside a numeric range: {s!r}")
        if str(t.get_fontweight()) in ("bold", "700"):
            assert s in allowed_bold, f"{name}: bold outside headline/panel letters: {s!r}"
        assert float(t.get_fontsize()) >= FLOOR_PT, (
            f"{name}: {t.get_fontsize()} pt under the {FLOOR_PT} pt floor: {s!r}")


def margin_gate(fig, name, min_mm=4.0):
    """No ink within 4 mm of any sheet edge, asserted on the rendered extent."""
    fig.canvas.draw()
    bb = fig.get_tightbbox(fig.canvas.get_renderer())
    W, H = fig.get_size_inches()
    m = min_mm * MM
    bad = []
    if bb.x0 < m:
        bad.append(f"left {bb.x0 / MM:.1f} mm")
    if bb.y0 < m:
        bad.append(f"bottom {bb.y0 / MM:.1f} mm")
    if W - bb.x1 < m:
        bad.append(f"right {(W - bb.x1) / MM:.1f} mm")
    if H - bb.y1 < m:
        bad.append(f"top {(H - bb.y1) / MM:.1f} mm")
    assert not bad, f"{name}: content closer than {min_mm} mm to an edge: {', '.join(bad)}"


def save_polished(fig, name, drawn, headline_lines=(), min_mm=4.0):
    """Prose gate, overlap gate, margin gate, journal preflight (height free for
    supplementary sheets), then PDF into Supplementary_Figures_Polished and the
    300 dpi PNG into _workfiles/polished_pngs. The only writer in this package."""
    _prose_gate(fig, name, headline_lines)
    nooverlap.gate(fig, name)
    margin_gate(fig, name, min_mm=min_mm)
    problems = preflight_figure(fig, JSTYLE)
    real = [p for p in problems if "height" not in p.lower()]
    assert not real, f"{name}: {real}"

    pdf = f"{POLISH}/{name}.pdf"
    fig.savefig(pdf)
    fig.savefig(f"{PNGS}/{name}.png", dpi=300)
    plt.close(fig)
    json.dump(drawn, open(f"{WORK}/{name}_polished_drawn_values.json", "w"), indent=1,
              default=float)
    for p in (pdf, f"{PNGS}/{name}.png", f"{WORK}/{name}_polished_drawn_values.json"):
        assert os.path.exists(p) and os.path.getsize(p) > 0, p
        print("wrote", p)
    return pdf
