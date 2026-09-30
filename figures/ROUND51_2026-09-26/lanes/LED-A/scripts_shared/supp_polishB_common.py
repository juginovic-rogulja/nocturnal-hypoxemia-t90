"""ROUND 49 lane LED-A copy (2026-09-26): TITLE_PT, TICK_PT, ANN_PT + 1 pt. Original: FINAL_FIGURES_2026-08-14/_scripts/supp_polishB_common.py.

Design-quality layer for the polished supplementary batch B
(efig14_polish.py, efig15_polish.py, efig16_polish.py, efig18_polish.py,
efig22_polish.py, efig23_polish.py).

The content, data wiring, numeric assertions and text budget of the shipped
Supplementary_Figures sheets are SETTLED. This module changes only how the settled marks
are drawn, to the one design system every polished sheet shares:

  canvas   171 mm wide exactly (savefig.bbox standard), height free per sheet,
           no ink within 4 mm of any edge (asserted).
  type     Arial via splitstyle rc(): axis and panel titles 11 pt regular, ticks, row
           labels and legends 10 pt, load-bearing annotations 9.5 pt, floor 9 pt,
           bold only for the headline, panel letters and group headers.
  ink      #1a1d21 for text, spines, reference rules.
  colour   #1f4257 primary, #4d7387 secondary, #8a9099 comparison, #ccd1d6 de-emphasis,
           #aeb5bc only in its two documented roles (external cohort, pale controls),
           at most 3 data colours per panel outside the documented dose-ramp exception.
  axes     0.8 pt spines, 3 pt outward ticks, 4-6 clean majors, explicit log ticks with
           NullLocator minors, units in every axis title.
  forests  uniform row pitch, CI 1.7 pt round caps, 6 pt white-edged markers, one
           aligned label gutter, HR = 1 rules 0.9 pt dashed ink behind the data,
           optional #eef0f1 0.5 pt vertical grid on wide forests.

splitstyle.save() is never called and legend_capture is neutralized exactly as the
existing polish scripts do, so nothing writes into the frozen New_Figures tree. Output
goes only to Supplementary_Figures_Polished (PDF, same base names as the shipped set)
plus the 300 dpi PNG in _workfiles/polished_pngs and a drawn-values JSON in _workfiles.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import atexit
import json
import os
import sys

import numpy as np

sys.dont_write_bytecode = True

sys.path.insert(0, f"{paths.FIGURE_ROOT}/New_Figures/"
                   "_NOT_THESE_workfiles/scripts")
from splitstyle import (rc, save, panel_label, axis_covers, sep_ok, luma, plt,      # noqa: F401,E402
                        BLUE, BLUE_MID, GREY, GREY_PALE, NEG_CONTROLS, PANEL)
import legend_capture as _lc                                                        # noqa: E402
import nooverlap                                                                    # noqa: E402
import matplotlib.text as mtext                                                     # noqa: E402
import matplotlib.ticker as mticker                                                 # noqa: E402
assert callable(save)  # imported per the house pattern, deliberately never called here

# restore the text methods legend_capture patched and disarm its atexit flush, which
# would otherwise append into the frozen New_Figures/FIGURE_LEGENDS.md
from matplotlib.axes import Axes as _Axes                                           # noqa: E402
from matplotlib.figure import Figure as _Figure                                     # noqa: E402
_Figure.text = _lc._orig_text
_Axes.text = _lc._orig_ax_text
_Axes.annotate = _lc._orig_ax_annotate
atexit.unregister(_lc.flush)

ROOT = paths.FIGURE_ROOT
NUM = f"{paths.NUMBERS_DIR}"
OUT = f"{paths.FIGURE_ROOT}/FINAL_FIGURES_2026-08-14"
POLISH = f"{OUT}/Supplementary_Figures_Polished"
WORK = f"{OUT}/_workfiles"
PNGS = f"{WORK}/polished_pngs"
for _d in (POLISH, WORK, PNGS):
    os.makedirs(_d, exist_ok=True)

# ---------------------------------------------------------------- the shared design system
MM = 1 / 25.4
FIGW = 171 * MM                       # every polished sheet is exactly this wide
EDGE = 0.18                           # nothing inside 0.18 in (4.6 mm) of any canvas edge
INK = "#1a1d21"                       # ink for text, spines, reference rules
GRID = "#eef0f1"                      # faint forest grid and row bands
PRIMARY, SECOND, COMPARE, DEEMPH = "#1f4257", "#4d7387", "#8a9099", "#ccd1d6"
PALE2 = "#aeb5bc"                     # documented exception roles only
assert (PRIMARY, SECOND, COMPARE, DEEMPH) == (BLUE, BLUE_MID, GREY, GREY_PALE)

TITLE_PT, TICK_PT, ANN_PT, FLOOR_PT = 12.0, 11.0, 10.5, 9.0    # round 49 (lane LED-A copy): 11.0, 10.0, 9.5 + 1 pt, the 9 pt floor unchanged

CI_LW = 1.7                           # CI lines, round caps
MEDGE = 0.8                           # white marker edge
S_CIRCLE, S_SQUARE, S_DIAMOND, S_TRI = 28.0, 23.0, 26.0, 30.0   # 6 pt equal-weight areas
MS_CIRCLE, MS_SQUARE, MS_DIAMOND, MS_TRI = 6.0, 5.4, 5.7, 6.3   # legend handle diameters


def rc_polish():
    r = dict(rc())
    r.update({"font.size": TICK_PT, "axes.labelsize": TITLE_PT, "axes.titlesize": TITLE_PT,
              "xtick.labelsize": TICK_PT, "ytick.labelsize": TICK_PT,
              "legend.fontsize": TICK_PT,
              "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
              "xtick.major.size": 3.0, "ytick.major.size": 3.0,
              "axes.edgecolor": INK, "text.color": INK, "axes.labelcolor": INK,
              "xtick.color": INK, "ytick.color": INK,
              "savefig.bbox": "standard"})
    return r


def headline(fig, text, y_in=None):
    """The sheet's one bold sentence, top left, inside the 4 mm margin."""
    W, H = fig.get_size_inches()
    y = 1.0 - (EDGE if y_in is None else y_in) / H
    return fig.text(EDGE / W, y, text, fontsize=TITLE_PT,
                    fontweight="bold", ha="left", va="top", color=INK, linespacing=1.3)


def letter(fig, s, x_in, y_in):
    """A bold panel letter placed in inches from the sheet's top left corner."""
    W, H = fig.get_size_inches()
    fig.text(x_in / W, 1.0 - y_in / H, s, fontsize=PANEL, fontweight="bold",
             ha="left", va="top", color=INK)


def logx_clean(ax, ticks):
    """Explicit log majors, NullLocator minors. The only log-axis form in this package."""
    ax.set_xscale("log")
    ax.set_xticks(ticks)
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.xaxis.set_minor_locator(mticker.NullLocator())


def logy_clean(ax, ticks):
    ax.set_yscale("log")
    ax.set_yticks(ticks)
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.yaxis.set_minor_locator(mticker.NullLocator())


def bare_y(ax):
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)


def ref_rule(ax, x=1.0):
    """The HR = 1 reference: 0.9 pt dashed ink, behind the data."""
    ax.axvline(x, color=INK, lw=0.9, ls=(0, (4, 3)), zorder=1)


def vgrid(ax, ticks):
    """Faint vertical grid at the tick positions on wide forests, behind everything."""
    for t in ticks:
        ax.axvline(t, color=GRID, lw=0.5, zorder=0.3)


def text_w(fig, s, size=TICK_PT, weight="normal"):
    """Rendered width of a string, in inches, at the sheet's real dpi."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    t = fig.text(0.5, 0.5, s, fontsize=size, fontweight=weight)
    w = t.get_window_extent(renderer=r).width / fig.dpi
    t.remove()
    return w


def solid_tint(hexcol, frac):
    """The solid colour equal to hexcol drawn at alpha=frac over white. PDFs record raw
    fills, so translucent bands fail the delivered-file palette audit: draw the mix."""
    h = hexcol.lstrip("#")
    rgb = [int(h[i:i + 2], 16) for i in (0, 2, 4)]
    mixed = [round(frac * c + (1 - frac) * 255) for c in rgb]
    return "#" + "".join(f"{c:02x}" for c in mixed)


# ---------------------------------------------------------------- the gates
def _dash_ok(s):
    """En dash is legal only inside a numeric range (0–1%, 1.19–1.44)."""
    for i, ch in enumerate(s):
        if ch == "–":
            if not (0 < i < len(s) - 1 and s[i - 1].isdigit() and s[i + 1].isdigit()):
                return False
    return True


def text_gate(fig, name, bold_ok=()):
    """Banned punctuation, the 9 pt floor, and bold only where the sheet licenses it."""
    bold_ok = set(bold_ok)
    for t in fig.findobj(mtext.Text):
        s = t.get_text()
        if not s or not s.strip() or not t.get_visible():
            continue
        for bad, nm in (("—", "em dash"), (";", "semicolon"),
                        ("×", "multiplication sign")):
            assert bad not in s, f"{name}: banned {nm} in on-sheet text: {s!r}"
        assert _dash_ok(s), f"{name}: en dash outside a numeric range: {s!r}"
        assert float(t.get_fontsize()) >= FLOOR_PT, (
            f"{name}: {t.get_fontsize()} pt under the {FLOOR_PT} pt floor: {s!r}")
        if str(t.get_fontweight()) in ("bold", "700"):
            assert s in bold_ok, f"{name}: bold outside the licensed set: {s!r}"


def margin_gate(fig, name, min_mm=4.0):
    """No ink within 4 mm of any sheet edge, asserted on every rendered text box and on
    the figure's tight bounding box."""
    boxes, _ = nooverlap._boxes(fig)
    Wpx, Hpx = fig.get_size_inches() * fig.dpi
    m = (min_mm - 0.3) * MM * fig.dpi          # 0.3 mm of numeric slack
    tight = [f"{b['where']}: {b['text'][:44]}" for b in boxes
             if b["bbox"].x0 < m or b["bbox"].y0 < m
             or b["bbox"].x1 > Wpx - m or b["bbox"].y1 > Hpx - m]
    assert not tight, f"{name}: text inside the {min_mm} mm edge margin: {tight}"
    bb = fig.get_tightbbox(fig.canvas.get_renderer())
    W, H = fig.get_size_inches()
    mi = (min_mm - 0.5) * MM
    bad = []
    if bb.x0 < mi:
        bad.append(f"left {bb.x0 / MM:.1f} mm")
    if bb.y0 < mi:
        bad.append(f"bottom {bb.y0 / MM:.1f} mm")
    if W - bb.x1 < mi:
        bad.append(f"right {(W - bb.x1) / MM:.1f} mm")
    if H - bb.y1 < mi:
        bad.append(f"top {(H - bb.y1) / MM:.1f} mm")
    assert not bad, f"{name}: ink closer than {min_mm} mm to an edge: {', '.join(bad)}"


def save_polished(fig, name, drawn, bold_ok=()):
    """All gates, then PDF (same base name as the shipped sheet) + 300 dpi PNG + JSON."""
    W = fig.get_size_inches()[0]
    assert abs(W - FIGW) < 1e-6, f"{name}: sheet is {W:.3f} in wide, not 171 mm"
    text_gate(fig, name, bold_ok)
    nooverlap.gate(fig, name)
    margin_gate(fig, name)
    pdf = f"{POLISH}/{name}.pdf"
    fig.savefig(pdf)
    fig.savefig(f"{PNGS}/{name}.png", dpi=300)
    plt.close(fig)
    json.dump(drawn, open(f"{WORK}/{name}_polish_drawn_values.json", "w"), indent=1,
              default=float)
    for p in (pdf, f"{PNGS}/{name}.png", f"{WORK}/{name}_polish_drawn_values.json"):
        assert os.path.exists(p) and os.path.getsize(p) > 0, p
        print("  wrote", p, f"({os.path.getsize(p):,} bytes)")
    return pdf


def endash(s):
    """Numeric ranges take the manuscript-standard en dash. Only a hyphen with digits
    on both sides is touched."""
    out = []
    for i, ch in enumerate(s):
        if ch == "-" and 0 < i < len(s) - 1 and s[i - 1].isdigit() and s[i + 1].isdigit():
            out.append("–")
        else:
            out.append(ch)
    return "".join(out)
