"""
Design-quality layer for the five polished Figure 2 panels (figure2A_polish.py ...
figure2E_polish.py).

The content, data wiring and the decluttered text budget of the split sheets are SETTLED.
This module changes only how the settled marks are drawn, to one design system shared by
every panel:

  ink #1a1d21 for text, spines and reference rules, 0.8 pt spines, 3 pt outward ticks,
  Arial with axis titles 11 pt regular, ticks and row labels 10 pt, legends 10 pt,
  load-bearing annotations 9.5 pt, bold reserved for 10 pt group headers,
  CI lines 1.7 pt with round caps, 6 pt markers with a 0.8 pt white edge,
  HR = 1 rules 0.9 pt dashed ink behind the data, a faint #eef0f1 vertical grid at the
  tick positions on every forest, and one shared left plot-area edge (2.12 in) with one
  shared left-aligned label gutter for the three forest panels.

Data plumbing, palette roles, the legend_capture neutralisation and the frozen-tree
discipline all come from figure2_split_common, which is imported (not modified). Output
goes to Main_Figures_Polished only. The split sheets and their scripts are not touched.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import json

import figure2_split_common as SC
from figure2_split_common import (rc_split, die, bare_y, logx_clean, logy_clean,   # noqa: F401
                                  axis_covers, BAND, BAND_LABELS, BAND_MARK, BLUE,
                                  BLUE_MID, GREY, PALE, plt, NEG_CONTROLS,
                                  FLOOR_PER_SD, FLOOR_SET_BY, NUM, SV, W_IN, WORK,
                                  BASE, TICKF, SMALLF, JSTYLE, nooverlap,
                                  preflight_figure, _prose_gate)

import os

INK = "#1a1d21"                  # the design system's ink, replaces splitstyle's #15181c
GRID = "#eef0f1"                 # faint forest grid, behind everything

# 2026-08-14 subfolder layout: the polished delivery files PDFs under per-figure
# subfolders (Main_Figures_Polished/Figure2/) and rasters under _workfiles/polished_pngs,
# matching the folder state Alen was shown. Same files, new addresses.
OUT = f"{paths.FIGURE_ROOT}/FINAL_FIGURES_2026-08-14/Main_Figures_Polished/Figure2"
PNGS = f"{SC.WORK}/polished_pngs"
os.makedirs(OUT, exist_ok=True)
os.makedirs(PNGS, exist_ok=True)

# ---------------------------------------------------------------- shared geometry, inches
EDGE = 0.18                      # nothing inside 0.18 in (4.6 mm) of any canvas edge
GUT_HEAD_X = 0.18                # group headers, bold, far left
GUT_COND_X = 0.38                # condition names, left-aligned, one shared gutter
AX_X0 = 2.12                     # the one left plot-area edge for every forest panel
AX_RIGHT = 0.18
FOREST_W = W_IN - AX_X0 - AX_RIGHT
COND_PAD_PT = (AX_X0 - GUT_COND_X) * 72.0    # tick pad that left-aligns names at GUT_COND_X

# ---------------------------------------------------------------- shared mark weights
CI_LW = 1.7                      # CI lines, round caps
MEDGE = 0.8                      # white marker edge
S_CIRCLE = 28.0                  # scatter areas for a 6 pt circle ...
S_SQUARE = 23.0                  # ... an equal-weight square ...
S_TRI = 30.0                     # ... and an equal-weight triangle
MS_CIRCLE = 6.0                  # legend handle diameters, points
MS_SQUARE = 5.4


def rc_polish():
    r = dict(rc_split())
    r.update({
        "axes.linewidth": 0.8, "axes.edgecolor": INK,
        "text.color": INK, "axes.labelcolor": INK,
        "xtick.color": INK, "ytick.color": INK,
        "xtick.major.width": 0.8, "ytick.major.width": 0.8,
        "xtick.major.size": 3.0, "ytick.major.size": 3.0,
    })
    return r


def ref_rule(ax, x=1.0):
    """The HR = 1 reference: 0.9 pt dashed ink, behind the data."""
    ax.axvline(x, color=INK, lw=0.9, ls=(0, (4, 3)), zorder=1)


def ref_rule_h(ax, y=1.0):
    ax.axhline(y, color=INK, lw=0.9, ls=(0, (4, 3)), zorder=1)


def vgrid(ax, ticks):
    """Faint vertical grid at the tick positions, identical on every forest panel."""
    for t in ticks:
        ax.axvline(t, color=GRID, lw=0.5, zorder=0.3)


def left_labels(ax, ys, labels):
    """Row labels left-aligned in the shared gutter, as tick labels (10 pt ink)."""
    ax.set_yticks(list(ys))
    ax.set_yticklabels(labels, fontsize=TICKF, color=INK, ha="left")
    ax.tick_params(axis="y", pad=COND_PAD_PT, length=0)


def fit_gutter(fig, texts, fontsize=TICKF, weight="normal", x0=GUT_COND_X, limit=AX_X0):
    """Assert every gutter string actually ends before the shared plot edge."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    for s in texts:
        t = fig.text(0.5, 0.5, s, fontsize=fontsize, fontweight=weight)
        w = t.get_window_extent(renderer=r).width / fig.dpi
        t.remove()
        assert x0 + w <= limit - 0.05, (
            f"gutter label {s!r} is {w:.2f} in wide and would cross the shared "
            f"plot edge at {limit:.2f} in")


def save_polished(fig, name, drawn):
    """The split gates (prose, overlap, journal preflight) plus the 4 mm edge rule,
    then PDF + 300 dpi PNG into Main_Figures_Polished only."""
    _prose_gate(fig, name)

    # canvas is fixed (savefig.bbox standard): every string keeps 4 mm of breathing room
    boxes, _ = nooverlap._boxes(fig)
    Wpx, Hpx = fig.get_size_inches() * fig.dpi
    m = (4.0 - 0.3) / 25.4 * fig.dpi          # 4 mm with 0.3 mm of numeric slack
    _tight = [f"{b['where']}: {b['text'][:44]}" for b in boxes
              if b["bbox"].x0 < m or b["bbox"].y0 < m
              or b["bbox"].x1 > Wpx - m or b["bbox"].y1 > Hpx - m]
    assert not _tight, f"{name}: text inside the 4 mm edge margin: {_tight}"

    nooverlap.gate(fig, name)
    problems = preflight_figure(fig, JSTYLE)
    real = [p for p in problems if "height" not in p.lower()]
    assert not real, f"{name}: {real}"

    pdf = f"{OUT}/{name}.pdf"
    fig.savefig(pdf)
    fig.savefig(f"{PNGS}/{name}.png", dpi=300)      # 2026-08-14: rasters to polished_pngs
    plt.close(fig)
    json.dump(drawn, open(f"{WORK}/{name}_polish_drawn_values.json", "w"), indent=1)
    print("wrote", pdf)
    print("wrote", f"{PNGS}/{name}.png")
    return pdf
