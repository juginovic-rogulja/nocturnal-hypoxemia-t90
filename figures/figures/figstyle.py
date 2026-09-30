"""
Shared figure style. Every figure in this paper is drawn through here so that a reader moving
between them never has to relearn what a color or a shape means.

Built on prism-plotter's manuscript style, with the choices this paper needs added on top:

  * One color per idea, held constant everywhere. Oxygen is always the same blue. The
    apnea-hypopnea index is always the same red. Context and negative controls are always grey.
  * L-shaped axes, outward ticks, no box, no gridlines. A gridline is only added where a reader
    genuinely has to read a value off the axis rather than compare two bars.
  * Type sizes fixed at the printed size, not scaled at save time, so 7 pt is 7 pt on paper.
  * A right-hand numeric column on every forest plot. A reviewer should never have to measure
    a whisker with a ruler to find out what the interval was.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../..")))  # repository root, where paths.py lives
import paths
import sys

sys.path.insert(0, paths.PRISM_PLOTTER_SRC)  # engine 2.0, 2026-08-13

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from prism_plotter.style import manuscript_style, register_manuscript_font

register_manuscript_font()

# ---- the paper's fixed vocabulary of color ------------------------------------------------
# Muted and cool throughout. Deep petrol blue for the exposure, slate for the comparator, warm
# grey for context, deep teal and muted brick for corrected and residual. Nothing saturated,
# nothing that reads as a warning label, and everything legible in greyscale because the
# lightness values are separated as well as the hues.
OXY = "#22485f"        # oxygenation, the exposure this paper is about
OXY_L = "#7e9cad"      # a lighter tint of the same, for a secondary series
AHI = "#5c6670"        # the apnea-hypopnea index, the comparator, slate
GREY = "#8a9099"       # context, and anything deliberately uninformative
GREY_L = "#ccd1d6"
INK = "#1a1d21"        # text and axes
ACCENT = "#7d6a52"     # a third series where one is unavoidable, bronze rather than orange
GOOD = "#3f6f63"       # treated, corrected, resolved
BAD = "#8a4a4a"        # untreated, residual, failed

FAMILY = {"oximetry": OXY, "autonomic/HRV": "#4a5f7a", "EEG microstructure": "#6b5b7b",
          "EEG band power": GREY_L, "architecture/duration": ACCENT}
COHORT = {"BDSP": OXY, "SHHS": "#4a5f7a", "MrOS": ACCENT}
BANDCOL = [GOOD, "#7fa79b", GREY, ACCENT, BAD]

# A category of nocturnal hypoxemia must carry the same color in every figure, whether that
# figure splits the range into 4 categories or 5. Keying on the label rather than on position
# is what guarantees it. ">10%" is the same muted brick everywhere, in the main figures, the
# survival curves, and the treatment panels.
BAND_COLOR = {
    "0":       "#2f5c50",   # exactly zero, the deepest teal
    "<1%":     GOOD,
    "0-1%":    GOOD,
    ">0-1":    GOOD,
    ">0-1%":   GOOD,
    "1-5%":    "#7fa79b",
    ">1-5":    "#7fa79b",
    ">1-5%":   "#7fa79b",
    "5-10%":   ACCENT,
    ">5-10":   ACCENT,
    ">5-10%":  ACCENT,
    ">10%":    BAD,
    ">10":     BAD,
}


def band_color(label):
    """Color for a category of nocturnal hypoxemia, keyed by its printed label."""
    k = str(label).strip()
    if k in BAND_COLOR:
        return BAND_COLOR[k]
    k2 = k.replace(" ", "").replace("%", "")
    for key, col in BAND_COLOR.items():
        if key.replace(" ", "").replace("%", "") == k2:
            return col
    return GREY

# printed type sizes
TITLE, LABEL, TICK, ANNOT, PANEL = 9.0, 8.0, 7.5, 7.0, 10.5


def rc():
    """rcParams for every figure in the paper, sourced from the manuscript style."""
    st = manuscript_style()
    r = dict(st.rc())
    r.update({
        "figure.dpi": 300, "savefig.dpi": 300,
        "font.size": TICK, "axes.labelsize": LABEL, "axes.titlesize": TITLE,
        "xtick.labelsize": TICK, "ytick.labelsize": TICK, "legend.fontsize": ANNOT,
        "axes.linewidth": 0.9, "axes.edgecolor": INK, "text.color": INK,
        "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
        "xtick.major.width": 0.9, "ytick.major.width": 0.9,
        "xtick.major.size": 3.2, "ytick.major.size": 3.2,
        "xtick.direction": "out", "ytick.direction": "out",
        "axes.spines.top": False, "axes.spines.right": False,
        "legend.frameon": False, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })
    return r


def new(w, h):
    plt.rcParams.update(rc())
    return plt.subplots(figsize=(w, h))


def panel_label(ax, letter, dx=-0.12, dy=1.03):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=PANEL,
            fontweight="bold", va="bottom", ha="left")


def bare_y(ax):
    """Drop the y spine and ticks. Used where the y axis is a list of labels, not a scale."""
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)


def save(fig, path):
    fig.savefig(f"{path}.pdf")
    fig.savefig(f"{path}.tif", dpi=300, pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(f"{path}.png", dpi=200)
    plt.close(fig)
    return path


def forest(ax, labels, hr, lo, hi, colors=None, ref=1.0, xlabel="", numcol=True,
           logx=True, group_breaks=None, marker_size=22):
    """
    Horizontal point and interval plot.

    group_breaks is a dict of {row index: heading text}. Headings are inserted above that row
    so that related outcomes read as a block rather than as a flat list.
    """
    n = len(labels)
    y = np.arange(n)[::-1].astype(float)
    if group_breaks:
        shift = 0.0
        newy = []
        for i in range(n):
            if i in group_breaks:
                shift += 1.0
            newy.append(y[i] - shift)
        y = np.array(newy)
    colors = colors or [OXY] * n

    ax.axvline(ref, color=INK, lw=0.8, ls=(0, (4, 3)), zorder=1)
    for i in range(n):
        ax.plot([lo[i], hi[i]], [y[i], y[i]], color=colors[i], lw=1.5,
                solid_capstyle="round", zorder=2)
        ax.scatter([hr[i]], [y[i]], s=marker_size, color=colors[i], zorder=3,
                   edgecolors="white", linewidths=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    bare_y(ax)
    if logx:
        ax.set_xscale("log")
        # a log axis defaults to labelling minor ticks in scientific notation, which is
        # unreadable at this range. Label only the major ticks we set explicitly.
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(
            lambda v, _: f"{v:g}"))
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.tick_params(axis="x", which="minor", length=1.8)
    ax.set_xlabel(xlabel)
    ax.set_ylim(min(y) - 0.9, max(y) + 0.9)

    if group_breaks:
        xmin = ax.get_xlim()[0]
        for i, txt in group_breaks.items():
            ax.text(xmin, y[i] + 0.95, txt, fontsize=ANNOT, fontweight="bold",
                    va="center", ha="left", color=INK)
    if numcol:
        xr = ax.get_xlim()
        xtext = xr[1] * (1.10 if logx else 1.02)
        for i in range(n):
            ax.text(xtext, y[i], f"{hr[i]:.2f} ({lo[i]:.2f}-{hi[i]:.2f})",
                    fontsize=ANNOT, va="center", ha="left", color=INK)
        ax.set_xlim(xr)
    return y
