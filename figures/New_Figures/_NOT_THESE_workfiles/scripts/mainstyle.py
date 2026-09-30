"""
Style for the revised main figures.

Built on the paper's existing figstyle, which is itself built on prism_plotter's manuscript
style, so type sizes, spine treatment and the save path stay identical to the figures already in
the manuscript. What this module adds is the thing Alen asked for twice: colour that separates.

The rule applied here
---------------------
Two series that a reader could confuse must differ in hue, in lightness and in shape. Lightness
is measured, not eyeballed: SEP_OK asserts that every pair of colours used inside one panel is at
least MIN_LUMA_GAP apart on the BT.601 luma scale, which is what a greyscale print or a
photocopier collapses to. Every series also carries its own marker shape and its own line dash,
so the panels survive being printed in black and white with no colour at all.

The palette is the paper's own. Deep petrol blue is the exposure. Warm grey is context and
anything deliberately uninformative. Bronze is the third series where a third is unavoidable.
Nothing new was invented, the values were only pulled apart in lightness.

The four bands of nocturnal hypoxemia are an ordered exposure, so they are drawn as one
monotone pale-to-deep blue ramp rather than as four unrelated hues. An ordinal quantity that
reads as ordinal in greyscale cannot produce the confusion Alen flagged.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../..")))  # repository root, where paths.py lives
import paths
import sys

sys.path.insert(0, f"{paths.FIGURE_ROOT}/figures")

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from figstyle import (rc, new, save, forest, panel_label, bare_y, INK, GREY, GREY_L,
                      OXY, ACCENT, GOOD, BAD, TITLE, LABEL, TICK, ANNOT, PANEL)
import legend_capture  # grey commentary off the sheet, into the legends

# ---------------------------------------------------------------- colour, with lightness held apart
BLUE_D = "#1f4257"     # the house dark blue, the exposure and the primary cohort.
                       # Was #12303f until 2026-08-07, which was a second blue for the same
                       # role as splitstyle's. One role, one colour.
BRONZE = "#7d6a52"     # the paper's existing third series
STEEL = "#8a9099"      # the house grey. Was #96a5b0, which is the same role as GREY_M and
                       # only 7 CIEDE2000 from it, so the two read as one colour anyway.
                       # Do not use STEEL and GREY_M as two series in one panel.
GREY_M = "#8a9099"     # context, negative controls
GREY_P = "#ccd1d6"     # the palest context fill

# three series that appear together in one panel, in fixed order
SERIES3 = [BLUE_D, BRONZE, STEEL]
MARK3 = ["o", "^", "s"]
DASH3 = ["-", (0, (5, 2)), (0, (1.4, 1.6))]

# the ordered exposure, pale to deep, monotone in greyscale
BANDS = ["0-1%", "1-5%", "5-10%", ">10%"]
BAND_RAMP = {"0-1%": "#d5dadd", "1-5%": "#97acb8", "5-10%": "#4d7387", ">10%": "#1f4257"}
# the same four values as splitstyle.BAND. They were within 4 CIEDE2000 of each other and
# encoded the same four bands, so they are now one ramp.

MIN_LUMA_GAP = 45.0


def luma(hexcol):
    """BT.601 luma, 0 to 255. What a greyscale print of this colour looks like."""
    h = hexcol.lstrip("#")
    r, g, b = (int(h[i:i + 2], 16) for i in (0, 2, 4))
    return 0.299 * r + 0.587 * g + 0.114 * b


def sep_report(name, colours):
    """Smallest greyscale gap between any two colours drawn in the same panel."""
    ls = [luma(c) for c in colours]
    worst, pair = 1e9, None
    for i in range(len(ls)):
        for j in range(i + 1, len(ls)):
            d = abs(ls[i] - ls[j])
            if d < worst:
                worst, pair = d, (colours[i], colours[j])
    return {"panel": name, "colours": list(colours), "luma": [round(v, 1) for v in ls],
            "min_gap": round(worst, 1), "closest_pair": list(pair), "passes": bool(worst >= MIN_LUMA_GAP)}


def sep_ok(name, colours, log=None):
    r = sep_report(name, colours)
    if log is not None:
        log.append(r)
    assert r["passes"], (f"{name}: greyscale gap {r['min_gap']} between {r['closest_pair']}, "
                         f"needs {MIN_LUMA_GAP}")
    return r


def band_col(label):
    return BAND_RAMP[str(label).strip()]


def stars(p):
    """Significance as asterisks. The paper never prints a P value on a figure."""
    if p is None or not np.isfinite(p):
        return ""
    return "***" if p < .001 else ("**" if p < .01 else ("*" if p < .05 else ""))


def ci(hr, lo, hi, dec=2):
    """A confidence interval in JAMA form, plain ASCII hyphen, never an en dash."""
    return f"{hr:.{dec}f} (95% CI, {lo:.{dec}f}-{hi:.{dec}f})"


def logticks(ax, ticks):
    ax.set_xscale("log")
    ax.set_xticks(ticks)
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.tick_params(axis="x", which="minor", length=1.8)


def forest_rows(ax, y, hr, lo, hi, colours, markers=None, ms=24, lw=1.5, hollow=None):
    """Point and interval rows, one marker shape per series so colour is never load-bearing."""
    markers = markers or ["o"] * len(y)
    hollow = hollow or [False] * len(y)
    for i in range(len(y)):
        ax.plot([lo[i], hi[i]], [y[i], y[i]], color=colours[i], lw=lw,
                solid_capstyle="round", zorder=2)
        ax.scatter([hr[i]], [y[i]], s=ms, marker=markers[i], zorder=3,
                   facecolors="white" if hollow[i] else colours[i],
                   edgecolors=colours[i] if hollow[i] else "white", linewidths=0.7)


def numcol(ax, y, hr, lo, hi, x=None, dec=2, size=None):
    """The right-hand numeric column. A reviewer should never measure a whisker with a ruler."""
    xr = ax.get_xlim()
    logx = ax.get_xscale() == "log"
    xt = x if x is not None else (xr[1] * 1.10 if logx else xr[1] + 0.02 * (xr[1] - xr[0]))
    for i in range(len(y)):
        ax.text(xt, y[i], f"{hr[i]:.{dec}f} ({lo[i]:.{dec}f}-{hi[i]:.{dec}f})",
                fontsize=size or (ANNOT - 0.5), va="center", ha="left", color=INK)
    ax.set_xlim(xr)


# ---------------------------------------------------------------- axis coverage
# Shared with splitstyle.axis_covers, which carries the full explanation. Imported rather than
# copied so the two sets of sheets cannot drift apart on what counts as clipped.
from splitstyle import axis_covers, AXIS_TOL_FRAC          # noqa: E402,F401

NEWFIG = f"{paths.FIGURE_ROOT}/New_Figures"


def _dest(kind):
    """
    Resolve where a figure of this kind belongs.

    The delivery folder has been renamed twice in one afternoon, Main -> 1_MAIN_FIGURES and
    _archive_raster -> _NOT_THESE_rasters_and_scripts, so the destination is looked up by
    pattern rather than hardcoded. A rename does not silently write a figure into a folder
    Alen is being told to ignore.
    """
    import os
    want = {"Main": ("main_figures", "main"),
            "Supplementary": ("supplementary_figures", "supplementary"),
            "raster": ("rasters_and_scripts", "archive_raster")}[kind]
    here = sorted(d for d in os.listdir(NEWFIG) if os.path.isdir(f"{NEWFIG}/{d}"))
    for pat in want:
        for d in here:
            key = d.lower().lstrip("_0123456789")
            if key.startswith(pat) or key == pat:
                return f"{NEWFIG}/{d}"
    fallback = {"Main": "1_MAIN_FIGURES", "Supplementary": "2_SUPPLEMENTARY_FIGURES",
                "raster": "_NOT_THESE_rasters_and_scripts"}[kind]
    os.makedirs(f"{NEWFIG}/{fallback}", exist_ok=True)
    return f"{NEWFIG}/{fallback}"


def save3(fig, name, kind="Main"):
    """
    Vector PDF, 300 dpi LZW TIFF, and a PNG for screen review, written into the layout Alen
    asked for on 2026-08-07: the PDF alone in the main or supplementary folder, the rasters out
    of sight in the working folder.
    """
    import os
    legend_capture.take(name)
    name = os.path.basename(name)
    pdf_dir, ras_dir = _dest(kind), _dest("raster")
    fig.savefig(f"{pdf_dir}/{name}.pdf")
    fig.savefig(f"{ras_dir}/{name}.tif", dpi=300, pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(f"{ras_dir}/{name}.png", dpi=200)
    plt.close(fig)
    return f"{pdf_dir}/{name}.pdf"
