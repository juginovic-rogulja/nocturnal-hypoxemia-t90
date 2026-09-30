"""
Style for the decompressed, split main-figure set (Figure01 ... Figure19).

Round 2 of Alen's review reversed the size brief. The earlier pass squeezed every figure to
JAMA's 7.25 in column at an 8 pt floor, and the result was dense. His instruction now is that
readability beats figure count and beats the width gate: split anything cramped into its own
numbered main figure, let the panels breathe, and he will compact later.

So this module does three things the old one did not.

  1. It sets a larger type scale. 9 pt is the floor here, not 8, and axis labels and titles sit
     at 10 to 11.5 pt. Nothing on any sheet is at the old floor.
  2. It drops the width constraint. Sheets are sized to the content, 7.5 to 11.5 in wide.
  3. It keeps the palette. Warm grey is the comparator, deep petrol blue is the measure of
     interest, and the four bands of nocturnal oxygen run one pale-to-deep ramp so an ordered
     exposure reads as ordered in greyscale.

THE ONE ADDED HUE
-----------------
Alen raised, but did not mandate, a third colour to encode a grouping. It is used in exactly
one place: BRONZE marks the Sleep Heart Health Study in the two external-validation sheets
(Figure14, Figure15), where three cohorts share one axis and two steps of grey were doing work
that a hue does better. Every other sheet is grey and blue only. Greyscale separation of the
three cohort colours is asserted at build time by sep_ok().

No individual data points. No printed P values, asterisks only. Confidence intervals are error
bars. Ns live in legends and titles, never inside the plot area.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../..")))  # repository root, where paths.py lives
import paths
import os
import sys

# 2026-08-13: engine moved to prism-plotter 2.0 (journal sizing, preflight, print
# mark weights). manuscript_style/register_manuscript_font are identical between
# versions, so the swap itself changes no pixel; the 2.0 features are opted into.
sys.path.insert(0, paths.PRISM_PLOTTER_SRC)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from prism_plotter.style import manuscript_style, register_manuscript_font

register_manuscript_font()

# ---------------------------------------------------------------- palette
BLUE = "#1f4257"        # the measure of interest, everywhere
BLUE_MID = "#4d7387"    # a second blue where one is unavoidable
BLUE_LT = "#8ba9b8"
GREY = "#8a9099"        # the comparator, everywhere
GREY_PALE = "#ccd1d6"   # deliberately uninformative, negative controls
GREY_MID = "#6b747d"
INK = "#15181c"
RULE = "#b9bec4"
BRONZE = "#7d6a52"      # the ONE added hue, external cohorts only

# the ordered exposure, pale to deep, monotone in greyscale
BAND_LABELS = ["0-1%", "1-5%", "5-10%", ">10%"]
BAND = {"0-1%": "#d5dadd", "1-5%": "#97acb8", "5-10%": "#4d7387", ">10%": BLUE}
# the top band was #1b3c4e, 2.3 CIEDE2000 from BLUE, which is a second dark blue nobody can
# see is different. It is BLUE now.
BAND_MARK = {"0-1%": "o", "1-5%": "o", "5-10%": "s", ">10%": "^"}
BAND_TEXT = {"0-1%": INK, "1-5%": INK, "5-10%": "white", ">10%": "white"}

# ---------------------------------------------------------------- type, in printed points
# 9 pt is the floor on every sheet in this set. The old set floored at 8.
TITLE, LABEL, TICK, ANNOT, SMALL, PANEL = 11.5, 10.0, 9.5, 9.5, 9.0, 13.0

MIN_LUMA_GAP = 45.0

# The negative-control panel and the confounding floor, settled 2026-08-07. Read from
# numbers/negcontrols_final.json so a figure can never print a number the analysis does not
# hold. Fracture and osteoarthritis are BOTH out, for the same reason: a plausible causal
# path runs from the exposure to each of them, which disqualifies a control however clean
# its estimate looks. Do not restate the body-mass-index argument for either; it is the
# reasoning Alen explicitly overturned.
import json as _json
import legend_capture  # grey commentary off the sheet, into the legends
_NC = _json.load(open(f"{paths.NUMBERS_DIR}/"
                      "negcontrols_final.json"))
NEG_CONTROLS = list(_NC["panel_labels"])
NEG_DROPPED = ["Fracture", "Osteoarthritis"]
FLOOR_PER_SD = _NC["confounding_floor_per_sd"]          # read from the file, never typed (v8 F1 fix)
FLOOR_SET_BY = _NC["confounding_floor_driver"]          # the control that sets it
assert NEG_CONTROLS == ["Back pain", "Cataract", "Glaucoma",
                        "Contact dermatitis", "Hemorrhoids"], NEG_CONTROLS
# v8 F1 fix (2026-09-12, integrity gate 6): the floor is checked against the file's own control table, not a typed
# number: it must be the largest per-SD control hazard ratio in negcontrols_final.json and the driver that control.
# negcontrols_final.json is a v7-era static input (built 2026-08-08 by make_negcontrols_final.py, no provenance
# sidecar, no v8 chain step rewrites it: audit check 4 STALE_NO_WRITER). The v8 primary table
# (numbers/bdsp_diseases_v3.csv, step 111) can move the control estimates; regenerate the json after step 111.
_NC_HR = {c["label"]: float(c["per_sd"]["hr"]) for c in _NC["controls"].values()}
assert len(_NC_HR) == int(_NC["n_controls"]) == len(NEG_CONTROLS), (_NC_HR, _NC["n_controls"], NEG_CONTROLS)
assert FLOOR_PER_SD == max(_NC_HR.values()) and _NC_HR[FLOOR_SET_BY] == FLOOR_PER_SD, \
    (FLOOR_PER_SD, FLOOR_SET_BY, _NC_HR)


def rc():
    r = dict(manuscript_style().rc())
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
        "legend.frameon": False, "savefig.bbox": "tight", "savefig.pad_inches": 0.06,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })
    return r


# ---------------------------------------------------------------- greyscale separation
def luma(hexcol):
    """BT.601 luma, 0 to 255. What a greyscale print of this colour collapses to."""
    h = hexcol.lstrip("#")
    r, g, b = (int(h[i:i + 2], 16) for i in (0, 2, 4))
    return 0.299 * r + 0.587 * g + 0.114 * b


def sep_ok(name, colours):
    ls = [luma(c) for c in colours]
    worst, pair = 1e9, None
    for i in range(len(ls)):
        for j in range(i + 1, len(ls)):
            d = abs(ls[i] - ls[j])
            if d < worst:
                worst, pair = d, (colours[i], colours[j])
    assert worst >= MIN_LUMA_GAP, (f"{name}: greyscale gap {worst:.1f} between {pair}, "
                                   f"needs {MIN_LUMA_GAP}")
    return worst


# ---------------------------------------------------------------- small helpers
def panel_label(ax, letter, dx=-0.10, dy=1.04):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=PANEL, fontweight="bold",
            va="bottom", ha="left", color=INK)


def bare_y(ax):
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)


def logx(ax, ticks):
    ax.set_xscale("log")
    ax.set_xticks(ticks)
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.tick_params(axis="x", which="minor", length=0)


def logy(ax, ticks):
    ax.set_yscale("log")
    ax.set_yticks(ticks)
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.tick_params(axis="y", which="minor", length=0)


def stars(p):
    """Significance as asterisks. A P value is never printed on a figure in this paper."""
    if p is None or not np.isfinite(p):
        return ""
    return "***" if p < .001 else ("**" if p < .01 else ("*" if p < .05 else ""))


def headline(fig, text, sub=None, x=0.008, y=0.995):
    """The sentence a reader should leave the sheet with, set above everything else."""
    t = fig.text(x, y, text, fontsize=TITLE, fontweight="bold", ha="left", va="top",
                 color=INK, linespacing=1.3)
    if sub:
        fig.text(x, y - 0.001, "", fontsize=1)
    return t


# ---------------------------------------------------------------- axis coverage
# Three sheets have now shipped with a value printed in the text and its mark cut off by the
# axis: Figure23 (obesity hypoventilation at 3.48 on an axis ending at 2.35), Figure01 (T90 at
# 0.020339 on an axis ending at 0.0196) and Figure13 (obesity hypoventilation at 5.13 on an axis
# ending at 4.80). Every one of them passed both delivery gates, because the gates only test a
# row that prints a parenthesised interval and only fire when the row draws nothing at all. A
# truncated mark draws something, so nothing fired.
#
# The fix belongs at build time, where the drawn values are still in hand. Call this immediately
# after the limits are set, with everything the panel draws, and the build stops rather than
# writing a sheet that lies.
AXIS_TOL_FRAC = 1e-9        # the limits are set from these values, so no slack is warranted


def axis_covers(ax, values, which="x", label="", tol_frac=AXIS_TOL_FRAC):
    """
    Assert the axis range holds every value the panel draws on it.

    values may hold None and NaN, which are skipped: a row with no upper bound draws no upper
    bound. Everything else has to sit inside the limits, on the axis's own scale, so a log axis
    is tested in logs and a value at the very edge is accepted.
    """
    lo, hi = (ax.get_xlim() if which == "x" else ax.get_ylim())
    lo, hi = (lo, hi) if lo <= hi else (hi, lo)
    scale = (ax.get_xscale() if which == "x" else ax.get_yscale())
    vals = [float(v) for v in np.ravel(np.asarray(values, dtype=object))
            if v is not None and np.isfinite(np.asarray(v, dtype=float))]
    if not vals:
        return (lo, hi)
    if scale == "log":
        assert lo > 0, f"{label or 'axis'}: log {which} axis with a non-positive limit {lo}"
        f = np.log
        vals = [v for v in vals if v > 0]
    else:
        def f(v):
            return v
    flo, fhi = f(lo), f(hi)
    tol = tol_frac * abs(fhi - flo)
    bad = [v for v in vals if f(v) < flo - tol or f(v) > fhi + tol]
    assert not bad, (f"{label or 'axis'}: {len(bad)} drawn value(s) fall outside the {which} "
                     f"axis [{lo:.6g}, {hi:.6g}], so the mark is cut off and the reader sees a "
                     f"printed number with nothing beside it: {sorted(bad)[:6]}")
    return (lo, hi)


# ---------------------------------------------------------------- output routing
NEWFIG = f"{paths.FIGURE_ROOT}/New_Figures"
MAIN = f"{NEWFIG}/1_MAIN_FIGURES"
SUPP = f"{NEWFIG}/2_SUPPLEMENTARY_FIGURES"
RAST = f"{NEWFIG}/_NOT_THESE_rasters_and_scripts"


def save(fig, name, kind="Main"):
    """
    PDF into the view folder, rasters into the working folder. No new subfolders, ever:
    Alen has reorganised this twice because agents scattered files.
    """
    d = {"Main": MAIN, "Supplementary": SUPP}[kind]
    legend_capture.take(name)
    os.makedirs(d, exist_ok=True)
    fig.savefig(f"{d}/{name}.pdf")
    fig.savefig(f"{RAST}/{name}.png", dpi=170)
    plt.close(fig)
    return f"{d}/{name}.pdf"
