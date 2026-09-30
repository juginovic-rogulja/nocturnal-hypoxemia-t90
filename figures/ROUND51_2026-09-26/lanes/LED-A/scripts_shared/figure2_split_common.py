"""ROUND 49 lane LED-A copy (2026-09-26): BASE, TICKF, SMALLF + 1 pt. Original: FINAL_FIGURES_2026-08-14/_scripts/figure2_split_common.py.

Shared plumbing for the five split, decluttered Figure 2 panels (figure2A_split.py ...
figure2E_split.py).

Alen's verdict on the composed Figure 2 was that it is insanely busy, so each panel becomes
its own standalone sheet with the text stripped back to what a reader needs to decode the
marks: no headline, no panel letters, no on-sheet commentary, no per-row printed estimate
columns. Structural labels stay (organ groups, band labels, condition names, axis labels,
minimal series legends).

House style is the frozen splitstyle module, applied through plt.rc_context only.
splitstyle.save() is never called: it writes into the frozen New_Figures tree. save_split()
below is the only door, and it writes into FINAL_FIGURES_2026-08-14/Main_Figures_Split only.

Two splitstyle side effects are undone at import, exactly as t90final.py undoes them:
legend_capture's monkey-patched text methods are restored (grey structural labels must draw
what they are told to draw) and its atexit flush into the frozen New_Figures/FIGURE_LEGENDS.md
is unregistered.

Every drawn number is read from the same frozen result files as the composite
(make_final_fig2_t90_and_sleep_duration.py) and asserted against the same pinned values.
A build that disagrees with its source aborts. The composite and its outputs are not touched.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import atexit
import json
import os
import sys

SCRIPTS = (f"{paths.FIGURE_ROOT}/New_Figures/"
           "_NOT_THESE_workfiles/scripts")
sys.path.insert(0, SCRIPTS)

from splitstyle import (rc, bare_y, logx, logy, axis_covers, sep_ok,            # noqa: E402,F401
                        BAND, BAND_LABELS, BAND_MARK, BLUE, BLUE_MID, GREY, GREY_PALE,
                        GREY_MID, INK, RULE, plt, NEG_CONTROLS, FLOOR_PER_SD, FLOOR_SET_BY)
import legend_capture as _lc                                                    # noqa: E402
import nooverlap                                                                # noqa: E402
from prism_plotter.journal import preflight_figure, journal_style               # noqa: E402
from matplotlib.ticker import NullLocator                                       # noqa: E402

# restore the patched text methods and disarm the frozen-folder flush (t90final.py pattern)
from matplotlib.axes import Axes as _Axes                                       # noqa: E402
from matplotlib.figure import Figure as _Figure                                 # noqa: E402
_Figure.text = _lc._orig_text
_Axes.text = _lc._orig_ax_text
_Axes.annotate = _lc._orig_ax_annotate
atexit.unregister(_lc.flush)

ROOT = paths.FIGURE_ROOT
NUM = f"{paths.NUMBERS_DIR}"
SV = paths.SV_ROOT
OUT = f"{paths.FIGURE_ROOT}/FINAL_FIGURES_2026-08-14/Main_Figures_Split"
WORK = f"{paths.FIGURE_ROOT}/FINAL_FIGURES_2026-08-14/_workfiles"
os.makedirs(OUT, exist_ok=True)
os.makedirs(WORK, exist_ok=True)

# 2026-08-14 audit fix carried over from the composite: not-significant markers are the house
# pale #ccd1d6, which holds 25 greyscale tones against the control grey.
PALE = GREY_PALE

W_IN = 171.0 / 25.4          # every split sheet is JAMA double column wide, heights vary
JSTYLE = journal_style("jama", columns=2)

# standalone type scale: base 11 pt, ticks and names 10 pt, small annotations 9.5, floor 9
BASE, TICKF, SMALLF = 12.0, 11.0, 10.5      # round 49 (lane LED-A copy): 11.0, 10.0, 9.5 + 1 pt on every text element
FONT_FLOOR = 9.0


def rc_split():
    r = dict(rc())
    r.update({
        "savefig.bbox": "standard",       # fixed canvas, 171 mm stays 171 mm
        "font.size": BASE, "axes.labelsize": BASE, "axes.titlesize": BASE,
        "xtick.labelsize": TICKF, "ytick.labelsize": TICKF, "legend.fontsize": TICKF,
    })
    return r


def die(cond, msg):
    assert cond, msg


def logx_clean(ax, ticks):
    logx(ax, ticks)
    ax.xaxis.set_minor_locator(NullLocator())


def logy_clean(ax, ticks):
    logy(ax, ticks)
    ax.yaxis.set_minor_locator(NullLocator())


def _all_strings(fig):
    out = []
    for t in fig.texts:
        out.append(t)
    for ax in fig.axes:
        out += list(ax.texts)
        out += ax.get_xticklabels() + ax.get_yticklabels()
        out += [ax.xaxis.label, ax.yaxis.label, ax.title]
        lg = ax.get_legend()
        if lg is not None:
            out += list(lg.get_texts())
    for lg in getattr(fig, "legends", []):
        out += list(lg.get_texts())
    return [t for t in out if t.get_text().strip()]


def _prose_gate(fig, name):
    """No em-dashes, no semicolons, nothing under the 9 pt floor, on any sheet."""
    for t in _all_strings(fig):
        s = t.get_text()
        assert "—" not in s, f"{name}: em-dash on sheet in {s!r}"
        assert ";" not in s, f"{name}: semicolon on sheet in {s!r}"
        assert t.get_fontsize() >= FONT_FLOOR, (
            f"{name}: {t.get_fontsize()} pt is under the {FONT_FLOOR} pt floor in {s!r}")


def save_split(fig, name, drawn):
    """Prose gate, canvas-cut gate, overlap gate, journal preflight, then PDF + 300 dpi PNG
    into Main_Figures_Split, plus the drawn-values provenance JSON into _workfiles."""
    _prose_gate(fig, name)

    # the canvas is fixed (savefig.bbox standard), so text past an edge is CUT, and the
    # overlap gate does not look for that. This does, with 2 px of slack.
    boxes, _ = nooverlap._boxes(fig)
    Wpx, Hpx = fig.get_size_inches() * fig.dpi
    _cut = [f"{b['where']}: {b['text'][:44]}" for b in boxes
            if b["bbox"].x0 < -2 or b["bbox"].y0 < -2
            or b["bbox"].x1 > Wpx + 2 or b["bbox"].y1 > Hpx + 2]
    assert not _cut, f"{name}: text cut by the canvas edge: {_cut}"

    nooverlap.gate(fig, name)
    problems = preflight_figure(fig, JSTYLE)
    real = [p for p in problems if "height" not in p.lower()]
    assert not real, f"{name}: {real}"

    pdf = f"{OUT}/{name}.pdf"
    fig.savefig(pdf)
    fig.savefig(f"{OUT}/{name}.png", dpi=300)
    plt.close(fig)
    json.dump(drawn, open(f"{WORK}/{name}_drawn_values.json", "w"), indent=1)
    print("wrote", pdf)
    print("wrote", f"{OUT}/{name}.png")
    return pdf
