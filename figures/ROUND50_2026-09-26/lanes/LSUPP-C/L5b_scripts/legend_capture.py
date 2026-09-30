"""
Move the grey commentary off the figures and into the figure legends.

Alen, 2026-08-08: "lots of these figures have some grey font around them. Any sort of details
like that will go into figure legends. So clear all the figures of any sort of grey font that's
kind of everywhere."

The sheets carried two different things in grey and only one of them is commentary:

  * figure-level text, drawn with fig.text() at the page margin. Always a caption, a subtitle
    or a footnote. This is what has to move.
  * axes-level text, drawn with ax.text() inside a panel. Rank numbers, the asleep-minus-awake
    difference column, negative-control tick labels, per-bar values. This is DATA and must
    stay on the sheet.

So the rule is drawn from where the text lives, not from its wording: Figure.text is captured
and suppressed, Axes.text is left alone. Nothing is deleted. Every captured string is written
to FIGURE_LEGENDS.md in the order it appeared, under the name of the sheet it came off, so it
can be pasted into the manuscript's legend for that figure.

Import this from a style module, not from each generator, so one import covers every script.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import atexit
import os

from matplotlib.axes import Axes
from matplotlib.figure import Figure

# Every grey the house palette draws commentary in. Data greys are never routed through
# Figure.text, so this list does not need to distinguish them.
GREY_TEXT = {
    "#8a9099",  # GREY / GREY_M, the house grey
    "#6b747d",  # GREY_MID
    "#5a6169",  # the darker note grey
    "#ccd1d6",  # GREY_PALE
    "#a0a6ad",
    "#7a828b",
}

_BANK = {}          # figure name -> [captured strings], filled by note() via save()
_PENDING = []       # captured strings for the figure currently being built

ROOT = paths.FIGURE_ROOT
OUT = f"{paths.FIGURE_ROOT}/ROUND50_2026-09-26/lanes/LSUPP-C/L5b_work/FIGURE_LEGENDS_lane_capture.md"   # V14_L5b: lane work folder, never New_Figures

_orig_text = Figure.text


def _norm(c):
    """A matplotlib colour to a lowercase hex string, or None if it is not a simple colour."""
    if isinstance(c, str):
        return c.lower() if c.startswith("#") else None
    try:
        from matplotlib.colors import to_hex
        return to_hex(c).lower()
    except Exception:
        return None


def _patched_text(self, x, y, s, *args, **kwargs):
    """Capture grey figure-level text instead of drawing it."""
    if _norm(kwargs.get("color")) in GREY_TEXT:
        t = " ".join(str(s).split())
        if t:
            _PENDING.append(t)
        # Return a real but invisible artist so callers that keep the handle still work.
        kwargs = dict(kwargs, alpha=0.0)
        art = _orig_text(self, x, y, "", *args, **kwargs)
        return art
    return _orig_text(self, x, y, s, *args, **kwargs)


Figure.text = _patched_text

# ---------------------------------------------------------------------------------------
# Axes-level grey. Some sheets put their commentary inside a panel with ax.text or
# ax.annotate rather than at the page margin, so those have to move too. The danger is that
# the SAME call draws data: rank numbers, the asleep-minus-awake column, "Cell 1", "top 10",
# asterisks. Those are one or two tokens; commentary is a sentence. So prose moves and short
# labels stay.
#
# Tick labels are deliberately NOT touched. Several sheets grey out their negative-control row
# names, and "Back pain (negative control)" is four words, so a word count alone would eat it.
# Tick labels are not created by Axes.text, so patching only these two methods leaves them be.
_PROSE_MIN_WORDS = 3

_orig_ax_text = Axes.text
_orig_ax_annotate = Axes.annotate


def _is_prose(s):
    t = " ".join(str(s).split())
    return len(t.split()) >= _PROSE_MIN_WORDS and any(ch.isalpha() for ch in t)


def _patched_ax_text(self, x, y, s, *args, **kwargs):
    if _norm(kwargs.get("color")) in GREY_TEXT and _is_prose(s):
        _PENDING.append(" ".join(str(s).split()))
        return _orig_ax_text(self, x, y, "", *args, **dict(kwargs, alpha=0.0))
    return _orig_ax_text(self, x, y, s, *args, **kwargs)


def _patched_ax_annotate(self, text, *args, **kwargs):
    if _norm(kwargs.get("color")) in GREY_TEXT and _is_prose(text):
        _PENDING.append(" ".join(str(text).split()))
        return _orig_ax_annotate(self, "", *args, **dict(kwargs, alpha=0.0))
    return _orig_ax_annotate(self, text, *args, **kwargs)


Axes.text = _patched_ax_text
Axes.annotate = _patched_ax_annotate


def take(name):
    """Called from the save helpers: bank whatever was captured for this sheet and reset."""
    global _PENDING
    if _PENDING:
        _BANK[name] = list(_PENDING)
    _PENDING = []


def flush():
    """Append this run's captured legends to FIGURE_LEGENDS.md, replacing any earlier entry."""
    if not _BANK:
        return
    existing = {}
    order = []
    if os.path.exists(OUT):
        cur = None
        for line in open(OUT):
            if line.startswith("## "):
                cur = line[3:].strip()
                if cur not in existing:
                    existing[cur] = []
                    order.append(cur)
            elif cur and line.strip():
                existing[cur].append(line.rstrip())
    for k, v in _BANK.items():
        existing[k] = v
        if k not in order:
            order.append(k)
    with open(OUT, "w") as f:
        f.write("# Figure legend text, lifted off the sheets\n\n")
        f.write("Every sentence here used to be printed in grey on the figure itself. Alen asked "
                "for the sheets to carry\nonly the data, with this detail moving into the legends. "
                "Nothing was discarded: each block below is\nverbatim what came off the sheet "
                "named above it, in the order it appeared.\n\n")
        for k in sorted(order):
            f.write(f"## {k}\n\n")
            for line in existing[k]:
                f.write(line + "\n\n")
    _BANK.clear()


atexit.register(flush)
