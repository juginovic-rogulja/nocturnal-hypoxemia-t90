#!$T90_PY
"""Round 49, lane L5: the panel a generator copy imports rc() and BAND_LABELS from New_Figures/_NOT_THESE_workfiles/scripts/splitstyle.py,
whose import now stops on an unrelated assertion (its negative-control list is pinned to a pre-v8 negcontrols_final.json). The copy takes
the two things it used from there verbatim (splitstyle.rc() over prism_plotter.manuscript_style().rc() with register_manuscript_font(),
BAND_LABELS) and drops the legend_capture monkeypatch undo (legend_capture is no longer imported, so Figure.text is the original)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
L5 = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/L5"
P4 = f"{L5}/Main_Fig5/scripts/04_build_panel_a_r49.py"; s = open(P4).read()
old = '''sys.path.insert(0, "$T90_FIGURE_ROOT/New_Figures/"
                   "_NOT_THESE_workfiles/scripts")
from splitstyle import (rc, BAND_LABELS, plt)  # noqa: E402
import legend_capture as _lc  # noqa: E402
import nooverlap  # noqa: E402
import matplotlib.text as mtext  # noqa: E402
from matplotlib.path import Path as MplPath  # noqa: E402
from matplotlib.patches import PathPatch, Rectangle  # noqa: E402

from matplotlib.axes import Axes as _Axes  # noqa: E402
from matplotlib.figure import Figure as _Figure  # noqa: E402
_Figure.text = _lc._orig_text
_Axes.text = _lc._orig_ax_text
_Axes.annotate = _lc._orig_ax_annotate
atexit.unregister(_lc.flush)
'''
new = '''sys.path.insert(0, "$T90_FIGURE_ROOT/New_Figures/"
                   "_NOT_THESE_workfiles/scripts")
# round 49: splitstyle.py no longer imports (its negative-control assertion is pinned to a pre-v8 file); rc() and BAND_LABELS are
# taken from it verbatim (splitstyle.py lines 36-43, 57, 65, 97-114) over the same prism_plotter manuscript style
sys.path.insert(0, "$T90_PRISM_PLOTTER_SRC")
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from prism_plotter.style import manuscript_style, register_manuscript_font  # noqa: E402
register_manuscript_font()
BAND_LABELS = ["0-1%", "1-5%", "5-10%", ">10%"]
_SS_TITLE, _SS_LABEL, _SS_TICK, _SS_ANNOT, _SS_INK = 11.5, 10.0, 9.5, 9.5, "#15181c"


def rc():
    r = dict(manuscript_style().rc())
    r.update({
        "figure.dpi": 300, "savefig.dpi": 300,
        "font.size": _SS_TICK, "axes.labelsize": _SS_LABEL, "axes.titlesize": _SS_TITLE,
        "xtick.labelsize": _SS_TICK, "ytick.labelsize": _SS_TICK, "legend.fontsize": _SS_ANNOT,
        "axes.linewidth": 0.9, "axes.edgecolor": _SS_INK, "text.color": _SS_INK,
        "axes.labelcolor": _SS_INK, "xtick.color": _SS_INK, "ytick.color": _SS_INK,
        "xtick.major.width": 0.9, "ytick.major.width": 0.9,
        "xtick.major.size": 3.2, "ytick.major.size": 3.2,
        "xtick.direction": "out", "ytick.direction": "out",
        "axes.spines.top": False, "axes.spines.right": False,
        "legend.frameon": False, "savefig.bbox": "tight", "savefig.pad_inches": 0.06,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })
    return r


import nooverlap  # noqa: E402
import matplotlib.text as mtext  # noqa: E402
from matplotlib.path import Path as MplPath  # noqa: E402
from matplotlib.patches import PathPatch, Rectangle  # noqa: E402
'''
assert old in s, "generator import block not found"
s = s.replace(old, new); open(P4, "w").write(s); print("patched", P4)
