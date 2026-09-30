#!$T90_PY
"""Round 54 (2026-09-29), lane L2, Main_Fig2: shared paths, style and helpers for the three panel builders, the composer and the
verifier. Copy of ROUND49_2026-09-26/lanes/L2/Main_Fig2/scripts/fig2_common.py (byte copy beside as fig2_common_R49_ORIGINAL.py)
with the round-54 changes only:
  PLUS = 4.0: pt(size) adds four points to the round-37 base sizes (10 -> 14 ticks, row labels, keys, values; 11 -> 15 axis titles,
  column headers, panel titles; 9 -> 13 small notes and the asterisk key; 14 -> 18 letters and the title strip, stamped by the composer);
  SIZE names the classes so that every builder passes a base size by class, never a literal;
  V30 (ROUND52 NEW_FINAL_SET_V30/Main_Fig2.pdf, byte-identical to the round-49 L2 output) is the sheet the census runs against;
  the geometry is no longer pinned to the V13 probe: work/layout.json (written by 01_build_a_r54.py) carries the round-54 page height,
  the a/b split and the b/c split. The V13 probe caches (work/r37_base_geometry.json, r37_base_text.json) are still read for the
  block membership of panel a and for the style constants the round-37 builders measured; the V13 sheet itself is never opened."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import hashlib, json, math, os, re, subprocess, sys
import matplotlib
matplotlib.use("Agg")
from matplotlib import font_manager as fm
for _f in ("Arial.ttf", "Arial Bold.ttf", "Arial Italic.ttf"):
    fm.fontManager.addfont(f"{paths.FONT_DIR}/{_f}")
import matplotlib.pyplot as plt
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND37_2026-09-14/figures/_common")
import v14lib
from v14lib import hydrated, sha256, blocks, sidecar_v8, hr_ci, q_text, halfup, int_comma

PLUS = 4.0                                  # round 54: every base size four points larger (round 49: 1.0)
def pt(size): return round(float(size) + PLUS, 3)
# base sizes by class (round-37 values): the builders pass these, pt() adds PLUS
SIZE = dict(tick=10.0, label=10.0, key=10.0, value=10.0, head=11.0, title=11.0, note=9.0, letter=14.0)
TARGET = {k: pt(v) for k, v in SIZE.items()}   # tick/label/key/value 14, head/title 15, note 13, letter 18

LANE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORK = os.path.join(LANE, "work"); VER = os.path.join(LANE, "verify"); CROPS = os.path.join(VER, "crops")
for _d in (WORK, VER, CROPS): os.makedirs(_d, exist_ok=True)
T90 = paths.FIGURE_ROOT
NUM = f"{paths.NUMBERS_DIR}"
OLD = f"{paths.V8_ROOT}/compare/old_numbers/numbers"       # v7 snapshot (record only, not run in round 54)
SHEET = f"{paths.FIGURE_ROOT}/ROUND36_2026-09-10/L9_assembly/NEW_FINAL_SET_V13/Main_Fig2.pdf"   # the V13 design baseline (record only, never opened here)
V30 = f"{paths.FIGURE_ROOT}/ROUND52_2026-09-26/figures/NEW_FINAL_SET_V30/Main_Fig2.pdf"          # the sheet reproduced at the round-54 sizes
LOCK = f"{paths.FIGURE_ROOT}/ROUND30_2026-09-08/_render.lock"
PAGE_W = 968.66                            # unchanged
PAGE_H_V30 = 1184.76
INK, GRID, BLUE, PALE, GREY = "#1a1d21", "#eef0f1", "#0288d1", "#ccd1d6", "#8a9099"
LADDER = ["#b3dcf2", "#7cc0e9", "#3f9fd8", "#0288d1"]
RC = {"font.family": "Arial", "font.size": pt(10.0), "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
      "xtick.major.size": 3.0, "ytick.major.size": 3.0, "xtick.direction": "out", "ytick.direction": "out",
      "axes.edgecolor": INK, "text.color": INK, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
      "axes.spines.top": False, "axes.spines.right": False, "axes.facecolor": "#ffffff",
      "figure.facecolor": "none", "savefig.facecolor": "none", "savefig.bbox": "standard", "savefig.pad_inches": 0.0,
      "pdf.fonttype": 42, "figure.dpi": 300, "savefig.dpi": 300, "mathtext.default": "regular", "axes.unicode_minus": False,
      "lines.solid_capstyle": "round"}
# Arial vertical metrics used for glyph boxes and vertical centring: cap height 0.716, descender 0.212 (of the size)
CAP, DESC = 0.716, 0.212
def centre_dy(size): return 0.5 * CAP * size      # baseline sits this far below the vertical centre of a capital


def sidecar(path):
    """v8.1 provenance gate (v14lib): the sidecar must name data_frozen_v8_2026-09 or a v8.1 STATUS ok row."""
    return sidecar_v8(path)


def layout():
    """work/layout.json written by 01_build_a_r54.py: page height, the a/b split, the b/c split, the panel boxes."""
    return json.load(open(hydrated(os.path.join(WORK, "layout.json"))))


def page_h():
    return float(layout()["page_h"])


def load_geometry():
    return json.load(open(hydrated(os.path.join(WORK, "r37_base_geometry.json"))))


def load_text():
    return json.load(open(hydrated(os.path.join(WORK, "r37_base_text.json"))))


def spans_in(text, region):
    """Spans of the V13 sheet's text layer by region: a (x<485, y<739), b (x>=485, y<739), c (y>=739)."""
    out = []
    for s in text["spans"]:
        x, y = s["origin"]
        r = "c" if y >= 739 else ("a" if x < 485 else "b")
        if r == region: out.append(s)
    return out


def hexcol(c):
    return None if c is None else "#" + "".join(f"{int(round(v * 255)):02x}" for v in c[:3])


def page_figure(H=None):
    """A transparent page-sized figure with an overlay axes whose data units are page points, y downward."""
    H = page_h() if H is None else H
    fig = plt.figure(figsize=(PAGE_W / 72.0, H / 72.0)); fig.patch.set_alpha(0.0)
    ov = fig.add_axes([0, 0, 1, 1]); ov.set_xlim(0, PAGE_W); ov.set_ylim(H, 0); ov.axis("off"); ov.patch.set_visible(False)
    return fig, ov


_WFIG = None
def text_width_pt(fig, s, size, weight="normal", style="normal"):
    global _WFIG
    if fig is None:
        if _WFIG is None:
            with plt.rc_context(RC): _WFIG = plt.figure(figsize=(2, 2))
        fig = _WFIG
    r = fig.canvas.get_renderer()
    t = fig.text(0.5, 0.5, s, fontsize=size, fontweight=weight, fontstyle=style, fontfamily="Arial")
    w = t.get_window_extent(renderer=r).width / fig.dpi * 72.0
    t.remove(); return w


def bh(pvals):
    """Benjamini-Hochberg adjusted P values (statsmodels fdr_bh, the paper's convention)."""
    from statsmodels.stats.multitest import multipletests
    return multipletests(list(pvals), method="fdr_bh")[1]


def hr_ci_text(hr, lo, hi):
    """HALF UP on the value given (v14lib.hr_ci, Decimal ROUND_HALF_UP)."""
    return hr_ci(hr, lo, hi)


def exact_half_2dp(v):
    """True when a stored value sits exactly on a half at two decimals (the only case a full-precision refit could move the digit)."""
    s = repr(float(v))
    return bool(re.fullmatch(r"-?\d+\.\d\d5", s))


def fit_log_axis(tick_x, tick_v):
    """x = a + b ln(v) through the tick marks; returns a, b and the largest residual."""
    import numpy as np
    lv = np.log(np.asarray(tick_v, float)); x = np.asarray(tick_x, float)
    b, a = np.polyfit(lv, x, 1); res = float(np.abs(a + b * lv - x).max())
    return float(a), float(b), res


def drawn_rec(text, x, baseline, ha="left", size=10.0, weight="normal", source_file="static:V13", source_key="", source_value="", rule="text", panel="", **extra):
    """One v14lib DRAWN_RECORD."""
    d = dict(text=str(text), x=round(float(x), 3), baseline=round(float(baseline), 3), ha=ha, size=size, weight=weight,
             source_file=source_file, source_key=source_key, source_value=source_value, rule=rule, panel=panel)
    d.update(extra); return d
