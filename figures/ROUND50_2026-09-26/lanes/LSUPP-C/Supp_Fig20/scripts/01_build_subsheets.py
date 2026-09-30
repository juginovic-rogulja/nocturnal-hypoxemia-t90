#!/usr/bin/env python3
"""Lane V7_L5b_PAP_SUPP, Supp 20 step 01: rebuild the eFigure 12B, 12C and 12D sub-sheets from a drawn-values JSON.

This is FINAL_FIGURES_2026-08-14/_scripts/redesign_eFigure12_supp_polish.py (2026-08-26 10:26, the builder of the
polished sub-sheets that the round-30 Supp_Fig20 base is composed of), with these lane changes and nothing else:
  1. SRC is a command-line choice: --old reads the 2026-08-07 eFigure12_drawn_values.json (the base sheet's
     numbers, used only to prove that this builder reproduces the base), --new reads this lane's v7 JSON
     (work/eFigure12_drawn_values_v7.json, stamped). Outputs go to work/subsheets_<mode>/, never to FINAL_FIGURES.
  2. Colours sampled from the base sheet's live render (RECOMPOSE_R30_RULES colour law, the LA_fig2 rule): the
     primary blue is #0288d1 on the base (the polish palette's #1f4257 is not on the sheet). Grey #8a9099, pale grey
     #ccd1d6, ink #1a1d21, bands #eef0f1, zone #ebedef are unchanged.
  3. No family letters: the base sheet keeps its letters (a, b, c, d) in its own strips, and the compose step drops
     only the data regions of the new sub-sheets into the base (flat cut-out). family_letter() is a no-op.
  4. 12C carries the base's later design patches: both column headers right-aligned at one x (the base's position,
     read from its text layer), and the two-entry legend centred under the axes (lower center) instead of lower left.
     The 12B legend and headers are unchanged in position (the base's 12B legend strip is kept from the base).
  5. Every literal result pin of the polish script (the 7 and 0 counts, the starred list, FLOORS_EXPECTED, the
     0.004 to 0.135 fall range) is replaced by a check against the files: counts must agree with the JSON's own
     B_counts, the recomputed floors must agree with attack_summary.json (both v7 outputs of step 46) within 5e-4,
     and the starred list, the falls and every printed value are RECORDED to work/subsheets_<mode>/build_record.json
     for the proof and the report, never asserted to a typed number. Row counts are recorded and compared with the
     base sheet's row counts by the compose step, which decides the cut-out geometry.
  6. 12A is not rebuilt: its twelve rows are identical to the base to full float precision (drawn_values_old_vs_new),
     so the base's panel a is kept and read back against the v7 file by the verify step.
The design layer (type sizes, pitches, gutters, markers, gates, save gates) is verbatim.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import atexit, glob, json, os, re, sys
sys.dont_write_bytecode = True
import numpy as np, pandas as pd

SCRIPTS = f"{paths.FIGURE_ROOT}/ROUND50_2026-09-26/lanes/LSUPP-C/L5b_scripts"   # R40: the lane copies of the design modules   # the lane's repointed design modules (splitstyle, legend_capture, nooverlap)
sys.path.insert(0, SCRIPTS)
from splitstyle import rc, save, panel_label, GREY, GREY_PALE, plt   # noqa: E402  (BLUE replaced below, rule 2)
from splitstyle import axis_covers, stars                             # noqa: E402,F401
assert callable(save)
import legend_capture as _lc                                          # noqa: E402
import nooverlap                                                      # noqa: E402
from matplotlib import ticker as mticker                              # noqa: E402
from matplotlib.axes import Axes as _Axes                             # noqa: E402
from matplotlib.figure import Figure as _Figure                       # noqa: E402
from matplotlib.lines import Line2D                                   # noqa: E402
from matplotlib.transforms import blended_transform_factory           # noqa: E402

_Figure.text = _lc._orig_text
_Axes.text = _lc._orig_ax_text
_Axes.annotate = _lc._orig_ax_annotate
atexit.unregister(_lc.flush)

LANE20 = f"{paths.FIGURE_ROOT}/ROUND50_2026-09-26/lanes/LSUPP-C/Supp_Fig20"   # R40: this lane's sheet folder (work/ records copied from the round-37 lane)
MODE = "new" if "--new" in sys.argv else ("old" if "--old" in sys.argv else None)
assert MODE, "usage: 01_build_subsheets.py --old | --new"
SRC = (f"{LANE20}/work/eFigure12_drawn_values_v8_2.json" if MODE == "new" else
       f"{LANE20}/work/eFigure12_drawn_values_v7_V13.json")   # --old = the values the V13 sheet prints (positive control of this builder)
OUT = f"{LANE20}/work/subsheets_{MODE}"; PNGS = OUT; WORK = OUT
os.makedirs(OUT, exist_ok=True)
st = os.stat(SRC); assert st.st_size > 0 and st.st_blocks > 0, ("evicted", SRC)

BLUE = "#0288d1"          # rule 2: the blue that is on the base sheet (81 fills, 60 CI strokes), not the polish #1f4257
assert GREY.lower() == "#8a9099" and GREY_PALE.lower() == "#ccd1d6"

# ---------------------------------------------------- the design system, portrait copy (verbatim)
W = 183.0 / 25.4
H_MAX = 247.0 / 25.4
INK = "#1a1d21"
PT_PLUS = 1.0   # R49: every text one point larger, geometry unchanged
LABF, TCKF, ANNF, FLOOR = 11.0 + PT_PLUS, 10.0 + PT_PLUS, 10.0 + PT_PLUS, 9.0
PTLAB = 9.5   # R49 exception: the 16 scatter-point labels of panel d keep 9.5 pt (at 10.5 pt the no-overlap solver finds no clear place for 'Death from any cause' among its 20 candidate offsets; the alternative, farther candidates with leader lines, would re-lay the whole cloud)
STARF = 10.0 + PT_PLUS
ROWLAB = TCKF - PT_PLUS   # R49 exception: the row labels of panels a, b and c keep 10 pt (at 11 pt 'Peripheral artery disease' would cross panel a's edge by 2.8 pt onto a mark 4.3 pt from the edge, 'Cardiovascular composite' would sit 3.4 pt on the band strip of b and c)
CI_LW, MEDGE = 1.7, 0.8
S_CIRCLE, S_SQUARE, MS_CIRCLE, MS_SQUARE = 31.0, 25.5, 6.3, 5.7
EDGE = 0.18
GUT, AX_X0, ROW_IN = 0.18, 1.89, 0.34
RIGHT = W - EDGE

# rule 4: the base sheet's 12C header anchors, read from its text layer (Supp_Fig20/work/base_text.json):
# 'Healthy subgroup, HR (95% CI)' bbox x1 493.103, baseline 750.809; 'P interaction' bbox x1 492.598, baseline 765.935;
# the 12C sub-sheet is placed at (14.0, 737.26). Both right edges sit at one x within 0.5 pt, so one anchor.
_BT = json.load(open(f"{LANE20}/work/base_text.json"))["spans"]
_h1 = [s for s in _BT if s["text"] == "Healthy subgroup, HR (95% CI)"]; _h2 = [s for s in _BT if s["text"] == "P interaction"]
assert len(_h1) == 1 and len(_h2) == 1, (len(_h1), len(_h2))
_C_OFF = (14.0, 737.26)
C_HDR_RIGHT_IN = (max(_h1[0]["bbox"][2], _h2[0]["bbox"][2]) - _C_OFF[0]) / 72.0
C_HDR_Y_IN = ((_h1[0]["origin"][1] - _C_OFF[1]) / 72.0, (_h2[0]["origin"][1] - _C_OFF[1]) / 72.0)   # baselines from the top


def rc_polish():
    r = dict(rc())
    r.update({
        "savefig.bbox": "standard",
        "font.size": TCKF, "axes.labelsize": LABF, "axes.titlesize": LABF,
        "xtick.labelsize": TCKF, "ytick.labelsize": TCKF, "legend.fontsize": TCKF,
        "axes.linewidth": 0.8, "axes.edgecolor": INK, "text.color": INK,
        "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
        "xtick.major.width": 0.8, "ytick.major.width": 0.8,
        "xtick.major.size": 3.0, "ytick.major.size": 3.0,
    })
    return r


def frac(fig, x_in, y_in, w_in, h_in):
    Wc, Hc = fig.get_size_inches()
    return [x_in / Wc, y_in / Hc, w_in / Wc, h_in / Hc]


def logx(ax, ticks):
    ax.set_xscale("log")
    ax.set_xticks(ticks)
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.xaxis.set_minor_locator(mticker.NullLocator())


def bare_y(ax):
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)


def bands(ax, n):
    for i in range(n):
        if i % 2 == 1:
            yy = n - 1 - i
            ax.axhspan(yy - 0.5, yy + 0.5, color="#eef0f1", lw=0, zorder=0)


def refline(ax, x=1.0):
    ax.axvline(x, color=INK, lw=0.9, ls=(0, (4, 3)), zorder=1)


def left_rows(ax, ys, labels, ax_x0=AX_X0, size=TCKF):
    ax.set_yticks(list(ys))
    ax.set_yticklabels(labels, fontsize=size, color=INK, ha="left")
    ax.tick_params(axis="y", pad=(ax_x0 - GUT) * 72.0, length=0)


def measure(fig, s, size, weight="normal"):
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    t = fig.text(0.5, 0.5, s, fontsize=size, fontweight=weight)
    w = t.get_window_extent(renderer=r).width / fig.dpi
    t.remove()
    return w


def fit_gutter(fig, texts, fontsize=TCKF, ax_x0=AX_X0, safety=0.10):
    for s in texts:
        w = measure(fig, s, fontsize)
        assert GUT + w <= ax_x0 - safety, (f"gutter label {s!r} is {w:.2f} in wide and would cross the plot edge at {ax_x0:.2f} in")


def _pts_up(ax, pts):
    bb = ax.get_position()
    return pts / 72.0 / (bb.height * ax.figure.get_size_inches()[1])


def family_letter(fig, lt, y_in=0.17):
    """Rule 3: no letters on the lane sub-sheets, the base sheet keeps its own."""
    RECORD.setdefault("letters_suppressed", []).append(lt)


def plain_ci(ax, lo, hi, hr, y, col, marker="o", size=S_CIRCLE):
    ax.plot([lo, hi], [y, y], color=col, lw=CI_LW, solid_capstyle="round", zorder=3)
    ax.scatter([hr], [y], s=size, marker=marker, color=col, zorder=5, edgecolors="white", linewidths=MEDGE)


def fmt_p(p):
    return "<0.001" if p < 0.001 else f"{p:.3f}"


_BOLD_OK = re.compile(r"^([a-g]|<0\.001|\d\.\d{3})$")


def _all_strings(fig):
    out = list(fig.texts)
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


def save_polished(fig, name, gatelog=None):
    """Prose, floor, edge, overlap and portrait-limit gates, then PDF + 300 dpi PNG (verbatim)."""
    for t in _all_strings(fig):
        s = t.get_text()
        for bad, nm in (("—", "em dash"), (";", "semicolon"), ("×", "multiplication sign")):
            assert bad not in s, f"{name}: banned {nm} on sheet in {s!r}"
        if str(t.get_fontweight()) in ("bold", "700"):
            assert _BOLD_OK.match(s), f"{name}: bold outside letters/values: {s!r}"
        assert t.get_fontsize() >= FLOOR, (f"{name}: {t.get_fontsize()} pt under the {FLOOR} pt floor in {s!r}")
    boxes, _ = nooverlap._boxes(fig)
    Wpx, Hpx = fig.get_size_inches() * fig.dpi
    m = (4.0 - 0.3) / 25.4 * fig.dpi
    # Rule 4b (relaunch 5): the two 12C column headers sit 13.5 and 28.7 pt below the sub-sheet's top edge,
    # inside this sheet's own 4 mm design margin. That is where the ROUND-30 BASE SHEET puts them, as page-level
    # stamps outside the sub-sheet, and on the composed page they are 750.8 pt from the page top, nowhere near a
    # margin. They are drawn here so they carry real Arial with real spaces (a TextWriter stamp writes the space
    # glyph back as a non-breaking space and breaks the word multiset), and they are exempt from this gate only.
    EDGE_EXEMPT = {"Healthy subgroup, HR (95% CI)", "P interaction", "REM window only\nHR (95% CI)"}   # R49: the two-line REM column header of 12A at 11 pt ends 3.2 mm from this sub-sheet's right edge (on the composed page 46 pt from panel b's labels), nothing moved
    tight = [f"{b['where']}: {b['text'][:44]}" for b in boxes
             if b["text"] not in EDGE_EXEMPT
             and (b["bbox"].x0 < m or b["bbox"].y0 < m or b["bbox"].x1 > Wpx - m or b["bbox"].y1 > Hpx - m)]
    assert not tight, f"{name}: text inside the 4 mm edge margin: {tight}"
    nooverlap.gate(fig, name, gatelog)
    w_mm, h_mm = fig.get_size_inches() * 25.4
    assert abs(w_mm - 183.0) < 0.5, f"{name}: canvas is {w_mm:.1f} mm, not 183 mm"
    assert h_mm <= 247.0 + 0.5, f"{name}: {h_mm:.0f} mm tall, over the 247 mm portrait limit"
    pdf = f"{OUT}/{name}.pdf"
    fig.savefig(pdf)
    fig.savefig(f"{PNGS}/{name}.png", dpi=300)
    plt.close(fig)
    print("wrote", pdf, f"({w_mm:.0f} x {h_mm:.1f} mm)")
    RECORD.setdefault("sheets", {})[name] = {"pdf": pdf, "w_mm": round(w_mm, 2), "h_mm": round(h_mm, 2)}
    return pdf


# ---------------------------------------------------------------- printed-number audit (verbatim)
SHOWN = []
RECORD = {"mode": MODE, "src": SRC, "blue": BLUE}


def shown(text, value, source, what):
    v, s = float(value), float(source)
    assert np.isfinite(v) and abs(v - s) < 1e-12, f"{what}: shown {v!r} != source {s!r}"
    t = str(text)
    assert ";" not in t and "—" not in t, f"{what}: banned punctuation in {t!r}"
    SHOWN.append({"what": what, "text": t, "value": v})
    return t


D = json.load(open(SRC))
GATE, WROTE = [], []
N3_KEY = "N3, too little N3 sleep in most patients to interpret"


def n3_key_legend(fig, y_in, h_in, x_in=GUT):
    fig.legend(handles=[Line2D([], [], ls="none", marker="x", ms=6.3, mew=1.4, color=GREY_PALE, label=N3_KEY)],
               loc="lower left", bbox_to_anchor=(x_in / W, y_in / h_in), fontsize=ANNF, frameon=False, handletextpad=0.5)


# =============================================== 12A, the windows track each other (ported from redesign_eFigure12_supp_polish.py lines 337 to 421, v8.2: panel a is rebuilt because its values change)
WINDOWS_A = [("wake", "Wake"), ("sleep", "Whole sleep"), ("nrem", "NREM"), ("n1", "N1"),
             ("n2", "N2"), ("n3", "N3"), ("rem", "REM (95% CI)")]


def panel_head(ax, letter, title=None, size=ANNF, x=0.0, lines_for_lift=None):
    """polish round 4: letter=None draws the window title alone above the panel."""
    if title:
        ax.annotate(title, xy=(x, 1.0), xycoords="axes fraction", xytext=(0, 6), textcoords="offset points",
                    ha="left", va="bottom", fontsize=size, color=INK, linespacing=1.18)


def fig_a():
    a = pd.DataFrame(D["A"]).iloc[::-1].reset_index(drop=True)      # the JSON holds bottom-to-top, top row first here
    n = len(a)
    assert n == 12, n
    GAPX = 0.14
    AX_X0_A = 1.84
    ROW_IN_A = 0.3055
    PW = (RIGHT - AX_X0_A - 3 * GAPX) / 4.0
    AX_H = n * ROW_IN_A
    BOT, INTER, TOPBAND = 0.86, 0.80, 0.72
    H = BOT + AX_H + INTER + AX_H + TOPBAND
    RECORD["A"] = {"n": n, "rows_top_to_bottom": list(a.outcome),
                   "printed": [{"outcome": r.outcome, "rem_hr": f"{r.t90_rem_hr:.2f} ({r.t90_rem_lo:.2f}-{r.t90_rem_hi:.2f})",
                                "src": {w: float(getattr(r, f"t90_{w}_hr")) for w, _l in WINDOWS_A} | {"t90_rem_lo": float(r.t90_rem_lo), "t90_rem_hi": float(r.t90_rem_hi)}} for r in a.itertuples()]}
    with plt.rc_context(rc_polish()):
        fig = plt.figure(figsize=(W, H))
        fit_gutter(fig, list(a.outcome), fontsize=ROWLAB, ax_x0=AX_X0_A)   # R49 exception (ROWLAB)
        y = np.arange(n, dtype=float)[::-1]
        y0_by_bank = {0: BOT + AX_H + INTER, 1: BOT}
        rem_ax = None
        for j, (w, lab) in enumerate(WINDOWS_A):
            bank, col_j = (0, j) if j < 4 else (1, j - 4)
            ax = fig.add_axes(frac(fig, AX_X0_A + col_j * (PW + GAPX), y0_by_bank[bank], PW, AX_H))
            bands(ax, n)
            refline(ax)
            col = BLUE if w == "rem" else (GREY_PALE if w == "n3" else GREY)
            vals = a[f"t90_{w}_hr"].astype(float).values
            drawn = list(vals)
            if w == "rem":
                rem_ax = ax
                for yy, lo, hi in zip(y, a.t90_rem_lo, a.t90_rem_hi):
                    ax.plot([lo, hi], [yy, yy], color=col, lw=CI_LW, solid_capstyle="round", zorder=3)
                drawn += list(a.t90_rem_lo) + list(a.t90_rem_hi)
            if w == "n3":
                ax.scatter(vals, y, s=33, marker="x", color=col, linewidths=1.4, zorder=4)
            else:
                ax.scatter(vals, y, s=S_CIRCLE, marker="o", color=col, zorder=4, edgecolors="white", linewidths=MEDGE)
            if col_j == 0:
                left_rows(ax, y, list(a.outcome), ax_x0=AX_X0_A, size=ROWLAB)   # R49 exception
            else:
                ax.set_yticks(y)
                ax.set_yticklabels([])
                ax.tick_params(axis="y", length=0)
            bare_y(ax)
            logx(ax, [1.0, 1.5, 2.0])
            ax.set_xlim(0.88, 2.30)
            ax.set_ylim(-0.7, n - 0.3)
            axis_covers(ax, drawn, "x", f"eFigure12A {lab} window")
            panel_head(ax, None, lab, size=ANNF, lines_for_lift=1)
        colx_in = AX_X0_A + 3 * (PW + GAPX) + 0.04
        tr = blended_transform_factory(fig.transFigure, rem_ax.transData)
        for yy, hr_, lo, hi in zip(y, a.t90_rem_hr, a.t90_rem_lo, a.t90_rem_hi):
            fig.text(colx_in / W, yy, f"{hr_:.2f} ({lo:.2f}-{hi:.2f})", transform=tr, fontsize=ANNF, va="center", ha="left", color=INK)
        fig.text(colx_in / W, n - 0.05, "REM window only\nHR (95% CI)", transform=tr, fontsize=ANNF, va="bottom", ha="left", color=INK, linespacing=1.25)
        fig.text((AX_X0_A + RIGHT) / 2.0 / W, (BOT - 0.34) / H, "Hazard ratio per 1 SD of T90 in the window", ha="center", va="center", fontsize=LABF, color=INK)
        n3_key_legend(fig, 0.05, H)
        family_letter(fig, "a")
        WROTE.append(save_polished(fig, "eFigure12A_windows_track_each_other", GATE))
    shown("12", 12, len(a), "12A n conditions")


# =============================================== 12B, REM against whole sleep, paired (verbatim except rules 3, 5)
def fig_b():
    b = pd.DataFrame(D["B"]).sort_values("rem_hr", ascending=False).reset_index(drop=True)
    n = len(b)
    assert n == D["B_counts"]["n_outcomes"], (n, D["B_counts"])
    n_rem_stars = int((b.lr_p < .05).sum())
    n_slp_keep = int((b.slp_p < .05).sum())
    assert n_rem_stars == D["B_counts"]["rem_adds_p05"], (n_rem_stars, D["B_counts"])
    assert n_slp_keep == D["B_counts"]["sleep_keeps_p05"], (n_slp_keep, D["B_counts"])
    RECORD["B"] = {"n": n, "rem_adds_p05": n_rem_stars, "sleep_keeps_p05": n_slp_keep, "rows_top_to_bottom": list(b.label),
                   "printed": [{"label": r.label, "hr": f"{r.rem_hr:.2f} ({r.rem_lo:.2f}-{r.rem_hi:.2f})", "p": fmt_p(r.lr_p), "bold": bool(r.lr_p < .05),
                                "src": {"rem_hr": r.rem_hr, "rem_lo": r.rem_lo, "rem_hi": r.rem_hi, "lr_p": r.lr_p}} for r in b.itertuples()]}
    y = np.arange(n, dtype=float)[::-1]
    COL_HR_IN = RIGHT - 1.52
    COL_P_IN = RIGHT - 0.44
    DW = COL_HR_IN - 0.14 - AX_X0          # the axes width: geometry, unchanged
    COL_HR_IN = COL_HR_IN - 0.08           # R49 nudge: the HR column (values and header) 5.8 pt to the left, so the 11 pt strings keep a 5.5 pt gap to the P column (V26: 6.8 pt at 10 pt); the column still starts 21 pt after the last interval end
    AX_H = n * ROW_IN
    BOT, TOPBAND, LEG_H = 0.64, 0.46, 0.50
    H = LEG_H + BOT + AX_H + TOPBAND
    with plt.rc_context(rc_polish()):
        fig = plt.figure(figsize=(W, H))
        fit_gutter(fig, list(b.label), fontsize=ROWLAB)   # R49 exception
        ax = fig.add_axes(frac(fig, AX_X0, LEG_H + BOT, DW, AX_H))
        bands(ax, n)
        refline(ax)
        tr = blended_transform_factory(fig.transFigure, ax.transData)
        for i in range(n):
            r = b.iloc[i]
            plain_ci(ax, r.rem_lo, r.rem_hi, r.rem_hr, y[i] + 0.20, BLUE)
            plain_ci(ax, r.slp_lo, r.slp_hi, r.slp_hr, y[i] - 0.20, GREY, marker="s", size=S_SQUARE)
            fig.text(COL_HR_IN / W, y[i], f"{r.rem_hr:.2f} ({r.rem_lo:.2f}-{r.rem_hi:.2f})", transform=tr, fontsize=ANNF, va="center", ha="left", color=INK)
            fig.text(COL_P_IN / W, y[i], fmt_p(r.lr_p), transform=tr, fontsize=ANNF, va="center", ha="left", color=INK,
                     fontweight="bold" if r.lr_p < .05 else "normal")
        left_rows(ax, y, list(b.label), size=ROWLAB)   # R49 exception
        bare_y(ax)
        logx(ax, [0.6, 0.8, 1, 1.5, 2, 2.5])
        ax.set_xlim(0.55, 2.92)
        ax.set_ylim(-0.7, n - 0.3)
        axis_covers(ax, np.concatenate([b[c] for c in ("rem_hr", "rem_lo", "rem_hi", "slp_hr", "slp_lo", "slp_hi")]), "x", "eFigure12B")
        ax.set_xlabel("Adjusted hazard ratio per 1 SD, both terms in one model", fontsize=LABF, labelpad=5)
        fig.text(COL_HR_IN / W, n - 0.05, "REM T90, HR (95% CI)", transform=tr, fontsize=ANNF, va="bottom", ha="left", color=INK)
        fig.text(COL_P_IN / W, n - 0.75, "P", transform=tr, fontsize=ANNF, va="bottom", ha="left", color=INK)
        fig.legend(handles=[
            Line2D([], [], color=BLUE, marker="o", ls="-", lw=CI_LW, ms=MS_CIRCLE, markeredgecolor="white", markeredgewidth=MEDGE,
                   label="REM T90, adjusted for whole-sleep T90"),
            Line2D([], [], color=GREY, marker="s", ls="-", lw=CI_LW, ms=MS_SQUARE, markeredgecolor="white", markeredgewidth=MEDGE,
                   label="Whole-sleep T90, adjusted for REM T90")],
            loc="lower center", bbox_to_anchor=((AX_X0 + DW / 2.0) / W, EDGE / H), ncol=1, fontsize=TCKF, frameon=False, handletextpad=0.6, labelspacing=0.45, handlelength=2.2, numpoints=1)   # R40: centred under the axes (where the V13 strip held it), one marker per entry, Arial
        family_letter(fig, "b")
        WROTE.append(save_polished(fig, "eFigure12B_rem_vs_whole_sleep", GATE))
    shown(str(n), n, D["B_counts"]["n_outcomes"], "12B n conditions")
    shown(str(n_rem_stars), n_rem_stars, D["B_counts"]["rem_adds_p05"], "12B rem adds count")
    shown(str(n_slp_keep), n_slp_keep, D["B_counts"]["sleep_keeps_p05"], "12B sleep keeps count")


# =============================================== 12C, healthy subgroup, paired (rules 3, 4, 5)
def fig_c():
    c = pd.DataFrame(D["C"]).sort_values("g_hr", ascending=False).reset_index(drop=True)
    n = len(c)
    starred = sorted(c.disease[c.p_interaction < .05])
    RECORD["C"] = {"n": n, "starred": starred, "rows_top_to_bottom": list(c.disease),
                   "printed": [{"disease": r.disease, "hr": f"{r.h_hr:.2f} ({r.h_lo:.2f}-{r.h_hi:.2f})", "p": fmt_p(r.p_interaction), "bold": bool(r.p_interaction < .05),
                                "src": {"h_hr": r.h_hr, "h_lo": r.h_lo, "h_hi": r.h_hi, "p_interaction": r.p_interaction, "g_hr": r.g_hr, "g_lo": r.g_lo, "g_hi": r.g_hi}} for r in c.itertuples()],
                   "header_right_in": round(C_HDR_RIGHT_IN, 4), "header_baselines_from_top_in": [round(v, 4) for v in C_HDR_Y_IN]}
    y = np.arange(n, dtype=float)[::-1]
    COL_HR_IN = RIGHT - 1.96
    COL_P_IN = RIGHT - 0.78
    DW = COL_HR_IN - 0.14 - AX_X0
    AX_H = n * ROW_IN
    BOT, TOPBAND, LEG_H = 0.64, 0.46, 0.50
    H = LEG_H + BOT + AX_H + TOPBAND
    with plt.rc_context(rc_polish()):
        fig = plt.figure(figsize=(W, H))
        fit_gutter(fig, list(c.disease), fontsize=ROWLAB)   # R49 exception
        ax = fig.add_axes(frac(fig, AX_X0, LEG_H + BOT, DW, AX_H))
        bands(ax, n)
        refline(ax)
        tr = blended_transform_factory(fig.transFigure, ax.transData)
        for i in range(n):
            r = c.iloc[i]
            plain_ci(ax, r.h_lo, r.h_hi, r.h_hr, y[i] + 0.20, BLUE)
            plain_ci(ax, r.g_lo, r.g_hi, r.g_hr, y[i] - 0.20, GREY, marker="s", size=S_SQUARE)
            fig.text(COL_HR_IN / W, y[i], f"{r.h_hr:.2f} ({r.h_lo:.2f}-{r.h_hi:.2f})", transform=tr, fontsize=ANNF, va="center", ha="left", color=INK)
            fig.text(COL_P_IN / W, y[i], fmt_p(r.p_interaction), transform=tr, fontsize=ANNF, va="center", ha="left", color=INK,
                     fontweight="bold" if r.p_interaction < .05 else "normal")
        left_rows(ax, y, list(c.disease), size=ROWLAB)   # R49 exception
        bare_y(ax)
        logx(ax, [0.6, 0.8, 1, 1.5, 2])
        ax.set_xlim(0.55, 2.45)
        ax.set_ylim(-0.7, n - 0.3)
        axis_covers(ax, np.concatenate([c[k] for k in ("g_hr", "g_lo", "g_hi", "h_hr", "h_lo", "h_hi")]), "x", "eFigure12C")
        ax.set_xlabel("Hazard ratio per 1 SD of REM T90", fontsize=LABF, labelpad=5)
        # rule 4a (relaunch 5): the two column headers are NOT drawn inside the sub-sheet. On the round-30 base
        # sheet they are page-level stamps sitting 13.5 and 28.7 pt below the sub-sheet's top edge, inside the
        # sub-sheet's own 4 mm design margin, so drawing them here trips this sheet's own edge gate (it did at
        # 04:31). The compose step stamps them on the page at the base sheet's own origins, font and size, from
        # this record, which is what the base sheet itself does.
        # rule 4: the two column headers right-aligned at the base sheet's anchor, in the top band, baselines as
        # on the base. Exempt from this sheet's 4 mm edge gate (see rule 4b in save_polished): on the base sheet
        # they are page-level stamps at exactly these origins.
        RECORD["C_headers"] = [{"text": "Healthy subgroup, HR (95% CI)", "right_x_in": round(C_HDR_RIGHT_IN, 4),
                                "baseline_from_top_in": round(C_HDR_Y_IN[0], 4), "size": ANNF},
                               {"text": "P interaction", "right_x_in": round(C_HDR_RIGHT_IN, 4),
                                "baseline_from_top_in": round(C_HDR_Y_IN[1], 4), "size": ANNF}]
        for txt, yb in (("Healthy subgroup, HR (95% CI)", C_HDR_Y_IN[0]), ("P interaction", C_HDR_Y_IN[1])):
            fig.text(C_HDR_RIGHT_IN, H - yb, txt, transform=fig.dpi_scale_trans, fontsize=ANNF, va="baseline", ha="right", color=INK)
        # rule 4: the legend centred under the axes (the base moved it there from the lower left), same distance from the bottom
        fig.legend(handles=[
            Line2D([], [], color=BLUE, marker="o", ls="-", lw=CI_LW, ms=MS_CIRCLE, markeredgecolor="white", markeredgewidth=MEDGE, label="Healthy subgroup"),
            Line2D([], [], color=GREY, marker="s", ls="-", lw=CI_LW, ms=MS_SQUARE, markeredgecolor="white", markeredgewidth=MEDGE, label="General cohort")],
            loc="lower center", bbox_to_anchor=((AX_X0 + DW / 2.0) / W, EDGE / H), ncol=1, fontsize=TCKF, frameon=False, handletextpad=0.6, labelspacing=0.45, handlelength=2.2, numpoints=1)   # R40: one marker per entry (explicit)
        family_letter(fig, "c")
        WROTE.append(save_polished(fig, "eFigure12C_healthy_subgroup", GATE))
    shown(str(n), n, len(c), "12C n conditions")
    shown(str(len(starred)), len(starred), int((c.p_interaction < .05).sum()), "12C interaction P under 0.05")


# =============================================== 12D (rules 3, 5)
FLOOR_ORDER = [("all", "Whole recording"), ("wake", "Wake"), ("sleep", "Whole sleep"), ("nrem", "NREM"), ("n1", "N1"), ("n2", "N2"), ("n3", "N3"), ("rem", "REM")]
_NC = json.load(open(f"{paths.NUMBERS_DIR}/negcontrols_final.json"))
KEEP_CONTROLS = list(_NC["panel_keys"])            # v8.1 panel from the file: alopecia, hernia_ing, glaucoma, dermatitis_contact, haemorrhoids
DROPPED_CONTROLS = ["fracture", "osteoarthritis", "back_pain", "cataract"]   # the v7-era controls no longer in the panel
STAGE_SRC = f"{paths.SV_ROOT}/stage_specific"
CHANCE = 0.25
CLAIM_EDGE = 1.25
assert CLAIM_EDGE == 1.25, "author display convention, no frozen source, change on purpose"
ZONE = "#ebedef"
SLEEPONLY_WINDOWS = ["t90_sleep_hr", "t90_nrem_hr", "t90_n2_hr", "t90_rem_hr"]
SLEEPONLY_STAGES = ("sleep", "nrem", "n2", "rem")
CORE_SRC = f"{STAGE_SRC}/stage_general_common_core.csv"
SUMM_SRC = f"{STAGE_SRC}/attack_summary.json"
for _p in (CORE_SRC, SUMM_SRC, f"{STAGE_SRC}/attack_controls_full.csv"):
    _s = os.stat(_p); assert _s.st_size > 0 and _s.st_blocks > 0, ("evicted", _p)
    assert os.path.exists(_p + ".provenance.json"), ("no sidecar", _p)


def control_table():
    frozen = pd.read_csv(f"{STAGE_SRC}/attack_controls_full.csv")
    v5 = pd.read_csv(f"{STAGE_SRC}/attack_controls_full_negpanel_v5.csv")   # step 1463, the v8.1 control panel refit
    RECORD["controls_v5_equal_full"] = bool(set(v5.control) == set(frozen.control) and all(abs(float(v5.set_index("control").loc[k, f"{w}_hr"]) - float(frozen.set_index("control").loc[k, f"{w}_hr"])) < 5e-4 for k in v5.control for w, _l in FLOOR_ORDER))
    if set(KEEP_CONTROLS) <= set(frozen.control):
        return frozen[frozen.control.isin(KEEP_CONTROLS)], "attack_controls_full.csv"
    table = pd.read_csv(f"{STAGE_SRC}/attack_controls_full_negpanel_v5.csv")
    shared = sorted(set(table.control) & set(frozen.control))
    assert shared, "no shared control to check the refit against"
    a, f = table.set_index("control"), frozen.set_index("control")
    for k in shared:
        for w, _lab in FLOOR_ORDER:
            assert abs(float(a.loc[k, f"{w}_hr"]) - float(f.loc[k, f"{w}_hr"])) < 5e-4, (k, w)
    return table[table.control.isin(KEEP_CONTROLS)], "attack_controls_full_negpanel_v5.csv"


def recompute_floors():
    """floor = max(worst single control HR, upper 95% of the pooled control). Rule 5: pinned against the v7
    attack_summary.json floors (same step 46), not against typed numbers."""
    keep, src = control_table()
    assert len(keep) == len(KEEP_CONTROLS), (src, sorted(keep.control))
    assert not set(DROPPED_CONTROLS) & set(keep.control), src
    out = {}
    for k, _lab in FLOOR_ORDER:
        w = 1.0 / keep[f"{k}_se"] ** 2
        cf = float((keep[f"{k}_coef"] * w).sum() / w.sum())
        pooled_hi = float(np.exp(cf + 1.96 * np.sqrt(1.0 / w.sum())))
        worst = float(keep[f"{k}_hr"].max())
        out[k] = max(worst, pooled_hi)
    summ_floors = json.load(open(SUMM_SRC))["floor"]
    for k, v in out.items():
        assert abs(v - float(summ_floors[k])) < 5e-4, (k, v, summ_floors[k])
    RECORD["D_floors"] = {"controls_from": src, "recomputed": out, "attack_summary": {k: float(summ_floors[k]) for k, _ in FLOOR_ORDER}}
    return out, src


def fig_d():
    s = pd.DataFrame(D["D_scatter"])
    core = pd.read_csv(CORE_SRC).set_index("key")
    for _i, r in s.iterrows():
        row = core.loc[r["key"]]
        four = [float(row[c]) for c in SLEEPONLY_WINDOWS]
        assert max(four) == float(row["win_hr_sleeponly"]), (r["key"], four, row["win_hr_sleeponly"])
        assert round(float(r["hr"]), 4) == round(float(row["win_hr_sleeponly"]), 4), (r["key"], r["hr"], row["win_hr_sleeponly"])
        assert r["observed"] == str(row["win_stage_sleeponly"]), (r["key"], r["observed"], row["win_stage_sleeponly"])
        assert r["observed"] in SLEEPONLY_STAGES, (r["key"], r["observed"])
    floors, floor_src = recompute_floors()
    falls = {k: float(D["D_floors"][k]) - floors[k] for k, _lab in FLOOR_ORDER}
    RECORD["D"] = {"n_points": len(s), "points": s.to_dict("records"), "floors_drawn": {k: round(floors[k], 6) for k, _ in FLOOR_ORDER},
                   "floors_printed": {}, "floor_values_3dp": {k: f"{floors[k]:.3f}" for k, _ in FLOOR_ORDER}, "json_D_floors": D["D_floors"],   # R50: nothing printed at the bar ends
                   "falls_json_minus_drawn": {k: round(v, 6) for k, v in falls.items()},
                   "fall_range": [round(min(falls.values()), 3), round(max(falls.values()), 3)]}
    A_X0, A_H = 1.06, 4.72
    AX_X0_D = 1.30
    B_H = len(FLOOR_ORDER) * ROW_IN
    BOT = 0.92
    B_Y0 = BOT
    TITLE_BAND = ((ANNF - PT_PLUS) * 1.25 + 6.0) / 72.0   # R49: the band that sets panel d's axes and the sub-sheet height stays at the round-40 value (geometry pinned)
    A_Y0 = B_Y0 + B_H + 1.02 - TITLE_BAND
    H = A_Y0 + A_H + 0.72 - TITLE_BAND
    with plt.rc_context(rc_polish()):
        fig = plt.figure(figsize=(W, H))
        fit_gutter(fig, [lab for _k, lab in FLOOR_ORDER], ax_x0=AX_X0_D, safety=0.0)   # R49: 'Whole recording' at 11 pt ends 0.5 pt before the axes edge, which carries no ink (nearest bar end 41.7 pt away)
        axA = fig.add_axes(frac(fig, A_X0, A_Y0, RIGHT - A_X0, A_H))
        axB = fig.add_axes(frac(fig, AX_X0_D, B_Y0, RIGHT - AX_X0_D, B_H))
        isrem = (s.observed == "rem").values
        axA.axvspan(0.95, CLAIM_EDGE, color=ZONE, lw=0, zorder=0)
        axA.axvline(CLAIM_EDGE, color=INK, lw=0.9, ls=(0, (1.5, 1.8)), zorder=1)
        axA.axhline(CHANCE, color=INK, lw=0.9, ls=(0, (4, 3)), zorder=1)
        axA.scatter(s.hr[isrem], s.p_obs[isrem], s=S_CIRCLE + 9, marker="o", color=BLUE, zorder=4, edgecolors="white", linewidths=0.7)
        axA.scatter(s.hr[~isrem], s.p_obs[~isrem], s=S_CIRCLE, marker="o", facecolors="white", edgecolors=GREY, linewidths=1.2, zorder=4)
        axA.set_xlim(0.95, 1.86)
        axA.set_ylim(-0.02, 1.13)
        axA.set_xticks([1.0, 1.2, 1.4, 1.6, 1.8])
        axA.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
        axis_covers(axA, s.hr, "x", "eFigure12D a")
        axis_covers(axA, list(s.p_obs) + [CHANCE], "y", "eFigure12D a")
        axA.set_xlabel("Largest observed hazard ratio across the four sleep windows", fontsize=LABF)
        axA.set_ylabel("Probability the same window has the largest\nhazard ratio under patient-level resampling", fontsize=LABF)
        # R40 (Alen, 2026-09-18): the printed "shaded zone, below 1.25: no window-versus-window claim is made" leaves (the shading stays; the rule goes to the figure legend)
        axA.text(1.848, CHANCE + 0.015, "chance", fontsize=ANNF, color=INK, va="bottom", ha="right")
        axA.legend(handles=[Line2D([], [], ls="none", marker="o", ms=6.7, mfc=BLUE, mec="white", mew=0.7, label="REM has the largest hazard ratio"),
                            Line2D([], [], ls="none", marker="o", ms=6.3, mfc="white", mec=GREY, mew=1.2, label="Another sleep stage has the largest hazard ratio")],   # R40: "window" -> "sleep stage"
                   loc="lower right", bbox_to_anchor=(1.0, 0.03), fontsize=TCKF, handletextpad=0.5, labelspacing=0.45, frameon=False, numpoints=1)
        fig.canvas.draw()
        from matplotlib.transforms import Bbox
        _xe = axA.transData.transform((CLAIM_EDGE, 0.0))[0]
        _y0 = axA.transAxes.transform((0.0, 0.0))[1]
        _y1 = axA.transAxes.transform((0.0, 1.0))[1]
        _edge_rule = Bbox([[_xe - 2.5, _y0], [_xe + 2.5, _y1]])
        _xf = axA.transAxes.transform((0.0, 0.0))[0]
        _left_frame = Bbox([[_xf - 2.5, _y0], [_xf + 2.5, _y1]])
        _labels_a = ["Urinary tract\ninfection" if lb == "Urinary tract infection" else lb for lb in s.label]
        assert [lb.replace("\n", " ") for lb in _labels_a] == list(s.label)
        nooverlap.place_labels(axA, s.hr.values, s.p_obs.values, _labels_a, PTLAB, INK, extra_avoid=[_edge_rule, _left_frame])

        y = np.arange(len(FLOOR_ORDER), dtype=float)[::-1]
        for (k, lab), yy in zip(FLOOR_ORDER, y):
            col = BLUE if k == "rem" else (GREY_PALE if k == "n3" else GREY)
            axB.plot([1.0, floors[k]], [yy, yy], color=col, lw=2.4, solid_capstyle="round", zorder=2)
            axB.scatter([floors[k]], [yy], s=33 if k == "n3" else S_CIRCLE, marker="x" if k == "n3" else "o", color=col, zorder=4,
                        edgecolors=None if k == "n3" else "white", linewidths=1.4 if k == "n3" else MEDGE)
            # R50 (Alen, 2026-09-26, item 6): the value at the end of each bar is NOT printed any more (the bars and the axis carry it);
            # the floor still equals the attack_summary floor within 5e-4 (recompute_floors) and is recorded, not shown
            assert abs(round(floors[k], 3) - round(float(RECORD["D_floors"]["attack_summary"][k]), 3)) < 1e-9, (k, floors[k])
        refline(axB)
        left_rows(axB, y, [lab for _k, lab in FLOOR_ORDER], ax_x0=AX_X0_D)
        bare_y(axB)
        # v8.2: the axis start is data-driven (the V13 axis started at 0.998; a floor below that extends the axis to the left, with a 0.95 tick)
        _xlo = 0.998 if min(floors.values()) >= 0.998 else min(0.998, min(floors.values()) - 0.03)
        RECORD["D_axis"] = {"xlim": [round(_xlo, 4), 1.245], "v13_xlim": [0.998, 1.245], "extended": bool(_xlo < 0.998)}
        axB.set_xlim(_xlo, 1.245)
        axB.set_ylim(-0.7, len(FLOOR_ORDER) - 0.3)
        axB.set_xticks(([0.95] if _xlo < 0.95 + 1e-9 else []) + [1.00, 1.05, 1.10, 1.15, 1.20])
        axB.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:.2f}"))
        axis_covers(axB, [1.0] + list(floors.values()), "x", "eFigure12D b")
        axB.set_xlabel("Level the negative controls reach, hazard ratio per 1 SD", fontsize=LABF)
        n3_key_legend(fig, 0.05, H)
        family_letter(fig, "d")
        WROTE.append(save_polished(fig, "eFigure12D_ranking_worth_and_floor", GATE))
    shown(str(len(s)), len(s), len(D["D_scatter"]), "12D n scatter conditions")


if __name__ == "__main__":
    fig_a()
    fig_b()
    fig_c()
    if "--skip-d" not in sys.argv: fig_d()   # --old positive control: panel d reads the CURRENT core csv, so it is skipped there
    else: RECORD["D_skipped"] = "positive control run on the V13 (v7) values: panel d is data-gated on the current core csv and skipped"
    RECORD.update({"figures": [os.path.basename(p) for p in WROTE], "printed_numbers_audited": SHOWN, "overlap_gate": GATE})
    json.dump(RECORD, open(f"{WORK}/build_record.json", "w"), indent=1, default=float)
    print("\n".join(WROTE))
    print(f"mode {MODE}: A n {RECORD['A']['n']} | B n {RECORD['B']['n']} ({RECORD['B']['rem_adds_p05']} of {RECORD['B']['n']} REM adds, {RECORD['B']['sleep_keeps_p05']} sleep keeps) | "
          f"C n {RECORD['C']['n']} starred {RECORD['C']['starred']} | D {RECORD.get('D', {}).get('n_points', 'skipped')} points floors {RECORD.get('D', {}).get('floor_values_3dp')}")
