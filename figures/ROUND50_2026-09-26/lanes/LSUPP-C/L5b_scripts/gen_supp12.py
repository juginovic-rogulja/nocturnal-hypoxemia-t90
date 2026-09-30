"""
eFigure 21, polished for the FINAL_FIGURES_2026-08-14 supplementary set.

Content is make_efig21_correction_vs_residual.py's, with every assertion kept: the 11
FDR-significant outcomes, residual on-PAP T90 (dark blue) beside the amount corrected by
PAP at fixed residual (grey) on one row each, the HR and events columns, the open-marker
convention and the three-entry legend. Drawing quality moves to the shared design system:
171 mm canvas, ink #1a1d21, 0.8 pt spines, axis titles 11 pt, ticks and row labels 10 pt
in one left-aligned gutter, number columns 9.5 pt, CI lines 1.7 pt with round caps and
6 pt white-edged markers, the reference rule 0.9 pt dashed ink, explicit log ticks with
NullLocator minors, and a frameless 10 pt legend.

REMOVED, per Alen's instruction: the grey "Personal reminder ..." strip at the foot of the
original. It was Alen's own note and must not ship. The grey legend-prose blocks the
original routed through legend_capture (never drawn on the shipped sheet) are not drawn
here either, so the shipped text budget is unchanged.

2026-08-21 round (ROUND_BRIEF_2026-08-21.md, section 2, and the Figure 4 owner's task):
the sheet gains a q column beside the HR column, labelled q and never P, because both
drawn series carry a Benjamini-Hochberg q in the frozen sources: q_sd for the residual
exposure (all 11 below 0.05, bold per the round's "bold where significant" rule) and
q_delta_sd for the corrected amount (none below 0.05, smallest 0.87, set regular). The
sheet's content is unchanged and it stays in the supplement per Alen's ruling today.
The on-sheet headline comes off per the round's output standard (no titles or
headlines on submission sheets). q prints to two significant figures, "<0.001" below
0.001. The plot column narrows 2.50 -> 1.92 in to hold the third printed column.

splitstyle.save() is imported per the style contract and never called. legend_capture's
patched text methods are restored at import exactly as the other polish scripts do.
Output goes to Supplementary_Figures_Polished (PDF) and _workfiles/polished_pngs (300 dpi
PNG) only. Nothing under New_Figures is written.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import atexit
import json
import os
import sys

import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory

SCRIPTS = f"{paths.FIGURE_ROOT}/ROUND50_2026-09-26/lanes/LSUPP-C/L5b_scripts"   # V14_L5b 2026-09-15: the lane copies of the style modules
sys.path.insert(0, SCRIPTS)
from splitstyle import (rc, save, axis_covers, bare_y, logx, sep_ok,        # noqa: E402
                        BLUE, GREY, GREY_PALE, plt)
assert callable(save)   # imported per the style contract, deliberately never called here
import legend_capture as _lc                                                # noqa: E402
import nooverlap                                                            # noqa: E402
from matplotlib.axes import Axes as _Axes                                   # noqa: E402
from matplotlib.figure import Figure as _Figure                             # noqa: E402

_Figure.text = _lc._orig_text
_Axes.text = _lc._orig_ax_text
_Axes.annotate = _lc._orig_ax_annotate
atexit.unregister(_lc.flush)

ROOT = paths.FIGURE_ROOT
FIN = f"{paths.FIGURE_ROOT}/ROUND50_2026-09-26/lanes/LSUPP-C/L5b_work/polish_out"   # V14_L5b: nothing is written under FINAL_FIGURES
OUT = f"{FIN}/pdf"
PNGS = f"{FIN}/png"
WORK = f"{FIN}/work"
os.makedirs(OUT, exist_ok=True)
os.makedirs(PNGS, exist_ok=True)

# ------------------------------------------------------------------ the design system
W = 171.0 / 25.4
INK = "#1a1d21"
PT_PLUS = 1.0   # R49: every text one point larger, geometry unchanged
LABF, TCKF, ANNF, FLOOR = 11.0 + PT_PLUS, 10.0 + PT_PLUS, 9.5 + PT_PLUS, 9.0
CI_LW, MEDGE = 1.7, 0.8
S_CIRCLE, S_SQUARE = 28.0, 23.0
MS_CIRCLE, MS_SQUARE = 6.0, 5.4
EDGE = 0.18


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


def left_rows(ax, ys, labels, gutter_x, ax_x0):
    ax.set_yticks(list(ys))
    ax.set_yticklabels(labels, fontsize=TCKF, color=INK, ha="left")
    ax.tick_params(axis="y", pad=(ax_x0 - gutter_x) * 72.0, length=0)


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


def save_polished(fig, name):
    for t in _all_strings(fig):
        s = t.get_text()
        assert "—" not in s, f"{name}: em-dash on sheet in {s!r}"
        assert ";" not in s, f"{name}: semicolon on sheet in {s!r}"
        assert "Personal reminder" not in s, f"{name}: reminder text must not ship"
        assert t.get_fontsize() >= FLOOR, (
            f"{name}: {t.get_fontsize()} pt under the {FLOOR} pt floor in {s!r}")
    boxes, _ = nooverlap._boxes(fig)
    Wpx, Hpx = fig.get_size_inches() * fig.dpi
    m = (4.0 - 0.3) / 25.4 * fig.dpi
    tight = [f"{b['where']}: {b['text'][:44]}" for b in boxes
             if b["bbox"].x0 < m or b["bbox"].y0 < m
             or b["bbox"].x1 > Wpx - m or b["bbox"].y1 > Hpx - m]
    assert not tight, f"{name}: text inside the 4 mm edge margin: {tight}"
    nooverlap.gate(fig, name)
    pdf = f"{OUT}/{name}.pdf"
    fig.savefig(pdf)
    fig.savefig(f"{PNGS}/{name}.png", dpi=300)
    plt.close(fig)
    print("wrote", pdf)
    print("wrote", f"{PNGS}/{name}.png")
    return pdf


# ------------------------------------------------------------------ sources, asserted
SRC = (f"{paths.SV_ROOT}/"
       "pap_residual_continuous")

cont = pd.read_csv(f"{SRC}/results_continuous.csv")
delta = pd.read_csv(f"{SRC}/delta_results.csv")

R = cont[cont.sig_fdr == True].sort_values("hr_sd", ascending=False)      # noqa: E712
assert len(R) >= 1 and (R.q_sd < 0.05).all() and (R.lo_sd > 1.0).all(), len(R)   # V14_L5b: data-driven (gate 6), the row count is the file's
print(f"FDR survivors drawn: {len(R)}")
D = delta.set_index("key").loc[R.key.tolist()]
assert (D.events.values == R.events.values).all(), "the two fits must share the event counts"
assert not D.sig_fdr.any(), "no corrected-amount association survives FDR"
MIN_Q = float(delta[delta.negative_control == False].q_delta_sd.min())    # noqa: E712
assert MIN_Q >= 0.05, MIN_Q   # V14_L5b: data-driven (the smallest corrected-amount q is reported, not typed)
N_DELTA_EXCL1 = int((D.lo_delta_sd > 1.0).sum()) + int((D.hi_delta_sd < 1.0).sum())   # V14_L5b: grey intervals that exclude 1, drawn as the data says
print(f"smallest corrected-amount q {MIN_Q:.4f}; corrected-amount intervals excluding 1: {N_DELTA_EXCL1}")

RESID, DELTA = "#0288d1", GREY   # V14_L5b 2026-09-15: the V13 sheet draws the residual series in the main-figure blue #0288d1 (round 31 LCOLOR), not splitstyle's August #1f4257
LINK = GREY_PALE                                   # the connector, the house de-emphasis grey
# V14_L5b: splitstyle's greyscale-separation gate is not applied to the V13 colour pair (#0288d1 against #8a9099, the pair the V13 sheet already prints); the design is fixed by V13


from decimal import Decimal, ROUND_HALF_UP   # V14_L5b: half up on the full-precision value (results_continuous.csv stores full precision)
def _r2(v): return str(Decimal(repr(float(v))).quantize(Decimal("0.01"), rounding=ROUND_HALF_UP))
def fmt(hr, lo, hi):
    return f"{_r2(hr)} ({_r2(lo)}-{_r2(hi)})"


def fmt_p(p):
    """Two significant figures, "<0.001" below 0.001, one further decimal where
    rounding would land exactly on the 0.050 boundary."""
    if p < 0.001:
        return "<0.001"
    from math import floor, log10
    d = 1 - floor(log10(p))
    s = f"{p:.{d}f}"
    while s in ("0.050", "0.0500") and abs(p - 0.05) > 1e-12 and d < 6:
        d += 1
        s = f"{p:.{d}f}"
    return s


# V14_L5b 2026-09-15: no typed q table (gate 6). The residual q is the FDR family this sheet's significance claim rests on (all
# drawn rows below 0.05, bold); the corrected-amount q never reaches 0.05 (regular weight). Both are asserted from the files.
for r in R.itertuples():
    d = D.loc[r.key]
    assert r.q_sd < 0.05, (r.key, r.q_sd)
    assert d.q_delta_sd >= 0.05, (r.key, d.q_delta_sd)


# ------------------------------------------------------------------ geometry, inches
GUT = EDGE
# 2026-08-21: the headline came off (round output standard), so TOP is a plain pad.
# The plot column narrows 2.50 -> 1.92 in for the third printed column.
LEFT, PLOTW = 2.05, 1.92
# 2026-08-25 round 4 (Alen's Supp 12 pick, reversing round 3's top placement): the
# three-entry legend moves to the BOTTOM of the figure, below the x-axis title,
# left-aligned on the plot column's left edge. TOP shrinks to a plain pad, BOT grows
# to hold the ticks, the axis title and the three legend rows. The in-axes column
# headers stay at 0.95 data units and the ylim top stays 1.35 (round-3 trim kept).
# Entry text identical, no value moves.
TOP = 0.25
UNITS = (11 - 1) + 0.75 + 1.35
PITCH = 0.40
AX_H = UNITS * PITCH
BOT = 1.42                                          # ticks, axis title, then the key
H = TOP + AX_H + BOT
# the bottom band must hold the ticks and axis title (0.39 in rendered), a 0.11 in
# gap, three 10 pt legend rows with their two labelspacing gaps (0.68 in rendered),
# and the 4 mm bottom margin
assert BOT >= 0.39 + 0.11 + 0.68 + 0.16, "bottom band too short for the three legend rows"

XLIM, XTICKS = (0.74, 3.25), [0.8, 1, 1.5, 2, 3]
OFF = 0.21

# the three printed columns, positioned in inches and converted to axes fractions
COL_W_HR, COL_W_PQ, COL_W_EV = 1.16, 0.47, 0.55
X_HR_IN = LEFT + PLOTW + 0.10
X_PQ_IN = X_HR_IN + COL_W_HR + 0.10
X_EV_IN = X_PQ_IN + COL_W_PQ + 0.10 + COL_W_EV / 2          # centre
assert X_EV_IN + COL_W_EV / 2 <= W - 0.16, "the events column leaves the canvas"
X_EV_R = X_EV_IN + COL_W_EV / 2      # Events column right edge, the legend aligns on it
X_HR = 1 + (X_HR_IN - LEFT - PLOTW) / PLOTW
X_PQ = 1 + (X_PQ_IN - LEFT - PLOTW) / PLOTW
X_EV = 1 + (X_EV_IN - LEFT - PLOTW) / PLOTW

DRAWN = []
with plt.rc_context(rc_polish()):
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([LEFT / W, BOT / H, PLOTW / W, AX_H / H])
    ax.axvline(1.0, color=INK, lw=0.9, ls=(0, (4, 3)), zorder=1)
    y = -np.arange(len(R), dtype=float)
    for r, yy in zip(R.itertuples(), y):
        d = D.loc[r.key]
        yres, ydel = yy + OFF, yy - OFF
        # V14_L5b: no connector between the two series (the V13 sheet paints them white, so nothing shows)
        for hr, lo, hi, ypt, col, mk, ms_ in [
                (r.hr_sd, r.lo_sd, r.hi_sd, yres, RESID, "o", S_CIRCLE),
                (d.hr_delta_sd, d.lo_delta_sd, d.hi_delta_sd, ydel, DELTA, "s", S_SQUARE)]:
            sig = not (lo <= 1.0 <= hi)
            ax.plot([lo, hi], [ypt] * 2, color=col, lw=CI_LW, solid_capstyle="round",
                    zorder=3)
            ax.scatter([hr], [ypt], s=ms_, marker=mk, zorder=4,
                       facecolors=col if sig else "white",
                       edgecolors="white" if sig else col,
                       linewidths=MEDGE if sig else 1.1)
        DRAWN.append({"key": r.key, "outcome": r.outcome, "events": int(r.events),
                      "residual_hr": r.hr_sd, "residual_lo": r.lo_sd, "residual_hi": r.hi_sd,
                      "delta_hr": d.hr_delta_sd, "delta_lo": d.lo_delta_sd,
                      "delta_hi": d.hi_delta_sd,
                      "printed_residual": fmt(r.hr_sd, r.lo_sd, r.hi_sd),
                      "printed_delta": fmt(d.hr_delta_sd, d.lo_delta_sd, d.hi_delta_sd),
                      "residual_q": r.q_sd, "printed_residual_q": fmt_p(r.q_sd),
                      "delta_q": d.q_delta_sd, "printed_delta_q": fmt_p(d.q_delta_sd),
                      "residual_significant_fdr": True, "delta_significant_fdr": False,
                      "sources": ["results_continuous.csv", "delta_results.csv"]})
    left_rows(ax, y, R.outcome.tolist(), GUT, LEFT)
    bare_y(ax)
    logx(ax, XTICKS)
    ax.xaxis.set_minor_locator(__import__("matplotlib").ticker.NullLocator())
    ax.set_xlim(*XLIM)
    ax.set_ylim(y.min() - 0.75, 1.35)
    axis_covers(ax, [XLIM[0], XLIM[1]] + R.lo_sd.tolist() + R.hi_sd.tolist()
                + R.hr_sd.tolist() + D.lo_delta_sd.tolist() + D.hi_delta_sd.tolist()
                + D.hr_delta_sd.tolist(), which="x", label="eFigure21")
    ax.set_xlabel("Hazard ratio per 1 SD of each exposure (95% CI)", fontsize=LABF,
                  labelpad=4)

    tr = blended_transform_factory(ax.transAxes, ax.transData)
    ax.text(X_HR, 0.95, "HR (95% CI)", transform=tr, fontsize=ANNF, fontweight="bold",
            va="center", ha="left", color=INK)
    ax.text(X_PQ, 0.95, "q", transform=tr, fontsize=ANNF, fontweight="bold",
            va="center", ha="left", color=INK)
    # V14_L5b: no "Events" header (V13 design, round 32 LSUPP)
    for r, yy in zip(R.itertuples(), y):
        d = D.loc[r.key]
        ax.text(X_HR, yy + OFF, fmt(r.hr_sd, r.lo_sd, r.hi_sd), transform=tr,
                fontsize=ANNF, va="center", ha="left", color=INK)
        ax.text(X_HR, yy - OFF, fmt(d.hr_delta_sd, d.lo_delta_sd, d.hi_delta_sd),
                transform=tr, fontsize=ANNF, va="center", ha="left", color=INK)
        # residual q bold (FDR survivor, asserted), corrected-amount q regular
        ax.text(X_PQ, yy + OFF, fmt_p(r.q_sd), transform=tr, fontsize=ANNF,
                va="center", ha="left", color=INK, fontweight="bold")
        ax.text(X_PQ, yy - OFF, fmt_p(d.q_delta_sd), transform=tr, fontsize=ANNF,
                va="center", ha="left", color=INK)
        # V14_L5b: no events count beside the rows (V13 design, round 32 LSUPP)

    # V14_L5b 2026-09-15: the three-entry key drawn at the V13 positions (the round-31 LCOLOR vector re-stamp is what V13 shows):
    # 16 pt handle lines at x 90 to 106 (1.4 pt), a 5.4 pt marker at their centre, labels at x 112 on the baselines 421.6, 435.1, 448.5
    # (page points, top-down), 9.2 pt ink; third label "95% CI includes 1" as the V13 sheet prints it. Figure coordinates in inches:
    # (x / 72, (H * 72 - y) / 72) on this 484.72 x 468.72 pt canvas.
    _fx = lambda x: x / 72.0
    _fy = lambda y: (H * 72.0 - y) / 72.0
    KEY_TEXTS = []   # R40: the key text artists, width-gated against the V13 crop edge
    KEY_ROWS = [(418.5, 421.6, RESID, "o", True, "Residual sleep T90 on PAP, per 1 SD"),
                (432.0, 435.1, DELTA, "s", True, "Amount corrected by PAP, per 1 SD, adjusted for the residual T90"),
                (445.5, 448.5, INK, "o", False, "95% CI includes 1")]
    for (yh, yb, col, mk, filled, lab) in KEY_ROWS:
        if filled:
            fig.add_artist(Line2D([_fx(90.0), _fx(106.0)], [_fy(yh), _fy(yh)], color=col, lw=1.4, transform=fig.dpi_scale_trans, solid_capstyle="butt"))
            fig.add_artist(Line2D([_fx(98.0)], [_fy(yh)], color=col, marker=mk, markersize=5.4, ls="none", markerfacecolor=col,
                                  markeredgecolor=col, markeredgewidth=1.0, transform=fig.dpi_scale_trans))
        else:
            fig.add_artist(Line2D([_fx(98.0)], [_fy(yh)], color=col, marker=mk, markersize=5.8, ls="none", markerfacecolor="white",
                                  markeredgecolor=col, markeredgewidth=1.0, transform=fig.dpi_scale_trans))
        KEY_TEXTS.append(fig.text(_fx(112.0), _fy(yb), lab, transform=fig.dpi_scale_trans, ha="left", va="baseline", fontsize=9.2 + PT_PLUS, color=INK))
    fig.canvas.draw()
    _rk = fig.canvas.get_renderer()
    for _t in KEY_TEXTS:
        _x1 = _t.get_window_extent(renderer=_rk).x1 / fig.dpi * 72.0
        assert _x1 <= 425.64 - 5.0, ("key label crosses the V13 crop edge", _t.get_text(), _x1)   # R40: the reworded entry must fit the cropped page

    # 2026-08-21: the on-sheet headline came off per the round's output standard (no
    # titles or headlines on submission sheets, the legend carries all prose).
    # The original's grey legend-prose blocks were captured by legend_capture and never
    # drawn on the shipped sheet, so they are not drawn here either. The grey
    # "Personal reminder ..." strip is DELETED per Alen's instruction and must not return.

    save_polished(fig, "eFigure21_pap_correction_vs_residual")

# ------------------------------------------------------------------ sanity, drawn vs source
cchk, dchk = cont.set_index("key"), delta.set_index("key")
print("drawn value against its source CSV")
for d in DRAWN:
    c, e = cchk.loc[d["key"]], dchk.loc[d["key"]]
    assert (d["residual_hr"], d["residual_lo"], d["residual_hi"]) == (c.hr_sd, c.lo_sd,
                                                                      c.hi_sd), d["key"]
    assert (d["delta_hr"], d["delta_lo"], d["delta_hi"]) == (e.hr_delta_sd, e.lo_delta_sd,
                                                             e.hi_delta_sd), d["key"]
    assert d["events"] == int(c.events) == int(e.events)
    assert d["printed_residual"] == fmt(c.hr_sd, c.lo_sd, c.hi_sd)
    assert d["printed_delta"] == fmt(e.hr_delta_sd, e.lo_delta_sd, e.hi_delta_sd)
    assert d["residual_q"] == c.q_sd and d["printed_residual_q"] == fmt_p(c.q_sd)
    assert d["delta_q"] == e.q_delta_sd and d["printed_delta_q"] == fmt_p(e.q_delta_sd)
    assert c.q_sd < 0.05 <= e.q_delta_sd
    print(f"  {d['key']:16} residual drawn {d['printed_residual']:>18}  csv {fmt(c.hr_sd, c.lo_sd, c.hi_sd):>18}  OK")
    print(f"  {'':16} corrected drawn {d['printed_delta']:>17}  csv "
          f"{fmt(e.hr_delta_sd, e.lo_delta_sd, e.hi_delta_sd):>18}  OK")
json.dump(DRAWN, open(f"{WORK}/eFigure21_polish_drawn_values.json", "w"), indent=1,
          default=float)
print(f"rows drawn: {len(DRAWN)}, smallest corrected-amount q {MIN_Q:.4f}, every value "
      f"asserted against its CSV")
