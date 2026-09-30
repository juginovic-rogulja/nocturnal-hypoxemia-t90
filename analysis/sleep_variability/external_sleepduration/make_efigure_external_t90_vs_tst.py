"""
Draft eFigure, the sleep-duration analogue of main Figure 5D: nocturnal oxygen and total
sleep time entered together, on one per-1-SD scale, for every external outcome.

WHAT IT SHOWS. Twelve cohort-outcome pairs, 7 in the Sleep Heart Health Study and 5 in the
Osteoporotic Fractures in Men Study. Each pair carries two estimates from ONE mutually
adjusted Cox model: nocturnal oxygen adjusted for total sleep time, and total sleep time
adjusted for nocturnal oxygen. Both exposures are rank-based inverse-normal transformed
within cohort, and within site in MrOS, so both hazard ratios are per 1 SD, and both sit on
one axis. Higher total sleep time means longer sleep, so a hazard ratio below 1 on the grey
series means longer sleep with lower risk.

DESIGN. Conventions copied from the delivered Nature Medicine package, chiefly
NATURE_MED_VERSION/_scripts (supp_polish_common.py for the sheet system) and
FINAL_FIGURES_2026-08-14/_scripts/figure5D_polish.py for this exact forest:
  canvas    171 mm wide exactly, savefig.bbox standard, height free, 4 mm margin gate
  type      Arial, axis title 11 pt, ticks and legend 10 pt, printed column 9.5 pt,
            floor 9 pt, bold only for the headline
  colour    ink #1a1d21, house blue #1f4257 for oxygen, house grey #8a9099 for total
            sleep time, greyscale separation asserted at build time
  axes      top and right spines off, left spine hidden because the y axis is categorical
            and the row labels sit in the gutter (bare_y, the convention of every forest
            in the package), 0.8 pt bottom spine, 3 pt outward ticks, log x with explicit
            ticks and no minors, faint grid at the ticks
  marks     CI lines 1.7 pt round caps, 6 pt markers with 0.8 pt edges, filled where the
            estimate survives Benjamini-Hochberg correction across the 12 pairs within its
            exposure family, open where it does not, HR = 1 rule 0.9 pt dashed behind
Single panel, so there is no panel letter to case.

PROVENANCE. Every drawn and printed number is read from results_persd.csv (written by
run_external_sleepduration_persd.py, whose 48-estimate positive control against the frozen
numbers/external.json must pass before it writes anything) and is additionally asserted
against a pinned copy below, so a silent drift in the analysis aborts the sheet rather than
redrawing it.

This script writes ONLY into this folder. It calls neither splitstyle.save nor
legend_capture.flush, so nothing in New_Figures, FINAL_FIGURES, JAMA_VERSION or
NATURE_MED_VERSION is touched.

Outputs: eFigure_external_t90_vs_tst.pdf, eFigure_external_t90_vs_tst.png (300 dpi),
         eFigure_external_t90_vs_tst_drawn_values.json
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import atexit
import json
import os
import sys

import numpy as np
import pandas as pd

SCRIPTS = (f"{paths.FIGURE_ROOT}/New_Figures/"
           "_NOT_THESE_workfiles/scripts")
sys.path.insert(0, SCRIPTS)
sys.dont_write_bytecode = True

from splitstyle import rc, BLUE, GREY, bare_y, logx, sep_ok, axis_covers, plt  # noqa: E402
import legend_capture as _lc                                                  # noqa: E402
import nooverlap                                                              # noqa: E402
from prism_plotter.journal import preflight_figure, journal_style             # noqa: E402
from matplotlib.axes import Axes as _Axes                                     # noqa: E402
from matplotlib.figure import Figure as _Figure                               # noqa: E402
from matplotlib.lines import Line2D                                           # noqa: E402
from matplotlib.ticker import NullLocator                                     # noqa: E402
import matplotlib.text as _mtext                                              # noqa: E402

# restore the text methods legend_capture patched and disarm its atexit flush, which would
# otherwise append into the frozen New_Figures/FIGURE_LEGENDS.md
_Figure.text = _lc._orig_text
_Axes.text = _lc._orig_ax_text
_Axes.annotate = _lc._orig_ax_annotate
atexit.unregister(_lc.flush)

HERE = os.path.dirname(os.path.abspath(__file__))
NAME = "eFigure_external_t90_vs_tst"

MM = 1 / 25.4
W = 171.0 * MM                       # the package width, exactly
EDGE = 0.16                          # left and right margin, inches
INK = "#1a1d21"
GRID = "#eef0f1"
HEAD_PT, TITLE_PT, TICK_PT, ANN_PT, FLOOR_PT = 11.5, 11.0, 10.0, 9.5, 9.0
CI_LW, MEDGE, S_MARK = 1.7, 0.8, 36.0

assert BLUE == "#1f4257" and GREY == "#8a9099", (BLUE, GREY)
sep_ok("oxygen against total sleep time", [BLUE, GREY])

# ---------------------------------------------------------------- the pinned estimates
# (cohort, outcome, oxygen (hr, lo, hi, q), total sleep time (hr, lo, hi, q))
PIN = [
    ("shhs", "Heart failure", (1.191, 1.086, 1.306, 0.0013), (0.983, 0.898, 1.075, 0.9267)),
    ("shhs", "Myocardial infarction", (1.089, 0.961, 1.234, 0.2148),
     (1.119, 0.992, 1.263, 0.3990)),
    ("shhs", "Stroke", (1.127, 0.977, 1.300, 0.1347), (1.025, 0.896, 1.173, 0.9267)),
    ("shhs", "Cardiovascular composite", (1.171, 1.100, 1.247, 0.0000),
     (0.991, 0.933, 1.053, 0.9267)),
    ("shhs", "Coronary heart disease", (1.121, 1.037, 1.211, 0.0098),
     (1.026, 0.952, 1.106, 0.8683)),
    ("shhs", "Cardiovascular death", (1.169, 1.044, 1.310, 0.0141),
     (0.953, 0.853, 1.064, 0.8683)),
    ("shhs", "Atrial fibrillation", (1.021, 0.906, 1.152, 0.7294),
     (0.996, 0.889, 1.116, 0.9743)),
    ("mros", "Death from any cause", (1.088, 1.035, 1.144, 0.0027),
     (0.942, 0.900, 0.987, 0.1348)),
    ("mros", "Cardiovascular death", (1.153, 1.062, 1.252, 0.0027),
     (0.940, 0.870, 1.015, 0.4512)),
    ("mros", "Pulmonary death", (1.223, 1.043, 1.435, 0.0229),
     (0.904, 0.776, 1.052, 0.5743)),
    ("mros", "Stroke death", (1.060, 0.846, 1.327, 0.6682), (1.003, 0.819, 1.230, 0.9743)),
    ("mros", "Cancer death", (1.110, 0.994, 1.238, 0.0955), (1.037, 0.938, 1.147, 0.8683)),
]
ROWLAB = {"Heart failure": "Heart failure",
          "Myocardial infarction": "Myocardial infarction",
          "Stroke": "Stroke",
          "Cardiovascular composite": "CVD composite",
          "Coronary heart disease": "Coronary heart disease",
          "Cardiovascular death": "CVD death",
          "Atrial fibrillation": "Atrial fibrillation",
          "Death from any cause": "Death, any cause",
          "Pulmonary death": "Pulmonary death",
          "Stroke death": "Stroke death",
          "Cancer death": "Cancer death"}
COHLAB = {"shhs": "SHHS", "mros": "MrOS"}

# ---------------------------------------------------------------- read and verify
R = pd.read_csv(os.path.join(HERE, "results_persd.csv"))
FIT = R[R.model == "S2_persd_adjT90"]
assert len(FIT) == 24, len(FIT)


def fetch(coh, out, term):
    sub = FIT[(FIT.cohort == coh) & (FIT.outcome == out) & (FIT.term == term)]
    assert len(sub) == 1, (coh, out, term, len(sub))
    r = sub.iloc[0]
    return float(r.hr), float(r.lo), float(r.hi), float(r.q), int(r.events), int(r.n)


ROWS, DRAWN = [], []
for coh, out, pt90, ptst in PIN:
    got_t = fetch(coh, out, "t90_z")
    got_u = fetch(coh, out, "tst_z")
    for got, pinned, who in ((got_t, pt90, "oxygen"), (got_u, ptst, "total sleep time")):
        hr, lo, hi, q = got[:4]
        assert (round(hr, 3), round(lo, 3), round(hi, 3)) == pinned[:3], \
            f"{coh} {out} {who}: CSV {(hr, lo, hi)} != pinned {pinned[:3]}"
        assert round(q, 4) == pinned[3], f"{coh} {out} {who}: q {q} != pinned {pinned[3]}"
    assert got_t[4] == got_u[4] and got_t[5] == got_u[5], (coh, out)
    ROWS.append({"cohort": coh, "outcome": out, "events": got_t[4], "n": got_t[5],
                 "t90": got_t[:4], "tst": got_u[:4]})

# the claim the sheet makes, asserted before anything is drawn
N_T90_FDR = sum(1 for r in ROWS if r["t90"][3] < 0.05)
N_TST_FDR = sum(1 for r in ROWS if r["tst"][3] < 0.05)
N_TST_RAW = int((FIT[FIT.term == "tst_z"].p < 0.05).sum())
assert (N_T90_FDR, N_TST_FDR, N_TST_RAW) == (7, 0, 1), (N_T90_FDR, N_TST_FDR, N_TST_RAW)
COH_N = {"shhs": 5802, "mros": 2911}
for r in ROWS:
    assert r["n"] <= COH_N[r["cohort"]], r
assert int(R[R.model == "T0_t90_alone_frozen"].shape[0]) == 12

# ---------------------------------------------------------------- geometry
ROW_IN = 0.56
OFF = 0.21                    # sub-row offset, in row units, as in figure5D_polish
AX_X0 = 1.86                  # gutter holds "SHHS" over the longest outcome name
GUT_X = EDGE
COL_FS, COL_GAP = ANN_PT, 0.12
COL_HEAD = "HR (95% CI)"
XTICKS = [0.8, 0.9, 1.0, 1.2, 1.5]
XLIM = (0.75, 1.52)

yD = np.arange(len(ROWS), dtype=float)[::-1]
yD[7:] -= 0.60                # a breath between the two cohorts
SPAN = (yD.max() + 0.72) - (yD.min() - 0.72)
AX_H = round(ROW_IN * SPAN, 2)

HEADLINE = ("Nocturnal oxygen holds in the external cohorts with total sleep time in the "
            "model")
FOOT = ("Both exposures are rank-based inverse-normal transformed within cohort, and within "
        "site in the Osteoporotic Fractures in Men Study, so each hazard ratio is per 1 SD "
        "and a value below 1 on the grey series means longer sleep with lower risk. Both "
        "estimates in a row come from one Cox model carrying the covariates of the frozen "
        "external-validation table, a natural cubic spline on age with 4 df plus sex in the "
        "Sleep Heart Health Study, stratified by clinical site in the Osteoporotic Fractures "
        "in Men Study. Filled markers denote q<0.05 after Benjamini-Hochberg correction "
        "across the 12 cohort-outcome pairs within each exposure family, open markers do "
        "not. Sleep Heart Health Study n = 5,802, Osteoporotic Fractures in Men Study "
        "n = 2,911, the second on follow-up from the sleep visit.")

BOT_AXIS = 0.52               # x axis title and tick labels
HEAD_BAND = 0.30              # the printed-column header
LEG_BAND = 0.64               # the two series and the open-marker key
TOP_BAND = 0.46               # the bold headline

RC = dict(rc())
RC["savefig.bbox"] = "standard"
RC.update({"font.size": TICK_PT, "axes.labelsize": TITLE_PT, "xtick.labelsize": TICK_PT,
           "ytick.labelsize": TICK_PT, "legend.fontsize": TICK_PT,
           "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
           "xtick.major.size": 3.0, "ytick.major.size": 3.0,
           "axes.edgecolor": INK, "text.color": INK, "axes.labelcolor": INK,
           "xtick.color": INK, "ytick.color": INK})


def hr_ci(hr, lo, hi):
    return f"{hr:.2f} ({lo:.2f}-{hi:.2f})"


with plt.rc_context(RC):
    # ---- measure the footnote so the canvas height is derived, never guessed
    probe = plt.figure(figsize=(W, 4.0))
    rend = probe.canvas.get_renderer()

    def width_in(fig, s, size):
        t = fig.text(0.5, 0.5, s, fontsize=size)
        w = t.get_window_extent(renderer=fig.canvas.get_renderer()).width / fig.dpi
        t.remove()
        return w

    FOOT_W = W - 2 * EDGE
    words, line, foot_lines = FOOT.split(), "", []
    for wd in words:
        trial = (line + " " + wd).strip()
        if width_in(probe, trial, ANN_PT) > FOOT_W and line:
            foot_lines.append(line)
            line = wd
        else:
            line = trial
    foot_lines.append(line)
    FOOT_TXT = "\n".join(foot_lines)
    COL_W = max(width_in(probe, s, COL_FS)
                for s in [COL_HEAD] + [hr_ci(*r[k][:3]) for r in ROWS
                                       for k in ("t90", "tst")])
    LAB_W = max(width_in(probe, f"{COHLAB[r['cohort']]}\n{ROWLAB[r['outcome']]}", TICK_PT)
                for r in ROWS)
    plt.close(probe)

    FOOT_H = len(foot_lines) * ANN_PT * 1.35 / 72.0 + EDGE + 0.10
    BOT = BOT_AXIS + FOOT_H
    H = round(BOT + AX_H + HEAD_BAND + LEG_BAND + TOP_BAND, 2)
    assert AX_X0 - GUT_X >= LAB_W + 0.10, (AX_X0 - GUT_X, LAB_W)

    COL_X_IN = W - EDGE
    AX_W = round(COL_X_IN - COL_W - COL_GAP - AX_X0, 3)
    assert AX_W >= 2.6, AX_W

    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([AX_X0 / W, BOT / H, AX_W / W, AX_H / H])

    for t in XTICKS:
        if t != 1.0:
            ax.axvline(t, color=GRID, lw=0.5, zorder=0.4)
    ax.axvline(1.0, color=INK, lw=0.9, ls=(0, (4, 3)), zorder=1)

    SUBROWS = []
    for r, yy in zip(ROWS, yD):
        for key, col, mk, off in (("tst", GREY, "s", -OFF), ("t90", BLUE, "o", OFF)):
            hr, lo, hi, q = r[key]
            sig = q < 0.05
            ax.plot([lo, hi], [yy + off] * 2, color=col, lw=CI_LW,
                    solid_capstyle="round", zorder=2)
            ax.scatter([hr], [yy + off], s=S_MARK, marker=mk, zorder=3,
                       facecolors=col if sig else "white",
                       edgecolors="white" if sig else col, linewidths=MEDGE)
            SUBROWS.append((yy + off, hr, lo, hi))
            DRAWN.append({"cohort": r["cohort"], "outcome": r["outcome"], "series": key,
                          "hr": hr, "lo": lo, "hi": hi, "q": q, "filled": bool(sig),
                          "printed": hr_ci(hr, lo, hi)})

    ax.set_yticks(yD)
    ax.set_yticklabels([f"{COHLAB[r['cohort']]}\n{ROWLAB[r['outcome']]}" for r in ROWS],
                       fontsize=TICK_PT, ha="left")
    ax.tick_params(axis="y", pad=(AX_X0 - GUT_X) * 72.0)
    bare_y(ax)
    logx(ax, XTICKS)
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xlim(*XLIM)
    ax.set_ylim(yD.min() - 0.72, yD.max() + 0.72)
    axis_covers(ax, [v for r in ROWS for k in ("t90", "tst") for v in r[k][:3]], "x", NAME)
    ax.set_xlabel("Hazard ratio per 1 SD (95% CI)", fontsize=TITLE_PT, labelpad=4)

    trCol = plt.matplotlib.transforms.blended_transform_factory(
        fig.dpi_scale_trans, ax.transData)
    for yy, hr, lo, hi in SUBROWS:
        fig.text(COL_X_IN, yy, hr_ci(hr, lo, hi), transform=trCol, fontsize=COL_FS,
                 color=INK, ha="right", va="center")
    fig.text(COL_X_IN, BOT + AX_H + 0.10, COL_HEAD, transform=fig.dpi_scale_trans,
             fontsize=COL_FS, color=INK, ha="right", va="bottom")

    fig.legend(handles=[
        Line2D([], [], color=BLUE, marker="o", ms=6.0, lw=CI_LW, markeredgecolor="white",
               markeredgewidth=MEDGE,
               label="Nocturnal oxygen, adjusted for total sleep time"),
        Line2D([], [], color=GREY, marker="s", ms=5.8, lw=CI_LW, markeredgecolor="white",
               markeredgewidth=MEDGE,
               label="Total sleep time, adjusted for nocturnal oxygen"),
        Line2D([], [], color=INK, marker="o", ms=6.0, lw=0, markerfacecolor="white",
               markeredgecolor=INK, markeredgewidth=MEDGE,
               label="Open marker, not significant after correction")],
        loc="upper left", bbox_to_anchor=(EDGE / W, (H - TOP_BAND) / H), ncol=1,
        fontsize=TICK_PT, frameon=False, handletextpad=0.55, labelspacing=0.40,
        handlelength=1.8, borderaxespad=0.0)

    fig.text(EDGE / W, (H - 0.16) / H, HEADLINE, fontsize=HEAD_PT, fontweight="bold",
             ha="left", va="top", color=INK)
    fig.text(EDGE / W, EDGE / H, FOOT_TXT, fontsize=ANN_PT, ha="left", va="bottom",
             color=INK, linespacing=1.35)

    # ---------------------------------------------------------------- gates
    nooverlap.gate(fig, NAME)
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    for t in fig.findobj(_mtext.Text):
        s = t.get_text()
        if not s or not s.strip() or not t.get_visible():
            continue
        for bad, what in (("—", "em dash"), (";", "semicolon"),
                          ("×", "multiplication sign")):
            assert bad not in s, f"{NAME}: banned {what} in {s[:50]!r}"
        for i, ch in enumerate(s):
            if ch == "–":
                assert 0 < i < len(s) - 1 and s[i - 1].isdigit() and s[i + 1].isdigit(), \
                    f"{NAME}: en dash outside a numeric range: {s[:50]!r}"
        if str(t.get_fontweight()) in ("bold", "700"):
            assert s == HEADLINE, f"{NAME}: bold outside the headline: {s[:50]!r}"
        assert float(t.get_fontsize()) >= FLOOR_PT, \
            f"{NAME}: {t.get_fontsize()} pt under the {FLOOR_PT} pt floor in {s[:50]!r}"

    bb = fig.get_tightbbox(rend)
    Wi, Hi = fig.get_size_inches()
    m = 4.0 * MM
    edges = {"left": bb.x0, "bottom": bb.y0, "right": Wi - bb.x1, "top": Hi - bb.y1}
    bad = {k: round(v / MM, 2) for k, v in edges.items() if v < m}
    assert not bad, f"{NAME}: content closer than 4 mm to an edge: {bad}"

    problems = preflight_figure(fig, journal_style("jama", columns=2))
    real = [p for p in problems if "height" not in p.lower()]
    assert not real, f"{NAME}: {real}"

    pdf = os.path.join(HERE, f"{NAME}.pdf")
    png = os.path.join(HERE, f"{NAME}.png")
    fig.savefig(pdf)
    fig.savefig(png, dpi=300)
    plt.close(fig)

json.dump({"figure": NAME,
           "source": "results_persd.csv, model S2_persd_adjT90",
           "canvas_mm": [round(W / MM, 1), round(H / MM, 1)],
           "xticks": XTICKS, "xlim": list(XLIM),
           "filled_means": "q<0.05, Benjamini-Hochberg across the 12 pairs within family",
           "counts_out_of_12": {"oxygen_fdr_significant": N_T90_FDR,
                                "total_sleep_time_fdr_significant": N_TST_FDR,
                                "total_sleep_time_raw_p_significant": N_TST_RAW},
           "rows": DRAWN},
          open(os.path.join(HERE, f"{NAME}_drawn_values.json"), "w"), indent=1)

for p in (pdf, png):
    assert os.path.exists(p) and os.path.getsize(p) > 0, p
print(f"wrote {pdf}")
print(f"wrote {png}")
print(f"{len(DRAWN)} drawn values, every one read from results_persd.csv and pinned")
print(f"canvas {W / MM:.1f} x {H / MM:.1f} mm")
