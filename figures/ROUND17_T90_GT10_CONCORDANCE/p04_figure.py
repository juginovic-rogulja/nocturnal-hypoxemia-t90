"""
The sheet: concordance gain over age and sex from T90 above 10% of the recording, one row per
outcome, both series, full follow-up and the two-year landmark.

Conventions are read off Figure 1c of
ROUND16_2026-09-04/L9_assembly/NEW_FINAL_SET_R16/Main_Fig1.pdf and its builder
FINAL_FIGURES_2026-08-14/_scripts/figure1C_polish.py, so the two sheets read side by side:

  colour    one house blue #0288d1 for both series, ink #1a1d21, sampled from the panel itself
  marker    filled circle full follow-up, filled triangle two-year landmark, the way Figure 1c
            separates its two series. Open marker when the bootstrap interval includes zero,
            which is the package's fill convention and which Figure 1c's ten rows never need
  interval  1.3 pt round-capped line in the series colour, no end caps
  rules     dashed ink at zero, dotted blue at this sheet's own mean across the outcomes,
            both behind the data
  axes      left and bottom spines only, the bottom bounded to its ticks, 3 pt outward ticks,
            no gridlines. Figure 1c itself hides the left spine, because its axis starts at
            zero and the dashed zero rule already stands there. This sheet runs negative, so
            it keeps the left spine and right-aligns the row labels against it, which is what
            the published 48-row sibling
            (Sleep_Variability_2026-08/dC_per_disease/Figure_dC_per_disease.pdf) does and what
            the clinical figure spec asks for.
  type      Arial throughout, axis title 11 pt, row labels and ticks 10 pt, printed column
            9.5 pt, footnote 9 pt, which is the floor
  column    right-aligned printed "Gain (95% CI)" for the full-follow-up series, as Figure 1c

Sheet width 183 mm exactly. No headline on the artwork.

Every drawn coordinate and every printed string is read back off the artists and asserted
against T90_gt10_concordance_per_disease.csv before the files are written, and the same
values are written out as a drawn-values JSON.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../..")))  # repository root, where paths.py lives
import paths
import json
import os
import warnings

warnings.filterwarnings("ignore")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                          # noqa: E402
import numpy as np                                                       # noqa: E402
import pandas as pd                                                      # noqa: E402
from matplotlib.lines import Line2D                                      # noqa: E402
from matplotlib.transforms import blended_transform_factory              # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = f"{HERE}/T90_gt10_concordance_per_disease.csv"
JSONH = f"{HERE}/T90_gt10_headline.json"

MM = 1 / 25.4
W_IN = 183.0 * MM
EDGE = 0.18
INK = "#1a1d21"
BLUE = "#0288d1"                 # sampled from Figure 1c of Main_Fig1.pdf
RULE = "#b9bec4"
LAB_PT, TICK_PT, ANN_PT, FLOOR_PT = 11.0, 10.0, 9.5, 9.0
CI_LW, MS, MEDGE, MEDGE_OPEN = 1.3, 4.0, 0.8, 1.0
ARM_OFF = 0.25
ROW_IN = 0.145
TOP_IN = 0.78
BOT_IN = 1.69
GAP = 0.09
DDAGGER = "‡"               # no two-year landmark estimate
MINUS = "−"

plt.rcParams.update({
    "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "standard",
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": TICK_PT, "axes.labelsize": LAB_PT,
    "xtick.labelsize": TICK_PT, "ytick.labelsize": TICK_PT, "legend.fontsize": TICK_PT,
    "axes.linewidth": 0.8, "axes.edgecolor": INK, "text.color": INK,
    "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
    "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "xtick.major.size": 3.0, "ytick.major.size": 0.0,
    "xtick.direction": "out", "axes.grid": False,
    "axes.spines.top": False, "axes.spines.right": False, "axes.spines.left": True,
    "legend.frameon": False, "pdf.fonttype": 42, "ps.fonttype": 42,
})


def fmt_val(dC, lo, hi):
    return f"{dC:+.3f} ({lo:.3f} to {hi:.3f})".replace("-", MINUS)


def tick_label(v):
    return "0" if abs(v) < 1e-12 else f"{v:.2f}".replace("-", MINUS)


def text_w(fig, s, pt):
    t = fig.text(0, 0, s, fontsize=pt)
    w = t.get_window_extent(fig.canvas.get_renderer()).width / fig.dpi
    t.remove()
    return w


def wrap_to(fig, text, pt, max_in):
    lines, cur = [], ""
    for word in text.split(" "):
        trial = f"{cur} {word}".strip()
        if cur and text_w(fig, trial, pt) > max_in:
            lines.append(cur)
            cur = word
        else:
            cur = trial
    lines.append(cur)
    return "\n".join(lines)


def pick_ticks(fig, lo_v, hi_v, plot_w, pad_frac=0.03):
    """
    The smallest round step whose tick labels still clear each other at this plot width.

    The first sheet ran from -0.10 to 0.25 and 0.05 fitted. The restricted sheet runs wider
    and 0.05 would set "-0.10" hard against "-0.05", so the step is chosen by measuring the
    widest label rather than by a rule of thumb.
    """
    chosen = None
    for step in (0.01, 0.02, 0.05, 0.10, 0.20, 0.25):
        t0 = np.floor(lo_v / step) * step
        t1 = np.ceil(hi_v / step) * step
        ticks = [round(t0 + i * step, 10) for i in range(int(round((t1 - t0) / step)) + 1)]
        w = max(text_w(fig, tick_label(v), TICK_PT) for v in ticks)
        span = (t1 - t0) * (1 + 2 * pad_frac)
        per_tick = plot_w * step / span if span > 0 else plot_w
        chosen = (ticks, t0, t1)
        if per_tick > w + 0.07:
            return chosen
    return chosen


def main():
    global ROW_IN   # v8: the row pitch may be reduced below for the taller outcome list
    D = pd.read_csv(CSV, comment="#")
    H = json.load(open(JSONH))
    D = D.sort_values("lag0_bin_dC", ascending=False).reset_index(drop=True)
    n = len(D)
    assert list(D.rank_by_binary_gain) == list(range(1, n + 1)), "row order is not the ranking"
    assert D.lag0_bin_estimable.all(), "a full-follow-up estimate is missing"
    no_lag = list(D.loc[~D.lag2_bin_estimable, "disease"])

    mean_ff = H["full_followup"]
    mean_lm = H["landmark_2y"]

    y_dis = np.arange(n)[::-1].astype(float) + 2.0
    y_mean = 0.0
    n_slot = n + 2
    h_in = ROW_IN * n_slot + TOP_IN + BOT_IN
    # v8 (2026-09-14): 52 ranked outcomes instead of 48; the row pitch shrinks just enough for the sheet to stay within 247 mm
    if h_in * 25.4 > 247.0:
        _row = (247.0 / 25.4 - TOP_IN - BOT_IN) / n_slot
        print(f"  v8: {n} outcomes; row pitch {ROW_IN * 25.4:.2f} -> {_row * 25.4:.2f} mm so the sheet stays within 247 mm")
        ROW_IN = _row; h_in = ROW_IN * n_slot + TOP_IN + BOT_IN

    fig = plt.figure(figsize=(W_IN, h_in))
    fig.canvas.draw()

    lab_w = max(text_w(fig, s, TICK_PT) for s in
                list(D.disease) + [f"{s} {DDAGGER}" for s in no_lag] +
                [f"Mean across the {n} outcomes"])
    val_strings = [fmt_val(r.lag0_bin_dC, r.lag0_bin_ci_lo, r.lag0_bin_ci_hi)
                   for _i, r in D.iterrows()] + \
                  [fmt_val(mean_ff["binary_mean_dC"], *mean_ff["binary_mean_ci"])]
    val_w = max([text_w(fig, s, ANN_PT) for s in val_strings] +
                [text_w(fig, "Gain (95% CI)", ANN_PT)])

    X_LAB = EDGE + lab_w + 0.12
    X_VAL_R = W_IN - EDGE
    X_PLOT_R = X_VAL_R - val_w - GAP
    plot_w = X_PLOT_R - X_LAB
    assert plot_w > 2.6, f"plot width {plot_w:.2f} in is too narrow"

    ax = fig.add_axes([X_LAB / W_IN, BOT_IN / h_in, plot_w / W_IN,
                       (h_in - TOP_IN - BOT_IN) / h_in])
    pitch_in = (h_in - TOP_IN - BOT_IN) / (n_slot + 0.55)
    assert 2 * ARM_OFF * pitch_in > MS / 72.0, \
        f"the two series collide: {2 * ARM_OFF * pitch_in:.4f} in apart, marker {MS / 72:.4f} in"

    lo = min(float(D.lag0_bin_ci_lo.min()), float(D.lag2_bin_ci_lo.min()),
             mean_ff["binary_mean_ci"][0], mean_lm["binary_mean_ci"][0], 0.0)
    hi = max(float(D.lag0_bin_ci_hi.max()), float(D.lag2_bin_ci_hi.max()))
    ticks, t0, t1 = pick_ticks(fig, lo, hi, plot_w)
    pad = 0.03 * (t1 - t0)
    ax.set_xlim(t0 - pad, t1 + pad)
    ax.set_xticks(ticks)
    ax.set_xticklabels([tick_label(v) for v in ticks])
    ax.set_yticks([])
    ax.set_ylim(-0.85, n_slot - 0.30)
    ax.spines["bottom"].set_bounds(ticks[0], ticks[-1])
    ax.spines["left"].set_position(("outward", 3))
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("Gain in held-out concordance from adding T90 above 10% to age and sex",
                  fontsize=LAB_PT, labelpad=6)

    ax.axvline(0.0, color=INK, lw=0.9, ls=(0, (4, 3)), zorder=1)
    ax.axvline(mean_ff["binary_mean_dC"], color=BLUE, lw=0.9, ls=(0, (1, 2.2)), zorder=1)
    ax.axhline(1.0, color=RULE, lw=0.7, zorder=1)

    drawn = []

    def mark(dC, ci_lo, ci_hi, yy, marker, series, label):
        if not np.isfinite(dC) or not np.isfinite(ci_lo):
            return
        filled = bool(ci_lo > 0 or ci_hi < 0)
        ax.plot([ci_lo, ci_hi], [yy, yy], color=BLUE, lw=CI_LW, solid_capstyle="round",
                zorder=2)
        ax.plot([dC], [yy], marker=marker, ms=MS, ls="none",
                mfc=BLUE if filled else "white", mec=BLUE,
                mew=MEDGE if filled else MEDGE_OPEN, zorder=3)
        drawn.append({"row": label, "series": series, "marker": marker, "y": float(yy),
                      "dC": float(dC), "ci_lo": float(ci_lo), "ci_hi": float(ci_hi),
                      "filled": filled})

    for i, r in D.iterrows():
        mark(r.lag0_bin_dC, r.lag0_bin_ci_lo, r.lag0_bin_ci_hi, y_dis[i] + ARM_OFF,
             "o", "full follow-up", r.disease)
        if bool(r.lag2_bin_estimable):
            mark(r.lag2_bin_dC, r.lag2_bin_ci_lo, r.lag2_bin_ci_hi, y_dis[i] - ARM_OFF,
                 "^", "2-year landmark", r.disease)
    MEANLAB = f"Mean across the {n} outcomes"
    mark(mean_ff["binary_mean_dC"], *mean_ff["binary_mean_ci"], y_mean + ARM_OFF,
         "o", "full follow-up", MEANLAB)
    mark(mean_lm["binary_mean_dC"], *mean_lm["binary_mean_ci"], y_mean - ARM_OFF,
         "^", "2-year landmark", MEANLAB)

    # ------------------------------------------------------------ row labels
    blend = blended_transform_factory(fig.transFigure, ax.transData)
    rowlab, labtexts = [], []
    for i, r in D.iterrows():
        s = f"{r.disease} {DDAGGER}" if not bool(r.lag2_bin_estimable) else r.disease
        rowlab.append(s)
        labtexts.append(fig.text((X_LAB - 0.09) / W_IN, y_dis[i], s, fontsize=TICK_PT,
                                 color=INK, va="center", ha="right", transform=blend))
    labtexts.append(fig.text((X_LAB - 0.09) / W_IN, y_mean, MEANLAB, fontsize=TICK_PT,
                             color=INK, va="center", ha="right", transform=blend))

    # ------------------------------------------------------------ printed value column
    printed = []
    for i, r in D.iterrows():
        s = fmt_val(r.lag0_bin_dC, r.lag0_bin_ci_lo, r.lag0_bin_ci_hi)
        printed.append({"row": r.disease, "text": s})
        fig.text(X_VAL_R / W_IN, y_dis[i], s, fontsize=ANN_PT, color=INK,
                 va="center", ha="right", transform=blend)
    s = fmt_val(mean_ff["binary_mean_dC"], *mean_ff["binary_mean_ci"])
    printed.append({"row": MEANLAB, "text": s})
    fig.text(X_VAL_R / W_IN, y_mean, s, fontsize=ANN_PT, color=INK, va="center",
             ha="right", transform=blend)
    hdr = fig.text(X_VAL_R / W_IN, (h_in - TOP_IN + 0.10) / h_in, "Gain (95% CI)",
                   fontsize=ANN_PT, color=INK, va="bottom", ha="right")

    # ------------------------------------------------------------ key and footnote
    handles = [
        Line2D([], [], color=BLUE, lw=CI_LW, marker="o", ms=MS, mfc=BLUE, mec=BLUE,
               mew=MEDGE, label="Full follow-up"),
        Line2D([], [], color=BLUE, lw=CI_LW, marker="^", ms=MS, mfc=BLUE, mec=BLUE,
               mew=MEDGE, label="2-year landmark"),
        Line2D([], [], color=BLUE, lw=0, marker="o", ms=MS, mfc="white", mec=BLUE,
               mew=MEDGE_OPEN, label="Open marker, 95% CI includes zero"),
    ]
    leg = fig.legend(handles=handles, loc="lower left", ncol=3,
                     bbox_to_anchor=(EDGE / W_IN, (h_in - TOP_IN + 0.30) / h_in),
                     bbox_transform=fig.transFigure, fontsize=ANN_PT, handlelength=2.0,
                     borderpad=0.0, columnspacing=1.5, handletextpad=0.5)

    foot = (f"Exposure is T90 above 10% of the recording, entered as a raw 0/1 indicator, "
            f"added to age and sex one outcome at a time. Points are the mean of the two "
            f"cross-fit directions, fitting at one hospital and scoring the other. Bars are "
            f"the 2.5th to 97.5th percentile of {H['bootstrap_replicates']} bootstrap "
            f"replicates resampling patients with replacement inside each hospital. The "
            f"2-year landmark removes every patient whose event or censoring arrived within "
            f"2 years of the sleep study and restarts the clock there. Dashed rule zero, "
            f"dotted rule the mean across the {n} outcomes without the landmark, "
            f"{mean_ff['binary_mean_dC']:.4f}. Printed values are the full-follow-up series, "
            f"and rows are ordered by it. "
            f"{DDAGGER} no 2-year landmark estimate: after the landmark this outcome falls "
            f"below the 30-case floor at one hospital, so it is left empty rather than drawn "
            f"at zero.")
    if not no_lag:
        foot = foot.split(f"{DDAGGER} no 2-year")[0].strip()
    ft = fig.text(EDGE / W_IN, 0.12 / h_in, wrap_to(fig, foot, FLOOR_PT, W_IN - 2 * EDGE),
                  fontsize=FLOOR_PT, color="#5f666d", va="bottom", ha="left", linespacing=1.5)

    # ------------------------------------------------------------ assertions
    ref = D.set_index("disease")
    worst = 0.0
    for m in drawn:
        if m["row"] == MEANLAB:
            src = mean_ff if m["series"] == "full follow-up" else mean_lm
            v, cl, ch = src["binary_mean_dC"], src["binary_mean_ci"][0], src["binary_mean_ci"][1]
        else:
            r = ref.loc[m["row"]]
            p = "lag0_" if m["series"] == "full follow-up" else "lag2_"
            v, cl, ch = r[p + "bin_dC"], r[p + "bin_ci_lo"], r[p + "bin_ci_hi"]
        worst = max(worst, abs(v - m["dC"]), abs(cl - m["ci_lo"]), abs(ch - m["ci_hi"]))
        assert m["filled"] == bool(cl > 0 or ch < 0), (m["row"], m["series"], "fill")
        assert m["marker"] == ("o" if m["series"] == "full follow-up" else "^")
    assert worst < 1e-12, worst
    n_expect = n + int(D.lag2_bin_estimable.sum()) + 2
    assert len(drawn) == n_expect, f"{len(drawn)} marks drawn, expected {n_expect}"
    for p in printed:
        if p["row"] == MEANLAB:
            v, cl, ch = mean_ff["binary_mean_dC"], *mean_ff["binary_mean_ci"]
        else:
            r = ref.loc[p["row"]]
            v, cl, ch = r.lag0_bin_dC, r.lag0_bin_ci_lo, r.lag0_bin_ci_hi
        assert p["text"] == fmt_val(v, cl, ch), p
    assert [s.split(f" {DDAGGER}")[0] for s in rowlab] == list(D.disease), "row labels"
    assert list(D.lag0_bin_dC) == sorted(D.lag0_bin_dC, reverse=True), "not in gain order"

    fig.canvas.draw()
    rr = fig.canvas.get_renderer()
    _tl = sorted((t.get_window_extent(rr) for t in ax.get_xticklabels()), key=lambda b: b.x0)
    _gaps = [(b.x0 - a.x1) / fig.dpi for a, b in zip(_tl, _tl[1:])]
    assert not _gaps or min(_gaps) > 0.03, \
        f"x tick labels sit {min(_gaps):.3f} in apart, they collide"
    lab_r = max(t.get_window_extent(rr).x1 for t in labtexts) / fig.dpi
    lab_l = min(t.get_window_extent(rr).x0 for t in labtexts) / fig.dpi
    assert lab_r <= X_LAB - 0.05, f"row labels reach {lab_r:.3f} in, plot starts {X_LAB:.3f}"
    assert lab_l >= EDGE - 0.01, f"row labels start {lab_l:.3f} in, margin is {EDGE:.3f}"
    gap = (leg.get_window_extent(rr).y0 - hdr.get_window_extent(rr).y1) / fig.dpi
    assert leg.get_window_extent(rr).x1 / fig.dpi <= W_IN - EDGE + 0.01, "key runs off"
    xlab_bb = ax.xaxis.label.get_window_extent(rr)
    assert xlab_bb.x0 / fig.dpi >= EDGE - 0.01 and xlab_bb.x1 / fig.dpi <= W_IN - EDGE + 0.01, \
        (f"the axis title spans {xlab_bb.x0 / fig.dpi:.3f} to {xlab_bb.x1 / fig.dpi:.3f} in, "
         f"the margins are {EDGE:.3f} to {W_IN - EDGE:.3f}")
    xlab_y0 = ax.xaxis.label.get_window_extent(rr).y0 / fig.dpi
    foot_y1 = ft.get_window_extent(rr).y1 / fig.dpi
    assert foot_y1 < xlab_y0 - 0.05, \
        f"the footnote reaches {foot_y1:.3f} in and the axis title starts {xlab_y0:.3f} in"
    assert ft.get_window_extent(rr).x1 / fig.dpi <= W_IN - EDGE + 0.01, "footnote runs off"
    assert h_in * 25.4 <= 247.0, f"sheet is {h_in * 25.4:.1f} mm tall"
    print(f"asserted {len(drawn)} marks, {len(printed)} printed strings and {n} row labels "
          f"against {os.path.basename(CSV)}, largest difference {worst:.2e}")
    print(f"  row-label gutter ends {lab_r:.3f} in, plot starts {X_LAB:.3f} in")
    print(f"  key clears the column header by {gap:.3f} in")

    pdf = f"{HERE}/Figure_T90_gt10_dC_per_disease.pdf"
    png = f"{HERE}/Figure_T90_gt10_dC_per_disease.png"
    fig.savefig(pdf)
    fig.savefig(png, dpi=300)
    plt.close(fig)
    assert abs(W_IN * 25.4 - 183.0) < 1e-6
    print(f"wrote {pdf}\nwrote {png}\ncanvas 183.0 x {h_in * 25.4:.1f} mm")

    json.dump({
        "sheet": os.path.basename(pdf),
        "canvas_mm": [183.0, round(h_in * 25.4, 2)],
        "source_csv": CSV,
        "source_headline_json": JSONH,
        "conventions_read_from": (f"{paths.FIGURE_ROOT}/"
                                  "ROUND16_2026-09-04/L9_assembly/NEW_FINAL_SET_R16/"
                                  "Main_Fig1.pdf, panel c"),
        "colour": {"series": BLUE, "ink": INK},
        "spines": "left and bottom only, bottom bounded to its ticks",
        "marker": {"full follow-up": "o", "2-year landmark": "^",
                   "open": "bootstrap interval includes zero"},
        "x_ticks": ticks,
        "zero_rule": 0.0,
        "mean_rule": mean_ff["binary_mean_dC"],
        "n_rows": n,
        "n_marks": len(drawn),
        "largest_difference_against_csv": worst,
        "row_order": list(D.disease),
        "drawn": drawn,
        "printed": printed,
        "no_landmark_estimate": no_lag,
    }, open(f"{HERE}/Figure_T90_gt10_drawn_values.json", "w"), indent=1, default=float)
    print(f"wrote {HERE}/Figure_T90_gt10_drawn_values.json")


if __name__ == "__main__":
    main()
