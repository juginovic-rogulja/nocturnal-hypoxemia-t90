"""
The figure: gain in held-out discrimination from T90, one row per disease, sorted descending,
with the bootstrap interval, and the 2-year landmark arm beside it.

Design is the Nature Medicine layer used in T90_Manuscript/NATURE_MED_VERSION/_scripts, held
here as literal constants rather than by importing that tree, so nothing can write into the
frozen figure folders.

  canvas   171 mm wide exactly, height set by the rows
  type     Arial, axis title 11 pt, ticks and row labels 10 pt, printed columns 9.5 pt,
           floor 9 pt
  colour   ink #1a1d21, house primary #1f4257 for the paper's windowing, comparison grey
           #8a9099 for the 2-year landmark
  axes     left and bottom spines only at 0.8 pt, 3 pt outward ticks
  marks    interval lines 1.4 pt round caps, markers with a 0.8 pt white edge, the zero rule
           0.9 pt dashed ink and the paper's average 0.9 pt dotted grey, both behind the data

Two channels, deliberately kept apart. Colour is the arm, blue without the landmark and grey
with it. Fill is the verdict on zero, solid when the interval clears it and open when it does
not. Both survive a greyscale print, which eleven organ hues would not, so organ system is a
printed column instead.

Every drawn coordinate is read back off the artists and asserted against dC_per_disease.csv.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import os
import warnings

warnings.filterwarnings("ignore")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                      # noqa: E402
import numpy as np                                                   # noqa: E402
import pandas as pd                                                  # noqa: E402
from matplotlib.lines import Line2D                                  # noqa: E402
import sys as _sys_v8; _sys_v8.path.insert(0, paths.ANALYSIS_DIR)  # v8 sweep 2026-09-12
from cohort_spec import N_RANKED_OUTCOMES  # v8 sweep 2026-09-12
__PUB_DC_V7__ = float(__import__("pandas").read_csv(f"{paths.NUMBERS_DIR}/ranking_v3.csv", comment="#").set_index("feature").loc["spo2_pct_below_90", "dC"])   # v7: T90 gain from the CURRENT ranking, was the August literal 0.02033875497249791


HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "dC_per_disease.csv")
PUB_DC = __PUB_DC_V7__          # the paper's average across the 48, for the reference rule

MM = 1 / 25.4
W_IN = 171.0 * MM
INK = "#1a1d21"
BLUE = "#1f4257"
GREY = "#8a9099"
RULE = "#b9bec4"
TITLE_PT, TICK_PT, ANN_PT, FLOOR_PT = 11.0, 10.0, 9.5, 9.0
CI_LW, MEDGE = 1.4, 0.8
MS = 4.6                              # marker diameter in points
OFF = 0.22                            # vertical offset of the two arms, in rows

# the width budget, in inches, left to right across the 171 mm canvas
X_LAB = 1.87          # disease names, widest 1.606 in at 10 pt, plus a 4.6 mm margin
X_PLOT = 3.73         # the data
X_GAP = 0.07
X_ORGAN = 0.88        # organ system column, widest 0.762 in at 9 pt
X_RIGHT = 0.18        # no ink inside 4.6 mm of the right edge
ROW_IN = 0.175        # 4.4 mm per row
TOP_IN = 0.40
BOT_IN = 2.40
DAGGER = "†"     # marks a negative control in the row label

plt.rcParams.update({
    "figure.dpi": 300, "savefig.dpi": 300,
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": TICK_PT, "axes.labelsize": TITLE_PT, "axes.titlesize": TITLE_PT,
    "xtick.labelsize": TICK_PT, "ytick.labelsize": TICK_PT, "legend.fontsize": ANN_PT,
    "axes.linewidth": 0.8, "axes.edgecolor": INK, "text.color": INK,
    "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
    "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "xtick.major.size": 3.0, "ytick.major.size": 0.0,
    "xtick.direction": "out", "ytick.direction": "out",
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "pdf.fonttype": 42, "ps.fonttype": 42,
})


def wrap_to(fig, text, pt, max_in):
    """Greedy wrap on measured Arial widths, so a caption can never run off the canvas."""
    r = fig.canvas.get_renderer()

    def w(s):
        t = fig.text(0, 0, s, fontsize=pt)
        bb = t.get_window_extent(r)
        t.remove()
        return bb.width / fig.dpi

    lines, cur = [], ""
    for word in text.split(" "):
        trial = f"{cur} {word}".strip()
        if cur and w(trial) > max_in:
            lines.append(cur)
            cur = word
        else:
            cur = trial
    lines.append(cur)
    return "\n".join(lines)


def main():
    raw = pd.read_csv(CSV)
    dis = raw[raw.row_type == "disease"]
    S = dis[dis.arm == "standard"].sort_values("dC", ascending=False).reset_index(drop=True)
    has_lag = (dis.arm == "lag2").any()
    L = dis[dis.arm == "lag2"].set_index("outcome") if has_lag else None
    n = len(S)

    mean_rows = raw[(raw.row_type == "mean_of_outcomes") & (raw.outcome == "__mean_all__")]
    mS = mean_rows[mean_rows.arm == "standard"].iloc[0]
    mL = mean_rows[mean_rows.arm == "lag2"].iloc[0] if has_lag else None

    # row 0 is the mean, rows 2..n+1 are the diseases, so a blank slot separates them
    y_dis = np.arange(n)[::-1].astype(float) + 2.0
    y_mean = 0.0
    n_slot = n + 2

    h_in = ROW_IN * n_slot + TOP_IN + BOT_IN
    # v8 (2026-09-13): the label column is sized to the widest disease name at draw time (the v8 outcome list carries longer
    # names, e.g. the added arrhythmia outcome); the plot column gives up the same width so the 171 mm canvas is unchanged.
    global X_LAB, X_PLOT
    _rowlab_pre = [f"{d} {DAGGER}" if c == "yes" else d for d, c in zip(S.disease, S.negative_control.fillna(""))]
    _tmp = plt.figure(figsize=(W_IN, 1.0)); _tmp.canvas.draw(); _r = _tmp.canvas.get_renderer()
    _widest = max(_tmp.text(0, 0, lab, fontsize=TICK_PT).get_window_extent(_r).width for lab in _rowlab_pre + ["Mean across outcomes"]) / _tmp.dpi
    plt.close(_tmp)
    _need = _widest + 4.6 * MM + 0.12   # widest label plus the 4.6 mm margin and the outward tick
    if _need > X_LAB:
        print(f"  label column widened {X_LAB:.3f} -> {_need:.3f} in for the widest name ({_widest:.3f} in); plot column narrowed by the same")
        X_PLOT = X_PLOT - (_need - X_LAB); X_LAB = _need
    fig = plt.figure(figsize=(W_IN, h_in))
    ax = fig.add_axes([X_LAB / W_IN, BOT_IN / h_in, X_PLOT / W_IN,
                       (h_in - TOP_IN - BOT_IN) / h_in])

    ax.axvline(0.0, color=INK, lw=0.9, ls=(0, (4, 3)), zorder=1)
    ax.axvline(PUB_DC, color=GREY, lw=0.9, ls=(0, (1, 2.2)), zorder=1)

    drawn = []

    def series(vals, y, colour):
        """One arm: interval line then marker, solid when the interval clears zero."""
        for yy, r in zip(y, vals):
            if r is None or not np.isfinite(r["dC"]) or not np.isfinite(r["ci_lo"]):
                continue
            sig = bool(r["excludes_zero"])
            ax.plot([r["ci_lo"], r["ci_hi"]], [yy, yy], color=colour, lw=CI_LW,
                    solid_capstyle="round", zorder=2)
            ax.plot([r["dC"]], [yy], marker="o", ms=MS, mfc=colour if sig else "white",
                    mec=colour if sig else colour, mew=MEDGE if sig else 1.1,
                    color="none", zorder=3)
            drawn.append({"disease": r["disease"], "arm": r["arm"], "y": float(yy),
                          "dC": float(r["dC"]), "ci_lo": float(r["ci_lo"]),
                          "ci_hi": float(r["ci_hi"]), "filled": sig})

    series([S.iloc[i].to_dict() for i in range(n)], y_dis + OFF, BLUE)
    if has_lag:
        series([L.loc[o].to_dict() if o in L.index else None for o in S.outcome],
               y_dis - OFF, GREY)
    series([mS.to_dict()], np.array([y_mean + OFF]), BLUE)
    if has_lag:
        series([mL.to_dict()], np.array([y_mean - OFF]), GREY)

    rowlab = [f"{s} {DAGGER}" if c == "yes" else s
              for s, c in zip(S.disease, S.negative_control.fillna(""))]
    ax.set_yticks(list(y_dis) + [y_mean])
    ax.set_yticklabels(rowlab + ["Mean across outcomes"], fontsize=TICK_PT)
    ax.set_ylim(-0.85, n_slot - 0.30)

    lo = float(min(dis.ci_lo.min(), 0.0))
    hi = float(dis.ci_hi.max())
    pad = 0.045 * (hi - lo)
    ax.set_xlim(lo - pad, hi + pad)
    ax.set_xticks([-0.05, 0.0, 0.05, 0.10, 0.15, 0.20, 0.25])
    ax.set_xticklabels(["−0.05", "0", "0.05", "0.10", "0.15", "0.20", "0.25"])
    ax.set_xlabel("Gain in held-out concordance from adding T90 to age and sex", labelpad=5)
    ax.tick_params(axis="y", length=0, pad=3)
    ax.spines["left"].set_position(("outward", 3))
    ax.spines["bottom"].set_bounds(-0.05, 0.25)
    ax.axhline(1.0, color=RULE, lw=0.7, zorder=1)          # the mean sits below this rule

    # ------------------------------------------------------------ printed column
    fx_org = (X_LAB + X_PLOT + X_GAP) / W_IN
    blend = matplotlib.transforms.blended_transform_factory(fig.transFigure, ax.transData)
    for i in range(n):
        fig.text(fx_org, y_dis[i], S.organ_group[i], fontsize=FLOOR_PT, color=GREY,
                 va="center", ha="left", transform=blend)
    fig.text(fx_org, (h_in - TOP_IN + 0.05) / h_in, "Organ system", fontsize=FLOOR_PT,
             color=INK, va="bottom", ha="left", style="italic")

    # ------------------------------------------------------------ legend and footnote
    handles = [
        Line2D([], [], color=BLUE, lw=CI_LW, marker="o", ms=MS, mec="white", mew=MEDGE,
               label="No landmark"),
        Line2D([], [], color=GREY, lw=CI_LW, marker="o", ms=MS, mec="white", mew=MEDGE,
               label="2-year landmark"),
        Line2D([], [], color=BLUE, lw=CI_LW, marker="o", ms=MS, mfc="white", mec=BLUE,
               mew=1.1, label="open marker, interval includes zero"),
    ]
    leg = fig.legend(handles=handles, loc="lower left", ncol=3,
                     bbox_to_anchor=(X_RIGHT / W_IN, 1.66 / h_in),
                     bbox_transform=fig.transFigure, fontsize=ANN_PT, handlelength=2.0,
                     borderpad=0.0, columnspacing=1.6)

    nrep = int(S.n_replicates.max())
    miss = [r.disease for _i, r in S.iterrows()
            if has_lag and r.outcome in L.index and not np.isfinite(L.loc[r.outcome].dC)]
    foot = ("Points are the mean of the two cross-fit directions, fitting at one hospital and "
            f"scoring the other. Bars are the 2.5th to 97.5th percentile of {nrep} bootstrap "
            "replicates resampling patients with replacement inside each hospital. The 2-year "
            "landmark removes every patient whose event or censoring arrived within 2 years of "
            "the sleep study and restarts the clock there, which is the rule the paper's "
            "published lag ladder uses. Dashed rule zero, dotted rule the paper's average "
            f"across all {N_RANKED_OUTCOMES} diseases, {PUB_DC:.4f}. {DAGGER} negative control. Only the row "
            "below the rule averages anything across diseases, over the "
            f"{int(mS.n_outcomes_in_mean)} outcomes estimable without the landmark and the "
            f"{int(mL.n_outcomes_in_mean) if mL is not None else 0} estimable with it, and "
            "the mean is taken inside each bootstrap replicate so the correlation between "
            "diseases is carried.")
    if miss:
        foot += (" " + ", ".join(miss) + " has fewer than 30 incident cases at one hospital "
                 "after the landmark, which is the ranking's own floor, so it carries no "
                 "2-year estimate.")
    ft = fig.text(X_RIGHT / W_IN, 0.10 / h_in, wrap_to(fig, foot, FLOOR_PT, W_IN - 2 * X_RIGHT),
                  fontsize=FLOOR_PT, color=GREY, va="bottom", ha="left", linespacing=1.5)

    # ------------------------------------------------------------ assert the drawn values
    dd = pd.DataFrame(drawn)
    key = raw.set_index(["outcome", "arm"]) if False else None
    ref = pd.concat([dis, raw[raw.outcome == "__mean_all__"]])
    ref = ref.set_index(["disease", "arm"])
    worst = 0.0
    for _i, r in dd.iterrows():
        w = ref.loc[(r.disease, r.arm)]
        worst = max(worst, abs(w.dC - r.dC), abs(w.ci_lo - r.ci_lo), abs(w.ci_hi - r.ci_hi))
        assert bool(w.excludes_zero) == bool(r.filled), (r.disease, r.arm)
    assert worst < 1e-12, worst
    n_expect = int(np.isfinite(S.dC).sum()) + 1
    if has_lag:
        n_expect += int(np.isfinite(dis[dis.arm == "lag2"].dC).sum()) + 1
    assert len(dd) == n_expect, f"{len(dd)} marks drawn, expected {n_expect}"
    assert [t.get_text() for t in ax.get_yticklabels()][:n] == rowlab, "row order"
    assert [s.split(f" {DAGGER}")[0] for s in rowlab] == list(S.disease), "row labels"
    org = [t.get_text() for t in fig.texts if t.get_text() in set(S.organ_group)]
    assert len(org) == n, f"{len(org)} organ strings, expected {n}"
    print(f"asserted {len(dd)} drawn marks and intervals, {n} row labels and {n} organ labels "
          f"against {os.path.basename(CSV)}, largest difference {worst:.2e}")

    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    gap = (leg.get_window_extent(r).y0 - ft.get_window_extent(r).y1) / fig.dpi
    assert gap > 0.04, f"legend sits {gap:.3f} in above the footnote, they collide"
    print(f"  legend clears the footnote by {gap:.3f} in")
    for what, arts, x0, x1 in (
            ("row labels", ax.get_yticklabels(), 0.0, X_LAB),
            ("organ", [t for t in fig.texts if t.get_text() in set(S.organ_group)],
             X_LAB + X_PLOT + X_GAP, W_IN - X_RIGHT),
            ("footnote", [ft], 0.0, W_IN)):
        l = min(a.get_window_extent(r).x0 for a in arts) / fig.dpi
        rr = max(a.get_window_extent(r).x1 for a in arts) / fig.dpi
        assert l >= x0 - 0.01 and rr <= x1 + 0.01, \
            f"{what} spans {l:.3f} to {rr:.3f} in, budget {x0:.3f} to {x1:.3f}"
        print(f"  {what:<12} {l:.3f} to {rr:.3f} in, budget {x0:.3f} to {x1:.3f}")

    pdf = os.path.join(HERE, "Figure_dC_per_disease.pdf")
    png = os.path.join(HERE, "Figure_dC_per_disease.png")
    fig.savefig(pdf)
    fig.savefig(png, dpi=300)
    plt.close(fig)
    assert abs(W_IN * 25.4 - 171.0) < 1e-6
    print(f"wrote {pdf}\nwrote {png}\ncanvas 171.0 x {h_in * 25.4:.1f} mm")
    dd.to_csv(os.path.join(HERE, "_work", "figure_drawn_values.csv"), index=False)


if __name__ == "__main__":
    main()
