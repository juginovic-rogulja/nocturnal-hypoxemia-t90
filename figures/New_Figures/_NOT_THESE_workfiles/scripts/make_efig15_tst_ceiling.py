"""
eFigure15  the reliability ceiling on sleep duration, and where it actually lands.

ALEN'S FRAMING, 2026-08-08, which replaces the earlier one. The argument is not about repeat
measurement. It is that people sleep differently at home than in the laboratory. So: we knew one
laboratory night may not represent habitual sleep, we corrected for that, and even under perfect
measurement total sleep time reaches only rank 84 of 197, and never better than rank 31 across the
ten scenarios of the published ladder.

THE OPEN PROBLEM THIS SHEET HAS TO SETTLE
-----------------------------------------
An eleventh scenario reaches rank 2, which would put sleep duration inside the top 10 and
contradict the Discussion clause. It is `classical_attenuation_lower_95_any_window`, and it stacks
two things:

  1. the CLASSICAL attenuation law, exponent 1.0, and
  2. an ICC of 0.171, which is not an estimate of reliability but the lower 95% bound of the single
     weakest of eight overlapping gap windows, from 114 pairs two years or more apart.

Neither survives contact with the data.

  On (1), the attenuation exponent was MEASURED in this cohort rather than assumed, on the three
  measures that have repeat nights: apnea-hypopnea index 0.536, oxygen nadir 0.546, T90 0.549, mean
  0.544. That is the square-root law. The classical law is not merely generous, it is the wrong
  law, and panel C shows the measurement.

  On (2), the point estimate for that window is 0.434 and the all-pairs estimate is 0.482 with a
  95% interval of [0.381, 0.567]. 0.171 sits far outside that interval, and it is reached only by
  taking the minimum over eight overlapping windows and then the lower bound of that minimum.

So the sheet draws all seventeen fitted scenarios, marks the one that reaches rank 2, and states
the defensible range: under every attenuation law this cohort's own repeat nights support, total
sleep time at perfect measurement lands no better than rank 24 of 197. The claim that it never
enters the top 10 holds. The clause as written, "never better than 31", does not: the measured-law
family reaches rank 24, and 24 is the number the Discussion should carry.

Sources
  Sleep_Variability_2026-08/tst_ceiling/attack_summary.json      the scenario ladders A1, A1b, A1c
  Sleep_Variability_2026-08/tst_ceiling/attack_results.csv       the same values, flat
  Sleep_Variability_2026-08/tst_ceiling/summary.json             the effect-to-gain map
  Sleep_Variability_2026-08/tst_ceiling/effect_vs_gain.csv       E_meanabs of every measurement
  Sleep_Variability_2026-08/tst_ceiling/ranking_reproduced.csv   the 197-measurement ranking
  Sleep_Variability_2026-08/tst_ceiling/_iccs.json               the repeat-night ICCs
  Sleep_Variability_2026-08/habitual_sleep/habitual_vs_lab_tst_correlation.csv
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../..")))  # repository root, where paths.py lives
import paths
import json
import sys
import textwrap

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import sys as _sys_v8; _sys_v8.path.insert(0, paths.ANALYSIS_DIR)  # v8 sweep 2026-09-12
from cohort_spec import N_MEASURES  # v8 sweep 2026-09-12

sys.path.insert(0, f"{paths.FIGURE_ROOT}/New_Figures/"
                   "_NOT_THESE_workfiles/scripts")
from splitstyle import (rc, save, panel_label, bare_y, axis_covers, luma, BLUE, BLUE_MID,
                        GREY, GREY_PALE, GREY_MID, INK, RULE, TITLE, LABEL, TICK, ANNOT,
                        SMALL)

TC = f"{paths.SV_ROOT}/tst_ceiling"
HAB = f"{paths.SV_ROOT}/habitual_sleep"
WORK = f"{paths.FIGURE_ROOT}/New_Figures/_NOT_THESE_workfiles"
NUM = paths.NUMBERS_DIR

A = json.load(open(f"{TC}/attack_summary.json"))
S = json.load(open(f"{TC}/summary.json"))
ICC = json.load(open(f"{TC}/_iccs.json"))
RK = pd.read_csv(f"{TC}/ranking_reproduced.csv")
EVG = pd.read_csv(f"{TC}/effect_vs_gain.csv")
COR = pd.read_csv(f"{HAB}/habitual_vs_lab_tst_correlation.csv")

N_MEAS = len(RK)
SLOPE = float(S["map_primary"]["slope"])
INTERCEPT = float(S["map_primary"]["intercept"])
E_NOW = float(EVG[EVG.feature == "TST_min"].E_meanabs.iloc[0])
TARGET10 = float(S["reproduced"]["dC_rank10"])
ICC_LAB = float(ICC["TST"]["icc_raw"])
N_PAIRS = int(ICC["TST"]["n_pairs"])
EXP_MEASURED = float(A["A1c"]["exponent_used"])
assert N_MEAS == N_MEASURES, (N_MEAS, ICC_LAB)  # v8 sweep: the ICC is read from _iccs.json, no longer compared to a typed literal

PROV = []


def rank_at(icc, power):
    """
    Where total sleep time would rank if its measurement error were removed entirely.

    Identical arithmetic to tst_ceiling/attack.py: the observed mean absolute log hazard is
    disattenuated by (1/ICC)**power, mapped to a concordance gain by the quadratic effect-to-gain
    map fitted across all 197 measurements, and looked up in the reproduced ranking.
    """
    e = E_NOW * (1.0 / np.asarray(icc, float)) ** power
    dc = SLOPE * e ** 2 + INTERCEPT
    if np.ndim(dc) == 0:
        return float(dc), int((RK.dC > dc).sum() + 1)
    return dc, np.array([(RK.dC > d).sum() + 1 for d in dc], int)


# the build refuses to run if it cannot reproduce the frozen ladder exactly
for tag, v in A["A1"]["ceilings"].items():
    pw = 1.0 if tag.startswith("classical") else 0.5
    dc, rk = rank_at(v["icc"], pw)
    assert abs(dc - v["dC_ceiling"]) < 1e-12 and rk == v["rank_ceiling"], (tag, dc, rk, v)
for tag, v in A["A1c"]["ceilings"].items():
    dc, rk = rank_at(v["icc"], EXP_MEASURED)
    assert abs(dc - v["dC_ceiling"]) < 1e-12 and rk == v["rank_ceiling"], (tag, dc, rk, v)
print(f"reproduced all {len(A['A1']['ceilings']) + len(A['A1c']['ceilings'])} frozen ceilings "
      f"to machine precision")

# ------------------------------------------------------------------ the scenario ladder
LAWS = [("sqrt", 0.5, "Square-root law", BLUE, "o"),
        ("measured", EXP_MEASURED, f"Measured law, {EXP_MEASURED:.3f}", BLUE_MID, "s"),
        ("classical", 1.0, "Classical law, 1.0, refuted", GREY_PALE, "^")]

SCEN = [  # (row label, ICC, the tags that use it, a note)
    ("Observed, 415 repeat nights, ICC 0.482", ICC_LAB),
    ("All 1,664 clean pairs, ICC 0.414", 0.414),
    ("Best single window, ICC 0.497", A["A1"]["ceilings"]["best_window_point_estimate"]["icc"]),
    ("Tightest window, under 90 days, ICC 0.452",
     A["A1"]["ceilings"]["tightest_window_under_90d"]["icc"]),
    ("Weakest window, ICC 0.424", A["A1"]["ceilings"]["lowest_window_point_estimate"]["icc"]),
    ("Upper 95% bound of the best window, ICC 0.647",
     A["A1"]["ceilings"]["upper_95_of_best_window"]["icc"]),
    ("Lower 95% bound, all pairs, ICC 0.381",
     A["A1"]["ceilings"]["lower_95_all_pairs_most_generous"]["icc"]),
    ("Lower 95% bound of the weakest window, ICC 0.171",
     A["A1"]["ceilings"]["lower_95_any_window_most_generous"]["icc"]),
]
# which law was actually fitted at which ICC, so nothing is drawn that was not run
FITTED = {"sqrt": {round(v["icc"], 6) for k, v in A["A1"]["ceilings"].items()
                   if not k.startswith("classical")},
          "measured": {round(v["icc"], 6) for v in A["A1c"]["ceilings"].values()},
          "classical": {round(v["icc"], 6) for k, v in A["A1"]["ceilings"].items()
                        if k.startswith("classical")}}

LADDER = []
for lab, icc in SCEN:
    for key, pw, _, _, _ in LAWS:
        if round(icc, 6) not in FITTED[key]:
            continue
        dc, rk = rank_at(icc, pw)
        LADDER.append({"scenario": lab, "icc": icc, "law": key, "power": pw,
                       "dC_ceiling": dc, "rank": rk, "top10": bool(dc > TARGET10)})
        PROV.append({"panel": "B", "label": f"{lab} / {key}", "icc": icc, "power": pw,
                     "dC": dc, "rank": rk,
                     "source": "tst_ceiling/attack_summary.json A1 and A1c ceilings"})
LAD = pd.DataFrame(LADDER)
assert len(LAD) == 17, len(LAD)

DEFENSIBLE = LAD[LAD.law != "classical"]
BEST_DEF = int(DEFENSIBLE["rank"].min())
WORST_DEF = int(DEFENSIBLE["rank"].max())
TOP10 = LAD[LAD.top10]
print(f"defensible scenarios (square-root and measured law): n = {len(DEFENSIBLE)}, "
      f"rank {BEST_DEF} to {WORST_DEF}")
print(f"scenarios reaching the top 10: {len(TOP10)} -> "
      f"{list(zip(TOP10.scenario, TOP10.law, TOP10['rank']))}")
# v8 (lane R1, 2026-09-13): the typed v7 conclusion "exactly one scenario reaches the top 10 and it is on the classical
# law" is a result, so it is no longer asserted (integrity gate 6). It is printed old-versus-new (rule 11) and the
# figure names whatever scenarios the table puts in the top 10. Every other check above and below is kept.
V7_TOP10 = [("Lower 95% bound of the weakest window, ICC 0.171", "classical", 2)]   # the v7 finding, typed for the comparison only
V8_TOP10 = [(r.scenario, r.law, int(r["rank"])) for _, r in TOP10.iterrows()]
TOP10_DRAMATIC = sorted((s, l) for s, l, _ in V7_TOP10) != sorted((s, l) for s, l, _ in V8_TOP10)
print(f"TOP10 FINDING old-versus-new: v7 = {V7_TOP10} | v8 = {V8_TOP10} | "
      + ("DRAMATIC: the set of scenarios reaching the top 10 changed" if TOP10_DRAMATIC else "same scenarios (ranks may differ)"))
LAWTXT = {"sqrt": "square-root law", "measured": f"measured law, {EXP_MEASURED:.3f}", "classical": "classical law, refuted"}

# the reliability each law would need before perfect measurement reached the top 10
BREAK = {}
for key, pw, *_ in LAWS:
    fold = float(np.sqrt((TARGET10 - INTERCEPT) / SLOPE) / E_NOW)
    BREAK[key] = float(fold ** (-1.0 / pw))
print("reliability that would be required to reach the top 10:",
      {k: round(v, 4) for k, v in BREAK.items()})

# ------------------------------------------------------------------ the home construct
HAB_A2 = A["A2"]["variants"]
HOME = [("Nocturnal oxygen, for scale", float(S["reproduced"]["T90_dC"]),
         int(S["reproduced"]["T90_rank"]), BLUE),
        ("Laboratory total sleep time, whole cohort",
         HAB_A2["lab TST_min, full cohort (reference)"]["dC"],
         HAB_A2["lab TST_min, full cohort (reference)"]["rank_if_inserted"], GREY),
        ("Laboratory total sleep time, the same people",
         HAB_A2["lab TST_min, SAME subsample (control)"]["dC"],
         HAB_A2["lab TST_min, SAME subsample (control)"]["rank_if_inserted"], GREY),
        ("Stated habitual sleep at home",
         HAB_A2["habitual_min (stated home sleep, all tiers)"]["dC"],
         HAB_A2["habitual_min (stated home sleep, all tiers)"]["rank_if_inserted"], GREY_MID),
        ("Stated habitual sleep under 6 hours",
         HAB_A2["stated home sleep under 6 h, binary"]["dC"],
         HAB_A2["stated home sleep under 6 h, binary"]["rank_if_inserted"], GREY_MID)]
for lab, dc, rk, _ in HOME:
    PROV.append({"panel": "D", "label": lab, "icc": None, "power": None, "dC": dc, "rank": rk,
                 "source": "tst_ceiling/attack_summary.json A2 variants, "
                           "summary.json reproduced"})

R_POOL = COR[COR.subset == "all_tiers"].iloc[0]
R_CLEAN = COR[COR.subset == "avg_nightly_hours_field"].iloc[0]
ICC_HOME_RANDOM = float(A["A2b"]["icc_two_random_mentions"]) if "icc_two_random_mentions" in \
    A["A2b"] else float(A["A2b"].get("icc_random_two", 0.472939))

# ==============================================================================================
# The sheet
#   A  rank against assumed reliability, one curve per attenuation law
#   B  the seventeen fitted scenarios, and the one that reaches rank 2
#   C  which attenuation law the data actually show
#   D  the home measure itself, measured rather than assumed
# ==============================================================================================
from matplotlib.transforms import blended_transform_factory  # noqa: E402

plt.rcParams.update(rc())
W, H = 14.1, 10.30
fig = plt.figure(figsize=(W, H))

BOT, HLOW, GAPV, HTOP = 1.62, 2.62, 1.12, 3.50
LA, WA = 1.42, 4.70                 # panel A
LB, WB = 4.02, 2.95                 # panel B, labels to its left
LC, WC = 9.35, 3.30                 # panel C
LD, WD = 9.35, 3.30                 # panel D

axA = fig.add_axes([LA / W, (BOT + HLOW + GAPV) / H, WA / W, HTOP / H])
axB = fig.add_axes([LB / W, BOT / H, WB / W, HLOW / H])
axC = fig.add_axes([LC / W, (BOT + HLOW + GAPV) / H, WC / W, HTOP / H])
axD = fig.add_axes([LD / W, BOT / H, WD / W, HLOW / H])

# ---------------------------------------------------------------- Panel A, rank vs reliability
AY_TOP, AY_BOT = 0.94, 206.0
BX_LO, BX_HI = 1.0, 220.0
grid = np.linspace(0.10, 1.0, 1400)
axA.axhspan(AY_TOP, 10.0, color=GREY_PALE, alpha=0.55, lw=0, zorder=0)
for key, pw, lab, col, mk in LAWS:
    _, rk = rank_at(grid, pw)
    axA.plot(grid, rk, color=col, lw=2.2 if key != "classical" else 1.8,
             ls="-" if key == "sqrt" else ((0, (5, 2)) if key == "measured" else (0, (2, 1.8))),
             zorder=3 if key != "classical" else 2)

GUIDES = [(ICC_LAB, "0.482, one laboratory night, the observed value"),
          (0.414, "0.414, all 1,664 clean pairs"),
          (ICC_HOME_RANDOM, "0.473, stated habitual sleep, its own reliability"),
          (float(A["A1"]["ceilings"]["lower_95_any_window_most_generous"]["icc"]),
           "0.171, lower 95% bound of the weakest gap window")]
for x, lab in GUIDES:
    axA.axvline(x, color=INK, lw=0.8, ls=(0, (1, 2.6)), zorder=1)
axA.set_yscale("log")
axA.set_ylim(AY_BOT, AY_TOP)
axA.set_yticks([1, 2, 5, 10, 24, 50, 100, N_MEASURES])
axA.set_yticklabels(["1", "2", "5", "10", "24", "50", "100", str(N_MEASURES)], fontsize=TICK)
axA.yaxis.set_minor_formatter(plt.NullFormatter())
axA.tick_params(axis="y", which="minor", length=0)
axA.set_xlim(0.094, 1.006)
axA.set_xticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
axA.set_xlabel("Assumed reliability of one laboratory night of total sleep time",
               fontsize=LABEL, labelpad=4)
axA.set_ylabel(f"Rank at perfect measurement, of {N_MEAS}", fontsize=LABEL, labelpad=4)
axA.text(0.985, 10.0, "top 10", fontsize=SMALL, ha="right", va="bottom", color=GREY_MID)
axA.text(0.022, 0.975, "Reliabilities this cohort actually supplies\n"
         + "\n".join(g[1] for g in GUIDES),
         transform=axA.transAxes, fontsize=SMALL, ha="left", va="top", color=INK,
         linespacing=1.42)
axA.legend(handles=[Line2D([], [], color=c, lw=2.2, marker=m, ms=6.2,
                           markeredgecolor=INK if k == "classical" else "white",
                           ls="-" if k == "sqrt" else ((0, (5, 2)) if k == "measured"
                                                       else (0, (2, 1.8))), label=l)
                    for k, _, l, c, m in LAWS],
           loc="lower left", bbox_to_anchor=(0.005, 0.005), fontsize=SMALL, frameon=False,
           handlelength=2.8, labelspacing=0.35)
panel_label(axA, "A", dx=-0.16, dy=1.035)
axA.set_title("Correct for measurement error, and duration still lands nowhere",
              fontsize=LABEL, fontweight="bold", loc="left", pad=9)

# ---------------------------------------------------------------- Panel B, the fitted scenarios
labels = [s for s, _ in SCEN]
yB = np.arange(len(labels), dtype=float)[::-1]
ymap = {s: y for s, y in zip(labels, yB)}
axB.axvspan(BX_LO, 10.0, color=GREY_PALE, alpha=0.55, lw=0, zorder=0)
LAWDY = {"sqrt": 0.20, "measured": 0.0, "classical": -0.20}
for _, r in LAD.iterrows():
    col, mk = next((c, m) for k, _, _, c, m in LAWS if k == r.law)
    axB.scatter([r["rank"]], [ymap[r.scenario] + LAWDY[r.law]], s=52, marker=mk, zorder=3,
                facecolors=col, edgecolors=INK if r.law == "classical" else "white",
                linewidths=0.8)
axB.set_xscale("log")
axB.set_xlim(BX_LO, BX_HI)
axB.set_xticks([1, 2, 5, 10, 24, 50, 100, N_MEASURES])
axB.set_xticklabels(["1", "2", "5", "10", "24", "50", "100", str(N_MEASURES)], fontsize=TICK)
axB.xaxis.set_minor_formatter(plt.NullFormatter())
axB.tick_params(axis="x", which="minor", length=0)
axB.set_yticks(yB)
axB.set_yticklabels(labels, fontsize=TICK)
bare_y(axB)
axB.set_ylim(yB.min() - 0.75, yB.max() + 0.95)
axB.set_xlabel(f"Rank at perfect measurement, of {N_MEAS}", fontsize=LABEL, labelpad=4)
axis_covers(axB, LAD["rank"].tolist(), "x", "eFigure15 panel B")
axB.text(0.005, 0.995, "Marker and colour as in A", transform=axB.transAxes, fontsize=SMALL,
         ha="left", va="top", color=GREY_MID)
# every scenario the table puts in the top 10 is named from the table (v8, lane R1 2026-09-13): law and rank, never typed
for _, bad in TOP10.iterrows():
    axB.annotate(f"reaches the top 10 (rank {int(bad['rank'])}),\non the {LAWTXT[bad.law]}",
                 xy=(bad["rank"], ymap[bad.scenario] + LAWDY[bad.law]),
                 xytext=(1.12, ymap[bad.scenario] + 0.72),
                 fontsize=SMALL, color=INK, ha="left", va="center", linespacing=1.3,
                 arrowprops=dict(arrowstyle="-", color=INK, lw=0.8,
                                 connectionstyle="angle3,angleA=0,angleB=90"))
if len(TOP10) == 0:
    axB.text(1.12, yB.max() + 0.72, "no fitted scenario reaches the top 10", fontsize=SMALL, color=INK,
             ha="left", va="center")
axB.text(0.995, 0.02, f"every defensible scenario, rank {BEST_DEF} to {WORST_DEF}",
         transform=axB.transAxes, fontsize=SMALL, ha="right", va="bottom", color=GREY_MID)
panel_label(axB, "B", dx=-1.34, dy=1.07)
axB.set_title(f"All {len(LAD)} fitted scenarios", fontsize=LABEL, fontweight="bold", loc="left",
              pad=9)

# ---------------------------------------------------------------- Panel C, the measured law
FEAT = {"AHI": "Apnea-hypopnea index", "spo2_nadir_corrected": "Oxygen nadir",
        "spo2_pct_below_90": "Time below 90%"}
MK = {"AHI": "o", "spo2_nadir_corrected": "s", "spo2_pct_below_90": "^"}
rs = np.linspace(0.15, 1.0, 400)
axC.plot(rs, rs ** 0.5, color=BLUE, lw=2.0, ls="-", zorder=2)
axC.plot(rs, rs ** 1.0, color=GREY_PALE, lw=1.8, ls=(0, (2, 1.8)), zorder=2)
for f, meta in A["A1b"]["per_feature"].items():
    e = meta["E_by_reliability"]
    x = np.array(sorted(float(k) for k in e))
    base = e[[k for k in e if abs(float(k) - 1.0) < 1e-9][0]]
    y = np.array([e[[k for k in e if abs(float(k) - xi) < 1e-9][0]] / base for xi in x])
    axC.plot(x, y, color=GREY, lw=1.0, ls="-", alpha=0.7, zorder=3)
    axC.scatter(x, y, s=42, marker=MK[f], facecolors="white", edgecolors=GREY,
                linewidths=1.3, zorder=4)
    PROV.append({"panel": "C", "label": f"{FEAT[f]} exponent", "icc": None,
                 "power": meta["exponent"], "dC": None, "rank": None,
                 "source": "tst_ceiling/attack_summary.json A1b per_feature"})
axC.set_xlim(0.13, 1.03)
axC.set_ylim(0.13, 1.10)
axC.set_xlabel("Reliability the measurement was degraded to", fontsize=LABEL, labelpad=4)
axC.set_ylabel("Log hazard kept, as a share of the whole", fontsize=LABEL, labelpad=4)
axC.legend(handles=[Line2D([], [], color=BLUE, lw=2.0, label="Square-root law, 0.5"),
                    Line2D([], [], color=GREY_PALE, lw=1.8, ls=(0, (2, 1.8)),
                           label="Classical law, 1.0")]
                   + [Line2D([], [], color=GREY, lw=1.0, marker=MK[f], ms=5.6,
                             markerfacecolor="white", markeredgecolor=GREY,
                             label=f"{FEAT[f]}, {A['A1b']['per_feature'][f]['exponent']:.3f}")
                      for f in FEAT],
           loc="upper left", bbox_to_anchor=(0.005, 0.995), fontsize=SMALL, frameon=False,
           handlelength=2.2, labelspacing=0.35)
for _a, _b in ((BLUE, GREY), (BLUE, GREY_PALE), (GREY, GREY_PALE),
               (BLUE, BLUE_MID), (BLUE_MID, GREY_PALE)):
    assert abs(luma(_a) - luma(_b)) >= 45.0, (_a, _b, abs(luma(_a) - luma(_b)))
panel_label(axC, "C", dx=-0.16, dy=1.035)
axC.set_title("Which attenuation law the data actually show", fontsize=LABEL, fontweight="bold",
              loc="left", pad=9)

# ---------------------------------------------------------------- Panel D, the home measure
yD = np.arange(len(HOME), dtype=float)[::-1]
axD.axvline(0.0, color=INK, lw=0.9, ls=(0, (4, 3)), zorder=1)
for (lab, dc, rk, col), yy in zip(HOME, yD):
    axD.barh(yy, dc, height=0.56, color=col, edgecolor="none", zorder=2)
    xt = dc + 0.0008 if dc >= 0 else 0.0008
    axD.text(xt, yy, f"rank {rk} of {N_MEAS}", va="center", ha="left", fontsize=SMALL,
             color=INK, zorder=3)
axD.set_yticks(yD)
axD.set_yticklabels([h[0] for h in HOME], fontsize=TICK)
bare_y(axD)
axD.set_xlim(-0.0135, 0.0255)
axD.set_xticks([-0.01, 0.0, 0.01, 0.02])
axD.set_xticklabels(["-0.010", "0", "+0.010", "+0.020"], fontsize=TICK)
axD.set_ylim(yD.min() - 0.72, yD.max() + 0.72)
axD.set_xlabel("Increase in held-out concordance over an age-and-sex model", fontsize=LABEL,
               labelpad=4)
axis_covers(axD, [h[1] for h in HOME], "x", "eFigure15 panel D")
panel_label(axD, "D", dx=-0.62, dy=1.07)
axD.set_title("The home measure itself, measured rather than assumed", fontsize=LABEL,
              fontweight="bold", loc="left", pad=9)

# ---------------------------------------------------------------- titles and footer
fig.text(0.006, 1 - 0.22 / H,
         "People sleep differently at home than in the laboratory. Correct for that, and sleep "
         "duration still lands nowhere",
         fontsize=TITLE, fontweight="bold", ha="left", va="top", color=INK)
fig.text(0.006, 1 - 0.56 / H, textwrap.fill(
    f"One laboratory night of total sleep time has a repeat-night reliability of {ICC_LAB:.3f} "
    f"({N_PAIRS} treatment-free pairs). A and B ask where total sleep time would rank among all "
    f"{N_MEAS} measurements if that error were removed entirely, under each attenuation law. C "
    f"shows which law the cohort's own repeat nights support. D is the home construct itself: "
    f"stated habitual sleep correlates with the laboratory night at r = {R_CLEAN.pearson_r:+.3f} "
    f"on the clean template field (n = {int(R_CLEAN.n):,}) and r = {R_POOL.pearson_r:+.3f} pooled "
    f"(n = {int(R_POOL.n):,}), and it is fitted here rather than assumed.", 216),
    fontsize=SMALL, ha="left", va="top", color=GREY, linespacing=1.42)

# footer sentences derived from the table (v8, lane R1 2026-09-13), never typed
DEF_TOP_TXT = ("It never enters the top 10, and it never enters the top 20. " if BEST_DEF > 20 else
               "It never enters the top 10. " if BEST_DEF > 10 else
               f"It enters the top 10 in {int((DEFENSIBLE['rank'] <= 10).sum())} of these scenarios. ")
if len(TOP10) == 0:
    TOP10_TXT = f"No scenario of the {len(LAD)} reaches the top 10. The classical law is refuted either way: "
else:
    TOP10_TXT = (f"{len(TOP10)} scenario{'s' if len(TOP10) > 1 else ''} of the {len(LAD)} "
                 f"reach{'es' if len(TOP10) == 1 else ''} the top 10: "
                 + "; ".join(f"{r.scenario[0].lower() + r.scenario[1:]} on the {LAWTXT[r.law]} (exponent {r.power:g}) at rank {int(r['rank'])}"
                             for _, r in TOP10.iterrows()) + ". ")
    _CL = TOP10[(TOP10.law == "classical") & (TOP10.icc.round(6) == round(A["A1"]["ceilings"]["lower_95_any_window_most_generous"]["icc"], 6))]
    if len(_CL):
        TOP10_TXT += (f"The classical-law entry uses an ICC of {float(_CL.iloc[0]['icc']):.3f}, the lower 95% bound of the weakest of "
                      f"eight overlapping gap windows, which sits outside the whole-sample interval of "
                      f"[{A['A1']['icc_TST_ci'][0]:.3f}, {A['A1']['icc_TST_ci'][1]:.3f}]. Both halves of it are refuted: ")
    else:
        TOP10_TXT += "The classical law is refuted either way: "
foot = (
    f"WHAT THE RANGE IS. Across the {len(DEFENSIBLE)} scenarios built on an attenuation law this "
    f"cohort supports, the square-root law and the exponent of {EXP_MEASURED:.3f} measured in C, "
    f"total sleep time at perfect measurement lands between rank {BEST_DEF} and rank {WORST_DEF} "
    f"of {N_MEAS}. " + DEF_TOP_TXT + TOP10_TXT
    + f"the exponent measured on the three measurements that have repeat nights is "
    f"{A['A1b']['per_feature']['AHI']['exponent']:.3f}, "
    f"{A['A1b']['per_feature']['spo2_nadir_corrected']['exponent']:.3f} and "
    f"{A['A1b']['per_feature']['spo2_pct_below_90']['exponent']:.3f}, mean "
    f"{EXP_MEASURED:.3f}, and reaching the top 10 under that exponent would need a reliability of "
    f"{BREAK['measured']:.3f} or worse, against a measured {ICC_LAB:.3f}. The Discussion clause "
    f"should therefore read no better than rank {BEST_DEF}, not {int(DEFENSIBLE[DEFENSIBLE.law == 'sqrt']['rank'].min())}. "
    f"D carries one omission by design: stated habitual sleep restricted to the highest-quality "
    f"records scores {HAB_A2['habitual_min, tier A only']['dC']:+.4f} and also ranks "
    f"{HAB_A2['habitual_min, tier A only']['rank_if_inserted']}, but it is fitted on a single "
    f"condition in {HAB_A2['habitual_min, tier A only']['n_patients']:,} people, so it is stated "
    f"here rather than drawn on the same axis as the stable fits."
)
fig.text(0.006, 0.055 / H, textwrap.fill(foot, 214), fontsize=SMALL, ha="left", va="bottom",
         color=GREY, linespacing=1.42)

out = save(fig, "eFigure15_duration_reliability_ceiling", kind="Supplementary")

json.dump({"n_measures": N_MEAS, "icc_lab_night": ICC_LAB, "n_repeat_pairs": N_PAIRS,
           "measured_attenuation_exponent": EXP_MEASURED,
           "exponent_per_feature": {k: v["exponent"] for k, v in A["A1b"]["per_feature"].items()},
           "ladder": LAD.to_dict("records"),
           "defensible_rank_range": [BEST_DEF, WORST_DEF],
           "n_scenarios_reaching_top10": int(len(TOP10)),
           "scenarios_reaching_top10": [{"scenario": s, "law": l, "rank": r} for s, l, r in V8_TOP10],
           "top10_finding_old_vs_new": {"v7": [{"scenario": s, "law": l, "rank": r} for s, l, r in V7_TOP10],
                                        "v8": [{"scenario": s, "law": l, "rank": r} for s, l, r in V8_TOP10],
                                        "dramatic": bool(TOP10_DRAMATIC)},
           "reliability_needed_for_top10": BREAK,
           "home": [{"label": l, "dC": d, "rank": r} for l, d, r, _ in HOME],
           "lab_vs_habitual_r": {"clean_template_field": [float(R_CLEAN.pearson_r),
                                                          int(R_CLEAN.n)],
                                 "pooled": [float(R_POOL.pearson_r), int(R_POOL.n)]},
           "tier_a_omitted_from_panel_D": {
               "dC": HAB_A2["habitual_min, tier A only"]["dC"],
               "rank": HAB_A2["habitual_min, tier A only"]["rank_if_inserted"],
               "n_patients": HAB_A2["habitual_min, tier A only"]["n_patients"],
               "n_outcomes_fitted": HAB_A2["habitual_min, tier A only"]["n_out"]},
           "reliability_of_the_home_measure": {
               "first_two_mentions": A["A2b"]["icc_first_two"],
               "two_random_mentions": A["A2b"]["icc_random_two"],
               "n": A["A2b"]["n_patients_with_repeats"]},
           "icc_all_pairs_95ci": A["A1"]["icc_TST_ci"],
           "icc_by_gap_window": A["A1"]["by_gap"]},
          open(f"{WORK}/eFigure15_drawn_values.json", "w"), indent=1)
# Same reason as Figure27: freeze the inputs where the manuscript's own audit can see them.
json.dump({"source": "Sleep_Variability_2026-08/tst_ceiling/attack_summary.json",
           "n_measures": N_MEAS, "icc_one_lab_night": ICC_LAB, "n_repeat_pairs": N_PAIRS,
           "icc_all_pairs_95ci": A["A1"]["icc_TST_ci"],
           "icc_by_gap_window": A["A1"]["by_gap"],
           "attenuation_exponent": {k: v["exponent"] for k, v in A["A1b"]["per_feature"].items()},
           "attenuation_exponent_mean": EXP_MEASURED,
           "effect_to_gain_map": {"slope": SLOPE, "intercept": INTERCEPT,
                                  "E_meanabs_TST": E_NOW, "dC_rank10": TARGET10},
           "ladder": LAD.to_dict("records"),
           "defensible_rank_range": [BEST_DEF, WORST_DEF],
           "reliability_needed_for_top10": BREAK,
           "home_measure": {l: {"dC": d, "rank": r} for l, d, r, _ in HOME},
           "home_measure_tier_a_only": HAB_A2["habitual_min, tier A only"],
           "reliability_of_the_home_measure": {"first_two": A["A2b"]["icc_first_two"],
                                               "two_random": A["A2b"]["icc_random_two"],
                                               "n": A["A2b"]["n_patients_with_repeats"]},
           "lab_vs_habitual_r": {"clean_template_field":
                                 {"r": float(R_CLEAN.pearson_r), "n": int(R_CLEAN.n)},
                                 "pooled": {"r": float(R_POOL.pearson_r), "n": int(R_POOL.n)}}},
          open(f"{NUM}/tst_ceiling.json", "w"), indent=1)
pd.DataFrame(PROV).to_csv(f"{WORK}/eFigure15_provenance.csv", index=False)
print(f"\nwrote {out}")
print(f"{len(PROV)} plotted values, every one from tst_ceiling or habitual_sleep")
