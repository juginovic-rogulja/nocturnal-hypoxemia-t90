"""Supp Fig 7, lane V14_L3b_DURATION_SUPP, round 37 (2026-09-15): re-plot the three panels at the V13 sheet's geometry from
the v8.1 tst_ceiling outputs (steps 141 to 143, T90_FULL_NIGHTS_ONLY=1, 15,551 full diagnostic nights) and numbers/ranking_v3.csv
(step 105). Repointed copy of the round-30 builder (01_build_PRE_V8_1.py beside it); the V13 design is the round-31 LTEXT
sheet: the six prose spans inside the panels ("415 untreated same-person pairs", "a median 476 days apart", "observed reliability
0.482", "rank 50 and 45 at perfect", "measurement", "rank 12 to 68 across all 14 supported scenarios") are NOT drawn; the statistic
labels stay ("median .. min", "IQR ..", the ICC line, "rank 1", "rank N", "top 10", the key).

v8.1 changes folded in: the number of measurements comes from the files (141, no literal); the rank axes end at N_MEAS with the
top limit scaled by the round-30 ratio (206 / 197); the "Total sleep time, as measured" rank is the paper's canonical
ranking_v3.csv rank (step 105); the ceiling frame's own rank (ranking_reproduced.csv, step 141) is logged beside it (the two files
rank on different refits under v8.1, see the report). The pairs file Sleep_Variability_2026-08/data/untreated_wide.parquet is the
Aug 7 pilot extraction step 141 reads for the ICC (no sidecar, unchanged, flagged as in round 30).

Writes work/sheet.pdf (the whole sheet, letters excluded), work/expected_values.json, work/build_log.txt."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import json
import os
import sys

import numpy as np
import pandas as pd

LANE = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LSUPP-B"   # R49: SD = LANE/Supp_Fig07 = this sheet folder
sys.path.insert(0, f"{LANE}/Supp_Fig07/scripts")
from l3b_common import (NUM, SV, INK, LABF, TCKF, ANNF, FLOOR, PT_PLUS, MEDGE, S_CIRCLE, S_SQUARE, rc_polish,  # noqa: E402
                        measure, text_gate, sidecar, sha256, jload, hydrated, plt)
import matplotlib.ticker as mticker  # noqa: E402

SHEET = "Supp_Fig07"
SD = f"{LANE}/{SHEET}"
WORK = f"{SD}/work"
os.makedirs(WORK, exist_ok=True)
LOG = open(f"{WORK}/build_log.txt", "w")


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    LOG.write(s + "\n")


# ------------------------------------------------------------------ colours sampled from the base sheet (V13 = the round-30 drawing layer)
BAR = "#ccd1d6"
IQR_BAND = "#e8edef"
TOP10_BAND = "#e6e8ea"
RANGE_BAR = "#afc0c9"
GREEN = "#298d32"
GREEN_MID = "#39c445"
GREEN_LT = "#79d475"
GREEN_TXT = "#5dcf66"
BLUE_MID = "#3f9fd8"

# ------------------------------------------------------------------ sources, gated (v8.1 sidecars)
TC = f"{SV}/tst_ceiling"
SRC = {}
for f in ("attack_summary.json", "summary.json", "_iccs.json", "ranking_reproduced.csv", "effect_vs_gain.csv"):
    sh, sc = sidecar(f"{TC}/{f}")
    SRC[f] = {"path": f"{TC}/{f}", "sha256": sh, "step": sc["step"], "mtime": sc["output"]["mtime_local"]}
    log(f, "sha256", sh, "step", sc["step"]["id"], "rc", sc["step"].get("rc"), "stamped", sc["output"]["mtime_local"])
sh, sc = sidecar(f"{NUM}/ranking_v3.csv")
SRC["ranking_v3.csv"] = {"path": f"{NUM}/ranking_v3.csv", "sha256": sh, "step": sc["step"], "mtime": sc["output"]["mtime_local"]}
log("ranking_v3.csv sha256", sh, "step", sc["step"]["id"], "stamped", sc["output"]["mtime_local"])
sh43, sc43 = sidecar(f"{NUM}/tst_ceiling.json")
SRC["tst_ceiling.json"] = {"path": f"{NUM}/tst_ceiling.json", "sha256": sh43, "step": sc43["step"], "mtime": sc43["output"]["mtime_local"]}
PAIRS = f"{SV}/data/untreated_wide.parquet"
SRC["untreated_wide.parquet"] = {"path": PAIRS, "sha256": sha256(PAIRS), "sidecar": None,
                                 "note": "Aug 7 pilot extraction, no sidecar, read by step 141 model.py for the ICC; unchanged since round 30 (same sha256)"}
log("untreated_wide.parquet sha256", SRC["untreated_wide.parquet"]["sha256"], "(no sidecar, static input, flagged)")

A = jload(f"{TC}/attack_summary.json")
S = jload(f"{TC}/summary.json")
ICC = jload(f"{TC}/_iccs.json")
RK = pd.read_csv(hydrated(f"{TC}/ranking_reproduced.csv"))
EVG = pd.read_csv(hydrated(f"{TC}/effect_vs_gain.csv"))
T43 = jload(f"{NUM}/tst_ceiling.json")
V3 = pd.read_csv(hydrated(f"{NUM}/ranking_v3.csv"), comment="#")

w = pd.read_parquet(hydrated(PAIRS))
DIFF = (w.tst_min_n1 - w.tst_min_n2).abs()
N_PAIRS = int(np.isfinite(DIFF).sum())
MED = float(DIFF.median())
Q1, Q3 = float(DIFF.quantile(.25)), float(DIFF.quantile(.75))
ICC_LAB = float(ICC["TST"]["icc_raw"])
ICC_CI = [float(v) for v in A["A1"]["icc_TST_ci"]]
GAP_MED = float(A["A1"]["gap_median_days"])
assert N_PAIRS == int(ICC["TST"]["n_pairs"]) == int(A["A1"]["n_pairs"]) == int(T43["n_repeat_pairs"]), (N_PAIRS, ICC["TST"]["n_pairs"], A["A1"]["n_pairs"])
assert abs(ICC_LAB - float(A["A1"]["icc_TST_all_pairs"])) < 1e-12 and abs(ICC_LAB - float(T43["icc_one_lab_night"])) < 1e-12
assert all(abs(a - b) < 1e-12 for a, b in zip(ICC_CI, T43["icc_all_pairs_95ci"]))


def icc_oneway(a, b):
    x = np.vstack([np.asarray(a, float), np.asarray(b, float)])
    n = x.shape[1]
    gm = x.mean()
    msb = 2.0 * ((x.mean(axis=0) - gm) ** 2).sum() / (n - 1)
    msw = ((x - x.mean(axis=0)) ** 2).sum() / n
    return (msb - msw) / (msb + msw)


ICC_RECOMP = icc_oneway(w.tst_min_n1, w.tst_min_n2)
log(f"pairs {N_PAIRS}, median diff {MED:.2f}, IQR {Q1:.2f} to {Q3:.2f}, ICC file {ICC_LAB:.6f}, ICC recomputed one-way {ICC_RECOMP:.6f}, gap median {GAP_MED}")

# ------------------------------------------------------------------ the ceiling arithmetic, verified (no literal counts)
N_MEAS = len(RK)
SLOPE = float(S["map_primary"]["slope"])
INTERCEPT = float(S["map_primary"]["intercept"])
E_NOW = float(EVG[EVG.feature == "TST_min"].E_meanabs.iloc[0])
TARGET10 = float(S["reproduced"]["dC_rank10"])
EXP_MEASURED = float(A["A1c"]["exponent_used"])
assert N_MEAS == int(S["n_measures"]) == int(T43["n_measures"]) == len(V3), (N_MEAS, S["n_measures"], T43["n_measures"], len(V3))
assert abs(EXP_MEASURED - float(T43["attenuation_exponent_mean"])) < 1e-12


def rank_at(icc, power):
    e = E_NOW * (1.0 / np.asarray(icc, float)) ** power
    dc = SLOPE * e ** 2 + INTERCEPT
    if np.ndim(dc) == 0:
        return float(dc), int((RK.dC > dc).sum() + 1)
    return dc, np.array([(RK.dC > v).sum() + 1 for v in dc], int)


n_rep = 0
for tag, v in A["A1"]["ceilings"].items():
    pw = 1.0 if tag.startswith("classical") else 0.5
    dc, rk = rank_at(v["icc"], pw)
    assert abs(dc - v["dC_ceiling"]) < 1e-12 and rk == v["rank_ceiling"], (tag, dc, rk, v)
    n_rep += 1
for tag, v in A["A1c"]["ceilings"].items():
    dc, rk = rank_at(v["icc"], EXP_MEASURED)
    assert abs(dc - v["dC_ceiling"]) < 1e-12 and rk == v["rank_ceiling"], (tag, dc, rk, v)
    n_rep += 1
log(f"reproduced all {n_rep} ceilings of attack_summary.json to machine precision")

SQRT_R = sorted(v["rank_ceiling"] for k, v in A["A1"]["ceilings"].items() if not k.startswith("classical"))
MEAS_R = sorted(v["rank_ceiling"] for v in A["A1c"]["ceilings"].values())
CLAS_R = sorted(v["rank_ceiling"] for k, v in A["A1"]["ceilings"].items() if k.startswith("classical"))
N_DEF = len(SQRT_R) + len(MEAS_R)
BEST_DEF, WORST_DEF = min(SQRT_R + MEAS_R), max(SQRT_R + MEAS_R)
assert [BEST_DEF, WORST_DEF] == list(T43["defensible_rank_range"]), (BEST_DEF, WORST_DEF, T43["defensible_rank_range"])
N_TOP10_ALL = sum(r <= 10 for r in SQRT_R + MEAS_R + CLAS_R)
OBS_KEY = [k for k in A["A1"]["ceilings"] if k.startswith("published_")]
assert len(OBS_KEY) == 1 and OBS_KEY[0] in A["A1c"]["ceilings"], OBS_KEY
OBS_KEY = OBS_KEY[0]
assert abs(float(A["A1"]["ceilings"][OBS_KEY]["icc"]) - ICC_LAB) < 1e-9, "the 'published' ceiling is not the observed ICC"
_, R_OBS_SQRT = rank_at(ICC_LAB, 0.5)
_, R_OBS_MEAS = rank_at(ICC_LAB, EXP_MEASURED)
assert R_OBS_SQRT == A["A1"]["ceilings"][OBS_KEY]["rank_ceiling"]
assert R_OBS_MEAS == A["A1c"]["ceilings"][OBS_KEY]["rank_ceiling"]
R_TST_RK = int(RK[RK.feature == "TST_min"]["rank"].iloc[0])
R_T90_RK = int(RK[RK.feature == "spo2_pct_below_90"]["rank"].iloc[0])
R_TST = int(V3[V3.feature == "TST_min"]["rank"].iloc[0])            # the paper's canonical rank (step 105)
R_T90 = int(V3[V3.feature == "spo2_pct_below_90"]["rank"].iloc[0])
assert R_TST_RK == int(S["reproduced"]["TST_rank"]) and R_T90_RK == int(S["reproduced"]["T90_rank"])
N_RANK_DISAGREE = int((V3.merge(RK, on="feature", suffixes=("_v3", "_rk")).eval("rank_v3 != rank_rk")).sum())
log(f"canonical ranks (ranking_v3.csv): TST {R_TST}, T90 {R_T90}; ceiling-frame ranks (ranking_reproduced.csv): TST {R_TST_RK}, T90 {R_T90_RK}; "
    f"{N_RANK_DISAGREE} of {N_MEAS} ranks differ between the two files (spearman_vs_frozen {S['reproduced'].get('spearman_vs_frozen')})")
BREAK_MEAS = float((np.sqrt((TARGET10 - INTERCEPT) / SLOPE) / E_NOW) ** (-1.0 / EXP_MEASURED))
EXPS = {k: v["exponent"] for k, v in A["A1b"]["per_feature"].items()}
log(f"at observed reliability sqrt {R_OBS_SQRT} measured {R_OBS_MEAS}, supported range {BEST_DEF} to {WORST_DEF} over {N_DEF}, classical {CLAS_R}, "
    f"exponents {EXPS}, mean {EXP_MEASURED:.4f}, reliability needed for top 10 (measured) {BREAK_MEAS:.4f}")

EXPECTED = {"sheet": SHEET, "sources": SRC, "values": [], "facts_for_legend": {
    "n_pairs": N_PAIRS, "median_diff_min": MED, "iqr": [Q1, Q3], "icc": ICC_LAB, "icc_ci": ICC_CI, "gap_median_days": GAP_MED,
    "exponents": EXPS, "exponent_mean": EXP_MEASURED, "rank_tst_canonical_ranking_v3": R_TST, "rank_tst_ceiling_frame_ranking_reproduced": R_TST_RK,
    "rank_t90": R_T90, "rank_t90_ceiling_frame": R_T90_RK, "n_ranks_differ_between_ranking_v3_and_ranking_reproduced": N_RANK_DISAGREE,
    "rank_obs_sqrt": R_OBS_SQRT, "rank_obs_meas": R_OBS_MEAS, "supported_range": [BEST_DEF, WORST_DEF], "n_supported": N_DEF,
    "sqrt_ranks": SQRT_R, "measured_ranks": MEAS_R, "classical_ranks": CLAS_R,
    "n_scenarios_reaching_top10_all": N_TOP10_ALL, "n_scenarios_all": len(SQRT_R) + len(MEAS_R) + len(CLAS_R),
    "reliability_needed_for_top10_measured": BREAK_MEAS, "icc_recomputed_oneway": ICC_RECOMP, "n_measures": N_MEAS,
    "prose_not_drawn_V13_design": ["%d untreated same-person pairs" % N_PAIRS, "a median %.0f days apart" % GAP_MED,
                                   "observed reliability %.3f" % ICC_LAB, "rank %d and %d at perfect" % (R_OBS_SQRT, R_OBS_MEAS), "measurement",
                                   "rank %d to %d across all %d supported scenarios" % (BEST_DEF, WORST_DEF, N_DEF)]}}


def expect(panel, text, source, key, value, rule="text"):
    EXPECTED["values"].append(dict(panel=panel, text=text, source=source, key=key, value=value, rule=rule))


# ------------------------------------------------------------------ geometry (the builder's, re-measured on the base sheet in round 30)
# R40 follow-up (audit): the 24 tick marked the measured projection's rank at the observed reliability; with one projection the tick is the
# square-root projection's rank at the observed reliability, R_OBS_SQRT (23), read from the files, never typed
RANK_TICKS = [t for t in (1, 10, R_OBS_SQRT, 50, 100) if np.log10(N_MEAS / t) >= 0.2] + [N_MEAS]   # the round-30 ticks (24 -> R_OBS_SQRT) plus the end tick; a tick within 0.2 decades of the end tick (100 when N is 141) would overprint it and is dropped
RANK_TOP = 206.0 / 197.0 * N_MEAS        # the round-30 axis end (206 for 197 measurements), scaled to the v8.1 count
RANK_LO = 0.94
W_IN = 171 / 25.4
EDGE = 0.18
AXL = 0.82
AXR = W_IN - 0.30
TOP_PAD = 0.17
LTR_H = 0.25
GAP = 0.36
H_A, H_B, H_C = 2.18, 2.36, 1.55
XL_H = 0.48
yA = TOP_PAD + LTR_H
yB = yA + H_A + XL_H + GAP + LTR_H
yC = yB + H_B + XL_H + GAP + LTR_H
H = yC + H_C + XL_H + 0.14
assert abs(W_IN * 72 - 484.72) < 0.02 and abs(H * 72 - 670.32) < 0.05, (W_IN * 72, H * 72)


def rank_axis(ax, which="x"):
    axis = ax.xaxis if which == "x" else ax.yaxis
    (ax.set_xscale if which == "x" else ax.set_yscale)("log")
    (ax.set_xticks if which == "x" else ax.set_yticks)(RANK_TICKS)
    axis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:g}"))
    axis.set_minor_locator(mticker.NullLocator())


def bare_y(ax):
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)


with plt.rc_context(rc_polish()):
    fig = plt.figure(figsize=(W_IN, H))
    axA = fig.add_axes([AXL / W_IN, 1 - (yA + H_A) / H, (AXR - AXL) / W_IN, H_A / H])
    axB = fig.add_axes([AXL / W_IN, 1 - (yB + H_B) / H, (AXR - AXL) / W_IN, H_B / H])

    # ---------------------------------------------------- a, the two nights disagree
    BINS = np.arange(0, 345, 15)
    assert DIFF.max() < BINS[-1], float(DIFF.max())
    counts, _ = np.histogram(DIFF, bins=BINS)
    assert int(counts.sum()) == N_PAIRS
    Y_TOP_A = 76 if counts.max() <= 76 else int(np.ceil(counts.max() * 1.05))
    axA.bar(BINS[:-1], counts, width=15, align="edge", color=BAR, edgecolor="white", linewidth=0.6, zorder=2)
    axA.axvspan(Q1, Q3, color=IQR_BAND, lw=0, zorder=1)
    axA.axvline(MED, color=GREEN, lw=1.8, zorder=3)
    axA.text(MED + 8, 73.5, f"median {MED:.1f} min", fontsize=ANNF, color=INK, ha="left", va="top")
    axA.text(MED + 8, 65.5, f"IQR {Q1:.1f} to {Q3:.1f}", fontsize=ANNF, color=INK, ha="left", va="top")
    # V13 design (round 31): the two prose lines of the block are gone, the ICC statistic keeps its line (the block's third line)
    axA.text(0.985, 0.70, f"\n\nICC {ICC_LAB:.3f} (95% CI, {ICC_CI[0]:.3f}-{ICC_CI[1]:.3f})",
             transform=axA.transAxes, fontsize=ANNF, color=INK, ha="right", va="top", linespacing=1.55)
    expect("a", f"median {MED:.1f} min", PAIRS, "median |tst_min_n1 - tst_min_n2| (no sidecar, static pilot file)", MED)
    expect("a", f"IQR {Q1:.1f} to {Q3:.1f}", PAIRS, "quartiles of |tst_min_n1 - tst_min_n2|", [Q1, Q3])
    expect("a", f"ICC {ICC_LAB:.3f} (95% CI, {ICC_CI[0]:.3f}-{ICC_CI[1]:.3f})", f"{TC}/_iccs.json",
           "TST/icc_raw and attack_summary.json A1/icc_TST_ci", [ICC_LAB] + ICC_CI)
    axA.set_xlim(0, 330)
    axA.set_xticks([0, 60, 120, 180, 240, 300])
    axA.set_ylim(0, Y_TOP_A)
    axA.set_yticks([0, 25, 50, 75])
    axA.set_xlabel("Between-night difference in total sleep time, min", fontsize=LABF, labelpad=5)
    axA.set_ylabel("Pairs of nights", fontsize=LABF)
    EXPECTED["hist_counts"] = counts.tolist()

    # ---------------------------------------------------- b, remove the error entirely
    grid = np.linspace(0.10, 1.0, 1400)
    _, r_sqrt = rank_at(grid, 0.5)
    _, r_meas = rank_at(grid, EXP_MEASURED)
    axB.axhspan(RANK_LO, 10.0, color=TOP10_BAND, lw=0, zorder=0)
    axB.text(0.985, 0.865, "top 10", transform=axB.transAxes, fontsize=FLOOR, color=GREEN_TXT, ha="right", va="top")
    # R40 (Alen, 2026-09-18): panel b shows ONE solid line, the square-root projection; the measured projection (dashed, blue), its
    # legend entry, the legend itself and the blue square at the observed reliability leave; the green dot at the observed reliability
    # (rank 23) and the dotted reliability line stay. The measured exponent goes to the figure legend (facts_for_legend).
    axB.plot(grid, r_sqrt, color=GREEN, lw=2.0, zorder=3)
    axB.axvline(ICC_LAB, color=INK, lw=0.9, ls=(0, (1, 2.6)), zorder=1)
    axB.scatter([ICC_LAB], [R_OBS_SQRT], s=S_CIRCLE, color=GREEN_MID, zorder=4, edgecolors="white", linewidths=MEDGE)
    # V13 design (round 31): no annotation beside the observed-reliability dots
    rank_axis(axB, "y")
    axB.set_ylim(RANK_TOP, RANK_LO)
    axB.set_xlim(0.094, 1.006)
    axB.set_xticks([0.2, 0.4, 0.6, 0.8, 1.0])
    axB.set_xlabel("Assumed reliability of one laboratory night", fontsize=LABF, labelpad=5)
    axB.set_ylabel(f"Rank at perfect measurement, of {N_MEAS}", fontsize=LABF)
    # R49: the rotated title overhangs its axis at both ends; at 12 pt its top end (y 251) ran into the 14 pt letter b (bottom y 258, x 13 to 22).
    # Nudged 12 pt down its axis (the one element, size kept): top end at y 263, 5 pt clear of the letter; the bottom end stays clear of the x tick row.
    YLAB_B_NUDGE_PT = 12.0
    axB.yaxis.label.set_y(0.5 - YLAB_B_NUDGE_PT / (H_B * 72.0))
    EXPECTED["r49_panel_b_ylabel_nudge_pt"] = YLAB_B_NUDGE_PT
    # R40: no legend (one line; the curve is named in the figure legend)
    for v in (float(r_sqrt.min()), float(r_sqrt.max()), R_OBS_SQRT):
        assert RANK_LO <= v <= RANK_TOP, v
    EXPECTED["r40_panel_b"] = {"drawn": "square-root projection (exponent 0.5) solid green, dotted line at the observed reliability, green dot at rank R_OBS_SQRT",
                               "removed": ["measured projection line (dashed blue, exponent %.3f)" % EXP_MEASURED, "legend (two entries)", "blue square at the observed reliability (rank %d)" % R_OBS_MEAS],
                               "rank_obs_sqrt": R_OBS_SQRT, "rank_obs_meas": R_OBS_MEAS, "exponent_measured": EXP_MEASURED,
                               "rank_tick_followup": {"was": 24, "now": R_OBS_SQRT, "where": "panel b y axis and panel c x axis; panel c keeps one marker at R_OBS_SQRT on the supported-range bar"}}
    expect("b", f"Rank at perfect measurement, of {N_MEAS}", f"{TC}/summary.json", "n_measures", N_MEAS)
    expect("b", f"{N_MEAS}", f"{TC}/summary.json", "n_measures", N_MEAS, rule="int")   # the last y tick

    # ---------------------------------------------------- c, where duration lands
    rows = [("Nocturnal oxygen (T90), as measured", "point", R_T90, GREEN_MID, "o"),
            ("Total sleep time, as measured", "point", R_TST, GREEN_LT, "s"),
            ("Total sleep time, error removed entirely", "band", (BEST_DEF, WORST_DEF), BLUE_MID, "s")]
    # R49: the axis left edge is pinned at the V26 value (the gutter measured at the round-40 size, TCKF - PT_PLUS): the marks of c do not move.
    # The row labels (right-anchored at the axis) take TCKF when the longest still starts at or right of EDGE, else the column keeps its current size.
    gut = max(measure(fig, r[0], TCKF - PT_PLUS) for r in rows)
    AXL_C = EDGE + gut + 0.14
    gut_plus = max(measure(fig, r[0], TCKF) for r in rows)
    TCKF_C = TCKF if AXL_C - 4.0 / 72.0 - gut_plus >= EDGE - 1e-9 else TCKF - PT_PLUS
    axC = fig.add_axes([AXL_C / W_IN, 1 - (yC + H_C) / H, (AXR - AXL_C) / W_IN, H_C / H])
    log(f"panel c axis left edge {AXL_C * 72:.2f} pt (base 195.69, pinned); row labels at {TCKF_C} pt (at {TCKF} the longest would start {(AXL_C - 4.0 / 72.0 - gut_plus) * 72:.1f} pt from the page edge, margin EDGE {EDGE * 72:.1f} pt)")
    EXPECTED["r49_panel_c"] = {"axis_left_pt": AXL_C * 72, "row_label_size": TCKF_C, "kept": TCKF_C != TCKF, "longest_label_start_at_plus_pt": (AXL_C - 4.0 / 72.0 - gut_plus) * 72}
    yCq = np.arange(len(rows), dtype=float)[::-1]
    axC.axvspan(RANK_LO, 10.0, color=TOP10_BAND, lw=0, zorder=0)
    axC.text(9.4, yCq.max() + 0.60, "top 10", fontsize=FLOOR, color=GREEN_TXT, ha="right", va="top")
    for (lab, kind, v, col, mk), y in zip(rows, yCq):
        if kind == "point":
            axC.scatter([v], [y], s=S_CIRCLE if mk == "o" else S_SQUARE, marker=mk, color=col, zorder=4,
                        edgecolors="white", linewidths=MEDGE)
            if v > 100:
                axC.text(v * 0.90, y, f"rank {v}", fontsize=ANNF, color=INK, va="center", ha="right")
            else:
                axC.text(v * 1.16, y, f"rank {v}", fontsize=ANNF, color=INK, va="center", ha="left")
            expect("c", f"rank {v}", f"{NUM}/ranking_v3.csv", f"feature={'spo2_pct_below_90' if mk == 'o' else 'TST_min'}:rank", v)
        else:
            lo, hi = v
            axC.plot([lo, hi], [y, y], color=RANGE_BAR, lw=6.5, solid_capstyle="round", zorder=3)
            axC.scatter([R_OBS_SQRT], [y], s=S_SQUARE, marker="o", color=BLUE_MID, zorder=4,
                        edgecolors="white", linewidths=0.7)   # R40 follow-up: one marker, the square-root rank at the observed reliability (the 24 marker of the measured projection leaves)
            # V13 design (round 31): no sentence above the range bar
    axC.set_yticks(yCq)
    axC.set_yticklabels([r[0] for r in rows], fontsize=TCKF_C)   # R49
    bare_y(axC)
    axC.tick_params(axis="y", pad=4.0)
    rank_axis(axC, "x")
    axC.set_xlim(RANK_LO, RANK_TOP)
    axC.set_ylim(-0.62, len(rows) - 0.20)
    axC.set_xlabel(f"Rank among the {N_MEAS} sleep measurements", fontsize=LABF, labelpad=5)
    expect("c", f"Rank among the {N_MEAS} sleep measurements", f"{TC}/summary.json", "n_measures", N_MEAS)
    expect("c", f"{N_MEAS}", f"{TC}/summary.json", "n_measures", N_MEAS, rule="int")   # the last x tick
    for t in RANK_TICKS[:-1]:
        expect("b", f"{t}", "static:V13 rank-axis tick (position follows the rescaled axis)", "rank axis tick", f"{t}")
        expect("c", f"{t}", "static:V13 rank-axis tick (position follows the rescaled axis)", "rank axis tick", f"{t}")
    expect("c", "top 10", "static:V13 label at data x 9.4 of the rescaled rank axis", "top-10 strip label", "top 10")
    EXPECTED["rank_ticks"] = RANK_TICKS

    text_gate(fig)
    out = f"{WORK}/sheet.pdf"
    fig.savefig(out)
    plt.close(fig)
    log("wrote", out)

EXPECTED["page"] = (W_IN * 72, H * 72)
EXPECTED["n_meas"] = N_MEAS
EXPECTED["rank_axis_top"] = RANK_TOP
EXPECTED["r49_text_sizes"] = dict(pt_plus=PT_PLUS, LABF=LABF, TCKF=TCKF, ANNF=ANNF, FLOOR=FLOOR, letters=13.0 + PT_PLUS)
json.dump(EXPECTED, open(f"{WORK}/expected_values.json", "w"), indent=1, default=float)
log("wrote", f"{WORK}/expected_values.json", "with", len(EXPECTED["values"]), "expected printed values")
LOG.close()
