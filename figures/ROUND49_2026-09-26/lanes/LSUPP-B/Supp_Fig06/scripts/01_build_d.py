"""Supp Fig 6 part d, lane V14_L3b_DURATION_SUPP, round v8.2 (2026-09-16): rebuilt from the re-extracted step-139 outputs.

The owner ruled (v8.2) that no panel may keep old data, so the V13 copy of part d (v7 values) is replaced by a fresh draw of the
round-30 design (01_build_d_PRE_V8_1.py lineage, FF/redesign_eFigure13_supp_polish.py fig_d):
  column 1  the observed short-with-low against short-with-preserved contrast, results.csv (step 138, v8.1, 15,551 full nights)
  column 2  the 1-year landmark re-fit, attack_landmark_vs_primary.csv with attack.json A4_selection (step 139, v8.2, 2026-09-15 19:30)
  column 3  the smallest hazard ratio detectable at 80% power, attack_min_detectable_hr.csv (step 139, v8.2), with the dotted median
            over every cell2_vs_cell1 row of that file (the count is read from the file, never typed)
Every step-139 source must carry a provenance sidecar dated after 2026-09-15 19:18 (the v8.2 re-extraction); older files are refused.
Cohort proof: step 139 ran on the all-nights frame (attack.json rebuild.cohort_n) while step 138 ran on the full diagnostic nights,
but the split nights carry no sleep amount and fall in no cell, so the cell sizes (and the per-row events) must be identical between
the two files: asserted below before anything is drawn.

Geometry is the round-30 design: 171 mm part, GUT 0.18 in, AX_X0 1.92 in, ROW_IN 0.30 in, three key rows of 0.24 in, TOPBAND 1.06 in,
ticks 0.5, 1, 2, 4, 8 on an axis 0.25 to 10.8. The axis end is data-driven: when the dotted median (or a drawn diamond) lies beyond
10.8 the range steps to 0.25 to 24 with the tick 16 added (declared in the record, the CHANGES csv and the report).
The row pitch is not reduced. The part height follows the row count and the page is re-sized by 02_compose.py.

Writes work/part_d.pdf, work/expected_values_d.json, work/build_log_d.txt, and merges the part-d record into work/expected_values.json
(kept_from_v13 emptied, the pre-merge file kept as work/expected_values_PRE_V8_2.json)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import json
import os
import shutil
import sys

import numpy as np
import pandas as pd

LANE = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LSUPP-B"   # R49: SD = LANE/Supp_Fig06 = this sheet folder
sys.path.insert(0, f"{LANE}/Supp_Fig06/scripts")
from l3b_common import (SV, INK, LABF, TCKF, ANNF, PT_PLUS, GREY, BLUE, BLUE_MID, CI_LW, MEDGE, S_CIRCLE,  # noqa: E402
                        rc_polish, measure, text_gate, sidecar, jload, hydrated, plt)
from matplotlib import ticker as mticker  # noqa: E402

SHEET = "Supp_Fig06"
SD = f"{LANE}/{SHEET}"
WORK = f"{SD}/work"
os.makedirs(WORK, exist_ok=True)
LOG = open(f"{WORK}/build_log_d.txt", "w")
V8_2_CUT = "2026-09-15 19:18:00"


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    LOG.write(s + "\n")


BAND = "#eef0f1"
STARF = 10.0 + PT_PLUS   # R49 (unused in d)
KEPT = []   # R49: gutter labels kept at their current size
SRC = f"{SV}/dur_x_oxygen_healthy"

# ------------------------------------------------------------------ sources and their sidecars (step 138 = the observed column, step 139 = v8.2)
SOURCES = {}
for f, step in (("results.csv", "138"), ("summary.json", "138"), ("attack.json", "139"),
                ("attack_min_detectable_hr.csv", "139"), ("attack_landmark_vs_primary.csv", "139"),
                ("attack_landmark1y.csv", "139")):
    sh, sc = sidecar(f"{SRC}/{f}")
    assert sc["step"]["id"] == step and sc["step"]["rc"] == 0, (f, sc["step"])
    mt = sc["output"]["mtime_local"]
    if step == "139":
        assert mt >= V8_2_CUT, f"{f}: step-139 output dated {mt}, before the v8.2 re-extraction cut {V8_2_CUT}: refused"
    SOURCES[f] = {"path": f"{SRC}/{f}", "sha256": sh, "step": sc["step"], "mtime": mt}
    log(f"{f:32s} step {step} rc 0 stamped {mt} sha256 {sh}")

res = pd.read_csv(hydrated(f"{SRC}/results.csv"))
summ = jload(f"{SRC}/summary.json")
atk = jload(f"{SRC}/attack.json")
mdhr = pd.read_csv(hydrated(f"{SRC}/attack_min_detectable_hr.csv"))
lmvp = pd.read_csv(hydrated(f"{SRC}/attack_landmark_vs_primary.csv"))
lm1y = pd.read_csv(hydrated(f"{SRC}/attack_landmark1y.csv"))
BASE_TEXT = json.load(open(hydrated(f"{WORK}/base_text.json")))          # the V13 sheet's text layer (written by 01_build.py)
BASE_STRINGS = {s["text"] for s in BASE_TEXT}

# ------------------------------------------------------------------ cohort and cell proof: step 139 (all nights) against step 138 (full nights)
CCS = summ["cell_counts_side_by_side"]
RB = atk["rebuild"]
for key, akey in (("healthy_tst300", "cells_tst300"), ("healthy_tst240", "cells_tst240")):
    assert RB[akey] == CCS[key]["n_by_cell"], (key, RB[akey], CCS[key]["n_by_cell"])
COHORT_NOTE = {"step_138_cohort_n": summ["healthy_definition"]["n_cohort"], "step_138_healthy_n": summ["healthy_definition"]["n_healthy"],
               "step_139_cohort_n": RB["cohort_n"], "step_139_healthy_n": RB["healthy_n"],
               "cells_identical": True, "attack_self_flags": {k: RB[k] for k in ("matches_model_py_4822", "reproduces_model_py_cells")},
               "note": "step 139 (attack.py) ran on the all-nights frame and step 138 on the full diagnostic nights; the split nights carry "
                       "no sleep amount and fall in no cell, so both cuts' cell sizes are identical (asserted) and the per-row events are "
                       "asserted equal below. The attack's own flags matches_model_py_4822 (a v7 literal) and reproduces_model_py_cells "
                       "are false because they compare against the all-nights healthy count and a typed v7 value: a chain-script item, not a data difference."}
log("cohort note:", json.dumps(COHORT_NOTE))

# ------------------------------------------------------------------ the rows, as the round-30 builder formed them
res300 = res[res.analysis == "healthy_tst300"]
d = res300[~res300.negative_control].dropna(subset=["p_cell2_vs_cell1_hr"]).copy()
md_all = mdhr[mdhr.contrast == "cell2_vs_cell1"].copy()
md = md_all.set_index("key")
lm = lmvp.set_index("key")
d["mdhr80"] = d.key.map(md.mdhr80)
d["lm_hr"] = d.key.map(lm.cell2_vs_cell1_hr_L)
d["lm_lo"] = d.key.map(lm.cell2_vs_cell1_lo)
d["lm_hi"] = d.key.map(lm.cell2_vs_cell1_hi)
d = d.sort_values("p_cell2_vs_cell1_hr", ascending=False).reset_index(drop=True)
n = len(d)
assert n > 0 and d.mdhr80.notna().all(), "a drawn row has no smallest detectable hazard ratio"
# the per-row events of the power file equal the observed contrast's (the same cells)
for r in d.itertuples():
    assert int(md.loc[r.key, "ev_a"]) == int(r.ev_cell2) and int(md.loc[r.key, "ev_ref"]) == int(r.ev_cell1), (r.key, md.loc[r.key, ["ev_a", "ev_ref"]].tolist(), r.ev_cell2, r.ev_cell1)
    # the per-row at-risk cells of the power file are the observed contrast's own (n_cell2, n_cell1 exclude the prevalent cases of that outcome)
    assert int(md.loc[r.key, "n_a"]) == int(r.n_cell2) and int(md.loc[r.key, "n_ref"]) == int(r.n_cell1), (r.key, md.loc[r.key, ["n_a", "n_ref"]].tolist(), r.n_cell2, r.n_cell1)

# the 1-year landmark values in attack_landmark_vs_primary.csv and attack_landmark1y.csv are the same numbers
l1 = lm1y.set_index("key")
chk = [(k, float(l1.loc[k, "cell2_vs_cell1_hr"]), float(lm.loc[k, "cell2_vs_cell1_hr_L"]))
       for k in d.key if k in l1.index and np.isfinite(l1.loc[k, "cell2_vs_cell1_hr"])]
assert all(abs(a - b) < 5e-4 for _k, a, b in chk), [c for c in chk if abs(c[1] - c[2]) >= 5e-4]

pc = summ["summary"]["healthy_tst300"]["cell2_vs_cell1"]
lc = atk["attacks"]["A4_selection"]["landmark_1y"]["cell2_vs_cell1"]
lc_cmp = atk["attacks"]["A4_selection"]["landmark_1y"].get("cell2_vs_cell1_primary_for_comparison", {})
n_lm = int(np.isfinite(d.lm_hr.astype(float)).sum())
assert n_lm == lc["n"], (n_lm, lc["n"])
assert n == pc["n_conditions_tested"], (n, pc["n_conditions_tested"])
if lc_cmp:
    assert lc_cmp["n"] == pc["n_conditions_tested"] and lc_cmp["n_sig"] == pc["n_significant"], (lc_cmp, pc)   # the attack's own copy of the observed column agrees with step 138

MDH = md_all.mdhr80.astype(float)
assert np.isfinite(MDH).all(), "a non-finite smallest detectable hazard ratio among the contrast rows"
med = float(MDH.median())
N_ALL = int(len(MDH))
N_ALL_CONTROLS = int(md_all.negative_control.sum())
med_str = f"{med:.2f}"
log(f"rows {n} (V13 24), landmark rows {n_lm} (V13 17), contrast rows in the power file {N_ALL} (V13 50, of which {N_ALL_CONTROLS} negative controls), "
    f"median smallest detectable hazard ratio {med} prints {med_str} (V13 5.97)")
log(f"drawn diamonds: median {float(np.median(d.mdhr80)):.3f}, max {float(np.max(d.mdhr80)):.3f}, "
    f"dotted line at {med:.3f}, so the line sits right of {int((d.mdhr80 < med).sum())} of {n} diamonds")

# ------------------------------------------------------------------ design system (efig1213_polish_common, the round-30 part d)
W = 171.0 / 25.4
GUT = 0.18
AX_X0 = 1.92
ROW_IN = 0.30
TOPBAND = 1.06
RIGHT = W - 0.18
PAIR_AXIS = ("Hazard ratio, short sleep with low oxygen against\n"
             "short sleep with preserved nocturnal oxygenation (95% CI)")
KEY_COLORS_D = "blue = all follow-up, grey = 1-year landmark re-fit"
KEY_DIAMOND = "open diamond = the smallest hazard ratio the data could detect, larger = thinner data"
KEY_DOTTED = f"dotted line = median of those ratios over all {N_ALL} outcomes, not the {n} drawn ({med_str})"
for s in (PAIR_AXIS.split("\n")[0], PAIR_AXIS.split("\n")[1], KEY_COLORS_D, KEY_DIAMOND,
          "Observed contrast", "1-year landmark re-fit", "Smallest hazard ratio", "detectable", "at 80% power"):
    assert s in BASE_STRINGS, s            # the static strings of the V13 part d (the numbers in the heads and the dotted key are data)

EXPECTED = {"sheet": SHEET, "panel": "d", "sources": SOURCES, "cohort_note": COHORT_NOTE, "values": []}
SUMM_PATH, ATK_PATH, MD_PATH, RES_PATH = f"{SRC}/summary.json", f"{SRC}/attack.json", f"{SRC}/attack_min_detectable_hr.csv", f"{SRC}/results.csv"


def expect(text, source, key, value, rule="text"):
    EXPECTED["values"].append(dict(panel="d", text=text, source=source, key=key, value=value, rule=rule))


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


def bands(ax, m):
    for i in range(m):
        if i % 2 == 1:
            yy = m - 1 - i
            ax.axhspan(yy - 0.5, yy + 0.5, color=BAND, lw=0, zorder=0)


def refline(ax, x=1.0):
    ax.axvline(x, color=INK, lw=0.9, ls=(0, (4, 3)), zorder=1)


def left_rows(ax, ys, labels, sizes=None):
    ax.set_yticks(list(ys))
    ax.set_yticklabels(labels, fontsize=TCKF, color=INK, ha="left")
    if sizes is not None:   # R49: per-label size (a label kept at its current size when the gutter rule fails at +1 pt)
        for t, sz in zip(ax.get_yticklabels(), sizes):
            t.set_fontsize(sz)
    ax.tick_params(axis="y", pad=(AX_X0 - GUT) * 72.0, length=0)


def fit_gutter(fig, texts):
    """R49: per-label sizes. A label fits the gutter rule at TCKF, or keeps its current size (TCKF - PT_PLUS) when it does not
    (the brief: never re-wrap, never move a data mark, keep that one element at its current size). Returns the sizes."""
    sizes = []
    for s in texts:
        if GUT + measure(fig, s, TCKF) <= AX_X0 - 0.10:
            sizes.append(TCKF)
        else:
            assert GUT + measure(fig, s, TCKF - PT_PLUS) <= AX_X0 - 0.10, f"gutter label too wide even at the current size: {s!r}"
            sizes.append(TCKF - PT_PLUS); KEPT.append(s)
    return sizes


def panel_head(ax, title, size=ANNF):
    ax.annotate(title, xy=(0.0, 1.0), xycoords="axes fraction", xytext=(0, 6), textcoords="offset points",
                ha="left", va="bottom", fontsize=size, color=INK, linespacing=1.18)


def ci_row(ax, lo, hi, hr, y, col, xlo, xhi, marker="o", size=S_CIRCLE):
    stop, start = xhi / 1.22, xlo * 1.22
    tip, ltip = xhi / 1.03, xlo * 1.03
    lo_d, hi_d = max(lo, start), min(hi, stop)
    ax.plot([lo_d, hi_d], [y, y], color=col, lw=CI_LW, solid_capstyle="round", zorder=3)
    # R49: the arrowhead is a data mark. annotate sizes it by the text size (mutation_scale defaults to font.size, now TCKF = 11), so it is
    # pinned to the V26 value (10.0 = the round-40 TCKF): heads stay 4.0 pt, identical to V26
    if hi > stop:
        ax.annotate("", xy=(tip, y), xytext=(stop, y), arrowprops=dict(arrowstyle="-|>", color=col, lw=1.2, shrinkA=0, shrinkB=0, mutation_scale=TCKF - PT_PLUS))
    if lo < start:
        ax.annotate("", xy=(ltip, y), xytext=(start, y), arrowprops=dict(arrowstyle="-|>", color=col, lw=1.2, shrinkA=0, shrinkB=0, mutation_scale=TCKF - PT_PLUS))
    if start <= hr <= stop:
        ax.scatter([hr], [y], s=size, marker=marker, color=col, zorder=5, edgecolors="white", linewidths=MEDGE)
    return dict(clipped_hi=bool(hi > stop), clipped_lo=bool(lo < start), marker_drawn=bool(start <= hr <= stop))


# ================================================================== part d
XLO = 0.25
AXIS_LADDER = [(10.8, [0.5, 1, 2, 4, 8]), (24.0, [0.5, 1, 2, 4, 8, 16]), (48.0, [0.5, 1, 2, 4, 8, 16, 32])]
need = max(med, float(d.mdhr80.max()))                     # the dotted median and every diamond must sit inside the axis
XHI, TICKS = next((x, t) for x, t in AXIS_LADDER if need <= x)
AXIS_NOTE = {"v13": [0.25, 10.8], "new": [XLO, XHI], "ticks_v13": [0.5, 1, 2, 4, 8], "ticks_new": TICKS,
             "why": f"the dotted median {med:.3f} and the largest drawn diamond {float(d.mdhr80.max()):.3f} must sit inside the axis; the V13 range ends at 10.8"}
log("axis:", json.dumps(AXIS_NOTE))
PW, GAPX = (RIGHT - AX_X0 - 2 * 0.30) / 3.0, 0.30
AX_H = n * ROW_IN
KEYROW = 0.24
BOT = 0.94            # R40: no key rows under the axis title (V16: 0.94 + 3 * KEYROW); the top-anchored content is unchanged
Hd = BOT + AX_H + TOPBAND
log(f"part d: {n} rows (V13 24), height {Hd * 72:.2f} pt (V13 714.24), pitch kept at {ROW_IN} in")

c_obs, c_all, c_fdr = pc["n_significant"], pc["n_conditions_tested"], pc["n_significant_fdr"]
l_sig, l_all, l_fdr = lc["n_sig"], lc["n"], lc["n_fdr"]
HEADS = [f"Observed contrast\n{c_obs} of {c_all} significant,\n{c_fdr} after FDR",
         f"1-year landmark re-fit\n{l_sig} of {l_all} significant,\n{l_fdr} after FDR",
         "Smallest hazard ratio\ndetectable\nat 80% power"]
expect(f"{c_obs} of {c_all} significant,", SUMM_PATH, "summary/healthy_tst300/cell2_vs_cell1", {"n_significant": c_obs, "n_conditions_tested": c_all}, rule="sig_of_tested")
expect(f"{c_fdr} after FDR", SUMM_PATH, "summary/healthy_tst300/cell2_vs_cell1", {"n_significant_fdr": c_fdr}, rule="fdr_after")
expect(f"{l_sig} of {l_all} significant,", ATK_PATH, "attacks/A4_selection/landmark_1y/cell2_vs_cell1", {"n_sig": l_sig, "n": l_all}, rule="lm_sig_of_n")
expect(f"{l_fdr} after FDR", ATK_PATH, "attacks/A4_selection/landmark_1y/cell2_vs_cell1", {"n_fdr": l_fdr}, rule="lm_fdr_after")
# R40: the dotted-line key no longer prints (its median, counts and the drawn-row count go to the figure legend, see EXPECTED["r40_keys_removed_d"])
for dis in d.disease:
    expect(dis, RES_PATH, "row set: healthy_tst300, non-control, the short-with-low against short-with-preserved contrast fitted, sorted by its hazard ratio", dis)
for t in TICKS:
    for _col in range(3):
        expect(f"{t:g}", "static:axis", f"x tick of part d (axis {XLO} to {XHI}, data-driven range)", f"{t:g}")


def build_d():
    y = np.arange(n, dtype=float)[::-1]
    EXPECTED["panelD"] = {"n_rows": n, "n_landmark": n_lm, "median_mdhr80_all": med, "n_all_contrast_rows": N_ALL, "n_all_controls": N_ALL_CONTROLS,
                          "primary_counts": {"n_significant": c_obs, "n_conditions_tested": c_all, "n_significant_fdr": c_fdr},
                          "landmark_counts": {"n_sig": l_sig, "n": l_all, "n_fdr": l_fdr,
                                              "n_retained": atk["attacks"]["A4_selection"]["landmark_1y"]["n_retained"]},
                          "rows": [], "row_order": list(d.disease), "axis": AXIS_NOTE,
                          "geometry": {"row_in": ROW_IN, "keyrow": KEYROW, "bot": BOT, "topband": TOPBAND,
                                       "height_pt": Hd * 72.0, "v13_height_pt": 714.24, "xlim": [XLO, XHI], "ticks": TICKS}}
    with plt.rc_context(rc_polish()):
        fig = plt.figure(figsize=(W, Hd))
        fig.canvas.draw()
        sizes_d = fit_gutter(fig, list(d.disease))   # R49
        axes = [fig.add_axes(frac(fig, AX_X0 + j * (PW + GAPX), BOT, PW, AX_H)) for j in range(3)]
        for ax in axes:
            bands(ax, n)
            logx(ax, TICKS)
            ax.set_xlim(XLO, XHI)
            refline(ax)
            ax.set_yticks(y)
            ax.set_yticklabels([])
            ax.tick_params(axis="y", length=0)
            bare_y(ax)
            ax.set_ylim(-0.7, n - 0.3)
        left_rows(axes[0], y, list(d.disease), sizes_d)
        for i in range(n):
            r = d.iloc[i]
            co = ci_row(axes[0], r.p_cell2_vs_cell1_lo, r.p_cell2_vs_cell1_hi, r.p_cell2_vs_cell1_hr, y[i], BLUE, XLO, XHI)
            cl = None
            if np.isfinite(r.lm_hr):
                cl = ci_row(axes[1], r.lm_lo, r.lm_hi, r.lm_hr, y[i], GREY, XLO, XHI)
            axes[2].scatter([r.mdhr80], [y[i]], s=34, marker="D", facecolors="white",
                            edgecolors=BLUE_MID, linewidths=1.2, zorder=5)
            EXPECTED["panelD"]["rows"].append(dict(
                disease=r.disease, key=r.key, ev_cell1=int(r.ev_cell1), ev_cell2=int(r.ev_cell2),
                hr=r.p_cell2_vs_cell1_hr, lo=r.p_cell2_vs_cell1_lo, hi=r.p_cell2_vs_cell1_hi, q=r.p_cell2_vs_cell1_q,
                lm_hr=r.lm_hr, lm_lo=r.lm_lo, lm_hi=r.lm_hi, mdhr80=r.mdhr80,
                observed=co, landmark=cl, diamond_drawn=bool(XLO <= r.mdhr80 <= XHI)))
        assert XLO <= med <= XHI, med
        axes[2].axvline(med, color=BLUE_MID, lw=1.0, ls=(0, (1.6, 1.8)), zorder=1)
        for ax, ttl in zip(axes, HEADS):
            panel_head(ax, ttl)
        # R40 (Alen, 2026-09-18): the three printed keys (colours, open diamond, dotted line) leave the sheet for the figure legend
        fig.text((AX_X0 + RIGHT) / 2.0 / W, 0.20 / Hd, PAIR_AXIS, ha="center", va="bottom",
                 fontsize=LABF, color=INK, linespacing=1.35)
        fig.canvas.draw()
        EXPECTED["r40_keys_removed_d"] = {"strings": [KEY_COLORS_D, KEY_DIAMOND, KEY_DOTTED], "median_mdhr80_all": med, "n_all": N_ALL, "n_drawn": n}
        text_gate(fig)
        out = f"{WORK}/part_d.pdf"
        fig.savefig(out)
        plt.close(fig)
        log("wrote", out, "size", [round(v * 72, 2) for v in fig.get_size_inches()])
        return out


build_d()
json.dump(EXPECTED, open(f"{WORK}/expected_values_d.json", "w"), indent=1, default=float)
log("wrote", f"{WORK}/expected_values_d.json", "with", len(EXPECTED["values"]), "expected printed values")

# ------------------------------------------------------------------ merge into the sheet record (kept_from_v13 emptied: nothing is kept any more)
EV = f"{WORK}/expected_values.json"
if not os.path.exists(f"{WORK}/expected_values_PRE_V8_2.json"):
    shutil.copy2(EV, f"{WORK}/expected_values_PRE_V8_2.json")
E = json.load(open(EV))
E["values"] = [v for v in E["values"] if v["panel"] != "d"] + EXPECTED["values"]
E["sources"].update(SOURCES)
E["panelD"] = EXPECTED["panelD"]
E["cohort_note_d"] = COHORT_NOTE
E["kept_from_v13"] = {}
E["part_d_v8_2"] = {"built_by": "01_build_d.py (v8.2)", "step_139_sidecars_after": V8_2_CUT, "part_file": f"{WORK}/part_d.pdf"}
E["r49_kept"]["gutter_labels_d"] = sorted(set(KEPT))
log("R49 part d labels kept at the current size:", sorted(set(KEPT)))
json.dump(E, open(EV, "w"), indent=1, default=float)
log("merged the part-d record into", EV, "(kept_from_v13 emptied; pre-merge copy expected_values_PRE_V8_2.json)")
LOG.close()
