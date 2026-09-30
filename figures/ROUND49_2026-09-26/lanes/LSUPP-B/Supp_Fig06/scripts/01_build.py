"""Supp Fig 6, lane V14_L3b_DURATION_SUPP, round 37 (2026-09-15): re-plot parts a, b, c at the V13 sheet's geometry from the
v8.1 step-138 outputs (Sleep_Variability_2026-08/dur_x_oxygen_healthy/, run on the 15,551 full diagnostic nights,
T90_FULL_NIGHTS_ONLY=1). Repointed copy of the round-30 builder (01_build_PRE_V8_1.py): the lane paths, the v8.1 sidecar
gate, the v8 frozen table through cohort_spec under the full-nights rule, the base text from the V13 sheet itself.
Part d (smallest detectable hazard ratio, 1-year landmark re-fit) is built by 01_build_d.py from the v8.2 step-139 outputs
(re-extracted 2026-09-15 19:30, the owner's v8.2 rule: no panel keeps old data); run 01_build_d.py after this script, it merges
its record into work/expected_values.json and empties kept_from_v13.

Geometry is the polished builders' (FF/split_efig13_polish.py for a, FF/redesign_eFigure13_supp_polish.py with
efig1213_polish_common for b and c): 171 mm parts, GUT 0.18 in, AX_X0 1.92 in, ROW_IN 0.30 in (b), PITCH 0.21 in (c),
TOPBAND 1.06 in, the same key rows. Strings are the V13 sheet's. Colours from the V13 drawing layer (#0288d1 low oxygen,
#8a9099 grey, #3f9fd8 normal-duration low-oxygen group, #eef0f1 bands and preserved squares).

Panel a's death rates per 1,000 person-years are formed exactly as the polished builder formed them: the v8 frozen table
through numbers/cohort_spec.apply_cohort (T90_FULL_NIGHTS_ONLY=1, the cohort step 138 ran on), the organ-disease-free
filter of dur_x_oxygen_healthy/model.py (ORGAN_GROUP blocks plus stroke and cancer) and its cell rule, with n and deaths per
cell asserted equal to summary.json before any rate is printed. GAP-9 note: the round-30 base sheet's panel a was drawn by
this same table-derived route (not from the frozen eFigure13_drawn_values.json of NF/split_efig13.py), so no GAP-9 artefact
feeds V14 either.

Writes work/part_a.pdf, part_b.pdf, part_c.pdf, work/base_text.json, work/expected_values.json, work/build_log.txt."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import json
import os
import sys

os.environ["T90_FULL_NIGHTS_ONLY"] = "1"      # step 138 ran under this rule: cohort_spec then returns the 15,551 full nights
import numpy as np
import pandas as pd

LANE = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LSUPP-B"   # R49: SD = LANE/Supp_Fig06 = this sheet folder
sys.path.insert(0, f"{LANE}/Supp_Fig06/scripts")
from l3b_common import (SV, T90, V8TAB, BASE, INK, LABF, TCKF, ANNF, PT_PLUS, CI_LW, MEDGE, S_CIRCLE, rc_polish, measure,  # noqa: E402
                        text_gate, sidecar, jload, hydrated, spans, plt)
from matplotlib import ticker as mticker  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

SHEET = "Supp_Fig06"
SD = f"{LANE}/{SHEET}"
WORK = f"{SD}/work"
os.makedirs(WORK, exist_ok=True)
LOG = open(f"{WORK}/build_log.txt", "w")


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    LOG.write(s + "\n")


BLUE, GREY, BLUE_MID, BAND = "#0288d1", "#8a9099", "#3f9fd8", "#eef0f1"
STARF = 10.0 + PT_PLUS   # R49
KEPT = []   # R49: strings kept at their current size (a placement rule fails at +1 pt), reported
KEPT_CELL = []
SRC = f"{SV}/dur_x_oxygen_healthy"
SOURCES = {}
for f in ("results.csv", "cells.csv", "comparison.csv", "summary.json"):
    sh, sc = sidecar(f"{SRC}/{f}")
    assert sc["step"]["id"] == "138" and sc["step"].get("rc", 0) == 0, sc["step"]
    SOURCES[f] = {"path": f"{SRC}/{f}", "sha256": sh, "step": sc["step"], "mtime": sc["output"]["mtime_local"]}
    log(f, "sha256", sh, "step 138 rc 0 stamped", sc["output"]["mtime_local"])
res = pd.read_csv(hydrated(f"{SRC}/results.csv"))
comp = pd.read_csv(hydrated(f"{SRC}/comparison.csv"))
summ = jload(f"{SRC}/summary.json")
HD = summ["healthy_definition"]
CCS = summ["cell_counts_side_by_side"]
# R40: the V13 text record copied from the round-37 lane (work/base_text.json), no PyMuPDF probe of a composed sheet
BASE_SP = json.load(open(hydrated(f"{WORK}/base_text.json"))); PAGE_W, PAGE_H = 1021.5, 1420.0
BASE_STRINGS = {s["text"] for s in BASE_SP}
assert (round(PAGE_W, 2), round(PAGE_H, 2)) == (1021.5, 1420.0), (PAGE_W, PAGE_H)
EXPECTED = {"sheet": SHEET, "sources": SOURCES, "values": [], "kept_from_v13": {}}


def expect(panel, text, source, key, value, rule="text"):
    EXPECTED["values"].append(dict(panel=panel, text=text, source=source, key=key, value=value, rule=rule))


# ------------------------------------------------------------------ design system (efig1213_polish_common)
W = 171.0 / 25.4
EDGE = GUT = 0.18
AX_X0 = 1.92
ROW_IN = 0.30
TOPBAND = 1.06
RIGHT = W - 0.18
RIGHT_STAR = W - 0.40
REF_AXIS = "Hazard ratio against normal-duration sleep\nwith preserved nocturnal oxygenation (95% CI)"
KEY_STARS = "* q < 0.05, ** q < 0.01, *** q < 0.001"
KEY_ARROW = "arrow = interval continues beyond the axis"
KEY_COLORS_B = "grey = preserved-oxygenation group, blue = low-oxygen groups"
KEY_COLORS_C = ("grey = full cohort", "blue = healthy subgroup")
for s in (KEY_STARS, KEY_ARROW, KEY_COLORS_B, *KEY_COLORS_C, "T90 ≤1%", "Low oxygen", "T90 above 10%", "4-hour cut",
          "5-hour cut, primary", "Reference group", "Short sleep", "Normal-duration sleep", "6 h or more"):
    assert s in BASE_STRINGS, s


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


def q_stars(q):
    if q is None or not np.isfinite(q):
        return ""
    return "***" if q < .001 else ("**" if q < .01 else ("*" if q < .05 else ""))


def finish(fig, name):
    text_gate(fig)
    out = f"{WORK}/part_{name}.pdf"
    fig.savefig(out)
    plt.close(fig)
    log("wrote", out, "size", [round(v * 72, 2) for v in fig.get_size_inches()])
    return out


# ================================================================== part a: the two 2 x 2 designs
sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
import cohort_spec  # noqa: E402
from cohort_spec import apply_cohort  # noqa: E402
from disease_definitions import DISEASES, ORGAN_GROUP  # noqa: E402

assert cohort_spec.FULL_NIGHTS_ONLY and cohort_spec.COHORT_N == HD["n_cohort"], (cohort_spec.COHORT_N, HD["n_cohort"])
GROUPS = ["cell1", "cell2", "cell3", "cell4"]
T90_LOW, T90_PRESERVED, TST_NORMAL = 10.0, 1.0, 360.0      # the model's definitions (dur_x_oxygen_healthy/model.py), proved below by the cell counts
V8 = f"{V8TAB}/t90_final.parquet"
sidecar_ok = os.path.exists(V8 + ".provenance.json") or os.path.exists(f"{V8TAB}/sha256.txt")
assert sidecar_ok, "v8 table without provenance"
_b = apply_cohort(pd.read_parquet(hydrated(V8)))
assert len(_b) == HD["n_cohort"], (len(_b), HD["n_cohort"])
# the healthy filter of dur_x_oxygen_healthy/model.py: ORGAN_GROUP blocks plus stroke and cancer (its two typed extras)
ORGAN_BLOCKS = [("cardiovascular", [k for k, g in ORGAN_GROUP.items() if g == "Cardiac"] + ["stroke_any", "cvd"]),
                ("pulmonary", [k for k, g in ORGAN_GROUP.items() if g == "Respiratory"]),
                ("renal", [k for k, g in ORGAN_GROUP.items() if g == "Kidney"]),
                ("hepatic", [k for k, g in ORGAN_GROUP.items() if g == "Liver"]),
                ("metabolic", [k for k, g in ORGAN_GROUP.items() if g == "Metabolic"]),
                ("cancer", ["cancer_any", "prostate_ca"])]
# the file's own record of the filter must be this list, label for label
for blk, ks in ORGAN_BLOCKS:
    assert HD["excluded_on"][blk] == [DISEASES[k][0] for k in ks], (blk, HD["excluded_on"][blk], [DISEASES[k][0] for k in ks])
assert HD["n_conditions_excluded_on"] == sum(len(ks) for _b_, ks in ORGAN_BLOCKS)
_keep = np.ones(len(_b), bool)
for _blk, _ks in ORGAN_BLOCKS:
    for _k in _ks:
        _keep &= ~(_b[f"{_k}_prevalent"] == 1).values
H = _b[_keep].copy()
assert len(H) == HD["n_healthy"], (len(H), HD["n_healthy"])
assert H.fu_valid.eq(1).all() and H.death_years.notna().all()
log("healthy subgroup reproduced from the v8 table under the full-nights rule:", len(H), "of", len(_b))


def assign(d, tst_short):
    low = d.spo2_pct_below_90 > T90_LOW
    pres = d.spo2_pct_below_90 <= T90_PRESERVED
    shrt = d.TST_min < tst_short
    norm = d.TST_min >= TST_NORMAL
    return pd.Series(np.select([shrt & pres, shrt & low, norm & low, norm & pres], GROUPS, default="excluded"), index=d.index)


RATE, PY, NBY, DBY = {}, {}, {}, {}
for key in ("healthy_tst240", "healthy_tst300"):
    thr = float(CCS[key]["tst_short_lt_min"])
    c = assign(H, thr)
    n = {g: int((c == g).sum()) for g in GROUPS}
    d = {g: int(H.loc[c == g, "death_incident"].sum()) for g in GROUPS}
    assert n == CCS[key]["n_by_cell"], (key, n, CCS[key]["n_by_cell"])
    assert d == CCS[key]["deaths_by_cell"], (key, d, CCS[key]["deaths_by_cell"])
    PY[key] = {g: float(H.loc[c == g, "death_years"].sum()) for g in GROUPS}
    RATE[key] = {g: 1000.0 * d[g] / PY[key][g] for g in GROUPS}
    NBY[key], DBY[key] = n, d
    log(key, "cut", thr, "n", n, "deaths", d, "rates", {g: round(RATE[key][g], 2) for g in GROUPS})
EXPECTED["panelA"] = {k: dict(n=NBY[k], deaths=DBY[k], person_years=PY[k], rate_per_1000=RATE[k], tst_short_lt_min=float(CCS[k]["tst_short_lt_min"])) for k in NBY}
EXPECTED["thresholds"] = dict(T90_LOW=T90_LOW, T90_PRESERVED=T90_PRESERVED, TST_NORMAL=TST_NORMAL, proved_by="cell n and deaths equal summary.json for both cuts")

DIAGRAMS = [("healthy_tst240", "4-hour cut"), ("healthy_tst300", "5-hour cut, primary")]
GUT_L = 1.56
AX_W = W - GUT_L - 0.30
AX_H = 1.86
LET_H = 0.46
LAB_DY = 0.06
BELOW = 0.46
INTER = 0.20
TOP_IN = 0.16
BOT_PAD = 0.22
Hc = TOP_IN + 2 * (LET_H + AX_H + BELOW) + INTER + BOT_PAD
assert abs(Hc * 72 - 442.08) < 0.05, Hc * 72
SUMM_PATH = f"{SRC}/summary.json"


def build_a():
    with plt.rc_context(rc_polish()):
        fig = plt.figure(figsize=(W, Hc))
        fig.canvas.draw()
        for i, (key, cut_lab) in enumerate(DIAGRAMS):
            y_top = TOP_IN + i * (LET_H + AX_H + BELOW + INTER) + LET_H
            ax = fig.add_axes([GUT_L / W, 1.0 - (y_top + AX_H) / Hc, AX_W / W, AX_H / Hc])
            ax.set_xlim(0, 2)
            ax.set_ylim(0, 2)
            ax.axis("off")
            for g, cx, cy in (("cell1", 0, 1), ("cell2", 1, 1), ("cell4", 0, 0), ("cell3", 1, 0)):
                low_oxy = cx == 1
                ax.add_patch(Rectangle((cx + 0.02, cy + 0.02), 0.96, 0.96, facecolor=BLUE if low_oxy else BAND,
                                       edgecolor=BLUE if low_oxy else GREY, lw=0.8))
                tc = "white" if low_oxy else INK
                n_txt, d_txt = f"n = {NBY[key][g]:,}", f"{DBY[key][g]} deaths"
                r_txt = f"{RATE[key][g]:.1f} deaths per 1,000 person-years"
                expect("a", n_txt, SUMM_PATH, f"cell_counts_side_by_side/{key}/n_by_cell/{g}", NBY[key][g], rule="n_eq_int_comma")
                expect("a", d_txt, SUMM_PATH, f"cell_counts_side_by_side/{key}/deaths_by_cell/{g}", DBY[key][g], rule="deaths_int")
                expect("a", r_txt, V8, f"1000 x deaths / sum(death_years) in {key} {g} (v8 table through cohort_spec, full nights, healthy filter)", RATE[key][g], rule="rate_1dp")
                lines = [n_txt, d_txt, r_txt]
                if g == "cell4":
                    lines = ["Reference group"] + lines
                step, start = 0.175, 0.5 + 0.175 * (len(lines) - 1) / 2.0
                for j, ln in enumerate(lines):
                    bold = (j == 0 and g == "cell4")
                    # R49: a cell line takes ANNF when it fits the cell rule, else the largest of ANNF - 0.5 and ANNF - PT_PLUS that fits (never smaller than now)
                    sz = next(s for s in (ANNF, ANNF - 0.5, ANNF - PT_PLUS) if measure(fig, ln, s, "bold" if bold else "normal") <= 0.96 * AX_W / 2.0 - 0.08)
                    if sz != ANNF: KEPT_CELL.append((ln, sz))
                    ax.text(cx + 0.5, cy + start - j * step, ln, ha="center", va="center", fontsize=sz, color=tc,
                            fontweight="bold" if bold else "normal")
                    assert measure(fig, ln, sz, "bold" if bold else "normal") <= 0.96 * AX_W / 2.0 - 0.08, ln
            for _y, _lab in ((1.5, "Short sleep"), (0.5, "Normal-duration sleep\n6 h or more")):
                ax.text(-0.03, _y, _lab, ha="right", va="center", fontsize=ANNF, color=INK, linespacing=1.2)
            ax.text(0.5, -0.06, "T90 ≤1%", ha="center", va="top", fontsize=ANNF, color=INK, linespacing=1.2)
            ax.text(1.5, -0.06, "Low oxygen\nT90 above 10%", ha="center", va="top", fontsize=ANNF, color=INK, linespacing=1.2)
            _lx = -(GUT_L - GUT) / AX_W
            ax.text(_lx, 1.0 + LAB_DY / AX_H, cut_lab, transform=ax.transAxes, fontsize=ANNF, color=INK, ha="left", va="bottom")
        return finish(fig, "a")


# ================================================================== part b: each group against the reference
CELLS = [("p_cell1_vs_cell4", "cell1", "Short sleep,\npreserved oxygenation", GREY),
         ("p_cell2_vs_cell4", "cell2", "Short sleep,\nlow oxygen", BLUE),
         ("p_cell3_vs_cell4", "cell3", "Normal-duration sleep,\nlow oxygen", BLUE_MID)]
res300 = res[res.analysis == "healthy_tst300"].copy()
RES_PATH = f"{SRC}/results.csv"


def build_b():
    b = res300[~res300.negative_control].dropna(subset=["p_cell1_vs_cell4_hr", "p_cell2_vs_cell4_hr", "p_cell3_vs_cell4_hr"]).copy()
    b = b.sort_values("p_cell2_vs_cell4_hr", ascending=False).reset_index(drop=True)
    n = len(b)
    y = np.arange(n, dtype=float)[::-1]
    XLO, XHI = 0.22, 18.0
    ticks = [0.25, 1, 2, 4, 8, 16]
    ns = CCS["healthy_tst300"]["n_by_cell"]
    PW, GAPX = (RIGHT_STAR - AX_X0 - 2 * 0.42) / 3.0, 0.42
    AX_H_b = n * ROW_IN
    KEYROW = 0.24
    BOT = 0.94            # R40: no key rows under the axis title (V16: 0.94 + 2 * KEYROW); the top-anchored content is unchanged
    Hb = BOT + AX_H_b + TOPBAND
    n_base = sum(1 for s in BASE_SP if s["size"] == 10.0 and s["font"] == "ArialMT" and abs(s["x"] - 535.68) < 0.5 and 40 < s["y"] < 530)
    log(f"part b: {n} rows (V13 {n_base}), height {Hb * 72:.2f} pt")
    EXPECTED["panelB"] = {"n_rows": n, "n_rows_v13": n_base, "rows": [], "stars": {}, "height_pt": Hb * 72}
    with plt.rc_context(rc_polish()):
        fig = plt.figure(figsize=(W, Hb))
        fig.canvas.draw()
        sizes_b = fit_gutter(fig, list(b.disease))   # R49
        for j, (pre, key, lab, col) in enumerate(CELLS):
            ax = fig.add_axes(frac(fig, AX_X0 + j * (PW + GAPX), BOT, PW, AX_H_b))
            bands(ax, n)
            logx(ax, ticks)
            ax.set_xlim(XLO, XHI)
            refline(ax)
            n_star = 0
            for i in range(n):
                r = b.iloc[i]
                clip = ci_row(ax, r[f"{pre}_lo"], r[f"{pre}_hi"], r[f"{pre}_hr"], y[i], col, XLO, XHI)
                st = q_stars(r[f"{pre}_q"])
                if st:
                    n_star += 1
                    ax.text(1.03, y[i], st, transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=STARF, color=INK)
                    expect("b", st, RES_PATH, f"analysis=healthy_tst300&disease={r.disease}:{pre}_q", float(r[f"{pre}_q"]), rule="stars_q")
                EXPECTED["panelB"]["rows"].append(dict(disease=r.disease, key=r.key, contrast=pre, hr=float(r[f"{pre}_hr"]), lo=float(r[f"{pre}_lo"]),
                                                       hi=float(r[f"{pre}_hi"]), q=float(r[f"{pre}_q"]), stars=st, **clip))
            EXPECTED["panelB"]["stars"][key] = n_star
            if j == 0:
                left_rows(ax, y, list(b.disease), sizes_b)
            else:
                ax.set_yticks(y)
                ax.set_yticklabels([])
                ax.tick_params(axis="y", length=0)
            bare_y(ax)
            ax.set_ylim(-0.7, n - 0.3)
            panel_head(ax, f"{lab}\nn = {ns[key]:,}")
            expect("b", f"n = {ns[key]:,}", SUMM_PATH, f"cell_counts_side_by_side/healthy_tst300/n_by_cell/{key}", ns[key], rule="n_eq_int_comma")
        for dis in b.disease:
            expect("b", dis, RES_PATH, "row set: healthy_tst300, non-control, the three pairwise contrasts fitted, sorted by the short-sleep low-oxygen hazard ratio", dis)
        # R40 (Alen, 2026-09-18): the three printed keys (stars, arrow, colours) leave the sheet for the figure legend; the stars stay on the plot
        axt = fig.text((AX_X0 + RIGHT_STAR) / 2.0 / W, 0.20 / Hb, REF_AXIS, ha="center", va="bottom", fontsize=LABF, color=INK, linespacing=1.35)
        fig.canvas.draw()
        EXPECTED["r40_keys_removed_b"] = [KEY_STARS, KEY_ARROW, KEY_COLORS_B]
        EXPECTED["panelB"]["row_order"] = list(b.disease)
        log("part b stars", EXPECTED["panelB"]["stars"], "rows", list(b.disease))
        return finish(fig, "b")


# ================================================================== part c: healthy against full cohort
COMP_PATH = f"{SRC}/comparison.csv"


def build_c():
    comp300 = comp[(comp.healthy_analysis == "healthy_tst300") & (~comp.negative_control)].copy()
    c2 = comp300[comp300.contrast == "cell2_vs_cell4"].dropna(subset=["ratio_healthy_over_full"])
    c3all = comp300[comp300.contrast == "cell3_vs_cell4"].dropna(subset=["ratio_healthy_over_full"])
    c3 = c3all[~c3all.key.isin(set(c2.key))]
    left = c2.sort_values("full_cohort_hr", ascending=False).reset_index(drop=True)
    right = c3.sort_values("full_cohort_hr", ascending=False).reset_index(drop=True)
    XLO, XHI = 0.22, 24.0
    ticks = [0.25, 1, 2, 4, 8, 16]
    n_full, n_heal = HD["n_cohort"], HD["n_healthy"]
    PW, GAPX = (RIGHT - AX_X0 - 0.42) / 2.0, 0.42
    PITCH = 0.21
    TOP_H, LOW_H = len(left) * PITCH, len(right) * PITCH
    BOT, INTER_C = 0.94, 1.14
    # v8.1: every condition with the normal-with-low contrast is also in the short-with-low block, so the lower block has NO row.
    # An empty block cannot be drawn (a zero-height axes has a singular transform): the lower pair is dropped and the part ends
    # under the top block; the two lower panel heads leave the sheet (declared word delta, flagged for the owner).
    RIGHT_BLOCK = len(right) > 0
    LOW_Y0 = BOT
    TOP_Y0 = LOW_Y0 + (LOW_H + INTER_C if RIGHT_BLOCK else 0.0)
    Hc3 = TOP_Y0 + TOP_H + TOPBAND
    log(f"part c: {len(left)} left rows (V13 25), {len(right)} right rows (V13 2), height {Hc3 * 72:.2f} pt" + ("" if RIGHT_BLOCK else " (lower block dropped: no row)"))
    EXPECTED["panelC"] = {"n_left": len(left), "n_right": len(right), "left": [], "right": [], "n_full": n_full, "n_healthy": n_heal, "height_pt": Hc3 * 72,
                          "right_block_drawn": RIGHT_BLOCK, "lower_block_note": "" if RIGHT_BLOCK else "no condition has the normal-duration-with-low-oxygen contrast outside the short-with-low block under v8.1: the lower pair of panels (Full cohort / Healthy subgroup, Normal-duration sleep with low oxygen) is not drawn"}
    with plt.rc_context(rc_polish()):
        fig = plt.figure(figsize=(W, Hc3))
        fig.canvas.draw()
        sizes_left = fit_gutter(fig, list(left.disease)); sizes_right = fit_gutter(fig, list(right.disease))   # R49
        panels = [(TOP_Y0, TOP_H, 0, left, "full", f"Full cohort, n = {n_full:,}\nShort sleep\nwith low oxygen", "left"),
                  (TOP_Y0, TOP_H, 1, left, "heal", f"Healthy subgroup, n = {n_heal:,}\nShort sleep\nwith low oxygen", "left")]
        if RIGHT_BLOCK:
            panels += [(LOW_Y0, LOW_H, 0, right, "full", "Full cohort\nNormal-duration sleep\nwith low oxygen", "right"),
                       (LOW_Y0, LOW_H, 1, right, "heal", "Healthy subgroup\nNormal-duration sleep\nwith low oxygen", "right")]
        for y0, hh, gj, dfr, side, ttl, blk in panels:
            ax = fig.add_axes(frac(fig, AX_X0 + gj * (PW + GAPX), y0, PW, hh))
            m = len(dfr)
            yy = np.arange(m, dtype=float)[::-1]
            bands(ax, m)
            logx(ax, ticks)
            ax.set_xlim(XLO, XHI)
            refline(ax)
            col = GREY if side == "full" else BLUE
            pre = "full_cohort" if side == "full" else "healthy"
            for i in range(m):
                r = dfr.iloc[i]
                clip = ci_row(ax, r[f"{pre}_lo"], r[f"{pre}_hi"], r[f"{pre}_hr"], yy[i], col, XLO, XHI)
                EXPECTED["panelC"][blk].append(dict(disease=r.disease, key=r.key, side=side, hr=float(r[f"{pre}_hr"]), lo=float(r[f"{pre}_lo"]), hi=float(r[f"{pre}_hi"]), **clip))
            if gj == 0:
                left_rows(ax, yy, list(dfr.disease), sizes_left if blk == "left" else sizes_right)
                for dis in dfr.disease:
                    expect("c", dis, COMP_PATH, f"row set {blk} block: healthy_tst300, non-control, contrast {'cell2_vs_cell4' if blk == 'left' else 'cell3_vs_cell4 not in the left block'}, sorted by full_cohort_hr", dis)
            else:
                ax.set_yticks(yy)
                ax.set_yticklabels([])
                ax.tick_params(axis="y", length=0)
            bare_y(ax)
            ax.set_ylim(-0.7, m - 0.3)
            panel_head(ax, ttl)
        expect("c", f"Full cohort, n = {n_full:,}", SUMM_PATH, "healthy_definition/n_cohort", n_full, rule="n_eq_int_comma_prefix")
        expect("c", f"Healthy subgroup, n = {n_heal:,}", SUMM_PATH, "healthy_definition/n_healthy", n_heal, rule="n_eq_int_comma_prefix")
        # R40: the two colour keys (grey = full cohort, blue = healthy subgroup) leave for the figure legend; the part height is unchanged
        axt = fig.text((AX_X0 + RIGHT) / 2.0 / W, 0.20 / Hc3, REF_AXIS, ha="center", va="bottom", fontsize=LABF, color=INK, linespacing=1.35)
        fig.canvas.draw()
        EXPECTED["r40_keys_removed_c"] = list(KEY_COLORS_C)
        EXPECTED["panelC"]["left_order"] = list(left.disease)
        EXPECTED["panelC"]["right_order"] = list(right.disease)
        return finish(fig, "c")


# ================================================================== part d: built by 01_build_d.py (v8.2); nothing is kept from V13
EXPECTED["kept_from_v13"] = {}
log("part d: built by 01_build_d.py from the v8.2 step-139 files (run it after this script); nothing kept from V13")

build_a()
build_b()
build_c()
EXPECTED["place_v13"] = {"a": (14.0, 14.0), "b": (522.72, 14.0), "c": (14.0, 691.76), "d": (522.72, 691.76)}
EXPECTED["r49_kept"] = {"gutter_labels_abc": sorted(set(KEPT)), "cell_lines": KEPT_CELL, "sizes": dict(pt_plus=PT_PLUS, LABF=LABF, TCKF=TCKF, ANNF=ANNF, STARF=STARF, letters=13.0 + PT_PLUS)}
log("R49 kept at the current size: gutter labels", sorted(set(KEPT)), "cell lines", KEPT_CELL)
EXPECTED["page_v13"] = (1021.50, 1420.0)
json.dump(EXPECTED, open(f"{WORK}/expected_values.json", "w"), indent=1, default=float)
log("wrote", f"{WORK}/expected_values.json", "with", len(EXPECTED["values"]), "expected printed values")
LOG.close()
