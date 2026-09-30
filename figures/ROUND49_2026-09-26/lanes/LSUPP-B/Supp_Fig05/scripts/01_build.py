"""Supp Fig 5, lane V14_L3b_DURATION_SUPP, round 37 (2026-09-15): re-plot the four parts (a to d) at the V13 sheet's geometry from
the v8.1 numbers. Repointed copy of the round-30 builder (01_build_PRE_V8_1.py) for parts a and b, and of the round-32 LSUPP
rebuild (lsupp_01_build_cd_PRE_V8_1.py, the V13 design) for parts c and d: no in-panel summaries, the x axes of c and d end
where the bars need them at the SAME points-per-unit scale as the V13 sheet (the round-30 axis ends 21.0 and 48.5 units over
the full forest width, read back from the V13 bars by the verifier).

Data (v8.1): numbers/one_exposure_v1.json (step 124: the T90 and apnea exposures on all 19,173 nights, the sleep exposures on the
15,551 full nights with a sleep value, restricted rows within the 10,122 preserved-oxygen nights) and
Sleep_Variability_2026-08/reorg_fits_efig4_lag2/joint_t90_tst_landmark.json + validation.json (step 252, 52 outcomes, the total
sleep time term complete-case on the full nights). Nothing is pinned to a literal result: the row set (top 14 by the T90 > 10 percent
hazard ratio among the non-circular, non-control conditions common to the four whole-cohort exposures), every count, median and
n is derived here from the files and written to work/expected_values.json for the verifier. Printed HR (95% CI) strings are half
up (Decimal) on the file's values; one_exposure_v1.json stores 3 decimals, so strings whose stored value sits on an exact half at
two decimals are listed in expected_values.json (exact_half_strings) for the report. The v8.1 controls are alopecia, inguinal
hernia, glaucoma, contact dermatitis, hemorrhoids (never drawn here: the panels show the 14 clinical rows).

Writes work/part_a.pdf .. part_d.pdf, work/expected_values.json, work/build_log.txt."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import json
import os
import sys
from decimal import Decimal

import numpy as np

LANE = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LSUPP-B/Supp_Fig05"   # R49: this lane's sheet folder (scripts/, work/, verify/ beside it)
sys.path.insert(0, f"{LANE}/scripts")
from l3b_common import (NUM, SV, INK, BLUE, ORANGE, ORANGE_PALE, GREEN, GREEN_MID, GREEN_PALE,  # noqa: E402
                        GREEN_BAND, LABF, TCKF, ANNF, PT_PLUS, CI_LW, MEDGE, S_CIRCLE, S_SQUARE, MS_CIRCLE,
                        MS_SQUARE, rc_polish, measure, wrap_measured, text_gate, sidecar, jload, spans, hydrated, BASE, plt)
import v14lib  # noqa: E402
from matplotlib import ticker as mticker  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.transforms import blended_transform_factory  # noqa: E402

SHEET = "Supp_Fig05"
SD = LANE   # round 38: outputs in the sheet folder
WORK = f"{SD}/work"
os.makedirs(WORK, exist_ok=True)
LOG = open(f"{WORK}/build_log.txt", "w")


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    LOG.write(s + "\n")


# ------------------------------------------------------------------ sources, gated (v8.1)
OE_PATH = f"{NUM}/one_exposure_v1.json"
LJ_PATH = f"{SV}/reorg_fits_efig4_lag2/joint_t90_tst_landmark.json"
VAL_PATH = f"{SV}/reorg_fits_efig4_lag2/validation.json"
SHA_OE, SC_OE = sidecar(OE_PATH)
SHA_LJ, SC_LJ = sidecar(LJ_PATH)
SHA_VAL, SC_VAL = sidecar(VAL_PATH)
assert SC_OE["step"].get("rc", 0) == 0 and SC_LJ["step"].get("rc", 0) == 0 and SC_VAL["step"].get("rc", 0) == 0
log("one_exposure_v1.json sha256", SHA_OE, "step", SC_OE["step"]["id"], "stamped", SC_OE["output"]["mtime_local"])
log("joint_t90_tst_landmark.json sha256", SHA_LJ, "step", SC_LJ["step"]["id"], "stamped", SC_LJ["output"]["mtime_local"])

D = jload(OE_PATH)
EX = D["exposures"]
assert D["normal_oxygen_definition"] == "spo2_pct_below_90 <= 1% of the recording"
COHORT_N = int(D["cohort_n"])
VAL = jload(VAL_PATH)
for gate in ("V1_lag0_vs_bdsp_diseases_v3", "V2_lag2_vs_lag_ladder"):
    assert len(VAL[gate]) == 3 and all(v["match"] for v in VAL[gate].values()), (
        f"landmark machinery validation failed at {gate}, panel d must not draw")
LJ = jload(LJ_PATH)
N_LJ = len(LJ["outcomes"])
assert LJ["provenance"]["n_cohort"] == COHORT_N, (LJ["provenance"]["n_cohort"], COHORT_N)
log(f"cohort {COHORT_N}, preserved-oxygen nights {D['normal_oxygen_n']}, joint landmark outcomes {N_LJ}")

# R49: the V13 text record from a Ghostscript census of the V13 sheet (work/base_text_v13.json), no PyMuPDF probe of a composed sheet
_BT = json.load(open(hydrated(f"{WORK}/base_text_v13.json"))); BASE_TEXT = _BT["spans"]; PAGE_W, PAGE_H = _BT["page"]
BASE_STRINGS = {s["text"].strip() for s in BASE_TEXT}


def hr_ci(hr, lo, hi):
    return v14lib.hr_ci(hr, lo, hi)


def on_half(x, nd=2):
    """True when the stored value sits exactly on a half at the printed precision (a full-precision twin could move the digit)."""
    d = Decimal(repr(float(x))).scaleb(nd)
    return d == d.to_integral_value() + Decimal("0.5") or (d - d.to_integral_value(rounding="ROUND_FLOOR")) == Decimal("0.5")


def bh_q(pmap):
    names = sorted(pmap, key=lambda k: pmap[k])
    m = len(names)
    q, prev = {}, 1.0
    for rank_from_top in range(m, 0, -1):
        nm = names[rank_from_top - 1]
        prev = min(prev, pmap[nm] * m / rank_from_top)
        q[nm] = prev
    return q


Q, POOL_N = {}, {}
for ename, e in EX.items():
    pool = {c: v["p"] for c, v in e["outcomes"].items() if not v["circular"] and not v["negative_control"]}
    POOL_N[ename] = len(pool)
    Q[ename] = bh_q(pool)
assert len(set(POOL_N.values())) == 1, POOL_N
N_POOL = list(POOL_N.values())[0]
log(f"BH pool per exposure: {N_POOL} non-circular, non-control conditions (v7 had 45)")

WHOLE_T90 = "Oxygen below 90% for more than 10% of the night"
WHOLE_SHORT = "Sleeping fewer than 5 hours"
AHI_ALL = "Severe sleep apnea, apnea-hypopnea index 30 or more"
AHI_NORM = "Severe sleep apnea, in normal oxygenation only"
SE_NORM = "Sleep efficiency below 85%, in normal oxygenation only"
N3_NORM = "N3 sleep below 5% of the night, in normal oxygenation only"
SHORT_NORM = "Sleeping fewer than 5 hours, in normal oxygenation only"

need = [WHOLE_T90, WHOLE_SHORT, AHI_ALL, AHI_NORM]
common = set(EX[WHOLE_T90]["outcomes"])
for e in need:
    common &= set(EX[e]["outcomes"])
common = {c for c in common if not EX[WHOLE_T90]["outcomes"][c]["circular"]
          and not EX[WHOLE_T90]["outcomes"][c]["negative_control"]}
CONDS = sorted(common, key=lambda c: -EX[WHOLE_T90]["outcomes"][c]["hr"])[:14]   # the builder's rule
assert len(CONDS) == 14, len(CONDS)
ROWS_IN = [c for c in CONDS if c not in BASE_STRINGS]
BASE_ROWS = [s["text"] for s in BASE_TEXT if abs(s["x"] - 26.96) < 0.5 and s["size"] == 10.0 and s["y"] < 380]
ROWS_OUT = [c for c in BASE_ROWS if c not in CONDS]
log("row set (top 14 by T90 hazard ratio):", CONDS)
log("rows in:", ROWS_IN, "rows out:", ROWS_OUT)

# panel d: recount from the per-outcome refit, never from its summary block
LARMS = [("t90", "0"), ("t90", "2"), ("tst", "0"), ("tst", "2")]
LCOUNT, LMED = {}, {}
for tag, lag in LARMS:
    hrs = [LJ["outcomes"][k][lag][f"{tag}_hr"] for k in LJ["outcomes"]]
    qs = [LJ["outcomes"][k][lag][f"{tag}_q"] for k in LJ["outcomes"]]
    LCOUNT[(tag, lag)] = sum(1 for h, qq in zip(hrs, qs) if qq < 0.05 and h > 1)
    LMED[(tag, lag)] = float(np.median(hrs))
log("panel d counts", LCOUNT, "medians", {k: round(v, 3) for k, v in LMED.items()})

# ------------------------------------------------------------------ geometry, inches (the builder's)
W = 171.0 / 25.4
EDGE = 0.18
GUT = EDGE
AX0 = 1.98
FW = W - AX0 - 0.20
DW = 2.76
COL_HR = AX0 + DW + 0.12
COL_Q = COL_HR + 1.16 + 0.02   # R49: the q column nudged 0.02 in (1.44 pt) right: at 10.5 pt the HR rule (worst case '8.88 (8.88-88.88)' ends 0.06 in before the q column) fell short by 1.11 pt; the q rule keeps 4.8 pt of room; text only, no data mark moves
XTICKS = [0.7, 1, 2, 4, 8]
V13_XLIM = (0.555, 13.0)          # the V13 (round-30) axis range of parts a and b
# v8.1: the axis range follows the data (the same rule ED_Fig02 and Supp_Fig04 applied this round): every drawn interval must sit
# inside it, the ticks stay the V13 ticks, and a limit is moved only when a value would fall off the V13 axis.
def _drawn_lo_hi():
    los, his = [], []
    for ename in (WHOLE_T90, WHOLE_SHORT, AHI_ALL, AHI_NORM):
        for c in CONDS:
            o = EX[ename]["outcomes"][c]; los.append(o["lo"]); his.append(o["hi"])
    return min(los), max(his)
_lo, _hi = _drawn_lo_hi()
XLIM = (min(V13_XLIM[0], round(_lo / 1.03, 4)), max(V13_XLIM[1], round(_hi * 1.03, 4)))
PITCH = 0.30
N_UNITS = 0.9 + 14 + 0.65
AH = N_UNITS * PITCH
C_PITCH = 0.40
C_GAP = 0.5
C_UNITS = 8.2
CH = C_UNITS * C_PITCH
DH = 4.30 * 0.52
TOP_M, XLBL, LEG, CXL, BOT_M = 0.34, 0.46, 0.48, 0.50, 0.24
H_AB = TOP_M + AH + XLBL + LEG + BOT_M
H_C = TOP_M + CH + CXL + BOT_M
H_D = TOP_M + DH + CXL + BOT_M
yAX = BOT_M + LEG + XLBL
yCX = BOT_M + CXL
PLACE = {"a": (14.0, 14.0), "b": (522.72, 14.0), "c": (14.0, 483.32), "d": (522.72, 483.32)}
PAGE = (round(PAGE_W, 2), round(PAGE_H, 2))
assert PAGE == (1021.5, 811.2), PAGE
assert abs(PLACE["a"][1] + H_AB * 72 - 459.32) < 0.02 and abs(PLACE["c"][1] + H_C * 72 - 797.24) < 0.02
assert abs(PLACE["d"][1] + H_D * 72 - 722.07) < 0.02
# the V13 (round-32) scale of parts c and d: the round-30 axes spanned the full forest width for 21.0 and 48.5 units
V13_XHI_C, V13_XHI_D = 21.0, 48.5
PER_UX_C, PER_UX_D = FW / V13_XHI_C, FW / V13_XHI_D

# ------------------------------------------------------------------ the sheet's strings
LAB_A_BLUE, LAB_A_GREEN = "T90 >10%", "Sleep <5 h"
LAB_B_DARK = f"Severe sleep apnea, preserved oxygenation (n exposed = {EX[AHI_NORM]['n_exposed']:,})"
LAB_B_PALE = f"All severe sleep apnea (n exposed = {EX[AHI_ALL]['n_exposed']:,})"
XLAB_AB = "Rate of new diagnosis (95% CI)"
XLAB_C = "Conditions with a raised rate of new diagnosis"
XLAB_D = f"Conditions with a raised rate (FDR q < 0.05), of {N_LJ}"
NOTE_C = f"preserved nocturnal oxygenation only, n = {D['normal_oxygen_n']:,}"
ORDER = [(WHOLE_T90, "T90 >10%", False, BLUE),
         (AHI_ALL, "All severe sleep apnea", False, ORANGE_PALE),
         (AHI_NORM, "Severe sleep apnea", True, ORANGE),
         (N3_NORM, "N3 sleep below 5%", True, GREEN_MID),
         (SHORT_NORM, "Sleep <5 h", True, GREEN_MID),
         (SE_NORM, "Sleep efficiency below 85%", True, GREEN_MID),
         (WHOLE_SHORT, "Sleeping fewer than 5 hours, whole cohort", False, GREEN_PALE)]
RESTRICTED = [r for _e, _l, r, _c in ORDER]
SUB = [i for i, r in enumerate(RESTRICTED) if r]
assert SUB == [2, 3, 4, 5], SUB
assert all(ORDER[i][0].endswith("in normal oxygenation only") for i in SUB)
DROWS = [("t90", "0", "Time below 90% saturation, per SD, all follow-up", BLUE),
         ("t90", "2", "Time below 90% saturation, per SD, 2-year landmark", BLUE),
         ("tst", "0", "Total sleep time, per SD, all follow-up", GREEN_PALE),
         ("tst", "2", "Total sleep time, per SD, 2-year landmark", GREEN_PALE)]
for s in (LAB_A_BLUE, LAB_A_GREEN, XLAB_AB, XLAB_C, "HR (95% CI)", "q") + tuple(l for _e, l, _r, _c in ORDER[:-1]):
    assert s in BASE_STRINGS, f"string not on the V13 sheet: {s!r}"


def labels_at_plus(fig, labels, max_in):
    """R49: each gutter label wrapped at the larger size; a label whose wrap would change against the V26 wrap (at TCKF - PT_PLUS) keeps
    its V26 wrap and its current size (the brief: never re-wrap, never move a data mark, keep that one element at its current size).
    Returns (wrapped strings, per-label sizes, kept list)."""
    out, sizes, kept = [], [], []
    for lab in labels:
        w_new = wrap_measured(fig, lab, max_in, TCKF)
        w_old = wrap_measured(fig, lab, max_in, TCKF - PT_PLUS)
        if w_new == w_old:
            out.append(w_new); sizes.append(TCKF)
        else:
            out.append(w_old); sizes.append(TCKF - PT_PLUS); kept.append(lab)
    return out, sizes, kept

EXPECTED = {"sheet": SHEET, "sources": {
    "one_exposure_v1.json": {"path": OE_PATH, "sha256": SHA_OE, "step": SC_OE["step"], "mtime": SC_OE["output"]["mtime_local"]},
    "joint_t90_tst_landmark.json": {"path": LJ_PATH, "sha256": SHA_LJ, "step": SC_LJ["step"], "mtime": SC_LJ["output"]["mtime_local"]},
    "validation.json": {"path": VAL_PATH, "sha256": SHA_VAL}},
    "row_set": CONDS, "rows_in": ROWS_IN, "rows_out": ROWS_OUT, "bh_pool_n": N_POOL, "cohort_n": COHORT_N, "n_lj_outcomes": N_LJ,
    "values": [], "exact_half_strings": []}


def expect(panel, text, source, key, value, rule="text"):
    EXPECTED["values"].append(dict(panel=panel, text=text, source=source, key=key, value=value, rule=rule))


def fmt_q(q):
    return "<0.001" if q < 0.001 else v14lib.halfup(q, 3)


def left_rows(ax, ys, labels, gutter_x, ax_x0, sizes=None):
    ax.set_yticks(list(ys))
    ax.set_yticklabels(labels, fontsize=TCKF, color=INK, ha="left")
    if sizes is not None:   # R49: per-label size (a label kept at its current size when its wrap would change)
        for t, sz in zip(ax.get_yticklabels(), sizes):
            t.set_fontsize(sz)
    ax.tick_params(axis="y", pad=(ax_x0 - gutter_x) * 72.0, length=0)


def bare_y(ax):
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)


def logx(ax, ticks):
    ax.set_xscale("log")
    ax.set_xticks(ticks)
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.xaxis.set_minor_locator(mticker.NullLocator())


def finish(fig, name):
    text_gate(fig)
    import matplotlib.text as mtext
    for t in fig.findobj(mtext.Text):
        assert "uncontrolled" not in t.get_text().lower(), t.get_text()
    out = f"{WORK}/part_{name}.pdf"
    fig.savefig(out)
    plt.close(fig)
    log("wrote", out)
    return out


# ------------------------------------------------------------------ parts a and b
def pair_panel(fig, ax, h_in, panel, e_prim, e_sec, col_prim, col_prim_fill, col_sec, col_sec_fill,
               lab_prim, lab_sec, leg_y):
    """Paired-offset forest. Filled marker = q < 0.05 within exposure (BH over the pool), open = not.
    The primary series (upper offset, circles) prints HR (95% CI) and q."""
    tr = blended_transform_factory(ax.transAxes, ax.transData)
    x_hr = (COL_HR - AX0) / DW
    x_q = (COL_Q - AX0) / DW
    y = np.arange(len(CONDS))[::-1].astype(float)
    for offs, ename, colour, fillc, mk, s_mk in ((0.20, e_prim, col_prim, col_prim_fill, "o", S_CIRCLE),
                                                 (-0.20, e_sec, col_sec, col_sec_fill, "s", S_SQUARE)):
        o = EX[ename]["outcomes"]
        for i, c in enumerate(CONDS):
            hr, lo, hi = o[c]["hr"], o[c]["lo"], o[c]["hi"]
            assert XLIM[0] < lo and hi < XLIM[1], (ename, c, lo, hi, "off the axis")
            qv = Q[ename][c]
            sig = qv < 0.05
            ax.plot([lo, hi], [y[i] + offs] * 2, color=colour, lw=CI_LW, solid_capstyle="round", zorder=2)
            ax.scatter([hr], [y[i] + offs], s=s_mk, marker=mk, zorder=3,
                       facecolors=fillc if sig else "white", edgecolors="white" if sig else colour,
                       linewidths=MEDGE if sig else 1.0)
            if offs > 0:
                hr_txt = hr_ci(hr, lo, hi)
                ax.text(x_hr, y[i], hr_txt, transform=tr, fontsize=ANNF, va="center", ha="left", color=INK)
                ax.text(x_q, y[i], fmt_q(qv), transform=tr, fontsize=ANNF, va="center", ha="left",
                        color=INK, fontweight="bold" if sig else "normal")
                expect(panel, hr_txt, OE_PATH, f"exposures/{ename}/outcomes/{c}", [hr, lo, hi], rule="hr_ci_2dp_dict")
                expect(panel, fmt_q(qv), OE_PATH, f"BH q over the {N_POOL} non-circular non-control conditions of exposures/{ename}, {c}", qv)
                for v in (hr, lo, hi):
                    if on_half(v): EXPECTED["exact_half_strings"].append(dict(panel=panel, exposure=ename, condition=c, stored=v, printed=hr_txt))
            EXPECTED.setdefault("markers", []).append(dict(panel=panel, exposure=ename, condition=c, hr=hr, lo=lo,
                                                           hi=hi, q=qv, filled=bool(sig)))
    ax.text(x_hr, len(CONDS) + 0.12, "HR (95% CI)", transform=tr, fontsize=ANNF, va="center", ha="left", color=INK)
    ax.text(x_q, len(CONDS) + 0.12, "q", transform=tr, fontsize=ANNF, va="center", ha="left", color=INK)
    ax.axvline(1.0, color=INK, lw=0.9, ls=(0, (4, 3)), zorder=1)
    left_rows(ax, y, CONDS, GUT, AX0)
    bare_y(ax)
    logx(ax, XTICKS)
    ax.set_xlim(*XLIM)
    ax.set_ylim(-0.9, len(CONDS) + 0.65)
    ax.set_xlabel(XLAB_AB, fontsize=LABF, labelpad=4)
    fig.legend(handles=[
        Line2D([], [], color=col_prim, marker="o", ls="-", lw=CI_LW, ms=MS_CIRCLE, markerfacecolor=col_prim_fill,
               markeredgecolor="white", markeredgewidth=MEDGE, label=lab_prim),
        Line2D([], [], color=col_sec, marker="s", ls="-", lw=CI_LW, ms=MS_SQUARE, markerfacecolor=col_sec_fill,
               markeredgecolor="white", markeredgewidth=MEDGE, label=lab_sec)],
        loc="upper left", bbox_to_anchor=(GUT / W, leg_y / h_in), ncol=1, fontsize=TCKF, frameon=False,
        handletextpad=0.6, labelspacing=0.45, handlelength=2.2)


def forest_widths(fig):
    for c in CONDS:
        assert GUT + measure(fig, c, TCKF) <= AX0 - 0.06, f"row label too wide: {c}"
    assert COL_HR + measure(fig, "8.88 (8.88-88.88)", ANNF) <= COL_Q - 0.06
    assert COL_Q + measure(fig, "<0.001", ANNF) <= W - EDGE + 0.005


def build_forest(name, e_prim, e_sec, col_prim, col_prim_fill, col_sec, col_sec_fill, lab_prim, lab_sec):
    fig = plt.figure(figsize=(W, H_AB))
    fig.canvas.draw()
    forest_widths(fig)
    ax = fig.add_axes([AX0 / W, yAX / H_AB, DW / W, AH / H_AB])
    pair_panel(fig, ax, H_AB, name, e_prim, e_sec, col_prim, col_prim_fill, col_sec, col_sec_fill,
               lab_prim, lab_sec, yAX - XLBL - 0.02)
    return finish(fig, name)


# ------------------------------------------------------------------ parts c and d (the V13 round-32 design: no in-panel summaries, shortened axes at the V13 scale)
counts, meds = [], []
for ename, _lab, _r, _c in ORDER:
    o = EX[ename]["outcomes"]
    counts.append(sum(1 for c in CONDS if o[c]["lo"] > 1))
    meds.append(float(np.median([o[c]["hr"] for c in CONDS])))
log("panel c counts", counts, "medians", [round(m, 3) for m in meds])
# R49: the note (and its bracket) left the sheet in round 40 (LNOTES, Alen 2026-09-18): not drawn, not expected
EXPECTED["r40_lnotes_removed"] = {"note": NOTE_C, "bracket": "three 0.8 pt strokes in the note colour beside the band", "band": "kept"}
expect("b", LAB_B_DARK, OE_PATH, f"exposures/{AHI_NORM}/n_exposed", EX[AHI_NORM]["n_exposed"])
expect("b", LAB_B_PALE, OE_PATH, f"exposures/{AHI_ALL}/n_exposed", EX[AHI_ALL]["n_exposed"])
expect("d", XLAB_D, LJ_PATH, "number of outcomes in the joint landmark file", N_LJ)
d_counts = [LCOUNT[(t, g)] for t, g, _l, _c in DROWS]
d_meds = [LMED[(t, g)] for t, g, _l, _c in DROWS]
EXPECTED["panelC"] = {e: dict(count=c, median=m, n_exposed=EX[e]["n_exposed"]) for (e, _l, _r, _c), c, m in zip(ORDER, counts, meds)}
EXPECTED["panelD"] = {f"{t}_lag{g}": dict(count=c, median=m) for (t, g, _l, _c), c, m in zip(DROWS, d_counts, d_meds)}
EXPECTED["n_exposed"] = {e: EX[e]["n_exposed"] for e in EX}
EXPECTED["normal_oxygen_n"] = D["normal_oxygen_n"]
XHI_C = float(np.ceil(((max(counts) + 0.32) * FW / (FW - 0.10)) * 2) / 2)    # the LSUPP rule (no annotation to make room for)
XHI_D = float(np.ceil(((max(d_counts) + 0.4) * FW / (FW - 0.10)) * 2) / 2)
AXW_C, AXW_D = XHI_C * PER_UX_C, XHI_D * PER_UX_D
D_TICKS = list(range(0, int(np.floor(XHI_D)) + 1, 10))
log(f"part c: x axis 0 to {XHI_C} at {PER_UX_C * 72:.4f} pt per unit (V13 scale), axes width {AXW_C:.4f} in")
log(f"part d: x axis 0 to {XHI_D} at {PER_UX_D * 72:.4f} pt per unit (V13 scale), axes width {AXW_D:.4f} in, ticks {D_TICKS}")


def c_row_positions():
    ys, cur = [], 0.0
    for i in range(len(ORDER) - 1, -1, -1):
        if i < len(ORDER) - 1:
            cur += 1.0 + (C_GAP if RESTRICTED[i] != RESTRICTED[i + 1] else 0.0)
        ys.append(cur)
    ys = ys[::-1]
    assert ys == [7.0, 6.0, 4.5, 3.5, 2.5, 1.5, 0.0], ys
    return np.array(ys)


GEOM = {}


def build_summary(name):
    fig = plt.figure(figsize=(W, H_C))
    fig.canvas.draw()
    ax = fig.add_axes([AX0 / W, yCX / H_C, AXW_C / W, CH / H_C])
    labs, sizes_c, kept_c = labels_at_plus(fig, [lab for _e, lab, _r, _c in ORDER], AX0 - GUT - 0.06)
    log("part c: labels kept at the current size (wrap would change at +1 pt):", kept_c)
    y = c_row_positions()
    ax.barh(y, counts, height=0.66, color=[c for _e, _l, _r, c in ORDER], linewidth=0, zorder=2)
    left_rows(ax, y, labs, GUT, AX0, sizes=sizes_c)
    bare_y(ax)
    ax.set_xlim(0, XHI_C)
    ax.set_xticks(range(0, len(CONDS) + 1, 2))
    ax.set_ylim(-0.6, C_UNITS - 0.6)
    ax.set_xlabel(XLAB_C, fontsize=LABF, labelpad=4)
    band_lo, band_hi = y[SUB[-1]] - 0.45, y[SUB[0]] + 0.45
    assert band_hi < y[SUB[0] - 1] - 0.33 and band_lo > y[SUB[-1] + 1] + 0.33
    ax.axhspan(band_lo, band_hi, color=GREEN_BAND, lw=0, zorder=0)
    per_ux = AXW_C / XHI_C
    assert abs(per_ux - PER_UX_C) < 1e-12
    # R49: the round-40 LNOTES step folded in (Alen, 2026-09-18): no bracket and no italic note beside the band; the band stays
    GEOM["c"] = dict(xhi=XHI_C, per_ux_in=per_ux, axes_in=[AX0, yCX, AXW_C, CH], band=[band_lo, band_hi], rows_y=list(map(float, y)),
                     kept_labels=kept_c, label_sizes=sizes_c)
    log("part c: band", [round(band_lo, 3), round(band_hi, 3)], "no bracket, no note (round 40)")
    return finish(fig, name)


def build_landmark(name):
    fig = plt.figure(figsize=(W, H_D))
    fig.canvas.draw()
    ax = fig.add_axes([AX0 / W, yCX / H_D, AXW_D / W, DH / H_D])
    d_labs, sizes_d, kept_d = labels_at_plus(fig, [lab for _t, _g, lab, _c in DROWS], AX0 - GUT - 0.06)
    if kept_d:   # R49: the four two-line labels of d form one column: when any would re-wrap, the column keeps its current size
        d_labs = [wrap_measured(fig, lab, AX0 - GUT - 0.06, TCKF - PT_PLUS) for _t, _g, lab, _c in DROWS]; sizes_d = [TCKF - PT_PLUS] * len(DROWS)
    log("part d: labels kept at the current size:", [lab for _t, _g, lab, _c in DROWS] if kept_d else [], "(re-wrap at +1 pt:", kept_d, ")")
    yd = np.arange(len(DROWS))[::-1].astype(float)
    ax.barh(yd, d_counts, height=0.62, color=[c for _t, _g, _l, c in DROWS], linewidth=0, zorder=2)
    left_rows(ax, yd, d_labs, GUT, AX0, sizes=sizes_d)
    bare_y(ax)
    ax.set_xlim(0, XHI_D)
    ax.set_xticks(D_TICKS)
    ax.set_ylim(-0.65, 3.65)
    ax.set_xlabel(XLAB_D, fontsize=LABF, labelpad=4)
    GEOM["d"] = dict(xhi=XHI_D, per_ux_in=AXW_D / XHI_D, axes_in=[AX0, yCX, AXW_D, DH], rows_y=list(map(float, yd)), ticks=D_TICKS,
                     kept_labels=([lab for _t, _g, lab, _c in DROWS] if kept_d else []), label_sizes=sizes_d)
    return finish(fig, name)


with plt.rc_context(rc_polish()):
    build_forest("a", WHOLE_T90, WHOLE_SHORT, BLUE, BLUE, GREEN, GREEN_MID, LAB_A_BLUE, LAB_A_GREEN)
    build_forest("b", AHI_NORM, AHI_ALL, ORANGE, ORANGE, ORANGE_PALE, ORANGE_PALE, LAB_B_DARK, LAB_B_PALE)
    build_summary("c")
    build_landmark("d")

EXPECTED["place"] = PLACE
EXPECTED["page"] = PAGE
EXPECTED["part_heights_pt"] = {"a": H_AB * 72, "b": H_AB * 72, "c": H_C * 72, "d": H_D * 72}
EXPECTED["labels"] = dict(a=[LAB_A_BLUE, LAB_A_GREEN], b=[LAB_B_DARK, LAB_B_PALE],
                          c=[l for _e, l, _r, _c in ORDER], d=[l for _t, _g, l, _c in DROWS])   # R49: no note
EXPECTED["xlabels"] = dict(ab=XLAB_AB, c=XLAB_C, d=XLAB_D)
EXPECTED["counts_c"] = counts
EXPECTED["counts_d"] = d_counts
EXPECTED["geometry"] = GEOM
EXPECTED["xlim_ab"] = dict(v13=list(V13_XLIM), new=list(XLIM), moved=[XLIM[0] != V13_XLIM[0], XLIM[1] != V13_XLIM[1]], drawn_lo=_lo, drawn_hi=_hi)
EXPECTED["v13_scale"] = dict(xhi_c=V13_XHI_C, xhi_d=V13_XHI_D, per_ux_c_in=PER_UX_C, per_ux_d_in=PER_UX_D)
EXPECTED["r49_text_sizes"] = dict(pt_plus=PT_PLUS, LABF=LABF, TCKF=TCKF, ANNF=ANNF, letters=13.0 + PT_PLUS)
EXPECTED["rounding_note"] = ("one_exposure_v1.json stores three decimals; printed HR (95% CI) strings are half up on those stored values "
                             "(no full-precision twin exists for v8.1); exact_half_strings lists the stored values that sit on a half at two decimals")
json.dump(EXPECTED, open(f"{WORK}/expected_values.json", "w"), indent=1, default=float)
log("wrote", f"{WORK}/expected_values.json", "with", len(EXPECTED["values"]), "expected printed values;",
    len(EXPECTED["exact_half_strings"]), "stored values on an exact half at 2 dp")
LOG.close()
