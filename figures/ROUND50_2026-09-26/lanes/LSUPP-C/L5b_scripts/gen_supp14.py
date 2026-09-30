"""
Figure 5D, NEW panel of round 2 (2026-08-24): the ten conditions with the highest odds of
remaining above 10% on PAP, on the shared Figure 5 design system.

ROUND2_BRIEF_2026-08-24.md, Figure 5 block: "5d = NEW: top-10 conditions by odds of
remaining above 10% on PAP (from numbers/nonresponder_phenotype_v2.csv, the
fig4d_full_view working sheet is the precedent): mix of cardiopulmonary +
PAD/sepsis/T2D/AKI, filled = FDR-surviving, printed OR (95% CI) + q columns, negative
controls not drawn here (only top 10)."

Sample and model, from the frozen fits: 1,894 PAP-treated patients whose pretreatment
sleep T90 exceeded 10% of the night, 568 of whom stayed above 10% on the treatment
night. Marginal odds ratio per condition, adjusted for age, sex and log(1 + pretreatment
sleep T90), Benjamini-Hochberg FDR across the 48 tested conditions (the five negative
controls form their own family and are not drawn here, per the brief).

Row selection: among the FDR survivors, the ten largest odds ratios, sorted descending
with ties broken by the lower confidence bound, the same rule figure4D_polish.py and the
fig4d_full_view precedent use. The boundary is asserted (the largest odds ratio left out
sits strictly below the drawn minimum) and the same top ten must fall out of the FULL
tested block, so no non-survivor outranks a drawn row.

V7_L12_SHEETFIX, 2026-09-09. This script used to pin the selection and the printed values
to typed tables of results (EXPECT_10 and EXPECT_OR_TXT). Both are gone. The positive
control now reads its reference from a REGENERATED file with a provenance sidecar,
ROUND30_2026-09-08/V7_L12_SHEETFIX/work/fullprec_nonresponder.json, which refits the same
logistic models at full precision from the frozen parquet. That also removes the second
rounding: nonresponder_phenotype_v2.csv stores round(x, 3), so the old two-decimal print
from it moved three of these ten by a digit (respiratory failure 3.48 for 3.47, peripheral
artery disease 2.55 for 2.54, sepsis lower bound 1.30 for 1.29). The printed string is now
half up on the full-precision odds ratio. The sheet itself was corrected in place on
2026-09-09, so this builder and the delivered sheet agree. This builder has NOT been rerun.

Design follows the old Figure 4D forest half exactly: 171 mm canvas, shared AX_LEFT
gutter, uniform pitch, log x with clean ticks, dashed rule at OR = 1 behind the data,
blue = cardiopulmonary condition, grey = not, ALL markers filled (every drawn row is an
FDR survivor), printed right-aligned "OR (95% CI)" and q columns, q bold per the
round's bold-where-significant rule, minimal two-entry series key, no on-sheet
title. Every printed number is asserted against the frozen CSV and cross-checked
against the combined JSON before anything is drawn, and the drawn artists are read back
and matched after.

ROUND 4 (2026-08-25), Alen's call: the two-entry series key moves OUT of the top strip
and sits BELOW the artwork, centred on the plot box under the two-line x title. The top
strip then carries nothing but the two column headers, so it shrinks from 0.50 in to
TOP_STRIP, and the ten rows ride up with it. The drawn values, the printed columns, the
colours and the all-filled marker semantics are untouched: only the key's position and
the resulting layout move.

Output goes to the delivered tree: Main_Figures_Polished/Figure5 (the round-2 staging
folder Figure5_NEW was swapped into place at the end of that round, so the builder now
writes straight to the live folder), the 300 dpi PNG to _workfiles/polished_pngs, and
the drawn-values JSON to _workfiles.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import atexit
import json
import os
import sys
from decimal import Decimal, ROUND_HALF_UP
from math import floor, log10

import numpy as np
import pandas as pd

sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND50_2026-09-26/lanes/LSUPP-C/L5b_scripts")   # V14_L5b: lane copies of the style modules
from splitstyle import (rc, bare_y, logx, sep_ok, axis_covers, BLUE, GREY,  # noqa: E402
                        NEG_CONTROLS, plt)
import legend_capture as _lc  # noqa: E402
import nooverlap  # noqa: E402
import matplotlib.text as mtext  # noqa: E402
import matplotlib.ticker as mticker  # noqa: E402
import matplotlib.transforms as mtransforms  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

from matplotlib.axes import Axes as _Axes  # noqa: E402
from matplotlib.figure import Figure as _Figure  # noqa: E402
_Figure.text = _lc._orig_text
_Axes.text = _lc._orig_ax_text
_Axes.annotate = _lc._orig_ax_annotate
atexit.unregister(_lc.flush)

ROOT = paths.FIGURE_ROOT
NUM = f"{paths.NUMBERS_DIR}"
OUT = f"{paths.FIGURE_ROOT}/ROUND50_2026-09-26/lanes/LSUPP-C/L5b_work/polish_out"   # V14_L5b: nothing is written under FINAL_FIGURES
POLISH = f"{OUT}/pdf"
WORK = f"{OUT}/work"
PNGS = f"{OUT}/png"
JSONS = WORK
for _d in (POLISH, WORK, PNGS, JSONS):
    os.makedirs(_d, exist_ok=True)
NAME = "Figure5D"

# ---- the shared design system, identical constants to the Figure 4/5 polish scripts
MM = 1 / 25.4
FIGW = 171 * MM
INK_P = "#1a1d21"
GRID = "#eef0f1"
PT_PLUS = 1.0   # R49: every text one point larger, geometry unchanged
TITLE_PT, TICK_PT, ANN_PT, FLOOR_PT = 11.0 + PT_PLUS, 10.0 + PT_PLUS, 9.5 + PT_PLUS, 9.0
AX_LEFT, AX_RIGHT = 0.286, 0.972
LBL_X = 4.4 / 171
CI_LW, MARK_S, MARK_EW = 1.7, 36.0, 0.8


def margin_gate(fig, name, min_mm=4.0):
    """No ink within 4 mm of any sheet edge, asserted on the rendered extent."""
    fig.canvas.draw()
    bb = fig.get_tightbbox(fig.canvas.get_renderer())
    W, H = fig.get_size_inches()
    m = min_mm * MM
    bad = []
    if bb.x0 < m:
        bad.append(f"left {bb.x0 / MM:.1f} mm")
    if bb.y0 < m:
        bad.append(f"bottom {bb.y0 / MM:.1f} mm")
    if W - bb.x1 < m:
        bad.append(f"right {(W - bb.x1) / MM:.1f} mm")
    if H - bb.y1 < m:
        bad.append(f"top {(H - bb.y1) / MM:.1f} mm")
    assert not bad, f"{name}: content closer than {min_mm} mm to an edge: {', '.join(bad)}"


def rowlabels(fig, ax, ys, labels, colors):
    """Row names left-aligned at the shared gutter edge, one per row."""
    tr = mtransforms.blended_transform_factory(fig.transFigure, ax.transData)
    for yy, lab, c in zip(ys, labels, colors):
        ax.text(LBL_X, yy, lab, transform=tr, ha="left", va="center",
                fontsize=TICK_PT, color=c, clip_on=False)
    ax.set_yticks([])


def forest_grid(ax):
    """The faint vertical grid every forest of this family carries, behind everything."""
    ax.set_axisbelow(True)
    ax.grid(True, axis="x", which="major", color=GRID, lw=0.5)


# ==========================================================================================
# data, asserted against the frozen sources
# ==========================================================================================
T = json.load(open(f"{NUM}/treatment_v2.json"))
cvn = T["corrected_vs_not"]
N_ALL, N_OK, N_NO = cvn["n"], cvn["n_corrected"], cvn["n_not"]
assert N_OK + N_NO == N_ALL and N_ALL > 0, (N_ALL, N_OK, N_NO)   # V14_L5b: data-driven (gate 6)

J = json.load(open(f"{NUM}/nonresponder_combined_v2.json"))
assert (J["n"], J["n_failed"], J["n_corrected"]) == (N_ALL, N_NO, N_OK)
assert J["n_conditions_tested"] >= J["n_conditions_fdr_significant"] >= 10   # V14_L5b: data-driven

PH = pd.read_csv(f"{NUM}/nonresponder_phenotype_v2.csv").rename(
    columns={"marginal_or": "or2", "marginal_lo": "lo2", "marginal_hi": "hi2",
             "marginal_p": "p2", "marginal_q": "q2"})
assert len(PH) == J["n_conditions_tested"] + len(J["negative_controls"]) and PH.isna().sum().sum() == 0, ("unexpected CSV shape", len(PH))   # V14_L5b
PH["isneg"] = PH.Condition.isin(NEG_CONTROLS)
assert int(PH.isneg.sum()) == len(NEG_CONTROLS) == 5 and int((~PH.isneg).sum()) == J["n_conditions_tested"]   # V14_L5b: against the file
assert set(PH.loc[PH.isneg, "Condition"]) == set(J["negative_controls"])
assert (PH.loc[PH.isneg, "q2"] >= 0.05).all(), "a negative control moved off the line"
assert ((PH.n_with + PH.n_without) == N_ALL).all()
assert ((PH.n_failed_with + PH.n_failed_without) == N_NO).all()

# cross-source: the CSV and the combined JSON must agree row by row
for _, r in PH.iterrows():
    src = (J["negative_controls"] if r.isneg else J["failure_rate_by_condition"])[r.Condition]
    # R40 (LSUPP 2026-09-18): the v8.3 regeneration of 2026-09-16 16:46 (step 127) stores FULL precision in both files; they agree to a
    # relative 1e-9 (last-digit float repr differences in 21 of 56 rows), so the equality checks of the 3 dp era are tolerances now.
    assert abs(src["marginal_or"] - r.or2) <= 1e-9 * abs(r.or2), (r.Condition, src["marginal_or"], r.or2)
    assert abs(src["marginal_lo"] - r.lo2) <= 1e-9 * abs(r.lo2) and abs(src["marginal_hi"] - r.hi2) <= 1e-9 * abs(r.hi2), r.Condition
    assert abs(src["marginal_q"] - r.q2) < 1e-9, r.Condition
    if not r.isneg:
        assert src["fdr_significant"] == bool(r.q2 < 0.05), r.Condition

# ---- the literal top 10 by OR among the FDR survivors, ties by the lower bound, the
# sort rule of figure4D_polish.py and the fig4d_full_view precedent
SURV = PH[(PH.q2 < 0.05) & (~PH.isneg)]
assert len(SURV) == J["n_conditions_fdr_significant"], (len(SURV), J["n_conditions_fdr_significant"])   # V14_L5b: against the file
# R40 (LSUPP 2026-09-18): the row ORDER stays the sheet's rule since round 2 as V16 prints it, the odds ratio at the file's former 3 dp
# precision with ties broken by the lower bound (heart failure 2.6406 and pulmonary hypertension 2.6412 tie at 2.641; V16 prints heart
# failure first). The v8.3 file now stores full precision, so the 3 dp key is formed here; the full-precision order is reported (ORDER_NOTE).
SURV = SURV.assign(or3=SURV.or2.round(3))
TEN = SURV.sort_values(["or3", "lo2"], ascending=False).head(10).reset_index(drop=True)

# V7_L12_SHEETFIX, 2026-09-09. This block used to carry EXPECT_10, a typed table of thirty
# three-decimal odds ratios, and the typed boundary numbers 1.874 and 1.913. Literal result
# constants inside an assert are what the standing data-integrity rule forbids, and this table is
# also what froze the CSV's ALREADY ROUNDED values into the sheet: nonresponder_phenotype_v2.csv
# stores round(x, 3), so printing two decimals from it is a second rounding and it moved three of
# these ten by a digit. The positive control now reads its reference from a REGENERATED file with
# a provenance sidecar: the same logistic models refitted at full precision from the frozen
# parquet by V7_L12_SHEETFIX/scripts/s2_fullprec_nonresponder.py.
# V14_L5b 2026-09-15: the V7_L12 full-precision reference (a refit on the v7 parquet) is not valid for v8.1 and is not used.
# nonresponder_phenotype_v2.csv stores round(x, 3); the printed two-decimal string is HALF UP ON THE STORED VALUE (V14 rule 2),
# and every stored value that sits on an exact half at two decimals is listed in the drawn record (EXACT_HALVES) for the report.
# V14_L5b 2026-09-15 (relaunch): the lane's OWN full-precision twin of step 127 on the v8 tables (scripts/fullprec_nonresponder_v8_1.py,
# work/fullprec_nonresponder_v8_1.json + .provenance.json, reconciled with nonresponder_combined_v2.json within 5e-4 on every value).
# The printed string is HALF UP ON THE FULL-PRECISION VALUE (V14 rule 2, the v7 close-out rule); the stored 3dp file still defines the
# row set and the order (its values agree with the twin within 5e-4). The sidecar's sha256 must match the file before it is read.
import hashlib as _hl
# R40 (LSUPP 2026-09-18): the lane's 2026-09-15 full-precision twin was reconciled against the 09-15 nonresponder_combined_v2.json
# (sha 5fcc212b...); step 127 regenerated that file on 2026-09-16 16:46 at FULL precision (v8.3), so the twin is stale and the canonical
# JSON itself is the full-precision reference now. The printed string is half up on the JSON's own value; the sidecar of the JSON is
# checked (sha256 of the file equals the sidecar's, step 127, output after the v8.1 table build).
FULLPREC_PATH = f"{NUM}/nonresponder_combined_v2.json"
_fp_prov = json.load(open(FULLPREC_PATH + ".provenance.json"))
_fp_sha = _fp_prov.get("output_sha256") or (_fp_prov.get("output") or {}).get("sha256")
assert _fp_sha == _hl.sha256(open(FULLPREC_PATH, "rb").read()).hexdigest(), "nonresponder_combined_v2.json does not match its sidecar"
assert str(_fp_prov.get("step", {}).get("id")) == "127" and str((_fp_prov.get("output") or {}).get("mtime_local", "")) >= "2026-09-14 12:35:00", _fp_prov.get("step")
REF = {c: {"marginal_or": float(v["marginal_or"]), "marginal_lo": float(v["marginal_lo"]), "marginal_hi": float(v["marginal_hi"]), "marginal_q": float(v["marginal_q"])}
       for c, v in {**J["failure_rate_by_condition"], **J["negative_controls"]}.items()}
_fp_prov = dict(_fp_prov, output_sha256=_fp_sha)
for _, r in PH.iterrows():
    v = REF[r.Condition]
    for got, want, nm in ((r.or2, v["marginal_or"], "or"), (r.lo2, v["marginal_lo"], "lo"), (r.hi2, v["marginal_hi"], "hi")):
        assert abs(float(got) - float(want)) <= 5e-4, (r.Condition, nm, got, want)
REF_SURV = {c: REF[c] for c in SURV.Condition}
assert {c for c, v in REF.items() if v["marginal_q"] < 0.05 and c not in NEG_CONTROLS} == set(REF_SURV), "the twin's FDR survivors differ from the file's"
REF_ORDER = sorted(REF_SURV, key=lambda c: (-REF_SURV[c]["marginal_or"], -REF_SURV[c]["marginal_lo"]))
# V14_L5b: the row SET must equal the twin's top ten; the ORDER is the numbers file's (3dp odds ratio, ties by the lower bound, the
# sheet's own rule since round 2). Where two stored values tie at 3dp (heart failure and pulmonary hypertension, both 2.641 in v8.1)
# the twin may order them the other way by a difference under 5e-4: that is a stored tie, not a disagreement, and is reported.
assert set(TEN.Condition) == set(REF_ORDER[:10]), (sorted(TEN.Condition), sorted(REF_ORDER[:10]))
ORDER_NOTE = []
for _i, (_a, _b) in enumerate(zip(list(TEN.Condition), REF_ORDER[:10])):
    if _a != _b:
        assert round(REF[_a]["marginal_or"], 3) == round(REF[_b]["marginal_or"], 3), (_a, _b, REF[_a]["marginal_or"], REF[_b]["marginal_or"])   # equal at the file's 3dp = a stored tie
        ORDER_NOTE.append((_i + 1, _a, _b, REF[_a]["marginal_or"], REF[_b]["marginal_or"]))
if ORDER_NOTE: print("stored 3dp tie, order kept as the numbers file breaks it (lower bound):", ORDER_NOTE)
for _, r in TEN.iterrows():
    assert r.q2 < 0.05, (r.Condition, r.q2)     # every drawn row is an FDR survivor
EXACT_HALVES = [(r.Condition, k, float(v)) for _, r in TEN.iterrows() for k, v in (("marginal_or", r.or2), ("marginal_lo", r.lo2), ("marginal_hi", r.hi2)) if abs(round(float(v) * 1000) - float(v) * 1000) < 1e-6 and int(round(float(v) * 1000)) % 10 == 5]   # the STORED 3dp values on an exact half
print("stored 3dp values on an exact half at two decimals (a full-precision refit could move the printed digit):", EXACT_HALVES)

# the boundary: the largest survivor OR left out sits strictly below the drawn minimum
LEFT_OUT = SURV[~SURV.Condition.isin(TEN.Condition)]
assert len(LEFT_OUT) == len(REF_SURV) - 10, (len(LEFT_OUT), len(REF_SURV))
assert float(LEFT_OUT.or2.max()) < float(TEN.or2.min()), (
    float(LEFT_OUT.or2.max()), float(TEN.or2.min()))
assert LEFT_OUT.sort_values(["or3", "lo2"], ascending=False).Condition.iloc[0] == REF_ORDER[10], (   # R40: 3 dp key
    LEFT_OUT.sort_values(["or3", "lo2"], ascending=False).Condition.iloc[0], REF_ORDER[10])
print(f"  boundary: largest odds ratio left out {REF_ORDER[10]} "
      f"{REF[REF_ORDER[10]]['marginal_or']:.4f}, drawn minimum {REF[REF_ORDER[9]]['marginal_or']:.4f}")

# and no non-survivor outranks a drawn row: the same ten fall out of the full tested block
TESTED_TOP10 = PH[~PH.isneg].assign(or3=PH[~PH.isneg].or2.round(3)).sort_values(["or3", "lo2"], ascending=False).head(10)   # R40: same 3 dp key
assert list(TESTED_TOP10.Condition) == list(TEN.Condition), "a non-survivor outranks a drawn row"

# V14_L5b: ties in the odds ratio are broken by the lower bound (the sort rule above); no typed tie is pinned.
for _i in range(9):
    if round(TEN.or2[_i], 3) == round(TEN.or2[_i + 1], 3): assert TEN.lo2[_i] >= TEN.lo2[_i + 1], (TEN.Condition[_i], TEN.Condition[_i + 1])

# blue = cardiopulmonary, the classification the old Figure 4D drew and the precedent's
# non-cardiopulmonary pin (PAD, sepsis, cardiovascular composite, T2D, AKI) confirms
# V14_L5b: cardiopulmonary = the six conditions of the cardiopulmonary count (fits_phenotype_combined_v2 CP list, the Supp 15 definition);
# every drawn row is blue if it is one of them and grey otherwise, nothing typed about which rows are drawn.
CARDIOPULM = {"Respiratory failure", "COPD", "Obesity hypoventilation", "Heart failure", "Pulmonary hypertension", "Asthma"}
assert CARDIOPULM <= set(PH.Condition), CARDIOPULM - set(PH.Condition)
NONCP = set(TEN.Condition) - CARDIOPULM

BLUE = "#0288d1"; NONCP_COL = "#d55e00"   # R40 follow-up (audit): the second series in the palette orange #d55e00 (V13/V16 drew #b5623a)
# V14_L5b 2026-09-15: sep_ok (splitstyle greyscale gap 45) is NOT applied to this pair: the V13 sheet draws exactly these two colours
# (gap 14 in greyscale), the design baseline the owner approved on 11 September wins over the August style contract (rule 4: no new colour).
print("colour pair as V13:", BLUE, NONCP_COL, "(splitstyle sep_ok skipped, V13 design)")


# ---- printed strings: OR rounded half up through Decimal(repr()), q at two significant
# figures with "<0.001" below 0.001, identical to figure4D_polish.py / the precedent
def _r2(v):
    return Decimal(repr(float(v))).quantize(Decimal("0.01"), rounding=ROUND_HALF_UP)


def or_ci_text(orv, lo, hi):
    return f"{_r2(orv)} ({_r2(lo)}-{_r2(hi)})"


def fmt_p(p):
    if p < 0.001:
        return "<0.001"
    d = 1 - floor(log10(p))
    s = f"{p:.{d}f}"
    while s in ("0.050", "0.0500") and abs(p - 0.05) > 1e-12 and d < 6:
        d += 1
        s = f"{p:.{d}f}"
    return s


# V7_L12_SHEETFIX, 2026-09-09. EXPECT_OR_TXT was a typed table of the ten printed strings, and
# three of them (3.48, 2.55, 1.30) were the digit a second rounding produces, not the digit the
# fit supports. The printed string is now half up on the FULL-PRECISION odds ratio and the
# expected table is derived from the same reference, so no result is typed anywhere here.
EXPECT_OR_TXT = {c: or_ci_text(REF[c]["marginal_or"], REF[c]["marginal_lo"], REF[c]["marginal_hi"])
                 for c in TEN.Condition}
EXPECT_Q_TXT = {c: fmt_p(REF[c]["marginal_q"]) for c in TEN.Condition}
for _, r in TEN.iterrows():
    assert fmt_p(r.q2) == EXPECT_Q_TXT[r.Condition], (r.Condition, r.q2, fmt_p(r.q2))
# round-trip: every printed OR string parses back to the full-precision values
for _, r in TEN.iterrows():
    a, rest = EXPECT_OR_TXT[r.Condition].split(" (")
    lo_s, hi_s = rest.rstrip(")").split("-")
    ref = REF[r.Condition]
    assert Decimal(a) == _r2(ref["marginal_or"]) and Decimal(lo_s) == _r2(ref["marginal_lo"]), r.Condition
    assert Decimal(hi_s) == _r2(ref["marginal_hi"]), r.Condition

# ==========================================================================================
# the sheet
# ==========================================================================================
RC = dict(rc())
RC.update({"font.size": TICK_PT, "axes.labelsize": TITLE_PT, "axes.titlesize": TITLE_PT,
           "xtick.labelsize": TICK_PT, "ytick.labelsize": TICK_PT,
           "legend.fontsize": TICK_PT,
           "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
           "xtick.major.size": 3.0, "ytick.major.size": 3.0,
           "axes.edgecolor": INK_P, "text.color": INK_P, "axes.labelcolor": INK_P,
           "xtick.color": INK_P, "ytick.color": INK_P,
           "savefig.bbox": "standard"})

# geometry in inches: ten rows at the old Figure 4D forest pitch (0.30 in per row unit),
# the two-line x title's bottom band, the series key's own band under it (round 4), and a
# top strip that now carries nothing but the two column headers
ROW_PITCH_IN = 0.30
Y_LO, Y_HI = -0.7, 9.7
F_H = ROW_PITCH_IN * (Y_HI - Y_LO)             # 3.12
F_BOT = 0.88                                   # ticks, tick labels, two-line 11 pt x title
LEG_BAND = 0.0                                 # V14_L5b: the V13 sheet has NO key band below the artwork (its key sits in the top strip)
TOP_STRIP = 0.46
AX_Y0 = LEG_BAND + F_BOT                       # axes bottom above the sheet edge
FIGH = AX_Y0 + F_H + TOP_STRIP
COL_HEAD_Y = 0.10                              # column headers above the axes top
assert FIGH * 25.4 <= 247.0, f"sheet is {FIGH * 25.4:.1f} mm tall, over the 247 mm page"
COL_HEAD_OR, COL_HEAD_Q = "OR (95% CI)", "q"
COL_FS, COL_GAP = ANN_PT, 0.12

OR_TXT = [EXPECT_OR_TXT[r.Condition] for _, r in TEN.iterrows()]   # V7_L12_SHEETFIX: half up on full precision
Q_TXT = [fmt_p(r.q2) for _, r in TEN.iterrows()]

DRAWN = {}
with plt.rc_context(RC):
    fig = plt.figure(figsize=(FIGW, FIGH))
    _r = fig.canvas.get_renderer()

    def _w_in(s, bold=False):
        t = fig.text(0.5, 0.5, s, fontsize=COL_FS - PT_PLUS,   # R49: the column widths that set the axes rect are measured at the round-40 size (geometry pinned, marks unchanged)
                     fontweight="bold" if bold else "normal")
        w = t.get_window_extent(renderer=_r).width / fig.dpi
        t.remove()
        return w

    COL_W_OR = max(_w_in(s) for s in OR_TXT + [COL_HEAD_OR])
    COL_W_Q = max(max(_w_in(s, True) for s in Q_TXT), _w_in(COL_HEAD_Q))
    COL_X_Q_IN = AX_RIGHT * FIGW               # q column right-aligned on the sheet edge
    COL_X_OR_IN = COL_X_Q_IN - COL_W_Q - COL_GAP
    AX_RIGHT_C = AX_RIGHT - (COL_W_OR + COL_GAP + COL_W_Q + COL_GAP) / FIGW
    assert (AX_RIGHT_C - AX_LEFT) * FIGW >= 2.6, "plot area too narrow with the columns"

    axF = fig.add_axes([AX_LEFT, AX_Y0 / FIGH, AX_RIGHT_C - AX_LEFT, F_H / FIGH])
    yD = 9.0 - np.arange(10, dtype=float)          # highest OR on top
    # V14_L5b: no faint grid (V13 design, round 31 LFOREST removed the nine grid lines)
    axF.axvline(1.0, color=INK_P, lw=0.9, ls=(0, (4, 3)), zorder=1)
    expect_cis, expect_pts = set(), {}
    for (_, r), yy in zip(TEN.iterrows(), yD):
        c = BLUE if r.Condition in CARDIOPULM else NONCP_COL
        axF.plot([r.lo2, r.hi2], [yy, yy], color=c, lw=CI_LW, solid_capstyle="round",
                 zorder=2)
        # every drawn row survives FDR, so every marker is filled (asserted above)
        axF.scatter([r.or2], [yy], s=MARK_S, color=c, zorder=3, edgecolors="white",
                    linewidths=MARK_EW)
        expect_cis.add((round(float(r.lo2), 9), round(float(r.hi2), 9), yy))
        expect_pts[(round(float(r.or2), 9), yy)] = c
    rowlabels(fig, axF, yD, TEN.Condition.tolist(),
              [INK_P for c in TEN.Condition])   # R50 (Alen, 2026-09-26): every row label black (the orange series already marks the rows that are not cardiopulmonary)
    bare_y(axF)
    logx(axF, [1, 2, 4, 8])
    axF.xaxis.set_minor_locator(mticker.NullLocator())
    axF.set_xlim(0.9, 10.0)
    axF.set_ylim(Y_LO, Y_HI)
    axis_covers(axF, TEN.lo2.tolist() + TEN.hi2.tolist() + TEN.or2.tolist(), "x",
                f"{NAME} forest")
    axF.set_xlabel("Odds ratio for staying above 10%\n(95% CI)", labelpad=7,
                   linespacing=1.3)

    # the two printed columns of the global convention, plus their headers. Every q is
    # bold because every drawn condition survives FDR (asserted above).
    BOLD_ARTISTS = []
    trCol = mtransforms.blended_transform_factory(fig.dpi_scale_trans, axF.transData)
    for s, yy in zip(OR_TXT, yD):
        axF.text(COL_X_OR_IN, yy, s, transform=trCol, ha="right", va="center",
                 fontsize=COL_FS, color=INK_P, clip_on=False)
    for s, yy in zip(Q_TXT, yD):
        t = axF.text(COL_X_Q_IN, yy, s, transform=trCol, ha="right", va="center",
                     fontsize=COL_FS, color=INK_P, clip_on=False, fontweight="bold")
        BOLD_ARTISTS.append(t)
    assert len(BOLD_ARTISTS) == 10, len(BOLD_ARTISTS)
    for _x, _h in ((COL_X_OR_IN, COL_HEAD_OR), (COL_X_Q_IN, COL_HEAD_Q)):
        fig.text(_x, AX_Y0 + F_H + COL_HEAD_Y, _h, transform=fig.dpi_scale_trans,
                 ha="right", va="bottom", fontsize=COL_FS, color=INK_P)

    # V14_L5b: the two-entry key drawn in the TOP strip at the V13 positions (V13 page coordinates, the figure sits at
    # (1.70, -0.57) on the 486.51 x 317.07 page, so figure inches = ((x - 1.70) / 72, (320.55 - y) / 72)): 15 pt handles with a
    # 6 pt white-edged marker at their centre, 10 pt ink labels, ArialMT.
    def _fx(x): return (x - 1.70) / 72.0
    def _fy(y): return (320.55 - y) / 72.0
    KEY_SHIFT_1, KEY_SHIFT_2 = -16.0, -8.0   # R49 nudge: entry 1 16 pt and entry 2 8 pt to the left, so the 11 pt labels keep 16 pt to the next handle and 12 pt to the column header (V26 15 and 18 pt at 10 pt; at +1 pt in place they were 8 and 4 pt)
    for (x0, x1, col, tx, lab) in ((136.64 + KEY_SHIFT_1, 151.64 + KEY_SHIFT_1, BLUE, 157.6 + KEY_SHIFT_1, "Cardiopulmonary"), (248.78 + KEY_SHIFT_2, 263.78 + KEY_SHIFT_2, NONCP_COL, 269.8 + KEY_SHIFT_2, "Not cardiopulmonary")):
        fig.add_artist(Line2D([_fx(x0), _fx(x1)], [_fy(18.0), _fy(18.0)], color=col, lw=CI_LW,
                              transform=fig.dpi_scale_trans, solid_capstyle="round"))   # R40: handle line without its end marker (one marker per entry)
        fig.add_artist(Line2D([_fx((x0 + x1) / 2)], [_fy(18.0)], color=col, marker="o", markersize=6, ls="none",
                              markeredgecolor="white", markeredgewidth=MARK_EW, transform=fig.dpi_scale_trans))
        fig.text(_fx(tx), _fy(21.5), lab, transform=fig.dpi_scale_trans, ha="left", va="baseline", fontsize=TICK_PT, color=INK_P)
    KEY = None

    # ---- read back every drawn artist and match it to the frozen values
    got_cis, n_ref = set(), 0
    for ln in axF.lines:
        xd, yd = np.asarray(ln.get_xdata(), float), np.asarray(ln.get_ydata(), float)
        if len(xd) == 2 and xd[0] == xd[1] == 1.0 and ln.is_dashed():
            n_ref += 1
        elif len(xd) == 2 and yd[0] == yd[1]:
            got_cis.add((round(xd[0], 9), round(xd[1], 9), yd[0]))
    assert n_ref == 1, "exactly one reference rule at 1.0"
    assert got_cis == expect_cis, "a confidence interval on the sheet is not the frozen one"
    got_pts = {}
    for c in axF.collections:
        off = np.asarray(c.get_offsets(), float)
        assert off.shape == (1, 2), "each marker is drawn by its own scatter call"
        got_pts[(round(off[0, 0], 9), off[0, 1])] = tuple(np.round(c.get_facecolor()[0][:3], 4))
    assert set(got_pts) == set(expect_pts), "a marker position is not the frozen estimate"
    for k, col in expect_pts.items():
        want = tuple(np.round(plt.matplotlib.colors.to_rgb(col), 4))
        assert got_pts[k] == want, f"marker fill at {k} is off palette (open marker?)"

    DRAWN["conditions"] = [{"condition": r.Condition, "or": r.or2, "lo": r.lo2,
                            "hi": r.hi2, "q": r.q2, "n_with": int(r.n_with),
                            "cardiopulmonary": r.Condition in CARDIOPULM,
                            "printed": EXPECT_OR_TXT[r.Condition],   # V14_L5b: half up on the full-precision twin
                            "printed_q": fmt_p(r.q2), "significant_fdr": True,
                            "marker_filled": True}
                           for _, r in TEN.iterrows()]
    DRAWN["selection"] = {"rule": f"top 10 by marginal OR among the {len(SURV)} FDR survivors, "
                                  "ties by lower bound",
                          "n_sample": N_ALL, "n_failed": N_NO,
                          "largest_or_left_out": float(LEFT_OUT.or2.max()),
                          "left_out_name": REF_ORDER[10], "exact_halves_3dp": EXACT_HALVES, "stored_tie_order_note": ORDER_NOTE, "fullprec_twin": FULLPREC_PATH, "fullprec_twin_sha256": _fp_prov["output_sha256"],
                          "n_fdr_survivors": int(len(SURV)), "n_tested": int(J["n_conditions_tested"])}

    # ---- gates, then write
    for t in fig.findobj(mtext.Text):
        s = t.get_text()
        if not s or not s.strip() or not t.get_visible():
            continue
        for bad, nm in (("—", "em dash"), ("–", "en dash"), (";", "semicolon"),
                        ("×", "multiplication sign")):
            assert bad not in s, f"banned {nm} in on-sheet text: {s!r}"
        if str(t.get_fontweight()) in ("bold", "700"):
            assert any(t is b for b in BOLD_ARTISTS), (
                f"bold outside the sanctioned significant q entries: {s!r}")
            assert s in set(Q_TXT), f"a bold string that is not a q value: {s!r}"
        assert float(t.get_fontsize()) >= FLOOR_PT, (
            f"{t.get_fontsize()} pt under the {FLOOR_PT} pt floor: {s!r}")

    # V14_L5b: the top strip carries the two column headers and the two-entry key (the V13 layout)
    fig.canvas.draw()
    _rk = fig.canvas.get_renderer()
    _ab = axF.get_window_extent(renderer=_rk)
    for t in fig.findobj(mtext.Text):
        s = t.get_text()
        if not s or not s.strip() or not t.get_visible():
            continue
        if t.get_window_extent(renderer=_rk).y1 > _ab.y1:
            assert s in (COL_HEAD_OR, COL_HEAD_Q, "Cardiopulmonary", "Not cardiopulmonary"), (
                f"the top strip carries only the column headers and the key, found {s!r}")

    nooverlap.gate(fig, NAME)
    margin_gate(fig, NAME)
    fig.savefig(f"{POLISH}/{NAME}.pdf")
    fig.savefig(f"{PNGS}/{NAME}.png", dpi=300)
    plt.close(fig)

json.dump(DRAWN, open(f"{JSONS}/{NAME}_polish_drawn_values.json", "w"), indent=1,
          default=float)

print("written:")
for p in (f"{POLISH}/{NAME}.pdf", f"{PNGS}/{NAME}.png",
          f"{JSONS}/{NAME}_polish_drawn_values.json"):
    assert os.path.exists(p) and os.path.getsize(p) > 0, p
    print("  ", p, f"{os.path.getsize(p):,} bytes")
print("drawn, in order:", ", ".join(f"{r.Condition} {_r2(r.or2)}" for _, r in TEN.iterrows()))
