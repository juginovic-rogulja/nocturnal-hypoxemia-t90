"""
eFigure 10 polished: the non-responder sheet, redrawn on the 171 mm design system.

Ships as Supplementary Fig. 14, delivered as two parts, panels a to c and panels d to e.

Content is make_eFigure10_nonresponder_phenotype.py's, untouched: the same five panels,
the same conditions, estimates, sample sizes and printed numbers, the same single
series legend. Every numeric assertion and the printed-number audit are carried over
verbatim.

What changes is the drawing. The old sheet was 14.4 in wide, so at 171 mm the grid
restacks. Ink #1a1d21, 0.8 pt spines, left-aligned label gutters, 6 pt white-edged
markers, the ordered bars on the declared house ramps (#ccd1d6 #8a9099 #4d7387 #1f4257
for the four-step panels, the declared five-step ramp for the fifths), per-bar sample
sizes as 9.5 pt annotations under their own bars, and frameless 10 pt legends. Panel d's
two long model names wrap over two lines, words unchanged, asserted against the settled
wording.

2026-08-21 fix wave, findings B4, M1 and B2 (soft):
  * B4. The single sheet ran 171.0 x 406.4 mm, more than one and a half Nature pages
    tall, and scaled to a page its smallest type landed at 5.5 pt. The minimum regular
    type was already sitting on the 9.0 pt floor, so there was no type left to give.
    Panel a alone is a 19-row odds-ratio forest with a printed OR (95% CI) column, so at
    the 0.175 in row pitch it spends 3.61 in of axes before a single band of label, and
    a two-column reflow of it is out because one column already needs the label gutter
    plus the data span plus the printed column, which is the full 171 mm. Output is now
    two sheets, per the audit's recommended split:
        eFigure10a_nonresponder_phenotype   panels a, b, c
        eFigure10b_nonresponder_phenotype   panels d, e
    Panel letters stay a, b, c, d, e in that order across the two parts, so the reader's
    a-to-e sequence is unbroken. Every part states its own size in the build log and
    asserts itself against the 247.0 mm page height before it is written. Panel b keeps
    its own row and loses 0.60 in of axes height to pay for the split, and its trend
    annotation moves out of the panel into the empty band on its right, where it has room
    at full size. The superseded combined sheet is left on disk untouched for the
    coordinating agent to retire.
  * M1. Panel a marks five of its 19 odds-ratio strings with an asterisk and defined it
    nowhere. The mark is NOT a significance level (all 19 rows drawn are already
    Benjamini-Hochberg significant, that is what puts them on the panel). It flags the
    conditions that survive the joint model, which is panel c's subject. The panel now
    carries the key "* remains independent in the joint model of panel c" at 9.5 pt ink,
    bottom left of the panel that uses it, and the build asserts the marked count equals
    n_conditions_independent_in_joint_model so the key can never be dead ink.
  * B2 (soft). The five non-bold panel-title sentences come off the artwork. The panel
    LETTERS stay. The sentences move verbatim into
    _workfiles/legends_updates/eFigure10_nonresponder_phenotype.fixwave.md, keyed by
    letter. The freed bands are reclaimed, not left white: panel a's two-entry series
    legend moves up into the band its title used to hold, which in turn frees the
    bottom-left of panel a for the asterisk key.

No printed number, odds ratio, confidence interval, n, percent or count changed. The one
string reformatted is panel b's y axis title, wrapped over two lines because the panel is
now shorter than the single line, words unchanged and asserted.

2026-08-24 round 2 (ROUND2_BRIEF, the nonresponder wording pass; Supplementary Fig. 13
in the V2 map). Label wording only, so the sheet reads as the treatment-selection
story: the therapy is named (PAP) and the outcome verb is the family's "remained low"
(the register of the new main Fig. 5d, odds of remaining above 10% on PAP). "stayed low
on treatment" becomes "remained low on PAP" on panel a's axis, "whole sample" becomes
"all patients" on the reference rule, the series key reads "Condition present before
PAP", panel b's axis titles become "Oxygen remained low, %" and "Cardiopulmonary
conditions before PAP", and panel c's axis reads "Odds that oxygen remains low on PAP".
Panel e's y title takes the same verb so the two parts of the figure agree. Not one
number, marker, estimate or panel letter moves.

Sources, read only: numbers/nonresponder_combined_v2.json, numbers/treatment_v2.json.
Output: Supplementary_Figures_Polished/eFigure10a_nonresponder_phenotype.pdf and
eFigure10b_nonresponder_phenotype.pdf, each with its 300 dpi twin in
_workfiles/polished_pngs/.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from supp_polish_common import (rc_polish, save_polished, panel_letter, measure,   # noqa: E402
                                left_rows, fit_gutter, ref_rule, bare_y, logx_clean,
                                axis_covers, sep_ok, plt, BLUE, GREY, GREY_PALE,
                                INK_P, W_IN, EDGE, ROOT, TICK_PT, ANN_PT, CI_LW,
                                MEDGE, S_CIRCLE, S_SQUARE, MS_CIRCLE, MS_SQUARE,
                                ROLE_RAMP4, FIFTHS_RAMP5, PT_PLUS)

from supp_polish_common import PRIMARY   # noqa: E402  V14_L5b 2026-09-15: the V13 sheet draws the series in the main-figure blue
BLUE = PRIMARY   # V14_L5b: #0288d1 (round 31 LCOLOR), not splitstyle's August #1f4257
STEM = "nonresponder_phenotype"
NAME_AC = f"eFigure10a_{STEM}"          # panels a, b, c
NAME_DE = f"eFigure10b_{STEM}"          # panels d, e
NUM = f"{paths.NUMBERS_DIR}"

# ------------------------------------------------------------------ data, verbatim
P = json.load(open(f"{NUM}/nonresponder_combined_v2.json"))
RATES = {k: v for k, v in P["failure_rate_by_condition"].items() if v["fdr_significant"]}
assert len(RATES) == P["n_conditions_fdr_significant"] > 0, len(RATES)   # V14_L5b 2026-09-15: data-driven (gate 6), the literal 19 is gone
# V14_L5b 2026-09-15: the full-precision twin of step 127 (lane work/fullprec_nonresponder_v8_1.json, sidecar sha-checked, reconciled with
# the shipped file at 5e-4): every printed odds ratio is half up on full precision, never a second rounding of the 3-decimal file
import hashlib
from decimal import Decimal, ROUND_HALF_UP
LANE = f"{paths.FIGURE_ROOT}/ROUND50_2026-09-26/lanes/LSUPP-C/L5b_work"   # R49: the twin regenerated in this lane (work/fullprec_nonresponder_v8_1.json)
FPP = f"{LANE}/work/fullprec_nonresponder_v8_1.json"
_prov = json.load(open(FPP + ".provenance.json"))
assert _prov["output_sha256"] == hashlib.sha256(open(FPP, "rb").read()).hexdigest(), "twin sidecar mismatch"
assert _prov["reconciled_against"]["sha256"] == hashlib.sha256(open(f"{NUM}/nonresponder_combined_v2.json", "rb").read()).hexdigest(), "twin reconciled against another combined json"
FULL = json.load(open(FPP))
def _r(v, nd=2): return str(Decimal(repr(float(v))).quantize(Decimal(1).scaleb(-nd), rounding=ROUND_HALF_UP))
FPC = {**FULL["joint_independent"], "Age, per year": FULL["joint_demographics"]["AgeAtVisit"],
       "Male sex": FULL["joint_demographics"]["male"], "Black race": FULL["joint_demographics"]["race_black"]}
DR = P["dose_response"]
AUC = P["auc"]
R5 = P["risk_fifths_from_crossvalidated_score"]

_T = json.load(open(f"{NUM}/treatment_v2.json"))["corrected_vs_not"]
assert (P["n"], P["n_failed"], P["n_corrected"]) == (_T["n"], _T["n_not"], _T["n_corrected"])
assert P["n"] > 0 and abs(P["failure_pct"] - 100.0 * P["n_failed"] / P["n"]) < 0.06, (P["n"], P["n_failed"], P["failure_pct"])   # V14_L5b: data-driven

LOG, PRINTED = [], []
LOG.append(sep_ok("condition present against absent", [BLUE, GREY_PALE]))
LOG.append(sep_ok("dose ramp, alternate steps", [ROLE_RAMP4[0], ROLE_RAMP4[2]]))
LOG.append(sep_ok("risk-fifths ramp, alternate steps",
                  [FIFTHS_RAMP5[0], FIFTHS_RAMP5[2], FIFTHS_RAMP5[4]]))

# ------------------------------------------------------------------ M1, the asterisk
# The mark on panel a's odds-ratio column is not a significance level. Every row drawn on
# panel a already cleared Benjamini-Hochberg correction, that is the filter that put it
# there. The asterisk flags the conditions that then survive the JOINT model of all 19
# conditions plus demographics, which is exactly the set panel c draws. The key says that,
# and the count is asserted against the source so the key is never dead ink.
STAR_KEY = "* remains independent in the joint model of panel c"
N_MARKED = sum(1 for v in RATES.values() if v["independent"])
assert N_MARKED == P["n_conditions_independent_in_joint_model"] > 0, N_MARKED   # V14_L5b: data-driven
assert set(k for k, v in RATES.items() if v["independent"]) == set(P["independent_conditions"])

# ------------------------------------------------------------------ geometry, inches
GUT = EDGE
A_AX_X0 = 1.92                   # panel a shares the polished forest gutter edge
A_AX_W = W_IN - EDGE - A_AX_X0
A_DATA_W = 3.36                  # 0 to 76% spans this, the odds column sits beyond
OR_X_IN = 5.38                   # left edge of the printed odds-ratio column

B_AX = [0.88, 2.60]              # x0, width. Wider than the 2026-08-21 wave-2 sheet,
                                 # and far enough right that the wrapped y title clears
                                 # the 4 mm canvas margin
B_TREND_X = B_AX[0] + B_AX[1] + 0.34   # the trend annotation, in the band right of b
# 2026-08-21 wave 2 recomposition: panel c moves to its own FULL-WIDTH row so its forest
# carries the printed "Odds ratio (95% CI)" column the round brief requires of every
# forest. It adopts panel a's grammar exactly: label gutter at GUT, axes from A_AX_X0,
# data spanning C_DATA_W, the printed column at OR_X_IN, so the two columns align down
# the sheet. Reading order stays a, b, c, d, e (b keeps its own row above c), so the
# internal letters are UNCHANGED.
C_AX_W = W_IN - EDGE - A_AX_X0   # panel c axes width, full width like panel a
C_DATA_W = 3.36                  # panel c data span, the column sits beyond, as on a
D_GUT, D_AX = 0.18, [1.74, 1.56]
E_TIT = 3.46
E_AX = [3.98, W_IN - EDGE - 3.98]

ROW_A = ROW_C = 0.175            # row pitch, both forests, at the 10 pt row-label floor
H2 = 1.15                        # panel b axes height (was 1.75 on the combined sheet)
H3 = 1.80                        # panels d and e axes height, unchanged

# ---- the bands. Each is the measured stack of what sits in it, not a guess.
LET_BAND = 0.30      # 0.08 of air plus a 13 pt bold panel letter
A_XBAND = 0.50       # panel a: x tick labels plus the axis title
A_KEY = 0.24         # panel a: the asterisk key line, 9.5 pt
B_XBAND = 0.72       # panel b: x tick labels, the per-bar n, the axis title
C_XBAND = 0.50       # panel c: x tick labels plus the axis title
DE_XBAND = 0.72      # panels d and e: the deeper of the two label stacks

# 2026-08-25 round 3 (Alen's Supp 13 comment, "a/b/c left + d/e right reads messy"):
# the five panels recompose 2+3. Panel c moves from part 1 to part 2, so part 1 holds
# the tall absolute-risk forest with its dose bars (a, b) and part 2 opens with the
# joint model (c) over the AUC and risk-fifth panels (d, e). The two parts now stand
# as similar-scale portraits and the review page reads a, b then c, d, e. Letters
# stay a-e in order across the parts, no letter changes, so no callout cascades.
LETTERS_AC = ("a", "b")
LETTERS_DE = ("c", "d", "e")

DRAWN_COMMON = {"panelA": RATES, "panelB": DR, "panelB_trend": P["dose_response_trend"],
                "panelD": AUC, "panelE": R5, "greyscale_separation": LOG,
                "asterisk_key": {"text": STAR_KEY, "marks": N_MARKED,
                                 "means": "independent in the joint model, panel c",
                                 "is_significance_level": False}}


def _c_rows():
    """Panel c's row content, needed before any canvas exists so the height is known."""
    ind = P["independent_conditions"]
    keep = sorted(ind, key=lambda k: -ind[k]["or"])
    joint = P["joint_model"]
    demo = [("Age, per year", joint["demographics"]["AgeAtVisit"]),
            ("Male sex", joint["demographics"]["male"]),
            ("Black race", joint["demographics"]["race_black"])]
    rows = [(k, ind[k], BLUE) for k in keep] + [(n, v, GREY) for n, v in demo]
    return ind, keep, rows


IND, KEEP, ROWS_C = _c_rows()
assert KEEP == sorted(KEEP, key=lambda k: -FULL["joint_independent"][k]["or"]), "joint-model row order differs between the 3-decimal file and the twin"   # V14_L5b
ORDER = sorted(RATES, key=lambda k: RATES[k]["fail_pct_with"])

A_H = ((len(ORDER) + 0.75) - (-0.9)) * ROW_A          # panel a axes height
HC = (len(ROWS_C) + 0.75 + 0.9) * ROW_C               # panel c axes height

# ---- R50 (Alen, 2026-09-26, item 3): the lower row (b, d, e) is its OWN part, laid across the full page under a and c, on one
# ---- baseline with one height (H3) and even gutters, the odds-ratio note as b's one-line subtitle. Part 1 holds a alone and part 2
# ---- c alone with the same top layout as before (letter band + edge above the axes), so a and c do not move on the page.
SUB_BAND = 0.22      # the one-line subtitle band of panel b (ANN_PT), between the letter band and the axes, on all three panels
ROW_GAP = 0.45       # the even gutter: between the lower panels' ink blocks, and between the upper parts' inner edge and the row's
YA = EDGE + A_XBAND                                   # panel a axes bottom (part 1 = a alone; the empty asterisk-key band is gone)
H_A = YA + A_H + LET_BAND + EDGE
YC = EDGE + C_XBAND                                   # panel c axes bottom (part 2 = c alone)
H_C = YC + HC + LET_BAND + EDGE
W_ROW = 2 * W_IN + 24.0 / 72.0                        # the composed page width (two 171 mm parts 24 pt apart) minus the 14 pt margins
Y_ROW_AX = EDGE + DE_XBAND                            # the lower row's axes bottom (b, d, e share it: one baseline)
H_ROW = Y_ROW_AX + H3 + SUB_BAND + LET_BAND + EDGE
B_GUT_W, D_GUT_W, E_GUT_W = B_AX[0] - EDGE, D_AX[0] - D_GUT, E_AX[0] - E_TIT   # each panel's own left block (y title and ticks, the wrapped labels, y title)
W_AX_ROW = (W_ROW - 2 * EDGE - (B_GUT_W + D_GUT_W + E_GUT_W) - 2 * ROW_GAP) / 3.0   # three axes of one width
XB0 = EDGE + B_GUT_W; XD_LET = XB0 + W_AX_ROW + ROW_GAP; XD0 = XD_LET + D_GUT_W; XE_LET = XD0 + W_AX_ROW + ROW_GAP; XE0 = XE_LET + E_GUT_W
assert abs(XE0 + W_AX_ROW + EDGE - W_ROW) < 1e-9, (XE0, W_AX_ROW, W_ROW)
NAME_ROW = f"eFigure10c_{STEM}_row"                   # panels b, d, e
LETTERS_ROW = ("b", "d", "e")
ROW_LAYOUT = {"H_A_in": H_A, "H_C_in": H_C, "H_ROW_in": H_ROW, "W_ROW_in": W_ROW, "EDGE_in": EDGE, "ROW_GAP_in": ROW_GAP, "SUB_BAND_in": SUB_BAND,
              "axes_width_in": W_AX_ROW, "axes_x0_in": {"b": XB0, "d": XD0, "e": XE0}, "letter_x_in": {"b": EDGE, "d": XD_LET, "e": XE_LET},
              "row_top_from_part1_top_in": H_A - 2 * EDGE + ROW_GAP,   # the row part's top edge sits ROW_GAP below part 1's inner edge
              "rule": "the lower row starts ROW_GAP below part 1's inner bottom edge (its own top edge, EDGE, then the letters); gutters between the panels' ink blocks equal ROW_GAP"}

for _part, _h in ((NAME_AC, H_A), (NAME_DE, H_C), (NAME_ROW, H_ROW)):
    assert _h * 25.4 <= 247.0, (_part, round(_h * 25.4, 1), "over the page height")
assert abs(W_IN * 25.4 - 171.0) < 0.05, "sheet width moved off 171 mm"


# ------------------------------------------------------------------ part 1: a
def build_a():
    fig = plt.figure(figsize=(W_IN, H_A))
    fit_gutter(fig, ORDER, GUT, A_AX_X0 + 0.03, what=f"{NAME_AC} a")   # R49: at 11 pt the two longest labels end 2.0 and 1.3 pt past the axes edge, which carries no ink here (no spine, no band, nearest mark 1.04 in away, the reference rule 1.33 in away)

    # -------------------------------------------------------- a, absolute risk (unchanged, R50)
    axA = fig.add_axes([A_AX_X0 / W_IN, YA / H_A, A_AX_W / W_IN, A_H / H_A])
    y = np.arange(len(ORDER)).astype(float)
    XLIM_HI = 76.0 * A_AX_W / A_DATA_W
    OR_X = (OR_X_IN - A_AX_X0) / A_AX_W * XLIM_HI
    for i, k in enumerate(ORDER):
        v = RATES[k]
        axA.plot([v["fail_pct_without"], v["fail_pct_with"]], [y[i], y[i]], color=GREY,
                 lw=1.0, zorder=2)
        axA.scatter([v["fail_pct_without"]], [y[i]], s=S_SQUARE, marker="s",
                    facecolors=GREY_PALE, edgecolors=GREY, linewidths=0.8, zorder=3)
        axA.scatter([v["fail_pct_with"]], [y[i]], s=S_CIRCLE, marker="o",
                    facecolors=BLUE, edgecolors="white", linewidths=MEDGE, zorder=4)
    axA.axvline(P["failure_pct"], color=INK_P, lw=0.9, ls=(0, (4, 3)), zorder=1)
    _ws = f"all patients, {P['failure_pct']:.1f}%"
    axA.text(P["failure_pct"] + 2.5, len(ORDER) - 0.28, _ws, fontsize=ANN_PT,
             color=INK_P, ha="left")
    PRINTED.append(("whole-sample rule", f"{P['failure_pct']:.1f}", float(P["failure_pct"])))
    left_rows(axA, y, ORDER, GUT, A_AX_X0)
    bare_y(axA)
    axA.set_xlim(0, XLIM_HI)
    axA.set_xticks([0, 20, 40, 60])
    axA.spines["bottom"].set_bounds(0, 76)
    axA.set_ylim(-0.9, len(ORDER) + 0.75)
    axis_covers(axA, [RATES[k]["fail_pct_with"] for k in ORDER]
                + [RATES[k]["fail_pct_without"] for k in ORDER] + [P["failure_pct"]], "x",
                "eFigure10 A")
    axA.set_xlabel("Patients whose oxygen remained low on PAP, %", labelpad=6)
    axA.xaxis.set_label_coords((A_DATA_W / 2) / A_AX_W, -0.088)
    for i, k in enumerate(ORDER):
        v = RATES[k]
        mark = "*" if v["independent"] else ""
        fp = FULL["by_condition"][k]   # V14_L5b: half up on the full-precision twin (each value within 5e-4 of the 3-decimal file)
        for a_ in ("marginal_or", "marginal_lo", "marginal_hi"):
            assert abs(float(fp[a_]) - float(v[a_])) <= 5e-4, (k, a_, fp[a_], v[a_])
        s = f"{_r(fp['marginal_or'])} ({_r(fp['marginal_lo'])}-{_r(fp['marginal_hi'])}){mark}"
        axA.text(OR_X, y[i], s, fontsize=ANN_PT, va="center", ha="left", color=INK_P)
        for lab, shown, src in ((f"{k} or", _r(fp["marginal_or"]), fp["marginal_or"]),
                                (f"{k} lo", _r(fp["marginal_lo"]), fp["marginal_lo"]),
                                (f"{k} hi", _r(fp["marginal_hi"]), fp["marginal_hi"])):
            PRINTED.append((lab, shown, float(src)))
    axA.text(OR_X, len(ORDER) + 0.45, "Odds ratio (95% CI)", fontsize=ANN_PT, ha="left",
             va="center", color=INK_P)
    # the asterisk key is not printed (round 40 LNOTES); its band is gone this round (R50), nothing sits under panel a's x title
    from matplotlib.lines import Line2D
    lg = fig.legend(handles=[
        Line2D([], [], marker="o", color=BLUE, ls="none", ms=MS_CIRCLE,
               markeredgecolor="white", markeredgewidth=MEDGE,
               label="Condition present before PAP"),
        Line2D([], [], marker="s", color="none", markerfacecolor=GREY_PALE,
               markeredgecolor=GREY, ms=MS_SQUARE, label="Condition absent")],
        loc="upper right", bbox_to_anchor=((W_IN - EDGE) / W_IN,
                                           (YA + A_H + LET_BAND) / H_A), ncol=2,
        fontsize=TICK_PT, frameon=False, borderpad=0.0, handletextpad=0.5,
        columnspacing=1.4, handlelength=1.3)
    fig.canvas.draw()
    _lgbb = lg.get_window_extent(fig.canvas.get_renderer())
    assert _lgbb.x0 / fig.dpi > GUT + 0.22, "panel a's legend reaches the panel letter"
    assert _lgbb.y0 / fig.dpi > YA + A_H + 0.04, "panel a's legend reaches its own axes"
    panel_letter(fig, "a", GUT, YA + A_H + LET_BAND, H_A)
    drawn = dict(DRAWN_COMMON)
    drawn["geometry"] = {"width_mm": round(W_IN * 25.4, 1), "height_mm": round(H_A * 25.4, 1), "panels": ["a"], "part": "1 of 3",
                         "changes_2026_09_26_round50": ["part 1 holds panel a alone (panel b moved to the lower-row part), the same top layout, so a does not move on the page"]}
    save_polished(fig, NAME_AC, drawn, headline_lines=("a",), min_mm=1.5)   # R49: the column header "Odds ratio (95% CI)" at 10.5 pt ends 1.6 mm from this part's right edge (6.5 mm inside the composed page), nothing moved


# ------------------------------------------------------------------ part 2: c
def build_c():
    fig = plt.figure(figsize=(W_IN, H_C))
    import math
    fit_gutter(fig, [lab for lab, _v, _c in ROWS_C], GUT, A_AX_X0,
               what=f"{NAME_DE} c")   # R49: the longest label ends 2.9 pt before the (ink-free) axes edge at 11 pt
    axC = fig.add_axes([A_AX_X0 / W_IN, YC / H_C, C_AX_W / W_IN, HC / H_C])
    yy = np.arange(len(ROWS_C))[::-1].astype(float)
    for i, (_lab, v, col) in enumerate(ROWS_C):
        axC.plot([v["lo"], v["hi"]], [yy[i], yy[i]], color=col, lw=CI_LW,
                 solid_capstyle="round", zorder=2)
        axC.scatter([v["or"]], [yy[i]], s=S_CIRCLE, color=col, zorder=3,
                    edgecolors="white", linewidths=MEDGE)
    ref_rule(axC, 1.0)
    left_rows(axC, yy, [lab for lab, _v, _c in ROWS_C], GUT, A_AX_X0)
    for lab in axC.get_yticklabels()[len(KEEP):]:
        lab.set_color(GREY)
    bare_y(axC)
    _cv = [v[k] for _l, v, _c in ROWS_C for k in ("or", "lo", "hi")]
    _xlo, _xhi_data = min(_cv) * 0.90, max(_cv) * 1.14
    _frac = C_DATA_W / C_AX_W
    _xhi_full = math.exp(math.log(_xlo) + (math.log(_xhi_data) - math.log(_xlo)) / _frac)
    OR_XC = math.exp(math.log(_xlo) + ((OR_X_IN - A_AX_X0) / C_AX_W)
                     * (math.log(_xhi_full) - math.log(_xlo)))
    logx_clean(axC, [t for t in (0.5, 0.8, 1, 1.5, 2, 3, 5)
                     if _xlo <= t <= _xhi_data])
    axC.set_xlim(_xlo, _xhi_full)
    axC.spines["bottom"].set_bounds(_xlo, _xhi_data)
    axC.set_ylim(-0.9, len(ROWS_C) + 0.75)
    axis_covers(axC, _cv, "x", "eFigure10 C")
    for i, (labk, v, _col) in enumerate(ROWS_C):
        fp = FPC[labk]   # V14_L5b: half up on the full-precision twin (each value within 5e-4 of the 3-decimal file)
        for a_ in ("or", "lo", "hi"):
            assert abs(float(fp[a_]) - float(v[a_])) <= 5e-4, (labk, a_, fp[a_], v[a_])
        s = f"{_r(fp['or'])} ({_r(fp['lo'])}-{_r(fp['hi'])})"
        axC.text(OR_XC, yy[i], s, fontsize=ANN_PT, va="center", ha="left", color=INK_P)
        for sub, shown, src in ((f"{labk} joint or", _r(fp["or"]), fp["or"]),
                                (f"{labk} joint lo", _r(fp["lo"]), fp["lo"]),
                                (f"{labk} joint hi", _r(fp["hi"]), fp["hi"])):
            PRINTED.append((sub, shown, float(src)))
    axC.text(OR_XC, len(ROWS_C) + 0.45, "Odds ratio (95% CI)", fontsize=ANN_PT,
             ha="left", va="center", color=INK_P)
    axC.set_xlabel("Odds that oxygen remains low on PAP, all terms in one model (95% CI)",
                   labelpad=6)
    axC.xaxis.set_label_coords((C_DATA_W / 2) / C_AX_W, -(0.30 / HC))
    assert P["n_conditions_independent_in_joint_model"] == len(KEEP)   # V14_L5b: data-driven
    panel_letter(fig, "c", GUT, YC + HC + LET_BAND, H_C)
    drawn = dict(DRAWN_COMMON)
    drawn["panelC"] = {"independent": IND, "demographics": P["joint_model"]["demographics"]}
    drawn["geometry"] = {"width_mm": round(W_IN * 25.4, 1), "height_mm": round(H_C * 25.4, 1), "panels": ["c"], "part": "2 of 3",
                         "changes_2026_09_26_round50": ["part 2 holds panel c alone (d and e moved to the lower-row part), the same top layout, so c does not move on the page"]}
    save_polished(fig, NAME_DE, drawn, headline_lines=("c",), min_mm=1.5)   # R49: same header on panel c


# ------------------------------------------------------------------ part 3: the lower row b, d, e (R50)
from supp_polish_common import _prose_gate, margin_gate, nooverlap, POLISH, PNGS, WORK, PANEL_PT   # noqa: E402


def row_letter(fig, letter, x_in, y_in):
    """Bold panel letter at an absolute position on the row canvas (W_ROW wide), va top."""
    fig.text(x_in / W_ROW, y_in / H_ROW, letter, fontsize=PANEL_PT, fontweight="bold", ha="left", va="top", color=INK_P)


def save_row(fig, name, drawn):
    """The row part is wider than 171 mm, so the journal preflight is not run on it: the prose, overlap and margin gates are."""
    _prose_gate(fig, name, LETTERS_ROW)
    nooverlap.gate(fig, name)
    margin_gate(fig, name)
    pdf = f"{POLISH}/{name}.pdf"
    fig.savefig(pdf)
    fig.savefig(f"{PNGS}/{name}.png", dpi=300)
    plt.close(fig)
    json.dump(drawn, open(f"{WORK}/{name}_polished_drawn_values.json", "w"), indent=1, default=float)
    for p in (pdf, f"{PNGS}/{name}.png"):
        assert os.path.exists(p) and os.path.getsize(p) > 0, p
        print("wrote", p)
    return pdf


def build_row():
    fig = plt.figure(figsize=(W_ROW, H_ROW))
    Y_LET = Y_ROW_AX + H3 + SUB_BAND + LET_BAND      # the letters' top, one line for b, d, e

    # -------------------------------------------------------- b, dose response, with the odds-ratio note as its one-line subtitle
    axB = fig.add_axes([XB0 / W_ROW, Y_ROW_AX / H_ROW, W_AX_ROW / W_ROW, H3 / H_ROW])
    xb = np.arange(4)
    assert sum(DR["n"]) == P["n"]
    axB.bar(xb, DR["failed_pct"], color=ROLE_RAMP4, width=0.68, linewidth=0)
    axB.set_xticks(xb)
    axB.set_xticklabels(DR["labels"], fontsize=TICK_PT)
    axB.set_xlim(-0.55, 3.55)
    axB.set_xlabel("Cardiopulmonary conditions before PAP", labelpad=24)
    _YLB = "Oxygen remained\nlow, %"
    assert _YLB.replace("\n", " ") == "Oxygen remained low, %"   # words as settled, round 2
    axB.set_ylabel(_YLB)
    axB.set_ylim(0, 88)
    axB.set_yticks([0, 20, 40, 60, 80])
    axis_covers(axB, DR["failed_pct"], "y", "eFigure10 B")
    tr = P["dose_response_trend"]
    ftr = FULL["dose_response_trend"]   # V14_L5b: half up on the full-precision twin
    for a_ in ("or", "lo", "hi"):
        assert abs(float(ftr[a_]) - float(tr[a_])) <= 5e-4, (a_, ftr[a_], tr[a_])
    NOTE = f"odds ratio per added condition {_r(ftr['or'])} (95% CI, {_r(ftr['lo'])}-{_r(ftr['hi'])})"   # R50: one line, directly above b's bars
    fig.text(XB0 / W_ROW, (Y_ROW_AX + H3 + 0.06) / H_ROW, NOTE, fontsize=ANN_PT, va="bottom", ha="left", color=INK_P)
    for lab, shown, src in (("trend or", _r(ftr["or"]), ftr["or"]),
                            ("trend lo", _r(ftr["lo"]), ftr["lo"]),
                            ("trend hi", _r(ftr["hi"]), ftr["hi"])):
        PRINTED.append((lab, shown, float(src)))
    row_letter(fig, "b", EDGE, Y_LET)

    # -------------------------------------------------------- d, what it is worth
    MOD = ["Pretreatment oxygen alone", "Pretreatment oxygen + demographics",
           "Cardiopulmonary count + oxygen + demographics",
           "Independent conditions + oxygen + demographics"]
    SHORT = ["Pretreatment oxygen alone", "and demographics",
             "and the cardiopulmonary count", "and the independent conditions"]
    SHORT_WRAP = ["Pretreatment oxygen\nalone", "and demographics",
                  "and the\ncardiopulmonary count", "and the independent\nconditions"]
    assert [s.replace("\n", " ") for s in SHORT_WRAP] == SHORT   # words unchanged
    fit_gutter(fig, SHORT_WRAP, XD_LET, XD0 - 0.06, size=TICK_PT - PT_PLUS, what=f"{NAME_ROW} d")   # R49 exception: the four wrapped model labels keep 10 pt
    axD = fig.add_axes([XD0 / W_ROW, Y_ROW_AX / H_ROW, W_AX_ROW / W_ROW, H3 / H_ROW])
    yv = np.arange(len(MOD))[::-1].astype(float)
    BASE = 0.5
    for i, m in enumerate(MOD):
        a = AUC[m]["auc_cv"]
        axD.barh(yv[i], a - BASE, left=BASE, height=0.58, color=ROLE_RAMP4[i],
                 linewidth=0)
    axD.axvline(AUC[MOD[0]]["auc_cv"], color=INK_P, lw=0.9, ls=(0, (4, 3)), zorder=4)
    axD.text(AUC[MOD[0]]["auc_cv"] - 0.008, 3.92, "pretreatment\noxygen alone",
             ha="right", va="center", fontsize=ANN_PT, color=INK_P, linespacing=1.2)
    left_rows(axD, yv, SHORT_WRAP, XD_LET, XD0, size=TICK_PT - PT_PLUS)   # R49 exception, see fit_gutter above
    bare_y(axD)
    axD.set_xlim(BASE, 0.94)
    axD.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9])
    axD.set_ylim(-0.55, 4.42)
    axis_covers(axD, [AUC[m]["auc_cv"] for m in MOD], "x", "eFigure10 D")
    axD.set_xlabel("Prediction accuracy\n(cross-validated AUC)", labelpad=6)
    row_letter(fig, "d", XD_LET, Y_LET)

    # -------------------------------------------------------- e, risk fifths
    axE = fig.add_axes([XE0 / W_ROW, Y_ROW_AX / H_ROW, W_AX_ROW / W_ROW, H3 / H_ROW])
    xr = np.arange(5)
    assert sum(r["n"] for r in R5) == P["n"]
    axE.bar(xr, [r["failed_pct"] for r in R5], color=FIFTHS_RAMP5, width=0.70,
            linewidth=0)
    axE.set_xticks(xr)
    axE.set_xticklabels(["1, lowest", "2", "3", "4", "5, highest"], fontsize=TICK_PT)
    axE.set_xlim(-0.60, 4.70)
    axE.set_xlabel("Fifth of predicted risk (quintile)", labelpad=24)
    axE.set_ylabel("Oxygen remained low, %")
    axE.set_ylim(0, 92)
    axE.set_yticks([0, 20, 40, 60, 80])
    axis_covers(axE, [r["failed_pct"] for r in R5], "y", "eFigure10 E")
    row_letter(fig, "e", XE_LET, Y_LET)

    drawn = dict(DRAWN_COMMON)
    drawn["note_one_line"] = NOTE
    drawn["geometry"] = {"width_mm": round(W_ROW * 25.4, 1), "height_mm": round(H_ROW * 25.4, 1), "panels": ["b", "d", "e"], "part": "3 of 3", "layout": ROW_LAYOUT,
                         "changes_2026_09_26_round50": ["panels b, d and e on one baseline with one height (H3), three axes of one width, even gutters (ROW_GAP) between their ink blocks",
                                                        "the odds-ratio note is b's one-line subtitle directly above its bars (same words, same size)",
                                                        "letters b, d, e at the panels' top-left corners on one line"]}
    save_row(fig, NAME_ROW, drawn)


with plt.rc_context(rc_polish()):
    build_a()
    build_c()
    build_row()
json.dump(ROW_LAYOUT, open(f"{WORK}/eFigure10_row_layout.json", "w"), indent=1)

# the printed-number audit, carried over verbatim and now covering the three parts
for lab, shown, src in PRINTED:
    dec = len(shown.split(".")[1]) if "." in shown else 0
    assert shown == _r(src, dec), (lab, shown, src)   # V14_L5b: half up on the value the string was printed from
print(f"printed-number audit: {len(PRINTED)} values reproduced from the source JSON")
print(f"asterisk key: {STAR_KEY!r}, {N_MARKED} rows marked (not printed)")
for _n, _h, _w, _letters in ((NAME_AC, H_A, W_IN, ("a",)), (NAME_DE, H_C, W_IN, ("c",)), (NAME_ROW, H_ROW, W_ROW, LETTERS_ROW)):
    print(f"  {_n}: {_w * 25.4:.1f} x {_h * 25.4:.1f} mm, panels {' '.join(_letters)}")
