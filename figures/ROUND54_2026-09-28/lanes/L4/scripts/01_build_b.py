#!/usr/bin/env python3
"""ROUND 49 (2026-09-26, lane L4): copy of the round-38 builder with every text one point larger through the size constants (DPT = 1.0):
row labels, block header and axis title 11 -> 12, legend and ticks 10 -> 11. EXCEPTION kept at 10 pt: the printed HR (95% CI) and q
columns with their headers (COL_DPT = 0.0, the column gap rules fail at 11 pt at the pinned edges). Nudged rules: the line breaks of
the row labels frozen at the V26 breaks, the gutter margin 0.05 -> 0.02 in, the axis-title allowance 0.30 -> 0.50 in. Round-38 docstring follows.
V14 lane L4 (round 37, 2026-09-15), Main Fig 4 panel b: the within-group paired forest, NO OR MILD apnea (AHI <15) beside
SEVERE apnea (AHI 30 or more), T90 above 10% against 10% or less inside each group, from the v8.1 numbers.
Copy of ROUND30_2026-09-08/V7_L4_APNEA/Main_Fig4/scripts/01_build_b.py (beside as 01_build_b_PRE_V8_1.py), repointed:
  numbers gated by the v8.1 sidecar (l4lib.sidecar = v14lib.sidecar_v8), no literal cohort counts (the file's cohort_n must equal
  crosstab_v2's n, the cohort n cohorts.json's), the negative-control panel read from numbers/negcontrols_final.json (controls the
  file could not fit in an arm are ABSENT and reported, not asserted away), the row set by the file's own lists (the round-30 rule),
  printed strings half up on the file's values (3 dp stored: the exact-half strings are listed in the drawn record), the V13 design:
  no pale vertical grid lines, no grey ties between the paired markers (rounds 31 to 36), the V13 panel box and plot box kept
  (the round-30 rule: the row pitch is what the rows need to fill the V13 box, so the axis rule stays on panel a's rule).
Outputs: work/panel_b_raw.pdf, work/panel_b_geometry.json, verify/Main_Fig4_b_drawn.json (DRAWN_RECORD fields per printed string).
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import json, os, re, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import l4lib as L
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND37_2026-09-14/figures/_common"); import v14lib as V

SHEET = "Main_Fig4"; SD = f"{L.LANE}/{SHEET}"; WORK = f"{SD}/work"; VER = f"{SD}/verify"
os.makedirs(WORK, exist_ok=True); os.makedirs(VER, exist_ok=True)
matplotlib, plt = L.mpl_setup()
from matplotlib.colors import to_rgba
from matplotlib.legend_handler import HandlerTuple
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
from matplotlib.transforms import blended_transform_factory

# ------------------------------------------------------------------ the numbers files, gated (v8.1)
SRC = f"{L.NUM}/fig3_within_stratum_v1.json"; SC = L.sidecar(SRC); assert SC["step"]["rc"] == 0, SC["step"]
XT_PATH = f"{L.NUM}/crosstab_v2.json"; SCX = L.sidecar(XT_PATH)
COH_PATH = f"{L.NUM}/cohorts.json"; SCC = L.sidecar(COH_PATH); COH = json.load(open(L.hydrated(COH_PATH)))
NC_PATH = f"{L.NUM}/negcontrols_final.json"; SCN = L.sidecar(NC_PATH); NCF = json.load(open(L.hydrated(NC_PATH)))
WS = json.load(open(L.hydrated(SRC))); FITS = WS["fits"]; ST = WS["strata"]; XT = json.load(open(L.hydrated(XT_PATH)))
N_COHORT = int(COH["bdsp"]["n_analysis"])
assert WS["provenance"]["cohort_n"] == XT["n"] and XT["n_cohort"] == N_COHORT and XT["n"] + XT["n_ahi_undefined"] == N_COHORT, (WS["provenance"]["cohort_n"], XT["n"], XT["n_cohort"], N_COHORT)
assert WS["provenance"]["contrast"].startswith("T90 above 10%"), WS["provenance"]["contrast"]
assert WS["drawn_left_group"] == "no_or_mild", WS["drawn_left_group"]
LEFT, RIGHT = "no_or_mild", "severe"
assert ST[LEFT]["rule"] == "AHI < 15" and ST[RIGHT]["rule"] == "AHI >= 30", (ST[LEFT], ST[RIGHT])
for g in (LEFT, RIGHT): assert ST[g]["n"] == ST[g]["n_low_o2"] + ST[g]["n_rest"], ST[g]
r3 = {r["category"]: r for r in XT["rows3"]}
assert ST[LEFT]["n"] == r3["No or mild (AHI <15)"]["n"] and ST[LEFT]["n_low_o2"] == r3["No or mild (AHI <15)"]["n_>10%"], (ST[LEFT], r3)
assert ST[RIGHT]["n"] == r3["Severe (AHI >=30)"]["n"] and ST[RIGHT]["n_low_o2"] == r3["Severe (AHI >=30)"]["n_>10%"], (ST[RIGHT], r3)

C_CONDS = list(WS["panel_c_rows_no_or_mild"]); C_CTRL = list(WS["panel_c_controls_no_or_mild"]); C_CTRL_ABSENT = list(WS["panel_c_controls_absent_no_or_mild"])
NEG_CONTROLS = list(NCF["panel_labels"]); assert len(NEG_CONTROLS) == NCF["n_controls"] == 5, NEG_CONTROLS
assert set(C_CTRL) | set(C_CTRL_ABSENT) == set(NEG_CONTROLS), (C_CTRL, C_CTRL_ABSENT, NEG_CONTROLS)
assert not set(C_CONDS) & set(NEG_CONTROLS) and len(C_CONDS) == len(set(C_CONDS))
QTHR = 0.05
for c in C_CONDS:
    r = FITS[c]; assert not r[LEFT]["unstable"] and not r[RIGHT]["unstable"], c
    assert r[LEFT]["q"] < QTHR, (c, r[LEFT]["q"])                      # the row set is anchored on the left group
    assert r[LEFT]["q_pool"] == r[RIGHT]["q_pool"] == "conditions", c
    assert r[LEFT]["hr"] > 1, (c, r[LEFT]["hr"])
for c in C_CTRL:
    r = FITS[c]; assert r[LEFT]["q_pool"] == r[RIGHT]["q_pool"] == "negative controls", c
    assert r[LEFT]["q"] >= QTHR and r[RIGHT]["q"] >= QTHR, (c, r[LEFT]["q"], r[RIGHT]["q"])
for c in C_CTRL_ABSENT:
    r = FITS[c]; assert r[LEFT]["unstable"] or r[RIGHT]["unstable"] or r[LEFT].get("q") is None, (c, "absent control is not unstable")
LOW_SUPPORT = [(c, a) for c in C_CONDS + C_CTRL for a in (LEFT, RIGHT) if FITS[c][a]["low_support"]]
assert sorted(LOW_SUPPORT) == sorted((x[0], x[1]) for x in WS["low_support_rows_no_or_mild"]), (LOW_SUPPORT, WS["low_support_rows_no_or_mild"])
hrs_left = [FITS[c][LEFT]["hr"] for c in C_CONDS]
ORDER_NOTE = "file order (panel_c_rows_no_or_mild), descending left-group HR" if hrs_left == sorted(hrs_left, reverse=True) else "file order (panel_c_rows_no_or_mild), NOT monotone in the left-group HR"


def arms(cond):
    a, s = FITS[cond][LEFT], FITS[cond][RIGHT]
    return ((a["hr"], a["lo"], a["hi"], a["q"]), (s["hr"], s["lo"], s["hi"], s["q"]))


C_ROWS = [(c, *arms(c)) for c in C_CONDS]; C_CTRL_ROWS = [(c, *arms(c)) for c in C_CTRL]

# ------------------------------------------------------------------ strings (the V13 wording)
XL_C = "Hazard ratio (95% CI), T90 above 10%\nversus 10% or less, within each apnea group"
LG_LEFT = f"No or mild apnea (AHI <15), {ST[LEFT]['n_low_o2']:,} of {ST[LEFT]['n']:,} in the low-oxygen arm"
LG_SEV = f"Severe apnea (AHI 30 or more), {ST[RIGHT]['n_low_o2']:,} of {ST[RIGHT]['n']:,}"
LG_OPEN = "Not significant"
BLOCK_HEAD = "Negative controls"; COL_HEAD = "HR (95% CI)"; QCOL_HEAD = "q"
# printed strings: half up on the file's values (the file stores 3 dp: an exact half at the printed precision is listed, not moved)
def is_half(v): s = repr(float(v)); return bool(re.fullmatch(r"-?\d+\.\d\d5", s))
COL_TXT, QCOL_TXT, QCOL_SIG, HALVES = [], [], [], []
for _c, _a, _s in C_ROWS + C_CTRL_ROWS:
    for _arm, (_hr, _lo, _hi, _q) in ((LEFT, _a), (RIGHT, _s)):
        COL_TXT.append(L.hr_ci(_hr, _lo, _hi)); QCOL_TXT.append(L.q_text(_q)); QCOL_SIG.append(bool(_q < QTHR))
        for nm, v in (("hr", _hr), ("lo", _lo), ("hi", _hi)):
            if is_half(v): HALVES.append(dict(condition=_c, arm=_arm, field=nm, stored=v, printed=L.r2(v)))
for _s in (XL_C, LG_LEFT, LG_SEV, LG_OPEN, BLOCK_HEAD, COL_HEAD, QCOL_HEAD, *COL_TXT, *QCOL_TXT, *C_CONDS, *C_CTRL):
    assert L.house_ok(_s), f"banned wording on the sheet: {_s!r}"

# ------------------------------------------------------------------ geometry, pinned to the V13 sheet (page points)
G13 = json.load(open(f"{WORK}/v13_geometry.json")); G = G13["panel_b"]
PX0, PY0, PX1, PY1 = G["panel_box"]                     # 484.859, 34.927, 969.585, 716.90 (the V13 box, unchanged since round 30)
PLOT_X0, PLOT_X1 = G["plot_x"]                          # 628.859, 833.385
PLOT_TOP, RULE_Y = G["plot_top"], G["axis_rule_y"]      # 96.848, 666.886 (= panel a's rule)
W_IN = (PX1 - PX0) / 72.0; H_IN = (PY1 - PY0) / 72.0
LEFT_PAD = 0.18; RIGHT_PAD = 0.18; COL_GAP = 0.12; QCOL_GAP = 0.10
X0 = (PLOT_X0 - PX0) / 72.0; AX_W = (PLOT_X1 - PLOT_X0) / 72.0
AX_Y0 = (PY1 - RULE_Y) / 72.0; AX_H = (RULE_Y - PLOT_TOP) / 72.0; TOP_STRIP = (PLOT_TOP - PY0) / 72.0
GAP_UNITS = 2.0; PAIR_OFF = 0.22; LEG_GAP = 0.07; LEG_ROW = 0.22; COLHEAD_Y = LEG_GAP
UNITS = (len(C_ROWS) - 1) + GAP_UNITS + (len(C_CTRL_ROWS) - 1) + 1.7
PITCH = AX_H / UNITS
N_ROWS_V13 = len([s for s in G["v13_strings"] if s["font"] == "ArialMT" and abs(s["size"] - 11.0) < 0.01 and s["origin"][0] < PLOT_X0 - 40 and s["text"] not in ("Negative controls",)])
UNITS_V13 = (N_ROWS_V13 - 2) + GAP_UNITS + 1.7 if N_ROWS_V13 else None      # for the record only (19 rows on V13: 14 conditions + 5 controls, the same stack formula)
C_TICKS_ALL = [0.5, 1.0, 2.0, 4.0, 8.0, 16.0]
DPT = 1.0        # round 49: every text element one point larger
COL_DPT = 0.0    # round 49 EXCEPTION: the two printed columns (HR (95% CI), q) and their headers stay 10 pt: at 11 pt the widest HR string
                 # (84.42 pt) exceeds the 77.31 pt its plot-box gap rule allows and the widest q string (33.93 pt) exceeds the 31.53 pt its
                 # column-gap rule allows at the pinned V13 column edges (moving the columns or the plot box is out of scope)
COL_FS, FS_ROW, FS_HEAD, FS_LEG, FS_XLAB, FS_TICK = 10.0 + COL_DPT, 11.0 + DPT, 11.0 + DPT, 10.0 + DPT, 11.0 + DPT, 10.0 + DPT
assert 2 * PAIR_OFF * PITCH >= COL_FS / 72.0, ("printed sub-rows would collide at this pitch", PITCH * 72)
assert (1 - 2 * PAIR_OFF) > 2 * PAIR_OFF * 1.2
INK9 = L.INK; GREY_LEFT = G["series_left_colour"]; ORANGE = L.ORANGE
COL_X_PAGE, QCOL_X_PAGE = G["col_right_x"], G["qcol_right_x"]   # 918.612, 956.621 (right edges of the two printed columns)
COL_X = (COL_X_PAGE - PX0) / 72.0; QCOL_X = (QCOL_X_PAGE - PX0) / 72.0
assert abs(QCOL_X - (W_IN - RIGHT_PAD)) < 0.02, (QCOL_X, W_IN - RIGHT_PAD)

RC = dict(L.RC); RC.update({"font.size": 11.0 + DPT, "axes.labelsize": 11.0 + DPT, "xtick.labelsize": 10.0 + DPT, "ytick.labelsize": 11.0 + DPT, "legend.fontsize": 10.0 + DPT,
                            "figure.facecolor": "none", "savefig.facecolor": "none", "axes.facecolor": "#ffffff"})


def fx(v): return v / W_IN
def fy(v): return v / H_IN
def r9(v): return round(float(v), 9)


def logx_strict(ax, ticks):
    ax.set_xscale("log"); ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}")); ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_minor_formatter(FuncFormatter(lambda v, _: "")); ax.tick_params(axis="x", which="minor", length=0)
    ax.tick_params(axis="x", labelsize=FS_TICK)


def stack(n_main, n_ctrl):
    y = [-float(i) for i in range(n_main)]; head = -(n_main - 1) - GAP_UNITS / 2.0
    y += [-(n_main - 1) - GAP_UNITS - float(i) for i in range(n_ctrl)]
    return np.array(y, float), head


def limits(values, pad=1.06):
    v = np.asarray(values, float); return float(v.min()) / pad, float(v.max()) * pad


def wrap_to(fig, text, fontsize, max_in):
    r = fig.canvas.get_renderer()
    def width(s):
        t = fig.text(0, 0, s, fontsize=fontsize); w = t.get_window_extent(renderer=r).width / fig.dpi; t.remove(); return w
    if "\n" in text:
        ws = [width(ln) for ln in text.split("\n")]; assert max(ws) <= max_in, ("authored line does not fit", text, ws, max_in); return text, ws
    lines, cur = [], ""
    for word in text.split():
        trial = f"{cur} {word}".strip()
        if not cur or width(trial) <= max_in: cur = trial
        else: lines.append(cur); cur = word
    lines.append(cur); return "\n".join(lines), [width(ln) for ln in lines]


def sheet_strings(fig, ax):
    lo, hi = ax.get_xlim(); out = [t for t in fig.texts if t.get_text().strip()] + [t for t in ax.texts if t.get_text().strip()]
    out += [t for t in ax.get_yticklabels() if t.get_text().strip()]
    for tick in ax.xaxis.get_major_ticks():
        if lo <= tick.get_loc() <= hi and tick.label1.get_text().strip(): out.append(tick.label1)
    for t in (ax.xaxis.label, ax.yaxis.label, ax.title):
        if t.get_text().strip(): out.append(t)
    for lg in getattr(fig, "legends", []): out += [t for t in lg.get_texts() if t.get_text().strip()]
    return out


with plt.rc_context(RC):
    fig = plt.figure(figsize=(W_IN, H_IN)); fig.patch.set_alpha(0.0); fig.canvas.draw(); _r = fig.canvas.get_renderer()
    def _w_in(s, bold=False):
        t = fig.text(0.5, 0.5, s, fontsize=COL_FS, fontweight="bold" if bold else "normal")
        w = t.get_window_extent(renderer=_r).width / fig.dpi; t.remove(); return w
    QCOL_W = max(max(_w_in(s, b) for s, b in zip(QCOL_TXT, QCOL_SIG)), _w_in(QCOL_HEAD)); COL_W = max(_w_in(s) for s in COL_TXT + [COL_HEAD])
    assert COL_X - COL_W - COL_GAP >= X0 + AX_W - 0.01, ("the widest HR string would run into the plot box", COL_W, COL_X - COL_GAP - (X0 + AX_W))
    assert QCOL_X - QCOL_W - QCOL_GAP >= COL_X - 0.01, ("the q column would run into the HR column", QCOL_W)
    ax = fig.add_axes([fx(X0), fy(AX_Y0), fx(AX_W), fy(AX_H)])
    rows = C_ROWS + C_CTRL_ROWS; y, y_head = stack(len(C_ROWS), len(C_CTRL_ROWS))
    # row labels wider than the gutter (X0 - LEFT_PAD, the V13 gutter) are set on two lines at the best space (the Fig 4d idiom,
    # "Cardiovascular / composite"); the rows the V13 sheet drew keep their single line
    GUTTER_IN = X0 - LEFT_PAD - 0.05
    def _w_lab(s_, size=FS_ROW - DPT):   # round 49: the line breaks stay as on V26 (the wrap rule evaluated at the V26 size 11 pt)
        t = fig.text(0.5, 0.5, s_, fontsize=size); w = t.get_window_extent(renderer=_r).width / fig.dpi; t.remove(); return w
    LABEL_DISPLAY = {}
    for cond_, _a, _s in rows:
        if _w_lab(cond_) <= GUTTER_IN: LABEL_DISPLAY[cond_] = cond_; continue
        parts = cond_.split(" "); best = None
        for k in range(1, len(parts)):
            l1, l2 = " ".join(parts[:k]), " ".join(parts[k:]); wmax = max(_w_lab(l1), _w_lab(l2))
            if wmax <= GUTTER_IN and (best is None or wmax < best[0]): best = (wmax, l1 + "\n" + l2)
        assert best is not None, ("row label does not fit the gutter even on two lines", cond_)
        LABEL_DISPLAY[cond_] = best[1]
    WRAPPED = {c: v for c, v in LABEL_DISPLAY.items() if "\n" in v}
    ax.axvline(1.0, color=INK9, lw=0.9, ls=(0, (4, 3)), zorder=1)
    drawn = []; expect_pts, expect_cis = {}, set(); DRAWN = []
    for (cond, a, s), yy in zip(rows, y):
        ya, ys = yy + PAIR_OFF, yy - PAIR_OFF          # no tie between the paired markers (V13 design, round 31)
        for (hr, lo, hi, q), col, mk, ms, ypt, arm in ((a, GREY_LEFT, "o", 36, ya, LEFT), (s, ORANGE, "D", 30, ys, RIGHT)):
            sig = q < QTHR
            ax.plot([lo, hi], [ypt] * 2, color=col, lw=1.7, solid_capstyle="round", zorder=3)
            ax.scatter([hr], [ypt], s=ms, marker=mk, zorder=4, facecolor=col if sig else "white", edgecolor="white" if sig else col, linewidths=0.8 if sig else 1.0)
            expect_cis.add((r9(lo), r9(hi), r9(ypt))); expect_pts[(r9(hr), r9(ypt))] = (to_rgba(col if sig else "white"), to_rgba("white" if sig else col))
            drawn += [hr, lo, hi]
            DRAWN.append(dict(condition=cond, arm=arm, hr=hr, lo=lo, hi=hi, q=q, sig=sig, y_data=ypt, control=cond in C_CTRL, key=f"fits/{cond}/{arm}", source=SRC))
    ax.set_yticks(y); ax.set_yticklabels([LABEL_DISPLAY[r[0]] for r in rows], fontsize=FS_ROW, color=INK9, linespacing=1.0)
    ax.yaxis.set_tick_params(length=0, pad=(X0 - LEFT_PAD) * 72.0)
    for t in ax.get_yticklabels(): t.set_ha("left")
    for sp in ("left", "top", "right"): ax.spines[sp].set_visible(False)
    C_LO, C_HI = limits(drawn); C_TICKS = [t for t in C_TICKS_ALL if C_LO <= t <= C_HI]
    logx_strict(ax, C_TICKS); ax.set_xlim(C_LO, C_HI); ax.set_ylim(float(y.min()) - 0.85, 0.85)
    ax.grid(False)                                       # V13 design: no pale vertical grid lines (round 31)
    XL_DRAWN, XL_W = wrap_to(fig, XL_C, FS_XLAB, AX_W + 0.50)   # round 49: allowance 0.30 -> 0.50 in (the authored second line is 238.8 pt at 12 pt, 226.1 allowed before; nothing sits beside it)
    assert XL_DRAWN.count("\n") == 1, ("the axis title must fit on two lines", XL_DRAWN, XL_W)
    ax.set_xlabel(XL_DRAWN, fontsize=FS_XLAB, labelpad=6)
    tr = blended_transform_factory(ax.transAxes, ax.transData)
    ax.text(-(X0 - LEFT_PAD) / AX_W, y_head, BLOCK_HEAD, transform=tr, fontsize=FS_HEAD, fontweight="bold", va="center", ha="left", color=INK9, clip_on=False, zorder=6)
    trCol = blended_transform_factory(fig.dpi_scale_trans, ax.transData)
    sub_y = [yy + off for yy in y for off in (PAIR_OFF, -PAIR_OFF)]
    for s_, qs, sig, yy in zip(COL_TXT, QCOL_TXT, QCOL_SIG, sub_y):
        fig.text(COL_X, yy, s_, transform=trCol, fontsize=COL_FS, color=INK9, ha="right", va="center")
        fig.text(QCOL_X, yy, qs, transform=trCol, fontsize=COL_FS, color=INK9, ha="right", va="center", fontweight="bold" if sig else "normal")
    for _x, _h in ((COL_X, COL_HEAD), (QCOL_X, QCOL_HEAD)):
        fig.text(_x, AX_Y0 + AX_H + COLHEAD_Y, _h, transform=fig.dpi_scale_trans, fontsize=COL_FS, color=INK9, ha="right", va="bottom")
    fig.legend(handles=[Line2D([], [], color=GREY_LEFT, marker="o", ls="-", lw=1.7, ms=6.0, markerfacecolor=GREY_LEFT, markeredgecolor="white", markeredgewidth=0.8, label=LG_LEFT),
                        Line2D([], [], color=ORANGE, marker="D", ls="-", lw=1.7, ms=5.5, markerfacecolor=ORANGE, markeredgecolor="white", markeredgewidth=0.8, label=LG_SEV)],
               loc="lower left", bbox_to_anchor=(fx(LEFT_PAD), fy(AX_Y0 + AX_H + LEG_GAP + LEG_ROW)), ncol=1, fontsize=FS_LEG, frameon=False, borderpad=0, borderaxespad=0.0, handletextpad=0.6, labelspacing=0.5, handlelength=1.6)
    fig.legend(handles=[(Line2D([], [], marker="o", ls="none", ms=6.0, markerfacecolor="white", markeredgecolor=GREY_LEFT, markeredgewidth=1.0),
                         Line2D([], [], marker="D", ls="none", ms=5.5, markerfacecolor="white", markeredgecolor=ORANGE, markeredgewidth=1.0))],
               labels=[LG_OPEN], handler_map={tuple: HandlerTuple(ndivide=None, pad=0.35)}, loc="lower left", bbox_to_anchor=(fx(LEFT_PAD), fy(AX_Y0 + AX_H + LEG_GAP)), ncol=1,
               fontsize=FS_LEG, frameon=False, borderpad=0, borderaxespad=0.0, handletextpad=0.6, handlelength=2.2)
    fig.canvas.draw(); rend = fig.canvas.get_renderer()
    # ---------------------------------------------------------------- read-back: every marker and interval is the fitted one, no tie, no grid
    got_cis, n_ref, n_other = set(), 0, 0
    for ln in ax.lines:
        xd, yd = np.asarray(ln.get_xdata(), float), np.asarray(ln.get_ydata(), float); assert len(xd) == 2
        if xd[0] == xd[1] == 1.0 and ln.is_dashed(): n_ref += 1
        elif yd[0] == yd[1]: got_cis.add((r9(xd[0]), r9(xd[1]), r9(yd[0])))
        else: n_other += 1
    assert n_ref == 1 and got_cis == expect_cis and n_other == 0, (n_ref, n_other, len(got_cis), len(expect_cis))
    assert not any(gl.get_visible() for gl in ax.get_xgridlines()) and not any(gl.get_visible() for gl in ax.get_ygridlines()), "grid lines drawn"
    got_pts = {}
    for c in ax.collections:
        off = np.asarray(c.get_offsets(), float); assert off.shape == (1, 2)
        got_pts[(r9(off[0, 0]), r9(off[0, 1]))] = (tuple(np.round(c.get_facecolor()[0], 4)), tuple(np.round(c.get_edgecolor()[0], 4)))
    assert set(got_pts) == set(expect_pts)
    for k, (efc, eec) in expect_pts.items(): assert got_pts[k] == (tuple(np.round(efc, 4)), tuple(np.round(eec, 4))), k
    assert len(got_pts) == 2 * len(rows) and len(got_cis) == 2 * len(rows)
    assert [t.get_text() for t in ax.get_yticklabels()] == [LABEL_DISPLAY[r[0]] for r in rows]
    EXPECT = sorted([LABEL_DISPLAY[r[0]] for r in rows] + [BLOCK_HEAD, XL_DRAWN, LG_LEFT, LG_SEV, LG_OPEN] + [f"{t:g}" for t in C_TICKS] + COL_TXT + QCOL_TXT + [COL_HEAD, QCOL_HEAD])
    GOT = sorted(t.get_text() for t in sheet_strings(fig, ax)); assert GOT == EXPECT, (GOT, EXPECT)
    for t in sheet_strings(fig, ax): assert t.get_fontsize() >= 9.0 and L.house_ok(t.get_text())
    for t in ax.get_yticklabels() + ax.texts:
        bb = t.get_window_extent(renderer=rend)
        assert bb.x0 / fig.dpi >= LEFT_PAD - 0.01 and bb.x1 / fig.dpi <= X0 - 0.02, ("row name outside the gutter", t.get_text(), bb.x1 / fig.dpi)   # round 49: margin 0.05 -> 0.02 in (Pulmonary hypertension ends 1.7 pt before the plot box at 12 pt)
    for t in sheet_strings(fig, ax):
        bb = t.get_window_extent(renderer=rend); assert bb.y0 / fig.dpi >= 0.02 and bb.y1 / fig.dpi <= H_IN - 0.02, ("text outside the panel box", t.get_text(), bb.y0 / fig.dpi, bb.y1 / fig.dpi)
    strings = []
    for t in sheet_strings(fig, ax):
        bb = t.get_window_extent(renderer=rend)
        strings.append(dict(text=t.get_text(), x0=round(PX0 + bb.x0 / fig.dpi * 72, 3), x1=round(PX0 + bb.x1 / fig.dpi * 72, 3),
                            y_top=round(PY1 - bb.y1 / fig.dpi * 72, 3), y_bot=round(PY1 - bb.y0 / fig.dpi * 72, 3), size=t.get_fontsize(), bold=t.get_fontweight() == "bold"))
    for d_, c in zip(DRAWN, ax.collections):
        off = c.get_offsets()[0]; px, py = ax.transData.transform(off)
        d_["x_page"] = round(PX0 + px / fig.dpi * 72, 3); d_["y_page"] = round(PY1 - py / fig.dpi * 72, 3)
        d_["hr_ci_text"] = L.hr_ci(d_["hr"], d_["lo"], d_["hi"]); d_["q_text"] = L.q_text(d_["q"])
    used = set()
    for c in ax.collections:
        for arr in (c.get_facecolor(), c.get_edgecolor()):
            for rowc in np.atleast_2d(arr): used.add(tuple(np.round(rowc[:3], 4)))
    for ln in ax.lines: used.add(tuple(np.round(matplotlib.colors.to_rgb(ln.get_color()), 4)))
    allowed = {tuple(np.round(matplotlib.colors.to_rgb(c), 4)) for c in (GREY_LEFT, ORANGE, INK9, "white")}
    assert not (used - allowed), used - allowed
    out_pdf = f"{WORK}/panel_b_raw.pdf"; fig.savefig(out_pdf, transparent=True); plt.close(fig)

# ------------------------------------------------------------------ the panel's own text layer -> DRAWN_RECORD rows (page coordinates)
met = L.sheet_font_metrics(json.load(open(f"{WORK}/base_text.json"))["spans"]); L.align_font_metrics(out_pdf, f"{WORK}/panel_b_m.pdf", met)
RB = L.spans_of(f"{WORK}/panel_b_m.pdf")
def page_of(s): return (s["origin"][0] + PX0, s["origin"][1] + PY0)
recs = []; taken = set()
def take(text, y_hint=None, x_hint=None, ha="left", **src):
    cands = [(i, s) for i, s in enumerate(RB) if i not in taken and s["text"].replace("\xa0", " ") == text]
    assert cands, ("drawn string not on the panel text layer", text)
    if y_hint is not None: cands.sort(key=lambda t: abs(page_of(t[1])[1] - y_hint) + (0 if x_hint is None else 0.01 * abs(page_of(t[1])[0] - x_hint)))
    i, s = cands[0]; taken.add(i); x, yb = page_of(s)
    recs.append(dict(text=text, x=round(x, 3), baseline=round(yb, 3), ha="left", size=s["size"], panel="b", **src)); return recs[-1]
for d_, ct, qt, sig in zip(DRAWN, COL_TXT, QCOL_TXT, QCOL_SIG):
    yp = d_["y_page"]
    take(ct, y_hint=yp, x_hint=COL_X_PAGE - 30, source_file=SRC, source_key=f"fits/{d_['condition']}/{d_['arm']}", source_value=[d_["hr"], d_["lo"], d_["hi"]], rule="arm_hr_ci")
    take(qt, y_hint=yp, x_hint=QCOL_X_PAGE - 10, source_file=SRC, source_key=f"fits/{d_['condition']}/{d_['arm']}", source_value=d_["q"], rule="arm_q")
for cond, a, s in rows:
    lst = "panel_c_controls_no_or_mild" if cond in C_CTRL else "panel_c_rows_no_or_mild"; idx = (C_CTRL if cond in C_CTRL else C_CONDS).index(cond)
    lines = LABEL_DISPLAY[cond].split("\n")
    for k, ln in enumerate(lines):
        take(ln, source_file=SRC, source_key=f"{lst}/{idx}", source_value=cond, rule="text" if len(lines) == 1 else f"wrap_line{k + 1}")
take(LG_LEFT, source_file=SRC, source_key="strata/no_or_mild", source_value=[ST[LEFT]["n_low_o2"], ST[LEFT]["n"]], rule="hdr_no_or_mild")
take(LG_SEV, source_file=SRC, source_key="strata/severe", source_value=[ST[RIGHT]["n_low_o2"], ST[RIGHT]["n"]], rule="hdr_severe")
for t in (LG_OPEN, BLOCK_HEAD, COL_HEAD, QCOL_HEAD): take(t, source_file="static:V13", source_key="V13 wording", source_value=t, rule="text")
for ln in XL_DRAWN.split("\n"): take(ln, source_file="static:V13", source_key="V13 axis title", source_value=ln, rule="text")
for t in C_TICKS: take(f"{t:g}", source_file="static:V13", source_key="axis tick (builder ticks inside the data range)", source_value=t, rule="text")
left = [s["text"] for i, s in enumerate(RB) if i not in taken and s["text"].strip()]; assert not left, ("panel strings not accounted for", left)   # whitespace-only spans (matplotlib legend spacing) carry no print

geom = dict(panel_box=[PX0, PY0, PX1, PY1], plot_box_page=[PLOT_X0, PLOT_TOP, PLOT_X1, RULE_Y], xlim=[C_LO, C_HI], ticks=C_TICKS, pitch_in=PITCH, pitch_pt=PITCH * 72,
            pitch_pt_v13=(AX_H / UNITS_V13 * 72) if UNITS_V13 else None, n_rows_v13=N_ROWS_V13, pair_off=PAIR_OFF, n_rows=len(rows), col_right_x=COL_X_PAGE, qcol_right_x=QCOL_X_PAGE,
            col_w_in=COL_W, qcol_w_in=QCOL_W, xlabel_lines=XL_DRAWN.split("\n"), xlabel_widths_in=XL_W, colours=dict(left=GREY_LEFT, severe=ORANGE, ink=INK9), order_note=ORDER_NOTE,
            design="V13 box kept (the round-30 rule: pitch = fill the box so the axis rule stays on panel a's rule), no grid lines, no ties")
json.dump(geom, open(f"{WORK}/panel_b_geometry.json", "w"), indent=1)
json.dump(dict(sheet=SHEET, panel="b", lane="V14_L4_APNEA", source=SRC, source_sha256=L.sha256(SRC), sidecar_step=SC["step"], v14_gate=SC["_v14_gate"],
               crosstab=dict(path=XT_PATH, sha256=L.sha256(XT_PATH), step=SCX["step"]["id"]), cohorts=dict(path=COH_PATH, n_analysis=N_COHORT), negcontrols_final=dict(path=NC_PATH, panel=NEG_CONTROLS),
               left_group=LEFT, right_group=RIGHT, strata=ST, rows_conditions=C_CONDS, rows_controls=C_CTRL, controls_absent=C_CTRL_ABSENT,
               header_left=LG_LEFT, header_severe=LG_SEV, key_open=LG_OPEN, xlabel=XL_C, xlabel_drawn=XL_DRAWN, ticks=C_TICKS, xlim=[C_LO, C_HI],
               low_support=LOW_SUPPORT, raw_p_only=WS["raw_p_only_no_or_mild"], dropped_unstable=WS["dropped_unstable_no_or_mild"], not_estimable=WS["not_estimable_no_or_mild"],
               values=DRAWN, strings=strings, drawn_records=recs, order_note=ORDER_NOTE, row_labels_wrapped=WRAPPED, stored_precision="3 dp in the numbers file, printed half up on the stored value",
               exact_half_strings=HALVES, geometry=geom), open(f"{VER}/{SHEET}_b_drawn.json", "w"), indent=1, ensure_ascii=False)
print(f"panel b: {len(C_ROWS)} conditions + {len(C_CTRL_ROWS)} controls (absent: {C_CTRL_ABSENT}), pitch {PITCH*72:.2f} pt (V13 {geom['pitch_pt_v13']:.2f} pt for {N_ROWS_V13} rows), xlim {C_LO:.3f} to {C_HI:.3f}, ticks {C_TICKS}")
print(f"  headers: {LG_LEFT!r} | {LG_SEV!r}; x label lines {XL_DRAWN.split(chr(10))} widths {[round(w,2) for w in XL_W]} in (plot {AX_W:.2f} in)")
print(f"  {ORDER_NOTE}; low support {LOW_SUPPORT}; exact halves at 3 dp: {len(HALVES)} {[(h['condition'], h['arm'], h['field'], h['stored']) for h in HALVES]}")
print(f"  open markers: {[(d['condition'], d['arm']) for d in DRAWN if not d['sig']]}")
print(f"  {len(recs)} drawn records; wrote {out_pdf}")
