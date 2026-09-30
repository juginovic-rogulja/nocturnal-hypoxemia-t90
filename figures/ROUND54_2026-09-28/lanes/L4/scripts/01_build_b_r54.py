#!/usr/bin/env python3
"""ROUND 54 (2026-09-29, lane L4): Main Fig 4 panel b at the round-54 sizes, re-laid on work/layout_r54.json. The round-49 builder's
data wiring is kept (fig3_within_stratum_v1.json gated by the v8.1 sidecar, the file's own row lists, crosstab and cohort cross-checks,
strings half up), the drawing is the page-pinned overlay idiom of panel a (page points, y down) instead of a matplotlib axes with a
legend object, so every baseline, column edge and handle is the layout spec's. Sizes: row labels, legend, tick labels, HR (95% CI) and q
columns 14 pt (the round-49 keep of the columns at 10 pt ends: the columns have their own edges in the spec), column headers and the two
axis-title lines 15 pt, the block header 14 pt bold. Marks as round 49: left group (no or mild) grey circles 6 pt, severe orange
diamonds 5.5 pt, filled with a 0.8 pt white edge when q < 0.05, open otherwise, 1.7 pt round-capped intervals, 0.9 pt dashed reference
rule at HR 1, no ties, no grid, the log axis on the data range of round 49 (limits(drawn, 1.06), asserted equal). Row labels on two
lines where the spec wraps them (Ventricular arrhythmia or cardiac arrest as on V30, plus the two the width forces, declared).
Outputs: work/panel_b_raw.pdf, work/panel_b_geometry.json, verify/Main_Fig4_b_drawn.json."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import json, math, os, re, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import l4lib as L, layout_r54 as LY
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND37_2026-09-14/figures/_common"); import v14lib as V

SHEET = "Main_Fig4"; SD = f"{L.LANE}/{SHEET}"; WORK = f"{SD}/work"; VER = f"{SD}/verify"; os.makedirs(WORK, exist_ok=True); os.makedirs(VER, exist_ok=True)
matplotlib, plt = L.mpl_setup()
from matplotlib.colors import to_rgb
S = LY.load(); B_ = S["b"]; FS = S["sizes"]

# ------------------------------------------------------------------ the numbers files, gated (v8.1), the round-49 wiring
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
    assert r[LEFT]["q"] < QTHR, (c, r[LEFT]["q"]); assert r[LEFT]["q_pool"] == r[RIGHT]["q_pool"] == "conditions", c; assert r[LEFT]["hr"] > 1, (c, r[LEFT]["hr"])
for c in C_CTRL:
    r = FITS[c]; assert r[LEFT]["q_pool"] == r[RIGHT]["q_pool"] == "negative controls", c; assert r[LEFT]["q"] >= QTHR and r[RIGHT]["q"] >= QTHR, (c, r[LEFT]["q"], r[RIGHT]["q"])
for c in C_CTRL_ABSENT:
    r = FITS[c]; assert r[LEFT]["unstable"] or r[RIGHT]["unstable"] or r[LEFT].get("q") is None, (c, "absent control is not unstable")
LOW_SUPPORT = [(c, a) for c in C_CONDS + C_CTRL for a in (LEFT, RIGHT) if FITS[c][a]["low_support"]]
assert sorted(LOW_SUPPORT) == sorted((x[0], x[1]) for x in WS["low_support_rows_no_or_mild"]), (LOW_SUPPORT, WS["low_support_rows_no_or_mild"])
hrs_left = [FITS[c][LEFT]["hr"] for c in C_CONDS]
ORDER_NOTE = "file order (panel_c_rows_no_or_mild), descending left-group HR" if hrs_left == sorted(hrs_left, reverse=True) else "file order (panel_c_rows_no_or_mild), NOT monotone in the left-group HR"
def arms(cond):
    a, s = FITS[cond][LEFT], FITS[cond][RIGHT]; return ((a["hr"], a["lo"], a["hi"], a["q"]), (s["hr"], s["lo"], s["hi"], s["q"]))
C_ROWS = [(c, *arms(c)) for c in C_CONDS]; C_CTRL_ROWS = [(c, *arms(c)) for c in C_CTRL]; rows = C_ROWS + C_CTRL_ROWS
assert [r["label"] for r in B_["rows"]] == [c for c, _, _ in rows], "the layout spec's row order is not the file's"

# ------------------------------------------------------------------ strings (the V13 wording)
XL_LINES = ["Hazard ratio (95% CI), T90 above 10%", "versus 10% or less, within each apnea group"]; XL_C = "\n".join(XL_LINES)
LG_LEFT = f"No or mild apnea (AHI <15), {ST[LEFT]['n_low_o2']:,} of {ST[LEFT]['n']:,} in the low-oxygen arm"
LG_SEV = f"Severe apnea (AHI 30 or more), {ST[RIGHT]['n_low_o2']:,} of {ST[RIGHT]['n']:,}"
LG_OPEN = "Not significant"; BLOCK_HEAD = "Negative controls"; COL_HEAD = "HR (95% CI)"; QCOL_HEAD = "q"
def is_half(v): return bool(re.fullmatch(r"-?\d+\.\d\d5", repr(float(v))))
COL_TXT, QCOL_TXT, QCOL_SIG, HALVES = [], [], [], []
for _c, _a, _s in rows:
    for _arm, (_hr, _lo, _hi, _q) in ((LEFT, _a), (RIGHT, _s)):
        COL_TXT.append(L.hr_ci(_hr, _lo, _hi)); QCOL_TXT.append(L.q_text(_q)); QCOL_SIG.append(bool(_q < QTHR))
        for nm, v in (("hr", _hr), ("lo", _lo), ("hi", _hi)):
            if is_half(v): HALVES.append(dict(condition=_c, arm=_arm, field=nm, stored=v, printed=L.r2(v)))
for _s in (XL_C, LG_LEFT, LG_SEV, LG_OPEN, BLOCK_HEAD, COL_HEAD, QCOL_HEAD, *COL_TXT, *QCOL_TXT, *C_CONDS, *C_CTRL): assert L.house_ok(_s), f"banned wording on the sheet: {_s!r}"

# ------------------------------------------------------------------ geometry from the layout spec (page points, y down)
PX0, PY0, PX1, PY1 = B_["box"]; W_PT, H_PT = PX1 - PX0, PY1 - PY0
LAB_X = B_["lab_x"]; COL_R = B_["col_r"]; QCOL_R = B_["qcol_r"]; AX_X0, AX_X1 = B_["plot_x0"], B_["plot_x1"]; RULE_Y = B_["rule_y"]; TICK_BASE = B_["tick_base"]
def limits(values, pad=1.06):
    v = np.asarray(values, float); return float(v.min()) / pad, float(v.max()) * pad
drawn_vals = [x for _c, a, s in rows for arm in (a, s) for x in arm[:3]]
C_LO, C_HI = limits(drawn_vals); assert abs(C_LO - S["xlim_b"][0]) < 1e-9 and abs(C_HI - S["xlim_b"][1]) < 1e-9, ("the data range differs from the layout spec's (round-49 values)", C_LO, C_HI, S["xlim_b"])
C_TICKS_ALL = [0.5, 1.0, 2.0, 4.0, 8.0, 16.0]; C_TICKS = [t for t in C_TICKS_ALL if C_LO <= t <= C_HI]; assert C_TICKS == list(S["ticks_b"])
def xv(v): return AX_X0 + (math.log(v) - math.log(C_LO)) / (math.log(C_HI) - math.log(C_LO)) * (AX_X1 - AX_X0)
REF_X = xv(1.0)
FS_ROW, FS_COL, FS_LEG, FS_TICK, FS_HEAD, FS_BLOCK, FS_XLAB = FS["row"], FS["col"], FS["key"], FS["tick"], FS["head"], FS["block"], FS["head"]
MARK_DY = S["mark_dy_em"] * FS_COL; LINE2 = 15.0
INK9 = L.INK; GREY_LEFT = json.load(open(f"{WORK}/v13_geometry.json"))["panel_b"]["series_left_colour"]; ORANGE = L.ORANGE
HANDLE = S["legend_handle_pt"]; HPAD = S["legend_pad_pt"]
MS = {LEFT: 6.0, RIGHT: 5.5}; MK = {LEFT: "o", RIGHT: "D"}; COLR = {LEFT: GREY_LEFT, RIGHT: ORANGE}   # ms 6.0 = scatter s 36, ms 5.5 = scatter s 30.25 (round 49: s 30)

RC = dict(L.RC); RC.update({"font.size": FS_ROW, "figure.facecolor": "none", "savefig.facecolor": "none"})
DRAWN, RECS = [], []
def T(x, y, s, size, ha="left", bold=False, **src):
    ov.text(x, y, s, fontsize=size, color=INK9, ha=ha, va="baseline", fontweight="bold" if bold else "normal")
    RECS.append(dict(text=s, x=round(x, 3), baseline=round(y, 3), ha=ha, size=size, panel="b", **src)); return ov.texts[-1]
with plt.rc_context(RC):
    fig = plt.figure(figsize=(W_PT / 72.0, H_PT / 72.0)); fig.patch.set_alpha(0.0)
    ov = fig.add_axes([0, 0, 1, 1]); ov.set_xlim(PX0, PX1); ov.set_ylim(PY1, PY0); ov.axis("off"); ov.patch.set_visible(False)
    rend = fig.canvas.get_renderer()
    def width_of(t): return t.get_window_extent(renderer=rend).width / fig.dpi * 72
    first_wake = B_["rows"][0]["wake"]
    ov.plot([REF_X, REF_X], [first_wake - MARK_DY - 12.0, RULE_Y], color=INK9, lw=0.9, ls=(0, (4, 3)), zorder=1, solid_capstyle="butt")
    ov.plot([AX_X0, AX_X1], [RULE_Y, RULE_Y], color=INK9, lw=0.8, zorder=4, solid_capstyle="projecting")
    for t in C_TICKS:
        ov.plot([xv(t), xv(t)], [RULE_Y, RULE_Y + LY.TICK_LEN], color=INK9, lw=0.8, zorder=4, solid_capstyle="butt")
        T(xv(t), TICK_BASE, f"{t:g}", FS_TICK, ha="center", source_file="static:V13", source_key="axis tick (builder ticks inside the data range)", source_value=t, rule="text")
    lst_of = lambda cond: ("panel_c_controls_no_or_mild", C_CTRL.index(cond)) if cond in C_CTRL else ("panel_c_rows_no_or_mild", C_CONDS.index(cond))
    for (cond, a, s), lr in zip(rows, B_["rows"]):
        assert lr["label"] == cond
        yl, ys = lr["wake"], lr["sleep"]; centre = (yl + ys) / 2 - MARK_DY; lines = lr["lines"]; lst, idx = lst_of(cond)
        if len(lines) == 1: T(LAB_X, centre + 0.36 * FS_ROW, cond, FS_ROW, source_file=SRC, source_key=f"{lst}/{idx}", source_value=cond, rule="text")
        else:
            assert len(lines) == 2 and " ".join(lines) == cond, (lines, cond)
            for k, ln in enumerate(lines): T(LAB_X, centre - 2.5 + k * LINE2, ln, FS_ROW, source_file=SRC, source_key=f"{lst}/{idx}", source_value=cond, rule=f"wrap_line{k + 1}", note="row label on two lines" + (" (as on V30)" if cond == "Ventricular arrhythmia or cardiac arrest" else " (round 54: a beside b does not fit the sheet width at 14 pt on one line)"))
        for (hr, lo, hi, q), arm, yb in ((a, LEFT, yl), (s, RIGHT, ys)):
            yc = yb - MARK_DY; sig = q < QTHR; col = COLR[arm]
            ov.plot([xv(lo), xv(hi)], [yc, yc], color=col, lw=1.7, solid_capstyle="round", zorder=3)
            ov.plot([xv(hr)], [yc], marker=MK[arm], ms=MS[arm], ls="None", color=col, markerfacecolor=col if sig else "white", markeredgecolor="white" if sig else col, markeredgewidth=0.8 if sig else 1.0, zorder=4)
            hs, qs = L.hr_ci(hr, lo, hi), L.q_text(q)
            T(COL_R, yb, hs, FS_COL, ha="right", source_file=SRC, source_key=f"fits/{cond}/{arm}", source_value=[hr, lo, hi], rule="arm_hr_ci")
            T(QCOL_R, yb, qs, FS_COL, ha="right", bold=sig, source_file=SRC, source_key=f"fits/{cond}/{arm}", source_value=q, rule="arm_q")
            DRAWN.append(dict(condition=cond, arm=arm, hr=hr, lo=lo, hi=hi, q=q, sig=sig, control=cond in C_CTRL, key=f"fits/{cond}/{arm}", source=SRC, x_page=round(xv(hr), 3), y_page=round(yc, 3), hr_ci_text=hs, q_text=qs))
    T(LAB_X, B_["block_base"], BLOCK_HEAD, FS_BLOCK, bold=True, source_file="static:V13", source_key="V13 wording", source_value=BLOCK_HEAD, rule="text")
    T(COL_R, B_["head_base"], COL_HEAD, FS_HEAD, ha="right", source_file="static:V13", source_key="V13 wording", source_value=COL_HEAD, rule="text")
    T(QCOL_R, B_["head_base"], QCOL_HEAD, FS_HEAD, ha="right", source_file="static:V13", source_key="V13 wording", source_value=QCOL_HEAD, rule="text")
    # legend: three lines, handles drawn by hand (line + filled marker, then the two open markers)
    lb = B_["legend_base"]; LX = B_["legend_x"]
    for k, (arm, lab, key, val) in enumerate(((LEFT, LG_LEFT, "strata/no_or_mild", [ST[LEFT]["n_low_o2"], ST[LEFT]["n"]]), (RIGHT, LG_SEV, "strata/severe", [ST[RIGHT]["n_low_o2"], ST[RIGHT]["n"]]))):
        mid = lb[k] - MARK_DY
        ov.plot([LX, LX + HANDLE], [mid, mid], color=COLR[arm], lw=1.7, solid_capstyle="round", zorder=3)
        ov.plot([LX + HANDLE / 2], [mid], marker=MK[arm], ms=MS[arm], ls="None", color=COLR[arm], markerfacecolor=COLR[arm], markeredgecolor="white", markeredgewidth=0.8, zorder=4)
        T(LX + HANDLE + HPAD, lb[k], lab, FS_LEG, source_file=SRC, source_key=key, source_value=val, rule="hdr_no_or_mild" if arm == LEFT else "hdr_severe")
    mid = lb[2] - MARK_DY
    ov.plot([LX + 3.0], [mid], marker="o", ms=MS[LEFT], ls="None", color=GREY_LEFT, markerfacecolor="white", markeredgecolor=GREY_LEFT, markeredgewidth=1.0, zorder=4)
    ov.plot([LX + 3.0 + 12.0], [mid], marker="D", ms=MS[RIGHT], ls="None", color=ORANGE, markerfacecolor="white", markeredgecolor=ORANGE, markeredgewidth=1.0, zorder=4)
    T(LX + HANDLE + HPAD + 6.0, lb[2], LG_OPEN, FS_LEG, source_file="static:V13", source_key="V13 wording", source_value=LG_OPEN, rule="text")
    XC = (AX_X0 + AX_X1) / 2
    for k, ln in enumerate(XL_LINES): T(XC, B_["title_base"][k], ln, FS_XLAB, ha="center", source_file="static:V13", source_key="V13 axis title", source_value=ln, rule="text")
    # read-back
    n_mk = sum(1 for ln in ov.lines if ln.get_marker() not in (None, "None", "")); n_ci = sum(1 for ln in ov.lines if ln.get_marker() in (None, "None", "") and ln.get_xdata()[0] != ln.get_xdata()[1] and ln.get_ydata()[0] == ln.get_ydata()[1] and ln.get_linewidth() == 1.7)
    assert n_mk == 2 * len(rows) + 4 and n_ci == 2 * len(rows) + 2, (n_mk, n_ci)
    used = {tuple(np.round(to_rgb(ln.get_color()), 4)) for ln in ov.lines} | {tuple(np.round(to_rgb(ln.get_markerfacecolor()), 4)) for ln in ov.lines if ln.get_marker() not in (None, "None", "")} | {tuple(np.round(to_rgb(ln.get_markeredgecolor()), 4)) for ln in ov.lines if ln.get_marker() not in (None, "None", "")}
    assert used <= {tuple(np.round(to_rgb(c), 4)) for c in (GREY_LEFT, ORANGE, INK9, "white")}, used
    for t in ov.texts: assert t.get_fontsize() >= 13.0 and L.house_ok(t.get_text())
    BB = {}; strings = []
    for t in ov.texts:
        bb = t.get_window_extent(renderer=rend); x0p, x1p = PX0 + bb.x0 / fig.dpi * 72, PX0 + bb.x1 / fig.dpi * 72; y0p, y1p = PY1 - bb.y1 / fig.dpi * 72, PY1 - bb.y0 / fig.dpi * 72
        assert x0p >= PX0 + 1.0 and x1p <= S["page_w"] - 1.0 and y0p >= PY0 + 1.0 and y1p <= PY1 - 1.0, ("text leaves panel b's page", t.get_text(), x0p, x1p, y0p, y1p)
        BB[(t.get_text(), round(t.get_position()[1], 3))] = (x0p, y0p, x1p, y1p); strings.append(dict(text=t.get_text(), x0=round(x0p, 3), x1=round(x1p, 3), y_top=round(y0p, 3), y_bot=round(y1p, 3), size=t.get_fontsize(), bold=t.get_fontweight() == "bold"))
    lab_right = max(v[2] for (s_, _), v in BB.items() if s_ in {ln for r in B_["rows"] for ln in r["lines"]} or s_ == BLOCK_HEAD)
    assert lab_right + LY.LAB_PLOT_GAP - 0.5 <= AX_X0, ("row labels run into the plot", lab_right, AX_X0)
    tip = max(xv(hi) for _c, a, s in rows for (_, _, hi, _) in (a, s)) + 0.85
    hr_left = min(v[0] for (s_, _), v in BB.items() if re.fullmatch(r"\d+\.\d\d \(\d+\.\d\d-\d+\.\d\d\)", s_))
    assert hr_left - tip >= 2.5, ("the HR column runs into the whiskers", hr_left, tip)
    hq = [min(v[0] for (s2, y2), v in BB.items() if y2 == rec["baseline"] and (s2.startswith("<0.") or re.fullmatch(r"0\.\d\d\d", s2))) - BB[(rec["text"], rec["baseline"])][2] for rec in RECS if rec["rule"] == "arm_hr_ci"]
    assert min(hq) >= 2.5, ("HR to q gap", min(hq))
    out_pdf = f"{WORK}/panel_b_raw.pdf"; fig.savefig(out_pdf, transparent=True); plt.close(fig)

# ------------------------------------------------------------------ the panel's own text layer -> every record at its origin
met = L.sheet_font_metrics(json.load(open(f"{WORK}/base_text.json"))["spans"]); L.align_font_metrics(out_pdf, f"{WORK}/panel_b_m.pdf", met)
RB = L.spans_of(f"{WORK}/panel_b_m.pdf"); taken = set()
for rec in RECS:
    cands = [(i, s) for i, s in enumerate(RB) if i not in taken and s["text"].replace("\xa0", " ") == rec["text"] and abs(s["origin"][1] + PY0 - rec["baseline"]) < 0.35]
    ha = rec["ha"]; best = None
    for i, s in cands:
        sx = PX0 + (s["origin"][0] if ha == "left" else (s["bbox"][2] if ha == "right" else (s["bbox"][0] + s["bbox"][2]) / 2))
        if abs(sx - rec["x"]) <= 0.6: best = i; break
    assert best is not None, ("drawn string not on the panel text layer at its origin", rec["text"], rec["x"], rec["baseline"], [(s["origin"], s["bbox"]) for _, s in cands][:2])
    taken.add(best)
left = [s["text"] for i, s in enumerate(RB) if i not in taken and s["text"].strip()]; assert not left, ("panel strings not accounted for", left)
WRAPPED = {r["label"]: "\n".join(r["lines"]) for r in B_["rows"] if len(r["lines"]) > 1}
geom = dict(panel_box=[PX0, PY0, PX1, PY1], plot_box_page=[AX_X0, first_wake - MARK_DY - 12.0, AX_X1, RULE_Y], xlim=[C_LO, C_HI], ticks=C_TICKS, sub_pitch=B_["sub_pitch"], pair_gap=B_["pair_gap"], n_rows=len(rows),
            col_right_x=COL_R, qcol_right_x=QCOL_R, rows=B_["rows"], head_base=B_["head_base"], block_base=B_["block_base"], legend_base=B_["legend_base"], title_base=B_["title_base"], tick_base=TICK_BASE, mark_dy=MARK_DY,
            xlabel_lines=XL_LINES, colours=dict(left=GREY_LEFT, severe=ORANGE, ink=INK9), order_note=ORDER_NOTE, widths=B_["widths"], wraps=B_["wraps"], design="round-54 re-lay on layout_r54.json, no grid lines, no ties")
json.dump(geom, open(f"{WORK}/panel_b_geometry.json", "w"), indent=1)
json.dump(dict(sheet=SHEET, panel="b", lane="L4 round 54", source=SRC, source_sha256=L.sha256(SRC), sidecar_step=SC["step"], v14_gate=SC["_v14_gate"],
               crosstab=dict(path=XT_PATH, sha256=L.sha256(XT_PATH), step=SCX["step"]["id"]), cohorts=dict(path=COH_PATH, n_analysis=N_COHORT), negcontrols_final=dict(path=NC_PATH, panel=NEG_CONTROLS),
               left_group=LEFT, right_group=RIGHT, strata=ST, rows_conditions=C_CONDS, rows_controls=C_CTRL, controls_absent=C_CTRL_ABSENT,
               header_left=LG_LEFT, header_severe=LG_SEV, key_open=LG_OPEN, xlabel=XL_C, xlabel_drawn=XL_C, ticks=C_TICKS, xlim=[C_LO, C_HI],
               low_support=LOW_SUPPORT, raw_p_only=WS["raw_p_only_no_or_mild"], dropped_unstable=WS["dropped_unstable_no_or_mild"], not_estimable=WS["not_estimable_no_or_mild"],
               values=DRAWN, strings=strings, drawn_records=RECS, order_note=ORDER_NOTE, row_labels_wrapped=WRAPPED, stored_precision="3 dp in the numbers file, printed half up on the stored value",
               exact_half_strings=HALVES, geometry=geom), open(f"{VER}/{SHEET}_b_drawn.json", "w"), indent=1, ensure_ascii=False)
print(f"panel b: {len(C_ROWS)} conditions + {len(C_CTRL_ROWS)} controls (absent: {C_CTRL_ABSENT}), sub-row pitch {B_['sub_pitch']} pt, xlim {C_LO:.3f} to {C_HI:.3f}, ticks {C_TICKS}; wraps {list(WRAPPED)}")
print(f"  {ORDER_NOTE}; low support {LOW_SUPPORT}; exact halves at 3 dp: {len(HALVES)}; open markers: {[(d['condition'], d['arm']) for d in DRAWN if not d['sig']]}; {len(RECS)} drawn records; wrote {out_pdf}")
