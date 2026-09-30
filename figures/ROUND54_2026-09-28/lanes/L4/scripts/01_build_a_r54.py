#!/usr/bin/env python3
"""ROUND 54 (2026-09-29, lane L4): Main Fig 4 panel a at the round-54 sizes, re-laid on work/layout_r54.json (layout_r54.py), the
round-49 builder's data wiring kept unchanged (the v8.2 step-121 files, the sidecar gates, the row rule, the strings half up, the drawn
record). Sizes: row labels, key, tick labels, HR (95% CI) and q columns 14 pt, column headers and the two caption lines 15 pt, the block
header 14 pt bold. Geometry: every baseline, column edge and the plot box come from the layout spec (sub-row pitch 14.5 pt, pair gap 2.5,
plot width from the tick-label rule, the three widest labels on two lines because a beside b does not fit the sheet width at 14 pt on one
line: declared in the record as wrap_line1 / wrap_line2). Marks as the generator: wake grey circles 6 pt, sleep blue diamonds 5.5 pt,
filled with a 0.8 pt white edge when q < 0.05, open otherwise, 1.7 pt round-capped intervals, 0.9 pt dashed reference rule at HR 1, no
ties, no grid, the log axis calibrated on the V13 range (the same xlim as round 49, so every mark keeps its value and its position on
the axis follows the new plot box). Outputs: work/panel_a_raw.pdf, work/panel_a_geometry.json, verify/Main_Fig4_a_drawn.json."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import json, math, os, re, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import l4lib as L, layout_r54 as LY
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND37_2026-09-14/figures/_common"); import v14lib as V

SHEET = "Main_Fig4"; SD = f"{L.LANE}/{SHEET}"; WORK = f"{SD}/work"; VER = f"{SD}/verify"; os.makedirs(WORK, exist_ok=True); os.makedirs(VER, exist_ok=True)
matplotlib, plt = L.mpl_setup()
from matplotlib.colors import to_rgb
V8_2_CUT = "2026-09-15 19:18:00"
S = LY.load(); A_ = S["a"]; FS = S["sizes"]

# ------------------------------------------------------------------ the numbers files, gated (v8.2: sidecar after 19:18 on 2026-09-15), the round-49 wiring
A_PATH = f"{L.NUM}/wake_sleep_healthy_A_full.csv"; J_PATH = f"{L.NUM}/wake_sleep_healthy.json"; NC_PATH = f"{L.NUM}/negcontrols_final.json"
GATES = {}
for p in (A_PATH, J_PATH):
    g = V.sidecar_v8(p); assert g["output_mtime"] >= V8_2_CUT, ("not the v8.2 re-extraction", p, g["output_mtime"]); GATES[p] = g
GATES[NC_PATH] = V.sidecar_v8(NC_PATH)
A = pd.read_csv(L.hydrated(A_PATH)); META = json.load(open(L.hydrated(J_PATH))); NCF = json.load(open(L.hydrated(NC_PATH)))
NEG_CONTROLS = list(NCF["panel_labels"]); assert len(NEG_CONTROLS) == NCF["n_controls"] == 5, NEG_CONTROLS
assert set(A.cond[A.neg]) == set(NEG_CONTROLS), ("the csv's controls are not the v8.1 panel", sorted(A.cond[A.neg]), NEG_CONTROLS)
assert not A.cond.duplicated().any()
N_A = int(META["strata"]["A_full"]["n"]); assert N_A == META["provenance"]["n_analysis_set"], (N_A, META["provenance"])
QTHR = float(META["provenance"]["q_threshold"]); assert QTHR == 0.05
DEATH = "Death from any cause"; assert DEATH in set(A.cond)
for c in ("wHR", "wlo", "whi", "wq", "sHR", "slo", "shi", "sq"): assert A[c].notna().all(), c

# ------------------------------------------------------------------ the row rule (the sheet's own since round 4)
d = A.copy(); d["gap"] = d.sHR - d.wHR
sig = d[(~d.neg) & (d.cond != DEATH) & ((d.wq < QTHR) | (d.sq < QTHR))]
TOP12 = sig.sort_values("sHR", ascending=False).head(12)
assert len(sig) >= 12, ("fewer than 12 significant conditions", len(sig))
main = pd.concat([TOP12, d[d.cond == DEATH]]).sort_values("gap", ascending=False)
ctrl = d[d.neg].sort_values("gap", ascending=False)
ROWS = list(main.itertuples()) + list(ctrl.itertuples())
assert len(main) == 13 and len(ctrl) == 5
DROPPED_SIG = sorted(set(sig.cond) - set(TOP12.cond), key=lambda c: -float(sig.set_index("cond").loc[c, "sHR"]))
N_SIG_W, N_SIG_S = int((A.wq < QTHR).sum()), int((A.sq < QTHR).sum())
assert [r["label"] for r in A_["rows"]] == [r.cond for r in ROWS], ("the layout spec's row order is not the row rule's", [r["label"] for r in A_["rows"]])

# ------------------------------------------------------------------ geometry from the layout spec (page points, y down)
PX0, PY0, PX1, PY1 = A_["box"]; W_PT, H_PT = PX1 - PX0, PY1 - PY0
LAB_X = A_["lab_x"]; COL_R = A_["col_r"]; QCOL_R = A_["qcol_r"]; AX_X0, AX_X1 = A_["plot_x0"], A_["plot_x1"]; RULE_Y = A_["rule_y"]
XLO, XHI = S["xlim_a"]; TICKS = list(S["ticks_a"]); TICK_BASE = A_["tick_base"]
def xv(v): return AX_X0 + (math.log(v) - math.log(XLO)) / (math.log(XHI) - math.log(XLO)) * (AX_X1 - AX_X0)
REF_X = xv(1.0)
FS_ROW, FS_COL, FS_KEY, FS_CAP, FS_TICK, FS_HEAD, FS_BLOCK = FS["row"], FS["col"], FS["key"], FS["head"], FS["tick"], FS["head"], FS["block"]
MARK_DY = S["mark_dy_em"] * FS_COL           # marker centre above the sub-row baseline
LINE2 = 15.0                                  # pitch of a two-line row label at 14 pt
GREY = "#8a9099"; BLUE = L.BLUE; INK = L.INK
MS = {"w": 6.0, "s": 5.5}; MK = {"w": "o", "s": "D"}; COL = {"w": GREY, "s": BLUE}
HANDLE = S["legend_handle_pt"]; HPAD = S["legend_pad_pt"]

# ------------------------------------------------------------------ strings
def hrci(r, p): return L.hr_ci(getattr(r, p + "HR"), getattr(r, p + "lo"), getattr(r, p + "hi"))
STR = {}; HALVES = []
def is_half(v): return bool(re.fullmatch(r"-?\d+\.\d\d5", repr(float(v))))
for r in ROWS:
    for p in ("w", "s"):
        STR[(r.cond, p)] = (hrci(r, p), L.q_text(float(getattr(r, p + "q"))), bool(float(getattr(r, p + "q")) < QTHR))
        for nm in ("HR", "lo", "hi"):
            if is_half(getattr(r, p + nm)): HALVES.append(dict(condition=r.cond, period=p, field=nm, stored=float(getattr(r, p + nm))))
vals = [float(getattr(r, p + k)) for r in ROWS for p in ("w", "s") for k in ("lo", "hi", "HR")]
assert XLO < min(vals) and max(vals) < XHI, ("a value leaves the axis range", XLO, min(vals), max(vals), XHI)
assert all(XLO < t < XHI for t in TICKS)
for s_ in ("Sleep-period T90", "Wake-period T90", "Not significant", "HR (95% CI)", "q", "Negative controls", "Hazard ratio (95% CI)", "per 1 SD higher T90", *[r.cond for r in ROWS]): assert L.house_ok(s_), s_

# ------------------------------------------------------------------ draw in page points (overlay axes, y downward)
RC = dict(L.RC); RC.update({"font.size": FS_ROW, "figure.facecolor": "none", "savefig.facecolor": "none"})
DRAWN, RECS = [], []
def T(x, y, s, size, ha="left", bold=False, **src):
    ov.text(x, y, s, fontsize=size, color=INK, ha=ha, va="baseline", fontweight="bold" if bold else "normal")
    RECS.append(dict(text=s, x=round(x, 3), baseline=round(y, 3), ha=ha, size=size, panel="a", **src))
with plt.rc_context(RC):
    fig = plt.figure(figsize=(W_PT / 72.0, H_PT / 72.0)); fig.patch.set_alpha(0.0)
    ov = fig.add_axes([0, 0, 1, 1]); ov.set_xlim(PX0, PX1); ov.set_ylim(PY1, PY0); ov.axis("off"); ov.patch.set_visible(False)
    # reference rule and axis
    first_wake = A_["rows"][0]["wake"]
    ov.plot([REF_X, REF_X], [first_wake - MARK_DY - 12.0, RULE_Y], color=INK, lw=0.9, ls=(0, (4, 3)), zorder=1, solid_capstyle="butt")
    ov.plot([AX_X0, AX_X1], [RULE_Y, RULE_Y], color=INK, lw=0.8, zorder=4, solid_capstyle="projecting")
    for t in TICKS:
        ov.plot([xv(t), xv(t)], [RULE_Y, RULE_Y + LY.TICK_LEN], color=INK, lw=0.8, zorder=4, solid_capstyle="butt")
        T(xv(t), TICK_BASE, f"{t:g}", FS_TICK, ha="center", source_file="static:V13", source_key="axis tick (ticks inside the V13 axis range)", source_value=t, rule="tick")
    # rows
    for r, lr in zip(ROWS, A_["rows"]):
        assert lr["label"] == r.cond
        yw, ys = lr["wake"], lr["sleep"]; centre = (yw + ys) / 2 - MARK_DY
        lines = lr["lines"]
        if len(lines) == 1: T(LAB_X, centre + 0.36 * FS_ROW, r.cond, FS_ROW, source_file=A_PATH, source_key=f"cond={r.cond}:cond", source_value=r.cond, rule="text")
        else:
            assert len(lines) == 2 and " ".join(lines) == r.cond, (lines, r.cond)
            T(LAB_X, centre - 2.5, lines[0], FS_ROW, source_file=A_PATH, source_key=f"cond={r.cond}:cond", source_value=r.cond, rule="wrap_line1", note="row label on two lines (round 54: a beside b does not fit the sheet width at 14 pt on one line)")
            T(LAB_X, centre - 2.5 + LINE2, lines[1], FS_ROW, source_file=A_PATH, source_key=f"cond={r.cond}:cond", source_value=r.cond, rule="wrap_line2", note="row label on two lines (round 54)")
        for p, yb in (("w", yw), ("s", ys)):
            yc = yb - MARK_DY; hr, lo, hi, q = (float(getattr(r, p + k)) for k in ("HR", "lo", "hi", "q")); sigp = q < QTHR; col = COL[p]
            ov.plot([xv(lo), xv(hi)], [yc, yc], color=col, lw=1.7, solid_capstyle="round", zorder=3)
            ov.plot([xv(hr)], [yc], marker=MK[p], ms=MS[p], ls="None", color=col, markerfacecolor=col if sigp else "white", markeredgecolor="white" if sigp else col, markeredgewidth=0.8 if sigp else 1.0, zorder=4)
            hs, qs, _ = STR[(r.cond, p)]
            T(COL_R, yb, hs, FS_COL, ha="right", source_file=A_PATH, source_key=f"cond={r.cond}:{p}HR,{p}lo,{p}hi", source_value=[hr, lo, hi], rule="hr_ci_2dp")
            T(QCOL_R, yb, qs, FS_COL, ha="right", bold=sigp, source_file=A_PATH, source_key=f"cond={r.cond}:{p}q", source_value=q, rule="q3")
            DRAWN.append(dict(condition=r.cond, period="wake" if p == "w" else "sleep", hr=hr, lo=lo, hi=hi, q=q, sig=sigp, control=bool(r.neg), x_page=round(xv(hr), 3), y_page=round(yc, 3), hr_ci_text=hs, q_text=qs, gap=float(r.gap)))
    T(LAB_X, A_["block_base"], "Negative controls", FS_BLOCK, bold=True, source_file="static:V13", source_key="V13 block header", source_value="Negative controls", rule="text")
    T(COL_R, A_["head_base"], "HR (95% CI)", FS_HEAD, ha="right", source_file="static:V13", source_key="V13 column header", source_value="HR (95% CI)", rule="text")
    T(QCOL_R, A_["head_base"], "q", FS_HEAD, ha="right", source_file="static:V13", source_key="V13 column header", source_value="q", rule="text")
    # key: line 1 = sleep handle + label, wake handle + label; line 2 = the two open markers + 'Not significant'
    kb1, kb2 = A_["key_base"]; kmid1 = kb1 - MARK_DY; kmid2 = kb2 - MARK_DY; x = LAB_X; KEY = {}
    for p, lab in (("s", "Sleep-period T90"), ("w", "Wake-period T90")):
        ov.plot([x, x + HANDLE], [kmid1, kmid1], color=COL[p], lw=1.7, solid_capstyle="round", zorder=3)
        ov.plot([x + HANDLE / 2], [kmid1], marker=MK[p], ms=MS[p], ls="None", color=COL[p], markerfacecolor=COL[p], markeredgecolor="white", markeredgewidth=0.8, zorder=4)
        T(x + HANDLE + HPAD, kb1, lab, FS_KEY, source_file="static:V13", source_key="V13 key", source_value=lab, rule="text"); KEY[p] = (x, x + HANDLE, x + HANDLE + HPAD)
        tw = ov.texts[-1].get_window_extent(renderer=fig.canvas.get_renderer()).width / fig.dpi * 72; x = x + HANDLE + HPAD + tw + 14.0
    ov.plot([LAB_X + 5.5], [kmid2], marker="D", ms=MS["s"], ls="None", color=BLUE, markerfacecolor="white", markeredgecolor=BLUE, markeredgewidth=1.0, zorder=4)
    ov.plot([LAB_X + 5.5 + 13.6], [kmid2], marker="o", ms=MS["w"], ls="None", color=GREY, markerfacecolor="white", markeredgecolor=GREY, markeredgewidth=1.0, zorder=4)
    T(LAB_X + HANDLE + HPAD + 9.0, kb2, "Not significant", FS_KEY, source_file="static:V13", source_key="V13 key (open markers)", source_value="Not significant", rule="text")
    # caption, two lines centred on the axis
    XC = (AX_X0 + AX_X1) / 2; cb1, cb2 = A_["cap_base"]
    T(XC, cb1, "Hazard ratio (95% CI)", FS_CAP, ha="center", source_file="static:V13", source_key="V13 axis caption", source_value="Hazard ratio (95% CI)", rule="text")
    T(XC, cb2, "per 1 SD higher T90", FS_CAP, ha="center", source_file="static:V13", source_key="V13 axis caption", source_value="per 1 SD higher T90", rule="text")
    # read-back of the matplotlib artists before saving
    n_mk = sum(1 for ln in ov.lines if ln.get_marker() not in (None, "None", "")); n_ci = sum(1 for ln in ov.lines if ln.get_marker() in (None, "None", "") and ln.get_xdata()[0] != ln.get_xdata()[1] and ln.get_ydata()[0] == ln.get_ydata()[1] and ln.get_linewidth() == 1.7)
    assert n_mk == 2 * len(ROWS) + 4 and n_ci == 2 * len(ROWS) + 2, (n_mk, n_ci)
    used = {tuple(np.round(to_rgb(ln.get_color()), 4)) for ln in ov.lines} | {tuple(np.round(to_rgb(ln.get_markerfacecolor()), 4)) for ln in ov.lines if ln.get_marker() not in (None, "None", "")} | {tuple(np.round(to_rgb(ln.get_markeredgecolor()), 4)) for ln in ov.lines if ln.get_marker() not in (None, "None", "")}
    assert used <= {tuple(np.round(to_rgb(c), 4)) for c in (GREY, BLUE, INK, "white")}, used
    for t in ov.texts: assert t.get_fontsize() >= 13.0 and L.house_ok(t.get_text())
    rend = fig.canvas.get_renderer(); BB = {}
    for t in ov.texts:
        bb = t.get_window_extent(renderer=rend); x0p, x1p = PX0 + bb.x0 / fig.dpi * 72, PX0 + bb.x1 / fig.dpi * 72; y0p, y1p = PY1 - bb.y1 / fig.dpi * 72, PY1 - bb.y0 / fig.dpi * 72
        assert x0p >= PX0 + 1.0 and x1p <= PX1 - 1.0 and y0p >= PY0 + 1.0 and y1p <= PY1 - 1.0, ("text leaves panel a's page", t.get_text(), x0p, x1p, y0p, y1p)
        BB[(t.get_text(), round(t.get_position()[1], 3))] = (x0p, y0p, x1p, y1p)
    # placement rules: labels end before the plot, the HR column starts after the widest whisker cap, HR to q at least 2.5 pt
    lab_right = max(v[2] for (s, _), v in BB.items() if s in {ln for r in A_["rows"] for ln in r["lines"]} or s == "Negative controls")
    assert lab_right + LY.LAB_PLOT_GAP - 0.5 <= AX_X0, ("row labels run into the plot", lab_right, AX_X0)
    tip = max(xv(float(getattr(r, p + "hi"))) for r in ROWS for p in ("w", "s")) + 0.85
    hr_left = min(v[0] for (s, _), v in BB.items() if re.fullmatch(r"\d+\.\d\d \(\d+\.\d\d-\d+\.\d\d\)", s))
    assert hr_left - tip >= 2.5, ("the HR column runs into the whiskers", hr_left, tip)
    hq = []
    for rec in RECS:
        if rec["rule"] == "hr_ci_2dp":
            qrec = [v for (s, y), v in BB.items() if y == rec["baseline"] and (s.startswith("<0.") or re.fullmatch(r"0\.\d\d\d", s))]; hrb = BB[(rec["text"], rec["baseline"])]
            hq.append(min(q[0] for q in qrec) - hrb[2])
    assert min(hq) >= 2.5, ("HR to q gap", min(hq))
    out_pdf = f"{WORK}/panel_a_raw.pdf"; fig.savefig(out_pdf, transparent=True); plt.close(fig)

# ------------------------------------------------------------------ the panel's own text layer: every record must be found at its origin (page frame = panel frame + box origin)
met = L.sheet_font_metrics(json.load(open(f"{WORK}/base_text.json"))["spans"]); L.align_font_metrics(out_pdf, f"{WORK}/panel_a_m.pdf", met)
RB = L.spans_of(f"{WORK}/panel_a_m.pdf"); taken = set()
for rec in RECS:
    cands = [(i, s) for i, s in enumerate(RB) if i not in taken and s["text"].replace("\xa0", " ") == rec["text"] and abs(s["origin"][1] + PY0 - rec["baseline"]) < 0.35]
    ha = rec["ha"]; best = None
    for i, s in cands:
        sx = PX0 + (s["origin"][0] if ha == "left" else (s["bbox"][2] if ha == "right" else (s["bbox"][0] + s["bbox"][2]) / 2))
        if abs(sx - rec["x"]) <= 0.6: best = i; break
    assert best is not None, ("drawn string not on the panel text layer at its origin", rec["text"], rec["x"], rec["baseline"], [(s["origin"], s["bbox"]) for _, s in cands][:2])
    taken.add(best)
left = [s["text"] for i, s in enumerate(RB) if i not in taken and s["text"].strip()]; assert not left, ("panel strings not accounted for", left)
geom = dict(panel_box=[PX0, PY0, PX1, PY1], axis_x=[AX_X0, AX_X1], rule_y=RULE_Y, ref_x=REF_X, xlim=[XLO, XHI], ticks=TICKS, label_x=LAB_X, rows=A_["rows"], sub_pitch=A_["sub_pitch"], pair_gap=A_["pair_gap"],
            head_base=A_["head_base"], block_base=A_["block_base"], col_right=COL_R, qcol_right=QCOL_R, key_base=A_["key_base"], cap_base=A_["cap_base"], tick_base=TICK_BASE, mark_dy=MARK_DY,
            colours=dict(wake=GREY, sleep=BLUE, ink=INK), marker_pt=MS, wraps=A_["wraps"], widths=A_["widths"], design="round-54 re-lay on layout_r54.json; markers and lines as figure3A_polish.py; no ties, no grid")
json.dump(geom, open(f"{WORK}/panel_a_geometry.json", "w"), indent=1)
json.dump(dict(sheet=SHEET, panel="a", lane="L4 round 54", round="v8.2 numbers, round-54 sizes", sources=dict(A_full=dict(path=A_PATH, sha256=GATES[A_PATH]["sha256"], step=GATES[A_PATH]["step"], output_mtime=GATES[A_PATH]["output_mtime"]),
               meta=dict(path=J_PATH, sha256=GATES[J_PATH]["sha256"], step=GATES[J_PATH]["step"], output_mtime=GATES[J_PATH]["output_mtime"]), negcontrols_final=dict(path=NC_PATH, panel=NEG_CONTROLS)),
               n_analysis_set=N_A, q_threshold=QTHR, n_sig_wake=N_SIG_W, n_sig_sleep=N_SIG_S, row_rule="significant for either period (q < 0.05, death excluded from the rule), the 12 largest sleep-period estimates, plus death, plus the controls; ordered by the asleep-minus-awake gap",
               rows_main=[r.cond for r in ROWS if not r.neg], rows_controls=[r.cond for r in ROWS if r.neg], significant_not_drawn=DROPPED_SIG, values=DRAWN, drawn_records=RECS, row_labels_wrapped=A_["wraps"],
               exact_half_strings=HALVES, stored_precision="full precision in the csv (printed half up 2 dp, q 3 dp)", geometry=geom), open(f"{VER}/{SHEET}_a_drawn.json", "w"), indent=1, ensure_ascii=False)
print(f"panel a: rows {[r.cond for r in ROWS]}; not drawn though significant: {DROPPED_SIG}; n {N_A:,}; sig wake {N_SIG_W} sleep {N_SIG_S}; ticks {TICKS}; xlim {XLO:.3f} to {XHI:.3f}; halves {len(HALVES)}; {len(RECS)} records; wraps {list(A_['wraps'])}; wrote {out_pdf}")
