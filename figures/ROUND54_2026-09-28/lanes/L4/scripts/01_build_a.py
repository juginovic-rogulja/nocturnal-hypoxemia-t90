#!/usr/bin/env python3
"""ROUND 49 (2026-09-26, lane L4): copy of the round-38 builder with every text one point larger (DPT = 1.0 added inside T(); sizes
11 -> 12 for row labels and captions, 10 -> 11 for the printed columns, key and ticks); positions, marks, colours and strings unchanged.
V14 lane L4, round v8.2 (2026-09-15 night), Main Fig 4 panel a: wake-period against sleep-period T90 per 1 SD in the full analysis
set, rebuilt from the RE-EXTRACTED step-121 files (the owner's v8.2 rule: no panel keeps old data). Until now the panel was the V13
render (group P). Data wiring = FINAL_FIGURES_2026-08-14/_scripts/figure3A_polish.py (the sheet's generator) with the trimmed row
rule the sheet has carried since round 4 (ROUND30 04_readback_a.py): the conditions significant for either period (FDR q < 0.05,
death excluded from the rule), the 12 with the largest sleep-period estimate, plus death, plus the negative controls the file
carries, ordered by the asleep-minus-awake gap (descending), controls the same. Design = the V13 panel a, pinned in page points:
row-label x and baselines, the two sub-rows of each pair (wake above, sleep below), the HR (95% CI) and q columns and their right
edges, the key rows, the log axis calibration (x = 240.604 + 147.35 ln HR, the V13 tick spacing and reference rule), the axis rule
at y 666.886 from x 151.68 to 356.16, tick labels at baseline 680.369, the two caption lines. Markers and lines as the generator:
wake grey circles 6 pt, sleep blue diamonds 5.5 pt, filled with a 0.8 pt white edge when q < 0.05, open with the series-colour ring
otherwise, 1.7 pt round-capped intervals, 0.9 pt dashed reference rule at HR 1, no ties, no grid (V13). Fonts Arial (Type 42).
Outputs: work/panel_a_raw.pdf, work/panel_a_geometry.json, verify/Main_Fig4_a_drawn.json (DRAWN_RECORD rows).
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import json, math, os, re, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import l4lib as L
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND37_2026-09-14/figures/_common"); import v14lib as V

SHEET = "Main_Fig4"; SD = f"{L.LANE}/{SHEET}"; WORK = f"{SD}/work"; VER = f"{SD}/verify"
matplotlib, plt = L.mpl_setup()
from matplotlib.colors import to_rgb
V8_2_CUT = "2026-09-15 19:18:00"
DPT = 1.0   # round 49: every text element one point larger, applied inside the text helper T() (the placement read-back uses the drawn artists)

# ------------------------------------------------------------------ the numbers files, gated (v8.2: sidecar after 19:18 today)
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

# ------------------------------------------------------------------ the row rule (the sheet's own since round 4, read back in round 30)
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

# ------------------------------------------------------------------ V13 pins (page points, from the V13 text layer and the round-30 walk)
IN = json.load(open(f"{WORK}/base_text.json")); SP = L.dedupe(IN["spans"])
pa = [s for s in SP if s["origin"][0] < 490 and 20 < s["origin"][1] < 717]
labs13 = sorted([s for s in pa if s["size"] == 11.0 and s["bbox"][0] < 100 and s["text"] != "Negative controls"], key=lambda s: s["origin"][1])
hr13 = [s for s in pa if re.fullmatch(r"\d+\.\d\d \(\d+\.\d\d-\d+\.\d\d\)", s["text"])]; q13 = [s for s in pa if re.fullmatch(r"<0\.001|0\.\d{3}", s["text"])]
assert len(labs13) == 18 and len(hr13) == 36 and len(q13) == 36, (len(labs13), len(hr13), len(q13))
LAB_X = float(np.median([s["origin"][0] for s in labs13])); FS_ROW = 11.0
Y_LAB = [s["origin"][1] for s in labs13]; PITCH = float(np.median(np.diff(Y_LAB)[:12])); Y_FIRST = Y_LAB[0]
HEAD13 = next(s for s in pa if s["text"] == "Negative controls"); HEAD_Y = HEAD13["origin"][1]
first_hr = sorted([s for s in hr13 if abs(s["origin"][1] - Y_FIRST) < 20], key=lambda s: s["origin"][1])
OFF_W, OFF_S = first_hr[0]["origin"][1] - Y_FIRST, first_hr[1]["origin"][1] - Y_FIRST      # -7.95 / +5.10: wake sub-row above, sleep below
COL_R = float(np.median([s["bbox"][2] for s in hr13])); QCOL_R = float(np.median([s["bbox"][2] for s in q13]))
HR_HEAD = next(s for s in pa if s["text"] == "HR (95% CI)"); Q_HEAD = next(s for s in pa if s["text"] == "q" and s["origin"][1] < 100)
KEY_S = next(s for s in pa if s["text"] == "Sleep-period T90"); KEY_W = next(s for s in pa if s["text"] == "Wake-period T90"); KEY_O = next(s for s in pa if s["text"] == "Not significant")
CAP1 = next(s for s in pa if s["text"] == "Hazard ratio (95% CI)"); CAP2 = next(s for s in pa if s["text"] == "per 1 SD higher T90")
ticks13 = {s["text"]: s for s in pa if s["text"] in ("0.6", "0.8", "1", "1.5", "2") and s["origin"][1] > 670}; TICK_BASE = float(np.median([s["origin"][1] for s in ticks13.values()]))
G13 = json.load(open(f"{WORK}/v13_geometry.json")); PB = G13["panel_b"]
RULE_Y = float(PB["rule_a_y"]); REF_X = float(PB["ref_rule_x_panel_a"])                 # 666.886, 240.604 (the round-30 walk of the V13 panel a)
tc = {float(k): (v["bbox"][0] + v["bbox"][2]) / 2 for k, v in ticks13.items()}
B_LN = (tc[2.0] - tc[0.6]) / (math.log(2.0) - math.log(0.6))                              # 147.35 pt per ln unit, the V13 tick spacing
def xv(v): return REF_X + B_LN * math.log(v)
assert max(abs(xv(v) - x) for v, x in tc.items()) < 0.6, {v: round(xv(v) - x, 2) for v, x in tc.items()}   # label centres sit within a glyph side bearing of the rule
AX_X0, AX_X1 = float(G13.get("panel_a_axis_x", [151.68, 356.16])[0]), float(G13.get("panel_a_axis_x", [151.68, 356.16])[1])   # measured on the V13 render (150 dpi, 0.5 pt)
XLO, XHI = math.exp((AX_X0 - REF_X) / B_LN), math.exp((AX_X1 - REF_X) / B_LN)
PANEL = [0.0, 0.0, 490.0, 716.9]; W_PT, H_PT = PANEL[2] - PANEL[0], PANEL[3] - PANEL[1]
# key handles (measured on the V13 render): sleep handle line x 14.72 to 32.0, wake handle ends 5.25 pt before its label, open row handles
HANDLE_LEN = 17.3; KEY_GAP = 5.25; KEY_MID = KEY_S["origin"][1] - 3.8; OPEN_MID = KEY_O["origin"][1] - 3.0
S_HANDLE = (KEY_S["origin"][0] - KEY_GAP - HANDLE_LEN, KEY_S["origin"][0] - KEY_GAP); W_HANDLE = (KEY_W["origin"][0] - KEY_GAP - HANDLE_LEN, KEY_W["origin"][0] - KEY_GAP)
OPEN_D_X, OPEN_O_X = 20.0, 33.6
GREY = "#8a9099"; BLUE = L.BLUE; INK = L.INK
MS = {"w": 6.0, "s": 5.5}; MK = {"w": "o", "s": "D"}; COL = {"w": GREY, "s": BLUE}
FS_COL, FS_KEY, FS_CAP, FS_TICK = 10.0, 10.0, 11.0, 10.0
TICKS_ALL = [0.5, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0, 4.0]

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
assert XLO < min(vals) and max(vals) < XHI, ("a value leaves the V13 axis", XLO, min(vals), max(vals), XHI)
TICKS = [t for t in TICKS_ALL if XLO < t < XHI]
for s_ in ("Sleep-period T90", "Wake-period T90", "Not significant", "HR (95% CI)", "q", "Negative controls", "Hazard ratio (95% CI)", "per 1 SD higher T90", *[r.cond for r in ROWS]): assert L.house_ok(s_), s_

# ------------------------------------------------------------------ draw in page points (overlay axes, y downward)
RC = dict(L.RC); RC.update({"font.size": 10.0 + DPT, "figure.facecolor": "none", "savefig.facecolor": "none"})
DRAWN, RECS = [], []
def T(x, y, s, size, ha="left", bold=False, **src):
    ov.text(x, y, s, fontsize=size + DPT, color=INK, ha=ha, va="baseline", fontweight="bold" if bold else "normal")   # round 49: +DPT
    RECS.append(dict(text=s, x=round(x, 3), baseline=round(y, 3), ha=ha, size=size + DPT, panel="a", **src))
with plt.rc_context(RC):
    fig = plt.figure(figsize=(W_PT / 72.0, H_PT / 72.0)); fig.patch.set_alpha(0.0)
    ov = fig.add_axes([0, 0, 1, 1]); ov.set_xlim(PANEL[0], PANEL[2]); ov.set_ylim(PANEL[3], PANEL[1]); ov.axis("off"); ov.patch.set_visible(False)
    # reference rule and axis
    ov.plot([REF_X, REF_X], [Y_FIRST + OFF_W - 12.0, RULE_Y], color=INK, lw=0.9, ls=(0, (4, 3)), zorder=1, solid_capstyle="butt")
    ov.plot([AX_X0, AX_X1], [RULE_Y, RULE_Y], color=INK, lw=0.8, zorder=4, solid_capstyle="projecting")
    for t in TICKS:
        ov.plot([xv(t), xv(t)], [RULE_Y, RULE_Y + 3.0], color=INK, lw=0.8, zorder=4, solid_capstyle="butt")
        T(xv(t), TICK_BASE, f"{t:g}", FS_TICK, ha="center", source_file="static:V13", source_key="axis tick (ticks inside the V13 axis range)", source_value=t, rule="tick")
    # rows
    y_lab = Y_FIRST; ROW_Y = {}
    for i, r in enumerate(ROWS):
        if i == 13: y_lab += PITCH   # the block gap before the controls (one extra pitch, V13: 59.26 pt between the last condition and the first control)
        ROW_Y[r.cond] = y_lab
        T(LAB_X, y_lab, r.cond, FS_ROW, source_file=A_PATH, source_key=f"cond={r.cond}:cond", source_value=r.cond, rule="text")
        for p, off in (("w", OFF_W), ("s", OFF_S)):
            yb = y_lab + off; yc = yb - 3.45; hr, lo, hi, q = (float(getattr(r, p + k)) for k in ("HR", "lo", "hi", "q")); sigp = q < QTHR; col = COL[p]
            ov.plot([xv(lo), xv(hi)], [yc, yc], color=col, lw=1.7, solid_capstyle="round", zorder=3)
            ov.plot([xv(hr)], [yc], marker=MK[p], ms=MS[p], ls="None", color=col, markerfacecolor=col if sigp else "white", markeredgecolor="white" if sigp else col, markeredgewidth=0.8 if sigp else 1.0, zorder=4)
            hs, qs, _ = STR[(r.cond, p)]
            T(COL_R, yb, hs, FS_COL, ha="right", source_file=A_PATH, source_key=f"cond={r.cond}:{p}HR,{p}lo,{p}hi", source_value=[hr, lo, hi], rule="hr_ci_2dp")
            T(QCOL_R, yb, qs, FS_COL, ha="right", bold=sigp, source_file=A_PATH, source_key=f"cond={r.cond}:{p}q", source_value=q, rule="q3")
            DRAWN.append(dict(condition=r.cond, period="wake" if p == "w" else "sleep", hr=hr, lo=lo, hi=hi, q=q, sig=sigp, control=bool(r.neg), x_page=round(xv(hr), 3), y_page=round(yc, 3), hr_ci_text=hs, q_text=qs, gap=float(r.gap)))
        y_lab += PITCH
    T(LAB_X, HEAD_Y, "Negative controls", FS_ROW, bold=True, source_file="static:V13", source_key="V13 block header", source_value="Negative controls", rule="text")
    T(HR_HEAD["bbox"][2], HR_HEAD["origin"][1], "HR (95% CI)", FS_COL, ha="right", source_file="static:V13", source_key="V13 column header", source_value="HR (95% CI)", rule="text")
    T(Q_HEAD["bbox"][2], HR_HEAD["origin"][1], "q", FS_COL, ha="right", source_file="static:V13", source_key="V13 column header", source_value="q", rule="text")
    # key
    for (x0, x1), p in ((S_HANDLE, "s"), (W_HANDLE, "w")):
        ov.plot([x0, x1], [KEY_MID, KEY_MID], color=COL[p], lw=1.7, solid_capstyle="round", zorder=3)
        ov.plot([(x0 + x1) / 2], [KEY_MID], marker=MK[p], ms=MS[p], ls="None", color=COL[p], markerfacecolor=COL[p], markeredgecolor="white", markeredgewidth=0.8, zorder=4)
    T(KEY_S["origin"][0], KEY_S["origin"][1], "Sleep-period T90", FS_KEY, source_file="static:V13", source_key="V13 key", source_value="Sleep-period T90", rule="text")
    T(KEY_W["origin"][0], KEY_W["origin"][1], "Wake-period T90", FS_KEY, source_file="static:V13", source_key="V13 key", source_value="Wake-period T90", rule="text")
    ov.plot([OPEN_D_X], [OPEN_MID], marker="D", ms=MS["s"], ls="None", color=BLUE, markerfacecolor="white", markeredgecolor=BLUE, markeredgewidth=1.0, zorder=4)
    ov.plot([OPEN_O_X], [OPEN_MID], marker="o", ms=MS["w"], ls="None", color=GREY, markerfacecolor="white", markeredgecolor=GREY, markeredgewidth=1.0, zorder=4)
    T(KEY_O["origin"][0], KEY_O["origin"][1], "Not significant", FS_KEY, source_file="static:V13", source_key="V13 key (open markers)", source_value="Not significant", rule="text")
    # caption, two lines centred on the axis
    XC = (AX_X0 + AX_X1) / 2
    T(XC, CAP1["origin"][1], "Hazard ratio (95% CI)", FS_CAP, ha="center", source_file="static:V13", source_key="V13 axis caption", source_value="Hazard ratio (95% CI)", rule="text")
    T(XC, CAP2["origin"][1], "per 1 SD higher T90", FS_CAP, ha="center", source_file="static:V13", source_key="V13 axis caption", source_value="per 1 SD higher T90", rule="text")
    # read-back of the matplotlib artists before saving
    n_mk = sum(1 for ln in ov.lines if ln.get_marker() not in (None, "None", "")); n_ci = sum(1 for ln in ov.lines if ln.get_marker() in (None, "None", "") and ln.get_xdata()[0] != ln.get_xdata()[1] and ln.get_ydata()[0] == ln.get_ydata()[1] and ln.get_linewidth() == 1.7)
    assert n_mk == 2 * len(ROWS) + 4 and n_ci == 2 * len(ROWS) + 2, (n_mk, n_ci)
    used = {tuple(np.round(to_rgb(ln.get_color()), 4)) for ln in ov.lines} | {tuple(np.round(to_rgb(ln.get_markerfacecolor()), 4)) for ln in ov.lines if ln.get_marker() not in (None, "None", "")} | {tuple(np.round(to_rgb(ln.get_markeredgecolor()), 4)) for ln in ov.lines if ln.get_marker() not in (None, "None", "")}
    assert used <= {tuple(np.round(to_rgb(c), 4)) for c in (GREY, BLUE, INK, "white")}, used
    for t in ov.texts: assert t.get_fontsize() >= 9.0 and L.house_ok(t.get_text())
    for t in ov.texts:
        bb = t.get_window_extent(renderer=fig.canvas.get_renderer()); x0p, x1p = bb.x0 / fig.dpi * 72, bb.x1 / fig.dpi * 72
        assert x0p >= 8 and x1p <= 484.86 + 0.5, ("text leaves panel a's box", t.get_text(), x0p, x1p)
    out_pdf = f"{WORK}/panel_a_raw.pdf"; fig.savefig(out_pdf, transparent=True); plt.close(fig)

# ------------------------------------------------------------------ the panel's own text layer: every record must be found at its origin (page frame = panel frame here)
met = L.sheet_font_metrics(IN["spans"]); L.align_font_metrics(out_pdf, f"{WORK}/panel_a_m.pdf", met)
RB = L.spans_of(f"{WORK}/panel_a_m.pdf"); taken = set()
for rec in RECS:
    cands = [(i, s) for i, s in enumerate(RB) if i not in taken and s["text"].replace("\xa0", " ") == rec["text"] and abs(s["origin"][1] - rec["baseline"]) < 0.35]
    ha = rec["ha"]; best = None
    for i, s in cands:
        sx = s["origin"][0] if ha == "left" else (s["bbox"][2] if ha == "right" else (s["bbox"][0] + s["bbox"][2]) / 2)
        if abs(sx - rec["x"]) <= 0.6: best = i; break
    assert best is not None, ("drawn string not on the panel text layer at its origin", rec["text"], rec["x"], rec["baseline"], [(s["origin"], s["bbox"]) for _, s in cands][:2])
    taken.add(best)
left = [s["text"] for i, s in enumerate(RB) if i not in taken and s["text"].strip()]; assert not left, ("panel strings not accounted for", left)
# the V13 strings of panel a (the OLD side of the delta) and the record
old_rows = {}
for lab in labs13:
    hs = sorted([s for s in hr13 if abs(s["origin"][1] - lab["origin"][1]) < 20], key=lambda s: s["origin"][1]); qs = sorted([s for s in q13 if abs(s["origin"][1] - lab["origin"][1]) < 20], key=lambda s: s["origin"][1])
    old_rows[lab["text"]] = dict(wake=(hs[0]["text"], qs[0]["text"], "Bold" in qs[0]["font"]), sleep=(hs[1]["text"], qs[1]["text"], "Bold" in qs[1]["font"]))
geom = dict(panel_box=PANEL, axis_x=[AX_X0, AX_X1], rule_y=RULE_Y, ref_x=REF_X, b_ln=B_LN, xlim=[XLO, XHI], ticks=TICKS, label_x=LAB_X, y_first=Y_FIRST, pitch=PITCH, off_wake=OFF_W, off_sleep=OFF_S,
            head_y=HEAD_Y, col_right=COL_R, qcol_right=QCOL_R, key=dict(sleep=S_HANDLE, wake=W_HANDLE, mid=KEY_MID, open_mid=OPEN_MID), colours=dict(wake=GREY, sleep=BLUE, ink=INK), marker_pt=MS, design="V13 panel a pins (row baselines, columns, key, axis calibration, ticks, captions); markers and lines as figure3A_polish.py; no ties, no grid")
json.dump(geom, open(f"{WORK}/panel_a_geometry.json", "w"), indent=1)
json.dump(dict(sheet=SHEET, panel="a", lane="V14_L4_APNEA", round="v8.2", sources=dict(A_full=dict(path=A_PATH, sha256=GATES[A_PATH]["sha256"], step=GATES[A_PATH]["step"], output_mtime=GATES[A_PATH]["output_mtime"]),
               meta=dict(path=J_PATH, sha256=GATES[J_PATH]["sha256"], step=GATES[J_PATH]["step"], output_mtime=GATES[J_PATH]["output_mtime"]), negcontrols_final=dict(path=NC_PATH, panel=NEG_CONTROLS)),
               n_analysis_set=N_A, q_threshold=QTHR, n_sig_wake=N_SIG_W, n_sig_sleep=N_SIG_S, row_rule="significant for either period (q < 0.05, death excluded from the rule), the 12 largest sleep-period estimates, plus death, plus the controls; ordered by the asleep-minus-awake gap",
               rows_main=[r.cond for r in ROWS if not r.neg], rows_controls=[r.cond for r in ROWS if r.neg], significant_not_drawn=DROPPED_SIG, values=DRAWN, drawn_records=RECS, v13_rows=old_rows,
               exact_half_strings=HALVES, stored_precision="full precision in the csv (printed half up 2 dp, q 3 dp)", geometry=geom), open(f"{VER}/{SHEET}_a_drawn.json", "w"), indent=1, ensure_ascii=False)
print(f"panel a: rows {[r.cond for r in ROWS]}; not drawn though significant: {DROPPED_SIG}; n {N_A:,}; sig wake {N_SIG_W} sleep {N_SIG_S}; ticks {TICKS}; xlim {XLO:.3f} to {XHI:.3f}; halves {len(HALVES)}; {len(RECS)} records; wrote {out_pdf}")
in13 = set(old_rows); now = {r.cond for r in ROWS}; print(f"  rows in: {sorted(now - in13)}; out: {sorted(in13 - now)}")
