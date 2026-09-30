#!$T90_PY
"""Supp_Fig21, round 40 (lane LSUPP): panel b deleted; panel a keeps everything and takes the "90% threshold" open-ring key that sat in
panel b; the letter a stays; the page box shrinks to panel a. Re-plot of the round-37 builder 01_build_s21.py (lane copy, same numbers:
the v8.1 ladder LADDER_SUMMARY.csv and ladder_peak.json, step 202, ranking_v3.csv step 105; the per-SD counts of panel b are no longer
drawn), geometry pinned to the V13 probe record (work/base_probe.json, copied). The key ring and its label sit in the free lower-right
corner of panel a (the lower-left, where panel b held them, is where the T80 marker sits). Page box: x 0 to the panel a spine plus the
sheet's right margin (V16: 484.72 - 475.2 = 9.52 pt after the panel b spine), height unchanged.
Writes work/Supp_Fig21_mpl.pdf, then 02 (fitz box + descriptor alignment as the round-37 composer) gives Supp_Fig21.pdf."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import json, os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
from matplotlib import font_manager as fm
L = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LSUPP-A"   # round 49 (Alen, 2026-09-26): every text one point larger, lane LSUPP-A (round-40 builder copied)
sys.path.insert(0, f"{L}/Supp_Fig21/scripts")
from s21_common import *          # noqa: F401,F403  (paths of the numbers files, sidecar gate, RUNGS, colours, fonts; LANE = this sheet folder)
for _f in (FONT_REG, FONT_BOLD, FONT_ITAL): fm.fontManager.addfont(hydrated(_f))
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

WORK = f"{L}/Supp_Fig21/work"; VER = f"{L}/Supp_Fig21/verify"; os.makedirs(WORK, exist_ok=True); os.makedirs(VER, exist_ok=True)
P = json.load(open(hydrated(f"{WORK}/base_probe.json")))          # the V13 probe record (round 37), never re-probed here
PAGE = P["page"]["rect"]; W_PT, H_PT = PAGE[2], PAGE[3]; assert (round(W_PT, 2), round(H_PT, 2)) == (484.72, 233.28), PAGE
D = P["drawings"]
axes_rects = [g["rect"] for g in D if g["type"] == "f" and g["fill"] == WHITE and abs((g["rect"][3] - g["rect"][1]) - 156.96) < 0.05]
assert len(axes_rects) == 2, axes_rects
AX_A, AX_B = sorted(axes_rects)
spans = P["spans"]
def ticklabels(axr, axis):
    if axis == "x":
        ls = [s for s in spans if s["size"] == 10.0 and s["text"].strip() and axr[0] - 12 < (s["bbox"][0] + s["bbox"][2]) / 2 < axr[2] + 12 and 184 < s["origin"][1] < 200]
        return sorted([(float(s["text"]), (s["bbox"][0] + s["bbox"][2]) / 2) for s in ls])
    ls = [s for s in spans if s["size"] == 10.0 and s["text"].strip() and abs(s["bbox"][2] - (axr[0] - 6.5)) < 1.0]
    return sorted([(float(s["text"]), s["origin"][1]) for s in ls])
def tickmarks(axr, axis):
    if axis == "x": return sorted({g["rect"][0] for g in D if g["type"] == "s" and abs(g["rect"][1] - axr[3]) < 0.01 and abs(g["rect"][3] - axr[3] - 3.0) < 0.01 and axr[0] <= g["rect"][0] <= axr[2]})
    return sorted({g["rect"][1] for g in D if g["type"] == "s" and abs(g["rect"][2] - axr[0]) < 0.01 and abs(g["rect"][0] - axr[0] + 3.0) < 0.01})
def limits(axr, axis):
    labs = ticklabels(axr, axis); marks = tickmarks(axr, axis); vals = [v for v, _ in labs]
    assert len(vals) == len(marks), (axis, labs, marks)
    if axis == "x":
        c1, c0 = np.polyfit(vals, marks, 1); res = float(np.abs(c0 + c1 * np.array(vals) - np.array(marks)).max()); lo, hi = (axr[0] - c0) / c1, (axr[2] - c0) / c1
    else:
        marks = sorted(marks, reverse=True); c1, c0 = np.polyfit(vals, marks, 1); res = float(np.abs(c0 + c1 * np.array(vals) - np.array(marks)).max()); lo, hi = (axr[3] - c0) / c1, (axr[1] - c0) / c1
    assert res < 0.01, (axis, res)
    return vals, (float(lo), float(hi)), marks
XT_A, XLIM_A, XM_A = limits(AX_A, "x"); YT_A, YLIM_A, YM_A = limits(AX_A, "y")
XLIM = (round(XLIM_A[0], 3), round(XLIM_A[1], 3)); YLIM_A = (round(YLIM_A[0], 6), round(YLIM_A[1], 6))
assert XLIM == (77.9, 96.6) and YLIM_A == (0.008, 0.0222) and YT_A == [0.008, 0.012, 0.016, 0.02] and XT_A == [80.0, 85.0, 88.0, 90.0, 92.0, 95.0], (XLIM, YLIM_A, YT_A, XT_A)
letters = {s["text"]: s for s in P["letters"]}; assert set(letters) == {"a", "b"} and letters["a"]["size"] == 13.0 and letters["a"]["font"] == "Arial-BoldMT"
def keyring():
    r = [g for g in D if g["type"] == "s" and g["kinds"] == "c" and g["color"] == BLUE and abs((g["rect"][2] - g["rect"][0]) - 10.488) < 0.02]; assert r, "key ring not found"
    return ((r[0]["rect"][0] + r[0]["rect"][2]) / 2, (r[0]["rect"][1] + r[0]["rect"][3]) / 2)
KEY_RING_V13 = keyring()      # (322.452, 160.38): 11.41 pt right of panel b's spine, 19.62 pt above its floor

# ------------------------------------------------------------------ the numbers (asserted before drawing; panel a only)
X = [80, 85, 88, 90, 92, 95]; PUBLISHED_X = 90
prov = {os.path.basename(p): sidecar(p) for p in (LADDER_CSV, PEAK_JSON, RANKING_CSV)}
S = pd.read_csv(hydrated(LADDER_CSV)); PEAK = json.load(open(hydrated(PEAK_JSON))); R = pd.read_csv(hydrated(RANKING_CSV), comment="#")
prim = S[S.threshold.isin(RUNGS)].copy(); assert len(prim) == 6, prim.threshold.tolist()
prim["x"] = prim.threshold.str[1:].astype(int); prim = prim.sort_values("x").reset_index(drop=True); assert prim.x.tolist() == X
DC = prim.dC_mean.astype(float).tolist(); NS = prim.n.astype(int).tolist()
row90 = R[R.feature == "spo2_pct_below_90"]; assert len(row90) == 1, row90
RULE = float(row90.dC.iloc[0]); assert int(row90["rank"].iloc[0]) == 1, row90.to_dict("records")
t90_dc = DC[RUNGS.index("T90")]; assert abs(RULE - t90_dc) < 1e-8, ("the T90 gain of ranking_v3 is not the ladder's T90 rung", RULE, t90_dc)
assert PEAK["peak_by_dC"] == RUNGS[int(np.argmax(DC))] and abs(PEAK["peak_dC"] - max(DC)) < 5e-10, (PEAK, DC)
i = int(np.argmax(DC)); assert all(DC[k] < DC[k + 1] for k in range(i)) and all(DC[k] > DC[k + 1] for k in range(i, 5)), DC
assert all(YLIM_A[0] < v < YLIM_A[1] for v in DC + [RULE]), (YLIM_A, DC, RULE)
# the round-37 record of the same numbers (positive control: the drawn values equal the record's to 1e-9)
R37 = json.load(open(hydrated(f"{WORK}/Supp_Fig21_drawn_R37.json")))
r37_a = [v for v in R37["values"] if v["panel"] == "a" and v["kind"] == "marker"]
assert [round(v["value"], 9) for v in r37_a] == [round(v, 9) for v in DC], "the ladder gains differ from the round-37 record"
assert abs([v for v in R37["values"] if v["kind"] == "dashed rule"][0]["value"] - RULE) < 1e-9
src_a = f"{LADDER_CSV} column dC_mean"; src_rule = f"{RANKING_CSV} row feature=spo2_pct_below_90 column dC"

# ------------------------------------------------------------------ geometry of the new page box: panel a plus the sheet's right margin
RIGHT_MARGIN = W_PT - AX_B[2]                 # 9.52 pt after the panel b spine on V16
NEW_W = round(AX_A[2] + RIGHT_MARGIN, 3)      # 236.16 + 9.52 = 245.68
# the key ring in panel a: the lower-right corner is free (the curve rises from T80 at the lower left to T90/T92 near the top and
# falls to T95 at mid height on the right); the ring keeps its V13 height above the floor (19.62 pt) and the label its 1.8 rung offset
KEY_DY_PT = AX_B[3] - KEY_RING_V13[1]         # 19.62 pt above the axes floor
KEY_LABEL_DX = 3.1 - 1.3                      # data units between the ring centre and the label start on V13 (panel b)
def data_x(px): return XLIM[0] + (px - AX_A[0]) / (AX_A[2] - AX_A[0]) * (XLIM[1] - XLIM[0])
def data_y(py): return YLIM_A[0] + (AX_A[3] - py) / (AX_A[3] - AX_A[1]) * (YLIM_A[1] - YLIM_A[0])
KEY_X_PT = 150.0                              # ring centre, page x (the label then spans about 166 to 228 pt, inside the spine at 236.16)
KEY_RING_XY = (data_x(KEY_X_PT), data_y(AX_A[3] - KEY_DY_PT)); KEY_LABEL_XY = (KEY_RING_XY[0] + KEY_LABEL_DX, KEY_RING_XY[1])

PT_PLUS = 1.0   # round 49: every text one point larger (axis titles 11 -> 12, tick labels 10 -> 11, the rule and key labels 9.5 -> 10.5)
LAB, TCK, ANN, PANEL = 11.0 + PT_PLUS, 10.0 + PT_PLUS, 9.5 + PT_PLUS, 13.0 + PT_PLUS
W_IN, H_IN = 171.0 / 25.4, 3.24
# round 49: at 12 pt the two-line y title, pushed left by the 11 pt tick labels, ran 2.05 pt past the page's left edge (the V26 sheet had 2.04 pt of
# room). Two nudges of that element restore the V26 clearance (1.95 pt): the y title 3 pt closer to its tick labels (labelpad 6 -> 3) and the
# y tick labels 1 pt closer to the axis (ytick.major.pad 3.5 -> 2.5). The x axis, the axes box and every data mark are untouched.
Y_LABELPAD, Y_TICKPAD = 3.0, 2.5
RC = {"font.family": "Arial", "font.size": TCK, "axes.labelsize": LAB, "xtick.labelsize": TCK, "ytick.labelsize": TCK, "ytick.major.pad": Y_TICKPAD, "legend.fontsize": 10.0,
      "axes.linewidth": 0.8, "axes.edgecolor": INK, "text.color": INK, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
      "xtick.major.width": 0.8, "ytick.major.width": 0.8, "xtick.major.size": 3.0, "ytick.major.size": 3.0, "xtick.direction": "out", "ytick.direction": "out",
      "axes.spines.top": False, "axes.spines.right": False, "axes.facecolor": WHITE, "figure.facecolor": WHITE, "savefig.facecolor": WHITE, "savefig.bbox": "standard",
      "savefig.pad_inches": 0.0, "figure.dpi": 300, "savefig.dpi": 300, "pdf.fonttype": 42, "ps.fonttype": 42, "axes.unicode_minus": False, "legend.frameon": False}
def frac_box(axr):
    x0, ytop, x1, ybot = axr; return [x0 / 72 / W_IN, (H_PT - ybot) / 72 / H_IN, (x1 - x0) / 72 / W_IN, (ybot - ytop) / 72 / H_IN]
def circle(ax, axr, xlim, ylim, xv, yv, diam_pt, face, edge, lw, z):
    w = diam_pt * (xlim[1] - xlim[0]) / (axr[2] - axr[0]); h = diam_pt * (ylim[1] - ylim[0]) / (axr[3] - axr[1])
    e = Ellipse((xv, yv), w, h, facecolor=face, edgecolor=edge, lw=lw, zorder=z, clip_on=False); ax.add_patch(e); return e
def page_xy(axr, xlim, ylim, xv, yv):
    x0, ytop, x1, ybot = axr; return (x0 + (xv - xlim[0]) / (xlim[1] - xlim[0]) * (x1 - x0), ybot - (yv - ylim[0]) / (ylim[1] - ylim[0]) * (ybot - ytop))
drawn = dict(sheet="Supp_Fig21", mode="round 49: +1 pt on the round-40 sheet (panel a only, the key ring in panel a, page box shrunk to panel a)", pt_plus=PT_PLUS, text_pt=dict(axis_title=LAB, tick=TCK, annotation=ANN), y_title_nudge=dict(labelpad_v26=6.0, labelpad=Y_LABELPAD, ytick_pad_v26=3.5, ytick_pad=Y_TICKPAD), provenance=prov, page_pt_v16=[W_PT, H_PT], page_pt_new=[NEW_W, H_PT],
             canvas_pt=[W_IN * 72, H_IN * 72], axes_a=AX_A, axes_b_v16=AX_B, xlim=XLIM, ylim_a=YLIM_A, key_ring_v16_page_xy=list(KEY_RING_V13), values=[], texts=[])
def rec(**kw): drawn["values"].append(kw)
def tx(**kw): drawn["texts"].append(kw)
with plt.rc_context(RC):
    fig = plt.figure(figsize=(W_IN, H_IN)); ax_a = fig.add_axes(frac_box(AX_A))
    ax_a.axhline(RULE, color=INK, lw=0.9, ls=(0, (4, 3)), zorder=1)
    ax_a.plot(X, DC, "-", color=BLUE, lw=1.7, solid_capstyle="round", zorder=2)
    for xv, yv in zip(X, DC): circle(ax_a, AX_A, XLIM, YLIM_A, xv, yv, 6.0, BLUE, WHITE, 0.8, 3)
    circle(ax_a, AX_A, XLIM, YLIM_A, PUBLISHED_X, DC[X.index(PUBLISHED_X)], 150 ** 0.5, "none", BLUE, 1.4, 4)
    ax_a.annotate(GAIN_LABEL, (XLIM[0] + 0.35, RULE), textcoords="offset points", xytext=(0, 5), fontsize=ANN, color=INK, ha="left", va="bottom", zorder=5)
    ax_a.set_xlim(*XLIM); ax_a.set_ylim(*YLIM_A); ax_a.set_xticks(X); ax_a.set_yticks(YT_A); ax_a.set_yticklabels([f"{t:.3f}" for t in YT_A])
    ax_a.set_xlabel("Saturation threshold defining\nthe exposure, %", fontsize=LAB, labelpad=4); ax_a.set_ylabel("Gain in held-out concordance\nover age and sex", fontsize=LAB, labelpad=Y_LABELPAD)
    # R40 follow-up (audit): no panel letter on a single-panel sheet (as ED 4, 8, 10)
    # the key of panel b, now in panel a (same ring size 110 ** 0.5 pt, same 1.4 pt blue stroke, same label, 9.5 pt)
    circle(ax_a, AX_A, XLIM, YLIM_A, KEY_RING_XY[0], KEY_RING_XY[1], 110 ** 0.5, "none", BLUE, 1.4, 4)
    ax_a.annotate(KEY_LABEL, KEY_LABEL_XY, fontsize=ANN, color=INK, ha="left", va="center", zorder=5)
    for t in ax_a.xaxis.get_major_ticks() + ax_a.yaxis.get_major_ticks(): t.tick1line.set_visible(False); t.tick2line.set_visible(False)
    dy = 3.0 * (YLIM_A[1] - YLIM_A[0]) / (AX_A[3] - AX_A[1]); dx = 3.0 * (XLIM[1] - XLIM[0]) / (AX_A[2] - AX_A[0])
    for xv in X: ax_a.plot([xv, xv], [YLIM_A[0], YLIM_A[0] - dy], color=INK, lw=0.8, solid_capstyle="butt", clip_on=False, zorder=1)
    for yv in ax_a.get_yticks(): ax_a.plot([XLIM[0], XLIM[0] - dx], [yv, yv], color=INK, lw=0.8, solid_capstyle="butt", clip_on=False, zorder=1)
    for k, xv in enumerate(X):
        px, py = page_xy(AX_A, XLIM, YLIM_A, xv, DC[k])
        rec(panel="a", kind="marker", label=f"T{xv} gain", threshold=xv, value=DC[k], page_xy=[round(px, 3), round(py, 3)], source=src_a, key=f"threshold={RUNGS[k]}", n=NS[k])
    px, py = page_xy(AX_A, XLIM, YLIM_A, XLIM[0], RULE)
    rec(panel="a", kind="dashed rule", label=GAIN_LABEL, value=RULE, page_y=round(py, 3), source=src_rule, key="spo2_pct_below_90.dC", text=GAIN_LABEL)
    px, py = page_xy(AX_A, XLIM, YLIM_A, PUBLISHED_X, DC[X.index(PUBLISHED_X)])
    rec(panel="a", kind="open ring", label=KEY_LABEL, threshold=PUBLISHED_X, value=DC[X.index(PUBLISHED_X)], page_xy=[round(px, 3), round(py, 3)], source="design: the ring marks the 90% threshold (T90)")
    px, py = page_xy(AX_A, XLIM, YLIM_A, *KEY_RING_XY)
    rec(panel="a", kind="key ring", label=KEY_LABEL, page_xy=[round(px, 3), round(py, 3)], v16_page_xy=[round(KEY_RING_V13[0], 3), round(KEY_RING_V13[1], 3)], source="design (the V16 key of panel b, moved into panel a)", text=KEY_LABEL)
    for xv in X: tx(panel="a", kind="x tick label", text=str(xv), value=xv, source=LADDER_CSV, key=f"threshold=T{xv}", rule="text")
    for t in YT_A: tx(panel="a", kind="y tick label", text=f"{t:.3f}", value=t, source="static:V13", key="fixed tick rule of the V13 sheet (limits 0.008 to 0.0222)", rule="text")
    tx(panel="a", kind="rule label", text=GAIN_LABEL, value=RULE, source=RANKING_CSV, key="feature=spo2_pct_below_90:dC (label rides 5 pt above the rule)", rule="text")
    tx(panel="a", kind="key label", text=KEY_LABEL, value=90, source="static:V16 (the key of panel b, moved)", rule="text")
    fig.canvas.draw()
    # the key must not touch the curve, the markers or the rule label: geometric check in page points
    rend = fig.canvas.get_renderer(); kb = [t for t in ax_a.texts if t.get_text() == KEY_LABEL][0].get_window_extent(renderer=rend)
    kx0, kx1 = kb.x0 / fig.dpi * 72, kb.x1 / fig.dpi * 72; ky_top, ky_bot = H_PT - kb.y1 / fig.dpi * 72, H_PT - kb.y0 / fig.dpi * 72
    ring_px = page_xy(AX_A, XLIM, YLIM_A, *KEY_RING_XY)
    assert kx1 <= AX_A[2] - 2.0, ("the key label crosses the panel a spine", kx1, AX_A[2])
    for k, xv in enumerate(X):
        mx, my = page_xy(AX_A, XLIM, YLIM_A, xv, DC[k])
        assert ((mx - ring_px[0]) ** 2 + (my - ring_px[1]) ** 2) ** 0.5 > 6.0 + 3.0 + 4.0, ("a marker touches the key ring", xv)
        assert not (kx0 - 4 <= mx <= kx1 + 4 and ky_top - 4 <= my <= ky_bot + 4), ("a marker sits under the key label", xv)
    # the polyline segments stay clear of the key box (distance of each segment to the box corners and centre)
    pts = [page_xy(AX_A, XLIM, YLIM_A, xv, DC[k]) for k, xv in enumerate(X)]
    def seg_dist(p, a, b):
        ax_, ay_ = a; bx_, by_ = b; px_, py_ = p; vx, vy = bx_ - ax_, by_ - ay_; L2 = vx * vx + vy * vy or 1.0
        t = max(0.0, min(1.0, ((px_ - ax_) * vx + (py_ - ay_) * vy) / L2)); return ((px_ - (ax_ + t * vx)) ** 2 + (py_ - (ay_ + t * vy)) ** 2) ** 0.5
    box_pts = [(kx0 - 6, ky_top), (kx1, ky_top), (kx0 - 6, ky_bot), (kx1, ky_bot), ring_px]
    for a_, b_ in zip(pts, pts[1:]):
        for q in box_pts: assert seg_dist(q, a_, b_) > 5.0, ("the curve passes through the key", a_, b_, q)
    rule_y = page_xy(AX_A, XLIM, YLIM_A, XLIM[0], RULE)[1]; assert ky_top > rule_y + 8, ("the key sits on the T90 gain rule", ky_top, rule_y)
    drawn["key_box_pt"] = [round(kx0, 2), round(ky_top, 2), round(kx1, 2), round(ky_bot, 2)]
    for t in fig.findobj(matplotlib.text.Text):
        s = t.get_text()
        if not t.get_visible() or not s.strip(): continue
        assert t.get_fontsize() >= 9.0, (s, t.get_fontsize())
        for bad in ("—", ";", "published"): assert bad not in s, s
    out = f"{WORK}/Supp_Fig21_mpl.pdf"
    fig.savefig(out, metadata={"Title": "Supp_Fig21", "Creator": "LSUPP-A round 49 Supp_Fig21 r40_build_s21.py +1 pt (matplotlib, Arial Type 42)"}); plt.close(fig)
drawn["output"] = out; drawn["new_page_box"] = [0, 0, NEW_W, H_PT]
json.dump(drawn, open(f"{VER}/Supp_Fig21_drawn.json", "w"), indent=1)
print(f"wrote {out}\n  DC {DC}\n  RULE {RULE}\n  key ring at page {drawn['values'][-1]['page_xy']} (V16 panel b {KEY_RING_V13}), label box {drawn['key_box_pt']}\n  new page box {drawn['new_page_box']}")
