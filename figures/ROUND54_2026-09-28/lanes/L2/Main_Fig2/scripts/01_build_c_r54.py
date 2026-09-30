#!$T90_PY
"""ROUND 54 (2026-09-29), lane L2: panel c of Main_Fig2 at the round-54 sizes (tick and band labels 14 pt, panel titles 15 pt, stars
14 pt bold, the band key 14 pt, the asterisk key 13 pt, the two captions 15 pt). Copy of ROUND49 .../L2/Main_Fig2/scripts/
01_build_c_r49.py (byte copy beside as 01_build_c_r49_R49_ORIGINAL.py): the data logic (the selection file's ten conditions in its order,
results_v2 graded, the figure2B y rule and tick rule, BH q across the non-control conditions, the star rule) is unchanged; the geometry
comes from work/layout.json (01_build_a_r54.py): the panel sits in the right slot below panel b as 3 columns x 4 rows of boxes
(V30: 2 rows x 5, full width), row-major order unchanged, the tenth box alone in the last row with the band key and the asterisk key
beside it, the x title centred under the block, the rotated y label at the left. Titles are centred over their box (a 15 pt title
left-aligned at the box would run off the sheet in the last column). Band width, marker size and edge, line widths, ribbon alpha,
tick length and the star offsets are the round-30 style values as in V30. usage: 01_build_c_r54.py"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import json, math, os, sys
import numpy as np, pandas as pd
import matplotlib.ticker as mticker
from matplotlib.patches import Rectangle, Polygon
from fig2_common import *

SRC = NUM
RES = hydrated(f"{SRC}/results_v2.json"); NCP_PATH = hydrated(f"{SRC}/negcontrols_v4_panel.csv")
SEL_PATH = hydrated(f"{NUM}/fig2c_selection_v8.json")
prov = {"results_v2.json": sidecar(RES), "negcontrols_v4_panel.csv": sidecar(NCP_PATH), "fig2c_selection_v8.json": sidecar(SEL_PATH)}
R = json.load(open(RES)); G = R["graded"]["outcomes"]; BAND_N = R["graded"]["band_n"]
NCP = pd.read_csv(NCP_PATH); SEL = json.load(open(SEL_PATH))
COH_N = int(json.load(open(hydrated(f"{SRC}/cohorts.json")))["bdsp"]["n_analysis"])
BANDS = ["0-1%", "1-5%", "5-10%", ">10%"]; BAND_LAB = ["0-1", "1-5", "5-10", ">10"]
conds_nc = [k for k, v in G.items() if k not in set(NCP.condition) and not v.get("negative_control")]
q_top = dict(zip(conds_nc, bh([G[k][">10%"]["p"] for k in conds_nc])))
assert sum(BAND_N.values()) == COH_N, (BAND_N, COH_N)
LAY = layout(); LC = LAY["c"]; PH = float(LAY["page_h"])
def tw(s, size, weight="normal"): return text_width_pt(None, s, size, weight)
S_TICK, S_TITLE, S_STAR, S_KEY, S_AKEY, S_CAP = SIZE["tick"], SIZE["head"], SIZE["label"], SIZE["key"], SIZE["note"], SIZE["title"]

GEO = load_geometry(); T = load_text(); C = [r for r in GEO["records"] if r["region"] == "c"]; SC = spans_in(T, "c")
def w_(b): return b[2] - b[0]
def h_(b): return b[3] - b[1]
ORDER = list(SEL["selected_labels"]); ORDER_SRC = f"{SEL_PATH} selected_labels ({SEL['rule']})"
assert len(ORDER) == 10 and all(c in G for c in ORDER), ORDER
assert ORDER == sorted(ORDER, key=lambda c: -G[c][">10%"]["hr"]), "the selection is not in top-band order"
assert all(abs(G[c][">10%"]["hr"] - e["top_band_hr"]) < 5e-4 for c, e in zip(ORDER, SEL["selected"])), "selection file and results_v2 graded disagree"
assert not set(ORDER) & set(NCP.condition) and "Obesity" not in ORDER and "Obesity hypoventilation" not in ORDER
# geometry: 3 x 4 grid from the planner
NCOL, NROW, BW4, HBOX, CGAP = int(LC["ncol"]), int(LC["nrow"]), float(LC["bw"]), float(LC["h_box"]), float(LC["col_gap"])
BW = BW4 / 4.0
BOXES = []
for n in range(len(ORDER)):
    r, j = divmod(n, NCOL); x0 = float(LC["box_x0"]) + j * (BW4 + CGAP); y0 = float(LC["first_box_y0"]) + r * float(LC["pitch_row"])
    BOXES.append((x0, y0, x0 + BW4, y0 + HBOX))
TITLE_DY, XL_DY, YT_PAD = float(LC["title_dy"]), float(LC["xl_dy"]), float(LC["yt_pad"])
YT_DY = centre_dy(pt(S_TICK))
def ylim_rule(v):
    his = [v[b]["hi"] for b in BANDS[1:]]; los = [v[b]["lo"] for b in BANDS[1:]]
    return min(0.72, min(los) * 0.9), max(1.05, max(his) * 1.22)
YLIM = {c: ylim_rule(G[c]) for c in ORDER}
def ymap(c, bx):
    lo_y, top = YLIM[c]; y0, y1 = bx[1], bx[3]
    return lambda v: y1 - (v - lo_y) / (top - lo_y) * (y1 - y0)
def yticks_rule(c):
    lo_y, top = YLIM[c]; cand = mticker.MaxNLocator(nbins=4, steps=[1, 2, 2.5, 5, 10]).tick_values(lo_y, top); span = top - lo_y
    return [float(t) for t in cand if lo_y + 0.06 * span <= t <= top]
# style: the round-30 record (the V13/V30 panel c is that build with text edits), cross-checked against the V13 geometry as round 49 did
r30 = json.load(open(hydrated(f"{paths.FIGURE_ROOT}/ROUND30_2026-09-08/V7_L2_T90/Main_Fig2/work/panel_c_drawn.json")))["style"]
MK_W, ME = float(r30["marker_w"]), float(r30["marker_edge"])
line = [r for r in C if r["op"] == "S" and hexcol(r["stroke"]) == "#ffffff"]; assert len(line) == 10; LINE_LW = line[0]["lw"]
dash = [r for r in C if r["dash"]]; assert len(dash) == 10; DASH_LW = dash[0]["lw"]
sp = [r for r in C if r["op"] == "S" and hexcol(r["stroke"]) == INK and not r["dash"]]; SPINE_LW = float(np.median([r["lw"] for r in sp]))
tk = [r for r in C if r["op"] == "S" and hexcol(r["stroke"]) == INK and r["npts"] == 2 and max(w_(r["bbox"]), h_(r["bbox"])) < 8]; TICK_LEN = float(np.median([max(w_(r["bbox"]), h_(r["bbox"])) for r in tk]))
for k, v in (("line_lw", LINE_LW), ("dash_lw", DASH_LW), ("spine_lw", SPINE_LW), ("tick_len", TICK_LEN)):
    assert abs(r30[k] - v) < 0.05, (k, r30[k], v, "V13 measurement disagrees with the round-30 record")
STAR_DX, STAR_DY = r30["star_dx"], r30["star_dy"]
sw = sorted([r for r in C if r["op"] == "f" and hexcol(r["fill"]) in LADDER and h_(r["bbox"]) < 20 and r["bbox"][1] > 1100], key=lambda r: r["bbox"][0]); assert len(sw) == 4
assert [hexcol(r["fill"]) for r in sw] == LADDER
SW_W, SW_H = float(np.mean([w_(r["bbox"]) for r in sw])), float(np.mean([h_(r["bbox"]) for r in sw]))
KEY_STAR = [s for s in SC if s["text"] == "***" and s["origin"][1] > 1100]; assert len(KEY_STAR) == 1; KEY_STAR = KEY_STAR[0]
KEY_TAIL = [s for s in SC if s["text"] == "q < 0.001"]; assert len(KEY_TAIL) == 1; KEY_TAIL = KEY_TAIL[0]
KEY_TAIL_GAP = KEY_TAIL["origin"][0] - KEY_STAR["bbox"][2]     # V13 gap between the star and its tail

for c in ORDER:
    v = G[c]; assert all(np.isfinite(v[b][k]) for b in BANDS[1:] for k in ("hr", "lo", "hi", "p")), c
    assert v["0-1%"]["hr"] == 1.0
    lo_y, top = YLIM[c]; assert all(lo_y < v[b]["lo"] and v[b]["hi"] < top for b in BANDS[1:]), (c, YLIM[c])
def star(q): return "***" if q < 0.001 else ("**" if q < 0.01 else ("*" if q < 0.05 else ""))
STARS = {c: star(q_top[c]) for c in ORDER}
KEY_SENTENCE = "q < 0.001" if all(s == "***" for s in STARS.values()) else None
if KEY_SENTENCE is None: print("NOTE: not every panel has top-band q below 0.001:", {c: (STARS[c], q_top[c]) for c in ORDER})
# the key sits in the free slots of the last row (columns 2 and 3), the x title under the block, the y label left of the block
last_r = (len(ORDER) - 1) // NCOL; free_x0 = float(LC["box_x0"]) + (len(ORDER) - last_r * NCOL) * (BW4 + CGAP) - YT_PAD - TICK_LEN - 2.0
last_y0 = float(LC["first_box_y0"]) + last_r * float(LC["pitch_row"]); KEY_Y1 = last_y0 + 0.42 * HBOX; KEY_Y2 = KEY_Y1 + 22.0
BLOCK_CX = float(LC["box_x0"]) + (NCOL * BW4 + (NCOL - 1) * CGAP) / 2
YLAB_CY = (BOXES[0][1] + BOXES[-1][3]) / 2

drawn = dict(panel="c", mode="new", source=RES, selection=SEL_PATH, provenance=prov, band_n=BAND_N, order=ORDER, order_source=ORDER_SRC, page_h=PH,
             ylim_rule_on_v8_1_data=YLIM, style=dict(band_w=BW, marker_w=MK_W, marker_edge=ME, line_lw=LINE_LW, dash_lw=DASH_LW, spine_lw=SPINE_LW, tick_len=TICK_LEN, ribbon_alpha=0.16,
                                                     title_dy=TITLE_DY, xlabel_dy=XL_DY, ytick_pad=YT_PAD, ytick_dy=YT_DY, star_dx=STAR_DX, star_dy=STAR_DY, swatch=[SW_W, SW_H]),
             top_band_q={c: float(q_top[c]) for c in ORDER}, n_conditions_bh=len(conds_nc), stars=STARS, key_sentence=KEY_SENTENCE,
             r54=dict(plus=PLUS, grid=[NCOL, NROW], box=[BW4, HBOX], box_v30=[143.78, 96.7], col_gap=CGAP, titles="centred over the box (V30: left-aligned at the box)", key_position="free slots of the last row", key_y=[KEY_Y1, KEY_Y2], free_x0=free_x0),
             panels=[], texts=[])
def T_(x, y, s, size, ha="left", weight="normal", rot=0, **rec):
    size = pt(size)   # round 54: +4.0 pt, drawn and recorded
    ov.text(x, y, s, fontsize=size, color=INK, ha=ha, va="baseline", fontweight=weight, rotation=rot, rotation_mode="anchor")
    drawn["texts"].append(drawn_rec(s, x, y, ha=ha, size=size, weight=weight, panel="c", rot=rot, **rec))
with plt.rc_context(RC):
    fig, ov = page_figure(PH)
    for c, bx in zip(ORDER, BOXES):
        f = ymap(c, bx); x0, y0, x1, y1 = bx; xs = [x0 + (j + 0.5) * BW for j in range(4)]
        for j in range(4): ov.add_patch(Rectangle((x0 + j * BW, y0), BW, y1 - y0, facecolor=LADDER[j], edgecolor="none", zorder=0))
        v = G[c]; hrs = [1.0] + [v[b]["hr"] for b in BANDS[1:]]; los = [v[b]["lo"] for b in BANDS[1:]]; his = [v[b]["hi"] for b in BANDS[1:]]
        poly = [(xs[j + 1], f(los[j])) for j in range(3)] + [(xs[j + 1], f(his[j])) for j in reversed(range(3))]
        ov.add_patch(Polygon(poly, closed=True, facecolor=INK, edgecolor="none", alpha=0.16, zorder=1))
        ov.plot([x0, x1], [f(1.0)] * 2, color=INK, lw=DASH_LW, ls=(0, (4, 3)), zorder=2, solid_capstyle="butt")
        ov.plot(xs, [f(h) for h in hrs], color="#ffffff", lw=LINE_LW, zorder=3, solid_joinstyle="round")
        ov.plot(xs, [f(h) for h in hrs], marker="o", ms=MK_W, ls="None", markerfacecolor=INK, markeredgecolor="white", markeredgewidth=ME, zorder=4)
        ov.plot([x0, x0], [y0, y1], color=INK, lw=SPINE_LW, zorder=5, solid_capstyle="projecting"); ov.plot([x0, x1], [y1, y1], color=INK, lw=SPINE_LW, zorder=5, solid_capstyle="projecting")
        for x in xs: ov.plot([x, x], [y1, y1 + TICK_LEN], color=INK, lw=SPINE_LW, zorder=5, solid_capstyle="butt")
        for j, lab in enumerate(BAND_LAB): T_(xs[j], y1 + XL_DY, lab, S_TICK, ha="center", source_key="band tick label (V13)")
        for t in yticks_rule(c):
            ov.plot([x0 - TICK_LEN, x0], [f(t), f(t)], color=INK, lw=SPINE_LW, zorder=5, solid_capstyle="butt")
            T_(x0 - YT_PAD, f(t) + YT_DY, f"{t:g}", S_TICK, ha="right", source_key=f"y tick of '{c}' (figure2B tick rule on the panel's own y limits)", source_value=f"{t:g}", rule="text", tick_value=t)
        T_((x0 + x1) / 2, y0 - TITLE_DY, c, S_TITLE, ha="center", source_file=RES, source_key=f"graded/outcomes/{c} (row from fig2c_selection_v8.json selected_labels)", source_value=c, rule="text")
        if STARS[c]: T_(xs[3] - STAR_DX, f(his[2]) - STAR_DY, STARS[c], S_STAR, ha="center", weight="bold", source_file=RES, source_key=f"stars from the BH q of graded.outcomes['{c}']['>10%'].p across the {len(conds_nc)} non-control conditions (computed, q = {float(q_top[c]):.2e})", source_value=STARS[c], rule="text", q_top=float(q_top[c]))
        drawn["panels"].append(dict(condition=c, box=list(bx), events=v["events"], ylim=YLIM[c], yticks=yticks_rule(c),
                                    bands=[dict(band=b, hr=hrs[j], lo=(None if j == 0 else los[j - 1]), hi=(None if j == 0 else his[j - 1]), x_page=round(xs[j], 3), y_hr=round(f(hrs[j]), 3),
                                                y_lo=(None if j == 0 else round(f(los[j - 1]), 3)), y_hi=(None if j == 0 else round(f(his[j - 1]), 3))) for j, b in enumerate(BANDS)],
                                    star=STARS[c], q_top=float(q_top[c])))
    T_(float(LC["ylabel_x"]), YLAB_CY, "Hazard ratio vs T90 0-1% (95% CI)", S_CAP, ha="center", rot=90, r54_anchor="centred on the vertical span of the boxes")
    T_(BLOCK_CX, float(LC["xtitle_base"]), "Time with oxygen saturation below 90%, % of the recording", S_CAP, ha="center", r54_anchor="centred on the block of boxes")
    # band key (14 pt) and the asterisk key (13 pt) in the free slots of the last row
    x = free_x0
    for j, lab in enumerate(BANDS):
        ov.add_patch(Rectangle((x, KEY_Y1 - centre_dy(pt(S_KEY)) - SW_H / 2), SW_W, SW_H, facecolor=LADDER[j], edgecolor="none", zorder=1))
        T_(x + SW_W + 5.0, KEY_Y1, lab, S_KEY, source_key="band key (V13 round-36 wording)")
        x += SW_W + 5.0 + tw(lab, pt(S_KEY)) + 12.0
    assert x <= PAGE_W, "band key runs off the sheet"
    T_(free_x0, KEY_Y2, "***", S_AKEY, weight="bold")
    T_(free_x0 + tw("***", pt(S_AKEY), "bold") + KEY_TAIL_GAP, KEY_Y2, KEY_SENTENCE or "q < 0.05", S_AKEY, r54_anchor=f"star right edge + V13 gap {KEY_TAIL_GAP:.2f} pt", source_key="star key: every drawn panel has top-band q below 0.001" if KEY_SENTENCE else "star key (not all panels below 0.001: SEE NOTE)")
    out = os.path.join(WORK, "panel_c.pdf"); fig.savefig(out, format="pdf"); plt.close(fig)
json.dump(drawn, open(os.path.join(WORK, "panel_c_drawn.json"), "w"), indent=1)
print(f"panel_c.pdf written: {len(ORDER)} panels {ORDER}, stars {set(STARS.values())}, key '{KEY_SENTENCE}', grid {NCOL} x {NROW}, box {BW4} x {HBOX}, band w {BW:.3f}, top-band q max {max(q_top[c] for c in ORDER):.2e} (BH across {len(conds_nc)})")
