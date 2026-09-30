#!$T90_PY
"""ROUND 54 (2026-09-29), lane L2: panel b of Main_Fig2 at the round-54 sizes (row labels, key, ticks, the control header 14 pt, the
caption 15 pt, 'no estimate' 13 pt). Copy of ROUND49 .../L2/Main_Fig2/scripts/01_build_b_r49.py (byte copy beside as
01_build_b_r49_R49_ORIGINAL.py): the data logic (top-12 rule, controls, lags, no-estimate slots, the row-set report) is unchanged; the
geometry comes from work/layout.json (01_build_a_r54.py, the planner): the panel sits in the right slot above panel c, the key in two
rows beside the letter (3 + 2 entries, V30 reading order), 18 row slots at a 27 pt pitch with the four lag sub-rows 6.44 pt apart (the
V30 spacing; marker sizes, line widths, colours, the x range 0.58 to 4.3 and the tick set unchanged), the axis about 256 pt wide.
usage: 01_build_b_r54.py"""
import json, math, os, sys
import numpy as np
import matplotlib.ticker as mticker
from matplotlib.transforms import blended_transform_factory
from fig2_common import *

SRC = NUM; TAG = ""
DATA = hydrated(f"{SRC}/lag_ladder.json"); prov = {"lag_ladder.json": sidecar(DATA)}
LAY = layout(); LB = LAY["b"]; PH = float(LAY["page_h"])
COH_N = int(json.load(open(hydrated(f"{SRC}/cohorts.json")))["bdsp"]["n_analysis"])
G = load_geometry(); T = load_text(); B = [r for r in G["records"] if r["region"] == "b"]; SB = spans_in(T, "b")
def w_(b): return b[2] - b[0]
def h_(b): return b[3] - b[1]
def tw(s, size, weight="normal", style="normal"): return text_width_pt(None, s, size, weight, style)
XTICKS = [0.7, 1.0, 1.5, 2.0, 3.0, 4.0]
XLO, XHI = 0.62, 4.30; XLO_V13, XHI_V13 = XLO, XHI; AXIS_EXTENDED = False
GUTTER_PAD = float(LB["gutter_pad"])
XLABEL = "Hazard ratio for a new diagnosis per 1 SD of T90 (95% CI)"; CTRL_HEAD = "Negative controls"
PRIMARY, COMPARE = BLUE, GREY
CI_LW, MEDGE = 1.7, 0.8
MARK = ["o", "s", "D", "^"]; SIZE_S = [28.0, 23.0, 26.0, 30.0]; MSIZE = [6.0, 5.4, 5.7, 6.3]
LAGS = ["0", "1", "2", "5"]; LAGNAME = ["All follow-up", "1-year landmark", "2-year landmark", "5-year landmark"]; QUIET_NAME = "5-year landmark, too few events"
TICK_PT, TITLE_PT, FLOOR_PT = pt(SIZE["tick"]), pt(SIZE["title"]), pt(SIZE["note"])   # 14, 15, 13
# the V13 key marks: line length and marker sizes kept
kln = [r for r in B if r["op"] == "S" and h_(r["bbox"]) < 0.05 and r["bbox"][0] < 640 and r["bbox"][1] < 110]; assert len(kln) == 4
KEY_LINE = float(np.mean([w_(r["bbox"]) for r in kln])); assert abs(KEY_LINE - 16.0) < 0.05, KEY_LINE
KEY_TXT_GAP, KEY_ENTRY_GAP = 5.5, 12.0       # V13: text origin 5.5 pt after the line end
key_rows = [[("o", MSIZE[0], False, LAGNAME[0]), ("s", MSIZE[1], False, LAGNAME[1]), ("D", MSIZE[2], False, LAGNAME[2])],
            [("^", MSIZE[3], False, LAGNAME[3]), ("^", MSIZE[3], True, QUIET_NAME)]]
KEY_ROWS = []
for ry, row in zip(LB["key_y"], key_rows):
    x = float(LB["key_x0"]); yc = ry - centre_dy(TICK_PT)
    for mk, ms, is_open, txt in row:
        line = None if is_open else (x, x + KEY_LINE)
        mx = x + KEY_LINE / 2; tx = x + KEY_LINE + KEY_TXT_GAP
        KEY_ROWS.append(dict(marker=mk, ms=ms, xy=(mx, yc), open=is_open, text=txt, tx=tx, ty=ry, line=line))
        x = tx + tw(txt, TICK_PT) + KEY_ENTRY_GAP
assert KEY_ROWS[-1]["tx"] + tw(QUIET_NAME, TICK_PT) <= PAGE_W - 4 and KEY_ROWS[2]["tx"] + tw(LAGNAME[2], TICK_PT) <= PAGE_W - 4, "key runs off the sheet"

# ------------------------------------------------------------------ data, asserted against the file (row set reported against the sheet)
J = json.load(open(DATA)); PC = J["per_condition"]
assert J["lags"] == [0, 1, 2, 5] and J["controls_n"] == 5 and J["cohort_n"] == COH_N, (J["lags"], J["controls_n"], J["cohort_n"], COH_N)
panel = [k for k, v in PC.items() if not v["negative_control"] and not v["circular"] and v["ci_lag0"][0] > 1]
panel.sort(key=lambda k: -PC[k]["hr_lag0"]); print(f"conditions with a lag-0 lower bound above 1 (non-control, non-circular): {len(panel)}")
TOP = panel[:12]
CTRL = sorted([k for k, v in PC.items() if v["negative_control"]], key=lambda k: -PC[k]["hr_lag0"]); assert len(CTRL) == 5
LABEL_SHORT = {"Ventricular arrhythmia or cardiac arrest": "Ventricular arrhythmia/arrest"}   # the V13 label column could not hold the full name (owner question)
LAB = {k: LABEL_SHORT.get(PC[k]["condition"], PC[k]["condition"]) for k in TOP + CTRL}; E5 = {k: int(PC[k]["events_by_lag"]["5"]) for k in TOP + CTRL}
noest = sorted((k, L) for k in TOP + CTRL for L in LAGS if PC[k]["hr_by_lag"][L] is None)
_los = [PC[k]["ci_by_lag"][L][0] for k in TOP + CTRL for L in LAGS if PC[k]["hr_by_lag"][L] is not None]; _his = [PC[k]["ci_by_lag"][L][1] for k in TOP + CTRL for L in LAGS if PC[k]["hr_by_lag"][L] is not None]
if min(_los) <= XLO or max(_his) >= XHI:
    XLO = min(XLO, 0.58); XHI = max(XHI, max(_his) * 1.02); AXIS_EXTENDED = True
    print(f"NOTE: x axis extended to {XLO:.3f}-{XHI:.3f} (V13 {XLO_V13:.3f}-{XHI_V13:.3f}), same ticks (as round 37 and V30)")
R49B = json.load(open(os.path.join(WORK, "r49_panel_b_drawn.json"))); assert [round(v, 6) for v in R49B["xlim"]] == [round(XLO, 6), round(XHI, 6)], (R49B["xlim"], XLO, XHI)
ROWS = TOP + [None] + CTRL; N = len(ROWS); assert N == 18 == int(LB["n_slots"])
base_labels = sorted(s["text"] for s in SB if s["font"] == "ArialMT" and abs(s["size"] - 10) < 0.01 and "(" in s["text"] and 110 < s["origin"][1] < 700)
new_labels = sorted(f"{LAB[k]} ({E5[k]:,})" for k in TOP + CTRL)
label_delta = dict(lost=sorted(set(base_labels) - set(new_labels)), gained=sorted(set(new_labels) - set(base_labels)))
base_noest = len([s for s in SB if s["text"] == "no estimate"])
print("row labels against the V13 sheet:", "identical set" if not label_delta["lost"] and not label_delta["gained"] else label_delta, f"| no-estimate slots: sheet {base_noest}, file {len(noest)} {noest}")
new_order = [f"{LAB[k]} ({E5[k]:,})" for k in TOP + CTRL]
# ------------------------------------------------------------------ geometry (planner values)
WIDEST = max(tw(f"{LAB[k]} ({E5[k]:,})", TICK_PT) for k in TOP + CTRL)
AX_L = float(LB["x0"]) + 6.0 + WIDEST + GUTTER_PAD; AX_R = float(LB["ax_r"]); AX_T, AX_B = float(LB["ax_t"]), float(LB["ax_b"])
PITCH = float(LB["pitch"]); assert abs((AX_B - AX_T) / N - PITCH) < 1e-6
SUB = float(LB["sub_row_spacing"]) / PITCH            # sub-row spacing in row units (6.44 pt)
OFFS = [1.5 * SUB, 0.5 * SUB, -0.5 * SUB, -1.5 * SUB]  # lag 0, 1, 2, 5 from the row centre (V30: 0.30, 0.10, -0.10, -0.30 at a 32.2 pt pitch, the same 6.44 pt)
assert AX_R - AX_L >= 200, (AX_L, AX_R)
def xpos(hr): return AX_L + (math.log(hr) - math.log(XLO)) / (math.log(XHI) - math.log(XLO)) * (AX_R - AX_L)
NOEST_W = tw("no estimate", FLOOR_PT, style="italic")
def marks_near(y_page, reach=8.5):
    """x extents (with the marker half width) of every estimable interval of any row whose sub-row centre lies within reach of y_page."""
    out = []
    for kk in TOP + CTRL:
        i = ROWS.index(kk); yc = N - i - 0.5
        for k2, L2 in enumerate(LAGS):
            h2 = PC[kk]["hr_by_lag"][L2]
            if h2 is None: continue
            yp = AX_T + (N - (yc + OFFS[k2])) * PITCH
            if abs(yp - y_page) <= reach: lo2, hi2 = PC[kk]["ci_by_lag"][L2]; out.append((min(xpos(lo2), xpos(h2) - 4.0), max(xpos(hi2), xpos(h2) + 4.0)))
    return out
NOEST_SIDES = {}
def noest_side(key, y_page):
    """'no estimate' sits at the left end of the axis (V30) unless a mark of a neighbouring sub-row lies under it, then at the right end."""
    left = (xpos(XLO * 1.02), xpos(XLO * 1.02) + NOEST_W); right = (xpos(XHI / 1.02) - NOEST_W, xpos(XHI / 1.02))
    near = marks_near(y_page)
    def clear(span): return all(span[1] < a - 2.0 or span[0] > b + 2.0 for a, b in near)
    side = "left" if clear(left) else ("right" if clear(right) else None)
    assert side is not None, ("no clear place for 'no estimate'", key, near)
    NOEST_SIDES[key] = side; return side
CAPTION_BASE = float(LB["caption_base"]); TICK_BASE = float(LB["tick_base"])
CAPTION_X = min((AX_L + AX_R) / 2, PAGE_W - 12.0 - tw(XLABEL, TITLE_PT) / 2)   # the 15 pt caption (388 pt) centred under the axis would leave the page: as far right as a 12 pt margin allows
drawn = []
with plt.rc_context({**RC, "axes.labelsize": TITLE_PT, "xtick.labelsize": TICK_PT, "ytick.labelsize": TICK_PT}):
    fig = plt.figure(figsize=(PAGE_W / 72.0, PH / 72.0)); fig.patch.set_alpha(0.0)
    ax = fig.add_axes([AX_L / PAGE_W, (PH - AX_B) / PH, (AX_R - AX_L) / PAGE_W, (AX_B - AX_T) / PH]); ax.set_facecolor("#ffffff")
    shade = True
    for i, key in enumerate(ROWS):
        if key is None: continue
        if shade: ax.axhspan(N - i - 1.0, N - i, color=GRID, lw=0, zorder=0)
        shade = not shade
    ib = ROWS.index(None)
    for ya, yb in ((0.0, N - ib - 1.0), (N - ib, float(N))): ax.plot([1.0, 1.0], [ya, yb], color=INK, lw=0.9, ls=(0, (4, 3)), zorder=1, solid_capstyle="butt")
    gut = blended_transform_factory(fig.dpi_scale_trans, ax.transData); ticks_y = []
    for i, key in enumerate(ROWS):
        y0 = N - i - 0.5
        if key is None:
            ax.text((AX_L - GUTTER_PAD) / 72.0, y0 - 0.15 - centre_dy(TICK_PT) / PITCH, CTRL_HEAD, fontsize=TICK_PT, color=INK, va="baseline", ha="right", fontweight="bold", transform=gut, clip_on=False, zorder=6); continue
        ticks_y.append((y0, f"{LAB[key]} ({E5[key]:,})"))
        for k, L in enumerate(LAGS):
            y = y0 + OFFS[k]; y_page = AX_T + (N - y) * PITCH
            hr = PC[key]["hr_by_lag"][L]
            if hr is None:
                side = noest_side(key, y_page)
                ax.text(XLO * 1.02 if side == "left" else XHI / 1.02, y - centre_dy(FLOOR_PT) / PITCH, "no estimate", fontsize=FLOOR_PT, color=COMPARE, va="baseline", ha=side, style="italic")
                drawn.append(dict(row=i, key=key, condition=LAB[key], lag=int(L), estimable=False, y_page=round(y_page, 3), side=side)); continue
            lo, hi = PC[key]["ci_by_lag"][L]; assert XLO < lo and hi < XHI, (key, L, lo, hi)
            quiet = (L == "5") and not bool(PC[key]["lag5_informative"])
            ax.plot([lo, hi], [y, y], color=PRIMARY, lw=CI_LW, zorder=2)
            ax.scatter([hr], [y], s=SIZE_S[k], marker=MARK[k], zorder=3, linewidths=MEDGE, facecolors="white" if quiet else PRIMARY, edgecolors=PRIMARY if quiet else "white")
            drawn.append(dict(row=i, key=key, condition=LAB[key], lag=int(L), estimable=True, hr=hr, lo=lo, hi=hi, quiet=quiet, marker=MARK[k],
                              x_page=round(xpos(hr), 3), x_lo=round(xpos(lo), 3), x_hi=round(xpos(hi), 3), y_page=round(y_page, 3), events=int(PC[key]["events_by_lag"][L])))
    ys, labs = zip(*ticks_y); ax.set_yticks(ys); ax.set_yticklabels(labs, fontsize=TICK_PT); ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0, pad=GUTTER_PAD)
    ax.set_xscale("log"); ax.set_xticks(XTICKS); ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:g}")); ax.xaxis.set_minor_locator(mticker.NullLocator())
    ax.set_xlim(XLO, XHI); ax.set_ylim(0.0, float(N))
    fig.text(CAPTION_X / PAGE_W, 1.0 - CAPTION_BASE / PH, XLABEL, fontsize=TITLE_PT, color=INK, ha="center", va="baseline")
    ov = fig.add_axes([0, 0, 1, 1]); ov.set_xlim(0, PAGE_W); ov.set_ylim(PH, 0); ov.axis("off"); ov.patch.set_visible(False)
    for r in KEY_ROWS:
        if r["line"]: ov.plot(list(r["line"]), [r["xy"][1]] * 2, color=PRIMARY, lw=CI_LW, solid_capstyle="projecting", zorder=2)
        ov.plot([r["xy"][0]], [r["xy"][1]], marker=r["marker"], ms=r["ms"], ls="None", zorder=3, markerfacecolor="white" if r["open"] else PRIMARY,
                markeredgecolor=PRIMARY if r["open"] else "white", markeredgewidth=MEDGE)
        ov.text(r["tx"], r["ty"], r["text"], fontsize=TICK_PT, color=INK, ha="left", va="baseline")
    out = os.path.join(WORK, "panel_b.pdf"); fig.savefig(out, format="pdf"); plt.close(fig)
# the matplotlib x tick labels sit at their own baseline: record it from the axes (tick length 3 + pad 3.5 + cap height)
XTICK_BASE = AX_B + 3.0 + 3.5 + CAP * TICK_PT
info = dict(panel="b", mode="new", source=DATA, provenance=prov, axes_box_page=[AX_L, AX_T, AX_R, AX_B], pitch=round(PITCH, 4), sub_row_pt=float(LB["sub_row_spacing"]), offsets=OFFS, xlim=[XLO, XHI], xlim_v13=[XLO_V13, XHI_V13], axis_extended=AXIS_EXTENDED, xticks=XTICKS,
            tick_x=[xpos(v) for v in XTICKS], caption_baseline=CAPTION_BASE, rows=[LAB[k] if k else CTRL_HEAD for k in ROWS], row_labels=new_order, base_row_labels_v13=sorted(base_labels), label_delta_v13=label_delta,
            noest=noest, base_noest=base_noest, key_rows=[dict(r, xy=list(r["xy"])) for r in KEY_ROWS], drawn=drawn,
            r54=dict(plus=PLUS, tick_pt=TICK_PT, title_pt=TITLE_PT, floor_pt=FLOOR_PT, widest_label=WIDEST, key_line=KEY_LINE, caption_x=CAPTION_X, caption_rule="min(axis centre, page width - 12 - half the caption width)", key_form="two rows beside the letter, 3 + 2 entries in the V30 reading order", noest_sides=NOEST_SIDES),
            model=J["provenance"]["model"], lag5_verdict=dict(n_informative=J["lag5_verdict"]["n_informative"], survive=J["lag5_verdict"]["survive_among_informative"]))
recs = []
for i, key in enumerate(ROWS):
    y0 = N - i - 0.5; base_y = AX_T + (N - y0) * PITCH + centre_dy(TICK_PT)
    if key is None:
        recs.append(drawn_rec(CTRL_HEAD, AX_L - GUTTER_PAD, AX_T + (N - (y0 - 0.15 - centre_dy(TICK_PT) / PITCH)) * PITCH, ha="right", size=TICK_PT, weight="bold", panel="b")); continue
    txt = f"{LAB[key]} ({E5[key]:,})"; cond = PC[key]["condition"]
    recs.append(drawn_rec(txt, AX_L - GUTTER_PAD, base_y, ha="right", size=TICK_PT, panel="b", source_file=DATA,
                          source_key=f"label: per_condition[{key}].condition + ' (' + events_by_lag[5] + ')'" + (" [display label shortened, owner question]" if cond != LAB[key] else ""),
                          source_value=txt, rule="text", condition=cond, events_5y=E5[key]))
    for k, L in enumerate(LAGS):
        if PC[key]["hr_by_lag"][L] is None:
            sd = NOEST_SIDES[key]
            recs.append(drawn_rec("no estimate", xpos(XLO * 1.02) if sd == "left" else xpos(XHI / 1.02), AX_T + (N - (y0 + OFFS[k])) * PITCH + centre_dy(FLOOR_PT), ha=sd, size=FLOOR_PT, panel="b", source_file=DATA, source_key=f"'no estimate' where per_condition[{key}].hr_by_lag[{L}] is null (at the {sd} end of the axis: the left end is under a neighbouring mark)" if sd == "right" else f"'no estimate' where per_condition[{key}].hr_by_lag[{L}] is null", source_value="no estimate", rule="text"))
recs.append(drawn_rec(XLABEL, CAPTION_X, CAPTION_BASE, ha="center", size=TITLE_PT, panel="b"))
for r in KEY_ROWS: recs.append(drawn_rec(r["text"], r["tx"], r["ty"], ha="left", size=TICK_PT, panel="b"))
for v in XTICKS: recs.append(drawn_rec(f"{v:g}", xpos(v), XTICK_BASE, ha="center", size=TICK_PT, panel="b", source_key="x tick (V30 axis)"))
info["texts"] = recs
json.dump(info, open(os.path.join(WORK, "panel_b_drawn.json"), "w"), indent=1)
print(f"panel_b.pdf written: {sum(1 for d in drawn if d['estimable'])} markers, {sum(1 for d in drawn if not d['estimable'])} no-estimate rows, axes {AX_L:.1f}-{AX_R:.1f} x {AX_T:.1f}-{AX_B:.1f}, pitch {PITCH:.2f} (sub-rows {LB['sub_row_spacing']} pt), top-12 {TOP}, controls {CTRL}")
