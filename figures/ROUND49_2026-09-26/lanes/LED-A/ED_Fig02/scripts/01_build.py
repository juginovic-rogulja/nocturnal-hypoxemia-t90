#!/usr/bin/env python3
"""ROUND 49 lane LED-A copy (2026-09-26): DELTA = +1 pt inside text() (ED2_DELTA=0 reproduces the round-38 build for the positive control),
one row label kept at 10 pt (KEEP_SIZE), the seven note strings of the two printed notes not drawn (round 40 LNOTES removal re-applied).
Original: ROUND38_2026-09-16/figures/lanes/ED_Fig02/scripts/01_build.py.
ED_Fig02, lane V7_L2_T90 (relaunch 3), step 1: the whole sheet re-plotted from the v7 numbers, pinned to the base sheet's geometry.

One matplotlib figure of the page size with one overlay axes whose data units are page points (y down). Every plot element is
drawn as an explicit patch or line in page points (no scatter, no markers, so the PDF has no XObjects, like the base sheet):
  panels a, b, c   x = a + PERLN ln(v) fitted to the base sheet's tick marks (work/base_geometry.json, positive control 0.001 pt),
                   19 rows at the base row centres, grey square (before) 3.744 pt above the row centre, blue interval and circle
                   (after) 3.744 pt below, square 4.796 pt with 0.8 pt grey edge, circle 5.292 pt with 0.8 pt white edge, interval 1.7 pt,
                   dashed reference rule at 1 (0.9 pt, dashes 3.6 on 2.7 off), bottom spine 0.8 pt with 3 pt ticks at 0.8 1 1.2 1.5 2,
                   the #eef0f1 control band at the base rectangle.
  text             every string at the base sheet's own baseline origin (ha left, va baseline), the base sizes 10, 9.5, 8.43, 8.0,
                   Arial regular and bold as Type 42; printed 1-year and 5-year values from the NEW file, NE and the non-informative
                   5-year values in grey #8a9099, the five control labels grey.
  key column       the visible key items as patches at the base rectangles (the hidden older copies are not reproduced).
  letters a, b, c  NOT drawn here: stamped by 02_compose.py with fitz TextWriter (Arial Bold 13 pt) at the base origins.
Data: ed2_common.derive(NUM) (the generator's row rule), asserted equal to work/data_new.json and inside the base x-limits before drawing.
v8.2 (2026-09-15 evening): panel c from the re-extracted sleep_unadjusted_v1.json (step 125, 19:29), all five controls fitted (no NE);
the row set follows the rule on the v8.2 files (the four new outcomes are now candidates). Backup of the v8.1 copy: 01_build_PRE_V8_2.py.
Writes work/sheet_mpl.pdf and verify/ED_Fig02_drawn.json.  usage: 01_build.py [--old]  (--old draws the pre-v7 snapshot for the eye)
"""
import json, os, sys, gc, math
DELTA = float(os.environ.get("ED2_DELTA", "1.0"))                  # round 49: +1 pt on every text element (0.0 for the OLD replot positive control)
KEEP_SIZE = {"Ventricular arrhythmia or cardiac arrest"}            # round 49 exception: at 11 pt this row label (185.8 -> 203.1 pt) would run into its own row's interval and square (195.3 to 222.8 pt); kept at 10 pt
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties, fontManager
from matplotlib.patches import Rectangle, Circle, Ellipse
from matplotlib.lines import Line2D
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ed2_common import *
from v14lib import halfup

OLDMODE = "--old" in sys.argv
SRC, TAG = (OLD, "_old") if OLDMODE else (NUM, "")
for f in (ARIAL, ARIAL_B): fontManager.addfont(f)
FP = {"regular": FontProperties(fname=ARIAL), "bold": FontProperties(fname=ARIAL_B)}
plt.rcParams.update({"pdf.fonttype": 42, "lines.scale_dashes": False, "savefig.bbox": "standard", "savefig.pad_inches": 0.0,
                     "figure.facecolor": "white", "savefig.facecolor": "white", "axes.unicode_minus": False, "pdf.compression": 6,
                     "path.simplify": False, "lines.solid_capstyle": "butt", "lines.dash_capstyle": "butt"})

G = load_geometry(); W, H = G["page"]
D = derive(SRC, "old" if OLDMODE else "new")
if not OLDMODE:
    ref = json.load(open(hydrated(os.path.join(WORK, "data_new.json"))))
    assert ref["rows"] == D["rows"] and ref["sha256"] == D["sha256"], "the build does not see the data step's files"
    assert D["c_missing"] == [], ("v8.2: every control must have a sleep-period estimate", D["c_missing"])
P = json.loads(json.dumps(G["panels"])); XLO, XHI = P["a"]["xlo"], P["a"]["xhi"]; TICKS_V13 = {n: list(P[n]["ticks_x"]) for n in P}
for n in "bc": assert abs(P[n]["xlo"] - XLO) < 1e-4 and abs(P[n]["xhi"] - XHI) < 1e-4 and abs(P[n]["perln"] - P["a"]["perln"]) < 1e-3, n
ROWC = G["row_centres"]; LANE_HALF = G["lane_half"]; SQ = G["square_side"]; CD = G["circle_diam"]; CIW = G["ci_w"]
assert len(ROWC) == 19 == len(D["rows"])
# ------------------------------------------------------------------ assertions before drawing
vals = [(rw["condition"], p, k, rw[p][k]) for rw in D["rows"] for p in "abc" if rw[p] is not None for k in ("before", "after", "lo", "hi")]
bad = [v for v in vals if not (XLO < v[3] < XHI)]
XLO_V13, XHI_V13 = XLO, XHI; AXIS_EXTENDED = False
if bad:
    # v8.1: an interval falls outside the V13 x range (inguinal hernia, panel b, lower bound 0.696 under 0.7034): the three panels share
    # one log mapping, extended to the left to the generator's own rule (0.96 x the smallest drawn value); the tick values stay, the
    # tick marks move right, the reference rules move with x(1). Declared in the record and the report.
    XLO = min(XLO, float(D["generator_xlim"][0])); XHI = max(XHI, float(D["generator_xlim"][1])) if any(v[3] >= XHI_V13 for v in bad) else XHI; AXIS_EXTENDED = True
    for n in "abc":
        p_ = P[n]; x0, x1 = p_["spine"][0], p_["spine"][2]
        p_["perln"] = (x1 - x0) / (math.log(XHI) - math.log(XLO)); p_["a"] = x0 - p_["perln"] * math.log(XLO)
        p_["ticks_x"] = [p_["a"] + p_["perln"] * math.log(v) for v in p_["tick_values"]]; p_["ref_x"] = p_["a"]; p_["xlo"], p_["xhi"] = XLO, XHI
    print(f"NOTE: x axis extended to {XLO:.4f}-{XHI:.4f} (V13 {XLO_V13:.4f}-{XHI_V13:.4f}) in all three panels, same tick values: {[v[:3] for v in bad]}")
    bad = [v for v in vals if not (XLO < v[3] < XHI)]
assert not bad, f"values outside the x-limits {XLO:.4f}..{XHI:.4f}: {bad}"
base_labels = [s["text"] for s in sorted([s for s in G["visible_spans"] if abs(s["origin"][0] - 12.96) < 0.05 and abs(s["size"] - 10.0) < 0.01 and s["text"] != "Negative controls"], key=lambda s: s["origin"][1])]
ROWSET_DELTA = dict(entering=sorted({rw["condition"] for rw in D["rows"]} - set(base_labels)), leaving=sorted(set(base_labels) - {rw["condition"] for rw in D["rows"]}))
print("row set against V13: entering", ROWSET_DELTA["entering"], "leaving", ROWSET_DELTA["leaving"], "(declared in the record)")
assert [rw["condition"] for rw in D["rows"][14:]] == D["controls"]
for rw in D["rows"]:
    assert rw["a"]["y1_text"] == halfup(rw['a']['y1'], 2)
    assert rw["a"]["y5_text"] == ("NE" if rw["a"]["y5"] is None else halfup(rw['a']['y5'], 2))
print(f"rows: {[rw['condition'] for rw in D['rows']]}")
print(f"drawn range {D['drawn_min']:.4f}..{D['drawn_max']:.4f} inside {XLO:.5f}..{XHI:.5f}")
# per-row baselines from the base sheet (index-aligned: row i of the new sheet takes row i's baselines)
lab_spans = sorted([s for s in G["visible_spans"] if abs(s["origin"][0] - 12.96) < 0.05 and abs(s["size"] - 10.0) < 0.01 and s["text"] != "Negative controls"], key=lambda s: s["origin"][1])
c1_spans = sorted([s for s in G["visible_spans"] if abs(s["origin"][0] - 248.69) < 0.05 and abs(s["size"] - 9.5) < 0.01 and s["text"] != "1-year"], key=lambda s: s["origin"][1])
c5_spans = sorted([s for s in G["visible_spans"] if abs(s["origin"][0] - 286.13) < 0.05 and abs(s["size"] - 9.5) < 0.01 and s["text"] != "5-year"], key=lambda s: s["origin"][1])
assert len(lab_spans) == len(c1_spans) == len(c5_spans) == 19
for i in range(19):   # the base baselines sit at a fixed offset from the row centres
    assert abs(lab_spans[i]["origin"][1] - ROWC[i] - 3.582) < 0.06, (i, lab_spans[i]["origin"][1] - ROWC[i])
    assert abs(c1_spans[i]["origin"][1] - ROWC[i] - 2.582) < 0.06 and abs(c5_spans[i]["origin"][1] - ROWC[i] - 2.582) < 0.06, i

# ------------------------------------------------------------------ the page
fig = plt.figure(figsize=(W / 72.0, H / 72.0), dpi=72)
ov = fig.add_axes([0, 0, 1, 1]); ov.set_xlim(0, W); ov.set_ylim(H, 0); ov.axis("off"); ov.patch.set_visible(False)
drawn = []
def rec(**kw): drawn.append(kw)
def text(x, y, s, size, color=INK, weight="regular", **meta):
    size = size + (0.0 if s in KEEP_SIZE else DELTA)                 # round 49: the +1 pt inside the helper; KEEP_SIZE strings stay at their V26 size
    ov.text(x, y, s, fontsize=size, color=color, fontproperties=FP[weight], ha="left", va="baseline", zorder=6)
    rec(kind="text", text=s, x=round(x, 3), baseline=round(y, 3), size=size, color=color, weight=weight, **meta)

SHA = D["sha256"]; SC = D.get("sidecars", {})
def src(file, key, value):
    return dict(source_file=f"{SRC}/{file}", key=key, value=value, sha256=SHA[file], sidecar=SC.get(file))

for n in "abc":
    p = P[n]; xf = xmap(p); x0, ys, x1 = p["spine"]
    if p.get("band"): ov.add_patch(Rectangle((p["band"][0], p["band"][1]), p["band"][2] - p["band"][0], p["band"][3] - p["band"][1], fc=BAND, ec="none", lw=0, zorder=0))   # V13: no control shading
    ov.add_line(Line2D([p["ref_x"], p["ref_x"]], [p["ref_y"][0], p["ref_y"][1]], color=INK, lw=p["ref_w"], dashes=(3.6, 2.7), zorder=1))
    ov.add_line(Line2D([x0, x1], [ys, ys], color=INK, lw=p["spine_w"], solid_capstyle="projecting", zorder=2))
    for tx in p["ticks_x"]: ov.add_line(Line2D([tx, tx], [ys, ys + p["tick_len"]], color=INK, lw=p["tick_w"], solid_capstyle="butt", zorder=2))
    for i, rw in enumerate(D["rows"]):
        v = rw[n]; yc = ROWC[i]; ysq, ycir = yc - LANE_HALF, yc + LANE_HALF
        if v is None:   # panel c, group P file without this (new) control: NE in grey at the reference rule, no marker
            text(p["ref_x"] - 6.5, ycir + 3.0, "NE", 9.5, color=GREY, what="panel c NE (no sleep-period estimate in the v7 file, group P)", row=i, condition=rw["condition"], panel=n)
            rec(kind="ne", panel=n, row=i, condition=rw["condition"]); continue
        keys = {"a": (f"lag_ladder.csv condition={rw['condition']} lag_years=0 hr", f"lag_ladder.csv condition={rw['condition']} lag_years=2 hr", f"lag_ladder.csv condition={rw['condition']} lag_years=2 lo", f"lag_ladder.csv condition={rw['condition']} lag_years=2 hi"),
                "b": (f"causal_tests_all.csv condition={rw['condition']} bmi_sub_hr", f"causal_tests_all.csv condition={rw['condition']} bmi_adj_hr", f"causal_tests_all.csv condition={rw['condition']} bmi_adj_lo", f"causal_tests_all.csv condition={rw['condition']} bmi_adj_hi"),
                "c": (f"sleep_unadjusted_v1.json outcomes.{rw['condition']}.sleep_alone.hr", f"sleep_unadjusted_v1.json outcomes.{rw['condition']}.sleep_adjusted_for_awake.hr", f"sleep_unadjusted_v1.json outcomes.{rw['condition']}.sleep_adjusted_for_awake.lo", f"sleep_unadjusted_v1.json outcomes.{rw['condition']}.sleep_adjusted_for_awake.hi")}[n]
        fil = {"a": "lag_ladder.csv", "b": "causal_tests_all.csv", "c": "sleep_unadjusted_v1.json"}[n]
        ov.add_line(Line2D([xf(v["lo"]), xf(v["hi"])], [ycir, ycir], color=BLUE, lw=CIW, solid_capstyle="round", zorder=3))
        rec(kind="interval", panel=n, row=i, condition=rw["condition"], x0=round(xf(v["lo"]), 3), x1=round(xf(v["hi"]), 3), y=round(ycir, 3), lo=v["lo"], hi=v["hi"], key_lo=keys[2], key_hi=keys[3], **src(fil, keys[2], v["lo"]))
        ov.add_patch(Rectangle((xf(v["before"]) - SQ / 2, ysq - SQ / 2), SQ, SQ, fc=PALE, ec=GREY, lw=G["square_edge_w"], zorder=4))
        rec(kind="square", panel=n, row=i, condition=rw["condition"], cx=round(xf(v["before"]), 3), cy=round(ysq, 3), side=SQ, **src(fil, keys[0], v["before"]))
        ov.add_patch(Circle((xf(v["after"]), ycir), CD / 2, fc=BLUE, ec="white", lw=G["circle_edge_w"], zorder=5))
        rec(kind="circle", panel=n, row=i, condition=rw["condition"], cx=round(xf(v["after"]), 3), cy=round(ycir, 3), diam=CD, **src(fil, keys[1], v["after"]))
        if n == "a":
            text(12.96, lab_spans[i]["origin"][1], rw["condition"], 10.0, color=GREY if rw["control"] else INK, what="row label", row=i, condition=rw["condition"])
            text(248.69, c1_spans[i]["origin"][1], rw["a"]["y1_text"], 9.5, what="1-year value", row=i, condition=rw["condition"],
                 **src("lag_ladder.csv", f"lag_ladder.csv condition={rw['condition']} lag_years=1 hr", rw["a"]["y1"]))
            text(286.13, c5_spans[i]["origin"][1], rw["a"]["y5_text"], 9.5, color=GREY if rw["a"]["y5_grey"] else INK, what="5-year value", row=i, condition=rw["condition"],
                 grey_reason=(None if not rw["a"]["y5_grey"] else ("estimable false" if not rw["a"]["y5_estimable"] else "lag5_informative false")),
                 **src("lag_ladder.csv", f"lag_ladder.csv condition={rw['condition']} lag_years=5 hr" + ("" if rw["a"]["y5_estimable"] else " (estimable false, prints NE)"), rw["a"]["y5"]))
# fixed strings at the base origins (design constants, identical words)
def fixed(txt, size, weight="regular", color=INK, **sel):
    s = span_at(G, txt, size, **sel); text(s["origin"][0], s["origin"][1], txt, s["size"], color=color, weight=weight, what="design constant (base sheet)")
fixed("Negative controls", 10.0, "bold"); fixed("1-year", 9.5); fixed("5-year", 9.5)
for n, p in P.items():
    for tx, tv, tx13 in zip(p["ticks_x"], p["tick_values"], G["panels"][n]["ticks_x"] if not AXIS_EXTENDED else TICKS_V13[n]):
        lab = f"{tv:g}"; s = span_at(G, lab, 10.0, xmin=tx13 - 12, xmax=tx13 + 12, ymin=425, ymax=436)
        text(s["origin"][0] + (tx - tx13), s["origin"][1], lab, 10.0, what="tick label" + (" (axis extended, moved with its tick)" if AXIS_EXTENDED else ""), panel=n, tick_value=tv)
for txt, xmin, xmax in (("Rate of new diagnosis", 130, 250), ("Rate of new diagnosis", 320, 440), ("Rate of new diagnosis", 455, 575),
                        ("per 1 SD (95% CI)", 130, 250), ("per 1 SD (95% CI)", 320, 440),
                        ("per 1 SD of sleep-period", 455, 575), ("time below 90%", 455, 575), ("saturation (95% CI)", 455, 575)):
    fixed(txt, 10.0, xmin=xmin, xmax=xmax)
for let in ("a", "b", "c"): fixed(let, 8.43, "bold", xmin=580, xmax=600)   # the key letters (x 587.2), at their V13 y (the round-31 edit moved the key column)
for txt in ("No landmark", "2-year landmark,", "with 95% CI", "Same patients,", "no adjustment", "Body mass index",
            "in the model", "Sleep oxygen", "alone", "Awake oxygen in", "the model"):   # V13 (round 31): no 'Negative / controls, shaded' key entry
    fixed(txt, 8.43, ymax=200)
NOTES_REMOVED = ("NE = no estimate, too few", "events remain after the", "5-year landmark", "Grey value = the 5-year", "re-fit has too few events",
                 "for that condition, not a", "negative finding")   # round 40 (LNOTES): the two printed notes left for the legend; not drawn (text-only removal re-applied)
if OLDMODE:                                                          # the OLD replot is the positive control against V13, which still carries the notes
    for txt in NOTES_REMOVED: fixed(txt, 8.0, ymax=275)
# the visible key markers, found by geometry on the V13 probe (key column x > 575): three pale squares, three blue handle lines,
# three blue circles, top to bottom = the a, b, c key rows. The V13 key draws the handle lines at 1.36 pt, squares' edges 0.8,
# circle rings 0.64 (the round-30 base's 80-percent key widths), the a and b handle lines inside a 7.28 pt clip; the V13 sheet
# (a round-31 stream edit of the round-30 build) keeps that construction, so the same widths are drawn.
KD = sorted(G["key_drawings"], key=lambda d: (d["rect"][1], d["rect"][0]))
def near_(a, b, tol=0.002): return a is not None and b is not None and all(abs(x - y) < tol for x, y in zip(a, b))
k_sq = sorted([d for d in KD if d["kinds"] == "re" and near_(d["fill"], [0.8, 0.8196, 0.8392])], key=lambda d: d["rect"][1])
k_ln = sorted([d for d in KD if d["kinds"] in ("l", "ll") and near_(d["color"], [0.0078, 0.5333, 0.8196])], key=lambda d: d["rect"][1])
k_ci = sorted([d for d in KD if d["kinds"] == "cccccccc" and near_(d["fill"], [0.0078, 0.5333, 0.8196])], key=lambda d: d["rect"][1])
assert len(k_sq) == len(k_ln) == len(k_ci) == 3, (len(k_sq), len(k_ln), len(k_ci))
KEY_LW, KEY_SQ_LW, KEY_CIRC_LW = 1.36, 0.8, 0.64
for j, (dsq, dln, dci) in enumerate(zip(k_sq, k_ln, k_ci)):
    r = dsq["rect"]; ov.add_patch(Rectangle((r[0], r[1]), r[2] - r[0], r[3] - r[1], fc=PALE, ec=GREY, lw=KEY_SQ_LW, zorder=4)); rec(kind="key square", rect=r, lw=KEY_SQ_LW)
    r = dln["rect"]; ln = Line2D([r[0], r[2]], [r[1], r[3]], color=BLUE, lw=KEY_LW, solid_capstyle="projecting", zorder=3); ov.add_line(ln)
    rec(kind="key line", rect=r, lw=KEY_LW, cap="projecting", v13_width=dln.get("width"))
    r = dci["rect"]; ov.add_patch(Ellipse(((r[0] + r[2]) / 2, (r[1] + r[3]) / 2), r[2] - r[0], r[3] - r[1], fc=BLUE, ec="white", lw=KEY_CIRC_LW, zorder=5)); rec(kind="key circle", rect=r, lw=KEY_CIRC_LW)
out_pdf = os.path.join(WORK, f"sheet_mpl{TAG}.pdf")
fig.savefig(out_pdf, format="pdf"); plt.close(fig); gc.collect()
meta = dict(builder=os.path.abspath(__file__), base_sheet=BASE, base_sha256=sha256(BASE), source_root=SRC, sha256=SHA, sidecars=SC, rows=[rw["condition"] for rw in D["rows"]], rowset_delta=ROWSET_DELTA, c_missing=D["c_missing"], excluded_by_sleep_file=D["excluded_by_sleep_file"],
            xlim=[XLO, XHI], xlim_v13=[XLO_V13, XHI_V13], axis_extended=AXIS_EXTENDED, delta_pt=DELTA, keep_size=sorted(KEEP_SIZE), notes_removed=list(NOTES_REMOVED), generator_xlim_on_these_numbers=D["generator_xlim"], note=("x range extended to the left (declared), same ticks" if AXIS_EXTENDED else "x-limits kept at the V13 sheet's"), panels_used=P,
            page=[W, H], n_text=sum(1 for d in drawn if d["kind"] == "text"), n_printed_values=sum(1 for d in drawn if d.get("what") in ("1-year value", "5-year value")))
json.dump(dict(meta=meta, drawn=drawn), open(os.path.join(VER, f"ED_Fig02_drawn{TAG}.json"), "w"), indent=1)
print("wrote", out_pdf, "and", os.path.join(VER, f"ED_Fig02_drawn{TAG}.json"), "| text items", meta["n_text"], "| printed values", meta["n_printed_values"])
