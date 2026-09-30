#!/usr/bin/env python3
"""Round 49 LED-B copy (2026-09-26): every text +1 pt (l5a14.Sheet), outputs into work/. Otherwise the round-37 builder. ED Fig 10 (single panel, flat sheet, nine organ-system bars), V14: re-plotted from scratch at the V13 geometry from the v8.1
treatment_v3.json (step 211, the comparison ADJUSTED for prevalent cardiopulmonary disease, FIX A1, decision 6) through the organ grouping of numbers/disease_definitions.py (the v8 module). V13 design: bars only,
no percentage labels beside the bars (round 31 LTEXT), "Group (responding of conditions)" row labels, ticks -10 to 40, axis title.
Lineage: ROUND30 V7_L5a_PAP_MAIN/ED_Fig10/scripts/run_ed10.py (read-back only; the panel was never rebuilt in round 30) and the
generator FINAL_FIGURES_2026-08-14/_scripts/r2_eFigureNEW_organ_rollup.py (design)."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from l5a14 import *

SHEET = "ED_Fig10"; D = f"{LANE}/{SHEET}"; W = f"{D}/work"; V = f"{D}/verify"; os.makedirs(f"{V}/crops", exist_ok=True)
BASE = f"{V13}/{SHEET}.pdf"; OUT = f"{W}/{SHEET}_built.pdf"   # round 49: the built sheet (V13 page size) goes to work/, the runner adds the title strip
ck = Checks(f"Lane V14_L5a_PAP_MAIN, {SHEET}, {time.strftime('%Y-%m-%d %H:%M')}. BASE (design) = V13 {BASE} (sha256 {sha256(BASE)[:16]}). NEW = v8.1 numbers with sidecars.")

# ---- 0 gates
sc_t = gate("treatment_v3"); sc_g = gate("gaps_v2"); dump({"treatment_v3": sc_t, "gaps_v2": sc_g}, f"{W}/sidecars.json")
ck.log("treatment_v3.json passes the v8.1 sidecar gate (step 211, inputs data_frozen_v8_2026-09, sha equal); the ADJUSTED fields are rolled up (FIX A1)", sc_t["v8_inputs"], f"out {sc_t['output_mtime']} sha {sc_t['sha256'][:16]}")
G = load_json(SPEC["gaps_v2"])
ck.log("gaps_v2.json v8.1 (step 118) and carries no organ roll-up block, so ED 10 derives from treatment_v3.json (adjusted) through disease_definitions.ORGAN_GROUP", sc_g["v8_inputs"] and not any("organ" in k.lower() for k in G), f"keys {sorted(G)}")

# ---- 1 the V13 geometry (flat sheet)
spans, info = text_spans(BASE); spans = collapse(spans); items = drawings(BASE)
dump({"info": info, "spans": spans}, f"{W}/base_text.json"); dump(items, f"{W}/base_drawings.json")
PW, PH = info["rect"][2], info["rect"][3]
ck.log("V13 sheet flat (0 XObjects), page 484.72 x 277.20", info["xobjects"] == 0 and abs(PW - 484.72) < 0.05 and abs(PH - 277.2) < 0.05, info)
labels = column(spans, 0, 20, size=10.0); assert len(labels) == 9, labels
LAB_X = float(np.mean([s["x"] for s in labels])); LAB_Y = [s["y"] for s in labels]
pitch = float(np.mean(np.diff(LAB_Y))); assert np.std(np.diff(LAB_Y)) < 0.05, np.diff(LAB_Y)
R = sorted([r for r in rects(items) if r["fill"] in (BLUE, GREY) and 8 <= r["h"] <= 20], key=lambda r: r["y0"]); assert len(R) == 9, len(R)
BAR_H = float(np.mean([r["h"] for r in R])); BAR_Y0 = [r["y0"] for r in R]
bar_off = [r["y0"] - s["y"] for r, s in zip(R, labels)]; assert max(bar_off) - min(bar_off) < 0.05, bar_off
BAR_OFF = float(np.mean(bar_off))
zero_x = sorted(set(round(r["x0"], 2) for r in R if r["fill"] == BLUE) | set(round(r["x1"], 2) for r in R if r["fill"] == GREY)); assert max(zero_x) - min(zero_x) < 0.05, zero_x
ZERO_X = float(np.mean(zero_x))
spine = [h for h in hlines(items, min_len=100) if h["stroke"] == INK]; assert len(spine) == 1, spine
SPINE = spine[0]; SPINE_Y = SPINE["y"]
zr = [v for v in vlines(items, max_len=300) if v["stroke"] == INK and (v["y1"] - v["y0"]) > 100]; assert len(zr) == 1, zr
ZR = zr[0]
ticks = [v for v in vlines(items, max_len=6.0) if v["stroke"] == INK and abs(v["y0"] - SPINE_Y) < 0.05]; ticks = sorted(ticks, key=lambda v: v["x"]); assert len(ticks) == 6, ticks
TICK_LEN = float(np.mean([v["y1"] - v["y0"] for v in ticks])); TICK_VALS = [-10, 0, 10, 20, 30, 40]
a, b, res = fit_axis(TICK_VALS, [v["x"] for v in ticks], logscale=False); assert float(np.abs(res).max()) < 0.05, res
ck.log(f"V13 geometry read: 9 row labels at x {LAB_X:.2f} pitch {pitch:.3f}, bars {BAR_H:.3f} pt tall at label y {BAR_OFF:+.3f}, zero x {ZERO_X:.3f}, linear axis {b:.4f} pt per percent (ticks -10..40, residual {np.abs(res).max():.3f}), spine y {SPINE_Y:.3f}", abs(a - ZERO_X) < 0.05, f"a {a:.3f}")
tick_labels = {t["text"]: t for t in spans if abs(t["y"] - 245.4) < 1.0}; assert set(tick_labels) == {"−10", "0", "10", "20", "30", "40"}, tick_labels
TICK_LAB_Y = float(np.mean([t["y"] for t in tick_labels.values()]))
title = find_span(spans, "Median rate lower where restored, %"); assert title is not None
axes_bg = [r for r in rects(items) if r["fill"] == WHITE and r["w"] > 300 and r["h"] > 200 and r["x0"] > 100]; assert len(axes_bg) == 1, axes_bg
AX = axes_bg[0]
old_labels = [s["text"] for s in labels]

# ---- 2 the v8.1 roll-up
grp, ungrouped, mapping = organ_rollup_adj(); dump({"groups": grp, "ungrouped": ungrouped, "mapping": mapping}, f"{W}/organ_rollup_v8_1.json")
T, rows = treatment_rows_adj(); n_ctrl = sum(1 for _, v in rows if v["negative_control"])
ck.log(f"roll-up from the v8.1 file (ADJUSTED comparison): {len(grp)} groups, {sum(g['conditions'] for g in grp)} conditions, {sum(g['controls_held_out'] for g in grp)} controls held out (the file carries {n_ctrl} controls), {sum(g['responding'] for g in grp)} responding, ungrouped {ungrouped}",
       len(grp) == 9 and sum(g["conditions"] for g in grp) + sum(g["controls_held_out"] for g in grp) + len(ungrouped) == len(rows) and sorted(ungrouped) == ["Cardiovascular composite", "Death from any cause"],
       f"{[(g['group'], g['conditions'], g['controls_held_out'], g['responding'], round(g['median_pct_lower'], 2)) for g in grp]}")
NEW_MAP = {c: g for g in grp for c in g["members"]}
new_out = [c for c in mapping if c in ("Ventricular arrhythmia or cardiac arrest", "Interstitial lung disease", "Polycythemia", "Parkinson's disease")]
ck.log("the v8 outcomes present in the PAP file have an organ group from disease_definitions.py (no typed map in this lane)", all(c in NEW_MAP for c in new_out), f"{[(c, NEW_MAP.get(c)) for c in new_out]}")
vals = [g["median_pct_lower"] for g in grp]
inside_v13 = all(AX["x0"] <= axis_x(a, b, v, False) <= AX["x1"] for v in vals)
# FIX A1: the adjusted Liver median (-16.6) falls below the V13 axis (-10 to 40): extend the tick range downward in steps of 10 and rescale
# it inside the V13 tick span (leftmost tick at the old -10 position, 40 at its old position), a declared axis change; the box, labels,
# bar heights, fonts and colours are the V13 design
AXIS_CHANGED = not inside_v13
if AXIS_CHANGED:
    lo_new = int(np.floor(min(vals) / 10.0) * 10); hi_new = int(np.ceil(max(vals) / 10.0) * 10) if max(vals) > 40 else 40
    TICK_VALS_NEW = list(range(lo_new, hi_new + 1, 10)); x_left, x_right = ticks[0]["x"], ticks[-1]["x"]
    b_new = (x_right - x_left) / (hi_new - lo_new); a_new = x_left - b_new * lo_new
else:
    TICK_VALS_NEW = TICK_VALS; a_new, b_new = a, b
TICK_X_NEW = [a_new + b_new * v for v in TICK_VALS_NEW]; ZERO_X_NEW = a_new
ck.log(f"axis: {'EXTENDED to ' + str(TICK_VALS_NEW[0]) + '..' + str(TICK_VALS_NEW[-1]) + ' and rescaled inside the V13 tick span (' + f'{b_new:.4f} pt per percent, zero x {ZERO_X_NEW:.2f}), a declared design change forced by the adjusted medians' if AXIS_CHANGED else 'the V13 ticks -10..40, unchanged'}; every median inside the axes box",
       all(AX["x0"] <= a_new + b_new * v <= AX["x1"] for v in vals), f"x range {a_new + b_new * min(vals):.1f} to {a_new + b_new * max(vals):.1f} in [{AX['x0']:.1f}, {AX['x1']:.1f}]; medians {[round(v, 1) for v in vals]}")

# ---- 3 draw
S = Sheet(PW, PH)
S.rect(AX["x0"], AX["y0"], AX["x1"], AX["y1"], fill=WHITE)
for i, g in enumerate(grp):
    y0 = LAB_Y[0] + i * pitch + BAR_OFF; xe = a_new + b_new * g["median_pct_lower"]
    if g["median_pct_lower"] >= 0: S.rect(ZERO_X_NEW, y0, xe, y0 + BAR_H, fill=BLUE)
    else: S.rect(xe, y0, ZERO_X_NEW, y0 + BAR_H, fill=GREY)
S.line(ZERO_X_NEW, ZR["y0"], ZERO_X_NEW, ZR["y1"], INK, ZR["width"], cap=0)
S.line(SPINE["x0"], SPINE_Y, SPINE["x1"], SPINE_Y, INK, SPINE["width"], cap=2)
TICK_W = ticks[0]["width"]
for v, tx in zip(TICK_VALS_NEW, TICK_X_NEW):
    S.line(tx, SPINE_Y, tx, SPINE_Y + TICK_LEN, INK, TICK_W, cap=0)
    lab = ("−" + str(-v)) if v < 0 else str(v)
    S.text(tx, TICK_LAB_Y, lab, 10.0, align="center", panel="single", source_key="axis tick" + (" (axis extended, FIX A1)" if AXIS_CHANGED else ""))
S.text((AX["x0"] + AX["x1"]) / 2, title["y"], title["text"], title["size"], align="center", panel="single", source_key="axis title")
for i, g in enumerate(grp):
    S.text(LAB_X, LAB_Y[0] + i * pitch, f"{g['group']} ({g['responding']} of {g['conditions']})", 10.0, panel="single", source_file=SPEC["treatment_v3"],
           source_key=f"outcomes[members of {g['group']}].adj_p < .05 and .adj_hr > 1 (responding) of {g['conditions']} (disease_definitions.ORGAN_GROUP)",
           source_value=[g["group"], g["responding"], g["conditions"]], rule="group_k_of_n", note="recomputed roll-up, members " + "; ".join(g["members"]))
S.save(OUT)
drawn = S.drawn; dump(drawn, f"{W}/{SHEET}_drawn.json")

# ---- 4 verify
sp2, info2 = text_spans(OUT); sp2 = collapse(sp2); items2 = drawings(OUT)
ck.log("new sheet: one page, V13 page box, flat", info2["rect"] == info["rect"] and info2["xobjects"] == 0, info2)
ck.log("fonts: Arial faces only", all("Arial" in f for f in info2["fonts"]), info2["fonts"])
static_ok = []
g = find_span(sp2, title["text"], y=title["y"], ytol=0.6); static_ok.append(g is not None and abs((g["bbox"][0] + g["bbox"][2]) / 2 - (AX["x0"] + AX["x1"]) / 2) < 0.6)   # round 49: centred on the axes box
exp_tick_labels = [("−" + str(-v)) if v < 0 else str(v) for v in TICK_VALS_NEW]
for lab, tx in zip(exp_tick_labels, TICK_X_NEW):
    g = find_span(sp2, lab, y=TICK_LAB_Y, ytol=0.6); static_ok.append(g is not None and abs((g["bbox"][0] + g["bbox"][2]) / 2 - tx) < 0.6)
ck.log(f"static strings: the axis title centred on the axes box (as V13); the {len(exp_tick_labels)} tick labels {exp_tick_labels} centred on the tick positions within 0.6 pt" + (" (axis extended and rescaled, FIX A1)" if AXIS_CHANGED else " (the V13 positions)"), all(static_ok), f"{sum(static_ok)} of {len(static_ok)}")
labs2 = column(sp2, 0, 20, size=10.0 + PT_PLUS)
exp_labels = [f"{g['group']} ({g['responding']} of {g['conditions']})" for g in grp]
ck.log("9 row labels read back in the recomputed order at the V13 origins", [s["text"] for s in labs2] == exp_labels and all(abs(s["y"] - y) < 0.3 and abs(s["x"] - LAB_X) < 0.3 for s, y in zip(labs2, LAB_Y)), [s["text"] for s in labs2])
R2 = sorted([r for r in rects(items2) if r["fill"] in (BLUE, GREY) and 8 <= r["h"] <= 20], key=lambda r: r["y0"])
ends = [(r["x1"] if g["median_pct_lower"] >= 0 else r["x0"]) for r, g in zip(R2, grp)]
a2, b2, res2 = fit_axis(vals, ends, logscale=False)
ck.log("9 bars read back: ends at the fitted positions of the medians (residual under 0.05 pt), the drawn zero x and scale, the V13 bar height, blue where positive and grey where not",
       len(R2) == 9 and float(np.abs(res2).max()) < 0.05 and abs(a2 - ZERO_X_NEW) < 0.05 and abs(b2 - b_new) < 0.001 and all(abs(r["h"] - BAR_H) < 0.02 for r in R2) and all((r["fill"] == BLUE) == (g["median_pct_lower"] >= 0) for r, g in zip(R2, grp)), f"b {b2:.4f} (V13 {b:.4f}) zero {a2:.2f} (V13 {ZERO_X:.2f}) max |res| {np.abs(res2).max():.3f}")
ck.log("no percentage label beside the bars (the V13 design of round 31)", not [s for s in sp2 if s["text"].endswith("%") and abs(s["size"] - 9.5) < 0.3], "")
W0, N0 = v14_multisets(v14_spans(BASE)); W1, N1 = v14_multisets(v14_spans(OUT))
old_tick_labels = [t["text"] for t in tick_labels.values()]
exp_removed = Counter(x for t in old_labels + old_tick_labels for x in re.findall(r"-?\d+\.?\d*", t)); exp_added = Counter(x for t in exp_labels + exp_tick_labels for x in re.findall(r"-?\d+\.?\d*", t))
ck.log("numeric-token multiset differs from V13 exactly by the row-label counts (old '(k of n)' out, new in) and the declared tick-label change", (N0 - N1) == (exp_removed - exp_added) and (N1 - N0) == (exp_added - exp_removed), f"removed {dict(N0 - N1)} added {dict(N1 - N0)}")
wr = (W0 - W1); wa = (W1 - W0)
ck.log("word multiset differs from V13 only by the row labels' words and the declared tick labels", all(w in " ".join(old_labels + old_tick_labels) for w in wr) and all(w in " ".join(exp_labels + exp_tick_labels) for w in wa), f"removed {dict(wr)} added {dict(wa)}")
colours = {hexcol(it["fill"]) for it in items2 if it["fill"]} | {hexcol(it["color"]) for it in items2 if it["color"]}
ck.log("colours within the V13 set", colours <= ({hexcol(it["fill"]) for it in items if it["fill"]} | {hexcol(it["color"]) for it in items if it["color"]}), f"{sorted(colours)}")
changes = []
# FIX A1: before = the UNADJUSTED V14 sheet (its readback json: groups with responding, conditions, median) snapshotted in work/
u14 = load_json(f"{R37LANE}/{SHEET}/work/{SHEET}_readback_UNADJUSTED_V14.json"); old_g = {g["group"]: g for g in u14["groups"]}; old_order = [g["group"] for g in u14["groups"]]   # round 49: the round-37 snapshot, read-only
for i, g in enumerate(grp):
    o = old_g[g["group"]]
    changes.append(dict(sheet=SHEET, panel="single", item=f"{g['group']}: row label (responding of conditions), unadjusted -> adjusted", before=f"{g['group']} ({o['responding']} of {o['conditions']})", after=f"{g['group']} ({g['responding']} of {g['conditions']})",
                        dramatic="count moves more than 20 percent" if (o["responding"] and abs(g["responding"] - o["responding"]) / o["responding"] > 0.2) or (o["responding"] == 0 and g["responding"] > 0) or (abs(g["conditions"] - o["conditions"]) / o["conditions"] > 0.2) else "", note=ADJ_NOTE + "; members: " + ", ".join(g["members"])))
    changes.append(dict(sheet=SHEET, panel="single", item=f"{g['group']}: median rate lower where restored, % (bar length), unadjusted -> adjusted", before=f"{o['median_pct_lower']:.1f}", after=f"{g['median_pct_lower']:.1f}",
                        dramatic=dramatic("pct", f"{o['median_pct_lower']:.1f}%", f"{g['median_pct_lower']:.1f}%"), note=""))
    changes.append(dict(sheet=SHEET, panel="single", item=f"{g['group']}: row position", before=str(old_order.index(g["group"]) + 1), after=str(i + 1), dramatic="rank changes" if old_order.index(g["group"]) != i else "", note="rows sorted by the adjusted median"))
if AXIS_CHANGED: changes.append(dict(sheet=SHEET, panel="single", item="axis range (ticks)", before=f"{TICK_VALS[0]} to {TICK_VALS[-1]}", after=f"{TICK_VALS_NEW[0]} to {TICK_VALS_NEW[-1]} (rescaled inside the V13 tick span)", dramatic="", note=f"the adjusted Liver median {min(vals):.1f} falls below the V13 axis; scale {b:.4f} -> {b_new:.4f} pt per percent"))
write_changes(f"{V}/CHANGES_{SHEET}.csv", changes)
n_rows, n_no = printed_values_csv(OUT, drawn, f"{V}/{SHEET}_printed_values.csv", static_from=BASE)
ck.log(f"printed values csv: {n_rows} rows, {n_no} NO", n_no == 0, f"{V}/{SHEET}_printed_values.csv")
render(OUT, 150, f"{W}/{SHEET}_built_150dpi.png"); render(BASE, 200, f"{W}/OLD_200dpi.png"); render(OUT, 200, f"{W}/NEW_200dpi.png"); render(BASE, 150, f"{W}/OLD_150dpi.png")
rd = raster_diff(f"{W}/OLD_150dpi.png", f"{W}/{SHEET}_built_150dpi.png", [(0, 12, PW, 270)], 150)
ck.log("raster diff at 150 dpi confined to the plot and axis region (a whole-sheet re-plot, Proof W: bars, labels, the re-set tick labels and title), 0 pixels in the margins", rd["shape_equal"] and rd["n_outside"] == 0, f"{rd}")
crop_pair(f"{W}/OLD_200dpi.png", f"{W}/NEW_200dpi.png", (0, 0, PW, PH), 200, f"{V}/crops/01_whole_sheet")
dump({"axis": {"a": a_new, "b": b_new, "zero_x": ZERO_X_NEW, "ticks": TICK_VALS_NEW, "v13_axis": {"a": a, "b": b, "zero_x": ZERO_X, "ticks": TICK_VALS}, "axis_changed": AXIS_CHANGED}, "groups": grp, "provenance": {"treatment_v3": sc_t, "gaps_v2": sc_g}, "adjusted": True, "note": ADJ_NOTE}, f"{V}/{SHEET}_readback.json")
RIGHT_LAB = max(s["bbox"][2] for s in labs2); ck.log(f"row labels at {10.0 + PT_PLUS:g} pt end at x {RIGHT_LAB:.2f}, clear of the axes box (x {AX['x0']:.2f}) by 4 pt", RIGHT_LAB + 4.0 <= AX["x0"], "")
dump({"pt_plus": PT_PLUS, "size_exceptions": {}, "nudges": {}, "not_drawn": [], "labels_right_x": RIGHT_LAB}, f"{W}/r49_record.json")
ok = ck.write(f"{W}/build_checks.txt"); print("RESULT ALL PASS" if ok else "RESULT FAIL")
