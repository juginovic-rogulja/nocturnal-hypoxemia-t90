#!/usr/bin/env python3
"""Round 49 LED-B copy (2026-09-26): every text +1 pt (l5a14.Sheet), the axis title centred on its V13 centre, the foot note not drawn (LNOTES round 40), outputs into work/. Otherwise the round-37 builder. ED Fig 8 (single panel, flat sheet, 49 outcomes still-low against restored), V14: re-plotted from scratch at the V13 geometry
from the v8.1 treatment_v3.json (step 211, the comparison ADJUSTED for prevalent cardiopulmonary disease, FIX A1, decision 6). Row rule (efig23_polish.py, round 30 read-back): all 49 corrected_vs_not outcomes sorted
by the hazard ratio descending, controls interleaved in grey; "% lower" printed where HR > 1; asterisk where P < 0.05 (not a
control). P VALUES ONLY, nothing added (owner's decision 7). V13 design (round 31 LFOREST): no pale grid lines, no grey connectors
between the reference circle and the marker; the hidden second copy of the key under a white rectangle (never printed) is not
reproduced (declared). Lineage: ROUND30 V7_L5a_PAP_MAIN/ED_Fig08/scripts/run_ed8.py (read-back only)."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from l5a14 import *

SHEET = "ED_Fig08"; D = f"{LANE}/{SHEET}"; W = f"{D}/work"; V = f"{D}/verify"; os.makedirs(f"{V}/crops", exist_ok=True)
BASE = f"{V13}/{SHEET}.pdf"; OUT = f"{W}/{SHEET}_built.pdf"   # round 49: the built sheet (V13 page size) goes to work/, the runner adds the title strip
ck = Checks(f"Lane V14_L5a_PAP_MAIN, {SHEET}, {time.strftime('%Y-%m-%d %H:%M')}. BASE (design) = V13 {BASE} (sha256 {sha256(BASE)[:16]}). NEW = treatment_v3.json v8.1 (step 211, ADJUSTED for prevalent cardiopulmonary disease, FIX A1).")

# ---- 0 gate
sc_t = gate("treatment_v3"); dump({"treatment_v3": sc_t}, f"{W}/sidecars.json")
ck.log("treatment_v3.json passes the v8.1 sidecar gate (step 211, inputs data_frozen_v8_2026-09, sha equal); the ADJUSTED fields adj_hr, adj_lo, adj_hi, adj_p, adj_pct_lower_if_corrected are drawn (FIX A1)", sc_t["v8_inputs"], f"out {sc_t['output_mtime']} sha {sc_t['sha256'][:16]}")

# ---- 1 the V13 geometry
spans, info = text_spans(BASE); spans = collapse(spans); items = drawings(BASE)
dump({"info": info, "spans": spans}, f"{W}/base_text.json"); dump(items, f"{W}/base_drawings.json")
PW, PH = info["rect"][2], info["rect"][3]
ck.log("V13 sheet flat (0 XObjects), page 484.72 x 656.755", info["xobjects"] == 0 and abs(PW - 484.72) < 0.05 and abs(PH - 656.755) < 0.05, info)
labels = column(spans, 0, 140, ymin=30, ymax=590, size=10.0); assert len(labels) == 49, len(labels)
LAB_R = float(np.mean([s["bbox"][2] for s in labels])); LAB_Y = [s["y"] for s in labels]; pitch = float(np.mean(np.diff(LAB_Y)))
assert np.std(np.diff(LAB_Y)) < 0.06 and abs(np.std([s["bbox"][2] for s in labels])) < 0.1, (np.diff(LAB_Y), LAB_R)
M = markers(items)
dat = sorted([m for m in M if m["shape"] == "circle" and not m["open"] and m["colour"] in (BLUE, CTRL_GREY) and abs(m["w"] - 5.29) < 0.05 and m["cx"] < 320], key=lambda m: m["cy"]); assert len(dat) == 49
ref = sorted([m for m in M if m["shape"] == "circle" and m["open"] and abs(m["w"] - 5.1) < 0.05 and m["cx"] < 320], key=lambda m: m["cy"]); assert len(ref) == 49
MARK_D = float(np.mean([m["w"] for m in dat])); MARK_EDGE = float(np.mean([m["widths"][0] for m in dat])); REF_D = float(np.mean([m["w"] for m in ref])); REF_EDGE = float(np.mean([m["widths"][0] for m in ref]))
OFF_M = float(np.mean([m["cy"] - s["y"] for m, s in zip(dat, labels)])); assert np.std([m["cy"] - s["y"] for m, s in zip(dat, labels)]) < 0.05
X_REF = float(np.median([m["cx"] for m in ref]))
ci = sorted([h for h in hlines(items, min_len=0.5) if h["stroke"] in (BLUE, CTRL_GREY) and h["width"] and h["width"] > 1.5 and h["x0"] < 320 and 30 <= h["y"] <= 590], key=lambda h: h["y"]); assert len(ci) == 49, len(ci)
CI_W = float(np.mean([h["width"] for h in ci]))
spine = [h for h in hlines(items, min_len=100) if h["stroke"] == INK]; assert len(spine) == 1; SPINE = spine[0]; SPINE_Y = SPINE["y"]
refline = [it for it in items if "line" in it and abs(it["line"][0] - it["line"][2]) < 0.05 and abs(it["line"][3] - it["line"][1]) > 100 and hexcol(it["color"]) == INK]; assert len(refline) == 1; REFL = refline[0]
ticks = sorted([v for v in vlines(items, max_len=6.0) if v["stroke"] == INK and abs(v["y0"] - SPINE_Y) < 0.05], key=lambda v: v["x"]); assert len(ticks) == 5
TICK_VALS = [0.25, 0.5, 1, 2, 4]; TICK_LEN = float(np.mean([v["y1"] - v["y0"] for v in ticks]))
a, b, res = fit_axis(TICK_VALS, [v["x"] for v in ticks], logscale=True); assert float(np.abs(res).max()) < 0.05, res
tick_labels = sorted([s for s in spans if abs(s["y"] - 601.29) < 0.6], key=lambda s: s["x"]); assert [s["text"] for s in tick_labels] == ["0.25", "0.5", "1", "2", "4"]
TICK_LAB_Y = float(np.mean([s["y"] for s in tick_labels]))
pcts = column(spans, 325, 340, ymin=30, ymax=590, size=9.5); stars = column(spans, 310, 320, ymin=30, ymax=590, size=11.0)
PCT_X = float(np.mean([p["x"] for p in pcts])); PCT_OFF = float(np.mean([p["y"] - min(labels, key=lambda s: abs(s["y"] - p["y"]))["y"] for p in pcts]))
STAR_X = float(np.mean([p["x"] for p in stars])); STAR_OFF = float(np.mean([p["y"] - min(labels, key=lambda s: abs(s["y"] - p["y"]))["y"] for p in stars]))
head_pct = find_span(spans, "% lower", ymin=20, ymax=30); title = find_span(spans, "Rate of new diagnosis, oxygen still low versus oxygen restored (95% CI)"); foot = find_span(spans, "* significant fall (P < 0.05)")
key_txt = [s for s in spans if s["x"] > 380 and 520 <= s["y"] <= 600]; assert len(key_txt) == 5, key_txt
hidden_key = [s for s in spans if s["x"] > 380 and 40 <= s["y"] <= 110]; assert len(hidden_key) == 5, hidden_key
ck.log(f"V13 geometry read: 49 labels right-aligned at x {LAB_R:.2f}, pitch {pitch:.3f}; markers {MARK_D:.2f} pt (edge {MARK_EDGE:.2f}) at label y {OFF_M:+.3f}; reference circles {REF_D:.2f} pt at x {X_REF:.3f}; CI {CI_W:.2f} pt; log axis {b:.3f} pt per ln unit; ticks 0.25 to 4; '% lower' at x {PCT_X:.2f} ({PCT_OFF:+.2f}), star at x {STAR_X:.2f} ({STAR_OFF:+.2f})",
       abs(a - X_REF) < 0.05 and head_pct is not None and title is not None and foot is not None, f"a {a:.3f}")
old_labels = [s["text"] for s in labels]; old_pcts = [p["text"] for p in pcts]; old_n_stars = len(stars)

# ---- 2 the v8.1 rows
T, rows = treatment_rows_adj(); cvn = T["corrected_vs_not"]
ck.log(f"treatment_v3.json v8.1 (adjusted view): {len(rows)} outcomes, n {cvn['n']:,} = {cvn['n_corrected']:,} restored + {cvn['n_not']:,} still above; controls in the file: {[k for k, v in rows if v['negative_control']]}",
       len(rows) == 49 and cvn["n"] == cvn["n_corrected"] + cvn["n_not"], "")
SIG = [k for k, v in rows if v["p"] < 0.05 and not v["negative_control"]]; N_GT1 = sum(1 for _, v in rows if v["hr"] > 1)
ck.log(f"recount on the v8.1 file: {len(SIG)} significant falls (P < 0.05, not a control), {N_GT1} rows with HR above 1 print '% lower'; controls significant: {[k for k, v in rows if v['negative_control'] and v['p'] < 0.05]}", not any(v["negative_control"] and v["p"] < 0.05 for _, v in rows), f"SIG {SIG}")
ck.log("adj_pct_lower_if_corrected equals 100 (1 - 1/adj_HR) within 0.2 on all rows", all(abs(v["pct_lower_if_corrected"] - 100 * (1 - 1 / v["hr"])) < 0.2 for _, v in rows), "")
xlo, xhi = exp((SPINE["x0"] - a) / b), exp((SPINE["x1"] - a) / b)
ck.log("every hazard ratio and interval end lies inside the V13 axis range (no tick change)", all(xlo <= v["lo"] and v["hi"] <= xhi for _, v in rows), f"axis {xlo:.3f} to {xhi:.3f}; data {min(v['lo'] for _, v in rows):.3f} to {max(v['hi'] for _, v in rows):.3f}")

# ---- 3 draw
S = Sheet(PW, PH)
S.line(REFL["line"][0], REFL["line"][1], REFL["line"][2], REFL["line"][3], INK, REFL["width"], cap=0, dashes=REFL.get("dashes"))
drawn_rows = []
for i, (k, v) in enumerate(rows):
    y = LAB_Y[0] + i * pitch; cy = y + OFF_M; col = CTRL_GREY if v["negative_control"] else BLUE
    S.line(axis_x(a, b, v["lo"]), cy, axis_x(a, b, v["hi"]), cy, col, CI_W, cap=1)
    S.circle(X_REF, cy, REF_D, WHITE, GREY, REF_EDGE)
    S.circle(axis_x(a, b, v["hr"]), cy, MARK_D, col, WHITE, MARK_EDGE)
    S.text(LAB_R, y, plab(k), 10.0, color=GREY if v["negative_control"] else INK, align="right", panel="single", source_file=SPEC["treatment_v3"], source_key=f"outcomes/{k} (row label, rows sorted by adj_hr descending)", source_value=k, rule="label", note="" if plab(k) == k else f"printed as the short form used on Fig 2 (file label: {k})")
    if v["hr"] > 1:
        S.text(PCT_X, y + PCT_OFF, pct_text(v["pct_lower_if_corrected"]), 9.5, color=GREY, panel="single", source_file=SPEC["treatment_v3"], source_key=f"outcomes/{k}/adj_pct_lower_if_corrected", source_value=v["pct_lower_if_corrected"], rule="pct_text")
    if v["p"] < 0.05 and not v["negative_control"]:
        S.text(STAR_X, y + STAR_OFF, "*", 11.0, panel="single", source_file=SPEC["treatment_v3"], source_key=f"outcomes/{k}/adj_p", source_value=v["p"], rule="star_p05", note=f"adjusted P {v['p']:.3g} < 0.05")
    drawn_rows.append(dict(row=k, hr=v["hr"], lo=v["lo"], hi=v["hi"], p=v["p"], pct=v["pct_lower_if_corrected"], control=v["negative_control"], x=axis_x(a, b, v["hr"]), y=cy, x_lo=axis_x(a, b, v["lo"]), x_hi=axis_x(a, b, v["hi"])))
S.line(SPINE["x0"], SPINE_Y, SPINE["x1"], SPINE_Y, INK, SPINE["width"], cap=2)
for t, val, lab in zip(ticks, TICK_VALS, tick_labels):
    S.line(t["x"], SPINE_Y, t["x"], SPINE_Y + TICK_LEN, INK, t["width"], cap=0)
    S.text(t["x"], TICK_LAB_Y, lab["text"], 10.0, align="center", panel="single", source_key="axis tick")
n_key = replay_items(S, items, 360, 520, 480, 600)
for s in key_txt + [head_pct]: S.text(s["x"], s["y"], s["text"], s["size"], color="#" + s["color"], panel="single", source_key="static V13 string")
TITLE_CX = (title["bbox"][0] + title["bbox"][2]) / 2   # round 49: the axis title is centred on its V13 centre (the page centre); the foot note is NOT drawn (removed by lane LNOTES, round 40)
S.text(TITLE_CX, title["y"], title["text"], title["size"], align="center", color="#" + title["color"], panel="single", source_key="static V13 string (centred)")
S.save(OUT); drawn = S.drawn; dump(drawn, f"{W}/{SHEET}_drawn.json"); dump(drawn_rows, f"{W}/{SHEET}_rows.json")

# ---- 4 verify
sp2, info2 = text_spans(OUT); sp2 = collapse(sp2); items2 = drawings(OUT)
ck.log("new sheet: one page, V13 page box, flat", info2["rect"] == info["rect"] and info2["xobjects"] == 0, info2)
ck.log("fonts: Arial faces only", all("Arial" in f for f in info2["fonts"]), info2["fonts"])
statics = key_txt + [head_pct]
ok_static = [find_span(sp2, s["text"], y=s["y"], ytol=0.6) is not None and abs(find_span(sp2, s["text"], y=s["y"], ytol=0.6)["x"] - s["x"]) < 0.6 for s in statics]
g = find_span(sp2, title["text"], y=title["y"], ytol=0.6); ok_static.append(g is not None and abs((g["bbox"][0] + g["bbox"][2]) / 2 - TITLE_CX) < 0.6)
for t, lab in zip(ticks, tick_labels):
    g = find_span(sp2, lab["text"], y=TICK_LAB_Y, ytol=0.6); ok_static.append(g is not None and abs((g["bbox"][0] + g["bbox"][2]) / 2 - t["x"]) < 0.6)
ck.log(f"static strings ({len(ok_static)}: key and column head left-aligned at the V13 x, the axis title centred on its V13 centre, tick labels centred on the ticks, within 0.6 pt) at +{PT_PLUS:g} pt; the foot note '* significant fall (P < 0.05)' not drawn (LNOTES round 40)", all(ok_static) and find_span(sp2, foot["text"]) is None, f"{sum(ok_static)} of {len(ok_static)}")
ck.log("the hidden second copy of the key (V13 ink text under a white rectangle, never printed) is not reproduced", not [s for s in sp2 if s["x"] > 380 and 40 <= s["y"] <= 110], "")
labs2 = column(sp2, 0, 140, ymin=30, ymax=590, size=10.0 + PT_PLUS)
ck.log("49 row labels read back in the file's order (hr descending), right-aligned at the V13 edge, on the V13 baselines, grey on the controls", [s["text"] for s in labs2] == [plab(k) for k, _ in rows] and all(abs(s["bbox"][2] - LAB_R) < 0.6 and abs(s["y"] - y) < 0.3 for s, y in zip(labs2, LAB_Y)) and all((s["color"] == GREY.lstrip("#")) == v["negative_control"] for s, (k, v) in zip(labs2, rows)), f"{len(labs2)}")
pc2 = column(sp2, 325, 340, ymin=30, ymax=590, size=9.5 + PT_PLUS); st2 = column(sp2, 310, 320, ymin=30, ymax=590, size=11.0 + PT_PLUS)
exp_pct = [(pct_text(v["pct_lower_if_corrected"]), LAB_Y[i] + PCT_OFF) for i, (k, v) in enumerate(rows) if v["hr"] > 1]
ck.log(f"'% lower' printed on exactly the {N_GT1} rows with HR above 1, equal to the file to one decimal, at the V13 offsets", len(pc2) == len(exp_pct) and all(p["text"] == t and abs(p["y"] - y) < 0.3 for p, (t, y) in zip(pc2, exp_pct)), f"{len(pc2)} strings")
ck.log(f"asterisk on exactly the {len(SIG)} rows with P < 0.05 (no control)", len(st2) == len(SIG) and all(abs(s["y"] - (LAB_Y[[k for k, _ in rows].index(kk)] + STAR_OFF)) < 0.3 for s, kk in zip(st2, SIG)), f"{len(st2)}")
# FIX A1: a data marker whose adjusted HR is close to 1 sits on its reference circle, so the two are clustered from SEPARATE item subsets (filled blue/grey data circles vs white-filled grey-edged reference circles) instead of one mixed pass
dat2 = sorted([m for m in markers([it for it in items2 if hexcol(it["fill"]) in (BLUE, CTRL_GREY)]) if m["shape"] == "circle" and abs(m["w"] - MARK_D) < 0.06 and m["cx"] < 320], key=lambda m: m["cy"])
ref2 = sorted([m for m in markers([it for it in items2 if hexcol(it["fill"]) == WHITE and hexcol(it["color"]) == GREY]) if m["shape"] == "circle" and abs(m["w"] - REF_D) < 0.06 and m["cx"] < 320], key=lambda m: m["cy"])
ci2 = sorted([h for h in hlines(items2, min_len=0.5) if h["stroke"] in (BLUE, CTRL_GREY) and h["width"] and h["width"] > 1.5 and h["x0"] < 320 and 30 <= h["y"] <= 590], key=lambda h: h["y"])
a2, b2, res2 = fit_axis([v["hr"] for _, v in rows], [m["cx"] for m in dat2], logscale=True) if len(dat2) == 49 else (0, 1, np.array([9]))
ci_dev = max(max(abs(h["x0"] - axis_x(a, b, v["lo"])), abs(h["x1"] - axis_x(a, b, v["hi"]))) for h, (k, v) in zip(ci2, rows)) if len(ci2) == 49 else 9
ck.log("49 data markers, 49 reference circles and 49 interval lines read back from the new drawings; marker x against ln(HR) of the file on the V13 axis (residual under 0.05 pt), interval ends within 0.05 pt of ln(lo), ln(hi), grey on exactly the controls",
       len(dat2) == 49 and len(ref2) == 49 and len(ci2) == 49 and float(np.abs(res2).max()) < 0.05 and abs(a2 - a) < 0.05 and abs(b2 - b) < 0.01 and ci_dev < 0.05 and all((m["colour"] == CTRL_GREY) == v["negative_control"] for m, (k, v) in zip(dat2, rows)) and all(abs(m["cx"] - X_REF) < 0.05 for m in ref2),
       f"markers {len(dat2)} ref {len(ref2)} ci {len(ci2)} max |res| {np.abs(res2).max():.3f} ci dev {ci_dev:.3f}")
ck.log("no pale #eef0f1 vertical grid line and no 1.0 pt grey connector on the sheet (the V13 design of round 31)", not [it for it in items2 if hexcol(it["color"]) == "#eef0f1"] and not [h for h in hlines(items2, min_len=0.0) if h["stroke"] == CTRL_GREY and h["width"] and abs(h["width"] - 1.0) < 0.05], "")
W0, N0 = v14_multisets(v14_spans(BASE)); W1, N1 = v14_multisets(v14_spans(OUT))
tok = lambda strs: Counter(x for t in strs for x in re.findall(r"-?\d+\.?\d*", t))
exp_removed = tok(old_pcts) + tok(old_labels) + tok([s["text"] for s in hidden_key]) + tok([foot["text"]]); exp_added = tok([t for t, _ in exp_pct]) + tok([plab(k) for k, _ in rows])
ck.log("numeric-token multiset differs from V13 exactly by the '% lower' strings, the row labels and the hidden key copy (out) against the new '% lower' strings and labels (in)", (N0 - N1) == (exp_removed - exp_added) and (N1 - N0) == (exp_added - exp_removed), f"removed {dict(N0 - N1)} added {dict(N1 - N0)}")
wtok = lambda strs: Counter(w for t in strs for w in t.split())
exp_wr = wtok(old_labels) + wtok(old_pcts) + wtok(["*"] * old_n_stars) + wtok([s["text"] for s in hidden_key]) + wtok([foot["text"]]); exp_wa = wtok([plab(k) for k, _ in rows]) + wtok([t for t, _ in exp_pct]) + wtok(["*"] * len(SIG))
ck.log("word multiset differs from V13 exactly by the row labels, '% lower' strings, asterisks and the hidden key copy", (W0 - W1) == (exp_wr - exp_wa) and (W1 - W0) == (exp_wa - exp_wr), f"removed {dict(W0 - W1)} added {dict(W1 - W0)}")
cols = {hexcol(it["fill"]) for it in items2 if it["fill"]} | {hexcol(it["color"]) for it in items2 if it["color"]}
ck.log("colours within the V13 set", cols <= ({hexcol(it["fill"]) for it in items if it["fill"]} | {hexcol(it["color"]) for it in items if it["color"]}), f"{sorted(cols)}")
# CHANGES (FIX A1): before = the UNADJUSTED V14 sheet snapshotted in work/ (its readback json and text layer), after = the adjusted sheet
u14 = load_json(f"{R37LANE}/{SHEET}/work/{SHEET}_readback_UNADJUSTED_V14.json"); old_rows = {r["row"]: r for r in u14["rows"]}   # round 49: the round-37 snapshot, read-only
usp, _ = unadjusted_v14_spans(SHEET); u_labels = [s["text"] for s in column(usp, 0, 140, ymin=30, ymax=590, size=10.0)]
u_order = [r["row"] for r in sorted(u14["rows"], key=lambda r: r["y"])]
changes = []; new_keys = [k for k, _ in rows]
for k in u_order:
    if k not in new_keys: changes.append(dict(sheet=SHEET, panel="single", item=f"{k}: row", before="drawn", after="absent", dramatic="row enters or leaves", note="not in treatment_v3.json"))
for i, (k, v) in enumerate(rows):
    o = old_rows.get(k)
    if o is None: changes.append(dict(sheet=SHEET, panel="single", item=f"{k}: row", before="absent", after=f"drawn at position {i + 1} as '{plab(k)}'", dramatic="row enters or leaves", note="")); continue
    oh = f"{o['hr']:.3f} ({o['lo']:.3f}-{o['hi']:.3f})"; nh = f"{v['hr']:.3f} ({v['lo']:.3f}-{v['hi']:.3f})"
    changes.append(dict(sheet=SHEET, panel="single", item=f"{k}: hazard ratio (marker, no numeral), unadjusted -> adjusted", before=oh, after=nh, dramatic=dramatic("hr", oh, nh), note=ADJ_NOTE))
    changes.append(dict(sheet=SHEET, panel="single", item=f"{k}: % lower", before=pct_text(o["pct"]) if o["hr"] > 1 else "", after=pct_text(v["pct_lower_if_corrected"]) if v["hr"] > 1 else "", dramatic=dramatic("pct", f"{o['pct']}%", f"{v['pct_lower_if_corrected']}%") if (o["hr"] > 1 and v["hr"] > 1) else ("sign flips" if (o["hr"] > 1) != (v["hr"] > 1) else ""), note=""))
    changes.append(dict(sheet=SHEET, panel="single", item=f"{k}: asterisk (P < 0.05)", before="*" if (o["p"] < 0.05 and not o["control"]) else "", after="*" if (v["p"] < 0.05 and not v["negative_control"]) else "", dramatic="q or P crosses 0.05" if (o["p"] < 0.05) != (v["p"] < 0.05) else "", note=f"P {o['p']:.3g} -> {v['p']:.3g} (adjusted)"))
    oi = u_order.index(k)
    if oi != i: changes.append(dict(sheet=SHEET, panel="single", item=f"{k}: row position", before=str(oi + 1), after=str(i + 1), dramatic="", note="rows sorted by the adjusted hazard ratio"))
write_changes(f"{V}/CHANGES_{SHEET}.csv", changes)
n_rows, n_no = printed_values_csv(OUT, drawn, f"{V}/{SHEET}_printed_values.csv", static_from=BASE)
ck.log(f"printed values csv: {n_rows} rows, {n_no} NO", n_no == 0, f"{V}/{SHEET}_printed_values.csv")
render(OUT, 150, f"{W}/{SHEET}_built_150dpi.png"); render(BASE, 200, f"{W}/OLD_200dpi.png"); render(OUT, 200, f"{W}/NEW_200dpi.png"); render(BASE, 150, f"{W}/OLD_150dpi.png")
rd = raster_diff(f"{W}/OLD_150dpi.png", f"{W}/{SHEET}_built_150dpi.png", [(0, 15, PW, 645)], 150)
ck.log("raster diff at 150 dpi confined to the content area (a whole-sheet re-plot, Proof W: the static strings are re-set in the full Arial face, so their pixels differ by anti-aliasing only), 0 pixels in the margins", rd["shape_equal"] and rd["n_outside"] == 0, f"{rd}")
crop_pair(f"{W}/OLD_200dpi.png", f"{W}/NEW_200dpi.png", (0, 20, 360, 250), 200, f"{V}/crops/01_top_rows"); crop_pair(f"{W}/OLD_200dpi.png", f"{W}/NEW_200dpi.png", (0, 250, 360, 610), 200, f"{V}/crops/02_lower_rows_and_axis"); crop_pair(f"{W}/OLD_200dpi.png", f"{W}/NEW_200dpi.png", (360, 30, 484.72, 610), 200, f"{V}/crops/03_key_column")
dump({"axis": {"a": a, "b": b, "x_ref": X_REF}, "rows": drawn_rows, "n_sig": len(SIG), "n_hr_gt1": N_GT1, "provenance": {"treatment_v3": sc_t}, "adjusted": True, "note": ADJ_NOTE}, f"{V}/{SHEET}_readback.json")
LEFT_INK = min(s["bbox"][0] for s in labs2); ck.log(f"row labels right-aligned at x {LAB_R:.2f}: the widest label at {10.0 + PT_PLUS:g} pt starts at x {LEFT_INK:.2f} (inside the page, no clipping; V13 gutter margin was {min(s['bbox'][0] for s in labels):.2f})", LEFT_INK > 0.0, "")
dump({"pt_plus": PT_PLUS, "size_exceptions": {}, "nudges": {}, "not_drawn": [foot["text"]], "title_centre_x": TITLE_CX, "left_ink_x": LEFT_INK}, f"{W}/r49_record.json")
ok = ck.write(f"{W}/build_checks.txt"); print("RESULT ALL PASS" if ok else "RESULT FAIL")
