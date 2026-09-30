#!$T90_PY
"""ED_Fig06, round 51 (lane LED-B, 2026-09-26), Alen's item 3 as the coordinator specified: panel b without the "New diagnoses / group size" column,
the freed width given to the forest axis (the HR (95% CI) column and the key stay at their V28 x). Rebuilt on a fresh page from the V28 sheet
(flat, 0 XObjects; read with PyMuPDF in this child under wd_run.sh):
 - every V28 drawing item outside panel b's plot region is replayed with fitz.Shape in its original order with its own fill, stroke, width,
   dashes, caps, joins, even-odd rule and opacity (the clipped stale full-page background of the old Chromium layer is skipped, invisible on V28);
 - panel b's forest is RE-PLOTTED on a wider log axis: the plot box keeps its right edge (356.8) and its y range and starts at x 209.0 instead of
   253.44 (the freed width), the axis range (the values at the box edges) is unchanged, so the ticks 0.5, 1, 2, 4, the dashed HR = 1 rule, the 24
   intervals and the 24 markers land at x' = a' + b' ln(v) with the V28 mark geometry (interval 1.7 pt round caps, grey square / blue circle with the
   0.8 pt white edge) from the values of results_crossed.csv (sidecar-gated, the same values the printed HR (95% CI) strings carry, checked against
   the V28 marks on the V28 axis); the overlay residue of V28 (old marks under white bands, the bands, the band segments of the rule) is not
   reproduced (invisible);
 - every V28 string except the 26 of the removed column is re-written with fitz.TextWriter at its V28 position, size, weight and colour, the four
   tick labels centred on the new tick x.
Text sizes and every printed number unchanged. usage: 01_r51_close_up_ed6.py <in V28.pdf> <work dir> <out.pdf>"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import fitz, json, os, sys, gc, re, hashlib, math
import pandas as pd
IN, WORK, OUT = sys.argv[1:4]; os.makedirs(WORK, exist_ok=True)
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND37_2026-09-14/figures/lanes/V14_L6_HABITUAL_EXTERNAL/_common"); import lane_common as LC
ARIAL = paths.FONT_ARIAL; ARIALB = paths.FONT_ARIAL_BOLD
BLUE, GREY, INK, WHITE = "#0288d1", "#8a9099", "#1a1d21", "#ffffff"
PLOT = fitz.Rect(246.0, 235.0, 360.0, 636.0); X0_NEW = 209.0
assert os.stat(IN).st_blocks > 0
FR, FB = fitz.Font(fontfile=ARIAL), fitz.Font(fontfile=ARIALB)
def hx(c): return None if c is None else "#%02x%02x%02x" % tuple(int(round(v * 255)) for v in c[:3])
src = fitz.open(IN); p = src[0]; W, H = p.rect.width, p.rect.height; assert len(p.get_xobjects()) == 0 and abs(W - 518.74) < 0.01 and abs(H - 666.5) < 0.01
spans = []
for b in p.get_text("rawdict")["blocks"]:
    if b["type"] != 0: continue
    for l in b["lines"]:
        for s in l["spans"]:
            txt = "".join(ch["c"] for ch in s["chars"]).replace("\xa0", " ").replace("\xad", "-")
            if not txt.strip(): continue
            spans.append(dict(text=txt, font=s["font"], size=round(s["size"], 3), color="#%06x" % s["color"], origin=[round(v, 3) for v in s["chars"][0]["origin"]], bbox=[round(v, 3) for v in s["bbox"]], dir=[round(v, 3) for v in l["dir"]]))
assert len(spans) == 116, len(spans)
EV_RE = re.compile(r"^\d+/[\d,]+$")
removed = [s for s in spans if (EV_RE.match(s["text"]) and s["origin"][1] > 235) or s["text"] in ("New diagnoses /", "group size")]; assert len(removed) == 26
ticklab = [s for s in spans if s["text"] in ("0.5", "1", "2", "4") and abs(s["origin"][1] - 626.876) < 0.5 and s["size"] == 11.0]; assert len(ticklab) == 4
kept = [s for s in spans if s not in removed and s not in ticklab]; assert len(kept) == 86
raw = p.get_drawings(); assert len(raw) == 284, len(raw)
# the clipped stale layer of the old Chromium export (invisible on V28): items that repeat a visible item 3.629 pt lower with the same x range and style
def same_style(a, b): return a["fill"] == b["fill"] and a["color"] == b["color"] and (a.get("width") or 0) == (b.get("width") or 0) and "".join(sorted(set(x[0] for x in a["items"]))) == "".join(sorted(set(x[0] for x in b["items"])))
STALE_BG = [k for k, it in enumerate(raw) if any(abs(it["rect"].x0 - o["rect"].x0) < 0.01 and abs(it["rect"].x1 - o["rect"].x1) < 0.01 and abs((it["rect"].y0 - o["rect"].y0) - 3.629) < 0.01 and same_style(it, o) for o in raw)]
assert len(STALE_BG) == 3 and any(raw[k]["rect"].width > 500 for k in STALE_BG), STALE_BG   # the page background, the reference cell's blue fill, its dashed white border
inplot = [k for k, it in enumerate(raw) if PLOT.contains(it["rect"])]
strad = [it["rect"] for it in raw if (not PLOT.contains(it["rect"])) and it["rect"].intersects(PLOT)]; assert all(r.width > 500 for r in strad), strad
# ---- the V28 forest: axis from the four tick rects, the visible marks (the last-drawn 24 intervals and 24 markers), the rule, the box
P = [raw[k] for k in inplot]
ticks = sorted([it for it in P if it["color"] is not None and hx(it["color"]) == INK and 2.5 < it["rect"].height < 3.5 and it["rect"].width < 1.0 and it["items"][0][0] == "l"], key=lambda it: it["rect"].x0); assert len(ticks) == 4, len(ticks)
spine = [it for it in P if it["items"][0][0] == "l" and len(it["items"]) == 1 and abs(it["items"][0][1].y - it["items"][0][2].y) < 0.01 and it["rect"].width > 50 and it["color"] is not None and hx(it["color"]) == INK and abs((it["width"] or 0) - 1.7) > 0.1]
assert len(spine) == 1, len(spine); SPINE = spine[0]; SPINE_Y = SPINE["items"][0][1].y; assert abs(SPINE_Y - 613.22) < 0.05 and abs(SPINE["rect"].x0 - 253.44) < 0.05 and abs(SPINE["rect"].x1 - 356.8) < 0.05
TX = [(it["rect"].x0 + it["rect"].x1) / 2 for it in ticks]; TV = [0.5, 1.0, 2.0, 4.0]
import numpy as np
b_, a_ = np.polyfit([math.log(v) for v in TV], TX, 1); assert max(abs(a_ + b_ * math.log(v) - x) for v, x in zip(TV, TX)) < 0.05
box = [it for it in P if it["fill"] is not None and hx(it["fill"]) == WHITE and it["color"] is None and it["rect"].height > 300]; assert len(box) == 1; BOX = box[0]["rect"]
rule = [it for it in P if it.get("dashes") not in (None, "[] 0") and it["color"] is not None and hx(it["color"]) == INK and it["rect"].height > 300]; assert len(rule) == 1; RULE = rule[0]
assert abs(RULE["rect"].x0 - (a_ + b_ * 0.0)) < 0.05
LO, HI = math.exp((BOX.x0 - a_) / b_), math.exp((BOX.x1 - a_) / b_)
b2 = (BOX.x1 - X0_NEW) / (math.log(HI) - math.log(LO)); a2 = X0_NEW - b2 * math.log(LO)
def xo(v): return a_ + b_ * math.log(v)
def xn(v): return a2 + b2 * math.log(v)
ci_all = [it for it in P if it["color"] is not None and hx(it["color"]) in (BLUE, GREY) and abs((it["width"] or 0) - 1.7) < 0.01 and it["rect"].width > 5]
mk_all = [it for it in P if it["fill"] is not None and hx(it["fill"]) in (BLUE, GREY) and it["rect"].width < 8 and it["rect"].height < 8]
assert len(ci_all) >= 48 and len(mk_all) >= 48, (len(ci_all), len(mk_all))
ci_vis, mk_vis = ci_all[-24:], mk_all[-24:]     # the topmost (last-drawn) marks are the visible ones (the round-37 overlay rule)
rows28 = []
for m in sorted(mk_vis, key=lambda it: it["rect"].y0):
    r = m["rect"]; cy = (r.y0 + r.y1) / 2; col = hx(m["fill"]); c = min((c for c in ci_vis if hx(c["color"]) == col), key=lambda c: abs(c["rect"].y0 - cy)); assert abs(c["rect"].y0 - cy) < 0.1
    rows28.append(dict(y=cy, col=col, kind="square" if m["items"][0][0] == "re" else "circle", size=r.width, edge_w=m["width"], edge=hx(m["color"]), hr=math.exp(((r.x0 + r.x1) / 2 - a_) / b_), lo=math.exp((c["rect"].x0 - a_) / b_), hi=math.exp((c["rect"].x1 - a_) / b_), ci_item=c, mk_item=m))
assert [r["col"] for r in rows28] == [GREY, BLUE] * 12
# ---- the values (the same files as the printed strings; round-37 gate)
SRC = f"{LC.SV}/habitual_outcomes_v2"; sc = LC.sidecar(LC.hydrated(f"{SRC}/results_crossed.csv"), min_mtime="2026-09-15 19:18:00")
crossed = pd.read_csv(f"{SRC}/results_crossed.csv"); cb = crossed[(crossed.model == "cells12") & (crossed["set"] == "primary")]
OUTC = [("cvd", "Cardiovascular composite"), ("death", "Death from any cause"), ("htn2", "Hypertension"), ("diabetes", "Type 2 diabetes"), ("hf", "Heart failure"), ("dementia", "Dementia")]
CELLS4 = [("<5h", "normal T90<=1%"), ("<5h", "low T90>10%"), (">=7h", "normal T90<=1%"), (">=7h", "low T90>10%")]
vals = []
for k, name in OUTC:
    dd = cb[cb.outcome_key == k].set_index("contrast")
    for hb, ob in CELLS4:
        r = dd.loc[f"{hb} | {ob}"]; vals.append(dict(outcome=name, cell=f"{hb} | {ob}", hr=float(r.hr), lo=float(r.lo95), hi=float(r.hi95), oxygen="normal" if ob.startswith("normal") else "low"))
assert len(vals) == 24
worst = 0.0
for r28, v in zip(rows28, vals):
    assert (r28["col"] == GREY) == (v["oxygen"] == "normal")
    worst = max(worst, abs(xo(v["hr"]) - xo(r28["hr"])), abs(xo(v["lo"]) - xo(r28["lo"])), abs(xo(v["hi"]) - xo(r28["hi"])))
assert worst < 0.05, ("the V28 marks do not sit at the file's values on the V28 axis", worst)
# ---- the fresh page
doc = fitz.open(); pg = doc.new_page(width=W, height=H); sh = pg.new_shape(); JOINS = []
def fin(it, **over):
    lc = it.get("lineCap"); lc = lc[0] if isinstance(lc, (tuple, list)) else (lc or 0); lj = int(it.get("lineJoin") or 0)
    closed = bool(it.get("closePath")) or all(x[0] in ("re", "qu") for x in it["items"])     # a rectangle or quad is a closed figure (V28 strokes the reference cell's quad closed)
    kw = dict(width=it.get("width") if it.get("width") is not None else 0, color=it.get("color"), fill=it.get("fill"), dashes=it.get("dashes"), lineCap=lc, lineJoin=lj, closePath=closed, even_odd=bool(it.get("even_odd")),
              fill_opacity=1 if it.get("fill_opacity") is None else it["fill_opacity"], stroke_opacity=1 if it.get("stroke_opacity") is None else it["stroke_opacity"]); kw.update(over)
    if kw["lineJoin"]: JOINS.append(kw["lineJoin"])       # PyMuPDF 1.27.2 writes the literal '{lineJoin} j' for a non-zero join: patched in the stream below
    sh.finish(**kw)
def draw_items(it, T=fitz.Matrix(1, 0, 0, 1, 0, 0)):
    for x in it["items"]:
        if x[0] == "l": sh.draw_line(fitz.Point(x[1]) * T, fitz.Point(x[2]) * T)
        elif x[0] == "re": sh.draw_rect(fitz.Rect(x[1]) * T)
        elif x[0] == "qu": sh.draw_quad(fitz.Quad(x[1]) * T)
        elif x[0] == "c": sh.draw_bezier(fitz.Point(x[1]) * T, fitz.Point(x[2]) * T, fitz.Point(x[3]) * T, fitz.Point(x[4]) * T)
        else: raise SystemExit(f"unknown item kind {x[0]}")
n_out = 0; plot_done = False; drawn_plot = []
def draw_plot():
    """The forest on the new axis, in the V28 layering: box, rule, then per row the interval and the marker."""
    global drawn_plot
    sh.draw_rect(fitz.Rect(X0_NEW, BOX.y0, BOX.x1, BOX.y1)); fin(box[0])
    l = RULE["items"][0]; sh.draw_line(fitz.Point(xn(1.0), l[1].y), fitz.Point(xn(1.0), l[2].y)); fin(RULE)
    for r28, v in zip(rows28, vals):
        c, m = r28["ci_item"], r28["mk_item"]; y = r28["y"]
        sh.draw_line(fitz.Point(xn(v["lo"]), y), fitz.Point(xn(v["hi"]), y)); fin(c)
        cx = xn(v["hr"]); s = r28["size"] / 2
        if r28["kind"] == "square": sh.draw_rect(fitz.Rect(cx - s, y - s, cx + s, y + s))
        else: sh.draw_circle(fitz.Point(cx, y), s)
        fin(m)
        drawn_plot.append(dict(outcome=v["outcome"], cell=v["cell"], y=round(y, 3), hr=v["hr"], lo=v["lo"], hi=v["hi"], x_hr=round(cx, 3), x_lo=round(xn(v["lo"]), 3), x_hi=round(xn(v["hi"]), 3), kind=r28["kind"], size=round(r28["size"], 3), colour=r28["col"], x28_hr=round(xo(r28["hr"]), 3)))
    sh.draw_line(fitz.Point(X0_NEW, SPINE_Y), fitz.Point(BOX.x1, SPINE_Y)); fin(SPINE)        # the bottom spine over the widened axis (V28: x 253.44 to 356.8 at y 613.22, 0.8 pt ink)
    for it, v in zip(ticks, TV):
        l = it["items"][0]; sh.draw_line(fitz.Point(xn(v), l[1].y), fitz.Point(xn(v), l[2].y)); fin(it)
for k, it in enumerate(raw):
    if k in STALE_BG: continue
    if k in inplot:
        if not plot_done: draw_plot(); plot_done = True
        continue
    draw_items(it); fin(it); n_out += 1
sh.commit()
# patch the join tokens
xr = pg.get_contents(); assert len(xr) == 1; cs = doc.xref_stream(xr[0]).decode("latin-1"); n_tok = cs.count("{lineJoin} j"); assert n_tok == len(JOINS), (n_tok, len(JOINS))
parts = cs.split("{lineJoin} j"); cs2 = "".join(parts[i] + (f"{JOINS[i]} j" if i < len(JOINS) else "") for i in range(len(parts))); doc.update_stream(xr[0], cs2.encode("latin-1"))
# text
writers = {}; placed = []
for s in kept + ticklab:
    rot = s["dir"] == [0.0, -1.0]; key = (s["color"], "Bold" in s["font"], rot); tw = writers.setdefault(key, fitz.TextWriter(pg.rect))
    if s in ticklab:
        v = float(s["text"]); wdt = FR.text_length(s["text"], fontsize=s["size"]); x = xn(v) - wdt / 2; placed.append(dict(text=s["text"], x=round(x, 3), y=s["origin"][1], centre=round(xn(v), 3), size=s["size"]))
    else: x = s["origin"][0]
    tw.append(fitz.Point(x, s["origin"][1]), s["text"], font=FB if key[1] else FR, fontsize=s["size"])
for (col, bold, rot), tw in writers.items():
    rgb = tuple(int(col.lstrip("#")[i:i + 2], 16) / 255 for i in (0, 2, 4))
    if rot:
        q = [s for s in kept if s["dir"] == [0.0, -1.0]]; assert len(q) == 1; q = q[0]; tw.write_text(pg, color=rgb, morph=(fitz.Point(q["origin"][0], q["origin"][1]), fitz.Matrix(90)))
    else: tw.write_text(pg, color=rgb)
n_cmap = 0; seen = set()
for f in pg.get_fonts(full=True):
    if f[3] not in ("Arial Regular", "Arial Bold") or f[0] in seen: continue
    seen.add(f[0]); tu = doc.xref_get_key(f[0], "ToUnicode")
    if tu[0] != "xref": continue
    tx = int(tu[1].split()[0]); st = doc.xref_stream(tx).decode("latin-1"); st2 = st.replace("<0003> <00a0>", "<0003> <0020>").replace("<0010> <00ad>", "<0010> <002d>")
    if st2 != st: doc.update_stream(tx, st2.encode("latin-1")); n_cmap += 1
doc.save(OUT, garbage=3, deflate=True); doc.close(); src.close(); gc.collect()
LAB_END = max(s["bbox"][2] for s in kept if abs(s["origin"][0] - 21.0) < 0.3 and s["size"] == 11.0)
rec = dict(sheet="ED_Fig06", round=51, input=IN, input_sha256=hashlib.sha256(open(IN, "rb").read()).hexdigest(), output=OUT, output_sha256=hashlib.sha256(open(OUT, "rb").read()).hexdigest(), page=[W, H],
           plot_region=[PLOT.x0, PLOT.y0, PLOT.x1, PLOT.y1], box_v28=[BOX.x0, BOX.y0, BOX.x1, BOX.y1], box_new=[X0_NEW, BOX.y0, BOX.x1, BOX.y1], axis_v28=dict(a=a_, b=b_, ticks_x=TX), axis_new=dict(a=a2, b=b2, ticks_x=[xn(v) for v in TV]), axis_range=[LO, HI],
           rule_x=[xo(1.0), xn(1.0)], removed=[dict(text=s["text"], x=s["origin"][0], y=s["origin"][1], size=s["size"], bbox=s["bbox"]) for s in removed], tick_labels=placed, rows=drawn_plot, row_label_right_edge=LAB_END, leftmost_ink_new=round(min(r["x_lo"] for r in drawn_plot), 3),
           n_items_v28=len(raw), n_items_in_plot_region=len(inplot), n_items_replayed_outside=n_out, n_plot_items_drawn=2 + 48 + 1 + 4, spine=dict(y=SPINE_Y, x_v28=[SPINE["rect"].x0, SPINE["rect"].x1], x_new=[X0_NEW, BOX.x1], width=SPINE["width"]), skipped_stale_layer=[list(raw[k]["rect"]) for k in STALE_BG], joins_patched=n_tok, cmaps_patched=n_cmap,
           values_source=dict(file=f"{SRC}/results_crossed.csv", sidecar_mtime=sc["_gate"]["output_mtime"], sha256=sc["_gate"]["sha256"]), v28_marks_vs_file_worst_pt=round(worst, 4))
json.dump(rec, open(f"{WORK}/r51_record.json", "w"), indent=1, ensure_ascii=False)
print(f"stale layer skipped: {[list(raw[k]['rect']) for k in STALE_BG]}"); print(f"axis: box x0 {BOX.x0:.2f} -> {X0_NEW}, ticks {[round(x, 2) for x in TX]} -> {[round(xn(v), 2) for v in TV]}, rule {xo(1.0):.2f} -> {xn(1.0):.2f}; {n_out} items replayed outside the plot, {len(inplot)} plot-region items replaced by {2 + 48 + 1 + 4} drawn (box, rule, spine, 24 intervals, 24 markers, 4 ticks); {len(removed)} strings dropped, {len(ticklab)} tick labels re-centred; joins patched {n_tok}; wrote {OUT}")
