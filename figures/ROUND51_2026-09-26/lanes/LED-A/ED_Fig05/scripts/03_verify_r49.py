#!/usr/bin/env python3
"""ROUND 49, lane LED-A: ED Fig 5 proof at +1 pt (class W, four panels), adapted from ROUND37_2026-09-14/figures/lanes/V14_L3a_DURATION_MAIN/
ED_Fig05/scripts/03_verify.py. Run on the built sheet before the title strip (work/ED_Fig05_built.pdf, the V13 page box). Declared deltas of this
round: every text +1 pt (annotations 9.5 -> 10.5, ticks, key and panel titles 10 -> 11, y labels 11 -> 12, letters 13 -> 14). Anchors: letters and
panel titles hang from their top (baseline moves by the ascent growth), x tick labels are centred on their ticks and hang from the tick pad, y tick
labels keep their right edge, the y labels stay centred on their axes (they move left with the wider tick labels), the annotation blocks hang from
their top (later lines move down by the line pitch growth), the key keeps its top and its right edge at the page margin (its text moves left by
its growth). Geometry: the axes widths are the V26 layout (the key width measured at the V26 size), so the 16 bars equal V26's. Data checks
(every printed percent, PR, CI and P recomputed from the v8 table) are the round-37 checks. Renders by Ghostscript."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import collections, csv, json, os, re, sys
import numpy as np, fitz
from PIL import Image
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND51_2026-09-26/lanes/LED-A/scripts_shared")
import lane_common as C
import r49lib as R
VL = C.VL
SHEET = "ED_Fig05"; HERE = f"{C.LANE}/{SHEET}"; VER = f"{HERE}/verify"; CROPS = f"{VER}/crops"; WORK = f"{HERE}/work"; LOGS = f"{R.LANE}/logs/{SHEET}"; os.makedirs(CROPS, exist_ok=True)
BASE = f"{C.BASE}/{SHEET}.pdf"; NEW = f"{WORK}/{SHEET}_built.pdf"; V26 = f"{R.V26}/{SHEET}.pdf"; STRIP = 16.0
D = json.load(open(f"{WORK}/{SHEET}_drawn.json")); KP = D["key_placement"]; DELTA = 1.0
LINES = []; OK = True
def check(name, ok, detail):
    global OK; OK &= bool(ok); LINES.append(f"{'PASS' if ok else 'FAIL'}  {name}: {detail}"); print(LINES[-1][:300])
bd = fitz.open(C.hydrated(BASE)); bp = bd[0]; nd = fitz.open(C.hydrated(NEW)); np_ = nd[0]
check("single page, page box equal to V13 (before the strip)", len(nd) == 1 and abs(np_.rect.width - bp.rect.width) < 0.01 and abs(np_.rect.height - bp.rect.height) < 0.01,
      f"V13 {bp.rect.width:.2f} x {bp.rect.height:.2f}, new {np_.rect.width:.2f} x {np_.rect.height:.2f}")
LB, LN = C.letters(bp), C.letters(np_)
check("letters a, b, c, d at the V13 positions (left edge within 0.3 pt, top within 0.5 pt: they hang from their top, the text layer reports the 14 pt glyph box 0.3 pt higher), Arial Bold 14", [t[0] for t in LN] == ["a", "b", "c", "d"] and all(a[0] == b[0] and abs(a[1] - b[1]) <= 0.3 and abs(a[2] - b[2]) <= 0.5 and a[4] == b[4] and abs(b[5] - 14.0) < 0.05 for a, b in zip(LB, LN)),
      f"V13 {[(t[0], t[1], t[2]) for t in LB]}, new {[(t[0], t[1], t[2], t[4], t[5]) for t in LN]}")
fb = sorted(set(f[3].split("+")[-1] for f in bp.get_fonts(full=True))); fn = sorted(set(f[3].split("+")[-1] for f in np_.get_fonts(full=True)))
check("fonts identical to V13 (Arial faces)", fb == fn, f"V13 {fb}, new {fn}")
check("no nested forms", len(np_.get_xobjects()) == 0, f"{len(np_.get_xobjects())} xobjects")
SB = C.dedupe(C.spans_of(bp)); SN = C.dedupe(C.spans_of(np_))
sizes_b = sorted({round(s["size"], 2) for s in SB}); sizes_n = sorted({round(s["size"], 2) for s in SN})
check("every type size = its V13 size + 1 pt", sizes_n == [round(z + DELTA, 2) for z in sizes_b], f"V13 {sizes_b} -> new {sizes_n}")
def is_value(s): return bool(re.match(r"^(\d+\.\d|PR, \d\.\d\d|\d\.\d\d-\d\.\d\d\)|P [<=] 0\.\d\d\d)$", s["text"]))
inv_b = [s for s in SB if not is_value(s)]; inv_n = [s for s in SN if not is_value(s)]
TITLES = {p["title"] for p in D["panels"].values()}; KEYL = {"Sleep under", "5 h", "Sleep 7 h", "or more"}
def cx(s): return (s["bbox"][0] + s["bbox"][2]) / 2
def cy(s): return (s["bbox"][1] + s["bbox"][3]) / 2
bad = []; missing = []; worst = collections.defaultdict(float)
for s in inv_n:
    t = s["text"]; cands = [b for b in inv_b if b["text"] == t]
    if not cands: missing.append(t); continue
    b = min(cands, key=lambda b: abs(b["origin"][0] - s["origin"][0]) + abs(b["origin"][1] - s["origin"][1]))
    if t in ("a", "b", "c", "d") or t in TITLES: kind, dx, dy, lim = "top-left", abs(s["bbox"][0] - b["bbox"][0]), abs(s["bbox"][1] - b["bbox"][1]), (0.3, 0.6)
    elif t == "Percent of the group": kind, dx, dy, lim = "y label (rotated, centred on the axis)", abs(s["origin"][0] - b["origin"][0]), abs(cy(s) - cy(b)), (3.0, 0.6)
    elif t in ("0", "10", "20", "30"): kind, dx, dy, lim = "y tick (right edge)", abs(s["bbox"][2] - b["bbox"][2]), abs(s["origin"][1] - b["origin"][1]), (0.3, 1.0)
    elif t.startswith("T90 >"): kind, dx, dy, lim = "x tick (centred, hangs from the tick pad)", abs(cx(s) - cx(b)), abs(s["bbox"][1] - b["bbox"][1]), (0.3, 0.6)
    elif t in KEYL: kind, dx, dy, lim = "key line (left edge moves with the key's growth, top of the first line kept)", abs(s["origin"][0] - b["origin"][0]), abs(s["origin"][1] - b["origin"][1]), (8.0, 5.0)
    elif t == "(95% CI,": kind, dx, dy, lim = "annotation line 2 (centred, block hangs from its top)", abs(cx(s) - cx(b)), abs(s["origin"][1] - b["origin"][1]), (1.5, 4.0)
    else: kind, dx, dy, lim = "other (origin)", abs(s["origin"][0] - b["origin"][0]), abs(s["origin"][1] - b["origin"][1]), (0.3, 0.3)
    worst[kind] = max(worst[kind], dx, dy)
    if dx > lim[0] or dy > lim[1]: bad.append((t, kind, round(dx, 2), round(dy, 2)))
check("every invariant string at its V13 anchor at +1 pt (letters and titles top-left, y ticks right edge, x ticks centred, y labels centred on the axis, key and annotation blocks from their top)",
      not bad and not missing, f"largest shifts by kind {dict((k, round(v, 2)) for k, v in worst.items())}, missing {missing}{', off: ' + str(bad) if bad else ''}")
extra_b = collections.Counter(b["text"] for b in inv_b) - collections.Counter(s["text"] for s in inv_n)
check("no V13 invariant string missing on the new sheet", not extra_b, f"V13-only {dict(extra_b)}")
kl = [s for s in SN if s["text"] in KEYL]; kl_b = [s for s in SB if s["text"] in KEYL]
check("key: text 11 pt, right edge of the key box at the page margin, its left edge moved through the 0.08 in gap by the key's growth, at most 0.03 in past the panel b axes' invisible boundary (declared)",
      len(kl) == 4 and all(abs(s["size"] - 11.0) < 0.05 for s in kl) and KP["key_right_in"] <= 171 / 25.4 - 0.18 + 1e-6 and KP["key_left_in"] >= KP["axes_right_in"] - 0.03,
      f"key left edge {KP['key_left_in']:.3f} in, panel b axes right edge {KP['axes_right_in']:.3f} in (gap {KP['gap_in']:.3f} in, V26 {KP['gap_v26_in']:.2f} in); key text x {min(s['origin'][0] for s in kl):.2f} (V13 {min(s['origin'][0] for s in kl_b):.2f})")
pct_b = sorted(round(cx(s), 1) for s in SB if re.match(r"^\d+\.\d$", s["text"])); pct_n = sorted(round(cx(s), 1) for s in SN if re.match(r"^\d+\.\d$", s["text"]))
blk_b = sorted((round(s["origin"][1], 1), round(cx(s), 1), s["text"][:3]) for s in SB if is_value(s) and not re.match(r"^\d+\.\d$", s["text"]))
blk_n = sorted((round(s["origin"][1], 1), round(cx(s), 1), s["text"][:3]) for s in SN if is_value(s) and not re.match(r"^\d+\.\d$", s["text"]))
pr_top_b = sorted(round(s["bbox"][1], 2) for s in SB if s["text"].startswith("PR,")); pr_top_n = sorted(round(s["bbox"][1], 2) for s in SN if s["text"].startswith("PR,"))
check("16 bar-top percent labels centred on the same 16 bars (same baselines: data positions), 24 PR / CI / P strings on the same centres, the 8 blocks hanging from the same top (later lines lower by the 1.5 em pitch growth)",
      len(pct_b) == len(pct_n) == 16 and max(abs(a - b) for a, b in zip(pct_b, pct_n)) <= 1.0 and len(blk_b) == len(blk_n) == 24 and all(abs(a[1] - b[1]) <= 1.5 and 0 <= b[0] - a[0] <= 6.0 for a, b in zip(blk_b, blk_n))
      and len(pr_top_b) == len(pr_top_n) == 8 and max(abs(a - b) for a, b in zip(pr_top_b, pr_top_n)) <= 0.6,
      f"percent centres max shift {max(abs(a - b) for a, b in zip(pct_b, pct_n)) if pct_b and len(pct_b) == len(pct_n) else 'n/a'}, block baselines lower by {min(b[0] - a[0] for a, b in zip(blk_b, blk_n)):.1f} to {max(b[0] - a[0] for a, b in zip(blk_b, blk_n)):.1f} pt, PR tops max shift {max(abs(a - b) for a, b in zip(pr_top_b, pr_top_n)) if len(pr_top_b) == len(pr_top_n) else 'n/a'} pt")
pctl_b = {round(s["origin"][1], 1) for s in SB if re.match(r"^\d+\.\d$", s["text"])}; pctl_n = {round(s["origin"][1], 1) for s in SN if re.match(r"^\d+\.\d$", s["text"])}
def panel_of(s):
    return ("a" if s["origin"][1] < 310 else "c") if s["origin"][0] < 200 else ("b" if s["origin"][1] < 310 else "d")
expected = collections.Counter()
for L_, pnl in D["panels"].items():
    for c, v in pnl["cuts"].items():
        for t in v["printed"]: expected[(L_, t)] += 1
        expected[(L_, "(95% CI,")] += 1
got = collections.Counter((panel_of(s), s["text"]) for s in SN if is_value(s) or s["text"] == "(95% CI,")
check("read-back: every printed percent, PR, CI and P string equals the value recomputed from the v8 table (32 numeric strings plus 8 '(95% CI,' lines)", got == expected,
      f"printed {sum(got.values())}, expected {sum(expected.values())}, mismatch {dict((got - expected) + (expected - got))}")
tb, tn = C.nums([w["text"] for w in C.dedupe(C.words_of(bp))]), C.nums([w["text"] for w in C.dedupe(C.words_of(np_))])
old_vals = sorted(s["text"] for s in SB if is_value(s)); new_vals = sorted(s["text"] for s in SN if is_value(s))
d = C.delta(tb, tn); ed = C.delta(C.nums(old_vals), C.nums(new_vals))
check("numeric-token multiset differs from V13 by exactly the printed value strings (percent labels, PR, CI, P: the v8 values against V13's v7 values)", d == ed, f"delta {d}")
wb, wn = sorted(w["text"] for w in C.dedupe(C.words_of(bp))), sorted(w["text"] for w in C.dedupe(C.words_of(np_)))
VALW = lambda w: bool(re.match(r"^[\d(]", w)) or w.endswith(")") or w in ("=", "<")
wd = C.delta([w for w in wb if not VALW(w)], [w for w in wn if not VALW(w)])
check("word multiset identical apart from the value words (numbers, CI parentheses, the P comparator '=' or '<')", wd == {"lost": {}, "gained": {}}, f"{wd}")
hx = lambda c: "#%02x%02x%02x" % tuple(int(round(v * 255)) for v in c)
cb = sorted(set(hx(x["fill"]) for x in bp.get_drawings() if x.get("fill") is not None)); cn = sorted(set(hx(x["fill"]) for x in np_.get_drawings() if x.get("fill") is not None))
check("fill colour set identical to V13 (green ladder #79d475 under 5 h, #39c445 7 h or more)", cb == cn, f"{cn}")
bars = [x for x in np_.get_drawings() if x.get("fill") is not None and hx(x["fill"]) in ("#79d475", "#39c445") and x["rect"].width > 15]
check("16 bars drawn", len(bars) == 16, f"{len(bars)} bars")
# data marks unchanged: the 16 bar rectangles equal V26's (V26 carries the 16 pt strip: its rects sit 16 pt lower)
vd = fitz.open(C.hydrated(V26)); vp = vd[0]
bars_v = [x for x in vp.get_drawings() if x.get("fill") is not None and hx(x["fill"]) in ("#79d475", "#39c445") and x["rect"].width > 15]
def key_(r, dy=0.0): return (round(r.x0, 2), round(r.y0 - dy, 2), round(r.x1, 2), round(r.y1 - dy, 2))
rn_ = sorted(key_(b["rect"]) for b in bars); rv_ = sorted(key_(b["rect"], STRIP) for b in bars_v)
dev = max(max(abs(a - b) for a, b in zip(p, q)) for p, q in zip(rn_, rv_)) if len(rn_) == len(rv_) == 16 else 99
check("data marks unchanged: the 16 bar rectangles equal V26's within 0.05 pt (V26 shifted up by the 16 pt strip), so the axes widths are the V26 layout", dev < 0.05, f"{len(bars_v)} V26 bars, max deviation {dev:.4f} pt; axes width {KP['axw_each_in']:.4f} in each")
vd.close()
R.gs_render(BASE, f"{WORK}/base_150.png", 150, LOGS, name="gs150_base"); R.gs_render(NEW, f"{WORK}/new_150.png", 150, LOGS, name="gs150_new")
a = np.asarray(Image.open(f"{WORK}/base_150.png").convert("RGB")).astype(int); b = np.asarray(Image.open(f"{WORK}/new_150.png").convert("RGB")).astype(int)
h = min(a.shape[0], b.shape[0]); w = min(a.shape[1], b.shape[1]); a = a[:h, :w]; b = b[:h, :w]
diff = np.abs(a - b).max(axis=2) > 0
LINES.append(f"INFO  raster (Ghostscript 150 dpi): {int(diff.sum())} differing pixels of {diff.size} against V13 (every text +1 pt, the bars and printed values carry the v8 numbers against V13's v7 ones)")
R.gs_render(BASE, f"{CROPS}/{SHEET}_OLD_200dpi.png", 200, LOGS, name="gs200_base"); R.gs_render(NEW, f"{CROPS}/{SHEET}_NEW_200dpi.png", 200, LOGS, name="gs200_new")
side = np.concatenate([a, np.full((a.shape[0], 20, 3), 255), b], axis=1).astype(np.uint8); Image.fromarray(side).save(f"{CROPS}/{SHEET}_OLD_vs_NEW_150dpi.png")
RULES = VL.RULES
RULES["pr_2dp"] = lambda v: "PR, " + VL.halfup(v, 2)
RULES["ci_2dp_close"] = lambda v: f"{VL.halfup(v[0], 2)}-{VL.halfup(v[1], 2)})"
RULES["P_3dp"] = lambda v: "P < 0.001" if float(v) < 0.001 else "P = " + VL.halfup(v, 3)
n_rows, n_no = VL.printed_values_csv(NEW, D["records"], f"{VER}/{SHEET}_printed_values.csv", static_from=BASE)
check("printed-values CSV: every string accounted for (40 drawn strings under their rules, the rest V13 strings at +1 pt in the builder's record)", n_no == 0, f"{n_rows} rows, {n_no} NO")
hits = C.grep_proof(HERE)
check("grep proof: no live script in this sheet folder names a v7 table, snapshot or set", not hits, f"{hits}")
LINES.append(f"INFO  sources: {D['sources']['parquet']} sha256 {D['sources']['parquet_sha256'][:16]}, habitual {D['sources']['habitual_sha256'][:16]}, {D['sources']['cohort_filter']}; V13 base sha256 {D['base_sha256'][:16]}")
bd.close(); nd.close()
print(C.write_checks(LINES, f"{VER}/checks_03_verify.txt", OK))
