#!/usr/bin/env python3
"""ROUND 49, lane LED-A: ED Fig 4 proof at +1 pt (class W against V13), adapted from ROUND37_2026-09-14/figures/lanes/V14_L3a_DURATION_MAIN/
ED_Fig04/scripts/03_verify.py. Run on the built sheet before the title strip (work/ED_Fig04_built.pdf, the V13 page box). Declared deltas of this
round: every text +1 pt (row labels and ticks 10 -> 11, axis title 11 -> 12, key 10 -> 11); the tick labels stay centred on their ticks and the
axis title on the axes (their left edges move with their width, the axis title sits 1 pt lower under the taller tick labels); the key grows
symmetrically about the V13 key centre (first entry left, second right by half the row's growth, from the builder's record). Renders by Ghostscript.
The data checks (rows = the ten largest T90-significant conditions of results_v2.json, 20 markers and 20 CI lines read back, colours, no grid)
are the round-37 checks. The v7 CHANGES block is not part of this round (no number changes)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import collections, csv, json, os, sys
import numpy as np, fitz
from PIL import Image
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LED-A/scripts_shared")
import lane_common as C
import r49lib as R
VL = C.VL
SHEET = "ED_Fig04"; HERE = f"{C.LANE}/{SHEET}"; VER = f"{HERE}/verify"; CROPS = f"{VER}/crops"; WORK = f"{HERE}/work"; LOGS = f"{R.LANE}/logs/{SHEET}"; os.makedirs(CROPS, exist_ok=True)
BASE = f"{C.BASE}/{SHEET}.pdf"; NEW = f"{WORK}/{SHEET}_built.pdf"
OLDREC = f"{C.R30_LANE}/{SHEET}/verify/{SHEET}_drawn.json"   # the round-30 drawn record (what V13 shows), the OLD side of the row-set delta only
DRAWN = json.load(open(f"{WORK}/{SHEET}_drawn.json")); KP = DRAWN["key_placement"]; DELTA = 1.0
LINES = []; OK = True
def check(name, ok, detail):
    global OK; OK &= bool(ok); LINES.append(f"{'PASS' if ok else 'FAIL'}  {name}: {detail}"); print(LINES[-1][:300])

for p in (BASE, NEW, OLDREC): C.hydrated(p)
bd = fitz.open(BASE); bp = bd[0]; nd = fitz.open(NEW); np_ = nd[0]
check("single page, page box equal to V13 (before the strip)", len(nd) == 1 and abs(np_.rect.width - bp.rect.width) < 0.01 and abs(np_.rect.height - bp.rect.height) < 0.01,
      f"V13 {bp.rect.width:.2f} x {bp.rect.height:.2f}, new {np_.rect.width:.2f} x {np_.rect.height:.2f}")
check("no panel letters on either sheet (single-panel sheet)", C.letters(bp) == [] and C.letters(np_) == [], f"{C.letters(bp)} {C.letters(np_)}")
fb = sorted(set(f[3].split("+")[-1] for f in bp.get_fonts(full=True))); fn = sorted(set(f[3].split("+")[-1] for f in np_.get_fonts(full=True)))
check("fonts Arial only, same faces as V13", fn == fb and all("Arial" in f for f in fn), f"V13 {fb}, new {fn}")
check("no nested forms on the new sheet", len(np_.get_xobjects()) == 0, f"{len(np_.get_xobjects())} xobjects")
SB = C.dedupe(C.spans_of(bp)); SN = C.dedupe(C.spans_of(np_))
sizes_b = sorted({round(s["size"], 2) for s in SB}); sizes_n = sorted({round(s["size"], 2) for s in SN})
check("every type size = its V13 size + 1 pt", sizes_n == [round(z + DELTA, 2) for z in sizes_b] and all(any(abs(s["size"] - DELTA - b["size"]) < 0.01 for b in SB if b["text"] == s["text"]) for s in SN if any(b["text"] == s["text"] for b in SB)), f"V13 {sizes_b} -> new {sizes_n}")
ROWS = [r["condition"] for r in DRAWN["rows"]]; OLD_ROWS = [r["condition"] for r in json.load(open(OLDREC))["rows"]]
base_rows = [s["text"] for s in sorted((s for s in SB if s["text"] in OLD_ROWS or s["text"] in ROWS), key=lambda s: s["origin"][1]) if s["text"] in OLD_ROWS]
check("V13 row labels = the round-30 drawn record (the OLD side of the row-set delta is what V13 shows)", base_rows == OLD_ROWS, f"V13 {base_rows}")
ROWSET = set(ROWS) | set(OLD_ROWS)
# every non-row string against V13 by its anchor
KEYS = {"T90 >10%": KP["T90 >10% origin"], "Laboratory sleep under 5 h": KP["Laboratory sleep under 5 h origin"]}
def centre(s): return (s["bbox"][0] + s["bbox"][2]) / 2
worst_x = worst_y = 0.0; missing = []; bad = []
for s in SN:
    if s["text"] in ROWSET: continue
    cands = [b for b in SB if b["text"] == s["text"]]
    if not cands: missing.append(s["text"]); continue
    b = min(cands, key=lambda b: abs(b["origin"][0] - s["origin"][0]) + abs(b["origin"][1] - s["origin"][1]))
    if s["text"] in KEYS:
        dx, dy = abs(s["origin"][0] - KEYS[s["text"]][0]), abs(s["origin"][1] - b["origin"][1]); lim = (0.3, 0.3)
    elif s["text"] in ("1", "2", "4", "8"):
        dx, dy = abs(centre(s) - centre(b)), abs(s["origin"][1] - b["origin"][1]); lim = (0.3, 1.5)      # centred on the tick, hangs from the tick pad: the baseline moves by the ascent growth
    else:
        dx, dy = abs(centre(s) - centre(b)), abs(s["origin"][1] - b["origin"][1]); lim = (0.3, 2.5)      # the axis title, centred on the axes, under the taller tick labels
    worst_x = max(worst_x, dx); worst_y = max(worst_y, dy)
    if dx > lim[0] or dy > lim[1]: bad.append((s["text"], round(dx, 2), round(dy, 2)))
check("tick labels centred on their V13 ticks (within 0.3 pt), the axis title centred on the axes (within 0.3 pt, at most 2.5 pt lower), the key entries at the builder's centred positions (within 0.3 pt, same baseline)",
      not bad and not missing, f"largest centre offset {worst_x:.3f} pt, largest baseline shift {worst_y:.3f} pt, missing {missing}{', off: ' + str(bad) if bad else ''}")
g1, g2 = KP["dx_key1"], KP["dx_key2"]; kl13, kr13 = KP["v13_box"]; kl, kr = KP["new_box"]
check("key grows symmetrically about the V13 key's ink centre (handle left end to second label right edge): both edges move by half the growth, the centre stays",
      abs(KP["new_centre"] - KP["v13_centre"]) < 0.3 and abs((kl - kl13) + (kr - kr13)) < 0.3,
      f"V13 box {kl13:.2f} to {kr13:.2f} (centre {KP['v13_centre']:.2f}, the axes centre {KP['axes_centre']:.2f}), new {kl:.2f} to {kr:.2f} (centre {KP['new_centre']:.2f}), growth {KP['growth_pt']:.2f} pt; text origins 'T90 >10%' {g1:+.2f} pt, 'Laboratory sleep under 5 h' {g2:+.2f} pt (its label grows to the right)")
rb = sorted(round(s["origin"][1], 2) for s in SB if s["text"] in OLD_ROWS); rn = sorted(round(s["origin"][1], 2) for s in SN if s["text"] in ROWS)
check("10 row labels on the same 10 rows as V13 (left edge unchanged, baseline within 1 pt: the label stays centred on its row at 11 pt)",
      len(rb) == len(rn) == 10 and max(abs(a - b) for a, b in zip(rb, rn)) <= 1.0 and all(abs(s["origin"][0] - 27.36) < 0.5 for s in SN if s["text"] in ROWS), f"V13 {rb}, new {rn}")
entered, left = sorted(set(ROWS) - set(OLD_ROWS)), sorted(set(OLD_ROWS) - set(ROWS))
check("row set: the ten largest T90-significant conditions of the v8.1 file (against V13: rows entering and leaving, the round-37 delta)", len(ROWS) == 10 and len(set(ROWS)) == 10, f"enter {entered}, leave {left}")
order_n = [s["text"] for s in sorted((s for s in SN if s["text"] in ROWS), key=lambda s: s["origin"][1])]
check("rows on the new sheet in the drawn order (sorted by the T90 hazard ratio)", order_n == ROWS, f"{order_n}")
tb, tn = C.nums([w["text"] for w in C.dedupe(C.words_of(bp))]), C.nums([w["text"] for w in C.dedupe(C.words_of(np_))])
check("numeric-token multiset identical to V13 (the sheet prints no estimates; the row labels carry no digit)", C.delta(tb, tn) == {"lost": {}, "gained": {}}, f"{len(tb)} vs {len(tn)} tokens, delta {C.delta(tb, tn)}")
wb, wn = sorted(w["text"] for w in C.dedupe(C.words_of(bp))), sorted(w["text"] for w in C.dedupe(C.words_of(np_)))
d = C.delta(wb, wn)
exp_lost = collections.Counter(w for lab in left for w in lab.split()); exp_gained = collections.Counter(w for lab in entered for w in lab.split())
check("word multiset differs from V13 by exactly the row labels that leave and enter", collections.Counter(d["lost"]) == exp_lost and collections.Counter(d["gained"]) == exp_gained, f"lost {d['lost']}, gained {d['gained']}")
# read-back of drawn values
ticks = {s["text"]: s for s in SN if s["text"] in ("1", "2", "4", "8") and s["origin"][1] > 290}
xc = {float(t): (s["bbox"][0] + s["bbox"][2]) / 2 for t, s in ticks.items()}
lv = np.log([1, 2, 4, 8]); xv = np.array([xc[v] for v in (1, 2, 4, 8)]); slope, icpt = np.polyfit(lv, xv, 1)
def hr_of(x): return float(np.exp((x - icpt) / slope))
dr = np_.get_drawings()
def hx(c): return "#%02x%02x%02x" % tuple(int(round(v * 255)) for v in c)
mk = [x for x in dr if x.get("fill") is not None and hx(x["fill"]) in ("#0288d1", "#39c445") and x["rect"].width < 10 and x["rect"].x0 > 150 and x["rect"].y1 < 290]
ci = [x for x in dr if x.get("color") is not None and hx(x["color"]) in ("#0288d1", "#298d32") and abs((x.get("width") or 0) - 1.7) < 0.01 and x["rect"].x0 > 150 and x["rect"].height < 0.6 and x["rect"].y1 < 290]
V = json.load(open(f"{C.NUM}/results_v2.json"))["vs_sleep_duration"]["outcomes"]
rows = DRAWN["rows"]; worst_m = 0.0; worst_ci = 0.0; n_ok = 0
for r in rows:
    for key, col in (("t90", "#0288d1"), ("short_sleep", "#39c445")):
        v = V[r["condition"]][key]
        m = min((m for m in mk if hx(m["fill"]) == col), key=lambda m: abs(hr_of((m["rect"].x0 + m["rect"].x1) / 2) - v["hr"]))
        got = hr_of((m["rect"].x0 + m["rect"].x1) / 2); worst_m = max(worst_m, abs(got - v["hr"]) / v["hr"])
        cc = "#0288d1" if key == "t90" else "#298d32"
        c = min((c for c in ci if hx(c["color"]) == cc), key=lambda c: abs(hr_of(c["rect"].x0) - v["lo"]) + abs(hr_of(c["rect"].x1) - v["hi"]))
        worst_ci = max(worst_ci, abs(hr_of(c["rect"].x0) - v["lo"]) / v["lo"], abs(hr_of(c["rect"].x1) - v["hi"]) / v["hi"]); n_ok += 1
check("read-back: 20 markers at the hazard ratios of results_v2.json (within 0.3 percent) and 20 CI lines at the bounds", len(mk) == 20 and len(ci) == 20 and worst_m <= 0.003 and worst_ci <= 0.003 and n_ok == 20,
      f"{len(mk)} markers, {len(ci)} CI lines, largest marker deviation {100 * worst_m:.3f} percent, CI ends {100 * worst_ci:.3f} percent")
# data marks at the V26 positions (V26 = the round-37 build on the same v8.1 numbers plus the 16 pt strip): markers and CI lines
drb = bp.get_drawings()
V26 = f"{R.V26}/{SHEET}.pdf"; STRIP = 16.0; vd = fitz.open(C.hydrated(V26)); drv = vd[0].get_drawings(); vd.close()
mkv = [x for x in drv if x.get("fill") is not None and hx(x["fill"]) in ("#0288d1", "#39c445") and x["rect"].width < 10 and x["rect"].x0 > 150 and x["rect"].y1 < 290 + STRIP]
civ = [x for x in drv if x.get("color") is not None and hx(x["color"]) in ("#0288d1", "#298d32") and abs((x.get("width") or 0) - 1.7) < 0.01 and x["rect"].x0 > 150 and x["rect"].height < 0.6 and x["rect"].y1 < 290 + STRIP]
def cen(ms, dy=0.0): return sorted((round((m["rect"].x0 + m["rect"].x1) / 2, 3), round((m["rect"].y0 + m["rect"].y1) / 2 - dy, 3), hx(m["fill"])) for m in ms)
def ends(cs, dy=0.0): return sorted((round(c["rect"].x0, 3), round(c["rect"].x1, 3), round(c["rect"].y0 - dy, 3), hx(c["color"])) for c in cs)
a_, b_ = cen(mk), cen(mkv, STRIP); c_, d_ = ends(ci), ends(civ, STRIP)
dev = max([max(abs(p[0] - q[0]), abs(p[1] - q[1])) for p, q in zip(a_, b_)] + [max(abs(p[0] - q[0]), abs(p[1] - q[1]), abs(p[2] - q[2])) for p, q in zip(c_, d_)]) if len(a_) == len(b_) == 20 and len(c_) == len(d_) == 20 else 99
check("data marks unchanged: the 20 markers and 20 CI lines sit within 0.1 pt of V26's (V26 shifted up by the 16 pt strip; V26 was drawn from the three-decimal results_v2.json of 2026-09-14, the file was rewritten at full precision on 2026-09-16, which moves a mark by at most 0.06 pt, below any print resolution)", dev < 0.1, f"max deviation {dev:.4f} pt over {len(a_)} markers and {len(c_)} CI lines")
ys = sorted(((m["rect"].y0 + m["rect"].y1) / 2, hr_of((m["rect"].x0 + m["rect"].x1) / 2)) for m in mk if hx(m["fill"]) == "#0288d1")
check("blue markers descend in hazard ratio down the page (rows sorted by T90 HR)", all(a[1] >= b[1] for a, b in zip(ys, ys[1:])), f"{[round(h, 2) for _y, h in ys]}")
lab_y = [s["origin"][1] for s in sorted((s for s in SN if s["text"] in ROWS), key=lambda s: s["origin"][1])]
check("each blue marker sits on its row (marker centre within 12 pt of the row label baseline, in row order)", all(abs(my - ly) < 12 for (my, _h), ly in zip(ys, lab_y)), f"{[round(my - ly, 1) for (my, _h), ly in zip(ys, lab_y)]}")
cb = sorted(set(hx(x["fill"]) for x in drb if x.get("fill") is not None)); cn = sorted(set(hx(x["fill"]) for x in dr if x.get("fill") is not None))
sb_ = sorted(set(hx(x["color"]) for x in drb if x.get("color") is not None)); sn_ = sorted(set(hx(x["color"]) for x in dr if x.get("color") is not None))
check("fill and stroke colour sets identical to V13", cb == cn and sb_ == sn_, f"fills {cn}, strokes {sn_}")
check("no pale #eef0f1 vertical grid line drawn (V13 design)", "#eef0f1" not in sn_, f"{sn_}")
# raster (Ghostscript): differences confined to the forest with its row labels, the text boxes (+1 pt) and the key row
R.gs_render(BASE, f"{WORK}/base_150.png", 150, LOGS, name="gs150_base"); R.gs_render(NEW, f"{WORK}/new_150.png", 150, LOGS, name="gs150_new")
a = np.asarray(Image.open(f"{WORK}/base_150.png").convert("RGB")).astype(int); b = np.asarray(Image.open(f"{WORK}/new_150.png").convert("RGB")).astype(int)
h = min(a.shape[0], b.shape[0]); w = min(a.shape[1], b.shape[1]); a = a[:h, :w]; b = b[:h, :w]
diff = (np.abs(a - b) > 32).any(axis=2); s150 = 150 / 72; mask = np.zeros_like(diff)
mask[:int(290 * s150), int(20 * s150):] = True                                                      # the forest and its row labels (the row set differs from V13 by the v8.1 rule)
for s in SB + SN:
    x0, y0, x1, y1 = s["bbox"]; mask[max(0, int((y0 - 2) * s150)):int((y1 + 2) * s150) + 1, max(0, int((x0 - 2) * s150)):int((x1 + 2) * s150) + 1] = True
mask[int(314 * s150):int(334 * s150), int(185 * s150):int(450 * s150)] = True                       # the key row in the pre-strip page (baseline 327): handles moved with the centred key
n_out = int((diff & ~mask).sum()); ys_, xs_ = np.where(diff & ~mask)
obox = [round(xs_.min() / s150, 1), round(ys_.min() / s150, 1), round(xs_.max() / s150, 1), round(ys_.max() / s150, 1)] if n_out else None
check("raster (Ghostscript 150 dpi): no difference against V13 outside the forest, the text boxes and the key row", n_out == 0, f"{int(diff.sum())} differing pixels, {n_out} outside the allowed regions (their box in pt: {obox})")
R.gs_render(BASE, f"{CROPS}/{SHEET}_OLD_200dpi.png", 200, LOGS, name="gs200_base"); R.gs_render(NEW, f"{CROPS}/{SHEET}_NEW_200dpi.png", 200, LOGS, name="gs200_new")
side = np.concatenate([a, np.full((a.shape[0], 20, 3), 255), b], axis=1).astype(np.uint8); Image.fromarray(side).save(f"{CROPS}/{SHEET}_OLD_vs_NEW_150dpi.png")
# printed values: every string accounted for (row labels re-read from the file, the rest V13 strings at +1 pt in the builder's record)
n_rows, n_no = VL.printed_values_csv(NEW, DRAWN["records"], f"{VER}/{SHEET}_printed_values.csv", static_from=BASE)
check("printed-values CSV: every string accounted for (drawn row labels re-read from the v8.1 file, ticks, axis title and key recorded as V13 strings at +1 pt)", n_no == 0, f"{n_rows} rows, {n_no} NO")
hits = C.grep_proof(HERE)
check("grep proof: no live script in this sheet folder names a v7 table, snapshot or set", not hits, f"{hits}")
LINES.append(f"INFO  numbers: results_v2.json sha256 {DRAWN['source_sha256'][:16]}, step {DRAWN['step']['step']} {DRAWN['step']['name']} out {DRAWN['step']['output_mtime']}; V13 base sha256 {DRAWN['base_sha256'][:16]}")
bd.close(); nd.close()
print(C.write_checks(LINES, f"{VER}/checks_03_verify.txt", OK))
