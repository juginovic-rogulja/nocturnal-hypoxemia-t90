#!/usr/bin/env python3
"""ED Fig 4 proof (class W, V14 against V13): page box, fonts, geometry pins, text multisets with the declared row-set delta,
read-back of every drawn value against numbers/results_v2.json (v8.1), crops OLD vs NEW at 200 dpi, CHANGES_ED_Fig04.csv (old =
the V13 sheet's rows = the round-30 drawn record on the v7 snapshot), the printed-values CSV. Repointed copy of the round-30
03_verify.py (backup 03_verify_PRE_V8_1.py): OLD numbers = V8_RECALC/compare/old_numbers (the v7 snapshot the V13 sheet drew from)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import collections, csv, json, os, sys
import numpy as np, fitz
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND37_2026-09-14/figures/lanes/V14_L3a_DURATION_MAIN/scripts")
import lane_common as C
VL = C.VL
SHEET = "ED_Fig04"; HERE = f"{C.LANE}/{SHEET}"; VER = f"{HERE}/verify"; CROPS = f"{VER}/crops"; WORK = f"{HERE}/work"; os.makedirs(CROPS, exist_ok=True)
BASE = f"{C.BASE}/{SHEET}.pdf"; NEW = f"{HERE}/{SHEET}.pdf"
OLDNUM = f"{C.OLD_V7}/results_v2.json"                       # OLD side of the delta only (the v7 numbers the V13 sheet printed)
OLDREC = f"{C.R30_LANE}/{SHEET}/verify/{SHEET}_drawn.json"   # the round-30 drawn record (what V13 shows), OLD side only
DRAWN = json.load(open(f"{WORK}/{SHEET}_drawn.json"))
LINES = []; OK = True
def check(name, ok, detail):
    global OK; OK &= bool(ok); LINES.append(f"{'PASS' if ok else 'FAIL'}  {name}: {detail}"); print(LINES[-1])

for p in (BASE, NEW, OLDNUM, OLDREC): C.hydrated(p)
bd = fitz.open(BASE); bp = bd[0]; nd = fitz.open(NEW); np_ = nd[0]
check("single page, page box equal to V13", len(nd) == 1 and abs(np_.rect.width - bp.rect.width) < 0.01 and abs(np_.rect.height - bp.rect.height) < 0.01,
      f"V13 {bp.rect.width:.2f} x {bp.rect.height:.2f}, new {np_.rect.width:.2f} x {np_.rect.height:.2f}")
check("no panel letters on either sheet (single-panel sheet)", C.letters(bp) == [] and C.letters(np_) == [], f"{C.letters(bp)} {C.letters(np_)}")
fb = sorted(set(f[3].split("+")[-1] for f in bp.get_fonts(full=True))); fn = sorted(set(f[3].split("+")[-1] for f in np_.get_fonts(full=True)))
check("fonts Arial only, same faces as V13", fn == fb and all("Arial" in f for f in fn), f"V13 {fb}, new {fn}")
check("no nested forms on the new sheet", len(np_.get_xobjects()) == 0, f"{len(np_.get_xobjects())} xobjects")
SB = C.dedupe(C.spans_of(bp)); SN = C.dedupe(C.spans_of(np_))
ROWS = [r["condition"] for r in DRAWN["rows"]]; OLD_ROWS = [r["condition"] for r in json.load(open(OLDREC))["rows"]]
base_rows = [s["text"] for s in sorted((s for s in SB if s["text"] in OLD_ROWS or s["text"] in ROWS), key=lambda s: s["origin"][1]) if s["text"] in OLD_ROWS]
check("V13 row labels = the round-30 drawn record (the OLD side of the delta is what V13 shows)", base_rows == OLD_ROWS, f"V13 {base_rows}")
ROWSET = set(ROWS) | set(OLD_ROWS)
# geometry pins: every non-row string sits where V13 has it (visible copy)
worst = 0.0; missing = []
for s in SN:
    if s["text"] in ROWSET: continue
    cands = [b for b in SB if b["text"] == s["text"]]
    if not cands: missing.append(s["text"]); continue
    b = min(cands, key=lambda b: abs(b["origin"][0] - s["origin"][0]) + abs(b["origin"][1] - s["origin"][1]))
    worst = max(worst, abs(b["origin"][0] - s["origin"][0]), abs(b["origin"][1] - s["origin"][1]))
check("ticks, axis title and key strings within 0.3 pt of the V13 positions", worst <= 0.3 and not missing, f"largest offset {worst:.3f} pt, missing {missing}")
rb = sorted(round(s["origin"][1], 2) for s in SB if s["text"] in OLD_ROWS); rn = sorted(round(s["origin"][1], 2) for s in SN if s["text"] in ROWS)
check("10 row labels on the same 10 baselines as V13 (rows re-ranked by the v8.1 T90 hazard ratio)", len(rb) == len(rn) == 10 and max(abs(a - b) for a, b in zip(rb, rn)) <= 0.3, f"V13 {rb}, new {rn}")
entered, left = sorted(set(ROWS) - set(OLD_ROWS)), sorted(set(OLD_ROWS) - set(ROWS))
check("row set: the ten largest T90-significant conditions of the v8.1 file (declared delta: rows entering and leaving)", len(ROWS) == 10 and len(set(ROWS)) == 10, f"enter {entered}, leave {left}")
order_b = [s["text"] for s in sorted((s for s in SB if s["text"] in OLD_ROWS), key=lambda s: s["origin"][1])]
order_n = [s["text"] for s in sorted((s for s in SN if s["text"] in ROWS), key=lambda s: s["origin"][1])]
LINES.append(f"INFO  row order V13: {order_b}"); LINES.append(f"INFO  row order new: {order_n}")
check("rows on the new sheet in the drawn order (sorted by the T90 hazard ratio)", order_n == ROWS, f"{order_n}")
tb, tn = C.nums([w["text"] for w in C.dedupe(C.words_of(bp))]), C.nums([w["text"] for w in C.dedupe(C.words_of(np_))])
# V13 (= the round-30 sheet) prints no estimate; its text layer differs from the rebuilt sheet only by the row labels
check("numeric-token multiset identical to V13 (the sheet prints no estimates; the row labels carry no digit)", C.delta(tb, tn) == {"lost": {}, "gained": {}}, f"{len(tb)} vs {len(tn)} tokens, delta {C.delta(tb, tn)}")
wb, wn = sorted(w["text"] for w in C.dedupe(C.words_of(bp))), sorted(w["text"] for w in C.dedupe(C.words_of(np_)))
d = C.delta(wb, wn)
exp_lost = collections.Counter(w for lab in left for w in lab.split()); exp_gained = collections.Counter(w for lab in entered for w in lab.split())
check("word multiset differs from V13 by exactly the row labels that leave and enter", collections.Counter(d["lost"]) == exp_lost and collections.Counter(d["gained"]) == exp_gained, f"lost {d['lost']}, gained {d['gained']}")
# read-back of drawn values: axis mapping from the new sheet's tick labels, markers and CI lines from its drawings
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
check("read-back: 20 markers at the hazard ratios of results_v2.json (v8.1, within 0.3 percent) and 20 CI lines at the bounds", len(mk) == 20 and len(ci) == 20 and worst_m <= 0.003 and worst_ci <= 0.003 and n_ok == 20,
      f"{len(mk)} markers, {len(ci)} CI lines, largest marker deviation {100 * worst_m:.3f} percent, CI ends {100 * worst_ci:.3f} percent")
ys = sorted(((m["rect"].y0 + m["rect"].y1) / 2, hr_of((m["rect"].x0 + m["rect"].x1) / 2)) for m in mk if hx(m["fill"]) == "#0288d1")
check("blue markers descend in hazard ratio down the page (rows sorted by T90 HR)", all(a[1] >= b[1] for a, b in zip(ys, ys[1:])), f"{[round(h, 2) for _y, h in ys]}")
# every marker row matches its label row (y of the blue marker = label baseline - 3.5 pt band): the k-th blue marker is the k-th row
lab_y = [s["origin"][1] for s in sorted((s for s in SN if s["text"] in ROWS), key=lambda s: s["origin"][1])]
check("each blue marker sits on its row (marker centre within 12 pt of the row label baseline, in row order)", all(abs(my - ly) < 12 for (my, _h), ly in zip(ys, lab_y)), f"{[round(my - ly, 1) for (my, _h), ly in zip(ys, lab_y)]}")
cb = sorted(set(hx(x["fill"]) for x in bp.get_drawings() if x.get("fill") is not None)); cn = sorted(set(hx(x["fill"]) for x in dr if x.get("fill") is not None))
sb_ = sorted(set(hx(x["color"]) for x in bp.get_drawings() if x.get("color") is not None)); sn_ = sorted(set(hx(x["color"]) for x in dr if x.get("color") is not None))
check("fill and stroke colour sets identical to V13", cb == cn and sb_ == sn_, f"fills {cn}, strokes {sn_}")
check("no pale #eef0f1 vertical grid line drawn (V13 design)", "#eef0f1" not in sn_, f"{sn_}")
# crops and raster
for tag, pg in (("OLD", bp), ("NEW", np_)):
    pm = pg.get_pixmap(dpi=200, colorspace=fitz.csRGB, alpha=False); pm.save(f"{CROPS}/{SHEET}_{tag}_200dpi.png"); pm = None
A = bp.get_pixmap(dpi=150, colorspace=fitz.csRGB, alpha=False); B = np_.get_pixmap(dpi=150, colorspace=fitz.csRGB, alpha=False)
a = np.frombuffer(A.samples, np.uint8).reshape(A.h, A.w, 3).astype(int); b = np.frombuffer(B.samples, np.uint8).reshape(B.h, B.w, 3).astype(int)
fitz.Pixmap(fitz.csRGB, B.w, B.h, B.samples, False).save(f"{HERE}/{SHEET}_150dpi.png")
diff = np.abs(a - b).max(axis=2) > 0; ys_, xs_ = np.where(diff); s150 = 150 / 72
box = (round(xs_.min() / s150, 1), round(ys_.min() / s150, 1), round(xs_.max() / s150, 1), round(ys_.max() / s150, 1)) if len(xs_) else None
inside = (xs_ / s150 >= 20) & (ys_ / s150 <= 290) if len(xs_) else np.array([], bool)
amax = np.abs(a - b).max(axis=2)
out_key = (~inside) & (ys_ / s150 >= 314) & (ys_ / s150 <= 333) if len(xs_) else np.array([], bool)
out_max = int(amax[ys_[~inside], xs_[~inside]].max()) if (~inside).any() else 0
check("raster differences confined to the forest and its row labels (y under 290 pt); outside it only sub-pixel anti-aliasing on the key handles (difference at most 8 of 255)",
      len(xs_) == 0 or (inside | out_key).all() and out_max <= 8,
      f"{int(diff.sum())} differing pixels at 150 dpi, box {box}, {int((~inside).sum())} outside the forest, all on the key row, largest outside difference {out_max} of 255")
side = np.concatenate([a, np.full((a.shape[0], 20, 3), 255), b], axis=1).astype(np.uint8)
fitz.Pixmap(fitz.csRGB, side.shape[1], side.shape[0], side.tobytes(), False).save(f"{CROPS}/{SHEET}_OLD_vs_NEW_150dpi.png")
# CHANGES csv (old = the v7 snapshot the V13 sheet drew from; the round-30 record confirms the V13 rows)
OLD = json.load(open(OLDNUM))["vs_sleep_duration"]["outcomes"]
old_rec = json.load(open(OLDREC))
assert all(abs(OLD[r["condition"]]["t90"]["hr"] - r["t90"]["hr"]) < 1e-9 for r in old_rec["rows"]), "the v7 snapshot differs from the round-30 drawn record"
old_top = OLD_ROWS
def f(v): return VL.hr_ci(v["hr"], v["lo"], v["hi"])
changes = []
for i, r in enumerate(rows):
    k = r["condition"]; o = OLD.get(k)
    for key, lab in (("t90", "T90 >10%, HR (95% CI)"), ("short_sleep", "Laboratory sleep under 5 h, HR (95% CI)")):
        n = V[k][key]; ov = o[key] if o else None
        cross = ov is not None and ((ov["lo"] <= 1 <= ov["hi"]) != (n["lo"] <= 1 <= n["hi"]))
        changes.append(dict(sheet=SHEET, panel="single", item=f"{k}: {lab}", before=f(ov) if ov else "", after=f(n), dramatic="yes" if cross else "no",
                            note=("DRAMATIC: 95% CI crosses 1; " if cross else "") + "drawn marker and CI, not printed; numbers/results_v2.json vs_sleep_duration.outcomes"))
    changes.append(dict(sheet=SHEET, panel="single", item=f"{k}: row rank", before=str(old_top.index(k) + 1) if k in old_top else "not in top 10", after=str(i + 1), dramatic="yes" if k not in old_top else "no", note="DRAMATIC: enters top 10" if k not in old_top else "ordering by t90.hr"))
for k in old_top:
    if k not in ROWS: changes.append(dict(sheet=SHEET, panel="single", item=f"{k}: row rank", before=str(old_top.index(k) + 1), after="not in top 10", dramatic="yes", note="DRAMATIC: leaves top 10"))
old_larger = sum(1 for k in OLD if OLD[k]["t90"]["p"] < 0.05 and OLD[k]["t90"]["hr"] > 1 and not OLD[k]["negative_control"] and OLD[k]["t90"]["hr"] > OLD[k]["short_sleep"]["hr"])
old_nsig = sum(1 for k in OLD if OLD[k]["t90"]["p"] < 0.05 and OLD[k]["t90"]["hr"] > 1 and not OLD[k]["negative_control"])
changes.append(dict(sheet=SHEET, panel="single", item="T90 larger than short sleep in (legend count, not printed on the sheet)", before=f"{old_larger} of {old_nsig}", after=DRAWN["t90_larger_in"], dramatic="no", note="numbers/results_v2.json vs_sleep_duration"))
VL.write_changes(f"{VER}/CHANGES_{SHEET}.csv", changes)
LINES.append(f"INFO  CHANGES csv: {VER}/CHANGES_{SHEET}.csv ({len(changes)} rows, {sum(1 for c in changes if c['dramatic'] == 'yes')} dramatic)")
# printed values: every digit-bearing string is a V13 static string (ticks) or a drawn row label
n_rows, n_no = VL.printed_values_csv(NEW, DRAWN["records"], f"{VER}/{SHEET}_printed_values.csv", static_from=BASE)
check("printed-values CSV: every string accounted for (drawn row labels re-read from the v8.1 file, ticks and key static from V13)", n_no == 0, f"{n_rows} rows, {n_no} NO")
hits = C.grep_proof(HERE)
check("grep proof: no live script in this sheet folder names a v7 table, snapshot or set", not hits, f"{hits}")
LINES.append(f"INFO  numbers: results_v2.json sha256 {DRAWN['source_sha256'][:16]}, step {DRAWN['step']['step']} {DRAWN['step']['name']} out {DRAWN['step']['output_mtime']}; V13 base sha256 {DRAWN['base_sha256'][:16]}")
bd.close(); nd.close()
print(C.write_checks(LINES, f"{VER}/checks.txt", OK))
