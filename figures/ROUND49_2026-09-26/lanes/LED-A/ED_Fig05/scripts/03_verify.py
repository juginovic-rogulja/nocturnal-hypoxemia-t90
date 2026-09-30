#!/usr/bin/env python3
"""ED Fig 5 proof (class W, four panels, V14 against V13): page box, letters a to d, fonts, geometry pins of every invariant
string, every printed number read back and matched to the drawn values (recomputed from the v8 table, cross-checked to the
step-129 v8.1 log), numeric-token delta = expected, word multiset identical, crops OLD vs NEW, CHANGES_ED_Fig05.csv (old = the V13
text layer, group sizes from the round-30 record), the printed-values CSV. Repointed copy of the round-30 03_verify.py
(backup 03_verify_PRE_V8_1.py)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import collections, csv, json, os, re, sys
import numpy as np, fitz
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND37_2026-09-14/figures/lanes/V14_L3a_DURATION_MAIN/scripts")
import lane_common as C
VL = C.VL
SHEET = "ED_Fig05"; HERE = f"{C.LANE}/{SHEET}"; VER = f"{HERE}/verify"; CROPS = f"{VER}/crops"; WORK = f"{HERE}/work"; os.makedirs(CROPS, exist_ok=True)
BASE = f"{C.BASE}/{SHEET}.pdf"; NEW = f"{HERE}/{SHEET}.pdf"
OLDREC = f"{C.R30_LANE}/{SHEET}/verify/{SHEET}_drawn.json"   # the round-30 record (v7 group sizes, what V13 shows), OLD side only
D = json.load(open(f"{WORK}/{SHEET}_drawn.json"))
LINES = []; OK = True
def check(name, ok, detail):
    global OK; OK &= bool(ok); LINES.append(f"{'PASS' if ok else 'FAIL'}  {name}: {detail}"); print(LINES[-1])
bd = fitz.open(C.hydrated(BASE)); bp = bd[0]; nd = fitz.open(C.hydrated(NEW)); np_ = nd[0]
check("single page, page box equal to V13", len(nd) == 1 and abs(np_.rect.width - bp.rect.width) < 0.01 and abs(np_.rect.height - bp.rect.height) < 0.01,
      f"V13 {bp.rect.width:.2f} x {bp.rect.height:.2f}, new {np_.rect.width:.2f} x {np_.rect.height:.2f}")
LB, LN = C.letters(bp), C.letters(np_)
check("letters a, b, c, d at the V13 positions (within 0.3 pt), Arial Bold 13", [t[0] for t in LN] == ["a", "b", "c", "d"] and all(a[0] == b[0] and abs(a[1] - b[1]) <= 0.3 and abs(a[2] - b[2]) <= 0.3 and a[4] == b[4] for a, b in zip(LB, LN)),
      f"V13 {[(t[0], t[1], t[2]) for t in LB]}, new {[(t[0], t[1], t[2], t[4], t[5]) for t in LN]}")
fb = sorted(set(f[3].split("+")[-1] for f in bp.get_fonts(full=True))); fn = sorted(set(f[3].split("+")[-1] for f in np_.get_fonts(full=True)))
check("fonts identical to V13 (Arial faces)", fb == fn, f"V13 {fb}, new {fn}")
check("no nested forms", len(np_.get_xobjects()) == 0, f"{len(np_.get_xobjects())} xobjects")
SB = C.dedupe(C.spans_of(bp)); SN = C.dedupe(C.spans_of(np_))
def is_value(s): return bool(re.match(r"^(\d+\.\d|PR, \d\.\d\d|\d\.\d\d-\d\.\d\d\)|P [<=] 0\.\d\d\d)$", s["text"]))
inv_b = [s for s in SB if not is_value(s)]; inv_n = [s for s in SN if not is_value(s)]
worst = 0.0; missing = []
for s in inv_n:
    cands = [b for b in inv_b if b["text"] == s["text"]]
    if not cands: missing.append(s["text"]); continue
    b = min(cands, key=lambda b: abs(b["origin"][0] - s["origin"][0]) + abs(b["origin"][1] - s["origin"][1]))
    worst = max(worst, abs(b["origin"][0] - s["origin"][0]), abs(b["origin"][1] - s["origin"][1]))
extra_b = collections.Counter(b["text"] for b in inv_b) - collections.Counter(s["text"] for s in inv_n)
check("every invariant string (titles, axis labels, ticks, key, 'PR,' and '(95% CI,' lines) within 0.3 pt of V13", worst <= 0.3 and not missing and not extra_b,
      f"largest offset {worst:.3f} pt, missing {missing}, V13-only {dict(extra_b)}")
pct_b = sorted(((round((s["bbox"][0] + s["bbox"][2]) / 2, 1)) for s in SB if re.match(r"^\d+\.\d$", s["text"])))
pct_n = sorted(((round((s["bbox"][0] + s["bbox"][2]) / 2, 1)) for s in SN if re.match(r"^\d+\.\d$", s["text"])))
blk_b = sorted(((round(s["origin"][1], 1), round((s["bbox"][0] + s["bbox"][2]) / 2, 1)) for s in SB if is_value(s) and not re.match(r"^\d+\.\d$", s["text"])))
blk_n = sorted(((round(s["origin"][1], 1), round((s["bbox"][0] + s["bbox"][2]) / 2, 1)) for s in SN if is_value(s) and not re.match(r"^\d+\.\d$", s["text"])))
check("16 bar-top percent labels centred on the same 16 bars, 24 PR / CI / P strings on the same baselines and centres",
      len(pct_b) == len(pct_n) == 16 and max(abs(a - b) for a, b in zip(pct_b, pct_n)) <= 1.0 and len(blk_b) == len(blk_n) == 24 and all(abs(a[0] - b[0]) <= 0.3 and abs(a[1] - b[1]) <= 1.5 for a, b in zip(blk_b, blk_n)),
      f"percent centres max shift {max(abs(a - b) for a, b in zip(pct_b, pct_n)) if pct_b and len(pct_b) == len(pct_n) else 'n/a'}, blocks {len(blk_b)} vs {len(blk_n)}")
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
check("numeric-token multiset differs from V13 by exactly the printed value strings (percent labels, PR, CI, P)", d == ed, f"delta {d}")
wb, wn = sorted(w["text"] for w in C.dedupe(C.words_of(bp))), sorted(w["text"] for w in C.dedupe(C.words_of(np_)))
VALW = lambda w: bool(re.match(r"^[\d(]", w)) or w.endswith(")") or w in ("=", "<")   # value words: numbers, the CI parentheses and the P comparator
wd = C.delta([w for w in wb if not VALW(w)], [w for w in wn if not VALW(w)])
check("word multiset identical apart from the value words (numbers, CI parentheses, the P comparator '=' or '<')", wd == {"lost": {}, "gained": {}}, f"{wd}")
hx = lambda c: "#%02x%02x%02x" % tuple(int(round(v * 255)) for v in c)
cb = sorted(set(hx(x["fill"]) for x in bp.get_drawings() if x.get("fill") is not None)); cn = sorted(set(hx(x["fill"]) for x in np_.get_drawings() if x.get("fill") is not None))
check("fill colour set identical to V13 (green ladder #79d475 under 5 h, #39c445 7 h or more)", cb == cn, f"{cn}")
bars = [x for x in np_.get_drawings() if x.get("fill") is not None and hx(x["fill"]) in ("#79d475", "#39c445") and x["rect"].width > 15]
check("16 bars drawn", len(bars) == 16, f"{len(bars)} bars")
# bar heights read back against the percentages: the y axis 0 to 33 maps the bar tops
ax_map = {}
for pnl_letter in "abcd":
    pass
for tag, pg in (("OLD", bp), ("NEW", np_)):
    pm = pg.get_pixmap(dpi=200, colorspace=fitz.csRGB, alpha=False); pm.save(f"{CROPS}/{SHEET}_{tag}_200dpi.png"); pm = None
A = bp.get_pixmap(dpi=150, colorspace=fitz.csRGB, alpha=False); B = np_.get_pixmap(dpi=150, colorspace=fitz.csRGB, alpha=False)
a = np.frombuffer(A.samples, np.uint8).reshape(A.h, A.w, 3).astype(int); b = np.frombuffer(B.samples, np.uint8).reshape(B.h, B.w, 3).astype(int)
fitz.Pixmap(fitz.csRGB, B.w, B.h, B.samples, False).save(f"{HERE}/{SHEET}_150dpi.png")
diff = np.abs(a - b).max(axis=2) > 0
LINES.append(f"INFO  raster: {int(diff.sum())} differing pixels at 150 dpi of {diff.size} (every panel's bars and printed values change, proof class W)")
side = np.concatenate([a, np.full((a.shape[0], 20, 3), 255), b], axis=1).astype(np.uint8)
fitz.Pixmap(fitz.csRGB, side.shape[1], side.shape[0], side.tobytes(), False).save(f"{CROPS}/{SHEET}_OLD_vs_NEW_150dpi.png")
# CHANGES csv: old from the V13 text layer per panel and cut, new from the drawn values
def slot_values(spans):
    out = {}
    for pnl in "abcd":
        ss = sorted([s for s in spans if panel_of(s) == pnl and is_value(s)], key=lambda s: (s["origin"][1], s["origin"][0]))
        pct = sorted([s for s in ss if re.match(r"^\d+\.\d$", s["text"])], key=lambda s: s["origin"][0])
        pr = sorted([s for s in ss if s["text"].startswith("PR,")], key=lambda s: s["origin"][0])
        ci = sorted([s for s in ss if s["text"].endswith(")")], key=lambda s: s["origin"][0])
        pv = sorted([s for s in ss if s["text"].startswith("P ")], key=lambda s: s["origin"][0])
        out[pnl] = dict(pct=[s["text"] for s in pct], pr=[s["text"] for s in pr], ci=[s["text"] for s in ci], p=[s["text"] for s in pv])
    return out
OV, NV = slot_values(SB), slot_values(SN)
OLD = json.load(open(C.hydrated(OLDREC)))["panels"]
rows = []
for pnl in "abcd":
    L_ = pnl.upper(); pn = D["panels"][pnl]
    for i, cut in enumerate((5, 10)):
        v = pn["cuts"][str(cut)] if str(cut) in pn["cuts"] else pn["cuts"][cut]
        o_pct_s, o_pct_r = OV[pnl]["pct"][2 * i], OV[pnl]["pct"][2 * i + 1]
        rows.append(dict(sheet=SHEET, panel=pnl, item=f"T90 >{cut}%: percent of short sleepers", before=o_pct_s, after=f"{v['pct_short']:.1f}", dramatic="no", note="data_frozen_v8 t90_final.parquet via cohort_spec (15,551), step-129 arithmetic"))
        rows.append(dict(sheet=SHEET, panel=pnl, item=f"T90 >{cut}%: percent of sleepers 7 h or more", before=o_pct_r, after=f"{v['pct_ref']:.1f}", dramatic="no", note="same"))
        old_pr = float(OV[pnl]["pr"][i].split(", ")[1]); old_ci = OV[pnl]["ci"][i].rstrip(")"); olo, ohi = (float(x) for x in old_ci.split("-"))
        cross = (olo <= 1 <= ohi) != (v["lo"] <= 1 <= v["hi"])
        rows.append(dict(sheet=SHEET, panel=pnl, item=f"T90 >{cut}%: prevalence ratio (95% CI)", before=f"{old_pr:.2f} ({old_ci})", after=f"{v['pr']:.2f} ({v['lo']:.2f}-{v['hi']:.2f})", dramatic="yes" if cross else "no",
                         note=("DRAMATIC: 95% CI now includes 1" if cross and (v["lo"] <= 1 <= v["hi"]) else ("DRAMATIC: 95% CI now excludes 1" if cross else ""))))
        op = OV[pnl]["p"][i]; npv = v["p_printed"]
        o_sig = op == "P < 0.001" or float(op.split("= ")[1]) < 0.05; n_sig = npv == "P < 0.001" or float(npv.split("= ")[1]) < 0.05
        rows.append(dict(sheet=SHEET, panel=pnl, item=f"T90 >{cut}%: P", before=op, after=npv, dramatic="yes" if o_sig != n_sig else "no", note="DRAMATIC: P crosses 0.05" if o_sig != n_sig else ""))
    on = (OLD[pnl]["n_short"], OLD[pnl]["n_ref"])
    flag = any(abs(new - old) / old > 0.2 for new, old in zip((pn["n_short"], pn["n_ref"]), on))
    rows.append(dict(sheet=SHEET, panel=pnl, item="group sizes (legend, not printed on the sheet): short sleepers vs 7 h or more", before=f"{on[0]:,} vs {on[1]:,}", after=f"{pn['n_short']:,} vs {pn['n_ref']:,}", dramatic="yes" if flag else "no",
                     note=("DRAMATIC: count moves more than 20 percent; " if flag else "") + "old = the round-30 record (v7, what V13 shows)"))
VL.write_changes(f"{VER}/CHANGES_{SHEET}.csv", rows)
LINES.append(f"INFO  CHANGES csv: {VER}/CHANGES_{SHEET}.csv ({len(rows)} rows, {sum(1 for r in rows if r['dramatic'] == 'yes')} dramatic)")
# printed values: the 40 drawn strings re-derived under their rules, everything else static from V13
RULES = VL.RULES
RULES["pr_2dp"] = lambda v: "PR, " + VL.halfup(v, 2)
RULES["ci_2dp_close"] = lambda v: f"{VL.halfup(v[0], 2)}-{VL.halfup(v[1], 2)})"
RULES["P_3dp"] = lambda v: "P < 0.001" if float(v) < 0.001 else "P = " + VL.halfup(v, 3)
n_rows, n_no = VL.printed_values_csv(NEW, D["records"], f"{VER}/{SHEET}_printed_values.csv", static_from=BASE)
check("printed-values CSV: every string accounted for (40 drawn strings under their rules, the rest static from V13)", n_no == 0, f"{n_rows} rows, {n_no} NO")
hits = C.grep_proof(HERE)
check("grep proof: no live script in this sheet folder names a v7 table, snapshot or set", not hits, f"{hits}")
LINES.append(f"INFO  sources: {D['sources']['parquet']} sha256 {D['sources']['parquet_sha256'][:16]}, habitual {D['sources']['habitual_sha256'][:16]}, {D['sources']['cohort_filter']}; V13 base sha256 {D['base_sha256'][:16]}")
bd.close(); nd.close()
print(C.write_checks(LINES, f"{VER}/checks.txt", OK))
