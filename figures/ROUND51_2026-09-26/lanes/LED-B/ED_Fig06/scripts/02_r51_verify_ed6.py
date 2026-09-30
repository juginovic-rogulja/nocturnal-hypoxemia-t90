#!$T90_PY
"""ED_Fig06 round 51 verify (child under wd_run.sh): the closed-up sheet against the V28 dump and the r51 record.
usage: 02_r51_verify_ed6.py <new.pdf> <v28_dump.json> <r51_record.json> <verify dir> <png_v28> <png_new>"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import fitz, json, os, sys, re, gc
import numpy as np, pandas as pd
from collections import Counter
NEW, DUMP, REC, VER, PNG_OLD, PNG_NEW = sys.argv[1:7]
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND37_2026-09-14/figures/lanes/V14_L6_HABITUAL_EXTERNAL/_common")
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND37_2026-09-14/figures/_common")
import lane_common as LC, v14lib as V
from PIL import Image
lines = []
def check(name, ok, detail=""): lines.append(f"{'PASS' if ok else 'FAIL'}  {name}" + (f": {detail}" if detail else "")); print(lines[-1][:400])
dump = json.load(open(DUMP)); rec = json.load(open(REC)); REG = rec["plot_region"]; import math
d = fitz.open(NEW); p = d[0]; W, H = p.rect.width, p.rect.height
check("page box equal to V28 (518.74 x 666.5)", abs(W - dump["page"][0]) < 0.01 and abs(H - dump["page"][1]) < 0.01, f"{W:.3f} x {H:.3f}")
new = []
for b in p.get_text("rawdict")["blocks"]:
    if b["type"] != 0: continue
    for l in b["lines"]:
        for s in l["spans"]:
            txt = "".join(ch["c"] for ch in s["chars"])
            if not txt.strip(): continue
            new.append(dict(text=txt, font=s["font"], size=round(s["size"], 3), color="#%06x" % s["color"], origin=[round(v, 3) for v in s["chars"][0]["origin"]], bbox=[round(v, 3) for v in s["bbox"]], dir=[round(v, 3) for v in l["dir"]]))
check("text layer clean: no NBSP and no soft hyphen", all("\xa0" not in s["text"] and "\xad" not in s["text"] for s in new), f"{len(new)} spans")
old = [dict(s, text=s["text"].replace("\xa0", " ").replace("\xad", "-")) for s in dump["spans"]]
rem = {(r["text"], round(r["x"], 2), round(r["y"], 2)) for r in rec["removed"]}; tl = {(t["text"], round(t["y"], 2)): t for t in rec["tick_labels"]}
used = [False] * len(new); miss = []; worst = 0.0
for s in old:
    key = (s["text"], round(s["origin"][0], 2), round(s["origin"][1], 2)); tkey = (s["text"], round(s["origin"][1], 2))
    if key in rem: continue
    ex = tl[tkey]["x"] if (tkey in tl and s["size"] == 11.0 and abs(s["origin"][1] - 626.876) < 0.5) else s["origin"][0]; hit = None
    for j, t in enumerate(new):
        if used[j] or t["text"] != s["text"] or abs(t["size"] - s["size"]) > 0.01 or (("Bold" in t["font"]) != ("Bold" in s["font"])) or t["color"] != s["color"] or t["dir"] != s["dir"]: continue
        if abs(t["origin"][1] - s["origin"][1]) < 0.05 and abs(t["origin"][0] - ex) < 0.05: hit = j; worst = max(worst, abs(t["origin"][0] - ex)); break
    if hit is None: miss.append((s["text"], s["origin"], "moved" if key in mvd else "kept"))
    else: used[hit] = True
check(f"every V28 span except the {len(rem)} removed is read back with the same text, size, weight, colour and baseline at its V28 x (0.05 pt), the four tick labels centred on the new tick x", not miss, f"largest offset {worst:.3f} pt, missing {miss[:6]}")
tl_ok = all(any(t["text"] == n["text"] and abs(n["origin"][1] - t["y"]) < 0.05 and abs((n["bbox"][0] + n["bbox"][2]) / 2 - t["centre"]) < 0.3 for n in new) for t in rec["tick_labels"])
check("the tick labels 0.5, 1, 2, 4 are centred on the new tick positions (0.3 pt)", tl_ok, f"{[(t['text'], t['centre']) for t in rec['tick_labels']]}")
extra = [(t["text"], t["origin"]) for j, t in enumerate(new) if not used[j]]
check("no span beyond those", not extra, f"{extra[:6]}")
check("the 26 removed strings are the 24 events/group-size strings and the two header lines, none of them left on the sheet", sorted(r["text"] for r in rec["removed"]) == sorted([r for r in (t["text"] for t in old) if re.fullmatch(r"\d+/[\d,]+", r)] + ["New diagnoses /", "group size"]) and not [t for t in new if re.fullmatch(r"\d+/[\d,]+", t["text"]) or t["text"] in ("New diagnoses /", "group size")], "")
fonts = sorted({f[3].split("+")[-1] for f in p.get_fonts(full=True)})
check("fonts: Arial faces only, no XObject, no image", all(f.replace(" ", "").lower().startswith("arial") for f in fonts) and len(p.get_xobjects()) == 0 and len(p.get_images()) == 0, f"{fonts}")
# drawings: outside the plot region every V28 item unchanged (the stale background skipped); inside: the forest read back against the file's values on the new axis
def style(it): return (tuple(round(v, 4) for v in it["fill"]) if it["fill"] else None, tuple(round(v, 4) for v in it["color"]) if it["color"] else None, round(it["width"] or 0, 3), str(it.get("dashes")), it["kinds"], it["n"])
def hx(c): return None if c is None else "#%02x%02x%02x" % tuple(int(round(v * 255)) for v in c[:3])
nd = []
for it in p.get_drawings():
    r = it["rect"]; nd.append(dict(rect=[r.x0, r.y0, r.x1, r.y1], fill=it.get("fill"), color=it.get("color"), width=it.get("width"), dashes=it.get("dashes"), kinds="".join(sorted(set(k[0] for k in it["items"]))), n=len(it["items"]), items=it["items"]))
def inside(r): return r[0] >= REG[0] and r[1] >= REG[1] and r[2] <= REG[2] and r[3] <= REG[3]                      # V28's plot region
REGN = [rec["box_new"][0] - 6.0, REG[1], REG[2], REG[3]]
def inside_new(r): return r[0] >= REGN[0] and r[1] >= REGN[1] and r[2] <= REGN[2] and r[3] <= REGN[3]                # the widened plot region on the new sheet
stale = [tuple(round(v, 3) for v in r) for r in rec["skipped_stale_layer"]]
od = [it for it in dump["drawings"] if tuple(round(v, 3) for v in it["rect"]) not in stale and not inside(it["rect"])]
n_in28 = sum(1 for it in dump["drawings"] if inside(it["rect"]))
taken = [False] * len(nd); unmatched = []; worst_r = 0.0
for it in od:
    hit = None
    for j, jt in enumerate(nd):
        if taken[j] or inside_new(jt["rect"]) or style(jt) != style(it): continue
        dr = max(abs(a - b) for a, b in zip(it["rect"], jt["rect"]))
        if dr <= 0.05: hit = j; worst_r = max(worst_r, dr); break
    if hit is None: unmatched.append((it["rect"], it["kinds"]))
    else: taken[hit] = True
extra_out = [jt["rect"] for j, jt in enumerate(nd) if not taken[j] and not inside_new(jt["rect"])]
check(f"line art outside the plot region: each of the {len(od)} V28 items (the 3 clipped stale-layer items skipped, declared) found once unchanged (same fill, stroke, width, dashes, path kinds, rect within 0.05 pt), nothing else outside", not unmatched and not extra_out, f"unmatched {len(unmatched)} {unmatched[:3]}, extra {len(extra_out)}, largest offset {worst_r:.3f}")
P = [jt for jt in nd if inside_new(jt["rect"])]; ax = rec["axis_new"]; a2, b2 = ax["a"], ax["b"]
def xn(v): return a2 + b2 * math.log(v)
boxes = [jt for jt in P if jt["fill"] and hx(jt["fill"]) == "#ffffff" and jt["color"] is None and jt["rect"][3] - jt["rect"][1] > 300]
rules = [jt for jt in P if jt["dashes"] not in (None, "[] 0") and jt["color"] and hx(jt["color"]) == "#1a1d21" and jt["rect"][3] - jt["rect"][1] > 300]
ticks = sorted([jt for jt in P if jt["color"] and hx(jt["color"]) == "#1a1d21" and 2.5 < jt["rect"][3] - jt["rect"][1] < 3.5 and jt["rect"][2] - jt["rect"][0] < 1.0 and jt["kinds"] == "l"], key=lambda jt: jt["rect"][0])
spines = [jt for jt in P if jt["kinds"] == "l" and jt["n"] == 1 and jt["rect"][3] - jt["rect"][1] < 0.01 and jt["rect"][2] - jt["rect"][0] > 50 and jt["color"] and hx(jt["color"]) == "#1a1d21" and abs((jt["width"] or 0) - 1.7) > 0.1]
cis = [jt for jt in P if jt["color"] and hx(jt["color"]) in ("#0288d1", "#8a9099") and abs((jt["width"] or 0) - 1.7) < 0.01 and jt["rect"][2] - jt["rect"][0] > 5]
mks = [jt for jt in P if jt["fill"] and hx(jt["fill"]) in ("#0288d1", "#8a9099") and jt["rect"][2] - jt["rect"][0] < 8 and jt["rect"][3] - jt["rect"][1] < 8]
check(f"inside the plot region: exactly one white plot box {[round(v, 2) for v in rec['box_new']]}, one dashed HR = 1 rule at x {xn(1.0):.2f}, one spine, 4 ticks at the new tick x, 24 intervals and 24 markers, nothing else ({n_in28} V28 items there, the overlay residue not reproduced)",
      len(boxes) == 1 and max(abs(a - b) for a, b in zip(boxes[0]["rect"], rec["box_new"])) < 0.05 and len(rules) == 1 and abs(rules[0]["rect"][0] - xn(1.0)) < 0.05 and len(spines) == 1 and len(ticks) == 4 and all(abs((t["rect"][0] + t["rect"][2]) / 2 - xn(v)) < 0.05 for t, v in zip(ticks, (0.5, 1, 2, 4))) and len(cis) == 24 and len(mks) == 24 and len(P) == 1 + 1 + 1 + 4 + 48,
      f"boxes {len(boxes)} rules {len(rules)} spines {len(spines)} ticks {len(ticks)} intervals {len(cis)} markers {len(mks)} items {len(P)}")
sp_ = spines[0] if spines else None; sw = rec["spine"]
check(f"the spine read back: one horizontal path at the axis baseline y {sw['y']:.2f} from x {sw['x_new'][0]:.2f} to {sw['x_new'][1]:.2f} (length {sw['x_new'][1] - sw['x_new'][0]:.2f} = the widened axis width; V28: {sw['x_v28'][0]:.2f} to {sw['x_v28'][1]:.2f}), {sw['width']:.1f} pt ink, the four ticks hanging from it (0.05 pt)",
      sp_ is not None and abs(sp_["rect"][1] - sw["y"]) < 0.05 and abs(sp_["rect"][0] - sw["x_new"][0]) < 0.05 and abs(sp_["rect"][2] - sw["x_new"][1]) < 0.05 and abs((sp_["rect"][2] - sp_["rect"][0]) - (rec["box_new"][2] - rec["box_new"][0])) < 0.05 and abs((sp_["width"] or 0) - sw["width"]) < 0.01 and all(abs(t["rect"][1] - sw["y"]) < 0.05 for t in ticks), f"{sp_['rect'] if sp_ else None}")
worst_m = 0.0; n_ok = 0
for r in rec["rows"]:
    m = [jt for jt in mks if abs((jt["rect"][1] + jt["rect"][3]) / 2 - r["y"]) < 0.1 and hx(jt["fill"]) == r["colour"]]; c = [jt for jt in cis if abs(jt["rect"][1] - r["y"]) < 0.1 and hx(jt["color"]) == r["colour"]]
    if len(m) == 1 and len(c) == 1:
        mr, cr = m[0]["rect"], c[0]["rect"]; worst_m = max(worst_m, abs((mr[0] + mr[2]) / 2 - xn(r["hr"])), abs(cr[0] - xn(r["lo"])), abs(cr[2] - xn(r["hi"])), abs((mr[2] - mr[0]) - r["size"])); n_ok += 1
check("the 24 intervals and markers read back at ln(lo), ln(hr), ln(hi) of results_crossed.csv on the new axis (0.05 pt), grey squares for preserved and blue circles for low oxygen with the V28 marker sizes", n_ok == 24 and worst_m < 0.05, f"{n_ok} rows, largest deviation {worst_m:.3f} pt")
check(f"the V28 marks sat at the file's values on the V28 axis (largest deviation {rec['v28_marks_vs_file_worst_pt']} pt): the same 24 hazard ratios and intervals are plotted", rec["v28_marks_vs_file_worst_pt"] < 0.05, "")
check(f"the widened axis: box x0 {rec['box_v28'][0]:.2f} -> {rec['box_new'][0]:.2f} (right edge {rec['box_new'][2]:.2f} and the value range {rec['axis_range'][0]:.3f} to {rec['axis_range'][1]:.3f} unchanged); the leftmost interval end at x {rec['leftmost_ink_new']:.1f} stays clear of the row labels (right edge {rec['row_label_right_edge']:.1f}) by 12 pt or more", rec["leftmost_ink_new"] - rec["row_label_right_edge"] >= 12.0 and abs(rec["box_new"][2] - rec["box_v28"][2]) < 0.01, "")
d.close(); gc.collect()
# raster proofs (Ghostscript 150 dpi renders made by the runner): identity outside the re-plotted region and the removed strings' boxes
A, B = (np.asarray(Image.open(x).convert("RGB")).astype(int) for x in (PNG_OLD, PNG_NEW)); S = 150 / 72
check("150 dpi renders of V28 and NEW have the same size", A.shape == B.shape, f"{A.shape}")
mask = np.zeros(A.shape[:2], bool); mask[int(REG[1] * S):int(REG[3] * S) + 1, int((rec["box_new"][0] - 6.0) * S):int(REG[2] * S) + 1] = True
for r in rec["removed"]:
    b0 = r["bbox"]; mask[max(0, int((b0[1] - 0.5) * S)):int((b0[3] + 0.5) * S) + 1, max(0, int((b0[0] - 0.5) * S)):int((b0[2] + 0.5) * S) + 1] = True
diff = (A != B).any(axis=2)
def boxes_of(m):
    ys, xs = np.nonzero(m)
    return [] if len(ys) == 0 else [(round(xs.min() / S, 1), round(ys.min() / S, 1), round(xs.max() / S, 1), round(ys.max() / S, 1))]
n_out = int((diff & ~mask).sum())
tmask = mask.copy()
for t in new:
    b0 = t["bbox"]; tmask[max(0, int((b0[1] - 0.5) * S)):int((b0[3] + 0.5) * S) + 1, max(0, int((b0[0] - 0.5) * S)):int((b0[2] + 0.5) * S) + 1] = True
n_out_t = int((diff & ~tmask).sum())
check("outside the re-plotted forest and the removed column, every differing pixel lies inside a glyph box of a re-written string (TextWriter packs strings of one colour into one text object and rounds their offsets to 1/1000 em: origins read back within 0.05 pt above); 0 differing pixels on any line art, mark, cell or margin", n_out_t == 0, f"{int(diff.sum())} differing pixels in all, {n_out} outside the forest and the removed column (all at glyph edges of re-written strings, box {boxes_of(diff & ~mask)}), {n_out_t} outside the glyph boxes too")
# values: the kept numbers equal the data files (round-37 gate) and V28
SRC = f"{LC.SV}/habitual_outcomes_v2"; MIN_MTIME = "2026-09-15 19:18:00"
sc_c = LC.sidecar(LC.hydrated(f"{SRC}/cells_crossed.csv"), min_mtime=MIN_MTIME); sc_r = LC.sidecar(LC.hydrated(f"{SRC}/results_crossed.csv"), min_mtime=MIN_MTIME)
cells = pd.read_csv(f"{SRC}/cells_crossed.csv"); crossed = pd.read_csv(f"{SRC}/results_crossed.csv")
A_ = cells[cells.design == "cells12"].set_index("cell"); counts = sorted(f"n = {int(n):,}" for n in A_.n)
OUTC = [("cvd", "Cardiovascular composite"), ("death", "Death from any cause"), ("htn2", "Hypertension"), ("diabetes", "Type 2 diabetes"), ("hf", "Heart failure"), ("dementia", "Dementia")]
CELLS4 = [("<5h", "normal T90<=1%"), ("<5h", "low T90>10%"), (">=7h", "normal T90<=1%"), (">=7h", "low T90>10%")]
cb = crossed[(crossed.model == "cells12") & (crossed["set"] == "primary")]; exp_hr = []
for k, name in OUTC:
    dd = cb[cb.outcome_key == k].set_index("contrast")
    for hb, ob in CELLS4:
        r = dd.loc[f"{hb} | {ob}"]; exp_hr.append(V.hr_ci(float(r.hr), float(r.lo95), float(r.hi95)))
got_counts = sorted(s["text"] for s in new if s["text"].startswith("n = ")); got_hr = [s["text"] for s in sorted(new, key=lambda s: s["origin"][1]) if re.fullmatch(r"\d\.\d\d \(\d\.\d\d-\d\.\d\d\)", s["text"])]
v28_hr = [s["text"] for s in sorted(old, key=lambda s: s["origin"][1]) if re.fullmatch(r"\d\.\d\d \(\d\.\d\d-\d\.\d\d\)", s["text"])]
check("every kept number identical to V28 and re-derived from the sidecar-gated data files: 12 group counts, 24 HR (95% CI) in row order (the events/group-size strings are the removed column)", got_counts == counts and got_hr == exp_hr == v28_hr, f"counts {got_counts == counts}, HR {got_hr == exp_hr}, HR = V28 {got_hr == v28_hr}")
alltext = " ".join(s["text"] for s in new)
check("house rules: no em dash, no semicolon on the sheet", "—" not in alltext and ";" not in alltext, "")
json.dump(dict(checks=lines, n_spans=len(new), fonts=fonts), open(f"{VER}/verify_record.json", "w"), indent=1, ensure_ascii=False)
open(f"{VER}/verify_lines.txt", "w").write("\n".join(lines) + "\n")
ok = all(l.startswith("PASS") for l in lines); print("VERIFY", "ALL PASS" if ok else "FAIL"); sys.exit(0 if ok else 1)
