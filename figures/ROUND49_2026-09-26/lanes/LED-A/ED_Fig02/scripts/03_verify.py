#!/usr/bin/env python3
"""ROUND 49 lane LED-A copy (2026-09-26): verifies the composed sheet before the title strip (work/ED_Fig02_composed.pdf, the V13 page box) with the
declared deltas of this round: every size +1 pt (the stretched key and caption strings Tf + 1 at the base x-scale, which PyMuPDF reports as 9.454 and
11.038), one row label kept at 10 pt, the two printed notes absent (round 40), the key-column raster window narrowed to the markers, the v7 CHANGES
block dropped (this round changes no number). Original: ROUND38_2026-09-16/figures/lanes/ED_Fig02/scripts/03_verify.py.
ED_Fig02, V14 lane L2 (round 37), step 3: Proof W of ED_Fig02.pdf against the V13 design baseline and the numbers files.
Rewritten from the round-30 verifier (byte copy beside as 03_verify_PRE_V8_1.py) for the V14 deltas: the row set (the generator's
rule on the v8.1 lag ladder and causal tests: chronic kidney disease enters, cellulitis leaves; the five v8.1 controls replace
back pain and cataract), panel c drawn from the group P v7 file with NE for the two controls it does not hold, the shared x range
extended to the left (inguinal hernia, panel b, 0.696), and every data string. Writes verify/checks.txt, values_agreement.txt,
expected_delta.json, CHANGES_ED_Fig02.csv, ED_Fig02_printed_values.csv, crops/ (from 02_compose.py)."""
import csv, gc, json, math, os, re, sys, time, collections
import numpy as np, pandas as pd, fitz
from PIL import Image
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ed2_common import *
from v14lib import printed_values_csv, Checks, write_changes, q_text

NEW = os.path.join(WORK, "ED_Fig02_composed.pdf"); OLDP = os.path.join(WORK, "ED_Fig02_OLDREPLOT.pdf")   # round 49: the composed sheet (the strip is stamped afterwards)
G = load_geometry(); W, H = G["page"]; BV = json.load(open(hydrated(os.path.join(WORK, "base_visible.json"))))
DRW = json.load(open(hydrated(os.path.join(VER, "ED_Fig02_drawn.json")))); DRO = json.load(open(hydrated(os.path.join(VER, "ED_Fig02_drawn_old.json"))))
DN = json.load(open(hydrated(os.path.join(WORK, "data_new.json")))); DO = json.load(open(hydrated(os.path.join(WORK, "data_old.json"))))
TOL = 0.3
DELTA = float(DRW["meta"]["delta_pt"]); KEEP = set(DRW["meta"]["keep_size"]); NOTES = list(DRW["meta"]["notes_removed"]); assert DELTA == 1.0
STRETCHED = {8.425: 9.454, 10.043: 11.038}                        # PyMuPDF size of the anisotropic strings at Tf + 1 (geometric mean of Tf and the x-scaled Tf)
def new_size(b, text=None): return b if text in KEEP else STRETCHED.get(round(b, 3), round(b + DELTA, 3))
ck = Checks(f"V14 lane L2, ED_Fig02, {time.strftime('%Y-%m-%d %H:%M')}. Base = V13 {BASE} (sha256 {sha256(BASE)[:16]}). Output {NEW} (sha256 {sha256(NEW)[:16]}). "
            f"Proof W (round v8.2): the whole sheet re-plotted; declared: row set (in {DRW['meta']['rowset_delta']['entering']}, out {DRW['meta']['rowset_delta']['leaving']}), panel c from the v8.2 "
            f"sleep_unadjusted_v1.json (step 125, all five controls fitted, c_missing {DRW['meta']['c_missing']}), x range {DRW['meta']['xlim_v13']} -> {DRW['meta']['xlim']} (same ticks).")
def check(ok, msg): return ck.log(msg, bool(ok))
def info(msg): ck.lines.append("INFO: " + msg); print("INFO: " + msg[:300], flush=True)

# ------------------------------------------------------------------ data and provenance
fresh = derive(NUM, "new")
check(fresh["rows"] == DN["rows"] and fresh["sha256"] == DN["sha256"], "the data step's cached work/data_new.json equals a fresh derivation from the numbers files (rows, values, sha256)")
for f in V8_FILES:
    sc = DN["sidecars"][f]
    check(sc["sha256"] == DN["sha256"][f] == DRW["meta"]["sha256"][f] and sc["v8_inputs"] and not sc["group_P"], f"{f}: v8 sidecar (step {sc['step']} {sc['name']}, out {sc['output_mtime']}), sha256 of the file read now = the builder's record = the sidecar's")
for f in V8_2_FILES:
    sc = DN["sidecars"][f]
    check(sc.get("output_mtime", "") >= V8_2_CUT, f"{f}: a v8.2 output (mtime {sc.get('output_mtime')} after {V8_2_CUT}); nothing kept from the v7 group P copy")
check(DN["c_missing"] == [] and DRW["meta"]["c_missing"] == [], "panel c: every drawn row, the five controls included, has a sleep-period estimate in the v8.2 file (no NE)")
check(DRW["meta"]["base_sha256"] == sha256(BASE), "the builder drew against the V13 sheet on disk now")
check(DN["conds"] + DN["controls"] == DRW["meta"]["rows"] == [r["condition"] for r in DN["rows"]], f"row rule (14 conditions ranked by the lag-0 hazard ratio among those present on all three tests, then the 5 controls) applied to the v8.1 files gives the drawn rows: {DRW['meta']['rows']}")
check(DN["offscale"] == ["Obesity hypoventilation"] and DO["offscale"] == ["Obesity hypoventilation"], f"the generator's off-scale exclusion (lag-0 hazard ratio above 2.3) removes the same condition on old and new files: {DN['offscale']}")
info(f"row set: entering {DRW['meta']['rowset_delta']['entering']}, leaving {DRW['meta']['rowset_delta']['leaving']}; conditions absent from the v7 sleep file and so never candidates: {DN['excluded_by_sleep_file']}")
check(set(DN["controls"]) == {"Alopecia", "Inguinal hernia", "Glaucoma", "Contact dermatitis", "Hemorrhoids"}, f"controls = the five v8.1 controls in the definitions order {DN['controls']}")
XLO, XHI = DRW["meta"]["xlim"]; P = DRW["meta"]["panels_used"]
allv = [(rw["condition"], p, k, rw[p][k]) for rw in DN["rows"] for p in "abc" if rw[p] is not None for k in ("before", "after", "lo", "hi")]
check(all(XLO < v[3] < XHI for v in allv), f"all {len(allv)} drawn values lie inside the x-limits {XLO:.4f} to {XHI:.4f} (V13 {DRW['meta']['xlim_v13'][0]:.4f} to {DRW['meta']['xlim_v13'][1]:.4f}, extended to the left, declared)")

# ------------------------------------------------------------------ the new sheet
doc = fitz.open(hydrated(NEW)); page = doc[0]
check(len(doc) == 1 and abs(page.rect.width - W) < 0.01 and abs(page.rect.height - H) < 0.01 and page.rotation == 0, f"page box equal within 0.01 pt: {page.rect.width:.3f} x {page.rect.height:.3f} (V13 {W:.3f} x {H:.3f}), 1 page, rotation 0")
check(len(page.get_xobjects()) == 0, f"single flat page: {len(page.get_xobjects())} XObjects (V13 0)")
fitz.TOOLS.mupdf_warnings(reset=True); page.get_pixmap(dpi=40); wmsg = fitz.TOOLS.mupdf_warnings(reset=True)
check(not wmsg, "MuPDF renders without warnings" + (f": {wmsg}" if wmsg else ""))
faces_new = sorted(set(f[3].split("+")[-1] for f in page.get_fonts(full=True))); faces_base = sorted(set(f["basefont"].split("+")[-1] for f in BV["fonts"]))
check(set(faces_new) <= {"ArialMT", "Arial-BoldMT", "Arial Bold"} and set(faces_new) == set(faces_base), f"fonts {faces_new} = the V13 set {faces_base}")
def spans_of(pg):
    """round 49: MuPDF joins two strings closer than about a third of an em into one span with a synthetic space ("1.2 1.5", the tick labels
    2.7 pt apart at 11 pt); a synthetic space has the width of the gap it bridges, a typed one the font's advance (0.278 em): split there."""
    out = []
    for b in pg.get_text("rawdict")["blocks"]:
        if b["type"] != 0: continue
        for l in b["lines"]:
            for s in l["spans"]:
                groups, cur = [], []
                for i, ch in enumerate(s["chars"]):
                    if ch["c"] == " " and i + 1 < len(s["chars"]) and abs((ch["bbox"][2] - ch["bbox"][0]) - 0.278 * s["size"]) > 0.15: groups.append(cur); cur = []; continue
                    cur.append(ch)
                groups.append(cur)
                for g in groups:
                    txt = "".join(ch["c"] for ch in g)
                    if not txt.strip(): continue
                    org = g[0]["origin"]; bb = [min(c["bbox"][0] for c in g), min(c["bbox"][1] for c in g), max(c["bbox"][2] for c in g), max(c["bbox"][3] for c in g)]
                    out.append(dict(text=txt, font=s["font"], size=round(s["size"], 3), bbox=[round(v, 3) for v in bb], origin=[round(org[0], 3), round(org[1], 3)], color="#%06x" % s["color"], asc=round(s["ascender"], 4), desc=round(s["descender"], 4)))
    return out
NS = spans_of(page); BS = [s for s in BV["spans"] if s.get("visible")]
check(collections.Counter((s["font"], s["asc"], s["desc"]) for s in NS if s["font"] in ("ArialMT", "Arial-BoldMT")) .keys() <= {("ArialMT", 0.728, -0.21), ("Arial-BoldMT", 1.01, -0.376)}, f"text-layer ascender/descender per face = the V13 faces: {dict(collections.Counter((s['font'], s['asc'], s['desc']) for s in NS))}")
check(set(s["color"] for s in NS) <= set(s["color"] for s in BS), f"text colours {sorted(set(s['color'] for s in NS))} within the V13 set")
sizes_new, sizes_base = sorted(set(s["size"] for s in NS)), sorted(set(s["size"] for s in BS))
sizes_exp = sorted({new_size(s["size"]) for s in BS if s["text"] not in NOTES} | ({10.0} if KEEP else set()))
check(len(sizes_new) == len(sizes_exp) and all(abs(a - b) < 0.002 for a, b in zip(sizes_new, sizes_exp)), f"type sizes {sizes_new} = V13 sizes {sizes_base} + 1 pt (stretched strings {STRETCHED}, kept label 10.0, note size 8.004 absent)")
LB = {s["text"]: s for s in BS if abs(s["size"] - 13.0) < 0.01}; LN = {s["text"]: s for s in NS if abs(s["size"] - 13.0 - DELTA) < 0.01}
for L in "abc":
    n, b = LN.get(L), LB.get(L)
    check(n and b and abs(n["origin"][0] - b["origin"][0]) < TOL and abs(n["origin"][1] - b["origin"][1]) < TOL and "Bold" in n["font"], f"letter {L}: {n and n['font']} {13 + DELTA:g} pt at {n and n['origin']} (V13 {b and b['origin']}, 13 pt)")
# fixed strings (everything but the data strings and the tick labels) at the V13 origins
def _sz(s, v): return abs(s["size"] - v) < 0.01 or abs(s["size"] - v - DELTA) < 0.01          # a V13 span or its +1 pt twin
def is_row_label(s): return abs(s["origin"][0] - 12.96) < 0.06 and _sz(s, 10.0) and s["text"] != "Negative controls"
def is_value(s): return _sz(s, 9.5) and (abs(s["origin"][0] - 248.69) < 0.06 or abs(s["origin"][0] - 286.13) < 0.06) and s["text"] not in ("1-year", "5-year")
def is_tick(s): return _sz(s, 10.0) and 425 < s["origin"][1] < 436
def is_ne_c(s): return s["text"] == "NE" and s["origin"][0] > 440 and s["origin"][1] < 425
def is_data(s): return is_row_label(s) or is_value(s) or is_tick(s) or is_ne_c(s)
fixed_n = [s for s in NS if not is_data(s)]; fixed_b = [s for s in BS if not is_data(s) and s["text"] not in NOTES]   # round 40: the seven note strings are absent by design
um = []; md = 0.0
for s in fixed_n:
    c = [b for b in fixed_b if b["text"] == s["text"] and abs(new_size(b["size"], b["text"]) - s["size"]) < 0.06]
    dm = min((max(abs(b["origin"][0] - s["origin"][0]), abs(b["origin"][1] - s["origin"][1])) for b in c), default=99.0)
    if dm >= TOL: um.append((s["text"], s["origin"], round(dm, 3)))
    else: md = max(md, dm)
check(not um and len(fixed_n) == len(fixed_b), f"every fixed string (letters, key letters and labels, headers, captions, 'Negative controls') sits within {TOL} pt of its V13 origin at +1 pt (max {md:.3f} pt, {len(fixed_n)} vs {len(fixed_b)} spans; the {len(NOTES)} note strings absent, round 40){': ' + str(um[:5]) if um else ''}")
check(all(abs(s["size"] - 10.0) < 0.01 for s in NS if s["text"] in KEEP) and all(abs(s["size"] - 11.0) < 0.01 for s in NS if is_row_label(s) and s["text"] not in KEEP), f"row labels 11 pt except the kept label {sorted(KEEP)} at 10 pt (at 11 pt it would run into its own row's square and interval)")
check(not [s for s in NS if s["text"] in NOTES], f"the two printed notes are not on the sheet (round 40 LNOTES removal re-applied): {NOTES}")
# tick labels: at the builder's positions (moved with the extended axis), same 15 strings
tk_n = sorted([s for s in NS if is_tick(s)], key=lambda s: s["origin"][0]); tk_b = sorted([s for s in BS if is_tick(s)], key=lambda s: s["origin"][0])
tk_rec = sorted([d for d in DRW["drawn"] if d["kind"] == "text" and d.get("what", "").startswith("tick label")], key=lambda d: d["x"])
dev_tk = max(max(abs(s["origin"][0] - d["x"]), abs(s["origin"][1] - d["baseline"])) for s, d in zip(tk_n, tk_rec)) if len(tk_n) == len(tk_rec) else 99
check(len(tk_n) == len(tk_b) == 15 and [s["text"] for s in tk_n] == [s["text"] for s in tk_b] and dev_tk < TOL, f"15 tick labels, the V13 strings, at the builder's positions on the extended axis (max deviation {dev_tk:.3f} pt; moved right by up to {max(abs(s['origin'][0] - b['origin'][0]) for s, b in zip(tk_n, tk_b)):.2f} pt, declared)")
# data strings on their slots
lab_n = sorted([s for s in NS if is_row_label(s)], key=lambda s: s["origin"][1]); lab_b = sorted([s for s in BS if is_row_label(s)], key=lambda s: s["origin"][1])
check([s["text"] for s in lab_n] == DRW["meta"]["rows"] and all(abs(n["origin"][1] - b["origin"][1]) < TOL for n, b in zip(lab_n, lab_b)), f"19 row labels in the builder's order on the V13 row slots: {[s['text'] for s in lab_n]}")
check(all(("#8a9099" == s["color"]) == (s["text"] in DN["controls"]) for s in lab_n), "control labels grey, condition labels ink")
v1_n = sorted([s for s in NS if is_value(s) and abs(s["origin"][0] - 248.69) < 0.06], key=lambda s: s["origin"][1]); v5_n = sorted([s for s in NS if is_value(s) and abs(s["origin"][0] - 286.13) < 0.06], key=lambda s: s["origin"][1])
check([s["text"] for s in v1_n] == [r["a"]["y1_text"] for r in DN["rows"]] and [s["text"] for s in v5_n] == [r["a"]["y5_text"] for r in DN["rows"]], f"printed 1-year and 5-year values = the v8.1 lag ladder half up on 2 decimals (NE where not estimable): 1-year {[s['text'] for s in v1_n]}; 5-year {[s['text'] for s in v5_n]}")
check(all((s["color"] == "#8a9099") == bool(r["a"]["y5_grey"]) for s, r in zip(v5_n, DN["rows"])), "5-year values grey exactly where the re-fit is not estimable or not informative")
ne_c = [s for s in NS if is_ne_c(s)]
check(len(ne_c) == 0, f"panel c: no NE string (v8.2: {len(ne_c)} NE strings on the sheet, expected 0)")

# ------------------------------------------------------------------ vectors: counts, styles, read-back
DW = page.get_drawings()
def near(a, b, tol=0.002): return a is not None and b is not None and all(abs(x - y) < tol for x, y in zip(a, b))
INKc = [0.1021, 0.114, 0.1289]; BLUEc = [0.008, 0.5332, 0.8203]; GREYc = [0.541, 0.5645, 0.5996]; PALEc = [0.8008, 0.8203, 0.8398]
top, bottom = G["plot_top"] - 2, G["plot_bottom"] + 1
sq = [d for d in DW if "".join(it[0] for it in d["items"]) == "re" and near(d.get("fill"), PALEc) and d["rect"].y0 > top and d["rect"].y1 < bottom and d["rect"].x0 < 575]
cc = [d for d in DW if "".join(it[0] for it in d["items"]) == "cccccccc" and near(d.get("fill"), BLUEc) and d["rect"].x0 < 575 and d["rect"].y0 > top]
ci = [d for d in DW if "".join(it[0] for it in d["items"]) == "l" and near(d.get("color"), BLUEc) and d["rect"].x0 < 575 and abs(d["width"] - 1.7) < 0.01]
n_expected = sum(1 for rw in DN["rows"] for p in "abc" if rw[p] is not None)
check(len(sq) == len(cc) == len(ci) == n_expected == 57, f"markers: {len(sq)} squares, {len(cc)} circles, {len(ci)} intervals (19 rows x 3 panels, every row with a value; V13 57)")
check(all(abs(d["rect"].width - G["square_side"]) < 0.05 for d in sq) and all(abs(d["rect"].width - G["circle_diam"]) < 0.05 for d in cc), f"marker sizes: squares {G['square_side']:.3f} pt, circles {G['circle_diam']:.3f} pt (V13)")
refs = [d for d in DW if d.get("dashes") and d["dashes"] not in ("[] 0",) and d["rect"].height > 100]
check(len(refs) == 3 and all(abs(d["rect"].x0 - P[n]["ref_x"]) < 0.05 for d, n in zip(sorted(refs, key=lambda d: d["rect"].x0), "abc")), f"three dashed reference rules at x(1) of the extended mapping {[round(P[n]['ref_x'], 2) for n in 'abc']} (V13 {[round(G['panels'][n]['ref_x'], 2) for n in 'abc']}, declared)")
spines = sorted([d for d in DW if near(d.get("color"), INKc) and abs(d["width"] - 0.8) < 0.01 and d["rect"].width > 50], key=lambda d: d["rect"].x0)
check(len(spines) == 3 and all(abs(d["rect"].x0 - G["panels"][n]["spine"][0]) < TOL and abs(d["rect"].x1 - G["panels"][n]["spine"][2]) < TOL and abs(d["rect"].y0 - G["panels"][n]["spine"][1]) < TOL for d, n in zip(spines, "abc")), "three bottom spines at the V13 boxes")
agree = ["ED_Fig02 V14: values read back from the vectors against the numbers files (percent of the axis span)"]; worst = 0.0; nread = 0
span_ln = math.log(XHI) - math.log(XLO)
for i, rw in enumerate(DN["rows"]):
    yc = G["row_centres"][i]; ysq, ycir = yc - G["lane_half"], yc + G["lane_half"]
    for n in "abc":
        v = rw[n]
        if v is None: continue
        p = P[n]; inv = lambda x: math.exp((x - p["a"]) / p["perln"])
        x0, x1 = p["spine"][0] - 1, p["spine"][2] + (73 if n == "a" else 1)
        s_ = [d for d in sq if abs((d["rect"].y0 + d["rect"].y1) / 2 - ysq) < 0.06 and x0 <= (d["rect"].x0 + d["rect"].x1) / 2 <= x1]
        c_ = [d for d in cc if abs((d["rect"].y0 + d["rect"].y1) / 2 - ycir) < 0.06 and x0 <= (d["rect"].x0 + d["rect"].x1) / 2 <= x1]
        i_ = [d for d in ci if abs(d["rect"].y0 - ycir) < 0.06 and x0 <= (d["rect"].x0 + d["rect"].x1) / 2 <= x1]
        assert len(s_) == 1 and len(c_) == 1 and len(i_) == 1, (rw["condition"], n, len(s_), len(c_), len(i_))
        got = (inv((s_[0]["rect"].x0 + s_[0]["rect"].x1) / 2), inv((c_[0]["rect"].x0 + c_[0]["rect"].x1) / 2), inv(i_[0]["rect"].x0), inv(i_[0]["rect"].x1))
        exp = (v["before"], v["after"], v["lo"], v["hi"])
        dev = max(abs(math.log(g) - math.log(e)) for g, e in zip(got, exp)) / span_ln * 100; worst = max(worst, dev); nread += 4
        agree.append(f"{n} {rw['condition']:28s} file {exp[0]:.3f} {exp[1]:.3f} ({exp[2]:.3f}-{exp[3]:.3f})  read {got[0]:.3f} {got[1]:.3f} ({got[2]:.3f}-{got[3]:.3f})  dev {dev:.3f}%")
open(os.path.join(VER, "values_agreement.txt"), "w").write("\n".join(agree) + "\n")
check(nread == 4 * n_expected == 228 and worst <= 0.5, f"NEW sheet: {nread} values (square = before, circle = after, interval ends) read back through the mapping agree with lag_ladder.csv (a), causal_tests_all.csv (b) and the v8.2 sleep_unadjusted_v1.json (c): max deviation {worst:.3f}% of the axis span")

# ------------------------------------------------------------------ positive control: the OLD replot against V13
docO = fitz.open(hydrated(OLDP)); OSP = spans_of(docO[0]); ODW = docO[0].get_drawings()
um = []; md = 0.0
for s in BS:
    c = [o for o in OSP if o["text"] == s["text"] and abs(o["size"] - s["size"]) < 0.06]
    dm = min((max(abs(o["origin"][0] - s["origin"][0]), abs(o["origin"][1] - s["origin"][1])) for o in c), default=99.0)
    if dm >= TOL: um.append((s["text"], s["origin"], round(dm, 3)))
    else: md = max(md, dm)
check(not um and len(OSP) == len(BS), f"positive control: the same builder on the v7 snapshot (work/ED_Fig02_OLDREPLOT.pdf) reproduces every visible V13 string at its origin within {TOL} pt (max {md:.3f} pt, {len(OSP)} vs {len(BS)} strings){': ' + str(um[:5]) if um else ''}")
bsq = [d for d in G["key_drawings"]]  # unused placeholder
base_dw = json.load(open(hydrated(os.path.join(WORK, "base_drawings_seq.json"))))
bcc = [d for d in base_dw if d["kinds"] == "cccccccc" and near(d["fill"], BLUEc) and d["rect"][0] < 575 and d["rect"][1] > top]
occ = [d for d in ODW if "".join(it[0] for it in d["items"]) == "cccccccc" and near(d.get("fill"), BLUEc) and d["rect"].x0 < 575 and d["rect"].y0 > top]
dm_pc = max(min(math.hypot((o["rect"].x0 + o["rect"].x1) / 2 - (b["rect"][0] + b["rect"][2]) / 2, (o["rect"].y0 + o["rect"].y1) / 2 - (b["rect"][1] + b["rect"][3]) / 2) for b in bcc) for o in occ) if occ and len(occ) == len(bcc) else 99
check(dm_pc < TOL, f"positive control: the OLD replot's {len(occ)} circles sit on V13's {len(bcc)} visible circles within {TOL} pt (max centre deviation {dm_pc:.3f} pt)")
docO.close()
imO = np.asarray(Image.open(os.path.join(WORK, "ED_Fig02_OLDREPLOT_150dpi.png")).convert("RGB")).astype(int); im13 = np.asarray(Image.open(hydrated(os.path.join(WORK, "base_150.png"))).convert("RGB")).astype(int) if os.path.exists(os.path.join(WORK, "base_150.png")) else None
if im13 is None:
    bdoc = fitz.open(hydrated(BASE)); pm = bdoc[0].get_pixmap(dpi=150, colorspace=fitz.csRGB, alpha=False); pm.save(os.path.join(WORK, "base_150.png")); pm = None; bdoc.close()
    im13 = np.asarray(Image.open(os.path.join(WORK, "base_150.png")).convert("RGB")).astype(int)
if imO.shape == im13.shape:
    frac = ((np.abs(imO - im13) > 96).any(axis=2)).mean() * 100
    check(frac < 1.0, f"positive control raster: OLD replot vs V13 at 150 dpi differ on {frac:.3f}% of pixels (font subsetting and anti-aliasing only)")
else: check(False, f"positive control raster shapes differ {imO.shape} vs {im13.shape}")
# NEW vs V13 raster: the key column (x > 578, y < 280) must be pixel-identical beyond the noise floor
imN = np.asarray(Image.open(os.path.join(WORK, "ED_Fig02_composed_150dpi.png")).convert("RGB")).astype(int)   # round 49: the composed sheet (the V13 page box; the delivered sheet carries the 16 pt strip)
s150 = 150 / 72.0; kx0, kx1, ky1 = int(598 * s150), int(614 * s150), int(178 * s150)      # the key handles, squares and circles (x 600.4 to 611.6, y 44 to 167); letters at 587, text from 616, the V13 notes from y 182 (absent here)
dkm = (np.abs(imN[:ky1, kx0:kx1] - im13[:ky1, kx0:kx1]) > 96).any(axis=2); dk = int(dkm.sum())
Image.fromarray(np.uint8(np.where(dkm[..., None], [255, 0, 0], imN[:ky1, kx0:kx1]))).resize((int(3 * (kx1 - kx0)), int(3 * ky1)), Image.NEAREST).save(os.path.join(WORK, "key_window_diff_150dpi_x3.png"))
info(f"150 dpi raster, key markers window (x 598 to 614, y < 178) against V13: {dk} of {dkm.size} pixels differ ({100 * dk / dkm.size:.3f}%): the ends of the a and b key handles, which V13 clips to 7.28 pt and the builder has drawn full width with projecting caps since round 30 (every delivered sheet since V14 carries the builder's handles; run_sheet.py checks the window against V26)")
dkm2 = (np.abs(imN[:ky1, kx0:kx1] - imO[:ky1, kx0:kx1]) > 96).any(axis=2); dk2 = int(dkm2.sum())
check(dk2 < 0.002 * dkm2.size, f"150 dpi raster: the key markers (x 598 to 614, y < 178) equal the OLD replot's (the same construction at the V26 sizes) beyond the anti-aliasing floor ({dk2} of {dkm2.size} pixels differ, {100 * dk2 / dkm2.size:.3f}%)")
# word / numeric deltas: observed = the data strings' delta
def wm(ss): return collections.Counter(w for t in ss for w in t.split())
def nm(ss): return collections.Counter(x for t in ss for x in re.findall(r"-?\d[\d,]*\.?\d*", t))
data_new = [s["text"] for s in NS if is_data(s) and not is_tick(s)]; data_base = [s["text"] for s in BS if is_data(s) and not is_tick(s)]
obs_wa, obs_wr = wm([s["text"] for s in NS]) - wm([s["text"] for s in BS]), wm([s["text"] for s in BS]) - wm([s["text"] for s in NS])
exp_wa, exp_wr = wm(data_new) - wm(data_base), (wm(data_base) - wm(data_new)) + wm(NOTES)      # round 40: the notes' words leave
obs_na, obs_nr = nm([s["text"] for s in NS]) - nm([s["text"] for s in BS]), nm([s["text"] for s in BS]) - nm([s["text"] for s in NS])
exp_na, exp_nr = nm(data_new) - nm(data_base), (nm(data_base) - nm(data_new)) + nm(NOTES)
# round 49: full multiset identities (NS = fixed + data_new; BS = fixed + data_base + notes), so a word that leaves with the notes and returns as a data string nets correctly
check(wm([s["text"] for s in NS]) + wm(NOTES) + wm(data_base) == wm([s["text"] for s in BS]) + wm(data_new), f"word multiset: new + the removed notes + the V13 data strings = V13 + the new data strings (fixed strings identical; data words added {sum(exp_wa.values())}, removed {sum((wm(data_base) - wm(data_new)).values())}, note words removed {sum(wm(NOTES).values())})")
check(nm([s["text"] for s in NS]) + nm(NOTES) + nm(data_base) == nm([s["text"] for s in BS]) + nm(data_new), f"numeric-token multiset: the same identity (data tokens added {sum(exp_na.values())}, removed {sum((nm(data_base) - nm(data_new)).values())}, note tokens removed {sum(nm(NOTES).values())})")
json.dump(dict(words_added=dict(exp_wa), words_removed=dict(exp_wr), numbers_added=dict(exp_na), numbers_removed=dict(exp_nr)), open(os.path.join(VER, "expected_delta.json"), "w"), indent=1)
# overlaps and ink inside the page
def overlaps(sp):
    out = []
    for i, a in enumerate(sp):
        for b in sp[i + 1:]:
            ax0, ay0, ax1, ay1 = a["bbox"]; bx0, by0, bx1, by1 = b["bbox"]
            if ax0 < bx1 - 0.5 and bx0 < ax1 - 0.5 and ay0 < by1 - 1 and by0 < ay1 - 1: out.append((a["text"], b["text"]))
    return out
ov_n, ov_b = overlaps(NS), overlaps(BS)
check(len(ov_n) <= len(ov_b), f"no text span overlaps another beyond V13's own ({len(ov_n)} overlapping pairs on the new sheet, {len(ov_b)} on V13){': ' + str(ov_n[:4]) if len(ov_n) > len(ov_b) else ''}")
ink = fitz.Rect()
for d in DW:
    if d.get("fill") is not None and near(d["fill"], [1.0, 1.0, 1.0]) and d.get("color") is None: continue   # the page background
    ink |= d["rect"]
for s in NS: ink |= fitz.Rect(s["bbox"])
check(ink.x0 > 2 and ink.y0 > 2 and ink.x1 < W - 2 and ink.y1 < H - 2, f"all ink inside the page with a 2 pt margin: {[round(v, 2) for v in ink]}")
check(all(os.path.exists(os.path.join(CROPS, f"{n}_OLD_vs_NEW_200dpi.png")) for n in ("panel_a", "panel_b", "panel_c", "key")), "200 dpi OLD (V13) vs NEW crops written for the three panels and the key")

# round 49: the v7-vs-v8 CHANGES block is not part of this round (no number changes); CHANGES_ED_Fig02.csv lists the typographic deltas (run_sheet.py)

# ------------------------------------------------------------------ printed values CSV
recs = []
for d in DRW["drawn"]:
    if d["kind"] != "text": continue
    what = d.get("what", ""); c = d.get("condition")
    if what == "row label": recs.append(dict(text=d["text"], x=d["x"], baseline=d["baseline"], ha="left", size=d["size"], source_file=f"{NUM}/lag_ladder.csv", source_key=f"condition={c}&lag_years=0.0:condition", source_value=c, rule="text", panel="a"))
    elif what == "1-year value": recs.append(dict(text=d["text"], x=d["x"], baseline=d["baseline"], ha="left", size=d["size"], source_file=f"{NUM}/lag_ladder.csv", source_key=f"condition={c}&lag_years=1.0:hr", source_value=d.get("value"), rule="2dp_halfup", panel="a"))
    elif what == "5-year value":
        if d["text"] == "NE": recs.append(dict(text="NE", x=d["x"], baseline=d["baseline"], ha="left", size=d["size"], source_file=f"{NUM}/lag_ladder.csv", source_key=f"NE: condition={c} lag_years=5.0 estimable False", source_value="NE", rule="text", panel="a"))
        else: recs.append(dict(text=d["text"], x=d["x"], baseline=d["baseline"], ha="left", size=d["size"], source_file=f"{NUM}/lag_ladder.csv", source_key=f"condition={c}&lag_years=5.0:hr", source_value=d.get("value"), rule="2dp_halfup", panel="a"))
    elif what.startswith("panel c NE"): recs.append(dict(text="NE", x=d["x"], baseline=d["baseline"], ha="left", size=d["size"], source_file=f"{NUM}/sleep_unadjusted_v1.json", source_key=f"NE: outcomes/{c} absent (group P v7 file, control added 14 September)", source_value="NE", rule="text", panel="c"))
    else: recs.append(dict(text=d["text"], x=d["x"], baseline=d["baseline"], ha="left", size=d["size"], source_file="static:V13", source_key=what, source_value=d["text"], rule="text", panel=d.get("panel", "")))
for L, s in LB.items(): recs.append(dict(text=L, x=s["origin"][0], baseline=s["origin"][1], ha="left", size=13.0, source_file="static:V13", source_key="panel letter", source_value=L, rule="text", panel=L))
# round 49: MuPDF's dict extraction joins the tick labels "1.2" and "1.5" (2.7 pt apart at 11 pt) into one span; the CSV matcher reads spans
# through v14lib.spans, so the sheet's spans are read from the rawdict and split at a typed-looking space that hides a real gap over 1 pt
import v14lib as _vl
def _spans_split(path, page_no=0):
    d = fitz.open(hydrated(path)); out = []
    for b in d[page_no].get_text("rawdict")["blocks"]:
        if b["type"] != 0: continue
        for l in b["lines"]:
            for s in l["spans"]:
                chars = s["chars"]; groups = []; cur = []
                for i, ch in enumerate(chars):
                    # a synthetic space: MuPDF gives it the width of the gap it bridges, a typed space the font's advance (0.278 em)
                    if ch["c"] == " " and i + 1 < len(chars) and abs((ch["bbox"][2] - ch["bbox"][0]) - 0.278 * s["size"]) > 0.15: groups.append(cur); cur = []; continue
                    cur.append(ch)
                groups.append(cur)
                for g in groups:
                    t = "".join(c["c"] for c in g).replace("\xa0", " ").strip()
                    if not t: continue
                    x0 = min(c["bbox"][0] for c in g); x1 = max(c["bbox"][2] for c in g); y0 = min(c["bbox"][1] for c in g); y1 = max(c["bbox"][3] for c in g)
                    out.append((t, g[0]["origin"][0], g[0]["origin"][1], s["font"], s["size"], s["color"], (x0, y0, x1, y1)))
    d.close(); gc.collect(); return out
_vl.spans = _spans_split
n_rows, n_no = printed_values_csv(NEW, recs, os.path.join(VER, "ED_Fig02_printed_values.csv"), static_from=BASE)
check(n_no == 0, f"printed values: {n_rows} strings on the sheet accounted for (drawn record with source, or identical to V13), {n_no} unaccounted (NO)")
doc.close(); ok = ck.write(os.path.join(VER, "checks_03_verify.txt")); print("RESULT ALL PASS" if ok else "RESULT FAIL")   # round 49: run_sheet.py merges this into checks.txt
