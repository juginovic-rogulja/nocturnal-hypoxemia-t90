#!/usr/bin/env python3
"""Supp_Fig04, V14 lane L2 (round 37), step 3: Proof W of Supp_Fig04.pdf against the V13 sheet and the v8.1 lag ladder (rewritten from the
round-30 verifier, byte copy beside as 03_verify_PRE_V8_1.py). Declared deltas: the row set (35 conditions, 18 in a and 17 + 5 controls in b,
the labels that enter and leave), the page height (+ the growth), every string under panel a moved by the panel-a growth, under panel b by
the panel-b growth, the c/d band and letters c, d by the larger, the panel c annotations and panel d counts. Writes verify/checks.txt,
values_agreement.txt, expected_delta.json, CHANGES_Supp_Fig04.csv, Supp_Fig04_printed_values.csv."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import csv, gc, json, math, os, re, sys, time, collections
import numpy as np, pandas as pd, fitz
from PIL import Image
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND37_2026-09-14/figures/_common")
from v14lib import V13, NUM, OLD_V7, hydrated, sha256, sidecar_v8, printed_values_csv, Checks, write_changes
LANE = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); WORK = os.path.join(LANE, "work"); VER = os.path.join(LANE, "verify"); CROPS = os.path.join(VER, "crops")
BASE = f"{V13}/Supp_Fig04.pdf"; NEW = os.path.join(LANE, "Supp_Fig04.pdf"); OLDP = os.path.join(WORK, "Supp_Fig04_OLD.pdf")
DATA = f"{NUM}/lag_ladder.json"; DATA_CSV = f"{NUM}/lag_ladder.csv"; OLD_DATA = f"{OLD_V7}/numbers/lag_ladder.json"
DR = json.load(open(hydrated(os.path.join(VER, "Supp_Fig04_drawn.json")))); DRO = json.load(open(hydrated(os.path.join(VER, "positive_control_drawn.json"))))
G = json.load(open(hydrated(os.path.join(WORK, "base_geometry.json")))); BT = json.load(open(hydrated(os.path.join(WORK, "base_text.json"))))
GR = DR["growth"]; TOL = 0.3; XLO, XHI = DR["xlim"]; XLO13, XHI13 = DR["xlim_v13"]
PT_PLUS = 1.0; REMOVED_R40 = list(DR.get("removed_r40_notes", []))   # round 49: every text one point larger; the two panel c key lines removed in round 40 (LNOTES) stay off
def sz(s, v): return abs(s["size"] - v) < 0.01 or abs(s["size"] - v - PT_PLUS) < 0.01   # a V13 size on the base, that size plus one on the new sheet
ck = Checks(f"Supp_Fig04 V14 lane L2 (round 37), Proof W: {NEW} (sha256 {sha256(NEW)[:16]}) against the V13 sheet {BASE} (sha256 {sha256(BASE)[:16]}) and numbers/lag_ladder.json (step 130). "
            f"Declared: rows a {GR['rows_a']} (V13 {GR['rows_a_v13']}), b {GR['rows_b']} (V13 {GR['rows_b_v13']}), growth a {GR['a']} b {GR['b']} band {GR['cd']} pt, page {DR['page_v13'][1]} -> {DR['page'][1]}, labels entering {DR['label_delta']['entering']}, leaving {DR['label_delta']['leaving']}, x limits {XLO} to {XHI} (V13 {XLO13} to {XHI13}, the 0.96 x minimum rule). {time.strftime('%Y-%m-%d %H:%M')}")
def check(ok, msg): return ck.log(msg, bool(ok))
def info(msg): ck.lines.append("INFO: " + msg); print("INFO: " + msg[:300], flush=True)
def spans_of(pdf):
    d = fitz.open(hydrated(pdf)); pg = d[0]; out = []
    for b in pg.get_text("rawdict")["blocks"]:
        if b["type"] != 0: continue
        for l in b["lines"]:
            for s in l["spans"]:
                txt = "".join(c["c"] for c in s["chars"]); org = s["chars"][0]["origin"]
                if txt.strip(): out.append(dict(text=txt, font=s["font"], size=round(s["size"], 3), bbox=[round(v, 3) for v in s["bbox"]], origin=[round(org[0], 3), round(org[1], 3)], color="#%06x" % s["color"], asc=round(s["ascender"], 4), desc=round(s["descender"], 4), dir=[round(v, 3) for v in l["dir"]]))
    words = [w[4] for w in pg.get_text("words")]; fonts = sorted({f[3].split("+")[-1] for f in pg.get_fonts(full=True)}); rect = [round(v, 3) for v in pg.rect]; nx = len(pg.get_xobjects()); n = len(d)
    d.close(); gc.collect(); return out, words, fonts, rect, nx, n
def collapse(sp, tol=0.25):
    keep = []
    for s in sp:
        if not any(k["text"] == s["text"] and abs(k["origin"][0] - s["origin"][0]) <= tol and abs(k["origin"][1] - s["origin"][1]) <= tol for k in keep): keep.append(s)
    return keep
def hx(c): return None if c is None else "#%02x%02x%02x" % tuple(int(round(v * 255)) for v in c)
def drawings_of(pdf):
    d = fitz.open(hydrated(pdf)); out = []
    for g in d[0].get_drawings():
        r = g["rect"]; out.append(dict(type=g["type"], rect=[r.x0, r.y0, r.x1, r.y1], fill=hx(g["fill"]), stroke=hx(g["color"]), width=g.get("width"), dashes=g.get("dashes"), kinds="".join(sorted({it[0] for it in g["items"]})), n_items=len(g["items"])))
    d.close(); gc.collect(); return out

# 1 provenance and the data re-derived
J = json.load(open(hydrated(DATA))); PC = J["per_condition"]; JO = json.load(open(hydrated(OLD_DATA))); PCO = JO["per_condition"]
sc = sidecar_v8(DATA); scc = sidecar_v8(DATA_CSV)
check(sc["sha256"] == DR["source"]["sha256"] and scc["sha256"] == DR["source"]["csv_sha256"] and sc["v8_inputs"], f"lag_ladder.json and .csv read now = the builder's record = the v8.1 sidecars (step {sc['step']} {sc['name']}, out {sc['output_mtime']})")
check(DR["panel_n"] == J["panel_n"] == 35 and DR["half"] == math.ceil(J["panel_n"] / 2) == 18, f"panel_n {J['panel_n']} (V13 {JO['panel_n']}), ceil(n/2) = {DR['half']} rows in a")
d0 = pd.read_csv(hydrated(DATA_CSV)); w0 = d0[d0.lag_years == 0]
panel = list(w0[(w0.lo > 1) & ~w0.negative_control & ~w0.circular].sort_values("hr", ascending=False).key); ctrl = list(w0[w0.negative_control].sort_values("hr", ascending=False).key)
rows_a = [PC[k]["condition"] for k in panel[:18]]; rows_b = [PC[k]["condition"] for k in panel[18:]] + ["Negative controls"] + [PC[k]["condition"] for k in ctrl]
check(XLO <= XLO13 and XHI >= XHI13 and DRO["xlim"] == [XLO13, XHI13], f"x limits: new {XLO} to {XHI}, V13 {XLO13} to {XHI13} (extended to the left only where a v8.1 interval runs below 0.62; the OLD replot keeps the V13 limits)")
check(DR["panels"]["a"]["rows"] == rows_a and DR["panels"]["b"]["rows"] == rows_b, f"row sets re-derived from the csv by the generator's rule equal the drawn rows: a {rows_a[:3]}... b ends {rows_b[-5:]}")
check(set(PC[k]["condition"] for k in ctrl) == {"Alopecia", "Inguinal hernia", "Glaucoma", "Contact dermatitis", "Hemorrhoids"}, f"controls are the five v8.1 controls in the lag-0 order {[PC[k]['condition'] for k in ctrl]}")

# 2 page, fonts, sizes, letters
ns, nw, nf, nrect, nnx, nn = spans_of(NEW); bs, bw, bf, brect, bnx, bn = spans_of(BASE)
check(nn == 1 and nnx == 0 and abs(nrect[2] - brect[2]) < 0.01 and abs(nrect[3] - (brect[3] + GR["cd"])) < 0.05, f"one flat page, width {nrect[2]} = V13, height {nrect[3]} = V13 {brect[3]} + {GR['cd']} (declared growth)")
fitz.TOOLS.mupdf_warnings(reset=True); dd_ = fitz.open(NEW); dd_[0].get_pixmap(dpi=40); w = fitz.TOOLS.mupdf_warnings(reset=True); dd_.close()
check(not w, "MuPDF renders without warnings" + (f": {w}" if w else ""))
check(set(nf) == set(bf) == {"ArialMT", "Arial-BoldMT", "Arial-ItalicMT"}, f"fonts {nf} = V13 {bf}")
check(collections.Counter((s["font"], s["asc"], s["desc"]) for s in ns).keys() <= collections.Counter((s["font"], s["asc"], s["desc"]) for s in bs).keys(), "text-layer ascender/descender per face within the V13 faces")
check(sorted({s["size"] for s in ns}) == sorted({round(s["size"] + PT_PLUS, 2) for s in bs}), f"type sizes {sorted({s['size'] for s in ns})} = V13 {sorted({s['size'] for s in bs})} plus one point (round 49)")
nsc, bsc = collapse(ns), collapse(bs)
LB = {s["text"]: s for s in bsc if abs(s["size"] - 13.0) < 0.01 and s["text"] in "abcd"}; LN = {s["text"]: s for s in nsc if abs(s["size"] - 13.0 - PT_PLUS) < 0.01 and s["text"] in "abcd"}
for L in "abcd":
    dy = GR["cd"] if L in "cd" else 0.0; n_, b_ = LN.get(L), LB.get(L)
    check(n_ and b_ and abs(n_["origin"][0] - b_["origin"][0]) < TOL and abs(n_["origin"][1] - (b_["origin"][1] + dy)) < TOL and "Bold" in n_["font"], f"letter {L}: 14 pt bold (V13 13 plus one) at {n_ and n_['origin']} (V13 {b_ and b_['origin']} + {dy})")

# 3 fixed strings at the V13 origin plus the region shift; data strings by the rule
def shift_y(x, y):
    if y >= 700.0: return GR["cd"]
    if x < 500.0 and y >= G["axes_box"]["a"][3] - 1.0: return GR["a"]
    if x >= 500.0 and y >= G["axes_box"]["b"][3] - 1.0: return GR["b"]
    return 0.0
def is_row_label(s): return sz(s, 10.0) and s["font"] == "ArialMT" and s["text"].endswith(")") and "(" in s["text"] and not s["text"].startswith("Hazard")
def is_c_ann(s): return sz(s, 9.5) and "/" in s["text"] and s["text"].replace("/", "").isdigit()
def is_d_count(s): return sz(s, 9.5) and s["text"].isdigit() and s["origin"][0] > 600 and s["origin"][1] > 700
def is_noest(s): return s["text"] == "No estimate"
def is_ctrl_head(s): return s["text"] == "Negative controls"
_YB = [DR["panels"][p]["axes_box"][3] + 13.5 for p in "ab"] + [G["axes_box"][p][3] + 13.5 for p in "ab"]
def is_xtick_ab(s): return sz(s, 10.0) and s["text"] in ("0.7", "1", "1.5", "2", "3", "4") and s["origin"][1] < 700 and any(abs(s["origin"][1] - y) < 3 for y in _YB)
def is_xtick_d(s): return sz(s, 10.0) and s["text"] in ("0", "5", "10", "15", "20", "25", "30") and s["origin"][1] > 700 and s["origin"][0] > 590
def is_data(s): return is_row_label(s) or is_c_ann(s) or is_d_count(s) or is_noest(s) or is_ctrl_head(s) or is_xtick_ab(s) or is_xtick_d(s)
fixed_n = [s for s in nsc if not is_data(s)]; fixed_b = [s for s in bsc if not is_data(s) and s["text"] not in REMOVED_R40]; um = []; md = 0.0
# round 49: a one-point larger string keeps its ANCHOR, not its box: the left origin (left-aligned strings), the centre (the axis titles, which the
# builder re-centres, and the matplotlib-centred tick labels) or the right edge (right-aligned tick labels); a rotated string keeps its x and its
# centre along y. The baseline of a top-aligned matplotlib tick label moves down by the ascent growth of one point (0.72 pt measured on the
# panel c x ticks), so the tolerance is 0.6 pt along the text and 1.0 pt across it.
TOL49 = 0.6; TOL49_Y = 1.0
def anchors(s):
    b = s["bbox"]; rot = abs(s["dir"][0]) < 0.5
    if rot: return [(s["origin"][0], (b[1] + b[3]) / 2)], True
    return [(s["origin"][0], s["origin"][1]), ((b[0] + b[2]) / 2, s["origin"][1]), (b[2], s["origin"][1])], False
for s in fixed_n:
    c = [b for b in fixed_b if b["text"] == s["text"] and abs(b["size"] + PT_PLUS - s["size"]) < 0.06 and b["font"] == s["font"]]
    best = 99.0
    for b in c:
        an, rot = anchors(s); ab, _ = anchors(b); dy = shift_y(*b["origin"])
        for (x1, y1), (x0, y0) in zip(an, ab):
            best = min(best, max(abs(x1 - x0), abs(y1 - (y0 + dy)) * TOL49 / TOL49_Y))
    if best >= TOL49: um.append((s["text"], s["origin"], round(best, 3)))
    else: md = max(md, best)
check(not um and len(fixed_n) == len(fixed_b), f"every fixed string (axis titles, keys, tick labels, panel d row labels, letters) keeps an anchor (left origin, centre or right edge) within {TOL49} pt along the text and {TOL49_Y} pt across it of its V13 anchor plus the declared region shift, one point larger (max {md:.3f} pt, {len(fixed_n)} vs {len(fixed_b)} spans, the two round-40 key lines declared off){': ' + str(um[:5]) if um else ''}")
check(all(not any(s["text"] == t for s in nsc) for t in REMOVED_R40) and len(REMOVED_R40) == 2, f"the two panel c key lines removed in round 40 are absent: {REMOVED_R40}")
xt_n = [s for s in nsc if is_xtick_ab(s)]; xt_bad = []
for s in xt_n:
    p = "a" if s["origin"][0] < 500 else "b"; box = DR["panels"][p]["axes_box"]; t = float(s["text"]); xc = box[0] + (math.log(t) - math.log(XLO)) / (math.log(XHI) - math.log(XLO)) * (box[2] - box[0])
    if abs((s["bbox"][0] + s["bbox"][2]) / 2 - xc) > TOL + 0.2 or abs(s["origin"][1] - (box[3] + 13.5)) > 3: xt_bad.append((p, s["text"], round((s["bbox"][0] + s["bbox"][2]) / 2, 2), round(xc, 2)))
check(len(xt_n) == 12 and not xt_bad, f"the 12 x tick labels of a and b sit centred on the ticks of the declared x mapping ({XLO} to {XHI}) under the grown axes{': ' + str(xt_bad[:4]) if xt_bad else ''}")
xd_n = sorted([s for s in nsc if is_xtick_d(s)], key=lambda s: s["origin"][0]); dbox = DR["panels"]["d"]["axes_box"]; dxl = DR["panels"]["d"]["xlim"][1]; xd_bad = []
for s in xd_n:
    xc = dbox[0] + float(s["text"]) / dxl * (dbox[2] - dbox[0])
    if abs((s["bbox"][0] + s["bbox"][2]) / 2 - xc) > TOL + 0.2 or abs(s["origin"][1] - (dbox[3] + 13.5)) > 3: xd_bad.append((s["text"], round((s["bbox"][0] + s["bbox"][2]) / 2, 2), round(xc, 2)))
check([s["text"] for s in xd_n] == [str(t) for t in [0, 5, 10, 15, 20, 25, 30] if t <= dxl] and not xd_bad, f"panel d x tick labels {[s['text'] for s in xd_n]} on the ticks of the axis to {dxl:.2f} (= 0.68 x panel_n, the generator's rule; V13 21.76){': ' + str(xd_bad) if xd_bad else ''}")
lab_n = sorted([s for s in nsc if is_row_label(s)], key=lambda s: (s["origin"][0] > 500, s["origin"][1])); exp_labels = [l["text"] for l in DR["panels"]["a"]["labels"]] + [l["text"] for l in DR["panels"]["b"]["labels"]]
check([s["text"] for s in lab_n] == exp_labels, f"{len(lab_n)} row labels in the rule's order: a {exp_labels[:2]} ... b {exp_labels[-2:]}")
pos_bad = []
for s, l in zip(lab_n, DR["panels"]["a"]["labels"] + DR["panels"]["b"]["labels"]):
    if abs(s["bbox"][2] - l["x_right"]) > TOL + 0.5 or abs((s["bbox"][1] + s["bbox"][3]) / 2 - l["y_center"]) > 3.0: pos_bad.append((s["text"], s["bbox"], l["x_right"], l["y_center"]))
check(not pos_bad, f"row labels right-aligned 4 pt left of the axes and centred on their rows (pitch {GR['pitch']} kept){': ' + str(pos_bad[:3]) if pos_bad else ''}")
hd = [s for s in nsc if is_ctrl_head(s)]; ib = DR["panels"]["b"]["row_keys"].index(None); yc = DR["panels"]["b"]["axes_box"][1] + (ib + 0.5) * GR["pitch"]
check(len(hd) == 1 and abs(hd[0]["origin"][1] - yc) < GR["pitch"], f"'Negative controls' header once, on its blank row ({hd[0]['origin'] if hd else None}, row centre {yc:.2f})")
noest_n = [s for s in nsc if is_noest(s)]; exp_noest = sum(1 for p in "ab" for r in DR["panels"][p]["rungs"] if not r["estimable"])
check(len(noest_n) == exp_noest == sum(1 for k in panel + ctrl for L in "0125" if PC[k]["hr_by_lag"][L] is None), f"{len(noest_n)} 'No estimate' strings = the rungs the file leaves unestimable")
c_ann = sorted([s for s in nsc if is_c_ann(s)], key=lambda s: s["origin"][0]); exp_ann = [a["text"] for a in DR["panels"]["c"]["annotations"]]
check([s["text"] for s in c_ann] == exp_ann == [f"{J['panel_n']}/{J['panel_n']}"] + [f"{J['by_lag'][L]['survive']}/{J['panel_n']}" for L in "125"], f"panel c annotations {[s['text'] for s in c_ann]} = survive/panel_n from the file")
d_cnt = sorted([s for s in nsc if is_d_count(s)], key=lambda s: s["origin"][1]); exp_cnt = [t["text"] for t in DR["panels"]["d"]["count_texts"]]
check([s["text"] for s in d_cnt] == exp_cnt, f"panel d counts {[s['text'] for s in d_cnt]} = the trend classes of the file (sum {sum(int(t) for t in exp_cnt)} = panel_n)")

# 4 vectors: every rung's marker and interval where the record says, values through the log mapping
DW = drawings_of(NEW)
def near_pt(a, b, tol=0.35): return abs(a[0] - b[0]) <= tol and abs(a[1] - b[1]) <= tol
fills = [g for g in DW if g["fill"] in ("#0288d1", "#ffffff") and g["rect"][2] - g["rect"][0] < 12 and g["rect"][3] - g["rect"][1] < 12 and g["type"] in ("f", "fs")]
cis = [g for g in DW if g["stroke"] == "#0288d1" and g["width"] is not None and abs(g["width"] - 1.7) < 0.01 and g["rect"][3] - g["rect"][1] < 0.05]
bad = []; nread = 0; worst = 0.0
for p in "ab":
    box = DR["panels"][p]["axes_box"]
    for r in DR["panels"][p]["rungs"]:
        if not r["estimable"]: continue
        c = [g for g in fills if near_pt(((g["rect"][0] + g["rect"][2]) / 2, (g["rect"][1] + g["rect"][3]) / 2), (r["x_page"], r["y_page"]))]
        l = [g for g in cis if abs(g["rect"][1] - r["y_page"]) < 0.35 and abs(g["rect"][0] - r["x_lo"]) < 0.35 and abs(g["rect"][2] - r["x_hi"]) < 0.35]
        if len(c) < 1 or len(l) != 1: bad.append((p, r["condition"], r["lag"], len(c), len(l))); continue
        want_face = "#ffffff" if r["quiet"] else "#0288d1"
        if not any(g["fill"] == want_face for g in c): bad.append((p, r["condition"], r["lag"], "face"))
        inv = lambda x: math.exp(math.log(XLO) + (x - box[0]) / (box[2] - box[0]) * (math.log(XHI) - math.log(XLO)))
        got = (inv((c[0]["rect"][0] + c[0]["rect"][2]) / 2), inv(l[0]["rect"][0]), inv(l[0]["rect"][2])); exp = (r["hr"], r["lo"], r["hi"])
        dev = max(abs(math.log(g) - math.log(e)) for g, e in zip(got, exp)) / (math.log(XHI) - math.log(XLO)) * 100; worst = max(worst, dev); nread += 3
check(not bad and worst <= 0.5, f"every estimable rung ({nread // 3}) has its marker (filled or open) and its 1.7 pt interval where the record puts them; values read back through the log mapping agree with lag_ladder.json within {worst:.3f}% of the axis span{': ' + str(bad[:4]) if bad else ''}")
xt = sorted({round(g["rect"][0], 2) for g in DW if g["stroke"] == "#1a1d21" and g["width"] is not None and abs(g["width"] - 0.8) < 0.01 and 2.9 < g["rect"][3] - g["rect"][1] < 3.1 and g["rect"][2] - g["rect"][0] < 0.05 and g["rect"][1] < 700})
check(len(xt) == 12 and all(any(abs(x - (DR["panels"][p]["axes_box"][0] + (math.log(t) - math.log(XLO)) / (math.log(XHI) - math.log(XLO)) * (DR["panels"][p]["axes_box"][2] - DR["panels"][p]["axes_box"][0]))) < 0.3 for p in "ab" for t in [0.7, 1, 1.5, 2, 3, 4]) for x in xt), f"12 x tick marks of a and b at the six hazard-ratio ticks of the V13 axes ({len(xt)})")
spine_a = [g for g in DW if g["stroke"] == "#1a1d21" and abs((g["width"] or 0) - 0.8) < 0.01 and g["rect"][3] - g["rect"][1] < 0.05 and abs(g["rect"][0] - DR["panels"]["a"]["axes_box"][0]) < 0.3 and abs(g["rect"][2] - DR["panels"]["a"]["axes_box"][2]) < 0.3]
spine_b = [g for g in DW if g["stroke"] == "#1a1d21" and abs((g["width"] or 0) - 0.8) < 0.01 and g["rect"][3] - g["rect"][1] < 0.05 and abs(g["rect"][0] - DR["panels"]["b"]["axes_box"][0]) < 0.3 and abs(g["rect"][2] - DR["panels"]["b"]["axes_box"][2]) < 0.3]
check(any(abs(g["rect"][1] - DR["panels"]["a"]["axes_box"][3]) < 0.3 for g in spine_a) and any(abs(g["rect"][1] - DR["panels"]["b"]["axes_box"][3]) < 0.3 for g in spine_b), f"bottom spines of a and b at the grown axes feet ({DR['panels']['a']['axes_box'][3]:.2f}, {DR['panels']['b']['axes_box'][3]:.2f}; V13 {G['axes_box']['a'][3]}, {G['axes_box']['b'][3]})")
bands = [g for g in DW if g["fill"] == "#eef0f1" and g["rect"][2] - g["rect"][0] > 100]
check(len(bands) == 9 + 11, f"row bands: {len(bands)} (9 in a for 18 rows, 11 in b for 22 rows around the blank header row; V13 8 + 11)")
c_box = DR["panels"]["c"]["axes_box"]; d_box = DR["panels"]["d"]["axes_box"]
cbars = sorted([g for g in DW if g["fill"] == "#8a9099" and c_box[0] - 1 <= g["rect"][0] and g["rect"][2] <= c_box[2] + 1 and c_box[1] - 1 <= g["rect"][1] and g["rect"][3] <= c_box[3] + 1 and g["rect"][2] - g["rect"][0] > 5], key=lambda g: g["rect"][0])
cvals = [(c_box[3] - g["rect"][1]) / (c_box[3] - c_box[1]) * 142.0 for g in cbars]
check(len(cbars) == 4 and max(abs(v - e) for v, e in zip(cvals, DR["panels"]["c"]["events_deleted_pct"])) < 0.71, f"panel c: 4 grey bars read back to median % events deleted {[round(v, 1) for v in cvals]} (file {DR['panels']['c']['events_deleted_pct']})")
cmk = sorted([g for g in DW if g["fill"] == "#0288d1" and g["stroke"] == "#ffffff" and c_box[0] <= (g["rect"][0] + g["rect"][2]) / 2 <= c_box[2] and c_box[1] <= (g["rect"][1] + g["rect"][3]) / 2 <= c_box[3] and 5 < g["rect"][2] - g["rect"][0] < 7.5], key=lambda g: g["rect"][0])
cy = [(c_box[3] - (g["rect"][1] + g["rect"][3]) / 2) / (c_box[3] - c_box[1]) * 142.0 for g in cmk]
check(len(cmk) == 4 and max(abs(v - e) for v, e in zip(cy, DR["panels"]["c"]["excess_retained_pct"])) < 0.71, f"panel c: 4 line markers read back to median % excess retained {[round(v, 1) for v in cy]} (file {DR['panels']['c']['excess_retained_pct']})")
dbars = sorted([g for g in DW if g["fill"] in ("#0288d1", "#3f9fd8", "#8a9099", "#aeb5bc") and d_box[0] - 1 <= g["rect"][0] and g["rect"][2] <= d_box[2] + 1 and d_box[1] - 1 <= g["rect"][1] and g["rect"][3] <= d_box[3] + 1 and g["rect"][3] - g["rect"][1] > 5], key=lambda g: g["rect"][1])
dvals = [(g["rect"][2] - d_box[0]) / (d_box[2] - d_box[0]) * DR["panels"]["d"]["xlim"][1] for g in dbars]
check(len(dbars) == 4 and max(abs(v - e) for v, e in zip(dvals, [DR["panels"]["d"]["counts"][k] for k in ["Stable", "Attenuating", "Collapsing", "Reversing"]])) < 0.05, f"panel d: 4 bars read back to the trend counts {[round(v, 2) for v in dvals]} (file {DR['panels']['d']['counts']}, axis to {DR['panels']['d']['xlim'][1]:.2f})")
colours = {g["fill"] for g in DW if g["fill"]} | {g["stroke"] for g in DW if g["stroke"]}
check(colours <= {"#1a1d21", "#eef0f1", "#0288d1", "#3f9fd8", "#8a9099", "#aeb5bc", "#ffffff"}, f"vector colours within the V13 set: {sorted(colours)}")
grid = [g for g in DW if g["stroke"] == "#eef0f1"]
check(not grid, f"no pale vertical grid strokes ({len(grid)})")

# 5 positive control: the OLD replot (v7 snapshot, growth 0) against V13
os_, ow, of_, orect, onx, on = spans_of(OLDP); osc = collapse(os_); um = []; md = 0.0
for s in osc:
    c = [b for b in bsc if b["text"] == s["text"] and abs(b["size"] - s["size"]) < 0.06]
    dm = min((max(abs(b["origin"][0] - s["origin"][0]), abs(b["origin"][1] - s["origin"][1])) for b in c), default=99.0)
    if dm >= TOL: um.append((s["text"], s["origin"], round(dm, 3)))
    else: md = max(md, dm)
check(not um and len(osc) == len(bsc) and DRO["growth"]["cd"] == 0.0, f"positive control: the same builder on the v7 snapshot (growth 0) reproduces every V13 string at its origin within {TOL} pt (max {md:.3f}, {len(osc)} vs {len(bsc)}){': ' + str(um[:5]) if um else ''}")
ODW = drawings_of(OLDP); BDW = drawings_of(BASE)
def centres(dw): return sorted(((g["rect"][0] + g["rect"][2]) / 2, (g["rect"][1] + g["rect"][3]) / 2) for g in dw if g["fill"] in ("#0288d1", "#ffffff") and 3 < g["rect"][2] - g["rect"][0] < 12 and 3 < g["rect"][3] - g["rect"][1] < 12)
oc, bc = centres(ODW), centres(BDW)
mdev = max(min(math.hypot(a[0] - b[0], a[1] - b[1]) for b in bc) for a in oc) if oc and bc else 99
check(mdev < 0.5, f"positive control: the OLD replot's {len(oc)} markers sit on V13's {len(bc)} markers within 0.5 pt (max centre deviation {mdev:.3f} pt)")

# 6 multisets: declared delta only
def wm(ss): return collections.Counter(w for t in ss for w in t.split())
def nm(ss): return collections.Counter(x for t in ss for x in re.findall(r"-?\d[\d,]*\.?\d*", t))
data_new = [s["text"] for s in nsc if is_data(s)]; data_base = [s["text"] for s in bsc if is_data(s)]
tick_new = [s["text"] for s in nsc if s["text"].isdigit() and s["origin"][1] > 700 and s["origin"][0] > 600 and not is_d_count(s)]
obs_wa, obs_wr = wm([s["text"] for s in nsc]) - wm([s["text"] for s in bsc]), wm([s["text"] for s in bsc]) - wm([s["text"] for s in nsc])
base_side_w = wm(data_base) + wm(REMOVED_R40)   # round 40: the two key lines leave with the data strings that leave (netted, a word can leave in a note and enter in a label)
exp_wa, exp_wr = wm(data_new) - base_side_w, base_side_w - wm(data_new)
check(obs_wa == exp_wa and obs_wr == exp_wr, f"word multiset: observed delta = the data strings' delta plus the two round-40 key lines (added {sum(exp_wa.values())}, removed {sum(exp_wr.values())}); every other fixed string identical" + ("" if (obs_wa == exp_wa and obs_wr == exp_wr) else f" [added obs-exp {dict(obs_wa - exp_wa)} exp-obs {dict(exp_wa - obs_wa)}; removed obs-exp {dict(obs_wr - exp_wr)} exp-obs {dict(exp_wr - obs_wr)}]"))
obs_na, obs_nr = nm([s["text"] for s in nsc]) - nm([s["text"] for s in bsc]), nm([s["text"] for s in bsc]) - nm([s["text"] for s in nsc])
base_side_n = nm(data_base) + nm(REMOVED_R40)
exp_na, exp_nr = nm(data_new) - base_side_n, base_side_n - nm(data_new)
check(obs_na == exp_na and obs_nr == exp_nr, f"numeric-token multiset: observed delta = the data strings' delta (added {sum(exp_na.values())}, removed {sum(exp_nr.values())})")
json.dump(dict(words_added=dict(exp_wa), words_removed=dict(exp_wr), numbers_added=dict(exp_na), numbers_removed=dict(exp_nr), growth=GR), open(os.path.join(VER, "expected_delta.json"), "w"), indent=1)
def overlaps(sp):
    out = []
    for i, a in enumerate(sp):
        for b in sp[i + 1:]:
            ax0, ay0, ax1, ay1 = a["bbox"]; bx0, by0, bx1, by1 = b["bbox"]
            if ax0 < bx1 - 0.5 and bx0 < ax1 - 0.5 and ay0 < by1 - 1 and by0 < ay1 - 1: out.append((a["text"], b["text"]))
    return out
ov_n, ov_b = overlaps(nsc), overlaps(bsc)
check(len(ov_n) <= len(ov_b), f"no text span overlaps another beyond V13's own ({len(ov_n)} pairs, V13 {len(ov_b)}){': ' + str(ov_n[:4]) if len(ov_n) > len(ov_b) else ''}")
ink = fitz.Rect()
for g in DW: ink |= fitz.Rect(g["rect"])
for s in nsc: ink |= fitz.Rect(s["bbox"])
check(ink.x0 > 2 and ink.y0 > 2 and ink.x1 < nrect[2] - 2 and ink.y1 < nrect[3] - 2, f"all ink inside the page with a 2 pt margin: {[round(v, 2) for v in ink]} on {nrect[2]} x {nrect[3]}")
check(all(os.path.exists(os.path.join(CROPS, f"panel_{p}_OLD_vs_NEW_200dpi.png")) for p in "abcd") and os.path.exists(os.path.join(LANE, "Supp_Fig04_150dpi.png")), "200 dpi OLD (V13) vs NEW crops for the four panels and the 150 dpi render")

# 7 CHANGES csv
rows = []
old_by = {}; new_by = {}
for p in "ab":
    for r in DRO["panels"][p]["rungs"]: old_by[(r["condition"], r["lag"])] = r
    for r in DR["panels"][p]["rungs"]: new_by[(r["condition"], r["lag"])] = r
old_conds = [c for c in DRO["panels"]["a"]["rows"] + DRO["panels"]["b"]["rows"] if c != "Negative controls"]; new_conds = [c for c in DR["panels"]["a"]["rows"] + DR["panels"]["b"]["rows"] if c != "Negative controls"]
for c in old_conds:
    if c not in new_conds: rows.append(dict(sheet="Supp_Fig04", panel="a/b", item=f"{c} (row leaves)", before=f"lag-0 {old_by[(c, 0)].get('hr')}", after="", dramatic="yes", note="control swap of 14 September" if c in ("Back pain", "Cataract") else "lag-0 lower bound no longer above 1 under v8.1"))
for c in new_conds:
    if c not in old_conds: rows.append(dict(sheet="Supp_Fig04", panel="a/b", item=f"{c} (row enters)", before="", after=f"lag-0 {new_by[(c, 0)].get('hr')}", dramatic="yes", note="new outcome or control under v8.1, or lag-0 lower bound now above 1")); continue
    for L in (0, 1, 2, 5):
        o, n = old_by.get((c, L)), new_by.get((c, L))
        if o is None or n is None: continue
        if o["estimable"] != n["estimable"]: rows.append(dict(sheet="Supp_Fig04", panel="a/b", item=f"{c} lag {L} estimable", before=str(o["estimable"]), after=str(n["estimable"]), dramatic="yes", note="")); continue
        if not n["estimable"]: continue
        sig_o, sig_n = (o["lo"] > 1 or o["hi"] < 1), (n["lo"] > 1 or n["hi"] < 1); dram = sig_o != sig_n or abs(n["hr"] - o["hr"]) / o["hr"] > 0.20 or o["quiet"] != n["quiet"]
        if abs(n["hr"] - o["hr"]) > 5e-4 or abs(n["lo"] - o["lo"]) > 5e-4 or abs(n["hi"] - o["hi"]) > 5e-4 or dram:
            rows.append(dict(sheet="Supp_Fig04", panel="a/b", item=f"{c} lag {L} hazard ratio (95% CI)", before=f"{o['hr']:.3f} ({o['lo']:.3f}-{o['hi']:.3f}){' open' if o['quiet'] else ''}", after=f"{n['hr']:.3f} ({n['lo']:.3f}-{n['hi']:.3f}){' open' if n['quiet'] else ''}", dramatic="yes" if dram else "no", note=("CI crosses 1 changed" if sig_o != sig_n else "") + (" open/filled changed" if o["quiet"] != n["quiet"] else "")))
for dl in DR.get("display_labels", []): rows.append(dict(sheet="Supp_Fig04", panel=dl["panel"], item="display label", before="", after=f"{dl['display']} for {dl['condition']}", dramatic="no", note=f"full label {dl['width_pt']} pt wide, gutter {dl['gutter_pt']} pt: the Main_Fig2 display label, owner question"))
rows.append(dict(sheet="Supp_Fig04", panel="a/b", item="rows (conditions with lag-0 lower bound above 1)", before=f"{JO['panel_n']} (16 + 16)", after=f"{J['panel_n']} ({DR['half']} + {J['panel_n'] - DR['half']})", dramatic="yes" if abs(J["panel_n"] - JO["panel_n"]) / JO["panel_n"] > 0.2 else "no", note=f"page grows {GR['cd']} pt, pitch {GR['pitch']} kept"))
for i, L in enumerate([0, 1, 2, 5]):
    rows.append(dict(sheet="Supp_Fig04", panel="c", item=f"landmark {L}: conditions surviving / panel", before=DRO["panels"]["c"]["annotations"][i]["text"], after=DR["panels"]["c"]["annotations"][i]["text"], dramatic="no", note=f"bar {DRO['panels']['c']['events_deleted_pct'][i]} -> {DR['panels']['c']['events_deleted_pct'][i]}, line {DRO['panels']['c']['excess_retained_pct'][i]} -> {DR['panels']['c']['excess_retained_pct'][i]}"))
for k in ["Stable", "Attenuating", "Collapsing", "Reversing"]:
    o, n = DRO["panels"]["d"]["counts"][k], DR["panels"]["d"]["counts"][k]
    rows.append(dict(sheet="Supp_Fig04", panel="d", item=f"{k} conditions", before=str(o), after=str(n), dramatic="yes" if o and abs(n - o) / o > 0.2 else "no", note=""))
write_changes(os.path.join(VER, "CHANGES_Supp_Fig04.csv"), rows); info(f"CHANGES_Supp_Fig04.csv: {len(rows)} rows, {sum(r['dramatic'] == 'yes' for r in rows)} dramatic")

# 8 printed values csv
recs = []
for s in nsc:
    t = s["text"]; base = dict(text=t, x=s["bbox"][0], baseline=s["origin"][1], ha="left", size=s["size"])
    if is_row_label(s):
        cond, ev = t.rsplit(" (", 1); ev = ev[:-1].replace(",", "")
        recs.append(dict(base, source_file=DATA_CSV, source_key=f"condition={cond}&lag_years=5.0:events", source_value=int(ev), rule="text", panel="a" if s["origin"][0] < 500 else "b", note="label = condition (events surviving the 5-year landmark)"))
        recs[-1]["rule"] = "text"; recs[-1]["source_value"] = t; recs[-1]["source_key"] = f"condition {cond}: condition and events at lag_years 5.0 ({ev})"
    elif is_c_ann(s):
        i = [x["text"] for x in DR["panels"]["c"]["annotations"]].index(t); L = ["0", "1", "2", "5"][i]
        recs.append(dict(base, source_file=DATA, source_key=("panel_n / panel_n" if L == "0" else f"by_lag.{L}.survive / panel_n"), source_value=t, rule="text", panel="c"))
    elif is_d_count(s): recs.append(dict(base, source_file=DATA, source_key=f"trends: count of '{[x['label'] for x in DR['panels']['d']['count_texts'] if x['text'] == t][0].lower()}' among the panel conditions", source_value=t, rule="text", panel="d"))
    elif is_noest(s): recs.append(dict(base, source_file=DATA, source_key="per_condition.<key>.hr_by_lag.<lag> is null", source_value=t, rule="text", panel="a" if s["origin"][0] < 500 else "b"))
    elif is_ctrl_head(s): recs.append(dict(base, source_file="static:V13", source_key="group header", source_value=t, rule="text", panel="b"))
    else: recs.append(dict(base, source_file="static:V13", source_key="V13 string (region shift declared)", source_value=t, rule="text", panel=""))
n_rows, n_no = printed_values_csv(NEW, recs, os.path.join(VER, "Supp_Fig04_printed_values.csv"), static_from=BASE)
check(n_no == 0, f"printed values: {n_rows} strings accounted for, {n_no} NO")
ok = ck.write(os.path.join(VER, "checks.txt")); print("RESULT ALL PASS" if ok else "RESULT FAIL")
