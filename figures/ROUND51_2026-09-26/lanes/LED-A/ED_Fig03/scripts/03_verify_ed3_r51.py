#!/usr/bin/env python3
"""ROUND 51, lane LED-A: 03_verify for ED_Fig03 with the larger boxes (item 1). Adapted from the round-49 verifier (03_verify_ed3_r49.py), run
on the built sheet before the title strip (work/ED_Fig03_built.pdf). Declared this round: the 20 axes boxes 85.16 x 60.83 pt (V28 72.3 x 48.66,
+17.8 and +25 percent), the horizontal gutter 33.6 pt (the widest y tick label plus 6.1 pt plus 2.5 mm), the row gutter 56 pt (V28 75.7),
the page 604.4 x 575.2 pt before the strip (V28 518.74 x 565.92), everything at its V28 text size with the V28 line breaks. Kept numbers:
every marker and ribbon value read back through each sheet's own y ticks equals V28's, every string equals V28's (the census in run_sheet.py).
The coordinator's rule: at least 2.5 mm between any panel's y tick labels and the box to its left (checked here). The data checks (rows by the
August rule, read-back against results_v2.json, asterisks by BH q, tick rule) are the round-37 checks."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import fitz, re, os, sys, gc, csv, collections, time, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ed3_common import *
from ed3_extract_r49 import extract, spans_of
from v14lib import printed_values_csv, Checks
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND51_2026-09-26/lanes/LED-A/scripts_shared")
import r49lib as R

NEW = f"{WORK}/ED_Fig03_built.pdf"; DRAWN = f"{VER}/ED_Fig03_drawn.json"; V28 = f"{R.V28}/ED_Fig03.pdf"; STRIP = 16.0
INK = (0x1a / 255, 0x1d / 255, 0x21 / 255); INK_INT = 0x1a1d21
LADDER = [(0xb3 / 255, 0xdc / 255, 0xf2 / 255), (0x7c / 255, 0xc0 / 255, 0xe9 / 255), (0x3f / 255, 0x9f / 255, 0xd8 / 255), (0x02 / 255, 0x88 / 255, 0xd1 / 255)]
HEX = ["#b3dcf2", "#7cc0e9", "#3f9fd8", "#0288d1"]; BAND_LABELS = {"0-1": "0-1%", "1-5": "1-5%", "5-10": "5-10%", ">10": ">10%"}
TOL = 0.3; MM = 72.0 / 25.4
for p in (BASE, NEW, DATA, PANEL_CSV, OLD_DATA, OLD_PANEL, DRAWN, V28): hydrated(p)
drawn = load_json(DRAWN); meta = drawn["meta"]; DELTA = float(drawn["delta_pt"]); KEEP = set(drawn["keep_size"]); R51 = drawn["round51"]; assert DELTA == 1.0 and not drawn["key2_drawn"]
TS, TK, KEYS = 9.5 + DELTA, 9.5 + DELTA, 9.0 + DELTA
KEY2_WORDS = "*** top band, q below 0.001"
ck = Checks(f"ED_Fig03 round 51 lane LED-A ({time.strftime('%Y-%m-%d %H:%M')}), Proof W with the larger boxes: {NEW} (sha256 {sha256(NEW)[:16]}) against V28 {V28} (sha256 {sha256(V28)[:16]}), the V13 design sheet "
            f"and numbers/results_v2.json graded + negcontrols_v4_panel.csv. Declared: boxes {R51['axw']:.2f} x {R51['axh']:.2f} pt (V28 {R51['axw_v28']:.2f} x {R51['axh_v28']:.2f}), gutter {R51['gap']:.2f} pt (V28 {R51['gap_v28']:.2f}), "
            f"row gutter {R51['rowgap']:.1f} pt (V28 {R51['rowgap_v28']:.1f}), page {R51['page'][0]:.2f} x {R51['page'][1]:.2f} pt before the strip (V28 {R51['page_v28_prestrip'][0]:.2f} x {R51['page_v28_prestrip'][1]:.2f}); text sizes, line breaks, y axes and key as V28.")
def check(ok, msg): return ck.log(msg, bool(ok))
def info(msg): ck.lines.append("INFO: " + msg); print("INFO: " + msg[:300], flush=True)
def near(c, t, tol=0.003): return c is not None and max(abs(a - b) for a, b in zip(c, t)) < tol

cur = extract(BASE); new = extract(NEW, title_sizes=(TS, 9.5), tick_size=TK); v28 = extract(V28, title_sizes=(TS, 9.5), tick_size=TK)
out = dict(graded(DATA)["outcomes"]); outo = dict(graded(OLD_DATA)["outcomes"]); Gn, Go = graded(DATA), graded(OLD_DATA)
check(sha256(DATA) == drawn["sources"][DATA]["sha256"] and sha256(PANEL_CSV) == drawn["sources"][PANEL_CSV]["sha256"], "the numbers files read now have the sha256 the builder recorded (v8.1 sidecars checked by the builder)")
check(drawn["sources"][DATA]["v8_inputs"] and drawn["sources"][PANEL_CSV]["v8_inputs"], f"sidecars name data_frozen_v8_2026-09: results_v2.json step {drawn['sources'][DATA]['step']} ({drawn['sources'][DATA]['output_mtime']}), negcontrols_v4_panel.csv step {drawn['sources'][PANEL_CSV]['step']} ({drawn['sources'][PANEL_CSV]['output_mtime']})")

# 1. page, fonts, sizes
W, H = new["page"]
check(abs(W - R51["page"][0]) < 0.01 and abs(H - R51["page"][1]) < 0.01, f"page (before the strip) {W:.2f} x {H:.2f} pt as the builder declares (V28 {R51['page_v28_prestrip'][0]:.2f} x {R51['page_v28_prestrip'][1]:.2f}: +{W - R51['page_v28_prestrip'][0]:.1f} pt wide, +{H - R51['page_v28_prestrip'][1]:.1f} pt tall)")
d = fitz.open(NEW); pg = d[0]
check(len(d) == 1 and len(pg.get_xobjects()) == 0, f"single flat page: {len(d)} page, {len(pg.get_xobjects())} XObjects")
fitz.TOOLS.mupdf_warnings(reset=True); pg.get_pixmap(dpi=50); w = fitz.TOOLS.mupdf_warnings(reset=True)
check(not w, "MuPDF renders without warnings" + (": " + w if w else ""))
fn = pg.get_fonts(full=True); names_new = sorted(f[3].split("+")[-1] for f in fn); d.close()
d = fitz.open(V28); fb = d[0].get_fonts(full=True); names_v28 = sorted(f[3].split("+")[-1] for f in fb); d.close()
check(set(names_new) <= set(names_v28) and all("Arial" in n for n in names_new), f"fonts: {names_new} (V28 {names_v28}, which also carries the strip's Arial Bold)")
check(all(f[2] == "TrueType" and f[5] == "WinAnsiEncoding" for f in fn), "fonts embedded as simple TrueType WinAnsi, as the V13 sheet's")
sizes_new = sorted({s["size"] for s in new["spans"]}); sizes_v28 = sorted({s["size"] for s in v28["spans"] if not s["text"].startswith("Extended")})
check(sizes_new == sizes_v28, f"type sizes {sizes_new} = V28's {sizes_v28} (round-49 values, the two kept 9.5 pt titles included)")

# 2. letters at the builder's positions
def letters(sp, z): return {s["text"]: s for s in sp if s["text"] in ("a", "b") and abs(s["size"] - z) < 0.01}
ln = letters(new["spans"], 13.0 + DELTA)
for L in "ab":
    ok = L in ln and "Bold" in ln[L]["font"] and abs(ln[L]["origin"][0] - drawn["letters"][L][0]) < TOL and abs(ln[L]["origin"][1] - drawn["letters"][L][1]) < TOL
    check(ok, f"letter {L}: {13 + DELTA:g} pt {ln.get(L, {}).get('font')} at origin {ln.get(L, {}).get('origin')} (builder {drawn['letters'][L]})")

# 3. row set, the new geometry, the growth against V28
rows_new = [p["title"] for p in new["panels"]]; rows_v13 = [p["title"] for p in cur["panels"]]; rows_v28 = [p["title"] for p in v28["panels"]]
rule_a, rule_b = august_rows(out, controls_of(PANEL_CSV))
check(len(new["panels"]) == 20 and len(v28["panels"]) == 20, f"20 panels ({len(new['panels'])} new, {len(v28['panels'])} V28)")
check(rows_new == rule_a + rule_b == meta["rows_drawn"] == rows_v28, "row set and order = the August rule on the v8.1 graded block = V28: " + " | ".join(rows_new))
dev_box = max(abs(pn["axes"][k] - v) for pn, rec in zip(new["panels"], drawn["panels"]) for k, v in zip(("left", "top", "right", "bottom"), rec["box"]))
check(dev_box < 0.05, f"axes boxes (spines) at the builder's positions within 0.05 pt: max deviation {dev_box:.3f} pt")
wr = [(pn["axes"]["right"] - pn["axes"]["left"]) / (pv["axes"]["right"] - pv["axes"]["left"]) for pn, pv in zip(new["panels"], v28["panels"])]
hr_ = [(pn["axes"]["bottom"] - pn["axes"]["top"]) / (pv["axes"]["bottom"] - pv["axes"]["top"]) for pn, pv in zip(new["panels"], v28["panels"])]
check(all(abs(r - R51["scale_w"]) < 0.006 for r in wr) and all(abs(r - R51["scale"]) < 0.006 for r in hr_), f"every box {min(wr):.3f} to {max(wr):.3f} times V28's width and {min(hr_):.3f} to {max(hr_):.3f} times its height (declared {R51['scale_w']:.3f} and {R51['scale']:.2f}; the extractor reads spines to 0.1 pt)")
gaps_h = [new["panels"][i + 1]["axes"]["left"] - new["panels"][i]["axes"]["right"] for r in range(4) for i in range(r * 5, r * 5 + 4)]
check(max(abs(g - R51["gap"]) for g in gaps_h) < 0.15, f"horizontal gutter {min(gaps_h):.2f} to {max(gaps_h):.2f} pt on every row (declared {R51['gap']:.2f}, V28 {R51['gap_v28']:.2f}; the widest y tick label {R51['ytick_w_max']:.2f} pt + 6.1 + 2.5 mm)")
rg = [new["panels"][r * 5 + 5]["axes"]["top"] - new["panels"][r * 5]["axes"]["bottom"] for r in range(2)]
check(max(abs(g - R51["rowgap"]) for g in rg) < 0.15, f"row gutter {min(rg):.2f} to {max(rg):.2f} pt between the panel a rows (declared {R51['rowgap']:.0f}, V28 {R51['rowgap_v28']:.1f})")
dev_xt = max(abs(a - (pn["axes"]["left"] + (pn["axes"]["right"] - pn["axes"]["left"]) / 4.0 * (k + 0.5))) for pn in new["panels"] for k, a in enumerate(pn["x_ticks"])) if all(len(pn["x_ticks"]) == 4 for pn in new["panels"]) else 99
check(dev_xt < TOL, f"x tick marks: 4 per panel at the band centres within {TOL} pt (max deviation {dev_xt:.3f} pt)")
# the coordinator's rule: y tick labels at least 2.5 mm from the box to their left
clear = []
for i, pn in enumerate(new["panels"]):
    col = i % 5; L = pn["axes"]["left"]
    labs = [s for s in new["spans"] if abs(s["size"] - TK) < 0.06 and "Bold" not in s["font"] and re.fullmatch(r"[\d.]+", s["text"]) and abs(s["bbox"][2] - (L - 7.32 / 1.2)) < 1.5 and pn["axes"]["top"] - 6 <= s["origin"][1] <= pn["axes"]["bottom"] + 6]
    if col == 0 or not labs: continue
    left_box = new["panels"][i - 1]["axes"]["right"]
    clear.append((pn["title"], min(s["bbox"][0] for s in labs) - left_box))
least = min(c for _, c in clear)
check(least >= 2.5 * MM - 0.05, f"every panel's y tick labels at least 2.5 mm from the box to its left: least clearance {least / MM:.2f} mm ({least:.2f} pt, {min(clear, key=lambda c: c[1])[0]}), {len(clear)} panels checked")
col0 = [s for s in new["spans"] if abs(s["size"] - TK) < 0.06 and re.fullmatch(r"[\d.]+", s["text"]) and abs(s["bbox"][2] - (new["panels"][0]["axes"]["left"] - 7.32 / 1.2)) < 1.5]
ylab = [s for s in new["spans"] if s["text"].startswith("Hazard ratio")]
info(f"first column: its y tick labels clear the rotated y label by {(min(s['bbox'][0] for s in col0) - max(s['bbox'][2] for s in ylab)) / MM:.2f} mm")
dr_t = drawn["text"]
def at(t, s): return abs(t["baseline"] - s["origin"][1]) < 0.35 and s["bbox"][0] - 0.35 <= t["x"] <= s["bbox"][2]
def rec_of(s):
    c = [t for t in dr_t if t["text"] == s["text"] and abs(t["x"] - s["origin"][0]) < 0.35 and at(t, s)]
    return c[0] if c else None
title_words = set(w for t in rows_new + rows_v13 for w in t.split())
def is_title(s): return abs(s["size"] - 9.5) < 0.06 and "Bold" not in s["font"] and bool(s["text"].split()) and s["text"].split()[0] in title_words and not re.fullmatch(r"[\d.]+", s["text"]) and s["text"] not in ("0-1", "1-5 5-10 >10", "1-5", "5-10", ">10")
def is_title_new(s): return (abs(s["size"] - TS) < 0.06 or abs(s["size"] - 9.5) < 0.06) and "Bold" not in s["font"] and bool(s["text"].split()) and s["text"].split()[0] in title_words and not re.fullmatch(r"[\d.]+", s["text"]) and s["text"] not in ("0-1", "1-5 5-10 >10", "1-5", "5-10", ">10")
fx_c = [s for s in cur["spans"] if not re.fullmatch(r"[\d.]+", s["text"]) and "*" not in s["text"] and not is_title(s) and s["text"] != "top band, q below 0.001"]
fx_n = [s for s in new["spans"] if not re.fullmatch(r"[\d.]+", s["text"]) and "*" not in s["text"] and not is_title_new(s)]
cx = R51["cx"]
def centred(text):
    sp = [s for s in new["spans"] if s["text"] == text]; return sp and abs((sp[0]["bbox"][0] + sp[0]["bbox"][2]) / 2 - cx) < 0.5
check(centred("Time with oxygen saturation below 90%, % of the recording") and centred("Negative controls"), f"the x title and 'Negative controls' centred on the panels' centre x {cx:.2f} (within 0.5 pt)")
yl = [s for s in new["spans"] if s["text"].startswith("Hazard ratio")]; bm = (new["panels"][0]["axes"]["top"] + new["panels"][15]["axes"]["bottom"]) / 2.0
check(len(yl) == 1 and abs((yl[0]["bbox"][1] + yl[0]["bbox"][3]) / 2 - bm) < 1.0 and abs(yl[0]["origin"][0] - 16.0) < TOL, f"the rotated y label at x 16.0, centred on the block of panels (mid {bm:.1f})")
key_n = [s for s in fx_n if (rec_of(s) or {}).get("what") == "key"]; key_c = [c for c in fx_c if 530 < c["origin"][1] < 548]
kw_n = collections.Counter(w for s in key_n for w in s["text"].split()); kw_c = collections.Counter(w for c in key_c for w in c["text"].split())
row_new = (min(s["bbox"][0] for s in key_n), max(s["bbox"][2] for s in key_n)); row_ink = (row_new[0] - 15.2, row_new[1])
sw_x = sorted(r_["x"] for r_ in dr_t if r_["what"] == "key" and r_["key"].startswith("band label")); wd_x = sorted(r_["x"] for r_ in dr_t if r_["what"] == "key" and r_["key"].startswith("band wording"))
check(len(key_n) == 8 and kw_n == kw_c and all(abs(s["size"] - KEYS) < 0.06 for s in key_n) and abs((row_ink[0] + row_ink[1]) / 2 - cx) < 0.5 and all(abs((b - a) - drawn["key_mean_dx"] + 15.2) < 0.05 for a, b in zip(sw_x, wd_x)),
      f"key row: the V13 key words at {KEYS:g} pt, wording {drawn['key_mean_dx'] - 15.2:.1f} pt right of its band label, ink {row_ink[0]:.1f} to {row_ink[1]:.1f} pt centred on x {cx:.2f}")
check(row_ink[0] > 2.0 and row_ink[1] < W - 2.0, f"the key row stays inside the page ({row_ink[0]:.1f} to {row_ink[1]:.1f} on {W:.2f} pt)")
title_bad = []
for p, rec in zip(new["panels"], drawn["panels"]):
    L, T = p["axes"]["left"], p["axes"]["top"]; lines = rec["title_lines"]; tsz = rec["title_size"]; pitch = rec["title_pitch"]
    for j, line in enumerate(lines):
        base_y = T - 4.0 - (len(lines) - 1 - j) * pitch
        sp = [s for s in new["spans"] if s["text"] == line and abs(s["size"] - tsz) < 0.06 and abs(s["origin"][0] - L) < TOL and abs(s["origin"][1] - base_y) < TOL]
        if not sp: title_bad.append((rec["title"], line))
    if (rec["title"] in KEEP) != (abs(tsz - 9.5) < 0.01): title_bad.append((rec["title"], "size", tsz))
same_breaks = [r["title_lines"] for r in drawn["panels"]] == [pv["title"] and r["title_lines"] for r, pv in zip(drawn["panels"], v28["panels"])]
v28_lines = {}
for pv in v28["panels"]:
    L = pv["axes"]["left"]; v28_lines[pv["title"]] = [s["text"] for s in sorted([s for s in v28["spans"] if abs(s["origin"][0] - L) < 0.6 and pv["axes"]["top"] - 40 < s["origin"][1] < pv["axes"]["top"] and (abs(s["size"] - TS) < 0.06 or abs(s["size"] - 9.5) < 0.06) and "Bold" not in s["font"]], key=lambda s: s["origin"][1])]
same_breaks = all(r["title_lines"] == v28_lines.get(r["title"]) for r in drawn["panels"])
check(not title_bad and same_breaks, f"titles: every line at its panel's left edge, 4 pt above the box, at the round-49 size and pitch, with the V28 line breaks (wrap limit kept at {R51['wrap_limit']:.2f} pt); wrapped {[r['title'] for r in drawn['panels'] if len(r['title_lines']) > 1]}{': ' + str(title_bad[:3]) if title_bad else ''}")
ytl_bad = []; NUD = {(n["title"], n["label"]): n for n in drawn.get("ytick_nudges", [])}
for pn in new["panels"]:
    L = pn["axes"]["left"]
    for y, lab in pn["y_ticks"]:
        if lab is None: ytl_bad.append((pn["title"], y, "no label")); continue
        n_ = NUD.get((pn["title"], lab))
        sp = [s for s in new["spans"] if s["text"] == lab and abs(s["size"] - TK) < 0.06 and (abs(s["bbox"][2] - (L - 7.32 / 1.2)) < TOL if n_ is None else abs(s["bbox"][0] - n_["x"]) < TOL) and abs(s["origin"][1] - 0.377 * TK - y) < TOL]
        if not sp: ytl_bad.append((pn["title"], y, lab))
check(not ytl_bad and not NUD, f"y tick labels {TK:g} pt, right-aligned 6.1 pt left of the spine and centred on their tick marks within {TOL} pt, no label nudged (the gutter holds the widest)")

# 4. style
band_ok = rib_ok = line_ok = mark_ok = ref_ok = star_ok = True; band_cols = set()
for p in new["panels"]:
    bf = p["band_fills"]
    if len(bf) != 4 or any(not near(bf[k][2], LADDER[k]) for k in range(4)): band_ok = False
    for b in bf: band_cols.add(b[2])
    if len(p["ribbons"]) != 1 or not near(p["ribbons"][0]["fill"], INK) or abs(p["ribbons"][0]["opacity"] - 0.16) > 0.005: rib_ok = False
    if len(p["line"]) != 1 or not near(p["line"][0]["color"], (1, 1, 1)) or len(p["line"][0]["pts"]) != 4: line_ok = False
    if len(p["markers"]) != 4 or any(not near(m["fill"], INK) or not near(m["stroke"], (1, 1, 1)) or abs(m["sw"] - 0.5) > 0.01 or abs(m["r"] - 2.0) > 0.05 for m in p["markers"]): mark_ok = False
    if len(p["ref"]) != 1 or not near(p["ref"][0]["color"], INK) or abs(p["ref"][0]["width"] - 0.9) > 0.01 or abs(p["ref_value"][0] - 1.0) > 0.0005 * (p["ymax"] - p["ymin"]): ref_ok = False
    if p["star"] and p["star_color"] != INK_INT: star_ok = False
hexes = sorted("#%02x%02x%02x" % tuple(round(c * 255) for c in col) for col in band_cols)
check(band_ok and hexes == sorted(HEX), f"band fills: four equal rects per panel, exactly {hexes} (the V13 four)")
check(rib_ok and line_ok and mark_ok and ref_ok and star_ok, "ribbon ink 16%, white 1.5 pt line through 4 points, 4 ink markers r 2.0 with 0.5 pt white ring, dashed ink reference at 1, ink asterisks (marker sizes and line weights unchanged)")

# 5. values: read back vs the numbers file, and identical to V28's read-back
rows = []; maxdev_json = maxdev_rib = maxdev_old = 0.0
for pn in new["panels"]:
    v = out[pn["title"]]; span = pn["ymax"] - pn["ymin"]
    hr_json = [1.0] + [v[k]["hr"] for k in BAND_KEYS[1:]]; dev_j = max(abs(a - b) for a, b in zip(pn["marker_hr"], hr_json)) / span * 100
    lo_hi = {}
    for x, val in pn["ribbon_values"][0]: lo_hi.setdefault(x, []).append(val)
    dev_r = 0.0
    for x, k in zip(sorted(lo_hi), BAND_KEYS[1:]):
        vs = sorted(lo_hi[x]); dev_r = max(dev_r, abs(vs[0] - v[k]["lo"]) / span * 100, abs(vs[-1] - v[k]["hi"]) / span * 100)
    maxdev_json = max(maxdev_json, dev_j); maxdev_rib = max(maxdev_rib, dev_r); rows.append((pn["title"], hr_json, pn["marker_hr"], dev_j, dev_r, pn["star"]))
for pc in cur["panels"]:
    vo = outo[pc["title"]]; hr_old = [1.0] + [vo[k]["hr"] for k in BAND_KEYS[1:]]
    maxdev_old = max(maxdev_old, max(abs(a - b) for a, b in zip(pc["marker_hr"], hr_old)) / (pc["ymax"] - pc["ymin"]) * 100)
check(maxdev_json <= 0.5, f"NEW sheet marker centres read back vs the v8.1 file: max deviation {maxdev_json:.3f}% of the axis span (limit 0.5%)")
check(maxdev_rib <= 0.5, f"NEW sheet ribbon vertices (95% CI) read back vs the v8.1 file: max deviation {maxdev_rib:.3f}% of span")
check(maxdev_old <= 0.5, f"positive control of the extractor: the V13 sheet's markers read back vs the v7 snapshot: max deviation {maxdev_old:.3f}% of span")
dv = max(abs(a - b) / (pn["ymax"] - pn["ymin"]) * 100 for pn, pv in zip(new["panels"], v28["panels"]) for a, b in zip(pn["marker_hr"], pv["marker_hr"]))
dl = max((abs(pn["ymin"] - pv["ymin"]) + abs(pn["ymax"] - pv["ymax"])) / (pn["ymax"] - pn["ymin"]) * 100 for pn, pv in zip(new["panels"], v28["panels"]))
dt = all(sorted(l for _, l in pn["y_ticks"] if l) == sorted(l for _, l in pv["y_ticks"] if l) for pn, pv in zip(new["panels"], v28["panels"]))
check(dv <= 0.1 and dl < 0.5 and dt, f"kept numbers: every marker value, y range and y tick set equals V28's read through each sheet's own ticks (max marker deviation {dv:.3f}% of span, y ranges within {dl:.2f}% of span, the extractor's 0.1 pt spine rounding)")
with open(f"{VER}/values_agreement.txt", "w") as f:
    f.write(f"ED_Fig03 round 51: plotted values read back from the vectors through each panel's y ticks\nnew: {NEW}\nV28: {V28}\nv8.1: {DATA} graded.outcomes\nold: {OLD_DATA}\n\n")
    for r in rows: f.write(f"{r[0]:42s} file {' '.join(f'{v:.3f}' for v in r[1]):28s} read {' '.join(f'{v:.3f}' for v in r[2]):28s} dev {r[3]:6.3f}% CI dev {r[4]:6.3f}%  star {r[5] or '-'}\n")

# 6. stars and the (absent) key sentence
q_all = bh_top(out); q_allo = bh_top(outo)
star_rule = ["***" if q_all[t] < 0.001 else "" for t in rows_new]
check([p["star"] for p in new["panels"]] == star_rule, f"asterisks = the rule (BH q of the top-band contrast across the {len(out)} top-band contrasts, *** below 0.001): " + "".join("*" if s else "-" for s in star_rule))
key2 = [s for s in new["spans"] if s["text"] == "top band, q below 0.001"]; all15 = all(q_all[t] < 0.001 for t in rows_new[:15])
check(len(key2) == 0 and all15, f"key sentence '{KEY2_WORDS}' absent (round 40); its rule still holds: every panel-a outcome has q < 0.001 (largest {max(q_all[t] for t in rows_new[:15]):.3g})")
check([p["star"] for p in new["panels"]] == [p["star"] for p in v28["panels"]], "asterisk pattern per slot identical to V28 (15 ***, 5 none)")

# 7. numeric tokens: every token matched to a source; the expected delta against V13
def tokens(sp):
    out_ = []
    for s in sp:
        for t in s["text"].split():
            if re.search(r"\d", t): out_.append((t, s))
    return out_
rb = []; unmatched_tok = []
for t, s in tokens(new["spans"]):
    x, y = s["origin"]; srcf = key = val = None
    if is_title_new(s): srcf, key, val = DATA, f"condition name graded.outcomes['{[r for r in rows_new if s['text'] in r][0] if any(s['text'] in r for r in rows_new) else s['text']}']", t
    elif re.fullmatch(r"[\d.]+", t) and abs(s["size"] - TK) < 0.06:
        pn = [p for p in new["panels"] if abs(s["bbox"][2] - (p["axes"]["left"] - 7.32 / 1.2)) < TOL and p["axes"]["top"] - 6 <= y <= p["axes"]["bottom"] + 6]
        if len(pn) == 1:
            ymin, ymax, ticks = ylim_ticks(out[pn[0]["title"]])
            if t in ticks: srcf, key, val = DATA, f"tick rule on ({ymin:.4f}, {ymax:.4f}) from graded.outcomes['{pn[0]['title']}'] lo/hi", float(t)
    elif t in BAND_LABELS and (abs(s["size"] - KEYS) < 0.06 or abs(s["size"] - TK) < 0.06): srcf, key, val = DATA, f"graded.band key '{BAND_LABELS[t]}'", BAND_LABELS[t]
    elif s["text"] in ("1% or less", "over 1 to 5%", "over 5 to 10%", "over 10%"): srcf, key, val = DATA, f"band wording for {BAND_KEYS[['1% or less', 'over 1 to 5%', 'over 5 to 10%', 'over 10%'].index(s['text'])]}", s["text"]
    elif s["text"].startswith("Hazard ratio vs"): srcf, key, val = DATA, "y label: reference band 0-1% and the 95% CI", t
    elif s["text"].startswith("Time with oxygen"): srcf, key, val = DATA, "x title: the 90% saturation threshold of T90", t
    if key is None: unmatched_tok.append((t, s["text"], s["origin"]))
    rb.append(dict(token=t, span=s["text"], x=x, y=y, size=s["size"], source=srcf, key=key, value=val))
with open(f"{VER}/readback_tokens.csv", "w", newline="") as f:
    wtr = csv.DictWriter(f, fieldnames=list(rb[0])); wtr.writeheader(); wtr.writerows(rb)
check(not unmatched_tok, f"READ-BACK: all {len(rb)} numeric tokens on the text layer matched to a source key (verify/readback_tokens.csv){': unmatched ' + str(unmatched_tok[:5]) if unmatched_tok else ''}")
tick_ok = all(sorted({l for _, l in p["y_ticks"] if l}, key=float) == sorted(ylim_ticks(out[p["title"]])[2], key=float) for p in new["panels"])
check(tick_ok, "y tick label set of every panel = the tick rule on the v8.1 limits: " + " | ".join(",".join(sorted({l for _, l in p['y_ticks'] if l}, key=float)) for p in new["panels"]))
tn, tc = collections.Counter(t for t, _ in tokens(new["spans"])), collections.Counter(t for t, _ in tokens(cur["spans"]))
old_ticks = collections.Counter(); new_ticks = collections.Counter()
for pc in cur["panels"]: old_ticks.update(ylim_ticks(outo[pc["title"]])[2])
for p in new["panels"]: new_ticks.update(ylim_ticks(out[p["title"]])[2])
title_tok_new = collections.Counter(t for r in rows_new for t in r.split() if re.search(r"\d", t)); title_tok_old = collections.Counter(t for r in rows_v13 for t in r.split() if re.search(r"\d", t))
exp_add, exp_rem = (new_ticks + title_tok_new) - (old_ticks + title_tok_old), (old_ticks + title_tok_old) - (new_ticks + title_tok_new) + collections.Counter({"0.001": 1})
got_add, got_rem = tn - tc, tc - tn
check(got_add == exp_add and got_rem == exp_rem, f"expected delta against V13: numeric tokens added {dict(got_add)} removed {dict(got_rem)} = tick rule on the v8.1 rows minus the V13 rows plus the title digits, plus the key sentence's 0.001 (round 40)")

# 8. word multiset: identical apart from the titles, the y tick labels and the key sentence (declared)
def collapse_(sp):
    keep = []
    for s in sp:
        if any(k["text"] == s["text"] and abs(k["size"] - s["size"]) < 0.06 and abs(k["origin"][0] - s["origin"][0]) < 0.25 and abs(k["origin"][1] - s["origin"][1]) < 0.25 for k in keep): continue
        keep.append(s)
    return keep
def words(sp):
    c = collections.Counter()
    for s in collapse_(sp): c.update(s["text"].split())
    return c
wn, wc_ = words(new["spans"]), words(cur["spans"])
title_w_new = collections.Counter(w for r in rows_new for w in r.split()); title_w_old = collections.Counter(w for r in rows_v13 for w in r.split())
wn_f = collections.Counter({k: v for k, v in wn.items() if not re.fullmatch(r"[\d.]+", k)}) - title_w_new
wc_f = collections.Counter({k: v for k, v in wc_.items() if not re.fullmatch(r"[\d.]+", k)}) - title_w_old - collections.Counter(w for w in KEY2_WORDS.split() if not re.fullmatch(r"[\d.]+", w))
check(wn_f == wc_f, f"word multiset identical to V13 apart from the titles, the y tick labels and the removed key sentence (declared): {sum(wn.values())} words")

# 9. overlaps, bounds, crowding, text layer vs drawn
d = fitz.open(NEW); pg = d[0]; sp = new["spans"]; boxes = [p["axes"] for p in new["panels"]]
def body(s):
    r = fitz.Rect(s["bbox"]); return fitz.Rect(r.x0, r.y0 + 0.15 * s["size"], r.x1, r.y1 - 0.2 * s["size"])
ov = [(sp[i]["text"], sp[j]["text"]) for i in range(len(sp)) for j in range(i + 1, len(sp)) if not (body(sp[i]) & body(sp[j])).is_empty and (body(sp[i]) & body(sp[j])).get_area() > 0.05]
check(not ov, f"no text span overlaps another ({len(ov)}{': ' + str(ov[:3]) if ov else ''})")
intr = []
for s in sp:
    r = fitz.Rect(s["bbox"]); r = fitz.Rect(r.x0, r.y0 + 0.2 * s["size"], r.x1, r.y1 - 0.2 * s["size"])
    for b in boxes:
        o = r & fitz.Rect(b["left"], b["top"], b["right"], b["bottom"])
        if not o.is_empty and o.get_area() > 0.05 and "*" not in s["text"]: intr.append((s["text"], b["left"], b["top"]))
check(not intr, f"no tick label or title runs into any axes box ({len(intr)}{': ' + str(intr[:3]) if intr else ''})")
ink_r = fitz.Rect()
for g in pg.get_drawings(): ink_r |= g["rect"]
for s in sp: ink_r |= fitz.Rect(s["bbox"])
check(ink_r.x0 > 2 and ink_r.y0 > 2 and ink_r.x1 < W - 2 and ink_r.y1 < H - 2, f"all ink inside the page with a 2 pt margin: {[round(v, 2) for v in ink_r]} on {round(W, 2)} x {round(H, 2)}")
wds = [w_ for w_ in pg.get_text("words") if w_[4] in ("0-1", "1-5", "5-10", ">10")]
gaps = [round(b[0] - a[2], 2) for a in wds for b in wds if a is not b and abs(a[3] - b[3]) < 0.5 and 0 < b[0] - a[2] < 20]
check(gaps and min(gaps) >= 2.5, f"neighbouring x tick labels never touch: least gap {min(gaps) if gaps else None} pt over {len(gaps)} pairs (centred on their ticks in the wider boxes)")
clear_bad = []; least_clear = 1e9
for p, rec in zip(new["panels"], drawn["panels"]):
    top_line = p["axes"]["top"] - 4.0 - (len(rec["title_lines"]) - 1) * rec["title_pitch"] - 0.72 * rec["title_size"]
    above = [s for s in sp if s["bbox"][3] <= p["axes"]["top"] and s["bbox"][1] < top_line and s["text"] not in rec["title_lines"] and s["bbox"][0] < p["axes"]["right"] and s["bbox"][2] > p["axes"]["left"]]
    if above: least_clear = min(least_clear, top_line - max(s["bbox"][3] for s in above))
    if any(s["bbox"][3] > top_line - 2.0 for s in above): clear_bad.append(rec["title"])
check(not clear_bad, f"every title clears the text above it (the previous row's x tick labels) by 2 pt or more: least {least_clear:.2f} pt with the {R51['rowgap']:.0f} pt row gutter{': ' + str(clear_bad) if clear_bad else ''}")
bad_pos = 0
for s in sp:
    if not s["text"].strip(): continue
    if any(t["text"] == s["text"] and abs(t["x"] - s["origin"][0]) < 0.35 and at(t, s) for t in dr_t): continue
    if all(any(t["text"] == w_ and at(t, s) for t in dr_t) for w_ in s["text"].split()): continue
    bad_pos += 1
check(bad_pos == 0, f"every text span on the sheet is at the position the builder logged ({bad_pos} not found)")
d.close()

# 10. printed values csv
recs = []
for s in sp:
    if not s["text"].strip(): continue
    ents = [t for t in dr_t if abs(t["baseline"] - s["origin"][1]) < 0.35 and s["bbox"][0] - 0.35 <= t["x"] <= s["bbox"][2] + 0.35]
    ents = sorted(ents, key=lambda t: t["x"]); e = ents[0] if ents else None
    if e is None and s["text"] == "Hazard ratio vs 0-1% (95% CI)": e = [t for t in dr_t if t["what"] == "ylabel"][0]
    if e is None: continue
    what = e["what"]; base = dict(text=s["text"], x=s["bbox"][0], baseline=s["origin"][1], ha="left", size=s["size"], panel=e.get("panel") or "")
    if what == "title":
        full = e["key"]
        if s["text"] == full: recs.append(dict(base, source_file=f"{NUM}/bdsp_diseases_v3.csv", source_key=f"disease={full}:disease", source_value=full, rule="text"))
        else: recs.append(dict(base, source_file=DATA, source_key=f"graded.outcomes key '{full}', {e['value']} of the wrapped title", source_value=s["text"], rule="text"))
    elif what == "ytick": recs.append(dict(base, source_file=DATA, source_key=e["key"], source_value=s["text"], rule="text"))
    elif what == "star": recs.append(dict(base, source_file=DATA, source_key=e["key"] + f" -> BH q {e['value']:.3g} < 0.001", source_value="***", rule="text"))
    elif what == "xtick": recs.append(dict(base, source_file=DATA, source_key="band labels (graded.band_n keys) " + s["text"], source_value=s["text"], rule="text"))
    else: recs.append(dict(base, source_file="static:V13", source_key=what, source_value=s["text"], rule="text"))
n_rows, n_no = printed_values_csv(NEW, recs, f"{VER}/ED_Fig03_printed_values.csv", static_from=BASE)
check(n_no == 0, f"printed values: {n_rows} strings on the sheet accounted for (drawn record with its source), {n_no} unaccounted (NO)")
ok = ck.write(f"{VER}/checks_03_verify.txt"); print("RESULT ALL PASS" if ok else "RESULT FAIL")
