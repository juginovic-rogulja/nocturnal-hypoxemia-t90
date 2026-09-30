#!/usr/bin/env python3
"""ED_Fig03, V14 lane L2 (round 37), step 03: Proof W of ED_Fig03.pdf against the V13 design baseline and the v8.1 numbers.
Repointed and extended copy of the round-30 verifier (byte copy beside as 03_verify_ed3_PRE_V8_1.py): the declared deltas are the row
set (the August rule on the v8.1 graded block, GAP-11 closed: entering / leaving rows and the five v8.1 controls), the titles that
follow it (wrapped by the generator's rule), the y tick labels (the tick rule on the new limits) and the asterisks (BH q across the
file's top-band contrasts). Writes verify/checks.txt, readback_tokens.csv, expected_delta.json, values_agreement.txt,
CHANGES_ED_Fig03.csv, ED_Fig03_printed_values.csv, crops/*.png (200 dpi OLD = V13 vs NEW per panel block)."""
import fitz, re, os, sys, gc, csv, collections, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ed3_common import *
from ed3_extract import extract, spans_of
from v14lib import printed_values_csv, Checks, write_changes

NEW = f"{LANE}/ED_Fig03.pdf"; DRAWN = f"{VER}/ED_Fig03_drawn.json"
INK = (0x1a / 255, 0x1d / 255, 0x21 / 255); INK_INT = 0x1a1d21
LADDER = [(0xb3 / 255, 0xdc / 255, 0xf2 / 255), (0x7c / 255, 0xc0 / 255, 0xe9 / 255), (0x3f / 255, 0x9f / 255, 0xd8 / 255), (0x02 / 255, 0x88 / 255, 0xd1 / 255)]
HEX = ["#b3dcf2", "#7cc0e9", "#3f9fd8", "#0288d1"]; BAND_LABELS = {"0-1": "0-1%", "1-5": "1-5%", "5-10": "5-10%", ">10": ">10%"}
TOL = 0.3
for p in (BASE, NEW, DATA, PANEL_CSV, OLD_DATA, OLD_PANEL, DRAWN): hydrated(p)
drawn = load_json(DRAWN); meta = drawn["meta"]
ck = Checks(f"ED_Fig03 V14 lane L2 (round 37), Proof W: {NEW} (sha256 {sha256(NEW)[:16]}) against the V13 sheet {BASE} (sha256 {sha256(BASE)[:16]}) and "
            f"numbers/results_v2.json graded (step 112) + negcontrols_v4_panel.csv (step 251); OLD side = the v7 snapshot {OLD_DATA}. "
            f"Declared: row set by the August rule (entering {meta['entering']}, leaving {meta['leaving']}; controls {meta['controls']}), titles, y tick labels, asterisks.")
def check(ok, msg): return ck.log(msg, bool(ok))
def info(msg): ck.lines.append("INFO: " + msg); print("INFO: " + msg[:300], flush=True)
def near(c, t, tol=0.003): return c is not None and max(abs(a - b) for a, b in zip(c, t)) < tol

cur = extract(BASE); new = extract(NEW)
out = dict(graded(DATA)["outcomes"]); outo = dict(graded(OLD_DATA)["outcomes"]); Gn, Go = graded(DATA), graded(OLD_DATA)
check(sha256(DATA) == drawn["sources"][DATA]["sha256"] and sha256(PANEL_CSV) == drawn["sources"][PANEL_CSV]["sha256"], "the numbers files read now have the sha256 the builder recorded (v8.1 sidecars checked by the builder)")
check(drawn["sources"][DATA]["v8_inputs"] and drawn["sources"][PANEL_CSV]["v8_inputs"], f"sidecars name data_frozen_v8_2026-09: results_v2.json step {drawn['sources'][DATA]['step']} ({drawn['sources'][DATA]['output_mtime']}), negcontrols_v4_panel.csv step {drawn['sources'][PANEL_CSV]['step']} ({drawn['sources'][PANEL_CSV]['output_mtime']})")
check(drawn["base_sheet_sha256"] == sha256(BASE), "the builder drew against the V13 sheet on disk now")

# 1. page, fonts
check(all(abs(a - b) < 0.01 for a, b in zip(new["mediabox"], cur["mediabox"])), f"page box equal within 0.01 pt: {[round(v, 3) for v in new['mediabox']]} (V13 {[round(v, 3) for v in cur['mediabox']]})")
d = fitz.open(NEW); pg = d[0]
check(len(d) == 1 and len(pg.get_xobjects()) == 0, f"single flat page: {len(d)} page, {len(pg.get_xobjects())} XObjects")
fitz.TOOLS.mupdf_warnings(reset=True); pg.get_pixmap(dpi=50); w = fitz.TOOLS.mupdf_warnings(reset=True)
check(not w, "MuPDF renders without warnings" + (": " + w if w else ""))
fn = pg.get_fonts(full=True); names_new = sorted(f[3].split("+")[-1] for f in fn); d.close()
d = fitz.open(BASE); fb = d[0].get_fonts(full=True); names_base = sorted(f[3].split("+")[-1] for f in fb); d.close()
check(names_new == names_base and all("Arial" in n for n in names_new), f"fonts: {names_new} (V13 {names_base})")
check(all(f[2] == "TrueType" and f[5] == "WinAnsiEncoding" for f in fn), "fonts embedded as simple TrueType WinAnsi, as the V13 sheet's")
sizes_new = sorted({s["size"] for s in new["spans"]}); sizes_base = sorted({s["size"] for s in cur["spans"]})
check(sizes_new == sizes_base, f"type sizes {sizes_new} (V13 {sizes_base})")

# 2. letters
def letters(sp): return {s["text"]: s for s in sp if s["text"] in ("a", "b") and s["size"] == 13.0}
ln, lc = letters(new["spans"]), letters(cur["spans"])
for L in "ab":
    ok = L in ln and L in lc and "Bold" in ln[L]["font"] and abs(ln[L]["origin"][0] - lc[L]["origin"][0]) < TOL and abs(ln[L]["origin"][1] - lc[L]["origin"][1]) < TOL and abs(ln[L]["bbox"][1] - lc[L]["bbox"][1]) < TOL
    check(ok, f"letter {L}: 13 pt {ln.get(L, {}).get('font')} at origin {ln.get(L, {}).get('origin')} (V13 {lc.get(L, {}).get('origin')})")

# 3. row set, geometry pins
rows_new = [p["title"] for p in new["panels"]]; rows_v13 = [p["title"] for p in cur["panels"]]
rule_a, rule_b = august_rows(out, controls_of(PANEL_CSV))
check(len(new["panels"]) == 20 and len(cur["panels"]) == 20, f"20 panels ({len(new['panels'])} new, {len(cur['panels'])} V13)")
check(rows_new == rule_a + rule_b == meta["rows_drawn"], "row set and order = the August rule on the v8.1 graded block (GAP-11 closed, declared): " + " | ".join(rows_new))
check(rows_v13 == meta["rows_v13"] and set(meta["entering"]) == set(rows_new) - set(rows_v13) and set(meta["leaving"]) == set(rows_v13) - set(rows_new), f"declared row delta: entering {meta['entering']}, leaving {meta['leaving']} (V13 rows recorded by the probe)")
dev_box = max(abs(pn["axes"][k] - pc["axes"][k]) for pn, pc in zip(new["panels"], cur["panels"]) for k in ("left", "right", "top", "bottom"))
check(dev_box < TOL, f"axes boxes (spines) within {TOL} pt of V13: max deviation {dev_box:.3f} pt")
dev_xt = max(abs(a - b) for pn, pc in zip(new["panels"], cur["panels"]) for a, b in zip(pn["x_ticks"], pc["x_ticks"])) if all(len(pn["x_ticks"]) == 4 == len(pc["x_ticks"]) for pn, pc in zip(new["panels"], cur["panels"])) else 99
check(dev_xt < TOL, f"x tick marks (4 per panel) within {TOL} pt: max deviation {dev_xt:.3f} pt")
title_words = set(w for t in rows_new + rows_v13 for w in t.split())
def is_title(s): return s["size"] == 9.5 and "Bold" not in s["font"] and bool(s["text"].split()) and s["text"].split()[0] in title_words and not re.fullmatch(r"[\d.]+", s["text"]) and s["text"] not in ("0-1", "1-5 5-10 >10", "1-5", "5-10", ">10")
def fixed_spans(sp):
    """text that does not depend on the data: everything but titles, y tick labels and asterisks"""
    return [s for s in sp if not re.fullmatch(r"[\d.]+", s["text"]) and "*" not in s["text"] and not is_title(s)]
fx_n, fx_c = fixed_spans(new["spans"]), fixed_spans(cur["spans"]); unmatched = []; maxdev = 0.0
for s in fx_n:
    cands = [c for c in fx_c if c["text"] == s["text"] and c["size"] == s["size"] and c["font"] == s["font"]]
    dm = min((max(abs(c["origin"][0] - s["origin"][0]), abs(c["origin"][1] - s["origin"][1])) for c in cands), default=99)
    if dm >= TOL: unmatched.append((s["text"], s["origin"], round(dm, 3)))
    maxdev = max(maxdev, dm if dm < 99 else 0)
check(not unmatched and len(fx_n) == len(fx_c), f"every x tick label, axis title, group title and key string sits within {TOL} pt of V13's (max {maxdev:.3f} pt, {len(fx_n)} vs {len(fx_c)} spans){': ' + str(unmatched[:4]) if unmatched else ''}")
# titles: each panel's title lines at the builder's positions, wrapped by the rule
title_bad = []
for p, rec in zip(new["panels"], drawn["panels"]):
    L, T = p["axes"]["left"], p["axes"]["top"]; lines = rec["title_lines"]
    for j, line in enumerate(lines):
        base_y = T - 4.0 - (len(lines) - 1 - j) * 10.0
        sp = [s for s in new["spans"] if s["text"] == line and s["size"] == 9.5 and abs(s["origin"][0] - L) < TOL and abs(s["origin"][1] - base_y) < TOL]
        if not sp: title_bad.append((rec["title"], line))
check(not title_bad and all(" ".join(r["title_lines"]) == r["title"] for r in drawn["panels"]), f"titles: every line of every title at its panel's left edge, 4 pt above the box, 10 pt line pitch, wrapped by the generator's rule; wrapped titles {[r['title'] for r in drawn['panels'] if len(r['title_lines']) > 1]}{': ' + str(title_bad[:3]) if title_bad else ''}")
ytl_bad = []
for pn in new["panels"]:
    L = pn["axes"]["left"]
    for y, lab in pn["y_ticks"]:
        if lab is None: ytl_bad.append((pn["title"], y, "no label")); continue
        sp = [s for s in new["spans"] if s["text"] == lab and s["size"] == 9.5 and abs(s["bbox"][2] - (L - 7.32 / 1.2)) < TOL and abs(s["origin"][1] - 0.377 * 9.5 - y) < TOL]
        if not sp: ytl_bad.append((pn["title"], y, lab))
check(not ytl_bad, f"y tick labels right-aligned 6.1 pt left of the spine and centred on their tick marks within {TOL} pt (data-driven set, declared): {ytl_bad[:3]}")

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
check(rib_ok and line_ok and mark_ok and ref_ok and star_ok, "ribbon ink 16%, white 1.5 pt line through 4 points, 4 ink markers r 2.0 with 0.5 pt white ring, dashed ink reference at 1, ink asterisks")

# 5. read-back of plotted values (through each panel's own y ticks) against the v8.1 file; the V13 sheet against the v7 snapshot
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
with open(f"{VER}/values_agreement.txt", "w") as f:
    f.write(f"ED_Fig03 V14: plotted values read back from the vectors through each panel's y ticks\nnew: {NEW}\nV13: {BASE}\nv8.1: {DATA} graded.outcomes\nold: {OLD_DATA}\n\n")
    for r in rows: f.write(f"{r[0]:42s} file {' '.join(f'{v:.3f}' for v in r[1]):28s} read {' '.join(f'{v:.3f}' for v in r[2]):28s} dev {r[3]:6.3f}% CI dev {r[4]:6.3f}%  star {r[5] or '-'}\n")

# 6. stars and the key sentence
q_all = bh_top(out); q_allo = bh_top(outo); q_nc = bh_top(out, [k for k in out if not out[k]["negative_control"]]); q15 = bh_top(out, rows_new[:15]); q5 = bh_top(out, rows_new[15:])
star_rule = ["***" if q_all[t] < 0.001 else "" for t in rows_new]
check([p["star"] for p in new["panels"]] == star_rule, f"asterisks = the rule (BH q of the top-band contrast across the {len(out)} top-band contrasts of results_v2.json graded, *** below 0.001): " + "".join("*" if s else "-" for s in star_rule))
agree = all(((q_all[k] < 0.001) == (q_nc.get(k, 1) < 0.001) == (q15.get(k, 1) < 0.001)) for k in rows_new[:15]) and all((q_all[k] < 0.001) == (q5[k] < 0.001) for k in rows_new[15:])
info(f"alternate families (non-controls, the 15 drawn, the 5 controls) give the same asterisks: {agree}; q max over the 15: all-family {max(q_all[k] for k in rows_new[:15]):.3g}, non-control family {max(q_nc[k] for k in rows_new[:15]):.3g}; q min over controls: all-family {min(q_all[k] for k in rows_new[15:]):.3g}, 5-family {min(q5[k] for k in rows_new[15:]):.3g}")
key2 = [s for s in new["spans"] if s["text"] == "top band, q below 0.001"]; all15 = all(q_all[t] < 0.001 for t in rows_new[:15])
check(len(key2) == 1 and all15, f"key sentence '*** top band, q below 0.001' (the V13 wording) printed once and true: every panel-a outcome has q < 0.001 (largest {max(q_all[t] for t in rows_new[:15]):.3g})")
check([p["star"] for p in new["panels"]] == [p["star"] for p in cur["panels"]], "asterisk pattern per slot identical to V13 (15 ***, 5 none)")

# 7. numeric tokens: every token matched to a source; the expected delta
def tokens(sp):
    out_ = []
    for s in sp:
        for t in s["text"].split():
            if re.search(r"\d", t): out_.append((t, s))
    return out_
rb = []; unmatched_tok = []
for t, s in tokens(new["spans"]):
    x, y = s["origin"]; srcf = key = val = None
    if is_title(s): srcf, key, val = DATA, f"condition name graded.outcomes['{[r for r in rows_new if s['text'] in r][0] if any(s['text'] in r for r in rows_new) else s['text']}']", t
    elif re.fullmatch(r"[\d.]+", t) and s["size"] == 9.5:
        pn = [p for p in new["panels"] if abs(s["bbox"][2] - (p["axes"]["left"] - 7.32 / 1.2)) < TOL and p["axes"]["top"] - 6 <= y <= p["axes"]["bottom"] + 6]
        if len(pn) == 1:
            ymin, ymax, ticks = ylim_ticks(out[pn[0]["title"]])
            if t in ticks: srcf, key, val = DATA, f"tick rule on ({ymin:.4f}, {ymax:.4f}) from graded.outcomes['{pn[0]['title']}'] lo/hi", float(t)
    elif t in BAND_LABELS and s["size"] in (9.0, 9.5): srcf, key, val = DATA, f"graded.band key '{BAND_LABELS[t]}'", BAND_LABELS[t]
    elif s["text"] in ("1% or less", "over 1 to 5%", "over 5 to 10%", "over 10%"): srcf, key, val = DATA, f"band wording for {BAND_KEYS[['1% or less', 'over 1 to 5%', 'over 5 to 10%', 'over 10%'].index(s['text'])]}", s["text"]
    elif s["text"] == "top band, q below 0.001" and t == "0.001": srcf, key, val = DATA, "asterisk threshold, BH q of the top-band contrast", 0.001
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
old_ticks = collections.Counter(); new_ticks = collections.Counter(); per_panel = []
for pc in cur["panels"]: old_ticks.update(ylim_ticks(outo[pc["title"]])[2])
for p in new["panels"]: new_ticks.update(ylim_ticks(out[p["title"]])[2]); per_panel.append(dict(title=p["title"], new_ticks=ylim_ticks(out[p["title"]])[2]))
title_tok_new = collections.Counter(t for r in rows_new for t in r.split() if re.search(r"\d", t)); title_tok_old = collections.Counter(t for r in rows_v13 for t in r.split() if re.search(r"\d", t))
exp_add, exp_rem = (new_ticks + title_tok_new) - (old_ticks + title_tok_old), (old_ticks + title_tok_old) - (new_ticks + title_tok_new)
got_add, got_rem = tn - tc, tc - tn
check(got_add == exp_add and got_rem == exp_rem, f"expected delta: numeric tokens added {dict(got_add)} removed {dict(got_rem)} = tick rule on NEW rows minus OLD rows plus the title digits (added {dict(exp_add)}, removed {dict(exp_rem)})")
json.dump(dict(rule="y tick labels and title digits are the only data-dependent numeric tokens; expected delta = tick rule(new rows, v8.1) - tick rule(V13 rows, v7 snapshot) + title digits", old_file=OLD_DATA, new_file=DATA,
               expected_added=dict(exp_add), expected_removed=dict(exp_rem), sheet_added=dict(got_add), sheet_removed=dict(got_rem), per_panel=per_panel, base_tokens=dict(tc), new_tokens=dict(tn)), open(f"{VER}/expected_delta.json", "w"), indent=1)

# 8. word multiset: identical apart from the titles and the y tick labels (declared)
def collapse_(sp):
    keep = []
    for s in sp:
        if any(k["text"] == s["text"] and k["size"] == s["size"] and abs(k["origin"][0] - s["origin"][0]) < 0.25 and abs(k["origin"][1] - s["origin"][1]) < 0.25 for k in keep): continue
        keep.append(s)
    return keep
def words(sp):
    c = collections.Counter()
    for s in collapse_(sp): c.update(s["text"].split())
    return c
wn, wc = words(new["spans"]), words(cur["spans"])
title_w_new = collections.Counter(w for r in rows_new for w in r.split()); title_w_old = collections.Counter(w for r in rows_v13 for w in r.split())
wn_f = collections.Counter({k: v for k, v in wn.items() if not re.fullmatch(r"[\d.]+", k)}) - title_w_new; wc_f = collections.Counter({k: v for k, v in wc.items() if not re.fullmatch(r"[\d.]+", k)}) - title_w_old
check(wn_f == wc_f, f"word multiset identical to V13 apart from the titles and the y tick labels (declared): {sum(wn.values())} words; title words added {dict(title_w_new - title_w_old)}, removed {dict(title_w_old - title_w_new)}")

# 9. overlaps, bounds, crowding, text layer vs drawn
d = fitz.open(NEW); pg = d[0]; W, H = pg.rect.width, pg.rect.height; sp = new["spans"]; boxes = [p["axes"] for p in new["panels"]]
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
gaps = [round(b[0] - a[2], 2) for a in wds for b in wds if a is not b and abs(a[3] - b[3]) < 0.5 and 0 < b[0] - a[2] < 12]
check(gaps and min(gaps) >= 2.5, f"neighbouring x tick labels never touch: least gap {min(gaps) if gaps else None} pt over {len(gaps)} pairs")
# the three-line titles must clear the row above (its x tick labels) by at least 2 pt
clear_bad = []
for p, rec in zip(new["panels"], drawn["panels"]):
    if len(rec["title_lines"]) < 2: continue
    top_line = p["axes"]["top"] - 4.0 - (len(rec["title_lines"]) - 1) * 10.0 - 0.72 * 9.5
    above = [s for s in sp if s["bbox"][3] <= p["axes"]["top"] and s["bbox"][1] < top_line and s["text"] not in rec["title_lines"] and s["bbox"][0] < p["axes"]["right"] and s["bbox"][2] > p["axes"]["left"]]
    if any(s["bbox"][3] > top_line - 2.0 for s in above): clear_bad.append(rec["title"])
check(not clear_bad, f"wrapped titles clear the text above them by 2 pt or more{': ' + str(clear_bad) if clear_bad else ''}")
dr_t = drawn["text"]; bad_pos = 0
def at(t, s): return abs(t["baseline"] - s["origin"][1]) < 0.35 and s["bbox"][0] - 0.35 <= t["x"] <= s["bbox"][2]
for s in sp:
    if not s["text"].strip(): continue
    if any(t["text"] == s["text"] and abs(t["x"] - s["origin"][0]) < 0.35 and at(t, s) for t in dr_t): continue
    if all(any(t["text"] == w_ and at(t, s) for t in dr_t) for w_ in s["text"].split()): continue
    bad_pos += 1
check(bad_pos == 0, f"every text span on the sheet is at the position the builder logged (merged x tick label spans matched word by word; {bad_pos} not found)")
d.close()

# 10. crops OLD (V13) vs NEW at 200 dpi per panel block
split = (new["panels"][14]["axes"]["bottom"] + 13.25 + new["panels"][15]["axes"]["top"] - 24.0 - 20.0) / 2.0
clips = dict(panel_a=fitz.Rect(0, 0, W, split), panel_b=fitz.Rect(0, split, W, H))
for tag, path in (("old", BASE), ("new", NEW)):
    dd = fitz.open(path); pg = dd[0]
    for pn, clip in clips.items(): pg.get_pixmap(dpi=200, clip=clip).save(f"{CROPS}/{pn}_{tag}_200dpi.png")
    dd.close(); gc.collect()
for pn in clips:
    a = fitz.Pixmap(f"{CROPS}/{pn}_old_200dpi.png"); b = fitz.Pixmap(f"{CROPS}/{pn}_new_200dpi.png")
    A = np.frombuffer(a.samples, dtype=np.uint8).reshape(a.height, a.width, a.n); B = np.frombuffer(b.samples, dtype=np.uint8).reshape(b.height, b.width, b.n)
    h = max(A.shape[0], B.shape[0])
    def pad(x): o = np.full((h, x.shape[1], x.shape[2]), 255, np.uint8); o[:x.shape[0]] = x; return o
    S = np.hstack([pad(A), np.full((h, 16, A.shape[2]), 200, np.uint8), pad(B)])
    fitz.Pixmap(fitz.csRGB if S.shape[2] == 3 else fitz.csGRAY, S.shape[1], S.shape[0], S.tobytes(), False).save(f"{CROPS}/{pn}_OLD_vs_NEW_200dpi.png")
    del a, b, A, B, S; gc.collect()
check(os.path.exists(f"{LANE}/ED_Fig03_150dpi.png") and all(os.path.exists(f"{CROPS}/{pn}_OLD_vs_NEW_200dpi.png") for pn in clips), "renders: ED_Fig03_150dpi.png and verify/crops/panel_{a,b}_OLD_vs_NEW_200dpi.png")

# 11. CHANGES csv (V14 columns): rows leaving, entering, per-band values, stars, band n
def crosses(lo, hi): return lo <= 1.0 <= hi
ch = []
for i, t in enumerate(rows_v13):
    if t not in rows_new:
        vo = outo[t]; ch.append(dict(sheet="ED_Fig03", panel="a" if i < 15 else "b", item=f"{t} (row leaves)", before=f"top band {vo['>10%']['hr']:.3f} ({vo['>10%']['lo']:.3f}-{vo['>10%']['hi']:.3f}), star {cur['panels'][i]['star'] or 'none'}", after="", dramatic="yes", note="control swap of 14 September" if i >= 15 else "August rule on the v8.1 graded block (GAP-11 closed)"))
for i, t in enumerate(rows_new):
    v = out[t]; panel = "a" if i < 15 else "b"
    if t not in rows_v13:
        ch.append(dict(sheet="ED_Fig03", panel=panel, item=f"{t} (row enters)", before="", after=f"top band {v['>10%']['hr']:.3f} ({v['>10%']['lo']:.3f}-{v['>10%']['hi']:.3f}), star {new['panels'][i]['star'] or 'none'}", dramatic="yes", note="new outcome or new control under v8.1" if t not in outo else "August rule on the v8.1 graded block")); continue
    vo = outo[t]; j = rows_v13.index(t)
    for b in BAND_KEYS[1:]:
        flip = crosses(vo[b]["lo"], vo[b]["hi"]) != crosses(v[b]["lo"], v[b]["hi"]); big = abs(v[b]["hr"] - vo[b]["hr"]) / vo[b]["hr"] > 0.20
        why = ([("95% CI now crosses 1" if crosses(v[b]["lo"], v[b]["hi"]) else "95% CI no longer crosses 1")] if flip else []) + (["hazard ratio moves more than 20 percent"] if big else [])
        if abs(v[b]["hr"] - vo[b]["hr"]) > 5e-4 or abs(v[b]["lo"] - vo[b]["lo"]) > 5e-4 or abs(v[b]["hi"] - vo[b]["hi"]) > 5e-4 or why:
            ch.append(dict(sheet="ED_Fig03", panel=panel, item=f"{t}, {b} band hazard ratio (95% CI), marker and ribbon", before=f"{vo[b]['hr']:.3f} ({vo[b]['lo']:.3f}-{vo[b]['hi']:.3f})", after=f"{v[b]['hr']:.3f} ({v[b]['lo']:.3f}-{v[b]['hi']:.3f})", dramatic="yes" if why else "no", note="; ".join(why) if why else ("slot moved" if i != j else "digits only")))
    qo, qn = q_allo[t], q_all[t]; why = (["q crosses 0.05"] if (qo < 0.05) != (qn < 0.05) else []) + (["q crosses 0.001, the asterisk threshold"] if (qo < 0.001) != (qn < 0.001) else [])
    if why or cur["panels"][j]["star"] != new["panels"][i]["star"]: ch.append(dict(sheet="ED_Fig03", panel=panel, item=f"{t}, top-band asterisk", before=f"{cur['panels'][j]['star'] or 'none'} (q {qo:.3g})", after=f"{new['panels'][i]['star'] or 'none'} (q {qn:.3g})", dramatic="yes" if why else "no", note="; ".join(why)))
    bt = sorted({l for _, l in cur["panels"][j]["y_ticks"] if l}, key=float); nt = sorted({l for _, l in new["panels"][i]["y_ticks"] if l}, key=float)
    if bt != nt or i != j: ch.append(dict(sheet="ED_Fig03", panel=panel, item=f"{t}, y tick labels (slot {j + 1} -> {i + 1})", before=" ".join(bt), after=" ".join(nt), dramatic="no", note="axis ticks follow the data-derived limits"))
for b in BAND_KEYS:
    mv = abs(Gn["band_n"][b] - Go["band_n"][b]) / Go["band_n"][b]
    ch.append(dict(sheet="ED_Fig03", panel="not printed", item=f"band n {b} (legend, not on this sheet)", before=Go["band_n"][b], after=Gn["band_n"][b], dramatic="yes" if mv > 0.20 else "no", note=f"count moves {mv * 100:.1f} percent"))
ch.append(dict(sheet="ED_Fig03", panel="key", item="key sentence", before="*** top band, q below 0.001", after="*** top band, q below 0.001", dramatic="no", note=f"re-derived: largest q over the 15 is {max(q_all[t] for t in rows_new[:15]):.3g}, smallest over the controls {min(q_all[t] for t in rows_new[15:]):.3g}; asterisk family {len(out)} top-band contrasts (was {len(outo)})"))
write_changes(f"{VER}/CHANGES_ED_Fig03.csv", ch)
info(f"CHANGES_ED_Fig03.csv: {len(ch)} rows, {sum(r['dramatic'] == 'yes' for r in ch)} dramatic")

# 12. printed values csv: one record per text span on the sheet, classified through the builder's log (merged x tick spans use their first word)
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
    elif what == "key" and (e.get("key") or "").startswith("max BH q"): recs.append(dict(base, source_file=DATA, source_key=e["key"] + f" = {e['value']:.3g} (< 0.001)", source_value=s["text"], rule="text"))
    else: recs.append(dict(base, source_file="static:V13", source_key=what, source_value=s["text"], rule="text"))
n_rows, n_no = printed_values_csv(NEW, recs, f"{VER}/ED_Fig03_printed_values.csv", static_from=BASE)
check(n_no == 0, f"printed values: {n_rows} strings on the sheet accounted for (drawn record with its source, or identical to V13), {n_no} unaccounted (NO)")
ok = ck.write(f"{VER}/checks.txt"); print("RESULT ALL PASS" if ok else "RESULT FAIL")
