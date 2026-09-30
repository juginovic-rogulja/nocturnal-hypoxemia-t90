#!/usr/bin/env python3
"""ED_Fig07 verify (V14 lane L4, round 37) against the V13 design baseline: page box, letters within 0.3 pt, fonts, every drawn string of
both panels read back at its origin, the static V13 strings (ticks, headers, keys, unchanged row labels) at their V13 origins, the
"HR (95% CI)" header at the V13 position, values re-derived from the numbers files, the axis and marker geometry read back from the
panel PDF (get_drawings on the small flat panel), word and numeric deltas declared, colours within the palette, 150 dpi render and
raster diff (whole sheet rebuilt: the diff is reported, the OLD vs NEW crops written at 200 dpi), the printed-values csv, CHANGES csv."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import collections, csv, gc, json, math, os, re, sys, time
import fitz, numpy as np
from PIL import Image
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ed7spec as S
L = S.L
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND37_2026-09-14/figures/_common"); import v14lib as V
from statsmodels.stats.multitest import multipletests

CL = json.load(open(f"{S.WORK}/compose_log.json")); assert CL["output_sha256"] == L.sha256(S.OUT_PDF), "the sheet on disk is not the composed one"
BT = json.load(open(L.hydrated(S.BASE_TEXT))); assert L.sha256(S.BASE_PDF) == BT["sha256"]
DA = json.load(open(S.DRAWN_A)); DB = json.load(open(S.DRAWN_B)); LAY = json.load(open(S.LAYOUT))
CK = L.Checks(f"ED_Fig07 V14 (lane V14_L4_APNEA, round 37, {time.strftime('%Y-%m-%d %H:%M')}): V13 base {S.BASE_PDF} sha {BT['sha256'][:16]}, output sha {CL['output_sha256'][:16]}. Both panels rebuilt from the v8.1 numbers in the V13 design.")
d = fitz.open(S.OUT_PDF); p = d[0]; W, H = p.rect.width, p.rect.height; spans = L.spans_of(p); words = L.words_of(p); fonts_out = sorted(set((f[3].split("+")[-1], f[2]) for f in p.get_fonts(full=True))); n_pages = d.page_count; d.close(); gc.collect()
SP = L.dedupe(spans); SI = L.dedupe(BT["spans"]); W_IN, H_IN = BT["page"]
def norm(t): return t.replace("\xa0", " ")
CK.check("one page, page box equal to V13 (484.72 x 900.63)", n_pages == 1 and abs(W - W_IN) < 0.01 and abs(H - H_IN) < 0.01, f"{W:.3f} x {H:.3f}")
lo = {t[0]: t for t in L.letters(spans)}; li = {t[0]: t for t in BT["letters"]}
CK.check("letters a and b at the V13 positions within 0.3 pt, Arial Bold 13", sorted(lo) == ["a", "b"] and all(abs(lo[k][1] - li[k][1]) <= 0.3 and abs(lo[k][2] - li[k][2]) <= 0.3 and lo[k][4] == "Arial-BoldMT" and lo[k][5] == 13.0 for k in "ab"), f"{[(k, lo[k][1], lo[k][2]) for k in sorted(lo)]}")
fi = set(tuple(x) for x in BT["fonts"]); fo = set(fonts_out)
CK.check("fonts within the V13 set (Arial faces)", (fo - fi) <= set(), f"gained {sorted(fo - fi)} lost {sorted(fi - fo)}")
def present(rec, tol=0.35):
    for s in SP:
        if norm(s["text"]) != norm(rec["text"]): continue
        if abs(s["origin"][1] - rec["baseline"]) > tol: continue
        ha = rec.get("ha", "left"); sx = s["origin"][0] if ha == "left" else (s["bbox"][2] if ha == "right" else (s["bbox"][0] + s["bbox"][2]) / 2)
        if abs(sx - rec["x"]) <= max(tol, 0.6): return True
    return False
RECS = {"a": DA["drawn_records"], "b": DB["drawn_records"]}
for pnl in "ab":
    miss = [(r["text"], r["x"], r["baseline"]) for r in RECS[pnl] if not present(r)]
    CK.check(f"panel {pnl}: every drawn string ({len(RECS[pnl])}) read back from the sheet's text layer at its origin (0.35 pt)", not miss, f"missing {miss[:6]}")
hdr = [s for s in SP if norm(s["text"]) == S.HDR_A and s["origin"][1] < 100]; hdr13 = [s for s in SI if norm(s["text"]) == S.HDR_A and s["origin"][1] < 100]
CK.check("'HR (95% CI)' column title of panel a at the V13 position (right edge and baseline within 0.3 pt), Arial 9.5", len(hdr) == 1 and len(hdr13) == 1 and abs(hdr[0]["bbox"][2] - hdr13[0]["bbox"][2]) <= 0.3 and abs(hdr[0]["origin"][1] - hdr13[0]["origin"][1]) <= 0.3 and abs(hdr[0]["size"] - 9.5) < 0.05, f"{[(h['bbox'][2], h['origin'][1]) for h in hdr]} vs V13 {[(h['bbox'][2], h['origin'][1]) for h in hdr13]}")
extra = [norm(s["text"]) for s in SP if not (len(s["text"]) == 1 and s["size"] == 13.0) and norm(s["text"]) != S.HDR_A and not any(norm(s["text"]) == norm(r["text"]) and abs(s["origin"][1] - r["baseline"]) < 0.35 for pnl in "ab" for r in RECS[pnl])]
CK.check("no string on the sheet beyond the drawn records, the letters and the column title", not extra, f"extra {extra[:8]}")
# static V13 strings (everything that is not a value or a row label) at their V13 origins
def sig(s): return (norm(s["text"]), round(s["origin"][0], 2), round(s["origin"][1], 2))
val_re = re.compile(r"^(\d+\.\d\d( \(\d+\.\d\d-\d+\.\d\d\))?\**|\(\d+\.\d\d-\d+\.\d\d\)|\*+)$")
labels13 = set(norm(s["text"]) for s in SI if s["font"] == "ArialMT" and abs(s["size"] - 10) < 0.01 and abs(s["bbox"][0] - S.LETTER_X) < 0.3 and s["origin"][1] < 320 and not s["text"].startswith("*"))
static13 = sorted(sig(s) for s in SI if not val_re.match(norm(s["text"])) and norm(s["text"]) not in labels13 and not (len(s["text"]) == 1 and s["size"] == 13.0))
static14 = sorted(sig(s) for s in SP if not val_re.match(norm(s["text"])) and norm(s["text"]) not in labels13 and not (len(s["text"]) == 1 and s["size"] == 13.0) and norm(s["text"]) not in set(DA["row_order_top_to_bottom"]))
CK.check("every static V13 string (ticks, headers, keys, column and group labels, the x title) at its V13 origin (0.01 pt)", static13 == static14, f"{len(static13)} V13, {len(static14)} V14; diff {L.delta([json.dumps(x, ensure_ascii=False) for x in static13], [json.dumps(x, ensure_ascii=False) for x in static14]) if static13 != static14 else 'none'}"[:500])
row13 = sorted(labels13); row14 = sorted(DA["row_order_top_to_bottom"])
CK.check("panel a row set = the v8.1 file's outcomes (7 conditions + the 5 v8.1 controls); rows entering/leaving DECLARED", set(row14) == set(DA["rows"][i]["label"] for i in range(len(DA["rows"]))) and len(row14) == 12, f"V13 {row13}; V14 {row14}; enter {sorted(set(row14) - set(row13))} leave {sorted(set(row13) - set(row14))}")
# values re-derived
SEV = json.load(open(L.hydrated(S.SRC_A))); V.sidecar_v8(S.SRC_A); X = json.load(open(L.hydrated(S.SRC_B))); V.sidecar_v8(S.SRC_B); C9 = X["cross9"]
bad_a = [r["label"] for r in DA["rows"] if r["estimate_text"] != L.hr_ci(SEV["outcomes"][r["label"]]["hr"], SEV["outcomes"][r["label"]]["lo"], SEV["outcomes"][r["label"]]["hi"]) or r["star_text"] != L.stars(SEV["outcomes"][r["label"]]["p"])]
CK.check("panel a: every estimate and star re-derived from severe_apnea_negctrl_v1.json (half up 2 dp, raw P stars)", not bad_a and len(DA["rows"]) == len(SEV["outcomes"]), f"{bad_a}")
CONDS9 = [c for c, v in C9.items() if not v["negative_control"]]; Q9 = {c: {} for c in CONDS9}
for cell in [f"{g}|{l}" for g in S.GROUPS for l in S.LEVELS]:
    if cell == S.REF: continue
    for c, q in zip(CONDS9, multipletests([C9[c]["cells"][cell]["p"] for c in CONDS9], method="fdr_bh")[1]): Q9[c][cell] = float(q)
bad_b = [(v["row"], v["column"]) for v in DB["values"] if v["row"] != S.REF and (v["hr_text"] != L.r2(C9[v["column"]]["cells"][v["row"]]["hr"]) + L.stars(Q9[v["column"]][v["row"]]) or v["ci_text"] != f"({L.r2(C9[v['column']]['cells'][v['row']]['lo'])}-{L.r2(C9[v['column']]['cells'][v['row']]['hi'])})" or v["fill"] != L.ramp_hex(C9[v["column"]]["cells"][v["row"]]["hr"]))]
CK.check("panel b: all 36 cells re-derived from alenfig3_cross_v1.json cross9 (half up 2 dp, BH stars within the cell, ramp fill)", not bad_b and len(DB["values"]) == 36, f"{bad_b[:4]}")
CK.check("panel b: cell sizes sum to the file's cohort n (18,271 in v8.1) = crosstab n", sum(DB["cell_sizes"].values()) == DB["cohort_n"], f"{sum(DB['cell_sizes'].values())}")
# panel a geometry read back from the flat panel PDF (small: get_drawings allowed)
GA = DA["geometry"]; dd = fitz.open(S.PANEL_A); DR = dd[0].get_drawings(); dd.close()
def segs():
    for dr in DR:
        for it in dr["items"]:
            if it[0] == "l": yield dr, it[1], it[2]
hz = [(dr, a, b) for dr, a, b in segs() if abs(a.y - b.y) < 1e-3]; spine = max(hz, key=lambda t: abs(t[2].x - t[1].x)); A_, B_ = GA["log_map"]["A"], GA["log_map"]["B"]
ci = [(dr, a, b) for dr, a, b in hz if abs(dr["width"] - GA["ci_lw"]) < 0.05]; rows_y = S.cluster([a.y for dr, a, b in ci], 0.05)
ok_ci = len(ci) == len(DA["rows"]) and all(any(abs(min(a.x, b.x) - r["ci_x"][0]) < 0.05 and abs(max(a.x, b.x) - r["ci_x"][1]) < 0.05 and abs(a.y - r["row_y_panel"]) < 0.05 for dr, a, b in ci) for r in DA["rows"])
CK.check("panel a: spine at the V13 axes (x 133.92 to 375.84, y 303.21 in panel coordinates), every interval read back at the fitted ln(lo), ln(hi) and its row y (0.05 pt)", abs(spine[1].y - GA["axes_panel"][3]) < 0.05 and abs(min(spine[1].x, spine[2].x) - GA["axes_panel"][0]) < 0.05 and ok_ci, f"spine y {spine[1].y:.3f}, {len(ci)} intervals, rows {len(rows_y)}")
circ = [dr for dr in DR if dr.get("fill") not in (None, (1.0, 1.0, 1.0)) and all(it[0] == "c" for it in dr["items"]) and len(dr["items"]) >= 4]
sq = [dr for dr in DR if dr.get("fill") not in (None, (1.0, 1.0, 1.0)) and len(dr["items"]) == 1 and dr["items"][0][0] == "re" and abs((dr["rect"].x1 - dr["rect"].x0) - (dr["rect"].y1 - dr["rect"].y0)) < 0.05 and (dr["rect"].x1 - dr["rect"].x0) > 3]
CK.check("panel a: markers = circles for conditions (#0288d1) and squares for controls (#8a9099), counts as the row set", len(circ) == DA["n_conditions"] and len(sq) == DA["n_controls"] and all(L.CM.rgb_to_hex(c["fill"]) == GA["colours"]["condition"] for c in circ) and all(L.CM.rgb_to_hex(c["fill"]) == GA["colours"]["control"] for c in sq), f"{len(circ)} circles, {len(sq)} squares")
mk_x = sorted(round((dr["rect"].x0 + dr["rect"].x1) / 2, 2) for dr in circ + sq); exp_x = sorted(round(r["marker_x"], 2) for r in DA["rows"])
CK.check("panel a: every marker centre at the fitted ln(HR) on the V13 log map (0.1 pt)", len(mk_x) == len(exp_x) and all(abs(a - b) < 0.1 for a, b in zip(mk_x, exp_x)), f"{list(zip(mk_x, exp_x))[:3]}")
# colours of both panels within the palette (+ the ramp)
RAMP_SWATCHES = {L.ramp_hex(math.sqrt(math.exp((k / 120) * math.log(3.25)) * math.exp(((k + 1) / 120) * math.log(3.25)))) for k in range(120)} | {v["fill"] for v in DB["values"]}   # the key's 120 swatches + every cell fill (asserted = ramp_hex(hr) by the builder)
def is_ramp(c): return c in RAMP_SWATCHES
bad_col = {}
for pnl, pdf in (("a", S.PANEL_A), ("b", S.PANEL_B)):
    dd = fitz.open(pdf); cols = collections.Counter()
    for dr in dd[0].get_drawings():
        for k in ("fill", "color"):
            if dr.get(k) is not None: cols[L.CM.rgb_to_hex(dr[k])] += 1
    dd.close(); off = {c: n for c, n in cols.items() if c not in V.PALETTE and not (pnl == "b" and is_ramp(c))}
    if off: bad_col[pnl] = off
CK.check("colours of both panels within the palette (blue, greys, ink, white, the orange and blue ladders, the red ramp)", not bad_col, f"{bad_col}")
# word and numeric deltas (whole sheet rebuilt): declared = the drawn records against the V13 strings
def toks(ss): return L.word_tokens([dict(text=norm(s["text"])) for s in ss]), L.num_tokens([dict(text=norm(s["text"])) for s in ss])
wi, ni = toks(SI); wo, no = toks(SP); wd = L.delta(wi, wo); nd = L.delta(ni, no)
json.dump(dict(words=wd, numbers=nd), open(f"{S.VER}/text_delta.json", "w"), indent=1, ensure_ascii=False)
new_all = L.num_tokens([dict(text=norm(r["text"])) for pnl in "ab" for r in RECS[pnl]] + [dict(text=S.HDR_A)]); exp = L.delta(ni, new_all)
CK.check("numeric-token delta on the sheet = the expected delta (drawn records against the V13 strings); word delta DECLARED in verify/text_delta.json", exp == nd, f"numbers lost {sum(nd['lost'].values())} gained {sum(nd['gained'].values())}; words lost {json.dumps(wd['lost'], ensure_ascii=False)[:200]} gained {json.dumps(wd['gained'], ensure_ascii=False)[:200]}")
# renders and crops
V.render(S.OUT_PDF, 150, S.OUT_PNG, heavy=False); OLD150 = f"{S.WORK}/OLD_V13_150dpi.png"; V.render(S.BASE_PDF, 150, OLD150, heavy=False)
rd = V.raster_diff(OLD150, S.OUT_PNG, [[0, 0, W, H]], 150)
CK.check("150 dpi renders of the same size; the rebuild is visible", rd["shape_equal"] and rd["n_diff"] > 100, f"{rd}")
for pnl, box in (("a", LAY["panel_a_rect"]), ("b", LAY["panel_b_rect"])):
    V.render(S.BASE_PDF, 200, f"{S.CROPS}/panel_{pnl}_OLD_V13_200dpi.png", heavy=False, clip=box); V.render(S.OUT_PDF, 200, f"{S.CROPS}/panel_{pnl}_NEW_200dpi.png", heavy=False, clip=box)
    a = Image.open(f"{S.CROPS}/panel_{pnl}_OLD_V13_200dpi.png"); b = Image.open(f"{S.CROPS}/panel_{pnl}_NEW_200dpi.png"); c = Image.new("RGB", (a.width + b.width + 20, max(a.height, b.height)), "white"); c.paste(a, (0, 0)); c.paste(b, (a.width + 20, 0)); c.save(f"{S.CROPS}/panel_{pnl}_OLD_vs_NEW_200dpi.png")
CK.check("200 dpi OLD vs NEW crops written for a and b", all(os.path.exists(f"{S.CROPS}/panel_{q}_OLD_vs_NEW_200dpi.png") for q in "ab"), "")
# printed values csv
V.RULES.update({"sev_hr_ci": lambda o: L.hr_ci(o["hr"], o["lo"], o["hi"]), "sev_stars_p": lambda o: L.stars(o["p"]),
                "cross9_hr_stars": lambda cell: next(L.r2(cell["hr"]) + L.stars(Q9[c][k]) for c in CONDS9 for k, v in C9[c]["cells"].items() if v == cell and not v["reference"]),
                "cell_ci_paren": lambda cell: f"({L.r2(cell['lo'])}-{L.r2(cell['hi'])})"})
recs = [dict(r, panel=pnl) for pnl in "ab" for r in RECS[pnl]] + [dict(text=S.HDR_A, x=hdr[0]["origin"][0] if hdr else 0, baseline=hdr[0]["origin"][1] if hdr else 0, ha="left", size=9.5, panel="a", source_file="static:V13", source_key="V13 column title (round 31)", source_value=S.HDR_A, rule="text")]
n_rows, n_no = V.printed_values_csv(S.OUT_PDF, recs, f"{S.VER}/{S.SHEET}_printed_values.csv", static_from=S.BASE_PDF)
rows_pv = list(csv.DictReader(open(f"{S.VER}/{S.SHEET}_printed_values.csv"))); cnt = collections.Counter(r["match"] for r in rows_pv)
CK.check(f"printed values csv: {n_rows} strings, {cnt.get('yes', 0)} re-derived from the source file, {cnt.get('yes_recorded', 0)} recorded, {cnt.get('static', 0)} static V13, {n_no} NO", n_no == 0, f"NO rows: {[(r['string'], r['note']) for r in rows_pv if r['match'] == 'NO'][:5]}")
# CHANGES csv
R30A = json.load(open(L.hydrated(f"{L.R30}/V7_L4_APNEA/ED_Fig07/verify/ED_Fig07_a_drawn.json"))); olda = {r["label"]: r for r in R30A["rows"]}; newa = {r["label"]: r for r in DA["rows"]}
v13_a_txt = set(norm(s["text"]) for s in SI); assert all(r["estimate_text"] in v13_a_txt for r in olda.values()), "the round-30 panel a record is not the V13 panel a"
rows = []
for k in sorted(set(olda) | set(newa)):
    o, n = olda.get(k), newa.get(k); rs = []
    if o is None: rs.append("row enters (v8.1 control panel)")
    elif n is None: rs.append("row leaves (control left the panel on 14 September)")
    else:
        if (o["star_text"] != "") != (n["star_text"] != ""): rs.append("P crosses 0.05")
        if (o["lo"] > 1 or o["hi"] < 1) != (n["lo"] > 1 or n["hi"] < 1): rs.append("CI crosses 1")
    rows.append(dict(sheet=S.SHEET, panel="a", item=f"{k} HR (95% CI), stars", before=(f"{o['estimate_text']} {o['star_text']}".strip() if o else "absent"), after=(f"{n['estimate_text']} {n['star_text']}".strip() if n else "absent"), dramatic="yes" if rs else "no", note="; ".join(rs)))
co, cn = R30A["cohort"], DA["cohort"]
for key in ("n", "n_reference", "n_low_oxygen"):
    mv = abs(cn[key] - co[key]) / co[key]; rows.append(dict(sheet=S.SHEET, panel="a", item=f"severe-apnea group {key} (not printed on the sheet, legend count)", before=f"{co[key]:,}", after=f"{cn[key]:,}", dramatic="yes" if mv > 0.2 else "no", note=f"moves {mv * 100:.1f}%"))
OLDB = json.load(open(L.hydrated(S.LF4_B_DRAWN))); oldb = {(v["row"], v["column"]): v for v in OLDB["values"]}
for v in DB["values"]:
    o = oldb.get((v["row"], v["column"])); rs = []
    if o and v["row"] != S.REF:
        if (o["stars"] != "") != (v["stars"] != ""): rs.append("q crosses 0.05")
        if (o["lo"] > 1 or o["hi"] < 1) != (v["lo"] > 1 or v["hi"] < 1): rs.append("CI crosses 1")
    rows.append(dict(sheet=S.SHEET, panel="b", item=f"{v['column']} | {v['row']}", before=(f"{o['hr_text']} {o['ci_text'] or ''}".strip() if o else "absent"), after=f"{v['hr_text']} {v['ci_text'] or ''}".strip(), dramatic="yes" if rs else "no", note="; ".join(rs)))
for cell, n_ in DB["cell_sizes"].items():
    o = OLDB["cell_sizes"].get(cell); mv = abs(n_ - o) / o if o else 1
    rows.append(dict(sheet=S.SHEET, panel="b", item=f"cell size {cell} (not printed on the sheet)", before=f"{o:,}" if o else "absent", after=f"{n_:,}", dramatic="yes" if mv > 0.2 else "no", note=f"moves {mv * 100:.1f}%"))
V.write_changes(S.CHANGES, rows)
CK.check(f"CHANGES csv written: {len(rows)} rows, {sum(1 for r in rows if r['dramatic'] == 'yes')} dramatic", os.path.exists(S.CHANGES), "")
ok = CK.write(S.CHECKS); print("RESULT ALL PASS" if ok else "RESULT FAIL")
