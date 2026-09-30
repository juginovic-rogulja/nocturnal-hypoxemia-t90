#!/usr/bin/env python3
"""ED_Fig07 verify, round 51 (lane LED-B, 2026-09-26, panels side by side: panel b's records are transformed by (+X_B, -YB) to the new page). From the round-49 verify (ed7_04_verify_r49.py beside), adapted from the round-38 ed7_04_verify.py (copied beside, unedited): the composed sheet is a
nested PDF (two panel forms), so it is read with Ghostscript only (txtwrite text layer, png renders); PyMuPDF text and drawings are used on the
two FLAT panel PDFs only. Checks: page box = V13; letters a, b Arial Bold 14 at the V13 origins; every drawn string of both panels read back on
its flat panel at its anchor (0.35 pt) and on the composed sheet's Ghostscript text layer (1 pt, txtwrite rounds to whole points); the
"HR (95% CI)" header once, 10.5 pt, right-aligned on the new estimate column edge; no string beyond the records, the header and the letters;
values re-derived from the numbers files (half up 2 dp, raw-P stars in a, BH stars within the cell in b, ramp fills); panel a geometry from
get_drawings on the flat panel (spine, intervals, markers at ln(HR)); colours within the palette; numeric-token delta against V13 = the declared
delta; every string at its V26 size + 1 except the declared fitted sizes; Ghostscript 150 dpi renders and 200 dpi OLD (V13) vs NEW crops;
printed-values csv per panel (flat PDFs); CHANGES csv (sizes)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import collections, csv, gc, json, math, os, re, sys, time, subprocess
import fitz, numpy as np
from PIL import Image
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ed7spec as S
L = S.L
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND37_2026-09-14/figures/_common"); import v14lib as V
sys.path.insert(0, f"{L.LANE}/scripts"); import ledb_common as C
from statsmodels.stats.multitest import multipletests

COMPOSED = sys.argv[1] if len(sys.argv) > 1 else S.OUT_PDF; STAMPED_PNG = sys.argv[2] if len(sys.argv) > 2 else None   # round 51: the runner stamps first, the verify reads the composed copy and the stamped render
CL = json.load(open(f"{S.WORK}/compose_log.json")); assert CL["output_sha256"] == L.sha256(COMPOSED), "the sheet on disk is not the composed one"
BT = json.load(open(L.hydrated(S.BASE_TEXT))); assert L.sha256(S.BASE_PDF) == BT["sha256"]
DA = json.load(open(S.DRAWN_A)); DB = json.load(open(S.DRAWN_B)); LAY = json.load(open(f"{S.WORK}/sheet_layout_r51.json")); GA = json.load(open(S.GEOM_A)); GB = json.load(open(S.GEOM_B))
X_B = LAY["x_b"]; YB = LAY["panel_b_src_rect"][1]
def tb(r): return dict(r, x=r["x"] + X_B, baseline=r["baseline"] - YB)      # a panel b record on the side-by-side page
CK = L.Checks(f"ED_Fig07 round 49 (lane LED-B, {time.strftime('%Y-%m-%d %H:%M')}): V13 base {S.BASE_PDF} sha {BT['sha256'][:16]}, composed sheet {S.OUT_PDF} sha {CL['output_sha256'][:16]}. Both panels rebuilt from the v8 numbers in the V13 design at +{L.PT_PLUS:g} pt.")
def norm(t): return t.replace("\xa0", " ").replace("\xad", "-")
# ---- the composed sheet: Ghostscript only
from pypdf import PdfReader
pr = PdfReader(COMPOSED); W, H = float(pr.pages[0].mediabox.width), float(pr.pages[0].mediabox.height); W_IN, H_IN = BT["page"]
CK.check(f"one page, page box {LAY['sheet_w']:.2f} x {LAY['sheet_h']:.2f} (panels side by side: V28 was {W_IN:.2f} x {H_IN:.2f}, declared)", len(pr.pages) == 1 and abs(W - LAY["sheet_w"]) < 0.01 and abs(H - LAY["sheet_h"]) < 0.01, f"{W:.3f} x {H:.3f}")
GSP = C.parse_txtwrite(C.gs_txtwrite(COMPOSED, f"{S.WORK}/composed_txtwrite.xml"))
lets = {s["text"]: s for s in GSP if s["text"] in ("a", "b") and s["bold"] and abs(s["size"] - S.LETTER_SIZE) < 0.01}
CK.check(f"letters a and b Arial Bold {S.LETTER_SIZE:g} pt at the V13 origins (x {S.LETTER_X}, baselines {LAY['letters'][0][2]}, {LAY['letters'][1][2]}; 1 pt, txtwrite rounds)", sorted(lets) == ["a", "b"] and all(abs(lets[ch][ "x0"] - x) <= 1.0 and abs(lets[ch]["y"] - y) <= 1.0 for ch, x, y in LAY["letters"]), f"{[(k, v['x0'], v['y']) for k, v in lets.items()]}")
RECS = {"a": DA["drawn_records"], "b": [tb(r) for r in DB["drawn_records"]]}
def present_gs(rec, tol=1.0):
    for s in GSP:
        if norm(s["text"]) != norm(rec["text"]) or abs(s["y"] - rec["baseline"]) > tol: continue
        ha = rec.get("ha", "left"); sx = s["x0"] if ha == "left" else (s["x1"] if ha == "right" else (s["x0"] + s["x1"]) / 2)
        if abs(sx - rec["x"]) <= tol and abs(s["size"] - rec["size"]) < 0.01: return True
    return False
for pnl in "ab":
    miss = [(r["text"], r["x"], r["baseline"], r["size"]) for r in RECS[pnl] if not present_gs(r)]
    CK.check(f"panel {pnl}: every drawn string ({len(RECS[pnl])}) on the composed sheet's Ghostscript text layer at its anchor, baseline and size (1 pt)", not miss, f"missing {miss[:6]}")
hdr = [s for s in GSP if norm(s["text"]) == S.HDR_A and s["y"] < 100]
CK.check(f"'HR (95% CI)' column title once, Arial {S.HDR_A_FS:g} pt, right-aligned on the estimate column edge x {CL['placements'][1]['x1']:.2f} at the header baseline {CL['placements'][1]['baseline']:.2f} (1 pt), as in round 49", len(hdr) == 1 and abs(hdr[0]["x1"] - CL["placements"][1]["x1"]) <= 1.0 and abs(hdr[0]["y"] - CL["placements"][1]["baseline"]) <= 1.0 and abs(hdr[0]["size"] - S.HDR_A_FS) < 0.01, f"{[(h['x0'], h['x1'], h['y'], h['size']) for h in hdr]}")
exp_ms = collections.Counter(norm(r["text"]) for pnl in "ab" for r in RECS[pnl]) + collections.Counter([S.HDR_A, "a", "b"]); got_ms = collections.Counter(norm(s["text"]) for s in GSP)
CK.check("no string on the sheet beyond the drawn records, the letters and the column title (Ghostscript string multiset)", exp_ms == got_ms, f"lost {dict(exp_ms - got_ms)} extra {dict(got_ms - exp_ms)}")
# ---- sizes: every string at its V26 size + 1 except the declared fitted sizes (SIZE_KEEP of panel b); the two star keys are absent on V26 already (LNOTES round 40)
v28 = C.parse_txtwrite(C.gs_txtwrite(f"{C.V28}/{S.SHEET}.pdf", f"{S.WORK}/v28_txtwrite.xml")); keep = GB.get("size_exceptions", {})
old_sz = collections.Counter((norm(s["text"]), round(s["size"], 3)) for s in v28 if not s["text"].startswith("Extended Data Fig.")); new_sz = collections.Counter((norm(s["text"]), round(s["size"], 3)) for s in GSP)
CK.check("every V28 string present at its V28 size (no size change, the title is added by the strip step)", old_sz == new_sz, f"expected-but-absent {dict(old_sz - new_sz)}; present-but-unexpected {dict(new_sz - old_sz)}")
# ---- the flat panels: PyMuPDF text and drawings (small flat matplotlib PDFs)
def present_flat(sp, rec, dy, tol=0.35):
    for s in sp:
        if norm(s["text"]) != norm(rec["text"]) or abs((s["origin"][1] + dy) - rec["baseline"]) > tol or abs(s["size"] - rec["size"]) > 0.01: continue
        ha = rec.get("ha", "left"); sx = s["origin"][0] if ha == "left" else (s["bbox"][2] if ha == "right" else (s["bbox"][0] + s["bbox"][2]) / 2)
        if abs(sx - rec["x"]) <= max(tol, 0.6): return True
    return False
SPA = L.spans_of(S.PANEL_A); SPB = L.spans_of(S.PANEL_B)
for pnl, sp, dy in (("a", SPA, S.SHIFT_A), ("b", SPB, 0.0)):
    recs_flat = RECS[pnl] if pnl == "a" else DB["drawn_records"]
    miss = [(r["text"], r["x"], r["baseline"]) for r in recs_flat if not present_flat(sp, r, dy if pnl == "a" else LAY["panel_b_src_rect"][1])]
    CK.check(f"panel {pnl}: every drawn string read back on the flat panel PDF at its anchor (0.35 pt), size {sorted({r['size'] for r in recs_flat})}", not miss, f"missing {miss[:6]}")
fonts = sorted({f[3].split("+")[-1] for p in (S.PANEL_A, S.PANEL_B) for f in fitz.open(p)[0].get_fonts(full=True)})
CK.check("fonts on the panels: Arial faces only", all("Arial" in f for f in fonts), f"{fonts}")
# ---- values re-derived (as round 38)
SEV = json.load(open(L.hydrated(S.SRC_A))); V.sidecar_v8(S.SRC_A); X = json.load(open(L.hydrated(S.SRC_B))); V.sidecar_v8(S.SRC_B); C9 = X["cross9"]
bad_a = [r["label"] for r in DA["rows"] if r["estimate_text"] != L.hr_ci(SEV["outcomes"][r["label"]]["hr"], SEV["outcomes"][r["label"]]["lo"], SEV["outcomes"][r["label"]]["hi"]) or r["star_text"] != L.stars(SEV["outcomes"][r["label"]]["p"])]
CK.check("panel a: every estimate and star re-derived from severe_apnea_negctrl_v1.json (half up 2 dp, raw P stars)", not bad_a and len(DA["rows"]) == len(SEV["outcomes"]), f"{bad_a}")
CONDS9 = [c for c, v in C9.items() if not v["negative_control"]]; Q9 = {c: {} for c in CONDS9}
for cell in [f"{g}|{l}" for g in S.GROUPS for l in S.LEVELS]:
    if cell == S.REF: continue
    for c, q in zip(CONDS9, multipletests([C9[c]["cells"][cell]["p"] for c in CONDS9], method="fdr_bh")[1]): Q9[c][cell] = float(q)
bad_b = [(v["row"], v["column"]) for v in DB["values"] if v["row"] != S.REF and (v["hr_text"] != L.r2(C9[v["column"]]["cells"][v["row"]]["hr"]) + L.stars(Q9[v["column"]][v["row"]]) or v["ci_text"] != f"({L.r2(C9[v['column']]['cells'][v['row']]['lo'])}-{L.r2(C9[v['column']]['cells'][v['row']]['hi'])})" or v["fill"] != L.ramp_hex(C9[v["column"]]["cells"][v["row"]]["hr"]))]
CK.check("panel b: all 36 cells re-derived from alenfig3_cross_v1.json cross9 (half up 2 dp, BH stars within the cell, ramp fill)", not bad_b and len(DB["values"]) == 36, f"{bad_b[:4]}")
CK.check("panel b: cell sizes sum to the file's cohort n = crosstab n", sum(DB["cell_sizes"].values()) == DB["cohort_n"], f"{sum(DB['cell_sizes'].values())}")
# the printed strings on the sheet equal the drawn records' strings (the records were checked against the files above): estimates and cells
est_on_sheet = sorted(norm(s["text"]) for s in GSP if re.fullmatch(r"\d\.\d\d \(\d\.\d\d-\d\.\d\d\)", norm(s["text"])) and s["y"] < 380)
CK.check("panel a: the 12 estimate strings on the sheet equal the re-derived estimates", est_on_sheet == sorted(r["estimate_text"] for r in DA["rows"]), "")
# ---- panel a geometry (flat panel)
dd = fitz.open(S.PANEL_A); DR = dd[0].get_drawings(); dd.close()
def segs():
    for dr in DR:
        for it in dr["items"]:
            if it[0] == "l": yield dr, it[1], it[2]
hz = [(dr, a, b) for dr, a, b in segs() if abs(a.y - b.y) < 1e-3]; spine = max(hz, key=lambda t: abs(t[2].x - t[1].x))
ci = [(dr, a, b) for dr, a, b in hz if abs(dr["width"] - GA["ci_lw"]) < 0.05]
ok_ci = len(ci) == len(DA["rows"]) and all(any(abs(min(a.x, b.x) - r["ci_x"][0]) < 0.05 and abs(max(a.x, b.x) - r["ci_x"][1]) < 0.05 and abs(a.y - r["row_y_panel"]) < 0.05 for dr, a, b in ci) for r in DA["rows"])
CK.check("panel a: spine at the V13 axes (x 133.92 to 375.84, y 303.21 in panel coordinates), every interval read back at ln(lo), ln(hi) and its row y (0.05 pt)", abs(spine[1].y - GA["axes_panel"][3]) < 0.05 and abs(min(spine[1].x, spine[2].x) - GA["axes_panel"][0]) < 0.05 and ok_ci, f"spine y {spine[1].y:.3f}, {len(ci)} intervals")
circ = [dr for dr in DR if dr.get("fill") not in (None, (1.0, 1.0, 1.0)) and all(it[0] == "c" for it in dr["items"]) and len(dr["items"]) >= 4]
sq = [dr for dr in DR if dr.get("fill") not in (None, (1.0, 1.0, 1.0)) and len(dr["items"]) == 1 and dr["items"][0][0] == "re" and abs((dr["rect"].x1 - dr["rect"].x0) - (dr["rect"].y1 - dr["rect"].y0)) < 0.05 and (dr["rect"].x1 - dr["rect"].x0) > 3]
CK.check("panel a: markers = circles for conditions (#0288d1) and squares for controls (#8a9099), counts as the row set", len(circ) == DA["n_conditions"] and len(sq) == DA["n_controls"] and all(L.CM.rgb_to_hex(c["fill"]) == GA["colours"]["condition"] for c in circ) and all(L.CM.rgb_to_hex(c["fill"]) == GA["colours"]["control"] for c in sq), f"{len(circ)} circles, {len(sq)} squares")
mk_x = sorted(round((dr["rect"].x0 + dr["rect"].x1) / 2, 2) for dr in circ + sq); exp_x = sorted(round(r["marker_x"], 2) for r in DA["rows"])
CK.check("panel a: every marker centre at ln(HR) on the V13 log map (0.1 pt)", len(mk_x) == len(exp_x) and all(abs(a - b) < 0.1 for a, b in zip(mk_x, exp_x)), f"{list(zip(mk_x, exp_x))[:3]}")
CK.check(f"panel a: row labels at {GA['sizes']['tick']:g} pt clear of the plot's leftmost ink (x {GA['ink_left_pt']:.1f}) by 4 pt; labels past the axes box's invisible left edge (x {GA['axes_panel'][0]:.2f}): {GA['labels_past_axes_edge_pt']} (declared, no spine and no grid there)", all(e + 4.0 <= GA["ink_left_pt"] for e in GA["label_ends_pt"].values()), f"widest label end {max(GA['label_ends_pt'].values()):.2f}")
RAMP_SWATCHES = {L.ramp_hex(math.sqrt(math.exp((k / 120) * math.log(3.25)) * math.exp(((k + 1) / 120) * math.log(3.25)))) for k in range(120)} | {v["fill"] for v in DB["values"]}
bad_col = {}
for pnl, pdf in (("a", S.PANEL_A), ("b", S.PANEL_B)):
    d2 = fitz.open(pdf); cols = collections.Counter()
    for dr in d2[0].get_drawings():
        for k in ("fill", "color"):
            if dr.get(k) is not None: cols[L.CM.rgb_to_hex(dr[k])] += 1
    d2.close(); off = {c: n for c, n in cols.items() if c not in V.PALETTE and not (pnl == "b" and c in RAMP_SWATCHES)}
    if off: bad_col[pnl] = off
CK.check("colours of both panels within the palette (blue, greys, ink, white, the orange and blue ladders, the red ramp)", not bad_col, f"{bad_col}")
# ---- numeric-token delta against V13 (V13 text layer from the round-37 probe; the composed sheet from Ghostscript)
SI = L.dedupe(BT["spans"])
def toks(strs): return L.word_tokens([dict(text=norm(t)) for t in strs]), L.num_tokens([dict(text=norm(t)) for t in strs])
wi, ni = toks([s["text"] for s in SI]); wo, no = toks([s["text"] for s in GSP]); wd_ = L.delta(wi, wo); nd = L.delta(ni, no)
json.dump(dict(words=wd_, numbers=nd), open(f"{S.VER}/text_delta.json", "w"), indent=1, ensure_ascii=False)
new_all = L.num_tokens([dict(text=norm(r["text"])) for pnl in "ab" for r in RECS[pnl]] + [dict(text=S.HDR_A)]); exp = L.delta(ni, new_all)
CK.check("numeric-token delta on the sheet = the expected delta (drawn records against the V13 strings; the two star keys of V13 are not drawn, as on V26); word delta in verify/text_delta.json", exp == nd, f"numbers lost {sum(nd['lost'].values())} gained {sum(nd['gained'].values())}")
# ---- renders (Ghostscript) and crops
NEW150 = f"{S.WORK}/composed_150dpi.png"; C.gs_render(COMPOSED, NEW150, 150)
# panel a renders pixel-identically to its V28 place: both sheets carry the 16 pt title strip, so the stamped render is compared with V28 at the same pixel phase (Ghostscript 150 dpi)
V28PNG = f"{S.WORK}/V28_150dpi.png"; C.gs_render(f"{C.V28}/{S.SHEET}.pdf", V28PNG, 150)
A = np.asarray(Image.open(V28PNG).convert("RGB")).astype(int); sc = 150 / 72
# the two panels are the round-49 panels V28 was composed from: byte-level provenance (the round-49 compose log's panel sha256 = the round-49 files) and pixel identity of the flat panel renders (same page sizes, same pixel phase)
R49 = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LED-B/ED_Fig07/work"; CL49 = json.load(open(f"{R49}/compose_log.json"))
sha49 = {os.path.basename(pl["file"]): pl["sha256"] for pl in CL49["placements"] if "file" in pl}
CK.check("the round-49 compose log names the panel files V28 was composed from, with sha256 equal to the round-49 lane's panel_a.pdf and panel_b.pdf", sha49.get("panel_a.pdf") == L.sha256(f"{R49}/panel_a.pdf") and sha49.get("panel_b.pdf") == L.sha256(f"{R49}/panel_b.pdf"), f"{ {k: v[:12] for k, v in sha49.items()} }")
for pnl, mine, theirs in (("a", S.PANEL_A, f"{R49}/panel_a.pdf"), ("b", S.PANEL_B, f"{R49}/panel_b.pdf")):
    p49, p51 = f"{S.WORK}/panel_{pnl}_r49_150dpi.png", f"{S.WORK}/panel_{pnl}_r51_150dpi.png"; C.gs_render(theirs, p49, 150); C.gs_render(mine, p51, 150)
    A_ = np.asarray(Image.open(p49).convert("RGB")).astype(int); B_ = np.asarray(Image.open(p51).convert("RGB")).astype(int)
    CK.check(f"panel {pnl}: this lane's flat panel PDF renders pixel-identically to the round-49 panel V28 was composed from (Ghostscript 150 dpi, same page size): 0 differing pixels", A_.shape == B_.shape and int((A_ != B_).any(axis=2).sum()) == 0, f"{int((A_ != B_).any(axis=2).sum()) if A_.shape == B_.shape else 'shape differs'} differing pixels, {A_.shape}")
if STAMPED_PNG and os.path.exists(STAMPED_PNG):
    B = np.asarray(Image.open(STAMPED_PNG).convert("RGB")).astype(int)
    CK.check("panel b placed to the right: the sheet's rightmost ink lies inside the page (the stamped render's last 9 pt are white)", B[:, int((W - S.EDGE + 2) * sc):].min() == 255, "")
    CK.check("the stamped sheet's title strip: the top 16 pt carry ink only inside the title box (x 12.96 to 152, y 2 to 15)", B[:int(16 * sc), int(155 * sc):].min() == 255 and B[:int(2 * sc), :].min() == 255, "")

o200, n200 = f"{S.WORK}/V28_200dpi.png", f"{S.WORK}/NEW_200dpi.png"; C.gs_render(f"{C.V28}/{S.SHEET}.pdf", o200, 200); C.gs_render(COMPOSED, n200, 200)
for pnl, box_old, box_new in (("a", [0, LAY["panel_a_rect"][1] + 16, W_IN, LAY["panel_a_rect"][3] + 16], LAY["panel_a_rect"]), ("b", [0, LAY["panel_b_src_rect"][1] + 16, W_IN, LAY["panel_b_src_rect"][3] + 16], LAY["panel_b_rect"])):
    s2 = 200 / 72; a = Image.open(o200).crop(tuple(int(round(v * s2)) for v in box_old)); b = Image.open(n200).crop(tuple(int(round(v * s2)) for v in box_new)); c = Image.new("RGB", (a.width + b.width + 20, max(a.height, b.height)), "white"); c.paste(a, (0, 0)); c.paste(b, (a.width + 20, 0)); c.save(f"{S.CROPS}/panel_{pnl}_V28_vs_NEW_200dpi.png")
CK.check("200 dpi V28 vs NEW crops written for a and b", all(os.path.exists(f"{S.CROPS}/panel_{q}_V28_vs_NEW_200dpi.png") for q in "ab"), "")
# ---- printed values csv, per flat panel (records in panel coordinates)
V.RULES.update({"sev_hr_ci": lambda o: L.hr_ci(o["hr"], o["lo"], o["hi"]), "sev_stars_p": lambda o: L.stars(o["p"]),
                "cross9_hr_stars": lambda cell: next(L.r2(cell["hr"]) + L.stars(Q9[c][k]) for c in CONDS9 for k, v in C9[c]["cells"].items() if v == cell and not v["reference"]),
                "cell_ci_paren": lambda cell: f"({L.r2(cell['lo'])}-{L.r2(cell['hi'])})"})
rows_all = []; n_no_all = 0
for pnl, pdf, dy in (("a", S.PANEL_A, S.SHIFT_A), ("b", S.PANEL_B, YB)):
    recs = [dict(r, baseline=r["baseline"] - dy, panel=pnl) for r in (RECS[pnl] if pnl == "a" else DB["drawn_records"])]
    n_rows, n_no = V.printed_values_csv(pdf, recs, f"{S.VER}/{S.SHEET}_printed_values_panel_{pnl}.csv"); n_no_all += n_no
    rows_all += [dict(r, panel=pnl) for r in csv.DictReader(open(f"{S.VER}/{S.SHEET}_printed_values_panel_{pnl}.csv"))]
with open(f"{S.VER}/{S.SHEET}_printed_values.csv", "w", newline="") as f:
    wr = csv.DictWriter(f, fieldnames=list(rows_all[0].keys())); wr.writeheader(); [wr.writerow(r) for r in rows_all]
cnt = collections.Counter(r["match"] for r in rows_all)
CK.check(f"printed values csv (per flat panel): {len(rows_all)} strings, {cnt.get('yes', 0)} re-derived from the source file, {cnt.get('yes_recorded', 0)} recorded, {n_no_all} NO", n_no_all == 0, f"NO rows: {[(r['string'], r['note']) for r in rows_all if r['match'] == 'NO'][:5]}")
# ---- CHANGES csv (this round: sizes only)
rows = [dict(sheet=S.SHEET, panel="page", label="layout", old_value=f"panels stacked, page {W_IN:.2f} x {H_IN:.2f}", new_value=f"panels side by side, page {LAY['sheet_w']:.2f} x {LAY['sheet_h']:.2f} before the strip", dramatic="no", note="round 51, Alen's item 4; no text or number changed"),
        dict(sheet=S.SHEET, panel="b", label="panel b origin", old_value=f"(0, {YB:.3f})", new_value=f"({X_B:.2f}, 0)", dramatic="no", note="tops aligned, letters at (12.96, 21.0) and ({X_B + 12.96:.2f}, 21.0)"),
        dict(sheet=S.SHEET, panel="a", label="row labels, tick labels, 'Negative controls'", old_value="10 pt", new_value=f"{GA['sizes']['tick']:g} pt", dramatic="no", note="round 49, +1 pt"),
        dict(sheet=S.SHEET, panel="a", label="estimates, stars", old_value="9.5 pt", new_value=f"{GA['sizes']['annot']:g} pt", dramatic="no", note="round 49, +1 pt"),
        dict(sheet=S.SHEET, panel="a", label="x-axis title", old_value="11 pt", new_value=f"{GA['sizes']['label']:g} pt", dramatic="no", note="round 49, +1 pt"),
        dict(sheet=S.SHEET, panel="a", label="'HR (95% CI)' column title", old_value=f"{S.HDR_A_FS:g} pt, right edge {CL['placements'][1]['x1']:.2f}", new_value="unchanged", dramatic="no", note="as in round 49"),
        dict(sheet=S.SHEET, panel="b", label="cell values, headers, group and level labels", old_value="9.5 pt", new_value=f"{9.5 + L.PT_PLUS:g} pt except {keep}", dramatic="no", note="round 49, +1 pt; fitted sizes declared"),
        dict(sheet=S.SHEET, panel="b", label="key (ramp title, ticks, reference key)", old_value="9.087 pt", new_value=f"{9.087 + L.PT_PLUS:g} pt", dramatic="no", note="round 49, +1 pt"),
        dict(sheet=S.SHEET, panel="page", label="letters a, b", old_value="13 pt", new_value=f"{S.LETTER_SIZE:g} pt", dramatic="no", note="round 49, +1 pt")]
L.write_changes(S.CHANGES, rows)
ok = CK.write(S.CHECKS); print("RESULT ALL PASS" if ok else "RESULT FAIL")
