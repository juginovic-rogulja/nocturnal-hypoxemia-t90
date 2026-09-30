#!$T90_PY
"""ROUND 49 (2026-09-26, lane L4), Main_Fig4 verify of the +1 pt vector compose (work/Main_Fig4_r49_vector.pdf).
Rules of the round: Ghostscript only on the composed sheet (pdfinfo, txtwrite, png16m), PyMuPDF only on the four flat panel pages.
Checks: (1) the compose log matches the files (output, band = the coordinator's band_4e_r49.pdf, the four panel pages against the build
manifest); (2) page size = V26's (pypdf mediabox, 0.05 pt); (3) the text census against the round-44 vector (the V26 source): same
multiset of line strings outside the band, every size exactly +1.0 except the declared keep box (panel b's printed columns), the band's
words the same and its sizes listed; (4) letters a to e and the title 14 pt at the V13 origins; (5) fonts Arial only; (6) every drawn
value re-derived from the numbers files (the round-38 idiom); (7) every drawn record string on the composed sheet's text layer at its
origin (1.6 pt); (8) every string inside its panel page, no string touching a neighbour (panels c and d cell text, all panels);
(9) geometry pins (c and d grids, panel b's plot box and column edges, the letter e clearance); (10) a fresh 150 dpi Ghostscript render.
Writes verify/checks.txt (ends RESULT ALL PASS), verify/census/census.json, verify/verify_record.json."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import collections, csv, fcntl, html, json, math, os, re, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import l4lib as L, lmain_lib as LM, census_r49 as CS
import fitz, pypdf
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND37_2026-09-14/figures/_common"); import v14lib as V
from statsmodels.stats.multitest import multipletests

SHEET = "Main_Fig4"; SD = f"{L.LANE}/{SHEET}"; WORK = f"{SD}/work"; VER = f"{SD}/verify"; os.makedirs(VER, exist_ok=True)
OUT = f"{WORK}/{SHEET}_r49_vector.pdf"; PNG150 = f"{SD}/{SHEET}_150dpi.png"
V26 = f"{paths.FIGURE_ROOT}/ROUND47_2026-09-24/figures/NEW_FINAL_SET_V26/{SHEET}.pdf"
OLDVEC = f"{paths.FIGURE_ROOT}/ROUND44_2026-09-21/lanes/LSKETCH/{SHEET}/work/{SHEET}_r44_vector.pdf"
BAND = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LSKETCH/work/build/band_4e_r49.pdf"
KEEP_BOX = (840.0, 34.927, 969.585, 716.9)    # panel b's printed columns (HR (95% CI), q, their headers): kept at 10 pt, see 01_build_b.py
G = json.load(open(f"{WORK}/v13_geometry.json")); CL = json.load(open(f"{WORK}/compose_log.json")); BM = json.load(open(f"{WORK}/build_manifest_r49.json"))
DA, DB, DC, DD = [json.load(open(f"{VER}/{SHEET}_{x}_drawn.json")) for x in "abcd"]
GEO = {x: json.load(open(f"{WORK}/panel_{x}_geometry.json")) for x in "abcd"}
BOX = {"a": G["panel_a"]["clip"], "b": G["panel_b"]["panel_box"], "c": G["panel_c"]["panel_box"], "d": GEO["d"]["page_box"]}
CK = L.Checks(f"Main_Fig4 round 49 (lane L4, {time.strftime('%Y-%m-%d %H:%M')}): every text one point larger. NEW = {OUT} (sha {L.sha256(OUT)[:16]}); "
              f"before = the round-44 vector {OLDVEC} (sha {L.sha256(OLDVEC)[:16]}, the source of V26's raster); band = {BAND} (sha {L.sha256(BAND)[:16]}).")
REC = {"sheet": SHEET, "out": OUT, "out_sha256": L.sha256(OUT), "old_vector": OLDVEC, "band": BAND, "band_sha256": L.sha256(BAND), "written": LM.now()}
E = 1.0   # the em unit in pt per pt of size
def norm(t): return t.replace("\xa0", " ")

# ---------------------------------------------------------------- (1) the compose log matches the files
CK.check("compose log: the sheet on disk is the composed one", CL["output_sha256"] == REC["out_sha256"] and CL["output"] == OUT, CL["output_sha256"][:16])
bp = [p for p in CL["placements"] if p["what"].startswith("band e")]; assert len(bp) == 1
CK.check("compose log: band e = the coordinator's round-49 band (sha256 of the file on disk), placed at the V13 band top", bp[0]["file"] == BAND and bp[0]["sha256"] == REC["band_sha256"] and abs(CL["band_top"] - G["band_e"]["rect"][1]) < 1e-6, f"{os.path.basename(BAND)} {REC['band_sha256'][:16]} top {CL['band_top']}")
okp = True; det = []
for x in "abcd":
    pp = [p for p in CL["placements"] if p["what"].startswith(f"panel {x} ")]; m = f"{WORK}/panel_{x}_m.pdf"; sha = L.sha256(m)
    okp &= len(pp) == 1 and pp[0]["sha256"] == sha == BM["panels"][x]["sha256"] and [round(v, 3) for v in pp[0]["rect"]] == [round(v, 3) for v in BOX[x]]; det.append(f"{x} {sha[:8]}")
CK.check("compose log: the four panel pages are this chain's builds (sha256 = build manifest = file) at their boxes (a, b, c the V13 boxes; d the V13 box widened 10 pt left)", okp, "; ".join(det))
CK.check("panel d page box = V13 box widened left by 10 pt only", GEO["d"]["v13_panel_box"] == G["panel_d"]["panel_box"] and abs(GEO["d"]["page_box"][0] - (G["panel_d"]["panel_box"][0] - 10.0)) < 1e-9 and GEO["d"]["page_box"][1:] == G["panel_d"]["panel_box"][1:], f"{GEO['d']['page_box']}")

# ---------------------------------------------------------------- (2) page size = V26
def mediabox(p):
    r = pypdf.PdfReader(L.hydrated(p)); mb = r.pages[0].mediabox; return len(r.pages), float(mb.width), float(mb.height)
n26, w26, h26 = mediabox(V26); nn, wn, hn = mediabox(OUT); ng, wg, hg = LM.gs_pdfinfo(OUT)
REC["page"] = dict(v26=[w26, h26], new=[wn, hn], new_gs=[wg, hg])
CK.check("page size = V26's (pypdf mediabox, 0.05 pt), one page (pypdf and gs PDFINFO agree)", n26 == nn == ng == 1 and abs(wn - w26) <= 0.05 and abs(hn - h26) <= 0.05 and abs(wg - wn) < 0.01 and abs(hg - hn) < 0.01, f"new {wn:.3f} x {hn:.3f}, V26 {w26:.3f} x {h26:.3f}")

# ---------------------------------------------------------------- (3) the text census (Ghostscript txtwrite, line level)
def REGION(s):   # the four panel regions of Main_Fig4 (panel d's widened page starts at 487.2; panel c's ink ends before 474)
    return ("a" if s["x0"] < 484.86 else "b") if s["y"] < 716.9 else ("c" if s["x0"] < 485.0 else "d")
rep, NSP, OSP = CS.census(OUT, OLDVEC, f"{VER}/census", CL["band_top"], [KEEP_BOX], region=REGION)
NSP_ALL = NSP; NSP = [s for s in NSP if not s["blank"]]
CK.check(f"census: the multiset of line strings outside the band is the same as on the round-44 vector ({rep['n_lines_outside_band']['old']} lines)", rep["strings_same"], f"lost {rep['lost'][:6]} gained {rep['gained'][:6]}")
CK.check(f"census: every size outside the band exactly +1.0 (sizes {rep['sizes_outside_band']['old']} -> {rep['sizes_outside_band']['new']}), EXCEPTION the keep box {KEEP_BOX} = panel b's printed columns kept at their size ({len(rep['kept_at_size'])} lines)", rep["sizes_ok"], f"missing {rep['size_missing_expected'][:6]} unexpected {rep['size_unexpected'][:6]}")
kept_sizes = collections.Counter(k[1] for k in rep["kept_at_size"])
CK.check("census: the kept lines are all 10 pt and all in panel b's column region (the declared exception, nothing else kept)", set(kept_sizes) == {10.0} and all(KEEP_BOX[0] <= k[2] and KEEP_BOX[1] <= k[3] <= KEEP_BOX[3] for k in rep["kept_at_size"]), f"{dict(kept_sizes)}")
CK.check("census: the band's words are the same as on the round-44 band (the coordinator's band: sizes listed, the two ramp-end labels on two lines)", rep["band"]["words_same"], f"sizes old {rep['band']['sizes']['old']} new {rep['band']['sizes']['new']}; lines old {len(rep['band']['lines_old'])} new {len(rep['band']['lines_new'])}")
CK.info(f"census: band lines before {[(t, s) for t, s, _, _ in rep['band']['lines_old']]}")
CK.info(f"census: band lines after  {[(t, s) for t, s, _, _ in rep['band']['lines_new']]}")
CK.info(f"census: fonts before {rep['fonts']['old']} after {rep['fonts']['new']}")

# ---------------------------------------------------------------- (4) letters and title, (5) fonts
lets = sorted([s for s in NSP if len(s["text"].strip()) == 1 and s["text"].strip() in "abcde" and s["size"] >= 12.5], key=lambda s: (s["y"], s["x0"]))
okl = [l["text"].strip() for l in lets] == list("abcde") and all(abs(l["x0"] - G["letters"][l["text"].strip()]["x"]) <= 1.0 and abs(l["y"] - G["letters"][l["text"].strip()]["baseline"]) <= 1.0 and abs(l["size"] - 14.0) < 0.01 and "Bold" in l["font"] for l in lets)
CK.check("letters a to e once each, Arial Bold 14 pt (were 13) at the V13 origins (1 pt, gs txtwrite)", okl, f"{[(l['text'].strip(), l['x0'], l['y'], l['size']) for l in lets]}")
ti = [s for s in NSP if norm(s["text"]).strip() == "Figure 4"]
CK.check("title 'Figure 4' once, Arial Bold 14 pt (was 13) at (14.17, 14.0)", len(ti) == 1 and abs(ti[0]["x0"] - G["letters"]["a"]["x"]) <= 1.0 and abs(ti[0]["y"] - LM.TITLE_BASELINE) <= 1.0 and abs(ti[0]["size"] - 14.0) < 0.01 and "Bold" in ti[0]["font"], f"{[(t['x0'], t['y'], t['size'], t['font']) for t in ti]}")
CK.check("fonts: Arial faces only on the composed sheet", all("Arial" in f for f in rep["fonts"]["new"]), f"{rep['fonts']['new']}")

minsz = min(s["size"] for s in NSP)
CK.check("no text below 10 pt anywhere on the sheet (the coordinator's final gate), band included", minsz >= 10.0 - 1e-6, f"smallest {minsz} pt; sizes present {sorted(set(s['size'] for s in NSP))}")
# ---------------------------------------------------------------- (6) values against the numbers files (re-derived here, the round-38 idiom)
X = json.load(open(L.hydrated(f"{L.NUM}/alenfig3_cross_v1.json"))); C6 = X["cross6"]; XT = json.load(open(L.hydrated(f"{L.NUM}/crosstab_v2.json"))); WS = json.load(open(L.hydrated(f"{L.NUM}/fig3_within_stratum_v1.json")))
for pth in (f"{L.NUM}/alenfig3_cross_v1.json", f"{L.NUM}/crosstab_v2.json", f"{L.NUM}/fig3_within_stratum_v1.json"): V.sidecar_v8(pth)
CONDS6 = [c for c, v in C6.items() if not v["negative_control"]]; CTRL6 = [c for c, v in C6.items() if v["negative_control"]]
Q6 = {c: {} for c in C6}
for cell in [f"{g}|{l}" for g in ("No or mild", "Moderate", "Severe") for l in ("<=1%", ">10%")]:
    if cell == DD["reference"]: continue
    for pool in (CONDS6, CTRL6):
        for c, q in zip(pool, multipletests([C6[c]["cells"][cell]["p"] for c in pool], method="fdr_bh")[1]): Q6[c][cell] = float(q)
bad_d = [(d["row"], d["cell"]) for d in DD["values"] if d["cell"] != DD["reference"] and (d["hr_text"] != L.r2(C6[d["row"]]["cells"][d["cell"]]["hr"]) + L.stars(Q6[d["row"]][d["cell"]]) or d["ci_text"] != f"({L.r2(C6[d['row']]['cells'][d['cell']]['lo'])}-{L.r2(C6[d['row']]['cells'][d['cell']]['hi'])})" or d["fill"] != L.ramp_hex(C6[d["row"]]["cells"][d["cell"]]["hr"]))]
CK.check("panel d: all 30 cells re-derived from alenfig3_cross_v1.json cross6 (half up 2 dp, BH stars within the cell, ramp fill) equal the drawn record", not bad_d and len(DD["values"]) == 30, f"{bad_d[:4]}")
bad_c = [(d["row"], d["col"]) for d in DC["values"] if d["text"] != f"{int(XT['rows3'][[r['category'] for r in XT['rows3']].index(d['row'])]['n_' + {'≤1%': '0-1%', '>1–5%': '1-5%', '>5–10%': '5-10%', '>10%': '>10%'}[d['col']]]):,} ({L.pct1(d['n'], d['group_n'])})"]
CK.check("panel c: all 12 cells re-derived from crosstab_v2.json rows3 (count, percent half up 1 dp) equal the drawn record", not bad_c and len(DC["values"]) == 12, f"{bad_c[:4]}")
AFULL = {r["cond"]: r for r in csv.DictReader(open(L.hydrated(f"{L.NUM}/wake_sleep_healthy_A_full.csv")))}
GA = V.sidecar_v8(f"{L.NUM}/wake_sleep_healthy_A_full.csv"); GJ = V.sidecar_v8(f"{L.NUM}/wake_sleep_healthy.json")
CK.check("panel a sources are the v8.2 re-extraction (step 121 sidecars after 2026-09-15 19:18), sha as recorded", GA["output_mtime"] >= "2026-09-15 19:18:00" and GJ["output_mtime"] >= "2026-09-15 19:18:00" and GA["sha256"] == DA["sources"]["A_full"]["sha256"], f"{GA['output_mtime']} / {GJ['output_mtime']}")
bad_a = [(d["condition"], d["period"]) for d in DA["values"] if d["hr_ci_text"] != L.hr_ci(AFULL[d["condition"]][d["period"][0] + "HR"], AFULL[d["condition"]][d["period"][0] + "lo"], AFULL[d["condition"]][d["period"][0] + "hi"]) or d["q_text"] != L.q_text(float(AFULL[d["condition"]][d["period"][0] + "q"])) or d["sig"] != (float(AFULL[d["condition"]][d["period"][0] + "q"]) < 0.05)]
CK.check("panel a: all 36 HR (95% CI) and q strings re-derived from wake_sleep_healthy_A_full.csv (half up 2 dp, q 3 dp, bold when q < 0.05) equal the drawn record", not bad_a and len(DA["values"]) == 36, f"{bad_a[:4]}")
_sig = {c for c, r in AFULL.items() if r["neg"] == "False" and c != "Death from any cause" and (float(r["wq"]) < 0.05 or float(r["sq"]) < 0.05)}
_top12 = sorted(_sig, key=lambda c: -float(AFULL[c]["sHR"]))[:12]
CK.check("panel a: the row set is the sheet's rule (12 largest sleep-period estimates among the conditions significant for either period, plus death, plus the five controls)", set(DA["rows_main"]) == set(_top12) | {"Death from any cause"} and set(DA["rows_controls"]) == {c for c, r in AFULL.items() if r["neg"] == "True"} and len(DA["rows_controls"]) == 5, f"main {DA['rows_main']}, controls {DA['rows_controls']}")
CK.check("panel a: rows ordered by the asleep-minus-awake gap within each block (the generator's rule)", all(float(AFULL[a]["sHR"]) - float(AFULL[a]["wHR"]) >= float(AFULL[b]["sHR"]) - float(AFULL[b]["wHR"]) - 1e-12 for blk in (DA["rows_main"], DA["rows_controls"]) for a, b in zip(blk, blk[1:])), "")
bad_b = [(d["condition"], d["arm"]) for d in DB["values"] if d["hr_ci_text"] != L.hr_ci(WS["fits"][d["condition"]][d["arm"]]["hr"], WS["fits"][d["condition"]][d["arm"]]["lo"], WS["fits"][d["condition"]][d["arm"]]["hi"]) or d["q_text"] != L.q_text(WS["fits"][d["condition"]][d["arm"]]["q"])]
CK.check(f"panel b: all {len(DB['values'])} HR (95% CI) and q strings re-derived from fig3_within_stratum_v1.json equal the drawn record", not bad_b and len(DB["values"]) == 2 * (len(DB["rows_conditions"]) + len(DB["rows_controls"])), f"{bad_b[:4]}")
CK.check("panel b: the row set is the file's own (panel_c_rows_no_or_mild + the controls it could fit), controls absent listed", DB["rows_conditions"] == WS["panel_c_rows_no_or_mild"] and DB["rows_controls"] == WS["panel_c_controls_no_or_mild"] and DB["controls_absent"] == WS["panel_c_controls_absent_no_or_mild"], f"{len(DB['rows_conditions'])} + {len(DB['rows_controls'])} rows, absent {DB['controls_absent']}")

# ---------------------------------------------------------------- (7) every drawn record on the composed sheet's text layer at its origin (gs txtwrite chars)
ROWS = collections.defaultdict(list)
for s in NSP_ALL:
    for c, x0, x1 in s["chars"]: ROWS[s["y"]].append((x0, x1, c))
ROWSTR = {}
for y, chs in ROWS.items():
    chs.sort(); ROWSTR[y] = ("".join(c for _, _, c in chs), chs)
def present(rec, tol=1.6):
    t = norm(rec["text"]); ha = rec.get("ha", "left")
    for y, (txt, chs) in ROWSTR.items():
        if abs(y - rec["baseline"]) > tol: continue
        i = txt.find(t)
        while i >= 0:
            x0 = chs[i][0]; x1 = chs[i + len(t) - 1][1]; ax = x0 if ha == "left" else (x1 if ha == "right" else (x0 + x1) / 2)
            if abs(ax - rec["x"]) <= tol: return True
            i = txt.find(t, i + 1)
    return False
RECS = {"a": DA["drawn_records"], "b": DB["drawn_records"], "c": DC["drawn_records"], "d": DD["drawn_records"]}
for pnl in "abcd":
    miss = [(r["text"], r["x"], r["baseline"]) for r in RECS[pnl] if not present(r)]
    CK.check(f"panel {pnl}: every drawn record string ({len(RECS[pnl])}) found on the composed sheet's text layer at its origin (gs txtwrite, 1.6 pt)", not miss, f"missing {miss[:6]}")
sizes_rec = {p: sorted(collections.Counter(round(r["size"], 3) for r in RECS[p]).items()) for p in "abcd"}
CK.info(f"drawn-record sizes after: {sizes_rec} (before: a 10/11, b 10/11, c 10/11, d 9.5/9.087)")

# ---------------------------------------------------------------- (8) flat panel pages: strings inside the page, nothing touching a neighbour
ASC = 0.73; XH = 0.52; DSC = 0.21
from fontTools.ttLib import TTFont
_TT = {"Arial": TTFont(L.ARIAL), "Arial Bold": TTFont(L.ARIAL_BOLD)}; _GB = {}
def glyph_box(c, font):
    """(yMin, yMax) of the character's outline in em units from the Arial face named by the span font (bold when 'Bold' in the name)."""
    face = "Arial Bold" if "Bold" in font else "Arial"; key = (face, c)
    if key not in _GB:
        tt = _TT[face]; gname = tt.getBestCmap().get(ord(c)); g = tt["glyf"][gname] if gname else None; upm = tt["head"].unitsPerEm
        _GB[key] = (0.0, 0.0) if g is None or g.numberOfContours == 0 else (g.yMin / upm, g.yMax / upm)
    return _GB[key]
def ink(s):
    t = s["text"]; size = s["size"]; y = s["origin"][1]
    top = y - (ASC if any(ch.isupper() or ch.isdigit() or ch in "()bdfhklt*%/[]{}|'\"≤≥<>" for ch in t) else XH) * size
    bot = y + (DSC if any(ch in "gjpqy(),;[]{}/" for ch in t) else 0.0) * size
    return top, bot
VIS = {}; touch = []; mingap = {}
for x in "abcd":
    pdoc = fitz.open(f"{WORK}/panel_{x}_m.pdf"); pg = pdoc[0]; pr = pg.rect; bx = BOX[x]
    assert abs(pr.width - (bx[2] - bx[0])) < 0.02 and abs(pr.height - (bx[3] - bx[1])) < 0.02, (x, pr, bx)
    sp = [s for s in L.spans_of(pg) if s["text"].strip()]; pdoc.close()
    for s in sp: s["page_bbox"] = [s["bbox"][0] + bx[0], s["bbox"][1] + bx[1], s["bbox"][2] + bx[0], s["bbox"][3] + bx[1]]; s["origin"] = [s["origin"][0] + bx[0], s["origin"][1] + bx[1]]
    VIS[x] = sp
    outside = [(s["text"], [round(v, 2) for v in s["page_bbox"]]) for s in sp if s["page_bbox"][0] < bx[0] - 0.3 or s["page_bbox"][2] > bx[2] + 0.3 or s["page_bbox"][1] < bx[1] - 0.3 or s["page_bbox"][3] > bx[3] + 0.3]
    CK.check(f"panel {x}: every string inside its panel page (nothing clipped at the page edge; {len(sp)} spans, PyMuPDF on the flat page)", not outside, f"outside {outside[:4]}")
    # horizontal neighbours on one baseline
    gaps = []
    ss = sorted(sp, key=lambda s: (round(s["origin"][1], 1), s["page_bbox"][0]))
    for i, s in enumerate(ss):
        for t in ss[i + 1:]:
            if abs(t["origin"][1] - s["origin"][1]) > 0.6: continue
            if t["page_bbox"][0] >= s["page_bbox"][2] - 0.01: gaps.append((round(t["page_bbox"][0] - s["page_bbox"][2], 2), s["text"], t["text"])); break
            if s["page_bbox"][0] < t["page_bbox"][2] and t["page_bbox"][0] < s["page_bbox"][2] and not (s["text"].endswith(" ") or t["text"].startswith(" ")):
                # PyMuPDF splits a single string at kerning points into adjacent spans (no gap): those share the drawn record's text
                joined = any(norm(s["text"]) + norm(t["text"]) in norm(r["text"]) or norm(t["text"]) + norm(s["text"]) in norm(r["text"]) for r in RECS[x])
                if not joined: touch.append((x, s["text"], t["text"], round(t["page_bbox"][0] - s["page_bbox"][2], 2)))
    # a same-record kerning split shows as a gap near 0: exclude pairs that belong to one record
    real = [g for g in gaps if not any(norm(g[1]) + norm(g[2]) in norm(r["text"]) for r in RECS[x])]
    mingap[x] = min(real) if real else None
    # vertical neighbours, glyph by glyph (PyMuPDF rawdict chars of the flat page): for every pair of x-overlapping characters on two
    # baselines less than 2.2 em apart, the clearance between the upper character's ink bottom (descender only for g j p q y and the like)
    # and the lower character's ink top (cap or ascender height for tall glyphs, x-height otherwise); two-line labels, the two lines of a cell
    pdoc = fitz.open(f"{WORK}/panel_{x}_m.pdf"); pg = pdoc[0]; CH = []
    for b_ in pg.get_text("rawdict")["blocks"]:
        if b_["type"] != 0: continue
        for l_ in b_["lines"]:
            for s_ in l_["spans"]:
                for ch in s_["chars"]:
                    if ch["c"].strip(): CH.append(dict(c=ch["c"], x0=ch["bbox"][0] + bx[0], x1=ch["bbox"][2] + bx[0], y=ch["origin"][1] + bx[1], size=s_["size"], font=s_["font"]))
    pdoc.close()
    rows_ = collections.defaultdict(list)
    for ch in CH: rows_[round(ch["y"], 1)].append(ch)
    ys_ = sorted(rows_); vgaps = []
    for i_, y1 in enumerate(ys_):
        for y2 in ys_[i_ + 1:]:
            if y2 - y1 > 2.2 * max(c["size"] for c in rows_[y1] + rows_[y2]): break
            for c1 in rows_[y1]:
                bot1 = c1["y"] + (DSC * c1["size"] if c1["c"] in "gjpqy(),;[]{}/" else 0.0)
                bot1 = c1["y"] - glyph_box(c1["c"], c1["font"])[0] * c1["size"]   # the glyph's own ink bottom (yMin of the outline, Arial or Arial Bold)
                for c2 in rows_[y2]:
                    if c2["x0"] < c1["x1"] - 0.3 and c1["x0"] < c2["x1"] - 0.3:
                        top2 = c2["y"] - glyph_box(c2["c"], c2["font"])[1] * c2["size"]   # the glyph's own ink top (yMax)
                        vgaps.append((round(top2 - bot1, 2), f"{c1['c']} of the line at {y1}", f"{c2['c']} of the line at {y2}"))
    vmin = min(vgaps) if vgaps else None
    CK.check(f"panel {x}: no string touches a horizontal neighbour on its baseline (min gap {mingap[x][0] if mingap[x] else None} pt between {mingap[x][1:] if mingap[x] else None}, gate 2.5) and no stacked glyphs overlap in ink (min vertical clearance {vmin[0] if vmin else None} pt between {vmin[1:] if vmin else None}, glyph outlines, gate above 0)", not [t for t in touch if t[0] == x] and (mingap[x] is None or mingap[x][0] >= 2.5) and (vmin is None or vmin[0] > 0.0), f"touching {[t for t in touch if t[0] == x][:4]}")
    REC[f"gaps_{x}"] = dict(min_horizontal=mingap[x], min_vertical_ink=vmin, n_spans=len(sp))
# panels c and d cell text: every cell string within its own cell rectangle (never over a neighbour cell)
for x, D in (("c", DC), ("d", DD)):
    bad = []
    for d in D["values"]:
        r = d["rect"]; texts = [d["text"]] if x == "c" else [d["hr_text"]] + ([d["ci_text"]] if d["ci_text"] else [])
        for t in texts:
            hits = [s for s in VIS[x] if norm(s["text"]) == t and r[1] - 1 <= s["origin"][1] <= r[3] + 1 and r[0] - 30 <= s["origin"][0] <= r[2] + 30]
            if not hits or any(s["page_bbox"][0] < r[0] + 0.5 or s["page_bbox"][2] > r[2] - 0.5 for s in hits): bad.append((t, [round(v, 1) for v in r], [[round(v, 1) for v in s["page_bbox"]] for s in hits]))
    CK.check(f"panel {x}: every cell string sits inside its own cell with at least 0.5 pt to the cell edge ({len(D['values'])} cells at +1 pt)", not bad, f"{bad[:3]}")

# the band (the coordinator's): stacked labels checked the same way from the composed sheet's txtwrite characters (integer boxes, 0.5 pt)
bch = []
for s_ in NSP_ALL:
    if s_["y"] < CL["band_top"]: continue
    for c, x0, x1 in s_["chars"]:
        if c.strip(): bch.append(dict(c=c, x0=x0, x1=x1, y=s_["y"], size=s_["size"], font=s_["font"]))
brows = collections.defaultdict(list)
for ch in bch: brows[ch["y"]].append(ch)
bys = sorted(brows); bgaps = []
for i_, y1 in enumerate(bys):
    for y2 in bys[i_ + 1:]:
        if y2 - y1 > 2.2 * max(c["size"] for c in brows[y1] + brows[y2]): break
        for c1 in brows[y1]:
            bot1 = c1["y"] - glyph_box(c1["c"], c1["font"])[0] * c1["size"]
            for c2 in brows[y2]:
                if c2["x0"] < c1["x1"] - 0.3 and c1["x0"] < c2["x1"] - 0.3:
                    bgaps.append((round(c2["y"] - glyph_box(c2["c"], c2["font"])[1] * c2["size"] - bot1, 2), f"{c1['c']} at y {y1}", f"{c2['c']} at y {y2}"))
bmin = min(bgaps) if bgaps else None
CK.check(f"band e (the coordinator's band_4e_r49.pdf): no stacked glyphs overlap in ink on the composed sheet (min vertical clearance {bmin[0] if bmin else None} pt between {bmin[1:] if bmin else None}; txtwrite integer boxes, 0.5 pt)", bmin is None or bmin[0] > 0.0, f"{sorted(bgaps)[:4]}")
REC["gaps_band"] = dict(min_vertical_ink=bmin)
# ---------------------------------------------------------------- (9) geometry pins
GC, GD_, GB = DC["geometry"], DD["geometry"], DB["geometry"]
CK.check("panel d grid touches the right margin (panel b's q column right edge, V13 956.621) and c and d grid bottoms are level (V13 1006.941)", abs(GD_["grid_x1"] - G["panel_b"]["qcol_right_x"]) < 0.01 and abs(GC["grid_bottom"] - GD_["grid_bottom"]) < 0.01 and abs(GC["grid_bottom"] - G["panel_c"]["grid_bottom"]) < 0.01, f"d x1 {GD_['grid_x1']:.3f}, bottoms {GC['grid_bottom']:.3f} / {GD_['grid_bottom']:.3f}")
ag = [r for r in DC["drawn_records"] if r["text"] == "apnea group"]; tk = [r for r in DD["drawn_records"] if r["text"] == "1" and r["baseline"] > 1000]
CK.check("panels c and d end on the same baseline (c's 'apnea group' = d's ramp ticks, V13 1058.208)", ag and tk and abs(ag[0]["baseline"] - tk[0]["baseline"]) <= 0.05 and abs(ag[0]["baseline"] - G["panel_c"]["bottom_ink_baseline"]) < 0.05, f"c {ag[0]['baseline'] if ag else None} d {tk[0]['baseline'] if tk else None}")
CK.check("panel b: plot box and printed columns at the V13 pins (plot x 628.859 to 833.385, rule on panel a's 666.886, HR column right 918.612, q column right 956.621)", GB["plot_box_page"] == [G["panel_b"]["plot_x"][0], G["panel_b"]["plot_top"], G["panel_b"]["plot_x"][1], G["panel_b"]["axis_rule_y"]] and GB["col_right_x"] == G["panel_b"]["col_right_x"] and GB["qcol_right_x"] == G["panel_b"]["qcol_right_x"], f"{GB['plot_box_page']}")
hr_sp = [s for s in DB["strings"] if re.fullmatch(r"\d+\.\d\d \(\d+\.\d\d-\d+\.\d\d\)", norm(s["text"]))]; q_sp = [s for s in DB["strings"] if norm(s["text"]).startswith("<0.") or re.fullmatch(r"0\.\d\d\d", norm(s["text"]))]
CK.check("panel b: the HR (95% CI) column right-aligned on the V13 edge and the q column on its edge (within 0.5 pt), both 10 pt (the declared exception)", hr_sp and q_sp and max(abs(s["x1"] - G["panel_b"]["col_right_x"]) for s in hr_sp) < 0.5 and max(abs(s["x1"] - G["panel_b"]["qcol_right_x"]) for s in q_sp) < 0.5 and {s["size"] for s in hr_sp + q_sp} == {10.0}, f"{len(hr_sp)} HR strings, {len(q_sp)} q strings")
lo_, hi_ = DB["xlim"]; PX = G["panel_b"]["plot_x"]
def xv(v): return PX[0] + (math.log(v) - math.log(lo_)) / (math.log(hi_) - math.log(lo_)) * (PX[1] - PX[0])
tip = max(xv(v["hi"]) for v in DB["values"]) + 0.85; wide = min(hr_sp, key=lambda s: s["x0"])
hq = [min(t["x0"] for t in q_sp if abs(t["y_bot"] - s["y_bot"]) < 0.5) - s["x1"] for s in hr_sp if any(abs(t["y_bot"] - s["y_bot"]) < 0.5 for t in q_sp)]
CK.check("panel b: the printed columns clear the plot ink and each other (widest HR string to the right-most whisker cap, HR to q on every sub-row, both at least 2.5 pt)", wide["x0"] - tip >= 2.5 and min(hq) >= 2.5, f"whisker cap {tip:.2f} to '{wide['text']}' at {wide['x0']:.2f}: {wide['x0'] - tip:.2f} pt; HR to q min {min(hq):.2f} pt")
CK.check("panel b: row labels keep the V26 line breaks (only 'Ventricular arrhythmia or cardiac arrest' on two lines) and the row set is unchanged", list(DB["row_labels_wrapped"]) == ["Ventricular arrhythmia or cardiac arrest"] and GB["n_rows"] == 21, f"{DB['row_labels_wrapped']}")
ea = [s for s in VIS["a"] if s["size"] == 12.0 and s["page_bbox"][0] < 100 and s["text"] != "Negative controls"]; ha_ = [s for s in VIS["a"] if s["size"] == 11.0 and abs(s["page_bbox"][2] - DA["geometry"]["col_right"]) < 0.6]; qa = [s for s in VIS["a"] if s["size"] == 11.0 and abs(s["page_bbox"][2] - DA["geometry"]["qcol_right"]) < 0.6]
gap_lab = DA["geometry"]["axis_x"][0] - max(s["page_bbox"][2] for s in ea); gap_hr = min(s["page_bbox"][0] for s in ha_) - DA["geometry"]["axis_x"][1]
hq_a = [min(t["page_bbox"][0] for t in qa if abs(t["origin"][1] - s["origin"][1]) < 0.5) - s["page_bbox"][2] for s in ha_ if any(abs(t["origin"][1] - s["origin"][1]) < 0.5 for t in qa)]
CK.check("panel a at +1 pt: row labels (12 pt) end before the axis start, the HR column (11 pt) starts after the axis end, HR to q at least 2.5 pt on every sub-row", gap_lab >= 1.5 and gap_hr >= 1.0 and min(hq_a) >= 2.5, f"label to axis {gap_lab:.2f} pt, axis end to HR column {gap_hr:.2f} pt, HR to q min {min(hq_a):.2f} pt")
ink_bottom = max(ink(s)[1] for x in "cd" for s in VIS[x]); e_top = G["letters"]["e"]["baseline"] - 0.905 * 14.0; e_glyph_top = G["letters"]["e"]["baseline"] - XH * 14.0
CK.check(f"letter e (14 pt) clears the lowest ink of panels c and d (box top at least 6 pt above the ink; V13 pin 8 pt at 13 pt), band top below the letter's baseline", e_top - ink_bottom >= 6.0 and CL["band_top"] > G["letters"]["e"]["baseline"], f"e box top {e_top:.2f} (glyph top {e_glyph_top:.2f}), c/d ink bottom {ink_bottom:.2f}, clearance {e_top - ink_bottom:.2f} pt, band top {CL['band_top']}")
CK.check("page bottom = band bottom (band top + band height), V13 height kept", abs(CL["band_top"] + bp[0]["band_page"][1] - hn) < 0.01 and abs(hn - G["page"][1]) < 0.01, f"{CL['band_top']} + {bp[0]['band_page'][1]} = {hn}")

# ---------------------------------------------------------------- (10) a fresh 150 dpi Ghostscript render (under the round-30 render lock, one heavy render at a time)
if os.path.exists(PNG150): os.remove(PNG150)
t0 = time.time()
with L.render_lock():
    r = subprocess.run([LM.GS, "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=png16m", "-r150", f"-sOutputFile={PNG150}", OUT], capture_output=True, text=True)
from PIL import Image
im = Image.open(PNG150) if os.path.exists(PNG150) else None; ew, eh = round(wn * 150 / 72), round(hn * 150 / 72)
CK.check("150 dpi Ghostscript render written fresh (png16m, render lock), pixel size = page size at 150 dpi (1 px)", r.returncode == 0 and im is not None and abs(im.width - ew) <= 1 and abs(im.height - eh) <= 1 and "error" not in (r.stdout + r.stderr).lower(), f"{PNG150} {im.size if im else None} in {time.time() - t0:.1f} s (expected {ew} x {eh})")
REC["png150"] = dict(path=PNG150, size=list(im.size) if im else None, sha256=L.sha256(PNG150) if im else None)
REC["census"] = f"{VER}/census/census.json"; REC["checks"] = f"{VER}/checks.txt"
json.dump(REC, open(f"{VER}/verify_record.json", "w"), indent=1, ensure_ascii=False)
ok = CK.write(f"{VER}/checks.txt"); sys.exit(0 if ok else 1)
