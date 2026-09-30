#!$T90_PY
"""ROUND 54 (2026-09-29, lane L4): the verifier of Main_Fig4 (a | b + the sketch band as panel c, owner's option A) and of the standalone
counts panel. Ghostscript only on the composed sheets (PDFINFO, txtwrite, png16m), PyMuPDF only on the flat panel pages.
Main_Fig4 checks: (1) the compose log matches the files (output, the coordinator's band_4e_r54.pdf, the panel pages = build manifest);
(2) page box declared, one page, at most 1110 pt (the aim) and 1260 (the gate), pypdf and gs agree, the flat twin the same box;
(3) the txtwrite census against V30's vector twin (the round-49 vector, the source of V30's raster): every string of V30's panels a and
b present at its round-54 size (14 pt, 15 pt for the six column headers, captions and axis-title lines, 18 pt for the letters and the
title), the declared re-wraps (labels on two lines) re-joined before the comparison, V30's panels c and d strings declared absent (listed),
the band's words the same as the round-49 band's, nothing below 12 pt (the smallest reported); (4) letters a, b, c and the title 18 pt
Arial Bold at the spec origins, Arial faces only; (5) every drawn value re-derived from the numbers files (the round-49 idiom) and the
marks preserved (the same hr, lo, hi, q per condition and arm as the round-49 record, every marker at xv(hr) on the new axis);
(6) every drawn record string on the composed sheet's text layer at its origin (txtwrite chars, 1.6 pt); (7) no overlapping text: on the
flat panel pages no string touches a horizontal neighbour (2.5 pt) and no stacked glyphs overlap in ink (fontTools outlines), on the whole
sheet no two txtwrite line boxes intersect, the band checked from the composed sheet; (8) placement rules: labels clear the plots, the
HR columns clear the whiskers, HR to q at least 2.5 pt, the two axis rules on one line, the letter c clear of the a/b ink and above the
band; (9) a fresh 150 dpi Ghostscript render; (10) the flat twin's census equals the nested compose's.
Counts panel checks: page box, strings = V30's panel c strings at 14 pt (title 15), values re-derived from crosstab_v2.json, cells inside
their rectangles, no overlaps, render, flat twin census. Writes <deliverable>/verify/checks.txt (RESULT ALL PASS), census.json, verify_record.json."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import collections, csv, json, math, os, re, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import l4lib as L, lmain_lib as LM, census_r49 as CS, layout_r54 as LY
import fitz, pypdf
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND37_2026-09-14/figures/_common"); import v14lib as V
from fontTools.ttLib import TTFont

SHEET = "Main_Fig4"; SD = f"{L.LANE}/{SHEET}"; WORK = f"{SD}/work"; VER = f"{SD}/verify"; os.makedirs(VER, exist_ok=True); os.makedirs(f"{VER}/census", exist_ok=True)
PC = f"{L.LANE}/panel_counts_r54"; PCV = f"{PC}/verify"; os.makedirs(f"{PCV}/census", exist_ok=True)
OUT = f"{SD}/{SHEET}.pdf"; FLAT = f"{WORK}/{SHEET}_flat.pdf"; PNG150 = f"{SD}/{SHEET}_150dpi.png"
OUT_C = f"{PC}/panel_counts_r54.pdf"; FLAT_C = f"{PC}/work/panel_counts_r54_flat.pdf"; PNG_C = f"{PC}/panel_counts_r54_150dpi.png"
V30 = f"{paths.FIGURE_ROOT}/ROUND52_2026-09-26/figures/NEW_FINAL_SET_V30/{SHEET}.pdf"
OLDVEC = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/L4/{SHEET}/work/{SHEET}_r49_vector.pdf"      # V30's vector twin (V30 = the 600-dpi raster of this file, ROUND49/lanes/LRASTER)
R49V = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/L4/{SHEET}/verify"
S = LY.load(); FS = S["sizes"]; BAND = S["e"]["band"]; BAND_R49 = LY.BAND_R49
CL = json.load(open(f"{WORK}/compose_log.json")); BM = json.load(open(f"{WORK}/build_manifest_r54.json")); CLC = json.load(open(f"{PC}/work/compose_log.json"))
DA, DB, DC = [json.load(open(f"{VER}/{SHEET}_{x}_drawn.json")) for x in "abc"]; GEO = {x: json.load(open(f"{WORK}/panel_{x}_geometry.json")) for x in "abc"}
RA, RB, RC49 = [json.load(open(f"{R49V}/{SHEET}_{x}_drawn.json")) for x in "abc"]
def norm(t): return t.replace("\xa0", " ")
HEAD15 = {"HR (95% CI)", "q", "Hazard ratio (95% CI)", "per 1 SD higher T90", "Hazard ratio (95% CI), T90 above 10%", "versus 10% or less, within each apnea group"}
def size_class(text):
    t = norm(text).strip()
    if t == "Figure 4" or (len(t) == 1 and t in "abc"): return 18.0
    return 15.0 if t in HEAD15 else 14.0
CK = L.Checks(f"Main_Fig4 round 54 (lane L4, {time.strftime('%Y-%m-%d %H:%M')}): a | b + the sketch band as panel c at the round-54 sizes (owner's option A). NEW = {OUT} (sha {L.sha256(OUT)[:16]}); "
              f"V30 = {V30} (sha {L.sha256(V30)[:16]}), its vector twin = {OLDVEC} (sha {L.sha256(OLDVEC)[:16]}); band = {BAND} (sha {L.sha256(BAND)[:16]}).")
REC = {"sheet": SHEET, "out": OUT, "out_sha256": L.sha256(OUT), "v30": V30, "v30_sha256": L.sha256(V30), "old_vector": OLDVEC, "band": BAND, "band_sha256": L.sha256(BAND), "written": LM.now(), "layout": S["arrangement"]}

# ---------------------------------------------------------------- (1) the compose log matches the files
CK.check("compose log: the sheet on disk is the composed one", CL["output_sha256"] == REC["out_sha256"] and CL["output"] == OUT, CL["output_sha256"][:16])
bp = [p for p in CL["placements"] if p["what"].startswith("band 4e")]; assert len(bp) == 1
CK.check("compose log: the band is the coordinator's round-54 band_4e_r54.pdf (sha256 of the file on disk), placed at the spec's band top, its page box read from the file", bp[0]["file"] == BAND and bp[0]["sha256"] == REC["band_sha256"] and abs(CL["band_top"] - S["e"]["band_top"]) < 1e-6 and abs(bp[0]["band_page"][1] - S["e"]["band_h"]) < 0.01, f"{os.path.basename(BAND)} {REC['band_sha256'][:16]} top {CL['band_top']:.3f}, band {bp[0]['band_page']}")
okp = True; det = []
for x in "ab":
    pp = [p for p in CL["placements"] if p["what"].startswith(f"panel {x} ")]; m = f"{WORK}/panel_{x}_m.pdf"; sha = L.sha256(m)
    okp &= len(pp) == 1 and pp[0]["sha256"] == sha == BM["panels"][x]["sha256"] and [round(v, 3) for v in pp[0]["rect"]] == [round(v, 3) for v in S[x]["box"]]; det.append(f"{x} {sha[:8]}")
CK.check("compose log: panels a and b are this chain's builds (sha256 = build manifest = file) at their spec boxes, V30's panels c and d not placed (owner's option A)", okp and not [p for p in CL["placements"] if p["what"].startswith("panel c") or p["what"].startswith("panel d")], "; ".join(det))

# ---------------------------------------------------------------- (2) page box
def mediabox(p):
    r = pypdf.PdfReader(L.hydrated(p)); mb = r.pages[0].mediabox; return len(r.pages), float(mb.width), float(mb.height)
n30, w30, h30 = mediabox(V30); nn, wn, hn = mediabox(OUT); ng, wg, hg = LM.gs_pdfinfo(OUT); nf, wf, hf = LM.gs_pdfinfo(FLAT)
REC["page"] = dict(v30=[w30, h30], new=[wn, hn], new_gs=[wg, hg], flat=[wf, hf], spec=[S["page_w"], S["page_h"]])
CK.check(f"page box declared: {wn:.3f} x {hn:.3f} pt, one page, width = V30's {w30:.3f} (0.05 pt), height at most 1110 (the aim) and 1260 (the gate), pypdf and gs PDFINFO agree, = the spec's", nn == ng == 1 and abs(wn - w30) <= 0.05 and hn <= 1110.0 and hn <= 1260.0 and abs(wg - wn) < 0.01 and abs(hg - hn) < 0.01 and abs(hn - S["page_h"]) < 0.01, f"V30 was {w30:.3f} x {h30:.3f}")
CK.check("the flat twin (gs pdfwrite re-distillation, work/Main_Fig4_flat.pdf): one page, the same page box", nf == 1 and abs(wf - wn) <= 0.05 and abs(hf - hn) <= 0.05 and CL["flat_sha256"] == L.sha256(FLAT), f"{wf} x {hf}")

# ---------------------------------------------------------------- (3) the census (Ghostscript txtwrite, line level)
AB_SPLIT = S["a"]["box"][2]; AB_BOTTOM = S["ab_bottom"]; BAND_TOP = S["e"]["band_top"]
nsp = CS.spans(CS.txtwrite(OUT, f"{VER}/census/new_txtwrite.xml")); osp = CS.spans(CS.txtwrite(OLDVEC, f"{VER}/census/old_txtwrite.xml"))
def reg_old(s): return ("a" if s["x0"] < 484.86 else "b") if s["y"] < 716.9 else (("c" if s["x0"] < 485.0 else "d") if s["y"] < 1083.38 else "band")
def reg_new(s): return ("a" if s["x0"] < AB_SPLIT else "b") if s["y"] < BAND_TOP else "band"
nl = CS.lines(nsp, reg_new); ol = CS.lines(osp, reg_old)
for l in nl: l["region"] = reg_new(dict(x0=l["x0"], y=l["y"]))
for l in ol: l["region"] = reg_old(dict(x0=l["x0"], y=l["y"]))
LETTERS = lambda l: l["text"].strip() in ("a", "b", "c", "d", "e", "Figure 4") and l["size"] >= 12.5
old_ab = [l for l in ol if l["region"] in "ab" and not LETTERS(l)]; old_cd = [l for l in ol if l["region"] in ("c", "d") and not LETTERS(l)]; old_band = [l for l in ol if l["region"] == "band"]
new_ab = [l for l in nl if l["region"] in "ab" and not LETTERS(l)]; new_band = [l for l in nl if l["region"] == "band"]
# the declared re-wraps: labels on two lines in round 54 that were one line on V30 (the two lines are re-joined for the string comparison)
WRAPS_OLD = {("a", k) for k in (RA.get("row_labels_wrapped", {}) or {})} | {("b", k) for k in RB["row_labels_wrapped"]}
REJOIN = [(pnl, k, list(v)) for pnl, W_ in (("a", DA["row_labels_wrapped"]), ("b", {k: v.split("\n") for k, v in DB["row_labels_wrapped"].items()})) for k, v in W_.items() if (pnl, k) not in WRAPS_OLD]
def rejoin(lines_):
    out = collections.Counter(l["text"] for l in lines_)
    for pnl, lab, (l1, l2) in REJOIN:
        assert out[l1] >= 1 and out[l2] >= 1, ("declared wrap not found on the new sheet", pnl, lab, l1, l2); out[l1] -= 1; out[l2] -= 1; out[lab] += 1
    return +out
n_str = rejoin(new_ab); o_str = collections.Counter(l["text"] for l in old_ab)
lost = sorted((o_str - n_str).elements()); gained = sorted((n_str - o_str).elements())
CK.check(f"census: every string of V30's panels a and b ({len(old_ab)} lines) is on the new sheet (line multiset equal after re-joining the {len(REJOIN)} declared re-wraps {[(p_, k_) for p_, k_, _ in REJOIN]}; V30's own two-line label kept)", n_str == o_str, f"lost {lost[:8]} gained {gained[:8]}")
CK.check(f"census: V30's panels c and d strings declared ABSENT ({len(old_cd)} lines, owner's option A: the counts panel is delivered standalone, panel d dropped), none of them on the new sheet outside the band except strings that panels a and b also carry", not [l for l in nl if l["region"] in "ab" and l["text"] in {x["text"] for x in old_cd} and l["text"] not in o_str], f"absent strings listed in census.json")
bad_size = [(l["text"], l["size"], size_class(l["text"])) for l in nl if l["region"] in "ab" and abs(l["size"] - size_class(l["text"])) > 0.01]
CK.check("census: every line outside the band at its round-54 size (14 pt, the eight header, caption and axis-title lines, six distinct strings, 15 pt, letters and the title 18 pt)", not bad_size and all(abs(l["size"] - 18.0) < 0.01 for l in nl if LETTERS(l)), f"wrong {bad_size[:6]}")
old_sizes = sorted(collections.Counter(l["size"] for l in old_ab).items()); new_sizes = sorted(collections.Counter(l["size"] for l in new_ab).items())
CK.info(f"census: sizes of V30's a and b lines {old_sizes} -> new {new_sizes}")
bw_old, bw_new = CS.words(old_band), CS.words(new_band)
CK.check("census: the band's words are the same as on the round-49 band (the coordinator's round-54 band, sizes listed)", bw_old == bw_new, f"lost {sorted((bw_old - bw_new).elements())[:6]} gained {sorted((bw_new - bw_old).elements())[:6]}; sizes old {sorted(collections.Counter(l['size'] for l in old_band).items())} new {sorted(collections.Counter(l['size'] for l in new_band).items())}")
minsz = min(l["size"] for l in nl)
CK.check("no text below 12 pt anywhere on the sheet (the round-54 floor), band included: the smallest outside the band is 14.0, the band's 14-pt caption reads 13.995 in txtwrite (the sketch engine's scaling, as the round-49 band's 10.995 did)", minsz >= 12.0 and min(l["size"] for l in nl if l["region"] != "band") >= 14.0 - 1e-6 and minsz >= 13.99, f"smallest {minsz} pt; sizes present {sorted(set(l['size'] for l in nl))}")
fonts_new = sorted(set(s["font"] for s in nsp)); CK.check("fonts: Arial faces only on the composed sheet", all("Arial" in f for f in fonts_new), f"{fonts_new}")
json.dump(dict(new=OUT, old=OLDVEC, band_top=BAND_TOP, lines_new=len(nl), lines_old=len(ol), old_ab=len(old_ab), new_ab=len(new_ab), strings_same=(n_str == o_str), lost=lost, gained=gained, rewraps=REJOIN,
               absent_cd=sorted(l["text"] for l in old_cd), absent_letters=["d", "e"], sizes_old_ab=old_sizes, sizes_new_ab=new_sizes, band_words_same=(bw_old == bw_new),
               band_lines_old=[(l["text"], l["size"]) for l in sorted(old_band, key=lambda l: (l["y"], l["x0"]))], band_lines_new=[(l["text"], l["size"]) for l in sorted(new_band, key=lambda l: (l["y"], l["x0"]))], fonts=fonts_new, min_size=minsz),
          open(f"{VER}/census/census.json", "w"), indent=1, ensure_ascii=False)

# ---------------------------------------------------------------- (4) letters and title
lets = sorted([l for l in nl if len(l["text"].strip()) == 1 and l["text"].strip() in "abcde" and l["size"] >= 12.5], key=lambda l: (l["y"], l["x0"]))
okl = [l["text"].strip() for l in lets] == list("abc") and all(abs(l["x0"] - S["letters"][l["text"].strip()]["x"]) <= 1.0 and abs(l["y"] - S["letters"][l["text"].strip()]["base"]) <= 1.0 and abs(l["size"] - 18.0) < 0.01 for l in lets)
fonts_lets = all("Bold" in s["font"] for s in nsp if s["text"].strip() in ("a", "b", "c") and len(s["text"].strip()) == 1 and s["size"] >= 12.5)
CK.check("letters a, b, c once each (no d, no e), Arial Bold 18 pt at the spec origins (1 pt, gs txtwrite)", okl and fonts_lets, f"{[(l['text'].strip(), l['x0'], l['y'], l['size']) for l in lets]}")
ti = [l for l in nl if norm(l["text"]).strip() == "Figure 4"]
CK.check("title 'Figure 4' once, Arial Bold 18 pt at (14.17, baseline 18)", len(ti) == 1 and abs(ti[0]["x0"] - S["title"]["x"]) <= 1.0 and abs(ti[0]["y"] - 18.0) <= 1.0 and abs(ti[0]["size"] - 18.0) < 0.01 and all("Bold" in s["font"] for s in nsp if norm(s["text"]).strip() == "Figure 4"), f"{[(t['x0'], t['y'], t['size']) for t in ti]}")

# ---------------------------------------------------------------- (5) values re-derived from the numbers files, marks preserved
WS = json.load(open(L.hydrated(f"{L.NUM}/fig3_within_stratum_v1.json"))); V.sidecar_v8(f"{L.NUM}/fig3_within_stratum_v1.json")
AFULL = {r["cond"]: r for r in csv.DictReader(open(L.hydrated(f"{L.NUM}/wake_sleep_healthy_A_full.csv")))}
GA = V.sidecar_v8(f"{L.NUM}/wake_sleep_healthy_A_full.csv"); GJ = V.sidecar_v8(f"{L.NUM}/wake_sleep_healthy.json")
CK.check("panel a sources are the v8.2 re-extraction (step 121 sidecars after 2026-09-15 19:18), sha as recorded", GA["output_mtime"] >= "2026-09-15 19:18:00" and GJ["output_mtime"] >= "2026-09-15 19:18:00" and GA["sha256"] == DA["sources"]["A_full"]["sha256"], f"{GA['output_mtime']} / {GJ['output_mtime']}")
bad_a = [(d["condition"], d["period"]) for d in DA["values"] if d["hr_ci_text"] != L.hr_ci(AFULL[d["condition"]][d["period"][0] + "HR"], AFULL[d["condition"]][d["period"][0] + "lo"], AFULL[d["condition"]][d["period"][0] + "hi"]) or d["q_text"] != L.q_text(float(AFULL[d["condition"]][d["period"][0] + "q"])) or d["sig"] != (float(AFULL[d["condition"]][d["period"][0] + "q"]) < 0.05)]
CK.check("panel a: all 36 HR (95% CI) and q strings re-derived from wake_sleep_healthy_A_full.csv (half up 2 dp, q 3 dp, bold when q < 0.05) equal the drawn record", not bad_a and len(DA["values"]) == 36, f"{bad_a[:4]}")
_sig = {c for c, r in AFULL.items() if r["neg"] == "False" and c != "Death from any cause" and (float(r["wq"]) < 0.05 or float(r["sq"]) < 0.05)}; _top12 = sorted(_sig, key=lambda c: -float(AFULL[c]["sHR"]))[:12]
CK.check("panel a: the row set is the sheet's rule (12 largest sleep-period estimates among the conditions significant for either period, plus death, plus the five controls), ordered by the asleep-minus-awake gap within each block", set(DA["rows_main"]) == set(_top12) | {"Death from any cause"} and set(DA["rows_controls"]) == {c for c, r in AFULL.items() if r["neg"] == "True"} and all(float(AFULL[a]["sHR"]) - float(AFULL[a]["wHR"]) >= float(AFULL[b]["sHR"]) - float(AFULL[b]["wHR"]) - 1e-12 for blk in (DA["rows_main"], DA["rows_controls"]) for a, b in zip(blk, blk[1:])), f"main {DA['rows_main']}, controls {DA['rows_controls']}")
bad_b = [(d["condition"], d["arm"]) for d in DB["values"] if d["hr_ci_text"] != L.hr_ci(WS["fits"][d["condition"]][d["arm"]]["hr"], WS["fits"][d["condition"]][d["arm"]]["lo"], WS["fits"][d["condition"]][d["arm"]]["hi"]) or d["q_text"] != L.q_text(WS["fits"][d["condition"]][d["arm"]]["q"])]
CK.check(f"panel b: all {len(DB['values'])} HR (95% CI) and q strings re-derived from fig3_within_stratum_v1.json equal the drawn record", not bad_b and len(DB["values"]) == 2 * (len(DB["rows_conditions"]) + len(DB["rows_controls"])) == 42, f"{bad_b[:4]}")
CK.check("panel b: the row set is the file's own (panel_c_rows_no_or_mild + the controls it could fit), controls absent listed", DB["rows_conditions"] == WS["panel_c_rows_no_or_mild"] and DB["rows_controls"] == WS["panel_c_controls_no_or_mild"] and DB["controls_absent"] == WS["panel_c_controls_absent_no_or_mild"], f"{len(DB['rows_conditions'])} + {len(DB['rows_controls'])} rows, absent {DB['controls_absent']}")
# marks preserved: the same values as the round-49 record (V30's marks), every marker at xv(hr) on the new axis
def xv_a(v): g = GEO["a"]; return g["axis_x"][0] + (math.log(v) - math.log(g["xlim"][0])) / (math.log(g["xlim"][1]) - math.log(g["xlim"][0])) * (g["axis_x"][1] - g["axis_x"][0])
def xv_b(v): g = GEO["b"]; return g["plot_box_page"][0] + (math.log(v) - math.log(g["xlim"][0])) / (math.log(g["xlim"][1]) - math.log(g["xlim"][0])) * (g["plot_box_page"][2] - g["plot_box_page"][0])
o49a = {(d["condition"], d["period"]): d for d in RA["values"]}; o49b = {(d["condition"], d["arm"]): d for d in RB["values"]}
mk_a = [k for k in [(d["condition"], d["period"]) for d in DA["values"]] if k not in o49a or any(abs(o49a[k][f] - next(d for d in DA["values"] if (d["condition"], d["period"]) == k)[f]) > 1e-12 for f in ("hr", "lo", "hi", "q"))]
mk_b = [k for k in [(d["condition"], d["arm"]) for d in DB["values"]] if k not in o49b or any(abs(o49b[k][f] - next(d for d in DB["values"] if (d["condition"], d["arm"]) == k)[f]) > 1e-12 for f in ("hr", "lo", "hi", "q"))]
pos_a = max(abs(d["x_page"] - xv_a(d["hr"])) for d in DA["values"]); pos_b = max(abs(d["x_page"] - xv_b(d["hr"])) for d in DB["values"])
CK.check("marks preserved: every marker and whisker of panels a and b carries the same hr, lo, hi and q as V30's (the round-49 record), the same significance fill, the same series colours, every marker at xv(hr) on its new axis (0.01 pt), axis ranges and ticks as round 49", not mk_a and not mk_b and len(DA["values"]) == len(RA["values"]) and len(DB["values"]) == len(RB["values"]) and pos_a < 0.01 and pos_b < 0.01 and GEO["a"]["xlim"] == RA["geometry"]["xlim"] and GEO["a"]["ticks"] == RA["geometry"]["ticks"] and GEO["b"]["xlim"] == RB["geometry"]["xlim"] and GEO["b"]["ticks"] == RB["geometry"]["ticks"] and GEO["a"]["colours"] == RA["geometry"]["colours"] and GEO["b"]["colours"] == RB["geometry"]["colours"], f"value diffs a {mk_a[:3]} b {mk_b[:3]}, position residual a {pos_a:.4f} b {pos_b:.4f}")

# ---------------------------------------------------------------- (6) every drawn record on the composed sheet's text layer at its origin
def rowstr(spans_):
    rows = collections.defaultdict(list)
    for s in spans_:
        for c, x0, x1 in s["chars"]: rows[s["y"]].append((x0, x1, c))
    out = {}
    for y, chs in rows.items(): chs.sort(); out[y] = ("".join(c for _, _, c in chs), chs)
    return out
ROWSTR = rowstr(nsp)
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
for pnl, D in (("a", DA), ("b", DB)):
    miss = [(r["text"], r["x"], r["baseline"]) for r in D["drawn_records"] if not present(r)]
    CK.check(f"panel {pnl}: every drawn record string ({len(D['drawn_records'])}) found on the composed sheet's text layer at its origin (gs txtwrite, 1.6 pt)", not miss, f"missing {miss[:6]}")
CK.info(f"drawn-record sizes: a {sorted(collections.Counter(r['size'] for r in DA['drawn_records']).items())}, b {sorted(collections.Counter(r['size'] for r in DB['drawn_records']).items())}")

# ---------------------------------------------------------------- (7) no overlapping text (flat pages: glyph outlines; whole sheet: txtwrite line boxes)
_TT = {"Arial": TTFont(L.ARIAL), "Arial Bold": TTFont(L.ARIAL_BOLD)}; _GB = {}
def glyph_box(c, font):
    face = "Arial Bold" if "Bold" in font else "Arial"; key = (face, c)
    if key not in _GB:
        tt = _TT[face]; gname = tt.getBestCmap().get(ord(c)); g = tt["glyf"][gname] if gname else None; upm = tt["head"].unitsPerEm
        _GB[key] = (0.0, 0.0) if g is None or g.numberOfContours == 0 else (g.yMin / upm, g.yMax / upm)
    return _GB[key]
def flat_page_checks(x, pdf, box, recs, D=None, cells=None):
    pdoc = fitz.open(pdf); pg = pdoc[0]; pr = pg.rect; bx = box
    assert abs(pr.width - (bx[2] - bx[0])) < 0.02 and abs(pr.height - (bx[3] - bx[1])) < 0.02, (x, pr, bx)
    sp = [s for s in L.spans_of(pg) if s["text"].strip()]
    for s in sp: s["page_bbox"] = [s["bbox"][0] + bx[0], s["bbox"][1] + bx[1], s["bbox"][2] + bx[0], s["bbox"][3] + bx[1]]; s["origin"] = [s["origin"][0] + bx[0], s["origin"][1] + bx[1]]
    outside = [(s["text"], [round(v, 2) for v in s["page_bbox"]]) for s in sp if s["page_bbox"][0] < bx[0] - 0.3 or s["page_bbox"][2] > bx[2] + 0.3 or s["page_bbox"][1] < bx[1] - 0.3 or s["page_bbox"][3] > bx[3] + 0.3]
    gaps = []; touch = []; ss = sorted(sp, key=lambda s: (round(s["origin"][1], 1), s["page_bbox"][0]))
    for i, s in enumerate(ss):
        for t in ss[i + 1:]:
            if abs(t["origin"][1] - s["origin"][1]) > 0.6: continue
            if t["page_bbox"][0] >= s["page_bbox"][2] - 0.01: gaps.append((round(t["page_bbox"][0] - s["page_bbox"][2], 2), s["text"], t["text"])); break
            if s["page_bbox"][0] < t["page_bbox"][2] and t["page_bbox"][0] < s["page_bbox"][2] and not (s["text"].endswith(" ") or t["text"].startswith(" ")):
                joined = any(norm(s["text"]) + norm(t["text"]) in norm(r["text"]) or norm(t["text"]) + norm(s["text"]) in norm(r["text"]) for r in recs)
                if not joined: touch.append((s["text"], t["text"], round(t["page_bbox"][0] - s["page_bbox"][2], 2)))
    real = [g for g in gaps if not any(norm(g[1]) + norm(g[2]) in norm(r["text"]) for r in recs)]; mingap = min(real) if real else None
    CH = []
    for b_ in pg.get_text("rawdict")["blocks"]:
        if b_["type"] != 0: continue
        for l_ in b_["lines"]:
            for s_ in l_["spans"]:
                for ch in s_["chars"]:
                    if ch["c"].strip(): CH.append(dict(c=ch["c"], x0=ch["bbox"][0] + bx[0], x1=ch["bbox"][2] + bx[0], y=ch["origin"][1] + bx[1], size=s_["size"], font=s_["font"]))
    pdoc.close(); rows_ = collections.defaultdict(list)
    for ch in CH: rows_[round(ch["y"], 1)].append(ch)
    ys_ = sorted(rows_); vgaps = []
    for i_, y1 in enumerate(ys_):
        for y2 in ys_[i_ + 1:]:
            if y2 - y1 > 2.2 * max(c["size"] for c in rows_[y1] + rows_[y2]): break
            for c1 in rows_[y1]:
                bot1 = c1["y"] - glyph_box(c1["c"], c1["font"])[0] * c1["size"]
                for c2 in rows_[y2]:
                    if c2["x0"] < c1["x1"] - 0.3 and c1["x0"] < c2["x1"] - 0.3: vgaps.append((round(c2["y"] - glyph_box(c2["c"], c2["font"])[1] * c2["size"] - bot1, 2), f"{c1['c']} of the line at {y1}", f"{c2['c']} of the line at {y2}"))
    vmin = min(vgaps) if vgaps else None
    CK.check(f"panel {x}: every string inside its panel page (nothing clipped; {len(sp)} spans, PyMuPDF on the flat page)", not outside, f"outside {outside[:4]}")
    CK.check(f"panel {x}: no string touches a horizontal neighbour on its baseline (min gap {mingap[0] if mingap else None} pt between {mingap[1:] if mingap else None}, gate 2.5) and no stacked glyphs overlap in ink (min vertical clearance {vmin[0] if vmin else None} pt between {vmin[1:] if vmin else None}, glyph outlines, gate above 0)", not touch and (mingap is None or mingap[0] >= 2.5) and (vmin is None or vmin[0] > 0.0), f"touching {touch[:4]}")
    if cells is not None:
        bad = []
        for d in cells:
            r = d["rect"]; hits = [s for s in sp if norm(s["text"]) == d["text"] and r[1] - 1 <= s["origin"][1] <= r[3] + 1 and r[0] - 30 <= s["origin"][0] <= r[2] + 30]
            if not hits or any(s["page_bbox"][0] < r[0] + 0.5 or s["page_bbox"][2] > r[2] - 0.5 or s["page_bbox"][1] < r[1] + 0.5 or s["page_bbox"][3] > r[3] - 0.5 for s in hits): bad.append((d["text"], [round(v, 1) for v in r], [[round(v, 1) for v in s["page_bbox"]] for s in hits]))
        CK.check(f"panel {x}: every cell string sits inside its own cell with at least 0.5 pt to the cell edge ({len(cells)} cells)", not bad, f"{bad[:3]}")
    REC[f"gaps_{x}"] = dict(min_horizontal=mingap, min_vertical_ink=vmin, n_spans=len(sp)); return sp
VIS = {}
for x, D in (("a", DA), ("b", DB)): VIS[x] = flat_page_checks(x, f"{WORK}/panel_{x}_m.pdf", S[x]["box"], D["drawn_records"])
# the whole sheet: no two txtwrite line boxes intersect (lines joined per baseline, size and region as in the census; the integer boxes shrunk 0.5 pt)
def boxes_overlap(lines_, shrink=0.5):
    bad = []
    LB = [(l["x0"] + shrink, l["y"] - 0.716 * l["size"] + shrink, l["x1"] - shrink, l["y"] + 0.212 * l["size"] - shrink, l["text"]) for l in lines_]
    for i in range(len(LB)):
        for j in range(i + 1, len(LB)):
            a, b = LB[i], LB[j]
            if a[0] < b[2] and b[0] < a[2] and a[1] < b[3] and b[1] < a[3]: bad.append((a[4], b[4]))
    return bad
ov = boxes_overlap(nl)
CK.check(f"whole sheet: no two text line boxes intersect (gs txtwrite line boxes, cap height to descender, {len(nl)} lines, band included)", not ov, f"{ov[:6]}")
bch = []
for s_ in nsp:
    if s_["y"] < BAND_TOP: continue
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
                if c2["x0"] < c1["x1"] - 0.3 and c1["x0"] < c2["x1"] - 0.3: bgaps.append((round(c2["y"] - glyph_box(c2["c"], c2["font"])[1] * c2["size"] - bot1, 2), f"{c1['c']} at y {y1}", f"{c2['c']} at y {y2}"))
bmin = min(bgaps) if bgaps else None
CK.check(f"band (panel c, the coordinator's band_4e_r54.pdf): no stacked glyphs overlap in ink on the composed sheet (min vertical clearance {bmin[0] if bmin else None} pt between {bmin[1:] if bmin else None}; txtwrite integer boxes, 0.5 pt)", bmin is None or bmin[0] > 0.0, f"{sorted(bgaps)[:4]}")
REC["gaps_band"] = dict(min_vertical_ink=bmin)

# ---------------------------------------------------------------- (8) placement rules
for x, D, G, labs in (("a", DA, GEO["a"], {ln for r in S["a"]["rows"] for ln in r["lines"]} | {"Negative controls"}), ("b", DB, GEO["b"], {ln for r in S["b"]["rows"] for ln in r["lines"]} | {"Negative controls"})):
    sp = VIS[x]; ax0, ax1 = (G["axis_x"] if x == "a" else (G["plot_box_page"][0], G["plot_box_page"][2])); xvf = xv_a if x == "a" else xv_b
    lab_right = max(s["page_bbox"][2] for s in sp if norm(s["text"]) in labs); tip = max(xvf(d["hi"]) for d in D["values"]) + 0.85
    hr_sp = [s for s in sp if re.fullmatch(r"\d+\.\d\d \(\d+\.\d\d-\d+\.\d\d\)", norm(s["text"]))]; q_sp = [s for s in sp if norm(s["text"]).startswith("<0.") or re.fullmatch(r"0\.\d\d\d", norm(s["text"]))]
    hq = [min(t["page_bbox"][0] for t in q_sp if abs(t["origin"][1] - s["origin"][1]) < 0.5) - s["page_bbox"][2] for s in hr_sp if any(abs(t["origin"][1] - s["origin"][1]) < 0.5 for t in q_sp)]
    col_r, qcol_r = (G["col_right"], G["qcol_right"]) if x == "a" else (G["col_right_x"], G["qcol_right_x"])
    CK.check(f"panel {x}: row labels end at least 4 pt before the plot, the HR column starts at least 2.5 pt after the right-most whisker cap, HR to q at least 2.5 pt on every sub-row ({len(hr_sp)} HR and {len(q_sp)} q strings, right-aligned on their column edges within 0.5 pt)", ax0 - lab_right >= 4.0 and min(s["page_bbox"][0] for s in hr_sp) - tip >= 2.5 and min(hq) >= 2.5 and len(hr_sp) == len(q_sp) == len(D["values"]) and max(abs(s["page_bbox"][2] - col_r) for s in hr_sp) < 0.5 and max(abs(s["page_bbox"][2] - qcol_r) for s in q_sp) < 0.5, f"label to plot {ax0 - lab_right:.2f} pt, whisker cap {tip:.2f} to HR column {min(s['page_bbox'][0] for s in hr_sp) - tip:.2f} pt, HR to q min {min(hq):.2f} pt")
CK.check("the two axis rules on one line (panel a's pitch fills panel b's block, the round-30 rule) and every sub-row pitch at least 14.5 pt (b 15.0, a derived)", abs(GEO["a"]["rule_y"] - GEO["b"]["plot_box_page"][3]) < 1e-6 and S["a"]["sub_pitch"] >= 14.5 and S["b"]["sub_pitch"] >= 14.5, f"rules {GEO['a']['rule_y']:.3f} / {GEO['b']['plot_box_page'][3]:.3f}, pitches a {S['a']['sub_pitch']:.3f} b {S['b']['sub_pitch']:.3f}")
def ink_bottom(sp): return max(s["origin"][1] + (0.212 * s["size"] if any(ch in "gjpqy(),;[]{}/" for ch in s["text"]) else 0.0) for s in sp)
ab_ink = max(ink_bottom(VIS["a"]), ink_bottom(VIS["b"])); c_top = S["letters"]["c"]["base"] - 0.905 * 18.0
CK.check("letter c (18 pt) clears the lowest ink of panels a and b (its ascender box top at least 6 pt above the ink, the round-49 gate) and the band top sits below the letter's baseline, page bottom = band top + band height", c_top - ab_ink >= 6.0 and BAND_TOP > S["letters"]["c"]["base"] and abs(BAND_TOP + bp[0]["band_page"][1] - hn) < 0.01, f"c box top {c_top:.2f}, a/b ink bottom {ab_ink:.2f}, clearance {c_top - ab_ink:.2f} pt, band top {BAND_TOP:.3f}")
leg = [s for s in VIS["b"] if norm(s["text"]).startswith("No or mild apnea")]; CK.check("panel b legend: the first line ends inside the page (1.5 pt margin), the three lines above the plot, the key of panel a beside no data", leg and leg[0]["page_bbox"][2] <= S["page_w"] - 1.5 and max(s["page_bbox"][3] for s in VIS["b"] if s["origin"][1] <= S["b"]["head_base"] + 0.1) < S["b"]["y0"] - 0.716 * 14.0, f"legend ends {leg[0]['page_bbox'][2]:.2f}" if leg else "no legend")

# ---------------------------------------------------------------- (9) the 150 dpi render, (10) the flat twin's census
def render(pdf, png, w_pt, h_pt, name):
    if os.path.exists(png): os.remove(png)
    t0 = time.time()
    with L.render_lock(): r = subprocess.run([LM.GS, "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=png16m", "-r150", f"-sOutputFile={png}", pdf], capture_output=True, text=True, timeout=600)
    from PIL import Image
    im = Image.open(png) if os.path.exists(png) else None; ew, eh = round(w_pt * 150 / 72), round(h_pt * 150 / 72)
    CK.check(f"{name}: 150 dpi Ghostscript render written fresh (png16m, render lock), pixel size = page size at 150 dpi (1 px)", r.returncode == 0 and im is not None and abs(im.width - ew) <= 1 and abs(im.height - eh) <= 1 and "error" not in (r.stdout + r.stderr).lower(), f"{png} {im.size if im else None} in {time.time() - t0:.1f} s (expected {ew} x {eh})")
    return dict(path=png, size=list(im.size) if im else None, sha256=L.sha256(png) if im else None)
REC["png150"] = render(OUT, PNG150, wn, hn, "Main_Fig4")
fl = CS.lines(CS.spans(CS.txtwrite(FLAT, f"{VER}/census/flat_txtwrite.xml")), reg_new)
CK.check("the flat twin's txtwrite census carries the same line multiset (text and size) as the nested compose", collections.Counter((l["text"], l["size"]) for l in fl) == collections.Counter((l["text"], l["size"]) for l in nl), f"{len(fl)} vs {len(nl)} lines")
REC["census"] = f"{VER}/census/census.json"; REC["checks"] = f"{VER}/checks.txt"; json.dump(REC, open(f"{VER}/verify_record.json", "w"), indent=1, ensure_ascii=False)
ok_main = CK.write(f"{VER}/checks.txt")

# ================================================================ the counts panel (V30's panel c), standalone
CK2 = L.Checks(f"panel_counts_r54 round 54 (lane L4, {time.strftime('%Y-%m-%d %H:%M')}): V30's Main_Fig4 panel c (the apnea-by-T90 count grid) as a standalone page at the round-54 sizes, no letter, no title. NEW = {OUT_C} (sha {L.sha256(OUT_C)[:16]}); V30's strings from {R49V}/{SHEET}_c_drawn.json.")
REC2 = {"out": OUT_C, "out_sha256": L.sha256(OUT_C), "written": LM.now()}
CK2.check("compose log: the page on disk is the composed one, from this chain's panel page (sha256 = build manifest)", CLC["output_sha256"] == REC2["out_sha256"] and CLC["panel_sha256"] == BM["panels"]["c"]["sha256"] == L.sha256(f"{WORK}/panel_c_m.pdf"), CLC["output_sha256"][:16])
nc, wc, hc = mediabox(OUT_C); ngc, wgc, hgc = LM.gs_pdfinfo(OUT_C); nfc, wfc, hfc = LM.gs_pdfinfo(FLAT_C)
CK2.check(f"page box declared: {wc:.3f} x {hc:.3f} pt, one page, = the spec's, pypdf and gs agree, the flat twin the same box", nc == ngc == nfc == 1 and abs(wc - S["c"]["box"][2]) < 0.01 and abs(hc - S["c"]["box"][3]) < 0.01 and abs(wgc - wc) < 0.01 and abs(hgc - hc) < 0.01 and abs(wfc - wc) <= 0.05 and abs(hfc - hc) <= 0.05, f"flat {wfc} x {hfc}")
csp = CS.spans(CS.txtwrite(OUT_C, f"{PCV}/census/new_txtwrite.xml")); cl = CS.lines(csp)
c_old = collections.Counter(norm(r["text"]) for r in RC49["drawn_records"]); c_new = collections.Counter(l["text"] for l in cl)
CK2.check(f"census: the same strings as V30's panel c ({sum(c_old.values())} strings), no letter, no title", c_old == c_new, f"lost {sorted((c_old - c_new).elements())[:6]} gained {sorted((c_new - c_old).elements())[:6]}")
bad_cs = [(l["text"], l["size"]) for l in cl if abs(l["size"] - (15.0 if l["text"] == "T90, % of the recording" else 14.0)) > 0.01]
CK2.check("census: every string at its round-54 size (the title 15 pt, band labels, row labels, cells and the key 14 pt), nothing below 12", not bad_cs and min(l["size"] for l in cl) >= 14.0 - 1e-6, f"wrong {bad_cs[:5]}; sizes {sorted(set(l['size'] for l in cl))}")
CK2.check("fonts: Arial faces only", all("Arial" in f for f in set(s["font"] for s in csp)), f"{sorted(set(s['font'] for s in csp))}")
XT = json.load(open(L.hydrated(f"{L.NUM}/crosstab_v2.json"))); V.sidecar_v8(f"{L.NUM}/crosstab_v2.json")
bad_c = [(d["row"], d["col"]) for d in DC["values"] if d["text"] != f"{int(XT['rows3'][[r['category'] for r in XT['rows3']].index(d['row'])]['n_' + {'≤1%': '0-1%', '>1–5%': '1-5%', '>5–10%': '5-10%', '>10%': '>10%'}[d['col']]]):,} ({L.pct1(d['n'], d['group_n'])})"]
CK2.check("all 12 cells re-derived from crosstab_v2.json rows3 (count, percent half up 1 dp) equal the drawn record, and equal V30's cell strings", not bad_c and len(DC["values"]) == 12 and sorted(d["text"] for d in DC["values"]) == sorted(d["text"] for d in RC49["values"]), f"{bad_c[:4]}")
CK2.check("cell fills: the orange ladder rung of every cell as V30's (the same rung labels), text white on the darkest rung", [(d["rung"], d["text_color"]) for d in DC["values"]] == [(d["rung"], d["text_color"]) for d in RC49["values"]], "")
ROWSTR = rowstr(csp)
miss = [(r["text"], r["x"], r["baseline"]) for r in DC["drawn_records"] if not present(r)]
CK2.check(f"every drawn record string ({len(DC['drawn_records'])}) found on the page's text layer at its origin (gs txtwrite, 1.6 pt)", not miss, f"missing {miss[:6]}")
CK_main = CK; CK = CK2
VIS["c"] = flat_page_checks("c", f"{WORK}/panel_c_m.pdf", S["c"]["box"], DC["drawn_records"], cells=DC["values"])
ovc = boxes_overlap(cl); CK.check(f"no two text line boxes intersect on the page ({len(cl)} lines)", not ovc, f"{ovc[:4]}")
REC2["png150"] = render(OUT_C, PNG_C, wc, hc, "panel_counts_r54")
flc = CS.lines(CS.spans(CS.txtwrite(FLAT_C, f"{PCV}/census/flat_txtwrite.xml")))
CK.check("the flat twin's txtwrite census carries the same line multiset as the nested page", collections.Counter((l["text"], l["size"]) for l in flc) == collections.Counter((l["text"], l["size"]) for l in cl), f"{len(flc)} vs {len(cl)} lines")
json.dump(dict(new=OUT_C, strings_old=sorted(c_old.elements()), strings_new=sorted(c_new.elements()), sizes=sorted(collections.Counter(l["size"] for l in cl).items())), open(f"{PCV}/census/census.json", "w"), indent=1, ensure_ascii=False)
json.dump(REC2, open(f"{PCV}/verify_record.json", "w"), indent=1, ensure_ascii=False)
ok_c = CK.write(f"{PCV}/checks.txt")
sys.exit(0 if ok_main and ok_c else 1)
