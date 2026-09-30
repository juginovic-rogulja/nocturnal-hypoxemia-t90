#!$T90_PY
"""Round 54, lane L5: verify Main_Fig5 (the vector compose at 14/15/18 pt, re-laid). Proves, per R54_BRIEF.md and the coordinator's two
additions (height aim 1110, the pdfwrite twin for the census):
  1. page box declared: pypdf mediabox of the vector (and of the flat twin) = 968.94 x (band top + band height), the width = V30's, the height
     at most 1110 (aim) and 1260 (hard gate), the band top = the flat content bottom + 20;
  2. Ghostscript text census (work/census_r54.json, 07_census_r54.py): the new sheet (through the flat twin) against V30's Main_Fig5, whose
     text is the round-49 L5 vector's (the assembly chain V27 = V28 = V29 = V30 = the ROUND49 LRASTER raster of that vector): the same
     multiset of tokens, every token's size mapped by the round-54 size map (no exception), every string present at its new size, no text
     below 12 pt (the band's 14-pt caption censuses at 13.9965 through the band's page scale: declared), the pieces and the vector's own
     census agreeing with the twin;
  3. printed values re-derived from the numbers files: every string of panels b and c with a numbers source (labels, HR (95% CI), P, q) equals
     the value re-derived from its file by the printing rule and sits on the flat page at its position and size; the seven panel a counts
     re-derived from treatment_v2.json equal the generator's strings and are on the sheet at 14 pt;
  4. marks preserved: panel a's drawings identical to the Aug-25 panel (fills remapped) and to the bump-0 rebuild, placed by the same V13
     similarity, and a masked raster diff of the panel a region against V30; panels b and c: the flat build's read-back (dots at ln(HR) of
     the file values on the declared axes, whiskers from ln(lo) to ln(hi), colours and fills as V13) rechecked here from the flat page;
     the band = LSKETCH's band_5d_r54.pdf (sha256 of BANDS_READY.txt) at scale 1.0;
  5. no overlapping text boxes on the composed sheet (census boxes, integer bboxes shrunk by 0.5 pt) and, on the flat page, the build's
     exact-width check; text clear of marks (the flat build's check, the generator's overlap gate for panel a);
  6. the 150 dpi render (Ghostscript) of this vector, the deliverable copies (Main_Fig5.pdf = the vector, Main_Fig5_150dpi.png), the flat twin;
  7. wording: no em dash, semicolon, the word 'printed' or a banned word in any string of the sheet.
Crops at 300 dpi for the eye in verify/crops/. usage: 07b_verify_r54.py"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import collections, gc, json, os, re, shutil, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND54_2026-09-28/lanes/L5/Main_Fig5/scripts")
import l5_lib as L
import fitz, numpy as np
from scipy import ndimage
from PIL import Image
from pypdf import PdfReader
import l5a14 as A
from v14lib import resolve_source, RULES
RULES["label_line"] = lambda v: str(v)
SHEET = "Main_Fig5"; SD = f"{L.LANE}/{SHEET}"; WORK = f"{SD}/work"; VER = f"{SD}/verify"; os.makedirs(f"{VER}/crops", exist_ok=True)
V30PDF = f"{L.V30}/{SHEET}.pdf"; VEC = f"{WORK}/{SHEET}_r54_vector.pdf"; TWIN = f"{WORK}/{SHEET}_flat.pdf"; PNG150 = f"{SD}/{SHEET}_150dpi.png"; DELIV = f"{SD}/{SHEET}.pdf"
clog = json.load(open(f"{WORK}/compose_log.json")); REC = json.load(open(f"{WORK}/flat_build_record.json")); DR = json.load(open(f"{WORK}/{SHEET}_bc_drawn.json")); CEN = json.load(open(f"{WORK}/census_r54.json"))
MAP = clog["panel_a"]["map"]; W, H = clog["page"]; BAND_TOP = clog["band_top"]; PT = REC["sizes"]
assert clog["output_sha256"] == L.sha256(VEC), "compose log does not match the vector"
assert CEN["new_sources"]["vector_sha256"] == clog["output_sha256"], "census is not of this vector"
# deliverables: the vector copy and the 150 dpi render
shutil.copyfile(VEC, DELIV)
if not (os.path.exists(PNG150) and os.path.exists(PNG150[:-4] + "_render_record.json") and json.load(open(PNG150[:-4] + "_render_record.json"))["pdf_sha256"] == L.sha256(VEC)):
    import subprocess; subprocess.run([L.PY, f"{SD}/scripts/06b_render_r54.py", VEC, PNG150, "150", "render150_deliv"], check=True)
rr = json.load(open(PNG150[:-4] + "_render_record.json")); assert rr["pdf_sha256"] == L.sha256(VEC) and rr["png_sha256"] == L.sha256(PNG150), "the 150 dpi render is not of this vector"
ck = L.Checks(f"Lane L5 round 54 (labels, values, ticks, keys 14 pt; column heads, axis titles, group head 15; letters and title 18; band notes 15, labels and captions 14), {SHEET}, {L.now()}. OLD = V30 {V30PDF} (sha256 {L.sha256(V30PDF)[:16]}; the assembly manifests V27 = V28 = V29 = V30 trace it to the ROUND49 LRASTER 600-dpi raster of the round-49 L5 vector {L.R49_VECTOR}, sha256 {L.sha256(L.R49_VECTOR)[:16]}, whose text is the old side of the census). NEW = {VEC} (sha256 {L.sha256(VEC)[:16]}; delivered as {DELIV}): panel a rebuilt from its generator at 14/15 pt under the V13 similarity, panels b and c re-laid by the round-49/42 builder at the round-54 sizes (axes compressed to keep every clearance, columns re-placed, key and reference line re-placed), the round-54 band, letters and title 18 pt. Page {W} x {H}.")


def spans_flat(pdf):
    d = fitz.open(L.hydrated(pdf)); out = []
    for b in d[0].get_text("rawdict")["blocks"]:
        if b["type"] != 0: continue
        for l in b["lines"]:
            for s in l["spans"]:
                t = L.norm_text("".join(c["c"] for c in s["chars"]))
                if not t.strip(): continue
                o = s["chars"][0]["origin"]; bb = s["bbox"]
                out.append(dict(text=t.strip(), x0=bb[0], ox=o[0], y=o[1], x1=bb[2], y0=bb[1], y1=bb[3], size=round(s["size"], 2), font=s["font"]))
    d.close(); gc.collect(); return out


def resolve(sf, sk):
    if sf.endswith(".json"):
        obj = json.load(open(L.hydrated(sf))); parts = [p for p in sk.split("/") if p != ""]
        for p in parts[:-1]: obj = obj[p]
        leaf = parts[-1]
        if "," in leaf: return [obj[f] for f in leaf.split(",")]
        return obj[leaf]
    return resolve_source(sf, sk)


# 1 page box
r_new = PdfReader(VEC); r_old = PdfReader(V30PDF); r_twin = PdfReader(TWIN); bn = r_new.pages[0].mediabox; bo = r_old.pages[0].mediabox; bt = r_twin.pages[0].mediabox
wn, hn, wo, ho, wt, ht = float(bn.width), float(bn.height), float(bo.width), float(bo.height), float(bt.width), float(bt.height)
bh = clog["band_page"][1]; cb = REC["content_bottom"]
ck.log(f"page box declared: the vector {wn:.4f} x {hn:.4f} pt (pypdf mediabox, one page) = the flat twin {wt:.4f} x {ht:.4f}; width = V30's {wo:.4f}; height = band top {BAND_TOP} + band height {bh:.3f} (V30 {ho:.3f}: +{hn - ho:.2f} pt), at most the coordinator's aim {L.HEIGHT_AIM} and the hard gate {L.HEIGHT_HARD}; band top = the flat content bottom {cb:.2f} (panel b's second axis title line, descender included) + 20 pt rounded up to 0.1; letter d at ({clog['letter_d']['x']}, {clog['letter_d']['baseline']}) = the V16 offset from the band top, letter a at ({clog['letter_a']['x']}, {clog['letter_a']['baseline']}) = the V13 origin, title at {clog['title']['origin']} (baseline 14 -> 18 for the 18-pt title)",
       len(r_new.pages) == 1 and len(r_twin.pages) == 1 and abs(wn - wo) < 0.05 and abs(wn - wt) < 0.05 and abs(hn - ht) < 0.05 and abs(hn - (BAND_TOP + bh)) < 0.01 and hn <= L.HEIGHT_AIM and hn <= L.HEIGHT_HARD and BAND_TOP >= cb + 20 - 1e-6 and BAND_TOP < cb + 20.2 and abs(clog["letter_a"]["baseline"] - 41.76) < 0.01)
# 2 census
c0 = CEN["old_r49_census_vs_old_twin"]; c1 = CEN["strings_old_vs_new_twin"]; c2 = CEN["strings_new_pieces_vs_new_twin"]; s2 = CEN["sizes_old_vs_new_twin"]; s3 = CEN["sizes_old_vs_new_pieces"]; nv = CEN["new_vector_census"]
ck.log(f"Ghostscript text census (txtwrite -dTextFormat=0 through the pdfwrite twins): OLD = the round-49 lane's census of the round-49 vector = a fresh census of its pdfwrite twin (tokens and sizes equal: {c0['equal'] and c0['sizes_equal']}); NEW twin tokens = OLD tokens (multiset of {sum(collections.Counter(tok for t, *_ in CEN['strings']['new_twin'] for tok in t.split()).values())} tokens in {CEN['counts']['new_twin']} strings, lost {c1['lost']} gained {c1['gained']}); the new pieces (panel a page x S, flat page, band, stamps) = the twin (lost {c2['lost']} gained {c2['gained']}); the vector's own txtwrite (finished under the watchdog: {nv['finished']}) = the twin ({nv.get('equal_to_twin')}, sizes equal {nv.get('sizes_equal')})",
       c0["equal"] and c0["sizes_equal"] and c1["equal"] and c2["equal"] and nv.get("finished") and nv.get("equal_to_twin") and nv.get("sizes_equal"))
ck.log(f"census sizes: every token's old size maps to its new size by the round-54 map {s2['size_map']} ({s2['n_tokens_mapped']} tokens, exceptions {s2['exceptions']}; pieces: {s3['n_tokens_mapped']} tokens, exceptions {s3['exceptions']}); size histogram old {CEN['size_histogram']['old']} -> new {CEN['size_histogram']['new_twin']}; the smallest text on the sheet {CEN['min_size_new']} pt (the band's 14-pt caption 'treat sleep apnea' through the band's page scale 0.99975 = 13.9965, declared; everything else 14, 15 or 18)",
       not s2["exceptions"] and not s3["exceptions"] and s2["n_tokens_mapped"] >= 250 and CEN["min_size_new"] >= 13.99)
j = CEN["strings_old_vs_new_twin_joined"]
ck.info(f"joined strings differ only by how txtwrite joins neighbours at the new sizes (the V13 heads 'HR (95% CI)' and 'P' were one joined string, the 14-pt tick labels join into one): lost {j['lost']} gained {j['gained']}; tokens are compared")
per = CEN["per_string"]; old_only = [r["text"] for r in per if r["old"] and not r["new"]]; new_only = [r["text"] for r in per if r["new"] and not r["old"]]
ck.log(f"every string present at its new size: of the {len(per)} distinct joined strings, those in the old census only {old_only} and in the new only {new_only} are the joining artefacts above (their tokens are present); every new string is at 14, 15 or 18 pt (the band caption 13.9965)", all(any(abs(sz - t) <= 0.03 for t in (14.0, 15.0, 18.0)) for r in per for sz, _, _ in r["new"]) and set(old_only) <= {"HR (95% CI) P", "0.5", "1", "2", "4", "8", "1.5"} and set(new_only) <= {"P", "1 1.5 2", "0.5 1 2 4 8", "HR (95% CI)"})
ck.info(f"multiplicity of spans in the twin census (the flat page is placed twice, clips b and c, and txtwrite emits the clipped-away copy; identical spans are counted once): {CEN['span_multiplicity']['new_twin']}")
# 3 data: flat build checks, re-derivation
fb = open(f"{VER}/checks_flat_build.txt").read().strip().splitlines()
ck.log(f"the flat build's own checks (checks_flat_build.txt, {sum(l.startswith('PASS') for l in fb)} pass, {sum(l.startswith('FAIL') for l in fb)} fail): rows recomputed from the refit file, every drawn value inside the V13 axis value ranges, the axes solved for every clearance at 14 pt, markers and whiskers read back at the file values, P bold and dots filled exactly where P < 0.05, sizes by role, no overlapping text, text clear of marks, RESULT ALL PASS", fb[-1] == "RESULT ALL PASS" and REC["flat_sha256"] == L.sha256(REC["flat"]) and clog["flat_sha256"] == REC["flat_sha256"])
FL = spans_flat(REC["flat"]); bc = [d for d in DR if d.get("panel") in ("b", "c") and d.get("source_file", "").startswith("/")]
bad_derive = []; bad_find = []; n_ok = 0
for d in bc:
    try:
        if d["rule"] in ("label", "label_line") and d["source_file"].endswith(".json"):
            ptr = d["source_key"].split(" (")[0]; node = resolve(d["source_file"], ptr); key = ptr.split("/")[-1]
            ok = isinstance(node, dict) and "adj_pub_hr" in node and ((d["text"] == A.plab(key)) if d["rule"] == "label" else (d["text"] in key and key == d["source_value"])); got = key
        else:
            v = resolve(d["source_file"], d["source_key"]); got = RULES[d["rule"]](v) if d["rule"] != "label" else A.plab(str(v)); ok = d["text"] == got
        if not ok: bad_derive.append((d["text"], d["rule"], d["source_key"][:60], str(got)[:40]))
    except Exception as e: bad_derive.append((d["text"], d["rule"], d["source_key"][:60], f"{type(e).__name__}: {str(e)[:60]}"))
    x = d["x"]; ha = d["ha"]
    hit = [s for s in FL if s["text"] == L.norm_text(d["text"]).strip() and abs(s["y"] - d["baseline"]) <= 0.6 and abs((s["x0"] if ha == "left" else s["x1"] if ha == "right" else (s["x0"] + s["x1"]) / 2) - x) <= 0.8 and abs(s["size"] - d["size"]) < 0.05 and abs(d["size"] - PT[d["role"]]) < 1e-9]
    if not hit: bad_find.append((d["text"], d["panel"], x, d["baseline"]))
    else: n_ok += 1
ck.log(f"every printed string of panels b and c with a numbers source ({len(bc)} drawn records: labels, HR (95% CI), P and q) re-derived from its file (json pointer or csv lookup, the printing rule: hr_ci half up at 2 dp, fmt_p two significant figures, label) equals the drawn text, and is on the flat page at its position and at its round-54 size ({n_ok} found)", not bad_derive and not bad_find, f"derive {bad_derive[:6]} find {bad_find[:6]}")
# 3b the rows' values against the files, independently of the drawn records, and the marker positions against the declared axes
Rj = json.load(open(L.hydrated(A.SPEC["fig5b_refit_v3_json"]))); rows_b = REC["rows_b"]; ab, bb = REC["panel_b_axis"]["a"], REC["panel_b_axis"]["b"]
vals_ok = all(abs(r["hr"] - Rj["outcomes"][r["outcome"]]["adj_B_hr"]) < 1e-12 and abs(r["lo"] - Rj["outcomes"][r["outcome"]]["adj_B_lo"]) < 1e-12 and abs(r["hi"] - Rj["outcomes"][r["outcome"]]["adj_B_hi"]) < 1e-12 and abs(r["p"] - Rj["outcomes"][r["outcome"]]["adj_B_p"]) < 1e-12 for r in rows_b)
pos_ok = all(abs(r["x"] - (ab + bb * np.log(r["hr"]))) < 1e-6 and abs(r["x_lo"] - (ab + bb * np.log(r["lo"]))) < 1e-6 and abs(r["x_hi"] - (ab + bb * np.log(r["hi"]))) < 1e-6 for r in rows_b)
import pandas as pd
cc = pd.read_csv(L.hydrated(A.SPEC["results_continuous"])).set_index("key"); rows_c = REC["rows_c"]; ac, bcx = REC["panel_c_axis"]["a"], REC["panel_c_axis"]["b"]
vals_c = all(abs(r["hr"] - cc.loc[r["key"], "hr_sd"]) < 1e-12 and abs(r["lo"] - cc.loc[r["key"], "lo_sd"]) < 1e-12 and abs(r["hi"] - cc.loc[r["key"], "hi_sd"]) < 1e-12 and abs(r["q"] - cc.loc[r["key"], "q_sd"]) < 1e-12 for r in rows_c)
pos_c = all(abs(r["x"] - (ac + bcx * np.log(r["hr"]))) < 1e-6 and abs(r["x_lo"] - (ac + bcx * np.log(r["lo"]))) < 1e-6 and abs(r["x_hi"] - (ac + bcx * np.log(r["hi"]))) < 1e-6 for r in rows_c)
itemsF = A.drawings(REC["flat"]); MF = A.markers(itemsF); HL = A.hlines(itemsF, min_len=5)
dotsB = sorted([m for m in MF if m["shape"] == "circle" and abs(m["w"] - 2 * 2.85) < 0.06 and m["cx"] > 500 and m["cy"] > 100], key=lambda m: m["cy"]); sqB = sorted([m for m in MF if m["shape"] == "square" and m["cx"] > 500], key=lambda m: m["cy"]); dotsC = sorted([m for m in MF if m["shape"] == "circle" and abs(m["w"] - 5.29 * 0.95) < 0.06 and m["cx"] < 400], key=lambda m: m["cy"])
mk_b = [(m["cx"], m["cy"], m["open"], m["colour"]) for m in dotsB + sqB]; exp_b = [(r["x"], r["y"], not r["filled"], r["colour"]["mark"]) for r in rows_b]
mk_ok = len(mk_b) == 21 and all(any(abs(m[0] - e[0]) < 0.05 and abs(m[1] - e[1]) < 0.05 and m[2] == e[2] and m[3] == e[3] for m in mk_b) for e in exp_b)
mk_c = len(dotsC) == 11 and all(any(abs(m["cx"] - r["x"]) < 0.05 and abs(m["cy"] - r["y"]) < 0.05 and m["colour"] == A.BLUE for m in dotsC) for r in rows_c)
wh_ok = all(any(abs(l["y"] - r["y"]) < 0.05 and abs(l["x0"] - r["x_lo"]) < 0.05 and abs(l["x1"] - r["x_hi"]) < 0.05 and l["stroke"] == r["colour"]["line"] for l in HL) for r in rows_b) and all(any(abs(l["y"] - r["y"]) < 0.05 and abs(l["x0"] - r["x_lo"]) < 0.05 and abs(l["x1"] - r["x_hi"]) < 0.05 and l["stroke"] == A.BLUE for l in HL) for r in rows_c)
ck.log(f"marks preserved (panels b and c): the 21 b rows' hr, lo, hi and P equal L4_fig5b_refit_v3.json adj_B_* and the 11 c rows' hr_sd, lo_sd, hi_sd, q_sd equal results_continuous.csv (exact); on the flat page (get_drawings) every marker sits at a + b ln(HR) of the declared axes (b: a {ab:.3f}, b {bb:.3f} = {100 * REC['panel_b_axis']['scale_vs_v13']:.1f} percent of V13; c: a {ac:.3f}, b {bcx:.3f} = {100 * REC['panel_c_axis']['scale_vs_v13']:.1f} percent of V13) and every whisker runs from a + b ln(lo) to a + b ln(hi) at the marker's y, orange dots filled exactly where P < 0.05, grey squares for the controls, blue dots in c: {len(mk_b)} b markers, {len(dotsC)} c markers; the axes keep their value ranges and tick values, the whiskers their values",
       vals_ok and pos_ok and vals_c and pos_c and mk_ok and mk_c and wh_ok)
# panel a: the seven counts from treatment_v2.json, the generator's drawn strings, the sheet census
T = json.load(open(L.hydrated(A.SPEC["treatment_v2"]))); fw = T["migration"]["from_worst"]; cvn = T["corrected_vs_not"]
exp_a = [f"n = {fw['n']:,}", f"0–1%, n = {fw['to']['0-1%']:,}", f"1–5%, n = {fw['to']['1-5%']:,}", f"5–10%, n = {fw['to']['5-10%']:,}", f">10%, n = {fw['to']['>10%']:,}", f"n = {cvn['n_corrected']:,} ({round(100 * cvn['n_corrected'] / cvn['n'])}%)", f"n = {cvn['n_not']:,} ({round(100 * cvn['n_not'] / cvn['n'])}%)"]
PAD = clog["panel_a"]["drawn_values"]; gen_strings = [PAD["pre"]["printed"].split("\n")[1]] + [p["printed"] for p in PAD["post"]] + ["n = " + a["printed"].split("n = ")[-1] for a in PAD["arms"]]
sheet_a = [(L.norm_text(s[0]), s[1]) for s in CEN["strings"]["new_twin"] if s[4] < 385 and s[3] < 497]
ck.log(f"panel a: the seven counts re-derived from treatment_v2.json (sha256 {L.sha256(A.SPEC['treatment_v2'])[:16]}) {exp_a} equal the generator's drawn strings and are on the composed sheet's census (panel a region) at 14 pt; the captions at 15 pt; panel a strings on the sheet {sorted(set(sz for _, sz in sheet_a))} pt (V30 12.34 / 13.48)",
       all(L.norm_text(e) in [t for t, _ in sheet_a] for e in exp_a) and all(any(L.norm_text(e) == L.norm_text(g) for g in gen_strings) for e in exp_a) and all(sz == 14.0 for t, sz in sheet_a if any(L.norm_text(e) == t for e in exp_a)) and sorted(set(sz for _, sz in sheet_a)) == [14.0, 15.0, 18.0])
def hx(c): return None if c is None else "#%02x%02x%02x" % tuple(int(round(v * 255)) for v in c[:3])
def draw_keys(pdf, remap=None):
    d = fitz.open(pdf); out = sorted((tuple(round(v, 2) for v in it["rect"]), (remap or {}).get(hx(it.get("fill")), hx(it.get("fill"))), hx(it.get("color")), it.get("width")) for it in d[0].get_drawings()); d.close(); return out
REMAP = MAP["fill_remap"]
ck.log(f"marks preserved (panel a): the copied generator reproduces the Aug-25 panel (15 spans, {MAP['n_drawings']} drawings identical under the fill remap {REMAP}); the round-54 page's 16 drawing items (rects, fills, strokes, widths) identical to the bump-0 page's and to the Aug-25 panel's; placed under the V13 similarity S {MAP['S']:.5f}, T ({MAP['TX']:.3f}, {MAP['TY']:.3f}) fitted on 15 strings with residual {MAP['max_residual']:.3f} pt, clipped to {clog['clips']['a']} (V13 {clog['clips']['v13']['a']}: the a/c boundary 393 -> 385 for the 18-pt letter c); the arm text column moved from X_ARM 77.5 to {clog['panel_a']['x_arm']} (text only) so the 14-pt arm strings end at {REC['panel_a_r54_page']['right']:.1f}, 20 pt left of panel b's labels",
       MAP["generator_reproduces_aug25"] and MAP["max_residual"] < 0.05 and draw_keys(f"{WORK}/panel_a/Figure5A_r54.pdf") == draw_keys(f"{WORK}/panel_a/Figure5A_bump0.pdf") == draw_keys(MAP["old_pdf"], REMAP) and abs(clog["panel_a"]["sheet_sizes"]["blocks_arms"] - 14.0) < 1e-3 and abs(clog["panel_a"]["sheet_sizes"]["captions"] - 15.0) < 1e-3)
# band
ready = open(L.BANDS_READY).read() if os.path.exists(L.BANDS_READY) else ""
m = re.search(r"band_5d_r54\.pdf\s+([\d.]+) x ([\d.]+) pt\s+sha256=([0-9a-f]{64})", ready)
ck.log(f"band d = LSKETCH's band_5d_r54.pdf (BANDS_READY.txt: {m.group(1) if m else '?'} x {m.group(2) if m else '?'} pt, sha256 {m.group(3)[:16] if m else '?'}) = the composed band (sha256 {clog['band_sha256'][:16]}), page box read from the file {clog['band_page']}, placed at scale 1.0 full width at {clog['band_placed_rect']}; the band's strings at 15 (notes) and 14 (caption) pt", bool(m) and m.group(3) == clog["band_sha256"] and clog["band_is_r54"] and abs(float(m.group(2)) - clog["band_page"][1]) < 0.01 and clog["band_scale"] == 1.0)
# 5 overlaps on the composed sheet (census boxes)
strs = CEN["strings"]["new_twin"]
def cbox(s):
    t, sz, f, x0, y, x1 = s; return (x0 + 0.5, y - 0.716 * sz + 0.5, x1 - 0.5, y + 0.212 * sz - 0.5)
def inter(a, b): return a[0] < b[2] and b[0] < a[2] and a[1] < b[3] and b[1] < a[3]
ov = [(strs[i][0], strs[j][0]) for i in range(len(strs)) for j in range(i + 1, len(strs)) if inter(cbox(strs[i]), cbox(strs[j]))]
ck.log(f"no overlapping text boxes on the composed sheet: {len(strs)} census strings (integer txtwrite boxes shrunk by 0.5 pt, cap top to descender bottom), overlapping pairs {ov}; on the flat page the build checked the {REC['n_drawn']} drawn strings at their exact widths (no pair overlaps) and every string against every whisker, marker, key handle, reference line and spine", not ov)
ck.info("text clear of marks: panels b and c by the flat build's check (exact geometry); panel a by the generator's overlap gate (nooverlap.gate: text against text and text against data, passed for Figure5A_r54); the key and the reference line: the line starts 6 pt under the key")
# 6 renders, deliverables
im_new = np.asarray(Image.open(PNG150).convert("RGB")).astype(int)
ck.log(f"150 dpi render (Ghostscript png16m, {rr['seconds']} s): {rr['px'][0]} x {rr['px'][1]} px of this vector (sha256 {rr['pdf_sha256'][:16]}); deliverables {DELIV} (sha256 {L.sha256(DELIV)[:16]} = the vector) and {PNG150}; the flat twin {TWIN} (sha256 {L.sha256(TWIN)[:16]}, gs pdfwrite 1.7 /prepress, no downsampling, colours unchanged, fonts embedded, one page {wt:.2f} x {ht:.2f}) for the coordinator's census", im_new.shape[1] == rr["px"][0] and rr["dpi"] == 150 and L.sha256(DELIV) == L.sha256(VEC) and abs(im_new.shape[1] - W * 150 / 72) < 2 and abs(im_new.shape[0] - H * 150 / 72) < 2)
# 4b panel a region masked raster diff against V30 (the panel a marks did not move: the same similarity, the same drawings)
v30_150 = f"{WORK}/v30_150dpi_gs.png"; L.gs_render(V30PDF, v30_150, 150, f"{SHEET}_gs150_v30", max_s=1800); im_old = np.asarray(Image.open(v30_150).convert("RGB")).astype(int)
s150 = 150 / 72.0; hh = min(im_new.shape[0], im_old.shape[0]); ww = min(im_new.shape[1], im_old.shape[1]); mask = np.zeros((hh, ww), bool)
def box(x0, y0, x1, y1, pad=2.5):
    mask[max(int((y0 - pad) * s150), 0):min(int((y1 + pad) * s150) + 1, hh), max(int((x0 - pad) * s150), 0):min(int((x1 + pad) * s150) + 1, ww)] = True
for t, sz, f, x0, y, *rest in CEN["strings"]["old"] + CEN["strings"]["new_twin"]:
    wdt = A.Sheet(10, 10).width(t, sz, "Bold" in f) if len(t) else 0
    box(x0 - 1, y - 0.95 * sz, x0 + wdt + 2, y + 0.3 * sz)
region = np.zeros_like(mask); region[int(16 * s150):int(385 * s150), :int(497 * s150)] = True
d = (np.abs(im_new[:hh, :ww] - im_old[:hh, :ww]) > 40).any(axis=2); out = d & ~mask & region
g = np.abs(np.diff(im_old[:hh, :ww].sum(axis=2), axis=0, prepend=0)) + np.abs(np.diff(im_old[:hh, :ww].sum(axis=2), axis=1, prepend=0)); edge2 = ndimage.binary_dilation(g > 60, iterations=2)
far = out & ~edge2; lab, n = ndimage.label(far); sizes = sorted((int((lab == i).sum()) for i in range(1, n + 1)), reverse=True)
locs = [tuple(round(float(v) / s150, 1) for v in (o[1].start, o[0].start)) for o in ndimage.find_objects(lab)][:10]
ck.log(f"masked raster diff of the panel a region (x < 497, 16 < y < 385) at 150 dpi against V30 rendered by the same Ghostscript, outside the union of the old and new text boxes: {int(out.sum())} px differ, {int((out & edge2).sum())} within 2 px of an edge of the old image (the 600 dpi raster resampled), {int(far.sum())} farther away in {n} clusters (largest {sizes[0] if sizes else 0} px; isolated 1-2 px flips tolerated, a cluster of 3 or more or more than 20 such pixels fails): the flows, blocks and brackets did not move", (not sizes or sizes[0] <= 2) and int(far.sum()) <= 20, f"far clusters at (pt) {locs}")
Image.fromarray((np.stack([out * 255, mask * 60, np.zeros_like(out)], axis=2)).astype(np.uint8)).save(f"{VER}/raster_diff_panel_a_outside_text_150dpi.png")
ck.info("panels b and c and the band are re-laid (axes compressed, columns and key re-placed, the band top 2 pt lower): their raster difference against V30 is the layout, proven instead by the vector read-back above")
# 7 wording
alltxt = [s[0] for s in CEN["strings"]["new_twin"]]
bad_w = [t for t in alltxt if "—" in t or ";" in t or re.search(r"\bprinted\b", t, re.I) or re.search(r"\b(honest|honestly|straightforward|null|prespecified|pre-specified)\b", t, re.I)]
ck.log(f"wording: no em dash, semicolon, the word 'printed' or a banned word in any of the {len(alltxt)} strings of the sheet", not bad_w, f"{bad_w}")
# crops at 300 dpi for the eye
png300 = f"{WORK}/{SHEET}_300dpi_r54.png"
if not (os.path.exists(png300) and os.path.getmtime(png300) > os.path.getmtime(VEC)): L.gs_render(VEC, png300, 300, f"{SHEET}_gs300", max_s=1800)
im300 = Image.open(png300).convert("RGB"); Z = 300 / 72.0
crops = {"title_letter_a": (0, 0, 300, 60), "panel_a_arms": (360, 100, 500, 300), "panel_a_captions": (20, 300, 500, 385), "b_key_heads": (497, 40, 969, 100), "b_rows_top": (497, 100, 969, 260), "b_two_line_row": (497, 330, 969, 400), "b_controls_axis": (497, 560, 969, 740), "c_heads_rows": (20, 400, 497, 560), "c_rows_axis": (20, 560, 497, 740), "letter_d_band": (0, 730, 969, 948)}
for nm, (x0, y0, x1, y1) in crops.items():
    im300.crop((int(x0 * Z), int(y0 * Z), min(int(x1 * Z), im300.width), min(int(y1 * Z), im300.height))).save(f"{VER}/crops/{nm}_300dpi.png")
ck.info(f"crops at 300 dpi in verify/crops/: {sorted(os.listdir(f'{VER}/crops'))}")
ck.info(f"round-54 layout decisions (declared): panel b labels, key and letter b 9 pt right (506.45); panel b axis at {100 * REC['panel_b_axis']['scale_vs_v13']:.1f} percent of the V13 scale, reference line at x {REC['panel_b_axis']['ref_x']:.2f} (V13 687.35) starting at y {REC['panel_b_axis']['ref_top']:.2f} under the key (V13 50.84), HR column right edge {REC['panel_b_axis']['hr_right']:.2f} (V13 896.98), P column right edge {REC['panel_b_axis']['p_right']:.2f} (V13 933.20); panel c axis at {100 * REC['panel_c_axis']['scale_vs_v13']:.1f} percent, reference line at x {REC['panel_c_axis']['ref_x']:.2f} (V13 164.30), HR column right edge {REC['panel_c_axis']['hr_right']:.2f} (V13 435.63), q column right edge {REC['panel_c_axis']['q_right']:.2f} (V13 471.76); rows at their V13 y (pitches 24.670 and 21.354); the two-line label at the row y minus and plus {REC['labels']['two_line_nudge_pt']} pt; key entries {REC['key_geometry']['dy']} pt apart; markers at the cap centre of the 14-pt label; tick labels at baseline {REC['panel_b_axis']['tick_label_baseline']:.2f} (V13 687.49); axis titles at {REC['panel_c_axis']['title_baseline']:.2f} and {[round(v, 2) for v in REC['panel_b_axis']['title_baselines']]} (V13 702.97, 703.83, 715.93), centred on their axes; column heads right-aligned at their column edges; the band top {BAND_TOP} (V30 747.80); page height {H} (V30 946.23). Nothing kept below its round-54 size.")
ok = ck.write(f"{VER}/checks.txt")
json.dump(dict(vector=VEC, vector_sha256=L.sha256(VEC), deliverable=DELIV, twin=TWIN, twin_sha256=L.sha256(TWIN), png150=PNG150, png150_sha256=L.sha256(PNG150), page=[wn, hn], v30=V30PDF, v30_sha256=L.sha256(V30PDF), r49_vector=L.R49_VECTOR, r49_vector_sha256=L.sha256(L.R49_VECTOR), band=clog["band"], band_sha256=clog["band_sha256"], census=f"{WORK}/census_r54.json", flat=REC["flat"], flat_sha256=REC["flat_sha256"], panel_a=clog["panel_a"]["page"], panel_a_sha256=clog["panel_a"]["sha256"], written=L.now()), open(f"{VER}/verify_record.json", "w"), indent=1)
print("RESULT", "ALL PASS" if ok else "FAIL"); sys.exit(0 if ok else 1)
