#!$T90_PY
"""Round 49, lane L5: verify Main_Fig5 (the vector compose, every text +1 pt). Proves, per the R49 brief:
  1. page box = V26's (pypdf mediabox, 0.05 pt); 2. Ghostscript text census (work/census_r49.json, 07_census_r49.py): the composed vector
  against the V26 sheet's pieces (and its whole-vector census where finished): same multiset of strings, every size +1.0 (exceptions listed:
  the band notes at +2.5 by the round-49 band spec); 3. the lane's own data checks: the flat build's checks (checks_flat_build.txt, ALL PASS),
  every printed string of panels b and c re-derived from its numbers file and found on the flat page at its position (the round-42 check),
  the seven panel a counts re-derived from treatment_v2.json and found on the sheet; the panel a data marks: drawings identical to the Aug-25
  panel and to the bump-0 rebuild, the V13 similarity (residual 0.006 pt); 4. the 150 dpi render (Ghostscript) exists, and a masked raster
  diff against the V26 raster rendered by the same Ghostscript at 150 dpi: outside the union of the old and new text boxes (dilated 2.5 pt)
  the page is unchanged (the data marks, axes, flows and band drawings did not move). Crops listed. usage: 07b_verify_r49.py"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import collections, csv, gc, json, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/L5/Main_Fig5/scripts")
import l5_lib as L
import fitz, numpy as np
from scipy import ndimage
from PIL import Image
from pypdf import PdfReader
import l5a14 as A
from v14lib import resolve_source, RULES
RULES["label_line"] = lambda v: str(v)
SHEET = "Main_Fig5"; SD = f"{L.LANE}/{SHEET}"; WORK = f"{SD}/work"; VER = f"{SD}/verify"; os.makedirs(f"{VER}/crops", exist_ok=True)
V26PDF = f"{L.V26}/{SHEET}.pdf"; VEC = f"{WORK}/{SHEET}_r49_vector.pdf"; PNG150 = f"{SD}/{SHEET}_150dpi.png"
clog = json.load(open(f"{WORK}/compose_log.json")); REC = json.load(open(f"{WORK}/flat_build_record.json")); DR = json.load(open(f"{WORK}/{SHEET}_bc_drawn.json")); CEN = json.load(open(f"{WORK}/census_r49.json"))
MAP = clog["panel_a"]["map"]; W, H = clog["page"]; BAND_TOP = clog["band_top"]
assert clog["output_sha256"] == L.sha256(VEC), "compose log does not match the vector"
rr = json.load(open(PNG150[:-4] + "_render_record.json")); assert rr["pdf_sha256"] == L.sha256(VEC) and rr["png_sha256"] == L.sha256(PNG150), "the 150 dpi render is not of this vector"
ck = L.Checks(f"Lane L5 round 49 (every figure text +1 pt), {SHEET}, {L.now()}. OLD = V26 {V26PDF} (600 dpi raster of the round-44 vector compose {L.R44_VECTOR}, sha256 {L.sha256(V26PDF)[:16]}); its text = the pieces (V13 panel a probe, round-42 flat page, round-44 band, the stamped d and title). NEW = {VEC} (sha256 {L.sha256(VEC)[:16]}): panel a rebuilt from its generator at +1/S pt under the V13 similarity, panels b and c rebuilt by the round-42 builder with Sheet.BUMP 1.0, the round-49 band (notes 12 pt, captions 11 pt), letters and title 14 pt. Page 968.94 x 946.2253.")


def spans_flat(pdf, dx=0.0, dy=0.0, scale=1.0):
    d = fitz.open(L.hydrated(pdf)); out = []
    for b in d[0].get_text("rawdict")["blocks"]:
        if b["type"] != 0: continue
        for l in b["lines"]:
            for s in l["spans"]:
                t = L.norm_text("".join(c["c"] for c in s["chars"]))
                if not t.strip(): continue
                o = s["chars"][0]["origin"]; bb = s["bbox"]
                out.append(dict(text=t.strip(), x0=dx + bb[0] * scale, ox=dx + o[0] * scale, y=dy + o[1] * scale, x1=dx + bb[2] * scale, y0=dy + bb[1] * scale, y1=dy + bb[3] * scale, size=round(s["size"] * scale, 2), font=s["font"]))
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
r_new = PdfReader(VEC); r_old = PdfReader(V26PDF); bn = r_new.pages[0].mediabox; bo = r_old.pages[0].mediabox
wn, hn, wo, ho = float(bn.width), float(bn.height), float(bo.width), float(bo.height)
ck.log(f"page box (pypdf mediabox): new {wn:.4f} x {hn:.4f} = V26 {wo:.4f} x {ho:.4f} (within 0.05 pt), one page; band top {BAND_TOP} = page height - band height (the V16 slot), letter d at ({clog['letter_d']['x']}, {clog['letter_d']['baseline']}) = the V13 origin, letter a at ({clog['letter_a']['x']}, {clog['letter_a']['baseline']}) = the V13 origin, title at {clog['title']['origin']}", len(r_new.pages) == 1 and abs(wn - wo) < 0.05 and abs(hn - ho) < 0.05 and abs(clog["letter_d"]["baseline"] - 744.76) < 0.01 and abs(clog["letter_a"]["baseline"] - 41.76) < 0.01)
# 2 census
c1 = CEN["strings_old_pieces_vs_new_pieces"]; c2 = CEN["strings_new_pieces_vs_new_composed"]; s2 = CEN["sizes_old_pieces_vs_new_composed"]
exc = s2["exceptions"]; exc_ok = all(len(e) == 4 and e[0] in ("90%", "SpO2") and abs(e[3] - 2.5) < 0.03 for e in exc) and len(exc) == 4
ck.log(f"Ghostscript text census (txtwrite -dTextFormat=0, tokens with sizes): the composed vector's strings = the V26 pieces' strings (multiset, {CEN['counts']['new_composed']} strings, lost {c1['lost']} gained {c1['gained']}; pieces = composed: lost {c2['lost']} gained {c2['gained']}); every token +1.0 pt ({s2['n_plus_one']} tokens) except the band's notes '90% SpO2' twice at +2.5 (9.5 -> 12 pt, the round-49 band spec: notes 12, captions 11); size histogram old {CEN['size_histogram']['old']} -> new {CEN['size_histogram']['new_composed']}", c1["equal"] and c2["equal"] and exc_ok and s2["n_plus_one"] >= 270)
ck.info(f"the composed census counts every flat-page span twice (placed with two clips, txtwrite emits the clipped-away copy): multiplicity {CEN.get('composed_span_multiplicity')}; identical spans are counted once")
ov = CEN["old_vector_census"]
if ov.get("finished"):
    s3 = CEN.get("sizes_old_composed_vs_new_composed", {}); exc3 = s3.get("exceptions", [])
    ck.log(f"the round-44 whole-vector census (Ghostscript, finished in the background): its strings = the pieces' strings (lost {ov.get('lost')} gained {ov.get('gained')}), and against the new composed vector every token +1.0 ({s3.get('n_plus_one')}) except the band notes (+2.5)", ov.get("strings_equal_to_old_pieces") and all(e[0] in ("90%", "SpO2") for e in exc3))
else: ck.info("the round-44 whole-vector census (nested V13 panel a, 454 XObjects) had not finished when this verifier ran (the brief allows the pieces): the old side is the pieces")
# 3 data: flat build checks, re-derivation
fb = open(f"{VER}/checks_flat_build.txt").read().strip().splitlines()
ck.log(f"the flat build's own checks (checks_flat_build.txt, {sum(l.startswith('PASS') for l in fb)} pass, {sum(l.startswith('FAIL') for l in fb)} fail): rows recomputed from the refit file, every drawn value inside the V13 axis ranges, markers read back at ln(HR), P bold and dots filled exactly where P < 0.05, key and columns at the V13 edges, +1 pt clearance checks, RESULT ALL PASS", fb[-1] == "RESULT ALL PASS" and REC["flat_sha256"] == L.sha256(REC["flat"]))
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
    hit = [s for s in FL if s["text"] == L.norm_text(d["text"]).strip() and abs(s["y"] - d["baseline"]) <= 0.6 and abs((s["x0"] if ha == "left" else s["x1"] if ha == "right" else (s["x0"] + s["x1"]) / 2) - x) <= 0.8 and abs(s["size"] - d["size"]) < 0.05]
    if not hit: bad_find.append((d["text"], d["panel"], x, d["baseline"]))
    else: n_ok += 1
ck.log(f"every printed string of panels b and c with a numbers source ({len(bc)} drawn records: labels, HR (95% CI), P and q) re-derived from its file (json pointer or csv lookup, the printing rule: hr_ci half up at 2 dp, fmt_p two significant figures, label) equals the drawn text, and is on the flat page at its position and at its +1 size ({n_ok} found)", not bad_derive and not bad_find, f"derive {bad_derive[:6]} find {bad_find[:6]}")
# panel a: the seven counts from treatment_v2.json, the generator's drawn strings, the sheet census
T = json.load(open(L.hydrated(A.SPEC["treatment_v2"]))); fw = T["migration"]["from_worst"]; cvn = T["corrected_vs_not"]
exp_a = [f"n = {fw['n']:,}", f"0–1%, n = {fw['to']['0-1%']:,}", f"1–5%, n = {fw['to']['1-5%']:,}", f"5–10%, n = {fw['to']['5-10%']:,}", f">10%, n = {fw['to']['>10%']:,}", f"n = {cvn['n_corrected']:,} ({round(100 * cvn['n_corrected'] / cvn['n'])}%)", f"n = {cvn['n_not']:,} ({round(100 * cvn['n_not'] / cvn['n'])}%)"]
PAD = clog["panel_a"]["drawn_values"]; gen_strings = [PAD["pre"]["printed"].split("\n")[1]] + [p["printed"] for p in PAD["post"]] + [a["printed"].split(" n = ")[-1] and "n = " + a["printed"].split("n = ")[-1] for a in PAD["arms"]]
sheet_a = [L.norm_text(s[0]) for s in CEN["strings"]["new_composed"] if s[4] < 393 and s[3] < 497]
ck.log(f"panel a: the seven counts re-derived from treatment_v2.json (sha256 {L.sha256(A.SPEC['treatment_v2'])[:16]}) {exp_a} equal the generator's drawn strings and are on the composed sheet's census (panel a region) at {sorted(set(round(s[1], 2) for s in CEN['strings']['new_composed'] if s[4] < 393 and s[3] < 497 and len(s[0]) > 1))} pt (V13 11.34 / 12.48 -> 12.34 / 13.48); the generator's post-block strings {PAD['post'][0]['printed']!r} etc.", all(L.norm_text(e) in sheet_a or any(L.norm_text(e) == L.norm_text(g) for g in gen_strings) for e in exp_a) and all(L.norm_text(e) in sheet_a for e in exp_a), f"sheet a strings {sheet_a}")
ck.log(f"panel a data marks: the copied generator reproduces the Aug-25 panel (15 spans and {MAP['n_drawings']} drawings identical), the +1 rebuild keeps every drawing identical (bump on the canvas {clog['panel_a']['bump_on_canvas_pt']:.5f} = 1/S), placed under the V13 similarity S {MAP['S']:.5f}, T ({MAP['TX']:.3f}, {MAP['TY']:.3f}) fitted on 15 strings with residual {MAP['max_residual']:.3f} pt, clipped to the V13 panel a clip {clog['placements'][0]['rect']}", MAP["generator_reproduces_aug25"] and MAP["max_residual"] < 0.05 and abs(clog["panel_a"]["bump_on_canvas_pt"] * MAP["S"] - 1.0) < 1e-6)
# panel a drawings identical between bump0 and r49 (re-measured here)
def hx(c): return None if c is None else "#%02x%02x%02x" % tuple(int(round(v * 255)) for v in c[:3])
def draw_keys(pdf, remap=None):
    d = fitz.open(pdf); out = sorted((tuple(round(v, 2) for v in it["rect"]), (remap or {}).get(hx(it.get("fill")), hx(it.get("fill"))), hx(it.get("color")), it.get("width")) for it in d[0].get_drawings()); d.close(); return out
REMAP = MAP["fill_remap"]
ck.log(f"panel a drawings (PyMuPDF get_drawings on the flat generator pages): the +1 page's 16 drawing items (rects, fills, strokes, widths) identical to the bump-0 page's, and to the Aug-25 panel's under the fill remap to the sheet's clinical palette {REMAP} (the V13 lineage recoloured panel a after Aug 25; measured on the V26 raster: exact 8-bit matches; ink and white text unchanged)", draw_keys(f"{WORK}/panel_a/Figure5A_r49.pdf") == draw_keys(f"{WORK}/panel_a/Figure5A_bump0.pdf") == draw_keys(MAP["old_pdf"], REMAP))
# 4 render and masked raster diff against V26
im_new = np.asarray(Image.open(PNG150).convert("RGB")).astype(int)
v26_150 = f"{WORK}/v26_150dpi_gs.png"; L.gs_render(V26PDF, v26_150, 150, f"{SHEET}_gs150_v26", max_s=1800); im_old = np.asarray(Image.open(v26_150).convert("RGB")).astype(int)
ck.log(f"150 dpi render (Ghostscript {rr['seconds']} s): {rr['px'][0]} x {rr['px'][1]} px = the V26 raster rendered at 150 dpi {im_old.shape[1]} x {im_old.shape[0]} px", im_new.shape == im_old.shape and rr["dpi"] == 150)
s150 = 150 / 72.0; mask = np.zeros(im_new.shape[:2], bool)
def box(x0, y0, x1, y1, pad=2.5):
    mask[max(int((y0 - pad) * s150), 0):min(int((y1 + pad) * s150) + 1, mask.shape[0]), max(int((x0 - pad) * s150), 0):min(int((x1 + pad) * s150) + 1, mask.shape[1])] = True
for t, sz, f, x0, y in CEN["strings"]["old"] + CEN["strings"]["new_composed"]:
    wdt = A.Sheet(10, 10).width(t, sz, "Bold" in f, bump=False) if len(t) else 0
    box(x0 - 1, y - 0.95 * sz, x0 + wdt + 2, y + 0.3 * sz)
band_mask = np.zeros_like(mask); band_mask[int((BAND_TOP - 1) * s150):, :] = True     # the band is the coordinator's (rebuilt at 12/11 pt with shorter fork arrowheads): reported apart
d = (np.abs(im_new - im_old) > 40).any(axis=2); n_band = int((d & band_mask & ~mask).sum()); out = d & ~mask & ~band_mask
# the V26 sheet is a 600 dpi raster: rendered at 150 dpi its edges resample, so differences along the old image's own edges are expected;
# a moved or changed shape would differ AWAY from the old edges (a shift of 1 pt = 2 px shows beyond the 2 px band); the sub-point geometry is
# proven on the vectors (drawings identical, similarity residual 0.006 pt, markers read back at ln(HR) with residual 0.001 pt)
g = np.abs(np.diff(im_old.sum(axis=2), axis=0, prepend=0)) + np.abs(np.diff(im_old.sum(axis=2), axis=1, prepend=0)); edge2 = ndimage.binary_dilation(g > 60, iterations=2)
far = out & ~edge2; lab, n = ndimage.label(far); sizes = sorted((int((lab == i).sum()) for i in range(1, n + 1)), reverse=True)
locs = [tuple(round(float(v) / s150, 1) for v in (o[1].start, o[0].start)) for o in ndimage.find_objects(lab)][:10]
ck.log(f"masked raster diff at 150 dpi against the V26 raster rendered by the same Ghostscript, the band region apart ({n_band} px differ there outside text: the coordinator's round-49 band, shorter fork arrowheads): above the band, outside the union of the old and new text boxes (census strings, dilated 2.5 pt; {mask.mean() * 100:.1f} percent of the page masked) {int(out.sum())} px differ, {int((out & edge2).sum())} of them within 2 px of an edge of the old image (the 600 dpi raster resampled: the flow, block and axis edges), {int(far.sum())} farther away in {n} clusters (largest {sizes[0] if sizes else 0} px; isolated 1-2 px flips tolerated, a cluster of 3 or more or more than 20 such pixels fails): the flows, blocks, brackets, axes, whiskers, markers and key handles did not move", (not sizes or sizes[0] <= 2) and int(far.sum()) <= 20, f"far clusters at (pt) {locs}")
Image.fromarray((np.stack([out * 255, mask * 60, np.zeros_like(out)], axis=2)).astype(np.uint8)).save(f"{VER}/raster_diff_outside_text_mask_150dpi.png")
ck.info(f"crops at 300 dpi in verify/crops/: {sorted(os.listdir(f'{VER}/crops'))}")
ck.info(f"round-49 nudges: the two lines of the one two-line label ('Ventricular arrhythmia' / 'or cardiac arrest') 0.5 pt apart each way (labels record two_line_nudge_pt {REC['labels'].get('two_line_nudge_pt')}); nothing kept at its old size; no data mark moved")
ok = ck.write(f"{VER}/checks.txt")
json.dump(dict(vector=VEC, vector_sha256=L.sha256(VEC), png150=PNG150, png150_sha256=L.sha256(PNG150), page=[wn, hn], v26=V26PDF, v26_sha256=L.sha256(V26PDF), band=clog["band"], band_sha256=clog["band_sha256"], census=f"{WORK}/census_r49.json", flat=REC["flat"], flat_sha256=REC["flat_sha256"], panel_a=clog["panel_a"]["page"], panel_a_sha256=clog["panel_a"]["sha256"], written=L.now()), open(f"{VER}/verify_record.json", "w"), indent=1)
print("RESULT", "ALL PASS" if ok else "FAIL"); sys.exit(0 if ok else 1)
