#!$T90_PY
"""Round 42, LFIG5b FINAL: verify Main_Fig5, one contrast per condition (Still above 10% on PAP vs the reference; the 1 to 10% line dropped
everywhere), 21 single-line rows at one pitch over the V13 span, the P column, the V13 axis, a two-entry key; the LFIG5b verifier of the
two-series P build adapted. Earlier docstring: verify Main_Fig5, the P-column / V13-axis variant of the LFIG5 sheet (the LFIG5 verifier, verify_fig5_r42_LFIG5_REFERENCE.py,
with: the printed column = the uncorrected P (bold, filled dots where P < 0.05), the head 'P', no q in panel b; panel b's spine, tick labels
and title at the V13 y, measured on the render against panel c's spine; a raster diff against the LFIG5 deliverable confined to panel b). Ghostscript only on
the composed sheets (renders under the lane watchdog); PyMuPDF only on the flat pieces (the flat page, the band). Proves:
  1. page box = the V16 page (band at the V16 slot); 2. the 18-row selection recomputed from the refit file = the drawn rows and order;
  3. every printed string of panels b and c re-derived from its numbers file (json pointer / csv lookup + the printing rule) and found on the
  flat page at its position; bold q faces and filled dots exactly where q < 0.05; the key measured (labels' x, handle start, rows); the 'q'
  head; 4. panel c identical to the round-40 flat page; band text = V16; 5. raster diff at 150 dpi against V16 (same Ghostscript) confined to
  the title strip, panel b and the band: panel a, panel c, the letter d and the strip between pixel-identical; 6. OCR read-back (r38 idiom,
  300 dpi) of every record string on the 600 dpi deliverable, all numeric strings of b and c exact; 7. crops OLD vs NEW; 8. CHANGES csv
  (round-40 sheet -> round 42, dramatic flags where the printed q and the old P sit on different sides of 0.05).
usage: verify_fig5_r42.py [--skip-ocr]"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import collections, csv, gc, json, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND42_2026-09-19/lanes/LFIG5b/Main_Fig5/scripts")
import lfig5b_lib as L, gs_text as T
import fitz, numpy as np
from scipy import ndimage
from PIL import Image
import l5a14 as A                                    # fmt_p, hr_ci, plab, RULES (registered on import), SPEC, drawings, markers, resolve via v14lib
from v14lib import resolve_source, RULES
RULES["label_line"] = lambda v: str(v)
OCR_SCRIPT = f"{L.R38}/figures/scripts/r38_ocr_check.py"
SHEET = "Main_Fig5"; SKIP_OCR = "--skip-ocr" in sys.argv; N = 5
SD = f"{L.LANE}/{SHEET}"; WORK = f"{SD}/work"; VER = f"{SD}/verify"; os.makedirs(f"{VER}/crops", exist_ok=True)
V16PDF = f"{L.V16}/{SHEET}.pdf"; VEC = f"{SD}/{SHEET}_vector.pdf"; RAS = f"{SD}/{SHEET}.pdf"; PREV = f"{WORK}/r42/Main_Fig5.pdf"          # the LFIG5 deliverable (q columns), the "before" side
clog = json.load(open(f"{WORK}/compose_log.json")); rr = json.load(open(f"{WORK}/raster_record.json")); REC = json.load(open(f"{WORK}/flat_build_record.json"))
DR = json.load(open(f"{WORK}/{SHEET}_bc_drawn.json")); W, H = clog["page"]; BAND_TOP = clog["band_top"]; BAND = clog["band"]


def spans_flat(pdf, dx=0.0, dy=0.0, scale=1.0):
    d = fitz.open(L.hydrated(pdf)); out = []
    for b in d[0].get_text("rawdict")["blocks"]:
        if b["type"] != 0: continue
        for l in b["lines"]:
            for s in l["spans"]:
                t = L.norm_text("".join(c["c"] for c in s["chars"]))
                if not t.strip(): continue
                o = s["chars"][0]["origin"]; bb = s["bbox"]
                out.append(dict(text=t.strip(), x0=dx + bb[0] * scale, ox=dx + o[0] * scale, y=dy + o[1] * scale, x1=dx + bb[2] * scale, size=round(s["size"] * scale, 2), font=s["font"]))
    d.close(); gc.collect(); return out


def record_rows(path):
    rows = list(csv.DictReader(open(L.hydrated(path))))
    for r in rows: r["string"] = L.norm_text(r["string"]).strip()
    return rows


def resolve(sf, sk):
    """JSON pointer with an optional comma list of leaf fields (outcomes/<k>/f1,f2,f3), or the v14lib csv lookup."""
    if sf.endswith(".json"):
        obj = json.load(open(L.hydrated(sf))); parts = [p for p in sk.split("/") if p != ""]
        for p in parts[:-1]: obj = obj[p]
        leaf = parts[-1]
        if "," in leaf: return [obj[f] for f in leaf.split(",")]
        return obj[leaf]
    return resolve_source(sf, sk)


def chars(strings): return collections.Counter(ch for s in strings for ch in L.norm_text(s) if not ch.isspace())


def write_rows_csv(path, rows):
    with open(path, "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=["string", "x", "y", "panel"]); wr.writeheader()
        for r in rows: wr.writerow({k: r.get(k, "") for k in wr.fieldnames})


def star_glyphs(png, s, x, y, dpi, size_pt=10.5, digit_pt=10.0):
    from PIL import ImageFont
    Z = dpi / 72.0; font = ImageFont.truetype(L.ARIAL, 40); w = font.getlength(s) / 40 * size_pt * 1.08
    im = Image.open(png).convert("L"); box = (max(0, int((x - 2) * Z)), max(0, int((y - size_pt * 0.95) * Z)), int((x + w + 2) * Z), int((y + size_pt * 0.30) * Z)); a = np.asarray(im.crop(box)).astype(int)
    Hd = 0.716 * digit_pt * Z; best = None
    for name, mask in (("ink", a <= 40), ("white", a >= 245)):
        lab, n = ndimage.label(mask); cs = []
        for i, o in enumerate(ndimage.find_objects(lab)):
            h, wd = o[0].stop - o[0].start, o[1].stop - o[1].start; npx = int((lab[o] == i + 1).sum())
            if npx >= 4 and h < 0.9 * a.shape[0]: cs.append(dict(x0=o[1].start, x1=o[1].stop, y0=o[0].start, h=h, w=wd))
        digits = [k for k in cs if 0.8 * Hd <= k["h"] <= 1.2 * Hd]
        if len(digits) >= 3 and (best is None or len(digits) > len(best[2])): best = (name, cs, digits)
    if best is None: return None, "no digit run"
    name, cs, digits = best; Hh = max(k["h"] for k in digits); ytop = min(k["y0"] for k in digits); dx1 = max(k["x1"] for k in digits)
    st = [k for k in cs if k["x0"] >= dx1 - 1 and 0.25 * Hh <= k["h"] <= 0.6 * Hh and k["y0"] <= ytop + 0.35 * Hh]
    return len(st), f"{name} mask, {len(digits)} digits, {len(st)} asterisk glyphs"


def run_ocr(ck, pdf, rows, subset_panels, subset_name):
    png = f"{WORK}/{SHEET}_ocr300.png"; L.gs_render(pdf, png, 300, f"{SHEET}_gs300_ocr", max_s=1800)
    pv = f"{WORK}/{SHEET}_ocr_rows.csv"; write_rows_csv(pv, rows); tmp = f"{WORK}/{SHEET}_ocr_checks_tmp.txt"; open(tmp, "w").write("")
    rc, st = L.wd(f"{SHEET}_ocr", [L.PY, OCR_SCRIPT, SHEET, pdf, pv, WORK, tmp], max_s=2400)
    for line in L.wd_out(f"{SHEET}_ocr").splitlines():
        if line.startswith("INFO: R38 OCR"): ck.info(line[6:][:700])
    rb = list(csv.DictReader(open(f"{WORK}/{SHEET}_ocr_readback.csv")))
    dig = [r for r in rb if any(ch.isdigit() for ch in r["string"])]
    n_ok = sum(1 for r in dig if r["status"] in ("exact", "digits_exact")); n_fz = sum(1 for r in dig if r["status"] == "fuzzy")
    misses = [(r["string"], r["panel"], r["x"], r["y"], r["ocr"][:30]) for r in dig if r["status"] not in ("exact", "digits_exact", "fuzzy")]
    ck.log(f"OCR read-back (r38_ocr_check idiom, 300 dpi, passes 1 to 3) of {len(rows)} record strings on the deliverable: {len(dig)} digit-bearing, {n_ok} exact or digit-skeleton exact, {n_fz} label-like fuzzy, {len(misses)} misses (gate 95 percent)", rc == 0 and (n_ok + n_fz) >= 0.95 * len(dig), f"misses {misses[:10]}")
    NUMERIC = re.compile(r"^[\d.,()%*/<>=\-–−≤≥ ]+$")
    def gnorm(s): return re.sub(r"\s+", "", L.norm_text(s).replace("≤", "<").replace("≥", ">").replace("<=", "<").replace(">=", ">"))
    skel = lambda s: "".join(ch for ch in L.norm_text(s).replace("O", "0").replace("o", "0").replace("l", "1").replace("I", "1").replace("S", "5").replace("B", "8") if ch in "0123456789.()-<>")
    sub = [r for r in rb if r["panel"] in subset_panels]; sub_dig = [r for r in sub if any(ch.isdigit() for ch in r["string"])]
    sub_num = [r for r in sub_dig if NUMERIC.match(r["string"])]; sub_lab = [r for r in sub_dig if not NUMERIC.match(r["string"])]
    num_bad = []; num_how = collections.Counter()
    for r in sub_num:
        if r["status"] in ("exact", "digits_exact"): num_how[r["status"]] += 1; continue
        x, y = float(r["x"]), float(r["y"]); wpt = 0.55 * 10.5 * len(r["string"]) + 4
        got = L.ocr_crop(png, (x - 2, y - 9.5, x + wpt, y + 3), 300, psm=7, whitelist="0123456789.()-<>*")
        dsk = lambda s: re.sub(r"[^0-9]", "", s)
        if dsk(got) and dsk(got) == dsk(r["string"]) and (skel(got) == skel(r["string"]) or len(dsk(got)) == len(dsk(r["string"]))):
            n_st, info = star_glyphs(png, r["string"], x, y, 300); n_ocr = got.count("*")
            num_how["exact_recrop_stars_counted" if r["string"].count("*") in (n_st, n_ocr) else "digits_exact_recrop"] += 1; continue
        num_bad.append((r["string"], r["x"], r["y"], r["status"], r["ocr"][:30], got[:30]))
    ck.log(f"OCR read-back of the brief's subset ({subset_name}): {len(sub_num)} numeric strings, every one exact or digit-skeleton exact: {dict(num_how)}", len(sub_num) > 0 and not num_bad, f"{num_bad[:8]}")
    lab_bad = []; lab_how = collections.Counter()
    for r in sub_lab:
        if r["status"] in ("exact", "fuzzy"): lab_how[r["status"]] += 1; continue
        if gnorm(r["ocr"]) and gnorm(r["ocr"]) == gnorm(r["string"]): lab_how["glyph_equal"] += 1; continue
        x, y = float(r["x"]), float(r["y"]); wpt = 0.62 * 10.5 * len(r["string"]) + 6
        got = L.ocr_crop(png, (x - 3, y - 11, x + wpt, y + 4), 300, psm=7)
        if gnorm(got) == gnorm(r["string"]): lab_how["recrop_exact"] += 1; continue
        lab_bad.append((r["string"], r["x"], r["y"], r["ocr"][:30], got[:30]))
    ck.log(f"OCR read-back of the subset's {len(sub_lab)} label-like strings with a digit: {dict(lab_how)}, unread by tesseract {len(lab_bad)} (design labels proven on the text layer; tesseract has no glyph for ≤)", len(lab_bad) <= 2, f"unread {lab_bad[:8]}")
    tit = [r for r in rb if r["string"] == f"Figure {N}"]; ck.log(f"OCR read-back of the title 'Figure {N}': {tit[0]['status'] if tit else 'absent'} ('{tit[0]['ocr'][:20] if tit else ''}')", bool(tit) and tit[0]["status"] in ("exact", "fuzzy"))
    return rb


def raster_gate(ck, v16_150, new_150, allowed, label):
    rd = L.raster_diff(v16_150, new_150, allowed, 150)
    a = np.asarray(Image.open(v16_150).convert("RGB")).astype(int); b = np.asarray(Image.open(new_150).convert("RGB")).astype(int)
    h, w = min(a.shape[0], b.shape[0]), min(a.shape[1], b.shape[1]); d = (np.abs(a[:h, :w] - b[:h, :w]) > 32).any(axis=2); mask = np.zeros_like(d); s = 150 / 72.0
    for (x0, y0, x1, y1) in allowed: mask[max(int((y0 - 1.5) * s), 0):min(int((y1 + 1.5) * s) + 1, h), max(int((x0 - 1.5) * s), 0):min(int((x1 + 1.5) * s) + 1, w)] = True
    out = d & ~mask; lab, n = ndimage.label(out); sizes = sorted((int((lab == i).sum()) for i in range(1, n + 1)), reverse=True)
    locs = [tuple(round(float(v) / s, 1) for v in (o[1].start, o[0].start)) for o in ndimage.find_objects(lab)][:8]
    ck.log(f"raster diff at 150 dpi against {label}, confined to the declared boxes {allowed}: {rd['n_diff']} px differ, {rd['n_outside']} outside in {n} clusters (largest {sizes[0] if sizes else 0} px; isolated 1-2 px anti-aliasing flips tolerated, a cluster of 3 or more or more than 20 flips fails)", (not sizes or sizes[0] <= 2) and rd["n_outside"] <= 20, f"outside clusters at (pt) {locs}; {rd}")
    return rd


# ============================================================================ checks
ck = L.Checks(f"Lane LFIG5b round 42, {SHEET}, {L.now()}. INPUT = V16 {V16PDF} (600 dpi raster, sha256 {L.sha256(V16PDF)[:16]}); V16 record work/r37/coordinator_Main_Fig5_printed_values.csv (146 rows); the LFIG5 deliverable work/r42/Main_Fig5.pdf (sha256 {L.sha256(PREV)[:16]}, the 'before' side: q columns, axis one pitch up). NEW = {RAS} (600 dpi raster of {VEC}, sha256 {rr['sha256_pdf'][:16]}). Changes: one contrast per condition (the B series only, the 1 to 10% line dropped from the rows, the key and the columns), the printed column = the uncorrected P (bold, filled dots where P < 0.05, head P, no q in panel b), 21 single-line rows at one pitch 24.670 from the V13 first row to the V13 last row with the head in its V13 block, panel b's spine, tick labels and title at the V13 y (= panel c's axis), a two-entry key on the left; everything else as LFIG5.")
assert rr["vector_sha256"] == L.sha256(VEC) and rr["sha256_pdf"] == L.sha256(RAS), "raster record does not match the files"
n, w, h = L.gs_pdfinfo(RAS); n0, w0, h0 = L.gs_pdfinfo(V16PDF)
ck.log(f"page box: raster {w} x {h} = vector {W} x {H} = V16 {w0} x {h0} (band top {BAND_TOP} = page height - band height, the V16 slot; letter d at ({clog['letter_d']['x']}, {clog['letter_d']['baseline']}) = the V13 origin), one page, 600 dpi ({rr['width_px']} x {rr['height_px']} px, size check {rr['size_check_within_0.5pt']})", n == 1 and abs(w - w0) < 0.01 and abs(h - h0) < 0.01 and abs(w - W) < 0.01 and abs(h - H) < 0.01 and abs(clog["letter_d"]["baseline"] - 744.76) < 0.01 and rr["dpi"] == 600 and rr["size_check_within_0.5pt"])
ck.info("the composed sheet nests the V13 panel a (454 XObjects): no text extraction on it; the text proof is by the flat page (PyMuPDF) against the numbers files, the raster identity of the unchanged regions, and OCR of the deliverable")
# 2 the selection, recomputed from the refit file
RJ = json.load(open(L.hydrated(A.SPEC["fig5b_refit_v3_json"]))); OC = RJ["outcomes"]
nc = sorted(((k, v["adj_pub_hr"]) for k, v in OC.items() if not v["negative_control"]), key=lambda kv: -kv[1]); ctrl = sorted(((k, v["adj_pub_hr"]) for k, v in OC.items() if v["negative_control"]), key=lambda kv: -kv[1])
exp_rows = [k for k, _ in nc[:18]]; exp_ctrl = [k for k, _ in ctrl]
ck.log(f"18-row selection recomputed from L4_fig5b_refit_v3.json (sha256 {L.sha256(A.SPEC['fig5b_refit_v3_json'])[:16]}, built '{RJ['built']}'): {exp_rows} then the controls {exp_ctrl}; = the builder's rows and order; 19th ratio {nc[18][0]} {nc[18][1]:.4f} below the 18th {nc[17][1]:.4f}; every outcome in the file has at least 40 events (min {min(v['events'] for v in OC.values())})", REC["selection"]["rows"] == exp_rows and REC["selection"]["controls"] == exp_ctrl and min(v["events"] for v in OC.values()) >= 40 and nc[17][1] > nc[18][1])
# 3 every printed string of b and c re-derived from its file and found on the flat page
FL = spans_flat(REC["flat"]); bc = [d for d in DR if d.get("panel") in ("b", "c") and d.get("source_file", "").startswith("/")]
bad_derive = []; bad_find = []; n_ok = 0
for d in bc:
    try:
        if d["rule"] in ("label", "label_line") and d["source_file"].endswith(".json"):
            ptr = d["source_key"].split(" (")[0]; node = resolve(d["source_file"], ptr); key = ptr.split("/")[-1]      # the label IS the outcome key: the node must exist
            ok = isinstance(node, dict) and "adj_pub_hr" in node and ((d["text"] == A.plab(key)) if d["rule"] == "label" else (d["text"] in key and key == d["source_value"]))
            got = key
        else:
            v = resolve(d["source_file"], d["source_key"]); got = RULES[d["rule"]](v) if d["rule"] != "label" else A.plab(str(v)); ok = d["text"] == got
        if not ok: bad_derive.append((d["text"], d["rule"], d["source_key"][:60], str(got)[:40]))
    except Exception as e: bad_derive.append((d["text"], d["rule"], d["source_key"][:60], f"{type(e).__name__}: {str(e)[:60]}"))
    x = d["x"]; ha = d["ha"]
    hit = [s for s in FL if s["text"] == L.norm_text(d["text"]).strip() and abs(s["y"] - d["baseline"]) <= 0.6 and abs((s["x0"] if ha == "left" else s["x1"] if ha == "right" else (s["x0"] + s["x1"]) / 2) - x) <= 0.8]
    if not hit: bad_find.append((d["text"], d["panel"], x, d["baseline"]))
    else: n_ok += 1
ck.log(f"every printed string of panels b and c with a numbers source ({len(bc)} drawn records: labels, HR (95% CI), q) re-derived from its file (json pointer or csv lookup, the printing rule: hr_ci half up at 2 dp, fmt_p two significant figures, label) equals the drawn text, and is on the flat page at its position ({n_ok} found)", not bad_derive and not bad_find, f"derive {bad_derive[:6]} find {bad_find[:6]}")
qb = [d for d in bc if d["rule"] == "fmt_p" and d["panel"] == "b"]
def span_of(d): return next(s for s in FL if s["text"] == d["text"] and abs(s["y"] - d["baseline"]) <= 0.6 and abs(s["x1"] - d["x"]) <= 0.8)
bold_ok = all(("Bold" in span_of(d)["font"]) == (float(d["source_value"]) < 0.05) for d in qb) and all(d["source_key"].endswith("_p") for d in qb)
ck.log(f"panel b P strings ({len(qb)}, every source key an adj_*_p field): bold on the flat page exactly where P < 0.05 ({sum(float(d['source_value']) < 0.05 for d in qb)} bold); bold values {[d['text'] for d in qb if float(d['source_value']) < 0.05]}", bold_ok)
ck.log(f"no 'q' string in panel b (x > 497): {[s['text'] for s in FL if s['x0'] > 497 and s['text'] == 'q']}; panel c keeps its V13 q column (the BH survivors' column, unchanged, {sum(1 for s in FL if s['x0'] < 497 and s['text'] == 'q')} head)", not any(s["x0"] > 497 and s["text"] == "q" for s in FL))
items = A.drawings(REC["flat"]); M = A.markers(items); rows_b = REC["rows_b"]
dots = sorted([m for m in M if m["shape"] in ("circle", "square") and m["cx"] > 600 and m["cy"] > 100], key=lambda m: (m["cy"], m["cx"])); exp_d = sorted(rows_b, key=lambda r: (r["y"], r["x"]))
ck.log(f"panel b markers ({len(dots)} = one per row, 18 circles + 3 squares): filled exactly where P < 0.05 ({sum(r['filled'] for r in rows_b)} filled), at ln(HR) on the V13 axis; every drawn record is the B series", len(dots) == len(exp_d) == 21 and all(r["series"] == "B" for r in rows_b) and all((not m["open"]) == r["filled"] and abs(m["cx"] - r["x"]) < 0.3 for m, r in zip(dots, exp_d)))
key_labels = sorted([s for s in FL if s["text"] in ("Still above 10% on PAP", "1 to 10% on PAP", "Reference: 1% or less on PAP")], key=lambda s: s["y"])
lab_x = REC["key_geometry"]; handles = [it for it in items if "line" in it and abs(it["line"][0] - lab_x["h0"]) < 0.3 and it["line"][1] < 100]
ck.log(f"key measured on the flat page: two entries {[s['text'] for s in key_labels]} at x0 {[round(s['x0'], 2) for s in key_labels]} (labels' x 497.45 + the V13 handle 20.52 + gap 6.84 = {lab_x['label_x']:.2f}), baselines {[s['y'] for s in key_labels]} (the first two V13 key rows), {len(handles)} handles starting at x {lab_x['h0']:.2f} = the condition labels' x; no '1 to 10% on PAP' on the sheet; the block sits above the label column (first row baseline {min(s['y'] for s in FL if abs(s['x0'] - 497.45) < 0.5 and s['size'] > 10 and s['y'] > 100):.2f})", [s["text"] for s in key_labels] == ["Still above 10% on PAP", "Reference: 1% or less on PAP"] and all(abs(s["x0"] - lab_x["label_x"]) < 0.4 for s in key_labels) and len(handles) == 2 and all(abs(s["y"] - y) < 0.05 for s, y in zip(key_labels, (66.23, 77.18))) and not any("1 to 10%" in s["text"] for s in FL))
qh = [s for s in FL if s["text"] == "q" and s["y"] < 80 and s["x0"] > 900]; ph = [s for s in FL if s["text"] == "P" and s["y"] < 80 and s["x0"] > 900]
ck.log(f"column head 'P' at its V13 origin ({ph[0]['x0'] if ph else None}, {ph[0]['y'] if ph else None}); no 'q' head", len(ph) == 1 and not qh and abs(ph[0]["x0"] - 927.17) < 0.4 and abs(ph[0]["y"] - 66.23) < 0.05)
# the axis of panel b at the V13 y: the flat page's drawings (both spines) and the tick labels and title baselines
items_ax = A.drawings(REC["flat"]); spines = sorted(set(round(l["y"], 2) for l in A.hlines(items_ax, min_len=150) if abs(l["width"] - 0.76) < 0.02))
tl = sorted(set(round(s["y"], 2) for s in FL if s["text"] in ("0.5", "1", "2", "4", "8") and s["x0"] > 600)); ttl = sorted(round(s["y"], 2) for s in FL if s["text"] in ("Hazard ratio for a new diagnosis", "vs 1% or less on PAP (95% CI)"))
ck.log(f"panel b axis at the V13 positions on the flat page: spine lines {spines} (b and c on one line 673.85), b tick labels at {tl} (V13 687.49), b axis title at {ttl} (V13 703.83, 715.93)", spines == [673.85] and tl == [687.49] and ttl == [703.83, 715.93])
lab_rows = sorted(round(s["y"], 2) for s in FL if abs(s["x0"] - 497.45) < 0.5 and s["size"] > 10 and s["y"] > 100 and s["text"] not in ("Negative controls", "Ventricular arrhythmia", "or cardiac arrest"))
hr_rows = sorted(round(s["y"], 2) for s in FL if abs(s["x1"] - 896.98) < 0.6 and s["y"] > 100); pitch = REC["layout_b"]["pitch"]; head_y = [round(s["y"], 2) for s in FL if s["text"] == "Negative controls"]
exp_y = [round(116.09 + i * pitch, 2) for i in range(18)] + [round(116.09 + (20 + j) * pitch, 2) for j in range(3)]
ck.log(f"panel b rows: 21 single-line rows (one HR string each) at one pitch {pitch:.3f} pt from the V13 first row 116.09 to the V13 last row {hr_rows[-1]} (V13 658.84), the head at {head_y} (1.5 pitches below the 18th row {hr_rows[17]} and above the first control {hr_rows[18]})", len(hr_rows) == 21 and all(abs(a - b) < 0.05 for a, b in zip(hr_rows, exp_y)) and abs(hr_rows[-1] - 658.84) < 0.05 and len(head_y) == 1 and abs(head_y[0] - (hr_rows[17] + 1.5 * pitch)) < 0.05, f"rows {hr_rows}")
# 4 panel c unchanged; band text
F40 = spans_flat(f"{WORK}/r40/panels_bc_flat.pdf"); kk = lambda s: (s["text"], round(s["x0"], 1), round(s["y"], 1), s["size"], "Bold" in s["font"])
c_new = collections.Counter(kk(s) for s in FL if s["x0"] < 497); c_old = collections.Counter(kk(s) for s in F40 if s["x0"] < 497)
ck.log(f"panel c unchanged: {sum(c_new.values())} spans identical (text, x, baseline, size, face) to the round-40 flat page", c_new == c_old, f"delta {L.delta(c_old, c_new)}")
rec = record_rows(f"{WORK}/r37/coordinator_Main_Fig5_printed_values.csv"); a_rows = [r for r in rec if float(r["x"]) < 497 and float(r["y"]) < 393 and r["string"] != "Figure 5"]; c_rows = [r for r in rec if float(r["x"]) < 497 and 393 <= float(r["y"]) < 725]; d_rows = [r for r in rec if float(r["y"]) >= 725]
miss_c = [r["string"] for r in c_rows if not any(s["text"] == r["string"] and abs(s["y"] - float(r["y"])) <= 0.6 and (abs(s["x0"] - float(r["x"])) <= 1.6 or abs(s["ox"] - float(r["x"])) <= 1.6) for s in FL)]
ck.log(f"every V16 record string of panel c ({len(c_rows)}) on the flat page at its record position", not miss_c, f"missing {miss_c[:6]}")
BS = spans_flat(BAND, 0.0, BAND_TOP)
ck.log(f"band d (full width, four organs) text = the V16 band's strings {sorted(r['string'] for r in d_rows if r['string'] != 'd')}: {[(s['text'], round(s['x0'], 1), round(s['y'], 1)) for s in BS]}", sorted(s["text"] for s in BS) == sorted(r["string"] for r in d_rows if r["string"] != "d"))
# 5 raster diff and crops
allowed = [(0, 0, 140, 22), (497, 0, W, 725), (0, BAND_TOP - 0.5, W, H)]
v16_150 = f"{WORK}/v16_150dpi_gs.png"; L.gs_render(V16PDF, v16_150, 150, f"{SHEET}_gs150_v16", max_s=1800)
rd = raster_gate(ck, v16_150, f"{SD}/{SHEET}_150dpi.png", allowed, "the V16 sheet rendered by the same Ghostscript (panel a, panel c, the letter d and the strip between panel c and the band must be identical)")
prev_150 = f"{WORK}/prev_lfig5_150dpi_gs.png"; L.gs_render(PREV, prev_150, 150, f"{SHEET}_gs150_prev", max_s=1800)
rd2 = raster_gate(ck, prev_150, f"{SD}/{SHEET}_150dpi.png", [(497, 0, W, 725)], "the LFIG5 deliverable rendered by the same Ghostscript, confined to panel b (title, panel a, panel c, the band and the letter d identical to LFIG5)")
# the spine of panel b measured on the render against panel c's: the darkest row of each spine's x-range between y 640 and 700
im150 = np.asarray(Image.open(f"{SD}/{SHEET}_150dpi.png").convert("L")).astype(int); s150 = 150 / 72.0
def spine_y(x0, x1):
    band = im150[int(640 * s150):int(700 * s150), int(x0 * s150):int(x1 * s150)]; prof = band.mean(axis=1); r = int(np.argmin(prof)); return round((int(640 * s150) + r + 0.5) / s150, 2), round(float(prof[r]), 1)
yc, lc = spine_y(170, 350); yb, lb = spine_y(630, 810)
ck.log(f"spine y measured on the 150 dpi render (darkest row over the spine's x-range, y 640 to 700): panel c {yc} pt (row mean level {lc}), panel b {yb} pt ({lb}); difference {abs(yb - yc):.2f} pt (one 150 dpi pixel = 0.48 pt); the flat page draws both at 673.85", abs(yb - yc) < 0.5 and abs(yc - 673.85) < 0.6 and lc < 160 and lb < 160)
new_200 = f"{WORK}/new_200dpi.png"; prev_200 = f"{WORK}/prev_200dpi.png"; L.gs_render(RAS, new_200, 200, f"{SHEET}_gs200_new", max_s=1800); L.gs_render(PREV, prev_200, 200, f"{SHEET}_gs200_prev", max_s=1800)
for name, box in [("01_panel_b_LFIG5_vs_LFIG5b", (480, 20, W, 725)), ("02_axis_alignment_c_vs_b", (100, 630, W, 725)), ("03_key_left", (480, 40, 780, 100)), ("04_controls_head_gap", (480, 500, W, 700))]: ck.info(f"crop OLD (LFIG5) vs NEW at 200 dpi, {name} {box}: {os.path.basename(L.crop_pair(prev_200, new_200, box, 200, f'{VER}/crops/{name}', label_old='OLD (LFIG5, q)', label_new='NEW (LFIG5b, P)'))}")
# 6 OCR
rows = [dict(string="Figure 5", x=36.0, y=14.0, panel="title")] + [dict(string=r["string"], x=float(r["x"]), y=float(r["y"]), panel="a") for r in a_rows] + [dict(string=r["string"], x=float(r["x"]), y=float(r["y"]), panel="c") for r in c_rows]
rows += [dict(string=s["text"], x=round(s["x0"], 2), y=round(s["y"], 2), panel="b") for s in FL if s["x0"] >= 497 and s["text"] != "b"] + [dict(string=s["text"], x=round(s["x0"], 2), y=round(s["y"], 2), panel="band_d") for s in BS] + [dict(string="d", x=36.0, y=744.76, panel="letter")]
if not SKIP_OCR: run_ocr(ck, RAS, rows, ("b", "c"), "every printed string of panels b and c")
# 8 CHANGES: the LFIG5 sheet's panel b (q columns, the same 21 rows) against the new P columns
F42 = spans_flat(f"{WORK}/r42/panels_bc_flat.pdf")
old_lab = sorted([s for s in F42 if abs(s["x0"] - 497.45) < 0.5 and s["size"] > 10 and s["y"] > 100 and s["text"] not in ("Negative controls", "or cardiac arrest")], key=lambda s: s["y"]); old_hr = sorted([s for s in F42 if abs(s["x1"] - 896.98) < 0.6 and s["y"] > 100], key=lambda s: s["y"]); old_p = sorted([s for s in F42 if abs(s["x1"] - 933.2) < 0.6 and s["y"] > 100 and s["size"] < 10], key=lambda s: s["y"])
old = {("Ventricular arrhythmia or cardiac arrest" if s["text"] == "Ventricular arrhythmia" else s["text"]): (old_hr[2 * i]["text"], old_p[2 * i]["text"], old_hr[2 * i + 1]["text"], old_p[2 * i + 1]["text"]) for i, s in enumerate(old_lab)}
assert len(old_lab) == 21 and len(old_hr) == 42 and len(old_p) == 42, (len(old_lab), len(old_hr), len(old_p))
by_out = collections.defaultdict(dict)
for r in rows_b: by_out[r["outcome"]][r["series"]] = r
changes = []; new_order = REC["selection"]["rows"] + REC["selection"]["controls"]
for i, k in enumerate(new_order):
    o = old[k]; rb_ = by_out[k]["B"]
    changes.append(dict(sheet=SHEET, panel="b", item=f"{k}: HR (95% CI), Still above 10% on PAP", old=o[0], new=rb_["s_hr"], kind="unchanged" if o[0] == rb_["s_hr"] else "CHANGED", x="", y="", note=""))
    flip = (float(o[1].lstrip("<")) < 0.05) != (rb_["p"] < 0.05)
    changes.append(dict(sheet=SHEET, panel="b", item=f"{k}: Still above 10% on PAP, q -> P", old=f"q {o[1]}", new=f"P {rb_['s_p']}{' (bold, filled dot)' if rb_['p'] < 0.05 else ''}", kind="column q -> P (uncorrected)", x="", y="", note="dramatic: the printed value and the marker cross 0.05 (q at or above 0.05, P < 0.05)" if flip else ""))
    changes.append(dict(sheet=SHEET, panel="b", item=f"{k}: 1 to 10% on PAP line", old=f"{o[2]}, q {o[3]}", new="removed (one contrast per condition)", kind="removed", x="", y="", note="dramatic: a printed line leaves the sheet"))
changes += [dict(sheet=SHEET, panel="b", item="column head", old="q (LFIG5)", new="P at the V13 origin (927.17, 66.23)", kind="wording", x=927.17, y=66.23, note=""),
            dict(sheet=SHEET, panel="b", item="key", old="three entries (Still above 10% on PAP, 1 to 10% on PAP, Reference: 1% or less on PAP)", new="two entries (Still above 10% on PAP, Reference: 1% or less on PAP) at the first two V13 key rows, on the left over the labels", kind="layout", x=497.45, y=66.23, note=""),
            dict(sheet=SHEET, panel="b", item="x-axis line, tick labels, axis title", old="one pitch up (spine 650.26, LFIG5)", new="the V13 y (spine 673.85 = panel c's axis; tick labels 687.49; title 703.83 and 715.93)", kind="layout", x="", y=673.85, note=""),
            dict(sheet=SHEET, panel="b", item="rows", old="21 two-line rows at the V13 pitch 23.598 from 116.09 (LFIG5)", new=f"21 single-line rows at one pitch {REC['layout_b']['pitch']:.3f} from the V13 first row 116.09 to the V13 last row 658.84, the head in its V13 block (1.5 pitches from each neighbour)", kind="layout", x="", y="", note=""),
            dict(sheet=SHEET, panel="a, c, d, title", item="everything else", old="LFIG5", new="identical (raster diff against LFIG5 confined to panel b)", kind="unchanged", x="", y="", note="")]
L.write_changes(f"{VER}/CHANGES_{SHEET}.csv", changes); ck.info(f"CHANGES_{SHEET}.csv: {len(changes)} rows ({sum('dramatic' in c['note'] for c in changes)} flagged dramatic: rows entering, and P-versus-q flips across 0.05)")
ok = ck.write(f"{VER}/checks.txt")
json.dump(dict(vector_sha256=L.sha256(VEC), raster_sha256=rr["sha256_pdf"], tif_sha256=rr["sha256_tif"], png150_sha256=L.sha256(f"{SD}/{SHEET}_150dpi.png"), page=[w, h], dpi=600, band=BAND, band_sha256=clog["band_sha256"], band_top=BAND_TOP, rows=new_order, p_printed={k: by_out[k]["B"]["s_p"] for k in new_order}, spine_measured={"c": yc, "b": yb}, raster_diff_v16=rd, raster_diff_lfig5=rd2, written=L.now()), open(f"{VER}/verify_record.json", "w"), indent=1, default=str)
print("RESULT", "ALL PASS" if ok else "FAIL"); sys.exit(0 if ok else 1)
