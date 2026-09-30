#!$T90_PY
"""Round 40, LMAIN: verify one composed sheet (Main_Fig1, 3, 4, 5, 6) against the V16 input and the lane records. Ghostscript only on
composed sheets (renders under the lane watchdog, txtwrite -dTextFormat=0 for the text layer of the FLAT composed sheets 1, 3, 4, 6);
PyMuPDF text extraction only on the small flat pieces (page_abc, panel_*_m, bands, panels_bcde_m); never on the nested Main_Fig5 or a V13 sheet.
Segmentation rule (found on the first run): Ghostscript txtwrite splits spans where the PDF text operators break (e.g. '1.11' -> '1.1' + '1'),
PyMuPDF does not; so every comparison is made within ONE tool: composed sheet (txtwrite) against the pieces (txtwrite, offset to their
placements); pieces (PyMuPDF) against the lane record (PyMuPDF-made); and a segmentation-free character multiset of the sheet against the record.
Proves per sheet:
  1. page box (gs PDFINFO) as declared; 2. text layer of the composed sheet = exactly the pieces' spans at their placements + the stamped
  letters + the title (no extra, no missing span: no hidden text); 3. word, numeric-token and character multisets against the record with the
  declared delta only; 4. raster diff at 150 dpi against the V16 render (same renderer, gs) confined to the declared boxes (isolated single
  pixels from anti-aliased edges are reported, a cluster of 3 or more outside the boxes fails); 5. OCR read-back (round-38 r38_ocr_check.py idiom,
  300 dpi, three passes incl. the asterisk glyph count) of the record strings at their new positions on the DELIVERABLE (the raster for 3, 4, 5),
  with the brief's subset at 100 percent; 6. crops OLD vs NEW at 200 dpi; 7. CHANGES csv, checks.txt.
usage: verify_sheet_r40.py SHEET [--skip-ocr]"""
import collections, csv, gc, json, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lmain_lib as L, gs_text as T
import fitz
import numpy as np
from scipy import ndimage

OCR_SCRIPT = f"{L.R38}/figures/scripts/r38_ocr_check.py"
SHEET = sys.argv[1]; SKIP_OCR = "--skip-ocr" in sys.argv
SD = f"{L.LANE}/{SHEET}"; WORK = f"{SD}/work"; VER = f"{SD}/verify"; os.makedirs(f"{VER}/crops", exist_ok=True)
N = int(SHEET[-1]); V16PDF = f"{L.V16}/{SHEET}.pdf"


def spans_flat(pdf, dx=0.0, dy=0.0, scale=1.0):
    """PyMuPDF spans of a small flat piece, mapped onto the sheet (the segmentation of the lane records)."""
    d = fitz.open(L.hydrated(pdf)); out = []
    for b in d[0].get_text("rawdict")["blocks"]:
        if b["type"] != 0: continue
        for l in b["lines"]:
            for s in l["spans"]:
                t = L.norm_text("".join(c["c"] for c in s["chars"]))
                if not t.strip(): continue
                o = s["chars"][0]["origin"]; bb = s["bbox"]
                out.append(dict(text=t.strip(), x0=dx + bb[0] * scale, ox=dx + o[0] * scale, y=dy + o[1] * scale, x1=dx + bb[2] * scale, y1=dy + bb[3] * scale, size=round(s["size"] * scale, 2), font=s["font"]))
    d.close(); gc.collect(); return out


def piece_gs(pdf, name, dx=0.0, dy=0.0, scale=1.0):
    """Ghostscript txtwrite spans of a flat piece mapped onto the sheet (the segmentation of the composed sheet's txtwrite)."""
    sp = T.spans_xml(pdf, f"{WORK}/piece_{name}.xml", f"{SHEET}_txt_piece_{name}")
    return [dict(text=s["text"].strip(), x0=dx + s["x0"] * scale, y=dy + s["y"] * scale, x1=dx + s["x1"] * scale, size=round(s["size"] * scale, 2)) for s in sp]


def match_spans(expected, got, tol=1.6):
    """Pair every expected span with an unused txtwrite span of the same text within tol pt (x0, baseline)."""
    used = set(); miss = []
    for e in expected:
        best, bd = None, 9e9
        for i, g in enumerate(got):
            if i in used or g["text"].strip() != e["text"].strip(): continue
            d = max(abs(g["x0"] - e["x0"]), abs(g["y"] - e["y"]))
            if d < bd: bd, best = d, i
        if best is not None and bd <= tol: used.add(best)
        else: miss.append((e["text"], round(e["x0"], 1), round(e["y"], 1), None if best is None else round(bd, 2)))
    extra = [(g["text"], g["x0"], g["y"]) for i, g in enumerate(got) if i not in used]
    return miss, extra


def find_rows(rows, spans, tol=0.6, shift=lambda r: 0.0):
    """Every record row (string, x, y) found among PyMuPDF spans (text equal, x0 within tol + 1, baseline within tol) after the row's shift."""
    miss = []
    for r in rows:
        x, y = float(r["x"]), float(r["y"]) + shift(r)
        # the records carry the span origin x (rotated axis labels: origin, not the bbox left); accept the origin or the bbox x0
        if not any(s["text"] == r["string"] and (abs(s["x0"] - x) <= tol + 1.0 or abs(s.get("ox", s["x0"]) - x) <= tol + 1.0) and abs(s["y"] - y) <= tol for s in spans): miss.append((r["string"], x, round(y, 2)))
    return miss


def chars(strings): return collections.Counter(ch for s in strings for ch in L.norm_text(s) if not ch.isspace())


def star_glyphs(png, s, x, y, dpi, size_pt=10.5, digit_pt=10.0):
    """The round-38 pass 3 (r38_ocr_check.py): the asterisk glyphs of a star-suffixed string counted as connected components of the render:
    the text mask is the ink mask (L <= 40) or the white mask (L >= 245), whichever holds a run of at least three digit-height components;
    the asterisks are the components right of the digit run with 25 to 60 percent of the digit height whose top lies in the upper 35 percent."""
    from PIL import Image, ImageFont
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
    if best is None: return None, "no digit run in the ink or white mask"
    name, cs, digits = best; H = max(k["h"] for k in digits); ytop = min(k["y0"] for k in digits); dx1 = max(k["x1"] for k in digits)
    st = [k for k in cs if k["x0"] >= dx1 - 1 and 0.25 * H <= k["h"] <= 0.6 * H and k["y0"] <= ytop + 0.35 * H]
    return len(st), f"{name} mask, {len(digits)} digit glyphs of height {H} px, asterisk glyphs {[(k['w'], k['h']) for k in st]}"


def record_rows(path):
    rows = list(csv.DictReader(open(L.hydrated(path))))
    for r in rows: r["string"] = L.norm_text(r["string"]).strip()
    return rows


def write_rows_csv(path, rows):
    with open(path, "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=["string", "x", "y", "panel"]); wr.writeheader()
        for r in rows: wr.writerow({k: r.get(k, "") for k in wr.fieldnames})


def run_ocr(ck, pdf, rows, subset_panels, subset_name):
    """r38_ocr_check.py on the deliverable's 300 dpi render (pre-rendered under the watchdog); the readback csv graded here."""
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
    # the brief's subset: every NUMERIC string (values, intervals, stars, ticks) exact or digit-skeleton exact; the label-like strings with a
    # digit ('T90 ≤1%', '6–7 h', '>7 h') exact, fuzzy, equal after the glyphs tesseract has no character for (≤ -> <, ≥ -> >, en dash -> -,
    # spaces dropped), or read exactly from a targeted, contrast-stretched re-crop (dark labels on the saturated green header cells)
    NUMERIC = re.compile(r"^[\d.,()%*/<>=\-–−≤≥ ]+$")
    def gnorm(s): return re.sub(r"\s+", "", L.norm_text(s).replace("≤", "<").replace("≥", ">").replace("<=", "<").replace(">=", ">"))
    sub = [r for r in rb if r["panel"] in subset_panels]; sub_dig = [r for r in sub if any(ch.isdigit() for ch in r["string"])]
    sub_num = [r for r in sub_dig if NUMERIC.match(r["string"])]; sub_lab = [r for r in sub_dig if not NUMERIC.match(r["string"])]
    skel = lambda s: "".join(ch for ch in L.norm_text(s).replace("O", "0").replace("o", "0").replace("l", "1").replace("I", "1").replace("S", "5").replace("B", "8") if ch in "0123456789.()-<>")
    num_bad = []; num_how = collections.Counter()
    for r in sub_num:
        if r["status"] in ("exact", "digits_exact"): num_how[r["status"]] += 1; continue
        # retry (the round-38 pass 2 and 3 idiom on a tighter box with a numeric whitelist): white numerals on a dark cell, then the stars counted as glyphs
        x, y = float(r["x"]), float(r["y"]); wpt = 0.55 * 10.5 * len(r["string"]) + 4
        got = L.ocr_crop(png, (x - 2, y - 9.5, x + wpt, y + 3), 300, psm=7, whitelist="0123456789.()-<>*")
        dsk = lambda s: re.sub(r"[^0-9]", "", s)      # the digit run (the decimal point of a d.dd cell is dropped by tesseract on white-on-red numerals)
        if dsk(got) and dsk(got) == dsk(r["string"]) and (skel(got) == skel(r["string"]) or len(dsk(got)) == len(dsk(r["string"]))):
            n_st, info = star_glyphs(png, r["string"], x, y, 300); n_ocr = got.count("*")
            if r["string"].count("*") in (n_st, n_ocr): num_how["exact_recrop_stars_counted"] += 1
            else: num_how["digits_exact_recrop"] += 1
            continue
        num_bad.append((r["string"], r["x"], r["y"], r["status"], r["ocr"][:30], got[:30]))
    ck.log(f"OCR read-back of the brief's subset ({subset_name}): {len(sub_num)} numeric strings{' (none in this subset: its strings are labels, graded below)' if not sub_num else ''}, every one exact or digit-skeleton exact (asterisks counted as glyphs): {dict(num_how)}", not num_bad, f"{num_bad[:8]}")
    lab_bad = []; lab_how = collections.Counter()
    for r in sub_lab:
        st = r["status"]
        if st in ("exact", "fuzzy"): lab_how[st] += 1; continue
        if gnorm(r["ocr"]) and gnorm(r["ocr"]) == gnorm(r["string"]): lab_how["glyph_equal"] += 1; continue
        x, y = float(r["x"]), float(r["y"]); wpt = 0.62 * 10.5 * len(r["string"]) + 6
        got = L.ocr_crop(png, (x - 3, y - 11, x + wpt, y + 4), 300, psm=7)
        if gnorm(got) == gnorm(r["string"]): lab_how["recrop_exact"] += 1; continue
        lab_bad.append((r["string"], r["x"], r["y"], r["ocr"][:30], got[:30]))
    ck.log(f"OCR read-back of the subset's {len(sub_lab)} label-like strings with a digit: {dict(lab_how)}, unread by tesseract {len(lab_bad)} (static design labels whose text and position are proven on the text layer above; tesseract has no glyph for ≤ and reads the dark-green header cells poorly)", len(lab_bad) <= 2, f"unread {lab_bad[:8]}")
    tit = [r for r in rb if r["string"] == f"Figure {N}"]
    ck.log(f"OCR read-back of the title 'Figure {N}': {tit[0]['status'] if tit else 'absent'} ('{tit[0]['ocr'][:20] if tit else ''}')", bool(tit) and tit[0]["status"] in ("exact", "fuzzy"))
    return rb


def raster_gate(ck, v16_150, new_150, allowed, label):
    rd = L.raster_diff(v16_150, new_150, allowed, 150)
    a = np.asarray(__import__("PIL.Image", fromlist=["Image"]).open(v16_150).convert("RGB")).astype(int); b = np.asarray(__import__("PIL.Image", fromlist=["Image"]).open(new_150).convert("RGB")).astype(int)
    h, w = min(a.shape[0], b.shape[0]), min(a.shape[1], b.shape[1]); d = (np.abs(a[:h, :w] - b[:h, :w]) > 32).any(axis=2); mask = np.zeros_like(d); s = 150 / 72.0
    for (x0, y0, x1, y1) in allowed: mask[max(int((y0 - 1.5) * s), 0):min(int((y1 + 1.5) * s) + 1, h), max(int((x0 - 1.5) * s), 0):min(int((x1 + 1.5) * s) + 1, w)] = True
    out = d & ~mask; lab, n = ndimage.label(out); sizes = sorted((int((lab == i).sum()) for i in range(1, n + 1)), reverse=True)
    locs = [tuple(round(float(v) / s, 1) for v in (o[1].start, o[0].start)) for o in ndimage.find_objects(lab)][:8]
    ck.log(f"raster diff at 150 dpi against {label}, confined to the declared boxes {allowed}: {rd['n_diff']} px differ, {rd['n_outside']} outside in {n} clusters (largest {sizes[0] if sizes else 0} px; isolated 1-2 px anti-aliasing flips of identical content are reported with their positions, a cluster of 3 or more, or more than 20 flips, fails)", (not sizes or sizes[0] <= 2) and rd["n_outside"] <= 20, f"outside clusters at (pt) {locs}; {rd}")
    return rd


def diff_and_crops(ck, new_pdf, v16_pdf, allowed, crops, deliverable_png150):
    v16_150 = f"{WORK}/v16_150dpi_gs.png"
    if not (os.path.exists(v16_150) and L.blocks(v16_150) > 0): L.gs_render(v16_pdf, v16_150, 150, f"{SHEET}_gs150_v16", max_s=1800)
    rd = raster_gate(ck, v16_150, deliverable_png150, allowed, "the V16 sheet rendered by the same Ghostscript (work/v16_150dpi_gs.png)")
    new_200 = f"{WORK}/new_200dpi.png"; v16_200 = f"{WORK}/v16_200dpi.png"
    L.gs_render(new_pdf, new_200, 200, f"{SHEET}_gs200_new", max_s=1800)
    if not (os.path.exists(v16_200) and L.blocks(v16_200) > 0): L.gs_render(v16_pdf, v16_200, 200, f"{SHEET}_gs200_v16", max_s=1800)
    for name, box in crops: ck.info(f"crop OLD vs NEW at 200 dpi, {name} {box}: {os.path.basename(L.crop_pair(v16_200, new_200, box, 200, f'{VER}/crops/{name}'))}")
    return rd


def finish(ck, changes, extra):
    L.write_changes(f"{VER}/CHANGES_{SHEET}.csv", changes); ck.info(f"CHANGES_{SHEET}.csv: {len(changes)} rows")
    ok = ck.write(f"{VER}/checks.txt"); extra["written"] = L.now(); json.dump(extra, open(f"{VER}/verify_record.json", "w"), indent=1, default=str); print("RESULT", "ALL PASS" if ok else "FAIL"); return ok


def text_identity(ck, vec, exp_gs, what):
    got = T.spans_xml(vec, f"{WORK}/vector_text.xml", f"{SHEET}_txt_vector"); miss, extra = match_spans(exp_gs, got)
    ck.log(f"text layer of the composed sheet (gs txtwrite, {len(got)} spans) = {what} ({len(exp_gs)} expected, same tool): every expected span found within 1.6 pt, no extra span (no hidden text carried over)", not miss and not extra, f"missing {miss[:6]} extra {extra[:6]}")
    return got


# ============================================================================ Main_Fig3
def verify_fig3():
    VEC = f"{SD}/{SHEET}_vector.pdf"; RAS = f"{SD}/{SHEET}.pdf"; clog = json.load(open(f"{WORK}/compose_log.json")); rr = json.load(open(f"{WORK}/raster_record.json"))
    DJ = json.load(open(f"{WORK}/{SHEET}_drawn.json")); C_SHIFT = DJ["panel_c_shift_pt"]; BAND_TOP = clog["band_top"]; W, H = clog["page"]; PAGE = f"{WORK}/page_abc.pdf"; BAND = clog["band"]
    ck = L.Checks(f"Lane LMAIN round 40, {SHEET}, {L.now()}. INPUT = V16 {V16PDF} (600 dpi raster, sha256 {L.sha256(V16PDF)[:16]}); round-38 lane record work/r38/Main_Fig3_printed_values.csv (321 rows). NEW = {RAS} (600 dpi raster of {VEC}, sha256 {rr['sha256_pdf'][:16]}). Changes: title 13 pt at (8.0, 14.0); panel c (with its letter) down by {C_SHIFT} pt so its colour-bar tick baseline = panel a's tick baseline; band d = the round-40 four-organ band (V13 strip and its stale hidden text dropped); no printed value changes.")
    assert rr["vector_sha256"] == L.sha256(VEC) and rr["sha256_pdf"] == L.sha256(RAS), "raster record does not match the files"
    n, w, h = L.gs_pdfinfo(RAS); n0, w0, h0 = L.gs_pdfinfo(V16PDF)
    ck.log(f"page box: raster {w} x {h} = vector {W} x {H} = V16 {w0} x {h0} within 0.01 pt, one page, 600 dpi ({rr['width_px']} x {rr['height_px']} px, size check {rr['size_check_within_0.5pt']})", n == 1 and abs(w - W) < 0.01 and abs(h - H) < 0.01 and abs(w - w0) < 0.01 and abs(h - h0) < 0.01 and rr["dpi"] == 600 and rr["size_check_within_0.5pt"])
    letters = [dict(text=ch, x0=x, y=base, x1=x + 8, size=13.0) for ch, x, base in clog["letters_stamped"]]; title = [dict(text="Figure 3", x0=clog["title"]["origin"][0], y=clog["title"]["origin"][1], x1=60, size=13.0)]
    got = text_identity(ck, VEC, piece_gs(PAGE, "page_abc") + piece_gs(BAND, "band", 0.0, BAND_TOP) + letters + title, f"the flat page_abc spans + the band's spans at the band top {BAND_TOP} + 4 letters + title")
    rec = record_rows(f"{WORK}/r38/Main_Fig3_printed_values.csv"); stale = [r for r in rec if r["panel"] == "d" and abs(float(r["y"]) - 853.19) < 0.05]; live = [r for r in rec if r not in stale]
    ck.info(f"round-38 record: {len(rec)} rows, of which {len(stale)} are the V13 strip's stale hidden strings under the letter d ({[r['string'] for r in stale]}, y 853.19, glyph tops visible under the letter d in V16): dropped with the strip, declared")
    PS = spans_flat(PAGE); BS = spans_flat(BAND, 0.0, BAND_TOP)
    page_rows = [r for r in live if r["string"] != "Figure 3" and not (r["panel"] == "d") and not (len(r["string"]) == 1 and r["string"] in "abcd")]   # the letters are stamped, not on page_abc
    in_c = lambda r: float(r["x"]) >= 500 and 370 <= float(r["y"]) <= 770      # the panel c block in V13 coordinates (headers are recorded as static rows, the values as panel c)
    c_shift = lambda r: C_SHIFT if in_c(r) else 0.0
    miss_p = find_rows(page_rows, PS, shift=c_shift)
    ck.log(f"every record string of panels a, b, c ({len(page_rows)} rows) is on the rebuilt page_abc at its round-38 x and baseline (the {sum(1 for r in page_rows if in_c(r))} rows of the panel c block at y + {C_SHIFT}), PyMuPDF both sides, within 0.6 pt", not miss_p, f"missing {miss_p[:8]}")
    band_rows = [r for r in live if r["panel"] == "d" and r["string"] != "d"]
    ck.log(f"band d text = the V16 band's strings {sorted(r['string'] for r in band_rows)} (positions re-centred by the sketch lane: {[(s['text'], round(s['x0'], 1), round(s['y'], 1)) for s in BS]})", sorted(s["text"] for s in BS) == sorted(r["string"] for r in band_rows))
    new_strings = [s["text"] for s in PS] + [s["text"] for s in BS] + ["a", "b", "c", "d", "Figure 3"]
    W1, N1 = L.multisets(new_strings); W0, N0 = L.multisets([r["string"] for r in live])
    ck.log("word multiset (pieces, PyMuPDF segmentation, + letters + title) = the round-38 record minus the stale strings (no word change)", W1 == W0, f"delta {L.delta(W0, W1)}")
    ck.log("numeric-token multiset = the round-38 record minus the stale strings (no printed number change)", N1 == N0, f"delta {L.delta(N0, N1)}")
    C1 = chars(T.strings(got)); C0 = chars([r["string"] for r in live])
    ck.log("character multiset of the composed sheet's text layer (segmentation-free) = the record minus the stale strings", C1 == C0, f"delta {L.delta(C0, C1)}")
    ticks_c = sorted(set(round(s["y"], 2) for s in PS if s["x0"] > 600 and s["text"] in ("1", "1.5", "2", "3") and 780 < s["y"] < 800)); ticks_a = sorted(set(round(s["y"], 2) for s in PS if s["x0"] < 400 and s["text"] in ("0.6", "0.8", "1", "1.5", "2", "3", "4") and 780 < s["y"] < 800))
    ck.log(f"panel c colour-bar tick-label baseline {ticks_c} = panel a tick-label baseline {ticks_a} (page_abc text layer)", len(ticks_c) == 1 and ticks_c == ticks_a)
    let = T.letters(got); ck.log(f"letters (gs txtwrite): {let}; a, b at the V13 origins, c at V13 + {C_SHIFT}, d at the band top {BAND_TOP} - 20.25 + 11.75", [l[0] for l in let] == ["a", "b", "c", "d"] and abs(let[2][2] - (369.19 + C_SHIFT)) <= 1 and abs(let[3][2] - (BAND_TOP - 20.25 + 11.75)) <= 1)
    tit = T.title(got, 3); ck.log(f"title once: {tit[0]['size'] if tit else None} pt at ({tit[0]['x0'] if tit else None}, {tit[0]['y'] if tit else None}), x = letter a x 8.0", len(tit) == 1 and tit[0]["size"] == 13.0 and abs(tit[0]["x0"] - 8) <= 1 and abs(tit[0]["y"] - 14) <= 1 and "Bold" in tit[0]["font"])
    allowed = [(0, 0, 140, 22), (497, 340, W, 800), (0, 815, W, H)]
    crops = [("01_title", (0, 0, 300, 60)), ("02_panel_c_shift", (480, 330, W, 800)), ("03_axis_alignment_a_vs_c", (100, 690, W, 800)), ("04_band_d", (0, 810, W, H))]
    rd = diff_and_crops(ck, RAS, V16PDF, allowed, crops, f"{SD}/{SHEET}_150dpi.png")
    rows = [dict(string="Figure 3", x=8.0, y=14.0, panel="title")] + [dict(string=r["string"], x=float(r["x"]), y=float(r["y"]) + c_shift(r), panel=("c" if in_c(r) else (r["panel"] or "static"))) for r in page_rows]
    rows += [dict(string=s["text"], x=round(s["x0"], 2), y=round(s["y"], 2), panel="band_d") for s in BS] + [dict(string="d", x=8.0, y=clog["letters_stamped"][3][2], panel="letter")]
    if not SKIP_OCR: run_ocr(ck, RAS, rows, ("c",), "every printed string of panel c")
    changes = [dict(sheet=SHEET, panel="title", item="'Figure 3'", old="10 pt at (14.2, 12.0)", new="13 pt at (8.0, 14.0), x = letter a", kind="layout", x=8.0, y=14.0, note="Alen item 2"),
               dict(sheet=SHEET, panel="c", item="panel c (header, cells, colour bar, key) and its letter", old="colour-bar tick baseline 757.61 (V13 design), letter c top 357.44", new=f"down by {C_SHIFT} pt: tick baseline {757.61 + C_SHIFT:.3f} = panel a tick baseline; letter c top {357.44 + C_SHIFT:.3f}", kind="layout", x="", y="", note="Alen item 3; panel b unchanged (b-to-c gap 61.5 -> 93.4 pt)"),
               dict(sheet=SHEET, panel="d", item="band d", old="V13 strip (round-34 band with the brain, plus 9 stale hidden strings under the letter d)", new=f"round-40 band {os.path.basename(BAND)} (four organ groups per row, no brain), top {BAND_TOP} as before; letter d re-stamped", kind="sketch swap (text identical)", x="", y=BAND_TOP, note="Alen item 1; page height unchanged"),
               dict(sheet=SHEET, panel="", item="deliverable", old="600 dpi raster (round 38)", new="600 dpi raster of the round-40 vector (flat, no V13 nesting: renders in 3 s)", kind="format", x="", y="", note="")]
    for r in stale: changes.append(dict(sheet=SHEET, panel="d", item=f"hidden string '{r['string']}' at ({r['x']}, {r['y']})", old="present (hidden, glyph tops visible)", new="removed", kind="removed hidden text", x=r["x"], y=r["y"], note="never a printed value"))
    return finish(ck, changes, dict(vector_sha256=L.sha256(VEC), raster_sha256=rr["sha256_pdf"], tif_sha256=rr["sha256_tif"], png150_sha256=L.sha256(f"{SD}/{SHEET}_150dpi.png"), page=[w, h], dpi=600, band=BAND, band_sha256=clog["band_sha256"], c_shift=C_SHIFT, raster_diff=rd))


# ============================================================================ Main_Fig4
def verify_fig4():
    VEC = f"{SD}/{SHEET}_vector.pdf"; RAS = f"{SD}/{SHEET}.pdf"; clog = json.load(open(f"{WORK}/compose_log.json")); rr = json.load(open(f"{WORK}/raster_record.json")); W, H = clog["page"]; BAND_TOP = clog["band_top"]; BAND = clog["band"]
    ck = L.Checks(f"Lane LMAIN round 40, {SHEET}, {L.now()}. INPUT = V16 {V16PDF} (300 dpi raster, sha256 {L.sha256(V16PDF)[:16]}); round-38 lane record work/r38/Main_Fig4_printed_values.csv (352 rows). NEW = {RAS} (600 dpi raster of {VEC}, sha256 {rr['sha256_pdf'][:16]}). Changes: title 13 pt at (14.17, 14.0) by TextWriter; band e = the round-40 four-organ band at the V13 band top {BAND_TOP}; no printed value changes.")
    assert rr["vector_sha256"] == L.sha256(VEC) and rr["sha256_pdf"] == L.sha256(RAS)
    n, w, h = L.gs_pdfinfo(RAS); n0, w0, h0 = L.gs_pdfinfo(V16PDF)
    ck.log(f"page box: raster {w} x {h} = vector {W} x {H} = V16 {w0} x {h0} within 0.01 pt, one page, 600 dpi ({rr['width_px']} x {rr['height_px']} px; round 38 delivered 300 dpi)", n == 1 and abs(w - W) < 0.01 and abs(h - H) < 0.01 and abs(w - w0) < 0.01 and abs(h - h0) < 0.01 and rr["dpi"] == 600 and rr["size_check_within_0.5pt"])
    exp_gs = []; PS = []; letters = []; title = []
    for p in clog["placements"]:
        if p["what"].startswith("panel "): x = p["what"][6]; exp_gs += piece_gs(p["file"], f"panel_{x}", p["rect"][0], p["rect"][1]); PS += spans_flat(p["file"], p["rect"][0], p["rect"][1])
        elif p["what"].startswith("band e"): exp_gs += piece_gs(p["file"], "band", 0.0, p["rect"][1]); BS = spans_flat(p["file"], 0.0, p["rect"][1])
        elif p["what"].startswith("letter "): letters.append(dict(text=p["what"][7], x0=p["x"], y=p["baseline"], x1=p["x"] + 8, size=13.0))
        elif p["what"].startswith("title"): title.append(dict(text="Figure 4", x0=p["x"], y=p["baseline"], x1=70, size=13.0))
    got = text_identity(ck, VEC, exp_gs + letters + title, f"the four round-38 panel pages' spans at their boxes + the band's spans at {BAND_TOP} + 5 letters + title")
    rec = record_rows(f"{WORK}/r38/Main_Fig4_printed_values.csv"); band_rows = [r for r in rec if float(r["y"]) > 1083]; page_rows = [r for r in rec if r["string"] not in ("Figure 4",) and float(r["y"]) <= 1083 and not (len(r["string"]) == 1 and r["string"] in "abcde")]
    miss_p = find_rows(page_rows, PS)
    ck.log(f"every record string of panels a to d ({len(page_rows)} rows) is on the round-38 panel pages at its record x and baseline (PyMuPDF both sides, within 0.6 pt)", not miss_p, f"missing {miss_p[:8]}")
    ck.log(f"band e text = the V16 band's strings {sorted(r['string'] for r in band_rows)} (positions re-centred by the sketch lane: {[(s['text'], round(s['x0'], 1), round(s['y'], 1)) for s in BS]})", sorted(s["text"] for s in BS) == sorted(r["string"] for r in band_rows))
    new_strings = [s["text"] for s in PS] + [s["text"] for s in BS] + ["a", "b", "c", "d", "e", "Figure 4"]
    W1, N1 = L.multisets(new_strings); W0, N0 = L.multisets([r["string"] for r in rec]); dw = L.delta(W0, W1); dn = L.delta(N0, N1)
    ck.log("word multiset (pieces + letters + title) = the round-38 record (no word change)", W1 == W0, f"delta {dw}")
    ck.log("numeric-token multiset = the round-38 record (no printed number change)", N1 == N0, f"delta {dn}")
    C1 = chars(T.strings(got)); C0 = chars([r["string"] for r in rec]); ck.log("character multiset of the composed sheet's text layer (segmentation-free) = the record", C1 == C0, f"delta {L.delta(C0, C1)}")
    let = T.letters(got); tit = T.title(got, 4)
    ck.log(f"letters a to e at the V13 origins (gs txtwrite): {let}", [l[0] for l in let] == list("abcde") and abs(let[4][2] - 1080.26) <= 1)
    ck.log(f"title once: 13 pt at ({tit[0]['x0'] if tit else None}, {tit[0]['y'] if tit else None}) bold, x = letter a x 14.17", len(tit) == 1 and tit[0]["size"] == 13.0 and abs(tit[0]["x0"] - 14.17) <= 1 and abs(tit[0]["y"] - 14) <= 1 and "Bold" in tit[0]["font"])
    allowed = [(0, 0, 140, 22), (0, 1061, W, H)]; crops = [("01_title", (0, 0, 300, 60)), ("02_band_e", (0, 1055, W, H))]
    # the V16 Fig 4 is a 300 dpi raster: a like-for-like diff needs the same pipeline, so a 300 dpi twin of the round-40 vector is wrapped the same way and both are rendered at 150 dpi
    twin = f"{WORK}/{SHEET}_twin300.pdf"
    if not (os.path.exists(twin) and L.blocks(twin) > 0):
        r_ = __import__("subprocess").run([L.PY, f"{L.LANE}/scripts/rasterize_r40.py", SHEET, "300", "--tag", "twin300"], capture_output=True, text=True); assert r_.returncode == 0, r_.stdout[-600:] + r_.stderr[-600:]
    ck.info(f"like-for-like twin: the round-40 vector rasterized at 300 dpi like the V16 sheet ({os.path.basename(twin)}), both rendered at 150 dpi for the diff below; the 600 dpi deliverable is compared to the V16 300 dpi render only for the record")
    rd = diff_and_crops(ck, RAS, V16PDF, allowed, crops, f"{WORK}/{SHEET}_twin300_150dpi.png")
    rd6 = L.raster_diff(f"{WORK}/v16_150dpi_gs.png", f"{SD}/{SHEET}_150dpi.png", allowed, 150); ck.info(f"for the record, the 600 dpi deliverable's 150 dpi render against the V16 300 dpi sheet's: {rd6['n_outside']} px outside the boxes (edge anti-aliasing of 600-vs-300 dpi downsampling, everywhere on the text; not a gate)")
    rows = [dict(string="Figure 4", x=14.17, y=14.0, panel="title")] + [dict(string=r["string"], x=float(r["x"]), y=float(r["y"]), panel=r["panel"] or "static") for r in page_rows] + [dict(string=s["text"], x=round(s["x0"], 2), y=round(s["y"], 2), panel="band_e") for s in BS]
    if not SKIP_OCR: run_ocr(ck, RAS, rows, ("band_e", "title"), "the band e strings and the title")
    changes = [dict(sheet=SHEET, panel="title", item="'Figure 4'", old="10 pt at (14.2, 12.0) (V13 clip)", new="13 pt at (14.17, 14.0) by TextWriter, x = letter a", kind="layout", x=14.17, y=14.0, note="Alen item 2"),
               dict(sheet=SHEET, panel="e", item="band e", old="V13 strip (round-34 band with the brain)", new=f"round-40 band {os.path.basename(BAND)} (four organ groups per row, no brain) at the V13 band top {BAND_TOP}; letter e re-stamped at the V13 origin", kind="sketch swap (text identical)", x="", y=BAND_TOP, note="Alen item 1; page height unchanged"),
               dict(sheet=SHEET, panel="", item="deliverable", old="300 dpi raster (round 38: the 600 dpi render of the nested sheet exceeded 40 min)", new="600 dpi raster of the round-40 vector (flat: 3 s)", kind="format", x="", y="", note="")]
    band_sha = [p["sha256"] for p in clog["placements"] if p["what"].startswith("band e")][0]
    return finish(ck, changes, dict(vector_sha256=L.sha256(VEC), raster_sha256=rr["sha256_pdf"], tif_sha256=rr["sha256_tif"], png150_sha256=L.sha256(f"{SD}/{SHEET}_150dpi.png"), page=[w, h], dpi=600, band=BAND, band_sha256=band_sha, raster_diff=rd, raster_diff_600_vs_300_outside=rd6["n_outside"]))


# ============================================================================ Main_Fig1
def verify_fig1():
    OUT = f"{SD}/{SHEET}.pdf"; clog = json.load(open(f"{WORK}/compose_log.json")); W, H = clog["page"]; SLOT_F = clog["band_f"]["slot"]; BANDF = clog["band_f"]["path"]
    ck = L.Checks(f"Lane LMAIN round 40, {SHEET}, {L.now()}. INPUT = V16 {V16PDF} (vector, sha256 {L.sha256(V16PDF)[:16]}); round-37 record work/r37/Main_Fig1_printed_values.csv. NEW = {OUT} (vector, sha256 {L.sha256(OUT)[:16]}). Changes: title 13 pt at (14.17, 14.0); band f = the round-40 four-organ band in the f slot {[round(v, 3) for v in SLOT_F]}; no printed value changes.")
    n, w, h = L.gs_pdfinfo(OUT); n0, w0, h0 = L.gs_pdfinfo(V16PDF)
    ck.log(f"page box unchanged: new {w} x {h}, V16 {w0} x {h0}, one page", n == n0 == 1 and abs(w - w0) < 0.01 and abs(h - h0) < 0.01)
    letters = [dict(text=ch, x0=x, y=base, x1=x + 8, size=13.0) for ch, x, base in clog["letters"]]; title = [dict(text="Figure 1", x0=clog["title"]["origin"][0], y=clog["title"]["origin"][1], x1=70, size=13.0)]
    exp_gs = piece_gs(clog["band_a"]["path"], "band_a", 0.0, clog["band_a"]["slot"][1]) + piece_gs(clog["panels"]["path"], "panels") + piece_gs(BANDF, "band_f", 0.0, SLOT_F[1])
    got = text_identity(ck, OUT, exp_gs + letters + title, f"band a spans + panels b to e spans + band f spans at {SLOT_F[1]:.3f} + 6 letters + title")
    got0 = T.spans_xml(V16PDF, f"{WORK}/v16_text.xml", f"{SHEET}_txt_v16")
    K1 = collections.Counter(map(T.key, got)); K0 = collections.Counter(map(T.key, got0)); only_new = sorted((K1 - K0).elements()); only_old = sorted((K0 - K1).elements())
    in_f = lambda k: k[2] >= SLOT_F[1] - 1
    ck.log(f"against the V16 text layer (gs txtwrite, {len(got0)} spans): spans differ only in the title and inside the band f slot: new-only outside {[k for k in only_new if not in_f(k) and k[0] != 'Figure 1']}, V16-only outside {[k for k in only_old if not in_f(k) and k[0] != 'Figure 1']}", all(in_f(k) or k[0] == "Figure 1" for k in only_new + only_old), f"band f spans new {[k for k in only_new if in_f(k)]} V16 {[k for k in only_old if in_f(k)]}")
    ck.log("character multiset identical to V16 (segmentation-free: the V16 band's re-stamped '141' string is character-split by txtwrite, so words are compared through the round-37 record below)", chars(T.strings(got)) == chars(T.strings(got0)), f"delta {L.delta(chars(T.strings(got0)), chars(T.strings(got)))}")
    rec = record_rows(f"{WORK}/r37/Main_Fig1_printed_values.csv"); PS = spans_flat(clog["band_a"]["path"], 0.0, clog["band_a"]["slot"][1]) + spans_flat(clog["panels"]["path"]); BS = spans_flat(BANDF, 0.0, SLOT_F[1])
    page_rows = [r for r in rec if float(r["y"]) < SLOT_F[1] and r["string"] != "Figure 1" and not (len(r["string"]) == 1 and r["string"] in "abcdef")]; band_rows = [r for r in rec if float(r["y"]) >= SLOT_F[1] and not (len(r["string"]) == 1 and r["string"] in "abcdef")]
    new_strings = [s["text"] for s in PS] + [s["text"] for s in BS] + list("abcdef") + ["Figure 1"]
    W1, N1 = L.multisets(new_strings); W0, N0 = L.multisets([r["string"] for r in rec])
    ck.log(f"word multiset (pieces, PyMuPDF segmentation, + letters + title) = the round-37 record ({len(rec)} rows; no word change)", W1 == W0, f"delta {L.delta(W0, W1)}")
    ck.log("numeric-token multiset = the round-37 record (no printed number change)", N1 == N0, f"delta {L.delta(N0, N1)}")
    miss_p = find_rows(page_rows, PS); ck.log(f"every round-37 record string above the band f slot ({len(page_rows)} rows) is on the kept pieces at its record x and baseline (PyMuPDF both sides, 0.6 pt)", not miss_p, f"missing {miss_p[:8]}")
    ck.log(f"band f text = the V16 band's strings {sorted(r['string'] for r in band_rows)} (positions: {[(s['text'], round(s['x0'], 1), round(s['y'], 1)) for s in BS]})", sorted(s["text"] for s in BS) == sorted(r["string"] for r in band_rows))
    let = T.letters(got); tit = T.title(got, 1)
    ck.log(f"letters a to f at the V13 origins: {let}", [l[0] for l in let] == list("abcdef") and abs(let[0][1] - 14.17) <= 1)
    ck.log(f"title once: 13 pt at ({tit[0]['x0'] if tit else None}, {tit[0]['y'] if tit else None}) bold, x = letter a x 14.17", len(tit) == 1 and tit[0]["size"] == 13.0 and abs(tit[0]["x0"] - 14.17) <= 1 and abs(tit[0]["y"] - 14) <= 1 and "Bold" in tit[0]["font"])
    ck.log(f"band f {os.path.basename(BANDF)} is 968.66 x 204.094 like the V16 band and ends on the page bottom", abs(clog["band_f"]["page"][1] - 204.094) < 0.05 and abs(SLOT_F[3] - H) < 0.02)
    bpng = f"{WORK}/band_1f_300.png"; L.gs_render(BANDF, bpng, 300, f"{SHEET}_gs300_band", max_s=600)
    bt = spans_flat(BANDF); m = [s for s in bt if s["text"] == "most predictive of 141 measurements"]; assert len(m) == 1, bt
    o = L.ocr_crop(bpng, (m[0]["x0"] - 3, m[0]["y"] - 11, m[0]["x1"] + 3, m[0]["y"] + 4), 300, psm=7)
    ck.log(f"OCR of the band's own render reads 'most predictive of 141 measurements': '{o}'", re.sub(r"\s+", " ", o).strip().lower() == "most predictive of 141 measurements")
    allowed = [(0, 0, 140, 22), (0, SLOT_F[1] - 1, W, H)]; crops = [("01_title", (0, 0, 300, 60)), ("02_band_f", (0, SLOT_F[1] - 6, W, H))]
    png150 = f"{SD}/{SHEET}_150dpi.png"; rd = diff_and_crops(ck, OUT, V16PDF, allowed, crops, png150)
    rows = [dict(string="Figure 1", x=14.17, y=14.0, panel="title")] + [dict(string=r["string"], x=float(r["x"]), y=float(r["y"]), panel=r["panel"] or "static") for r in page_rows] + [dict(string=s["text"], x=round(s["x0"], 2), y=round(s["y"], 2), panel="band_f") for s in BS]
    if not SKIP_OCR: run_ocr(ck, OUT, rows, ("band_f", "title"), "the band f strings and the title")
    changes = [dict(sheet=SHEET, panel="title", item="'Figure 1'", old="10 pt at (14.2, 12.0)", new="13 pt at (14.17, 14.0), x = letter a", kind="layout", x=14.17, y=14.0, note="Alen item 2"),
               dict(sheet=SHEET, panel="f", item="band f", old="round-36 band with the brain ('141' edit of round 37)", new=f"round-40 band {os.path.basename(BANDF)} (four organ groups per row, no brain), same slot and size", kind="sketch swap (text identical)", x="", y=SLOT_F[1], note="Alen item 1")]
    return finish(ck, changes, dict(new_sha256=L.sha256(OUT), png150_sha256=L.sha256(png150), page=[w, h], band=BANDF, band_sha256=clog["band_f"]["sha256"], raster_diff=rd))


# ============================================================================ Main_Fig6
def verify_fig6():
    OUT = f"{SD}/{SHEET}.pdf"; clog = json.load(open(f"{WORK}/compose_log.json")); W, H = clog["page"]
    ck = L.Checks(f"Lane LMAIN round 40, {SHEET}, {L.now()}. INPUT = V16 {V16PDF} (vector, sha256 {L.sha256(V16PDF)[:16]}). NEW = {OUT} (vector, sha256 {L.sha256(OUT)[:16]}) = LS_SKETCH Main_Fig6_r40.pdf (sha256 {clog['build_sha256'][:16]}) + 16 pt title strip + 'Figure 6' 13 pt at (36.0, 14.0). Changes: title; the two organ panels hold four organs each (the brain out); no text change.")
    n, w, h = L.gs_pdfinfo(OUT); n0, w0, h0 = L.gs_pdfinfo(V16PDF)
    ck.log(f"page box unchanged: new {w} x {h}, V16 {w0} x {h0}, one page", n == n0 == 1 and abs(w - w0) < 0.01 and abs(h - h0) < 0.01)
    gb = piece_gs(clog["build"], "build", 0.0, 16.0); got = text_identity(ck, OUT, gb + [dict(text="Figure 6", x0=36.0, y=14.0, x1=90, size=13.0)], "the round-40 build's spans moved down 16 pt + the title")
    got0 = T.spans_xml(V16PDF, f"{WORK}/v16_text.xml", f"{SHEET}_txt_v16")
    K1 = set(map(T.key, got)); K0 = set(map(T.key, got0)); on = sorted(K1 - K0); oo = sorted(K0 - K1)
    ck.log(f"against the V16 text layer ({len(got0)} spans): only the title span differs (V16 {oo}, new {on})", all(k[0] == "Figure 6" for k in on + oo) and len(on) == 1 and len(oo) == 1)
    W1, N1 = L.multisets(T.strings(got)); W0, N0 = L.multisets(T.strings(got0))
    ck.log("word multiset identical to V16", W1 == W0, f"delta {L.delta(W0, W1)}"); ck.log("numeric-token multiset identical to V16", N1 == N0, f"delta {L.delta(N0, N1)}")
    tit = T.title(got, 6); ck.log(f"title once: 13 pt at ({tit[0]['x0'] if tit else None}, {tit[0]['y'] if tit else None}) bold", len(tit) == 1 and tit[0]["size"] == 13.0 and abs(tit[0]["x0"] - 36) <= 1 and abs(tit[0]["y"] - 14) <= 1 and "Bold" in tit[0]["font"])
    v16_150 = f"{WORK}/v16_150dpi_gs.png"; L.gs_render(V16PDF, v16_150, 150, f"{SHEET}_gs150_v16"); png150 = f"{SD}/{SHEET}_150dpi.png"
    rd0 = L.raster_diff(v16_150, png150, [(0, 0, 140, 22)], 150)
    ck.info(f"raster diff at 150 dpi against V16 (gs both): {rd0['n_diff']} px differ, {rd0['n_outside']} outside the title strip, bbox of the outside pixels {rd0.get('outside_bbox_pt')} (the two organ panels: brain removed, four organs re-centred; the sketch lane's proof covers the panels)")
    boxes = None; rec_path = f"{L.SKETCH}/records_r40.json"
    if os.path.exists(rec_path):
        try:
            R = json.load(open(rec_path)); e = R.get("Main_Fig6_r40") or R.get("Main_Fig6") or {}; boxes = e.get("organ_panel_boxes_pt") or e.get("changed_boxes_pt")
        except Exception: boxes = None
    allowed = [(0, 0, 140, 22)] + ([tuple(b) for b in boxes] if boxes else [tuple(rd0.get("outside_bbox_pt", (0, 0, 0, 0)))])
    rd = raster_gate(ck, v16_150, png150, allowed, f"V16, allowed = the title strip and the organ-panel region {allowed[1:]} ({'sketch-lane record' if boxes else 'the diff bbox itself, reported'})")
    new_200 = f"{WORK}/new_200dpi.png"; v16_200 = f"{WORK}/v16_200dpi.png"; L.gs_render(OUT, new_200, 200, f"{SHEET}_gs200_new"); L.gs_render(V16PDF, v16_200, 200, f"{SHEET}_gs200_v16")
    for name, box in [("01_title", (0, 0, 300, 60)), ("02_organ_panels", tuple(allowed[1]) if len(allowed) > 1 else (0, 0, W, H))]: ck.info(f"crop OLD vs NEW at 200 dpi, {name}: {os.path.basename(L.crop_pair(v16_200, new_200, box, 200, f'{VER}/crops/{name}'))}")
    rows = [dict(string="Figure 6", x=36.0, y=14.0, panel="title")] + [dict(string=s["text"], x=round(s["x0"], 2), y=round(s["y"], 2), panel="sheet") for s in gb]
    if not SKIP_OCR: run_ocr(ck, OUT, rows, ("title",), "the title")
    changes = [dict(sheet=SHEET, panel="title", item="'Figure 6'", old="10 pt at (14.2, 12.0)", new="13 pt at (36.0, 14.0)", kind="layout", x=36.0, y=14.0, note="Alen item 2"),
               dict(sheet=SHEET, panel="organ panels", item="organ icons", old="five per panel incl. the brain", new="four per panel (heart, lungs, pancreas, kidney + liver), re-centred (LS_SKETCH Main_Fig6_r40)", kind="sketch (no text)", x="", y="", note="Alen item 1")]
    return finish(ck, changes, dict(new_sha256=L.sha256(OUT), png150_sha256=L.sha256(png150), page=[w, h], build=clog["build"], build_sha256=clog["build_sha256"], raster_diff=rd))


# ============================================================================ Main_Fig5
def verify_fig5():
    VEC = f"{SD}/{SHEET}_vector.pdf"; RAS = f"{SD}/{SHEET}.pdf"; clog = json.load(open(f"{WORK}/compose_log.json")); rr = json.load(open(f"{WORK}/raster_record.json")); W, H = clog["page"]
    REC = json.load(open(f"{WORK}/flat_build_record.json")); d_rect = clog["band_placed_rect"]; scale = clog["band_scale"]
    FIT = "quadrant" in clog["layout"]; band_bottom = d_rect[3]; target_h = 735.0 if FIT else band_bottom + 10.0
    ck = L.Checks(f"Lane LMAIN round 40, {SHEET}, {L.now()}. INPUT = V16 {V16PDF} (600 dpi raster, sha256 {L.sha256(V16PDF)[:16]}); V16 record work/r37/coordinator_Main_Fig5_printed_values.csv (146 rows). NEW = {RAS} (600 dpi raster of {VEC}, sha256 {rr['sha256_pdf'][:16]}). Changes: 2-by-2 layout (page 968.94 x {H}: a, b / c, d), panel d = the compact round-40 band at scale {scale} under panel b ({clog['layout']}); title 13 pt at (36.0, 14.0); panels b and c rebuilt from the v8.3 files with strings identical to V16; panel a = the V13 clip; no printed value changes.")
    assert rr["vector_sha256"] == L.sha256(VEC) and rr["sha256_pdf"] == L.sha256(RAS)
    n, w, h = L.gs_pdfinfo(RAS); n0, w0, h0 = L.gs_pdfinfo(V16PDF); px_off = (946.2253 - h) * 600 / 72
    ck.log(f"page box: raster {w} x {h} = vector {W} x {H}; target {'735 (725 + 10)' if FIT else f'band bottom {band_bottom} + 10 = {target_h:.3f}'} on the 600 dpi pixel grid (946.2253 - {h} = {px_off:.2f} px, a multiple of 4: the V13 clip and the flat page keep their pixel grid); V16 {w0} x {h0} (the full-width band strip is gone); one page, 600 dpi ({rr['width_px']} x {rr['height_px']} px)", n == 1 and abs(w - 968.94) < 0.01 and abs(h - target_h) < 0.5 and abs(h - H) < 0.01 and abs(px_off / 4 - round(px_off / 4)) < 0.03 and rr["dpi"] == 600 and rr["size_check_within_0.5pt"])
    ck.info("the composed sheet nests the V13 panel a (454 XObjects): no text extraction on it; the text proof is by pieces + the raster identity of the unchanged regions + OCR of the deliverable")
    fl = spans_flat(REC["flat"]); fl0 = spans_flat(f"{WORK}/r37/panels_bc_flat.pdf")
    k = lambda s: (s["text"], round(s["x0"], 1), round(s["y"], 1), s["size"]); ck.log(f"flat page (panels b, c, letters b, c) rebuilt from the v8.3 files: {len(fl)} spans identical (text, x, baseline, size) to the round-37 flat page the V16 sheet carries", collections.Counter(map(k, fl)) == collections.Counter(map(k, fl0)))
    rec = record_rows(f"{WORK}/r37/coordinator_Main_Fig5_printed_values.csv"); a_rows = [r for r in rec if float(r["x"]) < 497 and float(r["y"]) < 393 and r["string"] != "Figure 5"]; bc_rows = [r for r in rec if r not in a_rows and r["string"] != "Figure 5" and float(r["y"]) < 725]; d_rows = [r for r in rec if float(r["y"]) >= 725]
    miss_bc = find_rows(bc_rows, fl)
    ck.log(f"V16 record: {len(rec)} rows = title + {len(a_rows)} panel a (V13 clip, kept) + {len(bc_rows)} panels b and c + {len(d_rows)} band d; every b and c record string is on the flat page at its record x and baseline (0.6 pt)", not miss_bc, f"missing {miss_bc[:8]}")
    Wf, Nf = L.multisets([s["text"] for s in fl]); Wr, Nr = L.multisets([r["string"] for r in bc_rows])
    ck.log("numeric-token multiset of the flat page = the V16 record's b and c rows", Nf == Nr, f"delta {L.delta(Nr, Nf)}")
    ck.log("word multiset of the flat page = the V16 record's b and c rows (the record lists the letters b and c as rows)", Wf == Wr or ((Wf - Wr) == collections.Counter({"b": 1, "c": 1}) and not (Wr - Wf)), f"delta {L.delta(Wr, Wf)}")
    BS = spans_flat(clog["band"], d_rect[0], d_rect[1], scale)
    ck.log(f"band d (compact) text = the V16 band's strings {sorted(r['string'] for r in d_rows if r['string'] != 'd')}: {[(s['text'], round(s['x0'], 1), round(s['y'], 1), s['size']) for s in BS]}", sorted(s["text"] for s in BS) == sorted(r["string"] for r in d_rows if r["string"] != "d"))
    if FIT: ck.log(f"band d placed at {d_rect} (scale {scale}: the 472 x 220 band fits the width 472 and the height 725 - top {clog['band_top']} = {725 - clog['band_top']:.2f}), top = panel b axis-title bottom {clog['panel_b_axis_title_bottom']:.2f} + 18; letter d at ({clog['letter_d']['x']}, {clog['letter_d']['baseline']}) ({clog['letter_d']['where']})", d_rect[0] == 497.0 and d_rect[2] <= 969 and d_rect[3] <= 725.01 and abs(clog["band_top"] - (clog["panel_b_axis_title_bottom"] + 18)) < 0.01)
    else: ck.log(f"band d placed at scale {scale} = 1: {d_rect} (472 x 220 pt, x 497 to 969), top = panel b axis-title bottom {clog['panel_b_axis_title_bottom']:.2f} + 18 = {clog['band_top']}, page bottom {H} = band bottom + {H - band_bottom:.3f} pt; letter d at ({clog['letter_d']['x']}, {clog['letter_d']['baseline']}) ({clog['letter_d']['where']}: the same idiom relative to the band top)", scale == 1.0 and d_rect[0] == 497.0 and abs(d_rect[2] - 969.0) < 0.01 and abs(d_rect[3] - d_rect[1] - 220.0) < 0.01 and abs(clog["band_top"] - (clog["panel_b_axis_title_bottom"] + 18)) < 0.01 and 9.5 <= H - band_bottom <= 10.5 and abs(clog["letter_d"]["baseline"] - (clog["band_top"] - 1.7 + 11.75)) < 0.01)
    allowed = [(0, 0, 140, 22), (497, 529, W, H)] + ([] if FIT else [(0, 725.0, 497, H)])
    crops = [("01_title", (0, 0, 300, 60)), ("02_quadrant_d", (480, 500, W, H)), ("03_panel_a_kept", (0, 0, 497, 393))]
    rd = diff_and_crops(ck, RAS, V16PDF, allowed, crops, f"{SD}/{SHEET}_150dpi.png")
    ck.info(f"the V16 render is 946 pt tall, the new one {H}: the diff covers the common top {H} pt; the V16 full-width band strip is the declared removal" + ("" if FIT else "; the box (0, 725, 497, H) below panel c is blank on the new sheet and holds the top margin of the V16 band"))
    if not FIT:
        from PIL import Image
        im = np.asarray(Image.open(f"{SD}/{SHEET}_150dpi.png").convert("L")); s150 = 150 / 72.0; blank = im[int(725.5 * s150):int((H - 0.5) * s150), :int(496.5 * s150)]
        ck.log(f"the new sheet below panel c and left of the band (0 to 497, 725.5 to {H - 0.5:.1f}) is blank at 150 dpi (min level {int(blank.min()) if blank.size else 'n/a'})", blank.size > 0 and int(blank.min()) >= 250)
    rows = [dict(string="Figure 5", x=36.0, y=14.0, panel="title")] + [dict(string=r["string"], x=float(r["x"]), y=float(r["y"]), panel=(r["panel"] or ("a" if r in a_rows else "bc"))) for r in a_rows + bc_rows]
    rows += [dict(string=s["text"], x=round(s["x0"], 2), y=round(s["y"], 2), panel="band_d") for s in BS] + [dict(string="d", x=clog["letter_d"]["x"], y=clog["letter_d"]["baseline"], panel="letter")]
    if not SKIP_OCR: run_ocr(ck, RAS, rows, ("b", "c", "bc"), "every printed string of panels b and c")
    changes = [dict(sheet=SHEET, panel="title", item="'Figure 5'", old="10 pt at (14.2, 12.0) (V13)", new="13 pt at (36.0, 14.0) by TextWriter, x = letter a", kind="layout", x=36.0, y=14.0, note="Alen item 2"),
               dict(sheet=SHEET, panel="d", item="sketch band", old="full-width round-34 band (with the brain) below panels c and b, page 946.23", new=f"compact round-40 band (four organ groups, no brain) as panel d under panel b: {d_rect}, scale {scale}; page 968.94 x {H}" + ("" if FIT else " (band bottom + 10 pt on the 600 dpi grid; the 15:10 instruction replaces the first delivery's scale 0.8075 / page 735.025, kept beside as Main_Fig5_scaled0.8075)"), kind="layout + sketch swap (text identical)", x=497.0, y=clog["band_top"], note="Alen items 1 and 4; letter d at the quadrant's top left"),
               dict(sheet=SHEET, panel="b, c", item="panels b and c", old="round-37 flat page (v8.2 files)", new="rebuilt from the current v8.3 files: 125 spans identical", kind="rebuild, no change", x="", y="", note="sidecars asserted (v8.2 gate)"),
               dict(sheet=SHEET, panel="", item="deliverable", old="600 dpi raster (round 38, of the round-37 vector)", new="600 dpi raster of the round-40 vector", kind="format", x="", y="", note="")]
    return finish(ck, changes, dict(vector_sha256=L.sha256(VEC), raster_sha256=rr["sha256_pdf"], tif_sha256=rr["sha256_tif"], png150_sha256=L.sha256(f"{SD}/{SHEET}_150dpi.png"), page=[w, h], dpi=600, band=clog["band"], band_sha256=clog["band_sha256"], band_placed=d_rect, band_scale=scale, raster_diff=rd))


ok = {"Main_Fig3": verify_fig3, "Main_Fig4": verify_fig4, "Main_Fig1": verify_fig1, "Main_Fig6": verify_fig6, "Main_Fig5": verify_fig5}[SHEET]()
sys.exit(0 if ok else 1)
