#!$T90_PY
"""Round 54 (2026-09-29), lane L3, Main_Fig3 compose (vector): the rebuilt flat page work/page_abc.pdf (panels a, b, c at the round-54
sizes, page W x BAND_TOP) placed 1:1 on [0, BAND_TOP), the round-54 band placed with its top at BAND_TOP (its page box is read, the
page height = BAND_TOP + band height), the title "Figure 3" and the letters a, b, c, d stamped by TextWriter in Arial Bold 18 pt at the
planner's positions (the letter baselines = top + 0.905 x 18). No V13 base placement, no PyMuPDF text extraction on any nested sheet
(get_text only on the flat page_abc.pdf and on the band). Then a Ghostscript pdfwrite re-distillation of the composed sheet
(work/Main_Fig3_flat.pdf, the round-44 gs_flatten recipe: one page, no downsampling, colours unchanged, fonts embedded) for the
coordinator's census, and the 150 dpi Ghostscript render.
usage: 02_compose_r54.py [--band PATH] [--out PATH]   (default band = ROUND54/lanes/LSKETCH/work/build/band_3d_r54.pdf)"""
import argparse, gc, json, os, subprocess, sys, time
import fitz
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lmain_lib as L

SHEET = "Main_Fig3"; SD = f"{L.LANE}/{SHEET}"; WORK = f"{SD}/work"
ap = argparse.ArgumentParser(); ap.add_argument("--band", default=f"{L.SKETCH}/band_3d_r54.pdf"); ap.add_argument("--out", default=f"{SD}/{SHEET}.pdf"); a = ap.parse_args()
BAND = L.hydrated(a.band); OUT = a.out; PAGE = L.hydrated(f"{WORK}/page_abc.pdf")
DJ = json.load(open(f"{WORK}/{SHEET}_drawn.json")); Y = DJ["plan"]["y"]; PC = DJ["plan"]["c"]
W_DECL, H_PAGE = (float(v) for v in DJ["page"]); BAND_TOP = float(DJ["band_top"]); assert abs(H_PAGE - BAND_TOP) < 1e-6, (H_PAGE, BAND_TOP)
W = 952.73; assert abs(W_DECL - W) < 0.01, W_DECL
ASC = 0.905                                                       # Arial Bold ascent: letter baseline = top + ASC x size


def spans_flat(pdf):
    d = fitz.open(pdf); out = []
    for b in d[0].get_text("rawdict")["blocks"]:
        if b["type"] != 0: continue
        for l in b["lines"]:
            for s in l["spans"]:
                t = "".join(c["c"] for c in s["chars"]).replace("\xa0", " ")
                if t.strip(): out.append(dict(text=t, bbox=[round(v, 3) for v in s["bbox"]], origin=[round(s["chars"][0]["origin"][0], 3), round(s["chars"][0]["origin"][1], 3)], size=round(s["size"], 2), font=s["font"]))
    d.close(); gc.collect(); return out


pd_ = fitz.open(PAGE); pp = pd_[0]; assert abs(pp.rect.width - W) < 0.05 and abs(pp.rect.height - H_PAGE) < 0.05, (pp.rect, W, H_PAGE)
PSP = spans_flat(PAGE); assert all(Y["title_base"] + 6.0 <= s["bbox"][1] and s["bbox"][3] <= BAND_TOP - 2.0 for s in PSP), "rebuilt page has text outside the panel region"
bd = fitz.open(BAND); bp = bd[0]; bw, bh = bp.rect.width, bp.rect.height
assert abs(bw - W) < 0.05, ("band width must equal the sheet width", bw, W)
BSP = spans_flat(BAND)
H_NEW = round(BAND_TOP + bh, 3)
assert H_NEW <= 1260.0, ("the sheet is taller than the round-54 gate", H_NEW)
out = fitz.open(); pg = out.new_page(width=W, height=H_NEW)
pg.show_pdf_page(fitz.Rect(0, 0, W, BAND_TOP), pd_, 0, clip=fitz.Rect(0, 0, W, BAND_TOP))
pg.show_pdf_page(fitz.Rect(0, BAND_TOP, W, BAND_TOP + bh), bd, 0)
font = fitz.Font(fontfile=L.ARIALB)
LETTERS = [("a", 8.0, Y["a_top"]), ("b", 8.0, Y["b_top"]), ("c", PC["x0"], Y["b_top"]), ("d", 8.0, Y["d_top"])]
tw = fitz.TextWriter(pg.rect); stamped = []
for ch, x, top in LETTERS:
    base = round(top + ASC * L.LETTER_SIZE, 3); assert base < BAND_TOP or ch != "d" or top + L.LETTER_SIZE * 1.117 < BAND_TOP, ("letter d must sit above the band", top, BAND_TOP)
    tw.append(fitz.Point(x, base), ch, font=font, fontsize=L.LETTER_SIZE); stamped.append((ch, x, base))
assert Y["d_top"] + L.LETTER_SIZE * 1.117 < BAND_TOP, ("letter d must sit above the band", Y["d_top"], BAND_TOP)
tw.write_text(pg, color=(0, 0, 0))
tt = fitz.TextWriter(pg.rect); TITLE_X = 8.0; tt.append(fitz.Point(TITLE_X, L.TITLE_BASELINE), "Figure 3", font=font, fontsize=L.TITLE_SIZE); tt.write_text(pg, color=L.INK_RGB)
os.makedirs(os.path.dirname(OUT), exist_ok=True)
out.save(OUT, garbage=3, deflate=True); out.close(); bd.close(); pd_.close(); gc.collect()
n, w, h = L.gs_pdfinfo(OUT); assert n == 1 and abs(w - W) < 0.01 and abs(h - H_NEW) < 0.01, (n, w, h, H_NEW)

# the flat re-distillation (Ghostscript pdfwrite, the round-44 gs_flatten recipe) for the coordinator's census, under the watchdog
FLAT = f"{WORK}/{SHEET}_flat.pdf"
if os.path.exists(FLAT): os.remove(FLAT)
rc, st = L.wd("flatten_r54", [L.GS, "-dSAFER", "-dBATCH", "-dNOPAUSE", "-dQUIET", "-sDEVICE=pdfwrite", "-dCompatibilityLevel=1.7", "-dPDFSETTINGS=/prepress",
                              "-dDownsampleColorImages=false", "-dDownsampleGrayImages=false", "-dDownsampleMonoImages=false", "-dAutoFilterColorImages=false",
                              "-dAutoFilterGrayImages=false", "-dColorImageFilter=/FlateEncode", "-dGrayImageFilter=/FlateEncode", "-dColorConversionStrategy=/LeaveColorUnchanged",
                              "-dEmbedAllFonts=true", "-dSubsetFonts=true", "-dPreserveAnnots=false", f"-sOutputFile={FLAT}", OUT], max_s=600)
assert rc == 0 and os.path.exists(FLAT) and os.path.getsize(FLAT) > 0, ("gs pdfwrite re-distillation failed", st, L.wd_out("flatten_r54")[-400:])
nf, wf, hf = L.gs_pdfinfo(FLAT); assert nf == 1 and abs(wf - W) < 0.01 and abs(hf - H_NEW) < 0.01, ("flat page box", nf, wf, hf)
# the 150 dpi render of the vector (Ghostscript png16m under the watchdog)
PNG150 = f"{SD}/{SHEET}_150dpi.png"; st150 = L.gs_render(OUT, PNG150, 150, "render_r54_150", max_s=600)

log = {"sheet": SHEET, "round": 54, "page_abc": PAGE, "page_abc_sha256": L.sha256(PAGE), "band": BAND, "band_sha256": L.sha256(BAND), "band_page": [round(bw, 3), round(bh, 3)], "band_text": [s["text"] for s in BSP],
       "placements": [{"what": "rebuilt page (a, b, c)", "rect": [0, 0, W, BAND_TOP]}, {"what": "band d (round 54)", "rect": [0, BAND_TOP, W, H_NEW]}],
       "band_top": BAND_TOP, "letters_stamped": stamped, "letter_size": L.LETTER_SIZE, "letter_ascent": ASC,
       "title": {"text": "Figure 3", "origin": [TITLE_X, L.TITLE_BASELINE], "size": L.TITLE_SIZE, "font": "Arial Bold (TextWriter)", "color": "#1a1d21"},
       "page": [w, h], "page_v30": [952.73, 1043.821], "bytes": os.path.getsize(OUT), "output": OUT, "output_sha256": L.sha256(OUT),
       "flat": FLAT, "flat_sha256": L.sha256(FLAT), "flat_page": [wf, hf], "flat_status": st, "flat_recipe": "gs pdfwrite CompatibilityLevel 1.7 PDFSETTINGS /prepress, no downsampling, FlateEncode images, LeaveColorUnchanged, EmbedAllFonts, SubsetFonts",
       "png150": PNG150, "png150_status": st150, "png150_sha256": L.sha256(PNG150), "written": L.now()}
json.dump(log, open(f"{WORK}/compose_log.json", "w"), indent=1)
print(f"wrote {OUT}: page {w} x {h} (gs PDFINFO), band {os.path.basename(BAND)} {bw:.2f} x {bh:.3f} at top {BAND_TOP}, letters {stamped} at {L.LETTER_SIZE:g} pt, title at ({TITLE_X}, {L.TITLE_BASELINE}) {L.TITLE_SIZE:g} pt")
print(f"flat re-distillation {FLAT}: {nf} page {wf} x {hf} ({os.path.getsize(FLAT)} bytes, {st}); render {PNG150} ({st150})")
