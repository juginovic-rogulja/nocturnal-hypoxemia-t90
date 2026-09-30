#!$T90_PY
"""Round 49 (2026-09-26), lane L3 copy of the round-40 compose: the letters a, b, c, d and the title "Figure 3" are stamped at 14 pt
(13 + 1, L.LETTER_SIZE and L.TITLE_SIZE) at the SAME origins (baselines), the band is the round-49 band (notes 12 pt, captions 11 pt,
same page size), the page_abc piece is the round-49 build (+1 pt text). Everything else as in round 40. Original docstring follows.
Round 40, LMAIN, Main_Fig3 compose (vector): NO V13 base placement any more. The sheet = the rebuilt flat page work/page_abc.pdf
(panels a, b, c; c shifted down by C_SHIFT) placed 1:1 on [16, 772 + DELTA), the round-40 band placed with its top at
772 + DELTA + 22.52 (the V13 gap between the panels' region end 772 and the V13 band top 794.52), the title "Figure 3" Arial Bold
13 pt at (8.0, 14.0) (x = the letter a's x0) by TextWriter, the letters a, b at the V13 (x, top + 11.75), c at the V13 top + C_SHIFT,
d at the V13 x and the band top - 20.25 (the V13 letter-d top 774.27 against the V13 band top 794.52) + 11.75. Page height = band top +
band height. The V13 letter-d strip (and its stale hidden text line) is gone. Memory rule: no PyMuPDF text extraction on any nested
sheet; the V13 letters come from the round-38 probe cache (v13_cache); get_text only on the flat page_abc.pdf and the band.
usage: 02_compose_r40.py [--band PATH] [--out PATH]   (default band = LS_SKETCH build band_3d_r40.pdf)"""
import argparse, gc, json, os, sys
import fitz
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # round 49: this lane's copies of lmain_lib and v13_cache
import lmain_lib as L, v13_cache

SHEET = "Main_Fig3"; SD = f"{L.LANE}/{SHEET}"; WORK = f"{SD}/work"
ap = argparse.ArgumentParser(); ap.add_argument("--band", default=f"{L.SKETCH}/band_3d_r49.pdf"); ap.add_argument("--out", default=f"{WORK}/{SHEET}_r49_vector.pdf"); a = ap.parse_args()
BAND = L.hydrated(a.band); OUT = a.out; PAGE = L.hydrated(f"{WORK}/page_abc.pdf")
DJ = json.load(open(f"{WORK}/{SHEET}_drawn.json")); DELTA = float(DJ["delta_pt"]); C_SHIFT = float(DJ["panel_c_shift_pt"])
Y_TITLE, Y_D = 16.0, 772.0; Y_D_NEW = Y_D + DELTA
V13_BAND_TOP, V13_LETTER_D_TOP = 794.52, 774.27        # round-30 compose log (band located at [794.52, 992.945]) and the V13 letter record
GAP_BAND = round(V13_BAND_TOP - Y_D, 3); D_OFF = round(V13_LETTER_D_TOP - V13_BAND_TOP, 3)   # 22.52 and -20.25
W, H13 = v13_cache.page(SHEET); LB = {t[0]: t for t in v13_cache.letters_lane(SHEET)}; assert set(LB) == set("abcd"), LB
assert abs(W - 952.73) < 0.01, W


def spans_flat(pdf):
    d = fitz.open(pdf); out = []
    for b in d[0].get_text("rawdict")["blocks"]:
        if b["type"] != 0: continue
        for l in b["lines"]:
            for s in l["spans"]:
                t = "".join(c["c"] for c in s["chars"]).replace("\xa0", " ")
                if t.strip(): out.append(dict(text=t, bbox=[round(v, 3) for v in s["bbox"]], origin=[round(s["chars"][0]["origin"][0], 3), round(s["chars"][0]["origin"][1], 3)], size=round(s["size"], 2), font=s["font"]))
    d.close(); gc.collect(); return out


pd_ = fitz.open(PAGE); pp = pd_[0]; assert abs(pp.rect.width - W) < 0.05 and abs(pp.rect.height - (H13 + DELTA)) < 0.05, (pp.rect, W, H13, DELTA)
PSP = spans_flat(PAGE); assert all(Y_TITLE <= s["bbox"][1] and s["bbox"][3] <= Y_D_NEW for s in PSP), "rebuilt page has text outside [16, 772 + DELTA)"
bd = fitz.open(BAND); bp = bd[0]; bw, bh = bp.rect.width, bp.rect.height
assert abs(bw - W) < 0.05, ("band width must equal the sheet width", bw, W)
BSP = spans_flat(BAND)
BAND_TOP = round(Y_D_NEW + GAP_BAND, 3); H_NEW = round(BAND_TOP + bh, 3)
d_top = round(BAND_TOP + D_OFF, 3); d_base = round(d_top + 11.75, 3)
assert d_top + L.LETTER_SIZE * 1.117 < BAND_TOP, ("letter d must sit above the band", d_top, BAND_TOP)   # 1.117 = Arial Bold ascent 0.905 + descent 0.212 (13 pt: 14.52)
out = fitz.open(); pg = out.new_page(width=W, height=H_NEW)
pg.show_pdf_page(fitz.Rect(0, Y_TITLE, W, Y_D_NEW), pd_, 0, clip=fitz.Rect(0, Y_TITLE, W, Y_D_NEW))
pg.show_pdf_page(fitz.Rect(0, BAND_TOP, W, BAND_TOP + bh), bd, 0)
font = fitz.Font(fontfile=L.ARIALB)
tw = fitz.TextWriter(pg.rect); stamped = []
for ch in ("a", "b", "c"):
    x, top = LB[ch][1], LB[ch][2] + (C_SHIFT if ch == "c" else 0.0); tw.append(fitz.Point(x, top + 11.75), ch, font=font, fontsize=L.LETTER_SIZE); stamped.append((ch, x, round(top + 11.75, 3)))
tw.append(fitz.Point(LB["d"][1], d_base), "d", font=font, fontsize=L.LETTER_SIZE); stamped.append(("d", LB["d"][1], d_base))
tw.write_text(pg, color=(0, 0, 0))
tt = fitz.TextWriter(pg.rect); TITLE_X = LB["a"][1]; tt.append(fitz.Point(TITLE_X, L.TITLE_BASELINE), "Figure 3", font=font, fontsize=L.TITLE_SIZE); tt.write_text(pg, color=L.INK_RGB)
out.save(OUT, garbage=3, deflate=True); out.close(); bd.close(); pd_.close(); gc.collect()
n, w, h = L.gs_pdfinfo(OUT); assert n == 1 and abs(w - W) < 0.01 and abs(h - H_NEW) < 0.01, (n, w, h, H_NEW)
log = {"sheet": SHEET, "page_abc": PAGE, "page_abc_sha256": L.sha256(PAGE), "band": BAND, "band_sha256": L.sha256(BAND), "band_page": [round(bw, 3), round(bh, 3)], "band_text": [s["text"] for s in BSP],
       "delta_pt": DELTA, "panel_c_shift_pt": C_SHIFT, "placements": [{"what": "rebuilt page (a, b, c; c shifted)", "rect": [0, Y_TITLE, W, Y_D_NEW]}, {"what": "band d (round 40)", "rect": [0, BAND_TOP, W, H_NEW]}],
       "band_top": BAND_TOP, "v13_gap_panels_to_band": GAP_BAND, "letter_d_offset_from_band_top": D_OFF, "letters_stamped": stamped, "letters_v13_cache": [list(LB[ch]) for ch in "abcd"],
       "title": {"text": "Figure 3", "origin": [TITLE_X, L.TITLE_BASELINE], "size": L.TITLE_SIZE, "font": "Arial Bold (TextWriter)", "color": "#1a1d21"}, "letter_size": L.LETTER_SIZE, "round": 49,
       "page": [w, h], "page_v13": [W, H13], "bytes": os.path.getsize(OUT), "output": OUT, "output_sha256": L.sha256(OUT), "written": L.now()}
json.dump(log, open(f"{WORK}/compose_log{'' if OUT.endswith('_vector.pdf') else '_' + os.path.basename(OUT)[:-4]}.json", "w"), indent=1)
print(f"wrote {OUT}: page {w} x {h} (gs PDFINFO), band {os.path.basename(BAND)} {bw:.2f} x {bh:.3f} at top {BAND_TOP}, letters {stamped} at {L.LETTER_SIZE:g} pt, title at ({TITLE_X}, {L.TITLE_BASELINE}) {L.TITLE_SIZE:g} pt")
