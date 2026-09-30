#!$T90_PY
"""Round 42, LFIG5b (the P-column variant; a copy of the LFIG5 composer, the layout unchanged), Main_Fig5 compose (vector): the V16 geometry with the round-40 full-width four-organ band (decision 4).
  a = the V13 sheet clipped to (0, 16, 497, 393) (panel a with its letter; the V13 title strip excluded, the title is re-stamped);
      the V13 page is nested (454 XObjects): placed by show_pdf_page only, never read with get_text/get_drawings;
  b = the flat page work/panels_bc_flat.pdf clipped to (497, 0, 969, 725) (18 + 3 rows, the key on the left, q columns; letter b on the page);
  c = the flat page clipped to (0, 393, 497, 725) (letter c on the page);
  d = band_5d_r40.pdf (968.94 x 198.425, ROUND40 LS_SKETCH, sha256 eff4d4b8...) full width at the old band slot: the V16 band top
      946.2253 - 198.425 = 747.800 (the round-34 composer located the V13 band at page height minus band height; the V16 strip that held
      the letter d began at 725), page height = band top + band height + the V16 foot gap (0) = 946.2253 (the V16 page); letter d at its
      V13 x (36.0) and the band top + the V13 offset (the V13 letter baseline 744.76 - the V13 band top 747.800 = -3.040);
  title 'Figure 5' Arial Bold 13 pt at (36.0, 14.0) by TextWriter (round 40).
--band-top N places the band top at N pt instead (the letter d follows with the same offset, the page = N + band height); the default is the
V16 slot, recorded as this lane's decision. usage: 06_compose_r42.py [--band PATH] [--out PATH] [--band-top N]"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import argparse, gc, json, os, sys
import fitz
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND42_2026-09-19/lanes/LFIG5b/scripts")
import lfig5b_lib as L

SHEET = "Main_Fig5"; SD = f"{L.LANE}/{SHEET}"; WORK = f"{SD}/work"
ap = argparse.ArgumentParser(); ap.add_argument("--band", default=f"{L.SKETCH}/band_5d_r40.pdf"); ap.add_argument("--out", default=f"{SD}/{SHEET}_vector.pdf"); ap.add_argument("--band-top", type=float, default=None); a = ap.parse_args()
BAND = L.hydrated(a.band); OUT = a.out; BAND_SHA = "5dd157bf999fe50f"; assert L.sha256(BAND).startswith(BAND_SHA), ("band_5d_r44.pdf expected", L.sha256(BAND)[:16])
BASE = L.hydrated(f"{L.V13}/{SHEET}.pdf"); BASE_SHA = "f741b08242ee438a"; assert L.sha256(BASE).startswith(BASE_SHA), "V13 Main_Fig5 changed"
REC = json.load(open(f"{WORK}/flat_build_record.json")); FLAT = L.hydrated(REC["flat"]); assert L.sha256(FLAT) == REC["flat_sha256"], "flat page changed since the build record"
PW, PH13 = REC["page_v13"]; assert abs(PW - 968.94) < 0.01 and abs(PH13 - 946.2253) < 0.01, (PW, PH13)
LET13 = REC["letters_v13"]; assert abs(LET13["a"]["x"] - 36.0) < 0.01 and abs(LET13["d"]["x"] - 36.0) < 0.01 and abs(LET13["d"]["y"] - 744.76) < 0.01, LET13
Y_TITLE = 16.0; A_CLIP = fitz.Rect(0, Y_TITLE, 497, 393); B_CLIP = fitz.Rect(497, 0, PW, 725); C_CLIP = fitz.Rect(0, 393, 497, 725)
bd = fitz.open(BAND); bp = bd[0]; bw, bh = bp.rect.width, bp.rect.height
assert abs(bw - PW) < 0.05 and abs(bh - 198.425) < 0.05, ("band page", bw, bh)
V16_BAND_TOP = round(PH13 - bh, 4)                       # 747.8003: the round-34 composer located the V13 band at H - Hb
D_OFF = round(LET13["d"]["y"] - V16_BAND_TOP, 4)          # -3.0403: the V13 letter-d baseline against the V13 band top
BAND_TOP = round(a.band_top, 4) if a.band_top is not None else V16_BAND_TOP; FOOT_GAP = 0.0
PH_NEW = round(BAND_TOP + bh + FOOT_GAP, 4)
d_letter = (LET13["d"]["x"], round(BAND_TOP + D_OFF, 3))


def spans_flat(pdf):
    d = fitz.open(pdf); out = []
    for b in d[0].get_text("rawdict")["blocks"]:
        if b["type"] != 0: continue
        for l in b["lines"]:
            for s in l["spans"]:
                t = "".join(c["c"] for c in s["chars"]).replace("\xa0", " ")
                if t.strip(): out.append(dict(text=t, bbox=[round(v, 3) for v in s["bbox"]], origin=[round(s["chars"][0]["origin"][0], 3), round(s["chars"][0]["origin"][1], 3)], size=round(s["size"], 2)))
    d.close(); gc.collect(); return out


FSP = spans_flat(FLAT)
flat_letters = {s["text"]: s for s in FSP if len(s["text"]) == 1 and s["text"] in "bc" and s["size"] > 12.5}; assert set(flat_letters) == {"b", "c"}, flat_letters
assert all(s["bbox"][3] <= 725.0 for s in FSP), ("flat page text below the b/c clips", [s for s in FSP if s["bbox"][3] > 725.0][:3])
tb = [s for s in FSP if s["text"].strip() == "vs 1% or less on PAP (95% CI)"]; assert len(tb) == 1; B_TITLE_BOTTOM = tb[0]["bbox"][3]
assert B_TITLE_BOTTOM < BAND_TOP - 20 and d_letter[1] - 11.75 > B_TITLE_BOTTOM, ("panel b runs into the band or the letter d", B_TITLE_BOTTOM, BAND_TOP, d_letter)
out = fitz.open(); pg = out.new_page(width=PW, height=PH_NEW); placed = []
src = fitz.open(BASE); assert abs(src[0].rect.width - PW) < 0.01 and abs(src[0].rect.height - PH13) < 0.01
pg.show_pdf_page(A_CLIP, src, 0, clip=A_CLIP); placed.append({"what": "panel a (V13 sheet clipped; title strip excluded)", "rect": list(A_CLIP)})
flat = fitz.open(FLAT)
for nm, clip in (("b", B_CLIP), ("c", C_CLIP)):
    pg.show_pdf_page(clip, flat, 0, clip=clip); placed.append({"what": f"panel {nm} (flat page clipped, with its letter)", "rect": list(clip)})
d_rect = fitz.Rect(0, BAND_TOP, PW, BAND_TOP + bh); pg.show_pdf_page(d_rect, bd, 0); placed.append({"what": "band d (round-40 full-width four-organ band)", "rect": [round(v, 4) for v in d_rect], "scale": 1.0})
font = fitz.Font(fontfile=L.ARIALB)
tw = fitz.TextWriter(pg.rect); tw.append(fitz.Point(*d_letter), "d", font=font, fontsize=13.0); tw.write_text(pg, color=(0, 0, 0))
tt = fitz.TextWriter(pg.rect); tt.append(fitz.Point(LET13["a"]["x"], L.TITLE_BASELINE), "Figure 5", font=font, fontsize=L.TITLE_SIZE); tt.write_text(pg, color=L.INK_RGB)
out.save(OUT, garbage=1, deflate=True); out.close(); src.close(); flat.close(); bd.close(); gc.collect()
n, w, h = L.gs_pdfinfo(OUT); assert n == 1 and abs(w - PW) < 0.01 and abs(h - PH_NEW) < 0.01, (n, w, h, PH_NEW)
log = {"sheet": SHEET, "layout": "V16 geometry: a (V13 clip), b and c (flat clips), full-width band at the V16 band slot", "base": BASE, "base_sha256": L.sha256(BASE), "flat": FLAT, "flat_sha256": L.sha256(FLAT), "band": BAND, "band_sha256": L.sha256(BAND), "band_page": [round(bw, 3), round(bh, 3)],
       "band_top": BAND_TOP, "band_top_rule": "V16 band top = page height - band height (747.8003)" if a.band_top is None else f"--band-top {a.band_top}", "foot_gap": FOOT_GAP, "band_placed_rect": [round(v, 4) for v in d_rect], "band_scale": 1.0,
       "letter_d": {"x": d_letter[0], "baseline": d_letter[1], "v13_offset_from_band_top": D_OFF, "where": "V13 x, band top + the V13 offset"}, "letters_b_c_on_flat_page": {k: v["origin"] for k, v in flat_letters.items()}, "letter_a": "V13 (inside the panel a clip)",
       "panel_b_axis_title_bottom": B_TITLE_BOTTOM, "title": {"text": "Figure 5", "origin": [LET13["a"]["x"], L.TITLE_BASELINE], "size": L.TITLE_SIZE, "font": "Arial Bold (TextWriter)", "color": "#1a1d21"}, "placements": placed, "page": [w, h], "page_v13": [PW, PH13],
       "output": OUT, "output_sha256": L.sha256(OUT), "bytes": os.path.getsize(OUT), "written": L.now()}
json.dump(log, open(f"{WORK}/compose_log{'' if OUT.endswith('_vector.pdf') else '_' + os.path.basename(OUT)[:-4]}.json", "w"), indent=1)
print(f"wrote {OUT}: page {w} x {h}; band {os.path.basename(BAND)} {bw:.2f} x {bh:.3f} at top {BAND_TOP} ({log['band_top_rule']}); letter d at {d_letter}; b axis-title bottom {B_TITLE_BOTTOM:.2f}; title at ({LET13['a']['x']}, {L.TITLE_BASELINE}) 13 pt")
