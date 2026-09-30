#!$T90_PY
"""Round 49, lane L5 (2026-09-26): Main_Fig5 compose (vector), the round-44 layout (compose_fig5_r44.py) with every text one point larger.
  a = the panel a flow page rebuilt by 04_build_panel_a_r49.py (text + 1/S pt on its canvas, drawings identical to the Aug-25 panel), placed by
      the similarity fitted on the V13 sheet (work/panel_a/map_v13.json: S 1.13435, TX -29.086, TY 37.718, residual 0.006 pt) and clipped to the
      V13 panel a clip (0, 16, 497, 393); the letter a re-stamped at 14 pt at its V13 origin (36.0, 41.76) in the V13 ink (the V13 sheet is no
      longer nested into this compose);
  b = the flat page work/panels_bc_flat.pdf (05_build_flat_r49.py, +1 pt) clipped to (497, 0, 969, 725) (letter b on the page, 14 pt);
  c = the flat page clipped to (0, 393, 497, 725) (letter c on the page, 14 pt);
  d = the round-49 band (lanes/LSKETCH/work/build/band_5d_r49.pdf, 968.94 x 198.425, notes 12 pt, captions 11 pt) full width at the V16 band
      slot: band top = page height 946.2253 - band height = 747.800 (the round-44 rule), the letter d at its V13 x (36.0) and the band top + the
      V13 offset (-3.040) at 14 pt;
  title 'Figure 5' Arial Bold 14 pt at (36.0, 14.0) by TextWriter (round 40 idiom, 13 -> 14).
Page 968.94 x 946.2253 (= V26). usage: 06_compose_r49.py [--band PATH] [--out PATH] [--band-top N]"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import argparse, gc, json, os, sys
import fitz
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/L5/scripts")
import l5_lib as L

SHEET = "Main_Fig5"; SD = f"{L.LANE}/{SHEET}"; WORK = f"{SD}/work"
ap = argparse.ArgumentParser(); ap.add_argument("--band", default=f"{L.SKETCH}/band_5d_r49.pdf"); ap.add_argument("--out", default=f"{WORK}/{SHEET}_r49_vector.pdf"); ap.add_argument("--band-top", type=float, default=None); a = ap.parse_args()
BAND = L.hydrated(a.band); OUT = a.out
BAND_SHAS = {"band_5d_r49.pdf": "8f9f5568dfd6307f", "band_5d_r44.pdf": "5dd157bf999fe50f"}   # r49 band rebuilt 14:45 (fork arrowheads 4.5 mm shorter); the 14:24 build was 086f89253e7dd7e1
assert L.sha256(BAND).startswith(BAND_SHAS[os.path.basename(BAND)]), ("band changed since this lane read it", L.sha256(BAND)[:16])
REC = json.load(open(f"{WORK}/flat_build_record.json")); FLAT = L.hydrated(REC["flat"]); assert L.sha256(FLAT) == REC["flat_sha256"], "flat page changed since the build record"
PW, PH13 = REC["page_v13"]; assert abs(PW - 968.94) < 0.01 and abs(PH13 - 946.2253) < 0.01, (PW, PH13)
LET13 = REC["letters_v13"]; assert abs(LET13["a"]["x"] - 36.0) < 0.01 and abs(LET13["a"]["y"] - 41.76) < 0.01 and abs(LET13["d"]["x"] - 36.0) < 0.01 and abs(LET13["d"]["y"] - 744.76) < 0.01, LET13
MAP = json.load(open(f"{WORK}/panel_a/map_v13.json")); S, TX, TY = MAP["S"], MAP["TX"], MAP["TY"]; assert MAP["max_residual"] < 0.05 and MAP["generator_reproduces_aug25"], MAP
PA = L.hydrated(f"{WORK}/panel_a/Figure5A_r49.pdf"); PA_DRAWN = json.load(open(f"{WORK}/panel_a/Figure5A_r49_flow_drawn_values.json"))
assert abs(PA_DRAWN["bump_on_canvas_pt"] * S - 1.0) < 1e-6, ("panel a bump is not +1.0 pt on the sheet", PA_DRAWN["bump_on_canvas_pt"], S)
Y_TITLE = 16.0; A_CLIP = fitz.Rect(0, Y_TITLE, 497, 393); B_CLIP = fitz.Rect(497, 0, PW, 725); C_CLIP = fitz.Rect(0, 393, 497, 725)
bd = fitz.open(BAND); bp = bd[0]; bw, bh = bp.rect.width, bp.rect.height
assert abs(bw - PW) < 0.05 and abs(bh - 198.425) < 0.05, ("band page", bw, bh)
V16_BAND_TOP = round(PH13 - bh, 4); D_OFF = round(LET13["d"]["y"] - V16_BAND_TOP, 4)
BAND_TOP = round(a.band_top, 4) if a.band_top is not None else V16_BAND_TOP; FOOT_GAP = 0.0
PH_NEW = round(BAND_TOP + bh + FOOT_GAP, 4); assert abs(PH_NEW - PH13) < 0.001, ("page height must stay the V26 height", PH_NEW, PH13)
d_letter = (LET13["d"]["x"], round(BAND_TOP + D_OFF, 3))
LETTER_PT, TITLE_PT = L.LETTER_SIZE, L.TITLE_SIZE; assert LETTER_PT == 14.0 and TITLE_PT == 14.0


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
flat_letters = {s["text"]: s for s in FSP if len(s["text"]) == 1 and s["text"] in "bc" and s["size"] > 12.5}; assert set(flat_letters) == {"b", "c"} and all(abs(v["size"] - LETTER_PT) < 0.05 for v in flat_letters.values()), flat_letters
assert all(s["bbox"][3] <= 725.0 for s in FSP), ("flat page text below the b/c clips", [s for s in FSP if s["bbox"][3] > 725.0][:3])
tb = [s for s in FSP if s["text"].strip() == "vs 1% or less on PAP (95% CI)"]; assert len(tb) == 1; B_TITLE_BOTTOM = tb[0]["bbox"][3]
assert B_TITLE_BOTTOM < BAND_TOP - 20 and d_letter[1] - 0.905 * LETTER_PT > B_TITLE_BOTTOM, ("panel b runs into the band or the letter d", B_TITLE_BOTTOM, BAND_TOP, d_letter)
# panel a: the generator page under the V13 similarity, clipped to the V13 panel a clip
pa_doc = fitz.open(PA); pa_rect = pa_doc[0].rect; assert len(pa_doc) == 1 and len(pa_doc[0].get_xobjects()) == 0
src_clip = fitz.Rect(max(pa_rect.x0, (A_CLIP.x0 - TX) / S), max(pa_rect.y0, (A_CLIP.y0 - TY) / S), min(pa_rect.x1, (A_CLIP.x1 - TX) / S), min(pa_rect.y1, (A_CLIP.y1 - TY) / S))
tgt = fitz.Rect(S * src_clip.x0 + TX, S * src_clip.y0 + TY, S * src_clip.x1 + TX, S * src_clip.y1 + TY)
assert A_CLIP.contains(tgt) and abs(tgt.width / src_clip.width - S) < 1e-6 and abs(tgt.height / src_clip.height - S) < 1e-6, (tgt, src_clip)
PSP = spans_flat(PA); pa_sheet = [dict(text=s["text"], x=round(S * s["origin"][0] + TX, 3), y=round(S * s["origin"][1] + TY, 3), size=round(S * s["size"], 3)) for s in PSP]
assert all(A_CLIP.x0 < s["x"] and S * (sp["bbox"][2]) + TX < A_CLIP.x1 - 4 and A_CLIP.y0 + 4 < S * sp["bbox"][1] + TY and S * sp["bbox"][3] + TY < A_CLIP.y1 - 2 for s, sp in zip(pa_sheet, PSP)), "panel a text outside the clip"
out = fitz.open(); pg = out.new_page(width=PW, height=PH_NEW); placed = []
pg.show_pdf_page(tgt, pa_doc, 0, clip=src_clip); placed.append({"what": "panel a (round-49 flow page, +1 pt, under the V13 similarity, clipped to the V13 panel a clip)", "rect": [round(v, 4) for v in tgt], "src_clip": [round(v, 4) for v in src_clip], "scale": S})
flat = fitz.open(FLAT)
for nm, clip in (("b", B_CLIP), ("c", C_CLIP)):
    pg.show_pdf_page(clip, flat, 0, clip=clip); placed.append({"what": f"panel {nm} (flat page clipped, with its letter at {LETTER_PT} pt)", "rect": list(clip)})
d_rect = fitz.Rect(0, BAND_TOP, PW, BAND_TOP + bh); pg.show_pdf_page(d_rect, bd, 0); placed.append({"what": f"band d ({os.path.basename(BAND)}, full width)", "rect": [round(v, 4) for v in d_rect], "scale": 1.0})
font = fitz.Font(fontfile=L.ARIALB); ink_a = tuple(int(LET13["a"]["color"][i:i + 2], 16) / 255 for i in (0, 2, 4))
ta = fitz.TextWriter(pg.rect); ta.append(fitz.Point(LET13["a"]["x"], LET13["a"]["y"]), "a", font=font, fontsize=LETTER_PT); ta.write_text(pg, color=ink_a)
tw = fitz.TextWriter(pg.rect); tw.append(fitz.Point(*d_letter), "d", font=font, fontsize=LETTER_PT); tw.write_text(pg, color=(0, 0, 0))
tt = fitz.TextWriter(pg.rect); tt.append(fitz.Point(LET13["a"]["x"], L.TITLE_BASELINE), "Figure 5", font=font, fontsize=TITLE_PT); tt.write_text(pg, color=L.INK_RGB)
out.save(OUT, garbage=1, deflate=True); out.close(); pa_doc.close(); flat.close(); bd.close(); gc.collect()
n, w, h = L.gs_pdfinfo(OUT); assert n == 1 and abs(w - PW) < 0.01 and abs(h - PH_NEW) < 0.01, (n, w, h, PH_NEW)
log = {"sheet": SHEET, "round": 49, "layout": "round-44 geometry: a (the rebuilt flow page under the V13 similarity), b and c (flat clips), the round-49 band at the V16 band slot; every text +1 pt",
       "panel_a": {"page": PA, "sha256": L.sha256(PA), "map": MAP, "src_clip": [round(v, 4) for v in src_clip], "target": [round(v, 4) for v in tgt], "text_on_sheet": pa_sheet, "bump_on_canvas_pt": PA_DRAWN["bump_on_canvas_pt"], "drawn_values": PA_DRAWN},
       "flat": FLAT, "flat_sha256": L.sha256(FLAT), "band": BAND, "band_sha256": L.sha256(BAND), "band_page": [round(bw, 3), round(bh, 3)],
       "band_top": BAND_TOP, "band_top_rule": "V16 band top = page height - band height (747.8003)" if a.band_top is None else f"--band-top {a.band_top}", "foot_gap": FOOT_GAP, "band_placed_rect": [round(v, 4) for v in d_rect], "band_scale": 1.0,
       "letter_a": {"x": LET13["a"]["x"], "baseline": LET13["a"]["y"], "size": LETTER_PT, "color": LET13["a"]["color"], "where": "V13 origin, re-stamped (TextWriter, Arial Bold)"},
       "letter_d": {"x": d_letter[0], "baseline": d_letter[1], "size": LETTER_PT, "v13_offset_from_band_top": D_OFF, "where": "V13 x, band top + the V13 offset"}, "letters_b_c_on_flat_page": {k: dict(origin=v["origin"], size=v["size"]) for k, v in flat_letters.items()},
       "panel_b_axis_title_bottom": B_TITLE_BOTTOM, "title": {"text": "Figure 5", "origin": [LET13["a"]["x"], L.TITLE_BASELINE], "size": TITLE_PT, "font": "Arial Bold (TextWriter)", "color": "#1a1d21"}, "placements": placed, "page": [w, h], "page_v13": [PW, PH13],
       "output": OUT, "output_sha256": L.sha256(OUT), "bytes": os.path.getsize(OUT), "written": L.now()}
json.dump(log, open(f"{WORK}/compose_log{'' if OUT.endswith('_r49_vector.pdf') else '_' + os.path.basename(OUT)[:-4]}.json", "w"), indent=1)
print(f"wrote {OUT}: page {w} x {h} ({os.path.getsize(OUT):,} bytes); panel a at {[round(v, 2) for v in tgt]} scale {S:.5f}; band {os.path.basename(BAND)} {bw:.2f} x {bh:.3f} at top {BAND_TOP}; letters a ({LET13['a']['x']}, {LET13['a']['y']}) d {d_letter} at {LETTER_PT} pt, b/c on the flat page at {[v['size'] for v in flat_letters.values()]} pt; b axis-title bottom {B_TITLE_BOTTOM:.2f}; title at ({LET13['a']['x']}, {L.TITLE_BASELINE}) {TITLE_PT} pt")
