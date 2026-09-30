#!$T90_PY
"""Round 54, lane L5 (2026-09-29): Main_Fig5 compose (vector), the round-49 compose (06_compose_r49.py, the round-44 layout) at the round-54
sizes with the re-laid panels and a page height that follows the content and the band:
  a = the panel a flow page rebuilt by 04_build_panel_a_r54.py (block and arm strings 14 pt, captions 15 pt on the sheet; drawings identical
      to the Aug-25 panel), placed by the similarity fitted on the V13 sheet (work/panel_a/map_v13.json) and clipped to (0, 16, 497, CLIP_AC)
      where CLIP_AC = 385 (V13 393) comes from the flat build record; the letter a re-stamped at 18 pt at its V13 origin (36.0, 41.76);
  b = the flat page work/panels_bc_flat.pdf (05_build_flat_r54.py) clipped to (497, 0, PW, CONTENT_BOTTOM + 2) (letter b on the page, 18 pt);
  c = the flat page clipped to (0, CLIP_AC, 497, CONTENT_BOTTOM + 2) (letter c on the page, 18 pt);
  d = the band (--band, default the round-54 LSKETCH band when it exists, else the round-49 band as a stand-in) full width, its page box
      read from the file: band top = the flat content bottom (panel b's second axis title line, descender included) + BAND_AIR 20 (the
      round-44 rule: the band top at least 20 pt under panel b's axis title), rounded up to 0.1 pt; page height = band top + band height;
      the letter d at its V13 x (36.0) and band top + the V13 offset (-3.040), 18 pt;
  title 'Figure 5' Arial Bold 18 pt at (36.0, 18.0) by TextWriter (round 40 idiom; the baseline 14 -> 18 keeps the cap top 5 pt inside).
Page 968.94 x (band top + band height), asserted at most 1110 (the coordinator's aim) and 1260 (the hard gate).
usage: 06_compose_r54.py [--band PATH] [--out PATH]"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import argparse, gc, json, math, os, sys
import fitz
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND54_2026-09-28/lanes/L5/scripts")
import l5_lib as L

SHEET = "Main_Fig5"; SD = f"{L.LANE}/{SHEET}"; WORK = f"{SD}/work"
ap = argparse.ArgumentParser(); ap.add_argument("--band", default=None); ap.add_argument("--out", default=f"{WORK}/{SHEET}_r54_vector.pdf"); a = ap.parse_args()
BAND = a.band or (f"{L.SKETCH}/band_5d_r54.pdf" if os.path.exists(f"{L.SKETCH}/band_5d_r54.pdf") else f"{L.SKETCH_R49}/band_5d_r49.pdf")
BAND = L.hydrated(BAND); OUT = a.out; BAND_IS_R54 = os.path.basename(BAND) == "band_5d_r54.pdf"
REC = json.load(open(f"{WORK}/flat_build_record.json")); FLAT = L.hydrated(REC["flat"]); assert L.sha256(FLAT) == REC["flat_sha256"], "flat page changed since the build record"
assert REC["round"] == 54 and REC["sizes"]["label"] == 14.0 and REC["sizes"]["letter"] == 18.0
PW, PH13 = REC["page_v13"]; assert abs(PW - 968.94) < 0.01 and abs(PH13 - 946.2253) < 0.01, (PW, PH13)
LET13 = REC["letters_v13"]; assert abs(LET13["a"]["x"] - 36.0) < 0.01 and abs(LET13["a"]["y"] - 41.76) < 0.01 and abs(LET13["d"]["x"] - 36.0) < 0.01 and abs(LET13["d"]["y"] - 744.76) < 0.01, LET13
MAP = json.load(open(f"{WORK}/panel_a/map_v13.json")); S, TX, TY = MAP["S"], MAP["TX"], MAP["TY"]; assert MAP["max_residual"] < 0.05 and MAP["generator_reproduces_aug25"], MAP
PA = L.hydrated(f"{WORK}/panel_a/Figure5A_r54.pdf"); PA_DRAWN = json.load(open(f"{WORK}/panel_a/Figure5A_r54_flow_drawn_values.json"))
assert L.sha256(PA) == REC["panel_a_r54_page"]["sha256"], "panel a page changed since the flat build read it"
CAN = PA_DRAWN["text_pt_on_canvas"]; assert abs(CAN["blocks_arms"] * S - 14.0) < 1e-6 and abs(CAN["captions"] * S - 15.0) < 1e-6, ("panel a sizes on the sheet are not 14 / 15", CAN, S)
Y_TITLE = 16.0; CLIP_AC = REC["clip_ac"]; CONTENT_BOTTOM = REC["content_bottom"]; BAND_AIR = 20.0
BC_BOTTOM = math.ceil((CONTENT_BOTTOM + 2.0) * 10) / 10
A_CLIP = fitz.Rect(0, Y_TITLE, 497, CLIP_AC); B_CLIP = fitz.Rect(497, 0, PW, BC_BOTTOM); C_CLIP = fitz.Rect(0, CLIP_AC, 497, BC_BOTTOM)
bd = fitz.open(BAND); bp = bd[0]; bw, bh = bp.rect.width, bp.rect.height
assert abs(bw - PW) < 0.5, ("band page width", bw, PW)
V16_D_OFF = round(LET13["d"]["y"] - (PH13 - 198.4250704), 4)          # the V13 letter d sat 3.040 pt above the V16 band top (band height 198.425)
BAND_TOP = math.ceil((CONTENT_BOTTOM + BAND_AIR) * 10) / 10
PH_NEW = round(BAND_TOP + bh, 4)
assert PH_NEW <= L.HEIGHT_HARD, ("page taller than the hard gate", PH_NEW);
d_letter = (LET13["d"]["x"], round(BAND_TOP + V16_D_OFF, 3))
LETTER_PT, TITLE_PT = L.LETTER_SIZE, L.TITLE_SIZE; assert LETTER_PT == 18.0 and TITLE_PT == 18.0


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
flat_letters = {s["text"]: s for s in FSP if len(s["text"]) == 1 and s["text"] in "bc" and s["size"] > 16}; assert set(flat_letters) == {"b", "c"} and all(abs(v["size"] - LETTER_PT) < 0.05 for v in flat_letters.values()), flat_letters
assert all(s["bbox"][3] <= BC_BOTTOM for s in FSP), ("flat page text below the b/c clips", [s for s in FSP if s["bbox"][3] > BC_BOTTOM][:3])
assert all(s["bbox"][1] >= CLIP_AC for s in FSP if s["bbox"][0] < 497), ("panel c text above the a/c clip boundary", [s for s in FSP if s["bbox"][0] < 497 and s["bbox"][1] < CLIP_AC][:3])
assert all(s["bbox"][0] >= 497 for s in FSP if s["bbox"][1] < CLIP_AC), ("panel b text left of its clip", [s for s in FSP if s["bbox"][1] < CLIP_AC and s["bbox"][0] < 497][:3])
tb = [s for s in FSP if s["text"].strip() == "vs 1% or less on PAP (95% CI)"]; assert len(tb) == 1; B_TITLE_BOTTOM = tb[0]["bbox"][3]
assert B_TITLE_BOTTOM <= BAND_TOP - BAND_AIR + 0.5 and d_letter[1] - 0.905 * LETTER_PT > B_TITLE_BOTTOM, ("panel b runs into the band or the letter d", B_TITLE_BOTTOM, BAND_TOP, d_letter)
# panel a: the generator page under the V13 similarity, clipped to the panel a clip
pa_doc = fitz.open(PA); pa_rect = pa_doc[0].rect; assert len(pa_doc) == 1 and len(pa_doc[0].get_xobjects()) == 0
src_clip = fitz.Rect(max(pa_rect.x0, (A_CLIP.x0 - TX) / S), max(pa_rect.y0, (A_CLIP.y0 - TY) / S), min(pa_rect.x1, (A_CLIP.x1 - TX) / S), min(pa_rect.y1, (A_CLIP.y1 - TY) / S))
tgt = fitz.Rect(S * src_clip.x0 + TX, S * src_clip.y0 + TY, S * src_clip.x1 + TX, S * src_clip.y1 + TY)
assert A_CLIP.contains(tgt) and abs(tgt.width / src_clip.width - S) < 1e-6 and abs(tgt.height / src_clip.height - S) < 1e-6, (tgt, src_clip)
PSP = spans_flat(PA); pa_sheet = [dict(text=s["text"], x=round(S * s["origin"][0] + TX, 3), y=round(S * s["origin"][1] + TY, 3), size=round(S * s["size"], 3), bbox=[round(S * s["bbox"][0] + TX, 3), round(S * s["bbox"][1] + TY, 3), round(S * s["bbox"][2] + TX, 3), round(S * s["bbox"][3] + TY, 3)]) for s in PSP]
assert all(A_CLIP.x0 < s["bbox"][0] and s["bbox"][2] < A_CLIP.x1 - 4 and A_CLIP.y0 + 4 < s["bbox"][1] and s["bbox"][3] < A_CLIP.y1 - 4 for s in pa_sheet), "panel a text outside the clip"
out = fitz.open(); pg = out.new_page(width=PW, height=PH_NEW); placed = []
pg.show_pdf_page(tgt, pa_doc, 0, clip=src_clip); placed.append({"what": "panel a (round-54 flow page, 14/15 pt on the sheet, under the V13 similarity, clipped to the panel a clip)", "rect": [round(v, 4) for v in tgt], "clip": list(A_CLIP), "src_clip": [round(v, 4) for v in src_clip], "scale": S})
flat = fitz.open(FLAT)
for nm, clip in (("b", B_CLIP), ("c", C_CLIP)):
    pg.show_pdf_page(clip, flat, 0, clip=clip); placed.append({"what": f"panel {nm} (flat page clipped, with its letter at {LETTER_PT} pt)", "rect": list(clip)})
d_rect = fitz.Rect(0, BAND_TOP, PW, BAND_TOP + bh); pg.show_pdf_page(d_rect, bd, 0); placed.append({"what": f"band d ({os.path.basename(BAND)}, full width, page box read from the file)", "rect": [round(v, 4) for v in d_rect], "scale": 1.0})
font = fitz.Font(fontfile=L.ARIALB); ink_a = tuple(int(LET13["a"]["color"][i:i + 2], 16) / 255 for i in (0, 2, 4))
ta = fitz.TextWriter(pg.rect); ta.append(fitz.Point(LET13["a"]["x"], LET13["a"]["y"]), "a", font=font, fontsize=LETTER_PT); ta.write_text(pg, color=ink_a)
tw = fitz.TextWriter(pg.rect); tw.append(fitz.Point(*d_letter), "d", font=font, fontsize=LETTER_PT); tw.write_text(pg, color=(0, 0, 0))
tt = fitz.TextWriter(pg.rect); tt.append(fitz.Point(LET13["a"]["x"], L.TITLE_BASELINE), "Figure 5", font=font, fontsize=TITLE_PT); tt.write_text(pg, color=L.INK_RGB)
out.save(OUT, garbage=1, deflate=True); out.close(); pa_doc.close(); flat.close(); bd.close(); gc.collect()
n, w, h = L.gs_pdfinfo(OUT); assert n == 1 and abs(w - PW) < 0.01 and abs(h - PH_NEW) < 0.01, (n, w, h, PH_NEW)
log = {"sheet": SHEET, "round": 54, "layout": "round-44 geometry re-laid at 14/15/18 pt: a (the rebuilt flow page under the V13 similarity), b and c (flat clips), the band at the content bottom + 20 pt; page height = band top + band height",
       "panel_a": {"page": PA, "sha256": L.sha256(PA), "map": MAP, "src_clip": [round(v, 4) for v in src_clip], "target": [round(v, 4) for v in tgt], "text_on_sheet": pa_sheet, "canvas_sizes": CAN, "sheet_sizes": {"blocks_arms": round(CAN["blocks_arms"] * S, 4), "captions": round(CAN["captions"] * S, 4)}, "x_arm": PA_DRAWN.get("x_arm"), "drawn_values": PA_DRAWN},
       "flat": FLAT, "flat_sha256": L.sha256(FLAT), "band": BAND, "band_sha256": L.sha256(BAND), "band_is_r54": BAND_IS_R54, "band_page": [round(bw, 3), round(bh, 3)],
       "band_top": BAND_TOP, "band_top_rule": f"flat content bottom {CONTENT_BOTTOM:.3f} (panel b's second axis title line, descender included) + {BAND_AIR} pt, rounded up to 0.1", "band_air": BAND_AIR, "band_placed_rect": [round(v, 4) for v in d_rect], "band_scale": 1.0,
       "clips": {"a": list(A_CLIP), "b": list(B_CLIP), "c": list(C_CLIP), "clip_ac": CLIP_AC, "bc_bottom": BC_BOTTOM, "v13": {"a": [0, 16, 497, 393], "b": [497, 0, PW, 725], "c": [0, 393, 497, 725]}},
       "letter_a": {"x": LET13["a"]["x"], "baseline": LET13["a"]["y"], "size": LETTER_PT, "color": LET13["a"]["color"], "where": "V13 origin, re-stamped (TextWriter, Arial Bold)"},
       "letter_d": {"x": d_letter[0], "baseline": d_letter[1], "size": LETTER_PT, "v16_offset_from_band_top": V16_D_OFF, "where": "V13 x, band top + the V16 offset"}, "letters_b_c_on_flat_page": {k: dict(origin=v["origin"], size=v["size"]) for k, v in flat_letters.items()},
       "panel_b_axis_title_bottom": B_TITLE_BOTTOM, "title": {"text": "Figure 5", "origin": [LET13["a"]["x"], L.TITLE_BASELINE], "size": TITLE_PT, "font": "Arial Bold (TextWriter)", "color": "#1a1d21", "v13_origin": [36.0, 14.0], "note": "baseline 14 -> 18 at 18 pt: the cap top 5.1 pt inside the page (round 49: 4.0)"}, "placements": placed, "page": [w, h], "page_v13": [PW, PH13], "height_aim": L.HEIGHT_AIM, "height_hard": L.HEIGHT_HARD,
       "output": OUT, "output_sha256": L.sha256(OUT), "bytes": os.path.getsize(OUT), "written": L.now()}
json.dump(log, open(f"{WORK}/compose_log{'' if OUT.endswith('_r54_vector.pdf') else '_' + os.path.basename(OUT)[:-4]}.json", "w"), indent=1)
print(f"wrote {OUT}: page {w} x {h} ({os.path.getsize(OUT):,} bytes; aim <= {L.HEIGHT_AIM}, hard <= {L.HEIGHT_HARD}); panel a at {[round(v, 2) for v in tgt]} scale {S:.5f}; band {os.path.basename(BAND)} {bw:.2f} x {bh:.3f} at top {BAND_TOP} (content bottom {CONTENT_BOTTOM:.2f} + {BAND_AIR}); letters a ({LET13['a']['x']}, {LET13['a']['y']}) d {d_letter} at {LETTER_PT} pt, b/c on the flat page at {[v['size'] for v in flat_letters.values()]} pt; b axis-title bottom {B_TITLE_BOTTOM:.2f}; title at ({LET13['a']['x']}, {L.TITLE_BASELINE}) {TITLE_PT} pt")
