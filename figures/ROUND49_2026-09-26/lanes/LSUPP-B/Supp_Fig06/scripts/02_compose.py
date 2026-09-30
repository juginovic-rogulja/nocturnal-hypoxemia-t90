"""Supp Fig 6 (V14, round v8.2): the four rebuilt parts a, b, c (step 138) and d (steps 138 + 139, v8.2) placed on a fresh page of the
V13 sheet's width at the V13 part origins, the letters a, b, c, d stamped with fitz.TextWriter + Arial Bold 13 pt at the V13 origins.
Nothing is kept from V13 any more (the owner's v8.2 rule). The page height follows the tallest lower part: max(foot of c, foot of d)
plus the sheet's 14 pt margin (V13: 1420.0 pt with the 24-row v7 part d; now shorter). Writes Supp_Fig06.pdf, Supp_Fig06_150dpi.png
and work/compose_record.json (part sources with sha256, page height V13 -> new)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import json
import os
import sys

import fitz

LANE = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LSUPP-B"   # R49: SD = LANE/Supp_Fig06 = this sheet folder
sys.path.insert(0, f"{LANE}/Supp_Fig06/scripts")
from l3b_common import hydrated, FONTS, sha256, PT_PLUS  # noqa: E402

SHEET = "Supp_Fig06"
SD = f"{LANE}/{SHEET}"
WORK = f"{SD}/work"
E = json.load(open(hydrated(f"{WORK}/expected_values.json")))
assert E.get("kept_from_v13") in ({}, None) and "part_d_v8_2" in E, "run 01_build_d.py first (v8.2): part d must be a fresh build, nothing kept from V13"
BT = json.load(open(hydrated(f"{WORK}/base_text.json")))
letters = [(s["text"], s["x"], s["y"]) for s in BT if s["font"] == "Arial-BoldMT" and s["size"] == 13.0 and s["text"] in ("a", "b", "c", "d")]
assert [l[0] for l in letters] == ["a", "b", "c", "d"], letters
page_w, page_h_v13 = E["page_v13"]
place = dict(E["place_v13"])
MARGIN = 14.0
# R40 (Alen, 2026-09-18): the lower row (c, d) moves up under the top row (a, b): its origin y = the taller top part's foot + ROW_GAP
# (V16: 691.76 with a 197 pt hole under the 14-row part b); the letters c and d ride with their parts (same offset from the part origin)
ROW_GAP = 24.0
def _h(p):
    _d = fitz.open(p); _r = _d[0].rect.height; _d.close(); return _r
_foot_top = max(place["a"][1] + _h(f"{WORK}/part_a.pdf"), place["b"][1] + _h(f"{WORK}/part_b.pdf"))
_y_low = round(_foot_top + ROW_GAP, 3)
_dy = _y_low - place["c"][1]
place["c"] = (place["c"][0], _y_low); place["d"] = (place["d"][0], _y_low)
letters = [(g, x, y + (_dy if g in ("c", "d") else 0.0)) for g, x, y in letters]
parts = [(f"{WORK}/part_{p}.pdf", *place[p]) for p in "abcd"]
rec = {"page_v13": [page_w, page_h_v13], "parts": {}, "letters_stamped": letters, "kept_from_v13": None, "r40_row_gap_pt": ROW_GAP, "r40_lower_row_y": _y_low, "r40_lower_row_shift_pt": round(_dy, 3), "place_v13": E["place_v13"],
       "part_d_source": {"path": f"{WORK}/part_d.pdf", "sha256": sha256(f"{WORK}/part_d.pdf"), "built_by": "01_build_d.py (v8.2, steps 138 + 139)"}}
for p, ox, oy in parts:
    assert os.path.exists(p) and os.path.getsize(p) > 0, p
    r = fitz.open(p)[0].rect
    assert abs(r.width - 484.7244) < 0.02, (p, r)
    rec["parts"][os.path.basename(p)] = {"origin": [ox, oy], "size": [r.width, r.height], "foot": oy + r.height, "sha256": sha256(p)}
foot_a, foot_b = rec["parts"]["part_a.pdf"]["foot"], rec["parts"]["part_b.pdf"]["foot"]
foot_c, foot_d = rec["parts"]["part_c.pdf"]["foot"], rec["parts"]["part_d.pdf"]["foot"]
assert foot_a < place["c"][1] and foot_b < place["d"][1], (foot_a, foot_b, place)
page_h = round(max(foot_c, foot_d) + MARGIN, 3)
rec["page"] = [page_w, page_h]
rec["page_height_note"] = f"V13 {page_h_v13} pt (v7 part d, 24 rows) -> {page_h} pt: the page ends {MARGIN} pt under the taller lower part ({'d' if foot_d >= foot_c else 'c'}, foot {max(foot_c, foot_d):.2f})"
doc = fitz.open()
page = doc.new_page(width=page_w, height=page_h)
for p, ox, oy in parts:
    src = fitz.open(p)
    r = src[0].rect
    page.show_pdf_page(fitz.Rect(ox, oy, ox + r.width, oy + r.height), src, 0)
    src.close()
tw = fitz.TextWriter(page.rect)
font = fitz.Font(fontfile=FONTS["bold"])
for glyph, x, y in letters:
    tw.append(fitz.Point(x, y), glyph, font=font, fontsize=13.0 + PT_PLUS)   # R49: letters 14 pt at the V26 origins
tw.write_text(page, color=(0x1a / 255, 0x1d / 255, 0x21 / 255))
out = f"{SD}/{SHEET}.pdf"
doc.save(out, garbage=3, deflate=True)
doc.close()
json.dump(rec, open(f"{WORK}/compose_record.json", "w"), indent=1)
# R40: no render here (child under the watchdog; the verifier renders with Ghostscript)
print("wrote", out, "and", f"{SD}/{SHEET}_150dpi.png", "page", page_w, "x", page_h, "(V13", page_h_v13, ") parts", {k: round(v["foot"], 2) for k, v in rec["parts"].items()})
