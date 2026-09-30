#!$T90_PY
"""Supp_Fig20 round 49 composer (PyMuPDF placement only; run as a child under the watchdog), the round-40 composer repointed:
the four rebuilt sub-sheets on a clean page of the V26 size at the round-37 solved offsets, letters a to d stamped at the V26
origins (base_text.json) in Arial Bold at 13 + 1 = 14 pt."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import json, os, sys
import fitz
L = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LSUPP-C"
HERE = f"{L}/Supp_Fig20"; NEW = f"{HERE}/work/subsheets_new"; OUT = f"{HERE}/Supp_Fig20.pdf"
ARIALB = paths.FONT_ARIAL_BOLD
PT_PLUS = 1.0
W, H = 1089.48, 1441.46
R30 = {"A": (14.0, 14.029), "B": (556.74, 13.976), "C": (14.0, 737.154), "D": (556.74, 737.26)}
STEMS = {"A": "eFigure12A_windows_track_each_other", "B": "eFigure12B_rem_vs_whole_sleep", "C": "eFigure12C_healthy_subgroup", "D": "eFigure12D_ranking_worth_and_floor"}
BT = json.load(open(f"{HERE}/work/base_text.json"))["spans"]
out = fitz.open(); page = out.new_page(width=W, height=H); placed = {}
for tag in "ABCD":
    p = f"{NEW}/{STEMS[tag]}.pdf"; d = fitz.open(p); r = d[0].rect; dx, dy = R30[tag]
    page.show_pdf_page(fitz.Rect(dx, dy, dx + r.width, dy + r.height), d, 0); d.close()
    placed[tag] = {"pdf": p, "page": [r.width, r.height], "offset": [dx, dy], "foot": dy + r.height}
assert placed["A"]["foot"] < 730 and placed["B"]["foot"] < 730, (placed["A"]["foot"], placed["B"]["foot"])
tw = fitz.TextWriter(page.rect); font = fitz.Font(fontfile=ARIALB); stamped = []
for glyph in ("a", "b", "c", "d"):
    src = [s for s in BT if s["text"] == glyph and s["size"] >= 12]; assert len(src) >= 1, glyph
    o = src[0]["origin"]; tw.append(fitz.Point(o[0], o[1]), glyph, font=font, fontsize=src[0]["size"] + PT_PLUS); stamped.append({"letter": glyph, "origin": o, "size": src[0]["size"] + PT_PLUS})
tw.write_text(page, color=(0x1a / 255, 0x1d / 255, 0x21 / 255))
out.set_metadata({"title": "Supp_Fig20", "creator": "LSUPP-C r49 01_build_subsheets.py (matplotlib, +1 pt) + r49_compose_s20.py (fitz placement)"})
out.save(OUT, garbage=3, deflate=True); out.close()
json.dump({"sheet": "Supp_Fig20", "page": [W, H], "placed": placed, "stamps": stamped, "note": "R49: the round-40 composition (sub-sheets whole at the round-37 solved offsets, nothing kept from V13), letters 14 pt at the V26 origins"}, open(f"{HERE}/work/compose_record.json", "w"), indent=1)
print("wrote", OUT, {k: round(v["foot"], 2) for k, v in placed.items()})
