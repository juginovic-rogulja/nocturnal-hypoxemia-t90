#!$T90_PY
"""Round 54 (2026-09-29, Alen's option A): Extended Data Fig. 7 gains panel c = the apnea-by-oxygenation counts cross-tab that was Figure 4c
(the round-49 build at the Extended Data text size, ROUND49 L4 panel_c_m.pdf, 11 and 12 pt, unchanged), placed under panel a in the left
column of V30's ED_Fig07 sheet, the letter c stamped as the letters a and b are (fitz.TextWriter, Arial Bold 14). V30's page is placed
whole and untouched at the top. PyMuPDF for placement only (show_pdf_page on flat and V30 pages), censuses and renders by Ghostscript."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import json, os, sys, hashlib, subprocess
import fitz
T90 = paths.FIGURE_ROOT; R54 = f"{paths.FIGURE_ROOT}/ROUND54_2026-09-28"; LANE = f"{R54}/lanes/LED7C"; SD = f"{LANE}/ED_Fig07"; W_ = f"{SD}/work"
V30 = f"{paths.FIGURE_ROOT}/ROUND52_2026-09-26/figures/NEW_FINAL_SET_V30/ED_Fig07.pdf"; PC = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/L4/Main_Fig4/work/panel_c_m.pdf"
ARIAL_BOLD = paths.FONT_ARIAL_BOLD; INK = (0x1a / 255, 0x1d / 255, 0x21 / 255)
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
assert sha(PC).startswith("274aa860b216adc4"), sha(PC)          # the round-49 panel c (compose_log of Main_Fig4 r49)
v = fitz.open(V30); pc = fitz.open(PC); W, H30 = v[0].rect.width, v[0].rect.height; wc, hc = pc[0].rect.width, pc[0].rect.height
assert abs(W - 894.61) < 0.05 and abs(H30 - 539.98) < 0.05 and abs(wc - 497.2) < 0.05 and abs(hc - 447.1) < 0.05, (W, H30, wc, hc)
Y0 = 372.0; CLIP_H = 352.0; LETTER = ("c", 12.96, Y0 + 13.0); H = round(Y0 + CLIP_H + 8.0, 2)   # panel a's ink ends at y 363 on V30's page, c's first ink (title baseline 28 inside its page) lands at 400, its last baseline (341) at 713
doc = fitz.open(); page = doc.new_page(width=W, height=H)
page.show_pdf_page(fitz.Rect(0, 0, W, H30), v, 0)
page.show_pdf_page(fitz.Rect(0, Y0, wc, Y0 + CLIP_H), pc, 0, clip=fitz.Rect(0, 0, wc, CLIP_H))
tw = fitz.TextWriter(page.rect); tw.append(fitz.Point(LETTER[1], LETTER[2]), LETTER[0], font=fitz.Font(fontfile=ARIAL_BOLD), fontsize=14.0); tw.write_text(page, color=INK)
doc.set_metadata({"title": "ED_Fig07", "creator": "LED7C r54 compose_ed7c.py: V30 ED_Fig07 + round-49 Figure 4 panel c under panel a, letter c"})
os.makedirs(W_, exist_ok=True); out = f"{SD}/ED_Fig07.pdf"; doc.save(out, garbage=3, deflate=True); doc.close(); v.close(); pc.close()
json.dump(dict(v30=V30, v30_sha256=sha(V30), panel_c=PC, panel_c_sha256=sha(PC), page=[W, H], v30_page=[W, H30], panel_c_rect=[0, Y0, wc, Y0 + CLIP_H], panel_c_clip=[0, 0, wc, CLIP_H], letter=LETTER, letter_size=14.0, out=out, out_sha256=sha(out)), open(f"{W_}/compose_log.json", "w"), indent=1)
print("wrote", out, W, "x", H)
