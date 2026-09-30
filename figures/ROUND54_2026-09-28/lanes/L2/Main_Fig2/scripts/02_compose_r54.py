#!$T90_PY
"""Round 54 (2026-09-29), lane L2, Main_Fig2 step 2: the clean sheet at the round-54 sizes. The round-49 composer (byte copy beside as
02_compose_r49_R49_ORIGINAL.py) re-run on THIS lane's panels (work/panel_{a,b,c}.pdf) with the round-54 layout from work/layout.json:
  regions: a = [0, split_x] x [0, page_h], b = [split_x, page_w] x [0, bc_y], c = [split_x, page_w] x [bc_y, page_h] (an L-shape, panel a
  one column on the left, b above c on the right);
  the panel letters a, b, c: Helvetica-Bold 18 pt at (14.2, 44), (split_x + 4, 44), (split_x + 4, bc_y + 14);
  the title "Figure 2": Arial Bold 18 pt at (14.2, 18) (cap top about 5 pt inside the page, as Figures 5 and 6 of round 54).
Font descriptors of the matplotlib panels aligned to the V13 sheet's own Arial descriptors from the round-37 probe cache
(work/r37_base_text.json), so the V13 sheet is never opened. Placement by PyMuPDF show_pdf_page of the flat panel pages (allowed);
the composed sheet is inspected with Ghostscript only."""
import gc, json, os, sys
import fitz
from fig2_common import *
from v14lib import align_font_metrics, sheet_font_metrics
import lib_r54 as R

OUT = os.path.join(LANE, "Main_Fig2.pdf")
LAY = layout(); PH, SPLIT, BC_Y = float(LAY["page_h"]), float(LAY["split_x"]), float(LAY["bc_y"])
REGION = {"a": fitz.Rect(0.0, 0.0, SPLIT, PH), "b": fitz.Rect(SPLIT, 0.0, PAGE_W, BC_Y), "c": fitz.Rect(SPLIT, BC_Y, PAGE_W, PH)}
LETTERS = {"a": (14.2, float(LAY["letter_y"])), "b": (SPLIT + 4.0, float(LAY["letter_y"])), "c": (SPLIT + 4.0, BC_Y + 14.0)}
LETTER_SIZE = float(LAY["letter_size"]); TITLE_XY = (14.2, float(LAY["title_y"])); TITLE_SIZE = pt(SIZE["letter"])
assert LETTER_SIZE == TITLE_SIZE == 18.0, (LETTER_SIZE, TITLE_SIZE)
BT = load_text()


def content_bbox(path):
    doc = fitz.open(path); page = doc[0]; r = fitz.Rect()
    for d in page.get_drawings(): r |= d["rect"]
    for b in page.get_text("dict")["blocks"]:
        if b["type"] == 0: r |= fitz.Rect(b["bbox"])
    doc.close(); gc.collect(); return r


def main():
    met = sheet_font_metrics(BT["spans"])          # {font: (asc, desc)} per mille, from the round-37 probe cache of the V13 sheet
    log = dict(sheet="Main_Fig2", round=54, page=[PAGE_W, PH], split_x=SPLIT, bc_y=BC_Y, metrics_used={}, panels={}, letters=LETTERS, letter_size=LETTER_SIZE,
               title=dict(text="Figure 2", origin=list(TITLE_XY), size=TITLE_SIZE, font="Arial Bold (TextWriter)", color="#1a1d21"))
    panels = {}
    for p in "abc":
        src = os.path.join(WORK, f"panel_{p}.pdf"); dst = os.path.join(WORK, f"panel_{p}_metrics.pdf")
        log["metrics_used"][p] = align_font_metrics(src, dst, met)
        bb = content_bbox(dst); reg = REGION[p]
        assert reg.contains(bb), (p, "panel content leaves its region", tuple(bb), tuple(reg))
        panels[p] = dst; log["panels"][p] = dict(file=dst, sha256=R.sha256(dst), content_bbox=[round(v, 3) for v in bb], region=[round(v, 3) for v in reg])
    doc = fitz.open(); page = doc.new_page(width=PAGE_W, height=PH)
    for p in "abc":
        pdoc = fitz.open(panels[p]); assert abs(pdoc[0].rect.height - PH) < 0.01 and abs(pdoc[0].rect.width - PAGE_W) < 0.01, (p, pdoc[0].rect)
        page.show_pdf_page(REGION[p], pdoc, 0, clip=REGION[p]); pdoc.close()   # target = clip = the region: 1:1
    for Lt, (x, base_y) in LETTERS.items():
        page.insert_text(fitz.Point(x, base_y), Lt, fontsize=LETTER_SIZE, fontname="hebo", color=(0, 0, 0))
    tw = fitz.TextWriter(page.rect, color=R.INK_RGB)
    tw.append(fitz.Point(*TITLE_XY), "Figure 2", font=fitz.Font(fontfile=R.ARIALB), fontsize=TITLE_SIZE)
    tw.write_text(page)
    doc.save(OUT, garbage=4, deflate=True); doc.close(); gc.collect()
    n, w, h = R.gs_pdfinfo(OUT)
    log.update(out=OUT, out_sha256=R.sha256(OUT), bytes=os.path.getsize(OUT), page_out_gs=[n, w, h], written=R.now())
    assert n == 1 and abs(w - PAGE_W) < 0.01 and abs(h - PH) < 0.01, (n, w, h)
    json.dump(log, open(os.path.join(WORK, "compose_log.json"), "w"), indent=1)
    print(f"wrote {OUT} ({os.path.getsize(OUT):,} bytes): page {w} x {h} (gs PDFINFO), regions a/b/c split x {SPLIT}, y {BC_Y}, letters hebo {LETTER_SIZE} at {LETTERS}, "
          f"title Arial Bold {TITLE_SIZE} pt at {TITLE_XY}; metrics {log['metrics_used']}; content bboxes {[log['panels'][p]['content_bbox'] for p in 'abc']}")


if __name__ == "__main__":
    main()
