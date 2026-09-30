#!$T90_PY
"""Lane LSUPP round 40: PyMuPDF for PLACEMENT only, run as a child process under wd_run.sh.
usage: fitz_place.py job.json
job: {"mode": "compose", "out": pdf, "page": [w, h], "parts": [{"pdf": path, "rect": [x0, y0, x1, y1], "clip": [..] | null}],
      "letters": [{"text": "a", "x": 14.0, "y": 22.99, "size": 13.0, "bold": true}], "metadata": {...}}
     {"mode": "box", "src": pdf, "out": pdf, "box": [x0, y0, x1, y1]}            set MediaBox = CropBox of page 1
     {"mode": "redact", "src": pdf, "out": pdf, "rects": [[x0, y0, x1, y1]], "keep_graphics": true}   text-only redaction
No get_text / get_pixmap / get_drawings anywhere in this file."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import json, sys
import fitz

FONTS = {"regular": paths.FONT_ARIAL, "bold": paths.FONT_ARIAL_BOLD}
job = json.load(open(sys.argv[1]))
mode = job["mode"]
if mode == "compose":
    doc = fitz.open(); page = doc.new_page(width=job["page"][0], height=job["page"][1])
    for p in job["parts"]:
        src = fitz.open(p["pdf"]); r = fitz.Rect(*p["rect"]); clip = fitz.Rect(*p["clip"]) if p.get("clip") else None
        page.show_pdf_page(r, src, p.get("page", 0), clip=clip); src.close()
    if job.get("letters"):
        tw = fitz.TextWriter(page.rect); fb = fitz.Font(fontfile=FONTS["bold"]); fr = fitz.Font(fontfile=FONTS["regular"])
        for l in job["letters"]:
            tw.append(fitz.Point(l["x"], l["y"]), l["text"], font=fb if l.get("bold", True) else fr, fontsize=l.get("size", 13.0))
        tw.write_text(page, color=(0x1a / 255, 0x1d / 255, 0x21 / 255))
    if job.get("metadata"): doc.set_metadata(job["metadata"])
    doc.save(job["out"], garbage=3, deflate=True); doc.close()
elif mode == "box":
    doc = fitz.open(job["src"]); assert len(doc) == 1, len(doc); pg = doc[0]
    r = fitz.Rect(*job["box"]); cur = pg.rect
    # snap to the page's own edges when the requested edge is within 0.01 pt (a rounded 233.28 must not exceed a stored 233.2799988)
    if abs(r.x1 - cur.width) < 0.01: r.x1 = cur.width
    if abs(r.y1 - cur.height) < 0.01: r.y1 = cur.height
    pg.set_mediabox(r); pg.set_cropbox(pg.mediabox)   # the cropbox = the stored mediabox (its float32 width rejects the rounded rect)
    if job.get("metadata"): doc.set_metadata(job["metadata"])
    doc.save(job["out"], garbage=3, deflate=True); doc.close()
elif mode == "redact":
    doc = fitz.open(job["src"]); pg = doc[0]
    for rect in job["rects"]: pg.add_redact_annot(fitz.Rect(*rect), fill=False)
    ok = pg.apply_redactions(images=0, graphics=0, text=0)   # text only: images and vector graphics untouched
    assert ok, "apply_redactions returned False"
    doc.save(job["out"], garbage=3, deflate=True); doc.close()
else:
    raise SystemExit(f"unknown mode {mode}")
print("fitz_place ok", mode, job["out"])
