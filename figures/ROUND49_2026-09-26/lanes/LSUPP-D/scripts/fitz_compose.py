#!$T90_PY
"""Lane LSUPP-D child: PyMuPDF for PLACEMENT only (the round-40 fitz_place.py compose mode). No get_text / get_pixmap / get_drawings.
usage: fitz_compose.py job.json   job = {out, page: [w, h], parts: [{pdf, rect, clip}], metadata}"""
import json, sys, fitz
job = json.load(open(sys.argv[1]))
doc = fitz.open(); page = doc.new_page(width=job["page"][0], height=job["page"][1])
for prt in job["parts"]:
    src = fitz.open(prt["pdf"]); r = fitz.Rect(*prt["rect"]); clip = fitz.Rect(*prt["clip"]) if prt.get("clip") else None
    page.show_pdf_page(r, src, prt.get("page", 0), clip=clip); src.close()
if job.get("metadata"):
    doc.set_metadata(job["metadata"])
doc.save(job["out"], garbage=3, deflate=True); doc.close()
print("compose ok", job["out"])
