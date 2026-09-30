#!$T90_PY
"""Lane LSUPP-D child (under the watchdog, FLAT single-page files only): the text layer as PyMuPDF sees it (rawdict spans with
font, size, colour, line direction, origin, bbox, chars) and the bounding boxes of the vector drawing items.
usage: fitz_dump.py in.pdf out.json [--no-drawings]"""
import json, sys, fitz
src, out = sys.argv[1], sys.argv[2]; want_dr = "--no-drawings" not in sys.argv
d = fitz.open(src); assert len(d) == 1, len(d); p = d[0]
rec = {"src": src, "page": [p.rect.width, p.rect.height], "xobjects": len(p.get_xobjects()), "fonts": [[f[0], f[3], f[4]] for f in p.get_fonts(full=True)], "spans": [], "drawings": []}
for b in p.get_text("rawdict")["blocks"]:
    if b["type"] != 0:
        continue
    for l in b["lines"]:
        for s in l["spans"]:
            txt = "".join(ch["c"] for ch in s["chars"])
            if not txt.strip():
                continue
            rec["spans"].append({"text": txt, "font": s["font"], "size": round(s["size"], 4), "flags": s["flags"], "color": "#%06x" % s["color"], "dir": [round(v, 4) for v in l["dir"]],
                                 "bbox": [round(v, 3) for v in s["bbox"]], "origin": [round(s["origin"][0], 3), round(s["origin"][1], 3)], "asc": round(s["ascender"], 4), "desc": round(s["descender"], 4),
                                 "chars": [[ch["c"], round(ch["origin"][0], 3), round(ch["origin"][1], 3)] for ch in s["chars"]]})
if want_dr:
    for it in p.get_drawings():
        r = it["rect"]
        rec["drawings"].append({"rect": [round(r.x0, 3), round(r.y0, 3), round(r.x1, 3), round(r.y1, 3)], "type": it.get("type"), "fill": ("#%02x%02x%02x" % tuple(int(round(v * 255)) for v in it["fill"])) if it.get("fill") is not None else None,
                                "color": ("#%02x%02x%02x" % tuple(int(round(v * 255)) for v in it["color"])) if it.get("color") is not None else None, "width": it.get("width"), "n_items": len(it.get("items", []))})
json.dump(rec, open(out, "w"), indent=0)
print("spans", len(rec["spans"]), "drawings", len(rec["drawings"]), "xobjects", rec["xobjects"], "fonts", len(rec["fonts"]))
