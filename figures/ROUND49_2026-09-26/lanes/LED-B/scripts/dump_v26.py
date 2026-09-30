#!$T90_PY
"""LED-B round 49: span dump (PyMuPDF rawdict) of a V26 ED sheet, plus the drawings on ED_Fig06 (a flat sheet, 0 XObjects; run as a child under wd_run.sh).
usage: dump_v26.py <Sheet> <in.pdf> <out.json> [--drawings]"""
import fitz, json, sys, os, hashlib, gc
S, IN, OUT = sys.argv[1:4]; want_dr = "--drawings" in sys.argv
assert os.stat(IN).st_blocks > 0, ("evicted", IN)
d = fitz.open(IN); assert d.page_count == 1; p = d[0]
assert len(p.get_xobjects()) == 0, "not a flat sheet"
spans = []
for b in p.get_text("rawdict")["blocks"]:
    if b["type"] != 0: continue
    for l in b["lines"]:
        for s in l["spans"]:
            txt = "".join(ch["c"] for ch in s["chars"])
            if not txt.strip(): continue
            spans.append(dict(text=txt, font=s["font"], size=round(s["size"], 3), color="#%06x" % s["color"], origin=[round(v, 3) for v in s["chars"][0]["origin"]],
                              bbox=[round(v, 3) for v in s["bbox"]], dir=[round(v, 3) for v in l["dir"]], asc=round(s.get("ascender", 0), 4), dsc=round(s.get("descender", 0), 4),
                              chars=[dict(c=ch["c"], x=round(ch["origin"][0], 3), bbox=[round(v, 3) for v in ch["bbox"]]) for ch in s["chars"]]))
rec = dict(sheet=S, input=IN, sha256=hashlib.sha256(open(IN, "rb").read()).hexdigest(), page=[p.rect.width, p.rect.height], mediabox=list(p.mediabox), cropbox=list(p.cropbox),
           fonts=sorted({(f[3], f[2]) for f in p.get_fonts(full=True)}), xobjects=len(p.get_xobjects()), images=len(p.get_images()), n_spans=len(spans), spans=spans)
if want_dr:
    items = []
    for it in p.get_drawings():
        r = it["rect"]; kinds = "".join(sorted(set(k[0] for k in it["items"])))
        e = dict(rect=[round(r.x0, 3), round(r.y0, 3), round(r.x1, 3), round(r.y1, 3)], type=it.get("type"), fill=it.get("fill"), color=it.get("color"), width=it.get("width"), dashes=it.get("dashes"),
                 lineCap=it.get("lineCap"), kinds=kinds, n=len(it["items"]), even_odd=it.get("even_odd"), fill_opacity=it.get("fill_opacity"), stroke_opacity=it.get("stroke_opacity"))
        if kinds == "l" and len(it["items"]) == 1:
            p0, p1 = it["items"][0][1], it["items"][0][2]; e["line"] = [round(p0.x, 3), round(p0.y, 3), round(p1.x, 3), round(p1.y, 3)]
        items.append(e)
    rec["drawings"] = items; rec["n_drawings"] = len(items)
d.close(); gc.collect()
json.dump(rec, open(OUT, "w"), indent=1, ensure_ascii=False, default=lambda o: list(o) if hasattr(o, "__iter__") else str(o))
print(S, "spans", len(spans), "drawings", rec.get("n_drawings"), "fonts", rec["fonts"], "page", rec["page"])
