#!$T90_PY
import fitz, sys, json, re
NEW = sys.argv[1]
d = fitz.open(NEW); p = d[0]
print("xrefs contents:", p.get_contents()); cs = b"".join(d.xref_stream(x) for x in p.get_contents()); print("content bytes", len(cs))
txt = cs.decode("latin-1")
print("q count", txt.count("\nq\n") + txt.count(" q\n"), "Q count", txt.count("\nQ\n") + txt.count(" Q\n"), "'W n' count", txt.count("W n"), "BT count", txt.count("BT"), "Tf count", txt.count("Tf"))
print("HEAD:", txt[:400].replace("\n", " | ")); print("TAIL:", txt[-600:].replace("\n", " | "))
print("words:", len(p.get_text("words")), "dict blocks:", len(p.get_text("dict")["blocks"]), "drawings:", len(p.get_drawings()))
print("fonts:", p.get_fonts(full=True)[:4]); print("rect", p.rect, "mediabox", p.mediabox, "cropbox", p.cropbox, "rotation", p.rotation)
