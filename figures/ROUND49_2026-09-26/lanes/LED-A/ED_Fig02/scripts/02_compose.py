#!/usr/bin/env python3
"""ROUND 49 lane LED-A copy (2026-09-26): the anisotropic key and caption strings get Tf + 1 at the same x-scale, the letters 14 pt, no PyMuPDF
renders (Ghostscript renders and crops are made by run_sheet.py). Original: ROUND38_2026-09-16/figures/lanes/ED_Fig02/scripts/02_compose.py.
ED_Fig02, lane V7_L2_T90 (relaunch 3), step 2: compose the sheet.
1. align the matplotlib Type 42 fonts' FontDescriptor Ascent and Descent to the base sheet's own Arial descriptors (ArialMT 728 / -210,
   Arial-BoldMT 1010 / -376, read from the base by 00b_probe_visible.py) so the text layer reports the same span boxes (round-28 idiom),
2. stamp the panel letters a, b, c with fitz.TextWriter and /System/Library/Fonts/Supplemental/Arial Bold.ttf at 13 pt at the base
   sheet's letter origins (baseline 27.72, the base's own Arial-BoldMT 13 pt letters),
3. save <Sheet>.pdf (one flat page), render ED_Fig02_150dpi.png, and cut the 200 dpi OLD against NEW crops per panel into verify/crops/.
usage: 02_compose.py [--old]  (composes work/sheet_mpl_old.pdf into work/ED_Fig02_OLDREPLOT.pdf for the eye)"""
import json, os, re, sys, gc
import fitz
DELTA = float(os.environ.get("ED2_DELTA", "1.0"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ed2_common import *

OLDMODE = "--old" in sys.argv
SRC_PDF = os.path.join(WORK, "sheet_mpl_old.pdf" if OLDMODE else "sheet_mpl.pdf")
OUT_PDF = os.path.join(WORK, "ED_Fig02_OLDREPLOT.pdf") if OLDMODE else os.path.join(WORK, "ED_Fig02_composed.pdf")   # round 49: the strip is stamped afterwards into LANE/ED_Fig02.pdf
G = load_geometry(); W, H = G["page"]
BASEV = json.load(open(hydrated(os.path.join(WORK, "base_visible.json"))))
# the base sheet's descriptors per face: regular text spans 728/-210, bold 1010/-376 (the digits subset 718/-12 is a base artifact, not reproduced)
MET = {"ArialMT": (728, -210), "Arial-BoldMT": (1010, -376)}
for f in BASEV["fonts"]:
    base = f["basefont"].split("+")[-1]
    if base == "Arial-BoldMT": assert (int(f["descriptor"]["Ascent"]), int(f["descriptor"]["Descent"])) == MET[base], f
    if base == "ArialMT": assert (int(f["descriptor"]["Ascent"]), int(f["descriptor"]["Descent"])) in (MET[base], (718, -12)), f

doc = fitz.open(hydrated(SRC_PDF)); page = doc[0]
assert len(doc) == 1 and len(page.get_xobjects()) == 0, (len(doc), len(page.get_xobjects()))
assert abs(page.rect.width - W) < 0.01 and abs(page.rect.height - H) < 0.01, page.rect
# RELAUNCH 4: the base sheet writes its 13 key labels and 3 key letters with Tf 8.1881 and a text matrix scaling x by 1.05873, its 7
# notes with Tf 7.7787 and 1.05874, and the panel-a caption 'Rate of new diagnosis' with Tf 10.0854 and 0.991532 (an anisotropic
# artifact of the old build: the text layer reports the geometric-mean sizes 8.425, 8.004 and 10.043). matplotlib typeset exactly
# those 24 strings uniformly at the mean sizes, 3 percent taller and 3 percent narrower than the base's. Rewrite those text objects
# (matplotlib writes each string as 'q 1 0 0 1 x y cm BT /Fn size Tf 0 0 Td [..] TJ ET Q') to the base's Tf and Tm so glyph height
# and advances equal the base sheet's, keeping the origin (the cm translation) and the text layer.
_bs = [s_ for s_ in BASEV["spans"] if s_.get("visible")]
# round 49: the same construction one point larger in glyph height: Tf + DELTA, the same x-scale (Ghostscript reports the Tf, so the census reads
# 8.188 -> 9.188 and 10.085 -> 11.085); the seven 8.004 note strings are not on the sheet (round 40), so their count is 0 at DELTA 1
def _r(v): return f"{v:.10f}".rstrip("0").rstrip(".")
STRETCH = {_r(8.425 + DELTA): (_r(8.1881 + DELTA), "1.05873", sum(1 for s_ in _bs if abs(s_["size"] - 8.425) < 0.01)),
           _r(8.004 + DELTA): (_r(7.7787 + DELTA), "1.05874", sum(1 for s_ in _bs if abs(s_["size"] - 8.004) < 0.01) if OLDMODE else 0),
           _r(10.043 + DELTA): (_r(10.0854 + DELTA), "0.991532", sum(1 for s_ in _bs if abs(s_["size"] - 10.043) < 0.01))}   # V13: 14, 7 (0 at +1: notes removed), 1
xref = page.get_contents()[0]; cs = doc.xref_stream(xref).decode("latin-1"); patched = {}
for mean, (tf, a, expect) in STRETCH.items():
    cs, k = re.subn(r"(/F\d+) " + re.escape(mean) + r" Tf\n0 0 Td\n", lambda m: f"{m.group(1)} {tf} Tf\n{a} 0 0 1 0 0 Tm\n", cs); patched[mean] = k
    assert k == expect, (mean, k, expect)
doc.update_stream(xref, cs.encode("latin-1")); print("text objects given the base sheet's Tf and x-scaled Tm:", patched)
done = {}
for f in page.get_fonts(full=True):
    base = f[3].split("+")[-1]; assert base in MET, (f, "not an Arial face the base sheet carries")
    desc = doc.xref_get_key(f[0], "DescendantFonts"); assert desc[0] == "array", (f, desc)
    cid = int(re.search(r"(\d+) 0 R", desc[1]).group(1)); fdx = int(doc.xref_get_key(cid, "FontDescriptor")[1].split()[0])
    before = (doc.xref_get_key(fdx, "Ascent")[1], doc.xref_get_key(fdx, "Descent")[1])
    doc.xref_set_key(fdx, "Ascent", str(MET[base][0])); doc.xref_set_key(fdx, "Descent", str(MET[base][1])); done[base] = dict(before=before, after=MET[base])
print("font metrics aligned:", done)
# letters by TextWriter at the base origins (Arial Bold 13 pt); the base letters are Arial-BoldMT 13 pt at baseline 27.72
letters = sorted([s for s in G["visible_spans"] if abs(s["size"] - 13.0) < 0.01], key=lambda s: s["origin"][0])
assert [s["text"] for s in letters] == ["a", "b", "c"] and all(s["font"] == "Arial-BoldMT" for s in letters), letters
font = fitz.Font(fontfile=ARIAL_B); tw = fitz.TextWriter(page.rect)
for s in letters: tw.append(fitz.Point(s["origin"][0], s["origin"][1]), s["text"], font=font, fontsize=13.0 + DELTA)     # round 49: 14 pt
ink = tuple(int(INK[i:i + 2], 16) / 255 for i in (1, 3, 5))
tw.write_text(page, color=ink)
# the TextWriter font (BaseFont "Arial Bold", Type0) gets the base sheet's Arial-BoldMT descriptor too, so the letters' span boxes match
for f in page.get_fonts(full=True):
    if f[3].split("+")[-1] == "Arial Bold":
        desc = doc.xref_get_key(f[0], "DescendantFonts"); assert desc[0] == "array", (f, desc)
        cid = int(re.search(r"(\d+) 0 R", desc[1]).group(1)); fdx = int(doc.xref_get_key(cid, "FontDescriptor")[1].split()[0])
        before = (doc.xref_get_key(fdx, "Ascent")[1], doc.xref_get_key(fdx, "Descent")[1])
        doc.xref_set_key(fdx, "Ascent", str(MET["Arial-BoldMT"][0])); doc.xref_set_key(fdx, "Descent", str(MET["Arial-BoldMT"][1])); done["Arial Bold (TextWriter)"] = dict(before=before, after=MET["Arial-BoldMT"])
print("font metrics aligned:", done)
doc.save(OUT_PDF, garbage=1, deflate=True); doc.close(); gc.collect()
# read back: letters, fonts, page
doc = fitz.open(OUT_PDF); page = doc[0]
assert len(doc) == 1 and len(page.get_xobjects()) == 0
spans = [dict(text="".join(ch["c"] for ch in s["chars"]), font=s["font"], size=round(s["size"], 3), origin=[round(s["chars"][0]["origin"][0], 3), round(s["chars"][0]["origin"][1], 3)],
              asc=round(s["ascender"], 4), desc=round(s["descender"], 4), color="#%06x" % s["color"], bbox=[round(v, 3) for v in s["bbox"]])
         for b in page.get_text("rawdict")["blocks"] if b["type"] == 0 for l in b["lines"] for s in l["spans"]]
print("fonts in the composed sheet:", [(f[3], f[2], f[5]) for f in page.get_fonts(full=True)])
for s in letters:
    m = [q for q in spans if q["text"] == s["text"] and abs(q["size"] - 13 - DELTA) < 0.01]
    print(f"letter {s['text']}: base origin {s['origin']} asc {s['asc']} desc {s['desc']} | new {[(q['origin'], q['font'], q['asc'], q['desc']) for q in m]}")
import collections
print("spans", len(spans), "sizes", collections.Counter(s["size"] for s in spans), "fonts", collections.Counter(s["font"] for s in spans))
print("asc/desc by font:", collections.Counter((s["font"], s["asc"], s["desc"]) for s in spans))
print("sample:", [(s["text"], s["origin"], s["size"]) for s in spans[:8]])
json.dump(spans, open(os.path.join(WORK, "new_text_old.json" if OLDMODE else "new_text.json"), "w"), indent=0)
print("wrote", OUT_PDF, "(renders and crops: run_sheet.py, Ghostscript)")
