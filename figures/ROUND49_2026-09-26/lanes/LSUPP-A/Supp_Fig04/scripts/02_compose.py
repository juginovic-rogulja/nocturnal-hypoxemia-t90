#!/usr/bin/env python3
"""Supp_Fig04, V14 lane L2, step 2: the finished sheet from the re-plotted page (repointed copy of the round-30 composer, byte copy beside as
02_compose_PRE_V8_1.py): the three matplotlib fonts' FontDescriptor Ascent and Descent aligned to the values the V13 sheet's text layer
reports for the same faces (ArialMT 728/-210, Arial-BoldMT 1010/-376, Arial-ItalicMT 997/-324), saved as Supp_Fig04.pdf at the folder
root (page width = V13, height = V13 + the declared growth). Renders Supp_Fig04_150dpi.png and the 200 dpi OLD (V13) vs NEW crops of the
four panels. The positive-control page (work/panels_OLD.pdf) is aligned the same way to work/Supp_Fig04_OLD.pdf for 03_verify."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import gc, json, os, re, sys
import fitz, numpy as np
from PIL import Image
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND37_2026-09-14/figures/_common")
from v14lib import V13, hydrated, RenderLock
LANE = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); WORK = os.path.join(LANE, "work"); VER = os.path.join(LANE, "verify"); CROPS = os.path.join(VER, "crops"); os.makedirs(CROPS, exist_ok=True)
BASE = f"{V13}/Supp_Fig04.pdf"
PANELS = os.path.join(WORK, "panels.pdf"); PANELS_OLD = os.path.join(WORK, "panels_OLD.pdf"); OUT = os.path.join(LANE, "Supp_Fig04.pdf"); OUT_OLD = os.path.join(WORK, "Supp_Fig04_OLD.pdf"); PNG150 = os.path.join(LANE, "Supp_Fig04_150dpi.png")
DR = json.load(open(hydrated(os.path.join(VER, "Supp_Fig04_drawn.json")))); GROW = DR["growth"]["cd"]; PAGE_H = DR["page"][1]
def sheet_font_metrics(pdf):
    doc = fitz.open(hydrated(pdf)); page = doc[0]; cnt = {}
    for b in page.get_text("rawdict")["blocks"]:
        if b["type"] != 0: continue
        for l in b["lines"]:
            for s in l["spans"]:
                k = (s["font"], int(round(s["ascender"] * 1000)), int(round(s["descender"] * 1000))); cnt[k] = cnt.get(k, 0) + 1
    doc.close(); gc.collect(); met = {}
    for (font, asc, dsc), n in sorted(cnt.items(), key=lambda kv: -kv[1]): met.setdefault(font, (asc, dsc))
    return met, cnt
def align_font_metrics(src, dst, met, height):
    doc = fitz.open(hydrated(src)); page = doc[0]; done = {}
    for f in page.get_fonts(full=True):
        xref, base = f[0], f[3].split("+")[-1]; assert base in met, (base, sorted(met)); asc, dsc = met[base]
        desc = doc.xref_get_key(xref, "DescendantFonts")
        if desc[0] == "array": cid = int(re.search(r"(\d+) 0 R", desc[1]).group(1)); fdx = int(doc.xref_get_key(cid, "FontDescriptor")[1].split()[0])
        else: fdx = int(doc.xref_get_key(xref, "FontDescriptor")[1].split()[0])
        doc.xref_set_key(fdx, "Ascent", str(asc)); doc.xref_set_key(fdx, "Descent", str(dsc)); done[base] = (asc, dsc)
    assert abs(page.rect.width - 1021.45) < 0.01 and abs(page.rect.height - height) < 0.02, (page.rect, height)
    assert len(page.get_xobjects()) == 0, "the re-plotted page must be flat"
    doc.save(dst, garbage=4, deflate=True); doc.close(); gc.collect(); return done
def render(pdf, dpi, clip=None):
    doc = fitz.open(hydrated(pdf)); page = doc[0]; pm = page.get_pixmap(dpi=dpi, colorspace=fitz.csRGB, alpha=False, clip=None if clip is None else fitz.Rect(*clip))
    img = Image.frombytes("RGB", (pm.width, pm.height), pm.samples); pm = None; doc.close(); gc.collect(); return img
met, cnt = sheet_font_metrics(BASE); print("V13 font descriptors (majority):", met)
done = align_font_metrics(PANELS, OUT, met, PAGE_H); print("aligned", done, "->", OUT, os.path.getsize(OUT), "bytes")
done_old = align_font_metrics(PANELS_OLD, OUT_OLD, met, 989.728); print("aligned", done_old, "->", OUT_OLD)
img = render(OUT, 150); img.save(PNG150); print("wrote", PNG150, img.size)
CROP = dict(a=(14.0, 14.0, 499.0, 560.0 + DR["growth"]["a"]), b=(522.0, 14.0, 1008.0, 700.0 + DR["growth"]["b"]), c=(14.0, 698.0 + GROW, 500.0, 976.0 + GROW), d=(522.0, 698.0 + GROW, 1008.0, 976.0 + GROW))
CROP13 = dict(a=(14.0, 14.0, 499.0, 560.0), b=(522.0, 14.0, 1008.0, 700.0), c=(14.0, 698.0, 500.0, 976.0), d=(522.0, 698.0, 1008.0, 976.0))
with RenderLock():
    for p in CROP:
        old = render(BASE, 200, CROP13[p]); new = render(OUT, 200, CROP[p])
        old.save(os.path.join(CROPS, f"panel_{p}_old_200dpi.png")); new.save(os.path.join(CROPS, f"panel_{p}_new_200dpi.png"))
        gap = 24; sbs = Image.new("RGB", (old.width + gap + new.width, max(old.height, new.height)), (255, 255, 255)); sbs.paste(old, (0, 0)); sbs.paste(new, (old.width + gap, 0))
        a = np.asarray(sbs).copy(); a[:, old.width + 10:old.width + 14, :] = (200, 200, 200); Image.fromarray(a).save(os.path.join(CROPS, f"panel_{p}_OLD_vs_NEW_200dpi.png")); print(f"crops panel {p}: {old.size} old | {new.size} new")
json.dump(dict(base=BASE, out=OUT, out_old=OUT_OLD, font_metrics_applied=done, crops=CROP, crops_v13=CROP13, png150=PNG150, page_height=PAGE_H), open(os.path.join(WORK, "compose_record.json"), "w"), indent=1)
