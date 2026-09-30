#!$T90_PY
"""ROUND 54 (2026-09-29, lane L4): the composes. (1) Main_Fig4 = panels a and b (this lane's font-metric-aligned pages work/panel_a_m.pdf
and panel_b_m.pdf, sha-gated against work/build_manifest_r54.json) at their spec boxes, letters a, b and c (Arial Bold 18 pt, TextWriter)
at the spec origins, the title 'Figure 4' (Arial Bold 18 pt) on baseline 18, the coordinator's band_4e_r54.pdf as panel c placed with
its top at the spec's band_top (page height = band top + the band's page-box height, read from the file). Owner's option A: V30's
panels c and d are not on the sheet. (2) The counts panel (V30's panel c) as a standalone single page panel_counts_r54/panel_counts_r54.pdf
(white page of the spec's box, no letter, no title). (3) Ghostscript pdfwrite re-distillations of both (the round-44 gs_flatten recipe)
for the coordinator's text census: work/Main_Fig4_flat.pdf and panel_counts_r54/work/panel_counts_r54_flat.pdf.
usage: 02_compose_r54.py [--band PATH]"""
import argparse, gc, json, os, subprocess, sys
import fitz
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lmain_lib as L, l4lib as L4, layout_r54 as LY

SHEET = "Main_Fig4"; SD = f"{L.LANE}/{SHEET}"; WORK = f"{SD}/work"; PC = f"{L.LANE}/panel_counts_r54"; os.makedirs(f"{PC}/work", exist_ok=True); os.makedirs(f"{PC}/verify", exist_ok=True)
S = LY.load(); FS = S["sizes"]
ap = argparse.ArgumentParser(); ap.add_argument("--band", default=S["e"]["band"]); a = ap.parse_args()
BAND = L.hydrated(a.band); OUT = f"{SD}/{SHEET}.pdf"; FLAT = f"{WORK}/{SHEET}_flat.pdf"; OUT_C = f"{PC}/panel_counts_r54.pdf"; FLAT_C = f"{PC}/work/panel_counts_r54_flat.pdf"
BM = json.load(open(f"{WORK}/build_manifest_r54.json")); W = S["page_w"]


def gs_flatten(src, out):
    """Re-distil a nested sheet into a single-page pdfwrite file (no downsampling, colours untouched, fonts embedded), the round-44 recipe."""
    r = subprocess.run([L.GS, "-dSAFER", "-dBATCH", "-dNOPAUSE", "-dQUIET", "-sDEVICE=pdfwrite", "-dCompatibilityLevel=1.7", "-dPDFSETTINGS=/prepress",
                        "-dDownsampleColorImages=false", "-dDownsampleGrayImages=false", "-dDownsampleMonoImages=false", "-dAutoFilterColorImages=false",
                        "-dAutoFilterGrayImages=false", "-dColorImageFilter=/FlateEncode", "-dGrayImageFilter=/FlateEncode", "-dColorConversionStrategy=/LeaveColorUnchanged",
                        "-dEmbedAllFonts=true", "-dSubsetFonts=true", "-dPreserveAnnots=false", f"-sOutputFile={out}", src], capture_output=True, text=True, timeout=600)
    assert r.returncode == 0 and os.path.exists(out) and os.path.getsize(out) > 0, ("gs pdfwrite failed", r.stderr[-400:])
    return out


def stamp(page, text, x, baseline, size):
    tw = fitz.TextWriter(page.rect); tw.append(fitz.Point(x, baseline), text, font=fitz.Font(fontfile=L.ARIALB), fontsize=size); tw.write_text(page, color=L.INK_RGB)


# ------------------------------------------------------------------ (1) Main_Fig4 = a | b + band (panel c)
E = S["e"]; bd = fitz.open(BAND); bp = bd[0]; bw, bh = bp.rect.width, bp.rect.height
assert abs(bw - W) < 0.05 and abs(bh - E["band_h"]) < 0.01, ("band page box differs from the spec's reading", bw, bh, E["band_h"])
H_NEW = round(E["band_top"] + bh, 3); assert abs(H_NEW - E["page_h"]) < 0.01
doc = fitz.open(); page = doc.new_page(width=W, height=H_NEW); log = {"placements": []}
for x in "ab":
    m = L.hydrated(f"{WORK}/panel_{x}_m.pdf"); assert L.sha256(m) == BM["panels"][x]["sha256"], (x, "panel page differs from the build manifest of this chain")
    box = fitz.Rect(*S[x]["box"]); pdoc = fitz.open(m); pr = pdoc[0].rect
    assert abs(pr.width - box.width) < 0.02 and abs(pr.height - box.height) < 0.02, (x, pr, box)
    page.show_pdf_page(box, pdoc, 0); pdoc.close(); log["placements"].append({"what": f"panel {x} (round-54 page)", "rect": list(box), "file": m, "sha256": L.sha256(m)})
for ch in "abc":
    le = S["letters"][ch]; stamp(page, ch, le["x"], le["base"], FS["letter"]); log["placements"].append({"what": f"letter {ch} (TextWriter Arial Bold {FS['letter']:g} pt)", "x": le["x"], "baseline": le["base"], "size": FS["letter"]})
page.show_pdf_page(fitz.Rect(0, E["band_top"], W, E["band_top"] + bh), bd, 0)
log["placements"].append({"what": "band 4e (round 54, the coordinator's band) as panel c, top = the spec's band_top", "rect": [0, E["band_top"], W, H_NEW], "file": BAND, "sha256": L.sha256(BAND), "band_page": [round(bw, 3), round(bh, 3)], "band_mediabox": E["band_mediabox"]})
T = S["title"]; stamp(page, T["text"], T["x"], T["base"], FS["title"]); log["placements"].append({"what": f"title '{T['text']}' (TextWriter Arial Bold {FS['title']:g} pt, baseline {T['base']})", "x": T["x"], "baseline": T["base"], "size": FS["title"]})
bd.close(); doc.save(OUT, garbage=4, deflate=True); doc.close(); gc.collect()
n, w, h = L.gs_pdfinfo(OUT); assert n == 1 and abs(w - W) < 0.01 and abs(h - H_NEW) < 0.01, (n, w, h, H_NEW)
gs_flatten(OUT, FLAT); nf, wf, hf = L.gs_pdfinfo(FLAT); assert nf == 1 and abs(wf - w) < 0.05 and abs(hf - h) < 0.05, ("the flat twin's page differs", nf, wf, hf)
log.update({"sheet": SHEET, "band": BAND, "band_top": E["band_top"], "letters": S["letters"], "title": T, "page": [w, h], "output": OUT, "output_sha256": L.sha256(OUT), "bytes": os.path.getsize(OUT),
            "flat": FLAT, "flat_sha256": L.sha256(FLAT), "flat_page": [wf, hf], "build_manifest": BM, "written": L.now(), "arrangement": S["arrangement"]})
json.dump(log, open(f"{WORK}/compose_log.json", "w"), indent=1)
print(f"wrote {OUT}: page {w} x {h} (gs PDFINFO), band {os.path.basename(BAND)} {bw:.2f} x {bh:.3f} at top {E['band_top']:.3f}, letters a b c {FS['letter']:g} pt, title on baseline {T['base']}; flat twin {FLAT} {wf} x {hf}")

# ------------------------------------------------------------------ (2) the counts panel, a standalone single page
mc = L.hydrated(f"{WORK}/panel_c_m.pdf"); assert L.sha256(mc) == BM["panels"]["c"]["sha256"]
cbox = fitz.Rect(*S["c"]["box"]); pdoc = fitz.open(mc); pr = pdoc[0].rect; assert abs(pr.width - cbox.width) < 0.02 and abs(pr.height - cbox.height) < 0.02, (pr, cbox)
doc = fitz.open(); page = doc.new_page(width=cbox.width, height=cbox.height); page.show_pdf_page(cbox, pdoc, 0); pdoc.close(); doc.save(OUT_C, garbage=4, deflate=True); doc.close(); gc.collect()
nc, wc, hc = L.gs_pdfinfo(OUT_C); assert nc == 1 and abs(wc - cbox.width) < 0.01 and abs(hc - cbox.height) < 0.01, (nc, wc, hc)
gs_flatten(OUT_C, FLAT_C); nfc, wfc, hfc = L.gs_pdfinfo(FLAT_C); assert nfc == 1 and abs(wfc - wc) < 0.05 and abs(hfc - hc) < 0.05
json.dump({"panel": "counts (V30's Main_Fig4 panel c)", "page": [wc, hc], "output": OUT_C, "output_sha256": L.sha256(OUT_C), "panel_page": mc, "panel_sha256": L.sha256(mc), "flat": FLAT_C, "flat_sha256": L.sha256(FLAT_C),
           "letters": "none", "title": "none", "written": L.now()}, open(f"{PC}/work/compose_log.json", "w"), indent=1)
print(f"wrote {OUT_C}: page {wc} x {hc}; flat twin {FLAT_C}")
