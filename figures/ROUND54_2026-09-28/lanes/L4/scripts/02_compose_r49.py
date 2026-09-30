#!$T90_PY
"""ROUND 49 (2026-09-26, lane L4): the round-40 composer with this lane's +1 pt panel pages (work/panel_*_m.pdf, sha-gated against
work/build_manifest_r49.json), letters and the title at 14 pt, the coordinator's band_4e_r49.pdf, panel d's page placed at its widened box.
Round-40 docstring follows.
Round 40, LMAIN, Main_Fig4 compose (vector): the round-38 composer (lanes/Main_Fig4/scripts/02_compose.py) with NO V13 base
placement: the title "Figure 4" is stamped by TextWriter (Arial Bold 13 pt at (14.17, 14.0), x = the letter a's x) instead of the V13
title clip, and band e is the round-40 band placed with its top at the V13 band-e top (G["band_e"]["rect"][1] = 1083.38, inside the
V13 strip_with_letter [1061.308, 1293.14] whose upper 22 pt held the letter e); the letter e is stamped at G["letters"]["e"]. Page height
= band top + band height (= the V13 height when the band keeps its 209.76 pt). Panels a, b, c, d = the round-38 font-metric-aligned
pages panel_*_m.pdf (sha-gated against the round-38 compose log), letters a-d at the V13 origins (l4lib.stamp_letter).
usage: 02_compose_r40.py [--band PATH] [--out PATH]"""
import argparse, gc, json, os, sys
import fitz
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # round 49: lmain_lib.py and l4lib.py are the lane copies beside this script
import lmain_lib as L, l4lib as L4

SHEET = "Main_Fig4"; SD = f"{L.LANE}/{SHEET}"; WORK = f"{SD}/work"; R38W = WORK   # round 49: the panels are this lane's +1 pt rebuilds in work/
LETTER_SIZE = 14.0   # round 49: panel letters 13 -> 14 (the title through L.TITLE_SIZE = 14)
ap = argparse.ArgumentParser(); ap.add_argument("--band", default=f"{L.SKETCH}/band_4e_r49.pdf"); ap.add_argument("--out", default=f"{WORK}/{SHEET}_r49_vector.pdf"); a = ap.parse_args()
BAND = L.hydrated(a.band); OUT = a.out
G = json.load(open(f"{R38W}/v13_geometry.json")); BM = json.load(open(f"{WORK}/build_manifest_r49.json")); GD = json.load(open(f"{WORK}/panel_d_geometry.json"))
W, H13 = G["page"]
panels = {}
for x in "abcd":
    m = L.hydrated(f"{R38W}/panel_{x}_m.pdf"); assert L.sha256(m) == BM["panels"][x]["sha256"], (x, "panel page differs from the build manifest of this chain")
    box = {"a": G["panel_a"]["clip"], "b": G["panel_b"]["panel_box"], "c": G["panel_c"]["panel_box"], "d": GD["page_box"]}[x]   # round 49: panel d's page is its V13 box widened left by page_pad_left
    panels[x] = (m, fitz.Rect(*box))
assert GD["v13_panel_box"] == G["panel_d"]["panel_box"] and abs((GD["page_box"][0] + GD["page_pad_left"]) - G["panel_d"]["panel_box"][0]) < 1e-9 and GD["page_box"][1:] == G["panel_d"]["panel_box"][1:], (GD["page_box"], G["panel_d"]["panel_box"])
A = fitz.Rect(*G["panel_a"]["clip"]); assert list(panels["a"][1]) == list(A)
E_RECT = G["band_e"]["rect"]; E_STRIP = G["band_e"]["strip_with_letter"]; BAND_TOP = float(E_RECT[1]); LE = G["letters"]["e"]
assert E_STRIP[1] < LE["top"] < BAND_TOP and abs(E_RECT[3] - H13) < 0.01, (E_STRIP, LE, E_RECT, H13)
bd = fitz.open(BAND); bp = bd[0]; bw, bh = bp.rect.width, bp.rect.height
assert abs(bw - W) < 0.05, ("band width must equal the sheet width", bw, W)
H_NEW = round(BAND_TOP + bh, 3)
doc = fitz.open(); page = doc.new_page(width=W, height=H_NEW); log = {"placements": []}
for x in "abcd":
    m, box = panels[x]; pdoc = fitz.open(m); pr = pdoc[0].rect
    assert abs(pr.width - box.width) < 0.02 and abs(pr.height - box.height) < 0.02, (x, pr, box)
    page.show_pdf_page(box, pdoc, 0); pdoc.close(); log["placements"].append({"what": f"panel {x} (round-49 +1 pt page)", "rect": list(box), "file": m, "sha256": L.sha256(m)})
for ch in "abcde":
    le = G["letters"][ch]; L4.stamp_letter(page, ch, le["x"], le["baseline"], LETTER_SIZE, color=(0, 0, 0)); log["placements"].append({"what": f"letter {ch} (TextWriter Arial Bold {LETTER_SIZE:g} pt, V13 origin)", "x": le["x"], "baseline": le["baseline"], "size": LETTER_SIZE})
page.show_pdf_page(fitz.Rect(0, BAND_TOP, W, BAND_TOP + bh), bd, 0); log["placements"].append({"what": "band e (round 49, the coordinator's band), top = the V13 band-e top", "rect": [0, BAND_TOP, W, H_NEW], "file": BAND, "sha256": L.sha256(BAND), "band_page": [round(bw, 3), round(bh, 3)]})
TITLE_X = G["letters"]["a"]["x"]
tw = fitz.TextWriter(page.rect); tw.append(fitz.Point(TITLE_X, L.TITLE_BASELINE), "Figure 4", font=fitz.Font(fontfile=L.ARIALB), fontsize=L.TITLE_SIZE); tw.write_text(page, color=L.INK_RGB)
log["placements"].append({"what": f"title 'Figure 4' (TextWriter Arial Bold {L.TITLE_SIZE:g} pt, x = letter a x, baseline 14.0)", "x": TITLE_X, "baseline": L.TITLE_BASELINE, "size": L.TITLE_SIZE})
bd.close(); doc.save(OUT, garbage=4, deflate=True); doc.close(); gc.collect()
n, w, h = L.gs_pdfinfo(OUT); assert n == 1 and abs(w - W) < 0.01 and abs(h - H_NEW) < 0.01, (n, w, h, H_NEW)
log.update({"sheet": SHEET, "band": BAND, "band_top": BAND_TOP, "v13_band_rect": E_RECT, "v13_strip_with_letter": E_STRIP, "letter_e": LE, "page": [w, h], "page_v13": [W, H13], "output": OUT, "output_sha256": L.sha256(OUT), "bytes": os.path.getsize(OUT), "written": L.now()})
json.dump(log, open(f"{WORK}/compose_log{'' if OUT.endswith('_vector.pdf') else '_' + os.path.basename(OUT)[:-4]}.json", "w"), indent=1)
log["build_manifest"] = BM; log["letter_size"] = LETTER_SIZE; log["title_size"] = L.TITLE_SIZE; json.dump(log, open(f"{WORK}/compose_log{'' if OUT.endswith('_vector.pdf') else '_' + os.path.basename(OUT)[:-4]}.json", "w"), indent=1)
print(f"wrote {OUT}: page {w} x {h} (gs PDFINFO; V13 {W} x {H13}), band {os.path.basename(BAND)} {bw:.2f} x {bh:.3f} at top {BAND_TOP}, title at ({TITLE_X}, {L.TITLE_BASELINE}) {L.TITLE_SIZE:g} pt, letters {LETTER_SIZE:g} pt, letter e at ({LE['x']}, {LE['baseline']})")
