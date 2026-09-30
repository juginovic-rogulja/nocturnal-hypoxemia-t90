#!$T90_PY
"""Lane LSUPP-D, round 50 (2026-09-26), brief item 5: Supp_Fig16 re-laid as band (top, left-aligned, unchanged), panels a and b
side by side (a left, b right, tops aligned), panel c below at the left, on a wider and shorter page.

Pieces: the round-49 flat sheet (ROUND49 lane work/Supp_Fig16_restamped.pdf, the body form of V27's Supp_Fig16: every string at
its round-49 size, band-free) cut into three text-exact pieces by text-only redaction of everything outside each panel's region
(so that the Ghostscript census of the composed sheet carries every string once, a clipped placement alone would report the
whole page three times), then placed by show_pdf_page with clips on a new page: piece a at its own position, pieces b and c
translated by whole 150 dpi pixels (dx, dy multiples of 0.48 pt) so that the render is a pixel-exact translation of V27's.
The band (band_s16_r49.pdf) is placed as in round 49. Letters travel with their panels (each sits at its panel's top-left corner).
usage: relayout_s16.py [page_width_pt]   (default: b follows a with the sheet's own 32 pt gutter; 969.12 gives the main-sheet width)"""
import json, os, shutil, sys
from collections import Counter
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lsd_common as C

SHEET = "Supp_Fig16"; D = f"{C.LANE}/{SHEET}"; W = f"{D}/work"; V = f"{D}/verify"
for d_ in (W, V): os.makedirs(d_, exist_ok=True)
SRC = C.hydrated(f"{C.V27}/{SHEET}.pdf"); FLAT0 = C.hydrated(f"{C.R49LANE}/{SHEET}/work/{SHEET}_restamped.pdf"); BAND = C.hydrated(C.BAND_R49)
assert C.sha256(SRC).startswith("6e8eb1c3172d7b35"), C.sha256(SRC)            # V27 Supp_Fig16 = the round-49 15:32 build
assert C.sha256(FLAT0).startswith("26a4b3d4f0f01a91"), C.sha256(FLAT0)         # its body form (round-49 restamped flat, band-free)
assert C.sha256(BAND).startswith("13e8b4b2042422cb"), C.sha256(BAND)           # the band of the 15:32 rebuild (title 11 pt)
PW0, PH0 = 496.8, 1131.7; BAND_H = 79.0
DPI = 150; S = DPI / 72.0; PX = 72.0 / DPI                                     # 0.48 pt = one pixel at 150 dpi
# panel regions on the round-49 sheet (from its dumps: no drawing, no white fill and no string crosses y = 418 or y = 818,
# panel a's ink ends at x = 271.5, letters a, b, c at y 92.4, 429.8, 834.2)
REG = {"a": [0.0, BAND_H, 285.0, 418.0], "b": [0.0, 418.0, PW0, 818.0], "c": [0.0, 818.0, PW0, PH0]}
A_RIGHT, GUTTER, B_LEFT_MARGIN = 271.5, 32.2, 14.0
DY = -round((429.8 - 92.4) / PX) * PX                                          # b's letter top onto a's letter top, whole pixels: -337.44
alt = float(sys.argv[1]) if len(sys.argv) > 1 else None
DX = (round((A_RIGHT + GUTTER - B_LEFT_MARGIN) / PX) * PX) if alt is None else (round((alt - PW0) / PX) * PX)   # 289.44 (or 472.32 for 969.12)
PW, PH = round(DX + PW0, 2), round(PH0 + DY, 2)
TAG = "" if alt is None else f"_w{int(round(PW))}"
NEW = f"{D}/{SHEET}.pdf" if alt is None else f"{W}/{SHEET}{TAG}.pdf"
ck = C.Checks(f"{SHEET}, lane LSUPP-D round 50 (item 5, {C.now()}): band top left, a and b side by side, c below. Source V27 {SRC} (sha256 {C.sha256(SRC, 16)}), body form {FLAT0} (sha256 {C.sha256(FLAT0, 16)}), band {BAND} (sha256 {C.sha256(BAND, 16)}). Page {PW0} x {PH0} -> {PW} x {PH} pt, dx {DX} dy {DY} (whole 150 dpi pixels).")
FLAT = f"{W}/{SHEET}_body.pdf"; shutil.copyfile(FLAT0, FLAT)
# ------------------------------------------------------------------ 1 text-exact pieces
OFF = {"a": (0.0, 0.0), "b": (DX, DY), "c": (0.0, DY)}
PIECES = {}
for k, r in REG.items():
    comp = [[0, 0, PW0, r[1]], [0, r[3], PW0, PH0]] + ([[r[2], r[1], PW0, r[3]]] if r[2] < PW0 else [])
    out = f"{W}/piece_{k}{TAG}.pdf"; job = {"src": FLAT, "out": out, "whole_page": False, "redact_rects": comp, "items": []}
    json.dump(job, open(f"{W}/piece_{k}{TAG}_job.json", "w")); C.wd(f"{SHEET}_piece_{k}{TAG}", [C.PY, f"{C.SCRIPTS}/fitz_restamp.py", f"{W}/piece_{k}{TAG}_job.json"])
    runs = C.gs_txtwrite(out, f"{W}/piece_{k}{TAG}_census.xml"); PIECES[k] = (out, runs)
    inside = all(r[0] - 1.5 <= s_["x0"] and s_["x1"] <= r[2] + 1.5 and r[1] <= s_["y"] <= r[3] for s_ in runs)
    ck.log(f"piece {k}: every remaining run lies inside its region {r} ({len(runs)} runs)", inside, f"outside: {[(s_['text'], s_['x0'], s_['y']) for s_ in runs if not (r[0] - 1.5 <= s_['x0'] and s_['x1'] <= r[2] + 1.5 and r[1] <= s_['y'] <= r[3])][:5]}")
# ------------------------------------------------------------------ 2 compose
parts = []
for k in ("a", "b", "c"):
    r = REG[k]; dx, dy = OFF[k]
    parts.append({"pdf": PIECES[k][0], "rect": [r[0] + dx, r[1] + dy, r[2] + dx, r[3] + dy], "clip": r})
nb_, bw, bh = C.page_box(BAND); parts.append({"pdf": BAND, "rect": [0, 0, bw, bh]})
job = {"out": NEW, "page": [PW, PH], "parts": parts, "metadata": {"title": SHEET, "creator": f"LSUPP-D r50 relayout_s16.py: V27 body pieces (a at place, b dx {DX} dy {DY}, c dy {DY}) + band_s16_r49 (13e8b4b2)"}}
json.dump(job, open(f"{W}/compose{TAG}_job.json", "w")); C.wd(f"{SHEET}_compose{TAG}", [C.PY, f"{C.SCRIPTS}/fitz_compose.py", f"{W}/compose{TAG}_job.json"])
# ------------------------------------------------------------------ 3 verify
nn, nw, nh = C.page_box(NEW)
ck.log("one page, page box as planned", nn == 1 and abs(nw - PW) < 0.05 and abs(nh - PH) < 0.05, f"{nw:.2f} x {nh:.2f} pt (was {PW0} x {PH0})")
half = lambda v: round(v * 2) / 2
ref = C.gs_txtwrite(SRC, f"{W}/census_ref.xml"); new = C.gs_txtwrite(NEW, f"{W}/census_new{TAG}.xml")
cref = Counter((s_["text"], half(s_["size"])) for s_ in ref); cnew = Counter((s_["text"], half(s_["size"])) for s_ in new)
ck.log(f"txtwrite census equal to V27's: same multiset of strings and sizes ({sum(cnew.values())} runs, band included)", cref == cnew, f"difference {dict((cref - cnew) + (cnew - cref))}")
# positions: every V27 run reappears translated by its piece's offset
def piece_of(s_):
    if s_["y"] < BAND_H: return "band"
    for k, r in REG.items():
        if r[1] <= s_["y"] <= r[3]: return k
    return None
unmatched = []; used = [False] * len(new)
for s_ in ref:
    k = piece_of(s_); dx, dy = OFF.get(k, (0.0, 0.0))
    hit = [j for j, t in enumerate(new) if not used[j] and t["text"] == s_["text"] and abs(half(t["size"]) - half(s_["size"])) < 0.011 and abs(t["x0"] - (s_["x0"] + dx)) <= 1.0 and abs(t["y"] - (s_["y"] + dy)) <= 1.0]
    if hit: used[hit[0]] = True
    else: unmatched.append((s_["text"], k, s_["x0"], s_["y"]))
ck.log("every V27 run reappears at its piece's translation within 1 pt (band and a at place, b by (dx, dy), c by dy)", not unmatched and all(used), f"unmatched {unmatched[:6]}")
band_own = Counter((s_["text"], half(s_["size"])) for s_ in C.gs_txtwrite(BAND, f"{W}/census_band.xml"))
band_new = Counter((s_["text"], half(s_["size"])) for s_ in new if s_["y"] < BAND_H)
ck.log("band region census equals band_s16_r49's own census", band_new == band_own, f"{sum(band_new.values())} runs")
letters = {s_["text"]: (s_["x0"], s_["y"]) for s_ in new if s_["text"] in ("a", "b", "c") and half(s_["size"]) == 14.0}
ck.log("letters a, b, c present once at 14 pt, a and b on one baseline, each at its panel's top-left corner", len(letters) == 3 and abs(letters["a"][1] - letters["b"][1]) <= 1.0 and abs(letters["b"][0] - (14.0 + DX)) <= 1.0 and abs(letters["c"][0] - 14.0) <= 1.0, f"{letters}")
# raster: the new render is V27's render cut into the pieces and translated by whole pixels, white elsewhere
C.gs_render(SRC, f"{W}/v27_{DPI}.png", DPI); C.gs_render(NEW, f"{W}/new{TAG}_{DPI}.png", DPI)
A = C.load_png(f"{W}/v27_{DPI}.png"); B = C.load_png(f"{W}/new{TAG}_{DPI}.png")
E = np.full(B.shape, 255, int)
def px(v): return int(round(v * S))
E[0:px(REG["a"][3]), 0:A.shape[1]] = A[0:px(REG["a"][3]), 0:A.shape[1]]                                   # band + piece a rows, at place (V27 has only panel a and white there)
for k in ("b", "c"):
    r = REG[k]; dx, dy = OFF[k]; y0, y1 = px(r[1]), px(r[3]); x0 = px(dx); yy = px(r[1] + dy)
    E[yy:yy + (y1 - y0), x0:x0 + A.shape[1]] = A[y0:y1, 0:A.shape[1]]
dmax = np.abs(B - E).max(axis=2); diff = dmax > 0; n_bad = int(diff.sum()); worst = int(dmax.max()) if n_bad else 0
ck.log(f"render at {DPI} dpi equals V27's render cut at y = 418 and 818 and translated by whole pixels (every mark, rule, string and colour identical, nothing else on the page): at most 20 pixels may differ by at most 12 of 255 (Ghostscript's anti-aliasing of the band's transparency groups and of dashed rules on a raster of another size)", n_bad <= 20 and worst <= 12, f"{n_bad} differing pixels, largest difference {worst} of 255")
if n_bad:
    ys, xs = np.nonzero(diff); ck.info("  differing pixels (pt, expected, got): " + "; ".join(f"({x / S:.1f}, {y / S:.1f}) {E[y, x].tolist()} {B[y, x].tolist()}" for y, x in zip(ys, xs)).replace(";", ","))
sc = min(523.3 / PW, 683.7 / PH)
ck.info(f"assembly scale in the figures document (box 523.3 x 683.7 pt): {sc:.3f} (V27 skyscraper 0.604, main sheets 0.540): an 11 pt string prints at {11 * sc:.1f} pt")
C.gs_render(NEW, f"{D}/{SHEET}{TAG}_{DPI}dpi.png" if alt is None else f"{W}/{SHEET}{TAG}_{DPI}dpi.png", DPI)
json.dump({"sheet": SHEET, "source_v27": SRC, "source_sha256": C.sha256(SRC), "body": FLAT0, "body_sha256": C.sha256(FLAT0), "band": BAND, "band_sha256": C.sha256(BAND), "regions": REG, "offsets": OFF,
           "page_before": [PW0, PH0], "page_after": [PW, PH], "new": NEW, "new_sha256": C.sha256(NEW), "scale_in_document": round(sc, 4), "written": C.now()}, open(f"{W}/build_record{TAG}.json", "w"), indent=1)
ck.write(f"{V}/checks{TAG}.txt")
