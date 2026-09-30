#!$T90_PY
"""Lane LSUPP-D, round 52 (2026-09-26), brief item 2: Supp_Fig16 as in round 50 (a and b side by side, c below) with the sketch band
moved from the left edge to the middle of the page (horizontally centred, whole 150 dpi pixels). Same pieces, same offsets for a, b
and c as the round-50 build (relayout_s16.py), the band translated by BDX. Verified against V27 (the round-50 logic, band offset
added) and against V29 = the round-50 sheet (identical outside the band rows, the band rows a pixel translation).
usage: relayout_s16_r52.py"""
import json, os, shutil, sys
from collections import Counter
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lsd_common as C

SHEET = "Supp_Fig16"; D = f"{C.LANE}/{SHEET}"; W = f"{D}/work"; V = f"{D}/verify"
for d_ in (W, V): os.makedirs(d_, exist_ok=True)
SRC = C.hydrated(f"{C.V27}/{SHEET}.pdf"); FLAT0 = C.hydrated(f"{C.R49LANE}/{SHEET}/work/{SHEET}_restamped.pdf"); BAND = C.hydrated(C.BAND_R49)
PREV = C.hydrated(f"{C.V29}/{SHEET}.pdf"); PREV_REC = json.load(open(f"{C.R50LANE}/{SHEET}/work/build_record.json"))
assert C.sha256(SRC).startswith("6e8eb1c3172d7b35"), C.sha256(SRC)            # V27 Supp_Fig16 = the round-49 15:32 build
assert C.sha256(FLAT0).startswith("26a4b3d4f0f01a91"), C.sha256(FLAT0)         # its body form (round-49 restamped flat, band-free)
assert C.sha256(BAND).startswith("13e8b4b2042422cb"), C.sha256(BAND)           # the band of the 15:32 rebuild (title 11 pt)
assert C.sha256(PREV) == PREV_REC["new_sha256"], "V29 Supp_Fig16 is not the round-50 lane output"
PW0, PH0 = 496.8, 1131.7; BAND_H = 79.0
DPI = 150; S = DPI / 72.0; PX = 72.0 / DPI                                     # 0.48 pt = one pixel at 150 dpi
REG = {"a": [0.0, BAND_H, 285.0, 418.0], "b": [0.0, 418.0, PW0, 818.0], "c": [0.0, 818.0, PW0, PH0]}
A_RIGHT, GUTTER, B_LEFT_MARGIN = 271.5, 32.2, 14.0
DY = -round((429.8 - 92.4) / PX) * PX                                          # -337.44, as in round 50
DX = round((A_RIGHT + GUTTER - B_LEFT_MARGIN) / PX) * PX                       # 289.44, as in round 50
PW, PH = round(DX + PW0, 2), round(PH0 + DY, 2)
assert [PW, PH] == PREV_REC["page_after"] and PREV_REC["offsets"]["b"] == [DX, DY], (PW, PH, PREV_REC["page_after"])
nb_, bw, bh = C.page_box(BAND); assert abs(bw - PW0) < 0.05 and abs(bh - BAND_H) < 0.05, (bw, bh)
BDX = round(((PW - bw) / 2) / PX) * PX                                         # the band's translation: centred on the page, whole pixels (144.96)
NEW = f"{D}/{SHEET}.pdf"
ck = C.Checks(f"{SHEET}, lane LSUPP-D round 52 (item 2, {C.now()}): band centred on the page (dx {BDX} pt), a and b side by side, c below, as round 50. Source V27 {SRC} (sha256 {C.sha256(SRC, 16)}), body form {FLAT0} (sha256 {C.sha256(FLAT0, 16)}), band {BAND} (sha256 {C.sha256(BAND, 16)}), previous sheet V29 {PREV} (sha256 {C.sha256(PREV, 16)}). Page {PW} x {PH} pt (unchanged from round 50), dx {DX} dy {DY}, band dx {BDX} (whole 150 dpi pixels).")
FLAT = f"{W}/{SHEET}_body.pdf"; shutil.copyfile(FLAT0, FLAT)
# ------------------------------------------------------------------ 1 text-exact pieces (as round 50)
OFF = {"a": (0.0, 0.0), "b": (DX, DY), "c": (0.0, DY), "band": (BDX, 0.0)}
PIECES = {}
for k, r in REG.items():
    comp = [[0, 0, PW0, r[1]], [0, r[3], PW0, PH0]] + ([[r[2], r[1], PW0, r[3]]] if r[2] < PW0 else [])
    out = f"{W}/piece_{k}.pdf"; job = {"src": FLAT, "out": out, "whole_page": False, "redact_rects": comp, "items": []}
    json.dump(job, open(f"{W}/piece_{k}_job.json", "w")); C.wd(f"{SHEET}_piece_{k}", [C.PY, f"{C.SCRIPTS}/fitz_restamp.py", f"{W}/piece_{k}_job.json"])
    runs = C.gs_txtwrite(out, f"{W}/piece_{k}_census.xml"); PIECES[k] = (out, runs)
    inside = all(r[0] - 1.5 <= s_["x0"] and s_["x1"] <= r[2] + 1.5 and r[1] <= s_["y"] <= r[3] for s_ in runs)
    ck.log(f"piece {k}: every remaining run lies inside its region {r} ({len(runs)} runs)", inside, f"outside: {[(s_['text'], s_['x0'], s_['y']) for s_ in runs if not (r[0] - 1.5 <= s_['x0'] and s_['x1'] <= r[2] + 1.5 and r[1] <= s_['y'] <= r[3])][:5]}")
# ------------------------------------------------------------------ 2 compose
parts = []
for k in ("a", "b", "c"):
    r = REG[k]; dx, dy = OFF[k]
    parts.append({"pdf": PIECES[k][0], "rect": [r[0] + dx, r[1] + dy, r[2] + dx, r[3] + dy], "clip": r})
parts.append({"pdf": BAND, "rect": [BDX, 0, BDX + bw, bh]})
job = {"out": NEW, "page": [PW, PH], "parts": parts, "metadata": {"title": SHEET, "creator": f"LSUPP-D r52 relayout_s16_r52.py: V27 body pieces (a at place, b dx {DX} dy {DY}, c dy {DY}) + band_s16_r49 (13e8b4b2) centred at dx {BDX}"}}
json.dump(job, open(f"{W}/compose_job.json", "w")); C.wd(f"{SHEET}_compose", [C.PY, f"{C.SCRIPTS}/fitz_compose.py", f"{W}/compose_job.json"])
# ------------------------------------------------------------------ 3 verify
nn, nw, nh = C.page_box(NEW)
ck.log("one page, page box unchanged from round 50", nn == 1 and abs(nw - PW) < 0.05 and abs(nh - PH) < 0.05, f"{nw:.2f} x {nh:.2f} pt")
ck.log("the band is centred on the page (its centre within half a pixel of the page's centre)", abs((BDX + bw / 2) - PW / 2) <= PX / 2, f"band x {BDX:.2f} to {BDX + bw:.2f} of {PW}, centre offset {(BDX + bw / 2) - PW / 2:.3f} pt")
half = lambda v: round(v * 2) / 2
ref = C.gs_txtwrite(SRC, f"{W}/census_ref.xml"); new = C.gs_txtwrite(NEW, f"{W}/census_new.xml")
cref = Counter((s_["text"], half(s_["size"])) for s_ in ref); cnew = Counter((s_["text"], half(s_["size"])) for s_ in new)
ck.log(f"txtwrite census equal to V27's (and so to V29's): same multiset of strings and sizes ({sum(cnew.values())} runs, band included)", cref == cnew, f"difference {dict((cref - cnew) + (cnew - cref))}")
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
ck.log("every V27 run reappears at its piece's translation within 1 pt (a at place, b by (dx, dy), c by dy, the band by its dx)", not unmatched and all(used), f"unmatched {unmatched[:6]}")
band_own = Counter((s_["text"], half(s_["size"])) for s_ in C.gs_txtwrite(BAND, f"{W}/census_band.xml"))
band_new = Counter((s_["text"], half(s_["size"])) for s_ in new if s_["y"] < BAND_H)
ck.log("band region census equals band_s16_r49's own census", band_new == band_own, f"{sum(band_new.values())} runs")
letters = {s_["text"]: (s_["x0"], s_["y"]) for s_ in new if s_["text"] in ("a", "b", "c") and half(s_["size"]) == 14.0}
ck.log("letters a, b, c present once at 14 pt, a and b on one baseline, each at its panel's top-left corner (unchanged)", len(letters) == 3 and abs(letters["a"][1] - letters["b"][1]) <= 1.0 and abs(letters["b"][0] - (14.0 + DX)) <= 1.0 and abs(letters["c"][0] - 14.0) <= 1.0, f"{letters}")
# raster against V29 (the round-50 sheet): identical outside the band rows, the band rows translated by whole pixels
C.gs_render(PREV, f"{W}/v29_{DPI}.png", DPI); C.gs_render(NEW, f"{W}/new_{DPI}.png", DPI)
P = C.load_png(f"{W}/v29_{DPI}.png"); B = C.load_png(f"{W}/new_{DPI}.png")
ck.log("render size equals V29's", P.shape == B.shape, f"{B.shape[1]} x {B.shape[0]} px")
def px(v): return int(round(v * S))
bh_px = px(BAND_H); bdx_px = px(BDX); bw_px = px(bw)
E = P.copy(); E[0:bh_px, :, :] = 255; E[0:bh_px, bdx_px:bdx_px + bw_px] = P[0:bh_px, 0:bw_px]
dmax = np.abs(B - E).max(axis=2); diff = dmax > 0; n_bad = int(diff.sum()); worst = int(dmax.max()) if n_bad else 0
below = int(diff[bh_px:, :].sum())
ck.log(f"render at {DPI} dpi equals V29's with the band rows translated by {bdx_px} pixels (every mark, rule, string and colour identical): at most 40 pixels may differ by at most 16 of 255 (anti-aliasing at the band's edges and the row boundary), none below the band rows", n_bad <= 40 and worst <= 16 and below == 0, f"{n_bad} differing pixels, largest difference {worst} of 255, {below} below the band rows")
if n_bad:
    ys, xs = np.nonzero(diff); ck.info("  differing pixels (pt, expected, got): " + "; ".join(f"({x / S:.1f}, {y / S:.1f}) {E[y, x].tolist()} {B[y, x].tolist()}" for y, x in list(zip(ys, xs))[:40]).replace(";", ","))
sc = min(523.3 / PW, 683.7 / PH)
ck.info(f"assembly scale in the figures document (box 523.3 x 683.7 pt): {sc:.3f}, unchanged from round 50")
C.gs_render(NEW, f"{D}/{SHEET}_{DPI}dpi.png", DPI)
json.dump({"sheet": SHEET, "source_v27": SRC, "source_sha256": C.sha256(SRC), "body": FLAT0, "body_sha256": C.sha256(FLAT0), "band": BAND, "band_sha256": C.sha256(BAND), "previous_v29": PREV, "previous_sha256": C.sha256(PREV),
           "regions": REG, "offsets": OFF, "band_dx": BDX, "page": [PW, PH], "new": NEW, "new_sha256": C.sha256(NEW), "scale_in_document": round(sc, 4), "written": C.now()}, open(f"{W}/build_record.json", "w"), indent=1)
ck.write(f"{V}/checks.txt")
