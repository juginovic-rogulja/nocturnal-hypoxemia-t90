#!$T90_PY
"""Round 49, lane L5: (1) prove the copied generator reproduces the Aug-25 panel (FINAL_FIGURES_2026-08-14/Main_Figures_Polished/Figure5/
Figure5A.pdf): same strings at the same origins and sizes, same drawing rectangles (PyMuPDF on these small flat pages); (2) fit the map
from the generator page to the V13 sheet from the cached V13 text probe (work/base_text.json, the round-37 read of the nested V13 sheet):
X = S * x + TX, Y = S * y + TY over the 15 panel a strings; (3) write work/panel_a/map_v13.json. usage: 04b_panel_a_map_r49.py <bump0 pdf>"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import json, sys, os
import fitz, numpy as np
L5 = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/L5"; W = f"{L5}/Main_Fig5/work"
OLD = f"{paths.FIGURE_ROOT}/FINAL_FIGURES_2026-08-14/Main_Figures_Polished/Figure5/Figure5A.pdf"
NEW = sys.argv[1] if len(sys.argv) > 1 else f"{W}/panel_a/Figure5A_bump0.pdf"


def spans(pdf):
    d = fitz.open(pdf); assert len(d) == 1 and len(d[0].get_xobjects()) == 0, ("flat page expected", pdf); out = []
    for b in d[0].get_text("rawdict")["blocks"]:
        if b["type"] != 0: continue
        for l in b["lines"]:
            for s in l["spans"]:
                t = "".join(c["c"] for c in s["chars"]).replace("\xa0", " ").strip()
                if t: out.append(dict(text=t, x=round(s["chars"][0]["origin"][0], 3), y=round(s["chars"][0]["origin"][1], 3), size=round(s["size"], 3), font=s["font"], color="%06x" % s["color"]))
    rect = [round(v, 3) for v in d[0].rect]; dr = [[round(v, 3) for v in it["rect"]] + [it.get("fill"), it.get("color"), it.get("width")] for it in d[0].get_drawings()]
    d.close(); return out, rect, dr


so, ro, do = spans(OLD); sn, rn, dn = spans(NEW)
key = lambda s: (s["text"], round(s["x"], 2), round(s["y"], 2), round(s["size"], 2), s["font"], s["color"])
same_text = sorted(map(key, so)) == sorted(map(key, sn)); same_rect = ro == rn
def drkey(d): return (tuple(round(v, 2) for v in d[:4]), str(d[5]), d[6])          # geometry, stroke and width (the fills are remapped, see below)
same_draw = sorted(map(drkey, do)) == sorted(map(drkey, dn))
def hx(c): return None if c is None else "#%02x%02x%02x" % tuple(int(round(v * 255)) for v in c[:3])
REMAP = {"#ccd1d6": "#b3dcf2", "#8a9099": "#7cc0e9", "#4d7387": "#3f9fd8", "#1f4257": "#0288d1"}   # the Aug-25 ramp -> the sheet's clinical palette (measured on V26)
fills_ok = sorted((tuple(round(v, 2) for v in d[:4]), REMAP.get(hx(d[4]), hx(d[4]))) for d in do) == sorted((tuple(round(v, 2) for v in d[:4]), hx(d[4])) for d in dn)
print(f"old page {ro}, new page {rn}: same {same_rect}; spans old {len(so)} new {len(sn)}: identical (text, origin 2dp, size, font, colour) {same_text}; drawings old {len(do)} new {len(dn)}: identical rects/strokes/widths {same_draw}; fills = the Aug-25 fills under the remap {REMAP}: {fills_ok}")
same_draw = same_draw and fills_ok
if not same_text:
    for a, b in zip(sorted(map(key, so)), sorted(map(key, sn))):
        if a != b: print("  DIFF", a, b)
if not same_draw:
    for a, b in zip(sorted(map(drkey, do)), sorted(map(drkey, dn))):
        if a != b: print("  DRAW DIFF", a, b)
# the V13 map
J = json.load(open(f"{W}/base_text.json")); V = [s for s in J["spans"] if s["x"] < 497 and 16 < s["y"] < 393 and not (len(s["text"]) == 1 and s["text"] == "a")]
pairs = []
for v in V:
    m = [s for s in sn if s["text"].replace("\xa0", " ") == v["text"]]
    assert len(m) == 1, (v["text"], m); pairs.append((m[0], v))
G = np.array([[p[0]["x"], p[0]["y"]] for p in pairs]); Vv = np.array([[p[1]["x"], p[1]["y"]] for p in pairs])
# least squares for S, TX, TY: X = S x + TX, Y = S y + TY
A = np.zeros((2 * len(pairs), 3)); bvec = np.zeros(2 * len(pairs))
A[0::2, 0] = G[:, 0]; A[0::2, 1] = 1; bvec[0::2] = Vv[:, 0]; A[1::2, 0] = G[:, 1]; A[1::2, 2] = 1; bvec[1::2] = Vv[:, 1]
(S, TX, TY), *_ = np.linalg.lstsq(A, bvec, rcond=None); res = A @ np.array([S, TX, TY]) - bvec
sizes = sorted(set((p[0]["size"], p[1]["size"]) for p in pairs))
print(f"V13 map: S = {S:.5f}, TX = {TX:.3f}, TY = {TY:.3f}; max residual {np.abs(res).max():.3f} pt over {len(pairs)} strings; sizes generator -> V13 {sizes} (S x size = {[round(a * S, 3) for a, _ in sizes]})")
for (g, v), r in zip(pairs, res.reshape(-1, 2)): print(f"   {v['text']:42s} gen ({g['x']:7.2f},{g['y']:7.2f}) {g['size']:5.2f} -> V13 ({v['x']:7.2f},{v['y']:7.2f}) {v['size']:5.2f}  res ({r[0]:+.3f},{r[1]:+.3f})")
gen_rect = rn; tgt = [S * gen_rect[0] + TX, S * gen_rect[1] + TY, S * gen_rect[2] + TX, S * gen_rect[3] + TY]
print(f"generator page {gen_rect} maps to the sheet rect {[round(v, 3) for v in tgt]}")
json.dump(dict(S=S, TX=TX, TY=TY, max_residual=float(np.abs(res).max()), n=len(pairs), generator_page=gen_rect, target_rect=tgt, sizes_gen_to_v13=sizes, bump_on_canvas_for_plus_1=1.0 / S,
               old_pdf=OLD, new_pdf=NEW, generator_reproduces_aug25=bool(same_text and same_rect and same_draw), fill_remap=REMAP, n_spans=len(sn), n_drawings=len(dn)), open(f"{W}/panel_a/map_v13.json", "w"), indent=1)
print("bump on the canvas for +1.0 pt on the sheet:", 1.0 / S)
