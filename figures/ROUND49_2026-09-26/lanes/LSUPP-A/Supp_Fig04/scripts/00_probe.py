#!/usr/bin/env python3
"""Supp_Fig04, V14 lane L2, step 0: probe the V13 sheet (repointed copy of the round-30 probe, byte copy beside as 00_probe_PRE_V8_1.py).
V13 = the round-30 re-plot with the round-31 LTEXT rebuild of panels c and d (bigger axes, 2 Form XObjects, no hidden panel-d copy).
Writes work/base_text.json, work/base_drawings.json, work/base_geometry.json, work/base_150.png."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import collections, fcntl, gc, json, os, sys
import numpy as np, fitz
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND37_2026-09-14/figures/_common")
from v14lib import V13, hydrated, sha256, LOCK_ROOT
LANE = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); WORK = os.path.join(LANE, "work"); os.makedirs(WORK, exist_ok=True)
BASE = f"{V13}/Supp_Fig04.pdf"; hydrated(BASE)
def hx(c): return None if c is None else "#%02x%02x%02x" % tuple(int(round(v * 255)) for v in c)
doc = fitz.open(BASE); page = doc[0]; W, H = page.rect.width, page.rect.height
print("page", W, H, "rotation", page.rotation, "xobjects", len(page.get_xobjects()))
spans = []
for b in page.get_text("rawdict")["blocks"]:
    if b["type"] != 0: continue
    for l in b["lines"]:
        for s in l["spans"]:
            txt = "".join(ch["c"] for ch in s["chars"]); org = s["chars"][0]["origin"]
            spans.append(dict(text=txt, font=s["font"], size=round(s["size"], 3), color=s["color"], flags=s["flags"], asc=round(s["ascender"], 4), desc=round(s["descender"], 4),
                              bbox=[round(v, 3) for v in s["bbox"]], origin=[round(org[0], 3), round(org[1], 3)], dir=[round(v, 3) for v in l["dir"]]))
words = [dict(text=w[4], bbox=[round(v, 3) for v in w[:4]]) for w in page.get_text("words")]
json.dump(dict(page=[W, H], spans=spans, words=words, base_sha256=sha256(BASE)), open(os.path.join(WORK, "base_text.json"), "w"), indent=0)
font_desc = collections.Counter((s["font"], s["asc"], s["desc"]) for s in spans); print("font descriptors:", font_desc.most_common())
dr = page.get_drawings(); recs = []
for d in dr:
    r = d["rect"]; items = d.get("items", [])
    recs.append(dict(type=d.get("type"), fill=None if d.get("fill") is None else [round(v, 4) for v in d["fill"]], color=None if d.get("color") is None else [round(v, 4) for v in d["color"]], width=d.get("width"),
                     dashes=d.get("dashes"), rect=[round(r.x0, 3), round(r.y0, 3), round(r.x1, 3), round(r.y1, 3)], w=round(r.width, 3), h=round(r.height, 3), n_items=len(items),
                     kinds=dict(collections.Counter(it[0] for it in items)), seqno=d.get("seqno")))
json.dump(recs, open(os.path.join(WORK, "base_drawings.json"), "w"), indent=0)
G = dict(page=[W, H], base_sha256=sha256(BASE), n_xobjects=len(page.get_xobjects()))
BLUE, INK, WHITE, BAND = "#0288d1", "#1a1d21", "#ffffff", "#eef0f1"
whites = [r for r in recs if r["fill"] and hx(r["fill"]) == WHITE and r["w"] > 100 and r["h"] > 100]
AX = dict(a=[179.6, 44.241, 485.764, 461.265], b=[688.324, 44.24, 994.488, 617.648], c=[84.56, 735.248, 485.764, 913.648], d=[603.4, 735.248, 994.488, 913.648])
for k, box in AX.items():
    hit = [r for r in whites if all(abs(a - b) < 0.05 for a, b in zip(r["rect"], box))]; assert hit, (k, box, [r["rect"] for r in whites])
G["axes_box"] = AX
bands = sorted([r for r in recs if r["fill"] and hx(r["fill"]) == BAND and r["w"] > 100], key=lambda r: (r["rect"][0], r["rect"][1]))
G["row_bands"] = dict(a=sorted({(r["rect"][1], r["rect"][3]) for r in bands if abs(r["rect"][0] - AX["a"][0]) < 0.05}), b=sorted({(r["rect"][1], r["rect"][3]) for r in bands if abs(r["rect"][0] - AX["b"][0]) < 0.05}))
pitch_a = np.diff([b[0] for b in G["row_bands"]["a"]]); pitch_b = np.diff([b[0] for b in G["row_bands"]["b"]])
G["row_pitch"] = dict(a=round(float(np.median(pitch_a)) / 2, 4), b=round(float(np.median(pitch_b)) / 2, 4), rows_a=16, rows_b=22)
assert abs(G["row_pitch"]["a"] - (AX["a"][3] - AX["a"][1]) / 16) < 0.01 and abs(G["row_pitch"]["b"] - (AX["b"][3] - AX["b"][1]) / 22) < 0.01, G["row_pitch"]
G["hr1_rules"] = [dict(x=r["rect"][0], y0=r["rect"][1], y1=r["rect"][3], width=r["width"], dashes=r["dashes"], seq=r["seqno"]) for r in recs if r["dashes"] and r["dashes"] not in ("[] 0", "[]0")]
G["axis_rules"] = [dict(x0=r["rect"][0], x1=r["rect"][2], y=r["rect"][1], width=r["width"], seq=r["seqno"]) for r in recs if r["color"] and hx(r["color"]) == INK and r["w"] > 100 and r["h"] < 0.01]
ticks = [r for r in recs if r["color"] and hx(r["color"]) == INK and abs((r["width"] or 0) - 0.8) < 0.01 and 2.9 < max(r["w"], r["h"]) < 3.1 and min(r["w"], r["h"]) < 0.01]
G["ticks_x"] = sorted({(round(r["rect"][0], 3), round(r["rect"][1], 3), round(r["rect"][3], 3)) for r in ticks if r["h"] > r["w"]})
G["ticks_y"] = sorted({(round(r["rect"][1], 3), round(r["rect"][0], 3), round(r["rect"][2], 3)) for r in ticks if r["w"] > r["h"]})
mk = [r for r in recs if r["fill"] is not None and 0 < r["w"] < 10 and 0 < r["h"] < 10]
G["markers"] = [dict(cx=round((r["rect"][0] + r["rect"][2]) / 2, 3), cy=round((r["rect"][1] + r["rect"][3]) / 2, 3), w=r["w"], h=r["h"], fill=hx(r["fill"]), edge=hx(r["color"]), edge_w=r["width"], kinds=r["kinds"], seq=r["seqno"]) for r in mk]
ci = [r for r in recs if r["color"] and abs((r["width"] or 0) - 1.7) < 0.01 and r["h"] < 0.01]
G["ci_lines"] = [dict(x0=r["rect"][0], x1=r["rect"][2], y=r["rect"][1], colour=hx(r["color"]), seq=r["seqno"]) for r in ci]
inbox = lambda r, box: r["rect"][0] >= box[0] - 0.5 and r["rect"][2] <= box[2] + 0.5 and r["rect"][1] >= box[1] - 0.5 and r["rect"][3] <= box[3] + 0.5
G["c_bars"] = [dict(rect=r["rect"], fill=hx(r["fill"]), seq=r["seqno"]) for r in recs if r["fill"] and hx(r["fill"]) not in (WHITE, BLUE) and inbox(r, AX["c"]) and r["w"] > 5 and r["h"] > 1]
G["c_line"] = [dict(rect=r["rect"], colour=hx(r["color"]), width=r["width"], seq=r["seqno"]) for r in recs if r["color"] and abs((r["width"] or 0) - 1.8) < 0.01]
G["c_markers"] = [m for m in G["markers"] if AX["c"][0] <= m["cx"] <= AX["c"][2] and AX["c"][1] <= m["cy"] <= AX["c"][3]]
G["d_bars"] = [dict(rect=r["rect"], fill=hx(r["fill"]), seq=r["seqno"]) for r in recs if r["fill"] and hx(r["fill"]) != WHITE and inbox(r, AX["d"]) and r["h"] > 5]
G["key_a"] = dict(markers=[m for m in G["markers"] if 505 < m["cy"] < 545 and m["cx"] < 320], lines=[c for c in G["ci_lines"] if 505 < c["y"] < 545 and c["x0"] < 320])
G["key_b"] = dict(markers=[m for m in G["markers"] if m["cx"] > 870 and m["cy"] < 140], lines=[c for c in G["ci_lines"] if c["x0"] > 870 and c["y"] < 140])
G["colours_vector"] = dict(fills=sorted(set(hx(r["fill"]) for r in recs if r["fill"])), strokes=sorted(set(hx(r["color"]) for r in recs if r["color"])))
G["font_descriptors"] = [dict(font=k[0], asc=k[1], desc=k[2], spans=v) for k, v in font_desc.most_common()]
G["letters_text_layer"] = [dict(text=s["text"], x=s["origin"][0], baseline=s["origin"][1], top=s["bbox"][1], bottom=s["bbox"][3], font=s["font"], size=s["size"]) for s in spans if len(s["text"]) == 1 and s["text"] in "abcd" and s["size"] > 12]
pix = page.get_pixmap(dpi=150, colorspace=fitz.csRGB, alpha=False); pix.save(os.path.join(WORK, "base_150.png"))
A = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.h, pix.w, 3).copy(); pix = None; s = 150 / 72.0
def ink_in(box, thr=200):
    x0, y0, x1, y1 = [int(round(v * s)) for v in box]; sub = A[y0:y1, x0:x1]; return int((sub.min(axis=2) < thr).sum())
for L in G["letters_text_layer"]:
    L["ink_px_150dpi"] = ink_in((L["x"] - 0.5, L["top"], L["x"] + 8.5, L["bottom"])); L["visible"] = L["ink_px_150dpi"] > 20
G["c_annotation_span_colour"] = ["#%06x" % s["color"] for s in spans if s["text"] in ("32/32", "30/32", "19/32")][:1]
G["no_estimate_span_colour"] = ["#%06x" % s["color"] for s in spans if s["text"] == "No estimate"][:1]
json.dump(G, open(os.path.join(WORK, "base_geometry.json"), "w"), indent=1, default=float)
doc.close(); gc.collect()
print("axes boxes", G["axes_box"]); print("row pitch", G["row_pitch"]); print("letters:", G["letters_text_layer"]); print("key a", len(G["key_a"]["markers"]), len(G["key_a"]["lines"]), "key b", len(G["key_b"]["markers"]), len(G["key_b"]["lines"]))
print("c bars", len(G["c_bars"]), "c markers", len(G["c_markers"]), "d bars", len(G["d_bars"]), "markers", len(G["markers"]), "ci lines", len(G["ci_lines"]))
