#!$T90_PY
"""Supp_Fig21, round 49 (LSUPP-A): the page box step of the round-40 chain (MediaBox = CropBox = panel a plus the right margin, PyMuPDF box
mode on the simple matplotlib page), then the positive control of the round-40 finish: every marker, the dashed rule and the T90 ring of
panel a at the round-37 record's page coordinates (0.05 pt), the key ring inside the axes 19.62 pt above the floor, no panel letter
(the V16 letter box blank on the 150 dpi render). Writes verify/checks.txt (the round-49 finish step appends the census block)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import json, os, subprocess, sys
import fitz, numpy as np
from PIL import Image
L = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LSUPP-A"
S = "Supp_Fig21"; D = f"{L}/{S}"; W = f"{D}/work"; VER = f"{D}/verify"; NEW = f"{D}/{S}.pdf"; GS = paths.GS
DR = json.load(open(f"{VER}/{S}_drawn.json")); R37 = json.load(open(f"{W}/Supp_Fig21_drawn_R37.json"))
box = DR["new_page_box"]; W_NEW, H_NEW = round(box[2], 3), round(box[3], 3)
doc = fitz.open(f"{W}/Supp_Fig21_mpl.pdf"); assert len(doc) == 1; pg = doc[0]; r = fitz.Rect(0, 0, W_NEW, H_NEW); cur = pg.rect
if abs(r.x1 - cur.width) < 0.01: r.x1 = cur.width
if abs(r.y1 - cur.height) < 0.01: r.y1 = cur.height
pg.set_mediabox(r); pg.set_cropbox(pg.mediabox); doc.set_metadata({"title": S, "creator": "LSUPP-A round 49 r40_build_s21.py +1 pt (matplotlib) + box"})
doc.save(NEW, garbage=3, deflate=True); doc.close()
lines = []
def log(name, ok, detail=""):
    lines.append(f"{'PASS' if ok else 'FAIL'}  {name}" + (f": {detail}" if detail else "")); print(lines[-1][:300])
def pts(rec, kind): return [v for v in rec["values"] if v["panel"] == "a" and v["kind"] == kind]
for kind in ("marker", "open ring"):
    a, b = pts(DR, kind), pts(R37, kind)
    log(f"panel a {kind}s: {len(a)} drawn at the round-37 page coordinates within 0.05 pt, same values", len(a) == len(b) and all(abs(x["page_xy"][0] - y["page_xy"][0]) <= 0.05 and abs(x["page_xy"][1] - y["page_xy"][1]) <= 0.05 and abs(x["value"] - y["value"]) < 1e-9 for x, y in zip(a, b)), f"new {[x['page_xy'] for x in a]} r37 {[y['page_xy'] for y in b]}")
ra, rb = pts(DR, "dashed rule")[0], pts(R37, "dashed rule")[0]
log("panel a dashed rule at the round-37 page y within 0.05 pt, same T90 gain", abs(ra["page_y"] - rb["page_y"]) <= 0.05 and abs(ra["value"] - rb["value"]) < 1e-9, f"new {ra['page_y']} r37 {rb['page_y']}")
kr = pts(DR, "key ring")[0]
log("the key ring sits inside panel a's axes, 19.62 pt above the floor as on V16, its label clear of the curve, markers, rule label and spine (builder asserts at 10.5 pt)", DR["axes_a"][0] < kr["page_xy"][0] < DR["axes_a"][2] and abs((DR["axes_a"][3] - kr["page_xy"][1]) - (DR["axes_b_v16"][3] - kr["v16_page_xy"][1])) < 0.05, f"ring {kr['page_xy']}, label box {DR['key_box_pt']}, V16 ring {kr['v16_page_xy']}")
log("text sizes drawn: axis titles 12, tick labels 11, rule and key labels 10.5 (V26 11, 10, 9.5)", DR["text_pt"] == dict(axis_title=12.0, tick=11.0, annotation=10.5), str(DR["text_pt"]))
png = f"{D}/{S}_150dpi.png"
subprocess.run([GS, "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=png16m", "-r150", "-dTextAlphaBits=4", "-dGraphicsAlphaBits=4", f"-sOutputFile={png}", NEW], check=True)
im = np.asarray(Image.open(png).convert("L")); z = 150 / 72; lb = [4.32, 0.5, 11.55, 18.52]
sub = im[int(lb[1] * z):int(lb[3] * z) + 1, int(lb[0] * z):int(lb[2] * z) + 1]
log("no panel letter: the V16 letter a box (4.32 to 11.55 x, 0.5 to 18.52 y) is blank white on the render", int((sub < 250).sum()) == 0, f"{int((sub < 250).sum())} non-white px of {sub.size}")
d2 = fitz.open(NEW); log("page box set to panel a plus the right margin (MediaBox = CropBox)", abs(d2[0].rect.width - W_NEW) < 0.01 and abs(d2[0].rect.height - H_NEW) < 0.01 and len(d2) == 1, f"{d2[0].rect}"); d2.close()
ok = all(l.startswith("PASS") for l in lines)
open(f"{VER}/checks.txt", "w").write("\n".join([f"{S} round 49 (LSUPP-A): +1 pt on the round-40 build, positive control against the round-37 record", ""] + lines + ["", "RESULT ALL PASS" if ok else "RESULT FAIL"]) + "\n")
print("RESULT ALL PASS" if ok else "RESULT FAIL"); sys.exit(0 if ok else 1)
