#!$T90_PY
"""Verifier for ED_Fig07 round 54: Ghostscript census = V30's spans (same text, size, position) + the round-49 panel c spans translated by (0, Y0)
+ one letter c at 14 pt; page box declared; c's ink clear of panel b's ink and below panel a's; the 150-dpi render identical to V30's on the
top 539.98 pt outside the letter's box; nothing below 9 pt (the Extended Data gate); no em dash, semicolon, banned word."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import json, re, html, subprocess, collections, os
import numpy as np
from PIL import Image
T90 = paths.FIGURE_ROOT; LANE = f"{paths.FIGURE_ROOT}/ROUND54_2026-09-28/lanes/LED7C"; SD = f"{LANE}/ED_Fig07"; W_ = f"{SD}/work"; V_ = f"{SD}/verify"; GS = paths.GS
os.makedirs(V_, exist_ok=True); log = json.load(open(f"{W_}/compose_log.json")); NEW = log["out"]; V30 = log["v30"]; PC = log["panel_c"]; Y0 = log["panel_c_rect"][1]
lines = []; ok = True
def check(name, good, detail=""):
    global ok; ok &= bool(good); lines.append(f"{'PASS' if good else 'FAIL'} {name}" + (f": {detail}" if detail else "")); print(lines[-1])
def census(pdf, xml):
    subprocess.run([GS, "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=txtwrite", "-dTextFormat=0", f"-sOutputFile={xml}", pdf], check=True, capture_output=True, timeout=300)
    raw = open(xml, encoding="utf8", errors="replace").read(); out = []
    for m in re.finditer(r'<span bbox="(\S+) (\S+) (\S+) (\S+)" font="([^"]*)" size="([^"]*)">(.*?)</span>', raw, re.S):
        chars = re.findall(r'<char bbox="(\S+) (\S+) (\S+) (\S+)" c="(.*?)"/>', m.group(7)); text = "".join(html.unescape(c[4]) for c in chars).strip()
        if text: out.append(dict(text=text, x0=round(float(m.group(1)), 1), y0=round(float(m.group(2)), 1), x1=round(float(m.group(3)), 1), y1=round(float(m.group(4)), 1), size=round(float(m.group(6)) * 2) / 2))
    assert out, f"empty census {pdf}"; return out
new = census(NEW, f"{V_}/new_txt.xml"); v30 = census(V30, f"{V_}/v30_txt.xml"); pc = census(PC, f"{V_}/panel_c_txt.xml")
exp = [(s["text"], s["size"], s["x0"], s["y0"]) for s in v30] + [(s["text"], s["size"], s["x0"], s["y0"] + Y0) for s in pc] + [("c", 14.0, 13.0, Y0 + 13.0)]
pool = [(s["text"], s["size"], s["x0"], s["y0"]) for s in new]; used = [False] * len(pool); missing = []
for e in exp:
    hit = [j for j, g in enumerate(pool) if not used[j] and g[0] == e[0] and g[1] == e[1] and abs(g[2] - e[2]) <= 1.5 and abs(g[3] - e[3]) <= 1.5]
    if hit: used[hit[0]] = True
    else: missing.append(e)
extra = [g for j, g in enumerate(pool) if not used[j]]
check("census = V30's spans + the round-49 panel c spans translated by (0, Y0) + the letter c at 14 pt (text, size, position within 1.5 pt)", not missing and not extra, f"missing {missing[:5]}, extra {extra[:5]}")
from pypdf import PdfReader
b = PdfReader(NEW).pages[0].mediabox; w, h = float(b.width), float(b.height)
check("page box declared: width = V30's, height = Y0 + clip + 8", abs(w - log["page"][0]) < 0.05 and abs(h - log["page"][1]) < 0.05 and abs(w - 894.61) < 0.05, f"{w:.2f} x {h:.2f} (V30 894.61 x 539.98)")
lt = [s for s in new if s["size"] == 14.0 and s["text"] in ("a", "b", "c")]; letter_c = [s for s in lt if s["text"] == "c"]
cink = [s for s in new if s["y0"] > Y0 and s["x0"] < 480 and s not in letter_c]; bink = [s for s in new if s["x0"] >= 484 and s["y0"] > 40]; aink = [s for s in new if s["x0"] < 484 and 40 < s["y0"] <= 365]
c_top = min(s["y0"] - 0.75 * s["size"] for s in cink); a_bot = max(s["y1"] for s in aink)
check("panel c's ink lies under panel a's (its first cap top at least 14 pt below a's last ink), the letter c at least 10 pt below, and left of panel b's ink by at least 12 pt", c_top >= a_bot + 14 and (letter_c[0]["y0"] - 10) >= a_bot + 10 and max(s["x1"] for s in cink) + 12 <= min(s["x0"] for s in bink), f"c cap top {c_top:.0f} and letter cap top {letter_c[0]['y0'] - 10:.0f} vs a bottom {a_bot:.0f}, c right {max(s['x1'] for s in cink):.0f} vs b left {min(s['x0'] for s in bink):.0f}")
check("letters a, b, c once each at 14 pt, a and c on one x (13), b at its column (498)", sorted(s["text"] for s in lt) == ["a", "b", "c"] and all(abs(s["x0"] - (498 if s["text"] == "b" else 13)) <= 1 for s in lt), str([(s["text"], s["x0"], s["y0"]) for s in lt]))
check("no text below 9 pt (Extended Data gate, V30's own 9.5-pt keeps unchanged), panel c at 11 and 12 pt", min(s["size"] for s in new) >= 9.0 and set(s["size"] for s in cink) <= {11.0, 12.0}, f"sizes {sorted(set(s['size'] for s in new))}")
def render(pdf, png): subprocess.run([GS, "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=png16m", "-r150", f"-sOutputFile={png}", pdf], check=True, capture_output=True, timeout=300); return np.asarray(Image.open(png).convert("RGB")).astype(int)
A = render(V30, f"{V_}/v30_150.png"); B = render(NEW, f"{SD}/ED_Fig07_150dpi.png"); S = 150 / 72
rows = int(539.98 * S); top = B[:rows, :A.shape[1]]; mask = np.ones(top.shape[:2], bool); mask[int(Y0 * S):, 0:int(497.2 * S) + 1] = False   # the region panel c and its letter now occupy (left column under panel a), white on V30
d = np.abs(A[:rows] - top).max(axis=2); strong = (d > 64) & mask
core = strong.copy()
for dy, dx in ((1, 0), (-1, 0), (0, 1), (0, -1)): core &= np.roll(strong, (dy, dx), axis=(0, 1))       # a differing pixel whose four neighbours also differ = a real change, not a glyph-edge phase
mean_abs = float(np.abs(A[:rows] - top)[mask].mean())
check("150-dpi render: the top 539.98 pt equals V30's outside panel c's region up to glyph-edge phase (no differing pixel survives a 3 x 3 erosion, mean abs difference under 2 of 255)", A.shape[1] == B.shape[1] and int(core.sum()) == 0 and mean_abs < 2.0, f"{int((d > 0).sum())} differing pixels, {int(core.sum())} surviving erosion, mean abs {mean_abs:.2f}")
txt = " ".join(s["text"] for s in new); bad = [w for w in ("honest", "straightforward", "prespecified", "pre-specified", "printed") if w in txt.lower()] + (["em dash"] if "—" in txt else []) + (["semicolon"] if ";" in txt else [])
check("no em dash, semicolon or banned word", not bad, str(bad))
lines.append("RESULT ALL PASS" if ok else "RESULT FAIL"); open(f"{V_}/checks.txt", "w").write(f"ED_Fig07 round 54 lane LED7C (option A): V30 sheet {os.path.basename(V30)} sha256 {log['v30_sha256'][:16]} + round-49 Figure 4 panel c (sha256 {log['panel_c_sha256'][:16]}) at (0, {Y0}), letter c 14 pt\n" + "\n".join(lines) + "\n"); print(lines[-1])
