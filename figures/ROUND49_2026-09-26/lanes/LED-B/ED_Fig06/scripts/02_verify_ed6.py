#!$T90_PY
"""ED_Fig06 round 49 verify (child under wd_run.sh): the re-set sheet against the V26 dump and the declared placements.
usage: 02_verify_ed6.py <new.pdf> <v26_dump.json> <reset_record.json> <work dir> <verify dir> <png_v26> <png_notext> <png_new>
Checks: page box; 116 spans read back, each at the declared anchor (0.3 pt), baseline (0.05), size +1, weight, colour; no extra span; fonts Arial
faces only, no XObject; the 284 V26 drawing items on the new sheet unchanged (rect 0.05, fill, stroke, width, dashes) and nothing else drawn;
raster (Ghostscript 150 dpi): V26 -> NOTEXT differs only inside the visible V26 text boxes (the clipped stale copies were invisible),
V26 -> NEW differs only inside the union of the V26 and new text boxes (marks, cells, lines untouched); values against the data files
(habitual_outcomes_v2 cells_crossed.csv and results_crossed.csv under the round-37 sidecar gate): 12 counts, 24 events/n, 24 HR (95% CI)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import fitz, json, os, sys, re, gc
import numpy as np, pandas as pd
from PIL import Image
NEW, DUMP, REC, WORK, VER, PNG_OLD, PNG_NOTEXT, PNG_NEW = sys.argv[1:9]
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND37_2026-09-14/figures/lanes/V14_L6_HABITUAL_EXTERNAL/_common")
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND37_2026-09-14/figures/_common")
import lane_common as LC, v14lib as V
lines = []
def check(name, ok, detail=""): lines.append(f"{'PASS' if ok else 'FAIL'}  {name}" + (f": {detail}" if detail else "")); print(lines[-1][:400])
dump = json.load(open(DUMP)); rec = json.load(open(REC)); placed = rec["placed"]
d = fitz.open(NEW); p = d[0]; W, H = p.rect.width, p.rect.height
check("page box equal to V26", abs(W - dump["page"][0]) < 0.01 and abs(H - dump["page"][1]) < 0.01, f"{W:.3f} x {H:.3f}")
new = []
for b in p.get_text("rawdict")["blocks"]:
    if b["type"] != 0: continue
    for l in b["lines"]:
        for s in l["spans"]:
            txt = "".join(ch["c"] for ch in s["chars"])
            if not txt.strip(): continue
            new.append(dict(text=txt, font=s["font"], size=round(s["size"], 3), color="#%06x" % s["color"], origin=[round(v, 3) for v in s["chars"][0]["origin"]], bbox=[round(v, 3) for v in s["bbox"]], dir=[round(v, 3) for v in l["dir"]]))
check("text layer clean: no NBSP and no soft hyphen in any span", all("\xa0" not in s["text"] and "\xad" not in s["text"] for s in new), f"{len(new)} spans")
used = [False] * len(new); miss = []; worst = 0.0
for q in placed:
    hit = None
    for j, s in enumerate(new):
        if used[j] or s["text"] != q["text"] or abs(s["size"] - q["size"]) > 0.01 or (("Bold" in s["font"]) != q["bold"]) or s["color"] != q["color"]: continue
        if q["rotated"]:
            if s["dir"] == [0.0, -1.0] and abs(s["origin"][0] - q["x"]) < 0.3 and abs(s["origin"][1] - q["baseline"]) < 0.3: hit = j; break
        elif abs(s["origin"][1] - q["baseline"]) < 0.05:
            w = s["bbox"][2] - s["bbox"][0]
            got = s["origin"][0] if q["kind"] == "left" else (s["bbox"][2] if q["kind"] == "right" else s["bbox"][0] + w / 2)
            if abs(got - q["anchor"]) < 0.3: worst = max(worst, abs(got - q["anchor"])); hit = j; break
    if hit is None: miss.append((q["text"], q["kind"], q["anchor"], q["baseline"]))
    else: used[hit] = True
check(f"every declared placement ({len(placed)}) read back: text, size +1, weight, colour, baseline (0.05 pt) and anchor (0.3 pt)", not miss, f"largest anchor offset {worst:.3f} pt; missing {miss[:6]}")
extra = [(s["text"], s["origin"]) for j, s in enumerate(new) if not used[j]]
check("no span beyond the declared placements", not extra, f"{extra[:6]}")
fonts = sorted({f[3].split("+")[-1] for f in p.get_fonts(full=True)})
check("fonts: Arial faces only (system Arial Regular and Arial Bold), no XObject, no image", all(f.replace(" ", "").lower().startswith("arial") for f in fonts) and len(p.get_xobjects()) == 0 and len(p.get_images()) == 0, f"{fonts}")
# sizes: every V26 visible span's size + 1 exactly (by text and count)
from collections import Counter
co = Counter((s["text"].replace("\xa0", " ").replace("\xad", "-"), round(s["size"] + 1.0, 3), "Bold" in s["font"]) for s in dump["spans"]); cn = Counter((s["text"], s["size"], "Bold" in s["font"]) for s in new)
check("every V26 visible string present once at exactly +1.0 pt with its weight (116 of 116), nothing else", co == cn, f"lost {dict(co - cn)} gained {dict(cn - co)}")
# drawings identity
def key(it): return (tuple(round(v, 2) for v in it["rect"]), tuple(round(v, 4) for v in it["fill"]) if it["fill"] else None, tuple(round(v, 4) for v in it["color"]) if it["color"] else None, round(it["width"] or 0, 3), str(it.get("dashes")), it["kinds"], it["n"])
nd = []
for it in p.get_drawings():
    r = it["rect"]; kinds = "".join(sorted(set(k[0] for k in it["items"])))
    nd.append(dict(rect=[r.x0, r.y0, r.x1, r.y1], fill=it.get("fill"), color=it.get("color"), width=it.get("width"), dashes=it.get("dashes"), kinds=kinds, n=len(it["items"])))
od = dump["drawings"]
def style(it): return (tuple(round(v, 4) for v in it["fill"]) if it["fill"] else None, tuple(round(v, 4) for v in it["color"]) if it["color"] else None, round(it["width"] or 0, 3), str(it.get("dashes")), it["kinds"], it["n"])
taken = [False] * len(nd); unmatched = []; worst_r = 0.0
for it in od:
    hit = None
    for j, jt in enumerate(nd):
        if taken[j] or style(jt) != style(it): continue
        dr = max(abs(a - b) for a, b in zip(it["rect"], jt["rect"]))
        if dr <= 0.05: hit = j; worst_r = max(worst_r, dr); break
    if hit is None: unmatched.append(it["rect"])
    else: taken[hit] = True
check(f"line art identical to V26: each of the {len(od)} V26 drawing items found once on the new sheet with the same fill, stroke, width, dashes and path kinds and its rect within 0.05 pt (the redaction pass re-serialises coordinates to 0.01 pt), nothing added", not unmatched and all(taken), f"unmatched {len(unmatched)}, extra {sum(1 for t in taken if not t)}, largest rect offset {worst_r:.3f} pt")
d.close(); gc.collect()
# raster proofs (Ghostscript renders made by the runner)
def arr(png): return np.asarray(Image.open(png).convert("RGB")).astype(int)
A, B, C = arr(PNG_OLD), arr(PNG_NOTEXT), arr(PNG_NEW); S = 150 / 72
check("150 dpi renders of V26, NOTEXT and NEW have the same size", A.shape == B.shape == C.shape, f"{A.shape}")
def boxes_old():
    out = []
    for s in dump["spans"]:
        x0, y0, x1, y1 = s["bbox"]; out.append((x0 - 0.5, y0 - 0.5, x1 + 0.5, y1 + 0.5))
    return out
def boxes_new():
    out = []
    for q in placed:
        if q["rotated"]:
            x, y = q["x"], q["baseline"]; out.append((x - 0.95 * q["size"] - 0.5, y - q["w_new"] - 0.5, x + 0.3 * q["size"] + 0.5, y + 0.5))
        else:
            x, y = q["x"], q["baseline"]; out.append((x - 0.5, y - 0.95 * q["size"] - 0.5, x + q["w_new"] + 0.5, y + 0.3 * q["size"] + 0.5))
    return out
def mask(boxes, shape):
    m = np.zeros(shape[:2], bool)
    for (x0, y0, x1, y1) in boxes: m[max(0, int(y0 * S)):int(y1 * S) + 1, max(0, int(x0 * S)):int(x1 * S) + 1] = True
    return m
diff_notext = np.abs(A - B).max(axis=2) > 0; m_old = mask(boxes_old(), A.shape)
check("V26 -> NOTEXT: pixels differ only inside the visible V26 text boxes (the clipped stale text copies had no pixels: removing them changes nothing)", int((diff_notext & ~m_old).sum()) == 0, f"{int(diff_notext.sum())} differing pixels, {int((diff_notext & ~m_old).sum())} outside the text boxes")
diff_new = np.abs(A - C).max(axis=2) > 0; m_all = m_old | mask(boxes_new(), A.shape)
check("V26 -> NEW: pixels differ only inside the union of the V26 and the new text boxes (cells, marks, lines, key handles untouched)", int((diff_new & ~m_all).sum()) == 0, f"{int(diff_new.sum())} differing pixels, {int((diff_new & ~m_all).sum())} outside")
# values against the data files (round-37 gate)
SRC = f"{LC.SV}/habitual_outcomes_v2"; MIN_MTIME = "2026-09-15 19:18:00"
sc_c = LC.sidecar(LC.hydrated(f"{SRC}/cells_crossed.csv"), min_mtime=MIN_MTIME); sc_r = LC.sidecar(LC.hydrated(f"{SRC}/results_crossed.csv"), min_mtime=MIN_MTIME)
cells = pd.read_csv(f"{SRC}/cells_crossed.csv"); crossed = pd.read_csv(f"{SRC}/results_crossed.csv")
A_ = cells[cells.design == "cells12"].set_index("cell"); counts = sorted(f"n = {int(n):,}" for n in A_.n)
HAB = ["<5h", "5-<6h", "6-<7h(ref)", ">=7h"]; O2 = ["normal T90<=1%", "intermediate 1-10%", "low T90>10%"]
OUTC = [("cvd", "Cardiovascular composite"), ("death", "Death from any cause"), ("htn2", "Hypertension"), ("diabetes", "Type 2 diabetes"), ("hf", "Heart failure"), ("dementia", "Dementia")]
CELLS4 = [("<5h", "normal T90<=1%"), ("<5h", "low T90>10%"), (">=7h", "normal T90<=1%"), (">=7h", "low T90>10%")]
cb = crossed[(crossed.model == "cells12") & (crossed["set"] == "primary")]; exp_ev = []; exp_hr = []
for k, name in OUTC:
    dd = cb[cb.outcome_key == k].set_index("contrast")
    for hb, ob in CELLS4:
        r = dd.loc[f"{hb} | {ob}"]; exp_ev.append(f"{int(r.events_group)}/{int(r.n_group):,}"); exp_hr.append(V.hr_ci(float(r.hr), float(r.lo95), float(r.hi95)))
got_counts = sorted(s["text"] for s in new if s["text"].startswith("n = ")); got_ev = [s["text"] for s in sorted(new, key=lambda s: s["origin"][1]) if re.fullmatch(r"\d+/[\d,]+", s["text"])]
got_hr = [s["text"] for s in sorted(new, key=lambda s: s["origin"][1]) if re.fullmatch(r"\d\.\d\d \(\d\.\d\d-\d\.\d\d\)", s["text"])]
check("printed values re-derived from the data files (sidecar-gated): 12 group counts, 24 events/group size in row order, 24 HR (95% CI) in row order", got_counts == counts and got_ev == exp_ev and got_hr == exp_hr, f"counts {got_counts == counts}, events {got_ev == exp_ev}, HR {got_hr == exp_hr}; sidecars {sc_c['_gate']['output_mtime']} / {sc_r['_gate']['output_mtime']}")
alltext = " ".join(s["text"] for s in new)
check("house rules: no em dash, no semicolon on the sheet", "—" not in alltext and ";" not in alltext, "")
json.dump(dict(checks=lines, n_spans=len(new), fonts=fonts, spans=new), open(f"{VER}/verify_record.json", "w"), indent=1, ensure_ascii=False)
open(f"{VER}/verify_lines.txt", "w").write("\n".join(lines) + "\n")
print("VERIFY", "ALL PASS" if all(l.startswith("PASS") for l in lines) else "FAIL"); sys.exit(0 if all(l.startswith("PASS") for l in lines) else 1)
