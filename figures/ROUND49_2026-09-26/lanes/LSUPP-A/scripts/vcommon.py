#!/usr/bin/env python3
"""Lane V7_L1_RANK, shared proof helpers (03_verify): letters, font faces, token multisets by region, read-back of placed strings
by position, raster diff outside a region, side-by-side crops."""
import os, re, json, gc
from collections import Counter
import fitz
import common as C

FACE = {"ArialMT": "Arial", "Arial Regular": "Arial", "Arial-BoldMT": "Arial Bold", "Arial Bold": "Arial Bold", "Arial-ItalicMT": "Arial Italic"}
NUM_RE = re.compile(r"[+−\-]?\d[\d,]*\.?\d*")


def face(fontname):
    return FACE.get(fontname.split("+")[-1], fontname)


def sheet_spans(pdf_path):
    doc = fitz.open(C.hydrated(pdf_path)); pg = doc[0]
    sp = C.spans(pg); rect = [round(pg.rect.width, 2), round(pg.rect.height, 2)]
    lets = C.letters(pg); doc.close(); gc.collect()
    return sp, rect, lets


def in_region(s, region):
    """region = (y0, y1) or (x0, y0, x1, y1) on the span's origin."""
    ox, oy = s["origin"]
    if len(region) == 2: return region[0] <= oy < region[1]
    return region[0] <= ox < region[2] and region[1] <= oy < region[3]


def tokens(spans, region=None, drop_letters=True):
    """Word and number multisets of the spans (words are whitespace tokens with no digit, numbers are the numeric tokens)."""
    words, nums = Counter(), Counter()
    for s in spans:
        if drop_letters and len(s["text"].strip()) == 1 and s["text"].strip() in "abcdefgh" and s["size"] >= 12.5: continue
        if region is not None and not in_region(s, region): continue
        for w in s["text"].replace("\xa0", " ").split():
            ns = NUM_RE.findall(w)
            if ns:
                for n in ns: nums[n.replace("−", "-")] += 1
                rest = NUM_RE.sub("", w).strip("(),%")
                if rest: words[rest] += 1
            else:
                words[w] += 1
    return words, nums


def diff_counter(a, b):
    """(lost from a, gained in b) as sorted lists with multiplicity."""
    lost = sorted(((a - b)).elements()); gained = sorted(((b - a)).elements())
    return lost, gained


def compare_letters(base_lets, new_lets, tol=0.3, shifts=None, plus=1.0):
    """Same letters in the same reading order, each within tol pt (a declared per-letter y shift may be passed). Round 49: the new letters are
    one point larger (13 -> 14), so the size check expects base + plus and the box is compared on its left x and its baseline-side top
    (the bbox top rises with the larger ascent, so the top is compared after removing the ascent growth)."""
    shifts = shifts or {}
    if [l[0] for l in base_lets] != [l[0] for l in new_lets]: return False, f"letters {[l[0] for l in base_lets]} vs {[l[0] for l in new_lets]}"
    worst = 0.0
    for b, n in zip(base_lets, new_lets):
        dy = shifts.get(b[0], 0.0)
        asc_growth = (n[5] - b[5]) * 0.905; dsc_growth = (n[5] - b[5]) * 0.212      # Arial ascent 0.905, descent 0.212 per point
        d = max(abs(b[1] - n[1]), abs(b[2] + dy - asc_growth - n[2]), abs(b[3] + dy + dsc_growth - n[3]))
        worst = max(worst, d)
        if face(b[4]) != face(n[4]) or abs(b[5] + plus - n[5]) > 0.05: return False, f"letter {b[0]} font {b[4]} {b[5]} vs {n[4]} {n[5]} (expected base + {plus:g})"
    return worst <= tol, f"worst letter offset {worst:.3f} pt"


def find_span(spans, text, x, baseline, ha="left", tol=0.6, size=None):
    """The span with this exact text whose anchor (origin x, bbox right or bbox centre by ha) and baseline sit within tol."""
    best = None
    for s in spans:
        if s["text"] != text: continue
        if size is not None and abs(s["size"] - size) > 0.05: continue
        ax = s["origin"][0] if ha == "left" else (s["bbox"][2] if ha == "right" else (s["bbox"][0] + s["bbox"][2]) / 2)
        d = max(abs(ax - x), abs(s["origin"][1] - baseline))
        if best is None or d < best[0]: best = (d, s)
    if best is None: return None, None
    return (best[1] if best[0] <= tol else None), best[0]


def raster_outside(base_pdf, new_pdf, keep_region, dpi=150, lock=False, pad_px=2, extra_boxes=()):
    """Render both at dpi, return the count of differing pixels OUTSIDE keep_region (page pt, x0 y0 x1 y1) and inside, and the
    bbox of the outside differences. Both PNGs are kept beside the new sheet for the eye."""
    import numpy as np
    from PIL import Image
    wd = os.path.dirname(new_pdf); nm = os.path.splitext(os.path.basename(new_pdf))[0]
    a_png = f"{wd}/work/_base_{dpi}.png"; b_png = f"{wd}/work/_new_{dpi}.png"
    import shutil as _sh
    _sh.copyfile(C.full_render_cached(base_pdf, dpi, lock=True), a_png); _sh.copyfile(C.full_render_cached(new_pdf, dpi, lock=True), b_png)   # round 37: one cached child-process render per sheet and dpi
    a = np.asarray(Image.open(a_png).convert("RGB")).astype(int); b = np.asarray(Image.open(b_png).convert("RGB")).astype(int)
    shapes = (list(a.shape), list(b.shape))
    if a.shape != b.shape:
        # a page box equal within 0.01 pt can still render one pixel taller or wider (ceil of the float box): align on the common size, refuse more than that
        if abs(a.shape[0] - b.shape[0]) > 1 or abs(a.shape[1] - b.shape[1]) > 1: return dict(shape_a=a.shape, shape_b=b.shape, outside=-1, inside=-1)
        hh, ww = min(a.shape[0], b.shape[0]), min(a.shape[1], b.shape[1]); a = a[:hh, :ww]; b = b[:hh, :ww]
    d = np.abs(a - b).sum(axis=2) > 30
    k = dpi / 72.0
    x0, y0, x1, y1 = keep_region
    m = np.zeros_like(d); m[max(0, int(y0 * k) - pad_px):int(np.ceil(y1 * k)) + pad_px, max(0, int(x0 * k) - pad_px):int(np.ceil(x1 * k)) + pad_px] = True
    for (ex0, ey0, ex1, ey1) in extra_boxes:          # the re-stamped letter boxes (the proof allows them: rebuilt region plus the letter boxes)
        m[max(0, int(ey0 * k) - pad_px):int(np.ceil(ey1 * k)) + pad_px, max(0, int(ex0 * k) - pad_px):int(np.ceil(ex1 * k)) + pad_px] = True
    outside = d & ~m; inside = d & m
    ys, xs = np.where(outside)
    bbox = None if len(ys) == 0 else [round(xs.min() / k, 1), round(ys.min() / k, 1), round(xs.max() / k, 1), round(ys.max() / k, 1)]
    return dict(outside=int(outside.sum()), inside=int(inside.sum()), outside_bbox_pt=bbox, shape=list(a.shape), rendered_shapes=shapes)


def geometry_records(pdf_path, region=None):
    """Deduplicated vector records (fill, stroke, lw, kinds, bbox) of a page, optionally inside a region (x0 y0 x1 y1)."""
    doc = fitz.open(C.hydrated(pdf_path)); pg = doc[0]; seen = set(); out = []
    for d in pg.get_drawings():
        r = d["rect"]; b = (round(r.x0, 2), round(r.y0, 2), round(r.x1, 2), round(r.y1, 2))
        if region is not None and (b[2] < region[0] or b[0] > region[2] or b[3] < region[1] or b[1] > region[3]): continue
        kinds = "".join(it[0] for it in d.get("items", []))
        key = (b, C.__dict__.get("col", None) and None, d.get("fill") and tuple(round(v, 3) for v in d["fill"][:3]), d.get("color") and tuple(round(v, 3) for v in d["color"][:3]), round(d.get("width") or 0, 3), kinds[:6])
        if key in seen: continue
        seen.add(key)
        out.append(dict(bbox=list(b), fill=None if d.get("fill") is None else "#%02x%02x%02x" % tuple(int(round(v * 255)) for v in d["fill"][:3]),
                        stroke=None if d.get("color") is None else "#%02x%02x%02x" % tuple(int(round(v * 255)) for v in d["color"][:3]),
                        lw=round(d.get("width") or 0, 3), kinds=kinds[:8], dashes=d.get("dashes")))
    doc.close(); gc.collect()
    return out


def has_segment(recs, x0, y0, x1, y1, stroke=None, tol=0.3, lw=None):
    """Is there a stroked record whose bbox matches this segment within tol (any lw unless given)?"""
    for r in recs:
        if stroke is not None and r["stroke"] != stroke: continue
        if lw is not None and abs(r["lw"] - lw) > 0.05: continue
        b = r["bbox"]
        if max(abs(b[0] - x0), abs(b[1] - y0), abs(b[2] - x1), abs(b[3] - y1)) <= tol: return True
    return False


# ----------------------------------------------------------------------------------------------- proof W helpers (relaunch 5)
def pseudo(strings, size=10.0):
    """Pseudo spans so a list of strings can go through the same tokenizer as a text layer."""
    return [dict(text=s, size=size, origin=[0.0, 0.0], font="ArialMT") for s in strings]


def residue_check(check, name, base_spans, new_spans, old_strings, new_strings, region=None):
    """Proof W multiset check: the base tokens minus the strings the OLD numbers file prints must equal the new tokens minus the strings
    the NEW file prints (words and numbers separately), and every old string must be on the base sheet."""
    bw, bn = tokens(base_spans, region=region); nw, nn = tokens(new_spans, region=region)
    ow, on = tokens(pseudo(old_strings)); xw, xn = tokens(pseudo(new_strings))
    miss_old = sorted((ow - bw).elements()) + sorted((on - bn).elements())
    rw, rn = bw - ow, bn - on; sw, sn = nw - xw, nn - xn
    ok = (not miss_old) and rw == sw and rn == sn
    check(ok, name, f"old strings missing on the base {miss_old[:10]}; static residue words lost/gained {diff_counter(rw, sw)}; numbers lost/gained {diff_counter(rn, sn)}")
    return dict(ok=ok, words_lost_gained=diff_counter(ow, xw), numbers_lost_gained=diff_counter(on, xn))


def has_rect(recs, x0, y0, x1, y1, fill=None, tol=0.3):
    for r in recs:
        if fill is not None and (r["fill"] or "").lower() != fill.lower(): continue
        b = r["bbox"]
        if max(abs(b[0] - x0), abs(b[1] - y0), abs(b[2] - x1), abs(b[3] - y1)) <= tol: return True
    return False


def has_dot(recs, cx, cy, fill, diam=5.29, tol=0.35):
    for r in recs:
        if (r["fill"] or "").lower() != fill.lower(): continue
        b = r["bbox"]; w = b[2] - b[0]; h = b[3] - b[1]
        if abs(w - diam) <= 0.7 and abs(h - diam) <= 0.7 and abs((b[0] + b[2]) / 2 - cx) <= tol and abs((b[1] + b[3]) / 2 - cy) <= tol: return True
    return False


def readback(check, name, spans, values, tol=0.6, kinds=("printed", "tick", "label")):
    fails = []; n = 0
    for v in values:
        if v["kind"] not in kinds: continue
        n += 1; s, d = find_span(spans, v["text"], v["x"], v["baseline"], ha=v.get("ha", "left"), tol=tol)
        if s is None: fails.append((v["text"], v["x"], v["baseline"], None if d is None else round(d, 2)))
    check(not fails, f"{name}: {n} placed strings found on the text layer at their positions (within {tol} pt)", f"fails {fails[:8]}")
    return fails


def fonts_check(check, bsp, nsp, plus=1.0, kept=()):
    """Round 49: the new sheet carries the base's font faces at every size plus one point (10 -> 11, 9 -> 10, 13 -> 14), plus the declared
    (face, size) pairs of elements kept at their current size because a placement rule fails at +1 (listed in the lane REPORT)."""
    bf = {(face(s["font"]), round(s["size"] + plus, 2)) for s in bsp} | set(kept); nf = {(face(s["font"]), s["size"]) for s in nsp}
    check(bf == nf, f"font faces identical and every size the base size plus {plus:g} pt" + (f" (declared kept: {sorted(kept)})" if kept else ""), f"base+{plus:g}-only {sorted(bf - nf)} new-only {sorted(nf - bf)}")
