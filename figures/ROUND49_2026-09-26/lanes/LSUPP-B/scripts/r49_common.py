"""Lane LSUPP-B (round 49, 2026-09-26): shared helpers for Supp_Fig05, 06, 07 and 18 at +1 pt.

Rules of the round folded in: the text layer of any sheet is read with Ghostscript's txtwrite device (never PyMuPDF get_text on a
composed sheet), every render is Ghostscript png16m run one at a time, page boxes come from pypdf, PyMuPDF is used for placement
only (show_pdf_page, TextWriter) and for get_drawings on FLAT matplotlib parts (never on a composed sheet)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import hashlib, html, json, os, re, subprocess, time
from collections import Counter

GS = paths.GS
T90 = paths.FIGURE_ROOT
R49 = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26"
LANE = f"{R49}/lanes/LSUPP-B"
V26 = f"{paths.FIGURE_ROOT}/ROUND47_2026-09-24/figures/NEW_FINAL_SET_V26"
PT_PLUS = 1.0          # round 49 (Alen, 2026-09-26): every text one point larger
NUMRE = re.compile(r"-?\d[\d,]*\.?\d*")


# ---------------------------------------------------------------- files
def blocks(path):
    out = subprocess.run(["stat", "-f", "%b %z", path], capture_output=True, text=True, check=True).stdout.split()
    return int(out[0]), int(out[1])


def hydrated(path, wait_s=600):
    """iCloud eviction gate: nonzero size with zero allocated blocks = evicted. brctl download and poll, never read it evicted."""
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    b, z = blocks(path)
    if z > 0 and b == 0:
        subprocess.run(["brctl", "download", path], capture_output=True)
        t0 = time.time()
        while blocks(path)[0] == 0:
            if time.time() - t0 > wait_s:
                raise RuntimeError(f"EVICTED and not downloaded in {wait_s} s: {path}")
            time.sleep(2)
    return path


def sha256(path, n=None):
    h = hashlib.sha256()
    with open(hydrated(path), "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest() if n is None else h.hexdigest()[:n]


def now():
    return time.strftime("%Y-%m-%d %H:%M:%S")


def page_box(pdf):
    """(n_pages, width, height) from pypdf's mediabox (no content interpretation)."""
    from pypdf import PdfReader
    r = PdfReader(hydrated(pdf))
    mb = r.pages[0].mediabox
    return len(r.pages), float(mb.width), float(mb.height)


# ---------------------------------------------------------------- Ghostscript: renders and the text layer
def gs_render(pdf, png, dpi=150, timeout=900):
    """png16m render of the one page, anti-aliased text and graphics, one process at a time (the caller serialises)."""
    hydrated(pdf)
    os.makedirs(os.path.dirname(png), exist_ok=True)
    r = subprocess.run([GS, "-dSAFER", "-dBATCH", "-dNOPAUSE", "-dQUIET", "-sDEVICE=png16m", f"-r{dpi}", "-dTextAlphaBits=4",
                        "-dGraphicsAlphaBits=4", f"-sOutputFile={png}", pdf], capture_output=True, text=True, timeout=timeout)
    assert r.returncode == 0 and os.path.exists(png) and os.path.getsize(png) > 0, (pdf, r.stderr[-500:])
    return png


_SPAN = re.compile(r'<span bbox="(\S+) (\S+) (\S+) (\S+)" font="([^"]*)" size="([^"]*)">(.*?)</span>', re.S)
_CHAR = re.compile(r'<char bbox="(\S+) (\S+) (\S+) (\S+)" c="(.*?)"/>')


def gs_text(pdf, out_xml, timeout=600):
    """Spans of the sheet from Ghostscript's txtwrite device (TextFormat 0): text, x, x1, y (baseline, top-down page points, integer
    precision), font, size. The device splits one PDF text run at kerning pairs, so merge_kerned() is applied before any comparison."""
    hydrated(pdf)
    os.makedirs(os.path.dirname(out_xml), exist_ok=True)
    r = subprocess.run([GS, "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=txtwrite", "-dTextFormat=0", f"-sOutputFile={out_xml}", pdf],
                       capture_output=True, text=True, timeout=timeout)
    assert r.returncode == 0, r.stderr[-500:]
    return parse_txtwrite(out_xml)


def font_base(name):
    return name.split("+")[-1].replace("-Identity-H", "")


def parse_txtwrite(path):
    raw = open(path, encoding="utf-8", errors="replace").read()
    out = []
    for m in _SPAN.finditer(raw):
        x0, y0, x1, y1 = [float(v) for v in m.group(1, 2, 3, 4)]
        chars = [(float(c[0]), float(c[2]), html.unescape(c[4])) for c in _CHAR.findall(m.group(7))]
        text = "".join(c[2] for c in chars).replace("\xa0", " ")
        if not text.strip():
            continue
        lead = len(text) - len(text.lstrip())
        trail = len(text) - len(text.rstrip())
        xs = chars[lead][0] if lead < len(chars) else x0
        xe = chars[len(chars) - 1 - trail][1] if trail < len(chars) else x1
        out.append(dict(text=text.strip(), raw=text, x=xs, x1=xe, y=y0, y1=y1, font=font_base(m.group(5)), size=float(m.group(6))))
    out.sort(key=lambda s: (s["y"], s["x"]))
    return out


def merge_kerned(spans, gap=1.5, dy=0.5):
    """txtwrite splits a text run at kerning pairs ('T' + 'ype 2 diabetes', '1.1' + '1 (0.84-1.47)'): glue spans on one baseline with
    the same face and size whose gap is at most `gap` pt (a trailing space in the first piece is kept: 'adjusted for ' + 'T90')."""
    out = []
    for s in sorted(spans, key=lambda s: (s["y"], s["x"])):
        if out:
            p = out[-1]
            if abs(p["y"] - s["y"]) <= dy and p["font"] == s["font"] and p["size"] == s["size"] and -4.0 <= s["x"] - p["x1"] <= gap:
                p["raw"] = p["raw"].rstrip("\n") + s["raw"]
                p["text"] = p["raw"].strip()
                p["x1"] = max(p["x1"], s["x1"])
                p["merged"] = p.get("merged", 1) + 1
                continue
        out.append(dict(s))
    return out


def text_box(s, pad=2.0):
    """Page-point box of a census span from its baseline and size (top-down): ascent 1.05 em, descent 0.35 em, `pad` pt around.
    A span rotated 90 degrees (txtwrite reports x0 == x1 and y0 > y1, the baseline running up the page) gets the box the other way:
    the glyphs extend left of the baseline x, from y1 (top) to y0 (bottom)."""
    y1 = s.get("y1", s["y"])
    if abs(s["x1"] - s["x"]) < 2.0 and abs(y1 - s["y"]) > 2.0:
        return [s["x"] - 1.05 * s["size"] - pad, min(s["y"], y1) - pad, s["x"] + 0.35 * s["size"] + pad, max(s["y"], y1) + pad]
    return [s["x"] - pad, s["y"] - 1.05 * s["size"] - pad, s["x1"] + pad, s["y"] + 0.35 * s["size"] + pad]


# ---------------------------------------------------------------- the census comparison (strings equal, sizes +PT_PLUS, exceptions declared)
def census_compare(old, new, exceptions=None, plus=PT_PLUS):
    """old, new: merged span lists. Every old (text, face, size) must appear in new as (text, face, size + plus), except the declared
    exceptions {(text, face, old_size): new_size} (a string kept at a smaller step). Returns (ok, detail dict)."""
    exceptions = exceptions or {}
    want = Counter()
    for s in old:
        k = (s["text"], s["font"], s["size"])
        want[(s["text"], s["font"], exceptions.get(k, s["size"] + plus))] += 1
    have = Counter((s["text"], s["font"], s["size"]) for s in new)
    missing = want - have
    extra = have - want
    strings_equal = Counter(s["text"] for s in old) == Counter(s["text"] for s in new)
    return (not missing and not extra), dict(strings_equal=strings_equal, missing=dict(missing), extra=dict(extra),
                                             n_old=len(old), n_new=len(new), exceptions_used={f"{k[0]!r} {k[1]} {k[2]} -> {v}": 0 for k, v in exceptions.items()})


def size_table(spans):
    return sorted(Counter((s["font"], s["size"]) for s in spans).items())


# ---------------------------------------------------------------- raster diff outside declared boxes (Ghostscript renders, PIL only)
def raster_diff(png_a, png_b, allowed_boxes_pt, dpi, thresh=32, dilate_pt=1.5):
    import numpy as np
    from PIL import Image
    a = np.asarray(Image.open(png_a).convert("RGB")).astype(int)
    b = np.asarray(Image.open(png_b).convert("RGB")).astype(int)
    if a.shape != b.shape:
        return dict(n_diff=None, n_outside=None, shape_equal=False, shapes=[a.shape, b.shape])
    d = (np.abs(a - b) > thresh).any(axis=2)
    mask = np.zeros_like(d)
    z = dpi / 72.0
    for (x0, y0, x1, y1) in allowed_boxes_pt:
        X0, Y0 = max(int((x0 - dilate_pt) * z), 0), max(int((y0 - dilate_pt) * z), 0)
        X1, Y1 = min(int((x1 + dilate_pt) * z) + 1, d.shape[1]), min(int((y1 + dilate_pt) * z) + 1, d.shape[0])
        mask[Y0:Y1, X0:X1] = True
    ys, xs = np.nonzero(d & ~mask)
    bbox = None if len(ys) == 0 else [round(xs.min() / z, 1), round(ys.min() / z, 1), round(xs.max() / z, 1), round(ys.max() / z, 1)]
    return dict(n_diff=int(d.sum()), n_outside=int((d & ~mask).sum()), shape_equal=True, outside_bbox_pt=bbox, n_masked_px=int(mask.sum()))


def crop_png(png, dpi, box_pt, out_png, scale=1):
    from PIL import Image
    im = Image.open(png).convert("RGB")
    z = dpi / 72.0
    x0, y0, x1, y1 = box_pt
    c = im.crop((max(0, int(x0 * z)), max(0, int(y0 * z)), min(im.width, int(x1 * z)), min(im.height, int(y1 * z))))
    if scale != 1:
        c = c.resize((c.width * scale, c.height * scale), Image.LANCZOS)
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    c.save(out_png)
    return c.size


# ---------------------------------------------------------------- the vector layer of a FLAT matplotlib part (never a composed sheet)
def flat_drawings(pdf_path, nd=2):
    """get_drawings() on a flat matplotlib PDF (one page, no form XObjects): one record per path with its rect, fill, stroke, width,
    dashes and item kinds. Refuses a page that carries XObjects (a composed sheet)."""
    import fitz
    doc = fitz.open(hydrated(pdf_path))
    page = doc[0]
    assert len(doc) == 1 and len(page.get_xobjects()) == 0, f"not a flat page: {pdf_path}"
    items = []
    for g in page.get_drawings():
        f, c, r = g.get("fill"), g.get("color"), g["rect"]
        fh = "#%02x%02x%02x" % tuple(int(round(v * 255)) for v in f) if f else None
        ch = "#%02x%02x%02x" % tuple(int(round(v * 255)) for v in c) if c else None
        kinds = "".join(sorted(set(it[0] for it in g["items"])))
        items.append(dict(rect=[round(v, nd) for v in (r.x0, r.y0, r.x1, r.y1)], fill=fh, stroke=ch,
                          width=None if g.get("width") is None else round(g.get("width"), 3), kinds=kinds, n=len(g["items"]),
                          dashes=g.get("dashes")))
    doc.close()
    return items


def drawings_key(d):
    return (tuple(d["rect"]), d["fill"], d["stroke"], d["width"], d["kinds"], d["n"], str(d["dashes"]))


def drawings_compare(old_items, new_items, allowed_boxes=()):
    """Multiset equality of the drawing records; records whose rect lies inside an allowed box (legend handles) are excluded."""
    def inside(d):
        x0, y0, x1, y1 = d["rect"]
        return any(bx0 - 0.5 <= x0 and x1 <= bx1 + 0.5 and by0 - 0.5 <= y0 and y1 <= by1 + 0.5 for bx0, by0, bx1, by1 in allowed_boxes)
    co = Counter(drawings_key(d) for d in old_items if not inside(d))
    cn = Counter(drawings_key(d) for d in new_items if not inside(d))
    return (co == cn), dict(n_old=sum(co.values()), n_new=sum(cn.values()), only_old=list((co - cn).keys())[:12], only_new=list((cn - co).keys())[:12],
                            excluded_old=sum(1 for d in old_items if inside(d)), excluded_new=sum(1 for d in new_items if inside(d)))


# ---------------------------------------------------------------- checks
class Checks:
    def __init__(self, header):
        self.lines = [header, ""]; self.n_pass = 0; self.n_fail = 0

    def log(self, what, ok, detail=""):
        self.lines.append(f"{'PASS' if ok else 'FAIL'}: {what}" + (f" [{detail}]" if detail != "" else ""))
        if ok: self.n_pass += 1
        else: self.n_fail += 1
        print(self.lines[-1][:600], flush=True)
        return bool(ok)

    def info(self, what):
        self.lines.append(f"INFO: {what}"); print(self.lines[-1][:600], flush=True)

    def write(self, path):
        ok = self.n_fail == 0
        self.lines += ["", f"{self.n_pass} pass, {self.n_fail} fail", "RESULT ALL PASS" if ok else "RESULT FAIL"]
        os.makedirs(os.path.dirname(path), exist_ok=True)
        open(path, "w").write("\n".join(self.lines) + "\n")
        return ok
