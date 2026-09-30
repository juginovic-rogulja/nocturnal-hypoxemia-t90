"""Lane LSUPP-C (round 49, 2026-09-26): shared helpers. Ghostscript only for text layers (txtwrite, TextFormat 0) and renders
(png16m), one job at a time under the round-38 watchdog wd_run.sh (WD_LOGDIR = this lane's logs/watchdog); PyMuPDF for
placement only, in a child process (fitz_place.py). Every file is eviction-checked before it is read."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import hashlib, html, json, os, re, subprocess, time
from collections import Counter

T90 = paths.FIGURE_ROOT
LANE = f"{paths.FIGURE_ROOT}/ROUND50_2026-09-26/lanes/LSUPP-C"   # R50
V26 = f"{paths.FIGURE_ROOT}/ROUND47_2026-09-24/figures/NEW_FINAL_SET_V26"
V27 = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/figures/NEW_FINAL_SET_V27"   # R50: the base set (every text one point larger)
WD = f"{paths.FIGURE_ROOT}/ROUND38_2026-09-16/figures/scripts/wd_run.sh"
GS = paths.GS
PY = paths.PY
WDLOG = f"{LANE}/logs/watchdog"
os.makedirs(WDLOG, exist_ok=True)
PT_PLUS = 1.0
NUMRE = re.compile(r"-?\d[\d,]*\.?\d*")


def blocks(path):
    out = subprocess.run(["stat", "-f", "%b %z", path], capture_output=True, text=True, check=True).stdout.split()
    return int(out[0]), int(out[1])


def hydrated(path, wait_s=600):
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    b, z = blocks(path)
    if z > 0 and b == 0:
        subprocess.run(["brctl", "download", path], capture_output=True); t0 = time.time()
        while blocks(path)[0] == 0:
            if time.time() - t0 > wait_s: raise RuntimeError(f"EVICTED and not downloaded in {wait_s} s: {path}")
            time.sleep(2)
    return path


def sha256(path, n=None):
    h = hashlib.sha256()
    with open(hydrated(path), "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""): h.update(chunk)
    return h.hexdigest() if n is None else h.hexdigest()[:n]


def now():
    return time.strftime("%Y-%m-%d %H:%M:%S")


def page_box(pdf):
    """pypdf mediabox of page 1 (the brief's rule), no content interpretation."""
    from pypdf import PdfReader
    r = PdfReader(hydrated(pdf)); assert len(r.pages) == 1, len(r.pages)
    mb = r.pages[0].mediabox
    return float(mb.width), float(mb.height)


# ---------------------------------------------------------------- watchdog
def wd_wait_free(logdir=WDLOG):
    lock = f"{logdir}/.running"
    while os.path.exists(lock):
        try: pid = int(open(lock).read().strip() or 0)
        except ValueError: pid = 0
        if pid and subprocess.run(["kill", "-0", str(pid)], capture_output=True).returncode == 0: time.sleep(2)
        else: break


def wd(name, cmd, max_s=900, rss_kb=4000000, logdir=WDLOG):
    os.makedirs(logdir, exist_ok=True); wd_wait_free(logdir)
    env = dict(os.environ, WD_LOGDIR=logdir)
    r = subprocess.run([WD, name, str(max_s), str(rss_kb), "--", *[str(c) for c in cmd]], capture_output=True, text=True, env=env)
    st = {}
    try:
        for kv in open(f"{logdir}/{name}.status").read().split():
            k, v = kv.split("=", 1); st[k] = v
    except FileNotFoundError:
        st = {"rc": str(r.returncode), "reason": "no_status_file"}
    if st.get("rc") != "0" or st.get("reason") != "finished":
        tail = ""
        try: tail = open(f"{logdir}/{name}.out").read()[-3000:]
        except FileNotFoundError: pass
        raise RuntimeError(f"watchdog {name}: {st} :: {r.stdout[-500:]} :: {tail}")
    return st


def gs_render(pdf, png, dpi, name, aa=True):
    hydrated(pdf); os.makedirs(os.path.dirname(png), exist_ok=True)
    bits = "4" if aa else "1"
    wd(name, [GS, "-dSAFER", "-dBATCH", "-dNOPAUSE", "-dQUIET", "-sDEVICE=png16m", f"-r{dpi}", f"-dTextAlphaBits={bits}", f"-dGraphicsAlphaBits={bits}", f"-sOutputFile={png}", pdf])
    assert os.path.exists(png) and os.path.getsize(png) > 0, png
    return png


# ---------------------------------------------------------------- Ghostscript text census
_SPAN = re.compile(r'<span bbox="(\S+) (\S+) (\S+) (\S+)" font="([^"]*)" size="([^"]*)">(.*?)</span>', re.S)
_CHAR = re.compile(r'<char bbox="(\S+) (\S+) (\S+) (\S+)" c="(.*?)"/>')


def gs_text(pdf, out_xml, name):
    hydrated(pdf); os.makedirs(os.path.dirname(out_xml), exist_ok=True)
    wd(name, [GS, "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=txtwrite", "-dTextFormat=0", f"-sOutputFile={out_xml}", pdf])
    return parse_txtwrite(out_xml)


def parse_txtwrite(path):
    """Spans: text (stripped), raw, x (first printing char), x1, y (baseline, top-down), font, size, chars [(x, c)], rotated."""
    raw = open(path, encoding="utf-8", errors="replace").read(); out = []
    for m in _SPAN.finditer(raw):
        x0, y0, x1, y1 = [float(v) for v in m.group(1, 2, 3, 4)]
        chars = [(float(c[0]), float(c[2]), html.unescape(c[4]).replace("\xa0", " ")) for c in _CHAR.findall(m.group(7))]
        text = "".join(c[2] for c in chars)
        if not text.strip(): continue
        lead = len(text) - len(text.lstrip()); xs = chars[lead][0] if lead < len(chars) else x0
        rotated = (x1 - x0) < 0.5 and len(chars) > 1
        out.append(dict(text=text.strip(), raw=text, x=xs, x1=x1, y=y0, font=m.group(5).split("+")[-1], size=round(float(m.group(6)), 2), chars=chars, rotated=rotated))
    out.sort(key=lambda s: (s["y"], s["x"]))
    return out


def lines_of(spans):
    """Join the spans of one baseline (1.5 pt) in x order into lines (a space where a span carries one at the boundary or the gap
    exceeds 0.35 em); rotated spans are their own lines. Returns [(line_text, [(char, size), ...])]."""
    rows = {}; out = []
    for s in sorted(spans, key=lambda q: (q["y"], q["x"])):
        if s["rotated"]:
            out.append((s["text"], [(c, s["size"]) for c in s["text"]])); continue
        key = next((k for k in rows if abs(k - s["y"]) <= 1.5), None)
        if key is None: key = s["y"]; rows[key] = []
        rows[key].append(s)
    for k, ss in rows.items():
        ss = sorted(ss, key=lambda q: q["x"]); chars = []; prev = None
        for q in ss:
            if prev is None: sep = ""
            elif prev["raw"].endswith(" ") or q["raw"].startswith(" "): sep = " "
            else: sep = "" if (q["x"] - prev["x1"]) < 0.35 * q["size"] else " "
            if sep: chars.append((" ", None))
            chars += [(c, q["size"]) for c in q["raw"]]; prev = q
        text = "".join(c for c, _s in chars)
        out.append((" ".join(text.split()), chars))
    return out


def strings_of(spans, split_em=0.36):
    """Strings as drawn: the spans of one baseline joined across kerning splits and carried spaces, but split where the gap to
    the next span exceeds split_em (a separate text artist, or another panel on a coincident baseline: a kerning split is far
    under 0.35 em and a space inside one artist is carried as a character). Rotated spans stand alone."""
    rows = {}; out = []
    for s in sorted(spans, key=lambda q: (q["y"], q["x"])):
        if s["rotated"]: out.append(s["text"]); continue
        key = next((k for k in rows if abs(k - s["y"]) <= 1.5), None)
        if key is None: key = s["y"]; rows[key] = []
        rows[key].append(s)
    for k, ss in rows.items():
        ss = sorted(ss, key=lambda q: q["x"]); cur = ""; prev = None
        for q in ss:
            if prev is None: cur = q["raw"]
            elif prev["raw"].endswith(" ") or q["raw"].startswith(" "): cur += q["raw"] if (prev["raw"].endswith(" ") or q["raw"].startswith(" ")) else " " + q["raw"]
            elif (q["x"] - prev["x1"]) < 0.35 * q["size"]: cur += q["raw"]
            elif (q["x"] - prev["x1"]) < split_em * q["size"]: cur += " " + q["raw"]
            else:
                out.append(" ".join(cur.split())); cur = q["raw"]
            prev = q
        if cur.strip(): out.append(" ".join(cur.split()))
    return [t for t in out if t]


def string_multiset(spans):
    return Counter(strings_of(spans))


def word_multiset(spans):
    """Counter of (word, size) over the joined lines; a word takes the size of its first character."""
    c = Counter()
    for _text, chars in lines_of(spans):
        word, size = "", None
        for ch, s in chars:
            if ch.isspace():
                if word: c[(word, size)] += 1
                word, size = "", None
            else:
                if not word: size = s
                word += ch
        if word: c[(word, size)] += 1
    return c


def char_multiset(spans):
    c = Counter()
    for s in spans:
        for _x, _x1, ch in s["chars"]:
            if not ch.isspace(): c[(ch, s["size"])] += 1
    return c


def line_multiset(spans):
    return Counter(t for t, _c in lines_of(spans))


def expected_plus_one(old, kept):
    """Expected NEW multiset: every (item, size) of OLD at size + 1, except the counts listed in kept (a Counter of (item, size))."""
    exp = Counter()
    for (w, s), n in old.items():
        k = min(n, kept.get((w, s), 0))
        if k: exp[(w, s)] += k
        if n - k: exp[(w, round(s + PT_PLUS, 2))] += n - k
    return exp


def kept_counters(kept_strings):
    """(word, size) and (char, size) Counters of the strings that keep their size: [(string, size, multiplicity)]."""
    w, c = Counter(), Counter()
    for text, size, mult in kept_strings:
        for _ in range(mult):
            for tok in text.split(): w[(tok, size)] += 1
            for ch in text:
                if not ch.isspace(): c[(ch, size)] += 1
    return w, c


def fmt_delta(a, b, n=12):
    """Pretty difference of two Counters: what a has beyond b and what b has beyond a."""
    return f"missing {sorted((a - b).items())[:n]}; extra {sorted((b - a).items())[:n]}"


# ---------------------------------------------------------------- raster helpers (PIL and numpy only)
def colour_mask(png, colours, tol=8):
    import numpy as np
    from PIL import Image
    a = np.asarray(Image.open(png).convert("RGB")).astype(np.int16); m = np.zeros(a.shape[:2], bool)
    for h in colours:
        r, g, b = [int(h[i:i + 2], 16) for i in (1, 3, 5)]
        m |= (np.abs(a[..., 0] - r) <= tol) & (np.abs(a[..., 1] - g) <= tol) & (np.abs(a[..., 2] - b) <= tol)
    return m


def marks_unchanged(old_png, new_png, colours, dpi, allowed_boxes_pt, tol=8, dilate_pt=1.0):
    """Pixels of the mark colours on two renders WITHOUT anti-aliasing must coincide outside the allowed boxes (legends)."""
    import numpy as np
    a, b = colour_mask(old_png, colours, tol), colour_mask(new_png, colours, tol)
    if a.shape != b.shape: return dict(shape_equal=False, shapes=[a.shape, b.shape], n_diff=None, n_outside=None)
    d = a ^ b; mask = np.zeros_like(d); z = dpi / 72.0
    for (x0, y0, x1, y1) in allowed_boxes_pt:
        X0, Y0 = max(int((x0 - dilate_pt) * z), 0), max(int((y0 - dilate_pt) * z), 0)
        X1, Y1 = min(int((x1 + dilate_pt) * z) + 1, d.shape[1]), min(int((y1 + dilate_pt) * z) + 1, d.shape[0])
        mask[Y0:Y1, X0:X1] = True
    ys, xs = np.nonzero(d & ~mask)
    bbox = None if len(ys) == 0 else [round(xs.min() / z, 1), round(ys.min() / z, 1), round(xs.max() / z, 1), round(ys.max() / z, 1)]
    return dict(shape_equal=True, n_old=int(a.sum()), n_new=int(b.sum()), n_diff=int(d.sum()), n_outside=int((d & ~mask).sum()), outside_bbox_pt=bbox)


def diff_components(old_png, new_png, colours, dpi, allowed_boxes_pt, tol=8, dilate_pt=1.0):
    """Connected components of the mark-colour pixels that differ between the two renders outside the allowed boxes:
    [(x0_pt, y0_pt, x1_pt, y1_pt, n_px)], so a moved mark (a large component) is told from rasterization flicker (1 to 2 px)."""
    import numpy as np
    from scipy import ndimage
    a, b = colour_mask(old_png, colours, tol), colour_mask(new_png, colours, tol); d = a ^ b; mask = np.zeros_like(d); z = dpi / 72.0
    for (x0, y0, x1, y1) in allowed_boxes_pt:
        X0, Y0 = max(int((x0 - dilate_pt) * z), 0), max(int((y0 - dilate_pt) * z), 0)
        X1, Y1 = min(int((x1 + dilate_pt) * z) + 1, d.shape[1]), min(int((y1 + dilate_pt) * z) + 1, d.shape[0])
        mask[Y0:Y1, X0:X1] = True
    lab, n = ndimage.label(d & ~mask); out = []
    for i, o in enumerate(ndimage.find_objects(lab)):
        out.append((round(o[1].start / z, 1), round(o[0].start / z, 1), round(o[1].stop / z, 1), round(o[0].stop / z, 1), int((lab[o] == i + 1).sum())))
    return out


def region_diff(png_old, png_new, dpi, box_pt, allowed_boxes_pt, thresh=32, dilate_pt=1.5):
    """The same page-point box on two anti-aliased renders (top-left anchored, so equal page coordinates even when the page
    heights differ): differing pixels (any channel beyond thresh) outside the allowed boxes, plus their bounding box in page pt."""
    import numpy as np
    from PIL import Image
    z = dpi / 72.0
    def crop(png):
        im = Image.open(png).convert("RGB"); x0, y0, x1, y1 = box_pt
        return np.asarray(im.crop((int(round(x0 * z)), int(round(y0 * z)), int(round(x1 * z)), int(round(y1 * z))))).astype(np.int16)
    a, b = crop(png_old), crop(png_new); h, w = min(a.shape[0], b.shape[0]), min(a.shape[1], b.shape[1]); a, b = a[:h, :w], b[:h, :w]
    d = (np.abs(a - b) > thresh).any(axis=2); mask = np.zeros_like(d)
    for (x0, y0, x1, y1) in allowed_boxes_pt:
        X0, Y0 = max(int((x0 - box_pt[0] - dilate_pt) * z), 0), max(int((y0 - box_pt[1] - dilate_pt) * z), 0)
        X1, Y1 = min(int((x1 - box_pt[0] + dilate_pt) * z) + 1, w), min(int((y1 - box_pt[1] + dilate_pt) * z) + 1, h)
        mask[Y0:Y1, X0:X1] = True
    out = d & ~mask; ys, xs = np.nonzero(out)
    bbox = None if len(ys) == 0 else [round(box_pt[0] + xs.min() / z, 1), round(box_pt[1] + ys.min() / z, 1), round(box_pt[0] + xs.max() / z, 1), round(box_pt[1] + ys.max() / z, 1)]
    return dict(n_diff=int(d.sum()), n_outside=int(out.sum()), size=[w, h], outside_bbox_pt=bbox)


def blobs(png, dpi, box_pt, colour_hex, min_h_px=8, tol=8):
    import numpy as np
    from PIL import Image
    from scipy import ndimage
    im = np.asarray(Image.open(png).convert("RGB")).astype(int); z = dpi / 72.0
    x0, y0, x1, y1 = [int(round(v * z)) for v in box_pt]; sub = im[y0:y1, x0:x1]
    r, g, b = [int(colour_hex[i:i + 2], 16) for i in (1, 3, 5)]
    mask = (np.abs(sub[..., 0] - r) <= tol) & (np.abs(sub[..., 1] - g) <= tol) & (np.abs(sub[..., 2] - b) <= tol)
    lab, n = ndimage.label(mask); out = []
    for o in ndimage.find_objects(lab): out.append((int(o[1].stop - o[1].start), int(o[0].stop - o[0].start)))
    return [wh for wh in out if wh[1] >= min_h_px], out


def colour_census(png, min_px=4):
    import numpy as np
    from PIL import Image
    a = np.asarray(Image.open(png).convert("RGB")).reshape(-1, 3)
    vals, counts = np.unique(a, axis=0, return_counts=True)
    return {"#%02x%02x%02x" % tuple(int(v) for v in row): int(c) for row, c in zip(vals, counts) if c >= min_px}


def side_by_side(png_a, png_b, out_png, gap=24):
    from PIL import Image
    a, b = Image.open(png_a).convert("RGB"), Image.open(png_b).convert("RGB"); h = max(a.height, b.height)
    out = Image.new("RGB", (a.width + b.width + gap, h), (255, 255, 255)); out.paste(a, (0, 0)); out.paste(b, (a.width + gap, 0)); out.save(out_png)


def crop_png(png, dpi, box_pt, out_png, scale=1):
    from PIL import Image
    im = Image.open(png).convert("RGB"); z = dpi / 72.0; x0, y0, x1, y1 = box_pt
    c = im.crop((max(0, int(x0 * z)), max(0, int(y0 * z)), min(im.width, int(x1 * z)), min(im.height, int(y1 * z))))
    if scale != 1: c = c.resize((c.width * scale, c.height * scale), Image.LANCZOS)
    os.makedirs(os.path.dirname(out_png), exist_ok=True); c.save(out_png); return c.size


# ---------------------------------------------------------------- checks writer
class Checks:
    def __init__(self, header):
        self.lines = [header, ""]; self.n_pass = 0; self.n_fail = 0

    def log(self, what, ok, detail=""):
        self.lines.append(f"{'PASS' if ok else 'FAIL'}: {what}" + (f" [{detail}]" if detail != "" else ""))
        if ok: self.n_pass += 1
        else: self.n_fail += 1
        print(self.lines[-1][:300], flush=True); return bool(ok)

    def info(self, what):
        self.lines.append(f"INFO: {what}"); print(self.lines[-1][:300], flush=True)

    def write(self, path):
        ok = self.n_fail == 0
        self.lines += ["", f"{self.n_pass} pass, {self.n_fail} fail", "RESULT ALL PASS" if ok else "RESULT FAIL"]
        os.makedirs(os.path.dirname(path), exist_ok=True); open(path, "w").write("\n".join(self.lines) + "\n"); return ok
