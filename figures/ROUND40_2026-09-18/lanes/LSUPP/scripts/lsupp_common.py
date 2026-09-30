"""Lane LSUPP (round 40, 2026-09-18): shared helpers for the eight supplementary sheets.

Rules of the round folded in:
  - every render is Ghostscript (png16m) run ONE at a time under the round-38 watchdog wd_run.sh (WD_LOGDIR = this lane's logs/watchdog);
  - the text layer of a sheet is read with Ghostscript's txtwrite device (TextFormat 0: spans with positions), never with PyMuPDF
    get_text / pdftotext, so composed sheets are never opened by PyMuPDF for extraction;
  - PyMuPDF is used for PLACEMENT only (show_pdf_page, TextWriter, mediabox, redaction) in child scripts run under the watchdog;
  - every file is checked for iCloud eviction (stat -f %b nonzero) before it is read;
  - OCR read-back = Ghostscript 300 dpi render + tesseract (the round-38 idiom).
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import csv, hashlib, html, json, os, re, subprocess, time
from collections import Counter

T90 = paths.FIGURE_ROOT
R40 = f"{paths.FIGURE_ROOT}/ROUND40_2026-09-18"
LANE = f"{R40}/lanes/LSUPP"
V16 = f"{paths.FIGURE_ROOT}/ROUND38_2026-09-16/figures/NEW_FINAL_SET_V16"
R37 = f"{paths.FIGURE_ROOT}/ROUND37_2026-09-14/figures"
R38 = f"{paths.FIGURE_ROOT}/ROUND38_2026-09-16/figures"
WD = f"{R38}/scripts/wd_run.sh"
GS = paths.GS
PY = paths.PY
TESS = paths.TESSERACT
WDLOG = f"{LANE}/logs/watchdog"
os.makedirs(WDLOG, exist_ok=True)
INK = "#1a1d21"
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
    """Page width and height in points from the PDF's own MediaBox (regex on the bytes, no content interpretation);
    falls back to PyMuPDF's page rect (opening a document reads no content stream)."""
    b = open(hydrated(pdf), "rb").read()
    m = re.findall(rb"/MediaBox\s*\[\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*\]", b)
    if m:
        x0, y0, x1, y1 = [float(v) for v in m[0]]
        return round(x1 - x0, 3), round(y1 - y0, 3)
    import fitz
    d = fitz.open(pdf); r = d[0].rect; n = len(d); d.close()
    assert n == 1, n
    return round(r.width, 3), round(r.height, 3)


# ---------------------------------------------------------------- the watchdog
def wd_wait_free(logdir=WDLOG):
    lock = f"{logdir}/.running"
    while os.path.exists(lock):
        try:
            pid = int(open(lock).read().strip() or 0)
        except ValueError:
            pid = 0
        if pid and subprocess.run(["kill", "-0", str(pid)], capture_output=True).returncode == 0:
            time.sleep(2)
        else:
            break


def wd(name, cmd, max_s=900, rss_kb=4000000, logdir=WDLOG):
    """Run ONE command under wd_run.sh (RSS and wall-time capped), waiting first for the lane's watchdog slot to be free.
    Raises when the command fails or the watchdog killed it; returns the parsed .status line."""
    os.makedirs(logdir, exist_ok=True)
    wd_wait_free(logdir)
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
        try:
            tail = open(f"{logdir}/{name}.out").read()[-3000:]
        except FileNotFoundError:
            pass
        raise RuntimeError(f"watchdog {name}: {st} :: {r.stdout[-500:]} :: {tail}")
    return st


def gs_render(pdf, png, dpi, name, max_s=900):
    """Ghostscript png16m render of page 1 (the only page), anti-aliased text and graphics, under the watchdog."""
    hydrated(pdf)
    os.makedirs(os.path.dirname(png), exist_ok=True)
    st = wd(name, [GS, "-dSAFER", "-dBATCH", "-dNOPAUSE", "-dQUIET", "-sDEVICE=png16m", f"-r{dpi}", "-dTextAlphaBits=4",
                    "-dGraphicsAlphaBits=4", f"-sOutputFile={png}", pdf], max_s=max_s)
    assert os.path.exists(png) and os.path.getsize(png) > 0, png
    return st


# ---------------------------------------------------------------- text layer through Ghostscript txtwrite
_SPAN = re.compile(r'<span bbox="(\S+) (\S+) (\S+) (\S+)" font="([^"]*)" size="([^"]*)">(.*?)</span>', re.S)
_CHAR = re.compile(r'<char bbox="(\S+) (\S+) (\S+) (\S+)" c="(.*?)"/>')


def gs_text(pdf, out_txt, name):
    """Spans of the sheet from Ghostscript's txtwrite device (TextFormat 0): text, x0, x1, y (baseline, top-down page points,
    integer precision), font, size. The device concatenates the characters of one PDF text run into a span."""
    hydrated(pdf)
    os.makedirs(os.path.dirname(out_txt), exist_ok=True)
    wd(name, [GS, "-dSAFER", "-dBATCH", "-dNOPAUSE", "-dQUIET", "-sDEVICE=txtwrite", "-dTextFormat=0", f"-sOutputFile={out_txt}", pdf])
    return parse_txtwrite(out_txt)


def parse_txtwrite(path):
    raw = open(path, encoding="utf-8", errors="replace").read()
    out = []
    for m in _SPAN.finditer(raw):
        x0, y0, x1, y1 = [float(v) for v in m.group(1, 2, 3, 4)]
        chars = [(float(c[0]), html.unescape(c[4])) for c in _CHAR.findall(m.group(7))]
        text = "".join(c[1] for c in chars).replace("\xa0", " ")
        if not text.strip():
            continue
        # trim leading/trailing blanks, keep the x of the first printing char
        lead = len(text) - len(text.lstrip())
        xs = chars[lead][0] if lead < len(chars) else x0
        out.append(dict(text=text.strip(), raw=text, x=xs, x1=x1, y=y0, font=m.group(5).split("+")[-1], size=float(m.group(6))))
    out.sort(key=lambda s: (s["y"], s["x"]))
    return out


def collapse(spans, tol=1.0):
    """Stacked duplicate copies (same text within tol pt) count once (the round-28 rule)."""
    keep = []
    for s in spans:
        if not any(s["text"] == k["text"] and abs(s["x"] - k["x"]) <= tol and abs(s["y"] - k["y"]) <= tol for k in keep):
            keep.append(s)
    return keep


def words_of(spans):
    return Counter(w for s in spans for w in s["text"].split())


def nums_of(spans):
    return Counter(t for s in spans for t in NUMRE.findall(s["text"]))


def multiset_delta(old, new):
    return dict(removed=dict(old - new), added=dict(new - old))


def counter_of_strings(strings):
    w, n = Counter(), Counter()
    for t in strings:
        w.update(t.split()); n.update(NUMRE.findall(t))
    return w, n


# ---------------------------------------------------------------- OCR (the round-38 idiom: gs 300 dpi + tesseract)
def ocr_page(png, psm=11):
    r = subprocess.run([TESS, png, "stdout", "--psm", str(psm), "-l", "eng"], capture_output=True, text=True)
    return r.stdout


def ocr_crop(png, box_pt, dpi, out_png, psm=7, up=3, pad_pt=3.0):
    """OCR of a targeted crop (page points -> pixels), upscaled; returns the read string."""
    from PIL import Image, ImageOps
    im = Image.open(png).convert("L")
    z = dpi / 72.0
    x0, y0, x1, y1 = box_pt
    c = im.crop((max(0, int((x0 - pad_pt) * z)), max(0, int((y0 - pad_pt) * z)), min(im.width, int((x1 + pad_pt) * z)), min(im.height, int((y1 + pad_pt) * z))))
    c = ImageOps.expand(c.resize((c.width * up, c.height * up), Image.LANCZOS), border=24, fill=255)
    c.save(out_png)
    r = subprocess.run([TESS, out_png, "stdout", "--psm", str(psm), "-l", "eng"], capture_output=True, text=True)
    return " ".join(r.stdout.split())


def norm_ocr(s):
    return re.sub(r"[^a-z0-9%<>=.()*+-]", "", s.lower().replace("—", "-").replace("–", "-").replace("−", "-"))


# ---------------------------------------------------------------- crops from Ghostscript renders (PIL only, no PyMuPDF pixmaps)
def crop_png(png, dpi, box_pt, out_png, label=None):
    from PIL import Image, ImageDraw
    im = Image.open(png).convert("RGB")
    z = dpi / 72.0
    x0, y0, x1, y1 = box_pt
    c = im.crop((max(0, int(x0 * z)), max(0, int(y0 * z)), min(im.width, int(x1 * z)), min(im.height, int(y1 * z))))
    if label:
        out = Image.new("RGB", (c.width, c.height + 22), (255, 255, 255)); out.paste(c, (0, 22))
        ImageDraw.Draw(out).text((4, 4), label, fill=(120, 120, 120)); c = out
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    c.save(out_png)
    return c.size


def side_by_side(png_a, png_b, out_png, gap=24):
    from PIL import Image
    a, b = Image.open(png_a).convert("RGB"), Image.open(png_b).convert("RGB")
    h = max(a.height, b.height)
    out = Image.new("RGB", (a.width + b.width + gap, h), (255, 255, 255))
    out.paste(a, (0, 0)); out.paste(b, (a.width + gap, 0))
    out.save(out_png)
    return out.size


def raster_diff(png_a, png_b, allowed_boxes_pt, dpi, thresh=32, dilate_pt=1.5):
    """Differing pixels (any channel beyond thresh) outside the allowed page-point boxes; shapes must match."""
    import numpy as np
    from PIL import Image
    a = np.asarray(Image.open(png_a).convert("RGB")).astype(int); b = np.asarray(Image.open(png_b).convert("RGB")).astype(int)
    if a.shape != b.shape:
        return dict(n_diff=None, n_outside=None, shape_equal=False, shapes=[a.shape, b.shape])
    d = (np.abs(a - b) > thresh).any(axis=2); mask = np.zeros_like(d); z = dpi / 72.0
    for (x0, y0, x1, y1) in allowed_boxes_pt:
        X0, Y0 = max(int((x0 - dilate_pt) * z), 0), max(int((y0 - dilate_pt) * z), 0)
        X1, Y1 = min(int((x1 + dilate_pt) * z) + 1, d.shape[1]), min(int((y1 + dilate_pt) * z) + 1, d.shape[0])
        mask[Y0:Y1, X0:X1] = True
    ys, xs = np.nonzero(d & ~mask)
    bbox = None if len(ys) == 0 else [round(xs.min() / z, 1), round(ys.min() / z, 1), round(xs.max() / z, 1), round(ys.max() / z, 1)]
    return dict(n_diff=int(d.sum()), n_outside=int((d & ~mask).sum()), shape_equal=True, outside_bbox_pt=bbox)


# ---------------------------------------------------------------- checks and CHANGES
class Checks:
    def __init__(self, header):
        self.lines = [header, ""]; self.n_pass = 0; self.n_fail = 0

    def log(self, what, ok, detail=""):
        self.lines.append(f"{'PASS' if ok else 'FAIL'}: {what}" + (f" [{detail}]" if detail != "" else ""))
        if ok: self.n_pass += 1
        else: self.n_fail += 1
        print(self.lines[-1][:400], flush=True); return bool(ok)

    def info(self, what):
        self.lines.append(f"INFO: {what}"); print(self.lines[-1][:400], flush=True)

    def write(self, path):
        ok = self.n_fail == 0
        self.lines += ["", f"{self.n_pass} pass, {self.n_fail} fail", "RESULT ALL PASS" if ok else "RESULT FAIL"]
        os.makedirs(os.path.dirname(path), exist_ok=True); open(path, "w").write("\n".join(self.lines) + "\n"); return ok


def write_changes(path, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=["sheet", "panel", "item", "before", "after", "dramatic", "note"], extrasaction="ignore"); wr.writeheader()
        for r in rows: wr.writerow({k: r.get(k, "") for k in wr.fieldnames})


def status_update(rows, note=""):
    """Rewrite the lane STATUS.md (percent per sheet)."""
    lines = ["# LSUPP STATUS (round 40)", f"Updated {now()}. {note}".rstrip(), "", "| sheet | percent | note |", "|---|---|---|"]
    for sheet, pct, n in rows:
        lines.append(f"| {sheet} | {pct} | {n} |")
    open(f"{LANE}/STATUS.md", "w").write("\n".join(lines) + "\n")


# ---------------------------------------------------------------- marker blobs on a non-anti-aliased render (one marker per key entry)
def blobs(png, dpi, box_pt, colour_hex, min_h_px=8, tol=8):
    """Connected components of one exact colour inside a page-point box of a render made WITHOUT anti-aliasing.
    Returns (components with height >= min_h_px as (w, h) pairs, all components). A 1.7 pt handle line is 3.5 px tall at 150 dpi,
    a 6 pt marker about 12 px, and the marker's white edge separates it from the line, so the tall components are the markers."""
    import numpy as np
    from PIL import Image
    from scipy import ndimage
    im = np.asarray(Image.open(png).convert("RGB")).astype(int); z = dpi / 72.0
    x0, y0, x1, y1 = [int(round(v * z)) for v in box_pt]; sub = im[y0:y1, x0:x1]
    r, g, b = [int(colour_hex[i:i + 2], 16) for i in (1, 3, 5)]
    mask = (np.abs(sub[..., 0] - r) <= tol) & (np.abs(sub[..., 1] - g) <= tol) & (np.abs(sub[..., 2] - b) <= tol)
    lab, n = ndimage.label(mask); out = []
    for o in ndimage.find_objects(lab):
        out.append((int(o[1].stop - o[1].start), int(o[0].stop - o[0].start)))
    return [wh for wh in out if wh[1] >= min_h_px], out
