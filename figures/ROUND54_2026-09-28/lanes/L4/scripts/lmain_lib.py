#!$T90_PY
"""Lane LMAIN (round 40, 2026-09-18): shared helpers. Every render of a composed sheet goes through Ghostscript under the round-38
watchdog (wd_run.sh, WD_LOGDIR = this lane's logs/watchdog, one command at a time); text layers of composed sheets come from
Ghostscript txtwrite, never from PyMuPDF (PyMuPDF only for placement in a child process, and for the small flat pieces).
Paths are absolute. Before reading any file: blocks() > 0 (0 blocks with a nonzero size = iCloud-evicted)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import csv, hashlib, json, os, re, subprocess, time, collections
from decimal import Decimal, ROUND_HALF_UP

T90 = paths.FIGURE_ROOT
R37 = f"{paths.FIGURE_ROOT}/ROUND37_2026-09-14"; R38 = f"{paths.FIGURE_ROOT}/ROUND38_2026-09-16"; R40 = f"{paths.FIGURE_ROOT}/ROUND40_2026-09-18"
R49 = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26"; R54 = f"{paths.FIGURE_ROOT}/ROUND54_2026-09-28"; LANE = f"{R54}/lanes/L4"; WD = f"{R38}/figures/scripts/wd_run.sh"; WD_LOGDIR = f"{LANE}/logs/watchdog"   # round 49: this lane
V16 = f"{R38}/figures/NEW_FINAL_SET_V16"; V16_PNG = f"{R38}/figures/word_png_v16"; V13 = f"{paths.FIGURE_ROOT}/ROUND36_2026-09-10/L9_assembly/NEW_FINAL_SET_V13"
SKETCH = f"{R54}/lanes/LSKETCH/work/build"; SKETCH_STATUS = f"{R54}/lanes/LSKETCH/BANDS_READY.txt"   # round 54: the coordinator builds the bands here
GS = paths.GS; PY = paths.PY; TESS = paths.TESSERACT
ARIAL = paths.FONT_ARIAL; ARIALB = paths.FONT_ARIAL_BOLD
INK_RGB = (0x1a / 255, 0x1d / 255, 0x21 / 255)
TITLE_SIZE, TITLE_BASELINE = 18.0, 18.0   # round 54: the sheet title 18 pt on baseline 18 (cap top about 5 pt inside the page, as Figures 5 and 6)


# ------------------------------------------------------------------ files
def blocks(path):
    return int(subprocess.run(["stat", "-f", "%b", path], capture_output=True, text=True).stdout.strip() or 0)


def hydrated(path, wait_s=600):
    """0 blocks with a nonzero size = evicted: brctl download and poll."""
    if not os.path.exists(path): raise FileNotFoundError(path)
    if os.path.isdir(path) or os.path.getsize(path) == 0 or blocks(path) > 0: return path
    subprocess.run(["brctl", "download", path], check=False); t0 = time.time()
    while blocks(path) == 0:
        if time.time() - t0 > wait_s: raise RuntimeError(f"still evicted after {wait_s}s: {path}")
        time.sleep(5)
    return path


def sha256(path):
    h = hashlib.sha256()
    with open(hydrated(path), "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""): h.update(chunk)
    return h.hexdigest()


def load_json(path): return json.load(open(hydrated(path)))
def dump_json(obj, path): json.dump(obj, open(path, "w"), indent=1, default=str); return path
def now(): return time.strftime("%Y-%m-%d %H:%M:%S")


# ------------------------------------------------------------------ watchdog, Ghostscript
def wd(name, cmd, max_s=900, rss_kb=4000000):
    """Run ONE command under wd_run.sh (this lane's log dir). Returns (rc, status line). Refuses while another wd command runs here."""
    os.makedirs(WD_LOGDIR, exist_ok=True)
    env = dict(os.environ, WD_LOGDIR=WD_LOGDIR)
    r = subprocess.run([WD, name, str(max_s), str(rss_kb), "--"] + list(cmd), env=env, capture_output=True, text=True)
    st = open(f"{WD_LOGDIR}/{name}.status").read().strip() if os.path.exists(f"{WD_LOGDIR}/{name}.status") else r.stdout.strip()
    return r.returncode, st


def wd_out(name): return open(f"{WD_LOGDIR}/{name}.out").read() if os.path.exists(f"{WD_LOGDIR}/{name}.out") else ""


def gs_render(pdf, png, dpi, name, max_s=900, first_page_only=True):
    """png16m render through Ghostscript under the watchdog. Returns the status line; asserts rc 0, a FRESH output file (the target is
    removed first and must exist afterwards: a Ghostscript 'Page drawing error' exits 0 and leaves no page, which would otherwise hand a
    stale render to the next step, the 15:09 Fig 5 lesson) and no 'Error' line in Ghostscript's output."""
    hydrated(pdf)
    cmd = [GS, "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE"] + (["-dFirstPage=1", "-dLastPage=1"] if first_page_only else []) + ["-sDEVICE=png16m", f"-r{dpi}", f"-sOutputFile={png}", pdf]
    last = None
    for attempt in (1, 2):      # one retry: the 15:09 'Page drawing error' on Fig 5 did not recur on the next run of the same file
        if os.path.exists(png): os.remove(png)
        t0 = time.time(); rc, st = wd(name if attempt == 1 else f"{name}_retry", cmd, max_s=max_s); out = wd_out(name if attempt == 1 else f"{name}_retry")
        fresh = rc == 0 and os.path.exists(png) and os.path.getsize(png) > 0 and os.path.getmtime(png) >= t0 - 1
        if fresh and "error" not in out.lower(): return st + (" (attempt 2)" if attempt == 2 else "")
        last = ("gs render failed, wrote no fresh file, or reported an error", name, attempt, st, out[-500:]); print("WARN", last, flush=True)
    raise AssertionError(last)


def gs_pdfinfo(pdf):
    r = subprocess.run([GS, "-q", "-dNODISPLAY", "-dNOSAFER", "-dBATCH", "-dNOPAUSE", "-dPDFINFO", hydrated(pdf)], capture_output=True, text=True)
    txt = r.stdout + r.stderr
    boxes = re.findall(r"MediaBox:\s*\[([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\]", txt)
    m = re.search(r"Pages:\s*(\d+)", txt); n = int(m.group(1)) if m else len(boxes)
    assert boxes, ("no MediaBox in gs PDFINFO output", txt[-400:])
    x0, y0, x1, y1 = (float(v) for v in boxes[0]); return n, round(x1 - x0, 3), round(y1 - y0, 3)


def gs_txt(pdf, txt, name, max_s=600):
    """Text layer of a (flat or lightly composed) PDF through Ghostscript txtwrite under the watchdog. Never for the nested V13 sheets."""
    hydrated(pdf)
    rc, st = wd(name, [GS, "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=txtwrite", f"-sOutputFile={txt}", pdf], max_s=max_s)
    assert rc == 0 and os.path.exists(txt), ("gs txtwrite failed", name, st)
    return open(txt, encoding="utf-8", errors="replace").read()


# ------------------------------------------------------------------ multisets
NUMRE = re.compile(r"\d[\d,]*\.?\d*")


def norm_text(s):
    return s.replace("\xa0", " ").replace("−", "-").replace("–", "-").replace("—", "-").replace("\xad", "-")


def multisets(strings):
    words = collections.Counter(); nums = collections.Counter()
    for s in strings:
        s = norm_text(str(s))
        for w in s.split(): words[w] += 1
        for t in NUMRE.findall(s): nums[t] += 1
    return words, nums


def delta(a, b):
    return {"lost": dict(a - b), "gained": dict(b - a)}


def txt_strings(txt):
    """Lines of a gs txtwrite output as strings (collapsed whitespace), empty lines dropped."""
    out = []
    for line in txt.splitlines():
        line = re.sub(r"\s+", " ", norm_text(line)).strip()
        if line: out.append(line)
    return out


# ------------------------------------------------------------------ raster diff and crops
def raster_diff(png_a, png_b, allowed_boxes_pt, dpi, thresh=32, dilate_pt=1.5):
    """Differing pixels (any channel > thresh) outside the allowed page-point boxes; the two renders must share a shape unless
    the pages differ in height (then the common top part is compared and the extra rows are reported)."""
    import numpy as np
    from PIL import Image
    a = np.asarray(Image.open(png_a).convert("RGB")).astype(int); b = np.asarray(Image.open(png_b).convert("RGB")).astype(int)
    same_shape = a.shape == b.shape
    h = min(a.shape[0], b.shape[0]); w = min(a.shape[1], b.shape[1])
    d = (np.abs(a[:h, :w] - b[:h, :w]) > thresh).any(axis=2); mask = np.zeros_like(d); s = dpi / 72.0
    for (x0, y0, x1, y1) in allowed_boxes_pt:
        X0, Y0 = max(int((x0 - dilate_pt) * s), 0), max(int((y0 - dilate_pt) * s), 0); X1, Y1 = min(int((x1 + dilate_pt) * s) + 1, w), min(int((y1 + dilate_pt) * s) + 1, h)
        mask[Y0:Y1, X0:X1] = True
    out = dict(n_diff=int(d.sum()), n_outside=int((d & ~mask).sum()), shape_equal=same_shape, shape_a=list(a.shape), shape_b=list(b.shape), compared_px=[h, w])
    ys, xs = np.nonzero(d & ~mask)
    if len(ys): out["outside_bbox_pt"] = [round(float(xs.min()) / s, 1), round(float(ys.min()) / s, 1), round(float(xs.max()) / s, 1), round(float(ys.max()) / s, 1)]
    return out


def crop_pair(old_png, new_png, rect_pt, dpi, out_stem, label_old="OLD (V16)", label_new="NEW (round 40)"):
    """Side by side OLD vs NEW crop from two full-page renders at the same dpi (PIL only, never a clipped re-render)."""
    from PIL import Image, ImageDraw
    S = dpi / 72.0; x0, y0, x1, y1 = [int(round(v * S)) for v in rect_pt]; ims = []
    for p in (old_png, new_png):
        im = Image.open(p).convert("RGB"); ims.append(im.crop((max(x0, 0), max(y0, 0), min(x1, im.width), min(y1, im.height))))
    h = max(i.height for i in ims) + 26; w = sum(i.width for i in ims) + 30
    canvas = Image.new("RGB", (w, h), (255, 255, 255)); dr = ImageDraw.Draw(canvas); x = 10
    for im, lab in zip(ims, (label_old, label_new)):
        dr.text((x, 4), lab, fill=(200, 0, 0)); canvas.paste(im, (x, 22)); x += im.width + 10
    canvas.save(out_stem + "_OLD_vs_NEW_200dpi.png"); return out_stem + "_OLD_vs_NEW_200dpi.png"


# ------------------------------------------------------------------ OCR
def ocr_words(png):
    """tesseract --psm 11 word boxes (pixels) of a render."""
    base = png[:-4] + "_ocr"
    r = subprocess.run([TESS, png, base, "--psm", "11", "-l", "eng", "tsv"], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-400:]
    words = []
    for row in csv.DictReader(open(base + ".tsv"), delimiter="\t", quoting=csv.QUOTE_NONE):
        if row.get("level") != "5" or not (row.get("text") or "").strip(): continue
        L, T, Wd, Hd = (int(row[k]) for k in ("left", "top", "width", "height"))
        words.append(dict(text=row["text"], x0=L, y0=T, x1=L + Wd, y1=T + Hd, cx=L + Wd / 2, cy=T + Hd / 2, conf=float(row["conf"])))
    return words


def ocr_crop(png, box_pt, dpi, psm=7, upscale=3, whitelist=None):
    """tesseract on a crop of a render around a page-point box (upscaled, padded); the first non-empty reading of: the crop as is,
    contrast-stretched, binarized (dark ink on a saturated fill, e.g. the green header cells), and the inverted crop (white on dark)."""
    from PIL import Image, ImageOps
    import tempfile
    Z = dpi / 72.0; im = Image.open(png).convert("L")
    x0, y0, x1, y1 = box_pt; c = im.crop((max(0, int(x0 * Z)), max(0, int(y0 * Z)), int(x1 * Z), int(y1 * Z)))
    if c.width < 4 or c.height < 4: return ""
    c = c.resize((c.width * upscale, c.height * upscale), Image.LANCZOS)
    variants = [c, ImageOps.autocontrast(c, cutoff=2), ImageOps.autocontrast(c, cutoff=2).point(lambda v: 255 if v > 128 else 0), ImageOps.invert(c)]
    cmd0 = [TESS, "", "stdout", "--psm", str(psm), "-l", "eng"] + (["-c", f"tessedit_char_whitelist={whitelist}"] if whitelist else [])
    for v in variants:
        v = ImageOps.expand(v, border=24, fill=255); p = tempfile.mktemp(suffix=".png", dir=f"{LANE}/logs"); v.save(p)
        cmd = list(cmd0); cmd[1] = p; r = subprocess.run(cmd, capture_output=True, text=True); os.remove(p)
        got = r.stdout.strip().replace("\n", " ")
        if got: return got
    return ""


# ------------------------------------------------------------------ checks and CHANGES
class Checks:
    def __init__(self, header): self.lines = [header, ""]; self.n_pass = 0; self.n_fail = 0
    def log(self, what, ok, detail=""):
        ok = bool(ok); self.lines.append(f"{'PASS' if ok else 'FAIL'}: {what}" + (f" [{detail}]" if detail != "" else ""))
        self.n_pass += ok; self.n_fail += (not ok); print(self.lines[-1][:400], flush=True); return ok
    def info(self, text): self.lines.append(f"INFO: {text}"); print(self.lines[-1][:400], flush=True)
    def write(self, path):
        ok = self.n_fail == 0
        self.lines += ["", f"{self.n_pass} pass, {self.n_fail} fail", "RESULT ALL PASS" if ok else "RESULT FAIL"]
        os.makedirs(os.path.dirname(path), exist_ok=True); open(path, "w").write("\n".join(self.lines) + "\n"); print(self.lines[-1]); return ok


CHANGES_COLS = ["sheet", "panel", "item", "old", "new", "kind", "x", "y", "note"]


def write_changes(path, rows):
    with open(path, "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=CHANGES_COLS, extrasaction="ignore"); wr.writeheader()
        for r in rows: wr.writerow({k: r.get(k, "") for k in CHANGES_COLS})
    return len(rows)


def status(sheet_rows, note=""):
    """Rewrite STATUS.md (percent per sheet)."""
    lines = ["# LMAIN STATUS (round 40)", f"Updated {now()} EDT.", "", "| sheet | percent | state |", "|---|---|---|"]
    for s, p, st in sheet_rows: lines.append(f"| {s} | {p} | {st} |")
    if note: lines += ["", note]
    open(f"{LANE}/STATUS.md", "w").write("\n".join(lines) + "\n")


def sketch_ready(name):
    """True when LS_SKETCH/STATUS.md carries 'BAND READY <name>' and the file exists with blocks."""
    if not os.path.exists(SKETCH_STATUS): return False
    txt = open(SKETCH_STATUS).read()
    ok = re.search(r"BAND READY\s+" + re.escape(name) + r"\b", txt) is not None
    p = f"{SKETCH}/{name}.pdf"
    return ok and os.path.exists(p) and blocks(p) > 0


def halfup(v, nd): return str(Decimal(repr(float(v))).quantize(Decimal(1).scaleb(-nd), rounding=ROUND_HALF_UP))
