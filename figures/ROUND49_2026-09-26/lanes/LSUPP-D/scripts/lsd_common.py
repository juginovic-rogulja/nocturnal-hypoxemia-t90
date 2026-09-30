"""Lane LSUPP-D (round 49, 2026-09-26): shared helpers. The six sheets of this lane (Supp 8, 9, 10, 16, 17, 19) are legacy
overlay sheets (Chromium or matplotlib bases of August 2026 with fitz overlays of rounds 29b to 38, grid surgery of round 31,
the key move of round 32 and the band placement of round 44), not page-pinned matplotlib builders, so the +1 pt is applied at
the PDF level: every VISIBLE string of the current sheet is removed (text-only redaction) and re-set with the sheet's own Arial
face one point larger on its own anchor (left edge, right edge, centre, or the centre of a rotated axis title), graphics untouched.

Rules of the round folded in: renders and text censuses are Ghostscript only (png16m, txtwrite); PyMuPDF runs ONLY in child
processes under the round-38 watchdog (wd_run.sh, RSS and wall-time capped) and only on FLAT single-page files (nested sheets are
re-distilled by Ghostscript pdfwrite first); every input is checked for iCloud eviction before it is read."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import hashlib, html, json, os, re, subprocess, time
from collections import Counter

T90 = paths.FIGURE_ROOT
R49 = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26"
LANE = f"{R49}/lanes/LSUPP-D"
SCRIPTS = f"{LANE}/scripts"
V26 = f"{paths.FIGURE_ROOT}/ROUND47_2026-09-24/figures/NEW_FINAL_SET_V26"
V15 = f"{paths.FIGURE_ROOT}/ROUND38_2026-09-16/figures/NEW_FINAL_SET_V15"          # vector twins of the V16 to V26 rasters of Supp 9 and 10
BAND_R49 = f"{R49}/lanes/LSKETCH/work/build/band_s16_r49.pdf"
BAND_R44 = f"{paths.FIGURE_ROOT}/ROUND44_2026-09-21/lanes/LSKETCH/work/build/band_s16_r44.pdf"
GS = paths.GS
PY = paths.PY
WD = f"{SCRIPTS}/wd_run.sh"
WDLOG = f"{LANE}/logs/watchdog"
INK = "#1a1d21"
os.makedirs(WDLOG, exist_ok=True)


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
    """Page count and the first page's MediaBox width and height (pypdf, no content interpretation)."""
    from pypdf import PdfReader
    r = PdfReader(hydrated(pdf)); p = r.pages[0]
    return len(r.pages), float(p.mediabox.width), float(p.mediabox.height)


# ---------------------------------------------------------------- watchdog and Ghostscript
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


def wd(name, cmd, max_s=600, rss_kb=4000000, logdir=WDLOG):
    """Run ONE command under wd_run.sh (RSS and wall-time capped). Raises when it fails or the watchdog killed it."""
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
    out = ""
    try:
        out = open(f"{logdir}/{name}.out").read()
    except FileNotFoundError:
        pass
    if st.get("rc") != "0" or st.get("reason") != "finished":
        raise RuntimeError(f"watchdog {name}: {st} :: {out[-3000:]}")
    return st, out


def gs_render(pdf, png, dpi=150):
    hydrated(pdf); os.makedirs(os.path.dirname(png), exist_ok=True)
    subprocess.run([GS, "-dSAFER", "-dBATCH", "-dNOPAUSE", "-dQUIET", "-sDEVICE=png16m", f"-r{dpi}", "-dTextAlphaBits=4", "-dGraphicsAlphaBits=4",
                    f"-sOutputFile={png}", pdf], check=True, capture_output=True)
    assert os.path.exists(png) and os.path.getsize(png) > 0, png
    return png


def gs_flatten(src, out):
    """Re-distil a nested sheet into a single-page pdfwrite file (no downsampling, colours untouched, fonts embedded), the round-44 recipe."""
    hydrated(src); os.makedirs(os.path.dirname(out), exist_ok=True)
    subprocess.run([GS, "-dSAFER", "-dBATCH", "-dNOPAUSE", "-dQUIET", "-sDEVICE=pdfwrite", "-dCompatibilityLevel=1.7", "-dPDFSETTINGS=/prepress",
                    "-dDownsampleColorImages=false", "-dDownsampleGrayImages=false", "-dDownsampleMonoImages=false", "-dAutoFilterColorImages=false",
                    "-dAutoFilterGrayImages=false", "-dColorImageFilter=/FlateEncode", "-dGrayImageFilter=/FlateEncode", "-dColorConversionStrategy=/LeaveColorUnchanged",
                    "-dEmbedAllFonts=true", "-dSubsetFonts=true", "-dPreserveAnnots=false", f"-sOutputFile={out}", src], check=True, capture_output=True)
    return out


_SPAN = re.compile(r'<span bbox="(\S+) (\S+) (\S+) (\S+)" font="([^"]*)" size="([^"]*)">(.*?)</span>', re.S)
_CHAR = re.compile(r'<char bbox="(\S+) (\S+) (\S+) (\S+)" c="(.*?)"/>')


def gs_txtwrite(pdf, xml):
    hydrated(pdf); os.makedirs(os.path.dirname(xml), exist_ok=True)
    subprocess.run([GS, "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=txtwrite", "-dTextFormat=0", f"-sOutputFile={xml}", pdf], check=True, capture_output=True)
    return parse_txtwrite(xml)


def parse_txtwrite(xml):
    """Spans of the Ghostscript txtwrite census: text, x0, x1, y (baseline, integer page points), font, size."""
    raw = open(xml, encoding="utf-8", errors="replace").read(); out = []
    for m in _SPAN.finditer(raw):
        x0, y0, x1, y1 = [float(v) for v in m.group(1, 2, 3, 4)]
        chars = [(float(c[0]), float(c[2]), html.unescape(c[4])) for c in _CHAR.findall(m.group(7))]
        text = "".join(c[2] for c in chars).replace("\xa0", " ")
        if not text.strip():
            continue
        lead = len(text) - len(text.lstrip()); trail = len(text) - len(text.rstrip())
        xs = chars[lead][0] if lead < len(chars) else x0
        xe = chars[len(chars) - 1 - trail][1] if trail < len(chars) else x1
        out.append(dict(text=text.strip(), raw=text, x0=xs, x1=xe, y=y0, font=m.group(5).split("+")[-1], size=float(m.group(6))))
    return out


def census_rows(spans, gap=1.5, ytol=1.5):
    """Normalised census rows. txtwrite splits one string into several runs (kerning, glyph-cluster and pdfwrite run breaks) and
    reports stacked copies and copies of earlier generations separately. A row = runs of one size on one baseline (within ytol)
    that ABUT in x (next.x0 - prev.x1 within +-gap), joined in x order; runs that overlap an earlier run start a new row (a copy,
    never glued onto the visible string); rows with the same text at the same x (1 pt) and baseline are stacked copies, collapsed
    to one. Returns (Counter of (text, size), the rows with their positions)."""
    items = sorted(spans, key=lambda s: (round(s["size"], 2), round(s["y"] / (2 * ytol)), s["x0"], s["x1"]))
    rows = []
    for s in items:
        r = rows[-1] if rows else None
        if r is not None and abs(r["size"] - s["size"]) <= 0.011 and abs(r["y"] - s["y"]) <= ytol and -gap <= s["x0"] - r["x1"] <= gap:
            r["parts"].append(s); r["x1"] = max(r["x1"], s["x1"])
        else:
            rows.append(dict(size=s["size"], y=s["y"], x0=s["x0"], x1=s["x1"], parts=[s]))
    out = Counter(); detail = []
    for r in rows:
        text = _join_parts(r["parts"]) if len(r["parts"]) > 1 else r["parts"][0]["text"]
        if any(d["text"] == text and abs(d["x0"] - r["x0"]) <= 1.0 and abs(d["y"] - r["y"]) <= ytol and abs(d["size"] - round(r["size"], 2)) <= 0.011 for d in detail):
            for d in detail:
                if d["text"] == text and abs(d["x0"] - r["x0"]) <= 1.0 and abs(d["y"] - r["y"]) <= ytol and abs(d["size"] - round(r["size"], 2)) <= 0.011:
                    d["n_copies"] += 1; break
            continue
        out[(text, round(r["size"], 2))] += 1
        detail.append(dict(text=text, size=round(r["size"], 2), y=r["y"], x0=r["x0"], x1=r["x1"], n_parts=len(r["parts"]), n_copies=1))
    return out, detail


def _join_parts(parts):
    """Join abutting runs in x order; a gap wider than a quarter of the size becomes one space unless a space is already there."""
    parts = sorted(parts, key=lambda p: p["x0"]); s = parts[0].get("raw", parts[0]["text"]); size = parts[0]["size"]
    for a, b in zip(parts, parts[1:]):
        t = b.get("raw", b["text"])
        if b["x0"] - a["x1"] > 0.25 * size and not s.endswith(" ") and not t.startswith(" "):
            s += " "
        s += t
    return s.strip()


def load_png(png):
    import numpy as np
    from PIL import Image
    return np.asarray(Image.open(png).convert("RGB")).astype(int)


class Checks:
    def __init__(self, header):
        self.lines = [header]; self.ok = True
    def log(self, name, ok, detail=""):
        self.ok &= bool(ok); self.lines.append(f"{'PASS' if ok else 'FAIL'}  {name}: {detail}"); print(self.lines[-1]); return ok
    def info(self, text):
        self.lines.append("INFO  " + text); print(self.lines[-1])
    def write(self, path):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        self.lines.append("RESULT " + ("ALL PASS" if self.ok else "FAIL: see FAIL lines"))
        open(path, "w").write("\n".join(self.lines) + "\n"); print(self.lines[-1]); return self.ok
