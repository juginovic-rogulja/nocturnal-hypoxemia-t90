#!/usr/bin/env python3
"""Lane V14_L1_RANK (round 37, v8.1 numbers, V13 design), shared helpers, repointed from ROUND30 V7_L1_RANK/scripts/common.py: paths, numbers files with sidecar checks, fonts, page-pinned text, text-layer
extraction, token multisets, font-metric alignment, letter stamping, raster diffs, crops, CHANGES rows."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import os, re, json, hashlib, gc, fcntl, sys
from collections import Counter
import fitz

T90 = paths.FIGURE_ROOT
R30 = f"{paths.FIGURE_ROOT}/ROUND30_2026-09-08"
LANE = f"{paths.FIGURE_ROOT}/ROUND54_2026-09-28/lanes/L1"   # round 54 lane L1 (2026-09-29): every text +3 pt on the LFIG1E panels (PT_PLUS 4.0), panels re-laid   # LFIG1E: the family key moved into the empty upper-middle of panel b, away from panel c   # round 47: panels b and c swapped, callouts 11 pt (Dragana, 2026-09-24)   # round 44c: limb movements folded into sleep timing and structure (Alen, 2026-09-21)
BASE = f"{paths.FIGURE_ROOT}/ROUND36_2026-09-10/L9_assembly/NEW_FINAL_SET_V13"   # V14: the V13 sheets are the design baseline
NUM = f"{paths.NUMBERS_DIR}"
SV = paths.SV_ROOT
OLDNUM = f"{paths.V8_ROOT}/compare/old_numbers"   # the v7 numbers snapshot: the OLD side of every delta only, never a live source
DATA_V8 = f"{paths.TABLES_DIR}"
ARIAL = paths.FONT_ARIAL
ARIALB = paths.FONT_ARIAL_BOLD
ARIALI = paths.FONT_ARIAL_ITALIC
RENDER_LOCK = f"{R30}/_render.lock"
INK = "#1a1d21"
GREY = "#8a9099"
GREY_PALE = "#ccd1d6"
BLUE = "#0288d1"
FAMILY_COL = {"Oxygenation": "#0288d1", "Heart rate and variability": "#8d2dfa", "Breathing events": "#d55e00",
              "Brain, microstructure": "#fa93a5", "Brain, spectral power": "#f2e279",
              "Sleep timing and structure": "#39c445"}          # sampled from the round-30 sheets (what shows); v8.1: the seventh family (PLM index) in the palette grey, an owner call
FAMILY_ORDER_TOP_DOWN = ["Oxygenation", "Heart rate and variability", "Breathing events", "Brain, microstructure",
                         "Brain, spectral power", "Sleep timing and structure"]
BANNED = ["—", ";", "strata", "stratum", "honest", "straightforward", "null", "pre-specified", "prespecified"]


def hydrated(path):
    """iCloud eviction gate: a file with size and 0 blocks is evicted, never read it."""
    st = os.stat(path)
    assert not (st.st_size > 0 and st.st_blocks == 0), f"EVICTED (0 blocks): {path}"
    return path


def sha256(path):
    h = hashlib.sha256()
    with open(hydrated(path), "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sidecar_ok(path, need_today=False, today=None):
    """v8.1 provenance gate (V14): the sidecar must exist, name data_frozen_v8_2026-09 inputs (or be a v8.1 step output dated at or
    after the v8.1 table build), carry rc 0 and match the file's sha256. Delegates to the shared v14lib.sidecar_v8; the R30 arguments
    need_today/today are kept for call compatibility and ignored (every v8.1 file is dated 2026-09-14 12:35 or later)."""
    import sys as _s
    _c = f"{paths.FIGURE_ROOT}/ROUND37_2026-09-14/figures/_common"
    if _c not in _s.path: _s.path.insert(0, _c)
    from v14lib import sidecar_v8
    info = sidecar_v8(hydrated(path))
    d = json.load(open(hydrated(path + ".provenance.json")))
    assert d.get("step", {}).get("rc", 0) == 0, (path, d.get("step"))
    ins = d.get("inputs", {}); paths = [v["path"] if isinstance(v, dict) else str(v) for v in (ins.values() if isinstance(ins, dict) else ins)]
    assert not any(("_SUPERSEDED" in p or "_PRE_V7" in p or "_PRE_V8" in p or "data_frozen_v7_2026-09" in p) for p in paths), paths
    return dict(sidecar=path + ".provenance.json", output_sha256=info["sha256"], stamped=info["output_mtime"], step=info["step"],
                script=d.get("script", {}).get("path"), status_ok=info["status_ok"])


def setup_matplotlib():
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import font_manager as fm
    for f in (ARIAL, ARIALB, ARIALI):
        fm.fontManager.addfont(f)
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "Arial", "font.sans-serif": ["Arial"],
        "axes.unicode_minus": True, "savefig.bbox": "standard", "savefig.pad_inches": 0,
        "axes.linewidth": 0.8, "axes.edgecolor": INK, "text.color": INK, "axes.labelcolor": INK,
        "xtick.color": INK, "ytick.color": INK, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
        "xtick.major.size": 3.0, "ytick.major.size": 3.0, "xtick.direction": "out", "ytick.direction": "out",
        "lines.solid_capstyle": "projecting", "figure.dpi": 72, "savefig.dpi": 72,
    })
    return plt


class Page:
    """Page-pinned matplotlib figure: coordinates in PDF points, top-down y, like the sheet's text layer."""
    def __init__(self, plt, W, H):
        self.W, self.H = W, H
        self.fig = plt.figure(figsize=(W / 72.0, H / 72.0))
        self.fig.patch.set_alpha(0.0)
        self.texts = []

    def axes(self, x0, y0, x1, y1):
        ax = self.fig.add_axes([x0 / self.W, 1 - y1 / self.H, (x1 - x0) / self.W, (y1 - y0) / self.H])
        ax.patch.set_alpha(0.0)
        return ax

    PT_PLUS = 4.0   # round 54 (Dragana, 2026-09-28): ticks, labels, keys, callouts 14 pt, axis titles 15, notes 13 (base sizes 10, 11, 9); round 49 was 1.0
    def text(self, x, baseline, s, size, ha="left", color=INK, weight="normal", rotation=0, linespacing=1.2, record=None):
        size = size + self.PT_PLUS
        t = self.fig.text(x / self.W, 1 - baseline / self.H, s, fontsize=size, ha=ha, va="baseline", color=color,
                          fontweight=weight, rotation=rotation, rotation_mode="anchor", linespacing=linespacing)
        self.texts.append(dict(text=s, x=x, baseline=baseline, size=size, ha=ha, color=color, rotation=rotation,
                               width=round(self.width(s, size - self.PT_PLUS, weight), 2), record=record))   # round 54: width kept for the overlap check
        return t

    def width(self, s, size, weight="normal"):
        size = size + self.PT_PLUS
        self.fig.canvas.draw()
        r = self.fig.canvas.get_renderer()
        t = self.fig.text(0.5, 0.5, s, fontsize=size, fontweight=weight)
        w = t.get_window_extent(renderer=r).width * 72.0 / self.fig.dpi
        t.remove()
        return w


def spans(page):
    out = []
    for b in page.get_text("rawdict")["blocks"]:
        if b.get("type", 0) != 0:
            continue
        for l in b.get("lines", []):
            for s in l["spans"]:
                txt = "".join(ch["c"] for ch in s["chars"])
                if not txt.strip():
                    continue
                x0 = min(ch["bbox"][0] for ch in s["chars"]); y0 = min(ch["bbox"][1] for ch in s["chars"])
                x1 = max(ch["bbox"][2] for ch in s["chars"]); y1 = max(ch["bbox"][3] for ch in s["chars"])
                out.append(dict(text=txt, font=s["font"], size=round(s["size"], 2),
                                origin=[round(s["chars"][0]["origin"][0], 2), round(s["chars"][0]["origin"][1], 2)],
                                bbox=[round(x0, 2), round(y0, 2), round(x1, 2), round(y1, 2)], color=s.get("color"),
                                ascender=s.get("ascender"), descender=s.get("descender")))
    return out


def words(page, clip=None):
    ws = page.get_text("words", clip=clip) if clip is not None else page.get_text("words")
    return [dict(text=w[4], bbox=[round(v, 2) for v in w[:4]]) for w in ws]


NUM_RE = re.compile(r"[+−\-]?\d[\d,]*\.?\d*")


def num_tokens(ws):
    toks = []
    for w in ws:
        for t in NUM_RE.findall(w["text"]):
            toks.append(t.replace("−", "-"))
    return Counter(toks)


def word_tokens(ws):
    return Counter(w["text"] for w in ws)


def letters(page):
    out = []
    for s in spans(page):
        if len(s["text"].strip()) == 1 and s["text"].strip() in "abcdefgh" and s["size"] >= 12.5 and "Bold" in s["font"]:
            out.append((s["text"].strip(), s["bbox"][0], s["bbox"][1], s["bbox"][3], s["font"], s["size"]))
    return sorted(out, key=lambda t: (t[2], t[1]))


def sheet_font_metrics(pdf_path):
    doc = fitz.open(hydrated(pdf_path)); met = {}
    for s in spans(doc[0]):
        if s["ascender"] is not None:
            met.setdefault(s["font"], (int(round(s["ascender"] * 1000)), int(round(s["descender"] * 1000))))
    doc.close(); gc.collect()
    return met


def align_font_metrics(src, dst, met):
    """Round-28 idiom: give the matplotlib Type 42 fonts the Ascent and Descent the sheet's own Arial descriptors carry,
    so the text layer reports the same span boxes. Glyphs, widths and positions untouched."""
    doc = fitz.open(src); page = doc[0]; done = {}
    for f in page.get_fonts(full=True):
        base = f[3].split("+")[-1]
        key = base if base in met else None
        if key is None:
            # matplotlib Type 42 names: ArialMT, Arial-BoldMT, Arial-ItalicMT
            alt = {"ArialMT": "ArialMT", "Arial-BoldMT": "Arial-BoldMT"}.get(base)
            key = alt if alt in met else None
        if key is None:
            continue
        asc, dsc = met[key]
        desc = doc.xref_get_key(f[0], "DescendantFonts")
        if desc[0] != "array":
            fdx_s = doc.xref_get_key(f[0], "FontDescriptor")
            if fdx_s[0] != "xref":
                continue
            fdx = int(fdx_s[1].split()[0])
        else:
            cid = int(re.search(r"(\d+) 0 R", desc[1]).group(1))
            fdx = int(doc.xref_get_key(cid, "FontDescriptor")[1].split()[0])
        doc.xref_set_key(fdx, "Ascent", str(asc)); doc.xref_set_key(fdx, "Descent", str(dsc)); done[base] = (asc, dsc)
    doc.save(dst, garbage=1, deflate=True); doc.close()
    return done


def stamp_letters(pg, items, size=13.0, asc=11.75):
    """Letters with fitz.TextWriter + Arial Bold (the composer idiom): items = [(ch, x, bbox_top)]."""
    font = fitz.Font(fontfile=ARIALB); tw = fitz.TextWriter(pg.rect)
    for ch, x, top in items:
        tw.append((x, top + asc), ch, font=font, fontsize=size)
    tw.write_text(pg, color=(0, 0, 0))


def render_png(pdf_path, png_path, dpi=150, clip=None, lock=False):
    """Render with a shared render lock only when asked (300 dpi or nested sheets)."""
    fh = None
    if lock:
        fh = open(RENDER_LOCK, "w"); fcntl.flock(fh, fcntl.LOCK_EX)
    try:
        doc = fitz.open(hydrated(pdf_path)); pg = doc[0]
        pm = pg.get_pixmap(dpi=dpi, colorspace=fitz.csRGB, alpha=False, clip=clip)
        pm.save(png_path); w, h = pm.width, pm.height
        pm = None; doc.close(); gc.collect()
    finally:
        if fh is not None:
            fcntl.flock(fh, fcntl.LOCK_UN); fh.close()
    return w, h


def raster_diff(png_a, png_b):
    """Return (mask of differing pixels as numpy bool array, shape). Both PNGs at the same dpi."""
    import numpy as np
    from PIL import Image
    a = np.asarray(Image.open(png_a).convert("RGB")).astype(int)
    b = np.asarray(Image.open(png_b).convert("RGB")).astype(int)
    h = min(a.shape[0], b.shape[0]); w = min(a.shape[1], b.shape[1])
    d = np.abs(a[:h, :w] - b[:h, :w]).sum(axis=2) > 30
    return d, (a.shape, b.shape)


def full_render_cached(pdf_path, dpi, lock=True):
    """One full-page render per (pdf sha256, dpi), in a child process under the render lock, cached under the lane's work/_render_cache
    (round 37: the nested V13 Main_Fig1 takes minutes per render, so it is rendered ONCE per dpi and cropped with PIL, never six times with a clip)."""
    import subprocess
    cache_dir = f"{LANE}/_render_cache"; os.makedirs(cache_dir, exist_ok=True)
    png = f"{cache_dir}/{sha256(pdf_path)[:16]}_{dpi}.png"
    if os.path.exists(png) and os.stat(png).st_blocks > 0 and os.path.getsize(png) > 0: return png
    code = "import fitz,sys\nd=fitz.open(sys.argv[1]);pm=d[0].get_pixmap(dpi=int(sys.argv[3]),colorspace=fitz.csRGB,alpha=False);pm.save(sys.argv[2]);print(pm.width,pm.height)\n"
    fh = None
    if lock: fh = open(RENDER_LOCK, "w"); fcntl.flock(fh, fcntl.LOCK_EX)
    try:
        tmp = png[:-4] + "_tmp.png"          # MuPDF picks the format from the extension
        r = subprocess.run([sys.executable, "-c", code, hydrated(pdf_path), tmp, str(dpi)], capture_output=True, text=True, timeout=3600)
        assert r.returncode == 0, r.stderr[-400:]
        os.replace(tmp, png)
    finally:
        if fh is not None: fcntl.flock(fh, fcntl.LOCK_UN); fh.close()
    return png


def crop_pair(old_pdf, new_pdf, clip, out_png, dpi=200, label_old="OLD (V13 base)", label_new="NEW (v8.1)"):
    """Side by side OLD vs NEW crop at 200 dpi, cut from one cached full-page render per sheet (PIL), never a clipped re-render."""
    from PIL import Image, ImageDraw
    ims = []
    for p in (old_pdf, new_pdf):
        full = Image.open(full_render_cached(p, dpi)); s = dpi / 72.0
        x0, y0, x1, y1 = clip; box = (max(int(x0 * s), 0), max(int(y0 * s), 0), min(int(x1 * s + 0.5), full.width), min(int(y1 * s + 0.5), full.height))
        ims.append(full.crop(box).convert("RGB"))
    h = max(i.height for i in ims) + 26; w = sum(i.width for i in ims) + 30
    canvas = Image.new("RGB", (w, h), (255, 255, 255)); dr = ImageDraw.Draw(canvas)
    x = 10
    for im, lab in zip(ims, (label_old, label_new)):
        dr.text((x, 4), lab, fill=(200, 0, 0)); canvas.paste(im, (x, 22)); x += im.width + 10
    canvas.save(out_png)
    return out_png


def write_changes(path, rows):
    import csv
    cols = ["sheet", "panel", "item", "old_value", "new_value", "source_file", "key_path", "dramatic", "note"]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in cols})


def banned_words(strings):
    bad = []
    for s in strings:
        low = s.lower()
        for b in BANNED:
            if b.lower() in low:
                if b in ("null",) and "null" in low and not re.search(r"\bnull\b", low):
                    continue
                bad.append((s, b))
    return bad


def checks_writer(path):
    lines = []
    def check(ok, name, detail=""):
        lines.append(f"{'PASS' if ok else 'FAIL'} {name}  {detail}")
        print(lines[-1][:400])
        return ok
    def finish():
        ok = all(l.startswith("PASS") for l in lines)
        lines.append("RESULT ALL PASS" if ok else "RESULT FAIL: " + "; ".join(l.split("  ")[0][5:] for l in lines if l.startswith("FAIL")))
        open(path, "w").write("\n".join(lines) + "\n")
        print(lines[-1])
        return ok
    return check, finish
