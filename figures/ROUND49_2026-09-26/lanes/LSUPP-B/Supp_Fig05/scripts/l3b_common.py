"""Lane V14_L3b_DURATION_SUPP, shared layer (round 37, 2026-09-15). Repointed copy of the round-30 l3b_common.py
(byte copy beside it as l3b_common_PRE_V8_1.py): BASE = the V13 set, the numbers gate = v14lib.sidecar_v8 (v8.1),
the OLD side of every delta = the v7 snapshot under V8_RECALC/compare/old_numbers (read only for old-vs-new,
never as a live source). The matplotlib design system, the composer and the text-layer probe are unchanged."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import hashlib, json, os, re, subprocess, sys, time
from collections import Counter

FIG = f"{paths.FIGURE_ROOT}/ROUND37_2026-09-14/figures"
sys.path.insert(0, f"{FIG}/_common")
import v14lib  # noqa: E402
from v14lib import hydrated, sha256, sidecar_v8, load_json as jload, RenderLock  # noqa: E402,F401

LANE = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LSUPP-B/Supp_Fig05"   # R49: this lane's sheet folder
BASE = v14lib.V13                    # the V13 design baseline (never edited)
T90 = v14lib.T90
NUM = v14lib.NUM
SV = v14lib.SV
OLDNUM = f"{v14lib.OLD_V7}/numbers"                                                   # OLD side only (v7 snapshot)
OLDEXT = paths.OLD_SNAPSHOT_SV  # OLD side only (v7 snapshot)
V8TAB = f"{paths.TABLES_DIR}"
FONTS = {"regular": paths.FONT_ARIAL,
         "bold": paths.FONT_ARIAL_BOLD,
         "italic": paths.FONT_ARIAL_ITALIC}

INK = "#1a1d21"
BLUE = "#0288d1"
BLUE_MID = "#3f9fd8"
GREY = "#8a9099"
GREY_PALE = "#ccd1d6"
ORANGE = "#d55e00"
ORANGE_PALE = "#eeb08f"
GREEN = "#298d32"
GREEN_MID = "#39c445"
GREEN_PALE = "#d5f1d0"
GREEN_BAND = "#eefaee"


def sidecar(path, require_date=None):
    """The v8.1 provenance gate (v14lib.sidecar_v8): sidecar present, inputs name data_frozen_v8_2026-09 (or a v8.1 STATUS ok row),
    the file's sha256 equals the sidecar's. require_date (YYYY-MM-DD) additionally pins the output mtime's date. Returns (sha256, sidecar)."""
    rec = sidecar_v8(path)
    d = jload(path + ".provenance.json")
    if require_date:
        m = d["output"]["mtime_local"]
        assert m.startswith(require_date), f"sidecar of {path} is dated {m}, need {require_date}"
    return rec["sha256"], d


# ------------------------------------------------------------------ matplotlib
import matplotlib  # noqa: E402
matplotlib.use("Agg")
from matplotlib import font_manager as fm  # noqa: E402
for _p in FONTS.values():
    fm.fontManager.addfont(_p)
import matplotlib.pyplot as plt  # noqa: E402

PT_PLUS = 1.0   # round 49 (Alen, 2026-09-26): every text one point larger; the geometry constants below are unchanged
LABF, TCKF, ANNF, FLOOR = 11.0 + PT_PLUS, 10.0 + PT_PLUS, 9.5 + PT_PLUS, 9.0 + PT_PLUS
CI_LW, MEDGE = 1.7, 0.8
S_CIRCLE, S_SQUARE = 28.0, 23.0
MS_CIRCLE, MS_SQUARE = 6.0, 5.4


def rc_polish():
    return {
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "standard",
        "figure.facecolor": "#ffffff", "axes.facecolor": "#ffffff", "savefig.facecolor": "#ffffff",
        "axes.spines.top": False, "axes.spines.right": False,
        "xtick.direction": "out", "ytick.direction": "out",
        "legend.frameon": False, "grid.alpha": 0.25, "grid.color": "#000000",
        "font.size": TCKF, "axes.labelsize": LABF, "axes.titlesize": LABF,
        "xtick.labelsize": TCKF, "ytick.labelsize": TCKF, "legend.fontsize": TCKF,
        "axes.linewidth": 0.8, "axes.edgecolor": INK, "text.color": INK,
        "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
        "xtick.major.width": 0.8, "ytick.major.width": 0.8,
        "xtick.major.size": 3.0, "ytick.major.size": 3.0,
    }


def measure(fig, s, fontsize, weight="normal", style="normal"):
    r = fig.canvas.get_renderer()
    t = fig.text(0.5, 0.5, s, fontsize=fontsize, fontweight=weight, fontstyle=style)
    w = t.get_window_extent(renderer=r).width / fig.dpi
    t.remove()
    return w


def wrap_measured(fig, s, max_in, fontsize):
    s = s.replace("\n", " ")
    if measure(fig, s, fontsize) <= max_in:
        return s
    words, lines, cur = s.split(), [], ""
    for wd in words:
        trial = f"{cur} {wd}".strip()
        if cur and measure(fig, trial, fontsize) > max_in:
            lines.append(cur)
            cur = wd
        else:
            cur = trial
    lines.append(cur)
    return "\n".join(lines)


def text_gate(fig):
    import matplotlib.text as mtext
    for t in fig.findobj(mtext.Text):
        s = t.get_text()
        if not s.strip() or not t.get_visible():
            continue
        assert "—" not in s and ";" not in s, f"banned punctuation: {s!r}"
        assert float(t.get_fontsize()) >= FLOOR, f"{t.get_fontsize()} pt under the floor: {s!r}"


# ------------------------------------------------------------------ fitz
import fitz  # noqa: E402


def compose(page_w, page_h, parts, letters, out_pdf, letter_size=13.0 + PT_PLUS):   # R49: letters 14 pt
    """Place part PDFs (path, x, y of their top-left in page points, scale 1) on a fresh page of the base sheet's
    width (the height may differ from the base when a panel grew) and stamp the letters with TextWriter + Arial Bold."""
    doc = fitz.open()
    page = doc.new_page(width=page_w, height=page_h)
    for path, ox, oy in parts:
        src = fitz.open(path)
        r = src[0].rect
        page.show_pdf_page(fitz.Rect(ox, oy, ox + r.width, oy + r.height), src, 0)
        src.close()
    if letters:
        tw = fitz.TextWriter(page.rect)
        font = fitz.Font(fontfile=FONTS["bold"])
        for glyph, x, y in letters:
            tw.append(fitz.Point(x, y), glyph, font=font, fontsize=letter_size)
        tw.write_text(page, color=(0x1a / 255, 0x1d / 255, 0x21 / 255))
    doc.save(out_pdf, garbage=3, deflate=True)
    doc.close()
    return out_pdf


def spans(pdf_path):
    """Text layer as one record per span: text, origin x and baseline y, x1, font, size, colour."""
    doc = fitz.open(hydrated(pdf_path))
    p = doc[0]
    out = []
    for b in p.get_text("rawdict")["blocks"]:
        for l in b.get("lines", []):
            for s in l["spans"]:
                chars = [c for c in s["chars"]]
                txt = "".join(c["c"] for c in chars).replace("\xa0", " ")
                if not txt.strip():
                    continue
                i0 = next(i for i, c in enumerate(chars) if c["c"].strip())
                i1 = max(i for i, c in enumerate(chars) if c["c"].strip())
                out.append(dict(text=txt.strip(), x=round(chars[i0]["origin"][0], 2),
                                y=round(chars[i0]["origin"][1], 2), x1=round(chars[i1]["bbox"][2], 2),
                                font=s["font"].split("+")[-1], size=round(s["size"], 2),
                                color="#%06x" % s["color"]))
    rect = (p.rect.width, p.rect.height)
    doc.close()
    out.sort(key=lambda s: (s["y"], s["x"]))
    return out, rect


def collapse_spans(sp, tol=0.25):
    """Stacked duplicate text copies (same string, origin within tol) count once, the round-28 rule."""
    out = []
    for s in sp:
        if not any(o["text"] == s["text"] and abs(o["x"] - s["x"]) <= tol and abs(o["y"] - s["y"]) <= tol for o in out):
            out.append(s)
    return out


def drawings(pdf_path):
    """get_drawings() is safe on these flat supplementary sheets (never on the nested main figures)."""
    doc = fitz.open(hydrated(pdf_path))
    dr = doc[0].get_drawings()
    items = []
    for g in dr:
        f, c, r = g.get("fill"), g.get("color"), g["rect"]
        fh = "#%02x%02x%02x" % tuple(int(round(v * 255)) for v in f) if f else None
        ch = "#%02x%02x%02x" % tuple(int(round(v * 255)) for v in c) if c else None
        kinds = "".join(sorted(set(it[0] for it in g["items"])))
        items.append(dict(rect=[round(v, 2) for v in (r.x0, r.y0, r.x1, r.y1)], fill=fh, stroke=ch,
                          width=g.get("width"), kinds=kinds, n=len(g["items"]), dashes=g.get("dashes")))
    doc.close()
    return items


def render(pdf_path, png_path, dpi=150, clip=None):
    """R49: Ghostscript png16m (r49_common.gs_render), never a PyMuPDF pixmap; whole pages only."""
    assert clip is None, "R49 renders whole pages with Ghostscript"
    sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LSUPP-B/scripts")
    import r49_common as _C
    _C.gs_render(pdf_path, png_path, dpi)
    return png_path


NUMRE = re.compile(r"-?\d[\d,]*\.?\d*")


def num_tokens(sp):
    return Counter(t for s in sp for t in NUMRE.findall(s["text"]))


def word_tokens(sp):
    return Counter(w for s in sp for w in s["text"].split())


def fmt_delta(c_old, c_new):
    rem, add = c_old - c_new, c_new - c_old
    return dict(removed=dict(rem), added=dict(add))


def drawn_record(panel, text, x, baseline, source_file, source_key, source_value, rule, ha="left", size=None):
    """One v14lib DRAWN_RECORD row."""
    return dict(panel=panel, text=text, x=round(float(x), 3), baseline=round(float(baseline), 3), ha=ha, size=size,
                source_file=source_file, source_key=source_key, source_value=source_value, rule=rule)
