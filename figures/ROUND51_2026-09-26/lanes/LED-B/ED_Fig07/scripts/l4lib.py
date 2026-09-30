#!/usr/bin/env python3
"""Round 49 LED-B copy (2026-09-26, +1 pt rule) of the round-38 copy. V14 lane V14_L4_APNEA copy (round 37, 2026-09-15) of the lane V7_L4_APNEA shared helpers (round 30, 2026-09-08); repointed: LANE, BASE = V13, sidecar = v8.1 gate. Original beside as l4lib_PRE_V8_1.py.: eviction gate, text layer, multisets, letter stamping, font
metric alignment, checks writer, CHANGES csv, provenance, half-up rounding, matplotlib rc. Import from any sheet script:
    import sys; sys.path.insert(0, "$T90_FIGURE_ROOT/ROUND30_2026-09-08/V7_L4_APNEA/scripts"); import l4lib as L
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import collections, csv, datetime, fcntl, gc, hashlib, json, os, re, subprocess, time
from decimal import Decimal, ROUND_HALF_UP
import fitz

T90 = paths.FIGURE_ROOT
R30 = f"{paths.FIGURE_ROOT}/ROUND30_2026-09-08"
LANE = f"{paths.FIGURE_ROOT}/ROUND51_2026-09-26/lanes/LED-B"   # round 51 (lane LED-B): this lane root (SD = LANE/<Sheet> = the sheet folder)
PT_PLUS = 1.0   # round 49 (Alen, 2026-09-26): every text one point larger
BASE = f"{paths.FIGURE_ROOT}/ROUND36_2026-09-10/L9_assembly/NEW_FINAL_SET_V13"   # V14: the V13 design baseline
NUM = f"{paths.NUMBERS_DIR}"
SV = paths.SV_ROOT
ARIAL = paths.FONT_ARIAL
ARIAL_BOLD = paths.FONT_ARIAL_BOLD
ARIAL_ITALIC = paths.FONT_ARIAL_ITALIC
INK = "#1a1d21"; BLUE = "#0288d1"; ORANGE = "#d55e00"
ORANGE_LADDER = ["#f6d6c2", "#eeb08f", "#e4864f", "#d55e00"]
BLUE_LADDER = ["#b3dcf2", "#7cc0e9", "#3f9fd8", "#0288d1"]
RENDER_LOCK = f"{R30}/_render.lock"


# ------------------------------------------------------------------ eviction gate
def blocks(path): return int(subprocess.run(["stat", "-f", "%b", path], capture_output=True, text=True).stdout.strip() or 0)


def hydrated(path, wait_s=900):
    """Return path once it is local. 0 blocks with a nonzero size = iCloud-evicted: brctl download and poll."""
    if not os.path.exists(path): raise FileNotFoundError(path)
    if os.path.isdir(path): return path
    if os.path.getsize(path) == 0 or blocks(path) > 0: return path
    subprocess.run(["brctl", "download", path], check=False)
    t0 = time.time()
    while blocks(path) == 0 and os.path.getsize(path) > 0:
        if time.time() - t0 > wait_s: raise TimeoutError(f"still evicted after {wait_s}s: {path}")
        time.sleep(5)
    return path


def sha256(path, n=None):
    h = hashlib.sha256(); h.update(open(hydrated(path), "rb").read()); return h.hexdigest() if n is None else h.hexdigest()[:n]


def sidecar(path):
    """V14 (round 37): the provenance sidecar beside a numbers file, asserted to match the file's sha256 and to prove a v8.1
    output (inputs name data_frozen_v8_2026-09, or a STATUS ok row at or after the v8.1 table build): v14lib.sidecar_v8."""
    import sys as _s
    _c = f"{paths.FIGURE_ROOT}/ROUND37_2026-09-14/figures/_common"
    if _c not in _s.path: _s.path.insert(0, _c)
    from v14lib import sidecar_v8
    gate = sidecar_v8(path)
    sc = json.load(open(hydrated(path + ".provenance.json"))); sc["_v14_gate"] = gate
    return sc


# ------------------------------------------------------------------ text layer
def spans_of(page_or_path):
    close = False
    if isinstance(page_or_path, str):
        doc = fitz.open(hydrated(page_or_path)); page = doc[0]; close = True
    else: page = page_or_path
    out = []
    for b in page.get_text("rawdict")["blocks"]:
        if b["type"] != 0: continue
        for l in b["lines"]:
            for s in l["spans"]:
                txt = "".join(ch["c"] for ch in s["chars"]); org = s["chars"][0]["origin"]
                c = s.get("color"); col = ("#%06x" % c) if isinstance(c, int) else str(c)
                out.append(dict(text=txt, font=s["font"], size=round(s["size"], 3), color=col, bbox=[round(v, 3) for v in s["bbox"]],
                                origin=[round(org[0], 3), round(org[1], 3)], asc=round(s.get("ascender", 0), 4), dsc=round(s.get("descender", 0), 4)))
    if close: doc.close(); gc.collect()
    return out


def words_of(page_or_path):
    close = False
    if isinstance(page_or_path, str):
        doc = fitz.open(hydrated(page_or_path)); page = doc[0]; close = True
    else: page = page_or_path
    out = [dict(text=w[4], bbox=[round(v, 3) for v in w[:4]]) for w in page.get_text("words")]
    if close: doc.close(); gc.collect()
    return out


def dedupe(items, tol=0.25):
    """Collapse stacked copies: same text within tol pt of an already kept item (round-28 rule)."""
    kept = collections.defaultdict(list); out = []
    for it in items:
        key = it["text"].replace("\xa0", " "); x, y = it["bbox"][0], it["bbox"][1]
        if any(abs(x - kx) <= tol and abs(y - ky) <= tol for kx, ky in kept[key]): continue
        kept[key].append((x, y)); out.append(it)
    return out


NUMRE = re.compile(r"\d[\d,]*\.?\d*")


def num_tokens(words):
    t = []
    for w in words: t += NUMRE.findall(w["text"].replace("\xa0", " ") if isinstance(w, dict) else w)
    return sorted(t)


def word_tokens(words):
    return sorted(w["text"].replace("\xa0", " ") if isinstance(w, dict) else w for w in words)


def delta(a, b):
    ca, cb = collections.Counter(a), collections.Counter(b)
    return dict(lost=dict(ca - cb), gained=dict(cb - ca))


def letters(spans):
    return sorted([(s["text"], round(s["bbox"][0], 2), round(s["bbox"][1], 2), round(s["origin"][1], 2), s["font"], round(s["size"], 1))
                   for s in spans if len(s["text"].strip()) == 1 and s["text"].strip() in "abcdefgh" and s["size"] >= 12], key=lambda t: (t[2], t[1]))


# ------------------------------------------------------------------ stamping and fonts
def stamp_letter(page, ch, x, baseline, size=13.0, color=(0, 0, 0)):
    """Panel letter with fitz.TextWriter + Arial Bold (never insert_text with a fontfile)."""
    tw = fitz.TextWriter(page.rect); tw.append(fitz.Point(x, baseline), ch, font=fitz.Font(fontfile=ARIAL_BOLD), fontsize=size)
    tw.write_text(page, color=color)


def sheet_font_metrics(spans):
    """Ascent and descent (per mille) the base sheet's text layer reports for each face."""
    met = {}
    for s in spans: met.setdefault(s["font"], (int(round(s["asc"] * 1000)), int(round(s["dsc"] * 1000))))
    return met


def align_font_metrics(src, dst, met):
    """Give the matplotlib fonts the FontDescriptor Ascent/Descent the same face carries on the base sheet (round 28)."""
    doc = fitz.open(src); page = doc[0]; done = {}
    for f in page.get_fonts(full=True):
        base = f[3].split("+")[-1]
        if base not in met: continue
        asc, dsc = met[base]
        desc = doc.xref_get_key(f[0], "DescendantFonts")
        if desc[0] != "array": continue
        cid = int(re.search(r"(\d+) 0 R", desc[1]).group(1))
        fdx = int(doc.xref_get_key(cid, "FontDescriptor")[1].split()[0])
        doc.xref_set_key(fdx, "Ascent", str(asc)); doc.xref_set_key(fdx, "Descent", str(dsc)); done[base] = (asc, dsc)
    doc.save(dst, garbage=1, deflate=True); doc.close()
    return done


def mpl_setup():
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import font_manager as fm
    for f in (ARIAL, ARIAL_BOLD, ARIAL_ITALIC): fm.fontManager.addfont(f)
    import matplotlib.pyplot as plt
    return matplotlib, plt


RC = {"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42, "text.color": INK, "axes.edgecolor": INK, "axes.labelcolor": INK,
      "xtick.color": INK, "ytick.color": INK, "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
      "xtick.major.size": 3.0, "ytick.major.size": 3.0, "xtick.direction": "out", "ytick.direction": "out", "axes.unicode_minus": False,
      "savefig.bbox": "standard", "savefig.pad_inches": 0.0, "figure.dpi": 300, "savefig.dpi": 300, "mathtext.default": "regular"}


# ------------------------------------------------------------------ rendering (one heavy render at a time across lanes)
class render_lock:
    def __enter__(self): self.f = open(RENDER_LOCK, "w"); fcntl.flock(self.f, fcntl.LOCK_EX); return self
    def __exit__(self, *a): fcntl.flock(self.f, fcntl.LOCK_UN); self.f.close()


def render(pdf, png, dpi=150, clip=None, heavy=False):
    def go():
        d = fitz.open(hydrated(pdf)); p = d[0]
        if clip is None: pm = p.get_pixmap(dpi=dpi, colorspace=fitz.csRGB, alpha=False)
        else:
            Z = dpi / 72; x0, y0, x1, y1 = [int(round(float(v) * Z)) for v in clip]
            pm = p.get_pixmap(dpi=dpi, clip=fitz.Rect(x0 / Z, y0 / Z, x1 / Z, y1 / Z), colorspace=fitz.csRGB, alpha=False)
        pm.save(png); w, h = pm.w, pm.h; pm = None; d.close(); gc.collect(); return w, h
    t0 = time.time()
    if heavy:
        with render_lock(): w, h = go()
    else: w, h = go()
    return dict(png=png, w=w, h=h, dpi=dpi, seconds=round(time.time() - t0, 1))


# ------------------------------------------------------------------ rounding and strings
def r2(v): return str(Decimal(repr(float(v))).quantize(Decimal("0.01"), rounding=ROUND_HALF_UP))
def r1(v): return str(Decimal(repr(float(v))).quantize(Decimal("0.1"), rounding=ROUND_HALF_UP))
def r3(v): return str(Decimal(repr(float(v))).quantize(Decimal("0.001"), rounding=ROUND_HALF_UP))
def hr_ci(hr, lo, hi): return f"{r2(hr)} ({r2(lo)}-{r2(hi)})"
def q_text(q): return "<0.001" if q < 0.001 else r3(q)
def stars(q): return "" if q is None else ("***" if q < 0.001 else "**" if q < 0.01 else "*" if q < 0.05 else "")
def pct1(num, den): return str((Decimal(num) * 100 / Decimal(den)).quantize(Decimal("0.1"), rounding=ROUND_HALF_UP))
BANNED = ("—", ";", "strata", "stratum", "honest", "straightforward", "null", "pre-specified", "prespecified")


def house_ok(s):
    low = s.lower(); return not any(b in (s if b in ("—", ";") else low) for b in BANNED)


# ------------------------------------------------------------------ checks and CHANGES
class Checks:
    def __init__(self, title):
        self.lines = [title, ""]; self.ok = True
    def check(self, name, ok, detail=""):
        ok = bool(ok); self.ok &= ok
        self.lines.append(f"{'PASS' if ok else 'FAIL'}  {name}" + (f": {detail}" if detail else "")); print(self.lines[-1], flush=True)
    def info(self, text): self.lines.append(f"INFO  {text}"); print(self.lines[-1], flush=True)
    def write(self, path):
        self.lines += ["", "RESULT " + ("ALL PASS" if self.ok else "FAILED: " + "; ".join(l[6:].split(":")[0] for l in self.lines if l.startswith("FAIL")))]
        open(path, "w").write("\n".join(self.lines) + "\n"); print(self.lines[-1]); return self.ok


CHANGES_HEAD = ["sheet", "panel", "label", "old_value", "new_value", "source_file", "key", "dramatic", "note"]


def write_changes(path, rows):
    with open(path, "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=CHANGES_HEAD, extrasaction="ignore"); wr.writeheader()
        for r in rows: wr.writerow({k: r.get(k, "") for k in CHANGES_HEAD})
    return len(rows)


def provenance(out_path, inputs, script, extra=None):
    rec = {"output": {"path": out_path, "sha256": sha256(out_path) if os.path.exists(out_path) else None, "bytes": os.path.getsize(out_path) if os.path.exists(out_path) else None},
           "inputs": {k: {"path": v, "sha256": sha256(v)} for k, v in inputs.items()}, "script": {"path": script, "sha256": sha256(script)},
           "written": datetime.datetime.now().isoformat(timespec="seconds"), "lane": "V14_L4_APNEA"}
    if extra: rec.update(extra)
    json.dump(rec, open(out_path + ".provenance.json", "w"), indent=1); return rec


def now(): return datetime.datetime.now().strftime("%H:%M")


# ------------------------------------------------------------------ the hazard-ratio red ramp (round 27 lane B, Figs 3c, 4d, ED 7b)
# Old ramp (R18 j3_fig4d.py): log-even stops between HR 1 and 3.25. Round 27 mapped every old fill by its CIELAB lightness
# t = (L_palest - L) / (L_palest - L_darkest) with palest #fdefe9 and darkest #d15641, onto the straight CIELAB line from
# #fdecea (t = 0) to #d32f2f (t = 1). A fresh cell is coloured by the same two steps, so it lands exactly where the recolour
# would have put it. Text: paper ink where it reaches 4.5:1 on the fill, else white where white reaches it, else the higher.
import math as _math
import colour_math as CM
RAMP_OLD_STOPS = [(1.00, "#fdefe9"), (1.25, "#fbe3dc"), (1.90, "#f2a98f"), (3.25, "#d15540")]
RAMP_MIN, RAMP_MAX = 1.0, 3.25
RAMP_PALE, RAMP_TOP = "#fdecea", "#d32f2f"
_L_PALE = CM.rgb_to_lab(CM.hex_to_rgb("#fdefe9"))[0]; _L_DARK = CM.rgb_to_lab(CM.hex_to_rgb("#d15641"))[0]


def ramp_old(hr):
    x = min(max(_math.log(hr), _math.log(RAMP_MIN)), _math.log(RAMP_MAX))
    for (h0, c0), (h1, c1) in zip(RAMP_OLD_STOPS[:-1], RAMP_OLD_STOPS[1:]):
        if x <= _math.log(h1) + 1e-12:
            t = (x - _math.log(h0)) / (_math.log(h1) - _math.log(h0)); a, b = CM.hex_to_rgb(c0), CM.hex_to_rgb(c1)
            return tuple(u + t * (v - u) for u, v in zip(a, b))
    return CM.hex_to_rgb(RAMP_OLD_STOPS[-1][1])


def ramp_t(hr):
    L = CM.rgb_to_lab(ramp_old(hr))[0]; return min(1.0, max(0.0, (_L_PALE - L) / (_L_PALE - _L_DARK)))


def ramp_hex(hr):
    t = ramp_t(hr); a, b = CM.rgb_to_lab(CM.hex_to_rgb(RAMP_PALE)), CM.rgb_to_lab(CM.hex_to_rgb(RAMP_TOP))
    return CM.rgb_to_hex(CM.lab_to_rgb(tuple(u + t * (v - u) for u, v in zip(a, b))))


def text_on(fill_hex, ink="#1a1d21"):
    c_ink = CM.contrast(CM.hex_to_rgb(ink), CM.hex_to_rgb(fill_hex)); c_w = CM.contrast((1, 1, 1), CM.hex_to_rgb(fill_hex))
    return ink if c_ink >= 4.5 else ("#ffffff" if c_w >= 4.5 else (ink if c_ink >= c_w else "#ffffff"))


def hex_rgb01(h): return CM.hex_to_rgb(h)
