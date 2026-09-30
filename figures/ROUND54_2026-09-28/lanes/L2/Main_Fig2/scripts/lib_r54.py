#!$T90_PY
"""Round 54 (2026-09-29), lane L2, Main_Fig2: Ghostscript-only helpers (copy of the round-49 lib_r49.py, byte copy beside as
lib_r49_R49_ORIGINAL.py, with the round-54 additions at the end: regions of the L-shaped sheet, the class-based census against V30,
glyph-box overlaps, the pdfwrite flatten). Round-49 docstring: Ghostscript-only helpers for the composed sheet (renders and text layers under the
round-38 watchdog wd_run.sh with THIS lane's log dir, one command at a time, the round-root render lock held during every render),
the round-49 text census (every string present, every size +1.0, exceptions by string), the checks writer and the CHANGES csv.
Adapted from ROUND40 lanes/LMAIN/scripts/lmain_lib.py and gs_text.py (byte copies beside as *_R40_ORIGINAL.py). PyMuPDF is never
used on a composed sheet here."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import csv, fcntl, hashlib, html, json, os, re, subprocess, time, collections

LANE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WD = os.path.join(LANE, "scripts", "wd_run.sh"); WD_LOGDIR = os.path.join(LANE, "logs", "watchdog")
GS = paths.GS; PY = paths.PY
LOCK_ROOT = f"{paths.FIGURE_ROOT}/ROUND30_2026-09-08/_render.lock"
ARIALB = paths.FONT_ARIAL_BOLD
INK_RGB = (0x1a / 255, 0x1d / 255, 0x21 / 255)


def blocks(path):
    return int(subprocess.run(["stat", "-f", "%b", path], capture_output=True, text=True).stdout.strip() or 0)


def hydrated(path, wait_s=600):
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


def now(): return time.strftime("%Y-%m-%d %H:%M:%S")


class RenderLock:
    """ONE heavy render at a time across lanes (the round-root lock every lane uses)."""
    def __enter__(self):
        self.f = open(LOCK_ROOT, "a"); fcntl.flock(self.f, fcntl.LOCK_EX); return self
    def __exit__(self, *a):
        try: fcntl.flock(self.f, fcntl.LOCK_UN); self.f.close()
        except Exception: pass


def wd(name, cmd, max_s=900, rss_kb=4000000):
    os.makedirs(WD_LOGDIR, exist_ok=True)
    env = dict(os.environ, WD_LOGDIR=WD_LOGDIR)
    r = subprocess.run([WD, name, str(max_s), str(rss_kb), "--"] + list(cmd), env=env, capture_output=True, text=True)
    st = open(f"{WD_LOGDIR}/{name}.status").read().strip() if os.path.exists(f"{WD_LOGDIR}/{name}.status") else r.stdout.strip()
    return r.returncode, st


def wd_out(name): return open(f"{WD_LOGDIR}/{name}.out").read() if os.path.exists(f"{WD_LOGDIR}/{name}.out") else ""


def gs_render(pdf, png, dpi, name, max_s=900):
    """png16m render through Ghostscript under the watchdog and the round-root render lock; asserts rc 0, a fresh file, no Error line."""
    hydrated(pdf)
    cmd = [GS, "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-dFirstPage=1", "-dLastPage=1", "-sDEVICE=png16m", f"-r{dpi}", f"-sOutputFile={png}", pdf]
    last = None
    for attempt in (1, 2):
        if os.path.exists(png): os.remove(png)
        t0 = time.time()
        with RenderLock():
            rc, st = wd(name if attempt == 1 else f"{name}_retry", cmd, max_s=max_s)
        out = wd_out(name if attempt == 1 else f"{name}_retry")
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


def pypdf_mediabox(pdf):
    from pypdf import PdfReader
    r = PdfReader(hydrated(pdf)); mb = r.pages[0].mediabox
    return len(r.pages), round(float(mb.width), 3), round(float(mb.height), 3)


# ------------------------------------------------------------------ text layer (gs txtwrite -dTextFormat=0)
def norm_text(s):
    return s.replace("\xa0", " ").replace("−", "-").replace("–", "-").replace("—", "-").replace("\xad", "-")


def spans_xml(pdf, xml, name, max_s=600):
    hydrated(pdf)
    rc, st = wd(name, [GS, "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=txtwrite", "-dTextFormat=0", f"-sOutputFile={xml}", pdf], max_s=max_s)
    assert rc == 0 and os.path.exists(xml), ("gs txtwrite failed", name, st)
    return parse(xml)


def parse(xml_path):
    """Every <span> as dict(text, x0, y, x1, y0, y1, font, size). Characters joined in order (numeric entities unescaped)."""
    xml = open(xml_path, encoding="utf-8", errors="replace").read(); out = []
    for m in re.finditer(r'<span bbox="([^"]*)" font="([^"]*)" size="([^"]*)">(.*?)</span>', xml, re.S):
        bb = [float(v) for v in m.group(1).split()]; chars = re.findall(r'c="([^"]*)"', m.group(4))
        txt = norm_text(html.unescape("".join(chars)))
        if not txt.strip(): continue
        out.append(dict(text=txt.strip(), x0=bb[0], y=bb[1], x1=bb[2], y1=bb[3], font=m.group(2).split("+")[-1].replace("-Identity-H", ""), size=round(float(m.group(3)), 2)))
    return out


def strings_of(spans, split_x=479.0, dy_c=780.76, abut=1.5):
    """The sheet's strings rebuilt from gs spans: spans on one baseline (same region, same size) sorted by x, MERGED when they abut
    (gap below abut pt: Ghostscript splits one string into pieces such as '1.1' + '1 (1.06-1.16)' or 'T' + 'ype 2 diabetes'), kept
    separate otherwise. Returns dicts(text, x0, x1, y, size, font, region)."""
    def region(x, y): return "c" if y >= dy_c else ("a" if x < split_x else "b")
    rows = collections.defaultdict(list)
    for s in spans: rows[(region(s["x0"], s["y"]), round(s["y"]), s["size"])].append(s)
    out = []
    for (reg, _, sz), v in rows.items():
        cur = None
        for s in sorted(v, key=lambda s: s["x0"]):
            if cur is not None and s["x0"] - cur["x1"] < abut:
                cur["text"] += s["text"]; cur["x1"] = max(cur["x1"], s["x1"])
            else:
                if cur is not None: out.append(cur)
                cur = dict(text=s["text"], x0=s["x0"], x1=s["x1"], y=s["y"], size=sz, font=s["font"], region=reg)
        if cur is not None: out.append(cur)
    return out


def lines_of(spans):
    """Spans joined per text line (same baseline, same size, x order) so that a different span split between two sheets does not count."""
    rows = collections.defaultdict(list)
    for s in spans: rows[(round(s["y"]), s["size"])].append((s["x0"], s["text"]))
    return collections.Counter(" ".join(t for _, t in sorted(v)) for v in rows.values())


def census(new_spans, old_spans, plus=1.0, exceptions=None):
    """Round-49 census: same multiset of strings (span level, line level as the fallback when the span split differs), every (string,
    size) of the old sheet reappears as (string, size + plus), except the declared exceptions {string: old size kept}."""
    exc = exceptions or {}
    n_str = collections.Counter(s["text"] for s in new_spans); o_str = collections.Counter(s["text"] for s in old_spans)
    same_strings = n_str == o_str
    lost = sorted((o_str - n_str).elements()); gained = sorted((n_str - o_str).elements())
    ln = collections.Counter(s["text"] for s in strings_of(new_spans)); lo = collections.Counter(s["text"] for s in strings_of(old_spans)); line_same = ln == lo
    m_lost = sorted((lo - ln).elements()); m_gained = sorted((ln - lo).elements())
    o_pairs = collections.Counter((s["text"], s["size"]) for s in old_spans); n_pairs = collections.Counter((s["text"], s["size"]) for s in new_spans)
    expected = collections.Counter()
    for (t, sz), k in o_pairs.items(): expected[(t, sz if (t in exc and abs(exc[t] - sz) < 0.01) else round(sz + plus, 2))] += k
    size_ok = expected == n_pairs
    # line-level size check (robust to a different span split): every old line (text, size) -> (text, size + plus)
    lo_pairs = collections.Counter(); ln_pairs = collections.Counter()
    for s in strings_of(old_spans): lo_pairs[(s["text"], s["size"])] += 1
    for s in strings_of(new_spans): ln_pairs[(s["text"], s["size"])] += 1
    exp_l = collections.Counter()
    for (t, sz), k in lo_pairs.items(): exp_l[(t, sz if (t in exc and abs(exc[t] - sz) < 0.01) else round(sz + plus, 2))] += k
    line_size_ok = exp_l == ln_pairs
    miss = sorted((expected - n_pairs).elements()); extra = sorted((n_pairs - expected).elements())
    lmiss = sorted((exp_l - ln_pairs).elements()); lextra = sorted((ln_pairs - exp_l).elements())
    fonts_o = collections.Counter(s["font"] for s in old_spans); fonts_n = collections.Counter(s["font"] for s in new_spans)
    sizes_o = sorted({s["size"] for s in old_spans}); sizes_n = sorted({s["size"] for s in new_spans})
    # per-size accounting old -> new
    per_size = {str(sz): dict(old_spans=sum(1 for s in old_spans if s["size"] == sz), new_size=round(sz + plus, 2), new_spans=sum(1 for s in new_spans if abs(s["size"] - sz - plus) < 0.01)) for sz in sizes_o}
    return dict(n_spans_new=len(new_spans), n_spans_v26=len(old_spans), same_strings=same_strings, lost=lost, gained=gained, line_level_same=line_same, merged_lost=m_lost, merged_gained=m_gained,
                n_strings_new=sum(ln.values()), n_strings_v26=sum(lo.values()),
                sizes_v26=sizes_o, sizes_new=sizes_n, every_size_plus_one=size_ok, every_line_size_plus_one=line_size_ok, missing_expected=miss, unexpected=extra,
                line_missing_expected=lmiss, line_unexpected=lextra, exceptions=exc, fonts_v26=dict(fonts_o), fonts_new=dict(fonts_n), per_size=per_size, plus=plus)


def lines_of_sized(spans):
    rows = collections.defaultdict(list)
    for s in spans: rows[(round(s["y"]), s["size"])].append((s["x0"], s["text"]))
    return collections.Counter((" ".join(t for _, t in sorted(v)), k[1]) for k, v in rows.items())


def text_overlaps(spans, same_line_tol=1.5, cross_tol=0.5):
    """Pairs of spans whose bboxes overlap: on one baseline only when the x overlap exceeds same_line_tol (gs splits one string into
    abutting spans), on different baselines when both overlaps exceed cross_tol."""
    out = []
    for i in range(len(spans)):
        a = spans[i]
        for j in range(i + 1, len(spans)):
            b = spans[j]
            ox = min(a["x1"], b["x1"]) - max(a["x0"], b["x0"]); oy = min(a["y1"], b["y1"]) - max(a["y0"], b["y0"])
            dy = abs(a["y"] - b["y"]); x0 = min(a["x0"], b["x0"]); ymin = min(a["y"], b["y"])
            if dy < 0.5:
                if ox > same_line_tol and oy > 0: out.append((a["text"], b["text"], round(ox, 2), round(oy, 2), round(dy, 2), x0, ymin))
            elif ox > cross_tol and oy > cross_tol: out.append((a["text"], b["text"], round(ox, 2), round(oy, 2), round(dy, 2), x0, ymin))
    return out


def with_boxes(spans):
    """Add y0 (top) and y1 (bottom) to gs spans: txtwrite -dTextFormat=0 gives (x0, baseline, x1, baseline), no vertical extent, so the
    full Arial box is used: top = baseline - 0.905 * size (ascender), bottom = baseline + 0.212 * size (descender)."""
    for s in spans: s["y0"] = s["y"] - 0.905 * s["size"]; s["y1"] = s["y"] + 0.212 * s["size"]
    return spans


# ------------------------------------------------------------------ raster crops
def crop_pair(old_png, new_png, rect_pt, dpi, out_path, label_old="OLD (V26)", label_new="NEW (round 49)"):
    from PIL import Image, ImageDraw
    S = dpi / 72.0; x0, y0, x1, y1 = [int(round(v * S)) for v in rect_pt]; ims = []
    for p in (old_png, new_png):
        im = Image.open(p).convert("RGB"); ims.append(im.crop((max(x0, 0), max(y0, 0), min(x1, im.width), min(y1, im.height))))
    h = max(i.height for i in ims) + 26; w = sum(i.width for i in ims) + 30
    canvas = Image.new("RGB", (w, h), (255, 255, 255)); dr = ImageDraw.Draw(canvas); x = 10
    for im, lab in zip(ims, (label_old, label_new)):
        dr.text((x, 4), lab, fill=(200, 0, 0)); canvas.paste(im, (x, 22)); x += im.width + 10
    canvas.save(out_path); return out_path


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


# ================================================================== round 54 additions
def region_of(x, y, split_x, bc_y):
    """Panel of a point on the round-54 sheet: a left of split_x, b right of it above bc_y, c right of it below bc_y."""
    return "a" if x < split_x else ("b" if y < bc_y else "c")


def strings_r54(spans, split_x, bc_y, abut=1.5, region_fn=None):
    """As strings_of, with the region function of the sheet given (V30: region_of(x, y, 479, 780.76) reproduces round 49's split)."""
    rf = region_fn or (lambda x, y: region_of(x, y, split_x, bc_y))
    rows = collections.defaultdict(list)
    for s in spans: rows[(rf(s["x0"], s["y"]), round(s["y"]), s["size"])].append(s)
    out = []
    for (reg, _, sz), v in rows.items():
        cur = None
        for s in sorted(v, key=lambda s: s["x0"]):
            if cur is not None and s["x0"] - cur["x1"] < abut:
                cur["text"] += s["text"]; cur["x1"] = max(cur["x1"], s["x1"])
            else:
                if cur is not None: out.append(cur)
                cur = dict(text=s["text"], x0=s["x0"], x1=s["x1"], y=s["y"], size=sz, font=s["font"], region=reg)
        if cur is not None: out.append(cur)
    return out


def census_r54(new_strings, old_strings, rule, joins):
    """Round-54 census at the merged-string level. rule(old_string) -> the size the string must have on the new sheet (None = no rule,
    a failure). joins: list of (old texts in order, new text): the old strings are replaced by the joined string before comparing.
    Returns the dict with the pass flags, the missing and unexpected (text, size) pairs, the per-class accounting."""
    old = [dict(s) for s in old_strings]
    joined = []
    for parts, new_text in joins:
        found = []
        for part in parts:
            c = [s for s in old if s["text"] == part and s not in found]
            assert c, ("join part not on the old sheet", part)
            found.append(min(c, key=lambda s: s["y"]))
        for s in found: old.remove(s)
        j = dict(found[0]); j["text"] = new_text; j["joined_from"] = parts; old.append(j); joined.append((parts, new_text))
    expected = collections.Counter(); no_rule = []; per_class = collections.Counter()
    for s in old:
        sz = rule(s)
        if sz is None: no_rule.append((s["text"], s["size"], s["region"])); continue
        expected[(s["text"], round(sz, 2))] += 1; per_class[(s["region"], s["size"], round(sz, 2))] += 1
    got = collections.Counter((s["text"], s["size"]) for s in new_strings)
    missing = sorted((expected - got).elements()); unexpected = sorted((got - expected).elements())
    o_txt = collections.Counter(s["text"] for s in old); n_txt = collections.Counter(s["text"] for s in new_strings)
    return dict(n_old=len(old_strings), n_old_after_joins=len(old), n_new=len(new_strings), joins=joined, no_rule=no_rule,
                same_strings=o_txt == n_txt, lost=sorted((o_txt - n_txt).elements()), gained=sorted((n_txt - o_txt).elements()),
                sizes_ok=(not missing and not unexpected and not no_rule), missing=missing, unexpected=unexpected,
                per_class={f"{k[0]}:{k[1]}->{k[2]}": v for k, v in sorted(per_class.items())},
                sizes_old=sorted({s["size"] for s in old_strings}), sizes_new=sorted({s["size"] for s in new_strings}),
                min_size_new=min(s["size"] for s in new_strings), fonts_old=dict(collections.Counter(s["font"] for s in old_strings)), fonts_new=dict(collections.Counter(s["font"] for s in new_strings)))


CAP, DESC = 0.716, 0.212
def with_glyph_boxes(spans):
    """y0 (top) and y1 (bottom) of the glyph box on the gs baseline: cap height 0.716 and descender 0.212 of the size (Arial)."""
    for s in spans: s["y0"] = s["y"] - CAP * s["size"]; s["y1"] = s["y"] + DESC * s["size"]
    return spans


def glyph_overlaps(spans, same_line_tol=1.5, tol=0.3):
    """Pairs of strings whose glyph boxes overlap by more than tol in both directions (same-baseline pieces of one string excluded)."""
    out = []
    for i in range(len(spans)):
        a = spans[i]
        for j in range(i + 1, len(spans)):
            b = spans[j]
            ox = min(a["x1"], b["x1"]) - max(a["x0"], b["x0"]); oy = min(a["y1"], b["y1"]) - max(a["y0"], b["y0"])
            if abs(a["y"] - b["y"]) < 0.5 and abs(a["size"] - b["size"]) < 0.01:
                if ox > same_line_tol and oy > 0: out.append((a["text"], b["text"], round(ox, 2), round(oy, 2), round(min(a["x0"], b["x0"]), 1), round(min(a["y"], b["y"]), 1)))
            elif ox > tol and oy > tol: out.append((a["text"], b["text"], round(ox, 2), round(oy, 2), round(min(a["x0"], b["x0"]), 1), round(min(a["y"], b["y"]), 1)))
    return out


def gs_flatten(pdf, out, name, max_s=600):
    """Single-page pdfwrite re-distillation (the round-44 recipe used in rounds 49, 50 and 52): prepress, no downsampling, Flate images, colours unchanged, fonts embedded, no annotations."""
    hydrated(pdf)
    # the round-44 recipe as in ROUND49/50/52 lanes/LSUPP-D/scripts/lsd_common.py gs_flatten, flag for flag
    cmd = [GS, "-dSAFER", "-dBATCH", "-dNOPAUSE", "-dQUIET", "-sDEVICE=pdfwrite", "-dCompatibilityLevel=1.7", "-dPDFSETTINGS=/prepress",
           "-dDownsampleColorImages=false", "-dDownsampleGrayImages=false", "-dDownsampleMonoImages=false", "-dAutoFilterColorImages=false",
           "-dAutoFilterGrayImages=false", "-dColorImageFilter=/FlateEncode", "-dGrayImageFilter=/FlateEncode", "-dColorConversionStrategy=/LeaveColorUnchanged",
           "-dEmbedAllFonts=true", "-dSubsetFonts=true", "-dPreserveAnnots=false", f"-sOutputFile={out}", pdf]
    if os.path.exists(out): os.remove(out)
    with RenderLock():
        rc, st = wd(name, cmd, max_s=max_s)
    o = wd_out(name)
    assert rc == 0 and os.path.exists(out) and os.path.getsize(out) > 0 and "error" not in o.lower(), ("gs pdfwrite failed", name, st, o[-400:])
    return st


def crop(png, rect_pt, dpi, out_path, label):
    from PIL import Image, ImageDraw
    S = dpi / 72.0; x0, y0, x1, y1 = [int(round(v * S)) for v in rect_pt]
    im = Image.open(png).convert("RGB"); c = im.crop((max(x0, 0), max(y0, 0), min(x1, im.width), min(y1, im.height)))
    canvas = Image.new("RGB", (c.width + 20, c.height + 26), (255, 255, 255)); dr = ImageDraw.Draw(canvas)
    dr.text((10, 4), label, fill=(200, 0, 0)); canvas.paste(c, (10, 22)); canvas.save(out_path); return out_path
