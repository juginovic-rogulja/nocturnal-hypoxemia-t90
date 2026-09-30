"""Lane V14_L6_HABITUAL_EXTERNAL (round 37, 2026-09-15), shared helpers: the round-30 L6 helper repointed to V14.
Changes from the R30 original (lane_common_PRE_V8_1.py beside this file): LANE = the V14 lane, BASE = the V13 set (the design
baseline), the sidecar gate requires the v8.1 tables (data_frozen_v8_2026-09) and an optional minimum output time (the habitual
files regenerated 09:53 today), OLDNUM/OLDSV point at the v7 snapshot under V8_RECALC compare/old_numbers and are the OLD side of a
delta only (never a live source). Eviction gate, text-layer readers, render lock, font-metric alignment, checks writer unchanged."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import collections, fcntl, gc, hashlib, json, os, re, subprocess, time
import fitz

T90 = paths.FIGURE_ROOT
R30 = f"{paths.FIGURE_ROOT}/ROUND30_2026-09-08"
FIG = f"{paths.FIGURE_ROOT}/ROUND37_2026-09-14/figures"
LANE = f"{FIG}/lanes/V14_L6_HABITUAL_EXTERNAL"
BASE = f"{paths.FIGURE_ROOT}/ROUND36_2026-09-10/L9_assembly/NEW_FINAL_SET_V13"      # the V13 design baseline
NUM = f"{paths.NUMBERS_DIR}"
SV = paths.SV_ROOT
FF = f"{paths.FIGURE_ROOT}/FINAL_FIGURES_2026-08-14/_scripts"
OLDNUM = f"{paths.V8_ROOT}/compare/old_numbers"                 # OLD side of a delta only (the v7 numbers snapshot), never a live source
OLDSV = f"{OLDNUM}/external$T90_SV_ROOT"
V8_1_TABLE_BUILD = "2026-09-14 12:35:00"
ARIAL = paths.FONT_ARIAL
ARIALB = paths.FONT_ARIAL_BOLD
LOCK = f"{R30}/_render.lock"
INK = "#1a1d21"


def blocks(path):
    return int(subprocess.run(["stat", "-f", "%b", path], capture_output=True, text=True).stdout.strip() or 0)


def hydrated(path, wait_s=900):
    """iCloud eviction gate: a file with size and 0 blocks is evicted. brctl download and poll (until loop)."""
    assert os.path.exists(path), f"missing: {path}"
    if os.path.getsize(path) == 0 or blocks(path) > 0:
        return path
    subprocess.run(["brctl", "download", path], check=False)
    t0 = time.time()
    while blocks(path) == 0:
        if time.time() - t0 > wait_s:
            raise RuntimeError(f"still evicted after {wait_s}s: {path}")
        time.sleep(2)
    return path


def sha256(path):
    h = hashlib.sha256()
    with open(hydrated(path), "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def md5(path):
    h = hashlib.md5()
    with open(hydrated(path), "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sidecar(path, need_date=None, min_mtime=None):
    """Load <path>.provenance.json and prove the file is a v8.1 output: inputs name data_frozen_v8_2026-09 (never the v7 tables),
    the output sha256 equals the file, the output time is at or after the v8.1 table build (2026-09-14 12:35) and, when given,
    after min_mtime (the habitual files regenerated at 09:53 on 2026-09-15). need_date (a run date) is accepted for compatibility."""
    sc = hydrated(path + ".provenance.json")
    d = json.load(open(sc))
    ins = json.dumps(d.get("inputs", {}))
    assert "data_frozen_v8_2026-09" in ins, f"sidecar of {path} does not name the v8 tables"
    assert "data_frozen_v7_2026-09" not in ins, f"sidecar of {path} names the v7 tables (group P or stale)"
    step = d.get("step", {})
    assert step.get("rc", 0) == 0, (path, step)
    out = d.get("output") or {}
    mt = str(out.get("mtime_local", ""))
    assert mt >= V8_1_TABLE_BUILD, f"sidecar of {path} is older than the v8.1 table build: {mt}"
    if min_mtime:
        assert mt > min_mtime, f"sidecar of {path} is dated {mt}, not after {min_mtime}"
    if need_date:
        assert mt[:10] == need_date, (path, mt)
    out_sha = out.get("sha256") or d.get("output_sha256")
    assert out_sha, f"sidecar of {path} carries no output sha256"
    assert out_sha == sha256(path), f"sidecar sha256 does not match the file: {path}"
    d["_gate"] = {"step_id": step.get("id"), "step_name": step.get("name"), "output_mtime": mt, "sha256": out_sha}
    return d


def spans_of(page):
    out = []
    for b in page.get_text("rawdict")["blocks"]:
        if b["type"] != 0:
            continue
        for l in b["lines"]:
            for s in l["spans"]:
                txt = "".join(ch["c"] for ch in s["chars"]).replace("\xa0", " ")   # Chromium exports carry NBSP (compose_v7 does the same)
                if not txt.strip():
                    continue
                org = s["chars"][0]["origin"]
                out.append(dict(text=txt, font=s["font"], size=round(s["size"], 3), color="#%06x" % s.get("color", 0),
                                bbox=[round(v, 3) for v in s["bbox"]], origin=[round(org[0], 3), round(org[1], 3)]))
    return out


def letters(page):
    return sorted([(s["text"].strip(), round(s["bbox"][0], 2), round(s["bbox"][1], 2), round(s["bbox"][3], 2), s["font"], round(s["size"], 1))
                   for s in spans_of(page) if len(s["text"].strip()) == 1 and s["text"].strip() in "abcdefgh" and s["size"] >= 12],
                  key=lambda t: (t[2], t[1]))


def dedupe(items, tol=0.25):
    """Collapse stacked duplicate copies (same text within tol pt), the round-28 rule."""
    kept = collections.defaultdict(list); out = []
    for it in items:
        key = it["text"]; x, y = it["bbox"][0], it["bbox"][1]
        if any(abs(x - kx) <= tol and abs(y - ky) <= tol for kx, ky in kept[key]):
            continue
        kept[key].append((x, y)); out.append(it)
    return out


def words_of(page):
    return [dict(text=w[4], bbox=[round(v, 3) for v in w[:4]]) for w in page.get_text("words")]


def nums(texts):
    t = []
    for w in texts:
        t += re.findall(r"\d[\d,]*\.?\d*", w)
    return sorted(t)


def delta(a, b):
    ca, cb = collections.Counter(a), collections.Counter(b)
    return dict(lost=dict(ca - cb), gained=dict(cb - ca))


class render_lock:
    """ONE heavy render at a time across lanes (fcntl.flock on the round-root lock)."""
    def __enter__(self):
        self.f = open(LOCK, "w"); fcntl.flock(self.f, fcntl.LOCK_EX); return self
    def __exit__(self, *a):
        fcntl.flock(self.f, fcntl.LOCK_UN); self.f.close()


def sheet_font_metrics(sheet_pdf):
    doc = fitz.open(sheet_pdf); page = doc[0]; met = {}
    for b in page.get_text("rawdict")["blocks"]:
        if b["type"] != 0:
            continue
        for l in b["lines"]:
            for s in l["spans"]:
                met.setdefault(s["font"], (int(round(s["ascender"] * 1000)), int(round(s["descender"] * 1000))))
    doc.close(); gc.collect()
    return met


def align_font_metrics(src, dst, sheet_pdf):
    """Give the matplotlib Type 42 fonts the Ascent and Descent the base sheet's text layer reports (round 28)."""
    met = sheet_font_metrics(sheet_pdf)
    doc = fitz.open(src); page = doc[0]; done = {}
    for f in page.get_fonts(full=True):
        base = f[3].split("+")[-1]
        if base not in met:
            continue
        asc, dsc = met[base]
        desc = doc.xref_get_key(f[0], "DescendantFonts")
        if desc[0] != "array":
            continue
        cid = int(re.search(r"(\d+) 0 R", desc[1]).group(1))
        fdx = int(doc.xref_get_key(cid, "FontDescriptor")[1].split()[0])
        doc.xref_set_key(fdx, "Ascent", str(asc)); doc.xref_set_key(fdx, "Descent", str(dsc)); done[base] = (asc, dsc)
    doc.save(dst, garbage=1, deflate=True); doc.close()
    return done


def render(pdf, png, dpi=150, clip=None, heavy=False):
    def _do():
        d = fitz.open(pdf); p = d[0]
        pm = p.get_pixmap(dpi=dpi, clip=fitz.Rect(*clip) if clip else None, colorspace=fitz.csRGB, alpha=False)
        pm.save(png); w, h = pm.width, pm.height; pm = None; d.close(); gc.collect(); return w, h
    if heavy:
        with render_lock():
            return _do()
    return _do()


def pix_array(page, dpi=150, clip=None):
    import numpy as np
    pm = page.get_pixmap(dpi=dpi, clip=fitz.Rect(*clip) if clip else None, colorspace=fitz.csRGB, alpha=False)
    a = np.frombuffer(pm.samples, np.uint8).reshape(pm.h, pm.w, 3).astype(int); pm = None
    return a


def side_by_side(a, b, path, gap=20):
    import numpy as np
    h = max(a.shape[0], b.shape[0])
    def pad(x):
        if x.shape[0] == h: return x
        return np.concatenate([x, np.full((h - x.shape[0], x.shape[1], 3), 255)], axis=0)
    side = np.concatenate([pad(a), np.full((h, gap, 3), 255), pad(b)], axis=1).astype(np.uint8)
    fitz.Pixmap(fitz.csRGB, side.shape[1], side.shape[0], side.tobytes(), False).save(path)


def hexrgb(h):
    return tuple(int(h[i:i + 2], 16) / 255 for i in (1, 3, 5))


def hx(c):
    return "#%02x%02x%02x" % tuple(int(round(v * 255)) for v in c)


def write_checks(lines, path, result_ok, fail_text=None):
    lines = list(lines) + ["RESULT " + ("ALL PASS" if result_ok else ("FAIL: " + (fail_text or "see FAIL lines")))]
    os.makedirs(os.path.dirname(path), exist_ok=True)
    open(path, "w").write("\n".join(lines) + "\n")
    return lines[-1]


def now():
    return time.strftime("%Y-%m-%d %H:%M")


def status_append(text):
    with open(f"{LANE}/STATUS.md", "a") as f:
        f.write(text.rstrip("\n") + "\n")
