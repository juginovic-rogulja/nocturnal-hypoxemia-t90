"""Lane V14_L3a_DURATION_MAIN (round 37, 2026-09-15), shared helpers: a repointed copy of the round-30 lane_common.py
(backup beside it: lane_common_PRE_V8_1.py). Changes: LANE = this V14 lane, BASE = the V13 set (the design baseline),
no V7 table folder, sidecar() = the v8.1 gate of figures/_common/v14lib.py (inputs name data_frozen_v8_2026-09 or a
STATUS ok row after the v8.1 table build), and load_ff() which installs this lane's repointed splitstyle copy before any
FINAL_FIGURES helper is imported (the original splitstyle asserts the v7 control panel and fails on the v8.1 json)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import collections, fcntl, gc, hashlib, json, os, re, subprocess, sys, time
import fitz

T90 = paths.FIGURE_ROOT
FIG = f"{paths.FIGURE_ROOT}/ROUND37_2026-09-14/figures"
LANE = f"{FIG}/lanes/V14_L3a_DURATION_MAIN"
SCRIPTS = f"{LANE}/scripts"
BASE = f"{paths.FIGURE_ROOT}/ROUND36_2026-09-10/L9_assembly/NEW_FINAL_SET_V13"        # V13 = the design baseline
NUM = f"{paths.NUMBERS_DIR}"
SV = paths.SV_ROOT
FF = f"{paths.FIGURE_ROOT}/FINAL_FIGURES_2026-08-14/_scripts"
V8 = f"{paths.TABLES_DIR}"
OLD_V7 = f"{paths.V8_ROOT}/compare/old_numbers/numbers"       # the v7 numbers snapshot: OLD side of a delta only
R30_LANE = f"{paths.FIGURE_ROOT}/ROUND30_2026-09-08/V7_L3a_DURATION_MAIN"             # the v7 records (OLD side of a delta only)
ARIAL = paths.FONT_ARIAL
ARIALB = paths.FONT_ARIAL_BOLD
LOCK = f"{paths.FIGURE_ROOT}/ROUND30_2026-09-08/_render.lock"
sys.path.insert(0, f"{FIG}/_common")
import v14lib as VL  # noqa: E402


def load_ff():
    """Install this lane's splitstyle copy (assertion data-driven) as the module every FINAL_FIGURES helper imports, then
    expose FF and the NF scripts folder on sys.path. Call before importing figure2_polish_common or supp_polishB_common."""
    NF = f"{paths.FIGURE_ROOT}/New_Figures/_NOT_THESE_workfiles/scripts"          # legend_capture, nooverlap (unchanged NF helpers)
    if "splitstyle" not in sys.modules:
        for p in (NF, SCRIPTS):                                      # SCRIPTS ends up first: this lane's splitstyle copy wins
            while p in sys.path: sys.path.remove(p)
            sys.path.insert(0, p)
        import splitstyle  # noqa: F401  (this lane's copy: scripts/splitstyle.py)
        assert os.path.dirname(os.path.abspath(sys.modules["splitstyle"].__file__)) == SCRIPTS, sys.modules["splitstyle"].__file__
    if FF not in sys.path: sys.path.insert(0, FF)


def blocks(path):
    return int(subprocess.run(["stat", "-f", "%b", path], capture_output=True, text=True).stdout.strip() or 0)


def hydrated(path, wait_s=600):
    """iCloud eviction gate: a file with size and 0 blocks is evicted. Download and poll (until loop)."""
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


def sidecar(path):
    """The v8.1 provenance gate (v14lib.sidecar_v8): returns dict(path, step, name, output_mtime, status_ok, sha256, v8_inputs)."""
    return VL.sidecar_v8(hydrated(path))


def spans_of(page):
    out = []
    for b in page.get_text("rawdict")["blocks"]:
        if b["type"] != 0:
            continue
        for l in b["lines"]:
            for s in l["spans"]:
                txt = "".join(ch["c"] for ch in s["chars"])
                if not txt.strip():
                    continue
                org = s["chars"][0]["origin"]
                out.append(dict(text=txt, font=s["font"], size=round(s["size"], 3), color=s.get("color"),
                                bbox=[round(v, 3) for v in s["bbox"]], origin=[round(org[0], 3), round(org[1], 3)]))
    return out


def letters(page):
    return sorted([(s["text"].strip(), round(s["bbox"][0], 2), round(s["bbox"][1], 2), round(s["bbox"][3], 2), s["font"], round(s["size"], 1))
                   for s in spans_of(page) if len(s["text"].strip()) == 1 and s["text"].strip() in "abcdefgh" and s["size"] >= 12],
                  key=lambda t: (t[2], t[1]))


def dedupe(items, tol=0.25):
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
        self.f = open(LOCK, "a"); fcntl.flock(self.f, fcntl.LOCK_EX); return self
    def __exit__(self, *a):
        fcntl.flock(self.f, fcntl.LOCK_UN); self.f.close()


def sheet_font_metrics(sheet_pdf):
    doc = fitz.open(hydrated(sheet_pdf)); page = doc[0]; met = {}
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
    """Render in a child process (memory returns to the OS), under the round-root lock when heavy."""
    w, h = VL.render(pdf, dpi, png, heavy=heavy, clip=clip)
    return w, h


def write_checks(lines, path, result_ok):
    lines = list(lines) + ["RESULT " + ("ALL PASS" if result_ok else "FAILED")]
    os.makedirs(os.path.dirname(path), exist_ok=True)
    open(path, "w").write("\n".join(lines) + "\n")
    return lines[-1]


def now():
    return time.strftime("%Y-%m-%d %H:%M")


def status(pct, msg):
    with open(f"{LANE}/STATUS.md", "a") as f:
        f.write(f"- {time.strftime('%H:%M')} ET, {pct} percent. {msg}\n")


def grep_proof(sheet_dir):
    """The live-path proof: no script in this sheet folder names a v7 table, snapshot or set (backups excepted)."""
    hits = []
    for root, _d, files in os.walk(f"{sheet_dir}/scripts"):
        for fn in files:
            if not fn.endswith(".py") or fn.endswith("_PRE_V8_1.py"): continue
            for i, line in enumerate(open(os.path.join(root, fn), errors="replace"), 1):
                if re.search(r"data_frozen_v7|NEW_FINAL_SET_R30|NEW_FINAL_SET_V7\b", line) and "OLD" not in line and "_PRE_V8_1" not in line:
                    hits.append(f"{fn}:{i}: {line.strip()[:120]}")
    return hits
