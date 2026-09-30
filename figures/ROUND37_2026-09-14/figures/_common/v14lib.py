#!/usr/bin/env python3
"""V14 figure lanes (round 37, 2026-09-15): shared helpers. Import with
    import sys; sys.path.insert(0, "$T90_FIGURE_ROOT/ROUND37_2026-09-14/figures/_common"); from v14lib import *

Provides: paths (T90, NUM, SV, V13, LOCK_ROOT, FIG), the eviction gate (blocks, hydrated), sha256, the v8.1 sidecar gate
(sidecar_v8), text-layer probes (spans, collapse, multisets, letters, numeric tokens), the shared render (render, RenderLock),
the checks writer (Checks), CHANGES csv writer (write_changes), the printed-values CSV (printed_values_csv) with a source
resolver for JSON pointers and CSV lookups, and the half-up printing helpers (hr_ci, fmt_halfup).
PyMuPDF only. Never get_drawings() on Figures 1, 3, 4, 5.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../..")))  # repository root, where paths.py lives
import paths
import csv, fcntl, gc, hashlib, json, os, re, subprocess, sys, time
from collections import Counter
from decimal import Decimal, ROUND_HALF_UP
import fitz

T90 = paths.FIGURE_ROOT
NUM = f"{paths.NUMBERS_DIR}"
SV = paths.SV_ROOT
V8 = f"{paths.V8_ROOT}"
V13 = f"{paths.FIGURE_ROOT}/ROUND36_2026-09-10/L9_assembly/NEW_FINAL_SET_V13"
OLD_V7 = f"{V8}/compare/old_numbers"          # v7 numbers snapshot (OLD side of a delta only)
OLD_V8_0 = f"{V8}/compare/old_numbers_v8_0"   # v8.0 snapshot (OLD side of a delta only)
FIG = f"{paths.FIGURE_ROOT}/ROUND37_2026-09-14/figures"
LOCK_ROOT = f"{paths.FIGURE_ROOT}/ROUND30_2026-09-08/_render.lock"
V8_1_TABLE_BUILD = "2026-09-14 12:35:00"      # the v8.1 tables' build time: a v8.1 output is dated at or after this
ORDER = [f"Main_Fig{i}" for i in range(1, 7)] + [f"ED_Fig{i:02d}" for i in range(1, 11)] + [f"Supp_Fig{i:02d}" for i in range(1, 22)]
PALETTE = {"#0288d1", "#b3dcf2", "#7cc0e9", "#3f9fd8", "#2e7d32", "#298d32", "#39c445", "#79d475", "#5dcf66", "#d5f1d0", "#d55e00",
           "#f6d6c2", "#eeb08f", "#e4864f", "#1a1d21", "#1a1d20", "#8a9099", "#ccd1d6", "#eef0f1", "#ffffff", "#000000", "#3c4046",
           "#fdecea", "#d32f2f", "#0187d1"}
DRAWN_RECORD = {"text": "the string exactly as drawn", "x": "page x of the text origin (pt)", "baseline": "page y of the baseline (pt, top-down)",
                "ha": "left|center|right", "size": "font size pt", "source_file": "absolute path of the numbers file (or 'static:V13')",
                "source_key": "JSON pointer 'a/b/0/c' or CSV lookup 'col=value:field[,field2]' or free text", "source_value": "the raw value(s) read",
                "rule": "printing rule: 2dp_halfup | 1dp_halfup | int_comma | pct_1dp | q3 | text | p3 | ...", "panel": "letter"}


# ---------------------------------------------------------------- files
def blocks(path):
    out = subprocess.run(["stat", "-f", "%b %z", path], capture_output=True, text=True, check=True).stdout.split()
    return int(out[0]), int(out[1])


def hydrated(path, wait_s=600):
    """iCloud eviction gate. A file with a nonzero size and zero allocated blocks is evicted: request it and poll (never read it evicted)."""
    if not os.path.exists(path): raise FileNotFoundError(path)
    b, z = blocks(path)
    if z > 0 and b == 0:
        subprocess.run(["brctl", "download", path], capture_output=True)
        t0 = time.time()
        while blocks(path)[0] == 0:
            if time.time() - t0 > wait_s: raise RuntimeError(f"EVICTED and not downloaded in {wait_s} s: {path}")
            time.sleep(2)
    return path


def sha256(path, n=None):
    h = hashlib.sha256()
    with open(hydrated(path), "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""): h.update(chunk)
    return h.hexdigest() if n is None else h.hexdigest()[:n]


def load_json(path): return json.load(open(hydrated(path)))


_PSV = None
def psv_steps():
    """id -> (group, name, declared outputs) from rerun_steps_v8.psv (field 6 = declared outputs)."""
    global _PSV
    if _PSV is None:
        _PSV = {}
        for l in open(hydrated(f"{V8}/rerun_steps_v8.psv")):
            if not l.strip() or l.startswith("#"): continue
            f = [x.strip() for x in l.rstrip("\n").split("|")]
            if len(f) < 7 or f[0] == "id": continue
            _PSV[f[0]] = (f[1], f[2], [p for p in f[5].split(";") if p])
    return _PSV


def status_ok_rows():
    """(id) -> latest ok end time from STATUS.tsv."""
    out = {}
    for l in open(hydrated(f"{V8}/logs/STATUS.tsv")):
        f = l.rstrip("\n").split("\t")
        if len(f) < 8 or f[0] == "id": continue
        if f[5] == "0" and f[7] in ("ok", "MANUAL"): out[f[0]] = max(out.get(f[0], ""), f[4])
    return out


def sidecar_v8(path):
    """The provenance sidecar beside a numbers file must exist and prove the file is a v8.1 output: inputs name data_frozen_v8_2026-09,
    or its step id has a STATUS ok row at or after the v8.1 table build. Returns a dict for the provenance record. Raises otherwise."""
    sc_path = path + ".provenance.json"
    if not os.path.exists(sc_path): raise SystemExit(f"NO SIDECAR: {path} (a printed number without a sidecar is not printable)")
    sc = load_json(sc_path); ins = sc.get("inputs", {})
    paths = [v.get("path", "") if isinstance(v, dict) else str(v) for v in (ins.values() if isinstance(ins, dict) else ins)]
    v8in = any("data_frozen_v8_2026-09" in p for p in paths); v7in = any("data_frozen_v7_2026-09" in p for p in paths)
    st = sc.get("step", {}) or {}; sid = str(st.get("id", "")); mtime = (sc.get("output", {}) or {}).get("mtime_local", "")
    ok_rows = status_ok_rows(); ok_time = ok_rows.get(sid, "")
    fresh = bool(ok_time and ok_time >= V8_1_TABLE_BUILD) or (mtime >= V8_1_TABLE_BUILD if mtime else False)
    if v7in and not v8in: raise SystemExit(f"V7 SIDECAR (group P or stale): {path} step {sid} inputs {paths[:3]}")
    if not (v8in or fresh): raise SystemExit(f"SIDECAR NOT v8.1: {path} step {sid} out {mtime} status ok {ok_time}")
    osha = sc.get("output_sha256") or (sc.get("output", {}) or {}).get("sha256")
    sha = sha256(path)
    if osha and osha != sha: raise SystemExit(f"SIDECAR SHA MISMATCH: {path} sidecar {osha[:16]} file {sha[:16]}")
    return dict(path=path, step=sid, name=st.get("name"), output_mtime=mtime, status_ok=ok_time, sha256=sha, v8_inputs=v8in)


# ---------------------------------------------------------------- text layer
def spans(path, page_no=0):
    """[(text, x, y, font, size, colour_int, bbox)] of every span with a non-empty text, origin in top-down page points."""
    d = fitz.open(hydrated(path)); out = []
    for b in d[page_no].get_text("dict")["blocks"]:
        for l in b.get("lines", []):
            for s in l["spans"]:
                t = s["text"].replace("\xa0", " ").strip()
                if t: out.append((t, s["origin"][0], s["origin"][1], s["font"], s["size"], s["color"], tuple(s["bbox"])))
    d.close(); gc.collect(); return out


def collapse(sp, tol=1.5):
    keep = []
    for rec in sp:
        t, x, y = rec[0], rec[1], rec[2]
        if not any(t == k[0] and abs(x - k[1]) < tol and abs(y - k[2]) < tol for k in keep): keep.append(rec)
    return keep


NUMRE = re.compile(r"-?\d[\d,]*\.?\d*")
def multisets(sp):
    c = collapse(sp)
    return Counter(w for r in c for w in r[0].split()), Counter(x for r in c for x in NUMRE.findall(r[0]))


def letters(sp):
    return sorted([(r[0], round(r[1], 2), round(r[2], 2)) for r in collapse(sp) if len(r[0]) == 1 and r[0] in "abcdefgh" and r[4] >= 12], key=lambda t: (t[2], t[1]))


def page_box(path):
    d = fitz.open(hydrated(path)); n = len(d); w, h = d[0].rect.width, d[0].rect.height; d.close(); gc.collect()
    return n, round(w, 3), round(h, 3)


def fonts_of(path):
    d = fitz.open(hydrated(path)); fs = sorted({f[3].split("+")[-1] for f in d[0].get_fonts(full=True)}); d.close(); gc.collect(); return fs


# ---------------------------------------------------------------- render
class RenderLock:
    """ONE heavy render at a time across lanes (the round-root lock every lane uses)."""
    def __enter__(self):
        self.f = open(LOCK_ROOT, "a"); fcntl.flock(self.f, fcntl.LOCK_EX); return self
    def __exit__(self, *a):
        try: fcntl.flock(self.f, fcntl.LOCK_UN); self.f.close()
        except Exception: pass


RENDER_CODE = "import fitz,sys\nd=fitz.open(sys.argv[1]);pm=d[0].get_pixmap(dpi=int(sys.argv[3]),colorspace=fitz.csRGB,alpha=False,clip=None if sys.argv[4]=='-' else fitz.Rect(*[float(v) for v in sys.argv[4].split(',')]));pm.save(sys.argv[2]);print(pm.width,pm.height)\n"
def render(pdf, dpi, out_png, heavy=None, clip=None, timeout=2400):
    """Render page 0 in a child process (memory returns to the OS). heavy defaults to: main figure or file over 1 MB or dpi >= 300."""
    hydrated(pdf)
    if heavy is None: heavy = os.path.basename(pdf).startswith("Main") or os.stat(pdf).st_size > 1_000_000 or dpi >= 300
    clip_s = "-" if clip is None else ",".join(str(v) for v in clip)
    def run():
        r = subprocess.run([sys.executable, "-c", RENDER_CODE, pdf, out_png, str(dpi), clip_s], capture_output=True, text=True, timeout=timeout)
        if r.returncode != 0: raise RuntimeError(f"render failed for {pdf}: {r.stderr.strip()[-400:]}")
        return [int(x) for x in r.stdout.split()]
    if heavy:
        with RenderLock(): return run()
    return run()


def raster_diff(png_a, png_b, allowed_boxes_pt, dpi, thresh=32, dilate_pt=1.5):
    """Pixels that differ (any channel > thresh) outside the allowed page-point boxes. Returns dict(n_diff, n_outside, shape_equal)."""
    import numpy as np
    from PIL import Image
    a = np.asarray(Image.open(png_a).convert("RGB")).astype(int); b = np.asarray(Image.open(png_b).convert("RGB")).astype(int)
    if a.shape != b.shape: return dict(n_diff=None, n_outside=None, shape_equal=False, shapes=[a.shape, b.shape])
    d = (np.abs(a - b) > thresh).any(axis=2); mask = np.zeros_like(d)
    s = dpi / 72.0
    for (x0, y0, x1, y1) in allowed_boxes_pt:
        X0, Y0 = max(int((x0 - dilate_pt) * s), 0), max(int((y0 - dilate_pt) * s), 0); X1, Y1 = min(int((x1 + dilate_pt) * s) + 1, d.shape[1]), min(int((y1 + dilate_pt) * s) + 1, d.shape[0])
        mask[Y0:Y1, X0:X1] = True
    return dict(n_diff=int(d.sum()), n_outside=int((d & ~mask).sum()), shape_equal=True)


# ---------------------------------------------------------------- checks and CHANGES
class Checks:
    def __init__(self, header):
        self.lines = [header, ""]; self.n_pass = 0; self.n_fail = 0
    def log(self, what, ok, detail=""):
        self.lines.append(f"{'PASS' if ok else 'FAIL'}: {what}" + (f" [{detail}]" if detail != "" else ""))
        if ok: self.n_pass += 1
        else: self.n_fail += 1
        print(self.lines[-1][:300], flush=True); return ok
    def write(self, path):
        ok = self.n_fail == 0
        self.lines += ["", f"{self.n_pass} pass, {self.n_fail} fail", "RESULT ALL PASS" if ok else "RESULT FAIL"]
        os.makedirs(os.path.dirname(path), exist_ok=True); open(path, "w").write("\n".join(self.lines) + "\n"); return ok


def write_changes(path, rows):
    """rows: dicts with sheet, panel, item, before, after, dramatic (yes/no), note."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=["sheet", "panel", "item", "before", "after", "dramatic", "note"], extrasaction="ignore"); wr.writeheader()
        for r in rows: wr.writerow({k: r.get(k, "") for k in wr.fieldnames})


# ---------------------------------------------------------------- printing rules (half up on the value given)
def halfup(v, nd):
    return str(Decimal(repr(float(v))).quantize(Decimal(1).scaleb(-nd), rounding=ROUND_HALF_UP))


def fmt_halfup(v, nd): return halfup(v, nd)
def hr_ci(hr, lo, hi, nd=2): return f"{halfup(hr, nd)} ({halfup(lo, nd)}-{halfup(hi, nd)})"
def int_comma(v): return f"{int(round(float(v))):,}"
def q_text(q):
    q = float(q)
    if q < 0.001: return "<0.001"
    if (q < 0.05) != (float(halfup(q, 3)) < 0.05): return halfup(q, 4)
    return halfup(q, 3)
def p_text(p):
    p = float(p)
    return "<0.001" if p < 0.001 else halfup(p, 3)
def p_2sig(p):
    """Two significant figures, '<0.001' below 0.001, one further decimal where rounding lands on the 0.050 boundary (the polish generators' fmt_p)."""
    from math import floor, log10
    p = float(p)
    if p < 0.001: return "<0.001"
    d = 1 - floor(log10(p)); s = f"{p:.{d}f}"
    while s in ("0.050", "0.0500") and abs(p - 0.05) > 1e-12 and d < 6:
        d += 1; s = f"{p:.{d}f}"
    return s
RULES = {"2dp_halfup": lambda v: halfup(v, 2), "1dp_halfup": lambda v: halfup(v, 1), "3dp_halfup": lambda v: halfup(v, 3), "int_comma": int_comma,
         "int": lambda v: str(int(round(float(v)))), "pct_1dp": lambda v: halfup(v, 1), "q3": q_text, "p3": p_text, "text": lambda v: str(v),
         "hr_ci_2dp": lambda v: hr_ci(*v), "hr_ci_1dp": lambda v: hr_ci(*v, nd=1), "q2sig": p_2sig, "p2sig": p_2sig,
         "2dp_plain": lambda v: f"{float(v):.2f}", "1dp_plain": lambda v: f"{float(v):.1f}", "3dp_plain": lambda v: f"{float(v):.3f}",
         "hr_ci_2dp_plain": lambda v: f"{float(v[0]):.2f} ({float(v[1]):.2f}-{float(v[2]):.2f})",
         # a JSON object with hr/lo/hi fields (one_exposure_v1.json outcomes and the like): the pointer resolves to the dict
         "hr_ci_2dp_dict": lambda v: hr_ci(v["hr"], v["lo"], v["hi"]), "hr_ci_1dp_dict": lambda v: hr_ci(v["hr"], v["lo"], v["hi"], nd=1),
         # lane V14_L6 (ED_Fig06): a group count printed "n = 1,623" from one CSV field, and "events/n" from two CSV fields (events_group, n_group)
         "n_eq_int_comma": lambda v: f"n = {int(round(float(v))):,}", "ev_n": lambda v: f"{int(round(float(v[0])))}/{int(round(float(v[1]))):,}",
         # lane V14_L5b (Supp 13, 14): odds ratio dicts of the nonresponder files (marginal_or/lo/hi) and the joint model (or/lo/hi)
         "or_ci_2dp_marginal": lambda v: hr_ci(v["marginal_or"], v["marginal_lo"], v["marginal_hi"]), "or_ci_2dp_orlohi": lambda v: hr_ci(v["or"], v["lo"], v["hi"]),
         "or_ci_2dp_marginal_star": lambda v: hr_ci(v["marginal_or"], v["marginal_lo"], v["marginal_hi"]) + "*"}

# additive rules for lane V14_L3b (Supp_Fig06 panel strings): "<n> deaths", "<rate> deaths per 1,000 person-years", BH stars from a q value
RULES.update({"deaths_int": lambda v: f"{int(round(float(v)))} deaths",
              "rate_1dp": lambda v: f"{halfup(v, 1)} deaths per 1,000 person-years",
              "stars_q": lambda q: ("***" if float(q) < 0.001 else "**" if float(q) < 0.01 else "*" if float(q) < 0.05 else "")})



# V14_L5b (2026-09-15, Supp_Fig13): two additive rules, the "all patients, 30.0%" reference-rule label and the panel b trend string
RULES.update({"all_patients_pct": lambda v: f"all patients, {halfup(v, 1)}%",
              "trend_or_ci_2dp": lambda v: f"{halfup(v['or'], 2)} (95% CI, {halfup(v['lo'], 2)}-{halfup(v['hi'], 2)})"})


def hr_ci_wald_p(hr, p, nd=2):
    """Supp_Fig20 panel b (V14_L5b, v8.2): the Wald interval recovered from a hazard ratio and its two-sided p (make_efig12_stage.py), printed like hr_ci."""
    import math
    from scipy.stats import norm
    coef = math.log(float(hr)); z = norm.isf(max(min(float(p), 1.0), 1e-300) / 2.0); se = abs(coef) / z if z > 0 else float("nan")
    return f"{float(hr):.{nd}f} ({math.exp(coef - 1.96 * se):.{nd}f}-{math.exp(coef + 1.96 * se):.{nd}f})"
RULES["hr_ci_wald_p_2dp"] = lambda v: hr_ci_wald_p(v[0], v[1])


# ---------------------------------------------------------------- source resolver
def _get_json(obj, pointer):
    for part in [p for p in pointer.split("/") if p != ""]:
        if isinstance(obj, list): obj = obj[int(part)]
        elif isinstance(obj, dict):
            if part in obj: obj = obj[part]
            else:
                # allow "key=value" selection inside a list of dicts handled above, or numeric keys
                obj = obj[int(part)] if part.lstrip("-").isdigit() and int(part) in obj else obj[part]
        else: raise KeyError(part)
    return obj


def resolve_source(source_file, key):
    """Read the raw value(s) named by key from source_file.
    JSON: key is a pointer 'a/b/0/c' (list elements by index) or 'a/b/[name=Heart failure]/hr' (select the list element whose field equals).
    CSV : key is 'col=value[&col2=value2]:field' or ':field' with several fields 'field1,field2,field3' -> list.
    Returns the value, a list for several fields, or raises."""
    hydrated(source_file)
    if source_file.endswith(".json"):
        obj = load_json(source_file); cur = obj
        for part in [p for p in key.split("/") if p != ""]:
            m = re.fullmatch(r"\[(\w+)=(.+)\]", part)
            if m and isinstance(cur, list):
                sel = [e for e in cur if str(e.get(m.group(1))) == m.group(2)]
                if len(sel) != 1: raise KeyError(f"{part}: {len(sel)} matches")
                cur = sel[0]
            elif isinstance(cur, list): cur = cur[int(part)]
            else: cur = cur[part]
        return cur
    if source_file.endswith(".csv") or source_file.endswith(".tsv"):
        import pandas as pd
        sel, fields = key.split(":", 1) if ":" in key else (key, "")
        df = pd.read_csv(source_file, comment="#") if source_file.endswith(".csv") else pd.read_csv(source_file, sep="\t", comment="#")
        for cond in [c for c in sel.split("&") if c]:
            col, val = cond.split("=", 1)
            df = df[df[col].astype(str) == val]
        if len(df) != 1: raise KeyError(f"{sel}: {len(df)} rows")
        fl = [f for f in fields.split(",") if f]
        vals = [df.iloc[0][f] for f in fl]
        vals = [v.item() if hasattr(v, "item") else v for v in vals]
        return vals[0] if len(vals) == 1 else vals
    if source_file.endswith(".parquet"):
        raise KeyError("parquet sources: give source_value and rule, resolution not automated")
    raise KeyError(f"unknown source type: {source_file}")


def printed_values_csv(sheet_pdf, drawn_records, out_csv, static_from=None, group_p_boxes=(), tol=0.6):
    """One row per text span on the sheet that carries a digit (plus every drawn text): string, source_file, key, match, x, y,
    source_value, rule, note. drawn_records: list of DRAWN_RECORD dicts (text, x, baseline, ha, source_file, source_key, source_value, rule).
    A span is matched to a drawn record by text and origin (x within tol after alignment, y within tol); 'match' is:
      yes            drawn text equal to the sheet span AND the rule applied to the value re-read from source_file/source_key gives the same string
      yes_recorded   drawn text equal to the span, source re-read not automated (no resolvable key or a parquet source): the recorded source_value under the rule gives the string
      static         not drawn by the lane and identical (text + origin within tol) to the V13 sheet (static_from = V13 sheet path)
      kept_v13_groupP  inside a group P box (page points) and identical to V13
      NO             anything else (a digit string the lane cannot account for, or a value that does not re-derive)
    Returns (n_rows, n_no)."""
    sp = collapse(spans(sheet_pdf)); base = collapse(spans(static_from)) if static_from else []
    base_idx = {}
    for r in base: base_idx.setdefault(r[0], []).append(r)
    dr_idx = {}
    for d in drawn_records: dr_idx.setdefault(str(d["text"]).strip(), []).append(d)
    rows = []; n_no = 0; used = set()
    def in_boxes(x, y): return any(b[0] - 1 <= x <= b[2] + 1 and b[1] - 1 <= y <= b[3] + 1 for b in group_p_boxes)
    for (t, x, y, font, size, col, bb) in sp:
        has_digit = any(ch.isdigit() for ch in t)
        cands = dr_idx.get(t, [])
        hit = None
        for d in cands:
            dx = float(d["x"]); ha = d.get("ha", "left")
            # compare the span's origin to the drawn anchor: left origin, or centre / right of the span bbox
            sx = bb[0] if ha == "left" else (bb[2] if ha == "right" else (bb[0] + bb[2]) / 2)
            if abs(sx - dx) <= max(tol, 0.02 * (bb[2] - bb[0])) + 1.0 and abs(y - float(d["baseline"])) <= tol + 0.6: hit = d; break
        if hit is None and cands:
            # fall back: text equal and baseline equal (or x equal: rotated text), the NEAREST unused record first (round 37 fix: a used
            # record must not swallow a second span of the same text, e.g. the x ticks "1", "50" of two panels on one baseline)
            pool = [d for d in cands if id(d) not in used and (abs(y - float(d["baseline"])) <= tol + 0.6 or abs(x - float(d["x"])) <= tol + 1.0)]
            if pool: hit = min(pool, key=lambda d: abs(((bb[0] + bb[2]) / 2) - float(d["x"])) + abs(y - float(d["baseline"])))
        if hit is not None:
            used.add(id(hit)); rule = hit.get("rule", "text"); sf = hit.get("source_file", ""); sk = hit.get("source_key", ""); sv = hit.get("source_value", "")
            match = "NO"; note = ""
            try:
                fn = RULES.get(rule)
                rv = None
                if rule == "keyname":                     # the printed string IS a JSON key: the pointer must resolve, the string = its last segment
                    try:
                        resolve_source(sf, sk); got = sk.rstrip("/").split("/")[-1]
                        match = "yes" if got == t else "NO"; note = "" if match == "yes" else f"pointer last segment {got!r}"
                    except Exception as e:
                        match = "NO"; note = f"pointer did not resolve: {type(e).__name__}: {str(e)[:80]}"
                    rows.append(dict(string=t, source_file=sf, key=sk, match=match, x=round(x, 2), y=round(y, 2), source_value="", rule=rule, note=note, panel=hit.get("panel", "")))
                    if match == "NO": n_no += 1
                    continue
                try:
                    rv = resolve_source(sf, sk) if sf and sk and not sf.startswith("static") else None
                except Exception as e:
                    rv = None; note = f"source not re-read: {type(e).__name__}: {str(e)[:80]}"
                if fn is not None and rv is not None:
                    got = fn(rv)
                    match = "yes" if str(got) == t or str(got) == str(hit["text"]).strip() else "NO"
                    if match == "NO": note = f"rule {rule} on re-read value {rv!r} gives {got!r}"
                elif fn is not None and sv not in ("", None):
                    got = fn(sv); match = "yes_recorded" if str(got) == t else "NO"
                    if match == "NO": note = (note + "; " if note else "") + f"rule {rule} on recorded value {sv!r} gives {got!r}"
                else:
                    match = "yes_recorded" if rule == "text" and str(hit["text"]).strip() == t else "NO"
                    if match == "NO": note = (note + "; " if note else "") + f"no rule/value ({rule})"
            except Exception as e:
                match = "NO"; note = f"{type(e).__name__}: {str(e)[:100]}"
            rows.append(dict(string=t, source_file=sf, key=sk, match=match, x=round(x, 2), y=round(y, 2), source_value=json.dumps(sv, default=str)[:120] if sv not in ("", None) else "", rule=rule, note=note, panel=hit.get("panel", "")))
        else:
            same_v13 = any(abs(x - r[1]) <= tol + 1.0 and abs(y - r[2]) <= tol + 0.6 for r in base_idx.get(t, []))
            if same_v13 and in_boxes(x, y): rows.append(dict(string=t, source_file="static:V13", key="kept group P panel", match="kept_v13_groupP", x=round(x, 2), y=round(y, 2), source_value="", rule="", note="", panel=""))
            elif same_v13: rows.append(dict(string=t, source_file="static:V13", key="identical to V13 at the same origin", match="static", x=round(x, 2), y=round(y, 2), source_value="", rule="", note="", panel=""))
            elif has_digit: rows.append(dict(string=t, source_file="", key="", match="NO", x=round(x, 2), y=round(y, 2), source_value="", rule="", note="digit-bearing string not in the drawn records and not on V13 at this origin", panel=""))
            else: rows.append(dict(string=t, source_file="", key="", match="NO", x=round(x, 2), y=round(y, 2), source_value="", rule="", note="string not in the drawn records and not on V13 at this origin", panel=""))
        if rows[-1]["match"] == "NO": n_no += 1
    # drawn records never found on the sheet
    for d in drawn_records:
        if id(d) not in used:
            rows.append(dict(string=str(d["text"]), source_file=d.get("source_file", ""), key=d.get("source_key", ""), match="NO", x=d.get("x", ""), y=d.get("baseline", ""), source_value="", rule=d.get("rule", ""), note="drawn record not found on the sheet's text layer (not read back)", panel=d.get("panel", ""))); n_no += 1
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    with open(out_csv, "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=["string", "source_file", "key", "match", "x", "y", "source_value", "rule", "note", "panel"]); wr.writeheader()
        for r in rows: wr.writerow(r)
    return len(rows), n_no


# ---------------------------------------------------------------- font metrics alignment (round-28 idiom)
def sheet_font_metrics(base_text_spans_dict):
    """asc/desc per mille per font name from a text dict of the base sheet: {font: (asc, desc)}. Give get_text('dict') spans with 'asc'/'desc'."""
    met = {}
    for s in base_text_spans_dict:
        met.setdefault(s["font"], (int(round(s["asc"] * 1000)), int(round(s["desc"] * 1000))))
    return met


def align_font_metrics(src, dst, met):
    """Set every Type 42 font descriptor's Ascent/Descent in a matplotlib PDF to the base sheet's own Arial descriptors."""
    doc = fitz.open(hydrated(src)); page = doc[0]; done = {}
    for f in page.get_fonts(full=True):
        base = f[3].split("+")[-1]
        if base not in met: continue
        asc, dsc = met[base]
        desc = doc.xref_get_key(f[0], "DescendantFonts")
        if desc[0] != "array": continue
        cid = int(re.search(r"(\d+) 0 R", desc[1]).group(1))
        fdx = int(doc.xref_get_key(cid, "FontDescriptor")[1].split()[0])
        doc.xref_set_key(fdx, "Ascent", str(asc)); doc.xref_set_key(fdx, "Descent", str(dsc)); done[base] = (asc, dsc)
    doc.save(dst, garbage=1, deflate=True); doc.close(); gc.collect(); return done


def text_dict_spans(path):
    """get_text('dict') spans with asc/desc, for sheet_font_metrics()."""
    d = fitz.open(hydrated(path)); out = []
    for b in d[0].get_text("dict")["blocks"]:
        for l in b.get("lines", []):
            for s in l["spans"]:
                if s["text"].strip(): out.append(dict(text=s["text"], font=s["font"], size=s["size"], origin=s["origin"], bbox=s["bbox"], asc=s["ascender"], desc=s["descender"], color=s["color"]))
    d.close(); gc.collect(); return out


if __name__ == "__main__":
    print("v14lib ok; V13 set", V13, os.path.isdir(V13))

# FIX A2 (lane V14_L2_T90, 2026-09-15): a full-precision checkpoint record {HR, lo, hi} (dict or [hr, lo, hi] list) printed half up at 2 dp
RULES["hr_ci_HRlohi"] = lambda v: hr_ci(v["HR"], v["lo"], v["hi"]) if isinstance(v, dict) else hr_ci(*v)

# ---------------------------------------------------------------- additive rules, lane V14_L3b_DURATION_SUPP (v8.2, Supp_Fig06 panel d heads)
# each takes the dict node the JSON pointer resolves to (summary.json summary/<analysis>/<contrast>, attack.json attacks/A4_selection/landmark_1y/<contrast>)
RULES.setdefault("sig_of_tested", lambda d: f"{int(d['n_significant'])} of {int(d['n_conditions_tested'])} significant,")
RULES.setdefault("fdr_after", lambda d: f"{int(d['n_significant_fdr'])} after FDR")
RULES.setdefault("lm_sig_of_n", lambda d: f"{int(d['n_sig'])} of {int(d['n'])} significant,")
RULES.setdefault("lm_fdr_after", lambda d: f"{int(d['n_fdr'])} after FDR")
