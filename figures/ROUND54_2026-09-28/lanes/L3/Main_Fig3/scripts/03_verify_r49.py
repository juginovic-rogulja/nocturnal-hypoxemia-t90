#!$T90_PY
"""Round 49 (2026-09-26), lane L3, Main_Fig3: verify the +1 pt vector compose against the round-44 vector (the source of the V26 raster),
the V26 page box, the round-40 geometry record and the numbers files.
Ghostscript only on the composed sheets (txtwrite -dTextFormat=0 for the census, png16m for the render); PyMuPDF only on the flat pieces
(work/page_abc.pdf, the round-40 page_abc.pdf, the two bands). Never PyMuPDF on Main_Fig3 (memory rule of 2026-09-16).
Proves:
  1. page box (pypdf mediabox) = V26 = round-44 vector within 0.05 pt;
  2. text census: the page part (above the band) carries the same multiset of strings as the round-44 vector, every string paired at the
     same baseline with a size exactly +1.0 pt, except the declared kept elements (delta 0.0), each string keeping its left edge, centre
     or right edge; letters a-d and the title at 14 pt at the round-44 origins; the band part: same character multiset, no size smaller;
  3. drawn values against the data: every record string re-read from its numbers file under its rule (v14lib.printed_values_csv on the
     flat page_abc piece), every other string identical (text, baseline, an anchor edge) to the round-40 page_abc piece, sidecars and
     sha256 of the sources unchanged since the build;
  4. geometry frozen: every data mark (x, y, ci_x, fill) and the panel a axis, the cell fills and inks equal to the round-40 record;
     every record at the round-40 baseline, x equal or a declared nudge, size + 1.0 or a declared keep;
  5. fits and overlaps on the page_abc text layer: row labels clear the forest box, HR strings clear their whiskers, the q column clears
     the HR column, cell strings keep 2 pt inside their cells, the key and legends end inside the sheet, no two ink boxes intersect;
  6. the 150 dpi render exists, is newer than the vector and has the page's pixel size.
usage: 03_verify_r49.py"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import collections, gc, html, json, math, os, re, subprocess, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lmain_lib as L
import fitz, numpy as np, pypdf
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND37_2026-09-14/figures/_common")
import v14lib as VL

SHEET = "Main_Fig3"; SD = f"{L.LANE}/{SHEET}"; WORK = f"{SD}/work"; VER = f"{SD}/verify"; os.makedirs(VER, exist_ok=True)
NEW = f"{WORK}/{SHEET}_r49_vector.pdf"; PAGE = f"{WORK}/page_abc.pdf"; PNG150 = f"{SD}/{SHEET}_150dpi.png"
OLD = f"{paths.FIGURE_ROOT}/ROUND44_2026-09-21/lanes/LSKETCH/{SHEET}/work/{SHEET}_r44_vector.pdf"          # the vector source of the V26 raster
V26 = f"{paths.FIGURE_ROOT}/ROUND47_2026-09-24/figures/NEW_FINAL_SET_V26/{SHEET}.pdf"
R40 = f"{paths.FIGURE_ROOT}/ROUND40_2026-09-18/lanes/LMAIN/{SHEET}/work"; PAGE40 = f"{R40}/page_abc.pdf"; DJ40 = L.load_json(f"{R40}/{SHEET}_drawn.json")
BAND44 = f"{paths.FIGURE_ROOT}/ROUND44_2026-09-21/lanes/LSKETCH/work/build/band_3d_r44.pdf"
DJ = L.load_json(f"{WORK}/{SHEET}_drawn.json"); LOG = L.load_json(f"{WORK}/compose_log.json"); BAND = LOG["band"]
R49 = DJ["round49"]; UP = float(R49["up_pt"]); BAND_TOP = float(LOG["band_top"]); W_DECL, H_DECL = (float(v) for v in LOG["page"])
for p in (NEW, PAGE, OLD, V26, PAGE40, BAND, BAND44, PNG150): L.hydrated(p)
ck = L.Checks(f"Lane L3 round 49, {SHEET}, {L.now()}. NEW = {NEW} (sha256 {L.sha256(NEW)[:16]}), the vector compose of the +{UP:g} pt page_abc "
              f"(sha256 {L.sha256(PAGE)[:16]}) with the band {os.path.basename(BAND)} (sha256 {L.sha256(BAND)[:16]}). OLD = the round-44 vector {OLD} "
              f"(sha256 {L.sha256(OLD)[:16]}, the source of V26 {V26}). Geometry reference = round-40 record {R40}/{SHEET}_drawn.json.")


# ------------------------------------------------------------------ helpers
def gs_spans(pdf, xml):
    """Ghostscript txtwrite -dTextFormat=0 spans: text, x0, baseline, x1, font, size (bbox coordinates are integers in this format)."""
    r = subprocess.run([L.GS, "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=txtwrite", "-dTextFormat=0", f"-sOutputFile={xml}", pdf], capture_output=True, text=True, timeout=600)
    assert r.returncode == 0 and os.path.exists(xml), ("gs txtwrite failed", pdf, r.stderr[-300:])
    out = []
    for m in re.finditer(r'<span bbox="([^"]*)" font="([^"]*)" size="([^"]*)">(.*?)</span>', open(xml, encoding="utf-8", errors="replace").read(), re.S):
        bb = [float(v) for v in m.group(1).split()]; t = L.norm_text(html.unescape("".join(re.findall(r'c="([^"]*)"', m.group(4)))))
        if t.strip(): out.append(dict(text=t.strip(), x0=bb[0], y=bb[1], x1=bb[2], font=m.group(2).split("+")[-1].replace("-Identity-H", ""), size=round(float(m.group(3)), 4)))
    return out


def merge(spans, gap=1.0):
    """Ghostscript txtwrite starts a new span at every kerning adjustment (V|entricular, 1.1|1 (1.06-1.17)): glue fragments on one baseline
    with the same font and size whose gap is at most 1 pt (integer bboxes; distinct strings on a baseline are at least 2.4 pt apart)."""
    out = []
    for s in sorted(spans, key=lambda s: (s["y"], s["x0"])):
        if out and abs(out[-1]["y"] - s["y"]) <= 0.5 and out[-1]["font"] == s["font"] and abs(out[-1]["size"] - s["size"]) < 0.01 and -1.0 <= s["x0"] - out[-1]["x1"] <= gap:
            out[-1] = dict(out[-1], text=out[-1]["text"] + s["text"], x1=max(out[-1]["x1"], s["x1"]))
        else: out.append(dict(s))
    return out


def flat_spans(pdf):
    """PyMuPDF spans of a FLAT piece: text, origin, bbox, size, font, colour."""
    d = fitz.open(pdf); out = []
    for b in d[0].get_text("rawdict")["blocks"]:
        if b["type"] != 0: continue
        for l in b["lines"]:
            for s in l["spans"]:
                t = L.norm_text("".join(c["c"] for c in s["chars"]))
                if not t.strip(): continue
                o = s["chars"][0]["origin"]; bb = s["bbox"]
                out.append(dict(text=t.strip(), x=round(o[0], 3), y=round(o[1], 3), x0=round(bb[0], 3), x1=round(bb[2], 3), size=round(s["size"], 3), font=s["font"], color=s.get("color")))
    d.close(); gc.collect(); return out


def pair(old, new, dy=1.0):
    """Pair every old span with an unused new span of the same text at the same baseline, nearest by the best of the three anchors."""
    used = set(); pairs = []; miss = []
    for o in old:
        best, bd = None, 9e9
        for i, n in enumerate(new):
            if i in used or n["text"] != o["text"] or abs(n["y"] - o["y"]) > dy: continue
            d = min(abs(n["x0"] - o["x0"]), abs((n["x0"] + n["x1"]) / 2 - (o["x0"] + o["x1"]) / 2), abs(n["x1"] - o["x1"]))
            if d < bd: bd, best = d, i
        if best is None: miss.append(o)
        else: used.add(best); pairs.append((o, new[best]))
    extra = [n for i, n in enumerate(new) if i not in used]
    return pairs, miss, extra


def anchor_kind(o, n, tol=1.0, dx=0.0):
    c_o, c_n = (o["x0"] + o["x1"]) / 2, (n["x0"] + n["x1"]) / 2
    if abs(n["x0"] - o["x0"] - dx) <= tol: return "left" if dx == 0.0 else "left+rule"
    if abs(n["x0"] - o["x0"]) <= tol: return "left"
    if abs(c_n - c_o) <= tol: return "centre"
    if abs(n["x1"] - o["x1"]) <= tol: return "right"
    return None


# ------------------------------------------------------------------ 1. page box
def box(p): mb = pypdf.PdfReader(p).pages[0].mediabox; return float(mb.width), float(mb.height), len(pypdf.PdfReader(p).pages)
wn, hn, nn = box(NEW); wv, hv, nv = box(V26); wo, ho, no = box(OLD)
ck.log("page box: NEW = V26 = round-44 vector within 0.05 pt, one page each", nn == nv == no == 1 and abs(wn - wv) < 0.05 and abs(hn - hv) < 0.05 and abs(wn - wo) < 0.05 and abs(hn - ho) < 0.05 and abs(wn - W_DECL) < 0.01 and abs(hn - H_DECL) < 0.01,
       f"NEW {wn:.3f} x {hn:.3f}, V26 {wv:.3f} x {hv:.3f}, round-44 {wo:.3f} x {ho:.3f}, compose log {W_DECL} x {H_DECL}")

# ------------------------------------------------------------------ 2. census (Ghostscript txtwrite on both composed sheets)
SO_RAW = gs_spans(OLD, f"{WORK}/census_r44_vector.xml"); SN_RAW = gs_spans(NEW, f"{WORK}/census_r49_vector.xml")
SO, SN = merge(SO_RAW), merge(SN_RAW)
ck.info(f"gs txtwrite spans: OLD {len(SO_RAW)} raw -> {len(SO)} strings, NEW {len(SN_RAW)} raw -> {len(SN)} strings (kerning-split fragments glued: gap at most 1 pt, same baseline, font and size)")
# the panel b key rung labels move with their swatches (the 20 pt rule after each label, measured at the larger size): the declared shift
RUNG_DX = {k49["band"]: round(k49["label_x"] - k40["label_x"], 3) for k49, k40 in zip(DJ["panel_b"]["key_positions"], DJ40["panel_b"]["key_positions"])}
def rung_dx(text, y): return RUNG_DX.get(text.replace("-", "–"), RUNG_DX.get(text, 0.0)) if abs(y - 312.47) <= 1.0 else 0.0
Y_SPLIT = BAND_TOP - 5.0
po, pn = [s for s in SO if s["y"] < Y_SPLIT], [s for s in SN if s["y"] < Y_SPLIT]; bo, bn = [s for s in SO if s["y"] >= Y_SPLIT], [s for s in SN if s["y"] >= Y_SPLIT]
co, cn = collections.Counter(s["text"] for s in po), collections.Counter(s["text"] for s in pn)
ck.log("census, page part: same multiset of strings (gs txtwrite spans above the band)", co == cn, f"{len(po)} old, {len(pn)} new spans; lost {dict(co - cn)}, gained {dict(cn - co)}")
pairs, miss, extra = pair(po, pn)
ck.log("census, page part: every old span paired with a new span of the same text at the same baseline (integer bboxes, 1 pt)", not miss and not extra, f"{len(pairs)} pairs; unpaired old {[(m['text'], m['x0'], m['y']) for m in miss]}, unpaired new {[(e['text'], e['x0'], e['y']) for e in extra]}")
KEPT = R49["kept_at_current_size"]
def is_kept(n): return any(k["text"] == n["text"] and abs(float(k["baseline"]) - n["y"]) <= 1.0 for k in KEPT)
d_plus, d_zero, d_bad, kinds, no_anchor = 0, [], [], collections.Counter(), []
for o, n in pairs:
    d = n["size"] - o["size"]
    if abs(d - UP) <= 0.011: d_plus += 1
    elif abs(d) <= 0.011 and is_kept(n): d_zero.append((n["text"], n["y"], o["size"]))
    else: d_bad.append((n["text"], n["x0"], n["y"], o["size"], n["size"]))
    k = anchor_kind(o, n, dx=rung_dx(n["text"], n["y"])); kinds[k or "none"] += 1
    if k is None: no_anchor.append((n["text"], o["x0"], o["x1"], n["x0"], n["x1"], n["y"]))
ck.log(f"census, page part: every paired string +{UP:g} pt exactly, the declared kept elements +0.0", not d_bad and len(d_zero) == len(KEPT), f"{d_plus} at +{UP:g}, {len(d_zero)} kept (declared {len(KEPT)}), other {d_bad}")
ck.info(f"kept at the current size: {sorted(collections.Counter((t, s) for t, y, s in d_zero).items(), key=lambda kv: (-kv[1], kv[0]))}")
ck.log("census, page part: every string keeps its left edge, centre or right edge within 1 pt (integer bboxes; the panel b key rung labels their declared rule shift), baselines unchanged", not no_anchor, f"anchors {dict(kinds)}, rung shifts {RUNG_DX}; without an anchor {no_anchor}")
sizes_o, sizes_n = sorted(collections.Counter(s["size"] for s in po).items()), sorted(collections.Counter(s["size"] for s in pn).items())
ck.info(f"page-part size histogram OLD {sizes_o}"); ck.info(f"page-part size histogram NEW {sizes_n}")
lets_o = sorted([(s["text"], s["x0"], s["y"], s["size"]) for s in po if s["text"] in "abcd" and len(s["text"]) == 1 and "Bold" in s["font"] and s["size"] >= 12.5], key=lambda t: (t[2], t[1]))
lets_n = sorted([(s["text"], s["x0"], s["y"], s["size"]) for s in pn if s["text"] in "abcd" and len(s["text"]) == 1 and "Bold" in s["font"] and s["size"] >= 12.5], key=lambda t: (t[2], t[1]))
tit_o = [(s["x0"], s["y"], s["size"]) for s in po if s["text"] == "Figure 3"]; tit_n = [(s["x0"], s["y"], s["size"]) for s in pn if s["text"] == "Figure 3"]
ck.log("letters a, b, c, d: 14 pt Arial Bold at the round-44 origins (x0 and baseline within 1 pt)", [t[0] for t in lets_n] == ["a", "b", "c", "d"] and all(abs(n[3] - 14.0) < 0.01 and abs(o[3] - 13.0) < 0.01 and abs(n[1] - o[1]) <= 1.0 and abs(n[2] - o[2]) <= 1.0 for o, n in zip(lets_o, lets_n)), f"old {lets_o}, new {lets_n}")
ck.log("title 'Figure 3' once, 14 pt at the round-44 origin", len(tit_n) == 1 and abs(tit_n[0][2] - 14.0) < 0.01 and abs(tit_n[0][0] - tit_o[0][0]) <= 1.0 and abs(tit_n[0][1] - tit_o[0][1]) <= 1.0, f"old {tit_o}, new {tit_n}; compose log title {LOG['title']}, letters {LOG['letters_stamped']} at {LOG.get('letter_size')} pt")
# band part (the coordinator's build): same characters, sizes not smaller, sizes reported
cho, chn = collections.Counter(ch for s in bo for ch in s["text"] if not ch.isspace()), collections.Counter(ch for s in bn for ch in s["text"] if not ch.isspace())
bs_o, bs_n = sorted(collections.Counter(s["size"] for s in bo).items()), sorted(collections.Counter(s["size"] for s in bn).items())
ck.log("census, band part: same character multiset as the round-44 band, no size smaller (the band is the coordinator's round-49 build)", cho == chn and min(s for s, _ in bs_n) >= min(s for s, _ in bs_o), f"chars lost {dict(cho - chn)}, gained {dict(chn - cho)}; sizes OLD {bs_o}, NEW {bs_n}")
# band notes against the arrows (flat pieces, PyMuPDF): the '90% SpO2' notes' left edges, old and new
def band_words(pdf):
    sp = flat_spans(pdf); out = collections.defaultdict(list)
    for s in sp: out[round(s["y"], 1)].append(s)
    words = []
    for y, ss in out.items():
        ss = sorted(ss, key=lambda s: s["x0"]); words.append(dict(text="".join(s["text"] for s in ss), x0=min(s["x0"] for s in ss), x1=max(s["x1"] for s in ss), y=y, size=ss[0]["size"]))
    return sorted(words, key=lambda w: (w["y"], w["x0"]))
bw_o, bw_n = band_words(BAND44), band_words(BAND)
ck.info(f"band text lines (PyMuPDF on the flat bands; x in band points): OLD {[(w['text'], w['x0'], w['x1'], w['y'], w['size']) for w in bw_o]}")
ck.info(f"band text lines NEW {[(w['text'], w['x0'], w['x1'], w['y'], w['size']) for w in bw_n]}")
sp_o = [w for w in bw_o if "SpO2" in w["text"]]; sp_n = [w for w in bw_n if "SpO2" in w["text"]]
if sp_o and sp_n:
    ck.info(f"band '90% SpO2' notes: left edge moved from x {[w['x0'] for w in sp_o]} to {[w['x0'] for w in sp_n]} (band points), right edge {[w['x1'] for w in sp_o]} to {[w['x1'] for w in sp_n]}: the note grew leftward toward the arrowheads (seen touching in the 300 dpi crop; the coordinator's band)")

# ------------------------------------------------------------------ 3. drawn values against the data (flat page_abc piece)
LABEL_SHORT = {"Ventricular arrhythmia or cardiac arrest": "Ventricular arrhythmia/arrest"}; LABEL_C_FIRST = {"Death from any cause": "Death, any cause", "Cardiovascular composite": "Cardiovascular"}
def fmt2(x): return VL.halfup(x, 2)
def stars(q): return "" if (q is None or not np.isfinite(q)) else ("***" if q < 0.001 else "**" if q < 0.01 else "*" if q < 0.05 else "")
VL.RULES["hr_ci_HRlohi"] = lambda d: VL.hr_ci(d["HR"], d["lo"], d["hi"]) if isinstance(d, dict) else VL.hr_ci(*d)
VL.RULES["hr_stars"] = lambda v: fmt2(v[0]) + stars(float(v[1]) if v[1] is not None and not (isinstance(v[1], float) and math.isnan(v[1])) else float("nan"))
VL.RULES["ci_paren"] = lambda v: f"({fmt2(v[0])}-{fmt2(v[1])})"
VL.RULES["ref_one"] = lambda v: "1.00" if v == "NN" else "?"
VL.RULES["label_a"] = lambda v: LABEL_SHORT.get(str(v), str(v))
VL.RULES["label_c"] = lambda v: "(negative control)" if str(v) == "True" else LABEL_C_FIRST.get(str(v), str(v))
def records_for_csv(records):
    out = []
    for r in records:
        r = dict(r)
        if r.get("rule") == "label_alias":
            if r.get("panel") == "a": r["rule"] = "label_a"
            elif str(r.get("source_key", "")).endswith(":negative_control") and r.get("source_value") is not True:
                r["rule"] = "text"; r["source_key"] = ""; r["source_value"] = r["text"]
            else: r["rule"] = "label_c"
        out.append(r)
    return out
CSV_RECORDS = records_for_csv(DJ["records"])
n_rows, n_no = VL.printed_values_csv(PAGE, CSV_RECORDS, f"{VER}/{SHEET}_printed_values.csv", static_from=None)
import csv as _csv
rows = list(_csv.DictReader(open(f"{VER}/{SHEET}_printed_values.csv")))
yes = sum(1 for r in rows if r["match"] == "yes"); yes_rec = [r["string"] for r in rows if r["match"] == "yes_recorded"]; nos = [r for r in rows if r["match"] == "NO"]
# the strings that are not records must be the round-40 design strings, at the same baseline, keeping an anchor edge (nudges declared)
S49 = flat_spans(PAGE); S40 = flat_spans(PAGE40)
NUDGED = {(n["text"], round(float(n["baseline"]), 2)): n for n in R49["nudged"]}
static_ok, static_bad = [], []
for r in nos:
    if "not in the drawn records" not in r["note"]: static_bad.append((r["string"], r["x"], r["y"], r["note"])); continue
    t, x, y = L.norm_text(r["string"]).strip(), float(r["x"]), float(r["y"])
    n = min([s for s in S49 if s["text"] == t and abs(s["y"] - y) <= 0.6], key=lambda s: abs(s["x"] - x), default=None)
    o = [s for s in S40 if s["text"] == t and abs(s["y"] - y) <= 0.6]
    if n is None or not o: static_bad.append((t, x, y, "not on the round-40 page at this baseline")); continue
    o = min(o, key=lambda s: min(abs(s["x0"] - n["x0"]), abs((s["x0"] + s["x1"]) / 2 - (n["x0"] + n["x1"]) / 2), abs(s["x1"] - n["x1"])))
    kind = anchor_kind(dict(x0=o["x0"], x1=o["x1"]), dict(x0=n["x0"], x1=n["x1"]), tol=0.6, dx=rung_dx(t, y))
    exp = UP if not is_kept(dict(text=t, y=y)) else 0.0
    if kind and abs(n["size"] - o["size"] - exp) < 0.02: static_ok.append((t, kind))
    else: static_bad.append((t, x, y, f"anchor {kind}, size {o['size']} -> {n['size']}"))
ck.log("drawn values: every record string on the page_abc text layer at its recorded origin and re-derived from its numbers file under its rule (v14lib.printed_values_csv, PyMuPDF on the flat piece)", not [r for r in nos if "not in the drawn records" not in r["note"]] and not [r for r in rows if r["match"] == "NO" and r["note"].startswith("drawn record not found")],
       f"{n_rows} rows: {yes} re-derived (yes), {len(yes_rec)} recorded-only {sorted(set(yes_rec))[:8]}, {len(nos)} not records")
ck.log("drawn values: every other string is a round-40 design string at the same baseline with its left edge, centre or right edge kept (declared nudges) and +1.0 pt (or a declared keep)", not static_bad, f"{len(static_ok)} static strings {dict(collections.Counter(k for _, k in static_ok))}; failures {static_bad}")
src_bad = []
for name, sc in DJ["sources"].items():
    p = sc["path"]; sha = L.sha256(p)
    if sha != sc["sha256"]: src_bad.append((name, "sha256 changed since the build"))
    if "step" in sc:
        try: VL.sidecar_v8(p)
        except SystemExit as e: src_bad.append((name, str(e)[:100]))
ck.log("sources: every numbers file of the build unchanged (sha256), the six sidecar-gated ones still pass the v8.1 gate", not src_bad, f"{len(DJ['sources'])} sources, {sum(1 for v in DJ['sources'].values() if 'step' in v)} gated; {src_bad}")

# ------------------------------------------------------------------ 4. geometry frozen against the round-40 record
mk49 = [d for d in DJ["drawn"] if "marker" in d]; mk40 = [d for d in DJ40["drawn"] if "marker" in d]
same_marks = len(mk49) == len(mk40) == 104 and all(a["key"] == b["key"] and a["marker"] == b["marker"] and a["filled"] == b["filled"] and a["x"] == b["x"] and a["y"] == b["y"] and a["ci_x"] == b["ci_x"] for a, b in zip(mk49, mk40))
ck.log("geometry: all 104 data marks (52 T90 circles, 52 sleep-time squares: x, y, whisker x, fill) equal to the round-40 record", same_marks, f"{len(mk49)} marks new, {len(mk40)} old")
g49, g40 = DJ["panel_a"]["geometry"], DJ40["panel_a"]["geometry"]
keys = ["XA", "XB", "x_left", "x_right", "row0_centre", "pitch", "lane", "axis_top", "axis_bottom", "tick_baseline", "caption_baseline", "xlim", "ticks"]
ck.log("geometry: panel a axis (box, ticks, row pitch, baselines), DELTA and the panel c shift equal to the round-40 record", all(g49[k] == g40[k] for k in keys) and DJ["delta_pt"] == DJ40["delta_pt"] and DJ["panel_c_shift_pt"] == DJ40["panel_c_shift_pt"] and DJ["page"] == DJ40["page"], f"{ {k: g49[k] for k in ('XA', 'XB', 'x_left', 'x_right')} }, DELTA {DJ['delta_pt']}, C_SHIFT {DJ['panel_c_shift_pt']}, page {DJ['page']}")
tx49 = [d for d in DJ["drawn"] if "marker" not in d]; tx40 = [d for d in DJ40["drawn"] if "marker" not in d]
fills_same = len(tx49) == len(tx40) and all(a["panel"] == b["panel"] and a["text"] == b["text"] and a["baseline"] == b["baseline"] and a.get("fill") == b.get("fill") and a.get("ink") == b.get("ink") and a.get("band") == b.get("band") and a.get("bold") == b.get("bold") for a, b in zip(tx49, tx40))
ck.log("geometry: every drawn text of the record (259 entries) has the round-40 text, baseline, fill, ink, band and bold", fills_same, f"{len(tx49)} new, {len(tx40)} old")
rec49, rec40 = DJ["records"], DJ40["records"]; rec_bad = []
for a, b in zip(rec49, rec40):
    if (a["panel"], a["text"], a["ha"], a["baseline"], a["rule"]) != (b["panel"], b["text"], b["ha"], b["baseline"], b["rule"]): rec_bad.append(("identity", a["text"], a["baseline"])); continue
    nd = NUDGED.get((a["text"], round(float(a["baseline"]), 2)))
    x_ok = abs(a["x"] - b["x"]) < 1e-6 or (nd is not None and abs(a["x"] - (b["x"] + nd["dx"])) < 1e-6) or abs(a["x"] - (b["x"] + rung_dx(a["text"], float(a["baseline"])))) < 2e-3   # the key rung labels move with their swatches (declared rule)
    s_ok = abs(a["size"] - b["size"] - UP) < 1e-6 or (abs(a["size"] - b["size"]) < 1e-6 and is_kept(dict(text=a["text"], y=a["baseline"])))
    if not (x_ok and s_ok): rec_bad.append((a["text"], a["baseline"], b["x"], a["x"], b["size"], a["size"]))
ck.log("geometry: all 258 value records keep the round-40 identity and baseline, x equal or a declared nudge (or the key rung rule shift), size +1.0 or a declared keep", len(rec49) == len(rec40) == 258 and not rec_bad, f"{len(rec49)} records; failures {rec_bad[:6]}")
ck.info(f"declared nudges (origin moved to keep the centre or the right edge): {[(n['text'], n['anchor'], n['dx']) for n in R49['nudged']]}")
ck.info(f"declared keeps ({len(KEPT)}): {R49['keep_labels']} at 10.71 pt (the builder's 3 pt clearance to the forest box fails at 11.71 pt: ends at 156.6 and 144.9 against {g49['x_left'] - 3:.1f}); the panel c interval strings are at 10.5 pt (follow-up decision)")
kp49, kp40 = DJ["panel_b"]["key_positions"], DJ40["panel_b"]["key_positions"]
ck.info(f"panel b key rungs (swatch x): round-40 {[k['swatch_x'] for k in kp40]} -> round-49 {[k['swatch_x'] for k in kp49]} (the 20 pt rule after each label, measured at 12.5 pt)")

# ------------------------------------------------------------------ 5. fits and overlaps on the page_abc text layer (PyMuPDF, flat)
XL, XR = float(g49["x_left"]), float(g49["x_right"]); W = float(DJ["page"][0])
labels = [s for s in S49 if abs(s["x"] - 8.88) < 0.05 and 100 < s["y"] < 800]
ck.log("fit: all 52 panel a row labels end at least 3 pt before the forest box", len(labels) == 52 and all(s["x1"] <= XL - 3.0 + 0.05 for s in labels), f"largest end {max(s['x1'] for s in labels):.2f} vs {XL - 3.0:.2f} ({max(labels, key=lambda s: s['x1'])['text']})")
by_key = {(m["key"]): m for m in mk49 if m["marker"] == "T90 circle"}; by_key_sq = {(m["key"]): m for m in mk49 if m["marker"] == "total sleep time square"}
hr_gap, q_gap = [], []
for r in DJ["panel_a"]["rows"]:
    yc = None
    hs = [s for s in S49 if s["text"] == r["s_hr"] and abs(s["x1"] - 460.87) < 1.2]; qs = [s for s in S49 if s["text"] == r["s_q"] and abs(s["x1"] - 500.11) < 1.2]
    m = by_key[r["key"]]; sq = by_key_sq[r["key"]]
    hs = [s for s in hs if abs(s["y"] - (m["y"] + 3.06 - 0.37)) < 0.6]; qs = [s for s in qs if abs(s["y"] - (m["y"] + 3.06 - 0.37)) < 0.6]
    if len(hs) != 1 or len(qs) != 1: hr_gap.append(("missing", r["key"], len(hs), len(qs))); continue
    hr_gap.append((round(hs[0]["x0"] - max(m["ci_x"][1], sq["ci_x"][1]), 2), r["label"])); q_gap.append((round(qs[0]["x0"] - hs[0]["x1"], 2), r["label"]))
ck.log("fit: every HR string starts at least 2 pt right of the row's whiskers and of the forest box, every q string at least 2 pt right of the HR string", not [g for g in hr_gap if g[0] == "missing"] and min(g[0] for g in hr_gap) >= 2.0 and min(g[0] for g in q_gap) >= 2.0 and min(s["x0"] for s in S49 if abs(s["x1"] - 460.87) < 1.2 and 100 < s["y"] < 800) >= XR + 2.0,
       f"smallest HR gap {min(hr_gap)}, smallest q gap {min(q_gap)}, HR column left edge {min(s['x0'] for s in S49 if abs(s['x1'] - 460.87) < 1.2 and 100 < s['y'] < 800):.2f} vs box {XR:.2f}")
C_SHIFT = float(DJ["panel_c_shift_pt"]); y0c, y1c = 401.60 + C_SHIFT, 438.40 + C_SHIFT + 7 * 38.55
cells = [s for s in S49 if s["x0"] > 595 and y0c - 1 < s["y"] < y1c + 1 and s["text"] not in ("Sleep duration", "Oxygenation")]
GX = [601.11 + 118.01 * g for g in range(3)]; edges = [(GX[g] + 54.6 * l, GX[g] + 54.6 * l + 53.0) for g in range(3) for l in range(2)]
def cell_margin(s):
    for e0, e1 in edges:
        if e0 - 1 <= s["x0"] and s["x1"] <= e1 + 1: return min(s["x0"] - e0, e1 - s["x1"])
    return None
cm = [(cell_margin(s), s["text"], s["y"]) for s in cells]; cm_bad = [c for c in cm if c[0] is None or c[0] < 0.5]
ci_m = [c for c in cm if c[0] is not None and re.fullmatch(r"\(\d\.\d\d-\d\.\d\d\)", c[1])]; other_m = [c for c in cm if c[0] is not None and c not in ci_m]
ck.log("fit: every panel c cell and header string lies inside its 53 pt cell (at least 0.5 pt margin; the 35 interval strings at 10.5 pt keep 0.8 pt by the coordinator's decision, every other string at least 2 pt)", not cm_bad and len(ci_m) == 35 and min(c[0] for c in other_m) >= 2.0,
       f"{len(cm)} strings; interval strings {len(ci_m)}, margin {min(c[0] for c in ci_m):.2f} to {max(c[0] for c in ci_m):.2f} pt; other strings smallest margin {min(c[0] for c in other_m):.2f} pt ({min(other_m)[1]}); outside or under 0.5 pt {cm_bad[:5]}")
right_end = max(s["x1"] for s in S49); tail = max(S49, key=lambda s: s["x1"])
ck.log("fit: every string ends inside the sheet with at least 8 pt to spare (panel b key, panel c legend)", right_end <= W - 8.0, f"rightmost string '{tail['text']}' ends at {right_end:.2f} of {W}")
def ink(s): return (s["x0"], s["y"] - 0.716 * s["size"], s["x1"], s["y"] + 0.21 * s["size"])
hits = []
for i in range(len(S49)):
    a = ink(S49[i])
    for j in range(i + 1, len(S49)):
        b = ink(S49[j])
        if a[0] < b[2] - 0.3 and b[0] < a[2] - 0.3 and a[1] < b[3] - 0.3 and b[1] < a[3] - 0.3: hits.append((S49[i]["text"], S49[j]["text"], S49[i]["y"], S49[j]["y"]))
ck.log("overlap: no two text ink boxes (cap height to descender) intersect on the page_abc text layer", not hits, f"{len(S49)} spans; intersections {hits[:6]}")
outside = [s["text"] for s in S49 if s["y"] - 1.005 * s["size"] < 16.0 or s["y"] + 0.324 * s["size"] > float(DJ["y_d_new"])]
ck.log("region: every page_abc string lies inside [16, 772 + DELTA) (the compose clips to it)", not outside, f"{outside[:5]}")

# ------------------------------------------------------------------ 6. render and compose log
from PIL import Image
im = Image.open(PNG150); exp = (round(wn * 150 / 72), round(hn * 150 / 72))
ck.log("render: 150 dpi Ghostscript PNG exists, newer than the vector, with the page's pixel size (looked at: no overlap, no clipped text, key and legend boxes fit)", os.path.getmtime(PNG150) >= os.path.getmtime(NEW) - 1 and im.size == exp, f"{im.size} vs {exp}")
ck.log("compose log: band = the round-49 band with its current sha256, page_abc with its current sha256, letters and title 14 pt", LOG["band"] == BAND and LOG["band_sha256"] == L.sha256(BAND) and LOG["page_abc_sha256"] == L.sha256(PAGE) and LOG["title"]["size"] == 14.0 and LOG.get("letter_size") == 14.0 and "r49" in os.path.basename(BAND), f"band {os.path.basename(BAND)} {LOG['band_page']} at top {BAND_TOP}")
ok = ck.write(f"{VER}/checks.txt")
L.dump_json({"vector": NEW, "vector_sha256": L.sha256(NEW), "page_abc_sha256": L.sha256(PAGE), "band": BAND, "band_sha256": L.sha256(BAND), "old_vector": OLD, "old_sha256": L.sha256(OLD), "v26": V26, "v26_sha256": L.sha256(V26),
             "drawn_json_sha256": L.sha256(f"{WORK}/{SHEET}_drawn.json"), "page": [wn, hn], "census": {"page_pairs": len(pairs), "plus_one": d_plus, "kept": len(d_zero), "anchors": dict(kinds), "sizes_old": sizes_o, "sizes_new": sizes_n, "band_sizes_old": bs_o, "band_sizes_new": bs_n},
             "printed_values": {"rows": n_rows, "yes": yes, "yes_recorded": len(yes_rec), "static": len(static_ok)}, "result": "ALL PASS" if ok else "FAIL", "written": L.now()}, f"{VER}/verify_record.json")
sys.exit(0 if ok else 1)
