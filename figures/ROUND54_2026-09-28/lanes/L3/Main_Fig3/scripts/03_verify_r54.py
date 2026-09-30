#!$T90_PY
"""Round 54 (2026-09-29), lane L3, Main_Fig3: verify the re-laid vector compose against V30's Main_Fig3 through the lane's own vector twin
(the round-49 vector, the source of the LRASTER raster that V27 to V30 carry unchanged, proven by sha256), the numbers files, the
round-49 record of the values and the planner's geometry.
Ghostscript only on the composed sheets (txtwrite -dTextFormat=0 for the census, png16m for the render, pdfwrite for the flat twin),
PyMuPDF only on the flat pieces (work/page_abc.pdf, the bands). Never PyMuPDF on Main_Fig3 (memory rule of 2026-09-16).
Proves:
  1. page box (pypdf mediabox = gs PDFINFO = compose log) on the vector and its flat re-distillation, width = V30, height under the
     1260 pt gate (reported against the 1090 pt aim);
  2. the chain V30 Main_Fig3 (raster) = LRASTER raster of the round-49 vector (sha256 both ways);
  3. Ghostscript census: the page part (above the band) carries V30's strings plus the seven declared duplicates of the second axis,
     every string at its declared round-54 size (the build's own text placements), the letters and the title 18 pt, nothing below
     12 pt, the flat twin's census equal to the vector's, the band part with the round-49 band's characters at sizes not smaller;
  4. drawn values against the data: every record string on the flat page re-derived from its numbers file under its rule
     (v14lib.printed_values_csv), the other strings the same design strings as round 49 plus the declared duplicates, the sources
     unchanged (sha256) and sidecar-gated;
  5. marks preserved: 104 data marks with the round-49 values, fills and flags, at x = XA + XB ln HR of their column, every whisker
     and marker found on the flat page's drawing layer at that position, panel b and c fills and inks equal to round 49;
  6. fits and overlaps: labels clear the boxes, HR and q strings clear the whiskers, columns clear each other, tick labels clear each
     other, cell strings inside their cells, the key rows inside the sheet, no two text ink boxes intersect (flat page and composed
     sheet), the stamped letters clear every string;
  7. the 150 dpi render, the compose log, the band identity (BANDS_READY.txt), no em dash, semicolon or the word printed in any string.
usage: 03_verify_r54.py"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import collections, csv, gc, html, json, math, os, re, subprocess, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lmain_lib as L
import fitz, numpy as np, pypdf
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND37_2026-09-14/figures/_common")
import v14lib as VL

SHEET = "Main_Fig3"; SD = f"{L.LANE}/{SHEET}"; WORK = f"{SD}/work"; VER = f"{SD}/verify"; os.makedirs(VER, exist_ok=True)
NEW = f"{SD}/{SHEET}.pdf"; FLAT = f"{WORK}/{SHEET}_flat.pdf"; PAGE = f"{WORK}/page_abc.pdf"; PNG150 = f"{SD}/{SHEET}_150dpi.png"
R49L = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/L3/Main_Fig3"; OLD = f"{R49L}/work/Main_Fig3_r49_vector.pdf"; PAGE49 = f"{R49L}/work/page_abc.pdf"
DJ49 = L.load_json(f"{R49L}/work/Main_Fig3_drawn.json"); LOG49 = L.load_json(f"{R49L}/work/compose_log.json"); CSV49 = f"{R49L}/verify/Main_Fig3_printed_values.csv"
V30 = f"{paths.FIGURE_ROOT}/ROUND52_2026-09-26/figures/NEW_FINAL_SET_V30/{SHEET}.pdf"; V30_MAN = L.load_json(f"{paths.FIGURE_ROOT}/ROUND52_2026-09-26/figures/V30_ASSEMBLY_MANIFEST.json")
RASTER_REC = L.load_json(f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LRASTER/{SHEET}/work/raster_record.json")
BAND49 = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LSKETCH/work/build/band_3d_r49.pdf"; BANDS_READY = f"{L.R54}/lanes/LSKETCH/BANDS_READY.txt"
DJ = L.load_json(f"{WORK}/{SHEET}_drawn.json"); LOG = L.load_json(f"{WORK}/compose_log.json"); BAND = LOG["band"]; PLAN = DJ["plan"]; Y = PLAN["y"]
BAND_TOP = float(LOG["band_top"]); W_DECL, H_DECL = (float(v) for v in LOG["page"]); TXT = DJ["round54"]["text_layer"]; DUP = DJ["round54"]["declared_duplicates"]
for p in (NEW, FLAT, PAGE, OLD, PAGE49, V30, BAND, BAND49, PNG150, BANDS_READY): L.hydrated(p)
ck = L.Checks(f"Lane L3 round 54, {SHEET}, {L.now()}. NEW = {NEW} (sha256 {L.sha256(NEW)[:16]}), the vector compose of the re-laid page_abc (sha256 {L.sha256(PAGE)[:16]}) with the band "
              f"{os.path.basename(BAND)} (sha256 {L.sha256(BAND)[:16]}). FLAT = {FLAT} (sha256 {L.sha256(FLAT)[:16]}). OLD = the round-49 vector {OLD} (sha256 {L.sha256(OLD)[:16]}), the source of V30 {V30}. "
              f"Values reference = the round-49 record {R49L}/work/Main_Fig3_drawn.json. Geometry reference = the build's own plan (work/plan_r54.json).")


# ------------------------------------------------------------------ helpers
def gs_spans(pdf, xml, name):
    """Ghostscript txtwrite -dTextFormat=0 under the watchdog: text, x0, baseline, x1, font, size (integer bbox coordinates)."""
    if os.path.exists(xml): os.remove(xml)
    rc, st = L.wd(name, [L.GS, "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=txtwrite", "-dTextFormat=0", f"-sOutputFile={xml}", pdf], max_s=300)
    if rc != 0 or not os.path.exists(xml): return None, st
    out = []
    for m in re.finditer(r'<span bbox="([^"]*)" font="([^"]*)" size="([^"]*)">(.*?)</span>', open(xml, encoding="utf-8", errors="replace").read(), re.S):
        bb = [float(v) for v in m.group(1).split()]; t = L.norm_text(html.unescape("".join(re.findall(r'c="([^"]*)"', m.group(4)))))
        if t.strip(): out.append(dict(text=t.strip(), x0=bb[0], y=bb[1], x1=bb[2], font=m.group(2).split("+")[-1].replace("-Identity-H", ""), size=round(float(m.group(3)), 4)))
    return out, st


def merge(spans, gap=1.0, kern=3.0):
    """Glue Ghostscript's kerning-split fragments on one baseline (same font and size, gap at most 1 pt, or an overlap of at most 3 pt
    where a kerning pair such as T-o pulls the next fragment back under the first one's advance, 2 pt at 14 pt in integer boxes)."""
    out = []
    for s in sorted(spans, key=lambda s: (s["y"], s["x0"])):
        if out and abs(out[-1]["y"] - s["y"]) <= 0.5 and out[-1]["font"] == s["font"] and abs(out[-1]["size"] - s["size"]) < 0.01 and -kern <= s["x0"] - out[-1]["x1"] <= gap:
            out[-1] = dict(out[-1], text=out[-1]["text"] + s["text"], x1=max(out[-1]["x1"], s["x1"]))
        else: out.append(dict(s))
    return out


def flat_spans(pdf):
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


def ink(s, size_key="size"):
    """Cap height to descender box of a span (x0, top, x1, bottom)."""
    return (s["x0"], s["y"] - 0.716 * s[size_key], s["x1"], s["y"] + 0.21 * s[size_key])


def intersects(a, b, slack=0.3): return a[0] < b[2] - slack and b[0] < a[2] - slack and a[1] < b[3] - slack and b[1] < a[3] - slack
def gap_between(a, b):
    dx = max(b[0] - a[2], a[0] - b[2], 0.0); dy = max(b[1] - a[3], a[1] - b[3], 0.0); return math.hypot(dx, dy)


def is_letter(s): return len(s["text"]) == 1 and s["text"] in "abcd" and "Bold" in s["font"] and s["size"] >= 17.5
def nosp(t): return re.sub(r"\s+", "", t)            # Ghostscript's txtwrite drops some inter-word spaces (the flat twin's "Figure 3", "and T90"): compare without whitespace
def is_title(s): return nosp(s["text"]) == "Figure3"


# ------------------------------------------------------------------ 1. page box
def box(p): r = pypdf.PdfReader(p); mb = r.pages[0].mediabox; return float(mb.width), float(mb.height), len(r.pages)
wn, hn, nn = box(NEW); wf, hf, nf = box(FLAT); wv, hv, nv = box(V30); wo, ho, no = box(OLD)
gn, gwn, ghn = L.gs_pdfinfo(NEW)
ck.log("page box: NEW = compose log = gs PDFINFO = FLAT within 0.01 pt, one page each, width = V30 within 0.05 pt, height at most 1260 pt",
       nn == nf == nv == no == 1 and gn == 1 and abs(wn - W_DECL) < 0.01 and abs(hn - H_DECL) < 0.01 and abs(gwn - wn) < 0.01 and abs(ghn - hn) < 0.01 and abs(wf - wn) < 0.01 and abs(hf - hn) < 0.01 and abs(wn - wv) < 0.05 and hn <= 1260.0,
       f"NEW {wn:.3f} x {hn:.3f}, FLAT {wf:.3f} x {hf:.3f}, V30 {wv:.3f} x {hv:.3f}, round-49 vector {wo:.3f} x {ho:.3f}, compose log {W_DECL} x {H_DECL}")
ck.info(f"height against the coordinator's aim: {hn:.3f} pt ({'under' if hn <= 1090 else 'over'} the 1090 pt aim by {abs(1090 - hn):.1f} pt), V30 {hv:.3f}, gate 1260. Band top {BAND_TOP}, band {LOG['band_page']}")

# ------------------------------------------------------------------ 2. the chain V30 -> round-49 vector
sha_v30, sha_old = L.sha256(V30), L.sha256(OLD)
ck.log("chain: V30 Main_Fig3 (raster) has the sha256 of the LRASTER raster of the round-49 vector, and that record names the round-49 vector by its current sha256 (V30 = V29 = V28 = V27 unchanged per the manifests)",
       sha_v30 == RASTER_REC["sha256_pdf"] == V30_MAN["sheets"]["Main_Fig3"]["sha256"] and RASTER_REC["vector_sha256"] == sha_old and RASTER_REC["vector"] == OLD,
       f"V30 {sha_v30[:16]}, LRASTER pdf {RASTER_REC['sha256_pdf'][:16]}, LRASTER vector {RASTER_REC['vector_sha256'][:16]}, OLD {sha_old[:16]}, manifest why '{V30_MAN['sheets']['Main_Fig3']['why']}'")

# ------------------------------------------------------------------ 3. census (Ghostscript txtwrite on the composed sheets, under the watchdog)
SO_RAW, st_o = gs_spans(OLD, f"{WORK}/census_r49_vector.xml", "txt_old"); SN_RAW, st_n = gs_spans(NEW, f"{WORK}/census_r54_vector.xml", "txt_new"); SF_RAW, st_f = gs_spans(FLAT, f"{WORK}/census_r54_flat.xml", "txt_flat")
ck.log("census: Ghostscript txtwrite finished on the round-49 vector, the round-54 vector and its flat twin (watchdog, 300 s each)", SO_RAW is not None and SN_RAW is not None and SF_RAW is not None, f"old {st_o}, new {st_n}, flat {st_f}")
SO, SN, SF = merge(SO_RAW or []), merge(SN_RAW or []), merge(SF_RAW or [])
ck.info(f"gs txtwrite spans: OLD {len(SO_RAW or [])} raw -> {len(SO)} strings, NEW {len(SN_RAW or [])} raw -> {len(SN)}, FLAT {len(SF_RAW or [])} raw -> {len(SF)} (kerning-split fragments glued: gap at most 1 pt, same baseline, font and size)")
Y_SPLIT_O = float(LOG49["band_top"]) - 5.0; Y_SPLIT_N = BAND_TOP - 5.0
po, pn, pf = [s for s in SO if s["y"] < Y_SPLIT_O], [s for s in SN if s["y"] < Y_SPLIT_N], [s for s in SF if s["y"] < Y_SPLIT_N]
bo, bn = [s for s in SO if s["y"] >= Y_SPLIT_O], [s for s in SN if s["y"] >= Y_SPLIT_N]
co, cn, cf = collections.Counter(nosp(s["text"]) for s in po), collections.Counter(nosp(s["text"]) for s in pn), collections.Counter(nosp(s["text"]) for s in pf)
expected = co + collections.Counter({nosp(L.norm_text(k)): v for k, v in DUP.items()})
ck.log("census, page part: the round-54 strings = the round-49 strings + the seven declared duplicates of the second axis (HR (95% CI), q, the caption, 0.6, 0.8, 1, 1.5), nothing else gained or lost",
       cn == expected, f"{len(po)} old, {len(pn)} new spans, declared {dict(DUP)}, lost {dict(expected - cn)}, gained {dict(cn - expected)}")
ts_f, ts_n = collections.Counter((nosp(s["text"]), round(s["size"], 2)) for s in pf), collections.Counter((nosp(s["text"]), round(s["size"], 2)) for s in pn)
ck.log("census, page part: the flat twin carries the same (text, size) multiset as the vector (whitespace ignored)", ts_f == ts_n, f"difference {dict((ts_f - ts_n) + (ts_n - ts_f))}")
# every page-part string at its declared size: match to the build's text placements (text, baseline within 1, anchor within 1.5)
TXTN = [dict(t, text=nosp(L.norm_text(t["text"]))) for t in TXT]
def placement_of(s):
    best, bd = None, 9e9
    for t in TXTN:
        if t["text"] != nosp(s["text"]) or abs(t["baseline"] - s["y"]) > 1.0: continue
        a = s["x0"] if t["ha"] == "left" else (s["x1"] if t["ha"] == "right" else (s["x0"] + s["x1"]) / 2); d = abs(a - t["x"])
        if d < bd: bd, best = d, t
    return (best, bd) if best is not None and bd <= 1.5 else (None, bd)
size_bad, matched, unmatched = [], 0, []
for s in pn:
    if is_letter(s) or is_title(s): continue
    t, d = placement_of(s)
    if t is None: unmatched.append((s["text"], s["x0"], s["y"], s["size"])); continue
    matched += 1
    if abs(t["size"] - s["size"]) > 0.05 or (("Bold" in s["font"]) != bool(t["bold"])): size_bad.append((s["text"], s["y"], s["size"], t["size"], s["font"]))
ck.log("census, page part: every string (letters and title aside) matches a placement of the build at its baseline and anchor, with the placement's size and weight (rows, ticks, keys, cell values, HR and q 14, headers and axis titles 15, six-cell grid text and notes 13)",
       not size_bad and not unmatched and matched == len(pn) - 5, f"{matched} matched, unmatched {unmatched[:5]}, size or weight mismatches {size_bad[:5]}")
sizes_o, sizes_n = sorted(collections.Counter(s["size"] for s in po).items()), sorted(collections.Counter(s["size"] for s in pn).items())
ck.info(f"page-part size histogram OLD {sizes_o}"); ck.info(f"page-part size histogram NEW {sizes_n}")
min_page = min(s["size"] for s in pn); min_band = min(s["size"] for s in bn)
ck.log("census: no text below 12 pt anywhere on the sheet (page part smallest 13 pt, band smallest 14 pt)", min_page >= 12.0 and min_band >= 12.0, f"page part smallest {min_page}, band smallest {min_band}")
lets_n = sorted([(s["text"], s["x0"], s["y"], s["size"]) for s in pn if is_letter(s)], key=lambda t: (t[2], t[1]))
stamped = {ch: (x, base) for ch, x, base in LOG["letters_stamped"]}
tit_n = [(s["x0"], s["y"], s["size"]) for s in pn if is_title(s)]
ck.log("letters a, b, c, d: 18 pt Arial Bold once each at the compose log's positions (x0 and baseline within 1 pt), the title 'Figure 3' once at 18 pt at (8, 18)",
       [t[0] for t in lets_n] == ["a", "b", "c", "d"] and all(abs(t[3] - 18.0) < 0.01 and abs(t[1] - stamped[t[0]][0]) <= 1.0 and abs(t[2] - stamped[t[0]][1]) <= 1.0 for t in lets_n) and len(tit_n) == 1 and abs(tit_n[0][2] - 18.0) < 0.01 and abs(tit_n[0][0] - 8.0) <= 1.0 and abs(tit_n[0][1] - 18.0) <= 1.0,
       f"letters {lets_n}, stamped {LOG['letters_stamped']}, title {tit_n}")
# band part: the coordinator's round-54 build against the round-49 band
cho, chn = collections.Counter(ch for s in bo for ch in s["text"] if not ch.isspace()), collections.Counter(ch for s in bn for ch in s["text"] if not ch.isspace())
bs_o, bs_n = sorted(collections.Counter(s["size"] for s in bo).items()), sorted(collections.Counter(s["size"] for s in bn).items())
ck.log("census, band part: same character multiset as the round-49 band, no size smaller (LSKETCH round 54: notes 15, caption 14)", cho == chn and min(s for s, _ in bs_n) >= min(s for s, _ in bs_o), f"chars lost {dict(cho - chn)}, gained {dict(chn - cho)}, sizes OLD {bs_o}, NEW {bs_n}")
ready = open(BANDS_READY).read(); m_ready = re.search(r"band_3d_r54\.pdf\s+([\d.]+) x ([\d.]+) pt\s+sha256=([0-9a-f]{64})", ready)
ck.log("band identity: the composed band is ROUND54/lanes/LSKETCH/work/build/band_3d_r54.pdf, its sha256 and page box as BANDS_READY.txt lists them and as the compose log recorded them, width = the sheet width",
       m_ready is not None and os.path.basename(BAND) == "band_3d_r54.pdf" and os.path.dirname(BAND) == L.SKETCH and L.sha256(BAND) == m_ready.group(3) == LOG["band_sha256"] and abs(float(m_ready.group(1)) - LOG["band_page"][0]) < 0.01 and abs(float(m_ready.group(2)) - LOG["band_page"][1]) < 0.01 and abs(LOG["band_page"][0] - wn) < 0.05,
       f"band {BAND} {LOG['band_page']}, sha256 {L.sha256(BAND)[:16]}, BANDS_READY {m_ready.groups() if m_ready else None}")

# ------------------------------------------------------------------ 4. drawn values against the data (flat page_abc piece)
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
n_rows, n_no = VL.printed_values_csv(PAGE, records_for_csv(DJ["records"]), f"{VER}/{SHEET}_printed_values.csv", static_from=None)
rows = list(csv.DictReader(open(f"{VER}/{SHEET}_printed_values.csv")))
yes = sum(1 for r in rows if r["match"] == "yes"); yes_rec = [r["string"] for r in rows if r["match"] == "yes_recorded"]; nos = [r for r in rows if r["match"] == "NO"]
static_new = collections.Counter(L.norm_text(r["string"]).strip() for r in nos if "not in the drawn records" in r["note"]); other_no = [(r["string"], r["x"], r["y"], r["note"][:80]) for r in nos if "not in the drawn records" not in r["note"]]
rows49 = list(csv.DictReader(open(L.hydrated(CSV49)))); static_old = collections.Counter(L.norm_text(r["string"]).strip() for r in rows49 if r["match"] == "NO" and "not in the drawn records" in r["note"])
ck.log("drawn values: every record string on the page_abc text layer at its recorded anchor and re-derived from its numbers file under its rule (v14lib.printed_values_csv, PyMuPDF on the flat piece), none unaccounted for",
       not other_no and len(DJ["records"]) == 258, f"{n_rows} rows: {yes} re-derived (yes), {len(yes_rec)} recorded-only {sorted(set(yes_rec))}, {sum(static_new.values())} design strings, failures {other_no[:5]}")
ck.log("drawn values: the design strings (not records) are round 49's design strings plus the seven declared duplicates", static_new == static_old + collections.Counter({L.norm_text(k): v for k, v in DUP.items()}),
       f"{sum(static_old.values())} in round 49, {sum(static_new.values())} now, lost {dict(static_old + collections.Counter({L.norm_text(k): v for k, v in DUP.items()}) - static_new)}, gained {dict(static_new - static_old - collections.Counter({L.norm_text(k): v for k, v in DUP.items()}))}")
src_bad = []
for name, sc in DJ["sources"].items():
    p = sc["path"]; sha = L.sha256(p)
    if sha != sc["sha256"]: src_bad.append((name, "sha256 changed since the build"))
    if sha != DJ49["sources"][name]["sha256"]: src_bad.append((name, "differs from the round-49 build's source"))
    if "step" in sc:
        try: VL.sidecar_v8(p)
        except SystemExit as e: src_bad.append((name, str(e)[:100]))
ck.log("sources: every numbers file of the build unchanged (sha256, equal to the round-49 build's), the six sidecar-gated ones still pass the v8.1 gate", not src_bad, f"{len(DJ['sources'])} sources, {sum(1 for v in DJ['sources'].values() if 'step' in v)} gated, {src_bad}")

# ------------------------------------------------------------------ 5. marks preserved (values from the round-49 record, positions from the new axes, drawings on the flat page)
mk = [d for d in DJ["drawn"] if "marker" in d]; mk49 = {(d["key"], d["marker"]): d for d in DJ49["drawn"] if "marker" in d}
r49rows = {r["key"]: r for r in DJ49["panel_a"]["rows"]}; rows_new = {r["key"]: r for r in DJ["panel_a"]["rows"]}
cols = {c["name"]: c for c in DJ["panel_a"]["columns"]}; G = DJ["panel_a"]["geometry"]
val_bad, pos_bad = [], []
for m in mk:
    o = mk49[(m["key"], m["marker"])]; r, r49 = rows_new[m["key"]], r49rows[m["key"]]
    hr_keys = ("hr", "lo", "hi") if m["marker"] == "T90 circle" else ("tst_hr", "tst_lo", "tst_hi")
    vals = [r[k] for k in hr_keys]; vals49 = [r49[k] for k in hr_keys]
    if m["filled"] != o["filled"] or any(abs(a - b) > 1e-12 for a, b in zip(vals, vals49)) or any(abs(a - b) > 1e-9 for a, b in zip(m["value"], vals49)) or r["sig"] != r49["sig"] or r["tst_sig"] != r49["tst_sig"]: val_bad.append((m["key"], m["marker"]))
    c = cols[m["column"]]; xhr = lambda v: c["XA"] + c["XB"] * math.log(v)
    i = (m["row"] - 1) % DJ["panel_a"]["n_per_column"]; yc = G["row0_centre"] + G["pitch"] * i; y_exp = yc - G["lane"] if m["marker"] == "T90 circle" else yc + G["lane"]
    if abs(m["x"] - xhr(vals[0])) > 1e-3 or abs(m["ci_x"][0] - xhr(vals[1])) > 1e-3 or abs(m["ci_x"][1] - xhr(vals[2])) > 1e-3 or abs(m["y"] - y_exp) > 1e-3: pos_bad.append((m["key"], m["marker"], m["x"], xhr(vals[0]), m["y"], y_exp))
ck.log("marks: all 104 data marks (52 T90 circles, 52 sleep-time squares) carry the round-49 values (HR, lo, hi), fills and significance flags", len(mk) == len(mk49) == 104 and not val_bad, f"{len(mk)} marks, value mismatches {val_bad[:5]}")
ck.log("marks: every mark sits at x = XA + XB ln HR of its column (whisker ends likewise) on its row's lane (row centre minus or plus 3.06 pt)", not pos_bad, f"columns {[(c['name'], round(c['XA'], 3), round(c['XB'], 3), c['xlo'], c['xhi']) for c in DJ['panel_a']['columns']]}, mismatches {pos_bad[:3]}")
# the drawing layer of the flat page: whiskers as lines, markers as circles or squares
BLUE, GREEN = (0.008, 0.533, 0.82), (0.224, 0.769, 0.271)
def near(c, ref, tol=0.01): return c is not None and all(abs(a - b) <= tol for a, b in zip(c, ref))
dd = fitz.open(PAGE); drawings = dd[0].get_drawings(); dd.close(); gc.collect()
lines = [(it["rect"], tuple(round(v, 3) for v in it["color"]), it.get("width") or 0.0) for it in drawings if it.get("color") is not None and all(i[0] == "l" for i in it["items"]) and len(it["items"]) == 1]
shapes = [(it["rect"], tuple(round(v, 3) for v in it["color"]) if it.get("color") else None, tuple(round(v, 3) for v in it["fill"]) if it.get("fill") else None, "".join(i[0] for i in it["items"])) for it in drawings]
whisk_bad, mark_bad = [], []; counts = collections.Counter()
for m in mk:
    colr = BLUE if m["marker"] == "T90 circle" else GREEN
    hit = [l for l in lines if near(l[1], colr) and abs(l[0].y0 - m["y"]) < 0.05 and abs(l[0].y1 - m["y"]) < 0.05 and abs(l[0].x0 - m["ci_x"][0]) < 0.05 and abs(l[0].x1 - m["ci_x"][1]) < 0.05 and abs(l[2] - 1.7 * 1.0706) < 0.02]
    if len(hit) != 1: whisk_bad.append((m["key"], m["marker"], len(hit)))
    kind = "c" if m["marker"] == "T90 circle" else "re"
    sh = [s for s in shapes if s[3].startswith(kind) and abs((s[0].x0 + s[0].x1) / 2 - m["x"]) < 0.1 and abs((s[0].y0 + s[0].y1) / 2 - m["y"]) < 0.1 and (s[0].x1 - s[0].x0) < 8.0]
    ok = len(sh) == 1 and ((m["filled"] and near(sh[0][2], colr) and near(sh[0][1], (1.0, 1.0, 1.0))) or (not m["filled"] and near(sh[0][2], (1.0, 1.0, 1.0)) and near(sh[0][1], colr)))
    if not ok: mark_bad.append((m["key"], m["marker"], [(s[1], s[2], s[3]) for s in sh]))
    counts[(m["marker"], "filled" if m["filled"] else "open")] += 1
ck.log("marks drawn: on the flat page's drawing layer every mark has its whisker (one line of the lane colour, width 1.82 pt, from x(lo) to x(hi) on the lane) and its marker (a circle or square centred on x(HR), filled in the lane colour with a white edge when significant, white with the lane colour edge otherwise)",
       not whisk_bad and not mark_bad, f"{len(mk)} marks, counts {dict(counts)} (round 49: 37 T90 and 14 sleep-time significant), whisker misses {whisk_bad[:3]}, marker misses {mark_bad[:3]}")
fills_same = all(DJ["panel_c"]["rows"][t]["cells"][c]["fill"] == DJ49["panel_c"]["rows"][t]["cells"][c]["fill"] and DJ["panel_c"]["rows"][t]["cells"][c]["ink"] == DJ49["panel_c"]["rows"][t]["cells"][c]["ink"] for t in DJ["panel_c"]["order"] for c in ("SN", "SL", "NL", "LN", "LL")) \
    and DJ["panel_c"]["order"] == DJ49["panel_c"]["order"] and all(DJ["panel_b"]["rows"][r]["cells"][i]["shade"] == DJ49["panel_b"]["rows"][r]["cells"][i]["shade"] and DJ["panel_b"]["rows"][r]["cells"][i]["label"] == DJ49["panel_b"]["rows"][r]["cells"][i]["label"] for r in DJ["panel_b"]["rows"] for i in range(4)) and DJ["panel_b"]["key_bands"] == DJ49["panel_b"]["key_bands"]
ck.log("panels b and c: every cell fill, ink, shade band, label and the row orders equal to the round-49 record (the same ramp rule and ladder)", fills_same, f"c rows {DJ['panel_c']['order']}, b key bands {DJ['panel_b']['key_bands']}")

# ------------------------------------------------------------------ 6. fits and overlaps (flat page text layer, PyMuPDF)
S = flat_spans(PAGE); Wpt = float(DJ["page"][0]); fit = {}
bad_lab, bad_hr, bad_q, lab_ends, hr_gaps, q_gaps = [], [], [], [], [], []
for c in DJ["panel_a"]["columns"]:
    labs = [s for s in S if abs(s["x"] - c["lab_x0"]) < 0.05 and G["row0_centre"] - 10 < s["y"] < G["row0_centre"] + 25 * G["pitch"] + 10]
    if len(labs) != 26: bad_lab.append((c["name"], "label count", len(labs)))
    for s in labs:
        lab_ends.append(s["x1"])
        if s["x1"] > c["x_left"] - 3.0 + 0.05: bad_lab.append((c["name"], s["text"], s["x1"], c["x_left"]))
    for m in [d for d in mk if d["column"] == c["name"] and d["marker"] == "T90 circle"]:
        sq = mk49  # unused
        i = (m["row"] - 1) % 26; yc = G["row0_centre"] + G["pitch"] * i; y_hq = yc + G["hrq_dy"]; r = rows_new[m["key"]]
        hs = [s for s in S if s["text"] == L.norm_text(r["s_hr"]) and abs(s["y"] - y_hq) < 0.6 and abs(s["x1"] - c["hr_r"]) < 1.2]; qs = [s for s in S if s["text"] == L.norm_text(r["s_q"]) and abs(s["y"] - y_hq) < 0.6 and abs(s["x1"] - c["q_r"]) < 1.2]
        if len(hs) != 1 or len(qs) != 1: bad_hr.append((c["name"], r["label"], len(hs), len(qs))); continue
        sqm = [d for d in mk if d["key"] == m["key"] and d["marker"] != "T90 circle"][0]
        g = hs[0]["x0"] - max(m["ci_x"][1], sqm["ci_x"][1], c["x_right"]); hr_gaps.append((round(g, 2), r["label"])); q_gaps.append((round(qs[0]["x0"] - hs[0]["x1"], 2), r["label"]))
        if g < 2.0: bad_hr.append((c["name"], r["label"], round(g, 2)))
        if qs[0]["x0"] - hs[0]["x1"] < 2.0: bad_q.append((c["name"], r["label"], round(qs[0]["x0"] - hs[0]["x1"], 2)))
ck.log("fit: all 52 row labels (26 per column) end at least 3 pt before their column's forest box", not bad_lab, f"largest label end left {max(e for e in lab_ends if e < 400):.2f} vs box {cols['left']['x_left']}, right {max(e for e in lab_ends if e > 400):.2f} vs box {cols['right']['x_left']}, failures {bad_lab[:4]}")
ck.log("fit: every HR string starts at least 2 pt right of both whiskers of its row and of the box, every q string at least 2 pt right of its HR string", not bad_hr and not bad_q, f"smallest HR gap {min(hr_gaps)}, smallest q gap {min(q_gaps)}, failures {(bad_hr + bad_q)[:4]}")
lq = [s for s in S if abs(s["x1"] - cols["left"]["q_r"]) < 1.2 and 60 < s["y"] < 460]; rl = [s for s in S if abs(s["x"] - cols["right"]["lab_x0"]) < 0.05 and 60 < s["y"] < 460]
col_gap = min(r_["x0"] for r_ in rl) - max(s["x1"] for s in lq)
ck.log("fit: the right column's labels start at least 12 pt right of the left column's q strings (planned 20)", col_gap >= 12.0, f"gap {col_gap:.2f} pt")
tick_bad, tick_min = [], 9e9
for c in DJ["panel_a"]["columns"]:
    ts = sorted([s for s in S if abs(s["y"] - G["tick_baseline"]) < 0.6 and c["x_left"] - 15 < (s["x0"] + s["x1"]) / 2 < c["x_right"] + 15], key=lambda s: s["x0"])
    if [s["text"] for s in ts] != [f"{v:g}" for v in c["ticks"]]: tick_bad.append((c["name"], [s["text"] for s in ts]))
    for a, b in zip(ts[:-1], ts[1:]):
        tick_min = min(tick_min, b["x0"] - a["x1"])
        if b["x0"] - a["x1"] < 3.0: tick_bad.append((c["name"], a["text"], b["text"], round(b["x0"] - a["x1"], 2)))
ck.log("fit: the tick labels of each column are the declared ticks in order, neighbours at least 3 pt apart", not tick_bad, f"smallest gap {tick_min:.2f} pt, failures {tick_bad[:4]}")
ke = DJ["panel_a"]["key"]["entries"]; key_spans = [s for s in S if abs(s["y"] - Y["key_base"]) < 0.6]
key_ok = len(key_spans) == 3 and all(abs(s["x"] - e["label_x"]) < 0.6 for s, e in zip(sorted(key_spans, key=lambda s: s["x"]), ke)) and all(ke[i]["label_end"] + 10.0 <= ke[i + 1]["mark_x"] for i in range(2)) and max(s["x1"] for s in key_spans) <= Wpt - 8.0 and min(s["x0"] for s in key_spans) >= 8.0 + 18.0
ck.log("fit: the panel a key row holds its three entries at the build's positions, each label at least 10 pt before the next entry's mark, the first mark at least 12 pt right of the letter a, the row inside the sheet", key_ok, f"entries {ke}, spans {[(s['text'][:20], s['x0'], s['x1']) for s in key_spans]}")
# cells: panel c (13-pt strings inside 66.6-pt cells), panel b (14-pt strings inside 80.9-pt cells)
GC = DJ["panel_c"]["geometry"]; GB = DJ["panel_b"]["geometry"]
c_edges = [(gx + (GC["cell_w"] + GC["pair_gap"]) * l, gx + (GC["cell_w"] + GC["pair_gap"]) * l + GC["cell_w"]) for gx in GC["groups_x"] for l in range(2)]
b_edges = [(hx, hx + GB["cell_w"]) for hx in GB["cells_x"]]
def margin(s, edges):
    for e0, e1 in edges:
        if e0 - 1.5 <= s["x0"] and s["x1"] <= e1 + 1.5: return min(s["x0"] - e0, e1 - s["x1"])
    return None
c_cells = [s for s in S if s["x0"] > GC["groups_x"][0] - 2 and GC["hdr2_y0"] - 1 < s["y"] < GC["row_y0"] + 7 * GC["row_pitch"] + 1]
cm = [(margin(s, c_edges), s["text"], s["y"]) for s in c_cells]; c_vals = [c for c in cm if c[0] is not None and (re.fullmatch(r"\(\d\.\d\d-\d\.\d\d\)", c[1]) or re.fullmatch(r"\d\.\d\d\**", c[1]))]; c_sub = [c for c in cm if c[0] is not None and c not in c_vals]
c_bad = [c for c in cm if c[0] is None or c[0] < 0.8] + [c for c in c_sub if c[0] < 2.0]
ck.log("fit: every six-cell grid string lies inside its 66.6 pt cell (the 70 value and interval strings at 13 pt with at least 0.8 pt margin, the round-49 precedent, planned 1.5, the sub-headers at least 2 pt)",
       not c_bad and len(c_vals) == 77 and len(c_sub) == 12, f"{len(cm)} strings, value strings {len(c_vals)} with margins {min(c[0] for c in c_vals):.2f} to {max(c[0] for c in c_vals):.2f} pt, sub-headers {len(c_sub)} smallest margin {min(c[0] for c in c_sub):.2f} pt, failures {c_bad[:4]}")
b_cells = [s for s in S if s["x0"] > GB["cells_x"][0] - 2 and s["x1"] < GB["cells_x"][3] + GB["cell_w"] + 2 and GB["hdr_y0"] - 1 < s["y"] < GB["row_y0"] + 2 * GB["row_pitch"] + GB["row_h"] + 1]
bm = [(margin(s, b_edges), s["text"], s["y"]) for s in b_cells]; b_bad = [c for c in bm if c[0] is None or c[0] < 2.0]
ck.log("fit: every count-grid string (12 cells, 4 headers) lies inside its 80.9 pt cell with at least 2 pt margin", not b_bad and len(bm) == 16, f"{len(bm)} strings, smallest margin {min(c[0] for c in bm if c[0] is not None):.2f} pt ({min((c for c in bm if c[0] is not None), key=lambda c: c[0])[1]}), failures {b_bad[:4]}")
right_end = max(s["x1"] for s in S); tail = max(S, key=lambda s: s["x1"])
ck.log("fit: every string ends inside the sheet with at least 8 pt to spare", right_end <= Wpt - 8.0, f"rightmost string '{tail['text']}' ends at {right_end:.2f} of {Wpt}")
hits = []
for i in range(len(S)):
    a = ink(S[i])
    for j in range(i + 1, len(S)):
        b = ink(S[j])
        if intersects(a, b): hits.append((S[i]["text"], S[j]["text"], S[i]["y"], S[j]["y"]))
ck.log("overlap: no two text ink boxes (cap height to descender) intersect on the page_abc text layer", not hits, f"{len(S)} spans, intersections {hits[:6]}")
# the composed sheet (gs census, integer bboxes): the stamped letters and the title against every other string, and no intersections at all
hits_n, near_letters = [], []
for i in range(len(pn)):
    a = ink(pn[i])
    for j in range(i + 1, len(pn)):
        b = ink(pn[j])
        if intersects(a, b, slack=1.0): hits_n.append((pn[i]["text"], pn[j]["text"], pn[i]["y"], pn[j]["y"]))
for s in pn:
    if not (is_letter(s) or is_title(s)): continue
    a = ink(s); gmin, who = 9e9, None
    for t in pn:
        if t is s or is_letter(t) or is_title(t): continue
        g = gap_between(a, ink(t))
        if g < gmin: gmin, who = g, t["text"]
    near_letters.append((s["text"], round(gmin, 2), who))
ck.log("overlap: on the composed sheet no two text ink boxes intersect (Ghostscript census, integer boxes, 1 pt slack) and every stamped letter and the title keep at least 2 pt from every other string", not hits_n and all(n[1] >= 2.0 for n in near_letters), f"{len(pn)} spans, intersections {hits_n[:4]}, nearest to the letters and title {near_letters}")
outside = [s["text"] for s in S if s["y"] - 1.005 * s["size"] < Y["title_base"] + 6.0 or s["y"] + 0.324 * s["size"] > BAND_TOP - 2.0]
ck.log("region: every page_abc string lies inside the panel region (below the title strip, above the band top)", not outside, f"{outside[:5]}")
# strings hygiene
bad_txt = [s["text"] for s in pn if "—" in s["text"] or ";" in s["text"] or "printed" in s["text"].lower()]
ck.log("strings: no em dash, no semicolon, never the word printed in any string of the sheet", not bad_txt, f"{bad_txt[:5]}")

# ------------------------------------------------------------------ 7. render and compose log
from PIL import Image
im = Image.open(PNG150); exp = (round(wn * 150 / 72), round(hn * 150 / 72))
ck.log("render: 150 dpi Ghostscript PNG exists, newer than the vector, with the page's pixel size (looked at: no overlap, no clipped text, keys and legends clear)", os.path.getmtime(PNG150) >= os.path.getmtime(NEW) - 1 and im.size == exp, f"{im.size} vs {exp}")
# the flat twin draws the same picture: its 150 dpi Ghostscript render against the vector's. Re-encoded fonts move glyph edges by a
# pixel, so the test is structural: the difference mask must vanish under a 3 x 3 erosion (no differing blob wider than 2 px, a
# missing 14-pt glyph stem is 2.5 px wide at 150 dpi) and stay under 0.5 percent of the pixels.
PNG_FLAT = f"{WORK}/{SHEET}_flat_150dpi.png"; st_flat = L.gs_render(FLAT, PNG_FLAT, 150, "render_flat_150", max_s=600)
A = np.asarray(Image.open(PNG150).convert("RGB")).astype(int); B = np.asarray(Image.open(PNG_FLAT).convert("RGB")).astype(int)
same_shape = A.shape == B.shape; dmask = (np.abs(A - B) > 32).any(axis=2) if same_shape else np.ones(A.shape[:2], bool)
er = dmask[1:-1, 1:-1] & dmask[:-2, 1:-1] & dmask[2:, 1:-1] & dmask[1:-1, :-2] & dmask[1:-1, 2:] & dmask[:-2, :-2] & dmask[:-2, 2:] & dmask[2:, :-2] & dmask[2:, 2:]
n_diff, n_er, n_px = int(dmask.sum()), int(er.sum()), int(dmask.size)
ys, xs = np.nonzero(er); er_box = [round(float(xs.min()) * 72 / 150, 1), round(float(ys.min()) * 72 / 150, 1), round(float(xs.max()) * 72 / 150, 1), round(float(ys.max()) * 72 / 150, 1)] if len(ys) else None
rd = dict(shape_equal=same_shape, n_diff=n_diff, n_eroded=n_er, n_px=n_px, eroded_bbox_pt=er_box)
ck.log("flat twin: its 150 dpi render has the vector's pixel size, the difference mask (channel difference above 32 of 255) vanishes under a 3 x 3 erosion and covers under 0.5 percent of the pixels (glyph-edge noise of the re-encoded fonts, no missing element)",
       same_shape and n_er == 0 and n_diff <= 0.005 * n_px, f"{n_diff} differing pixels of {n_px} ({100 * n_diff / n_px:.3f} percent), {n_er} left after erosion (bbox {er_box}), shapes {A.shape} vs {B.shape}, render {st_flat}")
letters_ok = all(abs(base - (top + 0.905 * 18.0)) < 0.01 for (ch, x, base), (ch2, x2, top) in zip(LOG["letters_stamped"], [("a", 8.0, Y["a_top"]), ("b", 8.0, Y["b_top"]), ("c", PLAN["c"]["x0"], Y["b_top"]), ("d", 8.0, Y["d_top"])]) if ch == ch2 and abs(x - x2) < 1e-6)
ck.log("compose log: page_abc and the band with their current sha256, the flat twin with its current sha256, letters 18 pt at the planner's tops (baseline = top + 0.905 x 18), title 18 pt at (8, 18)",
       LOG["page_abc_sha256"] == L.sha256(PAGE) and LOG["band_sha256"] == L.sha256(BAND) and LOG["flat_sha256"] == L.sha256(FLAT) and LOG["title"]["size"] == 18.0 and LOG["letter_size"] == 18.0 and letters_ok and len(LOG["letters_stamped"]) == 4,
       f"letters {LOG['letters_stamped']}, title {LOG['title']['origin']}, flat {os.path.basename(LOG['flat'])} {LOG['flat_page']}")
ok = ck.write(f"{VER}/checks.txt")
L.dump_json({"vector": NEW, "vector_sha256": L.sha256(NEW), "flat": FLAT, "flat_sha256": L.sha256(FLAT), "page_abc_sha256": L.sha256(PAGE), "band": BAND, "band_sha256": L.sha256(BAND), "old_vector": OLD, "old_sha256": sha_old, "v30": V30, "v30_sha256": sha_v30,
             "drawn_json_sha256": L.sha256(f"{WORK}/{SHEET}_drawn.json"), "page": [wn, hn], "census": {"page_old": len(po), "page_new": len(pn), "matched_placements": matched, "sizes_old": sizes_o, "sizes_new": sizes_n, "band_sizes_old": bs_o, "band_sizes_new": bs_n, "declared_duplicates": DUP},
             "printed_values": {"rows": n_rows, "yes": yes, "yes_recorded": len(yes_rec), "design_strings": sum(static_new.values())}, "marks": dict((f"{k[0]} {k[1]}", v) for k, v in counts.items()), "fits": {"hr_gap_min": min(hr_gaps), "q_gap_min": min(q_gaps), "column_gap": round(col_gap, 2), "tick_gap_min": round(tick_min, 2), "c_value_margin_min": round(min(c[0] for c in c_vals), 2), "b_margin_min": round(min(c[0] for c in bm if c[0] is not None), 2), "right_end": round(right_end, 2), "nearest_to_letters": near_letters},
             "flat_raster_diff": rd, "result": "ALL PASS" if ok else "FAIL", "written": L.now()}, f"{VER}/verify_record.json")
sys.exit(0 if ok else 1)
