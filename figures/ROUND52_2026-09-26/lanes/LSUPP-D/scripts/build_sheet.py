#!$T90_PY
"""Lane LSUPP-D (round 49): rebuild one sheet with every visible string one point larger. usage: build_sheet.py <Sheet>

Route (the six sheets are legacy overlay sheets, see lsd_common.py):
  1. base = the current sheet (V26; for Supp 9 and 10 the V15 vector twin of the V26 raster), re-distilled flat by Ghostscript
     pdfwrite when nested (Supp 16, 19), byte copy otherwise;
  2. the text layer as PyMuPDF sees it (child fitz_dump.py) and a no-text twin (child fitz_restamp.py with no items): the
     Ghostscript renders of base and twin differ only inside the PyMuPDF span boxes (proof that those spans are ALL the visible
     text); a span whose own exclusive area shows no glyph pixels is a covered copy (invisible) and is dropped;
  3. plan: every visible span gets its anchor rule (left, right, centre, rotated centre), size + 1.0, the same face, weight and colour;
  4. child fitz_restamp.py: whole-page text-only redaction, old font resources dropped, every planned string re-set;
     Supp 16: the band area is redacted too and the sheet is recomposed as base (y >= 79) + the coordinator's band_s16_r49.pdf;
  5. verification: page box = V26 (0.05 pt); Ghostscript txtwrite census new against the reference (visible strings equal, every
     size + 1.0, every dropped string a proven invisible copy); anchors kept within 0.3 pt (independent PyMuPDF read-back);
     raster diff at 150 dpi confined to the text boxes (zero pixels outside: graphics untouched); no two strings overlap;
     text-against-graphics intersections listed; fonts Arial faces only; 150 dpi render; checks.txt."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import json, os, re, shutil, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lsd_common as C

SHEET = sys.argv[1]
CFG = {"Supp_Fig08": dict(src=f"{C.V26}/Supp_Fig08.pdf", nested=False),
       "Supp_Fig09": dict(src=f"{C.V15}/Supp_Fig09.pdf", nested=False, v26=f"{C.V26}/Supp_Fig09.pdf"),
       "Supp_Fig10": dict(src=f"{C.V15}/Supp_Fig10.pdf", nested=False, v26=f"{C.V26}/Supp_Fig10.pdf"),
       "Supp_Fig16": dict(src=f"{C.V26}/Supp_Fig16.pdf", nested=True, band=C.BAND_R49, band_h=79.0),
       "Supp_Fig17": dict(src=f"{C.V26}/Supp_Fig17.pdf", nested=False),
       "Supp_Fig19": dict(src=f"{C.V26}/Supp_Fig19.pdf", nested=True)}[SHEET]
D = f"{C.LANE}/{SHEET}"; W = f"{D}/work"; V = f"{D}/verify"
for d_ in (W, V): os.makedirs(d_, exist_ok=True)
SRC = C.hydrated(CFG["src"]); V26REF = C.hydrated(CFG.get("v26", CFG["src"])); NEW = f"{D}/{SHEET}.pdf"
BAND = CFG.get("band"); BAND_H = CFG.get("band_h", 0.0)
ck = C.Checks(f"{SHEET}, lane LSUPP-D (round 49, {C.now()}): every visible string one point larger, graphics untouched. Source {SRC} (sha256 {C.sha256(SRC, 16)}), page-size reference {V26REF} (sha256 {C.sha256(V26REF, 16)}).")
n_pages, PW, PH = C.page_box(V26REF)
DPI = 150; S = DPI / 72.0

# ------------------------------------------------------------------ 1 base
BASE = f"{W}/{SHEET}_base.pdf"
if CFG["nested"]:
    C.gs_flatten(SRC, BASE); ck.info(f"base = Ghostscript pdfwrite re-distillation of the source (nested source), sha256 {C.sha256(BASE, 16)}")
else:
    shutil.copyfile(SRC, BASE); ck.info(f"base = byte copy of the source, sha256 {C.sha256(BASE, 16)}")
nb, bw, bh = C.page_box(BASE); assert nb == 1 and abs(bw - PW) < 0.05 and abs(bh - PH) < 0.05, (nb, bw, bh, PW, PH)

# ------------------------------------------------------------------ 2 text layer and visibility
C.wd(f"{SHEET}_dump_base", [C.PY, f"{C.SCRIPTS}/fitz_dump.py", BASE, f"{W}/base_spans.json"])
BD = json.load(open(f"{W}/base_spans.json")); SP = BD["spans"]
assert BD["xobjects"] == 0 or SHEET == "Supp_Fig16", ("base still nested", BD["xobjects"])
job = {"src": BASE, "out": f"{W}/{SHEET}_notext.pdf", "whole_page": True, "items": []}
json.dump(job, open(f"{W}/notext_job.json", "w")); C.wd(f"{SHEET}_notext", [C.PY, f"{C.SCRIPTS}/fitz_restamp.py", f"{W}/notext_job.json"])
C.gs_render(BASE, f"{W}/base_{DPI}.png", DPI); C.gs_render(f"{W}/{SHEET}_notext.pdf", f"{W}/notext_{DPI}.png", DPI)
A = C.load_png(f"{W}/base_{DPI}.png"); N0 = C.load_png(f"{W}/notext_{DPI}.png")
tdiff = np.abs(A - N0).max(axis=2) > 0
def mask_of(boxes, pad=1.0):
    m = np.zeros(tdiff.shape, bool)
    for b in boxes:
        m[max(0, int((b[1] - pad) * S)):int((b[3] + pad) * S) + 1, max(0, int((b[0] - pad) * S)):int((b[2] + pad) * S) + 1] = True
    return m
inside = mask_of([s["bbox"] for s in SP])
n_out = int((tdiff & ~inside).sum())
ck.log("removing every text object changes pixels only inside the PyMuPDF span boxes (those spans are all the visible text; every other string of the txtwrite census is invisible)", n_out == 0, f"{int(tdiff.sum())} glyph pixels at {DPI} dpi, {n_out} outside")
vis, hid, amb = [], [], []
for i, s in enumerate(SP):
    own = mask_of([s["bbox"]], pad=0.0); others = mask_of([t["bbox"] for j, t in enumerate(SP) if j != i], pad=0.3)
    excl = own & ~others
    if excl.sum() < 0.3 * own.sum():
        amb.append(i); continue
    (vis if (tdiff & excl).sum() > 0 else hid).append(i)
for i in amb:      # too much overlap with neighbours for the exclusive test: redact this span alone and look
    s = SP[i]; tight = [s["bbox"][0], s["origin"][1] - 0.70 * s["size"], s["bbox"][2], s["origin"][1] + 0.12 * s["size"]]
    job = {"src": BASE, "out": f"{W}/vis_test_{i}.pdf", "whole_page": False, "redact_rects": [tight], "items": []}
    json.dump(job, open(f"{W}/vis_test_{i}.json", "w")); C.wd(f"{SHEET}_vis_{i}", [C.PY, f"{C.SCRIPTS}/fitz_restamp.py", f"{W}/vis_test_{i}.json"])
    C.gs_render(f"{W}/vis_test_{i}.pdf", f"{W}/vis_test_{i}.png", DPI); B = C.load_png(f"{W}/vis_test_{i}.png")
    changed = int((np.abs(A - B).max(axis=2) > 0).sum())
    (vis if changed > 0 else hid).append(i)
    os.remove(f"{W}/vis_test_{i}.pdf"); os.remove(f"{W}/vis_test_{i}.png")
vis.sort(); hid.sort()
if BAND:
    band_spans = [i for i in vis if SP[i]["origin"][1] < BAND_H + 1.0]
    vis = [i for i in vis if i not in band_spans]
else:
    band_spans = []
json.dump({"visible": vis, "hidden": hid, "ambiguous_tested": amb, "band": band_spans, "spans": SP}, open(f"{W}/visibility.json", "w"), indent=0)
ck.info(f"PyMuPDF spans {len(SP)}: visible {len(vis)}, covered (no glyph pixel of their own, dropped) {len(hid)}" + (f", band text replaced by the coordinator's band {len(band_spans)}" if BAND else "") + f"; {len(amb)} decided by a single-span redaction test")
for i in hid:
    s = SP[i]; ck.info(f"  covered copy dropped: {s['text']!r} {s['size']} pt at ({s['origin'][0]:.1f}, {s['origin'][1]:.1f})")

# ------------------------------------------------------------------ 3 plan
def classify(s):
    t = s["text"]; ox, oy = s["origin"]; x0, y0, x1, y1 = s["bbox"]; sz = round(s["size"] * 2) / 2; bold = "Bold" in s["font"]
    if s["dir"][1] < -0.5: return "rot"
    num = re.fullmatch(r"[<>]?\d[\d.,]*", t) is not None
    val = re.fullmatch(r"\d\.\d\d \(\d\.\d\d-\d\.\d\d\)", t) is not None
    if SHEET in ("Supp_Fig08", "Supp_Fig09", "Supp_Fig10"):
        if abs(ox - 27.36) < 0.5: return "left"            # row labels
        if t.startswith("Short sleep"): return "left"        # key labels, handle to the left
        if num: return "center"                              # x tick labels
        if sz == 11: return "center"                         # axis title
    elif SHEET == "Supp_Fig17":
        if oy < 25: return "left"                            # key row, handles to the left
        if t in ("HR (95% CI)", "P"): return "right"         # headers share the columns' right edges
        if val or (num and ox > 440): return "right"         # HR (95% CI) and P columns
        if ox < 15: return "left"                            # row labels
        if num: return "center"                              # x tick labels
        if sz == 11: return "center"                         # axis title
    elif SHEET == "Supp_Fig19":
        if "Helvetica" in s["font"]: return "left"           # key lines (base-14 Helvetica 9 pt, re-set in the sheet's Arial)
        if t in ("HR (95% CI)", "P"): return "center"        # headers sit on their own centres, not on the columns' edges
        if val or (num and ox > 320): return "right"
        if ox < 15: return "left"                            # row labels and bold block titles
        if num: return "center"
        if sz == 11: return "center"
    elif SHEET == "Supp_Fig16":
        if sz == 13: return "left"                           # panel letters
        if t in ("HR (95% CI)", "P"): return "right"
        if val or (num and ox > 440 and 500 < oy < 760): return "right"    # panel b columns
        if num and ox < 40: return "right"                   # y tick labels of panels a and c
        if num: return "center"                              # x tick labels
        if sz == 11: return "center"                         # horizontal axis titles
        if bold and 880 < oy < 910: return "center"          # facet titles of panel c
        if t in ("Primary", "SHHS", "MrOS"): return "center" # category labels of panel c
        return "left"                                        # row labels and key labels (handles to the left)
    raise ValueError((SHEET, t, ox, oy))
ITEMS = []
for i in vis:
    s = SP[i]; osz = round(s["size"] * 2) / 2; al = classify(s)
    ITEMS.append({"text": s["text"], "old_size": osz, "size": round(osz + 1.0, 2), "bold": "Bold" in s["font"], "color": s["color"], "ox": s["origin"][0], "oy": s["origin"][1],
                  "align": al if al != "rot" else "rot", "rot": -90 if al == "rot" else 0, "old_bbox": s["bbox"], "old_font": s["font"], "old_size_exact": s["size"], "span_index": i})
json.dump(ITEMS, open(f"{W}/plan.json", "w"), indent=1)
ck.info("plan: " + ", ".join(f"{k} {v}" for k, v in sorted(__import__('collections').Counter(it['align'] for it in ITEMS).items())) + "; sizes " + ", ".join(f"{a} -> {b}" for a, b in sorted(set((it['old_size'], it['size']) for it in ITEMS))))

# ------------------------------------------------------------------ 4 restamp (and compose)
FLAT_NEW = f"{W}/{SHEET}_restamped.pdf" if BAND else NEW
job = {"src": BASE, "out": FLAT_NEW, "whole_page": True, "redact_rects": ([[0.0, 0.0, PW, BAND_H + 1.0]] if BAND else []), "items": ITEMS,
       "metadata": {"title": SHEET, "creator": f"LSUPP-D r49 build_sheet.py: text-only redaction + TextWriter re-set of every visible string at +1 pt on {os.path.basename(SRC)}" + (" + band_s16_r49" if BAND else "")}}
json.dump(job, open(f"{W}/restamp_job.json", "w")); st, out = C.wd(f"{SHEET}_restamp", [C.PY, f"{C.SCRIPTS}/fitz_restamp.py", f"{W}/restamp_job.json"]); ck.info(out.strip().splitlines()[-1])
DR = json.load(open(FLAT_NEW + ".drawn.json"))["drawn"]
if BAND:
    C.hydrated(BAND); nbp, bwp, bhp = C.page_box(BAND); assert nbp == 1 and abs(bwp - PW) < 0.05 and abs(bhp - BAND_H) < 0.05, (bwp, bhp)
    job = {"out": NEW, "page": [PW, PH], "parts": [{"pdf": FLAT_NEW, "rect": [0, BAND_H, PW, PH], "clip": [0, BAND_H, PW, PH]}, {"pdf": BAND, "rect": [0, 0, bwp, bhp]}],
           "metadata": {"title": SHEET, "creator": "LSUPP-D r49 build_sheet.py: re-set text at +1 pt on the flattened V26 base (y >= 79) + band_s16_r49 (coordinator, LSKETCH)"}}
    json.dump(job, open(f"{W}/compose_job.json", "w")); C.wd(f"{SHEET}_compose", [C.PY, f"{C.SCRIPTS}/fitz_compose.py", f"{W}/compose_job.json"])
    ck.info(f"band {BAND} (sha256 {C.sha256(BAND, 16)}, {bwp:.2f} x {bhp:.2f} pt) placed at the top, base clipped to y >= {BAND_H}")

# ------------------------------------------------------------------ 5 verify
nn, nw, nh = C.page_box(NEW)
ck.log("one page, page box equal to V26 within 0.05 pt", nn == 1 and abs(nw - PW) < 0.05 and abs(nh - PH) < 0.05, f"{nw:.3f} x {nh:.3f} pt (V26 {PW:.3f} x {nh and PH:.3f})")
# anchors: independent PyMuPDF read-back of the flat new sheet
C.wd(f"{SHEET}_dump_new", [C.PY, f"{C.SCRIPTS}/fitz_dump.py", FLAT_NEW, f"{W}/new_spans.json", "--no-drawings"])
ND = json.load(open(f"{W}/new_spans.json")); NS = ND["spans"]
worst = {"left": 0.0, "right": 0.0, "center": 0.0, "rot": 0.0}; missing = []; nfont_bad = []
for it, dr in zip(ITEMS, DR):
    if it["rot"]:
        c = [s for s in NS if s["text"] == it["text"] and s["dir"][1] < -0.5 and abs(s["origin"][0] - it["ox"]) < 0.3]
        if not c: missing.append(it["text"]); continue
        s = c[0]; oc = (it["old_bbox"][1] + it["old_bbox"][3]) / 2; nc = (s["bbox"][1] + s["bbox"][3]) / 2
        worst["rot"] = max(worst["rot"], abs(nc - oc), abs(s["origin"][0] - it["ox"]))
    else:
        c = [s for s in NS if s["text"] == it["text"] and abs(s["origin"][1] - it["oy"]) < 0.3 and abs(s["origin"][0] - dr["x0"]) < 0.3]
        if not c: missing.append(it["text"]); continue
        s = c[0]; ob, nb_ = it["old_bbox"], s["bbox"]
        if it["align"] == "left": worst["left"] = max(worst["left"], abs(nb_[0] - ob[0]))
        elif it["align"] == "right": worst["right"] = max(worst["right"], abs(nb_[2] - ob[2]))
        else: worst["center"] = max(worst["center"], abs((nb_[0] + nb_[2]) / 2 - (ob[0] + ob[2]) / 2))
    if ("Bold" in s["font"]) != it["bold"] or abs(s["size"] - it["size"]) > 0.01: nfont_bad.append((it["text"], s["font"], s["size"]))
ck.log("every planned string read back on its baseline at + 1.0 pt, same weight, on its own anchor within 0.3 pt (left edge, right edge, centre, or baseline x and centre of a rotated title)", not missing and not nfont_bad and max(worst.values()) <= 0.3,
       f"largest anchor offsets: " + ", ".join(f"{k} {v:.3f}" for k, v in worst.items()) + f"; missing {missing}; weight or size mismatches {nfont_bad}")
fonts_new = sorted(set(f[1] for f in ND["fonts"]))
ck.log("fonts of the re-set text: the sheet's Arial faces only (TextWriter embedded subsets)", all(f in ("Arial Regular", "Arial Bold") for f in fonts_new) and fonts_new, f"{fonts_new}")
# census: Ghostscript txtwrite, the new sheet against the reference (Supp 9 and 10: the V15 vector twin, V26 being a raster without text)
from collections import Counter
ref_sp = C.gs_txtwrite(SRC, f"{W}/census_ref.xml"); new_sp = C.gs_txtwrite(NEW, f"{W}/census_new.xml")
half = lambda v: round(v * 2) / 2
new_body = [s for s in new_sp if not (BAND and s["y"] < BAND_H)]; ref_body = [s for s in ref_sp if not (BAND and s["y"] < BAND_H)]
planned = Counter((it["text"].strip(), it["size"]) for it in ITEMS); got = Counter((s["text"], half(s["size"])) for s in new_body)
ck.log(f"txtwrite census of the new sheet = the planned strings, one run each at + 1.0 pt ({sum(planned.values())} runs)", planned == got, f"difference {dict((planned - got) + (got - planned))}")
VIS = [SP[i] for i in vis]
def piece_of(run):
    for s_ in VIS:
        if abs(s_["origin"][1] - run["y"]) <= 1.5 and abs(half(s_["size"]) - half(run["size"])) <= 0.011 and run["text"] in s_["text"] and s_["bbox"][0] - 1.5 <= run["x0"] and run["x1"] <= s_["bbox"][2] + 1.5:
            return s_
    return None
pieces = {}; ghosts = []
for r in ref_body:
    s_ = piece_of(r)
    if s_ is None: ghosts.append(r)
    else: pieces.setdefault(id(s_), []).append(r)
uncovered = [s_["text"] for s_ in VIS if id(s_) not in pieces]
multi = [(s_["text"], len(rs), sorted(set(r["text"] for r in rs), key=lambda t: -len(t))) for s_ in VIS for rs in [pieces.get(id(s_), [])] if len(rs) > 1]
def ghost_class(q):
    if q["x1"] < 0 or q["x0"] > PW or q["y"] < 0 or q["y"] > PH: return "off-page"
    hits = [SP[i] for i in hid if abs(SP[i]["origin"][1] - q["y"]) <= 1.5 and abs(SP[i]["origin"][0] - q["x0"]) <= 1.5 and q["text"] in SP[i]["text"]]
    return "covered copy (no glyph pixel of its own)" if hits else "clipped copy of an earlier generation (not in the PyMuPDF view, no pixel outside the visible spans)"
ck.log(f"every run of the reference census is a piece of a visible string ({len(ref_body) - len(ghosts)} runs over {len(VIS)} visible strings: run splits and stacked copies) or an invisible copy ({len(ghosts)} runs), every visible string has at least one reference run", not uncovered, f"visible strings without a reference run {uncovered}")
for text, n, distinct in multi: ck.info(f"  reference runs of one string, one run on the new sheet: {text!r} <- {n} runs {distinct}")
for g in ghosts: ck.info(f"  not carried over: {g['text']!r} {g['size']} pt at ({g['x0']:.0f}, {g['y']:.0f}): {ghost_class(g)}")
if BAND:
    band_own = Counter((s_["text"], half(s_["size"])) for s_ in C.gs_txtwrite(BAND, f"{W}/census_band.xml"))
    band_rows = Counter((s_["text"], half(s_["size"])) for s_ in new_sp if s_["y"] < BAND_H)
    ck.log("band region census equals the coordinator's band_s16_r49 own census, no r44 band run left", band_rows == band_own, f"band runs {sum(band_rows.values())}; r44 band runs of the reference dropped {len(ref_sp) - len(ref_body)}")
sizes = sorted(set((it["old_size"], it["size"]) for it in ITEMS))
ck.log("every re-set string is exactly 1.0 pt larger than its old size", all(abs(b - a - 1.0) < 1e-9 for a, b in sizes), f"{sizes}")
# raster: graphics untouched
C.gs_render(NEW, f"{W}/new_{DPI}.png", DPI); Bn = C.load_png(f"{W}/new_{DPI}.png")
assert Bn.shape == A.shape, (Bn.shape, A.shape)
diff = np.abs(A - Bn).max(axis=2) > 0
allowed = mask_of([SP[i]["bbox"] for i in vis] + [SP[i]["bbox"] for i in hid] + [d_["box"] for d_ in DR], pad=1.0)
if BAND:
    allowed[:int((BAND_H + 1.0) * S) + 1, :] = True
n_outside = int((diff & ~allowed).sum())
ck.log(f"raster diff at {DPI} dpi against the base confined to the old and new text boxes" + (" and the band region" if BAND else "") + " (zero differing pixels outside: every mark, rule, key handle and colour untouched)", n_outside == 0, f"{int(diff.sum())} differing pixels, {n_outside} outside")
if BAND:
    # the band placed alone on a blank page of the sheet's size by the same placement call: the composed band region must equal it
    job = {"out": f"{W}/band_alone.pdf", "page": [PW, PH], "parts": [{"pdf": BAND, "rect": [0, 0, bwp, bhp]}]}
    json.dump(job, open(f"{W}/band_alone_job.json", "w")); C.wd(f"{SHEET}_band_alone", [C.PY, f"{C.SCRIPTS}/fitz_compose.py", f"{W}/band_alone_job.json"])
    C.gs_render(f"{W}/band_alone.pdf", f"{W}/band_alone_{DPI}.png", DPI); Ba = C.load_png(f"{W}/band_alone_{DPI}.png"); hb = int(BAND_H * S)
    bdiff = int((np.abs(Bn[:hb] - Ba[:hb]).max(axis=2) > 0).sum())
    ck.log("the band region of the new sheet renders exactly as band_s16_r49.pdf placed alone on a blank page by the same call (band content and placement intact)", bdiff == 0, f"{bdiff} differing pixels in the top {hb} rows")
    C.gs_render(BAND, f"{W}/band_{DPI}.png", DPI); Bb = C.load_png(f"{W}/band_{DPI}.png"); hb2 = min(Bb.shape[0], hb)
    pd_ = np.abs(Bn[:hb2, :Bb.shape[1]] - Bb[:hb2]).max(axis=2)
    ck.info(f"band rendered as its own page against placed as a form: {int((pd_ > 0).sum())} pixels differ at {DPI} dpi (max {int(pd_.max())}), Ghostscript's transparency-group and icon anti-aliasing of the Skia band, judged on the render")
if CFG["nested"]:
    C.gs_render(SRC, f"{W}/src_{DPI}.png", DPI); Asrc = C.load_png(f"{W}/src_{DPI}.png"); fmask = np.abs(Asrc - A).max(axis=2) > 0
    if BAND: fmask[:int((BAND_H + 1.0) * S) + 1, :] = False
    fd = int(fmask.sum())
    ck.log("the flattened base renders as the nested source" + (" below the band" if BAND else "") + " (pdfwrite re-distillation lossless at 150 dpi)", fd == 0, f"{fd} differing pixels")
# overlaps, on glyph INK extents (per-glyph bounds of the Arial faces via fontTools), not on font boxes
try:
    from fontTools.ttLib import TTFont
    _FT = {}
    def glyph_extents(text, bold):
        key = "b" if bold else "r"
        if key not in _FT:
            f_ = TTFont(paths.FONT_ARIAL_BOLD if bold else paths.FONT_ARIAL); _FT[key] = (f_, f_["head"].unitsPerEm, f_.getBestCmap(), f_["glyf"])
        f_, upm, cmap, glyf = _FT[key]; ymax, ymin = -9.0, 9.0
        for c in text:
            if c == " " or ord(c) not in cmap: continue
            g = glyf[cmap[ord(c)]]
            if g.numberOfContours == 0: continue
            ymax = max(ymax, g.yMax / upm); ymin = min(ymin, g.yMin / upm)
        return (ymax, ymin) if ymax > -9 else (0.716, 0.0)
    ink_source = "fontTools per-glyph bounds"
except Exception:
    def glyph_extents(text, bold):
        up = 0.716 + (0.02 if any(c in "bdfhklt()[]{}|/" for c in text) else 0.0); dn = -0.21 if any(c in "gjpqy()[]{},;" for c in text) else 0.0
        return (up, dn)
    ink_source = "class heuristic (fontTools absent)"
def glyph_boxes(text, bold, size, x0, baseline, rot=0, ox=None, y_start=None):
    """Per-glyph ink boxes (page coordinates) of a string set with the TextWriter face: x from the accumulated advances and the
    glyph's own xMin/xMax, y from its yMin/yMax. Rotated (-90) strings run upward from (ox, y_start), ascender side to the left."""
    try:
        key = "b" if bold else "r"
        if key not in _FT:
            f_ = TTFont(paths.FONT_ARIAL_BOLD if bold else paths.FONT_ARIAL); _FT[key] = (f_, f_["head"].unitsPerEm, f_.getBestCmap(), f_["glyf"])
        f_, upm, cmap, glyf = _FT[key]; hmtx = f_["hmtx"]; pen = 0.0; out = []
        for c in text:
            gn = cmap.get(ord(c)); adv = hmtx[gn][0] / upm * size if gn else 0.0
            if gn and glyf[gn].numberOfContours != 0:
                g = glyf[gn]; gx0, gx1 = pen + g.xMin / upm * size, pen + g.xMax / upm * size; gy1, gy0 = g.yMax / upm * size, g.yMin / upm * size
                if rot == -90: out.append([ox - gy1, y_start - gx1, ox - gy0, y_start - gx0])
                else: out.append([x0 + gx0, baseline - gy1, x0 + gx1, baseline - gy0])
            pen += adv
        return out
    except Exception:
        ymax, ymin = glyph_extents(text, bold); w = len(text) * 0.5 * size
        return [[x0, baseline - ymax * size, x0 + w, baseline - ymin * size]]
def ink_box(it, dr):
    ymax, ymin = glyph_extents(it["text"], it["bold"]); sz = it["size"]
    if it["rot"]: return [dr["x0"] - ymax * sz, dr["box"][1], dr["x0"] - ymin * sz, dr["box"][3]]
    return [dr["box"][0], it["oy"] - ymax * sz, dr["box"][2], it["oy"] - ymin * sz]
def old_ink_box(it):
    ymax, ymin = glyph_extents(it["text"], it["bold"]); sz = it["old_size"]; ob = it["old_bbox"]
    if it["rot"]: return [it["ox"] - ymax * sz, ob[1], it["ox"] - ymin * sz, ob[3]]
    return [ob[0], it["oy"] - ymax * sz, ob[2], it["oy"] - ymin * sz]
boxes = [(it["text"], ink_box(it, d_)) for it, d_ in zip(ITEMS, DR)]
GB = [glyph_boxes(it["text"], it["bold"], it["size"], d_["box"][0], it["oy"], rot=it["rot"], ox=d_.get("x0"), y_start=d_.get("y_start")) for it, d_ in zip(ITEMS, DR)]
def inter(a, b, shrink=0.2):
    return not (a[2] - shrink <= b[0] + shrink or b[2] - shrink <= a[0] + shrink or a[3] - shrink <= b[1] + shrink or b[3] - shrink <= a[1] + shrink)
def gap(a, b):
    dx = max(a[0] - b[2], b[0] - a[2], 0.0); dy = max(a[1] - b[3], b[1] - a[3], 0.0); return (dx ** 2 + dy ** 2) ** 0.5
pairs = []; gaps = []
for i_ in range(len(boxes)):
    for j_ in range(i_ + 1, len(boxes)):
        if not inter(boxes[i_][1], boxes[j_][1], shrink=-1.0): continue          # string boxes further than 1 pt apart: no glyph can touch
        g = min(gap(ga, gb) for ga in GB[i_] for gb in GB[j_]) if GB[i_] and GB[j_] else 9.9
        if any(inter(ga, gb) for ga in GB[i_] for gb in GB[j_]): pairs.append((boxes[i_][0], boxes[j_][0]))
        gaps.append((round(g, 2), boxes[i_][0], boxes[j_][0]))
gaps.sort()
ck.log(f"no two re-set strings overlap (per-glyph ink boxes, {ink_source}, 0.2 pt tolerance) and every string lies inside the page", not pairs and all(b[0] >= -0.5 and b[1] >= -0.5 and b[2] <= PW + 0.5 and b[3] <= PH + 0.5 for _, b in boxes), f"overlapping pairs {pairs}; tightest glyph clearances {gaps[:6]}")
DRW = [d_ for d_ in BD["drawings"] if not (d_["fill"] == "#ffffff" and d_["color"] is None)]     # white fills are covers, not marks
newhits = []
for it, d_ in zip(ITEMS, DR):
    nb_, ob_ = ink_box(it, d_), old_ink_box(it)
    for g in DRW:
        r = g["rect"]
        if (r[2] - r[0]) > 0.9 * PW: continue
        if inter(nb_, r, shrink=0.0) and not inter(ob_, r, shrink=0.0):
            newhits.append((it["text"], [round(v, 1) for v in r], g["color"] or g["fill"]))
ck.info(f"drawing items newly touched by a larger string's ink box (judged on the render): {len(newhits)}" + ("".join(f"\n        {h}" for h in newhits[:40])))
alltext = " ".join(it["text"] for it in ITEMS)
ck.log("house rules on the re-set strings (no em dash, no semicolon), wording unchanged", "—" not in alltext and ";" not in alltext, "ok")
C.gs_render(NEW, f"{D}/{SHEET}_{DPI}dpi.png", DPI)
ck.log(f"{DPI} dpi render written (Ghostscript)", os.path.exists(f"{D}/{SHEET}_{DPI}dpi.png"), f"{D}/{SHEET}_{DPI}dpi.png")
json.dump({"sheet": SHEET, "source": SRC, "source_sha256": C.sha256(SRC), "page_size_reference": V26REF, "reference_sha256": C.sha256(V26REF), "base": BASE, "base_sha256": C.sha256(BASE),
           "band": BAND, "band_sha256": C.sha256(BAND) if BAND else None, "new": NEW, "new_sha256": C.sha256(NEW), "page": [PW, PH], "n_visible": len(vis), "n_hidden": len(hid), "n_ghost_rows": len(ghosts),
           "sizes": sizes, "written": C.now()}, open(f"{W}/build_record.json", "w"), indent=1)
ck.write(f"{V}/checks.txt")
