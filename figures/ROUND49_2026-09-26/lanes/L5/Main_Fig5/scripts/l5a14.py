#!/usr/bin/env python3
"""Round 49, lane L5 (2026-09-26): the round-42 copy of this library with the page helper Sheet drawing every text ONE POINT larger
(Sheet.BUMP = 1.0 inside text() and width(); width(..., bump=False) gives the V13 width for the geometry read-back of the V13 sheet).
Lane V14_L5a_PAP_MAIN (round 37, 2026-09-15): shared library. ONE spec module for every input path (integrity gate 7), the v8.1
sidecar gate (v14lib), the builders' printing rules (half up on the file value, the round-30 close-out rule), the V13 geometry
read-back on flat sheets (get_drawings allowed on ED 8, 9, 10 only), the PyMuPDF drawing helpers (fitz.Shape + fitz.TextWriter with
the system Arial faces, the ROUND31 LFOREST ED 9 idiom), drawn-record collection for v14lib.printed_values_csv, renders, crops.
Lineage: ROUND30_2026-09-08/V7_L5a_PAP_MAIN/_lib_l5a.py (copied beside as _lib_l5a_PRE_V8_1.py), never edited there."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import os, re, json, csv, sys, subprocess, time, shutil
from collections import Counter
from decimal import Decimal, ROUND_HALF_UP
from math import floor, log10, log, exp
import numpy as np
import pandas as pd
import fitz
from PIL import Image
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND37_2026-09-14/figures/_common")
from v14lib import (T90, NUM, SV, V13, FIG, hydrated, sha256, load_json, sidecar_v8, spans as v14_spans, collapse as v14_collapse,
                    multisets as v14_multisets, letters as v14_letters, page_box, fonts_of, render as v14_render, raster_diff, Checks,
                    write_changes, halfup, printed_values_csv, RULES)

LANE = f"{FIG}/lanes/V14_L5a_PAP_MAIN"
ARIAL = paths.FONT_ARIAL; ARIAL_BOLD = paths.FONT_ARIAL_BOLD
SPEC = {
    "treatment_v2": f"{NUM}/treatment_v2.json",                                            # step 117 (v8.1, 2026-09-14 13:38)
    "gaps_v2": f"{NUM}/gaps_v2.json",                                                      # step 118
    "fig5b_refit_json": f"{paths.FIGURE_ROOT}/ROUND12_2026-09-03/L4_fig5/L4_fig5b_refit.json",          # step 155 (v8.1, 2026-09-14 15:26)
    "fig5b_refit_csv": f"{paths.FIGURE_ROOT}/ROUND12_2026-09-03/L4_fig5/L4_fig5b_refit.csv",            # step 155
    "treatment_v3": f"{NUM}/treatment_v3.json",                                            # step 211 (v8.1, 2026-09-15 09:10): the comparison ADJUSTED for prevalent cardiopulmonary disease (decision 6), FIX A1
    "fig5b_refit_v3_json": f"{paths.FIGURE_ROOT}/ROUND12_2026-09-03/L4_fig5/L4_fig5b_refit_v3.json",    # step 212 (v8.1, 2026-09-15 09:10): the band refit with the same adjustment, FIX A1
    "results_continuous": f"{SV}/pap_residual_continuous/results_continuous.csv",         # step 147 (v8.1, 2026-09-14 15:19)
    "results_graded": f"{SV}/pap_residual_continuous/results_graded.csv",                 # step 147
    "disease_definitions_dir": paths.ANALYSIS_DIR,                                                        # numbers/disease_definitions.py (v8, ORGAN_GROUP)
    "cohorts": f"{NUM}/cohorts.json",                                                      # step 101
}
INK = "#1a1d21"; INK20 = "#1a1d20"; BLUE = "#0288d1"; BLUE_LT = "#b3dcf2"; MID = "#3f9fd8"; GREY = "#8a9099"; GREY_LT = "#ccd1d6"; CTRL_GREY = "#aeb5bc"
ORANGE = "#d55e00"; WHITE = "#ffffff"


# ------------------------------------------------------------------ colours
def rgb(h):
    h = h.lstrip("#"); return tuple(int(h[i:i + 2], 16) / 255 for i in (0, 2, 4))


def hexcol(c):
    return None if c is None else "#%02x%02x%02x" % tuple(int(round(v * 255)) for v in c[:3])


# ------------------------------------------------------------------ the builders' printing rules, half up on the file value
def fmt_p(p):
    """Two significant figures, "<0.001" below 0.001, one further decimal where rounding would land on 0.050 (the generators'
    rule, L4_build_fig5b_new.py / make_efig20 / figure4C), the digits half up on the stored value (round-30 close-out rule)."""
    p = float(p)
    if p < 0.001: return "<0.001"
    d = 1 - floor(log10(p)); s = halfup(p, d)
    while s in ("0.050", "0.0500") and abs(p - 0.05) > 1e-12 and d < 6:
        d += 1; s = halfup(p, d)
    return s


def r2(v): return halfup(v, 2)
def hr_ci(hr, lo, hi): return f"{r2(hr)} ({r2(lo)}-{r2(hi)})"
def pct_text(pct): return f"{halfup(pct, 1)}%"
def pct_text_minus(pct): return pct_text(pct).replace("-", "−")   # ED 10 prints the true minus glyph
def int_comma(v): return f"{int(round(float(v))):,}"
RULES.update({"fmt_p": fmt_p, "hr_ci": lambda v: hr_ci(*v), "pct_text": pct_text, "pct_text_minus": pct_text_minus, "int_comma": int_comma,
              "count_pct": lambda v: f"n = {int_comma(v[0])} ({round(100 * v[0] / v[1])}%)", "group_k_of_n": lambda v: f"{v[0]} ({v[1]} of {v[2]})",
              "events_triplet": lambda v: f"{int(v[0])}/{int(v[1])}/{int(v[2])}", "n_eq": lambda v: f"n = {int_comma(v)}", "band_n": lambda v: f"{v[0]}, n = {int_comma(v[1])}",
              "star_p05": lambda v: "*" if float(v) < 0.05 else "", "label": lambda v: PRINT_LABEL.get(str(v), str(v))})
# Row labels wider than the sheets' label gutters: the ranking family (lane V14_L2_T90, Fig 2 and ED 3) prints the v8 outcome
# "Ventricular arrhythmia or cardiac arrest" (172.8 pt at 10 pt Arial) as "Ventricular arrhythmia/arrest" (126.2 pt); the same
# short form is used here so the figures agree (declared in the report; the numbers files keep the long name).
PRINT_LABEL = {"Ventricular arrhythmia or cardiac arrest": "Ventricular arrhythmia/arrest"}
def plab(k): return PRINT_LABEL.get(k, k)


# ------------------------------------------------------------------ data
V8_2_CUT = "2026-09-15 19:18:00"   # round v8.2 (2026-09-15 evening): the group P inputs were re-extracted and the chain rerun; every source must be dated after this
def gate(key):
    sc = sidecar_v8(SPEC[key]); sc["key"] = key
    mt = sc.get("output_mtime") or ""
    if mt < V8_2_CUT: raise SystemExit(f"V8.2 GATE: {SPEC[key]} sidecar output mtime {mt!r} is before {V8_2_CUT} (stale for round v8.2)")
    sc["v8_2_gate"] = f"output mtime {mt} >= {V8_2_CUT}"; return sc


def treatment_rows():
    """The 49 corrected_vs_not outcomes of treatment_v2.json sorted by the hazard ratio descending (ED 8 and Fig 5b row rule)."""
    T = load_json(SPEC["treatment_v2"]); oc = T["corrected_vs_not"]["outcomes"]
    return T, sorted(oc.items(), key=lambda kv: -kv[1]["hr"])


def residual_rows():
    """Panel c of Fig 5 and panel a of ED 9: the Benjamini-Hochberg survivors of results_continuous.csv by hr_sd descending."""
    c = pd.read_csv(hydrated(SPEC["results_continuous"]))
    return c, c[c.sig_fdr == True].sort_values("hr_sd", ascending=False)  # noqa: E712


def graded_rows():
    """ED 9 panel b: the principal outcomes of results_graded.csv by hr_gt10 descending."""
    g = pd.read_csv(hydrated(SPEC["results_graded"]))
    return g, g[g.principal == True].sort_values("hr_gt10", ascending=False)  # noqa: E712


def fig5b_rows():
    """Fig 5b: the conditions significant in the two-group comparison (in_current_5b) and the controls, both by pub_hr descending."""
    R = load_json(SPEC["fig5b_refit_json"]); OC = R["outcomes"]
    SIG = sorted(((k, v) for k, v in OC.items() if v["in_current_5b"]), key=lambda kv: -kv[1]["pub_hr"])
    CTRL = sorted(((k, v) for k, v in OC.items() if v["negative_control"]), key=lambda kv: -kv[1]["pub_hr"])
    return R, SIG, CTRL


def organ_rollup():
    """ED 10 roll-up recomputed from treatment_v2.json with the grouping of numbers/disease_definitions.py (the v8 module, the
    make_etables_v2.py rule verbatim: controls held out, composite and mortality excluded, Responding = P < .05 and HR > 1,
    median of the per-outcome percentage)."""
    sys.path.insert(0, SPEC["disease_definitions_dir"])
    import disease_definitions as dd
    NAME = {k: v[0] for k, v in dd.DISEASES.items()}
    GROUP_BY_NAME = {v[0]: dd.ORGAN_GROUP[k] for k, v in dd.DISEASES.items()}
    GROUP_BY_NAME["Death from any cause"] = "Mortality"
    NEG_NAMES = {NAME[k] for k in dd.NEGATIVE_CONTROLS}
    T = load_json(SPEC["treatment_v2"]); cv = T["corrected_vs_not"]["outcomes"]
    unmapped = [c for c in cv if c not in GROUP_BY_NAME]
    assert not unmapped, f"outcomes without an organ group in disease_definitions.py: {unmapped}"
    grp = []
    for name in sorted(set(GROUP_BY_NAME.values())):
        if name in ("Composite", "Mortality"): continue
        negs = [c for c in cv if GROUP_BY_NAME.get(c) == name and c in NEG_NAMES]
        real = [(c, v) for c, v in cv.items() if GROUP_BY_NAME.get(c) == name and c not in NEG_NAMES]
        if not real: continue
        grp.append({"group": name, "conditions": len(real), "controls_held_out": len(negs),
                    "responding": sum(1 for _, v in real if v["p"] < .05 and v["hr"] > 1),
                    "median_pct_lower": float(np.median([v["pct_lower_if_corrected"] for _, v in real])),
                    "members": [c for c, _ in real], "control_members": negs})
    grp.sort(key=lambda g: -g["median_pct_lower"])
    ungrouped = [c for c in cv if GROUP_BY_NAME.get(c) in ("Composite", "Mortality")]
    mapping = {c: GROUP_BY_NAME[c] for c in cv}
    return grp, ungrouped, mapping


# ------------------------------------------------------------------ V13 text layer and drawings (flat sheets only)
def text_spans(pdf):
    d = fitz.open(hydrated(pdf)); p = d[0]; out = []
    for b in p.get_text("dict")["blocks"]:
        for l in b.get("lines", []):
            for s in l["spans"]:
                t = s["text"].replace("\xa0", " ").strip()
                if t: out.append({"text": t, "x": round(s["origin"][0], 2), "y": round(s["origin"][1], 2), "font": s["font"], "size": round(s["size"], 2),
                                  "color": "%06x" % s["color"], "bbox": [round(v, 2) for v in s["bbox"]]})
    info = {"rect": [round(v, 4) for v in p.rect], "xobjects": len(p.get_xobjects()), "fonts": sorted({f[3].split("+")[-1] for f in p.get_fonts(full=True)})}
    d.close(); return out, info


def collapse(spans, tol=1.5):
    keep = []
    for s in spans:
        if not any(s["text"] == k["text"] and abs(s["x"] - k["x"]) < tol and abs(s["y"] - k["y"]) < tol for k in keep): keep.append(s)
    return keep


def norm(t): return t.replace("−", "-").replace("–", "-").replace("—", "-").replace("\xad", "-").strip()


def find_span(spans, text, y=None, ytol=2.0, xmin=None, xmax=None, ymin=None, ymax=None):
    c = [s for s in spans if norm(s["text"]) == norm(text) and (xmin is None or s["x"] >= xmin) and (xmax is None or s["x"] <= xmax)
         and (ymin is None or s["y"] >= ymin) and (ymax is None or s["y"] <= ymax)]
    if y is not None: c = sorted([s for s in c if abs(s["y"] - y) <= ytol], key=lambda s: abs(s["y"] - y))
    return c[0] if c else None


def column(spans, xmin, xmax, ymin=None, ymax=None, size=None):
    out = [s for s in spans if xmin <= s["x"] <= xmax and (ymin is None or s["y"] >= ymin) and (ymax is None or s["y"] <= ymax) and (size is None or abs(s["size"] - size) < 0.3)]
    return sorted(out, key=lambda s: s["y"])


def drawings(pdf):
    """get_drawings() on FLAT sheets only (never Figs 1, 3, 4, 5)."""
    d = fitz.open(hydrated(pdf)); p = d[0]
    assert len(p.get_xobjects()) == 0, "get_drawings is only allowed on a flat sheet"
    items = []
    for it in p.get_drawings():
        r = it["rect"]; kinds = "".join(sorted(set(k[0] for k in it["items"])))
        rec = {"rect": [round(r.x0, 3), round(r.y0, 3), round(r.x1, 3), round(r.y1, 3)], "type": it.get("type"), "fill": it.get("fill"), "color": it.get("color"),
               "width": it.get("width"), "dashes": it.get("dashes"), "kinds": kinds, "n": len(it["items"])}
        if kinds == "l" and len(it["items"]) == 1:
            p0, p1 = it["items"][0][1], it["items"][0][2]; rec["line"] = [round(p0.x, 3), round(p0.y, 3), round(p1.x, 3), round(p1.y, 3)]
        items.append(rec)
    d.close(); return items


def markers(items, xmax=None, ymin=None, ymax=None, side_min=3.0, side_max=10.0):
    cl = []
    for it in items:
        x0, y0, x1, y1 = it["rect"]; w, h = x1 - x0, y1 - y0
        if not (side_min <= w <= side_max and side_min <= h <= side_max and abs(w - h) < 1.5): continue
        cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
        if (xmax is not None and cx > xmax) or (ymin is not None and cy < ymin) or (ymax is not None and cy > ymax): continue
        if "c" in it["kinds"]: shape = "circle"
        elif it["kinds"] == "re" or (it["kinds"] == "l" and it["n"] == 4): shape = "square"
        elif it["kinds"] == "l" and it["n"] in (2, 3): shape = "triangle"
        else: continue
        fill, stroke = hexcol(it["fill"]), hexcol(it["color"])
        for c in cl:
            if c["shape"] == shape and abs(c["cx"] - cx) < 0.35 and abs(c["cy"] - cy) < 0.35:
                c["items"] += 1; c["fills"].append(fill); c["strokes"].append(stroke); c["widths"].append(it["width"]); break
        else: cl.append({"cx": cx, "cy": cy, "w": w, "h": h, "shape": shape, "items": 1, "fills": [fill], "strokes": [stroke], "widths": [it["width"]]})
    out = []
    for c in cl:
        nf = [f for f in c["fills"] if f not in (None, "#ffffff")]; ns = [s for s in c["strokes"] if s not in (None, "#ffffff")]
        out.append({"cx": round(c["cx"], 3), "cy": round(c["cy"], 3), "w": round(c["w"], 3), "h": round(c["h"], 3), "shape": c["shape"], "open": not nf,
                    "colour": (nf or ns or [None])[0], "n_items": c["items"], "widths": c["widths"]})
    return out


def hlines(items, min_len=0.5):
    out = []
    for it in items:
        if "line" in it:
            x0, y0, x1, y1 = it["line"]
            if abs(y0 - y1) < 0.05 and abs(x1 - x0) >= min_len:
                out.append({"x0": round(min(x0, x1), 3), "x1": round(max(x0, x1), 3), "y": round(y0, 3), "width": it["width"], "stroke": hexcol(it["color"]), "dashes": it.get("dashes")})
    return out


def vlines(items, max_len=6.0):
    out = []
    for it in items:
        if "line" in it:
            x0, y0, x1, y1 = it["line"]
            if abs(x0 - x1) < 0.05 and 0.5 <= abs(y1 - y0) <= max_len:
                out.append({"x": round(x0, 3), "y0": round(min(y0, y1), 3), "y1": round(max(y0, y1), 3), "stroke": hexcol(it["color"]), "width": it["width"], "dashes": it.get("dashes")})
    return out


def rects(items):
    out = []
    for it in items:
        if it["kinds"] == "re" and it["n"] == 1:
            x0, y0, x1, y1 = it["rect"]; out.append({"x0": x0, "y0": y0, "x1": x1, "y1": y1, "fill": hexcol(it["fill"]), "stroke": hexcol(it["color"]), "w": round(x1 - x0, 3), "h": round(y1 - y0, 3)})
    return out


def tick_xs(items, y_lo, y_hi, length=(2.0, 4.0)):
    xs = sorted(v["x"] for v in vlines(items, max_len=length[1]) if length[0] <= (v["y1"] - v["y0"]) <= length[1] and y_lo <= v["y0"] <= y_hi)
    out = []
    for x in xs:
        if not out or abs(x - out[-1]) > 0.05: out.append(x)
    return out


def fit_axis(values, xs, logscale=True):
    f = np.log(values) if logscale else np.asarray(values, float)
    A = np.vstack([np.ones_like(f), f]).T
    (a, b), *_ = np.linalg.lstsq(A, np.asarray(xs, float), rcond=None)
    res = np.asarray(xs, float) - (a + b * f)
    return float(a), float(b), res


def axis_x(a, b, v, logscale=True): return a + b * (log(v) if logscale else v)


# ------------------------------------------------------------------ drawing on a fresh page: fitz.Shape + fitz.TextWriter (Arial faces)
class Sheet:
    """A single-page PDF drawn from scratch at the V13 page size. Every text goes through text() and is recorded as a drawn record
    (v14lib.DRAWN_RECORD) for the printed-values csv; every path through Shape helpers. Fonts: the system Arial and Arial Bold
    (TextWriter), the ToUnicode hyphen entry corrected to U+002D after saving (the ROUND31 LFOREST idiom)."""
    def __init__(self, width, height):
        self.doc = fitz.open(); self.page = self.doc.new_page(width=width, height=height)
        self.font = fitz.Font(fontfile=ARIAL); self.bold = fitz.Font(fontfile=ARIAL_BOLD)
        self.writers = {}; self.drawn = []; self.shape = fitz.Shape(self.page)
        # a white page background rectangle as the matplotlib sheets carry (the raster proof wants opaque white)
        self.shape.draw_rect(self.page.rect); self.shape.finish(color=None, fill=(1, 1, 1)); self.shape.commit()

    BUMP = 1.0   # round 49: every text one point larger than the size the builder asks for (placement rules see the larger text)

    def width(self, s, size, bold=False, bump=True): return (self.bold if bold else self.font).text_length(s, fontsize=size + (self.BUMP if bump else 0.0))

    def text(self, x, y, s, size, bold=False, color=INK, align="left", panel="", source_file="static:V13", source_key="", source_value="", rule="text", note=""):
        size_asked = size; size = size + self.BUMP
        w = self.width(s, size, bold, bump=False); x0 = x - w if align == "right" else (x - w / 2 if align == "center" else x)
        key = (color, bold)
        if key not in self.writers: self.writers[key] = fitz.TextWriter(self.page.rect)
        self.writers[key].append(fitz.Point(x0, y), s, font=self.bold if bold else self.font, fontsize=size)
        self.drawn.append(dict(text=s, x=round(x, 3), baseline=round(y, 3), ha=align, size=size, bold=bold, color=color, panel=panel,
                               source_file=source_file, source_key=source_key, source_value=source_value, rule=rule, note=note, x_left=round(x0, 3), width=round(w, 3), size_asked=size_asked))
        return x0, w

    def line(self, x0, y0, x1, y1, color, width, cap=1, dashes=None):
        self.shape.draw_line(fitz.Point(x0, y0), fitz.Point(x1, y1)); self.shape.finish(color=rgb(color), width=width, lineCap=cap, dashes=dashes, closePath=False)

    def circle(self, cx, cy, diameter, fill, edge, edge_w):
        self.shape.draw_circle(fitz.Point(cx, cy), diameter / 2); self.shape.finish(color=rgb(edge) if edge else None, fill=rgb(fill) if fill else None, width=edge_w)

    def square(self, cx, cy, side, fill, edge, edge_w):
        self.shape.draw_rect(fitz.Rect(cx - side / 2, cy - side / 2, cx + side / 2, cy + side / 2)); self.shape.finish(color=rgb(edge) if edge else None, fill=rgb(fill) if fill else None, width=edge_w, closePath=True)

    def triangle(self, cx, cy, w, h, fill, edge, edge_w):
        self.shape.draw_polyline([fitz.Point(cx, cy - h / 2), fitz.Point(cx - w / 2, cy + h / 2), fitz.Point(cx + w / 2, cy + h / 2)]); self.shape.finish(color=rgb(edge) if edge else None, fill=rgb(fill) if fill else None, width=edge_w, closePath=True)

    def rect(self, x0, y0, x1, y1, fill=None, edge=None, edge_w=0.0):
        self.shape.draw_rect(fitz.Rect(x0, y0, x1, y1)); self.shape.finish(color=rgb(edge) if edge else None, fill=rgb(fill) if fill else None, width=edge_w, closePath=True)

    def save(self, path):
        self.shape.commit()
        for (color, bold), tw in self.writers.items(): tw.write_text(self.page, color=rgb(color))
        tmp = path + ".tmp.pdf"; self.doc.save(tmp, garbage=1, deflate=True); self.doc.close()
        fix_tounicode(tmp, path); os.remove(tmp); return path


def fix_tounicode(src, dst):
    """MuPDF's reverse glyph map gives the hyphen glyph U+00AD in the ToUnicode CMap of TextWriter fonts: set it to U+002D."""
    d = fitz.open(src); p = d[0]; fixed = []; seen = set()
    for f in p.get_fonts(full=True):
        if f[3] in ("Arial Regular", "Arial Bold", "ArialMT", "Arial-BoldMT") and f[0] not in seen:
            seen.add(f[0]); ff = ARIAL_BOLD if "Bold" in f[3] else ARIAL; gid = fitz.Font(fontfile=ff).has_glyph(0x2D)
            tu = d.xref_get_key(f[0], "ToUnicode")
            if tu[0] != "xref": continue
            tux = int(re.match(r"(\d+) 0 R", tu[1]).group(1)); cm = d.xref_stream(tux)
            old_e = b"<%04x> <00ad>" % gid
            if old_e in cm: d.update_stream(tux, cm.replace(old_e, b"<%04x> <002d>" % gid)); fixed.append((f[3], gid))
    d.save(dst, garbage=1, deflate=True); d.close(); return fixed


# ------------------------------------------------------------------ renders, crops
def render(pdf, dpi, out_png, heavy=None, clip=None): return v14_render(pdf, dpi, out_png, heavy=heavy, clip=clip)


def crop_pair(old_png, new_png, rect_pt, dpi, out_stem):
    S = dpi / 72.0
    x0, y0, x1, y1 = [int(round(v * S)) for v in rect_pt]
    o = Image.open(old_png).convert("RGB"); n = Image.open(new_png).convert("RGB")
    o = o.crop((x0, y0, min(x1, o.width), min(y1, o.height))); n = n.crop((x0, y0, min(x1, n.width), min(y1, n.height)))
    o.save(out_stem + "_OLD.png"); n.save(out_stem + "_NEW.png")
    st = Image.new("RGB", (max(o.width, n.width), o.height + n.height + 6), (255, 0, 0)); st.paste(o, (0, 0)); st.paste(n, (0, o.height + 6)); st.save(out_stem + "_OLD_over_NEW.png")
    return out_stem + "_OLD_over_NEW.png"


def dump(obj, path):
    json.dump(obj, open(path, "w"), indent=1, default=lambda o: float(o) if isinstance(o, (np.floating,)) else (int(o) if isinstance(o, np.integer) else (bool(o) if isinstance(o, np.bool_) else str(o))))
    return path


def dramatic(kind, old, new):
    """A conclusion changes, q or P crosses 0.05, a hazard ratio crosses 1, a rank leaves or enters the top 10, a count moves more
    than 20 percent, a row enters or leaves. Computed, never assumed."""
    o, n = norm(str(old)), norm(str(new))
    if o == n: return ""
    try:
        if kind in ("p", "q"):
            fo = 0.0005 if o.startswith("<") else float(o); fn = 0.0005 if n.startswith("<") else float(n)
            return "q or P crosses 0.05" if (fo < 0.05) != (fn < 0.05) else ""
        if kind == "hr":
            fo, fn = float(o.split()[0]), float(n.split()[0]); lo_o = float(o.split("(")[1].split("-")[0]); lo_n = float(n.split("(")[1].split("-")[0])
            hi_o = float(o.split("-")[-1].rstrip(")")); hi_n = float(n.split("-")[-1].rstrip(")"))
            return "hazard ratio crosses 1" if (fo > 1) != (fn > 1) else ("CI crosses 1" if ((lo_o <= 1 <= hi_o) != (lo_n <= 1 <= hi_n)) else "")
        if kind == "count":
            fo, fn = float(re.sub(r"[^0-9.]", "", o)), float(re.sub(r"[^0-9.]", "", n))
            return "count moves more than 20 percent" if fo and abs(fn - fo) / fo > 0.2 else ""
        if kind == "pct":
            fo, fn = float(o.rstrip("%")), float(n.rstrip("%"))
            return "sign flips" if (fo > 0) != (fn > 0) else ("count moves more than 20 percent" if fo and abs(fn - fo) / abs(fo) > 0.2 else "")
        if kind == "row": return "row enters or leaves"
    except Exception:
        return "changed (unparsed)"
    return ""


# ------------------------------------------------------------------ generic replay of V13 drawing items (key glyphs, static furniture)
def replay_items(S, items, x0, y0, x1, y1, dy=0.0, skip_white_only=False):
    """Redraw every V13 drawing item whose rect lies inside the region, in the V13 order: circles, single lines (round caps),
    filled polylines (triangles from the bbox), rectangles. dy shifts the copy vertically. Returns the number replayed."""
    n = 0
    for it in items:
        r = it["rect"]
        if not (r[0] >= x0 - 0.5 and r[1] >= y0 - 0.5 and r[2] <= x1 + 0.5 and r[3] <= y1 + 0.5): continue
        fill, stroke, w = hexcol(it["fill"]), hexcol(it["color"]), it["width"] or 0.0
        if skip_white_only and fill == WHITE and stroke in (None, WHITE): continue
        if "c" in it["kinds"]:
            S.shape.draw_circle(fitz.Point((r[0] + r[2]) / 2, (r[1] + r[3]) / 2 + dy), (r[2] - r[0]) / 2); S.shape.finish(color=rgb(stroke) if stroke else None, fill=rgb(fill) if fill else None, width=w)
        elif "line" in it:
            l = it["line"]; S.shape.draw_line(fitz.Point(l[0], l[1] + dy), fitz.Point(l[2], l[3] + dy)); S.shape.finish(color=rgb(stroke) if stroke else None, width=w, lineCap=1, closePath=False)
        elif it["kinds"] == "l" and it["n"] in (2, 3):
            cx, cy = (r[0] + r[2]) / 2, (r[1] + r[3]) / 2 + dy; ww, hh = r[2] - r[0], r[3] - r[1]
            S.shape.draw_polyline([fitz.Point(cx, cy - hh / 2), fitz.Point(cx - ww / 2, cy + hh / 2), fitz.Point(cx + ww / 2, cy + hh / 2)]); S.shape.finish(color=rgb(stroke) if stroke else None, fill=rgb(fill) if fill else None, width=w, closePath=True)
        elif "re" in it["kinds"]:
            S.shape.draw_rect(fitz.Rect(r[0], r[1] + dy, r[2], r[3] + dy)); S.shape.finish(color=rgb(stroke) if stroke else None, fill=rgb(fill) if fill else None, width=w, closePath=True)
        else: continue
        n += 1
    return n


def static_texts(S, spans, region, dy=0.0, panel="", exclude=()):
    """Copy the V13 static strings inside a region (x0, y0, x1, y1) at their origins (left-aligned as V13 draws them), bold where V13 is bold."""
    n = 0
    for s in spans:
        if s["text"] in exclude: continue
        if region[0] <= s["x"] <= region[2] and region[1] <= s["y"] <= region[3]:
            S.text(s["x"], s["y"] + dy, s["text"], s["size"], bold="Bold" in s["font"], color="#" + s["color"], panel=panel, source_key="static V13 string"); n += 1
    return n


# ------------------------------------------------------------------ FIX A1 (2026-09-15): the ADJUSTED PAP comparison (owner's decision 6)
ADJ_NOTE = "adjusted for prevalent cardiopulmonary disease (treatment_v3.json step 211 / L4_fig5b_refit_v3.json step 212, decision 6)"


def treatment_rows_adj():
    """ED 8 (and the Fig 5b row rule): the 49 outcomes of treatment_v3.json with the ADJUSTED hazard ratio (adj_hr, adj_lo, adj_hi,
    adj_p, adj_pct_lower_if_corrected) exposed under the field names the builders print, sorted by the adjusted hazard ratio
    descending. Returns (T, rows) where T carries a corrected_vs_not view (n, n_corrected, n_not, outcomes) plus the file's summary."""
    T3 = load_json(SPEC["treatment_v3"]); oc = T3["outcomes"]; c = T3["counts"]
    adj = {}
    for k, v in oc.items():
        adj[k] = dict(key=v["key"], negative_control=bool(v["negative_control"]), events=v["events"], n=v["n"], hr=v["adj_hr"], lo=v["adj_lo"], hi=v["adj_hi"], p=v["adj_p"],
                      q=v.get("adj_q"), pct_lower_if_corrected=v["adj_pct_lower_if_corrected"], unadj_hr=v["hr"], unadj_p=v["p"], unadj_pct=v["pct_lower_if_corrected"])
    T = {"corrected_vs_not": {"n": c["n"], "n_corrected": c["n_corrected"], "n_not": c["n_not"], "outcomes": adj}, "summary": T3["summary"], "counts": c, "covariate": T3["covariate"]}
    return T, sorted(adj.items(), key=lambda kv: -kv[1]["hr"])


def fig5b_rows_adj():
    """Fig 5b on the ADJUSTED comparison: rows = the conditions with adj_p < 0.05 in the two-group comparison of treatment_v3.json
    (the builder's rule, in_current_5b, applied to the adjusted values) and the controls the refit file carries, both by the adjusted
    published hazard ratio (adj_pub_hr) descending; the two band series exposed as B_*/A_* from adj_B_*/adj_A_*."""
    T, rows = treatment_rows_adj(); sig_set = {k for k, v in rows if v["p"] < 0.05 and not v["negative_control"]}
    R = load_json(SPEC["fig5b_refit_v3_json"]); OC = R["outcomes"]
    def view(v):
        w = dict(v)
        for band in ("A", "B"):
            for f in ("hr", "lo", "hi", "p"): w[f"{band}_{f}"] = v[f"adj_{band}_{f}"]
            w[f"{band}_q"] = v.get(f"adj_{band}_q")
        w["pub_hr"] = v["adj_pub_hr"]; w["in_current_5b_adj"] = (not v["negative_control"]) and (w.get("outcome_label") in sig_set)
        return w
    OCV = {k: view(dict(v, outcome_label=k)) for k, v in OC.items()}
    SIG = sorted(((k, v) for k, v in OCV.items() if k in sig_set), key=lambda kv: -kv[1]["pub_hr"])
    CTRL = sorted(((k, v) for k, v in OCV.items() if v["negative_control"]), key=lambda kv: -kv[1]["pub_hr"])
    missing = sorted(sig_set - set(OC)); assert not missing, f"adjusted-significant conditions absent from the refit file: {missing}"
    return R, SIG, CTRL, sig_set


def organ_rollup_adj():
    """ED 10 roll-up on the ADJUSTED comparison: the organ_rollup() rule (controls held out, composite and mortality excluded,
    Responding = adj_p < .05 and adj_hr > 1, median of the adjusted per-outcome percentage lower where restored)."""
    sys.path.insert(0, SPEC["disease_definitions_dir"])
    import disease_definitions as dd
    NAME = {k: v[0] for k, v in dd.DISEASES.items()}
    GROUP_BY_NAME = {v[0]: dd.ORGAN_GROUP[k] for k, v in dd.DISEASES.items()}
    GROUP_BY_NAME["Death from any cause"] = "Mortality"
    NEG_NAMES = {NAME[k] for k in dd.NEGATIVE_CONTROLS}
    T, rows = treatment_rows_adj(); cv = T["corrected_vs_not"]["outcomes"]
    unmapped = [c for c in cv if c not in GROUP_BY_NAME]
    assert not unmapped, f"outcomes without an organ group in disease_definitions.py: {unmapped}"
    grp = []
    for name in sorted(set(GROUP_BY_NAME.values())):
        if name in ("Composite", "Mortality"): continue
        negs = [c for c in cv if GROUP_BY_NAME.get(c) == name and c in NEG_NAMES]
        real = [(c, v) for c, v in cv.items() if GROUP_BY_NAME.get(c) == name and c not in NEG_NAMES]
        if not real: continue
        grp.append({"group": name, "conditions": len(real), "controls_held_out": len(negs),
                    "responding": sum(1 for _, v in real if v["p"] < .05 and v["hr"] > 1),
                    "median_pct_lower": float(np.median([v["pct_lower_if_corrected"] for _, v in real])),
                    "members": [c for c, _ in real], "control_members": negs})
    grp.sort(key=lambda g: -g["median_pct_lower"])
    ungrouped = [c for c in cv if GROUP_BY_NAME.get(c) in ("Composite", "Mortality")]
    mapping = {c: GROUP_BY_NAME[c] for c in cv}
    return grp, ungrouped, mapping


def unadjusted_v14_spans(sheet):
    """The text layer of the unadjusted V14 sheet snapshotted before FIX A1 (work/<Sheet>_UNADJUSTED_V14.pdf): the 'before' side of CHANGES."""
    sp, info = text_spans(f"{LANE}/{sheet}/work/{sheet}_UNADJUSTED_V14.pdf"); return collapse(sp), info
