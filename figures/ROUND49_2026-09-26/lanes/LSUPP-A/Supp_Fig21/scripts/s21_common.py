#!/usr/bin/env python3
"""Supp_Fig21, V14 lane L2 (round 37): shared paths and helpers (repointed copy of the round-30 module, byte copy beside as
s21_common_PRE_V8_1.py). Base = the V13 sheet (the round-30 re-plot with the round-31 LTEXT wording "T90 gain" / "90% threshold").
Numbers (v8.1): Sleep_Variability_2026-08/t_ladder_up/LADDER_SUMMARY.csv + ladder_peak.json (step 202), ladder_persd.csv (step 200),
numbers/ranking_v3.csv + ranking_v3_percondition.csv (step 105). OLD side (positive control and the delta only): the v7 snapshot
OLD = V8_RECALC_2026-09-12/compare/old_numbers (never a live source).
PANEL B COUNT RULE (decided in this lane, reported): the axis title names the ranked outcomes ("Outcomes significant per 1 SD (of N)",
N = ranking_v3.csv nout), so the drawn counts are the ranked outcomes with BH q < 0.05 per threshold, read from ladder_persd.csv (the
file's own q) restricted to the 52 ranked keys of ranking_v3_percondition.csv. LADDER_SUMMARY's column n_fdr_sig_of_48_persd counts
the 53 non-control outcomes instead (it includes the four circular outcomes that are not ranked and excludes the three ranked
controls), so it differs by one or two at four rungs; both are written to the CHANGES csv for the owner."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import os, json, hashlib, re, math, sys
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND37_2026-09-14/figures/_common")
from v14lib import *
T90 = paths.FIGURE_ROOT
SVL = f"{SV}/t_ladder_up"
LANE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORK, VER, CROPS = f"{LANE}/work", f"{LANE}/verify", f"{LANE}/verify/crops"
for _d in (WORK, VER, CROPS): os.makedirs(_d, exist_ok=True)
BASE = f"{V13}/Supp_Fig21.pdf"
LOCK = LOCK_ROOT
LADDER_CSV = f"{SVL}/LADDER_SUMMARY.csv"; PEAK_JSON = f"{SVL}/ladder_peak.json"; PERSD_CSV = f"{SVL}/ladder_persd.csv"
RANKING_CSV = f"{NUM}/ranking_v3.csv"; PERCOND_CSV = f"{NUM}/ranking_v3_percondition.csv"
OLD_LADDER_CSV = f"{paths.OLD_SNAPSHOT_SV}/t_ladder_up/LADDER_SUMMARY.csv"
OLD_RANKING_CSV = f"{OLD_V7}/numbers/ranking_v3.csv"
FONT_REG = paths.FONT_ARIAL; FONT_BOLD = paths.FONT_ARIAL_BOLD; FONT_ITAL = paths.FONT_ARIAL_ITALIC
RUNGS = ["T80", "T85", "T88", "T90", "T92", "T95"]; PERSD_COL = {"T80": "t80", "T85": "t85", "T88": "t88_frozen", "T90": "t90_frozen", "T92": "t92", "T95": "t95"}
INK, BLUE, WHITE = "#1a1d21", "#0288d1", "#ffffff"
GAIN_LABEL, KEY_LABEL = "T90 gain", "90% threshold"          # the V13 (round-31) wording, no "published"


def sidecar(path):
    """v8.1 gate through v14lib (data_frozen_v8 inputs or a v8.1 STATUS row), plus the sha256 of the file read now."""
    d = sidecar_v8(path); return dict(file=path, step=d["step"], name=d["name"], ended=d["status_ok"], output_mtime=d["output_mtime"], sha256=d["sha256"], v8_inputs=d["v8_inputs"])


def ranked_keys():
    import pandas as pd
    pc = pd.read_csv(hydrated(PERCOND_CSV), comment="#"); return sorted(set(pc.outcome))


def ranked_counts(persd_csv, keys):
    """Ranked outcomes with q < 0.05 per rung, from the per-SD file's own q columns restricted to the ranked keys."""
    import pandas as pd
    P = pd.read_csv(hydrated(persd_csv)); sub = P[P.key.isin(keys)]
    assert len(sub) == len(keys), (len(sub), len(keys), sorted(set(keys) - set(P.key)))
    return {r: int((sub[f"{PERSD_COL[r]}_q"] < 0.05).sum()) for r in RUNGS}, {r: int((P[~P.negative_control][f"{PERSD_COL[r]}_q"] < 0.05).sum()) for r in RUNGS}, int((~P.negative_control).sum())


def hexcol(c):
    return None if c is None else "#%02x%02x%02x" % tuple(int(round(v * 255)) for v in c[:3])


def spans_of(pdf_path):
    import fitz, gc
    doc = fitz.open(hydrated(pdf_path)); pg = doc[0]; out = []
    for b in pg.get_text("rawdict")["blocks"]:
        if b["type"] != 0: continue
        for l in b["lines"]:
            for s in l["spans"]:
                txt = "".join(c["c"] for c in s["chars"])
                out.append(dict(text=txt, origin=[round(v, 3) for v in s["origin"]], bbox=[round(v, 3) for v in s["bbox"]], size=round(s["size"], 3), font=s["font"],
                                color="#%06x" % s["color"], dir=[round(v, 3) for v in l["dir"]], asc=round(s["ascender"], 4), desc=round(s["descender"], 4)))
    words = [dict(text=w[4], bbox=[round(v, 3) for v in w[:4]]) for w in pg.get_text("words")]
    fonts = [dict(xref=f[0], type=f[2], name=f[3], enc=f[5]) for f in pg.get_fonts(full=True)]
    info = dict(rect=[round(v, 4) for v in pg.rect], mediabox=[round(v, 4) for v in pg.mediabox], cropbox=[round(v, 4) for v in pg.cropbox], n_pages=len(doc), n_xobjects=len(pg.get_xobjects()), fonts=fonts)
    doc.close(); gc.collect(); return out, words, info


def collapse_spans(spans, tol=0.25):
    kept = []
    for s in spans:
        if not any(k["text"] == s["text"] and abs(k["origin"][0] - s["origin"][0]) <= tol and abs(k["origin"][1] - s["origin"][1]) <= tol for k in kept): kept.append(s)
    return kept


def drawings_of(pdf_path):
    import fitz, gc
    doc = fitz.open(hydrated(pdf_path)); pg = doc[0]; out = []
    for g in pg.get_drawings():
        r = g["rect"]
        out.append(dict(type=g["type"], rect=[round(v, 3) for v in r], fill=hexcol(g["fill"]), stroke=hexcol(g["color"]), width=None if g.get("width") is None else round(g["width"], 3),
                        dashes=g.get("dashes"), n_items=len(g["items"]), kinds="".join(sorted({it[0] for it in g["items"]}))))
    doc.close(); gc.collect(); return out


NUMTOK = re.compile(r"\d+(?:\.\d+)?")
def numeric_tokens(words):
    toks = []
    for w in words:
        t = w["text"]
        if re.fullmatch(r"[\(\[]?\d+(?:\.\d+)?[\)\]%,.]*", t): toks.append(NUMTOK.search(t).group(0))
    return toks


def word_tokens(words): return [w["text"] for w in words]
