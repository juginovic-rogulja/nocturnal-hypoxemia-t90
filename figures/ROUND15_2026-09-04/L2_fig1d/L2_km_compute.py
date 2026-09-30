#!/usr/bin/env python3
"""ROUND 15, lane L2 (ruling 5): the eight Kaplan-Meier panels of Figure 1d.

Every curve and every log-rank P on the rebuilt panel d comes from this file, computed on
the frozen cohort. Nothing is copied off the round-14 sheet.

Method is taken unchanged from the script that produced the four panels now on the page,
FINAL_FIGURES_2026-08-14/_scripts/figure1D_polish.py:

  cohort      numbers/cohort_spec.apply_cohort on data_frozen_v7_2026-09/t90_final.parquet, 19,173
  bands       T90 <=1% / 1-10% / >10%, three bands, from spo2_pct_below_90
  risk set    per outcome, prevalent cases dropped, follow-up time present and > 0
  curves      crude Kaplan-Meier cumulative incidence, no adjustment, clipped at 8 years
  test        lifelines multivariate_logrank_test on the three bands, 2 df
  P string    Nature form, floored at P < 0.001, round-trip asserted against the value

PART A runs first; since the v8 F2 fix (2026-09-12) it is a logged OLD-versus-NEW comparison, not a gate. It
recomputes the four panels already drawn (cirrhosis, heart failure, type 2 diabetes,
hypertension) and checks them against the round-14 page three ways:
  1  the printed log-rank string of each tile, read out of the PDF text layer
  2  the y tick label set of each tile, read out of the PDF text layer
  3  the plotted level of all three curves near 8 years, read out of the RENDERED PIXELS
     of the round-14 page and converted to data units with that tile's own tick spacing
Check 3 also pins the axes geometry: ylim = 1.14 * ymax, so each tile independently
implies an axes box height, and all four must agree.

PART B computes the four new panels. They are the four highest ranked conditions of
panel c on the same sheet that row one does not already draw and that admit a clean
4 to 6 major y axis under the sheet's tick rule.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import warnings; warnings.filterwarnings("ignore")
import json, sys
import numpy as np
import pandas as pd
import fitz
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test
from scipy.stats import chi2

ROOT = paths.FIGURE_ROOT
LANE = f"{paths.FIGURE_ROOT}/ROUND15_2026-09-04/L2_fig1d"
R14 = f"{LANE}/Main_Fig1_R14_SOURCE.pdf"
sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
from cohort_spec import apply_cohort, COHORT_N, T90_FINAL  # noqa: E402  # v8 sweep 2026-09-12; T90_FINAL: F1 fix
from disease_definitions import DISEASES                # noqa: E402

XMAX = 8.0
BAND_LABELS = [("<=1%", 0, "#b3dcf2"), ("1-10%", 1, "#3f9fd8"), (">10%", 2, "#0288d1")]

OLD = [("cirrhosis", "Cirrhosis"), ("hf", "Heart failure"),
       ("diabetes", "Type 2 diabetes"), ("htn2", "Hypertension")]
NEW = [("resp_failure", "Respiratory failure"), ("obesity", "Obesity"),
       ("copd2", "COPD"), ("pneumonia", "Pneumonia")]
# risk set and event count of the four row-one panels, read from the primary analysis table
# numbers/bdsp_diseases_v3.csv (run_primary_unpenalized.py, step 111, provenance sidecar; its n and events follow the
# same risk-set rule as panel() below), never typed. v8 F1 fix (2026-09-12, integrity gate 6): the v7 literals were
# cirrhosis 18813/199, hf 16587/1187, diabetes 14861/1404, htn2 9274/1829; the v8 table moves diabetes to 14894/1413.
_PRIM = pd.read_csv(f"{paths.NUMBERS_DIR}/bdsp_diseases_v3.csv").set_index("key")
EXPECT_OLD = {k: (int(_PRIM.loc[k, "n"]), int(_PRIM.loc[k, "events"])) for k, _ in OLD}

log_lines = []
def say(s=""):
    print(s)
    log_lines.append(s)


# ------------------------------------------------------------------ cohort and bands
# v7 (2026-09-07) ran PART A (the pixel self-test against the round-14 sheet) on the superseded April table the sheet
# was drawn from. v8 F1 fix (2026-09-12, integrity gate 7): a chain step reads DATA_DIR only, so both parts read the
# spec's T90_FINAL; PART A now compares the current table's curves with the round-14 sheet.
V7_TABLE = f"{paths.TABLES_DIR}/t90_final.parquet"
CTRL_TABLE = T90_FINAL      # v8 F1 fix: was the April data_frozen/ table; self-test input, never printed

def load_bands(path):
    b_ = apply_cohort(pd.read_parquet(path))
    assert len(b_) == COHORT_N
    t90_ = b_.spo2_pct_below_90
    assert t90_.notna().all()
    b_ = b_.assign(band=np.where(t90_ <= 1, 0, np.where(t90_ <= 10, 1, 2)))
    bn = {g: int((b_.band == g).sum()) for _, g, _ in BAND_LABELS}
    assert sum(bn.values()) == COHORT_N, bn
    return b_, bn

b_v7, BAND_N_V7 = load_bands(V7_TABLE)
b, BAND_N = load_bands(CTRL_TABLE)
print(f"T90 band sizes: current table {BAND_N_V7}; self-test table {BAND_N}")


def jama_p(p):
    if p < 0.001:
        return "P < 0.001"
    if p < 0.01:
        return f"P = {p:.3f}"
    if p > 0.99:
        return "P > 0.99"
    return f"P = {p:.2f}"


def assert_p_string(s, p):
    if "<" in s:
        assert p < 0.001, (s, p)
    elif ">" in s:
        assert p > 0.99, (s, p)
    else:
        v = float(s.split("=")[1].strip())
        assert round(p, len(s.split(".")[-1])) == v, (s, p)


def km_curve(T, E):
    km = KaplanMeierFitter().fit(T, E)
    cd = km.cumulative_density_
    tt = np.asarray(cd.index, dtype=float)
    yy = np.asarray(cd.iloc[:, 0], dtype=float) * 100.0
    keep = tt <= XMAX
    tt, yy = tt[keep], yy[keep]
    if tt[-1] < XMAX:
        tt = np.append(tt, XMAX)
        yy = np.append(yy, yy[-1])
    return tt, yy


def clean_yticks(ymax_lim, n_max=6):
    """4 to 6 clean majors from 0, largest at or under the axis top. From figure1D_polish."""
    for step in (1, 2, 5, 3, 10, 20):   # v8: step 3 after 5, only when 5 leaves fewer than 4 majors (axis tops in [12, 15) and [18, 19))
        top = int(np.floor(ymax_lim / step)) * step
        n = top // step + 1
        if 4 <= n <= n_max:   # v8: the docstring's 4-to-6 rule, enforced here instead of by the caller's assert
            return list(range(0, top + 1, step))
    raise AssertionError(ymax_lim)


def panel(key, label, expect=None):
    assert DISEASES[key][0] == label, (key, DISEASES[key][0])
    f = b[(b[f"{key}_prevalent"] == 0) & b[f"{key}_years"].notna() & (b[f"{key}_years"] > 0)]
    n, ev = len(f), int(f[f"{key}_incident"].sum())
    if expect is not None:
        assert (n, ev) == expect, (key, n, ev, expect)
    T = f[f"{key}_years"].values
    E = f[f"{key}_incident"].astype(int).values
    G = f.band.values
    res = multivariate_logrank_test(T, G, E)
    stat, p = float(res.test_statistic), float(res.p_value)
    pc = float(chi2.sf(stat, 2))
    assert abs(p - pc) <= 1e-12 + 1e-9 * max(p, pc), (p, pc)
    pstr = jama_p(p)
    assert_p_string(pstr, p)
    rows = []
    for blab, g, col in BAND_LABELS:
        m = G == g
        assert float(T[m].max()) >= XMAX, (key, blab, "risk set empty before the display edge")
        tt, yy = km_curve(T[m], E[m])
        rows.append({"band": blab, "colour": col, "n": int(m.sum()), "events": int(E[m].sum()),
                     "at_risk": [int((T[m] >= q).sum()) for q in (0, 2, 4, 6, 8)],
                     "ci8_pct": float(yy[-1]), "t": tt, "y": yy})
    assert sum(r["n"] for r in rows) == n
    assert sum(r["events"] for r in rows) == ev
    assert rows[0]["ci8_pct"] < rows[1]["ci8_pct"] < rows[2]["ci8_pct"], (
        key, [r["ci8_pct"] for r in rows])
    ymax = max(r["y"].max() for r in rows)
    ylim = ymax * 1.14
    yt = clean_yticks(ylim)
    assert 4 <= len(yt) <= 6, (key, yt)
    return {"key": key, "label": label, "n": n, "events": ev, "chi2": stat, "p": p,
            "p_printed": pstr, "ymax": float(ymax), "ylim": float(ylim), "yticks": yt,
            "bands": rows}


# ================================================================== PART A, self-test
say("=" * 78)
say("PART A  SELF-TEST: recompute the four panels already on the round-14 sheet")
say("=" * 78)
say(f"cohort {len(b):,}   bands  <=1% {BAND_N[0]:,}   1-10% {BAND_N[1]:,}   >10% {BAND_N[2]:,}")
say()

old_panels = [panel(k, lab, EXPECT_OLD[k]) for k, lab in OLD]

# ---- the round-14 page: text layer and pixels
doc = fitz.open(R14)
page = doc[0]
SPINE = {"cirrhosis": 528.84, "hf": 641.07, "diabetes": 753.29, "htn2": 865.52}
AXW, BOTTOM = 83.70, 726.97
YLAB_RIGHT_GAP, YLAB_BASE_DY = 4.84, 2.66
ASC = 0.905                              # Arial ascent fraction, baseline = bbox.y0 + ASC*size

spans = []
for bl in page.get_text("dict")["blocks"]:
    if bl["type"] != 0:
        continue
    for ln in bl["lines"]:
        for sp in ln["spans"]:
            spans.append(dict(text=sp["text"], size=round(sp["size"], 2), bbox=sp["bbox"],
                              origin=sp["origin"]))

S = 12                                    # 864 dpi readback of the drawn curves
clipr = fitz.Rect(490, 595, 960, 745)
pm = page.get_pixmap(matrix=fitz.Matrix(S, S), clip=clipr)
img = np.frombuffer(pm.samples, dtype=np.uint8).reshape(pm.height, pm.width, 3).astype(int)

def read_curve_y(xpt, rgb):
    """centre y, in page points, of the run of `rgb` pixels in the column at x = xpt."""
    cx = int(round((xpt - clipr.x0) * S))
    col = np.abs(img[:, cx, :] - np.array(rgb)).sum(axis=1) < 40
    runs, st = [], None
    for i, v in enumerate(col):
        if v and st is None:
            st = i
        if not v and st is not None:
            runs.append((st, i)); st = None
    if st is not None:
        runs.append((st, len(col)))
    if not runs:
        return None
    st, en = max(runs, key=lambda r: r[1] - r[0])
    return clipr.y0 + (st + en) / 2.0 / S

fails, checks = [], 0
boxheights = []
# v8 F2 fix 2026-09-12 (rule 11): PART A is a logged OLD-versus-NEW comparison, not a gate. The round-14 sheet
# (Main_Fig1_R14_SOURCE.pdf) was drawn from the April table; the v8 numbers are compared with it here, every check
# printed with old, new and difference and written to km_pixel_oldnew_v8.csv, and the run never exits on a
# difference: the sheet is redrawn from the v8 numbers in the recompose phase. The machinery asserts stay (P-string
# round trip, tick rule, counts against bdsp_diseases_v3.csv, ordered bands, a log-rank span found on the page).
oldnew = []


def record(panel_key, check, old, new, tol, unit=""):
    """old = what the round-14 sheet shows, new = the v8 computation; tol None means an exact string comparison."""
    global checks
    checks += 1
    diff = None
    if tol is not None and old is not None and new is not None:
        diff = float(new) - float(old)
    if tol is None:
        differs = str(old) != str(new)
    else:
        differs = diff is None or abs(diff) > tol
    oldnew.append({"panel": panel_key, "check": check, "old_round14_sheet": old, "new_v8": new,
                   "difference_new_minus_old": None if diff is None else round(diff, 4), "tolerance": tol,
                   "unit": unit, "differs": bool(differs)})
    say(f"  {'DIFF' if differs else 'same'}  {panel_key:10s} {check:40s} old {old!s:>20} new {new!s:>20} "
        f"diff {'n/a' if diff is None else f'{diff:+.4f}'} {unit}")
    if differs:
        fails.append(f"{panel_key} {check}: old {old!r}, new {new!r}")
    return differs


for pn in old_panels:
    key = pn["key"]
    spine = SPINE[key]

    # ---- check 1: the printed log-rank string
    hits = [s for s in spans if s["size"] == 9.5 and abs(s["bbox"][0] - spine) < 0.6
            and 615 < s["bbox"][1] < 635]
    assert len(hits) >= 1, (key, "no log-rank span found")
    printed = hits[0]["text"].strip()
    want = "Log-rank " + pn["p_printed"]
    record(key, "log-rank string", printed, want, None)

    # ---- check 2: the y tick label set
    labs = [s for s in spans if s["size"] == 10.0
            and abs(s["bbox"][2] - (spine - YLAB_RIGHT_GAP)) < 1.5 and 620 < s["bbox"][1] < 735]
    labs.sort(key=lambda s: -s["origin"][1])
    got = [int(s["text"]) for s in labs]
    record(key, "y tick set", str(got), str(pn["yticks"]), None)

    # tick spacing from the drawn labels, and the axes box it implies
    ys = [s["origin"][1] - YLAB_BASE_DY for s in labs]
    step = pn["yticks"][1] - pn["yticks"][0]
    ppu = -float(np.mean(np.diff(ys))) / step
    boxh = ppu * pn["ylim"]
    boxheights.append(boxh)
    record(key, "zero tick y (pt) vs bottom spine", round(float(ys[0]), 3), BOTTOM, 0.05, "pt")

    # ---- check 3: the plotted level of each curve near 8 years, read off the pixels
    xpt = spine + AXW * (7.75 / 8.0)
    tq = 7.75
    for r in pn["bands"]:
        rgb = tuple(int(r["colour"][i:i + 2], 16) for i in (1, 3, 5))
        ypt = read_curve_y(xpt, rgb)
        computed = float(np.interp(tq, r["t"], r["y"]))
        if ypt is None:
            record(key, f"curve level at {tq} yr, T90 {r['band']} (pct)", None, round(computed, 3), None, "pct")
            continue
        drawn = (BOTTOM - ypt) / ppu
        r["_readback"] = {"drawn_pct": round(drawn, 3), "computed_pct": round(computed, 3)}
        tol = max(0.06, 0.012 * max(1.0, computed))     # 1.2%, floored at 0.06 percentage points
        record(key, f"curve level at {tq} yr, T90 {r['band']} (pct)", round(drawn, 3), round(computed, 3),
               tol, "pct")

    say(f"{pn['label']:18s} n={pn['n']:6,}  events={pn['events']:5,}  chi2={pn['chi2']:8.2f}  "
        f"{pn['p_printed']:10s}  yticks={pn['yticks']}")
    for r in pn["bands"]:
        rb = r.get("_readback", {})
        say(f"    T90 {r['band']:6s} n={r['n']:5,} ev={r['events']:5,}  8-yr incidence "
            f"{r['ci8_pct']:6.2f}%   readback at 7.75 yr: page {rb.get('drawn_pct')} vs "
            f"computed {rb.get('computed_pct')}")

# the four tiles of the round-14 sheet each imply an axes box; compared with their mean and with each other, not gated
BOX_H = float(np.mean(boxheights))
for pn, bh in zip(old_panels, boxheights):
    record(pn["key"], "axes box height implied by the tile (pt)", round(float(bh), 3), round(BOX_H, 3), 0.05, "pt")
record("all", "axes box height spread across tiles (pt)", round(float(max(boxheights) - min(boxheights)), 4),
       0.0, 0.05, "pt")
BOX_TOP = BOTTOM - BOX_H

say()
say(f"implied axes box: height {BOX_H:.3f} pt, top {BOX_TOP:.3f} pt, bottom {BOTTOM} pt "
    f"(spread across the four tiles {max(boxheights) - min(boxheights):.4f} pt)")
n_diff = sum(1 for r_ in oldnew if r_["differs"])
say(f"checks compared: {checks}   differ: {n_diff}")
for f_ in fails:
    say("  DIFF  " + f_)
say()
OLDNEW_CSV = f"{LANE}/km_pixel_oldnew_v8.csv"
pd.DataFrame(oldnew).to_csv(OLDNEW_CSV, index=False)
say(f"PART A OLD-VERSUS-NEW written -> {OLDNEW_CSV}  ({len(oldnew)} rows, {n_diff} differ). The round-14 sheet was "
    f"drawn from the April table and is redrawn from the v8 numbers in the recompose phase, so a difference here is "
    f"expected and is not a gate (v8 F2 fix 2026-09-12).")
say()
b, BAND_N = b_v7, BAND_N_V7          # v7: everything from here on comes from the corrected table
say(f"switched to the v7 table for PART B: bands  <=1% {BAND_N[0]:,}   1-10% {BAND_N[1]:,}   >10% {BAND_N[2]:,}")
old_panels = [panel(k, lab, EXPECT_OLD[k]) for k, lab in OLD]   # v7: the four existing panels re-drawn from the corrected bands too
for pn in old_panels:
    say(f"{pn['label']:20s} (v7) n={pn['n']:6,}  events={pn['events']:5,}  {pn['p_printed']:10s}  " + "  ".join(f"{r['band']} {r['ci8_pct']:.2f}%" for r in pn["bands"]))

# ================================================================== PART B, new panels
say("=" * 78)
say("PART B  the four NEW panels")
say("=" * 78)
new_panels = [panel(k, lab) for k, lab in NEW]
for pn in new_panels:
    say(f"{pn['label']:20s} n={pn['n']:6,}  events={pn['events']:5,}  chi2={pn['chi2']:8.2f}  "
        f"{pn['p_printed']:10s}  ymax={pn['ymax']:.2f}  yticks={pn['yticks']}")
    for r in pn["bands"]:
        say(f"    T90 {r['band']:6s} n={r['n']:5,} ev={r['events']:5,}  8-yr incidence {r['ci8_pct']:6.2f}%"
            f"  at risk {r['at_risk']}")
say()

# ------------------------------------------------------------------ values JSON
def strip(pn):
    return {k: v for k, v in pn.items() if k != "bands"} | {
        "bands": [{k: v for k, v in r.items() if k not in ("t", "y")} |
                  {"curve_t": [round(float(x), 6) for x in r["t"]],
                   "curve_cum_incidence_pct": [round(float(x), 6) for x in r["y"]]}
                  for r in pn["bands"]]}

out = {
    "lane": "ROUND15 L2_fig1d",
    "source_data": V7_TABLE, "selftest_source_data": CTRL_TABLE,
    "method": "crude Kaplan-Meier cumulative incidence by T90 band, clipped at 8 years, "
              "multivariate log-rank on 3 bands (2 df); identical to figure1D_polish.py",
    "cohort_n": int(len(b)),
    "band_n": {lab: BAND_N[g] for lab, g, _ in BAND_LABELS},
    "band_colours": {lab: c for lab, _, c in BAND_LABELS},
    "xmax_years": XMAX,
    "axes_box_pt": {"width": AXW, "height": round(BOX_H, 3), "top": round(BOX_TOP, 3),
                    "bottom": BOTTOM, "spines": SPINE},
    "selftest": {"mode": "old-versus-new comparison against the round-14 sheet, logged, not a gate (v8 F2 fix)",
                 "checks": checks, "differences_n": n_diff, "failures": fails, "oldnew_csv": OLDNEW_CSV,
                 "box_heights_implied": [round(v, 4) for v in boxheights]},
    "row1_existing": [strip(p) for p in old_panels],
    "row2_new": [strip(p) for p in new_panels],
}
json.dump(out, open(f"{LANE}/L2_fig1d_values.json", "w"), indent=1)
say(f"written -> L2_fig1d_values.json")
open(f"{LANE}/L2_selftest_output.txt", "w").write("\n".join(log_lines) + "\n")
