#!/usr/bin/env python3
"""01_build for Supp_Fig01 (lane V14_L1_RANK, round 37: v8.1 numbers on the V13 design; repointed from the V7_L1_RANK builder). The nine-row cohort sleep profile re-plotted from
numbers/figure1_sleep_profile_v1.json (step 28, rerun 2026-09-08) with matplotlib pinned to the round-30 base sheet's page points
(base_text.json, probe_geom.py records). Every quantile is re-derived from the v7 frozen table through numbers/cohort_spec.apply_cohort
and compared at the stored precision, the efig22_polish.py gate. No row prints a value: the median, interquartile box and 5th to 95th
range are drawn, the row axes carry the three tick values. The sleep-efficiency row's axis moves from 40/70/100 to 20/60/100 because
its 5th percentile is now 37.5 (the builder's rule: the axis must contain the 5th to 95th percentile at round values); every other
row keeps its axis, asserted. Colours as the base sheet shows: #0288d1 for the oxygenation rows, #8a9099 for the rest, #ccd1d6 axes.
Writes work/sheet_raw.pdf, Supp_Fig01.pdf (font metrics aligned to the base), work/placed_texts.json, verify/Supp_Fig01_rows_drawn.json."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import os, sys, json
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C

S = "Supp_Fig01"; SD = f"{C.LANE}/{S}"; WORK = f"{SD}/work"; VER = f"{SD}/verify"
os.makedirs(WORK, exist_ok=True); os.makedirs(VER, exist_ok=True)
plt = C.setup_matplotlib()
base_text = json.load(open(f"{WORK}/base_text.json")); W, H = base_text["page"]
assert [W, H] == [484.72, 369.36], (W, H)
SP = base_text["spans"]

# ------------------------------------------------------------------ numbers, with their sidecars, re-derived from the v7 table
PROF_PATH = f"{C.NUM}/figure1_sleep_profile_v1.json"; PROF_SC = C.sidecar_ok(PROF_PATH, need_today=True); PROF_SHA = C.sha256(PROF_PATH)
PROF = json.load(open(C.hydrated(PROF_PATH)))
COH_PATH = f"{C.NUM}/cohorts.json"; COH_SC = C.sidecar_ok(COH_PATH); COH_SHA = C.sha256(COH_PATH); COH = json.load(open(C.hydrated(COH_PATH)))
sys.path.insert(0, paths.ANALYSIS_DIR)
from cohort_spec import apply_cohort, COHORT_N, DATA_DIR
N_COHORT = int(COH["bdsp"]["n_analysis"]); assert N_COHORT == COHORT_N == PROF["cohort_n"] == 19173
assert PROF["source_file"].endswith("data_frozen_v8_2026-09/t90_final.parquet") and DATA_DIR.endswith("data_frozen_v8_2026-09"), (PROF["source_file"], DATA_DIR)
C_ROWS = ["tst_h", "sleep_efficiency_pct", "n3_pct", "rem_pct", "arousal_index", "ahi", "spo2_mean", "spo2_nadir", "t90"]
C_UNIT = {"tst_h": "hours", "sleep_efficiency_pct": "% of time in bed", "n3_pct": "% of sleep", "rem_pct": "% of sleep", "arousal_index": "per hour", "ahi": "per hour",
          "spo2_mean": "% saturation", "spo2_nadir": "% saturation", "t90": "% of the night"}
# each row's own axis at round values that contain the 5th to 95th percentile (the builder's rule). v8.1: the arousal index (p95 82.8) leaves the
# V13 axis 0/40/80 for 0/50/100 and T90 (p95 42.5) leaves 0/20/40 for 0/25/50; every other row keeps its V13 axis, asserted below.
C_AXIS = {"tst_h": (2, 5, 8), "sleep_efficiency_pct": (20, 60, 100), "n3_pct": (0, 25, 50), "rem_pct": (0, 20, 40), "arousal_index": (0, 50, 100), "ahi": (0, 50, 100),
          "spo2_mean": (88, 94, 100), "spo2_nadir": (70, 85, 100), "t90": (0, 25, 50)}
AXIS_V13 = {"tst_h": (2, 5, 8), "sleep_efficiency_pct": (20, 60, 100), "n3_pct": (0, 25, 50), "rem_pct": (0, 20, 40), "arousal_index": (0, 40, 80), "ahi": (0, 50, 100),
            "spo2_mean": (88, 94, 100), "spo2_nadir": (70, 85, 100), "t90": (0, 20, 40)}                 # the axes the V13 sheet prints (its text layer)
AXIS_OLD = {k: v for k, v in AXIS_V13.items() if C_AXIS[k] != v}
M = PROF["measures"]; assert all(k in M for k in C_ROWS)
T90F = f"{DATA_DIR}/t90_final.parquet"; C.hydrated(T90F)
cols = sorted({M[k]["source_column"] for k in C_ROWS} | {"fu_valid", "oximetry_bad", "spo2_pct_below_90", "BDSPPatientID"})
COHORT = apply_cohort(pd.read_parquet(T90F, columns=cols)); assert len(COHORT) == N_COHORT
def _dp(x):
    s = repr(float(x)); return len(s.split(".")[1].rstrip("0")) if "." in s else 0
REDERIVED = {}
for k in C_ROWS:
    m = M[k]; s = (COHORT["TST_min"].dropna() / 60.0) if k == "tst_h" else COHORT[m["source_column"]].dropna()
    assert len(s) == m["n"], f"{k}: parquet n {len(s)}, stored {m['n']}"
    got = {}
    for nm, q in (("median", .50), ("q1", .25), ("q3", .75), ("p05", .05), ("p95", .95)):
        d = _dp(m[nm]); got[nm] = round(float(s.quantile(q)), d)
        assert got[nm] == round(float(m[nm]), d), (k, nm, got[nm], m[nm])
    REDERIVED[k] = dict(n=int(len(s)), **got)
    lo, mid, hi = C_AXIS[k]
    assert lo <= m["p05"] and m["p95"] <= hi and lo < mid < hi, f"{k}: axis {lo} to {hi} clips {m['p05']}-{m['p95']}"
    if k in AXIS_OLD:
        olo, omid, ohi = AXIS_OLD[k]; assert not (olo <= m["p05"] and m["p95"] <= ohi), f"{k}: the old axis {AXIS_OLD[k]} still contains the 5th to 95th percentile, no change needed"
    else:
        assert C_AXIS[k] == AXIS_V13[k], k
    assert (m["family"] == "Oxygenation") == (k in ("spo2_mean", "spo2_nadir", "t90")), k
    # v8.1: the V13 tick strings of every axis must be on the base sheet (the axes the sheet prints are the ones AXIS_V13 declares)
    for v in AXIS_V13[k]: assert any(sp["text"] == f"{v:g}" for sp in SP), (k, v)

# ================================================================== geometry from the base sheet
PITCH = (218.537 - 68.203) / 5.0; GAP = 0.55
def row_y(i): return 68.203 + PITCH * (i + (GAP if i >= 6 else 0.0))
assert abs(row_y(6) - 265.14) < 0.02 and abs(row_y(8) - 325.273) < 0.02
G = dict(ax_l=148.32, ax_r=455.04, tick_len=3.24, bar_h=7.56, cap_h=2.45, name_right=135.35, name_dy=-2.81, unit_dy=8.45, tick_dy=12.56,
         key=[("Oxygenation", C.BLUE, (145.675, 33.987, 150.966, 39.279), (157.59, 39.22)), ("All other measurements", C.GREY, (239.717, 33.987, 245.009, 39.279), (251.63, 39.22))])
def px(k, v):
    lo, _, hi = C_AXIS[k]; return G["ax_l"] + (G["ax_r"] - G["ax_l"]) * (v - lo) / (hi - lo)

page = C.Page(plt, W, H); page.fig.patch.set_alpha(1.0); page.fig.patch.set_facecolor("white")
ax = page.axes(0, 0, W, H); ax.set_xlim(0, W); ax.set_ylim(H, 0); ax.axis("off")
DRAWN = dict(sheet=S, panel="rows", sources={"figure1_sleep_profile_v1.json": dict(path=PROF_PATH, sha256=PROF_SHA, sidecar=PROF_SC), "cohorts.json": dict(path=COH_PATH, sha256=COH_SHA, sidecar=COH_SC),
                                             "t90_final.parquet (v8)": dict(path=T90F, use="quantiles re-derived through cohort_spec.apply_cohort at the stored precision, asserted equal (the sleep rows over the non-masked full nights)")}, values=[])
def rec(text, x, baseline, source, key, value, kind="printed", **extra):
    DRAWN["values"].append(dict(text=text, x=round(x, 2), baseline=round(baseline, 2), source=source, key=key, value=value, kind=kind, **extra))
def T(x, baseline, s, size, ha="left", color=C.INK, static=True):
    page.text(x, baseline, s, size, ha=ha, color=color, record=dict(panel="rows", static=static))

# the key
for lab, col, (x0, y0, x1, y1), (tx, tb) in G["key"]:
    ax.add_patch(plt.Rectangle((x0, y0), x1 - x0, y1 - y0, facecolor=col, edgecolor="none", linewidth=0, zorder=4)); T(tx, tb, lab, 10.0)
rows = []
for i, k in enumerate(C_ROWS):
    m, y = M[k], row_y(i); oxy = m["family"] == "Oxygenation"; col = C.BLUE if oxy else C.GREY; lo, mid, hi = C_AXIS[k]
    ax.plot([G["ax_l"], G["ax_r"]], [y, y], color=C.GREY_PALE, lw=1.0, zorder=2, solid_capstyle="butt")
    for v in (lo, mid, hi):
        ax.plot([px(k, v)] * 2, [y, y + G["tick_len"]], color=C.GREY_PALE, lw=1.0, zorder=2, solid_capstyle="butt")
        T(px(k, v), y + G["tick_dy"], f"{v:g}", 10.0, ha="center", color=C.GREY, static=(k not in AXIS_OLD)); rec(f"{v:g}", px(k, v), y + G["tick_dy"], "axis", f"{k} axis tick", v, kind="tick", ha="center")
    ax.plot([px(k, m["p05"]), px(k, m["p95"])], [y, y], color=col, lw=1.5, zorder=3, solid_capstyle="butt")
    for v in (m["p05"], m["p95"]): ax.plot([px(k, v)] * 2, [y - G["cap_h"], y + G["cap_h"]], color=col, lw=1.5, zorder=3, solid_capstyle="butt")
    ax.add_patch(plt.Rectangle((px(k, m["q1"]), y - G["bar_h"] / 2), max(px(k, m["q3"]) - px(k, m["q1"]), 0.004), G["bar_h"], facecolor=col, edgecolor="none", zorder=4))
    ax.plot([px(k, m["median"])] * 2, [y - G["bar_h"] / 2, y + G["bar_h"] / 2], color="white", lw=1.6, zorder=5, solid_capstyle="butt")
    T(G["name_right"], y + G["name_dy"], m["label"], 10.0, ha="right"); T(G["name_right"], y + G["unit_dy"], C_UNIT[k], 9.0, ha="right", color=C.GREY)
    rec(m["label"], G["name_right"], y + G["name_dy"], "figure1_sleep_profile_v1.json", f"measures.{k}.label", m["label"], kind="label", ha="right")
    geom = dict(row_y=round(y, 3), p05_x=round(px(k, m["p05"]), 3), p95_x=round(px(k, m["p95"]), 3), q1_x=round(px(k, m["q1"]), 3), q3_x=round(px(k, m["q3"]), 3), median_x=round(px(k, m["median"]), 3), colour=col)
    rec(f"{m['label']} drawn quantiles", px(k, m["median"]), y, "figure1_sleep_profile_v1.json", f"measures.{k}: median, q1, q3, p05, p95, n", {q: m[q] for q in ("median", "q1", "q3", "p05", "p95", "n")}, kind="drawn", **geom)
    rows.append(dict(row=i + 1, key=k, label=m["label"], unit=C_UNIT[k], family=m["family"], axis=C_AXIS[k], axis_old=AXIS_OLD.get(k, C_AXIS[k]), values={q: m[q] for q in ("median", "q1", "q3", "p05", "p95", "n")}, rederived=REDERIVED[k], **geom))
DRAWN["geometry"] = dict(G, pitch=PITCH, gap_rows=GAP, rows=rows)
bad = C.banned_words([t["text"] for t in page.texts]); assert not bad, bad
raw_pdf = f"{WORK}/sheet_raw.pdf"; page.fig.savefig(raw_pdf, facecolor="white"); plt.close(page.fig)
met = C.sheet_font_metrics(f"{C.BASE}/{S}.pdf"); done = C.align_font_metrics(raw_pdf, f"{SD}/{S}.pdf", met); print("font metrics aligned:", done)
json.dump(DRAWN, open(f"{VER}/{S}_rows_drawn.json", "w"), indent=1, default=float)
json.dump(page.texts, open(f"{WORK}/placed_texts.json", "w"), indent=0, default=float)
print(f"wrote {SD}/{S}.pdf")
for rw in rows: print(f"  {rw['label']:26s} axis {rw['axis']} median {rw['values']['median']} IQR {rw['values']['q1']}-{rw['values']['q3']} p05-p95 {rw['values']['p05']}-{rw['values']['p95']} n {rw['values']['n']}")
