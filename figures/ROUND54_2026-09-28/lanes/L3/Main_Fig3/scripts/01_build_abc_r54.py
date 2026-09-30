#!/usr/bin/env python3
"""ROUND 54 (2026-09-29, lane L3): Main_Fig3 panels a, b, c re-laid at the round-54 sizes (rows, ticks, keys, cell values, hazard-ratio
and q strings 14 pt, axis titles and column headers 15 pt, small notes and the six-cell grid's cell text 13 pt, letters and title 18 pt
stamped by 02_compose). The 52-row forest is two side-by-side columns of 26 rows (rows 1 to 26 left, 27 to 52 right, each with its own
axis, hazard-ratio column and q column, the key in one row beside the letter a), panels b and c move below the forest and sit side by
side, band d below them. Every string of V30 stays (same strings, larger) plus the strings a second axis needs ("HR (95% CI)", "q",
"Hazard ratio per 1 SD (95% CI)", "0.6", "0.8", "1", "1.5" once more each, declared). Data marks keep their values (asserted against
the round-49 record), their positions follow the new axes (left column over [0.58, 4.15] with the V30 ticks, right column over [0.58,
1.55] with the ticks 0.6, 0.8, 1, 1.5 since its 26 rows span 0.61 to 1.34, the two boxes at one scale: the left box twice the right).
Geometry comes from 00_plan_r54.PLAN (numbers only).
The DATA section (sources, sidecars, checkpoint guard, panel b recompute, panel c) is the round-49 builder's code, unchanged.
Output: work/page_abc.pdf (page W x BAND_TOP, font metrics aligned to the V13 sheet as before) and work/Main_Fig3_drawn.json."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import json, math, os, re, sys
from decimal import Decimal, ROUND_HALF_UP
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.patches import Rectangle
import fitz
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v13_cache
import lane_common as C
import ramp_rule as RR
from plan_loader import PLAN
VL = C.VL

SHEET = "Main_Fig3"; HERE = f"{C.LANE}/{SHEET}"; WORK = f"{HERE}/work"; VER = f"{HERE}/verify"
R49_DRAWN = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/L3/Main_Fig3/work/Main_Fig3_drawn.json"   # the round-49 record (read only): mark values, fills, strings
os.makedirs(WORK, exist_ok=True); os.makedirs(VER, exist_ok=True)
BASE = f"{C.BASE}/{SHEET}.pdf"
CSV = f"{C.NUM}/bdsp_diseases_v3.csv"
MATRIX = f"{C.SV}/dC_side_analyses/hr_matrix_141x52.csv"
MATRIX_CK = f"{C.SV}/dC_side_analyses/_work/hr_matrix_ck/spo2_pct_below_90.json"     # full precision of the same fit (step 168 checkpoint)
RES = f"{C.SV}/dur_x_oxygen_6cell_tst300/results.csv"; SUMM = f"{C.SV}/dur_x_oxygen_6cell_tst300/summary.json"
J1 = f"{paths.FIGURE_ROOT}/ROUND18_2026-09-04/LC_panels/NEW_Fig3_duration_desat_values.json"
PARQ = f"{C.V8}/t90_final.parquet"; COHP = f"{C.NUM}/cohorts.json"
for p in (BASE, CSV, MATRIX, MATRIX_CK, RES, SUMM, J1, PARQ, COHP): C.hydrated(p)
SC111 = C.sidecar(CSV); SC168 = C.sidecar(MATRIX); SC136 = C.sidecar(RES); SC136b = C.sidecar(SUMM); SC156 = C.sidecar(J1); SC101 = C.sidecar(COHP)
for f_ in (C.ARIAL, C.ARIALB): fm.fontManager.addfont(f_)
sys.path.insert(0, paths.ANALYSIS_DIR)
from disease_definitions import DISEASES, NEGATIVE_CONTROLS, SWAPPED_CONTROLS_2026_09_14      # noqa: E402
from cohort_spec import apply_cohort, COHORT_N, FULL_NIGHTS_ONLY                              # noqa: E402
assert not FULL_NIGHTS_ONLY and COHORT_N == 19173, "panel b is recomputed on all rows with the masked sleep amounts (step 156's frame)"

INK, BLUE, GREEN, GRIDC = "#1a1d21", "#0288d1", "#39c445", "#eef0f1"
BLUE_LADDER = ["#b3dcf2", "#7cc0e9", "#3f9fd8", "#0288d1"]
GREEN_LADDER = {"<5": "#d5f1d0", "5–15": "#a9e3a2", "15–35": "#79d475", "≥35": "#39c445"}   # the sheet's count-grid ladder
GREEN3 = ["#d5f1d0", "#79d475", "#39c445"]                                                     # panel c duration headers
BLUE_PRES, BLUE_LOW = "#b3dcf2", "#0288d1"
W, H13 = v13_cache.page("Main_Fig3"); BSP = C.dedupe(v13_cache.spans_lane("Main_Fig3"))   # the V13 text layer from the probe cache (diagnostics only)
assert abs(W - 952.73) < 0.01 and abs(W - PLAN["W"]) < 1e-9, (W, PLAN["W"])
DRAWN = []; RECORDS = []; TXT = []
def rec(panel, text, x, y, ha, size, source, key, value, **kw):
    DRAWN.append(dict(panel=panel, text=text, x=round(x, 3), baseline=round(y, 3), ha=ha, size=size, source=source, key=key, value=value, **kw))
def vrec(panel, text, x, y, ha, size, source_file, source_key, source_value, rule, note=""):
    RECORDS.append(dict(panel=panel, text=text, x=round(x, 3), baseline=round(y, 3), ha=ha, size=size, source_file=source_file, source_key=source_key, source_value=source_value, rule=rule, note=note))


def bh(p):
    p = np.asarray(p, float); m = len(p); o = np.argsort(p); ranks = np.empty(m); ranks[o] = np.arange(1, m + 1)
    q = p * m / ranks; qs = q[o]; qs = np.minimum.accumulate(qs[::-1])[::-1]; out = np.empty(m); out[o] = np.minimum(qs, 1.0); return out


def q_text(q): return VL.q_text(q)          # three decimals half up, a fourth where three would cross the 0.05 boundary
def fmt2_halfup(x): return VL.halfup(x, 2)
def stars(q): return "" if (q is None or not np.isfinite(q)) else ("***" if q < 0.001 else "**" if q < 0.01 else "*" if q < 0.05 else "")
LABEL_SHORT = {"Ventricular arrhythmia or cardiac arrest": "Ventricular arrhythmia/arrest"}   # the lane L2 convention on Fig 2a (owner question)
def label_of(k): return LABEL_SHORT.get(DISEASES[k][0], DISEASES[k][0]) if k in DISEASES else "Death from any cause"

# =============================================================== panel a data (round-49 code, unchanged)
D = pd.read_csv(CSV, comment="#").set_index("key")
M = pd.read_csv(MATRIX)
T = M[M.measurement == "spo2_pct_below_90"].set_index("outcome"); S = M[M.measurement == "TST_min"].set_index("outcome")
CK = {r["outcome"]: r for r in json.load(open(MATRIX_CK))}; assert all(r["measurement"] == "spo2_pct_below_90" for r in CK.values())
N_A = len(T); assert N_A == 52 and len(S) == 52 and sorted(T.index) == sorted(S.index) == sorted(CK), (N_A, len(S), len(CK))
KEYS = list(T.index); assert all(k in D.index for k in KEYS), [k for k in KEYS if k not in D.index]
_CK_GUARD = max(abs(float(T.loc[k, c]) - float(CK[k][c])) for k in KEYS for c in ("HR", "lo", "hi", "p"))
assert _CK_GUARD < 1e-9, ("round 38: the checkpoint is not the matrix's own values", _CK_GUARD)
print(f"checkpoint = matrix on every HR, lo, hi, p of the {len(KEYS)} outcomes (max |diff| {_CK_GUARD:.1e}; step 168 rerun 2026-09-16 17:28 resumed from the checkpoints)")
CTRL_IN_A = sorted(k for k in KEYS if k in NEGATIVE_CONTROLS)
print(f"panel a: negative controls among the {len(KEYS)} ranked outcomes: {CTRL_IN_A}")
for k in KEYS:
    assert all(round(T.loc[k, a], 3) == round(D.loc[k, b], 3) for a, b in (("HR", "t90_hr"), ("lo", "t90_lo"), ("hi", "t90_hi"))), (k, T.loc[k, ["HR", "lo", "hi"]].tolist(), D.loc[k, ["t90_hr", "t90_lo", "t90_hi"]].tolist())
    assert all(round(CK[k][a], 3) == round(D.loc[k, b], 3) for a, b in (("HR", "t90_hr"), ("lo", "t90_lo"), ("hi", "t90_hi"))), (k, "checkpoint drifted from the csv at 3 decimals")
    assert abs(CK[k]["p"] - float(T.loc[k, "p"])) <= 1e-12 * max(1.0, abs(CK[k]["p"])) or abs(CK[k]["p"] - float(T.loc[k, "p"])) < 1e-300, (k, "checkpoint p differs from the matrix")
    assert all(abs(S.loc[k, a] - D.loc[k, b]) <= 0.005 for a, b in (("HR", "tst_hr"), ("lo", "tst_lo"), ("hi", "tst_hi"))), (k, "TST lane drifted between the csv and the matrix")
TST_DRIFT = max(abs(float(S.loc[k, a]) - float(D.loc[k, b])) for k in KEYS for a, b in (("HR", "tst_hr"), ("lo", "tst_lo"), ("hi", "tst_hi")))
print(f"TST lane (markers only, drawn from the csv): largest csv-vs-matrix difference {TST_DRIFT:.4f}")
Dk = D.loc[KEYS]
q_t90_csv = pd.Series(bh(Dk.t90_p.values), index=KEYS); q_tst_csv = pd.Series(bh(Dk.tst_p.values), index=KEYS)
t90_sig = set(q_t90_csv.index[q_t90_csv < 0.05]); tst_sig = set(q_tst_csv.index[q_tst_csv < 0.05])
assert t90_sig == set(T.index[T.q < 0.05]) and tst_sig == set(S.index[S.q < 0.05]), "BH within the exposure across the 52: csv-derived and matrix q disagree on the significant set"
ORDER = sorted(KEYS, key=lambda k: (-CK[k]["HR"], label_of(k)))
LABELS = [label_of(k) for k in ORDER]; assert len(set(LABELS)) == N_A
base_labels = sorted({s["text"] for s in BSP if abs(s["origin"][0] - 8.88) < 0.05 and 100 < s["origin"][1] < 730 and s["size"] > 10.5})
ROWS_IN = sorted(set(LABELS) - set(base_labels)); ROWS_OUT = sorted(set(base_labels) - set(LABELS))
print(f"panel a rows: {N_A} (V13 had {len(base_labels)}); enter {ROWS_IN}; leave {ROWS_OUT}")
ROWS_A = []
for i, k in enumerate(ORDER):
    hr, lo, hi = (float(CK[k][c]) for c in ("HR", "lo", "hi")); q = float(T.loc[k, "q"])
    s_hr = VL.hr_ci(hr, lo, hi); s_q = q_text(q); s_q_csv = q_text(float(q_t90_csv[k]))
    assert s_q == s_q_csv, (k, "printed q differs between the matrix q and the csv-derived BH q", s_q, s_q_csv)
    ROWS_A.append(dict(key=k, label=label_of(k), hr=hr, lo=lo, hi=hi, q=q, q_csv=float(q_t90_csv[k]), sig=k in t90_sig, s_hr=s_hr, s_q=s_q,
                       tst_hr=float(Dk.loc[k, "tst_hr"]), tst_lo=float(Dk.loc[k, "tst_lo"]), tst_hi=float(Dk.loc[k, "tst_hi"]), tst_p=float(Dk.loc[k, "tst_p"]), tst_q=float(q_tst_csv[k]), tst_sig=k in tst_sig,
                       csv_t90=[float(Dk.loc[k, c]) for c in ("t90_hr", "t90_lo", "t90_hi")]))
assert sum(r["sig"] for r in ROWS_A) == len(t90_sig) and sum(r["tst_sig"] for r in ROWS_A) == len(tst_sig)
print(f"panel a: {N_A} rows, T90 BH-significant {len(t90_sig)} of {N_A}, total sleep time {len(tst_sig)} of {N_A} {sorted(label_of(k) for k in tst_sig)}; first rows {LABELS[:6]}")
tie_cells = [(k, c, float(Dk.loc[k, c])) for k in KEYS for c in ("t90_hr", "t90_lo", "t90_hi") if abs(round(Dk.loc[k, c] * 1000) % 10) == 5]
tie_differ = [(k, c, v, fmt2_halfup(CK[k][{"t90_hr": "HR", "t90_lo": "lo", "t90_hi": "hi"}[c]]), fmt2_halfup(v)) for k, c, v in tie_cells if fmt2_halfup(CK[k][{"t90_hr": "HR", "t90_lo": "lo", "t90_hi": "hi"}[c]]) != fmt2_halfup(v)]
print(f"panel a rounding: {len(tie_cells)} csv values sit at a 3-decimal tie, {len(tie_differ)} print differently half-up from the csv than from the full-precision checkpoint (printed from the checkpoint): {tie_differ}")

# =============================================================== panel b data (round-49 code, unchanged)
J = json.load(open(J1)); DUR = ["<5 h", "6–7 h", ">7 h"]; BANDS = ["≤1%", ">1–5%", ">5–10%", ">10%"]
b = apply_cohort(pd.read_parquet(PARQ)); assert len(b) == COHORT_N == 19173 == json.load(open(COHP))["bdsp"]["n_analysis"], len(b)
assert b.spo2_pct_below_90.notna().all()
ox = b.spo2_pct_below_90.map(lambda v: 0 if v <= 1 else 1 if v <= 5 else 2 if v <= 10 else 3); tst = b.TST_min
dur = pd.Series(np.where(tst < 300, DUR[0], np.where((tst >= 360) & (tst <= 420), DUR[1], np.where(tst > 420, DUR[2], "outside"))), index=b.index)
N_TST_MISSING = int(tst.isna().sum())
def pct1(n, d): return str((Decimal(n) * 100 / Decimal(d)).quantize(Decimal("0.1"), rounding=ROUND_HALF_UP))
def band_of(p): p = float(p); return "<5" if p < 5 else "5–15" if p < 15 else "15–35" if p < 35 else "≥35"
GRID_B = {}
for r_ in DUR:
    m = dur == r_; n = int(m.sum()); cells = [int(((ox == c) & m).sum()) for c in range(4)]; assert sum(cells) == n
    jr = J["rows"][r_]; assert jr["total"] == n and [c["n"] for c in jr["cells"]] == cells, (r_, n, cells, jr)
    pct = [pct1(c, n) for c in cells]; assert pct == [c["pct"] for c in jr["cells"]], (r_, pct, jr)
    GRID_B[r_] = dict(total=n, cells=[dict(band=BANDS[c], n=cells[c], pct=pct[c], label=f"{cells[c]:,} ({pct[c]})", shade=band_of(pct[c])) for c in range(4)])
    assert [c["label"] for c in GRID_B[r_]["cells"]] == [c["label"] for c in jr["cells"]]
min_pct = min(float(c["pct"]) for r_ in DUR for c in GRID_B[r_]["cells"])
UNDER5 = [(r_, c["band"], c["pct"]) for r_ in DUR for c in GRID_B[r_]["cells"] if float(c["pct"]) < 5.0]
print("panel b: direct recompute from the v8 table equals the step-156 file on all 12 cells and 3 totals;", {r_: GRID_B[r_]["total"] for r_ in DUR}, "smallest cell", min_pct, "under 5:", UNDER5, "TST missing (masked split nights)", N_TST_MISSING)

# =============================================================== panel c data (round-49 code, unchanged)
R = pd.read_csv(RES); R = R[R.analysis == "primary"].set_index("key"); SUM = json.load(open(SUMM))
CTRL_KEY = SWAPPED_CONTROLS_2026_09_14["back_pain"]; assert CTRL_KEY == "alopecia" and CTRL_KEY in NEGATIVE_CONTROLS
CTRL_LABEL = f"{DISEASES[CTRL_KEY][0]} (negative control)"
KEYS3 = {"COPD": "copd2", "Heart failure": "hf", "Type 2 diabetes": "diabetes", "Cardiovascular composite": "cvd", "Death, any cause": "death", "Dyslipidemia": "dyslipid", CTRL_LABEL: CTRL_KEY}
CELLS = ["SN", "SL", "NL", "LN", "LL"]; CELL_POS = {"SN": (0, 0), "SL": (0, 1), "NN": (1, 0), "NL": (1, 1), "LN": (2, 0), "LL": (2, 1)}
cell_n = SUM["cell_n"]; assert cell_n["SN"] == GRID_B["<5 h"]["cells"][0]["n"] and cell_n["SL"] == GRID_B["<5 h"]["cells"][3]["n"] and cell_n["NN"] == GRID_B["6–7 h"]["cells"][0]["n"] \
    and cell_n["NL"] == GRID_B["6–7 h"]["cells"][3]["n"] and cell_n["LN"] == GRID_B[">7 h"]["cells"][0]["n"] and cell_n["LL"] == GRID_B[">7 h"]["cells"][3]["n"], (cell_n, "six-cell sizes differ from the count grid")
assert SUM["discards"]["cohort"] == 15551
C_ROWS = {}
for title, k in KEYS3.items():
    r = R.loc[k]; assert (k == CTRL_KEY) == bool(r.negative_control), (title, r.negative_control)
    cells = {}
    for c in CELLS:
        hr, lo, hi, p = (float(r[f"{c}_{s}"]) for s in ("hr", "lo", "hi", "p")); q = float(r[f"{c}_q"]) if f"{c}_q" in r and np.isfinite(r[f"{c}_q"]) else float("nan")
        assert lo < hr < hi, (title, c, hr, lo, hi)
        cells[c] = dict(hr=hr, lo=lo, hi=hi, p=p, q=q, s_hr=fmt2_halfup(hr) + stars(q), s_ci=f"({fmt2_halfup(lo)}-{fmt2_halfup(hi)})", fill=RR.ramp_hex(hr), ink=RR.text_on(RR.ramp_hex(hr)))
    C_ROWS[title] = dict(key=k, cells=cells, row_max=max(v["hr"] for v in cells.values()), control=bool(r.negative_control), disease=str(r.disease))
C_ORDER = sorted([t for t in KEYS3 if not C_ROWS[t]["control"]], key=lambda t: -C_ROWS[t]["row_max"]) + [t for t in KEYS3 if C_ROWS[t]["control"]]
assert len(C_ORDER) == 7 and C_ORDER[-1] == CTRL_LABEL
base_c_order = [s["text"] for s in sorted([s for s in BSP if 518 < s["origin"][0] < 600 and 440 < s["origin"][1] < 700 and s["size"] > 9.4 and s["text"] in ("COPD", "Heart failure", "Type 2 diabetes", "Cardiovascular", "Death, any cause", "Dyslipidemia", "Back pain", "Alopecia")], key=lambda s: s["origin"][1])]
print("panel c: row order by the largest hazard ratio in the row (control last): new", C_ORDER, "| V13", base_c_order)
c_ties = [(t, c, C_ROWS[t]["cells"][c][s]) for t in KEYS3 for c in CELLS for s in ("hr", "lo", "hi") if abs(round(C_ROWS[t]["cells"][c][s] * 1000) % 10) == 5]
print(f"panel c rounding: {len(c_ties)} of 105 csv values sit at a 3-decimal tie (printed half-up on the stored value, the file is 3 dp): {c_ties}")
sn_ln_sig = [(t, c, C_ROWS[t]["cells"][c]["q"]) for t in KEYS3 for c in ("SN", "LN") if np.isfinite(C_ROWS[t]["cells"][c]["q"]) and C_ROWS[t]["cells"][c]["q"] < 0.05]
ctrl_excl = [(c, C_ROWS[CTRL_LABEL]["cells"][c]["s_hr"], C_ROWS[CTRL_LABEL]["cells"][c]["s_ci"]) for c in CELLS if not (C_ROWS[CTRL_LABEL]["cells"][c]["lo"] <= 1 <= C_ROWS[CTRL_LABEL]["cells"][c]["hi"])]
print(f"panel c: preserved-oxygen cells (SN, LN) clearing BH on the drawn rows: {sn_ln_sig or 'none'}; control row {CTRL_LABEL} cells whose CI excludes 1: {ctrl_excl or 'none'}")

# =============================================================== round 54: the values against the round-49 record (marks, fills, strings)
DJ49 = json.load(open(C.hydrated(R49_DRAWN)))
_r49 = {r["key"]: r for r in DJ49["panel_a"]["rows"]}
assert [r["key"] for r in ROWS_A] == [r["key"] for r in DJ49["panel_a"]["rows"]], "round 54: the row order differs from the round-49 record"
for r in ROWS_A:
    o = _r49[r["key"]]
    assert all(abs(r[c] - o[c]) < 1e-12 for c in ("hr", "lo", "hi", "q", "tst_hr", "tst_lo", "tst_hi")) and r["sig"] == o["sig"] and r["tst_sig"] == o["tst_sig"] and r["s_hr"] == o["s_hr"] and r["s_q"] == o["s_q"] and r["label"] == o["label"], (r["key"], "round 54: a panel a value differs from the round-49 record")
for t in C_ORDER:
    o = DJ49["panel_c"]["rows"][t]
    for c in CELLS:
        assert all(abs(C_ROWS[t]["cells"][c][s] - o["cells"][c][s]) < 1e-12 for s in ("hr", "lo", "hi")) and all(C_ROWS[t]["cells"][c][s] == o["cells"][c][s] for s in ("s_hr", "s_ci", "fill", "ink")), (t, c, "round 54: a panel c cell differs from the round-49 record")
assert C_ORDER == DJ49["panel_c"]["order"] and {r_: GRID_B[r_] for r_ in DUR} == {r_: DJ49["panel_b"]["rows"][r_] for r_ in DUR}, "round 54: panel b or the panel c order differs from the round-49 record"
print("round 54: every panel a value (52 x HR, lo, hi, q, sleep-time HR, lo, hi, both flags, both strings), every panel c cell (35 x hr, lo, hi, strings, fill, ink) and every panel b cell equals the round-49 record")

# =============================================================== geometry (page points, y down) from the planner
SZ = PLAN["sizes"]; Y = PLAN["y"]; PB = PLAN["b"]; PC = PLAN["c"]
S_A = 1.0706                                                   # the sheet's scale of the figure2C_polish design (marker sizes, line widths)
S_CIRC, S_SQ = (math.sqrt(28.0) * S_A) ** 2, (math.sqrt(23.0) * S_A) ** 2
LW_CI, LW_MEDGE, LW_OPEN = 1.7 * S_A, 0.8 * S_A, 1.1 * S_A
N_COL = PLAN["n_col"]; PITCH = Y["pitch"]; ROW0 = Y["row0"]; LANE = Y["lane"]; AX_TOP, AX_BOT = Y["ax_top"], Y["ax_bot"]
TICK_BASE, CAP_BASE, HDR_BASE = Y["tick_base"], Y["cap_base"], Y["hdr_base"]
LAB_DY = round(0.36 * SZ["row"], 3)                             # label baseline below the row centre (x-height centre on the row centre)
HRQ_DY = round(-LANE + 0.36 * SZ["hrq"] - 1.3, 3)               # HR and q strings on the T90 lane, as in V30 (x-height centre 1.3 pt above the blue whisker)
COLS = []
for name, rows in (("left", ROWS_A[:N_COL]), ("right", ROWS_A[N_COL:])):
    p = PLAN[name]; XL, XR, xlo, xhi = p["box0"], p["box1"], p["xlo"], p["xhi"]
    XB = (XR - XL) / math.log(xhi / xlo); XA = XL - XB * math.log(xlo)
    vals = [v for r in rows for v in (r["lo"], r["hi"], r["tst_lo"], r["tst_hi"])]
    assert xlo < min(vals) and max(vals) < xhi, (name, min(vals), max(vals), xlo, xhi)
    lab_end = max(p["lab_x0"] + fitz.Font(fontfile=C.ARIAL).text_length(r["label"], fontsize=SZ["row"]) for r in rows)
    assert lab_end + PLAN["clear_label"] <= XL + 1e-6, (name, "a row label would cross the 3 pt clearance to the box", lab_end, XL)
    COLS.append(dict(name=name, rows=rows, XA=XA, XB=XB, x_left=XL, x_right=XR, xlo=xlo, xhi=xhi, ticks=p["ticks"], lab_x0=p["lab_x0"], hr_r=p["hr_r"], q_r=p["q_r"], lab_end=lab_end))
    print(f"panel a {name} column: rows {rows[0]['label']} .. {rows[-1]['label']}, box {XL} to {XR:.2f}, x = {XA:.3f} + {XB:.3f} ln HR over [{xlo}, {xhi}], labels end at {lab_end:.1f}, HR right edge {p['hr_r']:.2f}, q right edge {p['q_r']:.2f}")
BAND_TOP = Y["band_top"]; H_PAGE = BAND_TOP                     # the page piece ends where the band starts
KEY_ROWS = [("blue", "T90, percent of the recording below 90% saturation"), ("green", "Total sleep time on the study night"), ("open", "Not significant")]

RC = {"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.bbox": "standard", "savefig.pad_inches": 0.0, "figure.dpi": 100,
      "text.color": INK, "axes.unicode_minus": False, "lines.scale_dashes": False, "lines.solid_capstyle": "butt", "lines.dash_capstyle": "butt"}
FONT_R = fitz.Font(fontfile=C.ARIAL); FONT_B = fitz.Font(fontfile=C.ARIALB); SHIFT = {}; RIGHT = []
def wtxt(s, size, bold=False): return (FONT_B if bold else FONT_R).text_length(s, fontsize=size)
for PASS in (1, 2):
  DRAWN.clear(); RECORDS.clear(); RIGHT.clear(); TXT.clear()
  with plt.rc_context(RC):
      fig = plt.figure(figsize=(W / 72.0, H_PAGE / 72.0)); ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(H_PAGE, 0); ax.axis("off")
      def txt(x, y, s, size, ha="left", bold=False, color=INK, z=5):
          """x is the anchor the design gives the string: the left origin (ha left), the centre (ha center) or the right edge (ha right).
          Right-aligned strings are drawn left-aligned at x minus their advance, corrected in pass 2 by the measured edge (SHIFT, keyed by
          text, size and weight, so a string drawn twice at one baseline gets one correction)."""
          f = FONT_B if bold else FONT_R
          TXT.append(dict(text=s, x=round(x, 3), baseline=round(y, 3), ha=ha, size=size, bold=bold))
          if ha == "right":
              key = (s, size, bold); adv = f.text_length(s, fontsize=size)
              RIGHT.append((key, x, y)); ax.text(x - adv - SHIFT.get(key, 0.0), y, s, fontsize=size, ha="left", va="baseline", color=color, fontweight="bold" if bold else "normal", zorder=z); return
          ax.text(x, y, s, fontsize=size, ha=ha, va="baseline", color=color, fontweight="bold" if bold else "normal", zorder=z)
      def rect(x, y, w, h, fc, ec="none", lw=0.0, z=1):
          ax.add_patch(Rectangle((x, y), w, h, facecolor=fc, edgecolor=ec, linewidth=lw, zorder=z))
      # ---------------- panel a: the key in one row beside the letter a
      kx = Y["key_x0"]; ky = Y["key_y"]; kb = Y["key_base"]; KEY_POS_A = []
      for kind, s in KEY_ROWS:
          if kind == "blue":
              ax.plot([kx, kx + 15.94], [ky, ky], color=BLUE, lw=LW_CI, solid_capstyle="round", zorder=2)
              ax.scatter([kx + 7.97], [ky], s=(6.0 * S_A) ** 2, marker="o", facecolors=BLUE, edgecolors="white", linewidths=LW_MEDGE, zorder=3); xl = kx + 21.3
          elif kind == "green":
              ax.plot([kx, kx + 15.94], [ky, ky], color=GREEN, lw=LW_CI, solid_capstyle="round", zorder=2)
              ax.scatter([kx + 7.97], [ky], s=(5.4 * S_A) ** 2, marker="s", facecolors=GREEN, edgecolors="white", linewidths=LW_MEDGE, zorder=3); xl = kx + 21.3
          else:
              ax.scatter([kx + 8.03], [ky], s=(6.0 * S_A) ** 2, marker="o", facecolors="white", edgecolors=GREEN, linewidths=LW_OPEN, zorder=3); xl = kx + 21.33
          txt(xl, kb, s, SZ["key"]); KEY_POS_A.append(dict(entry=kind, mark_x=round(kx, 3), label_x=round(xl, 3), label_end=round(xl + wtxt(s, SZ["key"]), 3)))
          kx = xl + wtxt(s, SZ["key"]) + 20.0
      assert kx - 20.0 <= W - PLAN["right_margin"], ("the panel a key runs off the sheet", kx)
      # ---------------- panel a: the two columns
      for ci, col in enumerate(COLS):
          XA, XB, XL, XR, rows = col["XA"], col["XB"], col["x_left"], col["x_right"], col["rows"]
          def xhr(v, XA=XA, XB=XB): return XA + XB * math.log(v)
          for v in col["ticks"]:
              if v != 1: ax.plot([xhr(v)] * 2, [AX_TOP, AX_BOT], color=GRIDC, lw=0.5 * S_A, zorder=0.3)
          ax.plot([xhr(1.0)] * 2, [AX_BOT, AX_TOP], color=INK, lw=0.9 * S_A, ls=(0, (4 * 0.9 * S_A, 3 * 0.9 * S_A)), zorder=1)
          ax.plot([xhr(col["xlo"]), xhr(col["xhi"])], [AX_BOT] * 2, color=INK, lw=0.8 * S_A, solid_capstyle="projecting", zorder=4)
          for v in col["ticks"]:
              ax.plot([xhr(v)] * 2, [AX_BOT, AX_BOT + 3.0 * S_A], color=INK, lw=0.8 * S_A, zorder=4)
              txt(xhr(v), TICK_BASE, f"{v:g}", SZ["tick"], ha="center"); rec("a", f"{v:g}", xhr(v), TICK_BASE, "center", SZ["tick"], "axis", f"tick ({col['name']} column)", v)
          txt((XL + XR) / 2, CAP_BASE, "Hazard ratio per 1 SD (95% CI)", SZ["cap"], ha="center")
          txt(col["hr_r"], HDR_BASE, "HR (95% CI)", SZ["hdr"], ha="right"); txt(col["q_r"], HDR_BASE, "q", SZ["hdr"], ha="right")
          for i, r in enumerate(rows):
              yc = ROW0 + PITCH * i; yb, yg = yc - LANE, yc + LANE; irow = ci * N_COL + i
              ax.plot([xhr(r["lo"]), xhr(r["hi"])], [yb, yb], color=BLUE, lw=LW_CI, solid_capstyle="round", zorder=2)
              ax.plot([xhr(r["tst_lo"]), xhr(r["tst_hi"])], [yg, yg], color=GREEN, lw=LW_CI, solid_capstyle="round", zorder=2)
              if r["sig"]: ax.scatter([xhr(r["hr"])], [yb], s=S_CIRC, marker="o", facecolors=BLUE, edgecolors="white", linewidths=LW_MEDGE, zorder=3)
              else: ax.scatter([xhr(r["hr"])], [yb], s=S_CIRC, marker="o", facecolors="white", edgecolors=BLUE, linewidths=LW_OPEN, zorder=3)
              if r["tst_sig"]: ax.scatter([xhr(r["tst_hr"])], [yg], s=S_SQ, marker="s", facecolors=GREEN, edgecolors="white", linewidths=LW_MEDGE, zorder=3)
              else: ax.scatter([xhr(r["tst_hr"])], [yg], s=S_SQ, marker="s", facecolors="white", edgecolors=GREEN, linewidths=LW_OPEN, zorder=3)
              y_lab, y_hq = yc + LAB_DY, yc + HRQ_DY
              txt(col["lab_x0"], y_lab, r["label"], SZ["row"]); txt(col["hr_r"], y_hq, r["s_hr"], SZ["hrq"], ha="right"); txt(col["q_r"], y_hq, r["s_q"], SZ["hrq"], ha="right", bold=r["sig"])
              rec("a", r["label"], col["lab_x0"], y_lab, "left", SZ["row"], MATRIX_CK, f"[outcome={r['key']}] (row {irow + 1}, sorted by the T90 hazard ratio, {col['name']} column)", None)
              rec("a", r["s_hr"], col["hr_r"], y_hq, "right", SZ["hrq"], MATRIX_CK, f"[outcome={r['key']}]: HR, lo, hi (full precision; csv t90_hr/lo/hi equal at 3 decimals)", [r["hr"], r["lo"], r["hi"]], csv=r["csv_t90"], sidecar=SC168["sha256"])
              rec("a", r["s_q"], col["q_r"], y_hq, "right", SZ["hrq"], MATRIX, f"measurement=spo2_pct_below_90&outcome={r['key']}:q (BH across the 52 within the exposure; csv-derived q {r['q_csv']:.6g})", r["q"], bold=r["sig"], sidecar=SC168["sha256"])
              vrec("a", r["label"], col["lab_x0"], y_lab, "left", SZ["row"], CSV, f"key={r['key']}:disease", r["label"], "label_alias", f"row {irow + 1}, {col['name']} column")
              vrec("a", r["s_hr"], col["hr_r"], y_hq, "right", SZ["hrq"], MATRIX_CK, f"[outcome={r['key']}]", {"HR": r["hr"], "lo": r["lo"], "hi": r["hi"]}, "hr_ci_HRlohi", "half up on the step-168 full-precision checkpoint; bdsp_diseases_v3.csv equal at 3 dp")
              vrec("a", r["s_q"], col["q_r"], y_hq, "right", SZ["hrq"], MATRIX, f"measurement=spo2_pct_below_90&outcome={r['key']}:q", r["q"], "q3", "BH across the 52 ranked outcomes within the exposure (= the csv-derived q)")
              DRAWN.append(dict(panel="a", marker="T90 circle", filled=r["sig"], x=round(xhr(r["hr"]), 3), y=round(yb, 3), ci_x=[round(xhr(r["lo"]), 3), round(xhr(r["hi"]), 3)], key=r["key"], column=col["name"], row=irow + 1, value=[r["hr"], r["lo"], r["hi"]]))
              DRAWN.append(dict(panel="a", marker="total sleep time square", filled=r["tst_sig"], x=round(xhr(r["tst_hr"]), 3), y=round(yg, 3), ci_x=[round(xhr(r["tst_lo"]), 3), round(xhr(r["tst_hi"]), 3)], key=r["key"], column=col["name"], row=irow + 1,
                                source=CSV, value=[r["tst_hr"], r["tst_lo"], r["tst_hi"]], q=r["tst_q"], sidecar=SC111["sha256"]))
      # ---------------- panel b (below the forest, left)
      BX0 = PB["x0"]; BCELL = PB["cell"]; HX = [BX0 + PB["gut"] + (BCELL + PB["gap_h"]) * c for c in range(4)]; B_CENTRE = (HX[0] + HX[3] + BCELL) / 2
      txt(B_CENTRE, Y["b_title_base"], "T90, % of the recording", SZ["b_title"], ha="center")
      for c in range(4):
          rect(HX[c], Y["b_hdr_y0"], BCELL, Y["b_hdr_h"], BLUE_LADDER[c]); txt(HX[c] + BCELL / 2, Y["b_hdr_base"], BANDS[c], SZ["b_hdr"], ha="center")
      B_LAB_R = BX0 + PB["lab_w"] + 1.0
      for r_i, r_ in enumerate(DUR):
          y0 = Y["b_row_y0"] + Y["b_row_pitch"] * r_i; yb_ = round(y0 + Y["b_row_h"] / 2 + 0.36 * SZ["b_cell"], 3)
          txt(B_LAB_R, yb_, r_, SZ["b_row"], ha="right")
          for c in range(4):
              cell = GRID_B[r_]["cells"][c]; fc = GREEN_LADDER[cell["shade"]]
              rect(HX[c], y0, BCELL, Y["b_row_h"], fc); txt(HX[c] + BCELL / 2, yb_, cell["label"], SZ["b_cell"], ha="center")
              rec("b", cell["label"], HX[c] + BCELL / 2, yb_, "center", SZ["b_cell"], PARQ, f"{r_} x {cell['band']}: n and percent of the duration group (direct recompute = step-156 file)", [cell["n"], cell["pct"]], fill=fc, band=cell["shade"], sidecar=SC156["sha256"])
              vrec("b", cell["label"], HX[c] + BCELL / 2, yb_, "center", SZ["b_cell"], J1, f"rows/{r_}/cells/{c}/label", cell["label"], "text", f"n {cell['n']} of {GRID_B[r_]['total']}, recomputed from the v8 table")
      txt(BX0, Y["b_key_base"], "Percent of the group", SZ["b_key"])
      KEY_BANDS = [b_ for b_ in J["shade_bands_pct"] if any(c["shade"] == b_ for r_ in DUR for c in GRID_B[r_]["cells"])]
      assert KEY_BANDS and all(b_ in GREEN_LADDER for b_ in KEY_BANDS), KEY_BANDS
      _kx = BX0 + wtxt("Percent of the group", SZ["b_key"]) + PB["key_lead_gap"]; KEY_POS = []
      for lab in KEY_BANDS:
          xl = _kx + 17.42; rect(_kx, Y["b_key_sw_y"], 12.42, 10.76, GREEN_LADDER[lab], ec="white", lw=0.6, z=2); txt(xl, Y["b_key_base"], lab, SZ["b_key"])
          KEY_POS.append(dict(band=lab, swatch_x=round(_kx, 3), label_x=round(xl, 3)))
          vrec("b", lab, xl, Y["b_key_base"], "left", SZ["b_key"], J1, f"shade_bands_pct/{J['shade_bands_pct'].index(lab)}", lab, "text", "key rung, drawn because a cell of the grid falls in this band")
          _kx = xl + wtxt(lab, SZ["b_key"]) + 20.0
      assert _kx - 20.0 <= PB["right"], ("the panel b key runs past the panel b block", _kx - 20.0, PB["right"])
      # ---------------- panel c (below the forest, right of b)
      CX0 = PC["x0"]; CC = PC["cell"]; PAIR = 2 * CC + PC["pair_gap"]; GX = [PC["cell_x0"] + (PAIR + PC["group_gap"]) * g for g in range(3)]; C_LAB_R = CX0 + PC["lab_w"]
      for g in range(3): rect(GX[g], Y["c_hdr1_y0"], PAIR, Y["c_hdr1_h"], GREEN3[g])
      y_h1 = round(Y["c_hdr1_y0"] + Y["c_hdr1_h"] / 2 + 0.36 * SZ["c_hdr"], 3)
      for g, s in enumerate(DUR): txt(GX[g] + PAIR / 2, y_h1, s, SZ["c_hdr"], ha="center")
      y_h2a, y_h2b = round(Y["c_hdr2_y0"] + 13.8, 3), round(Y["c_hdr2_y0"] + 27.3, 3)
      for g in range(3):
          rect(GX[g], Y["c_hdr2_y0"], CC, Y["c_hdr2_h"], BLUE_PRES); rect(GX[g] + CC + PC["pair_gap"], Y["c_hdr2_y0"], CC, Y["c_hdr2_h"], BLUE_LOW)
          txt(GX[g] + CC / 2, y_h2a, "Preserved", SZ["c_sub"], ha="center"); txt(GX[g] + CC + PC["pair_gap"] + CC / 2, y_h2a, "Low", SZ["c_sub"], ha="center")
          txt(GX[g] + CC / 2, y_h2b, "T90 ≤1%", SZ["c_sub"], ha="center"); txt(GX[g] + CC + PC["pair_gap"] + CC / 2, y_h2b, "T90 >10%", SZ["c_sub"], ha="center")
      txt(C_LAB_R, y_h1, "Sleep duration", SZ["c_rowhdr"], ha="right"); txt(C_LAB_R, round(Y["c_hdr2_y0"] + Y["c_hdr2_h"] / 2 + 0.36 * SZ["c_rowhdr"], 3), "Oxygenation", SZ["c_rowhdr"], ha="right")
      rect(GX[1], Y["c_row_y0"], CC, 7 * Y["c_row_pitch"] - (Y["c_row_pitch"] - Y["c_row_h"]), "white", ec=INK, lw=0.9, z=2)
      TWO = {"Cardiovascular composite": ("Cardiovascular", "composite"), CTRL_LABEL: (DISEASES[CTRL_KEY][0], "(negative control)")}
      L2_DY = 14.5; L2_A, L2_B = round(-(L2_DY + 0.716 * SZ["c_row"] + 0.212 * SZ["c_row"]) / 2 + 0.716 * SZ["c_row"], 3), None; L2_B = round(L2_A + L2_DY, 3)   # two 14-pt lines centred on the row
      C2_A, C2_B = round(-(Y["c_line_dy"] + 0.716 * SZ["c_cell"] + 0.212 * SZ["c_cell"]) / 2 + 0.716 * SZ["c_cell"], 3), None; C2_B = round(C2_A + Y["c_line_dy"], 3)   # two 13-pt lines centred on the row
      ONE_DY = round(0.36 * SZ["c_row"], 3); ONE_DY13 = round(0.36 * SZ["c_cell"], 3)
      for j, title in enumerate(C_ORDER):
          y0 = Y["c_row_y0"] + Y["c_row_pitch"] * j; yc = y0 + Y["c_row_h"] / 2; kk = C_ROWS[title]["key"]
          if title in TWO:
              txt(C_LAB_R, yc + L2_A, TWO[title][0], SZ["c_row"], ha="right"); txt(C_LAB_R, yc + L2_B, TWO[title][1], SZ["c_row"], ha="right")
              vrec("c", TWO[title][0], C_LAB_R, yc + L2_A, "right", SZ["c_row"], RES, f"analysis=primary&key={kk}:disease", C_ROWS[title]["disease"], "label_alias", f"row {j + 1}, first line")
              vrec("c", TWO[title][1], C_LAB_R, yc + L2_B, "right", SZ["c_row"], RES, f"analysis=primary&key={kk}:negative_control", bool(C_ROWS[title]["control"]), "label_alias", f"row {j + 1}, second line")
              rec("c", title, C_LAB_R, yc + L2_A, "right", SZ["c_row"], RES, f"key={kk} (row {j + 1}, rows by the largest hazard ratio, control last)", None)
          else:
              txt(C_LAB_R, yc + ONE_DY, title, SZ["c_row"], ha="right")
              vrec("c", title, C_LAB_R, yc + ONE_DY, "right", SZ["c_row"], RES, f"analysis=primary&key={kk}:disease", C_ROWS[title]["disease"], "label_alias", f"row {j + 1}")
              rec("c", title, C_LAB_R, yc + ONE_DY, "right", SZ["c_row"], RES, f"key={kk} (row {j + 1}, rows by the largest hazard ratio, control last)", None)
          txt(GX[1] + CC / 2, yc + ONE_DY13, "1.00", SZ["c_cell"], ha="center"); rec("c", "1.00", GX[1] + CC / 2, yc + ONE_DY13, "center", SZ["c_cell"], RES, "reference cell NN (6–7 h, T90 ≤1%)", 1.0)
          vrec("c", "1.00", GX[1] + CC / 2, yc + ONE_DY13, "center", SZ["c_cell"], SUMM, "thresholds/reference_cell", "NN", "ref_one", "the reference cell prints 1.00 by definition")
          for c in CELLS:
              g, l = CELL_POS[c]; x0 = GX[g] + (CC + PC["pair_gap"]) * l; cell = C_ROWS[title]["cells"][c]; xc = x0 + CC / 2
              rect(x0, y0, CC, Y["c_row_h"], cell["fill"], z=1)
              txt(xc, yc + C2_A, cell["s_hr"], SZ["c_cell"], ha="center", color=cell["ink"]); txt(xc, yc + C2_B, cell["s_ci"], SZ["c_cell"], ha="center", color=cell["ink"])
              rec("c", cell["s_hr"], xc, yc + C2_A, "center", SZ["c_cell"], RES, f"analysis=primary key={kk} {c}_hr with stars from {c}_q", [cell["hr"], cell["q"]], fill=cell["fill"], ink=cell["ink"], sidecar=SC136["sha256"])
              rec("c", cell["s_ci"], xc, yc + C2_B, "center", SZ["c_cell"], RES, f"analysis=primary key={kk} {c}_lo, {c}_hi", [cell["lo"], cell["hi"]], fill=cell["fill"], ink=cell["ink"], sidecar=SC136["sha256"])
              vrec("c", cell["s_hr"], xc, yc + C2_A, "center", SZ["c_cell"], RES, f"analysis=primary&key={kk}:{c}_hr,{c}_q", [cell["hr"], cell["q"]], "hr_stars", "half up on the stored 3 dp value; stars by the cell's BH q (none on the control row)")
              vrec("c", cell["s_ci"], xc, yc + C2_B, "center", SZ["c_cell"], RES, f"analysis=primary&key={kk}:{c}_lo,{c}_hi", [cell["lo"], cell["hi"]], "ci_paren", "half up on the stored 3 dp value")
      # panel c legend: the ramp bar with its title and ticks, the reference box, the asterisk key
      bar0, bar1 = PC["bar0"], PC["bar1"]; BAR_W = bar1 - bar0
      txt((bar0 + bar1) / 2, Y["c_leg_title_base"], "Hazard ratio", SZ["c_note"], ha="center")
      n_bar = 128; pitch_bar = 0.96758; rw = 1.71
      for k in range(n_bar):
          h0 = math.exp(math.log(RR.RAMP_MIN) + (k / n_bar) * (math.log(RR.RAMP_MAX) - math.log(RR.RAMP_MIN))); h1 = math.exp(math.log(RR.RAMP_MIN) + ((k + 1) / n_bar) * (math.log(RR.RAMP_MAX) - math.log(RR.RAMP_MIN)))
          rect(bar0 + pitch_bar * k, Y["c_bar_y0"], rw, Y["c_bar_h"], RR.ramp_hex(math.sqrt(h0 * h1)), z=1)
      rect(bar0, Y["c_bar_y0"], BAR_W, Y["c_bar_h"], "none", ec="#ccd1d6", lw=0.4, z=2)
      BAR_TICKS = []
      for v in (1, 1.5, 2, 3):
          xt = bar0 + math.log(v) / math.log(RR.RAMP_MAX) * BAR_W; txt(xt, Y["c_bar_tick_base"], f"{v:g}", SZ["c_note"], ha="center"); BAR_TICKS.append(dict(value=v, x=round(xt, 3)))
          rec("c", f"{v:g}", xt, Y["c_bar_tick_base"], "center", SZ["c_note"], "ramp", "colour-bar tick (log scale 1 to 3.25 over the bar)", v)
      rect(PC["ref_x"], Y["c_ref_y0"], 14.5, 10.74, "white", ec=INK, lw=0.9, z=2)
      txt(PC["ref_x"] + 19.5, Y["c_ref_base"], "Reference, 6–7 h and T90 ≤1%", SZ["c_note"]); txt(PC["ref_x"], Y["c_ast_base"], "* q < 0.05, ** q < 0.01, *** q < 0.001", SZ["c_note"])
      raw = f"{WORK}/page_abc_raw.pdf"; fig.savefig(raw, facecolor="white"); plt.close(fig)

  if PASS == 1:
    _d = fitz.open(raw); _sp = C.spans_of(_d[0]); _d.close(); _worst = 0.0; _seen = {}
    for (key, xr, yr) in RIGHT:
        m = [t for t in _sp if t["text"] == key[0] and abs(t["origin"][1] - yr) < 0.06 and abs(t["size"] - key[1]) < 0.05]
        assert len(m) >= 1, ("pass-1 right-aligned string not found on the text layer", key, xr, yr)
        m = min(m, key=lambda t: abs(t["bbox"][2] - xr)); d = m["bbox"][2] - xr
        if key in _seen: assert abs(_seen[key] - d) < 0.05, ("pass 1: one string, two different edge corrections", key, _seen[key], d)
        else: _seen[key] = d; SHIFT[key] = d
        _worst = max(_worst, abs(d))
    print(f"pass 1: {len(RIGHT)} right-aligned strings measured ({len(SHIFT)} distinct text/size/weight keys), largest edge correction {_worst:.3f} pt; redrawing")
V8_2_METRICS = {k: tuple(v) for k, v in json.load(open(f"{WORK}/r38/v8_2_snapshot/Main_Fig3_drawn.json"))["font_metrics_aligned"].items()}   # the V13 sheet's descriptors (Arial-BoldMT 905/-211, ArialMT 1005/-324)
C.sheet_font_metrics = lambda sheet_pdf: (V8_2_METRICS if sheet_pdf == BASE else (_ for _ in ()).throw(RuntimeError("no text extraction on " + sheet_pdf)))
aligned = C.align_font_metrics(raw, f"{WORK}/page_abc.pdf", BASE)
d = fitz.open(f"{WORK}/page_abc.pdf"); pg = d[0]; assert abs(pg.rect.width - W) < 0.05 and abs(pg.rect.height - H_PAGE) < 0.05, pg.rect
_sp = C.spans_of(pg)
outside = [s["text"] for s in _sp if s["bbox"][1] < Y["title_base"] + 6.0 or s["bbox"][3] > H_PAGE - 2.0]; assert not outside, ("text outside the panel region", outside)
d.close()
n_marks = sum(1 for x in DRAWN if "marker" in x); assert n_marks == 104
out = {"figure": SHEET, "round": 54, "page_pdf": f"{WORK}/page_abc.pdf", "page": [W, H_PAGE], "band_top": BAND_TOP, "font_metrics_aligned": aligned, "plan": PLAN,
       "sources": {"bdsp_diseases_v3.csv": SC111, "hr_matrix_141x52.csv": SC168, "hr_matrix_ck_spo2_pct_below_90.json": {"path": MATRIX_CK, "sha256": C.sha256(MATRIX_CK), "role": "full precision of the step-168 fit (checkpoint, no sidecar of its own; equal to the csv at 3 dp on all 52 x 3 values, asserted)"},
                   "results.csv": SC136, "summary.json": SC136b, "NEW_Fig3_duration_desat_values.json": SC156, "t90_final.parquet": {"path": PARQ, "sha256": C.sha256(PARQ)}, "cohorts.json": SC101},
       "base_sheet": BASE, "base_sha256": C.sha256(BASE), "r49_record": R49_DRAWN, "r49_record_sha256": C.sha256(R49_DRAWN),
       "panel_a": {"n_rows": N_A, "n_per_column": N_COL, "order": [r["key"] for r in ROWS_A], "labels": LABELS, "rows_in": ROWS_IN, "rows_out": ROWS_OUT, "t90_sig": sorted(t90_sig), "tst_sig": sorted(tst_sig), "n_t90_sig": len(t90_sig), "n_tst_sig": len(tst_sig), "rows": ROWS_A,
                   "columns": [{k: v for k, v in col.items() if k != "rows"} | {"rows": [r["key"] for r in col["rows"]], "x_of_hr": f"{col['XA']:.4f} + {col['XB']:.4f} ln HR over [{col['xlo']}, {col['xhi']}]"} for col in COLS],
                   "geometry": {"row0_centre": ROW0, "pitch": PITCH, "lane": LANE, "axis_top": AX_TOP, "axis_bottom": AX_BOT, "tick_baseline": TICK_BASE, "caption_baseline": CAP_BASE, "header_baseline": HDR_BASE, "label_dy": LAB_DY, "hrq_dy": HRQ_DY, "scale_of_design": S_A, "box_w": PLAN["box_w"], "box_total": PLAN["box_total"], "pt_per_ln": PLAN["pt_per_ln"], "tick_gap": PLAN["tick_gap"]},
                   "key": {"baseline": Y["key_base"], "mark_y": Y["key_y"], "entries": KEY_POS_A, "rule": "one row beside the letter a, entries 20 pt apart, marks as V30 (line 15.94 long, marker on its centre, label 21.3 after the line start)"},
                   "rounding": {"csv_ties": tie_cells, "csv_halfup_vs_checkpoint_differ": tie_differ}, "tst_lane_csv_vs_matrix_max_diff": TST_DRIFT},
       "panel_b": {"rows": GRID_B, "ladder": GREEN_LADDER, "header_ladder": BLUE_LADDER, "min_pct": min_pct, "under_5": UNDER5, "not_shown": J["not_shown"], "tst_missing_masked": N_TST_MISSING, "key_bands": KEY_BANDS, "key_positions": KEY_POS, "key_bands_v13": ["5–15", "15–35", "≥35"],
                   "geometry": {"x0": BX0, "cells_x": HX, "cell_w": BCELL, "hdr_y0": Y["b_hdr_y0"], "hdr_h": Y["b_hdr_h"], "row_y0": Y["b_row_y0"], "row_h": Y["b_row_h"], "row_pitch": Y["b_row_pitch"], "label_right": B_LAB_R, "title_centre": B_CENTRE, "key_baseline": Y["b_key_base"]}},
       "panel_c": {"order": C_ORDER, "base_order": base_c_order, "control": {"key": CTRL_KEY, "label": CTRL_LABEL, "rule": "straight swap back pain -> alopecia (V8_1_DECISIONS item 2)", "ci_excludes_1_cells": ctrl_excl}, "rows": {t: C_ROWS[t] for t in C_ORDER}, "cell_n": cell_n, "rounding_ties": c_ties, "sn_ln_bh_significant": sn_ln_sig,
                   "ramp": {"rule": RR.__doc__.split("\n")[0], "pale": RR.RAMP_PALE, "top": RR.RAMP_TOP, "old_stops": RR.RAMP_OLD_STOPS},
                   "geometry": {"x0": CX0, "label_right": C_LAB_R, "groups_x": GX, "cell_w": CC, "pair_gap": PC["pair_gap"], "group_gap": PC["group_gap"], "hdr1_y0": Y["c_hdr1_y0"], "hdr1_h": Y["c_hdr1_h"], "hdr2_y0": Y["c_hdr2_y0"], "hdr2_h": Y["c_hdr2_h"], "row_y0": Y["c_row_y0"], "row_h": Y["c_row_h"], "row_pitch": Y["c_row_pitch"],
                                "line_dy_cells": Y["c_line_dy"], "cell_line_offsets": [C2_A, C2_B], "label_two_line_offsets": [L2_A, L2_B], "bar": [bar0, bar1, Y["c_bar_y0"], Y["c_bar_h"]], "bar_ticks": BAR_TICKS, "ref_box": [PC["ref_x"], Y["c_ref_y0"], 14.5, 10.74]}},
       "round54": {"sizes": SZ, "layout": "two side-by-side forest columns of 26 rows (own axis, HR and q columns each), key in one row beside the letter a, panels b and c side by side below the forest, band d below them",
                   "declared_duplicates": {"HR (95% CI)": 1, "q": 1, "Hazard ratio per 1 SD (95% CI)": 1, "0.6": 1, "0.8": 1, "1": 1, "1.5": 1},
                   "declared_sizes": {"panel a row labels, ticks, key labels (V30 11.71)": 14.0, "panel a HR and q strings (V30 11.17)": 14.0, "panel a column headers HR (95% CI), q (V30 11.17)": 15.0, "panel a caption (V30 12.78)": 15.0,
                                      "panel b title (V30 13.65)": 15.0, "panel b column headers (V30 12.5)": 15.0, "panel b cells, row labels, key (V30 12.5)": 14.0,
                                      "panel c duration headers (V30 10.91)": 15.0, "panel c row headers Sleep duration, Oxygenation (V30 10.5)": 15.0, "panel c row labels (V30 10.5)": 14.0,
                                      "panel c cell values and 1.00 (V30 10.5)": 13.0, "panel c sub-headers Preserved, Low, T90 ≤1%, T90 >10% (V30 10.5)": 13.0, "panel c notes: Hazard ratio, bar ticks, Reference, asterisk key (V30 10.09)": 13.0},
                   "exceptions": ["panel c cell values, the reference 1.00 and the four sub-headers at 13 pt (the round-54 small-note size) rather than 14: at 14 pt the interval strings are 68.4 pt wide and the six-cell grid beside the count grid does not fit the sheet width (planner: 4 x b cell + 6 x c cell budget)",
                                  "the right forest column's axis runs over [0.58, 1.55] (ticks 0.6, 0.8, 1, 1.5) because its 26 rows span 0.61 to 1.34, and the two boxes share one scale (the left box twice the right, one ln unit the same length in both, about 80 percent of the V30 length)"],
                   "text_layer": TXT},
       "drawn": DRAWN, "records": RECORDS, "written": C.now()}
json.dump(out, open(f"{WORK}/{SHEET}_drawn.json", "w"), indent=1, default=float)
print(f"wrote {WORK}/page_abc.pdf ({W:.2f} x {H_PAGE:.3f} pt, band top {BAND_TOP}) and {WORK}/{SHEET}_drawn.json ({len(DRAWN)} drawn, {len(RECORDS)} value records, {len(TXT)} text placements)")
import collections
print("round 54 size histogram of the placed text:", sorted(collections.Counter(t["size"] for t in TXT).items()))
