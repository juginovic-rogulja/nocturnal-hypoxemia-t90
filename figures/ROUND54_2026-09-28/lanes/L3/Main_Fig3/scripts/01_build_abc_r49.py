#!/usr/bin/env python3
"""ROUND 49 (2026-09-26, lane L3) copy of the round-40 builder: every text element +UP (1.0) pt inside the txt() helper, for the text and for
the width the placement rules see (right-aligned edge correction, the panel b key rung spacing). Geometry, data marks and printed strings
unchanged: the forest box left edge X_LEFT stays the round-40 value (the longest-label rule is evaluated at the OLD size, asserted against
the round-40 record), so two row labels whose end at +1 pt would cross the 3 pt clearance to the box stay at their current size (KEEP_LABELS,
computed by the builder's own rule). The panel c interval strings go to 10.5 pt (51.4 pt wide in a 53 pt cell, 0.8 pt margins, the
coordinator's decision of the follow-up: no text below 10 pt). Hand-placed
left-origin strings that the design centres in a cell or ends at a column edge keep that centre or edge (anchor=, a nudge of the origin by
half or all of the width growth, recorded in the drawn json). Records carry the size and origin as drawn. Original docstring follows.
ROUND 40 (2026-09-18, lane LMAIN) copy of the round-38 builder: ONE change, panel c moves DOWN by C_SHIFT so that the baseline of its
colour-bar tick labels (1, 1.5, 2, 3) equals panel a's x-axis tick-label baseline (TICK_BASE); the shift is one delta on every panel c
element (header rects and texts, row labels, cells, reference box, colour bar, key). Panel b stays (moving c down only widens the b-to-c
gap, 61.5 to 93.4 pt). Outputs go to the LMAIN lane; the V8.2 metrics snapshot is read from the lane's copy. Original docstring follows.
Main Fig 3, V14 rebuild of the three data panels as ONE matplotlib page pinned to the V13 sheet's page points (figure = page,
one axes in page coordinates, y down). Repointed copy of ROUND30_2026-09-08/V7_L3a_DURATION_MAIN/Main_Fig3/scripts/01_build_abc.py
(backup 01_build_abc_PRE_V8_1.py). Panel a = T90 against total sleep time per 1 SD for the 52 ranked outcomes (figure2C_polish.py
wiring, the sheet carries that design at scale 1.0706), panel b = the sleep-duration count grid (ROUND18 j1 lineage), panel c = the
six-cell hazard grid (ROUND15 build_r15_fig3b lineage, round-27 red ramp and white-numeral rule).
V8.1 changes against the R30 copy: (1) panel a has 52 rows (the 48 ranked outcomes plus polycythemia, Parkinson's disease, interstitial
lung disease, ventricular arrhythmia or cardiac arrest): the row pitch is kept, the forest grows by 4 rows at the foot, the axis rule,
tick labels and caption move down by the same amount (DELTA = 4 x 12.719 pt) and the page grows by DELTA (02_compose shifts the kept
letter-d strip and band d down by DELTA); (2) the T90 series is bdsp_diseases_v3.csv (step 111, all 19,173 nights), the printed strings
half up on the full-precision checkpoint of the same fit (dC_side_analyses/_work/hr_matrix_ck/spo2_pct_below_90.json, written by step
168 together with hr_matrix_141x52.csv, asserted equal to the csv at 3 decimals on every value); the total sleep time series is the csv's
tst columns (fitted on the 15,551 full nights, markers only); q = Benjamini-Hochberg across the 52 within the exposure from the csv's p,
asserted equal to the matrix's q sets; (3) panel b from NEW_Fig3_duration_desat_values.json (step 156, all 19,173 rows with the masked
sleep amounts, so the three rows are the full nights) re-derived from the v8 table; (4) panel c from dur_x_oxygen_6cell_tst300 (step
136, 15,551) with the control row ALOPECIA (the straight swap for back pain of 14 September); the SN and LN cells are reported, not
asserted; (5) drawn records carry the v14lib DRAWN_RECORD fields. Output: work/page_abc.pdf (page W x (H + DELTA), font metrics
aligned to the V13 sheet) and work/Main_Fig3_drawn.json. The title strip, the letter-d strip and band d are NOT drawn here."""
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
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # round 49: this lane's scripts folder (lane_common, ramp_rule, v13_cache copies)
import v13_cache
import lane_common as C
import ramp_rule as RR
VL = C.VL

SHEET = "Main_Fig3"; HERE = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/L3/Main_Fig3"; WORK = f"{HERE}/work"; VER = f"{HERE}/verify"
UP = 1.0                                                       # round 49: every text element one point larger
R40_DRAWN = f"{paths.FIGURE_ROOT}/ROUND40_2026-09-18/lanes/LMAIN/Main_Fig3/work/Main_Fig3_drawn.json"   # the frozen geometry (read only)
AT = {}; NUDGE = []; KEEP = []                                 # round 49: what txt() drew per (text, given x, baseline); the nudged and the kept elements
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
# round 38: step 168 is RESUMABLE (one checkpoint per measurement in _work/hr_matrix_ck/, a rerun reuses them), so the checkpoint's mtime
# (2026-09-14 15:37) predates the 2026-09-16 17:28 rerun that assembled the full-precision matrix from it; the round-37 mtime guard
# (checkpoint within [-3 h, +1 h] of the matrix) is replaced by a content guard below (_CK_GUARD): every HR, lo, hi and p of the checkpoint
# equals the matrix row within 1e-9.
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
W, H = v13_cache.page("Main_Fig3"); BSP = C.dedupe(v13_cache.spans_lane("Main_Fig3"))   # round 38: the V13 text layer from the probe cache, no PyMuPDF text extraction on the nested base
assert abs(W - 952.73) < 0.01 and abs(H - 992.945) < 0.01, (W, H)
DRAWN = []; RECORDS = []
def _as_drawn(text, x, y, size=None):
    """round 49: the origin and size txt() actually drew for this string (nudged origin, +UP or kept size)."""
    d = AT.get((text, round(x, 2), round(y, 2)))
    return (d[0], d[1]) if d else (x, size)
def rec(panel, text, x, y, ha, source, key, value, **kw):
    x, size = _as_drawn(text, x, y)
    DRAWN.append(dict(panel=panel, text=text, x=round(x, 3), baseline=round(y, 3), ha=ha, source=source, key=key, value=value, **({"size": size} if size else {}), **kw))
def vrec(panel, text, x, y, ha, size, source_file, source_key, source_value, rule, note=""):
    x, size = _as_drawn(text, x, y, size)
    RECORDS.append(dict(panel=panel, text=text, x=round(x, 3), baseline=round(y, 3), ha=ha, size=size, source_file=source_file, source_key=source_key, source_value=source_value, rule=rule, note=note))


def bh(p):
    p = np.asarray(p, float); m = len(p); o = np.argsort(p); ranks = np.empty(m); ranks[o] = np.arange(1, m + 1)
    q = p * m / ranks; qs = q[o]; qs = np.minimum.accumulate(qs[::-1])[::-1]; out = np.empty(m); out[o] = np.minimum(qs, 1.0); return out


def q_text(q): return VL.q_text(q)          # three decimals half up, a fourth where three would cross the 0.05 boundary
def fmt2_halfup(x): return VL.halfup(x, 2)
def stars(q): return "" if (q is None or not np.isfinite(q)) else ("***" if q < 0.001 else "**" if q < 0.01 else "*" if q < 0.05 else "")
LABEL_SHORT = {"Ventricular arrhythmia or cardiac arrest": "Ventricular arrhythmia/arrest"}   # the lane L2 convention on Fig 2a (owner question): the full name is 185 pt at 10.71 pt and cannot fit the 129 pt gutter
def label_of(k): return LABEL_SHORT.get(DISEASES[k][0], DISEASES[k][0]) if k in DISEASES else "Death from any cause"

# =============================================================== panel a data
D = pd.read_csv(CSV, comment="#").set_index("key")
M = pd.read_csv(MATRIX)
T = M[M.measurement == "spo2_pct_below_90"].set_index("outcome"); S = M[M.measurement == "TST_min"].set_index("outcome")
CK = {r["outcome"]: r for r in json.load(open(MATRIX_CK))}; assert all(r["measurement"] == "spo2_pct_below_90" for r in CK.values())
N_A = len(T); assert N_A == 52 and len(S) == 52 and sorted(T.index) == sorted(S.index) == sorted(CK), (N_A, len(S), len(CK))
KEYS = list(T.index); assert all(k in D.index for k in KEYS), [k for k in KEYS if k not in D.index]
_CK_GUARD = max(abs(float(T.loc[k, c]) - float(CK[k][c])) for k in KEYS for c in ("HR", "lo", "hi", "p"))
assert _CK_GUARD < 1e-9, ("round 38: the checkpoint is not the matrix's own values", _CK_GUARD)
print(f"checkpoint = matrix on every HR, lo, hi, p of the {len(KEYS)} outcomes (max |diff| {_CK_GUARD:.1e}; step 168 rerun 2026-09-16 17:28 resumed from the checkpoints)")
CTRL_IN_A = sorted(k for k in KEYS if k in NEGATIVE_CONTROLS)      # the ranked set (>= 150 events) includes the qualifying controls, as the V13 panel did (glaucoma, back pain, cataract in v7)
print(f"panel a: negative controls among the {len(KEYS)} ranked outcomes: {CTRL_IN_A}")
for k in KEYS:      # the csv (step 111) is the published panel; the matrix (step 168) and its checkpoint are the same fit at full precision
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
XLO_V13, XHI = 0.74, 4.15
XLO = 0.58                                                     # V8.1: Parkinson's disease T90 HR 0.72 (0.61-0.86) enters the panel; the axis extends to 0.58 with a 0.6 tick (the lane L2 rule on Fig 2a), the forest box unchanged
_allx = [CK[k][c] for k in KEYS for c in ("HR", "lo", "hi")] + [Dk.loc[k, c] for k in KEYS for c in ("tst_hr", "tst_lo", "tst_hi")]
assert XLO < min(_allx) and max(_allx) < XHI, ("a value leaves the sheet's axis", min(_allx), max(_allx))
_lab_font = fitz.Font(fontfile=C.ARIAL); _lab_end = {label_of(k): 8.88 + _lab_font.text_length(label_of(k), fontsize=10.71) for k in KEYS}
LAB_END_MAX = max(_lab_end.values()); LAB_LONGEST = max(_lab_end.items(), key=lambda t: t[1])[0]
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

# =============================================================== panel b data
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

# =============================================================== panel c data
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

# =============================================================== geometry (page points, y down); V8.1: the forest grows by DELTA at the foot
XA_V13, XB_V13 = 180.46, 139.54                                # V13: x = 180.46 + 139.54 ln HR over [0.74, 4.15]
X_LEFT_V13, X_RIGHT = XA_V13 + XB_V13 * math.log(XLO_V13), XA_V13 + XB_V13 * math.log(XHI)  # the V13 forest box
X_LEFT = max(X_LEFT_V13, math.ceil((LAB_END_MAX + 3.0) * 10) / 10)                          # V8.1: the left edge clears the longest row label by 3 pt (the right edge is kept)
XB = (X_RIGHT - X_LEFT) / math.log(XHI / XLO); XA = X_LEFT - XB * math.log(XLO)             # the box over [0.58, 4.15]
_R40 = json.load(open(C.hydrated(R40_DRAWN)))["panel_a"]["geometry"]                        # round 49: geometry frozen, the rule above is evaluated at the OLD label size
assert abs(X_LEFT - float(_R40["x_left"])) < 1e-6 and abs(XA - float(_R40["XA"])) < 1e-6 and abs(XB - float(_R40["XB"])) < 1e-6, ("round 49: the forest box moved", X_LEFT, XA, XB, _R40)
print(f"panel a axis: box x {X_LEFT:.2f} to {X_RIGHT:.2f} (V13 {X_LEFT_V13:.2f} to {X_RIGHT:.2f}), longest label '{LAB_LONGEST}' ends at {LAB_END_MAX:.1f}; x = {XA:.3f} + {XB:.3f} ln HR")
def xhr(v): return XA + XB * math.log(v)
PITCH = 12.719; ROW0 = 116.74; LANE = 3.06; N_BASE = 48
DELTA = round((N_A - N_BASE) * PITCH, 3)
AX_TOP = ROW0 - 0.75 * PITCH; AX_BOT = ROW0 + (N_A - 0.25) * PITCH
H_NEW = H + DELTA; Y_D_NEW = 772.0 + DELTA
TICKS_A = [0.6, 0.8, 1, 1.5, 2, 3, 4]                          # V8.1: the 0.6 tick added (axis extended to 0.58)
S_A = 1.0706                                                   # the sheet's scale of the figure2C_polish design
FS_ROW, FS_COL, FS_CAP = 10.71, 10.17, 11.78
assert abs(FS_ROW - 10.71) < 1e-9, "the longest-label rule above measures at 10.71 (the row label size)"
KEEP_LABELS = [lab for lab in LABELS if 8.88 + _lab_font.text_length(lab, fontsize=FS_ROW + UP) > X_LEFT - 3.0]   # round 49: the builder's clearance rule at +UP
print(f"round 49: row labels that would cross the 3 pt clearance to the forest box at {FS_ROW + UP} pt, kept at {FS_ROW} pt: {KEEP_LABELS}")
S_CIRC, S_SQ = (math.sqrt(28.0) * S_A) ** 2, (math.sqrt(23.0) * S_A) ** 2
LW_CI, LW_MEDGE, LW_OPEN = 1.7 * S_A, 0.8 * S_A, 1.1 * S_A
KEY_ROWS = [(54.68, 58.25, "T90, percent of the recording below 90% saturation"), (69.02, 72.67, "Total sleep time on the study night"), (83.60, 86.97, "Not significant")]
TICK_BASE, CAP_BASE = 738.63 + DELTA, 754.69 + DELTA
C_TICK_BASE_V13 = 757.61                                       # panel c colour-bar tick labels (1, 1.5, 2, 3) baseline in the V13 design
C_SHIFT = round(TICK_BASE - C_TICK_BASE_V13, 3)                # round 40: panel c moves down so its colour-bar tick baseline = panel a's tick baseline
def cy(y): return y + C_SHIFT                                  # every panel c y goes through cy()

RC = {"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.bbox": "standard", "savefig.pad_inches": 0.0, "figure.dpi": 100,
      "text.color": INK, "axes.unicode_minus": False, "lines.scale_dashes": False, "lines.solid_capstyle": "butt", "lines.dash_capstyle": "butt"}
FONT_R = fitz.Font(fontfile=C.ARIAL); FONT_B = fitz.Font(fontfile=C.ARIALB); SHIFT = {}; RIGHT = []
for PASS in (1, 2):
  DRAWN.clear(); RECORDS.clear(); RIGHT.clear(); AT.clear(); NUDGE.clear(); KEEP.clear()
  with plt.rc_context(RC):
      fig = plt.figure(figsize=(W / 72.0, H_NEW / 72.0)); ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(H_NEW, 0); ax.axis("off")
      def txt(x, y, s, size, ha="left", bold=False, color=INK, z=5, anchor=None, up=True):
          """round 49: size + UP for the text and for every width a placement rule measures. anchor="center" or "right" keeps a hand-placed
          left-origin string's centre or right edge (the origin moves by half or all of the width growth, recorded in NUDGE). up=False keeps
          one element at its current size (recorded in KEEP). AT records the origin and size drawn for rec()/vrec()."""
          x_given, size0 = x, size; size = size + UP if up else size; f = FONT_B if bold else FONT_R
          if ha == "left" and anchor in ("center", "right") and size != size0:
              dw = f.text_length(s, fontsize=size) - f.text_length(s, fontsize=size0); dx = -dw / 2.0 if anchor == "center" else -dw
              NUDGE.append(dict(text=s, x_given=round(x, 3), x_drawn=round(x + dx, 3), baseline=round(y, 3), anchor=anchor, dx=round(dx, 3), size=size)); x = x + dx
          if not up: KEEP.append(dict(text=s, x=round(x, 3), baseline=round(y, 3), size=size, ha=ha))
          AT[(s, round(x_given, 2), round(y, 2))] = (round(x, 3), size)
          if ha == "right":
              key = (s, round(y, 2), size, bold); adv = f.text_length(s, fontsize=size)
              RIGHT.append((key, x)); ax.text(x - adv - SHIFT.get(key, 0.0), y, s, fontsize=size, ha="left", va="baseline", color=color, fontweight="bold" if bold else "normal", zorder=z); return
          ax.text(x, y, s, fontsize=size, ha=ha, va="baseline", color=color, fontweight="bold" if bold else "normal", zorder=z)
      def rect(x, y, w, h, fc, ec="none", lw=0.0, z=1):
          ax.add_patch(Rectangle((x, y), w, h, facecolor=fc, edgecolor=ec, linewidth=lw, zorder=z))
      # ---------------- panel a
      for v in TICKS_A:
          if v != 1: ax.plot([xhr(v)] * 2, [AX_TOP, AX_BOT], color=GRIDC, lw=0.5 * S_A, zorder=0.3)
      ax.plot([xhr(1.0)] * 2, [AX_BOT, AX_TOP], color=INK, lw=0.9 * S_A, ls=(0, (4 * 0.9 * S_A, 3 * 0.9 * S_A)), zorder=1)
      ax.plot([xhr(XLO), xhr(XHI)], [AX_BOT] * 2, color=INK, lw=0.8 * S_A, solid_capstyle="projecting", zorder=4)
      for v in TICKS_A:
          ax.plot([xhr(v)] * 2, [AX_BOT, AX_BOT + 3.0 * S_A], color=INK, lw=0.8 * S_A, zorder=4)
          txt(xhr(v), TICK_BASE, f"{v:g}", FS_ROW, ha="center"); rec("a", f"{v:g}", xhr(v), TICK_BASE, "center", "axis", "tick", v)
      txt(176.82, CAP_BASE, "Hazard ratio per 1 SD (95% CI)", FS_CAP, anchor="center")
      txt(460.86, 99.47, "HR (95% CI)", FS_COL, ha="right"); txt(500.10, 99.60, "q", FS_COL, ha="right")
      for i, r in enumerate(ROWS_A):
          yc = ROW0 + PITCH * i; yb, yg = yc - LANE, yc + LANE
          ax.plot([xhr(r["lo"]), xhr(r["hi"])], [yb, yb], color=BLUE, lw=LW_CI, solid_capstyle="round", zorder=2)
          ax.plot([xhr(r["tst_lo"]), xhr(r["tst_hi"])], [yg, yg], color=GREEN, lw=LW_CI, solid_capstyle="round", zorder=2)
          if r["sig"]: ax.scatter([xhr(r["hr"])], [yb], s=S_CIRC, marker="o", facecolors=BLUE, edgecolors="white", linewidths=LW_MEDGE, zorder=3)
          else: ax.scatter([xhr(r["hr"])], [yb], s=S_CIRC, marker="o", facecolors="white", edgecolors=BLUE, linewidths=LW_OPEN, zorder=3)
          if r["tst_sig"]: ax.scatter([xhr(r["tst_hr"])], [yg], s=S_SQ, marker="s", facecolors=GREEN, edgecolors="white", linewidths=LW_MEDGE, zorder=3)
          else: ax.scatter([xhr(r["tst_hr"])], [yg], s=S_SQ, marker="s", facecolors="white", edgecolors=GREEN, linewidths=LW_OPEN, zorder=3)
          txt(8.88, yc + 3.81, r["label"], FS_ROW, up=r["label"] not in KEEP_LABELS); txt(460.87, yc - 0.37, r["s_hr"], FS_COL, ha="right"); txt(500.11, yc - 0.37, r["s_q"], FS_COL, ha="right", bold=r["sig"])
          rec("a", r["label"], 8.88, yc + 3.81, "left", MATRIX_CK, f"[outcome={r['key']}] (row {i + 1}, sorted by the T90 hazard ratio)", None)
          rec("a", r["s_hr"], 460.87, yc - 0.37, "right", MATRIX_CK, f"[outcome={r['key']}]: HR, lo, hi (full precision; csv t90_hr/lo/hi equal at 3 decimals)", [r["hr"], r["lo"], r["hi"]], csv=r["csv_t90"], sidecar=SC168["sha256"])
          rec("a", r["s_q"], 500.11, yc - 0.37, "right", MATRIX, f"measurement=spo2_pct_below_90&outcome={r['key']}:q (BH across the 52 within the exposure; csv-derived q {r['q_csv']:.6g})", r["q"], bold=r["sig"], sidecar=SC168["sha256"])
          vrec("a", r["label"], 8.88, yc + 3.81, "left", FS_ROW, CSV, f"key={r['key']}:disease", r["label"], "label_alias", f"row {i + 1}")
          vrec("a", r["s_hr"], 460.87, yc - 0.37, "right", FS_COL, MATRIX_CK, f"[outcome={r['key']}]", {"HR": r["hr"], "lo": r["lo"], "hi": r["hi"]}, "hr_ci_HRlohi", "half up on the step-168 full-precision checkpoint; bdsp_diseases_v3.csv equal at 3 dp")
          vrec("a", r["s_q"], 500.11, yc - 0.37, "right", FS_COL, MATRIX, f"measurement=spo2_pct_below_90&outcome={r['key']}:q", r["q"], "q3", "BH across the 52 ranked outcomes within the exposure (= the csv-derived q)")
          DRAWN.append(dict(panel="a", marker="T90 circle", filled=r["sig"], x=round(xhr(r["hr"]), 3), y=round(yb, 3), ci_x=[round(xhr(r["lo"]), 3), round(xhr(r["hi"]), 3)], key=r["key"]))
          DRAWN.append(dict(panel="a", marker="total sleep time square", filled=r["tst_sig"], x=round(xhr(r["tst_hr"]), 3), y=round(yg, 3), ci_x=[round(xhr(r["tst_lo"]), 3), round(xhr(r["tst_hi"]), 3)], key=r["key"],
                            source=CSV, value=[r["tst_hr"], r["tst_lo"], r["tst_hi"]], q=r["tst_q"], sidecar=SC111["sha256"]))
      for (yk, ybase, s), kind in zip(KEY_ROWS, ("blue", "green", "open")):
          if kind == "blue":
              ax.plot([18.63, 34.57], [yk, yk], color=BLUE, lw=LW_CI, solid_capstyle="round", zorder=2)
              ax.scatter([26.60], [yk], s=(6.0 * S_A) ** 2, marker="o", facecolors=BLUE, edgecolors="white", linewidths=LW_MEDGE, zorder=3); txt(39.93, ybase, s, FS_ROW)
          elif kind == "green":
              ax.plot([18.63, 34.57], [yk, yk], color=GREEN, lw=LW_CI, solid_capstyle="round", zorder=2)
              ax.scatter([26.60], [yk], s=(5.4 * S_A) ** 2, marker="s", facecolors=GREEN, edgecolors="white", linewidths=LW_MEDGE, zorder=3); txt(39.93, ybase, s, FS_ROW)
          else:
              ax.scatter([26.66], [yk], s=(6.0 * S_A) ** 2, marker="o", facecolors="white", edgecolors=GREEN, linewidths=LW_OPEN, zorder=3); txt(39.96, ybase, s, FS_ROW)
      # ---------------- panel b
      txt(712.0, 55.19, "T90, % of the recording", 12.65, anchor="center")
      HX = [610.76 + 85.29 * c for c in range(4)]; HW = 78.09
      for c in range(4):
          rect(HX[c], 71.14, HW, 16.0, BLUE_LADDER[c])
      for c, (s, x0) in enumerate(zip(BANDS, (638.34, 717.03, 799.13, 890.82))): txt(x0, 82.13, s, 11.5, anchor="center")
      for r_i, r_ in enumerate(DUR):
          y0 = 98.14 + 63.245 * r_i; yb_ = 130.45 + 63.245 * r_i
          txt(599.73, yb_, r_, 11.5, ha="right")
          for c in range(4):
              cell = GRID_B[r_]["cells"][c]; fc = GREEN_LADDER[cell["shade"]]
              rect(HX[c], y0, HW, 58.64, fc); txt(HX[c] + HW / 2, yb_, cell["label"], 11.5, ha="center")
              rec("b", cell["label"], HX[c] + HW / 2, yb_, "center", PARQ, f"{r_} x {cell['band']}: n and percent of the duration group (direct recompute = step-156 file)", [cell["n"], cell["pct"]], fill=fc, band=cell["shade"], sidecar=SC156["sha256"])
              vrec("b", cell["label"], HX[c] + HW / 2, yb_, "center", 11.5, J1, f"rows/{r_}/cells/{c}/label", cell["label"], "text", f"n {cell['n']} of {GRID_B[r_]['total']}, recomputed from the v8 table")
      txt(562.81, 312.47, "Percent of the group", 11.5)
      # V8.1: the key draws a rung for every ladder band that a cell of the grid uses, in the ladder order of the step-156 file
      # (shade_bands_pct). V13 printed three rungs because no cell was under 5 percent; the >7 h x >5-10% cell is now 3.9 percent,
      # so the "<5" rung returns (the design's own rule, as Fig 4c keeps it). The V13 rhythm: swatch 12.42 x 10.76 at y 304.39 from
      # x 685.14, label 17.42 pt after the swatch's left edge, the next swatch 20.0 pt after the label's end (V13: 685.14/702.56,
      # 748.14/765.56, 817.54/834.96).
      KEY_BANDS = [b_ for b_ in J["shade_bands_pct"] if any(c["shade"] == b_ for r_ in DUR for c in GRID_B[r_]["cells"])]
      assert KEY_BANDS and all(b_ in GREEN_LADDER for b_ in KEY_BANDS), KEY_BANDS
      _kx = 685.14; KEY_POS = []
      for lab in KEY_BANDS:
          xl = _kx + 17.42; rect(_kx, 304.39, 12.42, 10.76, GREEN_LADDER[lab], ec="white", lw=0.6, z=2); txt(xl, 312.47, lab, 11.5)
          KEY_POS.append(dict(band=lab, swatch_x=round(_kx, 3), label_x=round(xl, 3)))
          vrec("b", lab, xl, 312.47, "left", 11.5, J1, f"shade_bands_pct/{J['shade_bands_pct'].index(lab)}", lab, "text", "key rung, drawn because a cell of the grid falls in this band")
          _kx = xl + FONT_R.text_length(lab, fontsize=11.5 + UP) + 20.0   # round 49: the rung rule sees the larger label
      assert _kx < W - 20, ("the key runs off the sheet", _kx)
      # ---------------- panel c
      GX = [601.11 + 118.01 * g for g in range(3)]
      for g in range(3): rect(GX[g], cy(376.60), 107.6, 24.0, GREEN3[g])
      for s, x0 in zip(DUR, (645.13, 760.52, 881.15)): txt(x0, cy(391.18), s, 9.91, anchor="center")
      for g in range(3):
          rect(GX[g], cy(401.60), 53.0, 34.0, BLUE_PRES); rect(GX[g] + 54.6, cy(401.60), 53.0, 34.0, BLUE_LOW)
          txt(605.96 + 118.01 * g, cy(414.61), "Preserved", 9.5, anchor="center"); txt(673.50 + 118.01 * g, cy(414.61), "Low", 9.5, anchor="center")
          txt(608.64 + 118.01 * g, cy(426.20), "T90 ≤1%", 9.5, anchor="center"); txt(660.43 + 118.01 * g, cy(426.20), "T90 >10%", 9.5, anchor="center")
      txt(530.85, cy(391.07), "Sleep duration", 9.5, anchor="right"); txt(538.77, cy(421.07), "Oxygenation", 9.5, anchor="right")
      rect(719.12, cy(438.40), 53.0, 267.95, "white", ec=INK, lw=0.9, z=2)
      TWO = {"Cardiovascular composite": ("Cardiovascular", "composite"), CTRL_LABEL: (DISEASES[CTRL_KEY][0], "(negative control)")}
      for j, title in enumerate(C_ORDER):
          y0 = cy(438.40) + 38.55 * j; yc = y0 + 18.325; kk = C_ROWS[title]["key"]
          if title in TWO:
              txt(592.08, yc - 3.995, TWO[title][0], 9.5, ha="right"); txt(592.08, yc + 7.595, TWO[title][1], 9.5, ha="right")
              vrec("c", TWO[title][0], 592.08, yc - 3.995, "right", 9.5, RES, f"analysis=primary&key={kk}:disease", C_ROWS[title]["disease"], "label_alias", f"row {j + 1}, first line")
              vrec("c", TWO[title][1], 592.08, yc + 7.595, "right", 9.5, RES, f"analysis=primary&key={kk}:negative_control", bool(C_ROWS[title]["control"]), "label_alias", f"row {j + 1}, second line")
          else:
              txt(592.08, yc + 1.805, title, 9.5, ha="right")
              vrec("c", title, 592.08, yc + 1.805, "right", 9.5, RES, f"analysis=primary&key={kk}:disease", C_ROWS[title]["disease"], "label_alias", f"row {j + 1}")
          rec("c", title, 592.08, yc + 1.805, "right", RES, f"key={kk} (row {j + 1}, rows by the largest hazard ratio, control last)", None)
          txt(745.62, yc + 1.805, "1.00", 9.5, ha="center"); rec("c", "1.00", 745.62, yc + 1.805, "center", RES, "reference cell NN (6–7 h, T90 ≤1%)", 1.0)
          vrec("c", "1.00", 745.62, yc + 1.805, "center", 9.5, SUMM, "thresholds/reference_cell", "NN", "ref_one", "the reference cell prints 1.00 by definition")
          for c in CELLS:
              g, l = CELL_POS[c]; x0 = GX[g] + 54.6 * l; cell = C_ROWS[title]["cells"][c]
              rect(x0, y0, 53.0, 36.65, cell["fill"], z=1)
              txt(x0 + 26.5, yc - 3.995, cell["s_hr"], 9.5, ha="center", color=cell["ink"]); txt(x0 + 26.5, yc + 7.595, cell["s_ci"], 9.5, ha="center", color=cell["ink"])   # round 49 follow-up: 10.5 pt, 51.4 pt wide in the 53 pt cell (0.8 pt margins, the coordinator's decision: no text below 10 pt)
              rec("c", cell["s_hr"], x0 + 26.5, yc - 3.995, "center", RES, f"analysis=primary key={kk} {c}_hr with stars from {c}_q", [cell["hr"], cell["q"]], fill=cell["fill"], ink=cell["ink"], sidecar=SC136["sha256"])
              rec("c", cell["s_ci"], x0 + 26.5, yc + 7.595, "center", RES, f"analysis=primary key={kk} {c}_lo, {c}_hi", [cell["lo"], cell["hi"]], fill=cell["fill"], ink=cell["ink"], sidecar=SC136["sha256"])
              vrec("c", cell["s_hr"], x0 + 26.5, yc - 3.995, "center", 9.5, RES, f"analysis=primary&key={kk}:{c}_hr,{c}_q", [cell["hr"], cell["q"]], "hr_stars", "half up on the stored 3 dp value; stars by the cell's BH q (none on the control row)")
              vrec("c", cell["s_ci"], x0 + 26.5, yc + 7.595, "center", 9.5, RES, f"analysis=primary&key={kk}:{c}_lo,{c}_hi", [cell["lo"], cell["hi"]], "ci_paren", "half up on the stored 3 dp value")
      txt(663.40, cy(733.61), "Hazard ratio", 9.09, anchor="center")
      n_bar = 128; bx0, pitch_bar, rw, by0, bh_ = 626.028, 0.96758, 1.71, cy(738.612), 9.914
      for k in range(n_bar):
          h0 = math.exp(math.log(RR.RAMP_MIN) + (k / n_bar) * (math.log(RR.RAMP_MAX) - math.log(RR.RAMP_MIN))); h1 = math.exp(math.log(RR.RAMP_MIN) + ((k + 1) / n_bar) * (math.log(RR.RAMP_MAX) - math.log(RR.RAMP_MIN)))
          rect(bx0 + pitch_bar * k, by0, rw, bh_, RR.ramp_hex(math.sqrt(h0 * h1)), z=1)
      rect(626.03, cy(738.61), 750.27 - 626.03, 748.53 - 738.61, "none", ec="#ccd1d6", lw=0.4, z=2)
      for s, x0 in (("1", 623.50), ("1.5", 662.33), ("2", 696.37), ("3", 739.00)): txt(x0, cy(C_TICK_BASE_V13), s, 9.09, anchor="center")
      rect(776.27, cy(726.80), 14.5, 10.74, "white", ec=INK, lw=0.9, z=2)
      txt(795.77, cy(735.44), "Reference, 6–7 h and T90 ≤1%", 9.09); txt(776.27, cy(749.35), "* q < 0.05, ** q < 0.01, *** q < 0.001", 9.09)
      raw = f"{WORK}/page_abc_raw.pdf"; fig.savefig(raw, facecolor="white"); plt.close(fig)

  if PASS == 1:
    _d = fitz.open(raw); _sp = C.spans_of(_d[0]); _d.close(); _worst = 0.0
    for (key, xr) in RIGHT:
        m = [t for t in _sp if t["text"] == key[0] and abs(t["origin"][1] - key[1]) < 0.06 and abs(t["size"] - key[2]) < 0.05]
        assert len(m) >= 1, ("pass-1 right-aligned string not found on the text layer", key)
        m = min(m, key=lambda t: abs(t["bbox"][2] - xr)); SHIFT[key] = SHIFT.get(key, 0.0) + (m["bbox"][2] - xr); _worst = max(_worst, abs(m["bbox"][2] - xr))
    print(f"pass 1: {len(RIGHT)} right-aligned strings measured, largest edge correction {_worst:.3f} pt; redrawing")
V8_2_METRICS = {k: tuple(v) for k, v in json.load(open(f"{WORK}/r38/v8_2_snapshot/Main_Fig3_drawn.json"))["font_metrics_aligned"].items()}   # round 38: the V13 sheet's descriptors as the round-37 build recorded them (Arial-BoldMT 905/-211, ArialMT 1005/-324)
C.sheet_font_metrics = lambda sheet_pdf: (V8_2_METRICS if sheet_pdf == BASE else (_ for _ in ()).throw(RuntimeError("no text extraction on " + sheet_pdf)))   # round 38: no get_text on the nested V13 base
aligned = C.align_font_metrics(raw, f"{WORK}/page_abc.pdf", BASE)
d = fitz.open(f"{WORK}/page_abc.pdf"); pg = d[0]; assert abs(pg.rect.width - W) < 0.05 and abs(pg.rect.height - H_NEW) < 0.05, pg.rect
outside = [s["text"] for s in C.spans_of(pg) if s["bbox"][1] < 16.0 or s["bbox"][3] > Y_D_NEW]; assert not outside, ("text outside the rebuilt region", outside)
d.close()
out = {"figure": SHEET, "page_pdf": f"{WORK}/page_abc.pdf", "page": [W, H_NEW], "page_v13": [W, H], "delta_pt": DELTA, "y_d_new": Y_D_NEW, "font_metrics_aligned": aligned,
       "panel_c_shift_pt": C_SHIFT, "panel_c_shift_rule": f"round 40: panel c down by {C_SHIFT} pt so its colour-bar tick baseline ({C_TICK_BASE_V13} + shift = {cy(C_TICK_BASE_V13):.3f}) equals panel a's tick baseline ({TICK_BASE:.3f}); panel b unchanged",
       "sources": {"bdsp_diseases_v3.csv": SC111, "hr_matrix_141x52.csv": SC168, "hr_matrix_ck_spo2_pct_below_90.json": {"path": MATRIX_CK, "sha256": C.sha256(MATRIX_CK), "role": "full precision of the step-168 fit (checkpoint, no sidecar of its own; equal to the csv at 3 dp on all 52 x 3 values, asserted)"},
                   "results.csv": SC136, "summary.json": SC136b, "NEW_Fig3_duration_desat_values.json": SC156, "t90_final.parquet": {"path": PARQ, "sha256": C.sha256(PARQ)}, "cohorts.json": SC101},
       "base_sheet": BASE, "base_sha256": C.sha256(BASE),
       "panel_a": {"n_rows": N_A, "order": [r["key"] for r in ROWS_A], "labels": LABELS, "rows_in": ROWS_IN, "rows_out": ROWS_OUT, "t90_sig": sorted(t90_sig), "tst_sig": sorted(tst_sig), "n_t90_sig": len(t90_sig), "n_tst_sig": len(tst_sig), "rows": ROWS_A,
                   "geometry": {"x_of_hr": f"{XA:.4f} + {XB:.4f} ln HR (V13: 180.46 + 139.54 ln HR over [0.74, 4.15]; extended to [0.58, 4.15] on the same box)", "XA": XA, "XB": XB, "x_left": X_LEFT, "x_right": X_RIGHT, "x_left_v13": X_LEFT_V13, "longest_label": [LAB_LONGEST, LAB_END_MAX], "row0_centre": ROW0, "pitch": PITCH, "lane": LANE, "axis_top": AX_TOP, "axis_bottom": AX_BOT, "scale_of_design": S_A, "xlim": [XLO, XHI], "ticks": TICKS_A, "tick_baseline": TICK_BASE, "caption_baseline": CAP_BASE},
                   "rounding": {"csv_ties": tie_cells, "csv_halfup_vs_checkpoint_differ": tie_differ}, "tst_lane_csv_vs_matrix_max_diff": TST_DRIFT},
       "panel_b": {"rows": GRID_B, "ladder": GREEN_LADDER, "header_ladder": BLUE_LADDER, "min_pct": min_pct, "under_5": UNDER5, "not_shown": J["not_shown"], "tst_missing_masked": N_TST_MISSING, "key_bands": KEY_BANDS, "key_positions": KEY_POS, "key_bands_v13": ["5–15", "15–35", "≥35"]},
       "panel_c": {"order": C_ORDER, "base_order": base_c_order, "control": {"key": CTRL_KEY, "label": CTRL_LABEL, "rule": "straight swap back pain -> alopecia (V8_1_DECISIONS item 2)", "ci_excludes_1_cells": ctrl_excl}, "rows": {t: C_ROWS[t] for t in C_ORDER}, "cell_n": cell_n, "rounding_ties": c_ties, "sn_ln_bh_significant": sn_ln_sig,
                   "ramp": {"rule": RR.__doc__.split("\n")[0], "pale": RR.RAMP_PALE, "top": RR.RAMP_TOP, "old_stops": RR.RAMP_OLD_STOPS}},
       "round49": {"up_pt": UP, "rule": "every text element + UP pt (text and measured width); geometry, data marks and strings unchanged",
                   "x_left_frozen": {"x_left": X_LEFT, "reference": R40_DRAWN, "note": "the longest-label rule evaluated at the old label size, asserted equal to the round-40 record"},
                   "keep_labels": KEEP_LABELS, "kept_at_current_size": KEEP, "nudged": NUDGE},
       "drawn": DRAWN, "records": RECORDS, "written": C.now()}
json.dump(out, open(f"{WORK}/{SHEET}_drawn.json", "w"), indent=1, default=float)
print(f"wrote {WORK}/page_abc.pdf ({W:.2f} x {H_NEW:.3f} pt, DELTA {DELTA} pt, panel c shifted down by C_SHIFT {C_SHIFT} pt) and {WORK}/{SHEET}_drawn.json ({len(DRAWN)} drawn, {len(RECORDS)} value records)")
print(f"round 49: +{UP} pt; kept at the current size: {len(KEEP)} elements ({sorted(set(k['size'] for k in KEEP))} pt); nudged origins: {len(NUDGE)} ({sorted(set(n['text'] for n in NUDGE))})")
