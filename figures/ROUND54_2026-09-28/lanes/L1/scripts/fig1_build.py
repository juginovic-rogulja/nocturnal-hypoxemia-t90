#!/usr/bin/env python3
"""Round 54 lane L1 (2026-09-29, Dragana: "way too small text"): every text +3 pt on the LFIG1E panels (ticks, row labels, keys, callouts 14 pt,
axis titles 15, log-rank lines 13, panel d values 14 and its column header 15, tile titles 14 as labels) and the panels re-laid so that every placement
rule holds with at least 2 pt of slack: panel b callouts and key stacked with an 8-pt gap budget, panel c labels at x 486 with the axes from 660, panel d
scaled from its label and value columns, panel e tiles at a 121-pt pitch. Panel extents stay inside 259.736 to 948.868 (the composed sheet grows only with
the bands). Base: LFIG1E. Round 47 (ROUND47 LFIG1C, 2026-09-24, Dragana): panels b and c SWAPPED (the rank-ordered bar chart of all 141 measurements is now panel b on the left, the family dot plot panel c on the right) and the five callouts at 11 pt (were 9.5). Round 44c copy (ROUND44 LFIG1B): six families (the PLM index drawn inside sleep timing and structure), otherwise the round-37 builder. Original: 01_build for Main_Fig1 panels b, c, d, e (lane V14_L1_RANK, round 37: v8.1 numbers on the V13 design; repointed from the V7_L1_RANK builder).
v8.1: 141 measurements in seven families (the PLM index is a family of its own, drawn in the palette grey #8a9099), the x axes of b and c
run 1 to 141, the callouts of c are placed by a data-driven rule (text clear of every bar and of each other), panel d has 52 outcomes
(the ten largest gains excluding obesity hypoventilation), the Kaplan-Meier tiles follow the v8 tick rule of L2_km_compute.py. Re-plots the four data panels from the v7 numbers
files with matplotlib pinned to the round-30 sheet's page points. Geometry comes from the base sheet's text layer
(base_text.json) and its deduplicated vector records (base_drawings.json, probe_geom.py): the round-30 sheet carries shifted
clipped copies of older panels inside its strips, so every constant below was taken from the record that agrees with the text
layer, never from a heuristic. Writes work/panels_bcde.pdf (page sized, transparent), work/panels_bcde_m.pdf (font metrics
aligned to the sheet), work/placed_texts.json (every string placed) and verify/Main_Fig1_<x>_drawn.json."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import os, sys, json
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C

S = "Main_Fig1"; SD = f"{C.LANE}/{S}"; WORK = f"{SD}/work"; VER = f"{SD}/verify"
os.makedirs(WORK, exist_ok=True); os.makedirs(VER, exist_ok=True)
plt = C.setup_matplotlib()
W, H = 968.66, 1152.96
base_text = json.load(open(f"{WORK}/base_text.json"))
assert base_text["page"] == [W, H], base_text["page"]
SP = base_text["spans"]


def span_of(text, near_y=None, near_x=None):
    c = [s for s in SP if s["text"] == text]
    if near_y is not None: c = [s for s in c if abs(s["origin"][1] - near_y) < 8]
    if near_x is not None: c = [s for s in c if abs(s["origin"][0] - near_x) < 8]
    assert len(c) >= 1, (text, near_y, near_x)
    return c[0]

# ------------------------------------------------------------------ numbers, with their sidecars
RK_PATH = f"{C.NUM}/ranking_v3.csv"; RK_SC = C.sidecar_ok(RK_PATH, need_today=True)
RK = pd.read_csv(C.hydrated(RK_PATH), comment="#"); N_MEAS = len(RK); assert N_MEAS == 141 and sorted(RK["rank"]) == list(range(1, N_MEAS + 1))
RK_SHA = C.sha256(RK_PATH)
FAM_PATH = f"{C.NUM}/measure_families.csv"; FAM_SHA = C.sha256(FAM_PATH); FAM_SC = C.sidecar_ok(FAM_PATH)
FAM = pd.read_csv(C.hydrated(FAM_PATH)); FAMILY = FAM.set_index("feature")["family"]
assert RK.feature.isin(FAMILY.index).all()
RK["family"] = FAMILY.loc[RK.feature].values; RK["g"] = RK.dC * 1000.0
MERGED = RK.family == "Limb movements"; assert int(MERGED.sum()) == 1 and RK.loc[MERGED, "feature"].item() == "plm_index"
RK.loc[MERGED, "family"] = "Sleep timing and structure"   # round 44c (Alen): the limb-movement measurement joins sleep timing and structure
assert set(RK.family) == set(C.FAMILY_COL) == set(C.FAMILY_ORDER_TOP_DOWN), (set(RK.family) ^ set(C.FAMILY_COL))
N_FAM = len(C.FAMILY_ORDER_TOP_DOWN); assert N_FAM == 6
assert (RK.g.rank(ascending=False, method="first").astype(int) == RK["rank"]).all()   # the file's rank is largest gain first
RKI = RK.set_index("feature")
DC_PATH = f"{C.SV}/dC_per_disease/dC_per_disease.csv"; DC_SC = C.sidecar_ok(DC_PATH); DC_SHA = C.sha256(DC_PATH)
PC_PATH = f"{C.SV}/dC_per_disease/_work/positive_control.json"; PC_SC = C.sidecar_ok(PC_PATH); PC = json.load(open(C.hydrated(PC_PATH)))
KM_PATH = f"{paths.FIGURE_ROOT}/ROUND15_2026-09-04/L2_fig1d/L2_fig1d_values.json"; KM_SC = C.sidecar_ok(KM_PATH); KM_SHA = C.sha256(KM_PATH)
KM = json.load(open(C.hydrated(KM_PATH)))
assert KM["source_data"].endswith("data_frozen_v8_2026-09/t90_final.parquet") and KM["cohort_n"] == 19173
# v8: the self-test block is an old-versus-new comparison against the round-14 sheet, logged and not a gate (L2_km_compute.py, F2 fix)
assert KM["selftest"]["mode"].startswith("old-versus-new comparison") and KM["selftest"]["checks"] == 29, KM["selftest"]

page = C.Page(plt, W, H)
PP = C.Page.PT_PLUS; assert PP == 4.0, PP
SZ_TICK, SZ_TITLE, SZ_NOTE = 10.0, 11.0, 9.0   # + PP = 14, 15, 13
CAP = 0.727   # Arial cap height per point (the round-37 constant)
GAP = 8.0     # round 54: the gap budget between stacked text boxes in panel b (the rule below asks for 6)
DRAWN = {k: dict(sheet=S, panel=k, sources={}, values=[]) for k in "bcde"}
DRAWN["b"]["sources"] = DRAWN["c"]["sources"] = {"ranking_v3.csv": dict(path=RK_PATH, sha256=RK_SHA, sidecar=RK_SC),
                                                 "measure_families.csv": dict(path=FAM_PATH, sha256=FAM_SHA, sidecar=FAM_SC)}
DRAWN["d"]["sources"] = {"dC_per_disease.csv": dict(path=DC_PATH, sha256=DC_SHA, sidecar=DC_SC),
                         "positive_control.json": dict(path=PC_PATH, sha256=C.sha256(PC_PATH), sidecar=PC_SC)}
DRAWN["e"]["sources"] = {"L2_fig1d_values.json": dict(path=KM_PATH, sha256=KM_SHA, sidecar=KM_SC)}


def rec(panel, text, x, baseline, source, key, value, kind="printed", **extra):
    DRAWN[panel]["values"].append(dict(text=text, x=round(x, 2), baseline=round(baseline, 2), source=source, key=key,
                                       value=value, kind=kind, **extra))


def T(panel, x, baseline, s, size, ha="left", color=C.INK, rotation=0, static=True):
    """Place a string and record it as a placed text (static strings are identical to the base sheet)."""
    page.text(x, baseline, s, size, ha=ha, color=color, rotation=rotation, record=dict(panel=panel, static=static))

# ================================================================== panel b: the family rank view
# axes patch [141.959, 271.635, 466.123, 536.594]; ticks at 141.959 (1) ... 466.123 (197); median bars 28.025 pt tall, lw 1.8;
# row centres 297.96 + 42.46 i (six families), dots 3.87 pt, jitter 0.19 row units (LA_build.py, RandomState(0))
B = dict(x0=660.0, x1=944.507, y0=271.635, y1=536.594, pad_rows=0.62, jit=0.19, med_half=0.33, dot_d=3.873, med_lw=1.8)   # round 54: x0 660 (was 624.43) so that the 14-pt family labels clear the dots   # round 47: right column (was 141.959..466.123)
LABEL_X_B = 486.0   # round 54: family labels at 486 (were 499.15), the longest (162.1 pt at 14 pt) ends 12 pt before the axes
TICK_BASE_BC, TITLE_BASE_BC = 552.5, 569.0   # round 54: x tick baselines 552.5 (were 550.25) and axis titles 569 (were 565.14) for the 14 / 15 pt text
pitch_pt_v13 = (B["y1"] - B["y0"]) / (5 + 2 * B["pad_rows"]); assert abs(pitch_pt_v13 - 42.46) < 0.02, pitch_pt_v13
pitch_pt = (B["y1"] - B["y0"]) / ((N_FAM - 1) + 2 * B["pad_rows"])          # v8.1: seven family rows inside the V13 box (36.60 pt, was 42.46 for six)
LAB_DY = 301.61 - (B["y0"] + B["pad_rows"] * pitch_pt_v13)               # the V13 family label baseline sits 3.65 pt under its row centre
assert abs(LAB_DY - 3.65) < 0.05, LAB_DY
LAB_DY = LAB_DY + (14.0 - 11.0) * CAP / 2   # round 54: the 14-pt label keeps its cap centre on the row (4.74 pt under the row centre)
axb = page.axes(B["x0"], B["y0"], B["x1"], B["y1"])
axb.set_xlim(1, N_MEAS); axb.set_ylim((N_FAM - 1) + 2 * B["pad_rows"], 0)
centre = {f: B["pad_rows"] + i for i, f in enumerate(C.FAMILY_ORDER_TOP_DOWN)}
rng = np.random.RandomState(0)
dots = []
for f in reversed(C.FAMILY_ORDER_TOP_DOWN):                       # bottom family first, the LA_build.py ORDER
    v = sorted(int(r) for r in RK.loc[RK.family == f, "rank"])
    jit = rng.uniform(-B["jit"], B["jit"], len(v))
    axb.scatter(np.array(v, float), centre[f] - jit, s=B["dot_d"] ** 2, color=C.FAMILY_COL[f], linewidths=0, zorder=3, clip_on=False)
    dots += [dict(rank=rk_, family=f, jit=float(j)) for rk_, j in zip(v, jit)]
assert len(dots) == N_MEAS
MED = {}
for f in C.FAMILY_ORDER_TOP_DOWN:
    m = float(np.median(RK.loc[RK.family == f, "rank"])); MED[f] = m
    axb.plot([m, m], [centre[f] - B["med_half"], centre[f] + B["med_half"]], color=C.INK, lw=B["med_lw"], zorder=4, solid_capstyle="projecting")
    rec("b", f"median rank bar {f}", B["x0"] + (m - 1) * (B["x1"] - B["x0"]) / (N_MEAS - 1.0), B["y0"] + centre[f] * pitch_pt,
        "ranking_v3.csv x measure_families.csv", f"median of rank over family == {f}", m, kind="drawn", n_measurements=int((RK.family == f).sum()))
for s_ in ("top", "right", "left"): axb.spines[s_].set_visible(False)
XT = [1, 50, 100, N_MEAS]          # v8.1: 1, 50, 100, 141 (the 150 tick of the 197 axis has no place on a 141 axis)
axb.set_xticks(XT); axb.set_yticks([])
axb.tick_params(axis="x", length=3.0, width=0.8, labelbottom=False); axb.tick_params(axis="y", length=0)
for v in XT:
    xt = B["x0"] + (v - 1) * (B["x1"] - B["x0"]) / (N_MEAS - 1.0)
    T("b", xt, TICK_BASE_BC, str(v), SZ_TICK, ha="center", static=(v != N_MEAS)); rec("b", str(v), xt, TICK_BASE_BC, "axis", "x tick", v, kind="tick")
TITLE_X_B = (LABEL_X_B + B["x1"]) / 2   # round 54: the 15-pt title (421.6 pt) cannot centre on the 284.5-pt axes inside the page, so it centres on the panel span (labels to axes end)
assert TITLE_X_B + page.width(f"Rank among the {N_MEAS} parameters (rank 1 holds the largest gain)", SZ_TITLE) / 2 < W - 8
T("b", TITLE_X_B, TITLE_BASE_BC, f"Rank among the {N_MEAS} parameters (rank 1 holds the largest gain)", SZ_TITLE, ha="center", static=False)
rec("b", f"Rank among the {N_MEAS} parameters (rank 1 holds the largest gain)", TITLE_X_B, TITLE_BASE_BC, "ranking_v3.csv", "count of rows (the measurement list, N_MEASURES)", N_MEAS, kind="label")
for f in C.FAMILY_ORDER_TOP_DOWN:
    T("b", LABEL_X_B, B["y0"] + centre[f] * pitch_pt + LAB_DY, f, SZ_TICK, static=False)
    rec("b", f, LABEL_X_B, B["y0"] + centre[f] * pitch_pt + LAB_DY, "measure_families.csv", f"family row {f} (row centre + {LAB_DY:.2f} pt)", f, kind="label")
assert LABEL_X_B + max(page.width(f, SZ_TICK) for f in C.FAMILY_ORDER_TOP_DOWN) <= B["x0"] - 8.0, "a family label reaches the dot axes"
DRAWN["b"]["dots"] = dots; DRAWN["b"]["geometry"] = dict(B, label_x=LABEL_X_B, title_x=TITLE_X_B, tick_base=TICK_BASE_BC, title_base=TITLE_BASE_BC, pitch_pt=pitch_pt, pitch_pt_v13=pitch_pt_v13, medians=MED, n_measures=N_MEAS, n_families=N_FAM, label_dy=LAB_DY, xticks=XT, family_n={f: int((RK.family == f).sum()) for f in C.FAMILY_ORDER_TOP_DOWN})

# ================================================================== panel c: the 197 gain bars
# axes patch [546.182, 271.634, 944.507, 536.594]; x ticks 547.043 (1) ... 943.647 (197): 2.0235 pt per rank, xlim (0.575, 197.425);
# y ticks 293.176 (20) ... 508.59 (0), 530.132 (-2): 10.7707 pt per unit, ylim (-2.6, 22.0)
Cg = dict(x0=66.0, x1=466.123, y0=271.634, y1=536.594, xlim=(1 - 0.425, N_MEAS + 0.425), ylim=(-2.6, 22.0))   # round 54: x0 66 (was 63.8) for the 14-pt y tick labels   # round 47: left column (was 546.182..944.507)
DX_C = Cg["x0"] - 546.182   # every absolute x of the old bar chart moves by this      # the same 0.85 bar width padding as the 197 axis
axc = page.axes(Cg["x0"], Cg["y0"], Cg["x1"], Cg["y1"]); axc.set_xlim(*Cg["xlim"]); axc.set_ylim(*Cg["ylim"])
def cx(r): return Cg["x0"] + (r - Cg["xlim"][0]) * (Cg["x1"] - Cg["x0"]) / (Cg["xlim"][1] - Cg["xlim"][0])
def cy(g): return Cg["y1"] - (g - Cg["ylim"][0]) * (Cg["y1"] - Cg["y0"]) / (Cg["ylim"][1] - Cg["ylim"][0])
assert abs(cx(1) - (Cg["x0"] + 0.425 * (Cg["x1"] - Cg["x0"]) / (N_MEAS - 1 + 0.85))) < 0.01 and abs(cy(0) - 508.59) < 0.01 and abs(cy(20) - 293.176) < 0.01   # the y scale is the V13 one
ORDER = RK.sort_values("rank"); A_RANK = ORDER["rank"].to_numpy(float); A_GAIN = ORDER["g"].to_numpy(float)
assert A_GAIN.max() < Cg["ylim"][1] and A_GAIN.min() > Cg["ylim"][0] + 0.1, (A_GAIN.max(), A_GAIN.min())   # v8.1: the tail reaches -2.43, inside the axes box (the -2 tick no longer bounds it, reported)
axc.bar(A_RANK, A_GAIN, width=0.85, linewidth=0, zorder=2, color=[C.FAMILY_COL[f] for f in ORDER.family])
axc.axhline(0.0, color=C.INK, lw=0.8, zorder=3)
for s_ in ("top", "right"): axc.spines[s_].set_visible(False)
axc.set_xticks(XT); axc.set_yticks([-2, 0, 5, 10, 15, 20])
axc.tick_params(length=3.0, width=0.8, labelbottom=False, labelleft=False)
for v in XT:
    T("c", cx(v), TICK_BASE_BC, str(v), SZ_TICK, ha="center", static=(v != N_MEAS)); rec("c", str(v), cx(v), TICK_BASE_BC, "axis", "x tick", v, kind="tick")
YTICK_X_C = Cg["x0"] - 6.5   # right edge of the y tick labels, 3.5 pt left of the 3-pt tick marks
for v in (-2, 0, 5, 10, 15, 20):
    s_ = str(v).replace("-", "−"); bl = cy(v) + (14.0 * CAP) / 2 - 0.4   # round 54: digits centred on the tick (0.4 pt above it, the V13 offset), derived from cy, not from the V13 text layer
    T("c", YTICK_X_C, bl, s_, SZ_TICK, ha="right"); rec("c", s_, YTICK_X_C, bl, "axis", "y tick", v, kind="tick")
T("c", (Cg["x0"] + Cg["x1"]) / 2, TITLE_BASE_BC, f"Rank among the {N_MEAS} parameters", SZ_TITLE, ha="center", static=False)
rec("c", f"Rank among the {N_MEAS} parameters", (Cg["x0"] + Cg["x1"]) / 2, TITLE_BASE_BC, "ranking_v3.csv", "count of rows (the measurement list, N_MEASURES)", N_MEAS, kind="label")
YT_LINE1_X, YT_LINE2_X = 19.5, 37.0   # round 54: the two rotated 15-pt title lines (were 24.72 and 37.58 at 12 pt), 17.5 pt apart, centred on the axes
_yc = (Cg["y0"] + Cg["y1"]) / 2
T("c", YT_LINE1_X, _yc + page.width("Gain in held-out concordance", SZ_TITLE) / 2, "Gain in held-out concordance", SZ_TITLE, rotation=90)
T("c", YT_LINE2_X, _yc + page.width("over age and sex, x1000", SZ_TITLE) / 2, "over age and sex, x1000", SZ_TITLE, rotation=90)
assert YT_LINE2_X + 0.21 * (SZ_TITLE + PP) + 2.0 <= YTICK_X_C - page.width("−2", SZ_TICK), "the rotated title descenders reach the y tick labels"
# the six-entry family key: 5.1 pt squares centred at x 822.05, y 277.85 + 16.156 i, text at x 831.16
KEY_BASE0_V30 = span_of("Oxygenation", near_x=831.16)["origin"][1]; KEY_PITCH = 16.156
# round 54 (Dragana, 2026-09-28: the key sat at the right edge of panel b beside panel c's row labels, "oxygen not aligned with other oxygen, but close"):
# the key moves into the empty upper-middle of panel b, its squares at x 205, texts at x 214.11 (the V30 square-to-text offset 9.11), first baseline 318,
# below the T90 callout box and above the sleep-efficiency and total-sleep-time callout boxes; checked below against every callout box and the bars.
# round 54: 14-pt key at an 18-pt pitch (was 16.156), squares 6.5 pt (were 5.1, scaled with the text, declared), text 10.5 pt right of the square centre (was 9.11),
# first baseline set below from the T90 callout box plus the 8-pt gap budget
KEY_PITCH = 18.0; KEY_SQ = 6.5; KEY_SQ_X = 205.0; KEY_TX = KEY_SQ_X + 10.5; KEY_DY_SQ = -(14.0 * CAP) / 2   # square centre on the cap centre of the 14-pt text
CAP_TICK = (SZ_TICK + PP) * CAP
CALLOUT_PT = SZ_TICK; CALLOUT_CAP = round((CALLOUT_PT + PP) * CAP, 2)   # round 54: callouts 14 pt like every label (round 49 held them at 11 because the AHI and total-sleep-time callouts overlapped at 12: re-placed below)
# the five named bars: rank from the file, layout in data units (the figure1B_polish idiom). Text slots keep the base sheet's
# baselines; the two right-hand slots swap because sleep efficiency (195) and total sleep time (194) are now adjacent bars.
CALL = {"spo2_pct_below_90": "T90 (time below 90% saturation)", "AHI": "Apnea-hypopnea index", "N3_pct": "Deep sleep (N3)",
        "sleep_efficiency_pct": "Sleep efficiency", "TST_min": "Total sleep time"}
R_ = {f: int(RKI.loc[f, "rank"]) for f in CALL}; G_ = {f: float(RKI.loc[f, "g"]) for f in CALL}
assert R_["spo2_pct_below_90"] == 1 and R_["N3_pct"] < R_["TST_min"] < R_["sleep_efficiency_pct"], R_
gT = G_["spo2_pct_below_90"]
def xdata(xp): return Cg["xlim"][0] + (xp - Cg["x0"]) * (Cg["xlim"][1] - Cg["xlim"][0]) / (Cg["x1"] - Cg["x0"])
def ydata(yp): return Cg["ylim"][0] + (Cg["y1"] - yp) * (Cg["ylim"][1] - Cg["ylim"][0]) / (Cg["y1"] - Cg["y0"])
# the V13 slots: T90 beside its bar, AHI in the open area left of the tail, the three sleep measures in a staircase of right-aligned labels at the
# V13 baselines (484.79, 457.86, 430.95) with vertical leaders; every label is then proved clear of every bar and of the other labels (rule below)
# round 54 stack (top down, an 8-pt gap between every pair of boxes, each box = text + 1.5 pt margins): the T90 callout centred on its dot (as V30), the family
# key (six rows at 18 pt), the AHI callout (left-aligned at the V30 x, above the staircase because at 14 pt it no longer fits beside total sleep time),
# then the staircase of right-aligned sleep callouts at their bars: sleep efficiency, total sleep time, deep sleep (V30 order). Every base is derived here.
BOX_H = CALLOUT_CAP + 3.0
T90_BASE = cy(gT) + 2.57
KEY_BASE0 = (T90_BASE + 1.5) + GAP + CAP_TICK + 1.5                       # key top box edge = T90 box bottom + GAP
KEY_BOTTOM = KEY_BASE0 + KEY_PITCH * (N_FAM - 1) + 2.5
AHI_BASE = KEY_BOTTOM + GAP + CALLOUT_CAP + 1.5
SE_BASE = AHI_BASE + BOX_H + GAP; TST_BASE = SE_BASE + BOX_H + GAP; N3_BASE = TST_BASE + BOX_H + GAP
AHI_TX = 614.22 + DX_C
LAY = {"spo2_pct_below_90": dict(tx=cx(7.4), base=T90_BASE, ha="left", leader=[(R_["spo2_pct_below_90"] + 1.4, gT), (7.4 - 1.2, gT)]),
       "AHI": dict(tx=AHI_TX, base=AHI_BASE, ha="left", leader=[(R_["AHI"] + 1.2, G_["AHI"] + 0.55), (xdata(AHI_TX - 1.7), ydata(AHI_BASE - CALLOUT_CAP / 2))]),
       "N3_pct": dict(tx=cx(R_["N3_pct"]), base=N3_BASE, ha="right", leader=[(R_["N3_pct"], G_["N3_pct"] + 0.55), (R_["N3_pct"], ydata(N3_BASE) - 0.49)]),
       "TST_min": dict(tx=cx(R_["TST_min"]), base=TST_BASE, ha="right", leader=[(R_["TST_min"], G_["TST_min"] + 0.55), (R_["TST_min"], ydata(TST_BASE) - 0.49)]),
       "sleep_efficiency_pct": dict(tx=cx(R_["sleep_efficiency_pct"]), base=SE_BASE, ha="right",
                                    leader=[(R_["sleep_efficiency_pct"], G_["sleep_efficiency_pct"] + 0.55), (R_["sleep_efficiency_pct"], ydata(SE_BASE) - 0.49)])}
KEY_MOVE = dict(dx=round(KEY_TX - (831.16 + (Cg["x1"] - 944.507)), 2), dy=round(KEY_BASE0 - KEY_BASE0_V30, 2))
for i, f in enumerate(C.FAMILY_ORDER_TOP_DOWN):
    yc = KEY_BASE0 + KEY_DY_SQ + KEY_PITCH * i
    kax = page.axes(KEY_SQ_X - 7, yc - 5, KEY_SQ_X + 7, yc + 5); kax.set_xlim(KEY_SQ_X - 7, KEY_SQ_X + 7); kax.set_ylim(yc + 5, yc - 5); kax.axis("off")
    kax.add_patch(plt.Rectangle((KEY_SQ_X - KEY_SQ / 2, yc - KEY_SQ / 2), KEY_SQ, KEY_SQ, facecolor=C.FAMILY_COL[f], edgecolor="none"))
    T("c", KEY_TX, KEY_BASE0 + KEY_PITCH * i, f, SZ_TICK, static=False); rec("c", f, KEY_TX, KEY_BASE0 + KEY_PITCH * i, "measure_families.csv", f"key entry {i + 1} ({f})", f, kind="label", moved_r54=KEY_MOVE)
# callout-placement rule: a label box (text width at 9.5 pt, cap height 6.9 pt, 1.5 pt margin) must not touch any bar, any other label or another
# callout's leader; a leader must not cross another label. Proved here, recorded in the drawn file, and the verifier repeats it on the text layer.
BAR_HALF = 0.85 / 2
def label_box(f, L):
    w = page.width(f"{CALL[f]}, rank {R_[f]}", CALLOUT_PT); x0 = L["tx"] if L["ha"] == "left" else L["tx"] - w
    return (x0 - 1.5, L["base"] - CALLOUT_CAP - 1.5, x0 + w + 1.5, L["base"] + 1.5)
def bar_box(r, g):
    return (cx(r - BAR_HALF), cy(max(g, 0.0)), cx(r + BAR_HALF), cy(min(g, 0.0)))
def overlap(a, b): return not (a[2] <= b[0] or b[2] <= a[0] or a[3] <= b[1] or b[3] <= a[1])
BOXES = {f: label_box(f, L) for f, L in LAY.items()}
def _sep(a, b): return max(b[0] - a[2], a[0] - b[2], b[1] - a[3], a[1] - b[3])
SLACK = {}
for f, bx in BOXES.items():
    hits = [int(r) for r, g in zip(A_RANK, A_GAIN) if overlap(bx, bar_box(r, g))]
    assert not hits, (f"callout {CALL[f]} touches the bars of ranks {hits[:6]}: move its slot")
    SLACK[f"{f} to bars"] = round(min(_sep(bx, bar_box(r, g)) for r, g in zip(A_RANK, A_GAIN)), 2)
    for f2, bx2 in BOXES.items():
        if f2 != f: assert not overlap(bx, bx2), (f"callouts {CALL[f]} and {CALL[f2]} overlap")
    SLACK[f"{f} to other callouts"] = round(min(_sep(bx, bx2) for f2, bx2 in BOXES.items() if f2 != f), 2)
    assert bx[2] <= Cg["x1"] + 0.5 and bx[0] >= Cg["x0"] - 0.5, (f, bx)   # inside the axes box
for f, L in LAY.items():
    (xa, ya), (xb_, yb) = L["leader"]; lb = (min(cx(xa), cx(xb_)) - 0.5, min(cy(ya), cy(yb)) - 0.5, max(cx(xa), cx(xb_)) + 0.5, max(cy(ya), cy(yb)) + 0.5)
    for f2, bx2 in BOXES.items():
        if f2 != f: assert not overlap(lb, bx2), (f"leader of {CALL[f]} crosses the label of {CALL[f2]}")
CALLOUT_BOXES = {CALL[f]: [round(v, 2) for v in b] for f, b in BOXES.items()}
KEY_BOXES = []
for i, f in enumerate(C.FAMILY_ORDER_TOP_DOWN):
    w = page.width(f, SZ_TICK); b0 = KEY_BASE0 + KEY_PITCH * i; KEY_BOXES.append([KEY_SQ_X - KEY_SQ / 2, b0 - CAP_TICK - 1.5, KEY_TX + w + 1.5, b0 + 2.5])
_gaps = {name: round(min(_sep(k, bx) for k in KEY_BOXES), 2) for name, bx in BOXES.items()}
_key_right = max(k[2] for k in KEY_BOXES); _key_bottom = max(k[3] for k in KEY_BOXES); _key_left = min(k[0] for k in KEY_BOXES)
_bar_top_under_key = min(cy(g) for r_, g in zip(A_RANK, A_GAIN) if _key_left - 2 <= cx(r_) <= _key_right + 2)   # the smallest y (highest bar top) under the key's x range
KEY_CHECK = dict(gaps_to_callouts=_gaps, key_right=round(_key_right, 2), axes_right=Cg["x1"], key_bottom=round(_key_bottom, 2), highest_bar_top_under_key=round(_bar_top_under_key, 2),
                 ok=bool(min(_gaps.values()) >= 6.0 and _key_right <= Cg["x1"] - 25.0 and _bar_top_under_key >= _key_bottom + 6.0 and min(k[1] for k in KEY_BOXES) >= Cg["y0"] + 4.0))
assert KEY_CHECK["ok"], KEY_CHECK
SLACK["key to callouts (rule 6)"] = round(min(_gaps.values()) - 6.0, 2); SLACK["key to axes right (rule 25)"] = round(Cg["x1"] - 25.0 - _key_right, 2)
SLACK["key to bars beneath (rule 6)"] = round(_bar_top_under_key - _key_bottom - 6.0, 2); SLACK["key below axes top (rule 4)"] = round(min(k[1] for k in KEY_BOXES) - Cg["y0"] - 4.0, 2)
assert min(SLACK.values()) >= 2.0, SLACK   # round 54: every rule holds with at least 2 pt to spare
print("round-54 key placement:", KEY_CHECK); print("round-54 panel b slack (pt):", SLACK)
for f, L in LAY.items():
    (xa, ya), (xb_, yb) = L["leader"]
    axc.plot([xa, xb_], [ya, yb], color=C.GREY_PALE, lw=0.8, zorder=4, solid_capstyle="round")
    axc.scatter([R_[f]], [G_[f]], s=5.48 ** 2, color=C.INK, zorder=5, edgecolors="white", linewidths=0.7, clip_on=False)
    s_ = f"{CALL[f]}, rank {R_[f]}"
    T("c", L["tx"], L["base"], s_, CALLOUT_PT, ha=L["ha"], static=False)
    rec("c", s_, L["tx"], L["base"], "ranking_v3.csv", f"rank where feature == {f}", R_[f], gain_x1000=round(G_[f], 4),
        dot_page_xy=[round(cx(R_[f]), 2), round(cy(G_[f]), 2)], ha=L["ha"])
DRAWN["c"]["bars"] = [dict(rank=int(r), feature=fe, family=fa, gain_x1000=round(float(g), 5), colour=C.FAMILY_COL[fa])
                      for r, fe, fa, g in zip(ORDER["rank"], ORDER.feature, ORDER.family, ORDER.g)]
DRAWN["c"]["geometry"] = dict(Cg, xticks=XT, n_measures=N_MEAS, callout_boxes=CALLOUT_BOXES, callout_ranks=R_, callout_gains={f: round(g, 4) for f, g in G_.items()}, key_pitch=KEY_PITCH, key_base0=KEY_BASE0, key_base0_v30=KEY_BASE0_V30, key_move_r54=KEY_MOVE, key_check_r54=KEY_CHECK, key_square=KEY_SQ, stack_r54=dict(t90=T90_BASE, key_bottom=KEY_BOTTOM, ahi=AHI_BASE, se=SE_BASE, tst=TST_BASE, n3=N3_BASE, gap=GAP), slack_r54=SLACK, ytick_x=YTICK_X_C, ytitle_x=[YT_LINE1_X, YT_LINE2_X], tick_base=TICK_BASE_BC, title_base=TITLE_BASE_BC); DRAWN["c"]["negative_gains"] = int((A_GAIN < 0).sum()); DRAWN["c"]["min_gain_x1000"] = float(A_GAIN.min())

# ================================================================== panel d: per-disease gain forest, top 10 excluding OHS
RAW = pd.read_csv(C.hydrated(DC_PATH)); DIS = RAW[RAW.row_type == "disease"]
STD = DIS[DIS.arm == "standard"].sort_values("dC", ascending=False).reset_index(drop=True); LAG = DIS[DIS.arm == "lag2"].set_index("outcome")
sys.path.insert(0, paths.ANALYSIS_DIR)
from cohort_spec import N_RANKED_OUTCOMES, COHORT_N_FULL
N_OUT = len(STD); assert N_OUT == len(LAG) == N_RANKED_OUTCOMES == 52, (N_OUT, len(LAG), N_RANKED_OUTCOMES)
MEANROW = RAW[(RAW.row_type == "mean_of_outcomes") & (RAW.outcome == "__mean_all__") & (RAW.arm == "standard")].iloc[0]
PUB_DC = float(MEANROW["dC_pubbasis"]); assert int(MEANROW["n_outcomes_in_mean"]) == N_OUT
g1 = PC["gate1_published_machinery"]; assert PC["pass"] is True and g1["pass"] is True and PC["cohort_n"] == COHORT_N_FULL == 15551 and PC["n_outcomes"] == N_OUT
t90_dc = float(RKI.loc["spo2_pct_below_90", "dC"])
assert abs(g1["published_dC"] - t90_dc) < 1e-9 and abs(g1["reproduced_dC"] - t90_dc) < 1e-9, (g1, t90_dc)
assert abs(PUB_DC - t90_dc) < 5e-7, (PUB_DC, t90_dc)          # the dotted rule = the ranking's T90 gain, the mean over 48
ROWS = STD[STD.disease != "Obesity hypoventilation"].head(10).disease.tolist(); assert len(ROWS) == 10
SIX = STD.set_index("disease")
REPLICATES = {}
for dz in ROWS:
    r = SIX.loc[dz]; assert bool(r.estimable) and pd.isna(r.negative_control) and 400 <= int(r.n_replicates) <= 500, (dz, r.n_replicates)   # v8.1: polycythemia carries 438 of 500 bootstrap replicates (recorded, reported)
    assert float(r.ci_lo) < float(r.dC) < float(r.ci_hi) and bool(r.excludes_zero) == (float(r.ci_lo) > 0)
    lg = LAG.loc[str(r.outcome)]
    if bool(lg.estimable):
        assert 400 <= int(lg.n_replicates) <= 500 and float(lg.ci_lo) < float(lg.dC) < float(lg.ci_hi) and bool(lg.excludes_zero) == (float(lg.ci_lo) > 0), (dz, lg.n_replicates)
    else:
        # v8.1: polycythemia (173 incident cases) has NO 2-year landmark estimate (0 replicates): the row keeps its full-follow-up arm and the
        # landmark slot stays empty (no marker, no interval); listed in the report for the owner and the legend lane
        assert int(lg.n_replicates) == 0 and pd.isna(lg.dC), (dz, lg.n_replicates, lg.dC)
    REPLICATES[dz] = dict(standard=int(r.n_replicates), lag2=int(lg.n_replicates), lag2_estimable=bool(lg.estimable))
MISSING_LANDMARK = [dz for dz, v in REPLICATES.items() if not v["lag2_estimable"]]
# geometry from the sheet (probe_geom.py, records that agree with the text layer): axes patch [140.803, 616.776, 360.387, 898.776],
# bottom spine 153.471 (0) to 311.826 (0.15), ticks 0/0.05/0.10/0.15 at 52.785 pt per 0.05, the zero tick 9 pt long, the others 3 pt,
# full-arm CI lines at y 621.787 + 28.2353 i, landmark-arm lines at 633.646 + 28.2353 i, row labels at baseline 630.3 + 28.2353 i
# round 54: the 14-pt row labels (up to 142 pt) and the 14-pt value column (146 pt) no longer fit beside the V13 axes patch, so the axes are laid out from
# the two columns: labels at x 12.25 (were 15.25), the leftmost interval end at least 6 pt right of the longest label, the value column right-aligned at 469
# (was 474.05) and at least 3 pt right of the widest interval, 665 pt per unit of gain (V13 1055.7, v8.1 908.6), ticks every 0.05 as before (0 to 0.20)
LABEL_X_D = 12.25; VAL_RIGHT = 469.0; SCALE = 665.0; LABEL_GAP = 6.0
HI_ALL = float(np.nanmax([max(float(SIX.loc[dz].ci_hi), float(LAG.loc[str(SIX.loc[dz].outcome)].ci_hi) if bool(LAG.loc[str(SIX.loc[dz].outcome)].estimable) else -1.0) for dz in ROWS]))
LO_ALL = float(np.nanmin([min(float(SIX.loc[dz].ci_lo), float(LAG.loc[str(SIX.loc[dz].outcome)].ci_lo) if bool(LAG.loc[str(SIX.loc[dz].outcome)].estimable) else 9.0) for dz in ROWS]))
MAX_LABEL_W = max(page.width(dz, SZ_TICK) for dz in ROWS)
ZERO_X = 169.3
D = dict(x0=ZERO_X + LO_ALL * SCALE - 4.0, x1=VAL_RIGHT - 145.0, y0=616.776, y1=898.776, zero_x=ZERO_X, x015=ZERO_X + 0.15 * SCALE, full_y0=621.787, lag_y0=633.646, pitch=28.2353,
         label_base0=630.3, val_right=VAL_RIGHT, tick_base=917.0, title_base=933.0, header_base=613.0, key_base=609.0, key_y=605.5)   # tick baseline 917 (was 912.43): the 14-pt digits sit under the 9-pt zero tick
SCALE_V13 = 1055.7
XLIM = ((D["x0"] - D["zero_x"]) / SCALE, (D["x1"] - D["zero_x"]) / SCALE)
XTICKS_D = [round(0.05 * k, 2) for k in range(int(np.floor(XLIM[1] / 0.05)) + 1)]; assert XTICKS_D == [0.0, 0.05, 0.1, 0.15, 0.2], XTICKS_D   # the V30 tick strings
AXIS_D = dict(scale_v13=SCALE_V13, scale=SCALE, rescaled=SCALE < SCALE_V13 - 1e-9, xticks=XTICKS_D, max_upper_bound=HI_ALL, min_lower_bound=LO_ALL)
row_c = [(D["full_y0"] + D["lag_y0"]) / 2 + D["pitch"] * i for i in range(10)]; arm_pt = (D["lag_y0"] - D["full_y0"]) / 2
ppu = D["pitch"]; vmax = 10 + (row_c[0] - D["y0"]) / ppu; vmin = vmax - (D["y1"] - D["y0"]) / ppu
axd = page.axes(D["x0"], D["y0"], D["x1"], D["y1"]); axd.set_xlim(*XLIM); axd.set_ylim(vmin, vmax)
def dx(v): return D["x0"] + (v - XLIM[0]) * (D["x1"] - D["x0"]) / (XLIM[1] - XLIM[0])
def dy(v): return D["y1"] - (v - vmin) * (D["y1"] - D["y0"]) / (vmax - vmin)
assert abs(dx(0) - D["zero_x"]) < 0.01 and abs(dy(10) - row_c[0]) < 0.01 and abs(dy(1) - row_c[9]) < 0.01 and abs(dx(0.15) - D["x015"]) < 0.01
axd.axvline(0.0, color=C.INK, lw=0.9, ls=(0, (4, 3)), zorder=1); axd.axvline(PUB_DC, color=C.BLUE, lw=0.9, ls=(0, (1, 2.2)), zorder=1)
for s_ in ("top", "right", "left"): axd.spines[s_].set_visible(False)
axd.spines["bottom"].set_bounds(0.0, XTICKS_D[-1]); axd.set_xticks(XTICKS_D); axd.set_yticks([])
axd.tick_params(axis="x", length=3.0, width=0.8, labelbottom=False); axd.tick_params(axis="y", length=0)
zax = page.axes(D["zero_x"] - 3, D["y1"], D["zero_x"] + 3, D["y1"] + 9); zax.set_xlim(D["zero_x"] - 3, D["zero_x"] + 3); zax.set_ylim(D["y1"] + 9, D["y1"]); zax.axis("off")
zax.plot([D["zero_x"], D["zero_x"]], [D["y1"], D["y1"] + 9], color=C.INK, lw=0.8, solid_capstyle="projecting", clip_on=False)   # the sheet's 9 pt zero tick
arm_units = arm_pt / ppu
def fmt_val(dC, lo, hi): return f"{dC:+.3f} ({lo:.3f} to {hi:.3f})".replace("-", "−")
for i, dz in enumerate(ROWS):
    r = SIX.loc[dz]; lg = LAG.loc[str(r.outcome)]; v = 10 - i
    for arm, rr, dyy, mk in (("standard", r, +arm_units, "o"), ("lag2", lg, -arm_units, "^")):
        if arm == "lag2" and not bool(rr.estimable):
            rec("d", f"{dz} lag2 marker", 0.0, dy(v + dyy), "dC_per_disease.csv", f"row_type==disease & arm==lag2 & disease=={dz}: NOT ESTIMABLE (n_replicates 0), no landmark marker drawn", dict(estimable=False), kind="drawn")
            continue
        dC, lo, hi, sig = float(rr.dC), float(rr.ci_lo), float(rr.ci_hi), bool(rr.excludes_zero)
        axd.plot([lo, hi], [v + dyy, v + dyy], color=C.BLUE, lw=1.7, solid_capstyle="round", zorder=2, clip_on=False)
        if mk == "o":
            axd.plot([dC], [v + dyy], marker="o", ms=6.0, mfc=C.BLUE if sig else "white", mec=C.BLUE, mew=0.8 if sig else 1.1, ls="none", zorder=3, clip_on=False)
        else:
            axd.plot([dC], [v + dyy], marker="^", ms=7.8, mfc=C.BLUE if sig else "white", mec="white" if sig else C.BLUE, mew=0.8 if sig else 1.1, ls="none", zorder=3, clip_on=False)
        rec("d", f"{dz} {arm} marker", dx(dC), dy(v + dyy), "dC_per_disease.csv", f"row_type==disease & arm=={arm} & disease=={dz}: dC, ci_lo, ci_hi, excludes_zero",
            dict(dC=dC, ci_lo=lo, ci_hi=hi, filled=sig), kind="drawn")
        assert dx(lo) >= LABEL_X_D + MAX_LABEL_W + LABEL_GAP, (dz, arm, lo, dx(lo))   # round 54: every interval end clear of the label column (derived, not a literal)
    printed = fmt_val(float(r.dC), float(r.ci_lo), float(r.ci_hi)); base = D["label_base0"] + D["pitch"] * i
    T("d", D["val_right"], base, printed, SZ_TICK, ha="right", static=False)   # round 54: values 14 pt (cell values), were 10.5
    rec("d", printed, D["val_right"], base, "dC_per_disease.csv", f"row_type==disease & arm==standard & disease=={dz}: dC (ci_lo to ci_hi), 3 decimals",
        dict(dC=float(r.dC), ci_lo=float(r.ci_lo), ci_hi=float(r.ci_hi)), ha="right")
    T("d", LABEL_X_D, base, dz, SZ_TICK, static=False); rec("d", dz, LABEL_X_D, base, "dC_per_disease.csv", f"top 10 by dC excluding Obesity hypoventilation, position {i + 1}", dz, kind="label")
VAL_LEFT = min(D["val_right"] - page.width(fmt_val(float(SIX.loc[dz].dC), float(SIX.loc[dz].ci_lo), float(SIX.loc[dz].ci_hi)), SZ_TICK) for dz in ROWS)
assert HI_ALL < XLIM[1] and dx(HI_ALL) < VAL_LEFT - 3.0, (HI_ALL, XLIM, dx(HI_ALL), VAL_LEFT)     # every interval inside the axes patch and clear of the value column
SLACK_D = {"widest interval to value column (rule 3)": round(VAL_LEFT - 3.0 - dx(HI_ALL), 2), "leftmost interval end to label column (rule 6)": round(dx(LO_ALL) - (LABEL_X_D + MAX_LABEL_W + LABEL_GAP), 2),
           "tick label gap": round(0.05 * SCALE - page.width("0.05", SZ_TICK), 2)}
assert min(SLACK_D.values()) >= 2.0, SLACK_D
GUTTER = dict(max_ci_hi=HI_ALL, x_of_max=round(dx(HI_ALL), 2), min_ci_lo=LO_ALL, x_of_min=round(dx(LO_ALL), 2), axes_right=D["x1"], value_text_left=round(VAL_LEFT, 2), label_right=round(LABEL_X_D + MAX_LABEL_W, 2), axis=AXIS_D, slack_r54=SLACK_D)
T("d", D["val_right"], D["header_base"], "Gain (95% CI)", SZ_TITLE, ha="right")   # round 54: the column header 15 pt (was 10.5), baseline 613 (was 616.71)
for v in XTICKS_D:
    s_ = "0" if v == 0 else f"{v:.2f}"
    T("d", dx(v), D["tick_base"], s_, SZ_TICK, ha="center", static=False); rec("d", s_, dx(v), D["tick_base"], "axis", "x tick", v, kind="tick")
TITLE_X_D = (dx(0.0) + dx(XTICKS_D[-1])) / 2   # round 54: the 15-pt title centred on the tick span (was left-aligned at 99.52)
T("d", TITLE_X_D, D["title_base"], "Gain in held-out concordance from adding T90 to age and sex", SZ_TITLE, ha="center")
assert TITLE_X_D - page.width("Gain in held-out concordance from adding T90 to age and sex", SZ_TITLE) / 2 >= LABEL_X_D
# the key: 16-pt handle lines (lw 1.7) with the markers at their centres, re-spaced for the 14-pt texts (V30: handles 15.253-31.253 and 115.159-131.159 at y 602.444)
KEY_D = []; kx0 = 10.0
for mk, lab in (("o", "Full follow-up"), ("^", "2-year landmark")):
    kx1 = kx0 + 16.0; tx = kx1 + 7.0; ky = D["key_y"]
    kax = page.axes(kx0 - 5, ky - 6, kx1 + 5, ky + 6); kax.set_xlim(kx0 - 5, kx1 + 5); kax.set_ylim(ky + 6, ky - 6); kax.axis("off")
    kax.plot([kx0, kx1], [ky, ky], color=C.BLUE, lw=1.7, solid_capstyle="round", zorder=2)
    kax.plot([(kx0 + kx1) / 2], [ky], marker=mk, ms=6.0 if mk == "o" else 7.8, mfc=C.BLUE, mec="white", mew=0.8, ls="none", zorder=3)
    T("d", tx, D["key_base"], lab, SZ_TICK); KEY_D.append(dict(handle=[kx0, kx1], y=ky, text_x=tx, label=lab)); kx0 = tx + page.width(lab, SZ_TICK) + 16.0
DRAWN["d"]["geometry"] = dict(D, xlim=XLIM, scale_pt_per_unit=SCALE, axis=AXIS_D, xticks=XTICKS_D, rows_centre_pt=row_c, arm_offset_pt=arm_pt, pub_dc_rule=PUB_DC, pub_dc_rule_x=dx(PUB_DC), rows=ROWS, n_outcomes=N_OUT, replicates=REPLICATES, gutter=GUTTER, missing_landmark=MISSING_LANDMARK, label_x=LABEL_X_D, title_x=TITLE_X_D, key=KEY_D)
rec("d", f"dotted rule (mean gain over the {N_OUT} conditions)", dx(PUB_DC), 0, "dC_per_disease.csv", "row_type==mean_of_outcomes & outcome==__mean_all__ & arm==standard: dC_pubbasis", PUB_DC, kind="drawn")

# ================================================================== panel e: eight Kaplan-Meier tiles
# round 54: the 13-pt log-rank line (110.4 pt) and the 14-pt "Respiratory failure" title (113.4 pt) exceed the V30 tile pitch of 112.23 and, left-aligned
# at the spine, the last tile's line would leave the page, so each tile's title and log-rank line are CENTRED on its axis (V30: left-aligned at the spine)
# and the four tiles of a row sit at a 118-pt pitch from spine 514.2 (the panel starts at x 475.9 with the rotated 15-pt axis title, the 14-pt y tick
# labels end 5.5 pt left of each spine), tile axes 83.707 wide and 93.75 tall as before, title 20.5 and log-rank line 5 pt above the tile top (16.25, 4.73)
E = dict(spines=[514.2 + 118.0 * k for k in range(4)], axw=83.707, box_h=93.75, bottom1=746.78, bottom2=898.78, spine_w=0.5946, tick=2.23,
         ylab_gap=5.5, ylab_dy=(14.0 * CAP) / 2 - 1.34, xlab_dy=12.5, title_up=20.5, logrank_up=5.0, curve_lw=1.05, pitch=118.0, title_pt=SZ_TICK + PP, logrank_pt=SZ_NOTE + PP,
         key_y=601.85, key_base=605.35, xtitle_base=933.0, ytitle_x=486.8, left=475.9, title_align="center")
assert E["spines"][-1] + E["axw"] + page.width("8", SZ_TICK) / 2 <= W - 3.0, "the last tile's 8 reaches the page edge"
TILE_TEXT_BOXES = []   # (x0, x1) of every tile title and log-rank line, for the neighbour gap rule
BAND_COL = {"<=1%": "#b3dcf2", "1-10%": "#5cb3e4", ">10%": "#0288d1"}     # the sheet's T90 ramp (what shows)
ZORD = {"<=1%": 3, "1-10%": 4, ">10%": 5}
tiles = [(p, 1) for p in KM["row1_existing"]] + [(p, 2) for p in KM["row2_new"]]
assert [p["label"] for p, _ in tiles] == ["Cirrhosis", "Heart failure", "Type 2 diabetes", "Hypertension", "Respiratory failure", "Obesity", "COPD", "Pneumonia"]
def clean_yticks(ymax_lim, n_max=6):
    """the v8 rule of ROUND15 L2_km_compute.py: 4 to 6 clean majors from 0, step 3 tried after 5 (axis tops in [12, 15) and [18, 19))"""
    for step in (1, 2, 5, 3, 10, 20):
        top = int(np.floor(ymax_lim / step)) * step; n = top // step + 1
        if 4 <= n <= n_max: return list(range(0, top + 1, step))
    raise AssertionError(ymax_lim)
for k, (pn, row) in enumerate(tiles):
    sx = E["spines"][k % 4]; bottom = E["bottom1"] if row == 1 else E["bottom2"]; top = bottom - E["box_h"]
    assert abs(pn["ylim"] - 1.14 * pn["ymax"]) < 1e-6 and pn["yticks"] == clean_yticks(pn["ylim"]), (pn["label"], pn["yticks"])
    assert pn["p_printed"] == ("P < 0.001" if pn["p"] < 0.001 else None), (pn["label"], pn["p"], pn["p_printed"])
    assert sum(b["n"] for b in pn["bands"]) == pn["n"] and sum(b["events"] for b in pn["bands"]) == pn["events"]
    ax = page.axes(sx, top, sx + E["axw"], bottom); ax.set_xlim(0, 8.0); ax.set_ylim(0, pn["ylim"])
    for b in pn["bands"]:
        t = np.asarray(b["curve_t"]); y = np.asarray(b["curve_cum_incidence_pct"])
        # the stored curve arrays are rounded to six decimals, so the end value and the maximum agree with ci8_pct and ymax to 1e-6
        assert t[-1] == 8.0 and t[0] == 0.0 and y[0] == 0.0 and abs(y[-1] - b["ci8_pct"]) < 1e-6 and y.max() <= pn["ymax"] + 1e-6 and (np.diff(y) >= -1e-9).all(), (pn["label"], b["band"])
        ax.step(t, y, where="post", color=BAND_COL[b["band"]], lw=E["curve_lw"], solid_capstyle="round", solid_joinstyle="round", zorder=ZORD[b["band"]])
        rec("e", f"{pn['label']} curve {b['band']}", sx + E["axw"], bottom - b["ci8_pct"] / pn["ylim"] * E["box_h"], "L2_fig1d_values.json",
            f"{'row1_existing' if row == 1 else 'row2_new'}[label=={pn['label']}].bands[band=={b['band']}]: curve_cum_incidence_pct at 8 years, n, events",
            dict(ci8_pct=b["ci8_pct"], n=b["n"], events=b["events"]), kind="drawn")
    for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
    ax.spines["left"].set_bounds(0, pn["yticks"][-1]); ax.spines["left"].set_linewidth(E["spine_w"]); ax.spines["bottom"].set_linewidth(E["spine_w"])
    ax.set_yticks(pn["yticks"]); ax.set_xticks([0, 2, 4, 6, 8])
    ax.tick_params(length=E["tick"], width=E["spine_w"], labelbottom=False, labelleft=False)
    for v in pn["yticks"]:
        yt = bottom - v / pn["ylim"] * E["box_h"]
        T("e", sx - E["ylab_gap"], yt + E["ylab_dy"], str(v), SZ_TICK, ha="right", static=False)
        rec("e", str(v), sx - E["ylab_gap"], yt + E["ylab_dy"], "L2_fig1d_values.json", f"{pn['label']}.yticks", v, kind="tick", ha="right")
    for v in (0, 2, 4, 6, 8):
        T("e", sx + E["axw"] * v / 8.0, bottom + E["xlab_dy"], str(v), SZ_TICK, ha="center")
    xc_ = sx + E["axw"] / 2
    T("e", xc_, top - E["title_up"], pn["label"], SZ_TICK, ha="center"); rec("e", pn["label"], xc_, top - E["title_up"], "L2_fig1d_values.json", "label", pn["label"], kind="label", ha="center")   # round 54: tile titles 14 pt (label class, declared: at 15 pt "Respiratory failure" needs a 130-pt pitch the page cannot give), centred on the tile
    s_ = "Log-rank " + pn["p_printed"]; T("e", xc_, top - E["logrank_up"], s_, SZ_NOTE, ha="center")   # round 54: 13 pt (was 10.5), centred on the tile
    for s2, sz in ((pn["label"], SZ_TICK), (s_, SZ_NOTE)):
        w2 = page.width(s2, sz); TILE_TEXT_BOXES.append(dict(tile=k, text=s2, x0=round(xc_ - w2 / 2, 2), x1=round(xc_ + w2 / 2, 2)))
        assert xc_ + w2 / 2 <= W - 3.0 and xc_ - w2 / 2 >= E["left"] + 8.0, (s2, xc_ - w2 / 2, xc_ + w2 / 2)   # inside the page, clear of the rotated axis title
    rec("e", s_, sx, top - E["logrank_up"], "L2_fig1d_values.json", f"{pn['label']}.p (chi2 {pn['chi2']:.2f}, 2 df)", dict(p=pn["p"], chi2=pn["chi2"], n=pn["n"], events=pn["events"]))
# key at the top of the panel: 15-pt handle lines (lw 1.9) as on the sheet, re-spaced for the 14-pt texts (V30: handles from 504.95, 587.15, 680.806 at y 597.83, texts 7 pt right)
KEY_E = []; kx0 = E["spines"][0] - 8.8
kax = page.axes(495, 590, 820, 612); kax.set_xlim(495, 820); kax.set_ylim(612, 590); kax.axis("off")
for s_, col in (("T90 ≤1%", "#b3dcf2"), ("T90 1–10%", "#5cb3e4"), ("T90 >10%", "#0288d1")):
    kx1 = kx0 + 15.0; kax.plot([kx0, kx1], [E["key_y"], E["key_y"]], color=col, lw=1.9, solid_capstyle="butt", zorder=2)
    T("e", kx1 + 7.0, E["key_base"], s_, SZ_TICK); KEY_E.append(dict(handle=[kx0, kx1], text_x=kx1 + 7.0, label=s_)); kx0 = kx1 + 7.0 + page.width(s_, SZ_TICK) + 16.0
for k in range(len(tiles) - 1):   # a tile's title and log-rank line end at least 6 pt before the next tile's begin (row neighbours only)
    if k % 4 == 3: continue
    a = [b for b in TILE_TEXT_BOXES if b["tile"] == k]; b_ = [b for b in TILE_TEXT_BOXES if b["tile"] == k + 1]
    gap = min(bb["x0"] for bb in b_) - max(aa["x1"] for aa in a); assert gap >= 6.0, (k, gap)
    E.setdefault("neighbour_gaps", []).append(round(gap, 2))
_xc = (E["spines"][0] + E["spines"][-1] + E["axw"]) / 2
T("e", _xc, E["xtitle_base"], "Years since the sleep study", SZ_TITLE, ha="center")   # round 54: centred on the tiles (V30 671.77 left-aligned)
_yc = (E["bottom1"] - E["box_h"] + E["bottom2"]) / 2
T("e", E["ytitle_x"], _yc + page.width("Cumulative incidence, %", SZ_TITLE) / 2, "Cumulative incidence, %", SZ_TITLE, rotation=90)   # centred on the two rows
assert E["ytitle_x"] + 0.21 * (SZ_TITLE + PP) + 3.0 <= E["spines"][0] - E["ylab_gap"] - max(page.width(str(v), SZ_TICK) for p_, _ in tiles for v in p_["yticks"]), "the rotated title reaches the tick labels"
assert E["left"] - VAL_RIGHT >= 5.0 and abs(E["ytitle_x"] - (SZ_TITLE + PP) * CAP - E["left"]) < 0.2, (E["left"], VAL_RIGHT, E["ytitle_x"])   # the panel starts 5 pt or more right of the panel d value column
DRAWN["e"]["geometry"] = dict(E, band_colours=BAND_COL, band_n=KM["band_n"], tick_sets={p["label"]: p["yticks"] for p, _ in tiles}, key=KEY_E, xtitle_x=_xc, tile_text_boxes=TILE_TEXT_BOXES)

# ================================================================== write
bad = C.banned_words([t["text"] for t in page.texts]); assert not bad, bad
PANELS_TOP, PANELS_BOTTOM = 259.736, 948.868   # the composed sheet clips the panels page to the slot between band a and band f (compose_r44b)
for t in page.texts:
    top = t["baseline"] - t["size"] * CAP if t["rotation"] == 0 else t["baseline"] - t["width"]
    assert PANELS_TOP + 2 <= top and t["baseline"] + 0.21 * t["size"] <= PANELS_BOTTOM - 2, (t["text"], top, t["baseline"])
SIZES = sorted({t["size"] for t in page.texts}); assert SIZES == [13.0, 14.0, 15.0], SIZES
out = f"{WORK}/panels_bcde.pdf"; page.fig.savefig(out, transparent=True); plt.close(page.fig)
met = {}
for s_ in base_text["spans"]:
    if s_.get("ascender") is not None: met.setdefault(s_["font"].split("+")[-1], (int(round(s_["ascender"] * 1000)), int(round(s_["descender"] * 1000))))
print("sheet font metrics (from the V13 text record):", met)
done = C.align_font_metrics(out, f"{WORK}/panels_bcde_m.pdf", met); print("aligned:", done)
for k in "bcde":
    json.dump(DRAWN[k], open(f"{VER}/{S}_{k}_drawn.json", "w"), indent=1, default=float)
json.dump(page.texts, open(f"{WORK}/placed_texts.json", "w"), indent=0, default=float)
print(f"wrote {out} and metrics-aligned twin; {sum(len(DRAWN[k]['values']) for k in 'bcde')} recorded values, {len(page.texts)} strings placed")
print("callout ranks:", R_, "| panel d rows:", ROWS, "| dotted rule", round(PUB_DC, 6), "| replicates", REPLICATES, "| gutter", GUTTER)
