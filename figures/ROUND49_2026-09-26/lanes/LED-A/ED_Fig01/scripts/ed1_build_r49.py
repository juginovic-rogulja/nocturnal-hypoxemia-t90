#!/usr/bin/env python3
"""ROUND 49, lane LED-A (2026-09-26): ED_Fig01 rebuilt at +1 pt. Copy of ROUND37_2026-09-14/figures/lanes/V14_L1_RANK/scripts/ed1_build.py
(the V14 builder, v8.1 numbers on the V13 design) with the later steps of the chain folded in, so that the result equals the V26 sheet
except for text size:
  +1 pt        every text through common.Page.text()/width() (DELTA = 1.0 in the copied helper; the wrap rule measures the larger text)
  round 40     the panel a x-axis title "Share of the cohort, %" set as two lines "Share of the" / "cohort, %" centred on the axis
               (x 79.59 to 147.99), first baseline 293.0, second 11.41 pt below (ed1_edit.py), here placed by the builder itself
  round 43     panel b 14.4 pt further right on a page 14.4 pt wider (widen_ed1.py split the page at x = 160): every element of panel b
               (letter b, names, leaders, dots, gains, families, axis, ticks, axis title) drawn at x + 14.4, the page 533.14 pt wide
  round 40     the 16 pt title strip and "Extended Data Fig. 1" (14 pt) are stamped afterwards by scripts_shared/stamp_title_r49.py
Everything else (numbers, sidecar gates, geometry, colours, data marks) is the round-37 builder's, unchanged.
Writes work/sheet_raw.pdf, work/sheet_built.pdf (font metrics aligned to the V13 base; the stamp step makes ED_Fig01.pdf), work/placed_texts.json,
verify/ED_Fig01_{a,b}_drawn.json, work/build_record.json."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import os, sys, json
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C

S = "ED_Fig01"; SD = f"{C.LANE}/{S}"; WORK = f"{SD}/work"; VER = f"{SD}/verify"
os.makedirs(WORK, exist_ok=True); os.makedirs(VER, exist_ok=True)
plt = C.setup_matplotlib()
base_text = json.load(open(f"{WORK}/base_text.json")); W0, H_V13 = base_text["page"]
assert [W0, H_V13] == [518.74, 717.92], (W0, H_V13)          # the V13 sheet (30 named rows)
H_BASE = 677.81                                                # the 28-row geometry every constant below is pinned to (the round-28 base the R30 builder measured)
assert abs(H_V13 - (H_BASE + 2 * (616.124 - 74.648) / 27.0)) < 0.02, "the V13 page is not the 28-row base grown by two pitches"
SP = base_text["spans"]
LET = {L["text"]: L for L in base_text["letters"]}; assert set(LET) == {"a", "b"}
XB = 14.4                                                      # round 43: panel b 14.4 pt further right (widen_ed1.py, split at x = 160)
SPLIT = 160.0
W = W0 + XB
assert C.DELTA == 1.0

# ------------------------------------------------------------------ numbers, with their sidecars
RK_PATH = f"{C.NUM}/ranking_v3.csv"; RK_SC = C.sidecar_ok(RK_PATH, need_today=True); RK_SHA = C.sha256(RK_PATH)
RK = pd.read_csv(C.hydrated(RK_PATH), comment="#"); N_MEASURES = len(RK); assert sorted(RK["rank"]) == list(range(1, N_MEASURES + 1))
sys.path.insert(0, paths.ANALYSIS_DIR)
from cohort_spec import N_MEASURES as N_MEASURES_SPEC
assert N_MEASURES == N_MEASURES_SPEC == 141, (N_MEASURES, N_MEASURES_SPEC)
FAM_PATH = f"{C.NUM}/measure_families.csv"; FAM_SHA = C.sha256(FAM_PATH); FAM_SC = C.sidecar_ok(FAM_PATH)
FAMILY = pd.read_csv(C.hydrated(FAM_PATH)).set_index("feature")["family"]; assert RK.feature.isin(FAMILY.index).all()
RK["family"] = FAMILY.loc[RK.feature].values; RK["g"] = RK.dC * 1000.0
assert (RK.g.rank(ascending=False, method="first").astype(int) == RK["rank"]).all()
RES_PATH = f"{C.NUM}/results_v2.json"; RES_SC = C.sidecar_ok(RES_PATH); RES_SHA = C.sha256(RES_PATH); RES = json.load(open(C.hydrated(RES_PATH)))
COH_PATH = f"{C.NUM}/cohorts.json"; COH_SC = C.sidecar_ok(COH_PATH); COH_SHA = C.sha256(COH_PATH); COH = json.load(open(C.hydrated(COH_PATH)))
N_COHORT = int(COH["bdsp"]["n_analysis"])
from cohort_spec import apply_cohort, COHORT_N, DATA_DIR
assert N_COHORT == COHORT_N == RES["cohort"]["n"] == 19173 and DATA_DIR.endswith("data_frozen_v8_2026-09"), (N_COHORT, COHORT_N, DATA_DIR)
T90F = f"{DATA_DIR}/t90_final.parquet"; C.hydrated(T90F)
t90 = apply_cohort(pd.read_parquet(T90F, columns=["BDSPPatientID", "fu_valid", "spo2_pct_below_90", "oximetry_bad"]))["spo2_pct_below_90"]
assert len(t90) == N_COHORT and t90.notna().all()

# ---- panel a bands: the file's counts and shares, re-derived from the v7 table and from each other
BANDS = [("≤1", lambda t: t <= 1, "0-1%"), (">1 to 5", lambda t: (t > 1) & (t <= 5), "1-5%"), (">5 to 10", lambda t: (t > 5) & (t <= 10), "5-10%"), (">10", lambda t: t > 10, ">10%")]
BAND_N, BAND_PCT = [], []
for lab, rule, key in BANDS:
    n_file = int(RES["cohort"]["t90_bands_n"][key]); n_re = int(rule(t90).sum()); pct_file = float(RES["cohort"]["t90_bands_pct"][key])
    assert n_re == n_file == int(RES["graded"]["band_n"][key]), (lab, n_re, n_file, RES["graded"]["band_n"][key])
    assert abs(100.0 * n_file / N_COHORT - pct_file) < 0.051, (lab, pct_file, 100.0 * n_file / N_COHORT)
    BAND_N.append(n_file); BAND_PCT.append(100.0 * n_file / N_COHORT)
assert sum(BAND_N) == N_COHORT and int(RES["cohort"]["t90_bands_n"][">10%"]) == int(COH["bdsp"]["t90_above_10_n"]) and max(BAND_PCT) < 60.0
BAND_COL = ["#b3dcf2", "#7cc0e9", "#3f9fd8", "#0288d1"]                      # the base sheet's T90 ramp (sampled from its vector records)

# ---- panel b rows: gain above 0.003, siteZ twins folded (the fold rule stated from the v7 data, every condition asserted)
TAIL_G = 3.0
raw = set(RK.feature); RKI = RK.set_index("feature")
pairs = [(f, f[:-6]) for f in RK.feature if f.endswith("_siteZ") and f[:-6] in raw]
twin_diff = {f: abs(float(RKI.loc[f, "dC"]) - float(RKI.loc[b, "dC"])) for f, b in pairs}
assert len(pairs) == 0, ("v8.1 ranking still carries siteZ twins", pairs[:5])
dup = [f for f, _ in pairs]; distinct = RK[~RK.feature.isin(dup)].copy(); assert len(distinct) == N_MEASURES
K = int((RK.g > TAIL_G).sum()); assert float(RK.sort_values("rank").g.iloc[K]) <= TAIL_G and float(RK.sort_values("rank").g.iloc[K - 1]) > TAIL_G
tail = distinct[distinct.g > TAIL_G].sort_values("g", ascending=False).reset_index(drop=True)
N_ROWS = len(tail); assert N_ROWS == K - sum(1 for f, b in pairs if float(RKI.loc[f, "g"]) > TAIL_G), (N_ROWS, K)
TOP5_FAMILIES = sorted(set(RK[RK["rank"] <= 5].family)); print("families in the top 5:", TOP5_FAMILIES)
NICE = {"spo2_nadir_corrected": "Lowest saturation", "spo2_pct_below_90": "Time below 90% saturation", "spo2_pct_below_88": "Time below 88% saturation",
        "spo2_mean": "Mean saturation", "odi3_total": "3% desaturation index", "odi4_total": "4% desaturation index", "nrem_hr_bpm": "Heart rate in non-REM",
        "wholenight_hr_bpm": "Heart rate, whole night", "wholenight_lf_hf_ratio": "LF to HF ratio, whole night", "nrem_lf_hf_ratio": "LF to HF ratio in non-REM",
        "nrem_lf_power": "LF power in non-REM", "wholenight_lf_power": "LF power, whole night", "wholenight_hf_power": "HF power, whole night",
        "wholenight_n_rr": "Heartbeat intervals counted, whole night", "nrem_rmssd": "RMSSD in non-REM", "wholenight_rmssd": "RMSSD, whole night",
        "F_spindle_n_spindles": "Spindle count, frontal", "F_spindle_density_per_min": "Spindle density, frontal", "F_SO_mean_dur_sec": "Slow-oscillation duration, frontal",
        "C_SO_mean_dur_sec": "Slow-oscillation duration, central", "C_spindle_n_spindles": "Spindle count, central", "C_spindle_density_per_min": "Spindle density, central",
        "C_coupling_strength_r": "Spindle to slow-oscillation coupling, central", "O_coupling_strength_r": "Spindle to slow-oscillation coupling, occipital",
        "O_coupling_overlap_frac": "Spindle to slow-oscillation overlap, occipital", "C_nrem_sigma_mean": "Sigma power in non-REM, central",
        "resp_events_any_label": "Respiratory events, count", "AHI": "Apnea-hypopnea index",
        "O_spindle_n_spindles": "Spindle count, occipital", "C_nrem_theta_kurt": "Theta power kurtosis in non-REM, central",
        "hypoxic_burden": "Hypoxic burden", "plm_index": "Periodic limb movement index", "sleep_onset_latency_min": "Sleep onset latency", "rem_latency_min": "REM latency"}
BASE_NAMES = {s["text"] for s in SP}
NEW_NAMES = None
missing = sorted(f for f in tail.feature if f not in NICE); assert not missing, f"no plain-English name for: {missing}"
print(f"{N_ROWS} rows with gain above 0.003 (V13 printed 30)")

# ================================================================== geometry from the base sheet (V13 coordinates; panel b drawn at x + XB)
PITCH = (616.124 - 74.648) / 27.0                                   # 20.0547 pt, the base's row pitch (28 dot centres)
GROW = (N_ROWS - 28) * PITCH; H = H_BASE + GROW
A = dict(x0=79.59, y0=33.122, x1=147.99, y1=263.522, xmax=60.0, ylim=(-0.72, 3.72), bar_h=0.60, ylab_right=73.08, ylab_dy=3.58, xtick_base=277.18,
         xlabel=(25.68, 293.0), ylabel_lines=[("Time below 90% saturation,", 20.84, 216.06), ("% of the recording", 32.25, 193.01)])
XLAB_LINES = [("Share of the", 293.0), ("cohort, %", 293.0 + (32.25 - 20.84))]   # round 40 (ed1_edit.py): two lines centred on the panel a axis, 11.41 pt pitch
XLAB_XC = (A["x0"] + A["x1"]) / 2.0
B = dict(x0=292.59, y0=57.602, x1=371.05, xlim=(2.85, 24.5), y_row0=74.648, name_x=168.15, one_dy=3.58, two_dy=(-0.99, 8.15), val_right=395.35, val_dy=2.58,
         fam_x=402.55, fam_dy=2.58, fam_head=(402.55, 51.17), xtick_base_old=646.83, title_old=("Gain in held-out concordance, in units of 0.001", 217.96, 661.69),
         name_max_pt=1.67 * 72.0, dot_s=28.0, leader_col="#b3dcf2")
B["y1"] = B["y0"] + (N_ROWS + 0.7) * PITCH; B["ylim"] = (-0.85, N_ROWS - 1 + 0.85)
assert abs((B["y0"] + (28 + 0.7) * PITCH) - 633.17) < 0.02                 # the base's axes bottom for 28 rows
ppu = (B["x1"] - B["x0"]) / (B["xlim"][1] - B["xlim"][0])
def bx(g): return B["x0"] + (g - B["xlim"][0]) * ppu                       # V13 x of a gain (panel b drawn at bx(g) + XB)
def row_y(i): return B["y_row0"] + PITCH * i
assert abs(bx(5) - 300.382) < 0.02 and abs(bx(20) - 354.742) < 0.02 and abs(bx(TAIL_G) - 293.134) < 0.02
A_ppu = (A["x1"] - A["x0"]) / A["xmax"]; A_pitch = (A["y1"] - A["y0"]) / (A["ylim"][1] - A["ylim"][0])
def band_y(k): return A["y0"] + (A["ylim"][1] - (3 - k)) * A_pitch                   # band k (0 = top) centre
assert abs(band_y(0) - 70.484) < 0.02 and abs(band_y(3) - 226.16) < 0.02 and abs(A_ppu * 20 - 22.8) < 0.02
assert all(v >= SPLIT for v in (B["x0"], B["name_x"], B["val_right"], B["fam_x"], B["title_old"][1], LET["b"]["origin"][0])) and all(v < SPLIT for v in (A["x1"], A["ylab_right"], XLAB_XC + 40))

page = C.Page(plt, W, H); page.fig.patch.set_alpha(1.0); page.fig.patch.set_facecolor("white")
DRAWN = {"a": dict(sheet=S, panel="a", sources={"results_v2.json": dict(path=RES_PATH, sha256=RES_SHA, sidecar=RES_SC), "cohorts.json": dict(path=COH_PATH, sha256=COH_SHA, sidecar=COH_SC),
                                                "t90_final.parquet (v8)": dict(path=T90F, sha256=None, use="band counts re-derived through cohort_spec.apply_cohort, asserted equal")}, values=[]),
         "b": dict(sheet=S, panel="b", sources={"ranking_v3.csv": dict(path=RK_PATH, sha256=RK_SHA, sidecar=RK_SC), "measure_families.csv": dict(path=FAM_PATH, sha256=FAM_SHA, sidecar=FAM_SC)}, values=[])}
def rec(panel, text, x, baseline, source, key, value, kind="printed", **extra):
    DRAWN[panel]["values"].append(dict(text=text, x=round(x, 2), baseline=round(baseline, 2), source=source, key=key, value=value, kind=kind, **extra))
def T(panel, x, baseline, s, size, ha="left", color=C.INK, rotation=0, static=True, weight="normal"):
    page.text(x, baseline, s, size, ha=ha, color=color, rotation=rotation, weight=weight, record=dict(panel=panel, static=static))

# ================================================================== panel a
axa = page.axes(A["x0"], A["y0"], A["x1"], A["y1"]); axa.set_xlim(0, A["xmax"]); axa.set_ylim(*A["ylim"])
ys = np.arange(4, dtype=float)[::-1]
axa.barh(ys, BAND_PCT, A["bar_h"], color=BAND_COL, linewidth=0, zorder=2)
for k, ((lab, _, key), n, p) in enumerate(zip(BANDS, BAND_N, BAND_PCT)):
    rec("a", f"bar {lab}", A["x0"] + p * A_ppu, band_y(k), "results_v2.json", f"cohort.t90_bands_pct[{key}], cohort.t90_bands_n[{key}]", dict(pct=p, n=n), kind="drawn",
        bar_rect=[round(A["x0"], 3), round(band_y(k) - A["bar_h"] * A_pitch / 2, 3), round(A["x0"] + p * A_ppu, 3), round(band_y(k) + A["bar_h"] * A_pitch / 2, 3)], colour=BAND_COL[k])
    T("a", A["ylab_right"], band_y(k) + A["ylab_dy"], lab, 10.0, ha="right"); rec("a", lab, A["ylab_right"], band_y(k) + A["ylab_dy"], "axis", "band label", lab, kind="label", ha="right")
for s_ in ("top", "right"): axa.spines[s_].set_visible(False)
axa.set_xticks([0, 20, 40, 60]); axa.set_yticks(list(ys)); axa.tick_params(length=3.0, width=0.8, labelbottom=False, labelleft=False)
for v in (0, 20, 40, 60):
    T("a", A["x0"] + v * A_ppu, A["xtick_base"], str(v), 10.0, ha="center"); rec("a", str(v), A["x0"] + v * A_ppu, A["xtick_base"], "axis", "x tick", v, kind="tick", ha="center")
for s_, y in XLAB_LINES:                                                    # round 40: two centred lines instead of the one left-aligned line
    T("a", XLAB_XC, y, s_, 11.0, ha="center"); rec("a", s_, XLAB_XC, y, "axis", "x-axis title line (round 40, two centred lines)", s_, kind="label", ha="center")
for s_, x, y in A["ylabel_lines"]: T("a", x, y, s_, 11.0, rotation=90)
DRAWN["a"]["geometry"] = dict(A, ppu=A_ppu, pitch=A_pitch, bands=[dict(label=l, key=k, n=n, pct=p, colour=c_) for (l, _, k), n, p, c_ in zip(BANDS, BAND_N, BAND_PCT, BAND_COL)], xlabel_lines=XLAB_LINES, xlabel_xc=XLAB_XC)

# ================================================================== panel b
def wrap2(text, delta=None):
    """The base sheet's wrap rule: at most two lines of 1.67 in, words unchanged, the longest first line that leaves a fitting second line.
    Round 49: measured at the drawn size (10 + 1 pt) so the rule sees the larger text; delta=0.0 measures at the V13 size (record only)."""
    if page.width(text, 10.0, delta=delta) <= B["name_max_pt"] + 1.0: return [text]
    words, best = text.split(), None
    for k in range(1, len(words)):
        l1, l2 = " ".join(words[:k]), " ".join(words[k:])
        if page.width(l1, 10.0, delta=delta) <= B["name_max_pt"] + 1.0 and page.width(l2, 10.0, delta=delta) <= B["name_max_pt"] + 1.0: best = [l1, l2]
    assert best is not None, text; assert " ".join(best) == text
    return best
NEW_NAMES = [f for f in tail.feature if not (NICE[f] in BASE_NAMES or all(l in BASE_NAMES for l in wrap2(NICE[f], delta=0.0)))]
print("names new to the V13 sheet:", [(f, NICE[f]) for f in NEW_NAMES])
axb = page.axes(B["x0"] + XB, B["y0"], B["x1"] + XB, B["y1"]); axb.set_xlim(*B["xlim"]); axb.set_ylim(*B["ylim"])
yy = np.arange(N_ROWS)[::-1].astype(float)
rows = []; WRAP_CHANGED = []
LINES_ALL = [wrap2(NICE[f]) for f in tail.feature]
# round 49 nudge: a two-line name at 11 pt on the V26 line pitch (9.14 pt) lets a descender of line 1 touch a capital of line 2 (row 22,
# 'Sigma power in' over 'non-REM, central'); the two lines are set 1.22 pt further apart, symmetric about the same row centre, wherever both
# neighbouring rows are one-line names (between two consecutive two-line rows the V26 pitch stays: a wider pitch there would bring the lower
# row's first line within 0.1 pt of the upper row's second line)
TWO_DY_WIDE = (B["one_dy"] - 5.18, B["one_dy"] + 5.18)
def two_dy_for(i):
    prev_two = i - 1 >= 0 and len(LINES_ALL[i - 1]) == 2; next_two = i + 1 < len(LINES_ALL) and len(LINES_ALL[i + 1]) == 2
    return B["two_dy"] if (prev_two or next_two) else TWO_DY_WIDE
for i, r in tail.iterrows():
    g = float(r.g); y = row_y(i); col = C.FAMILY_COL[r.family]; name = NICE[r.feature]; lines = wrap2(name); lines_v26 = wrap2(name, delta=0.0); two_dy = two_dy_for(i)
    if lines != lines_v26: WRAP_CHANGED.append(dict(name=name, v26=lines_v26, new=lines, width_at_11=round(page.width(name, 10.0), 2), limit=B["name_max_pt"] + 1.0))
    assert abs((B["y0"] + (B["ylim"][1] - yy[i]) * PITCH) - y) < 0.01
    axb.plot([TAIL_G, g], [yy[i], yy[i]], color=B["leader_col"], lw=0.7, ls=(0, (1, 1.6)), zorder=2)
    axb.scatter([g], [yy[i]], s=B["dot_s"], color=col, zorder=3, edgecolors="white", linewidths=0.8, clip_on=False)
    if len(lines) == 1: T("b", B["name_x"] + XB, y + B["one_dy"], lines[0], 10.0, static=False)
    else:
        for ln, dy in zip(lines, two_dy): T("b", B["name_x"] + XB, y + dy, ln, 10.0, static=False)
    s = f"{g:.1f}"; T("b", B["val_right"] + XB, y + B["val_dy"], s, 9.5, ha="right", static=False)
    T("b", B["fam_x"] + XB, y + B["fam_dy"], r.family, 9.5, color=col, static=False)
    rec("b", s, B["val_right"] + XB, y + B["val_dy"], "ranking_v3.csv", f"dC x 1000 where feature == {r.feature} (rank {int(r['rank'])}), one decimal", g, ha="right",
        dot_page_xy=[round(bx(g) + XB, 3), round(y, 3)], colour=col, row=i + 1, name=name, name_lines=lines, family=r.family, feature=r.feature)
    rec("b", lines[0], B["name_x"] + XB, y + (B["one_dy"] if len(lines) == 1 else two_dy[0]), "ranking_v3.csv x NICE", f"row {i + 1} name ({r.feature})", name, kind="label", lines=lines, two_dy=list(two_dy))
    rec("b", r.family, B["fam_x"] + XB, y + B["fam_dy"], "measure_families.csv", f"family where feature == {r.feature}", r.family, kind="label")
    rows.append(dict(row=i + 1, feature=r.feature, rank=int(r["rank"]), gain_x1000=g, printed=s, name=name, lines=lines, lines_v26=lines_v26, two_dy=list(two_dy), family=r.family, colour=col, new_name=r.feature in NEW_NAMES,
                     folded_twin=next((f for f, b in pairs if b == r.feature), None)))
# the wrap rule at the V13 size must reproduce the base sheet's own two-line set on the rows it already had; at +1 pt more names may wrap (declared)
base_two = {"Slow-oscillation duration, frontal", "Heartbeat intervals counted, whole night", "Spindle to slow-oscillation coupling, occipital", "Sigma power in non-REM, central",
            "Spindle to slow-oscillation overlap, occipital", "Spindle to slow-oscillation coupling, central", "Slow-oscillation duration, central"}
present_names = {rw["name"] for rw in rows}
got_two_v26 = {rw["name"] for rw in rows if len(rw["lines_v26"]) == 2 and not rw["new_name"]}; assert got_two_v26 == {n for n in base_two if n in present_names}, (got_two_v26 ^ {n for n in base_two if n in present_names})
got_two = {rw["name"] for rw in rows if len(rw["lines"]) == 2}; assert got_two >= got_two_v26, (got_two, got_two_v26)
print("names wrapped at +1 pt that were one line on V26 (declared):", [w["name"] for w in WRAP_CHANGED])
for s_ in ("top", "right", "left"): axb.spines[s_].set_visible(False)
axb.set_xticks([5, 10, 15, 20]); axb.set_yticks([]); axb.tick_params(axis="x", length=3.0, width=0.8, labelbottom=False); axb.tick_params(axis="y", length=0)
assert tail.g.max() < B["xlim"][1] and tail.g.min() > B["xlim"][0]
for v in (5, 10, 15, 20):
    T("b", bx(v) + XB, B["xtick_base_old"] + GROW, str(v), 10.0, ha="center"); rec("b", str(v), bx(v) + XB, B["xtick_base_old"] + GROW, "axis", "x tick", v, kind="tick", ha="center")
T("b", B["title_old"][1] + XB, B["title_old"][2] + GROW, B["title_old"][0], 11.0)
T("b", B["fam_head"][0] + XB, B["fam_head"][1], "Family", 9.5)
DRAWN["b"]["geometry"] = dict(B, pitch=PITCH, grow_pt=GROW, grow_vs_v13_pt=H - H_V13, page=[W, H], page_v13=[W0, H_V13], x_shift_panel_b=XB, split_x=SPLIT, n_rows=N_ROWS, rows=rows,
                              twins_folded_above_tail=[f for f, b in pairs if float(RKI.loc[f, "g"]) > TAIL_G], twin_max_abs_dC_diff=(max(twin_diff.values()) if twin_diff else 0.0),
                              new_names={f: NICE[f] for f in NEW_NAMES if f in set(tail.feature)}, nice_all=NICE, top5_families=TOP5_FAMILIES, n_measures=N_MEASURES, wrap_changed_at_plus1=WRAP_CHANGED, two_dy_wide=list(TWO_DY_WIDE))

# ================================================================== letters (matplotlib bold, the base sheet's own way; b at x + XB), write
T("letters", LET["a"]["origin"][0], LET["a"]["origin"][1], "a", 13.0, weight="bold")
T("letters", LET["b"]["origin"][0] + XB, LET["b"]["origin"][1], "b", 13.0, weight="bold")
bad = C.banned_words([t["text"] for t in page.texts]); assert not bad, bad
raw_pdf = f"{WORK}/sheet_raw.pdf"; page.fig.savefig(raw_pdf, facecolor="white"); plt.close(page.fig)
built_pdf = f"{WORK}/sheet_built.pdf"
met = C.sheet_font_metrics(f"{C.BASE}/{S}.pdf"); done = C.align_font_metrics(raw_pdf, built_pdf, met); print("font metrics aligned:", done)
for k in "ab": json.dump(DRAWN[k], open(f"{VER}/{S}_{k}_drawn.json", "w"), indent=1, default=float)
json.dump(page.texts, open(f"{WORK}/placed_texts.json", "w"), indent=0, default=float)
json.dump(dict(sheet=S, builder=os.path.abspath(__file__), delta_pt=C.DELTA, page_built=[W, H], page_v13=[W0, H_V13], x_shift_panel_b=XB, split_x=SPLIT, grow_pt=GROW, n_rows=N_ROWS,
               xlabel_lines=XLAB_LINES, wrap_changed_at_plus1=WRAP_CHANGED, sizes_requested_plus_delta=sorted({t["size"] for t in page.texts}), built_pdf=built_pdf, built_sha256=C.sha256(built_pdf),
               sources={k: v for d in DRAWN.values() for k, v in d["sources"].items()}), open(f"{WORK}/build_record.json", "w"), indent=1, default=float)
print(f"wrote {built_pdf} page {W} x {H:.3f} ({GROW:+.3f} pt against the 28-row base, {H - H_V13:+.3f} pt against V13, {N_ROWS} rows, panel b at x + {XB}); bands {list(zip(BAND_N, [round(p, 2) for p in BAND_PCT]))}")
for rw in rows: print(f"  row {rw['row']:2d} rank {rw['rank']:3d} {rw['printed']:>5s} {rw['name']!r} {'(two lines)' if len(rw['lines']) == 2 else ''} {rw['family']}{' NEW NAME' if rw['new_name'] else ''}{' WRAP CHANGED' if rw['lines'] != rw['lines_v26'] else ''}")
