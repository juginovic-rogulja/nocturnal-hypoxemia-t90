#!/usr/bin/env python3
"""01_build for ED_Fig01 (lane V14_L1_RANK, round 37: v8.1 numbers on the V13 design; repointed from the V7_L1_RANK builder). The whole sheet is re-plotted with matplotlib pinned to the round-30 base
sheet's page points (base_text.json, probe_geom.py records).
Panel a: the four T90 band shares (numbers/results_v2.json cohort.t90_bands_pct and _n, step 12; counts re-derived from the v7 frozen
table through numbers/cohort_spec.apply_cohort and asserted; numbers/cohorts.json n_analysis, step 01).
Panel b: every measurement whose held-out gain exceeds 0.003 in numbers/ranking_v3.csv (step 05, rerun 2026-09-08), one named row
per distinct measurement (a site z-scored twin of a measurement already listed is folded into it: in the v7 ranking the sixty twins
differ from their originals by under 1.4e-6 in dC, print the same one-decimal gain and sit at adjacent ranks, all asserted), sorted by
gain, the printed gain in units of 0.001, the family column in the family colours, the dotted leaders.
The row set grows from 28 to 30 named rows: the coordinator's standing answer is to grow the page by the needed height with the row
pitch kept (2 rows x 20.0547 pt = 40.11 pt), so the page box becomes 518.74 x 717.92 pt and everything below the rows (x axis, tick
labels, axis title) moves down by that amount. Letters, panel a, the Family header and the row origin stay where they are.
Two measurements enter the row set without a plain-English name on the old sheet and get one here in the sheet's own idiom:
O_spindle_n_spindles = "Spindle count, occipital", C_nrem_theta_kurt = "Theta power kurtosis in non-REM, central" (listed in the report).
Writes work/sheet_raw.pdf, ED_Fig01.pdf (font metrics aligned to the base), work/placed_texts.json, verify/ED_Fig01_{a,b}_drawn.json."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import os, sys, json
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C

S = "ED_Fig01"; SD = f"{C.LANE}/{S}"; WORK = f"{SD}/work"; VER = f"{SD}/verify"
os.makedirs(WORK, exist_ok=True); os.makedirs(VER, exist_ok=True)
plt = C.setup_matplotlib()
base_text = json.load(open(f"{WORK}/base_text.json")); W, H_V13 = base_text["page"]
assert [W, H_V13] == [518.74, 717.92], (W, H_V13)          # the V13 sheet (30 named rows)
H_BASE = 677.81                                                # the 28-row geometry every constant below is pinned to (the round-28 base the R30 builder measured)
assert abs(H_V13 - (H_BASE + 2 * (616.124 - 74.648) / 27.0)) < 0.02, "the V13 page is not the 28-row base grown by two pitches"
SP = base_text["spans"]
LET = {L["text"]: L for L in base_text["letters"]}; assert set(LET) == {"a", "b"}

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
# v8.1: the site z-scored twins left the list (197 -> 141, decision 3), so no fold is needed; the fold logic stays for the record and is asserted on the data
assert len(pairs) == 0, ("v8.1 ranking still carries siteZ twins", pairs[:5])
dup = [f for f, _ in pairs]; distinct = RK[~RK.feature.isin(dup)].copy(); assert len(distinct) == N_MEASURES
K = int((RK.g > TAIL_G).sum()); assert float(RK.sort_values("rank").g.iloc[K]) <= TAIL_G and float(RK.sort_values("rank").g.iloc[K - 1]) > TAIL_G
tail = distinct[distinct.g > TAIL_G].sort_values("g", ascending=False).reset_index(drop=True)
N_ROWS = len(tail); assert N_ROWS == K - sum(1 for f, b in pairs if float(RKI.loc[f, "g"]) > TAIL_G), (N_ROWS, K)
TOP5_FAMILIES = sorted(set(RK[RK["rank"] <= 5].family)); print("families in the top 5:", TOP5_FAMILIES)   # v7 asserted Oxygenation only; v8.1: NREM heart rate is rank 4 (reported, not asserted)
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
        # new in the v7 row set (no name on the old sheet), named in the sheet's idiom, listed in the lane report
        "O_spindle_n_spindles": "Spindle count, occipital", "C_nrem_theta_kurt": "Theta power kurtosis in non-REM, central",
        # new in v8.1 (decision 3 added the measurement): named in the sheet's idiom, listed in the lane report
        "hypoxic_burden": "Hypoxic burden", "plm_index": "Periodic limb movement index", "sleep_onset_latency_min": "Sleep onset latency", "rem_latency_min": "REM latency"}
BASE_NAMES = {s["text"] for s in SP}
NEW_NAMES = None          # set after wrap2 is defined: a name is on the V13 sheet as one span or as its two wrapped lines
missing = sorted(f for f in tail.feature if f not in NICE); assert not missing, f"no plain-English name for: {missing}"
print(f"{N_ROWS} rows with gain above 0.003 (V13 printed 30)")

# ================================================================== geometry from the base sheet
PITCH = (616.124 - 74.648) / 27.0                                   # 20.0547 pt, the base's row pitch (28 dot centres)
GROW = (N_ROWS - 28) * PITCH; H = H_BASE + GROW
A = dict(x0=79.59, y0=33.122, x1=147.99, y1=263.522, xmax=60.0, ylim=(-0.72, 3.72), bar_h=0.60, ylab_right=73.08, ylab_dy=3.58, xtick_base=277.18,
         xlabel=(25.68, 293.0), ylabel_lines=[("Time below 90% saturation,", 20.84, 216.06), ("% of the recording", 32.25, 193.01)])
B = dict(x0=292.59, y0=57.602, x1=371.05, xlim=(2.85, 24.5), y_row0=74.648, name_x=168.15, one_dy=3.58, two_dy=(-0.99, 8.15), val_right=395.35, val_dy=2.58,
         fam_x=402.55, fam_dy=2.58, fam_head=(402.55, 51.17), xtick_base_old=646.83, title_old=("Gain in held-out concordance, in units of 0.001", 217.96, 661.69),
         name_max_pt=1.67 * 72.0, dot_s=28.0, leader_col="#b3dcf2")
B["y1"] = B["y0"] + (N_ROWS + 0.7) * PITCH; B["ylim"] = (-0.85, N_ROWS - 1 + 0.85)
assert abs((B["y0"] + (28 + 0.7) * PITCH) - 633.17) < 0.02                 # the base's axes bottom for 28 rows
ppu = (B["x1"] - B["x0"]) / (B["xlim"][1] - B["xlim"][0])
def bx(g): return B["x0"] + (g - B["xlim"][0]) * ppu
def row_y(i): return B["y_row0"] + PITCH * i
assert abs(bx(5) - 300.382) < 0.02 and abs(bx(20) - 354.742) < 0.02 and abs(bx(TAIL_G) - 293.134) < 0.02
A_ppu = (A["x1"] - A["x0"]) / A["xmax"]; A_pitch = (A["y1"] - A["y0"]) / (A["ylim"][1] - A["ylim"][0])
def band_y(k): return A["y0"] + (A["ylim"][1] - (3 - k)) * A_pitch                   # band k (0 = top) centre
assert abs(band_y(0) - 70.484) < 0.02 and abs(band_y(3) - 226.16) < 0.02 and abs(A_ppu * 20 - 22.8) < 0.02

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
T("a", *A["xlabel"], "Share of the cohort, %", 11.0)
for s_, x, y in A["ylabel_lines"]: T("a", x, y, s_, 11.0, rotation=90)
DRAWN["a"]["geometry"] = dict(A, ppu=A_ppu, pitch=A_pitch, bands=[dict(label=l, key=k, n=n, pct=p, colour=c_) for (l, _, k), n, p, c_ in zip(BANDS, BAND_N, BAND_PCT, BAND_COL)])

# ================================================================== panel b
def wrap2(text):
    """The base sheet's wrap rule: at most two lines of 1.67 in, words unchanged, the longest first line that leaves a fitting second line."""
    # the base's one-line limit is 1.67 in measured on the old builder's 100 dpi renderer; 'Time below 90% saturation' sits exactly on it (120.24 pt),
    # so the fit test allows one 72 dpi pixel (1 pt) of measurement rounding
    if page.width(text, 10.0) <= B["name_max_pt"] + 1.0: return [text]
    words, best = text.split(), None
    for k in range(1, len(words)):
        l1, l2 = " ".join(words[:k]), " ".join(words[k:])
        if page.width(l1, 10.0) <= B["name_max_pt"] + 1.0 and page.width(l2, 10.0) <= B["name_max_pt"] + 1.0: best = [l1, l2]
    assert best is not None, text; assert " ".join(best) == text
    return best
NEW_NAMES = [f for f in tail.feature if not (NICE[f] in BASE_NAMES or all(l in BASE_NAMES for l in wrap2(NICE[f])))]
print("names new to the V13 sheet:", [(f, NICE[f]) for f in NEW_NAMES])
axb = page.axes(B["x0"], B["y0"], B["x1"], B["y1"]); axb.set_xlim(*B["xlim"]); axb.set_ylim(*B["ylim"])
yy = np.arange(N_ROWS)[::-1].astype(float)
rows = []
for i, r in tail.iterrows():
    g = float(r.g); y = row_y(i); col = C.FAMILY_COL[r.family]; name = NICE[r.feature]; lines = wrap2(name)
    assert abs((B["y0"] + (B["ylim"][1] - yy[i]) * PITCH) - y) < 0.01
    axb.plot([TAIL_G, g], [yy[i], yy[i]], color=B["leader_col"], lw=0.7, ls=(0, (1, 1.6)), zorder=2)
    axb.scatter([g], [yy[i]], s=B["dot_s"], color=col, zorder=3, edgecolors="white", linewidths=0.8, clip_on=False)
    if len(lines) == 1: T("b", B["name_x"], y + B["one_dy"], lines[0], 10.0, static=False)
    else:
        for ln, dy in zip(lines, B["two_dy"]): T("b", B["name_x"], y + dy, ln, 10.0, static=False)
    s = f"{g:.1f}"; T("b", B["val_right"], y + B["val_dy"], s, 9.5, ha="right", static=False)
    T("b", B["fam_x"], y + B["fam_dy"], r.family, 9.5, color=col, static=False)
    rec("b", s, B["val_right"], y + B["val_dy"], "ranking_v3.csv", f"dC x 1000 where feature == {r.feature} (rank {int(r['rank'])}), one decimal", g, ha="right",
        dot_page_xy=[round(bx(g), 3), round(y, 3)], colour=col, row=i + 1, name=name, name_lines=lines, family=r.family, feature=r.feature)
    rec("b", lines[0], B["name_x"], y + (B["one_dy"] if len(lines) == 1 else B["two_dy"][0]), "ranking_v3.csv x NICE", f"row {i + 1} name ({r.feature})", name, kind="label", lines=lines)
    rec("b", r.family, B["fam_x"], y + B["fam_dy"], "measure_families.csv", f"family where feature == {r.feature}", r.family, kind="label")
    rows.append(dict(row=i + 1, feature=r.feature, rank=int(r["rank"]), gain_x1000=g, printed=s, name=name, lines=lines, family=r.family, colour=col, new_name=r.feature in NEW_NAMES,
                     folded_twin=next((f for f, b in pairs if b == r.feature), None)))
# the wrap rule must reproduce the base sheet's own two-line set on the rows it already had
base_two = {"Slow-oscillation duration, frontal", "Heartbeat intervals counted, whole night", "Spindle to slow-oscillation coupling, occipital", "Sigma power in non-REM, central",
            "Spindle to slow-oscillation overlap, occipital", "Spindle to slow-oscillation coupling, central", "Slow-oscillation duration, central"}
present_names = {rw["name"] for rw in rows}
got_two = {rw["name"] for rw in rows if len(rw["lines"]) == 2 and not rw["new_name"]}; assert got_two == {n for n in base_two if n in present_names}, (got_two ^ {n for n in base_two if n in present_names})
for s_ in ("top", "right", "left"): axb.spines[s_].set_visible(False)
axb.set_xticks([5, 10, 15, 20]); axb.set_yticks([]); axb.tick_params(axis="x", length=3.0, width=0.8, labelbottom=False); axb.tick_params(axis="y", length=0)
assert tail.g.max() < B["xlim"][1] and tail.g.min() > B["xlim"][0]
for v in (5, 10, 15, 20):
    T("b", bx(v), B["xtick_base_old"] + GROW, str(v), 10.0, ha="center"); rec("b", str(v), bx(v), B["xtick_base_old"] + GROW, "axis", "x tick", v, kind="tick", ha="center")
T("b", B["title_old"][1], B["title_old"][2] + GROW, B["title_old"][0], 11.0)
T("b", *B["fam_head"], "Family", 9.5)
DRAWN["b"]["geometry"] = dict(B, pitch=PITCH, grow_pt=GROW, grow_vs_v13_pt=H - H_V13, page=[W, H], page_v13=[W, H_V13], n_rows=N_ROWS, rows=rows, twins_folded_above_tail=[f for f, b in pairs if float(RKI.loc[f, "g"]) > TAIL_G],
                              twin_max_abs_dC_diff=(max(twin_diff.values()) if twin_diff else 0.0), new_names={f: NICE[f] for f in NEW_NAMES if f in set(tail.feature)}, nice_all=NICE, top5_families=TOP5_FAMILIES, n_measures=N_MEASURES)

# ================================================================== letters (matplotlib bold, the base sheet's own way), write
for ch in "ab": T("letters", LET[ch]["origin"][0], LET[ch]["origin"][1], ch, 13.0, weight="bold")
bad = C.banned_words([t["text"] for t in page.texts]); assert not bad, bad
raw_pdf = f"{WORK}/sheet_raw.pdf"; page.fig.savefig(raw_pdf, facecolor="white"); plt.close(page.fig)
met = C.sheet_font_metrics(f"{C.BASE}/{S}.pdf"); done = C.align_font_metrics(raw_pdf, f"{SD}/{S}.pdf", met); print("font metrics aligned:", done)
for k in "ab": json.dump(DRAWN[k], open(f"{VER}/{S}_{k}_drawn.json", "w"), indent=1, default=float)
json.dump(page.texts, open(f"{WORK}/placed_texts.json", "w"), indent=0, default=float)
print(f"wrote {SD}/{S}.pdf page {W} x {H:.3f} ({GROW:+.3f} pt against the 28-row base, {H - H_V13:+.3f} pt against V13, {N_ROWS} rows); bands {list(zip(BAND_N, [round(p, 2) for p in BAND_PCT]))}")
for rw in rows: print(f"  row {rw['row']:2d} rank {rw['rank']:3d} {rw['printed']:>5s} {rw['name']!r} {'(two lines)' if len(rw['lines']) == 2 else ''} {rw['family']}{' NEW NAME' if rw['new_name'] else ''}")
