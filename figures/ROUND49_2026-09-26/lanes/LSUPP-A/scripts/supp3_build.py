#!/usr/bin/env python3
"""01_build for Supp_Fig03 (lane V14_L1_RANK, round 37: v8.1 numbers on the V13 design; repointed from the V7_L1_RANK builder).
v8.1: seven families (the PLM index is a family of its own), so panel a has twelve rows: the row pitch is kept, the panel grows by one
pitch at the foot and panel b, its letter and the page foot move down by the same amount (rules 3.3). The five oxygen measures of
panel b are the run_combo_v2 rule (the five largest oxygen gains in ranking_v3: T90, lowest saturation, T88, mean saturation, ODI3),
not the file's stale correlation_labels literal (NUMBERS_DEFECT: it still names the 4% desaturation index). The whole sheet is re-plotted from numbers/combination_v2.json (step 16,
rerun 2026-09-08) with matplotlib pinned to the round-30 base sheet's page points (base_text.json, probe_geom.py records).
Panel a: the eleven gain bars in the base sheet's row logic (best single, best 2 to 5 oxygen, the five one-per-family additions
sorted by gain, the six-family composite), the printed four-decimal gains, the dashed best-single rule.
Panel b: the 5 x 5 Spearman matrix among the five oxygen measurements. The file's correlation_labels list is in the retired
ranking_v2 order while its matrix is in the ranking_v3 order (proved by the in-lane Spearman in work/corr_check.json, see
NUMBERS_DEFECT_combination_v2_labels.md): the labels are taken from the ranking_v3 order, as the coordinator decided.
Colours are the base sheet's own (sampled from its vector records): #0288d1, #ccd1d6, #8a9099, ink #1a1d21, and the matrix ramp
that shows on the approved sheet (off-diagonal cells on the #e2e8ec to #1f4257 ramp, the diagonal at #0288d1).
Writes work/sheet_raw.pdf, Supp_Fig03.pdf (font metrics aligned to the base), work/placed_texts.json and verify/Supp_Fig03_{a,b}_drawn.json."""
import os, sys, json
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C

S = "Supp_Fig03"; SD = f"{C.LANE}/{S}"; WORK = f"{SD}/work"; VER = f"{SD}/verify"
os.makedirs(WORK, exist_ok=True); os.makedirs(VER, exist_ok=True)
plt = C.setup_matplotlib()
base_text = json.load(open(f"{WORK}/base_text.json")); W, H_V13 = base_text["page"]
assert [W, H_V13] == [484.72, 678.24], (W, H_V13)
SP = base_text["spans"]
def span_of(text, near_y=None):
    c = [s for s in SP if s["text"] == text]
    if near_y is not None: c = [s for s in c if abs(s["origin"][1] - near_y) < 8]
    assert len(c) == 1, (text, near_y, len(c)); return c[0]

# ------------------------------------------------------------------ numbers, with their sidecars
CMB_PATH = f"{C.NUM}/combination_v2.json"; CMB_SC = C.sidecar_ok(CMB_PATH); CMB_SHA = C.sha256(CMB_PATH)
CMB = json.load(open(C.hydrated(CMB_PATH)))
RK_PATH = f"{C.NUM}/ranking_v3.csv"; RK_SC = C.sidecar_ok(RK_PATH); RK_SHA = C.sha256(RK_PATH)
RK = pd.read_csv(C.hydrated(RK_PATH), comment="#"); N_MEASURES = len(RK); assert N_MEASURES == 141, N_MEASURES
FAMILY = pd.read_csv(C.hydrated(f"{C.NUM}/measure_families.csv")).set_index("feature")["family"]
CC = json.load(open(C.hydrated(f"{WORK}/corr_check.json")))
assert CC["sidecar"]["output_sha256"] == CMB_SC["output_sha256"] and CC["combination_v2_sha256"] == CMB_SHA, "corr_check.json is stale, rerun supp3_corr_check.py"

# ---- panel a rows, the base sheet's row logic (make_efig08_09_polish.py) with the v7 file
single = CMB["single"]; best_name = max(single, key=single.get); best_single = float(single[best_name])
assert best_name == "spo2_pct_below_90", best_name                      # T90 is the best single oxygen measurement (report if not)
c = CMB["combinations"]; ANCHOR = CMB["composite_terms"]["anchor"]; assert ANCHOR == best_name
TERMS = CMB["composite_terms"]["oxygen_plus_all_families"]; N_FAM = len(set(FAMILY.loc[RK.feature])); assert TERMS[0] == ANCHOR and len(TERMS) == N_FAM == 7, (TERMS, N_FAM)
FB = {v["feature"]: k for k, v in CMB["family_best"].items()}; assert set(TERMS) == set(FB), (TERMS, FB)
PLAIN = {"nrem_hr_bpm": "heart rate in non-REM", "F_spindle_n_spindles": "the spindle count", "C_nrem_sigma_mean": "non-REM sigma power",
         "resp_events_any_label": "the respiratory event count", "N2_min": "minutes of stage N2", "AHI": "the apnea-hypopnea index",
         "spo2_pct_below_90": "time below 90% saturation", "spo2_nadir_corrected": "the lowest saturation",
         # v8.1 family bests (named in the sheet's idiom, listed in the lane report)
         "REM_min": "minutes of REM sleep", "plm_index": "the limb movement index"}          # "periodic limb movement index" does not fit the label column (it would run into the bars)
OTHER = sorted([(f, float(c[f"oxygen_plus_{f}"])) for f in TERMS if f != ANCHOR], key=lambda kv: kv[1]); assert len(OTHER) == N_FAM - 1 == 6
for f, _ in OTHER: assert f in PLAIN, f"no plain-English name for {f}"
RKI = RK.set_index("feature")
for f, _ in OTHER: assert int(RKI.loc[f, "rank"]) == CMB["family_best"][FB[f]]["rank"], (f, "family_best rank disagrees with ranking_v3")
ITEMS = [("Best single oxygen measurement", best_single, C.BLUE, f"single.{best_name}"),
         ("Best 2 oxygen measurements", float(c["top2_oxygen"]), C.GREY_PALE, "combinations.top2_oxygen"),
         ("Best 3 oxygen measurements", float(c["top3_oxygen"]), C.GREY_PALE, "combinations.top3_oxygen"),
         ("Best 4 oxygen measurements", float(c["top4_oxygen"]), C.GREY_PALE, "combinations.top4_oxygen"),
         ("Best 5 oxygen measurements", float(c["top5_oxygen"]), C.GREY_PALE, "combinations.top5_oxygen")]
for f, v in OTHER: ITEMS.append((f"Oxygen and {PLAIN[f]}", v, C.GREY if v > best_single else C.GREY_PALE, f"combinations.oxygen_plus_{f}"))
ITEMS.append(("Oxygen and the best of every other family", float(c["oxygen_plus_all_families"]), C.GREY, "combinations.oxygen_plus_all_families"))
N_ITEMS = len(ITEMS); assert N_ITEMS == 5 + (N_FAM - 1) + 1 == 12 and abs(c["oxygen_plus_all_families"] - c["all_six_family_bests"]) < 1e-12
N_ITEMS_V13 = 11; GROW = (N_ITEMS - N_ITEMS_V13) * 20.339; H = H_V13 + GROW     # one row pitch (the base's 20.339 pt) per extra row, panel b and its letter move down by GROW

# ---- panel b matrix in the ranking_v3 order (the label defect), asserted against the in-lane Spearman
FIVE_NAME = {"spo2_pct_below_90": "Below 90%", "spo2_pct_below_88": "Below 88%", "spo2_nadir_corrected": "Lowest saturation",
             "odi3_total": "3% desaturations", "odi4_total": "4% desaturations", "spo2_mean": "Mean saturation"}
PLAIN8 = {"Below 90%": "Time below 90% saturation", "Below 88%": "Time below 88% saturation", "Lowest saturation": "Lowest saturation",
          "3% desaturations": "3% desaturation index", "4% desaturations": "4% desaturation index", "Mean saturation": "Mean saturation"}
# the run_combo_v2 rule (TOP5_OXY): the five oxygen measures with the largest gain in ranking_v3, in ranking order (data-driven, never typed)
ORDER = RK.assign(family=FAMILY.loc[RK.feature].values).query("family == 'Oxygenation'").sort_values("dC", ascending=False).feature.head(5).tolist()
assert ORDER == CC["ranking_v3_order"] == CC["five"] and all(f in FIVE_NAME for f in ORDER), (ORDER, CC["ranking_v3_order"])
assert CC["max_abs_diff_z_scale"]["under_ranking_v3_order"] < 1e-3, CC["max_abs_diff_z_scale"]
assert CC["file_labels_are_the_five"] is False, "the file's label list now names the five: re-check the label defect note"
assert CC["n_complete_cases"] == CMB["correlation_n"] == 15551
R = np.array(CMB["correlation_matrix"], float); assert R.shape == (5, 5) and np.allclose(R, R.T) and np.allclose(np.diag(R), 1.0) and (np.abs(R) <= 1).all()
# The file stores the matrix at three decimals. A cell whose third decimal is a 5 is a rounding tie for the two-decimal print, so the
# printed string is taken from the in-lane full-precision recomputation (the file's own method on the v7 tables, reproduced to within
# 0.0005 in every cell); everywhere else the two roundings agree, asserted. The tie cells are listed for the coordinator.
RF = np.array(CC["recomputed_z_scale_matrix_ranking_order_6dp"], float); assert np.abs(RF - R).max() < 5.1e-4, np.abs(RF - R).max()
PRINT = [[f"{RF[a, b]:.2f}" for b in range(5)] for a in range(5)]
TIES = [(a, b, float(R[a, b]), float(RF[a, b]), f"{R[a, b]:.2f}", PRINT[a][b]) for a in range(5) for b in range(5) if f"{R[a, b]:.2f}" != PRINT[a][b]]
for a, b, rf, rfull, sfile, sfull in TIES: assert abs(round(rf * 1000)) % 10 == 5, ("a two-decimal disagreement that is not a rounding tie", a, b, rf, rfull)
NAMES = [FIVE_NAME[f] for f in ORDER]; ROWNAMES = [PLAIN8[n] for n in NAMES]
assert R[ORDER.index("spo2_pct_below_90"), ORDER.index("spo2_nadir_corrected")] < 0, "T90 against the lowest saturation must be negative"
assert sorted(NAMES) != sorted(CMB["correlation_labels"]), "the file's label list is a stale literal that still names the 4% desaturation index (NUMBERS_DEFECT); it is not used"

# ================================================================== geometry from the base sheet (probe_geom.py, base_text.json)
A = dict(x0=207.36, y0=37.44, x1=463.124, y1=275.04, xmax=0.0345, pitch=20.339, c0=54.545, bar_h=12.61, val_dx=0.0006, lab_x=12.96,
         lab_dy=3.60, val_dy=2.585, tick_base=288.7, title_x=335.242, title_base=(302.7, 314.63))
A["y1"] += GROW; A["tick_base"] += GROW; A["title_base"] = (A["title_base"][0] + GROW, A["title_base"][1] + GROW)     # the panel grows at the foot by one row pitch
A["ppu"] = (A["x1"] - A["x0"]) / A["xmax"]                                   # 7413.4 pt per unit of gain
ylim_top = (N_ITEMS - 1) + (A["c0"] - A["y0"]) / A["pitch"]; ylim_bot = ylim_top - (A["y1"] - A["y0"]) / A["pitch"]
def ax_x(v): return A["x0"] + v * A["ppu"]
def row_c(i): return A["c0"] + A["pitch"] * i
assert abs(ax_x(0.01) - 281.495) < 0.02 and abs(ax_x(0.03) - 429.764) < 0.02, (ax_x(0.01), ax_x(0.03))
B = dict(x0=200.16, y0=375.84 + GROW, cell=40.32, lab_right=196.65, lab_dy=3.58, col_top=580.95 + GROW, col_dx=2.59, val_dy=2.59, white_over=0.80,
         ramp_lo="#e2e8ec", ramp_hi="#1f4257", vmin=0.45, vmax=1.0, diag=C.BLUE)      # panel b moves down by GROW with its letter
LET = {L["text"]: L for L in base_text["letters"]}; assert set(LET) == {"a", "b"}
LET_SHIFT = {"a": 0.0, "b": GROW}

from matplotlib import colors as mcolors
_CMAP = mcolors.LinearSegmentedColormap.from_list("house_blue", [B["ramp_lo"], B["ramp_hi"]]); _NORM = mcolors.Normalize(vmin=B["vmin"], vmax=B["vmax"])
def cell_fill(r):
    """The base sheet's cell colour rule (what shows): the old builder's 256-entry colormap ramp for every off-diagonal cell, the diagonal at #0288d1."""
    m = abs(r)
    if m >= 0.9995: return B["diag"]
    rgba = _CMAP(_NORM(m)); return "#%02x%02x%02x" % tuple(int(round(v * 255)) for v in rgba[:3])
def _close(h1, h2, tol=1):
    return all(abs(int(h1[i:i + 2], 16) - int(h2[i:i + 2], 16)) <= tol for i in (1, 3, 5))
# the sampled base cells must come out of the same rule within one 8-bit unit per channel (PDF colour rounding): 0.943 -> #335366, 0.978 -> #27485d,
# 0.691 -> #8c9fab, 0.745 -> #798f9c, 1.0 -> #0288d1 exactly
for r_, hx in ((0.943, "#335366"), (0.978, "#27485d"), (0.691, "#8c9fab"), (0.745, "#798f9c"), (0.903, "#425f71"), (0.666, "#95a7b1")):
    assert _close(cell_fill(r_), hx), (r_, cell_fill(r_), hx)
assert cell_fill(1.0) == "#0288d1"

page = C.Page(plt, W, H); page.fig.patch.set_alpha(1.0); page.fig.patch.set_facecolor("white")
DRAWN = {k: dict(sheet=S, panel=k, sources={"combination_v2.json": dict(path=CMB_PATH, sha256=CMB_SHA, sidecar=CMB_SC),
                                            "ranking_v3.csv": dict(path=RK_PATH, sha256=RK_SHA, sidecar=RK_SC, use="label order of the five oxygen measurements, family_best ranks")},
                 values=[]) for k in "ab"}
def rec(panel, text, x, baseline, source, key, value, kind="printed", **extra):
    DRAWN[panel]["values"].append(dict(text=text, x=round(x, 2), baseline=round(baseline, 2), source=source, key=key, value=value, kind=kind, **extra))
def T(panel, x, baseline, s, size, ha="left", color=C.INK, rotation=0, static=True, weight="normal", keep_size=False):
    page.text(x, baseline, s, size, ha=ha, color=color, rotation=rotation, weight=weight, record=dict(panel=panel, static=static), keep_size=keep_size)

# ================================================================== panel a
axa = page.axes(A["x0"], A["y0"], A["x1"], A["y1"]); axa.set_xlim(0, A["xmax"]); axa.set_ylim(ylim_bot, ylim_top)
ys = np.arange(len(ITEMS))[::-1].astype(float)
assert max(v for _, v, _, _ in ITEMS) < A["xmax"], "a bar leaves the axis"
axa.barh(ys, [v for _, v, _, _ in ITEMS], color=[cc for _, _, cc, _ in ITEMS], height=A["bar_h"] / A["pitch"], linewidth=0, zorder=2)
axa.axvline(best_single, color=C.INK, lw=0.9, ls=(0, (4, 3)), zorder=3)
rec("a", "dashed best-single rule", ax_x(best_single), 0, "combination_v2.json", f"single.{best_name} (max of single)", best_single, kind="drawn")
for i, (lab, v, col, key) in enumerate(ITEMS):
    yc = row_c(i); s = f"{v:.4f}"
    # the value with its opaque white pad (the base sheet's idiom: a text-layer string, never a stroke path effect)
    page.fig.text(ax_x(v + A["val_dx"]) / W, 1 - yc / H, s, fontsize=9.5 + C.Page.PT_PLUS, ha="left", va="center", color=C.INK,
                  bbox=dict(facecolor="white", edgecolor="none", pad=0.8), zorder=4)   # round 49: +1 pt (placed through fig.text, so the helper's +1 is added here)
    page.texts.append(dict(text=s, x=ax_x(v + A["val_dx"]), baseline=yc + A["val_dy"], size=9.5 + C.Page.PT_PLUS, ha="left", color=C.INK, record=dict(panel="a", static=False)))
    rec("a", s, ax_x(v + A["val_dx"]), yc + A["val_dy"], "combination_v2.json", key, v, bar_right_x=round(ax_x(v), 3), colour=col, row=i + 1)
    # round 49 exception: the twelve row labels of panel a stay at 10 pt. At 11 pt the composite label (204 pt wide from x 12.96) would run 9.6 pt into
    # the bars (x0 207.36) and the apnea-hypopnea label would break the 4 pt rule, no nudge can hold the longest one, so the whole column keeps its size.
    assert A["lab_x"] + page.width(lab, 10.0, keep_size=True) < A["x0"] - 4.0, (lab, "row label runs into the bars")
    T("a", A["lab_x"], yc + A["lab_dy"], lab, 10.0, static=(i < 5 or i == N_ITEMS - 1), keep_size=True); rec("a", lab, A["lab_x"], yc + A["lab_dy"], "combination_v2.json", f"row {i + 1} label ({key})", lab, kind="label")
for s_ in ("top", "right", "left"): axa.spines[s_].set_visible(False)
axa.set_xticks([0, 0.01, 0.02, 0.03]); axa.set_yticks([]); axa.tick_params(axis="x", length=3.0, width=0.8, labelbottom=False); axa.tick_params(axis="y", length=0)
for v, s_ in ((0, "0.00"), (0.01, "0.01"), (0.02, "0.02"), (0.03, "0.03")):
    T("a", ax_x(v), A["tick_base"], s_, 10.0, ha="center"); rec("a", s_, ax_x(v), A["tick_base"], "axis", "x tick", v, kind="tick", ha="center")
T("a", A["title_x"], A["title_base"][0], "Gain in prediction accuracy (held-out", 11.0, ha="center")
T("a", A["title_x"], A["title_base"][1], "concordance) over age and sex", 11.0, ha="center")
DRAWN["a"]["geometry"] = dict(A, ylim=[ylim_bot, ylim_top], rows=[dict(row=i + 1, label=l, value=v, colour=cc, key=k) for i, (l, v, cc, k) in enumerate(ITEMS)],
                              best_single=dict(feature=best_name, value=best_single), baseline=CMB["baseline"], n_outcomes=CMB["n_outcomes"], n_items=N_ITEMS, n_items_v13=N_ITEMS_V13,
                              grow_pt=GROW, page=[W, H], page_v13=[W, H_V13], letter_shift=LET_SHIFT, plain_names=PLAIN, families=N_FAM)

# ================================================================== panel b
n = 5; bx0, by0, cell = B["x0"], B["y0"], B["cell"]
axb = page.axes(bx0, by0, bx0 + n * cell, by0 + n * cell); axb.set_xlim(-0.5, n - 0.5); axb.set_ylim(n - 0.5, -0.5); axb.axis("off")
for a in range(n):
    for b in range(n):
        r = float(R[a, b]); fill = cell_fill(r)
        axb.add_patch(plt.Rectangle((b - 0.5, a - 0.5), 1, 1, facecolor=fill, edgecolor="none", linewidth=0))
        cx_ = bx0 + (b + 0.5) * cell; cy_ = by0 + (a + 0.5) * cell; s = PRINT[a][b].replace("-", "−")
        col = "white" if abs(r) > B["white_over"] else C.INK
        T("b", cx_, cy_ + B["val_dy"], s, 9.5, ha="center", color=col, static=False)
        rec("b", s, cx_, cy_ + B["val_dy"], "combination_v2.json", f"correlation_matrix[{a}][{b}] with rows and columns in the ranking_v3 order ({ORDER[a]} x {ORDER[b]})", r, ha="center",
            value_full_precision=float(RF[a, b]), rounding_tie=any(t[0] == a and t[1] == b for t in TIES), fill=fill, text_colour=col, cell_rect=[round(bx0 + b * cell, 3), round(by0 + a * cell, 3), round(bx0 + (b + 1) * cell, 3), round(by0 + (a + 1) * cell, 3)])
for a in range(n):
    T("b", B["lab_right"], by0 + (a + 0.5) * cell + B["lab_dy"], ROWNAMES[a], 10.0, ha="right")
    rec("b", ROWNAMES[a], B["lab_right"], by0 + (a + 0.5) * cell + B["lab_dy"], "ranking_v3.csv", f"row {a + 1} = {ORDER[a]} (ranking_v3 order)", ROWNAMES[a], kind="label", ha="right")
for b in range(n):
    w = page.width(NAMES[b], 10.0); x = bx0 + (b + 0.5) * cell + B["col_dx"]; y = B["col_top"] + w
    T("b", x, y, NAMES[b], 10.0, rotation=90)
    rec("b", NAMES[b], x, y, "ranking_v3.csv", f"column {b + 1} = {ORDER[b]} (ranking_v3 order)", NAMES[b], kind="label", rotated=True)
DRAWN["b"]["geometry"] = dict(B, order=ORDER, names=NAMES, rownames=ROWNAMES, matrix=R.tolist(), file_labels=CMB["correlation_labels"],
                              mean_abs_correlation=CMB["mean_abs_correlation"], corr_check=os.path.basename(f"{WORK}/corr_check.json"),
                              rounding_ties=[dict(row=a, col=b, file_value=rf, full_precision=rfull, printed_from_file=sf, printed=sp) for a, b, rf, rfull, sf, sp in TIES])

# ================================================================== letters (matplotlib bold, the base sheet's own way), write
for ch in "ab":
    T("letters", LET[ch]["origin"][0], LET[ch]["origin"][1] + LET_SHIFT[ch], ch, 13.0, weight="bold")
bad = C.banned_words([t["text"] for t in page.texts]); assert not bad, bad
raw = f"{WORK}/sheet_raw.pdf"; page.fig.savefig(raw, facecolor="white"); plt.close(page.fig)
met = C.sheet_font_metrics(f"{C.BASE}/{S}.pdf"); done = C.align_font_metrics(raw, f"{SD}/{S}.pdf", met); print("font metrics aligned:", done)
for k in "ab": json.dump(DRAWN[k], open(f"{VER}/{S}_{k}_drawn.json", "w"), indent=1, default=float)
json.dump(page.texts, open(f"{WORK}/placed_texts.json", "w"), indent=0, default=float)
print(f"wrote {SD}/{S}.pdf page {W} x {H:.3f} ({GROW:+.3f} pt against V13 for {N_ITEMS} rows); rows:", [(l, round(v, 5)) for l, v, _, _ in ITEMS]); print("matrix order:", NAMES); print(np.round(R, 3))
