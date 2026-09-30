#!/usr/bin/env python3
"""Supp_Fig04, V14 lane L2 (round 37), step 1: re-plot the four data panels from numbers/lag_ladder.json (step 130, v8.1), pinned to
the V13 sheet's probed geometry (repointed copy of the round-30 builder, byte copy beside as 01_build_PRE_V8_1.py; panels c and d at the
round-31 LTEXT boxes the V13 sheet carries).

Data wiring = efig14_polish.py (a, b) and split_lagladder_polish.py (c, d), every rule carried over: rows = the conditions with a lower
95% bound above 1 at lag 0 that are not a control and not circular, sorted by the lag-0 hazard ratio, ceil(n/2) in panel a and the rest
in panel b, then the five controls under "Negative controls"; labels "<condition> (<events surviving the 5-year landmark>)"; rungs 0 circle,
1 square, 2 diamond, 5 triangle, open when lag5_informative is false, "No estimate" where not estimable; panel c bars and line, "<survive>/<n>";
panel d trend classes. ROW COUNT under v8.1: panel_n = 35 (V13: 32), so panel a holds 18 rows and panel b 17 + header + 5 (V13: 16 and 22).
GROWTH RULE (V14 rules 3.3): the row pitch, fonts and markers are kept, panel a grows by 2 rows at the foot and panel b by 1 row, the strings
under each grow with their panel, the c and d band and the letters c, d move down by the larger growth, the page grows by the same amount.
usage: 01_build.py            -> work/panels.pdf, verify/Supp_Fig04_drawn.json (v8.1)
       01_build.py --old      -> work/panels_OLD.pdf, verify/positive_control_drawn.json (the v7 snapshot the V13 sheet drew: growth 0)"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import collections, csv, json, math, os, sys
import matplotlib
matplotlib.use("Agg")
from matplotlib import font_manager as fm
for _f in ("Arial.ttf", "Arial Bold.ttf", "Arial Italic.ttf"): fm.fontManager.addfont(f"{paths.FONT_DIR}/{_f}")
import matplotlib.pyplot as plt, matplotlib.ticker as mticker, numpy as np, pandas as pd
from matplotlib.markers import MarkerStyle
from matplotlib.patches import PathPatch
from matplotlib.transforms import Affine2D
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND37_2026-09-14/figures/_common")
from v14lib import NUM, OLD_V7, hydrated, sha256, sidecar_v8

OLDMODE = "--old" in sys.argv
LANE = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); WORK = os.path.join(LANE, "work"); VER = os.path.join(LANE, "verify"); os.makedirs(VER, exist_ok=True)
NUMDIR = f"{OLD_V7}/numbers" if OLDMODE else NUM
DATA = f"{NUMDIR}/lag_ladder.json"; DATA_CSV = f"{NUMDIR}/lag_ladder.csv"
OUT_PDF = os.path.join(WORK, "panels_OLD.pdf" if OLDMODE else "panels.pdf"); OUT_JSON = os.path.join(VER, "positive_control_drawn.json" if OLDMODE else "Supp_Fig04_drawn.json")
for p in (DATA, DATA_CSV, os.path.join(WORK, "base_geometry.json"), os.path.join(WORK, "base_text.json")): hydrated(p)
G = json.load(open(os.path.join(WORK, "base_geometry.json"))); BT = json.load(open(os.path.join(WORK, "base_text.json")))
SHA = sha256(DATA); SHA_CSV = sha256(DATA_CSV)
if OLDMODE: SRC = dict(file=DATA, sha256=SHA, csv=DATA_CSV, csv_sha256=SHA_CSV, sidecar=None, note="v7 snapshot, positive control only")
else:
    sc, scc = sidecar_v8(DATA), sidecar_v8(DATA_CSV)
    SRC = dict(file=DATA, sha256=SHA, csv=DATA_CSV, csv_sha256=SHA_CSV, sidecar=DATA + ".provenance.json", step=sc["step"], step_name=sc["name"], step_ended=sc["status_ok"], output_mtime=sc["output_mtime"], v8_inputs=sc["v8_inputs"], csv_step=scc["step"])

# ------------------------------------------------------------------ page geometry (pt, top-down), from the V13 probe
PAGE_W, PAGE_H13 = G["page"]; assert (round(PAGE_W, 3), round(PAGE_H13, 3)) == (1021.45, 989.728), G["page"]
AX13 = G["axes_box"]; PITCH = G["row_pitch"]["a"]; assert abs(PITCH - G["row_pitch"]["b"]) < 0.01, G["row_pitch"]
ROWS13 = dict(a=G["row_pitch"]["rows_a"], b=G["row_pitch"]["rows_b"])
XLO13, XHI13 = 0.62, 4.30; XTICKS = [0.7, 1.0, 1.5, 2.0, 3.0, 4.0]; RUNG_OFF = [0.36, 0.12, -0.12, -0.36]   # the V13 (generator) fixed limits; extended below by the 0.96 x minimum rule when a v8.1 interval runs off
XLABEL = "Hazard ratio for a new diagnosis per 1 SD of sleep T90 (95% CI)"; CTRL_HEAD = "Negative controls"; GUTTER_PAD = 4.0
LAGS = ["0", "1", "2", "5"]; LAGNAME = ["All follow-up", "1-year landmark", "2-year landmark", "5-year landmark"]
QUIET_1 = "5-year landmark, too few events"; QUIET_2 = ("5-year landmark,", "too few events")
KEY_LINE = "line with circles = share of the extra risk that remains at each landmark"; KEY_BARS = "bars = share of disease events the landmark deletes"
C_XLABEL, C_YLABEL, D_XLABEL = "Landmark, years of diagnoses deleted", "Median across conditions, %", "Conditions"; D_LABELS = ["Stable", "Attenuating", "Collapsing", "Reversing"]
INK, GRID, PRIMARY, SECOND, COMPARE, PALE2, WHITE = "#1a1d21", "#eef0f1", "#0288d1", "#3f9fd8", "#8a9099", "#aeb5bc", "#ffffff"
assert set(G["colours_vector"]["fills"]) <= {INK, GRID, PRIMARY, SECOND, COMPARE, PALE2, WHITE}, G["colours_vector"]
CI_LW, MEDGE = 1.7, 0.8; MARK = ["o", "s", "D", "^"]; SIZE = [28.0, 23.0, 26.0, 30.0]; MSIZE = [6.0, 5.4, 5.7, 6.3]
PT_PLUS = 1.0   # round 49 (Alen, 2026-09-26): every text one point larger, lane LSUPP-A (ticks and labels 10 -> 11, axis titles 11 -> 12, annotations 9.5 -> 10.5, "No estimate" 9 -> 10, letters 13 -> 14)
TICK_PT, TITLE_PT, ANN_PT, FLOOR_PT, LETTER_PT = 10.0 + PT_PLUS, 11.0 + PT_PLUS, 9.5 + PT_PLUS, 9.0 + PT_PLUS, 13.0 + PT_PLUS
RC = {"font.family": "Arial", "font.size": TICK_PT, "axes.labelsize": TITLE_PT, "axes.titlesize": TITLE_PT, "xtick.labelsize": TICK_PT, "ytick.labelsize": TICK_PT, "legend.fontsize": TICK_PT, "legend.frameon": False,
      "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8, "xtick.major.size": 3.0, "ytick.major.size": 3.0, "xtick.direction": "out", "ytick.direction": "out",
      "axes.edgecolor": INK, "text.color": INK, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK, "axes.spines.top": False, "axes.spines.right": False, "axes.facecolor": WHITE,
      "figure.facecolor": "none", "savefig.facecolor": "none", "savefig.bbox": "standard", "savefig.pad_inches": 0.0, "figure.dpi": 72, "savefig.dpi": 72, "pdf.fonttype": 42, "mathtext.default": "regular",
      "axes.unicode_minus": False, "lines.solid_capstyle": "projecting"}

# ------------------------------------------------------------------ data, asserted against the file and its csv twin
J = json.load(open(DATA)); PC = J["per_condition"]; PANEL_N = int(J["panel_n"])
assert J["lags"] == [0, 1, 2, 5] and J["controls_n"] == 5, (J["lags"], J["controls_n"])
panel_set = {k for k, v in PC.items() if (not v["negative_control"]) and (not v["circular"]) and v["ci_lag0"][0] > 1}
assert len(panel_set) == PANEL_N, (len(panel_set), PANEL_N)
_d = pd.read_csv(DATA_CSV); _w0 = _d[_d.lag_years == 0]
panel = list(_w0[(_w0.lo > 1) & ~_w0.negative_control & ~_w0.circular].sort_values("hr", ascending=False).key)
ctrl = list(_w0[_w0.negative_control].sort_values("hr", ascending=False).key)
assert set(panel) == panel_set and all(PC[panel[i]]["hr_lag0"] >= PC[panel[i + 1]]["hr_lag0"] for i in range(PANEL_N - 1))
_byhr = collections.defaultdict(list)
for k in panel: _byhr[PC[k]["hr_lag0"]].append(k)
TIES = [dict(hr_lag0=v, drawn_order=[PC[k]["condition"] for k in ks], lower_ci=[PC[k]["ci_lag0"][0] for k in ks]) for v, ks in _byhr.items() if len(ks) > 1]
assert len(ctrl) == 5 and set(ctrl) == set(J["controls"].keys()), (ctrl, sorted(J["controls"]))
assert len({PC[k]["hr_lag0"] for k in ctrl}) == 5, "a tie among the controls"
half = int(math.ceil(PANEL_N / 2))                                  # efig14_polish.py line 144
ROWS_A = panel[:half]; ROWS_B = panel[half:] + [None] + ctrl
LAB = {k: PC[k]["condition"] for k in panel + ctrl}; E5 = {k: int(PC[k]["events_by_lag"]["5"]) for k in panel + ctrl}
# row labels must fit the gutter left of each ladder (page or sheet edge to the axes left minus the 4 pt pad); a label that does not fit
# takes the display label Main_Fig2 uses for the same outcome (owner question, recorded), never a smaller type
import fitz as _fz; _ARIAL = _fz.Font(fontfile=paths.FONT_ARIAL)
LABEL_SHORT = {"Ventricular arrhythmia or cardiac arrest": "Ventricular arrhythmia/arrest"}
SHEET_LEFT = dict(a=0.0, b=522.724)
DISPLAY = {}; LABEL_NOTES = []
for _p, _rows in (("a", panel[:int(math.ceil(PANEL_N / 2))]), ("b", panel[int(math.ceil(PANEL_N / 2)):] + ctrl)):
    gutter = AX13[_p][0] - 4.0 - SHEET_LEFT[_p] - 2.0
    for k in _rows:
        full = f"{LAB[k]} ({E5[k]:,})"; w = _ARIAL.text_length(full, TICK_PT)   # round 49: the rule sees the 11 pt label
        if w <= gutter: DISPLAY[k] = LAB[k]
        else:
            assert LAB[k] in LABEL_SHORT, ("row label wider than the gutter and no display label defined", LAB[k], w, gutter)
            DISPLAY[k] = LABEL_SHORT[LAB[k]]; LABEL_NOTES.append(dict(condition=LAB[k], display=DISPLAY[k], width_pt=round(w, 1), gutter_pt=round(gutter, 1), panel=_p))
            assert _ARIAL.text_length(f"{DISPLAY[k]} ({E5[k]:,})", TICK_PT) <= gutter, (DISPLAY[k], gutter)
if LABEL_NOTES: print("display labels (owner question, as Main_Fig2):", LABEL_NOTES)
base_labels = {s["text"] for s in BT["spans"] if s["size"] == 10.0 and s["font"] == "ArialMT" and "(" in s["text"] and ")" in s["text"] and s["text"][-1] == ")" and not s["text"].startswith("Hazard")}
new_labels = {f"{DISPLAY[k]} ({E5[k]:,})" for k in panel + ctrl}
LABEL_DELTA = dict(entering=sorted(new_labels - base_labels), leaving=sorted(base_labels - new_labels))
rows_csv = list(csv.DictReader(open(DATA_CSV))); CSVI = {(r["key"], str(int(float(r["lag_years"])))): r for r in rows_csv}
for k in panel + ctrl:
    for L in LAGS:
        r = CSVI[(k, L)]; assert r["condition"] == LAB[k]; est = r["estimable"] == "True"; assert est == (PC[k]["hr_by_lag"][L] is not None), (k, L)
        if est: assert abs(float(r["hr"]) - PC[k]["hr_by_lag"][L]) < 5e-4 and abs(float(r["lo"]) - PC[k]["ci_by_lag"][L][0]) < 5e-4 and abs(float(r["hi"]) - PC[k]["ci_by_lag"][L][1]) < 5e-4, (k, L)
        assert int(r["events"]) == int(PC[k]["events_by_lag"][L]), (k, L)
        if L == "5": assert (r["lag5_informative"] == "True") == bool(PC[k]["lag5_informative"]), (k, L)
noest = {(k, L) for k in panel + ctrl for L in LAGS if PC[k]["hr_by_lag"][L] is None}
_los = [PC[k]["ci_by_lag"][L][0] for k in panel + ctrl for L in LAGS if PC[k]["ci_by_lag"][L]]; _his = [PC[k]["ci_by_lag"][L][1] for k in panel + ctrl for L in LAGS if PC[k]["ci_by_lag"][L]]
XLO = min(XLO13, round(0.96 * min(_los), 4)); XHI = max(XHI13, round(1.04 * max(_his), 4))   # V13 limits unless an interval runs off (ED_Fig02's rule of round 37): declared
print(f"x limits {XLO} to {XHI} (V13 {XLO13} to {XHI13}; smallest lower bound {min(_los)}, largest upper {max(_his)})")
for k in panel + ctrl:
    for L in LAGS:
        hr, ci = PC[k]["hr_by_lag"][L], PC[k]["ci_by_lag"][L]
        if hr is None: assert ci is None and (k, L) in noest
        else:
            assert ci[0] <= hr <= ci[1] and XLO < ci[0] and ci[1] < XHI, (k, L, hr, ci, "an interval would run off the axis")
            assert PC[k]["significant_by_lag"][L] == (ci[0] > 1.0 or ci[1] < 1.0), (k, L)
assert all(not PC[k]["lag5_informative"] for k in ctrl)
n_informative = sum(1 for k in panel if PC[k]["hr_by_lag"]["5"] is not None and PC[k]["lag5_informative"])
assert n_informative == J["by_lag"]["5"]["n_informative"] == J["lag5_verdict"]["n_informative"]
BL = J["by_lag"]
rem = [0.0] + [BL[L]["median_pct_events_removed"] for L in "125"]; ret = [100.0] + [100.0 * BL[L]["median_excess_retained"] for L in "125"]; surv = [PANEL_N] + [BL[L]["survive"] for L in "125"]
assert all(BL[L]["of"] == PANEL_N for L in "125") and all(0 <= s <= PANEL_N for s in surv) and all(0 <= v <= 100 for v in rem) and rem[1] < rem[2] < rem[3] and all(0 < v <= 100 for v in ret)
for L in "125": assert sum(1 for k in panel if PC[k]["hr_by_lag"][L] is not None and PC[k]["ci_by_lag"][L][0] > 1.0) == BL[L]["survive"], L
tr = [J["trends"][k]["trend"] for k in panel]; assert all(J["trends"][k]["trend"] == PC[k]["trend"] for k in panel)
assert all(t in ("stable", "stable, strengthening", "attenuating", "collapsing", "reversing") for t in tr), set(tr)
counts = [("Stable", sum(t.startswith("stable") for t in tr), PRIMARY), ("Attenuating", sum(t == "attenuating" for t in tr), SECOND), ("Collapsing", sum(t == "collapsing" for t in tr), COMPARE), ("Reversing", sum(t == "reversing" for t in tr), PALE2)]
assert [c[0] for c in counts] == D_LABELS and sum(c[1] for c in counts) == PANEL_N, counts
D_XLIM = PANEL_N * 0.68                                             # split_lagladder_polish.py: the count axis runs to 0.68 x the panel size

# ------------------------------------------------------------------ the growth rule: pitch kept, panels grow at the foot, the strings under them and the c/d band move
N_A, N_B = len(ROWS_A), len(ROWS_B)
GROW_A, GROW_B = round((N_A - ROWS13["a"]) * PITCH, 3), round((N_B - ROWS13["b"]) * PITCH, 3); GROW_CD = max(GROW_A, GROW_B, 0.0)
AX = dict(a=[AX13["a"][0], AX13["a"][1], AX13["a"][2], AX13["a"][3] + GROW_A], b=[AX13["b"][0], AX13["b"][1], AX13["b"][2], AX13["b"][3] + GROW_B],
          c=[AX13["c"][0], AX13["c"][1] + GROW_CD, AX13["c"][2], AX13["c"][3] + GROW_CD], d=[AX13["d"][0], AX13["d"][1] + GROW_CD, AX13["d"][2], AX13["d"][3] + GROW_CD])
PAGE_H = PAGE_H13 + GROW_CD
def shift_y(x, y):
    """Region rule for a V13 string origin: top rows unchanged, under panel a +GROW_A, under panel b +GROW_B, the c/d band +GROW_CD."""
    if y >= 700.0: return GROW_CD
    if x < 500.0 and y >= AX13["a"][3] - 1.0: return GROW_A
    if x >= 500.0 and y >= AX13["b"][3] - 1.0: return GROW_B
    return 0.0
print(f"rows a {N_A} (V13 {ROWS13['a']}), rows b {N_B} (V13 {ROWS13['b']}), pitch {PITCH}: growth a {GROW_A}, b {GROW_B}, c/d band {GROW_CD}, page {PAGE_H13} -> {PAGE_H}")
print("rows a:", [LAB[k] for k in ROWS_A]); print("rows b:", [LAB[k] if k else "---" for k in ROWS_B]); print("label delta:", LABEL_DELTA)
print("no estimate:", sorted(noest)); print("panel c surv", surv, "rem", rem, "ret", [round(v, 3) for v in ret]); print("panel d", [(c[0], c[1]) for c in counts], "xlim", D_XLIM)

# ------------------------------------------------------------------ helpers
def xpos(box, hr): return box[0] + (math.log(hr) - math.log(XLO)) / (math.log(XHI) - math.log(XLO)) * (box[2] - box[0])
def add_axes(fig, box):
    x0, y0, x1, y1 = box; return fig.add_axes([x0 / PAGE_W, (PAGE_H - y1) / PAGE_H, (x1 - x0) / PAGE_W, (y1 - y0) / PAGE_H])
def probed_origin(text, xmin, xmax, ymin, ymax, font=None):
    """Origin (x, baseline) of a V13 span (first copy) plus the region shift."""
    c = [s for s in BT["spans"] if s["text"] == text and xmin <= s["origin"][0] <= xmax and ymin <= s["origin"][1] <= ymax and (font is None or s["font"] == font)]
    assert c, (text, xmin, xmax, ymin, ymax)
    o = min(c, key=lambda s: s["origin"][0])["origin"]; return [o[0], o[1] + shift_y(o[0], o[1])]
def unit_marker_path(m): ms = MarkerStyle(m); return ms.get_path().transformed(ms.get_transform())
def marker_patch(fig, px, py, m, size_pt, face, edge, lw, z=5):
    t = Affine2D().scale(size_pt * fig.dpi / 72.0).translate(px, py); p = PathPatch(unit_marker_path(m), transform=t, facecolor=face, edgecolor=edge, linewidth=lw, zorder=z, clip_on=False); fig.add_artist(p); return p
PRINTED = []
def fixed_text(ov, x, y, s, size, kind, weight="normal", style="normal", color=INK, rotation=0, centre_keep=False):
    """Round 49: a string V13 centres (the axis titles) is placed by its V13 LEFT origin, so at +1 pt it would grow to one side and lose its centre.
    centre_keep shifts the origin back by half the width growth (Arial advance widths at the two sizes), along the text direction, so the V13 centre is kept."""
    shift = 0.0
    if centre_keep:
        shift = (_ARIAL.text_length(s, size) - _ARIAL.text_length(s, size - PT_PLUS)) / 2.0
        if rotation == 90: y = y + shift          # the text runs upward from its origin: the origin moves down by half the growth
        else: x = x - shift
    ov.text(x, y, s, fontsize=size, color=color, ha="left", va="baseline", fontweight=weight, fontstyle=style, rotation=rotation, rotation_mode="anchor")
    PRINTED.append(dict(text=s, x=x, baseline=y, size=size, weight=weight, style=style, kind=kind, placement=("V13 centre kept (origin moved %.2f pt for the +1 pt width growth)" % shift) if centre_keep else "V13 origin plus the region shift", source="base sheet wording", centre_shift_r49=round(shift, 3)))

def draw_ladder(fig, ax, rows, box, letter):
    n = len(rows); pitch = (box[3] - box[1]) / n; assert abs(pitch - PITCH) < 0.01, (letter, pitch, PITCH)
    drawn, labels = [], []; shade = True
    for i, key in enumerate(rows):
        if key is None: continue
        if shade: ax.axhspan(n - i - 1.0, n - i, color=GRID, lw=0, zorder=0)
        shade = not shade
    if None in rows:
        ib = rows.index(None)
        for ya, yb in ((0.0, n - ib - 1.0), (n - ib, float(n))): ax.plot([1.0, 1.0], [ya, yb], color=INK, lw=0.9, ls=(0, (4, 3)), zorder=1)
    else: ax.axvline(1.0, color=INK, lw=0.9, ls=(0, (4, 3)), zorder=1)
    ax.set_xscale("log"); ax.set_xticks(XTICKS); ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:g}")); ax.xaxis.set_minor_locator(mticker.NullLocator())
    ax.set_xlim(XLO, XHI); ax.set_ylim(0.0, float(n)); ticks = []
    for i, key in enumerate(rows):
        y0 = n - i - 0.5
        if key is None: continue
        lab = f"{DISPLAY[key]} ({E5[key]:,})"; ticks.append((y0, lab))
        labels.append(dict(text=lab, kind="row label", panel=letter, row=i, x_right=round(box[0] - GUTTER_PAD, 3), y_center=round(box[1] + (i + 0.5) * pitch, 3), placement="matplotlib y tick label, pad 4 pt, centred on the row",
                           source=DATA, keys=[f"per_condition.{key}.condition", f"per_condition.{key}.events_by_lag.5"], value=E5[key], sha256=SHA, condition=LAB[key], display=DISPLAY[key]))
        for k, L in enumerate(LAGS):
            y = y0 + RUNG_OFF[k]; y_page = box[1] + (n - y) * pitch; hr = PC[key]["hr_by_lag"][L]
            if hr is None:
                ax.text(XLO * 1.02, y, "No estimate", fontsize=FLOOR_PT, color=COMPARE, va="center", ha="left", style="italic")
                drawn.append(dict(panel=letter, row=i, key=key, condition=LAB[key], lag=int(L), estimable=False, y_page=round(y_page, 3), x_text_left=round(xpos(box, XLO * 1.02), 3), keys=[f"per_condition.{key}.hr_by_lag.{L}"], sha256=SHA))
                PRINTED.append(dict(text="No estimate", kind="no estimate", panel=letter, row=i, x=round(xpos(box, XLO * 1.02), 3), y_center=round(y_page, 3), size=FLOOR_PT, style="italic", placement="matplotlib text, left at HR 0.6324, centred on the rung",
                                    source=DATA, keys=[f"per_condition.{key}.hr_by_lag.{L}"], value=None, sha256=SHA))
                continue
            lo, hi = PC[key]["ci_by_lag"][L]; quiet = (L == "5") and not bool(PC[key]["lag5_informative"])
            ax.plot([lo, hi], [y, y], color=PRIMARY, lw=CI_LW, solid_capstyle="round", zorder=2)
            px, py = ax.transData.transform((hr, y)); marker_patch(fig, px, py, MARK[k], math.sqrt(SIZE[k]), WHITE if quiet else PRIMARY, PRIMARY if quiet else WHITE, MEDGE)
            drawn.append(dict(panel=letter, row=i, key=key, condition=LAB[key], lag=int(L), estimable=True, hr=hr, lo=lo, hi=hi, quiet=quiet, marker=MARK[k], marker_width_pt=round(math.sqrt(SIZE[k]), 3),
                              x_page=round(xpos(box, hr), 3), x_lo=round(xpos(box, lo), 3), x_hi=round(xpos(box, hi), 3), y_page=round(y_page, 3),
                              keys=[f"per_condition.{key}.hr_by_lag.{L}", f"per_condition.{key}.ci_by_lag.{L}", f"per_condition.{key}.lag5_informative"], sha256=SHA))
    ys, labs = zip(*ticks); ax.set_yticks(ys); ax.set_yticklabels(labs, fontsize=TICK_PT); ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0, pad=GUTTER_PAD); ax.set_xlabel("")
    for t in XTICKS: PRINTED.append(dict(text=f"{t:g}", kind="x tick label", panel=letter, x_center=round(xpos(box, t), 3), size=TICK_PT, placement="matplotlib x tick label", source="axis design (V13 sheet)"))
    return drawn, labels

def key_geometry(which):
    K = G["key_" + which]
    mk = [m for m in K["markers"] if m["fill"] in (PRIMARY, WHITE) and m["edge"] in (PRIMARY, WHITE) and m["h"] > 3.0]; ln = [l for l in K["lines"] if l["colour"] == PRIMARY]
    def collapse(items, keyf, tol=0.7):
        out = []
        for it in sorted(items, key=lambda i: -i["seq"]):
            if not any(abs(keyf(it)[0] - keyf(o)[0]) < tol and abs(keyf(it)[1] - keyf(o)[1]) < tol for o in out): out.append(it)
        return out
    mk = collapse(mk, lambda m: (m["cx"], m["cy"])); ln = collapse(ln, lambda l: (l["x0"], l["y"]))
    return sorted(mk, key=lambda m: (round(m["cx"]), m["cy"])), sorted(ln, key=lambda l: (round(l["x0"]), l["y"]))

# ------------------------------------------------------------------ draw
with plt.rc_context(RC):
    fig = plt.figure(figsize=(PAGE_W / 72.0, PAGE_H / 72.0)); fig.patch.set_alpha(0.0)
    out = dict(sheet=OUT_PDF, page=[PAGE_W, PAGE_H], page_v13=[PAGE_W, PAGE_H13], growth=dict(a=GROW_A, b=GROW_B, cd=GROW_CD, pitch=PITCH, rows_a=N_A, rows_b=N_B, rows_a_v13=ROWS13["a"], rows_b_v13=ROWS13["b"]),
               source=SRC, mode="OLD snapshot (positive control)" if OLDMODE else "v8.1", pt_plus=PT_PLUS, text_pt=dict(tick=TICK_PT, title=TITLE_PT, annotation=ANN_PT, floor=FLOOR_PT, letter=LETTER_PT), panels={}, panel_n=PANEL_N, half=half, label_delta=LABEL_DELTA,
               xlim=[XLO, XHI], xlim_v13=[XLO13, XHI13], display_labels=LABEL_NOTES, row_order_rule="efig14_polish.py: pandas sort_values('hr', ascending=False) on the lag-0 rows of lag_ladder.csv, conditions with lo > 1, not a control, not circular; ceil(n/2) in a; controls the same way", lag0_ties=TIES, n_lag5_informative=n_informative)
    for letter, rows in (("a", ROWS_A), ("b", ROWS_B)):
        ax = add_axes(fig, AX[letter]); drawn, labels = draw_ladder(fig, ax, rows, AX[letter], letter)
        out["panels"][letter] = dict(axes_box=AX[letter], axes_box_v13=AX13[letter], rows=[LAB[k] if k else CTRL_HEAD for k in rows], row_keys=[k for k in rows], row_pitch=round((AX[letter][3] - AX[letter][1]) / len(rows), 4),
                                     xlim=[XLO, XHI], xticks=XTICKS, rung_offsets=RUNG_OFF, rungs=drawn, labels=labels)
    axc = add_axes(fig, AX["c"]); xc = [0, 1, 2, 3]
    axc.bar(xc, rem, 0.50, color=COMPARE, linewidth=0, zorder=2); axc.plot(xc, ret, color=PRIMARY, lw=1.8, zorder=3)
    axc.set_xticks(xc); axc.set_xticklabels(["0", "1", "2", "5"]); axc.set_ylim(0, 142); axc.set_yticks([0, 25, 50, 75, 100]); axc.set_xlim(-0.6, 3.6)
    c_ann = []
    for i in range(4):
        yv = max(ret[i], rem[i]) + 6.0; s = f"{surv[i]}/{PANEL_N}"; axc.text(xc[i], yv, s, ha="center", fontsize=ANN_PT, color=PRIMARY)
        px, py = axc.transData.transform((xc[i], ret[i])); marker_patch(fig, px, py, "o", MSIZE[0], PRIMARY, WHITE, MEDGE)
        xpage = AX["c"][0] + (xc[i] + 0.6) / 4.2 * (AX["c"][2] - AX["c"][0]); ypage = lambda v: AX["c"][3] - v / 142.0 * (AX["c"][3] - AX["c"][1])
        c_ann.append(dict(text=s, kind="annotation", panel="c", landmark=[0, 1, 2, 5][i], x_center=round(xpage, 3), y_baseline=round(ypage(yv), 3), size=ANN_PT, color=PRIMARY, placement="matplotlib text centred on the landmark, baseline 6 units above max(line, bar)",
                          source=DATA, keys=["panel_n"] if i == 0 else [f"by_lag.{[0, 1, 2, 5][i]}.survive", "panel_n"], value=surv[i], sha256=SHA, line_value_pct=round(ret[i], 3), bar_value_pct=rem[i], x_marker=round(xpage, 3), y_marker=round(ypage(ret[i]), 3),
                          bar_rect=[round(xpage - 0.25 / 4.2 * (AX["c"][2] - AX["c"][0]), 3), round(ypage(rem[i]), 3), round(xpage + 0.25 / 4.2 * (AX["c"][2] - AX["c"][0]), 3), AX["c"][3]],
                          value_keys=[] if i == 0 else [f"by_lag.{[0, 1, 2, 5][i]}.median_pct_events_removed", f"by_lag.{[0, 1, 2, 5][i]}.median_excess_retained"]))
    PRINTED.extend(c_ann)
    for t in ["0", "1", "2", "5"]: PRINTED.append(dict(text=t, kind="x tick label", panel="c", size=TICK_PT, placement="matplotlib", source="axis design (V13 sheet)"))
    for t in ["0", "25", "50", "75", "100"]: PRINTED.append(dict(text=t, kind="y tick label", panel="c", size=TICK_PT, placement="matplotlib", source="axis design (V13 sheet)"))
    out["panels"]["c"] = dict(axes_box=AX["c"], axes_box_v13=AX13["c"], ylim=[0, 142], xlim=[-0.6, 3.6], landmarks=[0, 1, 2, 5], events_deleted_pct=rem, excess_retained_pct=[round(v, 3) for v in ret], survive=surv, panel_n=PANEL_N, annotations=c_ann)
    axd = add_axes(fig, AX["d"]); yd = [3, 2, 1, 0]
    axd.barh(yd, [c[1] for c in counts], 0.55, color=[c[2] for c in counts], linewidth=0, zorder=2); axd.set_ylim(-0.4525, 3.4525); axd.set_yticks([]); axd.spines["left"].set_visible(False)
    axd.set_xlim(0, D_XLIM); axd.set_xticks([t for t in [0, 5, 10, 15, 20, 25, 30] if t <= D_XLIM]); d_counts = []
    for i, (lab, nn, col) in enumerate(counts):
        axd.text(nn + 0.5, yd[i], str(nn), va="center", fontsize=ANN_PT, color=INK)
        xr = AX["d"][0] + nn / D_XLIM * (AX["d"][2] - AX["d"][0]); yc = AX["d"][1] + (3.4525 - yd[i]) / 3.905 * (AX["d"][3] - AX["d"][1])
        d_counts.append(dict(text=str(nn), kind="count", panel="d", label=lab, x_left=round(AX["d"][0] + (nn + 0.5) / D_XLIM * (AX["d"][2] - AX["d"][0]), 3), y_center=round(yc, 3), size=ANN_PT, placement="matplotlib text 0.5 units right of the bar end, centred on the bar",
                             source=DATA, keys=[f"trends.{k}.trend" for k in panel], rule=("trend startswith 'stable'" if lab == "Stable" else f"trend == '{lab.lower()}'"), value=nn, sha256=SHA,
                             bar_rect=[AX["d"][0], round(yc - 0.275 / 3.905 * (AX["d"][3] - AX["d"][1]), 3), round(xr, 3), round(yc + 0.275 / 3.905 * (AX["d"][3] - AX["d"][1]), 3)], colour=col))
    PRINTED.extend(d_counts)
    for t in [str(t) for t in [0, 5, 10, 15, 20, 25, 30] if t <= D_XLIM]: PRINTED.append(dict(text=t, kind="x tick label", panel="d", size=TICK_PT, placement="matplotlib", source="axis design (V13 sheet)"))
    out["panels"]["d"] = dict(axes_box=AX["d"], axes_box_v13=AX13["d"], xlim=[0, D_XLIM], ylim=[-0.4525, 3.4525], counts={c[0]: c[1] for c in counts}, colours={c[0]: c[2] for c in counts}, trend_by_condition={LAB[k]: J["trends"][k]["trend"] for k in panel}, count_texts=d_counts)
    ov = fig.add_axes([0, 0, 1, 1]); ov.set_xlim(0, PAGE_W); ov.set_ylim(PAGE_H, 0); ov.axis("off"); ov.patch.set_visible(False)
    xa = probed_origin(XLABEL, 100, 300, 480, 520); xb = probed_origin(XLABEL, 600, 800, 640, 700)
    fixed_text(ov, xa[0], xa[1], XLABEL, TITLE_PT, "axis title a", centre_keep=True); fixed_text(ov, xb[0], xb[1], XLABEL, TITLE_PT, "axis title b", centre_keep=True)
    hd13 = [s for s in BT["spans"] if s["text"] == CTRL_HEAD and s["font"] == "Arial-BoldMT"]; assert len(hd13) == 1, hd13
    ib13 = ROWS13["b"] - 5 - 1; hd_off = hd13[0]["origin"][1] - (AX13["b"][1] + (ib13 + 0.5) * PITCH)          # the header's baseline offset from the centre of its (blank) row on V13
    ib = ROWS_B.index(None); fixed_text(ov, hd13[0]["origin"][0], AX["b"][1] + (ib + 0.5) * PITCH + hd_off, CTRL_HEAD, TICK_PT, "group header b", weight="bold")
    cx_ = probed_origin(C_XLABEL, 100, 400, 900, 960); fixed_text(ov, cx_[0], cx_[1], C_XLABEL, TITLE_PT, "axis title c x", centre_keep=True)
    cy_ = probed_origin(C_YLABEL, 30, 80, 800, 900); fixed_text(ov, cy_[0], cy_[1], C_YLABEL, TITLE_PT, "axis title c y", rotation=90, centre_keep=True)
    dx_ = probed_origin(D_XLABEL, 700, 900, 900, 960); fixed_text(ov, dx_[0], dx_[1], D_XLABEL, TITLE_PT, "axis title d", centre_keep=True)
    for lab in D_LABELS:
        o = probed_origin(lab, 500, 600, 740, 900); fixed_text(ov, o[0], o[1], lab, TICK_PT, "row label d")
    # round 40 (LNOTES, 2026-09-18): the two printed key lines of panel c were removed from the sheet (the wording lives in the legend); round 49 rebuilds without them
    REMOVED_R40 = [KEY_LINE, KEY_BARS]; out["removed_r40_notes"] = REMOVED_R40
    def draw_key(which, texts, dy):
        mk, ln = key_geometry(which); assert len(mk) == 5 and len(ln) == 4, (which, len(mk), len(ln))
        order = [("o", MSIZE[0], False), ("s", MSIZE[1], False), ("D", MSIZE[2], False), ("^", MSIZE[3], False), ("^", MSIZE[3], True)]; exp_w = [6.0, 5.4, 8.062, 6.3, 6.3]; rec = []
        for m, (marker, ms, opn), w in zip(mk, order, exp_w):
            assert abs(m["w"] - w) < 0.05 and (m["fill"] == WHITE) == opn, (which, m, marker, w, opn)
            px, py = ov.transData.transform((m["cx"], m["cy"] + dy)); marker_patch(fig, px, py, marker, ms, WHITE if opn else PRIMARY, PRIMARY if opn else WHITE, MEDGE); rec.append(dict(x=m["cx"], y=m["cy"] + dy, marker=marker, ms=ms, open=opn))
        for l in ln: ov.plot([l["x0"], l["x1"]], [l["y"] + dy, l["y"] + dy], color=PRIMARY, lw=CI_LW, solid_capstyle="projecting", zorder=2)
        for t in texts: fixed_text(ov, t[1][0], t[1][1], t[0], TICK_PT, f"key text {which}")
        return dict(markers=rec, lines=[dict(x0=l["x0"], x1=l["x1"], y=l["y"] + dy) for l in ln], texts=[dict(text=t[0], x=t[1][0], baseline=t[1][1]) for t in texts], shift=dy)
    ta = [(LAGNAME[0], probed_origin(LAGNAME[0], 40, 100, 505, 545)), (LAGNAME[1], probed_origin(LAGNAME[1], 40, 100, 505, 545)), (LAGNAME[2], probed_origin(LAGNAME[2], 150, 200, 505, 545)),
          (LAGNAME[3], probed_origin(LAGNAME[3], 150, 200, 505, 545)), (QUIET_1, probed_origin(QUIET_1, 250, 320, 505, 545))]
    tb = [(LAGNAME[i], probed_origin(LAGNAME[i], 880, 940, 40, 140)) for i in range(4)] + [(QUIET_2[0], probed_origin(QUIET_2[0], 880, 940, 40, 140)), (QUIET_2[1], probed_origin(QUIET_2[1], 880, 940, 40, 140))]
    out["panels"]["a"]["key"] = draw_key("a", ta, GROW_A); out["panels"]["b"]["key"] = draw_key("b", tb, 0.0)
    letters = []
    for L in G["letters_text_layer"]:
        if not L["visible"]: continue
        if any(abs(L["x"] - o["x"]) < 0.5 and abs(L["baseline"] - o["baseline_v13"]) < 0.5 for o in letters): continue
        dy = GROW_CD if L["text"] in ("c", "d") else 0.0; letters.append(dict(text=L["text"], x=L["x"], baseline=L["baseline"] + dy, baseline_v13=L["baseline"], shift=dy))
    assert [l["text"] for l in sorted(letters, key=lambda l: (round(l["baseline"]), l["x"]))] == ["a", "b", "c", "d"], letters
    for l in letters: fixed_text(ov, l["x"], l["baseline"], l["text"], LETTER_PT, "letter", weight="bold")
    out["letters"] = letters
    fig.savefig(OUT_PDF, format="pdf"); plt.close(fig)
out["printed"] = PRINTED
json.dump(out, open(OUT_JSON, "w"), indent=1, default=float)
na = sum(1 for d in out["panels"]["a"]["rungs"] if d["estimable"]); nb = sum(1 for d in out["panels"]["b"]["rungs"] if d["estimable"])
print(f"wrote {OUT_PDF} ({os.path.getsize(OUT_PDF)} bytes) and {OUT_JSON}; markers a {na}, no estimate a {sum(1 for d in out['panels']['a']['rungs'] if not d['estimable'])}, markers b {nb}, no estimate b {sum(1 for d in out['panels']['b']['rungs'] if not d['estimable'])}, printed strings {len(PRINTED)}, letters {[(l['text'], l['baseline']) for l in letters]}")
