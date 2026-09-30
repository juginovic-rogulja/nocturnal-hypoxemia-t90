#!/usr/bin/env python3
"""V14 lane L4, round v8.2 (2026-09-15 night), Supp_Fig11: wake-period and sleep-period T90 in the same Cox model in the disease-free
stratum (C_organ_metab_free), REBUILT from the re-extracted step-121 files (numbers/wake_sleep_healthy_C_organ_metab_free.csv and
wake_sleep_healthy.json, sidecars 19:27 today) with the v8.1 control panel. Until now the sheet was the V13 copy (group P).
Data wiring = FINAL_FIGURES_2026-08-14/_scripts/make_eFigureNEW_wake_sleep_diseasefree.py (every row of the stratum, ordered by the
asleep-minus-awake gap, the controls the file carries after them, ordered the same; FDR q < 0.05 fills the marker). That generator
cannot run tonight (its splitstyle import asserts the pre-14-September control list and prism_plotter is not on the path), so the
sheet is re-plotted here in page points pinned to the V13 sheet's text layer (label x and baselines, sub-row offsets, column edges,
key block, tick calibration, rule y) with the generator's marks (wake grey circles 6 pt, sleep blue diamonds 5.5 pt, filled with a
0.8 pt white edge when q < 0.05, open with the series ring otherwise, 1.7 pt round-capped intervals, 0.9 pt dashed reference rule,
0.8 pt spine, 3 pt ticks), the V13 design of round 31 (no grid lines, no raster key). The page height follows the row count at the
V13 pitch (the generator's own rule). Fonts Arial (the V13 key strings were Helvetica: declared).
Outputs: Supp_Fig11/Supp_Fig11.pdf, _150dpi.png, work/{raw.pdf, geometry.json}, verify/{Supp_Fig11_drawn.json, checks.txt,
CHANGES_Supp_Fig11.csv, Supp_Fig11_printed_values.csv, crops/}.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import collections, csv, json, math, os, re, sys, time
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import l4lib as L
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND37_2026-09-14/figures/_common"); import v14lib as V
import fitz
from PIL import Image

SHEET = "Supp_Fig11"; SD = f"{L.LANE}/{SHEET}"; WORK = f"{SD}/work"; VER = f"{SD}/verify"; CROPS = f"{VER}/crops"
for d_ in (WORK, VER, CROPS): os.makedirs(d_, exist_ok=True)
BASE = L.hydrated(f"{L.BASE}/{SHEET}.pdf"); OUT = f"{SD}/{SHEET}.pdf"; PNG = f"{SD}/{SHEET}_150dpi.png"
matplotlib, plt = L.mpl_setup()
from matplotlib.colors import to_rgb
V8_2_CUT = "2026-09-15 19:18:00"

# ------------------------------------------------------------------ files, gated (v8.2)
C_PATH = f"{L.NUM}/wake_sleep_healthy_C_organ_metab_free.csv"; J_PATH = f"{L.NUM}/wake_sleep_healthy.json"; NC_PATH = f"{L.NUM}/negcontrols_final.json"
GATES = {}
for p in (C_PATH, J_PATH):
    g = V.sidecar_v8(p); assert g["output_mtime"] >= V8_2_CUT, ("not the v8.2 re-extraction", p, g["output_mtime"]); GATES[p] = g
GATES[NC_PATH] = V.sidecar_v8(NC_PATH)
C = pd.read_csv(L.hydrated(C_PATH)); META = json.load(open(L.hydrated(J_PATH))); NCF = json.load(open(L.hydrated(NC_PATH)))
NEG_CONTROLS = list(NCF["panel_labels"]); assert len(NEG_CONTROLS) == 5
assert set(C.cond[C.neg]) <= set(NEG_CONTROLS), sorted(C.cond[C.neg]); assert not C.cond.duplicated().any()
N_C = int(META["strata"]["C_organ_metab_free"]["n"]); QTHR = float(META["provenance"]["q_threshold"]); assert QTHR == 0.05
d = C.copy(); d["gap"] = d.sHR - d.wHR
main = d[~d.neg].sort_values("gap", ascending=False); ctrl = d[d.neg].sort_values("gap", ascending=False)
ROWS = list(main.itertuples()) + list(ctrl.itertuples()); N_MAIN, N_CTRL = len(main), len(ctrl)
CTRL_ABSENT = [c for c in NEG_CONTROLS if c not in set(C.cond)]
N_SIG_W, N_SIG_S = int((C.wq < QTHR).sum()), int((C.sq < QTHR).sum()); N_SIG_S_UP = int(((C.sq < QTHR) & (C.sHR > 1)).sum())

# ------------------------------------------------------------------ V13 pins
IN = L.spans_of(BASE); SI = L.dedupe(IN); W13, H13 = fitz.open(BASE)[0].rect.width, fitz.open(BASE)[0].rect.height
labs13 = sorted([s for s in SI if s["size"] == 10.0 and s["font"] == "ArialMT" and s["origin"][0] < 20], key=lambda s: s["origin"][1])
hr13 = [s for s in SI if re.fullmatch(r"\d+\.\d\d \(\d+\.\d\d-\d+\.\d\d\)", s["text"])]; q13 = [s for s in SI if re.fullmatch(r"<0\.001|0\.\d{3}", s["text"])]
head13 = next(s for s in SI if s["text"] == "Negative controls")
assert len(hr13) == 2 * len(labs13) == len(q13), (len(labs13), len(hr13), len(q13))
N13 = len(labs13); Y13 = [s["origin"][1] for s in labs13]; PITCH = float(np.median(np.diff(Y13)[:15])); Y_FIRST = Y13[0]; LAB_X = float(np.median([s["origin"][0] for s in labs13]))
first_hr = sorted([s for s in hr13 if abs(s["origin"][1] - Y_FIRST) < 15], key=lambda s: s["origin"][1]); OFF_W, OFF_S = first_hr[0]["origin"][1] - Y_FIRST, first_hr[1]["origin"][1] - Y_FIRST
COL_R = float(np.median([s["bbox"][2] for s in hr13])); QCOL_R = float(np.median([s["bbox"][2] for s in q13]))
HR_HEAD = next(s for s in SI if s["text"] == "HR (95% CI)"); Q_HEAD = next(s for s in SI if s["text"] == "q" and s["origin"][1] < 40)
KEY = {t: next(s for s in SI if s["text"] == t) for t in ("Sleep-period T90", "Wake-period T90", "Not significant")}
CAP = next(s for s in SI if s["text"].startswith("Hazard ratio (95% CI) per 1 SD"))
ticks13 = {s["text"]: s for s in SI if s["text"] in ("0.6", "0.8", "1", "1.5", "2") and s["origin"][1] > 600}; TICK_BASE = float(np.median([s["origin"][1] for s in ticks13.values()]))
tc = {float(k): (v["bbox"][0] + v["bbox"][2]) / 2 for k, v in ticks13.items()}
B_LN = (tc[2.0] - tc[0.6]) / (math.log(2.0) - math.log(0.6)); X1 = tc[1.0]
def xv(v): return X1 + B_LN * math.log(v)
# the axis rule from the V13 render (150 dpi): the long horizontal ink line above the tick labels
im13 = np.asarray(Image.open(L.hydrated(f"{paths.FIGURE_ROOT}/ROUND36_2026-09-10/L9_assembly/word_png_v13/{SHEET}.png")).convert("RGB")).astype(int); s150 = im13.shape[1] / W13
band = im13[int((TICK_BASE - 22) * s150):int((TICK_BASE - 4) * s150), :]; dark = band.sum(axis=2) < 300; rr = np.where(dark.sum(axis=1) > 150)[0]
RULE_Y = float(np.mean(rr) / s150 + TICK_BASE - 22); xs = np.where(dark[rr[0]])[0]; AX_X0, AX_X1 = float(xs.min() / s150), float(xs.max() / s150)
XLO, XHI = math.exp((AX_X0 - X1) / B_LN), math.exp((AX_X1 - X1) / B_LN)
# the block header sits one pitch under the last row's label, the controls start one pitch under the header, the axis a fixed gap under the last control
HEAD_GAP = head13["origin"][1] - Y13[N13 - 1 - 2]     # V13: 2 controls after the header (the V13 sheet's own count, read from its text layer)
n_ctrl13 = sum(1 for s in labs13 if s["origin"][1] > head13["origin"][1]); assert n_ctrl13 == 2, n_ctrl13
HEAD_GAP = head13["origin"][1] - Y13[N13 - n_ctrl13 - 1]; CTRL_GAP = Y13[N13 - n_ctrl13] - head13["origin"][1]; RULE_GAP = RULE_Y - Y13[-1]
CAP_GAP = CAP["origin"][1] - RULE_Y; TICK_GAP = TICK_BASE - RULE_Y; BOTTOM_GAP = H13 - CAP["origin"][1]
# new layout: the same pins, rows counted from the file
y_main = [Y_FIRST + i * PITCH for i in range(N_MAIN)]; y_head = y_main[-1] + HEAD_GAP; y_ctrl = [y_head + CTRL_GAP + i * PITCH for i in range(N_CTRL)]
Y_LAST = y_ctrl[-1] if y_ctrl else y_main[-1]; RULE_NEW = Y_LAST + RULE_GAP; TICK_NEW = RULE_NEW + TICK_GAP; CAP_NEW = RULE_NEW + CAP_GAP; H_NEW = CAP_NEW + BOTTOM_GAP
GREY, BLUE, INK = "#8a9099", L.BLUE, L.INK; MS = {"w": 6.0, "s": 5.5}; MK = {"w": "o", "s": "D"}; COL = {"w": GREY, "s": BLUE}
FS_ROW, FS_COL, FS_KEY, FS_CAP, FS_TICK = 10.0, 9.5, 9.2, 11.0, 10.0
TICKS = [t for t in (0.5, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0, 4.0) if XLO < t < XHI]
vals = [float(getattr(r, p + k)) for r in ROWS for p in ("w", "s") for k in ("lo", "hi", "HR")]; assert XLO < min(vals) and max(vals) < XHI, (XLO, min(vals), max(vals), XHI)
# key handles (V13 render): a 17 pt line ending 5.3 pt before the label, the marker at its centre; the open row: diamond then circle before the label
KEY_NUDGE = -11.0   # round 49: the key block (two handles, three markers, three strings) moved 11 pt left as one element, so that the key strings at 10.2 pt
                    # stay inside the 8 pt page margin (at the V26 position the two T90 strings, 87 pt wide, would end 2.4 pt beyond the page edge)
KEY_X = KEY["Sleep-period T90"]["origin"][0] + KEY_NUDGE; HANDLE = (KEY_X - 5.3 - 17.0, KEY_X - 5.3); KEY_MID = {t: KEY[t]["origin"][1] - 3.3 for t in KEY}
OPEN_D_X, OPEN_O_X = KEY_X - 5.3 - 17.0 + 3.0, KEY_X - 5.3 - 3.0

def hrci(r, p): return L.hr_ci(getattr(r, p + "HR"), getattr(r, p + "lo"), getattr(r, p + "hi"))
STR = {(r.cond, p): (hrci(r, p), L.q_text(float(getattr(r, p + "q"))), bool(float(getattr(r, p + "q")) < QTHR)) for r in ROWS for p in ("w", "s")}
for s_ in ("Sleep-period T90", "Wake-period T90", "Not significant", "HR (95% CI)", "q", "Negative controls", CAP["text"], *[r.cond for r in ROWS]): assert L.house_ok(s_), s_

# ------------------------------------------------------------------ draw (page points, y downward)
RECS, DRAWN = [], []
def T(x, y, s, size, ha="left", bold=False, **src):
    size = size + L.PT_PLUS   # round 49: every text one point larger (10 -> 11, 9.5 -> 10.5, 9.2 -> 10.2, 11 -> 12)
    ov.text(x, y, s, fontsize=size, color=INK, ha=ha, va="baseline", fontweight="bold" if bold else "normal")
    RECS.append(dict(text=s, x=round(x, 3), baseline=round(y, 3), ha=ha, size=size, panel="", **src))
RC = dict(L.RC); RC.update({"font.size": 10.0, "figure.facecolor": "none", "savefig.facecolor": "none"})
with plt.rc_context(RC):
    fig = plt.figure(figsize=(W13 / 72.0, H_NEW / 72.0)); fig.patch.set_alpha(0.0)
    ov = fig.add_axes([0, 0, 1, 1]); ov.set_xlim(0, W13); ov.set_ylim(H_NEW, 0); ov.axis("off"); ov.patch.set_visible(False)
    ov.plot([X1, X1], [Y_FIRST + OFF_W - 12.0, RULE_NEW], color=INK, lw=0.9, ls=(0, (4, 3)), zorder=1, solid_capstyle="butt")
    ov.plot([AX_X0, AX_X1], [RULE_NEW, RULE_NEW], color=INK, lw=0.8, zorder=4, solid_capstyle="projecting")
    for t in TICKS:
        ov.plot([xv(t), xv(t)], [RULE_NEW, RULE_NEW + 3.0], color=INK, lw=0.8, zorder=4, solid_capstyle="butt")
        T(xv(t), TICK_NEW, f"{t:g}", FS_TICK, ha="center", source_file="static:V13", source_key="axis tick (ticks inside the V13 axis range)", source_value=t, rule="tick")
    for r, y_lab in zip(ROWS, y_main + y_ctrl):
        T(LAB_X, y_lab, r.cond, FS_ROW, source_file=C_PATH, source_key=f"cond={r.cond}:cond", source_value=r.cond, rule="text")
        for p, off in (("w", OFF_W), ("s", OFF_S)):
            yb = y_lab + off; yc = yb - 3.3; hr, lo, hi, q = (float(getattr(r, p + k)) for k in ("HR", "lo", "hi", "q")); sigp = q < QTHR; col = COL[p]
            ov.plot([xv(lo), xv(hi)], [yc, yc], color=col, lw=1.7, solid_capstyle="round", zorder=3)
            ov.plot([xv(hr)], [yc], marker=MK[p], ms=MS[p], ls="None", color=col, markerfacecolor=col if sigp else "white", markeredgecolor="white" if sigp else col, markeredgewidth=0.8 if sigp else 1.0, zorder=4)
            hs, qs, _ = STR[(r.cond, p)]
            T(COL_R, yb, hs, FS_COL, ha="right", source_file=C_PATH, source_key=f"cond={r.cond}:{p}HR,{p}lo,{p}hi", source_value=[hr, lo, hi], rule="hr_ci_2dp")
            T(QCOL_R, yb, qs, FS_COL, ha="right", bold=sigp, source_file=C_PATH, source_key=f"cond={r.cond}:{p}q", source_value=q, rule="q3")
            DRAWN.append(dict(condition=r.cond, period="wake" if p == "w" else "sleep", hr=hr, lo=lo, hi=hi, q=q, sig=sigp, control=bool(r.neg), x_page=round(xv(hr), 3), y_page=round(yc, 3), hr_ci_text=hs, q_text=qs, gap=float(r.gap)))
    T(LAB_X, y_head, "Negative controls", FS_ROW, bold=True, source_file="static:V13", source_key="V13 block header", source_value="Negative controls", rule="text")
    T(HR_HEAD["bbox"][2], HR_HEAD["origin"][1], "HR (95% CI)", FS_COL, ha="right", source_file="static:V13", source_key="V13 column header", source_value="HR (95% CI)", rule="text")
    T(Q_HEAD["bbox"][2], HR_HEAD["origin"][1], "q", FS_COL, ha="right", source_file="static:V13", source_key="V13 column header", source_value="q", rule="text")
    for t_, p in (("Sleep-period T90", "s"), ("Wake-period T90", "w")):
        ov.plot(list(HANDLE), [KEY_MID[t_]] * 2, color=COL[p], lw=1.7, solid_capstyle="round", zorder=3)
        ov.plot([sum(HANDLE) / 2], [KEY_MID[t_]], marker=MK[p], ms=MS[p], ls="None", color=COL[p], markerfacecolor=COL[p], markeredgecolor="white", markeredgewidth=0.8, zorder=4)
        T(KEY_X, KEY[t_]["origin"][1], t_, FS_KEY, source_file="static:V13", source_key="V13 key", source_value=t_, rule="text")
    ov.plot([OPEN_D_X], [KEY_MID["Not significant"]], marker="D", ms=MS["s"], ls="None", color=BLUE, markerfacecolor="white", markeredgecolor=BLUE, markeredgewidth=1.0, zorder=4)
    ov.plot([OPEN_O_X], [KEY_MID["Not significant"]], marker="o", ms=MS["w"], ls="None", color=GREY, markerfacecolor="white", markeredgecolor=GREY, markeredgewidth=1.0, zorder=4)
    T(KEY_X, KEY["Not significant"]["origin"][1], "Not significant", FS_KEY, source_file="static:V13", source_key="V13 key (open markers)", source_value="Not significant", rule="text")
    T((AX_X0 + AX_X1) / 2, CAP_NEW, CAP["text"], FS_CAP, ha="center", source_file="static:V13", source_key="V13 axis caption", source_value=CAP["text"], rule="text")
    used = {tuple(np.round(to_rgb(ln.get_color()), 4)) for ln in ov.lines} | {tuple(np.round(to_rgb(ln.get_markerfacecolor()), 4)) for ln in ov.lines if ln.get_marker() not in (None, "None", "")} | {tuple(np.round(to_rgb(ln.get_markeredgecolor()), 4)) for ln in ov.lines if ln.get_marker() not in (None, "None", "")}
    assert used <= {tuple(np.round(to_rgb(c), 4)) for c in (GREY, BLUE, INK, "white")}, used
    rend = fig.canvas.get_renderer()
    for t in ov.texts:
        assert t.get_fontsize() >= 9.0 and L.house_ok(t.get_text()); bb = t.get_window_extent(renderer=rend)
        assert bb.x0 / fig.dpi * 72 >= 8 and bb.x1 / fig.dpi * 72 <= W13 - 8, ("text leaves the page", t.get_text())
    raw = f"{WORK}/{SHEET}_raw.pdf"; fig.savefig(raw, transparent=True); plt.close(fig)
met = L.sheet_font_metrics(L.spans_of(BASE) and [dict(s, asc=s.get("asc", 0.905), desc=s.get("desc", -0.212)) for s in L.spans_of(BASE)]) if False else None
# font metrics: the V13 sheet's own Arial descriptors (asc/desc from its text dict)
tdict = fitz.open(BASE)[0].get_text("dict"); met = {}
for b in tdict["blocks"]:
    for l in b.get("lines", []):
        for s in l["spans"]: met.setdefault(s["font"], (int(round(s["ascender"] * 1000)), int(round(s["descender"] * 1000))))
L.align_font_metrics(raw, f"{WORK}/{SHEET}_m.pdf", met)
# flat single page: the composed sheet IS the panel (white page background added so the sheet prints like V13)
doc = fitz.open(); page = doc.new_page(width=W13, height=H_NEW); pdoc = fitz.open(f"{WORK}/{SHEET}_m.pdf")
page.draw_rect(page.rect, color=None, fill=(1, 1, 1), overlay=False); page.show_pdf_page(page.rect, pdoc, 0); pdoc.close()
doc.save(OUT, garbage=4, deflate=True); doc.close()

# ------------------------------------------------------------------ verify
CK = L.Checks(f"{SHEET} V14 v8.2 (lane V14_L4_APNEA, {time.strftime('%Y-%m-%d %H:%M')}): REBUILT from the re-extracted step-121 files (sidecars {GATES[C_PATH]['output_mtime']} / {GATES[J_PATH]['output_mtime']}); V13 base {BASE} sha {L.sha256(BASE)[:16]}; output sha {L.sha256(OUT)[:16]}. Design = V13 pins, page height follows the row count at the V13 pitch.")
SP = L.dedupe(L.spans_of(OUT)); dO = fitz.open(OUT); pO = dO[0]; W, H = pO.rect.width, pO.rect.height
CK.check(f"one page, V13 width {W13:.2f}; height {H13:.2f} -> {H:.2f} (DECLARED: {N13} -> {N_MAIN + N_CTRL} rows at the V13 pitch {PITCH:.2f} pt)", len(dO) == 1 and abs(W - W13) < 0.01 and abs(H - (H13 - (N13 - N_MAIN - N_CTRL) * PITCH)) < 0.5, f"{W:.2f} x {H:.2f}")
fonts_o = {f[3].split("+")[-1] for f in pO.get_fonts(full=True)}; fonts_i = {f[3].split("+")[-1] for f in fitz.open(BASE)[0].get_fonts(full=True)}
CK.check("fonts within the V13 set (Arial faces; the V13 key strings were Helvetica: DECLARED)", fonts_o <= fonts_i, f"out {sorted(fonts_o)}, lost {sorted(fonts_i - fonts_o)}")
def present(rec, tol=0.35):
    for s in SP:
        if s["text"].replace("\xa0", " ") != rec["text"] or abs(s["origin"][1] - rec["baseline"]) > tol: continue
        sx = s["origin"][0] if rec["ha"] == "left" else (s["bbox"][2] if rec["ha"] == "right" else (s["bbox"][0] + s["bbox"][2]) / 2)
        if abs(sx - rec["x"]) <= 0.6: return True
    return False
miss = [r["text"] for r in RECS if not present(r)]; CK.check(f"every drawn string ({len(RECS)}) read back at its origin (0.35 pt)", not miss, f"{miss[:6]}")
extra = [s["text"] for s in SP if s["text"].strip() and not any(s["text"].replace("\xa0", " ") == r["text"] and abs(s["origin"][1] - r["baseline"]) < 0.35 for r in RECS)]
CK.check("no string on the sheet beyond the drawn records", not extra, f"{extra[:6]}")
CROWS = {r["cond"]: r for r in csv.DictReader(open(C_PATH))}
bad = [(d_["condition"], d_["period"]) for d_ in DRAWN if d_["hr_ci_text"] != L.hr_ci(CROWS[d_["condition"]][d_["period"][0] + "HR"], CROWS[d_["condition"]][d_["period"][0] + "lo"], CROWS[d_["condition"]][d_["period"][0] + "hi"]) or d_["q_text"] != L.q_text(float(CROWS[d_["condition"]][d_["period"][0] + "q"])) or d_["sig"] != (float(CROWS[d_["condition"]][d_["period"][0] + "q"]) < QTHR)]
CK.check(f"all {len(DRAWN)} HR (95% CI) and q strings re-derived from the csv (half up 2 dp, q 3 dp, bold when q < 0.05) equal the drawn record", not bad, f"{bad[:4]}")
CK.check("row set = every condition of the stratum plus the controls the file carries, each block ordered by the asleep-minus-awake gap", [r.cond for r in ROWS] == list(main.cond) + list(ctrl.cond) and set(main.cond) == set(C.cond[~C.neg]) and all(np.diff(list(main.gap)) <= 1e-12) and all(np.diff(list(ctrl.gap)) <= 1e-12), f"{N_MAIN} + {N_CTRL} rows, controls absent from the file {CTRL_ABSENT}")
# marks read back from the flat sheet (get_drawings allowed: 0 xobjects, small)
drs = pO.get_drawings(); grid = [d_ for d_ in drs if d_.get("color") and L.CM.rgb_to_hex(d_["color"]) == "#eef0f1"]
CK.check("no pale #eef0f1 grid strokes (V13 design, round 31)", not grid, f"{len(grid)}")
cols = collections.Counter(L.CM.rgb_to_hex(d_[k]) for d_ in drs for k in ("fill", "color") if d_.get(k) is not None)
CK.check("colours within the palette (ink, blue #0288d1, grey #8a9099, white)", set(cols) <= {INK, BLUE, GREY, "#ffffff"}, f"{dict(cols)}")
ci = [d_ for d_ in drs if d_.get("color") and abs(d_.get("width", 0) - 1.7) < 0.05 and len(d_["items"]) == 1 and d_["items"][0][0] == "l"]
got_ci = sorted((round(it["items"][0][1].x, 2), round(it["items"][0][2].x, 2), round(it["items"][0][1].y, 2)) for it in ci)
exp_ci = sorted((round(min(xv(d_["lo"]), xv(d_["hi"])), 2), round(max(xv(d_["lo"]), xv(d_["hi"])), 2), round(d_["y_page"], 2)) for d_ in DRAWN)
CK.check(f"the {len(DRAWN)} intervals read back from the sheet's vectors at the fitted ends (0.05 pt), plus the two key handles", len(got_ci) == len(exp_ci) + 2 and all(any(abs(a[0] - b[0]) < 0.05 and abs(a[1] - b[1]) < 0.05 and abs(a[2] - b[2]) < 0.05 for a in got_ci) for b in exp_ci), f"{len(got_ci)} vs {len(exp_ci)}")
W0, N0 = L.word_tokens([dict(text=s["text"]) for s in SI]), L.num_tokens([dict(text=s["text"]) for s in SI]); W1, N1 = L.word_tokens([dict(text=s["text"]) for s in SP]), L.num_tokens([dict(text=s["text"]) for s in SP])
wd, nd = L.delta(W0, W1), L.delta(N0, N1); json.dump(dict(words=wd, numbers=nd), open(f"{VER}/text_delta.json", "w"), indent=1, ensure_ascii=False)
exp_n = L.delta(N0, L.num_tokens([dict(text=r["text"]) for r in RECS]))
CK.check("numeric-token delta vs V13 = the expected delta from the drawn records (DECLARED in verify/text_delta.json); word delta = the row labels that enter or leave and the hidden V13 key copies", nd == exp_n, f"lost {sum(nd['lost'].values())} gained {sum(nd['gained'].values())}; words lost {sum(wd['lost'].values())} gained {sum(wd['gained'].values())}")
old_rows = {}
for lab in labs13:
    hs = sorted([s for s in hr13 if abs(s["origin"][1] - lab["origin"][1]) < 15], key=lambda s: s["origin"][1]); qs = sorted([s for s in q13 if abs(s["origin"][1] - lab["origin"][1]) < 15], key=lambda s: s["origin"][1])
    old_rows[lab["text"]] = dict(wake=(hs[0]["text"], qs[0]["text"], "Bold" in qs[0]["font"]), sleep=(hs[1]["text"], qs[1]["text"], "Bold" in qs[1]["font"]))
V.render(OUT, 150, PNG); old150 = f"{WORK}/OLD_V13_150dpi.png"
if not os.path.exists(old150): V.render(BASE, 150, old150)
V.render(BASE, 200, f"{CROPS}/{SHEET}_OLD_V13_200dpi.png"); V.render(OUT, 200, f"{CROPS}/{SHEET}_NEW_200dpi.png")
a_ = Image.open(f"{CROPS}/{SHEET}_OLD_V13_200dpi.png"); b_ = Image.open(f"{CROPS}/{SHEET}_NEW_200dpi.png"); c_ = Image.new("RGB", (a_.width + b_.width + 20, max(a_.height, b_.height)), "white"); c_.paste(a_, (0, 0)); c_.paste(b_, (a_.width + 20, 0)); c_.save(f"{CROPS}/{SHEET}_OLD_vs_NEW_200dpi.png")
CK.check("150 dpi render and 200 dpi OLD vs NEW crop written", os.path.exists(PNG) and os.path.exists(f"{CROPS}/{SHEET}_OLD_vs_NEW_200dpi.png"), "")
V.RULES["tick"] = lambda v: f"{float(v):g}"
n_rows, n_no = V.printed_values_csv(OUT, RECS, f"{VER}/{SHEET}_printed_values.csv", static_from=BASE, group_p_boxes=[])
cnt = collections.Counter(r["match"] for r in csv.DictReader(open(f"{VER}/{SHEET}_printed_values.csv")))
CK.check(f"printed values csv: {n_rows} strings, {cnt.get('yes', 0)} re-derived, {cnt.get('yes_recorded', 0)} recorded, {cnt.get('static', 0)} static, {n_no} NO", n_no == 0, f"{dict(cnt)}")
rows = []
NEWA = {(d_["condition"], d_["period"]): d_ for d_ in DRAWN}
for cond in sorted(set(old_rows) | {r.cond for r in ROWS}):
    for per in ("wake", "sleep"):
        o = old_rows.get(cond, {}).get(per); n = NEWA.get((cond, per)); rs = []
        if o is None: rs.append("row enters")
        elif n is None: rs.append("row leaves")
        else:
            if o[2] != n["sig"]: rs.append("q crosses 0.05")
            om = re.match(r"([\d.]+) \(([\d.]+)-([\d.]+)\)", o[0]); olo, ohi = float(om.group(2)), float(om.group(3))
            if (olo > 1 or ohi < 1) != (n["lo"] > 1 or n["hi"] < 1): rs.append("CI crosses 1")
        rows.append(dict(sheet=SHEET, panel="", item=f"{cond} | {per}-period T90 HR (95% CI), q", before=(f"{o[0]}, q {o[1]}" if o else "absent"), after=(f"{n['hr_ci_text']}, q {n['q_text']}" if n else "absent"), dramatic="yes" if rs else "no", note="; ".join(rs) + ("" if rs else "; v8.2 rebuild of the V13 (v7) sheet")))
rows.append(dict(sheet=SHEET, panel="", item="stratum n (legend)", before="1,981 (v7)", after=f"{N_C:,}", dramatic="no", note="wake_sleep_healthy.json strata.C_organ_metab_free.n"))
rows.append(dict(sheet=SHEET, panel="", item="page height", before=f"{H13:.2f}", after=f"{H:.2f}", dramatic="no", note=f"{N13} -> {N_MAIN + N_CTRL} rows at the V13 pitch"))
rows.append(dict(sheet=SHEET, panel="", item="controls absent from the stratum file", before="", after=", ".join(CTRL_ABSENT), dramatic="no", note="the generator draws the controls the file carries"))
V.write_changes(f"{VER}/CHANGES_{SHEET}.csv", rows)
CK.check(f"CHANGES csv written: {len(rows)} rows, {sum(1 for r in rows if r['dramatic'] == 'yes')} dramatic", True, "")
json.dump(dict(sheet=SHEET, lane="V14_L4_APNEA", round="v8.2", sources={os.path.basename(k): v for k, v in GATES.items()}, n_stratum=N_C, q_threshold=QTHR, n_sig_wake=N_SIG_W, n_sig_sleep=N_SIG_S, n_sig_sleep_up=N_SIG_S_UP,
               rows_main=[r.cond for r in ROWS if not r.neg], rows_controls=[r.cond for r in ROWS if r.neg], controls_absent=CTRL_ABSENT, values=DRAWN, drawn_records=RECS, v13_rows=old_rows,
               geometry=dict(page=[W13, H_NEW], page_v13=[W13, H13], key_nudge_pt=KEY_NUDGE, key_x=KEY_X, pt_plus=L.PT_PLUS, label_x=LAB_X, y_first=Y_FIRST, pitch=PITCH, off_wake=OFF_W, off_sleep=OFF_S, col_right=COL_R, qcol_right=QCOL_R, rule_y=RULE_NEW, rule_y_v13=RULE_Y, axis_x=[AX_X0, AX_X1], x1=X1, b_ln=B_LN, xlim=[XLO, XHI], ticks=TICKS),
               output_sha256=L.sha256(OUT)), open(f"{VER}/{SHEET}_drawn.json", "w"), indent=1, ensure_ascii=False)
ok = CK.write(f"{VER}/checks.txt"); print("RESULT ALL PASS" if ok else "RESULT FAIL")
print(f"rows {N_MAIN} + {N_CTRL} (V13 {N13}); in {sorted({r.cond for r in ROWS} - set(old_rows))}, out {sorted(set(old_rows) - {r.cond for r in ROWS})}; n {N_C:,}; page {H13:.2f} -> {H:.2f}; sig wake {N_SIG_W} sleep {N_SIG_S} (up {N_SIG_S_UP})")
