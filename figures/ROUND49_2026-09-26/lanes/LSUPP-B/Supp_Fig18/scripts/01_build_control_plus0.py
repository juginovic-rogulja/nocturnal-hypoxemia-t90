#!/usr/bin/env python3
"""Supp Fig 18, V14 rebuild at the V13 design (the round-32 LSUPP re-set: 968 pt sheet, the forest axis 5.0 in wide, the block of row
labels, axis and columns centred on the page; the round-33 letters a and b at the block's left edge). Repointed copy of
ROUND32_2026-09-10/LSUPP/Supp_Fig18/scripts/01_build.py (backup 01_build_PRE_V8_1.py), which itself follows the round-30 lane's
ROUND12 build_sheet_de.py code path. Part a (external T90 against total sleep time, SHHS and MrOS) is unchanged: every one of its
strings is asserted equal to the V13 sheet's. Part b: the three primary-cohort rows (heart failure, cardiovascular composite, death;
T90 adjusted for the apnea index and the index adjusted for T90) come from numbers/alenfig3_cross_v1.json ["mutual"] (step 123,
v8.1 sidecar); the external rows come from numbers/external.json (unchanged, no sidecar, external cohorts) printed half up on the
full-precision refit ROUND30_2026-09-08/V7_L19_SHEETFIX2/work/fullprec_external.json (sidecar-checked) and are asserted equal to
the V13 sheet's strings. V14 changes against the LSUPP copy: paths, the v8.1 sidecar gate, the primary strings are NOT asserted
equal to the base (they change by design), the letters are drawn at the block's left edge (the round-33 position, the V13 design),
P values print half up (v14lib.p_text) and the drawn record carries v14lib DRAWN_RECORD fields. Writes work/Supp_Fig18_regen.pdf,
Supp_Fig18.pdf (the regenerated sheet is the delivered sheet, as in round 32) and work/Supp_Fig18_drawn.json."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import json, os, re, sys, hashlib as _hl
from collections import Counter
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerTuple
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
import matplotlib.transforms as mtransforms
import fitz
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LSUPP-B/Supp_Fig18/scripts")   # R49: this lane's copy of lane_common.py
import lane_common as C
VL = C.VL
SHEET = "Supp_Fig18"; HERE = f"{C.LANE}/{SHEET}/work/control_plus0"; WORK = f"{HERE}/work"; VER = f"{HERE}/verify"   # CONTROL: outputs under work/control_plus0
os.makedirs(WORK, exist_ok=True); os.makedirs(VER, exist_ok=True)
BASE = f"{C.BASE}/{SHEET}.pdf"
A3P = f"{C.NUM}/alenfig3_cross_v1.json"; EXTP = f"{C.NUM}/external.json"; COHP = f"{C.NUM}/cohorts.json"
J6D = f"{paths.FIGURE_ROOT}/FINAL_FIGURES_2026-08-14/_workfiles/Figure6D_polish_drawn_values.json"    # part a strings (external, unchanged), layout only
for p in (BASE, A3P, EXTP, COHP, J6D): C.hydrated(p)
SC23 = C.sidecar(A3P); SC01 = C.sidecar(COHP)
for f_ in (C.ARIAL, C.ARIALB): fm.fontManager.addfont(f_)
AX_W = 5.0                                                     # the forest axis width, inches (round 32)
AZ, GREEN_LINE, GREEN_FILL, TERRA = "#0288d1", "#298d32", "#39c445", "#d55e00"
INK = "#1a1d20"
PT_PLUS = 0.0   # CONTROL build (round 49): the same chain at the round-37 sizes from the CURRENT files, to separate the data-file precision effect from the +1 pt effect
FS_LETTER, FS_XLAB, FS_TICK, FS_LAB, FS_LEG, FS_COL = 13.0 + PT_PLUS, 11.0 + PT_PLUS, 10.0 + PT_PLUS, 10.0 + PT_PLUS, 10.0 + PT_PLUS, 9.5 + PT_PLUS
# R49: the forest axis (and every mark) stays where the round-37 build put it: ax_x0 from that lane's drawn record, never typed
R37_DRAWN = f"{paths.FIGURE_ROOT}/ROUND37_2026-09-14/figures/lanes/V14_L3a_DURATION_MAIN/Supp_Fig18/work/Supp_Fig18_drawn.json"
_G37 = json.load(open(R37_DRAWN))["geometry_in"]; assert abs(_G37["a"]["ax_x0"] - _G37["b"]["ax_x0"]) < 1e-9 and _G37["a"]["ax_w"] == 5.0
AX_X0_V26 = float(_G37["a"]["ax_x0"])
A3 = json.load(open(A3P)); E = json.load(open(EXTP)); COH = json.load(open(COHP)); D = json.load(open(J6D))
FPXP = f"{paths.FIGURE_ROOT}/ROUND30_2026-09-08/V7_L19_SHEETFIX2/work/fullprec_external.json"
_sc = json.load(open(C.hydrated(FPXP + ".provenance.json"))); FPX_SHA = _hl.sha256(open(FPXP, "rb").read()).hexdigest()
assert _sc["output_sha256"] == FPX_SHA, "fullprec_external.json does not match its sidecar"
FPX = json.load(open(FPXP))
N_BDSP = COH["bdsp"]["n_analysis"]; assert N_BDSP == 19173 and E["shhs"]["n"] == 5802 and E["mros"]["n"] == 2911
# ---- part a rows (external, unchanged): from the 6D drawn values, ordered as the builder ordered them
rows_d = D["rows"]; assert len(rows_d) == 24
LOOK_D = {(r["cohort"], r["outcome"], r["exposure"]): r for r in rows_d}
def d_outcomes(cohort): return [r["outcome"] for r in rows_d if r["cohort"] == cohort and r["exposure"] == "t90_adj_tst"]
D_BLOCKS = [("shhs", "Sleep Heart\nHealth Study", d_outcomes("shhs")), ("mros", "Osteoporotic\nFractures in Men", d_outcomes("mros"))]
# ---- part b rows: primary from the v8.1 mutual block, external from external.json (printed at full precision)
MUT = A3["mutual"]; PRIMARY = ["Heart failure", "Cardiovascular composite", "Death from any cause"]; LOOK_E = {}
for o in PRIMARY:
    m = MUT[o]; assert not m.get("negative_control", False), o
    for exp, k in (("t90_adj_ahi", "z_t90"), ("ahi_adj_t90", "z_ahi")):
        v = m[k]; LOOK_E[("BDSP", o, exp)] = dict(hr=float(v["hr"]), lo=float(v["lo"]), hi=float(v["hi"]), p=float(v["p"]), ci_excludes_1=not (v["lo"] <= 1.0 <= v["hi"]), source=A3P, key_path=f"mutual/{o}/{k}", n=m.get("n"), events=m.get("events"))
ext_outs = {}
for ck in ("shhs", "mros"):
    outs = list(E[ck]["outcomes"])
    for o in outs:
        for exp in ("t90_adj_ahi", "ahi_adj_t90"):
            v = FPX[ck]["outcomes"][o][exp]; v3 = E[ck]["outcomes"][o][exp]
            assert all(abs(float(v[k]) - float(v3[k])) <= 5e-4 for k in ("hr", "lo", "hi")) and (v["n"], v["events"]) == (v3["n"], v3["events"]), (ck, o, exp)
            LOOK_E[(ck, o, exp)] = dict(hr=float(v3["hr"]), lo=float(v3["lo"]), hi=float(v3["hi"]), p=float(v["p"]), ci_excludes_1=not (v3["lo"] <= 1.0 <= v3["hi"]), print=(float(v["hr"]), float(v["lo"]), float(v["hi"])), source=f"{EXTP} at full precision from {FPXP}", key_path=f"{ck}/outcomes/{o}/{exp}")
    ext_outs[ck] = sorted(outs, key=lambda o: -LOOK_E[(ck, o, "t90_adj_ahi")]["hr"])
prim_sorted = sorted(PRIMARY, key=lambda o: -LOOK_E[("BDSP", o, "t90_adj_ahi")]["hr"])
E_BLOCKS = [("BDSP", "Primary cohort", prim_sorted), ("shhs", "Sleep Heart\nHealth Study", ext_outs["shhs"]), ("mros", "Osteoporotic\nFractures in Men", ext_outs["mros"])]
assert [len(o) for _c, _h, o in E_BLOCKS] == [3, 7, 5] and len(LOOK_E) == 30
for r in LOOK_E.values(): assert r["lo"] < r["hr"] < r["hi"] and r["ci_excludes_1"] == (r["p"] < 0.05), r
def hr_ci_text(hr, lo, hi): return VL.hr_ci(hr, lo, hi)
def p_text(p): return VL.p_text(p)
def ptxt(r): return hr_ci_text(*r["print"]) if "print" in r else hr_ci_text(r["hr"], r["lo"], r["hi"])
# ---- the V13 sheet's text layer: part a and the external rows of part b must be reproduced exactly; the primary rows change by design
bd = fitz.open(BASE); bp = bd[0]; BSP = C.spans_of(bp); bd.close()
RX_HR = re.compile(r"^\d\.\d\d \(\d\.\d\d-\d\.\d\d\)$"); RX_PQ = re.compile(r"^(<0\.001|\d\.\d\d\d)$")
def layer(y0, y1, rx, x0=0): return sorted(s["text"].strip() for s in BSP if y0 <= s["origin"][1] < y1 and s["origin"][0] >= x0 and rx.match(s["text"].strip()))
mine_a_hr = sorted(hr_ci_text(LOOK_D[(ck, o, e)]["hr"], LOOK_D[(ck, o, e)]["lo"], LOOK_D[(ck, o, e)]["hi"]) for ck, _h, outs in D_BLOCKS for o in outs for e in ("t90_adj_tst", "tst_adj_t90"))
mine_a_q = sorted(p_text(LOOK_D[(ck, o, e)]["q"]) for ck, _h, outs in D_BLOCKS for o in outs for e in ("t90_adj_tst", "tst_adj_t90"))
assert mine_a_hr == layer(0, 570, RX_HR), ("part a HR strings differ from the V13 sheet", Counter(mine_a_hr) - Counter(layer(0, 570, RX_HR)), Counter(layer(0, 570, RX_HR)) - Counter(mine_a_hr))
assert mine_a_q == layer(0, 570, RX_PQ), ("part a q strings differ from the V13 sheet", Counter(mine_a_q) - Counter(layer(0, 570, RX_PQ)))
# part b on the V13 sheet: the primary block sits under the head "Primary cohort" (three rows), the external blocks below it
prim_head = [s for s in BSP if s["text"].strip() == "Primary cohort"]; assert len(prim_head) == 1, prim_head
shhs_head = [s for s in BSP if s["text"].strip() == "Sleep Heart" and s["origin"][1] > 600]; assert len(shhs_head) == 1, shhs_head
Y_PRIM0, Y_PRIM1 = prim_head[0]["origin"][1], shhs_head[0]["origin"][1]
base_prim_hr = layer(Y_PRIM0, Y_PRIM1, RX_HR); base_prim_p = layer(Y_PRIM0, Y_PRIM1, RX_PQ); assert len(base_prim_hr) == 6 and len(base_prim_p) == 6, (base_prim_hr, base_prim_p)
mine_ext_hr = sorted(ptxt(r) for (ck, o, e), r in LOOK_E.items() if ck != "BDSP"); mine_ext_p = sorted(p_text(r["p"]) for (ck, o, e), r in LOOK_E.items() if ck != "BDSP")
assert mine_ext_hr == layer(Y_PRIM1, 1230, RX_HR), ("external HR strings of part b differ from the V13 sheet", Counter(mine_ext_hr) - Counter(layer(Y_PRIM1, 1230, RX_HR)), Counter(layer(Y_PRIM1, 1230, RX_HR)) - Counter(mine_ext_hr))
assert mine_ext_p == layer(Y_PRIM1, 1230, RX_PQ), ("external P strings of part b differ from the V13 sheet", Counter(mine_ext_p) - Counter(layer(Y_PRIM1, 1230, RX_PQ)))
mine_prim_hr = sorted(ptxt(r) for (ck, o, e), r in LOOK_E.items() if ck == "BDSP"); mine_prim_p = sorted(p_text(r["p"]) for (ck, o, e), r in LOOK_E.items() if ck == "BDSP")
print("part a (24 + 24 strings) and the 24 external rows of part b reproduce the V13 sheet; primary rows V13", base_prim_hr, base_prim_p, "-> v8.1", mine_prim_hr, mine_prim_p)
# ---- geometry, the ROUND12 builder's, with the block centred and the axis AX_W wide (round 32), letters at the block edge (round 33)
W_PT = 968.0; W = W_PT / 72.0; MARG = 14.0 / 72.0; ROW_IN = 0.43; OFF = 0.21
HEAD_BAND = 0.70; BOT_BAND = 0.58; PART_GAP = 0.30; LEG_ROW = 0.21
XLIM = (0.635, 1.86); XTICKS = [0.7, 0.8, 1, 1.2, 1.5, 1.8]
_all_x = [v for r in list(LOOK_D.values()) + list(LOOK_E.values()) for v in (r["hr"], r["lo"], r["hi"])]
assert XLIM[0] < min(_all_x) and max(_all_x) < XLIM[1], (min(_all_x), max(_all_x))
def layout(blocks):
    items, y = [], 0.0
    for gi, (ck, head, outs) in enumerate(blocks):
        if gi > 0: y -= 0.35
        items.append(("head", (ck, head), y)); y -= 0.95
        for o in outs: items.append(("row", (ck, o), y)); y -= 1.0
    return items
KEY_A = ("T90, adjusted for total sleep time", "Total sleep time, adjusted for T90", "Not significant")
KEY_B = ("Oxygen, adjusted for the apnea index", "Apnea index, adjusted for oxygen", "95% CI includes 1")
PARTS = [dict(letter="a", blocks=D_BLOCKS, look=LOOK_D, exp=("t90_adj_tst", "tst_adj_t90"), line_cols=(AZ, GREEN_LINE), fill_cols=(AZ, GREEN_FILL), marks=("o", "s"), sig_of=lambda r: r["significant_fdr"], pcol="q", pkey="q", keys=KEY_A, n_rows=24),
         dict(letter="b", blocks=E_BLOCKS, look=LOOK_E, exp=("t90_adj_ahi", "ahi_adj_t90"), line_cols=(AZ, TERRA), fill_cols=(AZ, TERRA), marks=("o", "s"), sig_of=lambda r: r["ci_excludes_1"], pcol="P", pkey="p", keys=KEY_B, n_rows=30)]
for P in PARTS:
    P["items"] = layout(P["blocks"]); P["y_top"] = 0.72; P["y_bot"] = min(y for _k, _p, y in P["items"]) - 0.72
    P["ax_h"] = ROW_IN * (P["y_top"] - P["y_bot"]); P["h"] = HEAD_BAND + P["ax_h"] + BOT_BAND
H = round(MARG + PARTS[0]["h"] + PART_GAP + PARTS[1]["h"] + MARG, 3)
RC = {"font.family": "Arial", "font.size": FS_TICK, "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.bbox": "standard", "figure.dpi": 100,
      "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8, "xtick.major.size": 3.0, "ytick.major.size": 0.0, "xtick.direction": "out",
      "axes.edgecolor": INK, "text.color": INK, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK, "legend.frameon": False}
ALL_ROW_KEYS = [pl for P in PARTS for k, pl, _y in P["items"] if k == "row"]
ALL_HEADS = [pl[1] for P in PARTS for k, pl, _y in P["items"] if k == "head"]; HEADLINES = {ln for h in ALL_HEADS for ln in h.split("\n")}
ALL_COL_TXT, ALL_PQ_TXT = [], []
for P in PARTS:
    for ck, o in [pl for k, pl, _y in P["items"] if k == "row"]:
        for exp in P["exp"]:
            r = P["look"][(ck, o, exp)]; ALL_COL_TXT.append(ptxt(r)); ALL_PQ_TXT.append(p_text(r[P["pkey"]]))
def draw_part(fig, P, part_top):
    FW, FH = fig.get_size_inches(); rend = fig.canvas.get_renderer()
    def w_in(s, fs, bold=False):
        t = fig.text(0.5, 0.5, s, fontsize=fs, fontweight="bold" if bold else "normal"); w = t.get_window_extent(renderer=rend).width / fig.dpi; t.remove(); return w
    look, (ex1, ex2), (l1, l2), (f1, f2), (m1, m2) = P["look"], P["exp"], P["line_cols"], P["fill_cols"], P["marks"]
    sig_of, pkey = P["sig_of"], P["pkey"]; items = P["items"]
    GUT_W = max([w_in(o, FS_LAB) for _c, o in ALL_ROW_KEYS] + [w_in(ln, FS_LAB, bold=True) for ln in HEADLINES]) + 0.03
    GUT_PAD, COL_GAP, PQ_GAP = 0.05, 0.14, 0.10
    COL_W = max(w_in(s, FS_COL) for s in ALL_COL_TXT + ["HR (95% CI)"]); PQ_W = max(w_in(s, FS_COL, bold=True) for s in ALL_PQ_TXT + ["q", "P"])
    block_w = GUT_W + GUT_PAD + AX_W + COL_GAP + COL_W + PQ_GAP + PQ_W
    x_block_centred = (W - block_w) / 2.0                           # the round-37 rule (the block centred on the page)
    ax_x0 = AX_X0_V26                                               # R49: the axis pinned at the round-37 position (no mark moves)
    x_block = ax_x0 - GUT_W - GUT_PAD                               # the label block grows to the left, the columns to the right
    col_x = ax_x0 + AX_W + COL_GAP + COL_W; pq_x = col_x + PQ_GAP + PQ_W
    ax_top = part_top - HEAD_BAND; ax_y0 = ax_top - P["ax_h"]
    ax = fig.add_axes([ax_x0 / FW, ax_y0 / FH, AX_W / FW, P["ax_h"] / FH])
    ax.axvline(1.0, color=INK, lw=0.9, ls=(0, (4, 3)), zorder=1)   # no pale vertical grid lines (round 31)
    trH = mtransforms.blended_transform_factory(fig.dpi_scale_trans, ax.transData)
    DRAWN = []; yticks, ylabels = [], []
    for kind, payload, yy in items:
        if kind == "head":
            fig.text(x_block, yy, payload[1], transform=trH, fontsize=FS_LAB, fontweight="bold", color=INK, ha="left", va="center"); continue
        ck, outc = payload; yticks.append(yy); ylabels.append(outc)
        for exp, lc, fc, mk, offs, s_pt2 in ((ex1, l1, f1, m1, OFF, 36), (ex2, l2, f2, m2, -OFF, 31)):
            r = look[(ck, outc, exp)]; hr, lo, hi, pv = r["hr"], r["lo"], r["hi"], r[pkey]; sig = bool(sig_of(r))
            ax.plot([lo, hi], [yy + offs] * 2, color=lc, lw=1.7, solid_capstyle="round", zorder=2)
            ax.scatter([hr], [yy + offs], s=s_pt2, marker=mk, zorder=3, facecolors=fc if sig else "white", edgecolors="white" if sig else lc, linewidths=0.8)
            fig.text(col_x, yy + offs, ptxt(r), transform=trH, fontsize=FS_COL, color=INK, ha="right", va="center")
            fig.text(pq_x, yy + offs, p_text(pv), transform=trH, fontsize=FS_COL, color=INK, ha="right", va="center", fontweight="bold" if sig else "normal")
            DRAWN.append({"cohort": ck, "outcome": outc, "exposure": exp, "hr": hr, "lo": lo, "hi": hi, pkey: pv, "significant": sig, "printed_hr_ci": ptxt(r), f"printed_{pkey}": p_text(pv), "fill": "filled" if sig else "open", "row_unit_y": yy + offs, "source": r.get("source"), "key_path": r.get("key_path")})
    ax.set_yticks(yticks); ax.set_yticklabels(ylabels, fontsize=FS_LAB, ha="left"); ax.tick_params(axis="y", pad=(GUT_W + GUT_PAD) * 72.0, length=0)
    for spn in ("top", "right", "left"): ax.spines[spn].set_visible(False)
    ax.set_xscale("log"); ax.xaxis.set_major_locator(FixedLocator(XTICKS)); ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.xaxis.set_minor_locator(NullLocator()); ax.tick_params(axis="x", labelsize=FS_TICK); ax.set_xlim(*XLIM); ax.set_ylim(P["y_bot"], P["y_top"])
    head_y = ax_top + 0.07
    fig.text(col_x / FW, head_y / FH, "HR (95% CI)", fontsize=FS_COL, color=INK, ha="right", va="bottom")
    fig.text(pq_x / FW, head_y / FH, P["pcol"], fontsize=FS_COL, color=INK, ha="right", va="bottom")
    fig.text((ax_x0 + AX_W / 2) / FW, (ax_y0 - BOT_BAND + 0.08) / FH, "Hazard ratio per 1 SD (95% CI)", fontsize=FS_XLAB, color=INK, ha="center", va="bottom")
    fig.text(x_block / FW, part_top / FH, P["letter"], fontsize=FS_LETTER, fontweight="bold", color=INK, ha="left", va="top")   # round 33: the letter at the block's left edge
    k1, k2, k_open = P["keys"]; ms_o, ms_s = 6.0, 5.6
    fig.legend(handles=[Line2D([], [], color=l1, marker=m1, ms=ms_o, lw=1.7, markerfacecolor=f1, markeredgecolor="white", markeredgewidth=0.8, label=k1),
                        Line2D([], [], color=l2, marker=m2, ms=ms_s, lw=1.7, markerfacecolor=f2, markeredgecolor="white", markeredgewidth=0.8, label=k2)],
               loc="upper left", bbox_to_anchor=(ax_x0 / FW, part_top / FH), ncol=2, fontsize=FS_LEG, frameon=False, handletextpad=0.5, columnspacing=1.3, handlelength=1.6, borderpad=0.0, borderaxespad=0.0)
    fig.legend(handles=[(Line2D([], [], marker=m1, ls="none", ms=ms_o, markerfacecolor="white", markeredgecolor=l1, markeredgewidth=1.0),
                         Line2D([], [], marker=m2, ls="none", ms=ms_s, markerfacecolor="white", markeredgecolor=l2, markeredgewidth=1.0))],
               labels=[k_open], handler_map={tuple: HandlerTuple(ndivide=None, pad=0.35)}, loc="upper left", bbox_to_anchor=(ax_x0 / FW, (part_top - LEG_ROW) / FH),
               ncol=1, fontsize=FS_LEG, frameon=False, handletextpad=0.6, handlelength=2.2, borderpad=0.0, borderaxespad=0.0)
    return dict(drawn=DRAWN, x_block=x_block, gut_w=GUT_W, ax_x0=ax_x0, ax_w=AX_W, ax_y0=ax_y0, ax_top=ax_top, col_x=col_x, pq_x=pq_x, part_top=part_top, block_w=block_w,
                x_block_centred_rule=x_block_centred, ax_x0_centred_rule=x_block_centred + GUT_W + GUT_PAD, col_w=COL_W, pq_w=PQ_W)
with plt.rc_context(RC):
    fig = plt.figure(figsize=(W, H))
    top_a = H - MARG; res_a = draw_part(fig, PARTS[0], top_a)
    top_b = top_a - PARTS[0]["h"] - PART_GAP; res_b = draw_part(fig, PARTS[1], top_b)
    raw = f"{WORK}/{SHEET}_regen_raw.pdf"; fig.savefig(raw); plt.close(fig)
regen = f"{WORK}/{SHEET}_regen.pdf"
aligned = C.align_font_metrics(raw, regen, BASE)
d = fitz.open(regen); pg = d[0]; assert abs(pg.rect.width - 968.0) < 0.6 and abs(pg.rect.height - H * 72) < 0.6, pg.rect
n_t90_excl = sum(1 for (ck, o, exp), r in LOOK_E.items() if exp == "t90_adj_ahi" and r["ci_excludes_1"]); n_ahi_incl = sum(1 for (ck, o, exp), r in LOOK_E.items() if exp == "ahi_adj_t90" and not r["ci_excludes_1"])
ahi_excl = [(ck, o) for (ck, o, exp), r in LOOK_E.items() if exp == "ahi_adj_t90" and r["ci_excludes_1"]]
# drawn records (v14lib) for the six primary HR strings, six P strings and three primary row labels, from the regenerated sheet's text layer
sp = C.spans_of(pg); d.close()
RECORDS = []
def find_span(text, y_pt, tol=1.0):
    c = [s for s in sp if s["text"].strip() == text and abs(s["origin"][1] - y_pt) < tol]; assert len(c) == 1, (text, y_pt, c); return c[0]
def y_pt_of(row_unit_y): return (H - (res_b["ax_y0"] + (row_unit_y - PARTS[1]["y_bot"]) / (PARTS[1]["y_top"] - PARTS[1]["y_bot"]) * PARTS[1]["ax_h"])) * 72.0
for r in res_b["drawn"]:
    if r["cohort"] != "BDSP": continue
    yc = y_pt_of(r["row_unit_y"])
    s_hr = min((s for s in sp if s["text"].strip() == r["printed_hr_ci"]), key=lambda s: abs((s["bbox"][1] + s["bbox"][3]) / 2 - yc)); assert abs((s_hr["bbox"][1] + s_hr["bbox"][3]) / 2 - yc) < 3, (r, s_hr)
    s_p = min((s for s in sp if s["text"].strip() == r["printed_p"]), key=lambda s: abs((s["bbox"][1] + s["bbox"][3]) / 2 - yc)); assert abs((s_p["bbox"][1] + s_p["bbox"][3]) / 2 - yc) < 3, (r, s_p)
    RECORDS.append(dict(text=r["printed_hr_ci"], x=s_hr["bbox"][2], baseline=s_hr["origin"][1], ha="right", size=s_hr["size"], panel="b", source_file=A3P, source_key=r["key_path"], source_value={"hr": r["hr"], "lo": r["lo"], "hi": r["hi"]}, rule="hr_ci_dict", note=f"{r['outcome']} {r['exposure']}"))
    RECORDS.append(dict(text=r["printed_p"], x=s_p["bbox"][2], baseline=s_p["origin"][1], ha="right", size=s_p["size"], panel="b", source_file=A3P, source_key=r["key_path"] + "/p", source_value=r["p"], rule="p3", note=f"{r['outcome']} {r['exposure']}"))
prim_y = {r["outcome"]: y_pt_of(r["row_unit_y"] - OFF) for r in res_b["drawn"] if r["cohort"] == "BDSP" and r["exposure"] == "t90_adj_ahi"}   # the row centre of each primary row
for o in prim_sorted:
    s_l = min((s for s in sp if s["text"].strip() == o), key=lambda s: abs((s["bbox"][1] + s["bbox"][3]) / 2 - prim_y[o])); assert abs((s_l["bbox"][1] + s_l["bbox"][3]) / 2 - prim_y[o]) < 4, (o, s_l, prim_y[o])
    RECORDS.append(dict(text=o, x=s_l["origin"][0], baseline=s_l["origin"][1], ha="left", size=s_l["size"], panel="b", source_file=A3P, source_key=f"mutual/{o}", source_value=o, rule="keyname", note="primary-cohort row label"))
out = {"figure": SHEET, "lane": "V14_L3a_DURATION_MAIN", "regenerated_sheet": regen, "base_sheet": BASE, "base_sha256": C.sha256(BASE),
       "sources": {"alenfig3_cross_v1.json": SC23, "external.json": {"path": EXTP, "sha256": C.sha256(EXTP), "sidecar": "none (external cohorts, unchanged)", "printed_from": {"path": FPXP, "sha256": FPX_SHA, "rule": "half up on the full-precision refit"}},
                   "cohorts.json": SC01, "Figure6D_polish_drawn_values.json": {"path": J6D, "sha256": C.sha256(J6D), "role": "part a strings (external, unchanged), asserted equal to V13"}},
       "colours": {"oxygen": AZ, "sleep_line": GREEN_LINE, "sleep_fill": GREEN_FILL, "apnea_index": TERRA, "ink": INK}, "font_metrics_aligned": aligned, "ax_w_in": AX_W,
       "pattern": {"t90_ci_excludes_1": n_t90_excl, "of": 15, "ahi_ci_includes_1": n_ahi_incl, "ahi_excludes_1": ahi_excl},
       "primary_v13_strings": {"hr": base_prim_hr, "p": base_prim_p}, "primary_v14_strings": {"hr": mine_prim_hr, "p": mine_prim_p},
       "rows_a": res_a["drawn"], "rows_b": res_b["drawn"], "xlim": XLIM, "xticks": XTICKS, "page_pt": [W_PT, round(H * 72, 3)],
       "geometry_in": {k: {kk: res[kk] for kk in ("x_block", "gut_w", "ax_x0", "ax_w", "ax_y0", "ax_top", "col_x", "pq_x", "part_top", "block_w", "x_block_centred_rule", "ax_x0_centred_rule", "col_w", "pq_w")} for k, res in (("a", res_a), ("b", res_b))},
       "r49": {"pt_plus": PT_PLUS, "sizes": dict(letter=FS_LETTER, xlab=FS_XLAB, tick=FS_TICK, lab=FS_LAB, leg=FS_LEG, col=FS_COL), "ax_x0_pinned_from": R37_DRAWN, "geometry_v26_in": _G37},
       "records": RECORDS}
json.dump(out, open(f"{WORK}/{SHEET}_drawn.json", "w"), indent=1, default=float)
import shutil; shutil.copyfile(regen, f"{HERE}/{SHEET}.pdf")
for r in res_b["drawn"]:
    if r["cohort"] == "BDSP": print(f"  {r['outcome']:26s} {r['exposure']:12s} {r['printed_hr_ci']}  P {r['printed_p']}  {r['fill']}  (n {LOOK_E[('BDSP', r['outcome'], r['exposure'])]['n']}, events {LOOK_E[('BDSP', r['outcome'], r['exposure'])]['events']})")
print(f"oxygen adjusted for the index excludes 1 in {n_t90_excl} of 15; the index adjusted for oxygen includes 1 in {n_ahi_incl} of 15 (excludes 1: {ahi_excl})")
print(f"wrote {HERE}/{SHEET}.pdf ({W * 72:.1f} x {H * 72:.2f} pt), letters at x {res_a['x_block'] * 72:.3f} pt (block edge)")
