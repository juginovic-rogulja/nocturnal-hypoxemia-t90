#!/usr/bin/env python3
"""ROUND 54 (2026-09-29, lane L4): Main Fig 4 panel c at the round-54 sizes, re-laid on work/layout_r54.json (the panel is a full-width
row under a | b, c beside d does not fit the sheet width at 14 pt). The round-49 data wiring is kept (numbers/crosstab_v2.json rows3,
v8.1 sidecar gate, the ladder rule, every value asserted against the file). Sizes: the title 15 pt, band labels, row labels, cells and
the key 14 pt. Cells widened to 96 pt (the widest cell string is 86.9 pt at 14 pt), rows 44 pt tall (one line of 14 pt), the key
row under the grid (title on two lines, four rungs with their labels). Colours as V13: the four blue header rectangles, the orange
ladder by percent of the group, white text on the darkest rung. Outputs: work/panel_c_raw.pdf, work/panel_c_geometry.json,
verify/Main_Fig4_c_drawn.json."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import l4lib as L, layout_r54 as LY
matplotlib, plt = L.mpl_setup()
from matplotlib.patches import Rectangle

SHEET = "Main_Fig4"; SD = f"{L.LANE}/{SHEET}"; WORK = f"{SD}/work"; VER = f"{SD}/verify"
S = LY.load(); C_ = S["c"]; FS = S["sizes"]
SRC = f"{L.NUM}/crosstab_v2.json"; SC = L.sidecar(SRC); assert SC["step"]["rc"] == 0, SC["step"]
COH = json.load(open(L.hydrated(f"{L.NUM}/cohorts.json"))); L.sidecar(f"{L.NUM}/cohorts.json"); N_COHORT = int(COH["bdsp"]["n_analysis"])
XT = json.load(open(L.hydrated(SRC))); R3 = XT["rows3"]; assert len(R3) == 3 and XT["n_cohort"] == N_COHORT, (XT["n_cohort"], N_COHORT)
assert [r["category"] for r in R3] == ["No or mild (AHI <15)", "Moderate (AHI 15 to <30)", "Severe (AHI >=30)"]
assert sum(r["n"] for r in R3) == XT["n"] == N_COHORT - XT["n_ahi_undefined"], (sum(r["n"] for r in R3), XT["n"], XT["n_ahi_undefined"])
R4 = {r["category"]: r for r in XT["rows"]}; BN = ["0-1%", "1-5%", "5-10%", ">10%"]
for bn in BN:
    assert R3[0][f"n_{bn}"] == R4["None (AHI <5)"][f"n_{bn}"] + R4["Mild (AHI 5 to <15)"][f"n_{bn}"], bn
    assert R3[1][f"n_{bn}"] == R4["Moderate (AHI 15 to <30)"][f"n_{bn}"] and R3[2][f"n_{bn}"] == R4["Severe (AHI >=30)"][f"n_{bn}"], bn
for r in R3: assert sum(r[f"n_{bn}"] for bn in BN) == r["n"], r
def rung(p): return 0 if p < 5 else 1 if p < 15 else 2 if p < 35 else 3
ROWS = [("No or mild", "(AHI <15)"), ("Moderate", "(AHI 15 to <30)"), ("Severe", "(AHI 30 or more)")]
COLS = ["≤1%", ">1–5%", ">5–10%", ">10%"]; RUNG_LAB = ["<5", "5–15", "15–35", "≥35"]
AZURE = L.BLUE_LADDER; LADDER = L.ORANGE_LADDER; INK = L.INK
for t in ("T90, % of the recording", *COLS, *RUNG_LAB, "Percent of the", "apnea group", *[l for r in ROWS for l in r]): assert L.house_ok(t), t

# ---------------------------------------------------------------- geometry from the layout spec (page points, y down)
BX0, BY0, BX1, BY1 = C_["box"]; COL_X = C_["col_x"]; CW = C_["cell_w"]; HY0, HY1 = C_["header_y"]; ROW_TOP0 = C_["row_top0"]; ROW_H = C_["row_h"]; PITCH = ROW_H + C_["row_gap"]
GRID_BOTTOM = C_["grid_bottom"]; LAB_X1 = C_["lab_x1"]; GRID_X0, GRID_X1 = COL_X[0], C_["grid_x1"]
FS_T, FS_C = FS["head"], FS["clab"]; XH = 0.36 * FS_C; L1, L2 = -3.5, 11.5      # single line: baseline = centre + x-height/2; two lines at 15 pt pitch around the centre
assert abs(ROW_TOP0 + 2 * PITCH + ROW_H - GRID_BOTTOM) < 1e-6
W, H = BX1 - BX0, BY1 - BY0
fig = plt.figure(figsize=(W / 72, H / 72)); fig.patch.set_alpha(0.0); ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(BX0, BX1); ax.set_ylim(BY1, BY0); ax.axis("off"); ax.patch.set_visible(False)
DR = []; STR = []; RECTS = []
def text(x, y, s, size, ha, color=INK, **src):
    ax.text(x, y, s, fontsize=size, ha=ha, va="baseline", color=color); STR.append(dict(text=s, x=round(x, 4), baseline=round(y, 4), size=size, ha=ha, color=color, panel="c", **src)); return ax.texts[-1]
def rect(x0, y0, w, h, fill, edge="none", lw=0, kind=""):
    ax.add_patch(Rectangle((x0, y0), w, h, facecolor=fill, edgecolor=edge, lw=lw)); RECTS.append(dict(kind=kind, rect=[round(x0, 3), round(y0, 3), round(x0 + w, 3), round(y0 + h, 3)], fill=fill, edge=edge, lw=lw))
RC = dict(L.RC); RC.update({"figure.facecolor": "none", "savefig.facecolor": "none"})
with plt.rc_context(RC):
    rend = fig.canvas.get_renderer()
    text((GRID_X0 + GRID_X1) / 2, C_["title_base"], "T90, % of the recording", FS_T, "center", source_file="static:V13", source_key="V13 title", source_value="T90, % of the recording", rule="text")
    hc = (HY0 + HY1) / 2
    for i, x in enumerate(COL_X):
        rect(x, HY0, CW, HY1 - HY0, AZURE[i], kind=f"T90 header {COLS[i]}"); text(x + CW / 2, hc + XH, COLS[i], FS_C, "center", source_file="static:V13", source_key="V13 band label", source_value=COLS[i], rule="text")
    for ri, r in enumerate(R3):
        y0 = ROW_TOP0 + ri * PITCH; y1 = y0 + ROW_H; cy = (y0 + y1) / 2
        text(LAB_X1, cy + L1, ROWS[ri][0], FS_C, "right", source_file="static:V13", source_key="V13 row label", source_value=ROWS[ri][0], rule="text")
        text(LAB_X1, cy + L2, ROWS[ri][1], FS_C, "right", source_file="static:V13", source_key="V13 row label", source_value=ROWS[ri][1], rule="text")
        for ci, bn in enumerate(BN):
            n = int(r[f"n_{bn}"]); p = L.pct1(n, r["n"]); assert abs(float(p) - float(r[f"pct_{bn}"])) < 0.051, (r["category"], bn, p, r[f"pct_{bn}"])
            k = rung(float(p)); fill = LADDER[k]; col = "#ffffff" if k == 3 else INK; s = f"{n:,} ({p})"
            rect(COL_X[ci], y0, CW, ROW_H, fill, kind=f"cell {ROWS[ri][0]}|{COLS[ci]}")
            text(COL_X[ci] + CW / 2, cy + XH, s, FS_C, "center", color=col, source_file=SRC, source_key=f"rows3/{ri}", source_value=[n, r["n"]], rule=f"count_pct_{bn}")
            DR.append(dict(row=r["category"], row_label=ROWS[ri][0], col=COLS[ci], n=n, pct=p, pct_file=r[f"pct_{bn}"], group_n=r["n"], text=s, fill=fill, rung=RUNG_LAB[k], text_color=col,
                           rect=[COL_X[ci], round(y0, 3), COL_X[ci] + CW, round(y1, 3)], key=f"rows3[{ri}].n_{bn} / n", source=SRC))
    # the key row: 'Percent of the' / 'apnea group' on two lines, then the four rungs with their labels centred on the two-line block
    kb1, kb2 = C_["key_base"]; kx = GRID_X0 + 20.0
    t1 = text(kx, kb1, "Percent of the", FS_C, "left", source_file="static:V13", source_key="V13 key title", source_value="Percent of the", rule="text")
    t2 = text(kx, kb2, "apnea group", FS_C, "left", source_file="static:V13", source_key="V13 key title", source_value="apnea group", rule="text")
    tw = max(t.get_window_extent(renderer=rend).width / fig.dpi * 72 for t in (t1, t2))
    kcy = (kb1 - 0.716 * FS_C + kb2) / 2; SW, SH = 16.0, 12.0; x = kx + tw + 14.0; KEY_RECTS = []
    for k in range(4):
        rect(x, kcy - SH / 2, SW, SH, LADDER[k], edge="white", lw=0.6, kind=f"key rung {RUNG_LAB[k]}"); KEY_RECTS.append(RECTS[-1])
        t = text(x + SW + 4.0, kcy + XH, RUNG_LAB[k], FS_C, "left", source_file="static:V13", source_key="V13 key rung label", source_value=RUNG_LAB[k], rule="text")
        x = x + SW + 4.0 + t.get_window_extent(renderer=rend).width / fig.dpi * 72 + 14.0
    BOTTOM_INK = max(s["baseline"] for s in STR)
    for t in ax.texts:
        bb = t.get_window_extent(renderer=rend); x0p, x1p = BX0 + bb.x0 / fig.dpi * 72, BX0 + bb.x1 / fig.dpi * 72; y0p, y1p = BY1 - bb.y1 / fig.dpi * 72, BY1 - bb.y0 / fig.dpi * 72
        assert x0p >= BX0 + 1.0 and x1p <= BX1 - 1.0 and y0p >= BY0 + 1.0 and y1p <= BY1 - 1.0, ("text leaves panel c's page", t.get_text(), x0p, x1p, y0p, y1p)
        assert t.get_fontsize() >= 13.0 and L.house_ok(t.get_text())
    fig.savefig(f"{WORK}/panel_c_raw.pdf", transparent=True); plt.close(fig)
# ---------------------------------------------------------------- proofs: values from the file, the V13 string set
G0 = json.load(open(f"{WORK}/v13_geometry.json"))["panel_c"]
R49 = json.load(open(f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/L4/Main_Fig4/verify/Main_Fig4_c_drawn.json"))["drawn_records"]   # V30's panel c strings (the round-49 build is the source of V30's raster)
assert sorted(s["text"] for s in STR) == sorted(s["text"] for s in R49), ("the string set differs from V30's", sorted(set(s["text"] for s in STR) ^ set(s["text"] for s in R49)))
old_by = {(o["row"], o["col"]): o for o in G0["v13_values"]}; assert len(old_by) == len(DR) == 12
CHANGES = [dict(row=d["row_label"], col=d["col"], before=old_by[(d["row"], d["col"])]["text"], after=d["text"]) for d in DR if d["text"] != old_by[(d["row"], d["col"])]["text"] or d["fill"] != old_by[(d["row"], d["col"])]["fill"]]
geom = dict(panel_box=[BX0, BY0, BX1, BY1], col_x=COL_X, cell_w=CW, header_y=[HY0, HY1], row_top0=ROW_TOP0, row_h=ROW_H, pitch=PITCH, grid_bottom=GRID_BOTTOM, grid_x=[GRID_X0, GRID_X1], label_x1=LAB_X1,
            key_base=C_["key_base"], rects=RECTS, key_rects=KEY_RECTS, bottom_ink_baseline=BOTTOM_INK, design="round-54 re-lay on layout_r54.json, a standalone page (no letter, no title strip)")
json.dump(geom, open(f"{WORK}/panel_c_geometry.json", "w"), indent=1)
json.dump(dict(sheet=SHEET, panel="c", lane="L4 round 54", source=SRC, source_sha256=L.sha256(SRC), sidecar_step=SC["step"], v14_gate=SC["_v14_gate"],
               n=XT["n"], n_cohort=XT["n_cohort"], n_ahi_undefined=XT["n_ahi_undefined"], spearman=XT["spearman"], rows=[dict(category=r["category"], n=r["n"], t90_median=r["t90_median"]) for r in R3],
               values=DR, strings=STR, drawn_records=STR, changes_vs_v13=CHANGES, geometry=geom), open(f"{VER}/{SHEET}_c_drawn.json", "w"), indent=1, ensure_ascii=False)
print(f"panel c: 3 rows x 4 bands (cells {CW:.1f} x {ROW_H:.1f} pt, grid x {GRID_X0:.1f}..{GRID_X1:.1f}, y {HY0:.1f}..{GRID_BOTTOM:.1f}, lowest baseline {BOTTOM_INK:.1f}); n {XT['n']:,} of {XT['n_cohort']:,}; {len(CHANGES)} of 12 cells differ from V13 ({len(STR)} strings); wrote {WORK}/panel_c_raw.pdf")
