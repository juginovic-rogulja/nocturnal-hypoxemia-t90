#!/usr/bin/env python3
"""ROUND 49 (2026-09-26, lane L4): copy of the round-38 builder with every text one point larger (DPT = 1.0 inside text(): title 11 -> 12,
all other strings 10 -> 11); positions, cells, colours and strings unchanged. Round-38 docstring follows.
V14 lane L4 (round 37), Main Fig 4 panel c: the apnea-by-T90 count grid (numbers/crosstab_v2.json key rows3, three apnea groups x
four T90 bands, orange ladder by percent of the group, the four blue header rectangles, the key "Percent of the apnea group") from
the v8.1 numbers, at EXACTLY the V13 geometry (the round-31 LF4 panel: rows 67.004 pt tall, grid bottom level with panel d's, the
key's lowest baseline level with panel d's ramp ticks; LF4/Main_Fig4/work/panel_c_geometry.json, proved to be the V13 panel by
00_probe_v13.py). Copy of ROUND31_2026-09-10/LF4_FIG4_ED7/Main_Fig4/scripts/01_build_c.py (beside as 01_build_c_PRE_V8_1.py),
repointed: v8.1 sidecar gate, no round-30 drawn-record identity assertions (the values change), every string position pinned to the
V13 drawn record, every value asserted against the file. Outputs: work/panel_c_raw.pdf, work/panel_c_geometry.json, verify/Main_Fig4_c_drawn.json."""
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import l4lib as L
matplotlib, plt = L.mpl_setup()
from matplotlib.patches import Rectangle

SHEET = "Main_Fig4"; SD = f"{L.LANE}/{SHEET}"; WORK = f"{SD}/work"; VER = f"{SD}/verify"
DPT = 1.0   # round 49: every text element one point larger, applied inside text()
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

# ---------------------------------------------------------------- geometry: the V13 panel's (LF4 geometry, every string at its V13 origin)
G13 = json.load(open(f"{WORK}/v13_geometry.json")); G0 = G13["panel_c"]; OS = {}
for s in G0["v13_strings"]: OS.setdefault(s["text"], s)
BX0, BY0, BX1, BY1 = G0["panel_box"]; COL_X = G0["col_x"]; CW = G0["cell_w"]; HY0, HY1 = G0["header_y"]; ROW_TOP0 = G0["row_top0"]; ROW_H = G0["row_h"]; PITCH = G0["pitch"]
GRID_BOTTOM = G0["grid_bottom"]; LAB_X1 = G0["label_x1"]; LAB_DY = G0["label_dy"]; CELL_DY = G0["cell_text_dy"]; KEY_RECTS = [r for r in G0["rects"] if r["kind"].startswith("key rung")]
assert len(KEY_RECTS) == 4 and abs(ROW_TOP0 + 2 * PITCH + ROW_H - GRID_BOTTOM) < 1e-6
for t in ("T90, % of the recording", *COLS, *RUNG_LAB, "Percent of the", "apnea group", *[l for r in ROWS for l in r]): assert t in OS, ("V13 string missing", t)
W, H = BX1 - BX0, BY1 - BY0
fig = plt.figure(figsize=(W / 72, H / 72)); fig.patch.set_alpha(0.0); ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(BX0, BX1); ax.set_ylim(BY1, BY0); ax.axis("off"); ax.patch.set_visible(False)
DR = []; STR = []; RECTS = []
def text(x, y, s, size, ha, color=INK, **src):
    ax.text(x, y, s, fontsize=size + DPT, ha=ha, va="baseline", color=color); STR.append(dict(text=s, x=round(x, 4), baseline=round(y, 4), size=size + DPT, ha=ha, color=color, panel="c", **src))   # round 49: +DPT
def rect(x0, y0, w, h, fill, edge="none", lw=0, kind=""):
    ax.add_patch(Rectangle((x0, y0), w, h, facecolor=fill, edgecolor=edge, lw=lw)); RECTS.append(dict(kind=kind, rect=[round(x0, 3), round(y0, 3), round(x0 + w, 3), round(y0 + h, 3)], fill=fill, edge=edge, lw=lw))
RC = dict(L.RC); RC.update({"figure.facecolor": "none", "savefig.facecolor": "none"})
with plt.rc_context(RC):
    t = OS["T90, % of the recording"]; text(t["x"], t["baseline"], t["text"], 11.0, "center", source_file="static:V13", source_key="V13 title", source_value=t["text"], rule="text")
    for i, x in enumerate(COL_X):
        rect(x, HY0, CW, HY1 - HY0, AZURE[i], kind=f"T90 header {COLS[i]}"); text(x + CW / 2, OS[COLS[i]]["baseline"], COLS[i], 10.0, "center", source_file="static:V13", source_key="V13 band label", source_value=COLS[i], rule="text")
    for ri, r in enumerate(R3):
        y0 = ROW_TOP0 + ri * PITCH; y1 = y0 + ROW_H; cy = (y0 + y1) / 2
        text(LAB_X1, cy + LAB_DY[0], ROWS[ri][0], 10.0, "right", source_file="static:V13", source_key="V13 row label", source_value=ROWS[ri][0], rule="text")
        text(LAB_X1, cy + LAB_DY[1], ROWS[ri][1], 10.0, "right", source_file="static:V13", source_key="V13 row label", source_value=ROWS[ri][1], rule="text")
        for ci, bn in enumerate(BN):
            n = int(r[f"n_{bn}"]); p = L.pct1(n, r["n"]); assert abs(float(p) - float(r[f"pct_{bn}"])) < 0.051, (r["category"], bn, p, r[f"pct_{bn}"])
            k = rung(float(p)); fill = LADDER[k]; col = "#ffffff" if k == 3 else INK; s = f"{n:,} ({p})"
            rect(COL_X[ci], y0, CW, ROW_H, fill, kind=f"cell {ROWS[ri][0]}|{COLS[ci]}")
            text(COL_X[ci] + CW / 2, cy + CELL_DY, s, 10.0, "center", color=col, source_file=SRC, source_key=f"rows3/{ri}", source_value=[n, r["n"]], rule=f"count_pct_{bn}")
            DR.append(dict(row=r["category"], row_label=ROWS[ri][0], col=COLS[ci], n=n, pct=p, pct_file=r[f"pct_{bn}"], group_n=r["n"], text=s, fill=fill, rung=RUNG_LAB[k], text_color=col,
                           rect=[COL_X[ci], round(y0, 3), COL_X[ci] + CW, round(y1, 3)], key=f"rows3[{ri}].n_{bn} / n", source=SRC))
    for k, r0 in enumerate(KEY_RECTS):
        x0, y0, x1, y1 = r0["rect"]; rect(x0, y0, x1 - x0, y1 - y0, LADDER[k], edge="white", lw=0.6, kind=f"key rung {RUNG_LAB[k]}")
        text(OS[RUNG_LAB[k]]["x"], OS[RUNG_LAB[k]]["baseline"], RUNG_LAB[k], 10.0, "left", source_file="static:V13", source_key="V13 key rung label", source_value=RUNG_LAB[k], rule="text")
    for kt in ("Percent of the", "apnea group"): text(OS[kt]["x"], OS[kt]["baseline"], kt, 10.0, "left", source_file="static:V13", source_key="V13 key title", source_value=kt, rule="text")
    BOTTOM_INK = max(s["baseline"] for s in STR); assert abs(BOTTOM_INK - G0["bottom_ink_baseline"]) < 1e-6, (BOTTOM_INK, G0["bottom_ink_baseline"])
    for t in ax.texts:   # round 49: every string inside the panel page (nothing clipped at the page edge)
        bb = t.get_window_extent(renderer=fig.canvas.get_renderer()); x0p, x1p = BX0 + bb.x0 / fig.dpi * 72, BX0 + bb.x1 / fig.dpi * 72
        assert x0p >= BX0 + 1.0 and x1p <= BX1 - 1.0, ("text leaves panel c's page", t.get_text(), x0p, x1p)
    fig.savefig(f"{WORK}/panel_c_raw.pdf", transparent=True); plt.close(fig)
# ---------------------------------------------------------------- proofs: geometry identical to V13 (strings at the V13 origins, sizes), values from the file
assert sorted((s["text"], round(s["size"] - DPT, 3), s["color"], s["ha"]) for s in STR if not s["rule"].startswith("count_pct")) == sorted((s["text"], s["size"], s["color"], s["ha"]) for s in G0["v13_strings"] if "(" not in s["text"] or "AHI" in s["text"])   # round 49: sizes are the V13 sizes + DPT
for s in STR:
    if s["rule"].startswith("count_pct"): continue
    assert abs(s["x"] - OS[s["text"]]["x"]) < 1e-6 and abs(s["baseline"] - OS[s["text"]]["baseline"]) < 1e-6, (s["text"], s["x"], OS[s["text"]]["x"])
old_by = {(o["row"], o["col"]): o for o in G0["v13_values"]}; assert len(old_by) == len(DR) == 12
CHANGES = []
for d in DR:
    o = old_by[(d["row"], d["col"])]
    if d["text"] != o["text"] or d["fill"] != o["fill"]: CHANGES.append(dict(row=d["row_label"], col=d["col"], before=o["text"], after=d["text"], rung_before=o["rung"], rung_after=d["rung"]))
geom = dict(G0, rects=RECTS, bottom_ink_baseline=BOTTOM_INK); geom.pop("v13_strings", None); geom.pop("v13_values", None)
json.dump(geom, open(f"{WORK}/panel_c_geometry.json", "w"), indent=1)
json.dump(dict(sheet=SHEET, panel="c", lane="V14_L4_APNEA", source=SRC, source_sha256=L.sha256(SRC), sidecar_step=SC["step"], v14_gate=SC["_v14_gate"], v13_drawn_record=G0["drawn_v13"],
               n=XT["n"], n_cohort=XT["n_cohort"], n_ahi_undefined=XT["n_ahi_undefined"], spearman=XT["spearman"], rows=[dict(category=r["category"], n=r["n"], t90_median=r["t90_median"]) for r in R3],
               values=DR, strings=STR, drawn_records=STR, changes_vs_v13=CHANGES, geometry=geom), open(f"{VER}/{SHEET}_c_drawn.json", "w"), indent=1, ensure_ascii=False)
print(f"panel c: 3 rows x 4 bands at the V13 geometry (rows {ROW_H:.3f} pt, grid bottom {GRID_BOTTOM:.3f}, lowest baseline {BOTTOM_INK:.3f}); n {XT['n']:,} of {XT['n_cohort']:,} ({XT['n_ahi_undefined']} without an AHI), Spearman {XT['spearman']}")
for c in CHANGES: print(f"  {c['row']:<11} {c['col']:<7} {c['before']:<14} -> {c['after']:<14} rung {c['rung_before']} -> {c['rung_after']}")
print(f"  {len(CHANGES)} of 12 cells change; wrote {WORK}/panel_c_raw.pdf")
