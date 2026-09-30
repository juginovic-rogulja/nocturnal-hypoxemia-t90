#!/usr/bin/env python3
"""Round 49 LED-B copy (2026-09-26): every text one point larger through text() and tw(), the labels that do not fit their column at +1 take the largest fitting size at or above their V26 size (SIZE_KEEP, declared), the star key not drawn (LNOTES round 40). Otherwise the round-38 copy. ED_Fig07 panel b (V14 lane L4, round 37): the nine-cell apnea-by-oxygen cross (numbers/alenfig3_cross_v1.json key cross9, three apnea
groups x three T90 levels, reference = No or mild apnea with T90 <=1%, the four conditions of Fig 4d as columns, BH stars within the
cell across the clinical conditions of cross9) from the v8.1 numbers, in EXACTLY the V13 design (the round-31 LF4 panel: Fig 4d's cell
size 58.067 x 40.154 pt on 4d's pitches, header rectangles, 9.5 pt Arial, the red ramp, the reference row white in a 0.9 pt ink box,
the key at 4d's offsets; every V13 string at its V13 origin, LF4/ED_Fig07/verify/ED_Fig07_b_drawn.json proved present on the V13 sheet).
Copy of ROUND31_2026-09-10/LF4_FIG4_ED7/ED_Fig07/scripts/01_build_b.py (beside as ed7_02_build_b_PRE_V8_1.py), repointed: v8.1 sidecar
gate, geometry from this lane's Main_Fig4 panel d record (= V13's), no round-30 identity assertions, values from the file.
Outputs: work/panel_b_raw.pdf, work/panel_b.pdf, work/panel_b_geometry.json, work/sheet_layout.json, verify/ED_Fig07_b_drawn.json."""
import json, math, os, re, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ed7spec as S
L = S.L
from statsmodels.stats.multitest import multipletests
matplotlib, plt = L.mpl_setup()
from matplotlib.patches import Rectangle
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties

SC = L.sidecar(S.SRC_B); assert SC["step"]["rc"] == 0 and str(SC["step"]["id"]) == "123", SC["step"]
X = json.load(open(L.hydrated(S.SRC_B))); C = X["cross9"]; CS = X["cell_sizes"]["cross9"]
XT = json.load(open(L.hydrated(f"{L.NUM}/crosstab_v2.json"))); L.sidecar(f"{L.NUM}/crosstab_v2.json")
GROUPS, LEVELS, REF = S.GROUPS, S.LEVELS, S.REF
assert X["ahi3_groups"]["rule"] == {"No or mild": "AHI < 15", "Moderate": "15 <= AHI < 30", "Severe": "AHI >= 30"}, X["ahi3_groups"]
CELLS = [f"{g}|{l}" for g in GROUPS for l in LEVELS]
assert sorted(CS) == sorted(CELLS) and sum(CS.values()) == X["provenance"]["cohort_n"] == sum(X["ahi3_groups"]["n"].values()) == XT["n"], (CS, X["provenance"]["cohort_n"], XT["n"])
CONDS = [c for c, v in C.items() if not v["negative_control"]]; CTRLS = [c for c, v in C.items() if v["negative_control"]]
COLS = S.COLS; assert all(c in CONDS for c in COLS), (COLS, CONDS)
for c in C:
    assert C[c]["reference"] == REF and list(C[c]["cells"].keys()) == CELLS and C[c]["cells"][REF]["reference"] and C[c]["cells"][REF]["hr"] == 1.0 and C[c]["cells"][REF]["p"] is None
Q = {c: {} for c in CONDS}
for cell in CELLS:
    if cell == REF: continue
    ps = [C[c]["cells"][cell]["p"] for c in CONDS]; assert all(p is not None for p in ps), cell
    for c, q in zip(CONDS, multipletests(ps, method="fdr_bh")[1]): Q[c][cell] = float(q)
def is_half(v): return bool(re.fullmatch(r"-?\d+\.\d\d5", repr(float(v))))

# ---------------------------------------------------------------- strings: the V13 sheet's own (LF4 record), each at its V13 origin
OLD = json.load(open(L.hydrated(S.LF4_B_DRAWN))); OS = {}
for s in OLD["strings"]: OS.setdefault(s["text"], s)
V13_POS = {(s["text"], round(s["x"], 3), round(s["baseline"], 3)) for s in OLD["strings"]}
GROUP_LAB, LEVEL_LAB, COL_LAB, LVL_A, LVL_O = S.GROUP_LAB, S.LEVEL_LAB, S.COL_LAB, S.LVL_A, S.LVL_O
REF_KEY, STAR_KEY, RAMP_TITLE, TICKS = S.REF_KEY, S.STAR_KEY_B, S.RAMP_TITLE, S.TICKS_B
ALL_STR = [REF_KEY, STAR_KEY, RAMP_TITLE, LVL_O, *LVL_A, *TICKS, *sum(GROUP_LAB.values(), []), *LEVEL_LAB.values(), *sum(COL_LAB.values(), [])]
for s in ALL_STR: assert s in OS, ("not a string of the V13 sheet", s); assert L.house_ok(s), s

# ---------------------------------------------------------------- geometry: 4d's cells, header rectangles and key (this lane's panel d record = V13's), transposed as LF4 did
D = json.load(open(L.hydrated(S.FIG4_D_GEOM))); DO = D["offsets"]; LAY13 = json.load(open(L.hydrated(S.LF4_LAYOUT)))
CELL_W, CELL_H, CELL_PITCH, ROW_PITCH = D["cell_w"], D["cell_h"], D["cell_pitch"], D["row_pitch"]
CELL_GAP_X = CELL_PITCH - CELL_W; ROW_GAP = ROW_PITCH - CELL_H; GROUP_GAP = D["group_pitch"] - D["group_w"]
HDR_GAP_1 = D["header2_y"][0] - D["header1_y"][1]; HDR_GAP_2 = D["row0_y"] - D["header2_y"][1]
FS = 9.5; KEY_FS = D["key"]["key_fs"]; INK = L.INK; LIGHT = "#ccd1d6"
ORANGE3 = ["#f6d6c2", "#e4864f", "#d55e00"]; ORANGE3_TEXT = [INK, INK, "#ffffff"]; AZURE3 = ["#b3dcf2", "#7cc0e9", "#0288d1"]
BT = json.load(open(L.hydrated(S.BASE_TEXT))); W, BASE_H = BT["page"]; MET = L.sheet_font_metrics(BT["spans"])
EDGE = S.EDGE; LETTER_X = S.LETTER_X; LETTER_BASELINE = S.LETTER_BASELINE; SHIFT_A = S.SHIFT_A
YB = LAY13["panel_b_rect"][1]; assert abs(YB - (SHIFT_A + 376.649 - S.PANEL_B_OVERLAP)) < 1e-9, YB
AX0 = LETTER_X; AW = 74.0; OX0 = AX0 + AW + HDR_GAP_1; OW = 40.0; GX0 = OX0 + OW + HDR_GAP_2
_FP = FontProperties(fname=L.ARIAL)
def tw(s, fs=FS): return float(TextPath((0, 0), s, size=SIZE_KEEP.get(s, fs + L.PT_PLUS), prop=_FP).get_extents().width)   # round 49: widths at the size the string is drawn
def tw_at(s, fs): return float(TextPath((0, 0), s, size=fs, prop=_FP).get_extents().width)
SIZE_KEEP = {}   # round 49: strings that do not fit their column at +1 pt take the largest size at or above their V26 size, in 0.5 pt steps, that fits (declared)
def fit_size(s, fs, maxw, what):
    for cand in (fs + L.PT_PLUS, fs + L.PT_PLUS / 2, fs):
        if tw_at(s, cand) <= maxw:
            if cand != fs + L.PT_PLUS: SIZE_KEEP[s] = cand
            return cand
    raise AssertionError((what, s, tw_at(s, fs), maxw))
for g in GROUPS:
    for ln in GROUP_LAB[g]: fit_size(ln, FS, AW - 2.0, "group label wider than its column")
for lv in LEVELS: fit_size(LEVEL_LAB[lv], FS, OW - 2.0, "level label wider than its column")
for ln in range(2):
    for a, b in zip(COLS[:-1], COLS[1:]): assert (tw(COL_LAB[a][ln]) + tw(COL_LAB[b][ln])) / 2 + 2.0 <= CELL_PITCH, ("column labels collide", COL_LAB[a][ln], COL_LAB[b][ln])
HB1 = 40.0; HLINE = 15.2; HEAD_PAD = 9.5
hb1 = YB + HB1; hb2 = hb1 + HLINE; GRID_TOP = hb2 + HEAD_PAD
GROUP_H = 3 * CELL_H + 2 * ROW_GAP; GROUP_PITCH = GROUP_H + GROUP_GAP
def row_top(i): return GRID_TOP + (i // 3) * GROUP_PITCH + (i % 3) * ROW_PITCH
GRID_BOTTOM = row_top(8) + CELL_H; GX1 = GX0 + 3 * CELL_PITCH + CELL_W
assert GX1 <= W - EDGE, (GX1, W - EDGE)
dgb = D["grid_bottom"]; DS = {s["text"]: s for s in json.load(open(L.hydrated(S.FIG4_D_DRAWN)))["strings"] if s["baseline"] > 1000}
kx0d, ky0d, kx1d, ky1d = D["key"]["ramp"]; rbx, rby, rbw, rbh = D["key"]["refbox"]
KX0 = AX0 + (kx0d - D["grid_x0"]); KX1 = KX0 + (kx1d - kx0d); KY0 = GRID_BOTTOM + (ky0d - dgb); KY1 = GRID_BOTTOM + (ky1d - dgb)
RBX = KX1 + (rbx - kx1d); RBY = GRID_BOTTOM + (rby - dgb); REF_TX = RBX + (DS[REF_KEY]["x"] - rbx); STAR_X = RBX + (DS[STAR_KEY]["x"] - rbx)
B_TITLE = GRID_BOTTOM + (DS[RAMP_TITLE]["baseline"] - dgb); B_TICK = GRID_BOTTOM + (DS["1"]["baseline"] - dgb); B_REF = GRID_BOTTOM + (DS[REF_KEY]["baseline"] - dgb); B_STAR = GRID_BOTTOM + (DS[STAR_KEY]["baseline"] - dgb)
assert REF_TX + tw(REF_KEY, KEY_FS) <= W - EDGE, ("reference key runs past the margin", REF_TX + tw(REF_KEY, KEY_FS))
BOTTOM_INK = max(B_TICK, B_STAR); SHEET_H = BOTTOM_INK + EDGE; HB_ = SHEET_H - YB
assert abs(SHEET_H - LAY13["sheet_h"]) < 0.01, ("the V13 sheet height must reproduce", SHEET_H, LAY13["sheet_h"])
RAMP_STOPS = [("1", 0.0), ("1.5", math.log(1.5) / math.log(3.25)), ("2", math.log(2) / math.log(3.25)), ("3", math.log(3) / math.log(3.25))]

fig = plt.figure(figsize=(W / 72, HB_ / 72)); fig.patch.set_alpha(0.0); ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(YB + HB_, YB); ax.axis("off"); ax.patch.set_visible(False)
DR = []; STR = []; RECTS = []; HALVES = []
ST = lambda t: dict(source_file="static:V13", source_key="V13 wording", source_value=t, rule="text")
def text(x, y, s, size, ha, color=INK, z=5, **src):
    size = SIZE_KEEP.get(s, size + L.PT_PLUS)   # round 49: +1 pt (or the declared fitted size)
    ax.text(x, y, s, fontsize=size, ha=ha, va="baseline", color=color, zorder=z); STR.append(dict(text=s, x=round(x, 3), baseline=round(y, 3), size=size, ha=ha, color=color, coords="page", panel="b", **src))
def rect(x0, y0, w, h, fill, edge="none", lw=0, z=None, kind=""):
    kw = dict(facecolor=fill, edgecolor=edge, lw=lw); kw.update({} if z is None else dict(zorder=z)); ax.add_patch(Rectangle((x0, y0), w, h, **kw))
    RECTS.append(dict(kind=kind, rect=[round(x0, 3), round(y0, 3), round(x0 + w, 3), round(y0 + h, 3)], fill=fill, edge=edge, lw=lw))
RC = dict(L.RC); RC.update({"figure.facecolor": "none", "savefig.facecolor": "none"})
with plt.rc_context(RC):
    text(AX0 + AW / 2, OS[LVL_A[0]]["baseline"], LVL_A[0], FS, "center", **ST(LVL_A[0])); text(AX0 + AW / 2, OS[LVL_A[1]]["baseline"], LVL_A[1], FS, "center", **ST(LVL_A[1]))
    text(OX0 + OW / 2, OS[LVL_O]["baseline"], LVL_O, FS, "center", **ST(LVL_O))
    for ci, c in enumerate(COLS):
        xc = GX0 + ci * CELL_PITCH + CELL_W / 2; l1, l2 = COL_LAB[c]; text(xc, hb1, l1, FS, "center", **ST(l1)); text(xc, hb2, l2, FS, "center", **ST(l2))
    for gi, g in enumerate(GROUPS):
        gy0 = row_top(3 * gi); gc = gy0 + GROUP_H / 2; col = ORANGE3_TEXT[gi]
        rect(AX0, gy0, AW, GROUP_H, ORANGE3[gi], kind=f"group header {g}")
        l1, l2 = GROUP_LAB[g]; text(AX0 + AW / 2, gc + DO["G_TXT_DY"] - 5.4, l1, FS, "center", color=col, **ST(l1)); text(AX0 + AW / 2, gc + DO["G_TXT_DY"] + 5.4, l2, FS, "center", color=col, **ST(l2))
        for li, lv in enumerate(LEVELS):
            i = 3 * gi + li; ry0 = row_top(i); cy = ry0 + CELL_H / 2; cell = f"{g}|{lv}"
            rect(OX0, ry0, OW, CELL_H, AZURE3[li], kind=f"level header {cell}"); text(OX0 + OW / 2, cy + DO["LAB1_DY"], LEVEL_LAB[lv], FS, "center", **ST(LEVEL_LAB[lv]))
            for ci, cond in enumerate(COLS):
                v = C[cond]["cells"][cell]; x0 = GX0 + ci * CELL_PITCH; xc = x0 + CELL_W / 2
                if cell == REF:
                    text(xc, ry0 + DO["REF_DY"] + DO["DYC"], "1.00", FS, "center", z=6, source_file=S.SRC_B, source_key=f"cross9/{cond}/cells/{cell}/hr", source_value=1.0, rule="2dp_halfup")
                    DR.append(dict(row=cell, column=cond, hr=1.0, lo=None, hi=None, p=None, q=None, n=v["n"], events=v["events"], hr_text="1.00", ci_text=None, stars="", fill="#ffffff", text_color=INK,
                                   rect=[round(x0, 3), round(ry0, 3), round(x0 + CELL_W, 3), round(ry0 + CELL_H, 3)], key=f"cross9[{cond}].cells[{cell}]", q_rule="reference cell, hazard ratio 1 by definition", source=S.SRC_B)); continue
                q = Q[cond][cell]; st = L.stars(q); fill = L.ramp_hex(v["hr"]); col = L.text_on(fill); hr_s = L.r2(v["hr"]) + st; ci_s = f"({L.r2(v['lo'])}-{L.r2(v['hi'])})"
                for nm in ("hr", "lo", "hi"):
                    if is_half(v[nm]): HALVES.append(dict(column=cond, cell=cell, field=nm, stored=v[nm], printed=L.r2(v[nm])))
                rect(x0, ry0, CELL_W, CELL_H, fill, kind=f"cell {cell}|{cond}")
                text(xc, ry0 + DO["HR_DY"] + DO["DYC"], hr_s, FS, "center", color=col, source_file=S.SRC_B, source_key=f"cross9/{cond}/cells/{cell}", source_value=[v["hr"], st], rule="cross9_hr_stars")
                text(xc, ry0 + DO["CI_DY"] + DO["DYC"], ci_s, FS, "center", color=col, source_file=S.SRC_B, source_key=f"cross9/{cond}/cells/{cell}", source_value=[v["lo"], v["hi"]], rule="cell_ci_paren")
                DR.append(dict(row=cell, column=cond, hr=v["hr"], lo=v["lo"], hi=v["hi"], p=v["p"], q=q, n=v["n"], events=v["events"], hr_text=hr_s, ci_text=ci_s, stars=st, fill=fill, text_color=col,
                               rect=[round(x0, 3), round(ry0, 3), round(x0 + CELL_W, 3), round(ry0 + CELL_H, 3)], key=f"cross9[{cond}].cells[{cell}]",
                               q_rule=f"Benjamini-Hochberg across the {len(CONDS)} clinical conditions of cross9 within the cell", source=S.SRC_B))
    rect(GX0, row_top(0), GX1 - GX0, CELL_H, "white", edge=INK, lw=0.9, z=3, kind="reference row box")
    n = 120
    for k in range(n):
        h0 = math.exp((k / n) * math.log(3.25)); h1 = math.exp(((k + 1) / n) * math.log(3.25))
        ax.add_patch(Rectangle((KX0 + k * (KX1 - KX0) / n, KY0), (KX1 - KX0) / n + 0.4, KY1 - KY0, facecolor=L.ramp_hex(math.sqrt(h0 * h1)), edgecolor="none"))
    rect(KX0, KY0, KX1 - KX0, KY1 - KY0, "none", edge=LIGHT, lw=0.4, kind="ramp frame")
    text((KX0 + KX1) / 2, B_TITLE, RAMP_TITLE, KEY_FS, "center", **ST(RAMP_TITLE))
    for lab, frac in RAMP_STOPS: text(KX0 + frac * (KX1 - KX0), B_TICK, lab, KEY_FS, "center", **ST(lab))
    rect(RBX, RBY, rbw, rbh, "white", edge=INK, lw=0.9, kind="reference key box")
    text(REF_TX, B_REF, REF_KEY, KEY_FS, "left", **ST(REF_KEY))   # round 49: the star key is NOT drawn (removed by lane LNOTES, round 40)
    for t in ax.texts: assert t.get_fontsize() >= 9.0
    fig.savefig(S.PANEL_B_RAW, transparent=True); plt.close(fig)
done = L.align_font_metrics(S.PANEL_B_RAW, S.PANEL_B, MET); assert "ArialMT" in done, done
# ---------------------------------------------------------------- proofs: values against the file; every static string at its V13 origin; read-back on the panel's text layer
for d in DR:
    v = C[d["column"]]["cells"][d["row"]]
    if d["row"] == REF: assert d["hr_text"] == "1.00" and v["hr"] == 1.0; continue
    assert d["hr_text"].rstrip("*") == L.r2(v["hr"]) and d["ci_text"] == f"({L.r2(v['lo'])}-{L.r2(v['hi'])})" and d["stars"] == L.stars(d["q"]) and d["fill"] == L.ramp_hex(v["hr"]), d
for s in STR:
    if s["source_file"] == "static:V13": assert (s["text"], round(s["x"], 3), round(s["baseline"], 3)) in V13_POS, ("static string not at a V13 origin", s["text"], s["x"], s["baseline"])
RB = L.spans_of(S.PANEL_B)
for s in STR:
    h = [r for r in RB if r["text"].replace("\xa0", " ") == s["text"] and abs(r["origin"][1] + YB - s["baseline"]) < 0.05]; assert h, ("string not read back at its baseline", s)
old_by = {(o["row"], o["column"]): o for o in OLD["values"]}; CHANGES = []
for d in DR:
    o = old_by.get((d["row"], d["column"]))
    if o is None or d["hr_text"] != o["hr_text"] or d["ci_text"] != o["ci_text"] or d["fill"] != o["fill"]:
        CHANGES.append(dict(row=d["row"], column=d["column"], before=(f"{o['hr_text']} {o['ci_text']}" if o else None), after=f"{d['hr_text']} {d['ci_text']}", q_before=(o["q"] if o else None), q_after=d["q"]))
geom = dict(panel_box_page=[0, YB, W, SHEET_H], panel_h=HB_, sheet_w=W, sheet_h=SHEET_H, header_baselines=[hb1, hb2], grid=dict(top=GRID_TOP, bottom=GRID_BOTTOM, gx0=GX0, gx1=GX1, cell_w=CELL_W, cell_h=CELL_H, cell_pitch=CELL_PITCH, row_pitch=ROW_PITCH,
            group_h=GROUP_H, group_pitch=GROUP_PITCH, row_tops=[round(row_top(i), 3) for i in range(9)], group_col=[AX0, AX0 + AW], level_col=[OX0, OX0 + OW], gaps=dict(cell_x=CELL_GAP_X, row=ROW_GAP, group=GROUP_GAP, hdr1=HDR_GAP_1, hdr2=HDR_GAP_2)),
            key=dict(ramp=[KX0, KY0, KX1, KY1], refbox=[RBX, RBY, rbw, rbh], ref_text_x=REF_TX, star_x=STAR_X, baselines=dict(title=B_TITLE, ticks=B_TICK, ref=B_REF, stars=B_STAR), key_fs=KEY_FS),
            bottom_ink=BOTTOM_INK, edge=EDGE, rects=RECTS, font_metrics_aligned=done, offsets=DO, fs=FS, pt_plus=L.PT_PLUS, size_exceptions=SIZE_KEEP, star_key_drawn=False, design_source=S.FIG4_D_GEOM)
json.dump(geom, open(S.GEOM_B, "w"), indent=1)
json.dump(dict(sheet_w=W, sheet_h=SHEET_H, base_h=BASE_H, panel_a_rect=[0, SHIFT_A, W, SHIFT_A + 376.649], panel_b_rect=[0, YB, W, SHEET_H], shift_a=SHIFT_A,
               letters=[["a", LETTER_X, LETTER_BASELINE], ["b", LETTER_X, YB + LETTER_BASELINE]], letter_size=S.LETTER_SIZE, edge=EDGE), open(S.LAYOUT, "w"), indent=1)
json.dump(dict(sheet=S.SHEET, panel="b", lane="V14_L4_APNEA", source=S.SRC_B, source_sha256=L.sha256(S.SRC_B), sidecar_step=SC["step"], v14_gate=SC["_v14_gate"], v13_drawn_record=S.LF4_B_DRAWN,
               cross="cross9", reference=REF, groups=GROUPS, levels=LEVELS, ahi_rule=X["ahi3_groups"]["rule"], group_n=X["ahi3_groups"]["n"], cell_sizes=CS, cohort_n=X["provenance"]["cohort_n"], columns=COLS, clinical_pool=CONDS, control_pool_not_drawn=CTRLS,
               star_rule=f"Benjamini-Hochberg across the {len(CONDS)} clinical conditions of cross9 within the cell", values=DR, strings=STR, drawn_records=STR, changes_vs_v13=CHANGES, exact_half_strings=HALVES,
               stored_precision="3 dp in the numbers file, printed half up on the stored value", geometry=geom, text_layer_panel=RB), open(S.DRAWN_B, "w"), indent=1, ensure_ascii=False)
print(f"panel b: 9 rows x 4 columns at the V13 geometry (cells {CELL_W:.3f} x {CELL_H:.3f}, grid y {GRID_TOP:.2f}..{GRID_BOTTOM:.2f}); sheet {W:.2f} x {SHEET_H:.2f} pt (V13 {LAY13['sheet_h']:.2f}); cells {CS}")
for d in DR: print(f"  {d['row']:<20} {d['column']:<26} {d['hr_text']:<9} {str(d['ci_text']):<13} q={None if d['q'] is None else round(d['q'], 4)} n {d['n']} ev {d['events']}")
print(f"  {len(CHANGES)} of 36 cells change vs V13; exact halves {len(HALVES)}: {[(h['column'], h['cell'], h['field'], h['stored']) for h in HALVES]}; wrote {S.PANEL_B}")
