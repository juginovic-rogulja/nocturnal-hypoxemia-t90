#!/usr/bin/env python3
"""ROUND 49 (2026-09-26, lane L4): copy of the round-38 builder with every text one point larger (DPT = 1.0 inside text(): 9.5 -> 10.5,
the key 9.087 -> 10.087); the panel page widened 10 pt to the left (PAGE_PAD_LEFT) so the right-aligned row labels are not clipped;
positions, cells, colours and strings unchanged. Round-38 docstring follows.
V14 lane L4 (round 37), Main Fig 4 panel d: the twin of Fig 3c from numbers/alenfig3_cross_v1.json key cross6 (three apnea groups x
two oxygen levels, reference = No or mild apnea with T90 <=1%, the (1, 10] nights set aside) from the v8.1 numbers, at EXACTLY the
V13 geometry (the round-31 LF4 wide panel: grid to the right margin 956.621, cells 58.067 x 40.154 pt, header rectangles, key at its
offsets; LF4/Main_Fig4/work/panel_d_geometry.json, proved to be the V13 panel by 00_probe_v13.py). Rows: the owner's round-30 pick
(Cardiovascular composite, Type 2 diabetes, Acute kidney injury, Heart failure) plus ONE negative control: back pain on V13, which
left the control panel on 14 September (straight swap back pain -> alopecia, V8_1_DECISIONS item 2), so the control row is Alopecia.
Stars: Benjamini-Hochberg q within each cell across the clinical conditions of cross6 (the sheet's rule), the control row across the
five controls in the cell. Copy of ROUND31_2026-09-10/LF4_FIG4_ED7/Main_Fig4/scripts/01_build_d.py (beside as 01_build_d_PRE_V8_1.py),
repointed: v8.1 sidecar gate, no round-30 identity assertions, positions pinned to the V13 drawn record, values from the file.
Also writes verify/FIG6_BH_CHECK_v8_1.json (BH across the severe-with-preserved column, the Fig 6 claim check).
Outputs: work/panel_d_raw.pdf, work/panel_d_geometry.json, verify/Main_Fig4_d_drawn.json."""
import json, math, os, re, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import l4lib as L
from statsmodels.stats.multitest import multipletests
matplotlib, plt = L.mpl_setup()
from matplotlib.patches import Rectangle

SHEET = "Main_Fig4"; SD = f"{L.LANE}/{SHEET}"; WORK = f"{SD}/work"; VER = f"{SD}/verify"
DPT = 1.0             # round 49: every text element one point larger, applied inside text()
PAGE_PAD_LEFT = 10.0  # round 49: the panel PAGE extends 10 pt further left (transparent) so the right-aligned row labels at 10.5 pt
                      # ("(negative control)" 81.2 pt, "Type 2 diabetes" 74.2 pt, right edge pinned at 571.141) are not clipped at the
                      # V13 page edge 497.2; no element moves (every string and rectangle keeps its page coordinates)
SRC = f"{L.NUM}/alenfig3_cross_v1.json"; SC = L.sidecar(SRC); assert SC["step"]["rc"] == 0, SC["step"]
NCF = json.load(open(L.hydrated(f"{L.NUM}/negcontrols_final.json"))); L.sidecar(f"{L.NUM}/negcontrols_final.json")
XT = json.load(open(L.hydrated(f"{L.NUM}/crosstab_v2.json"))); L.sidecar(f"{L.NUM}/crosstab_v2.json")
X = json.load(open(L.hydrated(SRC))); C = X["cross6"]; CS = X["cell_sizes"]["cross6"]
GROUPS = ["No or mild", "Moderate", "Severe"]; LEVELS = ["<=1%", ">10%"]; REF = "No or mild|<=1%"
assert X["ahi3_groups"]["rule"] == {"No or mild": "AHI < 15", "Moderate": "15 <= AHI < 30", "Severe": "AHI >= 30"}, X["ahi3_groups"]
CONDS = [c for c, v in C.items() if not v["negative_control"]]; CTRLS = [c for c, v in C.items() if v["negative_control"]]
assert set(CTRLS) == set(NCF["panel_labels"]) and len(CTRLS) == 5, (CTRLS, NCF["panel_labels"])
SWAP = NCF.get("decision_record", ""); CONTROL_ROW = "Alopecia"
assert CONTROL_ROW in CTRLS, (CONTROL_ROW, CTRLS)
ROWS = ["Cardiovascular composite", "Type 2 diabetes", "Acute kidney injury", "Heart failure", CONTROL_ROW]
assert all(r in C for r in ROWS) and C[CONTROL_ROW]["negative_control"] and all(not C[r]["negative_control"] for r in ROWS[:4])
for r in ROWS: assert C[r]["reference"] == REF and C[r]["cells"][REF]["reference"] and C[r]["cells"][REF]["hr"] == 1.0
CELLS = [f"{g}|{l}" for g in GROUPS for l in LEVELS]; assert list(C[ROWS[0]]["cells"].keys()) == CELLS
SET_ASIDE = X["cell_sizes"]["cross6_set_aside_t90_gt1_le10"]; assert sum(CS.values()) + SET_ASIDE == XT["n"] == X["provenance"]["cohort_n"], (sum(CS.values()), SET_ASIDE, XT["n"])
for r in ROWS:
    for cell in CELLS:
        if cell != REF: assert C[r]["cells"][cell]["p"] is not None, (r, cell, "cell not fitted")
Q = {r: {} for r in C}
for cell in CELLS:
    if cell == REF: continue
    ps = [C[c]["cells"][cell]["p"] for c in CONDS]; assert all(p is not None for p in ps), cell
    for c, q in zip(CONDS, multipletests(ps, method="fdr_bh")[1]): Q[c][cell] = float(q)
    pc = [C[c]["cells"][cell]["p"] for c in CTRLS]
    for c, q in zip(CTRLS, multipletests(pc, method="fdr_bh")[1]): Q[c][cell] = float(q)
# ---------------------------------------------------------------- the Fig 6 claim check: BH across the severe-with-preserved column
SEV_PRES = "Severe|<=1%"
bh6 = dict(cell=SEV_PRES, reference=REF, rule="Benjamini-Hochberg across the clinical conditions of cross6 within the cell (the sheet's rule)", source=SRC, source_sha256=L.sha256(SRC),
           conditions={c: dict(hr=C[c]["cells"][SEV_PRES]["hr"], lo=C[c]["cells"][SEV_PRES]["lo"], hi=C[c]["cells"][SEV_PRES]["hi"], p=C[c]["cells"][SEV_PRES]["p"], q=Q[c][SEV_PRES], n=C[c]["cells"][SEV_PRES]["n"], events=C[c]["cells"][SEV_PRES]["events"]) for c in CONDS},
           controls={c: dict(hr=C[c]["cells"][SEV_PRES]["hr"], p=C[c]["cells"][SEV_PRES]["p"], q=Q[c][SEV_PRES]) for c in CTRLS})
bh6["survivors_q_lt_0_05"] = sorted([c for c in CONDS if Q[c][SEV_PRES] < 0.05], key=lambda c: Q[c][SEV_PRES]); bh6["raw_p_lt_0_05"] = sorted(c for c in CONDS if C[c]["cells"][SEV_PRES]["p"] < 0.05)
bh6["verdict"] = ("sleep apnea alone (severe apnea with preserved oxygen): " + (f"{len(bh6['survivors_q_lt_0_05'])} of {len(CONDS)} conditions clear BH q < 0.05: {bh6['survivors_q_lt_0_05']}" if bh6["survivors_q_lt_0_05"] else f"no condition of {len(CONDS)} clears BH q < 0.05"))
json.dump(bh6, open(f"{VER}/FIG6_BH_CHECK_v8_1.json", "w"), indent=1)

# ---------------------------------------------------------------- geometry: the V13 panel's (LF4), every string at its V13 origin
G13 = json.load(open(f"{WORK}/v13_geometry.json")); G0 = G13["panel_d"]; OS = {}
for s in G0["v13_strings"]: OS.setdefault(s["text"], s)
BX0, BY0, BX1, BY1 = G0["panel_box"]; V13_BX0 = BX0; BX0 = BX0 - PAGE_PAD_LEFT; O = G0["offsets"]; GX0 = G0["grid_x0"]; GRID_X1 = G0["grid_x1"]
CELL_W, CELL_PITCH, GRP_W, GRP_PITCH, CELL_H, ROW_PITCH = G0["cell_w"], G0["cell_pitch"], G0["group_w"], G0["group_pitch"], G0["cell_h"], G0["row_pitch"]; DYC = O["DYC"]
H1Y0, H1Y1 = G0["header1_y"]; H2Y0, H2Y1 = G0["header2_y"]; ROW0 = G0["row0_y"]; LAB_X1 = G0["label_x1"]; GRID_BOTTOM = G0["grid_bottom"]
assert abs(ROW0 + (len(ROWS) - 1) * ROW_PITCH + CELL_H - GRID_BOTTOM) < 1e-6 and abs(GX0 + 2 * GRP_PITCH + GRP_W - GRID_X1) < 1e-6
K = G0["key"]; RAMP = K["ramp"]; REFBOX = K["refbox"]; KEY_FS = K["key_fs"]; RAMP_STOPS = [("1", 0.0), ("1.5", math.log(1.5) / math.log(3.25)), ("2", math.log(2) / math.log(3.25)), ("3", math.log(3) / math.log(3.25))]
INK = L.INK; ORANGE_HDR = ["#f6d6c2", "#e4864f", "#d55e00"]; BLUE_HDR = ["#b3dcf2", "#0288d1"]
GROUP_LAB = {"No or mild": ["No or mild", "(AHI <15)"], "Moderate": ["Moderate", "(AHI 15 to <30)"], "Severe": ["Severe", "(AHI 30 or more)"]}
LEVEL_LAB = {"<=1%": ["Preserved", "T90 ≤1%"], ">10%": ["Low", "T90 >10%"]}
ROW_LAB = {"Cardiovascular composite": ["Cardiovascular", "composite"], "Type 2 diabetes": ["Type 2 diabetes"], "Acute kidney injury": ["Acute kidney", "injury"], "Heart failure": ["Heart failure"], CONTROL_ROW: [CONTROL_ROW, "(negative control)"]}
REF_KEY = "Reference, no or mild apnea with T90 ≤1%"; STAR_KEY = "* q < 0.05, ** q < 0.01, *** q < 0.001"; RAMP_TITLE = "Hazard ratio"
for s in (REF_KEY, STAR_KEY, RAMP_TITLE, *sum(GROUP_LAB.values(), []), *sum(LEVEL_LAB.values(), []), *sum(ROW_LAB.values(), [])): assert L.house_ok(s), s
for t in (REF_KEY, STAR_KEY, RAMP_TITLE, "1", "1.5", "2", "3", "Sleep apnea", "Oxygenation", *sum(GROUP_LAB.values(), []), *sum(LEVEL_LAB.values(), []), "(negative control)"): assert t in OS, ("V13 string missing", t)
def is_half(v): return bool(re.fullmatch(r"-?\d+\.\d\d5", repr(float(v))))

W, H = BX1 - BX0, BY1 - BY0
fig = plt.figure(figsize=(W / 72, H / 72)); fig.patch.set_alpha(0.0); ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(BX0, BX1); ax.set_ylim(BY1, BY0); ax.axis("off"); ax.patch.set_visible(False)
DR = []; STR = []; RECTS = []; HALVES = []
def text(x, y, s, size, ha, color=INK, z=5, **src):
    ax.text(x, y, s, fontsize=size + DPT, ha=ha, va="baseline", color=color, zorder=z); STR.append(dict(text=s, x=round(x, 3), baseline=round(y, 3), size=size + DPT, ha=ha, color=color, panel="d", **src))   # round 49: +DPT
def rect(x0, y0, w, h, fill, edge="none", lw=0, z=None, kind=""):
    kw = dict(facecolor=fill, edgecolor=edge, lw=lw); kw.update({} if z is None else dict(zorder=z)); ax.add_patch(Rectangle((x0, y0), w, h, **kw))
    RECTS.append(dict(kind=kind, rect=[round(x0, 3), round(y0, 3), round(x0 + w, 3), round(y0 + h, 3)], fill=fill, edge=edge, lw=lw))
ST = lambda t: dict(source_file="static:V13", source_key="V13 wording", source_value=t, rule="text")
RC = dict(L.RC); RC.update({"figure.facecolor": "none", "savefig.facecolor": "none"})
with plt.rc_context(RC):
    h1c = (H1Y0 + H1Y1) / 2; h2c = (H2Y0 + H2Y1) / 2
    text(LAB_X1, h1c + O["H1_TXT_DY"], "Sleep apnea", 9.5, "right", **ST("Sleep apnea")); text(LAB_X1, h2c + O["H2_TXT_DY"], "Oxygenation", 9.5, "right", **ST("Oxygenation"))
    for gi, g in enumerate(GROUPS):
        gx = GX0 + gi * GRP_PITCH; col = "#ffffff" if gi == 2 else INK
        rect(gx, H1Y0, GRP_W, H1Y1 - H1Y0, ORANGE_HDR[gi], kind=f"group header {g}")
        l1, l2 = GROUP_LAB[g]; text(gx + GRP_W / 2, h1c + O["G_TXT_DY"] - 5.4, l1, 9.5, "center", color=col, **ST(l1)); text(gx + GRP_W / 2, h1c + O["G_TXT_DY"] + 5.4, l2, 9.5, "center", color=col, **ST(l2))
        for li, lv in enumerate(LEVELS):
            cx0 = gx + li * CELL_PITCH; rect(cx0, H2Y0, CELL_W, H2Y1 - H2Y0, BLUE_HDR[li], kind=f"level header {g}|{lv}")
            a, b = LEVEL_LAB[lv]; text(cx0 + CELL_W / 2, h2c + O["L1_DY"], a, 9.5, "center", **ST(a)); text(cx0 + CELL_W / 2, h2c + O["L2_DY"], b, 9.5, "center", **ST(b))
    for ri, r in enumerate(ROWS):
        y0 = ROW0 + ri * ROW_PITCH; y1 = y0 + CELL_H; cy = (y0 + y1) / 2; labs = ROW_LAB[r]
        # row labels are the file's own keys (data-driven), recorded with the key name; a two-line label records each line
        if len(labs) == 1: text(LAB_X1, cy + O["LAB1_DY"], labs[0], 9.5, "right", source_file=SRC, source_key="", source_value=labs[0], rule="text", note=f"row label = key of cross6 ({r})")
        else:
            text(LAB_X1, cy + O["LAB2_DY"][0], labs[0], 9.5, "right", source_file=SRC, source_key="", source_value=labs[0], rule="text", note=f"line 1 of the row label, key of cross6 ({r})")
            text(LAB_X1, cy + O["LAB2_DY"][1], labs[1], 9.5, "right", **(ST(labs[1]) if labs[1] == "(negative control)" else dict(source_file=SRC, source_key="", source_value=labs[1], rule="text", note=f"line 2 of the row label, key of cross6 ({r})")))
        for gi, g in enumerate(GROUPS):
            for li, lv in enumerate(LEVELS):
                cell = f"{g}|{lv}"; v = C[r]["cells"][cell]; cx0 = GX0 + gi * GRP_PITCH + li * CELL_PITCH; cxc = cx0 + CELL_W / 2
                if cell == REF:
                    text(cxc, y0 + O["REF_DY"] + DYC, "1.00", 9.5, "center", z=6, source_file=SRC, source_key=f"cross6/{r}/cells/{cell}/hr", source_value=1.0, rule="2dp_halfup")
                    DR.append(dict(row=r, cell=cell, hr=1.0, lo=None, hi=None, p=None, q=None, n=v["n"], events=v["events"], hr_text="1.00", ci_text=None, stars="", fill="#ffffff", text_color=INK,
                                   rect=[round(cx0, 3), round(y0, 3), round(cx0 + CELL_W, 3), round(y1, 3)], key=f"cross6[{r}].cells[{cell}]", source=SRC)); continue
                q = Q[r][cell]; st = L.stars(q); fill = L.ramp_hex(v["hr"]); col = L.text_on(fill); hr_s = L.r2(v["hr"]) + st; ci_s = f"({L.r2(v['lo'])}-{L.r2(v['hi'])})"
                for nm in ("hr", "lo", "hi"):
                    if is_half(v[nm]): HALVES.append(dict(row=r, cell=cell, field=nm, stored=v[nm], printed=L.r2(v[nm])))
                rect(cx0, y0, CELL_W, CELL_H, fill, kind=f"cell {r}|{cell}")
                text(cxc, y0 + O["HR_DY"] + DYC, hr_s, 9.5, "center", color=col, source_file=SRC, source_key=f"cross6/{r}/cells/{cell}", source_value=[v["hr"], st], rule="cross6_hr_stars", stars=st, q=q)
                text(cxc, y0 + O["CI_DY"] + DYC, ci_s, 9.5, "center", color=col, source_file=SRC, source_key=f"cross6/{r}/cells/{cell}", source_value=[v["lo"], v["hi"]], rule="cell_ci_paren")
                DR.append(dict(row=r, cell=cell, hr=v["hr"], lo=v["lo"], hi=v["hi"], p=v["p"], q=q, n=v["n"], events=v["events"], hr_text=hr_s, ci_text=ci_s, stars=st, fill=fill, text_color=col,
                               rect=[round(cx0, 3), round(y0, 3), round(cx0 + CELL_W, 3), round(y1, 3)], key=f"cross6[{r}].cells[{cell}]",
                               q_rule=(f"BH across the {len(CONDS)} clinical conditions in the cell" if r != CONTROL_ROW else "BH across the five negative controls in the cell"), source=SRC))
    rect(GX0, ROW0, CELL_W, GRID_BOTTOM - ROW0, "white", edge=INK, lw=0.9, z=3, kind="reference column box")
    n = 120; kx0, ky0, kx1, ky1 = RAMP
    for k in range(n):
        h0 = math.exp((k / n) * math.log(3.25)); h1 = math.exp(((k + 1) / n) * math.log(3.25))
        ax.add_patch(Rectangle((kx0 + k * (kx1 - kx0) / n, ky0), (kx1 - kx0) / n + 0.4, ky1 - ky0, facecolor=L.ramp_hex(math.sqrt(h0 * h1)), edgecolor="none"))
    rect(kx0, ky0, kx1 - kx0, ky1 - ky0, "none", edge="#ccd1d6", lw=0.4, kind="ramp frame")
    text((kx0 + kx1) / 2, OS[RAMP_TITLE]["baseline"], RAMP_TITLE, KEY_FS, "center", **ST(RAMP_TITLE))
    for lab, frac in RAMP_STOPS: text(kx0 + frac * (kx1 - kx0), OS[lab]["baseline"], lab, KEY_FS, "center", **ST(lab))
    rect(REFBOX[0], REFBOX[1], REFBOX[2], REFBOX[3], "white", edge=INK, lw=0.9, kind="reference key box")
    text(OS[REF_KEY]["x"], OS[REF_KEY]["baseline"], REF_KEY, KEY_FS, "left", **ST(REF_KEY)); text(OS[STAR_KEY]["x"], OS[STAR_KEY]["baseline"], STAR_KEY, KEY_FS, "left", **ST(STAR_KEY))
    BOTTOM_INK = max(s["baseline"] for s in STR); assert abs(BOTTOM_INK - G0["bottom_ink_baseline"]) < 1e-6 and BOTTOM_INK + 4 <= BY1, (BOTTOM_INK, G0["bottom_ink_baseline"])
    LEFTMOST = None
    for t in ax.texts:   # round 49: every string inside the (widened) panel page, nothing clipped at the page edge
        bb = t.get_window_extent(renderer=fig.canvas.get_renderer()); x0p, x1p = BX0 + bb.x0 / fig.dpi * 72, BX0 + bb.x1 / fig.dpi * 72
        assert x0p >= BX0 + 1.0 and x1p <= BX1 - 1.0, ("text leaves panel d's page", t.get_text(), x0p, x1p)
        if LEFTMOST is None or x0p < LEFTMOST[0]: LEFTMOST = (round(x0p, 2), t.get_text())
    fig.savefig(f"{WORK}/panel_d_raw.pdf", transparent=True); plt.close(fig)
# ---------------------------------------------------------------- proofs: values against the file; static strings at the V13 origins
for d in DR:
    v = C[d["row"]]["cells"][d["cell"]]
    if d["cell"] == REF: assert d["hr_text"] == "1.00" and v["hr"] == 1.0; continue
    assert d["hr_text"].rstrip("*") == L.r2(v["hr"]) and d["ci_text"] == f"({L.r2(v['lo'])}-{L.r2(v['hi'])})" and d["stars"] == L.stars(d["q"]) and d["fill"] == L.ramp_hex(v["hr"]), d
V13_POS = {(s["text"], round(s["x"], 3), round(s["baseline"], 3)) for s in G0["v13_strings"]}
for s in STR:
    if s["source_file"] == "static:V13" and s["text"] != "(negative control)":
        assert (s["text"], round(s["x"], 3), round(s["baseline"], 3)) in V13_POS, ("static string not at a V13 origin", s["text"], s["x"], s["baseline"])
old_by = {(o["row"], o["cell"]): o for o in G0["v13_values"]}; CHANGES = []
for d in DR:
    key_old = ("Back pain", d["cell"]) if d["row"] == CONTROL_ROW else (d["row"], d["cell"]); o = old_by.get(key_old)
    if o is None or d["hr_text"] != o["hr_text"] or d["ci_text"] != o["ci_text"] or d["fill"] != o["fill"]:
        CHANGES.append(dict(row=d["row"], row_v13=key_old[0], cell=d["cell"], before=(f"{o['hr_text']} {o['ci_text']}" if o else None), after=f"{d['hr_text']} {d['ci_text']}", q_before=(o["q"] if o else None), q_after=d["q"]))
geom = dict({k: v for k, v in G0.items() if k not in ("v13_strings", "v13_values", "rects", "v7_geometry")}, rects=RECTS, bottom_ink_baseline=BOTTOM_INK,
            page_box=[BX0, BY0, BX1, BY1], page_pad_left=PAGE_PAD_LEFT, v13_panel_box=[V13_BX0, BY0, BX1, BY1], leftmost_text=LEFTMOST)   # round 49: the page is the V13 box widened left by PAGE_PAD_LEFT
json.dump(geom, open(f"{WORK}/panel_d_geometry.json", "w"), indent=1)
json.dump(dict(sheet=SHEET, panel="d", lane="V14_L4_APNEA", source=SRC, source_sha256=L.sha256(SRC), sidecar_step=SC["step"], v14_gate=SC["_v14_gate"], v13_drawn_record=G0["drawn_v13"],
               cross="cross6", reference=REF, cell_sizes=CS, set_aside_t90_gt1_le10=SET_ASIDE, cohort_n=X["provenance"]["cohort_n"], rows=ROWS, control_row=CONTROL_ROW, control_row_v13="Back pain",
               control_swap=SWAP, clinical_pool=CONDS, control_pool=CTRLS, values=DR, strings=STR, drawn_records=STR, changes_vs_v13=CHANGES, exact_half_strings=HALVES,
               stored_precision="3 dp in the numbers file, printed half up on the stored value", fig6_bh_check=bh6["verdict"], geometry=geom), open(f"{VER}/{SHEET}_d_drawn.json", "w"), indent=1, ensure_ascii=False)
print(f"panel d: {len(ROWS)} rows x 6 cells at the V13 geometry (grid x {GX0:.2f}..{GRID_X1:.2f}, rows {ROW0:.2f}..{GRID_BOTTOM:.2f}); cells {CS}, set aside {SET_ASIDE}; control row {CONTROL_ROW} (V13: Back pain)")
for d in DR: print(f"  {d['row']:<26} {d['cell']:<18} {d['hr_text']:<9} {str(d['ci_text']):<13} q={None if d['q'] is None else round(d['q'], 4)} n {d['n']} ev {d['events']} fill {d['fill']} text {d['text_color']}")
print(f"  {len(CHANGES)} of 30 cells change vs V13; exact halves {len(HALVES)}: {[(h['row'], h['cell'], h['field'], h['stored']) for h in HALVES]}")
print(f"  FIG 6 CHECK ({SEV_PRES}): {bh6['verdict']}; raw P < 0.05: {bh6['raw_p_lt_0_05']}")
print("wrote", f"{WORK}/panel_d_raw.pdf")
