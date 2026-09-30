#!/usr/bin/env python3
"""ROUND 54 (2026-09-29, lane L4): Main Fig 4 panel d (the nine-cell cross, three apnea groups x two oxygen levels, five rows) at the
round-54 sizes, re-laid on work/layout_r54.json (a full-width row under panel c). The round-49 data wiring is kept: numbers/
alenfig3_cross_v1.json cross6 gated by the v8.1 sidecar, the owner's rows plus the Alopecia control row, Benjamini-Hochberg stars within
each cell (clinical conditions in one pool, the five controls in another), the CIELAB hazard-ratio ramp of round 27, text colour by
contrast, the Fig 6 claim check written as before. Sizes: headers, row labels, cell values 14 pt, the ramp title, its ticks and the
reference key 14 pt (keys), the asterisk key 13 pt (a small note). Cells widened from 58.07 to 80.79 pt (the widest cell string
(2.41-3.51) is 78.6 pt at 14 pt), 40 pt tall, headers 34 pt for two 14-pt lines. Outputs: work/panel_d_raw.pdf, work/panel_d_geometry.json,
verify/Main_Fig4_d_drawn.json, verify/FIG6_BH_CHECK_v8_1.json."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import json, math, os, re, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import l4lib as L, layout_r54 as LY
from statsmodels.stats.multitest import multipletests
matplotlib, plt = L.mpl_setup()
from matplotlib.patches import Rectangle

SHEET = "Main_Fig4"; SD = f"{L.LANE}/{SHEET}"; WORK = f"{SD}/work"; VER = f"{SD}/verify"
S = LY.load(); D_ = S["d"]; FS = S["sizes"]
SRC = f"{L.NUM}/alenfig3_cross_v1.json"; SC = L.sidecar(SRC); assert SC["step"]["rc"] == 0, SC["step"]
NCF = json.load(open(L.hydrated(f"{L.NUM}/negcontrols_final.json"))); L.sidecar(f"{L.NUM}/negcontrols_final.json")
XT = json.load(open(L.hydrated(f"{L.NUM}/crosstab_v2.json"))); L.sidecar(f"{L.NUM}/crosstab_v2.json")
X = json.load(open(L.hydrated(SRC))); C = X["cross6"]; CS = X["cell_sizes"]["cross6"]
GROUPS = ["No or mild", "Moderate", "Severe"]; LEVELS = ["<=1%", ">10%"]; REF = "No or mild|<=1%"
assert X["ahi3_groups"]["rule"] == {"No or mild": "AHI < 15", "Moderate": "15 <= AHI < 30", "Severe": "AHI >= 30"}, X["ahi3_groups"]
CONDS = [c for c, v in C.items() if not v["negative_control"]]; CTRLS = [c for c, v in C.items() if v["negative_control"]]
assert set(CTRLS) == set(NCF["panel_labels"]) and len(CTRLS) == 5, (CTRLS, NCF["panel_labels"])
SWAP = NCF.get("decision_record", ""); CONTROL_ROW = "Alopecia"; assert CONTROL_ROW in CTRLS, (CONTROL_ROW, CTRLS)
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
SEV_PRES = "Severe|<=1%"
bh6 = dict(cell=SEV_PRES, reference=REF, rule="Benjamini-Hochberg across the clinical conditions of cross6 within the cell (the sheet's rule)", source=SRC, source_sha256=L.sha256(SRC),
           conditions={c: dict(hr=C[c]["cells"][SEV_PRES]["hr"], lo=C[c]["cells"][SEV_PRES]["lo"], hi=C[c]["cells"][SEV_PRES]["hi"], p=C[c]["cells"][SEV_PRES]["p"], q=Q[c][SEV_PRES], n=C[c]["cells"][SEV_PRES]["n"], events=C[c]["cells"][SEV_PRES]["events"]) for c in CONDS},
           controls={c: dict(hr=C[c]["cells"][SEV_PRES]["hr"], p=C[c]["cells"][SEV_PRES]["p"], q=Q[c][SEV_PRES]) for c in CTRLS})
bh6["survivors_q_lt_0_05"] = sorted([c for c in CONDS if Q[c][SEV_PRES] < 0.05], key=lambda c: Q[c][SEV_PRES]); bh6["raw_p_lt_0_05"] = sorted(c for c in CONDS if C[c]["cells"][SEV_PRES]["p"] < 0.05)
bh6["verdict"] = ("sleep apnea alone (severe apnea with preserved oxygen): " + (f"{len(bh6['survivors_q_lt_0_05'])} of {len(CONDS)} conditions clear BH q < 0.05: {bh6['survivors_q_lt_0_05']}" if bh6["survivors_q_lt_0_05"] else f"no condition of {len(CONDS)} clears BH q < 0.05"))
json.dump(bh6, open(f"{VER}/FIG6_BH_CHECK_v8_1.json", "w"), indent=1)

# ---------------------------------------------------------------- geometry from the layout spec (page points, y down)
BX0, BY0, BX1, BY1 = D_["box"]; GX0 = D_["grid_x0"]; GRID_X1 = D_["grid_x1"]; CELL_W, CELL_PITCH, GRP_W, GRP_PITCH, CELL_H, ROW_PITCH = D_["cell_w"], D_["cell_pitch"], D_["group_w"], D_["group_pitch"], D_["cell_h"], D_["row_pitch"]
H1Y0, H1Y1 = D_["header1_y"]; H2Y0, H2Y1 = D_["header2_y"]; ROW0 = D_["row0_y"]; LAB_X1 = D_["lab_x1"]; GRID_BOTTOM = D_["grid_bottom"]
assert abs(ROW0 + (len(ROWS) - 1) * ROW_PITCH + CELL_H - GRID_BOTTOM) < 1e-6 and abs(GX0 + 2 * GRP_PITCH + GRP_W - GRID_X1) < 1e-6
RAMP = D_["ramp"]; REFBOX = D_["refbox"]; RAMP_STOPS = [("1", 0.0), ("1.5", math.log(1.5) / math.log(3.25)), ("2", math.log(2) / math.log(3.25)), ("3", math.log(3) / math.log(3.25))]
FS_C, FS_KEY, FS_NOTE = FS["cell"], FS["key"], FS["note"]; XH = 0.36 * FS_C; L1, L2 = -3.5, 11.5; HR_DY, CI_DY, REF_DY = 16.0, 31.0, 25.0
INK = L.INK; ORANGE_HDR = ["#f6d6c2", "#e4864f", "#d55e00"]; BLUE_HDR = ["#b3dcf2", "#0288d1"]
GROUP_LAB = {"No or mild": ["No or mild", "(AHI <15)"], "Moderate": ["Moderate", "(AHI 15 to <30)"], "Severe": ["Severe", "(AHI 30 or more)"]}
LEVEL_LAB = {"<=1%": ["Preserved", "T90 ≤1%"], ">10%": ["Low", "T90 >10%"]}
ROW_LAB = {"Cardiovascular composite": ["Cardiovascular", "composite"], "Type 2 diabetes": ["Type 2 diabetes"], "Acute kidney injury": ["Acute kidney", "injury"], "Heart failure": ["Heart failure"], CONTROL_ROW: [CONTROL_ROW, "(negative control)"]}
REF_KEY = "Reference, no or mild apnea with T90 ≤1%"; STAR_KEY = "* q < 0.05, ** q < 0.01, *** q < 0.001"; RAMP_TITLE = "Hazard ratio"
for s in (REF_KEY, STAR_KEY, RAMP_TITLE, *sum(GROUP_LAB.values(), []), *sum(LEVEL_LAB.values(), []), *sum(ROW_LAB.values(), [])): assert L.house_ok(s), s
def is_half(v): return bool(re.fullmatch(r"-?\d+\.\d\d5", repr(float(v))))

W, H = BX1 - BX0, BY1 - BY0
fig = plt.figure(figsize=(W / 72, H / 72)); fig.patch.set_alpha(0.0); ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(BX0, BX1); ax.set_ylim(BY1, BY0); ax.axis("off"); ax.patch.set_visible(False)
DR = []; STR = []; RECTS = []; HALVES = []
def text(x, y, s, size, ha, color=INK, z=5, **src):
    ax.text(x, y, s, fontsize=size, ha=ha, va="baseline", color=color, zorder=z); STR.append(dict(text=s, x=round(x, 3), baseline=round(y, 3), size=size, ha=ha, color=color, panel="d", **src))
def rect(x0, y0, w, h, fill, edge="none", lw=0, z=None, kind=""):
    kw = dict(facecolor=fill, edgecolor=edge, lw=lw); kw.update({} if z is None else dict(zorder=z)); ax.add_patch(Rectangle((x0, y0), w, h, **kw))
    RECTS.append(dict(kind=kind, rect=[round(x0, 3), round(y0, 3), round(x0 + w, 3), round(y0 + h, 3)], fill=fill, edge=edge, lw=lw))
ST = lambda t: dict(source_file="static:V13", source_key="V13 wording", source_value=t, rule="text")
RC = dict(L.RC); RC.update({"figure.facecolor": "none", "savefig.facecolor": "none"})
with plt.rc_context(RC):
    h1c = (H1Y0 + H1Y1) / 2; h2c = (H2Y0 + H2Y1) / 2
    text(LAB_X1, h1c + XH, "Sleep apnea", FS_C, "right", **ST("Sleep apnea")); text(LAB_X1, h2c + XH, "Oxygenation", FS_C, "right", **ST("Oxygenation"))
    for gi, g in enumerate(GROUPS):
        gx = GX0 + gi * GRP_PITCH; col = "#ffffff" if gi == 2 else INK
        rect(gx, H1Y0, GRP_W, H1Y1 - H1Y0, ORANGE_HDR[gi], kind=f"group header {g}")
        l1, l2 = GROUP_LAB[g]; text(gx + GRP_W / 2, h1c + L1, l1, FS_C, "center", color=col, **ST(l1)); text(gx + GRP_W / 2, h1c + L2, l2, FS_C, "center", color=col, **ST(l2))
        for li, lv in enumerate(LEVELS):
            cx0 = gx + li * CELL_PITCH; rect(cx0, H2Y0, CELL_W, H2Y1 - H2Y0, BLUE_HDR[li], kind=f"level header {g}|{lv}")
            a, b = LEVEL_LAB[lv]; text(cx0 + CELL_W / 2, h2c + L1, a, FS_C, "center", **ST(a)); text(cx0 + CELL_W / 2, h2c + L2, b, FS_C, "center", **ST(b))
    for ri, r in enumerate(ROWS):
        y0 = ROW0 + ri * ROW_PITCH; y1 = y0 + CELL_H; cy = (y0 + y1) / 2; labs = ROW_LAB[r]
        if len(labs) == 1: text(LAB_X1, cy + XH, labs[0], FS_C, "right", source_file=SRC, source_key="", source_value=labs[0], rule="text", note=f"row label = key of cross6 ({r})")
        else:
            text(LAB_X1, cy + L1, labs[0], FS_C, "right", source_file=SRC, source_key="", source_value=labs[0], rule="text", note=f"line 1 of the row label, key of cross6 ({r})")
            text(LAB_X1, cy + L2, labs[1], FS_C, "right", **(ST(labs[1]) if labs[1] == "(negative control)" else dict(source_file=SRC, source_key="", source_value=labs[1], rule="text", note=f"line 2 of the row label, key of cross6 ({r})")))
        for gi, g in enumerate(GROUPS):
            for li, lv in enumerate(LEVELS):
                cell = f"{g}|{lv}"; v = C[r]["cells"][cell]; cx0 = GX0 + gi * GRP_PITCH + li * CELL_PITCH; cxc = cx0 + CELL_W / 2
                if cell == REF:
                    text(cxc, y0 + REF_DY, "1.00", FS_C, "center", z=6, source_file=SRC, source_key=f"cross6/{r}/cells/{cell}/hr", source_value=1.0, rule="2dp_halfup")
                    DR.append(dict(row=r, cell=cell, hr=1.0, lo=None, hi=None, p=None, q=None, n=v["n"], events=v["events"], hr_text="1.00", ci_text=None, stars="", fill="#ffffff", text_color=INK,
                                   rect=[round(cx0, 3), round(y0, 3), round(cx0 + CELL_W, 3), round(y1, 3)], key=f"cross6[{r}].cells[{cell}]", source=SRC)); continue
                q = Q[r][cell]; st = L.stars(q); fill = L.ramp_hex(v["hr"]); col = L.text_on(fill); hr_s = L.r2(v["hr"]) + st; ci_s = f"({L.r2(v['lo'])}-{L.r2(v['hi'])})"
                for nm in ("hr", "lo", "hi"):
                    if is_half(v[nm]): HALVES.append(dict(row=r, cell=cell, field=nm, stored=v[nm], printed=L.r2(v[nm])))
                rect(cx0, y0, CELL_W, CELL_H, fill, kind=f"cell {r}|{cell}")
                text(cxc, y0 + HR_DY, hr_s, FS_C, "center", color=col, source_file=SRC, source_key=f"cross6/{r}/cells/{cell}", source_value=[v["hr"], st], rule="cross6_hr_stars", stars=st, q=q)
                text(cxc, y0 + CI_DY, ci_s, FS_C, "center", color=col, source_file=SRC, source_key=f"cross6/{r}/cells/{cell}", source_value=[v["lo"], v["hi"]], rule="cell_ci_paren")
                DR.append(dict(row=r, cell=cell, hr=v["hr"], lo=v["lo"], hi=v["hi"], p=v["p"], q=q, n=v["n"], events=v["events"], hr_text=hr_s, ci_text=ci_s, stars=st, fill=fill, text_color=col,
                               rect=[round(cx0, 3), round(y0, 3), round(cx0 + CELL_W, 3), round(y1, 3)], key=f"cross6[{r}].cells[{cell}]",
                               q_rule=(f"BH across the {len(CONDS)} clinical conditions in the cell" if r != CONTROL_ROW else "BH across the five negative controls in the cell"), source=SRC))
    rect(GX0, ROW0, CELL_W, GRID_BOTTOM - ROW0, "white", edge=INK, lw=0.9, z=3, kind="reference column box")
    n = 120; kx0, ky0, kx1, ky1 = RAMP
    for k in range(n):
        h0 = math.exp((k / n) * math.log(3.25)); h1 = math.exp(((k + 1) / n) * math.log(3.25))
        ax.add_patch(Rectangle((kx0 + k * (kx1 - kx0) / n, ky0), (kx1 - kx0) / n + 0.4, ky1 - ky0, facecolor=L.ramp_hex(math.sqrt(h0 * h1)), edgecolor="none"))
    rect(kx0, ky0, kx1 - kx0, ky1 - ky0, "none", edge="#ccd1d6", lw=0.4, kind="ramp frame")
    text((kx0 + kx1) / 2, D_["ramp_title_base"], RAMP_TITLE, FS_KEY, "center", **ST(RAMP_TITLE))
    for lab, frac in RAMP_STOPS: text(kx0 + frac * (kx1 - kx0), D_["ramp_tick_base"], lab, FS_KEY, "center", **ST(lab))
    rect(REFBOX[0], REFBOX[1], REFBOX[2], REFBOX[3], "white", edge=INK, lw=0.9, kind="reference key box")
    text(D_["ref_key_x"], D_["ref_key_base"], REF_KEY, FS_KEY, "left", **ST(REF_KEY)); text(D_["ref_key_x"], D_["star_key_base"], STAR_KEY, FS_NOTE, "left", **ST(STAR_KEY))
    BOTTOM_INK = max(s["baseline"] for s in STR); assert BOTTOM_INK + 4 <= BY1, (BOTTOM_INK, BY1)
    rend = fig.canvas.get_renderer()
    for t in ax.texts:
        bb = t.get_window_extent(renderer=rend); x0p, x1p = BX0 + bb.x0 / fig.dpi * 72, BX0 + bb.x1 / fig.dpi * 72; y0p, y1p = BY1 - bb.y1 / fig.dpi * 72, BY1 - bb.y0 / fig.dpi * 72
        assert x0p >= BX0 + 1.0 and x1p <= BX1 - 1.0 and y0p >= BY0 + 1.0 and y1p <= BY1 - 1.0, ("text leaves panel d's page", t.get_text(), x0p, x1p, y0p, y1p)
        assert t.get_fontsize() >= 13.0 and L.house_ok(t.get_text())
    fig.savefig(f"{WORK}/panel_d_raw.pdf", transparent=True); plt.close(fig)
# ---------------------------------------------------------------- proofs: values against the file, the V13 string set
for d in DR:
    v = C[d["row"]]["cells"][d["cell"]]
    if d["cell"] == REF: assert d["hr_text"] == "1.00" and v["hr"] == 1.0; continue
    assert d["hr_text"].rstrip("*") == L.r2(v["hr"]) and d["ci_text"] == f"({L.r2(v['lo'])}-{L.r2(v['hi'])})" and d["stars"] == L.stars(d["q"]) and d["fill"] == L.ramp_hex(v["hr"]), d
G0 = json.load(open(f"{WORK}/v13_geometry.json"))["panel_d"]
R49 = json.load(open(f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/L4/Main_Fig4/verify/Main_Fig4_d_drawn.json"))["drawn_records"]   # V30's panel d strings
assert sorted(s["text"] for s in STR) == sorted(s["text"] for s in R49), ("the string set differs from V30's", sorted(set(s["text"] for s in STR) ^ set(s["text"] for s in R49)))
old_by = {(o["row"], o["cell"]): o for o in G0["v13_values"]}; CHANGES = []
for d in DR:
    key_old = ("Back pain", d["cell"]) if d["row"] == CONTROL_ROW else (d["row"], d["cell"]); o = old_by.get(key_old)
    if o is None or d["hr_text"] != o["hr_text"] or d["ci_text"] != o["ci_text"] or d["fill"] != o["fill"]:
        CHANGES.append(dict(row=d["row"], row_v13=key_old[0], cell=d["cell"], before=(f"{o['hr_text']} {o['ci_text']}" if o else None), after=f"{d['hr_text']} {d['ci_text']}", q_before=(o["q"] if o else None), q_after=d["q"]))
geom = dict(panel_box=[BX0, BY0, BX1, BY1], page_box=[BX0, BY0, BX1, BY1], grid_x0=GX0, grid_x1=GRID_X1, group_w=GRP_W, group_pitch=GRP_PITCH, cell_w=CELL_W, cell_pitch=CELL_PITCH, cell_h=CELL_H, row_pitch=ROW_PITCH,
            header1_y=[H1Y0, H1Y1], header2_y=[H2Y0, H2Y1], row0_y=ROW0, grid_bottom=GRID_BOTTOM, label_x1=LAB_X1, key=dict(ramp=RAMP, refbox=REFBOX, ramp_title_base=D_["ramp_title_base"], ramp_tick_base=D_["ramp_tick_base"], ref_key_base=D_["ref_key_base"], star_key_base=D_["star_key_base"]),
            rects=RECTS, bottom_ink_baseline=BOTTOM_INK, design="round-54 build for the record only (not delivered, owner option A), cells 80.79 x 40 pt")
json.dump(geom, open(f"{WORK}/panel_d_geometry.json", "w"), indent=1)
json.dump(dict(sheet=SHEET, panel="d", lane="L4 round 54", source=SRC, source_sha256=L.sha256(SRC), sidecar_step=SC["step"], v14_gate=SC["_v14_gate"],
               cross="cross6", reference=REF, cell_sizes=CS, set_aside_t90_gt1_le10=SET_ASIDE, cohort_n=X["provenance"]["cohort_n"], rows=ROWS, control_row=CONTROL_ROW, control_row_v13="Back pain",
               control_swap=SWAP, clinical_pool=CONDS, control_pool=CTRLS, values=DR, strings=STR, drawn_records=STR, changes_vs_v13=CHANGES, exact_half_strings=HALVES,
               stored_precision="3 dp in the numbers file, printed half up on the stored value", fig6_bh_check=bh6["verdict"], geometry=geom), open(f"{VER}/{SHEET}_d_drawn.json", "w"), indent=1, ensure_ascii=False)
print(f"panel d: {len(ROWS)} rows x 6 cells (cells {CELL_W:.2f} x {CELL_H:.1f}, grid x {GX0:.1f}..{GRID_X1:.1f}, y {H1Y0:.1f}..{GRID_BOTTOM:.1f}, lowest baseline {BOTTOM_INK:.1f}); control row {CONTROL_ROW}; {len(CHANGES)} of 30 cells differ from V13; halves {len(HALVES)}; {len(STR)} strings")
print(f"  FIG 6 CHECK ({SEV_PRES}): {bh6['verdict']}; wrote {WORK}/panel_d_raw.pdf")
