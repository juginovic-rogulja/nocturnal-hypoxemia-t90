#!/usr/bin/env python3
"""Round 49 LED-B copy (2026-09-26): every text one point larger (LABF 12, TCKF 11, ANNF 10.5), the star key not drawn (LNOTES round 40), the row-label rule against the plot's leftmost ink. Otherwise the round-38 copy. ED_Fig07 panel a (V14 lane L4, round 37): the within-severe-apnea forest, T90 above 10% against 10% or less, re-plotted from
numbers/severe_apnea_negctrl_v1.json (step 126, v8.1 sidecar gated, cross-checked against extras2_v2.json) in the V13 design:
the round-30 re-implementation of the eFigure 9 block (Arial Type 42, ink #1a1d21) with every geometry constant READ FROM THE
ROUND-30 MEASUREMENT OF THE BASE FOREST (work record panel_a_geometry.json: axes box, log map from the tick marks, row pitch, marker
sizes and colours, line widths, dashed reference rule, estimate column x, label gutter), which the V13 sheet still carries byte for
byte 8 pt down the page (proved in 04_verify against the V13 text layer). Rows sorted by hazard ratio ascending, the negative
controls (the v8.1 five: alopecia, inguinal hernia, glaucoma, contact dermatitis, hemorrhoids) as the top block under the bold
"Negative controls" header, a 1.5-row gap, then the conditions. Stars from the raw two-sided P (the file carries no q), estimates
half up to 2 dp on the file's 3 dp values (exact halves listed). Copy of the round-30 01_build_a.py (beside as ed7_01_build_a_PRE_V8_1.py)
without the base-sheet re-measurement and the round-30 L12 full-precision block.
Outputs: work/panel_a_raw.pdf, work/panel_a.pdf (font metrics aligned to the V13 sheet), work/panel_a_geometry.json, verify/ED_Fig07_a_drawn.json."""
import gc, json, math, os, re, sys
import fitz, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ed7spec as S
L = S.L
matplotlib, plt = L.mpl_setup()
from matplotlib import ticker as mticker
from matplotlib.transforms import blended_transform_factory
os.makedirs(S.WORK, exist_ok=True); os.makedirs(S.CROPS, exist_ok=True)

# ---------------------------------------------------------------- numbers, gated (v8.1)
SC = L.sidecar(S.SRC_A); assert SC["step"]["rc"] == 0 and str(SC["step"]["id"]) == "126", SC["step"]
SEV = json.load(open(L.hydrated(S.SRC_A)))
NCF = json.load(open(L.hydrated(f"{L.NUM}/negcontrols_final.json"))); L.sidecar(f"{L.NUM}/negcontrols_final.json")
SCX = L.sidecar(S.SRC_A_XCHECK); XC = json.load(open(L.hydrated(S.SRC_A_XCHECK)))
xo = XC["severe_apnea"]["outcomes"]
assert set(k.replace(" [NEG]", "") for k in xo) == set(SEV["outcomes"]), (sorted(xo), sorted(SEV["outcomes"]))
for k, v in xo.items():
    lab = k.replace(" [NEG]", ""); w = SEV["outcomes"][lab]
    assert w["negative_control"] == k.endswith("[NEG]"), k
    for f in ("hr", "lo", "hi", "p", "events"): assert v[f] == w[f], ("extras2_v2 disagrees", k, f, v[f], w[f])
assert XC["severe_apnea"]["n"] == SEV["n"], (XC["severe_apnea"]["n"], SEV["n"])
assert SEV["n_reference"] + SEV["n_low_oxygen"] == SEV["n"], (SEV["n_reference"], SEV["n_low_oxygen"], SEV["n"])
assert L.pct1(SEV["n_reference"], SEV["n"]) == str(SEV["pct_reference"]) and L.pct1(SEV["n_low_oxygen"], SEV["n"]) == str(SEV["pct_low_oxygen"])
assert sorted(SEV["negative_controls_revised"]) == sorted(k for k, v in SEV["outcomes"].items() if v["negative_control"]) == sorted(NCF["panel_labels"]), (SEV["negative_controls_revised"], NCF["panel_labels"])

items = sorted(SEV["outcomes"].items(), key=lambda kv: kv[1]["hr"])
labs = [k for k, _ in items]; hr = [v["hr"] for _, v in items]; lo = [v["lo"] for _, v in items]; hi = [v["hi"] for _, v in items]
pv = [v["p"] for _, v in items]; neg = [v["negative_control"] for _, v in items]
assert all(isinstance(p, float) and math.isfinite(p) for p in pv), pv
assert all(l < h <= u or l <= h < u for l, h, u in zip(lo, hr, hi)), list(zip(labs, lo, hr, hi))
assert all("q" not in v for _, v in items), "a q appeared in the source, the key must change"
STARS = [L.stars(p) for p in pv]
neg_idx = [i for i, n in enumerate(neg) if n]; N_CTRL = len(neg_idx); N_COND = len(labs) - N_CTRL; assert N_CTRL >= 1 and N_COND >= 1
CONTIGUOUS_TOP = neg_idx == list(range(N_CTRL))
assert any(STARS), "no row is marked, the significance key would be dead ink"
for lvl in sorted(set(STARS) - {""}): assert f"{lvl} P <" in S.KEY_STARS_A, lvl
for s in (S.KEY_STARS_A, S.BLOCK_HEAD, S.XLAB_A, *labs): assert L.house_ok(s), s
def is_half(v): return bool(re.fullmatch(r"-?\d+\.\d\d5", repr(float(v))))
HALVES = [dict(label=labs[i], field=f, stored=v) for i in range(len(labs)) for f, v in (("hr", hr[i]), ("lo", lo[i]), ("hi", hi[i])) if is_half(v)]

# ---------------------------------------------------------------- the geometry: the round-30 measurement of the base forest (the V13 panel a)
G = json.load(open(L.hydrated(S.R30_GEOM_A))); BW, BH = G["page_w"], G["page_h"]
AX_L, AX_T, AX_R, SPINE_Y = G["axes_panel"]; XLO, XHI = G["xlim"]; A, B = G["log_map"]["A"], G["log_map"]["B"]; TICKS = G["ticks"]; tick_x = G["tick_x"]
UNIT = G["unit_pt"]; BLOCK_GAP = G["block_gap_units"]; EST_X = G["est_x"]; X_NUM = G["x_num"]; SPINE_LW = G["spine_lw"]; CI_LW = G["ci_lw"]; REF_LW = G["ref_lw"]; REF_DASH = tuple(G["ref_dash"])
S_CIRCLE, S_SQUARE = G["marker_area"]["circle"], G["marker_area"]["square"]; MEDGE = G["marker_edge"]; COL_COND, COL_CTRL, INK = G["colours"]["condition"], G["colours"]["control"], G["colours"]["ink"]
XLAB_DROP_IN, KEY_Y_IN = G["xlab_drop_in"], G["key_y_in"]; TICK_LEN = 3.0
BT = json.load(open(L.hydrated(S.BASE_TEXT))); MET = L.sheet_font_metrics(BT["spans"])
# the V13 sheet must still carry that forest: tick labels, estimate column and the star key at the round-30 positions (shifted by SHIFT_A)
SPV = L.dedupe(BT["spans"])
for t, x in zip(TICKS, tick_x):
    h = [s for s in SPV if s["text"] == f"{t:g}" and abs((s["bbox"][0] + s["bbox"][2]) / 2 - x) < 0.6 and abs(s["origin"][1] - (SPINE_Y + TICK_LEN + 3.5 + S.SHIFT_A + 7.156)) < 2.5]
    assert h, ("V13 tick label not at the round-30 position", t, x)
est13 = [s for s in SPV if abs(s["size"] - 9.5) < 0.01 and re.fullmatch(r"\d\.\d\d \(\d\.\d\d-\d\.\d\d\)", s["text"].replace("\xa0", " "))]
assert len(est13) == len(labs) and all(abs(s["origin"][0] - EST_X) < 0.05 for s in est13), (len(est13), EST_X)
W_IN, H_IN = BW / 72, BH / 72; LEFT = AX_L / 72; PLOTW = (AX_R - AX_L) / 72; AXH = (SPINE_Y - AX_T) / 72; BOT = (BH - SPINE_Y) / 72; TOP = AX_T / 72
GUT = S.LETTER_X / 72; UNIT_IN = UNIT / 72
Y_CTRL_TOP = (N_COND - 1) + BLOCK_GAP + (N_CTRL - 1); Y_HEAD = Y_CTRL_TOP + 1.0; YLO, YHI = -0.8, Y_HEAD + 0.6
assert abs((YHI - YLO) * UNIT - (SPINE_Y - AX_T)) < 0.1, ((YHI - YLO) * UNIT, SPINE_Y - AX_T)   # same row count as V13: the y design reproduces the axes height
ci_ = [i for i in range(len(labs)) if neg[i]]; cc_ = [i for i in range(len(labs)) if not neg[i]]
def y_of(i): return (Y_CTRL_TOP - ci_.index(i)) if neg[i] else (len(cc_) - 1 - cc_.index(i))
yy = np.array([y_of(i) for i in range(len(labs))], dtype=float)
for v in lo + hi + hr: assert XLO < v < XHI, ("value off the axis", v, XLO, XHI)
LABF, TCKF, ANNF = 11.0 + L.PT_PLUS, 10.0 + L.PT_PLUS, 9.5 + L.PT_PLUS   # round 49: 12, 11, 10.5
rc = dict(L.RC); rc.update({"font.size": TCKF, "axes.labelsize": LABF, "axes.titlesize": LABF, "xtick.labelsize": TCKF, "ytick.labelsize": TCKF, "legend.fontsize": TCKF,
                            "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False, "figure.facecolor": "none", "savefig.facecolor": "none"})
STR = []
def rec(text, x, y, size, ha, va, weight="normal", color=INK, **src): STR.append(dict(text=text, x=round(x, 3), y=round(y, 3), size=size, ha=ha, va=va, weight=weight, color=color, coords="panel", **src))
with plt.rc_context(rc):
    fig = plt.figure(figsize=(W_IN, H_IN)); fig.patch.set_alpha(0.0); fig.canvas.draw(); rend = fig.canvas.get_renderer()
    # round 49: at 11 pt the widest row label reaches past the axes box's invisible left edge (no left spine, no grid); the rule becomes
    # "clear of the plot's leftmost ink by 4 pt" (the reference rule at x = A and every interval start), the axes-edge crossing is declared
    INK_LEFT_IN = min([A + B * math.log(v) for v in lo] + [A]) / 72.0; LABEL_ENDS = {}
    for s in labs:
        t = fig.text(0.5, 0.5, s, fontsize=TCKF); wid = t.get_window_extent(renderer=rend).width / fig.dpi; t.remove(); LABEL_ENDS[s] = round((GUT + wid) * 72.0, 2)
        assert GUT + wid + 4.0 / 72.0 <= INK_LEFT_IN, ("row label runs into the plot's ink", s, wid * 72.0, INK_LEFT_IN * 72.0)
    LABELS_PAST_AXES_EDGE = {s: round(e - AX_L, 2) for s, e in LABEL_ENDS.items() if e > AX_L}
    ax = fig.add_axes([LEFT / W_IN, BOT / H_IN, PLOTW / W_IN, AXH / H_IN])
    ax.axvline(1.0, color=INK, lw=REF_LW, ls=(0, REF_DASH), zorder=1)
    cols = [COL_CTRL if n else COL_COND for n in neg]
    for i in range(len(labs)):
        ax.plot([lo[i], hi[i]], [yy[i], yy[i]], color=cols[i], lw=CI_LW, solid_capstyle="round", zorder=2)
        ax.scatter([hr[i]], [yy[i]], s=S_SQUARE if neg[i] else S_CIRCLE, marker="s" if neg[i] else "o", zorder=3, facecolors=cols[i], edgecolors="white", linewidths=MEDGE)
        if STARS[i]:
            ax.text(hi[i] * 1.04, yy[i], STARS[i], va="center", fontsize=ANNF, color=cols[i])
            rec(STARS[i], A + B * math.log(hi[i] * 1.04), AX_T + (YHI - yy[i]) * UNIT, ANNF, "left", "center", color=cols[i], source_file=S.SRC_A, source_key=f"outcomes/{labs[i]}", source_value=pv[i], rule="sev_stars_p")
    ax.set_yticks(list(yy)); ax.set_yticklabels(labs, fontsize=TCKF, color=INK, ha="left"); ax.tick_params(axis="y", pad=(LEFT - GUT) * 72.0, length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xscale("log"); ax.set_xticks(TICKS); ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.xaxis.set_minor_formatter(mticker.NullFormatter()); ax.tick_params(axis="x", which="minor", length=0); ax.xaxis.set_minor_locator(mticker.NullLocator())
    ax.set_xlim(XLO, XHI); ax.set_ylim(YLO, YHI); ax.grid(False)
    tr = blended_transform_factory(ax.transAxes, ax.transData)
    ax.text(-(LEFT - GUT) / PLOTW, Y_HEAD, S.BLOCK_HEAD, transform=tr, fontsize=TCKF, fontweight="bold", va="center", ha="left", color=INK, clip_on=False, zorder=6)
    rec(S.BLOCK_HEAD, S.LETTER_X, AX_T + (YHI - Y_HEAD) * UNIT, TCKF, "left", "center", weight="bold", source_file="static:V13", source_key="V13 wording", source_value=S.BLOCK_HEAD, rule="text")
    EST = []
    for i in range(len(labs)):
        e = L.hr_ci(hr[i], lo[i], hi[i]); EST.append(e)
        ax.text(X_NUM, yy[i], e, transform=tr, fontsize=ANNF, va="center", ha="left", color=INK)
        rec(e, EST_X, AX_T + (YHI - yy[i]) * UNIT, ANNF, "left", "center", source_file=S.SRC_A, source_key=f"outcomes/{labs[i]}", source_value=[hr[i], lo[i], hi[i]], rule="sev_hr_ci")
        rec(labs[i], S.LETTER_X, AX_T + (YHI - yy[i]) * UNIT, TCKF, "left", "center_baseline", source_file=S.SRC_A, source_key="", source_value=labs[i], rule="text", note="row label = key of outcomes")
    ax.set_xlabel(S.XLAB_A, fontsize=LABF, labelpad=4); ax.xaxis.set_label_coords(0.5, -XLAB_DROP_IN / AXH)
    rec(S.XLAB_A, (AX_L + AX_R) / 2, SPINE_Y + XLAB_DROP_IN * 72, LABF, "center", "top", source_file="static:V13", source_key="V13 wording", source_value=S.XLAB_A, rule="text")
    for t, x in zip(TICKS, tick_x): rec(f"{t:g}", x, SPINE_Y + TICK_LEN + 3.5, TCKF, "center", "top", source_file="static:V13", source_key="V13 axis tick", source_value=f"{t:g}", rule="text")
    # round 49: the star key is NOT drawn (removed from the delivered sheet by lane LNOTES, round 40; the legend carries the sentence)
    for tobj in list(fig.texts) + [t for a_ in fig.axes for t in a_.texts] + ax.get_xticklabels() + ax.get_yticklabels() + [ax.xaxis.label]:
        if tobj.get_text().strip(): assert tobj.get_fontsize() >= 9.0 and L.house_ok(tobj.get_text()), tobj.get_text()
    fig.savefig(S.PANEL_A_RAW, transparent=True); plt.close(fig)
done = L.align_font_metrics(S.PANEL_A_RAW, S.PANEL_A, MET); assert set(done) >= {"ArialMT", "Arial-BoldMT"}, done
# ---------------------------------------------------------------- read back the panel's own text layer -> drawn records in PAGE coordinates
RB = L.spans_of(S.PANEL_A); taken = set(); recs = []
def take(text, y_hint, x_hint=None, **src):
    cands = [(i, s) for i, s in enumerate(RB) if i not in taken and s["text"].replace("\xa0", " ") == text]
    assert cands, ("drawn string not on the panel text layer", text)
    cands.sort(key=lambda t: abs(t[1]["origin"][1] - y_hint) + (0 if x_hint is None else 0.01 * abs(t[1]["origin"][0] - x_hint)))
    i, s = cands[0]; taken.add(i)
    recs.append(dict(text=text, x=round(s["origin"][0], 3), baseline=round(s["origin"][1] + S.SHIFT_A, 3), ha="left", size=s["size"], panel="a", **src)); return recs[-1]
for r in STR:
    y_hint = r["y"] + (0 if r["va"] in ("center", "center_baseline") else (7.2 if r["va"] == "top" else -2.0))
    take(r["text"], y_hint, r["x"], source_file=r["source_file"], source_key=r["source_key"], source_value=r["source_value"], rule=r["rule"], **({"note": r["note"]} if "note" in r else {}))
left = [s["text"] for i, s in enumerate(RB) if i not in taken and s["text"].strip()]; assert not left, ("panel strings not accounted for", left)
DR_ROWS = []
for i in range(len(labs)):
    yl = AX_T + (YHI - yy[i]) * UNIT
    DR_ROWS.append(dict(label=labs[i], negative_control=neg[i], hr=hr[i], lo=lo[i], hi=hi[i], p=pv[i], events=items[i][1]["events"], n=items[i][1]["n"], estimate_text=EST[i], star_text=STARS[i],
                        star_rule="raw two-sided P: * <0.05, ** <0.01, *** <0.001", row_unit=float(yy[i]), row_y_panel=round(yl, 3), row_y_page=round(yl + S.SHIFT_A, 3), marker_x=round(A + B * math.log(hr[i]), 3),
                        ci_x=[round(A + B * math.log(lo[i]), 3), round(A + B * math.log(hi[i]), 3)], colour=cols[i], marker="square" if neg[i] else "circle", key=f"outcomes/{labs[i]}", source=S.SRC_A))
for d in DR_ROWS:
    v = SEV["outcomes"][d["label"]]; assert d["estimate_text"] == L.hr_ci(v["hr"], v["lo"], v["hi"]) and d["star_text"] == L.stars(v["p"]), d
geom = dict(G, panel_box_page=[0, S.SHIFT_A, BW, BH + S.SHIFT_A], shift=S.SHIFT_A, pt_plus=L.PT_PLUS, sizes=dict(label=LABF, tick=TCKF, annot=ANNF), label_ends_pt=LABEL_ENDS, labels_past_axes_edge_pt=LABELS_PAST_AXES_EDGE, ink_left_pt=round(INK_LEFT_IN * 72.0, 3), star_key_drawn=False, y_range=[YLO, YHI], y_head=Y_HEAD, y_ctrl_top=Y_CTRL_TOP, n_rows=len(labs), font_metrics_aligned=done, source_geometry=S.R30_GEOM_A)
json.dump(geom, open(S.GEOM_A, "w"), indent=1)
json.dump(dict(sheet=S.SHEET, panel="a", lane="V14_L4_APNEA", source=S.SRC_A, source_sha256=L.sha256(S.SRC_A), sidecar_step=SC["step"], v14_gate=SC["_v14_gate"],
               cross_check=dict(file=S.SRC_A_XCHECK, sha256=L.sha256(S.SRC_A_XCHECK), agreed_rows=len(xo)), exposure=SEV["exposure"], subgroup=SEV["subgroup"],
               cohort=dict(n=SEV["n"], n_reference=SEV["n_reference"], pct_reference=SEV["pct_reference"], n_low_oxygen=SEV["n_low_oxygen"], pct_low_oxygen=SEV["pct_low_oxygen"], n_normal_oxygen_le1=SEV["n_normal_oxygen"], pct_normal_oxygen_le1=SEV["pct_normal_oxygen"]),
               negative_controls=SEV["negative_controls_revised"], controls_contiguous_at_top=CONTIGUOUS_TOP, n_controls=N_CTRL, n_conditions=N_COND,
               row_order_top_to_bottom=[labs[i] for i in sorted(range(len(labs)), key=lambda i: -yy[i])], rows=DR_ROWS, strings=STR, drawn_records=recs,
               key_strings=dict(block_head=S.BLOCK_HEAD, star_key=S.KEY_STARS_A, x_label=S.XLAB_A), stored_precision="3 dp in the numbers file, printed half up on the stored value", exact_half_strings=HALVES,
               geometry=geom, text_layer_panel=RB), open(S.DRAWN_A, "w"), indent=1, ensure_ascii=False)
print(f"panel a: {len(labs)} rows ({N_CTRL} controls on top, contiguous {CONTIGUOUS_TOP}), n {SEV['n']:,} (reference {SEV['n_reference']:,}, low oxygen {SEV['n_low_oxygen']:,}), axes x {AX_L:.2f}..{AX_R:.2f}, xlim {XLO:.4f}..{XHI:.4f}, unit {UNIT:.3f} pt")
for d in DR_ROWS: print(f"  {d['label']:<26} {d['estimate_text']:<18} {d['star_text']:<4} p {d['p']:.2e} y_page {d['row_y_page']}")
print(f"  exact halves {len(HALVES)}: {[(h['label'], h['field'], h['stored']) for h in HALVES]}; {len(recs)} drawn records; wrote {S.PANEL_A}")
