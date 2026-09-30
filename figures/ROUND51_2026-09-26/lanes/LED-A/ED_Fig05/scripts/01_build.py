#!/usr/bin/env python3
"""ROUND 51 lane LED-A copy (2026-09-26, item 2): no values above the bars, pale and dark green, the prevalence ratio with its 95% CI as the second tick-label
line under each pair and the P value as asterisks above the pair. Original of this copy: ROUND49_2026-09-26/lanes/LED-A/ED_Fig05/scripts/01_build.py.
ROUND 49 lane LED-A copy (2026-09-26): +1 pt through the copied helpers (supp_polishB_common TITLE_PT/TICK_PT/ANN_PT, splitstyle PANEL) and PT_TITLE 11,
the axes widths kept at the V26 layout (KEY_W measured at 10 pt), output work/ED_Fig05_built.pdf (the strip is stamped afterwards).
Original: ROUND37_2026-09-14/figures/lanes/V14_L3a_DURATION_MAIN/ED_Fig05/scripts/01_build.py.
ED Fig 5 (four panels a to d), V14 rebuild on the v8.1 tables at the V13 design. Repointed copy of
ROUND30_2026-09-08/V7_L3a_DURATION_MAIN/ED_Fig05/scripts/01_build.py (backup 01_build_PRE_V8_1.py). Design = efig18_polish.py round 5
(FINAL_FIGURES_2026-08-14/_scripts) as it shows on the V13 sheet (byte-identical to V7). Data recomputed from
data_frozen_v8_2026-09/t90_final.parquet through numbers/cohort_spec.apply_cohort with T90_FULL_NIGHTS_ONLY=1 (the 15,551 full
diagnostic nights, the v8.1 rule of step 129) and the v2 habitual sleep extraction, exactly the arithmetic of
make_efig18_shortsleep_oxygen.py (step 129, rc 0 on 2026-09-14 14:03), and cross-checked against the eight cells that step printed
to its v8.1 log. Colours sampled from the V13 sheet (#79d475 sleep under 5 h, #39c445 sleep 7 h or more).
V14 changes against the R30 copy: paths (LANE, BASE = V13, PARQ = the v8 table, S129LOG = the v8.1 driver log), the cohort env flag,
the drawn record carries v14lib DRAWN_RECORD fields for every printed string."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import json, os, re, sys
os.environ["T90_FULL_NIGHTS_ONLY"] = "1"                      # the v8.1 cohort rule of step 129 (set before cohort_spec is imported)
import numpy as np, pandas as pd
from scipy.stats import norm
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND51_2026-09-26/lanes/LED-A/scripts_shared")   # round 49: the lane's patched helpers (+1 pt)
import lane_common as C
C.load_ff()
from supp_polishB_common import (FIGW, EDGE, INK, TITLE_PT, TICK_PT, ANN_PT, rc_polish, letter, text_w, axis_covers, plt)
import supp_polishB_common as SB
sys.path.insert(0, paths.ANALYSIS_DIR)
from cohort_spec import apply_cohort, COHORT_N, FULL_NIGHTS_ONLY
import fitz
assert FULL_NIGHTS_ONLY and COHORT_N == 15551, (FULL_NIGHTS_ONLY, COHORT_N)

SHEET = "ED_Fig05"; HERE = f"{C.LANE}/{SHEET}"; WORK = f"{HERE}/work"; VER = f"{HERE}/verify"
os.makedirs(WORK, exist_ok=True); os.makedirs(VER, exist_ok=True)
BASE = f"{C.BASE}/{SHEET}.pdf"
PARQ = f"{C.V8}/t90_final.parquet"
HAB = f"{C.SV}/habitual_sleep_v2/habitual_records_v2_full.parquet"
S129LOG = f"{paths.V8_ROOT}/logs/129_make_efig18_shortsleep_oxygen.log"
for p in (BASE, PARQ, HAB, S129LOG): C.hydrated(p)
PRIMARY, DEEMPH = "#79d475", "#298d32"          # round 51 (2b): pale green for sleep under 5 h (V28 #79d475 kept), dark green #298d32 for 7 h or more (V28 #39c445 was too close: luma gap 34, now 73)
SB.sep_ok("ED_Fig05 bar pair", [PRIMARY, DEEMPH])
import splitstyle
PR_LAYOUT = os.environ.get("ED5_PR_LAYOUT", "one_line")   # round 51 (2c): "one_line" = the brief (PR line as the second tick-label line, panels widened to fit), "two_lines" = the alternative at the V28 width
OUTD = os.environ.get("ED5_OUT_DIR", WORK); os.makedirs(OUTD, exist_ok=True)

CARD = ["copd2", "asthma", "resp_failure", "obesity_hypovent", "pulm_htn", "hf", "ihd", "mi", "afib", "stroke_any"]
REN = ["ckd", "cirrhosis", "cancer_any"]; MET = ["diabetes", "obesity", "nafld"]     # the step-129 script's own lists (organ and metabolic disease)
b = apply_cohort(pd.read_parquet(PARQ)); assert len(b) == 15551, len(b)
_h = pd.read_parquet(HAB)
_u = _h[(_h.veto == "") & _h.tier.isin(["A", "B", "C"]) & _h.hours_v2.between(0.5, 14)]
_u = _u[~_u.lab_echo.fillna(False).astype(bool)].rename(columns={"person_id": "BDSPPatientID"})
_per = (_u.groupby("BDSPPatientID").hours_v2.median() * 60).rename("hab").reset_index()
b = b.merge(_per, on="BDSPPatientID", how="left")
healthy = pd.Series(True, index=b.index)
for k in CARD + REN + MET:
    healthy &= (b[f"{k}_prevalent"] != 1)


def cell(frame, cut, col):
    f = frame[frame[col].notna()]
    s, r = f[f[col] < 300], f[f[col] >= 420]
    k1, n1 = int((s.spo2_pct_below_90 > cut).sum()), len(s)
    k0, n0 = int((r.spo2_pct_below_90 > cut).sum()), len(r)
    rr = (k1 / n1) / (k0 / n0)
    se = np.sqrt(1 / k1 - 1 / n1 + 1 / k0 - 1 / n0)
    p = float(2 * (1 - norm.cdf(abs(np.log(rr) / se))))
    return 100 * k1 / n1, 100 * k0 / n0, n1, n0, rr, rr * np.exp(-1.96 * se), rr * np.exp(1.96 * se), p, k1, k0


PANELS = [("Everyone, laboratory sleep", b, "TST_min", "A"), ("Disease-free, laboratory sleep", b[healthy], "TST_min", "B"),
          ("Everyone, home-reported sleep", b, "hab", "C"), ("Disease-free, home-reported sleep", b[healthy], "hab", "D")]
GOT = {(L, c): cell(frame, c, col) for title, frame, col, L in PANELS for c in (5, 10)}

# cross-check against the step-129 driver log (the same arithmetic run by the v8.1 chain on 2026-09-14 14:03)
log = open(S129LOG).read()
assert "T90_FULL_NIGHTS_ONLY=1" in log and "data_frozen_v8_2026-09" in log, "the step-129 log is not the v8.1 run"
rx = re.compile(r"^\s+([ABCD]) .{38}\s*T90>\s*(5|10)%\s+([\d.]+) vs\s+([\d.]+)\s+([\d.]+) \(([\d.]+)-([\d.]+)\)", re.M)
S129 = {(m.group(1), int(m.group(2))): tuple(float(m.group(i)) for i in range(3, 8)) for m in rx.finditer(log)}
assert len(S129) == 8, S129
for k, g in GOT.items():
    mine = (round(g[0], 1), round(g[1], 1), round(g[4], 2), round(g[5], 2), round(g[6], 2))
    assert mine == S129[k], (k, mine, S129[k])
print("step-129 cross-check: all 8 cells agree with the v8.1 driver log")


def fmt_p(p):
    assert np.isfinite(p) and 0.0 <= p <= 1.0, p
    return "P < 0.001" if p < 0.001 else "P = " + f"{p:.3f}"


PSTR = {k: fmt_p(g[7]) for k, g in GOT.items()}

W_IN = FIGW
HEAD_H = EDGE; TITLE_BAND = 0.86; AX_H = 1.95; BELOW_V28 = 1.18; ROW_GAP = 0.26; X0_LEFT = 0.60; MID = 0.66; KEY_GAP = 0.08; PT_TITLE = 11.0   # round 49: panel titles 10 -> 11 pt
BELOW = {"one_line": 0.55, "two_lines": 0.76}[PR_LAYOUT]   # round 51 (2c): the five-line block is gone, the space below the axes holds two (or three) tick-label lines
DELTA = 1.0; TICK_PT_V26 = TICK_PT - DELTA; assert TICK_PT == 11.0 and TITLE_PT == 12.0 and ANN_PT == 10.5, (TICK_PT, TITLE_PT, ANN_PT)
yrow = [HEAD_H, HEAD_H + TITLE_BAND + AX_H + BELOW + ROW_GAP]
H = yrow[1] + TITLE_BAND + AX_H + BELOW + 0.06 + EDGE
TICKLAB = ["T90 >5%", "T90 >10%"]
LEG_LABELS = ("Sleep under 5 h", "Sleep 7 h or more"); LEG_WRAPPED = ("Sleep under\n5 h", "Sleep 7 h\nor more")
DRAWN = {"figure": SHEET, "sources": {"parquet": PARQ, "parquet_sha256": C.sha256(PARQ), "habitual": HAB, "habitual_sha256": C.sha256(HAB),
                                      "cohort_filter": "numbers/cohort_spec.apply_cohort with T90_FULL_NIGHTS_ONLY=1 (15,551 full diagnostic nights)", "step129_log": S129LOG},
         "base_sheet": BASE, "base_sha256": C.sha256(BASE), "panels": {}, "colours": {"under_5h": PRIMARY, "7h_or_more": DEEMPH}, "records": []}
with plt.rc_context(rc_polish()):
    fig = plt.figure(figsize=(FIGW, H))
    # round 49 nudge: at 12 pt the y labels of the left column crossed the 4 mm margin gate; each y label keeps its V26 left edge (the page-margin
    # side) and grows toward its tick labels: labelpad reduced by the growth of the widest tick label ("30", 10 -> 11 pt) plus the label's own
    # height growth (11 -> 12 pt), both measured in matplotlib's own metric (the same rule on all four panels)
    def _text_wh(s, size, rotation=0):
        fig.canvas.draw(); r_ = fig.canvas.get_renderer(); t_ = fig.text(0.5, 0.5, s, fontsize=size, rotation=rotation); bb_ = t_.get_window_extent(renderer=r_); t_.remove(); return bb_.width / fig.dpi, bb_.height / fig.dpi
    _d_tick = text_w(fig, "30", TICK_PT) - text_w(fig, "30", TICK_PT_V26)
    _d_lab = _text_wh("Percent of the group", TITLE_PT, 90)[0] - _text_wh("Percent of the group", TITLE_PT - DELTA, 90)[0]
    YLABELPAD_V26 = float(plt.rcParams["axes.labelpad"]); YLABELPAD = YLABELPAD_V26 - 72.0 * (_d_tick + _d_lab)
    assert YLABELPAD >= 1.0, (YLABELPAD, _d_tick, _d_lab)
    # round 49: the axes widths are the V26 ones (KEY_W measured at the V26 key size, 10 pt); the key text is set at 11 pt and the key stays right-aligned
    # to the edge, so its left edge moves into the 0.08 in gap by its growth (recorded and declared)
    _key_line_w = max(text_w(fig, ln, TICK_PT_V26) for s in LEG_WRAPPED for ln in s.split("\n"))
    KEY_W = (0.1 + 1.2 + 0.6 + 0.1) * TICK_PT_V26 / 72.0 + _key_line_w
    AXW_EACH = (W_IN - EDGE - KEY_W - KEY_GAP - X0_LEFT - MID) / 2.0
    X0 = [X0_LEFT, X0_LEFT + AXW_EACH + MID]; AXW = [AXW_EACH, AXW_EACH]
    # round 51 (2c): the one-line "PR x.xx (lo-hi)" tick-label line must fit under each pair: the pair pitch (half the axes width) at least the
    # widest line plus 6 pt, otherwise the two lines of a panel would touch. At 11 pt the line is about 96.6 pt and the V28 pitch 75.4 pt, so
    # the axes widen (and the page with them, the key and margins as before). The two-line alternative keeps the V28 width.
    AXW_EACH_V28 = AXW_EACH
    PR_LINES = [f"PR {GOT[(L_, c_)][4]:.2f} ({GOT[(L_, c_)][5]:.2f}-{GOT[(L_, c_)][6]:.2f})" for L_ in "ABCD" for c_ in (5, 10)]
    PR_W = max(text_w(fig, s_, TICK_PT) for s_ in PR_LINES)
    PR_GAP_PT = 9.0                                                                 # the gap between the two PR lines of a panel (one and a half spaces)
    if PR_LAYOUT == "one_line": AXW_EACH = max(AXW_EACH_V28, 2.0 * (PR_W + PR_GAP_PT / 72.0))
    W_IN = FIGW + 2.0 * (AXW_EACH - AXW_EACH_V28); fig.set_size_inches(W_IN, H)
    X0 = [X0_LEFT, X0_LEFT + AXW_EACH + MID]; AXW = [AXW_EACH, AXW_EACH]
    assert PR_LAYOUT == "two_lines" or AXW_EACH / 2.0 >= PR_W + PR_GAP_PT / 72.0, (AXW_EACH, PR_W)
    for idx, (title, frame, col, L) in enumerate(PANELS):
        cx, cy = idx % 2, idx // 2
        x0, axw = X0[cx], AXW[cx]; ytop = yrow[cy]
        ax = fig.add_axes([x0 / W_IN, 1 - (ytop + TITLE_BAND + AX_H) / H, axw / W_IN, AX_H / H])
        got = [GOT[(L, 5)], GOT[(L, 10)]]; xs = np.arange(2)
        ax.bar(xs - 0.19, [g[0] for g in got], 0.36, color=PRIMARY, linewidth=0)
        ax.bar(xs + 0.19, [g[1] for g in got], 0.36, color=DEEMPH, linewidth=0)
        STARS = {}
        for i, g in enumerate(got):   # round 51 (2a, 2c): no values above the bars, the P value as asterisks centred above the pair (*** P < 0.001, ** < 0.01, * < 0.05, none at or above 0.05)
            cut = (5, 10)[i]; st = splitstyle.stars(g[7]); STARS[cut] = st
            if st: ax.text(i, max(g[0], g[1]) + 0.8, st, ha="center", va="baseline", fontsize=ANN_PT, fontweight="bold", color=INK)
        PRL = [f"PR {g[4]:.2f} ({g[5]:.2f}-{g[6]:.2f})" for g in got]        # round 51 (2c): the prevalence ratio with its 95% CI under each pair
        if PR_LAYOUT == "one_line": labs = [f"{TICKLAB[i]}\n{PRL[i]}" for i in range(2)]
        else: labs = [f"{TICKLAB[i]}\nPR {got[i][4]:.2f}\n({got[i][5]:.2f}-{got[i][6]:.2f})" for i in range(2)]
        ax.set_xticks(xs); ax.set_xticklabels(labs, fontsize=TICK_PT, linespacing=1.35)
        ax.set_xlim(-0.5, 1.5); ax.set_ylim(0, 33); ax.set_yticks([0, 10, 20, 30])
        ax.set_ylabel("Percent of the group", fontsize=TITLE_PT, labelpad=YLABELPAD)      # round 49 nudge (see above)
        axis_covers(ax, [g[0] for g in got] + [g[1] for g in got], which="y", label=f"{SHEET} {L}")
        tx = x0 - 0.42
        letter(fig, L.lower(), tx, ytop)
        ptitle = title
        fig.text(tx / W_IN, 1 - (ytop + 0.24) / H, ptitle, ha="left", va="top", fontsize=PT_TITLE, color=INK, linespacing=1.4)
        for _ln in ptitle.split("\n"):
            assert tx + text_w(fig, _ln, PT_TITLE) < (X0[1] - 0.52 if cx == 0 else W_IN - EDGE), _ln
        DRAWN["panels"][L.lower()] = {"title": title, "n_short": got[0][2], "n_ref": got[0][3],
                                      "cuts": {c: {"pct_short": GOT[(L, c)][0], "pct_ref": GOT[(L, c)][1], "k_short": GOT[(L, c)][8], "k_ref": GOT[(L, c)][9],
                                                   "pr": GOT[(L, c)][4], "lo": GOT[(L, c)][5], "hi": GOT[(L, c)][6], "p": GOT[(L, c)][7], "p_printed": PSTR[(L, c)],
                                                   "printed": ([f"PR {GOT[(L, c)][4]:.2f} ({GOT[(L, c)][5]:.2f}-{GOT[(L, c)][6]:.2f})"] if PR_LAYOUT == "one_line" else [f"PR {GOT[(L, c)][4]:.2f}", f"({GOT[(L, c)][5]:.2f}-{GOT[(L, c)][6]:.2f})"]) + ([STARS[c]] if STARS[c] else []),
                                                   "stars": STARS[c], "removed_v28": [f"{GOT[(L, c)][0]:.1f}", f"{GOT[(L, c)][1]:.1f}", f"PR, {GOT[(L, c)][4]:.2f}", "(95% CI,", f"{GOT[(L, c)][5]:.2f}-{GOT[(L, c)][6]:.2f})", PSTR[(L, c)]]}
                                               for c in (5, 10)}}
    from matplotlib.patches import Patch
    _KEY_TOP_IN = HEAD_H + TITLE_BAND
    def _make_key18(x_in):
        return fig.legend(handles=[Patch(facecolor=PRIMARY, edgecolor="none", label=LEG_WRAPPED[0]), Patch(facecolor=DEEMPH, edgecolor="none", label=LEG_WRAPPED[1])],
                          loc="upper left", bbox_to_anchor=(x_in / W_IN, 1 - _KEY_TOP_IN / H), ncol=1, fontsize=TICK_PT, frameon=False, handletextpad=0.6,
                          labelspacing=0.7, handlelength=1.2, borderpad=0.1, borderaxespad=0)
    _leg18 = _make_key18(X0[1] + AXW[1] + KEY_GAP); fig.canvas.draw()
    _lbb = _leg18.get_window_extent(renderer=fig.canvas.get_renderer()); _leg18.remove()
    _leg18 = _make_key18(X0[1] + AXW[1] + KEY_GAP + ((W_IN - EDGE) - _lbb.x1 / fig.dpi)); fig.canvas.draw()
    _lbb = _leg18.get_window_extent(renderer=fig.canvas.get_renderer())
    KEY_LEFT_IN = _lbb.x0 / fig.dpi; KEY_GAP_IN = KEY_LEFT_IN - (X0[1] + AXW[1])
    # round 49: the key box (11 pt) is 0.104 in wider than the V26 one (10 pt); right-aligned to the page margin as before, its left edge crosses the
    # panel b axes' invisible right boundary by up to 0.03 in (the nearest bar ends 0.14 in inside that boundary, the percent labels 0.18 in)
    assert KEY_LEFT_IN >= X0[1] + AXW[1] - 0.03 and _lbb.x1 / fig.dpi <= W_IN - EDGE + 1e-6, (KEY_LEFT_IN, X0[1] + AXW[1], KEY_GAP)
    assert abs((H - _lbb.y1 / fig.dpi) - _KEY_TOP_IN) <= 0.03 and H - _lbb.y0 / fig.dpi <= _KEY_TOP_IN + AX_H
    assert abs(fig.get_size_inches()[0] - W_IN) < 1e-6
    H_V28 = HEAD_H + 2 * (TITLE_BAND + AX_H + BELOW_V28) + ROW_GAP + 0.06 + EDGE
    DRAWN["layout"] = dict(pr_layout=PR_LAYOUT, w_in=W_IN, h_in=H, w_v28_in=FIGW, h_v28_in=H_V28, axw_each_in=AXW_EACH, axw_each_v28_in=AXW_EACH_V28, pr_line_w_in=PR_W, pr_gap_pt=PR_GAP_PT, pair_pitch_pt=AXW_EACH * 36.0, pair_pitch_v28_pt=AXW_EACH_V28 * 36.0, below_in=BELOW, below_v28_in=BELOW_V28, bar_width_pt=0.36 * AXW_EACH * 36.0, bar_width_v28_pt=0.36 * AXW_EACH_V28 * 36.0)
    DRAWN["ylabel_pad"] = dict(v26_pt=YLABELPAD_V26, new_pt=round(YLABELPAD, 3), tick_growth_pt=round(72 * _d_tick, 3), label_growth_pt=round(72 * _d_lab, 3))
    DRAWN["key_placement"] = dict(key_left_in=KEY_LEFT_IN, key_right_in=_lbb.x1 / fig.dpi, axes_right_in=X0[1] + AXW[1], gap_in=KEY_GAP_IN, gap_v26_in=KEY_GAP, key_top_in=H - _lbb.y1 / fig.dpi, axw_each_in=AXW_EACH, key_w_layout_in=KEY_W)
    SB.text_gate(fig, SHEET, {"a", "b", "c", "d", "*", "**", "***"}); SB.nooverlap.gate(fig, SHEET); SB.margin_gate(fig, SHEET)
    raw = f"{OUTD}/{SHEET}_raw.pdf"; fig.savefig(raw); plt.close(fig)
final = f"{OUTD}/{SHEET}_built.pdf"                      # the strip is stamped afterwards into HERE/ED_Fig05.pdf
DRAWN["font_metrics_aligned"] = C.align_font_metrics(raw, final, BASE)
# drawn records from the final sheet's text layer: the 32 value strings and the 8 '(95% CI,' lines, panel by panel
d = fitz.open(final); sp = C.spans_of(d[0]); d.close()
def panel_of(s): return ("a" if s["origin"][1] < 310 else "c") if s["origin"][0] < 200 else ("b" if s["origin"][1] < 310 else "d")
RX_PR1 = re.compile(r"^PR \d\.\d\d \(\d\.\d\d-\d\.\d\d\)$"); RX_PR2 = re.compile(r"^PR \d\.\d\d$"); RX_CI2 = re.compile(r"^\(\d\.\d\d-\d\.\d\d\)$"); RX_ST = re.compile(r"^\*{1,3}$")
def panel_of(s):   # round 51: the panels' boxes from the layout (the page may be wider than V28)
    xmid = (X0[0] + AXW[0] + X0[1]) / 2.0 * 72.0; ymid = (yrow[0] + TITLE_BAND + AX_H + BELOW + yrow[1]) / 2.0 * 72.0
    return ("a" if s["origin"][1] < ymid else "c") if s["origin"][0] < xmid else ("b" if s["origin"][1] < ymid else "d")
for pnl in "abcd":
    pn = DRAWN["panels"][pnl]; ss = sorted([s for s in sp if panel_of(s) == pnl], key=lambda s: (round(s["origin"][1], 1), s["origin"][0]))
    prs = sorted([s for s in ss if RX_PR1.match(s["text"].strip()) or RX_PR2.match(s["text"].strip())], key=lambda s: s["origin"][0])
    cis = sorted([s for s in ss if RX_CI2.match(s["text"].strip())], key=lambda s: s["origin"][0]); sts = sorted([s for s in ss if RX_ST.match(s["text"].strip())], key=lambda s: s["origin"][0])
    n_st = sum(1 for c in (5, 10) if pn["cuts"][c]["stars"])
    assert len(prs) == 2 and len(cis) == (0 if PR_LAYOUT == "one_line" else 2) and len(sts) == n_st, (pnl, len(prs), len(cis), len(sts), n_st)
    cxs = [(p["bbox"][0] + p["bbox"][2]) / 2 for p in prs]
    for i, cut in enumerate((5, 10)):
        v = pn["cuts"][cut]; base = dict(panel=pnl, source_file=PARQ, note=f"panel {pnl} T90 >{cut}%: recomputed from the v8 table (cohort_spec, 15,551) = step-129 log")
        def rec(s, key, val, rule): DRAWN["records"].append(dict(base, text=s["text"].strip(), x=(s["bbox"][0] + s["bbox"][2]) / 2, baseline=s["origin"][1], ha="center", size=s["size"], source_key=key, source_value=val, rule=rule))
        rec(prs[i], f"prevalence ratio with its 95% CI, T90 >{cut}%", [v["pr"], v["lo"], v["hi"]], "pr_ci_line" if PR_LAYOUT == "one_line" else "pr_2dp_first")
        if cis: rec(cis[i], f"95% CI of the prevalence ratio T90 >{cut}%", [v["lo"], v["hi"]], "ci_2dp_paren")
        st_i = [s for s in sts if abs((s["bbox"][0] + s["bbox"][2]) / 2 - cxs[i]) < abs((s["bbox"][0] + s["bbox"][2]) / 2 - cxs[1 - i])]
        assert len(st_i) == (1 if v["stars"] else 0), (pnl, cut, v["stars"], len(st_i))
        if st_i: rec(st_i[0], f"two-sided P of the prevalence ratio T90 >{cut}% as asterisks (P = {v['p']:.3g})", v["p"], "stars_p")
VALUE_RX = re.compile(r"^(PR \d\.\d\d.*|\(\d\.\d\d-\d\.\d\d\)|\*{1,3})$")
for s in sp:   # round 49: every other string (titles, letters, tick labels, y labels, key lines, '(95% CI,' lines) recorded as the V13 string at +1 pt for the printed-values CSV
    t = s["text"].strip()
    if VALUE_RX.match(t) or not t: continue
    DRAWN["records"].append(dict(text=t, x=s["origin"][0], baseline=s["origin"][1], ha="left", size=s["size"], panel=panel_of(s), source_file="static:V13", source_key="V13 string at +1 pt, grown about its V13 anchor", source_value=t, rule="text"))
json.dump(DRAWN, open(f"{OUTD}/{SHEET}_drawn.json", "w"), indent=1, default=float)
for L in "ABCD":
    for c in (5, 10):
        g = GOT[(L, c)]; print(f"  {L} T90>{c:>2}%  {g[0]:5.1f} vs {g[1]:5.1f}  n {g[2]:,} vs {g[3]:,}  PR {g[4]:.2f} ({g[5]:.2f}-{g[6]:.2f}) {PSTR[(L, c)]}")
print("wrote", final, f"{W_IN * 72:.2f} x {H * 72:.2f} pt", DRAWN["layout"])
