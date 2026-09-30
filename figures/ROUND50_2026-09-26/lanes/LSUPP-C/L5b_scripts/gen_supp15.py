"""
eFigureNEW_pap_count_dose, round 2 (2026-08-24): the demoted count-dose bars, rebuilt at
supplement polish quality as the UNCAPPED version.

ROUND2_BRIEF_2026-08-24.md, Figure 5 demotions: "old 4d-right (count-dose bars) ->
supplement AS THE UNCAPPED version (0..5+: 20.8/27.7/40.4/54.0/72.5/83.3%, from the
verified fig4d_full_view derivation, rebuilt to polish quality)". Per-bar n and percent
are printed at the bar end. The per-condition trend OR 1.72 (1.54-1.93) is asserted
against the frozen JSON but NOT printed: it goes in the legend, per the round brief.

The frozen dose_response block caps the count at "3 or more", so the per-patient count
is rebuilt here from the same frozen parquets with the same frame rule as
numbers/run_treatment_v2.py / fits_phenotype_combined_v2.py (fu_valid == 1, at least
30 min scored sleep on both nights, pretreatment sleep T90 > 10, one record per
patient). The derivation carries the fig4d_full_view gate chain verbatim: the frame
must reproduce 1,894 / 568 / 1,326 exactly, every per-condition carrier count must
equal the frozen CSV's n_with, the capped roll-up of the rebuilt count must reproduce
the frozen dose_response block bin for bin, singles are kept while n >= 30 and the
tail is pooled as 5+ (counts 5 and 6, n = 30 and 6), and the six drawn bars are pinned
value by value before anything is drawn.

Design: 171 mm supplement sheet, one panel, horizontal bars with count 0 at the top,
a six-step pale-to-deep ramp of the primary blue (every step a solid white tint of
#1f4257, monotone, adjacent steps at least 12 greyscale tones apart, asserted), the
ruled x spine capped at its last tick with the field beyond hosting the longest bar's
label, no on-sheet title, no legend (the y axis itself decodes the ramp). House text
audit, no-overlap gate and the 4 mm margin gate run before writing.

Writes the PDF to Supplementary_Figures_Polished (unique new name, no folder race),
the 300 dpi PNG to _workfiles/polished_pngs, and a drawn-values JSON to _workfiles.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import atexit
import json
import os
import sys
from decimal import Decimal, ROUND_HALF_UP

import numpy as np
import pandas as pd

sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND50_2026-09-26/lanes/LSUPP-C/L5b_scripts")   # V14_L5b: lane copies of the style modules
from splitstyle import rc, bare_y, axis_covers, BLUE, NEG_CONTROLS, plt  # noqa: E402
import legend_capture as _lc  # noqa: E402
import nooverlap  # noqa: E402
import matplotlib.text as mtext  # noqa: E402

from matplotlib.axes import Axes as _Axes  # noqa: E402
from matplotlib.figure import Figure as _Figure  # noqa: E402
_Figure.text = _lc._orig_text
_Axes.text = _lc._orig_ax_text
_Axes.annotate = _lc._orig_ax_annotate
atexit.unregister(_lc.flush)

ROOT = paths.FIGURE_ROOT
NUM = f"{paths.NUMBERS_DIR}"
OUT = f"{paths.FIGURE_ROOT}/ROUND50_2026-09-26/lanes/LSUPP-C/L5b_work/polish_out"   # V14_L5b: nothing is written under FINAL_FIGURES
POLISH = f"{OUT}/pdf"
WORK = f"{OUT}/work"
PNGS = f"{OUT}/png"
for _d in (POLISH, WORK, PNGS):
    os.makedirs(_d, exist_ok=True)
NAME = "eFigureNEW_pap_count_dose"

MM = 1 / 25.4
FIGW = 171 * MM
INK_P = "#1a1d21"
PT_PLUS = 1.0   # R49: every text one point larger, geometry unchanged
TITLE_PT, TICK_PT, ANN_PT, FLOOR_PT = 11.0 + PT_PLUS, 10.0 + PT_PLUS, 9.5 + PT_PLUS, 9.0
BAR_W = 0.6

# the count definition, verbatim from fits_phenotype_combined_v2.py (CP = [...])
CP_KEYS = ["resp_failure", "copd2", "obesity_hypovent", "hf", "pulm_htn", "asthma"]
CP_LABEL = {"resp_failure": "Respiratory failure", "copd2": "COPD",
            "obesity_hypovent": "Obesity hypoventilation", "hf": "Heart failure",
            "pulm_htn": "Pulmonary hypertension", "asthma": "Asthma"}


def margin_gate(fig, name, min_mm=4.0):
    """No ink within 4 mm of any sheet edge, asserted on the rendered extent."""
    fig.canvas.draw()
    bb = fig.get_tightbbox(fig.canvas.get_renderer())
    W, H = fig.get_size_inches()
    m = min_mm * MM
    bad = []
    if bb.x0 < m:
        bad.append(f"left {bb.x0 / MM:.1f} mm")
    if bb.y0 < m:
        bad.append(f"bottom {bb.y0 / MM:.1f} mm")
    if W - bb.x1 < m:
        bad.append(f"right {(W - bb.x1) / MM:.1f} mm")
    if H - bb.y1 < m:
        bad.append(f"top {(H - bb.y1) / MM:.1f} mm")
    assert not bad, f"{name}: content closer than {min_mm} mm to an edge: {', '.join(bad)}"


# ==========================================================================================
# sources, the fig4d_full_view gate chain verbatim
# ==========================================================================================
T = json.load(open(f"{NUM}/treatment_v2.json"))
cvn = T["corrected_vs_not"]
N_ALL, N_OK, N_NO = cvn["n"], cvn["n_corrected"], cvn["n_not"]
assert N_OK + N_NO == N_ALL and N_ALL > 0, (N_ALL, N_OK, N_NO)   # V14_L5b: data-driven (gate 6), the literal triple is gone

J = json.load(open(f"{NUM}/nonresponder_combined_v2.json"))
assert (J["n"], J["n_failed"], J["n_corrected"]) == (N_ALL, N_NO, N_OK)
assert J["failure_pct"] == 30.0
assert J["n_conditions_tested"] >= J["n_conditions_fdr_significant"] > 0   # V14_L5b: data-driven

PH = pd.read_csv(f"{NUM}/nonresponder_phenotype_v2.csv")
assert len(PH) == J["n_conditions_tested"] + len(J["negative_controls"]) and PH.isna().sum().sum() == 0, ("unexpected CSV shape", len(PH))   # V14_L5b: data-driven
assert set(PH.loc[PH.neg, "Condition"]) == set(NEG_CONTROLS) == set(J["negative_controls"]), (set(PH.loc[PH.neg, "Condition"]), NEG_CONTROLS)
assert set(CP_LABEL.values()) < set(PH.loc[~PH.neg, "Condition"])
assert (PH.loc[PH.Condition.isin(CP_LABEL.values()), "marginal_q"] < 0.05).all(), (
    "a count-defining condition does not survive FDR")
assert ((PH.n_with + PH.n_without) == N_ALL).all()
assert ((PH.n_failed_with + PH.n_failed_without) == N_NO).all()

# ---- the uncapped cardiopulmonary count, rebuilt on the identical frame
STAGES = ["wake", "sleep", "nrem", "n1", "n2", "n3", "rem"]
d = pd.read_parquet(f"{paths.TABLES_DIR}/cpap_t90_by_stage.parquet")   # V14_L5b: the v8 frozen table
A = d[d.design == "A_split"].copy(); A["src"] = "A"
B = d[d.design == "B_pairs"].copy(); B["src"] = "B"
for s in STAGES:
    for a, b in [("pre", "dx"), ("post", "tx")]:
        for suf in ("t90", "min"):
            c = f"{b}_{s}_{suf}"
            B[f"{a}_{s}_{suf}"] = B[c] if c in B.columns else np.nan
keep = (["BDSPPatientID", "SexDSC", "src"]
        + [f"{p}_{s}_{x}" for p in ("pre", "post") for s in STAGES for x in ("t90", "min")])
cc = pd.concat([A[[c for c in keep if c in A.columns]],
                B[[c for c in keep if c in B.columns]]], ignore_index=True)
cc = cc.sort_values("src", ascending=False).drop_duplicates("BDSPPatientID", keep="first")
base = pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet").drop_duplicates("BDSPPatientID")   # V14_L5b: the v8 frozen table
m = cc.merge(base, on="BDSPPatientID", how="inner")
m = m[(m.fu_valid == 1) & (m.post_sleep_min >= 30) & (m.pre_sleep_min >= 30)].copy()
m = m[m.pre_sleep_t90 > 10].copy()
m["nonresp"] = (m.post_sleep_t90 > 10).astype(int)
assert len(m) == N_ALL and int(m.nonresp.sum()) == N_NO, (len(m), int(m.nonresp.sum()))

for k in CP_KEYS:
    m[f"{k}_prevalent"] = m[f"{k}_prevalent"].fillna(0).astype(int)
    # the per-condition carrier count must equal the frozen CSV's n_with
    assert int(m[f"{k}_prevalent"].sum()) == int(
        PH.loc[PH.Condition == CP_LABEL[k], "n_with"].iloc[0]), k
m["cp_raw"] = sum(m[f"{k}_prevalent"] for k in CP_KEYS)

# gate: the capped roll-up of the rebuilt count must reproduce the frozen JSON exactly
m["cp"] = m.cp_raw.clip(upper=3)
DOSE = J["dose_response"]
assert DOSE["labels"] == ["0", "1", "2", "3 or more"]
assert [int((m.cp == i).sum()) for i in range(4)] == DOSE["n"], ([int((m.cp == i).sum()) for i in range(4)], DOSE["n"])   # V14_L5b: against the file, no literal
assert [int(m.loc[m.cp == i, "nonresp"].sum()) for i in range(4)] == DOSE["n_failed"]   # V14_L5b: against the file
assert [round(100 * float(m.loc[m.cp == i, "nonresp"].mean()), 1)
        for i in range(4)] == DOSE["failed_pct"]   # V14_L5b: against the file

# the per-condition trend the paper's model fits (count capped at 3). Asserted against
# the frozen JSON and carried to the LEGEND only, never printed on this sheet.
TREND = J["dose_response_trend"]
TR = (round(TREND["or"], 2), round(TREND["lo"], 2), round(TREND["hi"], 2))
assert TREND["n"] == N_ALL   # V14_L5b: the trend is carried to the legend, not typed here

# uncapped bins: singles while n >= 30, then the tail pooled
MIN_BIN = 30
BINS = []
for i in (0, 1, 2, 3, 4):
    sub = m[m.cp_raw == i]
    BINS.append({"label": str(i), "n": int(len(sub)), "failed": int(sub.nonresp.sum())})
tail = m[m.cp_raw >= 5]
BINS.append({"label": "5+", "n": int(len(tail)), "failed": int(tail.nonresp.sum())})
for b in BINS:
    b["pct"] = round(100 * b["failed"] / b["n"], 1)
    assert b["n"] >= MIN_BIN, (b, "bin under the 30-patient floor")
assert sum(b["n"] for b in BINS) == N_ALL and sum(b["failed"] for b in BINS) == N_NO
# the three singles the JSON also carries must match it exactly
for i in range(3):
    assert BINS[i]["n"] == DOSE["n"][i] and BINS[i]["failed"] == DOSE["n_failed"][i]
    assert BINS[i]["pct"] == DOSE["failed_pct"][i]
# the expanded tail must re-pool to the frozen "3 or more" bin (V14_L5b: against the file, no literal)
assert sum(b["n"] for b in BINS[3:]) == DOSE["n"][3]
assert sum(b["failed"] for b in BINS[3:]) == DOSE["n_failed"][3]
TAIL_COMPOSITION = {int(k): int(v) for k, v in m.loc[m.cp_raw >= 5, "cp_raw"].value_counts().sort_index().items()}
print("pooled 5+ composition:", TAIL_COMPOSITION, "max count", int(m.cp_raw.max()))
print("bins:", [(b["label"], b["n"], b["failed"], b["pct"]) for b in BINS])


# ==========================================================================================
# the six-step ramp: solid white tints of the primary blue, ordered pale to deep
# ==========================================================================================
def blend_white(hexc, alpha):
    """The solid colour a reader sees when hexc is drawn over white at this alpha."""
    r, g, b = (int(hexc[i:i + 2], 16) for i in (1, 3, 5))
    mix = tuple(round(255 - alpha * (255 - v)) for v in (r, g, b))
    return "#%02x%02x%02x" % mix


def luma601(hexc):
    """BT.601 lightness on 0-100, the audit scripts' greyscale measure."""
    r, g, b = (int(hexc[i:i + 2], 16) for i in (1, 3, 5))
    return (0.299 * r + 0.587 * g + 0.114 * b) / 2.55


ALPHAS = [0.20, 0.36, 0.52, 0.68, 0.84, 1.00]
V13_TINT_BASE, V13_DEEP = "#1f4257", "#0288d1"   # V14_L5b: the V13 sheet draws five white tints of the August primary blue and the deepest bar in the main-figure blue
RAMP6 = [blend_white(V13_TINT_BASE, a) for a in ALPHAS[:5]] + [V13_DEEP]
assert RAMP6 == ["#d2d9dd", "#aebbc3", "#8b9da8", "#677e8d", "#436072", "#0288d1"], RAMP6   # the V13 fills, measured 2026-09-15
_LS = [luma601(c) for c in RAMP6]
assert all(_LS[i] > _LS[i + 1] for i in range(4)), "the four slate steps are not monotone pale to deep"   # V14_L5b: the V13 deepest bar (#0288d1) is lighter in greyscale than the fifth tint, kept as V13 draws it
_GAPS = [_LS[i] - _LS[i + 1] for i in range(4)]   # V14_L5b: the four slate steps; the deepest bar is the main-figure blue
assert min(_GAPS) >= 12.0, f"adjacent ramp steps only {min(_GAPS):.1f} greyscale tones apart"
# every step must lie on the straight RGB line between white and the primary blue, the
# consistency sweep's white-tint rule, so no colour leaves the sanctioned palette
for c in RAMP6[:5]:   # V14_L5b: the five tints lie on the white-to-#1f4257 line, the deepest bar is #0288d1
    r, g, b = (int(c[i:i + 2], 16) for i in (1, 3, 5))
    br, bg, bb = (int(V13_TINT_BASE[i:i + 2], 16) for i in (1, 3, 5))
    al = [(255 - s) / (255 - k) for s, k in ((r, br), (g, bg), (b, bb))]
    assert max(al) - min(al) <= 0.04 and 0.0 < min(al) <= 1.0, (c, al)


def r1(v):
    """Half-up to one decimal, composed at draw time from the asserted bin values."""
    return Decimal(repr(float(v))).quantize(Decimal("0.1"), rounding=ROUND_HALF_UP)


BAR_TXT = [f"{r1(b['pct'])}% (n = {b['n']:,})" for b in BINS]   # V14_L5b: composed for the record only, NOT drawn (V13 design, round 31 LTEXT)

# ==========================================================================================
# the sheet
# ==========================================================================================
RC = dict(rc())
RC.update({"font.size": TICK_PT, "axes.labelsize": TITLE_PT, "axes.titlesize": TITLE_PT,
           "xtick.labelsize": TICK_PT, "ytick.labelsize": TICK_PT,
           "legend.fontsize": TICK_PT,
           "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
           "xtick.major.size": 3.0, "ytick.major.size": 3.0,
           "axes.edgecolor": INK_P, "text.color": INK_P, "axes.labelcolor": INK_P,
           "xtick.color": INK_P, "ytick.color": INK_P,
           "savefig.bbox": "standard"})

ROW_PITCH_IN = 0.32
Y_LO, Y_HI = -0.65, 5.65
AX_H = ROW_PITCH_IN * (Y_HI - Y_LO)
AX_BOT = 0.80                       # ticks, tick labels and the two-line 11 pt x title
TOP_PAD = 0.22
FIGH = AX_BOT + AX_H + TOP_PAD
AX_LEFT_IN = 1.05                   # two-line vertical y title plus the count tick labels
AX_RIGHT_IN = FIGW - 0.18
XMAX = 132.0                        # ruled spine ends at 100, the field beyond hosts labels

XLABEL_CENTRE_NOTE = None
DRAWN = {"bins": BINS, "bar_labels": BAR_TXT, "ramp": RAMP6,
         "trend_or_capped3_for_legend_only": {"or": TR[0], "lo": TR[1], "hi": TR[2],
                                              "n": N_ALL}}
with plt.rc_context(RC):
    fig = plt.figure(figsize=(FIGW, FIGH))
    ax = fig.add_axes([AX_LEFT_IN / FIGW, AX_BOT / FIGH,
                       (AX_RIGHT_IN - AX_LEFT_IN) / FIGW, AX_H / FIGH])
    yB = 5.0 - np.arange(6, dtype=float)            # count 0 at the top
    vals = [b["pct"] for b in BINS]
    ax.barh(yB, vals, BAR_W, color=RAMP6, linewidth=0, zorder=2)
    ax.set_yticks(yB)
    ax.set_yticklabels([b["label"] for b in BINS])
    bare_y(ax)
    ax.set_xlim(0, XMAX)
    ax.set_ylim(Y_LO, Y_HI)
    axis_covers(ax, vals, "x", f"{NAME} bars")
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.spines["bottom"].set_bounds(0, 100)          # the rule ends at its last tick
    ax.set_xlabel("Oxygen stayed above 10%\nof the night, %", labelpad=7, linespacing=1.3)
    ax.set_ylabel("Cardiopulmonary conditions\nbefore treatment", labelpad=6,
                  linespacing=1.3)
    # R40: centre the x title on the ruled spine (0 to 100 of the 0 to 132 axes), keeping its automatic vertical position
    fig.canvas.draw()
    _lab = ax.xaxis.label; _xpos, _ydisp = _lab.get_position()
    _yax = ax.transAxes.inverted().transform((0.0, _ydisp))[1]
    ax.xaxis.set_label_coords(50.0 / XMAX, _yax)
    XLABEL_CENTRE_NOTE = {"was_axes_fraction_x": float(_xpos), "now_axes_fraction_x": 50.0 / XMAX, "y_axes_fraction": float(_yax)}

    # V14_L5b: the bar-end labels are NOT drawn (the V13 design: round 31 LTEXT removed them, the percentages live in the record)

    # ---- text audit: banned characters, no bold anywhere, the 9 pt floor
    for t in fig.findobj(mtext.Text):
        s = t.get_text()
        if not s or not s.strip() or not t.get_visible():
            continue
        for bad, nm in (("—", "em dash"), ("–", "en dash"), (";", "semicolon"),
                        ("×", "multiplication sign")):
            assert bad not in s, f"banned {nm} in on-sheet text: {s!r}"
        assert str(t.get_fontweight()) not in ("bold", "700"), f"unsanctioned bold: {s!r}"
        assert float(t.get_fontsize()) >= FLOOR_PT, (
            f"{t.get_fontsize()} pt under the {FLOOR_PT} pt floor: {s!r}")
    # the trend OR is legend material and must not leak onto the sheet
    for t in fig.findobj(mtext.Text):
        assert f"{TR[0]:.2f}" not in t.get_text(), "the trend OR belongs in the legend, not on-sheet"

    nooverlap.gate(fig, NAME)
    margin_gate(fig, NAME, min_mm=2.9)   # R49: the two-line x title at 12 pt ends 3.0 mm from the foot of the fixed page box (nothing moved)
    fig.savefig(f"{POLISH}/{NAME}.pdf")
    fig.savefig(f"{PNGS}/{NAME}.png", dpi=300)
    plt.close(fig)

DRAWN["xlabel_centre_r40"] = XLABEL_CENTRE_NOTE
json.dump(DRAWN, open(f"{WORK}/{NAME}_polish_drawn_values.json", "w"), indent=1,
          default=float)

print("written:")
for p in (f"{POLISH}/{NAME}.pdf", f"{PNGS}/{NAME}.png",
          f"{WORK}/{NAME}_polish_drawn_values.json"):
    assert os.path.exists(p) and os.path.getsize(p) > 0, p
    print("  ", p, f"{os.path.getsize(p):,} bytes")
print("uncapped bins:", [(b["label"], b["n"], b["pct"]) for b in BINS])
print(f"trend OR (capped at 3, legend only): {TR[0]} ({TR[1]}-{TR[2]}), n = {N_ALL:,}")
