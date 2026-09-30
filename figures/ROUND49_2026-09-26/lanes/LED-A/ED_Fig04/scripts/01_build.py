#!/usr/bin/env python3
"""ROUND 49 lane LED-A copy (2026-09-26): +1 pt through the copied helpers (figure2_split_common BASE/TICKF/SMALLF, splitstyle), the key anchored by its
centre on the forest axes (the V13 key row centre) so that it grows symmetrically, output work/ED_Fig04_built.pdf (the strip is stamped afterwards).
Original: ROUND37_2026-09-14/figures/lanes/V14_L3a_DURATION_MAIN/ED_Fig04/scripts/01_build.py.
ED Fig 4 (single-panel sheet), V14 rebuild on the v8.1 numbers at the V13 design. Repointed copy of
ROUND30_2026-09-08/V7_L3a_DURATION_MAIN/ED_Fig04/scripts/01_build.py (backup: 01_build_PRE_V8_1.py). Design = efigure24_polish.py
(FINAL_FIGURES_2026-08-14/_scripts) as it shows on the V13 sheet (byte-identical to the V7 sheet), data = numbers/results_v2.json
["vs_sleep_duration"] (step 112, v8.1 sidecar: the T90 term on all 19,173 nights, the sleep-under-5-h term on the 15,551 full nights).
Colours sampled from the V13 sheet (blue #0288d1 for T90, green #39c445 markers with #298d32 CI lines for laboratory sleep). Key
centred under the axis title as on the V13 sheet (two-pass placement measured through the PDF text layer). Row set = the ten
largest T90-significant (P < 0.05, HR > 1) non-control conditions, the R30 builder's rule, controls from the v8.1 panel (this lane's
splitstyle copy reads numbers/negcontrols_final.json, step 205). Writes ED_Fig04.pdf, work/ED_Fig04_raw.pdf, work/ED_Fig04_drawn.json.
V14 changes against the R30 copy: paths (LANE, BASE = V13, sidecar gate v8.1), the drawn record carries the v14lib DRAWN_RECORD fields.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import json, os, sys
import numpy as np
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LED-A/scripts_shared")   # round 49: the lane's patched helpers (+1 pt)
import lane_common as C
C.load_ff()
import figure2_polish_common as F2
from figure2_polish_common import (rc_polish, logx_clean, bare_y, axis_covers, left_labels, fit_gutter, ref_rule,
                                   INK, plt, NEG_CONTROLS, W_IN, TICKF, AX_X0, FOREST_W, EDGE, CI_LW, MEDGE, S_CIRCLE,
                                   S_SQUARE, MS_CIRCLE, MS_SQUARE)
from matplotlib.lines import Line2D
import fitz

SHEET = "ED_Fig04"; HERE = f"{C.LANE}/{SHEET}"; WORK = f"{HERE}/work"; VER = f"{HERE}/verify"
os.makedirs(WORK, exist_ok=True); os.makedirs(VER, exist_ok=True)
SRC = f"{C.NUM}/results_v2.json"
BASE = f"{C.BASE}/{SHEET}.pdf"
C.hydrated(SRC); C.hydrated(BASE)
SC = C.sidecar(SRC); SHA = SC["sha256"]
assert sorted(NEG_CONTROLS) == sorted(["Alopecia", "Inguinal hernia", "Glaucoma", "Contact dermatitis", "Hemorrhoids"]), NEG_CONTROLS   # v8.1 panel (read from the json)
BLUE, GREEN_MK, GREEN_CI = "#0288d1", "#39c445", "#298d32"       # sampled from the V13 sheet's drawings

R = json.load(open(SRC))["vs_sleep_duration"]
VSO = R["outcomes"]
assert R["n"] == 19173, R["n"]
for k, v in VSO.items(): assert v["negative_control"] == (k in NEG_CONTROLS), (k, v["negative_control"])
DCONDS = [k for _, k in sorted(((v["t90"]["hr"], k) for k, v in VSO.items() if v["t90"]["p"] < 0.05 and v["t90"]["hr"] > 1), reverse=True)
          if k not in NEG_CONTROLS]
N_D = len(DCONDS)
N_OX_LARGER = sum(1 for k in DCONDS if VSO[k]["t90"]["hr"] > VSO[k]["short_sleep"]["hr"])
TOPC = DCONDS[:10]
assert len(TOPC) == 10 and len(set(TOPC)) == 10, TOPC
_hrs = [VSO[k]["t90"]["hr"] for k in TOPC]
assert all(a >= b for a, b in zip(_hrs, _hrs[1:])), "rows not sorted by the T90 hazard ratio"

# the V13 sheet's anchors (text layer): key strings and their visible positions
bd = fitz.open(BASE); bp = bd[0]; base_spans = C.spans_of(bp); bd.close()
key1 = [s for s in base_spans if s["text"] == "T90 >10%" and s["origin"][0] > 200]
key2 = [s for s in base_spans if s["text"] == "Laboratory sleep under 5 h" and s["origin"][0] > 250]
assert len(key1) == 1 and len(key2) == 1, (key1, key2)
KEY1_ORG, KEY2_ORG = key1[0]["origin"], key2[0]["origin"]
BASE_LABELS = sorted({s["text"] for s in base_spans if abs(s["size"] - 10.0) < 0.05 and s["origin"][0] < 150 and s["origin"][1] < 290})

C_H = 3.75; C_Y0 = 1.00; H_IN = C_Y0 + C_H + EDGE; OFF = 0.185; XTICKS_C = [1, 2, 4, 8]
XLIM = (0.68, 12.3); YLIM = (-0.62, 9.80)


def build(leg_anchor, out_pdf):
    DRAWN = {"figure": SHEET, "source": SRC, "source_sha256": SHA, "step": SC, "base_sheet": BASE, "base_sha256": C.sha256(BASE),
             "t90_larger_in": f"{N_OX_LARGER} of {N_D}", "n_t90_significant_conditions": N_D,
             "rows": [], "colours": {"t90": BLUE, "short_sleep_marker": GREEN_MK, "short_sleep_ci": GREEN_CI},
             "cohort_n": R["n"], "pct_t90_hi": R["pct_t90_hi"], "pct_short": R["pct_short"], "base_row_labels": BASE_LABELS, "records": []}
    with plt.rc_context(rc_polish()):
        fig = plt.figure(figsize=(W_IN, H_IN))
        fit_gutter(fig, TOPC)
        axC = fig.add_axes([AX_X0 / W_IN, C_Y0 / H_IN, FOREST_W / W_IN, C_H / H_IN])
        yC = np.arange(len(TOPC), dtype=float)[::-1]
        ref_rule(axC, 1.0)
        drawnC = []
        for k, yy in zip(TOPC, yC):
            v = VSO[k]; t, s = v["t90"], v["short_sleep"]
            axC.plot([t["lo"], t["hi"]], [yy + OFF] * 2, color=BLUE, lw=CI_LW, solid_capstyle="round", zorder=2)
            axC.scatter([t["hr"]], [yy + OFF], s=S_CIRCLE, color=BLUE, zorder=3, edgecolors="white", linewidths=MEDGE)
            axC.plot([s["lo"], s["hi"]], [yy - OFF] * 2, color=GREEN_CI, lw=CI_LW, solid_capstyle="round", zorder=2)
            axC.scatter([s["hr"]], [yy - OFF], s=S_SQUARE, marker="s", color=GREEN_MK, zorder=3, edgecolors="white", linewidths=MEDGE)
            drawnC += [t["hr"], t["lo"], t["hi"], s["hr"], s["lo"], s["hi"]]
            DRAWN["rows"].append({"condition": k, "row_y": float(yy), "t90": dict(t), "short_sleep": dict(s), "events": v["events"],
                                  "key_path": f"vs_sleep_duration.outcomes.{k}"})
        left_labels(axC, yC, TOPC)
        bare_y(axC)
        logx_clean(axC, XTICKS_C)
        axC.set_xlim(*XLIM); axC.set_ylim(*YLIM)
        axis_covers(axC, drawnC, "x", SHEET)
        axC.set_xlabel("Hazard ratio, both exposures in one model (95% CI)", labelpad=4)
        lg = fig.legend(handles=[
            Line2D([], [], color=BLUE, marker="o", ms=MS_CIRCLE, lw=CI_LW, markeredgecolor="white", markeredgewidth=MEDGE, label="T90 >10%"),
            Line2D([], [], color=GREEN_CI, marker="s", ms=MS_SQUARE, lw=CI_LW, markerfacecolor=GREEN_MK, markeredgecolor="white",
                   markeredgewidth=MEDGE, label="Laboratory sleep under 5 h")],
            loc="lower center", bbox_to_anchor=leg_anchor, ncol=2, fontsize=TICKF, frameon=False, handletextpad=0.5,      # round 49: anchored by its centre (see below)
            columnspacing=1.6, labelspacing=0.42, handlelength=1.5)
        # the gates of save_polished, without its canonical write location
        F2._prose_gate(fig, SHEET)
        F2.nooverlap.gate(fig, SHEET)
        problems = [p for p in F2.preflight_figure(fig, F2.JSTYLE) if "height" not in p.lower()]
        assert not problems, problems
        fig.savefig(out_pdf)
        plt.close(fig)
    return DRAWN


def _hx(c): return "#%02x%02x%02x" % tuple(int(round(v * 255)) for v in c)
def key_geom(pdf):
    """The key's two text origins and its ink box: left = the left end of the first handle line (blue, CI_LW, in the key band), right = the second label's right edge."""
    d = fitz.open(pdf); sp = C.spans_of(d[0]); dr = d[0].get_drawings(); d.close()
    k1 = [s for s in sp if s["text"] == "T90 >10%"]; k2 = [s for s in sp if s["text"] == "Laboratory sleep under 5 h"]
    assert len(k1) == 1 and len(k2) == 1, (k1, k2)
    yb = k1[0]["origin"][1]
    hl = [g for g in dr if g.get("color") is not None and _hx(g["color"]) == BLUE and abs((g.get("width") or 0) - CI_LW) < 0.01 and abs(g["rect"].y0 - yb) < 12 and g["rect"].width > 5 and g["rect"].x0 < k1[0]["origin"][0]]
    assert len(hl) == 1, hl
    return k1[0]["origin"], k2[0]["origin"], hl[0]["rect"].x0, k2[0]["bbox"][2]
def key_origins(pdf): return key_geom(pdf)[:2]


# round 49: the V13 key is centred under the forest axes (its row centre 312.5 pt = the axes centre 312.4 pt), so at +1 pt it grows about that
# centre: the legend is anchored "lower center" on the axes centre and its baseline iterated onto the V13 key baseline; the first entry moves
# left and the second right by half the row's growth (asserted symmetric, declared in the record)
# the V13 key's ink box (first handle's left end to the second label's right edge) and its centre, read from the V13 sheet itself
_, _, KL_13, KR_13 = key_geom(BASE); K13_C = (KL_13 + KR_13) / 2.0
anchor = [K13_C / (W_IN * 72), (EDGE + 0.06) / H_IN]
raw = f"{WORK}/{SHEET}_raw.pdf"
for it in range(6):
    DRAWN = build(tuple(anchor), raw)
    o1, o2, kl, kr = key_geom(raw)
    dx = K13_C - (kl + kr) / 2.0; dy = KEY1_ORG[1] - o1[1]
    if abs(dx) < 0.02 and abs(dy) < 0.02:
        break
    anchor = [anchor[0] + dx / (W_IN * 72), anchor[1] - dy / (H_IN * 72)]
o1, o2, kl, kr = key_geom(raw); c_new = (kl + kr) / 2.0
dx1, dx2 = o1[0] - KEY1_ORG[0], o2[0] - KEY2_ORG[0]
DRAWN["key_placement"] = {"anchor_fig_frac": anchor, "loc": "lower center on the V13 key's own ink centre", "T90 >10% origin": o1, "base": KEY1_ORG,
                          "Laboratory sleep under 5 h origin": o2, "base2": KEY2_ORG, "dx_key1": dx1, "dx_key2": dx2, "v13_box": [KL_13, KR_13], "v13_centre": K13_C,
                          "new_box": [kl, kr], "new_centre": c_new, "growth_pt": (kr - kl) - (KR_13 - KL_13), "axes_centre": AX_X0 * 72 + FOREST_W * 36}
assert abs(o1[1] - KEY1_ORG[1]) < 0.3 and abs(o2[1] - KEY2_ORG[1]) < 0.3, (o1, o2, KEY1_ORG, KEY2_ORG)
assert abs(c_new - K13_C) < 0.3 and abs((kl - KL_13) + (kr - KR_13)) < 0.3, ("the key does not grow symmetrically about the V13 key centre", KL_13, KR_13, kl, kr)

final = f"{WORK}/{SHEET}_built.pdf"                      # round 49: the strip is stamped afterwards into HERE/ED_Fig04.pdf
DRAWN["font_metrics_aligned"] = C.align_font_metrics(raw, final, BASE)
# drawn records (v14lib DRAWN_RECORD) from the final sheet's own text layer: the ten row labels (source = the results file's
# outcome key), the key and axis strings and the tick labels are the V13 strings (static). The sheet prints no estimate.
d = fitz.open(final); sp = C.spans_of(d[0]); d.close()
for s in sp:
    t = s["text"].strip()
    if t in TOPC:
        DRAWN["records"].append(dict(text=t, x=s["origin"][0], baseline=s["origin"][1], ha="left", size=s["size"], panel="single",
                                     source_file=SRC, source_key=f"vs_sleep_duration/outcomes/{t}", source_value=t, rule="keyname",
                                     note=f"row {TOPC.index(t) + 1} of the ten largest T90-significant conditions (T90 HR {VSO[t]['t90']['hr']:.3f})"))
    else:   # round 49: the V13 strings at +1 pt (tick labels centred on their ticks, the axis title and the key centred on the axes) recorded for the printed-values CSV
        DRAWN["records"].append(dict(text=t, x=s["origin"][0], baseline=s["origin"][1], ha="left", size=s["size"], panel="single",
                                     source_file="static:V13", source_key="V13 string at +1 pt, grown about its V13 anchor", source_value=t, rule="text"))
json.dump(DRAWN, open(f"{WORK}/{SHEET}_drawn.json", "w"), indent=1, default=float)
print(f"{SHEET}: rows {TOPC}")
print(f"T90 larger in {N_OX_LARGER} of {N_D}; key at {o1} (base {KEY1_ORG}); wrote {final}")
