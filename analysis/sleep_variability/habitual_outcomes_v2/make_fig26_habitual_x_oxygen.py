"""
Figure26_habitual_sleep_x_oxygen — habitual home sleep crossed with laboratory oxygenation.

The question the sheet answers is which of the two sleep measurements carries the risk. Habitual
home sleep is what the patient says they sleep, pulled out of clinical notes for 5,295 of the
19,173 patients. Oxygenation is measured in the laboratory on the night of the study. Crossing
them puts every patient in one of twelve cells and asks what moves.

  A  the 4 by 3 cell grid, N and events, so the sample behind every estimate is on the sheet
  B  the main result, six outcomes across four cells
  C  oxygen inside each habitual band, the same effect at every duration
  D  the check, habitual against laboratory sleep against laboratory T90 in the same
     5,295 patients, with the negative-control band drawn

Nothing here is refitted. Every number is read out of the finished result files in
Sleep_Variability_2026-08/habitual_outcomes and asserted against the value it is expected to
carry before it is drawn.

Colour carries one idea. Grey is normal laboratory oxygen, dark blue is low laboratory oxygen.
The two are separated in lightness as well as hue and each carries its own marker shape, so the
sheet survives a greyscale print. There are no individual data dots, no printed P values, and
the intervals are drawn as error bars.

Sources
  cells_crossed.csv          panel A, cell N and events, design cells12
  results_crossed.csv        panel B, model cells12 set primary
                             panel C, model o2_within_hab set primary
  summary_crossed.json       reference cell, within-band counts and medians, sub-5-event counts
  attack_power_cells.csv     minimum detectable HR per cell, control behaviour in the weak cell
  attack_headtohead.csv      panel D, per-SD head to head, 45 non-control non-circular outcomes
  attack_log.txt             the published head-to-head table this panel reproduces
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

sys.path.insert(0, f"{paths.FIGURE_ROOT}/New_Figures/"
                   "_NOT_THESE_workfiles/scripts")
from mainstyle import (rc, panel_label, sep_ok, stars, logticks, INK, BLUE_D, GREY_M,  # noqa
                       BRONZE, save3, NEWFIG)

SRC = Path(f"{paths.SV_ROOT}/"
           "habitual_outcomes_v2")
WORK = Path(f"{paths.FIGURE_ROOT}/New_Figures/"
            "_NOT_THESE_workfiles")

# ------------------------------------------------------------------ type, all at or above 9 pt
FLOOR = 9.0
TITLE, LABEL, TICK, ANNOT, TINY = 12.0, 10.5, 9.5, 9.5, 9.0
PANEL = 14.0

# ------------------------------------------------------------------------------------- colour
NORMAL, LOW = GREY_M, BLUE_D          # the one contrast the sheet is built on
MID = "#97acb8"                       # the intermediate oxygen column of panel A only
PALE = "#dfe3e6"                      # the normal-oxygen fill of panel A
RULE = "#c8ccd0"
sep_ok("normal against low oxygen", [NORMAL, LOW])
sep_ok("panel A oxygen ramp", [PALE, LOW])

MK_NORMAL, MK_LOW = "s", "o"          # shape carries the contrast when colour cannot

# ------------------------------------------------------------------------------------ reading
cells = pd.read_csv(SRC / "cells_crossed.csv")
crossed = pd.read_csv(SRC / "results_crossed.csv")
power = pd.read_csv(SRC / "attack_power_cells.csv")
h2h = pd.read_csv(SRC / "attack_headtohead.csv")
S = json.load(open(SRC / "summary_crossed.json"))

assert S["cohort_n"] == 19173 and S["habitual_n"] == 5295   # v2 (v1: 8711)

# Panel A draws the 8,711 with a habitual sleep report and panel D said "the same 8,709
# patients", so one sheet gave two sizes for one set. The head-to-head of ATTACK 6 is fitted on
# HABSET, which is that same 8,711, and 8,711 is the largest n in attack_headtohead.csv. The
# 8,709 is the habitual arm of attack_selection_smd.csv, which is the selection comparison
# against the 10,458 without a report, not the modelled set. It is read from the result file
# now rather than typed.
HEAD_TO_HEAD_N = int(h2h.n.max())
assert HEAD_TO_HEAD_N == S["habitual_n"] == 5295, HEAD_TO_HEAD_N   # v2 (v1: 8711)
assert S["reference_cell"] == "6-<7h(ref) | normal T90<=1%" and S["reference_cell_n"] == 402

# The sheet referred to the negative controls in four places and named none of them, so a
# reader could not tell which panel was drawn. The names are read from the settled record and
# checked against the fitted file, never typed in here.
NEGC = json.load(open(f"{paths.NUMBERS_DIR}/"
                      "negcontrols_final.json"))
NEG_LABELS = list(NEGC["panel_labels"])
assert sorted(h2h[h2h.control].disease.unique()) == sorted(NEG_LABELS), (
    "Figure26 is drawing a different negative-control panel from negcontrols_final.json")
assert "Fracture" not in NEG_LABELS and "Osteoarthritis" not in NEG_LABELS
NEG_SENTENCE = ("The settled negative-control panel of "
                + NEGC["decision_date"] + ", "
                + str(NEGC["n_controls"]) + " conditions: "
                + ", ".join(NEG_LABELS[:-1]).lower() + " and " + NEG_LABELS[-1].lower()
                + ". Fracture and osteoarthritis are out of the panel.")

HAB = ["<5h", "5-<6h", "6-<7h(ref)", ">=7h"]
HAB_LAB = ["Under 5 h", "5 to under 6 h", "6 to under 7 h", "7 h or more"]
O2 = ["normal T90<=1%", "intermediate 1-10%", "low T90>10%"]
O2_LAB = ["Normal oxygen\nT90 1% or less", "Intermediate\nT90 1 to 10%", "Low oxygen\nT90 over 10%"]
O2_FILL = [PALE, MID, LOW]
O2_TEXT = [INK, INK, "white"]

OUT = [("cvd", "Cardiovascular composite"), ("death", "Death from any cause"),
       ("htn2", "Hypertension"), ("diabetes", "Type 2 diabetes"),
       ("hf", "Heart failure"), ("dementia", "Dementia")]
EVCOL = [("ev_cvd", "CVD"), ("ev_death", "Death"), ("ev_htn2", "HTN"),
         ("ev_diabetes", "T2D"), ("ev_hf", "HF")]

REF_CELL = "6-<7h(ref) | normal T90<=1%"
WEAK_CELL = "6-<7h(ref) | low T90>10%"
WEAK_CELL_POWER = "6-<7h(ref) x low T90>10%"

DRAWN = {"figure": "Figure26_habitual_sleep_x_oxygen",
         "built_from": str(SRC),
         "cohort": {"cohort_n": S["cohort_n"], "habitual_n": S["habitual_n"],
                    "habitual_pct": round(100 * S["habitual_n"] / S["cohort_n"], 1),
                    "head_to_head_n": HEAD_TO_HEAD_N,
                    "source": "summary_crossed.json habitual_n, which is also the largest n "
                              "in attack_headtohead.csv. The 8,709 this sheet used to print "
                              "for panel D is the habitual arm of attack_selection_smd.csv, "
                              "a different quantity."},
         "panelA": {}, "panelB": {}, "panelC": {}, "panelD": {}, "caveats": {}}

# ---------------------------------------------------------------------------- panel A, the grid
A = cells[cells.design == "cells12"].set_index("cell")
assert len(A) == 12 and int(A.n.sum()) == 5295

gridA = {}
for hb in HAB:
    for ob in O2:
        key = f"{hb} | {ob}"
        row = A.loc[key]
        gridA[key] = {"n": int(row.n),
                      "events": {lab: int(row[c]) for c, lab in EVCOL},
                      "median_habitual_h": float(row.median_habitual_h),
                      "median_t90_pct": float(row.median_t90_pct),
                      "is_reference": bool(row.is_reference)}
assert gridA["<5h | normal T90<=1%"]["n"] == 662
assert list(gridA["<5h | normal T90<=1%"]["events"].values()) == [65, 19, 71, 57, 41]
assert gridA["<5h | low T90>10%"]["n"] == 257
assert list(gridA["<5h | low T90>10%"]["events"].values()) == [40, 18, 23, 34, 30]
assert gridA["5-<6h | normal T90<=1%"]["n"] == 301 and gridA["5-<6h | low T90>10%"]["n"] == 87
assert gridA[REF_CELL]["n"] == 402 and gridA[WEAK_CELL]["n"] == 96
assert gridA[">=7h | normal T90<=1%"]["n"] == 1622
assert list(gridA[">=7h | normal T90<=1%"]["events"].values()) == [121, 86, 145, 68, 61]
assert gridA[">=7h | low T90>10%"]["n"] == 329
assert list(gridA[">=7h | low T90>10%"]["events"].values()) == [49, 32, 32, 35, 39]
DRAWN["panelA"] = {"cells": gridA,
                   "reference_cell": REF_CELL,
                   "unreliable_cell": WEAK_CELL,
                   "source": "cells_crossed.csv, design == cells12"}

# the weak cell, measured rather than asserted
pw = power[power.cell == WEAK_CELL_POWER]
weak_mdhr = float(pw.mdhr80.median())
weak_ctrl = pw[pw.control]
weak = {"n": gridA[WEAK_CELL]["n"],
        "median_mdhr80": round(weak_mdhr, 2),
        "outcomes_under_5_events": int(S["cells_with_outcomes_under_5_events"][WEAK_CELL]),
        "fittable_negative_controls": int(weak_ctrl.hr.notna().sum()),
        "significant_negative_controls": int((weak_ctrl.p < .05).sum()),
        "significant_control_names": sorted(weak_ctrl.loc[weak_ctrl.p < .05, "disease"]),
        "source": "attack_power_cells.csv (mdhr80, controls); summary_crossed.json (<5 events)"}
assert weak["outcomes_under_5_events"] == 10 and weak["significant_negative_controls"] == 1
DRAWN["caveats"]["weak_cell"] = weak

# --------------------------------------------------------------------- panel B, the main result
CELLS4 = [("<5h", "normal T90<=1%"), ("<5h", "low T90>10%"),
          (">=7h", "normal T90<=1%"), (">=7h", "low T90>10%")]
CELLS4_LAB = ["Under 5 h, normal oxygen", "Under 5 h, low oxygen",
              "7 h or more, normal oxygen", "7 h or more, low oxygen"]

cb = crossed[(crossed.model == "cells12") & (crossed["set"] == "primary")]
panelB = {}
for key, name in OUT:
    d = cb[cb.outcome_key == key].set_index("contrast")
    rows = []
    for (hb, ob), lab in zip(CELLS4, CELLS4_LAB):
        r = d.loc[f"{hb} | {ob}"]
        rows.append({"cell": f"{hb} | {ob}", "label": lab, "n": int(r.n_group),
                     "events": int(r.events_group), "hr": float(r.hr),
                     "lo": float(r.lo95), "hi": float(r.hi95), "p": float(r.p),
                     "oxygen": "normal" if ob.startswith("normal") else "low"})
    panelB[name] = rows
# the four verified cardiovascular values, to two decimals as the brief carries them
_cvd = [round(r["hr"], 2) for r in panelB["Cardiovascular composite"]]
assert _cvd == [1.41, 2.48, 1.11, 2.52], _cvd
assert [round(r["hr"], 2) for r in panelB["Heart failure"]] == [1.71, 3.79, 1.13, 3.79]
assert [round(r["hr"], 2) for r in panelB["Death from any cause"]] == [0.65, 1.64, 1.17, 2.24]
assert [round(r["hr"], 2) for r in panelB["Hypertension"]] == [1.36, 1.83, 0.89, 1.64]
assert [round(r["hr"], 2) for r in panelB["Type 2 diabetes"]] == [1.78, 3.11, 0.89, 2.33]
assert [round(r["hr"], 2) for r in panelB["Dementia"]] == [1.13, 0.66, 1.41, 1.70]
DRAWN["panelB"] = {"reference": f"{REF_CELL}, N=402",
                   "estimates": panelB,
                   "source": "results_crossed.csv, model == cells12, set == primary"}

# ------------------------------------------------------- panel C, oxygen inside each habitual band
co = crossed[(crossed.model == "o2_within_hab") & (crossed["set"] == "primary")]
panelC, bandsum = {}, {}
for hb, hl in zip(HAB, HAB_LAB):
    d = co[co.hab_band == hb].set_index("outcome_key")
    rows = []
    for key, name in OUT:
        r = d.loc[key]
        rows.append({"outcome": name, "n_low": int(r.n_group), "events_low": int(r.events_group),
                     "hr": float(r.hr), "lo": float(r.lo95), "hi": float(r.hi95),
                     "p": float(r.p)})
    panelC[hb] = rows
    s = S["o2_within_band"][hb]
    bandsum[hb] = {"label": hl, "fitted": int(s["outcomes_fitted"]),
                   "significant": int(s["significant"]), "median_hr": round(s["median_hr"], 2)}
assert [bandsum[b]["significant"] for b in HAB] == [18, 9, 15, 26]
assert [bandsum[b]["fitted"] for b in HAB] == [41, 26, 32, 44]
assert [bandsum[b]["median_hr"] for b in HAB] == [1.55, 1.71, 2.10, 1.80]
DRAWN["panelC"] = {"contrast": "low T90 over 10% against normal T90 1% or less, inside the band",
                   "estimates": panelC, "band_summary": bandsum,
                   "source": "results_crossed.csv model == o2_within_hab set == primary; "
                             "summary_crossed.json o2_within_band"}

# the mirror test, habitual band against the reference at fixed oxygen
mirror = {ob: {"fitted": int(S["hab_within_o2"][ob]["outcomes_fitted"]),
               "significant": int(S["hab_within_o2"][ob]["significant"]),
               "median_hr": round(S["hab_within_o2"][ob]["median_hr"], 2)} for ob in O2}
assert mirror["normal T90<=1%"] == {"fitted": 132, "significant": 6, "median_hr": 1.08}
assert mirror["intermediate 1-10%"]["significant"] == 7
assert mirror["low T90>10%"]["significant"] == 7
DRAWN["caveats"]["mirror_test"] = {"values": mirror,
                                   "source": "summary_crossed.json hab_within_o2"}
DRAWN["caveats"]["across_all_outcomes"] = {
    "low_oxygen_cells_raised": int(S["low_o2_cells_significant"]),
    "low_oxygen_cells_fitted": int(S["low_o2_cells_fitted"]),
    "low_oxygen_median_hr": round(S["median_hr_low_o2_cells"], 2),
    "normal_oxygen_cells_raised": int(S["normal_o2_cells_significant"]),
    "normal_oxygen_cells_fitted": int(S["normal_o2_cells_fitted"]),
    "normal_oxygen_median_hr": round(S["median_hr_normal_o2_cells"], 2),
    "source": "summary_crossed.json"}
assert DRAWN["caveats"]["across_all_outcomes"]["low_oxygen_cells_raised"] == 71
assert DRAWN["caveats"]["across_all_outcomes"]["low_oxygen_median_hr"] == 1.85
assert DRAWN["caveats"]["across_all_outcomes"]["normal_oxygen_median_hr"] == 1.09

# --------------------------------------------------------------- panel D, the measurement check
EXPO = [("z_hab", "Habitual home sleep\nself-reported, per SD", NORMAL, MK_NORMAL),
        ("z_tst", "Laboratory total sleep time\nmeasured, per SD", NORMAL, "^"),
        ("z_t90", "Laboratory T90\nmeasured, per SD", LOW, MK_LOW)]
panelD = {}
for key, lab, _, _ in EXPO:
    x = h2h[h2h.exposure == key]
    ctrl = x[x.control]
    main = x[(~x.control) & (~x.circular)]
    panelD[key] = {"label": lab.replace("\n", ", "),
                   "n_outcomes": int(len(main)),
                   "significant": int((main.p < .05).sum()),
                   "surviving_bh": int((main.q < .05).sum()),
                   "control_lo": round(float(ctrl.hr.min()), 3),
                   "control_hi": round(float(ctrl.hr.max()), 3),
                   "control_floor_abs": round(float(ctrl.abs_hr.max()), 3),
                   "median_abs_hr": round(float(main.abs_hr.median()), 3),
                   "max_abs_hr": round(float(main.abs_hr.max()), 3)}
assert panelD["z_hab"]["significant"] == 10 and panelD["z_hab"]["surviving_bh"] == 5
assert panelD["z_tst"]["significant"] == 1 and panelD["z_tst"]["surviving_bh"] == 0
assert panelD["z_t90"]["significant"] == 33 and panelD["z_t90"]["surviving_bh"] == 32
assert all(v["n_outcomes"] == 45 for v in panelD.values())
assert panelD["z_hab"]["control_lo"] == 0.864 and panelD["z_hab"]["control_hi"] == 0.930
DRAWN["panelD"] = {"estimates": panelD, "n_patients": HEAD_TO_HEAD_N,
                   "source": "attack_headtohead.csv, control == False and circular == False; "
                             "reproduces attack_log.txt ATTACK 6"}

hb_ = h2h[h2h.exposure == "z_hab"]
hb_ctrl_min = float(hb_[hb_.control].hr.min())
hb_main = hb_[(~hb_.control) & (~hb_.circular)]
outside = hb_main[hb_main.hr < hb_ctrl_min]
DRAWN["caveats"]["confounding_band"] = {
    "negative_control_range": [round(hb_ctrl_min, 3),
                               round(float(hb_[hb_.control].hr.max()), 3)],
    "outcomes_outside_band": int(len(outside)),
    "outcomes_outside_band_after_bh": int((outside.q < .05).sum()),
    "names": sorted(outside.disease),
    "source": "attack_headtohead.csv, exposure == z_hab"}
assert DRAWN["caveats"]["confounding_band"]["outcomes_outside_band"] == 4
assert DRAWN["caveats"]["confounding_band"]["outcomes_outside_band_after_bh"] == 3
# ================================================================================ drawing
#
# Every axes is placed by hand in inches on the sheet rather than by a grid. A gridspec puts
# panels in the right boxes and then lets the labels that hang off the side of a panel run into
# the panel next door, which is exactly what has to not happen here. Placing the boxes and the
# text gutters explicitly means the sheet can be reasoned about with a ruler.
plt.rcParams.update(rc())
plt.rcParams.update({"font.size": TICK, "axes.labelsize": LABEL, "xtick.labelsize": TICK,
                     "ytick.labelsize": TICK, "legend.fontsize": ANNOT,
                     "axes.titlesize": TITLE})

W, H = 17.0, 22.5
fig = plt.figure(figsize=(W, H))


def AX(x0, y0, w, h):
    return fig.add_axes([x0 / W, y0 / H, w / W, h / H])


def FT(x, y, s, size=TINY, weight=None, style=None, va="top", ha="left", color=INK,
       ls=1.7, **kw):
    return fig.text(x / W, y / H, s, fontsize=size, fontweight=weight, style=style,
                    va=va, ha=ha, color=color, linespacing=ls, **kw)


def plab(x, y, letter):
    fig.text(x / W, y / H, letter, fontsize=PANEL, fontweight="bold", va="bottom", ha="left",
             color=INK)


def fmt(v):
    return f"{v:,}"


# ------------------------------------------------------------------------------------ panel A
plab(0.35, 21.52, "A")
FT(0.80, 21.52, "Where the 5,295 patients sit when habitual home sleep is crossed with "
                "laboratory oxygenation,\nand the events behind every cell",
   size=TITLE, va="bottom")

axA = AX(2.35, 17.95, 7.65, 3.35)
axA.set_xlim(0, 3)
axA.set_ylim(0, 4.55)
axA.axis("off")

for j, ol in enumerate(O2_LAB):
    axA.text(j + 0.5, 4.06, ol, ha="center", va="bottom", fontsize=ANNOT, color=INK,
             linespacing=1.4)
for i, hl in enumerate(HAB_LAB):
    axA.text(-0.035, 3.5 - i, hl, ha="right", va="center", fontsize=ANNOT, color=INK)
fig.text(0.55 / W, (17.95 + 2.94 / 2) / H, "Habitual home sleep", rotation=90, va="center",
         ha="center", fontsize=ANNOT, color=INK, style="italic")

for i, hb in enumerate(HAB):
    for j, ob in enumerate(O2):
        key = f"{hb} | {ob}"
        c = gridA[key]
        x0, y0, w, hgt = j + 0.03, (3 - i) + 0.06, 0.94, 0.88
        axA.add_patch(Rectangle((x0, y0), w, hgt, facecolor=O2_FILL[j], edgecolor="white",
                                lw=1.2, zorder=1))
        tc, top, mid = O2_TEXT[j], y0 + hgt, x0 + w / 2
        flag = key in (REF_CELL, WEAK_CELL)
        if key == REF_CELL:
            axA.add_patch(Rectangle((x0, y0), w, hgt, facecolor="none", edgecolor=INK,
                                    lw=2.2, zorder=4))
            axA.text(mid, top - 0.135, "REFERENCE", ha="center", va="center", fontsize=TINY,
                     color=INK, fontweight="bold", zorder=5)
        if key == WEAK_CELL:
            axA.add_patch(Rectangle((x0, y0), w, hgt, facecolor="none", edgecolor=BRONZE,
                                    lw=2.2, ls=(0, (3.2, 2.0)), zorder=4))
            axA.text(mid, top - 0.135, "UNRELIABLE", ha="center", va="center", fontsize=TINY,
                     color="#e2d6c4", fontweight="bold", zorder=5)
        ny = top - (0.35 if flag else 0.25)
        axA.text(mid, ny, f"N = {fmt(c['n'])}", ha="center", va="center", fontsize=TITLE,
                 color=tc, fontweight="bold", zorder=5)
        ev = c["events"]
        axA.text(mid, ny - 0.22,
                 f"CVD {ev['CVD']}     Death {ev['Death']}     HTN {ev['HTN']}",
                 ha="center", va="center", fontsize=TINY, color=tc, zorder=5)
        axA.text(mid, ny - 0.40, f"T2D {ev['T2D']}     HF {ev['HF']}", ha="center",
                 va="center", fontsize=TINY, color=tc, zorder=5)

FT(0.35, 17.62,
   "Events are incident cases of the cardiovascular composite, death from any cause, "
   "hypertension, type 2 diabetes and heart failure.\n"
   "Reference cell 6 to under 7 h with T90 at or below 1%, N = 402, outlined in black. The 6 to "
   f"under 7 h low-oxygen cell, N = {weak['n']}, is marked unreliable:\n"
   f"median minimum detectable hazard ratio {weak['median_mdhr80']:.2f}, "
   f"{weak['outcomes_under_5_events']} outcomes with fewer than 5 events, and "
   f"{weak['significant_negative_controls']} of its {weak['fittable_negative_controls']} "
   "fittable negative controls are themselves raised.")

# ------------------------------------------------------------------------------------ panel D
plab(10.55, 21.52, "D")
FT(11.00, 21.52, f"The same {HEAD_TO_HEAD_N:,} patients, three measurements side by side\n"
                 f"against the band the negative controls occupy",
   size=TITLE, va="bottom")

axD = AX(11.45, 18.45, 5.05, 2.60)
DX = (1.0, 3.85)
ypos = [2.4, 1.3, 0.2]
axD.set_xscale("log")
axD.set_xlim(*DX)
axD.set_ylim(-0.55, 3.05)
axD.set_xticks([1.0, 1.2, 1.5, 2.0, 3.0])
axD.set_xticklabels(["1.0", "1.2", "1.5", "2.0", "3.0"])
axD.xaxis.set_minor_formatter(plt.NullFormatter())
axD.set_xlabel("Hazard ratio per standard deviation, folded to the raised direction",
               fontsize=LABEL, labelpad=3)
axD.set_yticks([])
for s in ("top", "right", "left"):
    axD.spines[s].set_visible(False)

for (key, lab, col, mk), y in zip(EXPO, ypos):
    v = panelD[key]
    axD.add_patch(Rectangle((1.0, y - 0.235), v["control_floor_abs"] - 1.0, 0.47,
                            facecolor="#e6eaec", edgecolor=RULE, lw=0.8, zorder=1))
    axD.plot([v["median_abs_hr"], v["max_abs_hr"]], [y, y], color=col, lw=2.0,
             solid_capstyle="round", zorder=3)
    axD.scatter([v["median_abs_hr"]], [y], s=72, marker=mk, facecolors="white",
                edgecolors=col, linewidths=1.8, zorder=4)
    axD.scatter([v["max_abs_hr"]], [y], s=72, marker=mk, facecolors=col,
                edgecolors="white", linewidths=0.9, zorder=4)
    axD.text(1.008, y + 0.30, lab.replace("\n", ", "), fontsize=ANNOT, color=INK,
             va="bottom", ha="left")
    axD.text(1.008, y - 0.30,
             f"{v['significant']} of {v['n_outcomes']} outcomes raised, "
             f"{v['surviving_bh']} after Benjamini-Hochberg."
             f"  Negative controls {v['control_lo']:.2f} to {v['control_hi']:.2f}",
             fontsize=TINY, color=INK, va="top", ha="left")
axD.axvline(1.0, color=INK, lw=1.0, zorder=2)

FT(10.55, 17.62,
   "Open marker, the median across the 45 outcomes. Filled marker, the largest single outcome.\n"
   "Shaded band, the largest effect any negative control reaches in the same model, which is\n"
   "the size an effect has to clear before it means anything. The five are named in note 5.")

# ------------------------------------------------------------------------------------ panel B
plab(0.35, 16.60, "B")
FT(0.80, 16.60, "Risk against the reference cell. At both sleep durations the normal-oxygen "
                "cells sit on the no-difference line\nand the low-oxygen cells sit well to the "
                "right of it", size=TITLE, va="bottom")

axB = AX(5.00, 10.25, 9.00, 6.00)
BX = (0.45, 6.2)
yy, gaps, cur = [], [], 0.0
for name in [n for _, n in OUT]:
    cur -= 1.15
    gaps.append((cur, name))
    for r in panelB[name]:
        cur -= 1.0
        yy.append((cur, r))
cur -= 0.60
axB.set_xscale("log")
axB.set_xlim(*BX)
axB.set_ylim(cur, 0.30)
axB.set_yticks([])
axB.set_xticks([0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0])
axB.set_xticklabels(["0.5", "0.75", "1.0", "1.5", "2.0", "3.0", "5.0"])
axB.xaxis.set_minor_formatter(plt.NullFormatter())
axB.set_xlabel("Hazard ratio against 6 to under 7 h with normal oxygen (N = 402)",
               fontsize=LABEL, labelpad=3)
for s in ("top", "right", "left"):
    axB.spines[s].set_visible(False)
axB.axvline(1.0, color=INK, lw=1.1, zorder=2)

TR = axB.get_yaxis_transform()
X_LAB, X_HEAD, X_NUM = (0.78 - 5.00) / 9.00, (0.35 - 5.00) / 9.00, (14.18 - 5.00) / 9.00
for y, name in gaps:
    axB.text(X_HEAD, y + 0.14, name, transform=TR, fontsize=LABEL, color=INK,
             fontweight="bold", va="center", ha="left", clip_on=False)
    axB.plot([X_HEAD, 1.0], [y + 0.60, y + 0.60], transform=TR, color=RULE, lw=0.7,
             zorder=1, clip_on=False)
for y, r in yy:
    col = NORMAL if r["oxygen"] == "normal" else LOW
    mk = MK_NORMAL if r["oxygen"] == "normal" else MK_LOW
    axB.plot([r["lo"], r["hi"]], [y, y], color=col, lw=2.0, solid_capstyle="round", zorder=3)
    axB.scatter([r["hr"]], [y], s=62, marker=mk, color=col, edgecolors="white",
                linewidths=0.8, zorder=4)
    axB.text(X_LAB, y, f"{r['label']}  (N = {fmt(r['n'])}, {r['events']} events)",
             transform=TR, fontsize=ANNOT, color=INK, va="center", ha="left", clip_on=False)
    axB.text(X_NUM, y, f"{r['hr']:.2f} (95% CI, {r['lo']:.2f}-{r['hi']:.2f}) {stars(r['p'])}",
             transform=TR, fontsize=ANNOT, color=INK, va="center", ha="left", clip_on=False)
axB.text(X_NUM, 0.10, "Hazard ratio (95% CI)", transform=TR, fontsize=ANNOT, color=INK,
         va="center", ha="left", fontweight="bold", clip_on=False)

fig.legend(handles=[Line2D([], [], marker=MK_NORMAL, color=NORMAL, lw=2.2,
                           markeredgecolor="white", markersize=8,
                           label="Normal laboratory oxygen, T90 1% or less   "
                                 "(N = 662 at under 5 h and 1,622 at 7 h or more)"),
                    Line2D([], [], marker=MK_LOW, color=LOW, lw=2.2,
                           markeredgecolor="white", markersize=8,
                           label="Low laboratory oxygen, T90 over 10%   "
                                 "(N = 257 at under 5 h and 329 at 7 h or more)")],
           loc="upper left", frameon=False, fontsize=ANNOT, handletextpad=0.7,
           labelspacing=0.55, bbox_to_anchor=(0.80 / W, 9.55 / H), bbox_transform=fig.transFigure)
FT(0.80, 8.95, "Asterisks mark * P below .05, ** P below .01, *** P below .001. The two "
               "intermediate habitual bands and the intermediate oxygen column are in panels "
               "A and C.")

# ------------------------------------------------------------------------------------ panel C
plab(0.35, 8.30, "C")
FT(0.80, 8.30, "Low against normal laboratory oxygen inside each habitual band. The oxygen "
               "effect is the same size at every habitual sleep duration,\nso it is not a "
               "duration effect wearing an oxygen label", size=TITLE, va="bottom")

CX = (0.09, 11.0)
CW, CGAP, CX0, CY0, CHH = 3.25, 0.30, 2.75, 3.35, 4.50
for k, hb in enumerate(HAB):
    x0 = CX0 + k * (CW + CGAP)
    ax = AX(x0, CY0, CW, CHH)
    ax.set_xscale("log")
    ax.set_xlim(*CX)
    ax.set_ylim(-0.80, len(OUT) - 0.20)
    ax.set_xticks([0.25, 1.0, 4.0])
    ax.set_xticklabels(["0.25", "1.0", "4.0"])
    ax.xaxis.set_minor_formatter(plt.NullFormatter())
    ax.set_yticks([])
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.axvline(1.0, color=NORMAL, lw=2.2, zorder=2)
    ax.set_xlabel("Hazard ratio", fontsize=TICK, labelpad=2)
    ax.set_title(bandsum[hb]["label"], fontsize=LABEL, color=INK, pad=7)

    for i, r in enumerate(panelC[hb]):
        y = len(OUT) - 1 - i
        ax.plot([r["lo"], r["hi"]], [y, y], color=LOW, lw=2.0, solid_capstyle="round", zorder=3)
        ax.scatter([r["hr"]], [y], s=58, marker=MK_LOW, color=LOW, edgecolors="white",
                   linewidths=0.8, zorder=4)
        ax.text(0.012, y + 0.24,
                f"{r['hr']:.2f} (95% CI, {r['lo']:.2f}-{r['hi']:.2f}) {stars(r['p'])}",
                transform=ax.get_yaxis_transform(), fontsize=TINY, color=INK, va="bottom",
                ha="left")
        if k == 0:
            ax.text((0.35 - CX0) / CW, y, r["outcome"], transform=ax.get_yaxis_transform(),
                    fontsize=ANNOT, color=INK, va="center", ha="left", clip_on=False)
    b = bandsum[hb]
    FT(x0 + CW / 2, 2.75,
       f"{b['significant']} of {b['fitted']} outcomes raised\nmedian HR "
       f"{b['median_hr']:.2f}", size=ANNOT, ha="center", ls=1.5)

fig.legend(handles=[Line2D([], [], marker=MK_LOW, color=LOW, lw=2.2, markeredgecolor="white",
                           markersize=8,
                           label="Low oxygen against normal oxygen within the same habitual "
                                 "band   (N = 257, 87, 96 and 329 across the four bands)"),
                    Line2D([], [], color=NORMAL, lw=2.4,
                           label="Normal oxygen within the same band, the reference   "
                                 "(N = 662, 301, 402 and 1,622)")],
           loc="upper left", frameon=False, fontsize=ANNOT, handletextpad=0.7,
           labelspacing=0.55, bbox_to_anchor=(0.80 / W, 2.15 / H), bbox_transform=fig.transFigure)

# --------------------------------------------------------------------------------- footnotes
cb_ = DRAWN["caveats"]["confounding_band"]
mt = DRAWN["caveats"]["mirror_test"]["values"]
aa = DRAWN["caveats"]["across_all_outcomes"]
FT(0.80, 1.45,
   "1  The whole habitual-sleep effect sits inside the confounding band. All five negative "
   f"controls run against habitual sleep at {cb_['negative_control_range'][0]:.2f} to "
   f"{cb_['negative_control_range'][1]:.2f} per standard deviation, the same direction and "
   f"about the same size as the real outcomes. Only {cb_['outcomes_outside_band']} of 45\n"
   f"    outcomes fall outside that band and {cb_['outcomes_outside_band_after_bh']} survive "
   "Benjamini-Hochberg correction (panel D).\n"
   "2  The 6 to under 7 h low-oxygen cell is marked unreliable in panel A and its estimates "
   "should not be read on their own.\n"
   "3  Habitual sleep is self-reported and was extracted from clinical notes for "
   f"{S['habitual_n']:,} of {S['cohort_n']:,} patients "
   f"({100 * S['habitual_n'] / S['cohort_n']:.0f}%), skewed towards one hospital, so panel A "
   "is not a random sample of the cohort.\n"
   f"4  Moving the habitual band while oxygen is held fixed does almost nothing: "
   f"{mt['normal T90<=1%']['significant']} of {mt['normal T90<=1%']['fitted']} outcomes are "
   f"raised within normal oxygen (median HR {mt['normal T90<=1%']['median_hr']:.2f}), "
   f"{mt['intermediate 1-10%']['significant']} of {mt['intermediate 1-10%']['fitted']} within "
   f"intermediate and {mt['low T90>10%']['significant']} of\n"
   f"    {mt['low T90>10%']['fitted']} within low. Across all outcomes the three low-oxygen "
   f"cells are raised in {aa['low_oxygen_cells_raised']} of {aa['low_oxygen_cells_fitted']} "
   f"contrasts (median HR {aa['low_oxygen_median_hr']:.2f}) against "
   f"{aa['normal_oxygen_cells_raised']} of {aa['normal_oxygen_cells_fitted']} for the "
   f"normal-oxygen cells (median HR {aa['normal_oxygen_median_hr']:.2f}).\n"
   f"5  {NEG_SENTENCE} A negative control is a condition nocturnal oxygen cannot plausibly "
   "cause, so whatever it shows is confounding.")

# ------------------------------------------------------------------------------------- filing
RAS = str(SRC)   # v2: rasters filed with the v2 outputs, not the workfiles
# 2026-08-08. This sheet is cited in the manuscript as eFigure 16, not as a main figure, so it
# is filed with the supplementary set under that name. It used to write into the main folder
# under Figure26, which left the eFigure 16 callout pointing at a file that did not exist and
# put two different figures under the same Figure26 number.
name = "eFigure16_habitual_sleep_x_oxygen_v2"


def _supp_dir():
    import os
    for d in sorted(os.listdir(NEWFIG)):
        if not os.path.isdir(f"{NEWFIG}/{d}"):
            continue
        key = d.lower().lstrip("_0123456789")
        if key.startswith("supplementary_figures") or key.startswith("supplement"):
            return f"{NEWFIG}/{d}"
    raise SystemExit("no supplementary figures folder found")


pdf_path = f"{SRC}/{name}.pdf"   # v2: never into New_Figures, the main session gates that
fig.savefig(pdf_path)
fig.savefig(f"{RAS}/{name}.png", dpi=200)
fig.savefig(f"{RAS}/{name}.tif", dpi=300, pil_kwargs={"compression": "tiff_lzw"})
plt.close(fig)

json.dump(DRAWN, open(SRC / "Figure26_drawn_values_v2.json", "w"), indent=1)
print("wrote", pdf_path)
print("wrote", f"{RAS}/{name}.png")
print("wrote", SRC / "Figure26_drawn_values_v2.json")
