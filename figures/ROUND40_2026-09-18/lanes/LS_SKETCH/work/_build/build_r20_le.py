#!/usr/bin/env python3
"""Round 20, lane LB: the Figure 4 bottom sketch, the Figure 5 bottom sketch and Figure 6.

WHAT ROUND 20 CHANGES, and nothing else.

band_fig4_bottom  970.79 pt
  7  the title line is deleted
  8  the sleeping person is the SAME FACE as Figure 1's, closed eyes and a neutral mouth.
     Round 18 drew this one with "sleeper-bed", which smiles, so the two sheets differed.
  9  the outcome count is the truth for THIS analysis. The apnea-by-oxygen file
     numbers/alenfig3_cross_v1.json carries 14 outcomes, not 50: cardiovascular composite,
     heart failure, respiratory failure, pulmonary hypertension, COPD, type 2 diabetes,
     chronic kidney disease, acute kidney injury and death from any cause, plus the five
     negative controls (back pain, cataract, glaucoma, contact dermatitis, hemorrhoids).
     The heading therefore reads "14 diseases and death", in the same grammar Figures 1 and
     3 use for their own 50, and 50 is never printed on this sheet.

band_fig5_bottom  968.94 pt
 10  the title line is deleted
 11  the disease claim is replaced by the measured result. numbers/treatment_honest_v3.json:
     42 outcomes, 13 of them significantly lower once the oxygen was restored on PAP, median
     39.8 per cent lower among those 13 (16.9 per cent across all 42). The heading prints the
     exact 39.8, not a rounded 40, because the manuscript prints 39.8.

Main_Fig6         968.94 pt
 12  the title line is deleted
 13  the organ block sits on a grey background panel
 14  the two band labels are set in filled pills, each with a rule that lands on its own
     arrow, so neither floats any more
 15  the key is centred under the oxygen node and says what it means
 16  both sleeping people are the SAME FACE as Figure 1's
 17  both cross badges take their x from ONE constant, r20_common.ROAD_BADGE_X
 18  the Benjamini-Hochberg line is deleted from both roads. The correction is stated in the
     figure legend instead, so the artwork keeps only the short claim.

COLOUR, round 20: sleep apnea and breathing events are CHARCOAL #37474f on the ladder
#d3d7d8 #a5acb0 #737e84 #37474f with the dark text twin #29353b, replacing round 18's orange.

    python3 build_r20_le.py            # 150 dpi previews into _build/previews/
    python3 build_r20_le.py --final    # SVG + 300 dpi PNG + PDF into the lane root
"""
from __future__ import annotations

import argparse
import json

from l1_common import (APNEA, APNEA_D, AZ, AZ_D, AZ_T, CARD, GREEN_D, GREY, GREY_D, GREY_L, INK,
                       MM_PER_PT, ROOT, WHITE, airflow, finish, new, swatch_label,
                       text_width_mm)
from build_l1 import oxygen_node
from r18_family import HIGH, LOW, SCOPE, T90_HEADER, airflow_ramp, big_crossed_road, \
    organ_row, oxygen_fork
from r20_common import ROAD_BADGE_X, band_pill, severity_ladder
from roads_r18 import ROAD_APNEA, ROAD_DURATION

OUT = ROOT
CREDIT = "Icons: Servier Medical Art (CC BY 3.0)."
MASK = "cpap-mask-nostrap"                        # head + face + mask + hose, no strap
SLEEPER = "sleeper-faceB"                       # the SAME face as Figure 1 (round 20)

W4_PT, H4 = 970.79, 74.0                          # the Figure 4 page width
W5_PT, H5 = 968.94, 70.0                          # the Figure 5 page width, Fig 3's height
W6_PT, H6 = 968.94, 180.0                         # Figure 6 keeps the Figure 5 page width

# ---- the geometry the Figure 3 bottom sketch fixes, read out of its delivered PDF
PITCH = 28.5                                      # organ pitch, LD centres 205 .. 319
ORG_CX = 262.0                                    # centre of the organ column
HEAD_Y = 12.6                                     # the heading row, LD prints it at 10.7 top
TX, TW, TH = 110.0, 44.0, 15.0                    # the one oxygen trace pair
LX = TX + TW + 3.0

# the severity range at the two ends of the ramp, worded as Main_Fig4.pdf prints its groups
APNEA_NONE, APNEA_SEVERE = "No apnea (AHI <5)", "Severe (AHI 30 or more)"

# ---- what each sheet says over its organ column, and where the number comes from
SCOPE_F4 = "diseases and death"                   # R22: the 14 keys of alenfig3_cross_v1.json are 8 diseases + death + 5 NEGATIVE CONTROLS, and no printed table carries that 9-outcome set, so the band states no count
SCOPE_F5 = "lower risk of disease"   # R24: Alen wants only this phrase on Fig 5d (was "17 of 44 conditions lower, median 39.8% lower")   # numbers/treatment_v2.json corrected_vs_not (7 Aug), matches the manuscript text and the Fig 5b legend; the 4 Aug treatment_honest_v3.json summary (13 of 42) is superseded


def credit(fig):
    """Bottom right, where the Figure 3 bottom sketch puts it."""
    fig.label(CREDIT, fig.w - 4.0, fig.h - 5.2, pt=fig.theme.floor_pt, anchor="end",
              va="top", color=GREY)


def two_endings(fig, cy_top, cy_bot, *, from_x, scope=SCOPE):
    """The two endings, drawn as the Figure 3 bottom sketch draws them: the scope of the
    outcome once in dark grey over the column, coloured organs on the top branch, their dark
    twins on the bottom branch, and the five group names written ONCE, at the BOTTOM of the
    block under the dark-grey row."""
    fig.label(scope, ORG_CX, HEAD_Y, bold=True, anchor="middle", va="middle", color=GREY_D)
    for cy, dark in ((cy_top, False), (cy_bot, True)):
        fig.arrow((from_x, cy), (ORG_CX - 63.0, cy), color=GREY, lw=0.6)
        organ_row(fig, ORG_CX, cy, dark=dark, pitch=PITCH, names=dark)


def t90_header(fig, lab_r, *, text=T90_HEADER):
    return fig.label(text, (TX - 4.0 + lab_r) / 2, HEAD_Y, bold=True, color=AZ,
                     anchor="middle", va="middle")


# ================================================================ Figure 4, bottom sketch

def fig4_bottom(final, *, stem="band_fig4_bottom"):
    W = W4_PT * MM_PER_PT                          # 342.47 mm
    fig = new(W, H4)
    t = fig.theme

    # ---- the one source: every apnea severity, drawn once, in one breathing trace over a
    # wedge that widens AND darkens from no apnea to severe, through the charcoal ladder
    RX, RW = 4.0, 74.0                             # the ramp runs x 4 .. 78
    BX = RX + RW / 2                               # 41.0, the column's spine
    fig.place(SLEEPER, BX, 23.0, w=28.0)
    airflow_ramp(fig, RX, 38.0, RW, 8.5)           # this sleeper's breathing, all severities
    severity_ladder(fig, RX, 46.0, RW - 3.0, left_text=APNEA_NONE, right_text=APNEA_SEVERE,
                    tri_h=2.0)
    fig.label("sleep apnea severity", BX, 56.0, anchor="middle", va="top", bold=True,
              color=APNEA_D)

    # ---- the one oxygen fork, both branches out of that one source
    CY_TOP, CY_BOT = 23.5, 52.5                    # round-19 spacing, names below the grey
    oxygen_fork(fig, (84.0, 40.0), TX, CY_TOP - TH / 2, CY_BOT - TH / 2, TW, TH, label_x=LX,
                dip_centers=(0.16, 0.40, 0.64, 0.86))
    lab_r = LX + max(text_width_mm(LOW, t.label_pt, True),
                     text_width_mm(HIGH, t.label_pt, True))
    t90_header(fig, lab_r)

    two_endings(fig, CY_TOP, CY_BOT, from_x=lab_r + 2.5, scope=SCOPE_F4)
    credit(fig)
    return finish(fig, stem, OUT, final=final, attribution=None)


# ================================================================ Figure 5, bottom sketch

def fig5_bottom(final, *, stem="band_fig5_bottom"):
    W = W5_PT * MM_PER_PT                          # 341.82 mm
    fig = new(W, H5)
    t = fig.theme

    # ---- the one source: people whose sleep apnea is treated with PAP. NOTHING here says
    # or implies that every one of them desaturated before treatment, because the paper
    # itself shows that many of them did not. There is no before-PAP trace on this band.
    BX = 40.0
    fig.place(MASK, BX, 28.0, h=26.0, flip_h=True)      # ROUND 27: the head faces right
    fig.label("sleep apnea\ntreated with PAP", BX, 43.5, anchor="middle", va="top",
              bold=True, color=APNEA_D, leading=1.2)

    # ---- the one oxygen fork. On PAP the oxygen goes two ways and BOTH are drawn.
    CY_TOP, CY_BOT = 23.0, 48.5                    # round-19 spacing, names below the grey
    IMP, STILL = "oxygen improved", "oxygen still low"
    oxygen_fork(fig, (78.0, 34.0), TX, CY_TOP - TH / 2, CY_BOT - TH / 2, TW, TH, label_x=LX,
                low=IMP, high=STILL, dip_centers=(0.16, 0.40, 0.64, 0.86))
    lab_r = LX + max(text_width_mm(IMP, t.label_pt, True),
                     text_width_mm(STILL, t.label_pt, True))
    t90_header(fig, lab_r, text="time with oxygen below 90% on PAP")

    two_endings(fig, CY_TOP, CY_BOT, from_x=lab_r + 2.5, scope=SCOPE_F5)
    credit(fig)
    return finish(fig, stem, OUT, final=final, attribution=None)


# ================================================================ Figure 6, standalone

def fig6(final, *, stem="Main_Fig6"):
    W = W6_PT * MM_PER_PT                          # 341.82 mm
    fig = new(W, H6)
    t = fig.theme
    BX, BED_W = 42.0, 26.0
    Y_DUR, Y_PAP, Y_APNEA = 52.0, 95.0, 138.0
    AX = 70.0                                      # every input arrow starts here
    nx, ny, nw, nh = 112.0, 79.0, 58.0, 36.0       # the oxygen node
    ncy = ny + nh / 2                              # 97.0
    CY_TOP, CY_BOT = 56.0, 134.0
    ROAD_TOP, ROAD_BOT, ROAD_X = 17.0, 170.0, 336.0
    KEY_Y = 119.5                                  # the key, centred under the node
    PANEL = (195.5, 39.0, 132.5, 111.0)            # the grey panel behind the organ block

    # ---- the grey background panel, drawn FIRST so everything else sits on it
    fig.rect(*PANEL, rx=2.0, fill=CARD, stroke=GREY_L, lw=0.3)

    # ---- input 1: sleep duration
    k = BED_W / 20.0
    fig.place(SLEEPER, BX, Y_DUR, w=BED_W)
    fig.place("moon-green", BX + 15.5 * k, Y_DUR - 6.0 * k, w=7 * k)
    fig.place("clock-green", BX + 15.5 * k, Y_DUR + 3.5 * k, w=6 * k)
    fig.label("sleep duration\nshort or long", BX, Y_DUR + 11.0, anchor="middle", va="top",
              bold=True, color=GREEN_D, leading=1.15)

    # ---- input 2: PAP treatment
    fig.place(MASK, BX, Y_PAP, h=24.0)
    fig.label("PAP treatment", BX, Y_PAP + 14.0, anchor="middle", va="top", bold=True,
              color=INK)

    # ---- input 3: sleep apnea
    airflow(fig, BX, Y_APNEA - 10.0, 16 * k, 5.5 * k, APNEA)
    fig.place(SLEEPER, BX, Y_APNEA, w=BED_W)
    fig.label("sleep apnea", BX, Y_APNEA + 11.0, anchor="middle", va="top", bold=True,
              color=APNEA_D)

    # ---- all three act on one thing, the oxygen
    fig.arrow((AX, Y_DUR), (nx - 1.5, ncy - 10.0), color=AZ, lw=0.6)
    fig.arrow((AX, Y_PAP), (nx - 1.5, ncy), color=AZ, lw=0.6)
    fig.arrow((AX, Y_APNEA), (nx - 1.5, ncy + 10.0), color=AZ, lw=0.6)
    oxygen_node(fig, nx, ny, nw, nh, title=None, before_bed=True,
                trace_inset=(9.0, 7.0, 4.0, 9.5), lw=1.0)
    # the key, centred on the node and saying what it measures (Alen's round-20 ruling)
    swatch_label(fig, nx + nw / 2, KEY_Y, "time with oxygen saturation below 90%",
                 anchor="middle")

    # ---- TWO endings out of the oxygen box. At or below 1% of the night keeps the organs
    # in their own colours, above 1% darkens them. Same cut as Figures 1 and 3.
    fig.label(SCOPE, ORG_CX, 44.0, bold=True, anchor="middle", va="middle", color=GREY_D)
    bands = ((CY_TOP, False, LOW, AZ_T[2], -8.0, 64.0, AZ_T[0], AZ_D),
             (CY_BOT, True, HIGH, AZ, 8.0, 127.0, AZ_D, WHITE))
    for cy, dark, band, col, dy, band_y, pill_fill, pill_ink in bands:
        p0 = (nx + nw + 1.5, ncy + dy)
        p1 = (ORG_CX - 63.0, cy)
        fig.arrow(p0, p1, color=col, lw=1.1)
        organ_row(fig, ORG_CX, cy, dark=dark, pitch=PITCH, names=dark)
        # where the arrow is at the pill's own height, so the rule lands ON it
        on_arrow = p0[0] + (p1[0] - p0[0]) * (band_y - p0[1]) / (p1[1] - p0[1])
        band_pill(fig, 172.0, band_y, band, fill=pill_fill, text_color=pill_ink,
                  to_arrow=on_arrow, rule_color=col)

    # ---- the two roads that do not reach disease. Both end at the damaged organs, because
    # that is the ending each road fails to reach. ONE badge x for both, so they cannot
    # drift apart, and no correction line: the legend states the standard (round 20).
    big_crossed_road(fig, (AX + 8.0, ROAD_TOP), (ROAD_X, CY_BOT - 11.0),
                     badge_at=ROAD_BADGE_X, label=ROAD_DURATION)
    big_crossed_road(fig, (AX + 8.0, ROAD_BOT), (ROAD_X, CY_BOT + 11.0),
                     badge_at=ROAD_BADGE_X, label=ROAD_APNEA, label_dy=-7.0,
                     label_va="bottom")
    credit(fig)
    return finish(fig, stem, OUT, final=final, attribution=None)


SHEETS = {"band_fig4_bottom": fig4_bottom, "band_fig5_bottom": fig5_bottom,
          "Main_Fig6": fig6}

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--final", action="store_true")
    ap.add_argument("--only", default="")
    args = ap.parse_args()
    recs = [fn(args.final) for stem, fn in SHEETS.items()
            if not args.only or args.only in stem]
    (ROOT / "_build" / "records_r20_le.json").write_text(json.dumps(recs, indent=1))
    for r in recs:
        print(json.dumps({k: r[k] for k in ("stem", "w_pt", "h_mm", "audit", "min_pt")}))
