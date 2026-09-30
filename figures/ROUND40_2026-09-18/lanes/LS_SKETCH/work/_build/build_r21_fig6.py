#!/usr/bin/env python3
"""Round 21, lane LB: Figure 6.

  Main_Fig6   968.94 pt wide, 180 mm tall, both unchanged from round 20.

WHAT ROUND 21 CHANGES, and nothing else.

 4  THE TWO BAND LABELS ARE THE SAME BLUE AS FIGURE 1'S OXYGEN AZURE, #0187d1.
    Round 20 gave them two different fills, a pale #b3dcf2 pill under the top arrow and a
    dark #02598a pill under the bottom one. Both pills are now #0187d1, the exact hex
    Figure 1 prints for oxygen, with the label in white. The pill and its rule are kept
    exactly as round 20 built them, and the RULE still takes its colour from the arrow it
    lands on, so the two branches are still told apart by their arrows.

    The contrast this costs is in the report and in verify_r21.py, and it is not good news:
    white on #0187d1 is 3.90 to 1, which is under the 4.5 to 1 that normal text wants.
    #1a1d20 ink on the same blue is 4.34 to 1, also under. There is no text colour that
    clears 4.5 on this blue, so the choice is between the blue Alen named and the standard.
    He named the blue, so the blue is what is drawn, in white, which is the better of the
    two and the conventional one for a filled colour pill.

 5  THE GREY RECTANGLE BEHIND THE ORGANS IS SPLIT INTO TWO.
    Round 20 drew one 132.5 x 111 mm panel spanning both organ rows, most of it empty. It
    is now two rectangles, one behind the top row of coloured organs and one behind the
    bottom row of dark organs and its five names. They are a deliberate pair and they are
    IDENTICAL: same left edge 195.5, same right edge 328.0, same width 132.5 mm, same
    corner radius 2.0 mm, same height 29.254 mm, same fill and same hairline. Both x edges
    are round 20's own, so the branch arrows, the band pills and the two roads all keep the
    clearances round 20 tuned. The gap between them is 53.1 mm, which the fork sets.

    The height is the height the BOTTOM box needs, with 3.5 mm of clearance around its own
    content, and the top box is padded out to match. Equal boxes read as a pair at a glance,
    and the alternative, each box hugging its own row, was rendered and read side by side
    before this one was taken: it makes the top box a 20.5 mm bar against a 29.3 mm one,
    which reads as two different things rather than as two endings.

    THE ORGAN NAMES ARE WRAPPED at 20 mm so they fit inside their box. "heart and vessels"
    is 26.26 mm on one line, centred on the heart at x 205, so it used to begin at 191.87,
    which is 3.63 mm OUTSIDE a panel that starts at 195.5. Round 20's single panel had the
    same overhang and it was easy to miss on a panel 111 mm tall. On a box that hugs the row
    it is not. Wrapped, the name begins at 197.92, inside. This is the same two-line
    treatment Figure 1's top band already gives these five names. Widening the boxes to the
    left instead was tried and dropped: the band pill and its rule already run to x 192.4,
    so a box wide enough for the name would swallow them.

    The disease-count heading, "50 diseases and death", heads the PAIR, so it sits 3.2 mm
    above the top rectangle rather than inside it. It is the scope of both endings, not of
    the healthy row alone, and Figure 1's bottom band already prints it above its own organ
    block. Its y is read off the box that was actually drawn, never recomputed.

    Both branch arrows still enter their own rectangle through its left edge, which is what
    ties each ending to its panel.

Everything else on the sheet is untouched: the roads, the badges, the node, the key, the
three inputs, the two organ rows and the credit.

COLOUR LAW, unchanged from round 20: sleep apnea and breathing events are CHARCOAL #37474f
on the ladder #d3d7d8 #a5acb0 #737e84 #37474f with the dark text twin #29353b. Oxygen azure
#0187d1. Sleep duration green #39c445 with the ladder #c8eecb #92df99 #5dcf66 #298d32.
Ink #1a1d20. Arial, 9 pt floor. No em-dashes, no semicolons, and the banned words are never
printed.

    python3 build_r21_fig6.py            # 150 dpi preview into _build/previews/
    python3 build_r21_fig6.py --final    # SVG + 300 dpi PNG + PDF into the lane root
"""
from __future__ import annotations

import argparse
import json

from l1_common import (APNEA, APNEA_D, AZ, AZ_D, AZ_T, CARD, GREEN_D, GREY, GREY_D, GREY_L,
                       INK, MM_PER_PT, ROOT, WHITE, airflow, finish, new, swatch_label)
from build_l1 import oxygen_node
from r18_family import HIGH, LOW, ORGANS_LIVE, SCOPE, big_crossed_road, organ_row
from r20_common import ROAD_BADGE_X, band_pill
from roads_r18 import ROAD_APNEA, ROAD_DURATION

OUT = ROOT
CREDIT = "Icons: Servier Medical Art (CC BY 3.0)."
MASK = "cpap-mask-nostrap"                        # head + face + mask + hose, no strap
SLEEPER = "sleeper-faceB"                       # the SAME face as Figure 1

W6_PT, H6 = 968.94, 180.0                         # both unchanged from round 20

PITCH = 28.5                                      # organ pitch, LD centres 205 .. 319
ORG_CX = 262.0                                    # centre of the organ column

# ---------------------------------------------------------------- the two grey rectangles
PANEL_X, PANEL_W = 195.5, 132.5                   # round 20's own left and right edges
PANEL_RX = 2.0                                    # round 20's own corner radius
PANEL_PAD = 3.5                                   # the same clearance on every side of both
ICON_HALF = max(h for _n, ic in ORGANS_LIVE for _i, h, _d in ic) / 2      # 6.75 mm
# The bottom row's group names are wrapped so the widest of them, "heart and vessels" at
# 26.26 mm, stops overhanging the rectangle's left edge. Its centre is the heart at x 205
# and the rectangle starts at 195.5, so a single line of it would begin at 191.87, which is
# 3.63 mm OUTSIDE the panel. Wrapped it begins at 197.92, inside. 20.0 mm is the wrap width
# that breaks the two long names and leaves "metabolism" (17.13 mm) on one line. This is the
# same two-line treatment Figure 1's top band already gives these five names.
NAME_MAX_W = 20.0


def split_panels(fig, cy_top, cy_bot, bot_content_bottom, *, equal_height=False,
                 panel_x=PANEL_X, panel_w=PANEL_W, pad=PANEL_PAD, height=None):
    """The pair of grey rectangles, drawn as ONE function so they cannot drift apart. The
    top one wraps the coloured organ row, the bottom one wraps the dark row AND its five
    names. Same x, same width, same radius, same padding. Returns the two boxes and the gap
    between them, so the report can quote measured numbers rather than intentions."""
    top_h = 2 * (ICON_HALF + pad)
    bot_y0 = cy_bot - ICON_HALF - pad
    bot_h = bot_content_bottom + pad - bot_y0
    if height is not None:                         # both boxes forced to one height, the
        top_h = bot_h = height                     # bottom one anchored on its own content
        bot_y0 = bot_content_bottom + pad - height
    elif equal_height:
        top_h = bot_h                              # the top row is padded out to match
    top = (panel_x, cy_top - top_h / 2, panel_w, top_h)
    bot = (panel_x, bot_y0, panel_w, bot_h)
    for box in (top, bot):
        fig.rect(*box, rx=PANEL_RX, fill=CARD, stroke=GREY_L, lw=0.3)
    return {"top": [round(v, 3) for v in top], "bottom": [round(v, 3) for v in bot],
            "gap_mm": round(bot[1] - (top[1] + top[3]), 3)}


def credit(fig):
    fig.label(CREDIT, fig.w - 4.0, fig.h - 5.2, pt=fig.theme.floor_pt, anchor="end",
              va="top", color=GREY)


# ================================================================ Figure 6

def fig6(final, *, stem="Main_Fig6", name_max_w=NAME_MAX_W, equal_height=True,
         panel_x=PANEL_X, panel_w=PANEL_W, height=None):
    W = W6_PT * MM_PER_PT                          # 341.82 mm
    fig = new(W, H6)
    BX, BED_W = 42.0, 26.0
    Y_DUR, Y_PAP, Y_APNEA = 52.0, 95.0, 138.0
    AX = 70.0                                      # every input arrow starts here
    nx, ny, nw, nh = 112.0, 79.0, 58.0, 36.0       # the oxygen node
    ncy = ny + nh / 2                              # 97.0
    CY_TOP, CY_BOT = 56.0, 134.0
    ROAD_TOP, ROAD_BOT, ROAD_X = 17.0, 170.0, 336.0
    KEY_Y = 119.5                                  # the key, centred under the node

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
    swatch_label(fig, nx + nw / 2, KEY_Y, "time with oxygen saturation below 90%",
                 anchor="middle")

    # ---- TWO endings out of the oxygen box. At or below 1% of the night keeps the organs
    # in their own colours, above 1% darkens them. Same cut as Figures 1 and 3.
    # ROUND 21: both pills are #0187d1, Figure 1's oxygen azure, with the label in white.
    # The RULE keeps the arrow's own colour, so the two branches stay distinguishable.
    bands = ((CY_TOP, False, LOW, AZ_T[2], -8.0, 64.0),
             (CY_BOT, True, HIGH, AZ, 8.0, 127.0))
    bot_bottom = CY_BOT
    for cy, dark, band, col, dy, band_y in bands:
        p0 = (nx + nw + 1.5, ncy + dy)
        p1 = (ORG_CX - 63.0, cy)
        fig.arrow(p0, p1, color=col, lw=1.1)
        b = organ_row(fig, ORG_CX, cy, dark=dark, pitch=PITCH, names=dark,
                      name_max_w=name_max_w)
        if dark:
            bot_bottom = b
        # where the arrow is at the pill's own height, so the rule lands ON it
        on_arrow = p0[0] + (p1[0] - p0[0]) * (band_y - p0[1]) / (p1[1] - p0[1])
        band_pill(fig, 172.0, band_y, band, fill=AZ, text_color=WHITE,
                  to_arrow=on_arrow, rule_color=col)

    # ---- ROUND 21: the grey panel is TWO rectangles, one per ending. They are drawn last
    # so their exact heights can come from the real organ and name boxes, then moved to the
    # FRONT of the part list so they sit behind everything, which is where round 20's single
    # panel sat. fig._parts is the render order, first drawn is furthest back.
    n0 = len(fig._parts)
    panels = split_panels(fig, CY_TOP, CY_BOT, bot_bottom, equal_height=equal_height,
                          panel_x=panel_x, panel_w=panel_w, height=height)
    moved = fig._parts[n0:]
    del fig._parts[n0:]
    fig._parts[0:0] = moved

    # ---- the heading of the PAIR of rectangles, 3.2 mm above the TOP one. It is drawn
    # after the boxes so its y comes from the box that was actually built, never from a
    # second copy of the same arithmetic.
    SCOPE_Y = panels["top"][1] - 3.2
    fig.label(SCOPE, ORG_CX, SCOPE_Y, bold=True, anchor="middle", va="middle", color=GREY_D)

    # ---- the two roads that do not reach disease. ONE badge x for both, and no correction
    # line: the legend states the standard.
    big_crossed_road(fig, (AX + 8.0, ROAD_TOP), (ROAD_X, CY_BOT - 11.0),
                     badge_at=ROAD_BADGE_X, label=ROAD_DURATION)
    big_crossed_road(fig, (AX + 8.0, ROAD_BOT), (ROAD_X, CY_BOT + 11.0),
                     badge_at=ROAD_BADGE_X, label=ROAD_APNEA, label_dy=-7.0,
                     label_va="bottom")
    credit(fig)
    rec = finish(fig, stem, OUT, final=final, attribution=None)
    rec["panels"] = panels
    rec["scope_heading_y"] = round(SCOPE_Y, 3)
    return rec


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--final", action="store_true")
    ap.add_argument("--variants", action="store_true")
    args = ap.parse_args()
    if args.variants:
        for name, kw in (("V0_names_one_line", {"name_max_w": None}),
                         ("V1_names_wrapped", {}),
                         ("V2_equal_height", {"equal_height": True}),
                         ("V3_wide_one_line", {"name_max_w": None, "panel_x": 189.0,
                                               "panel_w": 138.8, "height": 30.0})):
            r = fig6(False, stem=f"Main_Fig6_{name}", **kw)
            print(name, json.dumps(r["panels"]))
        raise SystemExit
    rec = fig6(args.final)
    (ROOT / "_build" / "records_r21_fig6.json").write_text(json.dumps([rec], indent=1))
    print(json.dumps({k: rec[k] for k in ("stem", "w_pt", "h_mm", "audit", "min_pt",
                                          "panels", "scope_heading_y")}, indent=1))
