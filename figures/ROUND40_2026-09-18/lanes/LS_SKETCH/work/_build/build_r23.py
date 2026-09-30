#!/usr/bin/env python3
"""Round 23, lane LA: all SIX schematic sheets of the T90 paper.

  band_fig1_top      968.66 pt      band_fig4_bottom   970.79 pt
  band_fig1_bottom   968.66 pt      band_fig5_bottom   968.94 pt
  band_fig3_bottom   952.73 pt      Main_Fig6          968.94 pt

Every width above is exact and none of them moves. Every page height is the height round 22
delivered. Figures 1 top and 5 bottom are rebuilt unchanged so the set is one build.

WHAT ROUND 23 CHANGES, and nothing else.

 1  THE OUTCOME HEADING MOVES UNDER THE GREY ORGANS AS A RISK STATEMENT. On band_fig3_bottom
    and band_fig4_bottom the heading used to sit above the organ rows. It now sits below the
    dark-grey row and its five names, centred on the organ column, the way band_fig1_bottom
    already places its own line. The wording is a risk statement:
        band_fig3_bottom   "increased risk of 50 diseases and death"
        band_fig4_bottom   "increased risk of future disease"
        band_fig1_bottom   "increased risk of 50 diseases and death"  (wording only)
        band_fig5_bottom   unchanged
    Both moved bands are then re-balanced by lifting their whole content so the white above
    the first ink equals the white below the last ink. SHIFT3 and SHIFT4 are those lifts,
    solved on real ink at 600 dpi by --solve and pinned here so a plain build reproduces.
    The organ names stay where they were relative to the grey row. No type shrinks.

 2  GREEN AUDIT. Every green is #39c445 (fills) or #298d32 (text on white, icon outlines).
    Two icon overrides in _build/icons/symbols carried other greens and were corrected:
    the moon's rim shade (#298d32 at opacity 0.3 over #39c445, which composites to #34b43f)
    is deleted, and the clock face's #e6f8e8 at opacity 0.7 is written as the opaque tint
    #eefaef it composites to, which is #39c445 toward white at 0.914.

 3  MAIN_FIG6, six changes.
    a  the two pills sit at the SAME perpendicular offset from their arrows and the SAME
       distance along the arrow from the organ box. The oxygen node moves up 2 mm so its
       centre sits on the figure's own centre line, y = 95, where the three inputs and the
       two organ rows already centre. That makes the two branch arrows mirror images, so
       the pills are placed by one rule and the mirror does the rest.
    b  the ≤1% arrow and pill become the faint blue #7cc0e9 (Figure 2b's second rung), the
       pill's text in ink. The >1% arrow and pill keep azure #0187d1 with white text.
    c  the arrow tips TOUCH the organ boxes' left edge. The head is a triangle whose tip
       point is at p2 and whose stroke (half the line width, round joins) reaches 0.275 mm
       past it, so the tip point sits 0.275 mm left of the box edge and the ink stops on it.
    d  ONE dashed road, not two. It leaves the sleep-apnea input at the bottom left, runs up
       a spine to the left of the three inputs, picks up a feeder from the sleep-duration
       input, crosses the top of the sheet and drops onto the TOP edge of the healthy organ
       box. Both statements sit on the top leg. The bottom road is deleted.
    e  "increased risk of 50 diseases and death" sits under the GREY organ box. The heading
       over the healthy box is gone.
    f  the two grey rectangles, the azure key line under the oximetry box and the organ
       names under the grey row are all kept.
    The bottom road's deletion leaves the sheet bottom-heavy at 180 mm, so the content is
    shifted down by DY6 so the top white equals the bottom white. The height is unchanged.

COLOUR LAW. Sleep apnea and breathing events are CHARCOAL #37474f on the ladder #d3d7d8
#a5acb0 #737e84 #37474f with the dark text twin #29353b. Oxygen azure #0187d1. Green: fills
#39c445 or its tints toward white, text on white #298d32, icons fill #39c445 outline #298d32.
Ink #1a1d20. Arial, 9 pt floor. Nothing prints an em-dash, a semicolon, the two words about
groups, or the phrase about oxygen that Alen banned. No credit line on any sheet.

    python3 build_r23.py             # 150 dpi previews into _build/previews/
    python3 build_r23.py --final     # SVG + 300 dpi PNG + PDF into the lane root
    python3 build_r23.py --options   # the not-shipped road option into _build/previews/
    python3 build_r23.py --solve     # re-solve SHIFT3, SHIFT4 and DY6 on real ink
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from l1_common import (APNEA, APNEA_D, AZ, AZ_T, GREEN_D, GREY, GREY_D, INK, MM_PER_PT, ROOT,
                       WHITE, airflow, finish, new, swatch_label, text_width_mm)
from r14_common import trace_pair
from r18_common import (HIGH_1PCT, LOW_1PCT, ORGANS_BAD, ORGANS_OK, organ_names,
                        organ_row as organ_row_18)
from r18_family import HIGH, LOW, airflow_ramp, organ_row as organ_row_f, oxygen_fork
from r20_common import band_pill, severity_ladder
from build_l1 import oxygen_node
from build_r21_fig6 import PANEL_X, PANEL_W, NAME_MAX_W, split_panels
from bioglyph.arrows import _head

import build_r22 as B22                  # band_fig1_top, band_fig1_bottom (credits already off)
import build_r20_le as B45               # band_fig5_bottom (credit turned off by build_r22)

OUT = ROOT
PREVIEWS = ROOT / "_build" / "previews"
SCRATCH = ROOT / "_build" / "_solve"

# ---------------------------------------------------------------- the printed wording
RISK_50 = "increased risk of 50 diseases and death"
RISK_F6 = "increased risk of disease"   # R24: Fig 6 only; bands 1 and 3 keep RISK_50
RISK_F4 = "increased risk of future disease"
ROAD_DURATION = "short or long sleep alone: no added risk"
ROAD_APNEA = "sleep apnea alone: no added risk"
DURATIONS = "short, normal or long sleep"

FAINT = AZ_T[1]                          # "#7cc0e9", Figure 2b's second rung

# ---------------------------------------------------------------- widths and heights
W1_PT, W3_PT, W4_PT, W5_PT, W6_PT = 968.66, 952.73, 970.79, 968.94, 968.94
H3, H4, H6 = 70.0, 74.0, 180.0

SCOPE_GAP = 4.6                          # names' last ink line to the risk line (band_fig1_bottom's)

# the balance lifts, solved by --solve on real ink at 600 dpi and pinned here
SHIFT3 = 4.79                            # band_fig3_bottom, whole content up by this
SHIFT4 = 4.11                            # band_fig4_bottom, whole content up by this
DY6 = 4.63                               # Main_Fig6, whole content down by this

# ---------------------------------------------------------------- band_fig1_bottom wording
B22.SCOPE = RISK_50


# ================================================================ Figure 3, bottom band

def fig3_bottom(final, *, stem="band_fig3_bottom", shift=None, outdir=None):
    """Round 20's band with the heading moved under the grey row as a risk statement, and
    the whole content lifted by ``shift`` so the band is balanced on its page."""
    d = SHIFT3 if shift is None else shift
    W = W3_PT * MM_PER_PT
    fig = new(W, H3)
    t = fig.theme

    Y_OK, Y_BAD = 23.5 - d, 50.0 - d
    TX, TW, TH = 98.0, 44.0, 14.0
    Y_FLAT, Y_DIPS = Y_OK - TH / 2, Y_BAD - TH / 2
    LX = TX + TW + 3.2
    OCX, PITCH, NAME_GAP = 262.0, 28.5, 1.8
    T90_HEADER = "time with oxygen below 90%"

    bed = fig.place("sleeper-faceB", 36.0, 34.0 - d, w=54)
    fig.label(DURATIONS, 36.0, bed.y + bed.h + 0.6, bold=True, anchor="middle", va="top",
              color=GREEN_D)

    fig.label(T90_HEADER, (TX - 4.0 + LX + text_width_mm(HIGH_1PCT, t.label_pt, True)) / 2,
              13.0 - d, bold=True, color=AZ, anchor="middle", va="middle")
    for yt in (Y_FLAT, Y_DIPS):
        fig.arrow((62.0, 34.0 - d), (TX - 13.0, yt + TH * 0.5), color=GREY, lw=0.5)
    trace_pair(fig, TX, Y_FLAT, Y_DIPS, TW, TH, label_x=LX, lw=0.6,
               low=LOW_1PCT, high=HIGH_1PCT)
    brx = LX + text_width_mm(HIGH_1PCT, t.label_pt, True) + 4.0

    # the two endings. ROUND 23: no heading above, the risk line below the names.
    for yt in (Y_OK, Y_BAD):
        fig.arrow((brx, yt), (OCX - PITCH * 2 - 9.0, yt), color=GREY, lw=0.6)
    organ_row_18(fig, OCX, Y_OK, ORGANS_OK, pitch=PITCH, icon_scale=0.95)
    bad_bottom = organ_row_18(fig, OCX, Y_BAD, ORGANS_BAD, pitch=PITCH, icon_scale=0.95)[2]
    names_bottom = organ_names(fig, OCX, ORGANS_OK, bad_bottom + NAME_GAP, pitch=PITCH,
                               color=INK)
    risk_y = names_bottom + SCOPE_GAP
    lb = fig.label(RISK_50, OCX, risk_y, bold=True, anchor="middle", va="top", color=GREY_D)

    rec = finish(fig, stem, outdir or OUT, final=final, attribution=None)
    rec.update(shift=round(d, 3), names_bottom=round(names_bottom, 3),
               risk_y=round(risk_y, 3), risk_box=[round(lb.x, 3), round(lb.y, 3),
                                                  round(lb.x + lb.w, 3), round(lb.y + lb.h, 3)])
    return rec


# ================================================================ Figure 4, bottom band

def fig4_bottom(final, *, stem="band_fig4_bottom", shift=None, outdir=None):
    """Round 20's band with the heading moved under the grey row as a risk statement, in
    Alen's own words for this analysis, and the whole content lifted by ``shift``."""
    d = SHIFT4 if shift is None else shift
    W = W4_PT * MM_PER_PT
    fig = new(W, H4)
    t = fig.theme
    PITCH, ORG_CX, HEAD_Y = B45.PITCH, B45.ORG_CX, B45.HEAD_Y
    TX, TW, TH, LX = B45.TX, B45.TW, B45.TH, B45.LX

    RX, RW = 4.0, 74.0
    BX = RX + RW / 2
    fig.place(B45.SLEEPER, BX, 23.0 - d, w=28.0)
    airflow_ramp(fig, RX, 38.0 - d, RW, 8.5)
    severity_ladder(fig, RX, 46.0 - d, RW - 3.0, left_text=B45.APNEA_NONE,
                    right_text=B45.APNEA_SEVERE, tri_h=2.0)
    fig.label("sleep apnea severity", BX, 56.0 - d, anchor="middle", va="top", bold=True,
              color=APNEA_D)

    CY_TOP, CY_BOT = 23.5 - d, 52.5 - d
    oxygen_fork(fig, (84.0, 40.0 - d), TX, CY_TOP - TH / 2, CY_BOT - TH / 2, TW, TH,
                label_x=LX, dip_centers=(0.16, 0.40, 0.64, 0.86))
    lab_r = LX + max(text_width_mm(LOW, t.label_pt, True),
                     text_width_mm(HIGH, t.label_pt, True))
    fig.label(B45.T90_HEADER, (TX - 4.0 + lab_r) / 2, HEAD_Y - d, bold=True, color=AZ,
              anchor="middle", va="middle")

    # the two endings. ROUND 23: no heading above, the risk line below the names.
    from_x = lab_r + 2.5
    names_bottom = CY_BOT
    for cy, dark in ((CY_TOP, False), (CY_BOT, True)):
        fig.arrow((from_x, cy), (ORG_CX - 63.0, cy), color=GREY, lw=0.6)
        b = organ_row_f(fig, ORG_CX, cy, dark=dark, pitch=PITCH, names=dark)
        if dark:
            names_bottom = b
    risk_y = names_bottom + SCOPE_GAP
    lb = fig.label(RISK_F4, ORG_CX, risk_y, bold=True, anchor="middle", va="top",
                   color=GREY_D)

    rec = finish(fig, stem, outdir or OUT, final=final, attribution=None)
    rec.update(shift=round(d, 3), names_bottom=round(names_bottom, 3),
               risk_y=round(risk_y, 3), risk_box=[round(lb.x, 3), round(lb.y, 3),
                                                  round(lb.x + lb.w, 3), round(lb.y + lb.h, 3)])
    return rec


# ================================================================ Figure 6

def road(fig, pts, *, color=GREY, lw=1.0, dash="4.4 2.8", r=1.8):
    """One dashed road through ``pts`` with rounded corners and the engine's own triangle
    head at the last point. Drawn as ONE path so the dash pattern runs unbroken through
    every corner. Returns the head's trim so the caller can quote where the line stops."""
    (xa, ya), (xb, yb) = pts[-2], pts[-1]
    dd = math.hypot(xb - xa, yb - ya)
    ux, uy = (xb - xa) / dd, (yb - ya) / dd
    h_svg, trim = _head((xb, yb), (ux, uy), lw, "triangle", color, None)
    end = (xb - ux * trim, yb - uy * trim)
    d = f"M{pts[0][0]:.3f} {pts[0][1]:.3f}"
    for i in range(1, len(pts) - 1):
        p, c, n = pts[i - 1], pts[i], pts[i + 1]
        d1 = math.hypot(c[0] - p[0], c[1] - p[1])
        d2 = math.hypot(n[0] - c[0], n[1] - c[1])
        u1 = ((c[0] - p[0]) / d1, (c[1] - p[1]) / d1)
        u2 = ((n[0] - c[0]) / d2, (n[1] - c[1]) / d2)
        rr = min(r, d1 / 2, d2 / 2)
        a1 = (c[0] - u1[0] * rr, c[1] - u1[1] * rr)
        a2 = (c[0] + u2[0] * rr, c[1] + u2[1] * rr)
        d += f" L{a1[0]:.3f} {a1[1]:.3f} Q{c[0]:.3f} {c[1]:.3f} {a2[0]:.3f} {a2[1]:.3f}"
    d += f" L{end[0]:.3f} {end[1]:.3f}"
    fig.path(d, stroke=color, lw=lw, dash=dash)
    fig._strokes_mm.append(lw)
    fig.add(h_svg)
    return trim


def pill_geometry(text, p0, tip, band_y, gap, pt, pad_x=3.0):
    """Where a band pill sits so that its right edge stops ``gap`` mm short of the arrow's
    centre line at the pill's own height. Returns (cx, x_on_arrow, d_perp, d_along)."""
    w = text_width_mm(text, pt, True) + 2 * pad_x
    x_on = p0[0] + (tip[0] - p0[0]) * (band_y - p0[1]) / (tip[1] - p0[1])
    cx = x_on - gap - w / 2
    # the pill centre against the arrow: perpendicular offset, and distance from the tip
    # along the arrow's own direction
    dx, dy = tip[0] - p0[0], tip[1] - p0[1]
    L = math.hypot(dx, dy)
    ux, uy = dx / L, dy / L
    vx, vy = cx - tip[0], band_y - tip[1]
    d_along = -(vx * ux + vy * uy)
    d_perp = abs(vx * uy - vy * ux)
    return cx, x_on, d_perp, d_along, w


def fig6(final, *, stem="Main_Fig6", dy=None, roads="one_x", spine_x=18.0, feeder_x0=27.0,
         badge_xs=(90.0, 205.0), pill_dy=31.0, rule_gap=3.0, outdir=None):
    """roads='one_x'  one shared cross badge on the top leg, the two statements stacked under
                      it (SHIPPED)
       roads='two_x'  two cross badges on the top leg, one statement under each (the option
                      kept in _build/previews/)"""
    d = DY6 if dy is None else dy
    W = W6_PT * MM_PER_PT
    fig = new(W, H6)
    t = fig.theme

    def Y(y):
        return y + d

    BX, BED_W = 42.0, 26.0
    Y_DUR, Y_PAP, Y_APNEA = Y(52.0), Y(95.0), Y(138.0)
    AX = 70.0
    MID = Y(95.0)                                  # the figure's own centre line
    nw, nh = 58.0, 36.0
    nx, ny = 112.0, MID - nh / 2                   # ROUND 23: the node centres on MID
    ncy = ny + nh / 2
    CY_TOP, CY_BOT = Y(56.0), Y(134.0)
    KEY_Y = ny + nh + 4.5                          # the key keeps its 4.5 mm under the node
    ROAD_TOP = Y(17.0)
    ORG_CX, PITCH = B45.ORG_CX, B45.PITCH
    ARROW_LW = 1.1
    TIP_X = PANEL_X - ARROW_LW * 0.25              # the head's stroke reaches the box edge

    # ---- input 1: sleep duration
    k = BED_W / 20.0
    fig.place(B45.SLEEPER, BX, Y_DUR, w=BED_W)
    # ROUND 27: the green moon and the green clock that sat at BX + 15.5 k beside this bed are
    # REMOVED, so the input reads as the bed and its label alone, like the other two inputs.
    # The bed, the label, the road feeder and the arrow keep their coordinates. The two boxes
    # the icons occupied are recorded (MOON_BOX, CLOCK_BOX) so the diff can be confined to them.
    MOON_BOX = [BX + 15.5 * k - 3.5 * k, Y_DUR - 6.0 * k - 3.5 * k, 7 * k, 7 * k]
    CLOCK_BOX = [BX + 15.5 * k - 3.0 * k, Y_DUR + 3.5 * k - 3.0 * k, 6 * k, 6 * k]
    fig.label("sleep duration\nshort or long", BX, Y_DUR + 11.0, anchor="middle", va="top",
              bold=True, color=GREEN_D, leading=1.15)

    # ---- input 2: PAP treatment
    fig.place(B45.MASK, BX, Y_PAP, h=24.0, flip_h=True)  # ROUND 27: the head faces right
    fig.label("PAP treatment", BX, Y_PAP + 14.0, anchor="middle", va="top", bold=True,
              color=INK)

    # ---- input 3: sleep apnea
    airflow(fig, BX, Y_APNEA - 10.0, 16 * k, 5.5 * k, APNEA)
    fig.place(B45.SLEEPER, BX, Y_APNEA, w=BED_W)
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

    # ---- TWO endings out of the oxygen box, mirror images about MID. ROUND 23: the ≤1%
    # branch is the faint blue with ink in its pill, the >1% branch keeps azure and white.
    # Each pill's right edge stops rule_gap short of its own arrow at the pill's height, and
    # the mirror makes the two offsets identical.
    bands = ((CY_TOP, False, LOW, FAINT, -10.0, MID - pill_dy, INK),
             (CY_BOT, True, HIGH, AZ, 10.0, MID + pill_dy, WHITE))
    bot_bottom = CY_BOT
    pills, arrows = {}, {}
    for cy, dark, band, col, dyy, band_y, pill_ink in bands:
        p0 = (nx + nw + 1.5, ncy + dyy)
        p1 = (TIP_X, cy)
        fig.arrow(p0, p1, color=col, lw=ARROW_LW)
        b = organ_row_f(fig, ORG_CX, cy, dark=dark, pitch=PITCH, names=dark,
                        name_max_w=NAME_MAX_W)
        if dark:
            bot_bottom = b
        cx, x_on, d_perp, d_along, pw = pill_geometry(band, p0, p1, band_y, rule_gap,
                                                      t.label_pt)
        band_pill(fig, cx, band_y, band, fill=col, text_color=pill_ink, to_arrow=x_on,
                  rule_color=col)
        key = "low" if not dark else "high"
        arrows[key] = {"p0": [round(v, 3) for v in p0], "tip": [round(v, 3) for v in p1],
                       "color": col, "lw": ARROW_LW,
                       "tip_ink_x": round(p1[0] + ARROW_LW * 0.25, 3)}
        pills[key] = {"cx": round(cx, 3), "cy": round(band_y, 3), "w": round(pw, 3),
                      "x_on_arrow": round(x_on, 3), "rule_mm": round(x_on - (cx + pw / 2), 3),
                      "d_perp": round(d_perp, 3), "d_along": round(d_along, 3),
                      "fill": col, "text": pill_ink}

    # ---- the two grey rectangles, drawn last and moved to the back (round 21)
    n0 = len(fig._parts)
    panels = split_panels(fig, CY_TOP, CY_BOT, bot_bottom, equal_height=True)
    moved = fig._parts[n0:]
    del fig._parts[n0:]
    fig._parts[0:0] = moved

    # ---- ROUND 23 e: the risk statement under the GREY box, read off the box drawn
    bx, by, bw, bh = panels["bottom"]
    RISK_Y = by + bh + 3.2
    lb = fig.label(RISK_F6, ORG_CX, RISK_Y, bold=True, anchor="middle", va="top",
                   color=GREY_D)

    # ---- ROUND 23 d: ONE road. From the sleep-apnea input, up the spine left of the
    # inputs, a feeder from the sleep-duration input, across the top, down onto the TOP
    # edge of the healthy box. The head's stroke reaches 0.25 mm past its tip point.
    top_y0 = panels["top"][1]
    ROAD_LW = 1.0
    tip = (ORG_CX, top_y0 - ROAD_LW * 0.25)
    pts = [(feeder_x0, Y_APNEA), (spine_x, Y_APNEA), (spine_x, ROAD_TOP),
           (ORG_CX, ROAD_TOP), tip]
    trim = road(fig, pts, color=GREY, lw=ROAD_LW)
    fig.path(f"M{feeder_x0:.3f} {Y_DUR:.3f} L{spine_x:.3f} {Y_DUR:.3f}", stroke=GREY,
             lw=ROAD_LW, dash="4.4 2.8")
    R = 4.8
    if roads == "two_x":
        for bxx, text in zip(badge_xs, (ROAD_DURATION, ROAD_APNEA)):
            fig.badge(bxx, ROAD_TOP, R, sign="×", color=GREY, glyph_color=WHITE)
            fig.label(text, bxx, ROAD_TOP + 7.0, pt=t.title_pt, bold=True, anchor="middle",
                      va="top", color=GREY)
        badge_at = list(badge_xs)
    else:
        bxx = (spine_x + ORG_CX) / 2
        fig.badge(bxx, ROAD_TOP, R, sign="×", color=GREY, glyph_color=WHITE)
        fig.label(ROAD_DURATION + "\n" + ROAD_APNEA, bxx, ROAD_TOP + 7.0, pt=t.title_pt,
                  bold=True, anchor="middle", va="top", color=GREY, leading=1.2)
        badge_at = [bxx]

    rec = finish(fig, stem, outdir or OUT, final=final, attribution=None)
    rec.update(dy=round(d, 3), node=[nx, round(ny, 3), nw, nh], node_cy=round(ncy, 3),
               removed_icons={"moon-green": [round(v, 3) for v in MOON_BOX],
                              "clock-green": [round(v, 3) for v in CLOCK_BOX]},
               mid=round(MID, 3), key_y=round(KEY_Y, 3), panels=panels, pills=pills,
               arrows=arrows, panel_x=PANEL_X, tip_x=round(TIP_X, 3),
               risk_y=round(RISK_Y, 3),
               risk_box=[round(lb.x, 3), round(lb.y, 3), round(lb.x + lb.w, 3),
                         round(lb.y + lb.h, 3)],
               road={"points": [[round(a, 3), round(b, 3)] for a, b in pts],
                     "feeder": [[feeder_x0, round(Y_DUR, 3)], [spine_x, round(Y_DUR, 3)]],
                     "tip": [round(tip[0], 3), round(tip[1], 3)],
                     "tip_ink_y": round(tip[1] + ROAD_LW * 0.25, 3),
                     "top_box_y0": round(top_y0, 3), "head_trim": round(trim, 3),
                     "badges_x": badge_at, "badge_y": round(ROAD_TOP, 3), "kind": roads})
    return rec


# ================================================================ the six, in one place

def build_all(final):
    return [B22.fig1_top(final), B22.fig1_bottom(final),
            fig3_bottom(final), fig4_bottom(final), B45.fig5_bottom(final),
            fig6(final)]


# ================================================================ the balance solver

def solve(fn, key, start, **kw):
    """Lift or drop a sheet's content until the white above the first ink equals the white
    below the last ink, measured on real ink at 600 dpi. Two passes converge because the
    content moves rigidly."""
    from measure_r22 import whites
    SCRATCH.mkdir(parents=True, exist_ok=True)
    v = start
    for _ in range(3):
        rec = fn(True, stem=f"solve_{key}", outdir=SCRATCH, **{key: v})
        w = whites(SCRATCH / f"solve_{key}.pdf")
        err = w["top"] - w["bottom"]
        print(f"  {key}={v:.3f}  top {w['top']:.3f}  bottom {w['bottom']:.3f}  err {err:+.3f}")
        if abs(err) < 0.03:
            break
        v = v + err / 2 if key == "shift" else v - err / 2
    return round(v, 2)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--final", action="store_true")
    ap.add_argument("--only", default="")
    ap.add_argument("--options", action="store_true")
    ap.add_argument("--solve", action="store_true")
    args = ap.parse_args()

    if args.solve:
        print("band_fig3_bottom"); s3 = solve(fig3_bottom, "shift", SHIFT3)
        print("band_fig4_bottom"); s4 = solve(fig4_bottom, "shift", SHIFT4)
        print("Main_Fig6"); d6 = solve(fig6, "dy", DY6)
        print(json.dumps({"SHIFT3": s3, "SHIFT4": s4, "DY6": d6}))
        raise SystemExit

    if args.options:
        # NOT SHIPPED: the top road with TWO cross badges, one statement under each
        r = fig6(False, stem="opt_Main_Fig6_two_x", roads="two_x")
        print(json.dumps({k: r[k] for k in ("stem", "audit", "road")}))
        raise SystemExit

    recs = [r for r in build_all(args.final) if not args.only or args.only in r["stem"]]
    (ROOT / "_build" / "records_r23.json").write_text(json.dumps(recs, indent=1))
    for r in recs:
        print(json.dumps({k: r[k] for k in ("stem", "w_pt", "h_mm", "audit", "min_pt")}))
