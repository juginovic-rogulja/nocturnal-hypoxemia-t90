#!/usr/bin/env python3
"""Round-12 L1: FigureLab schematics for the T90 paper.

A. "all roads run through oxygen", five iterations at 183 mm (Alen picks one)
B. opening and closing bands for Figs 1, 3 and 5 at the exact page widths
C. home-recording band for the external-validation supplementary sheet (496.8 pt)

    python3 build_l1.py                 # 150 dpi previews into _build/previews/ + records.json
    python3 build_l1.py --final         # SVG + 300 dpi PNG + PDF into L1_schematics/
    python3 build_l1.py --only band_fig3   # sheets whose stem starts with the prefix
"""
from __future__ import annotations

import argparse
import json

from l1_common import (AZ, AZ_D, AZ_T, CARD, GREEN, GREEN_D, GREY, GREY_L, INK, MM_PER_PT,
                       ORGANS, PREVIEWS, ROOT, TERRA, WHITE, airflow, answer, crossed_road,
                       finish, house, new, night_axis, number_circle, organ_row_named,
                       rank_ladder, spo2_trace, swatch_label, text_width_mm)

OUT = ROOT
W_FIG1, W_FIG3, W_FIG5 = 968.7 * MM_PER_PT, 952.7 * MM_PER_PT, 968.9 * MM_PER_PT
W_HOME, H_HOME = 496.8 * MM_PER_PT, 79.0 * MM_PER_PT


# ==================================================================== A. all roads
# Shared building blocks: the left column of sleepers, the oxygen node, the organ
# card and the two crossed grey roads Alen liked in sketchC.

def sleeper_short(fig, bx, y, *, bed_w=20.0, label_dy=9.0, pt=None):
    fig.place("sleeper-bed", bx, y, w=bed_w)
    k = bed_w / 20.0
    fig.place("moon-green", bx + 15.5 * k, y - 6.0 * k, w=7 * k)
    fig.place("clock-green", bx + 15.5 * k, y + 3.5 * k, w=6 * k)
    return fig.label("short sleep", bx, y + label_dy, anchor="middle", va="top", bold=True,
                     color=GREEN_D, pt=pt)


def sleeper_apnea(fig, bx, y, *, bed_w=20.0, label_dy=9.0, pt=None):
    k = bed_w / 20.0
    airflow(fig, bx, y - 10.0 * k, 14 * k, 5 * k, TERRA)
    fig.place("sleeper-bed", bx, y, w=bed_w)
    return fig.label("sleep apnea", bx, y + label_dy, anchor="middle", va="top", bold=True,
                     color=TERRA, pt=pt)


def pap_person(fig, bx, y, *, w=16.0, label_dy=10.0, text="on PAP treatment"):
    fig.place("cpap-mask", bx, y, w=w)
    return fig.label(text, bx, y + label_dy, anchor="middle", va="top", bold=True, color=INK)


def oxygen_node(fig, nx, ny, nw, nh, *, title="low oxygenation during sleep", before_bed=False,
                night=True, title_dy=3.0, trace_inset=(8.0, 3.0, 4.0, 9.0), lw=0.9):
    """White node with an azure border and the dipping overnight trace inside.
    before_bed=True draws a grey awake segment (crossed) before the sleep segment."""
    fig.rect(nx, ny, nw, nh, rx=3, fill=WHITE, stroke=AZ, lw=lw)
    li, ti, ri, bi = trace_inset
    tx, ty, tw, th = nx + li, ny + ti, nw - li - ri, nh - ti - bi
    t = fig.theme
    if before_bed:
        split = tx + tw * 0.32
        pts = [(tx + (split - tx) * i / 40, ty + th * 0.27 + 0.35 * ((i % 7) - 3) / 3) for i in range(41)]
        fig.path("M" + " L".join(f"{px:.2f} {py:.2f}" for px, py in pts), stroke=GREY, lw=0.55)
        fig.line(split, ty - 1.0, split, ty + th + 0.5, stroke=GREY_L, lw=0.35, dash="1.0 0.8")
        spo2_trace(fig, tx, ty, tw, th, dips=True, color=AZ, ref_at=0.5,
                   dip_centers=(0.18, 0.36, 0.55, 0.73, 0.91), label90="left", x_start=split)
        fig.badge(tx + (split - tx) / 2 - 6.4, ty - 3.0, 1.5, sign="×", color=GREY, glyph_color=WHITE)
        fig.label("awake", tx + (split - tx) / 2 - 3.8, ty - 2.9, pt=t.floor_pt, color=GREY, va="middle")
        fig.label("asleep", split + (tx + tw - split) / 2, ty - 2.9, pt=t.floor_pt, color=AZ,
                  anchor="middle", va="middle")
    else:
        spo2_trace(fig, tx, ty, tw, th, dips=True, color=AZ, ref_at=0.5, label90="left")
    if night:
        night_axis(fig, tx, ty + th + 1.2, tw, pt=t.floor_pt)
    lb = None
    if title:
        lb = fig.label(title, nx + nw / 2, ny + nh + title_dy, anchor="middle", va="top", bold=True,
                       color=AZ)
    return lb


def organ_card(fig, x, y, w, h, *, title="future disease", pitch=15.0, icon_h=8.5, organ_y=None,
               title_dy=2.5, name_pt=None):
    fig.card(x, y, w, h, fill=CARD)
    cm = x + w / 2
    if title:
        fig.label(title, cm, y + title_dy, anchor="middle", va="top", bold=True, color=INK)
    oy = organ_y if organ_y is not None else y + h * 0.42
    return organ_row_named(fig, cm, oy, pitch=pitch, icon_h=icon_h, max_w=pitch - 1.0, name_pt=name_pt)


def iter1(final, *, before_bed=False, rank_badge=False, stem="iteration_1"):
    """PAP as a third input road on the left that forks: oxygen improved (check,
    fewer diseases) or oxygen still low (into the node, on to the organs)."""
    W, H = 183.0, 108.0
    fig = new(W, H)
    t = fig.theme
    BX = 24.0
    Y_SHORT, Y_PAP, Y_APNEA = 19.0, 56.0, 91.0
    if before_bed:
        nx, ny, nw, nh = 71.0, 40.0, 45.0, 28.0
        cxb, cyb, cwb, chb, pitch = 123.0, 37.0, 58.0, 32.0, 14.5
        CHK = 56.5
    else:
        nx, ny, nw, nh = 75.0, 40.0, 38.0, 27.0
        cxb, cyb, cwb, chb, pitch = 121.0, 37.0, 60.0, 32.0, 15.0
        CHK = 60.0
    ncy = ny + nh / 2
    cmid = cxb + cwb / 2

    sleeper_short(fig, BX, Y_SHORT)
    pap_person(fig, BX, Y_PAP)
    sleeper_apnea(fig, BX, Y_APNEA)

    fig.arrow((45.0, Y_SHORT), (nx + 6.0, ny - 1.0), color=AZ, lw=0.4)
    fig.arrow((45.0, Y_APNEA), (nx - 1.0, ny + nh - 5.0), color=AZ, lw=0.4)
    F = (36.0, Y_PAP)
    fig.line(33.0, Y_PAP, F[0], F[1], stroke=AZ, lw=0.4)
    fig.arrow(F, (nx - 1.0, ncy + 3.5), color=AZ, lw=0.4)
    fig.label("oxygen still low", (F[0] + nx) / 2, ncy - 1.2, pt=t.note_pt, anchor="middle",
              va="bottom", color=AZ)
    fig.arrow(F, (CHK - 4.0, 40.0), color=AZ_T[2], lw=0.4, elbow="v")
    fig.badge(CHK, 40.0, 2.6, sign="check", color=AZ, glyph_color=WHITE)
    fig.label("oxygen improved", CHK - 13.0, 37.6, pt=t.note_pt, anchor="middle", va="bottom", color=AZ)
    fig.label("fewer diseases", CHK, 43.6, pt=t.note_pt, anchor="middle", va="top", color=AZ)

    title = "low oxygenation while asleep" if before_bed else "low oxygenation during sleep"
    oxygen_node(fig, nx, ny, nw, nh, title=title, before_bed=before_bed,
                trace_inset=(8.0, 7.0, 4.0, 8.0) if before_bed else (8.0, 3.0, 4.0, 9.0))
    yy = ny + nh + 7.5
    if before_bed:
        fig.label("oxygen before bed: no link to disease", nx + nw / 2, yy, pt=t.note_pt,
                  anchor="middle", va="top", color=GREY)
        yy += 5.5
    swatch_label(fig, nx + nw / 2, yy + 1.6, "time below 90% (T90)" if rank_badge else "time below 90%",
                 anchor="middle")
    if rank_badge:
        yy += 6.0
        rank_ladder(fig, nx + nw / 2 - 26.0, yy - 0.5, 4.0, 6.0)
        fig.label("T90 ranks 1 of 197 sleep measurements", nx + nw / 2 - 19.5, yy + 2.5,
                  pt=t.note_pt, va="middle", color=AZ)

    fig.arrow((nx + nw + 1.0, ncy), (cxb - 1.0, ncy), color=AZ, lw=1.0)
    organ_card(fig, cxb, cyb, cwb, chb, pitch=pitch, organ_y=cyb + 13.5)

    crossed_road(fig, (46.0, 7.0), (cmid, cyb - 1.0), badge_at=98.5,
                 label="short sleep alone: no extra disease")
    crossed_road(fig, (46.0, 104.0), (cmid, cyb + chb + 1.0), badge_at=98.5,
                 label="sleep apnea alone: no extra disease", label_dy=-3.6, label_va="bottom")
    return finish(fig, stem, OUT, final=final)


def iter2(final):
    """Right-hand column: the organ card on top, the two PAP outcomes under it."""
    W, H = 183.0, 118.0
    fig = new(W, H)
    t = fig.theme
    BX = 24.0
    Y_SHORT, Y_ANY, Y_APNEA = 19.0, 46.0, 74.0
    nx, ny, nw, nh = 71.0, 34.0, 38.0, 27.0
    ncy = ny + nh / 2
    cxb, cyb, cwb, chb = 118.0, 31.0, 63.0, 32.0
    cmid = cxb + cwb / 2

    sleeper_short(fig, BX, Y_SHORT)
    fig.place("sleeper-bed", BX, Y_ANY, w=20)
    fig.label("no sleep apnea", BX, Y_ANY + 9.0, anchor="middle", va="top", bold=True, color=GREY)
    sleeper_apnea(fig, BX, Y_APNEA)
    for y_from, y_to in ((Y_SHORT, ncy - 6), (Y_ANY, ncy), (Y_APNEA, ncy + 6)):
        fig.arrow((45.0, y_from), (nx - 1.0, y_to), color=AZ, lw=0.4)

    oxygen_node(fig, nx, ny, nw, nh)
    swatch_label(fig, nx + nw / 2, ny + nh + 9.0, "time below 90%", anchor="middle")
    fig.arrow((nx + nw + 1.0, ncy), (cxb - 1.0, ncy), color=AZ, lw=1.0)
    organ_card(fig, cxb, cyb, cwb, chb, pitch=15.5, organ_y=cyb + 13.5)

    crossed_road(fig, (46.0, 7.0), (cmid, cyb - 1.0), badge_at=98.5,
                 label="short sleep alone: no extra disease")
    crossed_road(fig, (46.0, 90.0), (cxb + 4.0, cyb + chb + 1.0), badge_at=84.0,
                 label="sleep apnea alone: no extra disease")

    # PAP outcomes under the organs (header right-aligned so the climbing road clears it)
    px = cxb
    fig.label("on PAP treatment", cxb + cwb, 70.0, anchor="end", va="top", bold=True, color=INK)
    rows = ((83.5, False, "oxygen improved", "fewer diseases"),
            (105.0, True, "oxygen still low", "still reaches the organs"))
    LX = px + 31.0
    for ry, dips, l1, l2 in rows:
        fig.place("cpap-mask", px + 6.0, ry, w=10)
        spo2_trace(fig, px + 13.0, ry - 6.0, 13.5, 12.0, dips=dips, color=AZ if dips else AZ_T[2],
                   ref_at=0.55, dip_centers=(0.2, 0.5, 0.8), label90=None, lw=0.5)
        fig.label(l1, LX, ry - 4.6, pt=t.note_pt, va="top", color=AZ)
        if dips:
            fig.label(l2, LX, ry + 1.0, pt=t.note_pt, va="top", color=AZ, max_w=33.0, leading=1.15)
        else:
            answer(fig, LX, ry + 3.0, l2, sign="check", color=AZ, r=1.9, pt=t.note_pt,
                   anchor="start", bold=False, text_color=AZ)
    fig.arrow((px + 28.7, 99.0), (px + 28.7, cyb + chb + 1.0), color=AZ, lw=0.6)
    return finish(fig, "iteration_2", OUT, final=final)


def iter3(final):
    """Compact: three inputs stacked, one oxygen gate with two exits."""
    W, H = 183.0, 88.0
    fig = new(W, H)
    t = fig.theme
    BX = 22.0
    Y_SHORT, Y_PAP, Y_APNEA = 12.0, 40.0, 66.0
    nx, ny, nw, nh = 66.0, 30.0, 36.0, 25.0
    ncy = ny + nh / 2
    cxb, cyb, cwb, chb = 120.0, 34.0, 61.0, 32.0
    cmid = cxb + cwb / 2

    sleeper_short(fig, BX, Y_SHORT, bed_w=16, label_dy=6.5)
    pap_person(fig, BX, Y_PAP, w=13, label_dy=7.5)
    sleeper_apnea(fig, BX, Y_APNEA, bed_w=16, label_dy=6.5)
    for y_from, y_to in ((Y_SHORT, ny + 3.0), (Y_PAP, ncy), (Y_APNEA, ny + nh - 4.0)):
        fig.arrow((40.0, y_from), (nx - 1.0, y_to), color=AZ, lw=0.4)

    oxygen_node(fig, nx, ny, nw, nh, title="how low does oxygen fall during sleep?", night=False,
                trace_inset=(8.0, 3.0, 3.0, 4.0))
    swatch_label(fig, nx + nw / 2, ny + nh + 9.5, "time below 90%", anchor="middle")

    fig.arrow((nx + nw + 1.0, ny + 1.0), (cxb + 1.0, 19.0), color=AZ_T[2], lw=0.5)
    fig.label("preserved, or improved on PAP", cxb + 2.0, 12.0, pt=t.note_pt, va="bottom", color=AZ)
    answer(fig, cxb + 2.0, 19.0, "no extra disease", sign="check", color=AZ, r=2.4, anchor="start")
    fig.arrow((nx + nw + 1.0, ncy + 4.0), (cxb - 1.0, cyb + chb / 2), color=AZ, lw=1.0)
    organ_card(fig, cxb, cyb, cwb, chb, title="oxygen stays low: future disease", pitch=15.0,
               organ_y=cyb + 13.5)

    crossed_road(fig, (40.0, 4.0), (cxb + cwb - 5.0, cyb - 1.0), badge_at=80.0,
                 label="short sleep alone: no extra disease")
    crossed_road(fig, (40.0, 84.0), (cmid, cyb + chb + 1.0), badge_at=95.0,
                 label="sleep apnea alone: no extra disease", label_dy=-3.6, label_va="bottom")
    return finish(fig, "iteration_3", OUT, final=final)


def fig5_closing(final):
    """The concluding schematic: iteration 1 re-laid at the Fig 5 page width."""
    W, H = W_FIG5, 98.0
    fig = new(W, H)
    t = fig.theme
    BX = 34.0
    Y_SHORT, Y_PAP, Y_APNEA = 18.0, 50.0, 84.0
    nx, ny, nw, nh = 150.0, 33.0, 52.0, 32.0
    ncy = ny + nh / 2
    cxb, cyb, cwb, chb = 240.0, 29.0, 96.0, 40.0
    cmid = cxb + cwb / 2

    sleeper_short(fig, BX, Y_SHORT, bed_w=26, label_dy=10.5)
    pap_person(fig, BX, Y_PAP, w=20, label_dy=11.5)
    sleeper_apnea(fig, BX, Y_APNEA, bed_w=26, label_dy=10.5)

    fig.arrow((62.0, Y_SHORT), (nx - 1.0, ny + 4.0), color=AZ, lw=0.5)
    fig.arrow((62.0, Y_APNEA), (nx - 1.0, ny + nh - 6.0), color=AZ, lw=0.5)
    F = (64.0, Y_PAP)
    fig.line(45.0, Y_PAP, F[0], F[1], stroke=AZ, lw=0.5)
    fig.arrow(F, (nx - 1.0, ncy + 4.0), color=AZ, lw=0.5)
    fig.label("oxygen still low", (F[0] + nx) / 2, ncy - 1.5, pt=t.note_pt, anchor="middle",
              va="bottom", color=AZ)
    fig.arrow(F, (96.0, 35.0), color=AZ_T[2], lw=0.5, elbow="v")
    fig.badge(100.5, 35.0, 3.0, sign="check", color=AZ, glyph_color=WHITE)
    fig.label("oxygen improved", 80.0, 32.2, pt=t.note_pt, anchor="middle", va="bottom", color=AZ)
    fig.label("fewer diseases", 105.5, 35.1, pt=t.note_pt, va="middle", color=AZ)

    oxygen_node(fig, nx, ny, nw, nh, trace_inset=(9.0, 3.5, 4.0, 10.0), lw=1.0)
    swatch_label(fig, nx + nw / 2, ny + nh + 9.5, "time below 90%", anchor="middle")
    fig.arrow((nx + nw + 1.5, ncy), (cxb - 1.5, ncy), color=AZ, lw=1.3)
    organ_card(fig, cxb, cyb, cwb, chb, pitch=23.0, icon_h=11.0, organ_y=cyb + 17.0)

    crossed_road(fig, (62.0, 6.0), (cmid, cyb - 1.5), badge_at=176.0,
                 label="short sleep alone: no extra disease", r=2.8)
    crossed_road(fig, (62.0, 94.0), (cmid, cyb + chb + 1.5), badge_at=176.0,
                 label="sleep apnea alone: no extra disease", label_dy=-3.8, label_va="bottom", r=2.8)
    return finish(fig, "band_fig5_bottom", OUT, final=final, attribution=None)


# ==================================================================== B. opening bands

def fig1_opening(final):
    W, H = W_FIG1, 66.0
    fig = new(W, H)
    t = fig.theme
    CAP_Y = 59.0
    number_circle(fig, 1, 6.0, 8.0)
    fig.place("sleeper-bed", 36.0, 33.0, w=46)
    fig.place("psg-head", 75.0, 25.0, w=17)
    fig.label("brain-wave sensors", 75.0, 35.0, pt=t.note_pt, anchor="middle", va="top", color=INK)
    fig.place("oximeter-azure", 75.0, 46.0, w=15)
    fig.label("finger oxygen sensor", 75.0, 52.0, pt=t.note_pt, anchor="middle", va="top", color=AZ_D)
    fig.label("One night of sleep recording in 19,173 adults", 50.0, CAP_Y, anchor="middle",
              va="top", color=INK)
    fig.block_arrow(97.0, 27.0, 12.0, 12.0)

    number_circle(fig, 2, 116.0, 8.0)
    fig.card(112.0, 12.0, 102.0, 42.0, fill=CARD)
    cells = ((116.0, 14.0, "blood oxygen", AZ, "spo2"), (165.0, 14.0, "breathing", TERRA, "air"),
             (116.0, 34.0, "brain waves", INK, "eeg"), (165.0, 34.0, "sleep stages", INK, "hyp"))
    for x0, y0, name, col, kind in cells:
        fig.label(name, x0, y0, pt=t.note_pt, va="top", color=col)
        tx, ty, tw, th = x0, y0 + 4.5, 44.0, 12.0
        if kind == "spo2":
            spo2_trace(fig, tx, ty, tw, th, dips=True, color=AZ, ref_at=0.55, label90=None, lw=0.5,
                       dip_centers=(0.15, 0.36, 0.58, 0.80))
        elif kind == "air":
            airflow(fig, tx + tw / 2, ty + th / 2, tw, th * 0.7, TERRA, lw=0.5)
        elif kind == "eeg":
            fig.trace("eeg", tx, ty, tw, th, color=GREY, lw=0.45)
        else:
            fig.hypnogram(tx, ty + 1.0, tw, th - 2.0, color=GREY, rem_color=AZ_T[2], lw=0.5)
    fig.label("197 measurements from that one night", 163.0, CAP_Y, anchor="middle", va="top",
              color=INK)
    fig.block_arrow(218.0, 27.0, 12.0, 12.0)

    number_circle(fig, 3, 237.0, 8.0)
    fig.place("calendar", 248.0, 28.0, w=16)
    fig.label("years later", 248.0, 38.5, pt=t.note_pt, anchor="middle", va="top", color=INK)
    organ_row_named(fig, 300.0, 27.0, pitch=22.0, icon_h=12.0, max_w=21.0)
    fig.label("Which measurement best predicts future disease?", 290.0, CAP_Y, anchor="middle",
              va="top", color=INK)
    return finish(fig, "band_fig1_top", OUT, final=final, attribution=None)


def fig1_closing(final):
    W, H = W_FIG1, 66.0
    fig = new(W, H)
    t = fig.theme
    fig.label("197 sleep measurements ranked by how well they predict disease", 8.0, 4.0,
              va="top", bold=True, color=INK)
    X0, CW = 8.0, 12.0
    rows = [(13.0, "1", AZ), (19.5, "2", AZ_T[2]), (26.0, "3", AZ_T[1]), (32.5, "4", AZ_T[0])]
    for y, rk, col in rows:
        fig.rect(X0, y - 2.4, CW, 4.8, rx=1.0, fill=col)
        fig.label(rk, X0 + CW / 2, y + 0.1, pt=t.note_pt, bold=True, anchor="middle", va="middle",
                  color=WHITE if col in (AZ, AZ_T[2]) else INK)
    top = "time below 90% (T90)"
    fig.label(top, X0 + CW + 3.0, 13.1, bold=True, color=AZ, va="middle")
    fig.label("rank 1 of 197", X0 + CW + 3.0 + text_width_mm(top, t.label_pt, True) + 3.0, 13.1,
              pt=t.note_pt, color=GREY, va="middle")
    fig.bracket(X0 + CW + 2.5, 17.0, 35.0, side="right", color=GREY, lw=0.45, tip=1.8)
    fig.label("the next ranks: other measures of oxygen during sleep", X0 + CW + 7.5, 26.1,
              pt=t.note_pt, color=AZ_D, va="middle", max_w=60.0, leading=1.15)
    for yy in (39.5, 42.5, 45.5):
        fig.circle(X0 + CW / 2, yy, 0.45, fill=GREY)
    fig.rect(X0, 52.0 - 2.4, CW, 4.8, rx=1.0, fill=GREEN)
    fig.label("192", X0 + CW / 2, 52.1, pt=t.note_pt, bold=True, anchor="middle", va="middle", color=WHITE)
    tst = "total sleep time"
    fig.label(tst, X0 + CW + 3.0, 52.1, bold=True, color=GREEN_D, va="middle")
    fig.label("rank 192 of 197", X0 + CW + 3.0 + text_width_mm(tst, t.label_pt, True) + 3.0, 52.1,
              pt=t.note_pt, color=GREY, va="middle")

    fig.line(118.0, 8.0, 118.0, 58.0, stroke=GREY_L, lw=0.3)
    TX, TY, TW, TH = 138.0, 12.0, 60.0, 24.0
    spo2_trace(fig, TX, TY, TW, TH, dips=True, color=AZ, ref_at=0.5, label90="left")
    night_axis(fig, TX, TY + TH + 1.5, TW)
    swatch_label(fig, TX + TW / 2, 46.5, "time below 90% (T90)", anchor="middle")
    fig.arrow((203.0, 24.0), (222.0, 24.0), color=AZ, lw=1.2)
    organ_row_named(fig, 270.0, 22.0, pitch=24.0, icon_h=12.0, max_w=23.0)
    fig.label("How low oxygen falls during sleep predicts future disease better than any other "
              "sleep measurement", 230.0, 54.0, pt=t.title_pt, bold=True, anchor="middle", va="top",
              color=INK, max_w=200.0)
    return finish(fig, "band_fig1_bottom", OUT, final=final, attribution=None)


def fig3_opening(final):
    """Two people whose sleep length and oxygen disagree: a short sleeper with
    preserved oxygenation, a long sleeper with low oxygenation. Which one matters?"""
    W, H = W_FIG3, 66.0
    fig = new(W, H)
    t = fig.theme
    fig.label("sleep length", 100.0, 3.0, anchor="middle", va="top", bold=True, color=GREEN_D)
    fig.label("oxygen during sleep", 187.0, 3.0, anchor="middle", va="top", bold=True, color=AZ)
    rows = ((22.0, "short sleeper", 20.0, True, False), (50.0, "long sleeper", 52.0, False, True))
    for y, name, bar_w, short, dips in rows:
        fig.place("sleeper-bed", 38.0, y, w=28)
        if short:
            fig.place("moon-green", 15.0, y - 5.0, w=8)
            fig.place("clock-green", 15.0, y + 5.0, w=7)
        else:
            fig.place("moon-green", 15.0, y, w=8)
        fig.label(name, 38.0, y + 10.5, pt=t.note_pt, anchor="middle", va="top", bold=True, color=GREEN_D)
        fig.rect(74.0, y - 3.0, bar_w, 6.0, rx=1.2, fill=GREEN)
        fig.label("a few hours" if short else "many hours", 74.0 + bar_w + 3.0, y + 0.1, pt=t.note_pt,
                  color=GREEN_D, va="middle")
        spo2_trace(fig, 162.0, y - 9.0, 50.0, 18.0, dips=dips, color=AZ if dips else AZ_T[2],
                   ref_at=0.5 if dips else 0.6, label90="left", dip_centers=(0.15, 0.36, 0.58, 0.80))
        fig.label("low oxygenation" if dips else "preserved oxygenation", 187.0, y + 10.5,
                  pt=t.note_pt, anchor="middle", va="top", bold=True, color=AZ)
    fig.line(150.0, 11.0, 150.0, 60.0, stroke=GREY_L, lw=0.3)
    fig.bracket(219.0, 11.0, 60.0, side="right", color=GREY, lw=0.5)
    fig.circle(230.0, 35.5, 4.2, fill=GREY)
    fig.label("?", 230.0, 35.7, pt=t.letter_pt, bold=True, anchor="middle", va="middle", color=WHITE)
    fig.arrow((235.5, 35.5), (246.0, 35.5), color=GREY, lw=0.6)
    organ_row_named(fig, 292.0, 22.0, pitch=22.0, icon_h=11.0, max_w=21.0)
    fig.label("Which one matters for future disease?", 292.0, 47.0, pt=t.title_pt, bold=True,
              anchor="middle", va="top", color=INK, max_w=88.0)
    return finish(fig, "band_fig3_top", OUT, final=final, attribution=None)


def fig3_closing(final):
    W, H = W_FIG3, 66.0
    fig = new(W, H)
    t = fig.theme
    for y, dips in ((17.0, False), (49.0, True)):
        fig.place("moon-green", 12.0, y - 5.0, w=7)
        fig.place("clock-green", 12.0, y + 4.5, w=6)
        fig.place("sleeper-bed", 34.0, y, w=28)
        fig.label("short sleep", 34.0, y + 10.5, pt=t.note_pt, anchor="middle", va="top", bold=True,
                  color=GREEN_D)
        spo2_trace(fig, 68.0, y - 9.0, 50.0, 18.0, dips=dips, color=AZ if dips else AZ_T[2],
                   ref_at=0.5 if dips else 0.6, label90="left", dip_centers=(0.15, 0.36, 0.58, 0.80))
        fig.label("low oxygenation" if dips else "preserved oxygenation", 93.0, y + 10.5,
                  pt=t.note_pt, anchor="middle", va="top", bold=True, color=AZ)
        fig.arrow((122.0, y), (136.0, y), color=AZ, lw=1.0)
        if not dips:
            answer(fig, 139.0, y, "no extra disease", sign="check", color=AZ, r=2.8, anchor="start")
        else:
            organ_row_named(fig, 174.0, y - 2.0, pitch=21.0, icon_h=10.5, max_w=20.0)
            fig.label("1.53 times the risk", 222.0, y - 1.9, bold=True, color=AZ, va="middle")
    fig.line(255.0, 8.0, 255.0, 58.0, stroke=GREY_L, lw=0.3)
    fig.label("Oxygen matters.", 296.0, 18.0, pt=t.title_pt, bold=True, anchor="middle", va="top",
              color=AZ)
    fig.label("Sleep length does not,", 296.0, 26.0, pt=t.title_pt, bold=True, anchor="middle",
              va="top", color=INK)
    fig.label("once oxygen is known.", 296.0, 32.0, pt=t.title_pt, bold=True, anchor="middle",
              va="top", color=INK)
    return finish(fig, "band_fig3_bottom", OUT, final=final, attribution=None)


def fig5_opening(final):
    W, H = W_FIG5, 66.0
    fig = new(W, H)
    t = fig.theme
    fig.place("cpap-mask", 30.0, 28.0, w=30)
    fig.label("sleep apnea treated with PAP", 30.0, 46.0, anchor="middle", va="top", bold=True,
              color=INK)
    fig.label("before PAP", 97.0, 6.0, pt=t.note_pt, anchor="middle", va="top", bold=True, color=GREY)
    spo2_trace(fig, 72.0, 12.0, 50.0, 20.0, dips=True, color=AZ, ref_at=0.5, label90="left")
    night_axis(fig, 72.0, 33.5, 50.0)
    swatch_label(fig, 97.0, 42.5, "time below 90%", anchor="middle")
    fig.block_arrow(128.0, 16.0, 14.0, 12.0)
    fig.label("starts PAP", 135.0, 30.0, pt=t.note_pt, anchor="middle", va="top", color=INK)
    fig.label("on PAP", 180.0, 2.0, pt=t.note_pt, anchor="middle", va="top", bold=True, color=GREY)
    spo2_trace(fig, 155.0, 8.0, 50.0, 16.0, dips=False, color=AZ_T[2], ref_at=0.6, label90="left")
    fig.label("oxygen improved", 208.0, 16.1, pt=t.note_pt, va="middle", color=AZ)
    spo2_trace(fig, 155.0, 30.0, 50.0, 16.0, dips=True, color=AZ, ref_at=0.5, label90="left",
               dip_centers=(0.15, 0.36, 0.58, 0.80))
    fig.label("oxygen still low", 208.0, 38.1, pt=t.note_pt, va="middle", color=AZ)
    fig.bracket(233.0, 8.0, 46.0, side="right", color=GREY, lw=0.5)
    fig.circle(243.5, 27.0, 4.2, fill=GREY)
    fig.label("?", 243.5, 27.2, pt=t.letter_pt, bold=True, anchor="middle", va="middle", color=WHITE)
    fig.arrow((249.0, 27.0), (258.0, 27.0), color=GREY, lw=0.6)
    organ_row_named(fig, 296.0, 22.0, pitch=21.0, icon_h=11.0, max_w=20.0)
    fig.label("Does disease risk change when oxygen improves?", 296.0, 47.0, pt=t.title_pt,
              bold=True, anchor="middle", va="top", color=INK, max_w=92.0)
    return finish(fig, "band_fig5_top", OUT, final=final, attribution=None)


# ==================================================================== C. home recording (496.8 x 79 pt)

def home_band(final):
    W, H = W_HOME, H_HOME                                  # 175.3 x 27.9 mm
    fig = new(W, H)
    t = fig.theme
    # two houses, a sleeper in each, cohort and size under
    for cx, name, n in ((11.5, "SHHS", "5,802 people"), (33.5, "MrOS", "2,911 people")):
        house(fig, cx, 10.0, 20.0, 17.0)
        fig.place("sleeper-bed", cx, 13.4, w=12)
        fig.label(name, cx, 19.4, pt=t.note_pt, anchor="middle", va="top", bold=True, color=INK)
        fig.label(n, cx, 23.0, pt=t.note_pt, anchor="middle", va="top", color=INK)
    # the home oximeter
    fig.place("oximeter-azure", 51.5, 6.8, w=10)
    fig.label("oxygen sensor\nat home", 51.5, 12.0, pt=t.note_pt, anchor="middle", va="top",
              color=AZ_D, leading=1.15)
    # the signal it records
    TX, TY, TW, TH = 64.5, 2.0, 19.0, 12.0
    spo2_trace(fig, TX, TY, TW, TH, dips=True, color=AZ, ref_at=0.5, label90="left",
               dip_centers=(0.15, 0.36, 0.58, 0.80))
    night_axis(fig, TX, TY + TH + 1.5, TW, pt=t.floor_pt)
    swatch_label(fig, TX + TW / 2, 23.6, "time below 90%", anchor="middle")
    # about 12 years later, the same diseases
    fig.arrow((86.0, 8.0), (112.0, 8.0), color=AZ, lw=1.0)
    fig.label("about 12 years later", 97.5, 10.2, pt=t.note_pt, anchor="middle", va="top", color=INK)
    _, bottom = organ_row_named(fig, 143.75, 4.6, pitch=16.5, icon_h=6.5, max_w=16.0, name_dy=1.5)
    fig.label("The same oxygen signal predicts disease at home", 143.75, bottom + 1.6, bold=True,
              anchor="middle", va="top", color=INK, max_w=66.0, leading=1.12)
    return finish(fig, "band_supp18_top", OUT, final=final, attribution=None)


# ==================================================================== drivers

SHEETS = {
    "iteration_1": lambda f: iter1(f),
    "iteration_2": iter2,
    "iteration_3": iter3,
    "iteration_4": lambda f: iter1(f, before_bed=True, stem="iteration_4"),
    "iteration_5": lambda f: iter1(f, before_bed=True, rank_badge=True, stem="iteration_5"),
    "band_fig1_top": fig1_opening, "band_fig1_bottom": fig1_closing,
    "band_fig3_top": fig3_opening, "band_fig3_bottom": fig3_closing,
    "band_fig5_top": fig5_opening, "band_fig5_bottom": fig5_closing,
    "band_supp18_top": home_band,
}

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--final", action="store_true")
    ap.add_argument("--only", default="")
    args = ap.parse_args()
    records = []
    for stem, fn in SHEETS.items():
        if args.only and not stem.startswith(args.only):
            continue
        records.append(fn(args.final))
    if not args.only:
        (ROOT / "_build" / "records.json").write_text(json.dumps(records, indent=1))
