#!/usr/bin/env python3
"""Round 20, lane LB: the three schematic bands for Figures 1 and 3.

  band_fig1_top      968.66 pt   the study band
  band_fig1_bottom   968.66 pt   the concluding fork of Figure 1
  band_fig3_bottom   952.73 pt   the Figure 3 sketch

WHAT ROUND 20 CHANGES, and nothing else.

FIG 1 TOP
 1  ONE LINE, y = 35. The middle card runs 12 to 58, so its centre is 35. The sleeping
    person, both grey block arrows and the organ row are now all centred on that number,
    which is what "the bed must sit in line with the grey arrow" asks for once it is applied
    to the whole band. The bed rises about 5 mm and the brain-wave callout drops about 4 mm,
    from a label centred on 14.6 to one centred on 19.0.
 2  "brain-wave sensors" becomes "electroencephalography (EEG)".
 3  pulse oximetry is ADDED, right of the bed, on the same line, named under the icon. The
    probe faces left so the finger runs back towards the sleeper.
 4  part 3 gets the SAME grey block arrow as the part 1 to part 2 transition, in place of the
    thin ink arrow. The organ row is centred on the middle card. "50 diseases and death"
    moves ABOVE the organs. The closing question names the count.

FIG 1 BOTTOM
 5  the band now says that T90 came first out of 197 measurements: the header over the trace
    carries the abbreviation, and a small rank ladder under it carries the claim.

FIG 3 BOTTOM
 6  the title line is deleted. Nothing else moves, so the band still drops into its parent
    figure at exactly 952.73 pt and exactly 70 mm.

COLOUR, round 20: sleep apnea and breathing events are CHARCOAL #37474f on the ladder
#d3d7d8 #a5acb0 #737e84 #37474f with the dark text twin #29353b, replacing round 18's orange.
Oxygen azure #0187d1, sleep duration green #39c445 with the dark twin #298d32.

    python3 build_r20_bands.py            # 150 dpi previews into _build/previews/
    python3 build_r20_bands.py --final    # the three bands into the lane root
"""
from __future__ import annotations

import argparse
import json

from l1_common import (APNEA, APNEA_D, AZ, AZ_T, CARD, GREEN_D, GREY, GREY_D, INK, MM_PER_PT, ROOT,
                       airflow, finish, new, numbered_caption, spo2_trace, text_width_mm)
from r14_common import emg_trace, trace_pair
from r18_common import (HIGH_1PCT, LOW_1PCT, ORGANS_BAD, ORGANS_OK, SCOPE, eeg_on_head,
                        fork, organ_names, organ_row)
from r20_common import ONE_LINE, block_arrow_on_line, rank_device, sensor_callout

OUT = ROOT
CREDIT = "Icons: Servier Medical Art (CC BY 3.0)."

W1_PT = 968.66                                   # Figure 1 page width
W3_PT = 952.73                                   # Figure 3 page width
H1_TOP, H1_BOT, H3_BOT = 66.0, 72.0, 70.0

CAP1 = "one night of sleep"
CAP2 = "197 measurements"                        # Main_Fig1: "Rank among the 197 measurements"
CAP3 = "Which measurement best predicts 50 future diseases?"
EEG_LABEL = "electroencephalography (EEG)"
OXI_LABEL = "pulse oximetry"

# Figure 1's bottom band, round 20: the header names the abbreviation and the ladder carries
# the rank. 197 is the number Main_Fig1 prints on its own axis, and T90 came first.
T90_HEADER_1B = "time with oxygen below 90% (T90)"
RANK_CLAIM = "the most predictive of 197 measurements"

MSG1 = "Of everything measured overnight, oxygen is what drives long-term outcomes."
DURATIONS = "short, normal or long sleep"

# the sleeping person's own head, in the TIGHT icon's units. The artwork is unchanged, so
# the head is still the circle at (30, 34.5) r 8 of the original 122 x 82 drawing, which in
# the cropped 4 14 114 62 box sits at (26, 20.5) r 8.
_TIGHT_HEAD = (26.0, 20.5, 8.0)
_TIGHT_W = 114.0


def tight_head_frame(placed):
    """(cx, cy, r) of the sleeper's head, in mm, for the tight-box icon."""
    k = placed.w / _TIGHT_W
    return (placed.x + _TIGHT_HEAD[0] * k, placed.y + _TIGHT_HEAD[1] * k, _TIGHT_HEAD[2] * k)


def credit(fig, y=None):
    fig.label(CREDIT, fig.w - 3.0, fig.h - 4.6 if y is None else y, pt=fig.theme.floor_pt,
              anchor="end", va="top", color=GREY)


# ================================================================ Figure 1, top band

def fig1_top(final, *, breathing=APNEA, stem="band_fig1_top"):
    W, H = W1_PT * MM_PER_PT, H1_TOP
    fig = new(W, H)
    t = fig.theme
    CAL_Y = 8.0                                  # the line every numbered caption sits on

    # (x, y, name, TEXT colour, TRACE colour, kind). The apnea family sets its type in the
    # dark twin #29353b and draws its marks in #37474f, the way the green family already
    # sets type in #298d32 and draws in #39c445.
    #
    # ROUND 20 consequence, worth Alen's eye: charcoal #37474f and ink #1a1d20 sit at a
    # contrast ratio of 1.77 to 1 against each other, so a charcoal breathing trace beside
    # an INK heart-rhythm trace would have looked like the same colour and the card would
    # have taught the reader nothing about the family. The heart-rhythm TRACE therefore
    # joins the other three unfamilied channels in grey, which leaves exactly two coloured
    # channels in the card, azure for oxygen and charcoal for breathing. Round 18 did not
    # need this because its apnea colour was orange. Every LABEL keeps its ink, because
    # grey #8a9099 at 9.5 pt does not clear 4.5 to 1 on white.
    cells = ((107.0, 14.0, "blood oxygen", AZ, AZ, "spo2"),
             (139.5, 14.0, "breathing", APNEA_D, breathing, "air"),
             (172.0, 14.0, "heart rhythm (ECG)", INK, GREY, "ecg"),
             (107.0, 35.0, "brain waves", INK, GREY, "eeg"),
             (139.5, 35.0, "sleep structure", INK, GREY, "hyp"),
             (172.0, 35.0, "leg movements", INK, GREY, "emg"))

    # ---------------------------------------------------------------- part 1
    # ONE picture of the recording happening: the person in bed with the electrodes on their
    # own head, and the pulse oximeter beside them. Everything sits on ONE_LINE.
    numbered_caption(fig, 1, 2.7, CAL_Y, CAP1)
    bed = fig.place("sleeper-faceB-tight", 30.0, ONE_LINE, w=50.0)
    hx, hy, hr = tight_head_frame(bed)
    gx, gy = eeg_on_head(fig, hx, hy, hr, gather=(hx - hr * 0.95, hy - hr * 2.4))
    lx, ly = 19.5, 19.0                          # the EEG name, clear above the bed's box
    fig.line(gx + 0.6, gy - 0.5, lx - 1.2, ly + 0.6, stroke=GREY_D, lw=0.32)
    fig.label(EEG_LABEL, lx, ly, pt=t.note_pt, va="middle", color=INK)
    sensor_callout(fig, "oximeter-azure-l", 70.5, ONE_LINE, 24.0, OXI_LABEL, color=INK)
    block_arrow_on_line(fig, 88.0)

    # ---------------------------------------------------------------- part 2 (untouched)
    numbered_caption(fig, 2, 104.7, CAL_Y, CAP2)
    fig.card(104.0, 12.0, 100.0, 46.0, fill=CARD)
    TW, TH = 28.0, 11.0
    for x0, y0, name, tcol, mcol, kind in cells:
        fig.label(name, x0, y0, pt=t.note_pt, va="top", color=tcol)
        tx, ty = x0, y0 + 5.0
        if kind == "spo2":
            spo2_trace(fig, tx, ty, TW, TH, dips=True, color=mcol, ref_at=0.55,
                       label90=None, lw=0.5, dip_centers=(0.15, 0.36, 0.58, 0.80))
        elif kind == "air":
            airflow(fig, tx + TW / 2, ty + TH / 2, TW, TH * 0.7, mcol, lw=0.5)
        elif kind == "ecg":
            fig.trace("ecg", tx, ty, TW, TH, color=mcol, lw=0.5)
        elif kind == "eeg":
            fig.trace("eeg", tx, ty, TW, TH, color=mcol, lw=0.45)
        elif kind == "emg":
            emg_trace(fig, tx, ty, TW, TH, color=mcol, lw=0.45)
        else:
            fig.hypnogram(tx, ty + 1.0, TW, TH - 2.0, color=mcol, rem_color=GREY_D, lw=0.5)

    # ---------------------------------------------------------------- part 3
    PITCH, X0, SCALE = 19.5, 247.5, 0.85
    ICON_HALF = max(h for _n, ic in ORGANS_OK for _i, h, _d in ic) * SCALE / 2   # 5.74 mm
    GX = X0 + PITCH * 2
    numbered_caption(fig, 3, 208.7, CAL_Y, CAP3)
    block_arrow_on_line(fig, 217.6)
    fig.label("years", 223.6, ONE_LINE + 6.0 + 1.4, pt=t.note_pt, anchor="middle", va="top",
              color=GREY)
    # the outcome heading now sits ABOVE the organs (Alen's round-20 ruling)
    fig.label(SCOPE, GX, ONE_LINE - ICON_HALF - 2.2, bold=True, anchor="middle",
              va="bottom", color=GREY_D)
    organ_row(fig, GX, ONE_LINE, ORGANS_OK, pitch=PITCH, icon_scale=SCALE,
              names_y=ONE_LINE + ICON_HALF + 1.8, max_w=PITCH - 0.5)

    credit(fig)
    return finish(fig, stem, OUT, final=final, attribution=None)


# ================================================================ Figure 1, bottom band

def fig1_bottom(final, *, stem="band_fig1_bottom"):
    W, H = W1_PT * MM_PER_PT, H1_BOT
    fig = new(W, H)
    t = fig.theme

    Y_OK, Y_BAD = 23.5, 51.5                     # the two organ rows, round-19 spacing
    CY = (Y_OK + Y_BAD) / 2                      # the fork's stem, level with the recording
    TX, TW, TH = 16.0, 72.0, 26.0                # the oxygen measurement
    TY = CY - TH / 2
    FORK_X, ROW_X = TX + TW + 4.0, 162.0
    OCX, PITCH = 248.0, 37.0                     # organ block centre and pitch
    SCOPE_Y, NAME_GAP = 13.0, 1.8                # the outcome heading, names under the grey

    # ---- the message, top left, in his own words
    fig.label(MSG1, 6.0, 4.4, pt=t.title_pt, bold=True, color=INK, va="top")

    # ---- ROUND 20: the rank. Of the 197 measurements this night carries, the time below
    # 90% came first, so the header names the abbreviation and the ladder carries the claim.
    rank_device(fig, TX + TW / 2, 13.0, RANK_CLAIM, pt=t.note_pt)

    # ---- one header line: what is measured on the left, what is at stake on the right
    fig.label(T90_HEADER_1B, TX + TW / 2, TY - 5.5, bold=True, color=AZ, anchor="middle",
              va="middle")
    fig.label(SCOPE, OCX, SCOPE_Y, bold=True, anchor="middle", va="middle", color=GREY_D)

    # ---- the one thing measured: blood oxygen across the night
    spo2_trace(fig, TX, TY, TW, TH, dips=True, color=AZ, ref_at=0.5, label90="left", lw=0.6,
               dip_centers=(0.18, 0.49, 0.79))
    fig.line(TX, TY + TH + 1.4, TX + TW, TY + TH + 1.4, stroke=GREY, lw=0.3)
    fig.label("one night of sleep", TX + TW / 2, TY + TH + 2.6, pt=t.note_pt, color=GREY,
              anchor="middle", va="top")

    # ---- the fork, at ONE percent of the night
    corner = fork(fig, FORK_X, CY, ROW_X, Y_OK, Y_BAD)
    lab_x = (corner + ROW_X) / 2
    fig.label(LOW_1PCT, lab_x, Y_OK - 2.6, bold=True, color=AZ_T[2], anchor="middle",
              va="baseline")
    fig.label(HIGH_1PCT, lab_x, Y_BAD - 2.6, bold=True, color=AZ, anchor="middle",
              va="baseline")

    # ---- the two endings: the same five organ groups, normal colour and dark grey, with
    # the five names written ONCE at the bottom of the block, under the dark-grey row
    organ_row(fig, OCX, Y_OK, ORGANS_OK, pitch=PITCH)
    bad_bottom = organ_row(fig, OCX, Y_BAD, ORGANS_BAD, pitch=PITCH)[2]
    organ_names(fig, OCX, ORGANS_OK, bad_bottom + NAME_GAP, pitch=PITCH, color=INK)

    credit(fig)
    return finish(fig, stem, OUT, final=final, attribution=None)


# ================================================================ Figure 3, bottom band

def fig3_bottom(final, *, stem="band_fig3_bottom"):
    W, H = W3_PT * MM_PER_PT, H3_BOT
    fig = new(W, H)
    t = fig.theme

    Y_OK, Y_BAD = 23.5, 50.0                     # the two oxygen outcomes and their organs
    TX, TW, TH = 98.0, 44.0, 14.0
    Y_FLAT, Y_DIPS = Y_OK - TH / 2, Y_BAD - TH / 2
    LX = TX + TW + 3.2
    OCX, PITCH, NAME_GAP = 262.0, 28.5, 1.8
    T90_HEADER = "time with oxygen below 90%"

    # ROUND 20: the title line is deleted. Everything else keeps its round-19 position, so
    # the page rectangle is unchanged and the band still drops into its parent figure.

    # ---- one bed, one person, one label that covers all three durations
    bed = fig.place("sleeper-faceB", 36.0, 34.0, w=54)
    fig.label(DURATIONS, 36.0, bed.y + bed.h + 0.6, bold=True, anchor="middle", va="top",
              color=GREEN_D)

    # ---- for those sleepers, two oxygen outcomes
    fig.label(T90_HEADER, (TX - 4.0 + LX + text_width_mm(HIGH_1PCT, t.label_pt, True)) / 2,
              13.0, bold=True, color=AZ, anchor="middle", va="middle")
    for yt in (Y_FLAT, Y_DIPS):
        fig.arrow((62.0, 34.0), (TX - 13.0, yt + TH * 0.5), color=GREY, lw=0.5)
    trace_pair(fig, TX, Y_FLAT, Y_DIPS, TW, TH, label_x=LX, lw=0.6,
               low=LOW_1PCT, high=HIGH_1PCT)
    brx = LX + text_width_mm(HIGH_1PCT, t.label_pt, True) + 4.0

    # ---- the two endings
    fig.label(SCOPE, OCX, 13.0, bold=True, anchor="middle", va="middle", color=GREY_D)
    for yt in (Y_OK, Y_BAD):
        fig.arrow((brx, yt), (OCX - PITCH * 2 - 9.0, yt), color=GREY, lw=0.6)
    organ_row(fig, OCX, Y_OK, ORGANS_OK, pitch=PITCH, icon_scale=0.95)
    bad_bottom = organ_row(fig, OCX, Y_BAD, ORGANS_BAD, pitch=PITCH, icon_scale=0.95)[2]
    organ_names(fig, OCX, ORGANS_OK, bad_bottom + NAME_GAP, pitch=PITCH, color=INK)

    credit(fig)
    return finish(fig, stem, OUT, final=final, attribution=None)


SHEETS = {"band_fig1_top": fig1_top, "band_fig1_bottom": fig1_bottom,
          "band_fig3_bottom": fig3_bottom}

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--final", action="store_true")
    ap.add_argument("--only", default="")
    args = ap.parse_args()
    recs = [fn(args.final) for stem, fn in SHEETS.items()
            if not args.only or args.only in stem]
    (ROOT / "_build" / "records_r20.json").write_text(json.dumps(recs, indent=1))
    for r in recs:
        print(json.dumps({k: r[k] for k in ("stem", "w_pt", "h_mm", "audit", "min_pt")}))
