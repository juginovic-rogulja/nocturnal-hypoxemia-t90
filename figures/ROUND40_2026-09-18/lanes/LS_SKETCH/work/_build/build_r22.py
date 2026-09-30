#!/usr/bin/env python3
"""Round 22, lane LA: all SIX schematic sheets of the T90 paper.

ROUND 27 (2026-09-06), FIG 1 TOP ONLY, everything below this block is round 22 as it was.
 A  EVERY TILE NAMES HOW ITS SIGNAL IS RECORDED, in parentheses after the name, the way the
    heart tile already did: blood oxygen (pulse oximetry), breathing (airflow and belts),
    heart rhythm (ECG), brain waves (EEG), sleep structure (staging), leg movements (leg EMG).
    A label goes on two lines (name, then the technique) when one line would overrun its
    tile, and stays on one line when it fits. Two lines at the theme's 9.5 pt need 6.97 mm
    of room above the trace ink and at the 9 pt floor still 6.60 mm, while the room under
    the round-26 anchors is 5.56 to 7.29 mm, so the six labels are set at the 9 pt floor
    (LABEL_PT) and their anchors are lifted by LABEL_LIFT per row. The traces, the card and
    everything else keep their round-26 coordinates. See TILE_LABELS.
 B  THREE SCALP ELECTRODES ON THE FACE-B SLEEPER, no text: ink discs on the visible crown
    (the hair side, which points at the headboard), thin leads gathered into one bundle
    above the pillow that runs off toward the headboard. eeg_on_faceB() below. The round-22
    option device (eeg_on_head, for the old forward-facing head) is left as it was.

  band_fig1_top      968.66 pt      band_fig4_bottom   970.79 pt
  band_fig1_bottom   968.66 pt      band_fig5_bottom   968.94 pt
  band_fig3_bottom   952.73 pt      Main_Fig6          968.94 pt

Every width above is exact and none of them moves. Every page height is the height round 20
or round 21 delivered, so all six still drop into their parent figures without shifting one
number on the page.

WHAT ROUND 22 CHANGES, and nothing else.

ALL SIX SHEETS
 1  THE ICON CREDIT LINE IS DELETED. Every sheet printed "Icons: Servier Medical Art
    (CC BY 3.0)." in a bottom corner. It is gone from all six. Alen checked Servier's own
    licensing page: their images are CC BY 4.0, not 3.0, so the line was wrong as well as
    unwanted, and the correct attribution now lives in the manuscript's acknowledgements in
    Servier's required wording for adapted images. The artwork carries no credit line.

    The three round-20 and round-21 builders are imported and their own ``credit`` function
    is replaced by a no-op, so this round removes the line WITHOUT changing one other
    coordinate on Figures 3, 4, 5 and 6. Nothing else in those files is touched.

FIG 1 TOP
 2  PART 1 IS THE BED AND THE SLEEPING PERSON, and nothing else. The
    "electroencephalography (EEG)" label goes, its leader line goes, and the three scalp
    electrodes and their lead bundle go with them, because "just the bed and the sleeping
    person" is what part 1 is now for and an unlabelled electrode set is exactly the kind of
    mark Alen has to ask about. Part 1 is labelled by the numbered caption it already
    carries, "one night of sleep".

    The recorded signals are NOT hinted at again in part 1. Part 2 of the same band already
    draws blood oxygen, breathing, heart rhythm, brain waves, sleep structure and leg
    movements, six named channels inside one card, so the hint would repeat what the next
    picture states outright. One option that WOULD carry the hint at almost no cost, the
    electrodes kept but silent, is rendered into _build/previews/ and is not shipped.

 3  THE CLOSING QUESTION CARRIES THE DEATHS: part 3's caption becomes "Which measurement
    best predicts 50 future diseases and death?" and the separate "50 diseases and death"
    heading over the organ row is DELETED, because the question now states it and the sheet
    was saying it twice.

FIG 1 BOTTOM
 4  ROOM FOR THE PANEL LETTER. The page composer stamps a 13 pt bold letter f at the top
    left of this band and it came out sliced through the middle. The band now keeps a clear
    rectangle at its top left, LETTER_ZONE_W x LETTER_ZONE_H, with no ink of its own in it,
    and verify_r22.py measures that the rectangle is empty on the delivered PDF.
 5  THE RANK LINE IS SHORTER: "the most predictive of 197 measurements" becomes
    "most predictive of 197 measurements".
 6  THE SORTED-BAR MARK IS DELETED. Three little bars, longest and azure on top, then two
    shorter greys, sat to the left of that line. It was drawn to say "top of a ranking"
    without printing a digit. It did not earn its place: Alen had to ask what it was, and
    the words beside it already say the same thing in full.
 7  "50 diseases and death" MOVES FROM THE TOP OF THE ORGAN BLOCK TO THE BOTTOM, under the
    dark-grey damaged row and its five names, so it reads with the damaged outcome instead
    of heading the pair. The block is then re-balanced against a page that no longer has a
    credit line at the foot of it.

COLOUR LAW, unchanged. Sleep apnea and breathing events are CHARCOAL #37474f on the ladder
#d3d7d8 #a5acb0 #737e84 #37474f with the dark text twin #29353b. Oxygen azure #0187d1. Sleep
duration green #39c445 on #c8eecb #92df99 #5dcf66 #298d32. Ink #1a1d20. Arial, 9 pt floor.
Nothing prints an em-dash, a semicolon, the two words about groups, or the phrase about
oxygen that Alen banned.

    python3 build_r22.py             # 150 dpi previews into _build/previews/
    python3 build_r22.py --final     # SVG + 300 dpi PNG + PDF into the lane root
    python3 build_r22.py --options   # the not-shipped part-1 option and the layout options
"""
from __future__ import annotations

import argparse
import json

from l1_common import (APNEA, APNEA_D, AZ, AZ_T, CARD, GREY, GREY_D, INK, MM_PER_PT, ROOT,
                       airflow, finish, new, numbered_caption, spo2_trace)
from r14_common import emg_trace
from r18_common import (HIGH_1PCT, LOW_1PCT, ORGANS_BAD, ORGANS_OK, SCOPE, eeg_on_head,
                        fork, organ_names, organ_row)
from r20_common import ONE_LINE, block_arrow_on_line

import build_r20_bands as B20            # band_fig3_bottom
import build_r20_le as B45               # band_fig4_bottom, band_fig5_bottom
import build_r21_fig6 as B6              # Main_Fig6

OUT = ROOT

# ---------------------------------------------------------------- change 1, all six sheets
# The credit line is removed by neutralising each builder's own credit(), which is the whole
# of the change on Figures 3, 4, 5 and 6. Their geometry is not touched, so the only ink
# those three sheets lose is the line itself.
B20.credit = lambda fig, y=None: None
B45.credit = lambda fig: None
B6.credit = lambda fig: None

W1_PT = 968.66                                   # Figure 1 page width, both bands
H1_TOP, H1_BOT = 66.0, 72.0

CAP1 = "one night of sleep"
CAP2 = "197 measurements"
# ROUND 22: the deaths move INTO the question, and the heading over the organs goes.
CAP3 = "Which measurement best predicts 50 future diseases and death?"

T90_HEADER_1B = "time with oxygen below 90% (T90)"
RANK_CLAIM = "most predictive of 197 measurements"     # ROUND 22: was "the most predictive"

# ---------------------------------------------------------------- part 1 balance (round 21)
CARD_X = 104.0
BED_W = 50.0
SLEEPER_TIGHT = "sleeper-faceB-tight"           # ROUND 26: the Fig 1 sleeper, switch here
BED_CX = 30.0                                    # box 5.0 .. 55.0, the page's own margin
BED_INK_L, BED_INK_R = 0.42, 49.57               # the sleeper's ink inside its 50 mm box
ARROW_W = 12.0
EEG_DX = -10.5                                   # only the not-shipped option still uses it
EEG_LABEL = "electroencephalography (EEG)"       # only the not-shipped option still uses it


def arrow_centred(ink_right, next_edge, w=ARROW_W):
    """A grey block arrow sits at the CENTRE of the white gap between the two things it
    joins. Round 21's rule, and deleting the EEG callout does not move either edge it uses:
    the bed's ink is unchanged and the card is unchanged."""
    gap = (next_edge - ink_right - w) / 2.0
    return (ink_right + gap, gap)


ARROW_X, PART1_GAP = arrow_centred(BED_CX - BED_W / 2 + BED_INK_R, CARD_X)

# the sleeping person's own head in the tight icon's units, for the option render only
_TIGHT_HEAD = (26.0, 20.5, 8.0)
_TIGHT_W = 114.0


def tight_head_frame(placed):
    k = placed.w / _TIGHT_W
    return (placed.x + _TIGHT_HEAD[0] * k, placed.y + _TIGHT_HEAD[1] * k, _TIGHT_HEAD[2] * k)


# ---------------------------------------------------------------- change 4, the letter zone
# The composer stamps the panel letter in 13 pt Arial Bold with its left edge at 14.17 pt,
# which is 5.00 mm, and the glyph box is 14.51 pt tall, which is 5.12 mm. The band reserves a
# rectangle wider and taller than that at its own top left and keeps every mark out of it.
LETTER_ZONE_W = 12.0                             # mm, from the page's left edge
LETTER_ZONE_H = 8.0                              # mm, from the page's top edge

# ---------------------------------------------------------------- change 7, the bottom band
# Round 21 balanced this band against the credit line at its foot. The credit line is gone
# and "50 diseases and death" has moved to the foot of the organ block, so the balance is
# solved again, this time against the page itself: the white above the first ink equals the
# white below the last ink. BOTTOM_SHIFT is the lift that does it, measured by
# solve_bottom_shift() on real ink at 600 dpi and pinned here so a plain build is
# reproducible. It also has to leave LETTER_ZONE_H clear at the top, and it does.
BOTTOM_SHIFT = 4.72
SCOPE_GAP = 4.6                                  # names' last ink to the scope line's top


# ---------------------------------------------------------------- round 27, the tile labels
# Column pitch is 32.5 mm and the traces are 28.0 mm wide. A one-line label "fits" when its
# right edge stays inside the tile's own column with the card's 3 mm padding, which is where
# the round-26 heart label already ran to. The other four overrun and take two lines.
LABEL_PT = 9.0                                   # the 9 pt floor, see the docstring
LABEL_LEADING = 1.15                             # the band's own two-line leading (part 3)
LABEL_LIFT = (1.0, 2.0)                          # mm the anchors rise, row 1 and row 2
TILE_LABELS = ("blood oxygen\n(pulse oximetry)", "breathing\n(airflow and belts)",
               "heart rhythm (ECG)", "brain waves (EEG)", "sleep structure\n(staging)",
               "leg movements\n(leg EMG)")

# ---------------------------------------------------------------- round 27, the wires
WIRES = True                                     # the shipped band carries the electrodes
LEAD_LW = 1.0 * MM_PER_PT                        # 1.0 pt leads
DOT_R = 0.8 * MM_PER_PT                          # 1.6 pt discs (diameter)
DOT_HALO = 0.4 * MM_PER_PT                       # white ring outside each disc, so an ink
#                                                  disc reads on the dark hair (round 18 did
#                                                  the same with a white stroke)
WIRE_ANGLES = (125.0, 150.0, 175.0)              # on the hair crescent, y up, 180 = crown
WIRE_GATHER = (-1.70, -1.05)                     # the bundle's knot, in head radii from the centre
WIRE_EXIT = (-2.35, -1.45)                       # where the bundle stops, at the headboard
_FACEB_HEAD = (26.0, 20.5, 8.0)                  # face B keeps the neutral head's circle


def eeg_on_faceB(fig, cx, cy, r, *, angles=WIRE_ANGLES, gather=WIRE_GATHER, exit_at=WIRE_EXIT,
                 color=INK, lead=INK, lw=LEAD_LW, dot_r=DOT_R, halo=None):
    """ROUND 27. Three scalp electrodes on the face-B head, which is turned three-quarter
    toward the viewer with its crown at the headboard: the hair crescent runs from about 104
    to 232 degrees (y up), so the discs sit on it above the pillow line. Each lead leaves
    its disc along the scalp normal and bends into one knot above the pillow, and the bundle
    runs on from the knot toward the headboard. No text. Returns the device's ink box in mm
    (x0, y0, x1, y1) and the disc centres."""
    import math
    gx, gy = cx + gather[0] * r, cy + gather[1] * r
    ex_, ey_ = cx + exit_at[0] * r, cy + exit_at[1] * r
    xs, ys, discs = [gx, ex_], [gy, ey_], []
    for a in angles:
        t = math.radians(a)
        px, py = cx + r * math.cos(t), cy - r * math.sin(t)
        mx, my = cx + r * 1.5 * math.cos(t), cy - r * 1.5 * math.sin(t)
        fig.path(f"M{px:.2f} {py:.2f} Q{mx:.2f} {my:.2f} {gx:.2f} {gy:.2f}", stroke=lead, lw=lw)
        fig._strokes_mm.append(lw)
        xs += [px, mx]
        ys += [py, my]
        discs.append((px, py))
    fig.path(f"M{gx:.2f} {gy:.2f} L{ex_:.2f} {ey_:.2f}", stroke=lead, lw=lw)
    fig._strokes_mm.append(lw)
    fig.circle(gx, gy, lw * 0.9, fill=lead)                  # the knot
    halo = DOT_HALO if halo is None else halo
    for px, py in discs:                                     # discs last, over the leads
        if halo:
            fig.circle(px, py, dot_r + halo, fill="#ffffff")
        fig.circle(px, py, dot_r, fill=color)
    pad = max(dot_r + (halo or 0.0), lw) + 0.1
    box = (min(xs) - pad, min(ys) - pad, max(xs) + pad, max(ys) + pad)
    fig._eeg = [{"cx": round(px, 2), "cy": round(py, 2), "r": round(dot_r, 3),
                 "head": [round(cx, 2), round(cy, 2), round(r, 2)]} for px, py in discs]
    return box, discs, (gx, gy), (ex_, ey_)


# ================================================================ Figure 1, top band

def fig1_top(final, *, breathing=APNEA, stem="band_fig1_top", bed_cx=None, arrow_x=None,
             electrodes=False, eeg_label=False, label_pt=None, label_leading=None,
             label_lift=None, labels=None, wires=None, dot_r=None, halo=None, outdir=None):
    """ROUND 22: part 1 is the bed and the sleeping person. ``electrodes`` and ``eeg_label``
    exist only so the not-shipped round-22 option can be rendered from this same code.
    ROUND 27: the six tile labels carry the recording technique (see TILE_LABELS, LABEL_PT,
    LABEL_LIFT) and the face-B sleeper wears three silent scalp electrodes (``wires``).
    The keyword knobs default to the module constants and exist for candidate renders."""
    W, H = W1_PT * MM_PER_PT, H1_TOP
    fig = new(W, H)
    t = fig.theme
    CAL_Y = 8.0                                  # the line every numbered caption sits on
    bed_cx = BED_CX if bed_cx is None else bed_cx
    arrow_x = ARROW_X if arrow_x is None else arrow_x
    label_pt = LABEL_PT if label_pt is None else label_pt
    label_leading = LABEL_LEADING if label_leading is None else label_leading
    label_lift = LABEL_LIFT if label_lift is None else label_lift
    labels = TILE_LABELS if labels is None else labels
    wires = WIRES if wires is None else wires
    dot_r = DOT_R if dot_r is None else dot_r

    cells = ((107.0, 14.0, labels[0], AZ, AZ, "spo2"),
             (139.5, 14.0, labels[1], APNEA_D, breathing, "air"),
             (172.0, 14.0, labels[2], INK, GREY, "ecg"),
             (107.0, 35.0, labels[3], INK, GREY, "eeg"),
             (139.5, 35.0, labels[4], INK, GREY, "hyp"),
             (172.0, 35.0, labels[5], INK, GREY, "emg"))

    # ---------------------------------------------------------------- part 1
    numbered_caption(fig, 1, 2.7, CAL_Y, CAP1)
    bed = fig.place(SLEEPER_TIGHT, bed_cx, ONE_LINE, w=BED_W)
    wire_box, wire_discs = None, []
    if wires:                                    # ROUND 27, the shipped device
        k = bed.w / _TIGHT_W
        hx, hy, hr = (bed.x + _FACEB_HEAD[0] * k, bed.y + _FACEB_HEAD[1] * k,
                      _FACEB_HEAD[2] * k)
        wire_box, wire_discs, _knot, _exit = eeg_on_faceB(fig, hx, hy, hr, dot_r=dot_r, halo=halo)
    if electrodes:                               # the round-22 option, never shipped
        hx, hy, hr = tight_head_frame(bed)
        gx, gy = eeg_on_head(fig, hx, hy, hr, gather=(hx - hr * 0.95, hy - hr * 2.4))
        if eeg_label:
            lx, ly = bed_cx + EEG_DX, 19.0
            fig.line(gx + 0.6, gy - 0.5, lx - 1.2, ly + 0.6, stroke=GREY_D, lw=0.32)
            fig.label(EEG_LABEL, lx, ly, pt=t.note_pt, va="middle", color=INK)
    block_arrow_on_line(fig, arrow_x)

    # ---------------------------------------------------------------- part 2
    # ROUND 27: the labels change, the card and the traces do not.
    numbered_caption(fig, 2, 104.7, CAL_Y, CAP2)
    fig.card(104.0, 12.0, 100.0, 46.0, fill=CARD)
    TW, TH = 28.0, 11.0
    tile_labels = []
    for x0, y0, name, tcol, mcol, kind in cells:
        lift = label_lift[0] if y0 < 30.0 else label_lift[1]
        lb = fig.label(name, x0, y0 - lift, pt=label_pt, va="top", color=tcol,
                       leading=label_leading)
        tile_labels.append({"text": name, "pt": label_pt, "lines": name.count("\n") + 1,
                            "anchor": [x0, round(y0 - lift, 3)],
                            "box": [round(lb.x, 3), round(lb.y, 3), round(lb.x + lb.w, 3),
                                    round(lb.y + lb.h, 3)],
                            "trace_box": [x0, y0 + 5.0, x0 + TW, y0 + 5.0 + TH]})
        tx, ty = x0, y0 + 5.0                    # the trace keeps its round-26 place
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
    # ROUND 22: the question carries the deaths and the heading over the organs is gone.
    PITCH, X0, SCALE = 19.5, 247.5, 0.85
    ICON_HALF = max(h for _n, ic in ORGANS_OK for _i, h, _d in ic) * SCALE / 2   # 5.74 mm
    GX = X0 + PITCH * 2
    numbered_caption(fig, 3, 208.7, CAL_Y, CAP3)
    block_arrow_on_line(fig, 217.6)
    fig.label("years", 223.6, ONE_LINE + 6.0 + 1.4, pt=t.note_pt, anchor="middle", va="top",
              color=GREY)
    organ_row(fig, GX, ONE_LINE, ORGANS_OK, pitch=PITCH, icon_scale=SCALE,
              names_y=ONE_LINE + ICON_HALF + 1.8, max_w=PITCH - 0.5)

    rec = finish(fig, stem, outdir or OUT, final=final, attribution=None)
    rec["tile_labels"] = tile_labels
    rec["label_pt"] = label_pt
    rec["label_leading"] = label_leading
    rec["label_lift"] = list(label_lift)
    rec["wires"] = None if wire_box is None else {
        "box": [round(v, 3) for v in wire_box],
        "discs": [[round(a, 3), round(b, 3)] for a, b in wire_discs],
        "lead_lw_mm": round(LEAD_LW, 4), "dot_r_mm": round(dot_r, 4),
        "halo_mm": round(DOT_HALO if halo is None else halo, 4),
        "bed_box": [round(bed.x, 3), round(bed.y, 3), round(bed.w, 3), round(bed.h, 3)]}
    return rec


# ================================================================ Figure 1, bottom band

def fig1_bottom(final, *, stem="band_fig1_bottom", shift=None, scope_at="bottom"):
    """ROUND 22: no sorted-bar mark, a shorter rank line, the outcome scope moved to the
    foot of the organ block, and a clear panel-letter rectangle at the top left."""
    W, H = W1_PT * MM_PER_PT, H1_BOT
    fig = new(W, H)
    t = fig.theme
    d = BOTTOM_SHIFT if shift is None else shift

    Y_OK, Y_BAD = 23.5 - d, 51.5 - d             # the two organ rows
    CY = (Y_OK + Y_BAD) / 2                      # the fork's stem, level with the recording
    TX, TW, TH = 16.0, 72.0, 26.0                # the oxygen measurement
    TY = CY - TH / 2
    FORK_X, ROW_X = TX + TW + 4.0, 162.0
    OCX, PITCH = 248.0, 37.0                     # organ block centre and pitch
    NAME_GAP = 1.8
    RANK_Y = 13.0 - d                            # the rank line, over the trace

    # ---- the rank, in words only. ROUND 22 deletes the sorted-bar mark that used to sit to
    # the left of this line and shortens the line itself. The claim is centred on the same
    # axis as the header under it and as the trace under that, so the left column reads as
    # one stack instead of as a mark plus a caption.
    fig.label(RANK_CLAIM, TX + TW / 2, RANK_Y, pt=t.note_pt, bold=True, color=AZ,
              anchor="middle", va="middle")

    # ---- what is measured
    fig.label(T90_HEADER_1B, TX + TW / 2, TY - 5.5, bold=True, color=AZ, anchor="middle",
              va="middle")

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

    # ---- the two endings, and ROUND 22's scope line at the FOOT of the pair
    if scope_at == "top":                        # round 21's position, option render only
        fig.label(SCOPE, OCX, 13.0 - d, bold=True, anchor="middle", va="middle",
                  color=GREY_D)
    organ_row(fig, OCX, Y_OK, ORGANS_OK, pitch=PITCH)
    bad_bottom = organ_row(fig, OCX, Y_BAD, ORGANS_BAD, pitch=PITCH)[2]
    names_bottom = organ_names(fig, OCX, ORGANS_OK, bad_bottom + NAME_GAP, pitch=PITCH,
                               color=INK)
    scope_y = None
    if scope_at == "bottom":
        scope_y = names_bottom + SCOPE_GAP
        fig.label(SCOPE, OCX, scope_y, bold=True, anchor="middle", va="top", color=GREY_D)

    rec = finish(fig, stem, OUT, final=final, attribution=None)
    rec["scope_y"] = None if scope_y is None else round(scope_y, 3)
    rec["shift"] = round(d, 3)
    rec["letter_zone"] = [0.0, 0.0, LETTER_ZONE_W, LETTER_ZONE_H]
    return rec


# ================================================================ the six, in one place

def build_all(final):
    """Every sheet this lane ships. Figures 3, 4, 5 and 6 come straight out of the round-20
    and round-21 builders with their credit line neutralised and nothing else changed."""
    recs = [fig1_top(final), fig1_bottom(final),
            B20.fig3_bottom(final), B45.fig4_bottom(final), B45.fig5_bottom(final),
            B6.fig6(final)]
    return recs


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--final", action="store_true")
    ap.add_argument("--only", default="")
    ap.add_argument("--options", action="store_true")
    ap.add_argument("--shift", type=float, default=None)
    args = ap.parse_args()

    if args.options:
        # NOT SHIPPED. The one small device that would hint at the recording in part 1
        # without a word of new text: the three scalp electrodes and their lead bundle,
        # kept but silent. Rendered so Alen can look at it and say yes or no.
        fig1_top(False, stem="opt_fig1_top_electrodes_silent", electrodes=True)
        fig1_top(False, stem="opt_fig1_top_bed_only")
        # the bottom band with the scope heading left where round 21 had it, for comparison
        fig1_bottom(False, stem="opt_fig1_bottom_scope_on_top", scope_at="top")
        raise SystemExit

    if args.shift is not None:
        fig1_bottom(False, stem=f"probe_bottom_shift{args.shift:g}", shift=args.shift)
        raise SystemExit

    recs = [r for r in build_all(args.final)
            if not args.only or args.only in r["stem"]] if not args.only else \
        [r for r in build_all(args.final) if args.only in r["stem"]]
    (ROOT / "_build" / "records_r22.json").write_text(json.dumps(recs, indent=1))
    for r in recs:
        print(json.dumps({k: r[k] for k in ("stem", "w_pt", "h_mm", "audit", "min_pt")}))
