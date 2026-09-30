"""Round 14 additions to the shared L1 schematic toolkit (FigureLab engine).

New here, all drawn with primitives so nothing depends on an external icon:
  finger_probe   a fingertip inside a pulse-oximeter clip, with the NAIL visible, the clip
                 gripping the tip, the hinge and the cable at the closed end (Alen: the
                 round-13 probe "looks like a band-aid" and had no nail)
  emg_trace      leg-movement EMG: a quiet baseline with two bursts
  airflow_sev    breathing airflow with a chosen number of pauses, so the four sleep-apnea
                 severity groups can be told apart at a glance
  trace_pair     the two oxygen traces one group splits into (shared by the Fig 3 and Fig 4
                 top bands so both use exactly the same wording and geometry)
  fork_to_card   the bracket, the question mark, the arrow and the organ-group card, with
                 the card centred vertically on that arrow (Alen's round-14 ruling 3)

Colour law and type sizes come from l1_common. No em-dashes, no semicolons, 9 pt floor.
"""
from __future__ import annotations

import math
import random

from l1_common import (AZ, AZ_D, AZ_T, CARD, GREY, GREY_D, INK, TEAL, TERRA, WHITE,
                       number_circle, numbered_caption, spo2_trace, text_width_mm)

SKIN, SKIN_EDGE, SKIN_SH = "#f1c9a5", "#c89b72", "#dfac82"
NAIL, NAIL_EDGE = "#fbe0d0", "#c08a63"
PROBE_UNITS = 38.0          # local frame of finger_probe: x runs -19 .. +19

# the five organ groups in dark grey, one row, as Figs 1, 3 and 5 already draw them
ORGAN_GROUPS = [("heart and vessels", [("heart-grey", 13.5, 0.0)]),
                ("lungs", [("sv-lung-grey", 13.5, 0.0)]),
                ("metabolism", [("pancreas-grey", 11.5, 0.0)]),
                ("kidney and liver", [("kidney-grey", 11.5, -6.0), ("liver-grey", 10.0, 5.0)]),
                ("brain", [("sv-brain-2-grey", 12.0, 0.0)])]
SCOPE = "50 diseases and death"
LOW_T90, HIGH_T90 = "≤1% of the night", ">10% of the night"
T90_HEADER = "time with oxygen below 90%"


# ---------------------------------------------------------------- the finger probe

def finger_probe(fig, cx, cy, *, w=27.0, fh=7.0, top_jaw=-3.0, bot_jaw=6.0,
                 nail=(0.4, 10.2), bulge=1.15, base_u=19.0):
    """Fingertip inside a pulse-oximeter clip, in profile. The finger enters from the
    right and its tip sits in the clip. The hinge and the cable are at the closed left
    end, as on a real clip probe. The lower housing runs the length of the device and
    the upper jaw is short, so the NAIL stays visible on the dorsum and the shape reads
    as a finger rather than as a dressing. ``w`` is the true drawn width, so the printed
    size can be checked directly (the Fig 1 band prints at about 0.53 of authored size)."""
    k = w / PROBE_UNITS

    def X(u):
        return cx + u * k

    def Y(u):
        return cy + u * k

    fh2, tip_u = fh / 2, -7.0
    # ---- the finger: rounded tip on the left, square cut on the right into the hand
    fig.path(f"M{X(base_u):.2f} {Y(-fh2):.2f} L{X(tip_u + fh2):.2f} {Y(-fh2):.2f} "
             f"A{fh2 * k:.2f} {fh2 * k:.2f} 0 0 0 {X(tip_u + fh2):.2f} {Y(fh2):.2f} "
             f"L{X(base_u):.2f} {Y(fh2):.2f} Z", fill=SKIN, stroke=SKIN_EDGE, lw=0.42)
    fig.path(f"M{X(tip_u + 3.0):.2f} {Y(fh2 - 1.0):.2f} L{X(base_u):.2f} {Y(fh2 - 1.0):.2f}",
             stroke=SKIN_SH, lw=0.9 * k, opacity=0.5)
    fig.path(f"M{X(base_u - 3.2):.2f} {Y(-fh2 + 0.9):.2f} Q{X(base_u - 4.1):.2f} {cy:.2f} "
             f"{X(base_u - 3.2):.2f} {Y(fh2 - 0.9):.2f}", stroke=SKIN_EDGE, lw=0.32, opacity=0.85)
    # ---- the clip: hinge block, a short upper jaw and the longer lower housing
    hx0, hx1, jx0, jt = -12.0, -8.2, -10.6, 3.2
    fig.rect(X(jx0), Y(-fh2 - jt + 0.6), (top_jaw - jx0) * k, jt * k, rx=1.2 * k,
             fill=AZ, stroke=AZ_D, lw=0.42)
    fig.rect(X(jx0), Y(fh2 - 0.6), (bot_jaw - jx0) * k, jt * k, rx=1.2 * k,
             fill=AZ, stroke=AZ_D, lw=0.42)
    fig.rect(X(hx0), Y(-fh2 - jt + 0.4), (hx1 - hx0) * k, (fh + 2 * jt - 0.8) * k, rx=1.5 * k,
             fill=AZ, stroke=AZ_D, lw=0.42)
    fig.rect(X(jx0 + 0.9), Y(-fh2 - jt + 1.4), 4.4 * k, 1.2 * k, rx=0.55 * k,
             fill=AZ_T[0], opacity=0.95)                       # sensor window on the top jaw
    fig.circle(X(hx0 + 1.9), cy, 0.9 * k, fill=WHITE, opacity=0.55)
    # ---- the nail, starting where the upper jaw ends
    n0, n1, ny = top_jaw + nail[0], top_jaw + nail[1], -fh2
    fig.path(f"M{X(n0):.2f} {Y(ny + 0.6):.2f} "
             f"C{X(n0 + 0.4):.2f} {Y(ny - bulge):.2f} {X(n1 - 0.6):.2f} {Y(ny - bulge - 0.1):.2f} "
             f"{X(n1):.2f} {Y(ny + 0.5):.2f} "
             f"C{X(n1 - 2.0):.2f} {Y(ny + 1.7):.2f} {X(n0 + 1.8):.2f} {Y(ny + 1.8):.2f} "
             f"{X(n0):.2f} {Y(ny + 0.6):.2f} Z", fill=NAIL, stroke=NAIL_EDGE, lw=0.36)
    fig.path(f"M{X(n1 - 2.4):.2f} {Y(ny + 0.2):.2f} Q{X(n1 - 1.2):.2f} {Y(ny + 1.3):.2f} "
             f"{X(n1 - 0.3):.2f} {Y(ny + 0.8):.2f}", stroke=NAIL_EDGE, lw=0.3, opacity=0.75)
    # ---- the cable out of the hinge
    fig.path(f"M{X(hx0):.2f} {Y(2.0):.2f} C{X(hx0 - 4.5):.2f} {Y(3.2):.2f} "
             f"{X(hx0 - 4.2):.2f} {Y(7.2):.2f} {X(-19.0):.2f} {Y(8.0):.2f}",
             stroke=GREY_D, lw=0.6)
    return (X(-19.0), Y(-fh2 - jt), X(base_u), Y(fh2 + jt + 1.0))


# ---------------------------------------------------------------- new traces

def emg_trace(fig, x, y, w, h, *, bursts=((0.24, 0.09), (0.60, 0.12)), color=GREY,
              lw=0.4, seed=11):
    """Leg-movement EMG: a quiet baseline with two bursts of activity."""
    rng = random.Random(seed)
    mid, n, pts = y + h / 2, 420, []
    for i in range(n + 1):
        t = i / n
        amp = 0.09
        for c, wd in bursts:
            amp = max(amp, math.exp(-((t - c) / (wd / 2.0)) ** 2))
        pts.append((x + w * t, mid - rng.uniform(-1.0, 1.0) * amp * h * 0.46))
    fig.path("M" + " L".join(f"{px:.2f} {py:.2f}" for px, py in pts), stroke=color, lw=lw)


def airflow_sev(fig, cx, cy, w, h, *, pauses=(), cycles=6.0, color=TEAL, lw=0.55):
    """Breathing airflow with the given pause windows (fractions of the trace). No pause
    means regular breathing, more and longer pauses mean more severe sleep apnea."""
    x0, n, pts = cx - w / 2, 320, []
    for i in range(n + 1):
        t = i / n
        amp = 1.0
        for a, b in pauses:
            if a <= t <= b:
                amp = 0.0
            else:
                d = min(abs(t - a), abs(t - b))
                if d < 0.045:
                    amp = min(amp, d / 0.045)
        pts.append((x0 + w * t, cy - (h / 2) * amp * math.sin(2 * math.pi * cycles * t)))
    fig.path("M" + " L".join(f"{px:.2f} {py:.2f}" for px, py in pts), stroke=color, lw=lw)


# ---------------------------------------------------------------- the shared fork

def trace_pair(fig, x, y_flat, y_dips, w, h, *, label_x, label_pt=None, lw=0.6,
               dip_centers=(0.14, 0.34, 0.56, 0.78), low=LOW_T90, high=HIGH_T90):
    """The two oxygen traces one group splits into: little time below 90%, and much time
    below 90% shaded. The band labels carry the figure's own T90 cuts. Round 18 lets the
    caller pass the cuts, because the two new fork bands split at ONE percent, not ten."""
    spo2_trace(fig, x, y_flat, w, h, dips=False, color=AZ_T[2], ref_at=0.62, label90="left",
               lw=lw, label_pt=label_pt)
    spo2_trace(fig, x, y_dips, w, h, dips=True, color=AZ, ref_at=0.5, label90="left",
               dip_centers=dip_centers, lw=lw, label_pt=label_pt)
    a = fig.label(low, label_x, y_flat + h * 0.5, bold=True, color=AZ_T[2], va="middle")
    b = fig.label(high, label_x, y_dips + h * 0.5, bold=True, color=AZ, va="middle")
    return a, b


def fork_to_card(fig, *, brx, bracket_y0, bracket_y1, card_x, card_w, arrow_y,
                 card_h=42.0, organ_dy=22.5, pitch=28.0, qr=4.2, icon_scale=1.0,
                 scope_number=None, scope_mode="beside", scope_dy=2.4, scope_color=INK):
    """Bracket over every trace, a question mark, one arrow, and the organ-group card
    centred vertically on that arrow (Alen's round-14 ruling 3). The card's internal
    layout is the one the round-13 Fig 3 band already used, so the Fig 3 and Fig 4 bands
    are twins. Returns the card box.

    Round 16: `scope_number` turns "50 diseases and death" into a numbered callout in the
    azure circle device. scope_mode='beside' puts the circle to the LEFT of the words and
    centres the pair on the card (Alen's Fig 3 ruling). scope_mode='above' puts the circle
    ABOVE the words, both centred on the card (his Fig 4 ruling)."""
    t = fig.theme
    fig.bracket(brx, bracket_y0, bracket_y1, side="right", color=GREY, lw=0.5)
    qx = brx + 9.5
    fig.circle(qx, arrow_y, qr, fill=GREY)
    fig.label("?", qx, arrow_y + 0.2, pt=t.letter_pt, bold=True, anchor="middle", va="middle",
              color=WHITE)
    fig.arrow((qx + qr + 1.3, arrow_y), (card_x - 1.5, arrow_y), color=GREY, lw=0.6)
    cy = arrow_y - card_h / 2
    cm = card_x + card_w / 2
    fig.card(card_x, cy, card_w, card_h, fill=CARD)
    if scope_number is None:
        fig.label(SCOPE, cm, cy + scope_dy, bold=True, anchor="middle", va="top",
                  color=scope_color)
    elif scope_mode == "beside":
        numbered_caption(fig, scope_number, cm, cy + scope_dy + 1.63, SCOPE,
                         anchor="middle", color=scope_color)
    else:
        number_circle(fig, scope_number, cm, cy + scope_dy - 5.5)
        fig.label(SCOPE, cm, cy + scope_dy, bold=True, anchor="middle", va="top",
                  color=scope_color)
    oy = cy + organ_dy
    x0 = cm - pitch * (len(ORGAN_GROUPS) - 1) / 2
    name_dy = max(h for _, ic in ORGAN_GROUPS for _, h, _ in ic) * icon_scale / 2 + 1.8
    for i, (name, icons) in enumerate(ORGAN_GROUPS):
        x = x0 + pitch * i
        for icon, ih, dx in icons:
            fig.place(icon, x + dx * icon_scale, oy, h=ih * icon_scale)
        fig.label(name, x, oy + name_dy, pt=t.note_pt, anchor="middle", va="top", color=INK,
                  max_w=pitch - 2.0, leading=1.15)
    return (card_x, cy, card_w, card_h)


def t90_header(fig, tx, tw, label_x, y=1.2, va="top"):
    """The one azure heading over the trace column, centred on traces plus band labels.
    Round 16: va='middle' lets it sit on the same line as the numbered callout row."""
    cx = (tx - 4.0 + label_x + text_width_mm(HIGH_T90, fig.theme.label_pt, True)) / 2
    return fig.label(T90_HEADER, cx, y, bold=True, color=AZ, anchor="middle", va=va)
