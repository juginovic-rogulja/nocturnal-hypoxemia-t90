"""Round 18, lane LE: the shared grammar of the BOTTOM sketches on Figures 3, 4 and 5, and
of the standalone Figure 6.

Alen settled the family on the Figure 3 bottom sketch and asked for Figures 4 and 5 to carry
the same visual sentence:

    ONE source on the left, covering the whole range of the exposure, so the picture never
    repeats a row per group
        ->  ONE oxygen fork, a flat trace against a dipping trace
        ->  TWO endings, normally coloured organs against dark organs

so the reader sees that the exposure on the left does not decide the ending. Oxygen does.

Everything here is drawn with the round-16 toolkit and the LOCAL icon overrides in
_build/icons and _build/icons_grey. FigureLab's own library is never edited.

Colour law, ROUND 20: oxygen azure #0187d1 family, sleep apnea CHARCOAL #37474f on the ladder
#d3d7d8 #a5acb0 #737e84 #37474f with the dark text twin #29353b, sleep duration green #39c445
with the dark twin #298d32, ink #1a1d20, greys #8a9099 and #ccd1d6. Arial, 9 pt floor. No
em-dashes, no semicolons, and the words about groups and about oxygen that Alen banned are
never printed. The colours themselves live in l1_common, so this module follows the law by
importing APNEA rather than by naming a hex.
"""
from __future__ import annotations

import math

from l1_common import (APNEA, APNEA_T, AZ, AZ_T, GREY, GREY_L, INK, WHITE, night_axis,
                       spo2_trace, text_width_mm)

# ---------------------------------------------------------------- the printed wording
LOW = "≤1% of the night"
HIGH = ">1% of the night"
T90_HEADER = "time with oxygen below 90%"
SCOPE = "50 diseases and death"
GOOD_VERDICT = "no added risk"
BAD_VERDICT = "higher risk"

# ---------------------------------------------------------------- the two organ rows
# The dark twins are the local grey icons built in round 13. The live ones are FigureLab's
# own full-colour originals, which is what the grey twins were derived from, so the two rows
# are the same five organ groups in the same order and only the colour differs.
# ROUND 40 (2026-09-18, lane LS_SKETCH): the BRAIN leaves every sketch (no brain condition has a positive association
# in the results, dementia and Parkinson's run inverse). FOUR groups in this order. ORGANS_LIVE_R36 keeps the five-group
# list for the probes that pin the old panel geometry (Figure 6).
ORGANS_LIVE_R36 = [("heart and vessels", [("heart", 13.5, 0.0)]),
                   ("lungs", [("sv-lung", 13.5, 0.0)]),
                   ("metabolism", [("pancreas", 11.5, 0.0)]),
                   ("kidney and liver", [("kidney", 11.5, -6.0), ("liver", 10.0, 5.0)]),
                   ("brain", [("sv-brain-2", 12.0, 0.0)])]
ORGANS_LIVE = [g for g in ORGANS_LIVE_R36 if g[0] != "brain"]
assert [n for n, _ in ORGANS_LIVE] == ["heart and vessels", "lungs", "metabolism", "kidney and liver"]
ORGANS_DARK_R36 = [(name, [(f"{icon}-grey", h, dx) for icon, h, dx in icons])
                   for name, icons in ORGANS_LIVE_R36]
ORGANS_DARK = [(name, [(f"{icon}-grey", h, dx) for icon, h, dx in icons])
               for name, icons in ORGANS_LIVE]


def organ_row(fig, cx, cy, *, dark: bool, pitch=25.0, scale=1.0, names=True,
              name_dy=None, name_color=INK, name_max_w=None, groups=None):
    """The organ groups on one baseline, centred on cx (four groups since round 40). dark=True draws the grey
    twins, dark=False draws the same organs in their own colours. groups= overrides the list (round 40: the
    probes pass the five-group lists to pin the old geometry). Returns the y where the
    names end (or the icon bottom when names are off)."""
    if groups is None:
        groups = ORGANS_DARK if dark else ORGANS_LIVE
    t = fig.theme
    x0 = cx - pitch * (len(groups) - 1) / 2
    dy = name_dy if name_dy is not None else \
        max(h for _, ic in ORGANS_LIVE for _, h, _ in ic) * scale / 2 + 1.8
    bottom = cy + dy
    for i, (name, icons) in enumerate(groups):
        x = x0 + pitch * i
        for icon, ih, dx in icons:
            fig.place(icon, x + dx * scale, cy, h=ih * scale)
        if names:
            # ROUND 21: name_max_w lets a caller narrow the wrap width. Figure 6 needs it
            # so the widest group name stays inside the grey rectangle behind the row,
            # which is the same two-line treatment Figure 1's top band already gives the
            # same five names. Default unchanged.
            lb = fig.label(name, x, cy + dy, pt=t.note_pt, anchor="middle", va="top",
                           color=name_color,
                           max_w=pitch - 2.0 if name_max_w is None else name_max_w,
                           leading=1.15)
            bottom = max(bottom, lb.y + lb.h)
    return bottom


def verdict(fig, cx, y, text, *, good: bool, r=3.0, pt=None, color=INK, badge_color=AZ):
    """The one-line reading of an organ row, a badge and its words, centred on cx.
    good=True gets the check, good=False gets the plus, because the risk is added."""
    pt = pt or fig.theme.title_pt
    tw = text_width_mm(text, pt, True)
    total = 2 * r + 2.4 + tw
    x = cx - total / 2
    fig.badge(x + r, y, r, sign="check" if good else "+", color=badge_color,
              glyph_color=WHITE)
    return fig.label(text, x + 2 * r + 2.4, y + 0.1, pt=pt, bold=True, color=color,
                     va="middle")


# ---------------------------------------------------------------- the oxygen fork

def oxygen_fork(fig, src, tx, y_top, y_bot, tw, th, *, label_x, lw=0.6,
                low=LOW, high=HIGH, head_back=13.0, arrow_color=GREY,
                dip_centers=(0.14, 0.34, 0.56, 0.78), sub_top=None, sub_bot=None,
                sub_dy=6.0):
    """Two grey arrows out of one source into the two overnight-oxygen traces, the flat one
    on top and the dipping one below, each with its band label to the right. sub_top and
    sub_bot print a second, smaller line under a band label (Figure 5 uses it for the
    numeric cut under 'oxygen improved')."""
    t = fig.theme
    for y in (y_top, y_bot):
        fig.arrow(src, (tx - head_back, y + th * 0.5), color=arrow_color, lw=lw)
    spo2_trace(fig, tx, y_top, tw, th, dips=False, color=AZ_T[2], ref_at=0.62,
               label90="left", lw=lw)
    spo2_trace(fig, tx, y_bot, tw, th, dips=True, color=AZ, ref_at=0.5, label90="left",
               dip_centers=dip_centers, lw=lw)
    a = fig.label(low, label_x, y_top + th * 0.5, bold=True, color=AZ_T[2], va="middle")
    b = fig.label(high, label_x, y_bot + th * 0.5, bold=True, color=AZ, va="middle")
    if sub_top:
        fig.label(sub_top, label_x, y_top + th * 0.5 + sub_dy, pt=t.note_pt, color=AZ_T[2],
                  va="middle")
    if sub_bot:
        fig.label(sub_bot, label_x, y_bot + th * 0.5 + sub_dy, pt=t.note_pt, color=AZ,
                  va="middle")
    return a, b


def t90_header(fig, cx, y, *, va="middle", color=AZ):
    return fig.label(T90_HEADER, cx, y, bold=True, color=color, anchor="middle", va=va)


# ---------------------------------------------------------------- the apnea severity ramp

def airflow_ramp(fig, x, cy, w, h, *, color=APNEA, lw=0.6, events=7, breath_period=0.042,
                 pause_max=0.66, n=900):
    """ONE breathing trace that carries the whole sleep apnea severity range in a single
    line: regular breathing at the left end, then pauses that grow longer and arrive more
    often towards the right end. This is what replaces the four repeated severity rows, so
    the picture never says the severity groups behave differently from one another.

    The trace is built event by event. Every event is the same width, and within it the
    last `pause` fraction is a flat pause while the rest is breathing at a FIXED rate, so
    the breaths stay the same size and only the still time grows."""
    ev = 1.0 / events
    pts = []
    for i in range(n + 1):
        t = i / n
        k = min(events - 1, int(t / ev))
        sev = k / (events - 1)
        pause = pause_max * sev
        local = (t - k * ev) / ev                      # 0 .. 1 within this event
        breathe = 1.0 - pause
        if local < pause:
            amp, ph, cycles = 0.0, 0.0, 1.0            # the apnea, a flat pause
        else:
            u = (local - pause) / breathe              # 0 .. 1 across the breathing part
            amp = min(1.0, u / 0.10, (1.0 - u) / 0.10) if k else 1.0
            cycles = max(1.0, round(ev * breathe / breath_period))
            ph = u
        pts.append((x + w * t, cy - (h / 2) * amp * math.sin(2 * math.pi * cycles * ph)))
    fig.path("M" + " L".join(f"{px:.2f} {py:.2f}" for px, py in pts), stroke=color, lw=lw)


def severity_axis(fig, x, y, w, *, left_text, right_text, color=APNEA, pt=None,
                  tri_h=2.2):
    """A widening orange wedge under the ramp trace, from no apnea on the left to severe on
    the right, with the figure's own group wording at the two ends."""
    t = fig.theme
    pt = pt or t.note_pt
    fig.polygon([(x, y - 0.35), (x + w, y - tri_h), (x + w, y + tri_h), (x, y + 0.35)],
                fill=color, opacity=0.55)
    fig.polygon([(x + w, y - tri_h - 1.1), (x + w + 3.0, y), (x + w, y + tri_h + 1.1)],
                fill=color, opacity=0.85)
    # bold: #ff8f00 is a light hue and these pages print at about half authored size, so the
    # apnea words need the stroke weight to survive (the same reason the family green sets
    # its type in GREEN_D rather than in #57cc60)
    a = fig.label(left_text, x, y + tri_h + 2.0, pt=pt, color=color, va="top", bold=True)
    b = fig.label(right_text, x + w + 3.0, y + tri_h + 2.0, pt=pt, color=color,
                  anchor="end", va="top", bold=True)
    return a, b


# ---------------------------------------------------------------- the crossed roads

def big_crossed_road(fig, p_from, p_to, *, badge_at, label, qualifier=None, label_dy=7.0,
                     label_va="top", color=GREY, r=4.8, pt=None, lw=1.0, dash="4.4 2.8",
                     max_w=None, qual_gap=5.2):
    """Round 18: the same grey dashed elbow road as round 16, drawn much bigger and bolder,
    because Alen wants the two 'this alone does not do it' messages to read instantly. The
    caption is bold and one type size up, and the cross badge is nearly twice the radius.

    `qualifier` prints a second, smaller line naming the significance standard, so the bold
    line stays instant while the claim stays exact. label_va='top' stacks it under the bold
    line, label_va='bottom' stacks it above."""
    t = fig.theme
    fig.arrow(p_from, p_to, color=color, lw=lw, dash=dash, elbow="h")
    fig.badge(badge_at, p_from[1], r, sign="×", color=color, glyph_color=WHITE)
    y = p_from[1] + label_dy
    if qualifier and label_va == "bottom":
        fig.label(qualifier, badge_at, y, pt=t.note_pt, anchor="middle", va="bottom",
                  color=color)
        y -= qual_gap
    lb = fig.label(label, badge_at, y, pt=pt or t.title_pt, bold=True, anchor="middle",
                   va=label_va, color=color, max_w=max_w, leading=1.2)
    if qualifier and label_va == "top":
        fig.label(qualifier, badge_at, lb.y + lb.h + 1.2, pt=t.note_pt, anchor="middle",
                  va="top", color=color)
    return lb
