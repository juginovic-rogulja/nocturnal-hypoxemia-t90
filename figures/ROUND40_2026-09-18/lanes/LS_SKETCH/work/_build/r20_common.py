"""Round 20, lane LB: the pieces the round-20 rebuild adds on top of the round-18 toolkit.

Everything here exists because of a specific instruction from Alen in this round.

  ONE_LINE        the single horizontal line Figure 1's three parts are aligned on. The bed,
                  both grey block arrows, the middle card and the organ row all centre on it,
                  which is what "the bed must sit in line with the grey arrow" asks for once
                  it is applied to the whole band instead of to one element.
  sensor_callout  an icon plus its name, used for the pulse oximeter added to Figure 1 part 1
  rank_device     the sorted-bar mark that says T90 came first out of 197 measurements
  band_pill       a filled pill carrying a band label, with a short rule that lands on the
                  arrow it belongs to, so the label is anchored instead of floating
  severity_ladder the apnea severity wedge redrawn in the round-20 charcoal ladder, light to
                  dark, so the severity range is carried by the family's own four rungs
  ROAD_BADGE_X    ONE constant for both of Figure 6's cross badges, so they cannot drift

COLOUR LAW, round 20: sleep apnea and breathing events are CHARCOAL #37474f on the ladder
#d3d7d8 #a5acb0 #737e84 #37474f with the dark text twin #29353b. Oxygen azure #0187d1, sleep
duration green #39c445 with the ladder #c8eecb #92df99 #5dcf66 #298d32, ink #1a1d20, greys
#8a9099 and #ccd1d6. Arial, 9 pt floor. No em-dashes, no semicolons, and the words about
groups and about oxygen that Alen banned are never printed.
"""
from __future__ import annotations

from l1_common import APNEA_D, APNEA_L, AZ, INK, text_width_mm

# ---------------------------------------------------------------- Figure 1's one line
# The middle card of Figure 1 top runs y 12 .. 58, so its centre is 35. Every part of the
# band is aligned on that number: the sleeping person, the two grey block arrows and the
# organ row. Alen's ruling was written about the bed and the first arrow, and applying it to
# the whole band is what makes the three parts read as one sentence.
ONE_LINE = 35.0
ARROW_W, ARROW_H = 12.0, 12.0                     # the grey block arrow, part 1 to part 2

# ---------------------------------------------------------------- Figure 6's two badges
ROAD_BADGE_X = 176.0                              # both cross badges, one constant, no drift


def block_arrow_on_line(fig, x, *, w=ARROW_W, h=ARROW_H, y=ONE_LINE):
    """The grey block arrow of Figure 1, centred on the band's one line. Both transitions
    call this, so part 2 to part 3 is the same arrow as part 1 to part 2 (Alen's ruling)."""
    return fig.block_arrow(x, y - h / 2, w, h)


def sensor_callout(fig, icon, cx, cy, w, text, *, pt=None, gap=1.6, color=INK, max_w=None):
    """A sensor drawn once and named under it, centred on cx. Used for the pulse oximeter
    Alen asked to add to Figure 1 part 1. Returns (placed, label)."""
    p = fig.place(icon, cx, cy, w=w)
    lb = fig.label(text, cx, p.y + p.h + gap, pt=pt or fig.theme.note_pt, anchor="middle",
                   va="top", color=color, max_w=max_w, leading=1.15)
    return p, lb


def rank_device(fig, cx, cy, text, *, w=5.2, bar_h=1.0, gap_y=1.55, pt=None, color=AZ,
                other=None, text_color=AZ, gap=2.0, bars=(1.0, 0.68, 0.46)):
    """Three sorted bars with the longest on top in azure, then the claim beside them, the
    pair centred on cx. Round 20: Figure 1's bottom band has to say that the time below 90%
    came first out of 197 measurements, and Alen asked for it briefly or as a small visual
    device. A sorted-bar mark says "top of a ranking" without printing a digit, which a
    numbered circle would, and the numbered circles on this figure already mean step 1, 2, 3.
    """
    from l1_common import GREY_L
    pt = pt or fig.theme.note_pt
    other = other or GREY_L
    tw = text_width_mm(text, pt, True)
    total = w + gap + tw
    x = cx - total / 2
    y0 = cy - gap_y - bar_h / 2
    for i, frac in enumerate(bars):
        fig.rect(x, y0 + i * gap_y, w * frac, bar_h, rx=bar_h / 2,
                 fill=color if i == 0 else other)
    lb = fig.label(text, x + w + gap, cy, pt=pt, bold=True, color=text_color, va="middle")
    return (x, x + total, lb)


def band_pill(fig, cx, cy, text, *, fill, text_color, to_arrow=None, rule_color=None,
              pt=None, pad_x=3.0, pad_y=1.4, rule_lw=0.9, bold=True):
    """A band label set in a filled pill, with a short rule running from the pill's right
    edge onto the arrow it labels. Round 20: on Figure 6 the two band labels used to float
    with nothing to anchor them, and Alen asked for emphasis plus an anchor to the arrow.

    ``to_arrow`` is the x on the arrow, at this same y, that the rule should land on.
    Returns (x_left, x_right) of the pill."""
    pt = pt or fig.theme.label_pt
    th = fig.theme
    tw = text_width_mm(text, pt, bold)
    h = pt * 0.352778 * 1.0 + 2 * pad_y                     # cap-and-descender box plus pad
    w = tw + 2 * pad_x
    x, y = cx - w / 2, cy - h / 2
    if to_arrow is not None:
        # the rule is drawn in the ARROW's colour and runs onto its centre line, so the pill
        # reads as a label ON that arrow rather than as a mark near it
        fig.line(x + w, cy, to_arrow, cy, stroke=rule_color or fill, lw=rule_lw)
    fig.rect(x, y, w, h, rx=h / 2, fill=fill)
    fig.label(text, cx, cy, pt=pt, bold=bold, color=text_color, anchor="middle",
              va="middle")
    return (x, x + w)


def severity_ladder(fig, x, y, w, *, left_text, right_text, ladder=APNEA_L,
                    text_color=APNEA_D, pt=None, tri_h=2.2, steps=4):
    """The apnea severity wedge, redrawn in the round-20 charcoal ladder. The wedge widens
    from no apnea on the left to severe on the right AND darkens through the family's four
    rungs, so the severity range is carried by the colour as well as by the shape. The two
    end labels are set in the dark text twin, never in a rung."""
    t = fig.theme
    pt = pt or t.note_pt
    seg = w / steps
    for i in range(steps):
        x0, x1 = x + i * seg, x + (i + 1) * seg
        h0 = 0.35 + (tri_h - 0.35) * (i / steps)
        h1 = 0.35 + (tri_h - 0.35) * ((i + 1) / steps)
        fig.polygon([(x0, y - h0), (x1, y - h1), (x1, y + h1), (x0, y + h0)],
                    fill=ladder[i])
    fig.polygon([(x + w, y - tri_h - 1.1), (x + w + 3.0, y), (x + w, y + tri_h + 1.1)],
                fill=ladder[-1])
    a = fig.label(left_text, x, y + tri_h + 2.0, pt=pt, color=text_color, va="top",
                  bold=True)
    b = fig.label(right_text, x + w + 3.0, y + tri_h + 2.0, pt=pt, color=text_color,
                  anchor="end", va="top", bold=True)
    return a, b
