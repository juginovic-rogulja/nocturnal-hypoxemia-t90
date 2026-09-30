"""Round 18 additions to the shared L1 schematic toolkit (FigureLab engine).

New here:
  head_frame     where the sleeping person's own head is, in millimetres, for any placed
                 sleeper icon (the icon draws the head at (30, 34.5) r 8 of a 122 x 82 box)
  eeg_on_head    brain-wave electrodes and their leads drawn ON that head, so ONE picture
                 shows the recording happening instead of a person plus a separate wired
                 head (Alen's round-18 ruling for the Figure 1 band)
  organ_row      one row of the five organ groups, either in their NORMAL Servier colours
                 or as the dark-grey twins, so a healthy row and a damaged row can be set
                 one above the other and told apart at a glance
  organ_names    the five group names on their own line, so round 19 can print them once at
                 the BOTTOM of the block, under the dark-grey row, instead of between the
                 two rows (Alen's round-19 ruling)
  ORGANS_OK / ORGANS_BAD   those two rows' icon lists, same geometry, different colour

Round-18 cuts: the fork in both new bands splits at ONE percent of the night, not ten,
because Alen is explicit that even above 1% is already worse.

Colour law and type sizes come from l1_common. No em-dashes, no semicolons, 9 pt floor.
"""
from __future__ import annotations

from l1_common import GREY, GREY_D, INK, WHITE

# ---- the two cuts this round forks on. LOW is r14_common's wording, unchanged.
LOW_1PCT = "≤1% of the night"
HIGH_1PCT = ">1% of the night"
T90_HEADER = "time with oxygen below 90%"
SCOPE = "50 diseases and death"

# ---- the five organ groups, normal colour and dark-grey twin, identical geometry
#      (the grey twins are recolours of exactly these files, so the shapes match)
# ROUND 40 (2026-09-18, lane LS_SKETCH): the BRAIN leaves every sketch. Alen: no brain condition has a positive
# association in the results (dementia and Parkinson's run inverse). FOUR groups, in this order, top and bottom rows
# identical in count and geometry. The five-group list is kept as ORGANS_OK_R36 for the probes that pin the old geometry.
ORGANS_OK_R36 = [("heart and vessels", [("heart", 13.0, 0.0)]),
                 ("lungs", [("sv-lung", 13.0, 0.0)]),
                 ("metabolism", [("pancreas", 11.0, 0.0)]),
                 ("kidney and liver", [("kidney", 11.5, -5.4), ("liver", 9.8, 4.4)]),
                 ("brain", [("sv-brain-2", 12.0, 0.0)])]
ORGANS_OK = [g for g in ORGANS_OK_R36 if g[0] != "brain"]
assert [n for n, _ in ORGANS_OK] == ["heart and vessels", "lungs", "metabolism", "kidney and liver"]
ORGANS_BAD_R36 = [(name, [(f"{ic}-grey", h, dx) for ic, h, dx in icons])
                  for name, icons in ORGANS_OK_R36]
ORGANS_BAD = [(name, [(f"{ic}-grey", h, dx) for ic, h, dx in icons])
              for name, icons in ORGANS_OK]

# the sleeper icons (sleeper-bed, sleeper-neutral) share one 122 x 82 viewBox and draw the
# head as a circle at (30, 34.5) with radius 8 in that box.
_HEAD_U = (30.0, 34.5, 8.0)
_VIEW_W = 122.0


def head_frame(placed):
    """(cx, cy, r) of the sleeping person's own head, in millimetres."""
    k = placed.w / _VIEW_W
    return (placed.x + _HEAD_U[0] * k, placed.y + _HEAD_U[1] * k, _HEAD_U[2] * k)


def eeg_on_head(fig, cx, cy, r, *, angles=(58.0, 100.0, 142.0), gather=None,
                color=INK, lead=GREY_D, lw=0.32, disc=0.23):
    """Brain-wave electrodes sitting on the scalp of a person already lying in bed, with
    their leads gathered into one bundle above the crown. ``angles`` are measured from the
    +x axis with y UP, so 90 is the top of the head and 142 is the back of it (the sleeper
    faces right). The bundle rises well clear of the crown so the leads read as cables
    running to a recorder, not as a hairnet. Returns the gather point, which is where a
    label's leader should land."""
    import math
    gx, gy = gather if gather else (cx - r * 0.95, cy - r * 3.2)
    for a in angles:
        t = math.radians(a)
        ex, ey = cx + r * math.cos(t), cy - r * math.sin(t)
        # the lead leaves the electrode along the scalp normal, then bends to the bundle
        mx, my = cx + r * 1.55 * math.cos(t), cy - r * 1.55 * math.sin(t)
        fig.path(f"M{ex:.2f} {ey:.2f} Q{mx:.2f} {my:.2f} {gx:.2f} {gy:.2f}",
                 stroke=lead, lw=lw)
    fig._eeg = []
    for a in angles:                                  # discs last, so no lead crosses one
        t = math.radians(a)
        ex, ey = cx + r * math.cos(t), cy - r * math.sin(t)
        fig.circle(ex, ey, r * disc, fill=color, stroke=WHITE, lw=lw * 0.7)
        fig._eeg.append({"cx": round(ex, 2), "cy": round(ey, 2), "r": round(r * disc, 2),
                         "head": [round(cx, 2), round(cy, 2), round(r, 2)]})
    fig.circle(gx, gy, r * 0.17, fill=lead)           # the bundle's collar
    return (gx, gy)


def organ_row(fig, cx, cy, groups, *, pitch, icon_scale=1.0, names_y=None,
              name_pt=None, name_color=INK, max_w=None):
    """One row of the five organ groups centred on cx. With ``names_y`` the group names are
    printed on that line instead of under the icons, so a single set of names can serve two
    rows stacked one above the other. Returns (x_first, x_last, bottom_y)."""
    x0 = cx - pitch * (len(groups) - 1) / 2
    bottom = cy
    for i, (name, icons) in enumerate(groups):
        x = x0 + pitch * i
        for icon, ih, dx in icons:
            p = fig.place(icon, x + dx * icon_scale, cy, h=ih * icon_scale)
            bottom = max(bottom, p.y + p.h)
        if names_y is not None:
            lb = fig.label(name, x, names_y, pt=name_pt or fig.theme.note_pt,
                           anchor="middle", va="top", color=name_color,
                           max_w=max_w or pitch - 1.5, leading=1.15)
            bottom = max(bottom, lb.y + lb.h)
    return (x0, x0 + pitch * (len(groups) - 1), bottom)


def organ_names(fig, cx, groups, y, *, pitch, pt=None, color=INK, max_w=None):
    """The five group names on one line centred on cx, written ONCE for a whole two-row
    organ block. Round 19 moves that line from between the two rows to the BOTTOM of the
    block, under the dark-grey row, so the names read as the caption of the pair rather than
    as the caption of the healthy row. Returns the y where the names end."""
    x0 = cx - pitch * (len(groups) - 1) / 2
    bottom = y
    for i, (name, _icons) in enumerate(groups):
        lb = fig.label(name, x0 + pitch * i, y, pt=pt or fig.theme.note_pt, anchor="middle",
                       va="top", color=color, max_w=max_w or pitch - 1.5, leading=1.15)
        bottom = max(bottom, lb.y + lb.h)
    return bottom


def fork(fig, x0, y0, x1, y_up, y_dn, *, color=GREY, lw=0.6):
    """One stem out of the measurement, then two elbow arrows, up to the healthy row and
    down to the damaged row. Returns the two corner x positions the labels sit over."""
    fig.line(x0, y0, x0 + 4.0, y0, stroke=color, lw=lw)
    for yt in (y_up, y_dn):
        fig.arrow((x0 + 4.0, y0), (x1, yt), color=color, lw=lw, elbow="v")
    return x0 + 4.0
