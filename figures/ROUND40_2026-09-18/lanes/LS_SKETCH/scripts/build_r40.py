#!/usr/bin/env python3
"""Round 40, lane LS_SKETCH: every sketch band and Figure 6 WITHOUT THE BRAIN, a compact 5d band, the Supp 16 band with one heart.
Rebuilt from the round-36 build (band_1f, Main_Fig6: ROUND36_2026-09-10/LS6_SKETCH/scripts/build_r36.py) and the round-34 build
(band_3d, band_4e, band_5d, band_s16: ROUND34_2026-09-10/LS4_SKETCH/scripts/build_r34.py), their code paths copied here, with
Alen's round-40 comments and nothing else. The build tree is this lane's copy of the round-36 tree (work/_build) with the organ
lists of r18_common.py and r18_family.py cut to FOUR groups.

  band_1f_r40          968.66 x 204.09 pt   Figure 1 panel f   (reads "most predictive of 141 measurements", as delivered on V16)
  band_3d_r40          952.73 x 198.43 pt   Figure 3 panel d
  band_4e_r40          970.79 x 209.76 pt   Figure 4 panel e
  band_5d_r40          968.94 x 198.43 pt   Figure 5 panel d, full width (kept for the coordinator)
  band_5d_compact_r40  472.00 x H pt        Figure 5 panel d, COMPACT, for the right-half quadrant slot (H = H_COMPACT_PT)
  band_s16_r40         496.80 x  79.00 pt   Supp Fig 16 top band
  Main_Fig6_r40        968.94 x 510.24 pt   Figure 6 without the title strip (the coordinator adds "Figure 6" at 13 pt)

WHAT CHANGES (round 40), and nothing else:
  organs   the BRAIN leaves every organ row (no brain condition has a positive association; dementia and Parkinson's inverse).
           FOUR groups: heart, lungs, pancreas (metabolism), kidney + liver, in that order, at the SAME pitch, both rows the
           same. The INK of the four-group block is centred on the same row centre as before (the kidney + liver group is
           now the last one and its liver reaches further right than the brain did, so centring the grid of group centres
           would leave the ink 2.6 to 3.0 mm right of the centre; the grid is moved left by that offset, ink_offset()).
           The arrows to the organs keep their round-34/36 gap to the first icon (the block starts half a pitch further
           right, so each organ arrow lengthens by that much; the fork arrows of 1f lengthen the same way and their labels
           follow the round-35 rule: centred on the arrow's mid-length x, one pill height clear of the arrow's edge).
  1f       "most predictive of 141 measurements" (the V16 sheet carries this wording, re-stamped in round 37; the r36 SVG
           still says 197).
  5d       a COMPACT band (new) for a quadrant slot: width 472 pt, height H_COMPACT_PT (aim 220 to 260), the same elements
           as the full band (mask + "treat sleep apnea", grey fork to two stacked traces, "90% SpO2" + dashed line, grey
           arrows to the four coloured / four grey organs), icons at least 9 pt tall, text at least 8 pt, balanced on its
           page by the round-23 rule (white above the first ink = white below the last ink, solved on real ink at 600 dpi
           with Ghostscript, pinned in SHIFT_COMPACT).
  s16      TWO named icons: one heart "heart attack or failure" and the human figure "death", centred under the question.
  Fig 6    the two organ panels hold four organs each; the panel rectangles keep their round-36 geometry (pinned by a
           five-group probe), the four-group ink centred in each (PANEL_FIT = True would shrink the panels to the block).

Render: FigureLab exports SVG -> PDF/PNG through Playwright; the bundled Chromium is missing on this machine (18 Sep), so
Renderer is pointed at the installed Google Chrome (CHROME), under the shared render lock. PyMuPDF touches nothing but the
freshly exported flat band's page box (l1_common.exact_boxes), in this process, which runs under wd_run.sh.

    python3 build_r40.py --preview             # 150 dpi previews into work/previews/ (no PDF)
    python3 build_r40.py --final               # SVG + 300 dpi PNG + PDF into work/build/ (under the render lock)
    python3 build_r40.py --final --only band_1f_r40,band_s16_r40
    python3 build_r40.py --solve-compact       # re-solve the compact band's balance shift on real ink (prints the value)
"""
from __future__ import annotations
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths

import argparse
import hashlib
import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path

LANE = Path(f"{paths.FIGURE_ROOT}/ROUND40_2026-09-18/lanes/LS_SKETCH")
WORK = LANE / "work"
BUILD = WORK / "_build"                       # this lane's copy of the round-36 LS6_SKETCH build tree, organ lists cut to four
OUT = WORK / "build"
SCRATCH = WORK / "_solve"
LOCK = f"{paths.FIGURE_ROOT}/ROUND30_2026-09-08/_render.lock"
CHROME = paths.CHROME
GS = paths.GS
sys.path.insert(0, str(BUILD))

import l1_common as L                                                        # noqa: E402
L.PREVIEWS = WORK / "previews"
from l1_common import (AZ, CARD, GREEN_D, GREY, GREY_L, INK, MM_PER_PT, WHITE, airflow, finish,   # noqa: E402
                       house, new, organ_row_named, text_width_mm)
from r18_common import ORGANS_BAD, ORGANS_OK, ORGANS_OK_R36, organ_row as organ_row_18     # noqa: E402
from r18_family import ORGANS_DARK, ORGANS_LIVE, ORGANS_LIVE_R36, airflow_ramp, organ_row as organ_row_f   # noqa: E402
from r20_common import severity_ladder                                       # noqa: E402
from build_r21_fig6 import ICON_HALF, PANEL_PAD, PANEL_RX                    # noqa: E402
import build_r20_le as B45                                                   # noqa: E402
import build_r22 as B22                                                      # noqa: E402
import build_r23 as B23                                                      # noqa: E402
from build_r23 import FAINT, H6, W6_PT                                       # noqa: E402
import bioglyph.render as _render                                            # noqa: E402

# ---------------------------------------------------------------- the render backend: the installed Google Chrome through Playwright
RENDER_INFO = {"backend": "playwright chromium, executable_path = Google Chrome", "chrome": CHROME, "version": None}


def _renderer_init(self, force=None):
    from playwright.sync_api import sync_playwright
    self.kind = "chromium"
    self._pw = sync_playwright().start()
    self._browser = self._pw.chromium.launch(headless=True, executable_path=CHROME)
    RENDER_INFO["version"] = self._browser.version


_render.Renderer.__init__ = _renderer_init

# ---------------------------------------------------------------- wording (Alen's, or already on the sheets)
PRESERVED, LOW = "Preserved oxygen", "Low oxygen"      # the 1f and Figure 6 labels (plain bold blue text since round 36)
LABEL90 = "90% SpO2"
INK21 = "#1a1d21"                                      # the sheets' ink
AMBER = "#d55e00"                                      # sleep apnea, burnt orange
AMBER_L = ["#f6d6c2", "#eeb08f", "#e4864f", "#d55e00"]
RANK_CLAIM_R36 = B22.RANK_CLAIM                        # "most predictive of 197 measurements" (the r36 SVG)
RANK_CLAIM_R40 = RANK_CLAIM_R36.replace("197", "141")  # the V16 sheet's wording (round 37 re-stamp): "most predictive of 141 measurements"
assert RANK_CLAIM_R40 == "most predictive of 141 measurements"
S16_TITLE = "Does the oxygen signal predict future disease at home?"
S16_ORGANS_R34 = [("heart", "heart failure"), ("heart", "heart attack"), ("sv-brain-2", "stroke"), ("human-grey", "death")]
S16_ORGANS_R40 = [("heart", "heart attack or failure"), ("human-grey", "death")]
S16_PITCH_R40 = 24.0                                   # the two icons, centred under the question (r34: four at 16.5)
S16_NAME_MAX_W_R40 = 18.0                              # wraps "heart attack or failure" as "heart attack" / "or failure" (r34: 16.0)
THRESHOLD_DASH = True
PANELS_FLUSH = True
PANEL_FIT = False                                      # Figure 6: False keeps the round-36 panel rectangles, True shrinks them to the four-group block

# ---------------------------------------------------------------- the label rule (round 36, kept)
LABEL_BLUE = AZ                                        # "#0187d1", the trace's blue
LABEL_BOLD = True
PILL_PAD_X, PILL_PAD_Y = 3.0, 1.4                      # the former pill's pads: reproduce the pill box in the record and the clearance rule

# ---------------------------------------------------------------- the one arrow rule (round 34, kept)
ARROW_COLOR = GREY                                     # "#8a9099"
ARROW_LW = 0.8
HEAD_LEN = max(2.3, 3.4 * ARROW_LW)                    # 2.72 mm

# ---------------------------------------------------------------- the dips convention (round 31, kept)
DIPS_LOW = (0.09, 0.25, 0.41, 0.57, 0.73, 0.89)
DIPS_PRESERVED = (0.30, 0.68)
DEEP_DEPTH, DEEP_SIGMA = 0.68, 0.035

# ---------------------------------------------------------------- fork-band geometry (round 32, kept)
TW_R31 = 38.0
LABEL90_GAP = 1.2
LABEL90_W = text_width_mm(LABEL90, 9.5, False)         # 16.21 mm
HEAD_CLEAR = 1.5
HEAD_BACK = LABEL90_GAP + LABEL90_W + HEAD_CLEAR       # 18.9 mm before the trace
ORGAN_ARROW_GAP = 4.0
TX_3D, TX_4E5D = 118.668, 130.568                      # the round-32 trace-pair x
ROW_X_1F_R36 = 162.0                                   # the round-36 fork tips' x on 1f (the gap to the first icon is kept, the x moves)

# ---------------------------------------------------------------- the balance shifts (rounds 34 to 36, kept; the compact band solved tonight)
SHIFT = {"band_1f_r40": 1.63, "band_3d_r40": 1.96, "band_4e_r40": 2.98, "band_5d_r40": 0.43, "Main_Fig6_r40": 0.0}
H_COMPACT_PT = 220.0                                   # the compact band's height (brief: at most 300, aim 220 to 260)
SHIFT_COMPACT = -0.71                                  # solved by --solve-compact 2026-09-18 14:20 (work/solve_compact_r40.json): the solver oscillates by one 600 dpi pixel (0.04 mm) between -0.699 and -0.72; whites 15.33 / 15.37 mm


# ================================================================ shared drawing (rounds 32 to 36, unchanged)

def trace(fig, x, y, w, h, *, color, dips, ref_at=0.5, lw=0.6, label=LABEL90, label_pt=None, label_bold=False,
          label_gap=LABEL90_GAP, ref_lw=0.3, ref_dash="1.4 1.0", label_color=INK21, ref_color=INK21, label_side="left"):
    y90 = y + h * ref_at
    n = 220
    pts = []
    for i in range(n + 1):
        t = i / n
        py = y + h * 0.24 + h * 0.03 * math.sin(t * 2 * math.pi * 7)
        for c in dips:
            py += h * DEEP_DEPTH * math.exp(-((t - c) / DEEP_SIGMA) ** 2)
        pts.append((x + w * t, py))
    fig.line(x, y90, x + w, y90, stroke=ref_color, lw=ref_lw, dash=ref_dash if THRESHOLD_DASH else None)
    if label_side == "left":
        lb = fig.label(label, x - label_gap, y90, pt=label_pt or fig.theme.note_pt, color=label_color, anchor="end", va="middle", bold=label_bold)
    else:
        lb = fig.label(label, x + w + label_gap, y90, pt=label_pt or fig.theme.note_pt, color=label_color, anchor="start", va="middle", bold=label_bold)
    seg, lenses = [], 0
    for px, py in pts + [(x + w, y90 - 1)]:
        if py > y90:
            seg.append((px, py))
        elif seg:
            fig.polygon([(seg[0][0], y90)] + seg + [(seg[-1][0], y90)], fill=color, opacity=0.35)
            seg = []
            lenses += 1
    fig.path("M" + " L".join(f"{px:.2f} {py:.2f}" for px, py in pts), stroke=color, lw=lw)
    return {"y90": round(y90, 3), "n_dips": len(dips), "lenses": lenses, "color": color,
            "trough_below_90_mm": round(max(py - y90 for _, py in pts), 3),
            "ink_top": round(min(py for _, py in pts), 3), "ink_bottom": round(max(py for _, py in pts), 3),
            "label_box": [round(lb.x, 3), round(lb.y, 3), round(lb.x + lb.w, 3), round(lb.y + lb.h, 3)],
            "label_pt": label_pt or fig.theme.note_pt, "label_color": label_color, "ref_color": ref_color, "ref_lw": ref_lw,
            "ref_dash": ref_dash if THRESHOLD_DASH else None, "label_side": label_side,
            "trace_box": [round(x, 3), round(y, 3), round(x + w, 3), round(y + h, 3)]}


def trace_pair(fig, tx, y_top, y_bot, tw, th, *, lw=0.6):
    top = trace(fig, tx, y_top, tw, th, color=AZ, dips=DIPS_PRESERVED, lw=lw)
    bot = trace(fig, tx, y_bot, tw, th, color=AZ, dips=DIPS_LOW, lw=lw)
    return {"top": top, "bottom": bot, "labels": {}}


def arrow(fig, p0, p1):
    fig.arrow(p0, p1, color=ARROW_COLOR, lw=ARROW_LW)
    return {"p0": [round(v, 3) for v in p0], "tip": [round(v, 3) for v in p1], "color": ARROW_COLOR, "lw": ARROW_LW,
            "head_len_mm": round(HEAD_LEN, 2), "length_mm": round(math.hypot(p1[0] - p0[0], p1[1] - p0[1]), 2)}


def fork_arrows(fig, src, tx, y_top, y_bot, th):
    tips = [(tx - HEAD_BACK, y + th * 0.5) for y in (y_top, y_bot)]
    recs = [arrow(fig, src, t) for t in tips]
    return {"src": [round(v, 3) for v in src], "tips": [r["tip"] for r in recs], "tips_x": round(tx - HEAD_BACK, 3),
            "color": ARROW_COLOR, "lw": ARROW_LW, "head_len_mm": round(HEAD_LEN, 2), "length_mm": [r["length_mm"] for r in recs],
            "angle_deg": [round(math.degrees(math.atan2(abs(t[1] - src[1]), t[0] - src[0])), 2) for t in tips]}


def box(p):
    return [round(p.x, 3), round(p.y, 3), round(p.x + p.w, 3), round(p.y + p.h, 3)]


def pill_h(fig):
    return fig.theme.label_pt * 0.352778 + 2 * PILL_PAD_Y      # 6.328 mm


def pill_w(fig, text):
    return text_width_mm(text, fig.theme.label_pt, True) + 2 * PILL_PAD_X


def seg_dist(px, py, a, b):
    ax, ay = a
    bx, by = b
    dx, dy = bx - ax, by - ay
    t = max(0.0, min(1.0, ((px - ax) * dx + (py - ay) * dy) / (dx * dx + dy * dy)))
    return math.hypot(px - (ax + t * dx), py - (ay + t * dy))


def box_seg_clearance(bx, a, b):
    x0, y0, x1, y1 = bx
    corners = [(x0, y0), (x1, y0), (x0, y1), (x1, y1)]
    d = min(seg_dist(px, py, a, b) for px, py in corners)
    for ex, ey in (a, b):
        cx_, cy_ = min(max(ex, x0), x1), min(max(ey, y0), y1)
        d = min(d, math.hypot(ex - cx_, ey - cy_))
    return d


def plain_label(fig, text, cx, cy, *, was_fill, was_text_color):
    pt = fig.theme.label_pt
    lb = fig.label(text, cx, cy, pt=pt, bold=LABEL_BOLD, color=LABEL_BLUE, anchor="middle", va="middle")
    ph, pw = pill_h(fig), pill_w(fig, text)
    return {"text": text, "cx": round(cx, 3), "cy": round(cy, 3), "pt": pt, "bold": LABEL_BOLD, "color": LABEL_BLUE,
            "box": box(lb), "w": round(lb.w, 3), "h": round(lb.h, 3), "pill": False, "rect": False, "rule": False,
            "was_pill": {"fill": was_fill, "text_color": was_text_color, "w": round(pw, 3), "h": round(ph, 3),
                         "box": [round(cx - pw / 2, 3), round(cy - ph / 2, 3), round(cx + pw / 2, 3), round(cy + ph / 2, 3)]}}


def label_beside_arrow(fig, text, p0, p1, *, side, was_fill, was_text_color, y_from=None):
    """Round 35's placement (the label's nearest corner one pill height clear of the arrow's visible edge, centred on the
    arrow's mid-length x), drawn as plain text since round 36. ROUND 40: the fork arrows of 1f are longer and shallower,
    so with y_from = the round-36 tip the label keeps the y the rule gave it on the round-36 arrow (its V16 y) and only its
    x follows the new arrow's mid-length; the clearance from the new arrow is recomputed and must still be one pill height."""
    ph = pill_h(fig)
    pw = pill_w(fig, text)

    def rule(p1_):
        cx_ = (p0[0] + p1_[0]) / 2
        x0, x1 = cx_ - pw / 2, cx_ + pw / 2
        slope = (p1_[1] - p0[1]) / (p1_[0] - p0[0])
        theta_ = math.atan(abs(slope))
        gap_v = (ph + ARROW_LW / 2) / math.cos(theta_)
        line_y = lambda x: p0[1] + slope * (x - p0[0])
        if side == "above":
            cy_ = min(line_y(x0), line_y(x1)) - gap_v - ph / 2
        else:
            cy_ = max(line_y(x0), line_y(x1)) + gap_v + ph / 2
        return cx_, cy_, theta_

    cx, cy_new, theta = rule(p1)
    cy = cy_new
    if y_from is not None:
        cx_r36, cy, theta_r36 = rule(y_from)
    rec = plain_label(fig, text, cx, cy, was_fill=was_fill, was_text_color=was_text_color)
    clear = box_seg_clearance(rec["box"], p0, p1) - ARROW_LW / 2
    rec.update(side=side, on_arrow=False, arrow_mid=[round(cx, 3), round((p0[1] + p1[1]) / 2, 3)], arrow_angle_deg=round(math.degrees(theta), 2),
               clearance_from_arrow_edge_mm=round(clear, 3), pill_clearance_was_mm=round(ph, 3))
    if y_from is not None:
        rec.update(y_kept_from_r36_arrow=True, r36_centre=[round(cx_r36, 3), round(cy, 3)], r36_arrow_angle_deg=round(math.degrees(theta_r36), 2),
                   cy_the_rule_would_give_on_the_new_arrow=round(cy_new, 3))
    return rec


# ================================================================ the four-group organ rows (round 40)

def row_extents(groups, pitch, scale, *, family=False):
    """Place one row on a probe figure (grid centred on 150) and return the ink extents about the grid centre, the offset of
    the ink centre from the grid centre, and every icon's box. Nothing is rendered."""
    probe = new(400, 60)
    if family:
        organ_row_f(probe, 150.0, 30.0, dark=False, pitch=pitch, scale=scale, names=False, groups=groups)
    else:
        organ_row_18(probe, 150.0, 30.0, groups, pitch=pitch, icon_scale=scale)
    x0 = min(p.x for p in probe._placed)
    x1 = max(p.x + p.w for p in probe._placed)
    return {"n_groups": len(groups), "n_icons": len(probe._placed), "pitch": pitch, "scale": scale,
            "left_of_centre": round(150.0 - x0, 3), "right_of_centre": round(x1 - 150.0, 3),
            "ink_width": round(x1 - x0, 3), "ink_offset": round((x0 + x1) / 2 - 150.0, 3),
            "icons": [{"name": p.name, "dx": round(p.x + p.w / 2 - 150.0, 3), "w": round(p.w, 2), "h": round(p.h, 2)} for p in probe._placed]}


def organ_rows(fig, ocx, ys, *, pitch, scale=1.0, family=False):
    """Both rows (coloured at ys[0], grey at ys[1]) with their INK centred on ocx. Returns the geometry for the record."""
    if family:
        ext = row_extents(ORGANS_LIVE, pitch, scale, family=True)
        old = row_extents(ORGANS_LIVE_R36, pitch, scale, family=True)
    else:
        ext = row_extents(ORGANS_OK, pitch, scale)
        old = row_extents(ORGANS_OK_R36, pitch, scale)
    gcx = ocx - ext["ink_offset"]
    n0 = len(fig._placed)
    for y, dark in ((ys[0], False), (ys[1], True)):
        if family:
            organ_row_f(fig, gcx, y, dark=dark, pitch=pitch, scale=scale, names=False)
        else:
            organ_row_18(fig, gcx, y, ORGANS_BAD if dark else ORGANS_OK, pitch=pitch, icon_scale=scale)
    placed = fig._placed[n0:]
    rows = {}
    for y, key in ((ys[0], "top"), (ys[1], "bottom")):
        ic = [p for p in placed if abs(p.y + p.h / 2 - y) < 0.01]
        rows[key] = {"y": round(y, 3), "n_icons": len(ic), "names": [p.name for p in ic],
                     "x_first": round(min(p.x for p in ic), 3), "x_last": round(max(p.x + p.w for p in ic), 3)}
    assert rows["top"]["n_icons"] == rows["bottom"]["n_icons"] == 5, rows       # four groups, five icons (kidney + liver is one group)
    assert rows["top"]["x_first"] == rows["bottom"]["x_first"] and rows["top"]["x_last"] == rows["bottom"]["x_last"], rows
    assert not any("brain" in p.name for p in placed), [p.name for p in placed]
    return {"row_centre": round(ocx, 3), "grid_centre": round(gcx, 3), "ink_centre": round(ocx, 3), "ink_offset_applied": round(-ext["ink_offset"], 3),
            "pitch": pitch, "scale": scale, "groups": [n for n, _ in (ORGANS_LIVE if family else ORGANS_OK)],
            "n_groups": 4, "rows": rows, "ink_x": [rows["top"]["x_first"], rows["top"]["x_last"]],
            "r36_five_groups": {"grid_centre": round(ocx, 3), "ink_x": [round(ocx - old["left_of_centre"], 3), round(ocx + old["right_of_centre"], 3)],
                                "ink_offset_unapplied": old["ink_offset"]},
            "first_icon_left_r36": round(ocx - old["left_of_centre"], 3)}


# ================================================================ Figure 1, panel f

def band_1f(final, *, stem="band_1f_r40", outdir=OUT, shift=None):
    """build_r36.band_1f with the four-group rows, the fork tips following the block (same gap to the first icon), the labels
    by the round-35 rule on the longer arrows, and the V16 wording "141"."""
    W, H = B22.W1_PT * MM_PER_PT, B22.H1_BOT
    fig = new(W, H)
    t = fig.theme
    d = SHIFT[stem] if shift is None else shift
    Y_OK, Y_BAD = 23.5 - d, 51.5 - d
    CY = (Y_OK + Y_BAD) / 2
    TX, TW, TH = 21.0, 67.0, 26.0
    TY = CY - TH / 2
    FORK_X = TX + TW + 4.0
    OCX, PITCH = 248.0, 37.0
    r1 = fig.label(RANK_CLAIM_R40, TX + TW / 2, TY - 5.5, pt=t.note_pt, bold=True, color=AZ, anchor="middle", va="middle")
    tr = trace(fig, TX, TY, TW, TH, color=AZ, dips=(0.18, 0.49, 0.79), lw=0.6)
    fig.line(TX, TY + TH + 1.4, TX + TW, TY + TH + 1.4, stroke=GREY, lw=0.3)
    r3 = fig.label("one night of sleep", TX + TW / 2, TY + TH + 2.6, pt=t.note_pt, color=GREY, anchor="middle", va="top")
    org = organ_rows(fig, OCX, (Y_OK, Y_BAD), pitch=PITCH)
    gap_r36 = org["first_icon_left_r36"] - ROW_X_1F_R36                   # 7.125 mm, the round-36 tip-to-first-icon gap
    ROW_X = org["ink_x"][0] - gap_r36
    src, tips = (FORK_X, CY), [(ROW_X, Y_OK), (ROW_X, Y_BAD)]
    tips_r36 = [(ROW_X_1F_R36, Y_OK), (ROW_X_1F_R36, Y_BAD)]          # the round-36 fork (its labels' y is kept)
    fa = [arrow(fig, src, tp) for tp in tips]
    labels = {"preserved": label_beside_arrow(fig, PRESERVED, src, tips[0], side="above", was_fill=FAINT, was_text_color=INK, y_from=tips_r36[0]),
              "low": label_beside_arrow(fig, LOW, src, tips[1], side="below", was_fill=AZ, was_text_color=WHITE, y_from=tips_r36[1])}
    for k, p in labels.items():
        assert p["clearance_from_arrow_edge_mm"] >= p["pill_clearance_was_mm"] - 0.01, (k, p)
        assert 0.0 < p["box"][1] and p["box"][3] < H and 0.0 < p["box"][0] and p["box"][2] < W, ("label inside the band", k, p["box"])
        assert p["box"][2] < org["ink_x"][0] - 5.0, ("label clear of the organ icons", k, p["box"], org["ink_x"])
        rb = box(r1)
        assert not (p["box"][0] < rb[2] and p["box"][2] > rb[0] and p["box"][1] < rb[3] and p["box"][3] > rb[1]), ("label clear of the rank line", k)
        assert p["box"][0] > TX + TW + 4.0, ("label right of the trace and fork", k, p["box"])
    assert org["ink_x"][1] < W - 5.0, ("organ block inside the band", org["ink_x"], W)
    rec = finish(fig, stem, outdir, final=final, attribution=None)
    rec.update(shift=round(d, 3), trace=tr, fork={"src": [round(v, 3) for v in src], "tips": [a["tip"] for a in fa], "color": ARROW_COLOR, "lw": ARROW_LW,
                                                    "head_len_mm": round(HEAD_LEN, 2), "length_mm": [a["length_mm"] for a in fa],
                                                    "tips_x": round(ROW_X, 3), "tips_x_r36": ROW_X_1F_R36, "gap_to_first_icon_mm": round(gap_r36, 3)},
               labels=labels, label_blue=LABEL_BLUE, letter_zone=[0.0, 0.0, B22.LETTER_ZONE_W, B22.LETTER_ZONE_H], rows_y=[round(Y_OK, 3), round(Y_BAD, 3)],
               row_x=round(ROW_X, 3), organs=org, left_block={"rank": box(r1), "night": box(r3)}, rank_claim=RANK_CLAIM_R40,
               deleted=["brain"], wording_delta={RANK_CLAIM_R36: RANK_CLAIM_R40},
               moved="round 40: the brain leaves both rows (four groups, ink centred on x 248, pitch 37); the fork tips move right with the block (same 7.1 mm gap to the first icon) and the two labels keep their V16 y, their x follows the arrows' mid-length (round-35 rule); '197' reads '141' as on the V16 sheet; nothing else (shift 1.63 kept)")
    return rec


# ================================================================ Figure 3, panel d

def band_3d(final, *, stem="band_3d_r40", outdir=OUT, shift=None):
    d = SHIFT[stem] if shift is None else shift
    W = B23.W3_PT * MM_PER_PT
    fig = new(W, B23.H3)
    Y_OK, Y_BAD = 23.5 - d, 50.0 - d
    TW, TH = TW_R31, 14.0
    Y_FLAT, Y_DIPS = Y_OK - TH / 2, Y_BAD - TH / 2
    OCX, PITCH = 262.0, 28.5
    TX = TX_3D
    bed = fig.place("sleeper-faceB", 36.0, 34.0 - d, w=54)
    dl = fig.label(B23.DURATIONS, 36.0, bed.y + bed.h + 0.6, bold=True, anchor="middle", va="top", color=GREEN_D)
    fa = fork_arrows(fig, (62.0, 34.0 - d), TX, Y_FLAT, Y_DIPS, TH)
    tp = trace_pair(fig, TX, Y_FLAT, Y_DIPS, TW, TH)
    brx = TX + TW + ORGAN_ARROW_GAP
    org = organ_rows(fig, OCX, (Y_OK, Y_BAD), pitch=PITCH, scale=0.95)
    tip_r34 = OCX - PITCH * 2 - 9.0                                          # 196.0, the round-34 tip
    gap_r34 = org["first_icon_left_r36"] - tip_r34                          # 4.37 mm
    tip_x = org["ink_x"][0] - gap_r34
    oa = [arrow(fig, (brx, yt), (tip_x, yt)) for yt in (Y_OK, Y_BAD)]
    assert org["ink_x"][1] < W - 5.0, (org["ink_x"], W)
    rec = finish(fig, stem, outdir, final=final, attribution=None)
    rec.update(shift=round(d, 3), traces=tp, fork=fa, organ_arrows=oa, organ_arrow=[round(brx, 3), round(tip_x, 3)], organ_arrow_r34=[round(brx, 3), round(tip_r34, 3)],
               gap_to_first_icon_mm=round(gap_r34, 3), tx=round(TX, 3), organs=org,
               bed_box=box(bed), duration_label_box=box(dl), deleted=["brain"],
               moved=f"round 40: the brain leaves both rows (four groups, ink centred on x {OCX}, pitch {PITCH}, scale 0.95); the organ arrows lengthen to the block (same {gap_r34:.2f} mm gap to the first icon); nothing else (shift {d} kept)")
    return rec


# ================================================================ Figure 4, panel e

def band_4e(final, *, stem="band_4e_r40", outdir=OUT, shift=None):
    class R4:
        W4_PT, H4 = B23.W4_PT, B23.H4; BED_W = 54.0; RX, RW = 4.0, 74.0; BX = RX + RW / 2
        BED_INK_TOP = 9.6; RAMP_CY, RAMP_H = 43.1, 8.5; LADDER_Y = 51.1; SEV_LABEL_Y = 59.3
    meas = json.load(open(f"{paths.FIGURE_ROOT}/ROUND29_2026-09-08/LC_fig4/work/bed_measure.json"))
    ink_frac_y = tuple(meas["ink_frac_y"])
    d = SHIFT[stem] if shift is None else shift
    dd = d - B23.SHIFT4
    W = R4.W4_PT * MM_PER_PT
    fig = new(W, R4.H4)
    PITCH, ORG_CX = B45.PITCH, B45.ORG_CX
    TW, TH = TW_R31, B45.TH
    TX = TX_4E5D
    bed_h = R4.BED_W * 82.0 / 122.0
    bed_cy = R4.BED_INK_TOP - ink_frac_y[0] * bed_h + bed_h / 2 - dd
    bed = fig.place(B45.SLEEPER, R4.BX, bed_cy, w=R4.BED_W)
    airflow_ramp(fig, R4.RX, R4.RAMP_CY - dd, R4.RW, R4.RAMP_H, color=AMBER)
    severity_ladder(fig, R4.RX, R4.LADDER_Y - dd, R4.RW - 3.0, left_text=B45.APNEA_NONE, right_text=B45.APNEA_SEVERE,
                    tri_h=2.0, ladder=AMBER_L, text_color=AMBER)
    sev = fig.label("sleep apnea severity", R4.BX, R4.SEV_LABEL_Y - dd, anchor="middle", va="top", bold=True, color=AMBER)
    CY_TOP, CY_BOT = 23.5 - d, 52.5 - d
    fa = fork_arrows(fig, (84.0, 40.0 - d), TX, CY_TOP - TH / 2, CY_BOT - TH / 2, TH)
    tp = trace_pair(fig, TX, CY_TOP - TH / 2, CY_BOT - TH / 2, TW, TH)
    from_x = TX + TW + ORGAN_ARROW_GAP
    org = organ_rows(fig, ORG_CX, (CY_TOP, CY_BOT), pitch=PITCH, family=True)
    tip_r34 = ORG_CX - 63.0                                                 # 199.0, the round-34 tip
    gap_r34 = org["first_icon_left_r36"] - tip_r34                          # 0.94 mm
    tip_x = org["ink_x"][0] - gap_r34
    oa = [arrow(fig, (from_x, cy), (tip_x, cy)) for cy in (CY_TOP, CY_BOT)]
    assert org["ink_x"][1] < W - 5.0, (org["ink_x"], W)
    rec = finish(fig, stem, outdir, final=final, attribution=None)
    rec.update(shift=round(d, 3), traces=tp, fork=fa, organ_arrows=oa, organ_arrow=[round(from_x, 3), round(tip_x, 3)], organ_arrow_r34=[round(from_x, 3), round(tip_r34, 3)],
               gap_to_first_icon_mm=round(gap_r34, 3), tx=round(TX, 3), organs=org,
               bed_box_mm=[round(bed.x, 3), round(bed.y, 3), round(bed.w, 3), round(bed.h, 3)], sev_label_box=box(sev), deleted=["brain"],
               moved=f"round 40: the brain leaves both rows (four groups, ink centred on x {ORG_CX}, pitch {PITCH}); the organ arrows lengthen to the block (same {gap_r34:.2f} mm gap to the first icon); nothing else (shift {d} kept)")
    return rec


# ================================================================ Figure 5, panel d (full width)

def band_5d(final, *, stem="band_5d_r40", outdir=OUT, shift=None):
    d = SHIFT[stem] if shift is None else shift
    W = B45.W5_PT * MM_PER_PT
    fig = new(W, B45.H5)
    TW, TH = TW_R31, B45.TH
    TX = TX_4E5D
    BX = 40.0
    mask = fig.place(B45.MASK, BX, 28.0 - d, h=26.0, flip_h=True)
    tl = fig.label("treat sleep apnea", BX, 43.5 - d, anchor="middle", va="top", bold=True, color=AMBER)
    CY_TOP, CY_BOT = 23.0 - d, 48.5 - d
    fa = fork_arrows(fig, (78.0, 34.0 - d), TX, CY_TOP - TH / 2, CY_BOT - TH / 2, TH)
    tp = trace_pair(fig, TX, CY_TOP - TH / 2, CY_BOT - TH / 2, TW, TH)
    from_x = TX + TW + ORGAN_ARROW_GAP
    org = organ_rows(fig, B45.ORG_CX, (CY_TOP, CY_BOT), pitch=B45.PITCH, family=True)
    tip_r34 = B45.ORG_CX - 63.0
    gap_r34 = org["first_icon_left_r36"] - tip_r34
    tip_x = org["ink_x"][0] - gap_r34
    oa = [arrow(fig, (from_x, cy), (tip_x, cy)) for cy in (CY_TOP, CY_BOT)]
    assert org["ink_x"][1] < W - 5.0, (org["ink_x"], W)
    rec = finish(fig, stem, outdir, final=final, attribution=None)
    rec.update(shift=round(d, 3), traces=tp, fork=fa, organ_arrows=oa, organ_arrow=[round(from_x, 3), round(tip_x, 3)], organ_arrow_r34=[round(from_x, 3), round(tip_r34, 3)],
               gap_to_first_icon_mm=round(gap_r34, 3), tx=round(TX, 3), organs=org,
               mask_box=box(mask), treat_label_box=box(tl), deleted=["brain"],
               moved=f"round 40: the brain leaves both rows (four groups, ink centred on x {B45.ORG_CX}, pitch {B45.PITCH}); the organ arrows lengthen to the block (same {gap_r34:.2f} mm gap to the first icon); nothing else (shift {d} kept)")
    return rec


# ================================================================ Figure 5, panel d, COMPACT (new, round 40)

C = dict(W_PT=472.0,            # the right half of the sheet, x 497 to 969
         MASK_CX=18.5, MASK_H=19.0, LABEL_GAP=2.5,   # the sleeper with the CPAP mask and "treat sleep apnea" under it
         SRC_DX=2.0,            # the fork source, this far right of the mask's box
         RISE=18.0,             # the two rows at CY -/+ RISE
         FORK_RUN=18.0,         # horizontal run from the fork source to the tips
         TW=22.0, TH=13.0,      # the traces
         ORG_RUN=8.0,           # the organ arrows' run (from trace right edge + ORGAN_ARROW_GAP)
         ORG_GAP=1.0,           # tip to the first icon (the full band's 0.94)
         PITCH=16.0, SCALE=0.70)  # the organ rows (heart 9.45 mm tall, liver 7.0: all above the 9 pt = 3.2 mm floor)


def band_5d_compact(final, *, stem="band_5d_compact_r40", outdir=OUT, shift=None, h_pt=None):
    """The 5d band for a quadrant slot: the same elements as the full band, laid out at 472 pt width. Balanced by the
    round-23 rule (SHIFT_COMPACT, solved on real ink at 600 dpi with Ghostscript)."""
    d = SHIFT_COMPACT if shift is None else shift
    W = C["W_PT"] * MM_PER_PT
    H = (H_COMPACT_PT if h_pt is None else h_pt) * MM_PER_PT
    fig = new(W, H)
    t = fig.theme
    CY = H / 2 + d
    # the mask + label block centred on CY
    lab_h = new(W, H).label("treat sleep apnea", 0, 0, anchor="middle", va="top", bold=True, color=AMBER).h
    block_h = C["MASK_H"] + C["LABEL_GAP"] + lab_h
    mask_cy = CY - block_h / 2 + C["MASK_H"] / 2
    mask = fig.place(B45.MASK, C["MASK_CX"], mask_cy, h=C["MASK_H"], flip_h=True)
    tl = fig.label("treat sleep apnea", C["MASK_CX"], mask.y + mask.h + C["LABEL_GAP"], anchor="middle", va="top", bold=True, color=AMBER)
    src = (mask.x + mask.w + C["SRC_DX"], CY)
    CY_TOP, CY_BOT = CY - C["RISE"], CY + C["RISE"]
    TX = src[0] + C["FORK_RUN"] + HEAD_BACK
    TW, TH = C["TW"], C["TH"]
    fa = fork_arrows(fig, src, TX, CY_TOP - TH / 2, CY_BOT - TH / 2, TH)
    tp = trace_pair(fig, TX, CY_TOP - TH / 2, CY_BOT - TH / 2, TW, TH)
    from_x = TX + TW + ORGAN_ARROW_GAP
    tip_x = from_x + C["ORG_RUN"]
    ext = row_extents(ORGANS_LIVE, C["PITCH"], C["SCALE"], family=True)
    first_left = tip_x + C["ORG_GAP"]
    ocx = first_left + ext["left_of_centre"] + ext["ink_offset"]           # the ink centre such that the first icon sits at first_left
    org = organ_rows(fig, ocx, (CY_TOP, CY_BOT), pitch=C["PITCH"], scale=C["SCALE"], family=True)
    assert abs(org["ink_x"][0] - first_left) < 0.01, (org["ink_x"], first_left)
    oa = [arrow(fig, (from_x, cy), (tip_x, cy)) for cy in (CY_TOP, CY_BOT)]
    min_icon_h = min(p.h for p in fig._placed if p.name != B45.MASK)
    assert min_icon_h >= 9.0 * MM_PER_PT, ("icons at least 9 pt tall", min_icon_h)
    assert min(fig._fonts_pt) >= 8.0, ("text at least 8 pt", min(fig._fonts_pt))
    assert org["ink_x"][1] < W - 2.5 and tl.x > 2.5, ("block inside the band with a margin", org["ink_x"], tl.x, W)
    assert tl.y > fa["tips"][1][1] - 20, "label under the mask"
    # the label must clear the lower fork arrow
    clear = box_seg_clearance(box(tl), src, tuple(fa["tips"][1])) - ARROW_LW / 2
    assert clear > 2.0, ("'treat sleep apnea' clear of the lower fork arrow", clear)
    rec = finish(fig, stem, outdir, final=final, attribution=None)
    rec.update(shift=round(d, 3), h_pt=round(H / MM_PER_PT, 2), w_pt=round(W / MM_PER_PT, 2), params=C, cy=round(CY, 3), traces=tp, fork=fa,
               organ_arrows=oa, organ_arrow=[round(from_x, 3), round(tip_x, 3)], gap_to_first_icon_mm=C["ORG_GAP"], tx=round(TX, 3), organs=org,
               mask_box=box(mask), treat_label_box=box(tl), label_clearance_from_lower_arrow_mm=round(clear, 3), min_icon_h_mm=round(min_icon_h, 3),
               min_icon_h_pt=round(min_icon_h / MM_PER_PT, 2), deleted=["brain"],
               moved="round 40: NEW compact layout of the 5d band for the quadrant slot (472 pt wide); same elements as the full band, four organ groups, balanced on its page")
    return rec


# ================================================================ Figure 6

def fig6(final, *, stem="Main_Fig6_r40", outdir=OUT, dy=None, box_pad=8.0):
    d = SHIFT[stem] if dy is None else dy
    W = W6_PT * MM_PER_PT
    fig = new(W, H6)
    t = fig.theme
    CX = W / 2
    MID = H6 / 2 + d
    BED_W = 26.0
    BED_H = BED_W * 82.0 / 122.0
    PITCH_IN = 43.0
    TW6, TH6 = 34.0, 19.5
    ORG_ARROW = 20.0
    PANEL_H = 2 * (ICON_HALF + PANEL_PAD)
    L_PCX = 72.0
    R_PCX = W - L_PCX
    PITCH = 25.0
    LGAP = 0.7
    ph = pill_h(fig)
    PILL_CLEAR = ph

    probe = new(W, H6)
    lab_h = probe.label("sleep apnea", CX, 0.0, anchor="middle", va="top", bold=True, color=AMBER).h
    c_top, c_bot = -PITCH_IN - BED_H / 2, PITCH_IN + 11.0 + lab_h
    Y_PAP = MID - (c_top + c_bot) / 2
    Y_DUR, Y_APNEA = Y_PAP - PITCH_IN, Y_PAP + PITCH_IN
    k = BED_W / 20.0
    b1 = fig.place(B45.SLEEPER, CX, Y_DUR, w=BED_W)
    l1 = fig.label("sleep duration", CX, Y_DUR + 11.0, anchor="middle", va="top", bold=True, color=GREEN_D)
    m2 = fig.place(B45.MASK, CX, Y_PAP, h=24.0, flip_h=True)
    l2 = fig.label("PAP treatment", CX, Y_PAP + 14.0, anchor="middle", va="top", bold=True, color=INK)
    airflow(fig, CX, Y_APNEA - 10.0, 16 * k, 5.5 * k, AMBER)
    af_box = [CX - 8 * k, Y_APNEA - 10.0 - 2.75 * k, CX + 8 * k, Y_APNEA - 10.0 + 2.75 * k]
    b3 = fig.place(B45.SLEEPER, CX, Y_APNEA, w=BED_W)
    l3 = fig.label("sleep apnea", CX, Y_APNEA + 11.0, anchor="middle", va="top", bold=True, color=AMBER)
    parts = [box(p) for p in (b1, l1, m2, l2, b3, l3)] + [af_box]
    cx0, cx1 = min(p[0] for p in parts), max(p[2] for p in parts)
    cy0, cy1 = min(p[1] for p in parts), max(p[3] for p in parts)
    half_w = max(CX - cx0, cx1 - CX) + box_pad
    inbox = [round(CX - half_w, 3), round(cy0 - box_pad, 3), round(2 * half_w, 3), round(cy1 - cy0 + 2 * box_pad, 3)]
    RX0, RX1 = inbox[0], inbox[0] + inbox[2]
    assert abs((inbox[1] + inbox[3] / 2) - MID) < 0.02, ("the rectangle is centred on the middle line", inbox, MID)

    TRACE_TOP = MID - TH6 / 2
    TXL, TXR = L_PCX - TW6 / 2, R_PCX - TW6 / 2
    trL = trace(fig, TXL, TRACE_TOP, TW6, TH6, color=AZ, dips=DIPS_PRESERVED, lw=1.0, label_pt=t.note_pt, label_bold=True,
                label_gap=LGAP, ref_lw=0.5, ref_dash="1.4 1.2", label_side="right")
    trR = trace(fig, TXR, TRACE_TOP, TW6, TH6, color=AZ, dips=DIPS_LOW, lw=1.0, label_pt=t.note_pt, label_bold=True,
                label_gap=LGAP, ref_lw=0.5, ref_dash="1.4 1.2", label_side="left")
    TIP_L = trL["label_box"][2] + HEAD_CLEAR
    TIP_R = trR["label_box"][0] - HEAD_CLEAR
    inputs = [dict(arrow(fig, (RX0, MID), (TIP_L, MID)), side="left"), dict(arrow(fig, (RX1, MID), (TIP_R, MID)), side="right")]

    A_Y0 = TRACE_TOP + TH6 + 1.5
    if PANELS_FLUSH:
        PANEL_Y0 = inbox[1] + inbox[3] - PANEL_H
        A_Y1 = PANEL_Y0 - ARROW_LW * 0.25
        ORG_ARROW = round(A_Y1 - A_Y0, 3)
    else:
        A_Y1 = A_Y0 + ORG_ARROW
        PANEL_Y0 = A_Y1 + ARROW_LW * 0.25
    PANEL_CY = PANEL_Y0 + PANEL_H / 2
    oaL = arrow(fig, (L_PCX, A_Y0), (L_PCX, A_Y1))
    oaR = arrow(fig, (R_PCX, A_Y0), (R_PCX, A_Y1))

    # ---- the two organ panels: FOUR groups each, ink centred on the panel's centre (round 40)
    extL = row_extents(ORGANS_LIVE, PITCH, 1.0, family=True)
    old5 = row_extents(ORGANS_LIVE_R36, PITCH, 1.0, family=True)
    gcx_off = extL["ink_offset"]
    n0 = len(fig._placed)
    organ_row_f(fig, L_PCX - gcx_off, PANEL_CY, dark=False, pitch=PITCH, names=False)
    organ_row_f(fig, R_PCX - gcx_off, PANEL_CY, dark=True, pitch=PITCH, names=False)
    placed = fig._placed[n0:]
    left_icons = [p for p in placed if p.x < W / 2]
    right_icons = [p for p in placed if p.x > W / 2]
    assert len(left_icons) == len(right_icons) == 5 and not any("brain" in p.name for p in placed), [p.name for p in placed]
    half_new = max(max(L_PCX - p.x for p in left_icons), max(p.x + p.w - L_PCX for p in left_icons),
                   max(R_PCX - p.x for p in right_icons), max(p.x + p.w - R_PCX for p in right_icons)) + PANEL_PAD
    half_r36 = max(old5["left_of_centre"], old5["right_of_centre"]) + PANEL_PAD     # 59.39, the round-36 panels (grid centred, five groups)
    half = half_new if PANEL_FIT else half_r36
    panels = {"left": [round(L_PCX - half, 3), round(PANEL_Y0, 3), round(2 * half, 3), round(PANEL_H, 3)],
              "right": [round(R_PCX - half, 3), round(PANEL_Y0, 3), round(2 * half, 3), round(PANEL_H, 3)]}
    organs = {"left": {"row_centre": L_PCX, "grid_centre": round(L_PCX - gcx_off, 3), "ink_offset_applied": round(-gcx_off, 3), "pitch": PITCH,
                       "ink_x": [round(min(p.x for p in left_icons), 3), round(max(p.x + p.w for p in left_icons), 3)], "names": [p.name for p in left_icons]},
              "right": {"row_centre": round(R_PCX, 3), "grid_centre": round(R_PCX - gcx_off, 3), "ink_offset_applied": round(-gcx_off, 3), "pitch": PITCH,
                        "ink_x": [round(min(p.x for p in right_icons), 3), round(max(p.x + p.w for p in right_icons), 3)], "names": [p.name for p in right_icons]},
              "panel_half_r36": round(half_r36, 3), "panel_half_fit": round(half_new, 3), "panel_fit": PANEL_FIT,
              "pad_left_inside_panel_mm": round(min(p.x for p in left_icons) - (L_PCX - half), 3), "pad_right_inside_panel_mm": round((L_PCX + half) - max(p.x + p.w for p in left_icons), 3)}
    for key, icons, pc in (("left", left_icons, L_PCX), ("right", right_icons, R_PCX)):
        assert abs((min(p.x for p in icons) + max(p.x + p.w for p in icons)) / 2 - pc) < 0.01, ("ink centred on the panel", key)
        assert min(p.x for p in icons) >= pc - half + PANEL_PAD - 0.01 and max(p.x + p.w for p in icons) <= pc + half - PANEL_PAD + 0.01, ("icons inside the panel with the pad", key)

    n0 = len(fig._parts)
    for bx_ in (inbox, panels["left"], panels["right"]):
        fig.rect(*bx_, rx=PANEL_RX, fill=CARD, stroke=GREY_L, lw=0.3)
    moved = fig._parts[n0:]
    del fig._parts[n0:]
    fig._parts[0:0] = moved

    PY = (A_Y0 + A_Y1) / 2
    pwL, pwR = pill_w(fig, PRESERVED), pill_w(fig, LOW)
    cxL = L_PCX - ARROW_LW / 2 - PILL_CLEAR - pwL / 2
    cxR = R_PCX + ARROW_LW / 2 + PILL_CLEAR + pwR / 2
    labels = {}
    for key, text, cx_, fill, tc, ar, side in (("preserved", PRESERVED, cxL, FAINT, INK, oaL, "left"), ("low", LOW, cxR, AZ, WHITE, oaR, "right")):
        rec_ = plain_label(fig, text, cx_, PY, was_fill=fill, was_text_color=tc)
        clear = box_seg_clearance(rec_["box"], tuple(ar["p0"]), tuple(ar["tip"])) - ARROW_LW / 2
        rec_.update(beside_arrow_x=ar["p0"][0], side=side, on_arrow=False, clearance_from_arrow_edge_mm=round(clear, 3), pill_clearance_was_mm=round(ph, 3))
        labels[key] = rec_

    for p in parts:
        assert p[0] >= RX0 + box_pad - 0.01 and p[2] <= RX1 - box_pad + 0.01 and p[1] >= inbox[1] + box_pad - 0.01 and p[3] <= inbox[1] + inbox[3] - box_pad + 0.01, (p, inbox)
    assert TIP_L < RX0 - 25.0 and TIP_R > RX1 + 25.0, ("input arrows at least 25 mm of run", TIP_L, RX0, TIP_R, RX1)
    assert abs(trL["y90"] - MID) < 0.01 and abs(trR["y90"] - MID) < 0.01
    assert abs((TXL + TW6 / 2) + (TXR + TW6 / 2) - W) < 0.02, "traces mirror about the page centre"
    for key, p in labels.items():
        assert p["clearance_from_arrow_edge_mm"] >= ph - 0.01, (key, p)
        assert p["box"][1] > trL["trace_box"][3] + 2.0 and p["box"][3] < PANEL_Y0 - 2.0, ("label between the trace and the panel", key, p["box"])
    assert labels["preserved"]["box"][2] < L_PCX and labels["low"]["box"][0] > R_PCX
    assert labels["preserved"]["box"][2] < trL["label_box"][0] or labels["preserved"]["box"][3] < trL["label_box"][1], "left label clear of '90% SpO2'"
    assert panels["left"][0] > 5.0 and panels["right"][0] + panels["right"][2] < W - 5.0, ("panels inside the page margins", panels)
    if PANELS_FLUSH:
        assert abs(panels["left"][1] + panels["left"][3] - (inbox[1] + inbox[3])) < 0.01, ("panels flush with the rectangle's bottom", panels, inbox)
    assert inbox[1] < TRACE_TOP and inbox[1] < min(p["box"][1] for p in labels.values()), "the rectangle is the first ink"
    assert abs(inbox[1] - (H6 - inbox[1] - inbox[3])) < 0.02 or d != 0.0, ("rectangle centred on the page at dy 0", inbox)

    rec = finish(fig, stem, outdir, final=final, attribution=None)
    rec.update(dy=round(d, 3), cx=round(CX, 3), mid=round(MID, 3), trace_top=round(TRACE_TOP, 3), box_pad=box_pad,
               inputs_y=[round(Y_DUR, 3), round(Y_PAP, 3), round(Y_APNEA, 3)], input_boxes={"duration": box(b1), "PAP": box(m2), "apnea": box(b3)},
               input_labels={"duration": box(l1), "PAP": box(l2), "apnea": box(l3)}, airflow_box=[round(v, 3) for v in af_box],
               inputs_box=inbox, inputs_content=[round(cx0, 3), round(cy0, 3), round(cx1, 3), round(cy1, 3)],
               traces={"left": trL, "right": trR}, tips_x=[round(TIP_L, 3), round(TIP_R, 3)],
               input_arrows=inputs, organ_arrows={"left": oaL, "right": oaR}, panels=panels, panel_cy=round(PANEL_CY, 3), organs=organs,
               labels=labels, label_blue=LABEL_BLUE,
               arrow_lw=ARROW_LW, arrow_color=ARROW_COLOR, threshold_dash=THRESHOLD_DASH, panels_flush=PANELS_FLUSH, org_arrow_mm=ORG_ARROW,
               road=None, badge=None, statements=None, deleted=["brain"],
               moved="round 40: the brain leaves both organ panels (four groups each, ink centred on the panel centre, pitch 25, the round-36 panel rectangles kept); nothing else (dy 0.0 kept)")
    return rec


# ================================================================ Supp Fig 16, top band

def band_s16(final, *, stem="band_s16_r40", outdir=OUT):
    """build_r34.band_s16 with TWO named icons centred under the question. Nothing else."""
    W, H = 496.8 * MM_PER_PT, 79.0 * MM_PER_PT
    fig = new(W, H)
    t = fig.theme
    AXIS = 14.0
    HX = 19.0
    house(fig, HX, 10.0, 20.0, 17.0)
    fig.place("sleeper-faceB", HX, 13.4, w=12)
    l1 = fig.label("SHHS + MrOS", HX, 19.4, pt=t.note_pt, anchor="middle", va="top", bold=True, color=INK)
    l2 = fig.label("5,802 + 2,911 people", HX, 23.0, pt=t.note_pt, anchor="middle", va="top", color=INK)
    TX, TW, TH = 61.0, 19.0, 12.0
    TY = AXIS - TH * 0.5
    tr = trace(fig, TX, TY, TW, TH, color=AZ, dips=(0.15, 0.36, 0.58, 0.80), lw=0.55)
    lbl_left = tr["label_box"][0]
    a1 = arrow(fig, (HX + 10.0 + 1.5, AXIS), (lbl_left - 2.0, AXIS))
    a2 = arrow(fig, (86.0, AXIS), (112.0, AXIS))
    yl = fig.label("years later", 99.0, AXIS + 2.2, pt=t.note_pt, anchor="middle", va="top", color=INK)
    OCX = 143.75
    ttl = fig.label(S16_TITLE.replace("predict future", "predict\nfuture"), OCX, 1.5, pt=10.0, bold=True, anchor="middle", va="top", color=INK21, leading=1.12)
    n0 = len(fig._placed)
    placed, bottom = organ_row_named(fig, OCX, AXIS, pitch=S16_PITCH_R40, icon_h=6.5, max_w=S16_NAME_MAX_W_R40, name_dy=1.5, organs=list(S16_ORGANS_R40))
    icons = fig._placed[n0:]
    assert [p.name for p in icons] == ["heart", "human-grey"], [p.name for p in icons]
    assert abs((icons[0].x + icons[0].w / 2 + icons[1].x + icons[1].w / 2) / 2 - OCX) < 0.01, "the two icons centred under the question"
    assert bottom < H - 0.5, ("names inside the band", bottom, H)
    # the names' line count: "heart attack or failure" as two lines, "death" as one (from the label boxes)
    from bioglyph.qc import label_boxes
    names = {tx: bx for tx, bx in label_boxes(fig) if tx in ("heart attack or failure", "death")}
    assert set(names) == {"heart attack or failure", "death"}, names
    rec = finish(fig, stem, outdir, final=final, attribution=None)
    rec.update(trace=tr, house_cx=HX, axis=AXIS, arrows={"house_to_trace": a1, "years_later": a2}, names_bottom=round(bottom, 3),
               title_box=box(ttl), title_pt=10.0, years_box=box(yl), organs=S16_ORGANS_R40, organs_r34=S16_ORGANS_R34, pitch=S16_PITCH_R40, name_max_w=S16_NAME_MAX_W_R40,
               icon_boxes={p.name: box(p) for p in icons}, name_boxes={k: [round(v, 3) for v in b] for k, b in names.items()},
               death_icon="human-grey", labels={"SHHS + MrOS": box(l1), "people": box(l2)},
               deleted=["heart failure", "heart attack", "stroke", "sv-brain-2 (stroke icon)", "the second heart"], added=["heart attack or failure"],
               moved=f"round 40: the four named icons become two, one heart 'heart attack or failure' and the human figure 'death', pitch {S16_PITCH_R40}, centred under the question; nothing else")
    return rec


SHEETS = {"band_1f_r40": band_1f, "band_3d_r40": band_3d, "band_4e_r40": band_4e, "band_5d_r40": band_5d,
          "band_5d_compact_r40": band_5d_compact, "band_s16_r40": band_s16, "Main_Fig6_r40": fig6}


def sha256(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def whites_gs(pdf, dpi=600):
    """The round-23 measuring stick with Ghostscript: white above the first ink and below the last ink, on real ink at 600 dpi."""
    import numpy as np
    from PIL import Image
    png = Path(pdf).with_suffix(f".{dpi}.png")
    subprocess.run([GS, "-dSAFER", "-dBATCH", "-dNOPAUSE", "-q", "-sDEVICE=png16m", f"-r{dpi}", "-dTextAlphaBits=4", "-dGraphicsAlphaBits=4",
                    f"-sOutputFile={png}", str(pdf)], check=True)
    a = np.asarray(Image.open(png).convert("RGB")).min(axis=2)
    ink = a < 245
    ys = np.nonzero(ink.any(axis=1))[0]
    k = 25.4 / dpi
    return {"top": round(float(ys.min() * k), 3), "bottom": round(float((ink.shape[0] - ys.max() - 1) * k), 3), "page_h_mm": round(ink.shape[0] * k, 3)}


def solve_compact(start):
    SCRATCH.mkdir(parents=True, exist_ok=True)
    v, log = start, []
    for _ in range(5):
        band_5d_compact(True, stem="solve_band_5d_compact", outdir=SCRATCH, shift=v)
        w = whites_gs(SCRATCH / "solve_band_5d_compact.pdf")
        err = w["top"] - w["bottom"]
        log.append({"shift": round(v, 3), "top": w["top"], "bottom": w["bottom"], "err": round(err, 3)})
        print(f"  compact shift={v:.3f}  top {w['top']:.3f}  bottom {w['bottom']:.3f}  err {err:+.3f}", flush=True)
        if abs(err) < 0.03:
            break
        v = v - err / 2
    return round(v, 2), log


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--final", action="store_true")
    ap.add_argument("--preview", action="store_true")
    ap.add_argument("--only", default="")
    ap.add_argument("--solve-compact", action="store_true")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    import fcntl
    lock = open(LOCK, "w")
    fcntl.flock(lock, fcntl.LOCK_EX)                 # every Chrome render of every lane, one at a time
    t0 = time.time()
    try:
        if args.solve_compact:
            v, log = solve_compact(SHIFT_COMPACT)
            (WORK / "solve_compact_r40.json").write_text(json.dumps({"value": v, "start": SHIFT_COMPACT, "log": log}, indent=1))
            print(json.dumps({"band_5d_compact_r40": v}))
        else:
            only = [s for s in args.only.split(",") if s]
            for stem, fn in SHEETS.items():
                if only and stem not in only:
                    continue
                r = fn(args.final)
                r["render"] = dict(RENDER_INFO)
                (WORK / f"records_r40_{stem}.json").write_text(json.dumps(r, indent=1, default=str))
                print(json.dumps({k: r[k] for k in ("stem", "w_pt", "h_mm", "audit", "min_pt", "texts", "numbers")}), flush=True)
    finally:
        fcntl.flock(lock, fcntl.LOCK_UN)
        lock.close()
    print(f"done in {time.time() - t0:.0f} s, render {RENDER_INFO}", flush=True)
