"""Shared drawn elements for the round-12 L1 schematics (FigureLab engine).
R13-C copy (2026-09-03): icons_grey root added, ORGANS_GREY, outputs under R13_schematic/.

Colour law 2026-09-01: oxygen azure, sleep duration green, sleep apnea
terracotta, greys for roads that do not reach disease. Type at theme sizes
only (11 / 10 / 9.5 pt, 9 pt floor). No em-dashes, no semicolons.
"""
from __future__ import annotations
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths

import math
import re
import sys
from pathlib import Path

import fitz
from PIL import Image

sys.path.insert(0, paths.FIGURELAB_DIR)
from bioglyph import Figure                                              # noqa: E402
from bioglyph.qc import audit, label_boxes                               # noqa: E402
from bioglyph.textmetrics import ASCENT, DESCENT, PT_MM, text_width_mm   # noqa: E402
from bioglyph.svg import el                                              # noqa: E402

HERE = Path(__file__).resolve().parent            # R13_schematic/_build/
ICON_DIRS = (HERE / "icons", HERE / "icons_grey")   # grey organ twins (make_grey_icons.py)
ROOT = HERE.parent                                # R13_schematic/
PREVIEWS = HERE / "previews"
MM_PER_PT = 0.352778

# ---------------------------------------------------------------- colour law
# ROUND 18: the azure and the green are pinned to the exact hexes Figure 1's own family key
# prints, so a sketch band and the panel it sits under name the same family in the same
# colour. Azure moves #0288d1 -> #0187d1 and the family green #39c445 -> #57cc60.
AZ = "#0187d1"                                    # oxygen / T90 (Fig 1 key: oxygenation)
AZ_T = ["#b3dcf2", "#7cc0e9", "#3f9fd8", "#0187d1"]
AZ_D = "#02598a"
GREEN = "#39c445"                                 # sleep duration (round-20 colour law)
GREEN_L = ["#c8eecb", "#92df99", "#5dcf66", "#298d32"]         # the ladder, light to dark
GREEN_D = "#298d32"                               # the dark twin, used for green TEXT only,
#                                                   #39c445 is too light to set type in
TERRA, TERRA_T = "#b5623a", "#d99b76"             # NOT sleep apnea any more (round 15)
# ROUND 15 colour law: sleep apnea leaves terracotta, which Alen says sits too close to
# the hazard-ratio red, and becomes purple. Purple is the only free hue: blue is oxygen,
# green is duration, red is the hazard ramp. Terracotta stays available for anything that
# was terracotta for a different reason, e.g. the breathing channel of a PSG montage.
PURPLE, PURPLE_T = "#7e57c2", "#c5b3e6"           # NOT sleep apnea any more (round 16)
# ROUND 16 colour law, settled by Alen and not to be re-litigated. In Figure 1 the family
# called "Breathing events" IS the apnea family, so it must be recoloured, and it cannot be
# purple because purple already means "Heart rate and variability" in the same key.
# SLEEP APNEA AND BREATHING EVENTS ARE DARK TEAL. Round 15's purple is replaced by this teal
# everywhere it stood for apnea, and so is the terracotta that stood for apnea before it.
TEAL, TEAL_T, TEAL_P = "#00695c", "#4d9b8f", "#b2dfdb"    # NOT sleep apnea any more (R18)
ORANGE, ORANGE_T, ORANGE_P = "#ff8f00", "#ffb74d", "#ffe0b2"   # NOT sleep apnea any more (R20)
# ROUND 20 colour law, settled by Alen and not to be re-litigated. Sleep apnea and breathing
# events leave round 18's orange and become CHARCOAL, on a four-rung ladder so a severity
# range can be drawn in the family itself. Type is set in the dark twin, never in a rung.
CHARCOAL = "#37474f"                              # sleep apnea, breathing events
CHARCOAL_L = ["#d3d7d8", "#a5acb0", "#737e84", "#37474f"]      # the ladder, light to dark
CHARCOAL_D = "#29353b"                            # the dark twin, for apnea TEXT
APNEA, APNEA_L, APNEA_D = CHARCOAL, CHARCOAL_L, CHARCOAL_D      # lane LE reads it by meaning
APNEA_T, APNEA_P = CHARCOAL_L[2], CHARCOAL_L[0]
GREY, GREY_L = "#8a9099", "#ccd1d6"
INK, WHITE, CARD = "#1a1d20", "#ffffff", "#f4f6f8"

ORGANS = [("heart", "heart failure"), ("sv-lung", "respiratory failure"),
          ("pancreas", "type 2 diabetes"), ("liver", "cirrhosis")]
ORGANS_GREY = [(f"{icon}-grey", name) for icon, name in ORGANS]   # organs under damage
GREY_D = "#3c4046"                                # dark-grey organ ink


def new(w: float, h: float) -> Figure:
    fig = Figure(w, h, style="real", icon_dirs=ICON_DIRS)
    fig._vboxes = []
    return fig


# ---------------------------------------------------------------- traces

def spo2_trace(fig, x, y, w, h, *, dips: bool, color: str, ref_at=0.5,
               dip_centers=(0.12, 0.29, 0.46, 0.62, 0.79, 0.93), label90="left",
               lw=0.55, label_pt=None, x_start=None,
               ref_lw=0.3, ref_dash="1.4 1.0", label_bold=False, label_gap=1.2):
    # ROUND 29 (lane LE): ref_lw / ref_dash weight the dashed 90 percent reference line,
    # label_bold sets the "90%" tick label bold, label_gap is its gap to the line's start.
    # Every default is the round-27 value, so every other sheet renders byte-identically.
    """Overnight blood-oxygen line with a dashed 90% reference. With dips the
    time spent below 90% is shaded in the line colour. x_start lets the line
    begin later than the reference line (used for the awake/asleep split)."""
    y90 = y + h * ref_at
    n = 220
    pts = []
    x0 = x if x_start is None else x_start
    w0 = w - (x0 - x)
    for i in range(n + 1):
        t = i / n
        if not dips:
            py = y + h * 0.26 + h * 0.05 * math.sin(t * 2 * math.pi * 3.3) \
                + h * 0.03 * math.sin(t * 2 * math.pi * 9.7 + 1)
        else:
            py = y + h * 0.24 + h * 0.03 * math.sin(t * 2 * math.pi * 7)
            for c in dip_centers:
                py += h * 0.68 * math.exp(-((t - c) / 0.035) ** 2)
        pts.append((x0 + w0 * t, py))
    fig.line(x, y90, x + w, y90, stroke=GREY, lw=ref_lw, dash=ref_dash)
    lp = label_pt or fig.theme.note_pt
    if label90 == "left":
        fig.label("90%", x - label_gap, y90, pt=lp, color=GREY, anchor="end", va="middle",
                  bold=label_bold)
    elif label90 == "right":
        fig.label("90%", x + w + label_gap, y90, pt=lp, color=GREY, va="middle",
                  bold=label_bold)
    if dips:
        seg = []
        for px, py in pts + [(x + w, y90 - 1)]:
            if py > y90:
                seg.append((px, py))
            elif seg:
                fig.polygon([(seg[0][0], y90)] + seg + [(seg[-1][0], y90)], fill=color, opacity=0.35)
                seg = []
    fig.path("M" + " L".join(f"{px:.2f} {py:.2f}" for px, py in pts), stroke=color, lw=lw)
    return y90


def night_axis(fig, x, y, w, *, text="night", color=GREY, pt=None, bold=False):
    """Thin baseline with the word 'night' under a trace. ROUND 29: bold= for the word."""
    fig.line(x, y, x + w, y, stroke=GREY, lw=0.3)
    return fig.label(text, x + w / 2, y + 1.0, pt=pt or fig.theme.note_pt, color=color,
                     anchor="middle", va="top", bold=bold)


def airflow(fig, cx, cy, w, h, color=TERRA, lw=0.55):
    """Breathing airflow with a pause in the middle: the sleep-apnea signature."""
    x0, n, pts = cx - w / 2, 240, []
    for i in range(n + 1):
        t = i / n
        if 0.36 < t < 0.66:
            py = cy
        elif t <= 0.36:
            py = cy - (h / 2) * math.sin(2 * math.pi * (t / 0.36) * 2.5)
        else:
            py = cy - (h / 2) * 0.85 * math.sin(2 * math.pi * ((t - 0.66) / 0.34) * 2.5)
        pts.append((x0 + w * t, py))
    fig.path("M" + " L".join(f"{px:.2f} {py:.2f}" for px, py in pts), stroke=color, lw=lw)


def swatch_label(fig, x, y, text, *, color=AZ, pt=None, anchor="start"):
    """Small shaded square plus its meaning, e.g. 'time below 90%'. With
    anchor='middle', x is the centre of the swatch+text pair."""
    pt = pt or fig.theme.note_pt
    tw = text_width_mm(text, pt, False)
    if anchor == "middle":
        x = x - (3.2 + 1.4 + tw) / 2
    fig.rect(x, y - 1.6, 3.2, 3.2, fill=color, opacity=0.35)
    fig.rect(x, y - 1.6, 3.2, 3.2, stroke=color, lw=0.3)
    return fig.label(text, x + 4.6, y, pt=pt, color=color, va="middle")


# ---------------------------------------------------------------- icons drawn by hand

def house(fig, cx, cy, w, h, *, stroke=GREY, roof=GREY_L):
    """Simple house outline: roof, chimney, body, door. (cx, cy) is the centre
    of the full w x h box. Returns the body box (x, y, w, h) for placing things inside."""
    top, bot = cy - h / 2, cy + h / 2
    roof_b = top + h * 0.36
    fig.rect(cx + w * 0.22, top + h * 0.10, w * 0.09, h * 0.2, fill=roof, stroke=stroke, lw=0.35)
    body = (cx - w * 0.42, roof_b - 0.4, w * 0.84, bot - roof_b + 0.4)
    fig.rect(*body, rx=0.8, fill=WHITE, stroke=stroke, lw=0.45)
    fig.polygon([(cx - w / 2, roof_b), (cx, top), (cx + w / 2, roof_b)],
                fill=roof, stroke=stroke, lw=0.45)
    return body


def rank_ladder(fig, x, y, w, h, *, rungs=4, top_color=AZ, other=AZ_T[0], stroke=GREY):
    """Tiny ranking ladder: two rails, `rungs` rungs, the top rung filled azure."""
    fig.line(x, y, x, y + h, stroke=stroke, lw=0.35)
    fig.line(x + w, y, x + w, y + h, stroke=stroke, lw=0.35)
    step = h / (rungs - 1) if rungs > 1 else h
    for i in range(rungs):
        ry = y + i * step
        fig.line(x, ry, x + w, ry, stroke=top_color if i == 0 else stroke, lw=0.9 if i == 0 else 0.35)
    fig.circle(x + w / 2, y, w * 0.28, fill=top_color)


# ---------------------------------------------------------------- text helpers

def vlabel(fig, text, x, y, *, pt=None, color=INK, bold=False):
    """Label rotated 90 degrees counter-clockwise, centred on (x, y)."""
    pt = fig.theme.label_pt if pt is None else pt
    size = pt * PT_MM
    wmm = text_width_mm(text, pt, bold)
    bx = x + (ASCENT - DESCENT) * size / 2
    fig.add(el("text", {"x": -y, "y": bx, "font_size": size, "fill": color,
                        "font_weight": "bold" if bold else None, "text_anchor": "middle",
                        "transform": "rotate(-90)"}, text))
    fig._fonts_pt.append(pt)
    box = (x - (ASCENT + DESCENT) * size / 2, y - wmm / 2, x + (ASCENT + DESCENT) * size / 2, y + wmm / 2)
    fig._vboxes.append((text, box))
    return box


def organ_row_named(fig, cx, cy, *, pitch, icon_h, name_dy=1.8, max_w=None, organs=ORGANS,
                    name_pt=None, name_color=INK):
    """Heart, lung, pancreas, liver on one optical baseline, the disease name
    written under each icon. Returns (placed icons, y where the names end)."""
    x0 = cx - pitch * (len(organs) - 1) / 2
    bottom = cy + icon_h / 2 + name_dy
    placed = []
    for i, (icon, name) in enumerate(organs):
        x = x0 + pitch * i
        placed.append(fig.place(icon, x, cy, h=icon_h))
        lb = fig.label(name, x, cy + icon_h / 2 + name_dy, pt=name_pt or fig.theme.note_pt,
                       color=name_color, anchor="middle", va="top",
                       max_w=max_w or pitch - 1.0, leading=1.15)
        bottom = max(bottom, lb.y + lb.h)
    return placed, bottom


def answer(fig, cx, y, text, *, sign, color, text_color=INK, bold=True, r=2.6, pt=None,
           anchor="middle"):
    """A badge and its one-line answer. anchor='middle' centres the pair on cx,
    anchor='start' puts the badge's left edge at cx."""
    pt = pt or fig.theme.label_pt
    tw = text_width_mm(text, pt, bold)
    total = 2 * r + 2.2 + tw
    x = cx - total / 2 if anchor == "middle" else cx
    fig.badge(x + r, y, r, sign=sign, color=color, glyph_color=WHITE)
    lb = fig.label(text, x + 2 * r + 2.2, y + 0.1, pt=pt, bold=bold, color=text_color, va="middle")
    return x, x + total, lb


def crossed_road(fig, p_from, p_to, *, badge_at=None, label=None, label_dy=3.6,
                 label_dx=4.0, label_va="top", color=GREY, r=2.4, pt=None):
    """Grey dashed elbow road (horizontal first) that does not reach disease on
    its own: a cross badge sits on the horizontal leg, the reason under it."""
    fig.arrow(p_from, p_to, color=color, lw=0.4, dash="1.6 1.2", elbow="h")
    bx = badge_at if badge_at is not None else (p_from[0] + p_to[0]) / 2
    fig.badge(bx, p_from[1], r, sign="×", color=color, glyph_color=WHITE)
    lb = None
    if label:
        lb = fig.label(label, bx + label_dx, p_from[1] + label_dy, pt=pt or fig.theme.note_pt,
                       anchor="middle", va=label_va, color=color)
    return lb


def number_circle(fig, n, cx, cy, r=3.3, color=AZ):
    fig.circle(cx, cy, r, fill=color)
    fig.label(str(n), cx, cy + 0.15, pt=fig.theme.label_pt, bold=True, color=WHITE,
              anchor="middle", va="middle")


def numbered_caption(fig, n, x, y, text, *, r=3.3, circle_color=AZ, color=AZ, pt=None,
                     gap=2.4, bold=True, anchor="start", max_w=None, leading=1.15):
    """Round 16: the numbered azure circle carries its caption BESIDE it, on one line, in
    the same azure, at the band's normal label size, so the numbers read as a titled
    sequence across the top of the band instead of as bare digits over a picture.

    anchor='start'  the pair's LEFT EDGE sits at x (the circle's centre is x + r)
    anchor='middle' the circle-plus-caption pair is centred on x
    Returns (x_left, x_right, label) so callers can assert the pair's own box.
    """
    pt = pt or fig.theme.label_pt
    tw = text_width_mm(text, pt, bold)
    total = 2 * r + gap + tw
    cx = x + r if anchor == "start" else x - total / 2 + r
    number_circle(fig, n, cx, y, r=r, color=circle_color)
    lb = fig.label(text, cx + r + gap, y + 0.1, pt=pt, bold=bold, color=color, va="middle",
                   max_w=max_w, leading=leading)
    return (cx - r, cx - r + total, lb)


# ---------------------------------------------------------------- QC and export

def _vbox_collisions(fig, tol=0.6):
    items = [("icon", p.name, (p.x, p.y, p.x + p.w, p.y + p.h)) for p in fig._placed]
    items += [("label", t, b) for t, b in label_boxes(fig)]
    out = []
    for text, vb in fig._vboxes:
        for kind, name, b in items:
            ox = min(vb[2], b[2]) - max(vb[0], b[0])
            oy = min(vb[3], b[3]) - max(vb[1], b[1])
            if ox > tol and oy > tol:
                out.append(f"rotated label '{text}' overlaps {kind} '{name}' by {min(ox, oy):.1f} mm")
    return out


def _bounds_problems(fig):
    """Labels that hang off the canvas (the engine only checks icons)."""
    out = []
    for text, (x0, y0, x1, y1) in label_boxes(fig):
        if x0 < -0.3 or y0 < -0.3 or x1 > fig.w + 0.3 or y1 > fig.h + 0.3:
            out.append(f"label '{text[:30]}' beyond canvas ({x0:.1f},{y0:.1f})-({x1:.1f},{y1:.1f})")
    return out


def exact_boxes(pdf_path: Path, png_path: Path, w_mm: float, h_mm: float, dpi: int) -> None:
    """Chromium pads the PDF page box by up to ~1.2 pt but draws the content at the
    requested size anchored top-left (measured 2026-09-03). Trim the MediaBox in PDF
    space (bottom strip) and set the CropBox so both boxes are exactly w x h pt with
    the content untouched. The PNG carries the same blank padding: crop it to the
    exact pixel size and stamp the dpi."""
    w_pt, h_pt = w_mm / MM_PER_PT, h_mm / MM_PER_PT
    if pdf_path.exists():
        doc = fitz.open(pdf_path)
        page = doc[0]
        mb = page.mediabox
        page.set_mediabox(fitz.Rect(0, mb.y1 - h_pt, w_pt, mb.y1))   # also resets the CropBox
        got = page.rect
        assert abs(got.width - w_pt) < 0.02 and abs(got.height - h_pt) < 0.02, (got, w_pt, h_pt)
        tmp = pdf_path.with_suffix(".tmp.pdf")
        doc.save(tmp, garbage=3, deflate=True)
        doc.close()
        tmp.replace(pdf_path)
    if png_path.exists():
        im = Image.open(png_path)
        W, H = round(w_mm / 25.4 * dpi), round(h_mm / 25.4 * dpi)
        im = im.crop((0, 0, min(W, im.width), min(H, im.height)))
        im.save(png_path, dpi=(dpi, dpi))


_NUM = re.compile(r"\d[\d,\.]*x?")


def finish(fig, stem, outdir, *, final: bool, attribution="auto", dpi=300):
    """Audit, print, export. Returns a record for the report."""
    probs = audit(fig) + _vbox_collisions(fig) + _bounds_problems(fig)
    texts = [t for t, _ in label_boxes(fig)]
    texts = [t for t in texts if not t.startswith("Icons:")]
    nums = sorted({m.group(0).rstrip(".") for t in texts for m in _NUM.finditer(t)
                   if m.group(0).rstrip(".") not in ("90", "2")})
    print(f"[{stem}] {fig.w:.1f} x {fig.h:.1f} mm  audit: {probs or 'clean'}")
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    PREVIEWS.mkdir(parents=True, exist_ok=True)
    if final:
        fig.export(outdir / f"{stem}.svg", outdir / f"{stem}.png", outdir / f"{stem}.pdf",
                   dpi=dpi, attribution=attribution)
        exact_boxes(outdir / f"{stem}.pdf", outdir / f"{stem}.png", fig.w, fig.h, dpi)
    fig.export(PREVIEWS / f"{stem}_preview.png", dpi=150, attribution=attribution)
    icons = [{"name": p.name, "cx": round(p.x + p.w / 2, 2), "cy": round(p.y + p.h / 2, 2),
              "w": round(p.w, 2), "h": round(p.h, 2)} for p in fig._placed]
    return {"stem": stem, "w_mm": round(fig.w, 1), "h_mm": round(fig.h, 1),
            "w_pt": round(fig.w / MM_PER_PT, 1), "audit": probs, "texts": texts,
            "numbers": nums, "min_pt": min(fig._fonts_pt) if fig._fonts_pt else None,
            "icons": icons, "electrodes": getattr(fig, "_eeg", []),
            "attribution": fig.attribution() if attribution else "(credit in legend)"}
