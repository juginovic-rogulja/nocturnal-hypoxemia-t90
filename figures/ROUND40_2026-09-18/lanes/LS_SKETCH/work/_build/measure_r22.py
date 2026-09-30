#!/usr/bin/env python3
"""Round 22, lane LA: the measuring stick. Every clearance the report quotes is read off a
rendered page here, at 600 dpi, from real ink. Nothing in the report is a number I typed.

  ink_box(pdf)             the ink's own bounding box, in mm
  top_white / bottom_white the white above the first ink and below the last ink
  zone_clear(pdf, box)     is a rectangle of the page free of ink
  corner_profile(pdf)      the first ink's y inside a series of left-hand strips
"""
from __future__ import annotations

import sys
from pathlib import Path

import fitz
import numpy as np

DPI = 600
K = 25.4 / DPI                                    # mm per pixel at 600 dpi
THRESH = 245                                      # anything not near-white is ink


def page_ink(pdf, dpi=DPI):
    d = fitz.open(str(pdf))
    p = d[0]
    pix = p.get_pixmap(dpi=dpi)
    a = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, pix.n)
    g = a[:, :, :3].min(axis=2) if pix.n >= 3 else a[:, :, 0]
    return (g < THRESH), (25.4 / dpi), p.rect


def ink_box(pdf, dpi=DPI):
    ink, k, rect = page_ink(pdf, dpi)
    ys, xs = np.nonzero(ink)
    return {"x0": xs.min() * k, "x1": (xs.max() + 1) * k,
            "y0": ys.min() * k, "y1": (ys.max() + 1) * k,
            "page_w_mm": rect.width * 25.4 / 72, "page_h_mm": rect.height * 25.4 / 72,
            "page_w_pt": rect.width, "page_h_pt": rect.height}


def whites(pdf, dpi=DPI):
    b = ink_box(pdf, dpi)
    return {"top": b["y0"], "bottom": b["page_h_mm"] - b["y1"],
            "left": b["x0"], "right": b["page_w_mm"] - b["x1"]}


def zone_clear(pdf, x0, y0, x1, y1, dpi=DPI):
    """Is the rectangle (mm) free of ink? Returns (clear, n_ink_pixels, first_ink_mm)."""
    ink, k, _ = page_ink(pdf, dpi)
    c0, c1 = int(round(x0 / k)), int(round(x1 / k))
    r0, r1 = int(round(y0 / k)), int(round(y1 / k))
    sub = ink[r0:r1, c0:c1]
    n = int(sub.sum())
    if n == 0:
        return True, 0, None
    ys, xs = np.nonzero(sub)
    return False, n, ((c0 + xs.min()) * k, (r0 + ys.min()) * k)


def corner_profile(pdf, widths=(4, 6, 8, 10, 12, 14), dpi=DPI):
    """For each left-hand strip of the given width in mm, the y of its first ink."""
    ink, k, _ = page_ink(pdf, dpi)
    out = {}
    for w in widths:
        sub = ink[:, :int(round(w / k))]
        ys, _ = np.nonzero(sub)
        out[w] = float(ys.min() * k) if len(ys) else None
    return out


def rows_with_ink(pdf, dpi=DPI):
    """The y of every ink row, so gaps between blocks can be found."""
    ink, k, _ = page_ink(pdf, dpi)
    return np.nonzero(ink.any(axis=1))[0] * k


def blocks(pdf, min_gap_mm=1.0, dpi=DPI):
    """The page's horizontal ink bands, split wherever there are min_gap_mm of white."""
    ys = rows_with_ink(pdf, dpi)
    if not len(ys):
        return []
    out, start, prev = [], ys[0], ys[0]
    for y in ys[1:]:
        if y - prev > min_gap_mm:
            out.append((float(start), float(prev)))
            start = y
        prev = y
    out.append((float(start), float(prev)))
    return out


if __name__ == "__main__":
    for f in sys.argv[1:]:
        b = ink_box(f)
        w = whites(f)
        print(Path(f).name)
        print(f"  page  {b['page_w_pt']:.4f} x {b['page_h_pt']:.4f} pt"
              f"   = {b['page_w_mm']:.3f} x {b['page_h_mm']:.3f} mm")
        print(f"  ink   x {b['x0']:.3f} .. {b['x1']:.3f}   y {b['y0']:.3f} .. {b['y1']:.3f}")
        print(f"  white top {w['top']:.3f}  bottom {w['bottom']:.3f}  "
              f"left {w['left']:.3f}  right {w['right']:.3f}")
        print(f"  corner {corner_profile(f)}")
