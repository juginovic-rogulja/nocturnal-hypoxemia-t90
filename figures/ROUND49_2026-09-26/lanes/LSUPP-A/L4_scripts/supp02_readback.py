"""
Read a Supp_Fig02 sheet's drawn values back from its PDF drawings (the base sheet or the new one).

The axes are calibrated from their own spines and tick marks (the tick marks' data values are the
known tick sets), then every marker centre, whisker end, median-line vertex and interquartile-polygon
vertex is converted to data units. Nothing here reads a numbers file: the result is what the sheet
shows, to be compared with the numbers file by the caller.

    import readback as RB
    r = RB.read_sheet(pdf, yticks=[88, 90, 92, 94, 96, 98])
"""
import gc

import fitz
import numpy as np

INK = (0x1a / 255, 0x1d / 255, 0x21 / 255)
LADDER = ["#b3dcf2", "#7cc0e9", "#3f9fd8", "#0288d1"]        # the base sheet's key handles, left to right
BAND_KEYS = ["le1", "1to5", "5to10", "gt10"]
OFFS = {"le1": -0.27, "1to5": -0.09, "5to10": 0.09, "gt10": 0.27}   # panel b x offsets of the four bands
STAGES = ["Wake", "N1", "N2", "N3", "REM"]
MS_A, MS_B = 4.2, 6.0                                        # marker diameters, panel a and panel b
LW_MEDIAN, LW_WHISKER, LW_AXIS, TICK_LEN = 1.6, 1.7, 0.8, 3.0


def hexc(c):
    return "#%02x%02x%02x" % tuple(int(round(v * 255)) for v in c)


def near(c, ref, tol=0.02):
    return c is not None and all(abs(a - b) <= tol for a, b in zip(c, ref))


def band_of(c):
    h = hexc(c) if c is not None else None
    return BAND_KEYS[LADDER.index(h)] if h in LADDER else None


def drawings(pdf):
    doc = fitz.open(pdf)
    page = doc[0]
    dr = page.get_drawings()
    rect = (page.rect.width, page.rect.height)
    doc.close()
    gc.collect()
    return dr, rect


def ink_segments(dr):
    """Unique single-segment ink strokes of the axis width: spines and tick marks."""
    out = set()
    for d in dr:
        if d["type"] not in ("s", "fs") or len(d["items"]) != 1 or d["items"][0][0] != "l":
            continue
        if not near(d.get("color"), INK) or abs((d.get("width") or 0) - LW_AXIS) > 0.01:
            continue
        p0, p1 = d["items"][0][1], d["items"][0][2]
        a, b = (p0.x, p0.y), (p1.x, p1.y)
        if a > b:
            a, b = b, a
        out.add((round(a[0], 3), round(a[1], 3), round(b[0], 3), round(b[1], 3)))
    return sorted(out)


class AxesBox:
    """An axes box in page coordinates (fitz: origin top left, y down)."""
    def __init__(self, x0, y0, x1, y1):
        self.x0, self.y0, self.x1, self.y1 = x0, y0, x1, y1

    def contains(self, x, y, pad=1.0):
        return self.x0 - pad <= x <= self.x1 + pad and self.y0 - pad <= y <= self.y1 + pad

    def as_list(self):
        return [round(self.x0, 3), round(self.y0, 3), round(self.x1, 3), round(self.y1, 3)]


def find_axes(segs):
    lefts = sorted([s for s in segs if abs(s[0] - s[2]) < 1e-3 and s[3] - s[1] > 50], key=lambda s: s[0])
    bottoms = sorted([s for s in segs if abs(s[1] - s[3]) < 1e-3 and s[2] - s[0] > 50], key=lambda s: s[0])
    assert len(lefts) == 2 and len(bottoms) == 2, (lefts, bottoms)
    axes = []
    for l, b in zip(lefts, bottoms):
        assert abs(l[0] - b[0]) < 0.01 and abs(l[3] - b[1]) < 0.01, (l, b)
        axes.append(AxesBox(l[0], l[1], b[2], l[3]))
    return axes


def calibrate(segs, ax, xticks, yticks):
    """Linear page-to-data maps from the outward tick marks against the known tick sets."""
    xt = sorted(s[0] for s in segs if abs(s[0] - s[2]) < 1e-3 and abs((s[3] - s[1]) - TICK_LEN) < 0.05
                and abs(s[1] - ax.y1) < 0.05 and ax.x0 - 1 <= s[0] <= ax.x1 + 1)
    yt = sorted((s[1] for s in segs if abs(s[1] - s[3]) < 1e-3 and abs((s[2] - s[0]) - TICK_LEN) < 0.05
                 and abs(s[2] - ax.x0) < 0.05 and ax.y0 - 1 <= s[1] <= ax.y1 + 1), reverse=True)
    assert len(xt) == len(xticks), (xt, xticks)
    assert len(yt) == len(yticks), (yt, yticks)
    px = np.polyfit(xt, xticks, 1)
    py = np.polyfit(yt, yticks, 1)
    rx = float(np.max(np.abs(np.polyval(px, xt) - np.asarray(xticks, float))))
    ry = float(np.max(np.abs(np.polyval(py, yt) - np.asarray(yticks, float))))
    return dict(x=lambda v: float(np.polyval(px, v)), y=lambda v: float(np.polyval(py, v)),
                x_ticks_page=[round(v, 3) for v in xt], y_ticks_page=[round(v, 3) for v in yt],
                fit_residual_data=(rx, ry),
                xlim=(float(np.polyval(px, ax.x0)), float(np.polyval(px, ax.x1))),
                ylim=(float(np.polyval(py, ax.y1)), float(np.polyval(py, ax.y0))))


def read_sheet(pdf, yticks, xticks_a=tuple(range(1, 11)), xticks_b=tuple(range(5))):
    dr, page = drawings(pdf)
    segs = ink_segments(dr)
    axa, axb = find_axes(segs)
    ca = calibrate(segs, axa, list(xticks_a), list(yticks))
    cb = calibrate(segs, axb, list(xticks_b), list(yticks))

    a_marker = {k: [None] * 10 for k in BAND_KEYS}
    b_marker = {k: [None] * 5 for k in BAND_KEYS}
    a_line = {k: [None] * 10 for k in BAND_KEYS}
    a_poly = {k: [[None, None] for _ in range(10)] for k in BAND_KEYS}
    b_whisker = {k: [[None, None] for _ in range(5)] for k in BAND_KEYS}
    legend = []
    stray = []

    for d in dr:
        t = d["type"]
        # markers: filled circles with a white edge in a band colour
        if t == "fs" and near(d.get("color"), (1.0, 1.0, 1.0)) and band_of(d.get("fill")) and len(d["items"]) >= 4:
            k = band_of(d["fill"])
            r = d["rect"]
            cx, cy = (r.x0 + r.x1) / 2, (r.y0 + r.y1) / 2
            if axa.contains(cx, cy) and abs(r.width - MS_A) < 0.3:
                xd, yd = ca["x"](cx), ca["y"](cy)
                i = int(round(xd))
                assert abs(xd - i) < 0.02 and 1 <= i <= 10, (k, xd)
                a_marker[k][i - 1] = yd
            elif axb.contains(cx, cy) and abs(r.width - MS_B) < 0.3:
                xd, yd = cb["x"](cx), cb["y"](cy)
                j = int(round(xd - OFFS[k]))
                assert abs(xd - OFFS[k] - j) < 0.02 and 0 <= j <= 4, (k, xd)
                b_marker[k][j] = yd
            elif cy < axa.y0:
                legend.append(("marker", k, round(cx, 3), round(cy, 3)))
            else:
                stray.append(("marker", k, round(cx, 3), round(cy, 3)))
            continue
        # whiskers: vertical strokes of the whisker width in a band colour
        if t == "s" and len(d["items"]) == 1 and d["items"][0][0] == "l" and abs((d.get("width") or 0) - LW_WHISKER) < 0.01 and band_of(d.get("color")):
            k = band_of(d["color"])
            p0, p1 = d["items"][0][1], d["items"][0][2]
            assert abs(p0.x - p1.x) < 1e-3, (k, p0, p1)
            xd = cb["x"](p0.x)
            j = int(round(xd - OFFS[k]))
            assert abs(xd - OFFS[k] - j) < 0.02 and 0 <= j <= 4, (k, xd)
            lo, hi = sorted([cb["y"](p0.y), cb["y"](p1.y)])
            b_whisker[k][j] = [lo, hi]
            continue
        # median lines: nine segments of the median width in a band colour inside panel a. The key handles
        # carry the same width and colour above the axes (one segment on the base sheet, two on a newer
        # matplotlib), so they are classed by position, not by segment count.
        if t == "s" and abs((d.get("width") or 0) - LW_MEDIAN) < 0.01 and band_of(d.get("color")):
            k = band_of(d["color"])
            r = d["rect"]
            if r.y1 < axa.y0:
                legend.append(("handle", k, round(r.x0, 3), round((r.y0 + r.y1) / 2, 3)))
                continue
            assert len(d["items"]) == 9 and all(it[0] == "l" for it in d["items"]), (k, len(d["items"]))
            pts = [d["items"][0][1]] + [it[2] for it in d["items"]]
            for p in pts:
                xd = ca["x"](p.x)
                i = int(round(xd))
                assert abs(xd - i) < 0.02 and 1 <= i <= 10, (k, xd)
                a_line[k][i - 1] = ca["y"](p.y)
            continue
        # interquartile polygons: translucent fills in a band colour
        if t == "f" and d.get("fill_opacity") is not None and d["fill_opacity"] < 0.5 and band_of(d.get("fill")):
            k = band_of(d["fill"])
            byx = {}
            for it in d["items"]:
                assert it[0] == "l", it[0]
                for p in (it[1], it[2]):
                    xd = ca["x"](p.x)
                    i = int(round(xd))
                    assert abs(xd - i) < 0.02 and 1 <= i <= 10, (k, xd)
                    byx.setdefault(i, []).append(ca["y"](p.y))
            assert sorted(byx) == list(range(1, 11)), sorted(byx)
            for i, ys in byx.items():
                a_poly[k][i - 1] = [min(ys), max(ys)]
            continue

    handles = sorted([h for h in legend if h[0] == "handle"], key=lambda h: h[2])
    ladder_order = []
    for h in handles:
        if h[1] not in ladder_order:
            ladder_order.append(h[1])
    return dict(pdf=pdf, page=[round(page[0], 3), round(page[1], 3)],
                axes_a=axa.as_list(), axes_b=axb.as_list(),
                calib_a={kk: vv for kk, vv in ca.items() if kk not in ("x", "y")},
                calib_b={kk: vv for kk, vv in cb.items() if kk not in ("x", "y")},
                a_marker=a_marker, a_line=a_line, a_poly=a_poly, b_marker=b_marker, b_whisker=b_whisker,
                legend_items=legend, legend_ladder_order=ladder_order, stray=stray, n_drawings=len(dr))
