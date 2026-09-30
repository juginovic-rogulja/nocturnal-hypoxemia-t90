"""PDF-level overlay helpers, the round-29b idiom (fitz Shape + TextWriter), for lane V7_L6_HABITUAL_EXTERNAL.

A data mark is replaced by (1) a white band painted over the old mark, (2) the rules that crossed the band redrawn
inside it (faint vertical grid, the dashed HR = 1 reference with its dash phase continued from the sheet's own line,
horizontal rules where a band touches one), (3) the new CI line and marker drawn with the base sheet's own geometry
(widths, caps, marker areas, colours read from the base drawing, never typed). A printed value is replaced by a
text-only redaction of the old span's box (no fill, graphics untouched) and a fitz.TextWriter stamp of the new string
with the sheet's Arial face at the measured size, colour and alignment. Nothing else on the sheet is touched.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import math
import fitz

ARIAL = paths.FONT_ARIAL
ARIALB = paths.FONT_ARIAL_BOLD
_FONTS = {}


def font(bold=False):
    key = "b" if bold else "r"
    if key not in _FONTS:
        _FONTS[key] = fitz.Font(fontfile=ARIALB if bold else ARIAL)
    return _FONTS[key]


def hexrgb(h):
    """A colour as an RGB tuple: a '#rrggbb' string, or a tuple passed through unchanged (the base drawing's own
    floats, so a re-drawn rule renders to the identical pixel values)."""
    if isinstance(h, (tuple, list)):
        return tuple(float(v) for v in h)
    return tuple(int(h[i:i + 2], 16) / 255 for i in (1, 3, 5))


def hx(c):
    return "#%02x%02x%02x" % tuple(int(round(v * 255)) for v in c)


class LogX:
    """x = a + b ln(v), fitted from the sheet's own tick positions {value: x}."""
    def __init__(self, ticks):
        vs = sorted(ticks); xs = [ticks[v] for v in vs]
        import numpy as np
        b, a = np.polyfit([math.log(v) for v in vs], xs, 1)
        self.a, self.b = float(a), float(b)
        self.resid = max(abs(a + b * math.log(v) - ticks[v]) for v in vs)
    def x(self, v):
        return self.a + self.b * math.log(v)
    def v(self, x):
        return math.exp((x - self.a) / self.b)


class LogY:
    """y = a + b ln(v) (page y grows downward, so b is negative), from {value: y}."""
    def __init__(self, ticks):
        vs = sorted(ticks); ys = [ticks[v] for v in vs]
        import numpy as np
        b, a = np.polyfit([math.log(v) for v in vs], ys, 1)
        self.a, self.b = float(a), float(b)
        self.resid = max(abs(a + b * math.log(v) - ticks[v]) for v in vs)
    def y(self, v):
        return self.a + self.b * math.log(v)
    def v(self, y):
        return math.exp((y - self.a) / self.b)


def snap_band(band, rule_xs, pad=1.2):
    """Widen a band so that no vertical rule is split by a band edge (a split rule would be half whited and then
    redrawn in full, double-painting its uncovered half). Any rule within pad of an edge is taken inside."""
    x0, y0, x1, y1 = band
    for rx in rule_xs:
        if x0 - pad <= rx <= x0 + pad: x0 = rx - pad
        if x1 - pad <= rx <= x1 + pad: x1 = rx + pad
    return (x0, y0, x1, y1)


def white_rect(page, rect):
    sh = page.new_shape(); sh.draw_rect(fitz.Rect(*rect)); sh.finish(fill=(1, 1, 1), color=None); sh.commit()


def fill_rect(page, rect, fill_hex):
    sh = page.new_shape(); sh.draw_rect(fitz.Rect(*rect)); sh.finish(fill=hexrgb(fill_hex), color=None); sh.commit()


def vline(page, x, y0, y1, color_hex, width, cap=0, dashes=None, phase=0.0):
    """A vertical rule from y0 to y1 (drawn in that direction, so a dash phase measured from y0 applies)."""
    sh = page.new_shape(); sh.draw_line(fitz.Point(x, y0), fitz.Point(x, y1))
    kw = dict(color=hexrgb(color_hex), width=width, lineCap=cap, fill=None, closePath=False)
    if dashes:
        kw["dashes"] = f"[{dashes[0]:.4f} {dashes[1]:.4f}] {phase:.4f}"
    sh.finish(**kw); sh.commit()


def hline(page, x0, x1, y, color_hex, width, cap=0, dashes=None, phase=0.0):
    sh = page.new_shape(); sh.draw_line(fitz.Point(x0, y), fitz.Point(x1, y))
    kw = dict(color=hexrgb(color_hex), width=width, lineCap=cap, fill=None, closePath=False)
    if dashes:
        kw["dashes"] = f"[{dashes[0]:.4f} {dashes[1]:.4f}] {phase:.4f}"
    sh.finish(**kw); sh.commit()


def polyline(page, pts, color_hex, width, cap=0):
    sh = page.new_shape(); sh.draw_polyline([fitz.Point(*p) for p in pts]); sh.finish(color=hexrgb(color_hex), width=width, lineCap=cap, fill=None, closePath=False); sh.commit()


def dashed_segment_v(page, x, y_top, y_bot, start_y, direction, color_hex, width, dashes):
    """Continue the sheet's dashed vertical reference line inside [y_top, y_bot]. start_y is the sheet line's own
    start point and direction +1 (drawn downward) or -1 (drawn upward), read from the base drawing item."""
    period = dashes[0] + dashes[1]
    if direction > 0:      # drawn top to bottom from start_y
        phase = (y_top - start_y) % period
        vline(page, x, y_top, y_bot, color_hex, width, cap=0, dashes=dashes, phase=phase)
    else:                  # drawn bottom to top from start_y
        phase = (start_y - y_bot) % period
        vline(page, x, y_bot, y_top, color_hex, width, cap=0, dashes=dashes, phase=phase)


def dashed_segment_h(page, y, x_left, x_right, start_x, direction, color_hex, width, dashes):
    period = dashes[0] + dashes[1]
    if direction > 0:
        phase = (x_left - start_x) % period
        hline(page, x_left, x_right, y, color_hex, width, cap=0, dashes=dashes, phase=phase)
    else:
        phase = (start_x - x_right) % period
        hline(page, x_right, x_left, y, color_hex, width, cap=0, dashes=dashes, phase=phase)


def ci_line(page, x0, x1, y, color_hex, width=1.7):
    """A confidence-interval line with round caps, the sheet's own mark."""
    sh = page.new_shape(); sh.draw_line(fitz.Point(x0, y), fitz.Point(x1, y)); sh.finish(color=hexrgb(color_hex), width=width, lineCap=1, fill=None, closePath=False); sh.commit()


def marker(page, cx, cy, kind, area_pt2, fill_hex, edge_hex, edge_w=0.8):
    """A matplotlib scatter marker: area s in pt^2 (circle diameter = sqrt(s), square side = sqrt(s)), edge centred."""
    sh = page.new_shape()
    if kind == "o":
        sh.draw_circle(fitz.Point(cx, cy), math.sqrt(area_pt2) / 2)
    elif kind == "s":
        h = math.sqrt(area_pt2) / 2
        sh.draw_rect(fitz.Rect(cx - h, cy - h, cx + h, cy + h))
    elif kind == "^":
        # matplotlib's triangle_up unit path (0,1),(-1,-1),(1,-1) scaled so its area equals s
        side = math.sqrt(area_pt2) * math.sqrt(2)         # for a 2 x 2 unit box the marker size is sqrt(s)
        hh = math.sqrt(area_pt2) / 2
        sh.draw_polyline([fitz.Point(cx, cy - hh), fitz.Point(cx - hh, cy + hh), fitz.Point(cx + hh, cy + hh), fitz.Point(cx, cy - hh)])
    else:
        raise ValueError(kind)
    sh.finish(fill=hexrgb(fill_hex) if fill_hex else None, color=hexrgb(edge_hex) if edge_hex else None, width=edge_w, closePath=True)
    sh.commit()


def redact_spans(page, spans, pad=0.0, tight=False):
    """Text-only redaction of the given spans' boxes (no fill drawn, graphics and images untouched). With tight=True the
    box is the glyph body only (baseline - 0.70 size to baseline + 0.12 size), so a neighbour 0.6 pt away is not hit."""
    for s in spans:
        if tight:
            r = fitz.Rect(s["bbox"][0], s["origin"][1] - 0.70 * s["size"], s["bbox"][2], s["origin"][1] + 0.12 * s["size"])
        else:
            r = fitz.Rect(*s["bbox"]); r = fitz.Rect(r.x0 - pad, r.y0 - pad, r.x1 + pad, r.y1 + pad)
        page.add_redact_annot(r, fill=False)
    if spans:
        assert page.apply_redactions(images=0, graphics=0, text=0), "apply_redactions returned False"


def drop_font_resource(page, basefont_prefix):
    """Remove an unused font resource (e.g. a base-14 Helvetica left behind by a redacted stamp) from the page's font
    dictionary. Returns the resource keys removed. Only call when no text on the page uses that face."""
    doc = page.parent; removed = []
    for f in page.get_fonts(full=True):
        if f[3].startswith(basefont_prefix):
            try:
                doc.xref_set_key(page.xref, f"Resources/Font/{f[4]}", "null"); removed.append(f[4])
            except Exception as e:
                raise RuntimeError(f"could not drop font resource {f[4]}: {e}")
    return removed


def fix_tounicode(doc):
    """The TextWriter-embedded Arial faces carry a ToUnicode CMap that maps the space glyph (gid 3) to U+00A0 and the
    hyphen glyph (gid 16) to U+00AD, so extracted text reads NBSP and soft hyphen. Rendering is unaffected. Repoint the
    two bfchar lines to U+0020 and U+002D on every 'Arial Regular' / 'Arial Bold' font of the document (the other
    lines, e.g. gid 0xA0 -> ae, are left alone). Returns the number of CMap streams patched."""
    done = set(); n = 0
    for pno in range(len(doc)):
        for f in doc[pno].get_fonts(full=True):
            xref, base = f[0], f[3]
            if base not in ("Arial Regular", "Arial Bold") or xref in done:
                continue
            done.add(xref)
            tu = doc.xref_get_key(xref, "ToUnicode")
            if tu[0] != "xref":
                continue
            tx = int(tu[1].split()[0]); s = doc.xref_stream(tx).decode("latin-1")
            s2 = s.replace("<0003> <00a0>", "<0003> <0020>").replace("<0010> <00ad>", "<0010> <002d>")
            if s2 != s:
                doc.update_stream(tx, s2.encode("latin-1")); n += 1
    return n


def text_width(text, size, bold=False):
    return font(bold).text_length(text, fontsize=size)


def stamp(page, text, x, baseline, size, color_hex, bold=False, align="left"):
    """TextWriter stamp with the sheet's Arial face. align: left (x = start), right (x = end), center (x = centre)."""
    w = text_width(text, size, bold)
    x0 = x if align == "left" else (x - w if align == "right" else x - w / 2)
    tw = fitz.TextWriter(page.rect); tw.append(fitz.Point(x0, baseline), text, font=font(bold), fontsize=size)
    tw.write_text(page, color=hexrgb(color_hex))
    return x0, w
