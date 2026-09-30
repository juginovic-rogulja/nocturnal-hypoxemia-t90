#!/usr/bin/env python3
"""Round 38, L5 SHEETS (2026-09-16): the V13 base sheets' text layers WITHOUT a PyMuPDF text extraction on a nested sheet.

Source: the coordinator's probe of the V13 set, ROUND37_2026-09-14/figures/_common/sheet_geometry_V13.json (written 2026-09-15 09:58
by probe_sheets_V13.py, one PyMuPDF get_text('dict') pass per V13 sheet, before the kernel panic of 2026-09-16 12:26). Per sheet it
holds every span as [text, origin_x, origin_y, font, size, colour_int, bbox] (text stripped, coordinates rounded to 2 decimals), the
letters, the font list, the page box. This module serves those spans in the lane formats:
  spans_lane(sheet)   -> the lane_common.spans_of dict form (text, font, size, color, bbox, origin)
  spans_v14(sheet)    -> the v14lib.spans 7-tuple form (text, x, y, font, size, colour_int, bbox)
  letters_lane(sheet) -> the lane_common.letters tuple form (text, bbox x0, bbox y0, bbox y1, font, size)
  fonts(sheet)        -> sorted base font names (subset prefixes removed)
  page(sheet)         -> (W, H)
The probe's sha256 is recorded by the caller (SOURCE_SCRIPTS_R38.md) so the text layer used is auditable."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import hashlib, json, os

GEOM = f"{paths.FIGURE_ROOT}/ROUND37_2026-09-14/figures/_common/sheet_geometry_V13.json"
_G = None


def _load():
    global _G
    if _G is None: _G = json.load(open(GEOM))
    return _G


def sha256():
    h = hashlib.sha256(); h.update(open(GEOM, "rb").read()); return h.hexdigest()


def raw(sheet): return _load()[f"{sheet}.pdf"]


def page(sheet): return tuple(raw(sheet)["page"])


def spans_lane(sheet):
    return [dict(text=t, font=f, size=round(float(sz), 3), color=int(c), bbox=[round(float(v), 3) for v in bb], origin=[round(float(x), 3), round(float(y), 3)])
            for t, x, y, f, sz, c, bb in raw(sheet)["spans"]]


def spans_v14(sheet):
    return [(t, float(x), float(y), f, float(sz), int(c), tuple(float(v) for v in bb)) for t, x, y, f, sz, c, bb in raw(sheet)["spans"]]


def letters_lane(sheet):
    return sorted([(s["text"].strip(), round(s["bbox"][0], 2), round(s["bbox"][1], 2), round(s["bbox"][3], 2), s["font"], round(s["size"], 1))
                   for s in spans_lane(sheet) if len(s["text"].strip()) == 1 and s["text"].strip() in "abcdefgh" and s["size"] >= 12], key=lambda t: (t[2], t[1]))


def fonts(sheet): return sorted(set(f.split("+")[-1] for f in raw(sheet)["fonts"]))


if __name__ == "__main__":
    import sys
    for sh in sys.argv[1:] or ["Main_Fig3", "Main_Fig4"]:
        print(sh, page(sh), len(spans_lane(sh)), "spans", letters_lane(sh), fonts(sh))
    print("probe sha256", sha256())
