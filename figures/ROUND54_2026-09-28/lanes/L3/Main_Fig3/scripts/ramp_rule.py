"""The hazard-ratio red ramp of Figs 3c, 4d and ED 7b (round 27 lane B), one implementation shared with lane V7_L4_APNEA
(copied from its l4lib.py on 2026-09-09 so the twin grids follow one rule). Old ramp (R15 build_r15_fig3b.py): RGB stops at
log-spaced hazard ratios 1, 1.25, 1.90, 3.25. Round 27 mapped every old fill by its CIELAB lightness position
t = (L_palest - L) / (L_palest - L_darkest), palest #fdefe9, darkest #d15641 (the top of the old gradient bar), onto the
straight CIELAB line from #fdecea (t = 0) to #d32f2f (t = 1). Positive control on the round-30 base Fig 3c: 35 of 35 cells
reproduced within one level per channel from their printed hazard ratios. Text: ink where it reaches 4.5:1 on the fill,
else white where white reaches it, else the higher of the two (round 27: 3.18*** and 3.08*** white, 2.76*** ink)."""
import math, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import colour_math as CM
RAMP_OLD_STOPS = [(1.00, "#fdefe9"), (1.25, "#fbe3dc"), (1.90, "#f2a98f"), (3.25, "#d15540")]
RAMP_MIN, RAMP_MAX = 1.0, 3.25
RAMP_PALE, RAMP_TOP = "#fdecea", "#d32f2f"
_L_PALE = CM.rgb_to_lab(CM.hex_to_rgb("#fdefe9"))[0]; _L_DARK = CM.rgb_to_lab(CM.hex_to_rgb("#d15641"))[0]


def ramp_old(hr):
    x = min(max(math.log(hr), math.log(RAMP_MIN)), math.log(RAMP_MAX))
    for (h0, c0), (h1, c1) in zip(RAMP_OLD_STOPS[:-1], RAMP_OLD_STOPS[1:]):
        if x <= math.log(h1) + 1e-12:
            t = (x - math.log(h0)) / (math.log(h1) - math.log(h0)); a, b = CM.hex_to_rgb(c0), CM.hex_to_rgb(c1)
            return tuple(u + t * (v - u) for u, v in zip(a, b))
    return CM.hex_to_rgb(RAMP_OLD_STOPS[-1][1])


def ramp_t(hr):
    L = CM.rgb_to_lab(ramp_old(hr))[0]; return min(1.0, max(0.0, (_L_PALE - L) / (_L_PALE - _L_DARK)))


def ramp_hex(hr):
    t = ramp_t(hr); a, b = CM.rgb_to_lab(CM.hex_to_rgb(RAMP_PALE)), CM.rgb_to_lab(CM.hex_to_rgb(RAMP_TOP))
    return CM.rgb_to_hex(CM.lab_to_rgb(tuple(u + t * (v - u) for u, v in zip(a, b))))


def text_on(fill_hex, ink="#1a1d21"):
    c_ink = CM.contrast(CM.hex_to_rgb(ink), CM.hex_to_rgb(fill_hex)); c_w = CM.contrast((1, 1, 1), CM.hex_to_rgb(fill_hex))
    return ink if c_ink >= 4.5 else ("#ffffff" if c_w >= 4.5 else (ink if c_ink >= c_w else "#ffffff"))


def contrasts(fill_hex, ink="#1a1d21"):
    return CM.contrast(CM.hex_to_rgb(ink), CM.hex_to_rgb(fill_hex)), CM.contrast((1, 1, 1), CM.hex_to_rgb(fill_hex))
