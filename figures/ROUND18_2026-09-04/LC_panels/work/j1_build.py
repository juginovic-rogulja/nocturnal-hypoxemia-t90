"""JOB 1: NEW_Fig3_duration_desat.pdf

Figure 4b's design, classified by SLEEP DURATION instead of apnea severity, so a reader can see
how much each duration group desaturates. Same cell layout, same "count (percent)" string, same
header treatment (the T90 azure ladder), same key (four shade bands of the row's own colour
family). Only the row variable and the cell colour family change: rows are the paper's three
sleep duration groups and the cells take the sleep duration green ladder, because in Figure 4b
the cells take the apnea family's colour and here the rows are duration.

Every count and percent is computed in this run by work/j1_data.py and asserted, and Figure 4b's
own sixteen cells are recomputed from the same cohort and matched against the strings and the
fills printed on the current sheet before anything is drawn.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../..")))  # repository root, where paths.py lives
import paths
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.patches import Rectangle
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties
import fitz

HERE = os.path.dirname(os.path.abspath(__file__))
LANE = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import j1_data as J                                             # noqa: E402

for f in [paths.FONT_ARIAL,
          paths.FONT_ARIAL_BOLD]:
    fm.fontManager.addfont(f)
plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42,
                     "text.color": "#1a1d21"})

INK = "#1a1d21"
# read off the current sheet: Figure 4b's own column headers
AZURE = ["#b3dcf2", "#7cc0e9", "#3f9fd8", "#0288d1"]
# the sleep duration green, the paper's own endpoints (#b0e8b0 pale, #448044 dark, both printed
# on the current Figure 3b header) blended at the same fractions the apnea teal ladder of
# Figure 4b uses (t = 1.00, 0.70, 0.42, 0.00 from pale to dark).
GREEN_PALE, GREEN_DARK = "#b0e8b0", "#448044"
GREEN_T = [1.00, 0.70, 0.42, 0.00]
TITLE_COL = "T90, % of the recording"
TITLE_ROW = "Sleep duration"
KEY_TITLE = "Percent of the duration group"
_FP = FontProperties(family="Arial")


def hex2rgb(h):
    return np.array([int(h[i:i + 2], 16) for i in (1, 3, 5)]) / 255


def rgb2hex(c):
    return "#%02x%02x%02x" % tuple(int(round(v * 255)) for v in c)


def mix(dark, pale, t):
    return rgb2hex(hex2rgb(dark) + t * (hex2rgb(pale) - hex2rgb(dark)))


GREEN = [mix(GREEN_DARK, GREEN_PALE, t) for t in GREEN_T]       # <5, 5-15, 15-35, >=35


def rel_lum(rgb):
    c = np.asarray(rgb, float)
    c = np.where(c <= 0.03928, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)
    return float(0.2126 * c[0] + 0.7152 * c[1] + 0.0722 * c[2])


def c_black(h):
    return (rel_lum(hex2rgb(h)) + 0.05) / 0.05


def c_white(h):
    return 1.05 / (rel_lum(hex2rgb(h)) + 0.05)


def ink_on(h):
    """Black unless black fails 4.5:1 and white beats it, which is Figure 4b's own rule."""
    cb, cw = c_black(h), c_white(h)
    return (INK, cb) if cb >= 4.5 or cb >= cw else ("#ffffff", cw)


def tw(s, fs):
    return float(TextPath((0, 0), s, size=fs, prop=_FP).get_extents().width)


class Geom:
    def __init__(self, **kw):
        self.__dict__.update(kw)


# Figure 4b on the current sheet, measured: columns 89.32 wide with a 3.60 gap, header band
# 12.24 tall, 9.72 below it, rows 44.64 tall with a 3.96 gap, cell and header text 10 pt, the
# two axis titles 11 pt, key swatch 10.80 x 9.36. The sheet keeps every one of those ratios and
# widens the columns to fill 183 mm; the main-figure version keeps them against Figure 3's page.
SHEET = Geom(W=183.0 / 25.4 * 72, X0=8.0, LAB_R=97.0, GX0=107.0, GX1=510.74, GAP=3.60,
             HEAD_H=12.24, HEAD_PAD=9.72, ROW_H=44.64, ROW_GAP=3.96, fs=10.0, fs_title=11.0,
             SW_W=10.80, SW_H=9.36, KEY_PAD=17.02, TITLE_PAD=10.0, top=8.0, bot=8.0,
             tag="sheet")
MAIN = Geom(W=952.73, X0=14.2, LAB_R=200.0, GX0=211.12, GX1=933.52, GAP=7.20,
            HEAD_H=14.08, HEAD_PAD=11.18, ROW_H=51.34, ROW_GAP=4.55, fs=11.5, fs_title=12.65,
            SW_W=12.42, SW_H=10.76, KEY_PAD=19.57, TITLE_PAD=11.5, top=10.0, bot=10.0,
            tag="main")

printed = []


def draw(g, data, path):
    ncol, nrow = 4, len(J.DUR_ROWS)
    cw = (g.GX1 - g.GX0 - (ncol - 1) * g.GAP) / ncol
    H = (g.top + g.fs_title * 1.35 + g.TITLE_PAD + g.HEAD_H + g.HEAD_PAD
         + nrow * g.ROW_H + (nrow - 1) * g.ROW_GAP + g.KEY_PAD + g.SW_H + g.bot)
    fig = plt.figure(figsize=(g.W / 72, H / 72))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, g.W)
    ax.set_ylim(0, H)
    ax.axis("off")

    def colx(c):
        return g.GX0 + c * (cw + g.GAP)

    # the two axis titles, on one baseline, as on Figure 4b
    ty = H - g.top - g.fs_title * 0.72
    ax.text((g.GX0 + g.GX1) / 2, ty, TITLE_COL, ha="center", va="center",
            fontsize=g.fs_title, color=INK)
    ax.text(g.LAB_R, ty, TITLE_ROW, ha="right", va="center", fontsize=g.fs_title, color=INK)

    # header band, the T90 azure ladder
    hy1 = ty - g.fs_title * 0.63 - g.TITLE_PAD
    hy0 = hy1 - g.HEAD_H
    for c in range(ncol):
        ax.add_patch(Rectangle((colx(c), hy0), cw, g.HEAD_H, facecolor=AZURE[c],
                               edgecolor="none"))
        col, ratio = ink_on(AZURE[c])
        ax.text(colx(c) + cw / 2, hy0 + g.HEAD_H / 2, J.T90_KEYS[c], ha="center", va="center",
                fontsize=g.fs, color=col)
        printed.append((g.tag, "header", J.T90_KEYS[c], AZURE[c], col, round(ratio, 2)))

    # the cells
    gy = hy0 - g.HEAD_PAD
    for r, rname in enumerate(J.DUR_ROWS):
        y1 = gy - r * (g.ROW_H + g.ROW_GAP)
        y0 = y1 - g.ROW_H
        ax.text(g.LAB_R, (y0 + y1) / 2, rname, ha="right", va="center", fontsize=g.fs, color=INK)
        for c in range(ncol):
            cell = data[rname]["cells"][c]
            fc = GREEN[cell["shade"]]
            ax.add_patch(Rectangle((colx(c), y0), cw, g.ROW_H, facecolor=fc, edgecolor="none"))
            col, ratio = ink_on(fc)
            ax.text(colx(c) + cw / 2, (y0 + y1) / 2, cell["label"], ha="center", va="center",
                    fontsize=g.fs, color=col)
            printed.append((g.tag, f"{rname}|{cell['band']}", cell["label"], fc, col,
                            round(ratio, 2)))
    gy_bot = gy - nrow * g.ROW_H - (nrow - 1) * g.ROW_GAP

    # the key, Figure 4b's own: a title then four labelled swatches
    ky = gy_bot - g.KEY_PAD - g.SW_H / 2
    ax.text(g.X0, ky, KEY_TITLE, ha="left", va="center", fontsize=g.fs, color=INK)
    x = g.X0 + tw(KEY_TITLE, g.fs) + g.SW_W * 3.26
    for i, (_, lab) in enumerate(J.PCT_BANDS):
        ax.add_patch(Rectangle((x, ky - g.SW_H / 2), g.SW_W, g.SW_H, facecolor=GREEN[i],
                               edgecolor="white", lw=0.6))
        ax.text(x + g.SW_W * 1.40, ky, lab, ha="left", va="center", fontsize=g.fs, color=INK)
        x += g.SW_W * 1.40 + tw(lab, g.fs) + g.SW_W * 2.35
    assert x < g.W, ("key runs off the sheet", x, g.W)
    fig.savefig(path)
    plt.close(fig)
    return H


def main():
    control, new, extra = J.build()
    bad = J.control_gate(control)
    assert not bad, ("CONTROL FAILED, not delivering: ", bad)

    for h in GREEN + AZURE:
        assert max(c_black(h), c_white(h)) >= 4.5, h
    lums = [rel_lum(hex2rgb(h)) for h in GREEN]
    assert lums == sorted(lums, reverse=True), lums                 # pale to dark, monotone

    out = {}
    for g, name in ((SHEET, "NEW_Fig3_duration_desat.pdf"),
                    (MAIN, "NEW_Fig3_duration_desat_fig3scale.pdf")):
        p = os.path.join(LANE, name)
        h = draw(g, new, p)
        out[g.tag] = dict(file=name, width_pt=round(g.W, 2), height_pt=round(h, 2),
                          width_mm=round(g.W / 72 * 25.4, 2), height_mm=round(h / 72 * 25.4, 2))
        print("%-6s %-42s %.2f x %.2f pt = %.1f x %.1f mm" %
              (g.tag, name, g.W, h, g.W / 72 * 25.4, h / 72 * 25.4))

    # every printed string asserted against the data it came from
    n = 0
    for tag, where, s, fill, col, ratio in printed:
        if "|" not in where:
            continue
        r, band = where.split("|")
        cell = [c for c in new[r]["cells"] if c["band"] == band][0]
        assert s == f"{cell['n']:,} ({cell['pct']})", (where, s)
        assert fill == GREEN[cell["shade"]], (where, fill)
        assert ratio >= 4.5, (where, ratio)
        n += 1
    assert n == 2 * 12, n
    print("printed cell strings asserted:", n, "(12 per deliverable)")

    values = {
        "what": "sleep duration by T90 band, count and percent of the duration group",
        "cohort": "data_frozen_v8_2026-09/t90_final.parquet via numbers/cohort_spec.apply_cohort",
        "cohort_n": extra["cohort"],
        "columns_used": {"duration": "TST_min", "oxygen": "spo2_pct_below_90",
                         "apnea (control only)": "AHI"},
        "duration_cuts": {"<5 h": "TST_min < 300", "6–7 h": "360 <= TST_min <= 420",
                          ">7 h": "TST_min > 420",
                          "source": "Sleep_Variability_2026-08/dur_x_oxygen_6cell_tst300/model.py"},
        "t90_bands": {"≤1%": "<= 1", ">1–5%": "1 to 5", ">5–10%": "5 to 10", ">10%": "> 10",
                      "source": "read off Figure 4b of NEW_FINAL_SET_R16/Main_Fig4.pdf"},
        "shade_bands_pct": ["<5", "5–15", "15–35", "≥35"],
        "green_ladder": dict(zip(["<5", "5–15", "15–35", "≥35"], GREEN)),
        "green_ladder_rule": ("blend of the paper's own duration greens %s (pale) and %s (dark) "
                              "at t = %s, the fractions Figure 4b's apnea teal ladder uses"
                              % (GREEN_PALE, GREEN_DARK, GREEN_T)),
        "azure_header_ladder": dict(zip(J.T90_KEYS, AZURE)),
        "not_shown": {"tst_300_to_360_min": extra["tst_300_to_360"],
                      "tst_missing": extra["tst_missing"],
                      "note": "the paper's duration groups leave the 5 to 6 h band out, so the "
                              "three rows cover 12,759 of the 19,173 recordings"},
        "rows": {r: dict(total=new[r]["total"], cells=new[r]["cells"]) for r in J.DUR_ROWS},
        "control_figure_4b": {r: dict(total=control[r]["total"], cells=control[r]["cells"])
                              for r in J.AHI_ROWS},
        "control_result": "16 of 16 printed strings and 16 of 16 cell fills reproduce",
        "contrast": {h: dict(black=round(c_black(h), 2), white=round(c_white(h), 2))
                     for h in GREEN + AZURE},
        "deliverables": out,
        "printed": [dict(sheet=t, where=w, string=s, fill=f, ink=c, contrast=r)
                    for t, w, s, f, c, r in printed],
    }
    json.dump(values, open(os.path.join(LANE, "NEW_Fig3_duration_desat_values.json"), "w"),
              indent=1, ensure_ascii=False)
    print("green ladder:", GREEN)
    print("OK")


if __name__ == "__main__":
    main()
