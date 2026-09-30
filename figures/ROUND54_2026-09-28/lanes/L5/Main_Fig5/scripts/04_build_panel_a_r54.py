"""
Round 54, lane L5 (2026-09-29): the Figure 5a flow generator (FINAL_FIGURES_2026-08-14/_scripts/figure5A_flow_polish.py, the panel the V13
sheet carries at scale S_V13 = 1.13435) rerun with the two text sizes given on the command line IN CANVAS POINTS (the round-54 sheet sizes
14 and 15 pt divided by S_V13: block and arm strings 12.342, captions 13.223), the arm text column X_ARM optional as a third number:
usage 04_build_panel_a_r54.py <tick_pt_canvas> <title_pt_canvas> <tag> [x_arm]. With 10 11 bump0 it reproduces the Aug-25 panel (the
fill remap apart). Writes only into this lane (work/panel_a/). Geometry in data units unchanged, so the flows, blocks and brackets (the
data marks) do not move. The round-49 and Aug-25 docstrings follow.
Round 49, lane L5 (2026-09-26): the Figure 5a flow generator rerun with every text size + BUMP where BUMP = 1/S_V13 pt on this canvas.
Figure5A, the pre/post FLOW render (round 3, 2026-08-25): the same 1,894 patients as
the on-PAP band bar chart, drawn as a paired pre-to-post flow. ROUND3_BRIEF_2026-08-25
MAINS (the item the brief labels "Fig 4a" after this sheet's round-1 workshop stem
Figure4A_prepost_VARIANT) and NUMBERING_MAP_V3_2026-08-25.md section 1f: final Fig. 5a
is a dual-render slot, both renders delivered, DEFAULT = this flow render if it gates
clean and reads clearly, with the band chart parked as the alternate. This builder is
the round-1 variant's lineage rebuilt at panel quality under the current numbering
(wider caption clearances, the final file name); the data asserts are carried verbatim.

WHAT IT DRAWS
  Left, one block: the pretreatment state. All 1,894 patients sit in the >10% band of
  sleep-period T90 on the diagnostic study, so the pre column is a single block in the
  band's own deep blue, its count printed on it.
  Right, four blocks: the on-PAP state, one block per T90 band on the strict role ramp
  (#ccd1d6, #8a9099, #4d7387, #1f4257, pale to deep), separated by white gaps, each
  with its band range and count printed on it.
  Between them, four flow bands, each in its destination band's colour, whose vertical
  thickness is the patient count it carries. The >10% stream runs deep blue from edge
  to edge: those 568 patients never left the band.
  Right of the stack, the two arms of Figure4A: restored to 10% or less (1,326, 70%)
  and still above 10% (568, 30%), bracketed exactly as on the band chart.

Every printed number is asserted against numbers/treatment_v2.json (migration,
from_worst and the migration matrix row of the >10% band) before anything is drawn.
The block and flow colours are the four steps of the declared house ramp, so every
adjacent pair is ramp-exempt in the colour gate and the order survives greyscale.

Design system as the other Figure 4 sheets: 171 mm canvas, Arial 11/10/9.5 with the
9 pt floor, ink #1a1d21, no titles, no bold, nothing within 4 mm of an edge. This
sheet is a flow diagram, so it carries no axes: every quantity a reader could take
off an axis is printed on the sheet instead.

splitstyle.save() is never called and legend_capture is disarmed exactly as in the
other polish scripts, so nothing writes into the frozen New_Figures tree. Output goes
only to Main_Figures_Polished (PDF) plus the 300 dpi PNG in _workfiles/polished_pngs
and a drawn-values JSON in _workfiles.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import atexit
import json
import os
import sys

import numpy as np

sys.path.insert(0, f"{paths.FIGURE_ROOT}/New_Figures/"
                   "_NOT_THESE_workfiles/scripts")
# round 49: splitstyle.py no longer imports (its negative-control assertion is pinned to a pre-v8 file); rc() and BAND_LABELS are
# taken from it verbatim (splitstyle.py lines 36-43, 57, 65, 97-114) over the same prism_plotter manuscript style
sys.path.insert(0, paths.PRISM_PLOTTER_SRC)
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from prism_plotter.style import manuscript_style, register_manuscript_font  # noqa: E402
register_manuscript_font()
BAND_LABELS = ["0-1%", "1-5%", "5-10%", ">10%"]
_SS_TITLE, _SS_LABEL, _SS_TICK, _SS_ANNOT, _SS_INK = 11.5, 10.0, 9.5, 9.5, "#15181c"


def rc():
    r = dict(manuscript_style().rc())
    r.update({
        "figure.dpi": 300, "savefig.dpi": 300,
        "font.size": _SS_TICK, "axes.labelsize": _SS_LABEL, "axes.titlesize": _SS_TITLE,
        "xtick.labelsize": _SS_TICK, "ytick.labelsize": _SS_TICK, "legend.fontsize": _SS_ANNOT,
        "axes.linewidth": 0.9, "axes.edgecolor": _SS_INK, "text.color": _SS_INK,
        "axes.labelcolor": _SS_INK, "xtick.color": _SS_INK, "ytick.color": _SS_INK,
        "xtick.major.width": 0.9, "ytick.major.width": 0.9,
        "xtick.major.size": 3.2, "ytick.major.size": 3.2,
        "xtick.direction": "out", "ytick.direction": "out",
        "axes.spines.top": False, "axes.spines.right": False,
        "legend.frameon": False, "savefig.bbox": "tight", "savefig.pad_inches": 0.06,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })
    return r


import nooverlap  # noqa: E402
import matplotlib.text as mtext  # noqa: E402
from matplotlib.path import Path as MplPath  # noqa: E402
from matplotlib.patches import PathPatch, Rectangle  # noqa: E402

ROOT = paths.FIGURE_ROOT
NUM = f"{paths.NUMBERS_DIR}"
TICK_ARG = float(sys.argv[1]) if len(sys.argv) > 1 else 10.0
TITLE_ARG = float(sys.argv[2]) if len(sys.argv) > 2 else 11.0
TAG = sys.argv[3] if len(sys.argv) > 3 else "bump0"
X_ARM_ARG = float(sys.argv[4]) if len(sys.argv) > 4 else 77.5
OUT = f"{paths.FIGURE_ROOT}/ROUND54_2026-09-28/lanes/L5/Main_Fig5/work/panel_a"
POLISH = OUT
WORK = OUT
PNGS = OUT
for _d in (POLISH, WORK, PNGS):
    os.makedirs(_d, exist_ok=True)
NAME = f"Figure5A_{TAG}"

# ---- the shared Figure 4 design system
MM = 1 / 25.4
FIGW = 171 * MM
INK_P = "#1a1d21"
PRIMARY, SECOND, COMPARE, DEEMPH = "#1f4257", "#4d7387", "#8a9099", "#ccd1d6"
# round 49: the sheet carries panel a in the clinical palette of the V13 lineage (measured on the V26 raster: the pre block and the >10%
# flow #0288d1, the 5-10% flow #3f9fd8, the 1-5% flow #7cc0e9, the 0-1% flow #b3dcf2, exact 8-bit matches of the v14lib PALETTE), the
# ink and the white text unchanged: the four fills are remapped one for one, everything else as the Aug-25 generator
PRIMARY, SECOND, COMPARE, DEEMPH = "#0288d1", "#3f9fd8", "#7cc0e9", "#b3dcf2"
RAMP = [DEEMPH, COMPARE, SECOND, PRIMARY]      # the ordered oxygen bands, pale to deep
TITLE_PT, TICK_PT, ANN_PT, FLOOR_PT = TITLE_ARG, TICK_ARG, 9.5, 9.0     # round 54: the canvas sizes from the command line (sheet size / S_V13)
# white text on the two deep steps, ink on the two pale steps, per band
BAND_TXT_COL = {DEEMPH: INK_P, COMPARE: INK_P, SECOND: "white", PRIMARY: "white"}


def margin_gate(fig, name, min_mm=4.0):
    """No ink within 4 mm of any sheet edge, asserted on the rendered extent."""
    fig.canvas.draw()
    bb = fig.get_tightbbox(fig.canvas.get_renderer())
    W, H = fig.get_size_inches()
    m = min_mm * MM
    bad = []
    if bb.x0 < m:
        bad.append(f"left {bb.x0 / MM:.1f} mm")
    if bb.y0 < m:
        bad.append(f"bottom {bb.y0 / MM:.1f} mm")
    if W - bb.x1 < m:
        bad.append(f"right {(W - bb.x1) / MM:.1f} mm")
    if H - bb.y1 < m:
        bad.append(f"top {(H - bb.y1) / MM:.1f} mm")
    assert not bad, f"{name}: content closer than {min_mm} mm to an edge: {', '.join(bad)}"


# ==========================================================================================
# data, asserted against the frozen source. Same pins as figure4A_polish.py plus the
# migration-matrix cross-check for the >10% row this sheet draws.
# ==========================================================================================
T = json.load(open(f"{NUM}/treatment_v2.json"))
mig, cvn = T["migration"], T["corrected_vs_not"]
N_ALL, N_OK, N_NO = cvn["n"], cvn["n_corrected"], cvn["n_not"]
assert (N_ALL, N_OK, N_NO) == (1894, 1326, 568), (N_ALL, N_OK, N_NO)
assert mig["labels"] == BAND_LABELS
BAND_TICKTXT = [b.replace("-", "–") for b in BAND_LABELS]
assert BAND_TICKTXT == ["0–1%", "1–5%", "5–10%", ">10%"], BAND_TICKTXT

FW = mig["from_worst"]
assert FW["n"] == N_ALL
DEST = [FW["to"][b] for b in BAND_LABELS]
assert DEST == [613, 498, 215, 568], DEST
assert sum(DEST[:3]) == N_OK and DEST[3] == N_NO
PCT_OK = round(100 * N_OK / N_ALL)
PCT_NO = round(100 * N_NO / N_ALL)
assert PCT_OK == 70 and PCT_NO == 30
DEST_PCT = [round(100 * v / N_ALL) for v in DEST]
assert DEST_PCT == [32, 26, 11, 30], DEST_PCT

# the migration matrix agrees with from_worst: its >10% row is this sheet's flow, its
# row sums are the pre-band totals and its column sums the post-band totals
MAT = mig["matrix"]
assert MAT[3] == DEST, (MAT[3], DEST)
assert [sum(row) for row in MAT] == mig["pre_totals"]
assert [sum(col) for col in zip(*MAT)] == mig["post_totals"]
assert mig["pre_totals"][3] == N_ALL

# ==========================================================================================
# geometry, in data units: x on 0..100, y in patients. The destination stack carries a
# white gap between blocks so the four flows separate; every block's HEIGHT is its
# count, and every count is printed, so no y axis is needed or drawn.
# ==========================================================================================
GAP = 70.0                                    # patients of white between the post blocks
X_SRC = (8.0, 22.0)                           # the pretreatment block
X_DST = (48.0, 72.0)                          # the on-PAP blocks
X_BRK = 74.6                                  # the arm bracket rule
X_BRK_TICK = 73.0                             # bracket end ticks point at the blocks
X_ARM = X_ARM_ARG                             # arm text, left-aligned (Aug-25: 77.5; round 54 may move the text left to clear panel b, text only)

# post blocks bottom to top: >10%, 5-10%, 1-5%, 0-1% (restored bands on top)
ORDER = [3, 2, 1, 0]                          # indices into BAND_LABELS / DEST / RAMP
d_lo, cur = {}, 0.0
for i in ORDER:
    d_lo[i] = cur
    cur += DEST[i] + GAP
TOP = cur - GAP                               # 1894 + 3 gaps
assert abs(TOP - (N_ALL + 3 * GAP)) < 1e-9

# the source block is centred against the stack; its slices follow the same order so
# the flows never cross
S_LO = 1.5 * GAP
s_lo, cur = {}, S_LO
for i in ORDER:
    s_lo[i] = cur
    cur += DEST[i]
assert abs(cur - (S_LO + N_ALL)) < 1e-9

Y_CAP = -170.0                                # column captions
Y_QTY = -400.0                                # the shared quantity line
FIGH = 4.30


def flow(ax, x0, y0a, y1a, x1, y0b, y1b, colour):
    """One flow band, cubic with horizontal tangents, no stroke."""
    xm0, xm1 = x0 + 0.42 * (x1 - x0), x0 + 0.58 * (x1 - x0)
    verts = [(x0, y0a),
             (xm0, y0a), (xm1, y0b), (x1, y0b),
             (x1, y1b),
             (xm1, y1b), (xm0, y1a), (x0, y1a),
             (x0, y0a)]
    codes = [MplPath.MOVETO,
             MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
             MplPath.LINETO,
             MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
             MplPath.CLOSEPOLY]
    ax.add_patch(PathPatch(MplPath(verts, codes), facecolor=colour, edgecolor="none",
                           lw=0, zorder=1))


# ==========================================================================================
# the sheet
# ==========================================================================================
RC = dict(rc())
RC.update({"font.size": TICK_PT, "axes.labelsize": TITLE_PT, "axes.titlesize": TITLE_PT,
           "xtick.labelsize": TICK_PT, "ytick.labelsize": TICK_PT,
           "legend.fontsize": TICK_PT,
           "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
           "xtick.major.size": 3.0, "ytick.major.size": 3.0,
           "axes.edgecolor": INK_P, "text.color": INK_P, "axes.labelcolor": INK_P,
           "xtick.color": INK_P, "ytick.color": INK_P,
           "savefig.bbox": "standard"})

DRAWN = {"pre": {}, "post": [], "arms": []}
with plt.rc_context(RC):
    fig = plt.figure(figsize=(FIGW, FIGH))
    ax = fig.add_axes([0.045, 0.045, 0.927, 0.912])
    ax.set_xlim(0, 100)
    ax.set_ylim(Y_QTY - 130, TOP + 30)
    ax.set_axis_off()

    # the pretreatment block: one band, the band's own deep blue
    ax.add_patch(Rectangle((X_SRC[0], S_LO), X_SRC[1] - X_SRC[0], N_ALL,
                           facecolor=PRIMARY, edgecolor="none", lw=0, zorder=2))
    src_txt = f"{BAND_TICKTXT[3]}\nn = {N_ALL:,}"
    ax.text(0.5 * (X_SRC[0] + X_SRC[1]), S_LO + 0.5 * N_ALL, src_txt, ha="center",
            va="center", fontsize=TICK_PT, color="white", linespacing=1.4)
    DRAWN["pre"] = {"band": BAND_LABELS[3], "n": N_ALL, "printed": src_txt}

    # the four flows and the four on-PAP blocks, each in its destination colour. The
    # flow ends run 1.2 units under the blocks (which draw above them in the same
    # colour), so no antialiasing seam shows where the two shapes meet.
    for i in ORDER:
        flow(ax, X_SRC[1] - 1.2, s_lo[i], s_lo[i] + DEST[i],
             X_DST[0] + 1.2, d_lo[i], d_lo[i] + DEST[i], RAMP[i])
        ax.add_patch(Rectangle((X_DST[0], d_lo[i]), X_DST[1] - X_DST[0], DEST[i],
                               facecolor=RAMP[i], edgecolor="none", lw=0, zorder=2))
        blk_txt = f"{BAND_TICKTXT[i]}, n = {DEST[i]}"
        ax.text(0.5 * (X_DST[0] + X_DST[1]), d_lo[i] + 0.5 * DEST[i], blk_txt,
                ha="center", va="center", fontsize=TICK_PT, color=BAND_TXT_COL[RAMP[i]])
        DRAWN["post"].append({"band": BAND_LABELS[i], "n": DEST[i],
                              "pct_of_1894": DEST_PCT[i], "printed": blk_txt})

    # the two arms of the band chart, bracketed beside the stack
    for y0, y1, lines, n in (
            (d_lo[2], TOP, ["Restored to", "10% or less", f"n = {N_OK:,} ({PCT_OK}%)"],
             N_OK),
            (0.0, DEST[3], ["Still above 10%", f"n = {N_NO} ({PCT_NO}%)"], N_NO)):
        ax.plot([X_BRK, X_BRK], [y0, y1], color=INK_P, lw=0.8, solid_capstyle="butt",
                zorder=3)
        ax.plot([X_BRK, X_BRK_TICK], [y0, y0], color=INK_P, lw=0.8, zorder=3)
        ax.plot([X_BRK, X_BRK_TICK], [y1, y1], color=INK_P, lw=0.8, zorder=3)
        arm_txt = "\n".join(lines)
        ax.text(X_ARM, 0.5 * (y0 + y1), arm_txt, ha="left", va="center",
                fontsize=TICK_PT, color=INK_P, linespacing=1.4)
        DRAWN["arms"].append({"n": n, "printed": arm_txt.replace("\n", " ")})

    # column captions and the shared quantity, the only 11 pt text on the sheet
    ax.text(0.5 * (X_SRC[0] + X_SRC[1]), Y_CAP, "Pretreatment\nsleep study",
            ha="center", va="center", fontsize=TITLE_PT, color=INK_P, linespacing=1.35)
    ax.text(0.5 * (X_DST[0] + X_DST[1]), Y_CAP, "On PAP", ha="center", va="center",
            fontsize=TITLE_PT, color=INK_P, linespacing=1.35)
    ax.text(0.5 * (X_SRC[0] + X_DST[1]), Y_QTY,
            "Time with oxygen below 90%, % of sleep", ha="center", va="center",
            fontsize=TITLE_PT, color=INK_P)

    # ---- gates, the Figure 4 standard, then write
    for t in fig.findobj(mtext.Text):
        s = t.get_text()
        if not s or not s.strip() or not t.get_visible():
            continue
        for bad, nm in (("—", "em dash"), (";", "semicolon"),
                        ("×", "multiplication sign")):
            assert bad not in s, f"banned {nm} in on-sheet text: {s!r}"
        # en dash is legal only inside a numeric range (the standard band form 0-1%)
        for i, ch in enumerate(s):
            if ch == "–":
                assert 0 < i < len(s) - 1 and s[i - 1].isdigit() and s[i + 1].isdigit(), (
                    f"en dash outside a numeric range in on-sheet text: {s!r}")
        assert str(t.get_fontweight()) not in ("bold", "700"), (
            f"bold is not allowed on this sheet: {s!r}")
        assert float(t.get_fontsize()) >= FLOOR_PT, (
            f"{t.get_fontsize()} pt under the {FLOOR_PT} pt floor: {s!r}")

    nooverlap.gate(fig, NAME)
    margin_gate(fig, NAME)
    fig.savefig(f"{POLISH}/{NAME}.pdf")
    fig.savefig(f"{PNGS}/{NAME}.png", dpi=300)
    plt.close(fig)

DRAWN["flow_total"] = sum(d["n"] for d in DRAWN["post"])
assert DRAWN["flow_total"] == N_ALL == DRAWN["pre"]["n"]
assert sum(a["n"] for a in DRAWN["arms"]) == N_ALL
DRAWN["bump_on_canvas_pt"] = TICK_PT - 10.0; DRAWN["text_pt_on_canvas"] = {"blocks_arms": TICK_PT, "captions": TITLE_PT}; DRAWN["x_arm"] = X_ARM
json.dump(DRAWN, open(f"{WORK}/{NAME}_flow_drawn_values.json", "w"), indent=1,
          default=float)

print("written:")
for p in (f"{POLISH}/{NAME}.pdf", f"{PNGS}/{NAME}.png",
          f"{WORK}/{NAME}_flow_drawn_values.json"):
    assert os.path.exists(p) and os.path.getsize(p) > 0, p
    print("  ", p, f"{os.path.getsize(p):,} bytes")
