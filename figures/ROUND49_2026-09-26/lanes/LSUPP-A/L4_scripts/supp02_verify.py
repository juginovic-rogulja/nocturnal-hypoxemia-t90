"""
Proof W for Supp_Fig02 (round 30, lane V7_L4_APNEA): values change, geometry pinned.

  1  page box 484.72 x 288.0 within 0.02, letters a and b at the base positions within 0.3 pt, Arial Bold 13
  2  geometry pins against the base sheet: both axes boxes (spines), the x tick marks, the x tick labels, the
     y-axis labels, the x labels and the key strings within 0.3 pt, fonts and sizes identical. The y tick set is
     data-driven (the builder's rule): its labels are pinned to the NEW calibration (right edge as the base, each
     label centred on its tick mark as the base centres its labels) and the change is declared.
  3  read-back: every marker centre, whisker end, median-line vertex and interquartile-polygon vertex of the
     NEW sheet, converted to data through the spines and tick marks, within 0.02 percentage points of the
     regenerated summary
  4  word multiset identical to the base after the round-28 dedupe, numeric-token multiset identical except
     for the declared y tick change (the expected delta is derived from the two tick sets, not typed)
  5  raster: Supp_Fig02_150dpi.png, OLD vs NEW crops of the key row, panel a and panel b at 200 dpi
  6  house words on every string of the sheet
  7  colours: the four band colours of the new sheet are the base sheet's, in the same key order
  8  counts: the drawn band sizes equal the summary's band_n_cohort and the v7 recompute, the stage panel n
     equals the summary's and n_ok minus the exclusion list's row count

Also writes verify/CHANGES_Supp_Fig02.csv (old values read from the base sheet's drawings, calibrated the same
way) and verify/readback_base.json, verify/readback_new.json.

    python3 verify_supp_fig02.py
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import collections
import csv
import datetime
import json
import os
import re
import sys

sys.dont_write_bytecode = True

T90 = paths.FIGURE_ROOT
R30 = f"{paths.FIGURE_ROOT}/ROUND30_2026-09-08"
LANE = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LSUPP-A"   # round 49: lane LSUPP-A
SHEET = f"{LANE}/Supp_Fig02"
sys.path.insert(0, f"{LANE}/L4_scripts")
import l4lib as L  # noqa: E402
import supp02_readback as RB  # noqa: E402
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND37_2026-09-14/figures/_common")
import v14lib as V  # noqa: E402
import fitz  # noqa: E402

NAME = "Supp_Fig02"
BASE = f"{L.BASE}/{NAME}.pdf"
NEW = f"{SHEET}/{NAME}.pdf"
DRAWN = f"{SHEET}/verify/{NAME}_drawn.json"
SUMMARY = f"{L.SV}/oxygen_profile/overnight_profile_summary.json"
EXCL = "oxyprofile_per_patient.csv:_v8_stage_collapsed"   # v8.2: the flag column
OLD_BUILDER = f"{paths.FIGURE_ROOT}/FINAL_FIGURES_2026-08-14/_scripts/efigNEW_overnight_oxygen_profile.py"
OLD_LEGEND_CHECK = f"{L.SV}/oxygen_profile/legend_check_output.txt"
CHECKS = f"{SHEET}/verify/checks.txt"
CHANGES = f"{SHEET}/verify/CHANGES_{NAME}.csv"
PNG150 = f"{SHEET}/{NAME}_150dpi.png"
CROPS = f"{SHEET}/verify/crops"
os.makedirs(CROPS, exist_ok=True)

BAND_LABEL = {"le1": "≤1%", "1to5": ">1 to 5%", "5to10": ">5 to 10%", "gt10": ">10%"}
BANDS = ["le1", "1to5", "5to10", "gt10"]
STAGES = ["Wake", "N1", "N2", "N3", "REM"]
BASE_YTICKS = None   # v8.2: read from the V13 sheet's own tick labels below (data-driven)
READ_TOL = 0.02                               # percentage points
PIN_TOL = 0.3                                 # points

PT_PLUS = 1.0                                 # round 49: every text one point larger than V13/V26
KEY_SHIFT_PT = -3.6                           # round 49: the shared key raised 3.6 pt as one block (declared, see the builder)
LAYOUT_TOL = 4.0                              # round 49: matplotlib re-lays out tick labels, axis titles and letters for the larger font (pts)
C = L.Checks(f"{NAME} proof W, lane LSUPP-A round 49 (+1 pt on the v8.2 build: V13 design, y axis pinned to the V13 tick set, key raised {-KEY_SHIFT_PT:g} pt), {datetime.datetime.now().strftime('%Y-%m-%d %H:%M')}")

# ------------------------------------------------------------------ inputs (gated)
sc = L.sidecar(SUMMARY)
C.check("summary sidecar: sha match, v8.1 tables named, step 153 on 2026-09-15 after 19:18 rc 0 (the v8.2 re-extraction)",
        str(sc["step"]["id"]) == "153" and sc["step"]["rc"] == 0 and str(sc["step"]["ended_local"]).startswith("2026-09-15") and sc["output"]["mtime_local"] >= "2026-09-15 19:18:00",
        f"{sc['step']['ended_local']} rc {sc['step']['rc']}")
S = json.load(open(L.hydrated(SUMMARY)))
D = json.load(open(L.hydrated(DRAWN)))
C.check("drawn JSON names the same summary (sha256)", D["sources"]["summary"]["sha256"] == sc["output"]["sha256"])
C.check("new PDF provenance sidecar present and current",
        os.path.exists(NEW + ".provenance.json") and json.load(open(NEW + ".provenance.json"))["output"]["sha256"] == L.sha256(NEW))

# ------------------------------------------------------------------ text layers
sb, sn = L.spans_of(BASE), L.spans_of(NEW)
wb, wn = L.dedupe(L.words_of(BASE)), L.dedupe(L.words_of(NEW))
base_ytick_texts = sorted({s["text"] for s in sb if s["font"] == "ArialMT" and s["size"] == 10.0 and s["bbox"][2] < 60}, key=int)
BASE_YTICKS = sorted(int(t) for t in set(base_ytick_texts))   # v8.2: the V13 sheet's own y tick set (86 to 98 by 4 on the delivered V13)
C.check("base sheet y tick labels read from its text layer form the V13 tick set (data-driven, not typed)", len(base_ytick_texts) >= 2 and [int(t) for t in sorted(set(base_ytick_texts), key=int)] == BASE_YTICKS, str(base_ytick_texts))
new_yticks = [float(t) for t in D["y_ticks"]]
new_ytick_texts = D["y_tick_labels"]

# ------------------------------------------------------------------ 1 page box and letters
doc = fitz.open(L.hydrated(NEW))
pw, ph = doc[0].rect.width, doc[0].rect.height
fonts_new = sorted((f[3].split("+")[-1], f[2]) for f in doc[0].get_fonts(full=True))
doc.close()
C.check("page box 484.72 x 288.0 within 0.02", abs(pw - 484.72) <= 0.02 and abs(ph - 288.0) <= 0.02, f"{pw:.3f} x {ph:.3f}")
C.check("fonts embedded: ArialMT and Arial-BoldMT, Type0", fonts_new == [("Arial-BoldMT", "Type0"), ("ArialMT", "Type0")], str(fonts_new))
lb, ln = L.letters(sb), L.letters(sn)
ok_letters = len(lb) == len(ln) == 2
for (tb, xb, yb, ob, fb, zb), (tn, xn, yn, on, fn, zn) in zip(lb, ln):
    ok_letters &= tb == tn and abs(xb - xn) <= PIN_TOL and abs(ob - on) <= PIN_TOL and fn == "Arial-BoldMT" and abs(zn - (13.0 + PT_PLUS)) < 0.01   # round 49: left edge and baseline (the bbox top rises with the larger ascent), 14 pt
C.check("letters a and b at the base positions within 0.3 pt (left edge, baseline), Arial-BoldMT 14 (V13 13 plus one point)", ok_letters, f"base {lb} new {ln}")

# ------------------------------------------------------------------ 3 (first) read-back of both sheets
rb = RB.read_sheet(BASE, BASE_YTICKS)
rn = RB.read_sheet(NEW, new_yticks)
json.dump(rb, open(f"{SHEET}/verify/readback_base.json", "w"), indent=1)
json.dump(rn, open(f"{SHEET}/verify/readback_new.json", "w"), indent=1)
C.check("axis calibration residuals under 0.001 data units on both sheets",
        max(max(rb["calib_a"]["fit_residual_data"]), max(rb["calib_b"]["fit_residual_data"]),
            max(rn["calib_a"]["fit_residual_data"]), max(rn["calib_b"]["fit_residual_data"])) < 1e-3,
        f"base a {rb['calib_a']['fit_residual_data']} b {rb['calib_b']['fit_residual_data']} new a {rn['calib_a']['fit_residual_data']} b {rn['calib_b']['fit_residual_data']}")
C.check("new sheet axis limits from the spines: y = drawn y_axis, x = drawn x limits (within 0.01)",
        all(abs(a - b) <= 0.01 for a, b in zip(rn["calib_a"]["ylim"], D["y_axis"])) and all(abs(a - b) <= 0.01 for a, b in zip(rn["calib_b"]["ylim"], D["y_axis"]))
        and all(abs(a - b) <= 0.01 for a, b in zip(rn["calib_a"]["xlim"], D["x_lim_a"])) and all(abs(a - b) <= 0.01 for a, b in zip(rn["calib_b"]["xlim"], D["x_lim_b"])),
        f"y a {rn['calib_a']['ylim']} b {rn['calib_b']['ylim']}, x a {rn['calib_a']['xlim']} b {rn['calib_b']['xlim']}")
C.check("no stray band-coloured marks outside the axes or the key on either sheet", not rb["stray"] and not rn["stray"], f"base {rb['stray']} new {rn['stray']}")

# ------------------------------------------------------------------ 2 geometry pins
C.check("panel a axes box (left and bottom spines) within 0.3 pt of the base", all(abs(a - b) <= PIN_TOL for a, b in zip(rb["axes_a"], rn["axes_a"])), f"base {rb['axes_a']} new {rn['axes_a']}")
C.check("panel b axes box (left and bottom spines) within 0.3 pt of the base", all(abs(a - b) <= PIN_TOL for a, b in zip(rb["axes_b"], rn["axes_b"])), f"base {rb['axes_b']} new {rn['axes_b']}")
C.check("x tick marks of both panels within 0.3 pt of the base",
        all(abs(a - b) <= PIN_TOL for a, b in zip(rb["calib_a"]["x_ticks_page"], rn["calib_a"]["x_ticks_page"]))
        and all(abs(a - b) <= PIN_TOL for a, b in zip(rb["calib_b"]["x_ticks_page"], rn["calib_b"]["x_ticks_page"])),
        f"a base {rb['calib_a']['x_ticks_page']} new {rn['calib_a']['x_ticks_page']}")


def is_ytick_label(s):
    sz = s["size"] in (10.0, 10.0 + PT_PLUS) and re.fullmatch(r"\d+", s["text"].strip()) is not None   # round 49: 10 pt on the base, 11 pt new, digits only (the rotated axis title is 11 pt at the same x)
    return s["font"] == "ArialMT" and sz and s["bbox"][2] < 60 or (s["font"] == "ArialMT" and sz and 300 < s["bbox"][0] < 330)


def nearest(span, pool):
    cands = [p for p in pool if p["text"] == span["text"]]
    if not cands:
        return None
    return min(cands, key=lambda p: abs(p["bbox"][0] - span["bbox"][0]) + abs(p["bbox"][1] - span["bbox"][1]))


# round 49: every string is one point larger, so matplotlib re-lays out the tick labels, the axis titles and the letters by fractions of a point
# to a few points, and the shared key (centred as a block) widens and was raised 3.6 pt. The pin therefore becomes: same text, same face,
# size + 1, within LAYOUT_TOL of the V13 position (the key strings: their declared shift in y, and the key block's centre within 0.5 pt).
pinned, worst, unmatched, font_mismatch, key_b, key_n, devs = 0, 0.0, [], [], [], [], []
for s in sb:
    if is_ytick_label(s) or not s["text"].strip():
        continue
    m = nearest(s, sn)
    if m is None:
        unmatched.append(s["text"])
        continue
    is_key = s["bbox"][1] < 50 and s["font"] == "ArialMT"
    # a one-point larger string keeps its ANCHOR, not its origin: the left origin (left-aligned), the centre (the axis titles, centred by
    # matplotlib) or the right edge; a rotated string (the y titles) keeps its x and its centre along y. The smallest anchor deviation counts.
    def anchors(t):
        b = t["bbox"]; rot = (b[3] - b[1]) > (b[2] - b[0]) and len(t["text"].strip()) > 2
        if rot: return [(t["origin"][0], (b[1] + b[3]) / 2)]
        return [(t["origin"][0], t["origin"][1]), ((b[0] + b[2]) / 2, t["origin"][1]), (b[2], t["origin"][1])]
    if is_key:   # the key entries move sideways with the widened, centre-anchored legend box: their x is proved by the block-centre check below, their y by the declared shift
        dev = abs(m["origin"][1] - (s["origin"][1] + KEY_SHIFT_PT)); key_b.append(s); key_n.append(m)
    else:
        dev = min(max(abs(x1 - x0), abs(y1 - y0)) for (x1, y1), (x0, y0) in zip(anchors(m), anchors(s)))
    worst = max(worst, dev); devs.append((s["text"], round(dev, 2)))
    pinned += 1
    if m["font"] != s["font"] or abs(m["size"] - (s["size"] + PT_PLUS)) > 0.01:
        font_mismatch.append((s["text"], s["font"], s["size"], m["font"], m["size"]))
C.check(f"every base string other than the y tick labels found on the new sheet with an anchor (left origin, centre or right edge, or x and centre for a rotated title) within {LAYOUT_TOL:g} pt of the base anchor (the matplotlib re-layout for the larger font; the key strings {-KEY_SHIFT_PT:g} pt higher, declared)",
        not unmatched and worst <= LAYOUT_TOL, f"{pinned} strings pinned, worst deviation {worst:.3f} pt, unmatched {unmatched}; deviations {sorted(devs, key=lambda d: -d[1])[:8]}")
kc_b = (min(k["bbox"][0] for k in key_b) + max(k["bbox"][2] for k in key_b)) / 2 if key_b else None; kc_n = (min(k["bbox"][0] for k in key_n) + max(k["bbox"][2] for k in key_n)) / 2 if key_n else None
def key_order(ks):
    """The key strings in reading order: the title row (the block's top) then the entry row, left to right."""
    t0 = min(k["bbox"][1] for k in ks); return [k["text"] for k in sorted(ks, key=lambda k: (k["bbox"][1] > t0 + 8, k["bbox"][0]))]
KEY_TEXT_OFFSET = (1.4 + 0.5) * PT_PLUS / 2.0   # the legend box stays anchored at the page centre; each entry's handle (1.4 em) and pad (0.5 em) grow with the font, so the strings' outer edges move right by half of that: 0.95 pt
C.check(f"the shared key stays centred as a block: its strings' outer edges sit {KEY_TEXT_OFFSET:.2f} pt right of the base's (within 0.3 pt: the handle and its pad, 1.9 em, one point larger, in a legend box anchored at the page centre), entries in the same order",
        key_b and abs((kc_n - kc_b) - KEY_TEXT_OFFSET) <= 0.3 and key_order(key_b) == key_order(key_n),
        f"key centre base {kc_b} new {kc_n}, {len(key_b)} strings")
C.check(f"fonts identical and every size one point larger for every pinned string (V13 sizes plus {PT_PLUS:g})", not font_mismatch, str(font_mismatch))
C.check("no string on the new sheet that is not on the base other than the declared y tick labels",
        sorted(s["text"] for s in sn if s["text"].strip() and nearest(s, sb) is None) == sorted(t for t in new_ytick_texts if t not in base_ytick_texts) * 2,
        str(sorted(s["text"] for s in sn if s["text"].strip() and nearest(s, sb) is None)))

# y tick labels: declared change, pinned to the new calibration the way the base pins its own
base_yl = [s for s in sb if is_ytick_label(s)]
new_yl = [s for s in sn if is_ytick_label(s)]
C.check("y tick label count: base = the V13 tick set per panel, new = the drawn tick set per panel",
        len(base_yl) == 2 * len(BASE_YTICKS) and len(new_yl) == 2 * len(new_yticks) and sorted({s["text"] for s in new_yl}, key=int) == new_ytick_texts,
        f"base {len(base_yl)} new {len(new_yl)} texts {sorted({s['text'] for s in new_yl}, key=int)}")


def label_offsets(spans, ticks_page_by_panel, ticks_data, panel_x):
    """Right edge and (label centre minus tick mark) offsets for the y tick labels of one panel."""
    out = []
    for s in spans:
        if not (panel_x[0] <= s["bbox"][0] <= panel_x[1]):
            continue
        v = float(s["text"])
        k = ticks_data.index(v)
        ty = ticks_page_by_panel[k]
        out.append((v, round(s["bbox"][2], 3), round((s["bbox"][1] + s["bbox"][3]) / 2 - ty, 3)))
    return out


# calibrate() lists the page y of the tick marks bottom-up, which is the tick set ascending
off_base = label_offsets(base_yl, rb["calib_a"]["y_ticks_page"], [float(t) for t in BASE_YTICKS], (0, 60)) + label_offsets(base_yl, rb["calib_b"]["y_ticks_page"], [float(t) for t in BASE_YTICKS], (300, 330))
off_new = label_offsets(new_yl, rn["calib_a"]["y_ticks_page"], new_yticks, (0, 60)) + label_offsets(new_yl, rn["calib_b"]["y_ticks_page"], new_yticks, (300, 330))
right_base = sorted({o[1] for o in off_base})
right_new = sorted({o[1] for o in off_new})
centre_base = sorted({o[2] for o in off_base})
centre_new = sorted({o[2] for o in off_new})
C.check("y tick labels (declared change): right edges as the base (within 0.3) and each label centred on its tick mark as the base does (within 0.3)",
        len(right_base) == len(right_new) == 2 and all(abs(a - b) <= PIN_TOL for a, b in zip(right_base, right_new))
        and max(abs(c - centre_base[0]) for c in centre_new) <= PIN_TOL,
        f"right edges base {right_base} new {right_new}, label centre minus tick base {centre_base} new {centre_new}")
Y_UNCHANGED = [int(v) for v in new_yticks] == BASE_YTICKS and abs(D["y_axis"][0] - BASE_YTICKS[0]) <= 0.01 and abs(D["y_axis"][1] - BASE_YTICKS[-1]) <= 0.01
C.check("y axis = the V13 design (range first to last V13 tick, the V13 tick set) or the change is declared in the drawn JSON",
        Y_UNCHANGED or D.get("y_axis_rule", "").startswith("DECLARED"), D.get("y_axis_rule", "no y_axis_rule in the drawn JSON"))
if Y_UNCHANGED:
    C.info(f"y axis unchanged from V13: {D['y_axis'][0]:g} to {D['y_axis'][1]:g}, ticks {[int(v) for v in new_yticks]} (no printed change); "
           f"the round-30 data rule alone would give {D.get('y_axis_data_rule')} ticks {D.get('y_ticks_data_rule')} on the data span {D.get('y_data_span')}")
else:
    C.info(f"DECLARED printed change: y ticks base {BASE_YTICKS} -> new {[int(v) for v in new_yticks]} (axis {D['y_axis'][0]:g} to {D['y_axis'][1]:g}): {D.get('y_axis_rule')}")

# ------------------------------------------------------------------ 3 read-back against the summary
def maxdev(pairs):
    return max(abs(a - b) for a, b in pairs) if pairs else 0.0


pairs_a_marker, pairs_a_line, pairs_a_q1, pairs_a_q3, pairs_b_marker, pairs_b_q1, pairs_b_q3 = [], [], [], [], [], [], []
missing = []
for k in BANDS:
    for i in range(10):
        r = S["deciles_by_band"][k][i]
        for pool, val, what in ((pairs_a_marker, rn["a_marker"][k][i], "marker"), (pairs_a_line, rn["a_line"][k][i], "line")):
            if val is None:
                missing.append((k, i + 1, what))
            else:
                pool.append((val, r["median"]))
        q = rn["a_poly"][k][i]
        if q[0] is None:
            missing.append((k, i + 1, "polygon"))
        else:
            pairs_a_q1.append((q[0], r["q1"]))
            pairs_a_q3.append((q[1], r["q3"]))
    for j, st in enumerate(STAGES):
        r = S["stages_by_band"][k][st]
        if rn["b_marker"][k][j] is None:
            missing.append((k, st, "marker"))
        else:
            pairs_b_marker.append((rn["b_marker"][k][j], r["median"]))
        w = rn["b_whisker"][k][j]
        if w[0] is None:
            missing.append((k, st, "whisker"))
        else:
            pairs_b_q1.append((w[0], r["q1"]))
            pairs_b_q3.append((w[1], r["q3"]))
C.check("every drawn element found on the new sheet: 40 + 40 + 40 panel a markers, lines, polygon columns and 20 + 20 panel b markers and whiskers",
        not missing and len(pairs_a_marker) == 40 and len(pairs_a_line) == 40 and len(pairs_a_q1) == 40 and len(pairs_b_marker) == 20 and len(pairs_b_q1) == 20,
        f"missing {missing}")
for name, pool in (("panel a marker centres = summary decile medians (40)", pairs_a_marker),
                   ("panel a median-line vertices = summary decile medians (40)", pairs_a_line),
                   ("panel a polygon lower edges = summary decile q1 (40)", pairs_a_q1),
                   ("panel a polygon upper edges = summary decile q3 (40)", pairs_a_q3),
                   ("panel b marker centres = summary stage medians (20)", pairs_b_marker),
                   ("panel b whisker lower ends = summary stage q1 (20)", pairs_b_q1),
                   ("panel b whisker upper ends = summary stage q3 (20)", pairs_b_q3)):
    C.check(name + " within 0.02 pp", maxdev(pool) <= READ_TOL, f"max abs deviation {maxdev(pool):.4f} pp over {len(pool)} values")
C.check("drawn JSON values equal the summary (medians and quartiles, 1e-9)",
        all(abs(D["panel_a_deciles"][k][f][i] - S["deciles_by_band"][k][i][f]) <= 1e-9 for k in BANDS for i in range(10) for f in ("median", "q1", "q3"))
        and all(abs(D["panel_b_stages"][k][f][j] - S["stages_by_band"][k][st][f]) <= 1e-9 for k in BANDS for j, st in enumerate(STAGES) for f in ("median", "q1", "q3"))
        and all(D["panel_a_deciles"][k]["n"][i] == S["deciles_by_band"][k][i]["n"] for k in BANDS for i in range(10))
        and all(D["panel_b_stages"][k]["n"][j] == S["stages_by_band"][k][st]["n"] for k in BANDS for j, st in enumerate(STAGES)))

# ------------------------------------------------------------------ 4 multisets
wt_b, wt_n = L.word_tokens(wb), L.word_tokens(wn)
dw = L.delta(wt_b, wt_n)
expected_lost = collections.Counter({str(t): 2 for t in BASE_YTICKS if float(t) not in new_yticks})
expected_gained = collections.Counter({t: 2 for t in new_ytick_texts if int(t) not in BASE_YTICKS})
C.check("word multiset identical to the base except the declared y tick labels (dedupe applied)",
        collections.Counter(dw["lost"]) == expected_lost and collections.Counter(dw["gained"]) == expected_gained, f"lost {dw['lost']} gained {dw['gained']}")
nt_b, nt_n = L.num_tokens(wb), L.num_tokens(wn)
dn = L.delta(nt_b, nt_n)
C.check("numeric-token multiset identical to the base except the declared y tick labels",
        collections.Counter(dn["lost"]) == expected_lost and collections.Counter(dn["gained"]) == expected_gained, f"lost {dn['lost']} gained {dn['gained']}")
C.check("no wording change: every non-numeric word of the base is on the new sheet with the same count",
        {w: c for w, c in collections.Counter(wt_b).items() if not w.replace(".", "").isdigit()} == {w: c for w, c in collections.Counter(wt_n).items() if not w.replace(".", "").isdigit()})

# ------------------------------------------------------------------ 6 house words
bad = [s["text"] for s in sn if not L.house_ok(s["text"]) or "×" in s["text"]]
C.check("house words on every string of the sheet (no em-dash, semicolon, multiplication sign, banned words)", not bad, str(bad))
C.check("panel letters lowercase bold, the only bold strings", [s["text"] for s in sn if s["font"] == "Arial-BoldMT"] == ["a", "b"])
C.check("text sizes on the new sheet = the V13 sizes plus one point (10, 11, 13 -> 11, 12, 14)", sorted({s["size"] for s in sn}) == sorted({round(v + PT_PLUS, 2) for v in {s["size"] for s in sb}}), f"new {sorted({s['size'] for s in sn})} base {sorted({s['size'] for s in sb})}")

# ------------------------------------------------------------------ 7 colours
C.check("key handle colours in the same order as the base (the light-blue T90 ladder)", rn["legend_ladder_order"] == rb["legend_ladder_order"] == BANDS,
        f"base {rb['legend_ladder_order']} new {rn['legend_ladder_order']}")
C.check("drawn JSON band colours are the base sheet's ladder", [D["band_colours"][k] for k in BANDS] == D["base_sheet_ladder"] == RB.LADDER)

# ------------------------------------------------------------------ 8 counts
import pandas as _pd
n_excl_rows = int(_pd.read_csv(L.hydrated(f"{L.SV}/oxygen_profile/oxyprofile_per_patient.csv"), usecols=["status", "_v8_stage_collapsed"]).query("status == 'ok'")["_v8_stage_collapsed"].astype(int).sum())   # v8.2: the flag column
C.check("band sizes drawn = summary band_n_cohort = v7 recompute, sum = cohort", D["band_n"] == S["band_n_cohort"] == D["band_n_v7_recomputed"] and sum(D["band_n"].values()) == S["cohort_n"] == D["cohort_n"],
        f"{D['band_n']} sum {sum(D['band_n'].values())}")
C.check("stage panel n = summary stage_panel_n = n_ok minus the exclusion list's row count",
        D["stage_panel_n"] == S["stage_panel_n"] == D["n_ok"] - n_excl_rows and D["stage_panel_excluded_n"] == S["stage_panel_excluded_n"] == n_excl_rows,
        f"stage panel {D['stage_panel_n']}, excluded {D['stage_panel_excluded_n']}, list rows {n_excl_rows}, n_ok {D['n_ok']}")
C.check("panel a frame = summary n_complete_decile_profile, per band = band_n_complete_profile",
        D["n_complete_decile_profile"] == S["n_complete_decile_profile"] and D["band_n_complete_profile"] == S["band_n_complete_profile"],
        f"{D['n_complete_decile_profile']} per band {D['band_n_complete_profile']}, collapsed nights kept in panel a {D['collapsed_nights_kept_in_panel_a']}")
C.check("per-patient file T90 = v7 table T90 for every record (max abs diff)", D["t90_per_patient_file_vs_v7_max_abs_diff"] <= 1e-9, f"{D['t90_per_patient_file_vs_v7_max_abs_diff']:.2e}")

# ------------------------------------------------------------------ 5 raster and crops
r150 = L.render(NEW, PNG150, 150)
C.check("150 dpi render written", os.path.exists(PNG150) and r150["w"] > 0, f"{r150['w']} x {r150['h']} px")
CLIPS = {"key_row": (100, 10, 400, 48), "panel_a": (8, 40, 285, 275), "panel_b": (285, 40, 484.72, 275)}
for tag, clip in CLIPS.items():
    L.render(BASE, f"{CROPS}/{tag}_OLD_200dpi.png", 200, clip=clip)
    L.render(NEW, f"{CROPS}/{tag}_NEW_200dpi.png", 200, clip=clip)
C.check("OLD vs NEW crops at 200 dpi written (key row, panel a, panel b)", all(os.path.getsize(f"{CROPS}/{t}_{v}_200dpi.png") > 0 for t in CLIPS for v in ("OLD", "NEW")))

# ------------------------------------------------------------------ CHANGES csv (old = the base sheet's drawings)
def f2(v):
    return L.r2(v)


rows = []
dramatic_rows = []
NOTE_OLD = "old value read from the base sheet's drawings through its axis calibration (the pre-v7 summary file was rewritten in place, the base drawings are the old record)"
order_changes = []
for i in range(10):
    old_order = sorted(BANDS, key=lambda k: -rb["a_marker"][k][i])
    new_order = sorted(BANDS, key=lambda k: -S["deciles_by_band"][k][i]["median"])
    if old_order != new_order:
        order_changes.append(("a", f"decile {i + 1}", old_order, new_order))
for j, st in enumerate(STAGES):
    old_order = sorted(BANDS, key=lambda k: -rb["b_marker"][k][j])
    new_order = sorted(BANDS, key=lambda k: -S["stages_by_band"][k][st]["median"])
    if old_order != new_order:
        order_changes.append(("b", st, old_order, new_order))
for k in BANDS:
    for i in range(10):
        r = S["deciles_by_band"][k][i]
        old_m, old_q = rb["a_marker"][k][i], rb["a_poly"][k][i]
        drama = abs(r["median"] - old_m) > 0.5 or any(o[0] == "a" and o[1] == f"decile {i + 1}" for o in order_changes)
        rows.append(dict(sheet=NAME, panel="a", label=f"{BAND_LABEL[k]} decile {i + 1} median", old_value=f2(old_m), new_value=f2(r["median"]),
                         source_file=SUMMARY, key=f"deciles_by_band.{k}[{i}].median", dramatic="yes" if drama else "no",
                         note=f"IQR old {f2(old_q[0])} to {f2(old_q[1])}, new {f2(r['q1'])} to {f2(r['q3'])}, new n {r['n']}. {NOTE_OLD}"))
    for j, st in enumerate(STAGES):
        r = S["stages_by_band"][k][st]
        old_m, old_w = rb["b_marker"][k][j], rb["b_whisker"][k][j]
        drama = abs(r["median"] - old_m) > 0.5 or any(o[0] == "b" and o[1] == st for o in order_changes)
        rows.append(dict(sheet=NAME, panel="b", label=f"{BAND_LABEL[k]} {st} median (IQR)", old_value=f"{f2(old_m)} ({f2(old_w[0])} to {f2(old_w[1])})",
                         new_value=f"{f2(r['median'])} ({f2(r['q1'])} to {f2(r['q3'])})", source_file=SUMMARY, key=f"stages_by_band.{k}.{st}",
                         dramatic="yes" if drama else "no",
                         note=f"new n {r['n']} on the ok records minus the {D['stage_panel_excluded_n']} collapsed-staging nights, old frame all ok records. {NOTE_OLD}"))
# stage ordering per band
for k in BANDS:
    old_lo = STAGES[int(min(range(5), key=lambda j: rb["b_marker"][k][j]))]
    old_hi = STAGES[int(max(range(5), key=lambda j: rb["b_marker"][k][j]))]
    new_lo = min(STAGES, key=lambda st: S["stages_by_band"][k][st]["median"])
    new_hi = max(STAGES, key=lambda st: S["stages_by_band"][k][st]["median"])
    rows.append(dict(sheet=NAME, panel="b", label=f"{BAND_LABEL[k]} lowest stage", old_value=old_lo, new_value=new_lo, source_file=SUMMARY,
                     key=f"stages_by_band.{k}.*.median", dramatic="yes" if old_lo != new_lo else "no",
                     note="stage with the lowest median, old from the base drawings"))
    rows.append(dict(sheet=NAME, panel="b", label=f"{BAND_LABEL[k]} highest stage", old_value=old_hi, new_value=new_hi, source_file=SUMMARY,
                     key=f"stages_by_band.{k}.*.median", dramatic="yes" if old_hi != new_hi else "no",
                     note="stage with the highest median, old from the base drawings"))
# band order across bands
rows.append(dict(sheet=NAME, panel="a and b", label="band order (highest to lowest median) at every decile and stage",
                 old_value="unchanged" if not order_changes else ", ".join(f"{p} {w}: {' > '.join(BAND_LABEL[x] for x in o)}" for p, w, o, n in order_changes),
                 new_value="≤1% > >1 to 5% > >5 to 10% > >10% everywhere" if not order_changes else ", ".join(f"{p} {w}: {' > '.join(BAND_LABEL[x] for x in n)}" for p, w, o, n in order_changes),
                 source_file=SUMMARY, key="deciles_by_band, stages_by_band", dramatic="yes" if order_changes else "no",
                 note="a band's curve changing order against another band anywhere is the dramatic rule"))
# band sizes: old from the 2026-08-25 builder's BAND_N_PUBLISHED literal (read from that file, not typed here)
old_line = [ln_ for ln_ in open(L.hydrated(OLD_BUILDER)) if ln_.startswith("BAND_N_PUBLISHED")][0]
old_band_n = json.loads(old_line.split("=", 1)[1].strip().replace("'", '"'))
for k in BANDS:
    o, n = old_band_n[k], S["band_n_cohort"][k]
    rows.append(dict(sheet=NAME, panel="text legend (band sizes)", label=f"band size {BAND_LABEL[k]}", old_value=f"{o:,}", new_value=f"{n:,}",
                     source_file=SUMMARY, key=f"band_n_cohort.{k}", dramatic="yes" if abs(n - o) / o > 0.2 else "no",
                     note=f"{100 * (n - o) / o:+.1f}%. old: BAND_N_PUBLISHED literal of the 2026-08-25 builder (the 2026-08-21 banding). new: recomputed from {D['sources']['v7_table']['path'].split('/T90_Manuscript/')[-1]} through cohort_spec.apply_cohort, asserted equal to the per-patient file and the summary, on all {S['n_ok']:,} ok records"))
rows.append(dict(sheet=NAME, panel="text legend (band sizes)", label="band sizes sum", old_value=f"{sum(old_band_n.values()):,}", new_value=f"{sum(S['band_n_cohort'].values()):,}",
                 source_file=SUMMARY, key="band_n_cohort", dramatic="no", note="all ok records, both rounds"))
# frames
old_dec_n = None
for ln_ in open(L.hydrated(OLD_LEGEND_CHECK)):
    if "records with all ten spans" in ln_:
        old_dec_n = int(ln_.split("legend")[1].split()[0].replace(",", ""))
rows.append(dict(sheet=NAME, panel="a", label="decile-frame n (records with all ten decile means)", old_value=f"{old_dec_n:,}" if old_dec_n else "not on disk",
                 new_value=f"{S['n_complete_decile_profile']:,}", source_file=SUMMARY, key="n_complete_decile_profile",
                 dramatic="yes" if old_dec_n and abs(S["n_complete_decile_profile"] - old_dec_n) / old_dec_n > 0.2 else "no",
                 note=f"old from oxygen_profile/legend_check_output.txt (2026-08-21 legend check, records with all ten spans). Collapsed-staging nights stay in panel a: {D['collapsed_nights_kept_in_panel_a']} of them have a complete profile"))
rows.append(dict(sheet=NAME, panel="b", label="stage-frame n", old_value=f"{S['n_ok']:,}", new_value=f"{S['stage_panel_n']:,}", source_file=SUMMARY, key="stage_panel_n",
                 dramatic="yes" if abs(S["stage_panel_n"] - S["n_ok"]) / S["n_ok"] > 0.2 else "no",
                 note=f"old: all ok records (the 2026-08-25 builder computed the stage panel on every ok record). new: ok minus the {S['stage_panel_excluded_n']} nights of {EXCL} (decision 6)"))
if not Y_UNCHANGED:
    rows.append(dict(sheet=NAME, panel="a and b", label="y axis (printed ticks)", old_value=f"{BASE_YTICKS[0]} to {BASE_YTICKS[-1]}, ticks {' '.join(str(t) for t in BASE_YTICKS)}",
                     new_value=f"{D['y_axis'][0]:g} to {D['y_axis'][1]:g}, ticks {' '.join(new_ytick_texts)}", source_file=DRAWN, key="y_axis, y_ticks", dramatic="no",
                     note=D.get("y_axis_rule", "declared printed change (the builder's data rule)")))
n_rows = L.write_changes(CHANGES, rows)
dramatic_rows = [r for r in rows if r["dramatic"] == "yes"]
C.info(f"CHANGES csv: {n_rows} rows, {len(dramatic_rows)} dramatic")
for r in dramatic_rows:
    C.info(f"DRAMATIC {r['panel']} | {r['label']} | old {r['old_value']} | new {r['new_value']}")
moves = sorted([(abs(S['deciles_by_band'][k][i]['median'] - rb['a_marker'][k][i]), f"a {BAND_LABEL[k]} decile {i + 1}") for k in BANDS for i in range(10)]
               + [(abs(S['stages_by_band'][k][st]['median'] - rb['b_marker'][k][j]), f"b {BAND_LABEL[k]} {st}") for k in BANDS for j, st in enumerate(STAGES)], reverse=True)
C.info("largest median moves (pp): " + ", ".join(f"{w} {v:.2f}" for v, w in moves[:5]))
band_sizes_text = " / ".join(f"{S['band_n_cohort'][k]:,}" for k in BANDS)
C.info(f"legend sentence: {D['legend_sentence_b']}. Band sizes {band_sizes_text} on all {S['n_ok']:,} ok records")

# v8.2 (V14): the printed-values csv, one row per string of the sheet. Supp 2 prints no data value (medians and bands are drawn, not
# printed): every string is an axis tick (the builder's data-driven y ticks, recorded with rule tick), a label, a key entry or a letter.
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND37_2026-09-14/figures/_common"); import v14lib as _V
_recs = []
for _s in L.dedupe(L.spans_of(NEW)):
    _t = _s["text"].replace("\xa0", " ").strip()
    if not _t: continue
    if re.fullmatch(r"\d+", _t) and (_s["bbox"][2] < 60 or 300 < _s["bbox"][0] < 330) and float(_t) in [float(v) for v in D["y_ticks"]]:
        _recs.append(dict(text=_t, x=_s["origin"][0], baseline=_s["origin"][1], ha="left", size=_s["size"], source_file=DRAWN, source_key=f"y_ticks/{[float(v) for v in D['y_ticks']].index(float(_t))}", source_value=float(_t), rule="tick", panel="a" if _s["bbox"][2] < 60 else "b"))
    else: _recs.append(dict(text=_t, x=_s["origin"][0], baseline=_s["origin"][1], ha="left", size=_s["size"], source_file="static:V13", source_key="V13 string (no data value is printed on this sheet)", source_value=_t, rule="text"))
_V.RULES["tick"] = lambda v: f"{float(v):g}"
_n, _no = _V.printed_values_csv(NEW, _recs, f"{SHEET}/verify/{NAME}_printed_values.csv", static_from=BASE, group_p_boxes=[])
C.check(f"printed values csv: {_n} strings ({len(_recs)} recorded: ticks, labels, key, letters; no data value is printed on this sheet), {_no} NO", _no == 0, f"{SHEET}/verify/{NAME}_printed_values.csv")
_yt = sorted({int(r["text"]) for r in _recs if r["rule"] == "tick"})
C.check("the y tick strings on the sheet are the builder's y_ticks", _yt == sorted(int(v) for v in D["y_ticks"]), f"{_yt} vs {D['y_ticks']}")

ok = C.write(CHECKS)
sys.exit(0 if ok else 1)
