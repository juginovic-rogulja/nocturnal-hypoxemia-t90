#!/usr/bin/env python3
"""Round 49 LED-B copy (2026-09-26, +1 pt rule): every text one point larger (l5a14.Sheet); the q/P column moves right so that the wider HR (95% CI) strings keep the V13 gap; the legend moves right with it, is drawn at 10 pt (9 + 1) with its third entry on three lines; row labels take the largest size at or above their V26 size (0.5 pt steps) that stays 4 pt clear of the row's own marks and the reference rule; the titles are centred on their V13 centres (panel b's title floored at the label x). Otherwise the round-40 builder. ED Fig 9, round 40 (lane LED, 2026-09-18): the round-37 builder (V14_L5a_PAP_MAIN/scripts/09_build_ed9.py, copied beside as
09_build_ed9_ORIGINAL_COPY.py, sha256 in ORIGINALS_sha256.txt) with ONE design change, Alen's item 7 / brief item 2: the V13 key that sat
beside panel a (13 strings at x 438 and marker glyphs replayed from the V13 drawing items at x 415 to 437, whose replayed stroke widths of
6.7 to 14.3 pt hid the square under a white stroke, blobbed the reference circle and filled the "open" square) is dropped, and a legend is
drawn INSIDE panel b (top right of b) by the same Shape helpers and the same marker geometry and colours as the plotted series: a filled
square on a line (the >5 to 10% series), a filled triangle on a line (the >10% series), the open grey reference circle, and an open square
beside an open triangle for "the interval includes 1". Numbers: unchanged (same v8 files, same printing rules, same row geometry; proved
against the delivered V16 sheet by verify_sheet.py). Output: work/ED_Fig09_rebuilt.pdf at the V16 page size (the 16 pt title strip is added
afterwards by stamp_title.py), work/legend_record.json, work/build/ (drawn records, rows, build_checks.txt, printed-values csv)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from l5a14 import *

SHEET = "ED_Fig09"; MYLANE = LANE; R40LANE = f"{paths.FIGURE_ROOT}/ROUND40_2026-09-18/lanes/LED"   # round 49: this lane; the round-40 lane read-only
D = f"{MYLANE}/{SHEET}"; W = f"{D}/work"; WB = f"{W}/build"; os.makedirs(WB, exist_ok=True)
BASE = f"{V13}/{SHEET}.pdf"; V16 = f"{paths.FIGURE_ROOT}/ROUND38_2026-09-16/figures/NEW_FINAL_SET_V16/{SHEET}.pdf"; OUT = f"{W}/{SHEET}_rebuilt.pdf"
R37 = load_json(f"{R37LANE}/{SHEET}/verify/{SHEET}_readback.json")["provenance"]
ck = Checks(f"Lane LED-B round 49 (+{PT_PLUS:g} pt), {SHEET} rebuild, {time.strftime('%Y-%m-%d %H:%M')}. BASE (design) = V13 {BASE} (sha256 {sha256(BASE)[:16]}); printed-number reference = V16 {V16} (sha256 {sha256(V16)[:16]}). NEW = results_continuous.csv and results_graded.csv (step 147).")
sc_c = gate("results_continuous"); sc_g = gate("results_graded"); dump({"continuous": sc_c, "graded": sc_g}, f"{WB}/sidecars.json")
ck.log("results_continuous.csv and results_graded.csv pass the v8 sidecar gate (step 147, inputs data_frozen_v8_2026-09, sha equal)", sc_c["v8_inputs"] and sc_g["v8_inputs"], f"out {sc_c['output_mtime']} / {sc_g['output_mtime']}")
ck.log("the number files are the ones the V16 sheet was built from (sha256 equal to the round-37 lane's readback provenance)", sc_c["sha256"] == R37["continuous"]["sha256"] and sc_g["sha256"] == R37["graded"]["sha256"], f"{sc_c['sha256'][:16]} {sc_g['sha256'][:16]}")

# ---- 1 V13 geometry (verbatim from the round-37 builder)
spans, info = text_spans(BASE); spans = collapse(spans); items = drawings(BASE)
dump({"info": info, "spans": spans}, f"{WB}/base_text.json"); dump(items, f"{WB}/base_drawings.json")
PW, PH = info["rect"][2], info["rect"][3]
ck.log("V13 sheet flat (0 XObjects), page 518.74 x 730.60", info["xobjects"] == 0 and abs(PW - 518.74) < 0.05 and abs(PH - 730.6) < 0.05, info)
letters = sorted([s for s in spans if len(s["text"]) == 1 and s["text"] in "ab" and s["size"] > 12.5], key=lambda s: s["y"]); assert [l["text"] for l in letters] == ["a", "b"]
M = markers(items)
labA = column(spans, 10, 16, ymin=30, ymax=200, size=10.0); assert len(labA) == 11
A_Y = [s["y"] for s in labA]; A_PITCH = float(np.mean(np.diff(A_Y))); A_X = float(np.mean([s["x"] for s in labA]))
datA = sorted([m for m in M if m["shape"] == "circle" and not m["open"] and m["colour"] == BLUE and abs(m["w"] - 5.29) < 0.05 and m["cy"] < 200 and m["cx"] < 250], key=lambda m: m["cy"]); assert len(datA) == 11
OFF_A = float(np.mean([m["cy"] - s["y"] for m, s in zip(datA, labA)])); MARK_D = float(np.mean([m["w"] for m in datA])); MARK_EDGE = float(np.mean([m["widths"][0] for m in datA]))
hrA = column(spans, 250, 258, ymin=40, ymax=200); qA = column(spans, 320, 335, ymin=40, ymax=200); assert len(hrA) == 11 and len(qA) == 11, (len(hrA), len(qA))
HR_X = float(np.mean([s["x"] for s in hrA])); Q_X = float(np.mean([s["x"] for s in qA])); COL_OFF_A = float(np.mean([h["y"] - s["y"] for h, s in zip(hrA, labA)]))
ciA = sorted([h for h in hlines(items, min_len=0.5) if h["stroke"] == BLUE and h["width"] and h["width"] > 1.5 and h["x0"] < 250 and 30 <= h["y"] <= 195], key=lambda h: h["y"]); assert len(ciA) == 11; CI_W = float(np.mean([h["width"] for h in ciA]))
spines = sorted([h for h in hlines(items, min_len=100) if h["stroke"] == INK], key=lambda h: h["y"]); assert len(spines) == 2; SPA, SPB = spines
refls = sorted([it for it in items if "line" in it and abs(it["line"][0] - it["line"][2]) < 0.05 and abs(it["line"][3] - it["line"][1]) > 100 and hexcol(it["color"]) == INK], key=lambda it: min(it["line"][1], it["line"][3])); assert len(refls) == 2; REFA, REFB = refls
ticksA = sorted([v for v in vlines(items, max_len=6.0) if v["stroke"] == INK and abs(v["y0"] - SPA["y"]) < 0.05], key=lambda v: v["x"]); ticksB = sorted([v for v in vlines(items, max_len=6.0) if v["stroke"] == INK and abs(v["y0"] - SPB["y"]) < 0.05], key=lambda v: v["x"])
assert len(ticksA) == 4 and len(ticksB) == 6; TICK_LEN = float(np.mean([v["y1"] - v["y0"] for v in ticksA + ticksB]))
TVA, TVB = [1, 1.5, 2, 3], [0.3, 0.5, 1, 2, 4, 8]
aA, bA, rA = fit_axis(TVA, [v["x"] for v in ticksA]); aB, bB, rB = fit_axis(TVB, [v["x"] for v in ticksB]); assert max(np.abs(rA).max(), np.abs(rB).max()) < 0.05
tlA = sorted([s for s in spans if abs(s["y"] - 198.08) < 0.6], key=lambda s: s["x"]); tlB = sorted([s for s in spans if abs(s["y"] - 701.06) < 0.6], key=lambda s: s["x"]); assert [s["text"] for s in tlA] == ["1", "1.5", "2", "3"] and [s["text"] for s in tlB] == ["0.3", "0.5", "1", "2", "4", "8"]
TL_OFF = float(np.mean([s["y"] for s in tlB])) - SPB["y"]; TL_OFF_A = float(np.mean([s["y"] for s in tlA])) - SPA["y"]
titleA = find_span(spans, "Hazard ratio per 1 SD of residual sleep T90 on PAP (95% CI)"); titleB = find_span(spans, "Hazard ratio against residual T90 at or below 5% of the night (95% CI)"); assert titleA and titleB
headA = [find_span(spans, "HR (95% CI)", ymin=25, ymax=40), find_span(spans, "q", ymin=25, ymax=40)]; headB = [find_span(spans, "HR (95% CI)", ymin=235, ymax=245), find_span(spans, "P", ymin=235, ymax=245)]; assert all(headA) and all(headB)
key_txt = [s for s in spans if s["x"] > 430 and s["y"] < 170]; assert len(key_txt) == 13, len(key_txt)      # the V13 key beside panel a: read for the record, NOT drawn (round 40)
labB = column(spans, 10, 16, ymin=250, ymax=700, size=10.0); assert len(labB) == 17; B_Y = [s["y"] for s in labB]; B_PITCH = float(np.mean(np.diff(B_Y))); assert abs(B_PITCH - 26.0) < 0.01
sq = sorted([m for m in M if m["shape"] == "square" and m["cy"] > 245 and abs(m["w"] - 4.8) < 0.1 and m["cx"] < 250], key=lambda m: m["cy"]); tr = sorted([m for m in M if m["shape"] == "triangle" and m["cy"] > 245 and abs(m["w"] - 5.48) < 0.1 and m["cx"] < 250], key=lambda m: m["cy"]); refB = sorted([m for m in M if m["shape"] == "circle" and m["open"] and m["cy"] > 245 and abs(m["w"] - 5.1) < 0.05 and m["cx"] < 250], key=lambda m: m["cy"])
assert len(sq) == len(tr) == len(refB) == 17
OFF_SQ = float(np.mean([m["cy"] - s["y"] for m, s in zip(sq, labB)])); OFF_TR = float(np.mean([m["cy"] - s["y"] for m, s in zip(tr, labB)])); OFF_REF = float(np.mean([m["cy"] - s["y"] for m, s in zip(refB, labB)]))
SQ_SIDE = float(np.mean([m["w"] for m in sq])); TR_W = float(np.mean([m["w"] for m in tr])); TR_H = float(np.mean([m["h"] for m in tr])); REF_D = float(np.mean([m["w"] for m in refB])); REF_EDGE = float(np.mean([m["widths"][0] for m in refB]))
FILL_EDGE = 0.8; OPEN_EDGE = 1.1
X_REF_B = float(np.median([m["cx"] for m in refB])); X_REF_A = REFA["line"][0]
hrB = column(spans, 250, 258, ymin=250, ymax=700); pB = column(spans, 320, 335, ymin=250, ymax=700); assert len(hrB) == 34 and len(pB) == 34
OFF_MID = float(np.mean([hrB[2 * i]["y"] - s["y"] for i, s in enumerate(labB)])); OFF_TOP = float(np.mean([hrB[2 * i + 1]["y"] - s["y"] for i, s in enumerate(labB)]))
ck.log(f"V13 geometry read: a 11 rows pitch {A_PITCH:.3f} (markers {MARK_D:.2f} pt at {OFF_A:+.2f}, columns at x {HR_X:.2f} and {Q_X:.2f}, {COL_OFF_A:+.2f}); b 17 rows pitch {B_PITCH:.2f} (square {SQ_SIDE:.2f} at {OFF_SQ:+.2f}, triangle {TR_W:.2f}x{TR_H:.2f} at {OFF_TR:+.2f}, reference {REF_D:.2f} at {OFF_REF:+.2f}, strings at {OFF_MID:+.2f} and {OFF_TOP:+.2f}); axes a {bA:.3f} and b {bB:.3f} pt per ln unit",
       abs(aA - X_REF_A) < 0.05 and abs(aB - X_REF_B) < 0.05 and abs(OFF_SQ - OFF_REF + 7.3) < 0.05 and abs(OFF_TR - OFF_REF - 7.3) < 0.05, f"aA {aA:.3f} aB {aB:.3f}")
old_labA = [s["text"] for s in labA]; old_labB = [s["text"] for s in labB]; old_hrA = [s["text"] for s in hrA]; old_qA = [s["text"] for s in qA]; old_hrB = [s["text"] for s in hrB]; old_pB = [s["text"] for s in pB]

# ---- 2 v8 rows (verbatim)
cont, A = residual_rows(); graded, B = graded_rows()
ck.log(f"panel a: {len(A)} Benjamini-Hochberg survivors (sig_fdr) of {len(cont)} outcomes, all q < 0.05 and lower CI above 1; controls in the file {list(cont[cont.negative_control == True].outcome)}", (A.q_sd < 0.05).all() and (A.lo_sd > 1.0).all(), list(A.outcome))  # noqa: E712
ck.log(f"panel b: {len(B)} principal outcomes of results_graded.csv (V13 had 17): the page grows by {(len(B) - 17) * B_PITCH:.1f} pt", len(B) >= 1, list(B.outcome))
NB = len(B); DELTA = (NB - 17) * B_PITCH; PH2 = PH + DELTA
xloA, xhiA = exp((SPA["x0"] - aA) / bA), exp((SPA["x1"] - aA) / bA); xloB, xhiB = exp((SPB["x0"] - aB) / bB), exp((SPB["x1"] - aB) / bB)
ck.log("every value lies inside the V13 axis ranges (no tick change)", all(xloA <= r.lo_sd and r.hi_sd <= xhiA for r in A.itertuples()) and all(xloB <= min(r.lo_5to10, r.lo_gt10) and max(r.hi_5to10, r.hi_gt10) <= xhiB for r in B.itertuples()), f"a {xloA:.3f}-{xhiA:.3f}, b {xloB:.3f}-{xhiB:.3f}")

# ---- 3 draw
S = Sheet(PW, PH2)
SIZE_KEEP = {}   # round 49: (panel, label) -> asked size below 10.0 where a label at 11 pt would come within 4 pt of its row's marks or the reference rule
def fit_label(panel, label, limit_x):
    for cand in (10.0, 9.5, 9.0):
        if A_X + S.width(label, cand) + 4.0 <= limit_x:
            if cand != 10.0: SIZE_KEEP[f"{panel}:{label}"] = cand + PT_PLUS
            return cand
    raise AssertionError(("row label reaches the marks even at its V26 size", panel, label, A_X + S.width(label, 9.0), limit_x))
GAP13 = Q_X - max(s["bbox"][2] for s in hrA + hrB)   # the V13 gap between the HR (95% CI) column's right edge and the q/P column
hr_texts = [hr_ci(r.hr_sd, r.lo_sd, r.hi_sd) for r in A.itertuples()] + [t for r in B.itertuples() for t in (hr_ci(r.hr_5to10, r.lo_5to10, r.hi_5to10), hr_ci(r.hr_gt10, r.lo_gt10, r.hi_gt10))]
Q_X_NEW = round(HR_X + max(S.width(t, 9.5) for t in hr_texts) + GAP13, 2)
ck.log(f"round 49: the q/P column moves from x {Q_X:.2f} to {Q_X_NEW:.2f} so that the widest HR (95% CI) string at {9.5 + PT_PLUS:g} pt keeps the V13 gap of {GAP13:.2f} pt (a text nudge, no mark moves)", Q_X_NEW > Q_X, "")
for l in letters: S.text(l["x"], l["y"], l["text"], 13.0, bold=True, panel=l["text"], source_key="panel letter")
# panel a (the key beside a is NOT drawn: round 40)
S.line(REFA["line"][0], REFA["line"][1], REFA["line"][2], REFA["line"][3], INK, REFA["width"], cap=0, dashes=REFA.get("dashes"))
rowsA = []
for i, r in enumerate(A.itertuples()):
    y = A_Y[0] + i * A_PITCH; cy = y + OFF_A
    S.line(axis_x(aA, bA, r.lo_sd), cy, axis_x(aA, bA, r.hi_sd), cy, BLUE, CI_W, cap=1); S.circle(axis_x(aA, bA, r.hr_sd), cy, MARK_D, BLUE, WHITE, MARK_EDGE)
    S.text(A_X, y, plab(r.outcome), fit_label("a", plab(r.outcome), min(axis_x(aA, bA, r.lo_sd), X_REF_A)), panel="a", source_file=SPEC["results_continuous"], source_key=f"key={r.key}:outcome", source_value=r.outcome, rule="label")
    S.text(HR_X, y + COL_OFF_A, hr_ci(r.hr_sd, r.lo_sd, r.hi_sd), 9.5, panel="a", source_file=SPEC["results_continuous"], source_key=f"key={r.key}:hr_sd,lo_sd,hi_sd", source_value=[r.hr_sd, r.lo_sd, r.hi_sd], rule="hr_ci")
    S.text(Q_X_NEW, y + COL_OFF_A, fmt_p(r.q_sd), 9.5, bold=bool(r.q_sd < 0.05), panel="a", source_file=SPEC["results_continuous"], source_key=f"key={r.key}:q_sd", source_value=r.q_sd, rule="fmt_p")
    rowsA.append(dict(key=r.key, outcome=r.outcome, hr=r.hr_sd, lo=r.lo_sd, hi=r.hi_sd, q=r.q_sd, x=axis_x(aA, bA, r.hr_sd), y=cy))
S.line(SPA["x0"], SPA["y"], SPA["x1"], SPA["y"], INK, SPA["width"], cap=2)
def tick_sizes(panel, ticks, labs):
    """Round 49: each tick label at 10 + 1 pt unless it would come within 2.5 pt of a neighbour at that size, then the pair keeps 10 pt (declared)."""
    sz = {lab["text"]: 10.0 for lab in labs}
    for (t0, l0), (t1, l1) in zip(zip(ticks, labs), zip(ticks[1:], labs[1:])):
        gap = (t1["x"] - S.width(l1["text"], 10.0) / 2) - (t0["x"] + S.width(l0["text"], 10.0) / 2)
        if gap < 2.5: sz[l0["text"]] = 9.0; sz[l1["text"]] = 9.0; SIZE_KEEP[f"{panel}:tick {l0['text']}"] = 10.0; SIZE_KEEP[f"{panel}:tick {l1['text']}"] = 10.0
    return sz
TSA = tick_sizes("a", ticksA, tlA); TSB = tick_sizes("b", ticksB, tlB)
for t, lab in zip(ticksA, tlA): S.line(t["x"], SPA["y"], t["x"], SPA["y"] + TICK_LEN, INK, t["width"], cap=0); S.text(t["x"], lab["y"], lab["text"], TSA[lab["text"]], align="center", panel="a", source_key="axis tick")
def cx13(s): return (s["bbox"][0] + s["bbox"][2]) / 2
S.text(headA[0]["x"], headA[0]["y"], headA[0]["text"], headA[0]["size"], bold=True, panel="a", source_key="static V13 string")
S.text(Q_X_NEW, headA[1]["y"], headA[1]["text"], headA[1]["size"], bold=True, panel="a", source_key="static V13 string (q header on the moved column)")
S.text(cx13(titleA), titleA["y"], titleA["text"], titleA["size"], align="center", panel="a", source_key="static V13 string (centred on its V13 centre)")
# panel b (verbatim)
S.line(REFB["line"][0], min(REFB["line"][1], REFB["line"][3]), REFB["line"][2], SPB["y"] + DELTA, INK, REFB["width"], cap=0, dashes=REFB.get("dashes"))
rowsB = []
for i, r in enumerate(B.itertuples()):
    y = B_Y[0] + i * B_PITCH; y_sq, y_tr, y_ref = y + OFF_SQ, y + OFF_TR, y + OFF_REF
    for hr, lo, hi, yy, colr, kind in ((r.hr_5to10, r.lo_5to10, r.hi_5to10, y_sq, MID, "square"), (r.hr_gt10, r.lo_gt10, r.hi_gt10, y_tr, BLUE, "triangle")):
        S.line(axis_x(aB, bB, lo), yy, axis_x(aB, bB, hi), yy, colr, CI_W, cap=1)
    S.circle(X_REF_B, y_ref, REF_D, WHITE, GREY, REF_EDGE)
    for hr, lo, hi, yy, colr, kind in ((r.hr_5to10, r.lo_5to10, r.hi_5to10, y_sq, MID, "square"), (r.hr_gt10, r.lo_gt10, r.hi_gt10, y_tr, BLUE, "triangle")):
        sig = not (lo <= 1.0 <= hi); cx = axis_x(aB, bB, hr)
        if kind == "square": S.square(cx, yy, SQ_SIDE, colr if sig else WHITE, WHITE if sig else colr, FILL_EDGE if sig else OPEN_EDGE)
        else: S.triangle(cx, yy, TR_W, TR_H, colr if sig else WHITE, WHITE if sig else colr, FILL_EDGE if sig else OPEN_EDGE)
        rowsB.append(dict(key=r.key, band=kind, hr=hr, lo=lo, hi=hi, x=cx, y=yy, filled=sig))
    S.text(A_X, y, plab(r.outcome), fit_label("b", plab(r.outcome), min(axis_x(aB, bB, min(r.lo_5to10, r.lo_gt10)), X_REF_B)), panel="b", source_file=SPEC["results_graded"], source_key=f"key={r.key}:outcome", source_value=r.outcome, rule="label", note="" if plab(r.outcome) == r.outcome else f"printed as the short form used on Fig 2 (file label: {r.outcome})")
    S.text(HR_X, y + OFF_MID, hr_ci(r.hr_5to10, r.lo_5to10, r.hi_5to10), 9.5, panel="b", source_file=SPEC["results_graded"], source_key=f"key={r.key}:hr_5to10,lo_5to10,hi_5to10", source_value=[r.hr_5to10, r.lo_5to10, r.hi_5to10], rule="hr_ci")
    S.text(HR_X, y + OFF_TOP, hr_ci(r.hr_gt10, r.lo_gt10, r.hi_gt10), 9.5, panel="b", source_file=SPEC["results_graded"], source_key=f"key={r.key}:hr_gt10,lo_gt10,hi_gt10", source_value=[r.hr_gt10, r.lo_gt10, r.hi_gt10], rule="hr_ci")
    S.text(Q_X_NEW, y + OFF_MID, fmt_p(r.p_5to10), 9.5, bold=bool(r.p_5to10 < 0.05), panel="b", source_file=SPEC["results_graded"], source_key=f"key={r.key}:p_5to10", source_value=r.p_5to10, rule="fmt_p")
    S.text(Q_X_NEW, y + OFF_TOP, fmt_p(r.p_gt10), 9.5, bold=bool(r.p_gt10 < 0.05), panel="b", source_file=SPEC["results_graded"], source_key=f"key={r.key}:p_gt10", source_value=r.p_gt10, rule="fmt_p")
S.line(SPB["x0"], SPB["y"] + DELTA, SPB["x1"], SPB["y"] + DELTA, INK, SPB["width"], cap=2)
for t, lab in zip(ticksB, tlB): S.line(t["x"], SPB["y"] + DELTA, t["x"], SPB["y"] + DELTA + TICK_LEN, INK, t["width"], cap=0); S.text(t["x"], lab["y"] + DELTA, lab["text"], TSB[lab["text"]], align="center", panel="b", source_key="axis tick")
S.text(headB[0]["x"], headB[0]["y"], headB[0]["text"], headB[0]["size"], bold=True, panel="b", source_key="static V13 string")
S.text(Q_X_NEW, headB[1]["y"], headB[1]["text"], headB[1]["size"], bold=True, panel="b", source_key="static V13 string (P header on the moved column)")
TITLE_B_CX = max(cx13(titleB), A_X + S.width(titleB["text"], titleB["size"]) / 2)   # centred on its V13 centre unless that puts its left edge left of the label x (then floored there)
S.text(TITLE_B_CX, titleB["y"] + DELTA, titleB["text"], titleB["size"], align="center", panel="b", source_key="static V13 string (moved with the axis, centred)")

# ---- 3b the legend (round 40): inside panel b, top right of b, the plotted helpers, geometry and colours
P_RIGHT = max(Q_X_NEW + S.width(t, 9.5, bold=True) for t in [t for r in B.itertuples() for t in (fmt_p(r.p_5to10), fmt_p(r.p_gt10))] + [fmt_p(r.q_sd) for r in A.itertuples()])   # right edge of the moved q/P column at +1 pt (bold is the wider face)
LEG_X_LINE0 = float(np.ceil(P_RIGHT + 6.5)); LEG_X_LINE1 = LEG_X_LINE0 + 20.0; LEG_X_MARK = (LEG_X_LINE0 + LEG_X_LINE1) / 2; LEG_X_TEXT = LEG_X_LINE1 + 6.0; LEG_X_RIGHT = PW - 4.5; LEG_GAP = 5.0   # round 49: the V13 offsets (line 20 pt, text 6 pt right of it) from the moved column
LEG_Y0 = headB[0]["y"]        # the first legend baseline on panel b's column-header line
# the line breaks are set by hand at the phrase boundaries ("|"): a number phrase such as "5 to 10%" is never split across lines
ENTRIES = [("Residual T90 above|5 to 10% of the night", "square_on_line"), ("Residual T90 above 10%|of the night", "triangle_on_line"),
           ("Residual T90|at or below 5%,|the reference at 1", "open_circle"), ("Open square or triangle:|the interval includes 1", "open_pair")]   # round 49: the third entry on three lines (at 10 pt its first line no longer fits the column), wording unchanged
def wrap(text, size, maxw):
    out = text.split("|"); assert all(S.width(l_, size) <= maxw for l_ in out), (text, [round(S.width(l_, size), 1) for l_ in out], maxw); return out
def fits(size): return all(S.width(l_, size) <= LEG_X_RIGHT - LEG_X_TEXT for t, _ in ENTRIES for l_ in t.split("|"))
LEG_SIZE = next((c for c in (9.5, 9.0, 8.41) if fits(c)), None); assert LEG_SIZE is not None
LEG_PITCH = round((LEG_SIZE + PT_PLUS) * 1.18, 2)   # round 49: the pitch rule on the drawn size
y = LEG_Y0; leg_lines = []; leg_handles = []
for text, kind in ENTRIES:
    ls = wrap(text, LEG_SIZE, LEG_X_RIGHT - LEG_X_TEXT); ys = [y + i * LEG_PITCH for i in range(len(ls))]
    cy = (ys[0] + ys[-1]) / 2 - 0.36 * (LEG_SIZE + PT_PLUS)
    if kind == "square_on_line": S.line(LEG_X_LINE0, cy, LEG_X_LINE1, cy, MID, CI_W, cap=1); S.square(LEG_X_MARK, cy, SQ_SIDE, MID, WHITE, FILL_EDGE)
    elif kind == "triangle_on_line": S.line(LEG_X_LINE0, cy, LEG_X_LINE1, cy, BLUE, CI_W, cap=1); S.triangle(LEG_X_MARK, cy, TR_W, TR_H, BLUE, WHITE, FILL_EDGE)
    elif kind == "open_circle": S.circle(LEG_X_MARK, cy, REF_D, WHITE, GREY, REF_EDGE)
    else: S.square(LEG_X_MARK - 5.0, cy, SQ_SIDE, WHITE, MID, OPEN_EDGE); S.triangle(LEG_X_MARK + 5.0, cy, TR_W, TR_H, WHITE, BLUE, OPEN_EDGE)
    leg_handles.append(dict(kind=kind, cx=LEG_X_MARK, cy=round(cy, 3), line=[LEG_X_LINE0, LEG_X_LINE1] if "line" in kind else None))
    for l_, yy in zip(ls, ys):
        S.text(LEG_X_TEXT, yy, l_, LEG_SIZE, panel="b", source_key="legend (round 40, inside panel b)"); leg_lines.append(dict(text=l_, x=LEG_X_TEXT, y=round(yy, 3), size=LEG_SIZE + PT_PLUS, bold=False, entry=text.replace("|", " ")))
    y = ys[-1] + LEG_PITCH + LEG_GAP
LEG_BOX = [LEG_X_LINE0 - 6.0, LEG_Y0 - LEG_SIZE - PT_PLUS - 2.0, LEG_X_RIGHT + 4.0, y - LEG_PITCH - LEG_GAP + 4.0]
ck.log(f"legend inside panel b: {len(ENTRIES)} entries, {len(leg_lines)} lines at {LEG_SIZE + PT_PLUS} pt (pitch {LEG_PITCH}), text from x {LEG_X_TEXT} to at most {LEG_X_RIGHT}, handles at x {LEG_X_LINE0} to {LEG_X_LINE1}, first baseline on the b header line y {LEG_Y0}; clear of the P column (right edge {P_RIGHT:.1f}) and inside the page",
       LEG_X_LINE0 - 6.0 > P_RIGHT and LEG_X_RIGHT + 4.0 < PW and LEG_BOX[3] < B_Y[0] + 5 * B_PITCH, f"box {[round(v, 1) for v in LEG_BOX]}")
ck.log("legend handles use the plotted marker geometry and colours (same helpers: square side, triangle w x h, reference diameter and edge, CI line width, fill/open edges)", True,
       f"square {SQ_SIDE:.2f} {MID} white edge {FILL_EDGE}; triangle {TR_W:.2f}x{TR_H:.2f} {BLUE}; reference {REF_D:.2f} white with {GREY} edge {REF_EDGE:.2f}; open pair edges {OPEN_EDGE}; line {CI_W:.2f}")
S.save(OUT); drawn = S.drawn; dump(drawn, f"{WB}/{SHEET}_drawn.json"); dump({"a": rowsA, "b": rowsB}, f"{WB}/{SHEET}_rows.json")
v16sp = load_json(f"{R40LANE}/verify/v16_spans/{SHEET}_spans.json")["spans"]
old_key = [dict(text=s["text"].replace("\xa0", " "), x=s["origin"][0], y=s["origin"][1]) for s in sorted(v16sp, key=lambda s: s["origin"][1]) if s["origin"][0] > 430 and s["origin"][1] < 170]
assert len(old_key) == 13 and [k["text"] for k in old_key] == [s["text"] for s in sorted(key_txt, key=lambda s: s["y"])], (old_key, key_txt)
dump(dict(size=LEG_SIZE + PT_PLUS, pitch=LEG_PITCH, gap=LEG_GAP, x_text=LEG_X_TEXT, x_line=[LEG_X_LINE0, LEG_X_LINE1], y0=LEG_Y0, lines=leg_lines, handles=leg_handles, legend_box=[round(v, 2) for v in LEG_BOX],
          old_key_spans=old_key, old_key_box=[412.0, 26.0, 512.0, 172.0], entries=[t.replace("|", " ") for t, _ in ENTRIES]), f"{W}/legend_record.json")

# ---- 4 verify (the round-37 checks, the key dropped from the statics; the V13 multiset delta now also carries key out / legend in)
sp2, info2 = text_spans(OUT); sp2 = collapse(sp2); items2 = drawings(OUT)
ck.log(f"new sheet: one page, width equal to V13, height {PH:.2f} -> {PH2:.2f} ({DELTA:+.1f} pt, {NB} rows in b), flat", abs(info2["rect"][2] - PW) < 0.01 and abs(info2["rect"][3] - PH2) < 0.01 and info2["xobjects"] == 0, info2)
ck.log("fonts: Arial faces only", all("Arial" in f for f in info2["fonts"]), info2["fonts"])
L2 = sorted([s for s in sp2 if len(s["text"]) == 1 and s["text"] in "ab" and s["size"] > 12.5], key=lambda s: s["y"])
ck.log("letters a and b Arial Bold 13 pt at the V13 origins within 0.3 pt", [l["text"] for l in L2] == ["a", "b"] and all(abs(l["x"] - o["x"]) < 0.3 and abs(l["y"] - o["y"]) < 0.3 and "Bold" in l["font"] for l, o in zip(L2, letters)), [(l["text"], l["x"], l["y"]) for l in L2])
def at_anchor(text, y, kind, ax, tol=0.6):
    g = find_span(sp2, text, y=y, ytol=tol)
    if g is None: return False
    got = g["x"] if kind == "left" else (g["bbox"][0] + g["bbox"][2]) / 2
    return abs(got - ax) < tol
ok_static = [at_anchor(headA[0]["text"], headA[0]["y"], "left", headA[0]["x"]), at_anchor(headA[1]["text"], headA[1]["y"], "left", Q_X_NEW), at_anchor(titleA["text"], titleA["y"], "centre", cx13(titleA)),
             at_anchor(headB[0]["text"], headB[0]["y"], "left", headB[0]["x"]), at_anchor(headB[1]["text"], headB[1]["y"], "left", Q_X_NEW)] + [at_anchor(lab["text"], lab["y"], "centre", t["x"]) for t, lab in zip(ticksA, tlA)]
moved = [at_anchor(lab["text"], lab["y"] + DELTA, "centre", t["x"]) for t, lab in zip(ticksB, tlB)] + [at_anchor(titleB["text"], titleB["y"] + DELTA, "centre", TITLE_B_CX)]
ck.log(f"static strings: 'HR (95% CI)' heads at the V13 x, the q and P heads on the moved column (x {Q_X_NEW:.2f}), panel a's title centred on its V13 centre, a's tick labels centred on the ticks; b's tick labels and axis title moved by exactly {DELTA:+.1f} pt (the title centred at x {TITLE_B_CX:.2f}, V13 centre {cx13(titleB):.2f})", all(ok_static) and all(moved), f"{sum(ok_static)} of {len(ok_static)}, moved {sum(moved)} of {len(moved)}")
ck.log("the V13 key beside panel a is gone (no string at x > 430 above y 170)", not [s for s in sp2 if s["x"] > 430 and s["y"] < 170], "")
labA2 = column(sp2, 10, 16, ymin=30, ymax=200); hrA2 = column(sp2, 250, 258, ymin=40, ymax=200); qA2 = column(sp2, Q_X_NEW - 3, Q_X_NEW + 3, ymin=40, ymax=200)
ck.log("panel a: 11 labels, HR (95% CI) and q strings read back in the file's order at the V13 offsets, q bold where q < 0.05", [s["text"] for s in labA2] == [plab(o) for o in A.outcome] and [norm(s["text"]) for s in hrA2] == [hr_ci(r.hr_sd, r.lo_sd, r.hi_sd) for r in A.itertuples()] and [s["text"] for s in qA2] == [fmt_p(r.q_sd) for r in A.itertuples()] and all(("Bold" in s["font"]) == (r.q_sd < 0.05) for s, r in zip(qA2, A.itertuples())) and all(abs(s["y"] - y) < 0.3 for s, y in zip(labA2, A_Y)), f"{len(labA2)} {len(hrA2)} {len(qA2)}")
labB2 = column(sp2, 10, 16, ymin=250, ymax=PH2 - 40); hrB2 = column(sp2, 250, 258, ymin=250, ymax=PH2 - 40); pB2 = column(sp2, Q_X_NEW - 3, Q_X_NEW + 3, ymin=250, ymax=PH2 - 40)
expB_hr = [t for r in B.itertuples() for t in (hr_ci(r.hr_5to10, r.lo_5to10, r.hi_5to10), hr_ci(r.hr_gt10, r.lo_gt10, r.hi_gt10))]; expB_p = [t for r in B.itertuples() for t in (fmt_p(r.p_5to10), fmt_p(r.p_gt10))]
expB_bold = [b_ for r in B.itertuples() for b_ in (r.p_5to10 < 0.05, r.p_gt10 < 0.05)]
ck.log(f"panel b: {NB} labels on the 26.0 pt pitch, both bands' HR (95% CI) and P strings read back from results_graded.csv in order, P bold where P < 0.05", [s["text"] for s in labB2] == [plab(o) for o in B.outcome] and [norm(s["text"]) for s in hrB2] == expB_hr and [s["text"] for s in pB2] == expB_p and all(("Bold" in s["font"]) == e for s, e in zip(pB2, expB_bold)) and all(abs(s["y"] - (B_Y[0] + i * B_PITCH)) < 0.3 for i, s in enumerate(labB2)), f"{len(labB2)} {len(hrB2)} {len(pB2)}")
ck.log("no events counts and no note on the sheet (V13 design)", not [s for s in sp2 if s["text"] == "Events" or re.fullmatch(r"\d+/\d+/\d+", s["text"])], "")
M2 = markers(items2); datA2 = sorted([m for m in M2 if m["shape"] == "circle" and not m["open"] and m["colour"] == BLUE and abs(m["w"] - MARK_D) < 0.06 and m["cy"] < 200 and m["cx"] < 250], key=lambda m: m["cy"])
sq2 = sorted([m for m in M2 if m["shape"] == "square" and m["cy"] > 245 and abs(m["w"] - SQ_SIDE) < 0.1 and m["cx"] < 250], key=lambda m: m["cy"]); tr2 = sorted([m for m in M2 if m["shape"] == "triangle" and m["cy"] > 245 and abs(m["w"] - TR_W) < 0.1 and m["cx"] < 250], key=lambda m: m["cy"]); ref2 = sorted([m for m in M2 if m["shape"] == "circle" and m["open"] and m["cy"] > 245 and abs(m["w"] - REF_D) < 0.06 and m["cx"] < 250], key=lambda m: m["cy"])
aA2, bA2, rA2 = fit_axis(list(A.hr_sd), [m["cx"] for m in datA2]) if len(datA2) == 11 else (0, 0, np.array([9]))
valsB = [v for r in B.itertuples() for v in (r.hr_5to10, r.hr_gt10)]; xsB = [m["cx"] for pair in zip(sq2, tr2) for m in pair] if len(sq2) == len(tr2) == NB else []
aB2, bB2, rB2 = fit_axis(valsB, xsB) if xsB else (0, 0, np.array([9]))
fill_ok = len(sq2) == NB and all(m["open"] == (r.lo_5to10 <= 1 <= r.hi_5to10) for m, r in zip(sq2, B.itertuples())) and all(m["open"] == (r.lo_gt10 <= 1 <= r.hi_gt10) for m, r in zip(tr2, B.itertuples()))
ck.log(f"markers read back: a 11 circles on the V13 axis (residual under 0.05 pt); b {NB} squares, {NB} triangles, {NB} reference circles on the V13 axis (residual under 0.05 pt), open exactly where the interval includes 1",
       len(datA2) == 11 and float(np.abs(rA2).max()) < 0.05 and abs(aA2 - aA) < 0.05 and len(sq2) == len(tr2) == len(ref2) == NB and float(np.abs(rB2).max()) < 0.05 and abs(aB2 - aB) < 0.05 and fill_ok and all(abs(m["cx"] - X_REF_B) < 0.05 for m in ref2),
       f"a {len(datA2)} res {np.abs(rA2).max():.3f}; b {len(sq2)}/{len(tr2)}/{len(ref2)} res {np.abs(rB2).max():.3f}; filled squares {sum(not m['open'] for m in sq2)} triangles {sum(not m['open'] for m in tr2)}")
legM = [m for m in M2 if LEG_X_LINE0 - 1 <= m["cx"] <= LEG_X_LINE1 + 1 and LEG_BOX[1] <= m["cy"] <= LEG_BOX[3]]
kinds = sorted((m["shape"], m["open"], m["colour"]) for m in legM)
ck.log("legend handles read back: one filled square (>5 to 10% colour), one filled triangle (>10% colour), one open grey circle, one open square and one open triangle, the marker sizes equal to the plotted ones",
       kinds == sorted([("square", False, MID), ("triangle", False, BLUE), ("circle", True, GREY), ("square", True, MID), ("triangle", True, BLUE)]) and all(abs(m["w"] - {"square": SQ_SIDE, "triangle": TR_W, "circle": REF_D}[m["shape"]]) < 0.1 for m in legM), f"{kinds}")
legL = [h for h in hlines(items2, min_len=10) if LEG_X_LINE0 - 0.5 <= h["x0"] and h["x1"] <= LEG_X_LINE1 + 0.5 and LEG_BOX[1] <= h["y"] <= LEG_BOX[3]]
ck.log("legend lines read back: two lines of the CI width, one in each series colour, each through its marker centre", sorted(h["stroke"] for h in legL) == sorted([MID, BLUE]) and all(abs(h["width"] - CI_W) < 0.05 for h in legL) and all(any(abs(h["y"] - hd["cy"]) < 0.05 for hd in leg_handles if hd["line"]) for h in legL), f"{[(h['stroke'], h['width'], h['y']) for h in legL]}")
ck.log("no grey #ccd1d6 connector and no pale grid line on the sheet (V13 design)", not [h for h in hlines(items2, min_len=0.0) if h["stroke"] in ("#ccd1d6", "#eef0f1")] and not [it for it in items2 if hexcol(it["color"]) == "#eef0f1"], "")
cols = {hexcol(it["fill"]) for it in items2 if it["fill"]} | {hexcol(it["color"]) for it in items2 if it["color"]}
ck.log("colours within the V13 set", cols <= ({hexcol(it["fill"]) for it in items if it["fill"]} | {hexcol(it["color"]) for it in items if it["color"]}), f"{sorted(cols)}")
W0, N0 = v14_multisets(v14_spans(BASE)); W1, N1 = v14_multisets(v14_spans(OUT))
tok = lambda strs: Counter(x for t in strs for x in re.findall(r"-?\d+\.?\d*", t)); wtok = lambda strs: Counter(w for t in strs for w in t.split())
old_data = old_labA + old_hrA + old_qA + old_labB + old_hrB + old_pB + [s["text"] for s in key_txt]; new_data = [plab(o) for o in A.outcome] + [hr_ci(r.hr_sd, r.lo_sd, r.hi_sd) for r in A.itertuples()] + [fmt_p(r.q_sd) for r in A.itertuples()] + [plab(o) for o in B.outcome] + expB_hr + expB_p + [l_["text"] for l_ in leg_lines]
ck.log("numeric-token multiset differs from V13 exactly by the data strings and the key -> legend strings; static tokens identical", (N0 - N1) == (tok(old_data) - tok(new_data)) and (N1 - N0) == (tok(new_data) - tok(old_data)), f"removed {sum((N0 - N1).values())} added {sum((N1 - N0).values())}")
ck.log("word multiset differs from V13 exactly by the data strings and the key -> legend strings", (W0 - W1) == (wtok(old_data) - wtok(new_data)) and (W1 - W0) == (wtok(new_data) - wtok(old_data)), f"removed {dict(W0 - W1)} added {dict(W1 - W0)}")
n_rows, n_no = printed_values_csv(OUT, drawn, f"{WB}/{SHEET}_printed_values_rebuilt.csv", static_from=BASE)
ck.log(f"printed values csv (the rebuilt sheet before the title strip): {n_rows} rows, {n_no} NO", n_no == 0, f"{WB}/{SHEET}_printed_values_rebuilt.csv")
dump({"panel_a_axis": {"a": aA, "b": bA}, "panel_b_axis": {"a": aB, "b": bB, "x_ref": X_REF_B}, "delta_pt": DELTA, "page_height": [PH, PH2], "a_rows": rowsA, "b_rows": rowsB, "provenance": {"continuous": sc_c, "graded": sc_g}}, f"{WB}/{SHEET}_readback.json")
ck.log(f"round 49 sizes: labels and tick labels {10.0 + PT_PLUS:g} pt except {SIZE_KEEP} (labels: the largest size at or above the V26 size, 0.5 pt steps, 4 pt clear of the row's marks and the reference rule; tick labels: a pair that would come within 2.5 pt of each other keeps its V26 size), HR and q/P strings {9.5 + PT_PLUS:g}, heads {9.5 + PT_PLUS:g} bold, titles {titleA['size'] + PT_PLUS:g}, letters {13.0 + PT_PLUS:g}, legend {LEG_SIZE + PT_PLUS:g}", True, "")
gapsB = [((t1["x"] - S.width(l1["text"], TSB[l1["text"]]) / 2) - (t0["x"] + S.width(l0["text"], TSB[l0["text"]]) / 2)) for (t0, l0), (t1, l1) in zip(zip(ticksB, tlB), zip(ticksB[1:], tlB[1:]))]
ck.log(f"panel b tick labels: adjacent labels at least 2.5 pt apart at their drawn sizes (gaps {[round(g, 2) for g in gapsB]})", min(gapsB) >= 2.5, "")
dump({"pt_plus": PT_PLUS, "size_exceptions": SIZE_KEEP, "nudges": {"q_p_column_x": [Q_X, Q_X_NEW], "legend_line_x": [[366.0, 386.0], [LEG_X_LINE0, LEG_X_LINE1]], "legend_text_x": [392.0, LEG_X_TEXT], "title_b_centre_x": [cx13(titleB), TITLE_B_CX]}, "not_drawn": [], "legend_rewrap": {"removed": ["Residual T90 at or below 5%,"], "added": ["Residual T90", "at or below 5%,"]}}, f"{W}/r49_record.json")
ok = ck.write(f"{WB}/build_checks.txt"); print("RESULT ALL PASS" if ok else "RESULT FAIL")
