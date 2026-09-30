#!$T90_PY
"""Supp_Fig13, round 50 (lane LSUPP-C, Alen's item 3): the lower row tidied. gen_supp13.py (this lane's copy of the round-49 builder,
re-laid: part 1 = panel a, part 2 = panel c, part 3 = the lower row b, d, e on one baseline with one height, three axes of one width,
even gutters, the odds-ratio note as b's one-line subtitle) builds the three parts; they are composed by PyMuPDF in a child under the
watchdog (a at (14, 14), c at (522.72, 14) as before, the row across the page below them, the page taller); verified against V27 with
r50_verify: same strings at the same sizes except the note's two lines becoming one, the regions of a and c pixel-identical, the row's
values equal to the round-49 record, letters and baselines proved from the text layer."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import json, os, subprocess, sys
from collections import Counter
from decimal import Decimal, ROUND_HALF_UP
L = f"{paths.FIGURE_ROOT}/ROUND50_2026-09-26/lanes/LSUPP-C"
sys.path.insert(0, f"{L}/scripts")
import r49_common as C
import r50_verify as V
S = "Supp_Fig13"; D = f"{L}/{S}"; W = f"{D}/work"; SC = f"{L}/L5b_scripts"; WK = f"{L}/L5b_work"; os.makedirs(W, exist_ok=True)
OLD = f"{C.V27}/{S}.pdf"; NEW = f"{D}/{S}.pdf"
NUM = paths.NUMBERS_DIR
SRC = {f: C.sha256(f"{NUM}/{f}") for f in ("nonresponder_phenotype_v2.csv", "nonresponder_combined_v2.json", "treatment_v2.json", "negcontrols_final.json")}
for f in SRC: C.hydrated(f"{NUM}/{f}.provenance.json")
R49L = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LSUPP-C"
R49 = json.load(open(f"{R49L}/Supp_Fig13/work/build_record.json"))["sources_sha256"]
FP = f"{WK}/work/fullprec_nonresponder_v8_1.json"; FPP = json.load(open(FP + ".provenance.json"))
assert FPP["output_sha256"] == C.sha256(FP), "twin sidecar mismatch"
assert FPP["reconciled_against"]["sha256"] == SRC["nonresponder_combined_v2.json"], "twin reconciled against another combined json"
FULL = json.load(open(FP)); P = json.load(open(f"{NUM}/nonresponder_combined_v2.json"))
if "--no-gen" not in sys.argv:
    r = subprocess.run([C.PY, f"{SC}/gen_supp13.py"], cwd=WK, capture_output=True, text=True)
    open(f"{W}/gen_supp13.log", "w").write(r.stdout + "\n--- stderr ---\n" + r.stderr); assert r.returncode == 0, r.stderr[-3000:]
P1 = f"{WK}/polish_out/pdf/eFigure10a_nonresponder_phenotype.pdf"; P2 = f"{WK}/polish_out/pdf/eFigure10b_nonresponder_phenotype.pdf"; P3 = f"{WK}/polish_out/pdf/eFigure10c_nonresponder_phenotype_row.pdf"
LAY = json.load(open(f"{WK}/polish_out/work/eFigure10_row_layout.json"))
OFF1, OFF2, MARGIN = (14.0, 14.0), (522.72, 14.0), 14.0
w1, h1 = C.page_box(P1); w2, h2 = C.page_box(P2); w3, h3 = C.page_box(P3)
assert abs(w1 - 484.72) < 0.2 and abs(w2 - 484.72) < 0.2, (w1, w2)
PW = OFF2[0] + w2 + MARGIN; assert abs(w3 - (PW - 2 * MARGIN)) < 0.2, (w3, PW)
Y_ROW = OFF1[1] + LAY["row_top_from_part1_top_in"] * 72.0; PH = Y_ROW + h3 + MARGIN
w27, h27 = C.page_box(OLD); assert abs(PW - w27) < 0.05, (PW, w27)
job = {"mode": "compose", "out": NEW, "page": [PW, PH], "parts": [{"pdf": P1, "rect": [OFF1[0], OFF1[1], OFF1[0] + w1, OFF1[1] + h1]}, {"pdf": P2, "rect": [OFF2[0], OFF2[1], OFF2[0] + w2, OFF2[1] + h2]}, {"pdf": P3, "rect": [MARGIN, Y_ROW, MARGIN + w3, Y_ROW + h3]}], "letters": [],
       "metadata": {"title": S, "creator": "LSUPP-C r50 gen_supp13.py (matplotlib) + fitz_place compose"}}
json.dump(job, open(f"{W}/place_job.json", "w")); C.wd(f"{S}_place", [C.PY, f"{L}/scripts/fitz_place.py", f"{W}/place_job.json"])
def _r(v, nd=2): return str(Decimal(repr(float(v))).quantize(Decimal(1).scaleb(-nd), rounding=ROUND_HALF_UP))
RATES = {k: v for k, v in P["failure_rate_by_condition"].items() if v["fdr_significant"]}; ORDER = sorted(RATES, key=lambda k: RATES[k]["fail_pct_with"])
drawn = Counter()
for k in ORDER:
    fp = FULL["by_condition"][k]; drawn[k] += 1; drawn[f"{_r(fp['marginal_or'])} ({_r(fp['marginal_lo'])}-{_r(fp['marginal_hi'])}){'*' if RATES[k]['independent'] else ''}"] += 1
drawn[f"all patients, {P['failure_pct']:.1f}%"] += 1
tr = FULL["dose_response_trend"]; TS = f"{_r(tr['or'])} (95% CI, {_r(tr['lo'])}-{_r(tr['hi'])})"; NOTE = f"odds ratio per added condition {TS}"
drawn[NOTE] += 1
KEEP = sorted(P["independent_conditions"], key=lambda k: -P["independent_conditions"][k]["or"])
FPC = {**FULL["joint_independent"], "Age, per year": FULL["joint_demographics"]["AgeAtVisit"], "Male sex": FULL["joint_demographics"]["male"], "Black race": FULL["joint_demographics"]["race_black"]}
for k in KEEP + ["Age, per year", "Male sex", "Black race"]:
    drawn[k] += 1; drawn[f"{_r(FPC[k]['or'])} ({_r(FPC[k]['lo'])}-{_r(FPC[k]['hi'])})"] += 1
DV3 = json.load(open(f"{WK}/polish_out/work/eFigure10c_nonresponder_phenotype_row_polished_drawn_values.json"))
DV49 = json.load(open(f"{R49L}/L5b_work/polish_out/work/eFigure10a_nonresponder_phenotype_polished_drawn_values.json"))
gen_log = open(f"{W}/gen_supp13.log").read()
REG_A = [0.0, 0.0, 512.0, 383.0]; REG_C = [516.0, 0.0, PW, 205.0]
def extra(ck, ctx):
    ck.log("numbers files unchanged since the round-49 build (sha256 equal the round-49 build record)", SRC == R49, str({k: v[:12] for k, v in SRC.items()}))
    ck.log("the builder's own data gates passed (rc 0: printed-number audit of every value against its source, row sets from the files, twin within 5e-4 of the json)", "printed-number audit: 91 values reproduced" in gen_log, "gen_supp13.log")
    ck.log("the lower row's values equal the round-49 record (panel b bins, the trend, panel d AUCs, panel e fifths)", DV3["panelB"] == DV49["panelB"] and DV3["panelB_trend"] == DV49["panelB_trend"] and DV3["panelD"] == DV49["panelD"] and DV3["panelE"] == DV49["panelE"], "")
    sp, spo = ctx["sp_new"], ctx["sp_old"]
    def letters(s_): return {t["text"]: (round(t["x"], 1), round(t["y"], 1), t["size"]) for t in s_ if len(t["text"]) == 1 and t["text"] in "abcde" and t["size"] >= 12}
    Lo, Ln = letters(spo), letters(sp)
    ck.log("letters a and c at the V27 origins (0.5 pt), 14 pt", all(abs(Lo[k][0] - Ln[k][0]) <= 0.5 and abs(Lo[k][1] - Ln[k][1]) <= 0.5 and Ln[k][2] == 14.0 for k in "ac"), f"V27 {Lo} NEW {Ln}")
    ck.log("letters b, d and e on one line at 14 pt, at the panels' top-left corners (x = the row's edge, d's label gutter, e's title gutter)", all(k in Ln for k in "bde") and abs(Ln["b"][1] - Ln["d"][1]) <= 0.5 and abs(Ln["d"][1] - Ln["e"][1]) <= 0.5 and all(Ln[k][2] == 14.0 for k in "bde")
           and abs(Ln["b"][0] - (MARGIN + LAY["letter_x_in"]["b"] * 72)) <= 1.0 and abs(Ln["d"][0] - (MARGIN + LAY["letter_x_in"]["d"] * 72)) <= 1.0 and abs(Ln["e"][0] - (MARGIN + LAY["letter_x_in"]["e"] * 72)) <= 1.0, f"NEW {Ln}, layout x {LAY['letter_x_in']}")
    xt = [t for t in sp if t["y"] > Y_ROW and t["text"] in ("3 or more", "0.5", "0.7", "0.9", "1, lowest", "5, highest")]   # x tick labels only (b's and e's y tick "0" sit on the axes bottom line itself)
    ys = sorted({round(t["y"]) for t in xt})
    ck.log("one baseline: the x tick labels of b, d and e share one baseline (1 pt), so the three axes share their bottom", len(xt) >= 6 and max(ys) - min(ys) <= 1.0, f"tick baselines {ys} ({len(xt)} labels: {[t['text'] for t in xt]})")
    ax = LAY["axes_x0_in"]; wax = LAY["axes_width_in"]; g = LAY["ROW_GAP_in"]
    ck.log("even gutters: the gap from b's axes to d's label block and from d's axes to e's title block both equal ROW_GAP, one axes width and one height for b, d, e (the builder's layout record)",
           abs((LAY["letter_x_in"]["d"] - (ax["b"] + wax)) - g) < 1e-6 and abs((LAY["letter_x_in"]["e"] - (ax["d"] + wax)) - g) < 1e-6, f"axes width {wax:.3f} in, gutter {g:.2f} in, axes x0 {ax}")
    ck.log("the odds-ratio note is one line directly above b's bars (left-aligned with b's axes, one string)", any(t["text"].startswith("odds ratio per added condition") and abs(t["x"] - (MARGIN + ax["b"] * 72)) <= 1.5 and t["y"] > Y_ROW for t in sp), str([(round(t["x"], 1), round(t["y"], 1), t["text"]) for t in sp if t["text"].startswith("odds ratio")]))
    ck.info(f"parts: a {w1:.2f} x {h1:.2f} at {OFF1}, c {w2:.2f} x {h2:.2f} at {OFF2}, row {w3:.2f} x {h3:.2f} at (14, {Y_ROW:.2f}); page {PW:.3f} x {PH:.3f} (V27 {w27:.3f} x {h27:.3f})")
spec = dict(sheet=S, old_pdf=OLD, new_pdf=NEW, out_dir=D, page_box=[PW, PH], page_box_note=f"declared: the tidier lower row needs a taller page, V27 {h27:.2f} -> {PH:.2f} pt (width unchanged)",
            removed=[("odds ratio per added condition", 10.5), (TS, 10.5)], added=[(NOTE, 10.5)],
            note="item 3: b, d, e on one baseline with one height and even gutters, the note as b's one-line subtitle, a and c unchanged",
            regions=[dict(name="panel a (part 1)", box=REG_A, allowed=[], note="a unchanged, same page coordinates"), dict(name="panel c (part 2)", box=REG_C, allowed=[], note="c unchanged, same page coordinates")],
            drawn=drawn, present=["Odds ratio (95% CI)", "Condition present before PAP", "Condition absent", "Cardiopulmonary conditions before PAP", "Fifth of predicted risk (quintile)"], absent=["remains independent"], sources=SRC, extra=extra,
            notes=[f"builder: {SC}/gen_supp13.py (edits in scripts/R50_EDITS.md), parts {P1}, {P2}, {P3}", f"layout record: {json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in LAY.items() if k != 'rule'})}"])
ok = V.run(spec)
json.dump({"sources_sha256": SRC, "twin": {"path": FP, "sha256": C.sha256(FP)}, "generator": f"{SC}/gen_supp13.py", "parts": [P1, P2, P3], "layout": LAY, "page": [PW, PH], "row_top": Y_ROW, "new": NEW, "new_sha256": C.sha256(NEW), "written": C.now()}, open(f"{W}/build_record.json", "w"), indent=1)
sys.exit(0 if ok else 1)
