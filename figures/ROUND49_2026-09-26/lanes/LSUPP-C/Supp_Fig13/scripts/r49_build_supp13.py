#!$T90_PY
"""Supp_Fig13, round 49 (lane LSUPP-C): the V26 sheet (round-37 v8.2 build, round-40 LNOTES note removal) rebuilt with every text one
point larger. gen_supp13.py (lane copy: the twin regenerated in this lane, PT_PLUS through supp_polish_common, the asterisk key not
printed, panel d's four wrapped model labels kept at 10 pt, the panel a gutter rule relaxed by 2.2 pt where the edge carries no ink,
the margin gate at 1.5 mm for the column header) builds the two 171 mm parts; they are composed as the round-30/37 sheet holds them
(part 1 at (14, 14), part 2 at (522.72, 14), page = the taller part + 28) by PyMuPDF in a child under the watchdog; verified against
V26 with r49_verify."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import json, os, subprocess, sys
from collections import Counter
from decimal import Decimal, ROUND_HALF_UP
L = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LSUPP-C"
sys.path.insert(0, f"{L}/scripts")
import r49_common as C
import r49_verify as V
S = "Supp_Fig13"; D = f"{L}/{S}"; W = f"{D}/work"; SC = f"{L}/L5b_scripts"; WK = f"{L}/L5b_work"; os.makedirs(W, exist_ok=True)
OLD = f"{C.V26}/{S}.pdf"; NEW = f"{D}/{S}.pdf"
NUM = paths.NUMBERS_DIR
SRC = {f: C.sha256(f"{NUM}/{f}") for f in ("nonresponder_phenotype_v2.csv", "nonresponder_combined_v2.json", "treatment_v2.json", "negcontrols_final.json")}
for f in SRC: C.hydrated(f"{NUM}/{f}.provenance.json")
R40 = json.load(open(f"{paths.FIGURE_ROOT}/ROUND40_2026-09-18/lanes/LSUPP/Supp_Fig14/work/build_record.json"))["sources_sha256"]
FP = f"{WK}/work/fullprec_nonresponder_v8_1.json"; FPP = json.load(open(FP + ".provenance.json"))
assert FPP["output_sha256"] == C.sha256(FP), "twin sidecar mismatch"
assert FPP["reconciled_against"]["sha256"] == SRC["nonresponder_combined_v2.json"], "twin reconciled against another combined json"
FULL = json.load(open(FP)); P = json.load(open(f"{NUM}/nonresponder_combined_v2.json"))
if "--no-gen" not in sys.argv:
    r = subprocess.run([C.PY, f"{SC}/gen_supp13.py"], cwd=WK, capture_output=True, text=True)
    open(f"{W}/gen_supp13.log", "w").write(r.stdout + "\n--- stderr ---\n" + r.stderr); assert r.returncode == 0, r.stderr[-3000:]
P1 = f"{WK}/polish_out/pdf/eFigure10a_nonresponder_phenotype.pdf"; P2 = f"{WK}/polish_out/pdf/eFigure10b_nonresponder_phenotype.pdf"
OFF1, OFF2, MARGIN = (14.0, 14.0), (522.72, 14.0), 14.0
w1, h1 = C.page_box(P1); w2, h2 = C.page_box(P2); assert abs(w1 - 484.72) < 0.2 and abs(w2 - 484.72) < 0.2, (w1, w2)
PW = OFF2[0] + w2 + MARGIN; PH = max(h1, h2) + 2 * MARGIN
w26, h26 = C.page_box(OLD); assert abs(PW - w26) < 0.05 and abs(PH - h26) < 0.05, (PW, PH, w26, h26)
job = {"mode": "compose", "out": NEW, "page": [PW, PH], "parts": [{"pdf": P1, "rect": [OFF1[0], OFF1[1], OFF1[0] + w1, OFF1[1] + h1]}, {"pdf": P2, "rect": [OFF2[0], OFF2[1], OFF2[0] + w2, OFF2[1] + h2]}], "letters": [],
       "metadata": {"title": S, "creator": "LSUPP-C r49 gen_supp13.py (matplotlib, +1 pt) + fitz_place compose"}}
json.dump(job, open(f"{W}/place_job.json", "w")); C.wd(f"{S}_place", [C.PY, f"{L}/scripts/fitz_place.py", f"{W}/place_job.json"])
# the drawn strings (the builder's printing rules on the twin, the strings V26 prints)
def _r(v, nd=2): return str(Decimal(repr(float(v))).quantize(Decimal(1).scaleb(-nd), rounding=ROUND_HALF_UP))
RATES = {k: v for k, v in P["failure_rate_by_condition"].items() if v["fdr_significant"]}; ORDER = sorted(RATES, key=lambda k: RATES[k]["fail_pct_with"])
drawn = Counter()
for k in ORDER:
    fp = FULL["by_condition"][k]; drawn[k] += 1; drawn[f"{_r(fp['marginal_or'])} ({_r(fp['marginal_lo'])}-{_r(fp['marginal_hi'])}){'*' if RATES[k]['independent'] else ''}"] += 1
drawn[f"all patients, {P['failure_pct']:.1f}%"] += 1
tr = FULL["dose_response_trend"]; drawn[f"{_r(tr['or'])} (95% CI, {_r(tr['lo'])}-{_r(tr['hi'])})"] += 1
KEEP = sorted(P["independent_conditions"], key=lambda k: -P["independent_conditions"][k]["or"])
FPC = {**FULL["joint_independent"], "Age, per year": FULL["joint_demographics"]["AgeAtVisit"], "Male sex": FULL["joint_demographics"]["male"], "Black race": FULL["joint_demographics"]["race_black"]}
for k in KEEP + ["Age, per year", "Male sex", "Black race"]:
    drawn[k] += 1; drawn[f"{_r(FPC[k]['or'])} ({_r(FPC[k]['lo'])}-{_r(FPC[k]['hi'])})"] += 1
KEPT = [(t, 10.0, 1) for t in ("Pretreatment oxygen", "alone", "and demographics", "and the", "cardiopulmonary count", "and the independent", "conditions")]
NOTE = "* remains independent in the joint model of panel c"
gen_log = open(f"{W}/gen_supp13.log").read()
def extra(ck, ctx):
    ck.log("numbers files unchanged since the round-40 build (sha256 of the four numbers files equal the round-40 build record)", SRC == R40, str({k: v[:12] for k, v in SRC.items()}))
    ck.log("full-precision twin regenerated in this lane and reconciled against the current v8.3 json (sidecar sha equal, every value within 5e-4, 0 printed strings move)", FPP["reconciled_against"]["sha256"] == SRC["nonresponder_combined_v2.json"] and "0 of 21 move" in open(f"{WK}/fullprec_r49_run.log").read(), f"{FP} sha {C.sha256(FP, 12)}")
    ck.log("the builder's own data gates passed (rc 0: printed-number audit of every value against its source, row sets from the files, twin within 5e-4 of the json)", "printed-number audit: 91 values reproduced" in gen_log, "gen_supp13.log")
    ck.log("round-40 LNOTES removal kept: the asterisk key sentence is not on the sheet", not any("remains independent" in ln for ln in ctx["lines_new"]), "")
    def letters(sp): return sorted([(s["text"], round(s["x"], 2), round(s["y"], 2), s["size"]) for s in sp if len(s["text"]) == 1 and s["text"] in "abcde" and s["size"] >= 12], key=lambda t: t[0])
    Lo, Ln = letters(ctx["sp_old"]), letters(ctx["sp_new"])
    ck.log("the five panel letters are on NEW at 14 pt at the V26 x origins (0.5 pt), baselines within 2 pt (the letters are top-anchored and grow downward)", [t[0] for t in Lo] == [t[0] for t in Ln] == ["a", "b", "c", "d", "e"] and all(abs(a[1] - b[1]) <= 0.5 and abs(a[2] - b[2]) <= 2.0 and b[3] == 14.0 for a, b in zip(Lo, Ln)), f"V26 {Lo} NEW {Ln}")
    comps = C.diff_components(ctx["R"]["old_noaa"], ctx["R"]["new_noaa"], spec["mark_colours"], 150, spec["legend_boxes"])
    ck.log("every differing mark pixel outside the legend box is rasterization flicker (a component of at most 2 px at an interval end or marker edge), not a moved mark: V26 drew the panel a and c intervals at the 15 Sep json's 3 dp values, this build at the 16 Sep v8.3 full-precision values (shifts under 0.05 px)", all(c[4] <= 2 for c in comps), f"{len(comps)} components, largest {max([c[4] for c in comps] or [0])} px, at {[c[:4] for c in comps][:12]}")
    ck.info(f"parts: 1 = {w1:.2f} x {h1:.2f} pt at {OFF1}, 2 = {w2:.2f} x {h2:.2f} pt at {OFF2}, page {PW:.3f} x {PH:.3f}")
sp_old_tmp = None
spec = dict(sheet=S, old_pdf=OLD, new_pdf=NEW, out_dir=D, kept=KEPT,
            note="every text +1 pt (13/11/10/9.5 -> 14/12/11/10.5) except panel d's four wrapped model labels (kept at 10 pt: at 11 pt the longest would end 1.1 pt before the bars that start at that edge); panel a's two longest row labels end 2.0 and 1.3 pt past an axes edge that carries no ink (no spine, no band, nearest mark 1.04 in); the column header 'Odds ratio (95% CI)' at 10.5 pt ends 1.6 mm from the part's right edge (6.5 mm inside the composed page)",
            mark_colours=["#0288d1", "#3f9fd8", "#ccd1d6", "#e2e7ea", "#b6c1c8", "#8a9aa5", "#4d6577"], marks_tolerance=20, drawn=drawn, present=["Odds ratio (95% CI)", "Condition present before PAP", "Condition absent"], absent=[NOTE, "remains independent"], sources=SRC, extra=extra,
            notes=[f"builder: {SC}/gen_supp13.py (edits in scripts/R49_EDITS.md), parts {P1}, {P2}", "the grey #8a9099 is not in the marks check (panel c's demographic labels are grey text); the blue series and every ramp step are"])
# legend box of panel a (its handles move with the 11 pt legend): from both text layers, computed inside run through a closure
def legend_boxes_from(sp_old, sp_new):
    b = [V.box_around(sp, lambda t: t.startswith("Condition present") or t.startswith("Condition absent") or t in ("PAP", "absent"), pad_left=40, pad_right=6, pad_y=8) for sp in (sp_old, sp_new)]
    b = [x for x in b if x]; return [[min(x[0] for x in b), min(x[1] for x in b), max(x[2] for x in b), max(x[3] for x in b)]] if b else []
_run = V.run
def run_with_boxes(spec):
    sp_old = C.gs_text(OLD, f"{W}/{S}_V26_txtwrite_pre.xml", f"{S}_txt_old_pre"); sp_new = C.gs_text(NEW, f"{W}/{S}_NEW_txtwrite_pre.xml", f"{S}_txt_new_pre")
    spec["legend_boxes"] = legend_boxes_from(sp_old, sp_new); print("legend boxes", spec["legend_boxes"]); return _run(spec)
ok = run_with_boxes(spec)
json.dump({"sources_sha256": SRC, "twin": {"path": FP, "sha256": C.sha256(FP)}, "generator": f"{SC}/gen_supp13.py", "parts": [P1, P2], "page": [PW, PH], "new": NEW, "new_sha256": C.sha256(NEW), "written": C.now()}, open(f"{W}/build_record.json", "w"), indent=1)
sys.exit(0 if ok else 1)
