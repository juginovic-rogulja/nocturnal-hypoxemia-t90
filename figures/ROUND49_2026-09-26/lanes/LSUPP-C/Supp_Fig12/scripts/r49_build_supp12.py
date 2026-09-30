#!$T90_PY
"""Supp_Fig12, round 49 (lane LSUPP-C): the V26 sheet (round 40) rebuilt from its data builder with every text one point larger.
gen_supp12.py (lane copy: paths, PT_PLUS on LABF/TCKF/ANNF and the 9.2 pt key labels, geometry constants untouched) is run plain;
the 171 mm output is cropped to the V26 page box (425.64 x 468.72) by PyMuPDF in a child under the watchdog (the round-33/40
rule); verified against V26 with r49_verify."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import json, os, subprocess, sys
from collections import Counter
L = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LSUPP-C"
sys.path.insert(0, f"{L}/scripts")
import r49_common as C
import r49_verify as V
S = "Supp_Fig12"; D = f"{L}/{S}"; W = f"{D}/work"; SC = f"{L}/L5b_scripts"; WK = f"{L}/L5b_work"; os.makedirs(W, exist_ok=True)
OLD = f"{C.V26}/{S}.pdf"; NEW = f"{D}/{S}.pdf"
SV = f"{paths.SV_ROOT}/pap_residual_continuous"
SRC = {f: C.sha256(f"{SV}/{f}") for f in ("results_continuous.csv", "delta_results.csv")}
for f in SRC: C.hydrated(f"{SV}/{f}.provenance.json")
R40 = json.load(open(f"{paths.FIGURE_ROOT}/ROUND40_2026-09-18/lanes/LSUPP/Supp_Fig12/work/build_record.json"))["sources_sha256"]
if "--no-gen" not in sys.argv:
    r = subprocess.run([C.PY, f"{SC}/gen_supp12.py"], cwd=WK, capture_output=True, text=True)
    open(f"{W}/gen_supp12.log", "w").write(r.stdout + "\n--- stderr ---\n" + r.stderr); assert r.returncode == 0, r.stderr[-3000:]
GEN = f"{WK}/polish_out/pdf/eFigure21_pap_correction_vs_residual.pdf"; pw, ph = C.page_box(GEN); assert abs(pw - 484.72) < 0.2, pw
w26, h26 = C.page_box(OLD); assert abs(w26 - 425.64) < 0.05 and abs(h26 - ph) < 0.5, (w26, h26, ph)
job = {"mode": "box", "src": GEN, "out": NEW, "box": [0, 0, w26, ph], "metadata": {"title": S, "creator": "LSUPP-C r49 gen_supp12.py (matplotlib, +1 pt) + fitz_place box"}}
json.dump(job, open(f"{W}/place_job.json", "w")); C.wd(f"{S}_place", [C.PY, f"{L}/scripts/fitz_place.py", f"{W}/place_job.json"])
DV = json.load(open(f"{WK}/polish_out/work/eFigure21_polish_drawn_values.json")); assert len(DV) == 11, len(DV)
drawn = Counter()
for d in DV:
    for k in ("outcome", "printed_residual", "printed_residual_q", "printed_delta", "printed_delta_q"): drawn[d[k]] += 1
OLD_LABEL = "Amount corrected by PAP, per 1 SD, residual held fixed"; NEW_LABEL = "Amount corrected by PAP, per 1 SD, adjusted for the residual T90"
gen_log = open(f"{W}/gen_supp12.log").read()
def extra(ck, ctx):
    ck.log("numbers files unchanged since the round-40 build (sha256 of results_continuous.csv and delta_results.csv equal the round-40 build record)", SRC == R40, str({k: v[:12] for k, v in SRC.items()}))
    ck.log("the builder's own data gates passed (rc 0: every drawn value asserted against its CSV, FDR rows from the file, q from the file)", "every value asserted against its CSV" in gen_log and "rows drawn: 11" in gen_log, "gen_supp12.log")
    sp = ctx["sp_new"]; key = [s for s in sp if s["text"].startswith("Amount corrected")]; assert key, "key label 2 not found"
    row = [s for s in sp if abs(s["y"] - key[0]["y"]) <= 1.5 and s["x"] >= 100]; x1 = max(s["x1"] for s in row)
    ck.log("the reworded key entry at 10.2 pt still ends before the V26 crop edge (425.64 - 5 pt, the round-40 gate)", x1 <= w26 - 5.0, f"ends at x {x1:.1f} pt of {w26}")
    q = [s for s in sp if s["text"] in ("<0.001", "0.0039", "0.0018") or (s["x"] > 370 and s["size"] > 10)]; x1q = max(s["x1"] for s in q) if q else 0
    ck.log("the q column at 10.5 pt ends before the crop edge", x1q <= w26 - 2.0, f"ends at x {x1q:.1f} pt of {w26}")
spec = dict(sheet=S, old_pdf=OLD, new_pdf=NEW, out_dir=D, kept=[], note="every text +1 pt (11/10/9.5/9.2 -> 12/11/10.5/10.2), page cropped to the V26 box, key at the V26 positions, marks unchanged",
            mark_colours=["#0288d1", "#8a9099"], legend_boxes=[], drawn=drawn, present=[NEW_LABEL, "Residual sleep T90 on PAP, per 1 SD", "95% CI includes 1"], absent=[OLD_LABEL], sources=SRC, extra=extra,
            notes=[f"builder: {SC}/gen_supp12.py (edits in scripts/R49_EDITS.md), generated {GEN}", "round-40 steps re-applied: the key wording of round 40 is in the generator copy, the round-33 crop to 425.64 pt by fitz_place box mode"])
ok = V.run(spec)
json.dump({"sources_sha256": SRC, "generator": f"{SC}/gen_supp12.py", "generated": GEN, "new": NEW, "new_sha256": C.sha256(NEW), "written": C.now()}, open(f"{W}/build_record.json", "w"), indent=1)
sys.exit(0 if ok else 1)
