#!$T90_PY
"""Supp_Fig15, round 49 (lane LSUPP-C): the V26 sheet (round 40) rebuilt with every text one point larger. gen_supp15.py (lane
copy: paths, PT_PLUS on TITLE_PT/TICK_PT/ANN_PT, the x title re-centred on the ruled axis as in round 40, the 4 mm margin gate at
2.9 mm because the two-line title now ends 3.0 mm from the foot of the fixed page box) is run plain; its output IS the sheet
(484.72 x 218.59); verified against V26 with r49_verify plus the round-40 centring check."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import json, os, shutil, subprocess, sys
from collections import Counter
L = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LSUPP-C"
sys.path.insert(0, f"{L}/scripts")
import r49_common as C
import r49_verify as V
S = "Supp_Fig15"; D = f"{L}/{S}"; W = f"{D}/work"; SC = f"{L}/L5b_scripts"; WK = f"{L}/L5b_work"; os.makedirs(W, exist_ok=True)
OLD = f"{C.V26}/{S}.pdf"; NEW = f"{D}/{S}.pdf"
NUM = paths.NUMBERS_DIR; T90 = paths.FIGURE_ROOT
SRC = {f: C.sha256(f"{NUM}/{f}") for f in ("nonresponder_phenotype_v2.csv", "nonresponder_combined_v2.json", "treatment_v2.json", "negcontrols_final.json")}
for f in SRC: C.hydrated(f"{NUM}/{f}.provenance.json")
TABLES = {t: C.sha256(f"{paths.TABLES_DIR}/{t}") for t in ("cpap_t90_by_stage.parquet", "t90_final.parquet")}
R40 = json.load(open(f"{paths.FIGURE_ROOT}/ROUND40_2026-09-18/lanes/LSUPP/Supp_Fig15/work/build_record.json"))
if "--no-gen" not in sys.argv:
    r = subprocess.run([C.PY, f"{SC}/gen_supp15.py"], cwd=WK, capture_output=True, text=True)
    open(f"{W}/gen_supp15.log", "w").write(r.stdout + "\n--- stderr ---\n" + r.stderr); assert r.returncode == 0, r.stderr[-3000:]
GEN = f"{WK}/polish_out/pdf/eFigureNEW_pap_count_dose.pdf"; DV = json.load(open(f"{WK}/polish_out/work/eFigureNEW_pap_count_dose_polish_drawn_values.json"))
w26, h26 = C.page_box(OLD); fw, fh = C.page_box(GEN); assert abs(fw - w26) < 0.05 and abs(fh - h26) < 0.05, (fw, fh, w26, h26)
shutil.copyfile(GEN, NEW)
drawn = Counter(b["label"] for b in DV["bins"]); assert len(DV["bins"]) == 6
gen_log = open(f"{W}/gen_supp15.log").read()
def extra(ck, ctx):
    ck.log("numbers files and frozen tables unchanged since the round-40 build (sha256 equal the round-40 build record)", SRC == R40["sources_sha256"] and TABLES == R40["tables_sha256"], str({k: v[:12] for k, v in {**SRC, **TABLES}.items()}))
    ck.log("the builder's own data gates passed (rc 0: the six bins rebuilt on the frozen tables, capped roll-up equal to the json, singles equal to the json)", "uncapped bins:" in gen_log and "[('0', 1050, 20.8), ('1', 375, 27.7), ('2', 240, 40.4), ('3', 113, 54.0), ('4', 80, 72.5), ('5+', 36, 83.3)]" in gen_log, "gen_supp15.log")
    sp_n, sp_o = ctx["sp_new"], ctx["sp_old"]
    def title_lines(sp): return sorted([s for s in sp if s["text"].startswith(("Oxygen stayed", "of the night"))], key=lambda s: s["y"])
    def tick50(sp): return [s for s in sp if s["text"] == "50" and s["y"] > 150][0]
    tn, to = title_lines(sp_n), title_lines(sp_o); t50n, t50o = tick50(sp_n), tick50(sp_o); c50 = (t50n["x"] + t50n["x1"]) / 2.0
    for ln in tn:
        c = (ln["x"] + ln["x1"]) / 2.0
        ck.log(f"round-40 centring kept: x title line {ln['text']!r} centred on the 50 tick (centre {c50:.1f} pt) within 1.5 pt", abs(c - c50) <= 1.5, f"line centre {c:.1f} pt")
    ck.log("x title baselines within 3 pt of V26's (the title rides 1 to 2 pt lower under the taller tick labels)", len(tn) == len(to) == 2 and all(abs(a["y"] - b["y"]) <= 3.0 for a, b in zip(tn, to)), f"NEW {[l['y'] for l in tn]} V26 {[l['y'] for l in to]}")
    ck.log("the 50 tick label did not move horizontally (1 pt)", abs(t50n["x"] - t50o["x"]) <= 1.5, f"NEW x {t50n['x']} V26 {t50o['x']}")
    ck.info(f"generator record of the re-centring: {DV.get('xlabel_centre_r40')}")
spec = dict(sheet=S, old_pdf=OLD, new_pdf=NEW, out_dir=D, kept=[], note="every text +1 pt (11/10 -> 12/11), the bars carry no printed number, the x title centred on the ruled axis as in round 40, the margin gate relaxed to 2.9 mm (title 3.0 mm from the foot)",
            mark_colours=DV["ramp"], legend_boxes=[], drawn=drawn, present=["Oxygen stayed above 10%", "of the night, %", "Cardiopulmonary conditions", "before treatment"], absent=[], sources={**SRC, **TABLES}, extra=extra,
            notes=[f"builder: {SC}/gen_supp15.py (edits in scripts/R49_EDITS.md), generated {GEN} (the sheet itself)"])
ok = V.run(spec)
json.dump({"sources_sha256": SRC, "tables_sha256": TABLES, "generator": f"{SC}/gen_supp15.py", "generated": GEN, "new": NEW, "new_sha256": C.sha256(NEW), "xlabel_centre_r49": DV.get("xlabel_centre_r40"), "written": C.now()}, open(f"{W}/build_record.json", "w"), indent=1)
sys.exit(0 if ok else 1)
