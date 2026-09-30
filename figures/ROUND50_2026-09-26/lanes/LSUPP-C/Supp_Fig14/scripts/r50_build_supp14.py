#!$T90_PY
"""Supp_Fig14, round 50 (lane LSUPP-C, Alen's item 4): every row label black. gen_supp14.py (this lane's copy of the round-49 builder,
one edit: the label colours) is run plain, the figure composed onto the V27 page box (486.51 x 317.07, figure at (1.70, -0.57)) by
PyMuPDF in a child under the watchdog, verified against V27 with r50_verify: same strings at the same sizes, the whole sheet
pixel-identical outside the five recoloured label boxes, marks identical, the grey gone, the round-49 key checks kept."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import json, os, subprocess, sys
from collections import Counter
L = f"{paths.FIGURE_ROOT}/ROUND50_2026-09-26/lanes/LSUPP-C"
sys.path.insert(0, f"{L}/scripts")
import r49_common as C
import r50_verify as V
S = "Supp_Fig14"; D = f"{L}/{S}"; W = f"{D}/work"; SC = f"{L}/L5b_scripts"; WK = f"{L}/L5b_work"; os.makedirs(W, exist_ok=True)
OLD = f"{C.V27}/{S}.pdf"; NEW = f"{D}/{S}.pdf"
NUM = paths.NUMBERS_DIR
SRC = {f: C.sha256(f"{NUM}/{f}") for f in ("nonresponder_phenotype_v2.csv", "nonresponder_combined_v2.json", "treatment_v2.json", "negcontrols_final.json")}
for f in SRC: C.hydrated(f"{NUM}/{f}.provenance.json")
R49 = json.load(open(f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LSUPP-C/Supp_Fig14/work/build_record.json"))["sources_sha256"]
if "--no-gen" not in sys.argv:
    r = subprocess.run([C.PY, f"{SC}/gen_supp14.py"], cwd=WK, capture_output=True, text=True)
    open(f"{W}/gen_supp14.log", "w").write(r.stdout + "\n--- stderr ---\n" + r.stderr); assert r.returncode == 0, r.stderr[-3000:]
GEN = f"{WK}/polish_out/pdf/Figure5D.pdf"; fw, fh = C.page_box(GEN); assert abs(fw - 484.72) < 0.2, fw
w27, h27 = C.page_box(OLD); assert abs(w27 - 486.51) < 0.05 and abs(h27 - 317.07) < 0.05, (w27, h27)
OFFX, OFFY = 1.70, -0.57
job = {"mode": "compose", "out": NEW, "page": [w27, h27], "parts": [{"pdf": GEN, "rect": [OFFX, OFFY, OFFX + fw, OFFY + fh]}], "letters": [], "metadata": {"title": S, "creator": "LSUPP-C r50 gen_supp14.py (matplotlib) + fitz_place compose"}}
json.dump(job, open(f"{W}/place_job.json", "w")); C.wd(f"{S}_place", [C.PY, f"{L}/scripts/fitz_place.py", f"{W}/place_job.json"])
DV = json.load(open(f"{WK}/polish_out/work/Figure5D_polish_drawn_values.json")); assert len(DV["conditions"]) == 10
drawn = Counter()
for c in DV["conditions"]:
    for k in ("condition", "printed", "printed_q"): drawn[c[k]] += 1
NONCP = [c["condition"] for c in DV["conditions"] if not c["cardiopulmonary"]]; assert len(NONCP) == 5, NONCP
KEY_SHIFT_1, KEY_SHIFT_2 = -16.0, -8.0
KEY_BLUE = [131.0 + KEY_SHIFT_1, 10.0, 157.5 + KEY_SHIFT_1, 26.0]; KEY_ORANGE = [243.0 + KEY_SHIFT_2, 10.0, 269.5 + KEY_SHIFT_2, 26.0]
sp27 = C.gs_text(OLD, f"{W}/{S}_V27_txtwrite_pre.xml", f"{S}_txt_old_pre")
LABEL_BOXES = []
for name in NONCP:
    ss = [s for s in sp27 if s["text"] == name and s["x"] < 60]; assert len(ss) == 1, (name, len(ss))
    s = ss[0]; LABEL_BOXES.append([s["x"] - 2.0, s["y"] - 12.0, s["x1"] + 3.0, s["y"] + 4.0])
gen_log = open(f"{W}/gen_supp14.log").read()
def extra(ck, ctx):
    ck.log("numbers files unchanged since the round-49 build (sha256 equal the round-49 build record)", SRC == R49, str({k: v[:12] for k, v in SRC.items()}))
    ck.log("the builder's own data gates passed (rc 0: row set and order from the files, every marker and interval read back against the frozen values)", "drawn, in order:" in gen_log, "gen_supp14.log")
    ck.log("the grey #8a9099 (the former label colour) is on V27 and gone from NEW, nothing else changes colour", "#8a9099" in ctx["cen_old"] and "#8a9099" not in ctx["cen_new"], f"NEW colours {sorted(ctx['cen_new'])}")
    for name, box, col in (("Cardiopulmonary", KEY_BLUE, "#0288d1"), ("Not cardiopulmonary", KEY_ORANGE, "#d55e00")):
        tall, allc = C.blobs(ctx["R"]["new_noaa"], 150, box, col, min_h_px=8)
        ck.log(f"round-40 key check kept: key entry '{name}' carries one marker blob of {col}", len(tall) == 1, f"tall components {tall}, all {allc}")
    ck.log("round-40 colour kept: the palette orange #d55e00 is painted and the old #b5623a is absent", "#d55e00" in ctx["cen_new"] and "#b5623a" not in ctx["cen_new"], "")
spec = dict(sheet=S, old_pdf=OLD, new_pdf=NEW, out_dir=D, removed=[], added=[], note="item 4: the five 'not cardiopulmonary' row labels black instead of grey, nothing else changes",
            regions=[dict(name="whole sheet", box=[0, 0, w27, h27], allowed=LABEL_BOXES, note="identical outside the five recoloured label boxes")],
            mark_colours=["#0288d1", "#d55e00"], legend_boxes=[], drawn=drawn, present=["Cardiopulmonary", "Not cardiopulmonary", "OR (95% CI)"] + NONCP, absent=[], sources=SRC, extra=extra,
            notes=[f"builder: {SC}/gen_supp14.py (edit in scripts/R50_EDITS.md), generated {GEN}", f"recoloured labels: {NONCP}"])
ok = V.run(spec)
json.dump({"sources_sha256": SRC, "generator": f"{SC}/gen_supp14.py", "generated": GEN, "new": NEW, "new_sha256": C.sha256(NEW), "written": C.now()}, open(f"{W}/build_record.json", "w"), indent=1)
sys.exit(0 if ok else 1)
