#!$T90_PY
"""Supp_Fig14, round 49 (lane LSUPP-C): the V26 sheet (round 40) rebuilt with every text one point larger. gen_supp14.py (lane
copy: paths, PT_PLUS on TITLE_PT/TICK_PT/ANN_PT, the printed-column widths that set the axes rect measured at the round-40 size so
every mark stays put) is run plain; the 171 mm figure is composed onto the V26 page box (486.51 x 317.07, figure at (1.70, -0.57),
the V13 idiom) by PyMuPDF in a child under the watchdog; verified against V26 with r49_verify plus the round-40 key-marker and
colour checks."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import json, os, subprocess, sys
from collections import Counter
L = f"{paths.FIGURE_ROOT}/ROUND50_2026-09-26/lanes/LSUPP-C"
sys.path.insert(0, f"{L}/scripts")
import r49_common as C
import r49_verify as V
S = "Supp_Fig14"; D = f"{L}/{S}"; W = f"{D}/work"; SC = f"{L}/L5b_scripts"; WK = f"{L}/L5b_work"; os.makedirs(W, exist_ok=True)
OLD = f"{C.V26}/{S}.pdf"; NEW = f"{D}/{S}.pdf"
NUM = paths.NUMBERS_DIR
SRC = {f: C.sha256(f"{NUM}/{f}") for f in ("nonresponder_phenotype_v2.csv", "nonresponder_combined_v2.json", "treatment_v2.json", "negcontrols_final.json")}
for f in SRC: C.hydrated(f"{NUM}/{f}.provenance.json")
R40 = json.load(open(f"{paths.FIGURE_ROOT}/ROUND40_2026-09-18/lanes/LSUPP/Supp_Fig14/work/build_record.json"))["sources_sha256"]
if "--no-gen" not in sys.argv:
    r = subprocess.run([C.PY, f"{SC}/gen_supp14.py"], cwd=WK, capture_output=True, text=True)
    open(f"{W}/gen_supp14.log", "w").write(r.stdout + "\n--- stderr ---\n" + r.stderr); assert r.returncode == 0, r.stderr[-3000:]
GEN = f"{WK}/polish_out/pdf/Figure5D.pdf"; fw, fh = C.page_box(GEN); assert abs(fw - 484.72) < 0.2, fw
w26, h26 = C.page_box(OLD); assert abs(w26 - 486.51) < 0.05 and abs(h26 - 317.07) < 0.05, (w26, h26)
OFFX, OFFY = 1.70, -0.57
job = {"mode": "compose", "out": NEW, "page": [w26, h26], "parts": [{"pdf": GEN, "rect": [OFFX, OFFY, OFFX + fw, OFFY + fh]}], "letters": [], "metadata": {"title": S, "creator": "LSUPP-C r49 gen_supp14.py (matplotlib, +1 pt) + fitz_place compose"}}
json.dump(job, open(f"{W}/place_job.json", "w")); C.wd(f"{S}_place", [C.PY, f"{L}/scripts/fitz_place.py", f"{W}/place_job.json"])
DV = json.load(open(f"{WK}/polish_out/work/Figure5D_polish_drawn_values.json")); assert len(DV["conditions"]) == 10
drawn = Counter()
for c in DV["conditions"]:
    for k in ("condition", "printed", "printed_q"): drawn[c[k]] += 1
KEY_SHIFT_1, KEY_SHIFT_2 = -16.0, -8.0   # the key nudge of gen_supp14.py (R49)
KEY_BLUE = [131.0 + KEY_SHIFT_1, 10.0, 157.5 + KEY_SHIFT_1, 26.0]; KEY_ORANGE = [243.0 + KEY_SHIFT_2, 10.0, 269.5 + KEY_SHIFT_2, 26.0]
KEY_STRIP = [104.0, 6.0, 280.0, 30.0]   # the key strip: its handles and markers moved 8 pt left (declared), the only mark-coloured pixels allowed to differ
gen_log = open(f"{W}/gen_supp14.log").read()
def extra(ck, ctx):
    ck.log("numbers files unchanged since the round-40 build (sha256 of the four numbers files equal the round-40 build record)", SRC == R40, str({k: v[:12] for k, v in SRC.items()}))
    ck.log("the builder's own data gates passed (rc 0: row set and order from the files, every marker and interval read back against the frozen values, printed strings half up on the v8.3 json)", "drawn, in order:" in gen_log, "gen_supp14.log")
    for name, box, col in (("Cardiopulmonary", KEY_BLUE, "#0288d1"), ("Not cardiopulmonary", KEY_ORANGE, "#d55e00")):
        tall, allc = C.blobs(ctx["R"]["new_noaa"], 150, box, col, min_h_px=8)
        ck.log(f"round-40 key check kept: key entry '{name}' carries one marker blob of {col}", len(tall) == 1, f"tall components {tall}, all {allc}")
    ck.log("round-40 colour kept: the palette orange #d55e00 is painted and the old #b5623a is absent", "#d55e00" in ctx["cen_new"] and "#b5623a" not in ctx["cen_new"], "")
    sp = ctx["sp_new"]; top = sorted([t for t in sp if t["y"] < 40], key=lambda t: t["x"])
    k1 = [t for t in top if t["text"].startswith("Cardiopulmonary")][0]; k2 = [t for t in top if t["text"].startswith("Not")][0]; h = [t for t in top if t["text"].startswith("OR")][0]
    ck.log("the nudged key at 11 pt keeps at least 10 pt to the column header and to the next handle", (h["x"] - k2["x1"]) >= 10.0 and (248.78 + KEY_SHIFT_2 - k1["x1"]) >= 10.0, f"key 1 ends {k1['x1']:.1f} (orange handle at {248.78 + KEY_SHIFT_2:.1f}), key 2 ends {k2['x1']:.1f}, header starts {h['x']:.1f}")
    ck.info("the axes rect and every mark stay where V26 has them: the printed-column widths that set the plot's right edge are measured at the round-40 size (gen_supp14.py edit), proved by the marks check above")
spec = dict(sheet=S, old_pdf=OLD, new_pdf=NEW, out_dir=D, kept=[], note="every text +1 pt (11/10/9.5 -> 12/11/10.5), geometry pinned, the two-entry key of the top strip nudged left (entry 1 by 16 pt, entry 2 by 8 pt: at 11 pt in place its labels were 8 pt from the next handle and 4 pt from the column header)",
            mark_colours=["#0288d1", "#d55e00"], legend_boxes=[KEY_STRIP], drawn=drawn, present=["Cardiopulmonary", "Not cardiopulmonary", "OR (95% CI)"], absent=[], sources=SRC, extra=extra,
            notes=[f"builder: {SC}/gen_supp14.py (edits in scripts/R49_EDITS.md), generated {GEN}", "round-40 steps re-applied in the generator copy: one marker per key entry, second series #d55e00, row order by the 3 dp odds ratio (heart failure before pulmonary hypertension)"])
ok = V.run(spec)
json.dump({"sources_sha256": SRC, "generator": f"{SC}/gen_supp14.py", "generated": GEN, "new": NEW, "new_sha256": C.sha256(NEW), "written": C.now()}, open(f"{W}/build_record.json", "w"), indent=1)
sys.exit(0 if ok else 1)
