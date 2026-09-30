#!$T90_PY
"""Supp_Fig20, round 50 (lane LSUPP-C, Alen's item 6): panel d's lower bars without the values printed at their ends. The sub-sheets are
rebuilt by this lane's copy of the round-49 builder (one edit: the eight floor values are not drawn), composed by r49_compose_s20.py
(PyMuPDF placement, child under the watchdog) at the round-37 offsets with the letters at 14 pt, verified against V27 with r50_verify:
the eight strings removed, everything else the same at the same size, the sheet pixel-identical outside the eight value boxes, marks
identical, the round-49 key checks kept. Delivered as the VECTOR sheet (the coordinator rasterizes)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import json, os, re, subprocess, sys
from collections import Counter
L = f"{paths.FIGURE_ROOT}/ROUND50_2026-09-26/lanes/LSUPP-C"
sys.path.insert(0, f"{L}/scripts")
import r49_common as C
import r50_verify as V
S = "Supp_Fig20"; D = f"{L}/{S}"; W = f"{D}/work"; os.makedirs(W, exist_ok=True)
OLD = f"{C.V27}/{S}.pdf"; NEW = f"{D}/{S}.pdf"
SS = f"{paths.SV_ROOT}/stage_specific"
DVJ = f"{W}/eFigure12_drawn_values_v8_2.json"; DV = json.load(open(DVJ))
SRC = {k: C.sha256(v["path"]) for k, v in DV["_provenance"]["inputs"].items()}; REC = {k: v["sha256"] for k, v in DV["_provenance"]["inputs"].items()}
for f in ("stage_general_common_core.csv", "attack_summary.json", "attack_controls_full.csv", "attack_controls_full_negpanel_v5.csv"): C.hydrated(f"{SS}/{f}.provenance.json")
if "--no-gen" not in sys.argv:
    r = subprocess.run([C.PY, f"{D}/scripts/01_build_subsheets.py", "--new"], cwd=W, capture_output=True, text=True)
    open(f"{W}/build_subsheets.log", "w").write(r.stdout + "\n--- stderr ---\n" + r.stderr); assert r.returncode == 0, r.stderr[-3000:]
BR = json.load(open(f"{W}/subsheets_new/build_record.json"))
R49BR = json.load(open(f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LSUPP-C/Supp_Fig20/work/subsheets_new/build_record.json"))
C.wd(f"{S}_compose", [C.PY, f"{D}/scripts/r49_compose_s20.py"])
CR = json.load(open(f"{W}/compose_record.json"))
drawn = Counter()
for r in BR["A"]["printed"]: drawn[r["outcome"]] += 2; drawn[r["rem_hr"]] += 1
for r in BR["B"]["printed"]: drawn[r["label"]] += 1; drawn[r["hr"]] += 1; drawn[r["p"]] += 1
for r in BR["C"]["printed"]: drawn[r["disease"]] += 1; drawn[r["hr"]] += 1; drawn[r["p"]] += 1
for r in DV["D_scatter"]:
    for part in (["Urinary tract", "infection"] if r["label"] == "Urinary tract infection" else [r["label"]]): drawn[part] += 1
assert BR["D"]["floors_printed"] == {}, BR["D"]["floors_printed"]
FLOORS = [BR["D"]["floor_values_3dp"][k] for k, _l in (("all", 0), ("wake", 0), ("sleep", 0), ("nrem", 0), ("n1", 0), ("n2", 0), ("n3", 0), ("rem", 0))]
sp27 = C.gs_text(OLD, f"{W}/{S}_V27_txtwrite_pre.xml", f"{S}_txt_old_pre")
VAL = sorted([s for s in sp27 if re.fullmatch(r"\d\.\d{3}", s["text"]) and s["x"] > 600 and 1150 < s["y"] < 1400], key=lambda s: s["y"])
assert [s["text"] for s in VAL] == FLOORS, ([s["text"] for s in VAL], FLOORS)
VALUE_BOXES = [[s["x"] - 2.0, s["y"] - 12.0, s["x1"] + 2.0, s["y"] + 4.0] for s in VAL]
REM = [(v, VAL[i]["size"]) for i, v in enumerate(FLOORS)]
KEYS = [("REM", "#0288d1", 560, 1100, 595, 625), ("Whole-sleep", "#8a9099", 560, 1100, 595, 625), ("Healthy subgroup", "#0288d1", 0, 556, 1300, 1350), ("General cohort", "#8a9099", 0, 556, 1300, 1350), ("REM has", "#0288d1", 560, 1100, 1060, 1095)]
def extra(ck, ctx):
    ck.log("stage_specific inputs unchanged since the drawn-values record (sha256 equal the record's provenance)", SRC == REC, str({k: v[:12] for k, v in SRC.items()}))
    same = BR["A"]["printed"] == R49BR["A"]["printed"] and BR["B"]["printed"] == R49BR["B"]["printed"] and BR["C"]["printed"] == R49BR["C"]["printed"] and BR["D"]["floor_values_3dp"] == R49BR["D"]["floors_printed"]
    ck.log("the builder's own gates passed (rc 0) and its record equals the round-49 record string for string (panels a, b, c and the floor values, now drawn as bar ends only)", same, f"floors {BR['D']['floor_values_3dp']}, D axis {BR.get('D_axis')}")
    sp = ctx["sp_new"]
    def letters(s_): return sorted([(t["text"], round(t["x"], 2), round(t["y"], 2), t["size"]) for t in s_ if len(t["text"]) == 1 and t["text"] in "abcd" and t["size"] >= 12], key=lambda t: t[0])
    Lo, Ln = letters(ctx["sp_old"]), letters(sp)
    ck.log("the four letters are on NEW at 14 pt at the V27 origins (0.5 pt)", Lo == Ln, f"V27 {Lo} NEW {Ln}")
    for prefix, col, x0, x1, y0, y1 in KEYS:
        c = sorted([t for t in sp if t["text"].startswith(prefix) and x0 <= t["x"] <= x1 and y0 <= t["y"] <= y1 and not t["text"].startswith("Healthy subgroup,")], key=lambda t: t["x"])
        if not c: ck.log(f"key entry '{prefix}' found on NEW", False, ""); continue
        t = c[0]; box = [t["x"] - 48, t["y"] - 10, t["x"] - 2, t["y"] + 4]
        tall, allc = C.blobs(ctx["R"]["new_noaa"], 150, box, col, min_h_px=8)
        ck.log(f"round-40 key check kept: key entry '{prefix}' carries one marker blob of {col}", len(tall) == 1, f"tall {tall}, all {allc}")
    ck.info(f"compose: sub-sheets at {[(k, v['offset']) for k, v in CR['placed'].items()]}, letters {CR['stamps']}")
spec = dict(sheet=S, old_pdf=OLD, new_pdf=NEW, out_dir=D, removed=REM, added=[], note="item 6: the eight values at the ends of panel d's lower bars are not printed, bars and axis unchanged, nothing else changes",
            regions=[dict(name="whole sheet", box=[0, 0, 1089.48, 1441.46], allowed=VALUE_BOXES, note="identical outside the eight value boxes")],
            mark_colours=["#0288d1", "#8a9099", "#ccd1d6"], legend_boxes=[], drawn=drawn, present=["Level the negative controls reach, hazard ratio per 1 SD"], absent=[], sources=SRC, extra=extra,
            notes=[f"builder: {D}/scripts/01_build_subsheets.py (edit in scripts/R50_EDITS.md), composer {D}/scripts/r49_compose_s20.py", "VECTOR sheet delivered (the coordinator rasterizes)"])
ok = V.run(spec)
json.dump({"sources_sha256": SRC, "drawn_values": DVJ, "drawn_values_sha256": C.sha256(DVJ), "builder": f"{D}/scripts/01_build_subsheets.py", "compose": CR, "removed": REM, "new": NEW, "new_sha256": C.sha256(NEW), "written": C.now()}, open(f"{W}/build_record.json", "w"), indent=1)
sys.exit(0 if ok else 1)
