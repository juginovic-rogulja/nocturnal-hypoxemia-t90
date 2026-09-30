#!$T90_PY
"""Supp_Fig20, round 49 (lane LSUPP-C): the V26 vector sheet (round 40) rebuilt with every text one point larger. 01_build_subsheets.py
(lane copy: paths, PT_PLUS on LABF/TCKF/ANNF/PTLAB, the row labels of panels a, b and c kept at 10 pt, panel d's title band pinned
at the round-40 value, the two right-growing column headers exempt from the sub-sheet edge gate) rebuilds the four sub-sheets from the
copied v8.2 record; r49_compose_s20.py (PyMuPDF placement, child under the watchdog) places them at the round-37 solved offsets on the
V26 page and stamps the letters at 14 pt at the V26 origins; verified against V26 with r49_verify plus the round-40 key checks.
Delivered as the VECTOR sheet (the coordinator rasterizes)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import json, os, subprocess, sys
from collections import Counter
L = f"{paths.FIGURE_ROOT}/ROUND50_2026-09-26/lanes/LSUPP-C"
sys.path.insert(0, f"{L}/scripts")
import r49_common as C
import r49_verify as V
S = "Supp_Fig20"; D = f"{L}/{S}"; W = f"{D}/work"; os.makedirs(W, exist_ok=True)
OLD = f"{C.V26}/{S}.pdf"; NEW = f"{D}/{S}.pdf"
SS = f"{paths.SV_ROOT}/stage_specific"
DVJ = f"{W}/eFigure12_drawn_values_v8_2.json"; DV = json.load(open(DVJ))
SRC = {k: C.sha256(v["path"]) for k, v in DV["_provenance"]["inputs"].items()}; REC = {k: v["sha256"] for k, v in DV["_provenance"]["inputs"].items()}
for f in ("stage_general_common_core.csv", "attack_summary.json", "attack_controls_full.csv", "attack_controls_full_negpanel_v5.csv"): C.hydrated(f"{SS}/{f}.provenance.json")
if "--no-gen" not in sys.argv:
    r = subprocess.run([C.PY, f"{D}/scripts/01_build_subsheets.py", "--new"], cwd=W, capture_output=True, text=True)
    open(f"{W}/build_subsheets.log", "w").write(r.stdout + "\n--- stderr ---\n" + r.stderr); assert r.returncode == 0, r.stderr[-3000:]
BR = json.load(open(f"{W}/subsheets_new/build_record.json"))
R40BR = json.load(open(f"{paths.FIGURE_ROOT}/ROUND40_2026-09-18/lanes/LSUPP/Supp_Fig20/work/subsheets_new/build_record.json"))
C.wd(f"{S}_compose", [C.PY, f"{D}/scripts/r49_compose_s20.py"])
CR = json.load(open(f"{W}/compose_record.json"))
drawn = Counter()
for r in BR["A"]["printed"]: drawn[r["outcome"]] += 2; drawn[r["rem_hr"]] += 1
for r in BR["B"]["printed"]: drawn[r["label"]] += 1; drawn[r["hr"]] += 1; drawn[r["p"]] += 1
for r in BR["C"]["printed"]: drawn[r["disease"]] += 1; drawn[r["hr"]] += 1; drawn[r["p"]] += 1
for k, v in BR["D"]["floors_printed"].items(): drawn[v] += 1
for r in DV["D_scatter"]:
    for part in (["Urinary tract", "infection"] if r["label"] == "Urinary tract infection" else [r["label"]]): drawn[part] += 1
KEPT = [(r["outcome"], 10.0, 2) for r in BR["A"]["printed"]] + [(r["label"], 10.0, 1) for r in BR["B"]["printed"]] + [(r["disease"], 10.0, 1) for r in BR["C"]["printed"]]
for r in DV["D_scatter"]:
    for part in (["Urinary tract", "infection"] if r["label"] == "Urinary tract infection" else [r["label"]]): KEPT.append((part, 9.5, 1))   # the scatter-point labels of panel d keep 9.5 pt
REM = ["shaded zone, below 1.25:", "no window-versus-window", "claim is made", "Another window has the largest hazard ratio"]
PRES = ["Another sleep stage has the largest hazard ratio", "REM T90, adjusted for whole-sleep T90", "Whole-sleep T90, adjusted for REM T90", "Healthy subgroup", "General cohort", "REM has the largest hazard ratio"]
KEYS = [("REM", "#0288d1", 560, 1100, 595, 625), ("Whole-sleep", "#8a9099", 560, 1100, 595, 625), ("Healthy subgroup", "#0288d1", 0, 556, 1300, 1350), ("General cohort", "#8a9099", 0, 556, 1300, 1350), ("REM has", "#0288d1", 560, 1100, 1060, 1095)]   # (first words of the entry, handle colour, x range, y range of the key zone)
def extra(ck, ctx):
    ck.log("stage_specific inputs unchanged since the drawn-values record (sha256 of core, panel, nested, bootstrap, summary and healthy files equal the record's provenance)", SRC == REC, str({k: v[:12] for k, v in SRC.items()}))
    same = BR["A"]["n"] == R40BR["A"]["n"] and BR["B"]["printed"] == R40BR["B"]["printed"] and BR["C"]["printed"] == R40BR["C"]["printed"] and BR["D"]["floors_printed"] == R40BR["D"]["floors_printed"] and BR["A"]["printed"] == R40BR["A"]["printed"]
    ck.log("the builder's own gates passed (rc 0) and its printed record equals the round-40 record string for string (panels a, b, c, d)", same, f"B {BR['B']['rem_adds_p05']} of {BR['B']['n']} REM adds, C starred {BR['C']['starred']}, D floors {BR['D']['floors_printed']}, D axis {BR.get('D_axis')}")
    sp = ctx["sp_new"]
    def letters(s_): return sorted([(t["text"], round(t["x"], 2), round(t["y"], 2), t["size"]) for t in s_ if len(t["text"]) == 1 and t["text"] in "abcd" and t["size"] >= 12], key=lambda t: t[0])
    Lo, Ln = letters(ctx["sp_old"]), letters(sp)
    ck.log("the four letters are on NEW at 14 pt at the V26 origins (x and baseline within 0.5 pt)", [t[0] for t in Lo] == [t[0] for t in Ln] == ["a", "b", "c", "d"] and all(abs(a[1] - b[1]) <= 0.5 and abs(a[2] - b[2]) <= 0.5 and b[3] == 14.0 for a, b in zip(Lo, Ln)), f"V26 {Lo} NEW {Ln}")
    for prefix, col, x0, x1, y0, y1 in KEYS:
        c = sorted([t for t in sp if t["text"].startswith(prefix) and x0 <= t["x"] <= x1 and y0 <= t["y"] <= y1 and not t["text"].startswith("Healthy subgroup,")], key=lambda t: t["x"])
        if not c: ck.log(f"key entry '{prefix}' found on NEW", False, ""); continue
        t = c[0]; box = [t["x"] - 48, t["y"] - 10, t["x"] - 2, t["y"] + 4]   # the round-40 handle box (the handle sits 7 to 31 pt left of the text)
        tall, allc = C.blobs(ctx["R"]["new_noaa"], 150, box, col, min_h_px=8)
        ck.log(f"round-40 key check kept: key entry '{prefix}' carries one marker blob of {col}", len(tall) == 1, f"box {[round(v, 1) for v in box]}, tall {tall}, all {allc}")
    ck.info(f"compose: sub-sheets at {[(k, v['offset']) for k, v in CR['placed'].items()]}, letters {CR['stamps']}")
spec = dict(sheet=S, old_pdf=OLD, new_pdf=NEW, out_dir=D, kept=KEPT,
            note="every text +1 pt (13/11/10/9.5 -> 14/12/11/10.5) except the row labels of panels a, b and c (kept at 10 pt: at 11 pt 'Peripheral artery disease' would cross panel a's edge by 2.8 pt onto a mark 4.3 pt from that edge and 'Cardiovascular composite' would sit 3.4 pt on the band strip of b and c); panel d's window labels are 11 pt; the 16 scatter-point labels of panel d keep 9.5 pt (the no-overlap solver has no clear place for 'Death from any cause' at 10.5 pt); the column headers of a and b grow into the sub-sheet margins (3.2 and 3.6 mm from the sub-sheet edge, 8.5 mm or more from the page edge)",
            mark_colours=["#0288d1", "#8a9099", "#ccd1d6"], marks_tolerance=0, drawn=drawn, present=PRES, absent=REM, sources=SRC, extra=extra,
            notes=[f"builder: {D}/scripts/01_build_subsheets.py (edits in scripts/R49_EDITS.md), composer {D}/scripts/r49_compose_s20.py", "VECTOR sheet delivered (the brief: the coordinator rasterizes)"])
def legend_boxes_from(sp_old, sp_new):
    """The legend zones whose handles move with the taller legends (declared): panel b's key, panel c's key, panel d's key and the
    N3 keys of panels a and d, each found by the first words of its entries inside its zone (txtwrite splits the entries at kerning)."""
    zones = [(("REM", "Whole-sleep"), 560, 1100, 595, 625), (("Healthy subgroup", "General cohort"), 0, 556, 1300, 1350),
             (("REM has", "Another"), 560, 1100, 1060, 1095), (("N3, too little",), 0, 556, 688, 706), (("N3, too little",), 560, 1100, 1402, 1420)]
    out = []
    for firsts, x0, x1, y0, y1 in zones:
        for sp in (sp_old, sp_new):
            zone = [s for s in sp if x0 <= s["x"] <= x1 and y0 <= s["y"] <= y1 and not s["text"].startswith("Healthy subgroup,")]
            b = V.box_around(zone, lambda t: t.startswith(firsts), pad_left=64, pad_right=6, pad_y=9)
            assert b, ("legend zone without text", firsts, x0, x1, y0, y1)
            out.append(b)
    return out
_run = V.run
def run_with_boxes(spec):
    sp_old = C.gs_text(OLD, f"{W}/{S}_V26_txtwrite_pre.xml", f"{S}_txt_old_pre"); sp_new = C.gs_text(NEW, f"{W}/{S}_NEW_txtwrite_pre.xml", f"{S}_txt_new_pre")
    spec["legend_boxes"] = legend_boxes_from(sp_old, sp_new); print("legend boxes", len(spec["legend_boxes"])); return _run(spec)
ok = run_with_boxes(spec)
json.dump({"sources_sha256": SRC, "drawn_values": DVJ, "drawn_values_sha256": C.sha256(DVJ), "builder": f"{D}/scripts/01_build_subsheets.py", "compose": CR, "new": NEW, "new_sha256": C.sha256(NEW), "written": C.now()}, open(f"{W}/build_record.json", "w"), indent=1)
sys.exit(0 if ok else 1)
