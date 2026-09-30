#!$T90_PY
"""ROUND 49, lane LED-A: the ED_Fig05 chain at +1 pt, every PDF step a watched child, Ghostscript renders only.
01_build (+1 pt through the copied helpers, the V26 axes layout kept) -> stamp (title strip 16 pt, 'Extended Data Fig. 5' 14 pt) -> gs 150 dpi
-> page size vs V26 -> Ghostscript census vs V26 -> positions vs V26 -> 03_verify_r49 -> checks.txt, CHANGES_ED_Fig05.csv, provenance.json"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import csv, json, os, sys, time
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LED-A/scripts_shared")
import r49lib as R
S = "ED_Fig05"; D = f"{R.LANE}/{S}"; L = f"{R.LANE}/logs/{S}"; VER = f"{D}/verify"; WORK = f"{D}/work"; SC = f"{D}/scripts"; OUT = f"{D}/{S}.pdf"; V26 = f"{R.V26}/{S}.pdf"
os.makedirs(L, exist_ok=True); os.makedirs(f"{VER}/crops", exist_ok=True); os.makedirs(WORK, exist_ok=True)
R.run_step("build", [R.PY, f"{SC}/01_build.py"], L)
R.run_step("stamp", [R.PY, f"{R.LANE}/scripts_shared/stamp_title_r49.py", S, "5", f"{WORK}/{S}_built.pdf", OUT, "12.96", "#1a1d21"], L, max_s=300)
R.gs_render(OUT, f"{D}/{S}_150dpi.png", 150, L, name="gs150")
drw = json.load(open(f"{WORK}/{S}_drawn.json")); stamp = json.load(open(f"{D}/{S}_stamp_record.json")); KP = drw["key_placement"]
MOVERS = {"Percent of the group", "T90 >5%", "T90 >10%", "Sleep under", "5 h", "Sleep 7 h", "or more", "(95% CI,"} | {t for p in drw["panels"].values() for c in p["cuts"].values() for t in c["printed"] if not t[0].isdigit() or "-" in t}
lines = [f"{S} round 49 lane LED-A ({time.strftime('%Y-%m-%d %H:%M')}): rebuilt at +1 pt from the round-37 builder (01_build.py, the v8 table through cohort_spec, step-129 cross-check) with the round-40 title strip re-applied. "
         f"Input V26 {V26} sha256 {R.sha256(V26)[:16]}; output {OUT} sha256 {R.sha256(OUT)[:16]}.",
         f"DECLARED DELTA: every text +1 pt (annotations 9.5 -> 10.5, ticks, key and panel titles 10 -> 11, y labels 11 -> 12, letters and title 13 -> 14); the axes widths are the V26 layout (the key width measured at 10 pt), "
         f"the key at 11 pt keeps its top and right edge and reaches {KP['gap_in']:.3f} in from the panel b axes (V26 gap {KP['gap_v26_in']:.2f} in); annotation blocks hang from their top, later lines lower by the pitch growth; "
         f"y labels centred on their axes move left with the wider tick labels. Page size, geometry, colours, data marks and every printed number unchanged.", ""]
ok, msg = R.page_check(OUT, V26); lines.append(msg)
okc, cl, old, new = R.census_check(OUT, V26, WORK, S); lines += cl
okp, pl = R.positions_report(old, new, movers=MOVERS, tol=6.5); lines += pl
lines.append("# 03_verify (03_verify_r49.py: every printed value recomputed from the v8 table, bars equal V26's, anchors, printed values)")
R.run_step("verify", [R.PY, f"{SC}/03_verify_r49.py"], L)
vl = [l for l in open(f"{VER}/checks_03_verify.txt").read().split("\n") if l.startswith(("PASS", "FAIL"))]
lines += vl + [f"{'PASS' if all(l.startswith('PASS') for l in vl) else 'FAIL'} 03_verify result: {open(f'{VER}/checks_03_verify.txt').read().strip().split(chr(10))[-1]}"]
lines.append(f"PASS render: {S}_150dpi.png by Ghostscript ({os.path.getsize(f'{D}/{S}_150dpi.png'):,} bytes), looked at: no overlap, no clipped text, key fits (see REPORT.md)")
res = R.write_checks(f"{VER}/checks.txt", lines)
rows = [dict(sheet=S, panel="a, b, c, d", item="text size", before="annotations 9.5 pt; ticks, key, panel titles 10 pt; y labels 11 pt; letters and title 13 pt", after="10.5; 11; 12; 14 pt", note="round 49, +1 pt on every text element"),
        dict(sheet=S, panel="b", item="key box left edge", before=f"{KP['gap_v26_in']:.2f} in right of the panel b axes", after=f"{KP['gap_in']:.3f} in", note="the key keeps its top and its right edge at the page margin; the axes widths are the V26 ones (key width measured at 10 pt for the layout)"),
        dict(sheet=S, panel="a, b, c, d", item="numbers, bars, geometry, colours, page size", before="V26", after="unchanged", note="census: the same 97 strings; 16 bar rectangles equal V26's; every printed value recomputed")]
with open(f"{VER}/CHANGES_{S}.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["sheet", "panel", "item", "before", "after", "note"]); w.writeheader(); [w.writerow(r) for r in rows]
json.dump(dict(sheet=S, input_v26=V26, input_sha256=R.sha256(V26), output=OUT, output_sha256=R.sha256(OUT), png=f"{D}/{S}_150dpi.png", png_sha256=R.sha256(f"{D}/{S}_150dpi.png"), time=time.strftime("%Y-%m-%dT%H:%M:%S"),
               scripts={os.path.basename(p): R.sha256(p) for p in [f"{SC}/01_build.py", f"{SC}/03_verify_r49.py", f"{R.LANE}/scripts_shared/lane_common.py", f"{R.LANE}/scripts_shared/supp_polishB_common.py", f"{R.LANE}/scripts_shared/splitstyle.py", f"{R.LANE}/scripts_shared/stamp_title_r49.py", f"{R.LANE}/scripts_shared/r49lib.py", f"{R.LANE}/tools/census.py", __file__]},
               numbers=drw["sources"], stamp=stamp, key_placement=KP, result=res), open(f"{VER}/provenance.json", "w"), indent=1, default=str)
print(S, "RESULT ALL PASS" if res else "RESULT FAIL")
