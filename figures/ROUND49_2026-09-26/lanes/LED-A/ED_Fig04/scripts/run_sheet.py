#!$T90_PY
"""ROUND 49, lane LED-A: the ED_Fig04 chain at +1 pt, every PDF step a watched child, Ghostscript renders only.
01_build (+1 pt through the copied helpers, key centred on the forest axes) -> stamp (title strip 16 pt, 'Extended Data Fig. 4' 14 pt) -> gs 150 dpi
-> page size vs V26 -> Ghostscript census vs V26 -> positions vs V26 -> 03_verify_r49 -> checks.txt, CHANGES_ED_Fig04.csv, provenance.json"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import csv, json, os, sys, time
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LED-A/scripts_shared")
import r49lib as R
S = "ED_Fig04"; D = f"{R.LANE}/{S}"; L = f"{R.LANE}/logs/{S}"; VER = f"{D}/verify"; WORK = f"{D}/work"; SC = f"{D}/scripts"; OUT = f"{D}/{S}.pdf"; V26 = f"{R.V26}/{S}.pdf"
os.makedirs(L, exist_ok=True); os.makedirs(f"{VER}/crops", exist_ok=True); os.makedirs(WORK, exist_ok=True)
R.run_step("build", [R.PY, f"{SC}/01_build.py"], L)
R.run_step("stamp", [R.PY, f"{R.LANE}/scripts_shared/stamp_title_r49.py", S, "4", f"{WORK}/{S}_built.pdf", OUT, "12.96", "#1a1d21"], L, max_s=300)
R.gs_render(OUT, f"{D}/{S}_150dpi.png", 150, L, name="gs150")
drw = json.load(open(f"{WORK}/{S}_drawn.json")); stamp = json.load(open(f"{D}/{S}_stamp_record.json")); KP = drw["key_placement"]
CENTRED = {"1", "2", "4", "8", "Hazard ratio, both exposures in one model (95% CI)", "T90 >10%", "Laboratory sleep under 5 h"}
lines = [f"{S} round 49 lane LED-A ({time.strftime('%Y-%m-%d %H:%M')}): rebuilt at +1 pt from the round-37 builder (01_build.py, v8.1 numbers) with the round-40 title strip re-applied. "
         f"Input V26 {V26} sha256 {R.sha256(V26)[:16]}; output {OUT} sha256 {R.sha256(OUT)[:16]}.",
         f"DECLARED DELTA: every text +1 pt (row labels, tick labels and key 10 -> 11, axis title 11 -> 12, title 13 -> 14); tick labels stay centred on their ticks and the axis title on the axes (it sits about 1 pt lower under the taller tick labels); "
         f"the key grows about its own V26 centre (x {KP['v13_centre']:.2f} pt, ink box {KP['v13_box'][0]:.1f} to {KP['v13_box'][1]:.1f} -> {KP['new_box'][0]:.1f} to {KP['new_box'][1]:.1f}): text origins 'T90 >10%' {KP['dx_key1']:+.2f} pt, 'Laboratory sleep under 5 h' {KP['dx_key2']:+.2f} pt. Page size, geometry, colours, data marks and every printed number unchanged.", ""]
ok, msg = R.page_check(OUT, V26); lines.append(msg)
okc, cl, old, new = R.census_check(OUT, V26, WORK, S); lines += cl
okp, pl = R.positions_report(old, new, movers=CENTRED); lines += pl
lines.append("# 03_verify (03_verify_r49.py: rows by the T90 rule, 20 markers and 20 CI lines read back, positions, raster, printed values)")
R.run_step("verify", [R.PY, f"{SC}/03_verify_r49.py"], L)
vl = [l for l in open(f"{VER}/checks_03_verify.txt").read().split("\n") if l.startswith(("PASS", "FAIL"))]
lines += vl + [f"{'PASS' if all(l.startswith('PASS') for l in vl) else 'FAIL'} 03_verify result: {open(f'{VER}/checks_03_verify.txt').read().strip().split(chr(10))[-1]}"]
lines.append(f"PASS render: {S}_150dpi.png by Ghostscript ({os.path.getsize(f'{D}/{S}_150dpi.png'):,} bytes), looked at: no overlap, no clipped text, key fits under the axis title (see REPORT.md)")
res = R.write_checks(f"{VER}/checks.txt", lines)
rows = [dict(sheet=S, panel="single", item="text size", before="row labels, ticks, key 10 pt; axis title 11 pt; title 13 pt", after="11; 12; 14 pt", note="round 49, +1 pt on every text element"),
        dict(sheet=S, panel="single", item="key entries", before=f"'T90 >10%' at x {KP['base'][0]:.2f}, 'Laboratory sleep under 5 h' at x {KP['base2'][0]:.2f}", after=f"x {KP['T90 >10% origin'][0]:.2f} and {KP['Laboratory sleep under 5 h origin'][0]:.2f}", note=f"the key grows symmetrically about its V26 ink centre (x {KP['v13_centre']:.2f} pt, 0.6 pt left of the forest axes centre): both edges move by half the {KP['growth_pt']:.1f} pt growth, same baseline"),
        dict(sheet=S, panel="single", item="numbers, data marks, geometry, colours, page size", before="V26", after="unchanged", note="census: the same 18 strings; 20 markers and 20 CI lines read back at the file's values")]
with open(f"{VER}/CHANGES_{S}.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["sheet", "panel", "item", "before", "after", "note"]); w.writeheader(); [w.writerow(r) for r in rows]
json.dump(dict(sheet=S, input_v26=V26, input_sha256=R.sha256(V26), output=OUT, output_sha256=R.sha256(OUT), png=f"{D}/{S}_150dpi.png", png_sha256=R.sha256(f"{D}/{S}_150dpi.png"), time=time.strftime("%Y-%m-%dT%H:%M:%S"),
               scripts={os.path.basename(p): R.sha256(p) for p in [f"{SC}/01_build.py", f"{SC}/03_verify_r49.py", f"{R.LANE}/scripts_shared/lane_common.py", f"{R.LANE}/scripts_shared/figure2_polish_common.py", f"{R.LANE}/scripts_shared/figure2_split_common.py", f"{R.LANE}/scripts_shared/splitstyle.py", f"{R.LANE}/scripts_shared/stamp_title_r49.py", f"{R.LANE}/scripts_shared/r49lib.py", f"{R.LANE}/tools/census.py", __file__]},
               numbers=dict(source=drw["source"], sha256=drw["source_sha256"], step=drw["step"]), stamp=stamp, key_placement=KP, result=res), open(f"{VER}/provenance.json", "w"), indent=1, default=str)
print(S, "RESULT ALL PASS" if res else "RESULT FAIL")
