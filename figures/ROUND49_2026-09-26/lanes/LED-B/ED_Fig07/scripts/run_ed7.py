#!$T90_PY
"""ED_Fig07 round 49 runner (lane LED-B): build a -> build b -> compose -> verify (children under wd_run.sh) -> title strip at 14 pt -> render -> census -> checks.txt."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import json, os, sys, time, shutil
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LED-B/scripts"); import ledb_common as C
S = "ED_Fig07"; N = 7; D = f"{C.LANE}/{S}"; WORK, VER, SCR = f"{D}/work", f"{D}/verify", f"{D}/scripts"; OUT = f"{D}/{S}.pdf"; V26 = f"{C.V26}/{S}.pdf"
os.makedirs(WORK, exist_ok=True); os.makedirs(f"{VER}/crops", exist_ok=True)
steps = sys.argv[1:] or ["build_a", "build_b", "compose", "verify", "stamp"]
if "build_a" in steps: C.wd("ed7_01_build_a_r49", 900, [C.PY, f"{SCR}/ed7_01_build_a.py"])
if "build_b" in steps: C.wd("ed7_02_build_b_r49", 900, [C.PY, f"{SCR}/ed7_02_build_b.py"])
if "compose" in steps: C.wd("ed7_03_compose_r49", 900, [C.PY, f"{SCR}/ed7_03_compose.py"])
if "verify" in steps: C.wd("ed7_04_verify_r49", 1800, [C.PY, f"{SCR}/ed7_04_verify_r49.py"])
COMPOSED = f"{D}/{S}.pdf"     # ed7spec OUT_PDF (the composed sheet at the V13 page size) is written by compose; the stamped sheet replaces it below
if "stamp" in steps:
    shutil.copyfile(COMPOSED, f"{WORK}/{S}_composed.pdf")
    C.wd(f"stamp_{S}", 300, [C.PY, f"{C.LANE}/scripts/stamp_title_r49.py", S, str(N), f"{WORK}/{S}_composed.pdf", OUT, str(C.V26_TITLE_X[S]), C.INK])
    C.gs_render(OUT, f"{D}/{S}_150dpi.png")
vlines = [l for l in open(f"{VER}/checks.txt").read().split("\n") if l.startswith(("PASS", "FAIL", "INFO"))]
keep = json.load(open(f"{WORK}/panel_b_geometry.json")).get("size_exceptions", {})
ex = dict(sizes=[dict(text=t, old=9.5, new=sz, why="does not fit its column at 10.5 pt: the largest size at or above 9.5 in 0.5 pt steps that fits (panel b, build_b fit_size)") for t, sz in keep.items()])
ok_c, clines = C.census(S, OUT, V26, VER, exceptions=ex)
ok_p, new_sz, old_sz = C.page_size_check(OUT, V26)
srec = json.load(open(f"{D}/{S}_stamp_record.json")); ga = json.load(open(f"{WORK}/panel_a_geometry.json"))
hdr = [f"{S} round 49 lane LED-B ({time.strftime('%Y-%m-%d %H:%M')}): rebuilt from the round-38 builders (copied, +1 pt) on the v8 numbers; V26 reference {V26} sha256 {C.sha256(V26)[:16]}; output {OUT} sha256 {C.sha256(OUT)[:16]}.",
       f"DECLARED DELTA: every text +1 pt (panel a labels 11, estimates and stars 10.5, axis title 12, 'HR (95% CI)' 10.5 right-aligned on the wider estimate column; panel b 10.5 except the fitted sizes {keep}; key 10.087; letters 14; title strip 16 pt with 'Extended Data Fig. 7' Arial Bold 14 at ({srec['origin'][0]}, {srec['origin'][1]})).",
       f"  the two star keys are not drawn (removed on V26 by lane LNOTES, round 40); panel a row labels past the axes box's invisible left edge: {ga['labels_past_axes_edge_pt']} (clear of the plot's leftmost ink at x {ga['ink_left_pt']:.1f} by 4 pt or more); numbers unchanged (re-derived from the numbers files below).", ""]
body = hdr + [f"{'PASS' if ok_p else 'FAIL'}  page size equals V26 (0.05 pt): {new_sz} vs {old_sz}", f"{'PASS' if abs(srec['page_out'][1] - srec['page_in'][1] - 16.0) < 0.01 and srec['size'] == 14.0 else 'FAIL'}  title strip: page {srec['page_in'][0]:.2f} x {srec['page_in'][1]:.3f} -> {srec['page_out'][0]:.2f} x {srec['page_out'][1]:.3f}, title '{srec['title']}' Arial Bold {srec['size']:g} pt at {srec['origin']}",
               "# builder and panel proofs (ed7_04_verify_r49.py, child under wd_run.sh, on the composed sheet before the strip)"] + vlines + ["# Ghostscript txtwrite census against V26 (ledb_common.census)"] + clines
C.finish_checks(f"{VER}/checks.txt", body)
json.dump(dict(sheet=S, v26=V26, v26_sha256=C.sha256(V26), output=OUT, output_sha256=C.sha256(OUT), png=f"{D}/{S}_150dpi.png", png_sha256=C.sha256(f"{D}/{S}_150dpi.png"), time=time.strftime("%Y-%m-%dT%H:%M:%S"),
               scripts={os.path.basename(x): C.sha256(x) for x in [f"{SCR}/ed7spec.py", f"{SCR}/l4lib.py", f"{SCR}/ed7_01_build_a.py", f"{SCR}/ed7_02_build_b.py", f"{SCR}/ed7_03_compose.py", f"{SCR}/ed7_04_verify_r49.py", f"{SCR}/run_ed7.py", f"{C.LANE}/scripts/stamp_title_r49.py", f"{C.LANE}/scripts/ledb_common.py"]},
               result=open(f"{VER}/checks.txt").read().rstrip().split("\n")[-1]), open(f"{VER}/provenance.json", "w"), indent=1)
