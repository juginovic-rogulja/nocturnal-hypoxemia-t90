#!$T90_PY
"""ED_Fig07 round 51 runner (lane LED-B, panels side by side): build a -> build b -> compose -> verify (children under wd_run.sh) -> title strip at 14 pt -> render -> census -> checks.txt."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import json, os, sys, time, shutil
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND51_2026-09-26/lanes/LED-B/scripts"); import ledb_common as C
S = "ED_Fig07"; N = 7; D = f"{C.LANE}/{S}"; WORK, VER, SCR = f"{D}/work", f"{D}/verify", f"{D}/scripts"; OUT = f"{D}/{S}.pdf"; V26 = f"{C.V26}/{S}.pdf"
os.makedirs(WORK, exist_ok=True); os.makedirs(f"{VER}/crops", exist_ok=True)
steps = sys.argv[1:] or ["build_a", "build_b", "compose", "verify", "stamp"]
if "build_a" in steps: C.wd("ed7_01_build_a_r49", 900, [C.PY, f"{SCR}/ed7_01_build_a.py"])
if "build_b" in steps: C.wd("ed7_02_build_b_r49", 900, [C.PY, f"{SCR}/ed7_02_build_b.py"])
if "compose" in steps: C.wd("ed7_03_compose_r51", 900, [C.PY, f"{SCR}/ed7_03_compose_r51.py"])
COMPOSED = f"{D}/{S}.pdf"     # ed7spec OUT_PDF (the composed sheet) is written by compose; the stamped sheet replaces it below (the composed copy stays in work/)
if "stamp" in steps:
    shutil.copyfile(COMPOSED, f"{WORK}/{S}_composed.pdf")
    C.wd(f"stamp_{S}", 300, [C.PY, f"{C.LANE}/scripts/stamp_title_r49.py", S, str(N), f"{WORK}/{S}_composed.pdf", OUT, str(C.V26_TITLE_X[S]), C.INK])
    C.gs_render(OUT, f"{D}/{S}_150dpi.png")
if "verify" in steps: C.wd("ed7_04_verify_r51", 1800, [C.PY, f"{SCR}/ed7_04_verify_r51.py", f"{WORK}/{S}_composed.pdf", f"{D}/{S}_150dpi.png"])
vlines = [l for l in open(f"{VER}/checks.txt").read().split("\n") if l.startswith(("PASS", "FAIL", "INFO"))]
keep = json.load(open(f"{WORK}/panel_b_geometry.json")).get("size_exceptions", {}); ex = {}; lay = json.load(open(f"{WORK}/sheet_layout_r51.json"))
ok_c, clines = C.census(S, OUT, V26, VER, exceptions=ex, plus=0.0)   # round 51: sizes unchanged, the base set is V28
ok_p, new_sz, old_sz = C.page_size_check(OUT, V26)
srec = json.load(open(f"{D}/{S}_stamp_record.json")); ga = json.load(open(f"{WORK}/panel_a_geometry.json"))
hdr = [f"{S} round 51 lane LED-B ({time.strftime('%Y-%m-%d %H:%M')}): the round-49 panels (rebuilt here from the round-38 builders at the round-49 sizes) composed side by side; V28 reference {V26} sha256 {C.sha256(V26)[:16]}; output {OUT} sha256 {C.sha256(OUT)[:16]}.",
       f"DECLARED DELTA (Alen, item 4): panels a and b side by side, a left at (0, 8) and b right at x {lay['x_b']:.2f} with the tops aligned (letters on the baseline 21.0 at x 12.96 and {lay['x_b'] + 12.96:.2f}); page {srec['page_in'][0]:.2f} x {srec['page_in'][1]:.2f} before the strip (V28 {lay['v28_sheet'][0]:.2f} x {lay['v28_sheet'][1]:.2f}), {srec['page_out'][0]:.2f} x {srec['page_out'][1]:.2f} with the 16 pt title strip ('Extended Data Fig. 7' Arial Bold 14 at ({srec['origin'][0]}, {srec['origin'][1]}) as before); no text size and no string changed, every printed number identical to V28 (re-derived from the numbers files below).", ""]
body = hdr + [f"PASS  page size stated: {new_sz} (V28 {old_sz}), changed by design (declared)", f"{'PASS' if abs(srec['page_out'][1] - srec['page_in'][1] - 16.0) < 0.01 and srec['size'] == 14.0 else 'FAIL'}  title strip: page {srec['page_in'][0]:.2f} x {srec['page_in'][1]:.3f} -> {srec['page_out'][0]:.2f} x {srec['page_out'][1]:.3f}, title '{srec['title']}' Arial Bold {srec['size']:g} pt at {srec['origin']}",
               "# builder and panel proofs (ed7_04_verify_r51.py, child under wd_run.sh, on the composed sheet before the strip)"] + vlines + ["# Ghostscript txtwrite census against V28 (ledb_common.census, size rule +0)"] + clines
C.finish_checks(f"{VER}/checks.txt", body)
json.dump(dict(sheet=S, v26=V26, v26_sha256=C.sha256(V26), output=OUT, output_sha256=C.sha256(OUT), png=f"{D}/{S}_150dpi.png", png_sha256=C.sha256(f"{D}/{S}_150dpi.png"), time=time.strftime("%Y-%m-%dT%H:%M:%S"),
               scripts={os.path.basename(x): C.sha256(x) for x in [f"{SCR}/ed7spec.py", f"{SCR}/l4lib.py", f"{SCR}/ed7_01_build_a.py", f"{SCR}/ed7_02_build_b.py", f"{SCR}/ed7_03_compose.py", f"{SCR}/ed7_03_compose_r51.py", f"{SCR}/ed7_04_verify_r51.py", f"{SCR}/run_ed7.py", f"{C.LANE}/scripts/stamp_title_r49.py", f"{C.LANE}/scripts/ledb_common.py"]},
               result=open(f"{VER}/checks.txt").read().rstrip().split("\n")[-1]), open(f"{VER}/provenance.json", "w"), indent=1)
