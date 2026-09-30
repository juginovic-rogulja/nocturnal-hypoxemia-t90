#!$T90_PY
"""ED_Fig07 round 52 runner (lane LED-B, the reference key on two lines): build a -> build b -> compose -> verify (children under wd_run.sh) -> title strip at 14 pt -> render -> census -> checks.txt."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import json, os, sys, time, shutil
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND52_2026-09-26/lanes/LED-B/scripts"); import ledb_common as C
S = "ED_Fig07"; N = 7; D = f"{C.LANE}/{S}"; WORK, VER, SCR = f"{D}/work", f"{D}/verify", f"{D}/scripts"; OUT = f"{D}/{S}.pdf"; V26 = f"{C.V26}/{S}.pdf"
os.makedirs(WORK, exist_ok=True); os.makedirs(f"{VER}/crops", exist_ok=True)
steps = sys.argv[1:] or ["build_a", "build_b", "compose", "verify", "stamp"]
if "build_a" in steps: C.wd("ed7_01_build_a_r49", 900, [C.PY, f"{SCR}/ed7_01_build_a.py"])
if "build_b" in steps: C.wd("ed7_02_build_b_r49", 900, [C.PY, f"{SCR}/ed7_02_build_b.py"])
if "compose" in steps: C.wd("ed7_03_compose_r52", 900, [C.PY, f"{SCR}/ed7_03_compose_r52.py"])
COMPOSED = f"{D}/{S}.pdf"     # ed7spec OUT_PDF (the composed sheet) is written by compose; the stamped sheet replaces it below (the composed copy stays in work/)
if "stamp" in steps:
    shutil.copyfile(COMPOSED, f"{WORK}/{S}_composed.pdf")
    C.wd(f"stamp_{S}", 300, [C.PY, f"{C.LANE}/scripts/stamp_title_r49.py", S, str(N), f"{WORK}/{S}_composed.pdf", OUT, str(C.V26_TITLE_X[S]), C.INK])
    C.gs_render(OUT, f"{D}/{S}_150dpi.png")
if "verify" in steps: C.wd("ed7_04_verify_r52", 1800, [C.PY, f"{SCR}/ed7_04_verify_r52.py", f"{WORK}/{S}_composed.pdf", f"{D}/{S}_150dpi.png"])
vlines = [l for l in open(f"{VER}/checks.txt").read().split("\n") if l.startswith(("PASS", "FAIL", "INFO"))]
gb = json.load(open(f"{WORK}/panel_b_geometry.json")); keep = gb.get("size_exceptions", {}); lay = json.load(open(f"{WORK}/sheet_layout_r51.json")); key = gb["key"]
ex = dict(removed=[dict(text="Reference, no or mild apnea with T90 ≤1%", count=1, why="the reference key set on two lines (Alen, round 52), wording unchanged")], added=[dict(text=t, count=1, why="the reference key set on two lines (Alen, round 52), wording unchanged") for t in key["ref_lines"]])
ok_c, clines = C.census(S, OUT, V26, VER, exceptions=ex, plus=0.0)   # round 52: sizes unchanged, the base set is V29
ok_p, new_sz, old_sz = C.page_size_check(OUT, V26)
srec = json.load(open(f"{D}/{S}_stamp_record.json")); ga = json.load(open(f"{WORK}/panel_a_geometry.json"))
hdr = [f"{S} round 52 lane LED-B ({time.strftime('%Y-%m-%d %H:%M')}): the round-51 build (panels side by side) with panel b's reference key on two lines; V29 reference {V26} sha256 {C.sha256(V26)[:16]}; output {OUT} sha256 {C.sha256(OUT)[:16]}.",
       f"DECLARED DELTA (Alen, round 52 item 1): the key line '{'Reference, no or mild apnea with T90 ≤1%'}' (right edge {lay['x_b'] + key['ref_right_v29']:.2f} on the page, past the grid's right edge {lay['x_b'] + key['grid_right']:.2f}) is set on two lines {key['ref_lines']} at the same size ({key['key_fs'] + 1:g} pt), same left edge (x {lay['x_b'] + key['ref_text_x']:.2f}), pitch {key['key_pitch']} pt, the reference box beside the first line, the key's right edge now {lay['x_b'] + key['key_right']:.2f}; nothing else moved (ramp key, panel a, letters, title strip as in round 51), page size unchanged ({srec['page_out'][0]:.2f} x {srec['page_out'][1]:.2f}), no text size changed, every printed number identical to V29.", ""]
body = hdr + [f"{'PASS' if ok_p else 'FAIL'}  page size equals V29 (0.05 pt): {new_sz} vs {old_sz}", f"{'PASS' if abs(srec['page_out'][1] - srec['page_in'][1] - 16.0) < 0.01 and srec['size'] == 14.0 else 'FAIL'}  title strip: page {srec['page_in'][0]:.2f} x {srec['page_in'][1]:.3f} -> {srec['page_out'][0]:.2f} x {srec['page_out'][1]:.3f}, title '{srec['title']}' Arial Bold {srec['size']:g} pt at {srec['origin']}",
               "# builder and panel proofs (ed7_04_verify_r52.py, child under wd_run.sh, on the composed sheet before the strip)"] + vlines + ["# Ghostscript txtwrite census against V29 (ledb_common.census, size rule +0, the split key declared)"] + clines
C.finish_checks(f"{VER}/checks.txt", body)
json.dump(dict(sheet=S, v26=V26, v26_sha256=C.sha256(V26), output=OUT, output_sha256=C.sha256(OUT), png=f"{D}/{S}_150dpi.png", png_sha256=C.sha256(f"{D}/{S}_150dpi.png"), time=time.strftime("%Y-%m-%dT%H:%M:%S"),
               scripts={os.path.basename(x): C.sha256(x) for x in [f"{SCR}/ed7spec.py", f"{SCR}/l4lib.py", f"{SCR}/ed7_01_build_a.py", f"{SCR}/ed7_02_build_b.py", f"{SCR}/ed7_03_compose.py", f"{SCR}/ed7_03_compose_r52.py", f"{SCR}/ed7_04_verify_r52.py", f"{SCR}/run_ed7.py", f"{C.LANE}/scripts/stamp_title_r49.py", f"{C.LANE}/scripts/ledb_common.py"]},
               result=open(f"{VER}/checks.txt").read().rstrip().split("\n")[-1]), open(f"{VER}/provenance.json", "w"), indent=1)
