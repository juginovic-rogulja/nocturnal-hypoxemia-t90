#!$T90_PY
"""ED_Fig08 / ED_Fig09 / ED_Fig10 round 49 runner (lane LED-B): build (child under wd_run.sh, the lane copy of the round-37 / round-40 builder at +1 pt)
-> title strip at 14 pt (child) -> Ghostscript 150 dpi render -> census against V26 -> checks.txt -> provenance.
usage: run_l5a.py <Sheet> [build] [stamp]"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import json, os, sys, time
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LED-B/scripts"); import ledb_common as C
S = sys.argv[1]; N = int(S[-2:]); steps = sys.argv[2:] or ["build", "stamp"]
D = f"{C.LANE}/{S}"; WORK, VER = f"{D}/work", f"{D}/verify"; OUT = f"{D}/{S}.pdf"; V26 = f"{C.V26}/{S}.pdf"; SCR = f"{C.LANE}/L5a_scripts"
BUILDER = {"ED_Fig08": f"{SCR}/08_build_ed8.py", "ED_Fig09": f"{SCR}/09_build_ed9_r49.py", "ED_Fig10": f"{SCR}/10_build_ed10.py"}[S]
BUILT = {"ED_Fig08": f"{WORK}/{S}_built.pdf", "ED_Fig09": f"{WORK}/{S}_rebuilt.pdf", "ED_Fig10": f"{WORK}/{S}_built.pdf"}[S]
BCHK = {"ED_Fig08": f"{WORK}/build_checks.txt", "ED_Fig09": f"{WORK}/build/build_checks.txt", "ED_Fig10": f"{WORK}/build_checks.txt"}[S]
os.makedirs(WORK, exist_ok=True); os.makedirs(f"{VER}/crops", exist_ok=True)
if "build" in steps: C.wd(f"build_{S}_r49", 1200, [C.PY, BUILDER])
if "stamp" in steps:
    C.wd(f"stamp_{S}", 300, [C.PY, f"{C.LANE}/scripts/stamp_title_r49.py", S, str(N), BUILT, OUT, str(C.V26_TITLE_X[S]), C.INK])
    C.gs_render(OUT, f"{D}/{S}_150dpi.png")
blines = [l for l in open(BCHK).read().split("\n") if l.startswith(("PASS", "FAIL", "INFO"))]
blines = [("PASS  " + l[6:]) if l.startswith("PASS: ") else (("FAIL  " + l[6:]) if l.startswith("FAIL: ") else l) for l in blines]
rec = json.load(open(f"{WORK}/r49_record.json"))
ex = dict(sizes=[dict(text=k.split(":", 1)[1].replace("tick ", ""), old=10.0, new=v, why=(f"tick label pair (panel {k.split(':')[0]}) that would come within 2.5 pt of each other at 11 pt keeps its V26 size" if "tick " in k else f"row label (panel {k.split(':')[0]}): the largest size at or above 10 in 0.5 pt steps that stays 4 pt clear of the row's marks and the reference rule")) for k, v in rec.get("size_exceptions", {}).items()],
          removed=[dict(text=t, count=1, why="legend re-wrapped at 10 pt (wording unchanged)") for t in rec.get("legend_rewrap", {}).get("removed", [])],
          added=[dict(text=t, count=1, why="legend re-wrapped at 10 pt (wording unchanged)") for t in rec.get("legend_rewrap", {}).get("added", [])])
ok_c, clines = C.census(S, OUT, V26, VER, exceptions=ex)
ok_p, new_sz, old_sz = C.page_size_check(OUT, V26)
srec = json.load(open(f"{D}/{S}_stamp_record.json"))
hdr = [f"{S} round 49 lane LED-B ({time.strftime('%Y-%m-%d %H:%M')}): rebuilt by the lane copy of {os.path.basename(BUILDER)} at +1 pt on the v8 numbers; V26 reference {V26} sha256 {C.sha256(V26)[:16]}; output {OUT} sha256 {C.sha256(OUT)[:16]}.",
       f"DECLARED DELTA: every text +1 pt (title strip 16 pt with '{srec['title']}' Arial Bold {srec['size']:g} at {srec['origin']}); size exceptions {rec.get('size_exceptions', {})}; text nudges {rec.get('nudges', {})}; not drawn {rec.get('not_drawn', [])}; legend re-wrap {rec.get('legend_rewrap', {})}; numbers unchanged (re-derived from the numbers files by the builder below).", ""]
body = hdr + [f"{'PASS' if ok_p else 'FAIL'}  page size equals V26 (0.05 pt): {new_sz} vs {old_sz}", f"{'PASS' if abs(srec['page_out'][1] - srec['page_in'][1] - 16.0) < 0.01 and srec['size'] == 14.0 else 'FAIL'}  title strip: page {srec['page_in'][0]:.2f} x {srec['page_in'][1]:.3f} -> {srec['page_out'][0]:.2f} x {srec['page_out'][1]:.3f}, title Arial Bold {srec['size']:g} pt at {srec['origin']}",
               f"# builder checks ({os.path.basename(BUILDER)}, child under wd_run.sh: V13 geometry, values re-derived from the sidecar-gated numbers files, read-back, printed-values csv)"] + blines + ["# Ghostscript txtwrite census against V26 (ledb_common.census)"] + clines
C.finish_checks(f"{VER}/checks.txt", body)
json.dump(dict(sheet=S, v26=V26, v26_sha256=C.sha256(V26), built=BUILT, built_sha256=C.sha256(BUILT), output=OUT, output_sha256=C.sha256(OUT), png=f"{D}/{S}_150dpi.png", png_sha256=C.sha256(f"{D}/{S}_150dpi.png"), time=time.strftime("%Y-%m-%dT%H:%M:%S"),
               scripts={os.path.basename(x): C.sha256(x) for x in [BUILDER, f"{SCR}/l5a14.py", f"{SCR}/run_l5a.py", f"{C.LANE}/scripts/stamp_title_r49.py", f"{C.LANE}/scripts/ledb_common.py"]}, result=open(f"{VER}/checks.txt").read().rstrip().split("\n")[-1]), open(f"{VER}/provenance.json", "w"), indent=1)
