#!$T90_PY
"""ED_Fig06 round 51 runner (lane LED-B): close-up edit (child) -> Ghostscript renders -> verify (child) -> census against V28 (sizes unchanged, 26 strings declared removed) -> checks.txt -> provenance."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import json, os, sys, time
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND51_2026-09-26/lanes/LED-B/scripts"); import ledb_common as C
S = "ED_Fig06"; D = f"{C.LANE}/{S}"; WORK, VER, SCR = f"{D}/work", f"{D}/verify", f"{D}/scripts"; OUT = f"{D}/{S}.pdf"; V28 = f"{C.V28}/{S}.pdf"; DUMP = f"{C.LANE}/verify/v28_spans/{S}_spans.json"
os.makedirs(WORK, exist_ok=True); os.makedirs(VER, exist_ok=True)
assert json.load(open(DUMP))["sha256"] == C.sha256(V28), "the V28 sheet changed since the dump"
C.wd(f"r51_closeup_{S}", 600, [C.PY, f"{SCR}/01_r51_close_up_ed6.py", V28, WORK, OUT])
C.gs_render(V28, f"{WORK}/{S}_V28_150dpi.png"); C.gs_render(OUT, f"{D}/{S}_150dpi.png")
try: C.wd(f"r51_verify_{S}", 900, [C.PY, f"{SCR}/02_r51_verify_ed6.py", OUT, DUMP, f"{WORK}/r51_record.json", VER, f"{WORK}/{S}_V28_150dpi.png", f"{D}/{S}_150dpi.png"])
except SystemExit as e: print("verify child failed:", e)
vlines = open(f"{VER}/verify_lines.txt").read().rstrip("\n").split("\n") if os.path.exists(f"{VER}/verify_lines.txt") else ["FAIL  verify child did not write its lines"]
rec = json.load(open(f"{WORK}/r51_record.json"))
from collections import Counter
ex = dict(removed=[dict(text=t, count=n, why="the 'New diagnoses / group size' column of panel b, removed (Alen, round 51 item 3)") for t, n in sorted(Counter(r["text"] for r in rec["removed"]).items())])
ok_c, clines = C.census(S, OUT, V28, VER, exceptions=ex, plus=0.0)
ok_p, new_sz, old_sz = C.page_size_check(OUT, V28)
hdr = [f"{S} round 51 lane LED-B ({time.strftime('%Y-%m-%d %H:%M')}): input V28 {V28} sha256 {C.sha256(V28)[:16]}; output {OUT} sha256 {C.sha256(OUT)[:16]}.",
       f"DECLARED DELTA (Alen, item 3, as the coordinator specified): panel b's 'New diagnoses / group size' column removed (24 events/group-size strings and the two header lines); the freed width given to the forest axis: plot box x0 {rec['box_v28'][0]:.2f} -> {rec['box_new'][0]:.2f} (right edge, y range and value range unchanged), ticks {[round(x, 1) for x in rec['axis_v28']['ticks_x']]} -> {[round(x, 1) for x in rec['axis_new']['ticks_x']]}, the 24 intervals and markers re-plotted from results_crossed.csv on the new axis, the four tick labels re-centred; the HR (95% CI) column, the key, the row labels and everything else at their V28 positions; page size unchanged; text sizes unchanged (round-49 values); every kept number identical to V28.", ""]
body = hdr + [f"{'PASS' if ok_p else 'FAIL'}  page size: {new_sz} (V28 {old_sz}), unchanged", "# text-layer, line-art, raster and value proofs (02_r51_verify_ed6.py, child under wd_run.sh)"] + vlines + ["# Ghostscript txtwrite census against V28 (ledb_common.census, size rule +0)"] + clines
C.finish_checks(f"{VER}/checks.txt", body)
json.dump(dict(sheet=S, round=51, input=V28, input_sha256=C.sha256(V28), output=OUT, output_sha256=C.sha256(OUT), png=f"{D}/{S}_150dpi.png", png_sha256=C.sha256(f"{D}/{S}_150dpi.png"), time=time.strftime("%Y-%m-%dT%H:%M:%S"),
               scripts={os.path.basename(x): C.sha256(x) for x in [f"{SCR}/01_r51_close_up_ed6.py", f"{SCR}/02_r51_verify_ed6.py", f"{SCR}/run_r51_ed6.py", f"{C.LANE}/scripts/ledb_common.py"]}, result=open(f"{VER}/checks.txt").read().rstrip().split("\n")[-1]), open(f"{VER}/provenance.json", "w"), indent=1)
