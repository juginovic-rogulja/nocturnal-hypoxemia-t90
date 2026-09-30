#!$T90_PY
"""ED_Fig06 round 49 runner (lane LED-B): reset (child) -> Ghostscript renders -> verify (child) -> census -> checks.txt -> provenance."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import json, os, sys, time, subprocess
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LED-B/scripts"); import ledb_common as C
S = "ED_Fig06"; D = f"{C.LANE}/{S}"; WORK, VER, SCR = f"{D}/work", f"{D}/verify", f"{D}/scripts"; OUT = f"{D}/{S}.pdf"; V26 = f"{C.V26}/{S}.pdf"; DUMP = f"{C.LANE}/verify/v26_spans/{S}_spans.json"
os.makedirs(WORK, exist_ok=True); os.makedirs(VER, exist_ok=True)
assert json.load(open(DUMP))["sha256"] == C.sha256(V26), "the V26 sheet changed since the dump"
C.wd(f"reset_{S}", 600, [C.PY, f"{SCR}/01_reset_text_ed6.py", V26, WORK, OUT])
C.gs_render(V26, f"{WORK}/{S}_V26_150dpi.png"); C.gs_render(f"{WORK}/{S}_NOTEXT.pdf", f"{WORK}/{S}_NOTEXT_150dpi.png"); C.gs_render(OUT, f"{D}/{S}_150dpi.png")
try: C.wd(f"verify_{S}", 900, [C.PY, f"{SCR}/02_verify_ed6.py", OUT, DUMP, f"{WORK}/reset_record.json", WORK, VER, f"{WORK}/{S}_V26_150dpi.png", f"{WORK}/{S}_NOTEXT_150dpi.png", f"{D}/{S}_150dpi.png"])
except SystemExit as e: print("verify child failed:", e)
vlines = open(f"{VER}/verify_lines.txt").read().rstrip("\n").split("\n") if os.path.exists(f"{VER}/verify_lines.txt") else ["FAIL  verify child did not write its lines"]
# census: the 55 clipped stale copies of panel a's text (Ghostscript lists them on V26, PyMuPDF and the render do not) are declared removed
from collections import Counter
old_gs = C.parse_txtwrite(C.gs_txtwrite(V26, f"{VER}/{S}_V26_txtwrite.xml")); vis = Counter(s["text"].replace("\xa0", " ").replace("\xad", "-") for s in json.load(open(DUMP))["spans"])
hidden = Counter(s["text"] for s in old_gs) - vis
ex = dict(removed=[dict(text=t, size=None, count=n, why="clipped stale copy of panel a's text on V26 (invisible: proved by the V26 -> NOTEXT raster identity), removed by the full-page text redaction") for t, n in sorted(hidden.items())])
ok_c, clines = C.census(S, OUT, V26, VER, exceptions=ex)
ok_p, new_sz, old_sz = C.page_size_check(OUT, V26)
rec = json.load(open(f"{WORK}/reset_record.json"))
hdr = [f"{S} round 49 lane LED-B ({time.strftime('%Y-%m-%d %H:%M')}): input V26 {V26} sha256 {C.sha256(V26)[:16]}; output {OUT} sha256 {C.sha256(OUT)[:16]}.",
       "DECLARED DELTA: every visible string re-set at its V26 size + 1.0 pt (title 13 -> 14, letters 13 -> 14, 11 -> 12, 10 -> 11, 9.5 -> 10.5) in the system Arial faces, same colour and weight, same baseline; anchors: centred strings on their V26 centre, right-aligned on their V26 right edge, left-aligned at their V26 x.",
       f"  nudges (text only, declared): panel b row labels x 27.36 -> {rec['nudges']['row_labels_x'][1]}; events/group-size column and its two header lines right-aligned at {rec['nudges']['events_right_edge'][1]} (was 240.4); smallest label-to-events gap on a row {rec['min_gap_pt']} pt; the column ends {rec['leftmost_plot_ink_x'] - rec['nudges']['events_right_edge'][1]:.1f} pt left of the plot's leftmost visible ink.",
       f"  removed: the {sum(hidden.values())} clipped stale copies of panel a's text that Ghostscript listed on V26 (invisible); no visible string removed, none added; numbers unchanged (values re-derived from the data files below).", ""]
body = hdr + [f"{'PASS' if ok_p else 'FAIL'}  page size equals V26 (0.05 pt): {new_sz} vs {old_sz}", "# text-layer, line-art, raster and value proofs (02_verify_ed6.py, child under wd_run.sh)"] + vlines + ["# Ghostscript txtwrite census against V26 (ledb_common.census)"] + clines
C.finish_checks(f"{VER}/checks.txt", body)
json.dump(dict(sheet=S, input=V26, input_sha256=C.sha256(V26), output=OUT, output_sha256=C.sha256(OUT), png=f"{D}/{S}_150dpi.png", png_sha256=C.sha256(f"{D}/{S}_150dpi.png"), time=time.strftime("%Y-%m-%dT%H:%M:%S"),
               scripts={os.path.basename(x): C.sha256(x) for x in [f"{SCR}/01_reset_text_ed6.py", f"{SCR}/02_verify_ed6.py", f"{SCR}/run_ed6.py", f"{C.LANE}/scripts/ledb_common.py"]}, result=open(f"{VER}/checks.txt").read().rstrip().split("\n")[-1]), open(f"{VER}/provenance.json", "w"), indent=1)
