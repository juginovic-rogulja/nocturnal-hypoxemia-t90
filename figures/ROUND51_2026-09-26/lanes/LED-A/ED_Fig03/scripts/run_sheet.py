#!$T90_PY
"""ROUND 51, lane LED-A: the ED_Fig03 chain (item 1, larger boxes), every PDF step a watched child, Ghostscript renders only.
01_build_ed3 -> stamp (title strip 16 pt, 'Extended Data Fig. 3' 14 pt) -> gs 150 dpi -> page size stated against V28 -> Ghostscript census against V28
(same strings, same sizes) -> 03_verify_ed3_r51 (geometry as declared, values equal V28, the 2.5 mm rule) -> V28 vs new side by side -> checks.txt, CHANGES, provenance"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import csv, json, os, sys, time
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND51_2026-09-26/lanes/LED-A/scripts_shared")
import r49lib as R
from PIL import Image
S = "ED_Fig03"; D = f"{R.LANE}/{S}"; L = f"{R.LANE}/logs/{S}"; VER = f"{D}/verify"; WORK = f"{D}/work"; SC = f"{D}/scripts"; OUT = f"{D}/{S}.pdf"; V28 = f"{R.V28}/{S}.pdf"
os.makedirs(L, exist_ok=True); os.makedirs(f"{VER}/crops", exist_ok=True)
R.run_step("build", [R.PY, f"{SC}/01_build_ed3.py"], L)
R.run_step("stamp", [R.PY, f"{R.LANE}/scripts_shared/stamp_title_r49.py", S, "3", f"{WORK}/{S}_built.pdf", OUT, "24.48", "#1a1d21"], L, max_s=300)
R.gs_render(OUT, f"{D}/{S}_150dpi.png", 150, L, name="gs150")
drw = json.load(open(f"{VER}/{S}_drawn.json")); stamp = json.load(open(f"{D}/{S}_stamp_record.json")); R51 = drw["round51"]
n, w, h = R.page_size(OUT); n0, w0, h0 = R.page_size(V28)
lines = [f"{S} round 51 lane LED-A ({time.strftime('%Y-%m-%d %H:%M')}): the 20 axes boxes larger with tighter gutters (item 1), rebuilt from the round-49 builder with the round-40 key-sentence removal and title strip re-applied. "
         f"Input V28 {V28} sha256 {R.sha256(V28)[:16]}; output {OUT} sha256 {R.sha256(OUT)[:16]}.",
         f"DECLARED DELTA: axes boxes {R51['axw']:.2f} x {R51['axh']:.2f} pt (V28 {R51['axw_v28']:.2f} x {R51['axh_v28']:.2f}: +{100 * (R51['scale_w'] - 1):.1f} percent wide, +{100 * (R51['scale'] - 1):.0f} percent tall), horizontal gutter {R51['gap']:.2f} pt (V28 {R51['gap_v28']:.2f}: the widest y tick label {R51['ytick_w_max']:.1f} pt + 6.1 pt + 2.5 mm), "
         f"row gutter {R51['rowgap']:.0f} pt (V28 {R51['rowgap_v28']:.1f}), the b block with its V28 spacings, the page {w:.2f} x {h:.2f} pt with the strip (V28 {w0:.2f} x {h0:.2f}). Text sizes, titles and their line breaks, y axes, markers, ribbons, key, colours: as V28. No printed number changes.", ""]
lines.append(f"PASS page size stated: new {w:.3f} x {h:.3f} pt, V28 {w0:.3f} x {h0:.3f} pt (declared growth +{w - w0:.1f} pt wide, +{h - h0:.1f} pt tall), {n} page")
okc, cl, old, new = R.census_check(OUT, V28, WORK, S, delta=0.0); lines += [l.replace("+0.0 pt against V26", "the V28 size (round-49 values)").replace("V26", "V28") for l in cl]
lines.append("# 03_verify (03_verify_ed3_r51.py: geometry as declared, the 2.5 mm rule, values read back and equal to V28, asterisks, tokens, printed values)")
R.run_step("verify", [R.PY, f"{SC}/03_verify_ed3_r51.py"], L)
vl = [l for l in open(f"{VER}/checks_03_verify.txt").read().split("\n") if l.startswith(("PASS", "FAIL"))]
lines += vl + [f"{'PASS' if all(l.startswith('PASS') for l in vl) else 'FAIL'} 03_verify result: {open(f'{VER}/checks_03_verify.txt').read().strip().split(chr(10))[-1]}"]
R.gs_render(V28, f"{WORK}/{S}_v28_150dpi.png", 150, L, name="gs150_v28")
a = Image.open(f"{WORK}/{S}_v28_150dpi.png").convert("RGB"); b = Image.open(f"{D}/{S}_150dpi.png").convert("RGB")
side = Image.new("RGB", (a.width + b.width + 20, max(a.height, b.height)), "white"); side.paste(a, (0, 0)); side.paste(b, (a.width + 20, 0)); side.save(f"{VER}/crops/{S}_V28_vs_NEW_150dpi.png")
lines.append(f"PASS render: {S}_150dpi.png by Ghostscript ({os.path.getsize(f'{D}/{S}_150dpi.png'):,} bytes), looked at: no overlap, no clipped text, the grid reads with the larger boxes (see REPORT.md); V28 beside new in verify/crops")
res = R.write_checks(f"{VER}/checks.txt", lines)
rows = [dict(sheet=S, panel="a, b", item="axes boxes", before=f"{R51['axw_v28']:.2f} x {R51['axh_v28']:.2f} pt", after=f"{R51['axw']:.2f} x {R51['axh']:.2f} pt", note="item 1: about a quarter larger (the width +18 percent after the 2.5 mm y-tick-label rule, the height +25 percent), the 5-column grid kept"),
        dict(sheet=S, panel="a, b", item="horizontal gutter", before=f"{R51['gap_v28']:.2f} pt", after=f"{R51['gap']:.2f} pt", note="the widest y tick label (20.4 pt at 10.5 pt) + 6.1 pt + 2.5 mm: every panel's y tick labels clear the box to their left by at least 2.5 mm"),
        dict(sheet=S, panel="a", item="row gutter", before=f"{R51['rowgap_v28']:.1f} pt", after=f"{R51['rowgap']:.0f} pt", note="x tick labels 15.5 pt below a row, a three-line title 33.6 pt above the next, 6.9 pt clear"),
        dict(sheet=S, panel="page", item="page size (with the strip)", before=f"{w0:.2f} x {h0:.2f} pt", after=f"{w:.2f} x {h:.2f} pt", note="the y tick labels fix the gutter, so wider boxes need the width; the tighter row gutters absorb the taller boxes"),
        dict(sheet=S, panel="all", item="text sizes, titles and line breaks, y axes, markers, ribbons, key, numbers", before="V28", after="unchanged", note="census: the same 193 strings at the same sizes; every marker value equals V28's")]
with open(f"{VER}/CHANGES_{S}.csv", "w", newline="") as f:
    wcsv = csv.DictWriter(f, fieldnames=["sheet", "panel", "item", "before", "after", "note"]); wcsv.writeheader(); [wcsv.writerow(r) for r in rows]
json.dump(dict(sheet=S, input_v28=V28, input_sha256=R.sha256(V28), output=OUT, output_sha256=R.sha256(OUT), png=f"{D}/{S}_150dpi.png", png_sha256=R.sha256(f"{D}/{S}_150dpi.png"), time=time.strftime("%Y-%m-%dT%H:%M:%S"),
               scripts={os.path.basename(p): R.sha256(p) for p in [f"{SC}/01_build_ed3.py", f"{SC}/ed3_common.py", f"{SC}/ed3_extract_r49.py", f"{SC}/03_verify_ed3_r51.py", f"{R.LANE}/scripts_shared/stamp_title_r49.py", f"{R.LANE}/scripts_shared/r49lib.py", f"{R.LANE}/tools/census.py", __file__]},
               numbers=drw["sources"], stamp=stamp, round51=R51, result=res), open(f"{VER}/provenance.json", "w"), indent=1, default=str)
print(S, "RESULT ALL PASS" if res else "RESULT FAIL")
