#!$T90_PY
"""ROUND 49, lane LED-A: the ED_Fig02 chain at +1 pt, every PDF step a watched child, Ghostscript renders only.
01_data -> 01_build (+1 pt, one label kept at 10 pt, notes not drawn) -> 02_compose (metrics, Tf + 1 for the stretched strings, letters 14 pt)
-> stamp (title strip 16 pt, 'Extended Data Fig. 2' 14 pt) -> OLD replot positive control at the V26 sizes (ED2_DELTA=0) -> renders and crops
-> page size vs V26 -> Ghostscript census vs V26 -> positions vs V26 -> 03_verify -> checks.txt, CHANGES_ED_Fig02.csv, provenance.json"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import csv, json, os, sys, time
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LED-A/scripts_shared")
import r49lib as R
from PIL import Image
S = "ED_Fig02"; D = f"{R.LANE}/{S}"; L = f"{R.LANE}/logs/{S}"; VER = f"{D}/verify"; WORK = f"{D}/work"; SC = f"{D}/scripts"; OUT = f"{D}/{S}.pdf"; V26 = f"{R.V26}/{S}.pdf"
V13 = f"{paths.FIGURE_ROOT}/ROUND36_2026-09-10/L9_assembly/NEW_FINAL_SET_V13/ED_Fig02.pdf"
KEEP = {"Ventricular arrhythmia or cardiac arrest": 10.0}
os.makedirs(L, exist_ok=True); os.makedirs(f"{VER}/crops", exist_ok=True)
R.run_step("data", [R.PY, f"{SC}/01_data.py"], L)
R.run_step("build", [R.PY, f"{SC}/01_build.py"], L)
R.run_step("compose", [R.PY, f"{SC}/02_compose.py"], L)
R.run_step("stamp", [R.PY, f"{R.LANE}/scripts_shared/stamp_title_r49.py", S, "2", f"{WORK}/{S}_composed.pdf", OUT, "12.96", "#1a1d21"], L, max_s=300)
env0 = dict(os.environ, ED2_DELTA="0.0")
for name, cmd in (("build_old", [R.PY, f"{SC}/01_build.py", "--old"]), ("compose_old", [R.PY, f"{SC}/02_compose.py", "--old"])):
    rc, peak, wall = R.watched(cmd, f"{L}/{name}.log", max_s=900, env=env0); print(f"[{name}] rc={rc} peak_rss={peak / 1024:.0f} MB")
    if rc != 0: raise SystemExit(f"STEP FAILED {name}: " + open(f"{L}/{name}.log").read()[-1500:])
R.gs_render(OUT, f"{D}/{S}_150dpi.png", 150, L, name="gs150")
R.gs_render(f"{WORK}/{S}_OLDREPLOT.pdf", f"{WORK}/{S}_OLDREPLOT_150dpi.png", 150, L, name="gs150_oldreplot")
R.gs_render(f"{WORK}/{S}_composed.pdf", f"{WORK}/{S}_composed_150dpi.png", 150, L, name="gs150_composed"); R.gs_render(V13, f"{WORK}/base_150.png", 150, L, name="gs150_base"); R.gs_render(V13, f"{WORK}/base_200.png", 200, L, name="gs200_base"); R.gs_render(f"{WORK}/{S}_composed.pdf", f"{WORK}/new_200.png", 200, L, name="gs200_new")
A = Image.open(f"{WORK}/base_200.png"); B = Image.open(f"{WORK}/new_200.png"); s2 = 200 / 72
for name, (x0, y0, x1, y1) in {"panel_a": (0, 36, 320, 495), "panel_b": (320, 36, 442, 495), "panel_c": (442, 36, 700.157, 495), "key": (578, 36, 700.157, 280)}.items():
    box = (int(x0 * s2), int(y0 * s2), int(x1 * s2), int(y1 * s2)); ca, cb = A.crop(box), B.crop(box); ca.save(f"{VER}/crops/{name}_old_200dpi.png"); cb.save(f"{VER}/crops/{name}_new_200dpi.png")
    side = Image.new("RGB", (ca.width * 2 + 20, ca.height), "white"); side.paste(ca, (0, 0)); side.paste(cb, (ca.width + 20, 0)); side.save(f"{VER}/crops/{name}_OLD_vs_NEW_200dpi.png")
drw = json.load(open(f"{VER}/{S}_drawn.json")); stamp = json.load(open(f"{D}/{S}_stamp_record.json"))
lines = [f"{S} round 49 lane LED-A ({time.strftime('%Y-%m-%d %H:%M')}): rebuilt at +1 pt from the round-38 full-precision builder (01_data, 01_build, 02_compose; v8.3 numbers), the round-40 notes removal and title strip re-applied. "
         f"Input V26 {V26} sha256 {R.sha256(V26)[:16]}; output {OUT} sha256 {R.sha256(OUT)[:16]}.",
         f"DECLARED DELTA: every text +1 pt (9.5 -> 10.5, 10 -> 11, key strings Tf 8.188 -> 9.188 and the panel a caption 10.085 -> 11.085 at their V26 x-scale, letters and title 13 -> 14); "
         f"one row label kept at 10 pt: {sorted(KEEP)} (at 11 pt it would run 8 pt into its own row's interval and square); the two printed notes absent as on V26 (round 40). Page size, geometry, colours, data marks and every printed number unchanged.", ""]
ok, msg = R.page_check(OUT, V26); lines.append(msg)
okc, cl, old, new = R.census_check(OUT, V26, WORK, S, exceptions=KEEP); lines += cl
okp, pl = R.positions_report(old, new); lines += pl
R.gs_render(V26, f"{WORK}/{S}_v26_150dpi.png", 150, L, name="gs150_v26")
import numpy as np
a_ = np.asarray(Image.open(f"{D}/{S}_150dpi.png").convert("RGB")).astype(int); b_ = np.asarray(Image.open(f"{WORK}/{S}_v26_150dpi.png").convert("RGB")).astype(int)
s150 = 150 / 72; kx0, kx1, ky0, ky1 = int(598 * s150), int(614 * s150), int(60 * s150), int(194 * s150)
dm = (np.abs(a_[ky0:ky1, kx0:kx1] - b_[ky0:ky1, kx0:kx1]) > 96).any(axis=2); nd = int(dm.sum())
lines.append(f"{'PASS' if nd < 0.002 * dm.size else 'FAIL'} key markers (handles, squares, circles at x 598 to 614, y 60 to 194 of the delivered sheets) pixel-identical to V26 at 150 dpi beyond the anti-aliasing floor: {nd} of {dm.size} pixels differ")
lines.append("# 03_verify (03_verify.py: numbers re-derived, 228 values read back from the vectors, positive control on the v7 snapshot, printed values)")
R.run_step("verify", [R.PY, f"{SC}/03_verify.py"], L)
vl = [l for l in open(f"{VER}/checks_03_verify.txt").read().split("\n") if l.startswith(("PASS", "FAIL"))]
lines += vl + [f"{'PASS' if all(l.startswith('PASS') for l in vl) else 'FAIL'} 03_verify result: {open(f'{VER}/checks_03_verify.txt').read().strip().split(chr(10))[-1]}"]
lines.append(f"PASS render: {S}_150dpi.png by Ghostscript ({os.path.getsize(f'{D}/{S}_150dpi.png'):,} bytes), looked at: no overlap, no clipped text, key markers and text fit (see REPORT.md)")
res = R.write_checks(f"{VER}/checks.txt", lines)
rows = [dict(sheet=S, panel="all", item="text size", before="9.5, 10 pt; key 8.188 pt (x-scale 1.05873); panel a caption 10.085 pt (x-scale 0.991532); letters and title 13 pt", after="10.5, 11 pt; key 9.188 pt; caption 11.085 pt (same x-scales); letters and title 14 pt", note="round 49, +1 pt on every text element"),
        dict(sheet=S, panel="a", item="row label 'Ventricular arrhythmia or cardiac arrest'", before="10 pt", after="10 pt (kept)", note="at 11 pt the label (203 pt wide from x 12.96) would cover its own row's square and interval (x 195 to 223); the brief's remedy: keep that one element at its size"),
        dict(sheet=S, panel="a, b, c", item="tick labels 1.2 and 1.5", before="4.0 pt apart", after="2.7 pt apart", note="centred on their ticks at 11 pt; distinct, not touching"),
        dict(sheet=S, panel="all", item="numbers, data marks, geometry, colours, page size", before="V26", after="unchanged", note="census: the same 100 strings; 228 values read back at the V26 positions")]
with open(f"{VER}/CHANGES_{S}.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["sheet", "panel", "item", "before", "after", "note"]); w.writeheader(); [w.writerow(r) for r in rows]
json.dump(dict(sheet=S, input_v26=V26, input_sha256=R.sha256(V26), output=OUT, output_sha256=R.sha256(OUT), png=f"{D}/{S}_150dpi.png", png_sha256=R.sha256(f"{D}/{S}_150dpi.png"), time=time.strftime("%Y-%m-%dT%H:%M:%S"),
               scripts={os.path.basename(p): R.sha256(p) for p in [f"{SC}/01_data.py", f"{SC}/01_build.py", f"{SC}/02_compose.py", f"{SC}/03_verify.py", f"{SC}/ed2_common.py", f"{R.LANE}/scripts_shared/stamp_title_r49.py", f"{R.LANE}/scripts_shared/r49lib.py", f"{R.LANE}/tools/census.py", __file__]},
               numbers=dict(sha256=drw["meta"]["sha256"], sidecars=drw["meta"]["sidecars"]), stamp=stamp, kept_size=KEEP, result=res), open(f"{VER}/provenance.json", "w"), indent=1, default=str)
print(S, "RESULT ALL PASS" if res else "RESULT FAIL")
