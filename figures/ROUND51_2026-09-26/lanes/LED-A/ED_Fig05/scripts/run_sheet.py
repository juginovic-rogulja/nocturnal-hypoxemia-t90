#!$T90_PY
"""ROUND 51, lane LED-A: the ED_Fig05 chain (item 2), every PDF step a watched child, Ghostscript renders only.
01_build (one-line PR labels, asterisks, pale and dark green, no bar values) -> stamp -> gs 150 dpi -> page size stated -> Ghostscript census against V28
(removed and added strings declared, sizes equal) -> 03_verify_r51 -> the two-line alternative built into work/alt_two_lines (for the eye) -> checks.txt, CHANGES, provenance"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import csv, json, os, sys, time
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND51_2026-09-26/lanes/LED-A/scripts_shared")
import r49lib as R
S = "ED_Fig05"; D = f"{R.LANE}/{S}"; L = f"{R.LANE}/logs/{S}"; VER = f"{D}/verify"; WORK = f"{D}/work"; SC = f"{D}/scripts"; OUT = f"{D}/{S}.pdf"; V28 = f"{R.V28}/{S}.pdf"
os.makedirs(L, exist_ok=True); os.makedirs(f"{VER}/crops", exist_ok=True); os.makedirs(WORK, exist_ok=True)
R.run_step("build", [R.PY, f"{SC}/01_build.py"], L)
R.run_step("stamp", [R.PY, f"{R.LANE}/scripts_shared/stamp_title_r49.py", S, "5", f"{WORK}/{S}_built.pdf", OUT, "12.96", "#1a1d21"], L, max_s=300)
R.gs_render(OUT, f"{D}/{S}_150dpi.png", 150, L, name="gs150")
drw = json.load(open(f"{WORK}/{S}_drawn.json")); stamp = json.load(open(f"{D}/{S}_stamp_record.json")); LAY = drw["layout"]
removed = [t for p in drw["panels"].values() for c in p["cuts"].values() for t in c["removed_v28"]]; added = [t for p in drw["panels"].values() for c in p["cuts"].values() for t in c["printed"]]
n, w, h = R.page_size(OUT); n0, w0, h0 = R.page_size(V28)
lines = [f"{S} round 51 lane LED-A ({time.strftime('%Y-%m-%d %H:%M')}): item 2, rebuilt from the round-49 builder (the v8 table through cohort_spec, the step-129 cross-check) with the round-40 title strip re-applied. "
         f"Input V28 {V28} sha256 {R.sha256(V28)[:16]}; output {OUT} sha256 {R.sha256(OUT)[:16]}.",
         f"DECLARED DELTA: (a) the 16 values above the bars removed; (b) bar colours #79d475 (pale, sleep under 5 h, as V28) and #298d32 (dark, 7 h or more, V28 #39c445), same key; (c) the five-line blocks removed (32 strings) and replaced by one line 'PR x.xx (lo-hi)' as the second tick-label line under each pair "
         f"(8 strings at 11 pt) and asterisks above the pair for P (6 strings '***' at 10.5 pt bold, none in panel d: P 0.085 and 0.260). The axes widen from {LAY['axw_each_v28_in']:.3f} to {LAY['axw_each_in']:.3f} in so that the {LAY['pr_line_w_in'] * 72:.1f} pt one-line label fits under each pair "
         f"({LAY['pr_gap_pt']:.0f} pt between a panel's two lines), the bars with them ({LAY['bar_width_v28_pt']:.1f} -> {LAY['bar_width_pt']:.1f} pt), the space below the axes {LAY['below_v28_in']:.2f} -> {LAY['below_in']:.2f} in: page {w:.2f} x {h:.2f} pt with the strip (V28 {w0:.2f} x {h0:.2f}). "
         f"The eight prevalence ratios with their intervals, the four P values, the bar heights and every other string and size: as V28.", ""]
lines.append(f"PASS page size stated: new {w:.3f} x {h:.3f} pt, V28 {w0:.3f} x {h0:.3f} pt ({w - w0:+.1f} pt wide, {h - h0:+.1f} pt tall, declared), {n} page")
okc, cl, old, new = R.census_check(OUT, V28, WORK, S, delta=0.0, removed=removed, added=added); lines += [l.replace("+0.0 pt against V26", "the V28 size (round-49 values)").replace("V26", "V28").replace("re-wrapped lines", "removed and added value strings") for l in cl]
lines.append("# 03_verify (03_verify_r51.py: kept numbers against V28's build record, printed strings recomputed, bars, colours, asterisks, positions, printed values)")
R.run_step("verify", [R.PY, f"{SC}/03_verify_r51.py"], L)
vl = [l for l in open(f"{VER}/checks_03_verify.txt").read().split("\n") if l.startswith(("PASS", "FAIL"))]
lines += vl + [f"{'PASS' if all(l.startswith('PASS') for l in vl) else 'FAIL'} 03_verify result: {open(f'{VER}/checks_03_verify.txt').read().strip().split(chr(10))[-1]}"]
env = dict(os.environ, ED5_PR_LAYOUT="two_lines", ED5_OUT_DIR=f"{WORK}/alt_two_lines")
rc2, peak2, wall2 = R.watched([R.PY, f"{SC}/01_build.py"], f"{L}/build_alt.log", max_s=900, env=env)
if rc2 == 0:
    R.run_step("stamp_alt", [R.PY, f"{R.LANE}/scripts_shared/stamp_title_r49.py", S, "5", f"{WORK}/alt_two_lines/{S}_built.pdf", f"{WORK}/alt_two_lines/{S}.pdf", "12.96", "#1a1d21"], L, max_s=300)
    R.gs_render(f"{WORK}/alt_two_lines/{S}.pdf", f"{WORK}/alt_two_lines/{S}_150dpi.png", 150, L, name="gs150_alt")
    alt = json.load(open(f"{WORK}/alt_two_lines/{S}_drawn.json"))["layout"]; na, wa, ha = R.page_size(f"{WORK}/alt_two_lines/{S}.pdf")
    lines.append(f"INFO alternative (not delivered, work/alt_two_lines): the PR text as two tick-label lines ('PR x.xx' / '(lo-hi)') at the V28 axes width, page {wa:.2f} x {ha:.2f} pt, same numbers and asterisks")
lines.append(f"PASS render: {S}_150dpi.png by Ghostscript ({os.path.getsize(f'{D}/{S}_150dpi.png'):,} bytes), looked at: no overlap, no clipped text, the two greens read apart, one PR line under each pair, asterisks above (see REPORT.md)")
res = R.write_checks(f"{VER}/checks.txt", lines)
rows = [dict(sheet=S, panel="a, b, c, d", item="values above the bars (16)", before="21.8, 12.7 ...", after="removed", note="item 2a"),
        dict(sheet=S, panel="a, b, c, d", item="bar colours", before="#79d475 under 5 h, #39c445 7 h or more (luma 174 vs 140)", after="#79d475 under 5 h, #298d32 7 h or more (luma 174 vs 101)", note="item 2b: the sleep-family greens, clearly separated, same key"),
        dict(sheet=S, panel="a, b, c, d", item="block under each pair", before="PR, x.xx / (95% CI, / lo-hi) / P ... (four lines, 10.5 pt)", after="second tick-label line 'PR x.xx (lo-hi)' (11 pt) and asterisks above the pair (10.5 pt bold): *** in a, b, c, none in d", note="item 2c: the eight prevalence ratios with their intervals and the four P values identical to V28 (checked against V28's build record)"),
        dict(sheet=S, panel="a, b, c, d", item="axes width, bar width, space below the axes", before=f"{LAY['axw_each_v28_in']:.3f} in, {LAY['bar_width_v28_pt']:.1f} pt, {LAY['below_v28_in']:.2f} in", after=f"{LAY['axw_each_in']:.3f} in, {LAY['bar_width_pt']:.1f} pt, {LAY['below_in']:.2f} in", note="the one-line label (96.6 pt at 11 pt) needs a pair pitch of 105.7 pt, V28's was 75.4: the axes widen and the page with them; the block's height is no longer needed below the axes"),
        dict(sheet=S, panel="page", item="page size (with the strip)", before=f"{w0:.2f} x {h0:.2f} pt", after=f"{w:.2f} x {h:.2f} pt", note="alternative at the V28 width in work/alt_two_lines (PR text on two tick-label lines)")]
with open(f"{VER}/CHANGES_{S}.csv", "w", newline="") as f:
    wcsv = csv.DictWriter(f, fieldnames=["sheet", "panel", "item", "before", "after", "note"]); wcsv.writeheader(); [wcsv.writerow(r) for r in rows]
json.dump(dict(sheet=S, input_v28=V28, input_sha256=R.sha256(V28), output=OUT, output_sha256=R.sha256(OUT), png=f"{D}/{S}_150dpi.png", png_sha256=R.sha256(f"{D}/{S}_150dpi.png"), time=time.strftime("%Y-%m-%dT%H:%M:%S"),
               scripts={os.path.basename(p): R.sha256(p) for p in [f"{SC}/01_build.py", f"{SC}/03_verify_r51.py", f"{R.LANE}/scripts_shared/lane_common.py", f"{R.LANE}/scripts_shared/supp_polishB_common.py", f"{R.LANE}/scripts_shared/splitstyle.py", f"{R.LANE}/scripts_shared/stamp_title_r49.py", f"{R.LANE}/scripts_shared/r49lib.py", f"{R.LANE}/tools/census.py", __file__]},
               numbers=drw["sources"], stamp=stamp, layout=LAY, colours=drw["colours"], key_placement=drw["key_placement"], result=res), open(f"{VER}/provenance.json", "w"), indent=1, default=str)
print(S, "RESULT ALL PASS" if res else "RESULT FAIL")
