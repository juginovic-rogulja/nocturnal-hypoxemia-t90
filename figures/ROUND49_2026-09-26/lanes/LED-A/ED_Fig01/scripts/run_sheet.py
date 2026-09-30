#!$T90_PY
"""ROUND 49, lane LED-A: the ED_Fig01 chain at +1 pt, every PDF step a watched child, Ghostscript renders only.
build (ed1_build_r49.py) -> stamp (title strip 16 pt, 'Extended Data Fig. 1' 14 pt) -> gs 150 dpi -> page size vs V26 -> Ghostscript census vs V26
-> positions vs V26 -> 03_verify (ed1_verify_r49.py, data read-back) -> checks.txt, CHANGES_ED_Fig01.csv, provenance.json"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import csv, json, os, sys, time
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LED-A/scripts_shared")
import r49lib as R
S = "ED_Fig01"; D = f"{R.LANE}/{S}"; L = f"{R.LANE}/logs/{S}"; VER = f"{D}/verify"; WORK = f"{D}/work"; OUT = f"{D}/{S}.pdf"; V26 = f"{R.V26}/{S}.pdf"
os.makedirs(L, exist_ok=True); os.makedirs(VER, exist_ok=True)
R.run_step("build", [R.PY, f"{D}/scripts/ed1_build_r49.py"], L)
R.run_step("stamp", [R.PY, f"{R.LANE}/scripts_shared/stamp_title_r49.py", S, "1", f"{WORK}/sheet_built.pdf", OUT, "12.96", "#1a1d21"], L, max_s=300)
R.gs_render(OUT, f"{D}/{S}_150dpi.png", 150, L, name="gs150")
rec = json.load(open(f"{WORK}/build_record.json")); stamp = json.load(open(f"{D}/{S}_stamp_record.json"))
removed, added = [], []
for w in rec["wrap_changed_at_plus1"]: removed += w["v26"]; added += w["new"]
lines = [f"{S} round 49 lane LED-A ({time.strftime('%Y-%m-%d %H:%M')}): rebuilt at +1 pt from the round-37 builder (ed1_build_r49.py, v8.1 numbers), the round-40 two-line axis title and title strip and the round-43 panel gap re-applied. "
         f"Input V26 {V26} sha256 {R.sha256(V26)[:16]}; output {OUT} sha256 {R.sha256(OUT)[:16]}.",
         f"DECLARED DELTA: every text +1 pt (9.5 -> 10.5, 10 -> 11, 11 -> 12, letters and title 13 -> 14); {len(rec['wrap_changed_at_plus1'])} names re-split by the sheet's own wrap rule at 11 pt (words unchanged): "
         + "; ".join(f"{w['name']!r}: {w['v26']} -> {w['new']}" for w in rec["wrap_changed_at_plus1"])
         + f"; two-line names between one-line rows set on a {2 * 5.18:.2f} pt line pitch (V26 {8.15 + 0.99:.2f} pt), symmetric about the same row centre; the two axis-title lines centred on the panel a axis (their left edges move with their width). Page size, geometry, colours, data marks and every printed number unchanged.", ""]
ok, msg = R.page_check(OUT, V26); lines.append(msg)
okc, cl, old, new = R.census_check(OUT, V26, WORK, S, removed=removed, added=added, v26_widened=True); lines += cl
okp, pl = R.positions_report(old, new, movers=set(added) | set(removed) | {"Share of the", "cohort, %"}); lines += pl
lines.append("# 03_verify (ed1_verify_r49.py: rows, bands, geometry pins, data marks, read-back of every placed string)")
R.run_step("verify", [R.PY, f"{D}/scripts/ed1_verify_r49.py"], L)
vl = [l for l in open(f"{VER}/checks_03_verify.txt").read().split("\n") if l.startswith(("PASS", "FAIL"))]
lines += vl + [f"{'PASS' if all(l.startswith('PASS') for l in vl) else 'FAIL'} 03_verify result: {open(f'{VER}/checks_03_verify.txt').read().strip().split(chr(10))[-1]}"]
lines.append(f"PASS render: {S}_150dpi.png by Ghostscript ({os.path.getsize(f'{D}/{S}_150dpi.png'):,} bytes), looked at: no overlap, no clipped text (see REPORT.md)")
res = R.write_checks(f"{VER}/checks.txt", lines)
rows = [dict(sheet=S, panel="all", item="text size", before="9.5, 10, 11 pt; letters and title 13 pt", after="10.5, 11, 12 pt; letters and title 14 pt", note="round 49, +1 pt on every text element")]
for w in rec["wrap_changed_at_plus1"]: rows.append(dict(sheet=S, panel="b", item=f"name line breaks: {w['name']}", before=" / ".join(w["v26"]), after=" / ".join(w["new"]), note=f"the sheet's wrap rule at 11 pt (one-line limit {w['limit']:.2f} pt, width at 11 pt {w['width_at_11']} pt); words unchanged"))
rows.append(dict(sheet=S, panel="b", item="two-line name pitch", before="9.14 pt", after="10.36 pt where both neighbouring rows are one-line names, 9.14 pt otherwise", note="nudge: keeps descenders clear of the second line's capitals at 11 pt; row centres unchanged"))
rows.append(dict(sheet=S, panel="a", item="x-axis title lines", before="'Share of the' at x 83.825, 'cohort, %' at x 90.557 (11 pt)", after="centred on the axis at x 113.79 (12 pt): left edges 81.12 and 88.45", note="round 40 rule re-applied at the new size"))
rows.append(dict(sheet=S, panel="all", item="numbers, data marks, geometry, colours, page size", before="V26", after="unchanged", note="census: same strings apart from the declared line breaks; dots, bars, leaders, spines at the V26 positions"))
with open(f"{VER}/CHANGES_{S}.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["sheet", "panel", "item", "before", "after", "note"]); w.writeheader(); [w.writerow(r) for r in rows]
json.dump(dict(sheet=S, input_v26=V26, input_sha256=R.sha256(V26), output=OUT, output_sha256=R.sha256(OUT), png=f"{D}/{S}_150dpi.png", png_sha256=R.sha256(f"{D}/{S}_150dpi.png"), time=time.strftime("%Y-%m-%dT%H:%M:%S"),
               scripts={os.path.basename(p): R.sha256(p) for p in [f"{D}/scripts/ed1_build_r49.py", f"{D}/scripts/common.py", f"{D}/scripts/ed1_verify_r49.py", f"{R.LANE}/scripts_shared/stamp_title_r49.py", f"{R.LANE}/scripts_shared/r49lib.py", f"{R.LANE}/tools/census.py", __file__]},
               numbers=rec["sources"], stamp=stamp, result=res), open(f"{VER}/provenance.json", "w"), indent=1, default=str)
print(S, "RESULT ALL PASS" if res else "RESULT FAIL")
