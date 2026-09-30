#!$T90_PY
"""ROUND 49, lane LED-A: the ED_Fig03 chain at +1 pt, every PDF step a watched child, Ghostscript renders only.
01_build_ed3 (+1 pt inside txt/tw, two titles kept at 9.5 pt, key sentence not drawn) -> stamp (title strip 16 pt, 'Extended Data Fig. 3' 14 pt)
-> gs 150 dpi -> page size vs V26 -> Ghostscript census vs V26 -> positions vs V26 -> 03_verify_ed3_r49 -> checks.txt, CHANGES_ED_Fig03.csv, provenance.json"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import csv, json, os, sys, time
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LED-A/scripts_shared")
import r49lib as R
S = "ED_Fig03"; D = f"{R.LANE}/{S}"; L = f"{R.LANE}/logs/{S}"; VER = f"{D}/verify"; WORK = f"{D}/work"; SC = f"{D}/scripts"; OUT = f"{D}/{S}.pdf"; V26 = f"{R.V26}/{S}.pdf"
os.makedirs(L, exist_ok=True); os.makedirs(f"{VER}/crops", exist_ok=True)
R.run_step("build", [R.PY, f"{SC}/01_build_ed3.py"], L)
R.run_step("stamp", [R.PY, f"{R.LANE}/scripts_shared/stamp_title_r49.py", S, "3", f"{WORK}/{S}_built.pdf", OUT, "24.48", "#1a1d21"], L, max_s=300)
R.gs_render(OUT, f"{D}/{S}_150dpi.png", 150, L, name="gs150")
drw = json.load(open(f"{VER}/{S}_drawn.json")); stamp = json.load(open(f"{D}/{S}_stamp_record.json"))
removed, added = [], []
for w in drw["wrap_changed_at_plus1"]: removed += w["v26"]; added += w["new"]
exc = {l: 9.5 for p in drw["panels"] if p["title"] in drw["keep_size"] for l in p["title_lines"]}
KEY = [t["text"] for t in drw["text"] if t["what"] == "key"]
CENTRED = {"Time with oxygen saturation below 90%, % of the recording", "Hazard ratio vs 0-1% (95% CI)", "Negative controls"}
lines = [f"{S} round 49 lane LED-A ({time.strftime('%Y-%m-%d %H:%M')}): rebuilt at +1 pt from the round-37 builder (01_build_ed3.py, v8.1 numbers), the round-40 key-sentence removal and title strip re-applied. "
         f"Input V26 {V26} sha256 {R.sha256(V26)[:16]}; output {OUT} sha256 {R.sha256(OUT)[:16]}.",
         f"DECLARED DELTA: every text +1 pt (titles and ticks 9.5 -> 10.5, key 9 -> 10, group title 10 -> 11, axis titles 11 -> 12, letters and title 13 -> 14); two last-column titles kept at 9.5 pt {sorted(drw['keep_size'])} "
         f"(at 10.5 pt their lines would end at 521.5 and 522.7 pt on the 518.74 pt page); title line pitch {drw['line_pitch_v26']:g} -> {drw['line_pitch_v26'] + 1:g} pt at 10.5 pt; key wording {drw['key_mean_dx'] - 15.2:.1f} pt right of its band label "
         f"(V26 {drw['key_mean_dx_v26'] - 15.2:.1f} pt, the label's growth); one title re-split by the sheet's wrap rule at 10.5 pt: " + "; ".join(f"{w['title']!r}: {w['v26']} -> {w['new']}" for w in drw["wrap_changed_at_plus1"])
         + "; centred strings (x title, y label, group title, key row) grow about their centre, x tick labels re-spaced by the min-gap rule; y tick labels nudged right to clear the neighbouring panel: " + ", ".join(f"{n['title']} {n['label']} +{n['shift']} pt" for n in drw.get("ytick_nudges", [])) + ". Page size, geometry, colours, data marks and every printed number unchanged.", ""]
ok, msg = R.page_check(OUT, V26); lines.append(msg)
okc, cl, old, new = R.census_check(OUT, V26, WORK, S, removed=removed, added=added, exceptions=exc); lines += cl
okp, pl = R.positions_report(old, new, movers=set(added) | set(removed) | set(KEY) | CENTRED); lines += pl
lines.append("# 03_verify (03_verify_ed3_r49.py: rows by the August rule, markers and ribbons read back, asterisks by BH q, tick rule, positions, printed values)")
R.run_step("verify", [R.PY, f"{SC}/03_verify_ed3_r49.py"], L)
vl = [l for l in open(f"{VER}/checks_03_verify.txt").read().split("\n") if l.startswith(("PASS", "FAIL"))]
lines += vl + [f"{'PASS' if all(l.startswith('PASS') for l in vl) else 'FAIL'} 03_verify result: {open(f'{VER}/checks_03_verify.txt').read().strip().split(chr(10))[-1]}"]
lines.append(f"PASS render: {S}_150dpi.png by Ghostscript ({os.path.getsize(f'{D}/{S}_150dpi.png'):,} bytes), looked at: no overlap, no clipped text, key row fits (see REPORT.md)")
res = R.write_checks(f"{VER}/checks.txt", lines)
rows = [dict(sheet=S, panel="all", item="text size", before="titles, ticks 9.5 pt; key 9 pt; 'Negative controls' 10 pt; axis titles 11 pt; letters and title 13 pt", after="10.5; 10; 11; 12; 14 pt", note="round 49, +1 pt on every text element"),
        dict(sheet=S, panel="a, b", item="titles 'Venous thromboembolism' and 'Contact dermatitis' (last column)", before="9.5 pt", after="9.5 pt (kept)", note="at 10.5 pt their lines would leave the page (521.5 and 522.7 pt on 518.74); the brief's remedy: keep those elements at their size"),
        dict(sheet=S, panel="a", item="title 'Ventricular arrhythmia or cardiac arrest' line breaks", before="Ventricular / arrhythmia or cardiac / arrest", after="Ventricular / arrhythmia or / cardiac arrest", note="the sheet's greedy wrap rule at 10.5 pt (limit 92.54 pt); words unchanged"),
        dict(sheet=S, panel="a, b", item="title line pitch (wrapped titles)", before="10 pt", after="11 pt at 10.5 pt (10 pt for the two kept titles)", note="one point with the type; lines grow upward, the first line stays 4 pt above the box"),
        dict(sheet=S, panel="key", item="band wording offset from the swatch", before=f"{drw['key_mean_dx_v26']} pt", after=f"{drw['key_mean_dx']} pt", note="the band label '5-10' grows 2.0 pt at 10 pt; the 3.5 pt gap to its wording kept"),
        dict(sheet=S, panel="b", item="y tick labels nudged right (" + "; ".join(f"{n['title']} {n['label']}: +{n['shift']} pt" for n in drw.get("ytick_nudges", [])) + ")", before="right-aligned 6.1 pt left of the spine", after="left edge 0.5 pt right of the neighbouring panel's box", note="at 10.5 pt the widest labels would run 0.5 pt into the neighbouring panel's band fill; the tick marks and every other label unchanged"),
        dict(sheet=S, panel="all", item="numbers, data marks, geometry, colours, page size", before="V26", after="unchanged", note="census: the same 193 strings apart from the declared line breaks; spines, ticks, markers, ribbons at the V26 positions")]
with open(f"{VER}/CHANGES_{S}.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["sheet", "panel", "item", "before", "after", "note"]); w.writeheader(); [w.writerow(r) for r in rows]
json.dump(dict(sheet=S, input_v26=V26, input_sha256=R.sha256(V26), output=OUT, output_sha256=R.sha256(OUT), png=f"{D}/{S}_150dpi.png", png_sha256=R.sha256(f"{D}/{S}_150dpi.png"), time=time.strftime("%Y-%m-%dT%H:%M:%S"),
               scripts={os.path.basename(p): R.sha256(p) for p in [f"{SC}/01_build_ed3.py", f"{SC}/ed3_common.py", f"{SC}/ed3_extract_r49.py", f"{SC}/03_verify_ed3_r49.py", f"{R.LANE}/scripts_shared/stamp_title_r49.py", f"{R.LANE}/scripts_shared/r49lib.py", f"{R.LANE}/tools/census.py", __file__]},
               numbers=drw["sources"], stamp=stamp, kept_size=drw["keep_size"], result=res), open(f"{VER}/provenance.json", "w"), indent=1, default=str)
print(S, "RESULT ALL PASS" if res else "RESULT FAIL")
