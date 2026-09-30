#!$T90_PY
"""Supp Fig 5, round 49 (+1 pt) verification against V26. Page box (pypdf); Ghostscript text census against V26 (same strings, every
size +1.0, the declared exceptions); read-back of every printed value re-derived here from the numbers files (the round-38 verifier's
derivation, matched against the census); the vector layer of the four flat parts against the round-38 parts (bracket and legend
handles excluded); masked raster diff of the composed sheets at 150 dpi (differences confined to text boxes and legend handles);
the 150 dpi render. Writes verify/checks.txt (ends RESULT ALL PASS), verify/census_*.json, verify/crops/*.png."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import json, os, re, sys
from collections import Counter
import numpy as np
L = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LSUPP-B"
sys.path.insert(0, f"{L}/scripts"); sys.path.insert(0, f"{L}/Supp_Fig05/scripts")
import r49_common as C
from l3b_common import NUM, SV, sidecar, jload, hydrated  # noqa: E402
import v14lib  # noqa: E402
S = "Supp_Fig05"; D = f"{L}/{S}"; W = f"{D}/work"; VER = f"{D}/verify"; os.makedirs(f"{VER}/crops", exist_ok=True)
OLD = f"{C.V26}/{S}.pdf"; NEW = f"{D}/{S}.pdf"
R38 = f"{paths.FIGURE_ROOT}/ROUND38_2026-09-16/figures/lanes/Supp_Fig05/work"
E = json.load(open(f"{W}/expected_values.json"))
ck = C.Checks(f"{S} round 49 (+{C.PT_PLUS} pt) verification against V26, {C.now()}")
# ---------------------------------------------------------------- 1 page box
n_new, w_new, h_new = C.page_box(NEW); n_old, w_old, h_old = C.page_box(OLD)
ck.log("one page, page box equals V26 (0.05 pt)", n_new == 1 and abs(w_new - w_old) < 0.05 and abs(h_new - h_old) < 0.05, f"V26 {w_old} x {h_old}, new {w_new} x {h_new}")
# ---------------------------------------------------------------- 2 text census
old = C.merge_kerned(C.gs_text(OLD, f"{W}/{S}_V26_txtwrite.xml")); new = C.merge_kerned(C.gs_text(NEW, f"{W}/{S}_NEW_txtwrite.xml"))
json.dump({"old": old, "new": new}, open(f"{VER}/census_spans.json", "w"), indent=0)
KEPT_C = E["geometry"]["c"]["kept_labels"]; KEPT_D = E["geometry"]["d"]["kept_labels"]
exceptions = {}
for lab in KEPT_C:
    exceptions[(lab, "ArialMT", 10.0)] = 10.0
for lab in KEPT_D:   # two-line labels: each line is one census string
    for line in [s["text"] for s in old if s["font"] == "ArialMT" and s["size"] == 10.0 and s["x"] > 530 and 520 < s["y"] < 660]:
        if line.replace("\n", " ") in lab or lab.startswith(line) or lab.endswith(line):
            exceptions[(line, "ArialMT", 10.0)] = 10.0
ok, det = C.census_compare(old, new, exceptions)
ck.log("census: the same multiset of strings as V26", det["strings_equal"], f"V26 {det['n_old']} spans, new {det['n_new']} (kerning splits merged)")
ck.log(f"census: every string one point larger than on V26 except the {len(exceptions)} declared strings kept at 10 pt", ok, f"missing {det['missing']} extra {det['extra']}")
ck.info(f"declared exceptions (kept at the current size, the wrap would change at +1 pt): {sorted(k[0] for k in exceptions)}")
ck.info(f"size table V26 {C.size_table(old)}; new {C.size_table(new)}")
letters_old = {s["text"]: (s["x"], s["y"]) for s in old if s["size"] == 13.0}; letters_new = {s["text"]: (s["x"], s["y"]) for s in new if s["size"] == 14.0}
ck.log("letters a to d at 14 pt at the V26 origins (1 pt, the census is integer)", set(letters_new) == {"a", "b", "c", "d"} and all(abs(letters_old[k][0] - letters_new[k][0]) <= 1 and abs(letters_old[k][1] - letters_new[k][1]) <= 1 for k in letters_old), f"{letters_old} -> {letters_new}")
fonts = {s["font"] for s in new}
ck.log("fonts are the Arial family only", fonts <= {"ArialMT", "Arial-BoldMT", "Arial-ItalicMT", "Arial Bold"}, str(fonts))
ck.log("no note, no italic face on the new sheet (round 40 LNOTES kept)", not any(s["font"] == "Arial-ItalicMT" for s in new) and not any("preserved nocturnal" in s["text"] for s in new))
# left-anchored strings keep their x (row labels, letters, column heads, legend text): report the largest shift per class
def cls(s):
    if s["size"] in (13.0, 14.0): return "letter"
    if s["font"] == "Arial-BoldMT" or re.fullmatch(r"\d+\.\d\d \(\d+\.\d\d-\d+\.\d\d\)|<0\.001|\d\.\d{3}|q|HR \(95% CI\)", s["text"]): return "column"
    if s["size"] in (11.0, 12.0): return "axis title"
    if re.fullmatch(r"[\d.]+", s["text"]): return "tick"
    return "label"
shift = {}
for s in old:
    cands = [n for n in new if n["text"] == s["text"] and abs(n["y"] - s["y"]) <= 3]
    if cands:
        n = min(cands, key=lambda n: abs(n["x"] - s["x"])); k = cls(s); shift[k] = max(shift.get(k, 0), abs(n["x"] - s["x"]))
ck.info(f"largest x shift by class (pt, integer census): {shift} (ticks and axis titles are centred, the q column was nudged 1.44 pt)")
# ---------------------------------------------------------------- 3 read-back from the numbers files (independent of the builder)
OE_PATH, LJ_PATH = f"{NUM}/one_exposure_v1.json", f"{SV}/reorg_fits_efig4_lag2/joint_t90_tst_landmark.json"
OE = jload(OE_PATH); sha_oe, sc_oe = sidecar(OE_PATH); LJ = jload(LJ_PATH); sha_lj, sc_lj = sidecar(LJ_PATH)
ck.log("sidecar sha256 of the sources unchanged since the build, both v8.1", sha_oe == E["sources"]["one_exposure_v1.json"]["sha256"] and sha_lj == E["sources"]["joint_t90_tst_landmark.json"]["sha256"], f"steps {sc_oe['step']['id']} ({sc_oe['output']['mtime_local']}), {sc_lj['step']['id']} ({sc_lj['output']['mtime_local']})")
EX = OE["exposures"]
WHOLE_T90 = "Oxygen below 90% for more than 10% of the night"; WHOLE_SHORT = "Sleeping fewer than 5 hours"
AHI_ALL = "Severe sleep apnea, apnea-hypopnea index 30 or more"; AHI_NORM = "Severe sleep apnea, in normal oxygenation only"
SE_NORM = "Sleep efficiency below 85%, in normal oxygenation only"; N3_NORM = "N3 sleep below 5% of the night, in normal oxygenation only"
SHORT_NORM = "Sleeping fewer than 5 hours, in normal oxygenation only"
ORDER = [WHOLE_T90, AHI_ALL, AHI_NORM, N3_NORM, SHORT_NORM, SE_NORM, WHOLE_SHORT]
def bh(pmap):
    names = sorted(pmap, key=lambda k: pmap[k]); m, q, prev = len(names), {}, 1.0
    for r in range(m, 0, -1):
        prev = min(prev, pmap[names[r - 1]] * m / r); q[names[r - 1]] = prev
    return q
common = set(EX[WHOLE_T90]["outcomes"])
for e in (WHOLE_T90, WHOLE_SHORT, AHI_ALL, AHI_NORM): common &= set(EX[e]["outcomes"])
common = {c for c in common if not EX[WHOLE_T90]["outcomes"][c]["circular"] and not EX[WHOLE_T90]["outcomes"][c]["negative_control"]}
conds = sorted(common, key=lambda c: -EX[WHOLE_T90]["outcomes"][c]["hr"])[:14]
Q = {e: bh({c: v["p"] for c, v in EX[e]["outcomes"].items() if not v["circular"] and not v["negative_control"]}) for e in EX}
printed = Counter()
for prim in (WHOLE_T90, AHI_NORM):
    for c in conds:
        o = EX[prim]["outcomes"][c]; printed[v14lib.hr_ci(o["hr"], o["lo"], o["hi"])] += 1; printed["<0.001" if Q[prim][c] < 0.001 else v14lib.halfup(Q[prim][c], 3)] += 1
counts_c = [sum(1 for c in conds if EX[e]["outcomes"][c]["lo"] > 1) for e in ORDER]
NEWD = {}
for tag in ("t90", "tst"):
    for lag in ("0", "2"):
        hrs = [LJ["outcomes"][k][lag][f"{tag}_hr"] for k in LJ["outcomes"]]; qs = [LJ["outcomes"][k][lag][f"{tag}_q"] for k in LJ["outcomes"]]
        NEWD[(tag, lag)] = sum(1 for h, q in zip(hrs, qs) if q < 0.05 and h > 1)
ck.log("verifier re-derives the builder's row set (top 14 by the T90 hazard ratio)", conds == E["row_set"], str(conds))
exp_strings = Counter(printed)
exp_strings[f"Severe sleep apnea, preserved oxygenation (n exposed = {EX[AHI_NORM]['n_exposed']:,})"] += 1
exp_strings[f"All severe sleep apnea (n exposed = {EX[AHI_ALL]['n_exposed']:,})"] += 1
exp_strings[f"Conditions with a raised rate (FDR q < 0.05), of {len(LJ['outcomes'])}"] += 1
for c in conds: exp_strings[c] += 2
new_strings = Counter(s["text"] for s in new)
missing = {k: v for k, v in exp_strings.items() if new_strings.get(k, 0) < v}
ck.log("read-back: every value re-derived from the numbers files is printed on the sheet the expected number of times (28 HR strings, 28 q strings, the two legend n, the panel d axis end, 14 rows twice)", not missing, f"missing {missing}")
ck.log("panel c counts and panel d counts re-derived equal the builder's record", counts_c == E["counts_c"] and [NEWD[k] for k in (("t90", "0"), ("t90", "2"), ("tst", "0"), ("tst", "2"))] == E["counts_d"], f"c {counts_c} d {[NEWD[k] for k in NEWD]}")
builder_texts = Counter(v["text"] for v in E["values"])
ck.log("builder's expected values all printed and within the independent re-derivation", all(new_strings.get(t, 0) >= n for t, n in builder_texts.items()) and set(builder_texts) <= set(exp_strings), str([t for t in builder_texts if t not in exp_strings]))
# ---------------------------------------------------------------- 4 the vector layer of the flat parts against the round-38 parts
dr_new = {p: C.flat_drawings(f"{W}/part_{p}.pdf") for p in "abcd"}; dr_old = {p: C.flat_drawings(f"{R38}/part_{p}.pdf") for p in "abcd"}
allowed = {"a": [], "b": [], "c": [], "d": []}
for p in "ab":   # legend handles (two rows under the axis title): the marker and line handles move with the legend's font size
    hs = [d for d in dr_old[p] + dr_new[p] if d["rect"][1] > (E["part_heights_pt"][p] - 0.8 * 72) and d["rect"][0] < 60]
    if hs: allowed[p].append([min(d["rect"][0] for d in hs) - 1, min(d["rect"][1] for d in hs) - 1, max(d["rect"][2] for d in hs) + 1, max(d["rect"][3] for d in hs) + 1])
brk = [d for d in dr_old["c"] if d["stroke"] == "#298d32" and d["width"] and abs(d["width"] - 0.8) < 0.05]   # the round-38 bracket (removed in round 40)
if brk: allowed["c"].append([min(d["rect"][0] for d in brk) - 1, min(d["rect"][1] for d in brk) - 1, max(d["rect"][2] for d in brk) + 1, max(d["rect"][3] for d in brk) + 1])
for p in "abcd":
    ok_d, det_d = C.drawings_compare(dr_old[p], dr_new[p], allowed[p])
    ck.log(f"part {p}: vector layer identical to the round-38 part (every path rect, fill, stroke, width) outside {len(allowed[p])} excluded box(es) ({'legend handles' if p in 'ab' else 'the removed bracket' if p == 'c' else 'none'})", ok_d, f"old {det_d['n_old']} new {det_d['n_new']} paths, excluded {det_d['excluded_old']}/{det_d['excluded_new']}; only old {det_d['only_old'][:3]} only new {det_d['only_new'][:3]}")
ck.log("part c: no bracket strokes on the new part (round 40 LNOTES kept)", not any(d["stroke"] == "#298d32" and d["width"] and abs(d["width"] - 0.8) < 0.05 for d in dr_new["c"]))
exp_marks = {}
M = E["markers"]
def nfill(exp, filled): return sum(1 for m in M if m["exposure"] == exp and m["filled"] == filled)
cnt = Counter((d["fill"], d["stroke"], d["kinds"]) for p in "ab" for d in dr_new[p])
exp_marks = {("#0288d1", "#ffffff", "c"): nfill(WHOLE_T90, True) + 1, ("#ffffff", "#0288d1", "c"): nfill(WHOLE_T90, False),
             ("#39c445", "#ffffff", "re"): nfill(WHOLE_SHORT, True) + 1, ("#ffffff", "#298d32", "re"): nfill(WHOLE_SHORT, False),
             ("#d55e00", "#ffffff", "c"): nfill(AHI_NORM, True) + 1, ("#ffffff", "#d55e00", "c"): nfill(AHI_NORM, False),
             ("#eeb08f", "#ffffff", "re"): nfill(AHI_ALL, True) + 1, ("#ffffff", "#eeb08f", "re"): nfill(AHI_ALL, False)}
bad = {k: (cnt.get(k, 0), v) for k, v in exp_marks.items() if cnt.get(k, 0) != v}
ck.log("marker fills (filled = q < 0.05) and series colours match the drawing layer of parts a and b, legend handles included", not bad, str(bad))
per_ux_c, per_ux_d = E["geometry"]["c"]["per_ux_in"] * 72, E["geometry"]["d"]["per_ux_in"] * 72
bars_c = [d for d in dr_new["c"] if d["kinds"] == "re" and d["stroke"] is None and d["fill"] in ("#0288d1", "#eeb08f", "#d55e00", "#39c445", "#d5f1d0")]
bars_d = [d for d in dr_new["d"] if d["kinds"] == "re" and d["stroke"] is None and d["fill"] in ("#0288d1", "#eeb08f", "#d55e00", "#39c445", "#d5f1d0")]
c_read = sorted(round((d["rect"][2] - d["rect"][0]) / per_ux_c, 2) for d in bars_c); d_read = sorted(round((d["rect"][2] - d["rect"][0]) / per_ux_d, 2) for d in bars_d)
ck.log("part c bars read back to the counts at the V13 scale", c_read == sorted(float(c) for c in counts_c if c > 0), f"read {c_read} counts {counts_c}")
ck.log("part d bars read back to the counts at the V13 scale", d_read == sorted(float(NEWD[k]) for k in NEWD if NEWD[k] > 0), f"read {d_read} counts {[NEWD[k] for k in NEWD]}")
ck.log("part c band present", any(d["fill"] == "#eefaee" for d in dr_new["c"]))
cols_old = {c for p in "abcd" for d in dr_old[p] for c in (d["fill"], d["stroke"]) if c}; cols_new = {c for p in "abcd" for d in dr_new[p] for c in (d["fill"], d["stroke"]) if c}
ck.log("no colour introduced (drawing colours within the round-38 parts')", cols_new <= cols_old, f"new-only {sorted(cols_new - cols_old)}")
# ---------------------------------------------------------------- 5 renders and the masked raster diff
C.gs_render(NEW, f"{D}/{S}_150dpi.png", 150); C.gs_render(OLD, f"{W}/{S}_V26_150dpi.png", 150)
boxes = [C.text_box(s, pad=2.5) for s in old + new]
PLACE = E["place"]
for p in "abcd":
    for b in allowed[p]:
        boxes.append([PLACE[p][0] + b[0], PLACE[p][1] + b[1], PLACE[p][0] + b[2], PLACE[p][1] + b[3]])
rd = C.raster_diff(f"{W}/{S}_V26_150dpi.png", f"{D}/{S}_150dpi.png", boxes, 150)
ck.log("150 dpi raster: every differing pixel lies inside a text box (V26 or new census) or a legend-handle box; data marks, axes and bands identical", rd["shape_equal"] and rd["n_outside"] == 0, str(rd))
for name, box in (("a_legend_foot", [0, 380, 535, 470]), ("b_columns", [850, 40, 1021, 200]), ("c_labels", [0, 500, 400, 800]), ("d_labels", [520, 500, 940, 720])):
    C.crop_png(f"{W}/{S}_V26_150dpi.png", 150, box, f"{VER}/crops/{S}_{name}_OLD.png"); C.crop_png(f"{D}/{S}_150dpi.png", 150, box, f"{VER}/crops/{S}_{name}_NEW.png")
ck.info("crops OLD (V26) and NEW in verify/crops/: a_legend_foot, b_columns, c_labels, d_labels")
ck.info(f"new sheet sha256 {C.sha256(NEW)}; V26 sha256 {C.sha256(OLD)}")
json.dump({"exceptions": [list(k) + [v] for k, v in exceptions.items()], "allowed_boxes_part": allowed, "census": det, "raster": rd, "x_shift": shift}, open(f"{VER}/{S}_verify_record.json", "w"), indent=1)
ok_all = ck.write(f"{VER}/checks.txt"); print("RESULT ALL PASS" if ok_all else "RESULT FAIL"); sys.exit(0 if ok_all else 1)
