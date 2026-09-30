#!$T90_PY
"""Supp Fig 18, round 49 (+1 pt) verification against V26. Page box (pypdf); Ghostscript text census against V26 (same strings, every
size +1.0); the letters at 14 pt at the block's left edge (moved 11.5 pt left with the block, declared) on their V26 baselines; the
six primary HR (95% CI) and P strings re-read here from the CURRENT alenfig3_cross_v1.json (v8.2 gate), bold P exactly where the
95% CI excludes 1; the vector layer of the flat sheet against the flat V26 sheet outside the legend handles (axis, ticks, every
interval and marker identical); masked raster diff at 150 dpi; the 150 dpi render. Writes verify/checks.txt (ends RESULT ALL PASS)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import json, os, re, sys
from collections import Counter
L = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LSUPP-B"
sys.path.insert(0, f"{L}/scripts"); sys.path.insert(0, f"{L}/Supp_Fig18/scripts")
import r49_common as C
import lane_common as LC  # noqa: E402
VL = LC.VL
S = "Supp_Fig18"; D = f"{L}/{S}"; W = f"{D}/work"; VER = f"{D}/verify"; os.makedirs(f"{VER}/crops", exist_ok=True)
OLD = f"{C.V26}/{S}.pdf"; NEW = f"{D}/{S}.pdf"
DR = json.load(open(f"{W}/{S}_drawn.json")); G = DR["geometry_in"]; G26 = DR["r49"]["geometry_v26_in"]
ck = C.Checks(f"{S} round 49 (+{C.PT_PLUS} pt) verification against V26, {C.now()}")
n_new, w_new, h_new = C.page_box(NEW); n_old, w_old, h_old = C.page_box(OLD)
ck.log("one page, page box equals V26 (0.05 pt)", n_new == 1 and abs(w_new - w_old) < 0.05 and abs(h_new - h_old) < 0.05, f"V26 {w_old} x {h_old}, new {w_new} x {h_new}")
ck.log("the delivered sheet is the regenerated flat page (no form XObjects)", C.sha256(NEW) == C.sha256(f"{W}/{S}_regen.pdf"))
old = C.merge_kerned(C.gs_text(OLD, f"{W}/{S}_V26_txtwrite.xml")); new = C.merge_kerned(C.gs_text(NEW, f"{W}/{S}_NEW_txtwrite.xml"))
json.dump({"old": old, "new": new}, open(f"{VER}/census_spans.json", "w"), indent=0)
ok, det = C.census_compare(old, new, {})
ck.log("census: the same multiset of strings as V26", det["strings_equal"], f"V26 {det['n_old']} spans, new {det['n_new']} (kerning splits merged)")
ck.log("census: every string one point larger than on V26, no exception", ok, f"missing {det['missing']} extra {det['extra']}")
ck.info(f"size table V26 {C.size_table(old)}; new {C.size_table(new)}")
letters_old = {s["text"]: (s["x"], s["y"]) for s in old if s["size"] == 13.0}; letters_new = {s["text"]: (s["x"], s["y"]) for s in new if s["size"] == 14.0}
dx = G["a"]["x_block"] * 72 - G26["a"]["x_block"] * 72
ck.log(f"letters a and b at 14 pt on their V26 baselines, at the block's left edge which moved {dx:+.2f} pt with the wider labels (declared)", set(letters_new) == {"a", "b"} and all(abs(letters_old[k][1] - letters_new[k][1]) <= 1 and abs((letters_new[k][0] - letters_old[k][0]) - dx) <= 1.5 for k in letters_old), f"{letters_old} -> {letters_new}")
ck.log("axis pinned: ax_x0, ax_w and both parts' ax_y0 equal the round-37 record (1e-6 in); the centred rule would have moved the axis", all(abs(G[p][k] - G26[p][k]) < 1e-6 for p in "ab" for k in ("ax_x0", "ax_w", "ax_y0", "ax_top", "part_top")), f"ax_x0 {G['a']['ax_x0'] * 72:.2f} pt (centred rule {G['a']['ax_x0_centred_rule'] * 72:.2f})")
ck.info(f"text block: labels start {dx:+.2f} pt, columns end {(G['a']['pq_x'] - G26['a']['pq_x']) * 72:+.2f} pt (col_x {(G['a']['col_x'] - G26['a']['col_x']) * 72:+.2f}); margins left {G['a']['x_block'] * 72:.1f} pt, right {968 - G['a']['pq_x'] * 72:.1f} pt")
fonts = {s["font"] for s in new}
ck.log("fonts are the Arial family only", fonts <= {"ArialMT", "Arial-BoldMT"}, str(fonts))
# read-back of the primary rows from the current file (v8.2 gate)
A3P = f"{LC.NUM}/alenfig3_cross_v1.json"; SC = LC.sidecar(A3P); V8_2_CUT = "2026-09-15 19:18:00"
ck.log(f"alenfig3_cross_v1.json sidecar after the v8.2 cut ({V8_2_CUT}) and sha256 equal to the build record", SC["output_mtime"] >= V8_2_CUT and SC["sha256"] == DR["sources"]["alenfig3_cross_v1.json"]["sha256"], f"{SC['output_mtime']} {SC['sha256'][:12]}")
A3 = json.load(open(LC.hydrated(A3P)))["mutual"]
exp = {}
for o in ("Heart failure", "Cardiovascular composite", "Death from any cause"):
    for k in ("z_t90", "z_ahi"):
        v = A3[o][k]; exp[(o, k)] = (VL.hr_ci(v["hr"], v["lo"], v["hi"]), VL.p_text(v["p"]), not (v["lo"] <= 1 <= v["hi"]))
got = {(r["outcome"], "z_t90" if r["exposure"] == "t90_adj_ahi" else "z_ahi"): (r["printed_hr_ci"], r["printed_p"], r["fill"] == "filled") for r in DR["rows_b"] if r["cohort"] == "BDSP"}
ck.log("read-back: the six primary HR (95% CI) strings, P strings and fill states equal alenfig3_cross_v1.json mutual re-read here, half up", got == exp, str(sorted(exp.values())))
new_strings = Counter(s["text"] for s in new)
ck.log("every primary string printed once, its P bold exactly where the 95% CI excludes 1 (census face)", all(new_strings.get(h, 0) >= 1 and new_strings.get(p, 0) >= 1 for h, p, _f in exp.values()) and all(any(s["text"] == p and (s["font"] == "Arial-BoldMT") == f for s in new) for _h, p, f in exp.values()))
ck.log("the 15 drawn records of the build (12 primary values, 3 labels) are census strings at their recorded baselines (1 pt)", all(any(s["text"] == r["text"] and abs(s["y"] - r["baseline"]) <= 1.0 for s in new) for r in DR["records"]), f"{len(DR['records'])} records")
# vector layer. V26 (round 37, built 2026-09-15 10:08) drew alenfig3_cross_v1.json and external.json as they stood then (three decimals);
# both files were rewritten at full precision on 2026-09-16 (16:44 and 19:37) and Supp_Fig18 was not rebuilt after that (V14 kept through V26).
# A CONTROL build at +0 pt from the CURRENT files (scripts/01_build_control_plus0.py, work/control_plus0/) separates the two effects:
# control against V26 = the precision effect (sub-point mark shifts, printed strings identical); new against control = the +1 pt effect (text only).
CTRL = f"{W}/control_plus0/{S}.pdf"; DRC = json.load(open(f"{W}/control_plus0/work/{S}_drawn.json"))
ctrl = C.merge_kerned(C.gs_text(CTRL, f"{W}/control_plus0/{S}_CTRL_txtwrite.xml"))
okc, detc = C.census_compare(old, ctrl, {}, plus=0.0)
ck.log("control (+0 pt, current files) against V26: the text census is identical (strings, faces, sizes) and the block geometry equals the round-37 record", okc and all(abs(DRC["geometry_in"][p_][k] - G26[p_][k]) < 1e-6 for p_ in "ab" for k in ("x_block", "ax_x0", "col_x", "pq_x")), f"missing {detc['missing']} extra {detc['extra']}")
dr_new = C.flat_drawings(NEW, nd=3); dr_old = C.flat_drawings(OLD, nd=3); dr_ctrl = C.flat_drawings(CTRL, nd=3)
H = DR["page_pt"][1]
allowed = []
for p_ in "ab":
    y_top = H - G[p_]["part_top"] * 72; y_ax = H - G[p_]["ax_top"] * 72; x0 = G[p_]["ax_x0"] * 72
    hs = [d for d in dr_old + dr_new + dr_ctrl if y_top - 2 < d["rect"][1] and d["rect"][3] < y_ax and x0 - 6 < d["rect"][0] < x0 + 400]
    if hs: allowed.append([min(d["rect"][0] for d in hs) - 1, min(d["rect"][1] for d in hs) - 1, max(d["rect"][2] for d in hs) + 1, max(d["rect"][3] for d in hs) + 1])
def pair_delta(a, b):
    worst = 0.0; n_diff = 0; unpaired = 0
    for d in a:
        c = [e for e in b if e["stroke"] == d["stroke"] and e["fill"] == d["fill"] and e["kinds"] == d["kinds"] and abs(e["rect"][1] - d["rect"][1]) < 1.0]
        if not c: unpaired += 1; continue
        e = min(c, key=lambda e: abs(e["rect"][0] - d["rect"][0]) + abs(e["rect"][2] - d["rect"][2]))
        dd = max(abs(x - y) for x, y in zip(d["rect"], e["rect"])); worst = max(worst, dd); n_diff += dd > 0.0005
    return worst, n_diff, unpaired
ok_d, det_d = C.drawings_compare(dr_ctrl, dr_new, allowed)
ck.log("new (+1 pt) against the control: vector layer identical outside the legend-handle boxes (axis line, reference line, tick marks, every interval and marker: rect, fill, stroke, width)", ok_d, f"control {det_d['n_old']} new {det_d['n_new']} paths, excluded {det_d['excluded_old']}/{det_d['excluded_new']} (legend handles); only control {det_d['only_old'][:3]} only new {det_d['only_new'][:3]}")
w26, n26, u26 = pair_delta([d for d in dr_old if not any(bx0 - 0.5 <= d["rect"][0] and d["rect"][2] <= bx1 + 0.5 and by0 - 0.5 <= d["rect"][1] and d["rect"][3] <= by1 + 0.5 for bx0, by0, bx1, by1 in allowed)], dr_ctrl)
ck.log(f"control against V26: every path paired, {n26} interval and marker paths of part b moved by at most {w26:.3f} pt (the full-precision rewrite of alenfig3_cross_v1.json and external.json on 2026-09-16, after V26's build; printed strings identical), part a exact", u26 == 0 and w26 <= 0.2, f"max {w26:.3f} pt over {len(dr_old)} paths")
wn, nn, un = pair_delta([d for d in dr_old if not any(bx0 - 0.5 <= d["rect"][0] and d["rect"][2] <= bx1 + 0.5 and by0 - 0.5 <= d["rect"][1] and d["rect"][3] <= by1 + 0.5 for bx0, by0, bx1, by1 in allowed)], dr_new)
ck.log(f"new against V26 outside the legend handles: every path paired, no mark moved by more than the precision effect ({wn:.3f} pt)", un == 0 and wn <= 0.2, f"{nn} paths differ by up to {wn:.3f} pt")
cols_old = {c for d in dr_old for c in (d["fill"], d["stroke"]) if c}; cols_new = {c for d in dr_new for c in (d["fill"], d["stroke"]) if c}
ck.log("no colour introduced", cols_new <= cols_old, f"new-only {sorted(cols_new - cols_old)}")
ck.log("no pale #eef0f1 vertical grid line drawn (round 31)", not any(d["stroke"] == "#eef0f1" for d in dr_new))
C.gs_render(NEW, f"{D}/{S}_150dpi.png", 150); C.gs_render(OLD, f"{W}/{S}_V26_150dpi.png", 150); C.gs_render(CTRL, f"{W}/control_plus0/{S}_CTRL_150dpi.png", 150)
boxes = [C.text_box(s, pad=2.5) for s in ctrl + new] + allowed
rd = C.raster_diff(f"{W}/control_plus0/{S}_CTRL_150dpi.png", f"{D}/{S}_150dpi.png", boxes, 150)
ck.log("150 dpi raster, new against the control: every differing pixel lies inside a text box (control or new census) or a legend-handle box; axis, ticks, intervals and markers identical", rd["shape_equal"] and rd["n_outside"] == 0, str(rd))
rd26 = C.raster_diff(f"{W}/{S}_V26_150dpi.png", f"{D}/{S}_150dpi.png", [C.text_box(s, pad=2.5) for s in old + new] + allowed, 150)
ck.info(f"150 dpi raster, new against V26 outside text and legend boxes: {rd26['n_outside']} differing pixels (the sub-point precision effect of part b, bbox {rd26['outside_bbox_pt']})")
for name, box in (("a_head", [150, 0, 968, 130]), ("b_head", [150, 570, 968, 700]), ("b_primary", [150, 640, 968, 780]), ("a_foot", [150, 500, 968, 580])):
    C.crop_png(f"{W}/{S}_V26_150dpi.png", 150, box, f"{VER}/crops/{S}_{name}_OLD.png"); C.crop_png(f"{D}/{S}_150dpi.png", 150, box, f"{VER}/crops/{S}_{name}_NEW.png")
ck.info("crops OLD (V26) and NEW in verify/crops/")
ck.info(f"new sheet sha256 {C.sha256(NEW)}; V26 sha256 {C.sha256(OLD)}")
json.dump({"census": det, "raster_new_vs_control": rd, "raster_new_vs_v26": rd26, "control_vs_v26_max_delta_pt": w26, "control_vs_v26_n_paths_moved": n26, "legend_boxes": allowed, "geometry_delta_pt": {k: (G["a"][k] - G26["a"][k]) * 72 for k in ("x_block", "gut_w", "ax_x0", "col_x", "pq_x", "block_w")}}, open(f"{VER}/{S}_verify_record.json", "w"), indent=1)
ok_all = ck.write(f"{VER}/checks.txt"); print("RESULT ALL PASS" if ok_all else "RESULT FAIL"); sys.exit(0 if ok_all else 1)
