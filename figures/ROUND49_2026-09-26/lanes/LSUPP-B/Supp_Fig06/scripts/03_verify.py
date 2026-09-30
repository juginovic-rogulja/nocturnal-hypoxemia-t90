#!$T90_PY
"""Supp Fig 6, round 49 (+1 pt) verification against V26. Page box (pypdf); Ghostscript text census against V26 (same strings, every
size +1.0, the declared exceptions); independent re-derivation of every printed value from the step-138 and step-139 files (the
round-40 verifier's derivations: panel a cells from the v8 table through cohort_spec, panel b rows and stars, panel c rows, part d
rows and heads) matched against the census; the vector layer of the four flat parts against the round-40 parts; masked raster diff
of the composed sheets at 150 dpi; the 150 dpi render. Writes verify/checks.txt (ends RESULT ALL PASS)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import json, os, re, sys
from collections import Counter
import numpy as np, pandas as pd
L = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LSUPP-B"
sys.path.insert(0, f"{L}/scripts"); sys.path.insert(0, f"{L}/Supp_Fig06/scripts")
import r49_common as C
from l3b_common import SV, T90, V8TAB, sidecar, jload, hydrated  # noqa: E402
S = "Supp_Fig06"; D = f"{L}/{S}"; W = f"{D}/work"; VER = f"{D}/verify"; os.makedirs(f"{VER}/crops", exist_ok=True)
OLD = f"{C.V26}/{S}.pdf"; NEW = f"{D}/{S}.pdf"
R40W = f"{paths.FIGURE_ROOT}/ROUND40_2026-09-18/lanes/LSUPP/Supp_Fig06/work"
E = json.load(open(f"{W}/expected_values.json")); CR = json.load(open(f"{W}/compose_record.json")); CR40 = json.load(open(f"{R40W}/compose_record.json"))
SRC = f"{SV}/dur_x_oxygen_healthy"; V8_2_CUT = "2026-09-15 19:18:00"
ck = C.Checks(f"{S} round 49 (+{C.PT_PLUS} pt) verification against V26, {C.now()}")
# ---------------------------------------------------------------- 1 page box and the compose record
n_new, w_new, h_new = C.page_box(NEW); n_old, w_old, h_old = C.page_box(OLD)
ck.log("one page, page box equals V26 (0.05 pt)", n_new == 1 and abs(w_new - w_old) < 0.05 and abs(h_new - h_old) < 0.05, f"V26 {w_old} x {h_old}, new {w_new} x {h_new}")
ck.log("compose record: part origins, part sizes and the lower-row y equal the round-40 record (the page height follows the parts)", all(CR["parts"][k]["origin"] == CR40["parts"][k]["origin"] and all(abs(a - b) < 0.01 for a, b in zip(CR["parts"][k]["size"], CR40["parts"][k]["size"])) for k in CR["parts"]) and CR["r40_lower_row_y"] == CR40["r40_lower_row_y"] and CR["page"] == CR40["page"], f"{ {k: v['origin'] for k, v in CR['parts'].items()} } lower row y {CR['r40_lower_row_y']}")
ck.log("letters stamped at the round-40 origins", CR["letters_stamped"] == CR40["letters_stamped"], str(CR["letters_stamped"]))
# ---------------------------------------------------------------- 2 text census
old = C.merge_kerned(C.gs_text(OLD, f"{W}/{S}_V26_txtwrite.xml")); new = C.merge_kerned(C.gs_text(NEW, f"{W}/{S}_NEW_txtwrite.xml"))
json.dump({"old": old, "new": new}, open(f"{VER}/census_spans.json", "w"), indent=0)
K = E["r49_kept"]
exceptions = {(lab, "ArialMT", 10.0): 10.0 for lab in set(K["gutter_labels_abc"]) | set(K["gutter_labels_d"])}
for ln, sz in K["cell_lines"]:
    exceptions[(ln, "ArialMT", 9.5)] = sz
ok, det = C.census_compare(old, new, exceptions)
ck.log("census: the same multiset of strings as V26", det["strings_equal"], f"V26 {det['n_old']} spans, new {det['n_new']} (kerning splits merged)")
ck.log(f"census: every string one point larger than on V26 except the declared exceptions ({len(exceptions)} strings)", ok, f"missing {det['missing']} extra {det['extra']}")
ck.info(f"declared exceptions: gutter labels kept at 10 pt (the gutter rule fails at 11 pt): {sorted(set(K['gutter_labels_abc']) | set(K['gutter_labels_d']))}; cell rate lines at 10.0 pt (+0.5, the cell rule fails at 10.5): {sorted(set(l for l, _s in K['cell_lines']))}")
ck.info(f"size table V26 {C.size_table(old)}; new {C.size_table(new)}")
letters_old = {s["text"]: (s["x"], s["y"]) for s in old if s["size"] == 13.0}; letters_new = {s["text"]: (s["x"], s["y"]) for s in new if s["size"] == 14.0}
ck.log("letters a to d at 14 pt at the V26 origins (1 pt, the census is integer)", set(letters_new) == {"a", "b", "c", "d"} and all(abs(letters_old[k][0] - letters_new[k][0]) <= 1 and abs(letters_old[k][1] - letters_new[k][1]) <= 1 for k in letters_old), f"{letters_old} -> {letters_new}")
fonts = {s["font"] for s in new}
ck.log("fonts are the Arial family only", fonts <= {"ArialMT", "Arial-BoldMT", "Arial Bold"}, str(fonts))
ck.log("no printed key strings (round 40 kept: stars, arrow, colour, diamond and dotted-line keys absent)", not any(s["text"].startswith(("* q <", "arrow =", "grey =", "blue =", "open diamond =", "dotted line =")) for s in new))
def cls(s):
    if s["size"] in (13.0, 14.0): return "letter"
    if s["text"] in ("*", "**", "***"): return "star"
    if s["size"] in (11.0, 12.0): return "axis title"
    if re.fullmatch(r"[\d.]+", s["text"]): return "tick"
    if s["size"] in (9.5, 10.5, 10.0) and (s["x"] < 470 or 640 < s["x"] < 1000) and s["y"] < 440: return "cell or head text"
    return "label or head"
shift = {}
for s in old:
    cands = [n for n in new if n["text"] == s["text"] and abs(n["y"] - s["y"]) <= 3]
    if cands:
        n = min(cands, key=lambda n: abs(n["x"] - s["x"])); k = cls(s); shift[k] = max(shift.get(k, 0), abs(n["x"] - s["x"]))
ck.info(f"largest x shift by class (pt, integer census): {shift} (centred and right-anchored strings move with their width, left-anchored labels and letters keep x)")
# ---------------------------------------------------------------- 3 sources and the independent re-derivation
SHA = {}
for f in ("results.csv", "cells.csv", "comparison.csv", "summary.json"):
    sh, sc = sidecar(f"{SRC}/{f}"); SHA[f] = (sh, sc)
ck.log("sidecar sha256 of every step-138 source unchanged since the build, all v8.1", all(SHA[f][0] == E["sources"][f]["sha256"] and SHA[f][1]["step"]["id"] == "138" for f in SHA), "; ".join(f"{f} {SHA[f][1]['output']['mtime_local']}" for f in SHA))
SHA9 = {}
for f in ("attack.json", "attack_min_detectable_hr.csv", "attack_landmark_vs_primary.csv", "attack_landmark1y.csv"):
    sh, sc = sidecar(f"{SRC}/{f}"); SHA9[f] = (sh, sc)
ck.log(f"every step-139 source (part d) is the v8.2 re-extraction (sidecar after {V8_2_CUT}), sha256 equal to the build record", all(SHA9[f][1]["step"]["id"] == "139" and SHA9[f][1]["step"]["rc"] == 0 and SHA9[f][1]["output"]["mtime_local"] >= V8_2_CUT and SHA9[f][0] == E["sources"][f]["sha256"] for f in SHA9), "")
ck.log("nothing on the sheet is kept from V13 (part d from 01_build_d.py v8.2)", E.get("kept_from_v13") in ({}, None) and CR.get("kept_from_v13") is None and "part_d_v8_2" in E)
res = pd.read_csv(hydrated(f"{SRC}/results.csv")); comp = pd.read_csv(hydrated(f"{SRC}/comparison.csv")); summ = jload(f"{SRC}/summary.json")
atk = jload(f"{SRC}/attack.json"); mdhr = pd.read_csv(hydrated(f"{SRC}/attack_min_detectable_hr.csv")); lmvp = pd.read_csv(hydrated(f"{SRC}/attack_landmark_vs_primary.csv"))
HD, CCS = summ["healthy_definition"], summ["cell_counts_side_by_side"]
os.environ["T90_FULL_NIGHTS_ONLY"] = "1"
sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
import cohort_spec  # noqa: E402
from cohort_spec import apply_cohort  # noqa: E402
from disease_definitions import DISEASES, ORGAN_GROUP  # noqa: E402
GROUPS = ["cell1", "cell2", "cell3", "cell4"]
tab = apply_cohort(pd.read_parquet(hydrated(f"{V8TAB}/t90_final.parquet")))
ck.log("v8 table through cohort_spec under the full-nights rule equals the file's cohort n", cohort_spec.FULL_NIGHTS_ONLY and len(tab) == HD["n_cohort"] == cohort_spec.COHORT_N, f"{len(tab)} vs {HD['n_cohort']}")
blocks = [("cardiovascular", [k for k, g in ORGAN_GROUP.items() if g == "Cardiac"] + ["stroke_any", "cvd"]), ("pulmonary", [k for k, g in ORGAN_GROUP.items() if g == "Respiratory"]),
          ("renal", [k for k, g in ORGAN_GROUP.items() if g == "Kidney"]), ("hepatic", [k for k, g in ORGAN_GROUP.items() if g == "Liver"]),
          ("metabolic", [k for k, g in ORGAN_GROUP.items() if g == "Metabolic"]), ("cancer", ["cancer_any", "prostate_ca"])]
ck.log("the healthy filter recorded in summary.json is the ORGAN_GROUP blocks plus stroke and cancer, label for label", all(HD["excluded_on"][b] == [DISEASES[k][0] for k in ks] for b, ks in blocks))
keep = np.ones(len(tab), bool)
for _b, ks in blocks:
    for k in ks: keep &= ~(tab[f"{k}_prevalent"] == 1).values
H = tab[keep]
ck.log("healthy subgroup re-derived from the v8 table equals the file (n_healthy)", len(H) == HD["n_healthy"], f"{len(H)} vs {HD['n_healthy']}")
def cells(d, cut):
    low, pres, shrt, norm = d.spo2_pct_below_90 > 10.0, d.spo2_pct_below_90 <= 1.0, d.TST_min < cut, d.TST_min >= 360.0
    return pd.Series(np.select([shrt & pres, shrt & low, norm & low, norm & pres], GROUPS, default="excluded"), index=d.index)
A = {}; exp_strings = Counter()
for key in ("healthy_tst240", "healthy_tst300"):
    c = cells(H, float(CCS[key]["tst_short_lt_min"]))
    n = {g: int((c == g).sum()) for g in GROUPS}; d = {g: int(H.loc[c == g, "death_incident"].sum()) for g in GROUPS}
    rate = {g: 1000.0 * d[g] / float(H.loc[c == g, "death_years"].sum()) for g in GROUPS}
    A[key] = (n, d, rate)
    for g in GROUPS:
        exp_strings[f"n = {n[g]:,}"] += 1; exp_strings[f"{d[g]} deaths"] += 1; exp_strings[f"{rate[g]:.1f} deaths per 1,000 person-years"] += 1
ck.log("panel a cell sizes and deaths re-derived from the v8 table equal summary.json (both cuts)", all(A[k][0] == CCS[k]["n_by_cell"] and A[k][1] == CCS[k]["deaths_by_cell"] for k in A), "; ".join(f"{k}: n {A[k][0]} deaths {A[k][1]}" for k in A))
ck.log("panel a rates re-derived agree with the builder's record (1e-9)", all(abs(A[k][2][g] - E["panelA"][k]["rate_per_1000"][g]) < 1e-9 for k in A for g in GROUPS))
res300 = res[res.analysis == "healthy_tst300"]
b = res300[~res300.negative_control].dropna(subset=["p_cell1_vs_cell4_hr", "p_cell2_vs_cell4_hr", "p_cell3_vs_cell4_hr"]).sort_values("p_cell2_vs_cell4_hr", ascending=False)
B_ROWS = list(b.disease)
ck.log("panel b row set re-derived (non-control, three contrasts fitted, sorted by the short-sleep low-oxygen hazard ratio) equals the build record", B_ROWS == E["panelB"]["row_order"], f"{len(B_ROWS)} rows")
def stars(q): return "" if q is None or not np.isfinite(q) else ("***" if q < .001 else ("**" if q < .01 else ("*" if q < .05 else "")))
star_expected = Counter(stars(r[f"{pre}_q"]) for _i, r in b.iterrows() for pre in ("p_cell1_vs_cell4", "p_cell2_vs_cell4", "p_cell3_vs_cell4") if stars(r[f"{pre}_q"]))
ns = CCS["healthy_tst300"]["n_by_cell"]
for key in ("cell1", "cell2", "cell3"): exp_strings[f"n = {ns[key]:,}"] += 1
for dis in B_ROWS: exp_strings[dis] += 1
comp300 = comp[(comp.healthy_analysis == "healthy_tst300") & (~comp.negative_control)]
c2 = comp300[comp300.contrast == "cell2_vs_cell4"].dropna(subset=["ratio_healthy_over_full"]).sort_values("full_cohort_hr", ascending=False)
c3 = comp300[(comp300.contrast == "cell3_vs_cell4")].dropna(subset=["ratio_healthy_over_full"]); c3 = c3[~c3.key.isin(set(c2.key))]
ck.log("panel c row sets re-derived (left = short-with-low fitted in both cohorts by full-cohort HR, right = none under v8.1) equal the build record", list(c2.disease) == E["panelC"]["left_order"] and len(c3) == 0 and not E["panelC"]["right_block_drawn"], f"left {len(c2)}, right {len(c3)}")
for dis in c2.disease: exp_strings[dis] += 1
exp_strings[f"Full cohort, n = {HD['n_cohort']:,}"] += 1; exp_strings[f"Healthy subgroup, n = {HD['n_healthy']:,}"] += 1
dd = res300[~res300.negative_control].dropna(subset=["p_cell2_vs_cell1_hr"]).copy()
md_all = mdhr[mdhr.contrast == "cell2_vs_cell1"]; md = md_all.set_index("key"); lm = lmvp.set_index("key")
dd["mdhr80"] = dd.key.map(md.mdhr80); dd["lm_hr"] = dd.key.map(lm.cell2_vs_cell1_hr_L)
dd = dd.sort_values("p_cell2_vs_cell1_hr", ascending=False).reset_index(drop=True)
D_ROWS = list(dd.disease); N_D = len(dd); N_LM = int(np.isfinite(dd.lm_hr.astype(float)).sum()); MED = float(md_all.mdhr80.astype(float).median())
pc = summ["summary"]["healthy_tst300"]["cell2_vs_cell1"]; lc = atk["attacks"]["A4_selection"]["landmark_1y"]["cell2_vs_cell1"]
ck.log("part d row set re-derived (non-control, short-with-low against short-with-preserved, by hazard ratio) equals the build record", D_ROWS == E["panelD"]["row_order"] and N_D == pc["n_conditions_tested"], f"{N_D} rows, landmark {N_LM}")
ck.log("part d: the step-139 cells equal the step-138 cells and the per-row events agree", atk["rebuild"]["cells_tst300"] == CCS["healthy_tst300"]["n_by_cell"] and all(int(md.loc[r.key, "ev_a"]) == int(r.ev_cell2) and int(md.loc[r.key, "ev_ref"]) == int(r.ev_cell1) for r in dd.itertuples()))
ck.log("part d landmark and power values in the build record equal the current files (1e-9)", all(abs(float(r["mdhr80"]) - float(md.loc[r["key"], "mdhr80"])) < 1e-9 for r in E["panelD"]["rows"]) and abs(E["panelD"]["median_mdhr80_all"] - MED) < 1e-9)
for t in (f"{pc['n_significant']} of {pc['n_conditions_tested']} significant,", f"{pc['n_significant_fdr']} after FDR", f"{lc['n_sig']} of {lc['n']} significant,", f"{lc['n_fdr']} after FDR"): exp_strings[t] += 1
for dis in D_ROWS: exp_strings[dis] += 1
for t in E["panelD"]["axis"]["ticks_new"]: exp_strings[f"{t:g}"] += 3
new_strings = Counter(s["text"] for s in new)
miss = {k: v for k, v in exp_strings.items() if new_strings.get(k, 0) < v}
ck.log("read-back: every re-derived value string (panel a counts, deaths, rates; b and c heads with n; d heads; row labels of b, c, d; d ticks) is printed at least the expected number of times", not miss, str(miss))
star_new = Counter(s["text"] for s in new if s["text"] in ("*", "**", "***"))
ck.log("panel b stars: the multiset of printed stars equals the re-derived stars", star_new == star_expected, f"printed {dict(star_new)} expected {dict(star_expected)}")
ck.log("builder's expected values all printed", all(new_strings.get(v["text"], 0) >= 1 for v in E["values"]), str([v["text"] for v in E["values"] if new_strings.get(v["text"], 0) < 1][:5]))
# row labels in order (top to bottom) per panel from the census
def labels_in(x0, y0, y1): return [s["text"] for s in sorted(new, key=lambda s: s["y"]) if abs(s["x"] - x0) <= 1.5 and y0 < s["y"] < y1 and s["size"] in (11.0, 10.0) and s["font"] == "ArialMT" and not re.fullmatch(r"[\d.]+", s["text"])]
XS = CR["place_v13"]["b"][0]; YL = CR["r40_lower_row_y"]
ck.log("panel b row labels on the sheet are the re-derived rows in order", labels_in(XS + 12.96, 40, 14 + CR["parts"]["part_b.pdf"]["size"][1] - 60) == B_ROWS)
ck.log("panel c row labels on the sheet are the re-derived rows in order", labels_in(26.96, YL + 60, YL + CR["parts"]["part_c.pdf"]["size"][1] - 60) == list(c2.disease))
ck.log("part d row labels on the sheet are the re-derived rows in order", labels_in(XS + 12.96, YL + 60, YL + CR["parts"]["part_d.pdf"]["size"][1] - 60) == D_ROWS)
# ---------------------------------------------------------------- 4 vector layers of the flat parts against the round-40 parts
dr_new = {p: C.flat_drawings(f"{W}/part_{p}.pdf") for p in "abcd"}; dr_old = {p: C.flat_drawings(f"{R40W}/part_{p}.pdf") for p in "abcd"}
for p in "abcd":
    ok_d, det_d = C.drawings_compare(dr_old[p], dr_new[p], [])
    ck.log(f"part {p}: vector layer identical to the round-40 part (every path rect, fill, stroke, width, dashes)", ok_d, f"old {det_d['n_old']} new {det_d['n_new']} paths; only old {det_d['only_old'][:3]} only new {det_d['only_new'][:3]}")
a_rect = [d for d in dr_new["a"] if d["kinds"] == "re" and d["fill"] in ("#0288d1", "#eef0f1")]
ck.log("panel a: two low-oxygen cells and two preserved cells per diagram", Counter(d["fill"] for d in a_rect) == Counter({"#0288d1": 4, "#eef0f1": 4}), str(Counter(d["fill"] for d in a_rect)))
nb = Counter(); [nb.update({d["stroke"]: 1}) for d in dr_new["b"] if d["kinds"] == "l" and d["width"] and abs(d["width"] - 1.7) < 0.05 and d["stroke"] in ("#8a9099", "#0288d1", "#3f9fd8")]
ck.log("panel b: one CI line per row per series (grey, blue, mid-blue)", all(v == len(B_ROWS) for v in nb.values()) and len(nb) == 3, f"{dict(nb)} rows {len(B_ROWS)}")
nd_ = Counter(); [nd_.update({d["stroke"]: 1}) for d in dr_new["d"] if d["kinds"] == "l" and d["width"] and abs(d["width"] - 1.7) < 0.05 and d["stroke"] in ("#8a9099", "#0288d1")]
diamonds = sum(1 for d in dr_new["d"] if d["fill"] == "#ffffff" and d["stroke"] == "#3f9fd8")
ck.log("part d drawing layer: one blue observed CI per row, one grey landmark CI per landmark row, one open diamond per row", nd_["#0288d1"] == N_D and nd_["#8a9099"] == N_LM and diamonds == N_D, f"{dict(nd_)} diamonds {diamonds}")
cols_old = {c for p in "abcd" for d in dr_old[p] for c in (d["fill"], d["stroke"]) if c}; cols_new = {c for p in "abcd" for d in dr_new[p] for c in (d["fill"], d["stroke"]) if c}
ck.log("no colour introduced", cols_new <= cols_old, f"new-only {sorted(cols_new - cols_old)}")
# ---------------------------------------------------------------- 5 renders and the masked raster diff
C.gs_render(NEW, f"{D}/{S}_150dpi.png", 150); C.gs_render(OLD, f"{W}/{S}_V26_150dpi.png", 150)
boxes = [C.text_box(s, pad=2.5) for s in old + new]
rd = C.raster_diff(f"{W}/{S}_V26_150dpi.png", f"{D}/{S}_150dpi.png", boxes, 150)
ck.log("150 dpi raster: every differing pixel lies inside a text box (V26 or new census); data marks, axes, bands, cells identical", rd["shape_equal"] and rd["n_outside"] == 0, str(rd))
for name, box in (("a_cells", [0, 40, 500, 250]), ("b_heads_stars", [640, 40, 1021, 200]), ("b_gutter", [520, 90, 700, 400]), ("c_heads_gutter", [0, 500, 480, 700]), ("d_heads_foot", [640, 500, 1021, 620]), ("d_foot", [640, 1000, 1021, 1096]), ("row_gap", [0, 420, 1021, 560])):
    C.crop_png(f"{W}/{S}_V26_150dpi.png", 150, box, f"{VER}/crops/{S}_{name}_OLD.png"); C.crop_png(f"{D}/{S}_150dpi.png", 150, box, f"{VER}/crops/{S}_{name}_NEW.png")
ck.info("crops OLD (V26) and NEW in verify/crops/")
ck.info(f"new sheet sha256 {C.sha256(NEW)}; V26 sha256 {C.sha256(OLD)}")
json.dump({"exceptions": [list(k) + [v] for k, v in exceptions.items()], "census": det, "raster": rd, "x_shift": shift}, open(f"{VER}/{S}_verify_record.json", "w"), indent=1)
ok_all = ck.write(f"{VER}/checks.txt"); print("RESULT ALL PASS" if ok_all else "RESULT FAIL"); sys.exit(0 if ok_all else 1)
