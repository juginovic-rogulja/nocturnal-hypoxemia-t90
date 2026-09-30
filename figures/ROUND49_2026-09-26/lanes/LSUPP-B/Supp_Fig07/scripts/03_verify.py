#!$T90_PY
"""Supp Fig 7, round 49 (+1 pt) verification against V26. Page box (pypdf); Ghostscript text census against V26 (same strings, every
size +1.0, the declared exceptions); independent re-derivation of every printed value from the v8.1 tst_ceiling files, the pairs
file and ranking_v3.csv (the round-37 verifier's derivation) matched against the census; the vector layer of the flat sheet
(work/sheet.pdf, letters excluded) against the round-40 flat sheet; masked raster diff of the composed sheets at 150 dpi; the
150 dpi render. Writes verify/checks.txt (ends RESULT ALL PASS)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import json, os, re, sys
from collections import Counter
import numpy as np, pandas as pd
L = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LSUPP-B"
sys.path.insert(0, f"{L}/scripts"); sys.path.insert(0, f"{L}/Supp_Fig07/scripts")
import r49_common as C
from l3b_common import NUM, SV, sidecar, jload, hydrated  # noqa: E402
S = "Supp_Fig07"; D = f"{L}/{S}"; W = f"{D}/work"; VER = f"{D}/verify"; os.makedirs(f"{VER}/crops", exist_ok=True)
OLD = f"{C.V26}/{S}.pdf"; NEW = f"{D}/{S}.pdf"
R40W = f"{paths.FIGURE_ROOT}/ROUND40_2026-09-18/lanes/LSUPP/Supp_Fig07/work"
E = json.load(open(f"{W}/expected_values.json")); N_MEAS = int(E["n_meas"]); LET = json.load(open(f"{W}/letters_v16.json"))
ck = C.Checks(f"{S} round 49 (+{C.PT_PLUS} pt) verification against V26, {C.now()}")
n_new, w_new, h_new = C.page_box(NEW); n_old, w_old, h_old = C.page_box(OLD)
ck.log("one page, page box equals V26 (0.05 pt)", n_new == 1 and abs(w_new - w_old) < 0.05 and abs(h_new - h_old) < 0.05, f"V26 {w_old} x {h_old}, new {w_new} x {h_new}")
old = C.merge_kerned(C.gs_text(OLD, f"{W}/{S}_V26_txtwrite.xml")); new = C.merge_kerned(C.gs_text(NEW, f"{W}/{S}_NEW_txtwrite.xml"))
json.dump({"old": old, "new": new}, open(f"{VER}/census_spans.json", "w"), indent=0)
PC = E["r49_panel_c"]; ROWS_C = ["Nocturnal oxygen (T90), as measured", "Total sleep time, as measured", "Total sleep time, error removed entirely"]
exceptions = {(t, "ArialMT", 10.0): PC["row_label_size"] for t in ROWS_C} if PC["kept"] else {}
ok, det = C.census_compare(old, new, exceptions)
ck.log("census: the same multiset of strings as V26", det["strings_equal"], f"V26 {det['n_old']} spans, new {det['n_new']} (kerning splits merged)")
ck.log(f"census: every string one point larger than on V26 except the {len(exceptions)} declared strings (panel c row labels kept at {PC['row_label_size']} pt)", ok, f"missing {det['missing']} extra {det['extra']}")
ck.info(f"panel c: axis left edge pinned at {PC['axis_left_pt']:.2f} pt (V26 195.69); at 11 pt the longest row label would start {PC['longest_label_start_at_plus_pt']:.1f} pt from the page edge (margin 13 pt), so the three labels keep {PC['row_label_size']} pt")
ck.info(f"size table V26 {C.size_table(old)}; new {C.size_table(new)}")
letters_old = {s["text"]: (s["x"], s["y"]) for s in old if s["size"] == 13.0}; letters_new = {s["text"]: (s["x"], s["y"]) for s in new if s["size"] == 14.0}
ck.log("letters a to c at 14 pt at the V26 origins (1 pt, the census is integer)", set(letters_new) == {"a", "b", "c"} and all(abs(letters_old[k][0] - letters_new[k][0]) <= 1 and abs(letters_old[k][1] - letters_new[k][1]) <= 1 for k in letters_old), f"{letters_old} -> {letters_new}")
ck.log("'top 10' labels at 10 pt (were 9, the type floor + 1)", sorted(s["size"] for s in new if s["text"] == "top 10") == [10.0, 10.0], str([s["size"] for s in new if s["text"] == "top 10"]))
fonts = {s["font"] for s in new}
ck.log("fonts are the Arial family only", fonts <= {"ArialMT", "Arial-BoldMT", "Arial Bold"}, str(fonts))
ck.log("the six round-31 prose strings and the round-40 legend strings stay absent", not any(t in s["text"] for s in new for t in ("Measured projection", "Square-root projection", "untreated same-person", "observed reliability", "supported scenarios")))
def cls(s):
    if s["size"] in (13.0, 14.0): return "letter"
    if s["size"] in (11.0, 12.0): return "axis title"
    if re.fullmatch(r"[\d.]+", s["text"]): return "tick"
    return "label"
shift = {}
for s in old:
    cands = [n for n in new if n["text"] == s["text"] and abs(n["y"] - s["y"]) <= 3]
    if cands:
        n = min(cands, key=lambda n: abs(n["x"] - s["x"])); k = cls(s); shift[k] = max(shift.get(k, 0), abs(n["x"] - s["x"]))
ck.info(f"largest x shift by class (pt, integer census): {shift}")
# ---------------------------------------------------------------- read-back from the files
TC = f"{SV}/tst_ceiling"
for f in ("attack_summary.json", "summary.json", "_iccs.json", "ranking_reproduced.csv", "effect_vs_gain.csv"):
    sh, sc = sidecar(f"{TC}/{f}"); assert sh == E["sources"][f]["sha256"], f
sh, sc = sidecar(f"{NUM}/ranking_v3.csv"); assert sh == E["sources"]["ranking_v3.csv"]["sha256"]
sh, sc = sidecar(f"{NUM}/tst_ceiling.json"); assert sh == E["sources"]["tst_ceiling.json"]["sha256"]
ck.log("sidecar sha256 of every source unchanged since the build, every sidecar v8.1", True)
A = jload(f"{TC}/attack_summary.json"); Sm = jload(f"{TC}/summary.json"); I = jload(f"{TC}/_iccs.json"); T43 = jload(f"{NUM}/tst_ceiling.json")
RK = pd.read_csv(hydrated(f"{TC}/ranking_reproduced.csv")); EVG = pd.read_csv(hydrated(f"{TC}/effect_vs_gain.csv")); V3 = pd.read_csv(hydrated(f"{NUM}/ranking_v3.csv"), comment="#")
w = pd.read_parquet(hydrated(f"{SV}/data/untreated_wide.parquet")); DIFF = (w.tst_min_n1 - w.tst_min_n2).abs()
SLOPE, INTER = Sm["map_primary"]["slope"], Sm["map_primary"]["intercept"]; E_NOW = float(EVG[EVG.feature == "TST_min"].E_meanabs.iloc[0]); EXP = A["A1c"]["exponent_used"]
def rank_at(icc, p):
    dc = SLOPE * (E_NOW * (1 / icc) ** p) ** 2 + INTER; return int((RK.dC > dc).sum() + 1)
icc = I["TST"]["icc_raw"]; ci = A["A1"]["icc_TST_ci"]
sq = sorted(v["rank_ceiling"] for k, v in A["A1"]["ceilings"].items() if not k.startswith("classical")); me = sorted(v["rank_ceiling"] for v in A["A1c"]["ceilings"].values())
n_meas = len(RK); assert n_meas == N_MEAS == int(Sm["n_measures"]) == int(T43["n_measures"]) == len(V3)
r_obs_sqrt = rank_at(icc, 0.5)
exp_strings = Counter({f"median {DIFF.median():.1f} min": 1, f"IQR {DIFF.quantile(.25):.1f} to {DIFF.quantile(.75):.1f}": 1, f"ICC {icc:.3f} (95% CI, {ci[0]:.3f}-{ci[1]:.3f})": 1,
                       f"rank {int(V3[V3.feature == 'spo2_pct_below_90']['rank'].iloc[0])}": 1, f"rank {int(V3[V3.feature == 'TST_min']['rank'].iloc[0])}": 1,
                       f"Rank at perfect measurement, of {n_meas}": 1, f"Rank among the {n_meas} sleep measurements": 1, f"{n_meas}": 2, "top 10": 2, **{f"{t}": 2 for t in E["rank_ticks"][:-1]}})
new_strings = Counter(s["text"] for s in new)
missing = {k: v for k, v in exp_strings.items() if new_strings.get(k, 0) < v}
ck.log("read-back: every value re-derived from the sources is printed (median, IQR, ICC, the two ranks from ranking_v3.csv, the count twice as the last rank tick and in both axis titles, the rank ticks twice, top 10 twice)", not missing, f"missing {missing}")
ck.log(f"rank tick {r_obs_sqrt} (the square-root projection's rank at the observed reliability, round-40 follow-up) is on both rank axes and no 24 or 100 tick", r_obs_sqrt in E["rank_ticks"] and new_strings.get(str(r_obs_sqrt), 0) == 2 and not any(s["text"] in ("24", "100") for s in new), str(E["rank_ticks"]))
ck.log("builder's expected values all printed and within the independent re-derivation", all(new_strings.get(v["text"], 0) >= 1 for v in E["values"]) and set(v["text"] for v in E["values"]) <= set(exp_strings), str([v["text"] for v in E["values"] if v["text"] not in exp_strings]))
ck.log("supported range of the ceiling frame equals tst_ceiling.json defensible_rank_range", [min(sq + me), max(sq + me)] == list(T43["defensible_rank_range"]), f"{[min(sq + me), max(sq + me)]}")
ck.log("observed-reliability ranks reproduce attack_summary.json", r_obs_sqrt == E["facts_for_legend"]["rank_obs_sqrt"] and rank_at(icc, EXP) == E["facts_for_legend"]["rank_obs_meas"], f"{r_obs_sqrt}, {rank_at(icc, EXP)}")
# ---------------------------------------------------------------- the vector layer of the flat sheet against the round-40 flat sheet
dr_new = C.flat_drawings(f"{W}/sheet.pdf"); dr_old = C.flat_drawings(f"{R40W}/sheet.pdf")
ok_d, det_d = C.drawings_compare(dr_old, dr_new, [])
ck.log("vector layer of the flat sheet identical to the round-40 sheet (every path rect, fill, stroke, width, dashes: histogram, bands, median line, curve, dotted line, dots, range bar, axes)", ok_d, f"old {det_d['n_old']} new {det_d['n_new']} paths; only old {det_d['only_old'][:3]} only new {det_d['only_new'][:3]}")
bars = [d for d in dr_new if d["fill"] == "#ccd1d6" and d["kinds"] == "re"]
ck.log("histogram bars: one per non-empty 15-min bin", len(bars) == sum(1 for c in E["hist_counts"] if c > 0), f"{len(bars)} bars")
ck.log("panel b: no measured-projection blue line (#3f9fd8 line), one green curve", not any(d["stroke"] == "#3f9fd8" and d["kinds"] == "l" and d["n"] > 10 for d in dr_new) and sum(1 for d in dr_new if d["stroke"] == "#298d32" and d["kinds"] == "l" and d["n"] > 10) == 1)
cols_old = {c for d in dr_old for c in (d["fill"], d["stroke"]) if c}; cols_new = {c for d in dr_new for c in (d["fill"], d["stroke"]) if c}
ck.log("no colour introduced", cols_new <= cols_old, f"new-only {sorted(cols_new - cols_old)}")
# ---------------------------------------------------------------- renders and the masked raster diff
C.gs_render(NEW, f"{D}/{S}_150dpi.png", 150); C.gs_render(OLD, f"{W}/{S}_V26_150dpi.png", 150)
boxes = [C.text_box(s, pad=2.5) for s in old + new]
rd = C.raster_diff(f"{W}/{S}_V26_150dpi.png", f"{D}/{S}_150dpi.png", boxes, 150)
ck.log("150 dpi raster: every differing pixel lies inside a text box (V26 or new census); data marks and axes identical", rd["shape_equal"] and rd["n_outside"] == 0, str(rd))
for name, box in (("a_stats", [100, 20, 484, 130]), ("b_axis", [0, 250, 484, 480]), ("c_labels", [0, 500, 484, 670])):
    C.crop_png(f"{W}/{S}_V26_150dpi.png", 150, box, f"{VER}/crops/{S}_{name}_OLD.png"); C.crop_png(f"{D}/{S}_150dpi.png", 150, box, f"{VER}/crops/{S}_{name}_NEW.png")
ck.info("crops OLD (V26) and NEW in verify/crops/")
ck.info(f"new sheet sha256 {C.sha256(NEW)}; V26 sha256 {C.sha256(OLD)}")
json.dump({"exceptions": [list(k) + [v] for k, v in exceptions.items()], "census": det, "raster": rd, "x_shift": shift, "panel_c": PC}, open(f"{VER}/{S}_verify_record.json", "w"), indent=1)
ok_all = ck.write(f"{VER}/checks.txt"); print("RESULT ALL PASS" if ok_all else "RESULT FAIL"); sys.exit(0 if ok_all else 1)
