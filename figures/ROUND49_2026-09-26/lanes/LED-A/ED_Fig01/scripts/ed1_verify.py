#!/usr/bin/env python3
"""03_verify for ED_Fig01 (lane V14_L1_RANK, round 37: v8.1 numbers against the V13 sheet; Proof W, whole sheet re-plotted). Page box (the declared 40.11 pt growth), letters, font faces, geometry pins on both
panels (spines, ticks, bars, leaders and dots at the positions the numbers put them, the static pins at the base sheet's positions), the
text-layer proof (base tokens minus the strings the OLD ranking file prints equal new tokens minus the strings the NEW file prints), the
read-back of every placed string against an independent re-derivation from the files, banned words, crops, the 150 dpi render and
verify/CHANGES_ED_Fig01.csv."""
import os, sys, json
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C, vcommon as V

S = "ED_Fig01"; SD = f"{C.LANE}/{S}"; WORK = f"{SD}/work"; VER = f"{SD}/verify"; CROPS = f"{VER}/crops"; os.makedirs(CROPS, exist_ok=True)
BASE_PDF = C.hydrated(f"{C.BASE}/{S}.pdf"); NEW_PDF = C.hydrated(f"{SD}/{S}.pdf")
check, finish = C.checks_writer(f"{VER}/checks.txt")
bsp, brect, blets = V.sheet_spans(BASE_PDF); nsp, nrect, nlets = V.sheet_spans(NEW_PDF)
DA = json.load(open(f"{VER}/{S}_a_drawn.json")); DB = json.load(open(f"{VER}/{S}_b_drawn.json")); GA = DA["geometry"]; GB = DB["geometry"]
GROW = GB["grow_pt"]; PITCH = GB["pitch"]

# ---------------------------------------------------------------- independent re-derivation from the files (old and new)
def rows_from(rk_path, fam_path):
    rk = pd.read_csv(C.hydrated(rk_path), comment="#"); fam = pd.read_csv(C.hydrated(fam_path)).set_index("feature")["family"]
    rk["family"] = fam.loc[rk.feature].values; rk["g"] = rk.dC * 1000.0; rki = rk.set_index("feature"); raw = set(rk.feature)
    pairs = [(f, f[:-6]) for f in rk.feature if f.endswith("_siteZ") and f[:-6] in raw]
    if pairs: assert max(abs(float(rki.loc[f, "dC"]) - float(rki.loc[b, "dC"])) for f, b in pairs) < 1e-5     # the v7 file carries 60 twins, the v8.1 file none
    distinct = rk[~rk.feature.isin([f for f, _ in pairs])]
    tail = distinct[distinct.g > 3.0].sort_values("g", ascending=False).reset_index(drop=True)
    return rk, tail
NICE = dict(GB["nice_all"]); NICE.update({r["feature"]: r["name"] for r in GB["rows"]})
RK_NEW, TAIL_NEW = rows_from(f"{C.NUM}/ranking_v3.csv", f"{C.NUM}/measure_families.csv")
OLD_FAM = f"{C.OLDNUM}/numbers/measure_families.csv" if os.path.exists(f"{C.OLDNUM}/numbers/measure_families.csv") else f"{C.NUM}/measure_families_PRE_V8.csv"   # the v7 family map (197 rows, siteZ twins) for the OLD side only
RK_OLD, TAIL_OLD = rows_from(f"{C.OLDNUM}/numbers/ranking_v3.csv", OLD_FAM)
OLD_NICE = dict(NICE); OLD_NICE.setdefault("AHI", "Apnea-hypopnea index")
assert all(f in OLD_NICE for f in TAIL_OLD.feature), sorted(set(TAIL_OLD.feature) - set(OLD_NICE))
RES_NEW = json.load(open(C.hydrated(f"{C.NUM}/results_v2.json"))); RES_OLD = json.load(open(C.hydrated(f"{C.OLDNUM}/numbers/results_v2.json")))
BAND_KEYS = ["0-1%", "1-5%", "5-10%", ">10%"]

# ---------------------------------------------------------------- 1 page box, 2 letters, 3 fonts
N_ROWS = GB["n_rows"]; DELTA_V13 = (N_ROWS - 30) * PITCH
check(abs(nrect[0] - brect[0]) < 0.01 and abs(nrect[1] - (677.81 + GROW)) < 0.02 and abs(brect[1] - 717.92) < 0.02 and abs(GROW - (N_ROWS - 28) * PITCH) < 1e-6,
      f"page box: width equal, height = the 28-row base plus {N_ROWS - 28} row pitches (V13 had 30 rows: {DELTA_V13:+.2f} pt against V13, the rules 3.3 growth rule)", f"{brect} -> {nrect}, grow {GROW:.3f} pt = {N_ROWS - 28} x {PITCH:.4f}")
ok, msg = V.compare_letters(blets, nlets); check(ok, "letters a, b within 0.3 pt, same face and size", msg)
V.fonts_check(check, bsp, nsp)

# ---------------------------------------------------------------- 4 rows: the build's row list equals the independent re-derivation
build_rows = [(r["feature"], r["printed"], r["family"]) for r in GB["rows"]]
rederived = [(r.feature, f"{r.g:.1f}", r.family) for _, r in TAIL_NEW.iterrows()]
check(build_rows == rederived and len(rederived) == N_ROWS, f"panel b row set, order, printed gains and families equal an independent re-derivation ({N_ROWS} rows)", f"{len(build_rows)} rows; first mismatch {next((i for i, (a, b) in enumerate(zip(build_rows, rederived)) if a != b), None)}")
band_ok = all(int(RES_NEW["cohort"]["t90_bands_n"][k]) == b["n"] and abs(float(RES_NEW["cohort"]["t90_bands_pct"][k]) - b["pct"]) < 0.051 for k, b in zip(BAND_KEYS, GA["bands"]))
check(band_ok and sum(b["n"] for b in GA["bands"]) == RES_NEW["cohort"]["n"], "panel a band counts and shares equal results_v2.json (sum = cohort n)", str([(b["label"], b["n"], round(b["pct"], 2)) for b in GA["bands"]]))

# ---------------------------------------------------------------- 5 geometry pins on the new sheet
G = V.geometry_records(NEW_PDF); INK = "#1a1d21"
pins = [("a left spine", 79.59, 33.122, 79.59, 263.522, INK), ("a bottom spine", 79.59, 263.522, 147.99, 263.522, INK)]
for v in (0, 20, 40, 60): pins.append((f"a x tick {v}", 79.59 + v * GA["ppu"], 263.522, 79.59 + v * GA["ppu"], 266.522, INK))
for k in range(4):
    cy = 33.122 + (3.72 - (3 - k)) * GA["pitch"]; pins.append((f"a y tick band {k}", 76.59, cy, 79.59, cy, INK))
yb = GB["y1"]; pins.append(("b bottom spine (moved down by the growth)", 292.59, yb, 371.05, yb, INK))
for v in (5, 10, 15, 20):
    x = 292.59 + (v - 2.85) * (371.05 - 292.59) / (24.5 - 2.85); pins.append((f"b x tick {v}", x, yb, x, yb + 3.0, INK))
missing = [p[0] for p in pins if not V.has_segment(G, *p[1:5], stroke=p[5])]
check(not missing and abs(yb - (633.17 + GROW)) < 0.02, f"static geometry pins ({len(pins)} spines and ticks at the base positions, panel b's axis {GROW:.2f} pt lower)", f"missing {missing}")
bar_miss = [b["label"] for v, b in zip(DA["values"], GA["bands"]) if v["kind"] == "drawn" and not V.has_rect(G, *v["bar_rect"], fill=v["colour"], tol=0.3)]
check(not bar_miss, "panel a bars: four rectangles at x0 = 79.59 with width = share x 1.14 pt and the T90 ramp colours", f"missing {bar_miss}")
dot_miss, lead_miss = [], []
for v in [v for v in DB["values"] if v["kind"] == "printed"]:
    cx, cy = v["dot_page_xy"]
    if not V.has_dot(G, cx, cy, v["colour"]): dot_miss.append((v["name"], cx, cy))
    if not V.has_segment(G, 293.134, cy, cx, cy, stroke="#b3dcf2", tol=0.35): lead_miss.append((v["name"], cx, cy))
check(not dot_miss and not lead_miss, f"panel b: {N_ROWS} dots at (gain, row) in the family colours and {N_ROWS} dotted leaders from 3.0 to the dot", f"dots missing {dot_miss[:5]} leaders missing {lead_miss[:5]}")

# ---------------------------------------------------------------- 6 text layer: base minus old strings == new minus new strings
old_strings = [OLD_NICE[r.feature] for _, r in TAIL_OLD.iterrows()] + [f"{r.g:.1f}" for _, r in TAIL_OLD.iterrows()] + [r.family for _, r in TAIL_OLD.iterrows()]
new_strings = [r["name"] for r in GB["rows"]] + [r["printed"] for r in GB["rows"]] + [r["family"] for r in GB["rows"]]
res = V.residue_check(check, "text layer: base tokens minus the OLD file's strings equal new tokens minus the NEW file's strings (static text identical)", bsp, nsp, old_strings, new_strings)
print("  words lost/gained:", res["words_lost_gained"]); print("  numbers lost/gained:", res["numbers_lost_gained"])
dup = [k for k, v in __import__("collections").Counter((s["text"], s["origin"][0], s["origin"][1]) for s in nsp).items() if v > 1]
check(not dup, "no duplicated spans on the new sheet", str(dup[:5]))

# ---------------------------------------------------------------- 7 read-back of every placed string, second lines of the wrapped names
V.readback(check, "read-back panel a", nsp, DA["values"]); V.readback(check, "read-back panel b", nsp, DB["values"])
second = []
for r in GB["rows"]:
    if len(r["lines"]) == 2:
        y = GB["y_row0"] + PITCH * (r["row"] - 1) + GB["two_dy"][1]
        s, d = V.find_span(nsp, r["lines"][1], GB["name_x"], y, ha="left", tol=0.6)
        if s is None: second.append((r["lines"][1], y, d))
check(not second, "second lines of the wrapped names at row centre + 8.15 pt", f"fails {second}")
present = {s["text"] for s in nsp}
miss = sorted(t for t in new_strings if t not in present and "\n" not in t and t in NICE.values() and all(len(r["lines"]) == 1 for r in GB["rows"] if r["name"] == t))
check(not miss, "every one-line name, printed gain and family string of the re-derived row set is on the sheet", f"missing {miss}")
bad = C.banned_words([s["text"] for s in nsp]); check(not bad, "no banned words on the sheet", str(bad))

# ---------------------------------------------------------------- crops, 150 dpi render, CHANGES
H_new = nrect[1]
C.crop_pair(BASE_PDF, NEW_PDF, (0, 0, 165, 310), f"{CROPS}/{S}_a_old_vs_new_200dpi.png", dpi=200)
C.crop_pair(BASE_PDF, NEW_PDF, (165, 0, 518.74, 370), f"{CROPS}/{S}_b_top_old_vs_new_200dpi.png", dpi=200)
C.crop_pair(BASE_PDF, NEW_PDF, (165, 340, 518.74, H_new), f"{CROPS}/{S}_b_bottom_old_vs_new_200dpi.png", dpi=200)
C.crop_pair(BASE_PDF, NEW_PDF, (0, 0, 518.74, H_new), f"{CROPS}/{S}_whole_old_vs_new_150dpi.png", dpi=150)
C.render_png(NEW_PDF, f"{SD}/{S}_150dpi.png", dpi=150)

rows = []
for k, b in zip(BAND_KEYS, GA["bands"]):
    o_n, o_p = int(RES_OLD["cohort"]["t90_bands_n"][k]), float(RES_OLD["cohort"]["t90_bands_pct"][k])
    rows.append(dict(sheet=S, panel="a", item=f"band {b['label']} share (drawn bar) and count", old_value=f"{o_p} % (n {o_n:,})", new_value=f"{b['pct']:.1f} % (n {b['n']:,})", source_file="numbers/results_v2.json",
                     key_path=f"cohort.t90_bands_pct[{k}], cohort.t90_bands_n[{k}]", dramatic="YES" if abs(b["n"] - o_n) / o_n > 0.2 else "", note=f"count {100 * (b['n'] - o_n) / o_n:+.1f}%"))
old_pos = {r.feature: i + 1 for i, r in TAIL_OLD.iterrows()}; new_pos = {r["feature"]: r["row"] for r in GB["rows"]}
old_g = {r.feature: f"{r.g:.1f}" for _, r in TAIL_OLD.iterrows()}; new_g = {r["feature"]: r["printed"] for r in GB["rows"]}
RKO = RK_OLD.set_index("feature"); RKN = RK_NEW.set_index("feature")
for f in list(dict.fromkeys(list(new_pos) + list(old_pos))):
    name = NICE.get(f, OLD_NICE.get(f))
    o = (f"row {old_pos[f]}: {old_g[f]} (rank {int(RKO.loc[f, 'rank'])})" if f in old_pos else (f"not in the row set (gain {RKO.loc[f, 'g']:.2f}, rank {int(RKO.loc[f, 'rank'])})" if f in RKO.index else "not a v7 measurement (added in v8, decision 3)"))
    n = (f"row {new_pos[f]}: {new_g[f]} (rank {int(RKN.loc[f, 'rank'])})" if f in new_pos else (f"not in the row set (gain {RKN.loc[f, 'g']:.2f}, rank {int(RKN.loc[f, 'rank'])})" if f in RKN.index else "not a v8.1 measurement (left the list, decision 3)"))
    note = ""
    if f not in old_pos: note = "enters the row set (gain above 0.003)" + (", new plain-English name on this sheet" if f in GB.get("new_names", {}) else "")
    elif f not in new_pos: note = "leaves the row set"
    elif f == "AHI": note = "notable: a named callout of Fig 1c"
    rows.append(dict(sheet=S, panel="b", item=f"{name} ({f})", old_value=o, new_value=n, source_file="numbers/ranking_v3.csv x measure_families.csv", key_path=f"dC x 1000 and rank where feature == {f}",
                     dramatic="notable" if note else "", note=note))
rows.append(dict(sheet=S, panel="b", item="row count", old_value=len(TAIL_OLD), new_value=len(GB["rows"]), source_file="numbers/ranking_v3.csv", key_path="count of distinct measurements with dC > 0.003 (siteZ twins folded)",
                 dramatic="notable", note=f"page {DELTA_V13:+.2f} pt against V13 (717.92) to {nrect[0]} x {nrect[1]} pt ({nrect[1] / 72 * 25.4:.1f} mm tall against the 247 mm ceiling), row pitch kept"))
top_old = RK_OLD.sort_values("rank").feature[:10].tolist(); top_new = RK_NEW.sort_values("rank").feature[:10].tolist()
rows.append(dict(sheet=S, panel="b", item="top 10 of the ranking", old_value=" | ".join(top_old), new_value=" | ".join(top_new), source_file="numbers/ranking_v3.csv", key_path="rank 1 to 10",
                 dramatic="no" if top_old == top_new else "YES", note="top 10 unchanged" if top_old == top_new else "top 10 changed"))
C.write_changes(f"{VER}/CHANGES_{S}.csv", rows)
print(f"  CHANGES rows: {len(rows)}, dramatic YES: {sum(1 for r in rows if r['dramatic'] == 'YES')}, notable: {sum(1 for r in rows if r['dramatic'] == 'notable')}")
finish()
