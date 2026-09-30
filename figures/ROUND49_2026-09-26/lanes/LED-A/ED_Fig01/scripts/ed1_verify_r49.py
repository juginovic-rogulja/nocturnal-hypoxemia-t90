#!/usr/bin/env python3
"""ROUND 49, lane LED-A: 03_verify for ED_Fig01 at +1 pt. Adapted copy of ROUND37_2026-09-14/figures/lanes/V14_L1_RANK/scripts/ed1_verify.py
(Proof W against the V13 sheet and the numbers files). Declared deltas of this round, applied to every pin: the 16 pt title strip (every
V13 object 16 pt lower, DY), panel b 14.4 pt further right (XB, round 43), every text +1 pt (DELTA), the panel a x-axis title as two centred
lines (round 40), the names the wrap rule re-splits at 11 pt (builder record), the 1.22 pt wider two-line pitch where the next row is one line.
The data checks (rows re-derived from ranking_v3.csv and measure_families.csv, bands from results_v2.json, dots and bars at the positions
the numbers put them, every placed string read back) are the round-37 checks, unchanged in substance. The OLD-vs-NEW CHANGES block of the
round-37 verifier (v7 snapshot) is not part of this round: CHANGES_ED_Fig01.csv lists this round's typographic deltas (written by run_sheet.py)."""
import os, sys, json
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C, vcommon as V

S = "ED_Fig01"; SD = f"{C.LANE}/{S}"; WORK = f"{SD}/work"; VER = f"{SD}/verify"; CROPS = f"{VER}/crops"; os.makedirs(CROPS, exist_ok=True)
BASE_PDF = C.hydrated(f"{C.BASE}/{S}.pdf"); NEW_PDF = C.hydrated(f"{SD}/{S}.pdf")
check, finish = C.checks_writer(f"{VER}/checks_03_verify.txt")
bsp, brect, blets = V.sheet_spans(BASE_PDF); nsp, nrect, nlets = V.sheet_spans(NEW_PDF)
DA = json.load(open(f"{VER}/{S}_a_drawn.json")); DB = json.load(open(f"{VER}/{S}_b_drawn.json")); GA = DA["geometry"]; GB = DB["geometry"]
GROW = GB["grow_pt"]; PITCH = GB["pitch"]; XB = GB["x_shift_panel_b"]; DY = 16.0; DELTA = C.DELTA
STAMP = json.load(open(f"{SD}/{S}_stamp_record.json")); assert abs(STAMP["strip_pt"] - DY) < 1e-9 and abs(STAMP["size"] - (13.0 + DELTA)) < 1e-9

# ---------------------------------------------------------------- independent re-derivation from the files
def rows_from(rk_path, fam_path):
    rk = pd.read_csv(C.hydrated(rk_path), comment="#"); fam = pd.read_csv(C.hydrated(fam_path)).set_index("feature")["family"]
    rk["family"] = fam.loc[rk.feature].values; rk["g"] = rk.dC * 1000.0; rki = rk.set_index("feature"); raw = set(rk.feature)
    pairs = [(f, f[:-6]) for f in rk.feature if f.endswith("_siteZ") and f[:-6] in raw]
    if pairs: assert max(abs(float(rki.loc[f, "dC"]) - float(rki.loc[b, "dC"])) for f, b in pairs) < 1e-5
    distinct = rk[~rk.feature.isin([f for f, _ in pairs])]
    tail = distinct[distinct.g > 3.0].sort_values("g", ascending=False).reset_index(drop=True)
    return rk, tail
NICE = dict(GB["nice_all"]); NICE.update({r["feature"]: r["name"] for r in GB["rows"]})
RK_NEW, TAIL_NEW = rows_from(f"{C.NUM}/ranking_v3.csv", f"{C.NUM}/measure_families.csv")
RES_NEW = json.load(open(C.hydrated(f"{C.NUM}/results_v2.json")))
BAND_KEYS = ["0-1%", "1-5%", "5-10%", ">10%"]

# ---------------------------------------------------------------- 1 page box, 2 letters, 3 fonts
N_ROWS = GB["n_rows"]
check(abs(nrect[0] - (brect[0] + XB)) < 0.01 and abs(nrect[1] - (677.81 + GROW + DY)) < 0.02 and abs(brect[1] - 717.92) < 0.02 and abs(GROW - (N_ROWS - 28) * PITCH) < 1e-6,
      f"page box: width = V13 + {XB} pt (round 43), height = the 28-row base plus {N_ROWS - 28} row pitches plus the {DY:g} pt strip (round 40)", f"{brect} -> {nrect}")
lets_ok = [l[0] for l in blets] == [l[0] for l in nlets] == ["a", "b"]
if lets_ok:
    for b, n in zip(blets, nlets):
        dx = XB if b[0] == "b" else 0.0
        lets_ok &= abs(b[1] + dx - n[1]) <= 0.3 and abs(b[2] + DY - n[2]) <= 0.3 + 1.0 and V.face(b[4]) == V.face(n[4]) and abs(n[5] - (b[5] + DELTA)) < 0.05
check(lets_ok, f"letters a, b: same face, {13 + DELTA:g} pt (V13 13), a at its V13 x, b at x + {XB}, both {DY:g} pt lower (the bbox top within 1.3 pt: the taller glyph)", f"V13 {[(l[0], l[1], l[2], l[5]) for l in blets]} new {[(l[0], l[1], l[2], l[5]) for l in nlets]}")
bf = {(V.face(s["font"]), round(s["size"] + DELTA, 2)) for s in bsp}; nf = {(V.face(s["font"]), round(s["size"], 2)) for s in nsp}
check(bf == nf, f"font faces identical to V13 and every size + {DELTA:g} pt", f"V13 + {DELTA:g}: {sorted(bf)} new: {sorted(nf)}")

# ---------------------------------------------------------------- 4 rows: the build's row list equals the independent re-derivation
build_rows = [(r["feature"], r["printed"], r["family"]) for r in GB["rows"]]
rederived = [(r.feature, f"{r.g:.1f}", r.family) for _, r in TAIL_NEW.iterrows()]
check(build_rows == rederived and len(rederived) == N_ROWS, f"panel b row set, order, printed gains and families equal an independent re-derivation ({N_ROWS} rows)", f"{len(build_rows)} rows; first mismatch {next((i for i, (a, b) in enumerate(zip(build_rows, rederived)) if a != b), None)}")
band_ok = all(int(RES_NEW["cohort"]["t90_bands_n"][k]) == b["n"] and abs(float(RES_NEW["cohort"]["t90_bands_pct"][k]) - b["pct"]) < 0.051 for k, b in zip(BAND_KEYS, GA["bands"]))
check(band_ok and sum(b["n"] for b in GA["bands"]) == RES_NEW["cohort"]["n"], "panel a band counts and shares equal results_v2.json (sum = cohort n)", str([(b["label"], b["n"], round(b["pct"], 2)) for b in GA["bands"]]))

# ---------------------------------------------------------------- 5 geometry pins on the new sheet (V13 pins + DY, panel b + XB)
G = V.geometry_records(NEW_PDF); INK = "#1a1d21"
pins = [("a left spine", 79.59, 33.122 + DY, 79.59, 263.522 + DY, INK), ("a bottom spine", 79.59, 263.522 + DY, 147.99, 263.522 + DY, INK)]
for v in (0, 20, 40, 60): pins.append((f"a x tick {v}", 79.59 + v * GA["ppu"], 263.522 + DY, 79.59 + v * GA["ppu"], 266.522 + DY, INK))
for k in range(4):
    cy = 33.122 + (3.72 - (3 - k)) * GA["pitch"] + DY; pins.append((f"a y tick band {k}", 76.59, cy, 79.59, cy, INK))
yb = GB["y1"] + DY; pins.append(("b bottom spine", 292.59 + XB, yb, 371.05 + XB, yb, INK))
for v in (5, 10, 15, 20):
    x = 292.59 + XB + (v - 2.85) * (371.05 - 292.59) / (24.5 - 2.85); pins.append((f"b x tick {v}", x, yb, x, yb + 3.0, INK))
missing = [p[0] for p in pins if not V.has_segment(G, *p[1:5], stroke=p[5])]
check(not missing and abs(yb - (633.17 + GROW + DY)) < 0.02, f"static geometry pins ({len(pins)} spines and ticks at the V13 positions + {DY:g} pt, panel b's + {XB} pt)", f"missing {missing}")
bar_miss = [b["label"] for v, b in zip(DA["values"], GA["bands"]) if v["kind"] == "drawn" and not V.has_rect(G, v["bar_rect"][0], v["bar_rect"][1] + DY, v["bar_rect"][2], v["bar_rect"][3] + DY, fill=v["colour"], tol=0.3)]
check(not bar_miss, "panel a bars: four rectangles at x0 = 79.59 with width = share x 1.14 pt and the T90 ramp colours (data marks unchanged)", f"missing {bar_miss}")
dot_miss, lead_miss = [], []
for v in [v for v in DB["values"] if v["kind"] == "printed"]:
    cx, cy = v["dot_page_xy"]; cy += DY
    if not V.has_dot(G, cx, cy, v["colour"]): dot_miss.append((v["name"], cx, cy))
    if not V.has_segment(G, 293.134 + XB, cy, cx, cy, stroke="#b3dcf2", tol=0.35): lead_miss.append((v["name"], cx, cy))
check(not dot_miss and not lead_miss, f"panel b: {N_ROWS} dots at (gain, row) in the family colours and {N_ROWS} dotted leaders from 3.0 to the dot (data marks unchanged, + {XB} pt with the panel)", f"dots missing {dot_miss[:5]} leaders missing {lead_miss[:5]}")

# ---------------------------------------------------------------- 6 text layer: V13 tokens minus its strings == new tokens minus the new strings (words unchanged by any re-wrap)
OLD_FAM = f"{C.OLDNUM}/numbers/measure_families.csv" if os.path.exists(f"{C.OLDNUM}/numbers/measure_families.csv") else f"{C.NUM}/measure_families_PRE_V8.csv"
RK_OLD, TAIL_OLD = rows_from(f"{C.OLDNUM}/numbers/ranking_v3.csv", OLD_FAM)
OLD_NICE = dict(NICE); OLD_NICE.setdefault("AHI", "Apnea-hypopnea index")
old_strings = [OLD_NICE[r.feature] for _, r in TAIL_OLD.iterrows()] + [f"{r.g:.1f}" for _, r in TAIL_OLD.iterrows()] + [r.family for _, r in TAIL_OLD.iterrows()]
new_strings = [r["name"] for r in GB["rows"]] + [r["printed"] for r in GB["rows"]] + [r["family"] for r in GB["rows"]] + [STAMP["title"]]
res = V.residue_check(check, f"text layer: V13 tokens minus the v7 strings equal new tokens minus the v8.1 strings and the title '{STAMP['title']}' (static words identical: the two-line axis title and the re-wrapped names keep their words)", bsp, nsp, old_strings, new_strings)
dup = [k for k, v in __import__("collections").Counter((s["text"], s["origin"][0], s["origin"][1]) for s in nsp).items() if v > 1]
check(not dup, "no duplicated spans on the new sheet", str(dup[:5]))

# ---------------------------------------------------------------- 7 read-back of every placed string (builder coordinates + DY)
def shifted(values): return [dict(v, baseline=v["baseline"] + DY) for v in values]
V.readback(check, "read-back panel a (+16 pt)", nsp, shifted(DA["values"])); V.readback(check, "read-back panel b (+16 pt)", nsp, shifted(DB["values"]))
second = []
for r in GB["rows"]:
    if len(r["lines"]) == 2:
        y = GB["y_row0"] + PITCH * (r["row"] - 1) + r["two_dy"][1] + DY
        s, d = V.find_span(nsp, r["lines"][1], GB["name_x"] + XB, y, ha="left", tol=0.6)
        if s is None: second.append((r["lines"][1], y, d))
check(not second, f"second lines of the {sum(1 for r in GB['rows'] if len(r['lines']) == 2)} two-line names at row centre + their pitch offset ({GB['two_dy'][1]} or {GB['two_dy_wide'][1]:.2f} pt)", f"fails {second}")
tl = [s for s in nsp if s["text"].replace("\xa0", " ") == STAMP["title"]]
check(len(tl) == 1 and abs(tl[0]["origin"][0] - STAMP["origin"][0]) < 0.05 and abs(tl[0]["origin"][1] - STAMP["origin"][1]) < 0.05 and abs(tl[0]["size"] - STAMP["size"]) < 0.01 and "Bold" in tl[0]["font"],
      f"title '{STAMP['title']}' once, Arial Bold {STAMP['size']:g} pt at ({STAMP['origin'][0]}, {STAMP['origin'][1]})", str([(t["text"], t["font"], t["size"], t["origin"]) for t in tl]))
xl = {s["text"]: s for s in nsp if s["text"] in ("Share of the", "cohort, %")}
xl_ok = len(xl) == 2 and all(abs((xl[t]["bbox"][0] + xl[t]["bbox"][2]) / 2 - GA["xlabel_xc"]) < 0.6 and abs(xl[t]["origin"][1] - (y + DY)) < 0.3 and abs(xl[t]["size"] - (11.0 + DELTA)) < 0.01 for t, y in GA["xlabel_lines"])
check(xl_ok, f"panel a x-axis title as two lines centred on the axis (x {GA['xlabel_xc']:.2f}), baselines {[round(y + DY, 2) for _, y in GA['xlabel_lines']]}, {11 + DELTA:g} pt (round 40 re-applied)", str([(t, s["origin"], s["size"]) for t, s in xl.items()]))
present = {s["text"] for s in nsp}
miss = sorted(t for t in new_strings if t not in present and "\n" not in t and t in NICE.values() and all(len(r["lines"]) == 1 for r in GB["rows"] if r["name"] == t))
check(not miss, "every one-line name, printed gain and family string of the re-derived row set is on the sheet", f"missing {miss}")
bad = C.banned_words([s["text"] for s in nsp]); check(not bad, "no banned words on the sheet", str(bad))
# every name that the wrap rule re-splits at 11 pt is declared in the builder record, and no one-line name at 11 pt exceeds the base's one-line limit
wc = {w["name"] for w in GB["wrap_changed_at_plus1"]}; resplit = {r["name"] for r in GB["rows"] if r["lines"] != r["lines_v26"]}
check(wc == resplit, f"re-wrapped names declared ({len(wc)}): {sorted(wc)}", "")
name_w = [(s["text"], round(s["bbox"][2] - s["bbox"][0], 2)) for s in nsp if abs(s["origin"][0] - (GB["name_x"] + XB)) < 0.3 and abs(s["size"] - (10 + DELTA)) < 0.01]
over = [(t, w) for t, w in name_w if w > GB["name_max_pt"] + 1.0 + 0.5]
check(not over and name_w, f"every name line at {10 + DELTA:g} pt within the base's one-line limit ({GB['name_max_pt'] + 1.0:.2f} pt, +0.5 measurement slack): widest {max(w for _, w in name_w):.2f} pt of {len(name_w)} lines, leaders start {293.134 + XB - (GB['name_x'] + XB) - max(w for _, w in name_w):.2f} pt right of the widest", f"over {over}")

# ---------------------------------------------------------------- crops (200 dpi, V13 base vs new, Ghostscript renders through the copied helper)
H_new = nrect[1]
C.crop_pair(BASE_PDF, NEW_PDF, (0, 0, 175, 340), f"{CROPS}/{S}_a_v13_vs_new_200dpi.png", dpi=200, label_old="V13 base", label_new="round 49 (+1 pt)")
C.crop_pair(BASE_PDF, NEW_PDF, (165, 0, 533.14, 380), f"{CROPS}/{S}_b_top_v13_vs_new_200dpi.png", dpi=200, label_old="V13 base", label_new="round 49 (+1 pt)")
C.crop_pair(BASE_PDF, NEW_PDF, (165, 340, 533.14, H_new), f"{CROPS}/{S}_b_bottom_v13_vs_new_200dpi.png", dpi=200, label_old="V13 base", label_new="round 49 (+1 pt)")
finish()
