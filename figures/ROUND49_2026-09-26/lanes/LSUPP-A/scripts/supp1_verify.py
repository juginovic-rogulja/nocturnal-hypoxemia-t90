#!/usr/bin/env python3
"""03_verify for Supp_Fig01 (lane V14_L1_RANK, round 37: v8.1 numbers against the V13 sheet; Proof W). Page box, fonts, no letters, geometry pins on every row (row axes, ticks, the drawn 5th to 95th line, the
interquartile box and the median tick at the positions the NEW file puts them, and the BASE sheet's rows at the positions the OLD file puts
them, which proves the base read back to the old numbers), the text-layer proof (only the sleep-efficiency ticks change), the read-back of
every placed string, banned words, crops, the 150 dpi render and verify/CHANGES_Supp_Fig01.csv."""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C, vcommon as V

S = "Supp_Fig01"; SD = f"{C.LANE}/{S}"; WORK = f"{SD}/work"; VER = f"{SD}/verify"; CROPS = f"{VER}/crops"; os.makedirs(CROPS, exist_ok=True)
BASE_PDF = C.hydrated(f"{C.BASE}/{S}.pdf"); NEW_PDF = C.hydrated(f"{SD}/{S}.pdf")
check, finish = C.checks_writer(f"{VER}/checks.txt")
bsp, brect, blets = V.sheet_spans(BASE_PDF); nsp, nrect, nlets = V.sheet_spans(NEW_PDF)
D = json.load(open(f"{VER}/{S}_rows_drawn.json")); G = D["geometry"]
PROF = json.load(open(C.hydrated(f"{C.NUM}/figure1_sleep_profile_v1.json"))); OPROF = json.load(open(C.hydrated(f"{C.OLDNUM}/numbers/figure1_sleep_profile_v1.json")))
AX_L, AX_R = G["ax_l"], G["ax_r"]
def px(axis, v): lo, _, hi = axis; return AX_L + (AX_R - AX_L) * (v - lo) / (hi - lo)

check(brect == nrect, "page box equal", f"{brect} vs {nrect}")
check(blets == [] and nlets == [], "no panel letters on either sheet", f"{blets} {nlets}")
V.fonts_check(check, bsp, nsp)

# ---------------------------------------------------------------- the build's values equal the file (independent read)
val_ok = all(all(abs(float(PROF["measures"][r["key"]][q]) - float(r["values"][q])) < 1e-9 for q in ("median", "q1", "q3", "p05", "p95", "n")) for r in G["rows"])
axis_ok = all(r["axis"][0] <= r["values"]["p05"] and r["values"]["p95"] <= r["axis"][2] for r in G["rows"])
check(val_ok and axis_ok and len(G["rows"]) == 9, "nine rows: drawn quantiles equal figure1_sleep_profile_v1.json, every axis contains the 5th to 95th percentile", "")
check(all(r["rederived"][q] == r["values"][q] for r in G["rows"] for q in ("median", "q1", "q3", "p05", "p95", "n")), "every quantile and n re-derived from the v8 frozen table at the stored precision (sleep rows over the full nights)", "")

# ---------------------------------------------------------------- geometry pins, new sheet against the new file, base sheet against the old file
def row_pins(recs, y, axis, m, col, label):
    miss = []
    if not V.has_segment(recs, AX_L, y, AX_R, y, stroke="#ccd1d6", tol=0.3): miss.append("row axis")
    for v in axis:
        if not V.has_segment(recs, px(axis, v), y, px(axis, v), y + G["tick_len"], stroke="#ccd1d6", tol=0.3): miss.append(f"tick {v}")
    if not V.has_segment(recs, px(axis, m["p05"]), y, px(axis, m["p95"]), y, stroke=col, tol=0.3): miss.append("p05 to p95 line")
    for q in ("p05", "p95"):
        if not V.has_segment(recs, px(axis, m[q]), y - G["cap_h"], px(axis, m[q]), y + G["cap_h"], stroke=col, tol=0.3): miss.append(f"{q} cap")
    if not V.has_rect(recs, px(axis, m["q1"]), y - G["bar_h"] / 2, px(axis, m["q3"]), y + G["bar_h"] / 2, fill=col, tol=0.3): miss.append("interquartile box")
    if not V.has_segment(recs, px(axis, m["median"]), y - G["bar_h"] / 2, px(axis, m["median"]), y + G["bar_h"] / 2, stroke="#ffffff", tol=0.3): miss.append("median tick")
    return [(label, x) for x in miss]
GN = V.geometry_records(NEW_PDF); GBASE = V.geometry_records(BASE_PDF)
miss_new, miss_base = [], []
for r in G["rows"]:
    col = r["colour"]; miss_new += row_pins(GN, r["row_y"], r["axis"], r["values"], col, r["label"])
    miss_base += row_pins(GBASE, r["row_y"], r["axis_old"], OPROF["measures"][r["key"]], col, r["label"])
check(not miss_new, "new sheet: every row's axis, ticks, 5th to 95th line with caps, interquartile box and median tick sit where the NEW file puts them (9 rows x 9 elements)", f"missing {miss_new[:8]}")
check(not miss_base, "base sheet: the same nine rows sit where the OLD file puts them on the old axes (the base reads back to the pre-v7 numbers)", f"missing {miss_base[:8]}")
key_ok = all(V.has_rect(GN, *box, fill=col, tol=0.3) for _, col, box, _ in [(k[0], k[1], k[2], k[3]) for k in G["key"]])
check(key_ok, "key squares at the base positions in #0288d1 and #8a9099", "")

# ---------------------------------------------------------------- text layer and read-back
old_strings = [f"{v:g}" for r in G["rows"] for v in r["axis_old"] if tuple(r["axis"]) != tuple(r["axis_old"])]; new_strings = [f"{v:g}" for r in G["rows"] for v in r["axis"] if tuple(r["axis"]) != tuple(r["axis_old"])]
res = V.residue_check(check, f"text layer: base minus the V13 ticks of the re-chosen axes {old_strings} equals new minus the new ticks {new_strings} (every other string identical)", bsp, nsp, old_strings, new_strings)
print("  numbers lost/gained:", res["numbers_lost_gained"])
V.readback(check, "read-back", nsp, D["values"])
bad = C.banned_words([s["text"] for s in nsp]); check(not bad, "no banned words on the sheet", str(bad))

# ---------------------------------------------------------------- crops, render, CHANGES
C.crop_pair(BASE_PDF, NEW_PDF, (0, 0, 484.72, 369.36), f"{CROPS}/{S}_whole_old_vs_new_150dpi.png", dpi=150)
C.crop_pair(BASE_PDF, NEW_PDF, (0, 20, 484.72, 240), f"{CROPS}/{S}_rows1to6_old_vs_new_200dpi.png", dpi=200)
C.crop_pair(BASE_PDF, NEW_PDF, (0, 240, 484.72, 369.36), f"{CROPS}/{S}_oxygen_rows_old_vs_new_200dpi.png", dpi=200)
C.render_png(NEW_PDF, f"{SD}/{S}_150dpi.png", dpi=150)
rows = []
for r in G["rows"]:
    o = OPROF["measures"][r["key"]]; n = r["values"]
    fmt = lambda m: f"median {m['median']:g} (IQR {m['q1']:g} to {m['q3']:g}), 5th to 95th {m['p05']:g} to {m['p95']:g}, n {m['n']:,}"
    rows.append(dict(sheet=S, panel="rows", item=f"{r['label']} (drawn)", old_value=fmt(o), new_value=fmt(n), source_file="numbers/figure1_sleep_profile_v1.json", key_path=f"measures.{r['key']}",
                     dramatic="YES" if abs(n["n"] - o["n"]) / o["n"] > 0.2 else "", note=f"median {100 * (n['median'] - o['median']) / o['median']:+.1f}%, n {n['n'] - o['n']:+,}"))
    if tuple(r["axis"]) != tuple(r["axis_old"]):
        rows.append(dict(sheet=S, panel="rows", item=f"{r['label']} axis ticks (printed)", old_value=" ".join(f"{v:g}" for v in r["axis_old"]), new_value=" ".join(f"{v:g}" for v in r["axis"]),
                         source_file="numbers/figure1_sleep_profile_v1.json", key_path=f"measures.{r['key']}.p05 = {n['p05']} (below the old axis floor of {r['axis_old'][0]})", dramatic="notable", note="axis re-chosen at round values that contain the 5th to 95th percentile"))
C.write_changes(f"{VER}/CHANGES_{S}.csv", rows)
print(f"  CHANGES rows: {len(rows)}, dramatic YES: {sum(1 for r in rows if r['dramatic'] == 'YES')}, notable: {sum(1 for r in rows if r['dramatic'] == 'notable')}")
finish()
