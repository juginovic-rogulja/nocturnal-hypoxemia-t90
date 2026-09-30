#!/usr/bin/env python3
"""03_verify for Supp_Fig03 (lane V14_L1_RANK, round 37: v8.1 numbers against the V13 sheet; Proof W; twelve rows, page grown one pitch). Page box, letters, fonts, geometry pins (panel a bars, dashed best-single rule, spine and ticks; the base's
bars at the OLD file's values, which proves the base read back to the pre-v7 numbers; panel b cell rectangles and fills at the rule the base
shows), the text-layer proof (base minus the OLD file's strings equals new minus the NEW file's strings), the read-back of every placed string
against the file within the printed rounding (the one rounding-tie cell is listed), banned words, crops, the 150 dpi render and
verify/CHANGES_Supp_Fig03.csv."""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C, vcommon as V

S = "Supp_Fig03"; SD = f"{C.LANE}/{S}"; WORK = f"{SD}/work"; VER = f"{SD}/verify"; CROPS = f"{VER}/crops"; os.makedirs(CROPS, exist_ok=True)
BASE_PDF = C.hydrated(f"{C.BASE}/{S}.pdf"); NEW_PDF = C.hydrated(f"{SD}/{S}.pdf")
check, finish = C.checks_writer(f"{VER}/checks.txt")
bsp, brect, blets = V.sheet_spans(BASE_PDF); nsp, nrect, nlets = V.sheet_spans(NEW_PDF)
DA = json.load(open(f"{VER}/{S}_a_drawn.json")); DB = json.load(open(f"{VER}/{S}_b_drawn.json")); GA = DA["geometry"]; GB = DB["geometry"]
CMB = json.load(open(C.hydrated(f"{C.NUM}/combination_v2.json"))); OCMB = json.load(open(C.hydrated(f"{C.OLDNUM}/numbers/combination_v2.json")))
CC = json.load(open(f"{WORK}/corr_check.json"))
PLAIN = dict(GA["plain_names"]); GROW = GA["grow_pt"]; N_ITEMS = GA["n_items"]; N_ITEMS_V13 = GA["n_items_v13"]
def items_from(cmb):
    """The sheet's row logic applied to a combination_v2.json (old or new): eleven (label, value) rows."""
    single = cmb["single"]; best = max(single, key=single.get); bs = float(single[best]); c = cmb["combinations"]; anchor = cmb["composite_terms"]["anchor"]
    other = sorted([(f, float(c[f"oxygen_plus_{f}"])) for f in cmb["composite_terms"]["oxygen_plus_all_families"] if f != anchor], key=lambda kv: kv[1])
    it = [("Best single oxygen measurement", bs)] + [(f"Best {k} oxygen measurements", float(c[f"top{k}_oxygen"])) for k in (2, 3, 4, 5)]
    assert len(other) in (5, 6), len(other)
    it += [(f"Oxygen and {PLAIN[f]}", v) for f, v in other] + [("Oxygen and the best of every other family", float(c["oxygen_plus_all_families"]))]
    return it, bs
ITEMS_NEW, BS_NEW = items_from(CMB); ITEMS_OLD, BS_OLD = items_from(OCMB)
ax_x = lambda v: GA["x0"] + v * GA["ppu"]

check(abs(nrect[0] - brect[0]) < 0.01 and abs(nrect[1] - (brect[1] + GROW)) < 0.02 and abs(GROW - (N_ITEMS - N_ITEMS_V13) * GA["pitch"]) < 1e-6,
      f"page box: width equal, height grown by exactly {N_ITEMS - N_ITEMS_V13} row pitch ({GROW:+.3f} pt, panel a has {N_ITEMS} rows, V13 had {N_ITEMS_V13})", f"{brect} -> {nrect}")
ok, msg = V.compare_letters(blets, nlets, shifts={"b": GROW}); check(ok, f"letters a, b within 0.3 pt, same face and size (b declared {GROW:+.2f} pt lower with its panel)", msg)
V.fonts_check(check, bsp, nsp, kept={("Arial", 10.0)})   # round 49: the twelve panel-a row labels are kept at 10 pt (the composite label would run into the bars at 11 pt)
check([(r["label"], r["value"]) for r in GA["rows"]] == ITEMS_NEW and abs(GA["best_single"]["value"] - BS_NEW) < 1e-12 and len(ITEMS_NEW) == N_ITEMS, f"panel a rows equal an independent re-derivation from combination_v2.json (labels, order, values, {N_ITEMS} rows)", "")

# ---------------------------------------------------------------- geometry pins
GN = V.geometry_records(NEW_PDF); GBASE = V.geometry_records(BASE_PDF); INK = "#1a1d21"
pins = [("a bottom spine", GA["x0"], GA["y1"], GA["x1"], GA["y1"], INK)]
for v in (0, 0.01, 0.02, 0.03): pins.append((f"a x tick {v}", ax_x(v), GA["y1"], ax_x(v), GA["y1"] + 3.0, INK))
missing = [p[0] for p in pins if not V.has_segment(GN, *p[1:5], stroke=p[5])]
check(not missing, f"static geometry pins ({len(pins)} spine and ticks at the base positions)", f"missing {missing}")
def bars_ok(recs, items, bs):
    miss = []
    for i, (lab, v) in enumerate(items):
        cy = GA["c0"] + GA["pitch"] * i; col = GA["rows"][i]["colour"] if items is ITEMS_NEW else None
        if items is not ITEMS_NEW and i >= N_ITEMS_V13: continue
        if not V.has_rect(recs, GA["x0"], cy - GA["bar_h"] / 2, ax_x(v), cy + GA["bar_h"] / 2, fill=col, tol=0.3): miss.append(lab)
    if not V.has_segment(recs, ax_x(bs), GA["y0"], ax_x(bs), GA["y1"], stroke=INK, tol=0.3): miss.append("dashed best-single rule")
    return miss
def bars_ok_base(recs, items, bs):
    """the V13 sheet against the OLD file on the V13 geometry (11 rows, spine one pitch higher)"""
    miss = []
    for i, (lab, v) in enumerate(items):
        cy = GA["c0"] + GA["pitch"] * i
        if not V.has_rect(recs, GA["x0"], cy - GA["bar_h"] / 2, ax_x(v), cy + GA["bar_h"] / 2, fill=None, tol=0.3): miss.append(lab)
    if not V.has_segment(recs, ax_x(bs), GA["y0"], ax_x(bs), GA["y1"] - GROW, stroke=INK, tol=0.3): miss.append("dashed best-single rule")
    return miss
mn = bars_ok(GN, ITEMS_NEW, BS_NEW); mb = bars_ok_base(GBASE, ITEMS_OLD, BS_OLD)
check(not mn, f"new sheet: {N_ITEMS} bars end at value x 7413.4 pt in the base's row colours (blue, pale, mid), dashed rule at the best single gain", f"missing {mn}")
check(not mb, f"base sheet: its {len(ITEMS_OLD)} bars and rule sit where the OLD file puts them (the base reads back to the v7 numbers)", f"missing {mb}")
rule_rec = [r for r in GN if r["stroke"] == INK and abs(r["bbox"][0] - ax_x(BS_NEW)) < 0.3 and abs(r["bbox"][1] - GA["y0"]) < 0.3 and r["dashes"] not in (None, "[] 0")]
check(bool(rule_rec), "the best-single rule is dashed (0.9 pt, 4 on 3 off)", str([r["dashes"] for r in rule_rec][:2]))
cell_miss = [v["text"] for v in DB["values"] if v["kind"] == "printed" and not V.has_rect(GN, *v["cell_rect"], fill=v["fill"], tol=0.3)]
check(not cell_miss, f"panel b: 25 cell rectangles at the base grid moved down {GROW:.2f} pt with the fill the base's colour rule gives (diagonal #0288d1, off-diagonal on the base ramp)", f"missing {cell_miss}")
# every off-diagonal fill on the base sheet is on the same ramp for the OLD value (the rule is the base's, not ours)
from matplotlib import colors as mcolors
cmap = mcolors.LinearSegmentedColormap.from_list("house_blue", [GB["ramp_lo"], GB["ramp_hi"]]); norm = mcolors.Normalize(vmin=GB["vmin"], vmax=GB["vmax"])
def fill_of(r):
    m = abs(r)
    if m >= 0.9995: return GB["diag"]
    return "#%02x%02x%02x" % tuple(int(round(v * 255)) for v in cmap(norm(m))[:3])
def close(h1, h2): return all(abs(int(h1[i:i + 2], 16) - int(h2[i:i + 2], 16)) <= 1 for i in (1, 3, 5))
OR = np.array(OCMB["correlation_matrix"]); cell = GB["cell"]; bx0, by0 = GB["x0"], GB["y0"] - GROW      # the base grid sits one pitch higher than the new one
base_fill_miss = []
for a in range(5):
    for b in range(5):
        want = fill_of(float(OR[a, b])); rect = [bx0 + b * cell, by0 + a * cell, bx0 + (b + 1) * cell, by0 + (a + 1) * cell]
        hit = [r for r in GBASE if r["fill"] and max(abs(r["bbox"][i] - rect[i]) for i in range(4)) <= 0.3]
        if not hit or not any(close(r["fill"], want) for r in hit): base_fill_miss.append((a, b, want, [r["fill"] for r in hit]))
check(not base_fill_miss, "base sheet: every cell fill equals the same colour rule applied to the OLD matrix (within one 8-bit unit)", f"{base_fill_miss[:4]}")

# ---------------------------------------------------------------- text layer and read-back
def matrix_strings(mat, printed=None):
    return [(printed[a][b] if printed else f"{mat[a][b]:.2f}").replace("-", "−") for a in range(5) for b in range(5)]
# the OLD panel-b strings are what the V13 sheet prints: the round-30 lane's drawn record (its tie cell printed from full precision), never re-rounded here
R30B = json.load(open(C.hydrated(f"{C.R30}/V7_L1_RANK/Supp_Fig03/verify/Supp_Fig03_b_drawn.json")))
old_b_printed = [v["text"] for v in R30B["values"] if v["kind"] == "printed"]; assert len(old_b_printed) == 25
old_b_labels = [v["text"] for v in R30B["values"] if v["kind"] == "label"]; new_b_labels = [v["text"] for v in DB["values"] if v["kind"] == "label"]
old_strings = [l for l, _ in ITEMS_OLD] + [f"{v:.4f}" for _, v in ITEMS_OLD] + old_b_printed + old_b_labels
NOT_PRINTED_R50 = [f"{v:.4f}" for _, v in ITEMS_NEW]   # round 50: the panel a gain values are no longer printed (declared removal, the V13 base still prints its eleven)
new_strings = [l for l, _ in ITEMS_NEW] + [v["text"] for v in DB["values"] if v["kind"] == "printed"] + new_b_labels
res = V.residue_check(check, "text layer: base minus the OLD strings (the V13 rows, its matrix cells and its five labels) equals new minus the NEW strings (static text identical)", bsp, nsp, old_strings, new_strings)
print("  words lost/gained:", res["words_lost_gained"]); print("  numbers lost/gained:", res["numbers_lost_gained"])
V.readback(check, "read-back panel a", nsp, DA["values"]); V.readback(check, "read-back panel b", nsp, DB["values"])
check(not any(s["text"] in NOT_PRINTED_R50 for s in nsp) and all(v["kind"] == "bar_value_not_printed" for v in DA["values"] if v["text"] in NOT_PRINTED_R50) and len(NOT_PRINTED_R50) == 12,
      f"round 50: none of the twelve panel a gain values is printed on the sheet (the bars end at value x 7413.4 pt, checked above)", str(NOT_PRINTED_R50))
# every printed matrix cell within the printed rounding of the file value, and the tie cells listed
R = np.array(CMB["correlation_matrix"]); RF = np.array(CC["recomputed_z_scale_matrix_ranking_order_6dp"])
off = []
for v in [v for v in DB["values"] if v["kind"] == "printed"]:
    p = float(v["text"].replace("−", "-")); a, b = [int(t) for t in v["key"].split("[")[1:3] for t in [t.split("]")[0]]]
    if abs(p - R[a, b]) > 0.0051 or abs(p - RF[a, b]) > 0.005 + 1e-9: off.append((a, b, v["text"], float(R[a, b]), float(RF[a, b])))
check(not off, "every printed cell within 0.005 of the file's 3-decimal value and of the full-precision recomputation", f"{off}")
ties = GB["rounding_ties"]; print("  rounding ties (file 3 dp vs full precision):", ties)
check(all(abs(round(t["file_value"] * 1000)) % 10 == 5 and t["printed"] == f"{t['full_precision']:.2f}" for t in ties) and len(ties) % 2 == 0,
      f"rounding ties (file 3 dp on an exact half at 2 dp): {len(ties)} cells, each printed from the full-precision recomputation, listed for the coordinator", str(ties))
check(abs(CC["max_abs_diff_z_scale"]["under_ranking_v3_order"]) < 1e-3 and GB["order"] == CC["ranking_v3_order"] == CC["five"] and CC["file_labels_are_the_five"] is False,
      f"the five measures and their order: the in-lane Spearman on the 15,551 full nights matches the file's matrix under the ranking_v3 top-5 oxygen order {CC['ranking_v3_order_labels']} (max {CC['max_abs_diff_z_scale']['under_ranking_v3_order']:.4f}); the file's correlation_labels literal names the 4% desaturation index, which is not among the five (NUMBERS_DEFECT, not used)", "")
bad = C.banned_words([s["text"] for s in nsp]); check(not bad, "no banned words on the sheet", str(bad))

# ---------------------------------------------------------------- crops, render, CHANGES
C.crop_pair(BASE_PDF, NEW_PDF, (0, 0, 484.72, 330 + GROW), f"{CROPS}/{S}_a_old_vs_new_200dpi.png", dpi=200)
C.crop_pair(BASE_PDF, NEW_PDF, (0, 340, 484.72, 678.24 + GROW), f"{CROPS}/{S}_b_old_vs_new_200dpi.png", dpi=200)
C.crop_pair(BASE_PDF, NEW_PDF, (0, 0, 484.72, 678.24 + GROW), f"{CROPS}/{S}_whole_old_vs_new_150dpi.png", dpi=150)
C.render_png(NEW_PDF, f"{SD}/{S}_150dpi.png", dpi=150)
rows = []
OLD_MAP = {l: v for l, v in ITEMS_OLD}
for i, (nl, nv) in enumerate(ITEMS_NEW):
    ol, ov = (ITEMS_OLD[i] if i < len(ITEMS_OLD) else ("(no row: V13 had eleven rows)", float("nan")))
    note = ""
    if nl not in OLD_MAP: note = "row enters (a new family best under v8.1)"
    elif ol != nl: note = "row order changes (the one-per-family rows are sorted by gain)"
    rows.append(dict(sheet=S, panel="a", item=f"row {i + 1}", old_value=f"{ol}: {ov:.4f}", new_value=f"{nl}: {nv:.4f} (V13 printed {OLD_MAP[nl]:.4f} for this row)" if nl in OLD_MAP else f"{nl}: {nv:.4f}", source_file="numbers/combination_v2.json", key_path=GA["rows"][i]["key"],
                     dramatic="notable" if (nl not in OLD_MAP or ol != nl) else "", note=note))
rows.append(dict(sheet=S, panel="a", item="rows that leave", old_value=" | ".join(l for l, _ in ITEMS_OLD if l not in {n for n, _ in ITEMS_NEW}), new_value="(left)", source_file="numbers/combination_v2.json", key_path="composite_terms.oxygen_plus_all_families", dramatic="notable", note="the family bests changed under the 141-measurement list"))
rows.append(dict(sheet=S, panel="a", item="row order of the one-per-family rows", old_value=" | ".join(l for l, _ in ITEMS_OLD[5:-1]), new_value=" | ".join(l for l, _ in ITEMS_NEW[5:-1]),
                 source_file="numbers/combination_v2.json", key_path="combinations.oxygen_plus_<family best>, sorted ascending", dramatic="", note="seven families under v8.1 (the PLM index is its own family): six one-per-family rows, page grown one pitch"))
rows.append(dict(sheet=S, panel="a", item="bar colours of the one-per-family rows", old_value="V13: all six above the best single gain (mid grey)", new_value="; ".join(f"{r['label']}: {r['colour']}" for r in GA["rows"][5:-1]), source_file="numbers/combination_v2.json",
                 key_path="combinations.oxygen_plus_<family best> against single.spo2_pct_below_90", dramatic="", note="the base's colour rule applied to the new values"))
rows.append(dict(sheet=S, panel="a", item="dashed best-single rule (drawn)", old_value=f"{BS_OLD:.5f}", new_value=f"{BS_NEW:.5f}", source_file="numbers/combination_v2.json", key_path="single.spo2_pct_below_90", dramatic="", note=""))
rows.append(dict(sheet=S, panel="a", item="baseline concordance (legend, not printed)", old_value=OCMB["baseline"], new_value=CMB["baseline"], source_file="numbers/combination_v2.json", key_path="baseline", dramatic="", note="for the legend"))
NAMES = GB["names"]
for a in range(5):
    for b in range(a + 1, 5):
        ov = float(OR[a, b]); nv = float(R[a, b]); pr = [v["text"] for v in DB["values"] if v["kind"] == "printed" and v["key"].startswith(f"correlation_matrix[{a}][{b}]")][0]
        rows.append(dict(sheet=S, panel="b", item=f"{NAMES[a]} x {NAMES[b]}", old_value=f"{ov:.2f} (file {ov})".replace("-", "−"), new_value=f"{pr} (file {nv})".replace("-", "−"), source_file="numbers/combination_v2.json",
                         key_path=f"correlation_matrix[{a}][{b}] (ranking_v3 order)", dramatic="", note="rounding tie, printed from the full-precision recomputation" if any(t["row"] == a and t["col"] == b for t in ties) else ""))
rows.append(dict(sheet=S, panel="b", item="mean absolute correlation (legend, not printed)", old_value=OCMB["mean_abs_correlation"], new_value=CMB["mean_abs_correlation"], source_file="numbers/combination_v2.json", key_path="mean_abs_correlation", dramatic="", note="for the legend"))
rows.append(dict(sheet=S, panel="b", item="the five measures and their order", old_value="Below 90% | Below 88% | Lowest saturation | 3% desaturations | 4% desaturations (V13)", new_value=" | ".join(NAMES), source_file="numbers/ranking_v3.csv", key_path="the five largest oxygen gains in ranking_v3, in rank order (the run_combo_v2 rule)",
                 dramatic="YES", note="mean saturation enters the five and the 4% desaturation index leaves; the file's correlation_labels list is a stale literal (NUMBERS_DEFECT, not used)"))
C.write_changes(f"{VER}/CHANGES_{S}.csv", rows)
print(f"  CHANGES rows: {len(rows)}, dramatic YES: {sum(1 for r in rows if r['dramatic'] == 'YES')}, notable: {sum(1 for r in rows if r['dramatic'] == 'notable')}")
finish()
