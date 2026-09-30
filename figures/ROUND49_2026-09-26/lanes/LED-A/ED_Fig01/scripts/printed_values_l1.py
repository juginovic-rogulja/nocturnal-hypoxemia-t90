#!/usr/bin/env python3
"""verify/<Sheet>_printed_values.csv for lane V14_L1_RANK (rules section 3.7): every text span of the delivered sheet that carries a
digit (plus every drawn string), matched to the builder's drawn records and re-read from the source file where the key is resolvable.
The round-30 style records (text, x, baseline, ha, source, key, value, kind) are converted to v14lib DRAWN_RECORD rows with a resolvable
CSV/JSON key and a printing rule where one exists; the rest carry the recorded value under the rule (match = yes_recorded).
usage: printed_values_l1.py <Sheet>"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import os, sys, json, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND37_2026-09-14/figures/_common")
import v14lib as L

S = sys.argv[1]; SD = f"{C.LANE}/{S}"; VER = f"{SD}/verify"
NEW_PDF = C.hydrated(f"{SD}/{S}.pdf"); BASE_PDF = C.hydrated(f"{C.BASE}/{S}.pdf")
SRC = {"ranking_v3.csv": f"{C.NUM}/ranking_v3.csv", "ranking_v3.csv x NICE": f"{C.NUM}/ranking_v3.csv", "ranking_v3.csv x measure_families.csv": f"{C.NUM}/ranking_v3.csv",
       "measure_families.csv": f"{C.NUM}/measure_families.csv", "results_v2.json": f"{C.NUM}/results_v2.json", "figure1_sleep_profile_v1.json": f"{C.NUM}/figure1_sleep_profile_v1.json",
       "combination_v2.json": f"{C.NUM}/combination_v2.json", "dC_per_disease.csv": f"{C.SV}/dC_per_disease/dC_per_disease.csv",
       "L2_fig1d_values.json": f"{paths.FIGURE_ROOT}/ROUND15_2026-09-04/L2_fig1d/L2_fig1d_values.json", "axis": "static:axis", "band": "static:band"}
# extra printing rules for this lane's quantities
L.RULES["gain_x1000_1dp"] = lambda v: f"{float(v) * 1000.0:.1f}"                 # the ED 1 gain column and the Fig 1 bars (one decimal, plain float as the builder prints)
L.RULES["gain_4dp"] = lambda v: f"{float(v):.4f}"
L.RULES["corr_2dp_minus"] = lambda v: f"{float(v):.2f}".replace("-", "−")
L.RULES["rank_callout"] = lambda v: str(int(v))
L.RULES["dc_ci_3dp"] = lambda v: f"{float(v['dC']):+.3f} ({float(v['ci_lo']):.3f} to {float(v['ci_hi']):.3f})".replace("-", "−")
L.RULES["logrank"] = lambda v: "Log-rank " + ("P < 0.001" if float(v["p"]) < 0.001 else f"P = {float(v['p']):.3f}")
L.RULES["g_int"] = lambda v: f"{float(v):g}"


def convert(panel, v):
    """One round-30 record -> one v14lib record with the best resolvable key."""
    src = v.get("source", ""); key = v.get("key", ""); val = v.get("value"); kind = v.get("kind"); text = str(v["text"])
    rec = dict(text=text, x=v["x"], baseline=v["baseline"], ha=v.get("ha", "left"), panel=panel, source_file=SRC.get(src, src), source_key=key, source_value=val, rule="text")
    if src in ("axis", "band") or kind == "tick":
        # tick labels are centred on the tick (x ticks) or right-aligned to the axis (y ticks): the drawn x is the anchor, not the left edge
        if kind == "tick": rec["ha"] = "center" if "x tick" in key else ("right" if "y tick" in key else rec["ha"])
        rec.update(source_file="static:V13", source_key=f"{kind}: {key}", rule="text", source_value=text); return rec
    m = re.search(r"where feature == (\S+)", key)
    if src.startswith("ranking_v3.csv") and m and ("dC x 1000" in key or "gain" in key.lower()) and kind == "printed":
        rec.update(source_file=SRC["ranking_v3.csv"], source_key=f"feature={m.group(1)}:dC", rule="gain_x1000_1dp"); return rec
    if src == "ranking_v3.csv" and m and key.startswith("rank where"):
        rec.update(source_key=f"feature={m.group(1)}:rank", rule="rank_callout" if not text.endswith(m.group(1)) else "text")
        if kind == "printed" and ", rank " in text:   # the Fig 1c callouts: "<name>, rank N": the re-read rank must equal the number at the end
            n = text.rsplit("rank ", 1)[-1]; rec.update(rule="text", source_value=text, note=f"rank {n} re-read separately"); rec["_rank_key"] = (f"feature={m.group(1)}:rank", n)
            rec["source_key"] = f"rank where feature == {m.group(1)}, printed as '<name>, rank N' (the rank is re-read below, the name is the label map)"
        return rec
    if src == "measure_families.csv" and m:
        rec.update(source_key=f"feature={m.group(1)}:family", rule="text"); return rec
    if src == "ranking_v3.csv x NICE" and kind == "label":
        rec.update(source_file="static:name map", source_key=key, rule="text", source_value=text); return rec
    if src == "figure1_sleep_profile_v1.json" and kind == "label" and key.startswith("measures."):
        rec.update(source_key=key.replace(".", "/"), rule="text"); return rec
    if src == "combination_v2.json" and kind == "printed" and key.startswith(("single.", "combinations.")):
        rec.update(source_key=key.replace(".", "/"), rule="gain_4dp"); return rec
    if src == "combination_v2.json" and kind == "label":
        rec.update(source_file="static:row logic", rule="text", source_value=text); return rec
    if src == "combination_v2.json" and kind == "printed" and key.startswith("correlation_matrix["):
        a, b = [int(t) for t in re.findall(r"\[(\d+)\]", key)[:2]]
        full = v.get("value_full_precision")
        if v.get("rounding_tie"): rec.update(source_file="static:in-lane full precision", source_key=f"corr_check.json recomputed[{a}][{b}] (rounding tie of the file's 3 dp)", rule="corr_2dp_minus", source_value=full)
        else: rec.update(source_key=f"correlation_matrix/{a}/{b}", rule="corr_2dp_minus")
        return rec
    if src == "ranking_v3.csv" and kind == "label" and ("ranking_v3 order" in key):
        rec.update(source_file="static:label map", rule="text", source_value=text); return rec
    if src == "dC_per_disease.csv" and kind == "printed" and isinstance(val, dict):
        dz = re.search(r"disease==(.+?): ", key); dz = dz.group(1) if dz else None
        rec.update(source_key=f"row_type=disease&arm=standard&disease={dz}:dC,ci_lo,ci_hi" if dz else key, rule="dc_ci_3dp", source_value=val)
        rec["_dc_row"] = dz; return rec
    if src == "dC_per_disease.csv" and kind == "label":
        rec.update(source_file="static:row set", rule="text", source_value=text); return rec
    if src == "L2_fig1d_values.json":
        if kind == "label": rec.update(source_file="static:tile title", rule="text", source_value=text); return rec
        if isinstance(val, dict) and "p" in val: rec.update(rule="logrank", source_value=val); return rec
        if kind == "tick": rec.update(rule="g_int", source_value=val); return rec
    rec.update(rule="text", source_value=text); return rec


records = []; tie_note = []
for fn in sorted(os.listdir(VER)):
    if not fn.endswith("_drawn.json"): continue
    d = json.load(open(f"{VER}/{fn}")); panel = d.get("panel", fn.split("_")[-2])
    for v in d["values"]:
        if v.get("kind") not in ("printed", "tick", "label"): continue
        records.append(convert(panel, v))
# every other string the builder placed (static captions, tick labels, letters, second lines of wrapped names) is accounted for from
# work/placed_texts.json: it was placed by the builder from the V13 geometry (moved with a grown or shrunk panel where the sheet says so)
placed = json.load(open(f"{SD}/work/placed_texts.json"))
def already(text, x, base, tol=0.3):
    return any(str(r["text"]) == text and abs(float(r["x"]) - x) <= tol and abs(float(r["baseline"]) - base) <= tol for r in records)
for t in placed:
    if already(str(t["text"]), float(t["x"]), float(t["baseline"])): continue
    rec_ = t.get("record") or {}
    records.append(dict(text=str(t["text"]), x=t["x"], baseline=t["baseline"], ha=t.get("ha", "left"), panel=rec_.get("panel", ""), source_file="static:V13 (placed by the builder at the V13 geometry)",
                        source_key="static string" + ("" if rec_.get("static", True) else " (declared change)"), source_value=str(t["text"]), rule="text"))
# Main_Fig1: the three band strings re-stamped with the v8.1 counts (compose_record.json), the title and the letters (static, at the V13 origins)
if S == "Main_Fig1":
    cr = json.load(open(f"{SD}/work/compose_record.json"))
    for band, y_off in (("band_a", cr["band_a"]["slot"][1]), ("band_f", cr["band_f"]["slot"][1])):
        for e in cr[band]["edits"]:
            n = e["new"]; txt = n["text"].replace("\xa0", " ")
            val = "141" if "141" in txt else "54"
            records.append(dict(text=txt, x=n["origin"][0], baseline=n["origin"][1] + y_off, ha="left", panel=band[-1], source_file=f"{C.NUM}/ranking_v3.csv" if val == "141" else f"{C.NUM}/bdsp_diseases_v3.csv",
                                source_key="row count = N_MEASURES (141)" if val == "141" else "53 non-control outcomes = 54 diseases and death (the v8 count, decision 11)", source_value=txt, rule="text",
                                note=f"band string re-stamped at the V13 origin: {e['old']['text']} -> {txt}"))
# resolve the special keys ourselves (callout ranks, dC rows) so the CSV carries a re-read where the generic resolver cannot
import pandas as pd
for r in records:
    if "_rank_key" in r:
        k, n = r.pop("_rank_key")
        got = L.resolve_source(SRC["ranking_v3.csv"], k); r["note"] = (r.get("note", "") + f"; rank re-read from ranking_v3.csv = {int(got)}").strip("; ")
        if str(int(got)) != n: r["rule"] = "MISMATCH"; r["source_value"] = f"re-read rank {int(got)} vs printed {n}"
    if "_dc_row" in r:
        dz = r.pop("_dc_row")
        if dz:
            df = pd.read_csv(SRC["dC_per_disease.csv"]); row = df[(df.row_type == "disease") & (df.arm == "standard") & (df.disease == dz)]
            assert len(row) == 1, dz
            r["source_value"] = dict(dC=float(row.dC.iloc[0]), ci_lo=float(row.ci_lo.iloc[0]), ci_hi=float(row.ci_hi.iloc[0])); r["source_key"] = f"row_type==disease & arm==standard & disease=={dz}: dC, ci_lo, ci_hi"
for r in records:
    if r["rule"] == "MISMATCH": r["rule"] = "text"; r["source_value"] = "MISMATCH"
out = f"{VER}/{S}_printed_values.csv"
n, n_no = L.printed_values_csv(NEW_PDF, records, out, static_from=BASE_PDF)
print(f"{S}: {n} rows, {n_no} NO -> {out}")
if n_no:
    import csv
    for row in csv.DictReader(open(out)):
        if row["match"] == "NO": print("   NO:", row["string"], "|", row["note"][:120])
