#!/usr/bin/env python3
"""ED_Fig02 relaunch 3, step 1 (data): the generator's row rule (FINAL_FIGURES_2026-08-14/_scripts/eFigure07_polish.py lines 182-194)
run on the OLD snapshot and on the NEW v7 files. Writes work/data_old.json and work/data_new.json with every drawn value and every
printed string, the sidecar facts and the sha256 of the files read. The OLD strings are the positive control against the base sheet.
Sources: numbers/lag_ladder.csv (+ .json twin, step 30), numbers/causal_tests_all.csv (step E2), numbers/sleep_unadjusted_v1.json (step 25),
numbers/disease_definitions.py for the control order."""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ed2_common import derive, NUM, OLD, WORK
old = derive(OLD, "old"); new = derive(NUM, "new")
json.dump(old, open(os.path.join(WORK, "data_old.json"), "w"), indent=1); json.dump(new, open(os.path.join(WORK, "data_new.json"), "w"), indent=1)

print("OLD rows:", old["conds"]); print("NEW rows:", new["conds"])
print("OLD generator xlim", old["generator_xlim"], "NEW", new["generator_xlim"], "drawn range new", new["drawn_min"], new["drawn_max"])
print("ranked lag0 hr NEW (first 18):", json.dumps(new["lag0_hr_ranked"]))
print("ranked lag0 hr OLD (first 18):", json.dumps(old["lag0_hr_ranked"]))
print("NE count old/new", old["n_ne"], new["n_ne"], " grey-5 count old/new", old["n_grey5"], new["n_grey5"])
print("sidecars:", json.dumps(new["sidecars"], indent=0))
# positive control: the OLD printed strings must all be on the base sheet's text layer, and the OLD row order must match it
T = json.load(open(os.path.join(WORK, "base_text.json")))
toks = [s["text"] for s in T["spans"]]
labels_base = [s["text"] for s in sorted([s for s in T["spans"] if abs(s["origin"][0] - 12.96) < 0.05 and s["size"] == 10.0 and s["text"] != "Negative controls"], key=lambda s: s["origin"][1])]
print("base row labels:", labels_base)
print("OLD rows == base labels:", [r["condition"] for r in old["rows"]] == labels_base)
print("NEW rows == base labels:", [r["condition"] for r in new["rows"]] == labels_base)
col1_base = [s["text"] for s in sorted([s for s in T["spans"] if abs(s["origin"][0] - 248.69) < 0.05 and s["size"] == 9.5 and s["text"] != "1-year"], key=lambda s: s["origin"][1])]
col5_base = [s["text"] for s in sorted([s for s in T["spans"] if abs(s["origin"][0] - 286.13) < 0.05 and s["size"] == 9.5 and s["text"] != "5-year"], key=lambda s: s["origin"][1])]
print("OLD 1-year strings == base column:", [r["a"]["y1_text"] for r in old["rows"]] == col1_base, [r["a"]["y1_text"] for r in old["rows"]], col1_base)
print("OLD 5-year strings == base column:", [r["a"]["y5_text"] for r in old["rows"]] == col5_base, [r["a"]["y5_text"] for r in old["rows"]], col5_base)
print("NEW 1-year:", [r["a"]["y1_text"] for r in new["rows"]]); print("NEW 5-year:", [r["a"]["y5_text"] for r in new["rows"]])
print("NEW grey 5-year rows:", [r["condition"] for r in new["rows"] if r["a"]["y5_grey"]]); print("OLD grey 5-year rows:", [r["condition"] for r in old["rows"] if r["a"]["y5_grey"]])
