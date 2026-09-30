"""
Step 233 fig2c_top10_rule (decision 13): the Figure 2c row set as a number file.

Rule (stated in the legend): the ten largest top-band hazard ratios (T90 above 10 percent against
1 percent or less) among the non-control outcomes, excluding obesity hypoventilation syndrome
(diagnosed partly on nocturnal oxygen) and obesity, read from results_v2.json ["graded"]. Every
row that is drawn today (ROUND14 R14_FIG2B_DRAWN.json row_order) and is not among the ten is
listed as leaving; those rows stay in Supp Table 5, which carries every outcome.

Output (cohort_spec.OUT_DIR): fig2c_selection_v8.json.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../..")))  # repository root, where paths.py lives
import paths
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from cohort_spec import NUMBERS_DIR, ROOT, sidecar, OUT_DIR  # noqa: E402
from disease_definitions import DISEASES, NEGATIVE_CONTROLS  # noqa: E402

TOP_BAND = ">10%"
EXCLUDE_KEYS = ["obesity_hypovent", "obesity"]
N_ROWS = 10
RESULTS = f"{NUMBERS_DIR}/results_v2.json"
DRAWN = f"{paths.FIGURE_ROOT}/ROUND14_2026-09-04/L3_fig2b/R14_FIG2B_DRAWN.json"

G = json.load(open(RESULTS))["graded"]
outcomes = G["outcomes"]
excl_labels = {DISEASES[k][0] for k in EXCLUDE_KEYS}
control_labels = {DISEASES[k][0] for k in NEGATIVE_CONTROLS}
missing = [k for k in EXCLUDE_KEYS if k not in DISEASES]
assert not missing, f"exclusion keys not in DISEASES: {missing}"

cand = []
for lab, v in outcomes.items():
    if v.get("negative_control") or lab in control_labels or lab in excl_labels:
        continue
    top = v.get(TOP_BAND, {})
    if "hr" not in top:
        continue
    cand.append({"outcome": lab, "top_band_hr": float(top["hr"]), "lo": float(top.get("lo", float("nan"))),
                 "hi": float(top.get("hi", float("nan"))), "p": top.get("p"), "events": v.get("events"),
                 "trend_p": v.get("trend_p")})
cand.sort(key=lambda r: -r["top_band_hr"])
chosen = cand[:N_ROWS]
for i, r in enumerate(chosen, 1):
    r["rank"] = i
chosen_labels = [r["outcome"] for r in chosen]

drawn_now = json.load(open(DRAWN))["row_order"] if os.path.exists(DRAWN) else []
leaving = [lab for lab in drawn_now if lab not in chosen_labels]
entering = [lab for lab in chosen_labels if lab not in drawn_now]
out = {"rule": (f"the {N_ROWS} largest hazard ratios in the top T90 band ({TOP_BAND} against 0-1%) among non-control outcomes, "
                f"excluding {sorted(excl_labels)}; the rows that leave stay in Supp Table 5"),
       "source": RESULTS, "band": TOP_BAND, "band_n": G.get("band_n"),
       "excluded": sorted(excl_labels), "n_candidates": len(cand),
       "selected": chosen, "selected_labels": chosen_labels,
       "drawn_today_source": DRAWN if drawn_now else None, "drawn_today": drawn_now,
       "leaving": leaving, "entering": entering,
       "next_five_below_the_cut": [r["outcome"] for r in cand[N_ROWS:N_ROWS + 5]]}
p = f"{OUT_DIR}/fig2c_selection_v8.json"
json.dump(out, open(p, "w"), indent=1)
sidecar(p, __file__, extra_inputs=[RESULTS] + ([DRAWN] if drawn_now else []),
        note=f"fig2c rule: {len(chosen)} rows, {len(leaving)} leave the current grid, {len(entering)} enter")
print(f"Figure 2c rows ({len(chosen)}): " + ", ".join(f"{r['outcome']} {r['top_band_hr']:.2f}" for r in chosen))
print(f"leave the current grid ({len(leaving)}): {leaving}")
print(f"enter ({len(entering)}): {entering}")
print("written:", p)
