"""
Step 240 outcomes_added_disclosure (decision 11): old against new for every condition whose
definition changed, as numbers for the "added after the initial analysis" disclosure.

Derived from the tables, never from a typed list:
  added      a condition of DISEASES whose <key>_incident column exists in the v8 table
             (cohort_spec.T90_FINAL) and not in the v7 table (cohort_spec.PREV_DIR)
  changed    a condition present in both whose incident or prevalent count differs
  unchanged  everything else (the 41 untouched conditions must be byte-identical in the build;
             here the counts are printed as a second, independent view of that gate)
For every added or changed condition: old and new incident and prevalent counts, and the per-SD
hazard ratio old (the snapshot of bdsp_diseases_v3.csv taken by step 001 under V8/compare/old_numbers,
else the v7 live file) and new (the regenerated numbers/bdsp_diseases_v3.csv).

Smoke on the v7 tables (PREV_DIR == DATA_DIR): nothing added, nothing changed, every count equal.

Output (cohort_spec.OUT_DIR): outcomes_v8_summary.json.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../..")))  # repository root, where paths.py lives
import paths
import json
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
NUMBERS_LOCAL = paths.NUMBERS_DIR
sys.path.insert(0, paths.ANALYSIS_DIR)
from cohort_spec import T90_FINAL, PREV_DIR, NUMBERS_DIR, OUT_DIR, ROOT, apply_cohort, sidecar  # noqa: E402
from disease_definitions import DISEASES, NEGATIVE_CONTROLS, CVD_COMPONENTS  # noqa: E402

OLD_SNAPSHOT = f"{paths.V8_ROOT}/compare/old_numbers/numbers/bdsp_diseases_v3.csv"
PREV_FINAL = f"{PREV_DIR}/t90_final.parquet"

new = apply_cohort(pd.read_parquet(T90_FINAL))
old_raw = pd.read_parquet(PREV_FINAL)
old = old_raw[(old_raw.fu_valid == 1) & old_raw.spo2_pct_below_90.notna() & (old_raw.oximetry_bad == 0)]
old = old[old.BDSPPatientID.isin(new.BDSPPatientID)]
print(f"new table {len(new):,} nights ({T90_FINAL}), old table {len(old):,} matched nights ({PREV_FINAL})")

keys = [k for k in DISEASES] + ["death"]
label = {k: DISEASES[k][0] for k in DISEASES}
label["death"] = "Death from any cause"


def counts(t, k):
    if f"{k}_incident" not in t.columns:
        return None
    return {"incident": int(t[f"{k}_incident"].fillna(0).sum()),
            "prevalent": int(t[f"{k}_prevalent"].fillna(0).sum()) if f"{k}_prevalent" in t.columns else None}


hr_new = pd.read_csv(f"{NUMBERS_DIR}/bdsp_diseases_v3.csv", comment="#").set_index("key") \
    if os.path.exists(f"{NUMBERS_DIR}/bdsp_diseases_v3.csv") else None
old_hr_src = OLD_SNAPSHOT if os.path.exists(OLD_SNAPSHOT) else f"{NUMBERS_DIR}/bdsp_diseases_v3.csv"
hr_old = pd.read_csv(old_hr_src, comment="#").set_index("key") if os.path.exists(old_hr_src) else None


def hr_of(tab, k):
    if tab is None or k not in tab.index:
        return None
    r = tab.loc[k]
    return {"hr": float(r.t90_hr), "lo": float(r.t90_lo), "hi": float(r.t90_hi), "p": float(r.t90_p), "events": int(r.events)}


added, changed, unchanged, missing_new = [], [], [], []
for k in keys:
    cn, co = counts(new, k), counts(old, k)
    rec = {"key": k, "label": label[k], "negative_control": k in NEGATIVE_CONTROLS,
           "old": co, "new": cn, "hr_old": hr_of(hr_old, k), "hr_new": hr_of(hr_new, k),
           "cvd_component": k in CVD_COMPONENTS}
    if cn is None:
        missing_new.append(rec)
    elif co is None:
        added.append(rec)
    elif cn != co:
        rec["incident_change"] = cn["incident"] - co["incident"]
        changed.append(rec)
    else:
        unchanged.append(rec)

out = {"new_table": T90_FINAL, "old_table": PREV_FINAL, "old_hr_source": old_hr_src,
       "n_nights_new": int(len(new)), "n_nights_old_matched": int(len(old)),
       "n_conditions": len(keys), "n_added": len(added), "n_changed": len(changed), "n_unchanged": len(unchanged),
       "n_missing_in_new": len(missing_new),
       "added_after_the_initial_analysis": [r["label"] for r in added],
       "definitions_fixed": [r["label"] for r in changed],
       "added": added, "changed": changed, "missing_in_new": missing_new,
       "unchanged_labels": [r["label"] for r in unchanged]}
p = f"{OUT_DIR}/outcomes_v8_summary.json"
json.dump(out, open(p, "w"), indent=1)
sidecar(p, __file__, extra_inputs=[PREV_FINAL] + [x for x in (OLD_SNAPSHOT, f"{NUMBERS_DIR}/bdsp_diseases_v3.csv") if os.path.exists(x)],
        note=f"outcomes summary: {len(added)} added, {len(changed)} changed, {len(unchanged)} unchanged")
print(f"added {len(added)}: {[r['label'] for r in added]}")
print(f"changed {len(changed)}: " + ", ".join(f"{r['label']} {r['old']['incident']}->{r['new']['incident']}" for r in changed))
print(f"unchanged {len(unchanged)}, missing in new {len(missing_new)}; written {p}")
