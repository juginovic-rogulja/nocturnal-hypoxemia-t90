#!/usr/bin/env python3
"""V7 (2026-09-07): rebuild numbers/negcontrols_v4_panel.csv from the regenerated bdsp_diseases_v3.csv (per-SD T90,
unpenalized, negative_control rows) and causal_tests_all.csv (BMI-adjusted HRs by condition). The legacy generator
(New_Figures/_NOT_THESE_workfiles/scripts/negcontrols.py) needs raw OMOP first-date columns that the frozen table no longer
carries. Positive control: run on the pre-V7 snapshot, this reproduces the old panel exactly.
usage: negcontrols_panel_v7.py [--v3 PATH --ct PATH --template PATH --out PATH]"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../..")))  # repository root, where paths.py lives
import paths
import argparse, pandas as pd, numpy as np
R = paths.NUMBERS_DIR
ap = argparse.ArgumentParser(); ap.add_argument("--v3", default=f"{R}/bdsp_diseases_v3.csv"); ap.add_argument("--ct", default=f"{R}/causal_tests_all.csv")
ap.add_argument("--template", default=f"{R}/negcontrols_v4_panel.csv"); ap.add_argument("--out", default=f"{R}/negcontrols_v4_panel.csv"); a = ap.parse_args()
v3 = pd.read_csv(a.v3); ct = pd.read_csv(a.ct); tpl = pd.read_csv(a.template)
c = v3[v3.negative_control == True].copy(); assert len(c) == 5, len(c)
p = pd.DataFrame({"key": c.key.values, "condition": c.disease.values, "hr": c.t90_hr.values, "lo": c.t90_lo.values, "hi": c.t90_hi.values, "events": c.events.astype(int).values})
p = p.merge(tpl[["key", "source"]], on="key", how="left")                      # the source label (retained / added) is a fixed design fact
b = ct.set_index("condition").reindex(p.condition.values)[["bmi_adj_hr", "bmi_adj_lo", "bmi_adj_hi"]]
p = pd.concat([p.reset_index(drop=True), b.reset_index(drop=True)], axis=1)[["key", "condition", "source", "hr", "lo", "hi", "events", "bmi_adj_hr", "bmi_adj_lo", "bmi_adj_hi"]]
p = p.set_index("key").loc[tpl.key].reset_index()                                # keep the template's row order
p.to_csv(a.out, index=False); print(p.to_string(index=False)); print(f"written -> {a.out}")
