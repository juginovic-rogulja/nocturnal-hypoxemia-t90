"""
eFigure 9 refitted onto the revised negative-control panel.

Two things had to change and only one of them was asked for.

  1. THE PANEL, settled 2026-08-07. Five controls: back pain, cataract, glaucoma, contact
     dermatitis, hemorrhoids. Fracture AND osteoarthritis are both out, for the same reason,
     which is that a plausible causal path exists from the exposure to each of them. Fracture:
     hypoxemia degrades alertness, alertness loss causes falls, falls cause fracture.
     Osteoarthritis: published work links it to hypoxia and to deoxygenation. That is the whole
     argument; the body-mass-index behaviour of either condition is not the reason and must not
     be given as the reason. An earlier version of this file reinstated osteoarthritis on the
     ground that it set the 1.133 floor. That block is superseded. The floor is 1.063, set by
     glaucoma. eFigure 9 carried osteoarthritis and back pain, so the flat row Alen liked is now
     five conditions wide instead of two. All five are read from the frozen outcome parquet,
     which carries the two replacements since numbers/add_newcontrols_v5.py.

  2. THE AXIS LABEL, which was wrong. The published eFigure 9 is labelled "per 1 SD of time
     with oxygen <90%". The exposure actually fitted in numbers/run_extras2_v2.py is the binary
     contrast lowox = (T90 > 10), severe apnea with low oxygen against severe apnea without.
     Those are different quantities and the figure was reporting one under the name of the
     other. The refit here keeps the exposure that was fitted and corrects the label.

Everything else is held: apnea-hypopnea index 30 or more, unpenalized Cox stratified by site,
age spline with one column dropped, sex, oximetry_bad == 0.

Writes numbers/severe_apnea_negctrl_v1.json.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../..")))  # repository root, where paths.py lives
import paths
import json
import sys
import warnings

import duckdb
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from patsy import dmatrix

warnings.filterwarnings("ignore")
ROOT = paths.FIGURE_ROOT
CACHE = paths.OMOP_CACHE_DIR
sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
from cohort_spec import apply_cohort, COHORT_N
from disease_definitions import DISEASES, NEGATIVE_CONTROLS

b = apply_cohort(pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet"))
assert len(b) == COHORT_N
b["t90"] = b.spo2_pct_below_90
b["male"] = (b.sex.astype(str).str.upper().str[0] == "M").astype(int)
sp = dmatrix("cr(a, df=4) - 1", {"a": b.AgeAtVisit.values}, return_type="dataframe").iloc[:, 1:]
for i in range(sp.shape[1]):
    b[f"age_s{i}"] = sp.iloc[:, i].values
ADJ = [f"age_s{i}" for i in range(sp.shape[1])] + ["male"]
assert len(ADJ) == 4

# ---- the two replacement controls -----------------------------------------------------------
# They used to be coded here on the fly from the OMOP cache, because they were not in the
# frozen outcome parquet. They are now, added by numbers/add_newcontrols_v5.py from the same
# ICD prefixes, so they are read like every other condition and cannot drift from the panel.
NEW = {k: (DISEASES[k][0], sorted(set(DISEASES[k][2]) | set(DISEASES[k][3])))
       for k in ("dermatitis_contact", "haemorrhoids")}
for k in NEW:
    assert f"{k}_incident" in b.columns, f"{k} missing: run numbers/add_newcontrols_v5.py"

# ---- validation: the two new controls must reproduce the per-SD fit of bdsp_diseases_v3.csv (v8) ----
panel = pd.read_csv(f"{paths.NUMBERS_DIR}/bdsp_diseases_v3.csv").rename(columns={"t90_hr": "hr"})   # v8: the regenerated per-SD table (step 111, sidecar; same unpenalized site-stratified model). negcontrols_v4_panel.csv is a v7-era static file with no v8 writer.
from scipy import stats as _st


def _rint(x):
    r = _st.rankdata(x)
    return _st.norm.ppf((r - 0.375) / (len(r) + 0.25))


b["z"] = b.groupby("site_id").spo2_pct_below_90.transform(_rint)
print("validation, the two new controls refitted per 1 SD on the whole cohort")
for k, (lab, _p) in NEW.items():
    f = b[(b[f"{k}_prevalent"] == 0) & b[f"{k}_years"].notna() & (b[f"{k}_years"] > 0)]
    d = pd.DataFrame({"T": f[f"{k}_years"], "E": f[f"{k}_incident"].astype(int),
                      "site": f.site_id, "z": f.z, **{c: f[c] for c in ADJ}}).dropna()
    c = CoxPHFitter(penalizer=0.0).fit(d, "T", "E", strata=["site"])
    hr = float(np.exp(c.params_["z"]))
    ref = panel[panel.key == k].iloc[0]
    print(f"   {lab:<20}refit {hr:.3f} ev {int(d.E.sum()):,}   frozen {ref.hr:.3f} "
          f"ev {int(ref.events):,}   "
          f"{'match' if abs(hr - ref.hr) < 0.002 and int(d.E.sum()) == int(ref.events) else 'MISMATCH'}")
    assert abs(hr - ref.hr) < 0.002 and int(d.E.sum()) == int(ref.events), lab

# ---- the severe-apnea subgroup --------------------------------------------------------------
sev = b[b.AHI.notna() & (b.AHI >= 30)].copy()
sev["lowox"] = (sev.t90 > 10).astype(int)
# TWO DIFFERENT COUNTS, and the figure used to print the wrong one in its title.
#   n_ref  is the REFERENCE GROUP of the contrast actually drawn, T90 at or below 10% of the
#          night. This is what the title must name, because every hazard ratio on the sheet is
#          measured against it.
#   n_norm is the much stricter "essentially normal oxygen" group, T90 at or below 1%. It is a
#          descriptive statistic about how many people with severe apnea keep normal oxygen. It
#          is NOT the reference group and the title said it was.
n_ref = int((sev.lowox == 0).sum())
n_low = int((sev.lowox == 1).sum())
n_norm = int((sev.t90 <= 1).sum())
print(f"\nsevere sleep apnea, apnea-hypopnea index 30 or more: n = {len(sev):,}, "
      f"reference group {n_ref:,} ({100 * n_ref / len(sev):.1f}%) at or below 10%, "
      f"of whom {n_norm:,} ({100 * n_norm / len(sev):.1f}%) were at or below 1%")

TEST = [("hf", "Heart failure"), ("cvd", "Cardiovascular composite"),
        ("resp_failure", "Respiratory failure"), ("diabetes", "Type 2 diabetes"),
        ("death", "Death from any cause"), ("copd2", "COPD"),
        ("aki", "Acute kidney injury")]
# imported, never listed: the panel cannot drift from the source of truth
CTRL = [(k, DISEASES[k][0]) for k in NEGATIVE_CONTROLS]

R = {"_built_by": "New_Figures/_working/scripts/fits_severe_apnea_negctrl.py",
     "exposure": "binary: time below 90% saturation above 10% of the night, against 10% or less",
     "exposure_is_not_per_sd": ("the published eFigure 9 axis label said per 1 SD; the fitted "
                                "exposure in run_extras2_v2.py is and always was this binary "
                                "contrast"),
     "subgroup": "apnea-hypopnea index 30 or more",
     "n": int(len(sev)),
     "n_reference": n_ref, "pct_reference": round(100 * n_ref / len(sev), 1),
     "n_low_oxygen": n_low, "pct_low_oxygen": round(100 * n_low / len(sev), 1),
     "n_normal_oxygen": n_norm,
     "pct_normal_oxygen": round(100 * n_norm / len(sev), 1),
     "_reference_note": ("n_reference is the comparison group of every hazard ratio on the "
                         "sheet, T90 at or below 10%. n_normal_oxygen is the stricter T90 at "
                         "or below 1% group and is descriptive only. The published eFigure 9 "
                         "title named n_normal_oxygen as if it were the reference group."),
     "negative_controls_revised": [lab for _k, lab in CTRL],
     "negative_controls_dropped": ["Fracture", "Osteoarthritis"],
     "outcomes": {}}

print(f"\n{'outcome':<28}{'ev':>6}   hazard ratio (95% CI)")
for k, lab in TEST + CTRL:
    if f"{k}_incident" not in sev.columns:
        print(f"   {lab:<25} not available")
        continue
    g = sev[(sev[f"{k}_prevalent"] == 0) & sev[f"{k}_years"].notna() & (sev[f"{k}_years"] > 0)]
    f = pd.DataFrame({"T": g[f"{k}_years"], "E": g[f"{k}_incident"].astype(int),
                      "x": g.lowox, "site": g.site_id,
                      **{c: g[c] for c in ADJ}}).dropna()
    if f.E.sum() < 60:
        print(f"   {lab:<25} {int(f.E.sum())} events, below the floor")
        continue
    c = CoxPHFitter(penalizer=0.0).fit(f, "T", "E", strata=["site"])
    lo, hi = np.exp(c.confidence_intervals_.loc["x"])
    hr, p = float(np.exp(c.params_["x"])), float(c.summary.loc["x", "p"])
    R["outcomes"][lab] = {"events": int(f.E.sum()), "n": int(len(f)),
                          "hr": hr, "lo": float(lo),  # v8.3 full precision (2026-09-16)
                          "hi": float(hi), "p": p,
                          "negative_control": (k, lab) in CTRL}
    print(f"{lab:<28}{int(f.E.sum()):>6}   {hr:.2f} (95% CI, {lo:.2f}-{hi:.2f})"
          f"{'  *' if p < .05 else ''}")

# ---- reproduce the frozen values for the conditions that were already there -----------------
old = json.load(open(f"{paths.NUMBERS_DIR}/extras2_v2.json"))["severe_apnea"]["outcomes"]
bad = [(lab, old[lab]["hr"], R["outcomes"][lab]["hr"]) for lab in R["outcomes"]
       if lab in old and abs(old[lab]["hr"] - R["outcomes"][lab]["hr"]) > 0.002]
print(f"\nvalidation against extras2_v2.json on the shared conditions: {len(bad)} mismatches")
assert not bad, bad

ctrl_hr = [v["hr"] for v in R["outcomes"].values() if v["negative_control"]]
R["control_range"] = [min(ctrl_hr), max(ctrl_hr)]
R["n_controls_significant"] = sum(1 for v in R["outcomes"].values()
                                  if v["negative_control"] and v["p"] < 0.05)
print(f"revised controls span {min(ctrl_hr):.2f} to {max(ctrl_hr):.2f}, "
      f"{R['n_controls_significant']} of {len(ctrl_hr)} significant")
json.dump(R, open(f"{paths.NUMBERS_DIR}/severe_apnea_negctrl_v1.json", "w"), indent=1)
print("wrote numbers/severe_apnea_negctrl_v1.json")
