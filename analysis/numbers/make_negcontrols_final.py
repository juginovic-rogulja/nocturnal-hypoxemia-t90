"""
Build numbers/negcontrols_final.json, the one file every downstream consumer should read for
the negative-control panel and the confounding floor.

Nothing here is typed in. The panel comes from disease_definitions.NEGATIVE_CONTROLS, the
primary per-SD estimates from bdsp_diseases_v3.csv, the landmark and body-mass-index tests from
causal_tests_all.csv, and the graded band estimates from results_v2.json. The floor is the
largest control hazard ratio on the primary per-SD scale, computed, not asserted.

The reason for a consolidated file rather than a constant: on 2026-08-07 a census found 87
scripts carrying a hardcoded control list or floor, four mutually incompatible panels across the
delivered figures, and three printed floors, none of them right. A count check does not catch
that, because every wrong panel also had five members. An identity check against this file does.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../..")))  # repository root, where paths.py lives
import paths
import hashlib
import json
import os
import sys
import time
from datetime import date

import pandas as pd

ROOT = paths.T90_ROOT
sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
from disease_definitions import (DISEASES, NEGATIVE_CONTROLS,          # noqa: E402
                                 NEGATIVE_CONTROL_LABELS, CONFOUNDING_FLOOR_PER_SD)
from cohort_spec import COHORT_N                                       # noqa: E402  v8 F2 fix

PATH = f"{paths.NUMBERS_DIR}/negcontrols_final.json"
INPUTS = {"bdsp_diseases_v3.csv": f"{paths.NUMBERS_DIR}/bdsp_diseases_v3.csv",
          "causal_tests_all.csv": f"{paths.NUMBERS_DIR}/causal_tests_all.csv",
          "results_v2.json": f"{paths.NUMBERS_DIR}/results_v2.json"}


def _describe(p):
    st = os.stat(p)
    assert st.st_blocks > 0, f"{p} is evicted (0 blocks)"
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return {"path": p, "bytes": st.st_size, "sha256": h.hexdigest(),
            "mtime_local": time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(st.st_mtime))}


_INPUTS = {k: _describe(v) for k, v in INPUTS.items()}

v3 = pd.read_csv(f"{paths.NUMBERS_DIR}/bdsp_diseases_v3.csv").set_index("key")
ct = pd.read_csv(f"{paths.NUMBERS_DIR}/causal_tests_all.csv").set_index("condition")
graded = json.load(open(f"{paths.NUMBERS_DIR}/results_v2.json"))["graded"]["outcomes"]

controls = {}
for k in NEGATIVE_CONTROLS:
    label = DISEASES[k][0]
    r = v3.loc[k]
    c = ct.loc[label] if label in ct.index else None
    rec = {
        "key": k,
        "label": label,
        "icd9_prefixes": DISEASES[k][2],
        "icd10_prefixes": DISEASES[k][3],
        "events": int(r.events),
        "at_risk": int(r.n),
        "per_sd": {"hr": float(r.t90_hr), "lo": float(r.t90_lo), "hi": float(r.t90_hi),
                   "p": float(r.t90_p)},
    }
    if c is not None:
        rec["landmark_2y"] = {"hr": float(c.lag_hr), "lo": float(c.lag_lo),
                              "hi": float(c.lag_hi), "events": int(c.lag_events)}
        rec["bmi_adjusted"] = {"hr": float(c.bmi_adj_hr), "lo": float(c.bmi_adj_lo),
                               "hi": float(c.bmi_adj_hi), "events": int(c.bmi_events),
                               "same_patients_unadjusted_hr": float(c.bmi_sub_hr)}
    if label in graded:
        rec["graded"] = graded[label]
    controls[k] = rec

hrs = {k: controls[k]["per_sd"]["hr"] for k in controls}
floor = max(hrs.values())
driver = max(hrs, key=hrs.get)
# v8 F2 fix 2026-09-12 (integrity gate 6): the floor and its driver are RESULTS and are no longer asserted against
# typed values (the v7 writer pinned floor == disease_definitions.CONFOUNDING_FLOOR_PER_SD and driver == "glaucoma").
# What is checked is definitional: every panel member is flagged negative_control in the primary table, its estimate
# is on the primary per-SD scale (exposure spo2_pct_below_90), and the driver is a panel member. Old versus new
# (rule 11) is printed and stored: the previous json and the module constant are compared, never asserted.
assert driver in NEGATIVE_CONTROLS, driver
assert all(bool(v3.loc[k, "negative_control"]) for k in NEGATIVE_CONTROLS), \
    v3.loc[NEGATIVE_CONTROLS, "negative_control"].to_dict()
assert all(v3.loc[k, "exposure"] == "spo2_pct_below_90" for k in NEGATIVE_CONTROLS), \
    v3.loc[NEGATIVE_CONTROLS, "exposure"].to_dict()
_prev = json.load(open(PATH)) if os.path.exists(PATH) else None
_prev_floor = float(_prev["confounding_floor_per_sd"]) if _prev else None
_prev_driver = _prev.get("confounding_floor_driver") if _prev else None
OLD_VS_NEW = {
    "previous_file": ({"built_on": _prev.get("_built_on"), "floor": _prev_floor, "driver": _prev_driver,
                       "controls_per_sd_hr": {k: _prev["controls"][k]["per_sd"]["hr"] for k in _prev["controls"]}}
                      if _prev else None),
    "module_constant_disease_definitions": CONFOUNDING_FLOOR_PER_SD,
    "new": {"floor": round(floor, 3), "driver": DISEASES[driver][0],
            "controls_per_sd_hr": {k: round(v, 3) for k, v in hrs.items()}},
    "floor_moved_up": bool(_prev and round(floor, 3) > _prev_floor),
    "dramatic": bool(_prev and (_prev_driver != DISEASES[driver][0]
                                or abs(round(floor, 3) - _prev_floor) / _prev_floor > 0.2)),
}

removed = {}
for k in ["fracture", "osteoarthritis"]:
    r = v3.loc[k]
    removed[k] = {
        "label": DISEASES[k][0],
        "status": "ordinary outcome, not a negative control",
        "per_sd": {"hr": float(r.t90_hr), "lo": float(r.t90_lo), "hi": float(r.t90_hi)},
        "events": int(r.events),
    }
removed["fracture"]["removal_reason"] = (
    "A plausible causal path runs from the exposure. Nocturnal hypoxemia degrades daytime "
    "alertness and cognition, alertness loss causes falls, and falls are the dominant mechanism "
    "of fracture in this age range. The body-mass-index behaviour is corroborating rather than "
    "disqualifying: fracture is the only outcome in the paper whose hazard ratio RISES on "
    "body-mass-index adjustment, from 1.144 to 1.223 in the same patients, because high body "
    "mass protects against fracture while also driving low nocturnal oxygen. That rise is to be "
    "explained, not used as a floor.")
removed["osteoarthritis"]["removal_reason"] = (
    "Alen's domain call: published literature links osteoarthritis to hypoxia and to "
    "deoxygenation, so a plausible causal path from the exposure exists. This is the whole "
    "reason. The earlier argument from its body-mass-index behaviour, 1.133 collapsing to 1.023, "
    "is EXPLICITLY OVERTURNED and must not be given as the rationale anywhere: the "
    "disqualifying issue is causation, not confounding.")
removed["thyroid_dis"] = {
    "label": "Thyroid disease",
    "status": "ordinary outcome, not a negative control",
    "removal_reason": ("Removed earlier, before 2026-08-07. Hypothyroidism can itself cause "
                       "disordered breathing during sleep."),
}

out = {
    "_what": ("The settled negative-control panel and the confounding floor. Read this file. Do "
              "not hardcode a control list or a floor number anywhere else."),
    "_built_by": "numbers/make_negcontrols_final.py",
    "_built_on": str(date.today()),
    "decision_date": "2026-08-07",
    "decided_by": "Alen Juginovic",
    "decision_record": "New_Figures/ALEN_COMMENTS_2026-08-07.md, final block",
    "qualifying_criterion": (
        "A negative control must have no plausible causal path FROM the exposure. This "
        "disqualifies a condition regardless of how clean its estimate looks, and it is the "
        "only criterion. Screening criteria applied to the replacements, all fixed before any "
        "model was fitted: at least 100 incident events in the 19,173-patient cohort, no "
        "plausible causal path from nocturnal hypoxemia, and not primarily a marker of "
        "healthcare contact intensity."),
    "n_controls": len(NEGATIVE_CONTROLS),
    "panel_keys": NEGATIVE_CONTROLS,
    "panel_labels": NEGATIVE_CONTROL_LABELS,
    "controls": controls,
    "confounding_floor_per_sd": round(floor, 3),
    "confounding_floor_driver": DISEASES[driver][0],
    "floor_definition": ("The largest negative-control hazard ratio on the paper's primary "
                         "per-SD scale: unpenalized Cox, stratified by site, adjusted for sex "
                         "and a 4-column natural cubic age spline with one column dropped so "
                         "the design is full rank, exposure rank-inverse-normalised within "
                         "site, oximetry_bad == 0, cohort 19,173."),
    "floor_history": [
        {"floor": 1.150, "driver": "Fracture", "panel_size": 5,
         "panel": ["Osteoarthritis", "Back pain", "Fracture", "Cataract", "Glaucoma"],
         "status": "superseded"},
        {"floor": 1.133, "driver": "Osteoarthritis", "panel_size": 4,
         "panel": ["Osteoarthritis", "Back pain", "Cataract", "Glaucoma"],
         "status": "superseded, fracture dropped alone"},
        {"floor": 1.063, "driver": "Glaucoma", "panel_size": 5,
         "panel": ["Back pain", "Cataract", "Glaucoma", "Contact dermatitis", "Hemorrhoids"],
         "status": "settled 2026-08-07 (v7 tables); the v7 panel as it was, recorded literally"},
        {"floor": 1.091, "driver": "Glaucoma", "panel_size": 5,
         "panel": ["Back pain", "Cataract", "Glaucoma", "Contact dermatitis", "Hemorrhoids"],
         "status": "v8.0 (2026-09-14 morning): the v7 panel refit on the v8 tables; cataract 1.067 and back pain 1.074 sat above 1"},
        {"floor": round(floor, 3), "driver": DISEASES[driver][0], "panel_size": len(NEGATIVE_CONTROLS),
         "panel": NEGATIVE_CONTROL_LABELS,
         "status": "v8.1 (2026-09-14, Alen's decision): cataract and back pain replaced by inguinal hernia and alopecia, "
                   "chosen from the pre-specified August candidate list; refit on the v8.1 tables (numbers/bdsp_diseases_v3.csv, "
                   "chain step 111), written by this script as chain step 205"},
    ],
    "old_vs_new": OLD_VS_NEW,
    "movement_is_downward_only": (
        "On the PRIMARY PER-SD SCALE the three v7 panel changes each moved the floor downward, 1.150 to 1.133 "
        "to 1.063, so on the v7 tables associations were promoted and none demoted. "
        + ((f"The v8 recalculation, with the panel unchanged, moved the floor "
            f"{'UP' if round(floor, 3) > _prev_floor else 'DOWN' if round(floor, 3) < _prev_floor else 'nowhere'} "
            f"from {_prev_floor} to {round(floor, 3)}; an association that cleared the earlier floor by less than "
            f"that move is not guaranteed to clear this one. The downward-only guarantee held for panel changes on "
            f"one set of tables, not across a recalculation. ") if _prev else "")
        + "This guarantee does NOT extend to the subgroup analyses, and assuming it does is a trap. See "
        "subgroup_floors_are_not_monotone below."),
    "subgroup_floors_are_not_monotone": {
        "_read_this": (
            "The duration-by-oxygen and stage-window analyses do not read the 1.063 floor. Each "
            "recomputes its own floor as a maximum over the five controls REFITTED inside that "
            "cell, using a TWO-SIDED rule, max(hr, 1/hr), the largest deviation from the null in "
            "either direction. Swapping two controls for two different ones therefore moves "
            "those floors in both directions, and on 2026-08-07 several rose."),
        "mechanism": (
            "Contact dermatitis carries 662 incident events cohort-wide against fracture's 846 "
            "and osteoarthritis's 2,092, and in the small cells it thins to single figures. In "
            "the low-oxygen short-sleep contrast it lands at hazard ratio 0.506 (95% CI, "
            "0.232-1.104) on 187 events spread across six cells with as few as 8 in one of them. "
            "A two-sided maximum reads 1/0.506 = 1.976 and turns that into the bar every real "
            "association must clear. A control estimate BELOW 1 is not evidence of upward "
            "confounding, so a two-sided rule lets sampling noise in a small control strike "
            "genuine findings."),
        "cells_where_the_floor_rose": {
            "dur_x_oxygen cell2_vs_cell1": {"old": 1.537, "old_driver": "Fracture",
                                            "new": 1.976, "new_driver": "Contact dermatitis",
                                            "struck": 2, "struck_now": 10},
            "dur_x_oxygen_6cell ox_short": {"old": 1.539, "new": 1.976,
                                            "new_driver": "Contact dermatitis"},
            "dur_x_oxygen_6cell SL": {"old": 1.415, "new": 1.495,
                                      "new_driver": "Contact dermatitis"},
            "dur_x_oxygen_healthy fourcell_cell2": {"old": 1.554, "new": 3.472,
                                                    "new_driver": "Hemorrhoids",
                                                    "note": "cell2 holds 99 patients"},
        },
        "recommendation": (
            "Do not change the floor rule without Alen. Two options, both defensible, and the "
            "choice is his. (a) Make the floor ONE-SIDED, max(hr), which is what the primary "
            "scale already does and what the word floor means. (b) Keep it two-sided but require "
            "a control to have a minimum number of events in the cell before it may set the "
            "floor; every control in these cells is already flagged unstable_lt5_events."),
    },
    "removed_from_panel": removed,
    "replacements_screened_from": 25,
    "replacement_screen": "numbers/negcontrols_v4.csv and numbers/negcontrols_v4.json",
    "source_of_truth": "numbers/disease_definitions.py, NEGATIVE_CONTROLS",
    "primary_estimates_from": "numbers/bdsp_diseases_v3.csv",
    "causal_tests_from": "numbers/causal_tests_all.csv",
    "cohort_n": int(COHORT_N),
    "_chain_step": "111b",
    "_built_at": time.strftime("%Y-%m-%d %H:%M:%S"),
    "_inputs": _INPUTS,
    "_input_vintage_note": (
        "confounding_floor_per_sd, confounding_floor_driver, controls.*.per_sd, events and at_risk come from "
        "bdsp_diseases_v3.csv (chain step 111, v8 tables) and controls.*.graded from results_v2.json (step 112). "
        "controls.*.landmark_2y and controls.*.bmi_adjusted come from causal_tests_all.csv, which has no v8 chain "
        f"step (file dated {_INPUTS['causal_tests_all.csv']['mtime_local']}); those two blocks are v7-era until it is "
        "regenerated. No chain script reads them (readers use panel_labels, confounding_floor_per_sd, "
        "confounding_floor_driver, controls.*.per_sd.hr, n_controls, decision_date)."),
    "ranking_outcome_set_note": (
        "Contact dermatitis (662 events) and hemorrhoids (1,445) clear the 150-event bar that "
        "defines the 197-measurement ranking's outcome set, so adding them to DISEASES would "
        "have taken it from 48 outcomes to 50 and moved every concordance gain behind Figure 1, "
        "eFigure 4, eFigure 5 and the total-sleep-time ceiling analysis. That is not implied by "
        "the negative-control decision, so the ranking is held at 48 through "
        "disease_definitions.RANKING_EXCLUDE. AWAITING ALEN'S CALL."),
}

path = PATH
json.dump(out, open(path, "w"), indent=1)
print("OLD VERSUS NEW (rule 11)")
print(f"  previous file   floor {_prev_floor} ({_prev_driver}), built {_prev.get('_built_on') if _prev else None}")
print(f"  module constant disease_definitions.CONFOUNDING_FLOOR_PER_SD = {CONFOUNDING_FLOOR_PER_SD}")
print(f"  new             floor {round(floor, 3)} ({DISEASES[driver][0]})   moved up: {OLD_VS_NEW['floor_moved_up']}"
      f"   dramatic (driver changed or >20%): {OLD_VS_NEW['dramatic']}")
if _prev:
    for k in hrs:
        _o = _prev["controls"].get(k, {}).get("per_sd", {}).get("hr")
        print(f"    {DISEASES[k][0]:<20} old {_o}   new {round(hrs[k], 3)}   "
              f"{'same' if _o == round(hrs[k], 3) else f'moved {round(hrs[k], 3) - _o:+.3f}' if _o is not None else 'new'}")
print(f"panel: {', '.join(NEGATIVE_CONTROL_LABELS)}")
for k in NEGATIVE_CONTROLS:
    p = controls[k]["per_sd"]
    print(f"  {controls[k]['label']:<20}{p['hr']:.3f} ({p['lo']:.3f}-{p['hi']:.3f})   "
          f"{controls[k]['events']:,} events")
print(f"floor {floor:.3f}, set by {DISEASES[driver][0]}")
print(f"written -> {path}")
