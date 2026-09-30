"""
Step 220 run_ph_all_v8 (decision 12): the proportional-hazards test for every outcome, and early
and late follow-up hazard ratios for the violators. Supplies the writer that ph_check_v2.json and
ph_timesplit_v2.json never had (v7 GROUP_E, blocker B7).

Model: the primary per-SD model of run_primary_unpenalized.py, exactly (age as cr(df=4) with the
first column dropped, sex, site strata, penalizer 0, the exposure as the within-site rank inverse
normal of cohort_spec.t90_column(), risk set free of the outcome at the study with follow-up > 0,
event floor 60). Test: lifelines proportional_hazard_test on the fitted model with the rank
time transform; the p value of the exposure term decides (violated at p < 0.05), the smallest p
over every term is reported beside it. Violators are refitted on 0 to 2 years and beyond 2 years
(left-truncated Cox, the same covariates), and on 0 to 1, 1 to 3 and 3 or more years as extra
columns; a split needs 30 events.

Positive control (smoke on the v7 tables): the three conditions ph_check_v2.json called violated
(respiratory failure, heart failure, COPD) must be violated here too, and their 0 to 2 and beyond 2
year hazard ratios must sit within 0.05 of ph_timesplit_v2.json. Those reference files were built
with a penalizer of 0.01, so the comparison also refits the violators at that penalizer and reports
both distances; the printed result is the unpenalized one. The reference is read from the file,
never typed.

Outputs (cohort_spec.OUT_DIR): ph_all_v8.json, ph_all_v8.csv (the supplementary table), and the two
legacy names in their old schema, ph_check_v2.json (label -> events, p, verdict, every outcome) and
ph_timesplit_v2.json (label -> "0 to 2 y", "beyond 2 y" blocks for every violator).
"""
import json
import os
import sys
import warnings

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test
from patsy import dmatrix
from scipy import stats

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from cohort_spec import (apply_cohort, T90_FINAL, OUT_DIR, NUMBERS_DIR, PRIMARY_EXPOSURE,  # noqa: E402
                         t90_column, exposure_tag, sidecar)
from disease_definitions import DISEASES, NEGATIVE_CONTROLS  # noqa: E402

EVENT_FLOOR, SPLIT_FLOOR, ALPHA = 60, 30, 0.05
XCOL = t90_column()
SFX = exposure_tag(XCOL) if XCOL != PRIMARY_EXPOSURE else ""

b = apply_cohort(pd.read_parquet(T90_FINAL))
b["male"] = (b.sex.astype(str).str.upper().str[0] == "M").astype(int)
sp = dmatrix("cr(a, df=4) - 1", {"a": b.AgeAtVisit.values}, return_type="dataframe").iloc[:, 1:]
for i in range(sp.shape[1]):
    b[f"age_s{i}"] = sp.iloc[:, i].values
ADJ = [f"age_s{i}" for i in range(sp.shape[1])] + ["male"]


def rint(x):
    x = pd.Series(x)
    out = pd.Series(np.nan, index=x.index)
    ok = x.notna()
    r = stats.rankdata(x[ok])
    out[ok] = stats.norm.ppf((r - 0.375) / (ok.sum() + 0.25))
    return out


b["z"] = b.groupby("site_id")[XCOL].transform(rint)
print(f"cohort {len(b):,}   exposure {XCOL}", flush=True)


def frame(key):
    g = b[(b[f"{key}_prevalent"] == 0) & b[f"{key}_years"].notna() & (b[f"{key}_years"] > 0)]
    return pd.DataFrame({"T": g[f"{key}_years"], "E": g[f"{key}_incident"].astype(int), "site": g.site_id,
                         "z": g.z, **{c: g[c] for c in ADJ}}).dropna()


def hr(m, x="z"):
    lo, hi = np.exp(m.confidence_intervals_.loc[x])
    return float(np.exp(m.params_[x])), float(lo), float(hi), float(m.summary.loc[x, "p"])


def split_fit(d, t0, t1, pen=0.0):
    """Cox restricted to (t0, t1], left-truncated at t0, same covariates."""
    e = d[d["T"] > t0].copy()
    e["E2"] = np.where(e["T"] > t1, 0, e["E"])
    e["T2"] = np.minimum(e["T"], t1)
    e["start"] = float(t0)
    e = e[e["T2"] > t0]
    if e.E2.sum() < SPLIT_FLOOR:
        return None
    m = CoxPHFitter(penalizer=pen).fit(e[["start", "T2", "E2", "site", "z"] + ADJ], "T2", "E2",
                                       entry_col="start", strata=["site"])
    h, lo, hi, p = hr(m)
    return {"events": int(e.E2.sum()), "n": int(len(e)), "hr": h, "lo": lo, "hi": hi, "p": p}  # v8.3 full precision (2026-09-16)


WINDOWS = {"split_0_2": (0, 2), "split_2_inf": (2, 1e9), "split_0_1": (0, 1), "split_1_3": (1, 3), "split_3_inf": (3, 1e9)}
rows = []
for key, (label, _neg, _, _) in list(DISEASES.items()) + [("death", ("Death from any cause", False, [], []))]:
    if f"{key}_incident" not in b.columns:
        continue
    d = frame(key)
    if d.E.sum() < EVENT_FLOOR:
        continue
    m = CoxPHFitter(penalizer=0.0).fit(d, "T", "E", strata=["site"])
    h, lo, hi, p = hr(m)
    try:
        pt = proportional_hazard_test(m, d, time_transform="rank")
        s = pt.summary
        p_x = float(s.loc[("z", "rank"), "p"]) if ("z", "rank") in s.index else float(s.loc["z", "p"])
        p_min = float(np.nanmin(s["p"].values))
        worst = str(s["p"].idxmin())
    except Exception as e:  # pragma: no cover
        p_x, p_min, worst = np.nan, np.nan, f"ERR {e}"
    rec = {"key": key, "disease": label, "negative_control": key in NEGATIVE_CONTROLS,
           "events": int(d.E.sum()), "n": int(len(d)),
           "hr_per_sd": h, "lo": lo, "hi": hi, "p": p,  # v8.3 full precision (2026-09-16)
           "schoenfeld_p_exposure": p_x, "schoenfeld_p_min_any_term": p_min, "schoenfeld_worst_term": worst,
           "violated": bool(p_x < ALPHA)}
    if rec["violated"]:
        for name, (t0, t1) in WINDOWS.items():
            r = split_fit(d, t0, t1)
            for k2 in ("events", "hr", "lo", "hi", "p"):
                rec[f"{name}_{k2}"] = r[k2] if r else np.nan
    rows.append(rec)
    print(f"{label:<30} ev={rec['events']:>5} HR={h:.3f} schoenfeld p={p_x:.4g} {'VIOLATED' if rec['violated'] else ''}", flush=True)

df = pd.DataFrame(rows)
nv = int(df.violated.sum())
print(f"\n{nv} of {len(df)} outcomes violate proportional hazards for the exposure at p < {ALPHA}")

# ------------------------------------------------------------------ positive control against the legacy files
control = {"reference_files": [], "checked": [], "note": ""}
ref_check, ref_split = f"{NUMBERS_DIR}/ph_check_v2.json", f"{NUMBERS_DIR}/ph_timesplit_v2.json"
if os.path.exists(ref_check) and os.path.exists(ref_split):
    RC, RS = json.load(open(ref_check)), json.load(open(ref_split))
    control["reference_files"] = [ref_check, ref_split]
    lab2key = {v[0]: k for k, v in DISEASES.items()}
    lab2key["Death"] = "death"
    lab2key["Death from any cause"] = "death"
    for lab, v in RC.items():
        key = lab2key.get(lab)
        mine = df[df.key == key]
        if key is None or mine.empty:
            control["checked"].append({"label": lab, "status": "not in this run"})
            continue
        mine = mine.iloc[0]
        entry = {"label": lab, "ref_verdict": v["verdict"], "verdict_here": "violated" if mine.violated else "satisfied",
                 "ref_events": v["events"], "events_here": int(mine.events), "ref_p": v["p"], "p_here": float(mine.schoenfeld_p_exposure)}
        if lab in RS:
            d = frame(key)
            dist_unpen, dist_pen = [], []
            for win, (t0, t1) in (("0 to 2 y", (0, 2)), ("beyond 2 y", (2, 1e9))):
                ref = RS[lab][win]["hr"]
                u = split_fit(d, t0, t1, 0.0)
                q = split_fit(d, t0, t1, 0.01)
                dist_unpen.append(abs(u["hr"] - ref) if u else np.nan)
                dist_pen.append(abs(q["hr"] - ref) if q else np.nan)
                entry[f"{win} ref"] = ref
                entry[f"{win} here unpenalized"] = u["hr"] if u else None
                entry[f"{win} here penalized 0.01"] = q["hr"] if q else None
            entry["max_abs_diff_unpenalized"] = float(np.nanmax(dist_unpen))
            entry["max_abs_diff_penalized_0.01"] = float(np.nanmax(dist_pen))
            entry["within_0.05"] = bool(min(entry["max_abs_diff_unpenalized"], entry["max_abs_diff_penalized_0.01"]) <= 0.05)
        control["checked"].append(entry)
    viol_ref = {lab for lab, v in RC.items() if v["verdict"] == "violated"}
    viol_here = {c["label"] for c in control["checked"] if c.get("verdict_here") == "violated"}
    control["reference_violators_reproduced"] = sorted(viol_ref & viol_here)
    control["reference_violators_missed"] = sorted(viol_ref - viol_here)
    control["split_within_0.05"] = [c["label"] for c in control["checked"] if c.get("within_0.05")]
    control["split_outside_0.05"] = [c["label"] for c in control["checked"] if "within_0.05" in c and not c["within_0.05"]]
    print("\npositive control against ph_check_v2.json / ph_timesplit_v2.json:")
    print(f"  reference violators reproduced: {control['reference_violators_reproduced']}  missed: {control['reference_violators_missed']}")
    for c in control["checked"]:
        if "within_0.05" in c:
            print(f"  {c['label']:<22} 0-2 y ref {c['0 to 2 y ref']:.3f} here {c['0 to 2 y here unpenalized']:.3f} (pen {c['0 to 2 y here penalized 0.01']:.3f})"
                  f"   2+ y ref {c['beyond 2 y ref']:.3f} here {c['beyond 2 y here unpenalized']:.3f} (pen {c['beyond 2 y here penalized 0.01']:.3f})"
                  f"   within 0.05: {c['within_0.05']}")
    control["note"] = ("the reference files were fitted with penalizer 0.01 on an earlier cohort; the printed result is the "
                       "unpenalized primary model, the penalized refit is the reproduction control")
else:
    control["note"] = "legacy reference files absent, no comparison"

# ------------------------------------------------------------------ outputs
csv_path, json_path = f"{OUT_DIR}/ph_all_v8{SFX}.csv", f"{OUT_DIR}/ph_all_v8{SFX}.json"
df.to_csv(csv_path, index=False)
out = {"exposure": XCOL, "model": "primary per-SD model (cr(age, df=4) minus one column, sex, site strata, penalizer 0)",
       "test": "lifelines proportional_hazard_test, time_transform rank, exposure term", "alpha": ALPHA,
       "event_floor": EVENT_FLOOR, "split_event_floor": SPLIT_FLOOR, "windows_years": {k: list(v) for k, v in WINDOWS.items()},
       "n_outcomes": int(len(df)), "n_violated": nv, "violators": df[df.violated].disease.tolist(),
       "outcomes": {r["disease"]: {k: (None if (isinstance(v, float) and np.isnan(v)) else v) for k, v in r.items() if k != "disease"}
                    for r in df.to_dict("records")},
       "positive_control": control}
json.dump(out, open(json_path, "w"), indent=1, default=float)
written = [csv_path, json_path]
if not SFX:
    # the legacy names, same schema as the April files, now written by code (blocker B7)
    chk = {r["disease"]: {"events": int(r["events"]), "p": round(float(r["schoenfeld_p_exposure"]), 4),
                          "verdict": "violated" if r["violated"] else "satisfied"} for r in rows}
    chk["Death"] = chk.pop("Death from any cause")
    ts = {}
    for r in rows:
        if r["violated"]:
            ts[r["disease"] if r["disease"] != "Death from any cause" else "Death"] = {
                "0 to 2 y": {"events": r["split_0_2_events"], "hr": r["split_0_2_hr"], "lo": r["split_0_2_lo"], "hi": r["split_0_2_hi"]},
                "beyond 2 y": {"events": r["split_2_inf_events"], "hr": r["split_2_inf_hr"], "lo": r["split_2_inf_lo"], "hi": r["split_2_inf_hi"]}}
    p1, p2 = f"{OUT_DIR}/ph_check_v2.json", f"{OUT_DIR}/ph_timesplit_v2.json"
    json.dump(chk, open(p1, "w"), indent=1, default=float)
    json.dump(ts, open(p2, "w"), indent=1, default=float)
    written += [p1, p2]
for p in written:
    sidecar(p, __file__, extra_inputs=[x for x in (ref_check, ref_split) if os.path.exists(x)],
            note=f"run_ph_all_v8: {len(df)} outcomes, {nv} violators, exposure {XCOL}")
print("written:", ", ".join(written))
