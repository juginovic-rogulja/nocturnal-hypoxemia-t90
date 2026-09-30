"""
Step 211 run_treatment_v3 (decisions 6 and 7): the PAP responder comparison of run_treatment_v2.py
with an adjustment for prevalent cardiopulmonary disease, and Benjamini-Hochberg q across the
non-control conditions.

The frame, the filters, the oximetry rule, the two groups (everyone started with sleep T90 above
10 percent; not corrected = still above 10 percent on PAP) and the unadjusted model (age, sex,
log1p pretreatment sleep T90, site strata, penalizer 0, 40-event floor) are those of
run_treatment_v2.py, reproduced here line for line. Two things are added:

  adjusted   the same Cox model plus cp_any, 1 when the patient carried any prevalent condition
             of cohort_spec.PREVALENT_CARDIOPULMONARY (D4 default: heart failure, ischemic heart
             disease, myocardial infarction, atrial fibrillation, stroke, peripheral artery disease,
             pulmonary hypertension, venous thromboembolism, COPD, asthma, respiratory failure,
             obesity hypoventilation; hypertension left out), the outcome under study excluded from
             its own covariate. No adherence data exist and none is modelled.
  q          Benjamini-Hochberg over the non-control conditions with at least 40 events, one family
             per contrast (unadjusted and adjusted); the negative controls corrected among
             themselves. The count of q < 0.05 is what replaces "sixteen conditions".

Positive control inside the script: with the adjustment off, every hazard ratio, CI bound (to the
3 dp the file stores) and p (to 1e-9) equals numbers/treatment_v2.json, and the group counts equal
its corrected_vs_not block. The script stops otherwise. Tables and paths come from cohort_spec only.

Output: cohort_spec.OUT_DIR/treatment_v3.json.
"""
import json
import os
import sys
import warnings

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from cohort_spec import (CPAP_STAGE_PARQUET, T90_FINAL, NUMBERS_DIR, OUT_DIR,  # noqa: E402
                         PREVALENT_CARDIOPULMONARY, sidecar)
from disease_definitions import DISEASES, NEGATIVE_CONTROLS  # noqa: E402

EVENT_FLOOR = 40
S = ["wake", "sleep", "nrem", "n1", "n2", "n3", "rem"]
d = pd.read_parquet(CPAP_STAGE_PARQUET)
A = d[d.design == "A_split"].copy(); A["src"] = "A"
B = d[d.design == "B_pairs"].copy(); B["src"] = "B"
P = d[d.design == "P_placebo"]
for s in S:
    for a, bb in [("pre", "dx"), ("post", "tx")]:
        for suf in ["t90", "min"]:
            c = f"{bb}_{s}_{suf}"
            B[f"{a}_{s}_{suf}"] = B[c] if c in B.columns else np.nan
keep = ["BDSPPatientID", "SexDSC", "src"] + [f"{p}_{s}_{x}" for p in ("pre", "post") for s in S for x in ("t90", "min")]
cc = pd.concat([A[[c for c in keep if c in A.columns]], B[[c for c in keep if c in B.columns]]], ignore_index=True)
cc = cc.sort_values("src", ascending=False).drop_duplicates("BDSPPatientID", keep="first")

base = pd.read_parquet(T90_FINAL)
m = cc.merge(base, on="BDSPPatientID", how="inner")
m = m[(m.fu_valid == 1) & (m.post_sleep_min >= 30) & (m.pre_sleep_min >= 30)].copy()
for _side in ("pre", "post"):
    _n, _mn = f"{_side}_sleep_nadir", f"{_side}_sleep_mean"
    if _n in m.columns and _mn in m.columns:
        m = m[~(m[_n] > m[_mn])]
    _t = f"{_side}_sleep_t90"
    if _t in m.columns and _mn in m.columns:
        m = m[~(m[_mn].between(50, 80) & (m[_t] > 70))]
m["male"] = (m.SexDSC.astype(str).str.upper().str[0] == "M").astype(int)
m["lbase"] = np.log1p(m.pre_sleep_t90)

sev = m[m.pre_sleep_t90 > 10].copy()
sev["notcorr"] = (sev.post_sleep_t90 > 10).astype(int)
CP = [k for k in PREVALENT_CARDIOPULMONARY if f"{k}_prevalent" in sev.columns]
cp_any_all = sev[[f"{k}_prevalent" for k in CP]].fillna(0).max(axis=1).astype(int)

T2 = json.load(open(f"{NUMBERS_DIR}/treatment_v2.json"))
CVN = T2["corrected_vs_not"]
counts = {"n": int(len(sev)), "n_corrected": int((1 - sev.notcorr).sum()), "n_not": int(sev.notcorr.sum()),
          "baseline_corrected": round(float(sev.loc[sev.notcorr == 0, "pre_sleep_t90"].median()), 1),
          "baseline_not": round(float(sev.loc[sev.notcorr == 1, "pre_sleep_t90"].median()), 1),
          "n_prevalent_cardiopulmonary_any": int(cp_any_all.sum()),
          "pct_prevalent_cardiopulmonary_corrected": round(100 * float(cp_any_all[sev.notcorr == 0].mean()), 1),
          "pct_prevalent_cardiopulmonary_not": round(100 * float(cp_any_all[sev.notcorr == 1].mean()), 1)}
for k in ("n", "n_corrected", "n_not", "baseline_corrected", "baseline_not"):
    assert counts[k] == CVN[k], f"{k}: here {counts[k]} against treatment_v2.json {CVN[k]}"
print(f"all started >10%: n={counts['n']:,}  corrected {counts['n_corrected']:,}  not {counts['n_not']:,}   "
      f"prevalent cardiopulmonary {counts['n_prevalent_cardiopulmonary_any']:,}  (counts equal treatment_v2.json)")


def bh(pvals):
    p = np.asarray(pvals, float)
    q = np.full(len(p), np.nan)
    ok = np.where(~np.isnan(p))[0]
    n = len(ok)
    if n == 0:
        return q
    order = ok[np.argsort(p[ok])]
    prev = 1.0
    for rank, i in enumerate(reversed(order), start=1):
        val = p[i] * n / (n - rank + 1)
        prev = min(prev, val)
        q[i] = min(prev, 1.0)
    return q


def est(c, x="x"):
    lo, hi = np.exp(c.confidence_intervals_.loc[x])
    h = float(np.exp(c.params_[x]))
    return {"hr": h, "lo": float(lo), "hi": float(hi), "p": float(c.summary.loc[x, "p"]),  # v8.3 full precision (2026-09-16)
            "pct_lower_if_corrected": 100 * (1 - 1 / h) if h > 0 else None}


rows, mism = [], []
for k in list(DISEASES.keys()) + ["death"]:
    if f"{k}_incident" not in sev.columns:
        continue
    lab = DISEASES[k][0] if k in DISEASES else "Death from any cause"
    neg = k in NEGATIVE_CONTROLS
    f = sev[(sev[f"{k}_prevalent"] == 0) & sev[f"{k}_years"].notna() & (sev[f"{k}_years"] > 0)]
    others = [f"{j}_prevalent" for j in CP if j != k]
    ff = pd.DataFrame({"T": f[f"{k}_years"], "E": f[f"{k}_incident"].astype(int), "x": f.notcorr,
                       "age": f.AgeAtVisit, "male": f.male, "lbase": f.lbase, "site": f.site_id,
                       "cp_any": f[others].fillna(0).max(axis=1).astype(int) if others else 0}).dropna()
    if ff.E.sum() < EVENT_FLOOR:
        assert lab not in CVN["outcomes"], f"{lab} is under the floor here but present in treatment_v2.json"
        continue
    try:
        c0 = CoxPHFitter(penalizer=0.0).fit(ff[["T", "E", "x", "age", "male", "lbase", "site"]], "T", "E", strata=["site"])
    except Exception:
        continue
    u = est(c0)
    pub = CVN["outcomes"].get(lab)
    # v8.3: reference now full precision, compare within 5e-4 (2026-09-16)
    ok = pub is not None and all(abs(u[k] - pub[k]) <= 5e-4 for k in ("hr", "lo", "hi")) \
        and int(ff.E.sum()) == pub["events"] and abs(u["p"] - pub["p"]) < 1e-9
    if not ok:
        mism.append({"outcome": lab, "here": u, "events": int(ff.E.sum()), "published": pub})
    c1 = CoxPHFitter(penalizer=0.0).fit(ff[["T", "E", "x", "age", "male", "lbase", "cp_any", "site"]], "T", "E", strata=["site"])
    adj = est(c1)
    cpe = est(c1, "cp_any")
    rows.append({"key": k, "outcome": lab, "negative_control": neg, "events": int(ff.E.sum()), "n": int(len(ff)),
                 "n_cp_any": int(ff.cp_any.sum()), "events_cp_any": int(ff.loc[ff.cp_any == 1, "E"].sum()),
                 "hr": u["hr"], "lo": u["lo"], "hi": u["hi"], "p": u["p"], "pct_lower_if_corrected": u["pct_lower_if_corrected"],
                 "adj_hr": adj["hr"], "adj_lo": adj["lo"], "adj_hi": adj["hi"], "adj_p": adj["p"],
                 "adj_pct_lower_if_corrected": adj["pct_lower_if_corrected"],
                 "cp_any_hr": cpe["hr"], "cp_any_lo": cpe["lo"], "cp_any_hi": cpe["hi"], "cp_any_p": cpe["p"]})

df = pd.DataFrame(rows)
print(f"positive control: {len(df)} outcomes fitted, {len(df) - len(mism)} reproduce treatment_v2.json within 5e-4, {len(mism)} differ")
for x in mism[:10]:
    print("   MISMATCH", x)
assert not mism, "the unadjusted contrast does not reproduce treatment_v2.json; stop"
assert len(df) == len(CVN["outcomes"]) and set(df.outcome) == set(CVN["outcomes"]), (len(df), len(CVN["outcomes"]))

for pc, qc in (("p", "q"), ("adj_p", "adj_q")):
    df[qc] = np.nan
    for flag in (False, True):
        idx = df.index[df.negative_control == flag]
        df.loc[idx, qc] = bh(df.loc[idx, pc].values)
nc = ~df.negative_control
summary = {"n_conditions_tested": int(nc.sum()), "n_negative_controls": int((~nc).sum()),
           "unadjusted": {"n_p_lt_0.05": int(((df.p < 0.05) & nc).sum()), "n_q_lt_0.05": int(((df.q < 0.05) & nc).sum()),
                          "n_hr_gt_1": int(((df.hr > 1) & nc).sum())},
           "adjusted": {"n_p_lt_0.05": int(((df.adj_p < 0.05) & nc).sum()), "n_q_lt_0.05": int(((df.adj_q < 0.05) & nc).sum()),
                        "n_hr_gt_1": int(((df.adj_hr > 1) & nc).sum())},
           "conditions_q_lt_0.05_unadjusted": df[(df.q < 0.05) & nc].sort_values("hr", ascending=False).outcome.tolist(),
           "conditions_q_lt_0.05_adjusted": df[(df.adj_q < 0.05) & nc].sort_values("adj_hr", ascending=False).outcome.tolist(),
           "median_abs_log_hr_shift_from_adjustment": float(np.median(np.abs(np.log(df.adj_hr) - np.log(df.hr))))}  # v8.3 full precision (2026-09-16)
R = {"designs": T2["designs"], "counts": counts,
     "covariate": {"name": "cp_any", "definition": f"any prevalent condition among {CP} at the sleep study, the outcome under study excluded",
                   "alternative_not_run": "the same set plus hypertension (about half the cohort)", "adherence": "no adherence data exist"},
     "multiple_testing": "Benjamini-Hochberg over the non-control conditions with at least 40 events, one family per contrast; controls corrected among themselves",
     "event_floor": EVENT_FLOOR, "summary": summary,
     "positive_control": {"reference": f"{NUMBERS_DIR}/treatment_v2.json", "n_outcomes": int(len(df)), "n_mismatch": len(mism)},
     "outcomes": {r["outcome"]: {k: (None if isinstance(v, float) and np.isnan(v) else v) for k, v in r.items() if k != "outcome"}
                  for r in df.to_dict("records")}}
out = f"{OUT_DIR}/treatment_v3.json"
json.dump(R, open(out, "w"), indent=2, default=float)
sidecar(out, __file__, extra_inputs=[f"{NUMBERS_DIR}/treatment_v2.json"],
        note=f"run_treatment_v3: {len(df)} outcomes, q<0.05 unadjusted {summary['unadjusted']['n_q_lt_0.05']}, adjusted {summary['adjusted']['n_q_lt_0.05']}")

print(f"\nnon-control conditions: p<0.05 {summary['unadjusted']['n_p_lt_0.05']} -> q<0.05 {summary['unadjusted']['n_q_lt_0.05']} (unadjusted); "
      f"p<0.05 {summary['adjusted']['n_p_lt_0.05']} -> q<0.05 {summary['adjusted']['n_q_lt_0.05']} (adjusted for prevalent cardiopulmonary disease)")
print(f"{'condition':<28}{'ev':>6}{'HR':>8}{'q':>9}{'adj HR':>9}{'adj q':>9}{'cp_any HR':>11}")
for _, r in df[nc].sort_values("adj_hr", ascending=False).head(20).iterrows():
    print(f"  {r.outcome:<26}{r.events:>6}{r.hr:>8.2f}{r.q:>9.3f}{r.adj_hr:>9.2f}{r.adj_q:>9.3f}{r.cp_any_hr:>11.2f}")
print("written:", out)
