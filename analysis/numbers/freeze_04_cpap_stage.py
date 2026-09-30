"""
Phase 1.6 to 1.10. Freeze the stage-resolved, treatment and residual-risk numbers.

Same model specification as freeze_03. Everything recomputed from source.
Writes numbers/cpap_stage.json.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../..")))  # repository root, where paths.py lives
import paths

# ----------------------------------------------------------------------------- superseded specification
# This script fits the stage-resolved and treatment numbers at penalizer=0.02 and 0.05, the specification that run_treatment_v2.py (step 117) replaced: the printed values
# come from run_treatment_v2.py (step 117) (unpenalized), and nothing in the repository reads numbers/cpap_stage.json. It is kept for the
# old-versus-new record of the v8 recalculation and stops here unless a reviewer asks for it explicitly.
import os as _guard_os
if _guard_os.environ.get("T90_RUN_SUPERSEDED") != "1":
    raise SystemExit("freeze_04_cpap_stage: superseded penalized specification, not a source of any printed number. "
                     "Set T90_RUN_SUPERSEDED=1 to run it anyway.")
import warnings
warnings.filterwarnings("ignore")
import json
import numpy as np
import pandas as pd
from scipy import stats
from patsy import dmatrix
from lifelines import CoxPHFitter
import sys
sys.path.insert(0, paths.ANALYSIS_DIR)
from disease_definitions import DISEASES, NEGATIVE_CONTROLS

OUT = paths.NUMBERS_DIR
R = {}
S = ["wake", "sleep", "nrem", "n1", "n2", "n3", "rem"]
MIN_MIN = 30.0
CUT = 10.0

d = pd.read_parquet(f"{paths.TABLES_DIR}/cpap_t90_by_stage.parquet")
A = d[d.design == "A_split"].copy(); A["src"] = "A"
B = d[d.design == "B_pairs"].copy(); B["src"] = "B"
P = d[d.design == "P_placebo"].copy()
for s in S:
    for a, bb in [("pre", "dx"), ("post", "tx")]:
        for suf in ["t90", "min"]:
            c = f"{bb}_{s}_{suf}"
            B[f"{a}_{s}_{suf}"] = B[c] if c in B.columns else np.nan

R["designs"] = {"split_night": int(len(A)), "diagnostic_titration_pairs": int(len(B)),
                "placebo_splits": int(len(P))}

# ---------------------------------------------------------- 1.7 CPAP effect by stage
print("1.7  CPAP effect on T90 by stage")
eff = {}
for s in S:
    g = A[(A[f"pre_{s}_min"] >= MIN_MIN) & (A[f"post_{s}_min"] >= MIN_MIN)]
    g = g[g[f"pre_{s}_t90"].notna() & g[f"post_{s}_t90"].notna()]
    if len(g) < 30:
        eff[s] = {"n": int(len(g)), "note": "too few with the stage measurable both sides"}
        continue
    pre, post = g[f"pre_{s}_t90"], g[f"post_{s}_t90"]
    w = stats.wilcoxon(pre, post) if len(g) > 10 else None
    eff[s] = {"n": int(len(g)),
              "pre_median": round(float(pre.median()), 2),
              "post_median": round(float(post.median()), 2),
              "pre_mean": round(float(pre.mean()), 2),
              "post_mean": round(float(post.mean()), 2),
              "relative_reduction_pct": (round(100 * (1 - post.median() / pre.median()), 1)
                                         if pre.median() > 0 else None),
              "fell_in_pct": round(100 * float((post < pre).mean()), 1),
              "wilcoxon_p": (float(w.pvalue) if w is not None else None)}
R["cpap_effect_by_stage"] = eff

# ---------------------------------------------------------- 1.8 residual distribution
print("1.8  residual distribution on treatment")
dist = {}
for s in ["rem", "sleep", "nrem", "n2"]:
    g = A[A[f"post_{s}_min"] >= MIN_MIN]
    v = g[f"post_{s}_t90"].dropna()
    if len(v) < 100:
        continue
    bands = [("exactly_0", -1, 1e-9), ("0_to_1", 1e-9, 1), ("1_to_5", 1, 5),
             ("5_to_10", 5, 10), ("above_10", 10, 1e9)]
    dist[s] = {"n": int(len(v)),
               "median": round(float(v.median()), 2), "mean": round(float(v.mean()), 2),
               "p90": round(float(v.quantile(.90)), 1), "p95": round(float(v.quantile(.95)), 1),
               "bands": {nm: {"n": int(((v > lo) & (v <= hi)).sum()),
                              "pct": round(100 * float(((v > lo) & (v <= hi)).mean()), 1)}
                         for nm, lo, hi in bands}}
R["residual_distribution"] = dist

# ---------------------------------------------------------- 1.9 residual risk
print("1.9  residual hypoxemia and incident disease")
keep = ["BDSPPatientID", "SexDSC", "src"] + \
       [f"{p}_{s}_{x}" for p in ("pre", "post") for s in S for x in ("t90", "min")]
cc = pd.concat([A[[c for c in keep if c in A.columns]],
                B[[c for c in keep if c in B.columns]]], ignore_index=True)
cc = cc.sort_values("src", ascending=False).drop_duplicates("BDSPPatientID", keep="first")

OUTC = [("hf", "Heart failure"), ("cvd", "Cardiovascular composite"),
        ("sepsis", "Sepsis"), ("resp_failure", "Respiratory failure"),
        ("aki", "Acute kidney injury"), ("pneumonia", "Pneumonia"),
        ("death", "Death from any cause"),
        # the settled negative-control panel, imported not listed
        ] + [(k, DISEASES[k][0]) for k in NEGATIVE_CONTROLS]
CARD = ["hf", "cvd", "ihd", "mi", "afib", "pulm_htn", "cardiac_arrest", "pad", "stroke_any"]
need = ["BDSPPatientID", "site_id", "fu_valid", "AgeAtVisit"]
for k, _ in OUTC:
    need += [f"{k}_incident", f"{k}_years", f"{k}_prevalent"]
need += [f"{c}_prevalent" for c in CARD]
av = set(pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet").columns)
base = pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet",
                       columns=[c for c in dict.fromkeys(need) if c in av])
m = cc.merge(base, on="BDSPPatientID", how="inner")
m = m[(m.fu_valid == 1) & (m.post_sleep_min >= MIN_MIN)].copy()
m["male"] = (m.SexDSC.astype(str).str.upper().str[0] == "M").astype(int)
m["lbase"] = np.log1p(m.pre_sleep_t90)
m["nr"] = (m.post_sleep_t90 > CUT).astype(int)
cardcols = [f"{c}_prevalent" for c in CARD if f"{c}_prevalent" in m.columns]
m["anycard"] = (m[cardcols].fillna(0).sum(axis=1) > 0).astype(int)

def fit(g, k, xcol="nr"):
    f = g[(g[f"{k}_prevalent"] == 0) & g[f"{k}_years"].notna() & (g[f"{k}_years"] > 0)]
    ff = pd.DataFrame({"T": f[f"{k}_years"], "E": f[f"{k}_incident"].astype(int),
                       "x": f[xcol], "age": f.AgeAtVisit, "male": f.male,
                       "lbase": f.lbase, "site": f.site_id}).dropna()
    if ff.E.sum() < 40:
        return None
    c = CoxPHFitter(penalizer=0.02).fit(ff, "T", "E", strata=["site"])
    lo, hi = np.exp(c.confidence_intervals_.loc["x"])
    return {"n": int(len(ff)), "events": int(ff.E.sum()),
            "hr": round(float(np.exp(c.params_["x"])), 3), "lo": round(float(lo), 3),
            "hi": round(float(hi), 3), "p": float(c.summary.loc["x", "p"])}

clean = m[m.anycard == 0]
R["residual_risk"] = {
    "n_total": int(len(m)), "n_nonresponder": int(m.nr.sum()),
    "pct_nonresponder": round(100 * float(m.nr.mean()), 1),
    "n_cardiac_free": int(len(clean)),
    "prevalent_cardiac_responders_pct": round(100 * float(m.loc[m.nr == 0, "anycard"].mean()), 1),
    "prevalent_cardiac_nonresponders_pct": round(100 * float(m.loc[m.nr == 1, "anycard"].mean()), 1),
    "baseline_t90_responders_median": round(float(m.loc[m.nr == 0, "pre_sleep_t90"].median()), 2),
    "baseline_t90_nonresponders_median": round(float(m.loc[m.nr == 1, "pre_sleep_t90"].median()), 2),
    "full": {}, "cardiac_free": {}}
for k, lab in OUTC:
    a, c = fit(m, k), fit(clean, k)
    if a:
        R["residual_risk"]["full"][lab] = a
    if c:
        R["residual_risk"]["cardiac_free"][lab] = c

# ---------------------------------------------------------- dose-response by band
print("1.8b dose-response by residual band")
BANDS = [(-1, 1e-9), (1e-9, 1), (1, 5), (5, 10), (10, 1e9)]
LBL = ["exactly 0", "0-1%", "1-5%", "5-10%", ">10%"]
dr = {}
for stage in ["sleep", "rem"]:
    g = m[m[f"post_{stage}_min"] >= MIN_MIN].copy() if f"post_{stage}_min" in m.columns else None
    if g is None or len(g) < 500:
        continue
    g["band"] = g[f"post_{stage}_t90"].apply(
        lambda v: next((i for i, (lo, hi) in enumerate(BANDS) if lo < v <= hi), np.nan))
    g = g[g.band.notna()]; g["band"] = g.band.astype(int)
    dr[stage] = {"band_n": {LBL[i]: int((g.band == i).sum()) for i in range(5)}, "outcomes": {}}
    for k, lab in OUTC:
        f = g[(g[f"{k}_prevalent"] == 0) & g[f"{k}_years"].notna() & (g[f"{k}_years"] > 0)]
        X = pd.get_dummies(f.band, prefix="b").astype(int)
        for i in range(5):
            if f"b_{i}" not in X:
                X[f"b_{i}"] = 0
        X = X[[f"b_{i}" for i in range(5)]].drop(columns=["b_0"])
        ff = pd.concat([pd.DataFrame({"T": f[f"{k}_years"].values,
                                      "E": f[f"{k}_incident"].astype(int).values,
                                      "age": f.AgeAtVisit.values, "male": f.male.values,
                                      "lbase": f.lbase.values, "site": f.site_id.values}),
                        X.reset_index(drop=True)], axis=1).dropna()
        if ff.E.sum() < 60:
            continue
        try:
            c = CoxPHFitter(penalizer=0.05).fit(ff, "T", "E", strata=["site"])
        except Exception:
            continue
        cells = {}
        for i in range(1, 5):
            nm = f"b_{i}"
            if nm in c.params_.index:
                lo, hi = np.exp(c.confidence_intervals_.loc[nm])
                cells[LBL[i]] = {"hr": round(float(np.exp(c.params_[nm])), 3),
                                 "lo": round(float(lo), 3), "hi": round(float(hi), 3),
                                 "p": float(c.summary.loc[nm, "p"])}
        # trend
        ft = pd.DataFrame({"T": f[f"{k}_years"].values,
                           "E": f[f"{k}_incident"].astype(int).values,
                           "band": f.band.values, "age": f.AgeAtVisit.values,
                           "male": f.male.values, "lbase": f.lbase.values,
                           "site": f.site_id.values}).dropna()
        try:
            ct = CoxPHFitter(penalizer=0.05).fit(ft, "T", "E", strata=["site"])
            cells["trend_p"] = float(ct.summary.loc["band", "p"])
        except Exception:
            pass
        dr[stage]["outcomes"][lab] = {"events": int(ff.E.sum()), **cells}
R["dose_response"] = dr

with open(f"{OUT}/cpap_stage.json", "w") as f:
    json.dump(R, f, indent=2)
print(f"\nwritten -> {OUT}/cpap_stage.json")
print("designs:", R["designs"])
print("residual risk n:", R["residual_risk"]["n_total"],
      " non-responders:", R["residual_risk"]["n_nonresponder"])
