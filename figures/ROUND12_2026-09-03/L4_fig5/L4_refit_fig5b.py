"""
L4 (round 12, 2026-09-03): Figure 5b refit against the <=1%-on-PAP reference.

Step 1, positive control. Rebuild the treatment cohort exactly as numbers/run_treatment_v2.py
does (same parquet files, same filters, same oximetry quality rule, same covariates, same
unpenalized Cox with site strata) and refit the published two-group contrast (residual sleep
T90 on PAP above 10% against restored to 10% or less). Every outcome must reproduce the
frozen numbers/treatment_v2.json hazard ratio and 95% CI to the third decimal before the
reference is changed. The group counts are checked against treatment_v2.json as well.

Step 2, the new exposure. Same cohort (everyone started above 10%), same adjustment
(log1p pretreatment sleep T90, age, sex, site strata), but the on-PAP exposure is now three
bands: reference = sleep T90 on PAP of 1% or less, A = above 1% up to 10%, B = above 10%.
One Cox fit per outcome with two dummies, HR, 95% CI and P for A and B, plus a 2-df
likelihood-ratio P for the exposure as a whole.

Note on age. The task brief describes the adjustment as an age spline (cr, df 4, one column
dropped). The published script adjusts for age as a single linear term ("age": AgeAtVisit),
and that is what the frozen 5b numbers reproduce from. The primary refit therefore keeps the
linear term so the new panel is on exactly the same footing as the current one (which stays
as Extended Data Fig. 8). The spline version is run as a sensitivity for both contrasts and
written beside it, never mixed into the drawn values.

v8 (2026-09-12, decisions 6 and 7), two flags, both off by default so step 155 reproduces the
round-12 output unchanged:
  --adjust-prevalent-cardiopulmonary   every contrast is refitted a second time with the
        decision 6 covariate: any prevalent condition among cohort_spec.PREVALENT_CARDIOPULMONARY
        (the outcome under study excluded from the set, its risk set is free of it anyway).
        Columns adj_pub_*, adj_A_*, adj_B_*, adj_p_lr_2df. No adherence data exist.
  --bh  Benjamini-Hochberg q per band across the non-control conditions (pub_q, A_q, B_q and
        the adj_ twins); the negative controls get their own q across themselves, so a control
        still has a q a figure can test.
  --out-dir DIR (default this folder), --suffix S (default "_v3" when either flag is on).
The paths come from cohort_spec (CPAP_STAGE_PARQUET, T90_FINAL, NUMBERS_DIR). The literal count
asserts of round 12 became reads of treatment_v2.json (gate 6: no result literal in an assert).

Outputs (nothing written into numbers/ or the frozen tables):
  L4_fig5b_refit{suffix}.csv     all outcomes, both contrasts, primary and spline sensitivity
  L4_fig5b_refit{suffix}.json    the same, keyed by outcome label, plus counts and the gate record
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import argparse
import json
import os
import sys
import warnings

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from patsy import dmatrix

warnings.filterwarnings("ignore")
ROOT = paths.FIGURE_ROOT
NUM = f"{paths.NUMBERS_DIR}"
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, paths.ANALYSIS_DIR)
from disease_definitions import DISEASES, NEGATIVE_CONTROLS  # noqa: E402
from cohort_spec import (CPAP_STAGE_PARQUET, T90_FINAL, NUMBERS_DIR,  # noqa: E402
                         PREVALENT_CARDIOPULMONARY, sidecar)

ap = argparse.ArgumentParser()
ap.add_argument("--adjust-prevalent-cardiopulmonary", action="store_true")
ap.add_argument("--bh", action="store_true")
ap.add_argument("--out-dir", default=HERE)
ap.add_argument("--suffix", default=None)
ARGS = ap.parse_args()
ADJUST, BH = ARGS.adjust_prevalent_cardiopulmonary, ARGS.bh
SFX = ARGS.suffix if ARGS.suffix is not None else ("_v3" if (ADJUST or BH) else "")
os.makedirs(ARGS.out_dir, exist_ok=True)

# ------------------------------------------------------------------ cohort, verbatim path
S = ["wake", "sleep", "nrem", "n1", "n2", "n3", "rem"]
d = pd.read_parquet(CPAP_STAGE_PARQUET)
A_ = d[d.design == "A_split"].copy(); A_["src"] = "A"
B_ = d[d.design == "B_pairs"].copy(); B_["src"] = "B"
for s in S:
    for a, b in [("pre", "dx"), ("post", "tx")]:
        for suf in ["t90", "min"]:
            c = f"{b}_{s}_{suf}"
            B_[f"{a}_{s}_{suf}"] = B_[c] if c in B_.columns else np.nan
keep = ["BDSPPatientID", "SexDSC", "src"] + [f"{p}_{s}_{x}" for p in ("pre", "post")
                                             for s in S for x in ("t90", "min")]
cc = pd.concat([A_[[c for c in keep if c in A_.columns]],
                B_[[c for c in keep if c in B_.columns]]], ignore_index=True)
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

# the age spline for the sensitivity: cr(df=4) on the analysis cohort, first column dropped
# (the unpenalized-fit convention of numbers/run_primary_unpenalized.py)
_sp = dmatrix("cr(a, df=4) - 1", {"a": m.AgeAtVisit.values}, return_type="dataframe").iloc[:, 1:]
for i in range(_sp.shape[1]):
    m[f"age_s{i}"] = _sp.iloc[:, i].values
SPL = [f"age_s{i}" for i in range(_sp.shape[1])]
assert len(SPL) == 3, SPL

# decision 6 covariate: any prevalent cardiopulmonary condition (D4 set), computed per outcome
# below so the outcome under study is excluded from its own covariate
CP = [k for k in PREVALENT_CARDIOPULMONARY if f"{k}_prevalent" in m.columns]
# the prevalent flags are read as they are (a missing flag never enters a risk set); the
# covariate fills missing with 0 inside its own computation only

# ------------------------------------------------------------------ bands and counts
sev = m[m.pre_sleep_t90 > 10].copy()
sev["notcorr"] = (sev.post_sleep_t90 > 10).astype(int)
sev["band3"] = np.select([sev.post_sleep_t90 <= 1, sev.post_sleep_t90 <= 10], [0, 1], 2)
sev["xA"] = (sev.band3 == 1).astype(int)
sev["xB"] = (sev.band3 == 2).astype(int)
BANDS4 = [(-1, 1), (1, 5), (5, 10), (10, 1e9)]
cnt4 = [int(((sev.post_sleep_t90 > lo) & (sev.post_sleep_t90 <= hi)).sum()) for lo, hi in BANDS4]
COUNTS = {"n_all": int(len(sev)), "n_restored_le10": int((sev.notcorr == 0).sum()),
          "n_still_gt10": int(sev.notcorr.sum()), "bands4_on_pap": cnt4,
          "n_ref_le1": int((sev.band3 == 0).sum()), "n_A_1_to_10": int((sev.band3 == 1).sum()),
          "n_B_gt10": int((sev.band3 == 2).sum()),
          "n_prevalent_cardiopulmonary_any": int(sev[[f"{k}_prevalent" for k in CP]].fillna(0).max(axis=1).sum()) if CP else None,
          "median_pre_t90_by_band3": {k: round(float(sev.loc[sev.band3 == i, "pre_sleep_t90"].median()), 1)
                                      for i, k in enumerate(("ref_le1", "A_1_to_10", "B_gt10"))},
          "median_post_t90_by_band3": {k: round(float(sev.loc[sev.band3 == i, "post_sleep_t90"].median()), 1)
                                       for i, k in enumerate(("ref_le1", "A_1_to_10", "B_gt10"))},
          "median_age_by_band3": {k: round(float(sev.loc[sev.band3 == i, "AgeAtVisit"].median()), 1)
                                  for i, k in enumerate(("ref_le1", "A_1_to_10", "B_gt10"))},
          "pct_male_by_band3": {k: round(100 * float(sev.loc[sev.band3 == i, "male"].mean()), 1)
                                for i, k in enumerate(("ref_le1", "A_1_to_10", "B_gt10"))}}
T = json.load(open(f"{NUMBERS_DIR}/treatment_v2.json"))
CVN = T["corrected_vs_not"]
assert (COUNTS["n_all"], COUNTS["n_restored_le10"], COUNTS["n_still_gt10"]) == (CVN["n"], CVN["n_corrected"], CVN["n_not"]), \
    (COUNTS, {k: CVN[k] for k in ("n", "n_corrected", "n_not")})
assert sum(cnt4) == COUNTS["n_all"] and COUNTS["n_ref_le1"] == cnt4[0] and COUNTS["n_B_gt10"] == cnt4[3], (cnt4, COUNTS)
print("counts ok against treatment_v2.json:", COUNTS)

# ------------------------------------------------------------------ fits
PUB = CVN["outcomes"]
N_PUB = len(PUB)
KEYS = list(DISEASES.keys()) + ["death"]


def frame(k):
    f = sev[(sev[f"{k}_prevalent"] == 0) & sev[f"{k}_years"].notna() & (sev[f"{k}_years"] > 0)]
    cols = {"T": f[f"{k}_years"], "E": f[f"{k}_incident"].astype(int), "x": f.notcorr,
            "xA": f.xA, "xB": f.xB, "age": f.AgeAtVisit, "male": f.male, "lbase": f.lbase,
            "site": f.site_id, "band3": f.band3}
    for c in SPL:
        cols[c] = f[c]
    others = [f"{j}_prevalent" for j in CP if j != k]
    cols["cp_any"] = f[others].fillna(0).max(axis=1).astype(int) if others else 0
    return pd.DataFrame(cols).dropna()


def cox(ff, cols):
    c = CoxPHFitter(penalizer=0.0).fit(ff[["T", "E", "site"] + cols], "T", "E", strata=["site"])
    return c


def est(c, name):
    lo, hi = np.exp(c.confidence_intervals_.loc[name])
    return {"hr": float(np.exp(c.params_[name])), "lo": float(lo), "hi": float(hi),
            "p": float(c.summary.loc[name, "p"])}


def bh(pvals):
    """Benjamini-Hochberg q values, in the order given; NaN stays NaN."""
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


from scipy import stats as st  # noqa: E402

rows, GATE = [], []
for k in KEYS:
    if f"{k}_incident" not in sev.columns:
        continue
    lab = DISEASES[k][0] if k in DISEASES else "Death from any cause"
    neg = k in NEGATIVE_CONTROLS
    ff = frame(k)
    if ff.E.sum() < 40:
        assert lab not in PUB, lab
        continue
    # --- step 1: the published contrast, same code path
    c2 = cox(ff, ["x", "age", "male", "lbase"])
    e2 = est(c2, "x")
    pub = PUB[lab]
    rep = (e2["hr"], e2["lo"], e2["hi"])  # v8.3 full precision (2026-09-16)
    # v8.3: reference now full precision, compare within 5e-4 (2026-09-16)
    ok = all(abs(a - b) <= 5e-4 for a, b in zip(rep, (pub["hr"], pub["lo"], pub["hi"]))) and int(ff.E.sum()) == pub["events"] \
        and abs(e2["p"] - pub["p"]) < 1e-9
    GATE.append({"outcome": lab, "events": int(ff.E.sum()), "events_published": pub["events"],
                 "hr_refit": rep[0], "lo_refit": rep[1], "hi_refit": rep[2], "p_refit": e2["p"],
                 "hr_published": pub["hr"], "lo_published": pub["lo"], "hi_published": pub["hi"],
                 "p_published": pub["p"], "match_3dp": bool(ok)})
    # --- step 2: the three-level exposure against <=1% on PAP
    c3 = cox(ff, ["xA", "xB", "age", "male", "lbase"])
    eA, eB = est(c3, "xA"), est(c3, "xB")
    c0 = cox(ff, ["age", "male", "lbase"])
    lr = 2 * (c3.log_likelihood_ - c0.log_likelihood_)
    p_lr = float(st.chi2.sf(lr, 2))
    # --- sensitivity: cr(df=4) age spline, one column dropped, both contrasts
    s2 = est(cox(ff, ["x"] + SPL + ["male", "lbase"]), "x")
    cs = cox(ff, ["xA", "xB"] + SPL + ["male", "lbase"])
    sA, sB = est(cs, "xA"), est(cs, "xB")
    ev = {g: int(ff.loc[ff.band3 == i, "E"].sum()) for i, g in enumerate(("ref", "A", "B"))}
    nn = {g: int((ff.band3 == i).sum()) for i, g in enumerate(("ref", "A", "B"))}
    row = {"outcome": lab, "key": k, "negative_control": neg, "n_at_risk": int(len(ff)),
           "events": int(ff.E.sum()),
           "n_ref": nn["ref"], "n_A": nn["A"], "n_B": nn["B"],
           "events_ref": ev["ref"], "events_A": ev["A"], "events_B": ev["B"],
           "A_hr": eA["hr"], "A_lo": eA["lo"], "A_hi": eA["hi"], "A_p": eA["p"],
           "B_hr": eB["hr"], "B_lo": eB["lo"], "B_hi": eB["hi"], "B_p": eB["p"],
           "p_lr_2df": p_lr,
           "pub_hr": pub["hr"], "pub_lo": pub["lo"], "pub_hi": pub["hi"], "pub_p": pub["p"],
           "refit_hr": e2["hr"], "refit_lo": e2["lo"], "refit_hi": e2["hi"], "refit_p": e2["p"],
           "spline_pub_hr": s2["hr"], "spline_pub_lo": s2["lo"], "spline_pub_hi": s2["hi"],
           "spline_pub_p": s2["p"],
           "spline_A_hr": sA["hr"], "spline_A_lo": sA["lo"], "spline_A_hi": sA["hi"],
           "spline_A_p": sA["p"],
           "spline_B_hr": sB["hr"], "spline_B_lo": sB["lo"], "spline_B_hi": sB["hi"],
           "spline_B_p": sB["p"],
           "in_current_5b": (not neg) and pub["p"] < 0.05}
    if ADJUST:
        # decision 6: the same three contrasts with the prevalent cardiopulmonary covariate
        row["n_cp_any"] = int(ff.cp_any.sum())
        a2 = est(cox(ff, ["x", "age", "male", "lbase", "cp_any"]), "x")
        ca = cox(ff, ["xA", "xB", "age", "male", "lbase", "cp_any"])
        aA, aB = est(ca, "xA"), est(ca, "xB")
        ca0 = cox(ff, ["age", "male", "lbase", "cp_any"])
        row.update({"adj_pub_hr": a2["hr"], "adj_pub_lo": a2["lo"], "adj_pub_hi": a2["hi"], "adj_pub_p": a2["p"],
                    "adj_A_hr": aA["hr"], "adj_A_lo": aA["lo"], "adj_A_hi": aA["hi"], "adj_A_p": aA["p"],
                    "adj_B_hr": aB["hr"], "adj_B_lo": aB["lo"], "adj_B_hi": aB["hi"], "adj_B_p": aB["p"],
                    "adj_p_lr_2df": float(st.chi2.sf(2 * (ca.log_likelihood_ - ca0.log_likelihood_), 2))})
    rows.append(row)

df = pd.DataFrame(rows)
assert len(df) == N_PUB and set(df.outcome) == set(PUB), (len(df), set(PUB) ^ set(df.outcome))
n_bad = sum(not g["match_3dp"] for g in GATE)
print(f"positive control: {len(GATE)} outcomes, {len(GATE) - n_bad} reproduce to 3 dp, {n_bad} differ")
for g in GATE:
    if not g["match_3dp"]:
        print("   MISMATCH", g)
assert n_bad == 0, "the published 5b numbers did not reproduce; stop"
n_neg_expected = sum(1 for k in NEGATIVE_CONTROLS if DISEASES[k][0] in PUB)
n_5b_expected = sum(1 for lab, v in PUB.items() if v["p"] < 0.05 and not v.get("negative_control", False))
assert df.negative_control.sum() == n_neg_expected and df.in_current_5b.sum() == n_5b_expected, \
    (int(df.negative_control.sum()), n_neg_expected, int(df.in_current_5b.sum()), n_5b_expected)

if BH:
    # decision 7: one q per band across the non-control conditions; the controls across themselves
    pcols = ["pub_p", "A_p", "B_p"] + (["adj_pub_p", "adj_A_p", "adj_B_p"] if ADJUST else [])
    for pc in pcols:
        qc = pc[:-2] + "_q"
        df[qc] = np.nan
        for flag in (False, True):
            idx = df.index[df.negative_control == flag]
            df.loc[idx, qc] = bh(df.loc[idx, pc].values)

df = df.sort_values("B_hr", ascending=False).reset_index(drop=True)
csv_path = f"{ARGS.out_dir}/L4_fig5b_refit{SFX}.csv"
json_path = f"{ARGS.out_dir}/L4_fig5b_refit{SFX}.json"
df.to_csv(csv_path, index=False)
built = "2026-09-03 L4_refit_fig5b.py" if not (ADJUST or BH) else \
    f"{pd.Timestamp.now():%Y-%m-%d} L4_refit_fig5b.py" + (" --adjust-prevalent-cardiopulmonary" if ADJUST else "") + (" --bh" if BH else "")
out = {"built": built, "counts": COUNTS,
       "model": {"exposure": "on-PAP sleep T90 band: ref <=1%, A (1,10], B >10",
                 "adjustment_primary": "log1p pretreatment sleep T90, age (linear, as the "
                                       "published run_treatment_v2.py), sex, site strata",
                 "adjustment_sensitivity": "cr(age, df=4) with the first column dropped in "
                                           "place of linear age",
                 "adjustment_prevalent_cardiopulmonary": (f"any prevalent condition among {CP} (the outcome "
                                                          "under study excluded), added to the primary adjustment") if ADJUST else None,
                 "multiple_testing": ("Benjamini-Hochberg per band across the non-control conditions; "
                                      "the negative controls corrected among themselves") if BH else None,
                 "fit": "lifelines CoxPHFitter penalizer 0, strata site_id",
                 "event_floor": 40, "risk_set": "free of the outcome at the sleep study, "
                                                "follow-up years > 0"},
       "positive_control": GATE,
       "outcomes": {r["outcome"]: {kk: (bool(v) if isinstance(v, (bool, np.bool_)) else v)
                                   for kk, v in r.items() if kk != "outcome"}
                    for r in df.to_dict("records")}}
json.dump(out, open(json_path, "w"), indent=1, default=float)
for p in (csv_path, json_path):
    sidecar(p, __file__, extra_inputs=[f"{NUMBERS_DIR}/treatment_v2.json"],
            note=f"L4_refit_fig5b adjust={ADJUST} bh={BH}, {len(df)} outcomes")

sig17 = df[df.in_current_5b].sort_values("pub_hr", ascending=False)
print(f"\n{'outcome':<28}{'ev':>5}  {'A: 1-10% vs <=1%':>24}  {'B: >10% vs <=1%':>24}   {'pub >10 vs <=10':>20}")
for _, r in pd.concat([sig17, df[df.negative_control]]).iterrows():
    fa = f"{r.A_hr:.2f} ({r.A_lo:.2f}-{r.A_hi:.2f}) p={r.A_p:.3g}"
    fb = f"{r.B_hr:.2f} ({r.B_lo:.2f}-{r.B_hi:.2f}) p={r.B_p:.3g}"
    fp = f"{r.pub_hr:.2f} ({r.pub_lo:.2f}-{r.pub_hi:.2f})"
    print(f"{r.outcome:<28}{r.events:>5}  {fa:>24}  {fb:>24}   {fp:>20}{'  NEG' if r.negative_control else ''}")
nc = ~df.negative_control
print("\nB significant (P<0.05), non-control:", int(((df.B_p < 0.05) & nc).sum()),
      " A significant:", int(((df.A_p < 0.05) & nc).sum()),
      f" B HR>1 of {int(nc.sum())}:", int(((df.B_hr > 1) & nc).sum()))
if BH:
    print("q < 0.05, non-control: pub", int(((df.pub_q < 0.05) & nc).sum()), " A", int(((df.A_q < 0.05) & nc).sum()),
          " B", int(((df.B_q < 0.05) & nc).sum()),
          (f"  adjusted: pub {int(((df.adj_pub_q < 0.05) & nc).sum())} A {int(((df.adj_A_q < 0.05) & nc).sum())} "
           f"B {int(((df.adj_B_q < 0.05) & nc).sum())}") if ADJUST else "")
print("written:", csv_path, json_path)
