"""
The 2-year-landmark arm for eFigure 4: the JOINT per-SD T90 + total-sleep-time model,
refit at lag 0 (all follow-up) and at the 2-year landmark, across the 48 ranked outcomes.

Ordered in the 2026-08-21 figure reorganization round (supplement pack A). The point the
panel has to carry: in one model holding both exposures, does per-SD T90 keep its
associations when the first two years of diagnoses and follow-up are deleted, and does
total sleep time stay flat.

MODEL, the paper's primary spec (ROUND_BRIEF section 0, identical to
New_Figures/_NOT_THESE_workfiles/scripts/lagladder.py and numbers/run_primary_unpenalized.py):
site-stratified Cox (Efron), unpenalized, age as a natural cubic spline with 4 df with the
first column dropped (rank trap), plus sex. Both exposures rank-inverse-normal transformed
within site, so each hazard ratio is per 1 SD. PSG date is time zero. oximetry_bad == 0.

LANDMARK, the published ladder's convention (numbers/lag_ladder.json provenance): for lag L,
anyone whose event or censoring arrived at or before L is removed, and the clock restarts
at L for the survivors.

POSITIVE CONTROLS, run before any joint fit and aborting on mismatch:
  V1  single-exposure T90 at lag 0 for heart failure, COPD and type 2 diabetes must
      reproduce numbers/bdsp_diseases_v3.csv (hr, lo, hi to 3 dp, n and events exactly).
  V2  the same three outcomes at the 2-year landmark must reproduce the published
      lag ladder rungs in numbers/lag_ladder.csv (hr, lo, hi to 3 dp, n and events exactly).
V1 proves cohort, eligibility, spline, transform and estimator. V2 proves the landmark
windowing. Only then is the joint model fitted.

FDR: BH within exposure within arm across the 48 ranked outcomes, the same convention the
frozen 197 x 48 matrix behind main Figure 2c uses (its q is BH within measurement across
the same 48, three of which are settled negative controls kept as ordinary rows there).

Outputs, all in this folder:
  joint_t90_tst_landmark.json   per-outcome estimates, both exposures, both arms, with q
  joint_t90_tst_landmark.csv    the same, long format
  validation.json               the V1/V2 reproduction, value by value
  provenance.json               inputs, spec, counts, date
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import hashlib
import json
import os
import sys
import time
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from patsy import dmatrix
from scipy import stats

T90ROOT = paths.T90_ROOT
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
from cohort_spec import COHORT_N, apply_cohort, N_MEASURES, N_RANKED_OUTCOMES   # noqa: E402  (v8: counts from the spec)
from disease_definitions import DISEASES                  # noqa: E402

LAGS = [0.0, 2.0]
MIN_EVENTS = 40                    # the ladder's fit floor
VAL_KEYS = ["hf", "copd2", "diabetes"]
MATRIX = (f"{paths.SV_ROOT}/"
          "dC_side_analyses/hr_matrix_141x52.csv")   # v8 matrix (141 measures x 52 ranked outcomes)

t0 = time.time()

# ------------------------------------------------------------------ the 48 ranked outcomes
M = pd.read_csv(MATRIX)
assert len(M) == N_MEASURES * N_RANKED_OUTCOMES and M.measurement.nunique() == N_MEASURES, (len(M), M.measurement.nunique())
OUTCOMES = sorted(M[M.measurement == "spo2_pct_below_90"].outcome)
assert len(OUTCOMES) == N_RANKED_OUTCOMES, len(OUTCOMES)

# ------------------------------------------------------------------ cohort, primary spec
b = apply_cohort(pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet"))
assert len(b) == COHORT_N
b["male"] = (b.sex.astype(str).str.upper().str[0] == "M").astype(int)

sp = dmatrix("cr(a, df=4) - 1", {"a": b.AgeAtVisit.values},
             return_type="dataframe").iloc[:, 1:]          # drop one column, rank trap
for i in range(sp.shape[1]):
    b[f"age_s{i}"] = sp.iloc[:, i].values
ADJ = [f"age_s{i}" for i in range(sp.shape[1])] + ["male"]


def rint(x):
    """Rank-inverse-normal, NaN-safe: missing values stay missing and drop at fit time."""
    x = np.asarray(x, dtype=float)
    out = np.full(len(x), np.nan)
    m = np.isfinite(x)
    r = stats.rankdata(x[m])
    out[m] = stats.norm.ppf((r - 0.375) / (m.sum() + 0.25))
    return out


b["z_t90"] = b.groupby("site_id").spo2_pct_below_90.transform(rint)
b["z_tst"] = b.groupby("site_id").TST_min.transform(rint)
n_tst_na = int(b.z_tst.isna().sum())
assert b.z_t90.notna().all()
print(f"cohort {len(b):,}, TST missing in {n_tst_na} recordings")


def base_frame(key, cols):
    """Eligibility, identical to the ladder: prevalent out, positive follow-up required."""
    yc, ec, pc = f"{key}_years", f"{key}_incident", f"{key}_prevalent"
    f = b[(b[pc] == 0) & b[yc].notna() & (b[yc] > 0)]
    return pd.DataFrame({"T": f[yc].astype(float), "E": f[ec].astype(int),
                         "site": f.site_id, **{c: f[c] for c in cols + ADJ}})


def landmark(base, L):
    """The ladder's convention: event or censoring at or before L removes the patient."""
    lm = base[base["T"] > L].copy()
    lm["T"] = lm["T"] - L
    return lm


def fit(d, exposures):
    d = d[["T", "E", "site"] + exposures + ADJ].dropna()
    n, e = int(len(d)), int(d["E"].sum())
    if e < MIN_EVENTS:
        return {"n": n, "events": e, "estimable": False}
    c = CoxPHFitter().fit(d, "T", "E", strata=["site"])
    out = {"n": n, "events": e, "estimable": True}
    for x in exposures:
        s = c.summary.loc[x]
        out[x] = {"hr": float(s["exp(coef)"]), "lo": float(s["exp(coef) lower 95%"]),
                  "hi": float(s["exp(coef) upper 95%"]), "p": float(s["p"])}
    return out


# ------------------------------------------------------------------ V1 and V2, abort on miss
REF0 = pd.read_csv(f"{paths.NUMBERS_DIR}/bdsp_diseases_v3.csv").set_index("key")
LAD = pd.read_csv(f"{paths.NUMBERS_DIR}/lag_ladder.csv")
LAD2 = LAD[LAD.lag_years == 2.0].set_index("key")
LAD0 = LAD[LAD.lag_years == 0.0].set_index("key")

validation = {"V1_lag0_vs_bdsp_diseases_v3": {}, "V2_lag2_vs_lag_ladder": {}}
for key in VAL_KEYS:
    base = base_frame(key, ["z_t90"])

    r = fit(landmark(base, 0.0), ["z_t90"])
    v = REF0.loc[key]
    rec = {"mine": {"n": r["n"], "events": r["events"],
                    "hr": round(r["z_t90"]["hr"], 3), "lo": round(r["z_t90"]["lo"], 3),
                    "hi": round(r["z_t90"]["hi"], 3)},
           "published": {"n": int(v.n), "events": int(v.events),
                         "hr": round(float(v.t90_hr), 3), "lo": round(float(v.t90_lo), 3),
                         "hi": round(float(v.t90_hi), 3)}}
    rec["match"] = rec["mine"] == rec["published"]
    validation["V1_lag0_vs_bdsp_diseases_v3"][key] = rec
    assert rec["match"], f"V1 FAILED for {key}: {rec}"
    # the ladder's own lag-0 row is the same fit, cross-checked while we are here
    lv = LAD0.loc[key]
    assert (int(lv.n), int(lv.events)) == (r["n"], r["events"]), (key, "ladder lag0 n/events")
    assert round(float(lv.hr), 3) == rec["mine"]["hr"], (key, "ladder lag0 hr")

    r2 = fit(landmark(base, 2.0), ["z_t90"])
    l2 = LAD2.loc[key]
    rec2 = {"mine": {"n": r2["n"], "events": r2["events"],
                     "hr": round(r2["z_t90"]["hr"], 3), "lo": round(r2["z_t90"]["lo"], 3),
                     "hi": round(r2["z_t90"]["hi"], 3)},
            "published": {"n": int(l2.n), "events": int(l2.events),
                          "hr": round(float(l2.hr), 3), "lo": round(float(l2.lo), 3),
                          "hi": round(float(l2.hi), 3)}}
    rec2["match"] = rec2["mine"] == rec2["published"]
    validation["V2_lag2_vs_lag_ladder"][key] = rec2
    assert rec2["match"], f"V2 FAILED for {key}: {rec2}"
    print(f"validated {key}: lag0 {rec['mine']['hr']} and lag2 {rec2['mine']['hr']} "
          "reproduce the published values exactly")

json.dump(validation, open(f"{HERE}/validation.json", "w"), indent=1)
print("V1 and V2 PASS, landmark machinery reproduces the published rungs\n")

# ------------------------------------------------------------------ the joint model, both arms
rows = []
for key in OUTCOMES:
    base = base_frame(key, ["z_t90", "z_tst"])
    for L in LAGS:
        r = fit(landmark(base, L), ["z_t90", "z_tst"])
        rec = {"key": key, "lag_years": L, "n": r["n"], "events": r["events"],
               "estimable": r["estimable"]}
        if r["estimable"]:
            for x, tag in (("z_t90", "t90"), ("z_tst", "tst")):
                for f_ in ("hr", "lo", "hi", "p"):
                    rec[f"{tag}_{f_}"] = r[x][f_]
        rows.append(rec)
    print(f"fitted {key}", flush=True)

d = pd.DataFrame(rows)
assert d.estimable.all(), d[~d.estimable]

# BH within exposure within arm across the 48 ranked outcomes
def bh(p):
    p = np.asarray(p, float)
    m = len(p)
    order = np.argsort(p)
    q = np.empty(m)
    prev = 1.0
    for rank_from_top in range(m, 0, -1):
        i = order[rank_from_top - 1]
        prev = min(prev, p[i] * m / rank_from_top)
        q[i] = prev
    return q


for tag in ("t90", "tst"):
    for L in LAGS:
        m = d.lag_years == L
        d.loc[m, f"{tag}_q"] = bh(d.loc[m, f"{tag}_p"].values)

d = d.sort_values(["key", "lag_years"])
d.to_csv(f"{HERE}/joint_t90_tst_landmark.csv", index=False)

# ------------------------------------------------------------------ summary and provenance
LABEL = {k: v[0] for k, v in DISEASES.items()}
LABEL["death"] = "Death from any cause"

summary = {}
for tag, nice in (("t90", "T90 per SD"), ("tst", "Total sleep time per SD")):
    for L in LAGS:
        s = d[d.lag_years == L]
        raised = s[(s[f"{tag}_q"] < 0.05) & (s[f"{tag}_hr"] > 1)]
        summary[f"{tag}_lag{int(L)}"] = {
            "exposure": nice, "arm": "all follow-up" if L == 0 else "2-year landmark",
            "n_fdr_raised": int(len(raised)), "of": len(OUTCOMES),
            "n_fdr_significant_either_direction": int((s[f"{tag}_q"] < 0.05).sum()),
            "median_hr": round(float(s[f"{tag}_hr"].median()), 3),
            "fdr_raised_outcomes": sorted(LABEL[k] for k in raised.key),
        }

t90_both = set(d[(d.lag_years == 0) & (d.t90_q < .05) & (d.t90_hr > 1)].key) \
    & set(d[(d.lag_years == 2) & (d.t90_q < .05) & (d.t90_hr > 1)].key)
summary["t90_fdr_raised_at_both_arms"] = {"n": len(t90_both),
                                          "outcomes": sorted(LABEL[k] for k in t90_both)}
s2 = d[d.lag_years == 2.0].set_index("key")
s0 = d[d.lag_years == 0.0].set_index("key")
rem = 100.0 * (1.0 - s2.events / s0.events)
summary["landmark_cost"] = {
    "median_pct_events_removed": round(float(rem.median()), 1),
    "min_pct_events_removed": round(float(rem.min()), 1),
    "max_pct_events_removed": round(float(rem.max()), 1),
}


def digest(p):
    h = hashlib.sha256(open(p, "rb").read()).hexdigest()[:16]
    return {"path": p, "sha256_16": h, "bytes": os.path.getsize(p)}


provenance = {
    "built_by": "reorg_fits_efig4_lag2/fit_joint_t90_tst_landmark.py",
    "built_on": time.strftime("%Y-%m-%d %H:%M"),
    "for": "eFigure 4 landmark panel, figure reorganization round 2026-08-21",
    "model": ("unpenalized site-stratified Cox (Efron), 4-df natural cubic age spline with "
              "the first column dropped, sex, z_t90 and z_tst together in one model, both "
              "rank-inverse-normalized within site, hazard ratios per 1 SD"),
    "landmark": ("for lag L, anyone whose event or censoring arrived at or before L is "
                 "removed and the clock restarts at L (numbers/lag_ladder.json convention)"),
    "fdr": (f"BH within exposure within arm across the {N_RANKED_OUTCOMES} ranked outcomes, the frozen "
            f"hr_matrix_{N_MEASURES}x{N_RANKED_OUTCOMES} convention"),
    "outcomes": OUTCOMES,
    "n_cohort": int(COHORT_N),
    "tst_missing_n": n_tst_na,
    "min_events_to_fit": MIN_EVENTS,
    "validated_against": [digest(f"{paths.NUMBERS_DIR}/bdsp_diseases_v3.csv"),
                          digest(f"{paths.NUMBERS_DIR}/lag_ladder.csv")],
    "inputs": [digest(f"{paths.TABLES_DIR}/t90_final.parquet"), digest(MATRIX)],
    "validation_keys": VAL_KEYS,
    "runtime_min": round((time.time() - t0) / 60.0, 1),
}
json.dump({"summary": summary, "outcomes": {
    k: {str(int(L)): {c: (float(r[c]) if isinstance(r[c], (int, np.floating, float))
                          else r[c])
                      for c in ("n", "events", "t90_hr", "t90_lo", "t90_hi", "t90_p",
                                "t90_q", "tst_hr", "tst_lo", "tst_hi", "tst_p", "tst_q")}
               for L, r in ((L, d[(d.key == k) & (d.lag_years == L)].iloc[0])
                            for L in LAGS)}
    for k in OUTCOMES}, "labels": LABEL, "provenance": provenance},
    open(f"{HERE}/joint_t90_tst_landmark.json", "w"), indent=1, default=float)
json.dump(provenance, open(f"{HERE}/provenance.json", "w"), indent=1)

print("\n" + "=" * 74)
for k, v in summary.items():
    if isinstance(v, dict) and "exposure" in v:
        print(f"{v['exposure']:<26} {v['arm']:<16} {v['n_fdr_raised']:>2} of {v['of']} "
              f"FDR-raised, median HR {v['median_hr']}")
print(f"T90 FDR-raised at BOTH arms: {summary['t90_fdr_raised_at_both_arms']['n']}")
print(f"2-year landmark removes a median "
      f"{summary['landmark_cost']['median_pct_events_removed']}% of events")
print(f"\nwrote {HERE}/joint_t90_tst_landmark.json  (runtime "
      f"{provenance['runtime_min']} min)")
