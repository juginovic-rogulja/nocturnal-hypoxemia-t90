"""
Is time below 90% saturation better expressed as MINUTES or as a PERCENTAGE of the night?

A percentage confounds severity with sleep duration: 10% of a 4-hour night is 24 minutes, 10%
of an 8-hour night is 48. If the biology responds to absolute hypoxic time, minutes should
discriminate better and the clinical threshold should be stated in minutes.

Method matches the paper's ranking exactly: each measure is added to a model containing age and
sex, fitted in 1 hospital and evaluated in the other, in both directions, and scored by the
average gain in held-out concordance across outcomes.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../..")))  # repository root, where paths.py lives
import paths
import warnings; warnings.filterwarnings("ignore")
import json, sys
import numpy as np, pandas as pd
from scipy import stats
from patsy import dmatrix
from lifelines import CoxPHFitter

ROOT = paths.T90_ROOT
sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
from disease_definitions import DISEASES, CIRCULAR
from cohort_spec import drop_split_nights_if_full   # v8.1

b = pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet")
b = b[(b.fu_valid == 1) & b.spo2_pct_below_90.notna() & (b.oximetry_bad == 0)].copy()
b = drop_split_nights_if_full(b)   # v8.1: recording minutes are a sleep amount, full diagnostic nights only
# drop the impossible-oximetry block found in the phenotype audit
bad = b.spo2_nadir_corrected > b.spo2_mean
print(f"cohort {len(b):,}, excluding {bad.sum()} artifact recordings -> {(~bad).sum():,}")
b = b[~bad].copy()

b["male"] = (b.sex.astype(str).str.upper().str[0] == "M").astype(int)
sp = dmatrix("cr(a, df=4) - 1", {"a": b.AgeAtVisit.values}, return_type="dataframe")
sp.columns = [f"age_s{i}" for i in range(4)]; sp.index = b.index
for c in sp.columns:
    b[c] = sp[c]
ADJ = [f"age_s{i}" for i in range(4)] + ["male"]

# the two candidate scales
b["t90_pct"] = b.spo2_pct_below_90
b["t90_min"] = b.spo2_pct_below_90 / 100.0 * b.recording_dur_min
b["t90_min_tst"] = b.spo2_pct_below_90 / 100.0 * b.TST_min
print(f"median minutes below 90%: {b.t90_min.median():.1f} over the recording, "
      f"{b.t90_min_tst.median():.1f} over sleep")
print(f"correlation between the percentage and the minutes: "
      f"{b.t90_pct.corr(b.t90_min, method='spearman'):.3f}")
print(f"sleep duration varies: TST median {b.TST_min.median():.0f} min, "
      f"IQR {b.TST_min.quantile(.25):.0f}-{b.TST_min.quantile(.75):.0f}\n")

def rint(x):
    r = stats.rankdata(x)
    return stats.norm.ppf((r - 0.375) / (len(r) + 0.25))

CAND = ["t90_pct", "t90_min", "t90_min_tst"]
for c in CAND:
    b[c + "__z"] = b.groupby("site_id")[c].transform(rint)

NAMES = {k: v[0] for k, v in DISEASES.items()}
OUT = [k for k in DISEASES if k not in CIRCULAR] + ["death"]

def heldout(cols):
    gains = []
    for tr, te in [("I0002", "I0006"), ("I0006", "I0002")]:
        for k in OUT:
            yc, ec = f"{k}_years", f"{k}_incident"
            if yc not in b.columns:
                continue
            g = b[(b[f"{k}_prevalent"] == 0) & b[yc].notna() & (b[yc] > 0)]
            A = g[g.site_id == tr][ADJ + cols + [yc, ec]].dropna()
            B = g[g.site_id == te][ADJ + cols + [yc, ec]].dropna()
            if A[ec].sum() < 40 or B[ec].sum() < 25:
                continue
            try:
                m0 = CoxPHFitter(penalizer=0.01).fit(A[ADJ + [yc, ec]], yc, ec)
                m1 = CoxPHFitter(penalizer=0.01).fit(A[ADJ + cols + [yc, ec]], yc, ec)
                from lifelines.utils import concordance_index as ci
                c0 = ci(B[yc], -m0.predict_partial_hazard(B), B[ec])
                c1 = ci(B[yc], -m1.predict_partial_hazard(B), B[ec])
            except Exception:
                continue
            gains.append(c1 - c0)
    return float(np.mean(gains)), len(gains)

R = {"n": int(len(b)), "scales": {}}
print(f"{'scale':34s}{'gain':>10}{'n outcomes':>12}")
print("-" * 58)
for c in CAND:
    g, k = heldout([c + "__z"])
    R["scales"][c] = {"gain": round(g, 5), "n_models": k}
    lab = {"t90_pct": "percentage of the recording",
           "t90_min": "minutes below 90% (recording)",
           "t90_min_tst": "minutes below 90% (sleep)"}[c]
    print(f"{lab:34s}{g:10.5f}{k:12d}")

# does adding sleep duration to the percentage help? if minutes were the right scale it should
g2, _ = heldout(["t90_pct__z", "TST_min"])
R["pct_plus_tst_gain"] = round(g2, 5)
print(f"{'percentage + total sleep time':34s}{g2:10.5f}")

# clinical thresholds on each scale, against death
print(f"\n{'threshold':34s}{'n above':>9}{'% cohort':>10}{'HR death':>10}{'95% CI':>16}")
print("-" * 80)
TH = ([("percentage", "t90_pct", v) for v in (5, 10, 20, 30)] +
      [("minutes", "t90_min", v) for v in (10, 20, 30, 45, 60, 90)])
R["thresholds"] = []
for lab, col, v in TH:
    b["_x"] = (b[col] > v).astype(int)
    g = b[(b.death_prevalent == 0) & b.death_years.notna() & (b.death_years > 0)]
    d = pd.DataFrame({"T": g.death_years, "E": g.death_incident.astype(int),
                      "x": g["_x"], "site": g.site_id,
                      **{c: g[c] for c in ADJ}}).dropna()
    try:
        m = CoxPHFitter(penalizer=0.01).fit(d, "T", "E", strata=["site"])
        s = m.summary.loc["x"]
        hr, lo, hi = s["exp(coef)"], s["exp(coef) lower 95%"], s["exp(coef) upper 95%"]
    except Exception:
        continue
    n = int(b["_x"].sum())
    R["thresholds"].append({"scale": lab, "cut": v, "n": n,
                            "pct": round(n / len(b) * 100, 1), "hr": round(float(hr), 3),
                            "lo": round(float(lo), 3), "hi": round(float(hi), 3)})
    print(f"{lab+' > '+str(v):34s}{n:9,}{n/len(b)*100:9.1f}%{hr:10.2f}   {lo:.2f}-{hi:.2f}")

json.dump(R, open(f"{paths.NUMBERS_DIR}/minutes_vs_pct_v2.json", "w"), indent=1)
print("\nwritten -> minutes_vs_pct_v2.json")
