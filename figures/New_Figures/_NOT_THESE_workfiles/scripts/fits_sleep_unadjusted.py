"""
The awake-oxygen adjustment test, fitted properly and with the comparator it was missing.

Two problems with what was on file.

  1. NO COMPARATOR. numbers/wake_sleep_healthy_A_full.csv holds the ADJUSTED estimate, sleep
     T90 with awake T90 in the same model. It does not hold sleep T90 on its own, so there is
     nothing to say how far the estimate moved on adjustment. A sentence like "diabetes goes
     from 1.41 to 1.31" needs a 1.41 that exists.

  2. THE FROZEN FIT IS RIDGE-PENALIZED. numbers/run_wake_sleep_healthy.py line 152 fits with
     penalizer=0.01. lifelines scales its penalty by the sample size, so on 6,803 patients that
     is an effective ridge near 70, and every hazard ratio in wake_sleep_healthy_A_full.csv is
     shrunk toward 1 by an amount that depends on the event count. The same script also carries
     all 4 age-spline columns, which is only survivable because of that penalty: unpenalized,
     the design is singular and the fit fails outright.

So both specifications are fitted here.

  penalized     penalizer=0.01, all 4 spline columns. Reproduces the frozen file exactly, which
                is the proof that the sample and conventions are right. Not the paper's
                primary specification and not what should be drawn.
  unpenalized   penalizer=0, one spline column dropped. The paper's primary specification, and
                what the redrawn eFigure 7 panel C and eFigure 11 use.

The sample is the CURRENT 19,173 cohort, entered through cohort_spec.apply_cohort like every
other analysis in the paper. It was the pre-ICD-fix cohort until 2026-08-07, which was a
standing-rule violation: that file predates the ICD-9 278.01/278.03 correction, the E888/E884
prefix-collision fix and the oximetry quality rule, and it has no column for contact dermatitis
or hemorrhoids, so two of the five settled negative controls could never be fitted here. The
merged wake-and-sleep sample is 6,800 on the current cohort, against 6,803 on the old one.

Writes numbers/sleep_unadjusted_v1.json.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../..")))  # repository root, where paths.py lives
import paths
import glob
import json
import sys
import warnings

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from patsy import dmatrix
from scipy import stats

warnings.filterwarnings("ignore")
ROOT = paths.FIGURE_ROOT
SHARDS = f"{paths.V8_ROOT}/X4_groupP/cpap_stage/out*.csv"   # v8.2: re-extracted, collapse-masked shards (X4_groupP)
sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
from disease_definitions import DISEASES, NEGATIVE_CONTROLS

MIN_EVENTS, MIN_WAKE_MIN = 60, 10.0

files = sorted(glob.glob(SHARDS))
assert files, f"no cpap_stage shards under {SHARDS}"   # v8.2: the shard count is a property of the fleet, not a result
raw = pd.concat([pd.read_csv(f, low_memory=False) for f in files], ignore_index=True)
st = raw.drop_duplicates("BDSPPatientID")
st = st[st.status == "ok"]

from cohort_spec import apply_cohort  # noqa: E402
b = apply_cohort(pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet"))
n_impossible = int((b.spo2_nadir_corrected > b.spo2_mean).sum())
b = b[~(b.spo2_nadir_corrected > b.spo2_mean)]
m = b.merge(st[["BDSPPatientID", "all_wake_t90", "all_sleep_t90", "all_wake_min"]],
            on="BDSPPatientID", how="inner")
m = m[m.all_wake_min >= MIN_WAKE_MIN].dropna(
    subset=["all_wake_t90", "all_sleep_t90"]).reset_index(drop=True)
m["male"] = (m.sex.astype(str).str.upper().str[0] == "M").astype(int)
sp = dmatrix("cr(a, df=4) - 1", {"a": m.AgeAtVisit.values}, return_type="dataframe")
sp.columns = [f"age_s{i}" for i in range(4)]
sp.index = m.index
for c in sp.columns:
    m[c] = sp[c]
ADJ4 = [f"age_s{i}" for i in range(4)] + ["male"]      # frozen: 4 columns, only safe penalized
ADJ3 = [f"age_s{i}" for i in range(1, 4)] + ["male"]   # primary: one column dropped
print(f"merged wake and sleep sample n = {len(m):,}, "
      f"{n_impossible} impossible-oximetry recordings dropped")

REF = pd.read_csv(f"{paths.NUMBERS_DIR}/wake_sleep_healthy_A_full.csv")
FROZ = json.load(open(f"{paths.NUMBERS_DIR}/wake_sleep_healthy.json"))["strata"]["A_full"]
assert len(m) == FROZ["n"], (len(m), FROZ["n"])


def rint(x):
    x = np.asarray(x, float)
    r = stats.rankdata(x)
    return stats.norm.ppf((r - 0.375) / (len(r) + 0.25))


m["wz"] = rint(m.all_wake_t90.values)
m["sz"] = rint(m.all_sleep_t90.values)
CONDS = list(DISEASES.items()) + [("death", ("Death from any cause", False, [], []))]


def fit(g, key, cols, pen):
    adj = ADJ4 if pen else ADJ3
    f = pd.DataFrame({"T": g[f"{key}_years"], "E": g[f"{key}_incident"].astype(int),
                      "site_id": g.site_id, **{c: g[c] for c in adj},
                      **{c: g[c] for c in cols}}).dropna()
    f = f[f["T"] > 0]
    if f.E.sum() < MIN_EVENTS:
        return None
    try:
        c = CoxPHFitter(penalizer=pen).fit(f, "T", "E", strata=["site_id"])
    except Exception:
        return None
    out = {"events": int(f.E.sum()), "n": int(len(f))}
    for x in cols:
        lo, hi = np.exp(c.confidence_intervals_.loc[x])
        out[x] = {"hr": round(float(np.exp(c.params_[x])), 4), "lo": round(float(lo), 4),
                  "hi": round(float(hi), 4), "p": float(c.summary.loc[x, "p"])}
    return out


# stratum C, free of organ and metabolic disease, the right panel of eFigure 11
CARDIOPULM = ["copd2", "asthma", "resp_failure", "obesity_hypovent", "pulm_htn",
              "hf", "ihd", "mi", "afib", "stroke_any"]
RENAL_HEP_CA = ["ckd", "cirrhosis", "cancer_any"]
METABOLIC = ["diabetes", "obesity", "nafld"]
mask = np.ones(len(m), bool)
for k in CARDIOPULM + RENAL_HEP_CA + METABOLIC:
    mask &= (m[f"{k}_prevalent"] == 0).values
mC = m[mask].reset_index(drop=True).copy()
mC["wz"] = rint(mC.all_wake_t90.values)
mC["sz"] = rint(mC.all_sleep_t90.values)
REFC = pd.read_csv(f"{paths.NUMBERS_DIR}/wake_sleep_healthy_C_organ_metab_free.csv")
FROZC = json.load(open(f"{paths.NUMBERS_DIR}/wake_sleep_healthy.json"))["strata"]["C_organ_metab_free"]
assert len(mC) == FROZC["n"], (len(mC), FROZC["n"])
print(f"organ and metabolic disease free stratum n = {len(mC):,}")

rows, bad, skipped = {}, [], []
for key, (label, _rawneg, _a, _c) in CONDS:
    neg = key in NEGATIVE_CONTROLS
    if f"{key}_incident" not in m.columns:
        continue
    g = m[m[f"{key}_prevalent"] == 0]
    pen_adj = fit(g, key, ["sz", "wz"], 0.01)
    if pen_adj is None:
        skipped.append(label)
        continue
    ref = REF[REF.cond == label]
    if len(ref):
        r = ref.iloc[0]
        if abs(pen_adj["sz"]["hr"] - r.sHR) > 1e-3 or abs(pen_adj["wz"]["hr"] - r.wHR) > 1e-3:
            bad.append((label, float(r.sHR), pen_adj["sz"]["hr"]))
    unp_adj = fit(g, key, ["sz", "wz"], 0.0)
    unp_una = fit(g, key, ["sz"], 0.0)
    if unp_adj is None or unp_una is None:
        skipped.append(label)
        continue
    rows[label] = {
        "events": unp_adj["events"], "at_risk": unp_adj["n"], "negative_control": bool(neg),
        "penalized_frozen": {"sleep_adjusted": pen_adj["sz"], "awake_adjusted": pen_adj["wz"]},
        "sleep_alone": unp_una["sz"],
        "sleep_adjusted_for_awake": unp_adj["sz"],
        "awake_adjusted_for_sleep": unp_adj["wz"],
        "shift": round(unp_adj["sz"]["hr"] - unp_una["sz"]["hr"], 4),
        "excess_retained": (round((unp_adj["sz"]["hr"] - 1) / (unp_una["sz"]["hr"] - 1), 3)
                            if unp_una["sz"]["hr"] > 1.0005 else None),
        "penalty_shrinkage": round(unp_adj["sz"]["hr"] - pen_adj["sz"]["hr"], 4),
    }
    gC = mC[mC[f"{key}_prevalent"] == 0]
    penC = fit(gC, key, ["sz", "wz"], 0.01)
    unpC = fit(gC, key, ["sz", "wz"], 0.0)
    if penC is not None:
        refC = REFC[REFC.cond == label]
        if len(refC) and abs(penC["sz"]["hr"] - float(refC.iloc[0].sHR)) > 1e-3:
            bad.append((label + " [C]", float(refC.iloc[0].sHR), penC["sz"]["hr"]))
    if unpC is not None:
        rows[label]["organ_free"] = {
            "events": unpC["events"], "at_risk": unpC["n"],
            "sleep_adjusted_for_awake": unpC["sz"], "awake_adjusted_for_sleep": unpC["wz"],
            "penalized_frozen_sleep": penC["sz"] if penC else None}

print(f"validation, the penalized refit against wake_sleep_healthy_A_full.csv: "
      f"{len(bad)} mismatches over {len(rows)} conditions")
assert not bad, bad[:5]

sig = [k for k, v in rows.items() if v["sleep_alone"]["lo"] > 1]
surv = [k for k in sig if rows[k]["sleep_adjusted_for_awake"]["lo"] > 1]
ret = [rows[k]["excess_retained"] for k in sig if rows[k]["excess_retained"] is not None]
shr = [v["penalty_shrinkage"] for v in rows.values()]
print(f"conditions significant on sleep oxygen alone: {len(sig)}")
print(f"   still significant after adjustment for awake oxygen: {len(surv)}")
print(f"   lost: {[k for k in sig if k not in surv]}")
print(f"median excess hazard retained on adjustment: {np.median(ret):.3f}")
print(f"median shrinkage the frozen ridge penalty imposes on the sleep estimate: "
      f"{np.median(shr):+.4f}, largest {max(shr, key=abs):+.4f}")

print(f"\n{'condition':<28}{'sleep alone':>13}{'adj for awake':>15}{'retained':>10}{'frozen':>9}")
for k in sorted(surv, key=lambda k: -rows[k]["sleep_alone"]["hr"])[:12]:
    v = rows[k]
    print(f"{k:<28}{v['sleep_alone']['hr']:>13.3f}{v['sleep_adjusted_for_awake']['hr']:>15.3f}"
          f"{v['excess_retained']:>10.2f}{v['penalized_frozen']['sleep_adjusted']['hr']:>9.3f}")

json.dump({"_built_by": "New_Figures/_NOT_THESE_workfiles/scripts/fits_sleep_unadjusted.py",
           "n": int(len(m)),
           "n_organ_metab_free": int(len(mC)),
           "cohort_file": "data_frozen_v8_2026-09/t90_final.parquet via cohort_spec.apply_cohort",
           "cohort_note": ("current 19,173 cohort, which gives 6,800 in this merged sample. "
                           "Until 2026-08-07 this script read t90_final_pre_icdfix.parquet and "
                           "returned 6,803, which violated the standing cohort rule and made "
                           "the two replacement negative controls unfittable here."),
           "impossible_oximetry_dropped": n_impossible,
           "exposure": "sleep-time T90, rank-inverse-normal within the sample, per 1 SD",
           "primary_specification": ("unpenalized Cox, stratified by site, natural cubic age "
                                     "spline with one column dropped, sex"),
           "frozen_specification": ("penalizer=0.01 with all 4 spline columns, which is what "
                                    "numbers/run_wake_sleep_healthy.py used and what "
                                    "wake_sleep_healthy_A_full.csv contains"),
           "validation": "the penalized refit reproduces wake_sleep_healthy_A_full.csv",
           "median_penalty_shrinkage": round(float(np.median(shr)), 4),
           "n_significant_sleep_alone": len(sig),
           "n_still_significant_after_awake_adjustment": len(surv),
           "lost_on_adjustment": [k for k in sig if k not in surv],
           "median_excess_retained": round(float(np.median(ret)), 3),
           "skipped": skipped,
           "outcomes": rows},
          open(f"{paths.NUMBERS_DIR}/sleep_unadjusted_v1.json", "w"), indent=1)
print("\nwrote numbers/sleep_unadjusted_v1.json")
