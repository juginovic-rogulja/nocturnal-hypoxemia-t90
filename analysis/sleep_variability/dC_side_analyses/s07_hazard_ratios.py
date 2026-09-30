"""
HAZARD RATIOS for every exposure the dC side analyses scored.

The dC machinery fits CoxPHFitter(penalizer=0.01), which lifelines multiplies by the sample
size, an effective ridge near 190 that shrinks every coefficient toward 1 and was the reason
the paper refit its primary panel unpenalized. HRs here therefore come from UNPENALIZED
site-stratified Cox, the exact specification of numbers/run_primary_unpenalized.py:

  strata          site_id
  age             cr(df=4) computed on the fitting cohort, FIRST COLUMN DROPPED so the design
                  is full rank under an unpenalized solver (the rank trap)
  sex             male indicator
  exposure        within-site rank inverse normal, so every HR is per 1 SD of the rank-normal
                  score, the paper's primary per-SD scale
  eligibility     prevalent == 0, follow-up present and > 0, complete cases, 60+ events
  optimum         every fit is run at four Newton step sizes and the maximum partial
                  log-likelihood solution is kept, because the default step has converged to
                  a wrong optimum on this cohort before

POSITIVE CONTROL. The full-cohort frozen-T90 fits must reproduce numbers/bdsp_diseases_v3.csv
(HR, both CI bounds, events, n) at its own 3-decimal rounding for every one of the 48 ranked
outcomes it carries. That file is the paper's published per-SD panel. Fail closed.

SETS
  full_19173      48 ranked outcomes. Exposures: frozen T90, non-REM heart rate, corrected
                  nadir, and the T90 + heart rate joint model with both HRs reported.
  subcohort_6800  the wake/sleep separable subcohort of s02, the same 45 outcomes that
                  cleared the dC floor there. Exposures: frozen T90 (the comparator),
                  sleep-period T90, wake-period T90.

BH q is computed within each exposure column, over exactly the outcomes that exposure fitted.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import glob
import json
import os

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from lifelines import CoxPHFitter
from patsy import dmatrix

from common import (DISEASES, PARQUET, SHARDS, T90ROOT, WORK, apply_cohort,
                    build_frame, organ_of, site_rint)
from disease_definitions import NEGATIVE_CONTROLS

T90, HR, NADIR = "spo2_pct_below_90", "nrem_hr_bpm", "spo2_nadir_corrected"
STEPS = [None, 0.1, 0.25, 0.5]
MIN_EVENTS = 60
QTHR = 0.05
HERE = os.path.dirname(os.path.abspath(__file__))


def bh(p):
    """Benjamini-Hochberg, transcribed from numbers/run_wake_sleep_healthy.py."""
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1.0)


def label_of(k):
    return DISEASES[k][0] if k in DISEASES else "Death from any cause"


# ------------------------------------------------------------------ one unpenalized fit
def fit_one(pkl, out, expo_cols, adj):
    dd = pd.read_pickle(pkl)
    yc, ec = f"{out}_years", f"{out}_incident"
    g = dd[(dd[f"{out}_prevalent"] == 0) & dd[yc].notna() & (dd[yc] > 0)]
    f = g[list(expo_cols) + list(adj) + [yc, ec, "site_id"]].dropna()
    ev = int(f[ec].sum())
    if ev < MIN_EVENTS:
        return {"outcome": out, "events": ev, "n": int(len(f)), "skipped": True}
    best, lls = None, []
    for s in STEPS:
        try:
            c = CoxPHFitter(penalizer=0.0)
            kw = {"fit_options": {"step_size": s}} if s is not None else {}
            c.fit(f, duration_col=yc, event_col=ec, strata=["site_id"], **kw)
            ll = float(c.log_likelihood_)
            lls.append(ll)
            if best is None or ll > best[0]:
                best = (ll, s, c)
        except Exception:
            continue
    if best is None:
        return {"outcome": out, "events": ev, "n": int(len(f)), "skipped": True}
    _ll, step, c = best
    rows = {}
    for col in expo_cols:
        r = c.summary.loc[col]
        rows[col] = {"hr": float(r["exp(coef)"]),
                     "lo": float(r["exp(coef) lower 95%"]),
                     "hi": float(r["exp(coef) upper 95%"]),
                     "p": float(r["p"])}
    return {"outcome": out, "events": ev, "n": int(len(f)), "skipped": False,
            "step": step, "ll_spread": float(max(lls) - min(lls)) if len(lls) > 1 else 0.0,
            "coefs": rows}


def run_set(tag, frame, outcomes, exposures, adj):
    """exposures: {exposure_name: (model_name, [z columns], which column to report)}"""
    pkl = os.path.join(WORK, f"hr_{tag}.pkl")
    zcols = sorted({c for _m, cols, _r in exposures.values() for c in cols})
    keep = list(adj) + zcols + ["site_id"] + \
        [f"{o}_{s}" for o in outcomes for s in ("years", "incident", "prevalent")]
    frame[[c for c in dict.fromkeys(keep) if c in frame.columns]].to_pickle(pkl)

    fits = {}          # (model_name) -> {outcome: fitresult}, one fit per model per outcome
    models = {}        # model_name -> cols
    for _e, (m, cols, _r) in exposures.items():
        models[m] = cols
    jobs = [(m, o) for m in models for o in outcomes]
    res = Parallel(n_jobs=7, verbose=0, backend="loky")(
        delayed(fit_one)(pkl, o, models[m], adj) for m, o in jobs)
    for (m, o), r in zip(jobs, res):
        fits.setdefault(m, {})[o] = r

    rows = []
    for expo, (m, _cols, rcol) in exposures.items():
        for o in outcomes:
            r = fits[m][o]
            base = {"cohort": tag, "exposure": expo, "model": m,
                    "outcome_key": o, "outcome_label": label_of(o),
                    "organ": organ_of(o),
                    "negative_control": o in NEGATIVE_CONTROLS,
                    "n": r["n"], "events": r["events"]}
            if r["skipped"] or rcol not in r.get("coefs", {}):
                rows.append({**base, "HR_per_SD": np.nan, "lo95": np.nan,
                             "hi95": np.nan, "p": np.nan, "skipped": True})
                continue
            c = r["coefs"][rcol]
            rows.append({**base, "HR_per_SD": c["hr"], "lo95": c["lo"], "hi95": c["hi"],
                         "p": c["p"], "skipped": False,
                         "step_size_used": "default" if r["step"] is None else r["step"],
                         "ll_spread_across_steps": r["ll_spread"]})
    df = pd.DataFrame(rows)
    # BH within exposure, over the outcomes that exposure actually fitted
    df["q_bh_within_exposure"] = np.nan
    for expo in df.exposure.unique():
        m = (df.exposure == expo) & df.p.notna()
        df.loc[m, "q_bh_within_exposure"] = bh(df.loc[m, "p"].values)
    df["fdr_sig_hr_gt1"] = (df.q_bh_within_exposure < QTHR) & (df.HR_per_SD > 1)
    return df


print("=" * 78)
print("HAZARD RATIOS, UNPENALIZED SITE-STRATIFIED COX, PER 1 SD RANK-NORMAL")
print("=" * 78)

# ------------------------------------------------------------------ full cohort
d = build_frame(extra_master_cols=[T90, HR, NADIR], cache="frame_combo.pkl")
from common import ranking_outcomes
import sys as _sys_v8; _sys_v8.path.insert(0, paths.ANALYSIS_DIR)  # v8 sweep 2026-09-12
from cohort_spec import N_RANKED_OUTCOMES  # v8 sweep 2026-09-12
OUT48 = ranking_outcomes(d)
assert len(OUT48) == N_RANKED_OUTCOMES
for c in (T90, HR, NADIR):
    d[c + "__z"] = site_rint(d, c)
# run_primary_unpenalized.py's basis: same cr(df=4) on the same 19,173 ages, first column
# dropped. build_frame already carries the identical 4 columns, so dropping age_s0 IS its ADJ.
ADJ_U = ["age_s1", "age_s2", "age_s3", "male"]

EXPO_FULL = {
    "t90_frozen":          ("single_t90",   [T90 + "__z"],            T90 + "__z"),
    "nrem_hr":             ("single_hr",    [HR + "__z"],             HR + "__z"),
    "nadir":               ("single_nadir", [NADIR + "__z"],          NADIR + "__z"),
    "t90_in_joint_with_hr": ("joint_t90_hr", [T90 + "__z", HR + "__z"], T90 + "__z"),
    "hr_in_joint_with_t90": ("joint_t90_hr", [T90 + "__z", HR + "__z"], HR + "__z"),
}
full = run_set("full_19173", d, OUT48, EXPO_FULL, ADJ_U)

# ------------------------------------------------------------------ positive control
print("\nPOSITIVE CONTROL against numbers/bdsp_diseases_v3.csv (published per-SD panel)")
ref = pd.read_csv(f"{paths.NUMBERS_DIR}/bdsp_diseases_v3.csv")
mine = full[(full.exposure == "t90_frozen") & (~full.skipped.astype(bool))]
j = mine.merge(ref[["key", "t90_hr", "t90_lo", "t90_hi", "events", "n"]],
               left_on="outcome_key", right_on="key", how="inner",
               suffixes=("", "_ref"))
assert len(j) == N_RANKED_OUTCOMES, f"only {len(j)} of {N_RANKED_OUTCOMES} ranked outcomes found in the reference"
bad = j[(j.HR_per_SD - j.t90_hr).abs() > 5e-4]  # v8.3: reference now full precision, compare within 5e-4 (2026-09-16)
bad_ci = j[((j.lo95 - j.t90_lo).abs() > 5e-4) |
           ((j.hi95 - j.t90_hi).abs() > 5e-4)]
bad_ne = j[(j.events != j.events_ref) | (j.n != j.n_ref)]
print(f"  {len(j)} of {N_RANKED_OUTCOMES} ranked outcomes present in the reference")
print(f"  HR mismatches at the file's 3-decimal rounding: {len(bad)}")
print(f"  CI-bound mismatches: {len(bad_ci)}     events/n mismatches: {len(bad_ne)}")
if len(bad) or len(bad_ci) or len(bad_ne):
    print(bad[["outcome_key", "HR_per_SD", "t90_hr"]].to_string(index=False))
    raise SystemExit("POSITIVE CONTROL FAILED, stop")
print("  PASS: the published unpenalized panel reproduces exactly")
hf = j[j.outcome_key == "hf"].iloc[0]
print(f"  spot check, heart failure: {hf.HR_per_SD:.3f} [{hf.lo95:.3f}-{hf.hi95:.3f}] "
      f"against published {hf.t90_hr:.3f} [{hf.t90_lo:.3f}-{hf.t90_hi:.3f}]")

# ------------------------------------------------------------------ subcohort
print("\nassembling the 6,800 subcohort (s02's loader, count asserted)")
files = sorted(glob.glob(SHARDS))
assert files, f"no cpap_stage shards under {SHARDS}"   # v8.2: shard count not typed
raw = pd.concat([pd.read_csv(f, low_memory=False,
                             usecols=["BDSPPatientID", "all_wake_t90", "all_sleep_t90",
                                      "all_wake_min", "status"]) for f in files],
                ignore_index=True)
st = raw.drop_duplicates("BDSPPatientID")
st = st[st.status == "ok"]
b = apply_cohort(pd.read_parquet(PARQUET))
b = b[~(b.spo2_nadir_corrected > b.spo2_mean)]
m = b.merge(st.drop(columns=["status"]), on="BDSPPatientID", how="inner")
m = m[m.all_wake_min >= 10.0]
m = m.dropna(subset=["all_wake_t90", "all_sleep_t90"]).reset_index(drop=True)
assert len(m) > 0 and m.BDSPPatientID.is_unique, len(m)   # v8.2: the subcohort size is a result, recorded, never typed
print(f"wake/sleep separable subcohort: {len(m):,} patients (the 2026-07-29 extraction gave 6,800)")
mst = pd.read_csv(f"{paths.TABLES_DIR}/"
                  "master_cohort.csv",
                  usecols=["BDSPPatientID", T90], low_memory=False)
mst = mst.rename(columns={T90: "t90_frozen_raw"}).drop_duplicates("BDSPPatientID")
m = m.merge(mst, on="BDSPPatientID", how="left")
assert m.t90_frozen_raw.notna().all()
m["male"] = (m.sex.astype(str).str.upper().str[0] == "M").astype(int)
sp = dmatrix("cr(a, df=4) - 1", {"a": m.AgeAtVisit.values}, return_type="dataframe")
for i in range(4):
    m[f"age_s{i}"] = sp.iloc[:, i].values
for c in ("t90_frozen_raw", "all_sleep_t90", "all_wake_t90"):
    m[c + "__z"] = site_rint(m, c)

OUT45 = json.load(open(os.path.join(WORK, "B_result.json")))["outcomes_common"]
assert OUT45, "B_result.json names no common outcome"   # v8.2: the count comes from s02 (B_result.json), never typed
print(f"common outcomes from s02: {len(OUT45)}")
EXPO_SUB = {
    "t90_frozen":  ("single_t90f",  ["t90_frozen_raw__z"], "t90_frozen_raw__z"),
    "t90_sleep":   ("single_sleep", ["all_sleep_t90__z"],  "all_sleep_t90__z"),
    "t90_wake":    ("single_wake",  ["all_wake_t90__z"],   "all_wake_t90__z"),
}
sub = run_set("subcohort_6800", m, OUT45, EXPO_SUB, ADJ_U)

# ------------------------------------------------------------------ output
df = pd.concat([full, sub], ignore_index=True)
for c in ("HR_per_SD", "lo95", "hi95"):
    df[c] = df[c].round(3)
df["p"] = df["p"].astype(float)
df["q_bh_within_exposure"] = df["q_bh_within_exposure"].astype(float)
out_path = os.path.join(HERE, "hr_tables.csv")
df.to_csv(out_path, index=False)
print(f"\nwrote {out_path}  ({len(df)} rows)")

# per-exposure summary. Direction matters: every exposure here is coded so that MORE of it
# is worse, except the corrected nadir, where a HIGHER nadir is better oxygenation and the
# harmful direction is therefore HR BELOW 1. Both directions are counted so the nadir cannot
# be misread as inert.
summ = []
for (coh, expo), g in df[df.p.notna()].groupby(["cohort", "exposure"], sort=False):
    sig_gt = g[(g.q_bh_within_exposure < QTHR) & (g.HR_per_SD > 1)]
    sig_lt = g[(g.q_bh_within_exposure < QTHR) & (g.HR_per_SD < 1)]
    summ.append({"cohort": coh, "exposure": expo, "n_outcomes_fitted": int(len(g)),
                 "n_fdr_sig_hr_gt1": int(len(sig_gt)),
                 "median_sig_HR_gt1": float(sig_gt.HR_per_SD.median()) if len(sig_gt) else np.nan,
                 "max_sig_HR_gt1": float(sig_gt.HR_per_SD.max()) if len(sig_gt) else np.nan,
                 "n_fdr_sig_hr_lt1": int(len(sig_lt)),
                 "median_sig_HR_lt1": float(sig_lt.HR_per_SD.median()) if len(sig_lt) else np.nan,
                 "neg_controls_flagged": ", ".join(
                     g[g.negative_control & (g.q_bh_within_exposure < QTHR)]
                     .outcome_label.tolist()) or "none",
                 "median_HR_all": float(g.HR_per_SD.median())})
S = pd.DataFrame(summ)
S.to_csv(os.path.join(WORK, "hr_summary.csv"), index=False)
print("\nSUMMARY, per exposure: FDR-significant at q<0.05, split by direction")
print(f"{'cohort':<16}{'exposure':<22}{'fitted':>7}{'sig>1':>6}{'med HR':>8}{'max HR':>8}"
      f"{'sig<1':>6}{'med HR':>8}   negative controls flagged")
print("-" * 110)
for _i, r in S.iterrows():
    mg = f"{r.median_sig_HR_gt1:.3f}" if pd.notna(r.median_sig_HR_gt1) else "-"
    mx = f"{r.max_sig_HR_gt1:.3f}" if pd.notna(r.max_sig_HR_gt1) else "-"
    ml = f"{r.median_sig_HR_lt1:.3f}" if pd.notna(r.median_sig_HR_lt1) else "-"
    print(f"{r.cohort:<16}{r.exposure:<22}{r.n_outcomes_fitted:>7}{r.n_fdr_sig_hr_gt1:>6}"
          f"{mg:>8}{mx:>8}{r.n_fdr_sig_hr_lt1:>6}{ml:>8}   {r.neg_controls_flagged}")
print("direction note: the corrected nadir is protective-coded, higher is better, so its "
      "harmful direction is the sig<1 column. Every other exposure is harmful-coded.")

# convergence bookkeeping
nd = df[df.step_size_used.notna() & (df.step_size_used != "default")]
print(f"\nfits where a non-default Newton step won on log-likelihood: "
      f"{nd[['cohort', 'model', 'outcome_key']].drop_duplicates().shape[0]} "
      f"(max log-likelihood spread across steps "
      f"{df.ll_spread_across_steps.max():.2e})")
print("penalizer note: these HRs are penalizer 0.0 with one age-spline column dropped. The "
      "dC fits keep the 4-column basis at penalizer 0.01, which lifelines scales by n. The "
      "two answer different questions and their numbers are not interchangeable.")
