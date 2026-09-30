"""
THE FULL MATRIX. Hazard ratios for all 197 ranked measurements against all 48 ranked
outcomes, 9,456 unpenalized site-stratified Cox fits.

Conventions are identical to s07_hazard_ratios.py, which is gated on reproducing the
published panel:
  unpenalized (penalizer 0.0), strata site_id
  age cr(df=4) on the 19,173 cohort with the FIRST column dropped (the rank trap), plus sex
  exposure as the within-site rank inverse normal, so every HR is per 1 SD rank-normal
  eligibility per outcome: prevalent == 0, follow-up present and > 0, complete cases
  60-event floor per fit, numbers/run_primary_unpenalized.py's line
  every fit at four Newton step sizes, maximum partial log-likelihood kept
  BH q within each measurement, across its fitted outcomes

MEASUREMENT LIST. Derived exactly as build_ranking_v3.py derives it, then ASSERTED equal to
ranking_v3.csv's 197 features, so the matrix and the ranking cannot drift apart. Any
per-fit skip (events under 60 after complete cases, or a failed fit) is recorded in the
long output with an empty HR and a reason, and counted in the log.

POSITIVE CONTROL. The spo2_pct_below_90 column of the finished matrix must reproduce
numbers/bdsp_diseases_v3.csv for all 48 ranked outcomes at its 3-decimal rounding, HR and
both CI bounds and events and n. Fail closed.

RESUMABLE. One checkpoint file per measurement in _work/hr_matrix_ck/. Rerunning skips
finished measurements.

Outputs (this folder): hr_matrix_197x48.csv, hr_matrix_summary.csv,
_work/hr_matrix_run.json (runtime and control bookkeeping).
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import json
import os
import time
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from lifelines import CoxPHFitter
from patsy import dmatrix
import sys as _sys_v8; _sys_v8.path.insert(0, paths.ANALYSIS_DIR)  # v8 sweep 2026-09-12
from cohort_spec import N_MEASURES, N_RANKED_OUTCOMES, SENSITIVITY_ONLY, RANKING_DROP_SUFFIX  # v8 sweep 2026-09-12

from common import (MASTER, PARQUET, T90ROOT, WORK, apply_cohort,
                    ranking_outcomes, site_rint)

STEPS = [None, 0.1, 0.25, 0.5]
MIN_EVENTS = 60
QTHR = 0.05
HERE = os.path.dirname(os.path.abspath(__file__))
CK = os.path.join(WORK, "hr_matrix_ck")
os.makedirs(CK, exist_ok=True)
PKL = os.path.join(WORK, "hr_matrix_frame.pkl")
ADJ_U = ["age_s1", "age_s2", "age_s3", "male"]

T0 = time.time()
print("=" * 78)
print(f"FULL HR MATRIX, {N_MEASURES} measurements x {N_RANKED_OUTCOMES} outcomes, unpenalized site-stratified Cox")
print("=" * 78, flush=True)

# ------------------------------------------------------------------ the 197, as the ranking
RK = pd.read_csv(f"{paths.NUMBERS_DIR}/ranking_v3.csv", comment="#")
RANKED = set(RK.feature)
assert len(RANKED) == N_MEASURES, len(RANKED)

base = apply_cohort(pd.read_parquet(PARQUET))
head = pd.read_csv(MASTER, nrows=0).columns.tolist()
FEAT = [c for c in head if c not in ("BDSPPatientID",)
        and not c.endswith(("_first_date", "_date"))]
raw = pd.read_csv(MASTER, usecols=["BDSPPatientID"] + FEAT, low_memory=False)
left = base[["BDSPPatientID", "site_id", "AgeAtVisit", "sex"] +
            [c for c in base.columns if c.endswith(("_incident", "_years", "_prevalent"))]]
raw = raw.drop(columns=[c for c in raw.columns
                        if c in set(left.columns) - {"BDSPPatientID"}])
FEAT = [c for c in FEAT if c in raw.columns]
d = left.merge(raw, on="BDSPPatientID", how="left")
d["male"] = (d.sex.astype(str).str.upper().str[0] == "M").astype(int)

num = [c for c in FEAT if pd.api.types.is_numeric_dtype(d[c])]
num = [c for c in num if d[c].notna().mean() >= 0.60 and d[c].nunique() > 10]
num = [c for c in num if c not in SENSITIVITY_ONLY and not c.endswith(RANKING_DROP_SUFFIX)]   # v8: the primary ranking's candidate rule (build_ranking_v3)
num = [c for c in num if c not in {"male", "AgeAtVisit"}
       and not c.startswith(("obs_", "follow_"))]
extra = sorted(set(num) - RANKED)
missing = sorted(RANKED - set(num))
assert not extra and not missing, (
    f"measurement list drifted from ranking_v3.csv: extra {extra}, missing {missing}")
num = [f for f in RK.feature]           # the ranking's own order, 197
print(f"measurements: {len(num)}, identical to ranking_v3.csv's set", flush=True)

OUT = ranking_outcomes(d)
assert len(OUT) == N_RANKED_OUTCOMES
print(f"outcomes: {len(OUT)}   fits: {len(num) * len(OUT):,} "
      f"(x{len(STEPS)} step sizes)", flush=True)

sp = dmatrix("cr(a, df=4) - 1", {"a": d.AgeAtVisit.values}, return_type="dataframe")
for i in range(4):
    d[f"age_s{i}"] = sp.iloc[:, i].values
for f in num:
    d[f + "__z"] = site_rint(d, f)

keep = ADJ_U + [f + "__z" for f in num] + ["site_id"] + \
    [f"{o}_{s}" for o in OUT for s in ("years", "incident", "prevalent")]
d[[c for c in dict.fromkeys(keep) if c in d.columns]].to_pickle(PKL)
print(f"frame pickled, {time.time() - T0:.0f}s so far", flush=True)

_CACHE = {}


def run_measurement(feat, outcomes):
    """All 48 outcomes for one measurement. Checkpointed, atomic, resumable."""
    ck = os.path.join(CK, f"{feat}.json")
    if os.path.exists(ck):
        try:
            return json.load(open(ck))
        except Exception:
            pass
    if "d" not in _CACHE:
        _CACHE["d"] = pd.read_pickle(PKL)
    dd = _CACHE["d"]
    z = feat + "__z"
    rows = []
    for out in outcomes:
        yc, ec = f"{out}_years", f"{out}_incident"
        g = dd[(dd[f"{out}_prevalent"] == 0) & dd[yc].notna() & (dd[yc] > 0)]
        f = g[[z] + ADJ_U + [yc, ec, "site_id"]].dropna()
        ev = int(f[ec].sum())
        if ev < MIN_EVENTS:
            rows.append({"measurement": feat, "outcome": out, "events": ev,
                         "n": int(len(f)), "skip_reason": "events_below_60"})
            continue
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
            rows.append({"measurement": feat, "outcome": out, "events": ev,
                         "n": int(len(f)), "skip_reason": "fit_failed"})
            continue
        _ll, step, c = best
        r = c.summary.loc[z]
        rows.append({"measurement": feat, "outcome": out,
                     "HR": float(r["exp(coef)"]),
                     "lo": float(r["exp(coef) lower 95%"]),
                     "hi": float(r["exp(coef) upper 95%"]),
                     "p": float(r["p"]), "events": ev, "n": int(len(f)),
                     "step": "default" if step is None else step,
                     "ll_spread": float(max(lls) - min(lls)) if len(lls) > 1 else 0.0})
    tmp = ck + ".part"
    json.dump(rows, open(tmp, "w"))
    os.replace(tmp, ck)
    return rows


done = [f for f in num if os.path.exists(os.path.join(CK, f"{f}.json"))]
if done:
    print(f"resuming, {len(done)} of {len(num)} measurements already checkpointed",
          flush=True)
t1 = time.time()
rr = Parallel(n_jobs=7, verbose=5, backend="loky")(
    delayed(run_measurement)(f, OUT) for f in num)
fit_minutes = (time.time() - t1) / 60
print(f"fitting took {fit_minutes:.1f} min", flush=True)

# ------------------------------------------------------------------ assemble
df = pd.DataFrame([r for chunk in rr for r in chunk])
assert len(df) == len(num) * len(OUT), (len(df), len(num) * len(OUT))
df["q"] = np.nan
for f in num:
    m = (df.measurement == f) & df.p.notna()
    if m.sum():
        p = df.loc[m, "p"].values
        n = len(p)
        o = np.argsort(p)
        q = np.empty(n)
        q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
        df.loc[m, "q"] = np.minimum(q, 1.0)

n_skip = int(df.p.isna().sum())
skips = df[df.p.isna()]
print(f"skipped fits: {n_skip} of {len(df)}", flush=True)
if n_skip:
    print(skips.groupby("skip_reason").size().to_string(), flush=True)

# ------------------------------------------------------------------ positive control
print("\nPOSITIVE CONTROL against numbers/bdsp_diseases_v3.csv", flush=True)
ref = pd.read_csv(f"{paths.NUMBERS_DIR}/bdsp_diseases_v3.csv")
mine = df[(df.measurement == "spo2_pct_below_90") & df.p.notna()]
j = mine.merge(ref[["key", "t90_hr", "t90_lo", "t90_hi", "events", "n"]],
               left_on="outcome", right_on="key", how="inner", suffixes=("", "_ref"))
assert len(j) == N_RANKED_OUTCOMES, len(j)
bad = j[((j.HR - j.t90_hr).abs() > 5e-4) |
        ((j.lo - j.t90_lo).abs() > 5e-4) |
        ((j.hi - j.t90_hi).abs() > 5e-4) |  # v8.3: reference now full precision, compare within 5e-4 (2026-09-16)
        (j.events != j.events_ref) | (j.n != j.n_ref)]
print(f"  {len(j)} of {N_RANKED_OUTCOMES} ranked outcomes, mismatches: {len(bad)}")
if len(bad):
    print(bad[["outcome", "HR", "t90_hr"]].to_string(index=False))
    raise SystemExit("POSITIVE CONTROL FAILED, matrix not written")
print("  PASS: the published unpenalized panel reproduces exactly", flush=True)

# ------------------------------------------------------------------ outputs
long_cols = ["measurement", "outcome", "HR", "lo", "hi", "p", "q", "events", "n",
             "skip_reason", "step", "ll_spread"]
for c in long_cols:
    if c not in df.columns:
        df[c] = np.nan
out_long = os.path.join(HERE, f"hr_matrix_{N_MEASURES}x{N_RANKED_OUTCOMES}.csv")
d3 = df.copy()
for c in ("HR", "lo", "hi"):
    pass  # v8.3: hr_matrix_141x52.csv keeps full precision (was d3[c].round(3)); Supp Table 7 prints its HRs at 2 dp
d3[long_cols].to_csv(out_long, index=False)
print(f"\nwrote {out_long}  ({len(d3)} rows)", flush=True)

summ = []
for f in num:
    g = df[(df.measurement == f) & df.p.notna()]
    sig_gt = g[(g.q < QTHR) & (g.HR > 1)]
    sig_lt = g[(g.q < QTHR) & (g.HR < 1)]
    summ.append({"measurement": f, "n_fitted": int(len(g)),
                 "n_sig_fdr_hr_gt1": int(len(sig_gt)),
                 "median_sig_HR_gt1": round(float(sig_gt.HR.median()), 3) if len(sig_gt) else np.nan,
                 "max_sig_HR_gt1": round(float(sig_gt.HR.max()), 3) if len(sig_gt) else np.nan,
                 "n_sig_fdr_hr_lt1": int(len(sig_lt)),
                 "median_sig_HR_lt1": round(float(sig_lt.HR.median()), 3) if len(sig_lt) else np.nan,
                 "n_sig_either": int(len(sig_gt) + len(sig_lt))})
S = pd.DataFrame(summ).merge(
    RK[["feature", "rank", "dC"]].rename(columns={"feature": "measurement",
                                                  "rank": "published_dC_rank",
                                                  "dC": "published_dC"}),
    on="measurement", how="left")
S = S.sort_values(["n_sig_either", "n_sig_fdr_hr_gt1"],
                  ascending=False).reset_index(drop=True)
S.insert(0, "sig_count_rank", range(1, len(S) + 1))
out_sum = os.path.join(HERE, "hr_matrix_summary.csv")
S.to_csv(out_sum, index=False)
print(f"wrote {out_sum}  ({len(S)} rows)", flush=True)

# ------------------------------------------------------------------ report block
print("\nTOP 10 BY FDR-SIGNIFICANT COUNT (either direction), against the published dC rank")
print(f"{'measurement':<30}{'sig':>5}{'>1':>4}{'<1':>4}{'med HR>1':>10}{'dC rank':>9}"
      f"{'dC':>10}")
print("-" * 74)
for _i, r in S.head(10).iterrows():
    mg = f"{r.median_sig_HR_gt1:.3f}" if pd.notna(r.median_sig_HR_gt1) else "-"
    print(f"{r.measurement:<30}{r.n_sig_either:>5}{r.n_sig_fdr_hr_gt1:>4}"
          f"{r.n_sig_fdr_hr_lt1:>4}{mg:>10}{int(r.published_dC_rank):>9}"
          f"{r.published_dC:>+10.5f}")

from scipy import stats as sst
rho, rp = sst.spearmanr(S.n_sig_either, -S.published_dC_rank)
rho2, rp2 = sst.spearmanr(S.n_sig_either, S.published_dC)
t90row = S[S.measurement == "spo2_pct_below_90"].iloc[0]
print(f"\nagreement, all {N_MEASURES}: Spearman(sig count, dC) {rho2:+.3f}, "
      f"Spearman(sig count, inverse dC rank) {rho:+.3f}")
print(f"T90 by significant count: position {int(t90row.sig_count_rank)} "
      f"({int(t90row.n_sig_either)} significant, {int(t90row.n_sig_fdr_hr_gt1)} above 1), "
      f"published dC rank 1")

nd = df[df.p.notna() & (df.step != "default")]
runtime_min = (time.time() - T0) / 60
print(f"\nnon-default step won on log-likelihood in {len(nd)} fits, "
      f"max log-likelihood spread across steps {df.ll_spread.max():.2e}")
print(f"total runtime {runtime_min:.1f} min (fitting {fit_minutes:.1f} min)")

json.dump({"runtime_min": round(runtime_min, 1), "fit_min": round(fit_minutes, 1),
           "n_measurements": len(num), "n_outcomes": len(OUT),
           "n_fits": int(len(df)), "n_skipped": n_skip,
           "skip_reasons": skips.groupby("skip_reason").size().to_dict() if n_skip else {},
           "positive_control_pass": True,
           "spearman_sigcount_vs_dC": round(float(rho2), 3),
           "max_ll_spread": float(df.ll_spread.max()),
           "t90_sig_either": int(t90row.n_sig_either),
           "t90_sig_count_position": int(t90row.sig_count_rank)},
          open(os.path.join(WORK, "hr_matrix_run.json"), "w"), indent=2)
print("run bookkeeping -> _work/hr_matrix_run.json", flush=True)
