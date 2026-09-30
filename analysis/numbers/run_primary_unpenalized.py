"""
Refit every primary association without the ridge penalty.

lifelines multiplies its elastic-net term by the sample size, so `penalizer=0.01` on 19,000
patients is an effective ridge of roughly 190, not 0.01. Every hazard ratio published from that
specification is shrunk toward the null by an amount that scales inversely with the event count,
which also reorders the panel. The penalty was never disclosed and the numbers do not reproduce
in standard survival software.

The age basis is a natural cubic spline with 4 degrees of freedom. Its four columns sum to one,
which is harmless in a Cox model because the partial likelihood has no intercept, and the
unpenalized fits converge. One column is nevertheless dropped so the design is full rank under
any solver.

v8 (2026-09-12): the table comes from cohort_spec (T90_FINAL) and the exposure column from
cohort_spec.t90_column(). With T90_COLUMN=spo2_pct_below_90_sleep (step 231, decision 2) the T90
column is the sleep-period T90 (pre-PAP sleep on split nights, missing on the 482 collapsed-staging
nights) and the outputs carry the _sleepT90 suffix; the AHI and TST columns are fitted on exactly
the rows of the primary run (they do not lose the rows where the sensitivity exposure is missing)
and must equal numbers/bdsp_diseases_v3.csv, which the script checks before writing.

Outputs numbers/primary_unpenalized{SFX}.json and numbers/bdsp_diseases_v3{SFX}.csv.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../..")))  # repository root, where paths.py lives
import paths
import json
import sys
import warnings

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from patsy import dmatrix
from scipy import stats

warnings.filterwarnings("ignore")
ROOT = paths.T90_ROOT
sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
from cohort_spec import (apply_cohort, T90_FINAL, OUT_DIR, NUMBERS_DIR, PRIMARY_EXPOSURE,
                         t90_column, exposure_tag, sidecar, cohort_tag)
from disease_definitions import DISEASES, NEGATIVE_CONTROLS

XCOL = t90_column()
SFX = exposure_tag(XCOL) if XCOL != PRIMARY_EXPOSURE else ""
b = apply_cohort(pd.read_parquet(T90_FINAL))
print(f"cohort {len(b):,}   exposure {XCOL}   present on {int(b[XCOL].notna().sum()):,}")

b["male"] = (b.sex.astype(str).str.upper().str[0] == "M").astype(int)
sp = dmatrix("cr(a, df=4) - 1", {"a": b.AgeAtVisit.values}, return_type="dataframe").iloc[:, 1:]
for i in range(sp.shape[1]):
    b[f"age_s{i}"] = sp.iloc[:, i].values
ADJ = [f"age_s{i}" for i in range(sp.shape[1])] + ["male"]


def rint(x):
    """Rank inverse normal over the non-missing values; missing stays missing (the primary T90
    has no missing value in the cohort, so this equals the original on the primary run)."""
    x = pd.Series(x)
    out = pd.Series(np.nan, index=x.index)
    ok = x.notna()
    r = stats.rankdata(x[ok])
    out[ok] = stats.norm.ppf((r - 0.375) / (ok.sum() + 0.25))
    return out


b["z"] = b.groupby("site_id")[XCOL].transform(rint)
# the apnea-hypopnea index and total sleep time are carried on the same basis, because the
# tables and figures compare all three side by side and a mixed penalty would not be comparable
b["z_ahi"] = b.groupby("site_id").AHI.transform(lambda x: rint(x.fillna(x.median())))
b["z_tst"] = b.groupby("site_id").TST_min.transform(rint)   # v8.1: no median fill, TST is missing by design on split nights, so the tst exposure fits on the nights with a real TST
EXPOSURES = {"t90": "z", "ahi": "z_ahi", "tst": "z_tst"}

rows = []
for key, (label, is_neg, _, _) in list(DISEASES.items()) + [("death", ("Death from any cause", False, [], []))]:
    yc, ec, pc = f"{key}_years", f"{key}_incident", f"{key}_prevalent"
    if yc not in b.columns:
        continue
    f = b[(b[pc] == 0) & b[yc].notna() & (b[yc] > 0)]
    d = pd.DataFrame({"T": f[yc], "E": f[ec].astype(int), "site": f.site_id,
                      **{c: f[c] for c in list(EXPOSURES.values()) + ADJ}})
    d = d.dropna(subset=["T", "E", "site", "z_ahi"] + ADJ)   # the primary row set; z may be missing (sensitivity); v8.1: z_tst missing on split nights, dropped per exposure below
    if d.E.sum() < 60:
        continue
    rec = {"key": key, "disease": label, "negative_control": key in NEGATIVE_CONTROLS,
           "events": int(d.E.sum()), "n": int(len(d))}
    for name, col in EXPOSURES.items():
        dd = d.dropna(subset=[col])
        if name == "t90":
            rec["t90_n"], rec["t90_events"] = int(len(dd)), int(dd.E.sum())
        # the penalized refit is kept only for the t90 exposure, to document the size of the bias
        for tag, pen in ((("", 0.0), ("_pen", 0.01)) if name == "t90" else (("", 0.0),)):
            try:
                keep = [col] + ADJ + ["T", "E", "site"]
                c = CoxPHFitter(penalizer=pen).fit(dd[keep], "T", "E", strata=["site"])
                r = c.summary.loc[col]
                rec[f"{name}_hr{tag}"] = float(r["exp(coef)"])  # v8.3 full precision (2026-09-16)
                rec[f"{name}_lo{tag}"] = float(r["exp(coef) lower 95%"])
                rec[f"{name}_hi{tag}"] = float(r["exp(coef) upper 95%"])
                rec[f"{name}_p{tag}"] = float(r["p"])
            except Exception:
                rec[f"{name}_hr{tag}"] = np.nan
    rows.append(rec)

df = pd.DataFrame(rows)
df["shift_pct"] = (df.t90_hr / df.t90_hr_pen - 1) * 100  # v8.3 full precision (2026-09-16)
df["exposure"] = XCOL

CTAG = cohort_tag()   # v8.1: '_fullnights' under T90_FULL_NIGHTS_ONLY (the supplementary table without split nights), else ''
if SFX and CTAG:
    print(f"note: cohort tag {CTAG}, the AHI/TST positive control against the primary file is skipped (different cohort)")
if SFX and not CTAG:
    # the sensitivity set changes only the exposure: the AHI and TST columns must equal the primary file
    prim = pd.read_csv(f"{NUMBERS_DIR}/bdsp_diseases_v3.csv", comment="#")
    j = df.merge(prim, on="key", suffixes=("", "_prim"))
    assert len(j) == len(df) == len(prim), f"row sets differ from the primary file: {len(df)} vs {len(prim)}"
    cols = ["ahi_hr", "ahi_lo", "ahi_hi", "tst_hr", "tst_lo", "tst_hi", "events", "n"]
    bad = [c for c in cols if not np.allclose(j[c].astype(float), j[c + "_prim"].astype(float), atol=5e-4, equal_nan=True)]
    assert not bad, f"AHI/TST columns differ from the primary file on {bad}"
    print("positive control: AHI, TST, events and n equal the primary bdsp_diseases_v3.csv on every row")

csv_path, json_path = f"{OUT_DIR}/bdsp_diseases_v3{SFX}{CTAG}.csv", f"{OUT_DIR}/primary_unpenalized{SFX}{CTAG}.json"
df.to_csv(csv_path, index=False)
sig = df[(df.t90_lo > 1) | (df.t90_hi < 1)]
neg = df[df.negative_control]
out = {"n": int(len(b)), "exposure": XCOL, "n_conditions": int(len(df)), "n_significant": int(len(sig)),
       "negative_control_range": [float(neg.t90_hr.min()), float(neg.t90_hr.max())],
       "median_shift_pct": float(df.shift_pct.median())}
json.dump(out, open(json_path, "w"), indent=1)
for p in (csv_path, json_path):
    sidecar(p, __file__, note=f"run_primary_unpenalized exposure {XCOL}, {len(df)} rows, cohort {len(b):,}{CTAG}")

print(f"\n{len(df)} outcomes, {len(sig)} significant")
print(f"median shift from removing the penalty: {df.shift_pct.median():+.1f}%")
print(f"negative controls now span {neg.t90_hr.min():.2f} to {neg.t90_hr.max():.2f}\n")
print(f"{'condition':26s}{'penalized':>11}{'unpenalized':>25}{'shift':>8}{'ev':>7}")
print("-" * 80)
for _, r in df.sort_values("t90_hr", ascending=False).head(14).iterrows():
    print(f"{r.disease:26s}{r.t90_hr_pen:11.2f}   {r.t90_hr:.2f} [{r.t90_lo:.2f}-{r.t90_hi:.2f}]"
          f"{r.shift_pct:11.1f}%{r.events:7,}")
print("\nnegative controls")
for _, r in neg.sort_values("t90_hr", ascending=False).iterrows():
    print(f"{r.disease:26s}{r.t90_hr_pen:11.2f}   {r.t90_hr:.2f} [{r.t90_lo:.2f}-{r.t90_hi:.2f}]"
          f"{r.shift_pct:11.1f}%{r.events:7,}")
