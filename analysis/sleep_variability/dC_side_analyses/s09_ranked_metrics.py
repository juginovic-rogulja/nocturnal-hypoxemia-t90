"""
RANKED_METRICS.csv, the one-sheet view Alen asked for.

One row per measurement, ordered by the published dC rank of ranking_v3.csv, joining the
197x48 HR matrix to the hand-classified labels of measure_families.csv. Built entirely from
finished outputs, computes nothing new.

Column definitions
  rank                    published dC rank, ranking_v3.csv
  measurement             the feature key, kept because two labels can collide (siteZ twins)
  label                   human definition from measure_families.csv
  family                  measure_families.csv family
  dC                      published discrimination gain, ranking_v3.csv
  n_sig_fdr_hr_gt1_of48   outcomes FDR-significant (q<0.05 within measurement) with HR > 1
  median_sig_HR           median HR across those HR > 1 significant outcomes
  strongest_HR            the FDR-significant HR farthest from 1 in EITHER direction,
                          max |log HR|, so protective-coded measures (nadir, mean
                          saturation, sleep efficiency) surface their real strongest hit
                          rather than a direction artifact
  strongest_HR_outcome    the outcome that HR belongs to, human label
  rank_by_sig_count       position when the 197 are re-ordered by FDR-significant count in
                          either direction (ties broken by count above 1, then dC)
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import os

import numpy as np
import pandas as pd

from common import DISEASES, T90ROOT
import sys as _sys_v8; _sys_v8.path.insert(0, paths.ANALYSIS_DIR)  # v8 sweep 2026-09-12
from cohort_spec import N_MEASURES, N_RANKED_OUTCOMES  # v8 sweep 2026-09-12

HERE = os.path.dirname(os.path.abspath(__file__))
QTHR = 0.05


def label_of(k):
    return DISEASES[k][0] if k in DISEASES else "Death from any cause"


M = pd.read_csv(os.path.join(HERE, f"hr_matrix_{N_MEASURES}x{N_RANKED_OUTCOMES}.csv"))
RK = pd.read_csv(f"{paths.NUMBERS_DIR}/ranking_v3.csv", comment="#")
FAM = pd.read_csv(f"{paths.NUMBERS_DIR}/measure_families.csv")

assert M.measurement.nunique() == N_MEASURES, M.measurement.nunique()
fitted = M[M.p.notna()].copy()

rows = []
for f, g in fitted.groupby("measurement"):
    sig = g[g.q < QTHR]
    sig_gt = sig[sig.HR > 1]
    rec = {"measurement": f,
           "n_sig_fdr_hr_gt1_of48": int(len(sig_gt)),
           "median_sig_HR": round(float(sig_gt.HR.median()), 3) if len(sig_gt) else np.nan,
           "n_sig_either": int(len(sig))}
    if len(sig):
        s = sig.loc[(sig.HR.apply(lambda h: abs(np.log(h)))).idxmax()]
        rec["strongest_HR"] = round(float(s.HR), 3)
        rec["strongest_HR_outcome"] = label_of(s.outcome)
    else:
        rec["strongest_HR"] = np.nan
        rec["strongest_HR_outcome"] = ""
    rows.append(rec)
H = pd.DataFrame(rows)

df = (RK[["rank", "feature", "dC"]]
      .rename(columns={"feature": "measurement"})
      .merge(FAM[["feature", "family", "definition"]]
             .rename(columns={"feature": "measurement", "definition": "label"}),
             on="measurement", how="left")
      .merge(H, on="measurement", how="left"))
assert df.label.notna().all(), df[df.label.isna()].measurement.tolist()
assert len(df) == N_MEASURES

by_sig = df.sort_values(["n_sig_either", "n_sig_fdr_hr_gt1_of48", "dC"],
                        ascending=False).reset_index(drop=True)
by_sig["rank_by_sig_count"] = range(1, len(by_sig) + 1)
df = df.merge(by_sig[["measurement", "rank_by_sig_count"]], on="measurement", how="left")

cols = ["rank", "measurement", "label", "family", "dC", "n_sig_fdr_hr_gt1_of48",
        "median_sig_HR", "strongest_HR", "strongest_HR_outcome", "rank_by_sig_count",
        "n_sig_either"]
df = df[cols].sort_values("rank").reset_index(drop=True)
df["dC"] = df.dC.round(6)
out = os.path.join(HERE, "RANKED_METRICS.csv")
df.to_csv(out, index=False)
print(f"wrote {out}  ({len(df)} rows)")
print("\ntop 12 by published dC rank")
print(df.head(12).to_string(index=False, max_colwidth=46))
print("\nspot checks")
t = df[df.measurement == "spo2_pct_below_90"].iloc[0]
print(f"  T90: rank {int(t['rank'])}, sig>1 {int(t.n_sig_fdr_hr_gt1_of48)} of {N_RANKED_OUTCOMES}, "
      f"median {t.median_sig_HR}, strongest {t.strongest_HR} ({t.strongest_HR_outcome}), "
      f"sig-count rank {int(t.rank_by_sig_count)}")
n = df[df.measurement == "spo2_nadir_corrected"].iloc[0]
print(f"  nadir: rank {int(n['rank'])}, sig>1 {int(n.n_sig_fdr_hr_gt1_of48)}, "
      f"strongest {n.strongest_HR} ({n.strongest_HR_outcome}) "
      f"<- protective-coded, strongest is below 1 by design")
