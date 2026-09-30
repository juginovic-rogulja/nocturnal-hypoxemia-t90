#!/usr/bin/env python3
"""Placebo work list, time-matched to the REAL Design-A switch-on distribution.

Each PAP-off night is cut at a minute drawn from the switch-on times actually observed
in the usable split nights, so the control carries the same back-loading of REM that the
treated nights carry by construction.
"""
import glob, numpy as np, pandas as pd
from pathlib import Path

HERE = Path(__file__).parent
STG = HERE.parent / "t90_by_stage"
rng = np.random.default_rng(20260729)

A = pd.concat([pd.read_csv(f, low_memory=False) for f in sorted(glob.glob(str(HERE / "outA_*.csv")))],
              ignore_index=True)
A = A[(A.status == "ok") & np.isfinite(A.cpap_start_min) & (A.pre_any_min >= 45)
      & (A.post_any_min >= 45) & (A.pap_frac_pre < 0.10) & (A.pap_frac_post >= 0.80)]
sw = A.cpap_start_min.values
print("usable split nights %d | switch-on min median %.0f IQR %.0f-%.0f"
      % (len(A), np.median(sw), np.percentile(sw, 25), np.percentile(sw, 75)))

d = pd.concat([pd.read_csv(f, low_memory=False) for f in sorted(glob.glob(str(STG / "out_*.csv")))],
              ignore_index=True)
d = d[(d.status == "ok") & d.cpap_frac.notna() & (d.cpap_frac < 0.02)
      & (d.tst_min >= 60) & (d.rec_min >= 150)].copy()
print("PAP-off nights available:", len(d), d.SiteID.value_counts().to_dict())

# match the SITE mix of the treated nights so the control is not a different hospital mix
mix = A.SiteID.value_counts(normalize=True)
N = 1600
take = []
for s, p in mix.items():
    pool = d[d.SiteID == s]
    k = min(int(round(p * N)), len(pool))
    take.append(pool.sample(n=k, random_state=11))
P = pd.concat(take, ignore_index=True)
P["split_min"] = np.clip(rng.choice(sw, size=len(P)), 45,
                         np.maximum(P.rec_min.values - 45, 46))
P["arm"] = "P"
cols = ["BDSPPatientID", "SiteID", "BIDSFolder", "SessionID", "AgeAtVisit", "SexDSC",
        "AHI", "TST_min", "spo2_pct_below_90", "split_min", "arm"]
P[cols].to_csv(HERE / "worklist_P.csv", index=False)
print("placebo work list:", len(P), P.SiteID.value_counts().to_dict(),
      "split_min median %.0f" % P.split_min.median())
