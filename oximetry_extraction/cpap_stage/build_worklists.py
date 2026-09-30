#!/usr/bin/env python3
"""Three work lists for the CPAP-by-stage question.

A  split nights      : cpap_frac strictly between 0.05 and 0.95 in the 27,658-session
                       per-stage sweep -> PAP came on part-way through.
P  placebo nights    : PAP never on (cpap_frac < 0.02) and long enough to split.
                       Each gets a split_min DRAWN FROM the treated switch-time
                       distribution so the control is time-matched, not fixed at 179.
B  diagnostic vs titration pairs : the 2,470 two-session patients already assembled by
                       pilots/cpap_t90 (pair_sessions.csv).
"""
import numpy as np, pandas as pd, glob
from pathlib import Path

HERE = Path(__file__).parent
STG = HERE.parent / "t90_by_stage"
CP = HERE.parent / "cpap_t90"
rng = np.random.default_rng(20260729)

d = pd.concat([pd.read_csv(f, low_memory=False) for f in sorted(glob.glob(str(STG / "out_*.csv")))],
              ignore_index=True)
d = d[d.status == "ok"].copy()
print("per-stage sweep sessions:", len(d))

cols = ["BDSPPatientID", "SiteID", "BIDSFolder", "SessionID", "AgeAtVisit", "SexDSC",
        "AHI", "TST_min", "spo2_pct_below_90"]

# ---------------------------------------------------------------- A: split nights
A = d[(d.cpap_frac > 0.15) & (d.cpap_frac < 0.92) & (d.tst_min >= 60)
      & (d.rec_min >= 120)].copy()
A["cpap_frac_prior"] = A.cpap_frac
A = A.sort_values("cpap_frac")
print("A candidates:", len(A), A.SiteID.value_counts().to_dict())
A[cols + ["cpap_frac_prior"]].assign(arm="A").to_csv(HERE / "worklist_A.csv", index=False)

# ---------------------------------------------------------------- P: placebo nights
# split time sampled from the observed PAP switch-on distribution of the 65 verified
# split nights; resampled again after run A with the real distribution if it differs.
prior = pd.read_csv(CP / "designB_splitnights.csv", low_memory=False)
sw = prior.cpap_start_min.dropna().values
sw = sw[(sw > 30) & (sw < 400)]
print("prior switch-on minutes: n=%d median=%.0f IQR %.0f-%.0f"
      % (sw.size, np.median(sw), np.percentile(sw, 25), np.percentile(sw, 75)))

P = d[(d.cpap_frac.notna()) & (d.cpap_frac < 0.02) & (d.tst_min >= 60)].copy()
P = P.sample(n=min(1400, len(P)), random_state=7)
draw = rng.choice(sw, size=len(P))
# keep both sides >= 45 min of recording
P["split_min"] = np.clip(draw, 45, np.maximum(P.rec_min.values - 45, 46))
P = P[P.rec_min >= 150]
print("P candidates:", len(P), P.SiteID.value_counts().to_dict())
P[cols + ["split_min"]].assign(arm="P").to_csv(HERE / "worklist_P.csv", index=False)

# ---------------------------------------------------------------- B: night pairs
B = pd.read_csv(CP / "pair_sessions.csv")
B["arm"] = "B"
print("B sessions:", len(B), B.BDSPPatientID.nunique(), "patients",
      B.SiteID.value_counts().to_dict())
B.to_csv(HERE / "worklist_B.csv", index=False)
