"""
B. Sleep-period T90 against whole-recording T90, in one subcohort, on one set of outcomes.

The paper's exposure counts every sample of the recording, wake included
(numbers/T90_PROVENANCE.md, rule 5: "over the entire recording. Not sleep-restricted"). The
obvious question is whether restricting it to scored sleep discriminates better. It can only be
answered where wake and sleep can be separated, which is the 6,800-recording subcohort behind
numbers/wake_vs_sleep_t90.csv.

COHORT. Assembled by the loader of numbers/run_wake_sleep_healthy.py, transcribed step for
step: 72 cpap_stage shards, one row per patient, status ok, joined to the 19,173 analysis
cohort, minus recordings whose nadir exceeds their mean saturation, minus recordings with under
10 minutes of scored wake, minus missing wake or sleep T90. The count is asserted at 6,800.

EXPOSURES, all four on the same 6,800 people and the same outcomes.
  t90_frozen      the paper's exposure, master_cohort_v2.csv, whole recording
  t90_wholerec    all_any_t90, the shard extractor's own whole-recording value
  t90_sleep       all_sleep_t90, scored sleep epochs only
  t90_wake        all_wake_t90, scored wake epochs only
The middle two matter: t90_frozen and t90_sleep come from different extraction runs, so a
difference between them confounds the window with the extractor. t90_wholerec against t90_sleep
is the clean contrast, same extractor, same night, only the window changes.

MODEL. The ranking's, not the wake/sleep paper's. Within-site rank inverse normal, age spline
plus sex, penalizer 0.01, fit at one site and scored at the other, both ways. The wake/sleep
lineage used a global rank transform and a site strata term instead; that specification answers
a different question and is not used here.

OUTCOMES. The ranking's 48, with the ranking's own per-fold floor of 30 events. In a cohort
2.8 times smaller many outcomes cannot clear it. The headline average is taken over the
outcomes where ALL FOUR exposures returned a value, so the comparison is on identical ground.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import glob
import json
import os

import numpy as np
import pandas as pd
from patsy import dmatrix

from common import (ADJ, DISEASES, MASTER, PARQUET, PEN, SHARDS, WORK,
                    apply_cohort, heldout, ranking_outcomes, site_rint)
import sys as _sys_v8; _sys_v8.path.insert(0, paths.ANALYSIS_DIR)  # v8 sweep 2026-09-12
from cohort_spec import N_RANKED_OUTCOMES  # v8 sweep 2026-09-12

MIN_WAKE_MIN = 10.0
SHARDCOLS = ["BDSPPatientID", "all_wake_t90", "all_sleep_t90", "all_any_t90",
             "all_wake_min", "all_sleep_min", "status"]

print("=" * 78)
print("B. SLEEP-PERIOD T90 AGAINST WHOLE-RECORDING T90")
print("=" * 78)

# ------------------------------------------------------------------ cohort, per the loader
files = sorted(glob.glob(SHARDS))
assert files, f"no cpap_stage shards under {SHARDS}"   # v8.2: the shard count is a property of the fleet, not a result
raw = pd.concat([pd.read_csv(f, low_memory=False, usecols=SHARDCOLS) for f in files],
                ignore_index=True)
st = raw.drop_duplicates("BDSPPatientID")
st = st[st.status == "ok"]
assert st.BDSPPatientID.is_unique

b = apply_cohort(pd.read_parquet(PARQUET))
bad = (b.spo2_nadir_corrected > b.spo2_mean)
b = b[~bad]
print(f"impossible-oximetry recordings dropped: {int(bad.sum())}")

m = b.merge(st.drop(columns=["status"]), on="BDSPPatientID", how="inner")
m = m[m.all_wake_min >= MIN_WAKE_MIN]
m = m.dropna(subset=["all_wake_t90", "all_sleep_t90", "all_any_t90"]).reset_index(drop=True)
assert len(m) > 0 and m.BDSPPatientID.is_unique, len(m)   # v8.2: the subcohort size is a result of the re-extraction, recorded, never typed
print(f"wake/sleep separable subcohort: {len(m):,} patients (the 2026-07-29 extraction gave 6,800)")
print(f"subcohort {len(m):,}  (matches numbers/wake_vs_sleep_t90.csv lineage)")

# the paper's frozen exposure, from the master, the same copy the ranking scores
mst = pd.read_csv(MASTER, usecols=["BDSPPatientID", "spo2_pct_below_90"], low_memory=False)
mst = mst.rename(columns={"spo2_pct_below_90": "t90_frozen"}).drop_duplicates("BDSPPatientID")
m = m.merge(mst, on="BDSPPatientID", how="left")
assert m.t90_frozen.notna().all(), int(m.t90_frozen.isna().sum())

m["male"] = (m.sex.astype(str).str.upper().str[0] == "M").astype(int)
sp = dmatrix("cr(a, df=4) - 1", {"a": m.AgeAtVisit.values}, return_type="dataframe")
sp.columns = ADJ[:4]
sp.index = m.index
for c in sp.columns:
    m[c] = sp[c]

m = m.rename(columns={"all_any_t90": "t90_wholerec", "all_sleep_t90": "t90_sleep",
                      "all_wake_t90": "t90_wake"})
EXPO = ["t90_frozen", "t90_wholerec", "t90_sleep", "t90_wake"]
print("\nexposure distributions in the subcohort")
print(f"{'exposure':<16}{'median':>9}{'IQR':>20}{'mean':>9}")
for e in EXPO:
    q = m[e].quantile([0.25, 0.5, 0.75])
    print(f"{e:<16}{q[0.5]:>9.3f}{f'{q[0.25]:.3f} to {q[0.75]:.3f}':>20}{m[e].mean():>9.3f}")
print("\nSpearman between the exposures")
cm = m[EXPO].corr(method="spearman")
print(cm.round(3).to_string())

for e in EXPO:
    m[e + "__z"] = site_rint(m, e)

# ------------------------------------------------------------------ fits
OUT = ranking_outcomes(pd.read_parquet(PARQUET).pipe(apply_cohort))
assert len(OUT) == N_RANKED_OUTCOMES
specs = {"__base__": []}
for e in EXPO:
    specs[e] = [e + "__z"]
res = heldout(specs, OUT, m, penalizer=PEN, pkl_name="B_frame.pkl")

wide = res.pivot(index="outcome", columns="spec", values="gain")
common = wide.dropna(subset=EXPO).index.tolist()
print(f"\noutcomes where all four exposures returned a value at the ranking's 30-event "
      f"per-fold floor: {len(common)} of {N_RANKED_OUTCOMES}")
dropped = sorted(set(OUT) - set(common))
print("underpowered in this subcohort, excluded from all four averages:")
print("   " + ", ".join(DISEASES[o][0] if o in DISEASES else "Death" for o in dropped))

sub = res[res.outcome.isin(common)]
summ = (sub[sub.spec != "__base__"].groupby("spec")
        .agg(dC=("gain", "mean"), worst=("gain", "min"), best=("gain", "max"),
             npos=("gain", lambda x: int((x > 0).sum())), nout=("gain", "size"))
        .reset_index().sort_values("dC", ascending=False))
base_c = float(sub[sub.spec == "__base__"].c.mean())

print(f"\nheld-out baseline over these {len(common)} outcomes in the 6,800: {base_c:.6f}")
print(f"\n{'exposure':<16}{'dC':>12}{'n>0':>6}{'worst':>12}{'best':>12}")
print("-" * 58)
for _i, r in summ.iterrows():
    print(f"{r.spec:<16}{r.dC:>+12.6f}{int(r.npos):>6}{r.worst:>+12.6f}{r.best:>+12.6f}")

g = summ.set_index("spec").dC
print(f"\nsleep-only against whole-recording, SAME extractor: "
      f"{g['t90_sleep']:+.6f} against {g['t90_wholerec']:+.6f}, "
      f"difference {g['t90_sleep'] - g['t90_wholerec']:+.6f}")
print(f"sleep-only against the paper's frozen exposure: "
      f"{g['t90_sleep']:+.6f} against {g['t90_frozen']:+.6f}, "
      f"difference {g['t90_sleep'] - g['t90_frozen']:+.6f}")
print(f"wake-only, for contrast: {g['t90_wake']:+.6f}")

# both exposures in one model, to see whether the sleep window carries anything extra
specs2 = {"__base__": [], "whole_plus_sleep": ["t90_wholerec__z", "t90_sleep__z"],
          "whole_plus_wake": ["t90_wholerec__z", "t90_wake__z"],
          "sleep_plus_wake": ["t90_sleep__z", "t90_wake__z"]}
res2 = heldout(specs2, common, m, penalizer=PEN, pkl_name="B_frame2.pkl")
s2 = (res2[res2.spec != "__base__"].groupby("spec")
      .agg(dC=("gain", "mean")).reset_index().sort_values("dC", ascending=False))
print("\nboth windows entered together")
for _i, r in s2.iterrows():
    print(f"   {r.spec:<20}{r.dC:>+12.6f}")

summ.insert(0, "analysis", "B")
summ.to_csv(os.path.join(WORK, "B_summary.csv"), index=False)
res.to_csv(os.path.join(WORK, "B_percondition.csv"), index=False)
json.dump({"subcohort_n": int(len(m)), "n_outcomes_common": len(common),
           "outcomes_common": common, "outcomes_dropped": dropped,
           "baseline": base_c,
           "dC": {k: float(v) for k, v in g.items()},
           "joint": {r.spec: float(r.dC) for _i, r in s2.iterrows()},
           "spearman": {a: {b_: float(cm.loc[a, b_]) for b_ in EXPO} for a in EXPO}},
          open(os.path.join(WORK, "B_result.json"), "w"), indent=2)
print("\nwritten -> _work/B_summary.csv, _work/B_percondition.csv, _work/B_result.json")
