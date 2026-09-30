"""
E. The ceiling if the ranking were scored only where oxygen plausibly acts.

The published ranking averages each measure's gain over all 48 outcomes, including bipolar
disorder, cataract and hearing loss. This re-averages the SAME frozen per-outcome gains, for
all 197 measures at once, over a restricted outcome universe, and reports where T90 lands.

THE EXACT SUBSET, primary. From ORGAN_GROUP in numbers/disease_definitions.py, the groups
Cardiac, Respiratory, Metabolic and Kidney, plus death. ORGAN_GROUP's vocabulary is "Cardiac"
for cardiovascular and "Kidney" for renal. 21 outcomes.

Two variants are reported beside it rather than instead of it, because two defensible readings
of the same instruction give different sets:
  extended   primary plus the cardiovascular composite, which ORGAN_GROUP files under
             "Composite", and plus stroke, which ORGAN_GROUP files under "Neuro/psych"
             although it is a cardiovascular event. 23 outcomes.
  plus-liver extended plus the Liver group, fatty liver disease and cirrhosis. Not requested,
             reported because it is the single largest omission from a hypoxia-pathway reading
             and leaving it out silently would flatter the primary number by comparison.
             25 outcomes.

No refit. This is arithmetic on numbers/ranking_v3_percondition.csv, which s00 reproduced.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import os

import numpy as np
import pandas as pd

from common import DISEASES, ORGAN_GROUP, OXYGEN_PATHWAY_GROUPS, T90ROOT, WORK, organ_of
import sys as _sys_v8; _sys_v8.path.insert(0, paths.ANALYSIS_DIR)  # v8 sweep 2026-09-12
from cohort_spec import N_RANKED_OUTCOMES  # v8 sweep 2026-09-12

FEAT = "spo2_pct_below_90"
P = pd.read_csv(f"{paths.NUMBERS_DIR}/ranking_v3_percondition.csv", comment="#")
RK = pd.read_csv(f"{paths.NUMBERS_DIR}/ranking_v3.csv", comment="#")
FAM = pd.read_csv(f"{paths.NUMBERS_DIR}/measure_families.csv")
dup = set(FAM[FAM.is_siteZ_duplicate.astype(bool)].feature)

ALL48 = sorted(P.outcome.unique())
assert len(ALL48) == N_RANKED_OUTCOMES, len(ALL48)

PRIMARY = sorted([o for o in ALL48 if organ_of(o) in OXYGEN_PATHWAY_GROUPS] + ["death"])
EXTENDED = sorted(set(PRIMARY) | {"cvd", "stroke_any"})
PLUSLIVER = sorted(set(EXTENDED) | {o for o in ALL48 if organ_of(o) == "Liver"})

print("=" * 78)
print("E. RESTRICTED RANKING CEILING")
print("=" * 78)
for nm, st in (("primary", PRIMARY), ("extended", EXTENDED), ("plus-liver", PLUSLIVER)):
    print(f"\n{nm}: {len(st)} outcomes")
    print("   " + ", ".join(
        f"{DISEASES[o][0] if o in DISEASES else 'Death'}" for o in st))

rows = []
for nm, st in (("published_all48", ALL48), ("primary", PRIMARY),
               ("extended", EXTENDED), ("plus_liver", PLUSLIVER)):
    sub = P[P.outcome.isin(st)]
    assert sub.groupby("feature").size().nunique() == 1, "ragged outcome coverage"
    r = (sub.groupby("feature")
            .agg(dC=("gain", "mean"), mC=("c", "mean"),
                 npos=("gain", lambda x: int((x > 0).sum())), nout=("gain", "size"))
            .reset_index().sort_values("dC", ascending=False))
    r["rank"] = range(1, len(r) + 1)
    nd = r[~r.feature.isin(dup)].copy()
    nd["rank_nodup"] = range(1, len(nd) + 1)
    r = r.merge(nd[["feature", "rank_nodup"]], on="feature", how="left")
    t = r[r.feature == FEAT].iloc[0]
    rows.append({"subset": nm, "n_outcomes": len(st), "n_measures": len(r),
                 "t90_dC": float(t.dC), "t90_rank": int(t["rank"]),
                 "t90_rank_nodup": int(t.rank_nodup),
                 "t90_npos": int(t.npos),
                 "runner_up": r.iloc[1].feature, "runner_up_dC": float(r.iloc[1].dC),
                 "lead_over_runner_up": float(t.dC - r.iloc[1].dC),
                 "top5": " | ".join(f"{a} {b:+.5f}" for a, b in
                                    zip(r.feature.head(5), r.dC.head(5)))})
    r.insert(0, "subset", nm)
    r.to_csv(os.path.join(WORK, f"E_ranking_{nm}.csv"), index=False)

S = pd.DataFrame(rows)
print("\n" + "=" * 78)
print(f"{'subset':<18}{'nOut':>5}{'T90 dC':>12}{'rank':>6}{'rank*':>7}{'n>0':>5}"
      f"{'lead over #2':>14}")
print("-" * 78)
for _i, r in S.iterrows():
    print(f"{r.subset:<18}{r.n_outcomes:>5}{r.t90_dC:>+12.6f}{r.t90_rank:>6}"
          f"{r.t90_rank_nodup:>7}{r.t90_npos:>5}{r.lead_over_runner_up:>+14.6f}")
print("rank* is the rank after dropping the site-z duplicate rows of measure_families.csv")

pub = S[S.subset == "published_all48"].iloc[0]
pri = S[S.subset == "primary"].iloc[0]
print(f"\nT90 goes from {pub.t90_dC:+.6f} over {N_RANKED_OUTCOMES} outcomes to {pri.t90_dC:+.6f} over the "
      f"{pri.n_outcomes} oxygen-pathway outcomes, a {pri.t90_dC / pub.t90_dC:.2f}-fold rise.")
print(f"It is rank {pri.t90_rank} of {pri.n_measures} in the restricted universe against "
      f"rank {pub.t90_rank} of {pub.n_measures} published.")

print("\ntop 5 in each universe")
for _i, r in S.iterrows():
    print(f"   {r.subset:<18}{r.top5}")

# the noise floor, for scale
import json
nf = json.load(open(f"{paths.NUMBERS_DIR}/ranking_v3_noisefloor.json"))
print(f"\nnoise floor for scale: a meaningless variable earns "
      f"{nf['dC_per_noise_variable_mean']:+.6f} on the published {N_RANKED_OUTCOMES}, 95% of them between "
      f"{nf['dC_95pct_range_of_a_useless_variable'][0]:+.6f} and "
      f"{nf['dC_95pct_range_of_a_useless_variable'][1]:+.6f}")

S.to_csv(os.path.join(WORK, "E_summary.csv"), index=False)
json.dump({"primary_subset": PRIMARY, "extended_subset": EXTENDED,
           "plus_liver_subset": PLUSLIVER,
           "organ_groups_used": OXYGEN_PATHWAY_GROUPS},
          open(os.path.join(WORK, "E_subsets.json"), "w"), indent=2)
print("\nwritten -> _work/E_summary.csv, _work/E_ranking_*.csv, _work/E_subsets.json")
