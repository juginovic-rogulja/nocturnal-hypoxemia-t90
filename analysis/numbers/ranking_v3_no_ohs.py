"""
Step 232 ranking_without_ohs (decision 4): the measure ranking with obesity hypoventilation
syndrome left out of the outcome set. Arithmetic only, on ranking_v3_percondition.csv: per measure
the mean held-out gain over the remaining outcomes, re-ranked. The with-OHS mean must equal
ranking_v3.csv's dC to 1e-9 first (the same file, the same arithmetic), so the only difference
between the two rankings is the one outcome removed.

Refreshes the Limitations sentence: "T90 remained the top-ranked measurement when this condition
was excluded" is printed as the rank T90 actually takes.

Outputs (cohort_spec.OUT_DIR): ranking_v3_no_ohs.csv, ranking_v3_no_ohs.json.
"""
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from cohort_spec import NUMBERS_DIR, OUT_DIR, PRIMARY_EXPOSURE, sidecar  # noqa: E402
from disease_definitions import DISEASES  # noqa: E402

OHS_KEY = "obesity_hypovent"
assert OHS_KEY in DISEASES, f"{OHS_KEY} is not a condition key"
PC, RK = f"{NUMBERS_DIR}/ranking_v3_percondition.csv", f"{NUMBERS_DIR}/ranking_v3.csv"
per = pd.read_csv(PC, comment="#")
rk = pd.read_csv(RK, comment="#")
assert OHS_KEY in set(per.outcome), f"{OHS_KEY} is not among the ranked outcomes of {PC}"

full = per.groupby("feature").gain.mean()
chk = rk.set_index("feature").dC
diff = float((full.reindex(chk.index) - chk).abs().max())
assert diff < 1e-9, f"the with-OHS mean recomputed from the per-outcome file differs from ranking_v3.csv by {diff:.2e}"

sub = per[per.outcome != OHS_KEY]
r2 = (sub.groupby("feature").agg(mC=("c", "mean"), dC=("gain", "mean"), worst=("gain", "min"),
                                 npos=("gain", lambda x: int((x > 0).sum())), nout=("gain", "size"))
      .reset_index().sort_values("dC", ascending=False))
r2.insert(0, "rank", range(1, len(r2) + 1))
r2 = r2.merge(rk[["feature", "rank", "dC"]].rename(columns={"rank": "rank_with_ohs", "dC": "dC_with_ohs"}), on="feature", how="left")
r2["rank_change"] = r2["rank"] - r2["rank_with_ohs"]

t = r2[r2.feature == PRIMARY_EXPOSURE].iloc[0]
top10_with = rk.sort_values("rank").head(10).feature.tolist()
top10_without = r2.sort_values("rank").head(10).feature.tolist()
out = {"outcome_excluded": OHS_KEY, "outcome_excluded_label": DISEASES[OHS_KEY][0],
       "n_outcomes_with": int(per.outcome.nunique()), "n_outcomes_without": int(sub.outcome.nunique()),
       "n_measures": int(len(r2)),
       "t90": {"rank_with_ohs": int(t.rank_with_ohs), "rank_without_ohs": int(t["rank"]),
               "dC_with_ohs": float(t.dC_with_ohs), "dC_without_ohs": float(t.dC),
               "top_ranked_without_ohs": bool(int(t["rank"]) == 1)},
       "top10_with_ohs": top10_with, "top10_without_ohs": top10_without,
       "n_top10_changed": len(set(top10_with) ^ set(top10_without)) // 2,
       "max_abs_rank_change": int(r2.rank_change.abs().max()),
       "with_ohs_mean_equals_ranking_v3_to": diff,
       "sentence": (f"T90 {'remained' if int(t['rank']) == 1 else 'did not remain'} the top-ranked measurement when this condition was excluded"
                    f" (rank {int(t['rank'])} of {len(r2)} without obesity hypoventilation, {int(t.rank_with_ohs)} with)")}
# v8.1 (14 Sept): the same arithmetic without respiratory failure (its definition rests on oxygen saturation too) and without both
RF_KEY = "resp_failure"
assert RF_KEY in DISEASES and RF_KEY in set(per.outcome), f"{RF_KEY} is not among the ranked outcomes"


def rank_without(keys):
    s2 = per[~per.outcome.isin(keys)]
    r = (s2.groupby("feature").agg(dC=("gain", "mean")).reset_index().sort_values("dC", ascending=False))
    r.insert(0, "rank", range(1, len(r) + 1))
    tt = r[r.feature == PRIMARY_EXPOSURE].iloc[0]
    return r, {"outcomes_excluded": list(keys), "n_outcomes_without": int(s2.outcome.nunique()), "t90_rank_without": int(tt["rank"]),
               "t90_dC_without": float(tt.dC), "top10_without": r.head(10).feature.tolist(), "top_ranked_without": bool(int(tt["rank"]) == 1)}


r_rf, without_rf = rank_without([RF_KEY])
r_both, without_both = rank_without([OHS_KEY, RF_KEY])
r2 = r2.merge(r_rf[["feature", "rank"]].rename(columns={"rank": "rank_without_resp_failure"}), on="feature", how="left")
r2 = r2.merge(r_both[["feature", "rank"]].rename(columns={"rank": "rank_without_ohs_and_resp_failure"}), on="feature", how="left")
out["without_resp_failure"] = without_rf
out["without_ohs_and_resp_failure"] = without_both
out["sentence_both"] = (f"T90 {'remained' if without_both['top_ranked_without'] else 'did not remain'} the top-ranked measurement when obesity hypoventilation "
                        f"and respiratory failure were both excluded (rank {without_both['t90_rank_without']} of {len(r2)}, gain {without_both['t90_dC_without']:.4f} against "
                        f"{float(t.dC_with_ohs):.4f} with them)")
csv_path, json_path = f"{OUT_DIR}/ranking_v3_no_ohs.csv", f"{OUT_DIR}/ranking_v3_no_ohs.json"
r2.to_csv(csv_path, index=False)
json.dump(out, open(json_path, "w"), indent=1)
for p in (csv_path, json_path):
    sidecar(p, __file__, extra_inputs=[PC, RK], note=out["sentence"])
print(out["sentence"]); print(out["sentence_both"])
print(f"top 10 changed by {out['n_top10_changed']} measure(s); largest rank move {out['max_abs_rank_change']}; written {csv_path}, {json_path}")
