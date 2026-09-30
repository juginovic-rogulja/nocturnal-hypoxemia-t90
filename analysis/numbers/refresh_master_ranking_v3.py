"""
Refresh the ranking keys of MASTER_v2.json from numbers/ranking_v3.csv.

MASTER_v2.json is the hand-assembled headline summary. Its ranking keys were copied out of
ranking_v2.csv, frozen 2026-08-02 on the pre-correction 19,383 cohort, and were never updated
when the oximetry quality rule and the two diagnosis-code fixes landed on 2026-08-04.

Touches only the ranking keys. The cohort, treatment, hidden, vs_sleep and graded blocks are not
read and not written. Nothing is deleted: the superseded values are kept in place under
_SUPERSEDED_ranking_v2_2026-08-02 so a reader can see which way each rank moved.

Companion to refresh_master_cohort.py, which does the same job for the cohort block.
"""
import json
import os

import pandas as pd
import os as _os_v8, sys as _sys_v8; _sys_v8.path.insert(0, _os_v8.path.dirname(_os_v8.path.abspath(__file__)))  # v8 sweep 2026-09-12
from cohort_spec import N_MEASURES  # v8 sweep 2026-09-12

HERE = os.path.dirname(os.path.abspath(__file__))
PATH = f"{HERE}/MASTER_v2.json"

RANK_KEYS = ["rank_spo2_pct_below_90", "rank_AHI", "rank_TST_min", "rank_arousal_index",
             "rank_sleep_efficiency_pct", "rank_N3_pct", "rank_nrem_hr_bpm"]
FEATURE_OF = {k: k[len("rank_"):] for k in RANK_KEYS}


def main():
    rk = pd.read_csv(f"{HERE}/ranking_v3.csv", comment="#")
    assert len(rk) == N_MEASURES, len(rk)
    M = json.load(open(PATH))

    cur = {k: M.get(k, {}).get("rank") for k in RANK_KEYS}
    want = {k: int(rk.loc[rk.feature == FEATURE_OF[k], "rank"].iloc[0]) for k in RANK_KEYS}
    if "_SUPERSEDED_ranking_v2_2026-08-02" in M and cur == want:
        print("already refreshed, nothing to do (v8: the unchanged file is rewritten so the rerun stamper sees a fresh output)")
        json.dump(M, open(PATH, "w"), indent=1); return 0
    if "_SUPERSEDED_ranking_v2_2026-08-02" in M:   # v7: the ranking moved again, keep the first superseded block, add this one
        M[f"_SUPERSEDED_ranking_v3_pre_v7"] = {k: M[k] for k in RANK_KEYS} | {"top5": M["top5"]}

    old = {"note": ("These are the ranking keys as MASTER_v2.json carried them before "
                    "2026-08-07. They come from numbers/ranking_v2.csv, frozen 2026-08-02 on "
                    "the 19,383 cohort, before the oximetry quality rule and before obesity "
                    "hypoventilation and falls were recoded. Kept, not deleted, so the move is "
                    "visible."),
           "top5": M["top5"], "n_measures": M["n_measures"]}
    for k in RANK_KEYS:
        old[k] = M[k]

    top5 = rk.sort_values("rank").head(5)
    M["top5"] = [{"rank": int(r["rank"]), "feature": r.feature, "dC": float(r.dC)}
                 for _i, r in top5.iterrows()]
    M["n_measures"] = int(len(rk))
    for k in RANK_KEYS:
        f = FEATURE_OF[k]
        row = rk[rk.feature == f]
        assert len(row) == 1, f
        M[k] = {"rank": int(row["rank"].iloc[0]), "gain": float(row.dC.iloc[0])}

    M.setdefault("_SUPERSEDED_ranking_v2_2026-08-02", old)
    M["_ranking_source"] = ("numbers/ranking_v3.csv, built 2026-08-07 on the 19,173 cohort of "
                            "cohort_spec.py. Read it with comment='#'.")
    json.dump(M, open(PATH, "w"), indent=1)

    for k in RANK_KEYS:
        print(f"  {k:<32}{old[k]['rank']:>4} -> {M[k]['rank']:<4}"
              f"  {old[k]['gain']:+.6f} -> {M[k]['gain']:+.6f}")
    print(f"\nwrote {PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
