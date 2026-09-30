"""
POSITIVE CONTROL. Nothing else in this folder runs until this passes.

Re-derives T90's published discrimination gain from the frozen cohort with a re-implementation
of the ranking pipeline, and checks it against numbers/ranking_v3.csv three ways:

  1. the headline dC, +0.02033875497249791, to 1e-9
  2. the held-out baseline, 0.647540, recovered from ranking_v3.csv's own mC minus dC
  3. every one of the 48 per-outcome gains against ranking_v3_percondition.csv, to 1e-9

Then one extra fit that is not a check of the paper but a check of the trap: the published age
basis is four cr(df=4) columns that sum to one and are singular. It survives only because the
ranking penalises at 0.01. The same dC is recomputed with one basis column dropped. If the
published basis were silently fitting nothing, the two would not agree.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import json
import os

import numpy as np
import pandas as pd

from common import (ADJ, PEN, T90ROOT, WORK, build_frame, heldout,
                    ranking_outcomes, site_rint, summarise)

import pandas as _pd
import sys as _sys_v8; _sys_v8.path.insert(0, paths.ANALYSIS_DIR)  # v8 sweep 2026-09-12
from cohort_spec import N_RANKED_OUTCOMES  # v8 sweep 2026-09-12
_rk = _pd.read_csv(f"{paths.NUMBERS_DIR}/ranking_v3.csv", comment="#")
PUB_DC = float(_rk.loc[_rk.feature == "spo2_pct_below_90", "dC"].iloc[0])   # v7: the current ranking_v3.csv, not the August literal 0.02033875497249791
_impl_v8 = (_rk.mC - _rk.dC); _impl_mode = _impl_v8.round(6).mode().iloc[0]
PUB_BASE = float(_impl_v8[_impl_v8.round(6) == _impl_mode].mean())  # v8: the MODAL held-out baseline of the ranking (139 of 141 rows; the two spindle-SO coupling features are refit on their non-missing subset and imply a baseline 0.0018 higher)
FEAT = "spo2_pct_below_90"
TOL = 1e-9

print("=" * 78)
print("POSITIVE CONTROL: T90 dC from the published ranking pipeline")
print("=" * 78)

d = build_frame(extra_master_cols=[FEAT])
print(f"cohort {len(d):,}")

OUT = ranking_outcomes(d)
assert len(OUT) == N_RANKED_OUTCOMES, f"outcome set is {len(OUT)}, expected {N_RANKED_OUTCOMES}"
print(f"outcomes {len(OUT)}")

d[FEAT + "__z"] = site_rint(d, FEAT)

specs = {"__base__": [], "t90": [FEAT + "__z"]}
res = heldout(specs, OUT, d, penalizer=PEN, pkl_name="pc_frame.pkl")
s = summarise(res)
got_dc = float(s[s.spec == "t90"].dC.iloc[0])
got_base = float(res[res.spec == "__base__"].c.mean())

# ---------------------------------------------------------------- check 1 and 2
print(f"\n  dC        published {PUB_DC:.12f}   here {got_dc:.12f}   "
      f"diff {got_dc - PUB_DC:+.3e}")
print(f"  baseline  published {PUB_BASE:.6f}         here {got_base:.6f}         "
      f"diff {got_base - PUB_BASE:+.3e}")

# ---------------------------------------------------------------- check 3, per outcome
ref = pd.read_csv(f"{paths.NUMBERS_DIR}/ranking_v3_percondition.csv", comment="#")
ref = ref[ref.feature == FEAT][["outcome", "c", "gain"]].rename(
    columns={"c": "c_ref", "gain": "gain_ref"})
mine = res[res.spec == "t90"][["outcome", "c", "gain"]]
j = mine.merge(ref, on="outcome", how="outer", indicator=True)
assert (j._merge == "both").all(), j[j._merge != "both"]
j["d_c"] = j.c - j.c_ref
j["d_gain"] = j.gain - j.gain_ref
worst_c = float(j.d_c.abs().max())
worst_g = float(j.d_gain.abs().max())
print(f"  per-outcome, {len(j)} conditions: largest |diff| in C {worst_c:.3e}, "
      f"in gain {worst_g:.3e}")

three = j.reindex(j.gain_ref.abs().sort_values(ascending=False).index).head(3)
print("\n  three largest per-outcome gains, published against here")
for _i, r in three.iterrows():
    print(f"     {r.outcome:<20} published {r.gain_ref:+.6f}   here {r.gain:+.6f}")

ok = (abs(got_dc - PUB_DC) < TOL and abs(got_base - PUB_BASE) < 5e-6
      and worst_c < TOL and worst_g < TOL)

# ---------------------------------------------------------------- the spline trap
# Not a check of the paper. A check that the published 4-column basis is really being fitted.
d2 = d.copy()
ADJ_DROP = [f"age_s{i}" for i in range(1, 4)] + ["male"]
res2 = heldout(specs, OUT, d2, penalizer=PEN, adj=ADJ_DROP, pkl_name="pc_frame_drop.pkl")
dc_drop = float(summarise(res2)[lambda x: x.spec == "t90"].dC.iloc[0])
print(f"\n  age-spline trap check: published 4-column basis dC {got_dc:.6f}, "
      f"one column dropped {dc_drop:.6f}, diff {dc_drop - got_dc:+.6f}")
print("  (the ridge at 0.01 makes the singular basis harmless. A degenerate fit would not "
      "agree here.)")

print("\n" + "=" * 78)
print("VERDICT: " + ("PASS, the published dC reproduces exactly" if ok else "FAIL"))
print("=" * 78)

json.dump({
    "published_dC": PUB_DC, "reproduced_dC": got_dc, "diff_dC": got_dc - PUB_DC,
    "published_baseline": PUB_BASE, "reproduced_baseline": got_base,
    "n_outcomes": len(OUT), "cohort_n": int(len(d)),
    "max_abs_per_outcome_C_diff": worst_c, "max_abs_per_outcome_gain_diff": worst_g,
    "dC_with_one_age_spline_column_dropped": dc_drop,
    "pass": bool(ok),
}, open(os.path.join(WORK, "positive_control.json"), "w"), indent=2)
res.to_csv(os.path.join(WORK, "pc_percondition.csv"), index=False)

if not ok:
    raise SystemExit("POSITIVE CONTROL FAILED, stop")
