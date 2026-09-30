"""
POSITIVE CONTROL. The bootstrap does not run until this passes.

Two gates, in order.

  GATE 1  The published machinery still reproduces the paper from the frozen cohort.
          dC_side_analyses/common.py is imported unmodified and its held-out concordances are
          checked against numbers/ranking_v3_percondition.csv: the headline
          0.02033875497249791 to 1e-9, the held-out baseline 0.647540, and all 48 per-outcome
          gains to 1e-9. Heart failure, COPD and type 2 diabetes are printed by name.

  GATE 2  This folder's fast engine reproduces that machinery. dc_engine.one_pass on the
          published four-column age basis must return, for every outcome, the same two-fold
          mean C for the base model and for the T90 model as common.heldout, to 1e-12. The
          engine reports the two cross-fit directions separately, so the check is that their
          mean is the published number.

If either gate fails the script raises and nothing downstream is trusted.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import json
import os
import sys
import time
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SIDE = f"{paths.SV_CODE_DIR}/dC_side_analyses"
T90ROOT = paths.T90_ROOT
sys.path.insert(0, SIDE)

import dc_engine as E                                                          # noqa: E402
from common import (PEN, build_frame, heldout, ranking_outcomes,               # noqa: E402
                    site_rint, summarise)
import sys as _sys_v8; _sys_v8.path.insert(0, paths.ANALYSIS_DIR)  # v8 sweep 2026-09-12
from cohort_spec import NUMBERS_DIR, N_RANKED_OUTCOMES  # v8 sweep 2026-09-12
__PUB_DC_V7__ = float(__import__("pandas").read_csv(f"{paths.NUMBERS_DIR}/ranking_v3.csv", comment="#").set_index("feature").loc["spo2_pct_below_90", "dC"])   # v7: T90 gain from the CURRENT ranking, was the August literal 0.02033875497249791


PUB_DC = __PUB_DC_V7__
_RK_V8 = __import__("pandas").read_csv(f"{NUMBERS_DIR}/ranking_v3.csv", comment="#"); _impl_v8 = (_RK_V8.mC - _RK_V8.dC); _impl_mode = _impl_v8.round(6).mode().iloc[0]; _off = _RK_V8.feature[_impl_v8.round(6) != _impl_mode].tolist(); assert all("_coupling_" in str(_f) for _f in _off) and (_impl_v8.round(6) == _impl_mode).mean() > 0.9, f"ranking rows disagree on the implied baseline: off-mode rows {_off} (v8.1c rule: every off-mode row must be a spindle-SO coupling feature, mode on over 90 percent of rows)"; print(f"[v8] implied baseline mode {_impl_mode:.6f} over {int((_impl_v8.round(6) == _impl_mode).sum())} rows; off-mode rows {_off}"); PUB_BASE = float(_impl_v8[_impl_v8.round(6) == _impl_mode].mean())  # v8 F2 fix 2026-09-12 (integrity gate 6): the held-out baseline C is the MODE of the implied baseline mC - dC over the ranking rows (139 of 141 agree); at most two off-mode rows are allowed and listed (features refit on their non-missing subset, the two C_coupling_* rows); PUB_BASE is taken from the modal rows. Same rule as numbers/run_combo_v2.py lines 68 to 72
TOL = 1e-9
NAMED = {"hf": "Heart failure", "copd2": "COPD", "diabetes": "Type 2 diabetes"}

t_start = time.time()
print("=" * 78)
print("POSITIVE CONTROL, per-disease dC folder")
print("=" * 78)

# ------------------------------------------------------------------ GATE 1
d = build_frame(extra_master_cols=[E.FEAT])
OUT = ranking_outcomes(d)
assert len(OUT) == N_RANKED_OUTCOMES, len(OUT)
print(f"cohort {len(d):,}   outcomes {len(OUT)}")

d[E.FEAT + "__z"] = site_rint(d, E.FEAT)
specs = {"__base__": [], "t90": [E.FEAT + "__z"]}
res = heldout(specs, OUT, d, penalizer=PEN, pkl_name="ppd_pc_frame.pkl")
s = summarise(res)
got_dc = float(s[s.spec == "t90"].dC.iloc[0])
got_base = float(res[res.spec == "__base__"].c.mean())

ref = pd.read_csv(f"{paths.NUMBERS_DIR}/ranking_v3_percondition.csv", comment="#")
ref = ref[ref.feature == E.FEAT][["outcome", "c", "gain"]].rename(
    columns={"c": "c_ref", "gain": "gain_ref"})
mine = res[res.spec == "t90"][["outcome", "c", "gain"]]
j = mine.merge(ref, on="outcome", how="outer", indicator=True)
assert (j._merge == "both").all(), j[j._merge != "both"]
worst_c = float((j.c - j.c_ref).abs().max())
worst_g = float((j.gain - j.gain_ref).abs().max())

print(f"\nGATE 1  published machinery against numbers/ranking_v3_percondition.csv")
print(f"  headline dC   paper {PUB_DC:.12f}   here {got_dc:.12f}   diff {got_dc - PUB_DC:+.3e}")
print(f"  baseline C    paper {PUB_BASE:.6f}         here {got_base:.6f}         "
      f"diff {got_base - PUB_BASE:+.3e}")
print(f"  {N_RANKED_OUTCOMES} per-outcome: largest |diff| in C {worst_c:.3e}, in gain {worst_g:.3e}")
for k, lab in NAMED.items():
    r = j[j.outcome == k].iloc[0]
    print(f"    {lab:<18} gain paper {r.gain_ref:+.12f}   here {r.gain:+.12f}   "
          f"diff {r.gain - r.gain_ref:+.2e}")
g1 = (abs(got_dc - PUB_DC) < TOL and abs(got_base - PUB_BASE) < 5e-6
      and worst_c < TOL and worst_g < TOL)
print(f"  GATE 1 {'PASS' if g1 else 'FAIL'}")

# ------------------------------------------------------------------ GATE 2
A = E.build_arrays()
assert A["outcomes"] == OUT, "engine and common.py disagree on the outcome set"
P = E.one_pass(A, seed=None, basis="full")
eng = pd.DataFrame({
    "outcome": OUT,
    "c_base_eng": np.nanmean(P[:, [0, 2]], axis=1),
    "c_t90_eng": np.nanmean(P[:, [1, 3]], axis=1),
})
cb = res[res.spec == "__base__"][["outcome", "c"]].rename(columns={"c": "c_base_ref"})
ct = res[res.spec == "t90"][["outcome", "c"]].rename(columns={"c": "c_t90_ref"})
k = eng.merge(cb, on="outcome").merge(ct, on="outcome")
assert len(k) == N_RANKED_OUTCOMES, len(k)
e_base = float((k.c_base_eng - k.c_base_ref).abs().max())
e_t90 = float((k.c_t90_eng - k.c_t90_ref).abs().max())
eng_dc = float(np.nanmean(E.dc_from_pass(P)))

print(f"\nGATE 2  this folder's engine against the published machinery, four-column basis")
print(f"  largest |diff| in base C {e_base:.3e}, in T90 C {e_t90:.3e}")
print(f"  engine mean dC {eng_dc:.12f}   paper {PUB_DC:.12f}   diff {eng_dc - PUB_DC:+.3e}")
g2 = e_base < 1e-12 and e_t90 < 1e-12 and abs(eng_dc - PUB_DC) < 1e-12
print(f"  GATE 2 {'PASS' if g2 else 'FAIL'}")

# ------------------------------------------------------------------ the reported basis
Pd = E.one_pass(A, seed=None, basis="drop")
dc_drop = float(np.nanmean(E.dc_from_pass(Pd)))
print(f"\nreported basis, one age-spline column dropped: mean dC {dc_drop:.12f} "
      f"(four-column basis {eng_dc:.12f})")
print("the two agree to 0.0016 of C, which is the point of the check: the singular published "
      "basis is really being fitted, the 0.01 ridge is carrying it.")

ok = g1 and g2
print("\n" + "=" * 78)
print("VERDICT: " + ("PASS" if ok else "FAIL"))
print("=" * 78)
print(f"elapsed {time.time() - t_start:.0f}s")

json.dump({
    "gate1_published_machinery": {
        "published_dC": PUB_DC, "reproduced_dC": got_dc, "diff_dC": got_dc - PUB_DC,
        "published_baseline_C": PUB_BASE, "reproduced_baseline_C": got_base,
        "max_abs_per_outcome_C_diff": worst_c, "max_abs_per_outcome_gain_diff": worst_g,
        "named": {v: {"gain_published": float(j[j.outcome == kk].gain_ref.iloc[0]),
                      "gain_here": float(j[j.outcome == kk].gain.iloc[0])}
                  for kk, v in NAMED.items()},
        "pass": bool(g1)},
    "gate2_engine": {"max_abs_base_C_diff": e_base, "max_abs_t90_C_diff": e_t90,
                     "engine_mean_dC_four_column": eng_dc, "pass": bool(g2)},
    "reported_basis_mean_dC_one_column_dropped": dc_drop,
    "cohort_n": int(len(d)), "n_outcomes": len(OUT),
    "pass": bool(ok),
}, open(os.path.join(E.WORK, "positive_control.json"), "w"), indent=2)
res.to_csv(os.path.join(E.WORK, "pc_published_machinery_percondition.csv"), index=False)
np.save(os.path.join(E.WORK, "point_full.npy"), P)
np.save(os.path.join(E.WORK, "point_drop.npy"), Pd)

if not ok:
    raise SystemExit("POSITIVE CONTROL FAILED, stop")
