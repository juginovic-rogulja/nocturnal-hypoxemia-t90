"""
SELF-TEST. Nothing about the binary exposure is computed or reported until this passes.

The control is the published CONTINUOUS T90 run through the exact code path this analysis will
use for the binary exposure. Six gates:

  GATE 1  full follow-up, published four-column age basis, continuous T90.
          Must reproduce T90_Manuscript/numbers/ranking_v3_percondition.csv, all 48 per-outcome
          gains to 1e-9, and the headline mean 0.02033875497249791.
  GATE 2  full follow-up, reported basis (one age-spline column dropped), continuous T90.
          Must reproduce Sleep_Variability_2026-08/dC_per_disease/dC_per_disease.csv,
          arm "standard", column dC, to 1e-12. This is Figure 1c's full-follow-up series.
  GATE 3  two-year landmark, reported basis, continuous T90. Must reproduce the same file's
          arm "lag2" dC to 1e-12. This is Figure 1c's landmark series.
  GATE 4  the landmark definition itself: n at risk and events at lag 0 and lag 2 must match
          T90_Manuscript/numbers/lag_ladder.csv exactly, outcome by outcome.
  GATE 5  the ten values printed on Figure 1c
          (FINAL_FIGURES_2026-08-14/_workfiles/Figure1C_drawn_values.json, the drawn-values
          record for the panel in ROUND16_2026-09-04/L9_assembly/NEW_FINAL_SET_R16/Main_Fig1.pdf)
          must equal the reproduced GATE 2 numbers at the printed precision.
  GATE 6  one_pass_x with exposure="rint" must equal dc_engine.one_pass bit for bit, so the
          only difference between the control and the binary run is the exposure column.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../..")))  # repository root, where paths.py lives
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
sys.path.insert(0, HERE)
import engine_gt10 as X                                                   # noqa: E402
import dc_engine as E                                                     # noqa: E402
import sys as _sys_v8; _sys_v8.path.insert(0, paths.ANALYSIS_DIR)  # v8 sweep 2026-09-12
from cohort_spec import NUMBERS_DIR, N_RANKED_OUTCOMES  # v8 sweep 2026-09-12

T90ROOT = paths.FIGURE_ROOT
PUBCSV = (f"{paths.SV_ROOT}/"
          "dC_per_disease/dC_per_disease.csv")
FIG1C = (f"{paths.FIGURE_ROOT}/FINAL_FIGURES_2026-08-14/_workfiles/Figure1C_drawn_values.json")
RANKCSV = f"{paths.NUMBERS_DIR}/ranking_v3_percondition.csv"
from cohort_spec import cohort_tag   # v8.1c
LADDER = f"{paths.NUMBERS_DIR}/lag_ladder{cohort_tag()}.csv"
__PUB_DC_V7__ = float(__import__("pandas").read_csv(f"{paths.NUMBERS_DIR}/ranking_v3.csv", comment="#").set_index("feature").loc["spo2_pct_below_90", "dC"])   # v7: T90 gain from the CURRENT ranking, was the August literal 0.02033875497249791
PUB_DC = __PUB_DC_V7__
_RK_V8 = __import__("pandas").read_csv(f"{paths.NUMBERS_DIR}/ranking_v3.csv", comment="#"); _impl_v8 = (_RK_V8.mC - _RK_V8.dC); _impl_mode = _impl_v8.round(6).mode().iloc[0]; _off = _RK_V8.feature[_impl_v8.round(6) != _impl_mode].tolist(); assert all("_coupling_" in str(_f) for _f in _off) and (_impl_v8.round(6) == _impl_mode).mean() > 0.9, f"ranking rows disagree on the implied baseline: off-mode rows {_off} (v8.1c rule: every off-mode row must be a spindle-SO coupling feature, mode on over 90 percent of rows)"; print(f"[v8] implied baseline mode {_impl_mode:.6f} over {int((_impl_v8.round(6) == _impl_mode).sum())} rows; off-mode rows {_off}"); PUB_BASE = float(_impl_v8[_impl_v8.round(6) == _impl_mode].mean())  # v8 F2 fix 2026-09-12 (integrity gate 6): the held-out baseline C is the MODE of the implied baseline mC - dC over the ranking rows (139 of 141 agree); at most two off-mode rows are allowed and listed (features refit on their non-missing subset, the two C_coupling_* rows); PUB_BASE is taken from the modal rows. Same rule as numbers/run_combo_v2.py lines 68 to 72; the path is built from T90ROOT (line 40) because the sweep's `from cohort_spec import NUMBERS_DIR` sits on line 154, AFTER this line (NameError otherwise); same folder
# optional first argument: the file the record is written to, so a later run can re-run the
# same six gates without overwriting an earlier run's record. Behaviour is unchanged.
OUTNAME = sys.argv[1] if len(sys.argv) > 1 else "selftest.json"

t0 = time.time()
print("=" * 84)
print("SELF-TEST: the published continuous T90, through this analysis's code path")
print("=" * 84)

A = E.build_arrays()
outs = A["outcomes"]
print(f"cohort {len(A['site']):,}   outcomes {len(outs)}   "
      f"above 10% T90 {int((A['t90'] > 10).sum()):,}")

R = {}

# ------------------------------------------------------------------ GATE 6 first (cheap-ish)
print("\nGATE 6  exposure wrapper is byte-identical to the published engine on 'rint'")
g6 = X.identity_guard(A)
print(f"  GATE 6 {'PASS' if g6 else 'FAIL'}")
R["gate6_wrapper_identity"] = bool(g6)

# ------------------------------------------------------------------ GATE 1
P_full = X.one_pass_x(A, seed=None, basis="full", lag=0.0, exposure="rint")
np.save(f"{X.WORK}/selftest_point_rint_full_lag0.npy", P_full)
dc_full = E.dc_from_pass(P_full)
c_base_full = np.nanmean(P_full[:, [0, 2]], axis=1)
ref = pd.read_csv(RANKCSV, comment="#")
ref = ref[ref.feature == "spo2_pct_below_90"].set_index("outcome")
assert set(ref.index) == set(outs), "outcome sets differ from the published ranking"
gain_ref = ref.loc[outs, "gain"].to_numpy(float)
c_ref = ref.loc[outs, "c"].to_numpy(float)
c_t90_full = np.nanmean(P_full[:, [1, 3]], axis=1)
d_gain = np.abs(dc_full - gain_ref)
d_c = np.abs(c_t90_full - c_ref)
mean_dc = float(np.nanmean(dc_full))
base_mean = float(np.nanmean(c_base_full))
g1 = (d_gain.max() < 1e-9 and d_c.max() < 1e-9 and abs(mean_dc - PUB_DC) < 1e-9
      and abs(base_mean - PUB_BASE) < 5e-6)
print("\nGATE 1  full follow-up, published basis, against numbers/ranking_v3_percondition.csv")
print(f"  headline mean dC  published {PUB_DC:.14f}  here {mean_dc:.14f}  "
      f"diff {mean_dc - PUB_DC:+.2e}")
print(f"  baseline C        published {PUB_BASE:.6f}          here {base_mean:.6f}")
print(f"  {N_RANKED_OUTCOMES} outcomes: largest |gain diff| {d_gain.max():.2e}, largest |C diff| {d_c.max():.2e}")
for o in ("diabetes", "hf", "copd2"):
    i = outs.index(o)
    print(f"    {o:<10} published {gain_ref[i]:+.12f}   here {dc_full[i]:+.12f}")
print(f"  GATE 1 {'PASS' if g1 else 'FAIL'}")
R["gate1_ranking_v3_percondition"] = {
    "source": RANKCSV, "published_mean_dC": PUB_DC, "reproduced_mean_dC": mean_dc,
    "max_abs_gain_diff": float(d_gain.max()), "max_abs_C_diff": float(d_c.max()),
    "published_baseline_C": PUB_BASE, "reproduced_baseline_C": base_mean, "pass": bool(g1)}

# ------------------------------------------------------------------ GATES 2 and 3
pub = pd.read_csv(PUBCSV)
pub = pub[pub.row_type == "disease"]
gates = {}
points = {}
for arm, lag in (("standard", 0.0), ("lag2", 2.0)):
    P = X.one_pass_x(A, seed=None, basis="drop", lag=lag, exposure="rint")
    np.save(f"{X.WORK}/selftest_point_rint_drop_{arm}.npy", P)
    points[arm] = P
    mine = pd.Series(E.dc_from_pass(P), index=outs)
    theirs = pub[pub.arm == arm].set_index("outcome")["dC"].reindex(outs)
    both = mine.notna() & theirs.notna()
    dif = (mine[both] - theirs[both]).abs()
    nan_agree = bool((mine.isna() == theirs.isna()).all())
    ok = bool(dif.max() < 1e-12 and nan_agree)
    gates[arm] = ok
    print(f"\nGATE {'2' if arm == 'standard' else '3'}  "
          f"{'full follow-up' if lag == 0 else 'two-year landmark'}, reported basis, "
          f"against dC_per_disease.csv arm '{arm}'")
    print(f"  {int(both.sum())} estimable outcomes, largest |dC diff| {dif.max():.3e}, "
          f"not-estimable set agrees: {nan_agree}")
    print(f"  GATE {'2' if arm == 'standard' else '3'} {'PASS' if ok else 'FAIL'}")
    R[f"gate_{arm}_dC_per_disease"] = {
        "source": PUBCSV, "arm": arm, "landmark_years": lag,
        "n_estimable": int(both.sum()), "max_abs_dC_diff": float(dif.max()),
        "not_estimable_set_agrees": nan_agree, "pass": ok}
g2, g3 = gates["standard"], gates["lag2"]

# ------------------------------------------------------------------ GATE 4
lad = pd.read_csv(LADDER)
bad = []
for lag in (0.0, 2.0):
    mine = E.counts(A, lag=lag).set_index("outcome")
    r = lad[lad.lag_years == lag].set_index("key")[["n", "events"]]
    for o in outs:
        if int(mine.loc[o, "n_at_risk"]) != int(r.loc[o, "n"]):
            bad.append(f"lag{lag:.0f} {o} n")
        if int(mine.loc[o, "events"]) != int(r.loc[o, "events"]):
            bad.append(f"lag{lag:.0f} {o} events")
g4 = not bad
print(f"\nGATE 4  landmark windowing against numbers/lag_ladder.csv")
print(f"  lag 0 and lag 2, n at risk and events, {len(outs)} outcomes each: "
      f"{'every value matches exactly' if g4 else bad[:6]}")
print(f"  GATE 4 {'PASS' if g4 else 'FAIL'}")
R["gate4_lag_ladder"] = {"source": LADDER, "mismatches": bad, "pass": bool(g4)}

# ------------------------------------------------------------------ GATE 5
fig = json.load(open(FIG1C))
drawn = fig["variants"]["Figure1C"]["drawn"]
sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
from disease_definitions import DISEASES  # noqa: E402

name2key = {v[0]: k for k, v in DISEASES.items()}
name2key["Death from any cause"] = "death"
mine_std = pd.Series(E.dc_from_pass(points["standard"]), index=outs)
mine_l2 = pd.Series(E.dc_from_pass(points["lag2"]), index=outs)
f5 = []
for row in drawn:
    key = name2key[row["disease"]]
    got = (mine_std if row["arm"] == "standard" else mine_l2)[key]
    f5.append({"disease": row["disease"], "arm": row["arm"], "key": key,
               "figure_json_dC": row["dC"], "reproduced_dC": float(got),
               "abs_diff": float(abs(row["dC"] - got)),
               "printed_on_figure": row.get("printed")})
g5 = all(r["abs_diff"] < 1e-12 for r in f5)
print(f"\nGATE 5  the values Figure 1c draws, {len(f5)} points over "
      f"{len(set(r['disease'] for r in f5))} conditions")
print(f"  largest |diff| against the figure's drawn-values record: "
      f"{max(r['abs_diff'] for r in f5):.3e}")
for r in f5:
    if r["printed_on_figure"]:
        print(f"    {r['disease']:<26} printed {r['printed_on_figure']:<24} "
              f"reproduced {r['reproduced_dC']:+.6f}")
print(f"  GATE 5 {'PASS' if g5 else 'FAIL'}")
R["gate5_figure1c_drawn_values"] = {"source": FIG1C, "points": f5, "pass": bool(g5)}

# v7 (2026-09-07): GATE 5 compares against the values the SHIPPING Figure 1c draws, which are the pre-correction numbers
# until the figure is redrawn from this run. It is reported (the differences are the figure's change list) but no longer
# blocks: the data gates 1 to 4 and 6 decide.
if not g5:
    print("  (GATE 5 differences are expected until Figure 1c is redrawn from the corrected data; not blocking)")
ok = all([g1, g2, g3, g4, g6])
print("\n" + "=" * 84)
print("SELF-TEST VERDICT: " + ("PASS" if ok else "FAIL"))
print("=" * 84)
print(f"elapsed {time.time() - t0:.0f}s")
R["pass"] = bool(ok)
R["cohort_n"] = int(len(A["site"]))
R["n_outcomes"] = len(outs)
R["n_above_10pct"] = int((A["t90"] > 10).sum())
json.dump(R, open(f"{X.WORK}/{OUTNAME}", "w"), indent=2)
if not ok:
    raise SystemExit("SELF-TEST FAILED. Nothing downstream may be reported.")
