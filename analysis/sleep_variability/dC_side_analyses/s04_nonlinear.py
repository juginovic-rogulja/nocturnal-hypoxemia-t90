"""
D. Does the paper's transform of T90 cost it anything?

The published exposure is the within-site rank inverse normal of T90. That throws away the
shape of the exposure entirely and keeps only the ordering within each hospital. Four
alternatives are fitted against it, everything else held fixed: same cohort, same 48 outcomes,
same age-and-sex base, same folds, same penalizer.

  rint_site      the published transform, the reference line
  raw_linear     T90 as recorded, untransformed, one column
  log1p_raw      log(1 + T90), untransformed afterwards
  log1p_sitez    log(1 + T90) then standardised within hospital, so it keeps the published
                 transform's site alignment and changes only its shape
  rcs3_raw       restricted cubic spline of raw T90, 3 knots at the 10th, 50th and 90th
                 percentiles, Harrell's placement, 2 basis columns after dropping one
  rcs5_raw       the same with 5 knots, 3 interior at the 10th, 50th and 90th percentiles
                 plus the two boundaries, 4 basis columns. Fitted because if curvature is
                 what the published transform is missing, more of it should help.
  rcs3_on_rint   the 3-knot spline applied to the published rank-inverse-normal score, which
                 asks whether curvature adds anything once the site alignment is already there
  bands4         the paper's four oxygen bands, 0-1, 1-5, 5-10 and over 10 percent, as 3
                 indicator columns with 0-1 percent as the reference

THE SPLINE RANK TRAP. Both spline bases here sum to one and are singular against each other and
against the age basis. One column is dropped from every T90 spline basis, per the standing rule.
The age basis stays at the published four columns because the positive control gates on
reproducing the published number, and s00 proved that basis is not fitting nothing.

A multi-column model is not free: three indicator columns get three chances to fit noise where
one score gets one. The noise floor of numbers/ranking_v3_noisefloor.json is printed for scale.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import json
import os

import numpy as np
import pandas as pd
from patsy import dmatrix

from common import (BANDS, PEN, T90ROOT, WORK, band_index, build_frame, heldout,
                    ranking_outcomes, site_rint, site_z, summarise)
import sys as _sys_v8; _sys_v8.path.insert(0, paths.ANALYSIS_DIR)  # v8 sweep 2026-09-12
from cohort_spec import N_RANKED_OUTCOMES  # v8 sweep 2026-09-12
__PUB_DC_V7__ = float(__import__("pandas").read_csv(f"{paths.NUMBERS_DIR}/ranking_v3.csv", comment="#").set_index("feature").loc["spo2_pct_below_90", "dC"])   # v7: T90 gain from the CURRENT ranking, was the August literal 0.02033875497249791


T90 = "spo2_pct_below_90"
print("=" * 78)
print("D. NON-LINEAR AND ALTERNATIVE TRANSFORMS OF T90")
print("=" * 78)

d = build_frame(extra_master_cols=[T90], cache="frame_main.pkl")
t = d[T90].astype(float)
OUT = ranking_outcomes(d)
assert len(OUT) == N_RANKED_OUTCOMES

d["t90_rint"] = site_rint(d, T90)
d["t90_raw"] = t
d["t90_log1p"] = np.log1p(t)
d["t90_log1p_sz"] = site_z(d.assign(_l=np.log1p(t)), "_l")

# --- restricted cubic splines.
# patsy's `knots=` argument means INTERIOR knots and adds the two boundaries on top, so a
# 3-knot spline is one interior knot with the outer two given as bounds. Getting that wrong
# silently fits a 5-knot basis.
def rcs(x, name, q=(0.10, 0.50, 0.90)):
    """k-knot restricted cubic spline, knots at the given quantiles, one column dropped."""
    s_ = pd.Series(x); kn = [float(s_.quantile(v)) for v in q]
    if len(set(np.round(kn, 9))) < len(kn):
        # v7: 12.9 percent of native-rate T90 is exactly 0 (the oximeter never dipped below 90), so the low quantiles
        # coincide. The outer knot stays at the minimum and the interior knots are the same quantiles of the values above it.
        pos = s_[s_ > s_.min()]
        kn = [float(s_.min())] + [float(pos.quantile(v)) for v in q[1:-1]] + [float(s_.quantile(q[-1]))]
        print(f"   {name}: knots recomputed on the values above the minimum (mass at {s_.min():g}: {(s_ == s_.min()).mean()*100:.1f} percent of nights)")
    assert len(set(np.round(kn, 9))) == len(kn), f"{name} knots collide: {kn}"
    B = dmatrix("cr(x, knots=k, lower_bound=lo, upper_bound=hi) - 1",
                {"x": np.asarray(x, float), "k": kn[1:-1], "lo": kn[0], "hi": kn[-1]},
                return_type="dataframe")
    assert B.shape[1] == len(kn), (name, B.shape[1], len(kn))
    B = B.iloc[:, 1:]                    # one column dropped, the rank trap
    cols = [f"{name}{i}" for i in range(B.shape[1])]
    print(f"   {name}: {len(kn)} knots at " + ", ".join(f"{v:.4f}" for v in kn) +
          f"  ->  {B.shape[1]} columns after dropping one of {len(kn)}")
    return cols, B


print("spline bases")
RCS, B3 = rcs(t.values, "t90_rcs3_")
for c, v in zip(RCS, B3.T.values):
    d[c] = v
RCS5, B5 = rcs(t.values, "t90_rcs5_", q=(0.02, 0.10, 0.50, 0.90, 0.98))
for c, v in zip(RCS5, B5.T.values):
    d[c] = v
RCSZ, BZ = rcs(d.t90_rint.values, "t90_zrcs3_")
for c, v in zip(RCSZ, BZ.T.values):
    d[c] = v

# --- the paper's four bands
bi = band_index(t.values)
assert not np.isnan(bi).any()
for i in range(1, 4):
    d[f"t90_band{i}"] = (bi == i).astype(float)
BAND = [f"t90_band{i}" for i in range(1, 4)]
print("band sizes: " + ", ".join(
    f"{BANDS[i][2]} n={int((bi == i).sum()):,}" for i in range(4))
    + "   (reference band 0-1%)")

specs = {
    "__base__": [],
    "rint_site": ["t90_rint"],
    "raw_linear": ["t90_raw"],
    "log1p_raw": ["t90_log1p"],
    "log1p_sitez": ["t90_log1p_sz"],
    "rcs3_raw": RCS,
    "rcs5_raw": RCS5,
    "rcs3_on_rint": RCSZ,
    "bands4": BAND,
}
res = heldout(specs, OUT, d, penalizer=PEN, pkl_name="D_frame.pkl")
s = summarise(res)
base = float(res[res.spec == "__base__"].c.mean())
pub = __PUB_DC_V7__
ref = float(s[s.spec == "rint_site"].dC.iloc[0])
assert abs(ref - pub) < 1e-9, (ref, pub)

print(f"\nheld-out baseline {base:.6f}   (ranking_v3.csv: 0.647540)")
print(f"published transform reproduces {ref:.6f} exactly\n")
print(f"{'transform':<16}{'cols':>5}{'dC':>12}{'vs published':>14}{'%':>9}"
      f"{'n>0':>5}{'worst':>12}")
print("-" * 76)
for _i, r in s.iterrows():
    n = len(specs[r.spec])
    print(f"{r.spec:<16}{n:>5}{r.dC:>+12.6f}{r.dC - ref:>+14.6f}"
          f"{100 * (r.dC / ref - 1):>+9.1f}{int(r.npos):>5}{r.worst:>+12.6f}")

best = s.iloc[0]
print(f"\nbest transform: {best.spec} at {best.dC:+.6f}, "
      f"{best.dC - ref:+.6f} against the published {ref:+.6f} "
      f"({100 * (best.dC / ref - 1):+.1f}%)")
beat = s[(s.dC > ref) & (s.spec != "rint_site")]
if len(beat):
    print("transforms that beat the published one: " +
          ", ".join(f"{r.spec} {r.dC - ref:+.6f}" for _i, r in beat.iterrows()))
else:
    print("no transform beats the published one.")

nf = json.load(open(f"{paths.NUMBERS_DIR}/ranking_v3_noisefloor.json"))
print(f"\nfor scale, the noise floor: one meaningless column earns "
      f"{nf['dC_per_noise_variable_mean']:+.6f} on this design, and the spread of a single "
      f"outcome's gain has sd {nf['single_outcome_gain_sd']:.4f}")

s.insert(0, "analysis", "D")
s.to_csv(os.path.join(WORK, "D_summary.csv"), index=False)
res.to_csv(os.path.join(WORK, "D_percondition.csv"), index=False)
json.dump({"baseline": base, "published_rint_dC": ref,
           "knot_quantiles_rcs3": [0.10, 0.50, 0.90],
           "knot_quantiles_rcs5": [0.02, 0.10, 0.50, 0.90, 0.98],
           "band_edges": [[b[0], b[1], b[2]] for b in BANDS],
           "dC": {r.spec: float(r.dC) for _i, r in s.iterrows()},
           "ncols": {k: len(v) for k, v in specs.items()},
           "best": str(best.spec), "best_dC": float(best.dC),
           "any_beats_published": bool(len(beat) > 0)},
          open(os.path.join(WORK, "D_result.json"), "w"), indent=2)
print("\nwritten -> _work/D_summary.csv, _work/D_percondition.csv, _work/D_result.json")
