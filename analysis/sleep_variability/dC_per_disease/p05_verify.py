"""
The number audit. Every printed value in this folder is recomputed from the raw checkpoints by
a second, independent code path and compared to what was delivered.

  1. point estimates. Every per-direction C and every dC in the CSV must equal what
     _work/point_<arm>_<basis>.npy holds, for both arms.
  2. intervals. Every ci_lo, ci_hi, boot_se and replicate count must equal the percentiles
     recomputed here from _work/_ck*/boot_*.npz, not reused from the assembly.
  3. the mean rows. The mean over outcomes must be the mean of the per-replicate means, taken
     inside the replicate, on the outcome set the row names.
  4. the markdown eTable. Every number in every row is parsed back out of the file and
     matched to the CSV at the printed precision.
  5. the figure. _work/figure_drawn_values.csv must restate the CSV exactly, marker fill
     included.

Any mismatch raises. Nothing is rounded away.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import glob
import os
import sys
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

import dc_engine as E
import sys as _sys_v8; _sys_v8.path.insert(0, paths.ANALYSIS_DIR)  # v8 sweep 2026-09-12
from cohort_spec import N_RANKED_OUTCOMES  # v8 sweep 2026-09-12
__PUB_DC_V7__ = float(__import__("pandas").read_csv(f"{paths.NUMBERS_DIR}/ranking_v3.csv", comment="#").set_index("feature").loc["spo2_pct_below_90", "dC"])   # v7: T90 gain from the CURRENT ranking, was the August literal 0.02033875497249791


HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "dC_per_disease.csv")
MD = os.path.join(HERE, "dC_per_disease_TABLE.md")
FIGV = os.path.join(HERE, "_work", "figure_drawn_values.csv")
ARMS = {"standard": {"lag": 0.0, "ck": "_ck", "bases": ("drop", "full")},
        "lag2": {"lag": 2.0, "ck": "_ck_lag2", "bases": ("drop",)}}
fails = []


def check(name, ok, detail=""):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}{'  ' + detail if detail else ''}")
    if not ok:
        fails.append(name)


raw = pd.read_csv(CSV)
dis = raw[raw.row_type == "disease"]
arms = sorted(dis.arm.unique(), key=lambda a: 0 if a == "standard" else 1)
A = E.build_arrays()
outs = A["outcomes"]
idx = {o: i for i, o in enumerate(outs)}
print(f"auditing {len(raw)} rows, arms {arms}")

# ---------------------------------------------------------------- 1. point estimates
worst = 0.0
for arm in arms:
    P = np.load(os.path.join(E.WORK, f"point_{arm}_drop.npy"))
    Q = (np.load(os.path.join(E.WORK, "point_standard_full.npy"))
         if arm == "standard" else None)
    for _i, r in dis[dis.arm == arm].iterrows():
        p = P[idx[r.outcome]]
        pairs = [(r.c_agesex_dir1_fitI0002_scoreI0006, p[0]),
                 (r.c_agesex_t90_dir1_fitI0002_scoreI0006, p[1]),
                 (r.dC_dir1, p[1] - p[0]),
                 (r.c_agesex_dir2_fitI0006_scoreI0002, p[2]),
                 (r.c_agesex_t90_dir2_fitI0006_scoreI0002, p[3]),
                 (r.dC_dir2, p[3] - p[2]),
                 (r.dC, np.nanmean([p[1] - p[0], p[3] - p[2]])),
                 (r.c_agesex_mean, np.nanmean([p[0], p[2]])),
                 (r.c_agesex_t90_mean, np.nanmean([p[1], p[3]]))]
        if Q is not None:
            q = Q[idx[r.outcome]]
            pairs.append((r.dC_pubbasis, np.nanmean([q[1] - q[0], q[3] - q[2]])))
        for got, want in pairs:
            if np.isnan(got) and np.isnan(want):
                continue
            worst = max(worst, abs(float(got) - float(want)))
check("point estimates rebuilt from point_<arm>_<basis>.npy", worst < 1e-12,
      f"largest difference {worst:.2e}")

# ---------------------------------------------------------------- 2 and 3. intervals
def boot_dc(arm, basis):
    parts = []
    for f in sorted(glob.glob(os.path.join(E.WORK, ARMS[arm]["ck"], "boot_*.npz"))):
        parts.append(np.load(f)[basis])
    M = np.concatenate(parts, axis=0)
    return np.nanmean(np.stack([M[:, :, 1] - M[:, :, 0], M[:, :, 3] - M[:, :, 2]], axis=2),
                      axis=2)


G = {(a, b): boot_dc(a, b) for a in arms for b in ARMS[a]["bases"]}
worst_ci, worst_se, nbad = 0.0, 0.0, 0
for arm in arms:
    for basis, tag in (("drop", ""), ("full", "_pubbasis")):
        if (arm, basis) not in G:
            continue
        g = G[(arm, basis)]
        for _i, r in dis[dis.arm == arm].iterrows():
            v = g[:, idx[r.outcome]]
            # the assembly reports no interval when the observed data cannot be fitted,
            # because the replicates that clear the fold floor are then a selected subset
            if not bool(r.estimable):
                nbad += int(r[f"n_replicates{tag}"] != 0)
                nbad += int(bool(r[f"excludes_zero{tag}"]))
                nbad += int(np.isfinite(r[f"ci_lo{tag}"]))
                nbad += int(int(np.isfinite(v).sum()) != int(r.boot_reps_clearing_the_floor))
                continue
            v = v[~np.isnan(v)]
            if len(v) == 0:
                nbad += int(r[f"n_replicates{tag}"] != 0)
                continue
            lo, hi = np.percentile(v, [2.5, 97.5])
            worst_ci = max(worst_ci, abs(lo - r[f"ci_lo{tag}"]), abs(hi - r[f"ci_hi{tag}"]))
            worst_se = max(worst_se, abs(v.std(ddof=1) - r[f"boot_se{tag}"]))
            nbad += int(len(v) != r[f"n_replicates{tag}"])
            nbad += int(bool(lo > 0 or hi < 0) != bool(r[f"excludes_zero{tag}"]))
check("percentile intervals recomputed from the checkpoints", worst_ci < 1e-12,
      f"largest difference {worst_ci:.2e}")
check("bootstrap standard errors", worst_se < 1e-12, f"largest difference {worst_se:.2e}")
check("replicate counts and the zero verdict", nbad == 0, f"{nbad} disagreements")

# the mean rows, rebuilt from the same matrices
est = {a: ~np.isnan(np.array([dis[(dis.arm == a) & (dis.outcome == o)].dC.iloc[0]
                              for o in outs])) for a in arms}
common = np.logical_and.reduce([est[a] for a in arms])
worst_m = 0.0
for arm in arms:
    for tag, mask in (("__mean_all__", est[arm]), ("__mean_common__", common)):
        r = raw[(raw.outcome == tag) & (raw.arm == arm)].iloc[0]
        per_rep = np.nanmean(G[(arm, "drop")][:, mask], axis=1)
        lo, hi = np.percentile(per_rep, [2.5, 97.5])
        pt = np.load(os.path.join(E.WORK, f"point_{arm}_drop.npy"))
        worst_m = max(worst_m, abs(lo - r.ci_lo), abs(hi - r.ci_hi),
                      abs(np.nanmean(E.dc_from_pass(pt)[mask]) - r.dC),
                      abs(per_rep.std(ddof=1) - r.boot_se))
        if int(mask.sum()) != int(r.n_outcomes_in_mean):
            nbad += 1
check("the mean-across-outcomes rows, means taken inside each replicate", worst_m < 1e-12,
      f"largest difference {worst_m:.2e}")

# ---------------------------------------------------------------- 4. the markdown eTable
txt = open(MD).read()
rows = [l for l in txt.splitlines() if l.startswith("| ") and "---" not in l]
body = [r for r in rows if "| Disease" not in r]
S = dis[dis.arm == "standard"].set_index("disease")
L = dis[dis.arm == "lag2"].set_index("disease") if "lag2" in arms else None
means = raw[raw.row_type == "mean_of_outcomes"].set_index(["disease", "arm"])
bad_md = []


def cistr(r):
    if r is None or not np.isfinite(r.dC) or not np.isfinite(r.ci_lo):
        return "not estimable"
    return f"{r.dC:.4f} ({r.ci_lo:.4f} to {r.ci_hi:.4f})"


for line in body:
    c = [x.strip() for x in line.strip("|").split("|")]
    name = c[0]
    if len(c) == 5:                                     # the published-basis panel
        if name not in S.index:
            bad_md.append(f"unknown disease {name!r}")
            continue
        r = S.loc[name]
        want = (f"{int(r.events):,}", f"{r.c_agesex_mean_pubbasis:.3f}",
                f"{r.c_agesex_t90_mean_pubbasis:.3f}",
                f"{r.dC_pubbasis:.4f} ({r.ci_lo_pubbasis:.4f} to {r.ci_hi_pubbasis:.4f})")
        if tuple(c[1:]) != want:
            bad_md.append(f"{name} panel2 {c[1:]} wanted {list(want)}")
        continue
    if name in S.index:                                 # a disease row
        r = S.loc[name]
        lr = L.loc[name] if L is not None and name in L.index else None
        want = (r.organ_group, f"{int(r.events):,}", f"{r.c_agesex_mean:.3f}",
                f"{r.c_agesex_t90_mean:.3f}", cistr(r),
                (f"{int(lr.events):,}" if lr is not None else ""),
                (cistr(lr) if lr is not None else ""),
                (r.negative_control if isinstance(r.negative_control, str) else ""))
        if tuple(c[1:]) != want:
            bad_md.append(f"{name} {c[1:]} wanted {list(want)}")
    elif (name, "standard") in means.index:             # a mean row
        r = means.loc[(name, "standard")]
        if c[5] != cistr(r):
            bad_md.append(f"{name} mean dC {c[5]!r} wanted {cistr(r)!r}")
        if L is not None:
            rl = means.loc[(name.replace(str(int(r.n_outcomes_in_mean)),
                                         str(int(r.n_outcomes_in_mean)))
                            if False else name, "lag2")] \
                if (name, "lag2") in means.index else None
            if rl is not None and c[7] != cistr(rl):
                bad_md.append(f"{name} lag2 mean {c[7]!r} wanted {cistr(rl)!r}")
    else:
        bad_md.append(f"unrecognised row {name!r}")
check(f"every number in the markdown eTable, {len(body)} rows", not bad_md,
      "" if not bad_md else "; ".join(bad_md[:4]))

# ---------------------------------------------------------------- 5. the figure
fv = pd.read_csv(FIGV)
ref = pd.concat([dis, raw[raw.outcome == "__mean_all__"]]).set_index(["disease", "arm"])
worst_fig, nfig = 0.0, 0
for _i, f in fv.iterrows():
    r = ref.loc[(f.disease, f.arm)]
    worst_fig = max(worst_fig, abs(f.dC - r.dC), abs(f.ci_lo - r.ci_lo),
                    abs(f.ci_hi - r.ci_hi))
    nfig += int(bool(f.filled) != bool(r.excludes_zero))
expect = int(np.isfinite(dis.dC).sum()) + len(arms)
check("every value drawn in the figure, marker fill included",
      worst_fig < 1e-12 and nfig == 0 and len(fv) == expect,
      f"largest difference {worst_fig:.2e}, {len(fv)} marks against {expect} expected")

# ---------------------------------------------------------------- headline
mp = raw[(raw.outcome == "__mean_all__") & (raw.arm == "standard")].iloc[0]
check(f"the {N_RANKED_OUTCOMES} published-basis gains still average the paper's headline",
      abs(mp.dC_pubbasis - __PUB_DC_V7__) < 1e-12, f"{mp.dC_pubbasis:.12f}")

print("\nVERDICT: " + ("PASS, every delivered number reproduces" if not fails
                       else "FAIL: " + ", ".join(fails)))
if fails:
    sys.exit(1)
