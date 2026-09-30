"""
Recompute the per-window negative-control fits for the SETTLED five-control panel, in the
stage-specific 9,392-patient PAP-off cohort.

attack_controls_full.csv on disk was written before the 2026-08-07 decision and holds
osteoarthritis and fracture, both of which are out of the panel, and does not hold contact
dermatitis or hemorrhoids, both of which are in. eFigure 12D's eight per-window confounding
floors are maxima over that file, so they are floors of a superseded panel.

This reproduces ATTACK 4's controls_for(b, STAGES) exactly: same parquet, same rank-inverse-
normal-within-site exposure on the 30-minute gated columns, same 4-column natural cubic age
spline with ONE COLUMN DROPPED, same site strata, same unpenalized Cox. The three controls that
survive in both panels (back pain, cataract, glaucoma) are refitted here too, purely as a check
that this reproduces the frozen file to the last digit.

Writes attack_controls_full_negpanel_v5.csv beside the original. The original is not touched.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../..")))  # repository root, where paths.py lives
import paths
import sys, warnings
import numpy as np, pandas as pd
from lifelines import CoxPHFitter
from patsy import dmatrix
from scipy import stats
warnings.filterwarnings("ignore")
np.random.seed(20260807)

HERE = f"{paths.SV_ROOT}/stage_specific"
TROOT = paths.FIGURE_ROOT
sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
from disease_definitions import NEGATIVE_CONTROLS

STAGES = ["all", "wake", "sleep", "nrem", "n1", "n2", "n3", "rem"]

b = pd.read_parquet(f"{HERE}/analysis.parquet")
assert (b.oximetry_bad == 0).all() and (b.fu_valid == 1).all() and b.papoff_strict.all()
print(f"cohort {len(b):,}")

sp = dmatrix("cr(a, df=4) - 1", {"a": b.AgeAtVisit.values}, return_type="dataframe").iloc[:, 1:]
for i in range(sp.shape[1]):
    b[f"age_s{i}"] = sp.iloc[:, i].values
ADJ = [f"age_s{i}" for i in range(sp.shape[1])] + ["male"]
assert len(ADJ) == 4, ADJ            # 3 spline columns + male: one column DROPPED


def rint(x):
    r = stats.rankdata(x)
    return stats.norm.ppf((r - 0.375) / (len(r) + 0.25))


def rint_site(frame, col):
    out = pd.Series(np.nan, index=frame.index)
    for _, g in frame.groupby("site_id"):
        v = g[col].dropna()
        if len(v) > 10:
            out.loc[v.index] = rint(v.values)
    return out


for s in STAGES:
    b[f"z_{s}"] = rint_site(b, f"t90_{s}_g30")


def frame_for(src, key, zcols):
    f = src[(src[f"{key}_prevalent"] == 0) & src[f"{key}_years"].notna()
            & (src[f"{key}_years"] > 0)]
    d = pd.DataFrame({"T": f[f"{key}_years"].values,
                      "E": f[f"{key}_incident"].astype(int).values,
                      "site": f.site_id.values}, index=f.index)
    for c in zcols + ADJ:
        d[c] = f[c].values
    return d


def fit(d, col):
    """Unpenalized site-stratified Cox. lifelines' default Newton step lands on a wrong
    optimum on some fits, so several step sizes are scanned and the maximum-likelihood
    solution is kept."""
    dd = d[[col] + ADJ + ["T", "E", "site"]].dropna()
    if dd.E.sum() < 20:
        return None
    best = None
    for step in (None, 0.9, 0.75, 0.5, 0.25, 0.1):
        try:
            c = CoxPHFitter(penalizer=0.0)
            kw = {} if step is None else {"step_size": step}
            c.fit(dd, "T", "E", strata=["site"], **kw)
        except Exception:
            continue
        ll = c.log_likelihood_
        if not np.isfinite(ll):
            continue
        if best is None or ll > best[0] + 1e-9:
            best = (ll, c)
    if best is None:
        return None
    r = best[1].summary.loc[col]
    return {"hr": float(r["exp(coef)"]), "coef": float(r["coef"]),
            "se": float(r["se(coef)"]), "p": float(r["p"]), "loglik": best[0]}


_frozen_controls = set(pd.read_csv(f"{HERE}/attack_controls_full.csv").control)
CHECK = [k for k in NEGATIVE_CONTROLS if k in _frozen_controls]   # v8.2: the controls in both the settled panel and the frozen file, data-driven (was a typed list)
assert CHECK, "no control shared between disease_definitions.NEGATIVE_CONTROLS and attack_controls_full.csv"
rows = []
for k in NEGATIVE_CONTROLS:
    d = frame_for(b, k, [f"z_{s}" for s in STAGES])
    rec = {"control": k, "events": int(d.E.sum())}
    line = f"{k:22s}{int(d.E.sum()):>6}"
    for s in STAGES:
        r = fit(d, f"z_{s}")
        if r:
            rec[f"{s}_hr"], rec[f"{s}_coef"] = r["hr"], r["coef"]
            rec[f"{s}_se"], rec[f"{s}_p"] = r["se"], r["p"]
            line += f"{r['hr']:>9.3f}" + ("*" if r["p"] < 0.05 else " ")
        else:
            line += f"{'.':>10}"
    rows.append(rec)
    print(line, flush=True)

out = pd.DataFrame(rows)
old = pd.read_csv(f"{HERE}/attack_controls_full.csv").set_index("control")
print("\nreproduction check against the frozen attack_controls_full.csv")
worst = 0.0
for k in CHECK:
    for s in STAGES:
        a = float(out.loc[out.control == k, f"{s}_hr"].iloc[0])
        c = float(old.loc[k, f"{s}_hr"])
        worst = max(worst, abs(a - c))
        flag = "" if abs(a - c) < 5e-4 else "   <-- DIFFERS"
        print(f"  {k:16s} {s:6s} refit {a:.6f}  frozen {c:.6f}  d={a - c:+.2e}{flag}")
print(f"max abs difference on the three shared controls: {worst:.3e}")

out.to_csv(f"{HERE}/attack_controls_full_negpanel_v5.csv", index=False)
print(f"\nwrote {HERE}/attack_controls_full_negpanel_v5.csv")

print("\nfloor per window = max(worst single control HR, upper 95% of the pooled control)")
for s in STAGES:
    w = 1.0 / out[f"{s}_se"] ** 2
    cf = float((out[f"{s}_coef"] * w).sum() / w.sum())
    hi = float(np.exp(cf + 1.96 * np.sqrt(1.0 / w.sum())))
    worst_hr = float(out[f"{s}_hr"].max())
    drv = str(out.loc[out[f"{s}_hr"].idxmax(), "control"])
    oldw = float(old[f"{s}_hr"].max()); olddrv = str(old[f"{s}_hr"].idxmax())
    ow = 1.0 / old[f"{s}_se"] ** 2
    ocf = float((old[f"{s}_coef"] * ow).sum() / ow.sum())
    ohi = float(np.exp(ocf + 1.96 * np.sqrt(1.0 / ow.sum())))
    print(f"  {s:6s} new {max(worst_hr, hi):.4f} (worst {worst_hr:.4f} {drv}, pooled_hi {hi:.4f})"
          f"   old {max(oldw, ohi):.4f} (worst {oldw:.4f} {olddrv}, pooled_hi {ohi:.4f})")
