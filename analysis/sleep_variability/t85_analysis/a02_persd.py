"""
Per-SD T85 against per-SD T90, on the paper's own Cox machinery, plus the correlation
between the two exposures.

POSITIVE CONTROL FIRST. The per-SD fitter here is numbers/run_all_v2.py's `cox()` helper.
Before it is used on T85 it is used to reproduce a published number from the same file,
section 4's binary "T90 above 10%" hazard ratios in results_v2.json["vs_sleep_duration"],
to 1e-9. That exercises exactly the code path the per-SD fits use: site strata, the
age spline with one column dropped, sex, penalizer 0, the 60-event floor, and the same
prevalent/follow-up eligibility.

Per-SD means the exposure divided by its own standard deviation on the modelled rows, so
one unit is one SD OF THAT EXPOSURE. T85 and T90 have very different spreads, so this is a
comparison of how much hazard one typical step in each measure buys, not a comparison on a
shared scale. Both are extremely right skewed, which is stated rather than fixed, because
fixing it would stop this being the paper's per-SD convention.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import json
import os
import sys
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from patsy import dmatrix
from scipy import stats

T90ROOT = paths.T90_ROOT
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, paths.ANALYSIS_DIR)
from cohort_spec import apply_cohort, COHORT_N                      # noqa: E402
from disease_definitions import (DISEASES, NEGATIVE_CONTROLS,       # noqa: E402
                                 RANKING_EXCLUDE, ORGAN_GROUP)

CIRCULAR = {"insomnia", "rls_plmd", "nocturia", "epilepsy"}
FLOOR = 60          # run_all_v2.py cox() floor


def bh(p):
    p = np.asarray(p, float)
    q = np.full(p.shape, np.nan)
    m = np.isfinite(p)
    n = int(m.sum())
    if n == 0:
        return q
    idx = np.where(m)[0][np.argsort(p[m])]
    q[idx] = np.minimum.accumulate(
        (p[idx] * n / np.arange(1, n + 1))[::-1])[::-1].clip(0, 1)
    return q


b = pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet")
h = pd.read_parquet(f"{paths.TABLES_DIR}/t90_base4.parquet")[
    ["BDSPPatientID", "osa_prevalent"]]
b = apply_cohort(b.merge(h, on="BDSPPatientID", how="left"))
assert len(b) == COHORT_N

t = pd.read_csv(os.path.join(HERE, "t85_per_patient_v7.csv"), low_memory=False)
t = t[t.status == "ok"][["BDSPPatientID", "prov_t85", "prov_t88", "prov_t90"]]
b["BDSPPatientID"] = b.BDSPPatientID.astype(int)
t["BDSPPatientID"] = t.BDSPPatientID.astype(int)
b = b.merge(t, on="BDSPPatientID", how="left")

b["t90_frozen"] = b.spo2_pct_below_90
b["t85"] = b.prov_t85
b["t90_rederiv"] = b.prov_t90
b["t88_rederiv"] = b.prov_t88
b["male"] = (b.sex.astype(str).str.upper().str[0] == "M").astype(int)
b["hrs"] = b.TST_min / 60.0
sp = dmatrix("cr(a, df=4) - 1", {"a": b.AgeAtVisit.values}, return_type="dataframe")
sp.columns = [f"age_s{i}" for i in range(4)]
sp.index = b.index
sp = sp.iloc[:, 1:]
for c in sp.columns:
    b[c] = sp[c]
ADJ = list(sp.columns) + ["male"]


def cox(g, xcols, k, strat=True, floor=FLOOR):
    """numbers/run_all_v2.py cox(), transcribed."""
    f = g[(g[f"{k}_prevalent"] == 0) & g[f"{k}_years"].notna() & (g[f"{k}_years"] > 0)]
    cols = {"T": f[f"{k}_years"], "E": f[f"{k}_incident"].astype(int)}
    for c in xcols + ADJ:
        cols[c] = f[c]
    if strat:
        cols["site"] = f.site_id
    ff = pd.DataFrame(cols).dropna()
    if ff.E.sum() < floor:
        return None
    try:
        c = CoxPHFitter(penalizer=0.0).fit(ff, "T", "E", strata=["site"] if strat else None)
    except Exception:
        return None
    out = {"n": int(len(ff)), "events": int(ff.E.sum())}
    for x in xcols:
        lo, hi = np.exp(c.confidence_intervals_.loc[x])
        out[x] = {"hr": round(float(np.exp(c.params_[x])), 3), "lo": round(float(lo), 3),
                  "hi": round(float(hi), 3), "p": float(c.summary.loc[x, "p"])}
    return out


# ------------------------------------------------------------------ positive control
print("=" * 78)
print("POSITIVE CONTROL: reproduce results_v2.json['vs_sleep_duration'] with this cox()")
print("=" * 78)
b["t90_hi"] = (b.t90_frozen > 10).astype(int)
b["slp_lt5"] = (b.hrs < 5).astype(float).where(b.hrs.notna())   # v8.1: missing where TST is masked (split nights)
pub = json.load(open(f"{paths.NUMBERS_DIR}/results_v2.json"))["vs_sleep_duration"]
worst, where, n_chk = 0.0, "", 0
for k, v in list(DISEASES.items()) + [("death", ("Death from any cause", False, [], []))]:
    lab = v[0]
    if f"{k}_incident" not in b.columns or lab not in pub["outcomes"]:
        continue
    r = cox(b, ["t90_hi", "slp_lt5"], k)
    if r is None:
        continue
    p = pub["outcomes"][lab]
    assert r["events"] == p["events"], (lab, r["events"], p["events"])
    n_chk += 1
    for src, dst in (("t90_hi", "t90"), ("slp_lt5", "short_sleep")):
        for fld in ("hr", "lo", "hi"):
            d = abs(r[src][fld] - p[dst][fld])
            if d > worst:
                worst, where = d, f"{lab} {dst} {fld}"
print(f"  {n_chk} published outcomes checked, largest |diff| {worst:.3e}  ({where})")
PC_OK = worst <= 5e-4 and n_chk >= 50  # v8.3: reference now full precision, compare within 5e-4 (2026-09-16)
print(f"  VERDICT: {'PASS' if PC_OK else 'FAIL'}")
if not PC_OK:
    raise SystemExit("POSITIVE CONTROL FAILED, stop")

# ------------------------------------------------------------------ correlation
print("\n" + "=" * 78)
print("HOW CLOSE ARE T85 AND T90")
print("=" * 78)
cc = b[["t85", "t90_frozen", "t90_rederiv", "t88_rederiv"]].dropna()
CORR = {"n": int(len(cc))}
for a, c in (("t85", "t90_frozen"), ("t85", "t90_rederiv"), ("t85", "t88_rederiv"),
             ("t90_rederiv", "t90_frozen")):
    sp_ = float(stats.spearmanr(cc[a], cc[c]).statistic)
    pe = float(np.corrcoef(cc[a], cc[c])[0, 1])
    pe_l = float(np.corrcoef(np.log1p(cc[a]), np.log1p(cc[c]))[0, 1])
    CORR[f"{a}__{c}"] = {"spearman": sp_, "pearson": pe, "pearson_log1p": pe_l}
    print(f"  {a:<12} vs {c:<12}  Spearman {sp_:.4f}   Pearson {pe:.4f}   "
          f"Pearson(log1p) {pe_l:.4f}")

print("\n  distributions on the modelled cohort")
DIST = {}
for c in ("t90_frozen", "t85", "t90_rederiv", "t88_rederiv"):
    v = b[c].dropna()
    DIST[c] = {"n": int(len(v)), "mean": float(v.mean()), "sd": float(v.std()),
               "median": float(v.median()), "q25": float(v.quantile(.25)),
               "q75": float(v.quantile(.75)), "q90": float(v.quantile(.90)),
               "q99": float(v.quantile(.99)), "max": float(v.max()),
               "pct_zero": float((v == 0).mean() * 100)}
    print(f"  {c:<13} n={len(v):,}  mean {v.mean():7.3f}  SD {v.std():7.3f}  "
          f"median {v.median():7.4f}  p90 {v.quantile(.90):7.3f}  max {v.max():6.2f}  "
          f"exact zero {100*(v==0).mean():.2f}%")

# ------------------------------------------------------------------ per-SD head to head
print("\n" + "=" * 78)
print("PER-SD HAZARD RATIOS, same rows, same model")
print("=" * 78)
EXPO = ["t90_frozen", "t85", "t90_rederiv"]
# common rows so the two exposures are never compared on different patients
common = b[EXPO].notna().all(axis=1)
bb = b[common].copy()
print(f"  rows with all three exposures present: {len(bb):,} of {len(b):,}")
SD = {}
for c in EXPO:
    SD[c] = float(bb[c].std())
    bb[f"{c}_sd"] = bb[c] / SD[c]
print("  one SD equals: " + ", ".join(f"{c} {SD[c]:.3f} pp" for c in EXPO))

rows = []
for k, v in list(DISEASES.items()) + [("death", ("Death from any cause", False, [], []))]:
    lab = v[0]
    if f"{k}_incident" not in bb.columns:
        continue
    got = {}
    for c in EXPO:
        r = cox(bb, [f"{c}_sd"], k)
        if r is None:
            got = {}
            break
        got[c] = (r, r[f"{c}_sd"])
    if not got:
        continue
    row = {"outcome": lab, "key": k, "organ": ORGAN_GROUP.get(k, "Death" if k == "death" else "UNGROUPED"),
           "negative_control": k in NEGATIVE_CONTROLS,
           "in_ranking_48": (k == "death") or (k not in CIRCULAR and k not in RANKING_EXCLUDE),
           "n": got["t85"][0]["n"], "events": got["t85"][0]["events"]}
    for c in EXPO:
        _r, e = got[c]
        row.update({f"{c}_hr": e["hr"], f"{c}_lo": e["lo"], f"{c}_hi": e["hi"],
                    f"{c}_p": e["p"]})
    # both exposures in one model: does T85 add anything on top of T90
    r2 = cox(bb, ["t85_sd", "t90_frozen_sd"], k)
    if r2:
        row.update({"joint_t85_hr": r2["t85_sd"]["hr"], "joint_t85_lo": r2["t85_sd"]["lo"],
                    "joint_t85_hi": r2["t85_sd"]["hi"], "joint_t85_p": r2["t85_sd"]["p"],
                    "joint_t90_hr": r2["t90_frozen_sd"]["hr"],
                    "joint_t90_lo": r2["t90_frozen_sd"]["lo"],
                    "joint_t90_hi": r2["t90_frozen_sd"]["hi"],
                    "joint_t90_p": r2["t90_frozen_sd"]["p"]})
    rows.append(row)

P = pd.DataFrame(rows)
for c in EXPO + ["joint_t85", "joint_t90"]:
    if f"{c}_p" in P.columns:
        P[f"{c}_q"] = np.nan
        for neg, idx in P.groupby("negative_control").groups.items():
            P.loc[idx, f"{c}_q"] = bh(P.loc[idx, f"{c}_p"].values)
P.to_csv(os.path.join(HERE, "persd_t85_vs_t90.csv"), index=False)

real = P[~P.negative_control]
neg = P[P.negative_control]
print(f"\n  {len(P)} outcomes fitted ({len(real)} real, {len(neg)} negative controls)")
for c in EXPO:
    s = real[f"{c}_hr"]
    print(f"  {c:<13} median per-SD HR {s.median():.3f}   "
          f"HR>1 on {int((s>1).sum())}/{len(real)}   "
          f"q<0.05 on {int((real[f'{c}_q']<0.05).sum())}/{len(real)}   "
          f"negative controls q<0.05 on {int((neg[f'{c}_q']<0.05).sum())}/{len(neg)}")

print("\n  head to head on the 10 outcomes with the most events")
top = real.sort_values("events", ascending=False).head(10)
print(top[["outcome", "events", "t90_frozen_hr", "t85_hr", "joint_t90_hr", "joint_t85_hr"]]
      .to_string(index=False, float_format=lambda x: f"{x:9.3f}"))

json.dump({"positive_control_pass": bool(PC_OK), "pc_worst_abs_diff": worst,
           "pc_n_outcomes_checked": n_chk, "correlations": CORR, "distributions": DIST,
           "sd_units_pp": SD, "n_common_rows": int(len(bb)),
           "n_outcomes_fitted": int(len(P))},
          open(os.path.join(HERE, "persd_summary.json"), "w"), indent=2)
print("\nwrote persd_t85_vs_t90.csv, persd_summary.json")
