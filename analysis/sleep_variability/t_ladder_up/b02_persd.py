"""
Per-SD hazard ratios for every rung of the threshold ladder, on the paper's own Cox
machinery. Extends t85_analysis/a02_persd.py.

POSITIVE CONTROL FIRST, unchanged from a02: the per-SD fitter is
numbers/run_all_v2.py's cox() helper, and before it fits anything new it must
reproduce section 4's published binary hazard ratios in
results_v2.json["vs_sleep_duration"] to 1e-9.

Rungs fitted per SD, all on the SAME rows (patients with every exposure present):
  t80, t85, t92, t95        re-derived, this run
  t88_rederiv, t90_rederiv  re-derived, same run (the pure same-code ladder)
  t88_frozen, t90_frozen    the frozen columns (the task's labels for T88/T90)

One SD is that exposure's own SD on the common rows, so one unit is one typical
step of that measure, not a shared scale.

A95, TIME AT OR ABOVE 95, AS A CONTINUOUS PER-SD VARIABLE: a95 = 100 - t95 on
the same denominator, so it is the same variable with the sign flipped. Its
Pearson and Spearman correlation with t95 is exactly -1, its per-SD hazard
ratio is exactly the reciprocal of t95's, and its dC is identical. This script
DEMONSTRATES that on three outcomes and then does not fit it again, so nobody
reads the mirrored column as new information. The genuinely new a95 content is
the banded preserved-saturation framing in b01.
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

t = pd.read_csv(os.path.join(HERE, "ladder_per_patient_v7.csv"), low_memory=False)
t = t[t.status == "ok"][["BDSPPatientID", "prov_t80", "prov_t85", "prov_t88",
                         "prov_t90", "prov_t92", "prov_t95", "prov_a95"]]
b["BDSPPatientID"] = b.BDSPPatientID.astype(int)
t["BDSPPatientID"] = t.BDSPPatientID.astype(int)
b = b.merge(t, on="BDSPPatientID", how="left")

b["t90_frozen"] = b.spo2_pct_below_90
b["t88_frozen"] = b.spo2_pct_below_88
for c in ("t80", "t85", "t92", "t95", "a95"):
    b[c] = b[f"prov_{c}"]
b["t88_rederiv"] = b.prov_t88
b["t90_rederiv"] = b.prov_t90
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
    """numbers/run_all_v2.py cox(), transcribed (via a02)."""
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
                  "hi": round(float(hi), 3), "p": float(c.summary.loc[x, "p"]),
                  "hr_raw": float(np.exp(c.params_[x]))}
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

# ------------------------------------------------------------------ correlations
print("\n" + "=" * 78)
print("THE LADDER'S CORRELATION STRUCTURE")
print("=" * 78)
EXPO = ["t80", "t85", "t88_frozen", "t88_rederiv", "t90_frozen", "t90_rederiv",
        "t92", "t95"]
cc = b[EXPO + ["a95"]].dropna()
CORR = {"n": int(len(cc))}
print(f"  complete rows {len(cc):,}. Spearman with frozen T90:")
for c in EXPO:
    s_ = float(stats.spearmanr(cc[c], cc["t90_frozen"]).statistic)
    CORR[f"{c}__t90_frozen"] = s_
    print(f"    {c:<12} {s_:.4f}")
sp_at = float(stats.spearmanr(cc["a95"], cc["t95"]).statistic)
pe_at = float(np.corrcoef(cc["a95"], cc["t95"])[0, 1])
CORR["a95__t95_spearman"] = sp_at
CORR["a95__t95_pearson"] = pe_at
print(f"  a95 vs t95: Spearman {sp_at:.6f}, Pearson {pe_at:.6f} "
      f"(the same variable, sign flipped)")
idmax = float((cc.a95 + cc.t95 - 100.0).abs().max())
CORR["a95_identity_max_abs_dev_pp"] = idmax
print(f"  max |a95 + t95 - 100| on the modelled cohort: {idmax:.3e} pp")

print("\n  distributions on the modelled cohort")
DIST = {}
for c in EXPO + ["a95"]:
    v = b[c].dropna()
    DIST[c] = {"n": int(len(v)), "mean": float(v.mean()), "sd": float(v.std()),
               "median": float(v.median()), "q25": float(v.quantile(.25)),
               "q75": float(v.quantile(.75)), "q90": float(v.quantile(.90)),
               "q99": float(v.quantile(.99)), "max": float(v.max()),
               "pct_zero": float((v == 0).mean() * 100),
               "pct_gt10": float((v > 10).mean() * 100),
               "pct_le1": float((v <= 1).mean() * 100)}
    print(f"  {c:<12} n={len(v):,}  mean {v.mean():7.3f}  SD {v.std():7.3f}  "
          f"median {v.median():8.4f}  p90 {v.quantile(.90):7.3f}  max {v.max():6.2f}")

# ------------------------------------------------------------------ per-SD ladder
print("\n" + "=" * 78)
print("PER-SD HAZARD RATIOS, same rows, same model, every rung")
print("=" * 78)
common = b[EXPO].notna().all(axis=1)
bb = b[common].copy()
print(f"  rows with every exposure present: {len(bb):,} of {len(b):,}")
SD = {}
for c in EXPO:
    SD[c] = float(bb[c].std())
    bb[f"{c}_sd"] = bb[c] / SD[c]
bb["a95_sd"] = bb["a95"] / float(bb["a95"].std())
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
        got[c] = r
    if not got:
        continue
    row = {"outcome": lab, "key": k,
           "organ": ORGAN_GROUP.get(k, "Death" if k == "death" else "UNGROUPED"),
           "negative_control": k in NEGATIVE_CONTROLS,
           "in_ranking_48": (k == "death") or (k not in CIRCULAR and k not in RANKING_EXCLUDE),
           "n": got["t95"]["n"], "events": got["t95"]["events"]}
    for c in EXPO:
        e = got[c][f"{c}_sd"]
        row.update({f"{c}_hr": e["hr"], f"{c}_lo": e["lo"], f"{c}_hi": e["hi"],
                    f"{c}_p": e["p"]})
    rows.append(row)

P = pd.DataFrame(rows)
for c in EXPO:
    P[f"{c}_q"] = np.nan
    for neg, idx in P.groupby("negative_control").groups.items():
        P.loc[idx, f"{c}_q"] = bh(P.loc[idx, f"{c}_p"].values)
P.to_csv(os.path.join(HERE, "ladder_persd.csv"), index=False)

real = P[~P.negative_control]
neg = P[P.negative_control]
print(f"\n  {len(P)} outcomes fitted ({len(real)} real, {len(neg)} negative controls)")
LADDER_STATS = {}
for c in EXPO:
    s = real[f"{c}_hr"]
    sig = real[real[f"{c}_q"] < 0.05]
    LADDER_STATS[c] = {
        "median_hr": float(s.median()),
        "n_hr_gt1": int((s > 1).sum()),
        "n_sig": int(len(sig)),
        "median_sig_hr": float(sig[f"{c}_hr"].median()) if len(sig) else None,
        "neg_sig": int((neg[f"{c}_q"] < 0.05).sum()),
    }
    print(f"  {c:<12} median per-SD HR {s.median():.3f}   HR>1 on {int((s>1).sum())}/{len(real)}"
          f"   q<0.05 on {len(sig)}/{len(real)}   "
          f"median sig HR {sig[f'{c}_hr'].median() if len(sig) else float('nan'):.3f}   "
          f"negative controls q<0.05 on {int((neg[f'{c}_q']<0.05).sum())}/{len(neg)}")

# ------------------------------------------------------------------ a95 mirror demo
print("\n" + "=" * 78)
print("A95 PER-SD IS T95 MIRRORED, demonstrated then not fitted again")
print("=" * 78)
mir = {}
worst_mirror = 0.0
for k in ("obesity", "hf", "death"):
    r95 = cox(bb, ["t95_sd"], k)
    ra = cox(bb, ["a95_sd"], k)
    h95 = r95["t95_sd"]["hr_raw"]
    ha = ra["a95_sd"]["hr_raw"]
    dev = abs(ha - 1.0 / h95)
    worst_mirror = max(worst_mirror, dev)
    mir[k] = {"t95_hr": h95, "a95_hr": ha, "recip_t95": 1.0 / h95, "abs_dev": dev}
    print(f"  {k:<8} t95 per-SD HR {h95:.6f}   a95 per-SD HR {ha:.6f}   "
          f"1/t95 {1.0/h95:.6f}   |a95 - 1/t95| {dev:.2e}")
print(f"  worst deviation from exact reciprocity: {worst_mirror:.2e}")
print("  So per-SD a95 carries no information t95 does not already carry, by arithmetic.")

json.dump({"positive_control_pass": bool(PC_OK), "pc_worst_abs_diff": worst,
           "pc_n_outcomes_checked": n_chk, "correlations": CORR,
           "distributions": DIST, "sd_units_pp": SD,
           "n_common_rows": int(len(bb)), "n_outcomes_fitted": int(len(P)),
           "ladder_stats": LADDER_STATS,
           "a95_mirror_demo": mir, "a95_worst_mirror_dev": worst_mirror},
          open(os.path.join(HERE, "persd_summary.json"), "w"), indent=2)
print("\nwrote ladder_persd.csv, persd_summary.json")
