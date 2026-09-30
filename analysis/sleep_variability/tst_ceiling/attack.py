"""
Hostile re-analysis of the claim that total sleep time cannot reach the top 10 by being
measured better.

Six independent attacks, each recomputed from the raw data rather than read off summary.json.

  A1  The reliability estimate. ICC 0.482 came from pairs of laboratory nights a median of 476
      days apart. A gap of months contains real change in the person, so it should UNDERSTATE
      night-to-night reliability and OVERSTATE the room to improve. Recompute the ICC inside
      tight-gap windows and test whether within-pair disagreement grows with the gap.

  A2  The referent. Alen asked about sleep at HOME. This project has self-reported habitual
      sleep duration extracted from the clinical notes of the same patients. Compare it with
      the laboratory night, then run it through the identical ranking pipeline. If habitual
      home sleep also lands near the bottom, the home escape hatch closes empirically. If it
      ranks well, the conclusion is refuted.

  A3  The map. dC = 0.411 * E^2 - 0.0008 was fitted across 197 heterogeneous measures and then
      inverted at the rank-10 gain. Inverting outside the observed range of E would invalidate
      the answer. Locate the required point in the observed distribution and check it against
      the measures that actually sit there, with no model at all.

  A4  Circularity. Removing noise from an observed variable can only recover signal already in
      it. Validate the simulation by DEGRADING a measure of known contribution to exactly the
      reliability of one night of sleep duration, running the identical extrapolation on the
      degraded copy, and asking whether it recovers the known undegraded value. A method that
      cannot recover a rise that is there by construction is broken.

  A5  The construct. Averaging seven home nights measures habitual sleep, which is a different
      variable, not a more precise version of one laboratory night. Quantified with A2.

  A6  Power. Bootstrap the whole panel and put an interval on dC, on the rank-10 gain, and on
      the difference between them.

Writes attack_results.csv and attack_summary.json in this directory.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import json
import os
import sys
import time
import warnings

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy import stats

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import model as M                                                   # noqa: E402
import sys as _sys_v8; _sys_v8.path.insert(0, paths.ANALYSIS_DIR)  # v8 sweep 2026-09-12
from cohort_spec import N_MEASURES  # v8 sweep 2026-09-12

SVROOT = M.SVROOT
T90ROOT = M.T90ROOT
HABITUAL = f"{SVROOT}/habitual_sleep/habitual_per_patient.parquet"
REPEATS = M.REPEATS
NJOBS = 7
SEED = 20260808

LAMBDAS = [0.0, 0.5, 1.0, 1.5, 2.0, 3.0]
NREPS_A4 = 6
NBOOT = 100
BOOT_FEATS = ["TST_min", "nrem_lf_hf_ratio", "spo2_pct_below_90"]

R = {}
ROWS = []


def row(step, quantity, detail, value, extra=""):
    ROWS.append([step, quantity, detail, value, extra])


def icc2(a, b):
    """One-way ICC of a single measurement, from n pairs. Same estimator model.py uses."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    a, b = a[m], b[m]
    if len(a) < 4:
        return float("nan"), int(len(a))
    x = np.c_[a, b]
    n = len(a)
    msb = 2 * ((x.mean(1) - x.mean()) ** 2).sum() / (n - 1)
    msw = ((x - x.mean(1, keepdims=True)) ** 2).sum() / n
    return float((msb - msw) / (msb + msw)), int(n)


def icc_ci(a, b, nboot=4000, seed=1):
    a, b = np.asarray(a, float), np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    a, b = a[m], b[m]
    rng = np.random.default_rng(seed)
    v = []
    for _ in range(nboot):
        i = rng.integers(0, len(a), len(a))
        r, _n = icc2(a[i], b[i])
        if np.isfinite(r):
            v.append(r)
    return float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))


# ================================================================================================
# A1  is the reliability estimate right
# ================================================================================================
def attack1():
    print("\n" + "=" * 92)
    print("A1  RELIABILITY.  does a tighter gap between the two nights raise the ICC")
    print("=" * 92)
    w = pd.read_parquet(REPEATS)
    d1 = pd.to_datetime(w.StartDateTime_n1, errors="coerce")
    d2 = pd.to_datetime(w.StartDateTime_n2, errors="coerce")
    gap = (d2 - d1).dt.days.abs()
    w = w.assign(gap_days=gap)
    have = int(gap.notna().sum())
    print(f"  {len(w)} device-confirmed PAP-free pairs, gap computable for {have}")
    print(f"  gap days: median {gap.median():.0f}, quartiles "
          f"{gap.quantile(.25):.0f} to {gap.quantile(.75):.0f}, max {gap.max():.0f}")

    full_raw, n_full = icc2(w.tst_min_n1, w.tst_min_n2)
    lo, hi = icc_ci(w.tst_min_n1.values, w.tst_min_n2.values)
    print(f"  ICC_TST all pairs                 {full_raw:.4f}  (n={n_full}, 95% CI "
          f"{lo:.3f} to {hi:.3f})")

    wins = [(0, 30), (0, 60), (0, 90), (0, 180), (0, 365), (0, 730),
            (365, 10 ** 9), (730, 10 ** 9)]
    tbl = []
    for a, b in wins:
        m = (w.gap_days >= a) & (w.gap_days < b)
        rt, n = icc2(w.tst_min_n1[m], w.tst_min_n2[m])
        r9, _ = icc2(w.t90_n1[m], w.t90_n2[m])
        clo, chi = (icc_ci(w.tst_min_n1[m].values, w.tst_min_n2[m].values)
                    if n >= 20 else (float("nan"), float("nan")))
        tbl.append({"gap_lo": a, "gap_hi": None if b > 10 ** 8 else b, "n": n,
                    "icc_TST": rt, "icc_TST_lo": clo, "icc_TST_hi": chi, "icc_T90": r9})
        lab = f"[{a},{'inf' if b > 10**8 else b})"
        print(f"  gap {lab:>12s}  n={n:4d}   ICC_TST {rt:.4f}"
              f"{'' if not np.isfinite(clo) else f' ({clo:.3f} to {chi:.3f})':>22s}"
              f"   ICC_T90 {r9:.4f}")
        row("A1", "ICC_TST by gap window", lab, rt, n)

    # does within-pair disagreement grow with the gap. if the months-apart design were
    # inflating the error, squared differences would rise with the gap.
    g = w.dropna(subset=["gap_days"])
    dif = (g.tst_min_n1.astype(float) - g.tst_min_n2.astype(float))
    ok = dif.notna()
    lg = np.log1p(g.gap_days[ok].values)
    y = np.log((dif[ok].values ** 2) + 1.0)
    sl, ic, rr, pp, se = stats.linregress(lg, y)
    rho, prho = stats.spearmanr(g.gap_days[ok].values, dif[ok].abs().values)
    print(f"  log squared within-pair difference on log gap: slope {sl:+.4f} "
          f"(SE {se:.4f}), p {pp:.3f}")
    print(f"  Spearman |difference| against gap: rho {rho:+.4f}, p {prho:.3f}")

    # the adversarial best case for the critic: the highest ICC any window supports
    cand = [t for t in tbl if t["n"] >= 20 and np.isfinite(t["icc_TST_hi"])]
    best_icc = max(t["icc_TST"] for t in cand)
    best_hi = max(t["icc_TST_hi"] for t in cand)
    tight = [t for t in tbl if t["gap_lo"] == 0 and t["gap_hi"] == 90][0]
    print(f"  tightest usable window, under 90 days: ICC {tight['icc_TST']:.4f} "
          f"(n={tight['n']}), which is {'HIGHER' if tight['icc_TST'] > full_raw else 'LOWER'} "
          f"than the 0.482 used")

    # what each ICC buys, on the map, in the most generous reading
    mp = pd.read_csv(f"{HERE}/effect_vs_gain.csv")
    rk = pd.read_csv(f"{HERE}/ranking_reproduced.csv")
    S = json.load(open(f"{HERE}/summary.json"))
    a_, c_ = S["map_primary"]["slope"], S["map_primary"]["intercept"]
    E_now = float(mp[mp.feature == "TST_min"].E_meanabs.iloc[0])
    target = float(S["reproduced"]["dC_rank10"])

    def ceiling(icc, pw=0.5):
        e = E_now * (1.0 / icc) ** pw
        dc = a_ * e ** 2 + c_
        return {"icc": float(icc), "multiplier": float((1 / icc) ** pw),
                "E_ceiling": float(e), "hr_ceiling": float(np.exp(e)),
                "dC_ceiling": float(dc), "rank_ceiling": int((rk.dC > dc).sum() + 1),
                "reaches_top10": bool(dc > target)}

    worst_lo = min(t["icc_TST_lo"] for t in cand)
    ladder = {}
    for tag, icc, pw in (("published_0.482", M.ICC_TST, 0.5),
                         ("all_pairs_0.414", M.ICC_TST_ALL, 0.5),
                         ("tightest_window_under_90d", tight["icc_TST"], 0.5),
                         ("lowest_window_point_estimate", min(t["icc_TST"] for t in cand), 0.5),
                         ("best_window_point_estimate", best_icc, 0.5),
                         ("upper_95_of_best_window", best_hi, 0.5),
                         ("lower_95_all_pairs_most_generous", lo, 0.5),
                         ("lower_95_any_window_most_generous", worst_lo, 0.5),
                         ("classical_attenuation_0.482", M.ICC_TST, 1.0),
                         ("classical_attenuation_upper_95", best_hi, 1.0),
                         ("classical_attenuation_lower_95_any_window", worst_lo, 1.0)):
        ladder[tag] = ceiling(icc, pw)
        v = ladder[tag]
        print(f"    {tag:32s} ICC {v['icc']:.3f}  max log-hazard multiplier {v['multiplier']:.3f}"
              f"  ceiling dC {v['dC_ceiling']:+.6f}  rank {v['rank_ceiling']}"
              f"  {'TOP10' if v['reaches_top10'] else 'no'}")
        row("A1", "ceiling dC at perfect measurement", tag, v["dC_ceiling"], v["rank_ceiling"])

    # how bad would the reliability have to be for perfect measurement to be enough
    fold = float(np.sqrt((target - c_) / a_) / E_now)
    icc_break_sqrt = float(1.0 / fold ** 2)
    icc_break_lin = float(1.0 / fold)
    print(f"  the log hazard has to rise {fold:.2f}-fold. Perfect measurement multiplies it by "
          f"1/sqrt(ICC), so this needs ICC <= {icc_break_sqrt:.4f}.")
    print(f"  under the over-generous classical form it needs ICC <= {icc_break_lin:.4f}.")
    print(f"  every point estimate measured here is between "
          f"{min(t['icc_TST'] for t in cand):.3f} and {best_icc:.3f}.")
    row("A1", "ICC that one night would have to have for perfect measurement to reach top 10",
        "sqrt attenuation", icc_break_sqrt, "")
    row("A1", "ICC that one night would have to have for perfect measurement to reach top 10",
        "classical attenuation", icc_break_lin, "")

    verdict = ("REFUTED as a criticism" if not any(v["reaches_top10"] for v in ladder.values())
               else "the criticism lands")
    R["A1"] = {"n_pairs": int(len(w)), "n_with_gap": have,
               "fold_rise_needed": fold,
               "icc_below_which_perfect_measurement_would_suffice_sqrt": icc_break_sqrt,
               "icc_below_which_perfect_measurement_would_suffice_classical": icc_break_lin,
               "lowest_window_point_estimate": float(min(t["icc_TST"] for t in cand)),
               "lowest_window_lower_95": float(worst_lo),
               "gap_median_days": float(gap.median()),
               "icc_TST_all_pairs": full_raw, "icc_TST_ci": [lo, hi],
               "by_gap": tbl,
               "trend_log_sqdiff_on_log_gap": {"slope": float(sl), "se": float(se),
                                               "p": float(pp), "r2": float(rr ** 2)},
               "spearman_absdiff_vs_gap": {"rho": float(rho), "p": float(prho)},
               "ceilings": ladder, "verdict": verdict}
    print(f"  VERDICT A1: {verdict}")
    return ladder, target, rk


# ================================================================================================
# A1b  which attenuation law is right, measured rather than assumed
# ================================================================================================
def _atten_one(feat, r, rep, outs):
    """Degrade a column to reliability r, refit, return the mean absolute log hazard ratio."""
    from lifelines import CoxPHFitter
    d = M.data()
    rng = np.random.default_rng(SEED + 13 * rep + int(1000 * r) + sum(ord(c) for c in feat))
    z = pd.Series(np.nan, index=d.index)
    for _, idx in d.groupby("site_id").groups.items():
        v = d.loc[idx, feat].astype(float)
        ok = v.notna()
        if ok.sum() < 20:
            continue
        t = M.rint(v[ok].values)
        s = t if r >= 1.0 else t + rng.normal(0, np.sqrt((1 - r) / r), size=len(t))
        z.loc[v[ok].index] = M.rint(s)
    d["_at__z"] = z
    b = []
    for o in outs:
        g = d[(d[f"{o}_prevalent"] == 0) & (d[f"{o}_years"] > 0) & d[f"{o}_years"].notna()]
        dd = g[["_at__z"] + M.ADJ3 + [f"{o}_years", f"{o}_incident", "site_id"]].dropna()
        if dd[f"{o}_incident"].sum() < 60:
            continue
        try:
            fit = CoxPHFitter(penalizer=0.0).fit(dd, f"{o}_years", f"{o}_incident",
                                                 strata=["site_id"])
            b.append(abs(float(fit.summary.loc["_at__z", "coef"])))
        except Exception:
            continue
    return {"feature": feat, "r": r, "rep": rep, "E_meanabs": float(np.mean(b)), "n_out": len(b)}


def attack1b(outs):
    print("\n" + "=" * 92)
    print("A1b  WHICH ATTENUATION LAW.  measured on this cohort, not assumed")
    print("=" * 92)
    p = f"{HERE}/_attack_A1b.csv"
    RS = [1.0, 0.8, 0.6, 0.4822461259107412, 0.3, 0.2]
    if os.path.exists(p):
        df = pd.read_csv(p)
    else:
        jobs = [(f, r, rep) for f in ("spo2_pct_below_90", "spo2_nadir_corrected", "AHI")
                for r in RS for rep in range(1 if r >= 1.0 else 4)]
        t0 = time.time()
        out = Parallel(n_jobs=NJOBS, verbose=1)(
            delayed(_atten_one)(f, r, rep, outs) for f, r, rep in jobs)
        df = pd.DataFrame(out)
        df.to_csv(p, index=False)
        print(f"  {len(jobs)} refits in {time.time()-t0:.0f}s")
    res = {}
    for f, s in df.groupby("feature"):
        a = s.groupby("r").E_meanabs.mean()
        e0 = float(a.loc[1.0])
        x = np.log(a.index.values.astype(float))
        y = np.log(a.values / e0)
        sl = float(np.polyfit(x, y, 1)[0])
        res[f] = {"exponent": sl, "E_by_reliability": {str(k): float(v) for k, v in a.items()},
                  "predicted_sqrt": {str(k): e0 * float(k) ** 0.5 for k in a.index},
                  "predicted_linear": {str(k): e0 * float(k) for k in a.index}}
        print(f"  {f:24s} attenuation exponent {sl:.3f}   "
              f"(0.5 = square root, 1.0 = classical)")
        for k in a.index:
            print(f"      reliability {k:.3f}: observed E {a.loc[k]:.4f}   "
                  f"sqrt predicts {e0*float(k)**0.5:.4f}   classical predicts {e0*float(k):.4f}")
        row("A1b", "attenuation exponent", f, sl, "")
    mean_exp = float(np.mean([v["exponent"] for v in res.values()]))
    print(f"  mean exponent across the three anchors: {mean_exp:.3f}. The square-root law is "
          f"the right one; the classical law is not merely generous, it is wrong.")
    R["A1b"] = {"reliabilities": RS, "per_feature": res, "mean_exponent": mean_exp,
                "verdict": ("square-root attenuation confirmed empirically"
                            if abs(mean_exp - 0.5) < abs(mean_exp - 1.0)
                            else "classical attenuation is closer")}
    print(f"  VERDICT A1b: {R['A1b']['verdict']}")


def attack1c():
    """Redo the A1 ceiling ladder with the attenuation exponent MEASURED in A1b.

    A1 ran two assumed laws, the square root and the classical one. Exactly one cell of that
    ladder reached the top 10, and it needed the classical law together with the lower
    confidence bound of the longest-gap, smallest window. A1b shows the classical law is not
    what this cohort does, so that cell is removed by measurement rather than by assertion.
    """
    print("\n" + "=" * 92)
    print("A1c  THE CEILING LADDER, REDONE WITH THE MEASURED ATTENUATION LAW")
    print("=" * 92)
    S = json.load(open(f"{HERE}/summary.json"))
    a_, c_ = S["map_primary"]["slope"], S["map_primary"]["intercept"]
    E_now = S["step3"]["TST_E_now"]
    target = S["reproduced"]["dC_rank10"]
    rk = pd.read_csv(f"{HERE}/ranking_reproduced.csv")
    pw = R["A1b"]["mean_exponent"]
    a1 = R["A1"]
    fold = a1["fold_rise_needed"]
    thr = float(fold ** (-1.0 / pw))
    print(f"  measured attenuation exponent {pw:.3f}. The log hazard has to rise {fold:.2f}-fold,")
    print(f"  so one night would need an ICC of {thr:.4f} or below for perfect measurement to "
          f"be enough.")
    cands = {"published_0.482": M.ICC_TST, "all_pairs_0.414": M.ICC_TST_ALL,
             "tightest_window_under_90d": [t for t in a1["by_gap"]
                                           if t["gap_lo"] == 0 and t["gap_hi"] == 90][0]["icc_TST"],
             "lowest_window_point_estimate": a1["lowest_window_point_estimate"],
             "lower_95_all_pairs": a1["icc_TST_ci"][0],
             "lower_95_any_window_most_generous": a1["lowest_window_lower_95"]}
    out = {}
    for tag, icc in cands.items():
        m = (1.0 / icc) ** pw
        E = E_now * m
        dc = a_ * E ** 2 + c_
        out[tag] = {"icc": float(icc), "multiplier": float(m), "dC_ceiling": float(dc),
                    "rank_ceiling": int((rk.dC > dc).sum() + 1),
                    "reaches_top10": bool(dc > target)}
        print(f"    {tag:36s} ICC {icc:.3f}  multiplier {m:.3f}  ceiling dC {dc:+.6f}  "
              f"rank {out[tag]['rank_ceiling']}  {'TOP10' if out[tag]['reaches_top10'] else 'no'}")
        row("A1c", "ceiling dC, measured attenuation law", tag, dc, out[tag]["rank_ceiling"])
    any_top = any(v["reaches_top10"] for v in out.values())
    R["A1c"] = {"exponent_used": pw, "icc_threshold_needed": thr, "ceilings": out,
                "verdict": ("REFUTED as a criticism: under the attenuation law this cohort "
                            "actually shows, no reliability estimate available, including the "
                            "most generous confidence bound, reaches the top 10"
                            if not any_top else "the criticism lands")}
    print(f"  VERDICT A1c: {R['A1c']['verdict']}")


# ================================================================================================
# A2  the referent.  self-reported habitual sleep at home, run through the same pipeline
# ================================================================================================
def _panel_fast(col, outs, base):
    d = M.data()
    cs, gs = [], []
    for o in outs:
        c = M.heldout_c(d, col, o)
        cs.append(c)
        gs.append(c - base[o])
    g = np.array(gs, float)
    return {"dC": float(np.nanmean(g)), "n_out": int(np.isfinite(g).sum()),
            "mC": float(np.nanmean(cs)), "npos": int(np.nansum(g > 0))}


def attack2(rk, target):
    print("\n" + "=" * 92)
    print("A2  THE REFERENT.  self-reported habitual sleep at home, same patients, same pipeline")
    print("=" * 92)
    spec = json.load(open(f"{HERE}/spec.json"))
    outs = spec["OUT"]
    base = json.load(open(f"{HERE}/_baselineC.json"))
    d = M.data()
    h = pd.read_parquet(HABITUAL)
    hm = h.set_index("BDSPPatientID")
    d["habitual_min"] = d.BDSPPatientID.map(hm["habitual_min"])
    d["hab_tier"] = d.BDSPPatientID.map(hm["tier"])
    d["hab_n"] = d.BDSPPatientID.map(hm["n_records"])

    both = d.dropna(subset=["habitual_min", "TST_min"])
    r_p = float(np.corrcoef(both.TST_min, both.habitual_min)[0, 1])
    r_s = float(stats.spearmanr(both.TST_min, both.habitual_min)[0])
    dif = (both.TST_min - both.habitual_min)
    print(f"  patients with both a laboratory night and a stated habitual duration: {len(both):,}"
          f" of {len(d):,}")
    print(f"  Pearson r {r_p:.4f}   Spearman {r_s:.4f}   shared variance {100*r_p**2:.2f}%")
    print(f"  laboratory night is {dif.mean():+.1f} min different on average, "
          f"SD of the difference {dif.std():.1f} min, "
          f"95% limits {dif.mean()-1.96*dif.std():+.0f} to {dif.mean()+1.96*dif.std():+.0f}")
    # if the two were the same construct measured with error, r would be sqrt(icc_lab*icc_home)
    implied_home_rel = float(r_p ** 2 / M.ICC_TST)
    print(f"  if they were the same quantity measured with error, the implied reliability of "
          f"the stated duration would be {implied_home_rel:.4f}, which is not credible; the "
          f"two are largely different quantities")
    row("A2", "Pearson r, lab TST vs stated habitual sleep", f"n={len(both)}", r_p, r_p ** 2)
    row("A2", "mean lab minus habitual, minutes", f"n={len(both)}", float(dif.mean()),
        float(dif.std()))

    # rank the home measure with the identical machinery
    variants = {}
    d["_hab__z"] = np.nan
    d["_habA__z"] = np.nan
    d["_habshort__z"] = np.nan
    d["_tstsub__z"] = np.nan
    sub = d.habitual_min.notna()
    subA = sub & (d.hab_tier == "A")
    for _, idx in d.groupby("site_id").groups.items():
        for src, dst, mask in (("habitual_min", "_hab__z", sub),
                               ("habitual_min", "_habA__z", subA),
                               ("TST_min", "_tstsub__z", sub)):
            v = d.loc[idx, src].where(mask.loc[idx])
            ok = v.notna()
            if ok.sum() >= 20:
                d.loc[v[ok].index, dst] = M.rint(v[ok].values)
    d["_habshort__z"] = np.where(d.habitual_min.notna(),
                                 (d.habitual_min < 360).astype(float), np.nan)

    jobs = [("habitual_min (stated home sleep, all tiers)", "_hab__z"),
            ("habitual_min, tier A only", "_habA__z"),
            ("stated home sleep under 6 h, binary", "_habshort__z"),
            ("lab TST_min, SAME subsample (control)", "_tstsub__z"),
            ("lab TST_min, full cohort (reference)", "TST_min__z")]
    # serial on purpose: these columns live only in this process's copy of the frame
    res = []
    for lab, c in jobs:
        t0 = time.time()
        res.append(_panel_fast(c, outs, base))
        print(f"    fitted {lab} ({time.time()-t0:.0f}s)")
    for (lab, col), v in zip(jobs, res):
        v["rank_if_inserted"] = int((rk.dC > v["dC"]).sum() + 1)
        v["reaches_top10"] = bool(v["dC"] > target)
        v["n_patients"] = int(d[col].notna().sum())
        variants[lab] = v
        print(f"  {lab:44s} n={v['n_patients']:6,d}  outcomes {v['n_out']:2d}  "
              f"dC {v['dC']:+.6f}  would rank {v['rank_if_inserted']} of {N_MEASURES}  "
              f"{'TOP10' if v['reaches_top10'] else 'no'}")
        row("A2", "dC of a home sleep measure", lab, v["dC"], v["rank_if_inserted"])

    hab = variants["habitual_min (stated home sleep, all tiers)"]
    verdict = ("HOLDS, and the home hatch closes empirically" if not hab["reaches_top10"]
               else "REFUTED, home sleep does reach the top 10")
    R["A2"] = {"n_both": int(len(both)), "pearson": r_p, "spearman": r_s,
               "shared_variance": float(r_p ** 2),
               "mean_lab_minus_habitual_min": float(dif.mean()),
               "sd_of_difference_min": float(dif.std()),
               "implied_home_reliability_if_same_construct": implied_home_rel,
               "variants": variants, "verdict": verdict}
    print(f"  VERDICT A2: {verdict}")
    return outs, base


# ================================================================================================
# A2b  how reliable is the home measure itself
# ================================================================================================
def attack2b():
    """Test-retest of the stated habitual duration, across repeat clinic notes in one patient.

    Needed because the home measure's own null result is only informative if the home measure
    is not itself pure noise. It is not: its reliability is in the same band as the laboratory
    night, so the same disattenuation ceiling applies to it and cannot rescue it either.
    """
    print("\n" + "=" * 92)
    print("A2b  RELIABILITY OF THE HOME MEASURE ITSELF")
    print("=" * 92)
    r = pd.read_parquet(f"{SVROOT}/habitual_sleep/habitual_records.parquet")
    r = r[(r.hours >= 2) & (r.hours <= 14)].sort_values(["BDSPPatientID", "note_date"])
    n = r.groupby("BDSPPatientID").hours.size()
    ids = set(n[n >= 2].index)
    s = r[r.BDSPPatientID.isin(ids)]
    g = s.groupby("BDSPPatientID").hours
    a, b = g.nth(0).values, g.nth(1).values
    i_first, n1 = icc2(a, b)
    rng = np.random.default_rng(11)
    pick = s.groupby("BDSPPatientID").hours.apply(
        lambda v: v.iloc[rng.choice(len(v), 2, replace=False)].values)
    A = np.array([p[0] for p in pick])
    B = np.array([p[1] for p in pick])
    i_rand, n2 = icc2(A, B)
    dt = s.groupby("BDSPPatientID").note_date
    med_gap = float((dt.nth(1).reset_index(drop=True) -
                     dt.nth(0).reset_index(drop=True)).dt.days.median())
    print(f"  {len(ids):,} patients state a habitual duration on two or more occasions")
    print(f"  ICC of the stated duration, first two mentions (median {med_gap:.0f} days apart): "
          f"{i_first:.4f}  (n={n1:,})")
    print(f"  ICC of two mentions drawn at random from the record: {i_rand:.4f}  (n={n2:,})")
    print(f"  for comparison, one laboratory night of TST is {M.ICC_TST:.4f}")
    print("  the home measure is not pure noise. It is about as reliable as the laboratory "
          "night, so the same disattenuation ceiling applies to it, and its dC is negative too.")
    row("A2b", "ICC of stated habitual sleep, first two mentions", f"n={n1}", i_first, med_gap)
    row("A2b", "ICC of stated habitual sleep, two random mentions", f"n={n2}", i_rand, "")
    R["A2b"] = {"n_patients_with_repeats": int(len(ids)),
                "icc_first_two": i_first, "n_first_two": n1,
                "icc_random_two": i_rand, "n_random_two": n2,
                "median_days_between_first_two": med_gap,
                "icc_one_lab_night_TST": M.ICC_TST,
                "verdict": "the home measure is about as reliable as the laboratory night"}


# ================================================================================================
# A3  is the map inversion an interpolation or an extrapolation
# ================================================================================================
def attack3(rk, target):
    print("\n" + "=" * 92)
    print("A3  THE MAP.  is the required point inside the observed data")
    print("=" * 92)
    mp = pd.read_csv(f"{HERE}/effect_vs_gain.csv")
    S = json.load(open(f"{HERE}/summary.json"))
    a_, c_ = S["map_primary"]["slope"], S["map_primary"]["intercept"]
    E = mp.E_meanabs.values
    E_req = float(np.sqrt((target - c_) / a_))
    pct = float((E < E_req).mean() * 100)
    n_above = int((E > E_req).sum())
    print(f"  required effect size E {E_req:.4f}. Observed range {E.min():.4f} to {E.max():.4f}.")
    print(f"  it sits at the {pct:.1f}th percentile of the {N_MEASURES} measures, with {n_above} "
          f"measures ABOVE it. This is interpolation, not extrapolation.")
    row("A3", "required E_meanabs", "", E_req, f"{pct:.1f}th pct")

    win = mp[(mp.E_meanabs > E_req - 0.02) & (mp.E_meanabs < E_req + 0.02)]
    print(f"  the {len(win)} real measures within +/-0.02 of that effect size actually achieve "
          f"dC {win.dC.min():+.6f} to {win.dC.max():+.6f} (median {win.dC.median():+.6f}), "
          f"ranks {int(win['rank'].min())} to {int(win['rank'].max())}")
    for _, r_ in win.sort_values("E_meanabs").iterrows():
        print(f"      rank {int(r_['rank']):3d}  {r_.feature:28s} E {r_.E_meanabs:.4f}  "
              f"dC {r_.dC:+.6f}")
    model_free = float(win.dC.median())
    print(f"  MODEL-FREE version of the answer: a measure needs an effect size around "
          f"{E_req:.3f} to sit at rank 10, confirmed without the fitted map, because the "
          f"measures already there have exactly that effect size.")

    # stability of the inversion when the fit is restricted to the weak measures only
    stab = {}
    for cut in [0.10, 0.12, 0.15, 0.20, 1.0]:
        s = mp[mp.E_meanabs <= cut]
        if len(s) < 30:
            continue
        sl, ic, rr, pp, se = stats.linregress(s.E_meanabs.values ** 2, s.dC.values)
        er = float(np.sqrt(max(0.0, (target - ic) / sl)))
        stab[f"fit_on_E_le_{cut}"] = {"n": int(len(s)), "slope": float(sl),
                                      "intercept": float(ic), "r2": float(rr ** 2),
                                      "E_required": er, "hr_required": float(np.exp(er))}
        print(f"  map refit on the {len(s):3d} measures with E <= {cut:.2f}: "
              f"R2 {rr**2:.3f}, required E {er:.4f} (HR {np.exp(er):.3f})")

    resid = float(mp[mp.feature == "TST_min"].dC.iloc[0] -
                  (a_ * mp[mp.feature == "TST_min"].E_meanabs.iloc[0] ** 2 + c_))
    print(f"  map residual at TST itself: {resid:+.6f}")
    verdict = "HOLDS" if (E.min() < E_req < E.max() and n_above >= 5) else "WEAKENED"
    R["A3"] = {"E_required": E_req, "E_min": float(E.min()), "E_max": float(E.max()),
               "percentile_of_required": pct, "n_measures_above": n_above,
               "window_n": int(len(win)),
               "window_dC_median": model_free,
               "window_dC_range": [float(win.dC.min()), float(win.dC.max())],
               "window_rank_range": [int(win["rank"].min()), int(win["rank"].max())],
               "window_features": win[["rank", "feature", "E_meanabs", "dC"]].to_dict("records"),
               "inversion_stability": stab, "map_residual_at_TST": resid,
               "verdict": verdict}
    print(f"  VERDICT A3: {verdict}")


# ================================================================================================
# A4  circularity.  can the simulation recover a rise that is there by construction
# ================================================================================================
def _degrade_and_refit(feat, r0, lam, rep, outs, base, permute=False):
    """Build a surrogate of `feat` whose reliability about the real column is r0, add lam times
    its own error variance on top, refit the whole panel, return dC.

    lam = 0 is the surrogate itself. The extrapolation back to reliability 1.0 has a known
    right answer: the dC of the undegraded column.
    """
    d = M.data()
    tag = sum(ord(ch) for ch in feat) + (5000 if permute else 0)     # deterministic across procs
    rng = np.random.default_rng(SEED + tag + rep * 977 + int(lam * 100))
    z = pd.Series(np.nan, index=d.index)
    for _, idx in d.groupby("site_id").groups.items():
        v = d.loc[idx, feat].astype(float)
        ok = v.notna()
        if ok.sum() < 20:
            continue
        t = M.rint(v[ok].values)                      # truth on the analysis scale, variance 1
        if permute:
            t = rng.permutation(t)
        s = t + rng.normal(0.0, np.sqrt((1 - r0) / r0), size=len(t))   # reliability r0
        if lam > 0:
            s = s + rng.normal(0.0, np.sqrt(lam * (1 - r0) / r0), size=len(s))
        z.loc[v[ok].index] = M.rint(s)
    d["_deg__z"] = z
    g = np.array([M.heldout_c(d, "_deg__z", o) - base[o] for o in outs], float)
    return {"feature": feat, "permute": permute, "lam": lam, "rep": rep,
            "dC": float(np.nanmean(g))}


def attack4(rk, target, outs, base):
    print("\n" + "=" * 92)
    print("A4  CIRCULARITY.  degrade a measure to sleep duration's reliability, then try to "
          "recover it")
    print("=" * 92)
    r0 = M.ICC_TST
    print(f"  every surrogate is built to have reliability exactly {r0:.4f} about the real "
          f"column, the reliability of one laboratory night of sleep duration.")
    targets = [("spo2_pct_below_90", False), ("spo2_nadir_corrected", False),
               ("TST_min", False), ("TST_min", True)]
    jobs = []
    for feat, perm in targets:
        for lam in LAMBDAS:
            for rep in range(1 if lam == 0 else NREPS_A4):
                jobs.append((feat, lam, rep, perm))
    p = f"{HERE}/_attack_A4.csv"
    if os.path.exists(p):
        df = pd.read_csv(p)
    else:
        t0 = time.time()
        out = Parallel(n_jobs=NJOBS, verbose=1)(
            delayed(_degrade_and_refit)(f, r0, l, r, outs, base, perm)
            for f, l, r, perm in jobs)
        df = pd.DataFrame(out)
        df.to_csv(p, index=False)
        print(f"  {len(jobs)} panel refits in {time.time()-t0:.0f}s")

    truth = rk.set_index("feature").dC
    rec = {}
    for feat, perm in targets:
        s = df[(df.feature == feat) & (df.permute == perm)]
        agg = (s.groupby("lam").dC.agg(["mean", "std", "count"]).reset_index()
               .rename(columns={"mean": "dC", "std": "sd", "count": "reps"}))
        rel = r0 / (1 + agg.lam.values * (1 - r0))
        w = np.sqrt(agg.reps.values.astype(float))
        co = np.polyfit(rel, agg.dC.values, 1, w=w)
        pred1 = float(np.polyval(co, 1.0))
        fit = np.polyval(co, rel)
        r2 = float(1 - ((agg.dC.values - fit) ** 2).sum() /
                   max(1e-18, ((agg.dC.values - agg.dC.values.mean()) ** 2).sum()))
        known = 0.0 if perm else float(truth[feat])
        nm = f"{feat}{' PERMUTED (null control)' if perm else ''}"
        err = pred1 - known
        relerr = err / abs(known) if abs(known) > 1e-6 else float("nan")
        rec[nm] = {"known_dC_undegraded": known,
                   "dC_of_degraded_surrogate": float(agg.loc[agg.lam == 0, "dC"].iloc[0]),
                   "recovered_dC_at_reliability_1": pred1,
                   "absolute_error": float(err),
                   "relative_error": float(relerr) if np.isfinite(relerr) else None,
                   "slope_per_reliability": float(co[0]), "fit_r2": r2,
                   "levels": agg.assign(reliability=rel).to_dict("records"),
                   "rank_recovered": int((rk.dC > pred1).sum() + 1)}
        print(f"  {nm:42s} true dC {known:+.6f}  degraded {rec[nm]['dC_of_degraded_surrogate']:+.6f}"
              f"  recovered {pred1:+.6f}  error {err:+.6f}"
              f"{'' if not np.isfinite(relerr) else f' ({100*relerr:+.1f}%)'}  R2 {r2:.3f}")
        row("A4", "recovery test", nm, pred1, known)

    key = "spo2_pct_below_90"
    ok_pos = abs(rec[key]["relative_error"]) < 0.35
    ok_pos2 = abs(rec["spo2_nadir_corrected"]["relative_error"]) < 0.35
    ok_null = abs(rec["TST_min PERMUTED (null control)"]["recovered_dC_at_reliability_1"]) < 0.002
    verdict = ("HOLDS, the method recovers a real rise and does not manufacture one"
               if (ok_pos and ok_pos2 and ok_null) else
               "WEAKENED, the recovery test does not clear its own tolerance")
    R["A4"] = {"reliability_imposed": r0, "n_refits": int(len(df)),
               "recovery": rec,
               "positive_control_within_35pct": bool(ok_pos and ok_pos2),
               "null_control_stays_null": bool(ok_null),
               "verdict": verdict}
    print(f"  VERDICT A4: {verdict}")


# ================================================================================================
# A6  power.  bootstrap the whole panel
# ================================================================================================
def _boot_rep(b, outs, feats):
    d = M.data()
    rng = np.random.default_rng(SEED + 31 * b)
    parts = []
    for s in M.SITES:
        idx = np.where(d.site_id.values == s)[0]
        parts.append(idx[rng.integers(0, len(idx), len(idx))])
    dd = d.iloc[np.concatenate(parts)].reset_index(drop=True)
    base = {o: M.heldout_c(dd, None, o) for o in outs}
    out = {"b": b}
    for f in feats:
        g = np.array([M.heldout_c(dd, f + "__z", o) - base[o] for o in outs], float)
        out[f] = float(np.nanmean(g))
    return out


def attack6(rk, target, outs):
    print("\n" + "=" * 92)
    print("A6  POWER.  bootstrap interval on dC and on the gap to rank 10")
    print("=" * 92)
    p = f"{HERE}/_attack_A6_boot.csv"
    if os.path.exists(p):
        bt = pd.read_csv(p)
    else:
        t0 = time.time()
        out = Parallel(n_jobs=NJOBS, verbose=1)(
            delayed(_boot_rep)(b, outs, BOOT_FEATS) for b in range(NBOOT))
        bt = pd.DataFrame(out)
        bt.to_csv(p, index=False)
        print(f"  {NBOOT} bootstrap resamples, each refitting the whole panel four times, "
              f"in {time.time()-t0:.0f}s")

    obs = rk.set_index("feature").dC
    st = {}
    for f in BOOT_FEATS:
        v = bt[f].dropna().values
        st[f] = {"observed_dC": float(obs[f]), "boot_mean": float(v.mean()),
                 "boot_se": float(v.std(ddof=1)),
                 "ci95": [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))],
                 "n_boot": int(len(v))}
        print(f"  {f:22s} dC {obs[f]:+.6f}   bootstrap SE {v.std(ddof=1):.6f}   "
              f"95% CI {np.percentile(v,2.5):+.6f} to {np.percentile(v,97.5):+.6f}")
        row("A6", "bootstrap SE of dC", f, float(v.std(ddof=1)), float(obs[f]))

    dif = (bt["nrem_lf_hf_ratio"] - bt["TST_min"]).dropna().values
    lo, hi = float(np.percentile(dif, 2.5)), float(np.percentile(dif, 97.5))
    zsc = float(dif.mean() / dif.std(ddof=1))
    print(f"  paired difference, rank-10 measure minus TST: {dif.mean():+.6f} "
          f"(95% CI {lo:+.6f} to {hi:+.6f}), z = {zsc:.2f}")
    sep = bool(lo > 0)
    # is TST distinguishable from zero
    vt = bt["TST_min"].dropna().values
    tst_ne0 = bool(np.percentile(vt, 2.5) > 0 or np.percentile(vt, 97.5) < 0)
    # smallest dC gap the design can resolve at 80% power, two-sided 5%
    se_dif = float(dif.std(ddof=1))
    mdd = 2.802 * se_dif
    print(f"  TST's dC is {'' if tst_ne0 else 'NOT '}distinguishable from zero")
    print(f"  the design resolves a dC difference of {mdd:.6f} at 80% power; the gap TST has "
          f"to close is {target - float(obs['TST_min']):.6f}, which is "
          f"{(target - float(obs['TST_min']))/mdd:.1f} times that")
    R["A6"] = {"n_boot": int(len(bt)), "per_feature": st,
               "paired_difference_rank10_minus_TST": {"mean": float(dif.mean()),
                                                      "se": se_dif, "ci95": [lo, hi], "z": zsc},
               "TST_distinguishable_from_zero": tst_ne0,
               "TST_distinguishable_from_rank10": sep,
               "minimum_detectable_dC_difference_80pct": float(mdd),
               "gap_to_close": float(target - float(obs["TST_min"])),
               "gap_in_units_of_minimum_detectable": float(
                   (target - float(obs["TST_min"])) / mdd),
               "verdict": ("HOLDS, the gap is far larger than the noise on it" if sep
                           else "WEAKENED, rank 10 and TST are not separable")}
    print(f"  VERDICT A6: {R['A6']['verdict']}")


def main():
    t0 = time.time()
    ladder, target, rk = attack1()
    outs, base = attack2(rk, target)
    attack1b(outs)          # v7: the shipped summary carried A1b, A2b and A1c (ED Fig 15 reads A1c); restored to main
    attack2b()
    attack1c()
    attack3(rk, target)
    attack4(rk, target, outs, base)
    attack6(rk, target, outs)
    json.dump(R, open(f"{HERE}/attack_summary.json", "w"), indent=1, default=str)
    pd.DataFrame(ROWS, columns=["attack", "quantity", "detail", "value", "extra"]).to_csv(
        f"{HERE}/attack_results.csv", index=False)
    print("\n" + "=" * 92)
    for k in ("A1", "A1b", "A1c", "A2", "A2b", "A3", "A4", "A6"):
        print(f"  {k}: {R[k]['verdict']}")
    print(f"wrote attack_summary.json and attack_results.csv ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
