"""
What it would take for total sleep time to reach the top 10 of the 197 measures.

Total sleep time sits at rank 190 of 197 in the frozen comparison, with a held-out concordance
gain of -0.001016, which is what a randomly permuted column costs. The question is what would
have to change for it to enter the top 10, and the question has two readings that must be
answered separately.

Reading one, measurement. One laboratory night measures habitual sleep duration with an ICC of
0.482 in this project's own repeat studies. Averaging k nights at home raises the reliability by
Spearman-Brown. How many nights would it take.

Reading two, biology. Independently of measurement, how much larger would the underlying
association with disease have to be.

The two readings are not the same, and the arithmetic that separates them is short. Measurement
error attenuates an association toward zero by a factor that is bounded by one. Perfect
measurement therefore MULTIPLIES an observed effect by a finite constant, 1/sqrt(0.482) = 1.44 on
the log-hazard scale. If the required multiple exceeds that constant, no number of home nights
suffices and the answer is a statement about the biology, not the measurement.

Five things are computed.

  Stage A, the ranking is reproduced from the frozen inputs so that the target, the effect-size
  map and the simulation all sit on one internally consistent scale. Specification copied from
  numbers/run_ranking_v2.py: rank-inverse-normal within site, natural cubic age spline with 4 df
  plus sex, Cox with penalizer 0.01, fitted at one hospital and scored at the other, both
  directions averaged, gain taken against an age-and-sex model fitted the same way, then averaged
  over the 48 outcomes. The frozen ranking_v2.csv is reported alongside as the authority.

  Stage B, a per-measure summary effect for all 197 measures on the analysis cohort. Standard
  model: site-stratified UNPENALIZED Cox, natural cubic age spline with 4 df of which ONE column
  is dropped so the design is full rank, plus sex, exposure rank-inverse-normal within site, so
  the coefficient is per 1 SD. An unpenalized fit on a rank-deficient design returns zeros
  without failing, so the age coefficients are checked to be non-zero. Each measure is summarised
  over the 48 outcomes as a root-mean-square log hazard ratio, both raw and with the sampling
  noise floor removed.

  Stage C, the map from effect size to concordance gain is fitted across the 197 measures in
  nine forms and inverted at the rank-10 gain. The best-fitting form is also the one with a
  mechanism behind it: adding any column to a held-out model costs a fixed amount in
  overfitting, which is the intercept, and buys an amount that scales with the SQUARE of its
  coefficient, because that is how a coefficient enters the information in a partial likelihood.

  Stage D, disattenuation and Spearman-Brown, to turn a required effect size into a required
  reliability and then into a number of home nights. Two attenuation forms are run, the square
  root of reliability and the classical regression form, because the second is more forgiving
  and still has to be beaten.

  Stage E, the direct simulation, which needs none of the above. Measurement noise is ADDED to
  the observed values at several levels, the ranking is refit end to end at each level, and the
  curve of concordance gain against noise is extrapolated back to the noise-free case. This is
  SIMEX. Noise is added two ways: in the measure's own units, which is right for sleep duration,
  where a second night differs by an additive number of minutes, and on the rank-inverse-normal
  scale the Cox model actually sees, which is right for T90, where additive noise in percentage
  points is not what a second night looks like on a zero-inflated skewed measure. Four
  extrapolants are reported. The primary one is linear in RELIABILITY, which is what the
  quadratic concordance map and square-root attenuation together imply, and which fits the T90
  curve with an R2 of 0.997. T90 is the positive control: if its gain rises steeply as
  measurement improves while sleep duration's does not move, the difference between them is
  signal, not noise.

Writes results.csv and summary.json in this directory. verify.py re-derives every headline
number from those files independently.
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
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from patsy import dmatrix
from scipy import stats
from scipy.optimize import curve_fit
import sys as _sys_v8; _sys_v8.path.insert(0, paths.ANALYSIS_DIR)  # v8 sweep 2026-09-12
from cohort_spec import N_MEASURES, N_RANKED_OUTCOMES, SENSITIVITY_ONLY, RANKING_DROP_SUFFIX  # v8 sweep 2026-09-12

warnings.filterwarnings("ignore")

T90ROOT = paths.T90_ROOT
SVROOT = paths.SV_ROOT
HERE = f"{SVROOT}/tst_ceiling"
PKL = f"{HERE}/rank_data.pkl"
MASTER = (f"{paths.TABLES_DIR}/"
          "master_cohort.csv")
CIRCULAR = {"insomnia", "rls_plmd", "nocturia", "epilepsy"}
SITES = ("I0002", "I0006")
NJOBS = 7

# reliability of one laboratory night, from this project's own repeat studies in the same
# untreated patients, data/untreated_wide.parquet
ICC_TST = 0.4822461259107412   # n = 415 device-confirmed PAP-free pairs
ICC_T90 = 0.4208377980350933   # n = 410 of the same pairs
ICC_TST_ALL = 0.414            # all 1,664 clean pairs, carried for sensitivity
ICC_T90_SHHS = 0.476           # 4,311 people, two polysomnograms about five years apart
REPEATS = f"{SVROOT}/data/untreated_wide.parquet"

# noise levels for the simulation, as multiples of the measurement-error variance already
# present. lam = -1 is the noise-free case, which is what the extrapolation targets.
LAMBDAS = [0.0, 0.5, 1.0, 1.5, 2.0, 3.0]
NREPS = 12
SEED = 20260807


# --------------------------------------------------------------------------------------------
# data
# --------------------------------------------------------------------------------------------
def rint(v):
    r = stats.rankdata(v, method="average")
    return stats.norm.ppf((r - 0.375) / (len(r) + 0.25))


def build():
    """Rebuild the ranking data frame exactly as numbers/run_ranking_v2.py does."""
    sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
    from disease_definitions import DISEASES
    from cohort_spec import drop_split_nights_if_full, ranked_outcome_keys   # v8.1

    base = pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet")
    _tbl = base
    base = base[(base.fu_valid == 1) & base.spo2_pct_below_90.notna() & (base.oximetry_bad == 0)]
    base = drop_split_nights_if_full(base)   # v8.1: the ranking family runs on the full diagnostic nights
    head = pd.read_csv(MASTER, nrows=0).columns.tolist()
    feat = [c for c in head if c != "BDSPPatientID" and not c.endswith(("_first_date", "_date"))]
    raw = pd.read_csv(MASTER, usecols=["BDSPPatientID"] + feat, low_memory=False)
    left = base[["BDSPPatientID", "site_id", "AgeAtVisit", "sex"] +
                [c for c in base.columns if c.endswith(("_incident", "_years", "_prevalent"))]]
    raw = raw.drop(columns=[c for c in raw.columns if c in set(left.columns) - {"BDSPPatientID"}])
    feat = [c for c in feat if c in raw.columns]
    d = left.merge(raw, on="BDSPPatientID", how="left")
    d["male"] = (d.sex.astype(str).str.upper().str[0] == "M").astype(int)

    num = [c for c in feat if pd.api.types.is_numeric_dtype(d[c])]
    num = [c for c in num if d[c].notna().mean() >= 0.60 and d[c].nunique() > 10]
    num = [c for c in num if c not in SENSITIVITY_ONLY and not c.endswith(RANKING_DROP_SUFFIX)]   # v8: the primary ranking's candidate rule (build_ranking_v3)
    num = [c for c in num if c not in {"male", "AgeAtVisit"}
           and not c.startswith(("obs_", "follow_"))]

    sp = dmatrix("cr(a, df=4) - 1", {"a": d.AgeAtVisit.values}, return_type="dataframe")
    sp.columns = [f"age_s{i}" for i in range(4)]
    sp.index = d.index
    for c in sp.columns:
        d[c] = sp[c]

    for f in num:
        z = pd.Series(np.nan, index=d.index)
        for _, idx in d.groupby("site_id").groups.items():
            v = d.loc[idx, f]
            ok = v.notna()
            if ok.sum() < 20:
                continue
            z.loc[v[ok].index] = rint(v[ok])
        d[f + "__z"] = z

    out = [k for k in DISEASES if k not in CIRCULAR and f"{k}_incident" in d.columns]
    out = ranked_outcome_keys(_tbl, out) + ["death"]      # v8.1: the ranked set is fixed on the all-nights cohort
    d.to_pickle(PKL)
    json.dump({"num": num, "OUT": out}, open(f"{HERE}/spec.json", "w"))
    return d, num, out


_CACHE = {}


def data():
    if "d" not in _CACHE:
        _CACHE["d"] = pd.read_pickle(PKL)
    return _CACHE["d"]


ADJ4 = [f"age_s{i}" for i in range(4)] + ["male"]          # ranking, penalized, full basis
ADJ3 = [f"age_s{i}" for i in range(3)] + ["male"]          # effects, unpenalized, one dropped


# --------------------------------------------------------------------------------------------
# Stage A, the ranking
# --------------------------------------------------------------------------------------------
def heldout_c(d, col, out):
    """Cross-site held-out concordance, both directions averaged. Ranking specification."""
    g = d[(d[f"{out}_prevalent"] == 0) & (d[f"{out}_years"] > 0) & d[f"{out}_years"].notna()]
    cols = ADJ4 + ([col] if col else [])
    res = []
    for tr, te in (SITES, SITES[::-1]):
        a = g[g.site_id == tr][cols + [f"{out}_years", f"{out}_incident"]].dropna()
        b = g[g.site_id == te][cols + [f"{out}_years", f"{out}_incident"]].dropna()
        if a[f"{out}_incident"].sum() < 30 or b[f"{out}_incident"].sum() < 30:
            res.append(np.nan)
            continue
        try:
            fit = CoxPHFitter(penalizer=0.01).fit(a, f"{out}_years", f"{out}_incident")
            p = fit.predict_partial_hazard(b)
            res.append(concordance_index(b[f"{out}_years"], -p, b[f"{out}_incident"]))
        except Exception:
            res.append(np.nan)
    return float(np.nanmean(res))


def rank_one(feat, outs):
    d = data()
    return {o: heldout_c(d, feat + "__z", o) for o in outs}


def stage_a(num, outs):
    p = f"{HERE}/_stageA_ranking.csv"
    if os.path.exists(p):
        return pd.read_csv(p)
    t0 = time.time()
    d = data()
    base = {o: heldout_c(d, None, o) for o in outs}
    json.dump(base, open(f"{HERE}/_baselineC.json", "w"))
    print(f"  baseline held-out C {np.mean(list(base.values())):.4f}  ({time.time()-t0:.0f}s)")
    res = Parallel(n_jobs=NJOBS, verbose=1)(delayed(rank_one)(f, outs) for f in num)
    rows = []
    for f, r in zip(num, res):
        g = np.array([r[o] - base[o] for o in outs], float)
        rows.append({"feature": f, "mC": float(np.nanmean([r[o] for o in outs])),
                     "dC": float(np.nanmean(g)), "worst": float(np.nanmin(g)),
                     "npos": int((g > 0).sum()), "nout": int(np.isfinite(g).sum())})
    rk = pd.DataFrame(rows).sort_values("dC", ascending=False).reset_index(drop=True)
    rk.insert(0, "rank", range(1, len(rk) + 1))
    rk.to_csv(p, index=False)
    print(f"  stage A done in {time.time()-t0:.0f}s")
    return rk


# --------------------------------------------------------------------------------------------
# Stage B, per-measure summary effect
# --------------------------------------------------------------------------------------------
def effect_one(feat, outs):
    d = data()
    col = feat + "__z"
    betas, ses, spline_ok = [], [], []
    for o in outs:
        g = d[(d[f"{o}_prevalent"] == 0) & (d[f"{o}_years"] > 0) & d[f"{o}_years"].notna()]
        dd = g[[col] + ADJ3 + [f"{o}_years", f"{o}_incident", "site_id"]].dropna()
        if dd[f"{o}_incident"].sum() < 60:
            continue
        try:
            fit = CoxPHFitter(penalizer=0.0).fit(dd, f"{o}_years", f"{o}_incident",
                                                 strata=["site_id"])
            r = fit.summary.loc[col]
            betas.append(float(r["coef"]))
            ses.append(float(r["se(coef)"]))
            spline_ok.append(bool(np.abs(fit.params_[[f"age_s{i}" for i in range(3)]]).max() > 1e-8))
        except Exception:
            continue
    b = np.array(betas)
    s = np.array(ses)
    if len(b) == 0:
        return None
    raw = float(np.sqrt((b ** 2).mean()))
    adj = float(np.sqrt(max(0.0, (b ** 2).mean() - (s ** 2).mean())))
    return {"feature": feat, "n_out": len(b), "E_rms": raw, "E_rms_adj": adj,
            "E_meanabs": float(np.abs(b).mean()), "max_hr": float(np.exp(np.abs(b).max())),
            "mean_chi2": float(((b / s) ** 2).mean()),
            "spline_nonzero_frac": float(np.mean(spline_ok))}


def stage_b(num, outs):
    p = f"{HERE}/_stageB_effects.csv"
    if os.path.exists(p):
        return pd.read_csv(p)
    t0 = time.time()
    res = Parallel(n_jobs=NJOBS, verbose=1)(delayed(effect_one)(f, outs) for f in num)
    df = pd.DataFrame([r for r in res if r])
    df.to_csv(p, index=False)
    print(f"  stage B done in {time.time()-t0:.0f}s")
    return df


# --------------------------------------------------------------------------------------------
# Stage E, the simulation
# --------------------------------------------------------------------------------------------
def simex_one(feat, icc, lam, rep, outs, scale):
    """Add measurement noise of lam times the error variance already present, then refit.

    scale "raw"  adds Gaussian noise in the measure's own units, which is the right model for
                 sleep duration, where night-to-night variation is additive in minutes.
    scale "rint" adds it on the rank-inverse-normal scale the Cox model actually sees, which is
                 the right model for a zero-inflated skewed measure such as T90, where additive
                 noise in percentage points is not what a second night looks like.
    """
    d = data()   # process-local, "_sim__z" is overwritten on every call
    rng = np.random.default_rng(SEED + rep * 1000 + int(lam * 100) + (7 if scale == "rint" else 0))
    x = d[feat].astype(float)
    z = pd.Series(np.nan, index=d.index)
    for _, idx in d.groupby("site_id").groups.items():
        v = x.loc[idx]
        ok = v.notna()
        if ok.sum() < 20:
            continue
        vv = v[ok].values
        if scale == "rint":
            vv = rint(vv)                       # variance 1 by construction
            sd = 1.0
        else:
            sd = vv.std(ddof=1)
        if lam > 0:
            vv = vv + rng.normal(0.0, sd * np.sqrt(lam * (1 - icc)), size=len(vv))
        z.loc[v[ok].index] = rint(vv)
    d["_sim__z"] = z
    if "base" not in _CACHE:
        _CACHE["base"] = json.load(open(f"{HERE}/_baselineC.json"))
    base = _CACHE["base"]
    g = np.array([heldout_c(d, "_sim__z", o) - base[o] for o in outs], float)
    return {"feature": feat, "scale": scale, "icc": icc, "lam": lam, "rep": rep,
            "dC": float(np.nanmean(g))}


def repeat_iccs():
    """ICC of one night, raw and on the rank-inverse-normal analysis scale, from the repeats."""
    p = f"{HERE}/_iccs.json"
    if os.path.exists(p):
        return json.load(open(p))
    w = pd.read_parquet(REPEATS)

    def icc(a, b):
        a, b = np.asarray(a, float), np.asarray(b, float)
        m = np.isfinite(a) & np.isfinite(b)
        a, b = a[m], b[m]
        x = np.c_[a, b]
        n = len(a)
        msb = 2 * ((x.mean(1) - x.mean()) ** 2).sum() / (n - 1)
        msw = ((x - x.mean(1, keepdims=True)) ** 2).sum() / n
        return float((msb - msw) / (msb + msw)), int(n)

    out = {}
    for nm, c1, c2 in (("TST", "tst_min_n1", "tst_min_n2"), ("T90", "t90_n1", "t90_n2")):
        a, b = w[c1].astype(float), w[c2].astype(float)
        m = a.notna() & b.notna()
        raw, n = icc(a[m], b[m])
        both = rint(np.r_[a[m].values, b[m].values])
        rt, _ = icc(both[:m.sum()], both[m.sum():])
        out[nm] = {"icc_raw": raw, "icc_rint": rt, "n_pairs": n}
    json.dump(out, open(p, "w"), indent=1)
    return out


def stage_e(outs, scale, iccs):
    p = f"{HERE}/_stageE_simex_{scale}.csv"
    if os.path.exists(p):
        return pd.read_csv(p)
    assert os.path.exists(f"{HERE}/_baselineC.json"), "run stage A first"
    t0 = time.time()
    key = "icc_raw" if scale == "raw" else "icc_rint"
    jobs = []
    for feat, nm in (("TST_min", "TST"), ("spo2_pct_below_90", "T90")):
        for lam in LAMBDAS:
            for rep in range(1 if lam == 0 else NREPS):
                jobs.append((feat, iccs[nm][key], lam, rep))
    res = Parallel(n_jobs=NJOBS, verbose=1)(
        delayed(simex_one)(f, i, l, r, outs, scale) for f, i, l, r in jobs)
    df = pd.DataFrame(res)
    df.to_csv(p, index=False)
    print(f"  stage E ({scale}) done in {time.time()-t0:.0f}s, {len(jobs)} refits of the panel")
    return df


def extrapolate(agg, icc):
    """Read dC off at zero measurement error, four ways.

    quadratic and linear in the noise multiplier are the textbook SIMEX extrapolants. The
    rational-linear one is the third textbook choice and is known to be unstable, so it is
    reported but not used. The fourth is the one with a mechanism: the concordance map is close
    to quadratic in the log hazard ratio, and the log hazard ratio is attenuated by the square
    root of reliability, so the gain should be close to LINEAR IN RELIABILITY. That form is
    fitted on reliability and read off at 1.0.
    """
    lam = agg.lam.values
    y = agg["dC"].values
    w = np.sqrt(agg["reps"].values.astype(float))
    rel = icc / (1 + lam * (1 - icc))
    q = np.polyfit(lam, y, 2, w=w)
    li = np.polyfit(lam, y, 1, w=w)
    lr = np.polyfit(rel, y, 1, w=w)
    try:
        rl, _ = curve_fit(lambda l, g0, g1, g2: (g0 + g1 * l) / (1 + g2 * l), lam, y,
                          p0=[y[0], 0.0, 1.0], maxfev=40000)
        rl_m1 = float((rl[0] - rl[1]) / (1 - rl[2]))
    except Exception:
        rl_m1 = float("nan")

    def r2(pred):
        return float(1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum())

    return {
        "quadratic_in_lambda": float(np.polyval(q, -1.0)),
        "linear_in_lambda": float(np.polyval(li, -1.0)),
        "rational_in_lambda": rl_m1,
        "linear_in_reliability": float(np.polyval(lr, 1.0)),
        "fit_r2_quadratic": r2(np.polyval(q, lam)),
        "fit_r2_linear": r2(np.polyval(li, lam)),
        "fit_r2_linear_in_reliability": r2(np.polyval(lr, rel)),
        "slope_per_lambda": float(li[0]),
        "slope_per_reliability": float(lr[0]),
        "_coef_quad": [float(v) for v in q],
        "_coef_rel": [float(v) for v in lr],
    }


# --------------------------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------------------------
def spearman_brown_k(target_r, icc):
    """Nights needed for a k-night average to reach reliability target_r."""
    if target_r >= 1.0:
        return np.inf
    return target_r * (1 - icc) / (icc * (1 - target_r))


def reliability_k(k, icc):
    return k * icc / (1 + (k - 1) * icc)


def lam_for_reliability(r, icc):
    """Noise multiplier that yields reliability r, given the observed icc. r = icc -> 0."""
    return (icc / r - 1.0) / (1.0 - icc)


def main():
    t0 = time.time()
    if os.path.exists(PKL) and os.path.exists(f"{HERE}/spec.json"):
        spec = json.load(open(f"{HERE}/spec.json"))
        num, outs = spec["num"], spec["OUT"]
        d = data()
    else:
        d, num, outs = build()
    print(f"cohort {len(d):,}   measures {len(num)}   outcomes {len(outs)}")
    R = {"cohort_n": int(len(d)), "n_measures": int(len(num)), "n_outcomes": int(len(outs)),
         "sites": list(SITES)}

    # ---- Step 1, the target ------------------------------------------------------------------
    frozen = pd.read_csv(f"{paths.NUMBERS_DIR}/ranking_v3.csv", comment="#")   # v8: the current ranking is the authority
    assert len(frozen) == N_MEASURES, len(frozen)
    fz = frozen.set_index("rank")
    tst_fz = frozen[frozen.feature == "TST_min"].iloc[0]
    target_fz = float(fz.loc[10, "dC"])
    R["frozen"] = {
        "n_rows": int(len(frozen)),
        "dC_rank5": float(fz.loc[5, "dC"]), "feature_rank5": str(fz.loc[5, "feature"]),
        "dC_rank10": target_fz, "feature_rank10": str(fz.loc[10, "feature"]),
        "dC_rank20": float(fz.loc[20, "dC"]), "feature_rank20": str(fz.loc[20, "feature"]),
        "dC_rank1": float(fz.loc[1, "dC"]), "feature_rank1": str(fz.loc[1, "feature"]),
        "TST_rank": int(tst_fz["rank"]), "TST_dC": float(tst_fz.dC), "TST_mC": float(tst_fz.mC),
        "T90_rank": int(frozen[frozen.feature == "spo2_pct_below_90"].iloc[0]["rank"]),
        "T90_dC": float(frozen[frozen.feature == "spo2_pct_below_90"].iloc[0].dC),
        "gain_TST_must_add": target_fz - float(tst_fz.dC),
    }
    print(f"\nSTEP 1  frozen target: rank 10 is {fz.loc[10,'feature']} at dC {target_fz:+.6f}")
    print(f"        rank 5 {fz.loc[5,'dC']:+.6f}   rank 20 {fz.loc[20,'dC']:+.6f}")
    print(f"        TST_min rank {int(tst_fz['rank'])} dC {tst_fz.dC:+.6f}, "
          f"must add {target_fz - tst_fz.dC:+.6f}")

    print(f"\nSTAGE A  reproducing the {N_MEASURES}-measure ranking")
    rk = stage_a(num, outs)
    rk.to_csv(f"{HERE}/ranking_reproduced.csv", index=False)
    rkx = rk.set_index("rank")
    target = float(rkx.loc[10, "dC"])
    tst_r = rk[rk.feature == "TST_min"].iloc[0]
    t90_r = rk[rk.feature == "spo2_pct_below_90"].iloc[0]
    m = frozen.set_index("feature").dC.reindex(rk.feature).values
    R["reproduced"] = {
        "dC_rank5": float(rkx.loc[5, "dC"]), "dC_rank10": target,
        "dC_rank20": float(rkx.loc[20, "dC"]),
        "feature_rank10": str(rkx.loc[10, "feature"]),
        "TST_rank": int(tst_r["rank"]), "TST_dC": float(tst_r.dC),
        "T90_rank": int(t90_r["rank"]), "T90_dC": float(t90_r.dC),
        "gain_TST_must_add": target - float(tst_r.dC),
        "spearman_vs_frozen": float(stats.spearmanr(rk.dC.values, m, nan_policy="omit")[0]),
        "pearson_vs_frozen": float(np.corrcoef(rk.dC.values, m)[0, 1]),
    }
    print(f"  reproduced: rank 10 dC {target:+.6f}, TST rank {int(tst_r['rank'])} "
          f"dC {tst_r.dC:+.6f}, T90 rank {int(t90_r['rank'])} dC {t90_r.dC:+.6f}")
    print(f"  agreement with frozen ranking: Spearman "
          f"{R['reproduced']['spearman_vs_frozen']:.3f}, Pearson "
          f"{R['reproduced']['pearson_vs_frozen']:.3f}")

    # ---- Step 2, the effect-size map ---------------------------------------------------------
    print("\nSTAGE B  per-measure summary effect, site-stratified unpenalized Cox")
    ef = stage_b(num, outs)
    print(f"  age-spline coefficients non-zero in "
          f"{100*ef.spline_nonzero_frac.mean():.1f}% of fits (rank-deficiency guard)")
    mp = rk.merge(ef, on="feature", how="inner")
    mp.to_csv(f"{HERE}/effect_vs_gain.csv", index=False)
    R["stageB_spline_nonzero_frac"] = float(ef.spline_nonzero_frac.mean())

    fits = {}
    for xcol in ["E_rms_adj", "E_rms", "E_meanabs"]:
        x, y = mp[xcol].values, mp.dC.values
        ok = np.isfinite(x) & np.isfinite(y)
        for form, xx in (("linear", x[ok]), ("quadratic", x[ok] ** 2), ("sqrt", np.sqrt(x[ok]))):
            sl, ic, r, p, se = stats.linregress(xx, y[ok])
            fits[f"{xcol}|{form}"] = {"slope": float(sl), "intercept": float(ic),
                                      "r2": float(r ** 2), "p": float(p), "se_slope": float(se),
                                      "n": int(ok.sum())}
    best = max(fits, key=lambda k: fits[k]["r2"])
    R["map_fits"] = fits
    R["map_best"] = best
    print(f"  best map: {best}  R2 {fits[best]['r2']:.3f}")
    for k in sorted(fits, key=lambda k: -fits[k]["r2"])[:6]:
        print(f"     {k:22s} R2 {fits[k]['r2']:.3f}  slope {fits[k]['slope']:+.5f} "
              f"intercept {fits[k]['intercept']:+.6f}")

    # primary map: the best-fitting one, which is quadratic in the mean absolute log hazard
    # ratio. The intercept is the held-out cost of adding any column at all, the slope is what
    # signal buys. A quadratic form is what theory predicts, because the information a predictor
    # adds to a partial likelihood scales with the square of its coefficient.
    PRI = best
    a, c = fits[PRI]["slope"], fits[PRI]["intercept"]
    xcol, form = PRI.split("|")
    R["map_primary"] = {"form": f"dC = slope * ({xcol})^{{{form}}} + intercept",
                        "xcol": xcol, "shape": form, **fits[PRI]}

    def invert(dc, key):
        """Effect size that a measure would need for a given concordance gain."""
        s, i = fits[key]["slope"], fits[key]["intercept"]
        xc, fm = key.split("|")
        t = (dc - i) / s
        if fm == "quadratic":
            return float(np.sqrt(t)) if t > 0 else float("nan")
        if fm == "sqrt":
            return float(t ** 2) if t > 0 else float("nan")
        return float(t)

    # ---- Step 3, invert -----------------------------------------------------------------------
    tst_e = mp[mp.feature == "TST_min"].iloc[0]
    t90_e = mp[mp.feature == "spo2_pct_below_90"].iloc[0]
    E_req = invert(target, PRI)
    E_now = float(tst_e[xcol])
    E_t90 = float(t90_e[xcol])
    alt = {k: {"E_required": invert(target, k),
               "hr_required": float(np.exp(invert(target, k))),
               "fold_change": invert(target, k) / float(tst_e[k.split("|")[0]]),
               "r2": fits[k]["r2"]}
           for k in fits if np.isfinite(invert(target, k))}
    R["step3"] = {
        "effect_summary_used": xcol, "map_shape": form,
        "TST_E_now": E_now, "TST_hr_per_sd_now": float(np.exp(E_now)),
        "TST_E_rms": float(tst_e.E_rms), "TST_E_rms_adj": float(tst_e.E_rms_adj),
        "TST_max_hr_any_outcome": float(tst_e.max_hr),
        "TST_max_hr_short_sleep_threshold_published": 1.18,
        "T90_E": E_t90, "T90_hr_per_sd": float(np.exp(E_t90)),
        "E_required": float(E_req), "hr_required_per_sd": float(np.exp(E_req)),
        "fold_change_log_hazard": float(E_req / E_now) if E_now > 0 else np.inf,
        "target_dC": target, "across_all_map_forms": alt,
        "fold_change_range_across_maps": [float(np.nanmin([v["fold_change"] for v in alt.values()])),
                                          float(np.nanmax([v["fold_change"] for v in alt.values()]))],
    }
    print(f"\nSTEP 3  TST now: {xcol} {E_now:.4f} (HR {np.exp(E_now):.3f} per 1 SD), "
          f"largest single-outcome HR {tst_e.max_hr:.3f}")
    print(f"        required for top 10: {E_req:.4f} (HR {np.exp(E_req):.3f} per 1 SD), "
          f"a {E_req/E_now:.2f}-fold rise in the log hazard")
    print(f"        for scale, T90 sits at {E_t90:.4f} (HR {np.exp(E_t90):.3f} per 1 SD)")
    print(f"        across all six map forms the required fold change is "
          f"{R['step3']['fold_change_range_across_maps'][0]:.2f} to "
          f"{R['step3']['fold_change_range_across_maps'][1]:.2f}")

    # ---- Step 4, reliability and nights -------------------------------------------------------
    # observed effect = true effect * attenuation(reliability). Primary attenuation is
    # sqrt(reliability). The classical regression form, attenuation = reliability, is carried
    # as a sensitivity because it is more forgiving and still has to be beaten.
    def dc_of(e):
        return float(a * (e ** 2 if form == "quadratic" else
                          np.sqrt(e) if form == "sqrt" else e) + c)

    step4 = {}
    for tag, pw, icc in (("sqrt", 0.5, ICC_TST), ("linear", 1.0, ICC_TST),
                         ("sqrt_all_pairs", 0.5, ICC_TST_ALL)):
        req_r = icc * (E_req / E_now) ** (1.0 / pw) if E_now > 0 else np.inf
        ceiling_E = E_now * (1.0 / icc) ** pw
        step4[tag] = {
            "icc_used": icc,
            "attenuation": f"observed = true * reliability**{pw}",
            "reliability_required": float(req_r),
            "achievable": bool(req_r <= 1.0),
            "times_over_the_maximum_possible": float(req_r),
            "nights_required": float(spearman_brown_k(req_r, icc)) if req_r <= 1.0 else None,
            "ceiling_E_at_perfect_reliability": float(ceiling_E),
            "ceiling_hr_at_perfect_reliability": float(np.exp(ceiling_E)),
            "ceiling_dC_via_map": dc_of(ceiling_E),
            "ceiling_rank_in_reproduced_ranking": int((rk.dC > dc_of(ceiling_E)).sum() + 1),
            "max_multiple_of_log_hazard": float((1.0 / icc) ** pw),
        }
    R["step4"] = step4
    sb = []
    for k in [1, 2, 3, 5, 7, 10, 14, 30, 100, 1000]:
        r_ = reliability_k(k, ICC_TST)
        e_ = E_now * np.sqrt(r_ / ICC_TST)
        sb.append({"nights": k, "reliability": float(r_), "E_attainable": float(e_),
                   "hr_attainable": float(np.exp(e_)), "dC_via_map": dc_of(e_),
                   "rank_via_map": int((rk.dC > dc_of(e_)).sum() + 1)})
    e_ = E_now / np.sqrt(ICC_TST)
    sb.append({"nights": None, "reliability": 1.0, "E_attainable": float(e_),
               "hr_attainable": float(np.exp(e_)), "dC_via_map": dc_of(e_),
               "rank_via_map": int((rk.dC > dc_of(e_)).sum() + 1)})
    R["step4_nights_table"] = sb
    for tag in step4:
        print(f"\nSTEP 4  reliability required ({tag}): "
              f"{step4[tag]['reliability_required']:.3f}"
              f"   {'ACHIEVABLE' if step4[tag]['achievable'] else 'IMPOSSIBLE, exceeds 1.0'}")
    print(f"        ceiling at perfect measurement: HR "
          f"{step4['sqrt']['ceiling_hr_at_perfect_reliability']:.3f} per 1 SD, "
          f"dC {step4['sqrt']['ceiling_dC_via_map']:+.6f} against a target of {target:+.6f}, "
          f"rank {step4['sqrt']['ceiling_rank_in_reproduced_ranking']} of {N_MEASURES}")

    # ---- Step 5, the other reading -------------------------------------------------------------
    E_true_now = E_now / np.sqrt(ICC_TST)
    R["step5"] = {
        "TST_true_E_implied": float(E_true_now),
        "TST_true_hr_implied": float(np.exp(E_true_now)),
        "E_true_required": float(E_req),
        "hr_true_required": float(np.exp(E_req)),
        "multiple_of_current_true_log_hazard": float(E_req / E_true_now),
        "multiple_of_current_observed_log_hazard": float(E_req / E_now),
        "for_reference_T90_true_hr": float(np.exp(E_t90 / np.sqrt(ICC_T90))),
    }
    print(f"\nSTEP 5  true association implied now: HR {np.exp(E_true_now):.3f} per 1 SD. "
          f"Required: HR {np.exp(E_req):.3f}, {E_req/E_true_now:.2f} times the current "
          f"true log hazard.")

    # ---- Step 6 and 7, the simulation ----------------------------------------------------------
    iccs = repeat_iccs()
    R["repeat_iccs"] = iccs
    print(f"\nreliability of one night: TST raw {iccs['TST']['icc_raw']:.4f} "
          f"rank-normal {iccs['TST']['icc_rint']:.4f} (n={iccs['TST']['n_pairs']}), "
          f"T90 raw {iccs['T90']['icc_raw']:.4f} rank-normal {iccs['T90']['icc_rint']:.4f} "
          f"(n={iccs['T90']['n_pairs']})")
    sim = {}
    for scale in ("raw", "rint"):
        print(f"\nSTAGE E  simulation on the {scale} scale, "
              f"refitting the whole {N_RANKED_OUTCOMES}-outcome panel at each noise level")
        sx = stage_e(outs, scale, iccs)
        sim[scale] = {}
        for feat, nm in (("TST_min", "TST"), ("spo2_pct_below_90", "T90")):
            icc = iccs[nm]["icc_raw" if scale == "raw" else "icc_rint"]
            s = sx[sx.feature == feat]
            agg = (s.groupby("lam").dC.agg(["mean", "std", "count"]).reset_index()
                   .rename(columns={"mean": "dC", "std": "sd", "count": "reps"}))
            agg["reliability"] = [icc / (1 + l * (1 - icc)) for l in agg.lam]
            ex = extrapolate(agg, icc)
            best_ex = ex["linear_in_reliability"]
            sweep = []
            for k in [1, 2, 3, 5, 10, 30]:
                r = reliability_k(k, icc)
                sweep.append({"nights": k, "reliability": float(r),
                              "lam": float(lam_for_reliability(r, icc)),
                              "dC_linear_in_reliability": float(np.polyval(ex["_coef_rel"], r)),
                              "dC_quadratic_in_lambda": float(
                                  np.polyval(ex["_coef_quad"], lam_for_reliability(r, icc)))})
            sweep.append({"nights": None, "reliability": 1.0, "lam": -1.0,
                          "dC_linear_in_reliability": best_ex,
                          "dC_quadratic_in_lambda": ex["quadratic_in_lambda"]})
            preds = [ex["quadratic_in_lambda"], ex["linear_in_lambda"],
                     ex["linear_in_reliability"]]
            sim[scale][nm] = {
                "icc": float(icc), "n_pairs": iccs[nm]["n_pairs"],
                "observed_dC": float(agg.loc[agg.lam == 0, "dC"].iloc[0]),
                "levels": agg.to_dict("records"),
                "extrapolated_dC_at_perfect_reliability": ex,
                "primary_extrapolant": "linear_in_reliability",
                "dC_at_perfect_reliability": float(best_ex),
                "dC_at_perfect_range_over_extrapolants": [float(np.min(preds)),
                                                          float(np.max(preds))],
                "rank_at_perfect_reliability": int((rk.dC > best_ex).sum() + 1),
                "reaches_top10_at_perfect": bool(max(preds) > target),
                "reliability_sweep": sweep,
            }
            v = sim[scale][nm]
            print(f"  {nm}: observed dC {v['observed_dC']:+.6f}  ->  at perfect measurement "
                  f"{best_ex:+.6f} (linear in reliability, fit R2 "
                  f"{ex['fit_r2_linear_in_reliability']:.3f})")
            print(f"       other extrapolants {ex['quadratic_in_lambda']:+.6f} (quadratic), "
                  f"{ex['linear_in_lambda']:+.6f} (linear), "
                  f"{ex['rational_in_lambda']:+.6f} (rational, unstable)")
            print(f"       dC lost per unit of added noise {ex['slope_per_lambda']:+.6f}   "
                  f"rank at perfect measurement {v['rank_at_perfect_reliability']} of {N_MEASURES}   "
                  f"{'REACHES top 10' if v['reaches_top10_at_perfect'] else 'still below top 10'}")
    R["step6_7_simulation"] = sim
    reaches = any(sim[s]["TST"]["reaches_top10_at_perfect"] for s in sim)
    R["verdict"] = ("no achievable measurement improvement suffices" if not reaches
                    else "achievable")

    # ---- results.csv ---------------------------------------------------------------------------
    rows = []
    rows.append(["step1", "frozen dC rank 1", fz.loc[1, "feature"], fz.loc[1, "dC"], ""])
    rows.append(["step1", "frozen dC rank 5", fz.loc[5, "feature"], fz.loc[5, "dC"], ""])
    rows.append(["step1", "frozen dC rank 10 (the target)", fz.loc[10, "feature"], target_fz, ""])
    rows.append(["step1", "frozen dC rank 20", fz.loc[20, "feature"], fz.loc[20, "dC"], ""])
    rows.append(["step1", "frozen dC TST_min", f"rank {int(tst_fz['rank'])}", tst_fz.dC, ""])
    rows.append(["step1", "gain TST must add (frozen)", "", target_fz - tst_fz.dC, ""])
    rows.append(["stepA", "reproduced dC rank 10 (the target)", rkx.loc[10, "feature"], target, ""])
    rows.append(["stepA", "reproduced dC TST_min", f"rank {int(tst_r['rank'])}", tst_r.dC, ""])
    rows.append(["stepA", "reproduced dC T90", f"rank {int(t90_r['rank'])}", t90_r.dC, ""])
    rows.append(["stepA", "gain TST must add (reproduced)", "", target - tst_r.dC, ""])
    rows.append(["stepA", "rank agreement with frozen, Spearman", "",
                 R["reproduced"]["spearman_vs_frozen"], ""])
    for k_ in sorted(fits, key=lambda z: -fits[z]["r2"]):
        rows.append(["step2", "map R2", k_, fits[k_]["r2"], ""])
    rows.append(["step2", "map slope (primary)", PRI, a, ""])
    rows.append(["step2", "map intercept, held-out cost of any column", PRI, c, ""])
    rows.append(["step3", "TST summary log HR now", xcol, E_now, np.exp(E_now)])
    rows.append(["step3", "TST largest single-outcome HR", "", np.log(tst_e.max_hr), tst_e.max_hr])
    rows.append(["step3", "TST published max HR, short-sleep threshold", "", np.log(1.18), 1.18])
    rows.append(["step3", "T90 summary log HR", xcol, E_t90, np.exp(E_t90)])
    rows.append(["step3", "log HR required for top 10", PRI, E_req, np.exp(E_req)])
    rows.append(["step3", "fold change in log hazard required", PRI, E_req / E_now, ""])
    for k_, v_ in alt.items():
        rows.append(["step3_maps", "fold change required", k_, v_["fold_change"],
                     v_["hr_required"]])
    for tag in step4:
        rows.append([f"step4_{tag}", "reliability required", "",
                     step4[tag]["reliability_required"], ""])
        rows.append([f"step4_{tag}", "nights required", "",
                     step4[tag]["nights_required"] if step4[tag]["nights_required"] else
                     float("inf"), ""])
        rows.append([f"step4_{tag}", "dC ceiling at perfect measurement", "",
                     step4[tag]["ceiling_dC_via_map"], ""])
        rows.append([f"step4_{tag}", "rank at perfect measurement", "",
                     step4[tag]["ceiling_rank_in_reproduced_ranking"], ""])
    for r_ in sb:
        rows.append(["step4_nights", f"nights {r_['nights']}", f"R={r_['reliability']:.3f}",
                     r_["dC_via_map"], r_["hr_attainable"]])
    rows.append(["step5", "true HR implied now", "", E_true_now, np.exp(E_true_now)])
    rows.append(["step5", "true HR required", "", E_req, np.exp(E_req)])
    rows.append(["step5", "multiple of current true log hazard", "", E_req / E_true_now, ""])
    for scale in sim:
        for nm in ("TST", "T90"):
            v = sim[scale][nm]
            for lv in v["levels"]:
                rows.append([f"step6_{nm}_{scale}", f"lambda {lv['lam']}",
                             f"reliability {lv['reliability']:.3f}", lv["dC"], lv["sd"]])
            for k_, v_ in v["extrapolated_dC_at_perfect_reliability"].items():
                if not k_.startswith("_"):
                    rows.append([f"step6_{nm}_{scale}", k_, "", v_, ""])
            rows.append([f"step6_{nm}_{scale}", "rank at perfect reliability", "",
                         v["rank_at_perfect_reliability"], ""])
            for s_ in v["reliability_sweep"]:
                rows.append([f"step6_{nm}_{scale}_sweep", f"nights {s_['nights']}",
                             f"R={s_['reliability']:.3f}", s_["dC_linear_in_reliability"], ""])
    out = pd.DataFrame(rows, columns=["step", "quantity", "detail", "value", "extra"])
    out.to_csv(f"{HERE}/results.csv", index=False)
    json.dump(R, open(f"{HERE}/summary.json", "w"), indent=1, default=str)
    print(f"\nVERDICT: {R['verdict']}")
    print(f"wrote results.csv and summary.json  ({time.time()-t0:.0f}s total)")


if __name__ == "__main__":
    main()
