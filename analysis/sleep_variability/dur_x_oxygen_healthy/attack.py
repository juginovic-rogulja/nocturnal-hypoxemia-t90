"""
Hostile review of the healthy-subgroup 2x2 discordance result.

Everything is recomputed from data_frozen_v7_2026-09/t90_final.parquet. Nothing is read back from
model.py's outputs, so a shared bug in that script cannot survive into this one. The healthy
definition, the cells and the model specification are rebuilt here from the same written rule.

Six attacks:
  A1  power, confidence-interval widths and minimum detectable hazard ratios
  A2  the whole-recording denominator, using the per-stage extraction
  A3  is "healthy" actually healthy
  A4  selection on the outcome, the low-oxygen healthy cells as a survivor stratum
  A5  negative controls per contrast, and what they strike
  A6  threshold sensitivity, TST at 4/5/6 h crossed with T90 at 5%/10%

Modelling rules obeyed throughout: site-stratified Cox, unpenalized (penalizer=0), sex covariate,
PSG date as time zero, and a 4-df natural cubic age spline with ONE COLUMN DROPPED. The script
asserts that no fit returned all-zero age coefficients.

Writes attack.json plus attack_*.csv in this folder.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import glob
import json
import sys
import time
import warnings

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from patsy import dmatrix
from scipy import stats

warnings.filterwarnings("ignore")

T90ROOT = paths.T90_ROOT
HERE = f"{paths.SV_ROOT}/dur_x_oxygen_healthy"
STAGEDIR = f"{paths.V8_ROOT}/X4_groupP/t90_by_stage"   # v8.2: re-extracted, collapse-masked shards (X4_groupP)
sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
from cohort_spec import apply_cohort, COHORT_N            # noqa: E402
from disease_definitions import (DISEASES, NEGATIVE_CONTROLS,   # noqa: E402
                                 CIRCULAR, ORGAN_GROUP)

T0 = time.time()
LOG = []


def log(m):
    print(m, flush=True)
    LOG.append(str(m))


def sec(t):
    log("\n" + "=" * 88)
    log(t)
    log("=" * 88)


# thresholds, primary specification
T90_LOW, T90_NORMAL = 10.0, 1.0
TST_NORMAL = 360.0
TST_PRIMARY = 300.0            # model.py's primary, chosen because 4 h leaves 46 in cell 2
TST_MATCHED = 240.0            # the full-cohort run's cut
CELLS = ["cell1", "cell2", "cell3", "cell4"]
CELLD = ["cell1", "cell2", "cell3"]
MIN_EVENTS = 20
MIN_CELL_EVENTS_PAIR = 3
ZA, ZB = 1.959964, 0.8416212   # two-sided alpha 0.05, 80% power

OUT = {"generated": time.strftime("%Y-%m-%d %H:%M"), "attacks": {}}

# ======================================================================================
# 0. COHORT, HEALTHY DEFINITION AND CELLS, REBUILT INDEPENDENTLY
# ======================================================================================
sec("0. INDEPENDENT REBUILD OF THE HEALTHY SUBGROUP AND THE CELLS")
b = apply_cohort(pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet"))
log(f"analysis cohort {len(b):,} (cohort_spec COHORT_N {COHORT_N:,})")

b["male"] = (b.sex.astype(str).str.upper().str[0] == "M").astype(float)
b.loc[~b.sex.astype(str).str.upper().str[0].isin(["M", "F"]), "male"] = np.nan

ORGAN_SYSTEMS = ("Cardiac", "Respiratory", "Kidney", "Liver", "Metabolic")
EXCLUDE_KEYS = [k for k, g in ORGAN_GROUP.items() if g in ORGAN_SYSTEMS]
EXCLUDE_KEYS += ["stroke_any", "cvd", "cancer_any", "prostate_ca"]
EXCLUDE_KEYS = list(dict.fromkeys(EXCLUDE_KEYS))
NOT_EXCLUDED = [k for k in DISEASES if k not in EXCLUDE_KEYS]

keep = np.ones(len(b), bool)
for k in EXCLUDE_KEYS:
    keep &= (b[f"{k}_prevalent"] == 0).values
h = b[keep].copy()
N_HEALTHY = len(h)
log(f"healthy (organ-disease-free) {N_HEALTHY:,} = {100 * N_HEALTHY / len(b):.1f}% of cohort")
log(f"  excluded on {len(EXCLUDE_KEYS)} conditions; sleep duration is NOT one of them; "
    f"insomnia, restless legs, nocturia and sleep apnoea all left in")
OUT["rebuild"] = {"cohort_n": int(len(b)), "healthy_n": int(N_HEALTHY),
                  "pct_of_cohort": round(100 * N_HEALTHY / len(b), 1),
                  "n_exclusion_conditions": len(EXCLUDE_KEYS),
                  "matches_model_py_4822": bool(N_HEALTHY == 4822)}


def assign_cells(d, tst_short, t90_low=T90_LOW, t90_norm=T90_NORMAL, tst_norm=TST_NORMAL):
    lo, no = d.spo2_pct_below_90 > t90_low, d.spo2_pct_below_90 <= t90_norm
    sh, ns = d.TST_min < tst_short, d.TST_min >= tst_norm
    return pd.Series(np.select([sh & no, sh & lo, ns & lo, ns & no], CELLS, default="excluded"),
                     index=d.index)


h["cell"] = assign_cells(h, TST_PRIMARY)
b["cell"] = assign_cells(b, TST_PRIMARY)
h240 = h.assign(cell=assign_cells(h, TST_MATCHED))
sizes = {k: int((h.cell == k).sum()) for k in CELLS}
sizes240 = {k: int((h240.cell == k).sum()) for k in CELLS}
log(f"primary cells (TST<300, T90>10 vs <=1): {sizes}")
log(f"matched  cells (TST<240):               {sizes240}")
OUT["rebuild"]["cells_tst300"] = sizes
OUT["rebuild"]["cells_tst240"] = sizes240
OUT["rebuild"]["reproduces_model_py_cells"] = bool(
    sizes == {"cell1": 681, "cell2": 109, "cell3": 181, "cell4": 1600})
log(f"reproduces model.py cell sizes: {OUT['rebuild']['reproduces_model_py_cells']}")


# ======================================================================================
# MODEL MACHINERY
# ======================================================================================
SPLINE_MAXABS = []


def agespline(d):
    sp = dmatrix("cr(a, df=4) - 1", {"a": d.AgeAtVisit.values}, return_type="dataframe")
    sp.columns = [f"age_s{i}" for i in range(4)]
    sp.index = d.index
    return sp.iloc[:, 1:]          # ONE COLUMN DROPPED, mandatory


OUTCOMES = [(k, v[0]) for k, v in DISEASES.items() if k not in CIRCULAR]
OUTCOMES.append(("death", "Death from any cause"))


def frame_for(d, key, extra=()):
    yc, ec, pc = f"{key}_years", f"{key}_incident", f"{key}_prevalent"
    if not all(c in d.columns for c in (yc, ec, pc)):
        return None
    f = d[(d[pc] == 0) & d[yc].notna() & (d[yc] > 0)].copy()
    if not len(f):
        return None
    sp = agespline(f)
    cols = {"T": f[yc], "E": f[ec].astype(int), "site": f.site_id, "male": f.male,
            "cell_lab": f.cell}
    for c in extra:
        cols[c] = f[c]
    for c in sp.columns:
        cols[c] = sp[c]
    return pd.DataFrame(cols).dropna()


def _fit(dd, terms):
    cov = terms + ["male"] + [c for c in dd.columns if c.startswith("age_s")]
    m = CoxPHFitter(penalizer=0.0).fit(dd[cov + ["T", "E", "site"]], "T", "E", strata=["site"])
    SPLINE_MAXABS.append(float(np.abs(
        m.params_[[c for c in m.params_.index if c.startswith("age_s")]]).max()))
    return m


def pairfit(dd, a, ref, min_events=MIN_EVENTS, min_cell=MIN_CELL_EVENTS_PAIR):
    sub = dd[dd.cell_lab.isin([a, ref])].copy()
    ea, er = int(sub[sub.cell_lab == a].E.sum()), int(sub[sub.cell_lab == ref].E.sum())
    na, nr = int((sub.cell_lab == a).sum()), int((sub.cell_lab == ref).sum())
    if ea + er < min_events or min(ea, er) < min_cell:
        return None
    sub["X"] = (sub.cell_lab == a).astype(float)
    try:
        m = _fit(sub, ["X"])
    except Exception:
        return None
    s = m.summary.loc["X"]
    hr, lo, hi = (float(s["exp(coef)"]), float(s["exp(coef) lower 95%"]),
                  float(s["exp(coef) upper 95%"]))
    return {"hr": round(hr, 3), "lo": round(lo, 3), "hi": round(hi, 3), "p": float(s["p"]),
            "se": float(s["se(coef)"]), "ev_a": ea, "ev_ref": er, "n_a": na, "n_ref": nr,
            "ci_ratio": round(hi / lo, 2) if lo > 0 else np.nan,
            "log_ci_width": round(float(np.log(hi) - np.log(lo)), 3)}


PAIRS = [("cell1", "cell4"), ("cell2", "cell4"), ("cell3", "cell4"), ("cell2", "cell1")]


def run_panel(d, tag, pairs=PAIRS, min_events=MIN_EVENTS):
    rows = []
    for key, label in OUTCOMES:
        dd = frame_for(d, key)
        if dd is None or not len(dd):
            continue
        ev = {k2: int(dd[dd.cell_lab == k2].E.sum()) for k2 in CELLS}
        nn = {k2: int((dd.cell_lab == k2).sum()) for k2 in CELLS}
        rec = {"analysis": tag, "key": key, "disease": label,
               "organ_group": ORGAN_GROUP.get(key, "Death"),
               "negative_control": key in NEGATIVE_CONTROLS,
               "used_as_healthy_exclusion": key in EXCLUDE_KEYS,
               "events_total": int(dd.E.sum()),
               **{f"n_{k2}": nn[k2] for k2 in CELLS}, **{f"ev_{k2}": ev[k2] for k2 in CELLS}}
        for a, ref in pairs:
            r = pairfit(dd, a, ref, min_events=min_events)
            pre = f"{a}_vs_{ref}"
            if r is None:
                rec[f"{pre}_hr"] = np.nan
                continue
            for kk, vv in r.items():
                rec[f"{pre}_{kk}"] = vv
        rows.append(rec)
    return pd.DataFrame(rows)


def fdr(p):
    p = np.asarray(p, float)
    ok = ~np.isnan(p)
    q = np.full(len(p), np.nan)
    pv = p[ok]
    n = len(pv)
    if n == 0:
        return q
    order = np.argsort(pv)
    ranked = pv[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(ranked, 1)
    q[ok] = out
    return q


def tally(df, pre, exclude_neg=True):
    s = df[~df.negative_control] if exclude_neg else df
    s = s[s[f"{pre}_hr"].notna()]
    if not len(s):
        return {"n": 0}
    p = s[f"{pre}_p"].values
    q = fdr(p)
    sig = ((s[f"{pre}_lo"] > 1) | (s[f"{pre}_hi"] < 1)).values
    return {"n": int(len(s)), "n_sig": int(sig.sum()), "expected_by_chance": round(.05 * len(s), 1),
            "n_fdr": int((q < .05).sum()), "median_hr": round(float(s[f"{pre}_hr"].median()), 3),
            "min_hr": round(float(s[f"{pre}_hr"].min()), 3),
            "max_hr": round(float(s[f"{pre}_hr"].max()), 3)}


def controls(df, pre):
    s = df[df.negative_control & df[f"{pre}_hr"].notna()]
    if not len(s):
        return {"n_controls": 0, "floor": None}
    hrs = s[f"{pre}_hr"].astype(float)
    floor = float(max(hrs.max(), 1 / hrs.min()))
    return {"n_controls": int(len(s)), "floor": round(floor, 3),
            "range": [round(float(hrs.min()), 3), round(float(hrs.max()), 3)],
            "n_significant_controls": int(((s[f"{pre}_lo"] > 1) | (s[f"{pre}_hi"] < 1)).sum()),
            "controls": {r.disease: [float(r[f"{pre}_hr"]), float(r[f"{pre}_lo"]),
                                     float(r[f"{pre}_hi"])] for _, r in s.iterrows()}}


# ======================================================================================
# A1. POWER
# ======================================================================================
sec("A1. POWER. THE MAIN THREAT.")
A1 = OUT["attacks"]["A1_power"] = {}

prim = run_panel(h, "healthy_tst300")
prim.to_csv(f"{HERE}/attack_primary_panel.csv", index=False)
m240 = run_panel(h240, "healthy_tst240")
m240.to_csv(f"{HERE}/attack_matched240_panel.csv", index=False)

A1["cell_sizes"] = sizes
A1["cell_sizes_tst240"] = sizes240
A1["smallest_cell"] = {"cell": "cell2", "n_tst300": sizes["cell2"], "n_tst240": sizes240["cell2"]}
A1["null_in_cell2_is_uninformative_at_tst240"] = bool(sizes240["cell2"] < 50)
log(f"cell 2 holds {sizes['cell2']} people at the 5-hour cut and {sizes240['cell2']} at the "
    f"4-hour cut. The 4-hour cell 2 is below 50, so any null there is uninformative.")

# --- every CI width in the smallest cell
wid = []
for _, r in prim.iterrows():
    for pre in ("cell2_vs_cell1", "cell2_vs_cell4"):
        if pd.isna(r.get(f"{pre}_hr", np.nan)):
            continue
        wid.append({"contrast": pre, "disease": r.disease,
                    "negative_control": bool(r.negative_control),
                    "hr": r[f"{pre}_hr"], "lo": r[f"{pre}_lo"], "hi": r[f"{pre}_hi"],
                    "ci_ratio_hi_over_lo": r[f"{pre}_ci_ratio"],
                    "log_ci_width": r[f"{pre}_log_ci_width"],
                    "ev_cell2": int(r.ev_cell2), "ev_ref": int(r[f"{pre}_ev_ref"]),
                    "significant": bool(r[f"{pre}_lo"] > 1 or r[f"{pre}_hi"] < 1)})
wid = pd.DataFrame(wid).sort_values("ci_ratio_hi_over_lo")
wid.to_csv(f"{HERE}/attack_ci_widths.csv", index=False)
for pre in ("cell2_vs_cell1", "cell2_vs_cell4"):
    s = wid[wid.contrast == pre]
    A1[f"ci_width_{pre}"] = {
        "n_estimable": int(len(s)),
        "median_ci_ratio": round(float(s.ci_ratio_hi_over_lo.median()), 2),
        "min_ci_ratio": round(float(s.ci_ratio_hi_over_lo.min()), 2),
        "max_ci_ratio": round(float(s.ci_ratio_hi_over_lo.max()), 2),
        "n_ci_ratio_over_5": int((s.ci_ratio_hi_over_lo > 5).sum()),
        "n_ci_ratio_over_10": int((s.ci_ratio_hi_over_lo > 10).sum()),
        "widest": s.nlargest(3, "ci_ratio_hi_over_lo")[
            ["disease", "hr", "lo", "hi", "ci_ratio_hi_over_lo"]].to_dict("records"),
        "narrowest": s.nsmallest(3, "ci_ratio_hi_over_lo")[
            ["disease", "hr", "lo", "hi", "ci_ratio_hi_over_lo"]].to_dict("records")}
    e = A1[f"ci_width_{pre}"]
    log(f"{pre}: {e['n_estimable']} estimable, median CI spans a factor of "
        f"{e['median_ci_ratio']}x (range {e['min_ci_ratio']}x to {e['max_ci_ratio']}x); "
        f"{e['n_ci_ratio_over_5']} span more than 5x, {e['n_ci_ratio_over_10']} more than 10x")

# --- minimum detectable HR at 80% power, Schoenfeld
def mdhr(events, n_a, n_ref):
    if events <= 0 or n_a <= 0 or n_ref <= 0:
        return np.nan
    p = n_a / (n_a + n_ref)
    v = events * p * (1 - p)
    if v <= 0:
        return np.nan
    return float(np.exp((ZA + ZB) / np.sqrt(v)))


mdrows = []
for _, r in prim.iterrows():
    for a, ref in PAIRS:
        pre = f"{a}_vs_{ref}"
        na, nr = int(r[f"n_{a}"]), int(r[f"n_{ref}"])
        ea, er = int(r[f"ev_{a}"]), int(r[f"ev_{ref}"])
        mdrows.append({"contrast": pre, "key": r.key, "disease": r.disease,
                       "negative_control": bool(r.negative_control),
                       "n_a": na, "n_ref": nr, "ev_a": ea, "ev_ref": er, "events": ea + er,
                       "mdhr80": round(mdhr(ea + er, na, nr), 3) if ea + er else np.nan,
                       "estimable": bool(not pd.isna(r.get(f"{pre}_hr", np.nan))),
                       "observed_hr": r.get(f"{pre}_hr", np.nan)})
md = pd.DataFrame(mdrows)
md.to_csv(f"{HERE}/attack_min_detectable_hr.csv", index=False)
HEAD = ["death", "hf", "ckd", "diabetes", "obesity", "htn2", "cvd", "ihd", "resp_failure",
        "copd2", "afib", "dementia", "stroke_any"]
for pre in ["cell2_vs_cell1", "cell1_vs_cell4", "cell2_vs_cell4", "cell3_vs_cell4"]:
    s = md[md.contrast == pre]
    A1[f"mdhr80_{pre}"] = {
        "median_over_all_outcomes": round(float(s.mdhr80.median()), 3),
        "headline": {r.disease: {"events": int(r.events), "mdhr80": r.mdhr80,
                                 "estimable": bool(r.estimable),
                                 "observed_hr": (None if pd.isna(r.observed_hr)
                                                 else float(r.observed_hr))}
                     for _, r in s[s.key.isin(HEAD)].iterrows()}}
    log(f"{pre}: median minimum detectable HR at 80% power "
        f"{A1[f'mdhr80_{pre}']['median_over_all_outcomes']}")
for k in ["death", "hf", "ckd", "diabetes", "htn2"]:
    r = md[(md.contrast == "cell2_vs_cell1") & (md.key == k)]
    if len(r):
        r = r.iloc[0]
        log(f"  cell2 vs cell1 {r.disease:26s} {int(r.events):4d} events, detectable only above "
            f"HR {r.mdhr80}")

# --- what the data cannot support: nulls that are not evidence of absence
uninf = []
for _, r in prim.iterrows():
    for pre in ("cell2_vs_cell1", "cell2_vs_cell4"):
        hr = r.get(f"{pre}_hr", np.nan)
        if pd.isna(hr):
            continue
        lo, hi = r[f"{pre}_lo"], r[f"{pre}_hi"]
        if lo <= 1 <= hi and (hi > 2.0 or lo < 0.5):
            uninf.append({"contrast": pre, "disease": r.disease, "hr": hr, "lo": lo, "hi": hi})
A1["nulls_that_are_not_evidence_of_absence"] = uninf
A1["n_uninformative_nulls"] = len(uninf)
log(f"{len(uninf)} of the non-significant cell-2 results have a confidence interval that still "
    f"admits a doubling of risk or a halving. Those are not evidence of absence.")

# --- power of the interaction test, the claim most at risk
ip = []
for key, label in OUTCOMES:
    dd = frame_for(h, key)
    if dd is None:
        continue
    dd = dd[dd.cell_lab != "excluded"]
    ev = {k2: int(dd[dd.cell_lab == k2].E.sum()) for k2 in CELLS}
    if min(ev.values()) < 3 or sum(ev.values()) < MIN_EVENTS:
        continue
    dd = dd.copy()
    dd["short"] = dd.cell_lab.isin(["cell1", "cell2"]).astype(float)
    dd["lowo2"] = dd.cell_lab.isin(["cell2", "cell3"]).astype(float)
    dd["ix"] = dd.short * dd.lowo2
    try:
        m = _fit(dd, ["short", "lowo2", "ix"])
    except Exception:
        continue
    s = m.summary.loc["ix"]
    ip.append({"disease": label, "hr": round(float(s["exp(coef)"]), 3),
               "lo": round(float(s["exp(coef) lower 95%"]), 3),
               "hi": round(float(s["exp(coef) upper 95%"]), 3), "p": float(s["p"]),
               "se": float(s["se(coef)"]),
               "ci_ratio": round(float(s["exp(coef) upper 95%"] / s["exp(coef) lower 95%"]), 1),
               "ev_cell2": ev["cell2"], "events": int(sum(ev.values())),
               "mdhr80_interaction": round(float(np.exp((ZA + ZB) * float(s["se(coef)"]))), 2)})
ipd = pd.DataFrame(ip)
ipd.to_csv(f"{HERE}/attack_interaction_power.csv", index=False)
A1["interaction_power"] = {
    "n_tested": int(len(ipd)),
    "n_nominally_significant": int(((ipd.lo > 1) | (ipd.hi < 1)).sum()),
    "median_ci_ratio": round(float(ipd.ci_ratio.median()), 1),
    "median_mdhr80": round(float(ipd.mdhr80_interaction.median()), 2),
    "min_mdhr80": round(float(ipd.mdhr80_interaction.min()), 2),
    "n_with_mdhr80_over_2": int((ipd.mdhr80_interaction > 2).sum())}
e = A1["interaction_power"]
log(f"interaction: {e['n_tested']} testable, median CI spans {e['median_ci_ratio']}x, and the "
    f"smallest interaction HR detectable at 80% power is {e['min_mdhr80']} (median "
    f"{e['median_mdhr80']}). {e['n_with_mdhr80_over_2']} of {e['n_tested']} cannot detect even a "
    f"doubling of the multiplicative effect.")


# ======================================================================================
# A2. THE WHOLE-RECORDING DENOMINATOR
# ======================================================================================
sec("A2. THE WHOLE-RECORDING DENOMINATOR")
A2 = OUT["attacks"]["A2_denominator"] = {}

files = sorted(glob.glob(f"{STAGEDIR}/out_*.csv"))
stg = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
raw_n = len(stg)
stg = stg[(stg.status == "ok") & stg.t90_sleep.notna() & stg.t90_wake.notna() &
          (stg.n_sleep > 0) & (stg.n_wake > 0)]
stg = stg.sort_values("n_sleep", ascending=False).drop_duplicates("BDSPPatientID")
SCOL = ["BDSPPatientID", "t90_all", "t90_wake", "t90_sleep", "n_sleep", "n_wake",
        "mean_wake", "mean_sleep", "tst_min", "rec_min"]
mg = h.merge(stg[SCOL], on="BDSPPatientID", how="inner", suffixes=("", "_stg"))
A2["stage_files"] = len(files)
A2["stage_rows_raw"] = int(raw_n)
A2["stage_rows_usable"] = int(len(stg))
A2["merged_into_healthy"] = int(len(mg))
A2["coverage_pct_of_healthy"] = round(100 * len(mg) / N_HEALTHY, 1)
A2["covered_by_cell"] = {k: int((mg.cell == k).sum()) for k in CELLS}
A2["coverage_pct_by_cell"] = {k: round(100 * (mg.cell == k).sum() / max(sizes[k], 1), 1)
                              for k in CELLS}
A2["sanity_rho_stage_t90_all_vs_frozen_t90"] = round(float(
    stats.spearmanr(mg.t90_all, mg.spo2_pct_below_90)[0]), 4)
log(f"per-stage extraction: {raw_n:,} rows, {len(stg):,} usable, merging onto {len(mg):,} of the "
    f"{N_HEALTHY:,} healthy ({A2['coverage_pct_of_healthy']}%)")
log(f"coverage inside the cells {A2['covered_by_cell']} = {A2['coverage_pct_by_cell']}% of each")
log(f"sanity check, whole-recording T90 in the extraction against the frozen T90, Spearman rho "
    f"{A2['sanity_rho_stage_t90_all_vs_frozen_t90']}")

A2["measured_by_cell"] = {}
for k in CELLS:
    s = mg[mg.cell == k]
    if len(s) < 5:
        A2["measured_by_cell"][k] = {"n": int(len(s)), "too_small": True}
        continue
    A2["measured_by_cell"][k] = {
        "n": int(len(s)),
        "sleep_t90_median": round(float(s.t90_sleep.median()), 3),
        "sleep_t90_iqr": [round(float(s.t90_sleep.quantile(.25)), 3),
                          round(float(s.t90_sleep.quantile(.75)), 3)],
        "wake_t90_median": round(float(s.t90_wake.median()), 3),
        "wake_t90_iqr": [round(float(s.t90_wake.quantile(.25)), 3),
                         round(float(s.t90_wake.quantile(.75)), 3)],
        "pct_wake_t90_over_10": round(100 * float((s.t90_wake > 10).mean()), 1),
        "pct_wake_t90_over_30": round(100 * float((s.t90_wake > 30).mean()), 1),
        "pct_wake_t90_exceeds_sleep_t90": round(100 * float((s.t90_wake > s.t90_sleep).mean()), 1),
        "mean_wake_spo2_median": round(float(s.mean_wake.median()), 2),
        "wake_minutes_median": round(float((s.rec_min - s.tst_min).median()), 1)}
    e = A2["measured_by_cell"][k]
    log(f"  {k} (n={e['n']}): measured ASLEEP T90 {e['sleep_t90_median']}%, AWAKE T90 "
        f"{e['wake_t90_median']}%, {e['pct_wake_t90_over_10']}% have awake T90 above 10%, "
        f"mean awake saturation {e['mean_wake_spo2_median']}%")

c1m, c2m = mg[mg.cell == "cell1"], mg[mg.cell == "cell2"]
A2["misclassification"] = {
    "cell1_n": int(len(c1m)),
    "cell1_pct_sleep_t90_over_1": round(100 * float((c1m.t90_sleep > 1).mean()), 1),
    "cell1_pct_sleep_t90_over_5": round(100 * float((c1m.t90_sleep > 5).mean()), 1),
    "cell1_pct_sleep_t90_over_10": round(100 * float((c1m.t90_sleep > 10).mean()), 1),
    "cell2_n": int(len(c2m)),
    "cell2_pct_sleep_t90_at_or_below_10": round(100 * float((c2m.t90_sleep <= 10).mean()), 1),
    "cell2_pct_sleep_t90_at_or_below_1": round(100 * float((c2m.t90_sleep <= 1).mean()), 1)}
log(f"  misclassification: {A2['misclassification']['cell1_pct_sleep_t90_over_1']}% of the "
    f"healthy 'normal oxygen' short sleepers had asleep T90 above 1%, and "
    f"{A2['misclassification']['cell2_pct_sleep_t90_at_or_below_10']}% of the healthy 'low "
    f"oxygen' short sleepers had asleep T90 at or below 10%")

# rebuild the cells on MEASURED ASLEEP T90 and refit
mgs = mg.copy()
mgs["cell"] = assign_cells(mgs.assign(spo2_pct_below_90=mgs.t90_sleep), TST_PRIMARY)
A2["sleep_denominated_cells"] = {"sizes": {k: int((mgs.cell == k).sum()) for k in CELLS}}
r_sd = run_panel(mgs, "sleep_denominated", min_events=15)
r_sd.to_csv(f"{HERE}/attack_sleep_denominated.csv", index=False)
mgo = mg.copy()   # same people, original whole-recording cells, so the comparison is like for like
r_wr = run_panel(mgo, "same_subset_whole_recording", min_events=15)
r_wr.to_csv(f"{HERE}/attack_same_subset_whole_recording.csv", index=False)
for pre in ["cell1_vs_cell4", "cell2_vs_cell1", "cell2_vs_cell4", "cell3_vs_cell4"]:
    A2["sleep_denominated_cells"][pre] = tally(r_sd, pre)
    A2.setdefault("same_subset_whole_recording", {})[pre] = tally(r_wr, pre)
log(f"  cells rebuilt on MEASURED ASLEEP T90: {A2['sleep_denominated_cells']['sizes']}")
for pre in ["cell1_vs_cell4", "cell2_vs_cell1", "cell3_vs_cell4"]:
    a, c = A2["sleep_denominated_cells"][pre], A2["same_subset_whole_recording"][pre]
    log(f"    {pre}: sleep-denominated {a.get('n_sig')} of {a.get('n')} sig, median HR "
        f"{a.get('median_hr')}  |  whole-recording on the same people {c.get('n_sig')} of "
        f"{c.get('n')}, median HR {c.get('median_hr')}")


# ======================================================================================
# A3. IS "HEALTHY" ACTUALLY HEALTHY
# ======================================================================================
sec("A3. IS 'HEALTHY' ACTUALLY HEALTHY")
A3 = OUT["attacks"]["A3_healthy_label"] = {}

bmi = pd.read_parquet(f"{paths.NUMBERS_DIR}/_bmi_derived.parquet")[["BDSPPatientID", "bmi"]]
hb = h.merge(bmi, on="BDSPPatientID", how="left")
# total prevalent burden. every condition, and separately only the ones NOT used to define health
prevcols_all = [f"{k}_prevalent" for k in DISEASES if f"{k}_prevalent" in h.columns]
prevcols_rem = [f"{k}_prevalent" for k in NOT_EXCLUDED if f"{k}_prevalent" in h.columns]
hb["n_prev_all"] = hb[prevcols_all].sum(axis=1)
hb["n_prev_remaining"] = hb[prevcols_rem].sum(axis=1)
A3["bmi_available_pct_of_healthy"] = round(100 * float(hb.bmi.notna().mean()), 1)
A3["by_cell"] = {}
for k in CELLS:
    s = hb[hb.cell == k]
    A3["by_cell"][k] = {
        "n": int(len(s)),
        "age_median": round(float(s.AgeAtVisit.median()), 1),
        "age_iqr": [round(float(s.AgeAtVisit.quantile(.25)), 1),
                    round(float(s.AgeAtVisit.quantile(.75)), 1)],
        "pct_male": round(100 * float(s.male.mean()), 1),
        "bmi_n": int(s.bmi.notna().sum()),
        "bmi_median": (None if s.bmi.notna().sum() < 5 else round(float(s.bmi.median()), 1)),
        "bmi_iqr": (None if s.bmi.notna().sum() < 5 else
                    [round(float(s.bmi.quantile(.25)), 1), round(float(s.bmi.quantile(.75)), 1)]),
        "pct_bmi_ge30": (None if s.bmi.notna().sum() < 5 else
                         round(100 * float((s.bmi >= 30).mean()), 1)),
        "ahi_median": round(float(s.AHI.median()), 1),
        "pct_ahi_ge15": round(100 * float((s.AHI >= 15).mean()), 1),
        "pct_ahi_ge30": round(100 * float((s.AHI >= 30).mean()), 1),
        "n_prev_remaining_mean": round(float(s.n_prev_remaining.mean()), 2),
        "n_prev_remaining_median": float(s.n_prev_remaining.median()),
        "pct_with_any_remaining_prev": round(100 * float((s.n_prev_remaining > 0).mean()), 1),
        "pct_prev_dementia": round(100 * float((s.dementia_prevalent == 1).mean()), 2),
        "pct_prev_depression": round(100 * float((s.depression_prevalent == 1).mean()), 1),
        "pct_prev_anaemia": round(100 * float((s.anaemia_prevalent == 1).mean()), 1),
        "followup_years_median": round(float(s.death_years.median()), 2)}
    e = A3["by_cell"][k]
    log(f"  {k:6s} n={e['n']:5d} age {e['age_median']:>5} {e['age_iqr']}  male {e['pct_male']:>5}%"
        f"  BMI {e['bmi_median']} (n={e['bmi_n']})  AHI {e['ahi_median']:>5}"
        f"  AHI>=15 {e['pct_ahi_ge15']:>5}%  other prevalent conditions "
        f"{e['n_prev_remaining_mean']:.2f}")

A3["age_gap_cell2_minus_cell4"] = round(A3["by_cell"]["cell2"]["age_median"] -
                                        A3["by_cell"]["cell4"]["age_median"], 1)
if A3["by_cell"]["cell2"]["bmi_median"] and A3["by_cell"]["cell4"]["bmi_median"]:
    A3["bmi_gap_cell2_minus_cell4"] = round(A3["by_cell"]["cell2"]["bmi_median"] -
                                            A3["by_cell"]["cell4"]["bmi_median"], 1)
    A3["bmi_gap_cell2_minus_cell1"] = round(A3["by_cell"]["cell2"]["bmi_median"] -
                                            A3["by_cell"]["cell1"]["bmi_median"], 1)
# formal tests
c1, c2, c4 = [hb[hb.cell == k] for k in ("cell1", "cell2", "cell4")]
A3["tests"] = {
    "age_c2_vs_c1_mannwhitney_p": float(stats.mannwhitneyu(c2.AgeAtVisit, c1.AgeAtVisit)[1]),
    "age_c2_vs_c4_mannwhitney_p": float(stats.mannwhitneyu(c2.AgeAtVisit, c4.AgeAtVisit)[1]),
    "ahi_c2_vs_c1_mannwhitney_p": float(stats.mannwhitneyu(c2.AHI.dropna(), c1.AHI.dropna())[1]),
    "nprev_c2_vs_c1_mannwhitney_p": float(
        stats.mannwhitneyu(c2.n_prev_remaining, c1.n_prev_remaining)[1])}
if c2.bmi.notna().sum() >= 5 and c1.bmi.notna().sum() >= 5:
    A3["tests"]["bmi_c2_vs_c1_mannwhitney_p"] = float(
        stats.mannwhitneyu(c2.bmi.dropna(), c1.bmi.dropna())[1])
    A3["tests"]["bmi_c2_vs_c4_mannwhitney_p"] = float(
        stats.mannwhitneyu(c2.bmi.dropna(), c4.bmi.dropna())[1])
pd.DataFrame(A3["by_cell"]).T.to_csv(f"{HERE}/attack_cell_confounders.csv")

# --- can the confounders be adjusted away? BMI is recorded on a tenth of the healthy subgroup,
# so it cannot. Apnoea severity can, and is the strongest available test of whether cell 2 is
# just sicker breathing rather than lower oxygen.
h_ahi = h[h.AHI.notna()].copy()
A3["ahi_adjustment"] = {"n_with_ahi": int(len(h_ahi)),
                        "pct_of_healthy": round(100 * len(h_ahi) / N_HEALTHY, 1),
                        "bmi_cannot_be_adjusted_n_in_cell2": A3["by_cell"]["cell2"]["bmi_n"]}
rows_a = []
for key, label in OUTCOMES:
    dd = frame_for(h_ahi, key, extra=("AHI",))
    if dd is None:
        continue
    for a, ref in [("cell2", "cell1"), ("cell1", "cell4")]:
        sub = dd[dd.cell_lab.isin([a, ref])].copy()
        ea = int(sub[sub.cell_lab == a].E.sum())
        er = int(sub[sub.cell_lab == ref].E.sum())
        if ea + er < MIN_EVENTS or min(ea, er) < MIN_CELL_EVENTS_PAIR:
            continue
        sub["X"] = (sub.cell_lab == a).astype(float)
        sub["log_ahi"] = np.log1p(sub.AHI)
        try:
            m0 = _fit(sub, ["X"])
            m1 = _fit(sub, ["X", "log_ahi"])
        except Exception:
            continue
        s0, s1 = m0.summary.loc["X"], m1.summary.loc["X"]
        rows_a.append({"contrast": f"{a}_vs_{ref}", "key": key, "disease": label,
                       "negative_control": key in NEGATIVE_CONTROLS,
                       "hr_unadj": round(float(s0["exp(coef)"]), 3),
                       "hr_ahi_adj": round(float(s1["exp(coef)"]), 3),
                       "lo_ahi_adj": round(float(s1["exp(coef) lower 95%"]), 3),
                       "hi_ahi_adj": round(float(s1["exp(coef) upper 95%"]), 3),
                       "pct_attenuation": round(100 * (1 - (float(s1["exp(coef)"]) - 1) /
                                                       (float(s0["exp(coef)"]) - 1)), 1)
                       if abs(float(s0["exp(coef)"]) - 1) > .02 else np.nan,
                       "events": ea + er})
ahiadj = pd.DataFrame(rows_a)
ahiadj.to_csv(f"{HERE}/attack_ahi_adjusted.csv", index=False)
for pre in ["cell2_vs_cell1", "cell1_vs_cell4"]:
    s = ahiadj[(ahiadj.contrast == pre) & ~ahiadj.negative_control]
    if not len(s):
        continue
    A3["ahi_adjustment"][pre] = {
        "n": int(len(s)),
        "median_hr_unadjusted": round(float(s.hr_unadj.median()), 3),
        "median_hr_ahi_adjusted": round(float(s.hr_ahi_adj.median()), 3),
        "n_sig_ahi_adjusted": int(((s.lo_ahi_adj > 1) | (s.hi_ahi_adj < 1)).sum()),
        "median_pct_attenuation": round(float(s.pct_attenuation.median()), 1),
        "headline": {r.disease: [r.hr_unadj, r.hr_ahi_adj, r.lo_ahi_adj, r.hi_ahi_adj]
                     for _, r in s[s.key.isin(HEAD)].iterrows()}}
    e = A3["ahi_adjustment"][pre]
    log(f"  {pre} adjusted for log apnoea-hypopnoea index: median HR {e['median_hr_unadjusted']} "
        f"-> {e['median_hr_ahi_adjusted']}, {e['n_sig_ahi_adjusted']} of {e['n']} still "
        f"significant, median attenuation {e['median_pct_attenuation']}%")
    for dz, v in e["headline"].items():
        log(f"      {dz:28s} {v[0]:.2f} -> {v[1]:.2f} [{v[2]:.2f}-{v[3]:.2f}]")


# ======================================================================================
# A4. SELECTION ON THE OUTCOME
# ======================================================================================
sec("A4. SELECTION ON THE OUTCOME. THE LOW-OXYGEN HEALTHY CELLS AS A SURVIVOR STRATUM")
A4 = OUT["attacks"]["A4_selection"] = {}

lowfull = int((b.spo2_pct_below_90 > T90_LOW).sum())
lowheal = int((h.spo2_pct_below_90 > T90_LOW).sum())
normfull = int((b.spo2_pct_below_90 <= T90_NORMAL).sum())
normheal = int((h.spo2_pct_below_90 <= T90_NORMAL).sum())
A4["survival_of_the_healthy_filter"] = {
    "low_oxygen_t90_over_10": {"full_cohort": lowfull, "healthy": lowheal,
                               "survival_pct": round(100 * lowheal / lowfull, 1)},
    "normal_oxygen_t90_le_1": {"full_cohort": normfull, "healthy": normheal,
                               "survival_pct": round(100 * normheal / normfull, 1)},
    "relative_risk_of_surviving": round((lowheal / lowfull) / (normheal / normfull), 3)}
e = A4["survival_of_the_healthy_filter"]
log(f"of the {lowfull:,} low-oxygen people in the full cohort, {lowheal:,} survive the healthy "
    f"filter ({e['low_oxygen_t90_over_10']['survival_pct']}%). Of the {normfull:,} "
    f"normal-oxygen people, {normheal:,} survive "
    f"({e['normal_oxygen_t90_le_1']['survival_pct']}%). A low-oxygen person is "
    f"{e['relative_risk_of_surviving']}x as likely to be called healthy.")

fullsizes = {k: int((b.cell == k).sum()) for k in CELLS}
A4["cell_survival"] = {k: {"full_cohort": fullsizes[k], "healthy": sizes[k],
                           "survival_pct": round(100 * sizes[k] / max(fullsizes[k], 1), 1)}
                       for k in CELLS}
for k in CELLS:
    e = A4["cell_survival"][k]
    log(f"  {k}: {e['healthy']:,} of {e['full_cohort']:,} survive = {e['survival_pct']}%")

# the T90 distribution inside the low-oxygen group, healthy vs excluded
lo_h = h[h.spo2_pct_below_90 > T90_LOW]
lo_x = b[(b.spo2_pct_below_90 > T90_LOW) & ~b.index.isin(h.index)]
A4["low_oxygen_healthy_vs_excluded"] = {
    "n_healthy": int(len(lo_h)), "n_excluded": int(len(lo_x)),
    "t90_median_healthy": round(float(lo_h.spo2_pct_below_90.median()), 2),
    "t90_median_excluded": round(float(lo_x.spo2_pct_below_90.median()), 2),
    "age_median_healthy": round(float(lo_h.AgeAtVisit.median()), 1),
    "age_median_excluded": round(float(lo_x.AgeAtVisit.median()), 1),
    "ahi_median_healthy": round(float(lo_h.AHI.median()), 1),
    "ahi_median_excluded": round(float(lo_x.AHI.median()), 1)}
e = A4["low_oxygen_healthy_vs_excluded"]
log(f"  within the low-oxygen group, the ones called healthy are younger "
    f"({e['age_median_healthy']} vs {e['age_median_excluded']}) with slightly lower T90 "
    f"({e['t90_median_healthy']}% vs {e['t90_median_excluded']}%) and lower AHI "
    f"({e['ahi_median_healthy']} vs {e['ahi_median_excluded']})")

# undiagnosed rather than well: do incident events cluster immediately after the sleep study?
early = []
for key, label in OUTCOMES:
    dd = frame_for(h, key)
    if dd is None:
        continue
    row = {"disease": label, "key": key}
    for k in CELLS:
        s = dd[(dd.cell_lab == k) & (dd.E == 1)]
        row[f"{k}_ev"] = int(len(s))
        row[f"{k}_pct_within_1y"] = (round(100 * float((s["T"] <= 1.0).mean()), 1)
                                     if len(s) >= 5 else np.nan)
        row[f"{k}_median_yrs"] = round(float(s["T"].median()), 2) if len(s) >= 5 else np.nan
    early.append(row)
early = pd.DataFrame(early)
early.to_csv(f"{HERE}/attack_event_timing.csv", index=False)
A4["diagnostic_catch_up"] = {
    k: {"n_outcomes_with_5plus_events": int(early[f"{k}_pct_within_1y"].notna().sum()),
        "median_pct_of_events_within_1y": round(float(early[f"{k}_pct_within_1y"].median()), 1),
        "median_time_to_event_yrs": round(float(early[f"{k}_median_yrs"].median()), 2)}
    for k in CELLS}
for k in CELLS:
    e = A4["diagnostic_catch_up"][k]
    log(f"  {k}: across outcomes with at least 5 events, a median of "
        f"{e['median_pct_of_events_within_1y']}% of incident diagnoses land inside the first year "
        f"(median time to event {e['median_time_to_event_yrs']} y)")

# and the direct test: refit the primary contrasts with the first year of follow-up discarded
LAND = 1.0
hl = h[h.death_years > LAND].copy()
A4["landmark_1y"] = {"n_retained": int(len(hl)),
                     "cells": {k: int((hl.cell == k).sum()) for k in CELLS}}
r_land = []
for key, label in OUTCOMES:
    dd = frame_for(h, key)
    if dd is None:
        continue
    dd = dd[dd["T"] > LAND].copy()
    dd["T"] = dd["T"] - LAND
    if not len(dd):
        continue
    rec = {"key": key, "disease": label, "negative_control": key in NEGATIVE_CONTROLS}
    for k2 in CELLS:
        rec[f"n_{k2}"] = int((dd.cell_lab == k2).sum())
        rec[f"ev_{k2}"] = int(dd[dd.cell_lab == k2].E.sum())
    for a, ref in [("cell1", "cell4"), ("cell2", "cell1")]:
        r = pairfit(dd, a, ref)
        pre = f"{a}_vs_{ref}"
        if r is None:
            rec[f"{pre}_hr"] = np.nan
            continue
        for kk, vv in r.items():
            rec[f"{pre}_{kk}"] = vv
    r_land.append(rec)
r_land = pd.DataFrame(r_land)
r_land.to_csv(f"{HERE}/attack_landmark1y.csv", index=False)
for pre in ["cell1_vs_cell4", "cell2_vs_cell1"]:
    A4["landmark_1y"][pre] = tally(r_land, pre)
    A4["landmark_1y"][f"{pre}_primary_for_comparison"] = tally(prim, pre)
    a, c = A4["landmark_1y"][pre], A4["landmark_1y"][f"{pre}_primary_for_comparison"]
    log(f"  landmark, first year discarded: {pre} {a.get('n_sig')} of {a.get('n')} significant, "
        f"median HR {a.get('median_hr')}  (primary was {c.get('n_sig')} of {c.get('n')}, "
        f"median HR {c.get('median_hr')})")

# is the loss of significance under the landmark attenuation or lost power? compare the point
# estimates and the event counts condition by condition.
cmpl = prim[["key", "disease", "cell2_vs_cell1_hr", "cell2_vs_cell1_ev_a",
             "cell2_vs_cell1_ev_ref"]].merge(
    r_land[["key", "cell2_vs_cell1_hr", "cell2_vs_cell1_lo", "cell2_vs_cell1_hi",
            "cell2_vs_cell1_ev_a", "cell2_vs_cell1_ev_ref"]], on="key", suffixes=("", "_L"))
cmpl = cmpl[cmpl.cell2_vs_cell1_hr.notna()]
cmpl["hr_ratio_landmark_over_primary"] = (cmpl.cell2_vs_cell1_hr_L /
                                          cmpl.cell2_vs_cell1_hr).round(3)
cmpl["cell2_events_retained_pct"] = (100 * cmpl.cell2_vs_cell1_ev_a_L /
                                     cmpl.cell2_vs_cell1_ev_a).round(1)
cmpl.to_csv(f"{HERE}/attack_landmark_vs_primary.csv", index=False)
A4["landmark_condition_by_condition"] = [
    {"disease": r.disease, "hr_primary": float(r.cell2_vs_cell1_hr),
     "hr_landmark": (None if pd.isna(r.cell2_vs_cell1_hr_L) else float(r.cell2_vs_cell1_hr_L)),
     "ev_cell2_primary": int(r.cell2_vs_cell1_ev_a),
     "ev_cell2_landmark": (None if pd.isna(r.cell2_vs_cell1_ev_a_L)
                           else int(r.cell2_vs_cell1_ev_a_L)),
     "not_estimable_after_landmark": bool(pd.isna(r.cell2_vs_cell1_hr_L))}
    for _, r in cmpl.iterrows()]
A4["landmark_median_cell2_events_retained_pct"] = round(float(
    cmpl.cell2_events_retained_pct.median()), 1)
log(f"  under the landmark, cell 2 keeps a median of "
    f"{A4['landmark_median_cell2_events_retained_pct']}% of its incident events")


# ======================================================================================
# A5. NEGATIVE CONTROLS
# ======================================================================================
sec("A5. NEGATIVE CONTROLS, PER CONTRAST, AND WHAT THEY STRIKE")
A5 = OUT["attacks"]["A5_negative_controls"] = {}
for pre in ["cell1_vs_cell4", "cell2_vs_cell4", "cell3_vs_cell4", "cell2_vs_cell1"]:
    A5[pre] = controls(prim, pre)
    e = A5[pre]
    log(f"  {pre}: {e['n_controls']} controls estimable, span {e.get('range')}, floor "
        f"{e.get('floor')}, {e.get('n_significant_controls')} of them significant")

struck = []
for pre in ["cell1_vs_cell4", "cell2_vs_cell4", "cell3_vs_cell4", "cell2_vs_cell1"]:
    fl = A5[pre].get("floor")
    if fl is None:
        continue
    s = prim[~prim.negative_control & prim[f"{pre}_hr"].notna()]
    s = s[(s[f"{pre}_lo"] > 1) | (s[f"{pre}_hi"] < 1)]
    for _, r in s.iterrows():
        hr = float(r[f"{pre}_hr"])
        eff = hr if hr >= 1 else 1 / hr
        struck.append({"contrast": pre, "disease": r.disease, "hr": hr,
                       "lo": float(r[f"{pre}_lo"]), "hi": float(r[f"{pre}_hi"]),
                       "floor": fl, "clears_floor": bool(eff > fl),
                       "verdict": "SURVIVES" if eff > fl else "STRUCK"})
struck = pd.DataFrame(struck)
struck.to_csv(f"{HERE}/attack_negcontrol_verdicts.csv", index=False)
A5["significant_results_vs_floor"] = struck.to_dict("records")
A5["n_struck"] = int((struck.verdict == "STRUCK").sum()) if len(struck) else 0
A5["n_survives"] = int((struck.verdict == "SURVIVES").sum()) if len(struck) else 0
log(f"  of {len(struck)} nominally significant non-control results across the four contrasts, "
    f"{A5['n_survives']} clear their floor and {A5['n_struck']} are struck")
for _, r in struck[struck.verdict == "STRUCK"].iterrows():
    log(f"    STRUCK  {r.contrast:18s} {r.disease:30s} HR {r.hr:.2f} against a floor of {r.floor}")


# ======================================================================================
# A6. THRESHOLD SENSITIVITY
# ======================================================================================
sec("A6. THRESHOLD SENSITIVITY, TST 4/5/6 h CROSSED WITH T90 5%/10%")
A6 = OUT["attacks"]["A6_thresholds"] = {"grid": []}
GRID = [(t, o) for t in (240.0, 300.0, 360.0) for o in (5.0, 10.0)]
for tst, o2 in GRID:
    tstnorm = max(TST_NORMAL, tst)
    d = h.assign(cell=assign_cells(h, tst, t90_low=o2, tst_norm=tstnorm))
    sz = {k: int((d.cell == k).sum()) for k in CELLS}
    r = run_panel(d, f"tst{int(tst)}_t90{int(o2)}",
                  pairs=[("cell1", "cell4"), ("cell2", "cell1"), ("cell2", "cell4"),
                         ("cell3", "cell4")])
    row = {"tst_short_lt_min": tst, "tst_normal_ge_min": tstnorm, "t90_low_gt_pct": o2,
           "t90_normal_le_pct": T90_NORMAL, "sizes": sz,
           "is_primary": bool(tst == TST_PRIMARY and o2 == T90_LOW),
           "cell1_vs_cell4": tally(r, "cell1_vs_cell4"),
           "cell2_vs_cell1": tally(r, "cell2_vs_cell1"),
           "cell2_vs_cell4": tally(r, "cell2_vs_cell4"),
           "cell3_vs_cell4": tally(r, "cell3_vs_cell4"),
           "negcontrol_floor_cell1_vs_cell4": controls(r, "cell1_vs_cell4").get("floor"),
           "negcontrol_floor_cell2_vs_cell1": controls(r, "cell2_vs_cell1").get("floor")}
    for kk in ["hf", "ckd", "obesity", "diabetes", "cvd", "death", "htn2"]:
        s = r[r.key == kk]
        if len(s) and not pd.isna(s.iloc[0].get("cell2_vs_cell1_hr", np.nan)):
            row[f"{kk}_c2vc1_hr"] = float(s.iloc[0]["cell2_vs_cell1_hr"])
    A6["grid"].append(row)
    a, c = row["cell1_vs_cell4"], row["cell2_vs_cell1"]
    log(f"  TST<{int(tst)} T90>{int(o2)}: cells {sz}  |  cell1vs4 {a.get('n_sig')}/{a.get('n')} "
        f"sig median HR {a.get('median_hr')}  |  cell2vs1 {c.get('n_sig')}/{c.get('n')} sig "
        f"median HR {c.get('median_hr')}")
pd.DataFrame([{**{k2: v for k2, v in g.items() if not isinstance(v, dict)},
               **{f"sz_{k2}": v for k2, v in g["sizes"].items()},
               **{f"c1v4_{k2}": v for k2, v in g["cell1_vs_cell4"].items()},
               **{f"c2v1_{k2}": v for k2, v in g["cell2_vs_cell1"].items()}}
              for g in A6["grid"]]).to_csv(f"{HERE}/attack_thresholds.csv", index=False)

c1v4 = [g["cell1_vs_cell4"] for g in A6["grid"]]
c2v1 = [g["cell2_vs_cell1"] for g in A6["grid"]]
A6["cell1_vs_cell4_conclusion_stable"] = all(
    (g.get("n_fdr", 0) == 0) for g in c1v4)
A6["cell1_vs_cell4_median_hr_range"] = [min(g["median_hr"] for g in c1v4),
                                        max(g["median_hr"] for g in c1v4)]
A6["cell2_vs_cell1_median_hr_range"] = [min(g["median_hr"] for g in c2v1 if g.get("n")),
                                        max(g["median_hr"] for g in c2v1 if g.get("n"))]
A6["cell2_vs_cell1_always_shows_excess"] = all(
    g["median_hr"] > 1 for g in c2v1 if g.get("n"))
A6["cell2_vs_cell1_n_grid_points_with_fdr_survivors"] = int(
    sum(1 for g in c2v1 if g.get("n_fdr", 0) > 0))
log(f"  cell1 vs cell4: zero FDR survivors at every one of the {len(GRID)} grid points: "
    f"{A6['cell1_vs_cell4_conclusion_stable']}; median HR ranges "
    f"{A6['cell1_vs_cell4_median_hr_range']}")
log(f"  cell2 vs cell1: median HR above 1 at every grid point: "
    f"{A6['cell2_vs_cell1_always_shows_excess']}; range "
    f"{A6['cell2_vs_cell1_median_hr_range']}; FDR survivors at "
    f"{A6['cell2_vs_cell1_n_grid_points_with_fdr_survivors']} of {len(GRID)} points")


# ======================================================================================
# GUARDS
# ======================================================================================
sec("GUARDS")
OUT["guards"] = {"n_fits": len(SPLINE_MAXABS),
                 "min_max_abs_age_spline_coef": round(float(np.min(SPLINE_MAXABS)), 4),
                 "n_fits_with_all_zero_age_coefs": int(np.sum(np.array(SPLINE_MAXABS) < 1e-8))}
log(f"{OUT['guards']['n_fits']:,} Cox fits, smallest max |age-spline coefficient| "
    f"{OUT['guards']['min_max_abs_age_spline_coef']}, fits returning all-zero age coefficients: "
    f"{OUT['guards']['n_fits_with_all_zero_age_coefs']}")
assert OUT["guards"]["n_fits_with_all_zero_age_coefs"] == 0, "an age spline collapsed to zero"

OUT["log"] = LOG
OUT["runtime_min"] = round((time.time() - T0) / 60, 1)
json.dump(OUT, open(f"{HERE}/attack.json", "w"), indent=1, default=str)
log(f"\nwrote attack.json in {OUT['runtime_min']} min")
