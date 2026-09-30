"""
Sleep duration crossed with nocturnal oxygenation, six cells.

The earlier 2x2 (../dur_x_oxygen/model.py) folded every sleeper of six hours or more into one
"normal" group, so long sleepers were invisible. This extends it to three duration strata by two
oxygenation strata, which is what Alen asked for: short sleep with good oxygen, short sleep with
bad oxygen, long sleep with good oxygen, long sleep with bad oxygen, all against a six-to-seven
hour normally oxygenated reference.

Cells, on the 19,173 analysis cohort:
    SN  short sleep  + normal oxygen   TST < 240 min,        T90 <= 1%
    SL  short sleep  + low oxygen      TST < 240 min,        T90 > 10%
    NN  normal sleep + normal oxygen   TST 360-420 min,      T90 <= 1%    REFERENCE
    NL  normal sleep + low oxygen      TST 360-420 min,      T90 > 10%
    LN  long sleep   + normal oxygen   TST > 420 min,        T90 <= 1%
    LL  long sleep   + low oxygen      TST > 420 min,        T90 > 10%

The 1-10% T90 band and the 240-360 min TST band are dropped so the contrasts are clean. Every
discarded count is reported.

TST is one laboratory night, not habitual sleep duration, so a short laboratory TST can mean a
short recording rather than little sleep. Every model is repeated on recordings of at least
360 minutes, and the cell profile reports recording length, sleep efficiency and the AHI so a
reader can see what these people actually are.

Modelling: site-stratified Cox, unpenalized (penalizer=0), natural cubic age spline with 4 df of
which one column is always dropped so the design is full rank, plus sex. PSG date is time zero.
Circular conditions are excluded from the outcome panel.

Writes cells.csv, results.csv, interaction.csv, summary.json in this directory.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import json
import os
import sys
import warnings

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from patsy import dmatrix
from scipy import stats
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")

T90ROOT = paths.T90_ROOT
HERE = f"{paths.SV_ROOT}/dur_x_oxygen_6cell"
os.makedirs(HERE, exist_ok=True)
sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
from cohort_spec import apply_cohort                          # noqa: E402
from disease_definitions import (DISEASES, NEGATIVE_CONTROLS,  # noqa: E402
                                 CIRCULAR, ORGAN_GROUP)

T90_LOW, T90_NORMAL = 10.0, 1.0            # low oxygen T90 > 10%, normal oxygen T90 <= 1%
TST_SHORT = 240.0                          # short sleep  < 4 h
TST_NORM_LO, TST_NORM_HI = 360.0, 420.0    # normal sleep 6 to 7 h, inclusive both ends
MIN_EVENTS = 20                            # per condition, summed over the six cells
REC_MIN = 360.0                            # sensitivity: recordings of at least 6 h
Z_A, Z_B = 1.959963985, 0.841621234        # alpha 0.05 two-sided, 80% power

OUT = {}
CELLS = ["SN", "SL", "NN", "NL", "LN", "LL"]
NONREF = ["SN", "SL", "NL", "LN", "LL"]
REF = "NN"
CELL_LABEL = {
    "SN": "short sleep (<4h), normal oxygen (T90<=1%)",
    "SL": "short sleep (<4h), low oxygen (T90>10%)",
    "NN": "normal sleep (6-7h), normal oxygen (T90<=1%)  REFERENCE",
    "NL": "normal sleep (6-7h), low oxygen (T90>10%)",
    "LN": "long sleep (>7h), normal oxygen (T90<=1%)",
    "LL": "long sleep (>7h), low oxygen (T90>10%)",
}
# within-duration-stratum oxygen contrasts: (name, low-oxygen cell, normal-oxygen cell)
OX_CONTRASTS = [("ox_short", "SL", "SN"), ("ox_normal", "NL", "NN"), ("ox_long", "LL", "LN")]


def log(m):
    print(m, flush=True)


# ======================================================================================
# 1. COHORT AND CELLS
# ======================================================================================
b = apply_cohort(pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet"))
log(f"cohort {len(b):,}")

b["male"] = (b.sex.astype(str).str.upper().str[0] == "M").astype(float)
b.loc[~b.sex.astype(str).str.upper().str[0].isin(["M", "F"]), "male"] = np.nan

lowO2 = b.spo2_pct_below_90 > T90_LOW
normO2 = b.spo2_pct_below_90 <= T90_NORMAL
shortS = b.TST_min < TST_SHORT
normS = (b.TST_min >= TST_NORM_LO) & (b.TST_min <= TST_NORM_HI)
longS = b.TST_min > TST_NORM_HI

b["cell"] = np.select(
    [shortS & normO2, shortS & lowO2, normS & normO2, normS & lowO2, longS & normO2, longS & lowO2],
    CELLS, default="excluded")

OUT["thresholds"] = {
    "t90_low_gt_pct": T90_LOW, "t90_normal_le_pct": T90_NORMAL,
    "tst_short_lt_min": TST_SHORT,
    "tst_normal_min": TST_NORM_LO, "tst_normal_max": TST_NORM_HI,
    "tst_long_gt_min": TST_NORM_HI,
    "reference_cell": REF, "penalizer": 0.0, "age_spline_df": 4, "age_spline_cols_used": 3,
}
OUT["discards"] = {
    "cohort": int(len(b)),
    "t90_middle_1_to_10_pct": int(((b.spo2_pct_below_90 > T90_NORMAL) &
                                   (b.spo2_pct_below_90 <= T90_LOW)).sum()),
    "tst_middle_240_to_360_min": int(((b.TST_min >= TST_SHORT) & (b.TST_min < TST_NORM_LO)).sum()),
    "tst_missing": int(b.TST_min.isna().sum()),
    "in_six_cells": int((b.cell != "excluded").sum()),
    "excluded_total": int((b.cell == "excluded").sum()),
}
OUT["duration_strata_totals"] = {
    "short_lt240": int(shortS.sum()), "normal_360_420": int(normS.sum()),
    "long_gt420": int(longS.sum()),
}
# how much the reference shrinks when normal sleep gains an upper bound
prior_ref = int(((b.TST_min >= TST_NORM_LO) & normO2).sum())
prior_nl = int(((b.TST_min >= TST_NORM_LO) & lowO2).sum())
OUT["reference_redefinition"] = {
    "prior_normal_rule": "TST >= 360 min, no upper bound",
    "prior_n_normal_normalO2": prior_ref, "prior_n_normal_lowO2": prior_nl,
    "new_n_normal_normalO2": int((b.cell == "NN").sum()),
    "new_n_normal_lowO2": int((b.cell == "NL").sum()),
    "moved_to_long_normalO2": prior_ref - int((b.cell == "NN").sum()),
    "moved_to_long_lowO2": prior_nl - int((b.cell == "NL").sum()),
}

counts = {k: int((b.cell == k).sum()) for k in CELLS}
OUT["cell_n"] = counts
log("cells: " + ", ".join(f"{k}={v:,}" for k, v in counts.items()))
log(f"discarded {OUT['discards']['excluded_total']:,} "
    f"({OUT['discards']['t90_middle_1_to_10_pct']:,} in the 1-10% T90 band, "
    f"{OUT['discards']['tst_middle_240_to_360_min']:,} in the 240-360 min TST band, "
    f"{OUT['discards']['tst_missing']} with no TST). "
    f"{OUT['discards']['in_six_cells']:,} kept.")
log(f"reference shrinks {prior_ref:,} -> {counts['NN']:,} once normal sleep is capped at 7 h")

f6 = b[b.cell != "excluded"].copy()

# ======================================================================================
# 2. CELL PROFILE  ->  cells.csv
# ======================================================================================
PREV_CHECK = ["copd2_prevalent", "hf_prevalent", "cancer_any_prevalent", "obesity_prevalent",
              "osa_prevalent", "dementia_prevalent", "diabetes_prevalent"]


def q(s, p):
    return round(float(s.quantile(p)), 2)


prof = []
for k in CELLS:
    s = f6[f6.cell == k]
    row = {
        "cell": k, "definition": CELL_LABEL[k], "n": int(len(s)),
        "pct_of_six_cells": round(100 * len(s) / len(f6), 1),
        "age_median": q(s.AgeAtVisit, .5), "age_q1": q(s.AgeAtVisit, .25),
        "age_q3": q(s.AgeAtVisit, .75),
        "pct_male": round(100 * float(s.male.mean()), 1),
        "t90_median": q(s.spo2_pct_below_90, .5), "t90_q1": q(s.spo2_pct_below_90, .25),
        "t90_q3": q(s.spo2_pct_below_90, .75),
        "tst_min_median": q(s.TST_min, .5), "tst_min_q1": q(s.TST_min, .25),
        "tst_min_q3": q(s.TST_min, .75),
        "recording_dur_min_median": q(s.recording_dur_min, .5),
        "recording_dur_min_q1": q(s.recording_dur_min, .25),
        "recording_dur_min_q3": q(s.recording_dur_min, .75),
        "pct_recording_lt360min": round(100 * float((s.recording_dur_min < REC_MIN).mean()), 1),
        "n_recording_ge360": int((s.recording_dur_min >= REC_MIN).sum()),
        "sleep_eff_pct_median": q(s.sleep_efficiency_pct, .5),
        "sleep_eff_pct_q1": q(s.sleep_efficiency_pct, .25),
        "sleep_eff_pct_q3": q(s.sleep_efficiency_pct, .75),
        "ahi_median": q(s.AHI, .5), "ahi_q1": q(s.AHI, .25), "ahi_q3": q(s.AHI, .75),
        "pct_ahi_ge15": round(100 * float((s.AHI >= 15).mean()), 1),
        "pct_ahi_ge30": round(100 * float((s.AHI >= 30).mean()), 1),
        "arousal_index_median": q(s.arousal_index, .5),
        "spo2_mean_median": q(s.spo2_mean, .5),
        "n_deaths": int(s.death_incident.sum()),
        "death_rate_pct": round(100 * float(s.death_incident.mean()), 2),
        "followup_years_median": q(s.death_years, .5),
    }
    for c in PREV_CHECK:
        if c in s.columns:
            row[c.replace("_prevalent", "_prev_pct")] = round(100 * float(s[c].mean()), 1)
    for site in sorted(f6.site_id.unique()):
        row[f"site_{site}"] = int((s.site_id == site).sum())
    prof.append(row)
cells = pd.DataFrame(prof)
cells.to_csv(f"{HERE}/cells.csv", index=False)
OUT["cell_profile"] = json.loads(cells.to_json(orient="records"))
log("\ncell profile")
log(cells[["cell", "n", "age_median", "pct_male", "tst_min_median", "recording_dur_min_median",
           "sleep_eff_pct_median", "t90_median", "ahi_median", "pct_ahi_ge15",
           "n_deaths"]].to_string(index=False))

# ======================================================================================
# 3. MODEL MACHINERY
# ======================================================================================
def agespline(d):
    """4-df natural cubic spline on age, first column dropped so the design is full rank.

    The four basis columns sum to one. A Cox partial likelihood has no intercept so the
    redundancy is harmless under a penalty, but an unpenalized fit on a rank-deficient design
    can return coefficients of zero without failing. One column is always dropped."""
    sp = dmatrix("cr(a, df=4) - 1", {"a": d.AgeAtVisit.values}, return_type="dataframe")
    sp.columns = [f"age_s{i}" for i in range(4)]
    sp.index = d.index
    return sp.iloc[:, 1:]


OUTCOMES = [(k, v[0]) for k, v in DISEASES.items() if k not in CIRCULAR]
OUTCOMES.append(("death", "Death from any cause"))


def frame_for(d, key):
    """Incident frame for one condition: no prevalent disease, positive follow-up."""
    yc, ec, pc = f"{key}_years", f"{key}_incident", f"{key}_prevalent"
    if not all(c in d.columns for c in (yc, ec, pc)):
        return None
    f = d[(d[pc] == 0) & d[yc].notna() & (d[yc] > 0)].copy()
    if not len(f):
        return None
    sp = agespline(f)
    cols = {"T": f[yc], "E": f[ec].astype(int), "site": f.site_id, "male": f.male,
            "cell_lab": f.cell}
    for c in sp.columns:
        cols[c] = sp[c]
    out = pd.DataFrame(cols).dropna()
    return out


def mdhr(events, n_cell, n_total):
    """Minimum detectable hazard ratio, 80% power at alpha 0.05, Schoenfeld.

    |log HR| = (z_a + z_b) / sqrt(D * p * (1-p)) with D the events and p the exposed fraction."""
    if events <= 0 or n_cell <= 0 or n_total <= 0 or n_cell >= n_total:
        return None
    p = n_cell / n_total
    return round(float(np.exp((Z_A + Z_B) / np.sqrt(events * p * (1 - p)))), 3)


def fit(dd, terms):
    cov = terms + ["male"] + [c for c in dd.columns if c.startswith("age_s")]
    return CoxPHFitter(penalizer=0.0).fit(dd[cov + ["T", "E", "site"]], "T", "E", strata=["site"])


def est(m, a, ref=None):
    """HR and 95% CI for term a, against the model baseline or against another term."""
    beta, V = m.params_, m.variance_matrix_
    if ref is None:
        d = float(beta[a])
        se = float(np.sqrt(V.loc[a, a]))
    else:
        d = float(beta[a] - beta[ref])
        se = float(np.sqrt(V.loc[a, a] + V.loc[ref, ref] - 2 * V.loc[a, ref]))
    z = d / se
    return (round(float(np.exp(d)), 3), round(float(np.exp(d - Z_A * se)), 3),
            round(float(np.exp(d + Z_A * se)), 3), float(2 * stats.norm.sf(abs(z))), se)


def run_panel(d, tag):
    """Every non-circular condition with at least MIN_EVENTS events across the six cells."""
    rows, skipped = [], []
    for key, label in OUTCOMES:
        dd = frame_for(d, key)
        if dd is None:
            skipped.append({"key": key, "reason": "columns absent"})
            continue
        if int(dd.E.sum()) < MIN_EVENTS:
            skipped.append({"key": key, "reason": f"{int(dd.E.sum())} events, below {MIN_EVENTS}"})
            continue
        ev = {k: int(dd[dd.cell_lab == k].E.sum()) for k in CELLS}
        nn = {k: int((dd.cell_lab == k).sum()) for k in CELLS}
        if ev[REF] == 0:
            skipped.append({"key": key, "reason": "reference cell has no events", "events": ev})
            continue
        # a cell with zero events has a maximum-likelihood coefficient of minus infinity, i.e.
        # zero hazard, so its members leave every risk set. Dropping them reproduces that fit
        # exactly and keeps the rest of the panel estimable.
        keep = [k for k in NONREF if ev[k] > 0]
        drop = [k for k in NONREF if ev[k] == 0]
        sub = dd[dd.cell_lab.isin(keep + [REF])].copy()
        for k in keep:
            sub[k] = (sub.cell_lab == k).astype(float)
        try:
            m = fit(sub, keep)
        except Exception as ex:
            skipped.append({"key": key, "reason": f"fit failed: {type(ex).__name__}: {ex}"})
            continue
        rec = {"analysis": tag, "key": key, "disease": label,
               "organ_group": ORGAN_GROUP.get(key, "Death"),
               "negative_control": key in NEGATIVE_CONTROLS,
               # obesity hypoventilation is defined partly by nocturnal desaturation, so an
               # association with T90 is mechanical. Kept in the panel, flagged, never headlined.
               "circular_flag": key == "obesity_hypovent",
               "n_total_six_cells": int(len(dd)), "events_total_six_cells": int(dd.E.sum()),
               "n_in_fit": int(len(sub)), "events_in_fit": int(sub.E.sum()),
               "cells_not_estimable": ",".join(drop)}
        for k in CELLS:
            rec[f"n_{k}"] = nn[k]
            rec[f"ev_{k}"] = ev[k]
            rec[f"mdhr_{k}"] = mdhr(int(dd.E.sum()), nn[k], int(len(dd))) if k != REF else None
        rec["min_cell_events"] = int(min(ev.values()))
        rec["unstable_lt5_events"] = bool(min(ev.values()) < 5)
        for k in NONREF:
            if k in keep:
                hr, lo, hi, p, _ = est(m, k)
                rec[f"{k}_hr"], rec[f"{k}_lo"] = hr, lo
                rec[f"{k}_hi"], rec[f"{k}_p"] = hi, p
            else:
                rec[f"{k}_hr"] = rec[f"{k}_lo"] = rec[f"{k}_hi"] = rec[f"{k}_p"] = np.nan
        # oxygen effect inside each duration stratum
        for name, hi_cell, lo_cell in OX_CONTRASTS:
            if hi_cell not in keep:
                rec[f"{name}_hr"] = rec[f"{name}_lo"] = np.nan
                rec[f"{name}_hi"] = rec[f"{name}_p"] = np.nan
                continue
            if lo_cell == REF:
                hr, lo, hi, p, _ = est(m, hi_cell)
            elif lo_cell not in keep:
                rec[f"{name}_hr"] = rec[f"{name}_lo"] = np.nan
                rec[f"{name}_hi"] = rec[f"{name}_p"] = np.nan
                continue
            else:
                hr, lo, hi, p, _ = est(m, hi_cell, lo_cell)
            rec[f"{name}_hr"], rec[f"{name}_lo"] = hr, lo
            rec[f"{name}_hi"], rec[f"{name}_p"] = hi, p
        # joint Wald test that all estimable non-reference cells equal the reference
        try:
            bvec = m.params_[keep].values
            Vk = m.variance_matrix_.loc[keep, keep].values
            chi2 = float(bvec @ np.linalg.solve(Vk, bvec))
            rec["global_chi2"] = round(chi2, 3)
            rec["global_df"] = len(keep)
            rec["global_p"] = float(stats.chi2.sf(chi2, len(keep)))
        except Exception:
            rec["global_chi2"] = rec["global_df"] = rec["global_p"] = np.nan
        rows.append(rec)
    df = pd.DataFrame(rows)
    if len(df):
        real = ~df.negative_control
        for k in NONREF + [c[0] for c in OX_CONTRASTS]:
            df[f"{k}_q"] = np.nan
            sel = real & df[f"{k}_p"].notna()
            if sel.sum():
                df.loc[sel, f"{k}_q"] = multipletests(df.loc[sel, f"{k}_p"],
                                                      method="fdr_bh")[1]
    log(f"\n[{tag}] {len(df)} conditions fitted, {len(skipped)} skipped")
    for s in skipped:
        log(f"    skip {s['key']}: {s['reason']}")
    return df, skipped


def floor_of(df):
    """Confounding floor per cell: the largest negative-control deviation from 1."""
    neg = df[df.negative_control]
    out = {}
    for k in NONREF + [c[0] for c in OX_CONTRASTS]:
        h = neg[f"{k}_hr"].dropna()
        if not len(h):
            out[k] = None
            continue
        dev = h.apply(lambda x: max(x, 1 / x))
        i = dev.idxmax()
        out[k] = {"max_abs_hr": round(float(dev.max()), 3),
                  "at_condition": str(neg.loc[i, "disease"]),
                  "hr_at": float(neg.loc[i, f"{k}_hr"]),
                  "ci_at": [float(neg.loc[i, f"{k}_lo"]), float(neg.loc[i, f"{k}_hi"])],
                  "range": [round(float(h.min()), 3), round(float(h.max()), 3)],
                  "n_controls": int(len(h)),
                  "n_controls_significant": int(((neg[f"{k}_lo"] > 1) |
                                                 (neg[f"{k}_hi"] < 1)).sum()),
                  "controls": {str(r.disease): [float(r[f"{k}_hr"]), float(r[f"{k}_lo"]),
                                                float(r[f"{k}_hi"])]
                               for _, r in neg.iterrows() if pd.notna(r[f"{k}_hr"])}}
    return out


# ======================================================================================
# 4. PRIMARY PANEL
# ======================================================================================
prim, prim_skip = run_panel(f6, "primary")
OUT["skipped_primary"] = prim_skip
FLOOR = floor_of(prim)
OUT["negative_control_floor_primary"] = FLOOR

# ======================================================================================
# 5. SENSITIVITY: RECORDINGS OF AT LEAST 6 HOURS
# ======================================================================================
f6r = f6[f6.recording_dur_min >= REC_MIN].copy()
OUT["long_recording_subset"] = {
    "rule": f"recording_dur_min >= {REC_MIN:.0f}",
    "n": int(len(f6r)),
    "by_cell": {k: int((f6r.cell == k).sum()) for k in CELLS},
    "pct_retained_by_cell": {k: round(100 * float((f6r.cell == k).sum() /
                                                  max((f6.cell == k).sum(), 1)), 1)
                             for k in CELLS},
}
log("\nlong-recording subset " + json.dumps(OUT["long_recording_subset"]["by_cell"]))
sens, sens_skip = run_panel(f6r, "recording_ge_360min")
OUT["skipped_sensitivity"] = sens_skip
OUT["negative_control_floor_sensitivity"] = floor_of(sens)

# ======================================================================================
# 6. FORMAL INTERACTION: DOES THE OXYGEN HR DEPEND ON SLEEP DURATION
# ======================================================================================
top10 = prim.sort_values("events_total_six_cells", ascending=False).head(10).key.tolist()
inter = []
for key in top10:
    dd = frame_for(f6, key)
    dd["dur_short"] = (dd.cell_lab.isin(["SN", "SL"])).astype(float)
    dd["dur_long"] = (dd.cell_lab.isin(["LN", "LL"])).astype(float)
    dd["lowo2"] = (dd.cell_lab.isin(["SL", "NL", "LL"])).astype(float)
    dd["short_x_lowo2"] = dd.dur_short * dd.lowo2
    dd["long_x_lowo2"] = dd.dur_long * dd.lowo2
    ev = {k: int(dd[dd.cell_lab == k].E.sum()) for k in CELLS}
    terms = ["dur_short", "dur_long", "lowo2", "short_x_lowo2", "long_x_lowo2"]
    row = {"key": key, "disease": str(prim.loc[prim.key == key, "disease"].iloc[0]),
           "events": int(dd.E.sum()),
           **{f"ev_{k}": ev[k] for k in CELLS}}
    if min(ev.values()) == 0:
        row.update({"fitted": False, "reason": "a cell has no events",
                    "interaction_joint_p": np.nan, "interaction_joint_chi2": np.nan})
        inter.append(row)
        continue
    m = fit(dd, terms)
    s = m.summary
    ip = ["short_x_lowo2", "long_x_lowo2"]
    bvec = m.params_[ip].values
    Vk = m.variance_matrix_.loc[ip, ip].values
    chi2 = float(bvec @ np.linalg.solve(Vk, bvec))
    row.update({
        "fitted": True,
        "main_short_hr": round(float(s.loc["dur_short", "exp(coef)"]), 3),
        "main_long_hr": round(float(s.loc["dur_long", "exp(coef)"]), 3),
        "main_lowo2_hr": round(float(s.loc["lowo2", "exp(coef)"]), 3),
        "main_lowo2_lo": round(float(s.loc["lowo2", "exp(coef) lower 95%"]), 3),
        "main_lowo2_hi": round(float(s.loc["lowo2", "exp(coef) upper 95%"]), 3),
        "short_x_lowo2_hr": round(float(s.loc["short_x_lowo2", "exp(coef)"]), 3),
        "short_x_lowo2_lo": round(float(s.loc["short_x_lowo2", "exp(coef) lower 95%"]), 3),
        "short_x_lowo2_hi": round(float(s.loc["short_x_lowo2", "exp(coef) upper 95%"]), 3),
        "short_x_lowo2_p": float(s.loc["short_x_lowo2", "p"]),
        "long_x_lowo2_hr": round(float(s.loc["long_x_lowo2", "exp(coef)"]), 3),
        "long_x_lowo2_lo": round(float(s.loc["long_x_lowo2", "exp(coef) lower 95%"]), 3),
        "long_x_lowo2_hi": round(float(s.loc["long_x_lowo2", "exp(coef) upper 95%"]), 3),
        "long_x_lowo2_p": float(s.loc["long_x_lowo2", "p"]),
        "interaction_joint_chi2": round(chi2, 3),
        "interaction_joint_df": 2,
        "interaction_joint_p": float(stats.chi2.sf(chi2, 2)),
        "ox_short_hr": float(prim.loc[prim.key == key, "ox_short_hr"].iloc[0]),
        "ox_normal_hr": float(prim.loc[prim.key == key, "ox_normal_hr"].iloc[0]),
        "ox_long_hr": float(prim.loc[prim.key == key, "ox_long_hr"].iloc[0]),
    })
    inter.append(row)
inter = pd.DataFrame(inter)
ok = inter.interaction_joint_p.notna()
inter["interaction_joint_q"] = np.nan
if ok.sum():
    inter.loc[ok, "interaction_joint_q"] = multipletests(inter.loc[ok, "interaction_joint_p"],
                                                         method="fdr_bh")[1]
inter.to_csv(f"{HERE}/interaction.csv", index=False)
OUT["interaction_top10"] = json.loads(inter.to_json(orient="records"))
log("\ninteraction: does the oxygen HR depend on sleep duration (top 10 by events)")
log(inter[["disease", "events", "ox_short_hr", "ox_normal_hr", "ox_long_hr",
           "interaction_joint_chi2", "interaction_joint_p"]].to_string(index=False))

# ======================================================================================
# 7. SUMMARIES
# ======================================================================================
def sigrows(df, k):
    r = df[~df.negative_control]
    return r[(r[f"{k}_lo"] > 1) | (r[f"{k}_hi"] < 1)]


def block(df, k, floor):
    r = df[~df.negative_control & df[f"{k}_hr"].notna()]
    rc = r[~r.circular_flag]
    s = sigrows(df, k)
    fl = floor.get(k) or {}
    thr = fl.get("max_abs_hr", 1.0)
    return {
        "label": CELL_LABEL.get(k, k),
        "n_conditions": int(len(r)),
        "n_significant": int(len(s)),
        "n_significant_excl_circular": int((~s.circular_flag).sum()),
        "n_expected_by_chance": round(0.05 * len(r), 1),
        "n_significant_fdr": int((r[f"{k}_q"] < 0.05).sum()),
        "median_hr": round(float(r[f"{k}_hr"].median()), 3) if len(r) else None,
        "median_hr_excl_circular": (round(float(rc[f"{k}_hr"].median()), 3) if len(rc) else None),
        "max_hr_excl_circular": (round(float(rc[f"{k}_hr"].max()), 3) if len(rc) else None),
        "hr_range": [round(float(r[f"{k}_hr"].min()), 3),
                     round(float(r[f"{k}_hr"].max()), 3)] if len(r) else None,
        "confounding_floor": thr,
        "n_clearing_floor": int(sum(max(x, 1 / x) > thr for x in s[f"{k}_hr"])),
        "significant": [{"disease": str(rr.disease), "hr": float(rr[f"{k}_hr"]),
                         "ci": [float(rr[f"{k}_lo"]), float(rr[f"{k}_hi"])],
                         "p": float(rr[f"{k}_p"]),
                         "q": (None if pd.isna(rr[f"{k}_q"]) else round(float(rr[f"{k}_q"]), 4)),
                         "events": int(rr.events_total_six_cells),
                         "clears_floor": bool(max(rr[f"{k}_hr"], 1 / rr[f"{k}_hr"]) > thr)}
                        for _, rr in s.sort_values(f"{k}_hr", ascending=False).iterrows()],
    }


OUT["cells_vs_reference"] = {k: block(prim, k, FLOOR) for k in NONREF}
OUT["oxygen_within_duration"] = {c[0]: block(prim, c[0], FLOOR) for c in OX_CONTRASTS}
OUT["oxygen_within_duration"]["ox_short"]["label"] = "low vs normal oxygen among short sleepers"
OUT["oxygen_within_duration"]["ox_normal"]["label"] = "low vs normal oxygen among 6-7 h sleepers"
OUT["oxygen_within_duration"]["ox_long"]["label"] = "low vs normal oxygen among long sleepers"

# side-by-side oxygen effect across the three duration strata
real = prim[~prim.negative_control]
OUT["oxygen_effect_across_strata"] = {
    "median_hr": {c[0]: (round(float(real[f"{c[0]}_hr"].median()), 3)
                         if real[f"{c[0]}_hr"].notna().any() else None) for c in OX_CONTRASTS},
    "n_significant": {c[0]: int(((real[f"{c[0]}_lo"] > 1) | (real[f"{c[0]}_hi"] < 1)).sum())
                      for c in OX_CONTRASTS},
    "n_estimable": {c[0]: int(real[f"{c[0]}_hr"].notna().sum()) for c in OX_CONTRASTS},
    "per_condition": [{"disease": str(r.disease),
                       "events": int(r.events_total_six_cells),
                       "ox_short": [r.ox_short_hr, r.ox_short_lo, r.ox_short_hi],
                       "ox_normal": [r.ox_normal_hr, r.ox_normal_lo, r.ox_normal_hi],
                       "ox_long": [r.ox_long_hr, r.ox_long_lo, r.ox_long_hi]}
                      for _, r in real.sort_values("events_total_six_cells",
                                                   ascending=False).iterrows()],
}
# paired comparison of the oxygen HR across strata, on conditions estimable in all three
paired = real.dropna(subset=["ox_short_hr", "ox_normal_hr", "ox_long_hr"])
if len(paired) >= 3:
    OUT["oxygen_effect_across_strata"]["paired"] = {
        "n_conditions": int(len(paired)),
        "median_ox_short": round(float(paired.ox_short_hr.median()), 3),
        "median_ox_normal": round(float(paired.ox_normal_hr.median()), 3),
        "median_ox_long": round(float(paired.ox_long_hr.median()), 3),
        "wilcoxon_short_vs_normal_p": float(stats.wilcoxon(np.log(paired.ox_short_hr),
                                                           np.log(paired.ox_normal_hr)).pvalue),
        "wilcoxon_long_vs_normal_p": float(stats.wilcoxon(np.log(paired.ox_long_hr),
                                                          np.log(paired.ox_normal_hr)).pvalue),
    }


def headline(k):
    s = (OUT["cells_vs_reference"] if k in NONREF else OUT["oxygen_within_duration"])[k]
    return {"n_significant": s["n_significant"],
            "n_significant_excl_circular": s["n_significant_excl_circular"],
            "n_conditions": s["n_conditions"],
            "n_expected_by_chance": s["n_expected_by_chance"],
            "n_significant_fdr": s["n_significant_fdr"],
            "median_hr": s["median_hr"],
            "median_hr_excl_circular": s["median_hr_excl_circular"],
            "confounding_floor": s["confounding_floor"],
            "n_clearing_floor": s["n_clearing_floor"],
            "top5": [{"disease": x["disease"], "hr": x["hr"], "ci": x["ci"]}
                     for x in s["significant"][:5]]}


OUT["answers"] = {
    "a_short_normal_oxygen_vs_reference": {"cell": "SN", "n": counts["SN"], **headline("SN")},
    "b_short_low_oxygen_vs_reference": {"cell": "SL", "n": counts["SL"], **headline("SL")},
    "c_long_normal_oxygen_vs_reference": {"cell": "LN", "n": counts["LN"], **headline("LN")},
    "d_long_low_oxygen_vs_reference": {"cell": "LL", "n": counts["LL"], **headline("LL")},
    "e_oxygen_effect_by_duration": {c[0]: headline(c[0]) for c in OX_CONTRASTS},
    "f_interaction": {
        "n_tested": int(inter.interaction_joint_p.notna().sum()),
        "min_p": (float(inter.interaction_joint_p.min())
                  if inter.interaction_joint_p.notna().any() else None),
        "max_p": (float(inter.interaction_joint_p.max())
                  if inter.interaction_joint_p.notna().any() else None),
        "n_p_lt_005": int((inter.interaction_joint_p < 0.05).sum()),
        "n_q_lt_005": int((inter.interaction_joint_q < 0.05).sum()),
    },
}

# minimum detectable HR per cell, on the largest outcome and the median across the panel
OUT["minimum_detectable_hr"] = {
    k: {"median_across_conditions": (round(float(prim[f"mdhr_{k}"].median()), 3)
                                     if prim[f"mdhr_{k}"].notna().any() else None),
        "best_case_death": (float(prim.loc[prim.key == "death", f"mdhr_{k}"].iloc[0])
                            if (prim.key == "death").any() else None),
        "worst_case": (round(float(prim[f"mdhr_{k}"].max()), 3)
                       if prim[f"mdhr_{k}"].notna().any() else None),
        "n_conditions_with_lt5_events": int((prim[f"ev_{k}"] < 5).sum()),
        "n_conditions_with_0_events": int((prim[f"ev_{k}"] == 0).sum())}
    for k in NONREF}

# sensitivity agreement
mg = prim.merge(sens, on="key", suffixes=("_p", "_s"))
OUT["sensitivity_agreement"] = {}
for k in NONREF + [c[0] for c in OX_CONTRASTS]:
    a = (mg[f"{k}_lo_p"] > 1) | (mg[f"{k}_hi_p"] < 1)
    c = (mg[f"{k}_lo_s"] > 1) | (mg[f"{k}_hi_s"] < 1)
    both = mg[f"{k}_hr_p"].notna() & mg[f"{k}_hr_s"].notna()
    OUT["sensitivity_agreement"][k] = {
        "n_paired": int(both.sum()),
        "sig_primary": int(a.sum()), "sig_sensitivity": int(c.sum()),
        "sig_both": int((a & c).sum()),
        "median_abs_pct_change_hr": (round(float(
            (100 * (mg.loc[both, f"{k}_hr_s"] / mg.loc[both, f"{k}_hr_p"] - 1)).abs().median()), 1)
            if both.sum() else None),
        "spearman_hr": (round(float(stats.spearmanr(mg.loc[both, f"{k}_hr_p"],
                                                    mg.loc[both, f"{k}_hr_s"]).statistic), 3)
                        if both.sum() > 2 else None),
        "lost": [str(x) for x in mg.loc[a & ~c, "disease_p"]],
        "gained": [str(x) for x in mg.loc[~a & c, "disease_p"]],
    }

res = pd.concat([prim, sens], ignore_index=True)
res.to_csv(f"{HERE}/results.csv", index=False)
json.dump(OUT, open(f"{HERE}/summary.json", "w"), indent=1, default=str)

# ======================================================================================
# 8. READOUT
# ======================================================================================
log("\n" + "=" * 95)
for k in NONREF:
    s = OUT["cells_vs_reference"][k]
    log(f"\n{k}  {s['label']}   n={counts[k]:,}")
    log(f"  {s['n_significant']} of {s['n_conditions']} significant "
        f"({s['n_expected_by_chance']} expected by chance), {s['n_significant_fdr']} after FDR, "
        f"median HR {s['median_hr']}, range {s['hr_range']}")
    log(f"  confounding floor {s['confounding_floor']}, {s['n_clearing_floor']} findings clear it")
    for x in s["significant"][:10]:
        log(f"    {x['disease'][:34]:34s} HR {x['hr']:.2f} [{x['ci'][0]:.2f}-{x['ci'][1]:.2f}] "
            f"p={x['p']:.2g} q={x['q']} floor={x['clears_floor']}")

log("\n" + "=" * 95)
log("OXYGEN EFFECT WITHIN EACH DURATION STRATUM")
for c in OX_CONTRASTS:
    s = OUT["oxygen_within_duration"][c[0]]
    log(f"\n{c[0]}  {s['label']}")
    log(f"  {s['n_significant']} of {s['n_conditions']} significant, "
        f"{s['n_significant_fdr']} after FDR, median HR {s['median_hr']}, range {s['hr_range']}, "
        f"floor {s['confounding_floor']}")
    for x in s["significant"][:10]:
        log(f"    {x['disease'][:34]:34s} HR {x['hr']:.2f} [{x['ci'][0]:.2f}-{x['ci'][1]:.2f}] "
            f"p={x['p']:.2g} q={x['q']} floor={x['clears_floor']}")

log("\nminimum detectable HR (80% power, alpha 0.05)")
for k in NONREF:
    m = OUT["minimum_detectable_hr"][k]
    log(f"  {k}: median {m['median_across_conditions']}, death {m['best_case_death']}, "
        f"worst {m['worst_case']}, conditions with <5 events {m['n_conditions_with_lt5_events']}, "
        f"with 0 events {m['n_conditions_with_0_events']}")

log("\nsensitivity, recordings of at least 6 hours")
for k in NONREF + [c[0] for c in OX_CONTRASTS]:
    a = OUT["sensitivity_agreement"][k]
    log(f"  {k}: significant {a['sig_primary']} -> {a['sig_sensitivity']}, both {a['sig_both']}, "
        f"median |change| {a['median_abs_pct_change_hr']}%, rank r {a['spearman_hr']}")

log("\nwrote cells.csv, results.csv, interaction.csv, summary.json")
