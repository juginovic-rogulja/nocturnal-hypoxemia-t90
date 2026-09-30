"""
Sleep duration crossed with nocturnal oxygenation, restricted to people who were free of organ
disease at the sleep study.

This repeats dur_x_oxygen/model.py inside a disease-free subgroup. The instruction that defines
the subgroup is that a short or a long night is NOT itself a disease, so sleep duration is never
an exclusion. Neither are the sleep disorders: insomnia, restless legs, nocturia and sleep apnoea
stay in, because excluding them would remove the exposure. Only ORGAN disease is excluded.

Healthy definition, built from ORGAN_GROUP in disease_definitions.py: no prevalent condition in
the Cardiac, Respiratory, Kidney, Liver or Metabolic groups, plus no prevalent stroke (grouped
Neuro/psych but a vascular disease, and excluded by both of the paper's earlier healthy
definitions) and no prevalent cancer (cancer_any, prostate_ca, grouped Other). Neuro/psych,
Infection and Sensory/MSK conditions are NOT exclusions, which keeps the five negative controls
available inside the subgroup.

Cells:
    cell1  short sleep  + normal oxygen    TST < SHORT, T90 <= 1%
    cell2  short sleep  + low oxygen       TST < SHORT, T90 > 10%
    cell3  normal sleep + low oxygen       TST >= 360 min, T90 > 10%
    cell4  normal sleep + normal oxygen    TST >= 360 min, T90 <= 1%   REFERENCE

SHORT is run at both 240 min (matching the full-cohort run exactly) and 300 min (more power once
organ disease is gone). Both are reported; the primary is chosen on cell occupancy, not on result.

Because cell 2 is thin in a disease-free population, every condition is fitted twice:
  - the four-cell factor, which needs events in all four cells and is directly comparable to the
    full-cohort run, and
  - pairwise two-cell models (1v4, 2v4, 3v4, 2v1), which need events only in the two cells being
    compared and therefore keep conditions the four-cell fit has to drop.

Modelling: site-stratified Cox, unpenalized, natural cubic age spline with 4 df of which ONE
column is dropped so the design is full rank, plus sex. PSG date is time zero. The spline
coefficients are checked to be non-zero, because an unpenalized fit on a rank-deficient design
returns zeros without failing.

TST here is one laboratory night, not habitual sleep duration, so every primary model is repeated
on recordings of at least 360 minutes.

Writes cells.csv, results.csv, interaction.csv, comparison.csv and summary.json in this directory.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import json
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
HERE = f"{paths.SV_ROOT}/dur_x_oxygen_healthy"
FULLRUN = f"{paths.SV_ROOT}/dur_x_oxygen"
sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
from cohort_spec import apply_cohort                          # noqa: E402
from disease_definitions import (DISEASES, NEGATIVE_CONTROLS,  # noqa: E402
                                 CIRCULAR, ORGAN_GROUP)

T90_LOW, T90_NORMAL = 10.0, 1.0        # low oxygen = T90 > 10%, normal oxygen = T90 <= 1%
TST_NORMAL = 360.0                     # normal sleep >= 6 h
TST_SHORT_MATCHED = 240.0              # matches the full-cohort run
TST_SHORT_POWERED = 300.0              # planned secondary, 5 h
MIN_EVENTS = 20                        # per condition, summed over the cells in the model
MIN_CELL_EVENTS_PAIR = 3               # a pairwise fit needs this many events in each cell
REC_MIN = 360.0                        # sensitivity: recordings of at least 6 h

OUT = {}
CELLS = ["cell1", "cell2", "cell3", "cell4"]
CELLD = ["cell1", "cell2", "cell3"]
CELL_LABEL = {
    "cell1": "short sleep, normal oxygen",
    "cell2": "short sleep, low oxygen",
    "cell3": "normal sleep, low oxygen",
    "cell4": "normal sleep, normal oxygen (reference)",
}


def log(m):
    print(m, flush=True)


# ======================================================================================
# 1. COHORT AND THE HEALTHY DEFINITION
# ======================================================================================
b = apply_cohort(pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet"))
log(f"analysis cohort {len(b):,}")

b["male"] = (b.sex.astype(str).str.upper().str[0] == "M").astype(float)
b.loc[~b.sex.astype(str).str.upper().str[0].isin(["M", "F"]), "male"] = np.nan

ORGAN_BLOCKS = [
    ("cardiovascular", [k for k, g in ORGAN_GROUP.items() if g == "Cardiac"] +
     ["stroke_any", "cvd"]),
    ("pulmonary", [k for k, g in ORGAN_GROUP.items() if g == "Respiratory"]),
    ("renal", [k for k, g in ORGAN_GROUP.items() if g == "Kidney"]),
    ("hepatic", [k for k, g in ORGAN_GROUP.items() if g == "Liver"]),
    ("metabolic", [k for k, g in ORGAN_GROUP.items() if g == "Metabolic"]),
    ("cancer", ["cancer_any", "prostate_ca"]),
]
EXCLUDE_KEYS = [k for _, ks in ORGAN_BLOCKS for k in ks]
assert len(EXCLUDE_KEYS) == len(set(EXCLUDE_KEYS))

# conditions deliberately NOT used as exclusions, listed so the definition is auditable
NOT_EXCLUDED = [k for k in DISEASES if k not in EXCLUDE_KEYS]

ladder, keep = [], np.ones(len(b), bool)
for block, keys in ORGAN_BLOCKS:
    for k in keys:
        prev = (b[f"{k}_prevalent"] == 1).values
        before = int(keep.sum())
        keep = keep & ~prev
        ladder.append({"block": block, "key": k, "condition": DISEASES[k][0],
                       "prevalent_in_cohort": int(prev.sum()),
                       "removed_here": before - int(keep.sum()),
                       "remaining": int(keep.sum())})
h = b[keep].copy()
N_HEALTHY = int(len(h))
log(f"healthy (organ-disease-free) {N_HEALTHY:,}  ({100 * N_HEALTHY / len(b):.1f}% of cohort)")

# the two healthy definitions the paper has used before, recomputed on THIS cohort so the
# comparison is like for like rather than against a number from a different filter
PAPER15 = ["hf", "cvd", "ihd", "mi", "afib", "stroke_any", "htn2", "diabetes", "ckd", "copd2",
           "asthma", "cancer_any", "dementia", "cirrhosis", "obesity_hypovent"]
WS_B = ["copd2", "asthma", "resp_failure", "obesity_hypovent", "pulm_htn", "hf", "ihd", "mi",
        "afib", "stroke_any", "ckd", "cirrhosis", "cancer_any"]
WS_C = WS_B + ["diabetes", "obesity", "nafld"]


def n_free(keys):
    m = np.ones(len(b), bool)
    for k in keys:
        m &= (b[f"{k}_prevalent"] == 0).values
    return int(m.sum())


OUT["healthy_definition"] = {
    "rule": "no prevalent disease in ORGAN_GROUP Cardiac, Respiratory, Kidney, Liver or "
            "Metabolic, plus no prevalent stroke and no prevalent cancer",
    "n_conditions_excluded_on": len(EXCLUDE_KEYS),
    "excluded_on": {blk: [DISEASES[k][0] for k in ks] for blk, ks in ORGAN_BLOCKS},
    "deliberately_not_excluded": [DISEASES[k][0] for k in NOT_EXCLUDED],
    "sleep_duration_is_not_an_exclusion": True,
    "sleep_disorders_not_excluded": ["Insomnia", "Restless legs syndrome", "Nocturia",
                                     "Sleep apnoea (AHI is never an exclusion)"],
    "n_cohort": int(len(b)),
    "n_healthy": N_HEALTHY,
    "pct_of_cohort": round(100 * N_HEALTHY / len(b), 1),
    "exclusion_ladder": ladder,
    "comparison_with_earlier_paper_definitions": {
        "this_run": N_HEALTHY,
        "paper_15_condition_rule_on_this_cohort": n_free(PAPER15),
        "paper_15_condition_rule_as_published": 6550,
        "wake_sleep_stratum_B_on_this_cohort": n_free(WS_B),
        "wake_sleep_stratum_C_on_this_cohort": n_free(WS_C),
        "wake_sleep_stratum_C_as_published": 1981,
    },
    "prevalent_dementia_remaining_in_healthy": int((h.dementia_prevalent == 1).sum()),
}
log("exclusion ladder")
for r in ladder:
    log(f"  {r['block']:16s} {r['condition']:28s} prevalent {r['prevalent_in_cohort']:6,d}  "
        f"removed {r['removed_here']:6,d}  remaining {r['remaining']:6,d}")
log(f"earlier definitions on this same cohort: 15-condition rule "
    f"{n_free(PAPER15):,} (published 6,550), wake/sleep stratum C {n_free(WS_C):,} "
    f"(published 1,981 on the 6,803 wake/sleep set)")


# ======================================================================================
# 2. CELLS AT BOTH SHORT-SLEEP THRESHOLDS
# ======================================================================================
def assign_cells(d, tst_short):
    lowO2 = d.spo2_pct_below_90 > T90_LOW
    normO2 = d.spo2_pct_below_90 <= T90_NORMAL
    shortS = d.TST_min < tst_short
    normS = d.TST_min >= TST_NORMAL
    return pd.Series(np.select(
        [shortS & normO2, shortS & lowO2, normS & lowO2, normS & normO2],
        CELLS, default="excluded"), index=d.index)


def discards(d, tst_short, cell):
    return {
        "n_in": int(len(d)),
        "t90_middle_1_to_10_pct": int(((d.spo2_pct_below_90 > T90_NORMAL) &
                                       (d.spo2_pct_below_90 <= T90_LOW)).sum()),
        "tst_middle_band_min": int(((d.TST_min >= tst_short) &
                                    (d.TST_min < TST_NORMAL)).sum()),
        "tst_missing": int(d.TST_min.isna().sum()),
        "in_four_cells": int((cell != "excluded").sum()),
        "excluded_total": int((cell == "excluded").sum()),
    }


cellcount = {}
for tag, thr in (("tst240", TST_SHORT_MATCHED), ("tst300", TST_SHORT_POWERED)):
    for pop, d in (("healthy", h), ("full_cohort", b)):
        c = assign_cells(d, thr)
        cellcount[f"{pop}_{tag}"] = {
            "tst_short_lt_min": thr,
            "n_by_cell": {k: int((c == k).sum()) for k in CELLS},
            "deaths_by_cell": {k: int(d.loc[c == k, "death_incident"].sum()) for k in CELLS},
            "discards": discards(d, thr, c),
        }
OUT["cell_counts_side_by_side"] = cellcount

log("\ncell sizes side by side (n / incident deaths)")
log(f"{'':26s}{'cell1':>18s}{'cell2':>18s}{'cell3':>18s}{'cell4':>18s}")
for k, v in cellcount.items():
    row = "".join(f"{v['n_by_cell'][c]:>13,d} /{v['deaths_by_cell'][c]:>4d}" for c in CELLS)
    log(f"{k:26s}{row}")

# Primary threshold: chosen on cell occupancy. Cell 2 is the cell every question depends on.
n_c2_240 = cellcount["healthy_tst240"]["n_by_cell"]["cell2"]
n_c2_300 = cellcount["healthy_tst300"]["n_by_cell"]["cell2"]
PRIMARY_THR = TST_SHORT_POWERED if n_c2_240 < 100 else TST_SHORT_MATCHED
OUT["primary_threshold"] = {
    "tst_short_lt_min": PRIMARY_THR,
    "rule": "the 4-hour cut is primary unless it leaves fewer than 100 people in cell 2, the "
            "cell every question depends on, in which case the planned 5-hour cut is primary",
    "cell2_n_at_240": n_c2_240, "cell2_n_at_300": n_c2_300,
    "chosen": "300 min (5 h)" if PRIMARY_THR == TST_SHORT_POWERED else "240 min (4 h)",
}
log(f"\nprimary short-sleep threshold: TST < {PRIMARY_THR:.0f} min "
    f"(cell 2 holds {n_c2_240} at 4 h and {n_c2_300} at 5 h)")

h240 = h.assign(cell=assign_cells(h, TST_SHORT_MATCHED))
h300 = h.assign(cell=assign_cells(h, TST_SHORT_POWERED))
f240 = h240[h240.cell != "excluded"].copy()
f300 = h300[h300.cell != "excluded"].copy()
FPRIM = f300 if PRIMARY_THR == TST_SHORT_POWERED else f240

# ======================================================================================
# 3. CELL PROFILE  ->  cells.csv
# ======================================================================================
PREV_CHECK = ["insomnia_prevalent", "rls_plmd_prevalent", "nocturia_prevalent",
              "depression_prevalent", "dementia_prevalent", "alopecia_prevalent"]   # v8.1: alopecia took back pain's slot
prof = []
for tag, f4 in (("tst240", f240), ("tst300", f300)):
    for k in CELLS:
        s = f4[f4.cell == k]
        row = {
            "threshold": tag, "tst_short_lt_min": (TST_SHORT_MATCHED if tag == "tst240"
                                                   else TST_SHORT_POWERED),
            "is_primary": bool((tag == "tst300") == (PRIMARY_THR == TST_SHORT_POWERED)),
            "cell": k, "definition": CELL_LABEL[k], "n": int(len(s)),
            "age_median": round(float(s.AgeAtVisit.median()), 1),
            "age_q1": round(float(s.AgeAtVisit.quantile(.25)), 1),
            "age_q3": round(float(s.AgeAtVisit.quantile(.75)), 1),
            "pct_male": round(100 * float(s.male.mean()), 1),
            "t90_median": round(float(s.spo2_pct_below_90.median()), 3),
            "t90_q1": round(float(s.spo2_pct_below_90.quantile(.25)), 3),
            "t90_q3": round(float(s.spo2_pct_below_90.quantile(.75)), 3),
            "tst_min_median": round(float(s.TST_min.median()), 1),
            "tst_min_q1": round(float(s.TST_min.quantile(.25)), 1),
            "tst_min_q3": round(float(s.TST_min.quantile(.75)), 1),
            "recording_dur_min_median": round(float(s.recording_dur_min.median()), 1),
            "recording_dur_min_q1": round(float(s.recording_dur_min.quantile(.25)), 1),
            "recording_dur_min_q3": round(float(s.recording_dur_min.quantile(.75)), 1),
            "pct_recording_lt360min": round(100 * float((s.recording_dur_min < REC_MIN).mean()), 1),
            "n_recording_ge360": int((s.recording_dur_min >= REC_MIN).sum()),
            "sleep_eff_pct_median": round(float(s.sleep_efficiency_pct.median()), 1),
            "sleep_eff_pct_q1": round(float(s.sleep_efficiency_pct.quantile(.25)), 1),
            "sleep_eff_pct_q3": round(float(s.sleep_efficiency_pct.quantile(.75)), 1),
            "ahi_median": round(float(s.AHI.median()), 2),
            "ahi_q1": round(float(s.AHI.quantile(.25)), 2),
            "ahi_q3": round(float(s.AHI.quantile(.75)), 2),
            "pct_ahi_ge15": round(100 * float((s.AHI >= 15).mean()), 1),
            "pct_ahi_ge30": round(100 * float((s.AHI >= 30).mean()), 1),
            "pct_ahi_lt5": round(100 * float((s.AHI < 5).mean()), 1),
            "arousal_index_median": round(float(s.arousal_index.median()), 2),
            "spo2_mean_median": round(float(s.spo2_mean.median()), 2),
            "odi3_median": round(float(s.odi3_total.median()), 2),
            "n_deaths": int(s.death_incident.sum()),
            "death_rate_pct": round(100 * float(s.death_incident.mean()), 2),
            "followup_years_median": round(float(s.death_years.median()), 2),
        }
        for c in PREV_CHECK:
            row[c.replace("_prevalent", "_prev_pct")] = round(100 * float(s[c].mean()), 1)
        for site, cnt in s.site_id.value_counts().items():
            row[f"site_{site}"] = int(cnt)
        prof.append(row)
cells = pd.DataFrame(prof)
cells.to_csv(f"{HERE}/cells.csv", index=False)
log("\ncell profile")
log(cells[["threshold", "cell", "n", "age_median", "pct_male", "tst_min_median",
           "recording_dur_min_median", "sleep_eff_pct_median", "t90_median",
           "ahi_median", "pct_ahi_ge15", "n_deaths"]].to_string(index=False))

OUT["cell_profile_sanity"] = {
    tag: {k: {"n": int(r.n), "ahi_median": float(r.ahi_median),
              "ahi_iqr": [float(r.ahi_q1), float(r.ahi_q3)],
              "pct_ahi_ge15": float(r.pct_ahi_ge15), "pct_ahi_lt5": float(r.pct_ahi_lt5),
              "recording_dur_min_median": float(r.recording_dur_min_median),
              "sleep_eff_pct_median": float(r.sleep_eff_pct_median),
              "pct_recording_lt360min": float(r.pct_recording_lt360min)}
          for k, r in cells[cells.threshold == tag].set_index("cell").iterrows()}
    for tag in ("tst240", "tst300")}


# ======================================================================================
# 4. MODEL MACHINERY
# ======================================================================================
SPLINE_COEFS = []


def agespline(d):
    """4-df natural cubic spline on age, first column dropped so the design is full rank.

    The four basis columns sum to one. A Cox partial likelihood has no intercept, so under a
    penalty the redundancy is harmless, but an unpenalized fit on a rank-deficient design can
    return coefficients of exactly zero without raising. One column is always dropped."""
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
    for k2 in CELLD:
        out[k2] = (out.cell_lab == k2).astype(float)
    return out


def _fit(dd, terms):
    cov = terms + ["male"] + [c for c in dd.columns if c.startswith("age_s")]
    m = CoxPHFitter(penalizer=0.0).fit(dd[cov + ["T", "E", "site"]], "T", "E", strata=["site"])
    SPLINE_COEFS.append(float(np.abs(
        m.params_[[c for c in m.params_.index if c.startswith("age_s")]]).max()))
    return m


def contrast(m, a, ref):
    """HR for cell a against cell ref, both non-baseline, from the fitted covariance."""
    V, beta = m.variance_matrix_, m.params_
    d = float(beta[a] - beta[ref])
    se = float(np.sqrt(V.loc[a, a] + V.loc[ref, ref] - 2 * V.loc[a, ref]))
    z = d / se
    return (round(float(np.exp(d)), 3), round(float(np.exp(d - 1.96 * se)), 3),
            round(float(np.exp(d + 1.96 * se)), 3), float(2 * stats.norm.sf(abs(z))))


PAIRS = [("cell1", "cell4"), ("cell2", "cell4"), ("cell3", "cell4"), ("cell2", "cell1")]


def pairfit(dd, a, ref):
    """Two-cell Cox: only the two cells being compared enter, so a condition with no events in
    a third cell is still estimable. Returns None if either cell is too thin."""
    sub = dd[dd.cell_lab.isin([a, ref])].copy()
    ea, er = int(sub[sub.cell_lab == a].E.sum()), int(sub[sub.cell_lab == ref].E.sum())
    if ea + er < MIN_EVENTS or min(ea, er) < MIN_CELL_EVENTS_PAIR:
        return None
    sub["X"] = (sub.cell_lab == a).astype(float)
    try:
        m = _fit(sub, ["X"])
    except Exception:
        return None
    s = m.summary.loc["X"]
    return {"hr": round(float(s["exp(coef)"]), 3),
            "lo": round(float(s["exp(coef) lower 95%"]), 3),
            "hi": round(float(s["exp(coef) upper 95%"]), 3),
            "p": float(s["p"]), "ev_a": ea, "ev_ref": er,
            "n_a": int((sub.cell_lab == a).sum()), "n_ref": int((sub.cell_lab == ref).sum()),
            "unstable": bool(min(ea, er) < 10)}


def run_panel(d, tag):
    """Every non-circular condition, four-cell model where possible, pairwise models always."""
    rows, skipped = [], []
    for key, label in OUTCOMES:
        dd = frame_for(d, key)
        if dd is None:
            skipped.append({"key": key, "disease": DISEASES.get(key, (key,))[0],
                            "reason": "columns absent"})
            continue
        ev = {k2: int(dd[dd.cell_lab == k2].E.sum()) for k2 in CELLS}
        nn = {k2: int((dd.cell_lab == k2).sum()) for k2 in CELLS}
        rec = {"analysis": tag, "key": key, "disease": label,
               "organ_group": ORGAN_GROUP.get(key, "Death"),
               "negative_control": key in NEGATIVE_CONTROLS,
               "used_as_healthy_exclusion": key in EXCLUDE_KEYS,
               "n_total": int(len(dd)), "events_total": int(dd.E.sum()),
               **{f"n_{k2}": nn[k2] for k2 in CELLS},
               **{f"ev_{k2}": ev[k2] for k2 in CELLS},
               "cells_with_zero_events": ",".join(k2 for k2 in CELLS if ev[k2] == 0),
               "unstable_lt5_events": bool(min(ev.values()) < 5)}

        # ---- pairwise two-cell models, always attempted
        any_pair = False
        for a, ref in PAIRS:
            r = pairfit(dd, a, ref)
            pre = f"p_{a}_vs_{ref}"
            if r is None:
                rec[f"{pre}_hr"] = np.nan
                continue
            any_pair = True
            rec.update({f"{pre}_hr": r["hr"], f"{pre}_lo": r["lo"], f"{pre}_hi": r["hi"],
                        f"{pre}_p": r["p"], f"{pre}_ev_a": r["ev_a"],
                        f"{pre}_ev_ref": r["ev_ref"], f"{pre}_unstable": r["unstable"]})

        # ---- four-cell model, comparable to the full-cohort run
        rec["fourcell_fitted"] = False
        if dd.E.sum() >= MIN_EVENTS and ev["cell4"] > 0 and min(ev[k2] for k2 in CELLD) > 0:
            try:
                m = _fit(dd, CELLD)
                rec["fourcell_fitted"] = True
                for k2 in CELLD:
                    s = m.summary.loc[k2]
                    rec[f"{k2}_hr"] = round(float(s["exp(coef)"]), 3)
                    rec[f"{k2}_lo"] = round(float(s["exp(coef) lower 95%"]), 3)
                    rec[f"{k2}_hi"] = round(float(s["exp(coef) upper 95%"]), 3)
                    rec[f"{k2}_p"] = float(s["p"])
                for a, ref in (("cell2", "cell1"), ("cell2", "cell3")):
                    hr, lo, hi, p = contrast(m, a, ref)
                    rec.update({f"{a}_vs_{ref}_hr": hr, f"{a}_vs_{ref}_lo": lo,
                                f"{a}_vs_{ref}_hi": hi, f"{a}_vs_{ref}_p": p})
                rec["product_c1_c3"] = round(rec["cell1_hr"] * rec["cell3_hr"], 3)
                rec["cell2_over_product"] = round(rec["cell2_hr"] / rec["product_c1_c3"], 3)
            except Exception as ex:
                rec["fourcell_reason"] = f"fit failed: {type(ex).__name__}"
        else:
            rec["fourcell_reason"] = (
                f"{int(dd.E.sum())} events, below {MIN_EVENTS}" if dd.E.sum() < MIN_EVENTS
                else "a cell has no events: " + rec["cells_with_zero_events"])

        if not rec["fourcell_fitted"] and not any_pair:
            skipped.append({"key": key, "disease": label,
                            "reason": rec.get("fourcell_reason", "no model estimable"),
                            "events_total": int(dd.E.sum()), "events_by_cell": ev})
            continue
        rows.append(rec)

    df = pd.DataFrame(rows)
    if len(df):
        real = ~df.negative_control
        for k2 in CELLD:
            if f"{k2}_p" in df.columns:
                df[f"{k2}_q"] = np.nan
                ok = real & df[f"{k2}_p"].notna()
                if ok.sum():
                    df.loc[ok, f"{k2}_q"] = multipletests(df.loc[ok, f"{k2}_p"],
                                                          method="fdr_bh")[1]
        for a, ref in PAIRS:
            pre = f"p_{a}_vs_{ref}"
            df[f"{pre}_q"] = np.nan
            ok = real & df[f"{pre}_hr"].notna()
            if ok.sum():
                df.loc[ok, f"{pre}_q"] = multipletests(df.loc[ok, f"{pre}_p"],
                                                       method="fdr_bh")[1]
    log(f"\n[{tag}] {len(df)} conditions with at least one estimable model, "
        f"{int(df.fourcell_fitted.sum()) if len(df) else 0} with the full four-cell fit, "
        f"{len(skipped)} with nothing estimable")
    return df, skipped


def floor_of(df, col_hr, col_lo, col_hi):
    """Confounding floor: the largest negative-control deviation from 1."""
    neg = df[df.negative_control & df[col_hr].notna()]
    if not len(neg):
        return None
    dev = neg[col_hr].apply(lambda x: max(x, 1 / x))
    i = dev.idxmax()
    return {"max_abs_hr": round(float(dev.max()), 3),
            "at_condition": str(neg.loc[i, "disease"]),
            "hr_at": float(neg.loc[i, col_hr]),
            "ci_at": [float(neg.loc[i, col_lo]), float(neg.loc[i, col_hi])],
            "range": [round(float(neg[col_hr].min()), 3), round(float(neg[col_hr].max()), 3)],
            "controls": {str(r.disease): [float(r[col_hr]), float(r[col_lo]), float(r[col_hi])]
                         for _, r in neg.iterrows()},
            "n_controls": int(len(neg)),
            "n_controls_significant": int(((neg[col_lo] > 1) | (neg[col_hi] < 1)).sum())}


def sig_block(df, hr, lo, hi, p, q, floor):
    r = df[(~df.negative_control) & df[hr].notna()]
    s = r[(r[lo] > 1) | (r[hi] < 1)]
    return {
        "n_conditions_tested": int(len(r)),
        "n_significant": int(len(s)),
        "n_expected_by_chance": round(0.05 * len(r), 1),
        "n_significant_fdr": int((r[q] < 0.05).sum()),
        "median_hr": round(float(r[hr].median()), 3),
        "hr_range": [round(float(r[hr].min()), 3), round(float(r[hr].max()), 3)],
        "n_clearing_negative_control_floor": (
            int(s[hr].apply(lambda x: max(x, 1 / x) > floor["max_abs_hr"]).sum())
            if floor else None),
        "significant": [{"disease": str(x.disease), "hr": float(x[hr]),
                         "ci": [float(x[lo]), float(x[hi])], "p": float(x[p]),
                         "q": round(float(x[q]), 4),
                         "clears_floor": (bool(max(x[hr], 1 / x[hr]) > floor["max_abs_hr"])
                                          if floor else None)}
                        for _, x in s.sort_values(hr, ascending=False).iterrows()],
    }


# ======================================================================================
# 5. PANELS
# ======================================================================================
panels, skips = {}, {}
for tag, f4 in (("healthy_tst240", f240), ("healthy_tst300", f300)):
    panels[tag], skips[tag] = run_panel(f4, tag)

PRIMTAG = "healthy_tst300" if PRIMARY_THR == TST_SHORT_POWERED else "healthy_tst240"
MATCHTAG = "healthy_tst240"
prim = panels[PRIMTAG]

# sensitivity: recordings of at least 6 hours, on the primary threshold
FPRIM_REC = FPRIM[FPRIM.recording_dur_min >= REC_MIN].copy()
OUT["long_recording_subset"] = {
    "rule": f"recording_dur_min >= {REC_MIN:.0f}",
    "threshold": PRIMTAG,
    "n": int(len(FPRIM_REC)),
    "by_cell": {k: int((FPRIM_REC.cell == k).sum()) for k in CELLS},
    "pct_retained_by_cell": {k: round(100 * float((FPRIM_REC.cell == k).sum() /
                                                  max((FPRIM.cell == k).sum(), 1)), 1)
                             for k in CELLS},
}
log("\nlong-recording subset " + json.dumps(OUT["long_recording_subset"]["by_cell"]))
panels["rec_ge360"], skips["rec_ge360"] = run_panel(FPRIM_REC, "rec_ge360")

# Age-overlap sensitivity. Once organ disease is excluded the reference cell is much younger
# than the low-oxygen cells (medians 35 against 65), so the age spline carries a heavy load and
# the negative controls in cell 1 all sit below 1. This repeats the primary panel inside an age
# band where all four cells are well populated, so the contrasts do not rest on extrapolation.
AGE_LO, AGE_HI = 35.0, 75.0
FPRIM_AGE = FPRIM[(FPRIM.AgeAtVisit >= AGE_LO) & (FPRIM.AgeAtVisit <= AGE_HI)].copy()
OUT["age_overlap_subset"] = {
    "rule": f"{AGE_LO:.0f} <= age <= {AGE_HI:.0f}",
    "reason": "the disease-free reference cell is much younger than the low-oxygen cells, so "
              "the age spline is doing heavy work and the cell-1 negative controls sit below 1",
    "age_median_by_cell_unrestricted": {k: float(FPRIM.loc[FPRIM.cell == k,
                                                           "AgeAtVisit"].median())
                                        for k in CELLS},
    "age_iqr_by_cell_unrestricted": {k: [float(FPRIM.loc[FPRIM.cell == k,
                                                         "AgeAtVisit"].quantile(.25)),
                                         float(FPRIM.loc[FPRIM.cell == k,
                                                         "AgeAtVisit"].quantile(.75))]
                                     for k in CELLS},
    "n": int(len(FPRIM_AGE)),
    "by_cell": {k: int((FPRIM_AGE.cell == k).sum()) for k in CELLS},
    "age_median_by_cell_restricted": {k: float(FPRIM_AGE.loc[FPRIM_AGE.cell == k,
                                                             "AgeAtVisit"].median())
                                      for k in CELLS},
}
log("age-overlap subset " + json.dumps(OUT["age_overlap_subset"]["by_cell"]))
panels["age_35_75"], skips["age_35_75"] = run_panel(FPRIM_AGE, "age_35_75")

OUT["skipped"] = skips
OUT["spline_check"] = {
    "n_fits": len(SPLINE_COEFS),
    "min_abs_max_age_coef": round(float(np.min(SPLINE_COEFS)), 6),
    "median_abs_max_age_coef": round(float(np.median(SPLINE_COEFS)), 4),
    "n_fits_with_all_zero_age_coefs": int(sum(1 for x in SPLINE_COEFS if x < 1e-8)),
    "verdict": ("age spline fits real coefficients in every model"
                if min(SPLINE_COEFS) > 1e-8 else "SOME FITS RETURNED ZERO AGE COEFFICIENTS"),
}
assert OUT["spline_check"]["n_fits_with_all_zero_age_coefs"] == 0, OUT["spline_check"]
log(f"\nage-spline check: {OUT['spline_check']['n_fits']} fits, smallest |max age coef| "
    f"{OUT['spline_check']['min_abs_max_age_coef']}, none zero")

# ======================================================================================
# 6. FLOORS AND HEADLINE SUMMARIES
# ======================================================================================
OUT["negative_control_floor"] = {}
OUT["summary"] = {}
for tag, df in panels.items():
    fl, sm = {}, {}
    for a, ref in PAIRS:
        pre = f"p_{a}_vs_{ref}"
        f = floor_of(df, f"{pre}_hr", f"{pre}_lo", f"{pre}_hi")
        fl[f"{a}_vs_{ref}"] = f
        sm[f"{a}_vs_{ref}"] = sig_block(df, f"{pre}_hr", f"{pre}_lo", f"{pre}_hi",
                                        f"{pre}_p", f"{pre}_q", f)
    for k2 in CELLD:
        if f"{k2}_hr" in df.columns and df[f"{k2}_hr"].notna().any():
            f = floor_of(df, f"{k2}_hr", f"{k2}_lo", f"{k2}_hi")
            fl[f"fourcell_{k2}"] = f
            sm[f"fourcell_{k2}"] = sig_block(df, f"{k2}_hr", f"{k2}_lo", f"{k2}_hi",
                                             f"{k2}_p", f"{k2}_q", f)
    OUT["negative_control_floor"][tag] = fl
    OUT["summary"][tag] = sm

# ======================================================================================
# 7. INTERACTION
# ======================================================================================
inter_rows = []
cand = prim[prim.fourcell_fitted].sort_values("events_total", ascending=False)
for _, r in cand.iterrows():
    dd = frame_for(FPRIM, r.key)
    dd["short"] = dd.cell_lab.isin(["cell1", "cell2"]).astype(float)
    dd["lowo2"] = dd.cell_lab.isin(["cell2", "cell3"]).astype(float)
    dd["short_x_lowo2"] = dd["short"] * dd["lowo2"]
    try:
        m = _fit(dd, ["short", "lowo2", "short_x_lowo2"])
    except Exception:
        continue
    s = m.summary
    inter_rows.append({
        "analysis": PRIMTAG, "key": r.key, "disease": r.disease,
        "events": int(dd.E.sum()), "ev_cell2": int(r.ev_cell2),
        "main_short_hr": round(float(s.loc["short", "exp(coef)"]), 3),
        "main_lowo2_hr": round(float(s.loc["lowo2", "exp(coef)"]), 3),
        "interaction_hr": round(float(s.loc["short_x_lowo2", "exp(coef)"]), 3),
        "interaction_lo": round(float(s.loc["short_x_lowo2", "exp(coef) lower 95%"]), 3),
        "interaction_hi": round(float(s.loc["short_x_lowo2", "exp(coef) upper 95%"]), 3),
        "interaction_p": float(s.loc["short_x_lowo2", "p"]),
        "cell2_hr": float(r.cell2_hr), "product_c1_c3": float(r.product_c1_c3),
        "cell2_exceeds_product": bool(r.cell2_hr > r.product_c1_c3),
    })
inter = pd.DataFrame(inter_rows)
if len(inter):
    inter["interaction_q"] = multipletests(inter.interaction_p, method="fdr_bh")[1]
inter.to_csv(f"{HERE}/interaction.csv", index=False)
OUT["interaction"] = {
    "analysis": PRIMTAG,
    "n_tested": int(len(inter)),
    "n_significant": int(((inter.interaction_lo > 1) | (inter.interaction_hi < 1)).sum())
    if len(inter) else 0,
    "n_significant_fdr": int((inter.interaction_q < 0.05).sum()) if len(inter) else 0,
    "p_range": [round(float(inter.interaction_p.min()), 4),
                round(float(inter.interaction_p.max()), 4)] if len(inter) else None,
    "median_interaction_hr": round(float(inter.interaction_hr.median()), 3) if len(inter) else None,
    "rows": json.loads(inter.to_json(orient="records")) if len(inter) else [],
}
if len(inter):
    log("\ninteraction, short by low oxygen, primary threshold")
    log(inter[["disease", "events", "interaction_hr", "interaction_lo", "interaction_hi",
               "interaction_p"]].to_string(index=False))

# ======================================================================================
# 8. HEALTHY AGAINST THE FULL COHORT  ->  comparison.csv
# ======================================================================================
full = pd.read_csv(f"{FULLRUN}/results.csv")
full = full[full.analysis == "primary"].set_index("key")

BIG, SMALL = 1.10, 1 / 1.10


def verdict(hh, ff):
    if pd.isna(hh) or pd.isna(ff) or ff == 0:
        return None, None
    ratio = hh / ff
    return round(float(ratio), 3), ("LARGER" if ratio > BIG else
                                    "SMALLER" if ratio < SMALL else "UNCHANGED")


comp_rows = []
for tag in ("healthy_tst240", "healthy_tst300"):
    df = panels[tag]
    for _, r in df.iterrows():
        if r.key not in full.index:
            continue
        fr = full.loc[r.key]
        base = {"healthy_analysis": tag, "matched_to_full_cohort": tag == MATCHTAG,
                "key": r.key, "disease": r.disease,
                "negative_control": bool(r.negative_control),
                "healthy_events": int(r.events_total),
                "full_cohort_events": int(fr.events_total)}
        for name, hcol, fcol in (
                ("cell1_vs_cell4", "p_cell1_vs_cell4", "cell1"),
                ("cell2_vs_cell4", "p_cell2_vs_cell4", "cell2"),
                ("cell3_vs_cell4", "p_cell3_vs_cell4", "cell3"),
                ("cell2_vs_cell1", "p_cell2_vs_cell1", "cell2_vs_cell1")):
            hh = r.get(f"{hcol}_hr", np.nan)
            ff = float(fr[f"{fcol}_hr"])
            ratio, v = verdict(hh, ff)
            row = dict(base)
            row.update({"contrast": name, "healthy_hr": (None if pd.isna(hh) else float(hh)),
                        "healthy_lo": (None if pd.isna(hh) else float(r[f"{hcol}_lo"])),
                        "healthy_hi": (None if pd.isna(hh) else float(r[f"{hcol}_hi"])),
                        "healthy_p": (None if pd.isna(hh) else float(r[f"{hcol}_p"])),
                        "full_cohort_hr": ff,
                        "full_cohort_lo": float(fr[f"{fcol}_lo"]),
                        "full_cohort_hi": float(fr[f"{fcol}_hi"]),
                        "ratio_healthy_over_full": ratio, "verdict": v,
                        "full_hr_inside_healthy_ci": (
                            None if pd.isna(hh) else
                            bool(float(r[f"{hcol}_lo"]) <= ff <= float(r[f"{hcol}_hi"])))})
            comp_rows.append(row)
comp = pd.DataFrame(comp_rows)
comp.to_csv(f"{HERE}/comparison.csv", index=False)

OUT["healthy_vs_full_cohort"] = {}
for tag in ("healthy_tst240", "healthy_tst300"):
    d = {}
    for name in ("cell1_vs_cell4", "cell2_vs_cell4", "cell3_vs_cell4", "cell2_vs_cell1"):
        s = comp[(comp.healthy_analysis == tag) & (comp.contrast == name) &
                 comp.verdict.notna() & ~comp.negative_control]
        if not len(s):
            d[name] = {"n_comparable": 0}
            continue
        d[name] = {
            "n_comparable": int(len(s)),
            "LARGER": int((s.verdict == "LARGER").sum()),
            "SMALLER": int((s.verdict == "SMALLER").sum()),
            "UNCHANGED": int((s.verdict == "UNCHANGED").sum()),
            "median_ratio": round(float(s.ratio_healthy_over_full.median()), 3),
            "median_hr_healthy": round(float(s.healthy_hr.median()), 3),
            "median_hr_full": round(float(s.full_cohort_hr.median()), 3),
            "n_full_hr_outside_healthy_ci": int((~s.full_hr_inside_healthy_ci).sum()),
            "per_condition": [{"disease": str(x.disease), "healthy_hr": float(x.healthy_hr),
                               "healthy_ci": [float(x.healthy_lo), float(x.healthy_hi)],
                               "full_hr": float(x.full_cohort_hr),
                               "ratio": float(x.ratio_healthy_over_full),
                               "verdict": str(x.verdict)}
                              for _, x in s.sort_values("ratio_healthy_over_full",
                                                        ascending=False).iterrows()],
        }
    OUT["healthy_vs_full_cohort"][tag] = d

# the four conditions the paper singled out in the wake-versus-sleep healthy analysis
WATCH = {"obesity": 1.64, "diabetes": 1.68, "hf": 1.85, "resp_failure": 2.16}
OUT["paper_watchlist"] = {}
for k, v in WATCH.items():
    e = {"paper_wake_sleep_healthy_hr": v}
    for tag in ("healthy_tst240", "healthy_tst300"):
        s = comp[(comp.healthy_analysis == tag) & (comp.key == k)]
        e[tag] = {str(x.contrast): {"healthy_hr": (None if x.healthy_hr is None or
                                                   pd.isna(x.healthy_hr) else float(x.healthy_hr)),
                                    "full_hr": float(x.full_cohort_hr),
                                    "verdict": (None if x.verdict is None or
                                                (isinstance(x.verdict, float) and pd.isna(x.verdict))
                                                else str(x.verdict))}
                  for _, x in s.iterrows()}
    OUT["paper_watchlist"][k] = e

# ======================================================================================
# 9. WRITE
# ======================================================================================
res = pd.concat([panels[t] for t in ("healthy_tst240", "healthy_tst300", "rec_ge360",
                                     "age_35_75")], ignore_index=True)
res.to_csv(f"{HERE}/results.csv", index=False)

OUT["dropped_for_too_few_events"] = {
    tag: {"nothing_estimable": [{"disease": s["disease"], "reason": s["reason"],
                                 "events": s.get("events_total")}
                                for s in skips[tag]],
          "no_fourcell_fit": [{"disease": str(r.disease),
                               "reason": str(r.get("fourcell_reason", "")),
                               "events_total": int(r.events_total),
                               "ev_by_cell": {c: int(r[f"ev_{c}"]) for c in CELLS}}
                              for _, r in panels[tag].iterrows() if not r.fourcell_fitted],
          "no_cell2_vs_cell1_pair": [str(r.disease) for _, r in panels[tag].iterrows()
                                     if pd.isna(r.get("p_cell2_vs_cell1_hr", np.nan))],
          "no_cell1_vs_cell4_pair": [str(r.disease) for _, r in panels[tag].iterrows()
                                     if pd.isna(r.get("p_cell1_vs_cell4_hr", np.nan))]}
    for tag in panels}

# sensitivity agreement, primary against each restricted panel
OUT["sensitivity_agreement"] = {}
for senstag in ("rec_ge360", "age_35_75"):
  mg = panels[PRIMTAG].merge(panels[senstag], on="key", suffixes=("_p", "_s"))
  OUT["sensitivity_agreement"][senstag] = {}
  for a, ref in PAIRS:
    pre = f"p_{a}_vs_{ref}"
    m2 = mg[mg[f"{pre}_hr_p"].notna() & mg[f"{pre}_hr_s"].notna()]
    if not len(m2):
        OUT["sensitivity_agreement"][senstag][f"{a}_vs_{ref}"] = {"n_paired": 0}
        continue
    ap = (m2[f"{pre}_lo_p"] > 1) | (m2[f"{pre}_hi_p"] < 1)
    cs = (m2[f"{pre}_lo_s"] > 1) | (m2[f"{pre}_hi_s"] < 1)
    OUT["sensitivity_agreement"][senstag][f"{a}_vs_{ref}"] = {
        "n_paired": int(len(m2)), "sig_primary": int(ap.sum()), "sig_sensitivity": int(cs.sum()),
        "sig_both": int((ap & cs).sum()),
        "median_hr_primary": round(float(m2[f"{pre}_hr_p"].median()), 3),
        "median_hr_sensitivity": round(float(m2[f"{pre}_hr_s"].median()), 3),
        "median_abs_pct_change_hr": round(float(
            (100 * (m2[f"{pre}_hr_s"] / m2[f"{pre}_hr_p"] - 1)).abs().median()), 1),
        "spearman_hr": (round(float(stats.spearmanr(m2[f"{pre}_hr_p"],
                                                    m2[f"{pre}_hr_s"]).statistic), 3)
                        if len(m2) > 2 else None),
        "lost": [str(x) for x in m2.loc[ap & ~cs, "disease_p"]],
        "gained": [str(x) for x in m2.loc[~ap & cs, "disease_p"]],
    }

json.dump(OUT, open(f"{HERE}/summary.json", "w"), indent=1, default=str)

# ======================================================================================
# 10. READ-OUT
# ======================================================================================
log("\n" + "=" * 95)
log(f"PRIMARY = {PRIMTAG}   healthy n = {N_HEALTHY:,}   four-cell n = {len(FPRIM):,}")
for name in ("cell1_vs_cell4", "cell2_vs_cell4", "cell3_vs_cell4", "cell2_vs_cell1"):
    s = OUT["summary"][PRIMTAG][name]
    fl = OUT["negative_control_floor"][PRIMTAG][name]
    log(f"\n{name}   {s['n_significant']} of {s['n_conditions_tested']} significant "
        f"(chance {s['n_expected_by_chance']}), {s['n_significant_fdr']} after FDR, "
        f"median HR {s['median_hr']}, range {s['hr_range'][0]}-{s['hr_range'][1]}")
    if fl:
        log(f"  negative-control floor {fl['max_abs_hr']} ({fl['at_condition']} "
            f"HR {fl['hr_at']}), {fl['n_controls']} controls span "
            f"{fl['range'][0]} to {fl['range'][1]}, {fl['n_controls_significant']} significant")
    for x in s["significant"][:15]:
        log(f"    {x['disease']:30s} HR {x['hr']:.2f} [{x['ci'][0]:.2f}-{x['ci'][1]:.2f}] "
            f"p={x['p']:.2g} q={x['q']:.3f} clears_floor={x['clears_floor']}")

log("\n" + "-" * 95)
log("healthy against the full cohort")
for tag in ("healthy_tst240", "healthy_tst300"):
    for name, d in OUT["healthy_vs_full_cohort"][tag].items():
        if not d.get("n_comparable"):
            continue
        log(f"  {tag:16s} {name:16s} n={d['n_comparable']:2d}  LARGER {d['LARGER']}, "
            f"SMALLER {d['SMALLER']}, UNCHANGED {d['UNCHANGED']}, median ratio "
            f"{d['median_ratio']}, median HR {d['median_hr_healthy']} vs "
            f"{d['median_hr_full']} in the full cohort")

for senstag, lab in (("rec_ge360", "recordings of at least 6 hours"),
                     ("age_35_75", "ages 35 to 75, where all four cells overlap")):
    log(f"\nsensitivity, {lab}")
    for k, a in OUT["sensitivity_agreement"][senstag].items():
        if not a.get("n_paired"):
            continue
        log(f"  {k:16s} significant {a['sig_primary']} -> {a['sig_sensitivity']}, "
            f"both {a['sig_both']}, median HR {a['median_hr_primary']} -> "
            f"{a['median_hr_sensitivity']}, median |change| {a['median_abs_pct_change_hr']}%")
    fl = OUT["negative_control_floor"][senstag]
    for k in ("cell1_vs_cell4", "cell2_vs_cell1"):
        if fl.get(k):
            log(f"  floor {k}: {fl[k]['max_abs_hr']} ({fl[k]['at_condition']}), "
                f"controls span {fl[k]['range'][0]} to {fl[k]['range'][1]}")
    for k in ("cell1_vs_cell4", "cell2_vs_cell1", "cell2_vs_cell4", "cell3_vs_cell4"):
        s = OUT["summary"][senstag].get(k)
        if s and s["n_conditions_tested"]:
            log(f"  {k:16s} {s['n_significant']}/{s['n_conditions_tested']} significant, "
                f"{s['n_significant_fdr']} FDR, median HR {s['median_hr']}")

log("\nwrote cells.csv, results.csv, interaction.csv, comparison.csv, summary.json")
