"""
Sleep duration crossed with nocturnal oxygenation: do the disease associations hold in each cell?

The question this answers is whether a short night carries any disease risk when the oxygen is
normal, and whether low oxygen still predicts disease in people who barely slept. The existing
discordance analysis (T90_Manuscript/numbers/run_phenotype_v2.py) crosses T90 with the AHI, not
with sleep duration, so this cross has never been run.

Cells, on the 19,173 analysis cohort:
    cell1  short sleep  + normal oxygen   TST < 240 min, T90 <= 1%
    cell2  short sleep  + low oxygen      TST < 240 min, T90 > 10%
    cell3  normal sleep + low oxygen      TST >= 360 min, T90 > 10%
    cell4  normal sleep + normal oxygen   TST >= 360 min, T90 <= 1%    REFERENCE

The 1-10% T90 band and the 240-360 min TST band are left out so the contrast is clean. The
discarded counts are reported.

TST here is one laboratory night. It is not habitual sleep duration, and a short laboratory TST
can mean a short recording rather than little sleep, which is why every model is repeated on
recordings of at least 360 minutes.

Modelling: site-stratified Cox, unpenalized, natural cubic age spline with 4 df of which one
column is dropped so the design is full rank, plus sex. PSG date is time zero.

Writes cells.csv and results.csv in this directory.
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
HERE = f"{paths.SV_ROOT}/dur_x_oxygen"
sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
from cohort_spec import apply_cohort                      # noqa: E402
from disease_definitions import (DISEASES, NEGATIVE_CONTROLS,  # noqa: E402
                                 CIRCULAR, ORGAN_GROUP)

T90_LOW, T90_NORMAL = 10.0, 1.0        # low oxygen = T90 > 10%, normal oxygen = T90 <= 1%
TST_SHORT, TST_NORMAL = 240.0, 360.0   # short sleep < 4 h, normal sleep >= 6 h
MIN_EVENTS = 20                        # per condition, summed over the four cells
REC_MIN = 360.0                        # sensitivity: recordings of at least 6 h

OUT = {}


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
normS = b.TST_min >= TST_NORMAL

b["cell"] = np.select(
    [shortS & normO2, shortS & lowO2, normS & lowO2, normS & normO2],
    ["cell1", "cell2", "cell3", "cell4"], default="excluded")

CELL_LABEL = {
    "cell1": "short sleep, normal oxygen",
    "cell2": "short sleep, low oxygen",
    "cell3": "normal sleep, low oxygen",
    "cell4": "normal sleep, normal oxygen (reference)",
}

OUT["thresholds"] = {"t90_low_gt_pct": T90_LOW, "t90_normal_le_pct": T90_NORMAL,
                     "tst_short_lt_min": TST_SHORT, "tst_normal_ge_min": TST_NORMAL}
OUT["discards"] = {
    "cohort": int(len(b)),
    "t90_middle_1_to_10_pct": int(((b.spo2_pct_below_90 > T90_NORMAL) &
                                   (b.spo2_pct_below_90 <= T90_LOW)).sum()),
    "tst_middle_240_to_360_min": int(((b.TST_min >= TST_SHORT) & (b.TST_min < TST_NORMAL)).sum()),
    "tst_missing": int(b.TST_min.isna().sum()),
    "in_four_cells": int((b.cell != "excluded").sum()),
    "excluded_total": int((b.cell == "excluded").sum()),
}
# alternative short-sleep thresholds, so the cut can be moved without rerunning the extraction
OUT["alternative_short_thresholds"] = {}
for thr, lab in ((180.0, "3h"), (240.0, "4h"), (300.0, "5h")):
    s = b.TST_min < thr
    OUT["alternative_short_thresholds"][lab] = {
        "tst_lt_min": thr, "n_short_total": int(s.sum()),
        "cell1_short_normalO2": int((s & normO2).sum()),
        "cell2_short_lowO2": int((s & lowO2).sum())}

log(f"cells: " + ", ".join(f"{k}={int((b.cell == k).sum()):,}" for k in CELL_LABEL))
log(f"discarded {OUT['discards']['excluded_total']:,} "
    f"({OUT['discards']['t90_middle_1_to_10_pct']:,} in the 1-10% T90 band, "
    f"{OUT['discards']['tst_middle_240_to_360_min']:,} in the 240-360 min TST band, "
    f"{OUT['discards']['tst_missing']} with no TST)")

f4 = b[b.cell != "excluded"].copy()

# ======================================================================================
# 2. CELL PROFILE  ->  cells.csv
# ======================================================================================
PREV_CHECK = ["copd2_prevalent", "hf_prevalent", "cancer_any_prevalent", "obesity_prevalent",
              "osa_prevalent", "dementia_prevalent"]
prof = []
for k in ["cell1", "cell2", "cell3", "cell4"]:
    s = f4[f4.cell == k]
    row = {
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
        "sleep_eff_pct_median": round(float(s.sleep_efficiency_pct.median()), 1),
        "sleep_eff_pct_q1": round(float(s.sleep_efficiency_pct.quantile(.25)), 1),
        "sleep_eff_pct_q3": round(float(s.sleep_efficiency_pct.quantile(.75)), 1),
        "ahi_median": round(float(s.AHI.median()), 2),
        "arousal_index_median": round(float(s.arousal_index.median()), 2),
        "spo2_mean_median": round(float(s.spo2_mean.median()), 2),
        "n_deaths": int(s.death_incident.sum()),
        "death_rate_pct": round(100 * float(s.death_incident.mean()), 2),
        "followup_years_median": round(float(s.death_years.median()), 2),
        "n_recording_ge360": int((s.recording_dur_min >= REC_MIN).sum()),
    }
    for c in PREV_CHECK:
        if c in s.columns:
            row[c.replace("_prevalent", "_prev_pct")] = round(100 * float(s[c].mean()), 1)
    for site, cnt in s.site_id.value_counts().items():
        row[f"site_{site}"] = int(cnt)
    prof.append(row)
cells = pd.DataFrame(prof)
cells.to_csv(f"{HERE}/cells.csv", index=False)
log("\ncell profile")
log(cells[["cell", "n", "age_median", "pct_male", "tst_min_median",
           "recording_dur_min_median", "sleep_eff_pct_median", "t90_median",
           "ahi_median", "n_deaths"]].to_string(index=False))

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


CELLD = ["cell1", "cell2", "cell3"]

OUTCOMES = [(k, v[0]) for k, v in DISEASES.items() if k not in CIRCULAR]
OUTCOMES.append(("death", "Death from any cause"))


def frame_for(d, key):
    """Incident frame for one condition: no prevalent disease, positive follow-up."""
    yc, ec, pc = f"{key}_years", f"{key}_incident", f"{key}_prevalent"
    if not all(c in d.columns for c in (yc, ec, pc)):
        return None
    f = d[(d[pc] == 0) & d[yc].notna() & (d[yc] > 0)].copy()
    sp = agespline(f)
    for c in sp.columns:
        f[c] = sp[c]
    cols = {"T": f[yc], "E": f[ec].astype(int), "site": f.site_id, "male": f.male,
            "cell_lab": f.cell}
    for c in sp.columns:
        cols[c] = f[c]
    out = pd.DataFrame(cols).dropna()
    for k2 in CELLD:
        out[k2] = (out.cell_lab == k2).astype(float)
    return out


def fit_cells(dd):
    """Site-stratified unpenalized Cox with the three non-reference cell indicators."""
    cov = CELLD + ["male"] + [c for c in dd.columns if c.startswith("age_s")]
    m = CoxPHFitter(penalizer=0.0).fit(dd[cov + ["T", "E", "site"]], "T", "E", strata=["site"])
    return m


def contrast(m, a, ref):
    """HR for cell a against cell ref, both non-baseline, from the fitted covariance."""
    V = m.variance_matrix_
    beta = m.params_
    d = float(beta[a] - beta[ref])
    v = float(V.loc[a, a] + V.loc[ref, ref] - 2 * V.loc[a, ref])
    se = np.sqrt(v)
    z = d / se
    return (round(float(np.exp(d)), 3), round(float(np.exp(d - 1.96 * se)), 3),
            round(float(np.exp(d + 1.96 * se)), 3), float(2 * stats.norm.sf(abs(z))))


def run_panel(d, tag):
    """Every non-circular condition with at least MIN_EVENTS events in the four cells."""
    rows, skipped = [], []
    for key, label in OUTCOMES:
        dd = frame_for(d, key)
        if dd is None:
            skipped.append({"key": key, "reason": "columns absent"})
            continue
        if dd.E.sum() < MIN_EVENTS:
            skipped.append({"key": key, "reason": f"{int(dd.E.sum())} events, below {MIN_EVENTS}"})
            continue
        ev = {k2: int(dd[dd.cell_lab == k2].E.sum()) for k2 in
              ["cell1", "cell2", "cell3", "cell4"]}
        nn = {k2: int((dd.cell_lab == k2).sum()) for k2 in
              ["cell1", "cell2", "cell3", "cell4"]}
        if ev["cell4"] == 0 or min(ev[k2] for k2 in CELLD) == 0:
            skipped.append({"key": key, "reason": "a cell has no events", "events": ev})
            continue
        try:
            m = fit_cells(dd)
        except Exception as ex:
            skipped.append({"key": key, "reason": f"fit failed: {type(ex).__name__}"})
            continue
        rec = {"analysis": tag, "key": key, "disease": label,
               "organ_group": ORGAN_GROUP.get(key, "Death"),
               "negative_control": key in NEGATIVE_CONTROLS,
               "circular_flag": key == "obesity_hypovent",
               "n_total": int(len(dd)), "events_total": int(dd.E.sum()),
               "n_cell1": nn["cell1"], "n_cell2": nn["cell2"],
               "n_cell3": nn["cell3"], "n_cell4": nn["cell4"],
               "ev_cell1": ev["cell1"], "ev_cell2": ev["cell2"],
               "ev_cell3": ev["cell3"], "ev_cell4": ev["cell4"],
               "unstable_lt5_events": bool(min(ev.values()) < 5)}
        for k2 in CELLD:
            s = m.summary.loc[k2]
            rec[f"{k2}_hr"] = round(float(s["exp(coef)"]), 3)
            rec[f"{k2}_lo"] = round(float(s["exp(coef) lower 95%"]), 3)
            rec[f"{k2}_hi"] = round(float(s["exp(coef) upper 95%"]), 3)
            rec[f"{k2}_p"] = float(s["p"])
        # the comparison Alen asked for: low oxygen among the short sleepers
        hr, lo, hi, p = contrast(m, "cell2", "cell1")
        rec.update({"cell2_vs_cell1_hr": hr, "cell2_vs_cell1_lo": lo,
                    "cell2_vs_cell1_hi": hi, "cell2_vs_cell1_p": p})
        hr, lo, hi, p = contrast(m, "cell2", "cell3")
        rec.update({"cell2_vs_cell3_hr": hr, "cell2_vs_cell3_lo": lo,
                    "cell2_vs_cell3_hi": hi, "cell2_vs_cell3_p": p})
        # multiplicative additivity: does cell2 exceed cell1 x cell3?
        rec["product_c1_c3"] = round(rec["cell1_hr"] * rec["cell3_hr"], 3)
        rec["cell2_over_product"] = round(rec["cell2_hr"] / rec["product_c1_c3"], 3)
        rows.append(rec)
    df = pd.DataFrame(rows)
    if len(df):
        real = ~df.negative_control
        for k2 in CELLD:
            df[f"{k2}_q"] = np.nan
            q = multipletests(df.loc[real, f"{k2}_p"], method="fdr_bh")[1]
            df.loc[real, f"{k2}_q"] = q
    log(f"\n[{tag}] {len(df)} conditions fitted, {len(skipped)} skipped")
    return df, skipped


def floor_of(df):
    """Confounding floor: the largest negative-control deviation from 1 in each cell."""
    neg = df[df.negative_control]
    out = {}
    for k2 in CELLD:
        h = neg[f"{k2}_hr"].dropna()
        if not len(h):
            continue
        dev = h.apply(lambda x: max(x, 1 / x))
        i = dev.idxmax()
        out[k2] = {"max_abs_hr": round(float(dev.max()), 3),
                   "at_condition": str(neg.loc[i, "disease"]),
                   "hr_at": float(neg.loc[i, f"{k2}_hr"]),
                   "ci_at": [float(neg.loc[i, f"{k2}_lo"]), float(neg.loc[i, f"{k2}_hi"])],
                   "range": [round(float(h.min()), 3), round(float(h.max()), 3)],
                   "n_controls": int(len(h)),
                   "n_controls_significant": int(((neg[f"{k2}_lo"] > 1) |
                                                  (neg[f"{k2}_hi"] < 1)).sum())}
    return out


# ======================================================================================
# 4. PRIMARY PANEL
# ======================================================================================
prim, prim_skip = run_panel(f4, "primary")
OUT["skipped_primary"] = prim_skip
OUT["negative_control_floor_primary"] = floor_of(prim)

# ======================================================================================
# 5. SENSITIVITY: RECORDINGS OF AT LEAST 6 HOURS
# ======================================================================================
f4r = f4[f4.recording_dur_min >= REC_MIN].copy()
OUT["long_recording_subset"] = {
    "rule": f"recording_dur_min >= {REC_MIN:.0f}",
    "n": int(len(f4r)),
    "by_cell": {k: int((f4r.cell == k).sum()) for k in CELL_LABEL},
    "pct_retained_by_cell": {k: round(100 * float((f4r.cell == k).sum() /
                                                  max((f4.cell == k).sum(), 1)), 1)
                             for k in CELL_LABEL},
}
log("\nlong-recording subset " + json.dumps(OUT["long_recording_subset"]["by_cell"]))
sens, sens_skip = run_panel(f4r, "recording_ge_360min")
OUT["skipped_sensitivity"] = sens_skip
OUT["negative_control_floor_sensitivity"] = floor_of(sens)

# ======================================================================================
# 6. FORMAL INTERACTION, TOP 10 BY EVENT COUNT
# ======================================================================================
top10 = prim.sort_values("events_total", ascending=False).head(10).key.tolist()
inter = []
for key in top10:
    dd = frame_for(f4, key)
    dd["short"] = (dd.cell_lab.isin(["cell1", "cell2"])).astype(float)
    dd["lowo2"] = (dd.cell_lab.isin(["cell2", "cell3"])).astype(float)
    dd["short_x_lowo2"] = dd["short"] * dd["lowo2"]
    cov = ["short", "lowo2", "short_x_lowo2", "male"] + \
          [c for c in dd.columns if c.startswith("age_s")]
    m = CoxPHFitter(penalizer=0.0).fit(dd[cov + ["T", "E", "site"]], "T", "E", strata=["site"])
    s = m.summary
    r = prim[prim.key == key].iloc[0]
    inter.append({
        "key": key, "disease": r.disease, "events": int(dd.E.sum()),
        "main_short_hr": round(float(s.loc["short", "exp(coef)"]), 3),
        "main_lowo2_hr": round(float(s.loc["lowo2", "exp(coef)"]), 3),
        "interaction_hr": round(float(s.loc["short_x_lowo2", "exp(coef)"]), 3),
        "interaction_lo": round(float(s.loc["short_x_lowo2", "exp(coef) lower 95%"]), 3),
        "interaction_hi": round(float(s.loc["short_x_lowo2", "exp(coef) upper 95%"]), 3),
        "interaction_p": float(s.loc["short_x_lowo2", "p"]),
        "cell2_hr": float(r.cell2_hr), "product_c1_c3": float(r.product_c1_c3),
        "cell2_exceeds_product": bool(r.cell2_hr > r.product_c1_c3),
    })
inter = pd.DataFrame(inter)
inter["interaction_q"] = multipletests(inter.interaction_p, method="fdr_bh")[1]
OUT["interaction_top10"] = json.loads(inter.to_json(orient="records"))
log("\ninteraction, top 10 by events")
log(inter[["disease", "events", "interaction_hr", "interaction_p",
           "cell2_hr", "product_c1_c3", "cell2_exceeds_product"]].to_string(index=False))

# ======================================================================================
# 7. WRITE
# ======================================================================================
res = pd.concat([prim, sens], ignore_index=True)
res.to_csv(f"{HERE}/results.csv", index=False)
inter.to_csv(f"{HERE}/interaction.csv", index=False)

# headline summaries
def sig(df, k2):
    r = df[(~df.negative_control)]
    return r[(r[f"{k2}_lo"] > 1) | (r[f"{k2}_hi"] < 1)]


OUT["summary_primary"] = {}
for k2 in CELLD:
    s = sig(prim, k2)
    OUT["summary_primary"][k2] = {
        "label": CELL_LABEL[k2], "n_conditions": int((~prim.negative_control).sum()),
        "n_significant": int(len(s)),
        "n_significant_fdr": int((prim.loc[~prim.negative_control, f"{k2}_q"] < 0.05).sum()),
        "hr_range": [round(float(prim.loc[~prim.negative_control, f"{k2}_hr"].min()), 3),
                     round(float(prim.loc[~prim.negative_control, f"{k2}_hr"].max()), 3)],
        "median_hr": round(float(prim.loc[~prim.negative_control, f"{k2}_hr"].median()), 3),
        "significant": [{"disease": r.disease, "hr": r[f"{k2}_hr"],
                         "ci": [r[f"{k2}_lo"], r[f"{k2}_hi"]], "p": r[f"{k2}_p"],
                         "q": round(float(r[f"{k2}_q"]), 4),
                         "clears_floor": bool(max(r[f"{k2}_hr"], 1 / r[f"{k2}_hr"]) >
                                              OUT["negative_control_floor_primary"][k2]["max_abs_hr"])}
                        for _, r in s.sort_values(f"{k2}_hr", ascending=False).iterrows()],
    }

# cell2 against cell1: does low oxygen still bite in the short sleepers
c21 = prim[(~prim.negative_control)]
OUT["cell2_vs_cell1"] = {
    "n_conditions": int(len(c21)),
    "n_significant": int(((c21.cell2_vs_cell1_lo > 1) | (c21.cell2_vs_cell1_hi < 1)).sum()),
    "median_hr": round(float(c21.cell2_vs_cell1_hr.median()), 3),
    "significant": [{"disease": r.disease, "hr": r.cell2_vs_cell1_hr,
                     "ci": [r.cell2_vs_cell1_lo, r.cell2_vs_cell1_hi],
                     "p": r.cell2_vs_cell1_p, "ev_cell1": int(r.ev_cell1),
                     "ev_cell2": int(r.ev_cell2)}
                    for _, r in c21[(c21.cell2_vs_cell1_lo > 1) | (c21.cell2_vs_cell1_hi < 1)]
                    .sort_values("cell2_vs_cell1_hr", ascending=False).iterrows()],
}

# does the sensitivity analysis change the answer
mg = prim.merge(sens, on="key", suffixes=("_p", "_s"))
OUT["sensitivity_agreement"] = {}
for k2 in CELLD:
    a = (mg[f"{k2}_lo_p"] > 1) | (mg[f"{k2}_hi_p"] < 1)
    c = (mg[f"{k2}_lo_s"] > 1) | (mg[f"{k2}_hi_s"] < 1)
    OUT["sensitivity_agreement"][k2] = {
        "n_paired": int(len(mg)),
        "sig_primary": int(a.sum()), "sig_sensitivity": int(c.sum()),
        "sig_both": int((a & c).sum()),
        "median_abs_pct_change_hr": round(float(
            (100 * (mg[f"{k2}_hr_s"] / mg[f"{k2}_hr_p"] - 1)).abs().median()), 1),
        "spearman_hr": round(float(stats.spearmanr(mg[f"{k2}_hr_p"],
                                                   mg[f"{k2}_hr_s"]).statistic), 3),
        "lost": [str(x) for x in mg.loc[a & ~c, "disease_p"]],
        "gained": [str(x) for x in mg.loc[~a & c, "disease_p"]],
    }

json.dump(OUT, open(f"{HERE}/summary.json", "w"), indent=1, default=str)

log("\n" + "=" * 90)
for k2 in CELLD:
    s = OUT["summary_primary"][k2]
    fl = OUT["negative_control_floor_primary"][k2]
    log(f"\n{k2}  {s['label']}")
    log(f"  {s['n_significant']} of {s['n_conditions']} conditions significant, "
        f"{s['n_significant_fdr']} after FDR. HR range {s['hr_range'][0]} to {s['hr_range'][1]}, "
        f"median {s['median_hr']}")
    log(f"  confounding floor {fl['max_abs_hr']} ({fl['at_condition']} "
        f"HR {fl['hr_at']}), controls span {fl['range'][0]} to {fl['range'][1]}")
    for x in s["significant"][:12]:
        log(f"    {x['disease']:30s} HR {x['hr']:.2f} [{x['ci'][0]:.2f}-{x['ci'][1]:.2f}] "
            f"p={x['p']:.2g} q={x['q']:.3f} clears_floor={x['clears_floor']}")

log(f"\ncell2 vs cell1 (low oxygen among short sleepers): "
    f"{OUT['cell2_vs_cell1']['n_significant']} of {OUT['cell2_vs_cell1']['n_conditions']} "
    f"significant, median HR {OUT['cell2_vs_cell1']['median_hr']}")
for x in OUT["cell2_vs_cell1"]["significant"][:12]:
    log(f"    {x['disease']:30s} HR {x['hr']:.2f} [{x['ci'][0]:.2f}-{x['ci'][1]:.2f}] "
        f"p={x['p']:.2g}  events {x['ev_cell1']} vs {x['ev_cell2']}")

log("\nsensitivity, recordings of at least 6 hours")
for k2 in CELLD:
    a = OUT["sensitivity_agreement"][k2]
    log(f"  {k2}: significant {a['sig_primary']} -> {a['sig_sensitivity']}, both {a['sig_both']}, "
        f"median |change| in HR {a['median_abs_pct_change_hr']}%, rank r {a['spearman_hr']}")

log("\nwrote cells.csv, results.csv, interaction.csv, summary.json")
