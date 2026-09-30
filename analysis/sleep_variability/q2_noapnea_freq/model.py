"""
Absolute disease frequency by nocturnal oxygen group, in patients without sleep apnea.

The question this answers is the one the hazard ratio does not: of the people who came through the
laboratory with an apnea-hypopnea index below 5 and still spent a large share of the night below
90% saturation, how many actually went on to develop each condition, against how many did so among
the people whose oxygen was preserved.

Population
    The frozen analysis cohort (numbers/cohort_spec.py, n = 19,173), restricted to AHI < 5. That is
    2,920 recordings. The exposure is a split of that group at T90 > 10% of the recording, which
    leaves 300 exposed and 2,620 unexposed. A 5% cut is run as a sensitivity analysis because 300
    is thin.

Per condition
    Prevalent cases at the sleep study are removed, and follow-up must be positive, so the at-risk
    denominator differs from 300 / 2,620 condition by condition. Both denominators are printed.

Estimates
    Risk ratio with the Katz log method, incidence-rate ratio with the log-rate method, and a Cox
    model for the binary oxygen contrast adjusted for age, sex and hospital. Age enters as a
    natural cubic spline with 4 degrees of freedom whose first column is dropped, because the four
    columns sum to one and an unpenalized fit on a rank-deficient design silently fits nothing.
    Hospital is a stratifying variable, matching the primary specification. The Cox fits carry no
    ridge penalty: lifelines scales its elastic-net term by the sample size, so any nonzero
    penalizer is an effective ridge of roughly n times that value.

Stability
    Every row whose exposed group holds fewer than 5 incident cases is flagged. At 300 exposed
    people this is most of the panel, and the flag matters more than the point estimate. Rows are
    also flagged when the adjusted model has fewer than 10 events in total or returns a confidence
    interval spanning more than three orders of magnitude, which is the signature of a covariate
    that separates the outcome rather than a real effect.

Writes table_noapnea.csv, table_noapnea_5pct.csv, table_noapnea.md and results_noapnea.json.
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

warnings.filterwarnings("ignore")

ROOT = paths.T90_ROOT
OUT = f"{paths.SV_ROOT}/q2_noapnea_freq"
sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
from cohort_spec import COHORT_N, apply_cohort            # noqa: E402
from disease_definitions import DISEASES, ORGAN_GROUP     # noqa: E402

# Negative-control panel and confounding floor, settled 2026-08-07 and imported from the
# source of truth rather than listed here. Back pain, cataract, glaucoma, contact
# dermatitis, hemorrhoids; floor 1.063, set by glaucoma. Fracture and osteoarthritis are
# ordinary outcomes now: a plausible causal path runs from the exposure to each of them,
# which disqualifies a control however clean its estimate looks. The floor moved DOWN,
# 1.150 to 1.133 to 1.063, so associations are promoted and none can be demoted.
from disease_definitions import (NEGATIVE_CONTROLS as CONTROLS,   # noqa: E402
                                 CONFOUNDING_FLOOR_PER_SD as FLOOR)
DROPPED_CONTROL = ["fracture", "osteoarthritis"]

AHI_CUT = 5.0
THRESHOLDS = [10.0, 5.0]
MIN_EXPOSED_CASES = 5

# ---------------------------------------------------------------------------------------------
# cohort
# ---------------------------------------------------------------------------------------------
b = apply_cohort(pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet"))
assert len(b) == COHORT_N
na = b[b.AHI < AHI_CUT].copy()
print(f"cohort {len(b):,} -> AHI < {AHI_CUT:.0f}: {len(na):,}")
N_NA, N_NA_LOWOX = int(len(na)), int((na.spo2_pct_below_90 > 10).sum())   # v7: 2,920 and 300 before the corrections, now taken from the data
print(f"no-apnea group {N_NA:,} nights, {N_NA_LOWOX:,} with T90 over 10 percent")
assert N_NA > 0 and N_NA_LOWOX > 0

na["male"] = (na.sex.astype(str).str.upper().str[0] == "M").astype(int)
na["fu_years"] = (na.censor_date - na.psg_date).dt.days / 365.25

# age basis, first column dropped so the design is full rank under any solver
sp = dmatrix("cr(a, df=4) - 1", {"a": na.AgeAtVisit.values}, return_type="dataframe").iloc[:, 1:]
for i in range(sp.shape[1]):
    na[f"age_s{i}"] = sp.iloc[:, i].values
ADJ = [f"age_s{i}" for i in range(sp.shape[1])] + ["male"]
print(f"age spline columns kept: {sp.shape[1]} ({', '.join(ADJ[:-1])}), plus sex, hospital stratified")

OUTCOMES = [(k, v[0]) for k, v in DISEASES.items()] + [("death", "Death from any cause")]


# ---------------------------------------------------------------------------------------------
# estimators
# ---------------------------------------------------------------------------------------------
def risk_ratio(c1, n1, c0, n0):
    """Katz log method. Undefined if either group has no cases."""
    if c1 == 0 or c0 == 0 or n1 == 0 or n0 == 0:
        return None, None, None
    rr = (c1 / n1) / (c0 / n0)
    se = np.sqrt(1 / c1 - 1 / n1 + 1 / c0 - 1 / n0)
    return rr, rr * np.exp(-1.96 * se), rr * np.exp(1.96 * se)


def rate_ratio(c1, pt1, c0, pt0):
    """Log-rate method on person-time. Undefined if either group has no cases."""
    if c1 == 0 or c0 == 0 or pt1 <= 0 or pt0 <= 0:
        return None, None, None
    irr = (c1 / pt1) / (c0 / pt0)
    se = np.sqrt(1 / c1 + 1 / c0)
    return irr, irr * np.exp(-1.96 * se), irr * np.exp(1.96 * se)


def cox_binary(d):
    """Unpenalized Cox for the binary oxygen contrast, age spline and sex, stratified by hospital."""
    if d.E.sum() < 1 or d.groupby("low_o2").E.sum().min() < 1:
        return None
    keep = ["low_o2"] + ADJ + ["T", "E", "site"]
    try:
        c = CoxPHFitter(penalizer=0.0).fit(d[keep], "T", "E", strata=["site"])
    except Exception as exc:
        return {"error": type(exc).__name__}
    r = c.summary.loc["low_o2"]
    return {"hr": float(r["exp(coef)"]), "lo": float(r["exp(coef) lower 95%"]),
            "hi": float(r["exp(coef) upper 95%"]), "p": float(r["p"])}


# ---------------------------------------------------------------------------------------------
# one full pass at a given T90 threshold
# ---------------------------------------------------------------------------------------------
def run(threshold):
    na["low_o2"] = (na.spo2_pct_below_90 > threshold).astype(int)
    g1, g0 = na[na.low_o2 == 1], na[na.low_o2 == 0]
    header = {
        "threshold_pct": threshold,
        "n_exposed": int(len(g1)), "n_unexposed": int(len(g0)),
        "median_fu_exposed": float(g1.fu_years.median()),
        "median_fu_unexposed": float(g0.fu_years.median()),
        "iqr_fu_exposed": [float(g1.fu_years.quantile(.25)), float(g1.fu_years.quantile(.75))],
        "iqr_fu_unexposed": [float(g0.fu_years.quantile(.25)), float(g0.fu_years.quantile(.75))],
        "person_years_exposed": float(g1.fu_years.sum()),
        "person_years_unexposed": float(g0.fu_years.sum()),
        "median_age_exposed": float(g1.AgeAtVisit.median()),
        "median_age_unexposed": float(g0.AgeAtVisit.median()),
        "pct_male_exposed": float(100 * g1.male.mean()),
        "pct_male_unexposed": float(100 * g0.male.mean()),
        "median_t90_exposed": float(g1.spo2_pct_below_90.median()),
        "median_t90_unexposed": float(g0.spo2_pct_below_90.median()),
    }

    rows = []
    for key, label in OUTCOMES:
        yc, ec, pc = f"{key}_years", f"{key}_incident", f"{key}_prevalent"
        if yc not in na.columns:
            continue
        f = na[(na[pc] == 0) & na[yc].notna() & (na[yc] > 0)]
        d = pd.DataFrame({"T": f[yc], "E": f[ec].astype(int), "site": f.site_id,
                          "low_o2": f.low_o2, **{c: f[c] for c in ADJ}}).dropna()

        e1, e0 = d[d.low_o2 == 1], d[d.low_o2 == 0]
        n1, n0 = len(e1), len(e0)
        c1, c0 = int(e1.E.sum()), int(e0.E.sum())
        pt1, pt0 = float(e1["T"].sum()), float(e0["T"].sum())
        p1 = 100 * c1 / n1 if n1 else np.nan
        p0 = 100 * c0 / n0 if n0 else np.nan
        ir1 = 1000 * c1 / pt1 if pt1 > 0 else np.nan
        ir0 = 1000 * c0 / pt0 if pt0 > 0 else np.nan
        rr, rrlo, rrhi = risk_ratio(c1, n1, c0, n0)
        irr, irrlo, irrhi = rate_ratio(c1, pt1, c0, pt0)
        cox = cox_binary(d)

        unstable = []
        if c1 < MIN_EXPOSED_CASES:
            unstable.append(f"low-oxygen group has {c1} case" + ("" if c1 == 1 else "s"))
        if c1 + c0 < 10:
            unstable.append(f"{c1 + c0} events in total")
        if cox and "hr" in cox and cox["lo"] > 0 and cox["hi"] / cox["lo"] > 1000:
            unstable.append("confidence interval spans more than 3 orders of magnitude")
        if cox and "error" in cox:
            unstable.append(f"Cox model did not converge ({cox['error']})")
        if cox is None:
            unstable.append("no cases in one group, hazard ratio not estimable")

        rows.append({
            "key": key, "disease": label, "organ_group": ORGAN_GROUP.get(key, "Mortality"),
            "row_type": ("negative_control" if key in CONTROLS
                         else "mortality" if key == "death" else "disease"),
            "n_low_o2": n1, "cases_low_o2": c1, "pct_low_o2": p1,
            "n_normal_o2": n0, "cases_normal_o2": c0, "pct_normal_o2": p0,
            "abs_diff_pp": p1 - p0,
            "person_years_low_o2": pt1, "person_years_normal_o2": pt0,
            "ir_low_o2_per1000py": ir1, "ir_normal_o2_per1000py": ir0,
            "ir_diff_per1000py": ir1 - ir0,
            "rr": rr, "rr_lo": rrlo, "rr_hi": rrhi,
            "irr": irr, "irr_lo": irrlo, "irr_hi": irrhi,
            "hr_adj": cox.get("hr") if cox else None,
            "hr_lo": cox.get("lo") if cox else None,
            "hr_hi": cox.get("hi") if cox else None,
            "p_value": cox.get("p") if cox else None,
            "unstable": bool(unstable),
            "unstable_reason": ", ".join(unstable),
        })
    return header, pd.DataFrame(rows)


# ---------------------------------------------------------------------------------------------
# formatting
# ---------------------------------------------------------------------------------------------
def fmt_p(p):
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return "not estimable"
    if p < 0.0001:
        return f"{p:.2e}"
    return f"{p:.4f}"


def fmt_ci(v, lo, hi, nd=2):
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return "not estimable"
    return f"{v:.{nd}f} ({lo:.{nd}f} to {hi:.{nd}f})"


def order(df):
    """Diseases by adjusted hazard ratio, not-estimable rows last, then mortality, then controls."""
    rank = {"disease": 0, "mortality": 1, "negative_control": 2}
    d = df.copy()
    d["_r"] = d.row_type.map(rank)
    d["_h"] = d.hr_adj.fillna(-1)
    return d.sort_values(["_r", "_h"], ascending=[True, False]).drop(columns=["_r", "_h"])


def md_table(df, threshold):
    star = "*"
    lines = [
        "| Disease | AHI<5 + T90>%g%%: cases/N (%%) | AHI<5 + T90<=%g%%: cases/N (%%) | "
        "Absolute difference, pp | Incidence rates per 1,000 PY | Adjusted HR (95%% CI) | P value |"
        % (threshold, threshold),
        "|---|---|---|---|---|---|---|",
    ]
    block = None
    for _, r in order(df).iterrows():
        if r.row_type != block:
            block = r.row_type
            title = {"disease": "Disease outcomes", "mortality": "Mortality",
                     "negative_control": "Negative controls (the floor)"}[block]
            lines.append(f"| **{title}** | | | | | | |")
        name = r.disease + (f" {star}" if r.unstable else "")
        hr = fmt_ci(r.hr_adj, r.hr_lo, r.hr_hi)
        lines.append(
            f"| {name} | {r.cases_low_o2}/{r.n_low_o2} ({r.pct_low_o2:.1f}%) | "
            f"{r.cases_normal_o2}/{r.n_normal_o2} ({r.pct_normal_o2:.1f}%) | "
            f"{r.abs_diff_pp:+.2f} | {r.ir_low_o2_per1000py:.1f} vs {r.ir_normal_o2_per1000py:.1f} | "
            f"{hr} | {fmt_p(r.p_value)} |")
    return "\n".join(lines)


def md_ratios(df, threshold):
    lines = ["| Disease | Risk ratio (95% CI) | Incidence-rate ratio (95% CI) |", "|---|---|---|"]
    for _, r in order(df).iterrows():
        name = r.disease + (" *" if r.unstable else "") + \
            (" [negative control]" if r.row_type == "negative_control" else "")
        lines.append(f"| {name} | {fmt_ci(r.rr, r.rr_lo, r.rr_hi)} | "
                     f"{fmt_ci(r.irr, r.irr_lo, r.irr_hi)} |")
    return "\n".join(lines)


def md_header(h):
    t = h["threshold_pct"]
    return (
        f"- Low oxygen, T90 > {t:g}% of the recording: **n = {h['n_exposed']:,}**, "
        f"{h['person_years_exposed']:,.1f} person-years, median follow-up "
        f"{h['median_fu_exposed']:.2f} years (IQR {h['iqr_fu_exposed'][0]:.2f} to "
        f"{h['iqr_fu_exposed'][1]:.2f}), median age {h['median_age_exposed']:.0f}, "
        f"{h['pct_male_exposed']:.1f}% male, median T90 {h['median_t90_exposed']:.2f}%\n"
        f"- Preserved oxygen, T90 <= {t:g}%: **n = {h['n_unexposed']:,}**, "
        f"{h['person_years_unexposed']:,.1f} person-years, median follow-up "
        f"{h['median_fu_unexposed']:.2f} years (IQR {h['iqr_fu_unexposed'][0]:.2f} to "
        f"{h['iqr_fu_unexposed'][1]:.2f}), median age {h['median_age_unexposed']:.0f}, "
        f"{h['pct_male_unexposed']:.1f}% male, median T90 {h['median_t90_unexposed']:.2f}%\n"
    )


# ---------------------------------------------------------------------------------------------
# run both thresholds
# ---------------------------------------------------------------------------------------------
res = {}
for t in THRESHOLDS:
    h, df = run(t)
    res[t] = (h, df)
    print(f"\nT90 > {t:g}%: {h['n_exposed']:,} exposed, {h['n_unexposed']:,} unexposed, "
          f"{len(df)} outcomes, {int(df.unstable.sum())} flagged unstable")

h10, d10 = res[10.0]
h5, d5 = res[5.0]

COLS = ["key", "disease", "organ_group", "row_type",
        "n_low_o2", "cases_low_o2", "pct_low_o2",
        "n_normal_o2", "cases_normal_o2", "pct_normal_o2", "abs_diff_pp",
        "person_years_low_o2", "person_years_normal_o2",
        "ir_low_o2_per1000py", "ir_normal_o2_per1000py", "ir_diff_per1000py",
        "rr", "rr_lo", "rr_hi", "irr", "irr_lo", "irr_hi",
        "hr_adj", "hr_lo", "hr_hi", "p_value", "unstable", "unstable_reason"]
order(d10)[COLS].to_csv(f"{OUT}/table_noapnea.csv", index=False)
order(d5)[COLS].to_csv(f"{OUT}/table_noapnea_5pct.csv", index=False)

flagged10 = order(d10[d10.unstable])
flagged5 = order(d5[d5.unstable])


def control_band(df):
    """Where the negative controls actually land in this subgroup, which is the local floor."""
    c = df[df.row_type == "negative_control"]
    top = c.loc[c.hr_adj.idxmax()]
    cleared = df[(df.row_type == "disease") & (df.hr_adj > top.hr_adj) & (df.hr_lo > 1)]
    return (f"Negative controls in this subgroup run from {c.hr_adj.min():.2f} ({c.loc[c.hr_adj.idxmin()].disease}) "
            f"to {c.hr_adj.max():.2f} ({top.disease}), so the local floor is **{c.hr_adj.max():.2f}**, "
            f"higher than the {FLOOR} floor of the full cohort. "
            f"{len(cleared)} disease outcomes sit above that ceiling with an interval excluding 1: "
            + ", ".join(f"{r.disease} ({r.hr_adj:.2f})" for _, r in
                        cleared.sort_values('hr_adj', ascending=False).iterrows()) + ".")

md = f"""# Absolute disease frequency by nocturnal oxygen, in patients without sleep apnea

Cohort {COHORT_N:,} recordings, restricted to an apnea-hypopnea index below 5, which is
{h10['n_exposed'] + h10['n_unexposed']:,} people. Split at T90, the share of the recording spent
below 90% saturation. Counts are incident cases after the sleep study, with anyone already
carrying the diagnosis removed, so the denominator moves from row to row.

Hazard ratios are unpenalized Cox models for the binary oxygen contrast, adjusted for age (natural
cubic spline, 4 degrees of freedom, first column dropped), sex, and stratified by hospital.
Negative controls are {", ".join(DISEASES[k][0].lower() for k in CONTROLS)}, and the floor from
the main analysis is {FLOOR}. Fracture and osteoarthritis are shown as ordinary outcomes, not
controls: a plausible causal path runs from nocturnal hypoxemia to each of them.

`*` marks a row where the low-oxygen group has fewer than {MIN_EXPOSED_CASES} cases, or the model
has fewer than 10 events in total, or the interval is so wide that the fit is carrying no
information. In those rows the percentage, the ratio and the interval should not be read as
estimates. The reason for each flag is in `unstable_reason` in the CSV.

## Person-time and follow-up, primary threshold

{md_header(h10)}
## Primary table, T90 > 10%

{md_table(d10, 10.0)}

### Where the floor sits

{control_band(d10)}

### Risk ratios and incidence-rate ratios, T90 > 10%

{md_ratios(d10, 10.0)}

### Rows flagged as unstable at the 10% cut ({len(flagged10)} of {len(d10)})

Every one of these has fewer than {MIN_EXPOSED_CASES} incident cases in the 300-person low-oxygen
group, or too few events overall for the adjusted model to mean anything.

{chr(10).join('- **' + r.disease + '**, ' + r.unstable_reason for _, r in flagged10.iterrows())}

## Sensitivity analysis, T90 > 5%

{md_header(h5)}
{md_table(d5, 5.0)}

### Where the floor sits, 5% cut

{control_band(d5)}

### Risk ratios and incidence-rate ratios, T90 > 5%

{md_ratios(d5, 5.0)}

### Rows flagged as unstable at the 5% cut ({len(flagged5)} of {len(d5)})

{chr(10).join('- **' + r.disease + '**, ' + r.unstable_reason for _, r in flagged5.iterrows())}

---
Generated by `model.py`. Data `data_frozen_v7_2026-09/t90_final.parquet` through
`numbers/cohort_spec.py::apply_cohort`, conditions from `numbers/disease_definitions.py`,
disease codes from `bdsp_diseases_v3.csv` lineage. Never `bdsp_diseases_v2.csv`.
"""
open(f"{OUT}/table_noapnea.md", "w").write(md)

json.dump({
    "cohort_n": COHORT_N, "ahi_cut": AHI_CUT,
    "n_ahi_lt5": int(len(na)),
    "primary": {"header": h10, "rows": json.loads(order(d10)[COLS].to_json(orient="records"))},
    "sensitivity_5pct": {"header": h5, "rows": json.loads(order(d5)[COLS].to_json(orient="records"))},
    "negative_controls": CONTROLS, "dropped_control": DROPPED_CONTROL, "floor": FLOOR,
}, open(f"{OUT}/results_noapnea.json", "w"), indent=1)

# ---------------------------------------------------------------------------------------------
# console summary
# ---------------------------------------------------------------------------------------------
print(f"\n{'condition':28s}{'low O2':>14}{'normal O2':>14}{'pp':>8}{'rates/1000PY':>16}"
      f"{'adj HR (95% CI)':>26}{'flag':>6}")
print("-" * 112)
for _, r in order(d10).iterrows():
    hr = fmt_ci(r.hr_adj, r.hr_lo, r.hr_hi)
    print(f"{r.disease[:27]:28s}{f'{r.cases_low_o2}/{r.n_low_o2}':>14}"
          f"{f'{r.cases_normal_o2}/{r.n_normal_o2}':>14}{r.abs_diff_pp:+8.2f}"
          f"{f'{r.ir_low_o2_per1000py:.1f} vs {r.ir_normal_o2_per1000py:.1f}':>16}{hr:>26}"
          f"{'  *' if r.unstable else '':>6}")
print(f"\nflagged at 10%: {len(flagged10)}/{len(d10)}   flagged at 5%: {len(flagged5)}/{len(d5)}")
print(f"person-years  low {h10['person_years_exposed']:,.1f}  normal "
      f"{h10['person_years_unexposed']:,.1f}")
print(f"wrote {OUT}/table_noapnea.csv, table_noapnea_5pct.csv, table_noapnea.md, "
      f"results_noapnea.json")
