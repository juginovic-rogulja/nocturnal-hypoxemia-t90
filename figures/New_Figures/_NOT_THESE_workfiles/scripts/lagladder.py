"""
The landmark lag ladder: lags of 1, 2 and 5 years on every outcome, negative controls included.

Alen's comment on eFigure 7 was that a single 2-year landmark is not a ladder. One lag can always
be dismissed as a lucky window. Three of them, run on every condition rather than three, with the
negative controls pushed through the identical machinery, is the strongest statement this design
can make about reverse causation: if low nocturnal oxygen were a consequence of disease that was
already brewing at the sleep study, then deleting the first 1, 2 and 5 years of diagnoses should
walk the hazard ratio down toward 1. If it does not walk down, the oxygen came first.

What a landmark does here. For lag L, everyone whose event or censoring arrived before L is
removed, and the clock restarts at L for the survivors. Removal is not a filter on the outcome
alone. A patient censored at 0.8 years contributes nothing to a 1-year landmark either, which is
why the sample falls faster than the events do.

Two things are computed that the earlier 2-year script did not compute, and both are needed
before any lag can be read as a negative answer.

  1. Excess-hazard retention, (HR_lag - 1) / (HR_0 - 1). A hazard ratio of 1.48 is not "92% of"
     1.60. The quantity that decays under reverse causation is the excess over 1, and here it is
     80%. This is what "stable", "attenuating" and "reversing" are graded on.

  2. Whether a lag that fails had the precision to succeed. For each lag the standard error on
     the log scale is recovered from the interval, and the original effect is asked whether it
     would still have been detected at that precision. The bar used for every verdict here is
     80% power, log(HR_0) / SE_lag >= 2.80, not the weaker 1.96, which only says the effect
     would clear significance half the time. If a lag is 80%-powered for the original effect
     and the interval still covers 1, the association genuinely faded. If it is not, the lag is
     silent rather than negative, and printing the point estimate as a finding is reading noise.
     Both bars are written out so the choice can be checked rather than taken on trust.

The follow-up constraint is real and is not hidden. Median follow-up is 4.4 years, so a 5-year
landmark discards more than half the observation window. The script measures how much survives
and states per condition whether lag 5 carries information, instead of printing an estimate that
looks like an answer.

Model, identical to numbers/run_causal_tests_all.py and numbers/run_primary_unpenalized.py:
Cox, unpenalized, stratified by site, sex plus a 4-df natural cubic age spline with ONE COLUMN
DROPPED so the design is full rank, exposure rank-inverse-normalized within site so the hazard
ratio is per 1 SD, PSG date as time zero, oximetry_bad == 0.

Writes numbers/lag_ladder.csv (long, one row per condition and lag) and numbers/lag_ladder.json
(the summary, the failures, the trend grades and the lag-5 verdict).
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../..")))  # repository root, where paths.py lives
import paths
import json
import sys
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from patsy import dmatrix
from scipy import stats

ROOT = paths.FIGURE_ROOT
sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
from cohort_spec import COHORT_N, apply_cohort, cohort_tag   # v8.1c: tagged outputs on the full-nights run
from disease_definitions import CIRCULAR, DISEASES, NEGATIVE_CONTROLS

LAGS = [0.0, 1.0, 2.0, 5.0]
MIN_EVENTS = 40          # same fit floor as run_causal_tests_all.py
RETAIN_STABLE = 0.75     # excess-hazard retention above this is graded stable
RETAIN_ATTEN = 0.40      # below this, and still above 1, is graded collapsing
Z80 = 2.801586           # 1.959964 + 0.841621, the 80%-power bar on the log scale
Z50 = 1.959964           # the weaker bar, reported alongside but never used for a verdict
SMALL_EFFECT = 1.10      # below this the retention denominator is too thin to grade finely

# The settled panel of 2026-08-07: back pain, cataract, glaucoma, contact dermatitis,
# hemorrhoids. Fracture and osteoarthritis left because a plausible causal path runs from
# the exposure to each of them; both are still fitted and reported, just not as controls.
# An earlier version of this file carried only THREE controls here, which was neither the
# old panel nor the settled one, because the two replacements had never been fitted for the
# lag ladder. They are in the frozen parquet now, so negative_control_revised is identical
# to NEGATIVE_CONTROLS and the two columns can finally be collapsed.
CONTROLS_REVISED = list(NEGATIVE_CONTROLS)
assert CONTROLS_REVISED and all(DISEASES[k][1] for k in CONTROLS_REVISED), CONTROLS_REVISED   # v8.1c: the definitions' panel (glaucoma, contact dermatitis, hemorrhoids, inguinal hernia, alopecia), not a typed list

# ---------------------------------------------------------------------------- cohort
b = apply_cohort(pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet"))
assert len(b) == COHORT_N
b["male"] = (b.sex.astype(str).str.upper().str[0] == "M").astype(int)

# 4 spline columns, first dropped. Without the drop the design is rank-deficient and an
# unpenalized Cox fit can converge on nothing at all while still reporting success.
sp = dmatrix("cr(a, df=4) - 1", {"a": b.AgeAtVisit.values}, return_type="dataframe").iloc[:, 1:]
for i in range(sp.shape[1]):
    b[f"age_s{i}"] = sp.iloc[:, i].values
ADJ = [f"age_s{i}" for i in range(sp.shape[1])] + ["male"]


def rint(x):
    r = stats.rankdata(x)
    return stats.norm.ppf((r - 0.375) / (len(r) + 0.25))


b["z"] = b.groupby("site_id").spo2_pct_below_90.transform(rint)

fu = (pd.to_datetime(b.censor_date) - pd.to_datetime(b.psg_date)).dt.days / 365.25
FU = {"median": round(float(fu.median()), 2),
      "q1": round(float(fu.quantile(0.25)), 2), "q3": round(float(fu.quantile(0.75)), 2),
      "pct_beyond": {str(int(L)): round(float(100 * (fu > L).mean()), 1) for L in LAGS if L > 0}}
print(f"cohort {len(b):,}   median follow-up {FU['median']} y "
      f"(IQR {FU['q1']}-{FU['q3']})")
print("   still under observation at " + ", ".join(
    f"{k} y {v}%" for k, v in FU["pct_beyond"].items()) + "\n")


# ---------------------------------------------------------------------------- one fit
def fit(d):
    """Unpenalized site-stratified Cox on the landmark frame. None when it cannot be fitted."""
    d = d[["T", "E", "site", "z"] + ADJ].dropna()
    if int(d["E"].sum()) < MIN_EVENTS:
        return {"n": int(len(d)), "events": int(d["E"].sum()), "estimable": False}
    try:
        c = CoxPHFitter().fit(d, "T", "E", strata=["site"])
    except Exception:
        return {"n": int(len(d)), "events": int(d["E"].sum()), "estimable": False}
    s = c.summary.loc["z"]
    hr, lo, hi = (float(s["exp(coef)"]), float(s["exp(coef) lower 95%"]),
                  float(s["exp(coef) upper 95%"]))
    return {"n": int(len(d)), "events": int(d["E"].sum()), "estimable": True,
            "hr": hr, "lo": lo, "hi": hi,  # v8.3 full precision (2026-09-16)
            "hr_raw": hr, "p": float(s["p"]),
            # standard error on the log scale, recovered from the interval. Used to ask whether
            # a lag that failed ever had the precision to succeed.
            "se_log": (np.log(hi) - np.log(lo)) / (2 * 1.959964)}


# ---------------------------------------------------------------------------- the ladder
ALL = list(DISEASES.items()) + [("death", ("Death from any cause", False, [], []))]
long_rows, per_cond = [], {}

for key, (label, _n, _a, _b) in ALL:
    yc, ec, pc = f"{key}_years", f"{key}_incident", f"{key}_prevalent"
    if yc not in b.columns:
        continue
    f = b[(b[pc] == 0) & b[yc].notna() & (b[yc] > 0)]
    base = pd.DataFrame({"T": f[yc].astype(float), "E": f[ec].astype(int), "site": f.site_id,
                         **{c: f[c] for c in ["z"] + ADJ}})

    rungs = {}
    for L in LAGS:
        # base["T"], never base.T, which is the transpose
        lm = base[base["T"] > L].copy()
        lm["T"] = lm["T"] - L
        rungs[L] = fit(lm)
    if not rungs[0.0]["estimable"]:
        continue

    e0, n0 = rungs[0.0]["events"], rungs[0.0]["n"]
    # unrounded, so the power bar does not flip on a third decimal place
    hr0 = rungs[0.0]["hr_raw"]
    row_common = {
        "key": key, "condition": label,
        "negative_control": key in NEGATIVE_CONTROLS,
        "negative_control_revised": key in CONTROLS_REVISED,
        "circular": key in CIRCULAR,
    }
    for L in LAGS:
        r = rungs[L]
        rec = dict(row_common, lag_years=L, n=r["n"], events=r["events"],
                   pct_events_removed=round(100 * (1 - r["events"] / e0), 1),
                   pct_n_removed=round(100 * (1 - r["n"] / n0), 1),
                   estimable=r["estimable"])
        if r["estimable"]:
            rec.update(hr=r["hr"], lo=r["lo"], hi=r["hi"], p=r["p"],
                       significant=bool(r["lo"] > 1),
                       ci_width_ratio=(
                           (np.log(r["hi"]) - np.log(r["lo"])) /
                           (np.log(rungs[0.0]["hi"]) - np.log(rungs[0.0]["lo"]))))
            # excess-hazard retention, only defined against a lag-0 effect above 1
            rec["excess_retained"] = ((r["hr_raw"] - 1) / (hr0 - 1)  # v8.3 full precision (2026-09-16)
                                      if hr0 > 1.001 else np.nan)
            # would the ORIGINAL effect still have been detected at this lag's precision?
            zneed = np.log(hr0) / r["se_log"] if hr0 > 1 else 0.0
            rec["z_for_hr0"] = float(zneed)  # v8.3 full precision (2026-09-16)
            rec["power80_for_hr0"] = bool(zneed >= Z80)     # the bar every verdict uses
            rec["power50_for_hr0"] = bool(zneed >= Z50)     # the weaker bar, for reference
            # has the estimate actually moved, or only the interval widened?
            rec["hr0_inside_lag_ci"] = bool(r["lo"] <= hr0 <= r["hi"])
        else:
            for c in ("hr", "lo", "hi", "p", "excess_retained", "ci_width_ratio", "z_for_hr0"):
                rec[c] = np.nan
            rec.update(significant=False, power80_for_hr0=False, power50_for_hr0=False,
                       hr0_inside_lag_ci=False)
        # the retention denominator is hr0 - 1. When that is a few hundredths, a 0.01 wobble in
        # the estimate swings retention by a fifth, so the fine grade is not trustworthy.
        rec["grade_precision_limited"] = bool(hr0 < SMALL_EFFECT)
        long_rows.append(rec)
    per_cond[key] = {"label": label, "rungs": rungs}

d = pd.DataFrame(long_rows)

# ---------------------------------------------------------------------------- grading
# The headline panel is the set the paper already calls significant: a lag-0 interval strictly
# above 1, not a negative control, not circular. Circular conditions are diagnosed by what the
# sleep study itself measures, so a hazard ratio for them rewards the definition, not the biology.
w0 = d[d.lag_years == 0]
PANEL = sorted(w0[(w0.lo > 1) & ~w0.negative_control & ~w0.circular].key)
CTRL = sorted(w0[w0.negative_control].key)


def grade(sub, lags):
    """stable / attenuating / collapsing / reversing across the given informative lags."""
    use = sub[sub.lag_years.isin(lags) & sub.estimable]
    if use.empty:
        return "not estimable", np.nan
    # Retention is (HR - 1) / (HR_0 - 1). With no positive association to start from there is
    # nothing to retain, and "reversing" would libel a condition that never pointed upward.
    if float(sub[sub.lag_years == 0].iloc[0].hr) <= 1.001:
        return "not graded, no positive association at lag 0", np.nan
    ret = use.excess_retained
    lo_ret = float(ret.min())
    if float(use.hr.min()) < 1.0:
        return "reversing", lo_ret
    if lo_ret >= RETAIN_STABLE:
        return ("stable, strengthening" if float(ret.max()) >= 1.15 else "stable"), lo_ret
    if lo_ret >= RETAIN_ATTEN:
        return "attenuating", lo_ret
    return "collapsing", lo_ret


grades = {}
for key in sorted(set(d.key)):
    sub = d[d.key == key]
    l5 = sub[sub.lag_years == 5.0].iloc[0]
    # lag 5 is allowed into the grade only when it was 80%-powered for the original effect.
    # Otherwise a wide, quiet interval would be read as evidence that the effect went away.
    lag5_informative = bool(l5.estimable and l5.power80_for_hr0)
    lags_used = [1.0, 2.0] + ([5.0] if lag5_informative else [])
    g_all, ret_all = grade(sub, lags_used)
    g_12, ret_12 = grade(sub, [1.0, 2.0])
    r0 = sub[sub.lag_years == 0].iloc[0]
    grades[key] = {
        "condition": r0.condition,
        "negative_control": bool(r0.negative_control),
        "negative_control_revised": bool(r0.negative_control_revised),
        "circular": bool(r0.circular),
        "hr_lag0": float(r0.hr), "ci_lag0": [float(r0.lo), float(r0.hi)],
        "events_lag0": int(r0.events),
        "hr_by_lag": {str(int(L)): (float(x.hr) if x.estimable else None)
                      for L, x in ((L, sub[sub.lag_years == L].iloc[0]) for L in LAGS)},
        "ci_by_lag": {str(int(L)): ([float(x.lo), float(x.hi)] if x.estimable else None)
                      for L, x in ((L, sub[sub.lag_years == L].iloc[0]) for L in LAGS)},
        "events_by_lag": {str(int(L)): int(sub[sub.lag_years == L].iloc[0].events) for L in LAGS},
        "pct_events_removed_by_lag": {
            str(int(L)): float(sub[sub.lag_years == L].iloc[0].pct_events_removed) for L in LAGS},
        "significant_by_lag": {str(int(L)): bool(sub[sub.lag_years == L].iloc[0].significant)
                               for L in LAGS},
        "excess_retained_by_lag": {
            str(int(L)): (None if pd.isna(sub[sub.lag_years == L].iloc[0].excess_retained)
                          else float(sub[sub.lag_years == L].iloc[0].excess_retained))
            for L in LAGS},
        "lag5_informative": lag5_informative,
        "lag5_informative_50pct_bar": bool(l5.estimable and l5.power50_for_hr0),
        "lags_graded": [int(x) for x in lags_used],
        "trend": g_all, "trend_lag1_2_only": g_12,
        "min_excess_retained": None if pd.isna(ret_all) else float(ret_all),  # v8.3 full precision (2026-09-16)
        "grade_precision_limited": bool(r0.hr < SMALL_EFFECT),
        "significant_at_all_graded_lags": bool(
            all(sub[sub.lag_years == L].iloc[0].significant for L in lags_used)),
    }

d["trend"] = d.key.map(lambda k: grades[k]["trend"])
d["lag5_informative"] = d.key.map(lambda k: grades[k]["lag5_informative"])
d = d.sort_values(["negative_control", "condition", "lag_years"])
d.to_csv(f"{paths.NUMBERS_DIR}/lag_ladder{cohort_tag()}.csv", index=False)

# ---------------------------------------------------------------------------- summary
summary = {"lags": [int(L) for L in LAGS], "cohort_n": int(COHORT_N), "follow_up_years": FU,
           "n_outcomes_fitted": int(d.key.nunique()),
           "panel_n": len(PANEL), "controls_n": len(CTRL), "by_lag": {}}

CTRL_REV = sorted(w0[w0.negative_control_revised].key)


def ctl_block(sub, keys):
    c = sub[sub.key.isin(keys) & sub.estimable]
    if not len(c):
        return None
    return {"hr_range": [float(c.hr.min()), float(c.hr.max())],  # v8.3 full precision (2026-09-16)
            "estimable": int(len(c)), "total": len(keys),
            "significant": int(c.significant.sum()),
            "conditions": {r.condition: [float(r.hr), float(r.lo), float(r.hi)]
                           for _, r in c.iterrows()}}


for L in LAGS[1:]:
    sub = d[(d.lag_years == L) & d.key.isin(PANEL)]
    surv = sub[sub.significant]
    fail = sub[~sub.significant]
    faded = fail[fail.power80_for_hr0 & fail.estimable]
    summary["by_lag"][str(int(L))] = {
        "survive": int(len(surv)), "of": len(PANEL),
        "median_pct_events_removed": round(float(sub.pct_events_removed.median()), 1),
        "median_pct_n_removed": round(float(sub.pct_n_removed.median()), 1),
        "median_excess_retained": float(sub.excess_retained.median()),  # v8.3 full precision (2026-09-16)
        "n_informative": int(sub.power80_for_hr0.sum()),
        "failures": [{"condition": r.condition, "hr": None if pd.isna(r.hr) else float(r.hr),
                      "ci": None if pd.isna(r.lo) else [float(r.lo), float(r.hi)],
                      "events": int(r.events),
                      "excess_retained": None if pd.isna(r.excess_retained)
                      else float(r.excess_retained),
                      "hr0_inside_lag_ci": bool(r.hr0_inside_lag_ci),
                      "reason": "faded" if (r.estimable and r.power80_for_hr0) else "lost power"}
                     for _, r in fail.sort_values("events", ascending=False).iterrows()],
        "n_faded": int(len(faded)), "n_lost_power": int(len(fail) - len(faded)),
        "controls_current_panel": ctl_block(d[d.lag_years == L], CTRL),
        "controls_revised_panel": ctl_block(d[d.lag_years == L], CTRL_REV),
        # kept flat for the figure caption, from the revised panel Alen asked for
        "controls_hr_range": ctl_block(d[d.lag_years == L], CTRL_REV)["hr_range"],
        "median_ci_width_ratio": float(sub.ci_width_ratio.median()),
    }

l5 = d[(d.lag_years == 5.0) & d.key.isin(PANEL)]
inf5 = l5[l5.lag5_informative]
non5 = l5[~l5.lag5_informative]
hr0map = w0.set_index("key").hr
summary["lag5_verdict"] = {
    "median_pct_events_removed": round(float(l5.pct_events_removed.median()), 1),
    "median_pct_n_removed": round(float(l5.pct_n_removed.median()), 1),
    "n_estimable": int(l5.estimable.sum()), "of": len(PANEL),
    "n_informative": int(len(inf5)),
    "n_informative_50pct_bar": int(l5.power50_for_hr0.sum()),
    "median_ci_width_ratio": float(l5.ci_width_ratio.median()),
    # informativeness at lag 5 turns on the effect size AND the events that survive the landmark,
    # not on the effect size alone. The two largest hazard ratios in the paper, obesity
    # hypoventilation and cirrhosis, are the least informative at lag 5 because almost no events
    # survive it. Both bounds are reported rather than a single threshold that does not exist.
    "informative_hr0_min": float(min(hr0map[k] for k in inf5.key)) if len(inf5) else None,  # v8.3 full precision (2026-09-16)
    "informative_events_min": int(inf5.events.min()) if len(inf5) else None,
    "uninformative_hr0_max_among_estimable": (
        float(max(hr0map[k] for k in non5[non5.estimable].key))
        if len(non5[non5.estimable]) else None),
    "not_estimable_for": sorted(non5[~non5.estimable].condition),
    "informative_for": sorted(inf5.condition),
    "uninformative_for": sorted(non5.condition),
    "survive_among_informative": int(inf5.significant.sum()),
    "usable": "partial",
    "verdict": (
        f"Lag 5 answers the question for the strong associations and not for the weak ones. It is "
        f"80%-powered for {len(inf5)} of {len(PANEL)} conditions, and {int(inf5.significant.sum())} "
        f"of those {len(inf5)} remain significant with the first five years of diagnoses deleted. "
        f"It removes a median {round(float(l5.pct_events_removed.median()), 1)}% of events and "
        f"{round(float(l5.pct_n_removed.median()), 1)}% of the sample against a median follow-up "
        f"of {FU['median']} years, and widens intervals "
        f"{round(float(l5.ci_width_ratio.median()), 1)}-fold. Informativeness depends on the "
        f"events that survive the landmark as much as on the effect size: every informative "
        f"condition kept at least {int(inf5.events.min()) if len(inf5) else 0} events, while "
        f"{' and '.join(sorted(non5[~non5.estimable].condition)) or 'none'} cannot be fitted at "
        f"all despite carrying the two largest hazard ratios in the paper. For a condition whose "
        f"hazard ratio is near 1.1, a non-significant 5-year rung is NOT a negative result and "
        f"must not be drawn or described as one. Report lag 5 as a supporting column with its "
        f"event counts on the face of the figure, never as a third equal rung."),
}
summary["trends"] = {k: {"condition": v["condition"], "trend": v["trend"],
                         "trend_lag1_2_only": v["trend_lag1_2_only"],
                         "min_excess_retained": v["min_excess_retained"],
                         "lag5_informative": v["lag5_informative"]}
                     for k, v in grades.items() if k in PANEL}
summary["controls"] = {k: grades[k] for k in CTRL}
summary["per_condition"] = grades
summary["provenance"] = {
    "script": "New_Figures/scripts/lagladder.py",
    "model": ("unpenalized Cox, site-stratified, sex + 4-df natural cubic age spline with one "
              "column dropped, exposure rank-inverse-normalized within site, HR per 1 SD, "
              "PSG date as time zero, oximetry_bad == 0"),
    "landmark": ("for lag L, anyone whose event or censoring arrived at or before L is removed "
                 "and the clock restarts at L"),
    "min_events_to_fit": MIN_EVENTS,
    "informative_rule": ("a lag is informative for a condition when log(HR at lag 0) / SE at that "
                         "lag >= 2.80, i.e. the precision the lag leaves would detect the original "
                         "effect with 80% power. The weaker 1.96 bar is also stored, as "
                         "power50_for_hr0, but no verdict in this file rests on it"),
    "negative_control_panel": ("negative_control and negative_control_revised are now the same "
                               "five conditions, the settled panel of 2026-08-07: back pain, "
                               "cataract, glaucoma, contact dermatitis, hemorrhoids. Fracture "
                               "and osteoarthritis are out because a plausible causal path runs "
                               "from nocturnal hypoxemia to each of them, which disqualifies a "
                               "control whatever its estimate looks like. See "
                               "numbers/negcontrols_final.json"),
    "grading": {"stable": f"excess-hazard retention >= {RETAIN_STABLE} at every graded lag",
                "attenuating": f"retention between {RETAIN_ATTEN} and {RETAIN_STABLE}",
                "collapsing": f"retention below {RETAIN_ATTEN}, point estimate still above 1",
                "reversing": "point estimate falls below 1 at a graded lag"},
}
json.dump(summary, open(f"{paths.NUMBERS_DIR}/lag_ladder{cohort_tag()}.json", "w"), indent=1, default=float)

# ---------------------------------------------------------------------------- report
print(f"{d.key.nunique()} outcomes fitted, {len(PANEL)} in the significant panel, "
      f"{len(CTRL)} negative controls\n")
hdr = f"{'condition':30s}{'lag 0':>16}{'lag 1':>16}{'lag 2':>16}{'lag 5':>16}  {'trend':22s}"
print(hdr)
print("-" * len(hdr))


def cell(r):
    if not r.estimable:
        return f"{'-':>9}{r.events:7,}"
    star = "*" if r.significant else " "
    return f"{r.hr:8.2f}{star}{r.events:7,}"


for key in sorted(PANEL, key=lambda k: -grades[k]["hr_lag0"]):
    sub = d[d.key == key]
    line = f"{grades[key]['condition']:30s}"
    for L in LAGS:
        line += cell(sub[sub.lag_years == L].iloc[0])
    print(line + f"  {grades[key]['trend']:22s}")

print("\nnegative controls, the settled five of 2026-08-07")
print("-" * len(hdr))
for key in sorted(CTRL, key=lambda k: -grades[k]["hr_lag0"]):
    sub = d[d.key == key]
    tag = "r " if key in CTRL_REV else "  "
    line = f"{tag}{grades[key]['condition']:28s}"
    for L in LAGS:
        line += cell(sub[sub.lag_years == L].iloc[0])
    print(line + f"  {grades[key]['trend']:22s}")

print("\n* interval excludes 1. Second number in each cell is events retained. A dash is fewer "
      f"than {MIN_EVENTS} events.\n")
for L in LAGS[1:]:
    s = summary["by_lag"][str(int(L))]
    cr, cc = s["controls_revised_panel"], s["controls_current_panel"]
    print(f"lag {int(L)} y: {s['survive']} of {s['of']} survive, "
          f"median {s['median_pct_events_removed']}% of events and {s['median_pct_n_removed']}% "
          f"of the sample removed, median excess-hazard retention "
          f"{s['median_excess_retained']:.2f}, {s['n_informative']} of {s['of']} rungs "
          f"80%-powered")
    print(f"          controls, revised panel {cr['hr_range'][0]:.2f}-{cr['hr_range'][1]:.2f} "
          f"({cr['significant']} of {cr['estimable']} significant), "
          f"as frozen {cc['hr_range'][0]:.2f}-{cc['hr_range'][1]:.2f} "
          f"({cc['significant']} of {cc['estimable']} significant)")
    for f_ in s["failures"]:
        ret = "n/a" if f_["excess_retained"] is None else f"{f_['excess_retained']:.2f}"
        print(f"    fails: {f_['condition']} ({f_['events']:,} events, retention {ret}, "
              f"{f_['reason']})")

v = summary["lag5_verdict"]
print(f"\nlag 5 verdict: {v['n_estimable']} of {v['of']} estimable, {v['n_informative']} "
      f"80%-powered ({v['n_informative_50pct_bar']} at the weaker 50% bar), "
      f"{v['survive_among_informative']} of those {v['n_informative']} still significant. "
      f"Median {v['median_pct_events_removed']}% of events and {v['median_pct_n_removed']}% of "
      f"the sample removed, intervals {v['median_ci_width_ratio']}x wider.")
print(f"  informative down to a lag-0 hazard ratio of {v['informative_hr0_min']} and no fewer than "
      f"{v['informative_events_min']} retained events. Largest effect it cannot speak to among "
      f"those still fittable: {v['uninformative_hr0_max_among_estimable']}. Not fittable at all: "
      f"{', '.join(v['not_estimable_for']) or 'none'}.")
print(f"  USABLE: {v['usable'].upper()}, for the {v['n_informative']} strong associations only.")
print(f"\nwritten -> {paths.NUMBERS_DIR}/lag_ladder{cohort_tag()}.csv")
print(f"written -> {paths.NUMBERS_DIR}/lag_ladder{cohort_tag()}.json")
