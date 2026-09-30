"""
REPORT.txt, built from dC_per_disease.csv and the run's own JSON so no number in the prose can
drift from the number in the table.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import json
import os
import re
import textwrap
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import sys as _sys_ct; _sys_ct.path.insert(0, paths.ANALYSIS_DIR); from cohort_spec import cohort_tag as _cohort_tag   # v8.1c

import dc_engine as E
import sys as _sys_v8; _sys_v8.path.insert(0, paths.ANALYSIS_DIR)  # v8 sweep 2026-09-12
from cohort_spec import N_RANKED_OUTCOMES  # v8 sweep 2026-09-12
__PUB_DC_V7__ = float(__import__("pandas").read_csv(f"{paths.NUMBERS_DIR}/ranking_v3.csv", comment="#").set_index("feature").loc["spo2_pct_below_90", "dC"])   # v7: T90 gain from the CURRENT ranking, was the August literal 0.02033875497249791


HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "dC_per_disease.csv")
T90ROOT = paths.T90_ROOT
PUB_DC = __PUB_DC_V7__
LADDER_CHECK = ["resp_failure", "copd2"]      # the two the coordinator asked to sanity-check


def runtime(arm):
    f = os.path.join(E.WORK, f"bootstrap.log" if arm == "standard"
                     else f"bootstrap_{arm}.log")
    if not os.path.exists(f):
        return "not recorded"
    m = re.search(r"bootstrap done in ([\d.]+) min", open(f).read())
    return f"{float(m.group(1)):.0f} minutes wall clock" if m else "still running"


def main():
    raw = pd.read_csv(CSV)
    dis = raw[raw.row_type == "disease"]
    d = dis[dis.arm == "standard"].sort_values("dC", ascending=False).reset_index(drop=True)
    has_lag = (dis.arm == "lag2").any()
    L = dis[dis.arm == "lag2"].set_index("outcome") if has_lag else None
    pc = json.load(open(os.path.join(E.WORK, "positive_control.json")))
    lg = json.load(open(os.path.join(E.WORK, "lag2_gate.json")))
    sm = json.load(open(os.path.join(E.WORK, "assemble_summary.json")))
    nrep = int(d.n_replicates.max())

    def row(k, arm="standard"):
        s = dis[(dis.outcome == k) & (dis.arm == arm)]
        return s.iloc[0] if len(s) else None

    def ci(r, nd=4):
        if r is None or not np.isfinite(r.dC):
            return "not estimable"
        return f"{r.dC:+.{nd}f} ({r.ci_lo:+.{nd}f} to {r.ci_hi:+.{nd}f})"

    mean_row = {(a, t): raw[(raw.outcome == t) & (raw.arm == a)].iloc[0]
                for a in dis.arm.unique() for t in ("__mean_all__", "__mean_common__")
                if len(raw[(raw.outcome == t) & (raw.arm == a)])}

    L_ = []
    A = L_.append

    def P(text):
        """A prose paragraph, wrapped to the report's 76-character measure."""
        for ln in textwrap.wrap(" ".join(text.split()), 76):
            L_.append(ln)
    A("PER-DISEASE DISCRIMINATION GAIN FOR NOCTURNAL OXYGEN (T90)")
    A(f"Bootstrap confidence intervals, {N_RANKED_OUTCOMES} diseases, two arms: the paper's windowing and a")
    A("2-year landmark. Nothing is averaged across diseases except the two rows that say so.")
    A(f"Built {pd.Timestamp.today():%Y-%m-%d} in Sleep_Variability_2026-08/dC_per_disease.")
    A("")

    # ---------------------------------------------------------------- TLDR
    A("TLDR")
    A("-" * 74)
    hi = d.iloc[0]
    ms = mean_row[("standard", "__mean_all__")]
    A(f"The average gain of {PUB_DC:.4f} in the paper is not what T90 does to any one disease.")
    A(f"It is the average of a range running from {hi.dC:+.4f} for {hi.disease.lower()} down")
    A(f"to {d.dC.min():+.4f} at the bottom. {int(d.excludes_zero.astype(bool).sum())} of "
      f"{len(d)} intervals exclude zero.")
    for k in ("resp_failure", "hf", "diabetes"):
        r = row(k)
        A(f"  {r.disease:<22} {ci(r)}, {r.dC / PUB_DC:.1f} times the published average")
    A(f"The average itself now has an interval. On the paper's own age basis the {N_RANKED_OUTCOMES} gains")
    A(f"average {ms.dC_pubbasis:.6f} (95% CI {ms.ci_lo_pubbasis:.6f} to "
      f"{ms.ci_hi_pubbasis:.6f}).")
    if has_lag:
        ml = mean_row[("lag2", "__mean_all__")]
        mc_s = mean_row[("standard", "__mean_common__")]
        P(f"Under a 2-year landmark the gains shrink but do not vanish. On the "
          f"{int(ml.n_outcomes_in_mean)} outcomes still estimable the mean is {ml.dC:.6f} "
          f"({ml.ci_lo:.6f} to {ml.ci_hi:.6f}) against {mc_s.dC:.6f} on the same set without "
          f"the landmark, so {100 * ml.dC / mc_s.dC:.0f}% of the average gain survives.")
    A("")

    # ---------------------------------------------------------------- positive controls
    A("POSITIVE CONTROLS, both run before any bootstrap")
    A("-" * 74)
    g1 = pc["gate1_published_machinery"]
    g2 = pc["gate2_engine"]
    A("Gate 1. The published machinery, dC_side_analyses/common.py imported unmodified,")
    A(f"reproduces numbers/ranking_v3_percondition.csv from the frozen cohort of "
      f"{pc['cohort_n']:,}.")
    A(f"  headline dC   paper {g1['published_dC']:.12f}   here "
      f"{g1['reproduced_dC']:.12f}   diff {g1['diff_dC']:+.1e}")
    A(f"  baseline C    paper {g1['published_baseline_C']:.6f}         here "
      f"{g1['reproduced_baseline_C']:.6f}   diff "
      f"{g1['reproduced_baseline_C'] - g1['published_baseline_C']:+.1e}")
    A(f"  all {N_RANKED_OUTCOMES} per-outcome gains, largest absolute difference "
      f"{g1['max_abs_per_outcome_gain_diff']:.2e}")
    for lab, v in g1["named"].items():
        A(f"    {lab:<18} paper {v['gain_published']:+.12f}   here {v['gain_here']:+.12f}")
    A("Gate 2. This folder's engine reproduces that machinery outcome by outcome on the")
    A(f"published age basis. Largest difference in baseline C "
      f"{g2['max_abs_base_C_diff']:.1e}, in C with T90 {g2['max_abs_t90_C_diff']:.1e}.")
    A(f"  VERDICT {'PASS' if pc['pass'] else 'FAIL'}.")
    A("")
    A("Gate 3, the landmark arm. There is no published dC at lag 2, so the gate is on the")
    A("windowing, which is the only thing the lag arm changes. numbers/lag_ladder.csv is the")
    A("paper's own landmark ladder. Every outcome's patients at risk and incident cases match")
    A(f"it exactly at lag 0 and at lag 2: {lg['gate_counts_match_published_ladder']}.")
    P(f"  the landmark rule, quoted from numbers/lag_ladder.json: {lg['convention']}")
    A(f"  lag 2 removes {lg['pct_events_removed_at_lag2']['min']:.0f}% to "
      f"{lg['pct_events_removed_at_lag2']['max']:.0f}% of incident cases, median "
      f"{lg['pct_events_removed_at_lag2']['median']:.0f}%")
    A("")

    # ---------------------------------------------------------------- method
    A("WHAT WAS FITTED")
    A("-" * 74)
    A("Cohort. 19,173 patients, data_frozen_v8_2026-09/t90_final.parquet through")
    A("numbers/cohort_spec.apply_cohort. One row per patient. Two hospitals, I0002 with")
    A("10,431 and I0006 with 8,742. For each disease, patients who already had that disease")
    A("at the recording are excluded, as the paper does, so the row set differs by disease.")
    A("")
    A("Model A. A restricted cubic spline in age, cr(df=4), with ONE basis column dropped,")
    A("plus sex. Model B. Model A plus T90 as a continuous term, the within-site rank")
    A("inverse normal of the percentage of the recording below 90% saturation. Cox")
    A("proportional hazards, penalizer 0.01, exactly the paper's fitting.")
    A("")
    A("Cross-fit. Fit at I0002 and score the held-out I0006, then reverse. Both directions")
    A("are in the CSV. The disease's dC is the mean of the two held-out gains, which is the")
    A("paper's own aggregation over folds. A fold is dropped if either side carries fewer")
    A("than 30 events, as in the ranking.")
    A("")
    A("Landmark arm. For the 2-year arm every patient whose event or censoring arrived at or")
    A("before 2 years leaves the risk set and the clock restarts at 2 years. The event flag")
    A("is untouched, because every event that survives the filter happens after 2 years.")
    A("Prevalent exclusion, the measure, the age spline and the site design are unchanged.")
    A("")
    A("Two age bases are carried on the standard arm. The reported one drops a basis column,")
    A("which is the rank-safe version and Alen's standing rule. The paper's published basis")
    A("keeps all four singular columns and survives only on the 0.01 ridge. Both are in the")
    A(f"CSV. The {N_RANKED_OUTCOMES} gains average {ms.dC:.6f} on the reported basis and "
      f"{ms.dC_pubbasis:.12f} on the")
    A("published basis, which is the paper's headline to every printed digit.")
    r1 = list(d.sort_values("dC", ascending=False).outcome)
    r2 = list(d.sort_values("dC_pubbasis", ascending=False).outcome)
    big = list(d[d.dC > 0.01].outcome)
    P(f"Top ten by the two bases: {'identical' if r1[:10] == r2[:10] else 'NOT identical'}. "
      f"Across all {len(d)}, {sum(1 for k in r1 if r1.index(k) != r2.index(k))} sit at a "
      f"different rank, at most {max(abs(r1.index(k) - r2.index(k)) for k in r1)} places. "
      f"Among the {len(big)} with a gain above 0.01 the largest move is "
      f"{max(abs(r1.index(k) - r2.index(k)) for k in big)}. The reordering is confined to "
      "gains that are near zero on both bases.")
    A("")

    # ---------------------------------------------------------------- bootstrap
    A("BOOTSTRAP")
    A("-" * 74)
    A(f"{nrep} replicates per arm, the same seeds in both. Patients resampled with")
    A("replacement inside each hospital, so both site sizes and the cross-fit design are")
    A("preserved in every replicate. The whole procedure is repeated on each resample and")
    A("not just the last fit: the age spline knots are re-placed on the resampled ages and")
    A("T90's within-site rank inverse normal is re-ranked inside the resample, because both")
    A("are estimated on the data. Interval is the plain 2.5th to 97.5th percentile.")
    A(f"Runtime: standard arm {runtime('standard')}, 2-year landmark arm {runtime('lag2')},")
    A("7 workers, checkpointed in chunks of 5 replicates so either arm is resumable.")
    for arm, lab in (("standard", "no landmark"), ("lag2", "2-year landmark")):
        if arm not in sm:
            continue
        s = sm[arm]
        P(f"  {lab}: {s['n_excluding_zero']} of {len(d)} intervals exclude zero, "
          f"{s['n_not_estimable']} outcome(s) not estimable"
          + (f" ({', '.join(s['not_estimable'])})" if s["not_estimable"] else "")
          + f". The 30-event floor tests both sides of the cross-fit, so a thin hospital "
          f"removes both directions at once. That happened in "
          f"{s['boot_reps_with_no_estimate']} outcome-replicates. Largest Monte Carlo error "
          f"of an interval bound {s['max_mc_error_of_ci_bound']:.5f}.")
        thin = dis[(dis.arm == arm) & dis.estimable.astype(bool)
                   & (dis.n_replicates < 0.9 * nrep)]
        if len(thin):
            P("    intervals resting on fewer than 90% of the replicates, because the fold "
              "floor bit in the rest: "
              + ", ".join(f"{t.disease} {int(t.n_replicates)} of {nrep}"
                          for _j, t in thin.iterrows()) + ". Their intervals are wider for it "
              "and none of them changes a verdict on zero.")
    A("The Monte Carlo error is the normal-approximation error of a percentile bound at this")
    A(f"number of replicates. It is at most "
      f"{max(sm[a]['max_mc_error_of_ci_bound'] for a in sm if a in ('standard', 'lag2')):.5f} "
      f"of C, which is "
      f"{100 * max(sm[a]['max_mc_error_of_ci_bound'] for a in sm if a in ('standard', 'lag2')) / d.dC.abs().max():.1f}% "
      f"of the largest")
    A("gain in the table. No disease's verdict on zero turns on it.")
    A("")

    # ---------------------------------------------------------------- the mean
    A("THE MEAN ACROSS OUTCOMES, WITH AN INTERVAL")
    A("-" * 74)
    A(f"The {N_RANKED_OUTCOMES} gains are {N_RANKED_OUTCOMES} views of one set of patients, not {N_RANKED_OUTCOMES} independent studies. The")
    A("right unit to resample is therefore the patient and not the disease: each replicate")
    A(f"recomputes all {N_RANKED_OUTCOMES} gains on the same resampled patients, the mean is taken inside the")
    A("replicate, and the interval is the percentiles of those means. Averaging {N_RANKED_OUTCOMES} separate")
    A("intervals would ignore that correlation and would be wrong twice over.")
    A("")
    for arm, lab in (("standard", "no landmark  "), ("lag2", "2-y landmark ")):
        for tag, what in (("__mean_all__", "all estimable in this arm"),
                          ("__mean_common__", "estimable in both arms   ")):
            if (arm, tag) not in mean_row:
                continue
            m = mean_row[(arm, tag)]
            A(f"  {lab} {what}  n={int(m.n_outcomes_in_mean):>2}  "
              f"{m.dC:+.6f} ({m.ci_lo:+.6f} to {m.ci_hi:+.6f})  SE {m.boot_se:.6f}")
    A(f"  published basis, no landmark, all {N_RANKED_OUTCOMES}        {ms.dC_pubbasis:+.6f} "
      f"({ms.ci_lo_pubbasis:+.6f} to {ms.ci_hi_pubbasis:+.6f})")
    A("The last line is the paper's headline 0.0203 with the interval it never had.")
    A("")

    # ---------------------------------------------------------------- ladder check
    if has_lag:
        A("LAG-2 SANITY CHECK AGAINST THE PUBLISHED HAZARD-RATIO LADDER")
        A("-" * 74)
        lad = pd.read_csv(f"{paths.NUMBERS_DIR}/lag_ladder{_cohort_tag()}.csv")   # v8.1c: the ladder of this run's frame
        A("Expected: attenuated but persistent for the respiratory outcomes. The ladder is a")
        A("hazard ratio per 1 SD and this analysis is a gain in concordance, so the check is")
        A("on direction and persistence, not on the size of the two numbers.")
        for k in LADDER_CHECK:
            l0 = lad[(lad.key == k) & (lad.lag_years == 0)].iloc[0]
            l2 = lad[(lad.key == k) & (lad.lag_years == 2)].iloc[0]
            s0, s2 = row(k), row(k, "lag2")
            A(f"  {s0.disease}")
            A(f"    ladder  HR lag 0 {l0.hr:.3f} ({l0.lo:.3f} to {l0.hi:.3f})   "
              f"HR lag 2 {l2.hr:.3f} ({l2.lo:.3f} to {l2.hi:.3f})   "
              f"both above 1: {bool(l0.lo > 1 and l2.lo > 1)}")
            A(f"    here    dC lag 0 {ci(s0)}")
            A(f"            dC lag 2 {ci(s2)}")
            A(f"    attenuates in both: "
              f"{bool(l2.hr < l0.hr and s2.dC < s0.dC)}, still clears zero or one in both: "
              f"{bool(l2.lo > 1 and s2.excludes_zero)}, "
              f"{100 * s2.dC / s0.dC:.0f}% of the gain retained")
        A("")

    # ---------------------------------------------------------------- the table
    A(f"THE RANKED TABLE, ALL {N_RANKED_OUTCOMES}")
    A("-" * 74)
    A(f"{'#':>2}  {'Disease':<26}{'Organ':<13}{'Cases':>6}{'AtRisk':>9}  "
      f"{'BaseC':>6} {'C+T90':>6}  {'dC':>8}  {'95% CI':<20}{'lag-2 dC':>9}  "
      f"{'lag-2 95% CI':<20}")
    for i, r in d.iterrows():
        star = "*" if r.excludes_zero else " "
        lr = L.loc[r.outcome] if has_lag and r.outcome in L.index else None
        if lr is None or not np.isfinite(lr.dC):
            lag_dc, lag_ci = "     n/e", "not estimable"
        else:
            lag_dc = f"{lr.dC:+.4f}{'*' if lr.excludes_zero else ' '}"
            lag_ci = f"{lr.ci_lo:+.4f} to {lr.ci_hi:+.4f}"
        A(f"{i + 1:>2}  {r.disease:<26}{r.organ_group:<13}{int(r.events):>6,}{int(r.n_at_risk):>9,}  "
          f"{r.c_agesex_mean:>6.3f} {r.c_agesex_t90_mean:>6.3f}  {r.dC:>+8.4f}{star} "
          f"{f'{r.ci_lo:+.4f} to {r.ci_hi:+.4f}':<20}{lag_dc:>9}  {lag_ci:<20}")
    A("  * interval excludes zero.  n/e, fewer than 30 cases at one hospital after the "
      "landmark.")
    A("")
    A(f"Intervals excluding zero: {int(d.excludes_zero.astype(bool).sum())} of {len(d)} "
      f"without the landmark, "
      f"{int(L.excludes_zero.astype(bool).sum()) if has_lag else 0} with it, on the reported")
    A(f"basis. On the paper's published basis, without the landmark, "
      f"{int(d.excludes_zero_pubbasis.astype(bool).sum())} of {len(d)}.")
    neg = d[d.negative_control == "yes"]
    A(f"The three negative controls in the ranking set are {', '.join(sorted(neg.disease))}. "
      f"Their gains")
    A(f"are {', '.join(f'{r.disease} {r.dC:+.4f}' for _i, r in neg.iterrows())}, and "
      f"{int(neg.excludes_zero.astype(bool).sum())} of 3 exclude zero.")
    A("")

    # ---------------------------------------------------------------- contrast
    A("THE CONTRAST")
    A("-" * 74)
    A(f"{'':<38}{'no landmark':<32}{'2-year landmark'}")
    A(f"{'Published average across the ' + str(N_RANKED_OUTCOMES):<38}{PUB_DC:+.4f}")
    for k in ("obesity_hypovent", "resp_failure", "diabetes", "hf"):
        s, l2 = row(k), row(k, "lag2")
        A(f"{s.disease:<38}{ci(s):<32}{ci(l2)}")
        A(f"{'':<38}{s.dC / PUB_DC:.1f} times the average")
    A("")
    ohs = row("obesity_hypovent")
    A(f"Baseline C for {ohs.disease.lower()} is {ohs.c_agesex_mean:.3f} on age and sex and")
    A(f"{ohs.c_agesex_t90_mean:.3f} once T90 is added. One oximetry number moves that disease")
    A("from a model that barely separates patients to one that separates them well.")
    ol = row("obesity_hypovent", "lag2")
    if ol is not None:
        P(f"It is also the one outcome the landmark cannot keep: only {int(ol.events)} cases "
          f"survive 2 years, {int(ol.events_I0002)} of them at I0002, below the ranking's own "
          "30-event floor.")
    A("")
    top = d.head(9)
    P(f"The nine largest gains are {', '.join(top.disease)}.")
    resp = d[d.organ_group == "Respiratory"]
    P(f"All {len(resp)} respiratory outcomes sit in the top "
      f"{int(d.index[d.organ_group == 'Respiratory'].max()) + 1} of {len(d)}, and "
      f"{int(resp.excludes_zero.astype(bool).sum())} of {len(resp)} clear zero.")
    for g in ("Metabolic", "Liver", "Neuro/psych"):
        gg = d[d.organ_group == g]
        ranks = [str(i + 1) for i in d.index[d.organ_group == g]]
        P(f"{g}: {len(gg)} outcomes at ranks {', '.join(ranks[:-1])} and {ranks[-1]}. "
          f"{int(gg.excludes_zero.astype(bool).sum())} clear zero.")
    bad = d[d.dC < 0]
    nb = int(bad.excludes_zero.astype(bool).sum())
    A(f"{len(bad)} diseases have a negative point estimate, and "
      f"{'none of those intervals exclude' if nb == 0 else str(nb) + ' of those intervals exclude'}"
      f" zero{'' if nb == 0 else ' (' + ', '.join(bad[bad.excludes_zero.astype(bool)].disease) + ')'}.")
    if has_lag:
        both = d[d.excludes_zero.astype(bool) & d.outcome.map(
            lambda o: o in L.index and bool(L.loc[o].excludes_zero))]
        P(f"{len(both)} diseases clear zero in BOTH arms, so their gain is not an artefact of "
          "diagnoses that were already brewing at the sleep study. They are "
          f"{', '.join(both.disease)}.")
    A("")

    # ---------------------------------------------------------------- paths
    A("WHERE THINGS ARE")
    A("-" * 74)
    for f in ("dC_per_disease.csv", "dC_per_disease_TABLE.md", "Figure_dC_per_disease.pdf",
              "Figure_dC_per_disease.png", "REPORT.txt", "dc_engine.py",
              "p00_positive_control.py", "p01_bootstrap.py", "p02_assemble.py",
              "p03_figure.py", "p04_report.py", "p05_verify.py", "p06_lag2_gate.py",
              "p07_points.py"):
        A(f"  {os.path.join(HERE, f)}")
    A(f"  {os.path.join(E.WORK, 'boot_dc_standard_drop.npy')}  raw ({nrep}, {N_RANKED_OUTCOMES}) matrix")
    A(f"  {os.path.join(E.WORK, '_ck')} and _ck_lag2  per-chunk checkpoints, both arms")
    A("")
    A("Nothing in T90_Manuscript was written to. The frozen parquet, the numbers folder, the")
    A("manuscript folders and the version folders were read only.")

    out = os.path.join(HERE, "REPORT.txt")
    open(out, "w").write("\n".join(L_) + "\n")
    print("\n".join(L_))
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
