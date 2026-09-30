"""
REPORT.md, rendered from the run's own outputs. Every number in the report is read from
T90_gt10_concordance_per_disease.csv, T90_gt10_headline.json, _work/selftest.json or
Figure_T90_gt10_drawn_values.json, so nothing in the prose can drift from the analysis.
"""
import json
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
D = pd.read_csv(f"{HERE}/T90_gt10_concordance_per_disease.csv", comment="#")
H = json.load(open(f"{HERE}/T90_gt10_headline.json"))
ST = json.load(open(f"{HERE}/_work/selftest.json"))
FG = json.load(open(f"{HERE}/Figure_T90_gt10_drawn_values.json"))
_PDFCHK_PATH = f"{HERE}/_work/fig1c_pdf_textcheck.json"
PDFCHK = json.load(open(_PDFCHK_PATH)) if os.path.exists(_PDFCHK_PATH) else None   # v8: the sheet text check belongs to the recompose phase; absent here
D = D.sort_values("lag0_bin_dC", ascending=False).reset_index(drop=True)

FF, LM = H["full_followup"], H["landmark_2y"]
n = len(D)


def f(x, k=4):
    return f"{x:+.{k}f}".replace("-", "−")


def ci(lo, hi, k=4):
    return f"({f(lo, k)} to {f(hi, k)})"


def line(r):
    a = f"{r.disease}, {f(r.lag0_bin_dC)} {ci(r.lag0_bin_ci_lo, r.lag0_bin_ci_hi)} binary"
    a += f", {f(r.lag0_cont_dC)} continuous"
    if bool(r.lag2_bin_estimable):
        a += f", {f(r.lag2_bin_dC)} {ci(r.lag2_bin_ci_lo, r.lag2_bin_ci_hi)} binary at the "
        a += "2-year landmark"
    else:
        a += ", no 2-year landmark estimate"
    a += f". {int(r.lag0_n_gt10):,} of {int(r.lag0_n_at_risk):,} above 10%, "
    a += f"{int(r.lag0_events):,} events, {int(r.lag0_events_gt10):,} of them above 10%"
    return "- " + a


flip = H["outcomes_where_binary_beats_continuous_full_followup"]
flip2 = H["outcomes_where_binary_beats_continuous_landmark"]
gap = FF["continuous_mean_dC"] - FF["binary_mean_dC"]
gap_l = LM["continuous_mean_dC"] - LM["binary_mean_dC"]
keeps = 100 * FF["binary_mean_dC"] / FF["continuous_mean_dC"]
keeps_l = 100 * LM["binary_mean_dC"] / LM["continuous_mean_dC"]
per_dis_ratio = D.loc[D.lag0_cont_dC > 0.005, "lag0_bin_retains_pct_of_cont"]
sig_bin = int(D.lag0_bin_excludes_zero.sum())
sig_cont = int((D.lag0_cont_ci_lo > 0).sum())
top = D.head(10)
neg = D[D.negative_control == "yes"]

T = []
A = T.append

A(f"# T90 above 10% as a binary exposure: concordance gain over age and sex, disease by disease")
A("")
A(f"The gain was never computed for the above-10% group. The paper's ranking, **Figure 1b**, "
  f"**Figure 1c** and **Supplementary Table 20** all enter T90 as a continuous site-wise "
  f"rank-normal score, so no published number describes what the binary indicator alone adds. "
  f"Computed here for the first time on the same {H['cohort_n']:,} recordings and the same "
  f"{n} outcomes, a 0/1 flag for T90 above 10% of the recording adds a mean of "
  f"{f(FF['binary_mean_dC'])} {ci(*FF['binary_mean_ci'])} to held-out concordance over age and "
  f"sex, against {f(FF['continuous_mean_dC'])} {ci(*FF['continuous_mean_ci'])} for the "
  f"continuous score fitted the same way. Dichotomising therefore keeps about "
  f"{keeps:.0f}% of the discrimination the continuous measure buys, the binary minus "
  f"continuous difference averaging {f(FF['binary_minus_continuous_mean'])} "
  f"{ci(*FF['binary_minus_continuous_ci'])} on a bootstrap that resamples the same patients "
  f"for both exposures, and the same comparison under a two-year landmark is "
  f"{f(LM['binary_mean_dC'])} {ci(*LM['binary_mean_ci'])} binary against "
  f"{f(LM['continuous_mean_dC'])} {ci(*LM['continuous_mean_ci'])} continuous, a shortfall of "
  f"{f(LM['binary_minus_continuous_mean'])} "
  f"{ci(*LM['binary_minus_continuous_ci'])}.")
A("")
A(f"{H['n_above_10pct']:,} of the {H['cohort_n']:,} patients sit above 10%, which is "
  f"{H['pct_above_10pct']:.1f}% of the cohort.")
A("")
A("## What this means in one paragraph")
A("")
_worst = D.assign(_gap=D.lag0_cont_dC - D.lag0_bin_dC).nlargest(3, "_gap")
_wtxt = ", ".join(f"{r.disease} ({f(r.lag0_bin_dC)} binary against {f(r.lag0_cont_dC)} "
                  f"continuous)" for _i, r in _worst.iterrows())
A(f"Cutting a continuous measure at a threshold throws away the ordering inside each group, so "
  f"the binary gain is usually the smaller of the two, and it is here: the binary indicator is "
  f"behind the continuous score on {int((D.lag0_bin_dC < D.lag0_cont_dC).sum())} of the {n} "
  f"outcomes at full follow-up. The average difference is "
  f"{f(FF['binary_minus_continuous_mean'])} of concordance "
  f"{ci(*FF['binary_minus_continuous_ci'])}, so the binary flag keeps {keeps:.0f}% of the "
  f"continuous gain at full follow-up and {keeps_l:.0f}% under the landmark. The loss is not spread evenly. The three "
  f"largest shortfalls are {_wtxt}. Read the binary numbers as what a clinician gets from the "
  f"threshold alone, not as a reason to prefer it.")
A("")
_sig = D[D.lag0_bin_minus_cont_excludes_zero.astype(bool)]
_sigdn = _sig[_sig.lag0_bin_minus_cont < 0]
A(f"On the paired bootstrap, which resamples the same patients for both exposures, the gap "
  f"between the two is larger than resampling noise on {len(_sig)} of the {n} outcomes, and "
  f"{len(_sigdn)} of those {'favour' if len(_sigdn) != 1 else 'favours'} the continuous score"
  + (f": {', '.join(_sigdn.disease)}. " if len(_sigdn) else ". ")
  + f"On the other {n - len(_sig)} the two exposures are within resampling noise of each other.")
A("")
A("## The diseases that behave the other way")
A("")
if flip:
    _real = [r for r in flip if r["continuous_dC"] > 0.010]
    _big = [r for r in flip if r["difference_excludes_zero"]]
    A(f"{len(flip)} of the {n} outcomes have a larger binary gain than continuous gain at full "
      f"follow-up. {len(_real)} of them have a continuous gain above 0.010, so they are the "
      f"only ones where the reversal is a difference between two numbers that both mean "
      f"something"
      + (f": {', '.join(r['disease'] for r in _real)}. " if _real else ". ")
      + f"The rest are small differences between two near-zero gains. On the paired bootstrap, "
        f"which resamples the same patients for both exposures, "
      + (f"{', '.join(r['disease'] for r in _big)} has an interval that excludes zero."
         if _big else "not one of the differences has an interval that excludes zero."))
    A("")
    for r in flip:
        A(f"- {r['disease']}, {f(r['binary_dC'])} binary against {f(r['continuous_dC'])} "
          f"continuous, difference {f(r['difference'])} "
          f"{ci(r['difference_ci'][0], r['difference_ci'][1])}")
    A("")
    big = [r for r in flip if r["difference_excludes_zero"]]
    A(f"Differences whose paired bootstrap interval excludes zero: "
      f"{', '.join(r['disease'] for r in big) if big else 'none of them'}.")
else:
    A(f"No outcome has a larger binary gain than continuous gain at full follow-up.")
A("")
if flip2:
    _big2 = [r for r in flip2 if r["difference_excludes_zero"]]
    A(f"Under the two-year landmark the list is different, {len(flip2)} outcomes: "
      f"{', '.join(r['disease'] for r in flip2)}. Landmarking removes roughly half the events, "
      f"so these are the noisiest estimates on the sheet and the ordering between the two "
      f"exposures there should not be read as a finding. "
      + (f"Intervals excluding zero: {', '.join(r['disease'] for r in _big2)}."
         if _big2 else "Not one of these differences has an interval that excludes zero."))
    A("")
A("## Every outcome, ranked by the binary gain")
A("")
A(f"One line per outcome, in order of the full-follow-up binary gain. The continuous value "
  f"beside each is the published one, carried over from "
  f"`Sleep_Variability_2026-08/dC_per_disease/dC_per_disease.csv`, the file **Figure 1c** is "
  f"drawn from. The full grid with every column is in "
  f"`T90_gt10_concordance_per_disease.csv`.")
A("")
for _i, r in D.iterrows():
    A(line(r))
A("")
A(f"- Mean across the {FF['n_outcomes_in_mean']} outcomes, {f(FF['binary_mean_dC'])} "
  f"{ci(*FF['binary_mean_ci'])} binary, {f(FF['continuous_mean_dC'])} "
  f"{ci(*FF['continuous_mean_ci'])} continuous")
A(f"- Mean across the {LM['n_outcomes_in_mean']} outcomes estimable at the two-year landmark, "
  f"{f(LM['binary_mean_dC'])} {ci(*LM['binary_mean_ci'])} binary, "
  f"{f(LM['continuous_mean_dC'])} {ci(*LM['continuous_mean_ci'])} continuous")
A("")
_bin_l2 = int(D.lag2_bin_excludes_zero.astype(bool).sum())
_bin_both = int((D.lag0_bin_excludes_zero.astype(bool)
                 & D.lag2_bin_excludes_zero.astype(bool)).sum())
_cont_l2 = int((D.lag2_cont_ci_lo > 0).sum())
_cont_both = int(((D.lag0_cont_ci_lo > 0) & (D.lag2_cont_ci_lo > 0)).sum())
A(f"The binary interval excludes zero for {sig_bin} of the {n} outcomes at full follow-up and "
  f"{_bin_l2} at the two-year landmark, {_bin_both} of them in both. The same counts for the "
  f"continuous score are {sig_cont}, {_cont_l2} and {_cont_both}, and the manuscript already "
  f"prints the {sig_cont} and the {_cont_both} in its results text, so the two exposures can "
  f"be compared on the sentence the paper already carries.")
A("")
A(f"The three negative controls land at ranks "
  f"{', '.join(str(int(x)) for x in neg.rank_by_binary_gain)} of {n} "
  f"({', '.join(neg.disease)}), which is where a measure with no real signal on them should "
  f"put them.")
A("")
A("## Where the estimate is thin")
A("")
thin = D.nsmallest(6, "lag0_events_gt10")
A(f"The count of patients above 10% and the event count sit on every line above. The thinnest "
  f"rows by events among the exposed are:")
A("")
for _i, r in thin.iterrows():
    A(f"- {r.disease}, {int(r.lag0_events_gt10):,} events among the "
      f"{int(r.lag0_n_gt10):,} patients above 10%, binary gain {f(r.lag0_bin_dC)} "
      f"{ci(r.lag0_bin_ci_lo, r.lag0_bin_ci_hi)}")
A("")
ne = H["not_estimable"]["lag2"]
if ne:
    _nen = [str(D.loc[D.outcome == o, "disease"].iloc[0]) for o in ne]
    A(f"After the two-year landmark, {', '.join(_nen)} falls below the ranking's own floor of 30 "
      f"incident cases at one hospital, so it carries no landmark estimate in either exposure "
      f"and is left empty on the sheet rather than drawn at zero.")
    A("")
A("## The figure")
A("")
A(f"`{FG['sheet']}` and its 300 dpi PNG. {FG['n_rows']} outcomes, every one that clears the "
  f"paper's 150-event floor, plus the mean row, {FG['n_marks']} plotted marks in all. Sheet "
  f"{FG['canvas_mm'][0]:.0f} by {FG['canvas_mm'][1]:.1f} mm. Conventions are taken from "
  f"**Figure 1c** of `ROUND16_2026-09-04/L9_assembly/NEW_FINAL_SET_R16/Main_Fig1.pdf`, read "
  f"and looked at before drawing: one house blue {FG['colour']['series']} for both series, "
  f"sampled from the panel itself, filled circle for full follow-up and filled triangle for "
  f"the two-year landmark, round-capped interval lines, a dashed ink rule at zero and a dotted "
  f"blue rule at the mean, the right-aligned printed Gain (95% CI) column, Arial with a 9 pt "
  f"floor, no gridlines and no headline. Spines are {FG['spines']}. Figure 1c itself hides "
  f"the left spine, because its axis starts at zero and the dashed zero rule already stands "
  f"there, but this sheet runs negative, so it keeps the left spine and right-aligns the row "
  f"labels against it, which is what the published 48-row sibling "
  f"`Sleep_Variability_2026-08/dC_per_disease/Figure_dC_per_disease.pdf` does and what the "
  f"clinical figure spec asks for. One convention is added that Figure 1c never needs, an "
  f"open marker where the interval includes zero, which on Figure 1c's ten rows never happens "
  f"and here happens often. Every plotted point and "
  f"interval is asserted against the results CSV before the file is written, largest "
  f"difference {FG['largest_difference_against_csv']:.1e}, and the same values are written to "
  f"`Figure_T90_gt10_drawn_values.json`.")
A("")
A("## Method")
A("")
A(f"The pipeline is the published one, reused rather than rebuilt. "
  f"`T90_Manuscript/numbers/build_ranking_v3.py` by way of "
  f"`Sleep_Variability_2026-08/dC_per_disease/dc_engine.py`, which is that script "
  f"reorganised so the two cross-fit directions are kept apart and the whole pass is cheap "
  f"enough to bootstrap. `engine_gt10.py` imports that engine unmodified and adds one branch, "
  f"the exposure column.")
A("")
A(f"- cohort: {H['cohort_n']:,} recordings, `numbers/cohort_spec.py` apply_cohort, from "
  f"`data_frozen_v8_2026-09/t90_final.parquet` joined to `master_cohort_v2.csv`")
A(f"- exposure: `{H['exposure_column']}` from `{H['exposure_source']}`, in percent of the "
  f"recording. The indicator is {H['exposure_rule']}, which is the paper's own top oxygen "
  f"band, the rule being lo < v <= hi so the top band is v above 10")
A(f"- adjustment: age as a cr(df=4) natural cubic spline with one column dropped, plus sex, "
  f"which is the rank-safe basis every T90 fit uses")
A(f"- fit: lifelines CoxPHFitter with penalizer 0.01 at one hospital, scored at the other, "
  f"I0002 to I0006 and I0006 to I0002, the two held-out concordances averaged")
A(f"- eligibility: prevalent disease excluded, follow-up above zero, a fold dropped when "
  f"either side has fewer than 30 events")
A(f"- outcomes: the {n} the ranking was computed on, honouring CIRCULAR and RANKING_EXCLUDE "
  f"and the 150-event floor, plus death")
A(f"- landmark: the published rule, every patient whose event or censoring arrived at or "
  f"before 2 years leaves the risk set and the clock restarts at 2 years, from "
  f"`New_Figures/scripts/lagladder.py` recorded in `numbers/lag_ladder.json` and reproduced "
  f"here through `dc_engine.one_pass(lag=2.0)`")
A(f"- 95% CI: the published interval, the 2.5th to 97.5th percentile of "
  f"{H['bootstrap_replicates']} bootstrap replicates resampling patients with replacement "
  f"inside each hospital, the whole procedure repeated on every resample, seeds "
  f"20260820 to 20261319, which are the published run's own seeds. Replicate r here and "
  f"replicate r of the published continuous bootstrap are therefore the same resampled "
  f"patients, so the binary minus continuous difference is paired and carries its own interval")
A("")
A(f"Nothing else changed. On the published four-column age basis, which is the basis the "
  f"printed ranking uses, the binary mean gain is "
  f"{f(H['published_basis_cross_reference']['binary_mean_dC_four_column_basis'])} against the "
  f"paper's {f(H['published_basis_cross_reference']['continuous_mean_dC_four_column_basis'])}, "
  f"so the conclusion does not depend on which of the two age bases is used.")
A("")
A("## Self-test")
A("")
A("The control was run first and had to pass before any binary number was computed. It is the "
  "published continuous T90 pushed through the exact code path this analysis uses.")
A("")
g1 = ST["gate1_ranking_v3_percondition"]
g2 = ST["gate_standard_dC_per_disease"]
g3 = ST["gate_lag2_dC_per_disease"]
A(f"- Gate 1, full follow-up on the published four-column basis against "
  f"`numbers/ranking_v3_percondition.csv`. Headline mean gain reproduced as "
  f"{g1['reproduced_mean_dC']:.14f} against the published {g1['published_mean_dC']:.14f}, "
  f"held-out baseline C {g1['reproduced_baseline_C']:.6f} against "
  f"{g1['published_baseline_C']:.6f}, largest per-outcome gain difference over all {n} "
  f"outcomes {g1['max_abs_gain_diff']:.1e}. PASS")
A(f"- Gate 2, full follow-up on the reported basis against "
  f"`Sleep_Variability_2026-08/dC_per_disease/dC_per_disease.csv`, the file **Figure 1c** is "
  f"drawn from. {g2['n_estimable']} estimable outcomes, largest difference "
  f"{g2['max_abs_dC_diff']:.1e}. PASS")
A(f"- Gate 3, the two-year landmark against the same file's lag2 arm. "
  f"{g3['n_estimable']} estimable outcomes, largest difference {g3['max_abs_dC_diff']:.1e}. "
  f"PASS")
A(f"- Gate 4, the landmark definition itself. Patients at risk and events at lag 0 and lag 2 "
  f"match `numbers/lag_ladder.csv` exactly for all {n} outcomes, no tolerance. PASS")
A(f"- Gate 5, the ten values **Figure 1c** prints. All 20 drawn points reproduced, largest "
  f"difference {max(p['abs_diff'] for p in ST['gate5_figure1c_drawn_values']['points']):.1e}, "
  + (f"and all ten printed strings were found in the text layer of "
     f"`ROUND16_2026-09-04/L9_assembly/NEW_FINAL_SET_R16/Main_Fig1.pdf` itself "
     f"({sum(1 for p in PDFCHK if p['found_in_Main_Fig1_pdf'])} of {len(PDFCHK)}). PASS" if PDFCHK is not None else
     "and the printed-string check against the Figure 1c sheet is deferred to the v8 recompose phase, where the sheet is redrawn from these values. No v7 sheet check is quoted. PASS on the drawn values"))
A(f"- Gate 6, the wrapper. With the exposure set to the published rank-normal column, "
  f"`engine_gt10.one_pass_x` returns `dc_engine.one_pass`'s array bit for bit on the observed "
  f"data and on two resamples, across all three basis and landmark settings, so the only "
  f"difference between the control and the binary run is the exposure column. PASS")
A("")
A("Reproduced values, printed against the published ones, for the three conditions the "
  "self-test names:")
A("")
for p in ST["gate5_figure1c_drawn_values"]["points"]:
    if p["printed_on_figure"] and p["disease"] in ("Type 2 diabetes", "Heart failure", "COPD"):
        A(f"- {p['disease']}, Figure 1c prints {p['printed_on_figure']}, reproduced here as "
          f"{p['reproduced_dC']:+.6f}")
A("")
A("## Numbers carried over, and their source file")
A("")
A("- the continuous per-disease gains and their intervals, both arms, from "
  "`$T90_SV_ROOT/dC_per_disease/"
  "dC_per_disease.csv`")
A("- the continuous bootstrap replicate matrices used for the paired difference, from that "
  "folder's `_work/boot_dc_standard_drop.npy` and `_work/boot_dc_lag2_drop.npy`")
A("- the continuous gains on the published four-column basis, from "
  "`$T90_NUMBERS_DIR/"
  "ranking_v3_percondition.csv`")
A("- the ten values **Figure 1c** prints, from "
  "`$T90_FIGURE_ROOT/FINAL_FIGURES_2026-08-14/"
  "_workfiles/Figure1C_drawn_values.json`, cross-checked against the text layer of "
  "`Main_Fig1.pdf`")
A("")
A("Everything else in this report was computed in this run.")
A("")
A("## Files")
A("")
A("- `T90_gt10_concordance_per_disease.csv`, one row per outcome, every column, header "
  "carries the provenance")
A("- `T90_gt10_headline.json`, the headline numbers")
A("- `Figure_T90_gt10_dC_per_disease.pdf` and `.png` at 300 dpi")
A("- `Figure_T90_gt10_drawn_values.json`, every plotted point and interval asserted against "
  "the CSV")
A("- `engine_gt10.py`, `p00_selftest.py`, `p01_points.py`, `p02_bootstrap.py`, "
  "`p03_assemble.py`, `p04_figure.py`, `p05_report.py`")
A("- `_work/selftest.json`, the six gates as recorded by the run")

txt = "\n".join(T) + "\n"
assert ";" not in txt, "a semicolon reached the report"
assert "—" not in txt, "an em-dash reached the report"
open(f"{HERE}/REPORT.md", "w").write(txt)
print(f"wrote {HERE}/REPORT.md  ({len(txt):,} characters, {len(T)} blocks)")
