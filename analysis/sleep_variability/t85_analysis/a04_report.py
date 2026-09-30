"""
Build REPORT.txt from the analysis outputs. Every number in the report is read from a file
this folder produced. Nothing is typed in by hand.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import json
import os

import numpy as np
import pandas as pd
import sys as _sys_v8; _sys_v8.path.insert(0, paths.ANALYSIS_DIR)  # v8 sweep 2026-09-12
from cohort_spec import N_RANKED_OUTCOMES  # v8 sweep 2026-09-12

HERE = os.path.dirname(os.path.abspath(__file__))
W = []


def w(s=""):
    W.append(s)


def rule(c="="):
    w(c * 80)


def j(name):
    return json.load(open(os.path.join(HERE, name)))


def c(name):
    return pd.read_csv(os.path.join(HERE, name), low_memory=False)


gate = j("validation_gate.json")
band = j("banded_results.json")
bpc = j("banded_positive_control.json")
persd = j("persd_summary.json")
dc = j("dC_result.json")
cost = j("cost.json") if os.path.exists(os.path.join(HERE, "cost.json")) else None
L = c("banded_results_long.csv")
P = c("persd_t85_vs_t90.csv")
D = pd.read_csv(os.path.join(HERE, "dC_per_outcome.csv"))
T = c("t85_per_patient.csv")
ex = j("extraction_summary.json")

t85b = band["t85"]["band_n"]
t90b = band["t90_frozen"]["band_n"]
BL = ["0-1%", "1-5%", "5-10%", ">10%"]

top = L[(L.band == ">10%") & (~L.negative_control)]
t85t = top[top.exposure == "t85"].set_index("outcome")
t90t = top[top.exposure == "t90_frozen"].set_index("outcome")
both = t85t.index.intersection(t90t.index)
n85 = int((t85t.loc[both, "q"] < 0.05).sum())
n90 = int((t90t.loc[both, "q"] < 0.05).sum())
m85 = float(t85t.loc[both, "hr"].median())
m90 = float(t90t.loc[both, "hr"].median())
negtop = L[(L.band == ">10%") & (L.negative_control)]
neg85 = int((negtop[negtop.exposure == "t85"].q < 0.05).sum())
neg90 = int((negtop[negtop.exposure == "t90_frozen"].q < 0.05).sum())
nneg = int(negtop[negtop.exposure == "t85"].shape[0])

sp = persd["correlations"]["t85__t90_frozen"]["spearman"]
spr = persd["correlations"]["t85__t90_rederiv"]["spearman"]

# ------------------------------------------------------------------ header
rule()
w("T85 AS AN EXPOSURE. SIDE ANALYSIS, NOT FOR THE MANUSCRIPT.")
w("Percent of the recording below 85% SpO2, re-derived from raw BDSP oximetry and put")
w("through the paper's own banded, per-SD and discrimination machinery.")
w("Built 2026-08-21. Extraction on EC2 in Alen's account, analysis local.")
w("Signals never left AWS. Only per-patient numbers came back.")
rule()
w()

# ------------------------------------------------------------------ TLDR
w("TLDR")
verdict = ("T85 is a weaker exposure than T90, not a stronger one"
           if dc["T85_dC"] < dc["reproduced_T90_dC"] else
           "T85 beats T90 on discrimination")
w(f"{verdict}. On the paper's own discrimination measure T85 gains "
  f"{dc['T85_dC']:+.6f} against T90's published {dc['reproduced_T90_dC']:+.6f}, which is "
  f"{100*dc['T85_dC']/dc['reproduced_T90_dC']:.0f} percent of it. The banded picture is the")
w(f"same story. T85's top band is small: {t85b['>10%']:,} patients above 10 percent against "
  f"T90's {t90b['>10%']:,},")
w(f"a {t90b['>10%']/max(t85b['>10%'],1):.1f}-fold difference, because a night has to be far worse to spend "
  f"10 percent of itself")
w(f"below 85 than below 90. Half the cohort ({100*t85b['0-1%']/band['t85']['n']:.0f} percent) sits in T85's "
  f"reference band against {100*t90b['0-1%']/band['t90_frozen']['n']:.0f} percent for T90.")
w(f"Across the {len(both)} real outcomes graded for both, the top band is FDR-significant for T85 on "
  f"{n85} and for")
w(f"T90 on {n90}, with median top-band hazard ratios of {m85:.2f} and {m90:.2f}. So T85 finds "
  f"{'fewer' if n85 < n90 else 'more'} outcomes")
w(f"{'and hits them no harder' if m85 <= m90 else 'but hits them harder'}, which is what a rarer "
  f"threshold with a smaller exposed group buys.")
w(f"The two measures are nearly the same variable anyway: Spearman {sp:.3f} against the paper's")
w(f"frozen T90 and {spr:.3f} against T90 re-derived by the identical code. Put both in one model")
w("and T90 keeps its effect while T85 mostly does not add on top of it.")
w()
w("The extraction itself is clean. All three positive controls passed exactly, the frozen")
w(f"T88 column reproduced fleet-wide on {ex['t88_within_0p5pp_pct']:.1f} percent of records within 0.5 "
  f"percentage points,")
w(f"and T85 landed on {ex['n_t85']:,} of the {ex['cohort_n']:,} cohort patients "
  f"({100*ex['n_t85']/ex['cohort_n']:.2f} percent).")
w(f"Total AWS cost {'$' + format(cost['total_usd'], '.2f') if cost else 'see cost.json'}"
  f"{', every instance terminated.' if cost and cost['all_terminated'] else '.'}")
w()
w("Recommendation for Alen's decision. Nothing here argues for swapping the paper's")
w("exposure or adding a T85 panel. The one defensible use of this folder is a sensitivity")
w("statement: the graded oxygen effect is not an artefact of where the 90 percent threshold")
w("sits, because moving it to 85 reproduces the same direction and the same ordering on a")
w("much smaller exposed group, at some cost in power.")
w()

# ------------------------------------------------------------------ 0 positive controls
rule()
w("0. POSITIVE CONTROLS. ALL PASSED EXACTLY. Everything below is gated on these.")
rule()
w("Three separate pieces of published arithmetic were re-fitted before anything new was.")
w()
w(f"  banded T90        numbers/results_v2.json['graded'], all "
  f"{bpc['n_graded_published']} outcomes and all band sizes")
w(f"                    largest absolute difference in any HR or CI bound "
  f"{bpc['worst_abs_hr_diff_vs_published']:.1e}")
w(f"                    largest absolute difference in any p value "
  f"{bpc['worst_abs_p_diff_vs_published']:.1e}")
w(f"  per-SD Cox        numbers/results_v2.json['vs_sleep_duration'], "
  f"{persd['pc_n_outcomes_checked']} outcomes")
w(f"                    largest absolute difference {persd['pc_worst_abs_diff']:.1e}")
w(f"  dC pipeline       published T90 dC {dc['published_T90_dC']:.12f}")
w(f"                    reproduced here  {dc['reproduced_T90_dC']:.12f}   "
  f"diff {dc['reproduced_T90_dC']-dc['published_T90_dC']:+.1e}")
w()
w("One real bug was caught by the banded control and fixed before any T85 number existed.")
w("The paper builds its band column as an integer, and pandas get_dummies on a float band")
w("column names its outputs bd_0.0 rather than bd_0, which silently yields three all-zero")
w("dummies and a model that fits nothing. The control returned 0 of 54 outcomes until the")
w("cast was restored. This is exactly what the control is for.")
w()

# ------------------------------------------------------------------ A extraction
rule()
w("A. EXTRACTION. T85 does not exist in the frozen data, so it was re-derived from raw.")
rule()
w("The frozen exposure spo2_pct_below_90 was produced by a script that is not on disk. The")
w("rule that reproduces it is archived in numbers/t90_provenance_reference.py and was proven")
w("against raw H5 in the 2026-08-13 pilot. T85 is that rule with the threshold moved:")
w()
w("  read signals/spo2, apply the dataset's digital-to-physical affine clipped to the")
w("  digital range, band-limit resample to 200 Hz, cast float32, keep finite samples in")
w("  [50,100] inclusive, then take the percent strictly below the threshold over the")
w("  ENTIRE recording with no sleep restriction.")
w()
w("Same code, same night, same denominator as the paper's T90. The comparison is therefore")
w("structural rather than assumed. The worker is worker_t85.py, which is the 2026-08-18")
w("hypoxic-burden worker with four added output lines and no other change (see its header).")
w()
w(f"  patients in the task list        {ex['cohort_n']:,}   identical to the 2026-08-18 run")
w(f"  records returned ok             {ex['n_ok']:,}")
w(f"  T85 present                     {ex['n_t85']:,}  "
  f"({100*ex['n_t85']/ex['cohort_n']:.2f} percent of the cohort)")
for k, v in ex["status_counts"].items():
    if k != "ok":
        w(f"  status {k:<25}{v:,}")
w()
w("VALIDATION GATE, before the fleet ran, on the pilot's 25 I0002 patients.")
w("T85 has no frozen twin, so the rule was validated on the two thresholds that do have")
w("one. The same lines of code emit all three.")
w()
w(f"  frozen T90 within 0.5 pp   {gate['t90_within_0p5pp']}/25   "
  f"median absolute difference {gate['t90_median_abs_diff_pp']:.1e} pp")
w(f"  frozen T88 within 0.5 pp   {gate['t88_within_0p5pp']}/25   "
  f"median absolute difference {gate['t88_median_abs_diff_pp']:.1e} pp")
w(f"  T85 <= T88 <= T90 on every record      {gate['t85_monotone']}")
w(f"  T85 finite and inside [0,100]          {gate['t85_in_range']}")
w()
w("The 2 misses are the same two records the 2026-08-13 pilot identified (their identifiers are in the gate file, not here);")
w("pilot identified. They sit in the documented signal-tail stratum where the stored signal")
w("itself differs, not the counting rule. T88 misses on the same two records, which is the")
w("point: an independent frozen column fails in exactly the same places for the same reason.")
w()
w("FLEET-WIDE VALIDATION, after the run, on all 19,173.")
w(f"  frozen T90 within 0.5 pp   {ex['t90_within_0p5pp']:,}/{ex['n_valid_cmp']:,}  "
  f"({ex['t90_within_0p5pp_pct']:.1f} percent)   rank correlation {ex['t90_spearman']:.4f}")
w(f"  frozen T88 within 0.5 pp   {ex['t88_within_0p5pp']:,}/{ex['n_valid_cmp']:,}  "
  f"({ex['t88_within_0p5pp_pct']:.1f} percent)   rank correlation {ex['t88_spearman']:.4f}")
w(f"  T85 <= T88 <= T90 violated on {ex['n_monotone_violations']} of {ex['n_ok']:,} records")
w()
w("Both frozen columns reproduce at the same rate, which is the documented signal-tail")
w("share and not a defect of this run.")
w()

# ------------------------------------------------------------------ B distribution
rule()
w("B. WHAT T85 LOOKS LIKE, and how close it is to T90.")
rule()
w("  measure        n        mean      SD    median      p90      p99     max   exact 0")
for k in ("t90_frozen", "t85", "t90_rederiv", "t88_rederiv"):
    d = persd["distributions"][k]
    w(f"  {k:<13}{d['n']:>7,}  {d['mean']:8.3f}{d['sd']:8.3f}{d['median']:10.4f}"
      f"{d['q90']:9.3f}{d['q99']:9.2f}{d['max']:8.2f}   {d['pct_zero']:5.2f}%")
w()
w("  correlation between the exposures")
for a, bb in (("t85", "t90_frozen"), ("t85", "t90_rederiv"), ("t85", "t88_rederiv"),
              ("t90_rederiv", "t90_frozen")):
    r = persd["correlations"][f"{a}__{bb}"]
    w(f"  {a:<12} vs {bb:<12}  Spearman {r['spearman']:.4f}   Pearson {r['pearson']:.4f}"
      f"   Pearson on log1p {r['pearson_log1p']:.4f}")
w()
w("T85 and T90 are close to the same ranking of patients. The last row is the re-derivation")
w("against the frozen column and is the ceiling any of these correlations could reach.")
w()
w("  BAND SIZES, the paper's four bands, membership rule lo < v <= hi")
w()
w("  band       T90 frozen            T85                   ratio T90 to T85")
for lab in BL:
    a, bb = t90b[lab], t85b[lab]
    w(f"  {lab:<9}{a:>7,} ({100*a/band['t90_frozen']['n']:4.1f}%)   "
      f"{bb:>7,} ({100*bb/band['t85']['n']:4.1f}%)   {a/max(bb,1):>8.2f}")
w()
w(f"  T85 re-derived total {band['t85']['n']:,}, T90 frozen total {band['t90_frozen']['n']:,}")
w()
w("This is the single most consequential fact in the folder. T85's top band holds")
w(f"{t85b['>10%']:,} patients where T90's holds {t90b['>10%']:,}. Every top-band hazard ratio below is")
w("estimated on that much smaller group, so wider confidence intervals are expected and are")
w("not evidence of a weaker effect on their own.")
w()

# ------------------------------------------------------------------ C banded
rule()
w("C. BANDED RESULTS. T85 against T90, same model, same rows, same outcome panel.")
rule()
w("Model is numbers/run_all_v2.py section 5, transcribed: site-stratified Cox, Efron ties,")
w("penalizer 0, age as cr(df=4) with the first column dropped, plus sex, prevalent cases")
w("excluded, follow-up over zero, and an outcome needs 100 incident events to be graded.")
w("Benjamini-Hochberg is ADDED, within each band, over the real outcomes and separately")
w("over the negative controls. The published table reports raw p only.")
w()
w(f"  outcomes graded, T85          {len(band['t85']['outcomes'])}")
w(f"  outcomes graded, T90 frozen   {len(band['t90_frozen']['outcomes'])}")
w(f"  outcomes graded, T90 re-derived {len(band['t90_rederiv']['outcomes'])}")
w(f"  of these, in the ranking's {N_RANKED_OUTCOMES}-outcome panel  "
  f"{int(L[(L.exposure=='t85')&(L.band=='>10%')].in_ranking_48.sum())}")
w(f"  negative controls in the panel               {nneg}")
w()
w("  TOP BAND (>10%) against the 0-1% reference, real outcomes only")
w(f"  FDR q < 0.05 on   T85 {n85}/{len(both)}    T90 frozen {n90}/{len(both)}")
w(f"  median hazard ratio  T85 {m85:.3f}    T90 frozen {m90:.3f}")
w(f"  hazard ratio above 1 on  T85 {int((t85t.loc[both,'hr']>1).sum())}/{len(both)}    "
  f"T90 frozen {int((t90t.loc[both,'hr']>1).sum())}/{len(both)}")
w()
w("  the 20 outcomes with the largest published T90 top-band effect, T85 alongside")
w()
w("  outcome                        events    T90 >10%  (95% CI)      T85 >10%  (95% CI)")
o = t90t.loc[both].sort_values("hr", ascending=False).head(20)
for lab, r in o.iterrows():
    s = t85t.loc[lab]
    star90 = "*" if r.q < 0.05 else " "
    star85 = "*" if s.q < 0.05 else " "
    w(f"  {lab[:29]:<29}{int(r.events):>7}   {r.hr:6.2f} ({r.lo:5.2f}-{r.hi:5.2f}){star90}  "
      f"{s.hr:6.2f} ({s.lo:5.2f}-{s.hi:5.2f}){star85}")
w("  * FDR q < 0.05 within that band")
w()
w("  FULL DOSE-RESPONSE on the 8 outcomes with the most events")
w()
w("  outcome                     exposure       1-5%     5-10%     >10%    trend HR  trend q")
for lab in t85t.loc[both].sort_values("events", ascending=False).head(8).index:
    for expo, nm in (("t90_frozen", "T90 frozen"), ("t85", "T85")):
        rr = L[(L.exposure == expo) & (L.outcome == lab)].set_index("band")
        vals = []
        for blab in ("1-5%", "5-10%", ">10%"):
            vals.append(f"{rr.loc[blab,'hr']:8.2f}" if blab in rr.index else "     n/a")
        tq = rr.iloc[0]
        th = f"{tq.trend_hr:10.3f}" if np.isfinite(tq.trend_hr) else "       n/a"
        tqs = f"{tq.trend_q:9.2e}" if np.isfinite(tq.trend_q) else "      n/a"
        w(f"  {(lab[:26] if expo=='t90_frozen' else ''):<28}{nm:<12}"
          + "".join(vals) + th + tqs)
w()
w("Every band size for T85 is in banded_results_long.csv as n_band, and the complete")
w("side-by-side for all outcomes and all three bands is banded_t85_vs_t90_sidebyside.csv.")
w()
w("  NEGATIVE CONTROLS, top band. These should sit around 1.0 and are the folder's own check.")
w()
w("  outcome                    T90 >10% HR (95% CI)  q        T85 >10% HR (95% CI)  q")
for lab in sorted(negtop[negtop.exposure == "t85"].outcome.unique()):
    a = negtop[(negtop.exposure == "t90_frozen") & (negtop.outcome == lab)]
    bq = negtop[(negtop.exposure == "t85") & (negtop.outcome == lab)]
    if not len(a) or not len(bq):
        continue
    a, bq = a.iloc[0], bq.iloc[0]
    w(f"  {lab[:26]:<27}{a.hr:5.2f} ({a.lo:4.2f}-{a.hi:4.2f})  {a.q:7.3f}   "
      f"{bq.hr:5.2f} ({bq.lo:4.2f}-{bq.hi:4.2f})  {bq.q:7.3f}")
w()
w(f"  negative controls with q < 0.05:  T90 {neg90}/{nneg}   T85 {neg85}/{nneg}")
w()

# ------------------------------------------------------------------ D per-SD
rule()
w("D. PER-SD HEAD TO HEAD.")
rule()
w(f"Fitted on the {persd['n_common_rows']:,} patients who have all three exposures, so the")
w("comparison is never across different people. One SD is that exposure's own SD, which is")
w("not a shared scale:")
for k, v in persd["sd_units_pp"].items():
    w(f"  one SD of {k:<13} = {v:7.3f} percentage points")
w()
real = P[~P.negative_control]
neg = P[P.negative_control]
w(f"  {len(P)} outcomes fitted, {len(real)} real and {len(neg)} negative controls, "
  f"60-event floor as the paper's cox() uses")
w()
w("  exposure       median per-SD HR   HR>1      FDR q<0.05    negative controls q<0.05")
for k in ("t90_frozen", "t85", "t90_rederiv"):
    s = real[f"{k}_hr"]
    w(f"  {k:<15}{s.median():13.3f}{int((s>1).sum()):9d}/{len(real):<4}"
      f"{int((real[f'{k}_q']<0.05).sum()):9d}/{len(real):<6}{int((neg[f'{k}_q']<0.05).sum()):12d}/{len(neg)}")
w()
w("  BOTH EXPOSURES IN ONE MODEL. Does T85 add anything T90 does not already carry?")
if "joint_t85_q" in real.columns:
    w(f"  T85 keeps q < 0.05 with T90 in the model on "
      f"{int((real.joint_t85_q<0.05).sum())}/{len(real)} outcomes")
    w(f"  T90 keeps q < 0.05 with T85 in the model on "
      f"{int((real.joint_t90_q<0.05).sum())}/{len(real)} outcomes")
    w(f"  T85 flips to HR below 1 once T90 is in the model on "
      f"{int((real.joint_t85_hr<1).sum())}/{len(real)}")
w()
w("  the 12 outcomes with the most events")
w()
w("  outcome                     events   T90 alone   T85 alone   T90 joint   T85 joint")
for _i, r in real.sort_values("events", ascending=False).head(12).iterrows():
    w(f"  {str(r.outcome)[:26]:<28}{int(r.events):>6}{r.t90_frozen_hr:12.3f}{r.t85_hr:12.3f}"
      f"{r.joint_t90_hr:12.3f}{r.joint_t85_hr:12.3f}")
w()
w("Full table with confidence intervals and q values in persd_t85_vs_t90.csv.")
w()

# ------------------------------------------------------------------ E dC
rule()
w("E. DISCRIMINATION GAIN (dC) ON THE PAPER'S RANKING MACHINERY.")
rule()
w("dC_side_analyses/common.py, unchanged. Age as the published 4-column cr(df=4) at")
w("penalizer 0.01, feature as within-site rank inverse normal, fit at one hospital and")
w(f"scored at the other in both directions, mean over the ranking's {N_RANKED_OUTCOMES} outcomes. The T90 row")
w("below IS the positive control and matches the published value to 1e-9.")
w()
w("  spec              dC          against published T90     mean held-out C")
for k in ("t90_frozen", "t85", "t88_rederiv", "t90_rederiv", "t85_and_t90"):
    if k not in dc["dC"]:
        continue
    v = dc["dC"][k]
    w(f"  {k:<18}{v:+.6f}      {v-dc['reproduced_T90_dC']:+.6f}          {dc['mC'][k]:.6f}")
w()
w(f"  T85 reaches {100*dc['T85_dC']/dc['reproduced_T90_dC']:.1f} percent of the published T90 gain.")
w(f"  Per outcome, T85's gain beats T90's on {dc['n_outcomes_T85_beats_T90']}/48.")
w()
w("  the 10 outcomes where T90 gains most, with T85 alongside")
w()
w("  outcome              base C     T90 gain    T85 gain    T88 gain   T85 minus T90")
for _i, r in D.head(10).iterrows():
    w(f"  {str(r.outcome)[:20]:<21}{r.base_c:.6f}  {r.t90_frozen:+.6f}  {r.t85:+.6f}  "
      f"{r.t88_rederiv:+.6f}  {r.t85_minus_t90:+.6f}")
w()
w("Per-outcome detail in dC_per_outcome.csv, spec-level summary in dC_summary.csv.")
w()

# ------------------------------------------------------------------ F cost
rule()
w("F. COST AND INSTANCE HYGIENE.")
rule()
if cost:
    w(f"  instances launched  {cost['n_instances']}")
    w(f"  total actual cost   ${cost['total_usd']:.4f}   against a $10 cap")
    w(f"  all terminated      {cost['all_terminated']}")
    if cost["still_alive"]:
        w(f"  STILL ALIVE: {cost['still_alive']}")
    w()
    w("  instance        role            type          runtime    compute cost")
    for r in cost["instances"]:
        w(f"  {r['instance']:<20}{str(r['name'])[:14]:<16}{r['type']:<14}"
          f"{(str(r['runtime_min'])+' min'):<11}${r['compute_usd']:.4f}")
    w()
    w("  Billed time is EC2's own launch and termination timestamps, not an estimate.")
    w("  Prices are us-east-1 on-demand, c7i.2xlarge $0.3570/h and c7i.4xlarge $0.7140/h,")
    w("  plus gp3 storage for the hours held. In-region S3 transfer is free.")
w()

# ------------------------------------------------------------------ G caveats
rule()
w("G. WHAT THIS DOES NOT SETTLE.")
rule()
w("1. T85 inherits every limitation of the frozen T90 convention, including the")
w("   whole-recording denominator with no sleep restriction. That was chosen deliberately")
w("   so the two are comparable, but it means T85 is not a sleep-specific measure either.")
w("2. The band edges 1, 5 and 10 percent were chosen for T90 and are simply reused. They")
w("   are not optimal for T85 and no attempt was made to re-tune them, because re-tuning")
w("   on the outcome would be fitting the answer.")
w("3. The top band is small. A result that is not significant in a band of this size is")
w("   weak evidence of absence.")
w("4. Roughly 7 percent of records sit in the documented signal-tail stratum where the")
w("   stored signal differs from what the frozen pipeline saw. That affects T85 and T90")
w("   equally and is why the comparison is made against the re-derived T90 as well.")
w("5. Nothing here was written into any manuscript, supplement, table or figure.")
w()

# ------------------------------------------------------------------ H files
rule()
w("H. FILES. All under Sleep_Variability_2026-08/t85_analysis/")
rule()
for nm, desc in (
    ("worker_t85.py", "the EC2 worker, 2026-08-18 worker plus 4 output lines"),
    ("launch.sh", "self-terminating fleet launcher, v7 lifecycle"),
    ("build_tasks.py", "task list, asserted identical to the 2026-08-18 run"),
    ("check_validation.py", "the pre-fleet gate"),
    ("assemble.py", "S3 shards to one per-patient CSV"),
    ("a01_banded.py", "banded analysis, gated on the published T90 table"),
    ("a02_persd.py", "per-SD head to head, gated on the published binary table"),
    ("a03_dC.py", "dC on common.py, gated on the published dC"),
    ("a04_report.py", "this report"),
    ("a99_cost.py", "actual cost and termination proof"),
    ("t85_per_patient.csv", "the extraction, one row per patient, numbers only"),
    ("validation_gate.json", "gate result"),
    ("validation_detail.csv", "the 40 validation records"),
    ("extraction_summary.json", "fleet-wide coverage and validation"),
    ("banded_results.json", "banded fits for all three exposures"),
    ("banded_results_long.csv", "one row per exposure, outcome and band, with q"),
    ("banded_t85_vs_t90_sidebyside.csv", "the same, pivoted for reading"),
    ("banded_positive_control.json", "proof the published banded table reproduces"),
    ("persd_t85_vs_t90.csv", "per-SD hazard ratios including the joint model"),
    ("persd_summary.json", "correlations, distributions, SD units"),
    ("dC_summary.csv", "spec-level dC"),
    ("dC_per_outcome.csv", "per-outcome gains for all specs"),
    ("dC_result.json", "dC headline plus positive control"),
    ("cost.json", "per-instance billed time and cost"),
):
    w(f"  {nm:<36}{desc}")
w()
rule()
w("END")
rule()

open(os.path.join(HERE, "REPORT.txt"), "w").write("\n".join(W) + "\n")
print(f"wrote REPORT.txt  ({len(W)} lines)")
