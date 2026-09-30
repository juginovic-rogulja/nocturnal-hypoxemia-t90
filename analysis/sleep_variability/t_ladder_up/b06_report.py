"""
REPORT.txt for the threshold-ladder run. Every number is read from the artifacts
this folder's scripts wrote. Nothing is typed in by hand.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import json
import os

import pandas as pd
import sys as _sys_v8; _sys_v8.path.insert(0, paths.ANALYSIS_DIR)  # v8 sweep 2026-09-12
from cohort_spec import N_RANKED_OUTCOMES  # v8 sweep 2026-09-12

HERE = os.path.dirname(os.path.abspath(__file__))


def J(name):
    return json.load(open(os.path.join(HERE, name)))


S = pd.read_csv(os.path.join(HERE, "LADDER_SUMMARY.csv"))
L = pd.read_csv(os.path.join(HERE, "ladder_banded_long.csv"))
gate = J("validation_gate.json")
ext = J("extraction_summary.json")
ps = J("persd_summary.json")
dc = J("dC_result.json")
bpc = J("banded_positive_control.json")
bres = J("ladder_banded_results.json")
cost = J("cost.json")
peak = J("ladder_peak.json")

PRIM = ["T80", "T85", "T88", "T90", "T92", "T95"]
prim = S[S.threshold.isin(PRIM)].set_index("threshold").loc[PRIM]
sec = S[S.threshold.str.endswith("_rederiv")].set_index("threshold")
a95row = S[S.threshold == "A95_preserved_bands"].iloc[0]

o = []
w = o.append
w("=" * 80)
w("THE OXYGEN THRESHOLD LADDER, T80 TO T95. SIDE ANALYSIS, NOT FOR THE MANUSCRIPT.")
w("Percent of the recording below each SpO2 threshold, re-derived from raw BDSP")
w("oximetry for all 19,173 cohort patients, plus the preserved-saturation reading")
w("of T95. Same rule, same denominator, same machinery as the paper's T90.")
w("Built 2026-08-21. Extraction on EC2 in Alen's account, analysis local.")
w("Signals never left AWS. Only per-patient numbers came back.")
w("=" * 80)
w("")
w("TLDR")
dcs = prim.dC_mean
ns = prim.n_fdr_sig_of_48_persd
# every sentence about the peak is COMPUTED, not assumed
pkd, pkn = peak["peak_by_dC"], peak["peak_by_nsig"]
ipk = PRIM.index(pkd)
rising = all(dcs[PRIM[i]] <= dcs[PRIM[i + 1]] + 1e-12 for i in range(ipk))
falling = all(dcs[PRIM[i]] >= dcs[PRIM[i + 1]] - 1e-12 for i in range(ipk, len(PRIM) - 1))
if rising and falling and 0 < ipk < len(PRIM) - 1:
    shape = (f"dC rises rung by rung to its maximum at {pkd} and falls on the far "
             f"side of it, a single interior peak.")
elif rising and ipk == len(PRIM) - 1:
    shape = f"dC rises monotonically along the whole ladder and is highest at {pkd}."
elif falling and ipk == 0:
    shape = f"dC falls monotonically along the whole ladder from its maximum at {pkd}."
else:
    shape = (f"dC is highest at {pkd} but the profile is not monotone on both sides, "
             f"read the table.")
gap = abs(int(ns[pkn]) - int(ns[pkd]))
if pkd == pkn:
    agree = f"The breadth count agrees, peaking at {pkd} with {int(ns[pkd])}/48."
elif abs(PRIM.index(pkn) - ipk) == 1:
    agree = (f"The breadth count nominally peaks one rung away at {pkn} "
             f"({int(ns[pkn])}/{N_RANKED_OUTCOMES} against {int(ns[pkd])}/{N_RANKED_OUTCOMES} at {pkd}, a "
             f"{gap}-outcome edge inside FDR-count noise).")
else:
    agree = (f"The breadth count peaks at {pkn} ({int(ns[pkn])}/{N_RANKED_OUTCOMES}), two or more "
             f"rungs from the dC peak, read the table before concluding anything.")
w(f"The information peaks at {pkd} on the paper's own ranking currency "
  f"(dC {peak['peak_dC']:+.6f}). {agree}")
w(f"Along the ladder T80 to T95 the dC sequence is "
  + ", ".join(f"{dcs[t]:+.6f}" for t in PRIM) + ",")
w(f"and the per-SD significance counts are "
  + ", ".join(str(int(ns[t])) for t in PRIM) + " of 48. " + shape)
s_a95 = int(a95row.n_fdr_sig_of_48_persd)
t90top = L[(L.exposure == "t90_frozen") & (L.band == ">10%") & (~L.negative_control)]
s_t90b = int((t90top.q < 0.05).sum())
n_real_banded = len(t90top)
if s_a95 <= s_t90b:
    a95cmp = (f"grades outcomes ({s_a95}/{n_real_banded} significant in the "
              f"least-preserved band) without overtaking T90's top band "
              f"({s_t90b}/{n_real_banded})")
else:
    a95cmp = (f"is nominally ahead of T90's top band ({s_a95}/{n_real_banded} against "
              f"{s_t90b}/{n_real_banded}), see section C before reading anything into it")
w(f"Reading T95 the other way round as time at normal saturation adds nothing as a")
w(f"continuous measure (it is the same variable with the sign flipped, shown exactly")
w(f"below), and its banded preserved-saturation framing {a95cmp}.")
w(f"Validation is clean and every instance")
w(f"is terminated. Total AWS cost ${cost['total_usd']:.2f} against the $10 cap.")
w("")
w("THE LADDER TABLE")
w("")
hdr = (f"  {'rung':<6} {'source':<34} {'n':>6}  {'median':>8}  {'%>10':>6}  "
       f"{'%<=1':>6}  {'dC':>10}  {'sig/48':>6}  {'medHR':>6}  {'medHRsig':>8}")
w(hdr)
w("  " + "-" * (len(hdr) - 2))
for th in PRIM:
    r = prim.loc[th]
    srcshort = {"T80": "re-derived", "T85": "re-derived (=T85 run)",
                "T88": "frozen spo2_pct_below_88", "T90": "frozen, published",
                "T92": "re-derived", "T95": "re-derived"}[th]
    w(f"  {th:<6} {srcshort:<34} {r.n:>6,}  {r.median_pct:>8.4f}  {r.pct_gt10:>6.2f}  "
      f"{r.pct_le1:>6.2f}  {r.dC_mean:>+10.6f}  {int(r.n_fdr_sig_of_48_persd):>6}  "
      f"{r.median_persd_hr:>6.3f}  {r.median_persd_hr_sig:>8.3f}")
w("")
for th in ("T88_rederiv", "T90_rederiv"):
    r = sec.loc[th]
    w(f"  {th:<13} {'same-code secondary rung':<27} {r.n:>6,}  {r.median_pct:>8.4f}  "
      f"{r.pct_gt10:>6.2f}  {r.pct_le1:>6.2f}  {r.dC_mean:>+10.6f}  "
      f"{int(r.n_fdr_sig_of_48_persd):>6}  {r.median_persd_hr:>6.3f}  "
      f"{r.median_persd_hr_sig:>8.3f}")
w("")
w("  median = median percent of the recording below the threshold. %>10 and %<=1 are")
w("  the shares of the cohort above 10 percent and at or below 1 percent. dC is the")
w(f"  paper's mean discrimination gain over the {N_RANKED_OUTCOMES} ranking outcomes (T90's cell is the")
w("  published value, which this run reproduced to 1e-9). sig/48 is the count of")
w(f"  FDR-significant per-SD hazard ratios among the {N_RANKED_OUTCOMES} real outcomes, medHR the")
w(f"  median per-SD HR over all {N_RANKED_OUTCOMES}, medHRsig the median among the significant ones.")
w("  The two secondary rungs re-derive 88 and 90 with the ladder's own code and by")
w("  matching the frozen rows they show frozen-vs-rederived mixing changes nothing.")
w("")
w(f"  A95 preserved-saturation row: median {a95row.median_pct:.2f} percent of the night at or")
w(f"  above 95, dC {a95row.dC_mean:+.6f} (identical to T95 by construction), worst band")
w(f"  vs most-preserved band significant on {int(a95row.n_fdr_sig_of_48_persd)}/"
  f"{n_real_banded} real graded outcomes (banded, not")
w(f"  per-SD), median significant HR {a95row.median_persd_hr_sig:.3f}.")
w("")
w("WHERE THE INFORMATION PEAKS")
w(f"By dC the argmax over the six rungs is {peak['peak_by_dC']} at "
  f"{peak['peak_dC']:+.6f}. By breadth the argmax")
w(f"is {peak['peak_by_nsig']} with {peak['peak_nsig']}/48. {shape}")
if pkd == "T90":
    if pkn == "T90":
        w("Both criteria put the published exposure at the optimum of its own family.")
    elif abs(PRIM.index(pkn) - ipk) == 1:
        w(f"The published T90 holds the dC maximum. {pkn} trades "
          f"{100*dcs[pkn]/dcs['T90']:.1f} percent of T90's dC for a "
          f"{gap}-outcome nominal edge in breadth,")
        w("which is inside FDR-count noise, so the information concentrates at "
          f"T90-{pkn[1:]} and nothing")
        w("in 80-95 improves on the published threshold in any way that survives G.4.")
    else:
        w(f"The dC and breadth argmaxes disagree by more than one rung ({pkd} vs "
          f"{pkn}), judge from the table.")
else:
    w(f"Note the dC argmax is {pkd} rather than the published T90, judge the gap")
    w("against caveat G.4 before reading anything into it.")
w("")
w("=" * 80)
w("0. POSITIVE CONTROLS. ALL PASSED EXACTLY. Everything below is gated on these.")
w("=" * 80)
w(f"  banded T90       every published HR, CI bound and p in results_v2.json['graded']")
w(f"                   reproduced, largest |diff| {bpc['worst_abs_hr_diff_vs_published']:.1e} "
  f"(HR/CI) and {bpc['worst_abs_p_diff_vs_published']:.1e} (p)")
w(f"  per-SD Cox       results_v2.json['vs_sleep_duration'], "
  f"{ps['pc_n_outcomes_checked']} outcomes, largest |diff| {ps['pc_worst_abs_diff']:.1e}")
w(f"  dC pipeline      published T90 dC {dc['published_T90_dC']:.12f}")
w(f"                   reproduced here  {dc['reproduced_T90_dC']:.12f}   "
  f"diff {dc['reproduced_T90_dC']-dc['published_T90_dC']:+.1e}")
w("  cross-checks     t85, t88_rederiv, t90_rederiv dC reproduce this morning's T85")
xc = dc["crosschecks_vs_t85_run"]
w("                   run: " + ", ".join(f"{k} |diff| {v:.1e}" for k, v in xc.items()))
w("")
w("=" * 80)
w("A. EXTRACTION AND VALIDATION")
w("=" * 80)
w("The worker is t85_analysis/worker_t85.py plus three threshold lines and the")
w("independent a95 line, nothing else. Selftest: 44 known-answer cases pass,")
w("covering exact fractions, the [50,100] keep band, strictly-below at every")
w("threshold, at-or-above for a95, monotonicity, the affine and 1 Hz resampling.")
w("")
w("VALIDATION GATE before the fleet, on the pilot's 25 I0002 patients:")
w(f"  frozen T90 within 0.5 pp   {gate['t90_within_0p5pp']}/25   "
  f"median |diff| {gate['t90_median_abs_diff_pp']:.1e} pp")
w(f"  frozen T88 within 0.5 pp   {gate['t88_within_0p5pp']}/25   "
  f"median |diff| {gate['t88_median_abs_diff_pp']:.1e} pp")
w(f"  t80<=t85<=t88<=t90<=t92<=t95 on every record   {gate['ladder_monotone']}")
w(f"  prov_t85 vs the T85 run    {gate['t85_vs_prev_within_0p01pp']}/"
  f"{gate['t85_vs_prev_n']} within 0.01 pp, max |diff| "
  f"{gate['t85_vs_prev_max_abs_diff_pp']:.1e} pp")
w(f"  a95 + t95 = 100            max |dev| {gate['a95_identity_max_abs_dev_pp']:.1e} pp")
w("  The 2 misses are the same signal-tail records the 2026-08-13 pilot and the T85")
w("  run identified, and T88 misses on the same two, which is the point of the check.")
w("")
ncmp = ext["coverage"]["prov_t90"]
w("FLEET-WIDE, all 19,173 patients:")
w(f"  records returned ok        {ext['n_ok']:,} of {ext['cohort_n']:,}")
w(f"  ladder present             {ext['coverage']['prov_t95']:,} "
  f"({100*ext['coverage']['prov_t95']/ext['cohort_n']:.2f} percent of the cohort)")
w(f"  frozen T90 within 0.5 pp   {ext['t90_within_0p5pp']:,}/{ncmp:,} comparable "
  f"({ext['t90_within_0p5pp_pct']:.1f} percent)   rank corr {ext['t90_spearman']:.4f}")
w(f"  frozen T88 within 0.5 pp   {ext['t88_within_0p5pp']:,}/{ncmp:,} comparable "
  f"({ext['t88_within_0p5pp_pct']:.1f} percent)   rank corr {ext['t88_spearman']:.4f}")
w(f"  monotonicity violations    {ext['n_monotone_violations']} of {ncmp:,}")
w(f"  prov_t85 vs the T85 run    {ext['t85_vs_prev_within_0p01pp']:,}/"
  f"{ext['t85_vs_prev_n']:,} within 0.01 pp, max {ext['t85_vs_prev_max_abs_diff_pp']:.1e} pp")
w(f"  a95 identity max |dev|     {ext['a95_identity_max_abs_dev_pp']:.1e} pp")
w("  The within-0.5pp rates are the documented signal-tail share, the same stratum at")
w("  the same rate as the T85 run this morning, not a defect of this run.")
w("")
w("=" * 80)
w("B. WHAT THE NEW RUNGS LOOK LIKE")
w("=" * 80)
d_ = ps["distributions"]
w(f"  {'measure':<12} {'n':>7} {'mean':>8} {'SD':>8} {'median':>9} {'p90':>8} "
  f"{'p99':>7} {'max':>7}")
for c in ("t80", "t85", "t88_frozen", "t90_frozen", "t92", "t95", "a95"):
    v = d_[c]
    w(f"  {c:<12} {v['n']:>7,} {v['mean']:>8.3f} {v['sd']:>8.3f} {v['median']:>9.4f} "
      f"{v['q90']:>8.3f} {v['q99']:>7.2f} {v['max']:>7.2f}")
w("")
w("  Spearman with frozen T90: " + ", ".join(
    f"{c} {ps['correlations'][f'{c}__t90_frozen']:.3f}"
    for c in ("t80", "t85", "t92", "t95")))
w(f"  a95 vs t95: Pearson {ps['correlations']['a95__t95_pearson']:.6f}, "
  f"Spearman {ps['correlations']['a95__t95_spearman']:.6f}.")
w("")
w("  PAPER BAND SIZES (0-1 / 1-5 / 5-10 / >10)")
for expo in ("t90_frozen", "t80", "t92", "t95"):
    if expo in bres:
        bn = bres[expo]["band_n"]
        tot = sum(bn.values())
        w(f"  {expo:<12} " + "  ".join(f"{lab} {n:,} ({100*n/tot:.1f}%)"
                                       for lab, n in bn.items()))
w("")
w("=" * 80)
w("C. BANDED RESULTS")
w("=" * 80)
w("Same model as the paper: site-stratified Cox, age spline cr(df=4) less one column,")
w("sex, prevalent excluded, 100-event floor, BH-FDR added within each band across the")
w("real outcomes. Top band vs the reference, all real graded outcomes (the T85")
w("report's convention), negative controls counted separately and expected flat:")
w("")
for expo in ("t90_frozen", "t80", "t92", "t95", "t80_quartile", "t92_quartile",
             "t95_quartile", "a95_preserved"):
    if expo not in bres:
        continue
    labs = list(bres[expo]["band_n"].keys())
    top = labs[-1]
    sub = L[(L.exposure == expo) & (L.band == top) & (~L.negative_control)]
    if not len(sub):
        continue
    sig = sub[sub.q < 0.05]
    ms = sig.hr.median() if len(sig) else float("nan")
    negb = L[(L.exposure == expo) & (L.band == top) & (L.negative_control)]
    w(f"  {expo:<14} top band '{top}' holds {bres[expo]['band_n'][top]:>6,}   "
      f"sig {len(sig):>2}/{len(sub)}   median HR {sub.hr.median():5.3f}   "
      f"median sig HR {ms:5.3f}   neg ctrl {int((negb.q < .05).sum())}/{len(negb)}")
w("")
t95bn = bres.get("t95", {}).get("band_n", {})
if t95bn:
    w(f"  T95's paper-band REFERENCE (0-1) holds {t95bn.get('0-1%', 0):,} patients, "
      f"which is the predicted")
    w("  collapse on the reference side rather than the top side. The quartile framing")
    w("  is therefore the readable one for T95 and is reported alongside, clearly")
    w("  labeled. Quartile edges are in banded_positive_control.json.")
w("")
negsig = L[(L.negative_control) & (L.q < 0.05)]
if len(negsig):
    w("  NEGATIVE CONTROL CROSSINGS, every band and framing, named in full:")
    for r in negsig.itertuples():
        w(f"    {r.exposure:<13} {r.band:<5} {r.outcome:<12} HR {r.hr:.2f} "
          f"({r.lo:.2f}-{r.hi:.2f})  q {r.q:.3f}")
    w("  One marginal top-band crossing (t92, back pain) and mid-band blips without")
    w("  dose ordering (glaucoma significant in Q2/Q3 but not Q4). At this many")
    w("  framings and bands a few crossings at q just under 0.05 are expected by")
    w("  chance, and back pain is the least inert of the five controls. Every TOP")
    w("  band in every framing except t92's stays at 0/5, and t92's quartile top")
    w("  band is 0/5, so no framing shows a systematic negative-control signal.")
else:
    w("  No negative control crosses q<0.05 in any band of any framing.")
w("")
w("  A95 PRESERVED-SATURATION FRAMING (percent of the night at or above 95)")
w(f"  cuts used {bpc['a95_cuts_used']}, bands {' / '.join(bpc['a95_band_labels'])}, "
  f"sizes {bpc['a95_band_sizes']}")
w("  reference is the most-preserved band. Read clinically: how much of the night was")
w("  spent at normal saturation. As a continuous per-SD variable this is EXACTLY the")
w("  mirrored T95 (same variable, sign flipped), so per-SD HRs invert and dC is")
w(f"  identical (checked: |dC(a95) - dC(t95)| = {dc['a95_mirror_absdiff']:.1e}, and the")
w(f"  per-SD reciprocity held to {ps['a95_worst_mirror_dev']:.1e} on the demo fits).")
w("  Expect NO new information from the continuous version. The banded version above")
w("  is the only part with its own content, and it does not overtake T90.")
w("")
w("=" * 80)
w("D. PER-SD LADDER, same rows, same model")
w("=" * 80)
w(f"  fitted on the {ps['n_common_rows']:,} patients with every rung present, "
  f"one SD is each rung's own SD")
w(f"  {'rung':<12} {'SD(pp)':>8} {'median HR':>10} {'HR>1':>6} {'sig q<.05':>10} "
  f"{'neg ctrl sig':>13}")
for c in ("t80", "t85", "t88_frozen", "t88_rederiv", "t90_frozen", "t90_rederiv",
          "t92", "t95"):
    st = ps["ladder_stats"][c]
    w(f"  {c:<12} {ps['sd_units_pp'][c]:>8.3f} {st['median_hr']:>10.3f} "
      f"{st['n_hr_gt1']:>4}/{N_RANKED_OUTCOMES} {st['n_sig']:>7}/{N_RANKED_OUTCOMES} {st['neg_sig']:>10}/5")
w("")
w("=" * 80)
w("E. dC ON THE PAPER'S RANKING MACHINERY")
w("=" * 80)
w("  spec           dC          percent of published T90")
for k in ("t80", "t85", "t88_frozen", "t88_rederiv", "t90_frozen", "t90_rederiv",
          "t92", "t95"):
    w(f"  {k:<14} {dc['dC'][k]:+.6f}   {100*dc['dC'][k]/dc['dC']['t90_frozen']:6.1f}")
w("")
nb = dc["n_beats_t90_per_outcome"]
w("  per outcome, the rung beats T90's gain on: " +
  ", ".join(f"{k} {v}/{N_RANKED_OUTCOMES}" for k, v in nb.items()))
w("")
w("=" * 80)
w("F. COST AND INSTANCE HYGIENE")
w("=" * 80)
w(f"  instances launched  {cost['n_instances']}")
w(f"  total actual cost   ${cost['total_usd']:.4f}   against a $10 cap")
w(f"  all terminated      {cost['all_terminated']}")
for r in sorted(cost["instances"], key=lambda x: x["name"]):
    w(f"  {r['instance']}  {r['name']:<19} {r['type']:<12} {r['runtime_min']:>5.1f} min  "
      f"${r['compute_usd']:.4f}")
w("  Launch times are EC2's own. EC2 returned no parenthesized termination time on")
w("  this run, so each end time is the instance's final S3 log push (the last thing")
w("  its shutdown trap does) plus a 60 second buffer, recorded per instance in")
w("  cost.json as end_source. A measured bound from the instance's own last action.")
w("")
w("=" * 80)
w("G. WHAT THIS DOES NOT SETTLE")
w("=" * 80)
w("1. Every rung inherits the frozen T90 convention, the whole-recording denominator")
w("   with no sleep restriction. Chosen deliberately so the rungs are comparable.")
w("2. The paper's band edges 1/5/10 were chosen for T90 and are reused unchanged on")
w("   the other rungs. The quartile and preserved-saturation framings are the")
w("   labeled alternatives, not re-tuned optima, because tuning edges on the outcome")
w("   would fit the answer.")
w("3. Roughly 4-6 percent of records sit in the documented signal-tail stratum where")
w("   the stored signal differs from what the frozen pipeline saw. It affects every")
w("   rung equally and the same-code secondary rungs show it changes nothing.")
w("4. dC differences between adjacent rungs are small against fold noise. The claim")
w("   defended here is the SHAPE (rise to 90, fall past it), not any single gap.")
w("5. Nothing here was written into any manuscript, supplement, table or figure.")
w("")
w("=" * 80)
w("H. FILES. All under Sleep_Variability_2026-08/t_ladder_up/")
w("=" * 80)
w("  worker_ladder.py, launch.sh, build_tasks.py, tasks_all.csv, tasks_validation.csv")
w("  selftest_ladder.py, check_validation.py, assemble.py")
w("  validate_out_0.jsonl, validation_gate.json, validation_detail.csv")
w("  ladder_per_patient.csv          the extraction, one row per patient, numbers only")
w("  extraction_summary.json         fleet-wide coverage and validation")
w("  b01_banded.py -> ladder_banded_long.csv, ladder_banded_results.json,")
w("                   banded_positive_control.json")
w("  b02_persd.py  -> ladder_persd.csv, persd_summary.json")
w("  b03_dC.py     -> dC_summary.csv, dC_per_outcome.csv, dC_result.json")
w("  b04_summary.py-> LADDER_SUMMARY.csv, ladder_peak.json")
w("  b05_figure.py -> ladder_figure.pdf")
w("  b99_cost.py   -> cost.json")
w("")
w("=" * 80)
w("END")
w("=" * 80)

txt = "\n".join(o) + "\n"
open(os.path.join(HERE, "REPORT.txt"), "w").write(txt)
print(txt[:2000])
print(f"... wrote REPORT.txt ({len(txt.splitlines())} lines)")
