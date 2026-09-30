"""
What moved between ranking_v2 (2026-08-02, 19,383 patients) and ranking_v3 (2026-08-07, 19,173).

Produces, from the two frozen files and nothing typed by hand:
  1. old and new rank and gain for every measure the manuscript names
  2. every measure whose rank moved by more than 5 places
  3. whether T90 is now first, tested rather than asserted, by pairing T90 against the corrected
     nadir on the same 48 outcomes and by resampling the outcomes
  4. how the negative-control outcomes behave under the top measures
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../..")))  # repository root, where paths.py lives
import paths
import warnings; warnings.filterwarnings("ignore")
import json
import numpy as np, pandas as pd
from scipy import stats
import os as _os_v8, sys as _sys_v8; _sys_v8.path.insert(0, _os_v8.path.dirname(_os_v8.path.abspath(__file__)))  # v8 sweep 2026-09-12
from cohort_spec import COHORT_N, N_RANKED_OUTCOMES  # v8 sweep 2026-09-12
from disease_definitions import NEGATIVE_CONTROLS, RANKING_EXCLUDE, SWAPPED_CONTROLS_2026_09_14   # v8.1: controls from the definitions

ROOT = paths.T90_ROOT
N = f"{paths.NUMBERS_DIR}"

NAMED = [
    ("spo2_pct_below_90", "T90, recording below 90%"),
    ("spo2_nadir_corrected", "corrected nadir saturation"),
    ("odi3_total", "ODI3"),
    ("odi4_total", "ODI4"),
    ("spo2_pct_below_88", "T88, recording below 88%"),
    ("spo2_mean", "mean SpO2"),
    ("nrem_hr_bpm", "non-REM heart rate"),
    ("wholenight_hr_bpm", "whole-night heart rate"),
    ("wholenight_lf_hf_ratio", "whole-night LF/HF ratio"),
    ("nrem_lf_hf_ratio", "non-REM LF/HF ratio"),
    ("AHI", "apnea-hypopnea index"),
    ("arousal_index", "arousal index"),
    ("N3_pct", "slow-wave sleep, percent"),
    ("N3_min", "slow-wave sleep, minutes"),
    ("sleep_efficiency_pct", "sleep efficiency"),
    ("TST_min", "total sleep time"),
]
NEGCTRL = [k for k in NEGATIVE_CONTROLS if k not in RANKING_EXCLUDE]      # v8.1: the settled controls that are also ranked (glaucoma, inguinal hernia, alopecia)
DROPPED_CTRL = ["fracture", "osteoarthritis"] + sorted(SWAPPED_CONTROLS_2026_09_14)   # scored as outcomes in v2, no longer controls (fracture, osteoarthritis, and the 14 Sept swap-outs)

v2 = pd.read_csv(f"{N}/ranking_v2.csv")
v3 = pd.read_csv(f"{N}/ranking_v3.csv", comment="#")
p2 = pd.read_csv(f"{N}/ranking_v2_percondition.csv")
p3 = pd.read_csv(f"{N}/ranking_v3_percondition.csv", comment="#")
# v8: the measure list is 141 (60 per-site twins dropped, 4 measures added) and the outcome list 52 (4 added), so the
# v2 comparison is made on the common measures and outcomes; the delta is printed and recorded, never hidden.
_common_f = sorted(set(v2.feature) & set(v3.feature)); _common_o = sorted(set(p2.outcome) & set(p3.outcome))
V8_SET_DELTA = {"measures_only_in_v2": sorted(set(v2.feature) - set(v3.feature)), "measures_only_in_v3": sorted(set(v3.feature) - set(v2.feature)),
                "outcomes_only_in_v2": sorted(set(p2.outcome) - set(p3.outcome)), "outcomes_only_in_v3": sorted(set(p3.outcome) - set(p2.outcome)),
                "common_measures": len(_common_f), "common_outcomes": len(_common_o)}
print("v8 set delta:", {k: (v if isinstance(v, int) else len(v)) for k, v in V8_SET_DELTA.items()})
v2 = v2[v2.feature.isin(_common_f)].copy(); v3 = v3[v3.feature.isin(_common_f)].copy()
p2 = p2[p2.outcome.isin(_common_o)].copy(); p3 = p3[p3.outcome.isin(_common_o)].copy()
assert set(v2.feature) == set(v3.feature) and set(p2.outcome) == set(p3.outcome), "restriction to the common sets failed"

m = v2.merge(v3, on="feature", suffixes=("_old", "_new"))
m["d_rank"] = m["rank_new"] - m["rank_old"]          # negative means the measure moved up
m["d_dC"] = m["dC_new"] - m["dC_old"]
# family labels come from the hand-classified lookup, never the retired keyword classifier
fam = pd.read_csv(f"{N}/measure_families.csv")[["feature", "family", "definition"]]
m = m.merge(fam, on="feature", how="left").sort_values("rank_new")
m.to_csv(f"{N}/ranking_v2_to_v3_moves.csv", index=False)

R = {"v8_set_delta": V8_SET_DELTA, "cohort_old": 19383, "cohort_new": COHORT_N,
     "baseline_heldout_C_old": float((v2.mC - v2.dC).mean()),
     "baseline_heldout_C_new": float((v3.mC - v3.dC).mean())}

# The retired root copy was an unlabelled 2026-08-04 rebuild on this same 19,173 cohort. It is an
# independent second run of the same procedure, so v3 should land on top of it.
RET = (f"{paths.FIGURE_ROOT}/New_Figures/_NOT_THESE_workfiles/RETIRED_root_ranking_v2_2026-08-04.csv")
try:
    if __import__("os").stat(RET).st_blocks == 0: raise FileNotFoundError(f"{RET} is iCloud-dataless (0 blocks), the read would hang")   # v8 guard
    ret = pd.read_csv(RET)
    c = v3.merge(ret, on="feature", suffixes=("_v3", "_ret"))
    dd = (c.dC_v3 - c.dC_ret).abs()
    R["retired_root_copy_check"] = {
        "max_abs_dC_difference": float(dd.max()),
        "n_measures_matching_to_1e9": int((dd < 1e-9).sum()), "n_measures": int(len(c)),
        "n_ranks_identical": int((c["rank_v3"] == c["rank_ret"]).sum())}
    print(f"cross-check against the retired 2026-08-04 rebuild on the same cohort: "
          f"max |dC difference| {dd.max():.2e}, ranks identical for "
          f"{int((c['rank_v3'] == c['rank_ret']).sum())} of {len(c)} measures")
except FileNotFoundError:
    pass

print("=" * 96)
print("1. THE MEASURES THE MANUSCRIPT NAMES")
print("=" * 96)
print(f"{'measure':<26}{'':<2}{'OLD rank':>9}{'NEW rank':>9}{'move':>7}   "
      f"{'OLD dC':>11}{'NEW dC':>11}{'change':>11}")
named_rows = []
for f, label in NAMED:
    r = m[m.feature == f]
    if not len(r):
        print(f"{label:<28}NOT IN THE RANKING"); continue
    r = r.iloc[0]
    mv = int(r.d_rank)
    print(f"{label:<28}{int(r['rank_old']):>9}{int(r['rank_new']):>9}"
          f"{('+' if mv > 0 else '')+str(mv):>7}   "
          f"{r.dC_old:>+11.6f}{r.dC_new:>+11.6f}{r.d_dC:>+11.6f}")
    named_rows.append(dict(feature=f, label=label, rank_old=int(r["rank_old"]),
                           rank_new=int(r["rank_new"]), move=mv,
                           dC_old=float(r.dC_old), dC_new=float(r.dC_new),
                           dC_change=float(r.d_dC)))
R["named_measures"] = named_rows
pd.DataFrame(named_rows).to_csv(f"{N}/ranking_v3_named_measures.csv", index=False)

print()
print("=" * 96)
print("2. IS T90 NOW FIRST?")
print("=" * 96)
top = m.nsmallest(6, "rank_new")[["feature", "rank_old", "rank_new", "dC_old", "dC_new"]]
print(top.to_string(index=False))

A, B = "spo2_pct_below_90", "spo2_nadir_corrected"


def paired(per, a, b, tag):
    x = per[per.feature == a].set_index("outcome")["gain"]
    y = per[per.feature == b].set_index("outcome")["gain"]
    common = sorted(set(x.index) & set(y.index))
    dif = (x[common] - y[common]).astype(float)
    se = dif.std(ddof=1) / np.sqrt(len(dif))
    t, p = stats.ttest_rel(x[common], y[common])
    w = stats.wilcoxon(x[common], y[common])
    rng = np.random.default_rng(20260807)
    idx = rng.integers(0, len(dif), size=(20000, len(dif)))
    bs = dif.values[idx].mean(axis=1)
    out = dict(tag=tag, n_outcomes=len(dif), mean_diff=float(dif.mean()), se=float(se),
               ci=[float(dif.mean() - 1.96 * se), float(dif.mean() + 1.96 * se)],
               t=float(t), p_paired_t=float(p), p_wilcoxon=float(w.pvalue),
               n_outcomes_A_beats_B=int((dif > 0).sum()),
               boot_p_A_ahead=float((bs > 0).mean()),
               boot_ci=[float(np.quantile(bs, 0.025)), float(np.quantile(bs, 0.975))])
    print(f"\n  {tag}: T90 minus corrected nadir, paired over {out['n_outcomes']} outcomes")
    print(f"    mean difference {out['mean_diff']:+.6f}  95% CI "
          f"{out['ci'][0]:+.6f} to {out['ci'][1]:+.6f}")
    print(f"    paired t p = {out['p_paired_t']:.3f}   Wilcoxon p = {out['p_wilcoxon']:.3f}   "
          f"T90 ahead on {out['n_outcomes_A_beats_B']} of {out['n_outcomes']}")
    print(f"    outcome bootstrap: T90 ahead in {100*out['boot_p_A_ahead']:.1f}% of resamples")
    return out


R["t90_vs_nadir_v2"] = paired(p2, A, B, "on v2, 19,383")
R["t90_vs_nadir_v3"] = paired(p3, A, B, "on v3, 19,173")


def rank_bootstrap(per, tag, nboot=4000):
    """Resample the 48 outcomes and re-rank, the same way the paper's 12th-to-137th was made."""
    piv = per.pivot_table(index="outcome", columns="feature", values="gain")
    outs = piv.index.to_numpy()
    rng = np.random.default_rng(20260807)
    feats = piv.columns.to_numpy()
    ranks = np.zeros((nboot, len(feats)), dtype=np.int32)
    vals = piv.to_numpy()
    for i in range(nboot):
        s = rng.integers(0, len(outs), len(outs))
        mean = np.nanmean(vals[s], axis=0)
        order = np.argsort(-mean)
        rk = np.empty(len(feats), dtype=np.int32)
        rk[order] = np.arange(1, len(feats) + 1)
        ranks[i] = rk
    res = {}
    for f in ("spo2_pct_below_90", "spo2_nadir_corrected", "odi3_total", "spo2_pct_below_88",
              "AHI", "arousal_index", "sleep_efficiency_pct", "TST_min", "N3_pct"):
        j = int(np.where(feats == f)[0][0])
        res[f] = dict(p_first=float((ranks[:, j] == 1).mean()),
                      lo=int(np.quantile(ranks[:, j], 0.025)),
                      hi=int(np.quantile(ranks[:, j], 0.975)),
                      median=int(np.median(ranks[:, j])))
    print(f"\n  {tag}: rank under outcome resampling, {nboot} resamples")
    for f, v in res.items():
        print(f"    {f:<24} median {v['median']:>4}   95% range {v['lo']} to {v['hi']}   "
              f"first in {100*v['p_first']:.1f}%")
    return res


R["rank_bootstrap_v2"] = rank_bootstrap(p2, "v2, 19,383")
R["rank_bootstrap_v3"] = rank_bootstrap(p3, "v3, 19,173")

print()
print("=" * 96)
print("3. MEASURES WHOSE RANK MOVED BY MORE THAN 5 PLACES")
print("=" * 96)
big = m[m.d_rank.abs() > 5].sort_values("d_rank")
print(f"{len(big)} of {len(m)} measures moved more than 5 places")
print(f"{'measure':<30}{'OLD':>5}{'NEW':>5}{'move':>7}{'OLD dC':>12}{'NEW dC':>12}{'change':>12}")
for _, r in big.iterrows():
    mv = int(r.d_rank)
    print(f"{r.feature:<30}{int(r['rank_old']):>5}{int(r['rank_new']):>5}"
          f"{('+' if mv > 0 else '')+str(mv):>7}"
          f"{r.dC_old:>+12.6f}{r.dC_new:>+12.6f}{r.d_dC:>+12.6f}")
big.to_csv(f"{N}/ranking_v3_big_movers.csv", index=False)
R["n_moved_gt5"] = int(len(big))
R["moved_gt5"] = [dict(feature=r.feature, rank_old=int(r["rank_old"]),
                       rank_new=int(r["rank_new"]), move=int(r.d_rank),
                       dC_old=float(r.dC_old), dC_new=float(r.dC_new))
                  for _, r in big.iterrows()]

print()
print("=" * 96)
print("3b. DOES THE BEST MEASURE OF ANY FAMILY CHANGE?")
print("=" * 96)
print("   run_combo_v2.py builds the combination figure from the top measure of each family,")
print("   so a family whose winner changes changes that figure, whatever the rank moves say.")
fb = []
for family, grp in m.groupby("family"):
    o = grp.loc[grp["rank_old"].idxmin()]
    n_ = grp.loc[grp["rank_new"].idxmin()]
    fb.append(dict(family=family, n_measures=len(grp), best_v2=o.feature, best_v3=n_.feature,
                   changed=bool(o.feature != n_.feature),
                   dC_v2_of_best_v2=float(o.dC_old), dC_v3_of_best_v3=float(n_.dC_new)))
fbd = pd.DataFrame(fb).sort_values("dC_v3_of_best_v3", ascending=False)
print(fbd.to_string(index=False))
fbd.to_csv(f"{N}/ranking_v3_family_best.csv", index=False)
R["family_best"] = fb
R["n_families_whose_best_changed"] = int(fbd.changed.sum())

print()
print("=" * 96)
print("4. NEGATIVE CONTROLS")
print("=" * 96)
ctrl = {}
for o in NEGCTRL + DROPPED_CTRL:
    row = {"is_settled_control": o in NEGCTRL}
    for f, label in NAMED[:8]:
        g2 = p2[(p2.feature == f) & (p2.outcome == o)]["gain"]
        g3 = p3[(p3.feature == f) & (p3.outcome == o)]["gain"]
        if len(g3):
            row[f] = {"v2": float(g2.iloc[0]) if len(g2) else None, "v3": float(g3.iloc[0])}
    ctrl[o] = row
    tag = "control" if o in NEGCTRL else "no longer a control"
    print(f"\n  {o}  ({tag})")
    for f, label in NAMED[:8]:
        if f in row:
            print(f"    {label:<28} v2 {row[f]['v2']:+.6f}   v3 {row[f]['v3']:+.6f}")
R["negative_controls"] = ctrl

# T90's gain on the controls against its gain on everything else
t3 = p3[p3.feature == A].set_index("outcome")["gain"]
inctrl = t3[[o for o in NEGCTRL if o in t3.index]]
rest = t3[[o for o in t3.index if o not in NEGCTRL + DROPPED_CTRL]]
print(f"\n  T90 mean gain on the {len(inctrl)} settled controls in this outcome set  {inctrl.mean():+.6f}")
print(f"  T90 mean gain on the other {len(rest)} outcomes                  {rest.mean():+.6f}")
R["t90_gain_controls_v3"] = float(inctrl.mean())
R["t90_gain_noncontrols_v3"] = float(rest.mean())

try:
    nf = json.load(open(f"{N}/ranking_v3_noisefloor.json"))
    R["noise_floor"] = nf
    print()
    print("=" * 96)
    print("5. WHAT A USELESS VARIABLE EARNS FOR FREE")
    print("=" * 96)
    print(f"  {nf['n_noise_variables']} random variables, same cohort, same {N_RANKED_OUTCOMES} outcomes")
    print(f"  mean gain per useless variable   {nf['dC_per_noise_variable_mean']:+.6f}")
    print(f"  spread across them, SD           {nf['dC_per_noise_variable_sd']:.6f}")
    print(f"  range                            {nf['dC_per_noise_variable_min']:+.6f} to "
          f"{nf['dC_per_noise_variable_max']:+.6f}")
    print(f"  a useless variable beats the baseline on {nf['npos_mean_of_48']:.1f} of {N_RANKED_OUTCOMES} outcomes")
    nfloor = nf["dC_per_noise_variable_max"]
    below = m[m.dC_new <= nfloor]
    print(f"  measures in v3 not beating the best of {nf['n_noise_variables']} random variables: "
          f"{len(below)} of {len(m)}, from rank {int(below['rank_new'].min())} down")
    R["n_measures_below_noise_max"] = int(len(below))
    R["noise_max_dC"] = float(nfloor)
except FileNotFoundError:
    print("\n  (noise floor not built yet, run ranking_v3_noisefloor.py)")

json.dump(R, open(f"{N}/ranking_v3_report.json", "w"), indent=2)
print(f"\nwritten: {N}/ranking_v3_report.json, ranking_v2_to_v3_moves.csv, "
      f"ranking_v3_named_measures.csv, ranking_v3_big_movers.csv")
