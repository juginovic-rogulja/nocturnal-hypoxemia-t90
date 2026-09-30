"""
Assemble the deliverable: one row per outcome, binary T90 above 10% beside the published
continuous T90, at full follow-up and under the two-year landmark.

Interval convention is the published one, copied from
Sleep_Variability_2026-08/dC_per_disease/p02_assemble.py: the plain 2.5th to 97.5th
percentile of the bootstrap distribution of the same quantity the point estimate reports,
500 replicates, patients resampled with replacement inside each hospital, the whole
procedure repeated on each resample.

Because this run uses the published seed base, replicate r here and replicate r of the
published continuous bootstrap are the same resampled patients, so the binary-minus-
continuous difference is paired and gets its own percentile interval.

Carried over, never recomputed, with their source named in the CSV header and in REPORT.md:
  cont_dC / cont_ci_*        Sleep_Variability_2026-08/dC_per_disease/dC_per_disease.csv
  cont_dC_pubbasis           T90_Manuscript/numbers/ranking_v3_percondition.csv
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../..")))  # repository root, where paths.py lives
import paths
import glob
import json
import os
import sys
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import engine_gt10 as X                                                   # noqa: E402
import dc_engine as E                                                     # noqa: E402

T90ROOT = paths.FIGURE_ROOT
PUBDIR = f"{paths.SV_ROOT}/dC_per_disease"
PUBCSV = f"{PUBDIR}/dC_per_disease.csv"
PUBBOOT = {"lag0": f"{PUBDIR}/_work/boot_dc_standard_drop.npy",
           "lag2": f"{PUBDIR}/_work/boot_dc_lag2_drop.npy"}
RANKCSV = f"{paths.NUMBERS_DIR}/ranking_v3_percondition.csv"
sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
from disease_definitions import DISEASES, ORGAN_GROUP, NEGATIVE_CONTROLS  # noqa: E402
import sys as _sys_v8; _sys_v8.path.insert(0, paths.ANALYSIS_DIR)  # v8 sweep 2026-09-12
from cohort_spec import COHORT_N  # v8 sweep 2026-09-12

ALPHA = 0.05
ARMS = {"lag0": {"lag": 0.0, "ck": "_ck_gt10_lag0", "pub_arm": "standard"},
        "lag2": {"lag": 2.0, "ck": "_ck_gt10_lag2", "pub_arm": "lag2"}}
NAME = {k: v[0] for k, v in DISEASES.items()}
NAME["death"] = "Death from any cause"
GROUP = dict(ORGAN_GROUP)
GROUP["death"] = "Mortality"
NEG = set(NEGATIVE_CONTROLS)
SELFTEST = json.load(open(f"{X.WORK}/selftest.json"))
assert SELFTEST["pass"], "the self-test did not pass, nothing may be assembled"


def interval(v):
    """p02_assemble.interval, unchanged: percentile interval plus its Monte Carlo error."""
    v = v[~np.isnan(v)]
    lo, hi = np.percentile(v, [100 * ALPHA / 2, 100 * (1 - ALPHA / 2)])
    sd = float(v.std(ddof=1))
    p = ALPHA / 2
    mc = sd * np.sqrt(p * (1 - p) / len(v)) / stats.norm.pdf(stats.norm.ppf(p))
    return dict(ci_lo=float(lo), ci_hi=float(hi), boot_se=sd, boot_mean=float(v.mean()),
                n_replicates=int(len(v)), mc_error_of_ci_bound=float(mc),
                excludes_zero=bool(lo > 0 or hi < 0))


def load_boot(arm):
    files = sorted(glob.glob(os.path.join(X.WORK, ARMS[arm]["ck"], "boot_*.npz")))
    parts, seeds = [], []
    for f in files:
        z = np.load(f)
        parts.append(z["drop"])
        seeds.append(z["seeds"])
    P = np.concatenate(parts, axis=0)
    s = np.concatenate(seeds)
    assert len(np.unique(s)) == len(s), f"duplicate seeds in arm {arm}"
    assert (np.diff(s) == 1).all() and s[0] == 20260820, "seed order is not the published one"
    return P, s


def dc_matrix(P):
    g = np.stack([P[:, :, 1] - P[:, :, 0], P[:, :, 3] - P[:, :, 2]], axis=2)
    return np.nanmean(g, axis=2)


A = E.build_arrays()
outs = A["outcomes"]
pub = pd.read_csv(PUBCSV)
pub_dis = pub[pub.row_type == "disease"]
rank = pd.read_csv(RANKCSV, comment="#")
rank = rank[rank.feature == "spo2_pct_below_90"].set_index("outcome")

point, boot, cont_boot, cnt = {}, {}, {}, {}
for arm, cfg in ARMS.items():
    point[arm] = np.load(f"{X.WORK}/point_gt10_drop_{arm}.npy")
    M, _s = load_boot(arm)
    boot[arm] = dc_matrix(M)
    np.save(f"{X.WORK}/boot_dc_gt10_{arm}.npy", boot[arm])
    cont_boot[arm] = np.load(PUBBOOT[arm])
    assert cont_boot[arm].shape == boot[arm].shape, (cont_boot[arm].shape, boot[arm].shape)
    cnt[arm] = X.exposed_counts(A, lag=cfg["lag"])
    print(f"{arm}: {boot[arm].shape[0]} binary replicates, "
          f"{cont_boot[arm].shape[0]} published continuous replicates")

point_full = np.load(f"{X.WORK}/point_gt10_full_lag0.npy")
dc_bin = {a: E.dc_from_pass(point[a]) for a in ARMS}
dc_bin_full = E.dc_from_pass(point_full)
est = {a: ~np.isnan(dc_bin[a]) for a in ARMS}

rows = []
for i, o in enumerate(outs):
    r = {"outcome": o, "disease": NAME[o], "organ_group": GROUP[o],
         "negative_control": "yes" if o in NEG else ""}
    for arm in ARMS:
        p = ARMS[arm]["pub_arm"]
        c = cnt[arm].loc[o]
        pr = pub_dis[(pub_dis.arm == p) & (pub_dis.outcome == o)].iloc[0]
        P = point[arm][i]
        pre = f"{arm}_"
        r[pre + "events"] = int(c.events)
        r[pre + "n_at_risk"] = int(c.n_at_risk)
        r[pre + "n_gt10"] = int(c.n_gt10)
        r[pre + "pct_gt10"] = round(float(c.pct_gt10), 3)
        r[pre + "events_gt10"] = int(c.events_gt10)
        r[pre + "events_le10"] = int(c.events_le10)
        r[pre + "events_gt10_I0002"] = int(c.events_gt10_I0002)
        r[pre + "events_gt10_I0006"] = int(c.events_gt10_I0006)
        r[pre + "bin_c_agesex_mean"] = float(np.nanmean([P[0], P[2]]))
        r[pre + "bin_c_agesex_t90gt10_mean"] = float(np.nanmean([P[1], P[3]]))
        r[pre + "bin_dC"] = float(dc_bin[arm][i])
        r[pre + "bin_estimable"] = bool(est[arm][i])
        iv = (interval(boot[arm][:, i]) if est[arm][i] else
              dict(ci_lo=np.nan, ci_hi=np.nan, boot_se=np.nan, boot_mean=np.nan,
                   n_replicates=0, mc_error_of_ci_bound=np.nan, excludes_zero=False))
        for k, v in iv.items():
            r[pre + "bin_" + k] = v
        r[pre + "bin_boot_reps_clearing_the_floor"] = int(np.isfinite(boot[arm][:, i]).sum())
        # the published continuous arm, carried over
        r[pre + "cont_dC"] = float(pr.dC) if np.isfinite(pr.dC) else np.nan
        r[pre + "cont_ci_lo"] = float(pr.ci_lo) if np.isfinite(pr.ci_lo) else np.nan
        r[pre + "cont_ci_hi"] = float(pr.ci_hi) if np.isfinite(pr.ci_hi) else np.nan
        r[pre + "cont_estimable"] = bool(pr.estimable)
        r[pre + "bin_minus_cont"] = (r[pre + "bin_dC"] - r[pre + "cont_dC"]
                                     if est[arm][i] and np.isfinite(pr.dC) else np.nan)
        if est[arm][i] and np.isfinite(pr.dC):
            d = boot[arm][:, i] - cont_boot[arm][:, i]
            dv = interval(d)
            r[pre + "bin_minus_cont_ci_lo"] = dv["ci_lo"]
            r[pre + "bin_minus_cont_ci_hi"] = dv["ci_hi"]
            r[pre + "bin_minus_cont_excludes_zero"] = dv["excludes_zero"]
            r[pre + "bin_retains_pct_of_cont"] = (100 * r[pre + "bin_dC"] / r[pre + "cont_dC"]
                                                  if r[pre + "cont_dC"] != 0 else np.nan)
        else:
            r[pre + "bin_minus_cont_ci_lo"] = np.nan
            r[pre + "bin_minus_cont_ci_hi"] = np.nan
            r[pre + "bin_minus_cont_excludes_zero"] = False
            r[pre + "bin_retains_pct_of_cont"] = np.nan
    r["lag0_bin_dC_pubbasis"] = float(dc_bin_full[i])
    r["lag0_cont_dC_pubbasis"] = float(rank.loc[o, "gain"])
    rows.append(r)

df = pd.DataFrame(rows).sort_values("lag0_bin_dC", ascending=False).reset_index(drop=True)
df.insert(0, "rank_by_binary_gain", range(1, len(df) + 1))

# ------------------------------------------------------------------ means across outcomes
common = est["lag0"] & est["lag2"]
H = {}
for arm in ARMS:
    m = est[arm]
    per_rep_b = np.nanmean(boot[arm][:, m], axis=1)
    per_rep_c = np.nanmean(cont_boot[arm][:, m], axis=1)
    per_rep_d = np.nanmean(boot[arm][:, m] - cont_boot[arm][:, m], axis=1)
    H[arm] = {
        "n_outcomes_in_mean": int(m.sum()),
        "binary_mean_dC": float(np.nanmean(dc_bin[arm][m])),
        "binary_mean_ci": [interval(per_rep_b)["ci_lo"], interval(per_rep_b)["ci_hi"]],
        "continuous_mean_dC": float(np.nanmean(
            pub_dis[pub_dis.arm == ARMS[arm]["pub_arm"]].set_index("outcome")
            .loc[[o for j, o in enumerate(outs) if m[j]], "dC"])),
        "continuous_mean_ci": [interval(per_rep_c)["ci_lo"], interval(per_rep_c)["ci_hi"]],
        "binary_minus_continuous_mean": float(np.nanmean(
            dc_bin[arm][m] - pub_dis[pub_dis.arm == ARMS[arm]["pub_arm"]]
            .set_index("outcome").loc[[o for j, o in enumerate(outs) if m[j]], "dC"].values)),
        "binary_minus_continuous_ci": [interval(per_rep_d)["ci_lo"],
                                       interval(per_rep_d)["ci_hi"]],
        "n_replicates": int(boot[arm].shape[0]),
    }
    H[arm]["binary_retains_pct_of_continuous"] = float(
        100 * H[arm]["binary_mean_dC"] / H[arm]["continuous_mean_dC"])

CSVOUT = f"{HERE}/T90_gt10_concordance_per_disease.csv"
HDR = [
    "# Concordance gain over age and sex for T90 ABOVE 10% of the recording, entered as a raw",
    "#   0/1 indicator, one disease at a time, beside the published continuous T90.",
    f"# built {pd.Timestamp.today().date()}   cohort 19,173   outcomes {len(outs)}   "
    f"above 10% {SELFTEST['n_above_10pct']:,}",
    "# exposure column: spo2_pct_below_90 > 10, from sleep-outcome-sandbox/feature_pipeline/",
    "#   data_frozen_v8_2026-09/master_cohort.csv, in percent of the recording",   # v8 repoint 2026-09-12 (report text)
    "# procedure: unchanged from T90_Manuscript/numbers/build_ranking_v3.py by way of",
    "#   Sleep_Variability_2026-08/dC_per_disease/dc_engine.py. Age spline cr(df=4) with one",
    "#   column dropped plus sex, CoxPHFitter(penalizer=0.01) fitted at I0002 and scored at",
    "#   I0006 then reversed, the two held-out concordances averaged, prevalent disease",
    "#   excluded, follow-up above zero, 30-event floor per fold.",
    "# lag0 = full follow-up. lag2 = two-year landmark, the published rule: every patient",
    "#   whose event or censoring arrived within 2 years leaves and the clock restarts there",
    "#   (numbers/lag_ladder.csv, reproduced exactly by this run's self-test).",
    "# 95% CI: 2.5th to 97.5th percentile of 500 bootstrap replicates, patients resampled",
    "#   with replacement inside each hospital, seeds 20260820..20261319, the published seeds.",
    "# cont_* columns are CARRIED OVER, not recomputed:",
    "#   lag0_cont_* and lag2_cont_* from Sleep_Variability_2026-08/dC_per_disease/"
    "dC_per_disease.csv",
    "#   lag0_cont_dC_pubbasis from T90_Manuscript/numbers/ranking_v3_percondition.csv",
    "# lag0_bin_dC_pubbasis is this run's binary gain on the published four-column age basis,",
    "#   the like-for-like partner of lag0_cont_dC_pubbasis.",
    "# readers must pass comment='#' to pandas.read_csv.",
]
with open(CSVOUT, "w") as fh:
    fh.write("\n".join(HDR) + "\n")
    df.to_csv(fh, index=False)
print(f"wrote {CSVOUT}  ({len(df)} rows, {len(df.columns)} columns)")

flip = df[(df.lag0_bin_dC > df.lag0_cont_dC) & df.lag0_bin_estimable]
flip2 = df[(df.lag2_bin_dC > df.lag2_cont_dC) & df.lag2_bin_estimable]
J = {
    "built": str(pd.Timestamp.today().date()),
    "question": "gain in held-out concordance over age and sex when the exposure is T90 above "
                "10% of the recording, as a binary indicator, one disease at a time",
    "exposure_column": "spo2_pct_below_90",
    "exposure_source": f"{paths.TABLES_DIR}/"
                       "master_cohort.csv",
    "exposure_rule": "1 when spo2_pct_below_90 > 10, else 0, entered raw",
    "cohort_n": COHORT_N,
    "n_above_10pct": SELFTEST["n_above_10pct"],
    "pct_above_10pct": round(100 * SELFTEST["n_above_10pct"] / COHORT_N, 2),
    "n_outcomes": len(outs),
    "bootstrap_replicates": int(boot["lag0"].shape[0]),
    "full_followup": H["lag0"],
    "landmark_2y": H["lag2"],
    "top5_binary_full_followup": [
        {"disease": r.disease, "dC": r.lag0_bin_dC,
         "ci": [r.lag0_bin_ci_lo, r.lag0_bin_ci_hi],
         "continuous_dC": r.lag0_cont_dC, "n_above_10pct": int(r.lag0_n_gt10),
         "events": int(r.lag0_events)}
        for _i, r in df.head(5).iterrows()],
    "outcomes_where_binary_beats_continuous_full_followup":
        [{"disease": r.disease, "binary_dC": r.lag0_bin_dC, "continuous_dC": r.lag0_cont_dC,
          "difference": r.lag0_bin_minus_cont,
          "difference_ci": [r.lag0_bin_minus_cont_ci_lo, r.lag0_bin_minus_cont_ci_hi],
          "difference_excludes_zero": bool(r.lag0_bin_minus_cont_excludes_zero)}
         for _i, r in flip.iterrows()],
    "outcomes_where_binary_beats_continuous_landmark":
        [{"disease": r.disease, "binary_dC": r.lag2_bin_dC, "continuous_dC": r.lag2_cont_dC,
          "difference": r.lag2_bin_minus_cont,
          "difference_ci": [r.lag2_bin_minus_cont_ci_lo, r.lag2_bin_minus_cont_ci_hi],
          "difference_excludes_zero": bool(r.lag2_bin_minus_cont_excludes_zero)}
         for _i, r in flip2.iterrows()],
    "not_estimable": {a: [outs[i] for i in range(len(outs)) if not est[a][i]] for a in ARMS},
    "published_basis_cross_reference": {
        "binary_mean_dC_four_column_basis": float(np.nanmean(dc_bin_full)),
        "continuous_mean_dC_four_column_basis": 0.02033875497249791,
        "continuous_source": RANKCSV},
    "selftest": {k: SELFTEST[k] for k in SELFTEST if k.startswith("gate") or k == "pass"},
}
JSONOUT = f"{HERE}/T90_gt10_headline.json"
json.dump(J, open(JSONOUT, "w"), indent=2, default=float)
print(f"wrote {JSONOUT}")

print(f"\nfull follow-up   binary mean dC {H['lag0']['binary_mean_dC']:+.6f} "
      f"({H['lag0']['binary_mean_ci'][0]:+.6f} to {H['lag0']['binary_mean_ci'][1]:+.6f})   "
      f"continuous {H['lag0']['continuous_mean_dC']:+.6f}")
print(f"2-year landmark  binary mean dC {H['lag2']['binary_mean_dC']:+.6f} "
      f"({H['lag2']['binary_mean_ci'][0]:+.6f} to {H['lag2']['binary_mean_ci'][1]:+.6f})   "
      f"continuous {H['lag2']['continuous_mean_dC']:+.6f}")
print(f"\nbinary above continuous, full follow-up: "
      f"{list(flip.disease) if len(flip) else 'none'}")
print(f"binary above continuous, landmark: {list(flip2.disease) if len(flip2) else 'none'}")
print("\ntop 12 by binary gain, full follow-up")
print(df.head(12)[["disease", "lag0_bin_dC", "lag0_bin_ci_lo", "lag0_bin_ci_hi",
                   "lag0_cont_dC", "lag0_n_gt10", "lag0_events"]].to_string(index=False))
