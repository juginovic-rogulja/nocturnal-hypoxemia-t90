"""
Assemble the per-disease table from the point estimates and the bootstrap checkpoints, for both
arms, and add the one quantity that IS an average: the mean over the 48 outcomes.

Writes
  dC_per_disease.csv          long, one row per disease per arm, plus the summary rows
  dC_per_disease_TABLE.md     the ranked eTable, both arms side by side
  _work/boot_dc_<arm>_<basis>.npy   the raw (n_rep, 48) bootstrap dC matrices

Confidence interval is the plain percentile interval, the 2.5th and 97.5th percentiles of the
bootstrap distribution of the same quantity the point estimate reports. No bias correction,
because the quantity is a difference of two concordances on the same held-out rows and the
bootstrap distribution is close to symmetric.

The mean across outcomes. The 48 gains are not independent. They are 48 views of one set of
patients, so a disease-level resample would be the wrong unit and averaging 48 separate
intervals would be wrong twice over. The patient resample already carries the dependence: each
replicate recomputes all 48 gains on the same resampled patients, their mean is taken inside
the replicate, and the interval is the percentiles of those means.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
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

import dc_engine as E

T90ROOT = paths.T90_ROOT
sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
from disease_definitions import DISEASES, ORGAN_GROUP, NEGATIVE_CONTROLS   # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ALPHA = 0.05
ARMS = {"standard": {"lag": 0.0, "ck": "_ck", "bases": ("drop", "full")},
        "lag2": {"lag": 2.0, "ck": "_ck_lag2", "bases": ("drop",)}}

NAME = {k: v[0] for k, v in DISEASES.items()}
NAME["death"] = "Death from any cause"
GROUP = dict(ORGAN_GROUP)
GROUP["death"] = "Mortality"
NEG = set(NEGATIVE_CONTROLS)


def load_boot(arm, basis):
    files = sorted(glob.glob(os.path.join(E.WORK, ARMS[arm]["ck"], "boot_*.npz")))
    if not files:
        return None, None
    parts, seeds = [], []
    for f in files:
        z = np.load(f)
        parts.append(z[basis])
        seeds.append(z["seeds"])
    P = np.concatenate(parts, axis=0)
    s = np.concatenate(seeds)
    assert len(np.unique(s)) == len(s), f"duplicate bootstrap seeds in arm {arm}"
    return P, s


def dc_matrix(P):
    """(n_rep, 48) averaged dC, the paper's nanmean over the two cross-fit directions."""
    g = np.stack([P[:, :, 1] - P[:, :, 0], P[:, :, 3] - P[:, :, 2]], axis=2)
    return np.nanmean(g, axis=2)


def interval(v):
    """Percentile interval, standard error and the Monte Carlo error of a bound."""
    v = v[~np.isnan(v)]
    lo, hi = np.percentile(v, [100 * ALPHA / 2, 100 * (1 - ALPHA / 2)])
    sd = float(v.std(ddof=1))
    p = ALPHA / 2
    mc = sd * np.sqrt(p * (1 - p) / len(v)) / stats.norm.pdf(stats.norm.ppf(p))
    return dict(ci_lo=float(lo), ci_hi=float(hi), boot_se=sd, boot_mean=float(v.mean()),
                n_replicates=int(len(v)), mc_error_of_ci_bound=float(mc),
                excludes_zero=bool(lo > 0 or hi < 0))


def main():
    A = E.build_arrays()
    outs = A["outcomes"]
    idx = {o: i for i, o in enumerate(outs)}

    pt, dcb, cnt, nrep, folds = {}, {}, {}, {}, {}
    for arm, cfg in ARMS.items():
        f = os.path.join(E.WORK, f"point_{arm}_drop.npy")
        if not os.path.exists(f):
            print(f"arm {arm}: no point estimate yet, skipped")
            continue
        got = {}
        for b in cfg["bases"]:
            M, _s = load_boot(arm, b)
            if M is None:
                break
            got[b] = (np.load(os.path.join(E.WORK, f"point_{arm}_{b}.npy")), M)
        if len(got) != len(cfg["bases"]):
            print(f"arm {arm}: checkpoints incomplete, arm skipped")
            continue
        pt[arm] = {}
        for b, (Pp, M) in got.items():
            pt[arm][b] = Pp
            # a direction is lost when either of its two concordances is missing, which the
            # 30-event fold floor causes. One direction lost still yields a dC, both lost
            # leaves the replicate with no estimate for that outcome.
            d1 = np.isnan(M[:, :, 0]) | np.isnan(M[:, :, 1])
            d2 = np.isnan(M[:, :, 2]) | np.isnan(M[:, :, 3])
            folds[(arm, b)] = {"one_direction": (d1 ^ d2).sum(axis=0),
                               "both_directions": (d1 & d2).sum(axis=0)}
            dcb[(arm, b)] = dc_matrix(M)
            np.save(os.path.join(E.WORK, f"boot_dc_{arm}_{b}.npy"), dcb[(arm, b)])
            nrep[(arm, b)] = M.shape[0]
        cnt[arm] = E.counts(A, lag=cfg["lag"]).set_index("outcome")
        print(f"arm {arm}: {nrep.get((arm, 'drop'), 0)} replicates on disk")

    arms = [a for a in ARMS if a in pt]
    point_dc = {a: E.dc_from_pass(pt[a]["drop"]) for a in arms}
    estimable = {a: ~np.isnan(point_dc[a]) for a in arms}
    common = np.logical_and.reduce([estimable[a] for a in arms])

    rows = []
    for arm in arms:
        for i, o in enumerate(outs):
            P = pt[arm]["drop"][i]
            r = {"row_type": "disease", "arm": arm, "landmark_years": ARMS[arm]["lag"],
                 "outcome": o, "disease": NAME[o], "organ_group": GROUP[o],
                 "negative_control": "yes" if o in NEG else "",
                 "events": int(cnt[arm].loc[o, "events"]),
                 "n_at_risk": int(cnt[arm].loc[o, "n_at_risk"]),
                 "events_I0002": int(cnt[arm].loc[o, "events_I0002"]),
                 "events_I0006": int(cnt[arm].loc[o, "events_I0006"]),
                 "n_I0002": int(cnt[arm].loc[o, "n_I0002"]),
                 "n_I0006": int(cnt[arm].loc[o, "n_I0006"]),
                 "c_agesex_dir1_fitI0002_scoreI0006": P[0],
                 "c_agesex_t90_dir1_fitI0002_scoreI0006": P[1],
                 "dC_dir1": P[1] - P[0],
                 "c_agesex_dir2_fitI0006_scoreI0002": P[2],
                 "c_agesex_t90_dir2_fitI0006_scoreI0002": P[3],
                 "dC_dir2": P[3] - P[2],
                 "c_agesex_mean": np.nanmean([P[0], P[2]]),
                 "c_agesex_t90_mean": np.nanmean([P[1], P[3]]),
                 "dC": float(point_dc[arm][i]),
                 "estimable": bool(estimable[arm][i])}
            # When the observed data cannot be fitted, the replicates that happen to clear
            # the fold floor are a selected subset of the bootstrap distribution, so no
            # interval is reported for that outcome. The count is kept as a diagnostic.
            v = dcb[(arm, "drop")][:, i]
            r.update(interval(v) if estimable[arm][i]
                     else dict(ci_lo=np.nan, ci_hi=np.nan, boot_se=np.nan, boot_mean=np.nan,
                               n_replicates=0, mc_error_of_ci_bound=np.nan,
                               excludes_zero=False))
            r["boot_reps_clearing_the_floor"] = int(np.isfinite(v).sum())
            r["boot_reps_losing_one_direction"] = int(folds[(arm, "drop")]["one_direction"][i])
            r["boot_reps_with_no_estimate"] = int(folds[(arm, "drop")]["both_directions"][i])
            if arm == "standard":
                Q = pt[arm]["full"][i]
                r["c_agesex_mean_pubbasis"] = np.nanmean([Q[0], Q[2]])
                r["c_agesex_t90_mean_pubbasis"] = np.nanmean([Q[1], Q[3]])
                r["dC_pubbasis"] = float(np.nanmean([Q[1] - Q[0], Q[3] - Q[2]]))
                for k, val in interval(dcb[("standard", "full")][:, i]).items():
                    r[k + "_pubbasis"] = val
            rows.append(r)

    # ------------------------------------------------------------ the mean across outcomes
    for arm in arms:
        for tag, mask, label in (
                ("__mean_all__", estimable[arm],
                 "Mean across every outcome estimable in that arm"),
                ("__mean_common__", common,
                 f"Mean across the {int(common.sum())} outcomes estimable in both arms")):
            M = dcb[(arm, "drop")][:, mask]
            per_rep = np.nanmean(M, axis=1)
            r = {"row_type": "mean_of_outcomes", "arm": arm,
                 "landmark_years": ARMS[arm]["lag"], "outcome": tag, "disease": label,
                 "organ_group": "", "negative_control": "",
                 "n_outcomes_in_mean": int(mask.sum()),
                 "min_outcomes_contributing_per_replicate":
                     int(np.isfinite(M).sum(axis=1).min()),
                 "events": int(cnt[arm].loc[[o for j, o in enumerate(outs) if mask[j]],
                                            "events"].sum()),
                 "dC": float(np.nanmean(point_dc[arm][mask]))}
            r.update(interval(per_rep))
            if arm == "standard":
                Mf = dcb[("standard", "full")][:, mask]
                r["dC_pubbasis"] = float(np.nanmean(
                    E.dc_from_pass(pt["standard"]["full"])[mask]))
                for k, val in interval(np.nanmean(Mf, axis=1)).items():
                    r[k + "_pubbasis"] = val
            rows.append(r)

    df = pd.DataFrame(rows)
    dis = df[df.row_type == "disease"].copy()
    order = (dis[dis.arm == "standard"].sort_values("dC", ascending=False).outcome.tolist())
    dis["_o"] = dis.outcome.map({o: i for i, o in enumerate(order)})
    dis = dis.sort_values(["_o", "arm"]).drop(columns="_o")
    df = pd.concat([dis, df[df.row_type != "disease"]], ignore_index=True)
    csv = os.path.join(HERE, "dC_per_disease.csv")
    df.to_csv(csv, index=False)
    print(f"wrote {csv}  ({len(df)} rows)")

    # ------------------------------------------------------------ the ranked eTable
    S = dis[dis.arm == "standard"].set_index("outcome")
    L = dis[dis.arm == "lag2"].set_index("outcome") if "lag2" in arms else None

    def fmt(r, tag=""):
        if not np.isfinite(r["dC" + tag]) or not np.isfinite(r["ci_lo" + tag]):
            return "not estimable"
        return f"{r['dC' + tag]:.4f} ({r['ci_lo' + tag]:.4f} to {r['ci_hi' + tag]:.4f})"

    body = {"Disease": [], "Organ-system group": [], "Incident cases": [], "Baseline C": [],
            "C with T90": [], "dC (95% CI)": [], "Cases after 2-y landmark": [],
            "dC at 2-y landmark (95% CI)": [], "Negative control": []}
    for o in order:
        s = S.loc[o]
        body["Disease"].append(s.disease)
        body["Organ-system group"].append(s.organ_group)
        body["Incident cases"].append(f"{int(s.events):,}")
        body["Baseline C"].append(f"{s.c_agesex_mean:.3f}")
        body["C with T90"].append(f"{s.c_agesex_t90_mean:.3f}")
        body["dC (95% CI)"].append(fmt(s))
        if L is not None:
            body["Cases after 2-y landmark"].append(f"{int(L.loc[o].events):,}")
            body["dC at 2-y landmark (95% CI)"].append(fmt(L.loc[o]))
        else:
            body["Cases after 2-y landmark"].append("")
            body["dC at 2-y landmark (95% CI)"].append("")
        body["Negative control"].append(s.negative_control)
    for tag in ("__mean_all__", "__mean_common__"):
        ms = df[(df.outcome == tag) & (df.arm == "standard")].iloc[0]
        body["Disease"].append(ms.disease)
        body["Organ-system group"].append("")
        body["Incident cases"].append("")
        body["Baseline C"].append("")
        body["C with T90"].append("")
        body["dC (95% CI)"].append(fmt(ms))
        if L is not None:
            ml = df[(df.outcome == tag) & (df.arm == "lag2")].iloc[0]
            body["Cases after 2-y landmark"].append("")
            body["dC at 2-y landmark (95% CI)"].append(fmt(ml))
        else:
            body["Cases after 2-y landmark"].append("")
            body["dC at 2-y landmark (95% CI)"].append("")
        body["Negative control"].append("")
    tab = pd.DataFrame(body)

    dp = dis[dis.arm == "standard"].sort_values("dC_pubbasis", ascending=False)
    tab2 = pd.DataFrame({
        "Disease": dp.disease,
        "Incident cases": dp.events.map("{:,}".format),
        "Baseline C": dp.c_agesex_mean_pubbasis.map("{:.3f}".format),
        "C with T90": dp.c_agesex_t90_mean_pubbasis.map("{:.3f}".format),
        "dC (95% CI)": [fmt(r, "_pubbasis") for _i, r in dp.iterrows()],
    })

    nex = int(S.excludes_zero.astype(bool).sum())
    nexl = int(L.excludes_zero.astype(bool).sum()) if L is not None else 0
    nrp = int(S.n_replicates.max())
    md = [
        "**eTable X.** Gain in held-out discrimination from adding nocturnal oxygen (T90) to "
        f"age and sex, one disease at a time. {len(order)} outcomes, the same set the paper "
        "ranks, in the frozen cohort of 19,173 patients. Baseline C is the held-out concordance "
        "of a Cox model with a restricted cubic spline in age and sex. C with T90 adds the "
        "within-site rank inverse normal of the percentage of the recording below 90% "
        "saturation. Models were fitted at one hospital and scored at the other, then the "
        "direction was reversed, and the two held-out values were averaged. The 2-year landmark "
        "columns repeat the whole procedure after removing every patient whose event or "
        "censoring arrived within 2 years of the sleep study and restarting the clock at 2 "
        "years, which is the rule the paper's published lag ladder uses. Confidence intervals "
        f"are the 2.5th and 97.5th percentiles of {nrp} bootstrap replicates resampling patients "
        "with replacement inside each hospital. Diseases are ordered by the size of the gain "
        "without the landmark. Only the last two rows average anything across diseases.",
        "",
        tab.to_markdown(index=False, disable_numparse=True),
        "",
        f"Gains whose interval excludes zero: {nex} of {len(order)} without the landmark, "
        f"{nexl} of {len(order)} at the 2-year landmark.",
        "",
        "---",
        "",
        "**eTable X, second panel.** The same 48 gains without the landmark, computed on the "
        "paper's published age basis with all four cr(df=4) columns kept. These point estimates "
        "are the ones already in numbers/ranking_v3_percondition.csv and are reproduced here to "
        "nine decimal places. The panel above drops one basis column, which is the rank-safe "
        "version of the same model and is the one this analysis reports.",
        "",
        tab2.to_markdown(index=False, disable_numparse=True),
        "",
    ]
    mdp = os.path.join(HERE, "dC_per_disease_TABLE.md")
    open(mdp, "w").write("\n".join(md))
    print(f"wrote {mdp}")

    summ = {"n_outcomes": len(order), "arms": arms,
            "n_replicates": {a: nrep.get((a, "drop"), 0) for a in arms}}
    for arm in arms:
        d2 = dis[dis.arm == arm]
        est = d2.estimable.astype(bool)          # object dtype, so ~ needs the cast
        summ[arm] = {
            "n_excluding_zero": int(d2.excludes_zero.astype(bool).sum()),
            "n_not_estimable": int((~est).sum()),
            "not_estimable": list(d2[~est].disease),
            "max_mc_error_of_ci_bound": float(d2.mc_error_of_ci_bound.max()),
            "boot_reps_losing_one_direction": int(d2.boot_reps_losing_one_direction.sum()),
            "boot_reps_with_no_estimate": int(d2.boot_reps_with_no_estimate.sum()),
            "outcomes_ever_losing_a_direction":
                int((d2.boot_reps_losing_one_direction > 0).sum()),
            "mean_all": df[(df.outcome == "__mean_all__") & (df.arm == arm)].iloc[0][
                ["dC", "ci_lo", "ci_hi", "boot_se", "n_outcomes_in_mean"]].to_dict(),
            "mean_common": df[(df.outcome == "__mean_common__") & (df.arm == arm)].iloc[0][
                ["dC", "ci_lo", "ci_hi", "boot_se", "n_outcomes_in_mean"]].to_dict(),
        }
    if "standard" in arms:
        summ["mean_dC_published_basis"] = float(
            df[(df.outcome == "__mean_all__") & (df.arm == "standard")].iloc[0].dC_pubbasis)
    json.dump(summ, open(os.path.join(E.WORK, "assemble_summary.json"), "w"), indent=2,
              default=float)

    out = dis[(dis.dC < dis.ci_lo) | (dis.dC > dis.ci_hi)]
    print(f"point estimate outside its own interval: {len(out)} of {len(dis)}"
          + ("" if out.empty else "  " + ", ".join(out.disease + " [" + out.arm + "]")))
    for arm in arms:
        m = df[(df.outcome == "__mean_all__") & (df.arm == arm)].iloc[0]
        print(f"{arm:<9} mean over {int(m.n_outcomes_in_mean)} outcomes "
              f"{m.dC:.6f} ({m.ci_lo:.6f} to {m.ci_hi:.6f})")
    print(f"\n{nex} of {len(order)} intervals exclude zero without the landmark, "
          f"{nexl} with it")


if __name__ == "__main__":
    main()
