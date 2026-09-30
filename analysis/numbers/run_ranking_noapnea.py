"""
The 197-measurement ranking, repeated in two restricted strata.

The reviewer objection this answers: oxygenation may top the parent ranking only because
84% of the parent cohort has sleep apnea. If oxygenation still ranks at the top among
people with no sleep apnea (AHI < 5), the objection does not hold.

Stratum 1: AHI < 5.
Stratum 2: AHI < 5 and free of prevalent organ disease (13 conditions).

Specification is copied from run_ranking_v2.py: age as a natural cubic spline with 4 df plus
sex as the adjustment set, each measurement added in turn as a rank-based inverse normal
score computed within hospital, Cox fitted in one hospital and scored in the other in both
directions, ranked by mean gain in the concordance statistic over the age-and-sex model.

Two departures from the parent, both deliberate and both reported:
  1. the impossible-oximetry block (spo2_nadir_corrected > spo2_mean) is dropped;
  2. the event floor is checked per direction and locked from the baseline model, so every
     measurement is scored on the same outcomes and the same directions.

Stages:
  A  positive control, reproduce frozen ranking_v2.csv rows on the parent cohort
  B  feasibility, incident events per outcome per hospital in each stratum
  C  the ranking itself, where feasible

Usage: python run_ranking_noapnea.py [A|B|C|all]
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../..")))  # repository root, where paths.py lives
import paths
import warnings; warnings.filterwarnings("ignore")
import json, sys, os, time
import numpy as np, pandas as pd
from scipy import stats
from patsy import dmatrix
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from joblib import Parallel, delayed

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from disease_definitions import DISEASES, CIRCULAR, RANKING_EXCLUDE
from cohort_spec import SENSITIVITY_ONLY, RANKING_DROP_SUFFIX, drop_split_nights_if_full   # v8.1   # v8: same candidate rule as build_ranking_v3

CIRC = set(CIRCULAR)
PARQUET = os.path.join(HERE, "..", "data_frozen_v8_2026-09", "t90_final.parquet")
MASTER = (f"{paths.TABLES_DIR}/"
          "master_cohort.csv")
FROZEN = os.path.join(HERE, "ranking_v3.csv")   # v8: the feature lock is the current ranking, not the frozen 197 file
OUTJSON = os.path.join(HERE, "ranking_noapnea.json")
SITES = ("I0002", "I0006")

# 13 conditions that define stratum 2's organ-disease exclusion
ORGAN_DZ = ["copd2", "asthma", "resp_failure", "obesity_hypovent", "pulm_htn", "hf", "ihd",
            "mi", "afib", "stroke_any", "ckd", "cirrhosis", "cancer_any"]

# floors. parent used 30 events in each hospital, both directions required.
MIN_TRAIN, MIN_TEST = 40, 25          # primary, as commissioned
PARENT_FLOOR = 30                     # sensitivity, parent's own rule
# hard degeneracy guard applied after a feature's own missingness is dropped
GUARD_TRAIN, GUARD_TEST = 25, 15

SHORTLIST = ["spo2_pct_below_90", "spo2_nadir_corrected", "odi3_total", "AHI", "TST_min",
             "arousal_index", "sleep_efficiency_pct", "N3_pct", "nrem_hr_bpm"]
REPORT = ["spo2_pct_below_90", "spo2_nadir_corrected", "AHI", "TST_min", "arousal_index",
          "sleep_efficiency_pct", "N3_pct", "odi3_total", "nrem_hr_bpm"]

RESULTS = {}


# ---------------------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------------------
def load_base(drop_artifact):
    """Parent cohort filter, optionally minus the impossible-oximetry block."""
    d = pd.read_parquet(PARQUET)
    b = d[(d.fu_valid == 1) & d.spo2_pct_below_90.notna()].copy()
    b = drop_split_nights_if_full(b)   # v8.1: the ranking family runs on the full diagnostic nights
    n_filtered = len(b)
    bad = (b.spo2_nadir_corrected > b.spo2_mean)
    n_bad = int(bad.sum())
    bad_by_site = b.loc[bad, "site_id"].value_counts().to_dict()
    if drop_artifact:
        b = b[~bad].copy()
    return b, dict(n_after_cohort_filter=n_filtered, n_impossible_oximetry=n_bad,
                   impossible_by_site={k: int(v) for k, v in bad_by_site.items()},
                   n_analysed=len(b))


def attach_features(base):
    """Merge the master feature table exactly as the parent script does."""
    head = pd.read_csv(MASTER, nrows=0).columns.tolist()
    FEAT = [c for c in head if c not in ("BDSPPatientID",)
            and not c.endswith(("_first_date", "_date"))]
    raw = pd.read_csv(MASTER, usecols=["BDSPPatientID"] + FEAT, low_memory=False)
    keep = ["BDSPPatientID", "site_id", "AgeAtVisit", "sex", "AHI_cohort"] + \
           [c for c in base.columns if c.endswith(("_incident", "_years", "_prevalent"))]
    base = base.copy()
    base["AHI_cohort"] = base["AHI"]          # stratifier from the frozen cohort file
    left = base[keep]
    raw = raw.drop(columns=[c for c in raw.columns
                            if c in set(left.columns) - {"BDSPPatientID"}])
    FEAT = [c for c in FEAT if c in raw.columns]
    d = left.merge(raw, on="BDSPPatientID", how="left")
    d["male"] = (d.sex.astype(str).str.upper().str[0] == "M").astype(int)
    return d, FEAT


def parent_feature_list(d, FEAT):
    """The parent's own screen: numeric, >=60% present, >10 distinct values."""
    num = [c for c in FEAT if pd.api.types.is_numeric_dtype(d[c])]
    num = [c for c in num if d[c].notna().mean() >= 0.60 and d[c].nunique() > 10]
    DROP = {"male", "AgeAtVisit", "AHI_cohort"}
    return [c for c in num if c not in DROP and not c.startswith(("obs_", "follow_"))
            and c not in SENSITIVITY_ONLY and not c.endswith(RANKING_DROP_SUFFIX)]   # v8: the parent's exclusions


def add_design(d, feats):
    """Age spline, then the within-hospital rank-based inverse normal score."""
    d = d.reset_index(drop=True)
    sp = dmatrix("cr(a, df=4) - 1", {"a": d.AgeAtVisit.values}, return_type="dataframe")
    sp.columns = [f"age_s{i}" for i in range(4)]
    sp.index = d.index
    for c in sp.columns:
        d[c] = sp[c]
    for f in feats:
        z = pd.Series(np.nan, index=d.index)
        for s, idx in d.groupby("site_id").groups.items():
            v = d.loc[idx, f]
            ok = v.notna()
            if ok.sum() < 20:
                continue
            r = stats.rankdata(v[ok], method="average")
            z.loc[v[ok].index] = stats.norm.ppf((r - 0.375) / (ok.sum() + 0.25))
        d[f + "__z"] = z
    return d


ADJ = [f"age_s{i}" for i in range(4)] + ["male"]


def outcome_list(d):
    OUT = [k for k in DISEASES if k not in CIRC and k not in RANKING_EXCLUDE
           and f"{k}_incident" in d.columns]
    OUT = [k for k in OUT if int(d[f"{k}_incident"].sum()) >= 150] + ["death"]
    return OUT


def at_risk(g0, o):
    return g0[(g0[f"{o}_prevalent"] == 0) & (g0[f"{o}_years"] > 0) & g0[f"{o}_years"].notna()]


def event_table(g0, OUT):
    rows = []
    for o in OUT:
        g = at_risk(g0, o)
        row = {"outcome": o}
        for s in SITES:
            gs = g[g.site_id == s]
            row[f"n_{s}"] = int(len(gs))
            row[f"ev_{s}"] = int(gs[f"{o}_incident"].sum())
        rows.append(row)
    return pd.DataFrame(rows)


def eligible_directions(ev_row, min_train, min_test):
    out = []
    for tr, te in ((SITES[0], SITES[1]), (SITES[1], SITES[0])):
        if ev_row[f"ev_{tr}"] >= min_train and ev_row[f"ev_{te}"] >= min_test:
            out.append((tr, te))
    return out


# ---------------------------------------------------------------------------------------
# the fit
# ---------------------------------------------------------------------------------------
def one(pkl, feat, out, dirs, guard):
    dd = pd.read_pickle(pkl)
    g = at_risk(dd, out)
    cols = ADJ + ([feat + "__z"] if feat else [])
    res = []
    for tr, te in dirs:
        A = g[g.site_id == tr][cols + [f"{out}_years", f"{out}_incident"]].dropna()
        B = g[g.site_id == te][cols + [f"{out}_years", f"{out}_incident"]].dropna()
        if guard and (A[f"{out}_incident"].sum() < GUARD_TRAIN
                      or B[f"{out}_incident"].sum() < GUARD_TEST):
            res.append(np.nan); continue
        if len(B) < 20:
            res.append(np.nan); continue
        try:
            c = CoxPHFitter(penalizer=0.01).fit(A, duration_col=f"{out}_years",
                                                event_col=f"{out}_incident")
            r = c.predict_partial_hazard(B)
            res.append(concordance_index(B[f"{out}_years"], -r, B[f"{out}_incident"]))
        except Exception:
            res.append(np.nan)
    if not res or np.all(np.isnan(res)):
        return {"feature": feat or "__BASE__", "outcome": out, "c": np.nan}
    return {"feature": feat or "__BASE__", "outcome": out, "c": float(np.nanmean(res))}


def run_ranking(d, feats, OUT, dirs_by_out, tag, guard=True, n_jobs=7):
    pkl = f"/tmp/rank_{tag}.pkl"
    d.to_pickle(pkl)
    jobs = [(None, o) for o in OUT] + [(f, o) for f in feats for o in OUT]
    print(f"[{tag}] {len(feats)} measures x {len(OUT)} outcomes -> {len(jobs)} model pairs")
    t0 = time.time()
    rr = Parallel(n_jobs=n_jobs, verbose=1, backend="loky")(
        delayed(one)(pkl, f, o, dirs_by_out[o], guard) for f, o in jobs)
    print(f"[{tag}] {time.time()-t0:.0f}s")
    s = pd.DataFrame(rr)
    bas = s[s.feature == "__BASE__"].set_index("outcome")["c"]
    f2 = s[s.feature != "__BASE__"].copy()
    f2["gain"] = f2.apply(lambda r: r.c - bas.get(r.outcome, np.nan), axis=1)
    f2 = f2[f2.gain.notna()]
    rank = (f2.groupby("feature")
              .agg(mC=("c", "mean"), dC=("gain", "mean"), worst=("gain", "min"),
                   npos=("gain", lambda x: int((x > 0).sum())), nout=("gain", "size"))
              .reset_index().sort_values("dC", ascending=False))
    rank.insert(0, "rank", range(1, len(rank) + 1))
    return rank, f2, bas


def positions(rank, keys):
    out = {}
    for k in keys:
        r = rank[rank.feature == k]
        out[k] = (None if not len(r) else
                  {"rank": int(r["rank"].iloc[0]), "of": int(len(rank)),
                   "dC": float(r.dC.iloc[0]), "mC": float(r.mC.iloc[0]),
                   "npos": int(r.npos.iloc[0]), "nout": int(r.nout.iloc[0])})
    return out


# =======================================================================================
def stage_A():
    """Positive control against the frozen full-cohort result."""
    print("\n=== STAGE A  positive control ===")
    base, prov = load_base(drop_artifact=False)
    print(f"parent cohort after filter: {len(base)}   (expected 19383)")
    d, FEAT = attach_features(base)
    num = parent_feature_list(d, FEAT)
    print(f"candidate measures: {len(num)}   (feature lock: {os.path.basename(FROZEN)})")
    frozen = pd.read_csv(FROZEN, comment="#")
    assert set(num) == set(frozen.feature), \
        f"feature list drift: {set(num) ^ set(frozen.feature)}"
    check = ["spo2_nadir_corrected", "spo2_pct_below_90", "odi3_total", "nrem_hr_bpm",
             "AHI", "arousal_index", "TST_min"]
    d = add_design(d, check)
    OUT = outcome_list(d)
    print(f"outcomes: {len(OUT)}   (frozen file has {int(frozen.nout.max())})")
    ev = event_table(d, OUT)
    # parent rule exactly: both directions, >=30 events each side, checked after dropna
    dirs = {o: [(SITES[0], SITES[1]), (SITES[1], SITES[0])] for o in OUT}
    pkl = "/tmp/rank_pc.pkl"; d.to_pickle(pkl)
    jobs = [(None, o) for o in OUT] + [(f, o) for f in check for o in OUT]
    rr = Parallel(n_jobs=7, verbose=1, backend="loky")(
        delayed(one_parent)(pkl, f, o) for f, o in jobs)
    s = pd.DataFrame(rr)
    bas = s[s.feature == "__BASE__"].set_index("outcome")["c"]
    f2 = s[s.feature != "__BASE__"].copy()
    f2["gain"] = f2.apply(lambda r: r.c - bas.get(r.outcome, np.nan), axis=1)
    got = f2.groupby("feature").agg(mC=("c", "mean"), dC=("gain", "mean")).reset_index()
    rows = []
    for _, r in got.iterrows():
        fr = frozen[frozen.feature == r.feature].iloc[0]
        rows.append(dict(feature=r.feature, frozen_rank=int(fr["rank"]),
                         frozen_dC=float(fr.dC), repro_dC=float(r.dC),
                         d_dC=float(r.dC - fr.dC),
                         frozen_mC=float(fr.mC), repro_mC=float(r.mC),
                         d_mC=float(r.mC - fr.mC)))
    pc = pd.DataFrame(rows).sort_values("frozen_rank")
    print(pc.to_string(index=False))
    worst = float(np.abs(pc.d_dC).max())
    print(f"\nlargest absolute difference in mean C gain: {worst:.2e}")
    print(f"baseline held-out C {bas.mean():.4f}  (parent printed 0.6474 territory)")
    RESULTS["stage_A_positive_control"] = {
        "n_cohort": len(base), "n_measures": len(num), "n_outcomes": len(OUT),
        "baseline_mean_C": float(bas.mean()),
        "max_abs_diff_dC": worst,
        "matched": bool(worst < 1e-9),
        "rows": pc.to_dict("records"),
    }
    RESULTS["provenance_parent"] = prov
    return worst


def one_parent(pkl, feat, out):
    """Byte-faithful copy of the parent's inner function, for the positive control."""
    dd = pd.read_pickle(pkl)
    g = dd[(dd[f"{out}_prevalent"] == 0) & (dd[f"{out}_years"] > 0) & dd[f"{out}_years"].notna()]
    cols = ADJ + ([feat + "__z"] if feat else [])
    res = []
    for tr, te in (("I0002", "I0006"), ("I0006", "I0002")):
        A = g[g.site_id == tr][cols + [f"{out}_years", f"{out}_incident"]].dropna()
        B = g[g.site_id == te][cols + [f"{out}_years", f"{out}_incident"]].dropna()
        if A[f"{out}_incident"].sum() < 30 or B[f"{out}_incident"].sum() < 30:
            res.append(np.nan); continue
        try:
            c = CoxPHFitter(penalizer=0.01).fit(A, duration_col=f"{out}_years",
                                                event_col=f"{out}_incident")
            r = c.predict_partial_hazard(B)
            res.append(concordance_index(B[f"{out}_years"], -r, B[f"{out}_incident"]))
        except Exception:
            res.append(np.nan)
    return {"feature": feat or "__BASE__", "outcome": out, "c": float(np.nanmean(res))}


def build_strata():
    base, prov = load_base(drop_artifact=True)
    d, FEAT = attach_features(base)
    num = parent_feature_list(d, FEAT)        # screened on the analysed cohort
    frozen = pd.read_csv(FROZEN, comment="#")
    feats = [f for f in frozen.feature if f in d.columns]   # lock to the ranking's feature list
    lost = [f for f in frozen.feature if f not in d.columns]
    s1 = d[d.AHI_cohort < 5].copy()
    prev = np.zeros(len(s1), dtype=bool)
    for k in ORGAN_DZ:
        prev |= (s1[f"{k}_prevalent"] == 1).values
    s2 = s1[~prev].copy()
    prov.update(dict(
        n_stratum1=len(s1), n_stratum2=len(s2),
        stratum1_by_site={s: int((s1.site_id == s).sum()) for s in SITES},
        stratum2_by_site={s: int((s2.site_id == s).sum()) for s in SITES},
        n_excluded_organ_disease=int(prev.sum()),
        n_measures_locked=len(feats), measures_lost=lost,
        parent_screen_on_analysed_cohort=len(num)))
    return d, s1, s2, feats, prov


def stage_B(d, s1, s2, feats):
    print("\n=== STAGE B  feasibility ===")
    OUT = outcome_list(d)
    print(f"outcomes carried from the parent run: {len(OUT)}")
    tabs, summ = {}, {}
    for name, g0 in (("full_minus_artifact", d), ("stratum1_no_apnea", s1),
                     ("stratum2_no_apnea_no_organ_dz", s2)):
        ev = event_table(g0, OUT)
        ev.to_csv(os.path.join(HERE, f"events_{name}.csv"), index=False)
        n_primary = int(sum(len(eligible_directions(r, MIN_TRAIN, MIN_TEST)) > 0
                            for _, r in ev.iterrows()))
        n_both = int(sum(len(eligible_directions(r, MIN_TRAIN, MIN_TEST)) == 2
                         for _, r in ev.iterrows()))
        n_parent = int(sum(len(eligible_directions(r, PARENT_FLOOR, PARENT_FLOOR)) == 2
                           for _, r in ev.iterrows()))
        print(f"{name:<32} n={len(g0):>6}  "
              f"floor40/25 any-direction {n_primary:>2}/{len(OUT)}  "
              f"both {n_both:>2}  parent30/30 both {n_parent:>2}")
        tabs[name] = ev
        summ[name] = dict(n=len(g0),
                          n_by_site={s: int((g0.site_id == s).sum()) for s in SITES},
                          outcomes_pass_40_25_any_direction=n_primary,
                          outcomes_pass_40_25_both_directions=n_both,
                          outcomes_pass_parent_30_30_both=n_parent,
                          events=ev.to_dict("records"))
    RESULTS["stage_B_feasibility"] = summ
    return OUT, tabs


def stage_C(s1, s2, feats, OUT, tabs):
    print("\n=== STAGE C  ranking ===")
    out = {}
    for name, g0, do_full in (("stratum1_no_apnea", s1, True),
                              ("stratum2_no_apnea_no_organ_dz", s2, False)):
        ev = tabs[name].set_index("outcome")
        dirs = {o: eligible_directions(ev.loc[o], MIN_TRAIN, MIN_TEST) for o in OUT}
        keep = [o for o in OUT if dirs[o]]
        print(f"\n--- {name}: {len(keep)} outcomes clear the floor")
        if len(keep) < 10:
            print("    fewer than 10 outcomes. the full ranking would be noise.")
            print("    running the pre-chosen shortlist only, flagged exploratory.")
            use, mode = [f for f in SHORTLIST if f in feats], "shortlist_exploratory"
        elif do_full:
            use, mode = feats, "full_197"
        else:
            use, mode = [f for f in SHORTLIST if f in feats], "shortlist"
        g = add_design(g0, use)          # rank transform recomputed within this stratum
        rank, per, bas = run_ranking(g, use, keep, dirs, name.split("_")[0])
        rank.to_csv(os.path.join(HERE, f"ranking_noapnea_{name}.csv"), index=False)
        per.to_csv(os.path.join(HERE, f"ranking_noapnea_{name}_percondition.csv"), index=False)
        print(f"baseline held-out C {bas.mean():.4f}")
        print(rank.head(20).to_string(index=False))
        pos = positions(rank, REPORT)
        for k, v in pos.items():
            if v:
                print(f"   {k:<24}rank {v['rank']:>3} of {v['of']}   gain {v['dC']:+.5f}"
                      f"   positive in {v['npos']}/{v['nout']}")
        out[name] = dict(mode=mode, n=len(g0), n_outcomes=len(keep), outcomes=keep,
                         baseline_mean_C=float(bas.mean()),
                         n_measures_ranked=int(len(rank)),
                         top20=rank.head(20).to_dict("records"),
                         key_positions=pos)

        # sensitivity: parent's own 30/30 both-direction rule
        d2 = {o: eligible_directions(ev.loc[o], PARENT_FLOOR, PARENT_FLOOR) for o in OUT}
        k2 = [o for o in OUT if len(d2[o]) == 2]
        if len(k2) >= 8 and set(k2) != set(keep):
            r2, _, b2 = run_ranking(g, use, k2, d2, name.split("_")[0] + "_sens")
            out[name]["sensitivity_parent_floor"] = dict(
                n_outcomes=len(k2), baseline_mean_C=float(b2.mean()),
                top10=[r.feature for _, r in r2.head(10).iterrows()],
                key_positions=positions(r2, REPORT))
            print(f"  sensitivity (parent 30/30, {len(k2)} outcomes) top 5: "
                  f"{[r.feature for _,r in r2.head(5).iterrows()]}")
    RESULTS["stage_C_ranking"] = out


def stage_D():
    """Robustness. Three questions the stratum-1 result has to survive."""
    print("\n=== STAGE D  robustness ===")
    out = {}
    b = pd.read_parquet(PARQUET)
    b = b[(b.fu_valid == 1) & b.spo2_pct_below_90.notna() & (b.oximetry_bad == 0)].copy()

    # D1. does dropping the impossible-oximetry block change the stratum-1 order?
    d1 = {}
    for tag, keep_art in (("with_artifact_rows", True), ("excl_artifact", False)):
        base = b if keep_art else b[~(b.spo2_nadir_corrected > b.spo2_mean)]
        dd, _ = attach_features(base.copy())
        s1 = dd[dd.AHI_cohort < 5].copy()
        ev = event_table(s1, outcome_list(dd)).set_index("outcome")
        dirs = {o: eligible_directions(ev.loc[o], MIN_TRAIN, MIN_TEST) for o in ev.index}
        keep = [o for o in ev.index if dirs[o]]
        g = add_design(s1, SHORTLIST)
        rank, _, _ = run_ranking(g, SHORTLIST, keep, dirs, "art_" + tag)
        d1[tag] = dict(n=len(s1), n_outcomes=len(keep), order=list(rank.feature),
                       dC={r.feature: float(r.dC) for _, r in rank.iterrows()})
        print(f"[{tag}] n={len(s1)} order={list(rank.feature)}")
    d1["order_identical"] = d1["with_artifact_rows"]["order"] == d1["excl_artifact"]["order"]
    out["D1_artifact_exclusion"] = d1

    # D2. how much of the stratum-1 gain rides on the 5 single-direction outcomes?
    per = pd.read_csv(os.path.join(HERE, "ranking_noapnea_stratum1_no_apnea_percondition.csv"))
    ev = pd.read_csv(os.path.join(HERE, "events_stratum1_no_apnea.csv")).set_index("outcome")
    both = [o for o in ev.index if len(eligible_directions(ev.loc[o], PARENT_FLOOR, PARENT_FLOOR)) == 2]
    single = sorted(set(per.outcome) - set(both))
    d2 = {"both_direction_outcomes": both, "single_direction_outcomes": single}
    for f in REPORT:
        x = per[per.feature == f].set_index("outcome")["gain"]
        d2[f] = dict(all19=float(x.mean()),
                     both14=float(x[[o for o in x.index if o in both]].mean()),
                     single5=float(x[[o for o in x.index if o in single]].mean()))
        print(f"  {f:<22} all19 {d2[f]['all19']:+.5f}  both14 {d2[f]['both14']:+.5f}  "
              f"single5 {d2[f]['single5']:+.5f}")
    out["D2_direction_dependence"] = d2

    # D3. is any measure distinguishable from any other? gain +/- SE across outcomes
    rows = []
    for f in per.feature.unique():
        x = per[(per.feature == f) & per.outcome.isin(both)]["gain"]
        se = x.std(ddof=1) / np.sqrt(len(x))
        rows.append(dict(feature=f, mean=float(x.mean()), se=float(se),
                         lo=float(x.mean() - 1.96 * se), hi=float(x.mean() + 1.96 * se)))
    st = pd.DataFrame(rows).sort_values("mean", ascending=False)
    st.to_csv(os.path.join(HERE, "ranking_noapnea_stratum1_uncertainty.csv"), index=False)
    out["D3_uncertainty"] = dict(
        n_outcomes=len(both),
        n_measures_ci_excludes_zero=int((st.lo > 0).sum()),
        measures_ci_excludes_zero=list(st[st.lo > 0].feature),
        key=st[st.feature.isin(REPORT)].to_dict("records"))
    print(f"  measures whose 95% CI excludes zero: {int((st.lo>0).sum())} of {len(st)} "
          f"-> {list(st[st.lo>0].feature)}")

    # D4. missingness of the oximetry measures inside stratum 1
    c = b[~(b.spo2_nadir_corrected > b.spo2_mean)]
    s1 = c[c.AHI < 5]
    out["D4_missingness_stratum1"] = {
        f: dict(n_missing=int(s1[f].isna().sum()), n=len(s1),
                by_site={k: int(v) for k, v in
                         s1.loc[s1[f].isna(), "site_id"].value_counts().items()})
        for f in ("spo2_nadir_corrected", "spo2_mean", "spo2_pct_below_90", "AHI")}
    RESULTS["stage_D_robustness"] = out


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what == "D":
        stage_D()
    if what in ("A", "all"):
        stage_A()
    if what in ("B", "C", "all"):
        d, s1, s2, feats, prov = build_strata()
        RESULTS["provenance_strata"] = prov
        print(json.dumps(prov, indent=2, default=str))
        OUT, tabs = stage_B(d, s1, s2, feats)
        if what in ("C", "all"):
            stage_C(s1, s2, feats, OUT, tabs)
    if os.path.exists(OUTJSON):
        old = json.load(open(OUTJSON))
        old.update(RESULTS); RESULTS = old
    json.dump(RESULTS, open(OUTJSON, "w"), indent=2, default=str)
    print(f"\nwrote {OUTJSON}")
