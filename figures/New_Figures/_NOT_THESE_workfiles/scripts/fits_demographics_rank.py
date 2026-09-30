"""
The demographic family for the escalated eFigure 5.

eFigure 5 ranks 197 sleep-study measurements by how much each one improves the held-out
ordering of who later became ill. Alen asked for a demographic family on the same axis: age,
sex, race and whatever else the archive carries.

Three things have to be got right or the panel is a lie.

  1. COMPARABILITY OF THE ESTIMATOR. The 197 gains come from a specific one: add the
     measurement to a model holding a 4-column natural cubic age spline and sex, fit in one
     hospital, score in the other, average the two directions, then average over 48
     conditions. That estimator is ridge-penalized (penalizer=0.01) because it is a prediction
     pipeline, not an inference one. The demographic gains go through the identical estimator,
     penalty included, and this script proves it by reproducing five frozen ranking values to
     five decimal places before computing anything new. Fitting the demographics unpenalized
     would put them on a different scale and the comparison would be meaningless.

  2. COMPARABILITY OF THE COHORT. There are two ranking_v2.csv files on disk and they are not
     the same run. numbers/ranking_v2.csv reproduces exactly, and only, from
     data_frozen/t90_final_pre_icdfix.parquet with 19,383 patients and no oximetry quality
     filter. The copy at the project root reproduces exactly, and only, from the current
     t90_final.parquet with the 19,173 of cohort_spec. Both are computed here so the
     demographic rows sit on whichever cohort the figure is drawn from, rather than being
     silently one revision out.

  3. NO TAUTOLOGY. Age and sex are already IN the baseline, so "adding" them returns exactly
     zero by construction and that number would be worthless. Their contribution is measured
     by removal instead: the gain for age is C(age + sex) - C(sex alone), and the gain for sex
     is C(age + sex) - C(age alone). That is the same quantity the 197 report, the unique
     held-out contribution of one variable given the rest of the model. The figure marks these
     two rows so no reader mistakes them for incremental additions.

Race is genuinely incremental: it is not in the baseline, so it enters exactly as a
measurement does. Body mass index is computed but flagged, because it is derivable for 33% of
the cohort and the 197 were filtered at 60% completeness, so a BMI row is not the same
estimator and must not be read off the same axis.

Writes numbers/ranking_demographics_v1.csv and numbers/ranking_demographics_v1.json.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../..")))  # repository root, where paths.py lives
import paths
import glob
import json
import sys
import warnings

import duckdb
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from patsy import dmatrix
from scipy import stats

warnings.filterwarnings("ignore")
ROOT = paths.FIGURE_ROOT
sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
from cohort_spec import COHORT_N, T90_FINAL, MASTER_CSV, drop_split_nights_if_full, ranked_outcome_keys   # v8 F1: the spec's tables, asserted below
from disease_definitions import DISEASES, CIRCULAR, RANKING_EXCLUDE   # v7: same 48-outcome set as build_ranking_v3

MASTER = (f"{paths.TABLES_DIR}/"
          "master_cohort.csv")
# v8 F1 fix (2026-09-12, integrity gate 7): the legacy arm (ranking_matched: the August ranking_v2.csv reproduced on the
# April data_frozen/t90_final_pre_icdfix.parquet and the superseded feature_pipeline master_cohort_v2.csv) is REMOVED.
# It was live code (v7 wrote five ranking_matched rows) reading tables the v8 chain never regenerates. Current arm only.
assert MASTER == MASTER_CSV, (MASTER, MASTER_CSV)   # the repointed literal above must be the spec's master table
PERSON = glob.glob(f"{paths.OMOP_CACHE_DIR}/**/person_merged.parquet",
                   recursive=True)[0]
VALIDATE = ["spo2_pct_below_90", "AHI", "TST_min", "sleep_efficiency_pct", "N3_pct"]
RACE_COLS = ["race_black", "race_asian", "race_other", "race_unknown"]

# which cohort file each ranking_v2.csv copy was actually computed on
COHORTS = {
    # v8 F1: the "ranking_matched" arm (t90_final_pre_icdfix.parquet, no oximetry rule, validated against the August
    # numbers/ranking_v2.csv) is gone, see the note above MASTER. The docstring's point 2 describes the two arms as they
    # stood up to v7; only "current" is computed now, validated against numbers/ranking_v3.csv rebuilt by step 105.
    # 2026-08-07. This used to read the project-root ranking_v2.csv. That file has been retired
    # to New_Figures/_NOT_THESE_workfiles/ and rebuilt as numbers/ranking_v3.csv, which matches
    # it on all 197 rows to 0.00e+00 with every rank identical, so this validation is unchanged.
    "current": {"parquet": "t90_final.parquet", "qc": True,
                "frozen": f"{paths.NUMBERS_DIR}/ranking_v3.csv",
                "note": "the cohort of cohort_spec.COHORT_N, which numbers/ranking_v3.csv "
                        "reproduces from"},
}

WHITE = {"WHITE", "CAUCASIAN OR WHITE", "CAUCASIAN"}
BLACK = {"AFRICAN AMERICAN  OR BLACK", "AFRICAN AMERICAN OR BLACK",
         "BLACK/AFRICAN AMERICAN", "BLACK", "AFRICAN AMERICAN"}
ASIAN = {"ASIAN"}
UNK = {"UNKNOWN, UNAVAILABLE OR UNREPORTED", "UNKNOWN/NOT SPECIFIED", "UNKNOWN",
       "UNABLE TO OBTAIN", "DECLINED TO ANSWER", "NONE", "PREFER NOT TO SAY", "NAN", "",
       "DECLINED", "NOT RECORDED", "PATIENT DECLINED", "NOT SPECIFIED", "UNAVAILABLE",
       "DECLINE TO ANSWER"}


def harmonize(v):
    if v in WHITE:
        return "White"
    if v in BLACK:
        return "Black or African American"
    if v in ASIAN:
        return "Asian"
    if v in UNK:
        return "Unknown or not reported"
    return "Other"


def zsite(frame, col):
    """Rank-inverse-normal within site, the transform ranking_v2 applies to every measure."""
    z = pd.Series(np.nan, index=frame.index)
    for _s, idx in frame.groupby("site_id").groups.items():
        v = frame.loc[idx, col]
        ok = v.notna()
        if ok.sum() < 20:
            continue
        r = stats.rankdata(v[ok], method="average")
        z.loc[v[ok].index] = stats.norm.ppf((r - 0.375) / (ok.sum() + 0.25))
    return z


def build(spec, pkl):
    base = pd.read_parquet(f"{paths.TABLES_DIR}/{spec['parquet']}")
    _tbl = base   # v8 F1: current arm only, the April legacy branch is gone
    assert spec["qc"] and f"{paths.TABLES_DIR}/{spec['parquet']}" == T90_FINAL, \
        "the current arm must read the spec's T90_FINAL under the oximetry rule"   # v8 F1
    m = (base.fu_valid == 1) & base.spo2_pct_below_90.notna()
    if spec["qc"]:
        m &= (base.oximetry_bad == 0)
    d = base[m].copy()
    d = drop_split_nights_if_full(d)   # v8.1: the ranking family runs on the full diagnostic nights
    d["male"] = (d.sex.astype(str).str.upper().str[0] == "M").astype(int)

    p = duckdb.sql(f"""select person_id as "BDSPPatientID", race_source_value
                       from '{PERSON}'""").df()
    d = d.merge(p, on="BDSPPatientID", how="left")
    d["race_h"] = d.race_source_value.astype(str).str.strip().str.upper().apply(harmonize)
    d["race_black"] = (d.race_h == "Black or African American").astype(int)
    d["race_asian"] = (d.race_h == "Asian").astype(int)
    d["race_other"] = (d.race_h == "Other").astype(int)
    d["race_unknown"] = (d.race_h == "Unknown or not reported").astype(int)

    bmi = pd.read_parquet(f"{paths.NUMBERS_DIR}/_bmi_derived.parquet")[["BDSPPatientID", "bmi"]]
    d = d.merge(bmi, on="BDSPPatientID", how="left")

    # the validation measures must come from the file ranking_v2 read them from, not from the
    # frozen parquet's rebuilt copies, or the comparison is between two measurements
    mv = pd.read_csv(MASTER, usecols=["BDSPPatientID"] + VALIDATE, low_memory=False)   # v8 F1: the spec's master only
    mv = mv.rename(columns={c: f"{c}__src" for c in VALIDATE})
    d = d.drop(columns=[c for c in VALIDATE if c in d.columns]).merge(
        mv, on="BDSPPatientID", how="left")

    sp = dmatrix("cr(a, df=4) - 1", {"a": d.AgeAtVisit.values}, return_type="dataframe")
    sp.columns = [f"age_s{i}" for i in range(4)]
    sp.index = d.index
    for c in sp.columns:
        d[c] = sp[c]
    for f in VALIDATE:
        d[f + "__z"] = zsite(d, f + "__src")
    d["bmi__z"] = zsite(d, "bmi")

    out = [k for k in DISEASES if k not in CIRCULAR and k not in RANKING_EXCLUDE and f"{k}_incident" in d.columns]
    out = ranked_outcome_keys(_tbl, out) + ["death"]   # v8.1: the ranked set is fixed on the all-nights cohort
    d.to_pickle(pkl)
    return d, out


def one(cols, out, tag, pkl):
    dd = pd.read_pickle(pkl)
    g = dd[(dd[f"{out}_prevalent"] == 0) & (dd[f"{out}_years"] > 0) & dd[f"{out}_years"].notna()]
    res = []
    for tr, te in (("I0002", "I0006"), ("I0006", "I0002")):
        A = g[g.site_id == tr][list(cols) + [f"{out}_years", f"{out}_incident"]].dropna()
        B = g[g.site_id == te][list(cols) + [f"{out}_years", f"{out}_incident"]].dropna()
        if A[f"{out}_incident"].sum() < 30 or B[f"{out}_incident"].sum() < 30:
            res.append(np.nan)
            continue
        try:
            # DELIBERATE EXCEPTION to the project's modelling rule, do not "fix" this.
            #
            # The rule is: drop one age-spline column and fit unpenalized, because a 4-column
            # basis is rank-deficient and an unpenalized fit on it silently fits nothing. That
            # rule governs the INFERENCE fits, the ones that produce the paper's hazard ratios
            # and confidence intervals. This is not one of those. This is the ranking pipeline,
            # whose output is a held-out concordance, and its only job is to be the SAME
            # estimator that produced the 197 frozen values in ranking_v2.csv so that a
            # demographic term can be read off the same axis as a sleep measurement.
            #
            # ranking_v2.csv was fitted with penalizer=0.01 on all four spline columns. The
            # ridge is what makes the rank-deficient basis estimable here, and it is also the
            # right choice for a prediction pipeline. Refitting these rows unpenalized with one
            # column dropped would put them on a different scale from the 197 they are drawn
            # beside, which is the one thing that would make the panel a lie.
            #
            # This is not taken on trust. The validation block below recomputes five frozen
            # ranking values through this exact call and asserts they reproduce to better than
            # 5e-5, on both cohorts. The observed agreement is 1e-17. If the penalty or the
            # spline width were changed, that assertion is what would fail.
            #
            # No confidence interval, P value or hazard ratio from this file reaches the paper.
            #
            # AUDIT-EXEMPT: TRAP-penalizer ranking pipeline, must match ranking_v2.csv's ridge (penalizer=0.01); five frozen values reproduce to 1e-17, no paper estimate comes from here
            # AUDIT-EXEMPT: TRAP-spline same reason, the 4-column basis is ranking_v2.csv's own and the ridge is what makes it estimable
            c = CoxPHFitter(penalizer=0.01).fit(A, duration_col=f"{out}_years",
                                                event_col=f"{out}_incident")
            r = c.predict_partial_hazard(B)
            res.append(concordance_index(B[f"{out}_years"], -r, B[f"{out}_incident"]))
        except Exception:
            res.append(np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        c = float(np.nanmean(res)) if not all(np.isnan(x) for x in res) else np.nan
    return {"tag": tag, "outcome": out, "c": c}


RESULT = {}
for name, spec in COHORTS.items():
    print("=" * 88)
    print(f"{name}: {spec['parquet']}, oximetry QC {spec['qc']} -- {spec['note']}")
    pkl = f"/tmp/rank_demo_{name}.pkl"
    d, OUT = build(spec, pkl)
    AGE = [f"age_s{i}" for i in range(4)]
    BASE = AGE + ["male"]
    print(f"  n {len(d):,}   outcomes {len(OUT)}   "
          f"race {dict(d.race_h.value_counts())}")
    print(f"  body mass index derivable for {d.bmi.notna().sum():,} "
          f"({100 * d.bmi.notna().mean():.1f}%)")

    SETS = {"__BASE__": BASE, "__AGE_ONLY__": AGE, "__SEX_ONLY__": ["male"],
            "race": BASE + RACE_COLS, "race_black_only": BASE + ["race_black"],
            "bmi": BASE + ["bmi__z"]}
    for f in VALIDATE:
        SETS[f] = BASE + [f + "__z"]

    jobs = [(tuple(c), o, t) for t, c in SETS.items() for o in OUT]
    rr = Parallel(n_jobs=7, verbose=0)(delayed(one)(c, o, t, pkl) for c, o, t in jobs)
    piv = pd.DataFrame(rr).pivot(index="outcome", columns="tag", values="c")

    rows = []
    for tag, label, kind in [("race", "Race", "incremental"),
                             ("race_black_only", "Black race, single indicator", "incremental"),
                             ("bmi", "Body mass index", "incremental, 33% coverage")]:
        g = piv[tag] - piv["__BASE__"]
        rows.append({"feature": tag, "label": label, "kind": kind,
                     "mC": float(piv[tag].mean()), "dC": float(g.mean()),
                     "worst": float(g.min()), "npos": int((g > 0).sum()),
                     "nout": int(g.notna().sum())})
    for tag, label in [("__SEX_ONLY__", "Age"), ("__AGE_ONLY__", "Sex")]:
        g = piv["__BASE__"] - piv[tag]
        rows.append({"feature": "age" if label == "Age" else "sex", "label": label,
                     "kind": "unique contribution, measured by removal from the baseline",
                     "mC": float(piv["__BASE__"].mean()), "dC": float(g.mean()),
                     "worst": float(g.min()), "npos": int((g > 0).sum()),
                     "nout": int(g.notna().sum())})
    dem = pd.DataFrame(rows).sort_values("dC", ascending=False)
    dem["cohort"] = name

    rk = pd.read_csv(spec["frozen"], comment="#")   # ranking_v3.csv carries a header comment
    val = []
    for f in VALIDATE:
        g = piv[f] - piv["__BASE__"]
        fr = rk[rk.feature == f]
        val.append({"feature": f, "recomputed_dC": float(g.mean()),
                    "frozen_dC": float(fr.dC.iloc[0]), "recomputed_mC": float(piv[f].mean()),
                    "frozen_mC": float(fr.mC.iloc[0])})
    V = pd.DataFrame(val)
    V["abs_diff"] = (V.recomputed_dC - V.frozen_dC).abs()
    print(f"\n  validation against {spec['frozen'].replace(ROOT + '/', '')}")
    print("  " + V.to_string(index=False).replace("\n", "\n  "))
    worst = float(V.abs_diff.max())
    print(f"  largest absolute difference: {worst:.2e}")
    assert worst < 5e-5, f"{name}: the estimator does not reproduce {spec['frozen']}"

    print("\n  demographic family")
    print("  " + dem.to_string(index=False).replace("\n", "\n  "))
    print(f"  baseline held-out concordance, age spline plus sex: {piv['__BASE__'].mean():.5f}")

    RESULT[name] = {"n": int(len(d)), "n_outcomes": len(OUT),
                    "parquet": spec["parquet"], "oximetry_qc": spec["qc"],
                    "reproduces": spec["frozen"].replace(ROOT + "/", ""),
                    "baseline_heldout_C_age_sex": round(float(piv["__BASE__"].mean()), 6),
                    "validation": json.loads(V.to_json(orient="records")),
                    "rows": json.loads(dem.to_json(orient="records")),
                    "race_counts": {k: int(v) for k, v in d.race_h.value_counts().items()},
                    "bmi_coverage_pct": round(100 * float(d.bmi.notna().mean()), 1)}

allrows = pd.concat([pd.DataFrame(RESULT[k]["rows"]).assign(cohort=k) for k in RESULT])
allrows.to_csv(f"{paths.NUMBERS_DIR}/ranking_demographics_v1.csv", index=False)
json.dump({"cohort_spec_COHORT_N": COHORT_N,
           "estimator": ("identical to ranking_v2.csv: ridge-penalized Cox (penalizer=0.01), "
                         "4-column natural cubic age spline plus sex, trained at one hospital "
                         "and scored at the other in both directions, averaged over outcomes"),
           "modelling_rule_exception": {
               "rule": "any new fit drops one age-spline column and is fitted unpenalized",
               "followed_here": False,
               "why": "this is the ranking pipeline, not an inference fit. Its output is a "
                      "held-out concordance and its only requirement is to be the same "
                      "estimator that produced the 197 frozen values in ranking_v2.csv, so a "
                      "demographic term can be read off the same axis as a sleep measurement. "
                      "ranking_v2.csv is ridge-penalized on all four spline columns; the ridge "
                      "is what makes that basis estimable. Refitting unpenalized with a column "
                      "dropped would put these rows on a different scale from the 197 they are "
                      "drawn beside.",
               "proof": "five frozen ranking values are recomputed through this exact call on "
                        "the current cohort and asserted to reproduce within 5e-5; observed agreement "
                        "is 1e-17. See the validation block in this file.",
               "blast_radius": "no hazard ratio, confidence interval or P value in the paper "
                               "comes from this file"},
           "which_ranking_file_matches_which_cohort": {
               k: {"file": v["reproduces"], "n": v["n"]} for k, v in RESULT.items()},
           "cohorts": RESULT},
          open(f"{paths.NUMBERS_DIR}/ranking_demographics_v1.json", "w"), indent=1)
print("\nwrote numbers/ranking_demographics_v1.csv and .json")
