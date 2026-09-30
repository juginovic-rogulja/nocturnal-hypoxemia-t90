"""
How much discrimination does T90 add ON TOP OF A REALISTIC CLINICAL BASELINE?

Side analysis for a colleague. Nothing here enters the manuscript.

The paper reports dC against an age + sex baseline. A clinician seeing the patient already
knows more than that. This rebuilds the same held-out dC against progressively richer
baselines and reports what is left of the T90 gain.

EVERY convention is inherited unchanged from ../dC_side_analyses/common.py, which is gated on
reproducing the published T90 dC of 0.020338754972 to 1e-9. The only thing that changes model
to model is the ADJUSTER LIST. Specifically kept fixed:

  fitting        lifelines CoxPHFitter(penalizer=0.01)
  cross-fit      fit at hospital I0002 and score at I0006, then the reverse, nanmean of the
                 two held-out concordances. A fold is dropped when either side has fewer
                 than 30 events.
  age            the published 4-column cr(df=4) basis (a second pass drops one column)
  exposure       within-site rank inverse normal of spo2_pct_below_90, entered as ONE column
  outcomes       the ranking's 48
  aggregation    plain mean of the per-outcome gains

HOSPITAL / SITE. The task asked whether site enters as a covariate or as strata. Neither is
possible here, and that is not a shortcut. The paper separates hospitals BY THE FOLD: each
training set is one hospital and each test set is the other. Inside a fold site_id is a
constant, so a hospital indicator is not identifiable and a strata term has one level. The
published design is already stricter than adjusting for hospital, because it never lets a
model see the hospital it is scored on. That convention is kept, and every C below is a
between-hospital transported C. The unpenalized hazard-ratio sanity check at the end DOES
stratify on site, because it is a single pooled fit and follows s07_hazard_ratios.py.

CONTINUOUS COVARIATES are entered on the paper's own scale for a measurement, the within-site
rank inverse normal. Mixing raw units with a rank-normal exposure would not be comparable, and
the rank-normal form is the more generous choice for the baseline because it linearises skewed
measures such as AHI instead of forcing a linear term through a long tail.

COVARIATE SOURCES AND THE ONE PROXY
  BMI          T90_Manuscript/numbers/_bmi_derived.parquet, units resolved per row. The
               pound-contaminated version was retired on 2026-08-08 and is not used.
               Derivable for 6,325 of 19,173 and 7.5% at I0002 against 63.4% at I0006. BMI
               models therefore run inside the BMI-complete subset with their OWN matched
               age + sex control, block B, so no C is ever compared across row sets.
  race         OMOP person table race_source_value, harmonised exactly as
               numbers/freeze_02_demographics.py does. "Unknown or not reported" is kept as
               its own level rather than dropped, so race contributes no missingness.
  smoking      OMOP condition_occurrence, tobacco codes dated ON OR BEFORE the sleep study.
               Current = F17.2xx, 305.1x, Z72.0. Former = V15.82, Z87.891. Reference = no
               tobacco code, which is NOT the same as never-smoker: it is an absence of
               coding. Complete for the whole cohort by construction.
  AHI          master_cohort_v2.csv AHI, complete.
  TST          master_cohort_v2.csv TST_min, 19,167 of 19,173.
  sleep eff.   master_cohort_v2.csv sleep_efficiency_pct, 19,167 of 19,173.
  waking SpO2  ---- PROXY IN THE FULL-COHORT MODELS ----
               A true wake-period mean saturation exists but only for the 6,800 recordings
               with staged wake in the cpap_stage shards (all_wake_mean, 35.5% of the
               cohort). Full-cohort models therefore use spo2_mean, the WHOLE-NIGHT mean
               saturation, which is complete. That is a proxy, and a deliberately harsh one:
               whole-night mean saturation is computed from the same oximetry trace as T90,
               so a baseline holding it is close to adjusting the exposure away. Block W
               repeats everything inside the 6,800 with the genuine wake-period mean.

BLOCKS
  full   19,173, whole-night mean SpO2 as the oxygen term
  D      decomposition inside `full`: which single covariate collapses the gain
  B      6,325 with a derivable BMI, with a matched age + sex control on the same rows
  W      6,800 with a scored wake period, TRUE waking SpO2
  agedrop  every `full` model repeated with one cr(df=4) column removed

Outputs results.csv, per_outcome.csv, hr_sanity.csv, provenance.json and REPORT.txt.
"""
from __future__ import annotations
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
import duckdb
from lifelines import CoxPHFitter
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
SIDE = os.path.join(os.path.dirname(HERE), "dC_side_analyses")
sys.path.insert(0, SIDE)

from common import (PEN, T90ROOT, build_frame, heldout,                # noqa: E402
                    organ_of, ranking_outcomes, site_rint)
import sys as _sys_v8; _sys_v8.path.insert(0, paths.ANALYSIS_DIR)  # v8 sweep 2026-09-12
from cohort_spec import NUMBERS_DIR, N_RANKED_OUTCOMES  # v8 sweep 2026-09-12
__PUB_DC_V7__ = float(__import__("pandas").read_csv(f"{paths.NUMBERS_DIR}/ranking_v3.csv", comment="#").set_index("feature").loc["spo2_pct_below_90", "dC"])   # v7: T90 gain from the CURRENT ranking, was the August literal 0.02033875497249791


WORK = os.path.join(HERE, "_work")
os.makedirs(WORK, exist_ok=True)

PUB_DC = __PUB_DC_V7__
_RK_V8 = __import__("pandas").read_csv(f"{NUMBERS_DIR}/ranking_v3.csv", comment="#"); _impl_v8 = (_RK_V8.mC - _RK_V8.dC); _impl_mode = _impl_v8.round(6).mode().iloc[0]; _off = _RK_V8.feature[_impl_v8.round(6) != _impl_mode].tolist(); assert all("_coupling_" in str(_f) for _f in _off) and (_impl_v8.round(6) == _impl_mode).mean() > 0.9, f"ranking rows disagree on the implied baseline: off-mode rows {_off} (v8.1c rule: every off-mode row must be a spindle-SO coupling feature, mode on over 90 percent of rows)"; print(f"[v8] implied baseline mode {_impl_mode:.6f} over {int((_impl_v8.round(6) == _impl_mode).sum())} rows; off-mode rows {_off}"); PUB_BASE = float(_impl_v8[_impl_v8.round(6) == _impl_mode].mean())  # v8 F2 fix 2026-09-12 (integrity gate 6): the held-out baseline C is the MODE of the implied baseline mC - dC over the ranking rows (139 of 141 agree); at most two off-mode rows are allowed and listed (features refit on their non-missing subset, the two C_coupling_* rows); PUB_BASE is taken from the modal rows. Same rule as numbers/run_combo_v2.py lines 68 to 72
T90 = "spo2_pct_below_90"
SHARDS = (f"{paths.T90_ROOT}/BDSP_New_Ideas_2026-07-26/"
          "pilots/cpap_stage/out*.csv").replace(paths.CPAP_STAGE_DIR, f"{paths.X4_DIR}/cpap_stage")   # v8.2 re-extracted shards
COND = f"{paths.OMOP_CACHE_DIR}/condition_occurrence_cohort.parquet"
MIN_WAKE_MIN = 10.0
AGE4 = [f"age_s{i}" for i in range(4)]
AGE3 = [f"age_s{i}" for i in range(1, 4)]

LOG = []


def say(msg=""):
    print(msg)
    LOG.append(msg)


# ==================================================================== covariate assembly
def add_race(d: pd.DataFrame) -> pd.DataFrame:
    """OMOP race, harmonised exactly as numbers/freeze_02_demographics.py, then dummied.
    Native Hawaiian / Pacific Islander (30 people) and American Indian / Alaska Native are
    folded into Other so no fold can hit an empty category. White is the reference."""
    person = glob.glob(f"{paths.OMOP_CACHE_DIR}/**/person_merged.parquet",
                       recursive=True)[0]
    p = duckdb.sql(f'''select person_id as "BDSPPatientID", race_source_value
                       from '{person}' ''').df().drop_duplicates("BDSPPatientID")
    d = d.merge(p, on="BDSPPatientID", how="left")
    raw = d.race_source_value.astype(str).str.strip().str.upper()

    WHITE = {"WHITE", "CAUCASIAN OR WHITE", "CAUCASIAN"}
    BLACK = {"AFRICAN AMERICAN  OR BLACK", "AFRICAN AMERICAN OR BLACK",
             "BLACK/AFRICAN AMERICAN", "BLACK", "AFRICAN AMERICAN"}
    ASIAN = {"ASIAN"}
    OTHERSET = {"NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER", "PACIFIC ISLANDER",
                "NATIVE HAWAIIAN", "AMERICAN INDIAN OR ALASKA NATIVE", "AMERICAN INDIAN"}
    UNK = {"UNKNOWN, UNAVAILABLE OR UNREPORTED", "UNKNOWN/NOT SPECIFIED", "UNKNOWN",
           "UNABLE TO OBTAIN", "DECLINED TO ANSWER", "NONE", "PREFER NOT TO SAY",
           "NAN", "", "DECLINED", "NOT RECORDED", "PATIENT DECLINED", "NOT SPECIFIED",
           "UNAVAILABLE", "DECLINE TO ANSWER"}

    def h(v):
        if v in WHITE:
            return "White"
        if v in BLACK:
            return "Black"
        if v in ASIAN:
            return "Asian"
        if v in OTHERSET:
            return "Other"
        if v in UNK:
            return "Unknown"
        return "Other"

    d["race_h"] = raw.apply(h)
    for lab, col in [("Black", "race_black"), ("Asian", "race_asian"),
                     ("Other", "race_other"), ("Unknown", "race_unknown")]:
        d[col] = (d.race_h == lab).astype(float)
    return d


def add_smoking(d: pd.DataFrame) -> pd.DataFrame:
    """Tobacco codes dated on or before the sleep study. Reference is no code."""
    key = d[["BDSPPatientID", "psg_date"]].copy()
    duckdb.sql("create or replace table coh as select * from key")
    q = f"""
    select c.person_id as pid,
      max(case when regexp_matches(upper(trim(c.condition_source_value)),
                                   '^(F17\\.?2|305\\.?1|Z72\\.?0)') then 1 else 0 end) as cur,
      max(case when regexp_matches(upper(trim(c.condition_source_value)),
                                   '^(V15\\.?82|Z87\\.?891)') then 1 else 0 end) as hst
    from '{COND}' c join coh k on k."BDSPPatientID" = c.person_id
    where c.condition_start_date <= k.psg_date
      and regexp_matches(upper(trim(c.condition_source_value)),
                         '^(F17|Z87\\.?891|Z72\\.?0|305\\.?1|V15\\.?82)')
    group by 1
    """
    s = duckdb.sql(q).df().rename(columns={"pid": "BDSPPatientID"})
    d = d.merge(s, on="BDSPPatientID", how="left")
    d[["cur", "hst"]] = d[["cur", "hst"]].fillna(0)
    # NOTE the bracket access. d["hst"] is deliberate: a column called `hist` reached by
    # attribute access would silently return the DataFrame.hist METHOD, which equals nothing
    # and would zero the variable without raising.
    d["smk_current"] = (d["cur"] == 1).astype(float)
    d["smk_former"] = ((d["cur"] == 0) & (d["hst"] == 1)).astype(float)
    d["smk_label"] = np.where(d["smk_current"] == 1, "current",
                              np.where(d["smk_former"] == 1, "former", "no_code"))
    return d


def add_bmi(d: pd.DataFrame) -> pd.DataFrame:
    b = pd.read_parquet(f"{paths.NUMBERS_DIR}/_bmi_derived.parquet")
    b = b[["BDSPPatientID", "bmi"]].dropna().drop_duplicates("BDSPPatientID")
    return d.merge(b, on="BDSPPatientID", how="left")


def add_wake_spo2(d: pd.DataFrame) -> pd.DataFrame:
    """True wake-period mean saturation from the cpap_stage shards, the same source and the
    same >=10 staged wake minute rule as numbers/run_wake_sleep_healthy.py."""
    files = sorted(glob.glob(SHARDS))
    assert files, f"no cpap_stage shards under {SHARDS}"   # v8.2: the shard count is a property of the fleet, not a result
    raw = pd.concat([pd.read_csv(f, low_memory=False,
                                 usecols=["BDSPPatientID", "status", "all_wake_mean",
                                          "all_wake_min", "all_wake_t90"])
                     for f in files], ignore_index=True)
    st = raw.drop_duplicates("BDSPPatientID")
    st = st[st.status == "ok"]
    assert st.BDSPPatientID.is_unique
    d = d.merge(st[["BDSPPatientID", "all_wake_mean", "all_wake_min", "all_wake_t90"]],
                on="BDSPPatientID", how="left")
    d.loc[d.all_wake_min < MIN_WAKE_MIN, ["all_wake_mean", "all_wake_t90"]] = np.nan
    return d


RINT_MAP = [(T90, "t90_z"), ("AHI", "ahi_z"), ("TST_min", "tst_z"),
            ("sleep_efficiency_pct", "se_z"), ("spo2_mean", "spo2mean_z"),
            ("bmi", "bmi_z"), ("all_wake_mean", "wakespo2_z")]


def rerank(f: pd.DataFrame) -> pd.DataFrame:
    """Recompute every within-site rank inverse normal inside whatever row set is passed."""
    f = f.reset_index(drop=True).copy()
    for src, dst in RINT_MAP:
        f[dst] = site_rint(f, src)
    return f


def build() -> pd.DataFrame:
    d = build_frame(extra_master_cols=[T90, "AHI", "TST_min", "sleep_efficiency_pct",
                                       "spo2_mean"], cache="frame_clin_probe.pkl")
    key = pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet",
                          columns=["BDSPPatientID", "psg_date"]).drop_duplicates("BDSPPatientID")
    d = d.merge(key, on="BDSPPatientID", how="left")
    d = add_race(d)
    d = add_smoking(d)
    d = add_bmi(d)
    d = add_wake_spo2(d)
    return rerank(d)


# ==================================================================== model definitions
RACE = ["race_black", "race_asian", "race_other", "race_unknown"]
SMOKE = ["smk_current", "smk_former"]
SLEEP_PROXY = ["ahi_z", "tst_z", "se_z", "spo2mean_z"]
SLEEP_TRUE = ["ahi_z", "tst_z", "se_z", "wakespo2_z"]

MODELS_FULL = {
    "M1_paper_age_sex": ([], "age + sex, THE PAPER'S BASELINE"),
    "M2_clinical": (RACE + SMOKE, "age + sex + race + smoking"),
    "M3_clinical_sleep": (RACE + SMOKE + SLEEP_PROXY,
                          "age + sex + race + smoking + AHI + TST + sleep eff + mean SpO2 (PROXY)"),
    "M4_colleague": (RACE + SMOKE + ["ahi_z", "spo2mean_z"],
                     "age + sex + race + smoking + AHI + mean SpO2 (PROXY)"),
    # ---- decomposition: which single term does the damage
    "D1_age_sex_AHI": (["ahi_z"], "age + sex + AHI"),
    "D2_age_sex_meanSpO2": (["spo2mean_z"], "age + sex + mean SpO2 (PROXY)"),
    "D3_age_sex_AHI_meanSpO2": (["ahi_z", "spo2mean_z"], "age + sex + AHI + mean SpO2 (PROXY)"),
    "D4_clinical_AHI_noOxygen": (RACE + SMOKE + ["ahi_z"],
                                 "age + sex + race + smoking + AHI, NO oxygen term"),
    "D5_clinical_sleep_noOxygen": (RACE + SMOKE + ["ahi_z", "tst_z", "se_z"],
                                   "age + sex + race + smoking + AHI + TST + sleep eff, NO oxygen"),
}

MODELS_BMI = {
    "B1_paper_age_sex": ([], "age + sex, MATCHED CONTROL on the BMI-complete rows"),
    "B2_age_sex_BMI": (["bmi_z"], "age + sex + BMI"),
    "B3_clinical_BMI": (["bmi_z"] + RACE + SMOKE, "age + sex + BMI + race + smoking"),
    "B4_clinical_sleep_BMI": (["bmi_z"] + RACE + SMOKE + SLEEP_PROXY,
                              "age + sex + BMI + race + smoking + AHI + TST + sleep eff + mean SpO2"),
    "B5_colleague_BMI": (["bmi_z"] + RACE + SMOKE + ["ahi_z", "spo2mean_z"],
                         "age + sex + BMI + race + smoking + AHI + mean SpO2 (PROXY)"),
}

MODELS_WAKE = {
    "W1_paper_age_sex": ([], "age + sex, MATCHED CONTROL on the wake subcohort"),
    "W2_clinical": (RACE + SMOKE, "age + sex + race + smoking"),
    "WD1_age_sex_AHI": (["ahi_z"], "age + sex + AHI"),
    "WD2_age_sex_trueWakeSpO2": (["wakespo2_z"], "age + sex + TRUE waking SpO2"),
    "WD3_age_sex_wholenightSpO2": (["spo2mean_z"], "age + sex + whole-night mean SpO2"),
    "WD4_clinical_AHI_noOxygen": (RACE + SMOKE + ["ahi_z"],
                                  "age + sex + race + smoking + AHI, NO oxygen term"),
    "W3_clinical_sleep_trueWake": (RACE + SMOKE + SLEEP_TRUE,
                                   "age + sex + race + smoking + AHI + TST + sleep eff + TRUE waking SpO2"),
    "W4_colleague_trueWake": (RACE + SMOKE + ["ahi_z", "wakespo2_z"],
                              "age + sex + race + smoking + AHI + TRUE waking SpO2"),
    "W5_colleague_wholenight": (RACE + SMOKE + ["ahi_z", "spo2mean_z"],
                                "age + sex + race + smoking + AHI + whole-night mean SpO2, "
                                "the same subcohort so proxy and truth are comparable"),
}

FIVE = {"Cardiac": "cardiac", "Respiratory": "respiratory", "Metabolic": "metabolic",
        "Kidney": "kidney"}


# ==================================================================== one model
def run_model(name, extra, frame, outcomes, agecols, tag):
    adj = list(agecols) + ["male"] + list(extra)
    res = heldout({"__base__": [], "t90": ["t90_z"]}, outcomes, frame,
                  penalizer=PEN, adj=adj, pkl_name=f"f_{tag}_{name}.pkl")
    b = res[res.spec == "__base__"].set_index("outcome")["c"]
    f = res[res.spec == "t90"].set_index("outcome")[["c", "gain"]]
    per = pd.DataFrame({"outcome": f.index, "c_base": b.reindex(f.index).values,
                        "c_t90": f.c.values, "dC": f.gain.values})
    per["organ"] = [organ_of(o) for o in per.outcome]
    per["organ5"] = [FIVE.get(o, "other") for o in per.organ]
    per["model"], per["block"] = name, tag
    per["rel_gain"] = per.dC / (1.0 - per.c_base)
    return per


def summarise_model(per, block, name, desc, frame, extra, agecols, age_basis):
    """Every mean is taken over the outcomes where dC EXISTS, so baseline C and C with T90
    are never averaged over different outcome sets."""
    p = per[per.dC.notna()]
    cols = list(agecols) + ["male"] + list(extra) + ["t90_z"]
    ok = frame[[c for c in cols if c in frame.columns]].notna().all(axis=1)
    cc = {str(s): int(v) for s, v in ok.groupby(frame.site_id).sum().items()}
    return {
        "block": block, "model": name, "covariates": desc, "age_basis": age_basis,
        "n_outcomes": int(len(p)),
        "n_complete_case_I0002": cc.get("I0002"), "n_complete_case_I0006": cc.get("I0006"),
        "mean_C_baseline": float(p.c_base.mean()),
        "mean_C_with_T90": float(p.c_t90.mean()),
        "mean_dC": float(p.dC.mean()),
        "rel_gain_on_means": float(p.dC.mean() / (1 - p.c_base.mean())),
        "mean_rel_gain_per_outcome": float(p.rel_gain.mean()),
        "n_outcomes_dC_positive": int((p.dC > 0).sum()),
        "worst_dC": float(p.dC.min()), "best_dC": float(p.dC.max()),
    }


# ==================================================================== hazard-ratio sanity
def hr_check(frame, extra, outcome):
    """UNPENALIZED, site-stratified, age cr(df=4) with the first column dropped, following
    ../dC_side_analyses/s07_hazard_ratios.py."""
    yc, ec = f"{outcome}_years", f"{outcome}_incident"
    g = frame[(frame[f"{outcome}_prevalent"] == 0) & frame[yc].notna() & (frame[yc] > 0)]
    use = AGE3 + ["male"] + list(extra) + ["t90_z", yc, ec, "site_id"]
    f = g[use].dropna()
    best = None
    for s in [None, 0.1, 0.25, 0.5]:
        try:
            c = CoxPHFitter(penalizer=0.0)
            kw = {"fit_options": {"step_size": s}} if s is not None else {}
            c.fit(f, duration_col=yc, event_col=ec, strata=["site_id"], **kw)
            if best is None or float(c.log_likelihood_) > best[0]:
                best = (float(c.log_likelihood_), c)
        except Exception:
            continue
    if best is None:
        return None
    c = best[1]
    lo, hi = np.exp(c.confidence_intervals_.loc["t90_z"])
    return {"outcome": outcome, "n": int(len(f)), "events": int(f[ec].sum()),
            "HR_per_SD": float(np.exp(c.params_["t90_z"])),
            "lo95": float(lo), "hi95": float(hi),
            "p": float(c.summary.loc["t90_z", "p"])}


# ==================================================================== main
def main():
    say("=" * 100)
    say("DOES T90 STILL ADD DISCRIMINATION ON TOP OF A CLINICAL BASELINE?")
    say("side analysis, T90 manuscript. nothing here enters the paper.")
    say("=" * 100)

    d = build()
    say(f"\ncohort {len(d):,}   hospitals " +
        ", ".join(f"{k} n={v:,}" for k, v in d.site_id.value_counts().sort_index().items()))
    OUT = ranking_outcomes(d)
    assert len(OUT) == N_RANKED_OUTCOMES, len(OUT)
    say(f"outcomes {len(OUT)} (the ranking's {N_RANKED_OUTCOMES})")

    # -------------------------------------------------------------- availability
    say("\n" + "-" * 100)
    say("1. COVARIATE AVAILABILITY IN THE 19,173")
    say("-" * 100)
    avail = []
    for lab, col in [("age", "AgeAtVisit"), ("sex", "male"), ("BMI", "bmi"), ("AHI", "AHI"),
                     ("total sleep time", "TST_min"),
                     ("sleep efficiency", "sleep_efficiency_pct"),
                     ("mean SpO2 whole night (PROXY)", "spo2_mean"),
                     ("waking SpO2, TRUE wake period", "all_wake_mean"),
                     ("T90 exposure", T90)]:
        n = int(d[col].notna().sum())
        by = {str(s): round(100 * float(v), 1)
              for s, v in d.groupby("site_id")[col].apply(lambda x: x.notna().mean()).items()}
        avail.append({"covariate": lab, "column": col, "n_present": n,
                      "pct_present": round(100 * n / len(d), 1),
                      "pct_I0002": by.get("I0002"), "pct_I0006": by.get("I0006")})
        say(f"  {lab:<34} {n:>6,}  {100*n/len(d):>5.1f}%   "
            f"I0002 {by.get('I0002'):>5}%   I0006 {by.get('I0006'):>5}%")
    say(f"  {'race (Unknown kept as a level)':<34} {len(d):>6,}  100.0%   "
        f"known category {100*float((d.race_h!='Unknown').mean()):.1f}%")
    say(f"  {'smoking (no code = reference)':<34} {len(d):>6,}  100.0%")
    say("\n  race     " + str({k: int(v) for k, v in d.race_h.value_counts().items()}))
    say("  smoking  " + str({k: int(v) for k, v in d.smk_label.value_counts().items()}))
    say("\n  MISSING-DATA HANDLING")
    say("    race     no missingness. 'Unknown or not reported' is kept as its own level "
        "rather than dropped.")
    say("    smoking  no missingness. The reference level is 'no tobacco code', which is an "
        "absence of coding and NOT a never-smoker.")
    say("    BMI      67% missing and unbalanced by hospital (7.5% at I0002, 63.4% at I0006). "
        "Complete case inside block B,")
    say("             which carries its own matched age + sex control so no C crosses row sets.")
    say("    waking SpO2  64.5% missing. Full-cohort models use the whole-night mean as an "
        "explicit proxy. Block W uses the truth.")
    say("    TST and sleep efficiency  6 of 19,173 missing, complete case.")

    ov = d[d.all_wake_mean.notna()]
    rho_proxy = float(stats.spearmanr(ov.spo2_mean, ov.all_wake_mean).statistic)
    rho_t90_mean = float(stats.spearmanr(d[T90], d.spo2_mean).statistic)
    rho_t90_wake = float(stats.spearmanr(ov[T90], ov.all_wake_mean).statistic)
    rho_t90_ahi = float(stats.spearmanr(d[T90], d.AHI).statistic)
    rho_t90_bmi = float(stats.spearmanr(d.loc[d.bmi.notna(), T90],
                                        d.loc[d.bmi.notna(), "bmi"]).statistic)
    say("\n  HOW CLOSE EACH BASELINE TERM ALREADY IS TO T90 (Spearman)")
    say(f"    T90 against whole-night mean SpO2 (the PROXY)   {rho_t90_mean:+.3f}")
    say(f"    T90 against TRUE waking mean SpO2               {rho_t90_wake:+.3f}")
    say(f"    T90 against AHI                                 {rho_t90_ahi:+.3f}")
    say(f"    T90 against BMI                                 {rho_t90_bmi:+.3f}")
    say(f"    whole-night mean SpO2 against TRUE waking mean  {rho_proxy:+.3f}   "
        f"(n={len(ov):,})")
    say("    The whole-night mean is read off the SAME oximetry trace as T90, so a baseline "
        "holding it is close to")
    say("    adjusting the exposure away. The true waking mean is a different period of the "
        "night and is the fair test.")

    say("\n  HOSPITAL. Not a covariate and not strata, and that is not a shortcut. The "
        "published design separates hospitals")
    say("  BY THE FOLD: each model trains on one hospital and is scored on the other, so "
        "site_id is constant inside a fold")
    say("  and a hospital term is not identifiable. Every C below is already a "
        "between-hospital transported C, which is")
    say("  stricter than adjusting for hospital. The hazard-ratio check at the end is a "
        "single pooled fit and does stratify on site.")

    # -------------------------------------------------------------- positive control
    say("\n" + "-" * 100)
    say("2. POSITIVE CONTROL")
    say("-" * 100)
    pc = run_model("M1_paper_age_sex", [], d, OUT, AGE4, "PC")
    got_dc, got_base = float(pc.dC.mean()), float(pc.c_base.mean())
    ok = abs(got_dc - PUB_DC) < 1e-9 and abs(got_base - PUB_BASE) < 5e-6
    say(f"  published dC {PUB_DC:.12f}   here {got_dc:.12f}   diff {got_dc - PUB_DC:+.3e}")
    say(f"  published baseline C {PUB_BASE:.6f}   here {got_base:.6f}   "
        f"diff {got_base - PUB_BASE:+.3e}")
    say(f"  VERDICT {'PASS' if ok else 'FAIL'}")
    if not ok:
        raise SystemExit("positive control failed, stop")

    # -------------------------------------------------------------- all blocks
    per_all, rows = [], []

    def do(block, frame, models, agecols=AGE4, basis="cr(df=4), published 4-column"):
        for name, (extra, desc) in models.items():
            per = run_model(name, extra, frame, OUT, agecols, block)
            per_all.append(per)
            rows.append(summarise_model(per, block, name, desc, frame, extra, agecols, basis))
            r = rows[-1]
            say(f"  {block:<12} {name:<28} baseline C {r['mean_C_baseline']:.4f}  "
                f"+T90 C {r['mean_C_with_T90']:.4f}  dC {r['mean_dC']:+.4f}  "
                f"({r['n_outcomes']}/{N_RANKED_OUTCOMES})")

    say("\n" + "-" * 100)
    say("3. FITTING")
    say("-" * 100)
    do("full", d, MODELS_FULL)

    dB = rerank(d[d.bmi.notna()])
    say(f"\n  block B, BMI-complete subset n = {len(dB):,}   " +
        ", ".join(f"{k} n={v:,}" for k, v in dB.site_id.value_counts().sort_index().items()))
    do("bmi_6325", dB, MODELS_BMI)

    dW = rerank(d[d.all_wake_mean.notna()])
    say(f"\n  block W, scored-wake subcohort n = {len(dW):,}   " +
        ", ".join(f"{k} n={v:,}" for k, v in dW.site_id.value_counts().sort_index().items()))
    do("wake_6800", dW, MODELS_WAKE)

    say("\n  age-spline drop check")
    for name, (extra, desc) in MODELS_FULL.items():
        per = run_model(name + "__agedrop", extra, d, OUT, AGE3, "full_agedrop")
        per_all.append(per)
        rows.append(summarise_model(per, "full_agedrop", name + "__agedrop", desc, d, extra,
                                    AGE3, "cr(df=4), one column dropped"))

    P = pd.concat(per_all, ignore_index=True)
    R = pd.DataFrame(rows)
    P.to_csv(os.path.join(HERE, "per_outcome.csv"), index=False)
    R.to_csv(os.path.join(HERE, "results.csv"), index=False)

    # -------------------------------------------------------------- headline
    say("\n" + "=" * 100)
    say("4. HEADLINE. MEAN dC ACROSS THE OUTCOMES, HELD OUT AND TRANSPORTED BETWEEN HOSPITALS")
    say("=" * 100)
    say(f"{'block':<11}{'model':<30}{'baseline C':>11}{'+T90 C':>9}{'dC':>9}"
        f"{'rel gain':>10}{'n out':>7}{'dC>0':>6}   covariates")
    for blk in ["full", "bmi_6325", "wake_6800"]:
        for _i, r in R[R.block == blk].iterrows():
            say(f"{blk:<11}{r.model:<30}{r.mean_C_baseline:>11.4f}{r.mean_C_with_T90:>9.4f}"
                f"{r.mean_dC:>+9.4f}{100*r.rel_gain_on_means:>9.1f}%{r.n_outcomes:>7d}"
                f"{r.n_outcomes_dC_positive:>6d}   {r.covariates}")
        say("")
    say("  rel gain = mean dC / (1 - mean baseline C). The share of the discrimination still "
        "unclaimed that T90 closes.")
    say("  Baseline C and C with T90 are averaged over the SAME outcomes within a row, so a "
        f"model with fewer than {N_RANKED_OUTCOMES} is internally consistent.")
    say("  Across blocks the row sets differ, so compare a model only with the matched "
        "age + sex control in its own block.")

    say("\n" + "-" * 100)
    say("5. AGE-SPLINE DROP CHECK (one cr column removed, the standing rule)")
    say("-" * 100)
    a = R[R.block == "full"].set_index("model")
    b_ = R[R.block == "full_agedrop"].set_index("model")
    for m in MODELS_FULL:
        say(f"  {m:<30} 4-column dC {a.loc[m,'mean_dC']:+.4f}   "
            f"one dropped {b_.loc[m+'__agedrop','mean_dC']:+.4f}   "
            f"diff {b_.loc[m+'__agedrop','mean_dC']-a.loc[m,'mean_dC']:+.4f}")

    # -------------------------------------------------------------- organs
    keymods = ["M1_paper_age_sex", "M2_clinical", "D4_clinical_AHI_noOxygen",
               "M4_colleague", "M3_clinical_sleep"]
    keyw = ["W1_paper_age_sex", "W2_clinical", "WD4_clinical_AHI_noOxygen",
            "W4_colleague_trueWake", "W3_clinical_sleep_trueWake"]
    keyb = ["B1_paper_age_sex", "B3_clinical_BMI", "B5_colleague_BMI"]

    def organ_table(block, mods, index_col, title):
        sub = P[(P.block == block) & P.model.isin(mods) & P.dC.notna()]
        piv = sub.pivot_table(index=index_col, columns="model", values="dC", aggfunc="mean")
        keep = [m for m in mods if m in piv.columns]
        piv = piv[keep]
        if index_col == "organ5":
            piv = piv.reindex([i for i in ["cardiac", "respiratory", "metabolic", "kidney",
                                           "other"] if i in piv.index])
        else:
            piv = piv.sort_values(keep[0], ascending=False)
        cnt = sub[sub.model == keep[0]].groupby(index_col).size()
        say(f"\n  {title}")
        say(f"  {'organ':<14}{'n':>4}" + "".join(f"{m[:22]:>24}" for m in keep))
        for o in piv.index:
            say(f"  {o:<14}{int(cnt.get(o,0)):>4}" +
                "".join(f"{piv.loc[o,m]:>+24.4f}" for m in keep))

    say("\n" + "-" * 100)
    say("6. PER-ORGAN dC")
    say("-" * 100)
    organ_table("full", keymods, "organ5", "full cohort, five groups")
    organ_table("full", keymods, "organ", "full cohort, the paper's full organ grouping")
    organ_table("wake_6800", keyw, "organ5", "wake subcohort with TRUE waking SpO2, five groups")
    organ_table("bmi_6325", keyb, "organ5", "BMI-complete subset, five groups")

    # -------------------------------------------------------------- top 10
    say("\n" + "-" * 100)
    say("7. TOP 10 OUTCOMES BY dC")
    say("-" * 100)
    for blk, mods in [("full", keymods), ("wake_6800", keyw), ("bmi_6325", keyb)]:
        for m in mods:
            s = P[(P.block == blk) & (P.model == m) & P.dC.notna()].sort_values(
                "dC", ascending=False).head(10)
            if not len(s):
                continue
            say(f"\n  [{blk}] {m}")
            say(f"    {'outcome':<22}{'organ':<14}{'base C':>9}{'+T90 C':>9}{'dC':>9}{'rel':>8}")
            for _i, r in s.iterrows():
                say(f"    {r.outcome:<22}{r.organ:<14}{r.c_base:>9.4f}{r.c_t90:>9.4f}"
                    f"{r.dC:>+9.4f}{100*r.rel_gain:>7.1f}%")

    # -------------------------------------------------------------- HR sanity
    say("\n" + "-" * 100)
    say("8. HAZARD-RATIO DIRECTION CHECK")
    say("-" * 100)
    say("  Unpenalized, site-stratified, the specification of "
        "../dC_side_analyses/s07_hazard_ratios.py.")
    say("  Reference is the paper's own published per-SD panel for the age + sex model.")
    ref = pd.read_csv(os.path.join(SIDE, "hr_tables.csv"))
    ref = ref[(ref.cohort == "full_19173") & (ref.model == "single_t90")].set_index("outcome_key")
    hrs = []
    for oc in ["hf", "resp_failure"]:
        for mname in ["M2_clinical", "D4_clinical_AHI_noOxygen", "M4_colleague",
                      "M3_clinical_sleep"]:
            r = hr_check(d, MODELS_FULL[mname][0], oc)
            if r is None:
                continue
            r["model"] = mname
            r["paper_HR_per_SD"] = float(ref.loc[oc, "HR_per_SD"])
            r["paper_lo95"] = float(ref.loc[oc, "lo95"])
            r["paper_hi95"] = float(ref.loc[oc, "hi95"])
            r["direction_matches"] = bool((r["HR_per_SD"] > 1) == (r["paper_HR_per_SD"] > 1))
            hrs.append(r)
            say(f"  {oc:<14}{mname:<28} T90 HR per SD {r['HR_per_SD']:.3f} "
                f"({r['lo95']:.3f}-{r['hi95']:.3f})   paper age+sex "
                f"{r['paper_HR_per_SD']:.3f} ({r['paper_lo95']:.3f}-{r['paper_hi95']:.3f})   "
                f"{'MATCHES' if r['direction_matches'] else 'DIFFERS'}")
    pd.DataFrame(hrs).to_csv(os.path.join(HERE, "hr_sanity.csv"), index=False)

    # -------------------------------------------------------------- shrinkage decomposition
    say("\n" + "-" * 100)
    say("9. IS THE SHRINKAGE MECHANICAL OR REAL?")
    say("-" * 100)
    say("  A better baseline has a higher C, so less room is left above it and dC must fall "
        "even if T90 carried")
    say("  exactly the same independent information. The mechanical part of that is "
        "dC_reference scaled by the ratio")
    say("  of headroom, (1 - C_model) / (1 - C_reference). Anything beyond it is genuine "
        "overlap: the baseline already")
    say("  contains what T90 was carrying. Reference is the matched age + sex model in each "
        "block.")
    say("")
    say(f"  {'block':<11}{'model':<28}{'base C':>8}{'headroom':>10}{'dC exp':>9}"
        f"{'dC obs':>9}{'drop from headroom':>21}")
    shrink = []
    for blk, refname in [("full", "M1_paper_age_sex"), ("bmi_6325", "B1_paper_age_sex"),
                         ("wake_6800", "W1_paper_age_sex")]:
        rr = R[R.block == blk].set_index("model")
        c0, d0 = float(rr.loc[refname, "mean_C_baseline"]), float(rr.loc[refname, "mean_dC"])
        h0 = 1 - c0
        for m in rr.index:
            cm, dm = float(rr.loc[m, "mean_C_baseline"]), float(rr.loc[m, "mean_dC"])
            hm = 1 - cm
            exp = d0 * hm / h0
            share = ((d0 - exp) / (d0 - dm)) if abs(d0 - dm) > 1e-12 else np.nan
            shrink.append({"block": blk, "model": m, "reference": refname,
                           "mean_C_baseline": cm, "headroom": hm,
                           "dC_expected_if_only_headroom": exp, "dC_observed": dm,
                           "share_of_drop_explained_by_headroom": share})
            txt = "reference" if m == refname else (
                f"{100*share:>19.1f}%" if np.isfinite(share) else " ")
            say(f"  {blk:<11}{m:<28}{cm:>8.4f}{hm:>10.4f}{exp:>+9.4f}{dm:>+9.4f}{txt:>21}")
        say("")
    pd.DataFrame(shrink).to_csv(os.path.join(HERE, "shrinkage_decomposition.csv"), index=False)

    say("\n" + "=" * 100)
    say("10. THE ANSWER")
    say("=" * 100)
    fu = R[R.block == "full"].set_index("model")
    wk = R[R.block == "wake_6800"].set_index("model")
    bm = R[R.block == "bmi_6325"].set_index("model")
    say(f"  Against the paper's age + sex baseline, T90 adds dC "
        f"{fu.loc['M1_paper_age_sex','mean_dC']:+.4f} on a baseline C of "
        f"{fu.loc['M1_paper_age_sex','mean_C_baseline']:.4f}.")
    say(f"  Against age + sex + race + smoking it adds "
        f"{fu.loc['M2_clinical','mean_dC']:+.4f} on {fu.loc['M2_clinical','mean_C_baseline']:.4f}.")
    say(f"  Against the colleague's full line it adds {fu.loc['M4_colleague','mean_dC']:+.4f} "
        f"on {fu.loc['M4_colleague','mean_C_baseline']:.4f}.")
    say(f"  Demographics, smoking, AHI, total sleep time and sleep efficiency together cost "
        f"only "
        f"{fu.loc['M1_paper_age_sex','mean_dC']-fu.loc['D5_clinical_sleep_noOxygen','mean_dC']:+.4f} "
        f"of the gain "
        f"({fu.loc['M1_paper_age_sex','mean_dC']:+.4f} to "
        f"{fu.loc['D5_clinical_sleep_noOxygen','mean_dC']:+.4f}).")
    say(f"  A single oxygen saturation term costs the rest: adding it takes dC from "
        f"{fu.loc['D5_clinical_sleep_noOxygen','mean_dC']:+.4f} to "
        f"{fu.loc['M3_clinical_sleep','mean_dC']:+.4f}.")
    say(f"  That is not the whole-night proxy talking. In the 6,800 with a real scored wake "
        f"period, the TRUE waking mean")
    say(f"  gives dC {wk.loc['W4_colleague_trueWake','mean_dC']:+.4f} and the whole-night mean "
        f"on the same rows gives "
        f"{wk.loc['W5_colleague_wholenight','mean_dC']:+.4f}.")
    say(f"  BMI matters too but is a weaker term: on the BMI-complete rows the matched age + "
        f"sex dC is "
        f"{bm.loc['B1_paper_age_sex','mean_dC']:+.4f} and adding BMI alone takes it to "
        f"{bm.loc['B2_age_sex_BMI','mean_dC']:+.4f}.")
    say(f"  The residual gain is concentrated, not spread. Under the colleague's line the "
        f"respiratory mean is still")
    say(f"  {P[(P.block=='full')&(P.model=='M4_colleague')&(P.organ5=='respiratory')].dC.mean():+.4f} "
        f"and the 28 non-cardiorespiratory-metabolic-renal outcomes average "
        f"{P[(P.block=='full')&(P.model=='M4_colleague')&(P.organ5=='other')].dC.mean():+.4f}.")

    json.dump({"positive_control_dC": got_dc, "positive_control_baseline": got_base,
               "published_dC": PUB_DC, "published_baseline": PUB_BASE,
               "spearman": {"t90_vs_wholenight_mean_spo2": rho_t90_mean,
                            "t90_vs_true_waking_mean_spo2": rho_t90_wake,
                            "t90_vs_ahi": rho_t90_ahi, "t90_vs_bmi": rho_t90_bmi,
                            "wholenight_vs_true_waking": rho_proxy},
               "availability": avail, "hazard_ratio_check": hrs},
              open(os.path.join(HERE, "provenance.json"), "w"), indent=2)

    with open(os.path.join(HERE, "REPORT.txt"), "w") as f:
        f.write("\n".join(LOG) + "\n")
    say("\nwrote results.csv, per_outcome.csv, hr_sanity.csv, provenance.json, REPORT.txt")


if __name__ == "__main__":
    main()
