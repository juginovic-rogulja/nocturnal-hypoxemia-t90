"""
THE SAME BASELINE LADDER AS run_clinical_baselines.py, BUT THE ADDED EXPOSURE IS TOTAL SLEEP TIME.

Side analysis. Nothing here enters the manuscript. No manuscript or frozen file is touched.

WHAT IS REUSED, UNCHANGED
  ../dC_side_analyses/common.py            build_frame, heldout, site_rint, organ_of,
                                           ranking_outcomes, PEN. Gated on reproducing the
                                           published T90 dC to 1e-9.
  ./run_clinical_baselines.py              add_race, add_smoking, add_bmi. Imported, not copied,
                                           so the covariate definitions cannot drift.

WHAT CHANGES
  the exposure column only. spo2_pct_below_90 -> TST_min, both entered as ONE column holding the
  within-site rank inverse normal of the measure, so every dC is "per 1 SD" on the paper's scale.

KEPT FIXED
  fitting        lifelines CoxPHFitter(penalizer=0.01)
  cross-fit      fit at hospital I0002 and score at I0006, then the reverse, nanmean of the two
                 held-out concordances. A fold is dropped when either side has under 30 events.
  age            the published 4-column cr(df=4) basis. A second pass drops one column.
  outcomes       the ranking's 48
  aggregation    plain mean of the per-outcome gains

TWO POSITIVE CONTROLS, BOTH BEFORE ANY LADDER RUNS
  TST   L1 must reproduce numbers/ranking_v3.csv row 192, dC -0.000956963972 on baseline 0.647540
  T90   the same L1 spec with T90 as the exposure must reproduce 0.020338754972 on 0.647540

THE TST-AS-COVARIATE ASYMMETRY
  The T90 ladder's rich rungs (M3_clinical_sleep, D5_clinical_sleep_noOxygen, B4_clinical_sleep_BMI)
  hold TST as a COVARIATE. That covariate cannot survive here because TST is the exposure. So the
  T90 side is re-run on the TST-free version of those covariate sets (suffix _noTSTcov below) and
  the comparison table quotes THOSE, so both exposures face an identical baseline. The original
  T90 rows are kept in the table for reference and flagged.

  Note the task's L4 and L4b coincide. L4 is "age + sex + race + smoking + AHI + sleep efficiency".
  L4b is "the T90 ladder's equivalent covariate set minus TST", and the T90 ladder's equivalent
  row D5_clinical_sleep_noOxygen is age + sex + race + smoking + AHI + TST + sleep efficiency.
  Removing TST from it gives exactly L4. They are reported as one row and the identity is stated
  rather than manufactured into a difference.

OXYGEN ROWS
  Every rung holding mean SpO2 is prefixed OX and flagged. Alen has said these are not the rows
  he cares about. They are run so the ladder is complete, not because they carry the argument.

Outputs results_tst.csv, per_outcome_tst.csv, provenance_tst.json, REPORT_TST.txt.
The comparison deliverables are built by make_comparison.py from these files plus results.csv.
"""
from __future__ import annotations
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths

import json
import os
import sys
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SIDE = os.path.join(os.path.dirname(HERE), "dC_side_analyses")
sys.path.insert(0, SIDE)
sys.path.insert(0, HERE)

from common import (PEN, T90ROOT, build_frame, heldout,                 # noqa: E402
                    organ_of, ranking_outcomes, site_rint)
from run_clinical_baselines import add_bmi, add_race, add_smoking       # noqa: E402
import sys as _sys_v8; _sys_v8.path.insert(0, paths.ANALYSIS_DIR)  # v8 sweep 2026-09-12
from cohort_spec import NUMBERS_DIR, N_RANKED_OUTCOMES  # v8 sweep 2026-09-12
__PUB_DC_V7__ = float(__import__("pandas").read_csv(f"{paths.NUMBERS_DIR}/ranking_v3.csv", comment="#").set_index("feature").loc["spo2_pct_below_90", "dC"])   # v7: T90 gain from the CURRENT ranking, was the August literal 0.02033875497249791


T90 = "spo2_pct_below_90"
TST = "TST_min"

# published, from T90_Manuscript/numbers/ranking_v3.csv
PUB_T90_DC = __PUB_DC_V7__
_RK_V8 = __import__("pandas").read_csv(f"{NUMBERS_DIR}/ranking_v3.csv", comment="#"); _impl_v8 = (_RK_V8.mC - _RK_V8.dC); _impl_mode = _impl_v8.round(6).mode().iloc[0]; _off = _RK_V8.feature[_impl_v8.round(6) != _impl_mode].tolist(); assert all("_coupling_" in str(_f) for _f in _off) and (_impl_v8.round(6) == _impl_mode).mean() > 0.9, f"ranking rows disagree on the implied baseline: off-mode rows {_off} (v8.1c rule: every off-mode row must be a spindle-SO coupling feature, mode on over 90 percent of rows)"; print(f"[v8] implied baseline mode {_impl_mode:.6f} over {int((_impl_v8.round(6) == _impl_mode).sum())} rows; off-mode rows {_off}"); PUB_T90_BASE = float(_impl_v8[_impl_v8.round(6) == _impl_mode].mean())  # v8 F2 fix 2026-09-12 (integrity gate 6): the held-out baseline C is the MODE of the implied baseline mC - dC over the ranking rows (139 of 141 agree); at most two off-mode rows are allowed and listed (features refit on their non-missing subset, the two C_coupling_* rows); PUB_T90_BASE is taken from the modal rows. Same rule as numbers/run_combo_v2.py lines 68 to 72
__rk_v7 = __import__("pandas").read_csv(f"{paths.NUMBERS_DIR}/ranking_v3.csv", comment="#").set_index("feature")   # v7: TST reference from the CURRENT ranking (was rank 192, dC -0.000957)
PUB_TST_DC = float(__rk_v7.loc["TST_min", "dC"])
PUB_TST_CWITH = float(__rk_v7.loc["TST_min", "mC"])
PUB_TST_RANK = int(__rk_v7.loc["TST_min", "rank"])

AGE4 = [f"age_s{i}" for i in range(4)]
AGE3 = [f"age_s{i}" for i in range(1, 4)]

LOG = []


def say(msg=""):
    print(msg, flush=True)
    LOG.append(msg)


# ==================================================================== frame
# The T90 script's RINT_MAP minus all_wake_mean. No wake block is requested here, so the 72
# cpap_stage shards are not read. site_rint is per column and independent, so dropping an unused
# column from the map cannot change any other column's transform.
RINT_MAP = [(T90, "t90_z"), ("AHI", "ahi_z"), (TST, "tst_z"),
            ("sleep_efficiency_pct", "se_z"), ("spo2_mean", "spo2mean_z"), ("bmi", "bmi_z")]


def rerank(f: pd.DataFrame) -> pd.DataFrame:
    """Recompute every within-site rank inverse normal inside whatever row set is passed."""
    f = f.reset_index(drop=True).copy()
    for src, dst in RINT_MAP:
        f[dst] = site_rint(f, src)
    return f


def build() -> pd.DataFrame:
    d = build_frame(extra_master_cols=[T90, "AHI", TST, "sleep_efficiency_pct", "spo2_mean"],
                    cache="frame_clin_probe.pkl")
    key = pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet",
                          columns=["BDSPPatientID", "psg_date"]).drop_duplicates("BDSPPatientID")
    d = d.merge(key, on="BDSPPatientID", how="left")
    d = add_race(d)
    d = add_smoking(d)
    d = add_bmi(d)
    return rerank(d)


# ==================================================================== models
RACE = ["race_black", "race_asian", "race_other", "race_unknown"]
SMOKE = ["smk_current", "smk_former"]

# (name, covariates beyond age+sex, description, is_oxygen_row, T90-ladder counterpart)
LADDER_FULL = [
    ("L1_age_sex", [],
     "age + sex, THE PAPER'S BASELINE", False, "M1_paper_age_sex"),
    ("L2_age_sex_AHI", ["ahi_z"],
     "age + sex + AHI", False, "D1_age_sex_AHI"),
    ("L3_age_sex_race_smoke", RACE + SMOKE,
     "age + sex + race + smoking", False, "M2_clinical"),
    ("LX_clinical_AHI", RACE + SMOKE + ["ahi_z"],
     "age + sex + race + smoking + AHI", False, "D4_clinical_AHI_noOxygen"),
    ("L4_clinical_AHI_SE", RACE + SMOKE + ["ahi_z", "se_z"],
     "age + sex + race + smoking + AHI + sleep efficiency (= L4 = L4b, no TST covariate)",
     False, "D5_clinical_sleep_noOxygen__noTSTcov"),
    ("OX1_age_sex_meanSpO2", ["spo2mean_z"],
     "age + sex + mean SpO2", True, "D2_age_sex_meanSpO2"),
    ("OX2_age_sex_AHI_meanSpO2", ["ahi_z", "spo2mean_z"],
     "age + sex + AHI + mean SpO2", True, "D3_age_sex_AHI_meanSpO2"),
    ("OX3_clinical_AHI_meanSpO2", RACE + SMOKE + ["ahi_z", "spo2mean_z"],
     "age + sex + race + smoking + AHI + mean SpO2", True, "M4_colleague"),
    ("OX4_clinical_AHI_SE_meanSpO2", RACE + SMOKE + ["ahi_z", "se_z", "spo2mean_z"],
     "age + sex + race + smoking + AHI + sleep efficiency + mean SpO2",
     True, "M3_clinical_sleep__noTSTcov"),
]

LADDER_BMI = [
    ("B1_age_sex", [],
     "age + sex, MATCHED CONTROL on the BMI-complete rows", False, "B1_paper_age_sex"),
    ("B2_age_sex_BMI", ["bmi_z"],
     "age + sex + BMI", False, "B2_age_sex_BMI"),
    ("B3_clinical_BMI", ["bmi_z"] + RACE + SMOKE,
     "age + sex + BMI + race + smoking", False, "B3_clinical_BMI"),
    ("B4_clinical_BMI_AHI_SE", ["bmi_z"] + RACE + SMOKE + ["ahi_z", "se_z"],
     "age + sex + BMI + race + smoking + AHI + sleep efficiency", False,
     "B4_clinical_sleep_BMI__noTSTcov_noOx"),
    ("OXB1_colleague_BMI", ["bmi_z"] + RACE + SMOKE + ["ahi_z", "spo2mean_z"],
     "age + sex + BMI + race + smoking + AHI + mean SpO2", True, "B5_colleague_BMI"),
    ("OXB2_clinical_sleep_BMI", ["bmi_z"] + RACE + SMOKE + ["ahi_z", "se_z", "spo2mean_z"],
     "age + sex + BMI + race + smoking + AHI + sleep efficiency + mean SpO2", True,
     "B4_clinical_sleep_BMI__noTSTcov"),
]

# T90 rows that must be produced fresh, because the covariate set is not in results.csv.
# Suffix __noTSTcov means the published T90 rung with its TST covariate removed.
T90_NEW_FULL = [
    ("D5_clinical_sleep_noOxygen__noTSTcov", RACE + SMOKE + ["ahi_z", "se_z"],
     "age + sex + race + smoking + AHI + sleep efficiency, NO oxygen, TST covariate REMOVED"),
    ("M3_clinical_sleep__noTSTcov", RACE + SMOKE + ["ahi_z", "se_z", "spo2mean_z"],
     "age + sex + race + smoking + AHI + sleep efficiency + mean SpO2, TST covariate REMOVED"),
]
T90_NEW_BMI = [
    ("B4_clinical_sleep_BMI__noTSTcov_noOx", ["bmi_z"] + RACE + SMOKE + ["ahi_z", "se_z"],
     "age + sex + BMI + race + smoking + AHI + sleep efficiency, TST covariate REMOVED"),
    ("B4_clinical_sleep_BMI__noTSTcov", ["bmi_z"] + RACE + SMOKE + ["ahi_z", "se_z", "spo2mean_z"],
     "age + sex + BMI + race + smoking + AHI + sleep eff + mean SpO2, TST covariate REMOVED"),
]
# T90 rows re-run purely to prove this script's machinery equals the one that wrote results.csv.
T90_VERIFY_FULL = [("M1_paper_age_sex", []), ("M2_clinical", RACE + SMOKE),
                   ("D5_clinical_sleep_noOxygen", RACE + SMOKE + ["ahi_z", "tst_z", "se_z"])]
T90_VERIFY_BMI = [("B1_paper_age_sex", [])]

FIVE = {"Cardiac": "cardiac", "Respiratory": "respiratory", "Metabolic": "metabolic",
        "Kidney": "kidney"}


# ==================================================================== one model
def run_model(name, extra, frame, outcomes, agecols, tag, expo):
    """expo is the single added exposure column, 'tst_z' or 't90_z'."""
    adj = list(agecols) + ["male"] + list(extra)
    res = heldout({"__base__": [], "expo": [expo]}, outcomes, frame,
                  penalizer=PEN, adj=adj, pkl_name=f"tstf_{tag}_{name}_{expo}.pkl")
    b = res[res.spec == "__base__"].set_index("outcome")["c"]
    f = res[res.spec == "expo"].set_index("outcome")[["c", "gain"]]
    per = pd.DataFrame({"outcome": f.index, "c_base": b.reindex(f.index).values,
                        "c_expo": f.c.values, "dC": f.gain.values})
    per["organ"] = [organ_of(o) for o in per.outcome]
    per["organ5"] = [FIVE.get(o, "other") for o in per.organ]
    per["model"], per["block"], per["exposure"] = name, tag, expo
    per["rel_gain"] = per.dC / (1.0 - per.c_base)
    return per


def summarise_model(per, block, name, desc, frame, extra, agecols, age_basis, expo,
                    oxygen, counterpart):
    p = per[per.dC.notna()]
    cols = list(agecols) + ["male"] + list(extra) + [expo]
    ok = frame[[c for c in cols if c in frame.columns]].notna().all(axis=1)
    cc = {str(s): int(v) for s, v in ok.groupby(frame.site_id).sum().items()}
    return {
        "block": block, "model": name, "exposure": expo, "covariates": desc,
        "oxygen_row": bool(oxygen), "t90_counterpart": counterpart, "age_basis": age_basis,
        "n_outcomes": int(len(p)),
        "n_complete_case_I0002": cc.get("I0002"), "n_complete_case_I0006": cc.get("I0006"),
        "mean_C_baseline": float(p.c_base.mean()),
        "mean_C_with_exposure": float(p.c_expo.mean()),
        "mean_dC": float(p.dC.mean()),
        "rel_gain_on_means": float(p.dC.mean() / (1 - p.c_base.mean())),
        "mean_rel_gain_per_outcome": float(p.rel_gain.mean()),
        "n_outcomes_dC_positive": int((p.dC > 0).sum()),
        "worst_dC": float(p.dC.min()), "best_dC": float(p.dC.max()),
    }


# ==================================================================== main
def main():
    say("=" * 100)
    say("HOW MUCH DISCRIMINATION DOES TOTAL SLEEP TIME ADD ON TOP OF A CLINICAL BASELINE?")
    say("The T90 ladder of run_clinical_baselines.py, rerun with TST_min as the added exposure.")
    say("Side analysis. Nothing here enters the paper.")
    say("=" * 100)

    d = build()
    say(f"\ncohort {len(d):,}   hospitals " +
        ", ".join(f"{k} n={v:,}" for k, v in d.site_id.value_counts().sort_index().items()))
    OUT = ranking_outcomes(d)
    assert len(OUT) == N_RANKED_OUTCOMES, len(OUT)
    say(f"outcomes {len(OUT)} (the ranking's {N_RANKED_OUTCOMES})")

    # -------------------------------------------------------------- availability
    say("\n" + "-" * 100)
    say("1. EXPOSURE AND COVARIATE AVAILABILITY IN THE 19,173")
    say("-" * 100)
    avail = []
    for lab, col in [("age", "AgeAtVisit"), ("sex", "male"), ("BMI", "bmi"), ("AHI", "AHI"),
                     ("TOTAL SLEEP TIME, the exposure", TST),
                     ("sleep efficiency", "sleep_efficiency_pct"),
                     ("mean SpO2 whole night", "spo2_mean"),
                     ("T90, the other exposure", T90)]:
        n = int(d[col].notna().sum())
        by = {str(s): round(100 * float(v), 1)
              for s, v in d.groupby("site_id")[col].apply(lambda x: x.notna().mean()).items()}
        avail.append({"covariate": lab, "column": col, "n_present": n,
                      "pct_present": round(100 * n / len(d), 1),
                      "pct_I0002": by.get("I0002"), "pct_I0006": by.get("I0006")})
        say(f"  {lab:<34} {n:>6,}  {100*n/len(d):>5.1f}%   "
            f"I0002 {by.get('I0002'):>5}%   I0006 {by.get('I0006'):>5}%")
    say(f"  {'race (Unknown kept as a level)':<34} {len(d):>6,}  100.0%")
    say(f"  {'smoking (no code = reference)':<34} {len(d):>6,}  100.0%")
    say("\n  TST is missing for 6 of 19,173, so every TST model is a complete case on 19,167.")
    say("  The BASELINE model never holds the exposure, so baseline C is identical for the T90")
    say("  and the TST run of the same rung. Only the with-exposure column differs.")

    tstmin = d[TST].describe()
    say(f"\n  TST_min raw   n {int(tstmin['count']):,}  mean {tstmin['mean']:.1f} min  "
        f"sd {tstmin['std']:.1f}  median {tstmin['50%']:.1f}  "
        f"IQR {d[TST].quantile(.25):.1f} to {d[TST].quantile(.75):.1f}")
    say("  The exposure enters as the within-site rank inverse normal of TST_min, so dC is "
        "per 1 SD of that.")

    # -------------------------------------------------------------- positive controls
    say("\n" + "-" * 100)
    say("2. POSITIVE CONTROLS, BOTH BEFORE ANY LADDER RUNS")
    say("-" * 100)

    pc_t = run_model("PC_TST", [], d, OUT, AGE4, "PC", "tst_z")
    t_dc, t_base, t_with = (float(pc_t.dC.mean()), float(pc_t.c_base.mean()),
                            float(pc_t.c_expo.mean()))
    ok_t = (abs(t_dc - PUB_TST_DC) < 1e-9 and abs(t_base - PUB_T90_BASE) < 5e-6
            and abs(t_with - PUB_TST_CWITH) < 1e-9)
    say(f"  TST, ranking_v3.csv rank {PUB_TST_RANK} of 197")
    say(f"    published dC          {PUB_TST_DC:.12f}   here {t_dc:.12f}   "
        f"diff {t_dc - PUB_TST_DC:+.3e}")
    say(f"    published C with TST  {PUB_TST_CWITH:.12f}   here {t_with:.12f}   "
        f"diff {t_with - PUB_TST_CWITH:+.3e}")
    say(f"    published baseline C  {PUB_T90_BASE:.6f}         here {t_base:.6f}   "
        f"diff {t_base - PUB_T90_BASE:+.3e}")
    say(f"    VERDICT {'PASS' if ok_t else 'FAIL'}")

    pc_o = run_model("PC_T90", [], d, OUT, AGE4, "PC", "t90_z")
    o_dc, o_base = float(pc_o.dC.mean()), float(pc_o.c_base.mean())
    ok_o = abs(o_dc - PUB_T90_DC) < 1e-9 and abs(o_base - PUB_T90_BASE) < 5e-6
    say(f"  T90")
    say(f"    published dC          {PUB_T90_DC:.12f}   here {o_dc:.12f}   "
        f"diff {o_dc - PUB_T90_DC:+.3e}")
    say(f"    published baseline C  {PUB_T90_BASE:.6f}         here {o_base:.6f}   "
        f"diff {o_base - PUB_T90_BASE:+.3e}")
    say(f"    VERDICT {'PASS' if ok_o else 'FAIL'}")

    if not (ok_t and ok_o):
        raise SystemExit("positive control failed, stop")

    # -------------------------------------------------------------- run everything
    per_all, rows = [], []

    def do(block, frame, spec_list, expo, agecols=AGE4,
           basis="cr(df=4), published 4-column", namesuffix=""):
        for item in spec_list:
            if len(item) == 5:
                name, extra, desc, oxy, cp = item
            else:
                name, extra = item[0], item[1]
                desc, oxy, cp = (item[2] if len(item) > 2 else name), False, ""
            nm = name + namesuffix
            per = run_model(nm, extra, frame, OUT, agecols, block, expo)
            per_all.append(per)
            rows.append(summarise_model(per, block, nm, desc, frame, extra, agecols, basis,
                                        expo, oxy, cp))
            r = rows[-1]
            say(f"  {block:<14}{expo:<7}{nm:<42} base C {r['mean_C_baseline']:.4f}  "
                f"+expo C {r['mean_C_with_exposure']:.4f}  dC {r['mean_dC']:+.5f}  "
                f"({r['n_outcomes']}/48)")

    say("\n" + "-" * 100)
    say("3. FITTING")
    say("-" * 100)
    say("  TST ladder, full cohort")
    do("full", d, LADDER_FULL, "tst_z")

    say("\n  T90 on the same full-cohort rungs. M1, M2 and D5 are re-runs that must equal "
        "results.csv.")
    say("  The __noTSTcov rows are new: the published rung with its TST covariate removed, so "
        "both exposures face one baseline.")
    do("full", d, [(n, e, f"T90 re-run of results.csv {n}", False, n) for n, e in T90_VERIFY_FULL],
       "t90_z", namesuffix="__RERUN")
    do("full", d, [(n, e, ds, "spo2mean_z" in e, n) for n, e, ds in T90_NEW_FULL], "t90_z")

    dB = rerank(d[d.bmi.notna()])
    say(f"\n  block B, BMI-complete subset n = {len(dB):,}   " +
        ", ".join(f"{k} n={v:,}" for k, v in dB.site_id.value_counts().sort_index().items()))
    do("bmi_6325", dB, LADDER_BMI, "tst_z")
    do("bmi_6325", dB,
       [(n, e, f"T90 re-run of results.csv {n}", False, n) for n, e in T90_VERIFY_BMI],
       "t90_z", namesuffix="__RERUN")
    do("bmi_6325", dB, [(n, e, ds, "spo2mean_z" in e, n) for n, e, ds in T90_NEW_BMI], "t90_z")

    say("\n  age-spline drop check, TST ladder, one cr column removed")
    do("full_agedrop", d, LADDER_FULL, "tst_z", agecols=AGE3,
       basis="cr(df=4), one column dropped", namesuffix="__agedrop")

    P = pd.concat(per_all, ignore_index=True)
    R = pd.DataFrame(rows)
    P.to_csv(os.path.join(HERE, "per_outcome_tst.csv"), index=False)
    R.to_csv(os.path.join(HERE, "results_tst.csv"), index=False)

    # -------------------------------------------------------------- verification
    say("\n" + "-" * 100)
    say("4. DOES THIS SCRIPT REPRODUCE results.csv FOR THE T90 ROWS IT RE-RAN?")
    say("-" * 100)
    old = pd.read_csv(os.path.join(HERE, "results.csv"))
    checks = []
    for blk, name in ([("full", n) for n, _ in T90_VERIFY_FULL] +
                      [("bmi_6325", n) for n, _ in T90_VERIFY_BMI]):
        o = old[(old.block == blk) & (old.model == name)].iloc[0]
        n_ = R[(R.block == blk) & (R.model == name + "__RERUN")].iloc[0]
        dd = float(n_.mean_dC) - float(o.mean_dC)
        db = float(n_.mean_C_baseline) - float(o.mean_C_baseline)
        checks.append({"block": blk, "model": name, "results_csv_dC": float(o.mean_dC),
                       "rerun_dC": float(n_.mean_dC), "diff_dC": dd,
                       "results_csv_baseline": float(o.mean_C_baseline),
                       "rerun_baseline": float(n_.mean_C_baseline), "diff_baseline": db,
                       "match": bool(abs(dd) < 1e-12 and abs(db) < 1e-12)})
        say(f"  {blk:<11}{name:<32} results.csv dC {float(o.mean_dC):+.6f}   "
            f"rerun {float(n_.mean_dC):+.6f}   diff {dd:+.2e}   "
            f"{'MATCH' if checks[-1]['match'] else 'DIFFERS'}")
    if not all(c["match"] for c in checks):
        raise SystemExit("re-run of the T90 rows did not match results.csv, stop")

    # -------------------------------------------------------------- headline
    say("\n" + "=" * 100)
    say("5. THE TST LADDER. MEAN dC ACROSS THE 48 OUTCOMES, HELD OUT AND TRANSPORTED "
        "BETWEEN HOSPITALS")
    say("=" * 100)
    say(f"{'rung':<32}{'baseline C':>11}{'+TST C':>9}{'dC':>10}{'rel gain':>10}"
        f"{'n out':>7}{'dC>0':>6}  oxygen  covariates")
    for blk in ["full", "bmi_6325"]:
        for _i, r in R[(R.block == blk) & (R.exposure == "tst_z")].iterrows():
            say(f"{r.model:<32}{r.mean_C_baseline:>11.4f}{r.mean_C_with_exposure:>9.4f}"
                f"{r.mean_dC:>+10.5f}{100*r.rel_gain_on_means:>9.2f}%{r.n_outcomes:>7d}"
                f"{r.n_outcomes_dC_positive:>6d}  {'YES' if r.oxygen_row else ' - ':<6}  "
                f"{r.covariates}")
        say("")
    say("  rel gain = mean dC / (1 - mean baseline C), the share of the still-unclaimed "
        "discrimination the exposure closes.")
    say("  Rows marked oxygen YES hold mean SpO2. Alen has said these are not the rows he "
        "cares about. Run for completeness.")
    say("  Block bmi_6325 carries its own matched age + sex control, so never compare a C "
        "across blocks.")

    say("\n" + "-" * 100)
    say("6. AGE-SPLINE DROP CHECK, TST (one cr column removed, the standing rule)")
    say("-" * 100)
    a = R[(R.block == "full") & (R.exposure == "tst_z")].set_index("model")
    b_ = R[R.block == "full_agedrop"].set_index("model")
    for name, _e, _d, _o, _c in LADDER_FULL:
        say(f"  {name:<32} 4-column dC {a.loc[name,'mean_dC']:+.5f}   "
            f"one dropped {b_.loc[name+'__agedrop','mean_dC']:+.5f}   "
            f"diff {b_.loc[name+'__agedrop','mean_dC']-a.loc[name,'mean_dC']:+.5f}")

    # -------------------------------------------------------------- organs
    say("\n" + "-" * 100)
    say("7. PER-ORGAN dC FOR THE L4 / L4b RUNG, T90 BESIDE TST")
    say("-" * 100)
    say("  L4 and L4b are the same covariate set: age + sex + race + smoking + AHI + sleep")
    say("  efficiency. L4b is defined as the T90 ladder's D5_clinical_sleep_noOxygen minus its")
    say("  TST covariate, and that removal lands exactly on L4.")

    def organ_table(pairs, title, index_col="organ5"):
        keys = {(b, m) for b, m, _l in pairs}
        inset = np.array([(bb, mm) in keys for bb, mm in zip(P.block, P.model)])
        sub = P[P.dC.notna().values & inset]
        piv = sub.pivot_table(index=index_col,
                              columns=["block", "model"], values="dC", aggfunc="mean")
        say(f"\n  {title}")
        order = ([i for i in ["cardiac", "respiratory", "metabolic", "kidney", "other"]
                  if i in piv.index] if index_col == "organ5" else list(piv.index))
        cnt = P[(P.block == pairs[0][0]) & (P.model == pairs[0][1]) & P.dC.notna()
                ].groupby(index_col).size()
        say(f"  {'organ':<14}{'n':>4}" + "".join(f"{lab:>26}" for _b, _m, lab in pairs))
        for o in order:
            say(f"  {o:<14}{int(cnt.get(o,0)):>4}" +
                "".join(f"{piv.loc[o,(b,m)]:>+26.5f}" for b, m, _l in pairs))

    L4PAIRS = [("full", "D5_clinical_sleep_noOxygen__RERUN", "T90 with TST covar"),
               ("full", "D5_clinical_sleep_noOxygen__noTSTcov", "T90 no TST covar"),
               ("full", "L4_clinical_AHI_SE", "TST (L4/L4b)")]
    organ_table(L4PAIRS, "L4 / L4b rung, five groups")
    organ_table(L4PAIRS, "L4 / L4b rung, the paper's full organ grouping", index_col="organ")
    organ_table([("full", "M1_paper_age_sex__RERUN", "T90 age+sex"),
                 ("full", "L1_age_sex", "TST age+sex")],
                "L1 age + sex rung, five groups")

    # -------------------------------------------------------------- top and bottom
    say("\n" + "-" * 100)
    say("8. WHERE TST HELPS AND WHERE IT HURTS, L1 AND L4")
    say("-" * 100)
    for m in ["L1_age_sex", "L4_clinical_AHI_SE"]:
        s = P[(P.block == "full") & (P.model == m) & P.dC.notna()].sort_values("dC",
                                                                              ascending=False)
        say(f"\n  [{m}] best 8")
        say(f"    {'outcome':<22}{'organ':<14}{'base C':>9}{'+TST C':>9}{'dC':>10}{'rel':>8}")
        for _i, r in s.head(8).iterrows():
            say(f"    {r.outcome:<22}{r.organ:<14}{r.c_base:>9.4f}{r.c_expo:>9.4f}"
                f"{r.dC:>+10.5f}{100*r.rel_gain:>7.2f}%")
        say(f"  [{m}] worst 8")
        for _i, r in s.tail(8).iterrows():
            say(f"    {r.outcome:<22}{r.organ:<14}{r.c_base:>9.4f}{r.c_expo:>9.4f}"
                f"{r.dC:>+10.5f}{100*r.rel_gain:>7.2f}%")

    # -------------------------------------------------------------- answer
    say("\n" + "=" * 100)
    say("9. THE ANSWER")
    say("=" * 100)
    ft = R[(R.block == "full") & (R.exposure == "tst_z")].set_index("model")
    fo = R[(R.block == "full") & (R.exposure == "t90_z")].set_index("model")
    bt = R[(R.block == "bmi_6325") & (R.exposure == "tst_z")].set_index("model")
    say(f"  Against age + sex, TST adds dC {ft.loc['L1_age_sex','mean_dC']:+.5f} on a baseline "
        f"C of {ft.loc['L1_age_sex','mean_C_baseline']:.4f}. It is a LOSS of discrimination, "
        f"and it reproduces")
    say(f"  the published rank {PUB_TST_RANK} of 197 exactly. T90 on the same baseline adds "
        f"{fo.loc['M1_paper_age_sex__RERUN','mean_dC']:+.5f}.")
    say(f"  Against age + sex + AHI, TST adds {ft.loc['L2_age_sex_AHI','mean_dC']:+.5f}.")
    say(f"  Against age + sex + race + smoking, TST adds "
        f"{ft.loc['L3_age_sex_race_smoke','mean_dC']:+.5f}.")
    say(f"  Against the L4 clinical line, TST adds {ft.loc['L4_clinical_AHI_SE','mean_dC']:+.5f} "
        f"while T90 on the IDENTICAL baseline adds "
        f"{fo.loc['D5_clinical_sleep_noOxygen__noTSTcov','mean_dC']:+.5f}.")
    say(f"  On the BMI-complete rows the matched age + sex control gives TST "
        f"{bt.loc['B1_age_sex','mean_dC']:+.5f}, adding BMI gives "
        f"{bt.loc['B2_age_sex_BMI','mean_dC']:+.5f}, the full clinical set with BMI gives "
        f"{bt.loc['B3_clinical_BMI','mean_dC']:+.5f}.")

    json.dump({"positive_control": {
                    "tst_dC": t_dc, "tst_published_dC": PUB_TST_DC,
                    "tst_C_with": t_with, "tst_published_C_with": PUB_TST_CWITH,
                    "tst_rank_in_ranking_v3": PUB_TST_RANK,
                    "t90_dC": o_dc, "t90_published_dC": PUB_T90_DC,
                    "baseline_C": t_base, "published_baseline_C": PUB_T90_BASE,
                    "tst_pass": bool(ok_t), "t90_pass": bool(ok_o)},
               "results_csv_rerun_check": checks,
               "model_spec": {"penalizer": PEN, "n_outcomes": N_RANKED_OUTCOMES,
                              "folds": "fit I0002 score I0006, then reverse, nanmean",
                              "min_events_per_fold": 30,
                              "age_basis": "cr(df=4) 4 columns, plus a one-column-dropped pass",
                              "exposure_transform": "within-site rank inverse normal, per 1 SD"},
               "availability": avail},
              open(os.path.join(HERE, "provenance_tst.json"), "w"), indent=2)

    with open(os.path.join(HERE, "REPORT_TST.txt"), "w") as f:
        f.write("\n".join(LOG) + "\n")
    say("\nwrote results_tst.csv, per_outcome_tst.csv, provenance_tst.json, REPORT_TST.txt")


if __name__ == "__main__":
    main()
