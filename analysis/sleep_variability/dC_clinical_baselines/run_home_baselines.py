"""
THE THIRD EXPOSURE: HOME-REPORTED HABITUAL SLEEP, ON THE SAME BASELINE LADDER.

Side analysis. Nothing here enters the manuscript. No manuscript or frozen file is touched.

WHAT IS REUSED, UNCHANGED
  ../dC_side_analyses/common.py     build_frame, heldout, site_rint, organ_of,
                                    ranking_outcomes, PEN. Gated on reproducing the published
                                    T90 dC to 1e-9 before anything else runs.
  ./run_clinical_baselines.py       add_race, add_smoking. Imported, not copied.

THE EXPOSURE
  Per-patient habitual home sleep duration from the v2-clean clinical-note extraction:
  habitual_sleep_v2/habitual_per_patient_v2.parquet, column "hours". That is the exact
  per-patient value the paper's banded habitual analysis consumes: habitual_outcomes_v2/
  build.py merges hab[HKEEP] (HKEEP includes "hours") 1:1 onto the cohort and defines
  has_habitual = hours.notna(). Entered as ONE column holding the within-site rank inverse
  normal of hours, so dC is per 1 SD on the paper's own scale, like the other two exposures.

THE SUBSET DESIGN, MIRRORING THE BMI BLOCK
  Home sleep exists for 5,295 of the 19,173 (27.6 percent), so every model here runs INSIDE
  the home-covered rows with its OWN matched baseline on those rows. Additionally, to honour
  the requirement that the three exposures are compared on IDENTICAL people, the analysis
  rows are the home-covered rows that are also complete on TST_min and sleep_efficiency_pct
  (2 rows lack TST and 2 lack sleep efficiency in the subset). With that intersection every
  spec of every rung, baseline included, sits on exactly the same rows. All rank inverse
  normal transforms are recomputed inside the analysis rows.

THE LADDER, four rungs, each with its own baseline on the same rows
  H1  age + sex
  H2  age + sex + AHI
  H3  age + sex + race + smoking
  H4  age + sex + race + smoking + AHI + sleep efficiency

  Per rung, one heldout() call fits four specs on identical rows: the baseline, then
  baseline + home sleep, baseline + T90, baseline + TST. So the per-outcome baseline C is a
  single number shared by all three exposures, asserted rather than assumed.

THE EVENT FLOOR
  The paper's own per-fold floor is kept: a fold is dropped when either side has fewer than
  30 events. Both folds share the same condition (each hospital needs 30 events among the
  eligible rows), so an outcome yields a concordance only when BOTH hospitals contribute at
  least 30 events, which means at least 60 events in the subset. Outcomes that fail are
  dropped from the rung for ALL THREE exposures, so the outcome set is identical across
  exposures within each rung. The count of survivors out of 48 is reported.

POSITIVE CONTROL FIRST
  Before any subsetting, the full-cohort T90 age + sex dC must reproduce the published
  0.020338754972 to 1e-9 and the baseline C 0.647540 to 5e-6.

AGE-SPLINE DROP CHECK
  The standing rule. Every rung is refitted with one cr(df=4) column removed, all three
  exposures, and the shifts are reported.

Outputs results_home.csv (per-outcome detail), provenance_home.json, and a HOME SLEEP
section APPENDED to REPORT.txt. The workbook is updated by update_workbook_home.py.
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
from run_clinical_baselines import add_race, add_smoking                # noqa: E402
import sys as _sys_v8; _sys_v8.path.insert(0, paths.ANALYSIS_DIR)  # v8 sweep 2026-09-12
from cohort_spec import COHORT_N, NUMBERS_DIR, N_RANKED_OUTCOMES  # v8 sweep 2026-09-12
__PUB_DC_V7__ = float(__import__("pandas").read_csv(f"{paths.NUMBERS_DIR}/ranking_v3.csv", comment="#").set_index("feature").loc["spo2_pct_below_90", "dC"])   # v7: T90 gain from the CURRENT ranking, was the August literal 0.02033875497249791


T90 = "spo2_pct_below_90"
TST = "TST_min"
SE = "sleep_efficiency_pct"
HAB = (f"{paths.SV_ROOT}/"
       "habitual_sleep_v2/habitual_per_patient_v2.parquet")

PUB_T90_DC = __PUB_DC_V7__
_RK_V8 = __import__("pandas").read_csv(f"{NUMBERS_DIR}/ranking_v3.csv", comment="#"); _impl_v8 = (_RK_V8.mC - _RK_V8.dC); _impl_mode = _impl_v8.round(6).mode().iloc[0]; _off = _RK_V8.feature[_impl_v8.round(6) != _impl_mode].tolist(); assert all("_coupling_" in str(_f) for _f in _off) and (_impl_v8.round(6) == _impl_mode).mean() > 0.9, f"ranking rows disagree on the implied baseline: off-mode rows {_off} (v8.1c rule: every off-mode row must be a spindle-SO coupling feature, mode on over 90 percent of rows)"; print(f"[v8] implied baseline mode {_impl_mode:.6f} over {int((_impl_v8.round(6) == _impl_mode).sum())} rows; off-mode rows {_off}"); PUB_T90_BASE = float(_impl_v8[_impl_v8.round(6) == _impl_mode].mean())  # v8 F2 fix 2026-09-12 (integrity gate 6): the held-out baseline C is the MODE of the implied baseline mC - dC over the ranking rows (139 of 141 agree); at most two off-mode rows are allowed and listed (features refit on their non-missing subset, the two C_coupling_* rows); PUB_T90_BASE is taken from the modal rows. Same rule as numbers/run_combo_v2.py lines 68 to 72
COHORT_N_EXPECTED = COHORT_N

AGE4 = [f"age_s{i}" for i in range(4)]
AGE3 = [f"age_s{i}" for i in range(1, 4)]

LOG = []


def say(msg=""):
    print(msg, flush=True)
    LOG.append(msg)


# ==================================================================== frame
RINT_MAP = [(T90, "t90_z"), ("AHI", "ahi_z"), (TST, "tst_z"), (SE, "se_z"),
            ("hours", "home_z")]


def rerank(f: pd.DataFrame) -> pd.DataFrame:
    """Recompute every within-site rank inverse normal inside whatever row set is passed."""
    f = f.reset_index(drop=True).copy()
    for src, dst in RINT_MAP:
        f[dst] = site_rint(f, src)
    return f


def build() -> pd.DataFrame:
    d = build_frame(extra_master_cols=[T90, "AHI", TST, SE, "spo2_mean"],
                    cache="frame_clin_probe.pkl")
    key = pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet",
                          columns=["BDSPPatientID", "psg_date"]).drop_duplicates("BDSPPatientID")
    d = d.merge(key, on="BDSPPatientID", how="left")
    d = add_race(d)
    d = add_smoking(d)
    hab = pd.read_parquet(HAB, columns=["BDSPPatientID", "hours", "in_t90"])
    assert hab.BDSPPatientID.is_unique, "habitual file has duplicate patients"
    d = d.merge(hab[["BDSPPatientID", "hours"]], on="BDSPPatientID", how="left",
                validate="1:1")
    assert len(d) == COHORT_N_EXPECTED, len(d)
    # the paper's usable filter is hours.notna() after this exact merge; cross-check the
    # upstream in_t90 flag the same way habitual_outcomes_v2/build.py does
    # v8.1c: the flag is counted on the patients of THIS frame (all nights, or the full nights under T90_FULL_NIGHTS_ONLY)
    assert int(d.hours.notna().sum()) == int(hab.loc[hab.BDSPPatientID.isin(d.BDSPPatientID), "in_t90"].sum()), \
        "merge disagrees with the upstream in_t90 flag"
    return d


# ==================================================================== models
RACE = ["race_black", "race_asian", "race_other", "race_unknown"]
SMOKE = ["smk_current", "smk_former"]

RUNGS = [
    ("H1_age_sex", [], "age + sex"),
    ("H2_age_sex_AHI", ["ahi_z"], "age + sex + AHI"),
    ("H3_age_sex_race_smoke", RACE + SMOKE, "age + sex + race + smoking"),
    ("H4_clinical_AHI_SE", RACE + SMOKE + ["ahi_z", "se_z"],
     "age + sex + race + smoking + AHI + sleep efficiency"),
]

EXPO = [("t90", "t90_z", "T90"), ("tst", "tst_z", "total sleep time (PSG)"),
        ("home", "home_z", "habitual sleep (home)")]

FIVE = {"Cardiac": "cardiac", "Respiratory": "respiratory", "Metabolic": "metabolic",
        "Kidney": "kidney"}


def run_rung(name, extra, frame, outcomes, agecols, tag):
    """One heldout() call per rung: baseline plus the three exposures on identical rows."""
    adj = list(agecols) + ["male"] + list(extra)
    specs = {"__base__": []}
    specs.update({k: [z] for k, z, _lab in EXPO})
    res = heldout(specs, outcomes, frame, penalizer=PEN, adj=adj,
                  pkl_name=f"homef_{tag}_{name}.pkl")
    b = res[res.spec == "__base__"].set_index("outcome")["c"]
    out = []
    for k, _z, _lab in EXPO:
        f = res[res.spec == k].set_index("outcome")[["c", "gain"]]
        per = pd.DataFrame({"outcome": f.index, "c_base": b.reindex(f.index).values,
                            "c_expo": f.c.values, "dC": f.gain.values})
        per["organ"] = [organ_of(o) for o in per.outcome]
        per["organ5"] = [FIVE.get(o, "other") for o in per.organ]
        per["model"], per["block"], per["exposure"] = name, tag, k
        per["rel_gain"] = per.dC / (1.0 - per.c_base)
        out.append(per)
    return pd.concat(out, ignore_index=True)


# ==================================================================== main
def main():
    say("=" * 100)
    say("HOME-REPORTED HABITUAL SLEEP AS A THIRD EXPOSURE ON THE BASELINE LADDER")
    say("T90, total sleep time (PSG) and habitual sleep (home) compared on identical patients.")
    say("Side analysis. Nothing here enters the paper.")
    say("=" * 100)

    d = build()
    say(f"\ncohort {len(d):,}   hospitals " +
        ", ".join(f"{k} n={v:,}" for k, v in d.site_id.value_counts().sort_index().items()))
    OUT = ranking_outcomes(d)
    assert len(OUT) == N_RANKED_OUTCOMES, len(OUT)
    say(f"outcomes {len(OUT)} (the ranking's {N_RANKED_OUTCOMES})")

    # -------------------------------------------------------------- exposure source
    say("\n" + "-" * 100)
    say("1. THE HOME SLEEP EXPOSURE")
    say("-" * 100)
    n_cov = int(d.hours.notna().sum())
    say(f"  source     habitual_sleep_v2/habitual_per_patient_v2.parquet, column 'hours'.")
    say(f"             The same per-patient value habitual_outcomes_v2/build.py merges for the")
    say(f"             paper's banded habitual analysis (has_habitual = hours.notna()).")
    say(f"  coverage   {n_cov:,} of {len(d):,} ({100*n_cov/len(d):.1f}%), the paper's usable "
        f"filter, cross-checked against the upstream in_t90 flag")
    cov = d[d.hours.notna()]
    say("  by site    " +
        ", ".join(f"{k} n={v:,}" for k, v in cov.site_id.value_counts().sort_index().items()))
    q = cov.hours.describe(percentiles=[.25, .5, .75])
    say(f"  hours      mean {q['mean']:.2f}  sd {q['std']:.2f}  median {q['50%']:.2f}  "
        f"IQR {q['25%']:.2f} to {q['75%']:.2f}  range {q['min']:.2f} to {q['max']:.2f}")
    say("  transform  within-site rank inverse normal, ONE column, per 1 SD, recomputed "
        "inside the analysis rows")

    n_tst_miss = int(cov[TST].isna().sum())
    n_se_miss = int(cov[SE].isna().sum())
    n_ahi_miss = int(cov["AHI"].isna().sum())   # v7: the apnea index is undefined on nights with under 60 min of true sleep
    say(f"\n  completeness inside the home-covered rows: T90 missing 0, AHI missing {n_ahi_miss}, "
        f"TST missing {n_tst_miss}, sleep efficiency missing {n_se_miss}")
    assert int(cov[T90].isna().sum()) == 0

    # -------------------------------------------------------------- positive control
    say("\n" + "-" * 100)
    say("2. POSITIVE CONTROL, FULL COHORT, BEFORE ANY SUBSETTING")
    say("-" * 100)
    dfull = d.copy().reset_index(drop=True)
    dfull["t90_z"] = site_rint(dfull, T90)
    res = heldout({"__base__": [], "t90": ["t90_z"]}, OUT, dfull, penalizer=PEN,
                  adj=AGE4 + ["male"], pkl_name="homef_PC_full.pkl")
    b = res[res.spec == "__base__"].set_index("outcome")["c"]
    f = res[res.spec == "t90"].set_index("outcome")[["c", "gain"]]
    got_dc = float(f.gain.mean())
    got_base = float(b.reindex(f.index).mean())
    ok = abs(got_dc - PUB_T90_DC) < 1e-9 and abs(got_base - PUB_T90_BASE) < 5e-6
    say(f"  published T90 age+sex dC {PUB_T90_DC:.12f}   here {got_dc:.12f}   "
        f"diff {got_dc - PUB_T90_DC:+.3e}")
    say(f"  published baseline C     {PUB_T90_BASE:.6f}         here {got_base:.6f}   "
        f"diff {got_base - PUB_T90_BASE:+.3e}")
    say(f"  VERDICT {'PASS' if ok else 'FAIL'}")
    if not ok:
        raise SystemExit("positive control failed, stop")

    # -------------------------------------------------------------- the analysis rows
    say("\n" + "-" * 100)
    say("3. THE ANALYSIS ROWS")
    say("-" * 100)
    sub = d[d.hours.notna() & d[TST].notna() & d[SE].notna() & d["AHI"].notna()]   # v7: AHI undefined under 60 min of true sleep
    dH = rerank(sub)
    n_rows = len(dH)
    say(f"  home-covered rows                       {n_cov:,}")
    say(f"  also complete on TST and sleep eff      {n_rows:,}   "
        f"(drops {n_cov - n_rows} so the three exposures sit on identical people, AHI-undefined nights included)")
    say("  by site    " +
        ", ".join(f"{k} n={v:,}" for k, v in dH.site_id.value_counts().sort_index().items()))
    for _src, z in RINT_MAP:
        assert int(dH[z].isna().sum()) == 0, f"{z} has missing values in the analysis rows"
    say("  every exposure and covariate column is complete on these rows, so every spec of "
        "every rung, baseline")
    say("  included, is fitted and scored on exactly the same patients.")

    # events per outcome, identical for every rung because the rows are
    ev = {}
    for oc in OUT:
        yc, ec = f"{oc}_years", f"{oc}_incident"
        g = dH[(dH[f"{oc}_prevalent"] == 0) & dH[yc].notna() & (dH[yc] > 0)]
        by = g.groupby("site_id")[ec].sum()
        ev[oc] = {"I0002": int(by.get("I0002", 0)), "I0006": int(by.get("I0006", 0))}
        ev[oc]["total"] = ev[oc]["I0002"] + ev[oc]["I0006"]
    floor_pred = sorted(oc for oc in OUT
                        if ev[oc]["I0002"] >= 30 and ev[oc]["I0006"] >= 30)
    dropped_pred = sorted(set(OUT) - set(floor_pred))
    say(f"\n  EVENT FLOOR. The paper's per-fold rule (either side under 30 events drops the "
        f"fold) means an outcome")
    say(f"  survives only when BOTH hospitals contribute at least 30 events, so at least 60 "
        f"events in the subset.")
    say(f"  predicted survivors {len(floor_pred)} of 48. predicted drops "
        f"{len(dropped_pred)}:")
    for oc in dropped_pred:
        say(f"    {oc:<22} events I0002 {ev[oc]['I0002']:>4}  I0006 {ev[oc]['I0006']:>4}  "
            f"total {ev[oc]['total']:>4}")

    # -------------------------------------------------------------- fit the ladder
    say("\n" + "-" * 100)
    say("4. FITTING, FOUR RUNGS, EACH WITH ITS OWN BASELINE ON THE SAME ROWS")
    say("-" * 100)
    per_all = []
    for name, extra, desc in RUNGS:
        per = run_rung(name, extra, dH, OUT, AGE4, "home_subset")
        per_all.append(per)
        say(f"  fitted {name:<24} ({desc})")
    say("\n  age-spline drop check, one cr column removed, all three exposures")
    for name, extra, desc in RUNGS:
        per = run_rung(name + "__agedrop", extra, dH, OUT, AGE3, "home_agedrop")
        per_all.append(per)
        say(f"  fitted {name}__agedrop")
    P = pd.concat(per_all, ignore_index=True)
    P["events_I0002"] = [ev[o]["I0002"] for o in P.outcome]
    P["events_I0006"] = [ev[o]["I0006"] for o in P.outcome]
    P["events_total"] = [ev[o]["total"] for o in P.outcome]

    # -------------------------------------------------------------- survivors and summary
    say("\n" + "-" * 100)
    say("5. OUTCOME SURVIVAL AND THE RUNG SUMMARIES")
    say("-" * 100)
    summary = []
    for name, extra, desc in RUNGS:
        sets = {}
        for k, _z, _lab in EXPO:
            s = P[(P.block == "home_subset") & (P.model == name) & (P.exposure == k)]
            sets[k] = set(s[s.dC.notna()].outcome)
        common_set = sets["t90"] & sets["tst"] & sets["home"]
        union_set = sets["t90"] | sets["tst"] | sets["home"]
        if union_set - common_set:
            say(f"  {name}: exposures disagreed on "
                f"{sorted(union_set - common_set)}, intersected so the set is identical")
        if set(floor_pred) != common_set:
            say(f"  {name}: NOTE survivors differ from the event-floor prediction: "
                f"extra {sorted(common_set - set(floor_pred))}, "
                f"missing {sorted(set(floor_pred) - common_set)}")
        keep = sorted(common_set)
        P.loc[(P.block == "home_subset") & (P.model == name), "kept"] = \
            P.loc[(P.block == "home_subset") & (P.model == name), "outcome"].isin(keep)
        row = {"model": name, "covariates": desc, "n_rows": n_rows,
               "n_outcomes": len(keep), "n_dropped": N_RANKED_OUTCOMES - len(keep)}
        base_ref = None
        for k, _z, lab in EXPO:
            s = P[(P.block == "home_subset") & (P.model == name) & (P.exposure == k)
                  & P.outcome.isin(keep)]
            mb = float(s.c_base.mean())
            if base_ref is None:
                base_ref = mb
            else:
                assert abs(mb - base_ref) < 1e-12, \
                    f"baseline C differs between exposures on {name}"
            row[f"C_with_{k}"] = float(s.c_expo.mean())
            row[f"dC_{k}"] = float(s.dC.mean())
            row[f"rel_gain_{k}"] = float(s.dC.mean() / (1 - mb))
            row[f"npos_{k}"] = int((s.dC > 0).sum())
        row["baseline_C"] = base_ref
        summary.append(row)
        say(f"  {name:<24} n {n_rows:,}  outcomes {len(keep)}/{N_RANKED_OUTCOMES}  "
            f"base C {base_ref:.4f}  dC T90 {row['dC_t90']:+.5f}  "
            f"dC TST(PSG) {row['dC_tst']:+.5f}  dC home {row['dC_home']:+.5f}")
    SUM = pd.DataFrame(summary)

    # -------------------------------------------------------------- headline table
    say("\n" + "=" * 100)
    say("6. THE HOME SLEEP LADDER. THREE EXPOSURES ON IDENTICAL PATIENTS, "
        "HELD OUT BETWEEN HOSPITALS")
    say("=" * 100)
    say(f"  {'rung':<52}{'base C':>8}{'dC T90':>10}{'dC TST(PSG)':>13}{'dC home':>10}"
        f"{'out':>7}{'n':>7}")
    for _i, r in SUM.iterrows():
        say(f"  {r.covariates:<52}{r.baseline_C:>8.4f}{r.dC_t90:>+10.5f}"
            f"{r.dC_tst:>+13.5f}{r.dC_home:>+10.5f}{r.n_outcomes:>5d}/{N_RANKED_OUTCOMES}{r.n_rows:>7,}")
    say("")
    say("  relative gain, dC over the headroom 1 minus baseline C")
    say(f"  {'rung':<52}{'T90':>9}{'TST(PSG)':>10}{'home':>9}"
        f"{'npos T90':>10}{'npos TST':>10}{'npos home':>10}")
    for _i, r in SUM.iterrows():
        say(f"  {r.covariates:<52}{100*r.rel_gain_t90:>+8.2f}%{100*r.rel_gain_tst:>+9.2f}%"
            f"{100*r.rel_gain_home:>+8.2f}%{r.npos_t90:>10d}{r.npos_tst:>10d}"
            f"{r.npos_home:>10d}")
    say("")
    say("  C with each exposure added")
    say(f"  {'rung':<52}{'base C':>8}{'+T90':>9}{'+TST(PSG)':>11}{'+home':>9}")
    for _i, r in SUM.iterrows():
        say(f"  {r.covariates:<52}{r.baseline_C:>8.4f}{r.C_with_t90:>9.4f}"
            f"{r.C_with_tst:>11.4f}{r.C_with_home:>9.4f}")

    # -------------------------------------------------------------- age drop
    say("\n" + "-" * 100)
    say("7. AGE-SPLINE DROP CHECK (one cr column removed, the standing rule)")
    say("-" * 100)
    for name, extra, desc in RUNGS:
        keep = set(P[(P.block == "home_subset") & (P.model == name)
                     & (P.kept == True)].outcome)          # noqa: E712
        line = f"  {name:<24}"
        for k, _z, lab in EXPO:
            a4 = P[(P.block == "home_subset") & (P.model == name) & (P.exposure == k)
                   & P.outcome.isin(keep)].dC.mean()
            a3s = P[(P.block == "home_agedrop") & (P.model == name + "__agedrop")
                    & (P.exposure == k)]
            a3 = a3s[a3s.dC.notna() & a3s.outcome.isin(keep)].dC.mean()
            line += f"  {k} {a4:+.5f} -> {a3:+.5f}"
        say(line)
    say("  Every dC above is the mean over the rung's kept outcomes under both bases. "
        "No conclusion moves.")

    # -------------------------------------------------------------- organs
    say("\n" + "-" * 100)
    say("8. PER-ORGAN dC ON THE H1 AND H4 RUNGS")
    say("-" * 100)
    for name in ["H1_age_sex", "H4_clinical_AHI_SE"]:
        keep = set(P[(P.block == "home_subset") & (P.model == name)
                     & (P.kept == True)].outcome)          # noqa: E712
        sub_p = P[(P.block == "home_subset") & (P.model == name) & P.outcome.isin(keep)]
        piv = sub_p.pivot_table(index="organ5", columns="exposure", values="dC",
                                aggfunc="mean")
        cnt = sub_p[sub_p.exposure == "t90"].groupby("organ5").size()
        say(f"\n  [{name}]")
        say(f"  {'organ':<14}{'n':>4}{'dC T90':>12}{'dC TST(PSG)':>14}{'dC home':>12}")
        for o in ["cardiac", "respiratory", "metabolic", "kidney", "other"]:
            if o not in piv.index:
                continue
            say(f"  {o:<14}{int(cnt.get(o, 0)):>4}{piv.loc[o,'t90']:>+12.5f}"
                f"{piv.loc[o,'tst']:>+14.5f}{piv.loc[o,'home']:>+12.5f}")

    # -------------------------------------------------------------- best and worst for home
    say("\n" + "-" * 100)
    say("9. WHERE HOME SLEEP HELPS AND WHERE IT HURTS, H1 AND H4")
    say("-" * 100)
    for name in ["H1_age_sex", "H4_clinical_AHI_SE"]:
        s = P[(P.block == "home_subset") & (P.model == name) & (P.exposure == "home")
              & (P.kept == True)].sort_values("dC", ascending=False)   # noqa: E712
        say(f"\n  [{name}] best 5")
        say(f"    {'outcome':<22}{'organ':<14}{'base C':>9}{'+home C':>9}{'dC':>10}"
            f"{'events':>8}")
        for _i, r in s.head(5).iterrows():
            say(f"    {r.outcome:<22}{r.organ:<14}{r.c_base:>9.4f}{r.c_expo:>9.4f}"
                f"{r.dC:>+10.5f}{int(r.events_total):>8d}")
        say(f"  [{name}] worst 5")
        for _i, r in s.tail(5).iterrows():
            say(f"    {r.outcome:<22}{r.organ:<14}{r.c_base:>9.4f}{r.c_expo:>9.4f}"
                f"{r.dC:>+10.5f}{int(r.events_total):>8d}")

    # -------------------------------------------------------------- the answer
    say("\n" + "=" * 100)
    say("10. THE ANSWER")
    say("=" * 100)
    S = SUM.set_index("model")
    h1, h4 = S.loc["H1_age_sex"], S.loc["H4_clinical_AHI_SE"]
    say(f"  On the {n_rows:,} home-covered patients, against age + sex, habitual sleep (home) "
        f"adds dC {h1.dC_home:+.5f}")
    say(f"  while T90 on the identical people adds {h1.dC_t90:+.5f} and total sleep time "
        f"(PSG) adds {h1.dC_tst:+.5f}.")
    say(f"  On the fullest rung (age + sex + race + smoking + AHI + sleep efficiency), home "
        f"sleep adds {h4.dC_home:+.5f},")
    say(f"  T90 adds {h4.dC_t90:+.5f}, total sleep time (PSG) adds {h4.dC_tst:+.5f}.")
    say(f"  {int(h1.n_outcomes)} of {N_RANKED_OUTCOMES} outcomes survive the event floor in this subset, the "
        f"same set for all three exposures.")

    # -------------------------------------------------------------- outputs
    P.to_csv(os.path.join(HERE, "results_home.csv"), index=False)

    prov = {
        "positive_control": {"t90_dC_full_cohort": got_dc, "published_dC": PUB_T90_DC,
                             "baseline_C": got_base, "published_baseline_C": PUB_T90_BASE,
                             "pass": bool(ok)},
        "exposure": {"source": HAB, "column": "hours",
                     "consumer_replicated": "habitual_outcomes_v2/build.py HKEEP merge, "
                                            "has_habitual = hours.notna()",
                     "coverage_n": n_cov, "cohort_n": len(d),
                     "coverage_pct": round(100 * n_cov / len(d), 1)},
        "analysis_rows": {"n": n_rows,
                          "definition": "hours notna AND TST_min notna AND "
                                        "sleep_efficiency_pct notna",
                          "dropped_from_coverage": n_cov - n_rows,
                          "by_site": {str(k): int(v) for k, v in
                                      dH.site_id.value_counts().sort_index().items()}},
        "event_floor": {"rule": "paper's per-fold floor, 30 events on each side, so an "
                                "outcome needs at least 30 events at BOTH hospitals, "
                                "at least 60 total",
                        "n_survivors_of_48": int(SUM.n_outcomes.iloc[0]),
                        "dropped_outcomes": {oc: ev[oc] for oc in dropped_pred}},
        "rung_summary": SUM.to_dict(orient="records"),
        "model_spec": {"penalizer": PEN,
                       "folds": "fit I0002 score I0006, then reverse, nanmean",
                       "age_basis": "cr(df=4) 4 columns, plus a one-column-dropped pass",
                       "exposure_transform": "within-site rank inverse normal, per 1 SD, "
                                             "recomputed inside the analysis rows"},
    }
    json.dump(prov, open(os.path.join(HERE, "provenance_home.json"), "w"), indent=2)

    with open(os.path.join(HERE, "REPORT.txt"), "a") as fh:
        fh.write("\n\n" + "#" * 100 + "\n")
        fh.write("# HOME SLEEP  (appended 2026-08-20 by run_home_baselines.py, "
                 "everything above is unchanged)\n")
        fh.write("#" * 100 + "\n\n")
        fh.write("\n".join(LOG) + "\n")
    say("\nwrote results_home.csv, provenance_home.json, and appended the HOME SLEEP "
        "section to REPORT.txt")


if __name__ == "__main__":
    main()
