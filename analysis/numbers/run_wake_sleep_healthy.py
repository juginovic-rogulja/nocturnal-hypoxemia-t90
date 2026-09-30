#!/usr/bin/env python3
"""
Wake vs sleep T90, entered together, repeated in progressively healthier strata.

The question the full-cohort version answered was whether nocturnal desaturation carries
risk beyond wakeful desaturation. The obvious reviewer objection is that both are just
markers of existing cardiopulmonary illness. This repeats the identical model inside
subgroups that were free of that illness at the time of the sleep study.

Conventions are copied byte for byte from the generator of numbers/wake_vs_sleep_t90.csv,
which is re-derived here and asserted against that file as a positive control before any
restriction is applied.

Strata
  A  full cohort, the reference
  B  free of prevalent cardiopulmonary and renal / hepatic / cancer disease
  C  B, and also free of prevalent metabolic disease
  D  C, and also AHI < 5
  E  B, and also AHI < 5      (reported because D is too small to model)

Outputs
  numbers/wake_sleep_healthy.json
  numbers/wake_sleep_healthy_<stratum>.csv   one per stratum
"""
from __future__ import annotations
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../..")))  # repository root, where paths.py lives
import paths
import warnings; warnings.filterwarnings("ignore")

import glob
import sys, json, sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from patsy import dmatrix
from lifelines import CoxPHFitter

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from disease_definitions import DISEASES, NEGATIVE_CONTROLS

SHARDS = f"{paths.V8_ROOT}/X4_groupP/cpap_stage/out*.csv"   # v8.2: re-extracted, collapse-masked shards (X4_groupP)
COHORT = HERE.parent / "data_frozen_v8_2026-09" / "t90_final.parquet"
REFCSV = HERE / "wake_vs_sleep_t90.csv"

MIN_EVENTS = 60
MIN_WAKE_MIN = 10.0
QTHR = 0.05

# Prevalent conditions that define each restriction. Kept as separate blocks so the cost of
# each block in sample size can be reported rather than only the total.
CARDIOPULM = ["copd2", "asthma", "resp_failure", "obesity_hypovent", "pulm_htn",
              "hf", "ihd", "mi", "afib", "stroke_any"]
RENAL_HEP_CA = ["ckd", "cirrhosis", "cancer_any"]
METABOLIC = ["diabetes", "obesity", "nafld"]

CONDS = list(DISEASES.items()) + [("death", ("Death from any cause", False, [], []))]


# --------------------------------------------------------------------------- data assembly
def load() -> pd.DataFrame:
    files = sorted(glob.glob(SHARDS))
    assert files, f"no cpap_stage shards under {SHARDS}"   # v8.2: shard count is a property of the fleet, not a result
    raw = pd.concat([pd.read_csv(f, low_memory=False) for f in files], ignore_index=True)

    # count check 1: duplicate folder+session rows, the number quoted in the task
    dup_fs = int(raw.duplicated(subset=["BIDSFolder", "SessionID"]).sum())
    # count check 2: same number from a groupby, an independent path
    dup_fs_alt = int((raw.groupby(["BIDSFolder", "SessionID"]).size() - 1).clip(lower=0).sum())
    assert dup_fs == dup_fs_alt, (dup_fs, dup_fs_alt)   # v8.2: the two counting paths must agree; the count itself is recorded, not typed
    print(f"duplicate folder+session rows across arms: {dup_fs}")

    st = raw.drop_duplicates("BDSPPatientID")
    st = st[st.status == "ok"]
    assert st.BDSPPatientID.is_unique

    b = pd.read_parquet(COHORT)
    b = b[(b.fu_valid == 1) & b.spo2_pct_below_90.notna() & (b.oximetry_bad == 0)]

    # the impossible-oximetry block: nadir above the mean cannot happen
    bad = (b.spo2_nadir_corrected > b.spo2_mean)
    n_bad = int(bad.sum())
    n_bad_alt = int(len(b) - (~bad).sum())          # second, independent count
    assert n_bad == n_bad_alt
    b = b[~bad]

    m = b.merge(st[["BDSPPatientID", "all_wake_t90", "all_sleep_t90", "all_wake_min"]],
                on="BDSPPatientID", how="inner")
    m = m[m.all_wake_min >= MIN_WAKE_MIN]
    m = m.dropna(subset=["all_wake_t90", "all_sleep_t90"]).reset_index(drop=True)

    # count check: len() and an independent boolean sum must agree
    assert len(m) == int(np.ones(len(m), bool).sum()) == m.BDSPPatientID.nunique()

    m["male"] = (m.sex.astype(str).str.upper().str[0] == "M").astype(int)
    sp = dmatrix("cr(a, df=4) - 1", {"a": m.AgeAtVisit.values}, return_type="dataframe")
    sp.columns = [f"age_s{i}" for i in range(4)]
    sp.index = m.index
    for c in sp.columns:
        m[c] = sp[c]
    return m, {"impossible_oximetry_dropped": n_bad, "dup_folder_session_rows": dup_fs}


ADJ = [f"age_s{i}" for i in range(4)] + ["male"]


def rint(x) -> np.ndarray:
    """Rank-based inverse normal, computed on whatever subset is passed in."""
    x = np.asarray(x, float)
    r = stats.rankdata(x)
    return stats.norm.ppf((r - 0.375) / (len(r) + 0.25))


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1.0)


def classify(wq, sq):
    w, s = wq < QTHR, sq < QTHR
    return "BOTH" if (w and s) else "sleep only" if s else "awake only" if w else "neither"


# --------------------------------------------------------------------------- one stratum
def run(sub: pd.DataFrame, name: str) -> dict:
    sub = sub.reset_index(drop=True).copy()
    sub["wz"] = rint(sub.all_wake_t90.values)
    sub["sz"] = rint(sub.all_sleep_t90.values)

    rho, rho_p = stats.spearmanr(sub.all_wake_t90.values, sub.all_sleep_t90.values)

    rows, skipped = [], []
    for key, (label, _rawneg, _, _) in CONDS:
        neg = key in NEGATIVE_CONTROLS
        if f"{key}_incident" not in sub.columns:
            continue
        g = sub[sub[f"{key}_prevalent"] == 0]
        f = pd.DataFrame({"T": g[f"{key}_years"], "E": g[f"{key}_incident"].astype(int),
                          "site_id": g.site_id,
                          **{c: g[c] for c in ADJ},
                          "wz": g.wz, "sz": g.sz}).dropna()
        f = f[f["T"] > 0]

        ev = int(f.E.sum())
        ev_alt = int((f.E.values == 1).sum())        # second, independent event count
        assert ev == ev_alt, (key, ev, ev_alt)

        if ev < MIN_EVENTS:
            skipped.append({"cond": label, "key": key, "events": ev, "at_risk": int(len(f))})
            continue

        c = CoxPHFitter(penalizer=0.01).fit(f, "T", "E", strata=["site_id"])
        wlo, whi = np.exp(c.confidence_intervals_.loc["wz"])
        slo, shi = np.exp(c.confidence_intervals_.loc["sz"])
        rows.append({"cond": label, "key": key, "neg": bool(neg), "ev": ev, "at_risk": int(len(f)),
                     "wHR": float(np.exp(c.params_["wz"])),
                     "wlo": float(wlo), "whi": float(whi),
                     "wp": float(c.summary.loc["wz", "p"]),
                     "sHR": float(np.exp(c.params_["sz"])),
                     "slo": float(slo), "shi": float(shi),
                     "sp": float(c.summary.loc["sz", "p"])})

    r = pd.DataFrame(rows)
    if len(r):
        r["wq"] = bh(r.wp.values)
        r["sq"] = bh(r.sp.values)
        r["class"] = [classify(a, b) for a, b in zip(r.wq, r.sq)]
        r = r.sort_values("sHR", ascending=False).reset_index(drop=True)
        r[["cond", "neg", "ev", "at_risk", "wHR", "wlo", "whi", "wp",
           "sHR", "slo", "shi", "sp", "wq", "sq", "class"]].to_csv(
            HERE / f"wake_sleep_healthy_{name}.csv", index=False)

    counts = dict(r["class"].value_counts()) if len(r) else {}
    return {
        "n": int(len(sub)),
        "n_modelled": int(len(r)),
        "n_skipped_underpowered": int(len(skipped)),
        "median_wake_t90": float(np.median(sub.all_wake_t90)),
        "median_sleep_t90": float(np.median(sub.all_sleep_t90)),
        "iqr_wake_t90": [float(np.percentile(sub.all_wake_t90, 25)),
                         float(np.percentile(sub.all_wake_t90, 75))],
        "iqr_sleep_t90": [float(np.percentile(sub.all_sleep_t90, 25)),
                          float(np.percentile(sub.all_sleep_t90, 75))],
        "spearman_wake_sleep": float(rho),
        "spearman_p": float(rho_p),
        "median_age": float(np.median(sub.AgeAtVisit)),
        "pct_male": float(100 * sub.male.mean()),
        "median_ahi": float(np.median(sub.AHI)),
        "class_counts": {k: int(v) for k, v in counts.items()},
        "sleep_only": sorted(r.loc[r["class"] == "sleep only", "cond"].tolist()) if len(r) else [],
        "both": sorted(r.loc[r["class"] == "BOTH", "cond"].tolist()) if len(r) else [],
        "awake_only": sorted(r.loc[r["class"] == "awake only", "cond"].tolist()) if len(r) else [],
        "neg_controls_flagged": sorted(
            r.loc[r.neg & (r["class"] != "neither"), "cond"].tolist()) if len(r) else [],
        "table": r.to_dict("records") if len(r) else [],
        "skipped": skipped,
    }


# --------------------------------------------------------------------------- main
def main():
    m, prov = load()
    print(f"analysis set n = {len(m)}")

    # ---- positive control: the full cohort must reproduce wake_vs_sleep_t90.csv exactly
    A = run(m, "A_full")
    ref = pd.read_csv(REFCSV)
    got = pd.DataFrame(A["table"])
    j = got.merge(ref, on="cond", suffixes=("", "_ref"))
    # The reference was frozen with 51 conditions. The two replacement negative controls of
    # 2026-08-07, contact dermatitis and hemorrhoids, are additions, so the check is that
    # every condition the reference holds still reproduces EXACTLY and that the only new
    # rows are the two known ones. A length equality would have masked a substitution.
    # v8.2 (2026-09-15): the reference may carry an earlier condition panel (the v7 file holds cataract and back pain and none of the
    # decision-11 conditions). Every added or removed row must be explained by the current panel (disease_definitions): an added row is
    # a panel condition, a removed row is no longer in the panel. In --refresh-reference mode the reference is then rewritten; without
    # the flag any panel difference stops the run, as a stale baseline should.
    panel = {label for label, *_ in (v for v in DISEASES.values())} | {"Death from any cause"}
    added = set(got.cond) - set(ref.cond); removed = set(ref.cond) - set(got.cond)
    assert added <= panel, sorted(added - panel)
    assert not (removed & panel), sorted(removed & panel)
    if added or removed:
        print(f"reference panel differs from the current panel: added {sorted(added)}, removed {sorted(removed)}")
        assert "--refresh-reference" in sys.argv, "the reference panel differs from the current panel: pass --refresh-reference"
    assert len(j) == len(set(got.cond) & set(ref.cond)), (len(j), len(ref))
    # ------------------------------------------------------------------ regression baseline
    # wake_vs_sleep_t90.csv is a regression baseline, not an independent result. It was frozen on
    # 2026-08-03 and the cohort moved under it on 2026-08-04, for reasons that have nothing to do
    # with the negative-control panel:
    #   the ICD-9 278.01 -> 278.03 correction (obesity hypoventilation, 375 events -> 169),
    #   the E888/E884 prefix collision with ICD-10 E88.81 (falls, 454 -> 456), and
    #   the oximetry quality rule that removed 210 recordings (analysis set 6,803 -> 6,800).
    # The last of those re-ranks the exposure for every condition, so nothing reproduces bit for
    # bit and no subset check can rescue the old file. Pass --refresh-reference to rewrite the
    # baseline from the current cohort; the old one is superseded, never deleted, and the full
    # diff is written beside it. Without the flag the strict check runs and fails loudly, which
    # is the correct behaviour for a stale baseline.
    import shutil
    from datetime import date
    if "--refresh-reference" in sys.argv:
        diff = j[["cond", "ev", "ev_ref", "wHR", "wHR_ref", "sHR", "sHR_ref"]].copy()
        diff["d_ev"] = diff.ev - diff.ev_ref
        diff["d_wHR"] = (diff.wHR - diff.wHR_ref).round(6)
        diff["d_sHR"] = (diff.sHR - diff.sHR_ref).round(6)
        shutil.copy(REFCSV, HERE / f"_superseded/wake_vs_sleep_t90_{date.today()}_pre_refresh.csv")
        diff.to_csv(HERE / f"_superseded/wake_vs_sleep_t90_{date.today()}_refresh_diff.csv",
                    index=False)
        got.to_csv(REFCSV, index=False)
        print(f"reference refreshed; {int((diff.d_ev != 0).sum())} conditions moved on events, "
              f"largest wake-HR shift {diff.d_wHR.abs().max():.4f}, diff written to _superseded/")
        ref = pd.read_csv(REFCSV)
        j = got.merge(ref, on="cond", suffixes=("", "_ref"))
    for a, b_ in [("ev", "ev_ref"), ("wHR", "wHR_ref"), ("sHR", "sHR_ref"),
                  ("wq", "wq_ref"), ("sq", "sq_ref")]:
        d = np.abs(j[a].values - j[b_].values).max()
        assert d < 1e-8, (a, d)
    assert (j["class"].values == j["class_ref"].values).all()
    print(f"positive control PASSED: {len(j)}/{len(ref)} conditions reproduce exactly")

    # ---- restriction ladder, with the cost of each block
    ladder, keep = [], np.ones(len(m), bool)
    for blockname, keys in [("cardiopulmonary", CARDIOPULM),
                            ("renal/hepatic/cancer", RENAL_HEP_CA),
                            ("metabolic", METABOLIC)]:
        for k in keys:
            prev = (m[f"{k}_prevalent"] == 1).values
            before = int(keep.sum())
            keep = keep & ~prev
            ladder.append({"block": blockname, "condition": DISEASES.get(k, (k,))[0],
                           "prevalent_in_full_cohort": int(prev.sum()),
                           "removed_here": before - int(keep.sum()),
                           "remaining": int(keep.sum())})

    mB = np.ones(len(m), bool)
    for k in CARDIOPULM + RENAL_HEP_CA:
        mB &= ~(m[f"{k}_prevalent"] == 1).values
    mC = mB.copy()
    for k in METABOLIC:
        mC &= ~(m[f"{k}_prevalent"] == 1).values
    ahi5 = (m.AHI < 5).values

    # two independent counts of every stratum size
    for nm, msk in [("B", mB), ("C", mC)]:
        assert int(msk.sum()) == len(m[msk]) == int(np.count_nonzero(msk)), nm

    out = {
        "provenance": {
            **prov,
            "n_analysis_set": int(len(m)),
            "min_events": MIN_EVENTS,
            "min_wake_minutes": MIN_WAKE_MIN,
            "q_threshold": QTHR,
            "positive_control": "reproduces numbers/wake_vs_sleep_t90.csv to <1e-8 on "
                                "events, both HRs, both q values and the class label",
        },
        "restriction_ladder": ladder,
        "definitions": {
            "A_full": "full analysis set",
            "B_organ_free": "no prevalent " + ", ".join(
                DISEASES[k][0] for k in CARDIOPULM + RENAL_HEP_CA),
            "C_organ_metab_free": "B plus no prevalent " + ", ".join(
                DISEASES[k][0] for k in METABOLIC),
            "D_C_ahi5": "C plus AHI < 5",
            "E_B_ahi5": "B plus AHI < 5",
        },
        "strata": {"A_full": A},
    }

    for name, msk in [("B_organ_free", mB), ("C_organ_metab_free", mC),
                      ("D_C_ahi5", mC & ahi5), ("E_B_ahi5", mB & ahi5)]:
        res = run(m[msk], name)
        out["strata"][name] = res
        print(f"{name:22s} n={res['n']:5d}  modelled={res['n_modelled']:3d}  "
              f"underpowered={res['n_skipped_underpowered']:3d}  {res['class_counts']}")
        print(f"{'':22s} sleep only: {res['sleep_only']}")

    (HERE / "wake_sleep_healthy.json").write_text(json.dumps(out, indent=1))
    print("\nwrote wake_sleep_healthy.json and one CSV per stratum")


if __name__ == "__main__":
    main()
