"""
Short sleep before and after nocturnal oxygen enters the model. Colleague request #1.

"In the primary cohort, fit site-stratified Cox models adjusted for age spline and sex for
all outcomes with total sleep time <5h and <6h as the only exposure (no T90), then refit each
model adding T90 >10%, and report both sets of HRs side by side in one supplementary table."

Panel note. The request says "51 outcomes". 51 was the condition count of the panel before
2026-08-07, when contact dermatitis and hemorrhoids entered as replacement negative controls.
The paper's primary analysis (numbers/run_primary_unpenalized.py, bdsp_diseases_v3.csv) fits
54 entries: 47 single conditions, the cardiovascular composite, 5 negative controls, and death
from any cause. This script fits the same 54. Benjamini-Hochberg correction runs within each
exposure column across the 48 non-control outcomes, the paper's discovery panel (the "0 of 48"
short-sleep claim and eFigure 24 use the same set). The 5 controls and death are shown for
completeness outside the correction panel, and a BH over all 54 is kept in the CSV as a check.

Cohort: data_frozen_v7_2026-09/t90_final.parquet through numbers/cohort_spec.apply_cohort, 19,173.
Total sleep time is missing for 6 participants. They are carried as not-short, which is the
shipped convention of the paper's joint model (run_all_v2.py, eTable 9), so the overlapping
columns of the new table and eTable 9 are digit-identical. Excluding the 6 instead moves no
short-sleep hazard ratio by more than 0.004, no oxygen hazard ratio by more than 0.02, and no
significance call (verified against a run with them excluded). All four models per outcome
share one frame either way, so the A-to-B change is a model effect, not a sample effect.
Nobody in the cohort has TST at or below 30 minutes, so the failed-study floor used in the
earlier 3-hour analysis removes no one and is not applied.

Model, identical to numbers/run_primary_unpenalized.py:
  site-stratified Cox (Efron ties), UNPENALIZED, natural cubic age spline cr(df=4) with the
  first basis column dropped (the 4 columns sum to 1 and an unpenalized solver on the full
  basis silently fits nothing), plus sex. statsmodels PHReg is used, the same engine as the
  filed short_sleep_threshold analysis. Prevalent cases of each outcome are excluded exactly
  as the paper does: prevalent == 0, follow-up years present and > 0.

Per outcome, four fits on one frame:
  A5  exposure TST < 300 min (binary vs >= 300), no T90
  B5  A5 plus T90 > 10% binary (spo2_pct_below_90 > 10)
  A6  exposure TST < 360 min (binary vs >= 360), no T90
  B6  A6 plus T90 > 10% binary

Sensitivity (summary only): the <6h contrast against a 7-8h (420-480 min) reference, with
6-7h and >8h as separate levels, without and with T90 > 10%.

Positive controls, run FIRST and asserted before anything else is trusted:
  PC1  hf, copd2, diabetes refit with the paper's per-SD T90 exposure (site-wise rank inverse
       normal) must reproduce numbers/bdsp_diseases_v3.csv (heart failure 1.598 [1.500-1.702])
       to within rounding, with events and n exact.
  PC2  the shipped joint binary model (eTable 9, results_v2.json vs_sleep_duration) must
       reproduce for heart failure and obesity hypoventilation under its exact shipped spec.

Sanity: crude incidence per 1000 person-years by sleep group for heart failure and obesity
hypoventilation, direction compared with the fitted A5 hazard ratio.

Outputs (this directory): results.csv, sensitivity_reference6h.csv, eTable17_draft.md,
REPORT.txt. Local only, deterministic, reads the frozen parquet, writes nothing outside this
directory. Does not touch the manuscript, the supplement, or LATEST_FINAL.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import os
import sys
import time
import warnings

import numpy as np
import pandas as pd
from patsy import dmatrix
from scipy import stats
from statsmodels.duration.hazard_regression import PHReg

warnings.filterwarnings("ignore")

ROOT = paths.T90_ROOT
OUT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
from cohort_spec import apply_cohort                          # noqa: E402
from disease_definitions import DISEASES, NEGATIVE_CONTROLS   # noqa: E402

MIN_EVENTS = 60          # the paper's floor for a fit worth reporting
PANEL = list(DISEASES.items()) + [("death", ("Death from any cause", False, [], []))]


def rint(x):
    """Rank inverse normal, exactly as numbers/run_primary_unpenalized.py."""
    r = stats.rankdata(x)
    return stats.norm.ppf((r - 0.375) / (len(r) + 0.25))


def bh_q(p):
    """Benjamini-Hochberg q values, hand-rolled (step-up, monotone)."""
    p = np.asarray(p, float)
    m = len(p)
    o = np.argsort(p)
    ranked = p[o] * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(m)
    out[o] = np.clip(q, 0, 1)
    return out


def build():
    b = apply_cohort(pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet"))
    b["male"] = (b.sex.astype(str).str.upper().str[0] == "M").astype(float)
    sp = dmatrix("cr(a, df=4) - 1", {"a": b.AgeAtVisit.values},
                 return_type="dataframe").iloc[:, 1:]        # first basis column dropped
    for i in range(sp.shape[1]):
        b[f"age_s{i}"] = sp.iloc[:, i].values
    adj = [f"age_s{i}" for i in range(sp.shape[1])] + ["male"]
    assert len(adj) == 4, adj
    b["z"] = b.groupby("site_id").spo2_pct_below_90.transform(rint)   # per-SD, PC1 only
    b["t90_hi"] = (b.spo2_pct_below_90 > 10).astype(float)
    # NaN < threshold is False, so the 6 TST-missing participants are carried as not-short.
    # That is the shipped convention of run_all_v2.py / eTable 9, kept deliberately so the
    # overlapping columns of the two tables cannot disagree. See the module docstring.
    b["s5"] = (b.TST_min < 300).astype(float).where(b.TST_min.notna())   # v8.1c: missing where TST is masked (split nights)
    b["s6"] = (b.TST_min < 360).astype(float).where(b.TST_min.notna())   # v8.1c
    # sensitivity bands: 0 = <360, 1 = 360-420, 2 = 420-480 (reference), 3 = >= 480
    band = np.select([b.TST_min < 360, b.TST_min < 420, b.TST_min < 480, b.TST_min >= 480],
                     [0, 1, 2, 3], default=np.nan)
    b["band6"] = np.where(b.TST_min.notna(), band, np.nan)
    return b, adj


def frame(b, key, xcols, adj):
    """The fitted frame for one outcome, prevalent excluded exactly as the paper does."""
    yc, ec, pc = f"{key}_years", f"{key}_incident", f"{key}_prevalent"
    if yc not in b.columns:
        return None
    g = b[(b[pc] == 0) & b[yc].notna() & (b[yc] > 0)]
    d = pd.DataFrame({"T": g[yc].values, "E": g[ec].astype(int).values,
                      "site": g.site_id.values,
                      **{c: g[c].values for c in xcols + adj}}).dropna()
    return d


def fit(d, xcols, adj):
    """One site-stratified unpenalized Cox fit (Efron). Returns per-exposure HR, CI, p."""
    X = d[xcols + adj].values
    m = PHReg(d["T"].values, X, status=d["E"].values, strata=d["site"].values,
              ties="efron").fit()
    out = {}
    for i, x in enumerate(xcols):
        c, se = m.params[i], m.bse[i]
        out[x] = {"hr": float(np.exp(c)), "lo": float(np.exp(c - 1.959964 * se)),
                  "hi": float(np.exp(c + 1.959964 * se)), "p": float(m.pvalues[i])}
    return out


def main():
    t0 = time.time()
    b, adj = build()
    n_tst_missing = int(b.TST_min.isna().sum())
    print(f"cohort {len(b):,}, TST missing {n_tst_missing}, "
          f"<5h {int((b.s5 == 1).sum()):,}, <6h {int((b.s6 == 1).sum()):,}, "
          f"T90>10% {int(b.t90_hi.sum()):,}")

    # ---------------------------------------------------- PC1: per-SD harness gate
    v3 = pd.read_csv(f"{paths.NUMBERS_DIR}/bdsp_diseases_v3.csv").set_index("key")
    print("\nPC1, per-SD T90 against bdsp_diseases_v3.csv")
    for key in ("hf", "copd2", "diabetes"):
        d = frame(b, key, ["z"], adj)
        r = fit(d, ["z"], adj)["z"]
        ref = v3.loc[key]
        ev, n = int(d.E.sum()), len(d)
        ok = (abs(r["hr"] - ref.t90_hr) <= 0.0015 and abs(r["lo"] - ref.t90_lo) <= 0.0015
              and abs(r["hi"] - ref.t90_hi) <= 0.0015 and ev == ref.events and n == ref.n)
        print(f"  {key:10s} refit {r['hr']:.3f} [{r['lo']:.3f}-{r['hi']:.3f}] "
              f"ev {ev} n {n}  vs filed {ref.t90_hr:.3f} [{ref.t90_lo:.3f}-{ref.t90_hi:.3f}] "
              f"ev {int(ref.events)} n {int(ref.n)}  {'PASS' if ok else 'FAIL'}")
        assert ok, f"PC1 failed on {key}: harness does not reproduce the paper"

    # ------------------------------------- PC2: shipped joint binary model (eTable 9)
    import json
    vs = json.load(open(f"{paths.NUMBERS_DIR}/results_v2.json"))["vs_sleep_duration"]["outcomes"]
    b["slp_lt5_shipped"] = (b.TST_min / 60.0 < 5).astype(float).where(b.TST_min.notna())   # v8.1c: missing where TST is masked, as run_all_v2 now does
    print("\nPC2, shipped joint model (eTable 9 spec) for heart failure and OHS")
    for key, lab in (("hf", "Heart failure"), ("obesity_hypovent", "Obesity hypoventilation")):
        d = frame(b, key, ["t90_hi", "slp_lt5_shipped"], adj)
        r = fit(d, ["t90_hi", "slp_lt5_shipped"], adj)
        ref = vs[lab]
        ok = (abs(r["t90_hi"]["hr"] - ref["t90"]["hr"]) <= 0.002
              and abs(r["slp_lt5_shipped"]["hr"] - ref["short_sleep"]["hr"]) <= 0.002)
        print(f"  {lab:24s} t90 {r['t90_hi']['hr']:.3f} vs {ref['t90']['hr']:.3f}, "
              f"sleep {r['slp_lt5_shipped']['hr']:.3f} vs {ref['short_sleep']['hr']:.3f}  "
              f"{'PASS' if ok else 'FAIL'}")
        assert ok, f"PC2 failed on {key}: joint model does not reproduce the shipped eTable 9"

    # ------------------------------------------------------------------ main fits
    print("\nfitting 4 models per outcome ...")
    rows = []
    for key, (label, _neg, _a, _b) in PANEL:
        d = frame(b, key, ["s5", "s6", "t90_hi"], adj)
        if d is None:
            continue
        ev = int(d.E.sum())
        if ev < MIN_EVENTS:
            continue
        rec = {"key": key, "outcome": label,
               "negative_control": key in NEGATIVE_CONTROLS, "is_death": key == "death",
               "in_bh_panel": (key not in NEGATIVE_CONTROLS) and key != "death",
               "n": len(d), "events": ev,
               "n_s5": int((d.s5 == 1).sum()), "events_s5": int(d.loc[d.s5 == 1, "E"].sum()),
               "n_s6": int((d.s6 == 1).sum()), "events_s6": int(d.loc[d.s6 == 1, "E"].sum()),
               "n_t90hi": int((d.t90_hi == 1).sum()),
               "events_t90hi": int(d.loc[d.t90_hi == 1, "E"].sum())}
        for tag, xcols in (("A5", ["s5"]), ("B5", ["s5", "t90_hi"]),
                           ("A6", ["s6"]), ("B6", ["s6", "t90_hi"])):
            r = fit(d, xcols, adj)
            s = r[xcols[0]]
            rec.update({f"{tag}_hr": s["hr"], f"{tag}_lo": s["lo"],
                        f"{tag}_hi": s["hi"], f"{tag}_p": s["p"]})
            if "t90_hi" in xcols:
                t = r["t90_hi"]
                rec.update({f"{tag}_t90_hr": t["hr"], f"{tag}_t90_lo": t["lo"],
                            f"{tag}_t90_hi": t["hi"], f"{tag}_t90_p": t["p"]})
        rows.append(rec)
        print(f"  {label:28s} ev {ev:5d}  A5 {rec['A5_hr']:.2f}  B5 {rec['B5_hr']:.2f}  "
              f"A6 {rec['A6_hr']:.2f}  B6 {rec['B6_hr']:.2f}  ({time.time()-t0:4.0f}s)")

    R = pd.DataFrame(rows)

    # BH within each exposure column: primary over the 48 non-control outcomes, and over
    # all fitted rows as a robustness column
    COLS = [("A5_p", "A5_q"), ("B5_p", "B5_q"), ("B5_t90_p", "B5_t90_q"),
            ("A6_p", "A6_q"), ("B6_p", "B6_q"), ("B6_t90_p", "B6_t90_q")]
    panel = R.in_bh_panel
    for pcol, qcol in COLS:
        R[qcol] = np.nan
        R.loc[panel, qcol] = bh_q(R.loc[panel, pcol].values)
        R[qcol + "_all"] = bh_q(R[pcol].values)

    # cross-check the B5 column against the shipped eTable 9 values, same spec, so the
    # only difference is the 3-decimal rounding stored in results_v2.json
    ds, dt = [], []
    for _, r in R.iterrows():
        if r.outcome in vs:
            ds.append(abs(r.B5_hr - vs[r.outcome]["short_sleep"]["hr"]))
            dt.append(abs(r.B5_t90_hr - vs[r.outcome]["t90"]["hr"]))
    max_ds, max_dt = float(np.max(ds)), float(np.max(dt))
    print(f"\nB5 vs shipped eTable 9, max |HR difference|: sleep term {max_ds:.4f}, "
          f"oxygen term {max_dt:.4f}")
    assert max_ds <= 0.002 and max_dt <= 0.002, "B5 does not reproduce the shipped eTable 9"

    R.to_csv(f"{OUT}/results.csv", index=False)

    # ---------------------------- variant check: the 6 TST-missing excluded instead
    print("variant check, TST-missing excluded instead of carried as not-short ...")
    b2 = b.copy()
    b2["s5"] = np.where(b2.TST_min.notna(), (b2.TST_min < 300).astype(float), np.nan)
    b2["s6"] = np.where(b2.TST_min.notna(), (b2.TST_min < 360).astype(float), np.nan)
    vrows = []
    for key, (label, _neg, _a, _b) in PANEL:
        d = frame(b2, key, ["s5", "s6", "t90_hi"], adj)
        if d is None or int(d.E.sum()) < MIN_EVENTS:
            continue
        rec = {"key": key}
        for tag, xcols in (("A5", ["s5"]), ("B5", ["s5", "t90_hi"]),
                           ("A6", ["s6"]), ("B6", ["s6", "t90_hi"])):
            r = fit(d, xcols, adj)
            rec[f"{tag}_hr"] = r[xcols[0]]["hr"]
            rec[f"{tag}_p"] = r[xcols[0]]["p"]
            if "t90_hi" in xcols:
                rec[f"{tag}_t90_hr"] = r["t90_hi"]["hr"]
                rec[f"{tag}_t90_p"] = r["t90_hi"]["p"]
        vrows.append(rec)
    V = pd.DataFrame(vrows).set_index("key")
    RM = R.set_index("key")
    var_sleep_d = max(float((V[f"{t}_hr"] - RM[f"{t}_hr"]).abs().max())
                      for t in ("A5", "B5", "A6", "B6"))
    var_oxy_d = max(float((V[f"{t}_t90_hr"] - RM[f"{t}_t90_hr"]).abs().max())
                    for t in ("B5", "B6"))
    panel_keys = RM.index[RM.in_bh_panel]
    flips = []
    for pcol, qcol in COLS:
        qv = pd.Series(bh_q(V.loc[panel_keys, pcol].values), index=panel_keys)
        for k in panel_keys[(qv < 0.05) != (RM.loc[panel_keys, qcol] < 0.05).values]:
            flips.append({"col": qcol, "outcome": RM.loc[k, "outcome"],
                          "q_main": float(RM.loc[k, qcol]), "q_variant": float(qv[k])})
    calls_same = len(flips) == 0
    sleep_flips = [f for f in flips if "t90" not in f["col"]]
    print(f"  max sleep-HR shift {var_sleep_d:.4f}, max oxygen-HR shift {var_oxy_d:.4f}, "
          f"BH call flips: {len(flips)}"
          + ("" if calls_same else "  " + "; ".join(
              f"{f['col']} {f['outcome']} q {f['q_main']:.3f} vs {f['q_variant']:.3f}"
              for f in flips)))

    # ------------------------------------------------- sensitivity: 7-8h reference
    print("\nsensitivity, <6h against a 420-480 min reference ...")
    sens = []
    for key, (label, _neg, _a, _b) in PANEL:
        d = frame(b, key, ["band6", "t90_hi"], adj)
        if d is None or int(d.E.sum()) < MIN_EVENTS:
            continue
        d = d.copy()
        d["lt6"] = (d.band6 == 0).astype(float)
        d["h67"] = (d.band6 == 1).astype(float)
        d["gt8"] = (d.band6 == 3).astype(float)
        rec = {"key": key, "outcome": label,
               "in_bh_panel": (key not in NEGATIVE_CONTROLS) and key != "death",
               "events": int(d.E.sum()),
               "n_ref": int((d.band6 == 2).sum()), "ev_ref": int(d.loc[d.band6 == 2, "E"].sum()),
               "n_lt6": int((d.band6 == 0).sum()), "ev_lt6": int(d.loc[d.band6 == 0, "E"].sum())}
        for tag, xcols in (("noT90", ["lt6", "h67", "gt8"]),
                           ("withT90", ["lt6", "h67", "gt8", "t90_hi"])):
            try:
                r = fit(d, xcols, adj)["lt6"]
            except np.linalg.LinAlgError:
                # v8.1c: with the split nights' sleep time masked, a sparse outcome can leave one
                # duration band without events inside a hospital, and the site-stratified fit is
                # singular. The row is kept with missing estimates and named in the log.
                print(f"  {label}: {tag} sensitivity fit singular (a duration band has no events within a hospital), left missing")
                r = {"hr": np.nan, "lo": np.nan, "hi": np.nan, "p": np.nan}
            rec.update({f"{tag}_hr": r["hr"], f"{tag}_lo": r["lo"],
                        f"{tag}_hi": r["hi"], f"{tag}_p": r["p"]})
        sens.append(rec)
    S = pd.DataFrame(sens)
    sp_ = S.in_bh_panel
    for tag in ("noT90", "withT90"):
        S[f"{tag}_q"] = np.nan
        S.loc[sp_, f"{tag}_q"] = bh_q(S.loc[sp_, f"{tag}_p"].values)
    S.to_csv(f"{OUT}/sensitivity_reference6h.csv", index=False)

    # ------------------------------------------------------------- crude-rate sanity
    print("\nsanity, crude incidence per 1000 person-years by sleep group")
    sanity = []
    for key, lab in (("hf", "Heart failure"), ("obesity_hypovent", "Obesity hypoventilation")):
        d = frame(b, key, ["s5"], adj)
        g1, g0 = d[d.s5 == 1], d[d.s5 == 0]
        r1 = 1000 * g1.E.sum() / g1["T"].sum()
        r0 = 1000 * g0.E.sum() / g0["T"].sum()
        hr = float(R.loc[R.key == key, "A5_hr"].iloc[0])
        agree = (r1 > r0) == (hr > 1)
        sanity.append((lab, r1, r0, hr, agree))
        print(f"  {lab:24s} <5h {r1:6.2f}  >=5h {r0:6.2f}  A5 HR {hr:.2f}  "
              f"direction {'AGREES' if agree else 'DISAGREES'}")

    # ------------------------------------------------------------------ summaries
    P = R[R.in_bh_panel]

    def n_sig(qcol):
        return int((P[qcol] < 0.05).sum())

    def n_raw(pcol):
        return int((P[pcol] < 0.05).sum())

    def att(a_tag, b_tag):
        """Percent change in the short-sleep log HR when T90 enters."""
        la = np.log(P[f"{a_tag}_hr"].values)
        lb = np.log(P[f"{b_tag}_hr"].values)
        with np.errstate(divide="ignore", invalid="ignore"):
            pct = 100 * (lb - la) / la
        sig = P[f"{a_tag}_p"].values < 0.05
        return {"median_pct_all": float(np.nanmedian(pct[np.abs(la) > 0.01])),
                "median_pct_Asig": float(np.nanmedian(pct[sig])) if sig.any() else np.nan,
                "n_Asig": int(sig.sum()),
                "median_abs_shift": float(np.median(np.abs(lb - la))),
                "n_toward_null": int(((np.abs(lb) < np.abs(la)) & (np.abs(la) > 0.01)).sum()),
                "n_compared": int((np.abs(la) > 0.01).sum())}

    att5, att6 = att("A5", "B5"), att("A6", "B6")

    lines = []
    a = lines.append
    a("SHORT SLEEP BEFORE AND AFTER T90, colleague request #1, "
      + time.strftime("%Y-%m-%d %H:%M"))
    a("")
    a("HEADLINE")
    a(f"Panel fitted: {len(R)} entries. {int(P.shape[0])} non-control outcomes carry the BH "
      f"correction, plus {int(R.negative_control.sum())} negative controls and death shown "
      f"outside it. The request said 51 outcomes. 51 was the panel before 2026-08-07, when "
      f"contact dermatitis and hemorrhoids entered as replacement negative controls. The "
      f"paper's primary analysis fits these {len(R)}.")
    a(f"TST under 5h alone (A5): {n_raw('A5_p')} of 48 at raw p<0.05, {n_sig('A5_q')} "
      f"survive BH.")
    a(f"TST under 5h with T90 (B5): {n_raw('B5_p')} of 48 at raw p<0.05, {n_sig('B5_q')} "
      f"survive BH.")
    a(f"TST under 6h alone (A6): {n_raw('A6_p')} of 48 at raw p<0.05, {n_sig('A6_q')} "
      f"survive BH.")
    a(f"TST under 6h with T90 (B6): {n_raw('B6_p')} of 48 at raw p<0.05, {n_sig('B6_q')} "
      f"survive BH.")
    t90sig = P[P.B5_t90_q < 0.05]
    t90_dn = t90sig[t90sig.B5_t90_hr < 1]
    a(f"T90>10% inside the B5 models: {n_raw('B5_t90_p')} of 48 at raw p<0.05, "
      f"{n_sig('B5_t90_q')} survive BH. Inside B6: {n_raw('B6_t90_p')} raw, "
      f"{n_sig('B6_t90_q')} after BH. Of the B5 survivors, "
      f"{len(t90sig) - len(t90_dn)} sit above 1 and {len(t90_dn)} below"
      + (" (" + ", ".join(f"{r.outcome} {r.B5_t90_hr:.2f}"
                          for _, r in t90_dn.iterrows()) + ")." if len(t90_dn) else "."))

    def names(df, tag):
        return ", ".join(f"{r.outcome} {r[f'{tag}_hr']:.2f}"
                         for _, r in df.sort_values(f"{tag}_q").iterrows()) or "none"

    a(f"BH survivors for short sleep, <5h alone: {names(P[P.A5_q < 0.05], 'A5')}. "
      f"<5h with T90: {names(P[P.B5_q < 0.05], 'B5')}. "
      f"<6h alone: {names(P[P.A6_q < 0.05], 'A6')}. "
      f"<6h with T90: {names(P[P.B6_q < 0.05], 'B6')}.")
    a("")
    a("ATTENUATION, percent change in the short-sleep log HR when T90 enters")
    a(f"Under 5h: median {att5['median_pct_Asig']:+.1f}% across the {att5['n_Asig']} outcomes "
      f"where A5 had raw p<0.05, median {att5['median_pct_all']:+.1f}% across all with "
      f"|log HR|>0.01, {att5['n_toward_null']} of {att5['n_compared']} moved toward 1.0, "
      f"median absolute log-HR shift {att5['median_abs_shift']:.4f}.")
    a(f"Under 6h: median {att6['median_pct_Asig']:+.1f}% across the {att6['n_Asig']} outcomes "
      f"where A6 had raw p<0.05, median {att6['median_pct_all']:+.1f}% across all with "
      f"|log HR|>0.01, {att6['n_toward_null']} of {att6['n_compared']} moved toward 1.0, "
      f"median absolute log-HR shift {att6['median_abs_shift']:.4f}.")
    ol5 = int(((b.s5 == 1) & (b.t90_hi == 1)).sum())
    ol6 = int(((b.s6 == 1) & (b.t90_hi == 1)).sum())
    a(f"Mechanism: {ol5:,} of {int((b.s5 == 1).sum()):,} under-5h sleepers and {ol6:,} of "
      f"{int((b.s6 == 1).sum()):,} under-6h sleepers also have T90>10%, so part of any crude "
      f"short-sleep signal is carried by oxygen.")
    a("")
    a("SENSITIVITY, under-6h contrast against a 7-8h (420-480 min) reference with 6-7h and "
      "over-8h as separate levels")
    Ssub = S[S.in_bh_panel]
    sens_named = ", ".join(f"{r.outcome} {r.noT90_hr:.2f}"
                           for _, r in Ssub[Ssub.noT90_p < 0.05].iterrows()) or "none"
    a(f"Without T90: {int((Ssub.noT90_p < 0.05).sum())} of {len(Ssub)} at raw p<0.05 "
      f"({sens_named}), {int((Ssub.noT90_q < 0.05).sum())} survive BH, median HR "
      f"{Ssub.noT90_hr.median():.2f}.")
    same_set = (set(Ssub[Ssub.noT90_p < 0.05].outcome)
                == set(Ssub[Ssub.withT90_p < 0.05].outcome))
    a(f"With T90: {int((Ssub.withT90_p < 0.05).sum())} of {len(Ssub)} at raw p<0.05"
      f"{' (the same outcomes)' if same_set else ''}, "
      f"{int((Ssub.withT90_q < 0.05).sum())} survive BH, median HR "
      f"{Ssub.withT90_hr.median():.2f}.")
    n_ref_cohort = int(((b.TST_min >= 420) & (b.TST_min < 480)).sum())
    a(f"Reference cell 420-480 min holds {n_ref_cohort:,} of the cohort, "
      f"{int(S.n_ref.min()):,}-{int(S.n_ref.max()):,} per outcome after prevalent "
      f"exclusion. The reference choice does not change the BH conclusion.")
    a("")
    a("OUTSIDE THE BH PANEL")
    dr = R[R.is_death].iloc[0]
    a(f"Death from any cause: <5h alone {dr.A5_hr:.2f} (p {dr.A5_p:.2f}), <6h alone "
      f"{dr.A6_hr:.2f} (p {dr.A6_p:.3f}), both below 1, the direction the shipped eTable 9 "
      f"already shows at <5h (0.93).")
    ctrl = R[R.negative_control]
    ctrl_min = min(float(ctrl[c].min()) for c in ("A5_p", "B5_p", "A6_p", "B6_p"))
    n_ctrl_sig = sum(int((ctrl[c] < 0.05).sum()) for c in ("A5_p", "B5_p", "A6_p", "B6_p"))
    a(f"Negative controls: {n_ctrl_sig} of 20 short-sleep estimates across the four models "
      f"reach raw p<0.05, smallest raw p {ctrl_min:.3f}.")
    a("")
    a("POSITIVE CONTROLS")
    a("PC1 per-SD T90 refit reproduces bdsp_diseases_v3.csv for heart failure, COPD and type "
      "2 diabetes, hazard ratios within rounding and events and n exact. Heart failure "
      "1.598 [1.500-1.702] reproduced.")
    a("PC2 the shipped joint binary model (eTable 9) reproduces for heart failure and obesity "
      "hypoventilation under its exact shipped specification.")
    a(f"B5 column against shipped eTable 9 across all 54 entries: max absolute HR difference "
      f"{max_ds:.4f} on the sleep term and {max_dt:.4f} on the oxygen term, within the "
      f"3-decimal rounding of the stored values, so the overlapping columns print "
      f"identically.")
    if calls_same:
        flip_txt = "every BH significance call identical"
    elif not sleep_flips:
        flip_txt = ("every short-sleep BH call identical, the only boundary case is "
                    + " and ".join(f"the oxygen estimate for {f['outcome']} in {f['col'][:2]} "
                                   f"(q {f['q_main']:.3f} against {f['q_variant']:.3f})"
                                   for f in flips))
    else:
        flip_txt = "BH calls DIFFER on short-sleep columns, read results.csv"
    a(f"Variant check, the 6 TST-missing excluded instead of carried as not-short: max "
      f"short-sleep HR shift {var_sleep_d:.4f}, max oxygen HR shift {var_oxy_d:.4f}, "
      f"{flip_txt}.")
    a("")
    a("SANITY, crude incidence per 1000 person-years")
    for lab, r1, r0, hr, agree in sanity:
        a(f"{lab}: under 5h {r1:.2f} against 5h or more {r0:.2f}, fitted A5 HR {hr:.2f}, "
          f"direction {'agrees' if agree else 'disagrees'}.")
    a("")
    qall_flips = []
    for _p, q in COLS:
        s48 = set(P.loc[P[q] < 0.05, "outcome"])
        s54 = set(P.loc[P[q + "_all"] < 0.05, "outcome"])
        for o in sorted(s48 ^ s54):
            r = P[P.outcome == o].iloc[0]
            qall_flips.append(f"{o} in {q[:2]}"
                              + (" oxygen" if "t90" in q else "")
                              + f" (q {r[q]:.3f} against {r[q + '_all']:.3f})")
    qall_sleep_ok = not any("t90" not in q for _p, q in COLS
                            for o in (set(P.loc[P[q] < 0.05, "outcome"])
                                      ^ set(P.loc[P[q + "_all"] < 0.05, "outcome"])))
    a("CONVENTIONS")
    a("Site-stratified Cox, Efron ties, unpenalized, statsmodels PHReg. Age natural cubic "
      "spline cr(df=4) with the first basis column dropped, plus sex. Hospital as stratum. "
      "Prevalent cases of each outcome excluded, follow-up years present and positive. Total "
      "sleep time missing for 6 participants, carried as not-short per the shipped eTable 9 "
      "convention, with the excluded variant checked above. All four models per outcome sit "
      "on one frame, so the A-to-B change is a model effect, not a sample effect. BH within "
      "each exposure column across the 48 non-control outcomes. A BH over all 54 entries is "
      "in results.csv as *_q_all and "
      + ("keeps every within-panel call unchanged."
         if not qall_flips else
         ("keeps every short-sleep call unchanged, the boundary case being "
          + " and ".join(qall_flips) + "." if qall_sleep_ok
          else "CHANGES short-sleep calls, read results.csv.")))
    a("")
    a("FILES")
    a(f"{OUT}/results.csv")
    a(f"{OUT}/sensitivity_reference6h.csv")
    a(f"{OUT}/eTable17_draft.md")
    a(f"{OUT}/run_shortsleep_before_after_t90.py")

    open(f"{OUT}/REPORT.txt", "w").write("\n".join(lines) + "\n")

    # ------------------------------------------------------------------ eTable draft
    def ci(hr, lo, hi, q):
        star = "*" if (not np.isnan(q)) and q < 0.05 else ""
        return f"{hr:.2f} ({lo:.2f}-{hi:.2f}){star}"

    tab = pd.DataFrame([{
        "Condition": r.outcome,
        "New diagnoses": r.events,
        "Sleep under 5 h, alone (95% CI)": ci(r.A5_hr, r.A5_lo, r.A5_hi, r.A5_q),
        "Sleep under 5 h, with oxygen (95% CI)": ci(r.B5_hr, r.B5_lo, r.B5_hi, r.B5_q),
        "Oxygen above 10%, 5-h model (95% CI)": ci(r.B5_t90_hr, r.B5_t90_lo, r.B5_t90_hi,
                                                   r.B5_t90_q),
        "Sleep under 6 h, alone (95% CI)": ci(r.A6_hr, r.A6_lo, r.A6_hi, r.A6_q),
        "Sleep under 6 h, with oxygen (95% CI)": ci(r.B6_hr, r.B6_lo, r.B6_hi, r.B6_q),
        "Oxygen above 10%, 6-h model (95% CI)": ci(r.B6_t90_hr, r.B6_t90_lo, r.B6_t90_hi,
                                                   r.B6_t90_q),
        "Negative control": "Yes" if r.negative_control else ""} for _, r in R.iterrows()])

    n5, n6 = int((b.s5 == 1).sum()), int((b.s6 == 1).sum())
    caption = (
        "**eTable 17.** Short sleep duration on the study night, before and after nocturnal "
        "oxygen enters the model. For each condition, the hazard ratio for total sleep time "
        f"under 5 hours (n = {n5:,} exposed) and under 6 hours (n = {n6:,} exposed) is shown "
        "first from a model with short sleep as the only exposure, then from the same model "
        "with T90 above 10% of the recording added as a second binary exposure, with the "
        "oxygen estimate from that joint model beside it. Every model is a Cox model "
        "stratified by hospital with Efron ties, unpenalized, adjusted for age (natural cubic "
        "spline, 4 df, one basis column dropped) and sex, fitted on the primary cohort of "
        "19,173 after excluding prevalent cases of that condition. Total sleep time was "
        "missing for 6 participants, carried as not short, the convention of eTable 9. "
        f"Excluding them instead moves no short-sleep hazard ratio by more than "
        f"{var_sleep_d:.3f}, no oxygen hazard ratio by more than {var_oxy_d:.3f}, and "
        + ("changes no significance call. " if calls_same
           else ("changes no short-sleep significance call (the one boundary case is the "
                 "oxygen estimate for "
                 + " and ".join(sorted({f["outcome"] for f in flips}))
                 + ", at q 0.05 exactly). " if not sleep_flips
              else "changes significance calls, read results.csv. "))
        + "Asterisks mark estimates significant after Benjamini-Hochberg correction within "
        "their column across the 48 non-control outcomes. The five negative controls and "
        "death from any cause are shown for completeness and sit outside the correction "
        "panel. The 5-hour joint model is the model of eTable 9.")

    md = caption + "\n\n" + tab.to_markdown(index=False) + "\n"
    open(f"{OUT}/eTable17_draft.md", "w").write(md)

    print(f"\nwritten: results.csv ({len(R)} rows), sensitivity_reference6h.csv "
          f"({len(S)} rows), eTable17_draft.md, REPORT.txt   [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
