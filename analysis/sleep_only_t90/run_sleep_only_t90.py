#!/usr/bin/env python3
"""Does T90 hold up when measured during SLEEP only?  (Alen, 2026-09-02)

Two questions, both inside the 6,800 participants whose recordings carry stage-resolved
oximetry (the wake-vs-sleep subset of Fig. 4a), so recording-based and sleep-only T90 are
compared in the SAME people:

  1. The four T90 categories (<=1, 1-5, 5-10, >10 % ) re-cut on sleep-only T90, graded
     hazard ratios versus the <=1 % group, fitted exactly as numbers/run_all_v2.py fits the
     published bands (Cox, age spline cr(df=4) with one column dropped, sex, stratified by
     hospital, unpenalized, >= 60 events).
  2. The held-out concordance gain over age and sex, computed exactly as
     numbers/build_ranking_v3.py computes the ranking (site-wise rank-normal score, Cox
     penalizer 0.01 fitted in one hospital and scored in the other, both ways, averaged),
     for recording T90, sleep-only T90, wake-only T90 and total sleep time.

Positive controls, asserted before anything is reported:
  - the published graded heart-failure bands reproduce from the full cohort
    (results_v2.json graded: 1.413 / 1.654 / 2.738);
  - the subset is the same 6,800 as numbers/wake_sleep_healthy_A_full.csv.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../..")))  # repository root, where paths.py lives
import paths
import warnings; warnings.filterwarnings("ignore")
import glob, json, os, sys, time
import numpy as np, pandas as pd
from scipy import stats
from patsy import dmatrix
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index

ROOT = paths.T90_ROOT
NUM = f"{paths.NUMBERS_DIR}"
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, paths.ANALYSIS_DIR)
from disease_definitions import DISEASES, NEGATIVE_CONTROLS

SHARDS = f"{paths.T90_ROOT}/cpap_stage/out*.csv"
T0 = time.time()
def say(*a): print(f"[{time.time()-T0:5.0f}s]", *a, flush=True)

# ------------------------------------------------------------------ cohort (run_all_v2 filter)
b = pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet")
b = b[(b.fu_valid == 1) & b.spo2_pct_below_90.notna() & (b.oximetry_bad == 0)].copy()
assert len(b) == 19173, len(b)
b["t90"] = b.spo2_pct_below_90
b["male"] = (b.sex.astype(str).str.upper().str[0] == "M").astype(int)
sp = dmatrix("cr(a, df=4) - 1", {"a": b.AgeAtVisit.values}, return_type="dataframe")
sp.columns = [f"age_s{i}" for i in range(4)]; sp.index = b.index
for c in sp.columns: b[c] = sp[c]
ADJ4 = [f"age_s{i}" for i in range(4)] + ["male"]        # ranking (ridge) keeps all four
ADJ3 = [f"age_s{i}" for i in range(1, 4)] + ["male"]     # unpenalized band fits drop one
BANDS = [(-1, 1, "0-1%"), (1, 5, "1-5%"), (5, 10, "5-10%"), (10, 1e9, ">10%")]
def band_of(v): return next(i for i, (lo, hi, _) in enumerate(BANDS) if lo < v <= hi)
CONDS = list(DISEASES.items()) + [("death", ("Death from any cause", False, [], []))]
def label(k): return CONDS[[c[0] for c in CONDS].index(k)][1][0]

def graded(g, col, k, min_events=60):
    """run_all_v2.py's graded fit, verbatim in spirit: dummies for bands 1..3 vs band 0."""
    f = g[(g[f"{k}_prevalent"] == 0) & g[f"{k}_years"].notna() & (g[f"{k}_years"] > 0)].copy()
    f["band"] = f[col].apply(band_of)
    if int(f[f"{k}_incident"].sum()) < min_events: return None
    X = pd.get_dummies(f.band, prefix="bd").astype(int)
    for i in range(4):
        if f"bd_{i}" not in X: X[f"bd_{i}"] = 0
    X = X[[f"bd_{i}" for i in range(4)]].drop(columns=["bd_0"])
    ff = pd.concat([pd.DataFrame({"T": f[f"{k}_years"].values, "E": f[f"{k}_incident"].astype(int).values,
                                  "site": f.site_id.values, **{c: f[c].values for c in ADJ3}}),
                    X.reset_index(drop=True)], axis=1).dropna()
    try: c = CoxPHFitter(penalizer=0.0).fit(ff, "T", "E", strata=["site"])
    except Exception: return None
    row = {"events": int(ff.E.sum()), "n": int(len(ff)),
           "band_n": [int((f.band == i).sum()) for i in range(4)]}
    for i in range(1, 4):
        nm = f"bd_{i}"; lo, hi = np.exp(c.confidence_intervals_.loc[nm])
        row[BANDS[i][2]] = {"hr": float(np.exp(c.params_[nm])), "lo": float(lo), "hi": float(hi),
                            "p": float(c.summary.loc[nm, "p"])}
    ft = pd.DataFrame({"T": f[f"{k}_years"].values, "E": f[f"{k}_incident"].astype(int).values,
                       "band": f.band.values, "site": f.site_id.values, **{c_: f[c_].values for c_ in ADJ3}}).dropna()
    try:
        ct = CoxPHFitter(penalizer=0.0).fit(ft, "T", "E", strata=["site"])
        row["trend_p"] = float(ct.summary.loc["band", "p"]); row["trend_hr"] = float(np.exp(ct.params_["band"]))
    except Exception:
        row["trend_p"] = np.nan; row["trend_hr"] = np.nan
    return row

# ------------------------------------------------------------------ positive control 1
say("positive control: published heart-failure bands from the full cohort")
pc = graded(b, "t90", "hf")
pub = json.load(open(f"{NUM}/results_v2.json"))["graded"]["outcomes"]["Heart failure"]
for lab in ("1-5%", "5-10%", ">10%"):
    assert abs(pc[lab]["hr"] - pub[lab]["hr"]) < 0.0015, (lab, pc[lab]["hr"], pub[lab]["hr"])
assert pc["events"] == pub["events"] == 1187, (pc["events"], pub["events"])
say(f"  reproduced: {pc['1-5%']['hr']:.3f} / {pc['5-10%']['hr']:.3f} / {pc['>10%']['hr']:.3f}, {pc['events']} events")

# ------------------------------------------------------------------ stage-resolved subset (run_wake_sleep_healthy.load)
files = sorted(glob.glob(SHARDS)); assert len(files) == 72, len(files)
raw = pd.concat([pd.read_csv(f, low_memory=False) for f in files], ignore_index=True)
st = raw.drop_duplicates("BDSPPatientID"); st = st[st.status == "ok"]
bb = b[~(b.spo2_nadir_corrected > b.spo2_mean)]
m = bb.merge(st[["BDSPPatientID", "all_wake_t90", "all_sleep_t90", "all_wake_min", "all_sleep_min"]],
             on="BDSPPatientID", how="inner")
m = m[m.all_wake_min >= 10.0].dropna(subset=["all_wake_t90", "all_sleep_t90"]).reset_index(drop=True)
ref = pd.read_csv(f"{NUM}/wake_sleep_healthy_A_full.csv")
n_ref = int(ref["n"].max()) if "n" in ref.columns else None
say(f"subset n = {len(m):,} (wake_sleep_healthy_A_full n column max = {n_ref})")
assert len(m) == 6800, len(m)   # the count the manuscript prints (Fig. 4a subset)
sites = m.site_id.value_counts().to_dict(); say("  sites:", sites)
say(f"  recording T90 median {m.t90.median():.2f}, sleep-only T90 median {m.all_sleep_t90.median():.2f}, "
    f"wake T90 median {m.all_wake_t90.median():.2f}; Spearman recording vs sleep {stats.spearmanr(m.t90, m.all_sleep_t90)[0]:.3f}")
band_rec = m.t90.apply(band_of); band_slp = m.all_sleep_t90.apply(band_of)
xt = pd.crosstab(band_rec, band_slp); xt.index = [x[2] for x in BANDS]; xt.columns = [x[2] for x in BANDS]
say("  band cross-tab (rows recording, cols sleep-only):\n" + xt.to_string())
agree = float((band_rec == band_slp).mean())

# ------------------------------------------------------------------ 1. bands, recording vs sleep-only, same people
say("1. graded bands in the subset, recording-based vs sleep-only")
rows = []
for k, (lab, is_nc, *_ ) in CONDS:
    if f"{k}_incident" not in m.columns: continue
    r_rec = graded(m, "t90", k); r_slp = graded(m, "all_sleep_t90", k)
    if r_rec is None or r_slp is None: continue
    for src, r in (("recording", r_rec), ("sleep_only", r_slp)):
        rows.append({"outcome": lab, "key": k, "negative_control": k in NEGATIVE_CONTROLS, "definition": src,
                     "events": r["events"], "n": r["n"], "n_band0": r["band_n"][0], "n_band1": r["band_n"][1],
                     "n_band2": r["band_n"][2], "n_band3": r["band_n"][3],
                     **{f"hr_{lb}": r[lb]["hr"] for lb in ("1-5%", "5-10%", ">10%")},
                     **{f"lo_{lb}": r[lb]["lo"] for lb in ("1-5%", "5-10%", ">10%")},
                     **{f"hi_{lb}": r[lb]["hi"] for lb in ("1-5%", "5-10%", ">10%")},
                     **{f"p_{lb}": r[lb]["p"] for lb in ("1-5%", "5-10%", ">10%")},
                     "trend_hr": r["trend_hr"], "trend_p": r["trend_p"]})
bands = pd.DataFrame(rows); bands.to_csv(f"{HERE}/bands_recording_vs_sleep_only_n6803.csv", index=False)
say(f"  {bands.key.nunique()} outcomes fitted both ways")

# ------------------------------------------------------------------ 2. held-out concordance gain (build_ranking_v3 procedure)
say("2. held-out concordance gain over age and sex, four exposures, same people")
d = m.copy()
EXPO = {"recording_T90": "t90", "sleep_only_T90": "all_sleep_t90", "wake_only_T90": "all_wake_t90", "total_sleep_time": "TST_min"}
for name, col in EXPO.items():
    z = pd.Series(np.nan, index=d.index)
    for s, idx in d.groupby("site_id").groups.items():
        v = d.loc[idx, col]; ok = v.notna()
        if ok.sum() < 20: continue
        r = stats.rankdata(v[ok], method="average")
        z.loc[v[ok].index] = stats.norm.ppf((r - 0.375) / (ok.sum() + 0.25))
    d[name + "__z"] = z
OUT = [k for k, _ in CONDS if f"{k}_incident" in d.columns and int(d[f"{k}_incident"].sum()) >= 60]
crec = []
for out in OUT:
    yc, ec = f"{out}_years", f"{out}_incident"
    g = d[(d[f"{out}_prevalent"] == 0) & (d[yc] > 0) & d[yc].notna()]
    side = {s: g[g.site_id == s] for s in ("I0002", "I0006")}
    res_by = {}
    for feat in [None] + list(EXPO):
        cols = ADJ4 + ([feat + "__z"] if feat else [])
        res = []
        for tr, te in (("I0002", "I0006"), ("I0006", "I0002")):
            A = side[tr][cols + [yc, ec]].dropna(); B = side[te][cols + [yc, ec]].dropna()
            if A[ec].sum() < 30 or B[ec].sum() < 30: res.append(np.nan); continue
            try:
                c = CoxPHFitter(penalizer=0.01).fit(A, duration_col=yc, event_col=ec)
                res.append(concordance_index(B[yc], -c.predict_partial_hazard(B), B[ec]))
            except Exception: res.append(np.nan)
        res_by[feat or "__BASE__"] = float(np.nanmean(res)) if not all(np.isnan(res)) else np.nan
    if np.isnan(res_by["__BASE__"]): continue
    crec.append({"outcome": label(out), "key": out, "negative_control": out in NEGATIVE_CONTROLS,
                 "events": int(g[ec].sum()), "baseline_C": res_by["__BASE__"],
                 **{f"gain_{n}": res_by[n] - res_by["__BASE__"] for n in EXPO}})
dc = pd.DataFrame(crec); dc.to_csv(f"{HERE}/dC_recording_vs_sleep_only_n6803.csv", index=False)
est = dc.dropna(subset=[f"gain_{n}" for n in EXPO])
summ = {n: {"mean_gain": float(est[f"gain_{n}"].mean()), "n_positive": int((est[f"gain_{n}"] > 0).sum()),
            "n_outcomes": int(len(est))} for n in EXPO}
say("  mean gain over age+sex, same outcomes, same people:")
for n, s in summ.items(): say(f"    {n:18s} {s['mean_gain']:+.4f}  positive in {s['n_positive']}/{s['n_outcomes']}")

# ------------------------------------------------------------------ report
KEY = [("hf", "heart failure"), ("resp_failure", "respiratory failure"), ("diabetes", "type 2 diabetes"),
       ("cirrhosis", "cirrhosis"), ("mi", "myocardial infarction"), ("copd2", "COPD"), ("death", "death")]
lines = ["# Sleep-only T90: do the categories and the concordance gain hold up?", "",
         f"Run {time.strftime('%Y-%m-%d %H:%M')}. Subset: {len(m):,} participants with stage-resolved oximetry "
         f"(the Fig. 4a wake-vs-sleep set), hospitals {sites}. Recording-based and sleep-only T90 compared in the same people.",
         f"Positive control: published heart-failure bands reproduced from the full cohort "
         f"({pc['1-5%']['hr']:.3f} / {pc['5-10%']['hr']:.3f} / {pc['>10%']['hr']:.3f}, 1,187 events).", "",
         f"Medians: recording T90 {m.t90.median():.2f} %, sleep-only T90 {m.all_sleep_t90.median():.2f} %, wake T90 {m.all_wake_t90.median():.2f} %. "
         f"Spearman recording vs sleep-only {stats.spearmanr(m.t90, m.all_sleep_t90)[0]:.3f}. Same band under both definitions: {agree*100:.1f} % of people.", "",
         "Band cross-tab (rows recording-based, columns sleep-only):", "", xt.to_markdown(), "",
         "## 1. Graded hazard ratios versus the <=1 % group (same people, same model)", "",
         "| Outcome | events | definition | 1-5 % | 5-10 % | >10 % | trend P |", "|---|---|---|---|---|---|---|"]
for k, lab in KEY:
    for src in ("recording", "sleep_only"):
        r = bands[(bands.key == k) & (bands.definition == src)]
        if r.empty: continue
        r = r.iloc[0]
        lines.append(f"| {lab} | {r.events} | {src} | {r['hr_1-5%']:.2f} ({r['lo_1-5%']:.2f}-{r['hi_1-5%']:.2f}) | "
                     f"{r['hr_5-10%']:.2f} ({r['lo_5-10%']:.2f}-{r['hi_5-10%']:.2f}) | {r['hr_>10%']:.2f} ({r['lo_>10%']:.2f}-{r['hi_>10%']:.2f}) | {r.trend_p:.1e} |")
# how many outcomes keep a significant top band under each definition
top = bands.pivot(index="key", columns="definition", values="p_>10%")
hr_top = bands.pivot(index="key", columns="definition", values="hr_>10%")
nc = bands.drop_duplicates("key").set_index("key").negative_control
real = top.index[~nc.loc[top.index]]
lines += ["", f"Across the {len(real)} non-control outcomes fitted both ways: top band (>10 %) P < 0.05 for "
          f"{int((top.loc[real,'recording']<0.05).sum())} outcomes recording-based and {int((top.loc[real,'sleep_only']<0.05).sum())} sleep-only; "
          f"median top-band HR {hr_top.loc[real,'recording'].median():.2f} recording-based vs {hr_top.loc[real,'sleep_only'].median():.2f} sleep-only. "
          f"Negative controls with top band P < 0.05: {int((top.loc[top.index[nc.loc[top.index]],'recording']<0.05).sum())} recording-based, "
          f"{int((top.loc[top.index[nc.loc[top.index]],'sleep_only']<0.05).sum())} sleep-only.", "",
          "## 2. Held-out concordance gain over age and sex (ranking procedure, same people)", "",
          "| exposure | mean gain | outcomes with a positive gain |", "|---|---|---|"]
for n, s in summ.items(): lines.append(f"| {n} | {s['mean_gain']:+.4f} | {s['n_positive']} of {s['n_outcomes']} |")
lines += ["", "Per-disease gains for the cover-letter diseases:", "",
          "| disease | events | recording T90 | sleep-only T90 | wake-only T90 | total sleep time |", "|---|---|---|---|---|---|"]
for k, lab in [("resp_failure", "respiratory failure"), ("cirrhosis", "cirrhosis"), ("diabetes", "type 2 diabetes"), ("mi", "myocardial infarction"), ("hf", "heart failure")]:
    r = dc[dc.key == k]
    if r.empty: lines.append(f"| {lab} | not estimable in the subset | | | | |"); continue
    r = r.iloc[0]
    lines.append(f"| {lab} | {r.events} | {r.gain_recording_T90:+.4f} | {r.gain_sleep_only_T90:+.4f} | {r.gain_wake_only_T90:+.4f} | {r.gain_total_sleep_time:+.4f} |")
lines += ["", "Files: bands_recording_vs_sleep_only_n6803.csv, dC_recording_vs_sleep_only_n6803.csv, RESULT.json.",
          "Not a re-ranking of 197 measurements (that needs the whole feature master in this subset); this is the four exposures that matter for the claim."]
open(f"{HERE}/SUMMARY.md", "w").write("\n".join(lines))
json.dump({"n_subset": int(len(m)), "sites": sites, "positive_control_hf_bands": {lab: pc[lab]["hr"] for lab in ("1-5%", "5-10%", ">10%")},
           "band_agreement": agree, "spearman_recording_sleep": float(stats.spearmanr(m.t90, m.all_sleep_t90)[0]),
           "dC_summary": summ}, open(f"{HERE}/RESULT.json", "w"), indent=1)
say("DONE ->", HERE)
