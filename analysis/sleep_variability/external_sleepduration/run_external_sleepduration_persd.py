"""
Colleague request 4, per-1-SD addendum. Total sleep time entered as the exposure in the
eTable 8 external-validation models, on the SAME transformation and scaling the external
T90 models use: a rank-based inverse-normal transform, so the exposure carries unit
variance and the hazard ratio is per 1 SD.

WHY THIS FILE EXISTS. run_external_sleepduration.py (same folder, 2026-08-19) already ran
sleep duration continuously per 1 hour and categorically. The colleague asked for the
per-1-SD scale, "the same transformation/scaling used for the existing external models".
In numbers/freeze_05_external.py that scaling is rint(), a rank-based inverse-normal
transform applied once to the whole SHHS analysis frame and within site in MrOS. This
script applies exactly that function to total sleep time and refits.

  SHHS   tst_z = rint(slpprdp)                     whole-cohort ranks, as for pctsa90h
  MrOS   tst_z = rint(POSLPRDP) within SITE        as for POPCSA90

Higher tst_z means LONGER sleep. A hazard ratio below 1 therefore means longer sleep sits
with lower risk. The mirror estimate per 1 SD SHORTER sleep is 1/HR, and is printed too.

MODELS, 12 cohort-outcome pairs (7 SHHS, 5 MrOS), covariates, exclusions, follow-up,
censoring and event definitions all copied from freeze_05_external.py, unchanged:
  S1  tst_z + covariates                  sleep duration alone, per 1 SD
  S2  tst_z + t90_z + covariates          mutually adjusted, T90 exactly as in eTable 8
Covariates: natural cubic spline on age, 4 df, first basis column dropped, plus sex in
SHHS. MrOS is stratified by clinical site. Fits are unpenalized maximum partial
likelihood with the frozen script's Newton step-size scan.

POSITIVE CONTROL FIRST. Before any sleep-duration model runs, the harness refits the 48
frozen eTable 8 estimates (T90 alone, AHI alone, and the mutually adjusted pair, for all
12 pairs) and compares them to numbers/external.json within 5e-4 on HR and both limits
(v8.3, both sides full precision; was equality at 3 decimals) and exactly on event counts. The run aborts if any of the 48 fails.

MULTIPLICITY. Hand-rolled Benjamini-Hochberg within each exposure family across the 12
cohort-outcome pairs. Four families of 12: sleep duration alone, sleep duration adjusted
for T90, T90 adjusted for sleep duration, and, for the like-for-like count, T90 alone
(its p values read from the frozen external.json, corrected here for comparability).

PER-HOUR BRIDGE. The per-hour fits already in results.csv are joined on, together with
the hours of sleep that one SD of tst_z spans, so the two scalings can be read against
each other rather than guessed at.

Reads only raw cohort files and the frozen numbers/external.json plus this folder's
results.csv. Writes only into this folder:
  results_persd.csv          positive control + every per-SD fit, with q values
  summary_persd.json         spec, counts, scaling bridge, FDR survivors
  eTable_persd_draft.md      supplement-style draft table, two panels
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import warnings
warnings.filterwarnings("ignore")
import json
import os
import subprocess
import sys
import numpy as np
import pandas as pd
from scipy import stats
from patsy import dmatrix
from lifelines import CoxPHFitter

HERE = os.path.dirname(os.path.abspath(__file__))
NUM = paths.NUMBERS_DIR
SH = paths.SHHS_DIR
STEPS = [None, 0.05, 0.1, 0.25, 0.5, 0.95]
EVMIN = {"shhs": 30, "mros": 25}
LABEL = {"shhs": "Sleep Heart Health Study", "mros": "Osteoporotic Fractures in Men Study"}


def mros_path(name):
    """Resolve an MrOS SAS file, identical resolver to run_external_sleepduration.py."""
    cands = [os.path.join(HERE, "data", name),
             f"{paths.MROS_DIR}/{name.split('.')[0].upper()}/{name}",
             os.path.expanduser(f"{paths.MROS_DIR}/{name}")]
    for c in cands:
        if os.path.exists(c):
            return c
    z = f"{paths.MROS_DIR}/{name.split('.')[0].upper()}.zip"
    if os.path.exists(z):
        os.makedirs(os.path.join(HERE, "data"), exist_ok=True)
        subprocess.run(["unzip", "-o", "-q", z, name, "-d", os.path.join(HERE, "data")],
                       check=True)
        return os.path.join(HERE, "data", name)
    raise FileNotFoundError(name)


# ---- fitting machinery, copied verbatim from freeze_05_external.py ----------------------
def fit_ml(f, strata=None):
    """Unpenalized Cox fit. Scan Newton step sizes, keep the maximum-likelihood solution."""
    best = None
    for st in STEPS:
        kw = {} if st is None else {"fit_options": {"step_size": st}}
        try:
            c = CoxPHFitter(penalizer=0.0).fit(f, "T", "E", strata=strata, **kw)
        except Exception:
            continue
        ll = float(c.log_likelihood_)
        if not np.isfinite(ll):
            continue
        if best is None or ll > best[0] + 1e-9:
            best = (ll, c, st)
    if best is None:
        raise RuntimeError("no Newton step size converged")
    return best[1]


def rint(v):
    """The frozen rank-based inverse-normal transform. Unit variance, so 1 unit = 1 SD."""
    z = pd.Series(np.nan, index=v.index)
    ok = v.notna()
    r = stats.rankdata(v[ok], method="average")
    z.loc[v[ok].index] = stats.norm.ppf((r - 0.375) / (ok.sum() + 0.25))
    return z


def spline(df, col, n=4):
    b = dmatrix(f"cr(a, df={n}) - 1", {"a": df[col].values}, return_type="dataframe")
    b.columns = [f"agesp{i}" for i in range(b.shape[1])]
    b.index = df.index
    b = b.iloc[:, 1:]
    for c in b.columns:
        df[c] = b[c]
    return list(b.columns)


def est(c, x, f):
    lo, hi = np.exp(c.confidence_intervals_.loc[x])
    return {"n": int(len(f)), "events": int(f.E.sum()),
            "hr": float(np.exp(c.params_[x])), "lo": float(lo), "hi": float(hi),
            "p": float(c.summary.loc[x, "p"])}


def bh(pvals):
    """Hand-rolled Benjamini-Hochberg q values."""
    p = np.asarray(pvals, dtype=float)
    m = len(p)
    order = np.argsort(p)
    q = np.empty(m)
    prev = 1.0
    for rank_from_top, idx in enumerate(order[::-1]):
        k = m - rank_from_top
        prev = min(prev, p[idx] * m / k)
        q[idx] = prev
    return q


# ============================================================ load SHHS
d = pd.read_csv(f"{SH}/shhs1-dataset-0.21.0.csv", low_memory=False)
cv = pd.read_csv(f"{SH}/shhs-cvd-summary-dataset-0.21.0.csv", low_memory=False)
cv = cv.drop(columns=[c for c in ["gender", "race", "age_s1"] if c in cv.columns])
s = d[["nsrrid", "pctsa90h", "ahi_a0h3a", "ahi_a0h4a", "age_s1", "gender",
       "bmi_s1", "smokstat_s1", "diasbp", "chol", "hdl", "htnderv_s1",
       "parrptdiab", "slpprdp"]].merge(cv, on="nsrrid", how="inner")
s["male"] = (s.gender == 1).astype(int)
ADJ = spline(s, "age_s1") + ["male"]
s["t90_z"] = rint(s.pctsa90h)
s["ahi_z"] = rint(s.ahi_a0h3a)
s["sleep_h"] = s.slpprdp / 60.0
s["tst_z"] = rint(s.slpprdp)                       # the eTable 8 scaling, applied to TST

SOUT = [("chf", "chf_date", "prev_chf", "Heart failure"),
        ("mi", "mi_date", "prev_mi", "Myocardial infarction"),
        ("stroke", "stk_date", "prev_stk", "Stroke"),
        ("any_cvd", None, None, "Cardiovascular composite"),
        ("any_chd", None, None, "Coronary heart disease"),
        ("cvd_death", "cvd_dthdt", None, "Cardiovascular death"),
        ("afibincident", None, "afibprevalent", "Atrial fibrillation")]


def sframe(k, dc, pc):
    g = s.copy()
    if pc and pc in g.columns:
        g = g[g[pc].fillna(0) == 0]
    ev = (g[k].fillna(0) > 0).astype(int).values
    t = g.censdate.astype(float).values
    if dc and dc in g.columns:
        t = np.where(ev == 1, g[dc].astype(float).values, t)
    g = g.assign(T=t / 365.25, E=ev)
    return g[(g["T"] > 0) & g["T"].notna()]


# ============================================================ load MrOS
pos = pd.read_sas(mros_path("POSFEB23.SAS7BDAT"), encoding="latin-1")
ef = pd.read_sas(mros_path("effeb24.sas7bdat"), encoding="latin-1")
vs = pd.read_sas(mros_path("vsfeb24.sas7bdat"), encoding="latin-1")
m = (pos[["ID", "POPCSA90", "PORDI4P", "POSLPRDP"]]
     .merge(ef[["ID", "SITE", "DADEAD", "DACARDIO", "DACANCER", "DASTROKE", "DAPULMON",
                "FUCDTIME", "FUVSDT"]], on="ID")
     .merge(vs[["ID", "VSAGE1", "HWBMI"]], on="ID", how="left"))
for c in m.columns:
    if c not in ("ID", "SITE"):
        m[c] = pd.to_numeric(m[c], errors="coerce")
m = m[m.FUVSDT.notna() & (m.FUVSDT > 0) & m.VSAGE1.notna() & m.POPCSA90.notna()].copy()
m["SITE"] = m.SITE.astype(str)
MADJ = spline(m, "VSAGE1")


def within_site(col):
    z = pd.Series(np.nan, index=m.index)
    for _site, idx in m.groupby("SITE").groups.items():
        z.loc[idx] = rint(m.loc[idx, col])
    return z


m["t90_z"] = within_site("POPCSA90")
m["ahi_z"] = within_site("PORDI4P")
m["sleep_h"] = m.POSLPRDP / 60.0
m["tst_z"] = within_site("POSLPRDP")               # within site, as POPCSA90 is

MOUT = [("DADEAD", "Death from any cause"), ("DACARDIO", "Cardiovascular death"),
        ("DAPULMON", "Pulmonary death"), ("DASTROKE", "Stroke death"),
        ("DACANCER", "Cancer death")]


def mframe(k, cols):
    f = pd.DataFrame({"T": (m.FUVSDT / 365.25).values,
                      "E": m[k].fillna(0).astype(int).values,
                      "site": m.SITE.values,
                      **{c: m[c].values for c in cols + MADJ}})
    return f.dropna()


# v8.3 (2026-09-16): the cohort sizes are read from the regenerated numbers/external.json (step 260, sidecar), not typed
with open(f"{NUM}/external.json") as fh:
    FROZEN = json.load(fh)
assert len(s) == FROZEN["shhs"]["n"], f"SHHS n drifted: {len(s)} vs frozen {FROZEN['shhs']['n']}"
assert len(m) == FROZEN["mros"]["n"], f"MrOS n drifted: {len(m)} vs frozen {FROZEN['mros']['n']}"
assert int(s.slpprdp.isna().sum()) == 0, "SHHS slpprdp has missing values"
assert int(m.POSLPRDP.isna().sum()) == 0, "MrOS POSLPRDP has missing values"
# the transform must be the frozen one, i.e. unit variance to 2 decimals
assert abs(float(s.tst_z.std(ddof=0)) - 1.0) < 0.02, float(s.tst_z.std(ddof=0))
assert abs(float(m.tst_z.std(ddof=0)) - 1.0) < 0.02, float(m.tst_z.std(ddof=0))

# ============================================================ positive control
# FROZEN (numbers/external.json) is loaded above, before the cohort-size asserts (v8.3)

control_rows = []
control_fail = 0


def check(cohort, lab, key, e):
    global control_fail
    fr = FROZEN[cohort]["outcomes"][lab].get(key)
    if fr is None:
        return
    # v8.3: reference now full precision, compare within 5e-4 (2026-09-16)
    ok = (abs(e["hr"] - fr["hr"]) <= 5e-4 and abs(e["lo"] - fr["lo"]) <= 5e-4
          and abs(e["hi"] - fr["hi"]) <= 5e-4 and e["events"] == fr["events"])
    if not ok:
        control_fail += 1
    control_rows.append({"cohort": cohort, "outcome": lab, "estimate": key,
                         "frozen_hr": fr["hr"], "refit_hr": e["hr"],  # v8.3 full precision (2026-09-16)
                         "frozen_lo": fr["lo"], "refit_lo": e["lo"],
                         "frozen_hi": fr["hi"], "refit_hi": e["hi"],
                         "frozen_events": fr["events"], "refit_events": e["events"],
                         "match": bool(ok)})


print("== POSITIVE CONTROL: refitting the 48 frozen eTable 8 estimates ==")
sfr = {}
for k, dc, pc, lab in SOUT:
    g = sframe(k, dc, pc)
    sfr[lab] = g
    for x in ["t90_z", "ahi_z"]:
        f = g[["T", "E", x] + ADJ].dropna()
        if f.E.sum() < EVMIN["shhs"]:
            continue
        c = fit_ml(f)
        check("shhs", lab, x.replace("_z", ""), est(c, x, f))
    f = g[["T", "E", "t90_z", "ahi_z"] + ADJ].dropna()
    if f.E.sum() >= EVMIN["shhs"]:
        c = fit_ml(f)
        check("shhs", lab, "t90_adj_ahi", est(c, "t90_z", f))
        check("shhs", lab, "ahi_adj_t90", est(c, "ahi_z", f))

for k, lab in MOUT:
    for x in ["t90_z", "ahi_z"]:
        f = mframe(k, [x])
        if f.E.sum() < EVMIN["mros"]:
            continue
        c = fit_ml(f, ["site"])
        check("mros", lab, x.replace("_z", ""), est(c, x, f))
    f = mframe(k, ["t90_z", "ahi_z"])
    if f.E.sum() >= EVMIN["mros"]:
        c = fit_ml(f, ["site"])
        check("mros", lab, "t90_adj_ahi", est(c, "t90_z", f))
        check("mros", lab, "ahi_adj_t90", est(c, "ahi_z", f))

cdf = pd.DataFrame(control_rows)
n_ok = int(cdf.match.sum())
print(f"  {n_ok}/{len(cdf)} frozen estimates reproduced within 5e-4 on HR, lo, hi "
      f"(full precision on both sides, v8.3), exact events")
if control_fail:
    print(cdf[~cdf.match].to_string(index=False))
    sys.exit("POSITIVE CONTROL FAILED, aborting before any sleep-duration model")

# ============================================================ per-SD sleep-duration models
rows = []


def add(cohort, lab, model, term, e):
    rows.append({"block": "sleep_duration_persd", "cohort": cohort,
                 "cohort_label": LABEL[cohort], "outcome": lab, "model": model,
                 "term": term, "n": e["n"], "events": e["events"],
                 "hr": e["hr"], "lo": e["lo"], "hi": e["hi"], "p": e["p"]})


def run_pair(cohort, lab, frame_fn, strata):
    evmin = EVMIN[cohort]
    f = frame_fn(["tst_z"])
    if f.E.sum() >= evmin:
        c = fit_ml(f, strata)
        add(cohort, lab, "S1_persd_alone", "tst_z", est(c, "tst_z", f))
    f = frame_fn(["tst_z", "t90_z"])
    if f.E.sum() >= evmin:
        c = fit_ml(f, strata)
        add(cohort, lab, "S2_persd_adjT90", "tst_z", est(c, "tst_z", f))
        add(cohort, lab, "S2_persd_adjT90", "t90_z", est(c, "t90_z", f))


for k, dc, pc, lab in SOUT:
    g = sfr[lab]
    run_pair("shhs", lab, lambda cols, g=g: g[["T", "E"] + cols + ADJ].dropna(), None)
for k, lab in MOUT:
    run_pair("mros", lab, lambda cols, k=k: mframe(k, cols), ["site"])

# T90 alone, from the frozen file, so the like-for-like count uses the same 12-test family
PAIRS = [("shhs", lab) for _, _, _, lab in SOUT] + [("mros", lab) for _, lab in MOUT]
for coh, lab in PAIRS:
    fr = FROZEN[coh]["outcomes"][lab]["t90"]
    rows.append({"block": "frozen_etable8", "cohort": coh, "cohort_label": LABEL[coh],
                 "outcome": lab, "model": "T0_t90_alone_frozen", "term": "t90_z",
                 "n": fr["n"], "events": fr["events"], "hr": fr["hr"], "lo": fr["lo"],
                 "hi": fr["hi"], "p": fr["p"]})

res = pd.DataFrame(rows)

FAMY = {("S1_persd_alone", "tst_z"): "tst_persd_alone",
        ("S2_persd_adjT90", "tst_z"): "tst_persd_adjT90",
        ("S2_persd_adjT90", "t90_z"): "t90_persd_adjTST",
        ("T0_t90_alone_frozen", "t90_z"): "t90_persd_alone_frozen"}
res["family"] = [FAMY[(r.model, r.term)] for r in res.itertuples()]
res["q"] = np.nan
for fam, idx in res.groupby("family").groups.items():
    assert len(idx) == 12, f"family {fam} has {len(idx)} members, expected 12"
    res.loc[idx, "q"] = bh(res.loc[idx, "p"].values)
res["hr_per_sd_shorter"] = 1.0 / res.hr
res["sig_raw"] = res.p < 0.05
res["sig_fdr"] = res.q < 0.05

# ============================================================ per-hour bridge
prev = pd.read_csv(os.path.join(HERE, "results.csv"))
prev = prev[prev.block == "sleep_duration"]
per_hour = {}
for mdl, key in [("M1_cont_alone", "alone"), ("M2_cont_adjT90", "adjT90")]:
    sub = prev[(prev.model == mdl) & (prev.term == "sleep_h")]
    for r in sub.itertuples():
        per_hour[(r.cohort, r.outcome, key)] = {"hr": float(r.hr), "lo": float(r.lo),
                                                "hi": float(r.hi), "p": float(r.p),
                                                "q": float(r.q)}
assert len(per_hour) == 24, len(per_hour)

# hours of sleep spanned by one SD of the transformed exposure, by ordinary least squares
# of hours on tst_z. This is the number that converts one scale into the other.
HRS_PER_SD = {}
for coh, frame in (("shhs", s), ("mros", m)):
    x = frame.tst_z.values
    y = frame.sleep_h.values
    slope = float(np.polyfit(x, y, 1)[0])
    HRS_PER_SD[coh] = {"hours_per_1sd_ols": round(slope, 4),
                       "raw_sd_hours": round(float(frame.sleep_h.std(ddof=1)), 4),
                       "mean_hours": round(float(frame.sleep_h.mean()), 4),
                       "median_hours": round(float(frame.sleep_h.median()), 4),
                       "iqr_hours": [round(float(frame.sleep_h.quantile(.25)), 4),
                                     round(float(frame.sleep_h.quantile(.75)), 4)]}

bridge = []
for coh, lab in PAIRS:
    for key, mdl in (("alone", "S1_persd_alone"), ("adjT90", "S2_persd_adjT90")):
        ph = per_hour[(coh, lab, key)]
        sd = res[(res.cohort == coh) & (res.outcome == lab) & (res.model == mdl)
                 & (res.term == "tst_z")].iloc[0]
        h = HRS_PER_SD[coh]["hours_per_1sd_ols"]
        bridge.append({"cohort": coh, "outcome": lab, "model": key,
                       "per_hour_hr": round(ph["hr"], 4), "per_hour_p": ph["p"],
                       "per_hour_q": ph["q"],
                       "hours_per_1sd": h,
                       "expected_per_sd_from_per_hour": round(ph["hr"] ** h, 4),
                       "observed_per_sd_hr": round(float(sd.hr), 4),
                       "observed_per_sd_p": float(sd.p), "observed_per_sd_q": float(sd.q),
                       "abs_diff_hr": round(abs(ph["hr"] ** h - float(sd.hr)), 4)})
bdf = pd.DataFrame(bridge)

# ============================================================ summary counts
def counts(fam):
    sub = res[res.family == fam]
    return {"n_outcomes": int(len(sub)), "raw_p_lt_05": int(sub.sig_raw.sum()),
            "fdr_q_lt_05": int(sub.sig_fdr.sum()),
            "ci_excludes_1": int((~((sub.lo <= 1.0) & (1.0 <= sub.hi))).sum())}


COUNTS = {f: counts(f) for f in ["t90_persd_alone_frozen", "t90_persd_adjTST",
                                 "tst_persd_alone", "tst_persd_adjTST"]
          if f in set(res.family)}
COUNTS["tst_persd_adjT90"] = counts("tst_persd_adjT90")
COUNTS.pop("tst_persd_adjTST", None)

print("\n== PER 1 SD, SLEEP DURATION (rank-inverse-normal, higher = longer sleep) ==")
for r in res[res.family.str.startswith("tst_")].itertuples():
    print(f"  {r.cohort:4s} {r.outcome:26s} {r.model:16s} HR {r.hr:6.3f} "
          f"({r.lo:6.3f}-{r.hi:6.3f}) p={r.p:.4f} q={r.q:.3f}")
print("\n== T90 PER 1 SD INSIDE THE MUTUALLY ADJUSTED MODEL ==")
for r in res[res.family == "t90_persd_adjTST"].itertuples():
    print(f"  {r.cohort:4s} {r.outcome:26s} HR {r.hr:6.3f} ({r.lo:6.3f}-{r.hi:6.3f}) "
          f"p={r.p:.4f} q={r.q:.3f}")
print("\n== SUMMARY COUNTS out of 12 ==")
for fam, c in COUNTS.items():
    print(f"  {fam:26s} raw p<0.05 {c['raw_p_lt_05']:2d}   FDR q<0.05 {c['fdr_q_lt_05']:2d}"
          f"   CI excludes 1 {c['ci_excludes_1']:2d}")

# T90 stability against the frozen eTable 8 value
shift = []
for coh, lab in PAIRS:
    fr = FROZEN[coh]["outcomes"][lab]["t90"]["hr"]
    o = res[(res.cohort == coh) & (res.outcome == lab)
            & (res.family == "t90_persd_adjTST")].iloc[0]
    shift.append(abs(float(o.hr) - fr))
MAXSHIFT = float(max(shift))
print(f"\nT90 stability: max |HR shift| against the frozen eTable 8 T90-alone value "
      f"across the 12 mutually adjusted models = {MAXSHIFT:.4f}")

# ============================================================ write outputs
full = pd.concat([cdf.assign(block="positive_control"), res,
                  bdf.assign(block="per_hour_bridge")], ignore_index=True)
full.to_csv(os.path.join(HERE, "results_persd.csv"), index=False)

summary = {
    "positive_control": {"n_estimates": int(len(cdf)), "n_match": n_ok,
                         "all_match": bool(control_fail == 0)},
    "spec": {"source": "numbers/freeze_05_external.py model spec, exposure swapped and "
                       "rescaled to the frozen rank-inverse-normal transform",
             "scaling": "rank-based inverse-normal, unit variance, hazard ratio per 1 SD",
             "direction": "higher tst_z = longer sleep, HR below 1 = longer sleep, "
                          "lower risk",
             "shhs_exposure": "rint(slpprdp), whole-cohort ranks",
             "mros_exposure": "rint(POSLPRDP) within SITE",
             "t90_terms": {"shhs": "rint(pctsa90h)",
                           "mros": "rint(POPCSA90) within SITE"},
             "fit": "Cox by maximum partial likelihood, penalizer 0, Newton step sizes "
                    "scanned, highest-likelihood solution kept",
             "age_basis": "natural cubic spline, 4 df, first basis column dropped",
             "shhs_covariates": "age spline + male",
             "mros_covariates": "age spline, stratified by site",
             "bh": "within each exposure family across the 12 cohort-outcome pairs, "
                   "4 families of 12"},
    "cohort_n": {"shhs": int(len(s)), "mros": int(len(m))},
    "sleep_duration_distribution": HRS_PER_SD,
    "summary_counts_out_of_12": COUNTS,
    "t90_max_hr_shift_vs_frozen": round(MAXSHIFT, 4),
    "fdr_survivors": res[res.q < 0.05][
        ["cohort", "outcome", "family", "hr", "lo", "hi", "p", "q"]].to_dict("records"),
    "per_hour_bridge_max_abs_diff": float(bdf.abs_diff_hr.max()),
}
with open(os.path.join(HERE, "summary_persd.json"), "w") as fh:
    json.dump(summary, fh, indent=2)


# ============================================================ supplement-style draft table
def cistr(r):
    return f"{r['hr']:.3f} ({r['lo']:.3f}-{r['hi']:.3f})"


def pstr(p):
    return "<.001" if p < 0.001 else f"{p:.3f}"


def pick(coh, lab, fam):
    return res[(res.cohort == coh) & (res.outcome == lab)
               & (res.family == fam)].iloc[0].to_dict()


def md_table(df):
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |",
             "|" + "|".join("---" for _ in cols) + "|"]
    for _, row in df.iterrows():
        lines.append("| " + " | ".join(str(row[c]) for c in cols) + " |")
    return "\n".join(lines)


pa = []
for coh, lab in PAIRS:
    a = pick(coh, lab, "tst_persd_alone")
    b = pick(coh, lab, "tst_persd_adjT90")
    pa.append({"Cohort": LABEL[coh], "Outcome": lab, "Events": int(a["events"]),
               "Sleep duration alone (95% CI)": cistr(a), "P value": pstr(a["p"]),
               "q value": f"{a['q']:.3f}",
               "Sleep duration, oxygen-adjusted (95% CI)": cistr(b),
               "Adjusted P value": pstr(b["p"]), "Adjusted q value": f"{b['q']:.3f}"})

pb = []
for coh, lab in PAIRS:
    fr = pick(coh, lab, "t90_persd_alone_frozen")
    o = pick(coh, lab, "t90_persd_adjTST")
    pb.append({"Cohort": LABEL[coh], "Outcome": lab, "Events": int(o["events"]),
               "Oxygen alone, frozen eTable 8 (95% CI)": cistr(fr),
               "P value": pstr(fr["p"]), "q value": f"{fr['q']:.3f}",
               "Oxygen, sleep-duration-adjusted (95% CI)": cistr(o),
               "Adjusted P value": pstr(o["p"]), "Adjusted q value": f"{o['q']:.3f}"})

C = COUNTS
# the caption's claim, asserted: the same outcomes survive, not merely the same number
SET_ALONE = {(r.cohort, r.outcome) for r in
             res[(res.family == "t90_persd_alone_frozen") & res.sig_fdr].itertuples()}
SET_ADJ = {(r.cohort, r.outcome) for r in
           res[(res.family == "t90_persd_adjTST") & res.sig_fdr].itertuples()}
assert SET_ALONE == SET_ADJ and len(SET_ALONE) == 7, (SET_ALONE, SET_ADJ)
TST_RAW_ALONE = {(r.cohort, r.outcome) for r in
                 res[(res.family == "tst_persd_alone") & res.sig_raw].itertuples()}
TST_RAW_ADJ = {(r.cohort, r.outcome) for r in
               res[(res.family == "tst_persd_adjT90") & res.sig_raw].itertuples()}
assert TST_RAW_ALONE == TST_RAW_ADJ == {("mros", "Death from any cause")}, TST_RAW_ALONE
CAPTION = (
    "**eTable 18B.** Total sleep time and nocturnal oxygen entered on the same per-1-SD "
    "scale in the external-validation Cox models of eTable 8, alone and mutually adjusted, "
    "in the Sleep Heart Health Study (n=5,802, 7 adjudicated cardiovascular outcomes) and "
    "the Osteoporotic Fractures in Men Study (n=2,911, 5 adjudicated causes of death). "
    "Total sleep time is the total sleep time of the baseline polysomnogram (SHHS variable "
    "slpprdp, MrOS variable POSLPRDP), rank-based inverse-normal transformed within cohort, "
    "and within site in MrOS, which is the identical transformation eTable 8 applies to "
    "nocturnal oxygen. Hazard ratios are therefore per 1 SD, and a hazard ratio below 1 "
    "means that longer sleep sits with lower risk. Oxygen is the percentage of sleep time "
    "below 90% saturation (SHHS pctsa90h, MrOS POPCSA90) on the same scale, exactly as in "
    "eTable 8. All models carry the eTable 8 covariates: a natural cubic spline on age with "
    "4 df and one basis column dropped, plus sex in SHHS, with stratification by clinical "
    "site in MrOS. Fits are by unpenalized maximum partial likelihood. The q values are "
    "Benjamini-Hochberg corrected within each exposure family across the 12 cohort-outcome "
    "pairs. Panel A gives total sleep time, Panel B the oxygen term of the same models "
    "beside its frozen eTable 8 value. Of the 12 cohort-outcome pairs, nocturnal oxygen "
    f"alone is significant in {C['t90_persd_alone_frozen']['fdr_q_lt_05']}, both at P<.05 "
    "and after correction, and the same "
    f"{C['t90_persd_adjTST']['fdr_q_lt_05']} remain significant once total sleep time is in "
    f"the model. Total sleep time is significant in "
    f"{C['tst_persd_alone']['raw_p_lt_05']} pair at P<.05 and in none after correction, "
    "and that one pair is the same with and without oxygen in the model, death from any "
    "cause in the Osteoporotic Fractures in Men Study.\n")

parts = [CAPTION,
         "\n**Panel A. Total sleep time, per 1 SD (higher = longer sleep)**\n",
         md_table(pd.DataFrame(pa)),
         "\n\n**Panel B. Nocturnal oxygen, per 1 SD, in the same models**\n",
         md_table(pd.DataFrame(pb)), ""]
with open(os.path.join(HERE, "eTable_persd_draft.md"), "w") as fh:
    fh.write("\n".join(parts))

print(f"\nwritten -> {HERE}/results_persd.csv, summary_persd.json, eTable_persd_draft.md")
