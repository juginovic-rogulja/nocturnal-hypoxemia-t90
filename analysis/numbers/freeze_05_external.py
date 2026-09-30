"""
Phase 1.5. Freeze the external validation numbers, matched to the primary-cohort method.

Matched: natural cubic spline on age (4 df), rank-based inverse normal exposure, sex where the
cohort has both sexes, site stratification where the cohort has sites.
Not matchable, and stated as limitations: SHHS and MrOS denominate T90 on total sleep time
while the primary cohort uses the whole recording; outcome ascertainment is physician
adjudication rather than diagnostic coding; SHHS publishes no field-centre variable; and
neither cohort recorded any incident non-cardiac condition, so the negative-control panel used
in the primary cohort has no counterpart.

Writes numbers/external.json.

v8.3 (2026-09-16, round 38 phase 1c): every result value is stored as the float the fit returned. The
August file stored round(x, 3) on hazard ratios and limits, round(x, 4) on the clinical-model concordance
and round(x, 2) on the follow-up years, so two printed digits sat on undecidable ties. No model, cohort
rule, seed or ordering changed; the superseded block is carried over unchanged.

Three corrections applied 2026-08-08 (decision B5):
  1  MrOS follow-up now runs on FUVSDT, days from the SLEEP visit, not FUCDTIME, which counts
     from the parent-cohort baseline roughly three and a half years earlier. The paper reports a
     12.2-year median from the sleep study, so the old clock did not match the stated design.
     The same 2,911 men and the same event counts carry over, only the time axis changes.
  2  The age basis returned all four spline columns, which sum to one and leave the design
     rank-deficient. One column is now dropped.
  3  Every fit was penalized, four at 0.01 and four at 0.05. lifelines scales its penalty by the
     sample size, so those were effective ridges in the tens to hundreds, shrinking every hazard
     ratio toward 1 by an amount that varies with the event count. All fits are now by maximum
     partial likelihood, with several Newton step sizes scanned and the highest-likelihood
     solution kept.
The previous contents of external.json are preserved under the key "superseded_fucdtime_penalized".
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../..")))  # repository root, where paths.py lives
import paths
import warnings
warnings.filterwarnings("ignore")
import json
import os
import numpy as np
import pandas as pd
from scipy import stats
from patsy import dmatrix
from lifelines import CoxPHFitter

OUT = paths.NUMBERS_DIR
SH = paths.SHHS_DIR
MR = os.path.expanduser(paths.MROS_DIR)
STEPS = [None, 0.05, 0.1, 0.25, 0.5, 0.95]
R = {}


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
    z = pd.Series(np.nan, index=v.index)
    ok = v.notna()
    r = stats.rankdata(v[ok], method="average")
    z.loc[v[ok].index] = stats.norm.ppf((r - 0.375) / (ok.sum() + 0.25))
    return z


def spline(df, col, n=4):
    """Natural cubic spline on age. One column is dropped: the n columns sum to one, so the
    full basis is rank-deficient and an unpenalized solver can silently return nothing.

    The basis columns are named agesp*, not age_s*. SHHS carries its own variable called
    age_s1, the man's age at visit 1, and the old naming overwrote it with the second basis
    column, so the clinical risk model at CLIN below was silently adjusting for a spline
    basis function instead of age. Fixed 2026-08-08.
    """
    b = dmatrix(f"cr(a, df={n}) - 1", {"a": df[col].values}, return_type="dataframe")
    b.columns = [f"agesp{i}" for i in range(b.shape[1])]
    b.index = df.index
    b = b.iloc[:, 1:]
    for c in b.columns:
        df[c] = b[c]
    return list(b.columns)


def res(c, x, f):
    lo, hi = np.exp(c.confidence_intervals_.loc[x])
    return {"n": int(len(f)), "events": int(f.E.sum()),
            "hr": float(np.exp(c.params_[x])), "lo": float(lo),  # v8.3 full precision (2026-09-16)
            "hi": float(hi), "p": float(c.summary.loc[x, "p"])}


# ============================================================ SHHS
d = pd.read_csv(f"{SH}/shhs1-dataset-0.21.0.csv", low_memory=False)
cv = pd.read_csv(f"{SH}/shhs-cvd-summary-dataset-0.21.0.csv", low_memory=False)
cv = cv.drop(columns=[c for c in ["gender", "race", "age_s1"] if c in cv.columns])
s = d[["nsrrid", "pctsa90h", "ahi_a0h3a", "ahi_a0h4a", "age_s1", "gender",
       "bmi_s1", "smokstat_s1", "diasbp", "chol", "hdl", "htnderv_s1",
       "parrptdiab"]].merge(cv, on="nsrrid", how="inner")
s["male"] = (s.gender == 1).astype(int)
ADJ = spline(s, "age_s1") + ["male"]
s["t90_z"] = rint(s.pctsa90h)
s["ahi_z"] = rint(s.ahi_a0h3a)

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


R["shhs"] = {"n": int(len(s)), "outcomes": {}}
for k, dc, pc, lab in SOUT:
    g = sframe(k, dc, pc)
    row = {}
    for x in ["t90_z", "ahi_z"]:
        f = g[["T", "E", x] + ADJ].dropna()
        if f.E.sum() < 30:
            continue
        c = fit_ml(f)
        row[x.replace("_z", "")] = res(c, x, f)
    f = g[["T", "E", "t90_z", "ahi_z"] + ADJ].dropna()
    if f.E.sum() >= 30:
        c = fit_ml(f)
        row["t90_adj_ahi"] = res(c, "t90_z", f)
        row["ahi_adj_t90"] = res(c, "ahi_z", f)
    # incremental value over a full clinical risk model
    CLIN = ["age_s1", "male", "bmi_s1", "smokstat_s1", "diasbp", "chol", "hdl",
            "htnderv_s1", "parrptdiab"]
    fc = g[["T", "E", "t90_z"] + CLIN].dropna()
    if fc.E.sum() >= 60:
        from lifelines.utils import concordance_index
        c0 = fit_ml(fc[["T", "E"] + CLIN])
        c1 = fit_ml(fc[["T", "E", "t90_z"] + CLIN])
        b0 = concordance_index(fc["T"], -c0.predict_partial_hazard(fc), fc.E)
        b1 = concordance_index(fc["T"], -c1.predict_partial_hazard(fc), fc.E)
        row["clinical_model_c"] = float(b0)  # v8.3 full precision (2026-09-16)
        row["clinical_model_c_with_t90"] = float(b1)
        row["clinical_model_gain"] = float(b1 - b0)
    # bands
    g2 = g.copy()
    g2["band"] = pd.cut(g2.pctsa90h, [-1, 1, 5, 10, 1e9], labels=[0, 1, 2, 3]).astype(float)
    g2 = g2[g2.band.notna()]
    X = pd.get_dummies(g2.band.astype(int), prefix="b").astype(int)
    for i in range(4):
        if f"b_{i}" not in X:
            X[f"b_{i}"] = 0
    X = X[[f"b_{i}" for i in range(4)]].drop(columns=["b_0"])
    fb = pd.concat([g2[["T", "E"] + ADJ].reset_index(drop=True),
                    X.reset_index(drop=True)], axis=1).dropna()
    if fb.E.sum() >= 60:
        try:
            cb = fit_ml(fb)
            row["bands"] = {}
            for i, nm in enumerate(["1-5%", "5-10%", ">10%"], start=1):
                col = f"b_{i}"
                if col in cb.params_.index:
                    lo, hi = np.exp(cb.confidence_intervals_.loc[col])
                    row["bands"][nm] = {"hr": float(np.exp(cb.params_[col])),  # v8.3 full precision (2026-09-16)
                                        "lo": float(lo), "hi": float(hi),
                                        "p": float(cb.summary.loc[col, "p"])}
        except Exception:
            pass
    R["shhs"]["outcomes"][lab] = row

# ============================================================ MrOS
pos = pd.read_sas(f"{MR}/POSFEB23.SAS7BDAT", encoding="latin-1")
ef = pd.read_sas(f"{MR}/effeb24.sas7bdat", encoding="latin-1")
vs = pd.read_sas(f"{MR}/vsfeb24.sas7bdat", encoding="latin-1")
m = (pos[["ID", "POPCSA90", "PORDI4P"]]
     .merge(ef[["ID", "SITE", "DADEAD", "DACARDIO", "DACANCER", "DASTROKE", "DAPULMON",
                "FUCDTIME", "FUVSDT"]], on="ID")
     .merge(vs[["ID", "VSAGE1", "HWBMI"]], on="ID", how="left"))
for c in m.columns:
    if c not in ("ID", "SITE"):
        m[c] = pd.to_numeric(m[c], errors="coerce")
# FUVSDT is days from the sleep visit. FUCDTIME counts from the parent-cohort baseline, about
# three and a half years earlier, and does not match the follow-up the paper reports.
m = m[m.FUVSDT.notna() & (m.FUVSDT > 0) & m.VSAGE1.notna() & m.POPCSA90.notna()].copy()
m["SITE"] = m.SITE.astype(str)
MADJ = spline(m, "VSAGE1")
z = pd.Series(np.nan, index=m.index)
for site, idx in m.groupby("SITE").groups.items():
    z.loc[idx] = rint(m.loc[idx, "POPCSA90"])
m["t90_z"] = z
z2 = pd.Series(np.nan, index=m.index)
for site, idx in m.groupby("SITE").groups.items():
    z2.loc[idx] = rint(m.loc[idx, "PORDI4P"])
m["ahi_z"] = z2

R["mros"] = {"n": int(len(m)), "sites": int(m.SITE.nunique()),
             "followup_var": "FUVSDT",
             "followup_years_median": float((m.FUVSDT / 365.25).median()),  # v8.3 full precision (2026-09-16): 2-dp store printed at 1 dp
             "followup_years_iqr": [float((m.FUVSDT / 365.25).quantile(.25)),
                                    float((m.FUVSDT / 365.25).quantile(.75))],
             "followup_years_median_superseded_FUCDTIME":
                 float((m.FUCDTIME / 365.25).median()),
             "outcomes": {}}
for k, lab in [("DADEAD", "Death from any cause"), ("DACARDIO", "Cardiovascular death"),
               ("DAPULMON", "Pulmonary death"), ("DASTROKE", "Stroke death"),
               ("DACANCER", "Cancer death")]:
    row = {}
    for x in ["t90_z", "ahi_z"]:
        f = pd.DataFrame({"T": (m.FUVSDT / 365.25).values,
                          "E": m[k].fillna(0).astype(int).values, x: m[x].values,
                          "site": m.SITE.values,
                          **{c: m[c].values for c in MADJ}}).dropna()
        if f.E.sum() < 25:
            continue
        c = fit_ml(f, ["site"])
        row[x.replace("_z", "")] = res(c, x, f)
    f = pd.DataFrame({"T": (m.FUVSDT / 365.25).values,
                      "E": m[k].fillna(0).astype(int).values,
                      "t90_z": m.t90_z.values, "ahi_z": m.ahi_z.values,
                      "site": m.SITE.values,
                      **{c: m[c].values for c in MADJ}}).dropna()
    if f.E.sum() >= 25:
        c = fit_ml(f, ["site"])
        row["t90_adj_ahi"] = res(c, "t90_z", f)
        row["ahi_adj_t90"] = res(c, "ahi_z", f)
    # bands
    mb = m.copy()
    mb["band"] = pd.cut(mb.POPCSA90, [-1, 1, 5, 10, 1e9], labels=[0, 1, 2, 3]).astype(float)
    mb = mb[mb.band.notna()]
    X = pd.get_dummies(mb.band.astype(int), prefix="b").astype(int)
    for i in range(4):
        if f"b_{i}" not in X:
            X[f"b_{i}"] = 0
    X = X[[f"b_{i}" for i in range(4)]].drop(columns=["b_0"])
    fb = pd.concat([pd.DataFrame({"T": (mb.FUVSDT / 365.25).values,
                                  "E": mb[k].fillna(0).astype(int).values,
                                  "site": mb.SITE.values,
                                  **{c: mb[c].values for c in MADJ}}),
                    X.reset_index(drop=True)], axis=1).dropna()
    if fb.E.sum() >= 40:
        try:
            cb = fit_ml(fb, ["site"])
            row["bands"] = {}
            for i, nm in enumerate(["1-5%", "5-10%", ">10%"], start=1):
                col = f"b_{i}"
                if col in cb.params_.index:
                    lo, hi = np.exp(cb.confidence_intervals_.loc[col])
                    row["bands"][nm] = {"hr": float(np.exp(cb.params_[col])),  # v8.3 full precision (2026-09-16)
                                        "lo": float(lo), "hi": float(hi),
                                        "p": float(cb.summary.loc[col, "p"])}
        except Exception:
            pass
    R["mros"]["outcomes"][lab] = row

R["spec"] = {"fit": "Cox by maximum partial likelihood, penalizer 0",
             "newton_step_sizes_scanned": [str(s) for s in STEPS],
             "age_basis": "natural cubic spline, 4 df, one column dropped",
             "mros_followup": "FUVSDT, days from the sleep visit",
             "corrected": "2026-08-08, decision B5"}

SUPKEY = "superseded_fucdtime_penalized"
if os.path.exists(f"{OUT}/external.json"):
    with open(f"{OUT}/external.json") as fh:
        prev = json.load(fh)
    # keep only the first superseded block, so the file does not grow a chain of copies
    R[SUPKEY] = prev.get(SUPKEY, {k: v for k, v in prev.items() if k != SUPKEY})

with open(f"{OUT}/external.json", "w") as f:
    json.dump(R, f, indent=2)

for coh in ("shhs", "mros"):
    print(f"\n=== {coh.upper()}  n={R[coh]['n']:,}")
    for lab, v in R[coh]["outcomes"].items():
        t = v.get("t90"); a = v.get("ahi")
        ta = v.get("t90_adj_ahi"); at = v.get("ahi_adj_t90")
        if not t:
            continue
        print(f"  {lab:<26} ev={t['events']:<5} T90 {t['hr']:.2f} ({t['lo']:.2f}-{t['hi']:.2f})"
              f"   AHI {a['hr']:.2f}" if a else "",
              f"  | adj: T90 {ta['hr']:.2f} p={ta['p']:.3f}, AHI {at['hr']:.2f} p={at['p']:.3f}"
              if ta and at else "")
print(f"\nwritten -> {OUT}/external.json")
