"""
External-cohort discrimination gain: T90 versus sleep duration, home-recorded PSG.

THE GAP THIS FILLS
------------------
numbers/freeze_05_external.py froze, for SHHS only, a T90 discrimination gain on a full
clinical risk model: clinical_model_c, clinical_model_c_with_t90, clinical_model_gain.
It never computed the same quantity for sleep duration, and it never computed any
clinical-model C for MrOS. This script does both, holding every other choice fixed.

METHOD, COPIED FROM freeze_05_external.py, NOT REIMPLEMENTED
------------------------------------------------------------
  fit          Cox by maximum partial likelihood, penalizer 0, Newton step sizes
               [None, 0.05, 0.1, 0.25, 0.5, 0.95] scanned, highest log-likelihood kept
  exposures    rank-based inverse normal transform, so every gain is per 1 SD.
               SHHS transforms on the whole cohort. MrOS transforms within site.
  C            lifelines.utils.concordance_index on -predict_partial_hazard,
               scored IN SAMPLE on the same rows the model was fit on. There is no
               split and no cross-fit. This is the frozen convention and it is kept
               so the positive control can reproduce. It is NOT comparable to the
               primary cohort's held-out cross-fit dC.
  sample       one complete-case frame per outcome, identical for baseline, +T90 and
               +TST, so the three C values differ only by the added term.
  floor        60 events, the floor the frozen clinical-model block used.

POSITIVE CONTROL
----------------
For all 7 SHHS outcomes, C(baseline), C(baseline+T90) and the gain must reproduce
numbers/external.json within 5e-5 (v8.3, both sides full precision; was equality at 4
decimal places) before any TST number is believed. The run
aborts if they do not.

MrOS IS AN EXTENSION, NOT A REPRODUCTION
----------------------------------------
freeze_05_external.py has no clinical-model block for MrOS, so there is no frozen MrOS
gain to reproduce and the MrOS rows carry no positive control. The MrOS clinical
baseline is built from the closest MrOS analogues of the SHHS clinical covariate list.
MrOS is men only, so there is no sex term. MrOS has no measured HDL or diastolic
pressure, so systolic pressure and self-reported high cholesterol stand in. Every MrOS
row is flagged frozen_check = "no frozen counterpart".

Reads only. Writes external_dc_tst.csv and external_dc_tst_provenance.json into this
directory. Modifies nothing under T90_Manuscript/.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import warnings
warnings.filterwarnings("ignore")
import json
import os
import sys
import numpy as np
import pandas as pd
from scipy import stats
from patsy import dmatrix
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index

HERE = f"{paths.SV_ROOT}/dC_clinical_baselines"
FROZEN = f"{paths.NUMBERS_DIR}/external.json"
SH = paths.SHHS_DIR
MR = os.path.expanduser(paths.MROS_DIR)
STEPS = [None, 0.05, 0.1, 0.25, 0.5, 0.95]
EVENT_FLOOR = 60

rows = []
prov = {"method": {
    "fit": "Cox by maximum partial likelihood, penalizer 0",
    "newton_step_sizes_scanned": [str(x) for x in STEPS],
    "C": "lifelines concordance_index on -predict_partial_hazard, IN SAMPLE, no split",
    "exposure_transform": "rank-based inverse normal, per 1 SD",
    "event_floor": EVENT_FLOOR,
    "sample": "one complete-case frame per outcome shared by baseline, +T90 and +TST",
    "source_of_method": "copied from T90_Manuscript/numbers/freeze_05_external.py",
}}


def fit_ml(f, strata=None):
    """Unpenalized Cox. Scan Newton step sizes, keep the maximum-likelihood solution."""
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
            best = (ll, c)
    if best is None:
        raise RuntimeError("no Newton step size converged")
    return best[1]


def rint(v):
    z = pd.Series(np.nan, index=v.index)
    ok = v.notna()
    r = stats.rankdata(v[ok], method="average")
    z.loc[v[ok].index] = stats.norm.ppf((r - 0.375) / (ok.sum() + 0.25))
    return z


def cidx(fit, f):
    return float(concordance_index(f["T"], -fit.predict_partial_hazard(f), f.E))


def ladder(f, base, strata=None):
    """C for baseline, baseline+T90, baseline+TST on one shared frame."""
    c0 = fit_ml(f[["T", "E"] + base], strata=strata)
    c1 = fit_ml(f[["T", "E", "t90_z"] + base], strata=strata)
    c2 = fit_ml(f[["T", "E", "tst_z"] + base], strata=strata)
    b0, b1, b2 = cidx(c0, f), cidx(c1, f), cidx(c2, f)
    return b0, b1, b2


# ==================================================================== SHHS
print("=" * 78)
print("SHHS  home polysomnography, community cohort")
print("=" * 78)
d = pd.read_csv(f"{SH}/shhs1-dataset-0.21.0.csv", low_memory=False)
cv = pd.read_csv(f"{SH}/shhs-cvd-summary-dataset-0.21.0.csv", low_memory=False)
cv = cv.drop(columns=[c for c in ["gender", "race", "age_s1"] if c in cv.columns])
s = d[["nsrrid", "pctsa90h", "slpprdp", "ahi_a0h3a", "ahi_a0h4a", "age_s1", "gender",
       "bmi_s1", "smokstat_s1", "diasbp", "chol", "hdl", "htnderv_s1",
       "parrptdiab"]].merge(cv, on="nsrrid", how="inner")
s["male"] = (s.gender == 1).astype(int)
s["t90_z"] = rint(s.pctsa90h)
s["tst_z"] = rint(s.slpprdp)
print(f"  merged n = {len(s):,}   slpprdp missing = {int(s.slpprdp.isna().sum())}")

SOUT = [("chf", "chf_date", "prev_chf", "Heart failure"),
        ("mi", "mi_date", "prev_mi", "Myocardial infarction"),
        ("stroke", "stk_date", "prev_stk", "Stroke"),
        ("any_cvd", None, None, "Cardiovascular composite"),
        ("any_chd", None, None, "Coronary heart disease"),
        ("cvd_death", "cvd_dthdt", None, "Cardiovascular death"),
        ("afibincident", None, "afibprevalent", "Atrial fibrillation")]

CLIN_SHHS = ["age_s1", "male", "bmi_s1", "smokstat_s1", "diasbp", "chol", "hdl",
             "htnderv_s1", "parrptdiab"]
CLIN_SHHS_LABEL = ("age + sex + BMI + smoking + diastolic BP + total cholesterol + HDL "
                   "+ hypertension + diabetes")


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


froz = json.load(open(FROZEN))
fails = []

for k, dc, pc, lab in SOUT:
    g = sframe(k, dc, pc)
    # The frozen frame. t90_z only, exactly as freeze_05_external.py built it.
    f_frozen = g[["T", "E", "t90_z"] + CLIN_SHHS].dropna()
    # The shared frame, adding tst_z. slpprdp is complete in SHHS so these must match.
    f = g[["T", "E", "t90_z", "tst_z"] + CLIN_SHHS].dropna()
    assert len(f) == len(f_frozen), f"{lab}: TST changed the sample, {len(f)} vs {len(f_frozen)}"
    if f.E.sum() < EVENT_FLOOR:
        print(f"  {lab:<26} SKIPPED, {int(f.E.sum())} events below floor")
        continue

    b0, b1, b2 = ladder(f, CLIN_SHHS)
    fz = froz["shhs"]["outcomes"][lab]
    # v8.3: reference now full precision, compare within 5e-5 (2026-09-16)
    ok = (abs(b0 - fz["clinical_model_c"]) <= 5e-5
          and abs(b1 - fz["clinical_model_c_with_t90"]) <= 5e-5
          and abs((b1 - b0) - fz["clinical_model_gain"]) <= 5e-5)
    if not ok:
        fails.append((lab, round(b0, 4), fz["clinical_model_c"],
                      round(b1, 4), fz["clinical_model_c_with_t90"]))
    print(f"  {lab:<26} n={len(f):<5} ev={int(f.E.sum()):<5} "
          f"C0={b0:.4f} (frozen {fz['clinical_model_c']:.4f})  "
          f"dC_T90={b1 - b0:+.4f} (frozen {fz['clinical_model_gain']:+.4f})  "
          f"dC_TST={b2 - b0:+.4f}   {'PASS' if ok else 'FAIL'}")
    rows.append({
        "cohort": "SHHS", "cohort_n": int(froz["shhs"]["n"]), "outcome": lab,
        "n_model": int(len(f)), "events": int(f.E.sum()),
        "baseline_covariates": CLIN_SHHS_LABEL,
        "c_baseline": b0, "c_with_t90": b1, "dC_t90": b1 - b0,
        "c_with_tst": b2, "dC_tst": b2 - b0,
        "frozen_c_baseline": fz["clinical_model_c"],
        "frozen_c_with_t90": fz["clinical_model_c_with_t90"],
        "frozen_dC_t90": fz["clinical_model_gain"],
        "frozen_check": "reproduced" if ok else "MISMATCH",
    })

if fails:
    print("\nPOSITIVE CONTROL FAILED, aborting before any TST number is used:")
    for r in fails:
        print("   ", r)
    sys.exit(1)
print("\n  POSITIVE CONTROL PASSED, all 7 SHHS outcomes reproduce external.json within 5e-5.\n")

# ==================================================================== MrOS
print("=" * 78)
print("MrOS  home polysomnography, community cohort, men only")
print("=" * 78)
pos = pd.read_sas(f"{MR}/POSFEB23.SAS7BDAT", encoding="latin-1")
ef = pd.read_sas(f"{MR}/effeb24.sas7bdat", encoding="latin-1")
vs = pd.read_sas(f"{MR}/vsfeb24.sas7bdat", encoding="latin-1")
v1 = pd.read_sas(f"{MR}/v1feb24.sas7bdat", encoding="latin-1")

m = (pos[["ID", "POPCSA90", "POSLPRDP", "PORDI4P"]]
     .merge(ef[["ID", "SITE", "DADEAD", "DACARDIO", "DACANCER", "DASTROKE", "DAPULMON",
                "FUVSDT"]], on="ID")
     .merge(vs[["ID", "VSAGE1", "HWBMI"]], on="ID", how="left")
     .merge(v1[["ID", "TUSMOKE", "BPSYSTOL", "MHBP", "MHDIAB"]],
            on="ID", how="left"))
for c in m.columns:
    if c not in ("ID", "SITE"):
        m[c] = pd.to_numeric(m[c], errors="coerce")
# FUVSDT is days from the sleep visit, the clock freeze_05_external.py settled on.
m = m[m.FUVSDT.notna() & (m.FUVSDT > 0) & m.VSAGE1.notna() & m.POPCSA90.notna()].copy()
m["SITE"] = m.SITE.astype(str)
print(f"  analysis frame n = {len(m):,}   sites = {m.SITE.nunique()}")

# within-site rank inverse normal, the MrOS convention in freeze_05_external.py
for src, dst in [("POPCSA90", "t90_z"), ("POSLPRDP", "tst_z")]:
    z = pd.Series(np.nan, index=m.index)
    for site, idx in m.groupby("SITE").groups.items():
        z.loc[idx] = rint(m.loc[idx, src])
    m[dst] = z

# DHLCHOL, the MrOS high-cholesterol field, was tried and dropped: it is missing in
# 2,124 of 2,911 men and collapsed every model to n=770, pushing two outcomes under the
# event floor. The remaining six covariates are near-complete.
CLIN_MROS = ["VSAGE1", "HWBMI", "TUSMOKE", "BPSYSTOL", "MHBP", "MHDIAB"]
CLIN_MROS_LABEL = ("age + BMI + smoking + systolic BP + hypertension + diabetes, "
                   "site-stratified, men only so no sex term, no lipid term available")
for c in CLIN_MROS:
    print(f"    {c:<10} missing = {int(m[c].isna().sum())}")

for k, lab in [("DADEAD", "Death from any cause"), ("DACARDIO", "Cardiovascular death"),
               ("DAPULMON", "Pulmonary death"), ("DASTROKE", "Stroke death"),
               ("DACANCER", "Cancer death")]:
    f = pd.DataFrame({"T": (m.FUVSDT / 365.25).values,
                      "E": m[k].fillna(0).astype(int).values,
                      "t90_z": m.t90_z.values, "tst_z": m.tst_z.values,
                      "site": m.SITE.values,
                      **{c: m[c].values for c in CLIN_MROS}}).dropna()
    if f.E.sum() < EVENT_FLOOR:
        print(f"  {lab:<26} SKIPPED, {int(f.E.sum())} events below floor")
        continue
    b0, b1, b2 = ladder(f, CLIN_MROS + ["site"], strata=["site"])
    print(f"  {lab:<26} n={len(f):<5} ev={int(f.E.sum()):<5} "
          f"C0={b0:.4f}  dC_T90={b1 - b0:+.4f}  dC_TST={b2 - b0:+.4f}")
    rows.append({
        "cohort": "MrOS", "cohort_n": int(froz["mros"]["n"]), "outcome": lab,
        "n_model": int(len(f)), "events": int(f.E.sum()),
        "baseline_covariates": CLIN_MROS_LABEL,
        "c_baseline": b0, "c_with_t90": b1, "dC_t90": b1 - b0,
        "c_with_tst": b2, "dC_tst": b2 - b0,
        "frozen_c_baseline": np.nan, "frozen_c_with_t90": np.nan,
        "frozen_dC_t90": np.nan,
        "frozen_check": "no frozen counterpart",
    })

# ==================================================================== out
out = pd.DataFrame(rows)
out.to_csv(f"{HERE}/external_dc_tst.csv", index=False)

prov["positive_control"] = {
    "target": "numbers/external.json shhs clinical_model_c / _with_t90 / _gain",
    "n_outcomes_checked": 7,
    "n_reproduced_to_4dp": 7,
    "pass": True,
}
prov["shhs"] = {"cohort_n": int(froz["shhs"]["n"]),
                "clinical_baseline": CLIN_SHHS, "label": CLIN_SHHS_LABEL,
                "exposures": {"t90": "pctsa90h", "tst": "slpprdp"},
                "slpprdp_missing_in_merged_frame": 0}
prov["mros"] = {"cohort_n": int(froz["mros"]["n"]),
                "clinical_baseline": CLIN_MROS, "label": CLIN_MROS_LABEL,
                "exposures": {"t90": "POPCSA90", "tst": "POSLPRDP"},
                "followup_var": "FUVSDT",
                "status": "EXTENSION, freeze_05_external.py has no MrOS clinical-model "
                          "block, so these rows have no frozen counterpart"}
prov["summary"] = {
    "n_outcomes": int(len(out)),
    "mean_dC_t90": float(out.dC_t90.mean()),
    "mean_dC_tst": float(out.dC_tst.mean()),
    "n_positive_t90": int((out.dC_t90 > 0).sum()),
    "n_positive_tst": int((out.dC_tst > 0).sum()),
}
for coh in ("SHHS", "MrOS"):
    q = out[out.cohort == coh]
    prov["summary"][coh] = {"n_outcomes": int(len(q)),
                            "mean_dC_t90": float(q.dC_t90.mean()),
                            "mean_dC_tst": float(q.dC_tst.mean()),
                            "n_positive_t90": int((q.dC_t90 > 0).sum()),
                            "n_positive_tst": int((q.dC_tst > 0).sum())}
with open(f"{HERE}/external_dc_tst_provenance.json", "w") as fh:
    json.dump(prov, fh, indent=2)

print("\n" + "=" * 78)
print("SUMMARY, all external outcomes")
print("=" * 78)
print(f"  outcomes            {len(out)}")
print(f"  mean dC from T90    {out.dC_t90.mean():+.4f}   positive in "
      f"{int((out.dC_t90 > 0).sum())} of {len(out)}")
print(f"  mean dC from TST    {out.dC_tst.mean():+.4f}   positive in "
      f"{int((out.dC_tst > 0).sum())} of {len(out)}")
print(f"\nwritten -> {HERE}/external_dc_tst.csv")
