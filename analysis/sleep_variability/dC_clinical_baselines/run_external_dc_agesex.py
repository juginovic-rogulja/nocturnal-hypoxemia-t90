"""
External-cohort discrimination gain on the BARE age + sex baseline: T90 versus sleep
duration, home-recorded PSG.

WHY
---
run_external_dc_tst.py computed the external dC on the frozen FULL clinical risk model.
The three-settings workbook compares that against primary-cohort dC values computed on an
age + sex baseline, so the external row sat on a much richer baseline than settings 1
and 3. This script computes the external dC on the bare baseline too, so all three
settings can be read on age + sex.

BASELINE: THE FROZEN EXTERNAL CONVENTION FOR AGE + SEX MODELS
-------------------------------------------------------------
numbers/freeze_05_external.py, the frozen external pipeline, adjusts its minimally
adjusted models by natural cubic spline on age, cr(df=4) with the first basis column
dropped, plus male in SHHS, and site strata in MrOS, which is men only and so has no sex
term. That is also the primary cohort's baseline functional form, the paper's age + sex
baseline. spline() below is copied verbatim from freeze_05_external.py. The linear age
term in run_external_dc_tst.py belongs to the frozen CLINICAL block only and is not the
pipeline's convention for age + sex models.

METHOD, REUSED, NOT REIMPLEMENTED
---------------------------------
This script IMPORTS run_external_dc_tst, which re-runs the whole frozen clinical
reproduction on import and aborts unless all 7 SHHS outcomes reproduce
numbers/external.json to 4 decimal places. Its loaded frames (s, m), helpers (sframe,
fit_ml, cidx, ladder, rint) and conventions (unpenalized Cox, Newton step scan,
EVENT_FLOOR 60, in-sample C, one shared complete-case frame per outcome) are then reused
unchanged.

POSITIVE CONTROL FOR THIS BASELINE
----------------------------------
external.json froze, per outcome, the minimally adjusted T90 hazard ratio computed on
exactly the frame this script builds: age spline + sex in SHHS, age spline + site strata
in MrOS. Before any C is computed, every outcome's frame must match the frozen n and
events, and the refit T90 HR must reproduce the frozen hr, lo and hi within 5e-4 (v8.3,
both sides full precision; was equality at 3 decimals).
12 of 12 must pass or the run aborts.

C IS STILL IN SAMPLE
--------------------
Same in-sample scoring as the frozen convention, kept and labeled. NOT comparable to the
held-out cross-fit dC of the primary cohort.

SANITY GATE
-----------
The bare baseline leaves more room, so the age + sex dC from T90 must be LARGER than the
clinical-baseline dC from run_external_dc_tst.py on cohort means and the combined mean.
The run exits nonzero if it is not, so the workbook step cannot ship an unexplained
number.

Writes external_dc_agesex.csv and external_dc_agesex_provenance.json into this
directory. Reads everything else. Modifies nothing under T90_Manuscript/.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import warnings
warnings.filterwarnings("ignore")
import json
import sys
import numpy as np
import pandas as pd
from patsy import dmatrix

HERE = f"{paths.SV_ROOT}/dC_clinical_baselines"
sys.path.insert(0, HERE)

print("Importing run_external_dc_tst: this re-runs the frozen clinical-baseline block and")
print("its positive control, which aborts the import on any mismatch with external.json.\n")
import run_external_dc_tst as base  # noqa: E402  executes fully on import

froz = json.load(open(base.FROZEN))


def spline(df, col, n=4):
    """Copied verbatim from numbers/freeze_05_external.py. cr(df=4) natural cubic
    spline, first basis column dropped, columns named agesp1..agesp3."""
    b = dmatrix(f"cr(a, df={n}) - 1", {"a": df[col].values}, return_type="dataframe")
    b.columns = [f"agesp{i}" for i in range(b.shape[1])]
    b.index = df.index
    b = b.iloc[:, 1:]
    for c in b.columns:
        df[c] = b[c]
    return list(b.columns)


def hrres(c, x, f):
    """Copied from freeze_05_external.py res(); v8.3: full precision, as res() now stores."""
    lo, hi = np.exp(c.confidence_intervals_.loc[x])
    return {"n": int(len(f)), "events": int(f.E.sum()),
            "hr": float(np.exp(c.params_[x])), "lo": float(lo),  # v8.3 full precision (2026-09-16)
            "hi": float(hi), "p": float(c.summary.loc[x, "p"])}


def hr_check(got, fz):
    # v8.3: reference now full precision, compare within 5e-4 (2026-09-16)
    ok = (got["n"] == fz["n"] and got["events"] == fz["events"]
          and abs(got["hr"] - fz["hr"]) <= 5e-4 and abs(got["lo"] - fz["lo"]) <= 5e-4
          and abs(got["hi"] - fz["hi"]) <= 5e-4 and abs(got["p"] - fz["p"]) < 1e-9)
    return "reproduced" if ok else "MISMATCH"


rows = []
hr_fail = []

# ==================================================================== SHHS
print("\n" + "=" * 78)
print("SHHS  age + sex baseline (age spline + male), home PSG")
print("=" * 78)
assert len(base.s) == froz["shhs"]["n"], "SHHS merged frame does not match frozen n"
SP_SH = spline(base.s, "age_s1")          # on the full merged frame, as freeze_05 did
BASE_SH = SP_SH + ["male"]
LAB_SH = ("age (natural cubic spline, cr df=4, first basis column dropped) + sex, "
          "the frozen external convention for age + sex models")

for k, dc, pc, lab in base.SOUT:
    g = base.sframe(k, dc, pc)
    fz = froz["shhs"]["outcomes"][lab]["t90"]
    # the frozen minimally adjusted frame, t90 only
    fh = g[["T", "E", "t90_z"] + BASE_SH].dropna()
    got = hrres(base.fit_ml(fh), "t90_z", fh)
    chk = hr_check(got, fz)
    if chk != "reproduced":
        hr_fail.append((lab, got, fz))
    # the shared ladder frame, tst added. slpprdp is complete in the merged frame.
    f = g[["T", "E", "t90_z", "tst_z"] + BASE_SH].dropna()
    assert len(f) == len(fh), f"{lab}: TST changed the sample, {len(f)} vs {len(fh)}"
    if f.E.sum() < base.EVENT_FLOOR:
        print(f"  {lab:<26} SKIPPED, {int(f.E.sum())} events below floor")
        continue
    b0, b1, b2 = base.ladder(f, BASE_SH)
    print(f"  {lab:<26} n={len(f):<5} ev={int(f.E.sum()):<5} "
          f"HR {got['hr']:.3f} (frozen {fz['hr']:.3f}) {chk:<10} "
          f"C0={b0:.4f}  dC_T90={b1 - b0:+.4f}  dC_TST={b2 - b0:+.4f}")
    rows.append({
        "row_type": "per_outcome", "cohort": "SHHS",
        "cohort_n": int(froz["shhs"]["n"]), "outcome": lab,
        "n_model": int(len(f)), "events": int(f.E.sum()),
        "baseline_covariates": LAB_SH,
        "c_baseline": b0, "c_with_t90": b1, "dC_t90": b1 - b0,
        "c_with_tst": b2, "dC_tst": b2 - b0,
        "hr_frame_n": got["n"], "hr_frame_events": got["events"],
        "t90_hr_frozen": fz["hr"], "t90_hr_recomputed": got["hr"],
        "frozen_hr_check": chk,
    })

# ==================================================================== MrOS
print("\n" + "=" * 78)
print("MrOS  age-only baseline (age spline, site-stratified), home PSG, men only")
print("=" * 78)
assert len(base.m) == froz["mros"]["n"], "MrOS analysis frame does not match frozen n"
SP_MR = spline(base.m, "VSAGE1")
LAB_MR = ("age only (natural cubic spline, cr df=4, first basis column dropped), "
          "site-stratified. MrOS is men only, so there is no sex term")
print("  MrOS is men only: the bare baseline is AGE ALONE, no sex term exists to add.")

MOUT = [("DADEAD", "Death from any cause"), ("DACARDIO", "Cardiovascular death"),
        ("DAPULMON", "Pulmonary death"), ("DASTROKE", "Stroke death"),
        ("DACANCER", "Cancer death")]

for k, lab in MOUT:
    fz = froz["mros"]["outcomes"][lab]["t90"]
    fh = pd.DataFrame({"T": (base.m.FUVSDT / 365.25).values,
                       "E": base.m[k].fillna(0).astype(int).values,
                       "t90_z": base.m.t90_z.values, "site": base.m.SITE.values,
                       **{c: base.m[c].values for c in SP_MR}}).dropna()
    got = hrres(base.fit_ml(fh, ["site"]), "t90_z", fh)
    chk = hr_check(got, fz)
    if chk != "reproduced":
        hr_fail.append((lab, got, fz))
    f = pd.DataFrame({"T": (base.m.FUVSDT / 365.25).values,
                      "E": base.m[k].fillna(0).astype(int).values,
                      "t90_z": base.m.t90_z.values, "tst_z": base.m.tst_z.values,
                      "site": base.m.SITE.values,
                      **{c: base.m[c].values for c in SP_MR}}).dropna()
    if f.E.sum() < base.EVENT_FLOOR:
        print(f"  {lab:<26} SKIPPED, {int(f.E.sum())} events below floor")
        continue
    b0, b1, b2 = base.ladder(f, SP_MR + ["site"], strata=["site"])
    print(f"  {lab:<26} n={len(f):<5} ev={int(f.E.sum()):<5} "
          f"HR {got['hr']:.3f} (frozen {fz['hr']:.3f}) {chk:<10} "
          f"C0={b0:.4f}  dC_T90={b1 - b0:+.4f}  dC_TST={b2 - b0:+.4f}")
    rows.append({
        "row_type": "per_outcome", "cohort": "MrOS",
        "cohort_n": int(froz["mros"]["n"]), "outcome": lab,
        "n_model": int(len(f)), "events": int(f.E.sum()),
        "baseline_covariates": LAB_MR,
        "c_baseline": b0, "c_with_t90": b1, "dC_t90": b1 - b0,
        "c_with_tst": b2, "dC_tst": b2 - b0,
        "hr_frame_n": got["n"], "hr_frame_events": got["events"],
        "t90_hr_frozen": fz["hr"], "t90_hr_recomputed": got["hr"],
        "frozen_hr_check": chk,
    })

if hr_fail:
    print("\nPOSITIVE CONTROL FAILED: the age-spline frames do not reproduce the frozen")
    print("minimally adjusted T90 hazard ratios in external.json. No C is trusted.")
    for lab, got, fz in hr_fail:
        print(f"  {lab}: got {got} frozen {fz}")
    sys.exit(1)
print(f"\n  POSITIVE CONTROL PASSED: all {len(rows)} outcomes reproduce the frozen "
      "minimally adjusted T90 HR (n, events, hr, lo, hi to 3 dp) from external.json.")

# ==================================================================== sanity vs clinical
out = pd.DataFrame(rows)
clin = pd.read_csv(f"{HERE}/external_dc_tst.csv")
cmap_t90 = {(r.cohort, r.outcome): r.dC_t90 for _, r in clin.iterrows()}
cmap_tst = {(r.cohort, r.outcome): r.dC_tst for _, r in clin.iterrows()}
out["dC_t90_clinical_baseline"] = [cmap_t90[(c, o)] for c, o in zip(out.cohort, out.outcome)]
out["dC_tst_clinical_baseline"] = [cmap_tst[(c, o)] for c, o in zip(out.cohort, out.outcome)]

print("\n" + "=" * 78)
print("SANITY: age + sex dC from T90 must exceed the clinical-baseline dC")
print("=" * 78)
viol = out[out.dC_t90 <= out.dC_t90_clinical_baseline]
for _, r in out.iterrows():
    flag = "" if r.dC_t90 > r.dC_t90_clinical_baseline else "   <-- NOT larger"
    print(f"  {r.cohort:<5} {r.outcome:<26} agesex {r.dC_t90:+.4f}  "
          f"clinical {r.dC_t90_clinical_baseline:+.4f}{flag}")
sane = True
for coh in ["SHHS", "MrOS"]:
    q = out[out.cohort == coh]
    a, c = q.dC_t90.mean(), q.dC_t90_clinical_baseline.mean()
    ok = a > c
    sane &= ok
    print(f"  {coh} mean: agesex {a:+.4f} vs clinical {c:+.4f}  "
          f"{'OK, larger' if ok else 'FAIL'}")
a, c = out.dC_t90.mean(), out.dC_t90_clinical_baseline.mean()
sane &= a > c
print(f"  Combined 12-outcome mean: agesex {a:+.4f} vs clinical {c:+.4f}  "
      f"{'OK, larger' if a > c else 'FAIL'}")
if len(viol):
    print(f"  Note: {len(viol)} individual outcome(s) not larger, listed above. The gate "
          "is on the means.")

# ==================================================================== means and CSV
mean_rows = []
for coh in ["SHHS", "MrOS"]:
    q = out[out.cohort == coh]
    mean_rows.append({
        "row_type": "cohort_mean", "cohort": coh, "cohort_n": int(q.cohort_n.iloc[0]),
        "outcome": f"MEAN of {len(q)} outcomes",
        "n_model": np.nan, "events": np.nan,
        "baseline_covariates": q.baseline_covariates.iloc[0],
        "c_baseline": q.c_baseline.mean(), "c_with_t90": q.c_with_t90.mean(),
        "dC_t90": q.dC_t90.mean(), "c_with_tst": q.c_with_tst.mean(),
        "dC_tst": q.dC_tst.mean(),
        "hr_frame_n": np.nan, "hr_frame_events": np.nan,
        "t90_hr_frozen": np.nan, "t90_hr_recomputed": np.nan,
        "frozen_hr_check": f"all {len(q)} reproduced",
        "dC_t90_clinical_baseline": q.dC_t90_clinical_baseline.mean(),
        "dC_tst_clinical_baseline": q.dC_tst_clinical_baseline.mean(),
    })
mean_rows.append({
    "row_type": "combined_mean", "cohort": "BOTH", "cohort_n": np.nan,
    "outcome": f"MEAN of {len(out)} outcomes, 7 SHHS + 5 MrOS",
    "n_model": np.nan, "events": np.nan,
    "baseline_covariates": "age + sex in SHHS, age only in MrOS (men only)",
    "c_baseline": out.c_baseline.mean(), "c_with_t90": out.c_with_t90.mean(),
    "dC_t90": out.dC_t90.mean(), "c_with_tst": out.c_with_tst.mean(),
    "dC_tst": out.dC_tst.mean(),
    "hr_frame_n": np.nan, "hr_frame_events": np.nan,
    "t90_hr_frozen": np.nan, "t90_hr_recomputed": np.nan,
    "frozen_hr_check": "all 12 reproduced",
    "dC_t90_clinical_baseline": out.dC_t90_clinical_baseline.mean(),
    "dC_tst_clinical_baseline": out.dC_tst_clinical_baseline.mean(),
})
full = pd.concat([out, pd.DataFrame(mean_rows)], ignore_index=True)

CSV = f"{HERE}/external_dc_agesex.csv"
with open(CSV, "w") as fh:
    fh.write(
        "# EXTERNAL COHORTS, AGE + SEX BASELINE: discrimination gain from T90 and from\n"
        "# home-PSG total sleep time. Sibling of external_dc_tst.csv, which holds the\n"
        "# same quantities on the frozen FULL clinical baseline.\n"
        "# Baseline: SHHS = age spline (cr df=4, first basis column dropped) + sex.\n"
        "#           MrOS = age spline only, site-stratified. MrOS is men only, so there\n"
        "#           is no sex term.\n"
        "# C convention: IN SAMPLE, no split, the frozen external convention, kept and\n"
        "#           labeled. NOT comparable to the primary cohort's held-out cross-fit.\n"
        "# Positive controls: run_external_dc_tst.py was re-executed on import and its 7\n"
        "#           SHHS clinical gains reproduced external.json to 4 dp. All 12 age-\n"
        "#           spline frames reproduce the frozen minimally adjusted T90 HR (n,\n"
        "#           events, hr, lo, hi to 3 dp) from external.json.\n"
        "# row_type: per_outcome, cohort_mean, combined_mean.\n"
        "# Generator: run_external_dc_agesex.py, 2026-08-20.\n")
    full.to_csv(fh, index=False)

prov = {
    "method": dict(base.prov["method"]),
    "baseline": {
        "shhs": LAB_SH, "mros": LAB_MR,
        "source_of_convention": "spline() copied verbatim from numbers/"
                                "freeze_05_external.py, its ADJ/MADJ convention",
        "note": "linear age in run_external_dc_tst.py belongs to the frozen CLINICAL "
                "block only; the pipeline's own age + sex convention is the spline",
    },
    "positive_controls": {
        "clinical_block_reproduction": "run_external_dc_tst.py re-executed on import; "
                                       "all 7 SHHS clinical gains reproduced "
                                       "external.json to 4 dp or the import aborts",
        "t90_hr_reproduction": {
            "target": "external.json per-outcome t90 block (minimally adjusted HR on "
                      "exactly this baseline's frame)",
            "n_checked": int(len(out)), "n_reproduced": int(len(out)), "pass": True},
    },
    "sanity_vs_clinical_baseline": {
        "rule": "age+sex dC_t90 mean must exceed clinical-baseline dC_t90 mean, "
                "per cohort and combined",
        "pass": bool(sane),
        "n_outcomes_not_larger": int(len(viol)),
        "outcomes_not_larger": [f"{r.cohort} {r.outcome}" for _, r in viol.iterrows()],
    },
    "summary": {},
}
for coh in ["SHHS", "MrOS"]:
    q = out[out.cohort == coh]
    prov["summary"][coh] = {
        "n_outcomes": int(len(q)),
        "mean_c_baseline": float(q.c_baseline.mean()),
        "mean_dC_t90": float(q.dC_t90.mean()), "mean_dC_tst": float(q.dC_tst.mean()),
        "n_positive_t90": int((q.dC_t90 > 0).sum()),
        "n_positive_tst": int((q.dC_tst > 0).sum()),
        "mean_dC_t90_clinical_baseline": float(q.dC_t90_clinical_baseline.mean()),
    }
prov["summary"]["combined"] = {
    "n_outcomes": int(len(out)),
    "mean_c_baseline": float(out.c_baseline.mean()),
    "mean_dC_t90": float(out.dC_t90.mean()), "mean_dC_tst": float(out.dC_tst.mean()),
    "n_positive_t90": int((out.dC_t90 > 0).sum()),
    "n_positive_tst": int((out.dC_tst > 0).sum()),
    "mean_dC_t90_clinical_baseline": float(out.dC_t90_clinical_baseline.mean()),
    "mean_dC_tst_clinical_baseline": float(out.dC_tst_clinical_baseline.mean()),
}
with open(f"{HERE}/external_dc_agesex_provenance.json", "w") as fh:
    json.dump(prov, fh, indent=2)

# ==================================================================== final table
print("\n" + "=" * 78)
print("FULL TABLE  age + sex baseline, in-sample C, per 1 SD of exposure")
print("=" * 78)
hdr = (f"  {'cohort':<6}{'outcome':<27}{'n':>6}{'events':>7}{'C(base)':>9}"
       f"{'C(+T90)':>9}{'dC_T90':>9}{'C(+TST)':>9}{'dC_TST':>9}")
print(hdr)
for _, r in out.iterrows():
    print(f"  {r.cohort:<6}{r.outcome:<27}{r.n_model:>6.0f}{r.events:>7.0f}"
          f"{r.c_baseline:>9.4f}{r.c_with_t90:>9.4f}{r.dC_t90:>+9.4f}"
          f"{r.c_with_tst:>9.4f}{r.dC_tst:>+9.4f}")
print()
for coh in ["SHHS", "MrOS"]:
    q = out[out.cohort == coh]
    print(f"  {coh} mean ({len(q)} outcomes):      dC_T90 {q.dC_t90.mean():+.4f}   "
          f"dC_TST {q.dC_tst.mean():+.4f}   C(base) {q.c_baseline.mean():.4f}")
print(f"  Combined mean (12 outcomes): dC_T90 {out.dC_t90.mean():+.4f}   "
      f"dC_TST {out.dC_tst.mean():+.4f}   C(base) {out.c_baseline.mean():.4f}")

print(f"\nwritten -> {CSV}")
print(f"written -> {HERE}/external_dc_agesex_provenance.json")
if not sane:
    print("\nSANITY GATE FAILED: investigate before updating the workbook.")
    sys.exit(2)
print("\nSANITY GATE PASSED: bare-baseline T90 gains exceed the clinical-baseline gains "
      "on every mean.")
