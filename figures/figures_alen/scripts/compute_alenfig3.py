"""
Numbers for AlenFigure3: low nocturnal oxygen is not the same thing as sleep apnea.

The frozen files carry the cross-tabulation (crosstab_v2.json) and the within-severe contrast
against a >10% cut (extras2_v2.json, phenotype_v2.json), but they do not carry a hazard ratio
for every cell of the apnea-by-oxygen cross, and they do not carry a mutual-adjustment fit of
oxygen against the apnea-hypopnea index on the current 19,173 cohort. Those are computed here
and written to numbers/alenfig3_cross_v1.json. Nothing else in numbers/ is touched.

Model, identical in every fit to the primary analysis (numbers/run_primary_unpenalized.py):
  * cohort  = data_frozen_v7_2026-09/t90_final.parquet filtered by cohort_spec.apply_cohort, n = 19,173
  * Cox     = site-stratified, UNPENALIZED (penalizer=0.0)
  * age     = natural cubic spline, cr(a, df=4), with one column dropped. The 4 columns sum to
              one, so an unpenalized solver silently fits nothing unless a column is removed.
  * sex     = male indicator
  * exposure for the mutual-adjustment fit = rank-inverse-normal transform of the measure,
              taken within site, exactly as the primary analysis does

Four blocks are written:
  cross8         8 cells, apnea category by oxygen dichotomised at 10% of the night below 90%,
                 reference = no apnea and oxygen at or below 10%. This is Panel A.
  cross12        12 cells, apnea category by oxygen in 3 levels (<=1%, >1-10%, >10%), reference
                 = no apnea with essentially normal oxygen. This is Panel B, and it is what
                 lets the severe-apnea-with-normal-oxygen cell be read against an outside
                 reference rather than against itself.
  within_severe  severe apnea only, oxygen >10% against oxygen <=1%. The direct contrast
                 behind the sentence the figure has to prove.
  mutual         oxygen and the apnea-hypopnea index entered in the same model.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import json
import sys
import warnings

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from patsy import dmatrix
from scipy import stats

warnings.filterwarnings("ignore")
ROOT = paths.FIGURE_ROOT
sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
from cohort_spec import COHORT_N, apply_cohort, N_MEASURES, N_RANKED_OUTCOMES  # noqa: E402  # v8 sweep 2026-09-12

OUT = f"{paths.NUMBERS_DIR}/alenfig3_cross_v1.json"

# the conditions the figure shows, plus the WHOLE settled negative-control panel of
# 2026-08-07, which must stay near 1.0. The list is imported, never retyped: this file used to
# carry osteoarthritis and back pain only, which is the superseded panel.
from disease_definitions import DISEASES, NEGATIVE_CONTROLS  # noqa: E402

OUTCOMES = [
    ("cvd", "Cardiovascular composite", False),
    ("hf", "Heart failure", False),
    ("resp_failure", "Respiratory failure", False),
    ("pulm_htn", "Pulmonary hypertension", False),
    ("copd2", "COPD", False),
    ("diabetes", "Type 2 diabetes", False),
    ("ckd", "Chronic kidney disease", False),
    ("aki", "Acute kidney injury", False),
    ("death", "Death from any cause", False),
] + [(k, DISEASES[k][0], True) for k in NEGATIVE_CONTROLS]
assert sum(1 for _, _, n in OUTCOMES if n) == 5
assert not {"Osteoarthritis", "Fracture"} & {lab for _, lab, n in OUTCOMES if n}

AHI_CATS = [
    ("None", lambda a: a < 5),
    ("Mild", lambda a: (a >= 5) & (a < 15)),
    ("Moderate", lambda a: (a >= 15) & (a < 30)),
    ("Severe", lambda a: a >= 30),
]

b = apply_cohort(pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet"))
assert len(b) == COHORT_N
b["t90"] = b.spo2_pct_below_90
b["male"] = (b.sex.astype(str).str.upper().str[0] == "M").astype(int)

sp = dmatrix("cr(a, df=4) - 1", {"a": b.AgeAtVisit.values}, return_type="dataframe").iloc[:, 1:]
for i in range(sp.shape[1]):
    b[f"age_s{i}"] = sp.iloc[:, i].values
ADJ = [f"age_s{i}" for i in range(sp.shape[1])] + ["male"]
assert len(ADJ) == 4, "expected 3 spline columns after dropping one, plus sex"


def rint(x):
    r = stats.rankdata(x)
    return stats.norm.ppf((r - 0.375) / (len(r) + 0.25))


b["z_t90"] = b.groupby("site_id").spo2_pct_below_90.transform(rint)
b["z_ahi"] = b.groupby("site_id").AHI.transform(lambda x: rint(x.fillna(x.median())))

# ---------------------------------------------------------------- cell labelling
ahi_cat = pd.Series(index=b.index, dtype=object)
for name, f in AHI_CATS:
    ahi_cat[f(b.AHI)] = name
b["ahi_cat"] = ahi_cat
N_AHI_UNDEFINED = int(b.ahi_cat.isna().sum())   # v7: the index is undefined on nights with under 60 min of true sleep
b = b[b.ahi_cat.notna()].copy()
print(f"nights without a defined AHI dropped from the apnea-group analyses: {N_AHI_UNDEFINED}")

b["ox2"] = np.where(b.t90 > 10, ">10%", "<=10%")
b["ox3"] = np.where(b.t90 <= 1, "<=1%", np.where(b.t90 <= 10, ">1-10%", ">10%"))
b["cell8"] = b.ahi_cat + "|" + b.ox2
b["cell12"] = b.ahi_cat + "|" + b.ox3

REF8, REF12 = "None|<=10%", "None|<=1%"
LEV8 = [f"{a}|{o}" for a, _ in AHI_CATS for o in ("<=10%", ">10%")]
LEV12 = [f"{a}|{o}" for a, _ in AHI_CATS for o in ("<=1%", ">1-10%", ">10%")]

# R30 (Alen 2026-09-08): THREE apnea groups, no and mild pooled (AHI under 15). Two crosses on them:
#   cross9  3 groups x 3 oxygen levels (<=1%, >1-10%, >10%), reference = no or mild apnea with T90 <=1%
#   cross6  3 groups x 2 oxygen levels (<=1% preserved, >10% low), the Figure 3c design; nights with T90 in
#           (1, 10] are set aside for that model exactly as the six-cell duration model sets aside its middle band
AHI3_CATS = [("No or mild", ("None", "Mild")), ("Moderate", ("Moderate",)), ("Severe", ("Severe",))]
b["ahi3"] = b.ahi_cat.map({"None": "No or mild", "Mild": "No or mild", "Moderate": "Moderate", "Severe": "Severe"})
assert b.ahi3.notna().all()
b["cell9"] = b.ahi3 + "|" + b.ox3
b["cell6"] = b.ahi3 + "|" + b.ox3
b6 = b[b.ox3 != ">1-10%"].copy()
LEV9 = [f"{a}|{o}" for a, _ in AHI3_CATS for o in ("<=1%", ">1-10%", ">10%")]
LEV6 = [f"{a}|{o}" for a, _ in AHI3_CATS for o in ("<=1%", ">10%")]
REF9 = REF6 = "No or mild|<=1%"


def cell_model(frame, cellcol, levels, ref, key):
    """Cox on a categorical cell variable. Returns one hazard ratio per non-reference cell."""
    yc, ec, pc = f"{key}_years", f"{key}_incident", f"{key}_prevalent"
    g = frame[(frame[pc] == 0) & frame[yc].notna() & (frame[yc] > 0)]
    X = pd.get_dummies(g[cellcol]).astype(int)
    for lv in levels:
        if lv not in X:
            X[lv] = 0
    X = X[levels]
    cols = {"T": g[yc].values, "E": g[ec].astype(int).values, "site": g.site_id.values}
    for c in ADJ:
        cols[c] = g[c].values
    d = pd.concat([pd.DataFrame(cols, index=g.index), X], axis=1).dropna()
    if d.E.sum() < 60:
        return None
    use = [lv for lv in levels if lv != ref and d[lv].sum() >= 25]
    fit = CoxPHFitter(penalizer=0.0).fit(d[["T", "E", "site"] + ADJ + use], "T", "E",
                                         strata=["site"])
    out = {"events": int(d.E.sum()), "n": int(len(d)), "reference": ref, "cells": {}}
    for lv in levels:
        n_lv, ev_lv = int((d[lv] == 1).sum()), int(d.loc[d[lv] == 1, "E"].sum())
        if lv == ref:
            out["cells"][lv] = {"hr": 1.0, "lo": None, "hi": None, "p": None,
                                "n": n_lv, "events": ev_lv, "reference": True}
        elif lv in use:
            r = fit.summary.loc[lv]
            out["cells"][lv] = {"hr": float(r["exp(coef)"]),  # v8.3 full precision (2026-09-16)
                                "lo": float(r["exp(coef) lower 95%"]),
                                "hi": float(r["exp(coef) upper 95%"]),
                                "p": float(r["p"]), "n": n_lv, "events": ev_lv,
                                "reference": False}
        else:
            out["cells"][lv] = {"hr": None, "n": n_lv, "events": ev_lv, "reference": False}
    return out


def simple_model(frame, xcols, key):
    """Cox with the named continuous or binary columns entered together."""
    yc, ec, pc = f"{key}_years", f"{key}_incident", f"{key}_prevalent"
    g = frame[(frame[pc] == 0) & frame[yc].notna() & (frame[yc] > 0)]
    d = pd.DataFrame({"T": g[yc], "E": g[ec].astype(int), "site": g.site_id,
                      **{c: g[c] for c in xcols + ADJ}}).dropna()
    if d.E.sum() < 60:
        return None
    fit = CoxPHFitter(penalizer=0.0).fit(d[["T", "E", "site"] + xcols + ADJ], "T", "E",
                                         strata=["site"])
    out = {"events": int(d.E.sum()), "n": int(len(d))}
    for x in xcols:
        r = fit.summary.loc[x]
        out[x] = {"hr": float(r["exp(coef)"]),  # v8.3 full precision (2026-09-16)
                  "lo": float(r["exp(coef) lower 95%"]),
                  "hi": float(r["exp(coef) upper 95%"]), "p": float(r["p"])}
    return out


R = {
    "provenance": {
        "cohort_n": int(len(b)),
        "source": "data_frozen_v8_2026-09/t90_final.parquet via numbers/cohort_spec.apply_cohort",
        "model": "site-stratified Cox, unpenalized, cr(a, df=4) age spline with one column "
                 "dropped, plus sex",
        "created_by": "figures_alen/scripts/compute_alenfig3.py",
        "note": "created because no frozen file carries a hazard ratio per cell of the "
                "apnea-by-oxygen cross, nor a mutual-adjustment fit on the 19,173 cohort",
    },
    "cell_sizes": {
        "cross8": {lv: int((b.cell8 == lv).sum()) for lv in LEV8},
        "cross12": {lv: int((b.cell12 == lv).sum()) for lv in LEV12},
        "cross9": {lv: int((b.cell9 == lv).sum()) for lv in LEV9},
        "cross6": {lv: int((b6.cell6 == lv).sum()) for lv in LEV6},
        "cross6_set_aside_t90_gt1_le10": int((b.ox3 == ">1-10%").sum()),
    },
    "ahi3_groups": {"rule": {"No or mild": "AHI < 15", "Moderate": "15 <= AHI < 30", "Severe": "AHI >= 30"},
                    "n": {a: int((b.ahi3 == a).sum()) for a, _ in AHI3_CATS},
                    "note": "R30, Alen 2026-09-08: no and mild apnea pooled; cross9 = 3 x 3 levels, cross6 = 3 x 2 levels (Fig 3c design)"},
    "cross8": {}, "cross12": {}, "cross9": {}, "cross6": {}, "within_severe": {}, "mutual": {},
}

sev = b[b.ahi_cat == "Severe"].copy()
sev["ox_gt10"] = (sev.t90 > 10).astype(int)
sev["ox_1to10"] = ((sev.t90 > 1) & (sev.t90 <= 10)).astype(int)
R["severe_n"] = {"n": int(len(sev)), "n_normal_oxygen_le1": int((sev.t90 <= 1).sum()),
                 "pct_normal_oxygen_le1": round(100 * float((sev.t90 <= 1).mean()), 1),
                 "n_low_oxygen_gt10": int((sev.t90 > 10).sum()),
                 "pct_low_oxygen_gt10": round(100 * float((sev.t90 > 10).mean()), 1)}

for key, label, neg in OUTCOMES:
    if f"{key}_incident" not in b.columns:
        print(f"  skip {label}, not in the frozen cohort")
        continue
    r8 = cell_model(b, "cell8", LEV8, REF8, key)
    r12 = cell_model(b, "cell12", LEV12, REF12, key)
    r9 = cell_model(b, "cell9", LEV9, REF9, key)
    r6 = cell_model(b6, "cell6", LEV6, REF6, key)
    rs = simple_model(sev, ["ox_1to10", "ox_gt10"], key)
    rm = simple_model(b, ["z_t90", "z_ahi"], key)
    if r8:
        R["cross8"][label] = {**r8, "negative_control": neg}
    if r12:
        R["cross12"][label] = {**r12, "negative_control": neg}
    if r9:
        R["cross9"][label] = {**r9, "negative_control": neg}
    if r6:
        R["cross6"][label] = {**r6, "negative_control": neg}
    if rs:
        R["within_severe"][label] = {**rs, "reference": "severe apnea, T90 <=1%",
                                     "negative_control": neg}
    if rm:
        R["mutual"][label] = {**rm, "negative_control": neg}
    print(f"  {label:28s} done")

# ---------------------------------------------------------------- ranking bands, Panel D
# No new model here. These are read straight out of the frozen ranking files and reduced to the
# span of the 5 oxygen measures, because a single rank is not defensible inside the no-apnea
# stratum and a band is.
OXY_FEATURES = ["spo2_nadir_corrected", "spo2_pct_below_90", "odi3_total",
                "spo2_pct_below_88", "odi4_total"]
NUM = f"{paths.NUMBERS_DIR}"


def band_from(df, source, stratum, n_outcomes, n_patients):
    r = {f: int(df.loc[df.feature == f, "rank"].iloc[0]) for f in OXY_FEATURES
         if (df.feature == f).any()}
    return {"stratum": stratum, "source": source, "n_measures": int(len(df)),
            "n_outcomes": n_outcomes, "n_patients": n_patients,
            "oxygen_ranks": r, "oxygen_band": [min(r.values()), max(r.values())],
            "t90_rank": r["spo2_pct_below_90"],
            "ahi_rank": int(df.loc[df.feature == "AHI", "rank"].iloc[0])}


full = pd.read_csv(f"{NUM}/ranking_v3.csv", comment="#")
noap = pd.read_csv(f"{NUM}/ranking_noapnea_stratum1_no_apnea.csv")
strict = pd.read_csv(f"{NUM}/ranking_noapnea_stratum1_uncertainty.csv")
strict = strict.sort_values("mean", ascending=False).reset_index(drop=True)
strict["rank"] = strict.index + 1

R["ranking_bands"] = [
    band_from(full, "numbers/ranking_v3.csv", "Whole cohort", N_RANKED_OUTCOMES, COHORT_N),
    band_from(noap, "numbers/ranking_noapnea_stratum1_no_apnea.csv",
              "No sleep apnea (AHI <5)", 19, 3017),
    band_from(strict, "numbers/ranking_noapnea_stratum1_uncertainty.csv",
              "No sleep apnea, stricter event rule", 14, 3017),
]

json.dump(R, open(OUT, "w"), indent=1)
print(f"\nwritten -> {OUT}")

print(f"\nRank among {N_MEASURES} sleep measures, oxygen family span against the apnea-hypopnea index")
for r in R["ranking_bands"]:
    print(f"  {r['stratum']:38s} oxygen {r['oxygen_band'][0]:>3}-{r['oxygen_band'][1]:<3}"
          f"  apnea index {r['ahi_rank']:>3}   ({r['n_outcomes']} conditions)")

# ------------------------------------------------------------------ console read-back
print(f"\ncohort {len(b):,}   cells of the cross, apnea category by oxygen\n")
hdr = f"{'':12s}" + "".join(f"{o:>12}" for o in ("<=10%", ">10%"))
print(hdr)
for a, _ in AHI_CATS:
    print(f"{a:12s}" + "".join(f"{R['cell_sizes']['cross8'][f'{a}|{o}']:12,}"
                               for o in ("<=10%", ">10%")))

print("\nPanel A, cardiovascular composite, reference = no apnea with oxygen at or below 10%")
c = R["cross8"]["Cardiovascular composite"]["cells"]
for lv in LEV8:
    v = c[lv]
    hr = "1.00 (reference)" if v.get("reference") else (
        f"{v['hr']:.2f} ({v['lo']:.2f}-{v['hi']:.2f})" if v.get("hr") else "not estimated")
    print(f"  {lv:20s} n={v['n']:6,}  events={v['events']:5,}   {hr}")

print("\nR30 cross6 (three apnea groups x preserved/low oxygen), reference = no or mild apnea with T90 <=1%")
for lab, v in R["cross6"].items():
    cells = v["cells"]
    print(f"  {lab:28s} " + "  ".join(f"{lv.split('|')[0][:8]}|{lv.split('|')[1]}: " + ("ref" if cells[lv].get("reference") else (f"{cells[lv]['hr']:.2f}" if cells[lv].get("hr") else "NE")) for lv in LEV6))

print("\nPanel B, severe apnea only, oxygen >10% against oxygen <=1%")
for lab, v in R["within_severe"].items():
    g = v["ox_gt10"]
    print(f"  {lab:28s} ev={v['events']:5,}   {g['hr']:.2f} ({g['lo']:.2f}-{g['hi']:.2f})")

print("\nMutual adjustment, oxygen and the apnea-hypopnea index in the same model")
for lab, v in R["mutual"].items():
    t, a = v["z_t90"], v["z_ahi"]
    print(f"  {lab:28s} oxygen {t['hr']:.2f} ({t['lo']:.2f}-{t['hi']:.2f})"
          f"   apnea index {a['hr']:.2f} ({a['lo']:.2f}-{a['hi']:.2f})")
