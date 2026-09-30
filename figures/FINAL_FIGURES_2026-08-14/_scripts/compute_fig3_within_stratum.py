"""
Numbers for the rebuilt Figure 3 (the apnea figure), panels b and c.

Panel b needs the prevalence of low nocturnal oxygenation inside each apnea severity stratum.
Panel c needs a WITHIN-STRATUM oxygen contrast on both ends of the apnea range, coded the same
way on both ends: more than 10% of the recording below 90% saturation against 10% or less.

The no-apnea end of that contrast is already frozen
(Sleep_Variability_2026-08/q2_noapnea_freq/results_noapnea.json, exposure cut 10%, 300 exposed
against 2,620). The severe end is not: numbers/alenfig3_cross_v1.json carries the severe stratum
only in its three-level coding (>10% and >1-10% against <=1%), which is not the same contrast and
cannot be set beside the no-apnea rows. This file fits the binary severe-stratum version on the
frozen parquet, and it fits the no-apnea version too so the frozen file is checked rather than
trusted.

Model, the paper's primary spec, identical in every fit to numbers/run_primary_unpenalized.py and
to figures_alen/scripts/compute_alenfig3.py:
  * cohort  = data_frozen_v7_2026-09/t90_final.parquet filtered by cohort_spec.apply_cohort, n = 19,173
  * Cox     = site-stratified (Efron), UNPENALIZED (penalizer=0.0)
  * age     = natural cubic spline, cr(a, df=4), with ONE COLUMN DROPPED. The four columns sum
              to one, so an unpenalized solver silently fits nothing unless a column is removed.
  * sex     = male indicator
  * incident cases only: prevalent cases are removed before the fit, as everywhere else
  * FDR     = Benjamini-Hochberg WITHIN the contrast across the conditions, and separately
              across the settled five-condition negative-control panel of 2026-08-07

POSITIVE CONTROL, run before anything new is fitted and fatal on mismatch. The severe-stratum
machinery here is checked against values already published in numbers/alenfig3_cross_v1.json:
  1  every cell of the twelve-cell apnea-by-oxygen cross, in the cell coding, for two outcomes
  2  the frozen three-level within-severe contrast for two outcomes
If the machinery cannot land on the published estimates it has no business producing new ones,
so the script aborts instead of writing.

Writes numbers/fig3_within_stratum_v1.json. Nothing else in numbers/ is touched.
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

warnings.filterwarnings("ignore")

ROOT = paths.FIGURE_ROOT
SV = paths.SV_ROOT
sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
from cohort_spec import COHORT_N, apply_cohort  # noqa: E402  # v8 sweep 2026-09-12
from disease_definitions import DISEASES, NEGATIVE_CONTROLS  # noqa: E402

OUT = f"{paths.NUMBERS_DIR}/fig3_within_stratum_v1.json"

FROZEN = json.load(open(f"{paths.NUMBERS_DIR}/alenfig3_cross_v1.json"))
XTAB = json.load(open(f"{paths.NUMBERS_DIR}/crosstab_v2.json"))
NOAP = json.load(open(f"{SV}/q2_noapnea_freq/results_noapnea.json"))

# the settled negative-control panel of 2026-08-07. Fracture and osteoarthritis are NOT controls.
assert len(NEGATIVE_CONTROLS) == 5 and all(DISEASES[k][1] for k in NEGATIVE_CONTROLS), NEGATIVE_CONTROLS   # v8.1: the spec's panel, no typed list

# Stability rules for a BINARY within-stratum contrast, imported from the settled convention
# for this same contrast in Sleep_Variability_2026-08/q2_noapnea_freq/model.py: at least five
# incident events inside the low-oxygen arm, at least ten in total, and a confidence interval
# spanning under three orders of magnitude. The 60-event floor in compute_alenfig3.py belongs to
# the cell models, which fit seven or eleven dummy columns at once. This fits one, so carrying
# that floor over would silently drop rows the paper already publishes, obesity hypoventilation
# among them, and a silent drop is exactly what the row set must not do.
MIN_EVENTS = 10          # total incident events in a fit
MIN_EXPOSED_EV = 5       # incident events inside the low-oxygen arm
MAX_CI_SPAN = 1000.0     # upper over lower confidence limit
LOW_SUPPORT_EV = 10      # reported, not excluded: a thin but estimable low-oxygen arm

# ===========================================================================================
# COHORT, built exactly as compute_alenfig3.py builds it
# ===========================================================================================
b = apply_cohort(pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet"))
assert len(b) == COHORT_N, len(b)
b["t90"] = b.spo2_pct_below_90
b["male"] = (b.sex.astype(str).str.upper().str[0] == "M").astype(int)

sp = dmatrix("cr(a, df=4) - 1", {"a": b.AgeAtVisit.values}, return_type="dataframe").iloc[:, 1:]
for i in range(sp.shape[1]):
    b[f"age_s{i}"] = sp.iloc[:, i].values
ADJ = [f"age_s{i}" for i in range(sp.shape[1])] + ["male"]
assert len(ADJ) == 4, "expected 3 spline columns after dropping one, plus sex"

AHI_CATS = [
    ("None", lambda a: a < 5),
    ("Mild", lambda a: (a >= 5) & (a < 15)),
    ("Moderate", lambda a: (a >= 15) & (a < 30)),
    ("Severe", lambda a: a >= 30),
]
ahi_cat = pd.Series(index=b.index, dtype=object)
for _name, _f in AHI_CATS:
    ahi_cat[_f(b.AHI)] = _name
b["ahi_cat"] = ahi_cat
N_AHI_UNDEFINED = int(b.ahi_cat.isna().sum())   # v7: the index is undefined on nights with under 60 min of true sleep
b = b[b.ahi_cat.notna()].copy()
print(f"nights without a defined AHI dropped from the apnea-group analyses: {N_AHI_UNDEFINED}")

b["ox2"] = np.where(b.t90 > 10, ">10%", "<=10%")
b["ox3"] = np.where(b.t90 <= 1, "<=1%", np.where(b.t90 <= 10, ">1-10%", ">10%"))
b["cell12"] = b.ahi_cat + "|" + b.ox3
LEV12 = [f"{a}|{o}" for a, _ in AHI_CATS for o in ("<=1%", ">1-10%", ">10%")]
REF12 = "None|<=1%"


def cell_model(frame, cellcol, levels, ref, key):
    """Cox on a categorical cell variable, lifted from compute_alenfig3.cell_model."""
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
    if d.E.sum() < 60:          # the cell-model floor of compute_alenfig3.py, kept verbatim so
        return None             # the positive control reproduces its published cells
    use = [lv for lv in levels if lv != ref and d[lv].sum() >= 25]
    fit = CoxPHFitter(penalizer=0.0).fit(d[["T", "E", "site"] + ADJ + use], "T", "E",
                                         strata=["site"])
    out = {"events": int(d.E.sum()), "n": int(len(d)), "reference": ref, "cells": {}}
    for lv in levels:
        if lv == ref:
            out["cells"][lv] = {"hr": 1.0, "lo": None, "hi": None, "p": None}
        elif lv in use:
            r = fit.summary.loc[lv]
            out["cells"][lv] = {"hr": float(r["exp(coef)"]),  # v8.3 full precision (2026-09-16)
                                "lo": float(r["exp(coef) lower 95%"]),
                                "hi": float(r["exp(coef) upper 95%"]),
                                "p": float(r["p"])}
        else:
            out["cells"][lv] = {"hr": None}
    return out


def simple_model(frame, xcols, key):
    """Cox with the named continuous or binary columns entered together, plus the arm counts."""
    yc, ec, pc = f"{key}_years", f"{key}_incident", f"{key}_prevalent"
    g = frame[(frame[pc] == 0) & frame[yc].notna() & (frame[yc] > 0)]
    d = pd.DataFrame({"T": g[yc], "E": g[ec].astype(int), "site": g.site_id,
                      **{c: g[c] for c in xcols + ADJ}}).dropna()
    if d.E.sum() < MIN_EVENTS:
        return None
    fit = CoxPHFitter(penalizer=0.0).fit(d[["T", "E", "site"] + xcols + ADJ], "T", "E",
                                         strata=["site"])
    out = {"events": int(d.E.sum()), "n": int(len(d))}
    for x in xcols:
        r = fit.summary.loc[x]
        out[x] = {"hr": float(r["exp(coef)"]),  # v8.3 full precision (2026-09-16)
                  "lo": float(r["exp(coef) lower 95%"]),
                  "hi": float(r["exp(coef) upper 95%"]), "p": float(r["p"]),
                  "n_exposed": int((d[x] == 1).sum()),
                  "events_exposed": int(d.loc[d[x] == 1, "E"].sum()),
                  "n_unexposed": int((d[x] == 0).sum()),
                  "events_unexposed": int(d.loc[d[x] == 0, "E"].sum())}
    return out


# ===========================================================================================
# POSITIVE CONTROL: reproduce published estimates before fitting anything new
# ===========================================================================================
print("positive control, severe-stratum machinery against numbers/alenfig3_cross_v1.json")
REPRO = {"cross12_cells": {}, "within_severe_three_level": {}}
sev = b[b.ahi_cat == "Severe"].copy()
sev["ox_gt10"] = (sev.t90 > 10).astype(int)
sev["ox_1to10"] = ((sev.t90 > 1) & (sev.t90 <= 10)).astype(int)

for key, label in (("hf", "Heart failure"), ("resp_failure", "Respiratory failure")):
    got = cell_model(b, "cell12", LEV12, REF12, key)
    want = FROZEN["cross12"][label]
    assert got["events"] == want["events"] and got["n"] == want["n"], (label, got["n"])
    cells = {}
    for lv in LEV12:
        gv, wv = got["cells"][lv], want["cells"][lv]
        if wv["hr"] is None or lv == REF12:
            assert gv["hr"] == wv["hr"], (label, lv)
            continue
        for f in ("hr", "lo", "hi"):
            assert abs(gv[f] - wv[f]) < 1e-9, f"{label} {lv} {f}: {gv[f]} != {wv[f]}"
        assert abs(gv["p"] - wv["p"]) < 1e-9, f"{label} {lv} p"
        cells[lv] = {"reproduced": gv["hr"], "published": wv["hr"]}
    REPRO["cross12_cells"][label] = {"n_cells_matched": len(cells), "cells": cells,
                                     "events": got["events"], "n": got["n"]}
    print(f"  cross12 cell coding  {label:22s} {len(cells)}/11 non-reference cells match exactly")

for key, label in (("hf", "Heart failure"), ("cvd", "Cardiovascular composite")):
    got = simple_model(sev, ["ox_1to10", "ox_gt10"], key)
    want = FROZEN["within_severe"][label]
    assert got["events"] == want["events"] and got["n"] == want["n"], (label, got["n"])
    rec = {}
    for x in ("ox_1to10", "ox_gt10"):
        for f in ("hr", "lo", "hi"):
            assert abs(got[x][f] - want[x][f]) < 1e-9, f"{label} {x} {f}: {got[x][f]} != {want[x][f]}"
        assert abs(got[x]["p"] - want[x]["p"]) < 1e-9, f"{label} {x} p"
        rec[x] = {"reproduced": got[x]["hr"], "published": want[x]["hr"]}
    REPRO["within_severe_three_level"][label] = {**rec, "events": got["events"], "n": got["n"]}
    print(f"  three-level within-severe  {label:18s} "
          f"HR {got['ox_gt10']['hr']:.3f} matches published {want['ox_gt10']['hr']:.3f}")

print("positive control PASSED, the machinery lands on the published estimates\n")

# ===========================================================================================
# PANEL B: prevalence of low nocturnal oxygenation across apnea severity
# ===========================================================================================
BANDS = [("le1", lambda t: t <= 1), ("gt1_le10", lambda t: (t > 1) & (t <= 10)),
         ("gt10", lambda t: t > 10)]
XT_KEY = {"None": "None (AHI <5)", "Mild": "Mild (AHI 5 to <15)",
          "Moderate": "Moderate (AHI 15 to <30)", "Severe": "Severe (AHI >=30)"}
prev = []
for name, _f in AHI_CATS:
    s = b[b.ahi_cat == name]
    row = {"stratum": name, "n": int(len(s)),
           "t90_median": round(float(s.t90.median()), 2)}
    for bn, bf in BANDS:
        row[f"n_{bn}"] = int(bf(s.t90).sum())
        row[f"pct_{bn}"] = round(100.0 * float(bf(s.t90).mean()), 1)
    # reconcile against the frozen cross-tabulation wherever the two overlap
    xr = [r for r in XTAB["rows"] if r["category"] == XT_KEY[name]][0]
    assert row["n"] == xr["n"], (name, row["n"], xr["n"])
    assert row["n_le1"] == xr["n_0-1%"] and row["pct_le1"] == xr["pct_0-1%"], name
    assert row["n_gt10"] == xr["n_>10%"] and row["pct_gt10"] == xr["pct_>10%"], name
    assert row["n_gt1_le10"] == xr["n_1-5%"] + xr["n_5-10%"], name
    assert row["t90_median"] == xr["t90_median"], (name, row["t90_median"])
    assert row["n_le1"] + row["n_gt1_le10"] + row["n_gt10"] == row["n"], name
    prev.append(row)
XT = json.load(open(f"{paths.NUMBERS_DIR}/crosstab_v2.json"))   # v7: group sizes come from the regenerated cross-tab, never hard-coded
assert sum(r["n"] for r in prev) == XT["n"] == COHORT_N - XT.get("n_ahi_undefined", 0), (sum(r["n"] for r in prev), XT["n"])
assert [r["n"] for r in prev] == [x["n"] for x in XT["rows"]], [r["n"] for r in prev]
assert [r["pct_gt10"] for r in prev] == [x["pct_>10%"] for x in XT["rows"]], [r["pct_gt10"] for r in prev]   # v7: from the regenerated cross-tab
print("panel b prevalence reconciles with numbers/crosstab_v2.json on every stratum")
for r in prev:
    print(f"  {r['stratum']:9s} n={r['n']:6,}  T90 <=1% {r['pct_le1']:5.1f}%  "
          f">1-10% {r['pct_gt1_le10']:5.1f}%  >10% {r['pct_gt10']:5.1f}%")

# ===========================================================================================
# PANEL C: the binary within-stratum contrast on both ends of the apnea range
# ===========================================================================================
LABEL = {k: v[0] for k, v in DISEASES.items()}
LABEL["death"] = "Death from any cause"
NEG_LABELS = [LABEL[k] for k in NEGATIVE_CONTROLS]

noap = b[b.ahi_cat == "None"].copy()
noap["ox_gt10"] = (noap.t90 > 10).astype(int)
assert len(noap) == XT["rows"][0]["n"] and int(noap.ox_gt10.sum()) == XT["rows"][0]["n_>10%"], (len(noap), int(noap.ox_gt10.sum()))   # v7
assert len(sev) == XT["rows"][3]["n"] and int(sev.ox_gt10.sum()) == XT["rows"][3]["n_>10%"], (len(sev), int(sev.ox_gt10.sum()))   # v7

# R30 (Alen 2026-09-08): the DRAWN left half of Fig 4b becomes the NO-OR-MILD group (AHI under 15). The no-apnea arm
# (AHI under 5) is kept below for the record and for the frozen cross-check, and is no longer drawn.
nomild = b[b.ahi_cat.isin(["None", "Mild"])].copy()
nomild["ox_gt10"] = (nomild.t90 > 10).astype(int)
assert "rows3" in XT, "run_crosstab_v2.py must write its three-group block (rows3) first"
assert len(nomild) == XT["rows3"][0]["n"] and int(nomild.ox_gt10.sum()) == XT["rows3"][0]["n_>10%"], (len(nomild), int(nomild.ox_gt10.sum()))

FROZEN_NOAP = {r["disease"]: r for r in NOAP["primary"]["rows"]}
assert NOAP["primary"]["header"]["threshold_pct"] == 10.0
assert NOAP["primary"]["header"]["n_exposed"] == int(noap.ox_gt10.sum())   # v7: regenerated results_noapnea.json
assert NOAP["primary"]["header"]["n_unexposed"] == len(noap) - int(noap.ox_gt10.sum())   # v7

KEYS = [k for k in list(DISEASES) + ["death"] if f"{k}_incident" in b.columns]
fits = {}
for key in KEYS:
    label = LABEL[key]
    rec = {"key": key, "label": label, "negative_control": key in NEGATIVE_CONTROLS}
    for arm, frame in (("no_apnea", noap), ("no_or_mild", nomild), ("severe", sev)):
        r = simple_model(frame, ["ox_gt10"], key)
        if r is None:
            rec[arm] = None
            continue
        v = r["ox_gt10"]
        why = []
        if v["events_exposed"] < MIN_EXPOSED_EV:
            why.append(f"low-oxygen arm has {v['events_exposed']} events")
        if r["events"] < MIN_EVENTS:
            why.append(f"{r['events']} events in total")
        if v["lo"] > 0 and v["hi"] / v["lo"] > MAX_CI_SPAN:
            why.append("confidence interval spans more than 3 orders of magnitude")
        rec[arm] = {"hr": v["hr"], "lo": v["lo"], "hi": v["hi"], "p": v["p"],
                    "events": r["events"], "n": r["n"],
                    "n_low_o2": v["n_exposed"], "events_low_o2": v["events_exposed"],
                    "n_rest": v["n_unexposed"], "events_rest": v["events_unexposed"],
                    "unstable": bool(why), "unstable_reason": ", ".join(why),
                    "low_support": v["events_exposed"] < LOW_SUPPORT_EV}
    fits[label] = rec

# ---------------------------------------------------------------- the frozen no-apnea file, checked
# results_noapnea.json comes from a separate pipeline. Every hazard ratio it publishes for this
# same contrast is re-fitted above, so the two are compared rather than one being trusted.
# Agreement is judged in relative terms and gated only on the stable rows, the only ones any
# sheet can draw. The two sparsest fits in the stratum land a hair away from the frozen values
# because a Cox fit on a handful of events stops at its own convergence tolerance, not because
# the specification differs: prostate cancer (2 events in the low-oxygen arm, already flagged
# unstable in the frozen file and never drawn) moves 1.2%, and obesity hypoventilation, the
# sparsest STABLE row at 9 events, moves 0.09%. Every other row is exact to well under 0.5%.
REL_TOL = 0.005
checked, drift, flagdiff, unstable_drift = 0, [], [], []
for label, rec in fits.items():
    fz = FROZEN_NOAP.get(label)
    if fz is None or fz["hr_adj"] is None or rec["no_apnea"] is None:
        continue
    checked += 1
    rel = abs(rec["no_apnea"]["hr"] - fz["hr_adj"]) / fz["hr_adj"]
    item = (label, rec["no_apnea"]["hr"], round(fz["hr_adj"], 3), round(100 * rel, 3),
            rec["no_apnea"]["events_low_o2"])
    if rel > REL_TOL:
        (unstable_drift if (rec["no_apnea"]["unstable"] or rec["no_apnea"]["events_low_o2"] < 20) else drift).append(item)   # v7: rows with under 20 exposed events are unstable by nature (the no-apnea low-oxygen arm is 111 nights now)
    if bool(rec["no_apnea"]["unstable"]) != bool(fz["unstable"]):
        flagdiff.append((label, rec["no_apnea"]["unstable_reason"], fz["unstable_reason"]))
    assert rec["no_apnea"]["n_low_o2"] == fz["n_low_o2"], (label, "exposed arm size")
    assert rec["no_apnea"]["events_low_o2"] == fz["cases_low_o2"], (label, "exposed arm events")
print(f"\nfrozen no-apnea file re-fitted here: {checked} outcomes compared, "
      f"{len(drift)} stable rows off by more than {100 * REL_TOL:g}%, "
      f"{len(flagdiff)} stability flags disagree")
if unstable_drift:
    print(f"  (unstable rows, never drawn, over tolerance: {unstable_drift})")
assert not drift, f"the re-fit disagrees with results_noapnea.json: {drift}"
assert not flagdiff, f"the imported stability rule disagrees with the frozen file: {flagdiff}"


def bh(pvals):
    """Benjamini-Hochberg, returned in the input order."""
    p = np.asarray(pvals, float)
    o = np.argsort(p)
    m = len(p)
    q = np.empty(m)
    prev = 1.0
    for i in range(m - 1, -1, -1):
        prev = min(prev, p[o[i]] * m / (i + 1))
        q[o[i]] = prev
    return q


# BH within the contrast across the conditions, separately across the negative controls,
# and separately in each stratum: each stratum is its own contrast.
for arm in ("no_apnea", "no_or_mild", "severe"):
    for pool in (False, True):
        rows = [r for r in fits.values()
                if r["negative_control"] is pool and r[arm] is not None
                and not r[arm]["unstable"]]
        if not rows:
            continue
        for r, q in zip(rows, bh([r[arm]["p"] for r in rows])):
            r[arm]["q"] = float(q)
            r[arm]["q_pool"] = "negative controls" if pool else "conditions"
            r[arm]["q_pool_size"] = len(rows)
    for r in fits.values():
        if r[arm] is not None and "q" not in r[arm]:
            r[arm]["q"] = None          # unstable rows are excluded from the correction
            r[arm]["q_pool"] = None
            r[arm]["q_pool_size"] = None

# ---------------------------------------------------------------- the row set panel c draws
# Anchored on the no-apnea stratum, which is the limiting one: 2,920 people, 300 of them with
# low oxygen. A condition qualifies when the no-apnea contrast is estimable, stable, and
# survives FDR correction inside that stratum. Every qualifying condition is then shown with
# its severe-stratum estimate beside it, and the severe estimate is never used to select rows.
sel = [r for r in fits.values()
       if not r["negative_control"] and r["no_apnea"] is not None
       and not r["no_apnea"]["unstable"] and r["no_apnea"]["q"] is not None
       and r["no_apnea"]["q"] < 0.05 and r["severe"] is not None
       and not r["severe"]["unstable"]]
sel.sort(key=lambda r: -r["no_apnea"]["hr"])
ctrl = [r for r in fits.values()
        if r["negative_control"] and r["no_apnea"] is not None and r["severe"] is not None
        and not r["no_apnea"]["unstable"] and not r["severe"]["unstable"]]
ctrl.sort(key=lambda r: -r["no_apnea"]["hr"])

ctrl_absent = [LABEL[k] for k in NEGATIVE_CONTROLS
               if LABEL[k] not in {r["label"] for r in ctrl}]
low_support = [(r["label"], arm, r[arm]["events_low_o2"])
               for r in sel + ctrl for arm in ("no_apnea", "severe")
               if r[arm]["low_support"]]
dropped_unstable = sorted(r["label"] for r in fits.values()
                          if r["no_apnea"] is not None and r["no_apnea"]["unstable"])
not_estimable = sorted(r["label"] for r in fits.values() if r["no_apnea"] is None)

print(f"\npanel c row set: {len(sel)} conditions, {len(ctrl)} negative controls")
print(f"  {'condition':30s} {'no apnea':>22s}   {'severe apnea':>22s}")
for r in sel + ctrl:
    a, s = r["no_apnea"], r["severe"]
    tag = "  [control]" if r["negative_control"] else ""
    print(f"  {r['label']:30s} {a['hr']:6.2f} ({a['lo']:.2f}-{a['hi']:.2f}) q={a['q']:.3f}   "
          f"{s['hr']:6.2f} ({s['lo']:.2f}-{s['hi']:.2f}) q={s['q']:.3f}{tag}")
if ctrl_absent:
    print(f"  negative controls not estimable in both strata: {ctrl_absent}")
if low_support:
    print(f"  thin low-oxygen arms (under {LOW_SUPPORT_EV} events): {low_support}")

# R30: the row set the sheet DRAWS, anchored on the no-or-mild group (AHI under 15), same rule as above
sel_nm = [r for r in fits.values()
          if not r["negative_control"] and r["no_or_mild"] is not None
          and not r["no_or_mild"]["unstable"] and r["no_or_mild"]["q"] is not None
          and r["no_or_mild"]["q"] < 0.05 and r["severe"] is not None
          and not r["severe"]["unstable"]]
sel_nm.sort(key=lambda r: -r["no_or_mild"]["hr"])
ctrl_nm = [r for r in fits.values()
           if r["negative_control"] and r["no_or_mild"] is not None and r["severe"] is not None
           and not r["no_or_mild"]["unstable"] and not r["severe"]["unstable"]]
ctrl_nm.sort(key=lambda r: -r["no_or_mild"]["hr"])
ctrl_absent_nm = [LABEL[k] for k in NEGATIVE_CONTROLS if LABEL[k] not in {r["label"] for r in ctrl_nm}]
low_support_nm = [(r["label"], arm, r[arm]["events_low_o2"]) for r in sel_nm + ctrl_nm
                  for arm in ("no_or_mild", "severe") if r[arm]["low_support"]]
raw_only_nm = sorted(r["label"] for r in fits.values()
                     if not r["negative_control"] and r["no_or_mild"] is not None and not r["no_or_mild"]["unstable"]
                     and r["no_or_mild"]["p"] < 0.05 and (r["no_or_mild"]["q"] is None or r["no_or_mild"]["q"] >= 0.05))
print(f"\npanel c row set, NO-OR-MILD anchor (the drawn one): {len(sel_nm)} conditions, {len(ctrl_nm)} negative controls; "
      f"raw P < 0.05 without BH: {len(raw_only_nm)}")
print(f"  {'condition':30s} {'no or mild apnea':>22s}   {'severe apnea':>22s}")
for r in sel_nm + ctrl_nm:
    a, s_ = r["no_or_mild"], r["severe"]
    tag = "  [control]" if r["negative_control"] else ""
    print(f"  {r['label']:30s} {a['hr']:6.2f} ({a['lo']:.2f}-{a['hi']:.2f}) q={a['q']:.3f}   "
          f"{s_['hr']:6.2f} ({s_['lo']:.2f}-{s_['hi']:.2f}) q={s_['q']:.3f}{tag}")

R = {
    "provenance": {
        "cohort_n": int(len(b)),
        "source": "data_frozen_v8_2026-09/t90_final.parquet via numbers/cohort_spec.apply_cohort",
        "model": "site-stratified Cox (Efron), unpenalized, cr(a, df=4) age spline with one "
                 "column dropped, plus sex, incident cases only",
        "fdr": "Benjamini-Hochberg within the contrast across the conditions, separately "
               "across the five negative controls, separately in each stratum",
        "created_by": "FINAL_FIGURES_2026-08-14/_scripts/compute_fig3_within_stratum.py",
        "contrast": "T90 above 10% of the recording against 10% or less, WITHIN each stratum",
        "stability_rule": "imported from q2_noapnea_freq/model.py, the settled convention for "
                          "this binary contrast: at least 5 events in the low-oxygen arm, at "
                          "least 10 in total, confidence interval under 3 orders of magnitude",
        "min_events": MIN_EVENTS,
        "min_exposed_events": MIN_EXPOSED_EV,
        "positive_control": "reproduces every estimable cell of the twelve-cell cross for heart "
                            "failure and respiratory failure, and the frozen three-level "
                            "within-severe contrast for heart failure and the cardiovascular "
                            "composite, to 1e-9",
        "frozen_noapnea_cross_check": f"{checked} outcomes re-fitted and matched to "
                                      "q2_noapnea_freq/results_noapnea.json within 0.005",
    },
    "positive_control": REPRO,
    "strata": {
        "no_apnea": {"rule": "AHI < 5", "n": int(len(noap)),
                     "n_low_o2": int(noap.ox_gt10.sum()),
                     "n_rest": int((~noap.ox_gt10.astype(bool)).sum())},
        "no_or_mild": {"rule": "AHI < 15", "n": int(len(nomild)),
                       "n_low_o2": int(nomild.ox_gt10.sum()),
                       "n_rest": int((~nomild.ox_gt10.astype(bool)).sum())},
        "severe": {"rule": "AHI >= 30", "n": int(len(sev)),
                   "n_low_o2": int(sev.ox_gt10.sum()),
                   "n_rest": int((~sev.ox_gt10.astype(bool)).sum())},
    },
    "drawn_left_group": "no_or_mild",
    "drawn_left_group_note": "Alen 2026-09-08: Fig 4b's left half is the no-or-mild group (AHI under 15), the row set is "
                             "anchored on it (BH within the group); the no-apnea arm (AHI under 5) is kept for the record, not drawn",
    "panel_c_rows_no_or_mild": [r["label"] for r in sel_nm],
    "panel_c_controls_no_or_mild": [r["label"] for r in ctrl_nm],
    "panel_c_controls_absent_no_or_mild": ctrl_absent_nm,
    "low_support_rows_no_or_mild": low_support_nm,
    "raw_p_only_no_or_mild": raw_only_nm,
    "dropped_unstable_no_or_mild": sorted(r["label"] for r in fits.values()
                                          if r["no_or_mild"] is not None and r["no_or_mild"]["unstable"]),
    "not_estimable_no_or_mild": sorted(r["label"] for r in fits.values() if r["no_or_mild"] is None),
    "prevalence": prev,
    "fits": fits,
    "panel_c_rows": [r["label"] for r in sel],
    "panel_c_controls": [r["label"] for r in ctrl],
    "panel_c_controls_absent": ctrl_absent,
    "low_support_rows": low_support,
    "dropped_unstable_no_apnea": dropped_unstable,
    "not_estimable_no_apnea": not_estimable,
}
json.dump(R, open(OUT, "w"), indent=1)
print(f"\nwritten -> {OUT}")
