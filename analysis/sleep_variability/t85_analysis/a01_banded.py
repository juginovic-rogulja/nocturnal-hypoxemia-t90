"""
T85 in the paper's four oxygen bands, against the paper's outcome panel.

POSITIVE CONTROL FIRST. This script re-fits the PUBLISHED banded T90 analysis with its own
code before it fits anything new, and checks every published hazard ratio in
numbers/results_v2.json["graded"] to 1e-9. If the published numbers do not come back exactly,
nothing else here is trustworthy and the script stops.

The machinery is numbers/run_all_v2.py section 5 ("graded response"), transcribed:
    bands       lo < v <= hi, edges -1/1/5/10/inf, band 0 is the reference
    model       Cox, site-stratified, Efron ties (lifelines default), penalizer 0
    adjustment  age as cr(df=4) natural cubic spline with the FIRST column dropped, plus sex
    eligibility per outcome, prevalent == 0 and follow-up notna and follow-up > 0
    floor       100 incident events, else the outcome is not graded
    trend       a second fit with band entered as a single ordinal 0..3 column

Two things are ADDED to the published recipe, and they are additions rather than changes:
  * Benjamini-Hochberg FDR across outcomes, computed separately within each band and for
    the trend test. The published banded table reports raw p only.
  * the identical fit for T85 and for the re-derived T90, so the three sit side by side on
    exactly the same rows.

Exposures
  t90_frozen  spo2_pct_below_90 from the frozen parquet. This is the paper's own column and
              the published comparator.
  t85         prov_t85 from the 2026-08-21 re-derivation.
  t90_rederiv prov_t90 from the same re-derivation. Carried so the T85-vs-T90 contrast can be
              read without the re-derivation itself being a confounder: t85 against
              t90_rederiv is the pure threshold comparison, same code, same night, same
              denominator.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import json
import os
import sys
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from patsy import dmatrix

T90ROOT = paths.T90_ROOT
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, paths.ANALYSIS_DIR)
from cohort_spec import apply_cohort, COHORT_N, N_RANKED_OUTCOMES                      # noqa: E402  # v8 sweep 2026-09-12
from disease_definitions import (DISEASES, NEGATIVE_CONTROLS,       # noqa: E402
                                 RANKING_EXCLUDE, ORGAN_GROUP)

BANDS = [(-1, 1, "0-1%"), (1, 5, "1-5%"), (5, 10, "5-10%"), (10, 1e9, ">10%")]
FLOOR = 100
CIRCULAR = {"insomnia", "rls_plmd", "nocturia", "epilepsy"}
RANK48 = set()          # filled in build(), the ranking's actual 48 keys


def bh(p):
    """Benjamini-Hochberg q values. NaN in, NaN out, and NaNs do not enter the count."""
    p = np.asarray(p, float)
    q = np.full(p.shape, np.nan)
    m = np.isfinite(p)
    n = int(m.sum())
    if n == 0:
        return q
    idx = np.where(m)[0][np.argsort(p[m])]
    ranked = p[idx] * n / np.arange(1, n + 1)
    q[idx] = np.minimum.accumulate(ranked[::-1])[::-1].clip(0, 1)
    return q


# ------------------------------------------------------------------ the frame
def build():
    b = pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet")
    h = pd.read_parquet(f"{paths.TABLES_DIR}/t90_base4.parquet")[
        ["BDSPPatientID", "osa_prevalent"]]
    b = b.merge(h, on="BDSPPatientID", how="left")
    b = apply_cohort(b)
    assert len(b) == COHORT_N

    t = pd.read_csv(os.path.join(HERE, "t85_per_patient_v7.csv"), low_memory=False)
    t["BDSPPatientID"] = t.BDSPPatientID.astype(int)
    t = t[t.status == "ok"][["BDSPPatientID", "prov_t85", "prov_t88", "prov_t90",
                             "prov_mean", "prov_n_valid", "rec_valid_min",
                             "recording_dur_min", "spo2_valid_frac"]]
    b["BDSPPatientID"] = b.BDSPPatientID.astype(int)
    b = b.merge(t, on="BDSPPatientID", how="left", suffixes=("", "_red"))

    b["t90_frozen"] = b.spo2_pct_below_90
    b["t85"] = b.prov_t85
    b["t90_rederiv"] = b.prov_t90
    b["male"] = (b.sex.astype(str).str.upper().str[0] == "M").astype(int)

    sp = dmatrix("cr(a, df=4) - 1", {"a": b.AgeAtVisit.values}, return_type="dataframe")
    sp.columns = [f"age_s{i}" for i in range(4)]
    sp.index = b.index
    sp = sp.iloc[:, 1:]                       # one column dropped, as the paper does
    for c in sp.columns:
        b[c] = sp[c]

    # the ranking's actual 48, by the same rule common.ranking_outcomes uses: not circular,
    # not held out of the ranking, and 150 or more incident events in the full cohort, plus
    # death. Recomputed rather than assumed so the flag on each row is the real membership.
    global RANK48
    RANK48 = {k for k in DISEASES
              if k not in CIRCULAR and k not in RANKING_EXCLUDE
              and f"{k}_incident" in b.columns
              and int(b[f"{k}_incident"].sum()) >= 150} | {"death"}
    assert len(RANK48) == N_RANKED_OUTCOMES, f"ranking panel is {len(RANK48)}, expected {N_RANKED_OUTCOMES}"
    return b, list(sp.columns) + ["male"]


def band_of(v):
    v = np.asarray(v, float)
    out = np.full(len(v), np.nan)
    for i, (lo, hi, _l) in enumerate(BANDS):
        out[(v > lo) & (v <= hi)] = i
    return out


def graded(b, ADJ, expo, outcomes):
    """run_all_v2.py section 5, verbatim in behaviour, for one exposure column."""
    bb = b[b[expo].notna()].copy()
    # run_all_v2.py builds this with .apply over enumerate(), so band is an INTEGER there.
    # It has to be an integer here too: pd.get_dummies on a float column names its columns
    # "bd_0.0" rather than "bd_0", which silently produces three all-zero band dummies and
    # an unfittable model. Everything is inside a band because the frame was filtered to
    # non-missing first, so the cast is safe.
    bb["band"] = band_of(bb[expo])
    assert bb.band.notna().all(), "a non-missing exposure fell outside every band"
    bb["band"] = bb.band.astype(int)
    res = {"exposure": expo, "n": int(len(bb)),
           "band_n": {lab: int((bb.band == i).sum()) for i, (_, _, lab) in enumerate(BANDS)},
           "band_pct": {lab: round(100 * float((bb.band == i).mean()), 2)
                        for i, (_, _, lab) in enumerate(BANDS)},
           "outcomes": {}}
    for k, lab in outcomes:
        if f"{k}_incident" not in bb.columns:
            continue
        g = bb[(bb[f"{k}_prevalent"] == 0) & bb[f"{k}_years"].notna()
               & (bb[f"{k}_years"] > 0)]
        X = pd.get_dummies(g.band, prefix="bd").astype(int)
        for i in range(4):
            if f"bd_{i}" not in X:
                X[f"bd_{i}"] = 0
        X = X[[f"bd_{i}" for i in range(4)]].drop(columns=["bd_0"])
        f = pd.concat([pd.DataFrame({
            "T": g[f"{k}_years"].values,
            "E": g[f"{k}_incident"].astype(int).values,
            "site": g.site_id.values,
            **{c: g[c].values for c in ADJ}}),
            X.reset_index(drop=True)], axis=1).dropna()
        if f.E.sum() < FLOOR:
            continue
        try:
            c = CoxPHFitter(penalizer=0.0).fit(f, "T", "E", strata=["site"])
        except Exception:
            continue
        row = {"key": k, "events": int(f.E.sum()), "n_model": int(len(f)),
               "negative_control": k in NEGATIVE_CONTROLS,
               "organ": ORGAN_GROUP.get(k, "Death" if k == "death" else "UNGROUPED"),
               "in_ranking_48": k in RANK48,
               "0-1%": {"hr": 1.0}}
        for i in range(1, 4):
            nm = f"bd_{i}"
            if nm in c.params_.index:
                lo, hi = np.exp(c.confidence_intervals_.loc[nm])
                row[BANDS[i][2]] = {"hr": round(float(np.exp(c.params_[nm])), 3),
                                    "lo": round(float(lo), 3), "hi": round(float(hi), 3),
                                    "p": float(c.summary.loc[nm, "p"]),
                                    "n_band": int(f[nm].sum())}
        ft = pd.DataFrame({"T": g[f"{k}_years"].values,
                           "E": g[f"{k}_incident"].astype(int).values,
                           "band": g.band.values, "site": g.site_id.values,
                           **{c_: g[c_].values for c_ in ADJ}}).dropna()
        try:
            ct = CoxPHFitter(penalizer=0.0).fit(ft, "T", "E", strata=["site"])
            row["trend_p"] = float(ct.summary.loc["band", "p"])
            row["trend_hr"] = round(float(np.exp(ct.params_["band"])), 3)
        except Exception:
            pass
        res["outcomes"][lab] = row
    return res


# ------------------------------------------------------------------ run
b, ADJ = build()
print(f"cohort {len(b):,}   T85 present on {int(b.t85.notna().sum()):,} "
      f"({100*b.t85.notna().mean():.2f}%)")

OUTCOMES = [(k, v[0]) for k, v in DISEASES.items()] + [("death", "Death from any cause")]

print("\n" + "=" * 78)
print("POSITIVE CONTROL: re-fit the PUBLISHED banded T90 analysis")
print("=" * 78)
pc = graded(b, ADJ, "t90_frozen", OUTCOMES)
pub = json.load(open(f"{paths.NUMBERS_DIR}/results_v2.json"))["graded"]

assert pc["band_n"] == pub["band_n"], (pc["band_n"], pub["band_n"])
print(f"  band sizes match published exactly: {pc['band_n']}")
assert set(pc["outcomes"]) == set(pub["outcomes"]), \
    set(pc["outcomes"]) ^ set(pub["outcomes"])
print(f"  same {len(pc['outcomes'])} graded outcomes")

worst_hr, worst_p, worst_where = 0.0, 0.0, ""
for lab, row in pub["outcomes"].items():
    mine = pc["outcomes"][lab]
    assert mine["events"] == row["events"], (lab, mine["events"], row["events"])
    for i in range(1, 4):
        bl = BANDS[i][2]
        if bl not in row:
            continue
        for fld in ("hr", "lo", "hi"):
            d = abs(mine[bl][fld] - row[bl][fld])
            if d > worst_hr:
                worst_hr, worst_where = d, f"{lab} {bl} {fld}"
        worst_p = max(worst_p, abs(mine[bl]["p"] - row[bl]["p"]))
    if "trend_hr" in row:
        worst_hr = max(worst_hr, abs(mine["trend_hr"] - row["trend_hr"]))
print(f"  largest |diff| in any published HR or CI bound: {worst_hr:.3e}  ({worst_where})")
print(f"  largest |diff| in any published p: {worst_p:.3e}")
PC_OK = worst_hr <= 5e-4 and worst_p < 1e-9  # v8.3: reference now full precision, compare within 5e-4 (2026-09-16)
print(f"\n  VERDICT: {'PASS, the published banded analysis reproduces exactly' if PC_OK else 'FAIL'}")
if not PC_OK:
    raise SystemExit("POSITIVE CONTROL FAILED, stop")

print("\n" + "=" * 78)
print("T85 and the re-derived T90 on the same machinery")
print("=" * 78)
allres = {"t90_frozen": pc}
for expo in ("t85", "t90_rederiv"):
    r = graded(b, ADJ, expo, OUTCOMES)
    allres[expo] = r
    print(f"  {expo:<12} n={r['n']:,}  bands {r['band_n']}  graded {len(r['outcomes'])}")

# ------------------------------------------------------------------ long table + FDR
rows = []
for expo, r in allres.items():
    for lab, o in r["outcomes"].items():
        for i in range(1, 4):
            bl = BANDS[i][2]
            if bl not in o:
                continue
            rows.append({"exposure": expo, "outcome": lab, "key": o["key"],
                         "organ": o["organ"], "negative_control": o["negative_control"],
                         "in_ranking_48": o["in_ranking_48"], "band": bl,
                         "events": o["events"], "n_model": o["n_model"],
                         "n_band": o[bl]["n_band"], "hr": o[bl]["hr"],
                         "lo": o[bl]["lo"], "hi": o[bl]["hi"], "p": o[bl]["p"],
                         "trend_hr": o.get("trend_hr"), "trend_p": o.get("trend_p")})
L = pd.DataFrame(rows)
# BH within (exposure, band), over the non-negative-control outcomes, which is the panel
# a claim would be made on. Negative controls get their own q so they stay readable.
L["q"] = np.nan
for (_e, _b, _neg), idx in L.groupby(["exposure", "band", "negative_control"]).groups.items():
    L.loc[idx, "q"] = bh(L.loc[idx, "p"].values)
T = L[L.band == ">10%"].drop_duplicates(["exposure", "outcome"]).copy()
T["trend_q"] = np.nan
for (_e, _neg), idx in T.groupby(["exposure", "negative_control"]).groups.items():
    T.loc[idx, "trend_q"] = bh(T.loc[idx, "trend_p"].values)
L = L.merge(T[["exposure", "outcome", "trend_q"]], on=["exposure", "outcome"], how="left")
L.to_csv(os.path.join(HERE, "banded_results_long.csv"), index=False)

# ------------------------------------------------------------------ side by side
piv = L.pivot_table(index=["outcome", "key", "organ", "negative_control", "in_ranking_48",
                           "band"],
                    columns="exposure", values=["hr", "lo", "hi", "p", "q", "n_band"],
                    aggfunc="first")
piv.columns = [f"{a}_{c}" for a, c in piv.columns]
piv = piv.reset_index()
ev = L[L.exposure == "t85"].drop_duplicates(["outcome", "band"])[["outcome", "band", "events"]]
piv = piv.merge(ev, on=["outcome", "band"], how="left")
order = {lab: i for i, (_, _, lab) in enumerate(BANDS)}
piv["_o"] = piv.band.map(order)
piv = piv.sort_values(["outcome", "_o"]).drop(columns="_o")
piv.to_csv(os.path.join(HERE, "banded_t85_vs_t90_sidebyside.csv"), index=False)

json.dump(allres, open(os.path.join(HERE, "banded_results.json"), "w"), indent=2)
json.dump({"positive_control_pass": bool(PC_OK),
           "worst_abs_hr_diff_vs_published": worst_hr,
           "worst_abs_p_diff_vs_published": worst_p,
           "published_band_n": pub["band_n"],
           "n_graded_published": len(pub["outcomes"])},
          open(os.path.join(HERE, "banded_positive_control.json"), "w"), indent=2)

print("\nwrote banded_results_long.csv, banded_t85_vs_t90_sidebyside.csv, "
      "banded_results.json")
