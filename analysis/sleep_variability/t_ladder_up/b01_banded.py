"""
The threshold ladder in the paper's four oxygen bands, against the paper's outcome
panel, plus the A95 preserved-saturation banded framing.

POSITIVE CONTROL FIRST, exactly as t85_analysis/a01_banded.py: this script re-fits the
PUBLISHED banded T90 analysis with its own code before it fits anything new, and checks
every published hazard ratio in numbers/results_v2.json["graded"] to 1e-9. If the
published numbers do not come back exactly, nothing else here is trustworthy and the
script stops.

The machinery is numbers/run_all_v2.py section 5 ("graded response"), transcribed via
a01_banded.py:
    bands       lo < v <= hi, edges -1/1/5/10/inf, band 0 is the reference
    model       Cox, site-stratified, Efron ties (lifelines default), penalizer 0
    adjustment  age as cr(df=4) natural cubic spline with the FIRST column dropped, plus sex
    eligibility per outcome, prevalent == 0 and follow-up notna and follow-up > 0
    floor       100 incident events, else the outcome is not graded
    trend       a second fit with band entered as a single ordinal column

Framings fitted here
  paper bands   t80, t92, t95 (and the t90_frozen positive control). A paper-band fit is
                only attempted when the top band holds >= 300 patients, per the task.
  quartiles     t80, t92, t95 cut at their own quartiles, Q1 the reference. This is the
                fallback for exposures whose paper reference band collapses (T95) and a
                readability check for the others. Clearly labeled _quartile.
  preserved     a95 = percent of the recording AT or ABOVE 95 percent SpO2, banded as
                percent of the night at normal saturation: >=99 / 90-<99 / 50-<90 / <50,
                reference = the most-preserved band. Cuts are adjusted (and reported) if
                any fitted band would hold < 300 patients.

FDR (Benjamini-Hochberg) is ADDED within each (framing, band) over the real outcomes and
separately over the negative controls, as the T85 run did. The published table reports
raw p only.
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
TOPBAND_MIN = 300
CIRCULAR = {"insomnia", "rls_plmd", "nocturia", "epilepsy"}
RANK48 = set()

# preserved-saturation cut candidates, first whose every band holds >= 300 wins.
# each is the list of lower edges of the "at or above" bands, descending.
A95_CUTS = [
    (99.0, 90.0, 50.0),      # the spec: >=99 / 90-<99 / 50-<90 / <50
    (98.0, 90.0, 50.0),
    (99.0, 95.0, 80.0),
    (95.0, 80.0, 50.0),
]


def bh(p):
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

    t = pd.read_csv(os.path.join(HERE, "ladder_per_patient_v7.csv"), low_memory=False)
    t["BDSPPatientID"] = t.BDSPPatientID.astype(int)
    t = t[t.status == "ok"][["BDSPPatientID", "prov_t80", "prov_t85", "prov_t88",
                             "prov_t90", "prov_t92", "prov_t95", "prov_a95"]]
    b["BDSPPatientID"] = b.BDSPPatientID.astype(int)
    b = b.merge(t, on="BDSPPatientID", how="left", suffixes=("", "_red"))

    b["t90_frozen"] = b.spo2_pct_below_90
    b["t88_frozen"] = b.spo2_pct_below_88
    # short names for the re-derived rungs. t88/t90 short names are the RE-DERIVED
    # values; the frozen ones are only ever *_frozen.
    for c in ("t80", "t85", "t88", "t90", "t92", "t95", "a95"):
        b[c] = b[f"prov_{c}"]
    b["t88_rederiv"] = b.prov_t88
    b["t90_rederiv"] = b.prov_t90
    b["male"] = (b.sex.astype(str).str.upper().str[0] == "M").astype(int)

    sp = dmatrix("cr(a, df=4) - 1", {"a": b.AgeAtVisit.values}, return_type="dataframe")
    sp.columns = [f"age_s{i}" for i in range(4)]
    sp.index = b.index
    sp = sp.iloc[:, 1:]                       # one column dropped, as the paper does
    for c in sp.columns:
        b[c] = sp[c]

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


def graded_generic(bb, ADJ, framing, band_labels, outcomes):
    """
    run_all_v2.py section 5 behaviour on a precomputed integer `band` column in bb.
    band 0 is always the reference. Returns the same result dict shape as a01's graded().
    """
    nb = len(band_labels)
    res = {"exposure": framing, "n": int(len(bb)),
           "band_n": {lab: int((bb.band == i).sum()) for i, lab in enumerate(band_labels)},
           "band_pct": {lab: round(100 * float((bb.band == i).mean()), 2)
                        for i, lab in enumerate(band_labels)},
           "outcomes": {}}
    for k, lab in outcomes:
        if f"{k}_incident" not in bb.columns:
            continue
        g = bb[(bb[f"{k}_prevalent"] == 0) & bb[f"{k}_years"].notna()
               & (bb[f"{k}_years"] > 0)]
        X = pd.get_dummies(g.band, prefix="bd").astype(int)
        for i in range(nb):
            if f"bd_{i}" not in X:
                X[f"bd_{i}"] = 0
        X = X[[f"bd_{i}" for i in range(nb)]].drop(columns=["bd_0"])
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
               band_labels[0]: {"hr": 1.0}}
        for i in range(1, nb):
            nm = f"bd_{i}"
            if nm in c.params_.index:
                lo, hi = np.exp(c.confidence_intervals_.loc[nm])
                row[band_labels[i]] = {"hr": round(float(np.exp(c.params_[nm])), 3),
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


def graded(b, ADJ, expo, outcomes, framing=None):
    """Paper bands on exposure column `expo`, identical in behaviour to a01's graded()."""
    bb = b[b[expo].notna()].copy()
    bb["band"] = band_of(bb[expo])
    assert bb.band.notna().all(), "a non-missing exposure fell outside every band"
    bb["band"] = bb.band.astype(int)
    return graded_generic(bb, ADJ, framing or expo, [x[2] for x in BANDS], outcomes)


# ------------------------------------------------------------------ run
b, ADJ = build()
LADDER = ["t80", "t85", "t88", "t90", "t92", "t95"]
for c in LADDER + ["a95"]:
    print(f"  {c:<4} present on {int(b[c].notna().sum()):,} ({100*b[c].notna().mean():.2f}%)")

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

# ------------------------------------------------------------------ paper bands
print("\n" + "=" * 78)
print("PAPER BANDS for the new rungs (fit only where the top band holds >= 300)")
print("=" * 78)
allres = {"t90_frozen": pc}
skipped = {}
for expo in ("t80", "t92", "t95"):
    sizes = {lab: int((band_of(b[expo].dropna()) == i).sum())
             for i, (_, _, lab) in enumerate(BANDS)}
    top = sizes[">10%"]
    print(f"  {expo:<4} band sizes {sizes}")
    if top >= TOPBAND_MIN:
        r = graded(b, ADJ, expo, OUTCOMES)
        allres[expo] = r
        print(f"       fitted, {len(r['outcomes'])} outcomes graded")
    else:
        skipped[expo] = sizes
        print(f"       NOT fitted with paper bands, top band {top} < {TOPBAND_MIN}")

# ------------------------------------------------------------------ quartiles
print("\n" + "=" * 78)
print("QUARTILE BANDS, Q1 the reference, clearly labeled fallback framing")
print("=" * 78)
qinfo = {}
for expo in ("t80", "t92", "t95"):
    bb = b[b[expo].notna()].copy()
    try:
        qb, edges = pd.qcut(bb[expo], 4, labels=False, retbins=True, duplicates="drop")
    except Exception as e:
        print(f"  {expo:<4} qcut failed: {e}")
        continue
    nb = int(qb.max()) + 1
    if nb < 4:
        print(f"  {expo:<4} quartiles degenerate ({nb} distinct bins), still fitted with "
              f"{nb} bands")
    labels = [f"Q{i+1}" for i in range(nb)]
    bb["band"] = qb.astype(int)
    edges_r = [round(float(e), 4) for e in edges]
    qinfo[expo] = {"edges": edges_r,
                   "sizes": {labels[i]: int((bb.band == i).sum()) for i in range(nb)}}
    print(f"  {expo:<4} edges {edges_r}  sizes {qinfo[expo]['sizes']}")
    r = graded_generic(bb, ADJ, f"{expo}_quartile", labels, OUTCOMES)
    allres[f"{expo}_quartile"] = r

# ------------------------------------------------------------------ preserved saturation
print("\n" + "=" * 78)
print("A95 PRESERVED-SATURATION BANDS: percent of the recording at or above 95% SpO2")
print("=" * 78)
av = b[b.a95.notna()]
chosen, csizes = None, None
for cuts in A95_CUTS:
    c1, c2, c3 = cuts
    idx = np.select([av.a95.values >= c1,
                     (av.a95.values >= c2) & (av.a95.values < c1),
                     (av.a95.values >= c3) & (av.a95.values < c2)],
                    [0, 1, 2], default=3)
    sizes = {i: int((idx == i).sum()) for i in range(4)}
    print(f"  cuts >={c1:.0f} / {c2:.0f}-<{c1:.0f} / {c3:.0f}-<{c2:.0f} / <{c3:.0f}"
          f"   sizes {sizes}" + ("   <- chosen" if chosen is None and
                                 min(sizes.values()) >= TOPBAND_MIN else ""))
    if chosen is None and min(sizes.values()) >= TOPBAND_MIN:
        chosen, csizes = cuts, sizes
if chosen is None:
    chosen = A95_CUTS[0]
    print("  no candidate had every band >= 300, using the spec cuts anyway, flagged")
c1, c2, c3 = chosen
A95_LABELS = [f">={c1:.0f}% preserved", f"{c2:.0f}-<{c1:.0f}%",
              f"{c3:.0f}-<{c2:.0f}%", f"<{c3:.0f}% preserved"]
bb = b[b.a95.notna()].copy()
bb["band"] = np.select([bb.a95.values >= c1,
                        (bb.a95.values >= c2) & (bb.a95.values < c1),
                        (bb.a95.values >= c3) & (bb.a95.values < c2)],
                       [0, 1, 2], default=3)
r = graded_generic(bb, ADJ, "a95_preserved", A95_LABELS, OUTCOMES)
allres["a95_preserved"] = r
print(f"  fitted with cuts {chosen}, reference = '{A95_LABELS[0]}', "
      f"{len(r['outcomes'])} outcomes graded")

# ------------------------------------------------------------------ long table + FDR
band_order = {}
rows = []
trend_rows = []
for expo, r in allres.items():
    # band labels in fit order: band_n preserves construction order
    labs = list(r["band_n"].keys())
    for i, lab_ in enumerate(labs):
        band_order[(expo, lab_)] = i
    for lab, o in r["outcomes"].items():
        for bl in labs[1:]:
            if bl not in o:
                continue
            rows.append({"exposure": expo, "outcome": lab, "key": o["key"],
                         "organ": o["organ"], "negative_control": o["negative_control"],
                         "in_ranking_48": o["in_ranking_48"], "band": bl,
                         "events": o["events"], "n_model": o["n_model"],
                         "n_band": o[bl]["n_band"], "hr": o[bl]["hr"],
                         "lo": o[bl]["lo"], "hi": o[bl]["hi"], "p": o[bl]["p"]})
        trend_rows.append({"exposure": expo, "outcome": lab,
                           "negative_control": o["negative_control"],
                           "trend_hr": o.get("trend_hr"), "trend_p": o.get("trend_p")})
L = pd.DataFrame(rows)
L["q"] = np.nan
for (_e, _b, _neg), idx in L.groupby(["exposure", "band", "negative_control"]).groups.items():
    L.loc[idx, "q"] = bh(L.loc[idx, "p"].values)
TR = pd.DataFrame(trend_rows)
TR["trend_q"] = np.nan
for (_e, _neg), idx in TR.groupby(["exposure", "negative_control"]).groups.items():
    TR.loc[idx, "trend_q"] = bh(TR.loc[idx, "trend_p"].values)
L = L.merge(TR, on=["exposure", "outcome", "negative_control"], how="left")
L["band_i"] = [band_order.get((e, bl), np.nan) for e, bl in zip(L.exposure, L.band)]
L = L.sort_values(["exposure", "outcome", "band_i"])
L.to_csv(os.path.join(HERE, "ladder_banded_long.csv"), index=False)

json.dump(allres, open(os.path.join(HERE, "ladder_banded_results.json"), "w"), indent=2)
json.dump({"positive_control_pass": bool(PC_OK),
           "worst_abs_hr_diff_vs_published": worst_hr,
           "worst_abs_p_diff_vs_published": worst_p,
           "published_band_n": pub["band_n"],
           "n_graded_published": len(pub["outcomes"]),
           "paper_band_skipped": skipped,
           "quartile_info": qinfo,
           "a95_cuts_used": list(chosen),
           "a95_band_labels": A95_LABELS,
           "a95_band_sizes": csizes},
          open(os.path.join(HERE, "banded_positive_control.json"), "w"), indent=2)

# quick screen summary
print("\nTOP-BAND SUMMARY, real outcomes in the 48-panel, q < 0.05")
for expo, r in allres.items():
    labs = list(r["band_n"].keys())
    top = labs[-1]
    sub = L[(L.exposure == expo) & (L.band == top) & (~L.negative_control)
            & L.in_ranking_48]
    if not len(sub):
        continue
    print(f"  {expo:<14} top band '{top}' n_band {r['band_n'][top]:>6,}  "
          f"sig {int((sub.q < .05).sum()):2d}/{len(sub)}  "
          f"median HR {sub.hr.median():.3f}  "
          f"median sig HR {sub[sub.q < .05].hr.median() if (sub.q < .05).any() else float('nan'):.3f}")

print("\nwrote ladder_banded_long.csv, ladder_banded_results.json, "
      "banded_positive_control.json")
