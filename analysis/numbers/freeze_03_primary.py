"""
Phase 1.3 to 1.11. Recompute every primary number the manuscript will state.

Nothing here is copied from an earlier working file. Every value is recomputed from the source
parquet so that the manuscript quotes a number that was produced by code sitting in the
repository, not a number transcribed from a conversation.

Writes numbers/primary.json.

Model specification, applied identically throughout unless stated:
  Cox proportional hazards, stratified by hospital, with a natural cubic spline on age (4 df)
  and sex. Exposure entered as a rank-based inverse-normal score computed within hospital, so
  hazard ratios are per 1 SD and are comparable across measures on different scales.
  Prevalent cases excluded per outcome. Follow-up must be positive.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../..")))  # repository root, where paths.py lives
import paths

# ----------------------------------------------------------------------------- superseded specification
# This script fits the primary associations at penalizer=0.01, the specification that run_primary_unpenalized.py (step 111) replaced: the printed values
# come from run_primary_unpenalized.py (step 111) (unpenalized), and nothing in the repository reads numbers/primary.json. It is kept for the
# old-versus-new record of the v8 recalculation and stops here unless a reviewer asks for it explicitly.
import os as _guard_os
if _guard_os.environ.get("T90_RUN_SUPERSEDED") != "1":
    raise SystemExit("freeze_03_primary: superseded penalized specification, not a source of any printed number. "
                     "Set T90_RUN_SUPERSEDED=1 to run it anyway.")
import warnings
warnings.filterwarnings("ignore")
import json
import numpy as np
import pandas as pd
from scipy import stats
from patsy import dmatrix
from lifelines import CoxPHFitter

OUT = paths.NUMBERS_DIR
R = {}

b = pd.read_parquet(f"{paths.TABLES_DIR}/t90_base4.parquet")
b = b[(b.fu_valid == 1) & (b.oximetry_bad == 0)].copy()
b["t90"] = b.spo2_pct_below_90
b = b[b.t90.notna()]
b["male"] = (b.sex.astype(str).str.upper().str[0] == "M").astype(int)

sp = dmatrix("cr(a, df=4) - 1", {"a": b.AgeAtVisit.values}, return_type="dataframe")
sp.columns = [f"age_s{i}" for i in range(sp.shape[1])]
sp.index = b.index
for c in sp.columns:
    b[c] = sp[c]
ADJ = [f"age_s{i}" for i in range(4)] + ["male"]


def rint(df, col):
    z = pd.Series(np.nan, index=df.index)
    for s, idx in df.groupby("site_id").groups.items():
        v = df.loc[idx, col]
        ok = v.notna()
        if ok.sum() < 20:
            continue
        r = stats.rankdata(v[ok], method="average")
        z.loc[v[ok].index] = stats.norm.ppf((r - 0.375) / (ok.sum() + 0.25))
    return z


b["t90_z"] = rint(b, "t90")
if "AHI" in b.columns:
    b["ahi_z"] = rint(b, "AHI")

import sys as _sys_dd; _sys_dd.path.insert(0, paths.ANALYSIS_DIR)
from disease_definitions import DISEASES, NEGATIVE_CONTROLS, RANKING_EXCLUDE   # v8.1c
DIS = [(k, v[0]) for k, v in DISEASES.items()] + [("death", "Death from any cause")]   # v8.1c: every condition of the definitions (57) plus death
# The panel is imported, never listed here. Settled 2026-08-07: back pain, cataract,
# glaucoma, contact dermatitis, hemorrhoids. Fracture and osteoarthritis stay as ordinary
# outcomes. See numbers/negcontrols_final.json.
import sys as _sys
_sys.path.insert(0, paths.ANALYSIS_DIR)
from disease_definitions import NEGATIVE_CONTROLS as _NC
NEG = set(NEGATIVE_CONTROLS)   # v8.1c


def cox(df, xcols, k, strat=True):
    inc, yrs, prev = f"{k}_incident", f"{k}_years", f"{k}_prevalent"
    if inc not in df.columns:
        return None
    g = df[(df[prev] == 0) & df[yrs].notna() & (df[yrs] > 0)]
    cols = {"T": g[yrs].values, "E": g[inc].astype(int).values}
    for c in xcols + ADJ:
        cols[c] = g[c].values
    if strat:
        cols["site"] = g.site_id.values
    f = pd.DataFrame(cols).dropna()
    if f.E.sum() < 40:
        return None
    try:
        c = CoxPHFitter(penalizer=0.01).fit(f, "T", "E",
                                            strata=["site"] if strat else None)
    except Exception:
        return None
    out = {"n": int(len(f)), "events": int(f.E.sum())}
    for x in xcols:
        lo, hi = np.exp(c.confidence_intervals_.loc[x])
        out[x] = {"hr": round(float(np.exp(c.params_[x])), 3),
                  "lo": round(float(lo), 3), "hi": round(float(hi), 3),
                  "p": float(c.summary.loc[x, "p"])}
    return out


# ---------------------------------------------------------------- 1.4 T90 vs AHI, 52 diseases
print("1.4  T90 and AHI across diseases")
rows = []
for k, lab in DIS:
    a = cox(b, ["t90_z"], k)
    if a is None:
        continue
    row = {"key": k, "disease": lab, "events": a["events"], "n": a["n"],
           "t90_hr": a["t90_z"]["hr"], "t90_lo": a["t90_z"]["lo"],
           "t90_hi": a["t90_z"]["hi"], "t90_p": a["t90_z"]["p"],
           "negative_control": k in NEG}
    if "ahi_z" in b.columns:
        c2 = cox(b, ["ahi_z"], k)
        if c2:
            row.update({"ahi_hr": c2["ahi_z"]["hr"], "ahi_lo": c2["ahi_z"]["lo"],
                        "ahi_hi": c2["ahi_z"]["hi"], "ahi_p": c2["ahi_z"]["p"]})
        j = cox(b, ["t90_z", "ahi_z"], k)
        if j:
            row.update({"t90_adj_ahi_hr": j["t90_z"]["hr"], "t90_adj_ahi_p": j["t90_z"]["p"],
                        "ahi_adj_t90_hr": j["ahi_z"]["hr"], "ahi_adj_t90_p": j["ahi_z"]["p"]})
    rows.append(row)
dis = pd.DataFrame(rows).sort_values("t90_hr", ascending=False)
dis.to_csv(f"{OUT}/bdsp_diseases.csv", index=False)
R["n_diseases_tested"] = int(len(dis))
R["n_diseases_significant"] = int((dis.t90_p < .05).sum())
R["diseases_top10"] = dis.head(10)[["disease", "events", "t90_hr", "t90_lo",
                                    "t90_hi", "t90_p"]].to_dict("records")
print(f"     {len(dis)} diseases, {int((dis.t90_p<.05).sum())} significant")

# ---------------------------------------------------------------- 1.11 healthy subset
print("1.11 healthy subset")
FREE = ["hf", "cvd", "ihd", "mi", "afib", "stroke_any", "htn2", "diabetes", "ckd", "copd2",
        "asthma", "cancer_any", "dementia", "cirrhosis", "obesity_hypovent"]
mask = np.ones(len(b), dtype=bool)
used = []
for c in FREE:
    if f"{c}_prevalent" in b.columns:
        mask &= (b[f"{c}_prevalent"] == 0).values
        used.append(c)
h = b[mask].copy()
h["hi10"] = (h.t90 > 10).astype(int)
b["hi10"] = (b.t90 > 10).astype(int)
R["healthy"] = {"n": int(len(h)), "pct_of_cohort": round(100 * len(h) / len(b), 1),
                "n_conditions_screened": len(used),
                "age_median": float(h.AgeAtVisit.median()),
                "male_pct": round(100 * h.male.mean(), 1),
                "t90_median": round(float(h.t90.median()), 2),
                "t90_above10_pct": round(100 * float((h.t90 > 10).mean()), 1),
                "results": {}}
for k, lab in ([("hf", "Heart failure"), ("copd2", "COPD"),
               ("resp_failure", "Respiratory failure"), ("aki", "Acute kidney injury"),
               ("cvd", "Cardiovascular composite"), ("diabetes", "Type 2 diabetes"),
               ("death", "Death from any cause"), ("osteoarthritis", "Osteoarthritis")]
              + [(k, DISEASES[k][0]) for k in NEGATIVE_CONTROLS if k not in RANKING_EXCLUDE]):   # v8.1c: the ranked controls (glaucoma, inguinal hernia, alopecia)
    a = cox(h, ["hi10"], k)
    c = cox(b, ["hi10"], k)
    if a:
        R["healthy"]["results"][lab] = {
            "healthy": {"events": a["events"], **a["hi10"]},
            "full_cohort": ({"events": c["events"], **c["hi10"]} if c else None)}
print(f"     n={len(h):,} ({R['healthy']['pct_of_cohort']}% of cohort)")

with open(f"{OUT}/primary.json", "w") as f:
    json.dump(R, f, indent=2)
print(f"\nwritten -> {OUT}/primary.json and bdsp_diseases.csv")
