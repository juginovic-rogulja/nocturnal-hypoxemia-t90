"""
The two replacement negative controls on the categorical oxygen scale.

results_v2.json holds the graded dose-response, oxygen category against rate of new diagnosis,
for all 52 frozen outcomes including the old five controls. Contact dermatitis and hemorrhoids
were coded after that freeze, so they have a per-SD estimate and no categorical one. The merged
eFigure 2 draws every control on the categorical scale, so the two new ones are fitted here on
exactly the specification results_v2 used.

Writes numbers/graded_newcontrols_v1.json.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../..")))  # repository root, where paths.py lives
import paths
import json
import sys
import warnings

import duckdb
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from patsy import dmatrix

warnings.filterwarnings("ignore")
ROOT = paths.FIGURE_ROOT
CACHE = paths.OMOP_CACHE_DIR
sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
from cohort_spec import apply_cohort

b = apply_cohort(pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet"))
b["t90"] = b.spo2_pct_below_90
b["male"] = (b.sex.astype(str).str.upper().str[0] == "M").astype(int)
sp = dmatrix("cr(a, df=4) - 1", {"a": b.AgeAtVisit.values}, return_type="dataframe").iloc[:, 1:]
for i in range(sp.shape[1]):
    b[f"age_s{i}"] = sp.iloc[:, i].values
ADJ = [f"age_s{i}" for i in range(sp.shape[1])] + ["male"]

BANDS = [(-1, 1, "0-1%"), (1, 5, "1-5%"), (5, 10, "5-10%"), (10, 1e9, ">10%")]
b["band"] = b.t90.apply(lambda v: next(i for i, (lo, hi, _l) in enumerate(BANDS) if lo < v <= hi))

NEW = {"dermatitis_contact": ("Contact dermatitis", ["692", "L23", "L24", "L25"]),
       "haemorrhoids": ("Hemorrhoids", ["455", "K64"])}
con = duckdb.connect()
con.execute(f"""create table c as
    select person_id, condition_start_date as dt,
           upper(replace(condition_source_value, '.', '')) as code
    from '{CACHE}/condition_occurrence_cohort.parquet'
    where condition_source_value is not null""")
for k, (_lab, pref) in NEW.items():
    where = " or ".join(f"code like '{p}%'" for p in pref)
    fd = con.execute(f"select person_id, min(dt) as fdt from c where {where} "
                     f"group by person_id").df()
    fd.columns = ["BDSPPatientID", f"{k}_first_date"]
    b = b.drop(columns=[c for c in (f"{k}_first_date",) if c in b.columns])   # v7: the frozen table now carries this column (add_newcontrols_v5), the merge would suffix it
    b = b.merge(fd, on="BDSPPatientID", how="left")
    f0 = pd.to_datetime(b[f"{k}_first_date"], errors="coerce")
    after = f0.notna() & (f0 > b.psg_date)
    b[f"{k}_prevalent"] = (f0.notna() & (f0 <= b.psg_date)).astype(int)
    b[f"{k}_incident"] = after.astype(int)
    end = np.where(after, f0.values, b.censor_date.values)
    b[f"{k}_years"] = (pd.to_datetime(end) - b.psg_date).dt.days / 365.25

OUT = {}
for k, (lab, _p) in NEW.items():
    g = b[(b[f"{k}_prevalent"] == 0) & b[f"{k}_years"].notna() & (b[f"{k}_years"] > 0)]
    X = pd.get_dummies(g.band, prefix="bd").astype(int)
    for i in range(4):
        if f"bd_{i}" not in X:
            X[f"bd_{i}"] = 0
    X = X[[f"bd_{i}" for i in range(4)]].drop(columns=["bd_0"])
    f = pd.concat([pd.DataFrame({"T": g[f"{k}_years"].values,
                                 "E": g[f"{k}_incident"].astype(int).values,
                                 "site": g.site_id.values,
                                 **{c: g[c].values for c in ADJ}}),
                   X.reset_index(drop=True)], axis=1).dropna()
    c = CoxPHFitter(penalizer=0.0).fit(f, "T", "E", strata=["site"])
    rec = {"events": int(f.E.sum()), "negative_control": True, "0-1%": {"hr": 1.0}}
    for i, name in enumerate(["1-5%", "5-10%", ">10%"], start=1):
        lo, hi = np.exp(c.confidence_intervals_.loc[f"bd_{i}"])
        rec[name] = {"hr": float(np.exp(c.params_[f"bd_{i}"])),  # v8.3 full precision (2026-09-16)
                     "lo": float(lo), "hi": float(hi),
                     "p": float(c.summary.loc[f"bd_{i}", "p"])}
    # the trend, oxygen category entered as an ordered score, matching results_v2
    ft = pd.DataFrame({"T": g[f"{k}_years"], "E": g[f"{k}_incident"].astype(int),
                       "site": g.site_id, "bandscore": g.band.astype(float),
                       **{c2: g[c2] for c2 in ADJ}}).dropna()
    ct = CoxPHFitter(penalizer=0.0).fit(ft, "T", "E", strata=["site"])
    rec["trend_hr"] = float(np.exp(ct.params_["bandscore"]))  # v8.3 full precision (2026-09-16)
    rec["trend_p"] = float(ct.summary.loc["bandscore", "p"])
    OUT[lab] = rec
    print(f"{lab:<22}ev {rec['events']:>5}   "
          f"1-5% {rec['1-5%']['hr']:.2f}  5-10% {rec['5-10%']['hr']:.2f}  "
          f">10% {rec['>10%']['hr']:.2f}   trend {rec['trend_hr']:.3f} "
          f"P {rec['trend_p']:.3f}")

json.dump({"_built_by": "New_Figures/_working/scripts/fits_graded_newcontrols.py",
           "_spec": "identical to the graded block of numbers/run_all_v2.py",
           "band_n": {lab: int((b.band == i).sum()) for i, (_a, _c, lab) in enumerate(BANDS)},
           "outcomes": OUT},
          open(f"{paths.NUMBERS_DIR}/graded_newcontrols_v1.json", "w"), indent=1)
print("wrote numbers/graded_newcontrols_v1.json")
