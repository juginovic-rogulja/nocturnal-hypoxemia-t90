#!/usr/bin/env python3
"""Full-precision reference for every value Supp Fig 13 and Supp Fig 14 print from the
non-responder analysis, regenerated here with a provenance sidecar.

numbers/nonresponder_combined_v2.json and numbers/nonresponder_phenotype_v2.csv store
round(x, 3), so a two-decimal print from either is a second rounding. This script refits the
same models, keeps full precision, and reports exactly which printed strings move. The sheet
verifiers read their reference from the file this writes, never from a literal.

Model spec lifted from New_Figures/_NOT_THESE_workfiles/scripts/fits_phenotype_combined_v2.py:
  marginal   sm.Logit(nonresp ~ <condition>_prevalent + AgeAtVisit + male + log1p(pre_sleep_t90))
  joint      the same on the FDR survivors together plus age, sex, race and baseline
  trend      the same on the clipped cardiopulmonary count

Writes work/fullprec_nonresponder.json and work/fullprec_nonresponder.json.provenance.json
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import glob
import hashlib
import json
import os
import subprocess
import sys
import warnings
from datetime import datetime

import duckdb
import numpy as np
import pandas as pd
import statsmodels.api as sm

warnings.filterwarnings("ignore")

ROOT = paths.FIGURE_ROOT
LANE = f"{paths.FIGURE_ROOT}/ROUND50_2026-09-26/lanes/LSUPP-C/L5b_work"   # R49: the twin regenerated in this lane, reconciled against the current v8.3 json
sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
from disease_definitions import DISEASES, NEGATIVE_CONTROLS   # noqa: E402
from rounding_v7 import half_up                               # noqa: E402

CP = ["resp_failure", "copd2", "obesity_hypovent", "hf", "pulm_htn", "asthma"]
STAGES = ["wake", "sleep", "nrem", "n1", "n2", "n3", "rem"]
from cohort_spec import CPAP_STAGE_PARQUET, T90_FINAL   # noqa: E402  V14_L5b: the v8 frozen tables through the spec module
IN_PARQ = [CPAP_STAGE_PARQUET, T90_FINAL]
assert "data_frozen_v8_2026-09" in CPAP_STAGE_PARQUET and "data_frozen_v8_2026-09" in T90_FINAL, IN_PARQ


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


# ------------------------------------------------------------------ the sample, verbatim
d = pd.read_parquet(IN_PARQ[0])
A = d[d.design == "A_split"].copy(); A["src"] = "A"
B = d[d.design == "B_pairs"].copy(); B["src"] = "B"
for s in STAGES:
    for a_, b_ in [("pre", "dx"), ("post", "tx")]:
        for suf in ("t90", "min"):
            c = f"{b_}_{s}_{suf}"
            B[f"{a_}_{s}_{suf}"] = B[c] if c in B.columns else np.nan
keep = (["BDSPPatientID", "SexDSC", "src"]
        + [f"{p}_{s}_{x}" for p in ("pre", "post") for s in STAGES for x in ("t90", "min")])
cc = pd.concat([A[[c for c in keep if c in A.columns]],
                B[[c for c in keep if c in B.columns]]], ignore_index=True)
cc = cc.sort_values("src", ascending=False).drop_duplicates("BDSPPatientID", keep="first")
base = pd.read_parquet(IN_PARQ[1]).drop_duplicates("BDSPPatientID")
m = cc.merge(base, on="BDSPPatientID", how="inner")
m = m[(m.fu_valid == 1) & (m.post_sleep_min >= 30) & (m.pre_sleep_min >= 30)].copy()
m = m[m.pre_sleep_t90 > 10].copy()
m["nonresp"] = (m.post_sleep_t90 > 10).astype(int)
SHIP0 = json.load(open(f"{paths.NUMBERS_DIR}/nonresponder_combined_v2.json"))   # V14_L5b: the sample size against the v8.1 file, no literal
assert len(m) == SHIP0["n"] and int(m.nonresp.sum()) == SHIP0["n_failed"], (len(m), int(m.nonresp.sum()), SHIP0["n"], SHIP0["n_failed"])
for k in DISEASES:
    c = f"{k}_prevalent"
    if c in m.columns:
        m[c] = m[c].fillna(0).astype(int)
m["male"] = (m.SexDSC.astype(str).str.upper().str[0] == "M").astype(int)
m["lbase"] = np.log1p(m.pre_sleep_t90)
m["cp_raw"] = sum(m[f"{k}_prevalent"] for k in CP)
m["cp"] = m.cp_raw.clip(upper=3)

PERSON = glob.glob(f"{paths.OMOP_CACHE_DIR}/**/person_merged.parquet",
                   recursive=True)[0]
p = duckdb.sql(f'select person_id as "BDSPPatientID", race_source_value from \'{PERSON}\'').df()
p = p.drop_duplicates("BDSPPatientID")
m = m.merge(p, on="BDSPPatientID", how="left")
raw = m.race_source_value.astype(str).str.strip().str.upper()
BLACK = {"AFRICAN AMERICAN  OR BLACK", "AFRICAN AMERICAN OR BLACK", "BLACK/AFRICAN AMERICAN",
         "BLACK", "AFRICAN AMERICAN"}
WHITE = {"WHITE", "CAUCASIAN OR WHITE", "CAUCASIAN"}
UNK = {"UNKNOWN, UNAVAILABLE OR UNREPORTED", "UNKNOWN/NOT SPECIFIED", "UNKNOWN",
       "UNABLE TO OBTAIN", "DECLINED TO ANSWER", "NONE", "PREFER NOT TO SAY", "NAN", "",
       "DECLINED", "NOT RECORDED", "PATIENT DECLINED", "NOT SPECIFIED", "UNAVAILABLE",
       "DECLINE TO ANSWER"}
m["race_black"] = raw.isin(BLACK).astype(int)
m["race_asian"] = (raw == "ASIAN").astype(int)
m["race_other"] = (~raw.isin(BLACK | WHITE | UNK) & (raw != "ASIAN")).astype(int)
m["race_unknown"] = raw.isin(UNK).astype(int)
RACE = ["race_black", "race_asian", "race_other", "race_unknown"]
DEMO = ["AgeAtVisit", "male"] + RACE


def logit(terms):
    X = sm.add_constant(m[terms].astype(float))
    f = pd.concat([X, m.nonresp], axis=1).dropna()
    r = sm.Logit(f.nonresp, f[X.columns]).fit(disp=0)
    ci = r.conf_int()
    return {t: {"or": float(np.exp(r.params[t])), "lo": float(np.exp(ci.loc[t, 0])),
                "hi": float(np.exp(ci.loc[t, 1])), "p": float(r.pvalues[t])} for t in terms}


def bh(pvals):
    n = len(pvals)
    order = np.argsort(pvals)
    q = np.empty(n, float)
    prev = 1.0
    for rank, i in enumerate(reversed(order), start=1):
        val = pvals[i] * n / (n - rank + 1)
        prev = min(prev, val)
        q[i] = min(prev, 1.0)
    return q


MIN_CELL = 20
label_of = {k: DISEASES[k][0] for k in DISEASES}
rows = []
for k in DISEASES:
    c = f"{k}_prevalent"
    if c not in m.columns:
        continue
    nw, nwo = int((m[c] == 1).sum()), int((m[c] == 0).sum())
    if nw < MIN_CELL or nwo < MIN_CELL:
        continue
    yw = m.loc[m[c] == 1, "nonresp"]
    if yw.sum() == 0 or yw.sum() == len(yw):
        continue
    v = logit([c, "AgeAtVisit", "male", "lbase"])[c]
    rows.append({"key": k, "Condition": label_of[k], "neg": k in NEGATIVE_CONTROLS,
                 "marginal_or": v["or"], "marginal_lo": v["lo"], "marginal_hi": v["hi"],
                 "marginal_p": v["p"]})
marg = pd.DataFrame(rows)
qs = []
for flag in (False, True):
    sub = marg[marg.neg == flag].copy()
    sub["marginal_q"] = bh(sub.marginal_p.values)
    qs.append(sub[["key", "marginal_q"]])
marg = marg.merge(pd.concat(qs, ignore_index=True), on="key", how="left")
SIG = [r.key for r in marg.itertuples() if not r.neg and r.marginal_q < 0.05]
joint = logit([f"{k}_prevalent" for k in SIG] + DEMO + ["lbase"])
indep = [k for k in SIG if joint[f"{k}_prevalent"]["p"] < 0.05]
trend = logit(["cp", "AgeAtVisit", "male", "lbase"])["cp"]

FP = {"by_condition": {r.Condition: {"marginal_or": r.marginal_or, "marginal_lo": r.marginal_lo,
                                     "marginal_hi": r.marginal_hi, "marginal_p": r.marginal_p,
                                     "marginal_q": r.marginal_q, "key": r.key}
                       for r in marg.itertuples()},
      "joint_independent": {label_of[k]: joint[f"{k}_prevalent"] for k in indep},
      "joint_demographics": {"AgeAtVisit": joint["AgeAtVisit"], "male": joint["male"],
                             "race_black": joint["race_black"]},
      "dose_response_trend": trend,
      "n": int(len(m)), "n_failed": int(m.nonresp.sum()),
      "n_conditions_tested": int((~marg.neg).sum()),
      "n_conditions_fdr_significant": len(SIG),
      "n_independent": len(indep)}

# ------------------------------------------------------------------ reconcile with the shipped file
SHIP = json.load(open(f"{paths.NUMBERS_DIR}/nonresponder_combined_v2.json"))
bad = []
missing = [lab for lab in FP["by_condition"] if lab not in SHIP["failure_rate_by_condition"]]
print("conditions refitted here but absent from the shipped json (negative controls are kept in a"
      f" separate block there): {missing}")
for lab, v in FP["by_condition"].items():
    if lab not in SHIP["failure_rate_by_condition"]:
        continue
    s = SHIP["failure_rate_by_condition"][lab]
    for a, b_ in (("marginal_or", "marginal_or"), ("marginal_lo", "marginal_lo"),
                  ("marginal_hi", "marginal_hi")):
        if abs(float(v[a]) - float(s[b_])) > 5e-4:
            bad.append((lab, a, v[a], s[b_]))
for lab, v in FP["by_condition"].items():   # V14_L5b: the negative controls live in their own block of the shipped file
    if lab in SHIP["negative_controls"]:
        s = SHIP["negative_controls"][lab]
        for a in ("marginal_or", "marginal_lo", "marginal_hi"):
            if abs(float(v[a]) - float(s[a])) > 5e-4: bad.append((lab, a, v[a], s[a]))
for lab, v in FP["joint_independent"].items():
    s = SHIP["independent_conditions"][lab]
    for a in ("or", "lo", "hi"):
        if abs(float(v[a]) - float(s[a])) > 5e-4: bad.append(("joint " + lab, a, v[a], s[a]))
for k, v in FP["joint_demographics"].items():
    s = SHIP["joint_model"]["demographics"][k]
    for a in ("or", "lo", "hi"):
        if abs(float(v[a]) - float(s[a])) > 5e-4: bad.append(("joint " + k, a, v[a], s[a]))
for a in ("or", "lo", "hi"):
    if abs(float(FP["dose_response_trend"][a]) - float(SHIP["dose_response_trend"][a])) > 5e-4: bad.append(("trend", a, FP["dose_response_trend"][a], SHIP["dose_response_trend"][a]))
assert not bad, bad
assert set(FP["joint_independent"]) == set(SHIP["independent_conditions"]), (set(FP["joint_independent"]) ^ set(SHIP["independent_conditions"]))
assert FP["n_conditions_fdr_significant"] == SHIP["n_conditions_fdr_significant"]
assert FP["n_independent"] == SHIP["n_conditions_independent_in_joint_model"]
print(f"reconciled with numbers/nonresponder_combined_v2.json: {len(FP['by_condition'])} conditions, "
      f"{FP['n_conditions_fdr_significant']} FDR survivors, {FP['n_independent']} independent, "
      "every value within 5e-4 of the shipped three decimals")

# ------------------------------------------------------------------ which printed strings move
print("\nSupp Fig 13 panel a, the FDR survivors, half up on the stored 3dp against half up on full precision")
moved13a = []
for lab, s in SHIP["failure_rate_by_condition"].items():
    if not s["fdr_significant"]:
        continue
    v = FP["by_condition"][lab]
    old = f"{half_up(s['marginal_or'], 2)} ({half_up(s['marginal_lo'], 2)}-{half_up(s['marginal_hi'], 2)})"
    new = f"{half_up(v['marginal_or'], 2)} ({half_up(v['marginal_lo'], 2)}-{half_up(v['marginal_hi'], 2)})"
    if old != new:
        moved13a.append((lab, old, new))
        print(f"  MOVES  {lab:28s} {old:22s} -> {new}")
print(f"  {len(moved13a)} of {sum(1 for s in SHIP['failure_rate_by_condition'].values() if s['fdr_significant'])} move")

print("\nSupp Fig 13 panel c, the joint model rows")
moved13c = []
for lab, s in SHIP["independent_conditions"].items():
    v = FP["joint_independent"][lab]
    old = f"{s['or']:.2f} ({s['lo']:.2f}-{s['hi']:.2f})"
    new = f"{half_up(v['or'], 2)} ({half_up(v['lo'], 2)}-{half_up(v['hi'], 2)})"
    if old != new:
        moved13c.append((lab, old, new))
        print(f"  MOVES  {lab:28s} {old:22s} -> {new}")
for nm, key in (("Age, per year", "AgeAtVisit"), ("Male sex", "male"), ("Black race", "race_black")):
    s = SHIP["joint_model"]["demographics"][key]
    v = FP["joint_demographics"][key]
    old = f"{s['or']:.2f} ({s['lo']:.2f}-{s['hi']:.2f})"
    new = f"{half_up(v['or'], 2)} ({half_up(v['lo'], 2)}-{half_up(v['hi'], 2)})"
    if old != new:
        moved13c.append((nm, old, new))
        print(f"  MOVES  {nm:28s} {old:22s} -> {new}")
s = SHIP["dose_response_trend"]
old = f"{s['or']:.2f} (95% CI, {s['lo']:.2f}-{s['hi']:.2f})"
new = (f"{half_up(trend['or'], 2)} (95% CI, {half_up(trend['lo'], 2)}-{half_up(trend['hi'], 2)})")
if old != new:
    moved13c.append(("panel b trend", old, new))
    print(f"  MOVES  {'panel b trend':28s} {old:22s} -> {new}")
print(f"  {len(moved13c)} move on panel b and c")

print("\nSupp Fig 14, the ten drawn rows, half up on the three decimals today against half up on full precision")
PH = pd.read_csv(f"{paths.NUMBERS_DIR}/nonresponder_phenotype_v2.csv")
surv = PH[(~PH.neg) & (PH.marginal_q < 0.05)].copy()
TEN = surv.sort_values(["marginal_or", "marginal_lo"], ascending=[False, False]).head(10)
moved14 = []
for r in TEN.itertuples():
    v = FP["by_condition"][r.Condition]
    old = (f"{half_up(r.marginal_or, 2)} ({half_up(r.marginal_lo, 2)}-{half_up(r.marginal_hi, 2)})")
    new = (f"{half_up(v['marginal_or'], 2)} ({half_up(v['marginal_lo'], 2)}-{half_up(v['marginal_hi'], 2)})")
    if old != new:
        moved14.append((r.Condition, old, new))
        print(f"  MOVES  {r.Condition:28s} {old:22s} -> {new}")
print(f"  {len(moved14)} of 10 move")

FP["moves"] = {"Supp_Fig13_panel_a": moved13a, "Supp_Fig13_panel_bc": moved13c,
               "Supp_Fig14": moved14}
os.makedirs(f"{LANE}/work", exist_ok=True)
OUTP = f"{LANE}/work/fullprec_nonresponder_v8_1.json"
json.dump(FP, open(OUTP, "w"), indent=1)
prov = {"output": OUTP, "output_sha256": sha256(OUTP), "written_local": datetime.now().isoformat(timespec="seconds"),
        "script": os.path.abspath(__file__), "script_sha256": sha256(os.path.abspath(__file__)),
        "inputs": [{"path": p_, "sha256": sha256(p_)} for p_ in IN_PARQ + [PERSON]],
        "reconciled_against": {"path": f"{paths.NUMBERS_DIR}/nonresponder_combined_v2.json",
                               "sha256": sha256(f"{paths.NUMBERS_DIR}/nonresponder_combined_v2.json"),
                               "tolerance": 5e-4},
        "model": "sm.Logit, spec lifted from fits_phenotype_combined_v2.py (v8.1 step 127), round(x, 3) left out", "step": {"id": "L5b-FP", "name": "fullprec_nonresponder_v8_1 (lane-side full-precision twin of step 127)"}, "output": {"mtime_local": datetime.now().strftime("%Y-%m-%d %H:%M:%S")},
        "python": sys.version.split()[0],
        "packages": subprocess.run([sys.executable, "-c",
                                    "import statsmodels,pandas,numpy;print(statsmodels.__version__,pandas.__version__,numpy.__version__)"],
                                   capture_output=True, text=True).stdout.strip()}
json.dump(prov, open(OUTP + ".provenance.json", "w"), indent=1)
print(f"\nwrote {OUTP}\nwrote {OUTP}.provenance.json")
