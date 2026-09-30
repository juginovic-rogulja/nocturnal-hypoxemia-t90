"""
Colleague request #4. Sleep duration as the exposure in the eTable 8 external-validation
models (SHHS and MrOS), alone and mutually adjusted with T90, for all 12 outcomes.

Model specification, cohort construction, covariates, follow-up, censoring and outcome
definitions are copied from numbers/freeze_05_external.py (the generator of the frozen
numbers/external.json behind eTable 8), with its three 2026-08-08 corrections already in
place: one age-spline column dropped, unpenalized maximum partial likelihood with a Newton
step-size scan, MrOS follow-up on FUVSDT. Nothing frozen is modified. This script only READS
the raw cohort files and WRITES into this folder.

POSITIVE CONTROL. Before any sleep-duration model runs, the harness refits the exact
eTable 8 estimates (T90 alone, AHI alone, mutually adjusted, 12 outcomes x 4 estimates)
and compares them to the frozen numbers/external.json within 5e-4 on HR, lo, hi (v8.3,
both sides full precision; was equality at 3 decimals) and exactly on event counts. The run aborts if any estimate fails to reproduce.

EXPOSURES (sleep duration = PSG total sleep time)
  SHHS  slpprdp   "Total Sleep Duration ... conventionally called total sleep time",
                  minutes, type II polysomnography (NSRR dictionary 0.21.0)
  MrOS  POSLPRDP  "TIME LIGHT OFF TO LIGHT ON SCORED AS SLEEP (MINS)" (POSFEB23 contents)
  Continuous: hours (minutes / 60), hazard ratio per 1-hour INCREASE in total sleep time.
  Categorical: <6 h, 6 to <7 h, 7 to <8 h (reference), >=8 h, half-open bins on minutes
  ( <360, 360-419.9, 420-479.9, >=480 ).

MODELS per outcome-cohort pair (12 pairs, same covariates as eTable 8: natural cubic
spline on age with 4 df and the first basis column dropped, plus sex in SHHS, site
stratification in MrOS)
  M1  sleep_h + covariates                       (continuous, alone)
  M2  sleep_h + t90_z + covariates               (continuous, mutually adjusted with T90)
  M3  lt6 + s6to7 + ge8 + covariates             (categorical, alone)
  M4  lt6 + s6to7 + ge8 + t90_z + covariates     (categorical, mutually adjusted with T90)
T90 enters exactly as in eTable 8: rank-based inverse-normal pctsa90h (SHHS) or
within-site rank-based inverse-normal POPCSA90 (MrOS).

MULTIPLICITY. Hand-rolled Benjamini-Hochberg within each exposure family across the 12
outcome-cohort pairs. Families: continuous alone, continuous T90-adjusted, each categorical
level alone, each categorical level T90-adjusted, T90 inside M2, T90 inside M4.

SANITY. Crude event rates per 1000 person-years by sleep category for SHHS cardiovascular
composite and MrOS death from any cause.

Writes: results.csv, summary.json, eTable18_draft.md, REPORT excerpts to stdout.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import warnings
warnings.filterwarnings("ignore")
import json
import os
import subprocess
import sys
import numpy as np
import pandas as pd
from scipy import stats
from patsy import dmatrix
from lifelines import CoxPHFitter

HERE = os.path.dirname(os.path.abspath(__file__))
NUM = paths.NUMBERS_DIR
SH = paths.SHHS_DIR
STEPS = [None, 0.05, 0.1, 0.25, 0.5, 0.95]
EVMIN = {"shhs": 30, "mros": 25}
LABEL = {"shhs": "Sleep Heart Health Study", "mros": "Osteoporotic Fractures in Men Study"}
CATS = ["lt6", "s6to7", "ge8"]
CATLAB = {"lt6": "Under 6 h", "s6to7": "6 to under 7 h", "ge8": "8 h or more",
          "ref": "7 to under 8 h"}


def mros_path(name):
    """Resolve an MrOS SAS file. Order: this folder's data/, the canonical extracted
    folders under $T90_MROS_DIR/, the tmp job folder the freeze script used (all three
    verified md5-identical 2026-08-19), else auto-extract from the canonical zip."""
    cands = [os.path.join(HERE, "data", name),
             f"{paths.MROS_DIR}/{name.split('.')[0].upper()}/{name}",
             os.path.expanduser(f"{paths.MROS_DIR}/{name}")]
    for c in cands:
        if os.path.exists(c):
            return c
    z = f"{paths.MROS_DIR}/{name.split('.')[0].upper()}.zip"
    if os.path.exists(z):
        os.makedirs(os.path.join(HERE, "data"), exist_ok=True)
        subprocess.run(["unzip", "-o", "-q", z, name, "-d", os.path.join(HERE, "data")],
                       check=True)
        return os.path.join(HERE, "data", name)
    raise FileNotFoundError(name)


# ---- fitting machinery, copied from freeze_05_external.py -------------------------------
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
    """Natural cubic spline on age, one basis column dropped (the 4 columns sum to one and
    the unpenalized fit is otherwise rank-deficient). Named agesp*, never age_s*."""
    b = dmatrix(f"cr(a, df={n}) - 1", {"a": df[col].values}, return_type="dataframe")
    b.columns = [f"agesp{i}" for i in range(b.shape[1])]
    b.index = df.index
    b = b.iloc[:, 1:]
    for c in b.columns:
        df[c] = b[c]
    return list(b.columns)


def est(c, x, f):
    lo, hi = np.exp(c.confidence_intervals_.loc[x])
    return {"n": int(len(f)), "events": int(f.E.sum()),
            "hr": float(np.exp(c.params_[x])), "lo": float(lo), "hi": float(hi),
            "p": float(c.summary.loc[x, "p"])}


def bh(pvals):
    """Hand-rolled Benjamini-Hochberg q values."""
    p = np.asarray(pvals, dtype=float)
    m = len(p)
    order = np.argsort(p)
    q = np.empty(m)
    prev = 1.0
    for rank_from_top, idx in enumerate(order[::-1]):
        k = m - rank_from_top
        prev = min(prev, p[idx] * m / k)
        q[idx] = prev
    return q


# ============================================================ load SHHS
d = pd.read_csv(f"{SH}/shhs1-dataset-0.21.0.csv", low_memory=False)
cv = pd.read_csv(f"{SH}/shhs-cvd-summary-dataset-0.21.0.csv", low_memory=False)
cv = cv.drop(columns=[c for c in ["gender", "race", "age_s1"] if c in cv.columns])
s = d[["nsrrid", "pctsa90h", "ahi_a0h3a", "ahi_a0h4a", "age_s1", "gender",
       "bmi_s1", "smokstat_s1", "diasbp", "chol", "hdl", "htnderv_s1",
       "parrptdiab", "slpprdp"]].merge(cv, on="nsrrid", how="inner")
s["male"] = (s.gender == 1).astype(int)
ADJ = spline(s, "age_s1") + ["male"]
s["t90_z"] = rint(s.pctsa90h)
s["ahi_z"] = rint(s.ahi_a0h3a)
s["sleep_h"] = s.slpprdp / 60.0
s["lt6"] = (s.slpprdp < 360).astype(int)
s["s6to7"] = ((s.slpprdp >= 360) & (s.slpprdp < 420)).astype(int)
s["ge8"] = (s.slpprdp >= 480).astype(int)

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


# ============================================================ load MrOS
pos = pd.read_sas(mros_path("POSFEB23.SAS7BDAT"), encoding="latin-1")
ef = pd.read_sas(mros_path("effeb24.sas7bdat"), encoding="latin-1")
vs = pd.read_sas(mros_path("vsfeb24.sas7bdat"), encoding="latin-1")
m = (pos[["ID", "POPCSA90", "PORDI4P", "POSLPRDP"]]
     .merge(ef[["ID", "SITE", "DADEAD", "DACARDIO", "DACANCER", "DASTROKE", "DAPULMON",
                "FUCDTIME", "FUVSDT"]], on="ID")
     .merge(vs[["ID", "VSAGE1", "HWBMI"]], on="ID", how="left"))
for c in m.columns:
    if c not in ("ID", "SITE"):
        m[c] = pd.to_numeric(m[c], errors="coerce")
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
m["sleep_h"] = m.POSLPRDP / 60.0
m["lt6"] = (m.POSLPRDP < 360).astype(int)
m["s6to7"] = ((m.POSLPRDP >= 360) & (m.POSLPRDP < 420)).astype(int)
m["ge8"] = (m.POSLPRDP >= 480).astype(int)

MOUT = [("DADEAD", "Death from any cause"), ("DACARDIO", "Cardiovascular death"),
        ("DAPULMON", "Pulmonary death"), ("DASTROKE", "Stroke death"),
        ("DACANCER", "Cancer death")]


def mframe(k, cols):
    f = pd.DataFrame({"T": (m.FUVSDT / 365.25).values,
                      "E": m[k].fillna(0).astype(int).values,
                      "site": m.SITE.values,
                      **{c: m[c].values for c in cols + MADJ}})
    return f.dropna()


# v8.3 (2026-09-16): the cohort sizes are read from the regenerated numbers/external.json (step 260, sidecar), not typed
with open(f"{NUM}/external.json") as fh:
    FROZEN = json.load(fh)
assert len(s) == FROZEN["shhs"]["n"], f"SHHS n drifted: {len(s)} vs frozen {FROZEN['shhs']['n']}"
assert len(m) == FROZEN["mros"]["n"], f"MrOS n drifted: {len(m)} vs frozen {FROZEN['mros']['n']}"

# ============================================================ positive control
# FROZEN (numbers/external.json) is loaded above, before the cohort-size asserts (v8.3)

control_rows = []
control_fail = 0


def check(cohort, lab, key, e):
    global control_fail
    fr = FROZEN[cohort]["outcomes"][lab].get(key)
    if fr is None:
        return
    # v8.3: reference now full precision, compare within 5e-4 (2026-09-16)
    ok = (abs(e["hr"] - fr["hr"]) <= 5e-4 and abs(e["lo"] - fr["lo"]) <= 5e-4
          and abs(e["hi"] - fr["hi"]) <= 5e-4 and e["events"] == fr["events"])
    if not ok:
        control_fail += 1
    control_rows.append({"cohort": cohort, "outcome": lab, "estimate": key,
                         "frozen_hr": fr["hr"], "refit_hr": e["hr"],  # v8.3 full precision (2026-09-16)
                         "frozen_lo": fr["lo"], "refit_lo": e["lo"],
                         "frozen_hi": fr["hi"], "refit_hi": e["hi"],
                         "frozen_events": fr["events"], "refit_events": e["events"],
                         "match": bool(ok)})


print("== POSITIVE CONTROL: refitting the eTable 8 estimates ==")
sfr = {}
for k, dc, pc, lab in SOUT:
    g = sframe(k, dc, pc)
    sfr[lab] = g
    for x in ["t90_z", "ahi_z"]:
        f = g[["T", "E", x] + ADJ].dropna()
        if f.E.sum() < EVMIN["shhs"]:
            continue
        c = fit_ml(f)
        check("shhs", lab, x.replace("_z", ""), est(c, x, f))
    f = g[["T", "E", "t90_z", "ahi_z"] + ADJ].dropna()
    if f.E.sum() >= EVMIN["shhs"]:
        c = fit_ml(f)
        check("shhs", lab, "t90_adj_ahi", est(c, "t90_z", f))
        check("shhs", lab, "ahi_adj_t90", est(c, "ahi_z", f))

for k, lab in MOUT:
    for x in ["t90_z", "ahi_z"]:
        f = mframe(k, [x])
        if f.E.sum() < EVMIN["mros"]:
            continue
        c = fit_ml(f, ["site"])
        check("mros", lab, x.replace("_z", ""), est(c, x, f))
    f = mframe(k, ["t90_z", "ahi_z"])
    if f.E.sum() >= EVMIN["mros"]:
        c = fit_ml(f, ["site"])
        check("mros", lab, "t90_adj_ahi", est(c, "t90_z", f))
        check("mros", lab, "ahi_adj_t90", est(c, "ahi_z", f))

cdf = pd.DataFrame(control_rows)
n_ok = int(cdf.match.sum())
print(f"  {n_ok}/{len(cdf)} frozen estimates reproduced within 5e-4 on HR, lo, hi "
      f"(full precision on both sides, v8.3), exact events")
if control_fail:
    print(cdf[~cdf.match].to_string(index=False))
    sys.exit("POSITIVE CONTROL FAILED, aborting before any sleep-duration model")

# ============================================================ sleep-duration models
rows = []


def add(cohort, lab, model, term, e, level_n=None, level_ev=None):
    rows.append({"block": "sleep_duration", "cohort": cohort, "cohort_label": LABEL[cohort],
                 "outcome": lab, "model": model, "term": term, "n": e["n"],
                 "events": e["events"], "level_n": level_n, "level_events": level_ev,
                 "hr": e["hr"], "lo": e["lo"], "hi": e["hi"], "p": e["p"]})


def run_pair(cohort, lab, frame_fn, strata):
    """frame_fn(cols) -> fitted frame with T, E, requested exposure columns, covariates."""
    evmin = EVMIN[cohort]
    # M1 continuous alone
    f = frame_fn(["sleep_h"])
    if f.E.sum() >= evmin:
        c = fit_ml(f, strata)
        add(cohort, lab, "M1_cont_alone", "sleep_h", est(c, "sleep_h", f))
    # M2 continuous + T90
    f = frame_fn(["sleep_h", "t90_z"])
    if f.E.sum() >= evmin:
        c = fit_ml(f, strata)
        add(cohort, lab, "M2_cont_adjT90", "sleep_h", est(c, "sleep_h", f))
        add(cohort, lab, "M2_cont_adjT90", "t90_z", est(c, "t90_z", f))
    # M3 categorical alone
    f = frame_fn(CATS)
    if f.E.sum() >= evmin:
        c = fit_ml(f, strata)
        for t in CATS:
            ln = int(f[t].sum())
            lev = int(f.loc[f[t] == 1, "E"].sum())
            add(cohort, lab, "M3_cat_alone", t, est(c, t, f), ln, lev)
    # M4 categorical + T90
    f = frame_fn(CATS + ["t90_z"])
    if f.E.sum() >= evmin:
        c = fit_ml(f, strata)
        for t in CATS:
            ln = int(f[t].sum())
            lev = int(f.loc[f[t] == 1, "E"].sum())
            add(cohort, lab, "M4_cat_adjT90", t, est(c, t, f), ln, lev)
        add(cohort, lab, "M4_cat_adjT90", "t90_z", est(c, "t90_z", f))


for k, dc, pc, lab in SOUT:
    g = sfr[lab]
    run_pair("shhs", lab, lambda cols, g=g: g[["T", "E"] + cols + ADJ].dropna(), None)
for k, lab in MOUT:
    run_pair("mros", lab, lambda cols, k=k: mframe(k, cols), ["site"])

res = pd.DataFrame(rows)

# family assignment for BH: one family per model x term, 12 pairs each
FAMY = {("M1_cont_alone", "sleep_h"): "dur_cont_alone",
        ("M2_cont_adjT90", "sleep_h"): "dur_cont_adjT90",
        ("M2_cont_adjT90", "t90_z"): "t90_in_cont_model",
        ("M3_cat_alone", "lt6"): "cat_lt6_alone",
        ("M3_cat_alone", "s6to7"): "cat_6to7_alone",
        ("M3_cat_alone", "ge8"): "cat_ge8_alone",
        ("M4_cat_adjT90", "lt6"): "cat_lt6_adjT90",
        ("M4_cat_adjT90", "s6to7"): "cat_6to7_adjT90",
        ("M4_cat_adjT90", "ge8"): "cat_ge8_adjT90",
        ("M4_cat_adjT90", "t90_z"): "t90_in_cat_model"}
res["family"] = [FAMY[(r.model, r.term)] for r in res.itertuples()]
res["q"] = np.nan
for fam, idx in res.groupby("family").groups.items():
    res.loc[idx, "q"] = bh(res.loc[idx, "p"].values)
res["unstable"] = [(r.level_events is not None and not pd.isna(r.level_events)
                    and r.level_events < 5) for r in res.itertuples()]

# ============================================================ crude-rate sanity checks
def crude(frame, tag):
    lines = []
    grp = np.select([frame.lt6 == 1, frame.s6to7 == 1, frame.ge8 == 1],
                    ["<6 h", "6 to <7 h", ">=8 h"], default="7 to <8 h")
    for gname in ["<6 h", "6 to <7 h", "7 to <8 h", ">=8 h"]:
        sub = frame[grp == gname]
        py = float(sub["T"].sum())
        ev = int(sub.E.sum())
        lines.append({"sanity": tag, "group": gname, "n": int(len(sub)), "events": ev,
                      "person_years": round(py, 1),
                      "rate_per_1000py": round(1000 * ev / py, 2) if py > 0 else None})
    return lines


san = []
g = sfr["Cardiovascular composite"][["T", "E", "lt6", "s6to7", "ge8"] + ADJ].dropna()
san += crude(g, "SHHS cardiovascular composite")
f = mframe("DADEAD", CATS)
san += crude(f, "MrOS death from any cause")

# ============================================================ write outputs
full = pd.concat([cdf.assign(block="positive_control"), res], ignore_index=True)
full.to_csv(os.path.join(HERE, "results.csv"), index=False)

summary = {"positive_control": {"n_estimates": int(len(cdf)), "n_match": n_ok,
                                "all_match": bool(control_fail == 0)},
           "spec": {"source": "numbers/freeze_05_external.py model spec, exposures swapped",
                    "fit": "Cox by maximum partial likelihood, penalizer 0, Newton step "
                           "sizes scanned, highest-likelihood solution kept",
                    "age_basis": "natural cubic spline, 4 df, first basis column dropped",
                    "shhs_covariates": "age spline + male",
                    "mros_covariates": "age spline, stratified by site",
                    "sleep_duration_vars": {"shhs": "slpprdp", "mros": "POSLPRDP"},
                    "t90_vars": {"shhs": "pctsa90h rank-inverse-normal",
                                 "mros": "POPCSA90 rank-inverse-normal within site"},
                    "categories_minutes": {"lt6": "<360", "6to7": "360-419.9",
                                           "ref_7to8": "420-479.9", "ge8": ">=480"},
                    "bh": "within each model x term family across the 12 outcome-cohort "
                          "pairs"},
           "cohort_n": {"shhs": int(len(s)), "mros": int(len(m))},
           "category_counts": {
               "shhs": {"lt6": int(s.lt6.sum()), "6to7": int(s.s6to7.sum()),
                        "7to8": int(((s.slpprdp >= 420) & (s.slpprdp < 480)).sum()),
                        "ge8": int(s.ge8.sum())},
               "mros": {"lt6": int(m.lt6.sum()), "6to7": int(m.s6to7.sum()),
                        "7to8": int(((m.POSLPRDP >= 420) & (m.POSLPRDP < 480)).sum()),
                        "ge8": int(m.ge8.sum())}},
           "sanity_crude_rates": san,
           "fdr_survivors": res[res.q < 0.05][
               ["cohort", "outcome", "family", "term", "hr", "p", "q"]
           ].to_dict("records")}
with open(os.path.join(HERE, "summary.json"), "w") as fh:
    json.dump(summary, fh, indent=2)

print("\n== SLEEP-DURATION RESULTS (sleep terms only) ==")
for r in res[res.term != "t90_z"].itertuples():
    star = "*" if r.unstable else " "
    print(f"  {r.cohort:4s} {r.outcome:26s} {r.model:14s} {r.term:6s} "
          f"HR {r.hr:6.3f} ({r.lo:6.3f}-{r.hi:6.3f}) p={r.p:.4f} q={r.q:.3f}{star}")
print("\n== T90 INSIDE THE MUTUALLY ADJUSTED MODELS ==")
for r in res[res.term == "t90_z"].itertuples():
    print(f"  {r.cohort:4s} {r.outcome:26s} {r.model:14s} "
          f"HR {r.hr:6.3f} ({r.lo:6.3f}-{r.hi:6.3f}) p={r.p:.4f} q={r.q:.3f}")
print("\n== CRUDE-RATE SANITY ==")
for x in san:
    print(f"  {x['sanity']:32s} {x['group']:10s} n={x['n']:5d} ev={x['events']:5d} "
          f"rate/1000py={x['rate_per_1000py']}")
surv = summary["fdr_survivors"]
print(f"\nFDR survivors across all families, sleep terms and T90 terms: {len(surv)}")
for x in surv:
    print(" ", x)


# ============================================================ eTable 18 draft
def cistr(r, star=False):
    if r is None or r["level_events"] == 0 or not np.isfinite(r["hi"]):
        return "—"
    txt = f"{r['hr']:.3f} ({r['lo']:.3f}-{r['hi']:.3f})"
    if star and r.get("level_events") is not None and r["level_events"] < 5:
        txt += "*"
    return txt


def pstr(r):
    if r is None or r["level_events"] == 0 or not np.isfinite(r["hi"]):
        return "—"
    return "<.001" if r["p"] < 0.001 else f"{r['p']:.3f}"


def pick(cohort, outcome, model, term):
    sub = res[(res.cohort == cohort) & (res.outcome == outcome)
              & (res.model == model) & (res.term == term)]
    if sub.empty:
        return None
    r = sub.iloc[0].to_dict()
    if r["level_events"] is None or pd.isna(r["level_events"]):
        r["level_events"] = -1  # sentinel: continuous term, never "0 events"
    return r


def md_table(df):
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |",
             "|" + "|".join("---" for _ in cols) + "|"]
    for _, row in df.iterrows():
        lines.append("| " + " | ".join(str(row[c]) for c in cols) + " |")
    return "\n".join(lines)


PAIRS = [("shhs", lab) for _, _, _, lab in SOUT] + [("mros", lab) for _, lab in MOUT]

pa = []
for coh, lab in PAIRS:
    a = pick(coh, lab, "M1_cont_alone", "sleep_h")
    b = pick(coh, lab, "M2_cont_adjT90", "sleep_h")
    o = pick(coh, lab, "M2_cont_adjT90", "t90_z")
    pa.append({"Cohort": LABEL[coh], "Outcome": lab, "Events": a["events"],
               "Sleep duration alone (95% CI)": cistr(a), "P value": pstr(a),
               "q value": f"{a['q']:.3f}",
               "Sleep duration, oxygen-adjusted (95% CI)": cistr(b),
               "Adjusted P value": pstr(b), "Adjusted q value": f"{b['q']:.3f}",
               "Oxygen per 1 SD, same model (95% CI)": cistr(o)})

pb = []
for coh, lab in PAIRS:
    for t in CATS:
        a = pick(coh, lab, "M3_cat_alone", t)
        b = pick(coh, lab, "M4_cat_adjT90", t)
        pb.append({"Cohort": LABEL[coh], "Outcome": lab, "Sleep group": CATLAB[t],
                   "n": int(a["level_n"]), "Events": int(a["level_events"]),
                   "Alone (95% CI)": cistr(a, star=True), "P value": pstr(a),
                   "Oxygen-adjusted (95% CI)": cistr(b, star=True),
                   "Adjusted P value": pstr(b)})

pc = []
for coh, lab in PAIRS:
    fr = FROZEN[coh]["outcomes"][lab]["t90"]
    o2 = pick(coh, lab, "M2_cont_adjT90", "t90_z")
    o4 = pick(coh, lab, "M4_cat_adjT90", "t90_z")
    pc.append({"Cohort": LABEL[coh], "Outcome": lab,
               "Oxygen alone, frozen eTable 8 (95% CI)":
                   f"{fr['hr']:.3f} ({fr['lo']:.3f}-{fr['hi']:.3f})",
               "With continuous sleep duration (95% CI)": cistr(o2),
               "P value, continuous model": pstr(o2),
               "With sleep duration categories (95% CI)": cistr(o4),
               "P value, categorical model": pstr(o4)})

qcat_min = float(res[res.family.str.startswith("cat_")].q.min())
CAPTION = (
    "**eTable 18.** Sleep duration entered as the exposure in the external-validation Cox "
    "models of eTable 8, alone and mutually adjusted with nocturnal oxygen, in the Sleep "
    "Heart Health Study (n=5,802, 7 adjudicated cardiovascular outcomes) and the "
    "Osteoporotic Fractures in Men Study (n=2,911, 5 adjudicated causes of death). Sleep "
    "duration is the total sleep time of the baseline polysomnogram (SHHS variable slpprdp, "
    "MrOS variable POSLPRDP, minutes divided by 60). Panel A enters it continuously, with "
    "the hazard ratio per 1-hour increase. Panel B enters it in categories of under 6 hours, "
    "6 to under 7 hours, and 8 hours or more, each against a reference of 7 to under 8 "
    "hours. Panel C shows the nocturnal-oxygen term from the same mutually adjusted models "
    "beside its frozen eTable 8 value. All models carry the eTable 8 covariates: a natural "
    "cubic spline on age with 4 df and one basis column dropped, plus sex in SHHS, with "
    "stratification by clinical site in MrOS. Oxygen is the percentage of sleep time below "
    "90% saturation (SHHS pctsa90h, MrOS POPCSA90), rank-based inverse-normal transformed, "
    "within site for MrOS, exactly as in eTable 8. Fits are by unpenalized maximum partial "
    "likelihood. The q values are Benjamini-Hochberg corrected within each exposure family "
    "across the 12 cohort-outcome pairs. No categorical sleep term survived correction "
    f"(smallest q={qcat_min:.2f}), so Panel B omits the q column, which is available in "
    "results.csv. "
    "An asterisk marks an estimate resting on fewer than 5 events in the exposed cell. An "
    "em-dash marks the one cell with no events (atrial fibrillation, 8 hours or more, "
    "SHHS), where the hazard ratio is not estimable.\n")

parts = [CAPTION,
         "\n**Panel A. Sleep duration, continuous per 1-hour increase**\n",
         md_table(pd.DataFrame(pa)),
         "\n\n**Panel B. Sleep duration in categories, reference 7 to under 8 hours**\n",
         md_table(pd.DataFrame(pb)),
         "\n\n**Panel C. The nocturnal-oxygen term in the mutually adjusted models**\n",
         md_table(pd.DataFrame(pc)), ""]
with open(os.path.join(HERE, "eTable18_draft.md"), "w") as fh:
    fh.write("\n".join(parts))

# T90 stability line for the report
shift = []
for coh, lab in PAIRS:
    fr = FROZEN[coh]["outcomes"][lab]["t90"]["hr"]
    for mdl in ["M2_cont_adjT90", "M4_cat_adjT90"]:
        o = pick(coh, lab, mdl, "t90_z")
        if o is not None:
            shift.append(abs(o["hr"] - fr))
print(f"\nT90 stability: max |HR shift| vs frozen T90-alone across all 24 mutually "
      f"adjusted models = {max(shift):.4f}")

print(f"\nwritten -> {HERE}/results.csv, summary.json, eTable18_draft.md")
