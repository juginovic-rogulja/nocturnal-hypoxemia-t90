"""
Hostile independent recomputation of the habitual-sleep analysis.

Nothing is read from analysis.parquet, results_habitual.csv or results_crossed.csv. The cohort,
the habitual bands, the oxygen bands and every hazard ratio are rebuilt from the two upstream
sources, t90_final.parquet and habitual_per_patient.parquet, so an error in the build step shows
up as a disagreement rather than being inherited.

Six attacks, in the order they were commissioned:

  1 selection      does the T90 result itself survive inside the 8,711, and does merely HAVING a
                   habitual value predict disease
  2 measurement    does the habitual answer flip across tier A, tier B, tier C and the template
                   field
  3 reverse cause  when was the note written relative to the study and to the diagnosis, and what
                   happens under a one- and two-year landmark
  4 note count     does the number of duration-bearing notes predict outcomes on its own, and does
                   adjusting for it move the habitual estimates
  5 power          minimum detectable hazard ratio per cell, and which nulls carry no information
  6 head to head   habitual versus laboratory total sleep time, per-SD and in cross-validated
                   held-out concordance computed with a fresh resampling scheme

Every model is an unpenalized Cox stratified by hospital with a natural cubic age spline
(cr(age, df=4), FIRST COLUMN DROPPED so the design is full rank) plus sex. A 4-column basis is
rank-deficient and an unpenalized fit silently returns nothing.

Writes attack_*.csv and attack_summary.json alongside this file, and prints the full transcript.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import json
import os
import sys
import time
import warnings

# This Mac runs several of these jobs at once. Each linear-algebra call otherwise grabs eight
# BLAS threads and the eight cores end up oversubscribed by a factor of a hundred, which makes a
# 0.6-second Cox fit take twenty. Pin to one thread per process before numpy is imported.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from patsy import dmatrix
from scipy import stats

warnings.filterwarnings("ignore")
pd.set_option("display.width", 200)

T90 = paths.T90_ROOT
HAB = f"{paths.SV_ROOT}/habitual_sleep_v2"
OUT = f"{paths.SV_ROOT}/habitual_outcomes_v2"
sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
from cohort_spec import apply_cohort, COHORT_N          # noqa: E402
from disease_definitions import (DISEASES, NEGATIVE_CONTROLS, CIRCULAR,   # noqa: E402
                                 CONFOUNDING_FLOOR_PER_SD)

T0 = time.time()
LOG = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.append(s)


def hdr(t):
    say("")
    say("=" * 100)
    say(t)
    say("=" * 100)


# ======================================================================================
# BUILD, from upstream only
# ======================================================================================
hdr("BUILD  (rebuilt from t90_final.parquet + habitual_per_patient.parquet, nothing inherited)")

b = apply_cohort(pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet"))
say(f"cohort after apply_cohort: {len(b):,}  (expected {COHORT_N:,})")

h = pd.read_parquet(f"{HAB}/habitual_per_patient_v2.parquet")
say(f"habitual_per_patient: {len(h):,} rows, {h.BDSPPatientID.nunique():,} unique patients")
h = h[["BDSPPatientID", "tier", "hours", "n_records", "n_tier", "n_pre_study",
       "note_date_first", "note_date_last", "any_third_party", "hours_min", "hours_max",
       "from_avg_nightly_hours_field"]].copy()
b = b.merge(h, on="BDSPPatientID", how="left", validate="one_to_one")
b["has_habitual"] = b.hours.notna().astype(int)
say(f"merged with a habitual value: {int(b.has_habitual.sum()):,} "
    f"({100*b.has_habitual.mean():.1f}%)   without: {int((1-b.has_habitual).sum()):,}")
with open(f"{OUT}/summary.json") as _fh:   # step 149 (build.py) output with sidecar: the build's own counts, never retyped. v8 F1 fix (2026-09-12), RECONSTRUCTED by lane R1 2026-09-13 from F1 STATUS + F1 smoke_attack.log after an iCloud rollback (F1's applier is evicted; documented behaviour, not F1's byte-identical text)
    _S_BUILD = json.load(_fh)
N_HAB_BUILD = int(_S_BUILD["step1"]["n_with_habitual"])
BAND_N_BUILD = [int(_S_BUILD["step2"]["n"][k]) for k in ["<5h", "5-<6h", "6-<7h(ref)", ">=7h"]]
assert int(b.has_habitual.sum()) == N_HAB_BUILD, f"habitual merge disagrees with the build ({N_HAB_BUILD:,})"

# habitual bands, left-closed exactly as specified
b["hab_cat"] = pd.cut(b.hours, [0, 5, 6, 7, 99], right=False,
                      labels=["<5h", "5-<6h", "6-<7h(ref)", ">=7h"])
say("habitual bands (left-closed): " +
    str({str(k): int(v) for k, v in b.hab_cat.value_counts().items()}))
assert list(b.hab_cat.value_counts().reindex(["<5h", "5-<6h", "6-<7h(ref)", ">=7h"])) == \
    BAND_N_BUILD, f"band counts disagree with the build {BAND_N_BUILD} (v1: 2693/907/1153/3958)"

# oxygen bands
b["o2_cat"] = pd.cut(b.spo2_pct_below_90, [-1, 1, 10, 1e9],
                     labels=["normal T90<=1%", "intermediate 1-10%", "low T90>10%"])
say("oxygen bands: " + str({str(k): int(v) for k, v in b.o2_cat.value_counts().items()}))

# timing of the note relative to the study. The upstream builder takes the median over the
# PRE-study notes whenever a patient has any, so "pre-study" means n_pre_study > 0. The stricter
# reading, every note in the patient's tier written before the study, is carried alongside
# because it is the one that actually rules out a post-study note contributing.
b["pre_study_only"] = ((b.has_habitual == 1) & (b.n_pre_study > 0)).astype(int)
b["pre_study_strict"] = ((b.has_habitual == 1) & (b.n_pre_study == b.n_tier)).astype(int)
b["post_study_only"] = ((b.has_habitual == 1) & (b.n_pre_study == 0)).astype(int)
say(f"pre-study (value taken from pre-study notes): {int(b.pre_study_only.sum()):,}   "
    f"post-study only: {int(b.post_study_only.sum()):,}   "
    f"strictly pre-study, no post-study note anywhere in the tier: "
    f"{int(b.pre_study_strict.sum()):,}")

# race and body mass index, for the selection table only. Race is harmonised from the same
# person file freeze_02_demographics.py uses; BMI from the frozen derivation.
try:
    import duckdb
    PERSON = f"{paths.OMOP_CACHE_DIR}/person_merged.parquet"
    p = duckdb.sql(f'select person_id as "BDSPPatientID", race_source_value '
                   f"from '{PERSON}'").df().drop_duplicates("BDSPPatientID")
    raw = p.race_source_value.astype(str).str.strip().str.upper()
    p["race_h"] = np.where(raw.isin(["WHITE", "CAUCASIAN OR WHITE", "CAUCASIAN"]), "White",
                           np.where(raw.str.contains("BLACK|AFRICAN AMERICAN", regex=True),
                                    "Black or African American", "Other/unknown"))
    b = b.merge(p[["BDSPPatientID", "race_h"]], on="BDSPPatientID", how="left")
    b["race_h"] = b.race_h.fillna("Other/unknown")
except Exception as e:                                            # pragma: no cover
    say(f"race merge unavailable ({e}); race dropped from the selection table")
    b["race_h"] = "Other/unknown"
try:
    bmi = pd.read_parquet(f"{paths.NUMBERS_DIR}/_bmi_derived.parquet")[["BDSPPatientID", "bmi"]]
    b = b.merge(bmi.drop_duplicates("BDSPPatientID"), on="BDSPPatientID", how="left")
except Exception as e:                                            # pragma: no cover
    say(f"BMI merge unavailable ({e})")
    b["bmi"] = np.nan
assert len(b) == COHORT_N
say(f"body mass index available for {int(b.bmi.notna().sum()):,} of {len(b):,} "
    f"({100*b.bmi.notna().mean():.1f}%)")

# adjustment set: age spline computed ONCE on the whole cohort so every subset uses one basis
b["male"] = (b.sex.astype(str).str.upper().str[0] == "M").astype(int)
sp = dmatrix("cr(a, df=4) - 1", {"a": b.AgeAtVisit.values}, return_type="dataframe").iloc[:, 1:]
for i in range(sp.shape[1]):
    b[f"age_s{i}"] = sp.iloc[:, i].values
ADJ = [f"age_s{i}" for i in range(sp.shape[1])] + ["male"]
say(f"age spline columns kept: {sp.shape[1]} of 4 (one dropped, unpenalized fit)")

# prevalent burden, a healthcare-contact proxy that exists for every patient
KEYS = [c[:-9] for c in b.columns if c.endswith("_incident")]
b["prevalent_count"] = b[[f"{k}_prevalent" for k in KEYS if f"{k}_prevalent" in b.columns]].sum(1)

OUTCOMES = [k for k in KEYS if f"{k}_years" in b.columns]
NONCIRC = [k for k in OUTCOMES if k not in CIRCULAR]
LABEL = {k: (DISEASES[k][0] if k in DISEASES else "Death from any cause") for k in OUTCOMES}
say(f"outcomes carried: {len(OUTCOMES)}   non-circular: {len(NONCIRC)}   "
    f"circular excluded: {CIRCULAR}")
say(f"negative controls: {NEGATIVE_CONTROLS}  (standing per-SD floor {CONFOUNDING_FLOOR_PER_SD})")


def rint(x):
    r = stats.rankdata(x)
    return stats.norm.ppf((r - 0.375) / (len(r) + 0.25))


def add_exposures(d):
    """Per-SD exposure columns, recomputed inside whatever set is passed in."""
    d = d.copy()
    d["z_t90"] = d.groupby("site_id").spo2_pct_below_90.transform(rint)
    d["z_tst"] = d.groupby("site_id").TST_min.transform(lambda x: rint(x.fillna(x.median())))
    d["z_ahi"] = d.groupby("site_id").AHI.transform(lambda x: rint(x.fillna(x.median())))
    if d.hours.notna().any():
        d["z_hab"] = (d.hours - d.hours.mean()) / d.hours.std()          # plain SD, as published
        d["z_hab_rint"] = d.groupby("site_id").hours.transform(
            lambda x: rint(x) if x.notna().all() else np.nan)            # like-for-like with T90
    return d


MINEV = 20


def frame(d, key, extra=()):
    """Rows at risk for one outcome with the survival columns attached."""
    yc, ec, pc = f"{key}_years", f"{key}_incident", f"{key}_prevalent"
    f = d[(d[pc] == 0) & d[yc].notna() & (d[yc] > 0)].copy()
    f["T"] = f[yc]
    f["E"] = f[ec].astype(int)
    return f


def cox(f, cols, strata=True, robust=False):
    """One unpenalized stratified Cox. Returns {col: (hr, lo, hi, p, se)} or None."""
    keep = list(cols) + ADJ + ["T", "E"] + (["site_id"] if strata else [])
    d = f[keep].dropna()
    if d.E.sum() < MINEV or len(d) < 50:
        return None
    if strata and d.groupby("site_id").E.sum().min() < 3:
        strata = False
        d = d.drop(columns=["site_id"])
    try:
        c = CoxPHFitter().fit(d, "T", "E", strata=(["site_id"] if strata else None),
                              robust=robust)
    except Exception:
        return None
    out = {"n": int(len(d)), "events": int(d.E.sum())}
    for col in cols:
        if col not in c.summary.index:
            continue
        r = c.summary.loc[col]
        out[col] = dict(hr=float(r["exp(coef)"]), lo=float(r["exp(coef) lower 95%"]),
                        hi=float(r["exp(coef) upper 95%"]), p=float(r["p"]),
                        se=float(r["se(coef)"]))
    return out


def bh(p):
    p = np.asarray(p, float)
    ok = ~np.isnan(p)
    q = np.full(len(p), np.nan)
    o = np.argsort(p[ok])
    v = p[ok][o]
    n = len(v)
    adj = np.minimum.accumulate((v * n / np.arange(1, n + 1))[::-1])[::-1]
    tmp = np.empty(n)
    tmp[o] = adj
    q[ok] = tmp
    return q


B = add_exposures(b)
HABSET = B[B.has_habitual == 1].copy()
HABSET = add_exposures(HABSET)                       # exposures re-standardised inside the subset
PRESET = B[B.pre_study_only == 1].copy()
PRESET = add_exposures(PRESET)
say(f"analysis sets: full {len(B):,}   habitual {len(HABSET):,}   pre-study {len(PRESET):,}")
say(f"habitual hours SD inside the habitual set: {HABSET.hours.std():.4f} h")

RES = {}

# ======================================================================================
# 1  SELECTION
# ======================================================================================
hdr("ATTACK 1  SELECTION.  Does T90 itself reproduce inside the 8,711, and does merely HAVING "
    "a habitual value predict disease?")

# ---- 1a standardised mean differences, recomputed
rows = []
cont = ["AgeAtVisit", "TST_min", "sleep_efficiency_pct", "AHI", "spo2_pct_below_90",
        "spo2_mean", "odi4_total", "odi3_total", "arousal_index", "recording_dur_min",
        "bmi", "prevalent_count", "N3_pct", "REM_pct"]
g1, g0 = B[B.has_habitual == 1], B[B.has_habitual == 0]
for v in cont:
    a, c = g1[v].dropna(), g0[v].dropna()
    sd = np.sqrt((a.var() + c.var()) / 2)
    rows.append(dict(var=v, kind="continuous", with_hab=a.mean(), without=c.mean(),
                     smd=(a.mean() - c.mean()) / sd if sd else np.nan,
                     n_with=len(a), n_without=len(c)))
binaries = {"male": B.male, "site_I0002": (B.site_id == "I0002").astype(int),
            "white": (B.race_h == "White").astype(int),
            "black": (B.race_h == "Black or African American").astype(int),
            "died": B.death_incident.astype(float),
            "any_prevalent": (B.prevalent_count > 0).astype(int),
            "t90_over_10pct": (B.spo2_pct_below_90 > 10).astype(int),
            "t90_under_1pct": (B.spo2_pct_below_90 <= 1).astype(int)}
for v, s in binaries.items():
    p1, p0 = s[B.has_habitual == 1].mean(), s[B.has_habitual == 0].mean()
    sd = np.sqrt((p1 * (1 - p1) + p0 * (1 - p0)) / 2)
    rows.append(dict(var=v, kind="binary", with_hab=p1, without=p0,
                     smd=(p1 - p0) / sd if sd else np.nan, n_with=int(B.has_habitual.sum()),
                     n_without=int((1 - B.has_habitual).sum())))
smd = pd.DataFrame(rows)
smd["abs_smd"] = smd.smd.abs()
smd = smd.sort_values("abs_smd", ascending=False)
smd.to_csv(f"{OUT}/attack_selection_smd.csv", index=False)
say(f"variables with |SMD| > 0.10: {int((smd.abs_smd > 0.10).sum())} of {len(smd)}")
say(smd.head(12).to_string(index=False, float_format=lambda x: f"{x:.3f}"))

# ---- 1b does having a habitual value predict outcomes, in the full cohort
rows = []
for k in OUTCOMES:
    f = frame(B, k)
    r = cox(f, ["has_habitual"])
    if r and "has_habitual" in r:
        rows.append(dict(key=k, disease=LABEL[k], control=k in NEGATIVE_CONTROLS,
                         n=r["n"], events=r["events"], **r["has_habitual"]))
sel = pd.DataFrame(rows)
sel["q"] = bh(sel.p.values)
sel.to_csv(f"{OUT}/attack_selection_ishabitual.csv", index=False)
say("")
say(f"HAVING a habitual note, as an exposure in the full {len(B):,}, age+sex adjusted, "
    f"site-stratified:")
say(f"  significant for {int((sel.p < .05).sum())} of {len(sel)} outcomes, "
    f"{int((sel.q < .05).sum())} survive Benjamini-Hochberg")
say(f"  median HR {sel.hr.median():.3f}, range {sel.hr.min():.2f} to {sel.hr.max():.2f}, "
    f"{int((sel.hr > 1).sum())} of {len(sel)} above 1")
say("  largest: " + ", ".join(
    f"{LABEL[r.key]} {r.hr:.2f} [{r.lo:.2f}-{r.hi:.2f}]"
    for r in sel.reindex(sel.hr.sub(1).abs().sort_values(ascending=False).index).head(6)
    .itertuples()))

# ---- 1c the decisive check: T90 per SD, full cohort vs inside the habitual subgroup
def per_sd_panel(d, expo, tag):
    rows = []
    for k in OUTCOMES:
        f = frame(d, k)
        r = cox(f, [expo])
        if r and expo in r:
            rows.append(dict(set=tag, exposure=expo, key=k, disease=LABEL[k],
                             control=k in NEGATIVE_CONTROLS, circular=k in CIRCULAR,
                             n=r["n"], events=r["events"], **r[expo]))
    out = pd.DataFrame(rows)
    if len(out):
        out["q"] = bh(out.p.values)
        out["mdhr80"] = np.exp(2.802 * out.se)
    return out


t90_full = per_sd_panel(B, "z_t90", "full_19173")
t90_hab = per_sd_panel(HABSET, "z_t90", "habitual_8711")
t90_no = per_sd_panel(B[B.has_habitual == 0].pipe(add_exposures), "z_t90", "no_habitual_10462")
t90_pre = per_sd_panel(PRESET, "z_t90", "pre_study_4631")
t90cmp = pd.concat([t90_full, t90_hab, t90_no, t90_pre])
t90cmp.to_csv(f"{OUT}/attack_t90_reproduction.csv", index=False)

m = t90_full.merge(t90_hab, on="key", suffixes=("_full", "_hab"))
m = m[~m.circular_full]
lf, lh = np.log(m.hr_full), np.log(m.hr_hab)
slope = np.polyfit(lf, lh, 1)[0]
say("")
say(f"T90 per SD, {len(m)} non-circular outcomes fitted in BOTH sets:")
say(f"  Pearson r of log HR  {np.corrcoef(lf, lh)[0,1]:.3f}   Spearman "
    f"{stats.spearmanr(lf, lh).statistic:.3f}")
say(f"  regression of subgroup log HR on full-cohort log HR: slope {slope:.3f} "
    f"(1.0 = identical magnitude)")
say(f"  same direction: {int((np.sign(lf) == np.sign(lh)).sum())} of {len(m)}")
say(f"  median HR full {m.hr_full.median():.3f} vs subgroup {m.hr_hab.median():.3f}")
say(f"  significant: full {int((m.p_full < .05).sum())} of {len(m)}, "
    f"subgroup {int((m.p_hab < .05).sum())} of {len(m)}")
say(f"  clearing the {CONFOUNDING_FLOOR_PER_SD} floor: full "
    f"{int(((m.hr_full > CONFOUNDING_FLOOR_PER_SD) & (m.p_full < .05)).sum())}, subgroup "
    f"{int(((m.hr_hab > CONFOUNDING_FLOOR_PER_SD) & (m.p_hab < .05)).sum())}")
HEAD = ["cvd", "death", "htn2", "diabetes", "hf", "dementia"]
for k in HEAD:
    r = m[m.key == k]
    if len(r):
        r = r.iloc[0]
        say(f"    {LABEL[k]:28s} full {r.hr_full:.3f} [{r.lo_full:.3f}-{r.hi_full:.3f}]"
            f"   subgroup {r.hr_hab:.3f} [{r.lo_hab:.3f}-{r.hi_hab:.3f}]  "
            f"(events {int(r.events_full)} -> {int(r.events_hab)})")
RES["attack1"] = dict(
    n_smd_over_010=int((smd.abs_smd > 0.10).sum()), n_smd_vars=len(smd),
    max_smd_var=str(smd.iloc[0]["var"]), max_smd=float(smd.iloc[0].abs_smd),
    has_habitual_significant=int((sel.p < .05).sum()), has_habitual_bh=int((sel.q < .05).sum()),
    has_habitual_median_hr=float(sel.hr.median()),
    t90_logHR_r=float(np.corrcoef(lf, lh)[0, 1]), t90_slope=float(slope),
    t90_sign_agree=int((np.sign(lf) == np.sign(lh)).sum()), t90_n=int(len(m)),
    t90_sig_full=int((m.p_full < .05).sum()), t90_sig_sub=int((m.p_hab < .05).sum()),
    t90_median_full=float(m.hr_full.median()), t90_median_sub=float(m.hr_hab.median()))

# ======================================================================================
# 2  MEASUREMENT QUALITY
# ======================================================================================
hdr("ATTACK 2  MEASUREMENT QUALITY.  Does the answer flip across tier A, tier B, tier C and the "
    "template field?")

SETS = {
    "primary_8711": HABSET,
    "tierA": HABSET[HABSET.tier == "A"],
    "tierB": HABSET[HABSET.tier == "B"],
    "tierC": HABSET[HABSET.tier == "C"],
    "template_field": HABSET[HABSET.from_avg_nightly_hours_field == 1],
    "pre_study": PRESET,
}
tier_rows = []
for tag, d in SETS.items():
    d = add_exposures(d)
    say(f"  {tag}: n={len(d):,}  median hours {d.hours.median():.2f}  SD {d.hours.std():.2f}")
    for k in OUTCOMES:
        f = frame(d, k)
        r = cox(f, ["z_hab"], strata=(d.site_id.nunique() > 1))
        if r and "z_hab" in r:
            tier_rows.append(dict(set=tag, n_set=len(d), key=k, disease=LABEL[k],
                                  control=k in NEGATIVE_CONTROLS, circular=k in CIRCULAR,
                                  n=r["n"], events=r["events"], **r["z_hab"]))
tier = pd.DataFrame(tier_rows)
tier["mdhr80"] = np.exp(2.802 * tier.se)
tier.to_csv(f"{OUT}/attack_tier.csv", index=False)

say("")
say(f"{'set':16s} {'n':>6s} {'fitted':>6s} {'sig':>4s} {'medHR':>7s} {'maxHR':>7s} "
    f"{'ctrl_med':>8s}")
for tag in SETS:
    t = tier[(tier.set == tag) & (~tier.circular)]
    if not len(t):
        continue
    nc = t[~t.control]
    ct = t[t.control]
    say(f"{tag:16s} {len(SETS[tag]):>6,} {len(t):>6d} {int((t.p<.05).sum()):>4d} "
        f"{nc.hr.median():>7.3f} {nc.hr.sub(1).abs().add(1).max():>7.3f} "
        f"{ct.hr.median() if len(ct) else np.nan:>8.3f}")

SURV = ["htn2", "diabetes", "obesity", "gerd"]
say("")
say("the four claimed survivors, and the controls, across measurement quality "
    "(HR per SD of habitual hours):")
piv = tier[tier.key.isin(SURV + NEGATIVE_CONTROLS)].pivot_table(
    index="key", columns="set", values="hr")
piv = piv.reindex([k for k in SURV + NEGATIVE_CONTROLS if k in piv.index])
piv.index = [LABEL[k] + ("  [CONTROL]" if k in NEGATIVE_CONTROLS else "") for k in piv.index]
say(piv[[c for c in ["primary_8711", "pre_study", "tierB", "tierA", "tierC",
                     "template_field"] if c in piv.columns]]
    .to_string(float_format=lambda x: f"{x:.3f}"))

# formal heterogeneity across tiers A / B / C, per outcome
het = []
for k in NONCIRC:
    sub = tier[(tier.key == k) & (tier.set.isin(["tierA", "tierB", "tierC"]))]
    if len(sub) < 2:
        continue
    w = 1 / sub.se.values ** 2
    lhr = np.log(sub.hr.values)
    mu = (w * lhr).sum() / w.sum()
    Q = float((w * (lhr - mu) ** 2).sum())
    df = len(sub) - 1
    het.append(dict(key=k, disease=LABEL[k], control=k in NEGATIVE_CONTROLS, k_tiers=len(sub),
                    Q=Q, df=df, p_het=float(stats.chi2.sf(Q, df)),
                    I2=max(0.0, 100 * (Q - df) / Q) if Q > 0 else 0.0,
                    pooled_hr=float(np.exp(mu))))
het = pd.DataFrame(het)
het["q_het"] = bh(het.p_het.values)
het.to_csv(f"{OUT}/attack_tier_heterogeneity.csv", index=False)
say("")
say(f"between-tier heterogeneity (Cochran Q on log HR across tiers A/B/C), {len(het)} outcomes: "
    f"{int((het.p_het < .05).sum())} at p<0.05, {int((het.q_het < .05).sum())} after "
    f"Benjamini-Hochberg (2.5 expected by chance)")
say("  the four survivors: " + ", ".join(
    f"{LABEL[r.key]} Q p={r.p_het:.2f}" for r in het[het.key.isin(SURV)].itertuples()))

# tier A versus tier B agreement across the whole panel
ab = tier[tier.set == "tierA"].merge(tier[tier.set == "tierB"], on="key", suffixes=("_A", "_B"))
ab = ab[~ab.circular_A]
say(f"  tier A vs tier B log HR: Pearson r {np.corrcoef(np.log(ab.hr_A), np.log(ab.hr_B))[0,1]:.3f}"
    f", same direction {int((np.sign(np.log(ab.hr_A))==np.sign(np.log(ab.hr_B))).sum())} "
    f"of {len(ab)}")
RES["attack2"] = dict(
    tier_sets={t: int(len(SETS[t])) for t in SETS},
    het_sig=int((het.p_het < .05).sum()), het_bh=int((het.q_het < .05).sum()), het_n=len(het),
    tierAB_r=float(np.corrcoef(np.log(ab.hr_A), np.log(ab.hr_B))[0, 1]),
    tierAB_sign_agree=int((np.sign(np.log(ab.hr_A)) == np.sign(np.log(ab.hr_B))).sum()),
    tierAB_n=int(len(ab)))

# ======================================================================================
# 3  REVERSE CAUSATION
# ======================================================================================
hdr("ATTACK 3  REVERSE CAUSATION.  When was the note written, relative to the study and to the "
    "diagnosis?")

hs = HABSET.copy()
hs["note_days"] = (pd.to_datetime(hs.note_date_first) -
                   pd.to_datetime(hs.psg_date).dt.tz_localize(None)).dt.days
say(f"days from sleep study to the FIRST note used: median {hs.note_days.median():.0f}, "
    f"quartiles {hs.note_days.quantile(.25):.0f} and {hs.note_days.quantile(.75):.0f}")
say(f"  notes used that begin AFTER the study: "
    f"{int((hs.note_days > 0).sum()):,} of {len(hs):,} ({100*(hs.note_days>0).mean():.1f}%)")

rows = []
for k in OUTCOMES:
    f = frame(hs, k)
    ev = f[f.E == 1].copy()
    if len(ev) < MINEV:
        continue
    ev["dx_days"] = ev["T"] * 365.25
    rows.append(dict(key=k, disease=LABEL[k], events=len(ev),
                     note_after_dx=int((ev.note_days > ev.dx_days).sum()),
                     pct_note_after_dx=100 * float((ev.note_days > ev.dx_days).mean()),
                     note_after_psg=int((ev.note_days > 0).sum()),
                     median_dx_days=float(ev.dx_days.median())))
tim = pd.DataFrame(rows).sort_values("pct_note_after_dx", ascending=False)
tim.to_csv(f"{OUT}/attack_note_timing.csv", index=False)
say("")
say(f"among INCIDENT CASES, the habitual note was written after the diagnosis in a median of "
    f"{tim.pct_note_after_dx.median():.1f}% of cases (range {tim.pct_note_after_dx.min():.1f} to "
    f"{tim.pct_note_after_dx.max():.1f}%) across {len(tim)} outcomes")
say(tim.head(8)[["disease", "events", "pct_note_after_dx"]].to_string(
    index=False, float_format=lambda x: f"{x:.1f}"))
for k in ["htn2", "diabetes", "obesity", "gerd", "cvd", "death"]:
    r = tim[tim.key == k]
    if len(r):
        say(f"    {LABEL[k]:28s} {r.iloc[0].pct_note_after_dx:5.1f}% of "
            f"{int(r.iloc[0].events)} cases had the exposure recorded after diagnosis")

# landmark: drop anyone whose event or censoring happened before L, restart the clock at L
lm_rows = []
for L in (0.0, 1.0, 2.0):
    for tag, dd in (("habitual_8711", HABSET), ("pre_study_4631", PRESET)):
        d = add_exposures(dd)
        for expo in ("z_hab", "z_t90"):
            for k in OUTCOMES:
                f = frame(d, k)
                if L > 0:
                    f = f[f["T"] > L].copy()
                    f["T"] = f["T"] - L
                r = cox(f, [expo], strata=(f.site_id.nunique() > 1))
                if r and expo in r:
                    lm_rows.append(dict(landmark_years=L, set=tag, exposure=expo, key=k,
                                        disease=LABEL[k], control=k in NEGATIVE_CONTROLS,
                                        circular=k in CIRCULAR, n=r["n"], events=r["events"],
                                        **r[expo]))
lm = pd.DataFrame(lm_rows)
lm["mdhr80"] = np.exp(2.802 * lm.se)
lm.to_csv(f"{OUT}/attack_landmark.csv", index=False)

say("")
say(f"{'set':16s} {'expo':7s} {'LM':>4s} {'fit':>4s} {'sig':>4s} {'medHR':>7s}  survivors "
    f"(hypertension / diabetes / obesity / reflux)")
for tag in ("habitual_8711", "pre_study_4631"):
    for expo in ("z_hab", "z_t90"):
        for L in (0.0, 1.0, 2.0):
            s = lm[(lm.set == tag) & (lm.exposure == expo) & (lm.landmark_years == L) &
                   (~lm.circular)]
            if not len(s):
                continue
            nc = s[~s.control]
            vals = " ".join(f"{s[s.key==k].hr.iloc[0]:.3f}" if len(s[s.key == k]) else "  .  "
                            for k in SURV)
            say(f"{tag:16s} {expo:7s} {L:>4.0f} {len(s):>4d} {int((s.p<.05).sum()):>4d} "
                f"{nc.hr.median():>7.3f}  {vals}")

# post-study versus pre-study, same model, directly compared
ps_rows = []
for k in OUTCOMES:
    d = add_exposures(HABSET)
    d["postflag"] = d.post_study_only
    f = frame(d, k)
    r = cox(f, ["z_hab"], strata=True)
    fpre = f[f.pre_study_only == 1]
    fpost = f[f.post_study_only == 1]
    rp = cox(fpre, ["z_hab"], strata=True)
    rq = cox(fpost, ["z_hab"], strata=True)
    if r and rp and rq and all("z_hab" in x for x in (r, rp, rq)):
        d1, d2 = np.log(rp["z_hab"]["hr"]), np.log(rq["z_hab"]["hr"])
        se = np.sqrt(rp["z_hab"]["se"] ** 2 + rq["z_hab"]["se"] ** 2)
        ps_rows.append(dict(key=k, disease=LABEL[k], control=k in NEGATIVE_CONTROLS,
                            circular=k in CIRCULAR,
                            hr_all=r["z_hab"]["hr"], hr_pre=np.exp(d1), hr_post=np.exp(d2),
                            ev_pre=rp["events"], ev_post=rq["events"],
                            p_diff=float(2 * stats.norm.sf(abs(d1 - d2) / se))))
ps = pd.DataFrame(ps_rows)
ps.to_csv(f"{OUT}/attack_pre_vs_post.csv", index=False)
say("")
say(f"pre-study-only vs post-study-only habitual HR, same model, {len(ps)} outcomes: "
    f"{int((ps.p_diff < .05).sum())} differ at p<0.05 "
    f"(2.5 expected). median HR pre {ps.hr_pre.median():.3f}, post {ps.hr_post.median():.3f}")
RES["attack3"] = dict(
    pct_notes_after_psg=float(100 * (hs.note_days > 0).mean()),
    median_pct_note_after_dx=float(tim.pct_note_after_dx.median()),
    max_pct_note_after_dx=float(tim.pct_note_after_dx.max()),
    prepost_differ=int((ps.p_diff < .05).sum()), prepost_n=int(len(ps)))

# ======================================================================================
# 4  NOTE COUNT
# ======================================================================================
hdr("ATTACK 4  THE NOTE IS NOT A RANDOM SAMPLE.  Does note count predict outcomes, and does "
    "adjusting for it move the habitual estimates?")

d = add_exposures(HABSET)
d["log_notes"] = np.log1p(d.n_records)
d["z_notes"] = (d.log_notes - d.log_notes.mean()) / d.log_notes.std()
d["z_prev"] = (d.prevalent_count - d.prevalent_count.mean()) / d.prevalent_count.std()
say(f"duration-bearing notes per patient: median {d.n_records.median():.0f}, "
    f"quartiles {d.n_records.quantile(.25):.0f} and {d.n_records.quantile(.75):.0f}, "
    f"max {d.n_records.max():.0f}")
say(f"correlation of log note count with habitual hours: "
    f"Pearson {np.corrcoef(d.log_notes, d.hours)[0,1]:+.4f}, "
    f"Spearman {stats.spearmanr(d.log_notes, d.hours).statistic:+.4f}")
say(f"correlation of log note count with prevalent condition count: "
    f"Spearman {stats.spearmanr(d.log_notes, d.prevalent_count).statistic:+.4f}")

nc_rows = []
for k in OUTCOMES:
    f = frame(d, k)
    r1 = cox(f, ["z_notes"])
    r2 = cox(f, ["z_hab"])
    r3 = cox(f, ["z_hab", "z_notes"])
    r4 = cox(f, ["z_hab", "z_notes", "z_prev"])
    if not (r1 and r2 and r3 and r4):
        continue
    nc_rows.append(dict(
        key=k, disease=LABEL[k], control=k in NEGATIVE_CONTROLS, circular=k in CIRCULAR,
        events=r2["events"], n=r2["n"],
        notes_hr=r1["z_notes"]["hr"], notes_lo=r1["z_notes"]["lo"], notes_hi=r1["z_notes"]["hi"],
        notes_p=r1["z_notes"]["p"],
        hab_hr=r2["z_hab"]["hr"], hab_p=r2["z_hab"]["p"],
        hab_hr_adj_notes=r3["z_hab"]["hr"], hab_p_adj_notes=r3["z_hab"]["p"],
        hab_hr_adj_notes_prev=r4["z_hab"]["hr"], hab_p_adj_notes_prev=r4["z_hab"]["p"]))
ncf = pd.DataFrame(nc_rows)
ncf["notes_q"] = bh(ncf.notes_p.values)
ncf["shift_pct"] = 100 * (np.log(ncf.hab_hr_adj_notes) / np.log(ncf.hab_hr) - 1)
ncf.to_csv(f"{OUT}/attack_notecount.csv", index=False)
nn = ncf[~ncf.circular]
say("")
say(f"NOTE COUNT AS AN EXPOSURE, per SD of log note count, {len(nn)} non-circular outcomes:")
say(f"  significant for {int((nn.notes_p < .05).sum())}, "
    f"{int((nn.notes_q < .05).sum())} survive Benjamini-Hochberg; median HR "
    f"{nn.notes_hr.median():.3f}, {int((nn.notes_hr > 1).sum())} of {len(nn)} above 1")
say(f"  largest: " + ", ".join(
    f"{LABEL[r.key]} {r.notes_hr:.2f}" for r in
    nn.reindex(nn.notes_hr.sub(1).abs().sort_values(ascending=False).index).head(6).itertuples()))
say(f"  note count is significant for {int((nn[nn.control].notes_p < .05).sum())} of "
    f"{int(nn.control.sum())} NEGATIVE CONTROLS "
    f"(HRs {', '.join(f'{x:.2f}' for x in nn[nn.control].notes_hr)})")
say("")
say("habitual HR before and after adjusting for note count, and for note count plus prevalent "
    "burden:")
for k in SURV + NEGATIVE_CONTROLS:
    r = ncf[ncf.key == k]
    if len(r):
        r = r.iloc[0]
        say(f"    {LABEL[k]:24s}{'  [CONTROL]' if k in NEGATIVE_CONTROLS else '           '} "
            f"{r.hab_hr:.3f} -> {r.hab_hr_adj_notes:.3f} -> {r.hab_hr_adj_notes_prev:.3f}")
say(f"  median absolute shift in log HR from the note-count adjustment: "
    f"{nn.shift_pct.abs().median():.1f}%")
RES["attack4"] = dict(
    notes_sig=int((nn.notes_p < .05).sum()), notes_bh=int((nn.notes_q < .05).sum()),
    notes_n=int(len(nn)), notes_median_hr=float(nn.notes_hr.median()),
    notes_sig_controls=int((nn[nn.control].notes_p < .05).sum()),
    median_shift_pct=float(nn.shift_pct.abs().median()),
    corr_notes_hours=float(np.corrcoef(d.log_notes, d.hours)[0, 1]))

# ======================================================================================
# 5  POWER
# ======================================================================================
hdr("ATTACK 5  POWER.  Minimum detectable hazard ratio per cell.  Which nulls are uninformative?")

d = add_exposures(HABSET)
d["cell"] = d.hab_cat.astype(str) + " x " + d.o2_cat.astype(str)
ref = "6-<7h(ref) x normal T90<=1%"
cells = [c for c in d.cell.unique() if c != ref]
say(f"reference cell {ref}: N={int((d.cell == ref).sum()):,}")
say("cell sizes: " + ", ".join(f"{c} {int((d.cell==c).sum()):,}" for c in sorted(cells)))

dum = pd.get_dummies(d.cell)
NAMEMAP = {}
for c in cells:
    nm = "C" + str(abs(hash(c)) % 100000)
    NAMEMAP[nm] = c
    d[nm] = dum[c].astype(int).values
cellcols = list(NAMEMAP)

pw_rows = []
for k in OUTCOMES:
    f = frame(d, k)
    r = cox(f, cellcols)
    if not r:
        continue
    for nm, c in NAMEMAP.items():
        if nm not in r:
            continue
        x = r[nm]
        ncell = int((f.cell == c).sum())
        ecell = int(f[f.cell == c].E.sum())
        pw_rows.append(dict(key=k, disease=LABEL[k], control=k in NEGATIVE_CONTROLS,
                            circular=k in CIRCULAR, cell=c, n_cell=ncell, events_cell=ecell,
                            hr=x["hr"], lo=x["lo"], hi=x["hi"], p=x["p"], se=x["se"],
                            mdhr80=float(np.exp(2.802 * x["se"])),
                            mdhr50=float(np.exp(1.96 * x["se"]))))
pw = pd.DataFrame(pw_rows)
pw["informative"] = np.where(pw.p < .05, "significant",
                             np.where(pw.mdhr80 < 1.5, "null, informative",
                                      np.where(pw.mdhr80 < 2.0, "null, weak", "UNINFORMATIVE")))
pw.to_csv(f"{OUT}/attack_power_cells.csv", index=False)
say("")
say(f"{'cell':38s} {'N':>6s} {'medMDHR80':>10s} {'medMDHR50':>10s} {'sig':>5s} {'unin':>5s}")
for c in sorted(cells, key=lambda x: -int((d.cell == x).sum())):
    s = pw[(pw.cell == c) & (~pw.circular)]
    if not len(s):
        continue
    say(f"{c:38s} {int((d.cell==c).sum()):>6,} {s.mdhr80.median():>10.2f} "
        f"{s.mdhr50.median():>10.2f} {int((s.p<.05).sum()):>5d} "
        f"{int((s.informative=='UNINFORMATIVE').sum()):>5d}")
nn = pw[~pw.circular]
say(f"total cell contrasts fitted {len(nn)}: significant {int((nn.p<.05).sum())}, "
    f"informative nulls {int((nn.informative=='null, informative').sum())}, "
    f"weak nulls {int((nn.informative=='null, weak').sum())}, "
    f"UNINFORMATIVE {int((nn.informative=='UNINFORMATIVE').sum())}")
say("MDHR80 = exp(2.802 x SE), the hazard ratio this contrast had 80% power to detect at "
    "alpha 0.05. MDHR50 = exp(1.96 x SE), the bare edge of significance.")
RES["attack5"] = dict(
    contrasts=int(len(nn)), significant=int((nn.p < .05).sum()),
    uninformative=int((nn.informative == "UNINFORMATIVE").sum()),
    informative_nulls=int((nn.informative == "null, informative").sum()),
    median_mdhr80=float(nn.mdhr80.median()))

# ======================================================================================
# 6  HEAD TO HEAD
# ======================================================================================
hdr("ATTACK 6  HEAD TO HEAD.  Habitual home sleep versus laboratory total sleep time, in the "
    "same patients.")

d = add_exposures(HABSET)
say(f"correlation of habitual hours with laboratory TST inside the modelled set: "
    f"Pearson {np.corrcoef(d.hours, d.TST_min.fillna(d.TST_min.median()))[0,1]:+.4f}, "
    f"Spearman {stats.spearmanr(d.hours, d.TST_min, nan_policy='omit').statistic:+.4f}")
say("median laboratory sleep hours by habitual band: " + ", ".join(
    f"{b} {g.TST_min.median()/60:.2f}" for b, g in d.groupby("hab_cat", observed=True)))

h2h_rows = []
for expo in ("z_hab", "z_hab_rint", "z_tst", "z_t90", "z_ahi"):
    for k in OUTCOMES:
        f = frame(d, k)
        r = cox(f, [expo])
        if r and expo in r:
            h2h_rows.append(dict(exposure=expo, key=k, disease=LABEL[k],
                                 control=k in NEGATIVE_CONTROLS, circular=k in CIRCULAR,
                                 n=r["n"], events=r["events"], **r[expo]))
h2h = pd.DataFrame(h2h_rows)
h2h["q"] = h2h.groupby("exposure").p.transform(lambda x: bh(x.values))
h2h["mdhr80"] = np.exp(2.802 * h2h.se)
h2h["abs_hr"] = np.where(h2h.hr < 1, 1 / h2h.hr, h2h.hr)
h2h.to_csv(f"{OUT}/attack_headtohead.csv", index=False)

say("")
say(f"{'exposure':12s} {'fit':>4s} {'sig':>4s} {'BH':>4s} {'>floor':>7s} {'medHR':>7s} "
    f"{'maxHR':>7s} {'ctrl band':>18s}")
floors = {}
for expo in ("z_hab", "z_hab_rint", "z_tst", "z_t90", "z_ahi"):
    s = h2h[(h2h.exposure == expo) & (~h2h.circular)]
    ct = s[s.control]
    nc = s[~s.control]
    floor = float(ct.abs_hr.max())
    floors[expo] = floor
    nsig = int((nc.p < .05).sum())
    nbh = int((nc.q < .05).sum())
    nfl = int(((nc.p < .05) & (nc.abs_hr > floor)).sum())
    say(f"{expo:12s} {len(nc):>4d} {nsig:>4d} {nbh:>4d} {nfl:>7d} {nc.abs_hr.median():>7.3f} "
        f"{nc.abs_hr.max():>7.3f} {ct.hr.min():.3f}-{ct.hr.max():.3f}     "
        f"(empirical floor {floor:.3f})")
say("'>floor' = significant AND larger than the biggest negative-control effect in the SAME "
    "model. 'ctrl band' = the range the five controls occupy, which for a protective-direction "
    "exposure is the yardstick, not the two-sided floor.")

# cross-validated held-out concordance, my own 5-fold scheme, not the published one
def cv_concordance(d, expo, k=5, seed=17):
    rng = np.random.RandomState(seed)
    rows = []
    for key in OUTCOMES:
        f = frame(d, key)
        cols = ADJ + ["T", "E", "site_id"] + ([expo] if expo else [])
        f = f[cols].dropna()
        if f.E.sum() < 60:
            continue
        idx = rng.permutation(len(f))
        folds = np.array_split(idx, k)
        cs, ws = [], []
        for i in range(k):
            te = f.iloc[folds[i]]
            tr = f.iloc[np.concatenate([folds[j] for j in range(k) if j != i])]
            if tr.E.sum() < 20 or te.E.sum() < 5:
                continue
            try:
                c = CoxPHFitter().fit(tr.drop(columns=[]), "T", "E", strata=["site_id"])
                pred = c.predict_partial_hazard(te)
            except Exception:
                continue
            for site, g in te.groupby("site_id"):
                if g.E.sum() < 3:
                    continue
                try:
                    ci = concordance_index(g["T"], -pred.loc[g.index], g.E)
                except Exception:
                    continue
                cs.append(ci)
                ws.append(g.E.sum())
        if cs:
            rows.append(dict(key=key, c=float(np.average(cs, weights=ws)),
                             events=int(f.E.sum())))
    return pd.DataFrame(rows)


say("")
say("cross-validated held-out concordance, 5-fold, fresh resampling (NOT the published "
    "cross-hospital protocol) ...")
base = cv_concordance(d, None).rename(columns={"c": "c_base"})
conc = base.copy()
for expo in ("z_hab", "z_tst", "z_t90", "z_ahi"):
    e = cv_concordance(d, expo).rename(columns={"c": f"c_{expo}"})
    conc = conc.merge(e[["key", f"c_{expo}"]], on="key", how="left")
    conc[f"gain_{expo}"] = conc[f"c_{expo}"] - conc.c_base
conc["disease"] = conc.key.map(LABEL)
conc.to_csv(f"{OUT}/attack_concordance.csv", index=False)
say(f"outcomes scored: {len(conc)}   mean baseline C (age + sex, site-stratified): "
    f"{conc.c_base.mean():.4f}")
for expo in ("z_hab", "z_tst", "z_t90", "z_ahi"):
    g = conc[f"gain_{expo}"].dropna()
    say(f"  {expo:8s} mean gain {g.mean():+.6f}   median {g.median():+.6f}   "
        f"improves {int((g>0).sum())} of {len(g)}   best {g.max():+.4f}  worst {g.min():+.4f}")
gh, gt = conc.gain_z_hab.dropna(), conc.gain_z_tst.dropna()
both = conc.dropna(subset=["gain_z_hab", "gain_z_tst"])
w = stats.wilcoxon(both.gain_z_hab, both.gain_z_tst)
say(f"  habitual vs laboratory TST, paired Wilcoxon over {len(both)} outcomes: "
    f"p={w.pvalue:.3f} (habitual better in {int((both.gain_z_hab>both.gain_z_tst).sum())})")
bt = conc.dropna(subset=["gain_z_hab", "gain_z_t90"])
w3 = stats.wilcoxon(bt.gain_z_hab, bt.gain_z_t90)
say(f"  habitual vs laboratory T90, paired Wilcoxon over {len(bt)} outcomes: "
    f"p={w3.pvalue:.2e} (T90 better in {int((bt.gain_z_t90>bt.gain_z_hab).sum())})")

RES["attack6"] = dict(
    floors=floors,
    per_sd={expo: dict(
        sig=int(((h2h.exposure == expo) & (~h2h.circular) & (~h2h.control) &
                 (h2h.p < .05)).sum()),
        bh=int(((h2h.exposure == expo) & (~h2h.circular) & (~h2h.control) & (h2h.q < .05)).sum()),
        median_abs_hr=float(h2h[(h2h.exposure == expo) & (~h2h.circular) &
                                (~h2h.control)].abs_hr.median()))
        for expo in ("z_hab", "z_hab_rint", "z_tst", "z_t90", "z_ahi")},
    concordance={expo: dict(mean_gain=float(conc[f"gain_{expo}"].mean()),
                            improved=int((conc[f"gain_{expo}"] > 0).sum()),
                            n=int(conc[f"gain_{expo}"].notna().sum()))
                 for expo in ("z_hab", "z_tst", "z_t90", "z_ahi")},
    wilcoxon_hab_vs_tst=float(w.pvalue), wilcoxon_hab_vs_t90=float(w3.pvalue))

# ======================================================================================
hdr("DONE")
RES["_runtime_s"] = round(time.time() - T0, 1)
RES["_built"] = time.strftime("%Y-%m-%d %H:%M:%S")
with open(f"{OUT}/attack_summary.json", "w") as fh:
    json.dump(RES, fh, indent=2, default=str)
with open(f"{OUT}/attack_log.txt", "w") as fh:
    fh.write("\n".join(LOG))
say(f"runtime {RES['_runtime_s']}s")
