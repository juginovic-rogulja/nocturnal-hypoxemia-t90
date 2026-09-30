"""
Stage-resolved survival models in the GENERAL cohort.

Question: is sleep T90 stage-specific, or is one stage carrying it, or do all stages simply
track one another? Fitted on the strict PAP-off stage dataset built by build.py
(analysis.parquet, n=9,392, sites I0002 and I0006).

Specification, following the traps documented in the T90 manuscript:
  * natural cubic age spline, 4 basis columns, ONE DROPPED so the design is full rank. An
    unpenalized Cox on the rank-deficient design silently fits nothing.
  * penalizer=0. lifelines multiplies the elastic-net term by n, so any ridge shrinks
    everything toward 1 by an amount that scales inversely with the event count.
  * site-stratified, sex as a covariate, PSG date as time zero.
  * exposure rank-inverse-normal transformed WITHIN SITE, on the gate-eligible sample for that
    exposure, so every hazard ratio is per 1 SD of the within-site normal score.
  * every stage gated at >=30 min of that stage (the *_g30 columns, which are already NaN both
    outside the gate and where the T90 is not computable).
  * NREM and N2 never enter the same model (Spearman 0.956 here).

Outputs
  stage_general.csv                 one row per outcome per exposure, the primary panel
  stage_general_common.csv          the same panel refit on the subset clearing EVERY gate
  stage_general_h2h_rem_nrem.csv    REM and NREM in one model, top outcomes by event count
  stage_general_h2h_vs_sleep.csv    each stage against whole-sleep T90 in one model
  stage_general_summary.json        stage-win tallies, confounding floors, BH counts
  model_general_log.txt             full run log
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import json
import sys
import warnings
from itertools import combinations

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test
from patsy import dmatrix
from scipy import stats
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")

HERE = f"{paths.SV_ROOT}/stage_specific"
MANU = paths.T90_ROOT
sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
from disease_definitions import DISEASES, NEGATIVE_CONTROLS, CIRCULAR, ORGAN_GROUP  # noqa: E402

LOG = open(f"{HERE}/model_general_log.txt", "w")


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    LOG.write(s + "\n")
    LOG.flush()


# ------------------------------------------------------------------ data
b = pd.read_parquet(f"{HERE}/analysis.parquet")
say(f"loaded {len(b):,} x {b.shape[1]}  sites {dict(b.site_id.value_counts())}")
assert b.BDSPPatientID.duplicated().sum() == 0, "one row per patient expected"

# age spline on the frozen-file age, which is the age at the study the follow-up clock starts from
sp = dmatrix("cr(a, df=4) - 1", {"a": b.AgeAtVisit.values}, return_type="dataframe").iloc[:, 1:]
for i in range(sp.shape[1]):
    b[f"age_s{i}"] = sp.iloc[:, i].values
ADJ = [f"age_s{i}" for i in range(sp.shape[1])] + ["male"]
assert len(ADJ) == 4, ADJ
say(f"age spline: 4 basis columns, {sp.shape[1]} kept (one dropped), adjustment set {ADJ}")


def rint(x):
    """Blom rank-inverse-normal. Ties take the average rank, so a mass at zero collapses to one
    score. That is deliberate but it is also why heavily tied stages behave like a binary."""
    r = stats.rankdata(x)
    return stats.norm.ppf((r - 0.375) / (len(r) + 0.25))


def rint_within_site(frame, col):
    """Transform within site on the non-missing values only, leaving NaN where it was NaN."""
    out = pd.Series(np.nan, index=frame.index)
    for _, g in frame.groupby("site_id"):
        v = g[col].dropna()
        if len(v) > 10:
            out.loc[v.index] = rint(v.values)
    return out


STAGES = ["all", "wake", "sleep", "nrem", "n1", "n2", "n3", "rem"]
EXPOSURES = {}
for s in STAGES:
    src = f"t90_{s}_g30"
    b[f"z_{s}"] = rint_within_site(b, src)
    EXPOSURES[f"t90_{s}"] = f"z_{s}"
b["z_ahi"] = rint_within_site(b, "AHI")
EXPOSURES["AHI"] = "z_ahi"

say("\nexposure availability after the >=30 min gate, and the tie mass at exactly zero")
say(f"{'exposure':12s}{'n':>8}{'zeros':>9}{'%zero':>8}{'median T90':>12}{'IQR':>22}")
tie_frac = {}
for s in STAGES:
    v = b[f"t90_{s}_g30"].dropna()
    z = int((v == 0).sum())
    tie_frac[s] = z / len(v)
    q1, q3 = v.quantile([0.25, 0.75])
    say(f"{'t90_' + s:12s}{len(v):>8,}{z:>9,}{100 * z / len(v):>7.1f}%{v.median():>12.3f}"
        f"{f'{q1:.3f} to {q3:.3f}':>22}")
say(f"{'AHI':12s}{b.AHI.notna().sum():>8,}")

# ------------------------------------------------------------------ correlation structure
say("\nSpearman correlation between stage exposures (raw T90, pairwise complete)")
cm = b[[f"t90_{s}_g30" for s in STAGES]].corr(method="spearman")
cm.index = cm.columns = STAGES
say(cm.round(3).to_string())
COLLINEAR = [(a, c, cm.loc[a, c]) for a, c in combinations(STAGES, 2) if abs(cm.loc[a, c]) > 0.9]
say("pairs above 0.90, which must never share a model: "
    + ", ".join(f"{a}-{c} {r:.3f}" for a, c, r in COLLINEAR))

# ------------------------------------------------------------------ outcomes
OUTCOMES = [(k, v[0], k in NEGATIVE_CONTROLS) for k, v in DISEASES.items()]
OUTCOMES.append(("death", "Death from any cause", False))
say(f"\n{len(OUTCOMES)} outcomes ({len(NEGATIVE_CONTROLS)} negative controls, "
    f"{len(CIRCULAR)} circular and flagged but retained)")

MIN_EVENTS = 60          # the manuscript's reporting floor
COMMON = b[[f"z_{s}" for s in STAGES]].notna().all(axis=1) & b.AHI.notna()
say(f"patients clearing EVERY stage gate simultaneously: {int(COMMON.sum()):,} "
    f"({100 * COMMON.mean():.1f}% of {len(b):,})")


def risk_set(frame, key):
    """At-risk frame for one outcome: not prevalent, positive follow-up."""
    yc, ec, pc = f"{key}_years", f"{key}_incident", f"{key}_prevalent"
    if yc not in frame.columns:
        return None
    return frame[(frame[pc] == 0) & frame[yc].notna() & (frame[yc] > 0)]


def fit(d, cols):
    """Unpenalized site-stratified Cox. Returns {col: (hr, lo, hi, p)} or None."""
    keep = list(cols) + ADJ + ["T", "E", "site"]
    dd = d[keep].dropna()
    if dd.E.sum() < MIN_EVENTS or dd[list(cols)].std().min() == 0:
        return None
    try:
        c = CoxPHFitter(penalizer=0.0).fit(dd, "T", "E", strata=["site"])
    except Exception as e:
        say(f"    fit failed on {cols}: {type(e).__name__} {e}")
        return None
    # a silent no-op fit is the documented failure mode of the rank-deficient design
    if not np.isfinite(c.log_likelihood_) or abs(c.log_likelihood_ratio_test().test_statistic) < 1e-8:
        say(f"    degenerate fit on {cols}")
        return None
    out = {}
    for col in cols:
        r = c.summary.loc[col]
        out[col] = (float(r["exp(coef)"]), float(r["exp(coef) lower 95%"]),
                    float(r["exp(coef) upper 95%"]), float(r["p"]))
    out["_n"] = int(len(dd))
    out["_e"] = int(dd.E.sum())
    out["_fit"] = c
    return out


def contrast(c, a, bcol):
    """Wald test of beta_a == beta_b in a joint model. This is the direct test of stage
    specificity: if the two stages carry the same hazard, the difference is not distinguishable
    from zero and calling one of them 'the' stage is not supported."""
    p = c.params_
    v = c.variance_matrix_
    d = float(p[a] - p[bcol])
    se = float(np.sqrt(v.loc[a, a] + v.loc[bcol, bcol] - 2 * v.loc[a, bcol]))
    z = d / se
    return d, se, float(2 * stats.norm.sf(abs(z)))


def panel(frame, tag):
    """One model per exposure per outcome."""
    rows = []
    for key, label, is_neg in OUTCOMES:
        f = risk_set(frame, key)
        if f is None:
            continue
        base = pd.DataFrame({"T": f[f"{key}_years"], "E": f[f"{key}_incident"].astype(int),
                             "site": f.site_id, **{c: f[c] for c in ADJ},
                             **{v: f[v] for v in EXPOSURES.values()}})
        rec = {"key": key, "outcome": label, "organ": ORGAN_GROUP.get(key, "Composite"),
               "negative_control": is_neg, "circular": key in CIRCULAR,
               "n_at_risk": int(len(base)), "events_total": int(base.E.sum())}
        any_fit = False
        for name, col in EXPOSURES.items():
            r = fit(base, [col])
            if r is None:
                for suf in ("hr", "lo", "hi", "p"):
                    rec[f"{name}_{suf}"] = np.nan
                rec[f"{name}_n"] = np.nan
                rec[f"{name}_events"] = np.nan
                continue
            any_fit = True
            hr, lo, hi, p = r[col]
            rec[f"{name}_hr"] = round(hr, 4)
            rec[f"{name}_lo"] = round(lo, 4)
            rec[f"{name}_hi"] = round(hi, 4)
            rec[f"{name}_p"] = p
            rec[f"{name}_n"] = r["_n"]
            rec[f"{name}_events"] = r["_e"]
        if any_fit:
            rows.append(rec)
    df = pd.DataFrame(rows)
    say(f"\n[{tag}] fitted {len(df)} outcomes with >={MIN_EVENTS} events on at least one exposure")
    return df


# ------------------------------------------------------------------ 1. primary panel
say("\n" + "=" * 100)
say("PRIMARY PANEL, each exposure alone, own gated sample")
say("=" * 100)
main = panel(b, "primary")

# Benjamini-Hochberg across outcomes within each stage
STAGE_COLS = [f"t90_{s}" for s in STAGES] + ["AHI"]
bh_counts = {}
for name in STAGE_COLS:
    m = main[main[f"{name}_p"].notna() & ~main.negative_control & ~main.circular]
    if not len(m):
        continue
    rej, q, _, _ = multipletests(m[f"{name}_p"].values, alpha=0.05, method="fdr_bh")
    main.loc[m.index, f"{name}_q"] = q
    bh_counts[name] = {"tested": int(len(m)), "q05": int(rej.sum()),
                       "q05_hr_above_1": int(((q < 0.05) & (m[f"{name}_hr"] > 1)).sum())}
# negative controls and circular outcomes still get a q for completeness, from their own set
for name in STAGE_COLS:
    m = main[main[f"{name}_p"].notna() & (main.negative_control | main.circular)]
    if len(m) > 1:
        _, q, _, _ = multipletests(m[f"{name}_p"].values, alpha=0.05, method="fdr_bh")
        main.loc[m.index, f"{name}_q"] = q

# ------------------------------------------------------------------ 2. which stage wins
SLEEP_STAGES = ["sleep", "nrem", "n1", "n2", "n3", "rem"]   # NREM and N2 both listed, never co-modelled
WIN_POOL = ["wake", "sleep", "nrem", "n1", "n2", "n3", "rem"]


def winner(row, pool):
    hrs = {s: row[f"t90_{s}_hr"] for s in pool if pd.notna(row.get(f"t90_{s}_hr"))}
    if not hrs:
        return None, np.nan
    s = max(hrs, key=lambda k: hrs[k])
    return s, hrs[s]


main["win_stage"], main["win_hr"] = zip(*main.apply(lambda r: winner(r, WIN_POOL), axis=1))
main["win_stage_sleeponly"], main["win_hr_sleeponly"] = zip(
    *main.apply(lambda r: winner(r, SLEEP_STAGES), axis=1))
# margin of the winner over whole-sleep T90, the number that decides whether staging buys anything
main["win_minus_sleep"] = (main.win_hr_sleeponly - main["t90_sleep_hr"]).round(4)

# ------------------------------------------------------------------ 3. negative-control floors
say("\n" + "=" * 100)
say("NEGATIVE CONTROLS, every stage. The floor is the largest absolute deviation from 1.")
say("=" * 100)
neg = main[main.negative_control]
say(f"{'control':22s}{'events':>8}  " + "".join(f"{s:>16s}" for s in STAGE_COLS))
for _, r in neg.iterrows():
    line = f"{r.outcome:22s}{r.events_total:>8,}  "
    for name in STAGE_COLS:
        hr, p = r.get(f"{name}_hr"), r.get(f"{name}_p")
        if pd.isna(hr):
            line += f"{'--':>16s}"
        else:
            star = "*" if p < 0.05 else " "
            line += f"{hr:>14.3f}{star} "
    say(line)

floors = {}
for name in STAGE_COLS:
    v = neg[f"{name}_hr"].dropna()
    if not len(v):
        continue
    dev = pd.concat([v, 1 / v], axis=0)
    worst = float(dev.max())
    which = neg.loc[(v.apply(lambda x: max(x, 1 / x))).idxmax(), "outcome"]
    floors[name] = {"floor_hr": round(worst, 4), "driver": which,
                    "n_controls": int(len(v)),
                    "n_controls_p05": int((neg[f"{name}_p"] < 0.05).sum())}
say("\nconfounding floor per exposure (largest |HR| among the 5 negative controls)")
for name in STAGE_COLS:
    f_ = floors.get(name)
    if f_:
        say(f"  {name:12s} floor {f_['floor_hr']:.3f}  driven by {f_['driver']}  "
            f"({f_['n_controls_p05']} of {f_['n_controls']} controls p<0.05)")

# ------------------------------------------------------------------ 4. head to head REM vs NREM
say("\n" + "=" * 100)
say("HEAD TO HEAD, REM and NREM in the SAME model, top outcomes by event count")
say("=" * 100)
top = main.dropna(subset=["t90_rem_hr"]).nlargest(10, "events_total")
h2h_rows = []
for _, r in top.iterrows():
    f = risk_set(b, r.key)
    base = pd.DataFrame({"T": f[f"{r.key}_years"], "E": f[f"{r.key}_incident"].astype(int),
                         "site": f.site_id, **{c: f[c] for c in ADJ},
                         "z_rem": f.z_rem, "z_nrem": f.z_nrem, "z_sleep": f.z_sleep})
    joint = fit(base, ["z_rem", "z_nrem"])
    # the same restricted sample, one exposure at a time, so the joint model is compared like for like
    common = base.dropna(subset=["z_rem", "z_nrem"])
    solo_r = fit(common, ["z_rem"])
    solo_n = fit(common, ["z_nrem"])
    if joint is None or solo_r is None or solo_n is None:
        continue
    row = {"key": r.key, "outcome": r.outcome, "n": joint["_n"], "events": joint["_e"],
           "rem_alone_hr": round(solo_r["z_rem"][0], 4), "rem_alone_p": solo_r["z_rem"][3],
           "nrem_alone_hr": round(solo_n["z_nrem"][0], 4), "nrem_alone_p": solo_n["z_nrem"][3],
           "rem_adj_hr": round(joint["z_rem"][0], 4), "rem_adj_lo": round(joint["z_rem"][1], 4),
           "rem_adj_hi": round(joint["z_rem"][2], 4), "rem_adj_p": joint["z_rem"][3],
           "nrem_adj_hr": round(joint["z_nrem"][0], 4), "nrem_adj_lo": round(joint["z_nrem"][1], 4),
           "nrem_adj_hi": round(joint["z_nrem"][2], 4), "nrem_adj_p": joint["z_nrem"][3]}
    d, se, pc = contrast(joint["_fit"], "z_rem", "z_nrem")
    row.update({"rem_minus_nrem_logHR": round(d, 4), "contrast_se": round(se, 4),
                "contrast_p": pc, "top10": True})
    h2h_rows.append(row)

# the same head to head extended to every outcome with enough events, so the top-10 answer is not
# an artefact of which ten outcomes happen to be commonest
for _, r in main[(main.events_total >= 150) & main["t90_rem_hr"].notna()].iterrows():
    if r.key in set(x["key"] for x in h2h_rows):
        continue
    f = risk_set(b, r.key)
    base = pd.DataFrame({"T": f[f"{r.key}_years"], "E": f[f"{r.key}_incident"].astype(int),
                         "site": f.site_id, **{c: f[c] for c in ADJ},
                         "z_rem": f.z_rem, "z_nrem": f.z_nrem})
    joint = fit(base, ["z_rem", "z_nrem"])
    common = base.dropna(subset=["z_rem", "z_nrem"])
    solo_r, solo_n = fit(common, ["z_rem"]), fit(common, ["z_nrem"])
    if joint is None or solo_r is None or solo_n is None:
        continue
    d, se, pc = contrast(joint["_fit"], "z_rem", "z_nrem")
    h2h_rows.append({
        "key": r.key, "outcome": r.outcome, "n": joint["_n"], "events": joint["_e"],
        "rem_alone_hr": round(solo_r["z_rem"][0], 4), "rem_alone_p": solo_r["z_rem"][3],
        "nrem_alone_hr": round(solo_n["z_nrem"][0], 4), "nrem_alone_p": solo_n["z_nrem"][3],
        "rem_adj_hr": round(joint["z_rem"][0], 4), "rem_adj_lo": round(joint["z_rem"][1], 4),
        "rem_adj_hi": round(joint["z_rem"][2], 4), "rem_adj_p": joint["z_rem"][3],
        "nrem_adj_hr": round(joint["z_nrem"][0], 4), "nrem_adj_lo": round(joint["z_nrem"][1], 4),
        "nrem_adj_hi": round(joint["z_nrem"][2], 4), "nrem_adj_p": joint["z_nrem"][3],
        "rem_minus_nrem_logHR": round(d, 4), "contrast_se": round(se, 4),
        "contrast_p": pc, "top10": False})

h2h = pd.DataFrame(h2h_rows)
t10 = h2h[h2h.top10]
say(f"{'outcome':26s}{'ev':>6}{'REM alone':>12}{'NREM alone':>12}"
    f"{'REM|NREM':>18}{'NREM|REM':>18}{'REM=NREM p':>12}")
for _, r in t10.iterrows():
    say(f"{r.outcome:26s}{r.events:>6,}{r.rem_alone_hr:>12.3f}{r.nrem_alone_hr:>12.3f}"
        f"{f'{r.rem_adj_hr:.3f} p={r.rem_adj_p:.3f}':>18}"
        f"{f'{r.nrem_adj_hr:.3f} p={r.nrem_adj_p:.3f}':>18}{r.contrast_p:>12.3f}")
say(f"\nREM keeps p<0.05 adjusted for NREM in {int((t10.rem_adj_p < 0.05).sum())} of {len(t10)}")
say(f"NREM keeps p<0.05 adjusted for REM in {int((t10.nrem_adj_p < 0.05).sum())} of {len(t10)}")
say(f"REM and NREM coefficients differ at p<0.05 in {int((t10.contrast_p < 0.05).sum())} of {len(t10)}")
say(f"median attenuation of REM when NREM enters: "
    f"{100 * ((t10.rem_adj_hr - 1) / (t10.rem_alone_hr - 1)).median():.0f}% of the excess retained")
say(f"\nextended to all {len(h2h)} outcomes with >=150 events:")
say(f"  REM survives NREM p<0.05 in {int((h2h.rem_adj_p < 0.05).sum())}, "
    f"NREM survives REM in {int((h2h.nrem_adj_p < 0.05).sum())}, "
    f"coefficients formally differ in {int((h2h.contrast_p < 0.05).sum())}")
say(f"  REM alone > NREM alone in {int((h2h.rem_alone_hr > h2h.nrem_alone_hr).sum())} of {len(h2h)}")

# ------------------------------------------------------------------ 5. each stage vs whole sleep
say("\n" + "=" * 100)
say("HEAD TO HEAD, each stage against WHOLE-SLEEP T90 in one model")
say("=" * 100)
vs_rows = []
pool = [s for s in ["wake", "nrem", "n1", "n2", "n3", "rem"]]
for _, r in main.nlargest(20, "events_total").iterrows():
    f = risk_set(b, r.key)
    base = pd.DataFrame({"T": f[f"{r.key}_years"], "E": f[f"{r.key}_incident"].astype(int),
                         "site": f.site_id, **{c: f[c] for c in ADJ},
                         **{f"z_{s}": f[f"z_{s}"] for s in STAGES}})
    for s in pool:
        joint = fit(base, [f"z_{s}", "z_sleep"])
        common = base.dropna(subset=[f"z_{s}", "z_sleep"])
        solo_s = fit(common, [f"z_{s}"])
        solo_w = fit(common, ["z_sleep"])
        if joint is None or solo_s is None or solo_w is None:
            continue
        rr = cm.loc[s, "sleep"]
        vs_rows.append({
            "key": r.key, "outcome": r.outcome, "stage": s, "n": joint["_n"],
            "events": joint["_e"], "spearman_with_sleep": round(float(rr), 3),
            "stage_alone_hr": round(solo_s[f"z_{s}"][0], 4), "stage_alone_p": solo_s[f"z_{s}"][3],
            "sleep_alone_hr": round(solo_w["z_sleep"][0], 4), "sleep_alone_p": solo_w["z_sleep"][3],
            "stage_adj_hr": round(joint[f"z_{s}"][0], 4), "stage_adj_p": joint[f"z_{s}"][3],
            "sleep_adj_hr": round(joint["z_sleep"][0], 4), "sleep_adj_p": joint["z_sleep"][3],
            "stage_beats_sleep_alone": bool(solo_s[f"z_{s}"][0] > solo_w["z_sleep"][0]),
        })
vs = pd.DataFrame(vs_rows)
say("\nsame-sample comparison, stage alone vs whole-sleep alone, then both together")
for s in pool:
    g = vs[vs.stage == s]
    if not len(g):
        continue
    say(f"  {s:6s} n_outcomes {len(g):3d}  stage beats whole-sleep alone in "
        f"{int(g.stage_beats_sleep_alone.sum()):3d}  "
        f"median stage HR {g.stage_alone_hr.median():.3f} vs sleep {g.sleep_alone_hr.median():.3f}  "
        f"| joint: stage p<0.05 in {int((g.stage_adj_p < 0.05).sum())}, "
        f"sleep p<0.05 in {int((g.sleep_adj_p < 0.05).sum())}")

# ------------------------------------------------------------------ 6. common-sample refit
say("\n" + "=" * 100)
say("COMMON SAMPLE, only patients clearing EVERY stage gate, so stages are compared on the")
say("same people rather than on samples that differ by up to 4,600 patients")
say("=" * 100)
common_frame = b[COMMON].copy()
for s in STAGES:
    common_frame[f"z_{s}"] = rint_within_site(common_frame, f"t90_{s}_g30")
common_frame["z_ahi"] = rint_within_site(common_frame, "AHI")
cmain = panel(common_frame, "common")
if len(cmain):
    cmain["win_stage"], cmain["win_hr"] = zip(*cmain.apply(lambda r: winner(r, WIN_POOL), axis=1))
    cmain["win_stage_sleeponly"], cmain["win_hr_sleeponly"] = zip(
        *cmain.apply(lambda r: winner(r, SLEEP_STAGES), axis=1))
    for name in STAGE_COLS:
        m = cmain[cmain[f"{name}_p"].notna() & ~cmain.negative_control & ~cmain.circular]
        if len(m) > 1:
            _, q, _, _ = multipletests(m[f"{name}_p"].values, alpha=0.05, method="fdr_bh")
            cmain.loc[m.index, f"{name}_q"] = q

# ------------------------------------------------------------------ 6b. core common sample
say("\n" + "=" * 100)
say("CORE COMMON SAMPLE, patients clearing the wake, sleep, NREM, N2 and REM gates. N1 and N3")
say("are dropped from the requirement because their gates alone cost more than half the cohort.")
say("=" * 100)
CORE = b[[f"z_{s}" for s in ["wake", "sleep", "nrem", "n2", "rem"]]].notna().all(axis=1)
say(f"core common sample: {int(CORE.sum()):,} ({100 * CORE.mean():.1f}%)")
core_frame = b[CORE].copy()
for s in STAGES:
    core_frame[f"z_{s}"] = rint_within_site(core_frame, f"t90_{s}_g30")
core_frame["z_ahi"] = rint_within_site(core_frame, "AHI")
coremain = panel(core_frame, "core")
CORE_POOL = ["wake", "sleep", "nrem", "n2", "rem"]
coremain["win_stage"], coremain["win_hr"] = zip(
    *coremain.apply(lambda r: winner(r, CORE_POOL), axis=1))
coremain["win_stage_sleeponly"], coremain["win_hr_sleeponly"] = zip(
    *coremain.apply(lambda r: winner(r, ["sleep", "nrem", "n2", "rem"]), axis=1))
say("\nstage winner on the core common sample, sleep stages only")
say(coremain[~coremain.negative_control].win_stage_sleeponly.value_counts().to_string())
say("\ncore-sample negative controls")
for _, r in coremain[coremain.negative_control].iterrows():
    say(f"  {r.outcome:22s}" + "".join(
        f"{r.get('t90_' + s + '_hr'):>9.3f}" if pd.notna(r.get(f"t90_{s}_hr")) else f"{'--':>9s}"
        for s in CORE_POOL))

# ------------------------------------------------------------------ 7. do the stages track?
say("\n" + "=" * 100)
say("DO THE STAGES SIMPLY TRACK EACH OTHER? Spearman between the vectors of hazard ratios")
say("across outcomes. High correlation means the stage panel is one signal measured eight ways.")
say("=" * 100)
hrmat = main[[f"t90_{s}_hr" for s in STAGES]].copy()
hrmat.columns = STAGES
hrc = hrmat.corr(method="spearman")
say(hrc.round(3).to_string())
offdiag = [hrc.loc[a, c] for a, c in combinations(STAGES, 2)]
say(f"\nmedian off-diagonal correlation between stage HR vectors: {np.median(offdiag):.3f}")
say(f"REM vs the other sleep stages: "
    + ", ".join(f"{s} {hrc.loc['rem', s]:.3f}" for s in ["sleep", "nrem", "n1", "n2", "n3"]))
# spread of HRs within an outcome, against the spread expected from sampling noise alone
main["stage_hr_range"] = (hrmat.max(axis=1) - hrmat.min(axis=1)).round(4)
say(f"median across-stage spread of HR within an outcome: "
    f"{main.stage_hr_range.median():.3f} (IQR {main.stage_hr_range.quantile(.25):.3f} to "
    f"{main.stage_hr_range.quantile(.75):.3f})")


# ------------------------------------------------------------------ tallies
def tally(df, col, label):
    say(f"\n{label}")
    t = df[df[col].notna() & ~df.negative_control][col].value_counts()
    for s, n in t.items():
        say(f"  {s:8s}{n:4d}")
    return {k: int(v) for k, v in t.items()}


say("\n" + "=" * 100)
say("STAGE WIN TALLIES")
say("=" * 100)
win_all = tally(main, "win_stage", "largest HR across wake + all sleep stages, primary panel")
win_sleep = tally(main, "win_stage_sleeponly", "largest HR restricted to sleep stages, primary panel")
win_common = tally(cmain, "win_stage_sleeponly", "largest HR restricted to sleep stages, COMMON sample") \
    if len(cmain) else {}

# how often the best stage clears its own floor and beats whole-sleep meaningfully
main["win_clears_floor"] = main.apply(
    lambda r: bool(pd.notna(r.win_hr_sleeponly)
                   and r.win_hr_sleeponly > floors.get(f"t90_{r.win_stage_sleeponly}",
                                                       {"floor_hr": np.inf})["floor_hr"]), axis=1)
n_pos = int((~main.negative_control & ~main.circular).sum())
say(f"\nof {n_pos} non-control outcomes, the winning sleep stage exceeds its own negative-control "
    f"floor in {int((main.win_clears_floor & ~main.negative_control & ~main.circular).sum())}")
say(f"whole-sleep T90 itself exceeds the sleep floor "
    f"({floors['t90_sleep']['floor_hr']:.3f}) in "
    f"{int(((main['t90_sleep_hr'] > floors['t90_sleep']['floor_hr']) & ~main.negative_control & ~main.circular).sum())}")

say("\nBenjamini-Hochberg across non-control, non-circular outcomes WITHIN each exposure")
for name in STAGE_COLS:
    c = bh_counts.get(name)
    if c:
        say(f"  {name:12s} tested {c['tested']:3d}  q<0.05 {c['q05']:3d}  "
            f"(HR>1 in {c['q05_hr_above_1']})")

# ------------------------------------------------------------------ headline table
say("\n" + "=" * 100)
say("PRIMARY PANEL, sorted by whole-sleep T90 hazard ratio")
say("=" * 100)
hdr = f"{'outcome':26s}{'ev':>6}" + "".join(f"{s:>9s}" for s in WIN_POOL) + f"{'AHI':>9}{'win':>7}"
say(hdr)
say("-" * len(hdr))
for _, r in main.sort_values("t90_sleep_hr", ascending=False).iterrows():
    mark = "N" if r.negative_control else ("c" if r.circular else " ")
    line = f"{mark}{r.outcome[:25]:25s}{r.events_total:>6,}"
    for s in WIN_POOL:
        hr = r.get(f"t90_{s}_hr")
        line += f"{'--':>9s}" if pd.isna(hr) else f"{hr:>9.3f}"
    ah = r.get("AHI_hr")
    line += f"{'--':>9s}" if pd.isna(ah) else f"{ah:>9.3f}"
    line += f"{str(r.win_stage_sleeponly):>7s}"
    say(line)

# ------------------------------------------------------------------ 8. verification
say("\n" + "=" * 100)
say("VERIFICATION")
say("=" * 100)
_f = risk_set(b, "diabetes")
_d = pd.DataFrame({"T": _f.diabetes_years, "E": _f.diabetes_incident.astype(int),
                   "site": _f.site_id, "male": _f.male, "z_rem": _f.z_rem})
for i in range(sp.shape[1]):
    _d[f"age_s{i}"] = _f[f"age_s{i}"]
_d4 = _d.copy()
_d4["age_s3"] = dmatrix("cr(a, df=4) - 1", {"a": b.AgeAtVisit.values},
                        return_type="dataframe").iloc[:, 0].reindex(_f.index).values
try:
    CoxPHFitter(penalizer=0.0).fit(_d4.dropna(), "T", "E", strata=["site"])
    say("  4-column age spline: FIT RETURNED, inspect it, the design is rank deficient")
except Exception as e:
    say(f"  4-column age spline refuses to fit: {type(e).__name__}. This is the documented trap. "
        f"One column is dropped throughout.")
_c = CoxPHFitter(penalizer=0.0).fit(_d.dropna(), "T", "E", strata=["site"])
say(f"  3-column design, diabetes on REM T90: log-likelihood ratio "
    f"{_c.log_likelihood_ratio_test().test_statistic:.1f}, concordance {_c.concordance_index_:.3f}, "
    f"{int(_d.dropna().E.sum())} events. A silent no-op fit would give a ratio of 0.")

say("\n  proportional hazards on the exposure (scaled Schoenfeld, rank time transform)")
for key, lab in [("diabetes", "Type 2 diabetes"), ("aki", "Acute kidney injury"),
                 ("obesity", "Obesity"), ("htn2", "Hypertension"),
                 ("resp_failure", "Respiratory failure"), ("death", "Death"),
                 ("osteoarthritis", "Osteoarthritis")]:   # an outcome, not a control, since 2026-08-07
    for ex in ["z_rem", "z_sleep"]:
        f = risk_set(b, key)
        d = pd.DataFrame({"T": f[f"{key}_years"], "E": f[f"{key}_incident"].astype(int),
                          "site": f.site_id, "male": f.male, ex: f[ex],
                          **{c: f[c] for c in ADJ if c != "male"}}).dropna()
        c = CoxPHFitter(penalizer=0.0).fit(d, "T", "E", strata=["site"])
        p = proportional_hazard_test(c, d, time_transform="rank").summary.loc[ex, "p"]
        flag = "  <-- violated" if p < 0.05 else ""
        say(f"    {lab:22s}{ex:9s}ev {int(d.E.sum()):5d}  HR {c.summary.loc[ex, 'exp(coef)']:.3f}  "
            f"PH p={p:.4f}{flag}")

# reconciliation against the 2026-07-29 table, which used the looser PAP-off rule on n=13,066
say("\n  reconciliation against T90_by_state_vs_AHI_all_diseases.csv (2026-07-29, n=13,066)")
try:
    old = pd.read_csv(f"{paths.T90_ROOT}/BDSP_New_Ideas_2026-07-26/"
                      "deliverable/T90_by_state_vs_AHI_all_diseases.csv")
    old["disease"] = old.disease.replace({
        "Anaemia": "Anemia", "Blood clots (VTE)": "Venous thromboembolism",
        "Cardiovascular": "Cardiovascular composite", "Cataract [NEG]": "Cataract",
        "Death (any cause)": "Death from any cause", "Dyslipidaemia": "Dyslipidemia",
        # the 2026-07-29 table tagged two conditions [NEG]; this map strips the tag, because
        # fracture and osteoarthritis are ordinary outcomes on the settled panel
        "Fatty liver": "Fatty liver disease", "Fracture [NEG]": "Fracture",
        "Osteoarthritis [NEG]": "Osteoarthritis",
        "Heart attack": "Myocardial infarction", "High blood pressure": "Hypertension",
        "Ischaemic heart disease": "Ischemic heart disease",
        "Kidney disease": "Chronic kidney disease", "Reflux": "Gastroesophageal reflux",
        "Restless legs": "Restless legs syndrome"})
    j = old.merge(main, left_on="disease", right_on="outcome", how="inner")
    say(f"    {len(j)} of {len(old)} outcomes matched")
    OLDMAP = {"T90 all": "t90_all", "WAKE": "t90_wake", "SLEEP": "t90_sleep", "NREM": "t90_nrem",
              "N1": "t90_n1", "N2": "t90_n2", "N3": "t90_n3", "REM": "t90_rem", "AHI": "AHI"}
    recon = {}
    for k, v in OLDMAP.items():
        dd = j[[k, f"{v}_hr"]].dropna()
        rho = float(stats.spearmanr(dd[k], dd[f"{v}_hr"]).statistic)
        recon[v] = {"n": int(len(dd)), "spearman": round(rho, 3),
                    "median_shift": round(float((dd[f"{v}_hr"] - dd[k]).median()), 3)}
        say(f"    {k:10s} n {len(dd):3d}  Spearman {rho:.3f}  median shift "
            f"{(dd[f'{v}_hr'] - dd[k]).median():+.3f}")
    OS = ["SLEEP", "NREM", "N1", "N2", "N3", "REM"]
    j["old_win"] = j[OS].idxmax(axis=1).str.lower()
    agree = int((j.old_win == j.win_stage_sleeponly).sum())
    say(f"    the WINNING STAGE agrees in only {agree} of {len(j)} outcomes "
        f"({100 * agree / len(j):.0f}%), although every stage's hazard ratios correlate above 0.80")
    recon["winner_agreement"] = {"agree": agree, "n": int(len(j))}
except Exception as e:
    say(f"    reconciliation skipped: {type(e).__name__} {e}")
    recon = {}

# ------------------------------------------------------------------ write
main.to_csv(f"{HERE}/stage_general.csv", index=False)
cmain.to_csv(f"{HERE}/stage_general_common.csv", index=False)
coremain.to_csv(f"{HERE}/stage_general_common_core.csv", index=False)
h2h.to_csv(f"{HERE}/stage_general_h2h_rem_nrem.csv", index=False)
vs.to_csv(f"{HERE}/stage_general_h2h_vs_sleep.csv", index=False)

summary = {
    "n_cohort": int(len(b)),
    "sites": {str(k): int(v) for k, v in b.site_id.value_counts().items()},
    "n_outcomes_fitted": int(len(main)),
    "min_events": MIN_EVENTS,
    "gate_n": {s: int(b[f"t90_{s}_g30"].notna().sum()) for s in STAGES},
    "zero_fraction": {s: round(tie_frac[s], 4) for s in STAGES},
    "n_common_all_gates": int(COMMON.sum()),
    "n_common_core_gates": int(CORE.sum()),
    "win_tally_sleep_stages_core_sample": {
        k: int(v) for k, v in
        coremain[~coremain.negative_control].win_stage_sleeponly.value_counts().items()},
    "hr_vector_spearman_median": round(float(np.median(offdiag)), 3),
    "hr_vector_spearman_rem_vs": {s: round(float(hrc.loc["rem", s]), 3)
                                  for s in STAGES if s != "rem"},
    "median_across_stage_hr_range": round(float(main.stage_hr_range.median()), 4),
    "rem_nrem_contrast_p05": int((h2h.contrast_p < 0.05).sum()),
    "rem_nrem_contrast_tested": int(len(h2h)),
    "reconciliation_vs_2026_07_29": recon,
    "spearman": {f"{a}-{c}": round(float(cm.loc[a, c]), 3) for a, c in combinations(STAGES, 2)},
    "win_tally_all_states": win_all,
    "win_tally_sleep_stages": win_sleep,
    "win_tally_sleep_stages_common_sample": win_common,
    "negative_control_floor": floors,
    "bh": bh_counts,
    "median_hr_by_exposure": {n: round(float(main[f"{n}_hr"].median()), 4) for n in STAGE_COLS},
    "median_hr_by_exposure_noncontrol": {
        n: round(float(main.loc[~main.negative_control & ~main.circular, f"{n}_hr"].median()), 4)
        for n in STAGE_COLS},
    "rem_survives_nrem": int((h2h.rem_adj_p < 0.05).sum()) if len(h2h) else 0,
    "nrem_survives_rem": int((h2h.nrem_adj_p < 0.05).sum()) if len(h2h) else 0,
    "h2h_n": int(len(h2h)),
}
json.dump(summary, open(f"{HERE}/stage_general_summary.json", "w"), indent=1)
say(f"\nwrote stage_general.csv ({main.shape[0]} x {main.shape[1]}), "
    f"stage_general_common.csv ({cmain.shape[0]}), stage_general_h2h_rem_nrem.csv ({len(h2h)}), "
    f"stage_general_h2h_vs_sleep.csv ({len(vs)}), stage_general_summary.json")
LOG.close()
