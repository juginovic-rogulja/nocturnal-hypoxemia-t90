"""
Hostile re-analysis of the claim that sleep T90's association with disease is STAGE-SPECIFIC.

Everything here is recomputed from analysis.parquet (and, for reliability, from the raw
per-stage pilot shards). Nothing is taken from stage_general.csv or stage_healthy.csv except
for a final reconciliation of point estimates.

Five attacks, each ending in HOLDS / WEAKENED / REFUTED:
  1 collinearity + bootstrap stability of the winning stage
  2 differential gating: does a stage win because its gate selected a sicker subgroup
  3 reliability ceiling: night-to-night ICC per stage, disattenuation
  4 negative controls recomputed per stage, confounding floor
  5 is whole-sleep T90 sufficient: nested LR test and cross-validated concordance gain

Modelling follows the documented traps: unpenalized Cox (penalizer=0), site strata, sex,
natural cubic age spline df=4 with ONE COLUMN DROPPED, exposure rank-inverse-normal within
site, HR per 1 SD, each stage gated at >=30 min via the *_g30 columns, NREM and N2 never in
one model.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import glob
import json
import os
import sys
import time
import warnings
from itertools import combinations

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from patsy import dmatrix
from scipy import stats

warnings.filterwarnings("ignore")
np.random.seed(20260807)

HERE = f"{paths.SV_ROOT}/stage_specific"
TROOT = paths.T90_ROOT
PILOT = f"{paths.V8_ROOT}/X4_groupP/t90_by_stage"   # v8.2: re-extracted, collapse-masked shards (X4_groupP)
sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
from disease_definitions import DISEASES, NEGATIVE_CONTROLS, CIRCULAR  # noqa: E402

SKIP_BOOT = os.environ.get("SKIP_BOOT") == "1"
LOG = open(f"{HERE}/attack_log{'_noboot' if SKIP_BOOT else ''}.txt", "w")
T0 = time.time()


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.write(s + "\n")
    LOG.flush()


def head(t):
    say("\n" + "=" * 100)
    say(t)
    say("=" * 100)


B_BOOT = int(os.environ.get("B_BOOT", 250))
STAGES = ["all", "wake", "sleep", "nrem", "n1", "n2", "n3", "rem"]
SLEEPSTAGES = ["sleep", "nrem", "n1", "n2", "n3", "rem"]
CORE = ["wake", "sleep", "nrem", "n2", "rem"]  # stages with a workable common sample

# ---------------------------------------------------------------- load
b = pd.read_parquet(f"{HERE}/analysis.parquet")
say(f"analysis.parquet {b.shape[0]:,} x {b.shape[1]}   sites {dict(b.site_id.value_counts())}")
assert b.BDSPPatientID.duplicated().sum() == 0
assert (b.oximetry_bad == 0).all() and (b.fu_valid == 1).all() and b.papoff_strict.all()

sp = dmatrix("cr(a, df=4) - 1", {"a": b.AgeAtVisit.values}, return_type="dataframe").iloc[:, 1:]
for i in range(sp.shape[1]):
    b[f"age_s{i}"] = sp.iloc[:, i].values
ADJ = [f"age_s{i}" for i in range(sp.shape[1])] + ["male"]
assert len(ADJ) == 4


def rint(x):
    r = stats.rankdata(x)
    return stats.norm.ppf((r - 0.375) / (len(r) + 0.25))


def rint_site(frame, col):
    out = pd.Series(np.nan, index=frame.index)
    for _, g in frame.groupby("site_id"):
        v = g[col].dropna()
        if len(v) > 10:
            out.loc[v.index] = rint(v.values)
    return out


for s in STAGES:
    b[f"z_{s}"] = rint_site(b, f"t90_{s}_g30")
b["z_ahi"] = rint_site(b, "AHI")

OUTCOMES = [(k, v[0], k in NEGATIVE_CONTROLS) for k, v in DISEASES.items()]
OUTCOMES += [("death", "Death from any cause", False)]
OUTCOMES = [(k, lab, neg) for k, lab, neg in OUTCOMES if f"{k}_years" in b.columns]
NEG = [k for k, _, n in OUTCOMES if n]
MAIN = [k for k, _, n in OUTCOMES if not n and k not in CIRCULAR]
say(f"{len(OUTCOMES)} outcomes, {len(NEG)} negative controls {NEG}, {len(MAIN)} main")

# ---------------------------------------------------------------- provenance spot check
head("PROVENANCE: analysis.parquet against the raw pilot shards")
sh = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"{PILOT}/out_*.csv"))], ignore_index=True)
say(f"shards {len(sh):,} rows before dedupe, {sh.BDSPPatientID.nunique():,} patients")
sh = sh.drop_duplicates(subset=["BIDSFolder", "SessionID"])
say(f"{len(sh):,} rows after dedupe on BIDSFolder+SessionID")
m = b.merge(sh, on=["BIDSFolder", "SessionID"], suffixes=("", "_raw"))
say(f"{len(m):,} of {len(b):,} analysis rows matched back to a shard row")
for s in ["all", "sleep", "rem", "n3"]:
    d = (m[f"t90_{s}"] - m[f"t90_{s}_raw"]).abs()
    say(f"  t90_{s:5s} max abs diff vs shard {np.nanmax(d):.6g}, "
        f"mismatches>1e-6: {int((d > 1e-6).sum())}")
say(f"cpap_frac==0 in the matched rows: {(m.cpap_frac_raw == 0).mean():.4f}; "
    f"has_cpap_chan: {m.has_cpap_chan_raw.mean():.4f}")

# ================================================================================
head("ATTACK 1  COLLINEARITY, AND IS THE WINNING STAGE STABLE")
# ================================================================================
raw = b[[f"t90_{s}_g30" for s in STAGES]].copy()
raw.columns = STAGES
sp_c = raw.corr(method="spearman")
pe_c = raw.corr(method="pearson")
say("\nSpearman, pairwise complete, gated values (the analysis exposures)")
say(sp_c.round(3).to_string())
say("\nPearson on the rank-inverse-normal analysis scale")
zc = b[[f"z_{s}" for s in STAGES]].corr()
zc.index = zc.columns = STAGES
say(zc.round(3).to_string())
say("\ncomplete-case correlation, patients clearing every gate "
    f"(n={int(raw.dropna().shape[0]):,})")
say(raw.dropna().corr(method="spearman").round(3).to_string())
say("\nkey pairs (Spearman, gated pairwise):")
for a, c in [("rem", "nrem"), ("rem", "sleep"), ("nrem", "sleep"), ("n2", "nrem"),
             ("rem", "all"), ("wake", "sleep"), ("n3", "rem"), ("n1", "rem")]:
    say(f"  {a:5s} vs {c:5s}  {sp_c.loc[a, c]:.3f}")
hi = [(a, c, sp_c.loc[a, c]) for a, c in combinations(STAGES, 2) if sp_c.loc[a, c] > 0.8]
say(f"pairs above 0.80: {len(hi)} of {len(list(combinations(STAGES, 2)))}  ->  "
    + ", ".join(f"{a}-{c} {r:.2f}" for a, c, r in hi))
say("\nnon-zero-only Spearman (drop the shared tie mass at zero; this is the correlation of the")
say("informative part of each pair, i.e. how much independent information a stage could carry)")
for a, c in [("rem", "nrem"), ("rem", "sleep"), ("nrem", "sleep")]:
    sub = raw[[a, c]].dropna()
    sub = sub[(sub[a] > 0) & (sub[c] > 0)]
    say(f"  {a} vs {c}: n={len(sub):,}  rho={stats.spearmanr(sub[a], sub[c])[0]:.3f}")

# ---- core common sample ----------------------------------------------------
core = b[np.logical_and.reduce([b[f"gate_{s}"].values for s in CORE])
         & b[[f"t90_{s}_g30" for s in CORE]].notna().all(axis=1).values].copy()
for s in CORE:
    core[f"z_{s}"] = rint_site(core, f"t90_{s}_g30")
say(f"\ncore common sample (clears wake, sleep, NREM, N2 and REM): n={len(core):,} "
    f"({100 * len(core) / len(b):.1f}% of {len(b):,})")
allg = b[b[[f"t90_{s}_g30" for s in STAGES]].notna().all(axis=1)]
say(f"clears ALL eight gates: n={len(allg):,} ({100 * len(allg) / len(b):.1f}%)")


def fit(d, cols, strata=True):
    """Unpenalized site-stratified Cox. Returns dict of coef/se/p/loglik/concordance."""
    keep = list(cols) + ADJ + ["T", "E"] + (["site"] if strata else [])
    dd = d[keep].dropna()
    if dd.E.sum() < 20:
        return None
    c = CoxPHFitter(penalizer=0.0)
    c.fit(dd, "T", "E", strata=(["site"] if strata else None))
    out = {"n": len(dd), "events": int(dd.E.sum()), "loglik": c.log_likelihood_,
           "c": c.concordance_index_, "model": c}
    for col in cols:
        r = c.summary.loc[col]
        out[col] = {"hr": float(r["exp(coef)"]), "coef": float(r["coef"]),
                    "se": float(r["se(coef)"]), "p": float(r["p"]),
                    "lo": float(r["exp(coef) lower 95%"]), "hi": float(r["exp(coef) upper 95%"])}
    return out


def frame_for(src, key, zcols):
    f = src[(src[f"{key}_prevalent"] == 0) & src[f"{key}_years"].notna()
            & (src[f"{key}_years"] > 0)]
    d = pd.DataFrame({"T": f[f"{key}_years"].values, "E": f[f"{key}_incident"].astype(int).values,
                      "site": f.site_id.values},
                     index=f.index)
    for c in zcols + ADJ:
        d[c] = f[c].values
    return d


# ---- observed per-stage panel on the core sample ---------------------------
say("\nObserved hazard ratios on the CORE COMMON SAMPLE (identical patients, all five windows)")
say(f"{'outcome':22s}{'ev':>6}" + "".join(f"{s:>10}" for s in CORE) + f"{'winner':>10}{'spread':>9}")
core_rows = []
for k, lab, neg in OUTCOMES:
    d = frame_for(core, k, [f"z_{s}" for s in CORE])
    if d.E.sum() < 60:
        continue
    rec = {"key": k, "label": lab, "neg": neg, "events": int(d.E.sum()), "n": len(d)}
    for s in CORE:
        r = fit(d, [f"z_{s}"])
        if r:
            rec[f"{s}_hr"] = r[f"z_{s}"]["hr"]
            rec[f"{s}_p"] = r[f"z_{s}"]["p"]
            rec[f"{s}_lo"] = r[f"z_{s}"]["lo"]
            rec[f"{s}_hi"] = r[f"z_{s}"]["hi"]
            rec[f"{s}_coef"] = r[f"z_{s}"]["coef"]
            rec[f"{s}_se"] = r[f"z_{s}"]["se"]
            rec[f"{s}_ll"] = r["loglik"]
            rec[f"{s}_c"] = r["c"]
    sl = [s for s in CORE if s != "wake"]
    hrs = {s: rec.get(f"{s}_hr") for s in sl if rec.get(f"{s}_hr") == rec.get(f"{s}_hr")}
    rec["winner"] = max(hrs, key=hrs.get) if hrs else None
    rec["spread"] = (max(hrs.values()) - min(hrs.values())) if hrs else np.nan
    core_rows.append(rec)
    say(f"{lab[:21]:22s}{rec['events']:>6}"
        + "".join(f"{rec.get(s + '_hr', float('nan')):>10.3f}" for s in CORE)
        + f"{str(rec['winner']):>10}{rec['spread']:>9.3f}"
        + ("   [negative control]" if neg else ""))
core_df = pd.DataFrame(core_rows)
core_df.to_csv(f"{HERE}/attack_core_panel.csv", index=False)

sl = [s for s in CORE if s != "wake"]
med_ciw = {s: np.nanmedian(core_df[f"{s}_hi"] - core_df[f"{s}_lo"]) for s in sl}
say(f"\nmedian across-stage spread within an outcome: {np.nanmedian(core_df.spread):.3f}")
say("median 95% CI width of a single stage estimate: "
    + ", ".join(f"{s} {v:.3f}" for s, v in med_ciw.items()))
say(f"outcomes where the spread exceeds the REM CI width: "
    f"{int((core_df.spread > med_ciw['rem']).sum())} of {len(core_df)}")

# ---- bootstrap the identity of the winner -----------------------------------
BOOT_KEYS = [r["key"] for r in core_rows if not r["neg"] and r["key"] not in CIRCULAR
             and r["events"] >= 150]
BOOT_KEYS = sorted(BOOT_KEYS, key=lambda k: -core_df.set_index("key").loc[k, "events"])[:16]
NINE = ["aki", "gout", "diabetes", "pneumonia", "obesity", "cellulitis", "ckd", "htn2",
        "dyslipid"]  # the nine REM claims said to survive everything
TEST_KEYS = BOOT_KEYS + [k for k in NINE
                         if k not in BOOT_KEYS and k in set(core_df.key)]
say(f"\nbootstrapping the winner for {len(BOOT_KEYS)} outcomes (>=150 events), "
    f"B={B_BOOT}, patient-level resampling of the core sample")
idx = np.arange(len(core))
boot_win = {k: [] for k in BOOT_KEYS}
boot_hr = {k: {s: [] for s in sl} for k in BOOT_KEYS}
t = time.time()
for bi in range(0 if not SKIP_BOOT else B_BOOT, B_BOOT):
    take = np.random.randint(0, len(core), len(core))
    cb = core.iloc[take].copy()
    cb.index = np.arange(len(cb))
    for s in sl:
        cb[f"z_{s}"] = rint_site(cb, f"t90_{s}_g30")
    for k in BOOT_KEYS:
        d = frame_for(cb, k, [f"z_{s}" for s in sl])
        hrs = {}
        for s in sl:
            try:
                r = fit(d, [f"z_{s}"])
            except Exception:
                r = None
            if r:
                hrs[s] = r[f"z_{s}"]["hr"]
                boot_hr[k][s].append(r[f"z_{s}"]["hr"])
        if len(hrs) == len(sl):
            boot_win[k].append(max(hrs, key=hrs.get))
    if bi % 25 == 0:
        say(f"   boot {bi}/{B_BOOT}  {time.time() - t:.0f}s")
say(f"bootstrap done in {time.time() - t:.0f}s")

say(f"\n{'outcome':22s}{'ev':>6}{'observed':>10}{'modal':>8}{'modal%':>8}{'P(obs wins)':>13}"
    f"{'P(REM wins)':>13}   full distribution")
stab = []
for k in BOOT_KEYS:
    w = pd.Series(boot_win[k])
    if not len(w):
        continue
    row = core_df.set_index("key").loc[k]
    vc = w.value_counts(normalize=True)
    stab.append({"key": k, "label": row.label, "events": int(row.events),
                 "observed": row.winner, "modal": vc.index[0], "modal_frac": vc.iloc[0],
                 "p_obs": float(vc.get(row.winner, 0.0)), "p_rem": float(vc.get("rem", 0.0)),
                 **{f"p_{s}": float(vc.get(s, 0.0)) for s in sl}})
    say(f"{row.label[:21]:22s}{int(row.events):>6}{str(row.winner):>10}{vc.index[0]:>8}"
        f"{100 * vc.iloc[0]:>7.0f}%{100 * vc.get(row.winner, 0):>12.0f}%"
        f"{100 * vc.get('rem', 0):>12.0f}%   "
        + " ".join(f"{s}:{100 * vc.get(s, 0):.0f}" for s in sl))
stab_df = pd.DataFrame(stab)
if len(stab_df):
    stab_df.to_csv(f"{HERE}/attack_winner_bootstrap.csv", index=False)
    say(f"\nmedian P(the observed winner wins again) = {stab_df.p_obs.median():.2f}")
    say(f"outcomes where the observed winner is reproduced in >50% of bootstraps: "
        f"{int((stab_df.p_obs > 0.5).sum())} of {len(stab_df)}")
    say(f"outcomes where the observed winner is reproduced in >80% of bootstraps: "
        f"{int((stab_df.p_obs > 0.8).sum())} of {len(stab_df)}")
    say(f"median P(REM wins) = {stab_df.p_rem.median():.2f}; "
        f"REM is the modal winner in {int((stab_df.modal == 'rem').sum())} of {len(stab_df)}")
    say("expected P(winner) under 'all five equivalent' = 0.25 with four sleep windows")

# ---- likelihood-ratio chi-square, a scale-free comparator --------------------
say("\nSCALE CHECK. Per-1-SD hazard ratios are not comparable across stages when the tie mass")
say("at zero differs (N3 78.7%, REM 42.9%, sleep 29.1%). A scale-free comparator is the 1-df")
say("likelihood-ratio chi-square of adding the exposure to the covariate-only model.")
say(f"{'outcome':22s}{'ev':>6}" + "".join(f"{s:>9}" for s in sl) + f"{'best':>8}")
lr_win = []
for k in TEST_KEYS:
    d = frame_for(core, k, [f"z_{s}" for s in sl])
    base = fit(d, [])
    chis = {}
    for s in sl:
        r = fit(d, [f"z_{s}"])
        if r and base:
            chis[s] = 2 * (r["loglik"] - base["loglik"])
    if not chis:
        continue
    bw = max(chis, key=chis.get)
    lr_win.append({"key": k, "best": bw, **{f"chi_{s}": chis[s] for s in chis}})
    lab = core_df.set_index("key").loc[k, "label"]
    say(f"{lab[:21]:22s}{int(core_df.set_index('key').loc[k, 'events']):>6}"
        + "".join(f"{chis.get(s, float('nan')):>9.1f}" for s in sl) + f"{bw:>8}")
lrw = pd.DataFrame(lr_win)
say("LR-chi2 winner tally: " + ", ".join(f"{s} {int((lrw.best == s).sum())}" for s in sl))
_hw = core_df.set_index("key").loc[TEST_KEYS, "winner"]
say("HR winner tally on the same outcomes: "
    + ", ".join(f"{s} {int((_hw == s).sum())}" for s in sl))
_agree = float((lrw.set_index("key").best == _hw.loc[lrw.key].values).mean())
say(f"LR-chi2 winner agrees with the HR winner in {_agree:.0%} of these outcomes")

say("\nSECOND SCALE CHECK. Drop the rank scale entirely: any hypoxaemia in the window versus")
say("none (T90 > 0 vs T90 == 0). This unit is identical across stages and immune to the tie mass.")
for s in sl:
    core[f"bin_{s}"] = (core[f"t90_{s}_g30"] > 0).astype(float)
    say(f"  exposed fraction, {s}: {100 * core[f'bin_{s}'].mean():.1f}%")
say(f"{'outcome':22s}{'ev':>6}" + "".join(f"{s:>9}" for s in sl) + f"{'winner':>8}")
bin_rows = []
for k in TEST_KEYS:
    d = frame_for(core, k, [f"bin_{s}" for s in sl])
    hh = {}
    for s in sl:
        r = fit(d, [f"bin_{s}"])
        if r:
            hh[s] = r[f"bin_{s}"]["hr"]
    if not hh:
        continue
    bw = max(hh, key=hh.get)
    bin_rows.append({"key": k, "winner_bin": bw, **{f"{s}_hr": hh[s] for s in hh}})
    lab = core_df.set_index("key").loc[k, "label"]
    say(f"{lab[:21]:22s}{int(core_df.set_index('key').loc[k, 'events']):>6}"
        + "".join(f"{hh.get(s, float('nan')):>9.3f}" for s in sl) + f"{bw:>8}")
binw = pd.DataFrame(bin_rows)
binw.to_csv(f"{HERE}/attack_binary_scale.csv", index=False)
say("binary-contrast winner tally: "
    + ", ".join(f"{s} {int((binw.winner_bin == s).sum())}" for s in sl))
_ag2 = float((binw.set_index("key").winner_bin == _hw.loc[binw.key].values).mean())
say(f"binary winner agrees with the per-1-SD winner in {_ag2:.0%} of these outcomes")

# ================================================================================
head("ATTACK 2  DIFFERENTIAL GATING AND DIFFERENTIAL MEASUREMENT ERROR")
# ================================================================================
COV = ["AgeAtVisit", "male", "AHI", "t90_all", "spo2_mean", "tst_min", "sleep_efficiency_pct",
       "arousal_index", "n_prevalent_excl"]
say("who is lost to each >=30 min gate (mean in gated vs excluded, standardised difference)")
say(f"{'stage':7s}{'gated':>8}{'lost':>7}  " + "".join(f"{c[:11]:>13}" for c in COV))
gate_rows = []
for s in STAGES:
    g = b[f"t90_{s}_g30"].notna()
    rec = {"stage": s, "gated": int(g.sum()), "lost": int((~g).sum())}
    line = f"{s:7s}{int(g.sum()):>8,}{int((~g).sum()):>7,}  "
    for c in COV:
        x1, x0 = b.loc[g, c], b.loc[~g, c]
        sd = np.sqrt((x1.var() + x0.var()) / 2)
        smd = (x1.mean() - x0.mean()) / sd if sd > 0 else np.nan
        rec[f"{c}_gated"], rec[f"{c}_lost"], rec[f"{c}_smd"] = x1.mean(), x0.mean(), smd
        line += f"{smd:>13.2f}"
    # event rate for a few outcomes
    for k in ["death", "diabetes", "htn2", "aki"]:
        e1 = b.loc[g & (b[f"{k}_prevalent"] == 0), f"{k}_incident"].mean()
        e0 = b.loc[~g & (b[f"{k}_prevalent"] == 0), f"{k}_incident"].mean()
        rec[f"{k}_rate_gated"], rec[f"{k}_rate_lost"] = e1, e0
    gate_rows.append(rec)
    say(line + "   (SMD, gated minus lost)")
gates = pd.DataFrame(gate_rows)
gates.to_csv(f"{HERE}/attack_gate_selection.csv", index=False)
say("\nincident rates, gated vs lost (prevalent-free denominator)")
say(f"{'stage':7s}" + "".join(f"{k:>22}" for k in ["death", "diabetes", "htn2", "aki"]))
for _, r in gates.iterrows():
    say(f"{r['stage']:7s}" + "".join(
        f"{r[k + '_rate_gated']:>11.3f}{r[k + '_rate_lost']:>11.3f}"
        for k in ["death", "diabetes", "htn2", "aki"]))

say("\nDoes the REM gate select a SICKER group (which would manufacture the REM signal) or a")
say("HEALTHIER one (which would attenuate it)? positive SMD = gated group higher")
r = gates.set_index("stage").loc["rem"]
for c in COV:
    say(f"  {c:22s} gated {r[c + '_gated']:.2f}  lost {r[c + '_lost']:.2f}  SMD {r[c + '_smd']:+.3f}")

say("\nMEASUREMENT-TIME TEST. Within REM-gated patients, split by minutes of REM scored.")
say("If REM's advantage is real signal diluted by short sampling, the HR should be LARGER")
say("where REM was measured for longer. If it is a selection artefact, it should not track time.")
say(f"{'outcome':22s}{'tertile':>10}{'REM min':>10}{'n':>7}{'ev':>6}{'HR':>8}{'95% CI':>18}")
remg = b[b.t90_rem_g30.notna()].copy()
tert = pd.qcut(remg.min_rem, 3, labels=["short", "mid", "long"])
remg["tert"] = tert
timerows = []
for k in ["htn2", "obesity", "diabetes", "aki", "pneumonia"]:
    for lv in ["short", "mid", "long"]:
        sub = remg[remg.tert == lv].copy()
        sub["z_rem"] = rint_site(sub, "t90_rem_g30")
        d = frame_for(sub, k, ["z_rem"])
        r2 = fit(d, ["z_rem"])
        if r2:
            e = r2["z_rem"]
            timerows.append({"key": k, "tertile": lv, "rem_min": sub.min_rem.median(),
                             "n": r2["n"], "events": r2["events"], **e})
            ci = "%.2f-%.2f" % (e["lo"], e["hi"])
            say(f"{k[:21]:22s}{lv:>10}{sub.min_rem.median():>10.0f}{r2['n']:>7}{r2['events']:>6}"
                f"{e['hr']:>8.3f}{ci:>18}")
pd.DataFrame(timerows).to_csv(f"{HERE}/attack_rem_duration.csv", index=False)

# ================================================================================
head("ATTACK 3  RELIABILITY CEILING: NIGHT-TO-NIGHT ICC PER STAGE")
# ================================================================================
say("Computed from the raw pilot shards: patients with >=2 strictly PAP-off recordings")
say("(cpap_frac==0 with a PAP channel present) that both clear the same >=30 min stage gate.")
sh2 = sh[(sh.cpap_frac == 0) & (sh.has_cpap_chan.astype(str).str.lower().isin(["true", "1"]))
         & (sh.status == "ok")].copy()
# site_id is only populated on the spine rows in the shards; SiteID is complete
sh2["site_id"] = sh2.SiteID.fillna(sh2.site_id)
say(f"sites in the repeat-night set: {dict(sh2.site_id.value_counts())}")
say(f"strictly PAP-off shard rows: {len(sh2):,} over {sh2.BDSPPatientID.nunique():,} patients")
icc_rows = []
for s in STAGES:
    col = f"t90_{s}"
    mn = "rec_min" if s == "all" else f"min_{s}"
    v = sh2[sh2[col].notna()].copy()
    if s != "all":
        v = v[v[mn] >= 30]
    cnt = v.groupby("BDSPPatientID").size()
    rep = v[v.BDSPPatientID.isin(cnt[cnt >= 2].index)].copy()
    if rep.BDSPPatientID.nunique() < 30:
        icc_rows.append({"stage": s, "patients": rep.BDSPPatientID.nunique()})
        continue
    # two nights per patient (first two chronological by SessionID)
    rep = rep.sort_values(["BDSPPatientID", "SessionID"]).groupby("BDSPPatientID").head(2)
    # analysis scale: rank-inverse-normal within site across all measured nights
    rep["z"] = np.nan
    for _, g in rep.groupby("site_id"):
        rep.loc[g.index, "z"] = rint(g[col].values)
    w = rep.pivot_table(index="BDSPPatientID", columns=rep.groupby("BDSPPatientID").cumcount(),
                        values="z")
    w = w.dropna()
    if len(w) < 30:
        icc_rows.append({"stage": s, "patients": len(w)})
        say(f"  t90_{s:5s} only {len(w)} repeat pairs clear the gate on both nights, not estimable")
        continue
    x1, x2 = w[0].values, w[1].values
    # one-way random-effects ICC(1,1) from the two-night ANOVA
    m = (x1 + x2) / 2
    msb = 2 * np.var(m, ddof=1)
    msw = np.mean((x1 - x2) ** 2) / 2
    icc = (msb - msw) / (msb + msw)
    r_p = stats.pearsonr(x1, x2)[0]
    # raw-scale ICC too
    wr = rep.pivot_table(index="BDSPPatientID", columns=rep.groupby("BDSPPatientID").cumcount(),
                         values=col).dropna()
    y1, y2 = wr[0].values, wr[1].values
    mr = (y1 + y2) / 2
    iccr = (2 * np.var(mr, ddof=1) - np.mean((y1 - y2) ** 2) / 2) / \
           (2 * np.var(mr, ddof=1) + np.mean((y1 - y2) ** 2) / 2)
    boot = []
    for _ in range(2000):
        t2 = np.random.randint(0, len(x1), len(x1))
        a, c2 = x1[t2], x2[t2]
        mm = (a + c2) / 2
        bb = 2 * np.var(mm, ddof=1)
        ww = np.mean((a - c2) ** 2) / 2
        boot.append((bb - ww) / (bb + ww))
    lo, hi2 = np.percentile(boot, [2.5, 97.5])
    icc_rows.append({"stage": s, "patients": len(w), "icc_z": icc, "icc_lo": lo, "icc_hi": hi2,
                     "icc_raw": iccr, "pearson_z": r_p})
    say(f"  t90_{s:5s} n={len(w):>4} pairs   ICC(analysis scale)={icc:.3f} "
        f"({lo:.3f} to {hi2:.3f})   ICC(raw)={iccr:.3f}   r={r_p:.3f}")
icc_df = pd.DataFrame(icc_rows)
icc_df.to_csv(f"{HERE}/attack_icc.csv", index=False)

say("\nDISATTENUATION. With classical measurement error, an exposure standardised on its own")
say("OBSERVED scale has its log HR attenuated by sqrt(ICC), so the corrected log HR is")
say("logHR / sqrt(ICC). A NOISIER stage is therefore penalised MORE, and a stage that wins")
say("despite being noisier is understated, not flattered.")
iccm = icc_df.set_index("stage")
PRIOR_ICC = {"all": 0.42, "sleep": 0.40, "rem": 0.31}  # project's own repeat-night work
use_icc, icc_src = {}, {}
for s in sl:
    v = iccm["icc_z"].get(s, np.nan) if "icc_z" in iccm.columns else np.nan
    if v != v:
        v = PRIOR_ICC.get(s, np.nan)
        icc_src[s] = "prior"
    else:
        icc_src[s] = "recomputed"
    use_icc[s] = v
say("ICC source: " + ", ".join(f"{s}={icc_src[s]}" for s in sl))
say("ICC used: " + ", ".join(f"{s} {use_icc[s]:.3f}" for s in sl if use_icc[s] == use_icc[s]))
say(f"\n{'outcome':22s}" + "".join(f"{s + ' obs':>11}{s + ' corr':>11}" for s in ["sleep", "rem"])
    + f"{'REM/sleep obs':>15}{'REM/sleep corr':>16}")
dis_rows = []
for k in TEST_KEYS:
    row = core_df.set_index("key").loc[k]
    if not (row.get("rem_hr") == row.get("rem_hr")):
        continue
    cs, cr = np.log(row.sleep_hr), np.log(row.rem_hr)
    ds = cs / np.sqrt(use_icc["sleep"]) if use_icc["sleep"] == use_icc["sleep"] else np.nan
    dr = cr / np.sqrt(use_icc["rem"]) if use_icc["rem"] == use_icc["rem"] else np.nan
    dis_rows.append({"key": k, "label": row.label, "sleep_obs": row.sleep_hr,
                     "sleep_corr": np.exp(ds), "rem_obs": row.rem_hr, "rem_corr": np.exp(dr),
                     "ratio_obs": cr / cs, "ratio_corr": dr / ds})
    say(f"{row.label[:21]:22s}{row.sleep_hr:>11.3f}{np.exp(ds):>11.3f}"
        f"{row.rem_hr:>11.3f}{np.exp(dr):>11.3f}{cr / cs:>15.2f}{dr / ds:>16.2f}")
dis = pd.DataFrame(dis_rows)
dis.to_csv(f"{HERE}/attack_disattenuated.csv", index=False)
if len(dis):
    say(f"\nmedian observed REM/sleep log-HR ratio {dis.ratio_obs.median():.3f}; "
        f"after disattenuation {dis.ratio_corr.median():.3f}")
    exp_ratio = np.sqrt(use_icc["rem"] / use_icc["sleep"])
    say(f"if REM and whole-sleep measured the SAME underlying quantity, the observed REM/sleep")
    say(f"log-HR ratio would be sqrt(ICC_rem/ICC_sleep) = {exp_ratio:.3f}; observed is "
        f"{dis.ratio_obs.median():.3f}")

# ================================================================================
head("ATTACK 4  NEGATIVE CONTROLS RECOMPUTED PER STAGE")
# ================================================================================


def controls_for(frame, tag, stages):
    rows = []
    say(f"\n[{tag}] n={len(frame):,}")
    say(f"{'control':22s}{'ev':>6}" + "".join(f"{s:>10}" for s in stages))
    for k in NEG:
        d = frame_for(frame, k, [f"z_{s}" for s in stages])
        if d.E.sum() < 40:
            continue
        rec = {"control": k, "events": int(d.E.sum())}
        line = f"{k[:21]:22s}{int(d.E.sum()):>6}"
        for s in stages:
            r = fit(d, [f"z_{s}"])
            if r:
                rec[f"{s}_hr"] = r[f"z_{s}"]["hr"]
                rec[f"{s}_coef"] = r[f"z_{s}"]["coef"]
                rec[f"{s}_se"] = r[f"z_{s}"]["se"]
                rec[f"{s}_p"] = r[f"z_{s}"]["p"]
                line += f"{r[f'z_{s}']['hr']:>9.3f}" + ("*" if r[f"z_{s}"]["p"] < 0.05 else " ")
            else:
                line += f"{'.':>10}"
        rows.append(rec)
        say(line)
    cd = pd.DataFrame(rows)
    say(f"{'FLOOR (max HR)':22s}{'':>6}" + "".join(
        f"{cd[f'{s}_hr'].max():>10.3f}" for s in stages))
    pooled = {}
    for s in stages:
        wgt = 1 / cd[f"{s}_se"] ** 2
        cf = (cd[f"{s}_coef"] * wgt).sum() / wgt.sum()
        se = np.sqrt(1 / wgt.sum())
        pooled[s] = (np.exp(cf), np.exp(cf - 1.96 * se), np.exp(cf + 1.96 * se))
    say(f"{'POOLED control':22s}{'':>6}" + "".join(f"{pooled[s][0]:>10.3f}" for s in stages))
    say(f"{'  pooled upper 95%':22s}{'':>6}" + "".join(f"{pooled[s][2]:>10.3f}" for s in stages))
    return cd, pooled


ctl_core, pooled_core = controls_for(core, "CORE COMMON SAMPLE", CORE)
ctl_full, pooled_full = controls_for(b, "FULL COHORT, each stage on its own gated sample", STAGES)
ctl_core.to_csv(f"{HERE}/attack_controls_core.csv", index=False)
ctl_full.to_csv(f"{HERE}/attack_controls_full.csv", index=False)

FLOOR = {s: max(ctl_full[f"{s}_hr"].max(), pooled_full[s][2]) for s in STAGES}
say("\nWorking floor per stage = max(worst single control HR, upper 95% of the pooled control):")
say("  " + ", ".join(f"{s} {FLOOR[s]:.3f}" for s in STAGES))

say("\nSTRIKE LIST. The nine REM claims that were said to survive everything, tested against the")
say("REM floor recomputed here, on the CORE COMMON SAMPLE (so the comparison is within-person).")
NINE = ["aki", "gout", "diabetes", "pneumonia", "obesity", "cellulitis", "ckd", "htn2",
        "dyslipid"]
cf_idx = core_df.set_index("key")
say(f"{'claim':22s}{'ev':>6}{'REM HR':>9}{'95% CI':>18}{'sleep HR':>10}{'vs floor':>10}"
    f"{'beats sleep':>13}")
strike = []
for k in NINE:
    if k not in cf_idx.index:
        say(f"{k:22s}{'not fitted on the core common sample':>60}")
        continue
    r = cf_idx.loc[k]
    ok_floor = r.rem_hr > FLOOR["rem"]
    ok_beat = r.rem_hr > r.sleep_hr
    strike.append({"key": k, "events": int(r.events), "rem_hr": r.rem_hr, "rem_lo": r.rem_lo,
                   "sleep_hr": r.sleep_hr, "clears_floor": bool(ok_floor),
                   "lo_clears_floor": bool(r.rem_lo > FLOOR["rem"]), "beats_sleep": bool(ok_beat)})
    say(f"{r.label[:21]:22s}{int(r.events):>6}{r.rem_hr:>9.3f}"
        f"{f'{r.rem_lo:.2f}-{r.rem_hi:.2f}':>18}{r.sleep_hr:>10.3f}"
        f"{('PASS' if ok_floor else 'FAIL'):>10}{('yes' if ok_beat else 'NO'):>13}")
pd.DataFrame(strike).to_csv(f"{HERE}/attack_strike_list.csv", index=False)

say("\nIS THE CONFOUNDING ITSELF STAGE-SPECIFIC? Joint model on each negative control,")
say("whole-sleep T90 and REM T90 together, core common sample. If REM carries EXTRA signal on")
say("outcomes that cannot be caused by nocturnal hypoxaemia, its excess is a confounding")
say("channel, not a physiological one.")
say(f"{'control':22s}{'ev':>6}{'sleep HR|REM':>14}{'p':>9}{'REM HR|sleep':>14}{'p':>9}")
jrows = []
for k in NEG:
    d = frame_for(core, k, ["z_sleep", "z_rem"])
    r = fit(d, ["z_sleep", "z_rem"])
    if not r:
        continue
    jrows.append({"control": k, "events": r["events"],
                  "sleep_hr": r["z_sleep"]["hr"], "sleep_p": r["z_sleep"]["p"],
                  "rem_hr": r["z_rem"]["hr"], "rem_p": r["z_rem"]["p"]})
    say(f"{k[:21]:22s}{r['events']:>6}{r['z_sleep']['hr']:>14.3f}{r['z_sleep']['p']:>9.3f}"
        f"{r['z_rem']['hr']:>14.3f}{r['z_rem']['p']:>9.3f}")
jr = pd.DataFrame(jrows)
jr.to_csv(f"{HERE}/attack_controls_joint.csv", index=False)
say(f"REM exceeds whole-sleep on {int((jr.rem_hr > jr.sleep_hr).sum())} of {len(jr)} negative "
    f"controls; REM p<0.05 net of sleep on {int((jr.rem_p < 0.05).sum())}")

say("\nOBESITY ADJUSTMENT. BMI is not in this dataset. Prevalent obesity (an ICD proxy) is.")
say("The reverse-causation work found 26 of 30 outcomes moved on BMI adjustment, so the REM")
say("claims are refitted with prevalent obesity added as a covariate.")
core["prev_ob"] = core.obesity_prevalent.astype(int)
say(f"prevalent obesity in the core sample: {100 * core.prev_ob.mean():.1f}%")
say(f"{'claim':22s}{'REM crude':>11}{'REM +ob':>10}{'% change':>10}{'sleep crude':>13}"
    f"{'sleep +ob':>11}{'% change':>10}")
obrows = []
for k in [x for x in NINE if x != "obesity"] + NEG:
    d = frame_for(core, k, ["z_sleep", "z_rem"])
    d["prev_ob"] = core.loc[d.index, "prev_ob"].values
    out = {"key": k}
    for s in ["rem", "sleep"]:
        a = fit(d, [f"z_{s}"])
        bm = fit(d.assign(**{}), [f"z_{s}", "prev_ob"])
        if a and bm:
            out[f"{s}_crude"] = a[f"z_{s}"]["hr"]
            out[f"{s}_adj"] = bm[f"z_{s}"]["hr"]
            out[f"{s}_pct"] = 100 * (np.log(bm[f"z_{s}"]["hr"]) / np.log(a[f"z_{s}"]["hr"]) - 1)
    obrows.append(out)
    if "rem_crude" in out:
        say(f"{k[:21]:22s}{out['rem_crude']:>11.3f}{out['rem_adj']:>10.3f}{out['rem_pct']:>9.1f}%"
            f"{out['sleep_crude']:>13.3f}{out['sleep_adj']:>11.3f}{out['sleep_pct']:>9.1f}%")
ob = pd.DataFrame(obrows)
ob.to_csv(f"{HERE}/attack_obesity_adjusted.csv", index=False)
say(f"median log-HR change on obesity adjustment: REM {ob.rem_pct.median():+.1f}%, "
    f"whole-sleep {ob.sleep_pct.median():+.1f}%")

# ================================================================================
head("ATTACK 5  IS WHOLE-SLEEP T90 SUFFICIENT? NESTED MODELS AND CROSS-VALIDATED C")
# ================================================================================
say("Model A: whole-sleep T90 + age spline + sex, site strata.")
say("Model B: A + the stage. LR test on 1 df. Concordance is reported in-sample AND")
say("5-fold cross-validated, because an in-sample gain is guaranteed to be positive.")


def cv_delta_c(d, cols_a, cols_b, folds=5):
    dd = d[list(set(cols_a + cols_b)) + ADJ + ["T", "E", "site"]].dropna()
    if dd.E.sum() < 60:
        return np.nan, np.nan
    ix = np.random.permutation(len(dd))
    parts = np.array_split(ix, folds)
    ca, cb = [], []
    for f in range(folds):
        te = dd.iloc[parts[f]]
        tr = dd.iloc[np.concatenate([parts[j] for j in range(folds) if j != f])]
        if te.E.sum() < 5 or tr.E.sum() < 30:
            continue
        for cols, acc in ((cols_a, ca), (cols_b, cb)):
            try:
                m = CoxPHFitter(penalizer=0.0).fit(tr[cols + ADJ + ["T", "E", "site"]], "T", "E",
                                                   strata=["site"])
                pr = m.predict_partial_hazard(te)
                acc.append(concordance_index(te["T"], -pr, te["E"]))
            except Exception:
                acc.append(np.nan)
    return np.nanmean(ca), np.nanmean(cb)


say(f"\n{'outcome':22s}{'ev':>6}{'stage':>7}{'LR chi2':>9}{'p':>10}{'dC in':>8}{'dC cv':>9}"
    f"{'sleep p|stage':>15}{'stage p|sleep':>15}")
nest = []
for k in TEST_KEYS:
    d = frame_for(core, k, [f"z_{s}" for s in sl])
    a = fit(d, ["z_sleep"])
    if not a:
        continue
    for s in ["rem", "nrem", "n2", "wake"]:
        if s == "wake":
            d2 = frame_for(core, k, ["z_sleep", "z_wake"])
        else:
            d2 = d
        try:
            bmod = fit(d2, ["z_sleep", f"z_{s}"])
        except Exception:
            bmod = None
        if not bmod:
            continue
        chi = 2 * (bmod["loglik"] - a["loglik"])
        p = stats.chi2.sf(max(chi, 0), 1)
        ca, cb2 = cv_delta_c(d2, ["z_sleep"], ["z_sleep", f"z_{s}"])
        nest.append({"key": k, "label": cf_idx.loc[k, "label"], "events": a["events"],
                     "stage": s, "lr_chi2": chi, "lr_p": p,
                     "c_a": a["c"], "c_b": bmod["c"], "dc_in": bmod["c"] - a["c"],
                     "c_a_cv": ca, "c_b_cv": cb2, "dc_cv": cb2 - ca,
                     "sleep_p_adj": bmod["z_sleep"]["p"], "sleep_hr_adj": bmod["z_sleep"]["hr"],
                     "stage_p_adj": bmod[f"z_{s}"]["p"], "stage_hr_adj": bmod[f"z_{s}"]["hr"]})
        say(f"{cf_idx.loc[k, 'label'][:21]:22s}{a['events']:>6}{s:>7}{chi:>9.1f}{p:>10.4f}"
            f"{bmod['c'] - a['c']:>8.4f}{cb2 - ca:>9.4f}{bmod['z_sleep']['p']:>15.4f}"
            f"{bmod[f'z_{s}']['p']:>15.4f}")
nest_df = pd.DataFrame(nest)
nest_df.to_csv(f"{HERE}/attack_nested.csv", index=False)
for s in ["rem", "nrem", "n2", "wake"]:
    q = nest_df[nest_df.stage == s]
    if not len(q):
        continue
    say(f"\n{s.upper()} added to whole-sleep T90 across {len(q)} outcomes:")
    say(f"  LR p<0.05 in {int((q.lr_p < 0.05).sum())} of {len(q)}; "
        f"after Benjamini-Hochberg q<0.05 in "
        f"{int((stats.false_discovery_control(q.lr_p.values) < 0.05).sum())}")
    say(f"  median in-sample dC {q.dc_in.median():+.4f}; median 5-fold CV dC {q.dc_cv.median():+.4f}; "
        f"CV dC positive in {int((q.dc_cv > 0).sum())} of {len(q)}")
    say(f"  whole-sleep T90 keeps p<0.05 with the stage in the model: "
        f"{int((q.sleep_p_adj < 0.05).sum())} of {len(q)}; "
        f"stage keeps p<0.05: {int((q.stage_p_adj < 0.05).sum())} of {len(q)}")
    say(f"  median absolute CV concordance: whole-sleep alone {q.c_a_cv.median():.4f}, "
        f"with the stage {q.c_b_cv.median():.4f}")

# reverse direction: does whole-sleep add to REM
say("\nREVERSE: does whole-sleep T90 add anything to a model that already has REM?")
rev = []
for k in TEST_KEYS:
    d = frame_for(core, k, ["z_sleep", "z_rem"])
    a = fit(d, ["z_rem"])
    bmod = fit(d, ["z_rem", "z_sleep"])
    if not (a and bmod):
        continue
    chi = 2 * (bmod["loglik"] - a["loglik"])
    ca, cb2 = cv_delta_c(d, ["z_rem"], ["z_rem", "z_sleep"])
    rev.append({"key": k, "lr_chi2": chi, "lr_p": stats.chi2.sf(max(chi, 0), 1),
                "dc_cv": cb2 - ca})
rev = pd.DataFrame(rev)
say(f"  whole-sleep added to REM: LR p<0.05 in {int((rev.lr_p < 0.05).sum())} of {len(rev)}; "
    f"median CV dC {rev.dc_cv.median():+.4f}")
say(f"  REM added to whole-sleep: LR p<0.05 in "
    f"{int((nest_df[nest_df.stage == 'rem'].lr_p < 0.05).sum())} of "
    f"{len(nest_df[nest_df.stage == 'rem'])}; median CV dC "
    f"{nest_df[nest_df.stage == 'rem'].dc_cv.median():+.4f}")

# ================================================================================
head("HEALTHY SUBGROUP: the same attacks where power allows")
# ================================================================================
hcore = core[core.healthy].copy()
for s in CORE:
    hcore[f"z_{s}"] = rint_site(hcore, f"t90_{s}_g30")
say(f"healthy in the full cohort {int(b.healthy.sum()):,}; healthy in the core common sample "
    f"{len(hcore):,}")
say("\nzero mass per stage in the healthy core sample")
for s in CORE:
    v = hcore[f"t90_{s}_g30"].dropna()
    say(f"  t90_{s:6s} n={len(v):,}  {100 * (v == 0).mean():.1f}% exactly zero  "
        f"median {v.median():.3f}")
hrows = []
say(f"\n{'outcome':22s}{'ev':>6}" + "".join(f"{s:>10}" for s in CORE) + f"{'winner':>10}")
for k, lab, neg in OUTCOMES:
    d = frame_for(hcore, k, [f"z_{s}" for s in CORE])
    if d.E.sum() < 40:
        continue
    rec = {"key": k, "label": lab, "neg": neg, "events": int(d.E.sum()), "n": len(d)}
    for s in CORE:
        r = fit(d, [f"z_{s}"])
        if r:
            rec[f"{s}_hr"] = r[f"z_{s}"]["hr"]
            rec[f"{s}_p"] = r[f"z_{s}"]["p"]
            rec[f"{s}_lo"] = r[f"z_{s}"]["lo"]
            rec[f"{s}_hi"] = r[f"z_{s}"]["hi"]
            rec[f"{s}_coef"] = r[f"z_{s}"]["coef"]
            rec[f"{s}_se"] = r[f"z_{s}"]["se"]
    hh = {s: rec.get(f"{s}_hr") for s in sl if rec.get(f"{s}_hr") == rec.get(f"{s}_hr")}
    rec["winner"] = max(hh, key=hh.get) if hh else None
    hrows.append(rec)
    say(f"{lab[:21]:22s}{rec['events']:>6}"
        + "".join(f"{rec.get(s + '_hr', float('nan')):>10.3f}" for s in CORE)
        + f"{str(rec['winner']):>10}" + ("   [negative control]" if neg else ""))
hdf = pd.DataFrame(hrows)
hdf.to_csv(f"{HERE}/attack_healthy_panel.csv", index=False)
say("healthy winner tally: " + ", ".join(
    f"{s} {int((hdf[~hdf.neg].winner == s).sum())}" for s in sl))
say(f"median 95% CI width, healthy REM {np.nanmedian(hdf.rem_hi - hdf.rem_lo):.3f} vs "
    f"general core REM {med_ciw['rem']:.3f}")

HKEYS = [r["key"] for r in hrows if not r["neg"] and r["key"] not in CIRCULAR
         and r["events"] >= 80]
HKEYS = sorted(HKEYS, key=lambda k: -hdf.set_index("key").loc[k, "events"])[:6]
say(f"\nbootstrapping the healthy winner for {HKEYS}, B={B_BOOT}")
hb = {k: [] for k in HKEYS}
for bi in range(0 if not SKIP_BOOT else B_BOOT, B_BOOT):
    take = np.random.randint(0, len(hcore), len(hcore))
    cb = hcore.iloc[take].copy()
    cb.index = np.arange(len(cb))
    for s in sl:
        cb[f"z_{s}"] = rint_site(cb, f"t90_{s}_g30")
    for k in HKEYS:
        d = frame_for(cb, k, [f"z_{s}" for s in sl])
        hh = {}
        for s in sl:
            try:
                r = fit(d, [f"z_{s}"])
            except Exception:
                r = None
            if r:
                hh[s] = r[f"z_{s}"]["hr"]
        if len(hh) == len(sl):
            hb[k].append(max(hh, key=hh.get))
hstab = []
say(f"{'outcome':22s}{'ev':>6}{'observed':>10}{'P(obs)':>9}{'P(REM)':>9}")
for k in HKEYS:
    w = pd.Series(hb[k])
    if not len(w):
        continue
    vc = w.value_counts(normalize=True)
    row = hdf.set_index("key").loc[k]
    hstab.append({"key": k, "observed": row.winner, "p_obs": float(vc.get(row.winner, 0)),
                  "p_rem": float(vc.get("rem", 0)), "events": int(row.events)})
    say(f"{row.label[:21]:22s}{int(row.events):>6}{str(row.winner):>10}"
        f"{100 * vc.get(row.winner, 0):>8.0f}%{100 * vc.get('rem', 0):>8.0f}%")
_hs = pd.DataFrame(hstab)
if len(_hs):  # never clobber a real bootstrap file on a SKIP_BOOT rerun
    _hs.to_csv(f"{HERE}/attack_healthy_bootstrap.csv", index=False)

say("\nhealthy negative controls, core common sample")
try:
    hctl, hpooled = controls_for(hcore, "HEALTHY CORE", CORE)
    hctl.to_csv(f"{HERE}/attack_controls_healthy.csv", index=False)
except Exception as e:
    say(f"  healthy controls not estimable: {e}")

say("\nhealthy nested test: does REM add to whole-sleep T90 in the disease-free subgroup?")
say(f"{'outcome':22s}{'ev':>6}{'LR chi2':>9}{'p':>9}{'dC cv':>9}{'sleep p|REM':>13}{'REM p|sleep':>13}")
hn = []
for k in HKEYS:
    d = frame_for(hcore, k, ["z_sleep", "z_rem"])
    a = fit(d, ["z_sleep"])
    bmod = fit(d, ["z_sleep", "z_rem"])
    if not (a and bmod):
        continue
    chi = 2 * (bmod["loglik"] - a["loglik"])
    ca, cb2 = cv_delta_c(d, ["z_sleep"], ["z_sleep", "z_rem"])
    hn.append({"key": k, "lr_chi2": chi, "lr_p": stats.chi2.sf(max(chi, 0), 1), "dc_cv": cb2 - ca,
               "sleep_p": bmod["z_sleep"]["p"], "rem_p": bmod["z_rem"]["p"],
               "sleep_hr": bmod["z_sleep"]["hr"], "rem_hr": bmod["z_rem"]["hr"]})
    say(f"{hdf.set_index('key').loc[k, 'label'][:21]:22s}{a['events']:>6}{chi:>9.1f}"
        f"{stats.chi2.sf(max(chi, 0), 1):>9.4f}{cb2 - ca:>9.4f}{bmod['z_sleep']['p']:>13.4f}"
        f"{bmod['z_rem']['p']:>13.4f}")
hn = pd.DataFrame(hn)
hn.to_csv(f"{HERE}/attack_healthy_nested.csv", index=False)
if len(hn):
    say(f"  REM adds at p<0.05 in {int((hn.lr_p < 0.05).sum())} of {len(hn)}; "
        f"median CV dC {hn.dc_cv.median():+.4f}")

# ================================================================================
head("RECONCILIATION against the delivered tables")
# ================================================================================
try:
    g = pd.read_csv(f"{HERE}/stage_general.csv")
    say(f"stage_general.csv {g.shape}")
    kk = [c for c in g.columns if c.endswith("_hr")][:12]
    say("columns: " + ", ".join(list(g.columns)[:14]))
    j = g.set_index("key") if "key" in g.columns else None
    if j is not None:
        for s in ["sleep", "rem"]:
            col = f"t90_{s}_hr"
            if col in j.columns:
                m2 = core_df.set_index("key")[[f"{s}_hr"]].join(j[[col]], how="inner").dropna()
                rho = stats.spearmanr(m2[f"{s}_hr"], m2[col])[0]
                say(f"  {s}: n={len(m2)} matched, Spearman {rho:.3f}, "
                    f"median delivered minus mine {np.median(m2[col] - m2[f'{s}_hr']):+.3f}")
except Exception as e:
    say(f"reconciliation skipped: {e}")

say(f"\nTOTAL RUNTIME {time.time() - T0:.0f}s")
json.dump({"floor": FLOOR, "icc": icc_df.to_dict("records"),
           "median_p_obs_winner": (float(stab_df.p_obs.median()) if len(stab_df) else None),
           "core_n": int(len(core)), "healthy_core_n": int(len(hcore))},
          open(f"{HERE}/attack_summary{'_noboot' if SKIP_BOOT else ''}.json", "w"),
          indent=1, default=str)
LOG.close()
