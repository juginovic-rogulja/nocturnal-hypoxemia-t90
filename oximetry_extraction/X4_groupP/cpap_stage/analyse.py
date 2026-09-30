#!/usr/bin/env python3
"""
Does CPAP fix REM hypoxaemia as well as it fixes the other stages?

DESIGN A  same night, PAP-off portion vs PAP-on portion (split nights)
DESIGN B  diagnostic night vs later titration night, same patient
CONTROL   PAP-off nights cut at a matched split time (the within-night time trend)

Every T90 is the project-standard  mean(spo2[50..100] < 90) * 100.
A stage is used for a patient only when it has >= MINSTAGE minutes on BOTH sides.
"""
from __future__ import annotations
import glob, os, sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).parent
MINSTAGE = 20.0            # minutes required on EACH side
STAGES = ["wake", "sleep", "n1", "n2", "n3", "rem"]
OUT = Path(os.environ.get("CPAP_STAGE_OUT") or "/tmp/cpap_t90_by_stage.parquet")   # v8.2: the one change against the pilot's analyse.py, rebuild_cpap_stage.py passes a path beside the shards
pd.set_option("display.width", 250)


def load(pat):
    fs = sorted(glob.glob(str(HERE / pat)))
    if not fs:
        return pd.DataFrame()
    return pd.concat([pd.read_csv(f, low_memory=False) for f in fs], ignore_index=True)


def paired(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    return np.asarray(a)[m], np.asarray(b)[m]


def block(df, pre, post, label, minutes=True):
    """Per-stage before/after table."""
    rows = []
    for s in STAGES:
        cb, ca = f"{pre}_{s}_t90", f"{post}_{s}_t90"
        mb, ma = f"{pre}_{s}_min", f"{post}_{s}_min"
        if cb not in df.columns:
            continue
        g = df[(df[mb] >= MINSTAGE) & (df[ma] >= MINSTAGE)]
        a, b = paired(g[cb].values, g[ca].values)
        if a.size < 8:
            rows.append(dict(stage=s, n=a.size))
            continue
        w = stats.wilcoxon(a, b)
        rel = (np.median(a) - np.median(b)) / np.median(a) * 100 if np.median(a) > 0 else np.nan
        # patient-level relative reduction, only where there was something to fix
        sub = a >= 1.0
        relpat = np.median((a[sub] - b[sub]) / a[sub] * 100) if sub.sum() >= 8 else np.nan
        rows.append(dict(
            stage=s, n=a.size,
            drop_lt20=int(((df[mb] < MINSTAGE) | (df[ma] < MINSTAGE)).sum()),
            min_before=float(np.median(df[mb][df[mb] > 0])),
            min_after=float(np.median(df[ma][df[ma] > 0])),
            t90_before=float(np.median(a)), t90_after=float(np.median(b)),
            iqr_before=f"{np.percentile(a,25):.2f}-{np.percentile(a,75):.2f}",
            iqr_after=f"{np.percentile(b,25):.2f}-{np.percentile(b,75):.2f}",
            abs_drop=float(np.median(a - b)),
            rel_drop_med=rel, rel_drop_pat=relpat,
            pct_improved=float(100 * np.mean(b < a)),
            p=w.pvalue))
    t = pd.DataFrame(rows)
    print(f"\n{'='*100}\n{label}\n{'='*100}")
    print(t.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    return t


def complete_case(df, pre, post, label, stages=("n2", "rem"), thr=MINSTAGE):
    """Restrict to patients who clear the minute gate in EVERY listed stage on BOTH sides,
    so the stage comparison is within the same people."""
    m = np.ones(len(df), bool)
    for s in stages:
        m &= (df[f"{pre}_{s}_min"] >= thr) & (df[f"{post}_{s}_min"] >= thr)
        m &= np.isfinite(df[f"{pre}_{s}_t90"]) & np.isfinite(df[f"{post}_{s}_t90"])
    g = df[m]
    print(f"\n{label}  complete cases with >={thr:.0f} min in {list(stages)} on BOTH sides: "
          f"{len(g)} of {len(df)}")
    if len(g) < 15:
        return pd.DataFrame(), g
    rows = []
    for s in stages:
        a = g[f"{pre}_{s}_t90"].values; b = g[f"{post}_{s}_t90"].values
        w = stats.wilcoxon(a, b)
        sub = a >= 1.0
        rows.append(dict(stage=s, n=len(g),
                         t90_before=float(np.median(a)), t90_after=float(np.median(b)),
                         abs_drop=float(np.median(a - b)),
                         rel_drop_med=float((np.median(a) - np.median(b)) / np.median(a) * 100)
                         if np.median(a) > 0 else np.nan,
                         rel_drop_pat=float(np.median((a[sub] - b[sub]) / a[sub] * 100))
                         if sub.sum() >= 8 else np.nan,
                         pct_improved=float(100 * np.mean(b < a)), p=w.pvalue))
    t = pd.DataFrame(rows)
    print(t.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    # direct paired contrast of the log-reduction, REM vs N2, in the same patient
    if set(stages) >= {"n2", "rem"}:
        dl = {}
        for s in ("n2", "rem"):
            dl[s] = (np.log1p(g[f"{pre}_{s}_t90"].values)
                     - np.log1p(g[f"{post}_{s}_t90"].values))
        w = stats.wilcoxon(dl["rem"], dl["n2"])
        print(f"  paired log1p reduction REM {np.median(dl['rem']):.3f} vs "
              f"N2 {np.median(dl['n2']):.3f}  diff {np.median(dl['rem']-dl['n2']):+.3f}  "
              f"p={w.pvalue:.3g}  (positive = REM improved MORE)")
    return t, g


def gate_sensitivity(df, pre, post, label):
    print(f"\n-- gate sensitivity ({label}): n and median relative reduction per stage --")
    rows = []
    for thr in (10, 20, 30, 45):
        r = {"gate_min": thr}
        for s in STAGES:
            mb, ma = f"{pre}_{s}_min", f"{post}_{s}_min"
            g = df[(df[mb] >= thr) & (df[ma] >= thr)]
            a, b = paired(g[f"{pre}_{s}_t90"].values, g[f"{post}_{s}_t90"].values)
            sub = a >= 1.0
            r[f"{s}_n"] = a.size
            r[f"{s}_rel"] = float(np.median((a[sub] - b[sub]) / a[sub] * 100)) if sub.sum() >= 8 else np.nan
        rows.append(r)
    print(pd.DataFrame(rows).to_string(index=False, float_format=lambda v: f"{v:.1f}"))


def ancova(df, pre, post, label, ref="n2"):
    """FLOOR-SAFE test.  Many stages hit T90 = 0 on treatment, so 'bigger log drop' is
    partly mechanical: REM starts higher, so it has more room to fall.  The question that
    survives a floor is: AT THE SAME STARTING HYPOXAEMIA, does one stage end up with more
    residual hypoxaemia than another?  Positive stage coefficient = fixed LESS well."""
    import statsmodels.formula.api as smf
    rows = []
    for s in ["n1", "n2", "n3", "rem"]:
        mb, ma = f"{pre}_{s}_min", f"{post}_{s}_min"
        g = df[(df[mb] >= MINSTAGE) & (df[ma] >= MINSTAGE)
               & np.isfinite(df[f"{pre}_{s}_t90"]) & np.isfinite(df[f"{post}_{s}_t90"])]
        rows.append(pd.DataFrame(dict(pid=g.BDSPPatientID.values, stage=s,
                                      y0=np.log1p(g[f"{pre}_{s}_t90"].values),
                                      y1=np.log1p(g[f"{post}_{s}_t90"].values))))
    d = pd.concat(rows, ignore_index=True)
    d["stage"] = pd.Categorical(d.stage, categories=[ref] + [x for x in ["n1", "n2", "n3", "rem"] if x != ref])
    m = smf.mixedlm("y1 ~ y0 + C(stage)", d, groups=d.pid).fit(reml=False)
    print(f"\n--- ANCOVA {label}: residual log1p(T90) on treatment, adjusted for the SAME "
          f"stage's pre-treatment value.  ref = {ref}.  positive coef = fixed LESS well ---")
    print(pd.DataFrame({"coef": m.params, "se": m.bse, "p": m.pvalues})
          .to_string(float_format=lambda v: f"{v:.4f}"))
    # same thing restricted to a matched starting range so it is not an extrapolation
    for lo, hi in [(1.0, 100.0), (5.0, 100.0)]:
        dd = d[(np.expm1(d.y0) >= lo) & (np.expm1(d.y0) <= hi)]
        if dd.stage.nunique() < 2 or len(dd) < 60:
            continue
        mm = smf.mixedlm("y1 ~ y0 + C(stage)", dd, groups=dd.pid).fit(reml=False)
        c = [i for i in mm.params.index if "rem" in i]
        if c:
            print(f"    restricted to pre-treatment T90 >= {lo:.0f}%  (n={len(dd)}): "
                  f"REM vs {ref} coef {mm.params[c[0]]:+.4f}  p={mm.pvalues[c[0]]:.3g}")
    return m


def severity_matched_placebo(P, A):
    """The untreated nights are far milder than the treated ones, so the raw placebo change
    is pinned at zero and cannot correct anything.  Restrict the control to nights whose
    FIRST-part T90 in that stage is at least as bad as the treated group's lower quartile."""
    print(f"\n{'='*100}\nSEVERITY-MATCHED CONTROL — untreated nights restricted to the treated "
          f"group's own baseline range, per stage\n{'='*100}")
    rows = []
    for s in STAGES:
        ta = A[(A[f"pre_{s}_min"] >= MINSTAGE) & (A[f"post_{s}_min"] >= MINSTAGE)]
        a0 = ta[f"pre_{s}_t90"].dropna().values
        if a0.size < 20:
            continue
        lo = float(np.percentile(a0, 25))
        g = P[(P[f"h1_{s}_min"] >= MINSTAGE) & (P[f"h2_{s}_min"] >= MINSTAGE)
              & (P[f"h1_{s}_t90"] >= lo) & np.isfinite(P[f"h2_{s}_t90"])]
        if len(g) < 15:
            rows.append(dict(stage=s, thr=lo, n_ctrl=len(g))); continue
        b0, b1 = g[f"h1_{s}_t90"].values, g[f"h2_{s}_t90"].values
        w = stats.wilcoxon(b0, b1)
        at = ta[ta[f"pre_{s}_t90"] >= lo]
        a1, a2 = at[f"pre_{s}_t90"].values, at[f"post_{s}_t90"].values
        rows.append(dict(stage=s, thr=lo, n_ctrl=len(g), n_tx=len(at),
                         ctrl_before=float(np.median(b0)), ctrl_after=float(np.median(b1)),
                         ctrl_rel=float(np.median((b0 - b1) / b0 * 100)),
                         ctrl_p=w.pvalue,
                         tx_before=float(np.median(a1)), tx_after=float(np.median(a2)),
                         tx_rel=float(np.median((a1 - a2) / a1 * 100)),
                         net_rel=float(np.median((a1 - a2) / a1 * 100)
                                       - np.median((b0 - b1) / b0 * 100))))
    print(pd.DataFrame(rows).to_string(index=False, float_format=lambda v: f"{v:.3f}"))


def long_form(df, pre, post, key, arm):
    """Long table for the mixed model: one row per patient x stage x period."""
    out = []
    for s in STAGES:
        cb, ca = f"{pre}_{s}_t90", f"{post}_{s}_t90"
        mb, ma = f"{pre}_{s}_min", f"{post}_{s}_min"
        if cb not in df.columns:
            continue
        g = df[(df[mb] >= MINSTAGE) & (df[ma] >= MINSTAGE)
               & np.isfinite(df[cb]) & np.isfinite(df[ca])]
        for per, col in (("before", cb), ("on", ca)):
            out.append(pd.DataFrame(dict(pid=g[key].values, arm=arm, stage=s,
                                         period=per, t90=g[col].values)))
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def mixed(L, label, ref="n2"):
    """log1p(T90) ~ period * stage, patient random effect."""
    try:
        import statsmodels.formula.api as smf
    except Exception as e:
        print("statsmodels missing:", e); return
    d = L.copy()
    d["y"] = np.log1p(d.t90)
    d["on"] = (d.period == "on").astype(float)
    d = d[d.stage.isin(["n1", "n2", "n3", "rem"])]
    d["stage"] = pd.Categorical(d.stage, categories=[ref] + [s for s in ["n1", "n2", "n3", "rem"] if s != ref])
    m = smf.mixedlm("y ~ on * C(stage)", d, groups=d.pid).fit(reml=False)
    print(f"\n--- MIXED MODEL {label}: log1p(T90) ~ PAP * stage, ref stage = {ref}, "
          f"patient random effect (n obs={len(d)}, patients={d.pid.nunique()}) ---")
    tab = pd.DataFrame({"coef": m.params, "se": m.bse, "z": m.tvalues, "p": m.pvalues})
    print(tab.to_string(float_format=lambda v: f"{v:.4f}"))
    # joint test of the interaction
    inter = [c for c in m.params.index if c.startswith("on:")]
    if inter:
        R = np.zeros((len(inter), len(m.params)))
        for i, c in enumerate(inter):
            R[i, list(m.params.index).index(c)] = 1
        w = m.wald_test(R, scalar=True)
        print(f"  JOINT stage-by-PAP interaction: chi2={float(w.statistic):.2f} "
              f"df={len(inter)} p={float(w.pvalue):.3g}")
    return m


# ====================================================================== DESIGN A
A = load("outA_*.csv")
if len(A):
    A = A[A.status == "ok"].copy()
    print(f"DESIGN A raw sessions read: {len(A)}   sites {A.SiteID.value_counts().to_dict()}")
    A["pre_min_tot"] = A["pre_any_min"]
    A["post_min_tot"] = A["post_any_min"]
    keep = (np.isfinite(A.cpap_start_min) & (A.pre_any_min >= 45) & (A.post_any_min >= 45)
            & (A.pap_frac_pre < 0.10) & (A.pap_frac_post >= 0.80))
    print("  usable split nights (switch found, >=45 min each side, "
          "PAP off before <10%% / on after >=80%%): %d" % keep.sum())
    print("  reasons excluded: no switch %d | side <45 min %d | PAP leaked pre %d | "
          "PAP not sustained post %d"
          % ((~np.isfinite(A.cpap_start_min)).sum(),
             ((A.pre_any_min < 45) | (A.post_any_min < 45)).sum(),
             (A.pap_frac_pre >= 0.10).sum(), (A.pap_frac_post < 0.80).sum()))
    A = A[keep].copy()
    print("  switch-on minute: median %.0f  IQR %.0f-%.0f"
          % (A.cpap_start_min.median(), A.cpap_start_min.quantile(.25),
             A.cpap_start_min.quantile(.75)))
    tA = block(A, "pre", "post", "DESIGN A — SAME NIGHT, PAP-off portion vs PAP-on portion")
    print("\nREM minutes before vs after the switch (treated nights): "
          "median %.1f vs %.1f  (ratio %.2f)"
          % (A.pre_rem_min.median(), A.post_rem_min.median(),
             A.post_rem_min.median() / max(A.pre_rem_min.median(), 1e-9)))
    print("REM share of staged time: before %.1f%%  after %.1f%%"
          % (100 * (A.pre_rem_min / A.pre_sleep_min.replace(0, np.nan)).median(),
             100 * (A.post_rem_min / A.post_sleep_min.replace(0, np.nan)).median()))
    gate_sensitivity(A, "pre", "post", "Design A")
    tAcc, Acc = complete_case(A, "pre", "post", "DESIGN A")
    LA = long_form(A, "pre", "post", "BDSPPatientID", "A")
    mixed(LA, "DESIGN A")
    if len(Acc) > 30:
        mixed(long_form(Acc, "pre", "post", "BDSPPatientID", "Acc"),
              "DESIGN A, complete cases only")
    ancova(A, "pre", "post", "DESIGN A")
else:
    A, tA, LA = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

# ====================================================================== CONTROL
P = load("outP_*.csv")
if len(P):
    P = P[(P.status == "ok")].copy()
    kp = ((P.cpap_frac.fillna(0) < 0.02) & (P.h1_any_min >= 45) & (P.h2_any_min >= 45))
    print(f"\nPLACEBO nights read {len(P)}, usable {int(kp.sum())}")
    P = P[kp].copy()
    print("  split minute used: median %.0f  IQR %.0f-%.0f"
          % (P.split_min_used.median(), P.split_min_used.quantile(.25),
             P.split_min_used.quantile(.75)))
    tP = block(P, "h1", "h2", "CONTROL — PAP NEVER ON, cut at a matched split time")
    print("\nREM minutes first vs second part (placebo nights): median %.1f vs %.1f"
          % (P.h1_rem_min.median(), P.h2_rem_min.median()))
    gate_sensitivity(P, "h1", "h2", "placebo")
    complete_case(P, "h1", "h2", "PLACEBO")
    LP = long_form(P, "h1", "h2", "BDSPPatientID", "P")
    mixed(LP, "PLACEBO (time trend only)")
    a, b = paired(P.h1_any_t90.values, P.h2_any_t90.values)
    w = stats.wilcoxon(a, b)
    print(f"\nPLACEBO, WHOLE-SEGMENT T90 (no stage stratification), n={a.size}: "
          f"first {np.median(a):.3f}% -> second {np.median(b):.3f}%, "
          f"paired median diff {np.median(b-a):+.3f}%, p={w.pvalue:.3g}, "
          f"fell in {100*np.mean(b<a):.1f}%   [this is the back-loading confound]")
    if len(A):
        severity_matched_placebo(P, A)
else:
    P, tP, LP = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

# ============================================== PLACEBO-CORRECTED DESIGN A
if len(tA) and len(tP):
    print(f"\n{'='*100}\nPLACEBO-CORRECTED DESIGN A  (treated change minus untreated time trend)\n{'='*100}")
    j = tA.merge(tP, on="stage", suffixes=("_tx", "_pl"))
    j["abs_drop_corrected"] = j.abs_drop_tx - j.abs_drop_pl
    j["rel_corrected"] = j.abs_drop_corrected / j.t90_before_tx * 100
    print(j[["stage", "n_tx", "n_pl", "t90_before_tx", "t90_after_tx", "abs_drop_tx",
             "abs_drop_pl", "abs_drop_corrected", "rel_drop_med_tx", "rel_corrected"]]
          .to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    if len(LA) and len(LP):
        LL = pd.concat([LA, LP], ignore_index=True)
        LL["pid"] = LL.arm + "_" + LL.pid.astype(str)
        LL["y"] = np.log1p(LL.t90)
        LL["on"] = (LL.period == "on").astype(float)
        LL["tx"] = (LL.arm == "A").astype(float)
        LL = LL[LL.stage.isin(["n1", "n2", "n3", "rem"])]
        LL["stage"] = pd.Categorical(LL.stage, categories=["n2", "n1", "n3", "rem"])
        try:
            import statsmodels.formula.api as smf
            m = smf.mixedlm("y ~ on * C(stage) * tx", LL, groups=LL.pid).fit(reml=False)
            print("\n--- TRIPLE-INTERACTION MODEL (treated vs placebo, stage, period) ---")
            k = [i for i in m.params.index if ("on:" in i and "tx" in i)]
            tt = pd.DataFrame({"coef": m.params, "se": m.bse, "p": m.pvalues})
            print(tt.loc[[i for i in m.params.index if "on" in i]]
                  .to_string(float_format=lambda v: f"{v:.4f}"))
        except Exception as e:
            print("triple model failed:", e)

# ====================================================================== DESIGN B
B = load("outB_*.csv")
tB = pd.DataFrame()
if len(B):
    B = B[B.status == "ok"].copy()
    B["is_dx"] = B.cpap_frac.fillna(0) < 0.02
    B["is_tx"] = B.cpap_frac >= 0.80
    pairs = []
    for pid, g in B.sort_values("rank").groupby("BDSPPatientID"):
        dg, tx = g[g.is_dx], g[g.is_tx]
        if len(dg) and len(tx):
            d0 = dg.iloc[0]
            later = tx[tx["rank"] > d0["rank"]]
            if len(later):
                t0 = later.iloc[0]
                rec = dict(BDSPPatientID=pid, SiteID=d0.SiteID, gap=d0.pair_gap)
                # v8.2 (2026-09-15): carry the collapse flag of each night into the pair (the stage values are already masked)
                rec["_v8_stage_collapsed_dx"] = int(d0.get("_v8_stage_collapsed", 0) or 0)
                rec["_v8_stage_collapsed_tx"] = int(t0.get("_v8_stage_collapsed", 0) or 0)
                rec["_v8_stage_collapsed"] = max(rec["_v8_stage_collapsed_dx"], rec["_v8_stage_collapsed_tx"])
                for s in STAGES + ["any"]:
                    for f in ("t90", "min"):
                        rec[f"dx_{s}_{f}"] = d0.get(f"all_{s}_{f}", np.nan)
                        rec[f"tx_{s}_{f}"] = t0.get(f"all_{s}_{f}", np.nan)
                pairs.append(rec)
    BP = pd.DataFrame(pairs)
    print(f"\n\nDESIGN B sessions read {len(B)}, patients {B.BDSPPatientID.nunique()}, "
          f"eligible diagnostic->titration pairs {len(BP)}")
    if len(BP):
        print("  median gap between nights: %.0f days" % BP.gap.median())
        tB = block(BP, "dx", "tx", "DESIGN B — DIAGNOSTIC NIGHT vs LATER TITRATION NIGHT")
        gate_sensitivity(BP, "dx", "tx", "Design B")
        tBcc, Bcc = complete_case(BP, "dx", "tx", "DESIGN B")
        LB = long_form(BP, "dx", "tx", "BDSPPatientID", "B")
        mixed(LB, "DESIGN B")
        if len(Bcc) > 30:
            mixed(long_form(Bcc, "dx", "tx", "BDSPPatientID", "Bcc"),
                  "DESIGN B, complete cases only")
        ancova(BP, "dx", "tx", "DESIGN B")
else:
    BP = pd.DataFrame()

# ====================================================================== RESIDUAL
print(f"\n{'='*100}\nRESIDUAL HYPOXAEMIA ON TREATMENT — which stage is still worst?\n{'='*100}")
for nm, df, tag in (("Design A on-PAP portion", A, "post"),
                    ("Design B titration night", BP, "tx")):
    if not len(df):
        continue
    r = []
    for s in STAGES:
        c, m = f"{tag}_{s}_t90", f"{tag}_{s}_min"
        if c not in df.columns:
            continue
        g = df[(df[m] >= MINSTAGE) & np.isfinite(df[c])]
        if len(g) < 8:
            continue
        r.append(dict(stage=s, n=len(g), median=float(np.median(g[c])),
                      p75=float(np.percentile(g[c], 75)),
                      mean=float(np.mean(g[c])),
                      pct_over_1=float(100 * np.mean(g[c] >= 1.0)),
                      pct_over_10=float(100 * np.mean(g[c] >= 10.0))))
    print(f"\n{nm}")
    print(pd.DataFrame(r).to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    # paired REM vs N2 on treatment
    if f"{tag}_rem_t90" in df.columns:
        g = df[(df[f"{tag}_rem_min"] >= MINSTAGE) & (df[f"{tag}_n2_min"] >= MINSTAGE)]
        a, b = paired(g[f"{tag}_rem_t90"].values, g[f"{tag}_n2_t90"].values)
        if a.size > 10:
            w = stats.wilcoxon(a, b)
            print(f"  paired REM vs N2 on treatment (n={a.size}): median {np.median(a):.3f} "
                  f"vs {np.median(b):.3f}, diff {np.median(a-b):+.3f}, p={w.pvalue:.3g}, "
                  f"REM worse in {100*np.mean(a>b):.1f}%")

# ====================================================================== SAVE
frames = []
for nm, df in (("A_split", A), ("P_placebo", P), ("B_pairs", BP)):
    if len(df):
        d = df.copy(); d["design"] = nm
        frames.append(d)
if frames:
    allf = pd.concat(frames, ignore_index=True)
    allf.columns = [str(c) for c in allf.columns]
    allf.to_parquet(OUT, index=False)
    print(f"\nwrote {OUT}  rows={len(allf)} cols={allf.shape[1]}")

gb = 0.0
for f in glob.glob(str(HERE / "out*_*.csv")):
    try:
        gb += pd.read_csv(f, usecols=["bytes"], low_memory=False).bytes.sum()
    except Exception:
        pass
print(f"\nS3 bytes across all reads: {gb/1e9:.2f} GB")
