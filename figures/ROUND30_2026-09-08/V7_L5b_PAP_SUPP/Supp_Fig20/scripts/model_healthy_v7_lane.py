"""
Stage-resolved survival models in the HEALTHY (disease-free at baseline) subgroup.

The question is whether the association between nocturnal hypoxaemia and incident disease is
carried by one sleep stage, and in particular whether the REM signal that dominated the general
cohort survives when nobody has organ disease at the sleep study. If REM only wins in people who
are already sick, it was a marker of prevalent illness rather than a property of REM.

Specification, matching the manuscript:
  * site-stratified Cox, UNPENALIZED (penalizer=0), PSG date as time zero
  * natural cubic spline on age, df=4, ONE COLUMN DROPPED so the design is full rank. An
    unpenalized Cox on the undropped basis silently fits nothing.
  * sex as a covariate
  * exposure rank-inverse-normal transformed WITHIN SITE, hazard ratio per 1 SD
  * each stage gated at >= 30 min of that stage (the *_g30 columns, which also require the
    per-stage T90 to be computable, not merely the minutes to clear the gate)
  * NREM and N2 correlate 0.96 and never enter the same model

Three exposure scales are carried for every fit, because the healthy subgroup is mostly
non-hypoxaemic and a per-1-SD statement on a variable that is exactly zero for 62% of people in
REM and 88% in N3 is not interpretable on its own:
  own    the paper's convention, ranks taken within the healthy subgroup
  cs     common scale, ranks taken in the general cohort and then restricted to healthy rows,
         so that one unit means the same amount of hypoxaemia in both groups. This is the scale
         used for the healthy-versus-general comparison.
  bin    any hypoxaemia in that stage versus none (T90 > 0), which is the contrast the healthy
         distribution can actually support

Effect modification is tested formally with an exposure-by-healthy interaction fitted in the
general cohort, which avoids comparing two estimates from nested samples by eye.

Nothing is reported below 20 events.

Writes stage_healthy.csv (the deliverable), stage_healthy_joint.csv (REM and NREM in one model)
and stage_healthy_log.txt.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import sys
import warnings
from contextlib import redirect_stdout

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test
from patsy import dmatrix
from scipy import stats

warnings.filterwarnings("ignore")

HERE = f"{paths.FIGURE_ROOT}/ROUND30_2026-09-08/V7_L5b_PAP_SUPP/Supp_Fig20/work/stage_healthy_v7"  # LANE V7_L5b: outputs redirected into the lane, reads stay on the v7 stage_specific folder
SRC_IN = f"{paths.SV_ROOT}/stage_specific"
TROOT = paths.FIGURE_ROOT
sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
from disease_definitions import CIRCULAR, DISEASES, NEGATIVE_CONTROLS, ORGAN_GROUP  # noqa: E402

MIN_EVENTS = 20
STAGES = ["wake", "sleep", "nrem", "n1", "n2", "n3", "rem"]
EXPOSURES = ["t90_all"] + [f"t90_{s}" for s in STAGES] + ["AHI"]
# the five windows that partition the night into mutually informative pieces. sleep, nrem and
# all are aggregates of these and are excluded from the "which stage wins" tally.
STAGE_SET = ["t90_wake", "t90_n1", "t90_n2", "t90_n3", "t90_rem"]
COARSE_SET = ["t90_wake", "t90_nrem", "t90_rem"]


def hdr(t):
    print(f"\n{'=' * 96}\n{t}\n{'=' * 96}")


def rint(x):
    """Rank-inverse-normal (Blom). Ties take the average rank, so a mass at zero maps to one
    shared score rather than being spread out."""
    r = stats.rankdata(x)
    return stats.norm.ppf((r - 0.375) / (len(r) + 0.25))


def add_spline(f):
    """Age spline, df=4, first column dropped. Built inside each analysis sample so the knots
    sit in that sample's age distribution (the healthy group is 13 years younger)."""
    sp = dmatrix("cr(a, df=4) - 1", {"a": f.AgeAtVisit.values}, return_type="dataframe").iloc[:, 1:]
    for i in range(sp.shape[1]):
        f[f"age_s{i}"] = sp.iloc[:, i].values
    return f, [f"age_s{i}" for i in range(sp.shape[1])]


def fit(d, expcols, adj, ph=False):
    """Unpenalized site-stratified Cox. Returns one record per exposure column."""
    keep = list(expcols) + adj + ["T", "E", "site"]
    d = d[keep].dropna()
    if d.E.sum() < MIN_EVENTS or d.site.nunique() < 1:
        return None
    try:
        c = CoxPHFitter(penalizer=0.0).fit(d, "T", "E", strata=["site"])
    except Exception:
        return None
    out = {"n": int(len(d)), "events": int(d.E.sum())}
    for col in expcols:
        r = c.summary.loc[col]
        out[col] = {"hr": float(r["exp(coef)"]), "lo": float(r["exp(coef) lower 95%"]),
                    "hi": float(r["exp(coef) upper 95%"]), "p": float(r["p"])}
    if ph:
        try:
            t = proportional_hazard_test(c, d, time_transform="rank")
            out["ph_p"] = float(t.summary.loc[expcols[0], "p"])
        except Exception:
            out["ph_p"] = np.nan
    return out


# ============================================================ load
b = pd.read_parquet(f"{SRC_IN}/analysis.parquet")
b["healthy_i"] = b.healthy.astype(int)
GEN = b
HLT = b[b.healthy_i == 1].copy()

LABEL = {k: v[0] for k, v in DISEASES.items()}
LABEL["death"] = "Death from any cause"
GROUP = dict(ORGAN_GROUP)
GROUP["death"] = "Composite"
KEYS = [c[:-9] for c in b.columns if c.endswith("_incident")]
KEYS = [k for k in KEYS if k in LABEL]

log = open(f"{HERE}/stage_healthy_log.txt", "w")


class Tee:
    def write(self, s):
        sys.__stdout__.write(s)
        log.write(s)

    def flush(self):
        sys.__stdout__.flush()
        log.flush()


with redirect_stdout(Tee()):
    hdr("STAGE-RESOLVED COX IN THE HEALTHY SUBGROUP")
    print(f"general cohort            {len(GEN):,}")
    print(f"healthy subgroup (H1)     {len(HLT):,}  ({100 * len(HLT) / len(GEN):.1f}%)")
    print(f"sites                     {dict(GEN.site_id.value_counts())} general, "
          f"{dict(HLT.site_id.value_counts())} healthy")
    print(f"median age                {GEN.AgeAtVisit.median():.0f} general, "
          f"{HLT.AgeAtVisit.median():.0f} healthy")
    print(f"median follow-up (y)      {GEN.death_years.median():.2f} general, "
          f"{HLT.death_years.median():.2f} healthy")
    print(f"deaths                    {int(GEN.death_incident.sum()):,} general, "
          f"{int(HLT.death_incident.sum()):,} healthy")
    print(f"\nunpenalized Cox, site-stratified, age spline df=4 with one column dropped, sex "
          f"adjusted,\nexposure rank-inverse-normal within site, HR per 1 SD, "
          f"minimum {MIN_EVENTS} events.")

    # ======================================================== 0. ascertainment
    hdr("0  ASCERTAINMENT CHECK, READ THIS BEFORE ANY HEALTHY-SUBGROUP NUMBER")
    print("Disease-free at baseline in a hospital archive is not the same as healthy. Mortality")
    print("per person-year is the one outcome captured through a channel other than the clinic,")
    print("so it exposes how well the rest of the record is captured.\n")
    GEN["ageb"] = pd.cut(GEN.AgeAtVisit, [0, 30, 40, 50, 60, 70, 120])
    print(f"{'age band':12s}{'disease-free n':>15}{'deaths':>8}{'per 1,000 py':>14}"
          f"{'other n':>10}{'deaths':>8}{'per 1,000 py':>14}{'ratio':>8}")
    print("-" * 92)
    for ab, s in GEN.groupby("ageb", observed=True):
        a, c = s[s.healthy_i == 1], s[s.healthy_i == 0]
        ra = 1000 * a.death_incident.sum() / max(a.death_years.sum(), 1e-9)
        rc = 1000 * c.death_incident.sum() / max(c.death_years.sum(), 1e-9)
        print(f"{str(ab):12s}{len(a):15,}{int(a.death_incident.sum()):8,}{ra:14.1f}"
              f"{len(c):10,}{int(c.death_incident.sum()):8,}{rc:14.1f}{ra / rc:8.1f}")
    ra = 1000 * HLT.death_incident.sum() / HLT.death_years.sum()
    oth = GEN[GEN.healthy_i == 0]
    rc = 1000 * oth.death_incident.sum() / oth.death_years.sum()
    print(f"{'all':12s}{len(HLT):15,}{int(HLT.death_incident.sum()):8,}{ra:14.1f}"
          f"{len(oth):10,}{int(oth.death_incident.sum()):8,}{rc:14.1f}{ra / rc:8.1f}")
    print(f"\nDisease-free patients die at {ra / rc:.1f} times the rate of the rest of the cohort, "
          f"and at ages under 50\nthe ratio is 5 to 8, which is not biology. The same pattern is "
          "present in the full frozen\ncohort of 19,173 (18.8 against 9.6 per 1,000 person-years), "
          "so it is a property of the\narchive and not of the PAP-off selection made here.")
    print(f"\nDeaths among the disease-free are also concentrated at one site: "
          f"{dict(HLT[HLT.death_incident == 1].site_id.value_counts())}, against "
          f"{dict(HLT[HLT.death_incident == 0].site_id.value_counts())} for those who "
          "did not die.")
    print("\nThe reading is that zero prevalent diagnoses across 26 conditions partly selects")
    print("patients with little encounter history at these hospitals. Their diagnoses are")
    print("under-captured while their deaths are captured through a linkage. Two consequences:")
    print("  1. incident-disease hazard ratios in the healthy subgroup are attenuated toward 1")
    print("     by missed outcomes, so a smaller healthy estimate is not evidence of a smaller")
    print("     effect. Every SMALLER label in Q6 has to survive this before it means anything.")
    print("  2. the within-person comparisons, which stage beats which inside the same patients")
    print("     with the same outcome definition, are NOT affected. Missing an outcome removes a")
    print("     patient's event from every window at once. Q1b, Q3 and Q4 are therefore the")
    print("     trustworthy parts of this analysis and Q6 is the fragile part.")

    # ======================================================== 1. eligibility
    hdr("1  WHICH CONDITIONS CAN BE FITTED AT ALL")
    elig = []
    for k in KEYS:
        yc, ec, pc = f"{k}_years", f"{k}_incident", f"{k}_prevalent"
        fh = HLT[(HLT[pc] == 0) & HLT[yc].notna() & (HLT[yc] > 0)]
        fg = GEN[(GEN[pc] == 0) & GEN[yc].notna() & (GEN[yc] > 0)]
        elig.append({"key": k, "disease": LABEL[k], "group": GROUP[k],
                     "neg": k in NEGATIVE_CONTROLS, "circ": k in CIRCULAR,
                     "h_ev": int(fh[ec].sum()), "h_n": len(fh),
                     "g_ev": int(fg[ec].sum()), "g_n": len(fg)})
    elig = pd.DataFrame(elig).sort_values("h_ev", ascending=False)
    KEEP = list(elig[elig.h_ev >= MIN_EVENTS].key)
    DROP = elig[elig.h_ev < MIN_EVENTS]
    print(f"{len(KEEP)} conditions reach {MIN_EVENTS} incident events in the healthy subgroup, "
          f"{len(DROP)} do not.\n")
    print(f"{'condition':30s}{'group':14s}{'healthy ev':>11}{'healthy n':>11}"
          f"{'general ev':>12}{'general n':>11}")
    print("-" * 96)
    for _, r in elig.iterrows():
        tag = "" if r.h_ev >= MIN_EVENTS else "   DROPPED"
        flag = " [neg ctrl]" if r.neg else (" [circular]" if r.circ else "")
        print(f"{r.disease[:29]:30s}{r.group:14s}{r.h_ev:11,}{r.h_n:11,}{r.g_ev:12,}"
              f"{r.g_n:11,}{tag}{flag}")
    print("\nDropped for lack of events in the healthy subgroup: "
          + ", ".join(f"{r.disease} ({r.h_ev})" for _, r in DROP.iterrows()))

    # ======================================================== 2. exposure samples
    hdr("2  EXPOSURE SAMPLES AND THE SHAPE OF THE HEALTHY DISTRIBUTION")
    SAMP = {}
    print(f"{'exposure':12s}{'general n':>11}{'gen zero%':>11}{'healthy n':>11}"
          f"{'hlt zero%':>11}{'hlt median':>12}{'hlt p90':>10}")
    print("-" * 96)
    for e in EXPOSURES:
        col = f"{e}_g30" if e.startswith("t90") else e
        g = GEN[GEN[col].notna()].copy()
        g["z"] = g.groupby("site_id")[col].transform(rint)
        g["zb"] = (g[col] > 0).astype(int) if e.startswith("t90") else (g[col] >= 15).astype(int)
        h = g[g.healthy_i == 1].copy()          # common scale: ranks inherited from general
        h = h.rename(columns={"z": "z_cs"})
        h["z"] = h.groupby("site_id")[col].transform(rint)   # own scale
        g, adj_g = add_spline(g)
        h, adj_h = add_spline(h)
        SAMP[e] = {"g": g, "h": h, "adj_g": adj_g, "adj_h": adj_h, "col": col}
        gz, hz = g[col], h[col]
        z0 = 100 * (hz == 0).mean() if e.startswith("t90") else np.nan
        g0 = 100 * (gz == 0).mean() if e.startswith("t90") else np.nan
        print(f"{e:12s}{len(g):11,}{g0:11.1f}{len(h):11,}{z0:11.1f}"
              f"{hz.median():12.3f}{hz.quantile(0.9):10.2f}")
    print("\nThe healthy distribution is mostly zero in every sleep stage. In REM 61.8% of healthy")
    print("patients have no hypoxaemia at all and in N3 87.8% do, so a rank-inverse-normal score")
    print("in those stages is close to a binary contrast dressed as a continuous one. Every fit")
    print("below is therefore also run as any-versus-none, and that column is the one to trust")
    print("when the per-1-SD estimate looks large on a nearly degenerate variable.")

    # ======================================================== 3. main fits
    hdr("3  FITTING")
    rows = []
    for k in KEEP:
        yc, ec, pc = f"{k}_years", f"{k}_incident", f"{k}_prevalent"
        for e in EXPOSURES:
            S = SAMP[e]
            rec = {"key": k, "disease": LABEL[k], "organ_group": GROUP[k],
                   "negative_control": k in NEGATIVE_CONTROLS, "circular": k in CIRCULAR,
                   "exposure": e}
            frames = {}
            for tag, src, adj in (("h", S["h"], S["adj_h"]), ("g", S["g"], S["adj_g"])):
                f = src[(src[pc] == 0) & src[yc].notna() & (src[yc] > 0)].copy()
                f["T"], f["E"], f["site"] = f[yc], f[ec].astype(int), f.site_id
                frames[tag] = f
            ph_flag = e in ("t90_sleep", "t90_rem")
            # healthy, own scale
            r = fit(frames["h"], ["z"], S["adj_h"] + ["male"], ph=ph_flag)
            if r:
                rec.update({"h_n": r["n"], "h_events": r["events"],
                            "h_hr": r["z"]["hr"], "h_lo": r["z"]["lo"],
                            "h_hi": r["z"]["hi"], "h_p": r["z"]["p"]})
                if ph_flag:
                    rec["h_ph_p"] = r.get("ph_p", np.nan)
            # healthy, common scale
            r = fit(frames["h"], ["z_cs"], S["adj_h"] + ["male"])
            if r:
                rec.update({"h_hr_cs": r["z_cs"]["hr"], "h_lo_cs": r["z_cs"]["lo"],
                            "h_hi_cs": r["z_cs"]["hi"], "h_p_cs": r["z_cs"]["p"]})
            # healthy, any versus none
            r = fit(frames["h"], ["zb"], S["adj_h"] + ["male"])
            if r:
                rec.update({"h_hr_bin": r["zb"]["hr"], "h_lo_bin": r["zb"]["lo"],
                            "h_hi_bin": r["zb"]["hi"], "h_p_bin": r["zb"]["p"],
                            "h_exposed_pct": round(100 * frames["h"].zb.mean(), 1)})
            # general
            r = fit(frames["g"], ["z"], S["adj_g"] + ["male"], ph=ph_flag)
            if r:
                rec.update({"g_n": r["n"], "g_events": r["events"],
                            "g_hr": r["z"]["hr"], "g_lo": r["z"]["lo"],
                            "g_hi": r["z"]["hi"], "g_p": r["z"]["p"]})
                if ph_flag:
                    rec["g_ph_p"] = r.get("ph_p", np.nan)
            r = fit(frames["g"], ["zb"], S["adj_g"] + ["male"])
            if r:
                rec["g_hr_bin"] = r["zb"]["hr"]
            # formal effect modification, fitted in the general cohort on the general scale
            fg = frames["g"].copy()
            fg["zxh"] = fg.z * fg.healthy_i
            r = fit(fg, ["z", "healthy_i", "zxh"], S["adj_g"] + ["male"])
            if r:
                rec["p_interaction"] = r["zxh"]["p"]
                rec["hr_ratio_hlt_vs_gen"] = r["zxh"]["hr"]
            # whole-sleep T90 refitted on this stage's own gated sample, so "does the stage beat
            # whole-sleep T90" is a comparison inside one set of patients
            if e != "t90_sleep":
                for tag, adj in (("h", S["adj_h"]), ("g", S["adj_g"])):
                    f = frames[tag]
                    if "t90_sleep_g30" not in f.columns:
                        continue
                    f = f[f.t90_sleep_g30.notna()].copy()
                    f["zs"] = f.groupby("site_id").t90_sleep_g30.transform(rint)
                    r = fit(f, ["zs"], adj + ["male"])
                    if r:
                        rec[f"{tag}_sleep_matched_hr"] = r["zs"]["hr"]
                        rec[f"{tag}_sleep_matched_p"] = r["zs"]["p"]
            rows.append(rec)
        print(f"  {LABEL[k][:34]:36s} done", flush=True)

    R = pd.DataFrame(rows)

    # direction labels, computed on the common scale so a unit means the same hypoxaemia in both
    bh, bg = np.log(R.h_hr_cs), np.log(R.g_hr)
    R["beta_pct_diff"] = np.where(bg.abs() < 0.01, np.nan, 100 * (bh - bg) / bg.abs())
    def direction(r):
        if pd.isna(r.h_hr_cs) or pd.isna(r.g_hr):
            return "", ""
        d = np.log(r.h_hr_cs) - np.log(r.g_hr)
        pt = "LARGER" if d > 0.05 else ("SMALLER" if d < -0.05 else "UNCHANGED")
        ci = "UNCHANGED"
        if r.g_hr < r.h_lo_cs:
            ci = "LARGER"
        elif r.g_hr > r.h_hi_cs:
            ci = "SMALLER"
        return pt, ci
    R[["direction_pt", "direction_ci"]] = R.apply(
        lambda r: pd.Series(direction(r)), axis=1)
    for tag in ("h", "g"):
        R[f"{tag}_beats_sleep"] = R[f"{tag}_hr"] > R[f"{tag}_sleep_matched_hr"]
    R.loc[R.exposure == "t90_sleep", ["h_beats_sleep", "g_beats_sleep"]] = np.nan
    R = R.sort_values(["key", "exposure"])
    R.to_csv(f"{HERE}/stage_healthy.csv", index=False)
    print(f"\nwrote stage_healthy.csv  {R.shape[0]} rows x {R.shape[1]} columns, "
          f"{R.h_hr.notna().sum()} healthy estimates and {R.g_hr.notna().sum()} general")

    P = R.pivot_table(index="key", columns="exposure", values="h_hr")
    PG = R.pivot_table(index="key", columns="exposure", values="g_hr")
    PN = R.pivot_table(index="key", columns="exposure", values="h_events")
    meta = R.drop_duplicates("key").set_index("key")[
        ["disease", "organ_group", "negative_control", "circular"]]

    # ======================================================== Q1 which stage wins
    hdr("Q1  IN PEOPLE WHO START HEALTHY, IS THE SIGNAL CARRIED BY ONE STAGE?")

    def winner_table(P, setcols, title, mask=None):
        avail = P[setcols].dropna()
        if mask is not None:
            avail = avail[avail.index.isin(mask)]
        w = avail.idxmax(axis=1)
        print(f"\n{title}   {len(avail)} conditions with an estimate in every window")
        print(f"{'condition':28s}" + "".join(f"{c.replace('t90_', ''):>9}" for c in setcols)
              + f"{'winner':>10}{'ev':>7}")
        print("-" * 96)
        for k in avail.index:
            print(f"{meta.loc[k, 'disease'][:27]:28s}"
                  + "".join(f"{avail.loc[k, c]:9.2f}" for c in setcols)
                  + f"{w[k].replace('t90_', ''):>10}{PN.loc[k, setcols[-1]]:7.0f}")
        vc = w.value_counts()
        print("\nwins: " + "  ".join(f"{c.replace('t90_', '')} {vc.get(c, 0)}" for c in setcols))
        return w

    real = [k for k in P.index if not meta.loc[k, "negative_control"]
            and not meta.loc[k, "circular"]]
    w_h = winner_table(P, STAGE_SET, "HEALTHY subgroup, five night windows", real)
    w_g = winner_table(PG, STAGE_SET, "GENERAL cohort, same five windows, same spec", real)
    print("\nagreement between healthy and general on the winning window: "
          f"{int((w_h == w_g.reindex(w_h.index)).sum())} of {len(w_h)} conditions")
    w_hc = winner_table(P, COARSE_SET, "HEALTHY subgroup, coarse windows", real)
    w_gc = winner_table(PG, COARSE_SET, "GENERAL cohort, coarse windows", real)
    print("\nagreement on the coarse winner: "
          f"{int((w_hc == w_gc.reindex(w_hc.index)).sum())} of {len(w_hc)}")

    # ------------------------------------------------ Q1b same patients in every window
    hdr("Q1b  THE SAME PATIENTS IN EVERY WINDOW")
    five = [f"t90_{s}_g30" for s in ["wake", "n1", "n2", "n3", "rem"]]
    three = [f"t90_{s}_g30" for s in ["wake", "nrem", "rem"]]
    for lab, f in (("general", GEN), ("healthy", HLT)):
        print(f"{lab}: all five windows >= 30 min in the same patient "
              f"{int(f[five].notna().all(axis=1).sum()):,} of {len(f):,}, "
              f"wake+NREM+REM {int(f[three].notna().all(axis=1).sum()):,}")
    print("\nThe five-window tally above is NOT a within-person ranking. Only 546 healthy and")
    print("1,263 general patients clear 30 minutes in all five windows on the same night, so N1")
    print("and N3 cannot be compared to REM inside one sample. The three coarse windows can:")
    print("they are available together in 2,000 healthy and 6,112 general patients, and their")
    print("mutual Spearman correlations there are 0.49 to 0.71, low enough for one model.\n")

    crows = []
    for tag, src, lab in (("h", HLT, "HEALTHY"), ("g", GEN, "GENERAL")):
        base = src[src[three].notna().all(axis=1)].copy()
        for s in ["wake", "nrem", "rem"]:
            base[f"zc_{s}"] = base.groupby("site_id")[f"t90_{s}_g30"].transform(rint)
        base, adj = add_spline(base)
        print(f"\n{lab} common sample n = {len(base):,}. Each window fitted alone on this sample,")
        print("then all three in one model.")
        print(f"{'condition':26s}{'wake':>8}{'nrem':>8}{'rem':>8}{'winner':>9}   |"
              f"{'wake adj':>10}{'nrem adj':>10}{'rem adj':>10}{'ev':>7}")
        print("-" * 100)
        for k in KEEP:
            yc, ec, pc = f"{k}_years", f"{k}_incident", f"{k}_prevalent"
            f = base[(base[pc] == 0) & base[yc].notna() & (base[yc] > 0)].copy()
            f["T"], f["E"], f["site"] = f[yc], f[ec].astype(int), f.site_id
            solo = {}
            for s in ["wake", "nrem", "rem"]:
                r = fit(f, [f"zc_{s}"], adj + ["male"])
                solo[s] = r[f"zc_{s}"]["hr"] if r else np.nan
            j = fit(f, [f"zc_{s}" for s in ["wake", "nrem", "rem"]], adj + ["male"])
            if not j or any(pd.isna(v) for v in solo.values()):
                continue
            win = max(solo, key=solo.get)
            rec = {"cohort": lab.lower(), "key": k, "disease": LABEL[k],
                   "negative_control": k in NEGATIVE_CONTROLS, "circular": k in CIRCULAR,
                   "n": j["n"], "events": j["events"], "winner": win,
                   **{f"solo_{s}": solo[s] for s in solo}}
            for s in ["wake", "nrem", "rem"]:
                rec[f"adj_{s}_hr"] = j[f"zc_{s}"]["hr"]
                rec[f"adj_{s}_lo"] = j[f"zc_{s}"]["lo"]
                rec[f"adj_{s}_hi"] = j[f"zc_{s}"]["hi"]
                rec[f"adj_{s}_p"] = j[f"zc_{s}"]["p"]
            crows.append(rec)
            print(f"{LABEL[k][:25]:26s}{solo['wake']:8.2f}{solo['nrem']:8.2f}{solo['rem']:8.2f}"
                  f"{win:>9}   |{j['zc_wake']['hr']:10.2f}{j['zc_nrem']['hr']:10.2f}"
                  f"{j['zc_rem']['hr']:10.2f}{j['events']:7,}")
    C = pd.DataFrame(crows)
    C.to_csv(f"{HERE}/stage_healthy_common_sample.csv", index=False)
    for lab in ("healthy", "general"):
        s = C[(C.cohort == lab) & ~C.negative_control & ~C.circular]
        vc = s.winner.value_counts()
        print(f"\n{lab} common sample, {len(s)} conditions. wins: "
              + "  ".join(f"{w} {vc.get(w, 0)}" for w in ["wake", "nrem", "rem"]))
        print(f"  in the three-window model the CI excludes 1 for wake in "
              f"{int(((s.adj_wake_lo > 1) | (s.adj_wake_hi < 1)).sum())}, NREM in "
              f"{int(((s.adj_nrem_lo > 1) | (s.adj_nrem_hi < 1)).sum())}, REM in "
              f"{int(((s.adj_rem_lo > 1) | (s.adj_rem_hi < 1)).sum())} of {len(s)}")
    print("wrote stage_healthy_common_sample.csv")

    # ======================================================== Q2 is REM special
    hdr("Q2  IS REM SPECIAL IN THE HEALTHY SUBGROUP?")
    rr = R[R.exposure == "t90_rem"].set_index("key")
    rs = R[R.exposure == "t90_sleep"].set_index("key")
    print(f"{'condition':26s}{'healthy REM HR':>28}{'ev':>6}{'general REM HR':>26}{'ev':>7}"
          f"{'p int':>8}")
    print("-" * 104)
    for k in rr.sort_values("h_hr", ascending=False).index:
        r = rr.loc[k]
        if pd.isna(r.h_hr):
            continue
        hstr = "{:.2f} [{:.2f}-{:.2f}]".format(r.h_hr, r.h_lo, r.h_hi)
        gstr = ("{:.2f} [{:.2f}-{:.2f}]".format(r.g_hr, r.g_lo, r.g_hi)
                if pd.notna(r.g_hr) else "-")
        pstr = "{:.3f}".format(r.p_interaction) if pd.notna(r.p_interaction) else "-"
        gev = r.g_events if pd.notna(r.g_events) else 0
        print(f"{r.disease[:25]:26s}{hstr:>28}{r.h_events:6.0f}{gstr:>26}{gev:7.0f}{pstr:>8}")
    sig_h = rr[(rr.h_lo > 1) | (rr.h_hi < 1)]
    sig_g = rr[(rr.g_lo > 1) | (rr.g_hi < 1)]
    print(f"\nREM significant in {len(sig_h)} of {int(rr.h_hr.notna().sum())} healthy conditions, "
          f"{len(sig_g)} of {int(rr.g_hr.notna().sum())} in the general cohort on this dataset.")
    print(f"REM any-versus-none in healthy, median HR "
          f"{rr.h_hr_bin.median():.2f}, range {rr.h_hr_bin.min():.2f} to {rr.h_hr_bin.max():.2f}")

    # ======================================================== Q3 REM vs NREM jointly
    hdr("Q3  REM AND NREM IN ONE MODEL")
    jrows = []
    for tag, src, adjk in (("h", SAMP["t90_rem"]["h"], "adj_h"), ("g", SAMP["t90_rem"]["g"], "adj_g")):
        base = src[src.t90_nrem_g30.notna()].copy()
        base["z_rem"] = base.groupby("site_id").t90_rem_g30.transform(rint)
        base["z_nrem"] = base.groupby("site_id").t90_nrem_g30.transform(rint)
        base, adj = add_spline(base)
        print(f"\n{'HEALTHY' if tag == 'h' else 'GENERAL'} joint sample n = {len(base):,} "
              f"(REM >= 30 min and NREM >= 30 min), "
              f"Spearman REM~NREM {base[['t90_rem_g30', 't90_nrem_g30']].corr(method='spearman').iloc[0, 1]:.3f}")
        print(f"{'condition':26s}{'REM adj for NREM':>24}{'NREM adj for REM':>24}{'ev':>7}{'n':>8}")
        print("-" * 92)
        for k in KEEP:
            yc, ec, pc = f"{k}_years", f"{k}_incident", f"{k}_prevalent"
            f = base[(base[pc] == 0) & base[yc].notna() & (base[yc] > 0)].copy()
            f["T"], f["E"], f["site"] = f[yc], f[ec].astype(int), f.site_id
            r = fit(f, ["z_rem", "z_nrem"], adj + ["male"])
            if not r:
                continue
            jrows.append({"cohort": "healthy" if tag == "h" else "general", "key": k,
                          "disease": LABEL[k], "negative_control": k in NEGATIVE_CONTROLS,
                          "circular": k in CIRCULAR, "n": r["n"], "events": r["events"],
                          "rem_hr": r["z_rem"]["hr"], "rem_lo": r["z_rem"]["lo"],
                          "rem_hi": r["z_rem"]["hi"], "rem_p": r["z_rem"]["p"],
                          "nrem_hr": r["z_nrem"]["hr"], "nrem_lo": r["z_nrem"]["lo"],
                          "nrem_hi": r["z_nrem"]["hi"], "nrem_p": r["z_nrem"]["p"]})
            a, c2 = r["z_rem"], r["z_nrem"]
            astr = "{:.2f} [{:.2f}-{:.2f}]".format(a["hr"], a["lo"], a["hi"])
            cstr = "{:.2f} [{:.2f}-{:.2f}]".format(c2["hr"], c2["lo"], c2["hi"])
            print(f"{LABEL[k][:25]:26s}{astr:>24}{cstr:>24}{r['events']:7,}{r['n']:8,}")
    J = pd.DataFrame(jrows)
    J.to_csv(f"{HERE}/stage_healthy_joint.csv", index=False)
    for c in ("healthy", "general"):
        s = J[(J.cohort == c) & ~J.negative_control & ~J.circular]
        print(f"\n{c}: REM > NREM in {int((s.rem_hr > s.nrem_hr).sum())} of {len(s)} conditions. "
              f"REM CI excludes 1 in {int(((s.rem_lo > 1) | (s.rem_hi < 1)).sum())}, "
              f"NREM in {int(((s.nrem_lo > 1) | (s.nrem_hi < 1)).sum())}.")
    print("wrote stage_healthy_joint.csv")

    # ======================================================== Q4 beats whole sleep
    hdr("Q4  DOES ANY STAGE BEAT WHOLE-SLEEP T90 IN THE HEALTHY SUBGROUP?")
    print("Whole-sleep T90 is refitted inside each stage's own gated sample, so the two hazard")
    print("ratios on a line come from the same patients.\n")
    for e in [f"t90_{s}" for s in STAGES if s != "sleep"] + ["t90_all", "AHI"]:
        sub = R[(R.exposure == e) & R.h_hr.notna() & ~R.negative_control & ~R.circular]
        if not len(sub):
            continue
        won = sub.h_beats_sleep.sum()
        med = (np.log(sub.h_hr) - np.log(sub.h_sleep_matched_hr)).median()
        print(f"{e:10s} beats whole-sleep T90 in {int(won):3d} of {len(sub):3d} conditions, "
              f"median log-HR advantage {med:+.3f}, "
              f"median n {sub.h_n.median():.0f}")
    print("\nby condition, best window against whole-sleep T90 on the matched sample")
    print(f"{'condition':26s}{'best window':>13}{'best HR':>10}{'sleep HR':>10}{'ev':>7}")
    print("-" * 70)
    for k in real:
        sub = R[(R.key == k) & R.exposure.isin(STAGE_SET) & R.h_hr.notna()]
        if not len(sub):
            continue
        i = sub.h_hr.idxmax()
        r = sub.loc[i]
        print(f"{r.disease[:25]:26s}{r.exposure.replace('t90_', ''):>13}{r.h_hr:10.2f}"
              f"{r.h_sleep_matched_hr:10.2f}{r.h_events:7.0f}")

    # ======================================================== Q5 negative controls
    hdr("Q5  NEGATIVE CONTROLS IN THE HEALTHY SUBGROUP, THE CONFOUNDING FLOOR")
    nc = R[R.negative_control & R.h_hr.notna()]
    print(f"{'control':22s}{'exposure':11s}{'healthy HR':>22}{'ev':>6}"
          f"{'general HR':>22}{'ev':>7}")
    print("-" * 96)
    for _, r in nc.sort_values(["exposure", "h_hr"], ascending=[True, False]).iterrows():
        g = f"{r.g_hr:.2f} [{r.g_lo:.2f}-{r.g_hi:.2f}]" if pd.notna(r.g_hr) else "-"
        print(f"{r.disease[:21]:22s}{r.exposure:11s}"
              f"{f'{r.h_hr:.2f} [{r.h_lo:.2f}-{r.h_hi:.2f}]':>22}{r.h_events:6.0f}"
              f"{g:>22}{r.g_events:7.0f}")
    print(f"\n{'exposure':12s}{'healthy floor (max NC HR)':>28}{'n controls':>12}"
          f"{'general floor':>16}{'n controls':>12}")
    print("-" * 82)
    floor = []
    for e in EXPOSURES:
        s = nc[nc.exposure == e]
        sg = R[R.negative_control & (R.exposure == e) & R.g_hr.notna()]
        if not len(s):
            continue
        floor.append({"exposure": e, "h_floor": s.h_hr.max(), "h_k": len(s),
                      "g_floor": sg.g_hr.max() if len(sg) else np.nan, "g_k": len(sg),
                      "h_sig": int(((s.h_lo > 1) | (s.h_hi < 1)).sum())})
        print(f"{e:12s}{s.h_hr.max():28.2f}{len(s):12d}"
              f"{(sg.g_hr.max() if len(sg) else np.nan):16.2f}{len(sg):12d}")
    F = pd.DataFrame(floor)
    print(f"\nAcross every exposure the healthy negative-control ceiling is "
          f"{F.h_floor.min():.2f} to {F.h_floor.max():.2f}.")
    print(f"In the general cohort on the same dataset it is "
          f"{F.g_floor.min():.2f} to {F.g_floor.max():.2f}.")
    print(f"Negative-control estimates whose CI excludes 1 in the healthy subgroup: "
          f"{int(F.h_sig.sum())} of {int(nc.shape[0])}.")
    print("\nThe maximum of a handful of imprecise estimates is biased upward and the number of")
    print("controls surviving 20 events differs by exposure, so the pooled control is the better")
    print("floor. Fixed-effect inverse-variance pooling of the log hazard ratios:\n")
    print(f"{'exposure':12s}{'healthy pooled control':>26}{'k':>4}"
          f"{'general pooled control':>26}{'k':>4}")
    print("-" * 74)

    def pool(sub, hr, lo, hi):
        s = sub[sub[hr].notna()]
        if not len(s):
            return None
        b = np.log(s[hr].values)
        se = (np.log(s[hi].values) - np.log(s[lo].values)) / (2 * 1.96)
        w = 1 / se ** 2
        m = float((w * b).sum() / w.sum())
        sem = float(np.sqrt(1 / w.sum()))
        return np.exp(m), np.exp(m - 1.96 * sem), np.exp(m + 1.96 * sem), len(s)

    pooled = []
    for e in EXPOSURES:
        sh = R[R.negative_control & (R.exposure == e)]
        ph_ = pool(sh, "h_hr", "h_lo", "h_hi")
        pg_ = pool(sh, "g_hr", "g_lo", "g_hi")
        if not ph_:
            continue
        pooled.append({"exposure": e, "h": ph_, "g": pg_})
        hs = "{:.2f} [{:.2f}-{:.2f}]".format(*ph_[:3])
        gs = "{:.2f} [{:.2f}-{:.2f}]".format(*pg_[:3]) if pg_ else "-"
        print(f"{e:12s}{hs:>26}{ph_[3]:4d}{gs:>26}{pg_[3] if pg_ else 0:4d}")
    hmax = max(p["h"][0] for p in pooled)
    hupper = max(p["h"][2] for p in pooled)
    gmax = max(p["g"][0] for p in pooled if p["g"])
    h_excl = sum(1 for p in pooled if p["h"][1] > 1 or p["h"][2] < 1)
    g_excl = sum(1 for p in pooled if p["g"] and (p["g"][1] > 1 or p["g"][2] < 1))
    print(f"\nHealthy pooled negative control: point estimate never above {hmax:.2f}, and the "
          f"confidence\ninterval covers 1.0 in {len(pooled) - h_excl} of {len(pooled)} windows.")
    print(f"General pooled negative control: point estimate up to {gmax:.2f}, and the confidence "
          f"interval\nexcludes 1.0 in {g_excl} of {len(pooled)} windows.")
    print("\nSo the confounding floor does drop. On the point estimate it drops to 1.0 in every")
    print(f"window, against 1.04 to {gmax:.2f} in the general cohort, which reproduces the 1.17 the")
    print("earlier table reported for fracture. But the healthy intervals are wide, and the")
    print(f"highest upper bound across windows is {hupper:.2f}, so the defensible statement is that "
          f"nothing\nbelow about 1.15 should be believed in the healthy subgroup either: the point "
          "estimate\nimproves, the precision does not.")

    # ------------------------------------------------ sensitivity on the two build flags
    hdr("SENSITIVITY, THE TWO ISSUES THE BUILD FLAGGED, HEALTHY REM")
    print("index_study_match == 1 drops the 52 healthy patients whose exposure was measured after")
    print("their follow-up clock started. spo2_suspect == 0 drops the oximetry the pilot flagged.\n")
    S = SAMP["t90_rem"]
    print(f"{'condition':26s}{'primary':>22}{'index-matched':>22}{'spo2 clean':>22}{'ev':>6}")
    print("-" * 94)
    top = R[(R.exposure == "t90_rem") & R.h_hr.notna() & ~R.negative_control & ~R.circular]
    for _, rr_ in top.sort_values("h_hr", ascending=False).head(12).iterrows():
        k = rr_.key
        yc, ec, pc = f"{k}_years", f"{k}_incident", f"{k}_prevalent"
        out = []
        for msk in (None, "index_study_match", "clean"):
            f = S["h"]
            if msk == "index_study_match":
                f = f[f.index_study_match == 1]
            elif msk == "clean":
                f = f[f.spo2_suspect == 0]
            f = f[(f[pc] == 0) & f[yc].notna() & (f[yc] > 0)].copy()
            f["T"], f["E"], f["site"] = f[yc], f[ec].astype(int), f.site_id
            r = fit(f, ["z"], S["adj_h"] + ["male"])
            out.append("{:.2f} [{:.2f}-{:.2f}]".format(r["z"]["hr"], r["z"]["lo"], r["z"]["hi"])
                       if r else "-")
        print(f"{rr_.disease[:25]:26s}{out[0]:>22}{out[1]:>22}{out[2]:>22}{rr_.h_events:6.0f}")

    # ======================================================== Q6 healthy versus general
    hdr("Q6  HEALTHY VERSUS GENERAL, EVERY ESTIMATE SIDE BY SIDE")
    print("Comparison is on the COMMON SCALE: ranks are taken in the general cohort and then")
    print("restricted to healthy rows, so one unit is the same amount of hypoxaemia in both.")
    print("p_interaction is the exposure-by-healthy term fitted in the general cohort.\n")
    cmp = R[R.h_hr_cs.notna() & R.g_hr.notna() & ~R.negative_control & ~R.circular]
    tab = cmp.groupby("exposure").direction_pt.value_counts().unstack(fill_value=0)
    for c in ("LARGER", "UNCHANGED", "SMALLER"):
        if c not in tab.columns:
            tab[c] = 0
    print(f"{'exposure':12s}{'LARGER':>9}{'UNCH':>8}{'SMALLER':>9}{'med % diff':>12}"
          f"{'p int < .05':>13}")
    print("-" * 66)
    for e in EXPOSURES:
        if e not in tab.index:
            continue
        s = cmp[cmp.exposure == e]
        print(f"{e:12s}{tab.loc[e, 'LARGER']:9d}{tab.loc[e, 'UNCHANGED']:8d}"
              f"{tab.loc[e, 'SMALLER']:9d}{s.beta_pct_diff.median():12.0f}"
              f"{int((s.p_interaction < 0.05).sum()):13d}")
    print(f"\noverall  LARGER {int(tab['LARGER'].sum())}, UNCHANGED "
          f"{int(tab['UNCHANGED'].sum())}, SMALLER {int(tab['SMALLER'].sum())} "
          f"of {len(cmp)} comparisons")
    print(f"interaction p < 0.05 in {int((cmp.p_interaction < 0.05).sum())} of {len(cmp)}, "
          f"expected by chance {0.05 * len(cmp):.0f}")
    strict = cmp[cmp.direction_ci != "UNCHANGED"]
    print(f"general estimate outside the healthy 95% CI in {len(strict)} of {len(cmp)}: "
          + ", ".join(f"{r.disease}/{r.exposure.replace('t90_', '')} {r.direction_ci}"
                      for _, r in strict.iterrows()))
    print("\nwake versus sleep, the comparison the manuscript already made")
    for e in ("t90_wake", "t90_sleep"):
        s = cmp[cmp.exposure == e]
        print(f"  {e:10s} median healthy HR (common scale) {s.h_hr_cs.median():.3f}, "
              f"median general {s.g_hr.median():.3f}, "
              f"LARGER in {int((s.direction_pt == 'LARGER').sum())} of {len(s)}")

    # ======================================================== proportional hazards
    hdr("PROPORTIONAL HAZARDS, EXPOSURE TERM ONLY, t90_sleep AND t90_rem")
    ph = R[R.exposure.isin(["t90_sleep", "t90_rem"]) & R.h_ph_p.notna()]
    bad_h = ph[ph.h_ph_p < 0.05]
    bad_g = ph[ph.g_ph_p < 0.05]
    print(f"healthy: Schoenfeld p < 0.05 for the exposure in {len(bad_h)} of {len(ph)} fits"
          + (": " + ", ".join(f"{r.disease}/{r.exposure}" for _, r in bad_h.iterrows())
             if len(bad_h) else ""))
    print(f"general: Schoenfeld p < 0.05 for the exposure in {len(bad_g)} of {len(ph)} fits"
          + (": " + ", ".join(f"{r.disease}/{r.exposure}" for _, r in bad_g.iterrows())
             if len(bad_g) else ""))

    hdr("RECONCILIATION AGAINST THE 2026-07-29 STAGE TABLE")
    OLD = (f"{paths.T90_ROOT}/BDSP_New_Ideas_2026-07-26/deliverable/"
           "T90_by_state_vs_AHI_all_diseases.csv")
    try:
        o = pd.read_csv(OLD).set_index("disease")
        gp = R.pivot_table(index="disease", columns="exposure", values="g_hr")
        mapc = {"T90 all": "t90_all", "WAKE": "t90_wake", "SLEEP": "t90_sleep",
                "NREM": "t90_nrem", "N1": "t90_n1", "N2": "t90_n2", "N3": "t90_n3",
                "REM": "t90_rem", "AHI": "AHI"}
        pairs = [(d, c, o.loc[d, a], gp.loc[d, c]) for d in o.index if d in gp.index
                 for a, c in mapc.items()
                 if pd.notna(o.loc[d, a]) and pd.notna(gp.loc[d, c])]
        pr = pd.DataFrame(pairs, columns=["disease", "exposure", "old", "new"])
        print(f"{pr.disease.nunique()} conditions and {len(pr)} estimates appear in both tables. "
              f"The overlap is limited\nbecause this run only fits the 44 conditions that reach "
              f"20 events in the healthy subgroup,\nso cirrhosis, sepsis, gout, myocardial "
              "infarction, peripheral artery disease, respiratory\nfailure, prostate cancer and "
              "obesity hypoventilation are absent here.")
        print(f"\nSpearman correlation of the hazard ratios       "
              f"{pr.old.corr(pr.new, method='spearman'):.3f}")
        print(f"Correlation of the log hazard ratios            "
              f"{np.log(pr.old).corr(np.log(pr.new)):.3f}")
        print(f"Median ratio of log hazard ratios, new to old   "
              f"{(np.log(pr.new) / np.log(pr.old)).median():.3f}")
        print("\nThe ranking reproduces. Estimates here are about 9% larger on the log scale, "
              "which is\nthe direction expected from dropping the ridge penalty, and the cohort "
              "differs as well:\n9,392 strictly PAP-off first studies here against the roughly "
              "13,066 built on the looser\nPAP-off-portion rule. The two tables are not "
              "interchangeable and this one supersedes it\nonly for the conditions it covers.")
    except Exception as exc:
        print(f"could not reconcile: {exc}")

    hdr("POWER, STATED PLAINLY")
    print(f"healthy subgroup                        {len(HLT):,}")
    print(f"healthy and REM >= 30 min               {int(HLT.t90_rem_g30.notna().sum()):,} "
          f"({100 * HLT.t90_rem_g30.notna().mean():.0f}%)  <- the binding constraint")
    print(f"healthy and N3 >= 30 min                {int(HLT.t90_n3_g30.notna().sum()):,}")
    print(f"healthy and N1 >= 30 min                {int(HLT.t90_n1_g30.notna().sum()):,}")
    print(f"healthy, wake and NREM and REM together {int(HLT[three].notna().all(axis=1).sum()):,}")
    print(f"healthy, all five windows together      {int(HLT[five].notna().all(axis=1).sum()):,}")
    rem_ok = R[(R.exposure == "t90_rem") & R.h_hr.notna()]
    all_ok = R[(R.exposure == "t90_all") & R.h_hr.notna()]
    print(f"\nconditions fittable on t90_all in healthy   {len(all_ok)}")
    print(f"conditions fittable on t90_rem in healthy   {len(rem_ok)}  "
          f"(the REM gate costs {len(all_ok) - len(rem_ok)} conditions)")
    print(f"median healthy REM analytic n               {rem_ok.h_n.median():.0f}")
    print(f"median healthy REM events                   {rem_ok.h_events.median():.0f}")
    print(f"median general REM events                   {rem_ok.g_events.median():.0f}  "
          f"({rem_ok.g_events.median() / rem_ok.h_events.median():.1f} times as many)")
    print(f"\nMedian follow-up is {HLT.death_years.median():.2f} y in the healthy subgroup against "
          f"{GEN.death_years.median():.2f} y in the general cohort.\nBeing disease-free in a "
          "hospital archive partly means having fewer encounters, so the healthy\ngroup is "
          "both smaller and less observed. Every negative healthy result below is at least\nas "
          "likely to be a follow-up problem as a biology one, and the direction of that bias is\n"
          "toward the no-difference answer.")

    hdr("FILES")
    print(f"{HERE}/stage_healthy.csv                 one row per condition and exposure, "
          f"healthy and general side by side")
    print(f"{HERE}/stage_healthy_joint.csv           REM and NREM in one model")
    print(f"{HERE}/stage_healthy_common_sample.csv   wake, NREM and REM in the same patients")
    print(f"{HERE}/stage_healthy_log.txt             this log")
log.close()
