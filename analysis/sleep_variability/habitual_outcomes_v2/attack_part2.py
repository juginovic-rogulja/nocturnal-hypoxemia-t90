"""
Continuation of attack.py: the cross-validated concordance of attack 6, and the summary file.

attack.py ran all six attacks and wrote every CSV except attack_concordance.csv. It was stopped
inside the last block, the five-fold held-out concordance, because this Mac was running six other
jobs at once, eight cores against a load average of 430, and a Cox fit that takes 0.6 s alone was
taking twenty. Its transcript is attack_run.out.

This file finishes that one block and assembles attack_summary.json from the CSVs attack.py
already wrote, refitting nothing else. Two changes from attack.py's version of the block, both
forced by the machine load and both stated in the transcript: three folds rather than five, and
the twelve outcomes with the most incident events rather than all fifty. Neither changes what the
comparison is, because every exposure is scored on the same folds and the same outcomes.

Model, unchanged: unpenalized Cox, stratified by hospital, natural cubic age spline
cr(age, df=4) with the FIRST COLUMN DROPPED, plus sex.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import json
import os
import sys
import time
import warnings

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
from cohort_spec import apply_cohort                     # noqa: E402
from disease_definitions import DISEASES, NEGATIVE_CONTROLS, CIRCULAR   # noqa: E402

T0 = time.time()
LOG = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.append(s)


b = apply_cohort(pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet"))
h = pd.read_parquet(f"{HAB}/habitual_per_patient_v2.parquet")[
    ["BDSPPatientID", "tier", "hours", "n_records", "n_tier", "n_pre_study"]]
b = b.merge(h, on="BDSPPatientID", how="left", validate="one_to_one")
b["has_habitual"] = b.hours.notna().astype(int)
assert int(b.has_habitual.sum()) == 5295   # v2 habitual set (v1 was 8711)
b["male"] = (b.sex.astype(str).str.upper().str[0] == "M").astype(int)
sp = dmatrix("cr(a, df=4) - 1", {"a": b.AgeAtVisit.values}, return_type="dataframe").iloc[:, 1:]
for i in range(sp.shape[1]):
    b[f"age_s{i}"] = sp.iloc[:, i].values
ADJ = [f"age_s{i}" for i in range(sp.shape[1])] + ["male"]
KEYS = [c[:-9] for c in b.columns if c.endswith("_incident")]
OUTCOMES = [k for k in KEYS if f"{k}_years" in b.columns]
NONCIRC = [k for k in OUTCOMES if k not in CIRCULAR]
LABEL = {k: (DISEASES[k][0] if k in DISEASES else "Death from any cause") for k in OUTCOMES}


def rint(x):
    r = stats.rankdata(x)
    return stats.norm.ppf((r - 0.375) / (len(r) + 0.25))


D = b[b.has_habitual == 1].copy()
D["z_t90"] = D.groupby("site_id").spo2_pct_below_90.transform(rint)
D["z_tst"] = D.groupby("site_id").TST_min.transform(lambda x: rint(x.fillna(x.median())))
D["z_ahi"] = D.groupby("site_id").AHI.transform(lambda x: rint(x.fillna(x.median())))
D["z_hab"] = (D.hours - D.hours.mean()) / D.hours.std()
say(f"habitual analysis set {len(D):,}, age spline columns {sp.shape[1]} of 4")


def frame(d, key):
    yc, ec, pc = f"{key}_years", f"{key}_incident", f"{key}_prevalent"
    f = d[(d[pc] == 0) & d[yc].notna() & (d[yc] > 0)].copy()
    f["T"] = f[yc]
    f["E"] = f[ec].astype(int)
    return f


ev = {k: int(frame(D, k).E.sum()) for k in NONCIRC}
CVSET = [k for k, v in sorted(ev.items(), key=lambda x: -x[1])[:12]]
say("")
say("CROSS-VALIDATED HELD-OUT CONCORDANCE, 3-fold, fresh resampling, NOT the published "
    "cross-hospital protocol.")
say("outcomes, the twelve with the most incident events: " +
    ", ".join(f"{LABEL[k]} {ev[k]}" for k in CVSET))


def cv_concordance(expo, k=3, seed=17):
    rng = np.random.RandomState(seed)
    out = {}
    for key in CVSET:
        f = frame(D, key)
        cols = ADJ + ["T", "E", "site_id"] + ([expo] if expo else [])
        f = f[cols].dropna()
        idx = rng.permutation(len(f))
        folds = np.array_split(idx, k)
        cs, ws = [], []
        for i in range(k):
            te = f.iloc[folds[i]]
            tr = f.iloc[np.concatenate([folds[j] for j in range(k) if j != i])]
            try:
                c = CoxPHFitter().fit(tr, "T", "E", strata=["site_id"])
                pred = c.predict_partial_hazard(te)
            except Exception:
                continue
            for site, g in te.groupby("site_id"):
                if g.E.sum() < 3:
                    continue
                try:
                    cs.append(concordance_index(g["T"], -pred.loc[g.index], g.E))
                    ws.append(g.E.sum())
                except Exception:
                    pass
        if cs:
            out[key] = float(np.average(cs, weights=ws))
    return out


base = cv_concordance(None)
conc = pd.DataFrame({"key": list(base), "disease": [LABEL[k] for k in base],
                     "events": [ev[k] for k in base], "c_base": list(base.values())})
for expo in ("z_hab", "z_tst", "z_t90", "z_ahi"):
    e = cv_concordance(expo)
    conc[f"c_{expo}"] = conc.key.map(e)
    conc[f"gain_{expo}"] = conc[f"c_{expo}"] - conc.c_base
conc.to_csv(f"{OUT}/attack_concordance.csv", index=False)
say(f"outcomes scored {len(conc)}, mean baseline C (age + sex, hospital-stratified) "
    f"{conc.c_base.mean():.4f}")
for expo in ("z_hab", "z_tst", "z_t90", "z_ahi"):
    g = conc[f"gain_{expo}"].dropna()
    say(f"  {expo:6s} mean gain {g.mean():+.6f}   median {g.median():+.6f}   improves "
        f"{int((g>0).sum())} of {len(g)}   best {g.max():+.4f}   worst {g.min():+.4f}")
wht = stats.wilcoxon(conc.gain_z_hab, conc.gain_z_tst)
wh9 = stats.wilcoxon(conc.gain_z_hab, conc.gain_z_t90)
say(f"  habitual vs laboratory total sleep time, paired Wilcoxon p={wht.pvalue:.3f} "
    f"(habitual larger in {int((conc.gain_z_hab>conc.gain_z_tst).sum())} of {len(conc)})")
say(f"  habitual vs laboratory T90, paired Wilcoxon p={wh9.pvalue:.5f} "
    f"(T90 larger in {int((conc.gain_z_t90>conc.gain_z_hab).sum())} of {len(conc)})")
say(conc[["disease", "events", "c_base", "gain_z_hab", "gain_z_tst", "gain_z_t90"]]
    .to_string(index=False, float_format=lambda x: f"{x:.4f}"))


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    v = p[o]
    n = len(v)
    adj = np.minimum.accumulate((v * n / np.arange(1, n + 1))[::-1])[::-1]
    q = np.empty(n)
    q[o] = adj
    return q


# ----------------------------------------------------------------- summary, from the CSVs
smd = pd.read_csv(f"{OUT}/attack_selection_smd.csv")
sel = pd.read_csv(f"{OUT}/attack_selection_ishabitual.csv")
sel["q"] = bh(sel.p.values)
tt = pd.read_csv(f"{OUT}/attack_t90_reproduction.csv")
mf = tt[tt.set == "full_19173"].merge(tt[tt.set == "habitual_8711"], on="key",
                                      suffixes=("_f", "_h"))
mf = mf[~mf.circular_f]
lf, lh = np.log(mf.hr_f), np.log(mf.hr_h)
het = pd.read_csv(f"{OUT}/attack_tier_heterogeneity.csv")
tier = pd.read_csv(f"{OUT}/attack_tier.csv")
tim = pd.read_csv(f"{OUT}/attack_note_timing.csv")
lm = pd.read_csv(f"{OUT}/attack_landmark.csv")
ps = pd.read_csv(f"{OUT}/attack_pre_vs_post.csv")
ncf = pd.read_csv(f"{OUT}/attack_notecount.csv")
pw = pd.read_csv(f"{OUT}/attack_power_cells.csv")
h2h = pd.read_csv(f"{OUT}/attack_headtohead.csv")

a = tier[tier.set == "tierA"].set_index("key")
bb = tier[tier.set == "tierB"].set_index("key")
j = a.join(bb, rsuffix="_B", how="inner")
j = j[~j.circular]

RES = {
    "attack1_selection": dict(
        n_smd_over_010=int((smd.abs_smd > 0.10).sum()), n_smd_vars=int(len(smd)),
        max_smd_var=str(smd.sort_values("abs_smd", ascending=False).iloc[0]["var"]),
        max_smd=float(smd.abs_smd.max()),
        has_habitual_significant=int((sel.p < .05).sum()),
        has_habitual_bh=int((sel.q < .05).sum()),
        has_habitual_median_hr=float(sel.hr.median()),
        has_habitual_n_outcomes=int(len(sel)),
        t90_logHR_r=float(np.corrcoef(lf, lh)[0, 1]),
        t90_slope=float(np.polyfit(lf, lh, 1)[0]),
        t90_sig_full=int((mf.p_f < .05).sum()), t90_sig_sub=int((mf.p_h < .05).sum()),
        t90_n=int(len(mf)), t90_median_full=float(mf.hr_f.median()),
        t90_median_sub=float(mf.hr_h.median())),
    "attack2_measurement": dict(
        het_sig=int((het.p_het < .05).sum()), het_bh=int((het.q_het < .05).sum()),
        het_n=int(len(het)),
        tierAB_r=float(np.corrcoef(np.log(j.hr), np.log(j.hr_B))[0, 1]),
        tierAB_sign_agree=int((np.sign(np.log(j.hr)) == np.sign(np.log(j.hr_B))).sum()),
        tierAB_n=int(len(j))),
    "attack3_reverse": dict(
        median_pct_note_after_dx=float(tim.pct_note_after_dx.median()),
        max_pct_note_after_dx=float(tim.pct_note_after_dx.max()),
        prepost_differ=int((ps.p_diff < .05).sum()), prepost_n=int(len(ps))),
    "attack4_notecount": dict(
        notes_sig=int((ncf[~ncf.circular].notes_p < .05).sum()),
        notes_q_sig=int((ncf[~ncf.circular].notes_q < .05).sum()),
        notes_n=int((~ncf.circular).sum()),
        notes_median_hr=float(ncf[~ncf.circular].notes_hr.median()),
        notes_sig_controls=int((ncf[ncf.control].notes_p < .05).sum())),
    "attack5_power": dict(
        contrasts=int((~pw.circular).sum()),
        significant=int(((~pw.circular) & (pw.p < .05)).sum()),
        no_information=int(((~pw.circular) & (pw.informative == "UNINFORMATIVE")).sum()),
        median_mdhr80=float(pw[~pw.circular].mdhr80.median())),
    "attack6_headtohead": dict(
        concordance={expo: dict(mean_gain=float(conc[f"gain_{expo}"].mean()),
                                improved=int((conc[f"gain_{expo}"] > 0).sum()),
                                n=int(conc[f"gain_{expo}"].notna().sum()))
                     for expo in ("z_hab", "z_tst", "z_t90", "z_ahi")},
        baseline_c=float(conc.c_base.mean()),
        wilcoxon_hab_vs_tst=float(wht.pvalue), wilcoxon_hab_vs_t90=float(wh9.pvalue)),
}
for L in (0.0, 1.0, 2.0):
    for expo in ("z_hab", "z_t90"):
        s = lm[(lm.set == "habitual_8711") & (lm.exposure == expo) &
               (lm.landmark_years == L) & (~lm.circular)]
        RES["attack3_reverse"][f"{expo}_LM{int(L)}_significant"] = int((s.p < .05).sum())
for expo in ("z_hab", "z_tst", "z_t90"):
    s = h2h[(h2h.exposure == expo) & (~h2h.circular)]
    ct, nc = s[s.control], s[~s.control]
    floor = float(ct.abs_hr.max())
    RES["attack6_headtohead"][expo] = dict(
        n=int(len(nc)), significant=int((nc.p < .05).sum()), bh=int((nc.q < .05).sum()),
        above_floor=int(((nc.p < .05) & (nc.abs_hr > floor)).sum()),
        median_abs_hr=float(nc.abs_hr.median()), max_abs_hr=float(nc.abs_hr.max()),
        control_floor=floor, control_min=float(ct.hr.min()), control_max=float(ct.hr.max()))

RES["_runtime_s"] = round(time.time() - T0, 1)
RES["_built"] = time.strftime("%Y-%m-%d %H:%M:%S")
with open(f"{OUT}/attack_summary.json", "w") as fh:
    json.dump(RES, fh, indent=2, default=str)
with open(f"{OUT}/attack_log_part2.txt", "w") as fh:
    fh.write("\n".join(LOG))
say("")
say(json.dumps(RES, indent=2, default=str))
say(f"runtime {RES['_runtime_s']}s")
