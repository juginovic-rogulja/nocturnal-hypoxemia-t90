"""
Step 013 checker (after the EC2 pilot on the 300 split nights + the 12 probe nights at zero cut):
  * split nights: pre-PAP whole-recording native T90 vs cpap_t90_by_stage.pre_any_t90 (Spearman > 0.99, median |diff| < 0.5),
    pre-PAP sleep minutes vs pre_sleep_min (|diff| <= 1 min on >= 95 percent), stage percentages sum to 100, SE <= 100,
    AHI x TST/60 integer, per-site distributions of every family beside the whole-night v7 values (reextract_FINAL);
  * probe nights (t_end empty): every shared column within 1e-6 of reextract_FINAL on the same platform (x86) as v7.
Writes logs/PILOT300_SPLIT.md with a sidecar. Exit 1 when a gate fails.
usage: check_pilot_split.py [--jsonl data/pilot300_split_FINAL.jsonl]
"""
from __future__ import annotations
import argparse, json, math, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import x1_paths as P
from provenance import write_sidecar
from compare_to_reextract_final import compare, TOL

FAMILIES = {"oxygen": ["spo2_pct_below_90", "spo2_pct_below_90_native", "spo2_mean", "spo2_nadir", "odi3_total", "odi4_total", "spo2_pct_below_88"],
            "architecture": ["TST_min", "sleep_efficiency_pct", "N1_pct", "N2_pct", "N3_pct", "REM_pct", "WASO_min", "sol_min", "rem_latency_min"],
            "indices": ["AHI", "arousal_index"],
            "brain": ["F_nrem_delta_mean", "C_nrem_delta_mean", "O_nrem_delta_mean", "F_rem_theta_mean", "F_nrem_sigma_mean"],
            "micro": ["F_spindle_density_per_min", "C_SO_density_per_min", "F_coupling_strength_r"],
            "heart": ["wholenight_rmssd", "nrem_rmssd", "wholenight_hr_bpm", "nrem_lf_hf_ratio"],
            "new": ["hypoxic_burden", "plm_index", "lm_index", "spo2_pct_below_90_sleep"]}


def q(s):
    s = pd.to_numeric(s, errors="coerce")
    return f"n={int(s.notna().sum())} med={s.median():.3g} [{s.quantile(.25):.3g},{s.quantile(.75):.3g}]"


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--jsonl", default=f"{P.X1_DATA}/pilot300_split_FINAL.jsonl")
    ap.add_argument("--report", default=f"{P.X1_LOGS}/PILOT300_SPLIT.md")
    a = ap.parse_args()
    P.assert_not_evicted(a.jsonl, P.REEXTRACT_FINAL, P.CPAP_STAGE_V7)
    df = pd.DataFrame([json.loads(l) for l in open(a.jsonl)]).drop_duplicates("BDSPPatientID", keep="last")
    L = [f"# PILOT300_SPLIT ({len(df)} records)", f"status: {df.status.astype(str).str[:40].value_counts().to_dict()}"]
    ok = df[df.status == "ok"].copy(); ok["is_split"] = ok["_v8_prepap_window"].astype(bool)
    fails = []
    cp = pd.read_parquet(P.CPAP_STAGE_V7, columns=["BDSPPatientID", "pre_any_t90", "pre_sleep_min", "cpap_start_min"])
    fin = pd.read_parquet(P.REEXTRACT_FINAL)
    sp = ok[ok.is_split].merge(cp, on="BDSPPatientID")
    from scipy.stats import spearmanr
    x = pd.to_numeric(sp.spo2_pct_below_90_native, errors="coerce"); y = sp.pre_any_t90; m = x.notna() & y.notna()
    rho = spearmanr(x[m], y[m]).correlation; med = float((x - y)[m].abs().median())
    L.append(f"\n## Gate A: pre-PAP native T90 vs pre_any_t90 on {int(m.sum())} split nights: Spearman {rho:.4f} (>0.99 needed), median|diff| {med:.3f} (<0.5), p95 {float((x-y)[m].abs().quantile(.95)):.3f}, max {float((x-y)[m].abs().max()):.3f}")
    if not (rho > 0.99 and med < 0.5):
        fails.append("gate A")
    d = (pd.to_numeric(sp.TST_min, errors="coerce") - sp.pre_sleep_min).abs(); frac = float((d <= 1).mean() * 100)
    L.append(f"## Gate B: pre-PAP TST vs pre_sleep_min: within 1 min on {frac:.1f}% (>=95 needed), median {float(d.median()):.3f}, max {float(d.max()):.2f}")
    if frac < 95:
        fails.append("gate B")
    L.append(sp.assign(dT90=(x - y).abs(), dTST=d).sort_values("dTST", ascending=False)[["BDSPPatientID", "site_id", "TST_min", "pre_sleep_min", "spo2_pct_below_90_native", "pre_any_t90", "_v8_window_end_min", "cpap_start_min"]].head(8).to_string(index=False))
    pct = sp[["N1_pct", "N2_pct", "N3_pct", "REM_pct"]].apply(pd.to_numeric, errors="coerce").sum(axis=1, min_count=1)
    n100 = int(((pct - 100).abs() < 1e-6).sum()); L.append(f"## Invariants: stage pct sum 100 on {n100} of {int(pct.notna().sum())} (TST>0)")
    se = pd.to_numeric(sp.sleep_efficiency_pct, errors="coerce"); L.append(f"- SE <= 100: {int((se <= 100 + 1e-9).sum())} of {int(se.notna().sum())}")
    ahi_n = pd.to_numeric(sp.AHI, errors="coerce") * pd.to_numeric(sp.TST_min, errors="coerce") / 60
    L.append(f"- AHI x TST/60 integer: {int((ahi_n.dropna() - ahi_n.dropna().round()).abs().lt(1e-6).sum())} of {int(ahi_n.notna().sum())}")
    L.append(f"- pre-PAP TST < 60 min (index undefined under D1): {int((pd.to_numeric(sp.TST_min, errors='coerce') < 60).sum())}; < 30: {int((pd.to_numeric(sp.TST_min, errors='coerce') < 30).sum())}")
    if n100 != int(pct.notna().sum()) or int((se <= 100 + 1e-9).sum()) != int(se.notna().sum()):
        fails.append("invariants")
    L.append("\n## Per-site distributions, pre-PAP (new) beside whole-night v7 (reextract_FINAL), split nights")
    v7 = fin[fin.BDSPPatientID.isin(sp.BDSPPatientID)]
    for fam, cols in FAMILIES.items():
        for c in cols:
            if c in sp.columns:
                L.append(f"- {fam}/{c}: I0002 new {q(sp[sp.site_id=='I0002'][c])} v7 {q(v7[v7.site_id=='I0002'][c]) if c in v7 else 'n/a'} | "
                         f"I0006 new {q(sp[sp.site_id=='I0006'][c])} v7 {q(v7[v7.site_id=='I0006'][c]) if c in v7 else 'n/a'}")
    pr = ok[~ok.is_split]
    L.append(f"\n## Gate C: {len(pr)} probe nights at zero cut vs reextract_FINAL on x86 (1e-6)")
    if len(pr):
        dd, nc, nm, ncol = compare(pr, fin[fin.BDSPPatientID.isin(pr.BDSPPatientID)], "reextract_FINAL_x86")
        L.append(f"- {ncol} columns, {nc} cells, {nm} within {TOL} ({100*nm/max(nc,1):.3f}%), {len(dd)} mismatches: {dd['class'].value_counts().to_dict() if len(dd) else {}}")
        if len(dd):
            dd = dd.assign(rel=dd.absdiff / np.maximum(dd.old.abs(), 1e-12))
            L.append(f"- mismatch relative diff: median {dd.rel.median():.2e} max {dd.rel.max():.2e}; columns: {sorted(dd.column.unique())[:30]}")
            dd.to_csv(f"{P.X1_LOGS}/probe12_x86_reproduction.csv", index=False)
        if nm != nc:
            fails.append("gate C (see relative diffs)")
    L.append(f"\n## Timing: extract_wall_sec median {pd.to_numeric(ok.extract_wall_sec, errors='coerce').median():.0f}, p95 {pd.to_numeric(ok.extract_wall_sec, errors='coerce').quantile(.95):.0f}; total_sec median {pd.to_numeric(ok.total_sec, errors='coerce').median():.0f}")
    L.append(f"\nGATES FAILED: {fails if fails else 'none'}")
    rep = "\n".join(L); open(a.report, "w").write(rep + "\n"); print(rep)
    write_sidecar(a.report, [a.jsonl, P.REEXTRACT_FINAL, P.CPAP_STAGE_V7], __file__, "013_check_pilot_split", index_path=f"{P.X1_LOGS}/PROVENANCE_INDEX.tsv")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
