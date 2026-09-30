"""
Step 015, assemble_split_prepap: data/split3622_FINAL.jsonl -> data/split_prepap_FINAL.parquet (3,622 rows; failures by id in
data/split_failures.csv), the D1 rules (AHI and arousal index NaN with _v8_index_undefined when pre-PAP true sleep < 60 min;
stage-gated EEG features NaN where a stage holds < 5 epochs is the pipeline's own rule; oxygen measures NaN when the pre-PAP
window holds < 5 valid minutes), row invariants, _v8_window_end_min, and logs/ASSEMBLE_SPLIT.md with old (whole night, v7
reextract_FINAL) against new (pre-PAP) medians per family and per site.
usage: assemble_split.py [--jsonl ...] [--manifest ...] [--out ...] [--expect-n 3622] [--log-suffix _pilot]
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import x1_paths as P
from provenance import write_sidecar
import assemble_common as AC

FAMILIES = {"oxygen": ["spo2_pct_below_90", "spo2_pct_below_90_native", "spo2_mean", "spo2_nadir", "odi3_total", "odi4_total", "spo2_pct_below_88"],
            "architecture": ["TST_min", "sleep_efficiency_pct", "N1_pct", "N2_pct", "N3_pct", "REM_pct", "WASO_min", "sol_min", "rem_latency_min", "recording_dur_min"],
            "indices": ["AHI", "arousal_index"],
            "EEG band power": ["F_nrem_delta_mean", "C_nrem_delta_mean", "O_nrem_delta_mean", "F_nrem_sigma_mean", "F_rem_theta_mean", "C_nrem_beta_kurt"],
            "spindles/SO": ["F_spindle_density_per_min", "C_spindle_density_per_min", "F_SO_density_per_min", "F_coupling_strength_r"],
            "heart": ["wholenight_rmssd", "nrem_rmssd", "wholenight_hr_bpm", "nrem_hr_bpm", "nrem_lf_hf_ratio"],
            "new": ["hypoxic_burden", "plm_index", "lm_index", "spo2_pct_below_90_sleep"]}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--jsonl", default=f"{P.X1_DATA}/split3622_FINAL.jsonl"); ap.add_argument("--manifest", default=P.MANIFEST_SPLIT)
    ap.add_argument("--out", default=f"{P.X1_DATA}/split_prepap_FINAL.parquet"); ap.add_argument("--expect-n", type=int, default=None); ap.add_argument("--log-suffix", default="")
    a = ap.parse_args()
    P.assert_not_evicted(a.jsonl, a.manifest, P.REEXTRACT_FINAL, P.T90_FINAL_V7)
    man = pd.read_csv(a.manifest)
    df = AC.load_jsonl(a.jsonl)
    coh = P.cohort_frame(["_v7_1_stage_collapsed"]); collapsed = set(coh[coh._v7_1_stage_collapsed.astype(bool)].BDSPPatientID)
    keep = [c for c in ("BDSPPatientID", "site_id", "bids_col", "SessionID", "_encoding_used", "t_end_sec", "cpap_start_min", "_v8_session_differs_from_v7") if c in man]
    df = man[keep].merge(df.drop(columns=[c for c in ("site_id", "_encoding_used") if c in df]), on="BDSPPatientID", how="left")
    df["status"] = df.status.fillna("missing_record")
    expect = a.expect_n or len(man)
    assert len(df) == expect and not df.BDSPPatientID.duplicated().any(), f"{len(df)} rows, expected {expect}"
    df = AC.apply_rules(df, collapsed)
    df["_v8_prepap_source"] = f"{os.path.basename(a.jsonl)} sha256:{AC.sha256(a.jsonl)[:16]}"
    df["_v8_window_end_min"] = pd.to_numeric(df.get("_v8_window_end_min"), errors="coerce")
    fails = df[~df.status.astype(str).str.startswith("ok")][["BDSPPatientID", "site_id", "status"]]
    fails.to_csv(f"{P.X1_DATA}/split_failures{a.log_suffix}.csv", index=False)
    ok = df[df.status.astype(str).str.startswith("ok")]
    win = pd.to_numeric(ok.get("_v8_window_end_min"), errors="coerce"); cs = pd.to_numeric(ok.get("cpap_start_min"), errors="coerce")
    L = [f"# ASSEMBLE_SPLIT ({len(df)} rows from {os.path.basename(a.jsonl)}; {len(fails)} failures in split_failures{a.log_suffix}.csv = {100*len(fails)/max(len(df),1):.2f}%)",
         f"status: {df.status.astype(str).str[:40].value_counts().to_dict()}",
         f"windowed rows among ok: {int(ok._v8_prepap_window.fillna(False).astype(bool).sum())} of {len(ok)}; window end vs cpap_start_min max |diff| {float((win - cs).abs().max()) if len(ok) else float('nan'):.4f} min",
         f"D1: pre-PAP TST < 60 min -> AHI/arousal index NaN: {int(df._v8_index_undefined.sum())} (plan expected about 478 of 3,622); < 30 min: {int((pd.to_numeric(df.TST_min, errors='coerce') < 30).sum())}; "
         f"oxygen window < 5 valid min: {int(df._v8_oxygen_window_short.sum())}; collapsed (D11): {int(df._v7_1_stage_collapsed.sum())}",
         "\n## Row invariants (ok rows)"]
    inv = AC.invariants(df)
    for name, npass, nchk, ids in inv:
        L.append(f"- {name}: {npass} of {nchk}{'' if npass == nchk else '  FAIL ids ' + str(ids)}")
    fin = pd.read_parquet(P.REEXTRACT_FINAL); fin = fin[fin.BDSPPatientID.isin(ok.BDSPPatientID)]
    L.append("\n## Old (whole night, v7 reextract_FINAL) against new (pre-PAP window) medians, per family and site")
    for fam, cols in FAMILIES.items():
        L.append(f"### {fam}")
        for c in cols:
            if c not in ok:
                continue
            for site in ("I0002", "I0006"):
                n = pd.to_numeric(ok[ok.site_id == site][c], errors="coerce"); o = pd.to_numeric(fin[fin.site_id == site][c], errors="coerce") if c in fin else pd.Series(dtype=float)
                L.append(f"- {c} {site}: new n={int(n.notna().sum())} median {n.median():.4g} IQR [{n.quantile(.25):.4g}, {n.quantile(.75):.4g}] | old n={int(o.notna().sum())} median {o.median() if len(o) else float('nan'):.4g}"
                         + (f" | missing new but old finite: {int((n.isna() & o.reindex(n.index).notna()).sum()) if len(o) else 0}" if c in fin else " | old n/a"))
    rep = "\n".join(L); open(f"{P.X1_LOGS}/ASSEMBLE_SPLIT{a.log_suffix}.md", "w").write(rep + "\n"); print(rep)
    df.to_parquet(a.out, index=False)
    idx = f"{P.X1_LOGS}/PROVENANCE_INDEX.tsv"
    for o in (a.out, f"{P.X1_LOGS}/ASSEMBLE_SPLIT{a.log_suffix}.md", f"{P.X1_DATA}/split_failures{a.log_suffix}.csv"):
        write_sidecar(o, [a.jsonl, a.manifest, P.REEXTRACT_FINAL, P.T90_FINAL_V7], __file__, "015_assemble_split", index_path=idx)
    fail_inv = [n for n, p, k, _ in inv if p != k]
    print(f"\nASSEMBLE_SPLIT: failures {len(fails)} ({100*len(fails)/max(len(df),1):.2f}%); invariant failures {fail_inv}; wrote {a.out}")
    return 0 if (not fail_inv and len(fails) < 0.01 * len(df)) else 1


if __name__ == "__main__":
    sys.exit(main())
