"""
Step 023, assemble_light: data/light19173_FINAL.jsonl -> data/light_FINAL.parquet (one row per cohort night; failures listed in
data/light_failures.csv), with the D1/D11 rules and the provenance columns, plus logs/ASSEMBLE_LIGHT.md (invariants and the
positive controls against v7 on non-split nights) and logs/CROSS_SITE_NEW_MEASURES.md (gate 4, per-site distribution of every
new measure with the convention note; a difference without a note is a STOP for the owner).
usage: assemble_light.py [--jsonl ...] [--manifest ...] [--out ...] [--expect-n 19173]
"""
from __future__ import annotations
import argparse, json, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import x1_paths as P
from provenance import write_sidecar
import assemble_common as AC

NOTES = {"plm_index": "B3: I0002 scores every limb movement (codes 1/2/3/4), I0006 far fewer rows per night; the series rule is applied to both; the ranking uses the within-site rank-normal score (D6)",
         "lm_index": "B3 as for plm_index",
         "hypoxic_burden": "D5: event table differs by site (I0002 resp_3 with the arousal hypopnea rule; I0006 single table with RERA rows); oximeter 25 vs 10 Hz; windows per site from the ensemble curve",
         "sol_min": "staging start differs by site (I0006 stages from sample 0 with an unscored code; I0002 from lights-off), see the integrity rules",
         "rem_latency_min": "same staging-start note as sol_min; collapsed nights (D11) are NaN",
         "spo2_pct_below_90_sleep": "oximeter grid differs by site (integer plateau at 90 at one site); sleep-coded native samples only"}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--jsonl", default=f"{P.X1_DATA}/light19173_FINAL.jsonl"); ap.add_argument("--manifest", default=P.MANIFEST_FULL_V8)
    ap.add_argument("--out", default=f"{P.X1_DATA}/light_FINAL.parquet"); ap.add_argument("--expect-n", type=int, default=None)
    ap.add_argument("--log-suffix", default="")
    a = ap.parse_args()
    P.assert_not_evicted(a.jsonl, a.manifest, P.REEXTRACT_FINAL, P.T90_FINAL_V7)
    man = pd.read_csv(a.manifest)
    df = AC.load_jsonl(a.jsonl)
    coh = P.cohort_frame(["_v7_1_stage_collapsed"]); collapsed = set(coh[coh._v7_1_stage_collapsed.astype(bool)].BDSPPatientID)
    expect = a.expect_n or len(man)
    df = man[["BDSPPatientID", "site_id", "bids_col", "SessionID", "_encoding_used", "t_end_sec"]].merge(df.drop(columns=[c for c in ("site_id", "_encoding_used") if c in df]), on="BDSPPatientID", how="left")
    df["status"] = df.status.fillna("missing_record")
    df = AC.apply_rules(df, collapsed)
    hb_sha = AC.sha256(P.HB_DEFINITION)[:16] if os.path.isfile(P.HB_DEFINITION) else None
    df["_v8_new_measure_source"] = f"{os.path.basename(a.jsonl)} sha256:{AC.sha256(a.jsonl)[:16]}"
    df["_v8_hb_definition_sha"] = hb_sha
    df["_v8_light_status"] = df.status
    df["_v8_window_end_min"] = pd.to_numeric(df.get("_v8_window_end_min"), errors="coerce")
    fails = df[~df.status.astype(str).str.startswith("ok")][["BDSPPatientID", "site_id", "status"]]
    fails.to_csv(f"{P.X1_DATA}/light_failures{a.log_suffix}.csv", index=False)
    assert len(df) == expect, f"{len(df)} rows, expected {expect}"
    assert not df.BDSPPatientID.duplicated().any()
    L = [f"# ASSEMBLE_LIGHT ({len(df)} rows from {os.path.basename(a.jsonl)}; {len(fails)} failures listed in light_failures{a.log_suffix}.csv)",
         f"status: {df.status.astype(str).str[:40].value_counts().to_dict()}",
         f"windowed (split) rows: {int(df._v8_prepap_window.fillna(False).astype(bool).sum())}; index undefined (D1, pre-PAP TST < 60): {int(df._v8_index_undefined.sum())}; "
         f"oxygen window short (< 5 valid min): {int(df._v8_oxygen_window_short.sum())}; collapsed staging (D11): {int(df._v7_1_stage_collapsed.sum())}",
         "\n## Positive controls on non-split ok nights (exact equality with v7 / reextract_FINAL)"]
    fin = pd.read_parquet(P.REEXTRACT_FINAL, columns=["BDSPPatientID", "sol_min", "rem_latency_min", "spo2_pct_below_90_native", "spo2_pct_below_90", "TST_min"])
    ns = df[df.status.astype(str).str.startswith("ok") & ~df._v8_prepap_window.fillna(False).astype(bool)].merge(fin, on="BDSPPatientID", suffixes=("", "_v7"))
    n_bad = 0
    for c in ("spo2_pct_below_90_native", "spo2_pct_below_90", "sol_min", "rem_latency_min", "TST_min"):
        x = pd.to_numeric(ns[c], errors="coerce"); y = pd.to_numeric(ns[f"{c}_v7"], errors="coerce")
        # D11 sets sol/rem NaN on collapsed nights in v8 while v7 kept a value: compare on non-collapsed nights only for those
        m = ~ns._v7_1_stage_collapsed if c in ("sol_min", "rem_latency_min") else pd.Series(True, index=ns.index)
        eq = ((x == y) | (x.isna() & y.isna()))[m]
        n_bad += int((~eq).sum())
        L.append(f"- {c}: {int(eq.sum())} of {int(m.sum())} identical; mismatched ids: {ns.loc[m & ~((x == y) | (x.isna() & y.isna())), 'BDSPPatientID'].tolist()[:10]}")
    L.append("\n## Row invariants (ok rows)")
    inv = AC.invariants(df)
    for name, npass, nchk, ids in inv:
        L.append(f"- {name}: {npass} of {nchk}{'' if npass == nchk else '  FAIL ids ' + str(ids)}")
    L.append(f"\n## D11 check: stage-gated new measures NaN on the {int(df._v7_1_stage_collapsed.sum())} collapsed nights: " +
             str({c: int(df.loc[df._v7_1_stage_collapsed, c].notna().sum()) for c in AC.NEW_MEASURES if c in df}) + " non-NaN (must all be 0)")
    L.append(f"\n## Timing: total_sec median {pd.to_numeric(df.total_sec, errors='coerce').median():.1f}, p95 {pd.to_numeric(df.total_sec, errors='coerce').quantile(.95):.1f}")
    rep = "\n".join(L); open(f"{P.X1_LOGS}/ASSEMBLE_LIGHT{a.log_suffix}.md", "w").write(rep + "\n"); print(rep)
    # cross-site (gate 4)
    ok = df[df.status.astype(str).str.startswith("ok")]
    C = [f"# CROSS_SITE_NEW_MEASURES ({len(ok)} ok rows). Gate 4: every per-site difference explained by physiology or the scoring convention, or STOP.", ""]
    for line in AC.per_site_table(ok, AC.NEW_MEASURES + ["spo2_pct_below_90_native", "hb_n_events_used", "lm_rows", "plm_n_series", "TST_min"]):
        c = line.split(":")[0][2:]
        C.append(line + (f"\n    note: {NOTES[c]}" if c in NOTES else ""))
    C.append("\nOwner action: confirm each note or STOP the measure (D6 for plm_index).")
    open(f"{P.X1_LOGS}/CROSS_SITE_NEW_MEASURES{a.log_suffix}.md", "w").write("\n".join(C) + "\n"); print("\n".join(C))
    df.to_parquet(a.out, index=False)
    idx = f"{P.X1_LOGS}/PROVENANCE_INDEX.tsv"
    inputs = [a.jsonl, a.manifest, P.REEXTRACT_FINAL, P.T90_FINAL_V7] + ([P.HB_DEFINITION] if hb_sha else [])
    for o in (a.out, f"{P.X1_LOGS}/ASSEMBLE_LIGHT{a.log_suffix}.md", f"{P.X1_LOGS}/CROSS_SITE_NEW_MEASURES{a.log_suffix}.md", f"{P.X1_DATA}/light_failures{a.log_suffix}.csv"):
        write_sidecar(o, inputs, __file__, "023_assemble_light", index_path=idx)
    fail_inv = [n for n, p, k, _ in inv if p != k]
    print(f"\nASSEMBLE_LIGHT: positive-control mismatches {n_bad}; invariant failures {fail_inv}; wrote {a.out}")
    return 0 if (n_bad == 0 and not fail_inv) else 1


if __name__ == "__main__":
    sys.exit(main())
