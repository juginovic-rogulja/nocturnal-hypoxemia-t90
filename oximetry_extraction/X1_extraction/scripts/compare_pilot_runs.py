"""Compare two light-pilot runs of the same nights (default 100/100 windows vs the frozen definition): every column except the
HB family must be identical (S3 re-read reproducibility); the HB family is summarised per site, old vs new (decision D5 effect).
usage: compare_pilot_runs.py <a.jsonl> <b.jsonl> [--out logs/PILOT_RUNS_COMPARE.md]"""
from __future__ import annotations
import argparse, json, math, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import x1_paths as P
from provenance import write_sidecar

HB = ("hypoxic_burden", "hb_", "runtime_sec", "total_sec", "download_sec", "source", "diag_upsampled_fft", "diag_light_index")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("a"); ap.add_argument("b"); ap.add_argument("--out", default=f"{P.X1_LOGS}/PILOT_RUNS_COMPARE.md")
    x = ap.parse_args()
    A = pd.DataFrame([json.loads(l) for l in open(x.a)]).drop_duplicates("BDSPPatientID", keep="last").set_index("BDSPPatientID").sort_index()
    B = pd.DataFrame([json.loads(l) for l in open(x.b)]).drop_duplicates("BDSPPatientID", keep="last").set_index("BDSPPatientID").sort_index()
    ids = A.index.intersection(B.index)
    cols = [c for c in A.columns if c in B.columns and not any(c.startswith(h) or c == h for h in HB)]
    diff = []
    for c in cols:
        a = A.loc[ids, c]; b = B.loc[ids, c]
        an = pd.to_numeric(a, errors="coerce"); bn = pd.to_numeric(b, errors="coerce")
        if an.notna().sum() or bn.notna().sum():
            eq = (an == bn) | (an.isna() & bn.isna())
            if a.dtype == object and an.isna().all():
                eq = (a.astype(str) == b.astype(str))
        else:
            eq = (a.astype(str) == b.astype(str))
        n_bad = int((~eq).sum())
        if n_bad:
            diff.append((c, n_bad))
    L = [f"# PILOT_RUNS_COMPARE: {os.path.basename(x.a)} vs {os.path.basename(x.b)} on {len(ids)} shared nights",
         f"- non-HB columns compared: {len(cols)}; columns with any difference: {len(diff)} {diff[:20]}"]
    for site in ("I0002", "I0006"):
        m = A.loc[ids, "site_id"] == site
        ha = pd.to_numeric(A.loc[ids][m].hypoxic_burden, errors="coerce"); hb = pd.to_numeric(B.loc[ids][m].hypoxic_burden, errors="coerce")
        ok = ha.notna() & hb.notna()
        rho = ha[ok].corr(hb[ok], method="spearman") if ok.sum() > 3 else np.nan
        L.append(f"- HB {site}: A median {ha.median():.2f} IQR [{ha.quantile(.25):.2f}, {ha.quantile(.75):.2f}] | B median {hb.median():.2f} IQR [{hb.quantile(.25):.2f}, {hb.quantile(.75):.2f}] | "
                 f"Spearman A~B {rho:.4f} | median B/A {float((hb[ok] / ha[ok].replace(0, np.nan)).median()):.3f} | source A {A.loc[ids][m].hb_window_source.iloc[0][:40]} | B {B.loc[ids][m].hb_window_source.iloc[0][:40]}")
    rep = "\n".join(L); open(x.out, "w").write(rep + "\n"); print(rep)
    write_sidecar(x.out, [x.a, x.b], __file__, "020b_compare_pilot_runs", index_path=f"{P.X1_LOGS}/PROVENANCE_INDEX.tsv")
    return 0 if not diff else 1


if __name__ == "__main__":
    sys.exit(main())
