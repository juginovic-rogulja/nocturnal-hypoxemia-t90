"""
Step 012 gate: the v8 extractor with t_end=None on the 12 probe nights against (1) the v7 raw extractor output (the
full_s0/full_s1 jsonl, the pure extract_night record) and (2) the frozen reextract_FINAL.parquet (which carries the v7
post-processing: AHI construction override, cells set missing by rule, siteZ copies, oxygen-pass columns).
Every column present on both sides is compared; a value matches when |new - old| <= 1e-6 or both are missing.
Output: logs/probe12_reproduction.csv (one row per night x column with a mismatch, classified) and a summary on stdout.
usage: compare_to_reextract_final.py <probe12.jsonl> [--raw-jsonl a.jsonl b.jsonl ...]
"""
from __future__ import annotations
import argparse, json, math, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import x1_paths as P
from provenance import write_sidecar

TOL = 1e-6
V7_POST = ("_siteZ", "_oxypass", "_aprilstyle")   # v7 post-processing suffixes, never produced by the extractor
V7_POST_COLS = {"spo2_nadir_corrected", "file_category", "comparable_file", "subgroup", "oxygen_pass_status", "spo2_frac_integer",
                "spo2_frac_exact_90", "spo2_frac_exact_89", "spo2_pct_below_88_native", "spo2_nadir_native", "spo2_p1_native",
                "download_sec", "file_mb", "extract_wall_sec", "total_sec", "runtime_sec", "source", "traceback", "stem", "site_id",
                "_encoding_used", "BDSPPatientID", "status"}
NEW_V8 = ("_v8_", "hb_", "hypoxic_burden", "plm_", "lm_", "limb_code_hist", "limb_event_map", "limb_dur_", "resp_code_hist",
          "resp_event_map", "_sleep", "sol_min_nm", "rem_latency_min_nm", "spo2_native_valid_min", "diag_prepap_", "diag_light_index")


def load_jsonl(paths, ids):
    rows = []
    for p in paths:
        with open(p) as f:
            for line in f:
                try:
                    d = json.loads(line)
                except Exception:
                    continue
                if int(d.get("BDSPPatientID", -1)) in ids:
                    rows.append(d)
    return pd.DataFrame(rows)


def compare(new: pd.DataFrame, ref: pd.DataFrame, label: str):
    out = []
    cols = [c for c in new.columns if c in ref.columns and c not in V7_POST_COLS and not c.endswith(V7_POST)
            and not any(c.startswith(s) or s in c for s in NEW_V8)]
    n_cells = 0; n_match = 0
    for pid in new.BDSPPatientID.astype(int):
        a = new[new.BDSPPatientID.astype(int) == pid].iloc[0]; b = ref[ref.BDSPPatientID.astype(int) == pid]
        if b.empty:
            out.append({"ref": label, "BDSPPatientID": pid, "column": "<night>", "new": None, "old": None, "class": "night_missing_in_ref"}); continue
        b = b.iloc[0]
        for c in cols:
            x, y = a[c], b[c]
            n_cells += 1
            xf = _tofloat(x); yf = _tofloat(y)
            if xf is not None and yf is not None:
                if (math.isnan(xf) and math.isnan(yf)) or (abs(xf - yf) <= TOL) or (math.isinf(xf) and xf == yf):
                    n_match += 1; continue
                if math.isnan(yf) and not math.isnan(xf):
                    cls = "old_missing_new_finite"
                elif math.isnan(xf) and not math.isnan(yf):
                    cls = "new_missing_old_finite"
                else:
                    cls = "both_finite_differ"
                out.append({"ref": label, "BDSPPatientID": pid, "column": c, "new": xf, "old": yf, "absdiff": abs(xf - yf) if not (math.isnan(xf) or math.isnan(yf)) else np.nan, "class": cls})
            else:
                if str(x) == str(y) or (x is None and (y is None or (isinstance(y, float) and math.isnan(y)))):
                    n_match += 1; continue
                out.append({"ref": label, "BDSPPatientID": pid, "column": c, "new": str(x)[:80], "old": str(y)[:80], "absdiff": np.nan, "class": "text_differs"})
    return pd.DataFrame(out), n_cells, n_match, len(cols)


def _tofloat(v):
    if v is None:
        return float("nan")
    if isinstance(v, (bool, np.bool_)):
        return float(v)
    if isinstance(v, (int, float, np.integer, np.floating)):
        return float(v)
    if isinstance(v, str):
        try:
            return float(v)
        except Exception:
            return None
    return None


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("jsonl"); ap.add_argument("--raw-jsonl", nargs="*", default=[])
    ap.add_argument("--out", default=f"{P.X1_LOGS}/probe12_reproduction.csv")
    a = ap.parse_args()
    P.assert_not_evicted(P.REEXTRACT_FINAL, a.jsonl, *a.raw_jsonl)
    new = pd.DataFrame([json.loads(l) for l in open(a.jsonl)])
    ids = set(new.BDSPPatientID.astype(int))
    print(f"new: {len(new)} nights, status: {new.status.astype(str).str[:30].value_counts().to_dict()}")
    fin = pd.read_parquet(P.REEXTRACT_FINAL); fin = fin[fin.BDSPPatientID.isin(ids)]
    frames = []
    if a.raw_jsonl:
        raw = load_jsonl(a.raw_jsonl, ids)
        print(f"raw v7 jsonl: {len(raw)} of {len(ids)} nights found")
        d, nc, nm, ncol = compare(new, raw, "v7_raw_jsonl"); frames.append(d)
        print(f"[v7_raw_jsonl] {ncol} columns compared, {nc} cells, {nm} within {TOL} ({100*nm/max(nc,1):.3f}%), {len(d)} mismatches")
        if len(d):
            print(d["class"].value_counts().to_dict())
            print(d.groupby("column").size().sort_values(ascending=False).head(25).to_string())
    d2, nc2, nm2, ncol2 = compare(new, fin, "reextract_FINAL"); frames.append(d2)
    print(f"[reextract_FINAL] {ncol2} columns compared, {nc2} cells, {nm2} within {TOL} ({100*nm2/max(nc2,1):.3f}%), {len(d2)} mismatches")
    if len(d2):
        print(d2["class"].value_counts().to_dict())
        print(d2.groupby(["class", "column"]).size().sort_values(ascending=False).head(40).to_string())
    allm = pd.concat(frames) if frames else pd.DataFrame()
    allm.to_csv(a.out, index=False)
    write_sidecar(a.out, [a.jsonl, P.REEXTRACT_FINAL] + a.raw_jsonl, __file__, "012_probe12_reproduction", index_path=f"{P.X1_LOGS}/PROVENANCE_INDEX.tsv")
    print("wrote", a.out)


if __name__ == "__main__":
    main()
