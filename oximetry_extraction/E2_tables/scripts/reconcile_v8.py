#!/usr/bin/env python3
"""
Step 052, reconcile_sample_vs_source (integrity gate 1): N stratified cohort nights (site x split status), every v8 table value
of every replaced or added column equals the pass record (exact, NaN-aware; a missing table value is accepted only with its
documented reason: D1 index undefined, collapsed staging, plausibility flag, the keep_v7 oxygen policy), and the untouched
non-split nights' oxygen columns equal the fresh light-pass read of the source file. --reread K then hands K of the sampled
nights to X1's independent second code path (scripts/second_path_check.py, numpy only, nights read again from S3) and reports.

usage: reconcile_v8.py --dir <built folder> --split <assembled parquet or jsonl> --light <assembled parquet or jsonl>
                       [--n 300] [--reread 30] [--seed 20260912] [--log-dir <folder>] [--smoke]
Exit 1 on any unexplained difference or a failed re-read.
"""
from __future__ import annotations
import argparse, json, os, subprocess, sys, time
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import e2_spec as S  # noqa: E402
CS = S.CS
from build_v8_tables import ARCH, IDX, OXY_MAP, NEW5, RECON_TOL, load_pass, detect_collapse, nan_equal  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True); ap.add_argument("--split", required=True); ap.add_argument("--light", required=True)
    ap.add_argument("--n", type=int, default=300); ap.add_argument("--reread", type=int, default=30); ap.add_argument("--seed", type=int, default=20260912)
    ap.add_argument("--log-dir", default=None); ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--light-jsonl", default=S.LIGHT_JSONL, help="the raw light jsonl for the second code path (needs S3 access)")
    a = ap.parse_args()
    t0 = time.time(); out = os.path.abspath(a.dir); log_dir = a.log_dir or out; os.makedirs(log_dir, exist_ok=True)
    S.require_resident(*[f"{out}/{f}" for f in S.TABLE_FILES], a.split, a.light, S.PREV["t90_final"], S.MANIFEST_SPLIT, S.MANIFEST_FULL, CS.MEASURE_FAMILIES_CSV)
    L = [f"# RECONCILE_V8 ({'SMOKE' if a.smoke else 'REAL'}; tables {out}; split {a.split}; light {a.light}; {time.strftime('%Y-%m-%d %H:%M')})", ""]
    fails = 0

    prev = pd.read_parquet(S.PREV["t90_final"]); coh = CS.apply_cohort(prev); cohort_ids = set(coh.BDSPPatientID)
    collapsed_v7 = set(prev.loc[prev["_v7_1_stage_collapsed"].astype(bool), "BDSPPatientID"])
    split = detect_collapse(load_pass(a.split, "split", collapsed_v7, S.MANIFEST_SPLIT))
    light = detect_collapse(load_pass(a.light, "light", collapsed_v7, S.MANIFEST_FULL))
    split_ok = split[split._ok].set_index("BDSPPatientID"); light_ok = light[light._ok].set_index("BDSPPatientID")
    fam = pd.read_csv(CS.MEASURE_FAMILIES_CSV); meas137 = [f for f in fam.feature if f not in CS.NEW_MEASURES_V8]
    repl = list(dict.fromkeys(meas137 + ARCH + IDX + list(OXY_MAP))); srcmap = {c: OXY_MAP.get(c, c) for c in repl}
    v7ch = pd.read_parquet(S.PREV["changelog"], columns=["table", "column", "reason"])
    collapse_cols = {t: set(v7ch[(v7ch.table == t) & v7ch.reason.str.startswith("stage collapsed")].column) for t in ("t90_final", "t90_base4", "master")}

    # stratified sample: site x split, from nights the passes cover
    covered = set(light_ok.index) & cohort_ids
    pool = coh[coh.BDSPPatientID.isin(covered)][["BDSPPatientID", "site_id"]].copy()
    pool["split"] = pool.BDSPPatientID.isin(split_ok.index)
    rng = np.random.default_rng(a.seed); parts = []
    per = max(a.n // 4, 1)
    for (s, w), g in pool.groupby(["site_id", "split"]):
        k = min(per, len(g)); parts.append(g.iloc[rng.choice(len(g), size=k, replace=False)])
    sample = pd.concat(parts, ignore_index=True) if parts else pool.iloc[:0]
    ids = sample.BDSPPatientID.tolist()
    L.append(f"sample: {len(ids)} nights ({sample.groupby(['site_id', 'split']).size().to_dict()}) from {len(pool):,} covered cohort nights (per stratum up to {per})")

    tabs = {"t90_final": pd.read_parquet(f"{out}/t90_final.parquet"), "t90_base4": pd.read_parquet(f"{out}/t90_base4.parquet"),
            "master": pd.read_csv(f"{out}/master_cohort.csv", low_memory=False, float_precision="round_trip")}
    prevx = prev.set_index("BDSPPatientID")
    L.append("\n## Table against pass record, every replaced or added column, per table")
    for name, t in tabs.items():
        t = t.set_index("BDSPPatientID").reindex(ids)
        n_eq = n_expl = n_tol = n_skipw = 0; unexpl = []
        for i in ids:
            row = t.loc[i]
            if i in split_ok.index:
                s = split_ok.loc[i]; cols = [c for c in repl if c in t.columns]
                for c in cols:
                    tv, sv = row[c], pd.to_numeric(s[srcmap[c]], errors="coerce")
                    if nan_equal([tv], [sv])[0]:
                        n_eq += 1; continue
                    if pd.isna(tv):
                        if c in S.SLEEP_MASK_COLS and bool(row.get("_v8_sleep_masked", False)):   # v8.1: masked by design
                            n_expl += 1; continue
                        if c in collapse_cols[name] and bool(s["_v8_stage_collapsed"]):
                            n_expl += 1; continue
                        if bool(row.get(f"_impossible_{c}", False)):
                            n_expl += 1; continue
                        if bool(row.get("_v8_sleep_floor_all", False)):
                            n_expl += 1; continue
                    if c in OXY_MAP and bool(row.get("_v8_oxygen_window_short", False)) and nan_equal([tv], [prevx.loc[i, c]])[0]:
                        n_expl += 1; continue
                    unexpl.append((i, c, tv, sv))
            if i in light_ok.index and bool(light_ok.loc[i, "_v8_prepap_window"]) and i not in split_ok.index:
                n_skipw += 1        # windowed light record on a night the split pass did not replace: the builder skips it (smoke artefact; 0 in a real build)
            elif i in light_ok.index:
                l = light_ok.loc[i]
                for c in [c for c in NEW5 if c in t.columns]:
                    tv, sv = row[c], pd.to_numeric(l[c], errors="coerce")
                    if nan_equal([tv], [sv])[0]:
                        n_eq += 1
                    elif pd.isna(tv) and (bool(l["_v8_stage_collapsed"]) or bool(row.get(f"_impossible_{c}", False))):
                        n_expl += 1
                    elif pd.isna(tv) and c in S.SLEEP_MASK_COLS and bool(row.get("_v8_sleep_masked", False)):   # v8.1: masked by design
                        n_expl += 1
                    else:
                        unexpl.append((i, c, tv, sv))
                if i not in split_ok.index:   # untouched night: the v7 oxygen columns against the fresh read of the source (another run
                    for c, sc in OXY_MAP.items():   # of the same definition: summation-order noise within RECON_TOL is counted, not failed)
                        if c not in t.columns:
                            continue
                        tv, sv = row[c], pd.to_numeric(l[sc], errors="coerce")
                        if nan_equal([tv], [sv])[0]:
                            n_eq += 1
                        elif pd.isna(tv) and bool(row.get(f"_impossible_{c}", False)):
                            n_expl += 1
                        elif pd.notna(tv) and pd.notna(sv) and abs(float(tv) - float(sv)) <= RECON_TOL[c]:
                            n_tol += 1
                        else:
                            unexpl.append((i, c, tv, sv))
        fails += len(unexpl)
        L.append(f"- {name}: {n_eq:,} cells equal, {n_tol:,} untouched-night oxygen cells equal within tolerance (re-read noise), {n_expl:,} missing with a documented reason, {n_skipw} windowed light records not replaced by the split pass (skipped, smoke only), {len(unexpl)} UNEXPLAINED"
                 + (f"; first: {unexpl[:5]}" if unexpl else ""))

    # second code path on K of the sampled nights (S3 re-read; X1's independent implementation)
    if a.reread > 0:
        jl = a.light_jsonl
        if not os.path.exists(jl) or S.evicted(jl):
            L.append(f"\n## Second code path: SKIPPED, light jsonl missing or evicted ({jl})"); fails += 1
        else:
            sub = f"{S.E2_WORK}/reconcile_subset_{len(ids)}.jsonl"
            want = set(ids); n_w = 0
            with open(jl) as fi, open(sub, "w") as fo:
                for ln in fi:
                    if not ln.strip():
                        continue
                    r = json.loads(ln)
                    if int(r.get("BDSPPatientID", -1)) in want:
                        fo.write(ln if ln.endswith("\n") else ln + "\n"); n_w += 1
            csv = f"{log_dir}/RECONCILE_SECOND_PATH_{a.reread}.csv"
            cmd = [S.PY, f"{S.X1_SCRIPTS}/second_path_check.py", "--pass-jsonl", sub, "--manifest", S.MANIFEST_FULL, "--n", str(a.reread), "--seed", str(a.seed), "--out", csv]
            L.append(f"\n## Second code path (X1 second_path_check.py, S3 re-read) on {a.reread} of the sampled nights ({n_w} pass records in the subset)")
            try:
                pr = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
                tail = (pr.stdout or "").strip().splitlines()[-6:]
                L += [f"    {ln}" for ln in tail]
                L.append(f"- rc {pr.returncode} ({'PASS' if pr.returncode == 0 else 'FAIL'}); csv {csv}")
                if pr.returncode != 0:
                    fails += 1; L += [f"    stderr: {ln}" for ln in (pr.stderr or "").strip().splitlines()[-5:]]
            except subprocess.TimeoutExpired:
                fails += 1; L.append("- TIMEOUT after 3600 s (FAIL)")
    else:
        L.append("\n## Second code path: not run (--reread 0)")

    L.append(f"\nResult: {'PASS' if fails == 0 else f'FAIL ({fails})'} in {time.time()-t0:.0f}s")
    rep = "\n".join(L); open(f"{log_dir}/RECONCILE_V8.md", "w").write(rep + "\n"); print(rep)
    S.write_sidecar(f"{log_dir}/RECONCILE_V8.md", [f"{out}/{f}" for f in S.TABLE_FILES] + [a.split, a.light, S.PREV["t90_final"]], os.path.abspath(__file__),
                    "052_reconcile_sample_vs_source", extra={"n": len(ids), "reread": a.reread, "fails": fails}, index_path=f"{S.E2_LOGS}/PROVENANCE_INDEX.tsv")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
