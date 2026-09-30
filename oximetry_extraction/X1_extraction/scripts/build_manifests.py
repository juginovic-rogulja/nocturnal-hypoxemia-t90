"""
Step 010, build_manifests_v8. Inputs (v7, read only): cpap_t90_by_stage (design A_split, cpap_start_min, BIDSFolder,
SessionID, SiteID, pre_any_t90, pre_sleep_min), t90_final through apply_cohort, manifest_full_19173.csv, reextract_FINAL
(spo2_native_fs). Outputs (X1/work + V8/smoke): manifest_split_3622.csv, manifest_full_19173_v8.csv, manifest_probe12.csv,
sample_300.csv, sample_300_split.csv, each with a provenance sidecar. Smoke checks inside: counts, unique ids, t_end_sec >= 2700.
The split file = the file the PAP switch-on was measured on (cpap table BIDSFolder + SessionID); where that session differs
from the v7 manifest's session the row is flagged _v8_session_differs_from_v7 (owner decision, default: the cpap table's file).
"""
from __future__ import annotations
import os, shutil, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import x1_paths as P
from provenance import write_sidecar

STEP = "010_build_manifests_v8"


def stratified_sample(df, strata_cols, n, seed):
    """Proportional allocation over the strata, at least one row per non-empty stratum, exact total n."""
    rng = np.random.default_rng(seed)
    g = df.groupby(strata_cols, dropna=False, sort=True)
    sizes = g.size()
    alloc = np.floor(sizes / sizes.sum() * n).astype(int)
    alloc[(alloc == 0) & (sizes > 0)] = 1
    # fix the total by largest remainders
    while alloc.sum() > n:
        k = (alloc - sizes / sizes.sum() * n).idxmax(); alloc[k] -= 1
    while alloc.sum() < n:
        rem = sizes / sizes.sum() * n - alloc
        rem[alloc >= sizes] = -np.inf
        k = rem.idxmax(); alloc[k] += 1
    parts = []
    for key, sub in g:
        k = alloc[key]
        if k > 0:
            parts.append(sub.sample(n=min(k, len(sub)), random_state=int(rng.integers(0, 2**31 - 1))))
    out = pd.concat(parts).sort_values("BDSPPatientID").reset_index(drop=True)
    assert len(out) == n, len(out)
    return out


def main():
    P.ensure_dirs()
    P.assert_not_evicted(P.CPAP_STAGE_V7, P.T90_FINAL_V7, P.MANIFEST_FULL_V7, P.REEXTRACT_FINAL)
    coh = P.cohort_frame(["site_id", "_v7_1_stage_collapsed"])
    ids = set(coh.BDSPPatientID.astype(int))
    m = pd.read_csv(P.MANIFEST_FULL_V7)
    assert len(m) == len(ids) and set(m.BDSPPatientID.astype(int)) == ids, "v7 manifest is not the cohort"
    c = pd.read_parquet(P.CPAP_STAGE_V7, columns=["BDSPPatientID", "SiteID", "BIDSFolder", "SessionID", "design",
                                                  "cpap_start_min", "pre_any_t90", "pre_sleep_min", "pap_src", "status", "fs_spo2"])
    s = c[(c.design == "A_split") & c.BDSPPatientID.isin(ids)].copy()
    assert not s.BDSPPatientID.duplicated().any(), "duplicate split ids"
    r = pd.read_parquet(P.REEXTRACT_FINAL, columns=["BDSPPatientID", "spo2_native_fs", "recording_dur_min"])
    s = s.merge(m.rename(columns={"SessionID": "SessionID_v7", "bids_col": "bids_col_v7"}), on="BDSPPatientID", how="left")
    s = s.merge(r, on="BDSPPatientID", how="left")
    s = s.merge(coh[["BDSPPatientID", "_v7_1_stage_collapsed"]], on="BDSPPatientID", how="left")
    assert (s.SiteID == s.site_id).all() and (s.BIDSFolder == s.bids_col_v7).all()
    s["_v8_session_differs_from_v7"] = (s.SessionID.astype(int) != s.SessionID_v7.astype(int))
    s["t_end_sec"] = (s.cpap_start_min * 60.0).round(3)
    split = pd.DataFrame({"site_id": s.site_id, "BDSPPatientID": s.BDSPPatientID.astype(int), "bids_col": s.BIDSFolder,
                          "SessionID": s.SessionID.astype(int), "_encoding_used": s._encoding_used,
                          "t_end_sec": s.t_end_sec, "spo2_native_fs": s.spo2_native_fs, "fs_spo2_cpap_table": s.fs_spo2,
                          "cpap_start_min": s.cpap_start_min, "pre_any_t90_v7": s.pre_any_t90, "pre_sleep_min_v7": s.pre_sleep_min,
                          "recording_dur_min_v7": s.recording_dur_min, "SessionID_v7": s.SessionID_v7.astype(int),
                          "_v8_session_differs_from_v7": s._v8_session_differs_from_v7,
                          "_v7_1_stage_collapsed": s._v7_1_stage_collapsed}).sort_values("BDSPPatientID").reset_index(drop=True)
    # smoke checks (counts from the plan; the cohort size itself comes from cohort_spec)
    assert len(split) == P.COHORT_N_EXPECTED_SPLIT, f"split rows {len(split)}"
    assert (split.t_end_sec >= 2700).all(), "a split row has t_end_sec < 2700"
    assert split.t_end_sec.notna().all()
    assert ((split.recording_dur_min_v7 * 60 > split.t_end_sec) | split._v8_session_differs_from_v7).all(), "t_end beyond the v7 recording"
    split.to_csv(P.MANIFEST_SPLIT, index=False)
    # full 19,173 manifest: v7 columns + t_end_sec (empty on non-split nights); split rows use the cpap table's session
    full = m.copy()
    full["t_end_sec"] = np.nan
    full = full.merge(split[["BDSPPatientID", "SessionID", "t_end_sec"]].rename(columns={"SessionID": "SessionID_split", "t_end_sec": "t_end_split"}),
                      on="BDSPPatientID", how="left")
    is_split = full.t_end_split.notna()
    full.loc[is_split, "t_end_sec"] = full.loc[is_split, "t_end_split"]
    full.loc[is_split, "SessionID"] = full.loc[is_split, "SessionID_split"].astype(int)
    full["is_split"] = is_split.astype(int)
    full = full.drop(columns=["SessionID_split", "t_end_split"]).merge(coh[["BDSPPatientID", "_v7_1_stage_collapsed"]], on="BDSPPatientID")
    assert len(full) == len(ids) and full.is_split.sum() == len(split)
    full = full.sort_values("BDSPPatientID").reset_index(drop=True)
    full.to_csv(P.MANIFEST_FULL_V8, index=False)
    # probe 12
    probe = full[full.BDSPPatientID.isin(P.PROBE12)].copy()
    probe["order"] = probe.BDSPPatientID.map({b: i for i, b in enumerate(P.PROBE12)})
    probe = probe.sort_values("order").drop(columns="order")
    assert len(probe) == 12 and probe.is_split.sum() == 0, "probe nights must be non-split"
    probe["t_end_sec"] = ""
    probe.to_csv(P.MANIFEST_PROBE12, index=False)
    # samples
    samp = stratified_sample(full, ["site_id", "is_split", "_encoding_used"], 300, P.SEED)
    samp.to_csv(P.SAMPLE_300, index=False)
    samp_split = stratified_sample(split, ["site_id", "_encoding_used"], 300, P.SEED)
    samp_split.to_csv(P.SAMPLE_300_SPLIT, index=False)
    for src in (P.SAMPLE_300, P.SAMPLE_300_SPLIT):
        shutil.copy2(src, os.path.join(P.V8_SMOKE, os.path.basename(src)))
    inputs = [P.CPAP_STAGE_V7, P.T90_FINAL_V7, P.MANIFEST_FULL_V7, P.REEXTRACT_FINAL, f"{P.NUMBERS}/cohort_spec.py"]
    idx = f"{P.X1_LOGS}/PROVENANCE_INDEX.tsv"
    for out in (P.MANIFEST_SPLIT, P.MANIFEST_FULL_V8, P.MANIFEST_PROBE12, P.SAMPLE_300, P.SAMPLE_300_SPLIT,
                f"{P.V8_SMOKE}/sample_300.csv", f"{P.V8_SMOKE}/sample_300_split.csv"):
        write_sidecar(out, inputs, __file__, STEP, index_path=idx)
    print(f"split {len(split)} (I0002 {int((split.site_id=='I0002').sum())}, I0006 {int((split.site_id=='I0006').sum())}); "
          f"session differs from v7 on {int(split._v8_session_differs_from_v7.sum())} rows: "
          f"{split[split._v8_session_differs_from_v7].BDSPPatientID.tolist()}")
    print(f"full {len(full)} (split {int(full.is_split.sum())}); probe {len(probe)}; sample_300 {len(samp)} "
          f"(split {int(samp.is_split.sum())}); sample_300_split {len(samp_split)}")
    print("t_end_sec on split rows: min %.1f median %.1f max %.1f" % (split.t_end_sec.min(), split.t_end_sec.median(), split.t_end_sec.max()))
    print("sample_300 strata:", samp.groupby(["site_id", "is_split", "_encoding_used"]).size().to_dict())
    print("spo2_native_fs on split:", split.spo2_native_fs.value_counts().to_dict(),
          "fs disagreement cpap-table vs reextract:", int((split.spo2_native_fs != split.fs_spo2_cpap_table).sum()))


if __name__ == "__main__":
    main()
