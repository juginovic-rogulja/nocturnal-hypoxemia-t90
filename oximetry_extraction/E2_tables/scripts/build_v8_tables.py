#!/usr/bin/env python3
"""
Step 050, build_v8_tables: the four v8 frozen T90 tables from the v7 tables (PREV_DIR, read only) plus the three v8 inputs.

  --split    X1 step 015 output (split3622_assembled.parquet; a jsonl pass is accepted, X1's assemble rules are then applied here)
  --light    X1 step 023 output (light19173_assembled.parquet or the jsonl)
  --outcomes X2 step 041 output (outcomes_v8.parquet)
  --out      the output folder (the real build: cohort_spec DATA_DIR = data_frozen_v8_2026-09; a smoke build: anywhere else)
  --smoke    the passes may cover a subset of nights (pilot data); completeness is reported instead of asserted
  --ids CSV  restrict the v8 replacement to these BDSPPatientIDs (the psv smoke form)

What it does, per table (t90_final, t90_base4, master_cohort; cpap_t90_by_stage covariates only):
  1. split nights (decision 1): every night-level measure column (the 137 kept measures, the sleep architecture columns, the
     indices, the six oxygen measures) replaced by the pre-PAP window value from --split. D1 default: AHI and arousal index
     undefined under 60 min of pre-PAP sleep (X1 applied it; counted here), oxygen measures need 5 valid minutes (policy
     --short-oxygen-window: keep_v7 = keep the v7 whole-night oxygen values on those nights, flagged; nan = set missing, the
     cohort then shrinks), --sleep-floor-all = the D1 alternative (every sleep-gated measure missing under 60 min).
  2. collapsed staging (v7.1 rule, addendum 2026-09-08): a pass record whose sleep holds ONE stage code (all-REM current
     files) is impossible; composition and stage-gated features are set missing (the v7.1 column set, read from the v7
     CHANGELOG), the five new measures too (D11). Detected from stage_epoch_hist in the pass, not from the v7 flag alone.
  3. all cohort nights: hypoxic_burden, plm_index, sol_min, rem_latency_min and spo2_pct_below_90_sleep from --light
     (into every table; QC columns hb_n_events_used, lm_index, plm_n_series ... into t90_final only, so the ranking's
     scan of the master header sees exactly N_MEASURES candidates).
  4. plausibility bounds (v7 step 5) on every replaced or added value: outside -> missing + _impossible_<col> flag.
  5. outcomes: every *_first_date/_prevalent/_incident/_years column of --outcomes into t90_final and t90_base4 (added
     where absent), death columns and the v8 provenance columns; the stale obesity_hypovent_first_date rebuilt from the
     hospital-record cache with the v8 list and the verbatim v7 rule (flags kept, agreement reported).
  6. the 60 _siteZ columns dropped from master_cohort and t90_base4 (decision 3; the swept scripts never name them).
  7. cpap_t90_by_stage: whole-night AHI, TST_min, spo2_pct_below_90 refreshed from the v8 t90_final as v7 did.
  8. provenance columns (v7 ones kept): _v8_prepap_window, _v8_window_end_sec, _v8_window_end_min, _v8_untreated_sleep_min,
     _v8_session_differs_from_v7, _v8_index_undefined, _v8_oxygen_window_short, _v8_stage_collapsed, _v8_measure_source,
     _v8_ahi_source, _v8_prepap_source, _v8_new_measure_source, _v8_outcome_list_version, _v8_light_status, _v8_hb_definition_sha.
Outputs: the four tables, CHANGELOG.parquet (one row per changed cell), CHANGELOG.md, PROVENANCE.md, sha256.txt, a
provenance sidecar per table, then invariants_v8.run on the folder (INVARIANTS_V8.md, PER_SITE_ALL_COLUMNS_V8.csv); rc 1 on
any invariant failure. Every untouched cell is asserted byte-identical to v7 (the CHANGELOG is the only difference).
No result literal in any assert: every reference count is read from the v7 tables, the manifests or X2's event_counts at run time.
"""
from __future__ import annotations
import argparse, json, os, shutil, sys, time, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import e2_spec as S  # noqa: E402
CS = S.CS
sys.path.insert(0, S.X1_SCRIPTS)
import assemble_common as AC  # noqa: E402  (X1's jsonl loader, D1/D11 rules; used only when a jsonl is passed)

OXY_MAP = {"spo2_pct_below_90": "spo2_pct_below_90_native", "spo2_pct_below_88": "spo2_pct_below_88_native",
           "spo2_mean": "spo2_mean_native", "spo2_nadir_corrected": "spo2_p1_native",
           "odi3_total": "odi3_total_native", "odi4_total": "odi4_total_native"}
ARCH = ["TST_min", "sleep_efficiency_pct", "N1_pct", "N2_pct", "N3_pct", "REM_pct", "N1_min", "N2_min", "N3_min", "REM_min",
        "WASO_min", "wake_min", "recording_dur_min"]
IDX = ["AHI", "arousal_index"]
NEW4 = list(CS.NEW_MEASURES_V8)
SLEEP_T90 = "spo2_pct_below_90_sleep"
NEW5 = NEW4 + [SLEEP_T90]
QC_NEW = ["hb_n_events_used", "hb_w_pre", "hb_w_post", "hb_window_source", "lm_index", "plm_n_series", "spo2_sleep_valid_min",
          "spo2_native_valid_min", "n_scored_epochs"]
COMPOSITION = ["N1_pct", "N2_pct", "N3_pct", "REM_pct", "N1_min", "N2_min", "N3_min", "REM_min"]
# gate-1 tolerances for the light pass against the v7 oxygen columns on untouched nights (percent of recording, percent saturation,
# events per hour): summation-order noise sits at 1e-15 / 1e-5 / 1e-7; a definitional change would exceed these
RECON_TOL = {"spo2_pct_below_90": 1e-3, "spo2_pct_below_88": 1e-3, "spo2_mean": 1e-2, "spo2_nadir_corrected": 0.5, "odi3_total": 1e-3, "odi4_total": 1e-3}
T0 = time.time()
LOG = []


def P(msg):
    LOG.append(msg)
    print(f"[{time.time()-T0:6.0f}s] {msg}", flush=True)


class Changes:
    """One row per changed cell: table, BDSPPatientID, column, old, new (floats where numeric, strings otherwise), reason."""

    def __init__(self):
        self.parts = []
        self.dropped = {}

    @staticmethod
    def _diff(old, new):
        o, n = pd.Series(np.asarray(old), dtype=object), pd.Series(np.asarray(new), dtype=object)
        try:
            eq = (o == n) | (o.isna() & n.isna())
        except Exception:
            eq = o.astype(str) == n.astype(str)
        return ~eq.to_numpy(dtype=bool)

    def add(self, table, ids, col, old, new, reason):
        d = self._diff(old, new)
        if not d.any():
            return 0
        o, n = np.asarray(old, dtype=object)[d], np.asarray(new, dtype=object)[d]
        num = pd.api.types.is_numeric_dtype(np.asarray(old)) and pd.api.types.is_numeric_dtype(np.asarray(new))
        if num:
            of, nf, os_, ns_ = pd.to_numeric(o, errors="coerce"), pd.to_numeric(n, errors="coerce"), None, None
        else:
            of, nf = np.nan, np.nan
            os_ = [None if pd.isna(v) else str(v) for v in o]
            ns_ = [None if pd.isna(v) else str(v) for v in n]
        self.parts.append(pd.DataFrame({"table": table, "BDSPPatientID": np.asarray(ids)[d].astype("int64"), "column": col,
                                        "old": of, "new": nf, "old_str": os_, "new_str": ns_, "reason": reason}))
        return int(d.sum())

    def frame(self):
        if not self.parts:
            return pd.DataFrame(columns=["table", "BDSPPatientID", "column", "old", "new", "old_str", "new_str", "reason"])
        f = pd.concat(self.parts, ignore_index=True)
        f["old"] = pd.to_numeric(f["old"], errors="coerce").astype(float); f["new"] = pd.to_numeric(f["new"], errors="coerce").astype(float)
        return f


CH = Changes()


def _b(x):
    """object/bool series with NaN -> numpy bool array (NaN = False), without pandas' downcasting warning"""
    x = pd.Series(x)
    return x.where(x.notna(), False).astype(bool).to_numpy()


def ensure_float(t, col):
    if col in t.columns and not pd.api.types.is_float_dtype(t[col]):
        if pd.api.types.is_bool_dtype(t[col]) or pd.api.types.is_integer_dtype(t[col]) or t[col].dtype == object:
            t[col] = pd.to_numeric(t[col], errors="coerce").astype(float)


def nan_equal(a, b):
    a, b = pd.Series(np.asarray(a), dtype=object), pd.Series(np.asarray(b), dtype=object)
    try:
        return ((a == b) | (a.isna() & b.isna())).to_numpy(dtype=bool)
    except Exception:
        return (a.astype(str) == b.astype(str)).to_numpy(dtype=bool)


# ---------------------------------------------------------------------------------------------------------------- inputs
def verify_prev_sha(paths):
    listed = {}
    for ln in open(S.PREV["sha256"]).read().strip().splitlines():
        sha, name = ln.split()
        listed[name] = sha
    bad = []
    for k in ("t90_final", "t90_base4", "master", "cpap"):
        name = os.path.basename(paths[k]); got = S.sha256(paths[k])
        if listed.get(name) != got:
            bad.append(f"{name}: sha256.txt {listed.get(name)} vs file {got}")
    return listed, bad


def load_pass(path, kind, collapsed_v7, manifest):
    """An assembled parquet as is; a jsonl through X1's loader + D1/D11 rules + the manifest columns (assemble_* logic)."""
    if path.endswith(".parquet"):
        df = pd.read_parquet(path)
        src_tag = f"{os.path.basename(path)} sha256:{S.sha256(path)[:16]}"
        if "status" not in df.columns and "_v8_light_status" in df.columns:
            df["status"] = df["_v8_light_status"]
    else:
        df = AC.load_jsonl(path)
        man = pd.read_csv(manifest)
        keep = [c for c in ("BDSPPatientID", "site_id", "SessionID", "_encoding_used", "t_end_sec", "cpap_start_min",
                            "_v8_session_differs_from_v7") if c in man.columns]
        man = man[man.BDSPPatientID.isin(df.BDSPPatientID)]
        df = man[keep].merge(df.drop(columns=[c for c in ("site_id", "_encoding_used") if c in df.columns]), on="BDSPPatientID", how="left")
        df["status"] = df.status.fillna("missing_record")
        df = AC.apply_rules(df, collapsed_v7, S.MIN_TST_FOR_INDEX, S.MIN_OXYGEN_VALID_MIN)
        src_tag = f"{os.path.basename(path)} sha256:{S.sha256(path)[:16]}"
        df["_v8_window_end_min"] = pd.to_numeric(df.get("_v8_window_end_min"), errors="coerce")
    for c in ("_v8_prepap_window", "_v8_index_undefined", "_v8_oxygen_window_short", "_v7_1_stage_collapsed", "_v8_session_differs_from_v7"):
        df[c] = _b(df[c]) if c in df.columns else False
    df["_src_tag"] = src_tag
    df["_ok"] = df.status.astype(str).str.startswith("ok")
    assert not df.BDSPPatientID.duplicated().any(), f"{kind}: duplicated BDSPPatientID"
    return df


def stage_codes_present(hist):
    if hist is None or (isinstance(hist, float) and np.isnan(hist)):
        return np.nan, None
    if isinstance(hist, str):
        try:
            hist = json.loads(hist)
        except Exception:
            return np.nan, None
    if not isinstance(hist, dict):
        return np.nan, None
    present = [int(k) for k, v in hist.items() if int(k) in S.STAGE_CODES_SLEEP and v and v > 0]
    return len(present), (present[0] if len(present) == 1 else None)


def detect_collapse(df):
    """Gate 2 of the 2026-09-08 addendum on the pass record: sleep present, one stage code. A windowed record under 60 min
    of sleep whose single code is not REM is a short window, flagged separately and kept."""
    n, code = zip(*[stage_codes_present(h) for h in df.get("stage_epoch_hist", pd.Series([None] * len(df)))]) if len(df) else ((), ())
    df["_v8_stage_codes_present"] = pd.to_numeric(pd.Series(n, index=df.index), errors="coerce")
    single_code = pd.Series(code, index=df.index, dtype=object)
    tst = pd.to_numeric(df.get("TST_min"), errors="coerce")
    one = df._ok & (tst > 0) & (df._v8_stage_codes_present == 1)
    windowed = df._v8_prepap_window
    short_win = one & windowed & (tst < S.COLLAPSE_MIN_TST_WINDOWED) & (single_code != 4)
    df["_v8_single_stage_short_window"] = short_win.to_numpy(dtype=bool)
    df["_v8_stage_collapsed"] = (one & ~short_win).to_numpy(dtype=bool)
    return df


# ---------------------------------------------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--split", default=S.first_existing(S.SPLIT_ASSEMBLED_CANDIDATES) or S.SPLIT_ASSEMBLED_CANDIDATES[0])
    ap.add_argument("--light", default=S.first_existing(S.LIGHT_ASSEMBLED_CANDIDATES) or S.LIGHT_ASSEMBLED_CANDIDATES[0])
    ap.add_argument("--outcomes", default=S.OUTCOMES_V8)
    ap.add_argument("--ids", default="", help="CSV with BDSPPatientID: restrict the v8 replacement to these nights")
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--sleep-floor-all", action="store_true", help="D1 alternative: every sleep-gated measure missing under 60 min of pre-PAP sleep")
    ap.add_argument("--short-oxygen-window", choices=("keep_v7", "nan"), default="keep_v7")
    ap.add_argument("--no-ohs-rebuild", action="store_true")
    ap.add_argument("--outcomes-into-master", action="store_true")
    ap.add_argument("--allow-recon-mismatch", action="store_true")
    ap.add_argument("--skip-invariants", action="store_true")
    a = ap.parse_args()

    out = os.path.abspath(a.out)
    if os.path.abspath(S.PREV_DIR) == out:
        sys.exit("refusing to write into PREV_DIR")
    if a.smoke and out == os.path.abspath(S.REAL_OUT):
        sys.exit("a smoke build must not write into the real DATA_DIR")
    if not a.smoke and out != os.path.abspath(S.REAL_OUT):
        P(f"NOTE: real mode writes to {out}, not cohort_spec DATA_DIR {S.REAL_OUT}")
    if os.path.isdir(out) and os.listdir(out) and not a.force:
        sys.exit(f"{out} exists and is not empty, pass --force")
    if CS.smoke_env_active():
        sys.exit(f"cohort_spec smoke variables are set ({CS.smoke_env_active()}); unset them, the builder uses the real cohort rule")
    os.makedirs(out, exist_ok=True)
    for f in os.listdir(out):
        if f.startswith("INVARIANTS_FAILED"):
            os.remove(os.path.join(out, f))
    P(f"build_v8_tables --out {out} {'SMOKE' if a.smoke else 'REAL'}; split {a.split}; light {a.light}; outcomes {a.outcomes}")

    inputs = [S.PREV[k] for k in ("t90_final", "t90_base4", "master", "cpap", "sha256", "changelog")] + \
             [a.split, a.light, a.outcomes, S.MANIFEST_SPLIT, S.MANIFEST_FULL, S.EVENT_COUNTS_V8, S.DEFS_V8, CS.MEASURE_FAMILIES_CSV]
    if not a.no_ohs_rebuild:
        inputs.append(S.OMOP_CACHE)
    S.require_resident(*inputs)
    bounds_path = S.pick(S.BOUNDS_CSV)
    bounds_src = "csv" if (os.path.exists(bounds_path) and not S.evicted(bounds_path)) else "v7 CHANGELOG reasons (bounds csv evicted)"
    listed, bad = verify_prev_sha(S.PREV)
    if bad:
        sys.exit("PREV_DIR tables do not match PREV_DIR/sha256.txt:\n  " + "\n  ".join(bad))
    P(f"v7 tables verified against {S.PREV['sha256']} (t90_final {listed['t90_final.parquet'][:16]})")

    # ------------------------------------------------------------------ v7 tables (read only)
    fz = pd.read_parquet(S.PREV["t90_final"]); b4 = pd.read_parquet(S.PREV["t90_base4"])
    # float_precision="round_trip": pandas' default parser is off by one ulp on some values; the v7 text must survive the rewrite
    ma = pd.read_csv(S.PREV["master"], low_memory=False, float_precision="round_trip"); cp = pd.read_parquet(S.PREV["cpap"])
    prev_fz, prev_b4, prev_ma, prev_cp = fz.copy(), b4.copy(), ma.copy(), cp.copy()
    n_rows_v7 = len(prev_fz)
    assert len(prev_b4) == n_rows_v7 and len(prev_ma) == n_rows_v7, "v7 tables differ in row count"
    assert set(prev_fz.BDSPPatientID) == set(prev_ma.BDSPPatientID) == set(prev_b4.BDSPPatientID)
    coh7 = CS.apply_cohort(prev_fz)
    cohort_ids = set(coh7.BDSPPatientID)
    collapsed_v7 = set(prev_fz.loc[prev_fz["_v7_1_stage_collapsed"].astype(bool), "BDSPPatientID"])
    v7ch = pd.read_parquet(S.PREV["changelog"], columns=["table", "column", "reason"])
    # the v7.1 patch labelled the master rows "master_cohort"; the v7.0 builder used "master"
    _lab = {"t90_final": ("t90_final",), "t90_base4": ("t90_base4",), "master": ("master", "master_cohort")}
    collapse_cols = {t: sorted(set(v7ch[v7ch.table.isin(_lab[t]) & v7ch.reason.str.startswith("stage collapsed")].column))
                     for t in ("t90_final", "t90_base4", "master")}
    man_split = pd.read_csv(S.MANIFEST_SPLIT)
    split_ids_all = set(man_split.BDSPPatientID)
    P(f"v7 read: rows {n_rows_v7:,}, cohort {len(cohort_ids):,}, split manifest {len(split_ids_all):,} nights "
      f"({len(split_ids_all - cohort_ids)} outside the cohort), v7.1 collapsed flag {len(collapsed_v7)}, "
      f"v7.1 column set t90_final {len(collapse_cols['t90_final'])} / master {len(collapse_cols['master'])} / base4 {len(collapse_cols['t90_base4'])}")

    # ------------------------------------------------------------------ measure columns (137 + architecture + indices + oxygen)
    fam = pd.read_csv(CS.MEASURE_FAMILIES_CSV)
    feats = fam.feature.tolist()
    assert not any(f.endswith(CS.RANKING_DROP_SUFFIX) for f in feats), "measure_families.csv still lists _siteZ twins"
    assert len(feats) == CS.N_MEASURES, f"measure_families.csv has {len(feats)} rows, cohort_spec.N_MEASURES = {CS.N_MEASURES}"
    meas137 = [f for f in feats if f not in NEW4]
    assert len(meas137) == CS.N_MEASURES - len(NEW4)
    repl = list(dict.fromkeys(meas137 + ARCH + IDX + list(OXY_MAP)))
    srcmap = {c: OXY_MAP.get(c, c) for c in repl}
    oxy_cols = list(OXY_MAP)
    sleep_gated = [c for c in repl if c not in oxy_cols and c != "recording_dur_min"]

    # ------------------------------------------------------------------ the two passes
    split = load_pass(a.split, "split", collapsed_v7, S.MANIFEST_SPLIT)
    light = load_pass(a.light, "light", collapsed_v7, S.MANIFEST_FULL)
    if a.ids:
        ids = pd.read_csv(a.ids); col = "BDSPPatientID" if "BDSPPatientID" in ids.columns else ids.columns[0]
        keep = set(ids[col]); split = split[split.BDSPPatientID.isin(keep)].copy(); light = light[light.BDSPPatientID.isin(keep)].copy()
        P(f"--ids {a.ids}: {len(keep):,} ids; split {len(split):,} rows, light {len(light):,} rows kept")
    split = split[split.BDSPPatientID.isin(split_ids_all)].copy()       # the split pass replaces only manifest split nights
    split = detect_collapse(split); light = detect_collapse(light)
    missing_src = [c for c in repl if srcmap[c] not in split.columns]
    if missing_src:
        sys.exit(f"the split pass lacks {len(missing_src)} measure columns: {missing_src}")
    split_ok = split[split._ok].copy(); split_fail = split[~split._ok]
    light_ok = light[light._ok].copy(); light_fail = light[~light._ok]
    if not a.smoke:
        miss_split = split_ids_all - set(split.BDSPPatientID); miss_light = cohort_ids - set(light.BDSPPatientID)
        if miss_split or miss_light:
            sys.exit(f"REAL build: split pass lacks {len(miss_split)} manifest nights, light pass lacks {len(miss_light)} cohort nights")
    P(f"split pass: {len(split):,} manifest nights present, ok {len(split_ok):,}, failed {len(split_fail):,} "
      f"{split_fail.status.astype(str).str[:30].value_counts().to_dict() if len(split_fail) else ''}; "
      f"light pass: {len(light):,} present, ok {len(light_ok):,}, failed {len(light_fail):,}")
    # windowed rows of the light pass that are not manifest split nights would be a manifest error
    stray = light_ok[light_ok._v8_prepap_window & ~light_ok.BDSPPatientID.isin(split_ids_all)]
    assert len(stray) == 0, f"{len(stray)} light rows windowed but not in the split manifest"

    # D1 counts on the split pass (X1 applied the rules; the builder counts and carries them)
    tst_w = pd.to_numeric(split_ok.TST_min, errors="coerce")
    n_under60 = int((tst_w < S.MIN_TST_FOR_INDEX).sum()); n_under30 = int((tst_w < 30).sum())
    n_noscored = int(((tst_w <= 0) | tst_w.isna() | (pd.to_numeric(split_ok.get("n_scored_epochs"), errors="coerce") == 0)).sum())
    n_oxyshort = int(split_ok._v8_oxygen_window_short.sum()); n_idx_undef = int(split_ok._v8_index_undefined.sum())
    n_coll_split = int(split_ok._v8_stage_collapsed.sum()); n_coll_light = int(light_ok._v8_stage_collapsed.sum())
    n_short_single = int(split_ok._v8_single_stage_short_window.sum())
    P(f"D1 on the pre-PAP window (ok split nights): true sleep under 60 min {n_under60:,} (under 30 min {n_under30:,}); indices undefined {n_idx_undef:,}; "
      f"no scored sleep epoch {n_noscored:,}; oxygen window under {S.MIN_OXYGEN_VALID_MIN:g} valid min {n_oxyshort:,} (policy {a.short_oxygen_window})")
    P(f"collapsed staging in the pass (one sleep code, TST > 0): split {n_coll_split:,} of {len(split_ok):,} ok; light {n_coll_light:,} of {len(light_ok):,} ok "
      f"(v7.1 flag marks {int(light_ok._v7_1_stage_collapsed.sum())} of the light rows); windowed single non-REM code under 60 min, kept: {n_short_single:,}")

    # ------------------------------------------------------------------ the replacement source frame (split nights)
    src = split_ok.set_index("BDSPPatientID").copy()
    for c in repl:
        src[srcmap[c]] = pd.to_numeric(src[srcmap[c]], errors="coerce")
    fzi = prev_fz.set_index("BDSPPatientID")
    n_floor = 0
    if a.sleep_floor_all:
        m = pd.to_numeric(src.TST_min, errors="coerce") < S.MIN_TST_FOR_INDEX
        for c in sleep_gated:
            src.loc[m, srcmap[c]] = np.nan
        n_floor = int(m.sum())
        P(f"D1 ALTERNATIVE applied (--sleep-floor-all): every sleep-gated measure missing on {n_floor:,} split nights under 60 min")
    # collapsed staging in the window: v7.1 rule (composition + stage-gated set), applied per table below through the src frame
    src_collapsed = src[src._v8_stage_collapsed].index
    # short oxygen window policy
    keep_oxy_ids = src[src._v8_oxygen_window_short].index if a.short_oxygen_window == "keep_v7" else pd.Index([])
    if len(keep_oxy_ids):
        for c in oxy_cols:
            src.loc[keep_oxy_ids, srcmap[c]] = fzi.loc[keep_oxy_ids, c].values     # a no-op replacement, recorded as kept
    src["_v8_measure_source"] = "pre-PAP window (split pass)"
    src.loc[keep_oxy_ids, "_v8_measure_source"] = f"pre-PAP window (split pass); oxygen family v7 whole-night (pre-PAP oxygen window under {S.MIN_OXYGEN_VALID_MIN:g} valid min)"

    def replace_split(t, name):
        idx = t.BDSPPatientID.isin(src.index).to_numpy()
        ids = t.loc[idx, "BDSPPatientID"].values; s = src.reindex(ids)
        ncell = 0; ncol = 0
        ccols = set(collapse_cols[name])
        for c in repl:
            if c not in t.columns:
                continue
            ensure_float(t, c)
            new = s[srcmap[c]].to_numpy(dtype=float)
            if c in ccols:
                new = np.where(s._v8_stage_collapsed.to_numpy(dtype=bool), np.nan, new)
            n = CH.add(name, ids, c, t.loc[idx, c].values, new, "pre-PAP window (decision 1)")
            t.loc[idx, c] = new; ncell += n; ncol += 1
        if "spo2_native_fs" in t.columns and "spo2_native_fs" in s.columns:
            ensure_float(t, "spo2_native_fs"); t.loc[idx, "spo2_native_fs"] = pd.to_numeric(s.spo2_native_fs, errors="coerce").values
        # collapsed windows: the v7.1 column set that is not among the replaced columns (siteZ twins etc.) -> missing too
        cidx = t.BDSPPatientID.isin(src_collapsed).to_numpy()
        for c in ccols - set(repl):
            if c in t.columns:
                ensure_float(t, c)
                CH.add(name, t.loc[cidx, "BDSPPatientID"].values, c, t.loc[cidx, c].values, np.full(int(cidx.sum()), np.nan), "stage collapsed in the pre-PAP window (v7.1 rule)")
                t.loc[cidx, c] = np.nan
        P(f"{name}: split replacement on {int(idx.sum()):,} nights, {ncol} columns, {ncell:,} cells changed; collapsed windows {int(cidx.sum()):,}")
        return t

    fz = replace_split(fz, "t90_final"); b4 = replace_split(b4, "t90_base4"); ma = replace_split(ma, "master")

    # ------------------------------------------------------------------ new measures (all cohort nights, from the light pass)
    lsrc = light_ok.set_index("BDSPPatientID")
    # a light record windowed on a split night is used only when that night's measures were replaced by the split pass (the same
    # window); otherwise the whole-night table and a windowed new measure would disagree. In the real build every windowed night
    # is replaced, so nothing is skipped; the smoke (subset passes) reports the skipped count.
    skip_w = lsrc._v8_prepap_window & ~lsrc.index.isin(src.index)
    n_skip_w = int(skip_w.sum()); lsrc = lsrc[~skip_w].copy()
    lsrc["_v8_light_tst_min"] = pd.to_numeric(lsrc.TST_min, errors="coerce")
    for c in NEW5 + QC_NEW:
        if c in lsrc.columns and c != "hb_window_source":
            lsrc[c] = pd.to_numeric(lsrc[c], errors="coerce")
    # D11 on the histogram-detected collapse: only the composition-gated measure (REM latency) goes missing; hypoxic burden, the PLM
    # index, sleep latency and sleep-period T90 are sleep-versus-wake quantities, which the 2026-09-08 addendum keeps on collapsed
    # nights ("sleep-versus-wake quantities may stay"). Coordinator default 2026-09-12; alternative (all five missing) in the report.
    D11_COLS = ["rem_latency_min"]
    for c in D11_COLS:
        lsrc.loc[lsrc._v8_stage_collapsed, c] = np.nan
    coll_new_ids = set(lsrc[lsrc._v8_stage_collapsed].index)
    # cross-check: the split pass carries the same new measures on the same window
    both = src.index.intersection(lsrc.index)
    xchk = {}
    for c in NEW5:
        if c in src.columns and len(both):
            x = pd.to_numeric(src.loc[both, c], errors="coerce"); y = lsrc.loc[both, c]
            xchk[c] = int((~((x == y) | (x.isna() & y.isna()))).sum())
    P(f"new measures: light ok nights in the cohort {int(lsrc.index.isin(cohort_ids).sum()):,} (windowed light records skipped because the split pass "
      f"does not cover the night: {n_skip_w:,}{', a smoke artefact' if a.smoke else ', MUST be 0 in a real build'}); D11 (REM latency only) missing on {len(coll_new_ids):,} collapsed nights; "
      f"split-vs-light disagreement on the {len(both):,} shared nights: {xchk}")
    if n_skip_w and not a.smoke:
        sys.exit("REAL build: windowed light records on nights the split pass did not replace")

    def add_new(t, name, cols):
        idx = t.BDSPPatientID.isin(lsrc.index).to_numpy() & t.BDSPPatientID.isin(cohort_ids).to_numpy()
        ids = t.loc[idx, "BDSPPatientID"].values; s = lsrc.reindex(ids)
        for c in cols:
            if c not in s.columns:
                continue
            if c not in t.columns:
                t[c] = np.nan if c != "hb_window_source" else None
            if c != "hb_window_source":
                ensure_float(t, c)
                new = s[c].to_numpy(dtype=float)
            else:
                new = s[c].astype(object).where(s[c].notna(), None).to_numpy()
            CH.add(name, ids, c, t.loc[idx, c].values, new, "new measure (decision 3 / 2)" if c in NEW5 else "new measure QC column")
            t[c] = t[c].astype(float) if c != "hb_window_source" else t[c]
            t.loc[idx, c] = new
        return t

    fz = add_new(fz, "t90_final", NEW5 + QC_NEW + ["_v8_light_tst_min"]); b4 = add_new(b4, "t90_base4", NEW5); ma = add_new(ma, "master", NEW5)
    P(f"new columns added: t90_final {NEW5 + [c for c in QC_NEW if c in lsrc.columns]}, t90_base4 and master {NEW5}")

    # ------------------------------------------------------------------ plausibility bounds on every replaced or added value
    if bounds_src == "csv":
        bounds = pd.read_csv(bounds_path)[["column", "lo", "hi"]]
    else:
        r = v7ch[v7ch.reason.str.startswith("impossible value outside")].drop_duplicates(["column", "reason"])
        bounds = pd.DataFrame({"column": r.column.values,
                               "lo": [float(s.split("[")[1].split(",")[0]) for s in r.reason],
                               "hi": [float(s.split(",")[1].split("]")[0]) for s in r.reason]}).drop_duplicates("column")
    unchecked_bounded = sorted({c[len("_impossible_"):] for c in prev_ma.columns if c.startswith("_impossible_")} - set(bounds.column))
    touched_ids = set(src.index) | set(lsrc.index)
    nimp = 0
    for _, r in bounds.iterrows():
        col, lo, hi = r["column"], float(r["lo"]), float(r["hi"])
        for t, name in ((fz, "t90_final"), (b4, "t90_base4"), (ma, "master")):
            if col not in t.columns:
                continue
            m = t.BDSPPatientID.isin(touched_ids) & t.BDSPPatientID.isin(cohort_ids) & t[col].notna() & ((t[col] < lo) | (t[col] > hi))
            flag = f"_impossible_{col}"
            if flag not in t.columns:
                t[flag] = False
            t[flag] = _b(t[flag]) | m.to_numpy()
            if m.any():
                nimp += CH.add(name, t.loc[m, "BDSPPatientID"].values, col, t.loc[m, col].values, np.full(int(m.sum()), np.nan), f"impossible value outside [{lo:g}, {hi:g}]")
                t.loc[m, col] = np.nan
    P(f"plausibility bounds ({bounds_src}, {len(bounds)} columns): {nimp:,} replaced/added cells outside the bounds set missing + flagged; "
      f"bounded columns without a bound on file: {unchecked_bounded or 'none'}")

    # ------------------------------------------------------------------ reconciliation of the light pass against v7 (non-split nights)
    ns = lsrc[~lsrc._v8_prepap_window & lsrc.index.isin(cohort_ids)]
    v7n = fzi.reindex(ns.index)
    # The v7 oxygen columns came from the v7 native oximeter pass (oxygen_full_FINAL); the light pass is the same definition in
    # another run (it equals reextract_FINAL exactly). Differences are summation-order noise (T90 1e-15, mean 1e-5, 1st percentile
    # 1e-7) and are counted; a definitional shift would exceed the tolerances below and stops the build.
    recon = {}; n_recon_bad = 0
    for c, sc in OXY_MAP.items():
        flag = f"_impossible_{c}"
        ok_rows = ~v7n[flag].astype(bool) if flag in v7n.columns else pd.Series(True, index=ns.index)
        x = pd.to_numeric(ns[sc], errors="coerce"); y = pd.to_numeric(v7n[c], errors="coerce")
        d = (x - y).abs().where(~(x.isna() & y.isna()), 0.0).where(~(x.isna() ^ y.isna()), np.inf)[ok_rows]
        tol = RECON_TOL[c]
        recon[c] = {"n": int(len(d)), "differ_1e-9": int((d > 1e-9).sum()), "max_abs": float(d[np.isfinite(d)].max()) if np.isfinite(d).any() else 0.0, "beyond_tol": int((d > tol).sum()), "tol": tol}
        n_recon_bad += recon[c]["beyond_tol"]
    tst_ok = ~v7n["_v7_stagefix_source"].astype(str).str.startswith("April") & ~v7n["_v7_1_stage_collapsed"].astype(bool) & ~ns._v8_stage_collapsed
    x = pd.to_numeric(ns.TST_min, errors="coerce"); y = v7n.TST_min
    d = (x - y).abs()[tst_ok]
    recon["TST_min (current-file staging vs v7, excl. April-record and collapsed nights; informational)"] = {"n": int(tst_ok.sum()), "differ_1e-9": int((d > 1e-9).sum()), "max_abs_min": float(d.max()) if len(d) else 0.0}
    P(f"reconciliation, light pass vs v7 t90_final on {len(ns):,} non-split cohort nights: " + "; ".join(f"{k}: {v}" for k, v in recon.items()))
    if n_recon_bad and not a.allow_recon_mismatch:
        sys.exit(f"RECONCILIATION MISMATCH: {n_recon_bad} oxygen cells of the light pass differ from v7 beyond tolerance on non-split nights (gate 1); pass --allow-recon-mismatch to continue")
    # the pre-PAP T90 against the cpap table's own pre-PAP T90 (the v7 pilot value carried in the split manifest)
    if "pre_any_t90_v7" in man_split.columns:
        cmp_ids = src.index.difference(keep_oxy_ids)
        mm = man_split.set_index("BDSPPatientID").reindex(cmp_ids)
        d = (pd.to_numeric(src.loc[cmp_ids, OXY_MAP["spo2_pct_below_90"]], errors="coerce") - pd.to_numeric(mm.pre_any_t90_v7, errors="coerce")).abs()
        P(f"pre-PAP T90 vs the cpap table's own pre_any_t90 on {int(d.notna().sum()):,} split nights (the {len(keep_oxy_ids)} keep_v7 nights excluded): "
          f"identical (1e-6) on {int((d < 1e-6).sum()):,}, max |diff| {float(d.max()) if d.notna().any() else float('nan'):.4f}")

    # ------------------------------------------------------------------ outcomes (t90_final, t90_base4)
    oc = pd.read_parquet(a.outcomes)
    assert len(oc) == n_rows_v7 and set(oc.BDSPPatientID) == set(prev_fz.BDSPPatientID), "outcomes_v8 rows differ from the v7 table"
    ocx = oc.set_index("BDSPPatientID")
    counts = pd.read_csv(S.EVENT_COUNTS_V8)
    copied_keys = counts.loc[counts.status.str.startswith("unchanged"), "key"].tolist()
    changed_keys = counts.loc[~counts.status.str.startswith("unchanged") & (counts.key != "death"), "key"].tolist()
    ohs_note = "not rebuilt (--no-ohs-rebuild)"
    ohs_flag_diff = None
    if not a.no_ohs_rebuild:
        import duckdb, importlib.util
        spec = importlib.util.spec_from_file_location("disease_definitions_v8", S.DEFS_V8); D8 = importlib.util.module_from_spec(spec); spec.loader.exec_module(D8)
        lab, neg, i9, i10 = D8.DISEASES["obesity_hypovent"]; excl = getattr(D8, "EXCLUDE", {}).get("obesity_hypovent", ())
        pref = sorted(set(i9) | set(i10)); where = " or ".join(f"code like '{p}%'" for p in pref)
        if excl:
            where = f"({where}) and not (" + " or ".join(f"code like '{e}%'" for e in excl) + ")"
        con = duckdb.connect()
        fd = con.execute(f"""select person_id as BDSPPatientID, min(condition_start_date) as first_date from
                             (select person_id, condition_start_date, upper(replace(condition_source_value, '.', '')) as code
                              from '{S.OMOP_CACHE}' where condition_source_value is not null) where {where} group by person_id""").df()
        con.close()
        fd["first_date"] = pd.to_datetime(fd.first_date, errors="coerce"); fd = fd.set_index("BDSPPatientID")
        new_fd = fd.first_date.reindex(ocx.index).astype("datetime64[us]")
        old_fd = ocx["obesity_hypovent_first_date"].astype("datetime64[us]")
        psg = ocx.psg_date
        prev_flag = (new_fd.notna() & psg.notna() & (new_fd <= psg)).astype(int); inc_flag = (new_fd.notna() & psg.notna() & (new_fd > psg)).astype(int)
        cohm = ocx.index.isin(cohort_ids)
        ohs_flag_diff = {"prevalent_cohort_rows_differ": int(((prev_flag != ocx.obesity_hypovent_prevalent) & cohm).sum()),
                         "incident_cohort_rows_differ": int(((inc_flag != ocx.obesity_hypovent_incident) & cohm).sum()),
                         "first_date_rows_changed_all": int((~((new_fd == old_fd) | (new_fd.isna() & old_fd.isna()))).sum()),
                         "first_date_rows_changed_cohort": int((~((new_fd == old_fd) | (new_fd.isna() & old_fd.isna())) & cohm).sum()),
                         "prefixes": pref, "exclude": list(excl)}
        ocx["obesity_hypovent_first_date"] = new_fd.values
        ohs_note = (f"obesity_hypovent_first_date rebuilt from the hospital-record cache with the v8 list {pref} (verbatim v7 rule: strip the point, "
                    f"upper-case, prefix match, earliest date per patient); flags kept from X2/v7; flags recomputed from the rebuilt dates differ from the "
                    f"kept flags on {ohs_flag_diff['prevalent_cohort_rows_differ']} (prevalent) and {ohs_flag_diff['incident_cohort_rows_differ']} (incident) cohort rows; "
                    f"first dates changed on {ohs_flag_diff['first_date_rows_changed_cohort']:,} cohort rows ({ohs_flag_diff['first_date_rows_changed_all']:,} of all rows)")
        P(ohs_note)
    out_cols = [c for c in ocx.columns if c.endswith(S.OUTCOME_SUFFIXES) or c.startswith("death_") or c.startswith("_v8_outcome")]

    def merge_outcomes(t, name):
        ids = t.BDSPPatientID.values; s = ocx.reindex(ids)
        n_rep = n_add = ncell = 0
        for c in out_cols:
            new = s[c].values
            if c in t.columns:
                if pd.api.types.is_datetime64_any_dtype(t[c]) or pd.api.types.is_datetime64_any_dtype(s[c]):
                    new = pd.to_datetime(s[c]).astype("datetime64[us]").values
                    old = pd.to_datetime(t[c]).astype("datetime64[us]").values
                    ncell += CH.add(name, ids, c, old, new, "outcome v8 (decision 5 / 11; OHS first date rebuilt)")
                    t[c] = new
                else:
                    if pd.api.types.is_float_dtype(s[c]) and not pd.api.types.is_float_dtype(t[c]):
                        ensure_float(t, c)
                    ncell += CH.add(name, ids, c, t[c].values, new, "outcome v8 (decision 5 / 11)")
                    t[c] = new
                n_rep += 1
            else:
                t[c] = new; n_add += 1
                if c.endswith(S.OUTCOME_SUFFIXES):
                    ncell += CH.add(name, ids, c, np.full(len(ids), np.nan, dtype=object), new, "outcome v8, new condition (decision 11)")
        P(f"{name}: outcomes merged, {n_rep} columns replaced, {n_add} added, {ncell:,} cells changed")
        return t

    fz = merge_outcomes(fz, "t90_final"); b4 = merge_outcomes(b4, "t90_base4")
    if a.outcomes_into_master:
        ma = merge_outcomes(ma, "master")
    # positive control (X2's gate re-asserted on the merged table): untouched conditions byte-identical to v7
    bad = []
    rebuilt = set() if a.no_ohs_rebuild else {"obesity_hypovent_first_date"}     # rebuilt on purpose; its flags stay byte-identical
    for k in copied_keys:
        for suf in S.OUTCOME_SUFFIXES:
            c = k + suf
            if c in prev_fz.columns and c not in rebuilt and not nan_equal(fz[c].values, prev_fz[c].values).all():
                bad.append(c)
    assert not bad, f"untouched outcome columns changed: {bad}"
    fz = fz.copy(); b4 = b4.copy(); ma = ma.copy()      # de-fragment after the column additions
    coh8 = CS.apply_cohort(fz)
    inc_ref = counts.set_index("key")["incident"]
    inc_bad = {k: (int(coh8[f"{k}_incident"].sum()), int(inc_ref[k])) for k in changed_keys if int(coh8[f"{k}_incident"].sum()) != int(inc_ref[k])}
    assert not inc_bad, f"incident counts in the merged table differ from X2's event_counts_v8.csv: {inc_bad}"
    P(f"outcome positive control: {len(copied_keys)} untouched conditions byte-identical to v7 in t90_final; {len(changed_keys)} changed/new conditions' incident counts equal event_counts_v8.csv")

    # ------------------------------------------------------------------ cpap_t90_by_stage: whole-night covariates from the v8 t90_final
    v8i = fz.set_index("BDSPPatientID"); idx = cp.BDSPPatientID.isin(v8i.index).to_numpy(); ids = cp.loc[idx, "BDSPPatientID"].values
    for c in ("AHI", "TST_min", "spo2_pct_below_90"):
        if c in cp.columns:
            ensure_float(cp, c)
            CH.add("cpap_t90_by_stage", ids, c, cp.loc[idx, c].values, v8i.loc[ids, c].values, "whole-night covariate from v8 t90_final")
            cp.loc[idx, c] = v8i.loc[ids, c].values
    cp["_v8_covariates_refreshed"] = idx
    P(f"cpap_t90_by_stage: whole-night AHI, TST, T90 refreshed on {int(idx.sum()):,} of {len(cp):,} rows (its per-stage columns unchanged)")

    # ------------------------------------------------------------------ provenance columns
    split_src_tag = split._src_tag.iloc[0] if len(split) else ""; light_src_tag = light._src_tag.iloc[0] if len(light) else ""
    hb_sha = light_ok["_v8_hb_definition_sha"].dropna().iloc[0] if "_v8_hb_definition_sha" in light_ok.columns and light_ok["_v8_hb_definition_sha"].notna().any() else \
        (S.sha256(S.HB_DEFINITION)[:16] if os.path.isfile(S.HB_DEFINITION) and not S.evicted(S.HB_DEFINITION) else None)
    defs_sha = S.sha256(S.DEFS_V8)[:12]
    spl_all = split.set_index("BDSPPatientID")

    def stamp(t, name, outcomes_here):
        ids = t.BDSPPatientID
        in_split_ok = ids.isin(src.index).to_numpy(); in_split_any = ids.isin(spl_all.index).to_numpy()
        s = src.reindex(ids.values); sa = spl_all.reindex(ids.values); ls = lsrc.reindex(ids.values)
        t["_v8_prepap_window"] = in_split_ok
        t["_v8_window_end_sec"] = np.where(in_split_any, pd.to_numeric(sa.get("t_end_sec"), errors="coerce"), np.nan)
        t["_v8_window_end_min"] = np.where(in_split_ok, pd.to_numeric(s.get("_v8_window_end_min"), errors="coerce"), np.nan)
        t["_v8_untreated_sleep_min"] = np.where(in_split_ok, pd.to_numeric(s.get("TST_min"), errors="coerce"), np.nan)
        t["_v8_session_differs_from_v7"] = in_split_any & _b(sa["_v8_session_differs_from_v7"])
        t["_v8_index_undefined"] = in_split_ok & _b(s["_v8_index_undefined"])
        t["_v8_oxygen_window_short"] = in_split_ok & _b(s["_v8_oxygen_window_short"])
        t["_v8_sleep_floor_all"] = bool(a.sleep_floor_all) & in_split_ok & _b(pd.to_numeric(s.get("TST_min"), errors="coerce") < S.MIN_TST_FOR_INDEX)
        t["_v8_stage_collapsed"] = (in_split_ok & _b(s["_v8_stage_collapsed"])) | ids.isin(coll_new_ids).to_numpy()
        t["_v8_single_stage_short_window"] = in_split_ok & _b(s["_v8_single_stage_short_window"])
        msrc = np.where(in_split_ok, s["_v8_measure_source"].astype(object),
                        np.where(in_split_any, "v7 whole-night (v8 split extraction failed: " + sa.status.astype(str).str[:40] + ")",
                                 np.where(ids.isin(cohort_ids), "v7 whole-night (not a split night)", "v7 (outside the analysis cohort)")))
        t["_v8_measure_source"] = msrc
        prev_ahi = t["_v7_ahi_source"].astype(str) if "_v7_ahi_source" in t.columns else pd.Series("v7", index=t.index)
        t["_v8_ahi_source"] = np.where(in_split_ok, np.where(t["_v8_index_undefined"], "undefined, pre-PAP true sleep under 60 min", "pre-PAP window count over pre-PAP true sleep (split pass)"), prev_ahi)
        t["_v8_prepap_source"] = np.where(in_split_ok, split_src_tag, "")
        t["_v8_new_measure_source"] = np.where(ids.isin(lsrc.index) & ids.isin(cohort_ids), light_src_tag, "")
        t["_v8_light_status"] = np.where(ids.isin(light.BDSPPatientID), light.set_index("BDSPPatientID").status.reindex(ids.values).astype(str), "not run")
        if name == "t90_final":
            t["_v8_hb_definition_sha"] = np.where(ids.isin(lsrc.index) & ids.isin(cohort_ids), hb_sha, None)
        if outcomes_here:
            t["_v8_outcome_list_version"] = defs_sha
        return t

    fz = stamp(fz, "t90_final", True); b4 = stamp(b4, "t90_base4", True); ma = stamp(ma, "master", a.outcomes_into_master)

    # ------------------------------------------------------------------ drop the _siteZ twins (master, base4)
    for t, name in ((ma, "master"), (b4, "t90_base4")):
        z = [c for c in t.columns if c.endswith(CS.RANKING_DROP_SUFFIX)]
        CH.dropped[name] = z
        t.drop(columns=z, inplace=True)
    P(f"_siteZ columns dropped: master {len(CH.dropped['master'])}, t90_base4 {len(CH.dropped['t90_base4'])}")

    # ------------------------------------------------------------------ v8.1 (Alen, 14 Sept): split nights contribute oxygen and event-rate measures only
    import importlib.util as _ilu
    masked_ids = set()

    def mask_split_sleep(t, name):
        m = t["_v8_prepap_window"].to_numpy(dtype=bool)
        ids = t.loc[m, "BDSPPatientID"].values; ncell = 0; done = []
        for c in S.SLEEP_MASK_COLS:
            if c not in t.columns:
                continue
            ensure_float(t, c)
            old = t.loc[m, c].values
            ncell += CH.add(name, ids, c, old, np.full(len(old), np.nan), S.SLEEP_MASK_REASON)
            t.loc[m, c] = np.nan; done.append(c)
        t["_v8_sleep_masked"] = m
        masked_ids.update(ids.tolist())
        P(f"{name}: v8.1 split-night sleep mask: {len(done)} columns, {ncell:,} cells set missing on {int(m.sum()):,} nights ({', '.join(done)})")
        return ncell

    n_mask = mask_split_sleep(fz, "t90_final") + mask_split_sleep(b4, "t90_base4") + mask_split_sleep(ma, "master")
    if "TST_min" in cp.columns:           # the whole-night covariate copy in cpap_t90_by_stage follows t90_final
        mcp = cp.BDSPPatientID.isin(masked_ids).to_numpy(); ensure_float(cp, "TST_min")
        n_mask += CH.add("cpap_t90_by_stage", cp.loc[mcp, "BDSPPatientID"].values, "TST_min", cp.loc[mcp, "TST_min"].values,
                         np.full(int(mcp.sum()), np.nan), S.SLEEP_MASK_REASON)
        cp.loc[mcp, "TST_min"] = np.nan
    # v8.1: the two swapped-out negative controls leave the outcome table (their v7 columns are not carried forward)
    _sp = _ilu.spec_from_file_location("disease_definitions_v8_swap", S.DEFS_V8); _dm = _ilu.module_from_spec(_sp); _sp.loader.exec_module(_dm)
    gone_cols = [k + suf for k in _dm.SWAPPED_CONTROLS_2026_09_14 for suf in S.OUTCOME_SUFFIXES]
    for t, name in ((fz, "t90_final"), (b4, "t90_base4"), (ma, "master")):
        g = [c for c in gone_cols if c in t.columns]
        CH.dropped[name] = CH.dropped.get(name, []) + g
        t.drop(columns=g, inplace=True)
    P(f"v8.1 swapped-out control columns dropped where present: {gone_cols}")

    # ------------------------------------------------------------------ untouched cells byte-identical to v7 (in memory)
    ch = CH.frame()
    for t, prev, name in ((fz, prev_fz, "t90_final"), (b4, prev_b4, "t90_base4"), (ma, prev_ma, "master"), (cp, prev_cp, "cpap_t90_by_stage")):
        assert len(t) == len(prev) and (t.BDSPPatientID.values == prev.BDSPPatientID.values).all(), f"{name}: row order or count changed"
        chn = ch[ch.table == name]
        bad = []
        for c in prev.columns:
            if c not in t.columns:
                if c not in CH.dropped.get(name, []):
                    bad.append((c, "column missing"))
                continue
            if c.startswith("_impossible_"):
                continue
            a_, b_ = t[c].values, prev[c].values
            if pd.api.types.is_datetime64_any_dtype(prev[c]) or pd.api.types.is_datetime64_any_dtype(t[c]):
                a_ = pd.to_datetime(t[c]).astype("datetime64[us]").values; b_ = pd.to_datetime(prev[c]).astype("datetime64[us]").values
            diff_rows = set(prev.BDSPPatientID.values[~nan_equal(a_, b_)])
            lc = chn.loc[chn.column == c, "BDSPPatientID"].value_counts()
            unlogged = diff_rows - set(lc.index)
            # a cell logged twice may net to its v7 value (pre-PAP value, then the bounds rule back to missing): a chain, accepted
            single_equal = {i for i in set(lc.index) - diff_rows if lc[i] < 2}
            if unlogged or single_equal:
                bad.append((c, f"{len(unlogged)} unlogged / {len(single_equal)} logged once but equal"))
        assert not bad, f"{name}: untouched cells changed or the CHANGELOG is incomplete: {bad[:10]}"
        n_same = sum(1 for c in prev.columns if c in t.columns and c not in set(chn.column))
        P(f"{name}: {n_same} of {prev.shape[1]} v7 columns byte-identical; {chn.column.nunique()} columns changed on {chn.BDSPPatientID.nunique():,} nights ({len(chn):,} cells); {len(t.columns) - prev.shape[1] + len(CH.dropped.get(name, []))} columns added")
    # the non-split cohort nights: every measure column byte-identical
    nonsplit = cohort_ids - split_ids_all
    for t, prev, name in ((fz, prev_fz, "t90_final"), (b4, prev_b4, "t90_base4"), (ma, prev_ma, "master")):
        m = t.BDSPPatientID.isin(nonsplit).to_numpy()
        badc = [c for c in repl if c in prev.columns and not nan_equal(t.loc[m, c].values, prev.loc[m, c].values).all()]
        assert not badc, f"{name}: measure columns changed on non-split nights: {badc}"
    P(f"non-split cohort nights ({len(nonsplit):,}): all {len([c for c in repl if c in prev_ma.columns])} measure columns byte-identical to v7 in every table")
    assert len(CS.apply_cohort(fz)) == len(cohort_ids) and len(CS.apply_cohort(b4)) == len(cohort_ids), "the analysis cohort changed"

    # ------------------------------------------------------------------ write
    files = {"t90_final.parquet": fz, "t90_base4.parquet": b4, "master_cohort.csv": ma, "cpap_t90_by_stage.parquet": cp}
    for f, t in files.items():
        p = f"{out}/{f}"
        (t.to_csv(p, index=False) if f.endswith(".csv") else t.to_parquet(p, index=False))
    ch.to_parquet(f"{out}/CHANGELOG.parquet", index=False)
    outs = {f: S.sha256(f"{out}/{f}") for f in files}
    with open(f"{out}/sha256.txt", "w") as f:
        for k, v in outs.items():
            f.write(f"{v}  {k}\n")
    summ = ch.groupby(["table", "column", "reason"]).size().reset_index(name="cells")
    by_reason = ch.groupby(["table", "reason"]).agg(cells=("column", "size"), columns=("column", "nunique"), nights=("BDSPPatientID", "nunique")).reset_index()
    in_sha = {p: S.sha256(p) for p in inputs if os.path.isfile(p)}
    decisions = [
        f"D1 (default, v7 rules): AHI and arousal index undefined on {n_idx_undef:,} split nights with under 60 min of pre-PAP sleep "
        f"(under 30 min: {n_under30:,}); no scored sleep epoch on {n_noscored:,}. ALTERNATIVE (--sleep-floor-all, {'APPLIED' if a.sleep_floor_all else 'not applied'}): "
        f"every sleep-gated measure missing on those {n_under60:,} nights.",
        f"D1 oxygen: {n_oxyshort:,} split nights have under {S.MIN_OXYGEN_VALID_MIN:g} valid pre-PAP oximeter minutes. Policy --short-oxygen-window {a.short_oxygen_window}: "
        + ("the six oxygen measures keep their v7 whole-night values on those nights (flag _v8_oxygen_window_short; the cohort stays intact). ALTERNATIVE nan: "
           "the six become missing and the analysis cohort loses those nights (cohort_spec.COHORT_N would change)." if a.short_oxygen_window == "keep_v7" else
           "the six oxygen measures are missing on those nights; the analysis cohort shrinks by that count."),
        f"Collapsed staging (addendum 2026-09-08 gate 2, D11): detected from the pass histogram on {n_coll_split:,} split and {n_coll_light:,} light nights "
        f"(the v7.1 flag alone marks {len(collapsed_v7)} rows); composition and the v7.1 stage-gated column set set missing on the collapsed split windows, "
        f"the five new measures missing on every collapsed night. Windowed single non-REM-code records under 60 min ({n_short_single:,}) are kept as short windows.",
        f"Outcomes: {ohs_note}.",
        "Outcome clock: v7 psg_date and censor_date (plan D2 default), carried by outcomes_v8.parquet.",
        f"master_cohort.csv carries the four new measures and {SLEEP_T90} only (the ranking scans its header for candidates: N_MEASURES = {CS.N_MEASURES}); "
        f"outcome columns {'were merged into master too (--outcomes-into-master)' if a.outcomes_into_master else 'are NOT in master (no live script reads them there; flag --outcomes-into-master)'}.",
        f"_siteZ twins dropped from master ({len(CH.dropped['master'])}) and t90_base4 ({len(CH.dropped['t90_base4'])}); t90_final never carried them.",
        f"Plausibility bounds from {bounds_src}; bounded columns without a bound on file: {unchecked_bounded or 'none'}.",
        f"Session flag: {int(fz['_v8_session_differs_from_v7'].sum())} split nights whose PAP switch-on session differs from the v7 feature session (X1 flag 1).",
    ]
    prov = [f"# PROVENANCE, T90 tables v8 (2026-09){' SMOKE BUILD, not the frozen tables' if a.smoke else ''}", "",
            f"Built {time.strftime('%Y-%m-%d %H:%M %Z')} by {os.path.abspath(__file__)} (sha256 {S.sha256(os.path.abspath(__file__))[:16]})",
            f"Mode: {'smoke' if a.smoke else 'real'}; --out {out}; --split {a.split}; --light {a.light}; --outcomes {a.outcomes}; --ids {a.ids or 'none'}",
            f"Spec: {os.path.abspath(CS.__file__)} (PREV_DIR {S.PREV_DIR}; COHORT_N {CS.COHORT_N}; N_MEASURES {CS.N_MEASURES}; N_RANKED_OUTCOMES {CS.N_RANKED_OUTCOMES})",
            "", "## Inputs (sha256)"] + [f"- {os.path.basename(p)}: {h}  ({p})" for p, h in in_sha.items()] + \
           ["", "## Outputs (sha256)"] + [f"- {k}: {v}" for k, v in outs.items()] + \
           ["", "## Decisions applied as defaults (gate 12), with the alternative"] + [f"- {d}" for d in decisions] + \
           ["", "## Log"] + [f"- {s}" for s in LOG] + \
           ["", "## Changed cells by table and reason", "", by_reason.to_string(index=False),
            "", "## Changed cells by table, column and reason", "", summ.to_string(index=False),
            "", "## v8.1 split-night sleep mask (Alen, 14 Sept 2026)", "",
            f"Split nights contribute oxygen and event-rate measures only: sleep amounts, stage composition, REM latency and event counts "
            f"({len(S.SLEEP_MASK_COLS)} columns: {', '.join(S.SLEEP_MASK_COLS)}) are set missing on the {int(fz['_v8_sleep_masked'].sum()):,} split nights "
            f"({n_mask:,} cells across the four tables, flag column _v8_sleep_masked). The ranking of measurements and every sleep-duration analysis run "
            f"on the full diagnostic nights (cohort_spec.FULL_NIGHTS_ONLY); the T90 hazard ratios and the PAP analyses keep every night. "
            f"Negative-control swap of the same date: the outcome columns of {list(_dm.SWAPPED_CONTROLS_2026_09_14)} are dropped, "
            f"{list(_dm.SWAPPED_CONTROLS_2026_09_14.values())} enter through X2 (disease_definitions_v8.SWAPPED_CONTROLS_2026_09_14).",
            "", "## Columns dropped", ""] + [f"- {k}: {v}" for k, v in CH.dropped.items()]
    open(f"{out}/PROVENANCE.md", "w").write("\n".join(prov) + "\n")
    md = [f"# CHANGELOG v7 -> v8 ({'SMOKE' if a.smoke else 'REAL'} build, {time.strftime('%Y-%m-%d %H:%M')})", "",
          "The v8 tables start from the v7 tables (data_frozen_v7_2026-09, read only) and change exactly the cells listed in CHANGELOG.parquet "
          "(table, BDSPPatientID, column, old, new, reason). In words:", "",
          f"1. Split-night patients ({len(split_ids_all):,} in the manifest, {len(split_ok):,} extracted): every night-level measure (the {len(meas137)} kept measures, sleep architecture, AHI and arousal index, "
          f"the six oxygen measures) now describes the untreated pre-PAP window, not the whole file (decision 1). {n_idx_undef:,} of them have their indices undefined (under 60 min of pre-PAP sleep, D1).",
          f"2. Four new measures on every cohort night (decision 3): {', '.join(NEW4)}; and the sleep-period T90 {SLEEP_T90} for the decision 2 sensitivity set. Missing on {len(coll_new_ids):,} collapsed-staging nights (D11).",
          f"3. Outcomes (decisions 5 and 11): {len(changed_keys)} conditions recomputed or new ({', '.join(changed_keys)}), {len(copied_keys)} untouched and byte-identical to v7; {ohs_note}.",
          f"4. The {len(CH.dropped['master'])} _siteZ twins dropped from master_cohort.csv and t90_base4.parquet (decision 3).",
          "5. cpap_t90_by_stage.parquet: whole-night AHI, TST_min and spo2_pct_below_90 refreshed from the v8 t90_final (as v7 did); its per-stage columns are its own.",
          f"6. v8.1 (14 Sept decision): split nights contribute oxygen and event-rate measures only; {len(S.SLEEP_MASK_COLS)} sleep-amount, composition, REM-latency and count columns set missing on {int(fz['_v8_sleep_masked'].sum()):,} split nights ({n_mask:,} cells, flag _v8_sleep_masked); the swapped-out controls' outcome columns dropped ({gone_cols}).",
          "", "## Cells by table and reason", "", by_reason.to_string(index=False), "", "## Decisions", ""] + [f"- {d}" for d in decisions]
    open(f"{out}/CHANGELOG.md", "w").write("\n".join(md) + "\n")
    for f in list(files) + ["CHANGELOG.parquet"]:
        S.write_sidecar(f"{out}/{f}", inputs, os.path.abspath(__file__), "050_build_v8_tables",
                        extra={"mode": "smoke" if a.smoke else "real", "n_rows": n_rows_v7, "cohort": len(cohort_ids), "split_ok": len(split_ok), "light_ok": len(light_ok),
                               "outputs_sha256": outs, "ohs": ohs_flag_diff}, index_path=f"{S.E2_LOGS}/PROVENANCE_INDEX.tsv")
    P(f"written -> {out} ({len(ch):,} changed cells) in {time.time()-T0:.0f}s")

    # ------------------------------------------------------------------ invariants (gates 3 and 4), stop on failure
    if a.skip_invariants:
        return 0
    import invariants_v8 as INV
    n_fail = INV.run(out, S.PREV_DIR, out, smoke=a.smoke)
    if n_fail:
        open(f"{out}/INVARIANTS_FAILED_{time.strftime('%Y%m%d_%H%M%S')}", "w").write(f"{n_fail} invariant failures, see INVARIANTS_V8.md\n")
        P(f"INVARIANTS FAILED: {n_fail} (see {out}/INVARIANTS_V8.md)")
        return 1
    P("invariants passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
