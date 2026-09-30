"""
Comparison helpers: map re-extracted columns onto master_cohort_v1 / t90_final columns and score
agreement per metric.

Tolerance (from the brief): a value matches when |new - old| <= max(0.01, 0.01 * |old|); NaN
matches NaN only. Architecture columns are compared against BOTH tables where they exist, and the
master's N1/N3 (known to be swapped relative to the frozen table) are compared against the frozen
values, which the raw legend supports.
"""
from __future__ import annotations
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../..")))  # repository root, where paths.py lives
import paths

import numpy as np
import pandas as pd

SB = paths.FEATURE_WORK
T90 = paths.T90_ROOT
AUDIT = f"{paths.T90_ROOT}/EncodingB_Audit_2026-09-06"

# master column -> re-extracted column (same name unless stated)
# The four feature families keep their April names. Architecture columns are decoded with the
# file legend; the master's wake_min held the code-0 (unscored) minutes on unflagged nights.
RENAME = {
    "wake_min": "unscored_min",          # what the April extractor put in wake_min (A nights)
}
# columns whose master values are compared to the frozen table's version (N1/N3 swap fix)
FROZEN_PREFERRED = ["N1_pct", "N3_pct"]
FROZEN_COLS = ["spo2_pct_below_90", "spo2_pct_below_88", "spo2_nadir_corrected", "odi3_total", "odi4_total",
               "spo2_mean", "AHI", "arousal_index", "N2_pct", "REM_pct", "N3_pct", "sleep_efficiency_pct",
               "N1_pct", "TST_min", "recording_dur_min"]
ARCH_COLS = ["recording_dur_min", "TST_min", "wake_min", "N1_min", "N2_min", "N3_min", "REM_min",
             "sleep_efficiency_pct", "N1_pct", "N2_pct", "N3_pct", "REM_pct", "WASO_min", "arousal_index", "AHI"]
NON_NUMERIC = {"ecg_channel_used"}


def load_master(cols=None):
    m = pd.read_csv(f"{SB}/outputs/master_cohort_v1.csv", low_memory=False, float_precision="round_trip",
                    usecols=(None if cols is None else (["BDSPPatientID"] + [c for c in cols if c != "BDSPPatientID"])))
    return m.set_index("BDSPPatientID")


def load_frozen():
    f = pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet")
    return f.set_index("BDSPPatientID")


def cohort_ids():
    f = pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet",
                        columns=["BDSPPatientID", "fu_valid", "spo2_pct_below_90", "oximetry_bad"])
    c = f[(f.fu_valid == 1) & f.spo2_pct_below_90.notna() & (f.oximetry_bad == 0)]
    assert len(c) == 19173, len(c)
    return c.BDSPPatientID.tolist()


def ranked_names():
    r = pd.read_csv(f"{paths.NUMBERS_DIR}/ranking_v3.csv", comment="#")
    return r.feature.tolist()


def feature_columns(master_cols=None):
    """The 122 April feature columns + architecture columns + spo2_nadir_corrected, in master order.
    Read from the CSV header so the index column does not matter."""
    header = pd.read_csv(f"{SB}/outputs/master_cohort_v1.csv", nrows=0).columns.tolist()
    idx = header.index("BDSPPatientID")
    raw = [c for c in header[:idx] if c not in NON_NUMERIC]
    return raw + ARCH_COLS + ["spo2_nadir_corrected"]


def match(new, old, rel=0.01, abs_=0.01):
    """Elementwise match per the brief's tolerance, NaN pattern exact."""
    new = np.asarray(new, dtype=float); old = np.asarray(old, dtype=float)
    both_nan = np.isnan(new) & np.isnan(old)
    one_nan = np.isnan(new) ^ np.isnan(old)
    ok = np.zeros(new.shape, dtype=bool)
    fin = ~np.isnan(new) & ~np.isnan(old)
    ok[fin] = np.abs(new[fin] - old[fin]) <= np.maximum(abs_, rel * np.abs(old[fin]))
    ok[both_nan] = True
    return ok, one_nan, both_nan


def site_constants(master):
    """Per-site mean and sd (ddof=1) of every raw feature that has a _siteZ twin, computed exactly as the
    master's siteZ columns were (verified to 1e-14 on the master itself)."""
    zcols = [c for c in master.columns if c.endswith("_siteZ")]
    rows = []
    for z in zcols:
        raw = z[:-6]
        for site, g in master.groupby("site_id"):
            x = g[raw].astype(float)
            rows.append({"feature": raw, "site_id": site, "mean": x.mean(), "sd": x.std(ddof=1), "n": int(x.notna().sum())})
    return pd.DataFrame(rows)


def add_siteZ(df, consts):
    """Standardise re-extracted raw features with the master's own site constants."""
    out = df.copy()
    for (feat), g in consts.groupby("feature"):
        if feat not in out.columns:
            continue
        z = np.full(len(out), np.nan)
        for _, r in g.iterrows():
            m = (out["site_id"] == r.site_id).to_numpy()
            z[m] = (out.loc[m, feat].astype(float).to_numpy() - r["mean"]) / r["sd"]
        out[f"{feat}_siteZ"] = z
    return out


def compare_table(new: pd.DataFrame, master: pd.DataFrame, frozen: pd.DataFrame, cols, flag_col="_encoding_used"):
    """Per-metric agreement, split by flag. `new` indexed by BDSPPatientID with re-extracted columns."""
    ids = new.index
    flags = master.loc[ids, flag_col].fillna("none")
    rows = []
    per_night = {}
    for col in cols:
        src = RENAME.get(col, col)
        if src not in new.columns:
            rows.append({"metric": col, "status": "not_computed"}); continue
        old_m = master.loc[ids, col].astype(float) if col in master.columns else pd.Series(np.nan, index=ids)
        old = old_m.copy()
        table = "master"
        if col in FROZEN_PREFERRED and col in frozen.columns:
            fz = frozen.reindex(ids)[col].astype(float)
            old = fz.where(fz.notna(), old_m); table = "frozen"
        nv = new[src].astype(float)
        ok, one_nan, both_nan = match(nv.to_numpy(), old.to_numpy())
        d = (nv - old).to_numpy()
        rec = {"metric": col, "compared_to": table, "source_column": src, "n": len(ids),
               "n_both_present": int((~np.isnan(nv.to_numpy()) & ~np.isnan(old.to_numpy())).sum()),
               "n_nan_mismatch": int(one_nan.sum()), "n_both_nan": int(both_nan.sum()),
               "frac_match_all": float(ok.mean())}
        for fl in ("A", "B", "none"):
            m = (flags == fl).to_numpy()
            rec[f"n_{fl}"] = int(m.sum())
            rec[f"frac_match_{fl}"] = float(ok[m].mean()) if m.any() else np.nan
        fin = ~np.isnan(d)
        rec["median_abs_diff"] = float(np.median(np.abs(d[fin]))) if fin.any() else np.nan
        rec["max_abs_diff"] = float(np.max(np.abs(d[fin]))) if fin.any() else np.nan
        rec["median_rel_diff"] = float(np.nanmedian(np.abs(d[fin]) / np.maximum(1e-12, np.abs(old.to_numpy()[fin])))) if fin.any() else np.nan
        if fin.any():
            worst = np.argsort(-np.abs(np.where(fin, d, 0)))[:3]
            rec["worst_nights"] = ";".join(f"{ids[i]}:{nv.iloc[i]:.6g}vs{old.iloc[i]:.6g}" for i in worst)
        rows.append(rec)
        per_night[col] = ok
    return pd.DataFrame(rows), pd.DataFrame(per_night, index=ids)
