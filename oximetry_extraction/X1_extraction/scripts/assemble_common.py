"""Shared pieces of the v8 assembly steps (015 split, 023 light): jsonl loading, the D1/D11 rules, row invariants, per-site tables."""
from __future__ import annotations
import hashlib, json, os
import numpy as np, pandas as pd

STAGE_GATED_NEW = ["sol_min", "rem_latency_min", "plm_index", "lm_index", "plm_index_excl_resp", "hypoxic_burden", "hb_area_pctmin", "spo2_pct_below_90_sleep",
                   "spo2_pct_below_88_sleep", "spo2_mean_sleep", "spo2_p1_sleep", "spo2_sleep_valid_min"]
OXYGEN_FAMILY = ["spo2_mean", "spo2_nadir", "spo2_pct_below_90", "spo2_pct_below_88", "odi3_total", "odi4_total", "spo2_p1_fft", "spo2_pct_below_90_native",
                 "spo2_mean_native", "spo2_pct_below_88_native", "spo2_nadir_native", "spo2_p1_native", "odi3_total_native", "odi4_total_native"]
INDEX_COLS = ["AHI", "AHI_rows", "AHI_r4", "AHI_r4_rows", "arousal_index", "arousal_index_rows"]
NEW_MEASURES = ["hypoxic_burden", "plm_index", "lm_index", "sol_min", "rem_latency_min", "spo2_pct_below_90_sleep"]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def load_jsonl(path):
    rows = [json.loads(l) for l in open(path) if l.strip()]
    df = pd.DataFrame(rows)
    # keep the last record per night (a retry after a failure appends a second row)
    df["_ok"] = df.status.astype(str).str.startswith("ok")
    df = df.sort_values(["BDSPPatientID", "_ok"]).drop_duplicates("BDSPPatientID", keep="last").drop(columns="_ok")
    for c in df.columns:
        if df[c].dtype == object:
            conv = pd.to_numeric(df[c], errors="coerce")
            if conv.notna().sum() >= 0.5 * df[c].notna().sum() and c not in ("status", "stem", "site_id", "_encoding_used", "source"):
                if df[c].dropna().map(lambda v: isinstance(v, (int, float))).all():
                    df[c] = conv
    return df


def apply_rules(df, collapsed_ids, min_tst_for_index=60.0, min_oxygen_valid_min=5.0):
    """D1 (index undefined under 60 min of pre-PAP sleep, oxygen undefined under 5 valid minutes) on windowed rows; D11
    (stage-gated new measures NaN on collapsed-staging nights) on every row. Returns the frame with the flag columns."""
    df = df.copy()
    win = df["_v8_prepap_window"].fillna(False).astype(bool) if "_v8_prepap_window" in df else pd.Series(False, index=df.index)
    tst = pd.to_numeric(df.get("TST_min"), errors="coerce")
    df["_v8_index_undefined"] = win & (tst < min_tst_for_index)
    for c in INDEX_COLS:
        if c in df:
            df.loc[df._v8_index_undefined, c] = np.nan
    vm = pd.to_numeric(df.get("spo2_native_valid_min"), errors="coerce")
    df["_v8_oxygen_window_short"] = win & (vm < min_oxygen_valid_min)
    for c in OXYGEN_FAMILY + ["hypoxic_burden", "hb_area_pctmin", "spo2_pct_below_90_sleep"]:
        if c in df:
            df.loc[df._v8_oxygen_window_short, c] = np.nan
    df["_v7_1_stage_collapsed"] = df.BDSPPatientID.isin(collapsed_ids)
    for c in ["rem_latency_min"]:   # D11 = REM latency only (coordinator default 2026-09-12: sleep-versus-wake quantities stay on collapsed nights per the 2026-09-08 addendum)
        if c in df:
            df.loc[df._v7_1_stage_collapsed, c] = np.nan
    return df


def invariants(df):
    """Row invariants on ok rows; returns a list of (name, n_pass, n_checked, offending ids)."""
    ok = df[df.status.astype(str).str.startswith("ok")]
    out = []
    pct = ok[["N1_pct", "N2_pct", "N3_pct", "REM_pct"]].apply(pd.to_numeric, errors="coerce").sum(axis=1, min_count=4)
    m = pct.notna(); out.append(("stage percentages sum to 100 (TST>0)", int(((pct[m] - 100).abs() < 1e-6).sum()), int(m.sum()), ok.loc[m & ((pct - 100).abs() >= 1e-6), "BDSPPatientID"].tolist()[:10]))
    se = pd.to_numeric(ok.sleep_efficiency_pct, errors="coerce"); m = se.notna()
    out.append(("sleep efficiency <= 100", int((se[m] <= 100 + 1e-9).sum()), int(m.sum()), ok.loc[m & (se > 100 + 1e-9), "BDSPPatientID"].tolist()[:10]))
    tst = pd.to_numeric(ok.TST_min, errors="coerce"); rec = pd.to_numeric(ok.recording_dur_min, errors="coerce"); m = tst.notna() & rec.notna()
    out.append(("TST <= recording (window) minutes", int((tst[m] <= rec[m] + 1e-6).sum()), int(m.sum()), ok.loc[m & (tst > rec + 1e-6), "BDSPPatientID"].tolist()[:10]))
    for idx in ("AHI", "arousal_index", "plm_index"):
        if idx in ok:
            n = pd.to_numeric(ok[idx], errors="coerce") * tst / 60.0; m = n.notna()
            out.append((f"{idx} x TST/60 is an integer count", int(((n[m] - n[m].round()).abs() < 1e-6).sum()), int(m.sum()), ok.loc[m & ((n - n.round()).abs() >= 1e-6), "BDSPPatientID"].tolist()[:10]))
    hb = pd.to_numeric(ok.hypoxic_burden, errors="coerce"); nev = pd.to_numeric(ok.hb_n_events_used, errors="coerce"); m = hb.notna()
    out.append(("hypoxic burden >= 0", int((hb[m] >= 0).sum()), int(m.sum()), ok.loc[m & (hb < 0), "BDSPPatientID"].tolist()[:10]))
    out.append(("hypoxic burden > 0 only when events were used (HB may be 0 with events: no dip below the pre-event baseline)", int((~((hb[m] > 0) & (nev[m] == 0))).sum()), int(m.sum()), ok.loc[m & ((hb > 0) & (nev == 0)), "BDSPPatientID"].tolist()[:10]))
    for c, lo, hi in (("spo2_pct_below_90_native", 0, 100), ("spo2_pct_below_90_sleep", 0, 100), ("spo2_mean_native", 50, 100), ("sol_min", 0, 1e4), ("rem_latency_min", 0, 1e4), ("plm_index", 0, 1e4), ("hypoxic_burden", 0, 1e5)):
        if c in ok:
            x = pd.to_numeric(ok[c], errors="coerce"); m = x.notna()
            out.append((f"{c} within [{lo}, {hi}]", int(((x[m] >= lo) & (x[m] <= hi)).sum()), int(m.sum()), ok.loc[m & ((x < lo) | (x > hi)), "BDSPPatientID"].tolist()[:10]))
    return out


def per_site_table(df, cols, label_new="new"):
    from scipy.stats import ks_2samp
    L = []
    for c in cols:
        if c not in df:
            continue
        a = pd.to_numeric(df[df.site_id == "I0002"][c], errors="coerce").dropna(); b = pd.to_numeric(df[df.site_id == "I0006"][c], errors="coerce").dropna()
        ks = ks_2samp(a, b).statistic if len(a) > 5 and len(b) > 5 else np.nan
        L.append(f"- {c}: I0002 n={len(a)} median {a.median():.3g} IQR [{a.quantile(.25):.3g}, {a.quantile(.75):.3g}] | I0006 n={len(b)} median {b.median():.3g} IQR [{b.quantile(.25):.3g}, {b.quantile(.75):.3g}] | KS {ks:.3f}")
    return L
