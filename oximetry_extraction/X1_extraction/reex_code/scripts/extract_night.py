"""
Re-extract every reproducible master_cohort_v1 column for ONE night from a local .h5.

The four feature families (brain, microstructure, heart, lung) are the ORIGINAL April 2026
functions, imported unchanged from reex_pipeline.features (byte-identical copies of the
feature_pipeline tarball, only np.trapz -> np.trapezoid). Only the h5 reader differs, to bridge
the August 2026 file layout (see reex_pipeline/core/h5_reader.py).

Architecture columns follow arch_worker.py's definitions but with the file's own stage legend
(1 N3, 2 N2, 3 N1, 4 REM, 5 Wake, 0 unscored) and, as the April extractor did, minutes counted at
sample level (samples / fs / 60). Epoch-mode variants are kept with suffix _ep.

Usage: python extract_night.py <h5 path> [--json out.json]
"""
from __future__ import annotations

import json
import sys
import time
import traceback
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from reex_pipeline.core.h5_reader import (HSPH5Reader, REGION_CHANNELS, ECG_CHANNEL_PRIORITY,
                                          has_usable_staging_from_epochs)
from reex_pipeline.features.brain_features import extract_brain_features
from reex_pipeline.features.heart_features import extract_heart_features, detect_r_peaks, rr_intervals_ms
from reex_pipeline.features.lung_features import extract_lung_features, spo2_summary, count_desats
from reex_pipeline.features.microstructure import extract_microstructure

# file legend (identical on every file opened so far, re-read from each file and asserted)
LEGEND = {"0": "Unknown/Unscored", "1": "N3", "2": "N2", "3": "N1", "4": "REM", "5": "Wake"}
CODE_N3, CODE_N2, CODE_N1, CODE_REM, CODE_WAKE, CODE_UNSCORED = 1, 2, 3, 4, 5, 0
SLEEP_CODES = (1, 2, 3, 4)


def _sig_stats(x: np.ndarray) -> dict:
    ax = np.abs(x[np.isfinite(x)])
    if ax.size == 0:
        return {"sd": np.nan, "p50_abs": np.nan, "p999_abs": np.nan, "max_abs": np.nan}
    return {"sd": float(np.std(x)), "p50_abs": float(np.percentile(ax, 50)),
            "p999_abs": float(np.percentile(ax, 99.9)), "max_abs": float(ax.max())}


def architecture(r: HSPH5Reader) -> dict:
    fs = r.fs
    stage = r._stage_array()                     # sample level, reader fs
    n = stage.size
    rec_min = n / fs / 60.0
    out = {"recording_dur_min": rec_min, "fs_hz": fs}
    # sample-level minutes per legend code (the April extractor's accounting, correct legend)
    mins = {c: float(np.sum(stage == c)) / fs / 60.0 for c in range(0, 10)}
    out.update({
        "N3_min": mins[CODE_N3], "N2_min": mins[CODE_N2], "N1_min": mins[CODE_N1],
        "REM_min": mins[CODE_REM], "wake_min": mins[CODE_WAKE], "unscored_min": mins[CODE_UNSCORED],
    })
    tst = out["N1_min"] + out["N2_min"] + out["N3_min"] + out["REM_min"]
    out["TST_min"] = tst
    out["sleep_efficiency_pct"] = 100.0 * tst / rec_min if rec_min > 0 else np.nan
    for k in ("N1", "N2", "N3", "REM"):
        out[f"{k}_pct"] = 100.0 * out[f"{k}_min"] / tst if tst > 0 else np.nan
    sleep_samples = np.where(np.isin(stage, SLEEP_CODES))[0]
    if sleep_samples.size:
        onset = sleep_samples[0]
        out["WASO_min"] = float(np.sum(stage[onset:] == CODE_WAKE)) / fs / 60.0
        out["sol_min"] = onset / fs / 60.0
        rem = np.where(stage[onset:] == CODE_REM)[0]
        out["rem_latency_min"] = rem[0] / fs / 60.0 if rem.size else np.nan
    else:
        out["WASO_min"] = np.nan; out["sol_min"] = np.nan; out["rem_latency_min"] = np.nan
    # epoch-mode variants (arch_worker.py accounting, correct legend)
    ep = r.stage_per_epoch()
    out["n_epochs"] = int(ep.size)
    out["n_scored_epochs"] = int(np.isin(ep, (1, 2, 3, 4, 5)).sum())
    n_sleep = int(np.isin(ep, SLEEP_CODES).sum())
    out["TST_min_ep"] = n_sleep * 0.5
    out["tib_min_ep"] = ep.size * 0.5
    for k, c in (("N1", CODE_N1), ("N2", CODE_N2), ("N3", CODE_N3), ("REM", CODE_REM)):
        out[f"{k}_pct_ep"] = 100.0 * float(np.sum(ep == c)) / n_sleep if n_sleep else np.nan
    out["wake_min_ep"] = float(np.sum(ep == CODE_WAKE)) * 0.5
    out["usable_staging_april_rule"] = bool(has_usable_staging_from_epochs(ep))
    hist = {int(k): int(v) for k, v in zip(*np.unique(ep, return_counts=True))}
    out["stage_epoch_hist"] = json.dumps(hist)
    # events
    tst_h = tst / 60.0
    for kind in ("arousal", "resp", "resp_3", "resp_4", "limb"):
        if r.has_table(kind) or (kind == "resp" and r.has_table("resp_3")):
            try:
                cnt = r.event_counts(kind)
            except KeyError:
                continue
            out[f"{kind}_rows"] = cnt["rows"]
            out[f"{kind}_rows_nonzero"] = cnt["rows_nonzero"]
            out[f"{kind}_rising"] = cnt["rising_edges"]
            out[f"{kind}_rising_incl_first"] = cnt["rising_edges_incl_first"]
    out["has_resp_table"] = bool(r.has_table("resp"))
    out["has_resp3_table"] = bool(r.has_table("resp_3"))
    # indices over true sleep hours (rising-edge counts, the count that reproduces the
    # I0006 A-night arousal index; row-count versions kept beside them)
    def idx(num):
        return float(num) / tst_h if (tst_h > 0 and num is not None and np.isfinite(num)) else np.nan
    out["arousal_index"] = idx(out.get("arousal_rising", np.nan))
    out["arousal_index_rows"] = idx(out.get("arousal_rows", np.nan))
    resp_key = "resp" if r.has_table("resp") else "resp_3"
    out["AHI"] = idx(out.get(f"{resp_key}_rising", np.nan))
    out["AHI_rows"] = idx(out.get(f"{resp_key}_rows", np.nan))
    out["AHI_source_table"] = resp_key
    if r.has_table("resp_4"):
        out["AHI_r4"] = idx(out.get("resp_4_rising", np.nan))
        out["AHI_r4_rows"] = idx(out.get("resp_4_rows", np.nan))
    return out


def lung_variants(r: HSPH5Reader) -> dict:
    """Extra oxygen columns: the nadir the April code produced (spo2_nadir) at reader fs, the same
    quantities at the channel's native rate (no upsampling), and resp counts per table."""
    out = {}
    if r.has_channel("spo2"):
        info = r.channel_info("spo2")
        out["spo2_native_fs"] = info["fs"]; out["spo2_unit"] = info["unit"]
        x = r.read_raw("spo2")
        if np.nanmedian(x[x > 0]) < 2.0:
            x = x * 100
        s = spo2_summary(x)
        out["spo2_nadir"] = s.get("spo2_nadir", np.nan)
        # candidate for the master's spo2_nadir_corrected: 1st percentile of the April-style trace
        # within [50, 100] (matched it exactly on 9 of 12 local nights)
        xv = x[(x >= 50.0) & (x <= 100.0)]
        out["spo2_p1_fft"] = float(np.percentile(xv, 1)) if xv.size else np.nan
        out["spo2_frac_below_50"] = float(np.mean(x < 50.0)) if x.size else np.nan
        # native-rate versions
        raw = r._f["signals"]["spo2"][:].ravel().astype(np.float64) * info["gain"] + info["offset"]
        raw = raw.astype(np.float32)
        if np.nanmedian(raw[raw > 0]) < 2.0:
            raw = raw * 100
        sn = spo2_summary(raw)
        out["spo2_pct_below_90_native"] = sn.get("spo2_pct_below_90", np.nan)
        out["spo2_mean_native"] = sn.get("spo2_mean", np.nan)
        nfs = int(round(info["fs"])) if np.isfinite(info["fs"]) and info["fs"] >= 1 else None
        if nfs:
            out["odi3_total_native"] = float(count_desats(raw, nfs, 3.0))
            out["odi4_total_native"] = float(count_desats(raw, nfs, 4.0))
    for kind in ("resp_3", "resp_4"):
        if r.has_table(kind):
            arr = r.annotations_array(kind)
            out[f"resp_events_any_label_{kind}"] = int(np.sum(np.diff((arr != 0).astype(np.int8)) == 1))
    return out


def diagnostics(r: HSPH5Reader) -> dict:
    out = {"duration_sec_attr": float(r._f.attrs.get("duration_sec", np.nan)),
           "n_samples": int(r.n_samples), "legend": json.dumps(r.stage_event_map(), sort_keys=True),
           "legend_matches": (r.stage_event_map() == LEGEND),
           "hypopnea_rule_source": str(r._f.attrs.get("hypopnea_rule_source", "")),
           "ann_groups": ",".join(sorted(r._f["annotations/expert_1"].keys()))}
    for region, (l, rr) in REGION_CHANNELS.items():
        for ch in (l, rr):
            if r.has_channel(ch):
                info = r.channel_info(ch)
                out[f"ch_{ch}_fs"] = info["fs"]; out[f"ch_{ch}_unit"] = info["unit"]
                out[f"ch_{ch}_gain"] = info.get("gain", np.nan); out[f"ch_{ch}_physmax"] = info.get("phys_max", np.nan)
                out[f"ch_{ch}_type"] = info["type"]
            else:
                out[f"ch_{ch}_fs"] = np.nan
        if r.has_channel(l) and r.has_channel(rr):
            st = _sig_stats(r.region_signal_uV(region))
            out.update({f"amp_{region}_{k}": v for k, v in st.items()})
    for c in ECG_CHANNEL_PRIORITY:
        if r.has_channel(c):
            info = r.channel_info(c)
            out["ecg_fs"] = info["fs"]; out["ecg_unit"] = info["unit"]; out["ecg_gain"] = info.get("gain", np.nan)
            ecg = r.read_channel_uV(c)
            out.update({f"amp_ecg_{k}": v for k, v in _sig_stats(ecg).items()})
            pk = detect_r_peaks(ecg, r.fs)
            rr_ms = np.diff(pk) / r.fs * 1000.0
            out["ecg_n_peaks"] = int(pk.size)
            out["ecg_rr_median_ms_all"] = float(np.median(rr_ms)) if rr_ms.size else np.nan
            out["ecg_rr_frac_kept_300_2000"] = float(np.mean((rr_ms >= 300) & (rr_ms <= 2000))) if rr_ms.size else np.nan
            out["ecg_rr_frac_ge_1900"] = float(np.mean(rr_ms >= 1900)) if rr_ms.size else np.nan
            break
    out["channels"] = ",".join(r.channels())
    return out


def extract(path: str | Path) -> dict:
    t0 = time.time()
    rec = {"source": str(path), "status": "ok"}
    try:
        with HSPH5Reader(path) as r:
            rec.update(diagnostics(r))
            # the four April families, unchanged functions
            brain = extract_brain_features(r)
            micro = extract_microstructure(r)
            heart = extract_heart_features(r)
            lung = extract_lung_features(r)
            rec["n_brain"], rec["n_micro"], rec["n_heart"], rec["n_lung"] = len(brain), len(micro), len(heart), len(lung)
            rec.update(brain); rec.update(micro); rec.update(heart); rec.update(lung)
            rec.update(lung_variants(r))
            rec.update(architecture(r))
            rec.update({f"diag_{k}": (json.dumps(v) if isinstance(v, (list, dict)) else v) for k, v in r.diag.items()})
    except Exception as e:
        rec["status"] = f"fail:{type(e).__name__}:{str(e)[:200]}"
        rec["traceback"] = traceback.format_exc()[-1500:]
    rec["runtime_sec"] = round(time.time() - t0, 2)
    return rec


if __name__ == "__main__":
    p = sys.argv[1]
    rec = extract(p)
    if "--json" in sys.argv:
        Path(sys.argv[sys.argv.index("--json") + 1]).write_text(json.dumps(rec, indent=1, default=str))
    keys = ["status", "runtime_sec", "recording_dur_min", "TST_min", "N1_min", "N2_min", "N3_min", "REM_min",
            "wake_min", "unscored_min", "sleep_efficiency_pct", "WASO_min", "arousal_index", "AHI",
            "spo2_pct_below_90", "spo2_mean", "odi3_total", "odi4_total", "resp_events_any_label",
            "F_nrem_delta_mean", "C_nrem_delta_mean", "wholenight_rmssd", "F_spindle_mean_amp_uV"]
    for k in keys:
        print(f"{k:28s} {rec.get(k)}")
