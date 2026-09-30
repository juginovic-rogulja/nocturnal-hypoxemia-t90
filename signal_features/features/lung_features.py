"""
Lung-age features — respiratory + SpO2.

Phase 0.5 scaffold: SpO2 summary, basic desaturation counts, AHI proxy.
Full spec (hypoxic burden via Azarbarzin method, position-stratified AHI,
respiratory rate variability) is Phase 2.
"""

from __future__ import annotations

from typing import Dict

import numpy as np

from ..core.h5_reader import HSPH5Reader


def spo2_summary(spo2: np.ndarray) -> Dict[str, float]:
    """Basic SpO2 summary — mean, nadir, time-below-90%, ODI3, ODI4."""
    spo2 = spo2.astype(np.float32)
    # drop impossible SpO2 values
    spo2 = spo2[(spo2 >= 50.0) & (spo2 <= 100.0)]
    if spo2.size == 0:
        return {k: float("nan") for k in ("spo2_mean", "spo2_nadir",
                                          "spo2_pct_below_90", "odi_events")}
    return {
        "spo2_mean": float(np.mean(spo2)),
        "spo2_nadir": float(np.min(spo2)),
        "spo2_pct_below_90": float(np.mean(spo2 < 90.0) * 100),
        "spo2_pct_below_88": float(np.mean(spo2 < 88.0) * 100),
    }


def count_desats(spo2: np.ndarray, fs: int, drop: float = 3.0) -> int:
    """Count episodes where SpO2 drops by >= `drop` % within a 30-sec window."""
    if spo2.size == 0:
        return 0
    window = 30 * fs
    events = 0
    i = 0
    while i < spo2.size - window:
        seg = spo2[i : i + window]
        seg = seg[(seg >= 50.0) & (seg <= 100.0)]
        if seg.size >= 0.5 * window:
            if seg.max() - seg.min() >= drop:
                events += 1
                i += window  # skip forward to avoid double-counting
                continue
        i += fs  # advance 1 sec
    return events


def extract_lung_features(reader: HSPH5Reader) -> Dict[str, float]:
    fs = reader.fs
    out: Dict[str, float] = {}
    if reader.has_channel("spo2"):
        spo2 = reader.read_raw("spo2")
        # BDSP spo2 sometimes stored raw as percentage (50–100), sometimes
        # scaled 0.5–1.0. Auto-detect and normalize to % scale.
        if np.nanmedian(spo2[spo2 > 0]) < 2.0:
            spo2 = spo2 * 100
        out.update(spo2_summary(spo2))
        out["odi3_total"] = float(count_desats(spo2, fs, drop=3.0))
        out["odi4_total"] = float(count_desats(spo2, fs, drop=4.0))
    else:
        for k in ("spo2_mean", "spo2_nadir", "spo2_pct_below_90",
                  "spo2_pct_below_88", "odi3_total", "odi4_total"):
            out[k] = float("nan")

    # Respiratory event count from annotations (non-zero markers = events)
    try:
        resp = reader.annotations_array("resp")
        out["resp_events_any_label"] = int(np.sum(np.diff((resp != 0).astype(np.int8)) == 1))
    except KeyError:
        out["resp_events_any_label"] = -1

    return out
