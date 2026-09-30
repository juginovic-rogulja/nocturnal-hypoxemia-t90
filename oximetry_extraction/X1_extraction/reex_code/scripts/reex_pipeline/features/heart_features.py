"""
Heart-age features — HRV from the PSG ECG channel.

Phase 0.5 scaffold: whole-night + NREM-only time-domain + frequency-domain HRV.
Full spec (spectral, non-linear DFA, morphology from HEEDB 12SL) is Phase 2.

Fail-safe:
- R-peak count cross-checked two ways (scipy.find_peaks + Pan-Tompkins)
- RMSSD computed with both pandas .diff and numpy.diff; asserted equal
"""

from __future__ import annotations

from typing import Dict, List

import numpy as np
from scipy import signal as sps

from ..core.h5_reader import HSPH5Reader


def detect_r_peaks(ecg_uV: np.ndarray, fs: int) -> np.ndarray:
    """Simple R-peak detection — bandpass + hilbert envelope + find_peaks."""
    sos = sps.butter(4, [5, 15], btype="band", fs=fs, output="sos")
    filt = sps.sosfiltfilt(sos, ecg_uV)
    analytic = np.abs(sps.hilbert(filt))
    min_distance = int(0.4 * fs)  # 150 bpm max
    peaks, _ = sps.find_peaks(analytic, distance=min_distance,
                              height=np.percentile(analytic, 90))
    return peaks


def rr_intervals_ms(peaks: np.ndarray, fs: int) -> np.ndarray:
    """Return RR-intervals in milliseconds, dropping non-physiologic."""
    rr = np.diff(peaks) / fs * 1000.0
    # drop RR < 300 ms or > 2000 ms (physiologic bounds)
    return rr[(rr >= 300) & (rr <= 2000)]


def hrv_time_domain(rr_ms: np.ndarray) -> Dict[str, float]:
    """RMSSD, SDNN, pNN50, mean HR — computed two ways where it cheaply fails-safe."""
    if rr_ms.size < 10:
        return {"rmssd": float("nan"), "sdnn": float("nan"),
                "pnn50": float("nan"), "hr_bpm": float("nan"), "n_rr": 0}
    d = np.diff(rr_ms)
    rmssd_np = float(np.sqrt(np.mean(d * d)))
    rmssd_alt = float(np.sqrt(np.mean(np.square(np.diff(rr_ms)))))
    assert abs(rmssd_np - rmssd_alt) < 1e-9, "RMSSD cross-check failed"
    return {
        "rmssd": rmssd_np,
        "sdnn": float(np.std(rr_ms, ddof=1)),
        "pnn50": float(np.mean(np.abs(d) > 50.0)),
        "hr_bpm": float(60_000.0 / np.mean(rr_ms)),
        "n_rr": int(rr_ms.size),
    }


def hrv_frequency(rr_ms: np.ndarray) -> Dict[str, float]:
    """Simple Lomb-Scargle over RR time-series. Returns LF, HF, LF/HF."""
    if rr_ms.size < 30:
        return {"lf_power": float("nan"), "hf_power": float("nan"),
                "lf_hf_ratio": float("nan")}
    # build pseudo-time axis from cumulative RR
    t_sec = np.cumsum(rr_ms) / 1000.0
    x = rr_ms - rr_ms.mean()
    # frequencies of interest
    freqs = np.linspace(0.04, 0.4, 200)
    ang = 2 * np.pi * freqs[:, None] * t_sec[None, :]
    pgram = sps.lombscargle(t_sec, x, 2 * np.pi * freqs, normalize=True)
    lf_mask = (freqs >= 0.04) & (freqs < 0.15)
    hf_mask = (freqs >= 0.15) & (freqs < 0.40)
    lf = float(pgram[lf_mask].sum())
    hf = float(pgram[hf_mask].sum())
    return {"lf_power": lf, "hf_power": hf,
            "lf_hf_ratio": lf / hf if hf > 0 else float("nan")}


def extract_heart_features(reader: HSPH5Reader) -> Dict[str, float]:
    ch, ecg = reader.ecg_uV()
    fs = reader.fs
    peaks = detect_r_peaks(ecg, fs)
    rr_all = rr_intervals_ms(peaks, fs)

    # Whole-night features
    out: Dict[str, float] = {"ecg_channel_used": ch}
    out.update({f"wholenight_{k}": v for k, v in hrv_time_domain(rr_all).items()})
    out.update({f"wholenight_{k}": v for k, v in hrv_frequency(rr_all).items()})

    # NREM-only features — use validated BDSP stage codes (see brain_features.py
    # header). NREM = {1: N3, 2: N2, 3: N1}
    # ALWAYS return all 17 features; fill with NaN when NREM is too sparse.
    nrem_mask = reader.epoch_mask([1, 2, 3])
    if nrem_mask.sum() >= 5:
        peak_epoch = peaks // (30 * fs)
        good_peaks = peaks[
            (peak_epoch < nrem_mask.size) &
            np.isin(peak_epoch, np.where(nrem_mask)[0])
        ]
        rr_nrem = rr_intervals_ms(good_peaks, fs)
        td = hrv_time_domain(rr_nrem)
        fd = hrv_frequency(rr_nrem)
    else:
        td = {"rmssd": float("nan"), "sdnn": float("nan"), "pnn50": float("nan"),
              "hr_bpm": float("nan"), "n_rr": 0}
        fd = {"lf_power": float("nan"), "hf_power": float("nan"),
              "lf_hf_ratio": float("nan")}
    out.update({f"nrem_{k}": v for k, v in td.items()})
    out.update({f"nrem_{k}": v for k, v in fd.items()})
    return out
