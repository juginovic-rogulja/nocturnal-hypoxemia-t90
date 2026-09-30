"""
Sleep microstructure features — spindles, slow oscillations, and spindle-SO coupling.

Phase 2 targets per SAP:
  - Spindle density, amplitude, frequency, duration (per region, per stage)
  - Slow oscillation density, amplitude, negative slope (per region)
  - Spindle-SO coupling: angle, overlap fraction, modulation strength

These are the features that dominated Sun 2024's brain-age clock and are
named as the expected top predictors in the analysis plan (SAP section 6.1).

Phase 0.5 scaffold:
  * Spindle detection via Moelle-style band-pass RMS threshold (10-16 Hz,
    z-score > 2 for > 0.5 sec).
  * SO detection via low-pass (0.5-1.25 Hz) zero-crossing detection with
    amplitude threshold.
  * Coupling: for each spindle peak, find nearest SO phase via Hilbert
    transform; report circular-mean angle + modulation strength (r).

All computed per region (F/C/O) during NREM-only.
"""

from __future__ import annotations

from typing import Dict

import numpy as np
from scipy import signal as sps

from ..core.h5_reader import HSPH5Reader, REGION_CHANNELS

STAGE_NREM = [1, 2, 3]  # validated 2026-04-22; see brain_features.py header


def detect_spindles(sig_uV: np.ndarray, fs: int, nrem_mask_sample: np.ndarray) -> Dict:
    """
    Moelle-like spindle detector:
      1. Bandpass 10-16 Hz
      2. Hilbert envelope
      3. z-score within NREM
      4. Event = z > 2 for >= 0.5 sec, <= 3 sec
    """
    sos = sps.butter(4, [10, 16], btype="band", fs=fs, output="sos")
    filt = sps.sosfiltfilt(sos, sig_uV)
    env = np.abs(sps.hilbert(filt))

    # z-score over NREM samples only
    nrem_env = env[nrem_mask_sample]
    if nrem_env.size < fs * 60:  # need at least 1 min of NREM
        return {"density_per_min": float("nan"), "mean_amp_uV": float("nan"),
                "mean_dur_sec": float("nan"), "mean_freq_Hz": float("nan"),
                "n_spindles": 0}
    mu, sd = nrem_env.mean(), nrem_env.std()
    if sd == 0:
        return {"density_per_min": float("nan"), "mean_amp_uV": float("nan"),
                "mean_dur_sec": float("nan"), "mean_freq_Hz": float("nan"),
                "n_spindles": 0}
    z = (env - mu) / sd

    # Threshold crossings within NREM only
    over = (z > 2.0) & nrem_mask_sample
    # Group contiguous samples into events
    diffs = np.diff(over.astype(np.int8))
    starts = np.where(diffs == 1)[0] + 1
    stops = np.where(diffs == -1)[0] + 1
    if over[0]: starts = np.r_[0, starts]
    if over[-1]: stops = np.r_[stops, over.size]

    durations = (stops - starts) / fs
    valid = (durations >= 0.5) & (durations <= 3.0)
    starts, stops, durations = starts[valid], stops[valid], durations[valid]

    if starts.size == 0:
        return {"density_per_min": 0.0, "mean_amp_uV": float("nan"),
                "mean_dur_sec": float("nan"), "mean_freq_Hz": float("nan"),
                "n_spindles": 0}

    # Peak amplitude per spindle
    amps = np.array([np.max(env[s:e]) for s, e in zip(starts, stops)])
    # Dominant frequency within spindle window
    freqs = []
    for s, e in zip(starts, stops):
        seg = filt[s:e]
        if seg.size < 8: continue
        fft_freqs = np.fft.rfftfreq(seg.size, d=1/fs)
        fft_pow = np.abs(np.fft.rfft(seg)) ** 2
        band = (fft_freqs >= 10) & (fft_freqs <= 16)
        if band.any():
            freqs.append(float(fft_freqs[band][np.argmax(fft_pow[band])]))

    nrem_minutes = nrem_mask_sample.sum() / fs / 60
    return {
        "density_per_min": starts.size / max(nrem_minutes, 1e-9),
        "mean_amp_uV": float(np.mean(amps)),
        "mean_dur_sec": float(np.mean(durations)),
        "mean_freq_Hz": float(np.mean(freqs)) if freqs else float("nan"),
        "n_spindles": int(starts.size),
        "_peak_samples": [int((s + e) // 2) for s, e in zip(starts, stops)],
    }


def detect_slow_oscillations(sig_uV: np.ndarray, fs: int,
                              nrem_mask_sample: np.ndarray) -> Dict:
    """
    SO detector (Staresina-style):
      1. Low-pass 1.25 Hz
      2. Zero-crossings (pos->neg transitions)
      3. Keep waves with period 0.8-2.0 sec and neg-peak-to-pos-peak amplitude > 75 uV
    """
    sos = sps.butter(4, [0.5, 1.25], btype="band", fs=fs, output="sos")
    filt = sps.sosfiltfilt(sos, sig_uV)

    # Find negative peaks (SO trough)
    neg_peaks, _ = sps.find_peaks(-filt, distance=int(0.8 * fs))
    # Keep only within NREM
    neg_peaks = neg_peaks[nrem_mask_sample[neg_peaks]]

    valid_peaks = []
    amps, durations, neg_slopes = [], [], []
    for p in neg_peaks:
        # window [-1, +1] sec around neg peak
        w0, w1 = max(0, p - fs), min(filt.size, p + fs)
        seg = filt[w0:w1]
        if seg.size < 2 * fs * 0.8: continue
        # Find the following positive peak within 2 sec
        after = seg[p - w0:]
        if after.size < int(0.4 * fs): continue
        pos_rel = np.argmax(after)
        if pos_rel == 0 or pos_rel > int(1.0 * fs): continue
        amp = after[pos_rel] - seg[p - w0]  # neg-peak-to-pos-peak
        if amp < 75.0: continue
        dur = pos_rel / fs
        slope = (after[pos_rel] - seg[p - w0]) / dur  # uV/sec
        valid_peaks.append(p)
        amps.append(amp)
        durations.append(dur)
        neg_slopes.append(slope)

    nrem_minutes = nrem_mask_sample.sum() / fs / 60
    if not valid_peaks:
        return {"density_per_min": 0.0, "mean_amp_uV": float("nan"),
                "mean_dur_sec": float("nan"), "mean_slope_uV_per_sec": float("nan"),
                "n_SOs": 0, "_peak_samples": []}
    return {
        "density_per_min": len(valid_peaks) / max(nrem_minutes, 1e-9),
        "mean_amp_uV": float(np.mean(amps)),
        "mean_dur_sec": float(np.mean(durations)),
        "mean_slope_uV_per_sec": float(np.mean(neg_slopes)),
        "n_SOs": len(valid_peaks),
        "_peak_samples": valid_peaks,
    }


def spindle_so_coupling(sig_uV: np.ndarray, fs: int,
                        spindle_peaks: list, so_peaks: list) -> Dict:
    """
    For each spindle, compute the SO phase at the spindle peak (if within 1 sec).
    Return circular mean (angle) and modulation strength (r = sqrt(mean(cos)^2 + mean(sin)^2)).
    """
    if not spindle_peaks or not so_peaks:
        return {"coupling_angle_rad": float("nan"),
                "coupling_strength_r": float("nan"),
                "coupling_overlap_frac": float("nan")}

    # SO phase via Hilbert on 0.5-1.25 Hz filtered signal
    sos = sps.butter(4, [0.5, 1.25], btype="band", fs=fs, output="sos")
    filt = sps.sosfiltfilt(sos, sig_uV)
    phase = np.angle(sps.hilbert(filt))

    so_peaks = np.array(so_peaks)
    overlaps = []
    for sp in spindle_peaks:
        distances = np.abs(so_peaks - sp)
        if distances.min() <= fs:  # within 1 sec of any SO
            overlaps.append(phase[sp])

    if not overlaps:
        return {"coupling_angle_rad": float("nan"),
                "coupling_strength_r": float("nan"),
                "coupling_overlap_frac": 0.0}

    angles = np.array(overlaps)
    r = np.sqrt(np.mean(np.cos(angles)) ** 2 + np.mean(np.sin(angles)) ** 2)
    theta = np.arctan2(np.mean(np.sin(angles)), np.mean(np.cos(angles)))
    return {
        "coupling_angle_rad": float(theta),
        "coupling_strength_r": float(r),
        "coupling_overlap_frac": len(overlaps) / len(spindle_peaks),
    }


def extract_microstructure(reader: HSPH5Reader) -> Dict[str, float]:
    fs = reader.fs
    nrem_epoch_mask = reader.epoch_mask(STAGE_NREM)
    samples_per_epoch = 30 * fs
    nrem_sample_mask = np.repeat(nrem_epoch_mask, samples_per_epoch)

    out: Dict[str, float] = {}
    for region in REGION_CHANNELS:
        if not all(reader.has_channel(c) for c in REGION_CHANNELS[region]):
            continue
        sig = reader.region_signal_uV(region)
        # pad/trim sample mask to signal length
        nsm = nrem_sample_mask[: sig.size]
        if nsm.size < sig.size:
            nsm = np.pad(nsm, (0, sig.size - nsm.size), constant_values=False)

        sp = detect_spindles(sig, fs, nsm)
        so = detect_slow_oscillations(sig, fs, nsm)
        cpl = spindle_so_coupling(sig, fs, sp.get("_peak_samples", []),
                                  so.get("_peak_samples", []))

        for k, v in sp.items():
            if not k.startswith("_"):
                out[f"{region}_spindle_{k}"] = v
        for k, v in so.items():
            if not k.startswith("_"):
                out[f"{region}_SO_{k}"] = v
        for k, v in cpl.items():
            out[f"{region}_{k}"] = v
    return out
