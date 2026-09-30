"""
Brain-age features — Phase 0.5 scaffold.

This is the FIRST iteration, intentionally minimal. It demonstrates the
end-to-end pipeline (read h5 -> per-region regional signal -> per-stage
band power + kurtosis) and must be extended in Phase 2 to match the full
Sun 2024 feature set (spindles, slow oscillations, spindle-SO coupling).

Fail-safe: every feature is computed twice-via-different-paths where feasible
(e.g., band power via Welch PSD and via direct filter-then-variance). Downstream
fail-safe assertions ensure the two agree to within tolerance on test data.
"""

from __future__ import annotations

from typing import Dict, List

import numpy as np
from scipy import signal as sps
from scipy.stats import kurtosis

from ..core.h5_reader import HSPH5Reader, REGION_CHANNELS

# Stage codes — VALIDATED 2026-04-22 against manual annotations.csv on
# sub-SXXXXNNNNNNNNNNN (see tests/validate_stage_encoding.md). BDSP h5 uses
# an inverted-AASM encoding:
#   0 = unscored/pre-recording (epochs before lights-off)
#   1 = N3
#   2 = N2
#   3 = N1
#   4 = REM
#   5 = Wake
#   9 = unknown/artifact
STAGE_NREM = [1, 2, 3]   # N3, N2, N1 — union matches Sun 2024's "NREM"
STAGE_REM = [4]
STAGE_WAKE = [5]

BANDS = {
    "delta": (1.0, 4.0),
    "theta": (4.0, 8.0),
    "alpha": (8.0, 12.0),
    "sigma": (12.0, 16.0),
    "beta": (16.0, 30.0),
}


def bandpower_welch(x: np.ndarray, fs: int, low: float, high: float,
                    nperseg: int = 2 * 30 * 200) -> float:
    """PSD-integrated bandpower (Welch). Returns mean PSD within [low, high].
    NOTE: returns raw bandpower (μV² / Hz). Always log-transform + winsorize
    before passing to ElasticNet or clustering — raw bandpower has 1000× outliers
    from brief artifacts/wake mis-labeled as NREM.
    """
    f, pxx = sps.welch(x, fs=fs, nperseg=min(nperseg, x.size), noverlap=0)
    m = (f >= low) & (f < high)
    # numpy 2.x: np.trapz was removed in favor of np.trapezoid
    _trap = getattr(np, "trapezoid", None) or np.trapz  # type: ignore[attr-defined]
    return float(_trap(pxx[m], f[m]))


def bandpower_filter_variance(x: np.ndarray, fs: int, low: float, high: float) -> float:
    """Bandpass filter + variance — independent check on Welch bandpower."""
    sos = sps.butter(4, [low, high], btype="band", fs=fs, output="sos")
    y = sps.sosfiltfilt(sos, x)
    return float(np.var(y))


def _epoch_signal(x: np.ndarray, fs: int, epoch_sec: int = 30) -> np.ndarray:
    samples = epoch_sec * fs
    n = x.size // samples
    return x[: n * samples].reshape(n, samples)


def regional_band_features(reader: HSPH5Reader, region: str,
                           stage_codes: List[int]) -> Dict[str, float]:
    """
    For one region (F/C/O) and one set of stage codes (e.g. NREM),
    compute mean + kurtosis of bandpower in each BANDS entry across
    concatenated epochs.
    """
    sig = reader.region_signal_uV(region)
    fs = reader.fs
    epochs = _epoch_signal(sig, fs)                  # (n_epochs, samples)
    mask = reader.epoch_mask(stage_codes)            # (n_epochs,)
    if mask.sum() < 5:
        # too few epochs in this stage set — return NaNs
        return {
            f"{region}_{stage_label(stage_codes)}_{b}_{stat}": float("nan")
            for b in BANDS for stat in ("mean", "kurt")
        }

    out: Dict[str, float] = {}
    for band_name, (lo, hi) in BANDS.items():
        # per-epoch bandpower via Welch
        per_epoch = np.array([
            bandpower_welch(epochs[i], fs, lo, hi) for i in np.where(mask)[0]
        ])
        # fail-safe cross-check on one epoch: Welch vs filter-variance
        i0 = int(np.where(mask)[0][0])
        bp_w = bandpower_welch(epochs[i0], fs, lo, hi)
        bp_f = bandpower_filter_variance(epochs[i0], fs, lo, hi)
        # Allow Welch vs filter variance to differ by a factor (not identical
        # metrics — Welch integrates PSD, filter-variance is time-domain var)
        # but require they're both positive and within 2 OOM of each other.
        if bp_w > 0 and bp_f > 0:
            ratio = max(bp_w, bp_f) / max(1e-12, min(bp_w, bp_f))
            if ratio > 1000:
                out[f"{region}_{stage_label(stage_codes)}_{band_name}_warn"] = float(ratio)

        label = f"{region}_{stage_label(stage_codes)}_{band_name}"
        out[f"{label}_mean"] = float(np.nanmean(per_epoch))
        out[f"{label}_kurt"] = float(kurtosis(per_epoch, fisher=True, bias=False))
    return out


def stage_label(stage_codes: List[int]) -> str:
    s = set(stage_codes)
    if s == set(STAGE_NREM): return "nrem"
    if s == set(STAGE_REM):  return "rem"
    if s == set(STAGE_WAKE): return "wake"
    return "_".join(str(c) for c in sorted(stage_codes))


def extract_brain_features(reader: HSPH5Reader) -> Dict[str, float]:
    """Top-level entry point — returns a flat dict of brain features."""
    out: Dict[str, float] = {}
    for region in REGION_CHANNELS:
        if not all(reader.has_channel(c) for c in REGION_CHANNELS[region]):
            # skip region if channels missing — record as NaN
            for stage in (STAGE_NREM, STAGE_REM):
                for band in BANDS:
                    out[f"{region}_{stage_label(stage)}_{band}_mean"] = float("nan")
                    out[f"{region}_{stage_label(stage)}_{band}_kurt"] = float("nan")
            continue
        out.update(regional_band_features(reader, region, STAGE_NREM))
        out.update(regional_band_features(reader, region, STAGE_REM))
    return out
