"""
HSP H5 reader — single source of truth for reading BDSP preprocessed HDF5 files.

BDSP preprocessed .h5 structure (discovered 2026-04-22 on sub-SXXXXNNNNNNNNNNN):
  attrs:
    meas_date: ISO timestamp
    sampling_rate: 200  (signals & annotations; always 200 Hz)
    unit_voltage: 'V'   (signals are in volts -> multiply by 1e6 for uV)
  /signals/<name>        float32 (N, 1) at 200 Hz
  /annotations/stage     float32 (N, 1) sample-level stage labels
  /annotations/arousal   float32 (N, 1) arousal events
  /annotations/resp      float32 (N, 1) respiratory events
  /annotations/limb      float32 (N, 1) limb movement events

Stage encoding (inferred, 2026-04-22 — will be validated against manual
annotations in a separate step before Phase 2):
  0 = Wake
  1 = N1
  2 = N2
  3 = N3
  4 = ambiguous (appears in data, not documented — treat as unknown)
  5 = REM
  9 = artifact / unscorable

This reader DOES NOT assume the stage mapping. Stage codes are returned as
floats; downstream code must apply a validated STAGE_MAP.
"""

from __future__ import annotations
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../..")))  # repository root, where paths.py lives
import paths

from pathlib import Path
from typing import Dict, List, Tuple

import h5py
import numpy as np


# Region aggregation for brain features (matches Sun 2024).
# Values are (left, right) channel pairs — code averages the pair into a
# single regional signal. Missing channels are skipped (returns NaN-marked).
REGION_CHANNELS: Dict[str, Tuple[str, str]] = {
    "F": ("f3-m2", "f4-m1"),
    "C": ("c3-m2", "c4-m1"),
    "O": ("o1-m2", "o2-m1"),
}

ECG_CHANNEL_PRIORITY: List[str] = ["ecg", "ecg-v1", "ecg-v2", "ecg-ll", "ecg-ra"]


class HSPH5Reader:
    """Read a single BDSP .h5 file. Lazy; opens on context-manager use."""

    def __init__(self, path: str | Path):
        self.path = Path(path)
        self._f: h5py.File | None = None

    def __enter__(self) -> "HSPH5Reader":
        self._f = h5py.File(self.path, "r")
        return self

    def __exit__(self, *exc):
        if self._f is not None:
            self._f.close()
            self._f = None

    # ------------------------------------------------------------------ attrs
    @property
    def fs(self) -> int:
        return int(self._f.attrs["sampling_rate"])

    @property
    def duration_sec(self) -> float:
        return self._stage_array().size / self.fs

    @property
    def meas_date(self) -> str:
        return str(self._f.attrs["meas_date"])

    # -------------------------------------------------------------- channels
    def channels(self) -> List[str]:
        return sorted(self._f["signals"].keys())

    def has_channel(self, name: str) -> bool:
        return name.lower() in self._f["signals"]

    def read_channel_uV(self, name: str) -> np.ndarray:
        """Return signal in microvolts (float32, 1-D)."""
        raw = self._f[f"signals/{name.lower()}"][:].flatten().astype(np.float32)
        return raw * 1e6  # V -> uV

    def read_raw(self, name: str) -> np.ndarray:
        """Return raw signal in native units (Volts or whatever)."""
        return self._f[f"signals/{name.lower()}"][:].flatten().astype(np.float32)

    # ---------------------------------------------------------- aggregations
    def region_signal_uV(self, region: str) -> np.ndarray:
        """Average of the left/right channel pair for a region (F/C/O). Returns uV."""
        left, right = REGION_CHANNELS[region]
        if not (self.has_channel(left) and self.has_channel(right)):
            raise KeyError(f"Missing channel for region {region}: need {left} and {right}")
        return 0.5 * (self.read_channel_uV(left) + self.read_channel_uV(right))

    def ecg_uV(self) -> Tuple[str, np.ndarray]:
        """Return (channel_name_used, signal_uV) from priority list."""
        for c in ECG_CHANNEL_PRIORITY:
            if self.has_channel(c):
                return c, self.read_channel_uV(c)
        raise KeyError("No ECG channel present")

    # ------------------------------------------------------------- stages
    def _stage_array(self) -> np.ndarray:
        return self._f["annotations/stage"][:].flatten().astype(np.float32)

    def stage_per_epoch(self, epoch_sec: int = 30) -> np.ndarray:
        """
        Downsample sample-level stage to per-epoch stage (one integer per epoch,
        via mode). Returns 1-D int array of length n_epochs.
        """
        stage = self._stage_array()
        samples_per_epoch = epoch_sec * self.fs
        n_epochs = stage.size // samples_per_epoch
        stage = stage[: n_epochs * samples_per_epoch].reshape(n_epochs, samples_per_epoch)
        # mode per row — use scipy to avoid numpy-bincount-fail on NaN
        out = np.zeros(n_epochs, dtype=np.int16)
        for i in range(n_epochs):
            vals, cts = np.unique(stage[i][~np.isnan(stage[i])], return_counts=True)
            out[i] = int(vals[cts.argmax()]) if vals.size else -1
        return out

    def epoch_mask(self, stage_codes: List[int], epoch_sec: int = 30) -> np.ndarray:
        """Boolean array of length n_epochs — True where stage matches."""
        st = self.stage_per_epoch(epoch_sec)
        return np.isin(st, stage_codes)

    # ------------------------------------------------------------- events
    def annotations_array(self, kind: str) -> np.ndarray:
        """Return sample-level annotation array for 'arousal', 'resp', or 'limb'."""
        return self._f[f"annotations/{kind}"][:].flatten().astype(np.float32)


def probe_h5(path: str | Path) -> Dict:
    """One-shot diagnostic — returns a dict of what's inside. Safe to call before analysis."""
    path = Path(path)
    with HSPH5Reader(path) as r:
        info = {
            "path": str(path),
            "size_mb": path.stat().st_size / 1e6,
            "fs": r.fs,
            "duration_hr": r.duration_sec / 3600,
            "n_channels": len(r.channels()),
            "channels": r.channels(),
            "has_ecg": any(r.has_channel(c) for c in ECG_CHANNEL_PRIORITY),
            "has_F_region": all(r.has_channel(c) for c in REGION_CHANNELS["F"]),
            "has_C_region": all(r.has_channel(c) for c in REGION_CHANNELS["C"]),
            "has_O_region": all(r.has_channel(c) for c in REGION_CHANNELS["O"]),
        }
        st = r.stage_per_epoch()
        vals, cts = np.unique(st, return_counts=True)
        info["stage_histogram_epochs"] = dict(zip(vals.tolist(), cts.tolist()))
    return info


def has_usable_staging(path: str | Path, min_scored_epochs: int = 200) -> bool:
    """
    Sanity check: does this h5 have usable AASM staging?
    - Returns False if >95% of epochs are stage=0 (unscored) or stage=9 (artifact)
    - Returns False if fewer than `min_scored_epochs` are scored (NREM+REM+Wake)
    Use this BEFORE feature extraction to skip patients that would waste compute.
    """
    with HSPH5Reader(path) as r:
        st = r.stage_per_epoch()
    unscored = int(((st == 0) | (st == 9)).sum())
    total = int(st.size)
    scored = total - unscored
    if total == 0:
        return False
    if scored < min_scored_epochs:
        return False
    if unscored / total > 0.95:
        return False
    return True


if __name__ == "__main__":
    import json, sys
    p = sys.argv[1] if len(sys.argv) > 1 else (
        f"{paths.PSG_DIR}/example/ses-1/eeg/"
        "sub-SXXXXNNNNNNNNNNN_ses-1.h5"
    )
    print(json.dumps(probe_h5(p), indent=2, default=str))
