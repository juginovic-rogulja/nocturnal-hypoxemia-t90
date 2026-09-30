"""
HSP H5 reader adapted to the AUGUST 2026 BDSP file layout.

Same public interface as the April 2026 feature_pipeline/core/h5_reader.py (HSPH5Reader,
REGION_CHANNELS, ECG_CHANNEL_PRIORITY, probe_h5, has_usable_staging) so the feature
modules (brain_features, microstructure, heart_features, lung_features) run unchanged.

What changed in the files between April and August 2026 and how it is bridged here:

  April layout                                  August layout (this reader)
  root attr sampling_rate = 200                 no root attr; per-channel attr fs (200 for EEG/ECG,
                                                 25 or 10 for spo2, 50 for belts, 1 for cpap_on)
  root attr unit_voltage = 'V'                  per-channel attr unit ('uV', '%', 'bpm', ...)
  /signals/<name> float32 (N,1), volts, 200 Hz  /signals/<name> int16 digital values with attrs
                                                 dig_min, dig_max, phys_min, phys_max (EDF scaling)
  /annotations/stage float32 (N,1) sample-level /annotations/expert_1/stage event table
                                                 (codes int16, starts, ends, durations in seconds)
                                                 + group attr event_map
  /annotations/{arousal,resp,limb} sample-level  /annotations/expert_1/{arousal,limb} event tables,
                                                 resp (I0006) or resp_3 + resp_4 (I0002)

Bridging rules (layout adaptation only, no feature definition is touched):
  * fs of the reader = the sampling rate of the EEG derivations (200 Hz on every file seen).
  * read_raw returns physical values: dig * gain + offset with gain = (phys_max - phys_min) /
    (dig_max - dig_min). A channel sampled below reader fs (SpO2 at 25 or 10 Hz) is brought to
    reader fs by an FFT resample (scipy.signal.resample), because that, and only that, reproduces
    the April oxygen values on the probe nights (T90, T88, ODI3, ODI4, mean): the April files'
    200 Hz SpO2 was an FFT upsample of the oximeter trace, ringing included. resample_mode='hold'
    gives the native step trace instead (used for the *_native columns).
    A channel sampled above reader fs is decimated by scipy.signal.resample_poly (flagged).
  * read_channel_uV converts the physical unit to microvolts (uV as is, mV x 1e3, V x 1e6).
  * _stage_array rebuilds the sample-level stage vector at reader fs from the event table:
    fill 0 (unscored), then paint each event's code over [start, end) in FILE order.
    Events whose start lies beyond the recording end are day-wrap artefacts of the August
    regeneration when start - 86400 falls inside the recording (unwrapped, duration replaced
    by the table's median positive duration when the stored duration is negative), otherwise
    they are dropped. Counts of both are kept for the report.
  * stage_per_epoch / epoch_mask are byte-for-byte the April code (mode per 30 s epoch on the
    recording-start grid).
  * annotations_array(kind) rasterises the matching event table to reader fs (code value on
    event samples, 0 elsewhere). 'resp' falls back to 'resp_3' when only resp_3/resp_4 exist.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Tuple

import h5py
import numpy as np


# Region aggregation for brain features (matches Sun 2024). Unchanged from April.
REGION_CHANNELS: Dict[str, Tuple[str, str]] = {
    "F": ("f3-m2", "f4-m1"),
    "C": ("c3-m2", "c4-m1"),
    "O": ("o1-m2", "o2-m1"),
}

ECG_CHANNEL_PRIORITY: List[str] = ["ecg", "ecg-v1", "ecg-v2", "ecg-ll", "ecg-ra"]

UNIT_TO_UV = {"uv": 1.0, "µv": 1.0, "microvolt": 1.0, "mv": 1e3, "v": 1e6}
DAY_SEC = 86400.0


class HSPH5Reader:
    """Read a single BDSP .h5 file (August 2026 layout). Lazy; opens on context-manager use."""

    def __init__(self, path: str | Path, target_fs: int | None = None, resample_mode: str = "fft"):
        self.path = Path(path)
        self._f: h5py.File | None = None
        self._target_fs = target_fs
        self.resample_mode = resample_mode   # 'fft' reproduces the April files, 'hold' is the native trace
        self._chan: Dict[str, dict] = {}
        self._stage_cache: np.ndarray | None = None
        self.diag: Dict[str, object] = {}

    # ------------------------------------------------------------------ open/close
    def __enter__(self) -> "HSPH5Reader":
        self._f = h5py.File(self.path, "r")
        self._index_channels()
        return self

    def __exit__(self, *exc):
        if self._f is not None:
            self._f.close()
            self._f = None

    def _index_channels(self):
        g = self._f["signals"]
        for k in g.keys():
            d = g[k]
            if not isinstance(d, h5py.Dataset):
                continue
            a = dict(d.attrs)
            info = {"n": int(d.shape[0]), "fs": float(a.get("fs", np.nan)),
                    "unit": str(a.get("unit", "")), "type": str(a.get("type", "")),
                    "dtype": str(d.dtype)}
            if all(x in a for x in ("dig_min", "dig_max", "phys_min", "phys_max")):
                dmin, dmax = float(a["dig_min"]), float(a["dig_max"])
                pmin, pmax = float(a["phys_min"]), float(a["phys_max"])
                gain = (pmax - pmin) / (dmax - dmin) if dmax != dmin else 1.0
                info.update(gain=gain, offset=pmin - dmin * gain, phys_min=pmin, phys_max=pmax,
                            dig_min=dmin, dig_max=dmax, scaled=True)
            else:
                info.update(gain=1.0, offset=0.0, scaled=False)
            self._chan[k.lower()] = info
        # synthesize a missing referential derivation from its two electrodes (the April files
        # carried every derivation; the August files sometimes omit one). Flagged in diag.
        for pair in REGION_CHANNELS.values():
            for der in pair:
                if der not in self._chan:
                    a, b = der.split("-")
                    if a in self._chan and b in self._chan and self._chan[a]["fs"] == self._chan[b]["fs"]:
                        self._chan[der] = {"n": min(self._chan[a]["n"], self._chan[b]["n"]), "fs": self._chan[a]["fs"],
                                           "unit": self._chan[a]["unit"], "type": "synthesized", "dtype": "float32",
                                           "gain": 1.0, "offset": 0.0, "scaled": False, "synth_from": (a, b)}
                        self.diag.setdefault("synthesized_derivations", []).append(der)
        # reader fs: the April files were ALL at 200 Hz (root attr sampling_rate = 200); every
        # channel is brought to 200 Hz here, so the feature code sees what it saw in April.
        if "sampling_rate" in self._f.attrs:
            self._fs = int(self._f.attrs["sampling_rate"])
        else:
            self._fs = int(self._target_fs or 200)
        eeg = [c for pair in REGION_CHANNELS.values() for c in pair if c in self._chan]
        self.diag["eeg_rates"] = sorted({self._chan[c]["fs"] for c in eeg if np.isfinite(self._chan[c]["fs"])})
        self.diag["fs"] = self._fs

    # ------------------------------------------------------------------ attrs
    @property
    def fs(self) -> int:
        return int(self._fs)

    @property
    def n_samples(self) -> int:
        """Length of the reader-fs time base: the EEG derivation length when present,
        else duration_sec * fs."""
        eeg = [c for pair in REGION_CHANNELS.values() for c in pair if c in self._chan
               and np.isfinite(self._chan[c]["fs"])]
        if eeg:
            return max(int(round(self._chan[c]["n"] * self.fs / self._chan[c]["fs"])) for c in eeg)
        d = self._f.attrs.get("duration_sec", None)
        if d is not None and np.isfinite(float(d)):
            return int(round(float(d) * self.fs))
        # fall back to the longest channel rescaled
        return max(int(round(v["n"] * self.fs / v["fs"])) for v in self._chan.values() if np.isfinite(v["fs"]))

    @property
    def duration_sec(self) -> float:
        return self._stage_array().size / self.fs

    @property
    def meas_date(self) -> str:
        return str(self._f.attrs["meas_date"])

    # -------------------------------------------------------------- channels
    def channels(self) -> List[str]:
        return sorted(self._chan.keys())

    def has_channel(self, name: str) -> bool:
        return name.lower() in self._chan

    def channel_info(self, name: str) -> dict:
        return dict(self._chan[name.lower()])

    def _to_reader_fs(self, x: np.ndarray, ch_fs: float, name: str) -> np.ndarray:
        if not np.isfinite(ch_fs) or abs(ch_fs - self.fs) < 1e-6:
            return x
        ratio = self.fs / ch_fs
        if ratio > 1:
            if self.resample_mode == "fft":
                # The April 2026 files carried every channel at 200 Hz. On the four probe nights
                # the April oxygen values (T90, T88, ODI3, ODI4, mean) are reproduced only when the
                # native 25 or 10 Hz trace is brought to 200 Hz with an FFT resample
                # (scipy.signal.resample), not by sample-and-hold or linear interpolation.
                # So the platform's April upsampling was an FFT resample and it is mimicked here.
                from scipy.signal import resample
                n_out = int(round(x.size * ratio))
                self.diag.setdefault("upsampled_fft", []).append(f"{name}:{ch_fs:g}->{self.fs}")
                y = resample(x.astype(np.float64), n_out).astype(np.float32)
            elif abs(ratio - round(ratio)) < 1e-9:
                self.diag.setdefault("upsampled_hold", []).append(f"{name}:{ch_fs:g}->{self.fs}")
                y = np.repeat(x, int(round(ratio)))
            else:
                self.diag.setdefault("upsampled_nearest", []).append(f"{name}:{ch_fs:g}->{self.fs}")
                idx = np.minimum((np.arange(int(round(x.size * ratio))) / ratio).astype(np.int64), x.size - 1)
                y = x[idx]
        else:
            n_out = int(round(x.size * ratio))
            if self.resample_mode == "fft":
                from scipy.signal import resample
                self.diag.setdefault("downsampled_fft", []).append(f"{name}:{ch_fs:g}->{self.fs}")
                y = resample(x.astype(np.float64), n_out).astype(np.float32)
            else:
                from scipy.signal import resample_poly
                from fractions import Fraction
                fr = Fraction(self.fs, int(round(ch_fs))).limit_denominator(1000)
                self.diag.setdefault("downsampled_poly", []).append(f"{name}:{ch_fs:g}->{self.fs}")
                y = resample_poly(x.astype(np.float64), fr.numerator, fr.denominator).astype(np.float32)
        n = self.n_samples
        if y.size > n:
            y = y[:n]
        return y

    def read_raw(self, name: str) -> np.ndarray:
        """Return the signal in its physical unit at reader fs (float32, 1-D)."""
        key = name.lower()
        info = self._chan[key]
        if "synth_from" in info:
            a, b = info["synth_from"]
            xa = self._read_phys(a); xb = self._read_phys(b)
            n = min(xa.size, xb.size)
            x = (xa[:n] - xb[:n]).astype(np.float32)
            return self._to_reader_fs(x, info["fs"], key)
        return self._to_reader_fs(self._read_phys(key), info["fs"], key)

    def _read_phys(self, key: str) -> np.ndarray:
        info = self._chan[key]
        d = self._f["signals"][key]
        raw = d[:].ravel()
        if info["scaled"] and not np.issubdtype(raw.dtype, np.floating):
            x = raw.astype(np.float64) * info["gain"] + info["offset"]
        elif info["scaled"]:
            x = raw.astype(np.float64) * info["gain"] + info["offset"]
        else:
            x = raw.astype(np.float64)
        return x.astype(np.float32)

    def read_channel_uV(self, name: str) -> np.ndarray:
        """Return signal in microvolts (float32, 1-D)."""
        info = self._chan[name.lower()]
        unit = info["unit"].strip().lower()
        factor = UNIT_TO_UV.get(unit, None)
        if factor is None:
            # April files were volts; an unlabelled voltage channel is treated as uV and flagged
            self.diag.setdefault("unit_unknown", []).append(f"{name}:{info['unit']}")
            factor = 1.0
        x = self.read_raw(name)
        if factor != 1.0:
            x = (x * factor).astype(np.float32)
        return x

    # ---------------------------------------------------------- aggregations
    def region_signal_uV(self, region: str) -> np.ndarray:
        left, right = REGION_CHANNELS[region]
        if not (self.has_channel(left) and self.has_channel(right)):
            raise KeyError(f"Missing channel for region {region}: need {left} and {right}")
        return 0.5 * (self.read_channel_uV(left) + self.read_channel_uV(right))

    def ecg_uV(self) -> Tuple[str, np.ndarray]:
        for c in ECG_CHANNEL_PRIORITY:
            if self.has_channel(c):
                return c, self.read_channel_uV(c)
        raise KeyError("No ECG channel present")

    # ------------------------------------------------------------- events
    def _event_table(self, kind: str):
        base = "annotations/expert_1"
        if f"{base}/{kind}" in self._f:
            g = self._f[f"{base}/{kind}"]
        elif kind == "resp" and f"{base}/resp_3" in self._f:
            g = self._f[f"{base}/resp_3"]
            self.diag["resp_source"] = "resp_3"
        else:
            raise KeyError(f"no event table {kind}")
        codes = g["codes"][:].astype(np.int64)
        starts = g["starts"][:].astype(np.float64)
        ends = g["ends"][:].astype(np.float64)
        durs = g["durations"][:].astype(np.float64)
        emap = g.attrs.get("event_map", None)
        try:
            emap = json.loads(emap) if emap is not None else None
        except Exception:
            emap = None
        return codes, starts, ends, durs, emap

    def rasterise(self, kind: str, n: int | None = None) -> np.ndarray:
        """Paint an event table onto a sample-level int16 array at reader fs (0 = no event)."""
        n = self.n_samples if n is None else n
        fs = self.fs
        D = n / fs
        codes, starts, ends, durs, _ = self._event_table(kind)
        out = np.zeros(n, dtype=np.int16)
        pos = durs[durs > 0]
        default_dur = {"stage": 30.0, "arousal": 3.0, "resp": 10.0, "resp_3": 10.0, "resp_4": 10.0, "limb": 0.5}.get(kind, 1.0)
        med_dur = float(np.median(pos)) if pos.size else default_dur
        n_wrapped = n_dropped = n_offgrid = n_negdur = 0
        for c, s, e, d in zip(codes, starts, ends, durs):
            if s > D + 1.0:
                if 0.0 <= s - DAY_SEC <= D:
                    s2 = s - DAY_SEC
                    e2 = e - DAY_SEC if (e - DAY_SEC) > s2 else s2 + (d if d > 0 else med_dur)
                    s, e = s2, e2
                    n_wrapped += 1
                else:
                    n_dropped += 1
                    continue
            elif e <= s:
                n_negdur += 1
                e = s + (d if d > 0 else med_dur)
            if kind == "stage" and abs((s % 30.0) - round(s % 30.0)) < 1e-6 and (round(s) % 30) != 0:
                n_offgrid += 1
            a = int(round(s * fs)); b = int(round(e * fs))
            a = max(a, 0); b = min(b, n)
            if b > a:
                out[a:b] = c
        self.diag[f"{kind}_n_events"] = int(len(codes))
        self.diag[f"{kind}_n_wrapped"] = n_wrapped
        self.diag[f"{kind}_n_dropped_beyond_end"] = n_dropped
        self.diag[f"{kind}_n_negdur_fixed"] = n_negdur
        if kind == "stage":
            self.diag["stage_n_offgrid"] = n_offgrid
            self.diag["stage_med_dur"] = med_dur
        return out

    def event_counts(self, kind: str) -> Dict[str, int]:
        """Row count, rising-edge count after rasterising (touching events merge), and
        rising-edge count of the event starts only (no merging), for the report."""
        codes, starts, ends, durs, _ = self._event_table(kind)
        ras = self.rasterise(kind)
        rising = int(np.sum(np.diff((ras != 0).astype(np.int8)) == 1)) + int(ras[0] != 0)
        rising_april_style = int(np.sum(np.diff((ras != 0).astype(np.int8)) == 1))
        return {"rows": int(len(codes)), "rows_nonzero": int(np.sum(codes != 0)),
                "rising_edges": rising_april_style, "rising_edges_incl_first": rising}

    def has_table(self, kind: str) -> bool:
        return f"annotations/expert_1/{kind}" in self._f

    # ------------------------------------------------------------- stages
    def stage_event_map(self) -> dict | None:
        try:
            return self._event_table("stage")[4]
        except KeyError:
            return None

    def _stage_array(self) -> np.ndarray:
        if self._stage_cache is None:
            self._stage_cache = self.rasterise("stage").astype(np.float32)
        return self._stage_cache

    def stage_per_epoch(self, epoch_sec: int = 30) -> np.ndarray:
        """April code, unchanged: per-epoch mode of the sample-level stage on the
        recording-start grid."""
        stage = self._stage_array()
        samples_per_epoch = epoch_sec * self.fs
        n_epochs = stage.size // samples_per_epoch
        stage = stage[: n_epochs * samples_per_epoch].reshape(n_epochs, samples_per_epoch)
        out = np.zeros(n_epochs, dtype=np.int16)
        for i in range(n_epochs):
            vals, cts = np.unique(stage[i][~np.isnan(stage[i])], return_counts=True)
            out[i] = int(vals[cts.argmax()]) if vals.size else -1
        return out

    def epoch_mask(self, stage_codes: List[int], epoch_sec: int = 30) -> np.ndarray:
        st = self.stage_per_epoch(epoch_sec)
        return np.isin(st, stage_codes)

    def annotations_array(self, kind: str) -> np.ndarray:
        """Sample-level annotation array for 'arousal', 'resp', or 'limb' (float32)."""
        return self.rasterise(kind).astype(np.float32)


def probe_h5(path: str | Path) -> Dict:
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


def has_usable_staging_from_epochs(st: np.ndarray, min_scored_epochs: int = 200) -> bool:
    """April rule: False if more than 95 percent of epochs are 0 or 9, or fewer than 200 scored."""
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


def has_usable_staging(path: str | Path, min_scored_epochs: int = 200) -> bool:
    with HSPH5Reader(path) as r:
        st = r.stage_per_epoch()
    return has_usable_staging_from_epochs(st, min_scored_epochs)
