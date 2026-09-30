"""
PrePapReader: the v7 HSPH5Reader with ONE addition, a time cut at t_end_sec (the first sustained PAP switch-on),
implemented in the three places every feature computation goes through (plan B4):
  * n_samples          -> min(full length, round(t_end * fs)); the stage vector, every rasterised event table and every
                          resampled signal are therefore cut on the recording-start grid (samples at or after t_end are
                          absent, which the architecture code counts as unscored/nothing).
  * read_raw           -> the reader-fs signal truncated at n_samples (same-rate channels were not truncated by the base
                          class; resampled ones already were). The FFT resample is applied to the WHOLE native trace first,
                          exactly as in the whole-night pass, then cut, so the April-style samples are the same samples.
  * _event_table       -> rows whose (day-wrap-corrected) start lies at or after t_end are dropped, ends are clipped at t_end,
                          so row counts (AHI_rows, limb_rows) and rising-edge counts see the window only.
  plus read_native(name): the physical trace at the channel's own rate, cut at round(t_end * fs_native) (the v7
  lung_variants read this dataset directly, bypassing read_raw).
With t_end_sec=None every override is a pass-through (step 012 proves it on the 12 probe nights to 1e-6).
Light mode (light=True) indexes only spo2, the EEG derivations (for the time base) and the ECG, for targeted s3fs reads.
"""
from __future__ import annotations
import os, sys
from pathlib import Path
import numpy as np
import h5py

HERE = os.path.dirname(os.path.abspath(__file__))
REEX = os.path.join(os.path.dirname(HERE), "reex_code", "scripts")
if REEX not in sys.path:
    sys.path.insert(0, REEX)
from reex_pipeline.core.h5_reader import HSPH5Reader, REGION_CHANNELS, ECG_CHANNEL_PRIORITY, DAY_SEC  # noqa: E402

LIGHT_CHANNELS = ["spo2"] + [c for pair in REGION_CHANNELS.values() for c in pair] + ECG_CHANNEL_PRIORITY
DEFAULT_DUR = {"stage": 30.0, "arousal": 3.0, "resp": 10.0, "resp_3": 10.0, "resp_4": 10.0, "limb": 0.5}


def clean_event_times(codes, starts, ends, durs, D, kind):
    """The per-row rules of HSPH5Reader.rasterise, returned as times (s) instead of a painted array:
    day-wrapped starts are unwrapped when start - 86400 lies inside [0, D], starts beyond that are dropped, an end at or
    before its start is replaced by start + duration (or the table's median positive duration, else the kind's default).
    Ends are NOT clipped at D here (rasterise clips by the array length); callers clip as they need."""
    codes = np.asarray(codes); starts = np.asarray(starts, float); ends = np.asarray(ends, float); durs = np.asarray(durs, float)
    pos = durs[durs > 0]
    med_dur = float(np.median(pos)) if pos.size else DEFAULT_DUR.get(kind, 1.0)
    oc, os_, oe = [], [], []
    for c, s, e, d in zip(codes, starts, ends, durs):
        if s > D + 1.0:
            if 0.0 <= s - DAY_SEC <= D:
                s2 = s - DAY_SEC
                e2 = e - DAY_SEC if (e - DAY_SEC) > s2 else s2 + (d if d > 0 else med_dur)
                s, e = s2, e2
            else:
                continue
        elif e <= s:
            e = s + (d if d > 0 else med_dur)
        oc.append(c); os_.append(s); oe.append(e)
    return np.asarray(oc), np.asarray(os_, float), np.asarray(oe, float)


class PrePapReader(HSPH5Reader):
    def __init__(self, path, t_end_sec=None, fileobj=None, light=False, **kw):
        super().__init__(path if path is not None else "s3fs-object", **kw)
        t = None if t_end_sec is None else float(t_end_sec)
        self.t_end_sec = t if (t is not None and np.isfinite(t)) else None
        self._fileobj = fileobj
        self._light = bool(light)
        self._n_cut = None
        self._n_full = None
        self.window_end_sec = None

    # ---------------------------------------------------------------- open
    def __enter__(self):
        self._f = h5py.File(self._fileobj if self._fileobj is not None else self.path, "r")
        self._index_channels()
        self._n_full = HSPH5Reader.n_samples.fget(self)
        if self.t_end_sec is not None:
            n_cut = int(round(self.t_end_sec * self.fs))
            self._n_cut = int(min(self._n_full, max(n_cut, 0)))
            self.window_end_sec = self._n_cut / self.fs
            self.diag["prepap_t_end_sec"] = self.t_end_sec
            self.diag["prepap_n_full"] = int(self._n_full)
            self.diag["prepap_n_cut"] = int(self._n_cut)
            self.diag["prepap_window_end_sec"] = self.window_end_sec
        return self

    def _index_channels(self):
        if not self._light:
            return super()._index_channels()
        # light mode: attrs of the needed channels only (each attrs read is one range request over s3fs)
        g = self._f["signals"]
        keys = {k.lower(): k for k in g.keys()}
        for lk in LIGHT_CHANNELS:
            if lk not in keys:
                continue
            d = g[keys[lk]]
            if not isinstance(d, h5py.Dataset):
                continue
            a = dict(d.attrs)
            info = {"n": int(d.shape[0]), "fs": float(a.get("fs", np.nan)), "unit": str(a.get("unit", "")),
                    "type": str(a.get("type", "")), "dtype": str(d.dtype)}
            if all(x in a for x in ("dig_min", "dig_max", "phys_min", "phys_max")):
                dmin, dmax = float(a["dig_min"]), float(a["dig_max"]); pmin, pmax = float(a["phys_min"]), float(a["phys_max"])
                gain = (pmax - pmin) / (dmax - dmin) if dmax != dmin else 1.0
                info.update(gain=gain, offset=pmin - dmin * gain, phys_min=pmin, phys_max=pmax, dig_min=dmin, dig_max=dmax, scaled=True)
            else:
                info.update(gain=1.0, offset=0.0, scaled=False)
            self._chan[lk] = info
        self._fs = int(self._f.attrs["sampling_rate"]) if "sampling_rate" in self._f.attrs else int(self._target_fs or 200)
        eeg = [c for pair in REGION_CHANNELS.values() for c in pair if c in self._chan]
        self.diag["eeg_rates"] = sorted({self._chan[c]["fs"] for c in eeg if np.isfinite(self._chan[c]["fs"])})
        self.diag["fs"] = self._fs
        self.diag["light_index"] = True

    # ---------------------------------------------------------------- cut points
    @property
    def n_samples(self) -> int:
        n = HSPH5Reader.n_samples.fget(self)
        return n if self._n_cut is None else int(min(n, self._n_cut))

    @property
    def n_samples_full(self) -> int:
        return int(self._n_full)

    def read_raw(self, name: str) -> np.ndarray:
        x = super().read_raw(name)
        return x if self._n_cut is None else x[: self._n_cut]

    def _event_table(self, kind: str):
        codes, starts, ends, durs, emap = super()._event_table(kind)
        if self._n_cut is None:
            return codes, starts, ends, durs, emap
        D_full = self._n_full / self.fs
        t = self.window_end_sec
        shift = np.where((starts > D_full + 1.0) & (starts - DAY_SEC >= 0.0) & (starts - DAY_SEC <= D_full), DAY_SEC, 0.0)
        keep = (starts - shift) < t
        ends2 = np.minimum(ends, t + shift)
        self.diag[f"prepap_{kind}_rows_in"] = int(len(codes))
        self.diag[f"prepap_{kind}_rows_kept"] = int(keep.sum())
        return codes[keep], starts[keep], ends2[keep], durs[keep], emap

    def read_native(self, name: str) -> np.ndarray:
        """Physical values at the channel's own rate (float32), cut at t_end. Same arithmetic as v7 lung_variants."""
        key = name.lower()
        info = self._chan[key]
        raw = self._f["signals"][key][:].ravel().astype(np.float64) * info["gain"] + info["offset"]
        if self._n_cut is not None and np.isfinite(info["fs"]) and info["fs"] > 0:
            raw = raw[: int(round(self.window_end_sec * info["fs"]))]
        return raw.astype(np.float32)

    def clean_events(self, kind: str):
        """(codes, start_s, end_s) of the kind's events after the reader's own wrap/negative-duration rules, in the
        window (ends clipped at the window end)."""
        codes, starts, ends, durs, _ = self._event_table(kind)
        D = self.n_samples / self.fs
        c, s, e = clean_event_times(codes, starts, ends, durs, D, kind)
        if s.size:
            e = np.minimum(e, D)
        return c, s, e

    def event_map(self, kind: str):
        try:
            return self._event_table(kind)[4]
        except KeyError:
            return None
