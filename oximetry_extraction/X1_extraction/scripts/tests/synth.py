"""Synthetic BDSP-layout h5 for the unit tests (August 2026 layout: int16 signals with EDF scaling attrs, event tables
under annotations/expert_1 with codes/starts/ends/durations and an event_map attr)."""
from __future__ import annotations
import json
import numpy as np
import h5py

LEGEND = {"0": "Unknown/Unscored", "1": "N3", "2": "N2", "3": "N1", "4": "REM", "5": "Wake"}
RESP_MAP = {"1": "Obstructive Apnea", "2": "Central Apnea", "9": "Hypopnea (unspecified)", "0": "No Event"}
LIMB_MAP = {"1": "Limb Movement", "3": "Left Leg", "4": "Right Leg", "0": "No Event"}


def _sig(g, name, x, fs, unit="uV", pmin=-3200.0, pmax=3200.0):
    x = np.asarray(x, float)
    dmin, dmax = -32768.0, 32767.0
    gain = (pmax - pmin) / (dmax - dmin)
    dig = np.clip(np.round((x - pmin) / gain + dmin), dmin, dmax).astype(np.int16)
    d = g.create_dataset(name, data=dig.reshape(-1, 1))
    d.attrs["fs"] = float(fs); d.attrs["unit"] = unit; d.attrs["type"] = "raw"
    d.attrs["dig_min"] = dmin; d.attrs["dig_max"] = dmax; d.attrs["phys_min"] = pmin; d.attrs["phys_max"] = pmax


def _table(g, name, codes, starts, ends, emap):
    t = g.create_group(name)
    codes = np.asarray(codes, np.int16); starts = np.asarray(starts, float); ends = np.asarray(ends, float)
    t.create_dataset("codes", data=codes); t.create_dataset("starts", data=starts); t.create_dataset("ends", data=ends)
    t.create_dataset("durations", data=ends - starts)
    t.attrs["event_map"] = json.dumps(emap)


def make_synthetic_h5(path, dur_sec=3600, fs=200, spo2_fs=25, seed=1, stage_epochs=None, resp=None, limb=None,
                      spo2=None, wrap_arousal=True):
    """dur_sec seconds; stage_epochs = list of codes per 30 s epoch (default: 10 wake, then N1 N2 N3 REM cycle);
    resp = list of (start, end, code); limb = list of (start, end, code); spo2 = native trace or None (synthetic)."""
    rng = np.random.default_rng(seed)
    n = dur_sec * fs
    t = np.arange(n) / fs
    with h5py.File(path, "w") as f:
        f.attrs["duration_sec"] = np.float32(dur_sec); f.attrs["meas_date"] = "2020-01-01 22:00:00"
        f.attrs["hypopnea_rule_source"] = "both"
        g = f.create_group("signals")
        for ch in ("f3-m2", "f4-m1", "c3-m2", "c4-m1", "o1-m2", "o2-m1"):
            x = 30 * np.sin(2 * np.pi * 1.0 * t + rng.uniform(0, 6)) + 10 * np.sin(2 * np.pi * 13 * t) + rng.normal(0, 5, n)
            _sig(g, ch, x, fs)
        ecg = np.zeros(n)
        beats = np.arange(0.5, dur_sec, 0.9)
        for b in beats:
            i = int(b * fs)
            if i + 4 < n:
                ecg[i:i + 4] = [200, 800, -300, 50]
        _sig(g, "ecg", ecg + rng.normal(0, 5, n), fs, pmin=-5000, pmax=5000)
        if spo2 is None:
            spo2 = np.full(dur_sec * spo2_fs, 96.0)
        _sig(g, "spo2", spo2, spo2_fs, unit="%", pmin=0.0, pmax=100.0)
        a = f.create_group("annotations/expert_1")
        n_ep = dur_sec // 30
        if stage_epochs is None:
            cyc = [5] * 10 + ([3, 3, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 4, 4, 4, 5] * 20)
            stage_epochs = cyc[:n_ep]
        stage_epochs = list(stage_epochs)[:n_ep]
        st = np.arange(n_ep) * 30.0
        _table(a, "stage", stage_epochs, st, st + 30.0, LEGEND)
        ar_s = np.array([600.0, 1200.0, 1800.0]); ar_s = ar_s[ar_s < dur_sec]
        if wrap_arousal:
            _table(a, "arousal", [1] * len(ar_s), ar_s + 86400.0, ar_s + 3.0, {"0": "No Arousal", "1": "Arousal"})
        else:
            _table(a, "arousal", [1] * len(ar_s), ar_s, ar_s + 3.0, {"0": "No Arousal", "1": "Arousal"})
        resp = resp or []
        _table(a, "resp_3", [c for _, _, c in resp], [s for s, _, _ in resp], [e for _, e, _ in resp], RESP_MAP)
        _table(a, "resp_4", [c for _, _, c in resp[::2]], [s for s, _, _ in resp[::2]], [e for _, e, _ in resp[::2]], RESP_MAP)
        limb = limb or []
        _table(a, "limb", [c for _, _, c in limb], [s for s, _, _ in limb], [e for _, e, _ in limb], LIMB_MAP)
    return path
