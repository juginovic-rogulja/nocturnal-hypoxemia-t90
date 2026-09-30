#!/usr/bin/env python3
"""
Split WAKE into its three physiological windows and compute T90 in each.

The existing extract.py pools every wake epoch, so its "wake T90" mixes the resting evening
baseline with post-apnoeic arousals during the night. Those are different things: the evening
window is the only one that measures a patient's oxygenation before sleep has done anything to
them, which is the quantity that separates a gas-exchange problem from an airway problem.

  evening  wake epochs BEFORE the first sleep sample
  waso     wake epochs BETWEEN the first and last sleep sample
  morning  wake epochs AFTER the last sleep sample

Each is computed over the whole recording and, separately, within the untreated portion before
the sustained PAP switch-on. All conventions are copied byte-identically from extract.py so the
numbers stay comparable: valid = spo2 in [50,100], T90 = mean(valid < 90) * 100, BDSP stage
codes 1=N3 2=N2 3=N1 4=REM 5=Wake, scaling from DATASET attrs only, SpO2 aligned to staging by
time and never by position.
"""
from __future__ import annotations
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../..")))  # repository root, where paths.py lives
import paths
import argparse, time
from pathlib import Path
import numpy as np, pandas as pd, s3fs, h5py

AP = f"{paths.BDSP_S3_ACCESS_POINT}/PSG/bids"
N3, N2, N1, REM, WAKE = 1, 2, 3, 4, 5
SLEEP_CODES = (N1, N2, N3, REM)
BLOCK = 32 * 1024


def phys(dset, raw):
    a = dset.attrs
    try:
        dmin, dmax = float(a["dig_min"]), float(a["dig_max"])
        pmin, pmax = float(a["phys_min"]), float(a["phys_max"])
        if dmax > dmin and pmax > pmin:
            return pmin + (raw.astype(np.float64) - dmin) * (pmax - pmin) / (dmax - dmin)
    except Exception:
        pass
    return raw.astype(np.float64)


def t90_of(spo2, mask):
    """T90, mean, nadir and minutes over a boolean mask, using only valid samples."""
    v = spo2[mask]
    v = v[(v >= 50) & (v <= 100)]
    if v.size == 0:
        return {"t90": np.nan, "mean": np.nan, "nadir": np.nan, "n": 0}
    return {"t90": float((v < 90.0).mean() * 100), "mean": float(v.mean()),
            "nadir": float(v.min()), "n": int(v.size)}


def one(fsys, site, sub, ses):
    key = f"{AP}/{site}/{sub}/ses-{ses}/eeg/{sub}_ses-{ses}.h5"
    out = {"status": "ok"}
    with fsys.open(key, "rb") as r:
        with h5py.File(r, "r") as h5:
            sig = h5["signals"]
            names = {k.lower(): k for k in sig.keys()}
            if "spo2" not in names:
                out["status"] = "no_spo2"; return out
            d = sig[names["spo2"]]
            fs = float(d.attrs.get("fs", 25.0))
            spo2 = phys(d, np.asarray(d[:]).ravel())
            n = spo2.size
            if n < 100 * fs:
                out["status"] = "too_short"; return out
            t = np.arange(n) / fs

            # PAP switch-on, same sustained-90%-over-30-min rule as extract.py
            switch = np.nan
            con, fsc = np.zeros(0, bool), 1.0
            if "cpap_on" in names:
                dc = sig[names["cpap_on"]]; fsc = float(dc.attrs.get("fs", 1.0))
                con = np.asarray(dc[:]).ravel().astype(np.float32) > 0.5
            elif "cpres" in names:
                dp = sig[names["cpres"]]; fsc = float(dp.attrs.get("fs", 25.0))
                p = phys(dp, np.asarray(dp[:]).ravel())
                w = max(int(round(30 * fsc)), 1)
                con = np.convolve((p > 3.0).astype(np.float32), np.ones(w) / w, "same") > 0.5
            if con.size and con.any():
                w = int(round(30 * 60 * fsc))
                if con.size > w:
                    run = np.convolve(con.astype(np.float32), np.ones(w) / w, "valid")
                    g = np.flatnonzero(run >= 0.90)
                    if g.size:
                        switch = float(g[0] / fsc / 60.0)
                if not np.isfinite(switch):
                    switch = float(np.flatnonzero(con)[0] / fsc / 60.0)
            out["cpap_start_min"] = switch

            # staging -> sample-level stage array
            try:
                g = h5["annotations/expert_1"]["stage"]
                ss = np.asarray(g["starts"][:], np.float64)
                cc = np.asarray(g["codes"][:], np.int32)
                try:
                    dur = np.asarray(g["ends"][:], np.float64) - ss
                except KeyError:
                    dur = np.asarray(g["durations"][:], np.float64)
            except Exception as e:
                out["status"] = f"no_stage:{type(e).__name__}"; return out
            dur = np.where(~np.isfinite(dur) | (dur <= 0) | (dur > 60), 30.0, dur)
            o = np.argsort(ss, kind="stable")
            ss, dur, cc = ss[o], dur[o], cc[o]
            stage = np.zeros(n, np.int8)
            for a, b, c in zip(np.rint(ss * fs).astype(np.int64),
                               np.rint((ss + dur) * fs).astype(np.int64), cc):
                if b <= 0 or a >= n or c < 1 or c > 5:
                    continue
                stage[max(a, 0):min(b, n)] = c

            sleepmask = np.isin(stage, SLEEP_CODES)
            if not sleepmask.any():
                out["status"] = "no_sleep_epochs"; return out
            si = np.flatnonzero(sleepmask)
            first, last = si[0], si[-1]
            idx = np.arange(n)
            wake = stage == WAKE
            WIN = {"evening": wake & (idx < first),
                   "waso":    wake & (idx >= first) & (idx <= last),
                   "morning": wake & (idx > last),
                   "sleep":   sleepmask}
            out["sleep_onset_min"] = float(first / fs / 60.0)
            out["final_wake_min"] = float(last / fs / 60.0)

            for w, m in WIN.items():
                for k, v in t90_of(spo2, m).items():
                    out[f"all_{w}_{k}"] = v
                out[f"all_{w}_min"] = float(m.sum() / fs / 60.0)
            if np.isfinite(switch) and switch > 0:
                pre = t < switch * 60.0
                for w, m in WIN.items():
                    for k, v in t90_of(spo2, m & pre).items():
                        out[f"pre_{w}_{k}"] = v
                    out[f"pre_{w}_min"] = float((m & pre).sum() / fs / 60.0)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--infile", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    todo = pd.read_csv(a.infile, low_memory=False)
    if a.limit:
        todo = todo.head(a.limit)
    outp = Path(a.out)
    rows, done = [], set()
    if outp.exists():
        prev = pd.read_csv(outp, low_memory=False)
        rows = prev.to_dict("records")
        done = set(zip(prev.BIDSFolder, prev.SessionID))
    fsys = s3fs.S3FileSystem(anon=False, default_block_size=BLOCK)
    t0 = time.time()
    for i, (_, r) in enumerate(todo.iterrows()):
        if (r.BIDSFolder, int(r.SessionID)) in done:
            continue
        rec = {c: r[c] for c in todo.columns}
        try:
            rec.update(one(fsys, r.SiteID, r.BIDSFolder, int(r.SessionID)))
        except Exception as e:
            rec.update({"status": f"ERR:{type(e).__name__}", "err": str(e)[:150]})
        rows.append(rec)
        if (i + 1) % 25 == 0:
            pd.DataFrame(rows).to_csv(outp, index=False)
            el = time.time() - t0
            print(f"  {i+1}/{len(todo)}  {el:.0f}s  {el/(i+1):.2f}s/rec  "
                  f"eta {(len(todo)-i-1)*el/(i+1)/60:.0f} min", flush=True)
    pd.DataFrame(rows).to_csv(outp, index=False)
    print(f"done {len(rows)} rows in {(time.time()-t0)/60:.1f} min -> {outp}")


if __name__ == "__main__":
    main()
