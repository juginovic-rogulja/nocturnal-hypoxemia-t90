#!/usr/bin/env python3
"""
Per-STAGE T90 on each side of a within-night split, and for whole nights.

Merges the two verified extractors:
  * pilots/cpap_t90/extract.py     -> PAP channel handling + sustained switch-on detection
  * pilots/t90_by_stage/extract.py -> staging -> sample-level stage array, per-state T90

For every session it returns, for each of three segmentations:
    all   : whole recording
    pre / post  : split at the first SUSTAINED PAP switch-on (Design A)
    h1  / h2    : split at `split_min` from the worklist (placebo control; default 179 min)
per state in {wake, sleep, nrem, n1, n2, n3, rem}:
    t90, valid-sample count, minutes, mean SpO2, nadir.

Definitions kept byte-identical to the existing work so numbers stay comparable:
    valid = spo2[(spo2 >= 50) & (spo2 <= 100)];  T90 = mean(valid < 90.0) * 100

Traps handled (all previously verified in this project)
  * channels.tsv mislabels EVERY channel type=EEG units=uV -> scaling comes ONLY from the
    DATASET attrs dig_min/dig_max/phys_min/phys_max.
      I0002 int16 over +-32768 mapped 0-100 %, fs 25 Hz | I0006 identity, 10 Hz | I0004 4 Hz
  * BDSP stage codes are 1=N3, 2=N2, 3=N1, 4=REM, 5=Wake.  NOT depth ordered.
  * I0002 stage `starts` are not sorted and some `durations` are corrupt negatives ->
    epoch length from ends-starts with a 30 s fallback when outside (0, 60].
  * SpO2 aligned to staging by TIME (index = round(sec*fs)), never by position.
  * cpap_on is a 1 Hz derived boolean; cpres is cmH2O.  Both aligned by time.
"""
from __future__ import annotations
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../..")))  # repository root, where paths.py lives
import paths
import argparse, time
from pathlib import Path

import numpy as np
import pandas as pd
import s3fs, h5py

AP = f"{paths.BDSP_S3_ACCESS_POINT}/PSG/bids"
HERE = Path(__file__).parent

N3, N2, N1, REM, WAKE = 1, 2, 3, 4, 5          # BDSP codes, NOT depth ordered
SLEEP_CODES = (N1, N2, N3, REM)
NREM_CODES = (N1, N2, N3)
STATES = {"wake": (WAKE,), "sleep": SLEEP_CODES, "nrem": NREM_CODES,
          "n1": (N1,), "n2": (N2,), "n3": (N3,), "rem": (REM,)}
BLOCK = 32 * 1024        # measured 0.51 MB/file, no wall-time penalty vs 256 KB here
DEFAULT_SPLIT = 179.0    # median PAP switch-on minute in the split-night pilot


def phys(dset, raw):
    a = dset.attrs
    try:
        dmin, dmax = float(a["dig_min"]), float(a["dig_max"])
        pmin, pmax = float(a["phys_min"]), float(a["phys_max"])
    except KeyError:
        return raw.astype(np.float32)
    if dmax == dmin:
        return raw.astype(np.float32)
    return (pmin + (raw.astype(np.float64) - dmin) * (pmax - pmin)
            / (dmax - dmin)).astype(np.float32)


def t90(x):
    v = x[(x >= 50.0) & (x <= 100.0)]
    if v.size < 100:
        return dict(t90=np.nan, n=int(v.size), mean=np.nan, nadir=np.nan)
    return dict(t90=float(np.mean(v < 90.0) * 100), n=int(v.size),
                mean=float(np.mean(v)), nadir=float(np.min(v)))


def _b(r):
    try:
        return int(r.cache.total_requested_bytes)
    except Exception:
        return np.nan


def seg_states(out, spo2, stage, segmask, fs, tag):
    """Per-state T90 within one time segment."""
    for nm, codes in STATES.items():
        m = np.isin(stage, codes) & segmask
        res = t90(spo2[m])
        out[f"{tag}_{nm}_t90"] = res["t90"]
        out[f"{tag}_{nm}_n"] = res["n"]
        out[f"{tag}_{nm}_min"] = float(m.sum() / fs / 60.0)
        out[f"{tag}_{nm}_mean"] = res["mean"]
        out[f"{tag}_{nm}_nadir"] = res["nadir"]
    r = t90(spo2[segmask])                       # unstaged, whole segment
    out[f"{tag}_any_t90"] = r["t90"]
    out[f"{tag}_any_min"] = float(segmask.sum() / fs / 60.0)
    out[f"{tag}_any_nadir"] = r["nadir"]


def one(fsys, site, sub, ses, split_min):
    key = f"{AP}/{site}/{sub}/ses-{ses}/eeg/{sub}_ses-{ses}.h5"
    out = {"status": "ok"}
    with fsys.open(key, "rb") as r:
        with h5py.File(r, "r") as h5:
            sig = h5["signals"]
            names = {k.lower(): k for k in sig.keys()}
            if "spo2" not in names:
                out["status"] = "no_spo2"; out["bytes"] = _b(r); return out
            d = sig[names["spo2"]]
            fs = float(d.attrs.get("fs", 25.0))
            spo2 = phys(d, np.asarray(d[:]).ravel())
            n = spo2.size
            out["fs_spo2"] = fs
            out["rec_min"] = n / fs / 60.0
            out["spo2_valid_frac"] = float(((spo2 >= 50) & (spo2 <= 100)).mean()) if n else np.nan
            out.update({f"all_any_{k}": v for k, v in t90(spo2).items()})
            if n < 100 * fs:
                out["status"] = "too_short"; out["bytes"] = _b(r); return out
            t = np.arange(n) / fs                       # seconds

            # ---------------- PAP -------------------------------------------
            papw = np.zeros(n, bool)
            out["has_cpap_chan"] = "cpap_on" in names
            out["has_cpres"] = "cpres" in names
            cf, src = np.nan, "none"
            if "cpap_on" in names:
                dc = sig[names["cpap_on"]]
                fsc = float(dc.attrs.get("fs", 1.0))
                con = np.asarray(dc[:]).ravel().astype(np.float32) > 0.5
                src = "cpap_on"
            elif "cpres" in names:
                dp = sig[names["cpres"]]
                fsc = float(dp.attrs.get("fs", 25.0))
                p = phys(dp, np.asarray(dp[:]).ravel())
                raw = p > 3.0                            # cmH2O, ambient near 0
                w = max(int(round(30 * fsc)), 1)
                # v8.2 (2026-09-15): integer count, the same rule (more than half the window above 3 cmH2O) without a floating-point tie
                con = 2 * np.convolve(raw.astype(np.int64), np.ones(w, np.int64), "same") > w
                src = "cpres"
            else:
                con = np.zeros(0, bool)
            out["pap_src"] = src
            switch = np.nan
            if con.size:
                cf = float(con.mean())
                ii = np.clip(np.floor(t * fsc).astype(np.int64), 0, con.size - 1)
                papw = con[ii]
                # first SUSTAINED switch-on: >=90 % on over the following 30 min
                if con.any():
                    w = int(round(30 * 60 * fsc))
                    if con.size > w:
                        # v8.2 (2026-09-15): integer count of on-samples per window; the rule (at least 90 percent on) is unchanged.
                        # The float form (mean >= 0.90) sat on an exact tie at 1620 of 1800 samples on nearly every sustained
                        # switch-on; the tie evaluates >= 0.90 on this Mac (the 2026-07-29 extraction) and < 0.90 on the EC2
                        # build (numpy accumulation order), a one-sample (1 s) shift of the switch-on on 95 percent of split nights.
                        run = np.convolve(con.astype(np.int64), np.ones(w, np.int64), "valid")
                        g = np.flatnonzero(run >= int(np.ceil(0.90 * w - 1e-9)))
                        if g.size:
                            switch = float(g[0] / fsc / 60.0)
                    if not np.isfinite(switch):
                        switch = float(np.flatnonzero(con)[0] / fsc / 60.0)
            out["cpap_frac"] = cf
            out["cpap_start_min"] = switch

            # ---------------- staging ---------------------------------------
            try:
                g = h5["annotations/expert_1"]["stage"]
                ss = np.asarray(g["starts"][:], np.float64)
                cc = np.asarray(g["codes"][:], np.int32)
                try:
                    dur = np.asarray(g["ends"][:], np.float64) - ss
                except KeyError:
                    dur = np.asarray(g["durations"][:], np.float64)
            except Exception as e:
                out["status"] = f"no_stage:{type(e).__name__}"; out["bytes"] = _b(r); return out
            dur = np.where(~np.isfinite(dur) | (dur <= 0) | (dur > 60), 30.0, dur)
            o = np.argsort(ss, kind="stable")
            ss, dur, cc = ss[o], dur[o], cc[o]
            out["n_stage_ep"] = int(cc.size)
            stage = np.zeros(n, np.int8)
            i0 = np.rint(ss * fs).astype(np.int64)
            i1 = np.rint((ss + dur) * fs).astype(np.int64)
            for a, b, c in zip(i0, i1, cc):
                if b <= 0 or a >= n or c < 1 or c > 5:
                    continue
                stage[max(a, 0):min(b, n)] = c
            sleepmask = np.isin(stage, SLEEP_CODES)
            out["tst_min"] = float(sleepmask.sum() / fs / 60.0)
            out["staged_min"] = float((stage > 0).sum() / fs / 60.0)
            if not sleepmask.any():
                out["status"] = "no_sleep_epochs"; out["bytes"] = _b(r); return out

            allm = np.ones(n, bool)
            seg_states(out, spo2, stage, allm, fs, "all")

            # ---- Design A: split at the sustained PAP switch-on --------------
            if np.isfinite(switch) and switch > 0:
                pre = t < switch * 60.0
                seg_states(out, spo2, stage, pre, fs, "pre")
                seg_states(out, spo2, stage, ~pre, fs, "post")
                out["pap_frac_pre"] = float(papw[pre].mean()) if pre.any() else np.nan
                out["pap_frac_post"] = float(papw[~pre].mean()) if (~pre).any() else np.nan

            # ---- placebo / fixed split -------------------------------------
            sp = float(split_min) if np.isfinite(split_min) else DEFAULT_SPLIT
            out["split_min_used"] = sp
            h1 = t < sp * 60.0
            if h1.any() and (~h1).any():
                seg_states(out, spo2, stage, h1, fs, "h1")
                seg_states(out, spo2, stage, ~h1, fs, "h2")
                out["pap_frac_h1"] = float(papw[h1].mean())
                out["pap_frac_h2"] = float(papw[~h1].mean())
        out["bytes"] = _b(r)
    return out


KEEP = ["BDSPPatientID", "SiteID", "BIDSFolder", "SessionID", "AgeAtVisit", "SexDSC",
        "arm", "rank", "pair_gap", "sdt", "split_min", "AHI", "TST_min",
        "spo2_pct_below_90", "cpap_frac_prior"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--infile", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--maxsec", type=float, default=0.0)
    args = ap.parse_args()

    todo = pd.read_csv(args.infile, low_memory=False)
    if args.nshards > 1:
        todo = todo.iloc[args.shard::args.nshards]
    if args.limit:
        todo = todo.head(args.limit)
    outp = Path(f"{args.out}_{args.shard}.csv")
    rows, done = [], set()
    if outp.exists():
        prev = pd.read_csv(outp, low_memory=False)
        rows = prev.to_dict("records")
        done = set(zip(prev.BIDSFolder, prev.SessionID))

    fsys = s3fs.S3FileSystem(anon=False, default_block_size=BLOCK)
    t0, nb = time.time(), 0.0
    for _, r in todo.iterrows():
        if (r.BIDSFolder, int(r.SessionID)) in done:
            continue
        rec = {c: r[c] for c in KEEP if c in todo.columns}
        try:
            rec.update(one(fsys, r.SiteID, r.BIDSFolder, int(r.SessionID),
                           float(r.get("split_min", np.nan)) if "split_min" in todo.columns
                           else np.nan))
        except Exception as e:
            rec.update({"status": f"ERR:{type(e).__name__}", "err": str(e)[:150]})
        b = rec.get("bytes", np.nan)
        nb += b if isinstance(b, (int, float, np.integer)) and np.isfinite(b) else 0
        rows.append(rec)
        if args.maxsec and (time.time() - t0) > args.maxsec:
            pd.DataFrame(rows).to_csv(outp, index=False)
            print(f"[s{args.shard}] TIMEBOX n={len(rows)} GB={nb/1e9:.3f}", flush=True)
            return
        if len(rows) % 50 == 0:
            pd.DataFrame(rows).to_csv(outp, index=False)
            print(f"[s{args.shard}] {len(rows)}/{len(todo)} {time.time()-t0:.0f}s "
                  f"GB={nb/1e9:.3f}", flush=True)
    pd.DataFrame(rows).to_csv(outp, index=False)
    print(f"[s{args.shard}] DONE n={len(rows)} {time.time()-t0:.0f}s GB={nb/1e9:.3f}", flush=True)


if __name__ == "__main__":
    main()
