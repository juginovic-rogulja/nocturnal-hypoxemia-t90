#!/usr/bin/env python3
"""
Overnight oxygen profile extraction worker. Runs on EC2 (us-east-1) under the
EC2BDSPAccess instance profile. Per the BDSP licence no signal data is written
anywhere: the only thing this emits is one row of numbers per patient.

2026-08-21 PROVENANCE OF THIS FILE
This is t_ladder_up/worker_ladder.py, the worker that processed all 19,173
patients earlier today (itself the 2026-08-18 hypoxic-burden worker), with the
hypoxic-burden and event-anchored blocks REMOVED and two blocks added. Nothing
about how the signal is read, mapped to physical units, cleaned, or staged has
been touched, so every number here sits on the same footing as the frozen T90
column and the ladder run.

  1. DECILES OF ELAPSED RECORDING TIME. The SpO2 sample grid is cut into ten
     equal spans of elapsed recording time and each span reports the mean of its
     valid samples, the count of them, the elapsed minutes it covers, and the
     percent of them below 90. The spans partition the whole recording exactly,
     so two identities are checkable per record and are gated downstream:
         sum of the ten counts  = rec_n   (the whole-recording valid count)
         count-weighted mean of the ten means = rec_mean
     The SpO2 channel spans the full recording on every record of this cohort
     (spo2_dur_min / recording_dur_min = 1.000 on all 19,042 with usable
     oximetry in the ladder run), so a decile of the SpO2 index is a decile of
     elapsed recording time. Both durations are emitted so that stays checkable.

  2. PER-STAGE SATURATION. The stage array is already built on the SpO2 sampling
     grid by stage_on_grid, so a stage label sits on every oximetry sample with
     no resampling of either series. For each of the five codes (1=N3, 2=N2,
     3=N1, 4=REM, 5=Wake) the worker emits the mean of the valid SpO2 samples in
     that stage, the count of them, the stage's minutes (all samples, the
     definition the ladder run used for n3_min etc., so it is reproducible
     across runs), the minutes with valid oximetry, and the percent below 90.
     `staged_valid_n` is emitted independently as (stage > 0 & valid).sum(), so
     "the five stage counts sum to the staged valid count" is a real check
     rather than an identity of the analysis code.

VALIDATION TARGETS, all frozen and all local
  rec_mean   the whole-recording mean of the cleaned signal, against the frozen
             spo2_mean column. This is the run's primary validation axis and the
             one the deliverable is built on.
  prov_t90   the frozen T90 rule verbatim (band-limited resample to 200 Hz,
             float32 cast, keep [50,100], percent strictly below 90 over the
             whole recording), against frozen spo2_pct_below_90.
  stage mins against the frozen N1_pct/N2_pct/N3_pct/REM_pct times TST_min, and
             against the ladder run's own n1_min..rem_min (same code, same data,
             so any difference at all is a defect).

Everything below the primitives is the ladder worker's signal handling verbatim:
    signals/spo2 -> per-dataset affine (dig_min,dig_max -> phys_min,phys_max,
    clipped to the digital range) -> keep finite samples in [50,100].
Stage comes from annotations/expert_1/stage (starts/codes/ends|durations), run-length
expanded, with the pilot's whole-day offset repair.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import json
import multiprocessing as mp
import os
import queue as queue_mod
import sys
import time
import traceback
from multiprocessing.pool import ThreadPool

import numpy as np
import pandas as pd

AP = f"{paths.BDSP_S3_ACCESS_POINT}/PSG/bids"
SPO2_LO, SPO2_HI = 50.0, 100.0
TARGET_FS = 200.0                  # the declared rate the frozen T90 rule resamples to
SLEEP_CODES = (1, 2, 3, 4)
WAKE_CODE = 5
STAGE_NAMES = ((1, "n3"), (2, "n2"), (3, "n1"), (4, "rem"), (5, "wake"))
N_DECILES = 10
DAY = 86400.0
RECORD_TIMEOUT_S = 420
WORKERS = int(os.environ.get("HB_WORKERS", "16"))
PROBE = os.environ.get("HB_PROBE", "0") == "1"


# ------------------------------------------------------------------ primitives
def affine(raw, a):
    """Digital -> physical mapping from dataset attrs. None when unusable."""
    try:
        dmin, dmax = float(a["dig_min"]), float(a["dig_max"])
        pmin, pmax = float(a["phys_min"]), float(a["phys_max"])
    except Exception:
        return None, {}
    meta = {"dig_min": dmin, "dig_max": dmax, "phys_min": pmin, "phys_max": pmax}
    if not np.isfinite([dmin, dmax, pmin, pmax]).all() or dmax == dmin:
        return None, meta
    r = np.clip(raw, min(dmin, dmax), max(dmin, dmax))
    return pmin + (r - dmin) * (pmax - pmin) / (dmax - dmin), meta


def best_day_offset(ss, rec_s, max_days=2):
    if ss.size == 0 or not np.isfinite(rec_s) or rec_s <= 0:
        return 0
    best_k, best_in = 0, int(((ss >= -1.0) & (ss <= rec_s)).sum())
    for k in range(1, max_days + 1):
        s = ss - k * DAY
        n_in = int(((s >= -1.0) & (s <= rec_s)).sum())
        if n_in > best_in:
            best_k, best_in = k, n_in
    return best_k


def stage_on_grid(g, rec_s, fs, n_grid, out):
    """Sample-level stage array on an arbitrary grid (here: the SpO2 grid)."""
    ss = np.asarray(g["starts"][:], np.float64)
    cc = np.asarray(g["codes"][:], np.int32)
    try:
        dur = np.asarray(g["ends"][:], np.float64) - ss
    except KeyError:
        dur = np.asarray(g["durations"][:], np.float64)
    if ss.size == 0:
        out["stage_status"] = "stage_empty"
        return None
    bad = ~np.isfinite(dur) | (dur <= 0)
    out["n_stage_dur_repaired"] = int(bad.sum())
    dur = np.where(bad, 30.0, dur)
    dur = np.minimum(dur, rec_s)
    k = best_day_offset(ss, rec_s)
    out["stage_day_offset_days"] = int(k)
    ss = ss - k * DAY
    inside = (ss >= -1.0) & (ss <= rec_s)
    out["stage_out_of_range_frac"] = float((~inside).mean())

    stage = np.zeros(int(n_grid), np.int8)
    o = np.argsort(ss, kind="stable")
    i0 = np.rint(ss[o] * fs).astype(np.int64)
    i1 = np.rint((ss[o] + dur[o]) * fs).astype(np.int64)
    cs, ins = cc[o], inside[o]
    shifted = np.zeros(ss.size, bool) if k == 0 else np.ones(ss.size, bool)
    placed = 0
    for pass_shifted in (False, True):
        for a, b, c, ok, sh in zip(i0, i1, cs, ins, shifted):
            if bool(sh) != pass_shifted or not ok:
                continue
            if b <= 0 or a >= n_grid or c < 1 or c > 5:
                continue
            lo, hi = max(int(a), 0), min(int(b), int(n_grid))
            if hi <= lo:
                continue
            if pass_shifted:
                free = stage[lo:hi] == 0
                if not free.any():
                    continue
                stage[lo:hi] = np.where(free, np.int8(c), stage[lo:hi])
            else:
                stage[lo:hi] = c
            placed += 1
    out["n_stage_rows_placed"] = placed
    if placed == 0:
        out["stage_status"] = "stage_all_out_of_range"
        return None
    return stage


def describe(x, prefix, out):
    if x.size == 0:
        return
    out[f"{prefix}_n"] = int(x.size)
    out[f"{prefix}_mean"] = float(x.mean())
    out[f"{prefix}_nadir"] = float(x.min())
    out[f"{prefix}_median"] = float(np.median(x))
    out[f"{prefix}_t90"] = float((x < 90.0).mean() * 100.0)


def decile_profile(x, valid, fs, out):
    """
    Mean SpO2 in each of ten equal spans of elapsed recording time.

    The cut points are the ten equal spans of the SpO2 sample index, which is
    elapsed recording time on this archive's files. Each span reports the mean
    over its VALID samples only ([50,100] after the affine, the frozen rule), so
    a span that is all dropout reports no mean rather than a fabricated one. The
    spans partition the recording exactly, so the counts sum to the
    whole-recording valid count and the count-weighted mean of the means is the
    whole-recording mean.
    """
    n = int(x.size)
    if n < N_DECILES:
        out["decile_status"] = "too_few_samples"
        return
    edges = np.rint(np.linspace(0, n, N_DECILES + 1)).astype(np.int64)
    got = 0
    for i in range(N_DECILES):
        lo, hi = int(edges[i]), int(edges[i + 1])
        tag = f"{i + 1:02d}"
        out[f"dec_min_{tag}"] = float((hi - lo) / fs / 60.0)
        v = valid[lo:hi]
        cnt = int(v.sum())
        out[f"dec_n_{tag}"] = cnt
        if cnt == 0:
            continue
        seg = x[lo:hi][v]
        out[f"dec_mean_{tag}"] = float(seg.mean())
        out[f"dec_t90_{tag}"] = float((seg < 90.0).mean() * 100.0)
        got += 1
    out["n_deciles_with_data"] = got
    out["decile_status"] = "ok" if got == N_DECILES else "partial"


def stage_profile(x, valid, stage, fs, out):
    """Mean SpO2, minutes and percent below 90 within each scored stage."""
    per_min = 1.0 / fs / 60.0
    scored = stage > 0
    out["staged_min"] = float(scored.sum() * per_min)
    out["staged_valid_n"] = int((scored & valid).sum())      # emitted independently
    for code, nm in STAGE_NAMES:
        m = stage == code
        out[f"{nm}_min"] = float(m.sum() * per_min)
        mv = m & valid
        cnt = int(mv.sum())
        out[f"{nm}_n"] = cnt
        out[f"{nm}_valid_min"] = float(cnt * per_min)
        if cnt == 0:
            continue
        seg = x[mv]
        out[f"{nm}_mean"] = float(seg.mean())
        out[f"{nm}_t90"] = float((seg < 90.0).mean() * 100.0)
    sl = np.isin(stage, SLEEP_CODES)
    out["tst_min"] = float(sl.sum() * per_min)
    m = sl & valid
    out["sleep_valid_min"] = float(m.sum() * per_min)
    out["sleep_spo2_coverage"] = float(m.sum() / sl.sum()) if sl.sum() else np.nan
    if m.any():
        describe(x[m], "sleep", out)


def walk(h5, out_names, limit=400):
    def visit(name, obj):
        if len(out_names) >= limit:
            return
        try:
            shp = getattr(obj, "shape", None)
            out_names.append(f"{name}{'' if shp is None else ' ' + str(tuple(shp))}")
        except Exception:
            out_names.append(name)
    h5.visititems(visit)


# ------------------------------------------------------------------ per record
def extract_one(task, q):
    import h5py
    import s3fs
    rec = {"BDSPPatientID": int(task["BDSPPatientID"]), "site": task["site"],
           "sub": task["sub"], "ses": str(task["ses"]), "status": "ok"}
    t0 = time.time()
    try:
        fsys = s3fs.S3FileSystem(anon=False, default_block_size=4 * 1024 * 1024)
        key = (f"{AP}/{task['site']}/{task['sub']}/ses-{task['ses']}/eeg/"
               f"{task['sub']}_ses-{task['ses']}.h5")
        with fsys.open(key, "rb") as fh:
            with h5py.File(fh, "r") as h5:
                rec["meas_date"] = str(h5.attrs.get("meas_date", ""))
                try:
                    file_fs = float(h5.attrs["sampling_rate"])
                except Exception:
                    file_fs = np.nan

                if PROBE:
                    names = []
                    walk(h5, names)
                    rec["h5_structure"] = names

                sig = h5.get("signals")
                names = {k.lower(): k for k in sig.keys()} if sig is not None else {}
                ch = names.get("spo2")
                if ch is None:
                    for lk, orig in names.items():
                        if any(t in lk for t in ("spo2", "sao2", "osat")):
                            ch = orig
                            break
                if ch is None:
                    rec["status"] = "no_spo2_channel"
                    q.put(rec)
                    return
                rec["spo2_channel"] = ch
                d = sig[ch]

                raw = np.asarray(d[:]).ravel().astype(np.float64)
                rec["spo2_n_samples"] = int(raw.size)
                try:
                    fs = float(d.attrs["fs"])
                except Exception:
                    fs = np.nan
                if not np.isfinite(fs) or fs <= 0:
                    fs = file_fs if (np.isfinite(file_fs) and file_fs > 0) else TARGET_FS
                    rec["spo2_fs_source"] = "file_attr_or_default"
                else:
                    rec["spo2_fs_source"] = "dataset_attr"
                rec["spo2_fs"] = float(fs)
                rec["spo2_dur_min"] = float(raw.size / fs / 60.0)

                x, meta = affine(raw, d.attrs)
                rec["affine_used"] = x is not None
                if x is None:
                    x = raw
                for kk, vv in meta.items():
                    rec[f"spo2_{kk}"] = vv
                del raw

                valid = np.isfinite(x) & (x >= SPO2_LO) & (x <= SPO2_HI)
                rec["spo2_valid_frac"] = float(valid.mean()) if x.size else np.nan
                xv = x[valid]
                describe(xv, "rec", rec)              # rec_mean validates against spo2_mean

                # ---- recording length from the longest channel
                rec_s = 0.0
                for kk in names.values():
                    dd = sig[kk]
                    try:
                        cfs = float(dd.attrs["fs"])
                    except Exception:
                        cfs = file_fs if np.isfinite(file_fs) else TARGET_FS
                    if cfs and np.isfinite(cfs) and cfs > 0:
                        rec_s = max(rec_s, dd.shape[0] / cfs)
                rec["recording_dur_min"] = float(rec_s / 60.0)

                # ---- frozen-T90 validation variant (resample to 200 Hz, whole recording)
                n_out = int(round(x.size * TARGET_FS / fs))
                if x.size >= 32 and 32 <= n_out <= 60_000_000:
                    if abs(TARGET_FS / fs - 1.0) < 1e-12:
                        y = x
                    else:
                        from scipy.signal import resample
                        y = resample(x, n_out)
                    y = np.asarray(y, dtype=np.float32).astype(np.float64)
                    yv = y[np.isfinite(y) & (y >= SPO2_LO) & (y <= SPO2_HI)]
                    if yv.size >= 10:
                        rec["prov_t90"] = float((yv < 90.0).mean() * 100.0)
                        rec["prov_mean"] = float(yv.mean())
                        rec["prov_n_valid"] = int(yv.size)
                    del y, yv

                # ---- deliverable 1: ten deciles of elapsed recording time
                decile_profile(x, valid, fs, rec)

                # ---- stage on the SpO2 grid
                stage = None
                if "annotations" in h5 and "stage" in h5["annotations"] and \
                        hasattr(h5["annotations"]["stage"], "shape"):
                    rec["stage_layout"] = "sample_array"
                    st = np.asarray(h5["annotations"]["stage"][:]).ravel()
                    if st.size and np.isfinite(file_fs) and file_fs > 0:
                        idx = np.clip((np.arange(x.size) * (file_fs / fs)).astype(np.int64),
                                      0, st.size - 1)
                        stage = np.nan_to_num(st[idx], nan=0).astype(np.int8)
                elif "annotations/expert_1" in h5 and "stage" in h5["annotations/expert_1"]:
                    rec["stage_layout"] = "expert_1_events"
                    stage = stage_on_grid(h5["annotations/expert_1"]["stage"],
                                          rec_s, fs, x.size, rec)
                else:
                    rec["stage_layout"] = "none_found"

                rec["rec_valid_min"] = float(valid.sum() / fs / 60.0)

                # ---- deliverable 2: saturation within each stage
                if stage is not None:
                    stage_profile(x, valid, stage, fs, rec)
                del x, valid, xv
    except Exception as e:
        rec["status"] = "error"
        rec["error"] = f"{type(e).__name__}: {e}"
        rec["trace"] = traceback.format_exc()[-500:]
    rec["elapsed_s"] = round(time.time() - t0, 2)
    q.put(rec)


def run_with_timeout(task):
    ctx = mp.get_context("spawn")
    q = ctx.Queue()
    p = ctx.Process(target=extract_one, args=(task, q))
    p.start()
    try:
        rec = q.get(timeout=RECORD_TIMEOUT_S)
        p.join(15)
        if p.is_alive():
            p.kill()
            p.join(5)
        return rec
    except queue_mod.Empty:
        pass
    if p.is_alive():
        p.terminate()
        p.join(10)
        if p.is_alive():
            p.kill()
            p.join(5)
        return {"BDSPPatientID": int(task["BDSPPatientID"]), "site": task["site"],
                "sub": task["sub"], "ses": str(task["ses"]),
                "status": "timeout_or_crash"}
    try:
        return q.get(timeout=5)
    except Exception:
        return {"BDSPPatientID": int(task["BDSPPatientID"]), "site": task["site"],
                "sub": task["sub"], "ses": str(task["ses"]), "status": "no_result"}


def main():
    inp = sys.argv[1]
    out_jsonl = sys.argv[2]
    shard = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    nshard = int(sys.argv[4]) if len(sys.argv) > 4 else 1

    tasks = pd.read_csv(inp)
    tasks = tasks.iloc[shard::nshard].to_dict("records")

    done = set()
    if os.path.exists(out_jsonl):
        with open(out_jsonl) as f:
            for line in f:
                try:
                    r = json.loads(line)
                    if r.get("status") in ("ok", "no_spo2_channel"):
                        done.add((int(r["BDSPPatientID"]), str(r["ses"])))
                except Exception:
                    pass
    todo = [t for t in tasks
            if (int(t["BDSPPatientID"]), str(t["ses"])) not in done]
    print(f"shard {shard}/{nshard}: {len(tasks)} assigned, {len(done)} already done, "
          f"{len(todo)} to run, {WORKERS} concurrent", flush=True)

    t0 = time.time()
    fh = open(out_jsonl, "a", buffering=1)
    n_ok = 0
    with ThreadPool(WORKERS) as tp:
        for i, rec in enumerate(tp.imap_unordered(run_with_timeout, todo)):
            fh.write(json.dumps(rec, default=str) + "\n")
            n_ok += int(rec.get("status") == "ok")
            if (i + 1) % 100 == 0 or i + 1 == len(todo):
                el = time.time() - t0
                rate = (i + 1) / el if el > 0 else 0
                print(f"[{i+1}/{len(todo)}] ok {n_ok}  {rate:.2f} rec/s  "
                      f"eta {(len(todo)-i-1)/max(rate,1e-9)/60:.1f} min", flush=True)
    fh.close()
    print(f"done in {(time.time()-t0)/60:.1f} min, ok {n_ok}/{len(todo)}", flush=True)


if __name__ == "__main__":
    main()
