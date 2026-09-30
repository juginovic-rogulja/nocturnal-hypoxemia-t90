#!/usr/bin/env python3
"""
Threshold-ladder extraction worker (T80/T92/T95). Runs on EC2 (us-east-1) under
the EC2BDSPAccess instance profile. Per the BDSP licence no signal data is
written anywhere: the only thing this emits is one row of numbers per patient.

2026-08-21 PROVENANCE OF THIS FILE
This is t85_analysis/worker_t85.py, the worker that processed all 19,173
patients earlier today (which itself is the 2026-08-18 hypoxic-burden worker
plus the T85 output lines), with exactly ONE addition and no other change:

  the frozen-T90 validation block also reports `prov_t80`, `prov_t92` and
  `prov_t95`, the same rule at thresholds 80, 92 and 95, plus `prov_a95`, the
  percent of valid samples AT OR ABOVE 95 (preserved saturation), computed
  independently so that a95 + t95 = 100 is a checkable identity per record.
  Same affine, same 200 Hz band-limited resample, same float32 cast, same
  [50,100] keep band, same whole-recording denominator. prov_t90 and prov_t88
  still validate the rule against the two frozen columns, and prov_t85 must
  reproduce the 2026-08-21 T85 run to 0.01 pp because it is the same code on
  the same data.

The prior t85 header, kept because everything it says still holds:
This is hypoxic_burden/worker_burden.py, the worker that processed all 19,173
patients on 2026-08-18, with exactly two additions and no other change:

  1. describe() also reports `{prefix}_t85`, percent strictly below 85.
  2. the frozen-T90 validation block also reports `prov_t88`, `prov_t85` and
     `prov_n_valid`.

`prov_t85` is the deliverable. It is produced by the SAME rule, on the SAME
whole-recording denominator, as the frozen `spo2_pct_below_90` column
(numbers/t90_provenance_reference.py): affine to physical, band-limited resample
to 200 Hz, float32 cast, keep finite samples in [50,100] inclusive, then percent
strictly below the threshold over the ENTIRE recording with no sleep
restriction. So T85 is comparable to the paper's T90 by construction, not by
assumption. `prov_t88` exists only because the frozen data carries
`spo2_pct_below_88`, which gives the run a second independent validation target.

Everything below this line is the 2026-08-18 worker verbatim.

Signal handling is transcribed from the 2026-08-13 raw pilot
(raw_pilot_rederive/worker_rederive.py), which reproduced the frozen T90 column
for 23/25 patients:
    signals/spo2 -> per-dataset affine (dig_min,dig_max -> phys_min,phys_max,
    clipped to the digital range) -> keep finite samples in [50,100].
    T90 provenance variant additionally band-limit-resamples to 200 Hz and casts
    float32 before filtering; that variant is what the frozen column equals, so it
    is carried here purely as the per-patient validation check.
Stage comes from annotations/expert_1/stage (starts/codes/ends|durations), run-length
expanded, codes 1=N3 2=N2 3=N1 4=REM 5=Wake, with the pilot's whole-day offset repair.
Unlike the pilot the stage array is built on the SpO2 sampling grid, so a stage label
sits on every oximetry sample without any resampling of either series.

BURDEN DEFINITIONS (all integrals use only samples that are both scored sleep and
valid oximetry, and are normalised by that same amount of time, so dropout cannot
inflate or deflate a rate):

  A. Primary, event-free area under a night-specific baseline
       B    = 95th percentile of the night's cleaned SpO2 (whole recording)
       area = sum over sleep samples of max(0, B - SpO2) / fs / 60      [%-min]
       HB   = area / (sleep-with-valid-SpO2 minutes / 60)              [%-min/h]
     hb_p95_recbase_pctmin_per_h. A sleep-only-baseline twin is also emitted.

  B. Fixed threshold, area below 90%
       area = sum over sleep samples of max(0, 90 - SpO2) / fs / 60
       HB   = area / (sleep-with-valid-SpO2 minutes / 60)
     hb_below90_pctmin_per_h. An 88% twin is emitted too.

  C. Event-anchored (only where respiratory event annotations exist)
       per event, pre-event baseline = mean cleaned SpO2 over the 100 s before onset,
       window = onset to onset + 90 s, area = integral of max(0, baseline - SpO2).
       HB   = sum of event areas / (sleep-with-valid-SpO2 hours)
     A fixed-window simplification of Azarbarzin's individual-ensemble window.

Whole-recording twins of A and B are emitted so records without usable staging are
not simply lost.
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
DAY = 86400.0
RECORD_TIMEOUT_S = 420
WORKERS = int(os.environ.get("HB_WORKERS", "16"))
PROBE = os.environ.get("HB_PROBE", "0") == "1"

EVENT_HINTS = ("apnea", "hypopnea", "resp", "desat", "event", "arousal", "sao2", "spo2")
RESP_HINTS = ("apnea", "hypopnea", "resp")
PRE_EVENT_S = 100.0
EVENT_WINDOW_S = 90.0


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


def area_stats(x, fs, minutes, prefix, out, baseline):
    """Integral of max(0, baseline - x) over the samples given, as %-min and %-min/h."""
    if x.size == 0 or not np.isfinite(baseline) or minutes <= 0:
        return
    deficit = baseline - x
    np.maximum(deficit, 0.0, out=deficit)
    pctmin = float(deficit.sum() / fs / 60.0)
    out[f"{prefix}_pctmin"] = pctmin
    out[f"{prefix}_pctmin_per_h"] = float(pctmin / (minutes / 60.0))


def describe(x, prefix, out):
    if x.size == 0:
        return
    out[f"{prefix}_n"] = int(x.size)
    out[f"{prefix}_mean"] = float(x.mean())
    out[f"{prefix}_nadir"] = float(x.min())
    out[f"{prefix}_p95"] = float(np.percentile(x, 95))
    out[f"{prefix}_median"] = float(np.median(x))
    out[f"{prefix}_t90"] = float((x < 90.0).mean() * 100.0)
    out[f"{prefix}_t88"] = float((x < 88.0).mean() * 100.0)
    out[f"{prefix}_t85"] = float((x < 85.0).mean() * 100.0)


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


def event_burden(h5, x_full, valid_full, fs, sleep_hours, rec_s, out):
    """
    Azarbarzin-style event-anchored burden, fixed-window simplification.
    Looks for respiratory event annotations under annotations/ (any expert group).
    """
    if "annotations" not in h5 or sleep_hours <= 0:
        return
    groups = []
    try:
        for gname, gobj in h5["annotations"].items():
            if hasattr(gobj, "keys"):
                for k in gobj.keys():
                    groups.append((f"annotations/{gname}/{k}", k.lower()))
            else:
                groups.append((f"annotations/{gname}", gname.lower()))
    except Exception:
        return
    resp = [(p, k) for p, k in groups if any(h in k for h in RESP_HINTS)]
    out["event_dataset_candidates"] = ";".join(sorted({p for p, _ in groups}))[:900]
    if not resp:
        out["event_burden_status"] = "no_respiratory_annotation"
        return

    # Two layouts exist in this archive: a single `resp` dataset, or `resp_3` plus `resp_4`,
    # which are the 3% and 4% desaturation hypopnea rules and therefore NEST. Concatenating
    # them would count the same event twice and would do so only for the files that carry
    # both, which is a between-file artefact, not physiology. Onsets are de-duplicated to
    # 0.1 s, so the union collapses to the more inclusive set. Per-dataset counts are kept
    # so the assumption is checkable after the fact rather than merely asserted.
    per_ds = {}
    for path, key in resp:
        try:
            g = h5[path]
            ss = np.asarray(g["starts"][:], np.float64)
            try:
                dd = np.asarray(g["ends"][:], np.float64) - ss
            except KeyError:
                dd = np.asarray(g["durations"][:], np.float64)
            per_ds[key] = (ss, dd)
            out[f"n_events_{key}"] = int(ss.size)
        except Exception:
            continue
    if not per_ds:
        out["event_burden_status"] = "respiratory_annotation_unreadable"
        return
    out["resp_layout"] = "+".join(sorted(per_ds))

    def score(starts_list):
        ss = np.concatenate(starts_list) if starts_list else np.zeros(0)
        if ss.size == 0:
            return None
        k = best_day_offset(ss, rec_s)
        ss = ss - k * DAY
        ss = ss[np.isfinite(ss) & (ss >= 0) & (ss <= rec_s)]
        if ss.size == 0:
            return None
        ss = np.unique(np.round(ss, 1))          # de-duplicate the nested rules
        n = x_full.size
        total, used = 0.0, 0
        pre_n = int(fs * 10)
        for i, s in enumerate(ss):
            a = int(round(s * fs))
            pre_lo = max(0, int(round((s - PRE_EVENT_S) * fs)))
            if a - pre_lo < pre_n:
                continue
            pv = valid_full[pre_lo:a]
            if pv.sum() < pre_n:
                continue
            base = float(x_full[pre_lo:a][pv].mean())
            # the window stops at the next onset, so back-to-back events in severe disease
            # cannot have their desaturations integrated two and three times over
            end_s = s + EVENT_WINDOW_S
            if i + 1 < ss.size:
                end_s = min(end_s, ss[i + 1])
            b = min(n, int(round(end_s * fs)))
            if b <= a:
                continue
            wv = valid_full[a:b]
            if not wv.any():
                continue
            d = base - x_full[a:b][wv]
            np.maximum(d, 0.0, out=d)
            total += float(d.sum() / fs / 60.0)
            used += 1
        return total, used, int(ss.size)

    r = score([v[0] for v in per_ds.values()])
    if r is None:
        out["event_burden_status"] = "no_events_in_range"
        out["n_resp_events"] = 0
        return
    total, used, n_ev = r
    out["n_resp_events"] = n_ev
    out["n_resp_events_scored"] = used
    out["hb_event_pctmin"] = total
    out["hb_event_pctmin_per_h"] = float(total / sleep_hours)
    out["event_burden_status"] = "ok"

    # the strict 4% rule alone, where the file carries it, as the sensitivity axis
    strict = [v[0] for k, v in per_ds.items() if k.endswith("_4")]
    if strict:
        r4 = score(strict)
        if r4 is not None:
            out["hb_event4_pctmin"] = r4[0]
            out["hb_event4_pctmin_per_h"] = float(r4[0] / sleep_hours)
            out["n_resp_events_4"] = r4[2]


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
                describe(xv, "rec", rec)

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
                        # T90 is the frozen column's own rule and is carried as the
                        # per-patient validation check. T88 has a frozen twin
                        # (spo2_pct_below_88) so it validates the run a second,
                        # independent way. T85 is the new exposure and is produced by
                        # the SAME rule on the SAME denominator, so it is comparable to
                        # the frozen T90 by construction rather than by assumption.
                        # a95 is the preserved-saturation complement of t95, emitted
                        # INDEPENDENTLY ((yv >= 95).mean(), not 100 - t95) so the
                        # identity a95 + t95 = 100 is checkable on every record
                        # rather than true by definition of the analysis code.
                        rec["prov_a95"] = float((yv >= 95.0).mean() * 100.0)
                        rec["prov_t95"] = float((yv < 95.0).mean() * 100.0)
                        rec["prov_t92"] = float((yv < 92.0).mean() * 100.0)
                        rec["prov_t90"] = float((yv < 90.0).mean() * 100.0)
                        rec["prov_t88"] = float((yv < 88.0).mean() * 100.0)
                        rec["prov_t85"] = float((yv < 85.0).mean() * 100.0)
                        rec["prov_t80"] = float((yv < 80.0).mean() * 100.0)
                        rec["prov_mean"] = float(yv.mean())
                        rec["prov_n_valid"] = int(yv.size)
                    del y, yv

                # ---- stage on the SpO2 grid
                stage = None
                if "annotations" in h5 and "stage" in h5["annotations"] and \
                        hasattr(h5["annotations"]["stage"], "shape"):
                    rec["stage_layout"] = "sample_array"
                    st = np.asarray(h5["annotations"]["stage"][:]).ravel()
                    if st.size and np.isfinite(file_fs) and file_fs > 0:
                        # legacy layout is on the file rate; map onto the SpO2 grid
                        idx = np.clip((np.arange(x.size) * (file_fs / fs)).astype(np.int64),
                                      0, st.size - 1)
                        stage = np.nan_to_num(st[idx], nan=0).astype(np.int8)
                elif "annotations/expert_1" in h5 and "stage" in h5["annotations/expert_1"]:
                    rec["stage_layout"] = "expert_1_events"
                    stage = stage_on_grid(h5["annotations/expert_1"]["stage"],
                                          rec_s, fs, x.size, rec)
                else:
                    rec["stage_layout"] = "none_found"

                per_min = 1.0 / fs / 60.0
                base_rec = rec.get("rec_p95", np.nan)

                # ---- whole-recording burden (fallback, no sleep restriction)
                rec_min = float(valid.sum() * per_min)
                rec["rec_valid_min"] = rec_min
                area_stats(xv, fs, rec_min, "hb_rec_p95", rec, base_rec)
                area_stats(xv, fs, rec_min, "hb_rec_below90", rec, 90.0)

                sleep_hours = 0.0
                if stage is not None:
                    sl = np.isin(stage, SLEEP_CODES)
                    rec["tst_min"] = float(sl.sum() * per_min)
                    rec["staged_min"] = float((stage > 0).sum() * per_min)
                    rec["wake_min"] = float((stage == WAKE_CODE).sum() * per_min)
                    for code, nm in ((1, "n3"), (2, "n2"), (3, "n1"), (4, "rem")):
                        rec[f"{nm}_min"] = float((stage == code).sum() * per_min)
                    m = sl & valid
                    sleep_min = float(m.sum() * per_min)
                    rec["sleep_valid_min"] = sleep_min
                    rec["sleep_spo2_coverage"] = (
                        float(m.sum() / sl.sum()) if sl.sum() else np.nan)
                    xs = x[m]
                    if xs.size:
                        describe(xs, "sleep", rec)
                        base_sleep = float(np.percentile(xs, 95))
                        area_stats(xs, fs, sleep_min, "hb_p95_recbase", rec, base_rec)
                        area_stats(xs, fs, sleep_min, "hb_p95_sleepbase", rec, base_sleep)
                        area_stats(xs, fs, sleep_min, "hb_below90", rec, 90.0)
                        area_stats(xs, fs, sleep_min, "hb_below88", rec, 88.0)
                        sleep_hours = sleep_min / 60.0
                    del xs

                # ---- event-anchored
                try:
                    event_burden(h5, x, valid, fs, sleep_hours, rec_s, rec)
                except Exception as e:
                    rec["event_burden_status"] = f"error:{type(e).__name__}"
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
