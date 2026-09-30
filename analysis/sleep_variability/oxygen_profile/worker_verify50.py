#!/usr/bin/env python3
"""
INDEPENDENT verification worker: 50 patients, a second implementation written from scratch.

This file deliberately does NOT import worker_oxyprofile. It shares the definitions (clean to
[50,100] after the digital-to-physical affine, ten equal spans of elapsed recording time, stage
codes 1=N3 2=N2 3=N1 4=REM 5=Wake, stage minutes counted over scored samples) and nothing else.
Where a mechanism could hide a bug it uses a different one on purpose:

  the affine        written as scale-and-shift rather than the ratio form
  the ten spans     each sample is LABELLED by integer division, (j * 10) // n, and the means
                    come from np.bincount sums. The fleet worker instead computes edge indices
                    with np.rint(np.linspace(...)) and slices. The two agree except possibly at
                    a single sample per boundary, which moves a decile mean by about 1e-4 of a
                    percentage point on a 600,000-sample recording.
  the hypnogram     events are sorted and painted with an explicit index loop of its own,
                    with the same whole-day offset repair rule found by scanning k in 0..2
  the aggregates    per-stage means come from np.bincount over the stage code, not from five
                    boolean masks

So a shared bug cannot self-confirm: an error in either implementation shows up as a numeric
disagreement in the comparison table, which is the gate.

Per the BDSP licence, signals never leave the account. This emits one row of numbers per
patient, exactly as the fleet worker does.

    python3 worker_verify50.py tasks_verify50.csv out.jsonl [shard] [nshard]
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import json
import os
import sys
import time
import traceback
from multiprocessing.pool import ThreadPool

import numpy as np
import pandas as pd

ACCESS_POINT = f"{paths.BDSP_S3_ACCESS_POINT}/PSG/bids"
KEEP_LO, KEEP_HI = 50.0, 100.0
NSPAN = 10
CODE_NAME = {1: "n3", 2: "n2", 3: "n1", 4: "rem", 5: "wake"}
SECONDS_PER_DAY = 86400.0
NTHREAD = int(os.environ.get("HB_WORKERS", "8"))


def to_physical(digital, attrs):
    """Digital counts to saturation in percent, as scale and shift."""
    dlo, dhi = float(attrs["dig_min"]), float(attrs["dig_max"])
    plo, phi = float(attrs["phys_min"]), float(attrs["phys_max"])
    if not all(np.isfinite([dlo, dhi, plo, phi])) or dhi == dlo:
        return None
    scale = (phi - plo) / (dhi - dlo)
    lo, hi = (dlo, dhi) if dlo <= dhi else (dhi, dlo)
    bounded = np.clip(digital, lo, hi)
    return (bounded - dlo) * scale + plo


def span_labels(n):
    """Which of the ten spans of elapsed recording time each sample belongs to."""
    return np.minimum((np.arange(n, dtype=np.int64) * NSPAN) // n, NSPAN - 1)


def group_mean(labels, values, keep, ngroup):
    """Mean of `values` per label, over the samples `keep` selects. Counts come back too."""
    lab = labels[keep]
    val = values[keep]
    counts = np.bincount(lab, minlength=ngroup).astype(np.int64)
    sums = np.bincount(lab, weights=val, minlength=ngroup)
    means = np.full(ngroup, np.nan, float)
    hit = counts > 0
    means[hit] = sums[hit] / counts[hit]
    return means, counts


def hypnogram(group, n_samples, rate, record_seconds):
    """
    Sample-level stage codes on the oximetry grid, painted from the scored events.

    Written independently of the fleet worker: events are sorted once, the whole-day offset is
    chosen by scanning candidates, and the painting is a single explicit loop that writes into
    free positions when the events had to be shifted (so the earliest scoring wins) and
    straight over the top when they did not.
    """
    starts = np.asarray(group["starts"][:], float)
    codes = np.asarray(group["codes"][:], int)
    if "ends" in group:
        lengths = np.asarray(group["ends"][:], float) - starts
    else:
        lengths = np.asarray(group["durations"][:], float)
    if starts.size == 0:
        return None, {"stage_note": "no_events"}
    lengths = np.where(np.isfinite(lengths) & (lengths > 0), lengths, 30.0)
    lengths = np.minimum(lengths, record_seconds)

    offsets = [k for k in range(3)]
    inside = [int(np.sum((starts - k * SECONDS_PER_DAY >= -1.0)
                         & (starts - k * SECONDS_PER_DAY <= record_seconds)))
              for k in offsets]
    k = int(np.argmax(inside))
    starts = starts - k * SECONDS_PER_DAY

    order = np.argsort(starts, kind="stable")
    stage = np.zeros(n_samples, np.int8)
    painted = 0
    for pos in order:
        code = int(codes[pos])
        if code < 1 or code > 5:
            continue
        begin = int(np.rint(starts[pos] * rate))
        stop = int(np.rint((starts[pos] + lengths[pos]) * rate))
        if starts[pos] < -1.0 or starts[pos] > record_seconds:
            continue
        begin = max(begin, 0)
        stop = min(stop, n_samples)
        if stop <= begin:
            continue
        if k > 0:
            window = stage[begin:stop]
            empty = window == 0
            if not empty.any():
                continue
            stage[begin:stop] = np.where(empty, np.int8(code), window)
        else:
            stage[begin:stop] = np.int8(code)
        painted += 1
    if painted == 0:
        return None, {"stage_note": "nothing_painted"}
    return stage, {"stage_note": "ok", "v_day_offset": k, "v_events_painted": painted}


def one_patient(task):
    import h5py
    import s3fs
    row = {"BDSPPatientID": int(task["BDSPPatientID"]), "site": task["site"],
           "sub": task["sub"], "ses": str(task["ses"]), "status": "ok"}
    started = time.time()
    try:
        store = s3fs.S3FileSystem(anon=False, default_block_size=4 * 1024 * 1024)
        path = (f"{ACCESS_POINT}/{task['site']}/{task['sub']}/ses-{task['ses']}/eeg/"
                f"{task['sub']}_ses-{task['ses']}.h5")
        with store.open(path, "rb") as handle:
            with h5py.File(handle, "r") as f:
                channels = dict(f["signals"].items())
                pick = None
                for name in channels:
                    if name.lower() == "spo2":
                        pick = name
                        break
                if pick is None:
                    for name in channels:
                        if "spo2" in name.lower() or "sao2" in name.lower():
                            pick = name
                            break
                if pick is None:
                    row["status"] = "no_spo2_channel"
                    return row
                dset = channels[pick]
                row["v_channel"] = pick

                try:
                    rate = float(dset.attrs["fs"])
                except Exception:
                    rate = float(f.attrs.get("sampling_rate", 200.0))
                if not np.isfinite(rate) or rate <= 0:
                    rate = float(f.attrs.get("sampling_rate", 200.0))
                row["v_fs"] = rate

                counts_raw = np.asarray(dset[:], float).reshape(-1)
                n = counts_raw.size
                row["v_n_samples"] = int(n)
                sat = to_physical(counts_raw, dset.attrs)
                if sat is None:
                    sat = counts_raw
                    row["v_affine"] = False
                else:
                    row["v_affine"] = True

                usable = np.isfinite(sat) & (sat >= KEEP_LO) & (sat <= KEEP_HI)
                row["v_rec_n"] = int(usable.sum())
                if usable.any():
                    row["v_rec_mean"] = float(sat[usable].mean())

                # --- ten spans of elapsed recording time
                labels = span_labels(n)
                means, counts = group_mean(labels, sat, usable, NSPAN)
                for i in range(NSPAN):
                    row[f"v_dec_n_{i+1:02d}"] = int(counts[i])
                    if counts[i] > 0:
                        row[f"v_dec_mean_{i+1:02d}"] = float(means[i])
                row["v_deciles_with_data"] = int((counts > 0).sum())

                # --- recording length, taken from the longest channel
                longest = 0.0
                for name, chan in channels.items():
                    try:
                        crate = float(chan.attrs["fs"])
                    except Exception:
                        crate = float(f.attrs.get("sampling_rate", 200.0))
                    if np.isfinite(crate) and crate > 0:
                        longest = max(longest, chan.shape[0] / crate)
                row["v_record_min"] = float(longest / 60.0)

                # --- stages
                if "annotations/expert_1" in f and "stage" in f["annotations/expert_1"]:
                    stage, note = hypnogram(f["annotations/expert_1"]["stage"], n, rate, longest)
                    row.update(note)
                else:
                    stage, row["stage_note"] = None, "no_expert_1_stage"
                if stage is not None:
                    smean, scount = group_mean(stage, sat, usable, 6)
                    scored = np.bincount(stage, minlength=6).astype(np.int64)
                    for code, nm in CODE_NAME.items():
                        row[f"v_{nm}_min"] = float(scored[code] / rate / 60.0)
                        row[f"v_{nm}_n"] = int(scount[code])
                        if scount[code] > 0:
                            row[f"v_{nm}_mean"] = float(smean[code])
                    row["v_tst_min"] = float(scored[1:5].sum() / rate / 60.0)
    except Exception as exc:
        row["status"] = "error"
        row["error"] = f"{type(exc).__name__}: {exc}"
        row["trace"] = traceback.format_exc()[-400:]
    row["v_elapsed_s"] = round(time.time() - started, 2)
    return row


def main():
    tasks = pd.read_csv(sys.argv[1]).to_dict("records")
    out_path = sys.argv[2]
    print(f"{len(tasks)} patients, {NTHREAD} concurrent", flush=True)
    started = time.time()
    with open(out_path, "a", buffering=1) as fh:
        with ThreadPool(NTHREAD) as pool:
            for i, row in enumerate(pool.imap_unordered(one_patient, tasks), start=1):
                fh.write(json.dumps(row, default=str) + "\n")
                if i % 10 == 0 or i == len(tasks):
                    print(f"[{i}/{len(tasks)}] {(time.time()-started)/60:.1f} min", flush=True)
    print(f"done in {(time.time()-started)/60:.1f} min", flush=True)


if __name__ == "__main__":
    main()
