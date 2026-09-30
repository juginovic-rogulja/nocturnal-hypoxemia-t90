"""
EC2 worker for the v8 FULL-family pass (the v7 worker.py pattern, robust version): manifest CSV with site_id, bids_col,
SessionID, BDSPPatientID, _encoding_used, t_end_sec (empty = whole night). Each night: aws s3 cp the .h5 from the BDSP access
point to tmpfs, run scripts/extract_v8.py --one in ITS OWN SUBPROCESS with a hard timeout, delete the file, append the record
to a JSONL. Failed nights are retried once at the end with a quarter of the concurrency. Partials to S3 every PARTIAL_EVERY
nights and every 10 minutes; the heartbeat comes from user-data.
Environment: OUTPUT_BUCKET, OUTPUT_PREFIX, RUN_TAG, N_WORKERS, MANIFEST, SHARD (i/n), WORK_DIR, OUT_DIR, NIGHT_TIMEOUT (s,
default 1500), PARTIAL_EVERY (250), HB_DEFINITION (path), FAMILIES (full).
"""
from __future__ import annotations
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import json, os, subprocess, sys, threading, time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent
EXTRACT = HERE.parent / "scripts" / "extract_v8.py"
PY = sys.executable
S3_BASE = paths.S3_PSG_BASE
OUTPUT_BUCKET = os.environ.get("OUTPUT_BUCKET", paths.S3_OUTPUT_BUCKET)
OUTPUT_PREFIX = os.environ.get("OUTPUT_PREFIX", "outputs/v8_2026-09/run")
RUN_TAG = os.environ.get("RUN_TAG", "run")
N_WORKERS = int(os.environ.get("N_WORKERS", str(os.cpu_count() or 4)))
WORK = Path(os.environ.get("WORK_DIR", "/opt/reex/work")); WORK.mkdir(parents=True, exist_ok=True)
OUT_DIR = Path(os.environ.get("OUT_DIR", "/opt/reex/out")); OUT_DIR.mkdir(parents=True, exist_ok=True)
PARTIAL_EVERY = int(os.environ.get("PARTIAL_EVERY", "250"))
NIGHT_TIMEOUT = int(os.environ.get("NIGHT_TIMEOUT", "1500"))
HB_DEF = os.environ.get("HB_DEFINITION", "")
FAMILIES = os.environ.get("FAMILIES", "full")
LOCK = threading.Lock()


def s3_upload(local, key):
    try:
        subprocess.run(["aws", "s3", "cp", str(local), f"s3://{OUTPUT_BUCKET}/{OUTPUT_PREFIX}/{key}", "--only-show-errors"], check=False, timeout=300)
    except Exception:
        pass


def _t_end(row):
    v = row.get("t_end_sec", "")
    try:
        f = float(v)
        return f if f == f else None
    except Exception:
        return None


def process(row):
    site, bids, ses = row["site_id"], row["bids_col"], int(row["SessionID"])
    stem = f"{bids}_ses-{ses}"; local = WORK / f"{stem}.h5"; jout = WORK / f"{stem}.json"
    url = f"{S3_BASE}/{site}/{bids}/ses-{ses}/eeg/{stem}.h5"
    t_end = _t_end(row)
    rec = {"BDSPPatientID": int(row["BDSPPatientID"]), "site_id": site, "stem": stem, "_encoding_used": row.get("_encoding_used", ""),
           "t_end_sec_manifest": t_end}
    t0 = time.time()
    try:
        r = None
        for attempt in range(3):
            try:
                r = subprocess.run(["aws", "s3", "cp", url, str(local), "--only-show-errors"], capture_output=True, text=True, timeout=600)
            except subprocess.TimeoutExpired:
                r = None
            if r is not None and r.returncode == 0 and local.exists():
                break
            time.sleep(5 * (attempt + 1))
        else:
            rec["status"] = "download_fail:" + (r.stderr[-200:] if r is not None else "timeout"); return rec
        rec["download_sec"] = round(time.time() - t0, 1); rec["file_mb"] = round(local.stat().st_size / 1e6, 1)
        t1 = time.time()
        cmd = [PY, str(EXTRACT), "--one", str(local), "--families", FAMILIES, "--json", str(jout)]
        if t_end is not None:
            cmd += ["--t-end-sec", str(t_end)]
        if HB_DEF:
            cmd += ["--hb-definition", HB_DEF]
        try:
            p = subprocess.run(cmd, capture_output=True, text=True, timeout=NIGHT_TIMEOUT,
                               env={**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"})
        except subprocess.TimeoutExpired:
            rec["status"] = f"extract_timeout_{NIGHT_TIMEOUT}s"; return rec
        if p.returncode != 0 or not jout.exists():
            rec["status"] = f"extract_exit_{p.returncode}:" + (p.stderr[-300:] if p.stderr else ""); return rec
        out = json.load(open(jout))
        if out.get("status") == "ok":
            out.pop("traceback", None)
        rec.update(out); rec["extract_wall_sec"] = round(time.time() - t1, 1)
    except Exception as e:
        rec["status"] = f"fail:{type(e).__name__}:{str(e)[:200]}"
    finally:
        for f in (local, jout):
            try:
                if f.exists():
                    f.unlink()
            except Exception:
                pass
    rec["total_sec"] = round(time.time() - t0, 1)
    return rec


def run_pass(rows, n_workers, fh, jsonl, counters, t0, label):
    last_partial = time.time()
    with ThreadPoolExecutor(max_workers=n_workers) as ex:
        futs = {ex.submit(process, r): r for r in rows}
        for i, fut in enumerate(as_completed(futs), 1):
            try:
                rec = fut.result()
            except Exception as e:
                r = futs[fut]; rec = {"BDSPPatientID": int(r["BDSPPatientID"]), "site_id": r["site_id"], "status": f"fail:{type(e).__name__}:{str(e)[:200]}"}
            with LOCK:
                fh.write(json.dumps(rec, default=str) + "\n"); fh.flush()
                counters["done"] += 1
                if str(rec.get("status", "")).startswith("ok"):
                    counters["ok"] += 1
                else:
                    counters["fail"] += 1; counters["failed_rows"].append(futs[fut])
            if i % 25 == 0 or i == len(rows):
                el = time.time() - t0
                print(f"[{label} {i}/{len(rows)}] {el/60:.1f} min ok={counters['ok']} fail={counters['fail']} {el/max(1,counters['done']):.2f} s/night wall", flush=True)
            if counters["done"] % PARTIAL_EVERY == 0 or time.time() - last_partial > 600:
                s3_upload(jsonl, f"{RUN_TAG}_partial.jsonl"); last_partial = time.time()


def main():
    df = pd.read_csv(os.environ.get("MANIFEST", "/opt/reex/manifest.csv"), low_memory=False)
    shard = os.environ.get("SHARD", "")
    if shard:
        i, n = [int(x) for x in shard.split("/")]; df = df.iloc[i::n].copy()
    done = set()
    resume = os.environ.get("RESUME_FROM", "")
    if resume and os.path.isfile(resume):
        for line in open(resume):
            try:
                d = json.loads(line)
                if str(d.get("status", "")).startswith("ok"):
                    done.add(int(d["BDSPPatientID"]))
            except Exception:
                pass
    rows = [r for r in df.to_dict("records") if int(r["BDSPPatientID"]) not in done]
    print(f"V8 {RUN_TAG}: {len(rows)} nights ({len(done)} resumed), {N_WORKERS} workers, families={FAMILIES}, shard={shard or 'all'}, timeout {NIGHT_TIMEOUT}s, hb_def={HB_DEF}", flush=True)
    jsonl = OUT_DIR / f"{RUN_TAG}.jsonl"
    if resume and os.path.isfile(resume) and str(jsonl) != resume:
        jsonl.write_text(open(resume).read())
    t0 = time.time(); counters = {"done": 0, "ok": 0, "fail": 0, "failed_rows": []}
    with open(jsonl, "a") as fh:
        run_pass(rows, N_WORKERS, fh, jsonl, counters, t0, "pass1")
        retry = list(counters["failed_rows"])
        if retry:
            print(f"retrying {len(retry)} failed nights with {max(2, N_WORKERS // 4)} workers", flush=True)
            counters["failed_rows"] = []
            run_pass(retry, max(2, N_WORKERS // 4), fh, jsonl, counters, t0, "retry")
    summary = {"run_tag": RUN_TAG, "n_in": len(rows), "n_records": counters["done"], "n_ok": counters["ok"], "n_fail_records": counters["fail"],
               "n_retried": len(retry), "n_still_failed": len(counters["failed_rows"]), "elapsed_min": (time.time() - t0) / 60,
               "n_workers": N_WORKERS, "sec_per_night_wall": (time.time() - t0) / max(1, len(rows)), "families": FAMILIES, "hb_definition": HB_DEF}
    (OUT_DIR / f"{RUN_TAG}_summary.json").write_text(json.dumps(summary, indent=2))
    s3_upload(jsonl, f"{RUN_TAG}_FINAL.jsonl"); s3_upload(OUT_DIR / f"{RUN_TAG}_summary.json", f"{RUN_TAG}_summary.json")
    print(f"DONE {RUN_TAG}: {summary}", flush=True)


if __name__ == "__main__":
    main()
