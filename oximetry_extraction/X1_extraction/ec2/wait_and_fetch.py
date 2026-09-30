"""
Poll an EC2 v8 run: heartbeat every 30 s, STALL after 10 min without nights_done moving, done when <tag>_FINAL.jsonl appears
(or the instance's TERMINATED marker appears without FINAL: then the LAST/partial jsonl is fetched and the run is marked
incomplete). Fetches FINAL jsonl, run/worker logs, summary and versions into X1/data/ with a provenance sidecar.
usage: wait_and_fetch.py <run_tag> [--max-wait-min N] [--stall-min 10]
"""
from __future__ import annotations
import argparse, json, os, subprocess, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))
import x1_paths as P
from provenance import write_sidecar


def s3(*args, timeout=120):
    return subprocess.run(["aws", "s3", *args], capture_output=True, text=True, timeout=timeout)


def exists(uri):
    return s3("ls", uri).returncode == 0 and bool(s3("ls", uri).stdout.strip())


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("tag"); ap.add_argument("--max-wait-min", type=float, default=360)
    ap.add_argument("--stall-min", type=float, default=10); ap.add_argument("--poll-sec", type=int, default=30)
    a = ap.parse_args()
    base = f"s3://{P.S3_BUCKET}/{P.S3_OUTPUT_PREFIX}/{a.tag}"
    t0 = time.time(); last_n = None; last_move = time.time(); status = "running"
    while True:
        if exists(f"{base}/{a.tag}_FINAL.jsonl"):
            status = "final"; break
        if exists(f"{base}/{a.tag}_TERMINATED.txt"):
            reason = s3("cp", f"{base}/{a.tag}_TERMINATED.txt", "-").stdout.strip()
            status = f"terminated_without_final:{reason}"; break
        hb = s3("cp", f"{base}/heartbeat_{a.tag}.txt", "-").stdout.strip()
        n = None
        for tok in hb.split():
            if tok.startswith("nights_done="):
                n = int(tok.split("=")[1])
        if n is not None and n != last_n:
            last_n = n; last_move = time.time()
        el = (time.time() - t0) / 60
        stall = (time.time() - last_move) / 60 > a.stall_min
        print(f"[{el:6.1f} min] {hb or 'no heartbeat yet'}{'   STALL' if stall else ''}", flush=True)
        if el > a.max_wait_min:
            status = "max_wait_exceeded"; break
        time.sleep(a.poll_sec)
    print("status:", status, flush=True)
    os.makedirs(P.X1_DATA, exist_ok=True)
    fetched = []
    for key, dst in ((f"{a.tag}_FINAL.jsonl", f"{a.tag}_FINAL.jsonl"), (f"{a.tag}_LAST.jsonl", f"{a.tag}_LAST.jsonl"),
                     (f"{a.tag}_partial.jsonl", f"{a.tag}_partial.jsonl"), (f"{a.tag}_run.log", f"{a.tag}_run.log"),
                     (f"{a.tag}_worker.log", f"{a.tag}_worker.log"), (f"{a.tag}_summary.json", f"{a.tag}_summary.json"),
                     (f"{a.tag}_versions.txt", f"{a.tag}_versions.txt")):
        if exists(f"{base}/{key}"):
            r = s3("cp", f"{base}/{key}", f"{P.X1_DATA}/{dst}", "--only-show-errors", timeout=1800)
            if r.returncode == 0:
                fetched.append(f"{P.X1_DATA}/{dst}")
    print("fetched:", fetched, flush=True)
    final = f"{P.X1_DATA}/{a.tag}_FINAL.jsonl"
    if os.path.isfile(final):
        n = sum(1 for _ in open(final)); ok = sum(1 for l in open(final) if '"status": "ok"' in l)
        print(f"FINAL: {n} records, {ok} ok", flush=True)
        write_sidecar(final, [f"{P.X1_EC2}/reex_code_v8.tar.gz"] if os.path.isfile(f"{P.X1_EC2}/reex_code_v8.tar.gz") else [], __file__,
                      f"ec2_{a.tag}", extra={"status": status, "s3_base": base}, index_path=f"{P.X1_LOGS}/PROVENANCE_INDEX.tsv")
    return 0 if status == "final" else 2


if __name__ == "__main__":
    sys.exit(main())
