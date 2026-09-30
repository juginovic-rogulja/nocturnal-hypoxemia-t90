"""Provenance sidecars: <output>.provenance.json with sha256 of the output, every input, the script, time and
environment, plus one line in X1/logs/PROVENANCE_INDEX.tsv. A printed number without a sidecar is not printable."""
from __future__ import annotations
import hashlib, json, os, platform, socket, sys, time
from pathlib import Path


def sha256_file(path, chunk=1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def _versions():
    out = {"python": sys.version.split()[0], "platform": platform.platform()}
    for m in ("numpy", "scipy", "pandas", "h5py", "s3fs", "pyarrow"):
        try:
            out[m] = __import__(m).__version__
        except Exception:
            out[m] = None
    return out


def write_sidecar(output, inputs, script, step, extra=None, index_path=None):
    output = str(output)
    rec = {"output": output, "output_sha256": sha256_file(output) if os.path.isfile(output) else None,
           "output_bytes": os.path.getsize(output) if os.path.isfile(output) else None,
           "step": step, "script": str(script), "script_sha256": sha256_file(script) if os.path.isfile(script) else None,
           "inputs": {str(p): (sha256_file(p) if os.path.isfile(p) else "DIR_OR_MISSING") for p in inputs},
           "time_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "host": socket.gethostname(),
           "env": _versions(), "extra": extra or {}}
    side = output + ".provenance.json"
    Path(side).write_text(json.dumps(rec, indent=1, default=str))
    if index_path:
        with open(index_path, "a") as f:
            f.write(f"{rec['time_utc']}\t{step}\t{output}\t{rec['output_sha256']}\t{os.path.basename(str(script))}\n")
    return side
