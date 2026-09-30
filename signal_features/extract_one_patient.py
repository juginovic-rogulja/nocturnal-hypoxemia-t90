"""
Run the full feature pipeline on ONE patient and write a feature row.

Usage:
    python -m feature_pipeline.extract_one_patient <path-to-.h5>

Outputs:
    feature_pipeline/outputs/features_<patient_id>_<timestamp>.json
    feature_pipeline/outputs/probe_<patient_id>.json   (h5 diagnostic)
"""

from __future__ import annotations
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..")))  # repository root, where paths.py lives
import paths

import json
import sys
import time
from datetime import datetime
from pathlib import Path

from .core.h5_reader import HSPH5Reader, probe_h5
from .features.brain_features import extract_brain_features
from .features.heart_features import extract_heart_features
from .features.lung_features import extract_lung_features
from .features.microstructure import extract_microstructure


def run(path: str) -> dict:
    t0 = time.time()
    probe = probe_h5(path)

    with HSPH5Reader(path) as r:
        brain = extract_brain_features(r)
        micro = extract_microstructure(r)
        heart = extract_heart_features(r)
        lung = extract_lung_features(r)

    rec = {
        "source": str(path),
        "meta": {
            "n_channels": probe["n_channels"],
            "duration_hr": probe["duration_hr"],
            "fs": probe["fs"],
            "stage_histogram_epochs": probe["stage_histogram_epochs"],
        },
        "brain": brain,
        "microstructure": micro,
        "heart": heart,
        "lung": lung,
        "runtime_sec": round(time.time() - t0, 2),
        "extracted_at": datetime.now().isoformat(timespec="seconds"),
    }
    return rec


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else (
        f"{paths.PSG_DIR}/example/ses-1/eeg/"
        "sub-SXXXXNNNNNNNNNNN_ses-1.h5"
    )
    rec = run(path)
    outdir = Path(paths.FEATURE_OUTPUTS_DIR)
    outdir.mkdir(exist_ok=True, parents=True)
    stem = Path(path).stem
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    fp = outdir / f"features_{stem}_{ts}.json"
    fp.write_text(json.dumps(rec, indent=2, default=str))
    print(f"Wrote {fp}")
    print(f"Runtime: {rec['runtime_sec']} sec")
    print(f"Brain: {len(rec['brain'])}  Micro: {len(rec['microstructure'])}  "
          f"Heart: {len(rec['heart'])}  Lung: {len(rec['lung'])}")
