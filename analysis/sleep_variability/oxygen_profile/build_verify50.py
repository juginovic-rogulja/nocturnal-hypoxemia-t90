"""
Draw the 50 patients for the independent verification pass.

Fixed seed, drawn uniformly at random from the full frozen 19,173-adult cohort (no
stratification, no exclusion of the awkward records, because the point of the pass is to catch
what a stratified sample would hide). The seed and the drawn IDs are written next to the task
file so the draw is reproducible.

    python3 build_verify50.py
"""
import json
import os

import pandas as pd
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SEED = 20260821
N_DRAW = 50

tasks = pd.read_csv(os.path.join(HERE, "tasks_all.csv"))
assert len(tasks) == 19173, len(tasks)

rng = np.random.RandomState(SEED)
idx = rng.choice(len(tasks), size=N_DRAW, replace=False)
draw = tasks.iloc[np.sort(idx)].reset_index(drop=True)
assert len(draw) == N_DRAW and draw.BDSPPatientID.nunique() == N_DRAW

draw.to_csv(os.path.join(HERE, "tasks_verify50.csv"), index=False)
json.dump({"seed": SEED, "n": N_DRAW, "drawn_from": "tasks_all.csv, 19,173 rows",
           "BDSPPatientID": [int(v) for v in draw.BDSPPatientID],
           "sites": {str(k): int(v) for k, v in draw.site.value_counts().items()}},
          open(os.path.join(HERE, "verify50_draw.json"), "w"), indent=2)
print(f"seed {SEED}   drew {len(draw)} patients   sites {draw.site.value_counts().to_dict()}")
print("wrote tasks_verify50.csv, verify50_draw.json")
