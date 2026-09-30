"""
LADDER_SUMMARY.csv: one row per threshold, the columns Alen asked for, built ONLY
from artifacts already on disk (no number is retyped by hand).

  threshold      T80 T85 T88 T90 T92 T95, plus the two re-derived secondary rungs
                 for 88/90 (pure same-code ladder) and the A95 banded framing row.
  source         which column the row is computed from. T88 and T90 are the FROZEN
                 columns per the task's labels. T90's dC cell is the published value
                 (the refit equals it to 1e-9, which is the positive control).
  median_pct, pct_gt10, pct_le1     distribution on the cohort.
  dC_mean                           the paper's discrimination gain, mean over 48.
  n_fdr_sig_of_48_persd             per-SD Cox, BH-FDR q<0.05 among the 48 real.
  median_persd_hr                   median per-SD HR over all 48.
  median_persd_hr_sig               median per-SD HR among the FDR-significant.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import json
import os

import numpy as np
import pandas as pd
import sys as _sys_v8; _sys_v8.path.insert(0, paths.ANALYSIS_DIR)  # v8 sweep 2026-09-12
from cohort_spec import N_RANKED_OUTCOMES  # v8 sweep 2026-09-12
__PUB_DC_V7__ = float(__import__("pandas").read_csv(f"{paths.NUMBERS_DIR}/ranking_v3.csv", comment="#").set_index("feature").loc["spo2_pct_below_90", "dC"])   # v7: T90 gain from the CURRENT ranking, was the August literal 0.02033875497249791


HERE = os.path.dirname(os.path.abspath(__file__))
PUB_DC = __PUB_DC_V7__

ps = json.load(open(os.path.join(HERE, "persd_summary.json")))
dc = json.load(open(os.path.join(HERE, "dC_result.json")))
bd = json.load(open(os.path.join(HERE, "ladder_banded_results.json")))
bpc = json.load(open(os.path.join(HERE, "banded_positive_control.json")))
L = pd.read_csv(os.path.join(HERE, "ladder_banded_long.csv"))

DIST, STATS, DC = ps["distributions"], ps["ladder_stats"], dc["dC"]

ROWS = [
    ("T80", "t80", "re-derived 2026-08-21", DC["t80"]),
    ("T85", "t85", "re-derived 2026-08-21 (= T85 run bit-for-bit)", DC["t85"]),
    ("T88", "t88_frozen", "frozen spo2_pct_below_88", DC["t88_frozen"]),
    ("T90", "t90_frozen", "frozen spo2_pct_below_90, published", PUB_DC),
    ("T92", "t92", "re-derived 2026-08-21", DC["t92"]),
    ("T95", "t95", "re-derived 2026-08-21", DC["t95"]),
    ("T88_rederiv", "t88_rederiv", "secondary: same-code ladder rung", DC["t88_rederiv"]),
    ("T90_rederiv", "t90_rederiv", "secondary: same-code ladder rung", DC["t90_rederiv"]),
]

rows = []
for name, col, src, dcv in ROWS:
    d_, s_ = DIST[col], STATS[col]
    rows.append({
        "threshold": name, "source": src, "n": d_["n"],
        "median_pct": round(d_["median"], 4),
        "pct_gt10": round(d_["pct_gt10"], 2),
        "pct_le1": round(d_["pct_le1"], 2),
        "dC_mean": round(dcv, 9),
        "n_fdr_sig_of_48_persd": s_["n_sig"],
        "median_persd_hr": round(s_["median_hr"], 3),
        "median_persd_hr_sig": (round(s_["median_sig_hr"], 3)
                                if s_["median_sig_hr"] is not None else np.nan),
    })

# the A95 preserved-saturation banded framing, one row. Its continuous dC equals
# T95's by the mirror identity (checked in b03), stated rather than re-listed.
# Banded significance is counted over ALL real graded outcomes (the T85 report's
# convention), denominator in the source string.
a95d = DIST["a95"]
alabs = list(bd["a95_preserved"]["band_n"].keys())
worst_lab = alabs[-1]
sub = L[(L.exposure == "a95_preserved") & (L.band == worst_lab)
        & (~L.negative_control)]
sig = sub[sub.q < 0.05]
rows.append({
    "threshold": "A95_preserved_bands",
    "source": (f"percent of night at/above 95, bands {' / '.join(alabs)}, "
               f"reference most-preserved. BANDED counts, worst band vs reference, "
               f"over the {len(sub)} real graded outcomes, not per-SD"),
    "n": bd["a95_preserved"]["n"],
    "median_pct": round(a95d["median"], 4),
    "pct_gt10": np.nan, "pct_le1": np.nan,
    "dC_mean": round(DC["a95"], 9),
    "n_fdr_sig_of_48_persd": int((sig).shape[0]),
    "median_persd_hr": round(float(sub.hr.median()), 3) if len(sub) else np.nan,
    "median_persd_hr_sig": round(float(sig.hr.median()), 3) if len(sig) else np.nan,
})

S = pd.DataFrame(rows)
S.to_csv(os.path.join(HERE, "LADDER_SUMMARY.csv"), index=False)
print(S.to_string(index=False))

# where does the information peak: judged on dC (the paper's ranking currency) and
# on the per-SD FDR count, over the six primary rungs.
prim = S[S.threshold.isin(["T80", "T85", "T88", "T90", "T92", "T95"])]
peak_dc = prim.loc[prim.dC_mean.idxmax()]
peak_ns = prim.loc[prim.n_fdr_sig_of_48_persd.idxmax()]
print(f"\nPEAK by dC:            {peak_dc.threshold}  (dC {peak_dc.dC_mean:.6f})")
print(f"PEAK by n significant: {peak_ns.threshold}  ({int(peak_ns.n_fdr_sig_of_48_persd)}/{N_RANKED_OUTCOMES})")
json.dump({"peak_by_dC": peak_dc.threshold, "peak_dC": float(peak_dc.dC_mean),
           "peak_by_nsig": peak_ns.threshold,
           "peak_nsig": int(peak_ns.n_fdr_sig_of_48_persd)},
          open(os.path.join(HERE, "ladder_peak.json"), "w"), indent=2)
print("\nwrote LADDER_SUMMARY.csv, ladder_peak.json")
