#!/usr/bin/env python3
"""
Step 030 measure_families_v8: numbers/measure_families.csv from 197 rows to N_MEASURES (141).

Decision 3 (2026-09-12): drop the 60 per-site standardized copies (the rows flagged
is_siteZ_duplicate, every one a Brain, spectral power twin ending in _siteZ), add hypoxic burden
(Oxygenation), the periodic limb movement index (new family, Limb movements), sleep onset latency
and REM latency (Sleep timing and structure, beside total sleep time and efficiency).

    build_measure_families_v8.py                 rewrite numbers/measure_families.csv, backup
                                                 numbers/measure_families_PRE_V8.csv (refuses a
                                                 second run when the backup exists and differs)
    build_measure_families_v8.py --out PATH      write elsewhere (the smoke run), no backup
    build_measure_families_v8.py --dry-run       print the counts, write nothing

The four new rows carry rank_ranking_v3 and dC_ranking_v3 empty until step 105 ranks them;
nothing reads those two columns in the chain (they are the 2026-08-07 snapshot values kept for
the record). Every count is checked against cohort_spec, never against a literal.
Input and output paths come from cohort_spec only.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import argparse
import os
import shutil
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, paths.ANALYSIS_DIR)
import cohort_spec as CS  # noqa: E402

# definition text from RERUN_PLAN_V8.md (D5, D6, step 020) and the extraction code (extract_night.py
# emits sol_min and rem_latency_min). The family scheme columns copy the sibling rows of the
# family joined (spo2_pct_below_90 for Oxygenation, TST_min for Sleep timing and structure).
NEW_ROWS = [
    {"feature": "hypoxic_burden", "family": "Oxygenation",
     "definition": ("Hypoxic burden, percent-minutes of desaturation per hour of sleep: per respiratory "
                    "event the area between the pre-event baseline saturation and the trace, native "
                    "oximeter rate, event windows frozen per site in extract_v8/hb_definition.json"),
     "sibling": "spo2_pct_below_90"},
    {"feature": "plm_index", "family": "Limb movements",
     "definition": ("Periodic limb movement index, series of 4 or more limb movements 5 to 90 s apart, "
                    "0.5 to 10 s long (AASM rule) per hour of sleep; scoring convention differs by "
                    "hospital, rank-normalised within site"),
     "sibling": None},
    {"feature": "sol_min", "family": "Sleep timing and structure",
     "definition": "Sleep onset latency, minutes from lights off (or recording start) to the first sleep epoch",
     "sibling": "TST_min"},
    {"feature": "rem_latency_min", "family": "Sleep timing and structure",
     "definition": "REM latency, minutes from sleep onset to the first REM epoch",
     "sibling": "TST_min"},
]
LIMB_SCHEME = {"family": "Limb movements", "family_make_v2": "Limb movements",
               "family_audit": "Limb movements", "family_etable": "Limb movements",
               "keyword_family_efig05": "Limb movements", "keyword_family_make_v2": "Limb movements",
               "keyword_family_make_supp_v2": "Limb movements"}


def build(src):
    t = pd.read_csv(src)
    assert t.feature.is_unique, "duplicate feature in the source file"
    dup = t.is_siteZ_duplicate.astype(bool)
    dropped = t[dup].feature.tolist()
    bad = [f for f in dropped if not f.endswith(CS.RANKING_DROP_SUFFIX)]
    assert not bad, f"is_siteZ_duplicate rows that do not end in {CS.RANKING_DROP_SUFFIX}: {bad[:5]}"
    also = [f for f in t[~dup].feature if f.endswith(CS.RANKING_DROP_SUFFIX)]
    assert not also, f"rows ending in {CS.RANKING_DROP_SUFFIX} not flagged as duplicates: {also[:5]}"
    keep = t[~dup].copy()
    rows = []
    for nr in NEW_ROWS:
        assert nr["feature"] not in set(t.feature), f"{nr['feature']} already in the file"
        if nr["sibling"]:
            base = t[t.feature == nr["sibling"]].iloc[0].to_dict()
            assert base["family"] == nr["family"], (nr, base["family"])
        else:
            base = {c: np.nan for c in t.columns}
            base.update(LIMB_SCHEME)
        base = dict(base)
        base.update({"feature": nr["feature"], "family": nr["family"], "definition": nr["definition"],
                     "resolved_by": "explicit (v8 decision 3, 2026-09-12)", "is_siteZ_duplicate": False,
                     "rank_ranking_v3": np.nan, "dC_ranking_v3": np.nan,
                     "keyword_matches": "explicit:v8"})
        rows.append(base)
    new = pd.DataFrame(rows, columns=t.columns)
    out = pd.concat([keep, new], ignore_index=True)
    out["rank_ranking_v3"] = out["rank_ranking_v3"].astype("Int64")
    assert set(CS.NEW_MEASURES_V8) == {r["feature"] for r in NEW_ROWS}, "cohort_spec.NEW_MEASURES_V8 disagrees with this script"
    assert len(out) == CS.N_MEASURES, f"{len(out)} rows, cohort_spec.N_MEASURES is {CS.N_MEASURES}"
    assert out.feature.is_unique
    assert not out.is_siteZ_duplicate.astype(bool).any()
    assert out.family.notna().all() and out.definition.notna().all()
    return t, out, dropped


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=None, help="output path (default: rewrite the live file with a backup)")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    src = CS.MEASURE_FAMILIES_CSV
    if not os.path.exists(src):
        sys.exit(f"missing {src}")
    t, out, dropped = build(src)
    fam_old = t.family.value_counts().to_dict()
    fam_new = out.family.value_counts().to_dict()
    print(f"source {src}: {len(t)} rows, {len(dropped)} siteZ duplicates dropped, {len(NEW_ROWS)} added -> {len(out)} rows")
    print("families before:", fam_old)
    print("families after: ", fam_new)
    if a.dry_run:
        return 0
    dst = a.out or src
    if dst == src:
        bak = os.path.join(os.path.dirname(src), "measure_families_PRE_V8.csv")
        cur = open(src, "rb").read()
        if os.path.exists(bak) and open(bak, "rb").read() != cur:
            sys.exit(f"REFUSED: backup {bak} exists and differs from the current file (second run?)")
        if not os.path.exists(bak):
            shutil.copy2(src, bak)
        print(f"backup {bak}")
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    out.to_csv(dst, index=False)
    side = CS.sidecar(dst, __file__, extra_inputs=[src] if dst != src else [],
                      note=f"step 030: {len(t)} -> {len(out)} rows; dropped {len(dropped)} siteZ twins; added {[r['feature'] for r in NEW_ROWS]}",
                      record={"families_before": fam_old, "families_after": fam_new, "dropped": dropped})
    print(f"wrote {dst}  ({len(out)} rows)  sidecar {side}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
