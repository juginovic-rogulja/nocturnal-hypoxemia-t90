"""
Cohort aggregates for the overnight oxygen profile sheet.

Reads the extraction's per-patient file, applies the paper's four T90 bands from the FROZEN
spo2_pct_below_90 column (not from anything this run computed, so the banding is the paper's),
and writes the medians and interquartile ranges the figure draws. The figure script recomputes
every one of these from the same per-patient file and aborts if any disagrees, so this file is
the contract between the extraction and the sheet rather than a second source of truth.

Inclusion for the decile panel: a record must carry all ten decile means. A span that is
entirely oximeter dropout has no mean and the worker reports none rather than fabricating one,
so a record missing a span would otherwise contribute to some deciles and not others and bend
the curve. Inclusion for the stage panel is per stage: a record contributes to a stage if it
has at least one valid oximetry sample scored in that stage.

    python3 b01_profile.py
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import json
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DEC = [f"{i + 1:02d}" for i in range(10)]
STAGES = (("wake", "Wake"), ("n1", "N1"), ("n2", "N2"), ("n3", "N3"), ("rem", "REM"))
BANDS = [("le1", "T90 1% or less", lambda t: t <= 1),
         ("1to5", "T90 above 1% to 5%", lambda t: (t > 1) & (t <= 5)),
         ("5to10", "T90 above 5% to 10%", lambda t: (t > 5) & (t <= 10)),
         ("gt10", "T90 above 10%", lambda t: t > 10)]
# Reference band counts are NOT a literal (integrity gate 6): they are recomputed here from the v7
# analysis table through the same cohort rule, so a stale frozen_reference.csv fails loudly.
import sys  # noqa: E402
T90ROOT = paths.T90_ROOT
V7_TABLE = f"{paths.TABLES_DIR}/t90_final.parquet"
sys.path.insert(0, paths.ANALYSIS_DIR)
from cohort_spec import apply_cohort, COHORT_N  # noqa: E402
_ref = apply_cohort(pd.read_parquet(V7_TABLE))
assert len(_ref) == COHORT_N, len(_ref)
BAND_N_FROZEN = {k: int(rule(_ref.spo2_pct_below_90).sum()) for k, _l, rule in BANDS}
assert sum(BAND_N_FROZEN.values()) == COHORT_N, BAND_N_FROZEN
# Rule 6 of RECOMPOSE_R30_RULES.md: the stage panel excludes the nights whose current BDSP file
# carries a collapsed staging (the extraction of 2026-08-21 read those files after the collapse).
d = pd.read_csv(os.path.join(HERE, "oxyprofile_per_patient.csv"), low_memory=False)
assert len(d) == COHORT_N, len(d)
# v8.2 (2026-09-15): the collapsed-staging nights are the ones the re-extraction itself flagged (_v8_stage_collapsed, rule D11 applied
# at assembly from the per-state minutes of every record), not a list carried over from the 2026-09-08 audit; that list is printed beside it.
assert "_v8_stage_collapsed" in d.columns, "oxyprofile_per_patient.csv carries no _v8_stage_collapsed column: assemble the masked fleet_v8_2 shards first"
COLLAPSED_IDS = set(d.loc[d._v8_stage_collapsed.fillna(0).astype(int) == 1, "BDSPPatientID"].astype(int))
COLLAPSED_LIST = "_v8_stage_collapsed column of oxyprofile_per_patient.csv (rule D11 at assembly, X4_groupP/mask_collapsed.py)"
_audit_ids = set(pd.read_csv(f"{paths.T90_ROOT}/EncodingB_Audit_2026-09-06/s3_dates/degenerate_staging_nights.csv").BDSPPatientID.astype(int))
print(f"collapsed-staging nights flagged by the extraction: {len(COLLAPSED_IDS):,}; on the 2026-09-08 audit list: {len(_audit_ids):,}; overlap {len(COLLAPSED_IDS & _audit_ids):,}")
ok = d[d.status == "ok"].copy()
print(f"records {len(d):,}   ok {len(ok):,}")

# the paper's bands, from the frozen exposure column
ok["band"] = None
for key, _lab, rule in BANDS:
    ok.loc[rule(ok.spo2_pct_below_90), "band"] = key
assert ok.band.notna().all(), int(ok.band.isna().sum())
band_n = {k: int((ok.band == k).sum()) for k, _l, _r in BANDS}
print("band sizes:", band_n)
assert band_n == BAND_N_FROZEN, (band_n, BAND_N_FROZEN)

full = ok[ok.n_deciles_with_data == 10].copy()
# stage panel frame: ok records minus the collapsed-staging nights (rule 6)
stg = ok[~ok.BDSPPatientID.astype(int).isin(COLLAPSED_IDS)].copy()
n_collapsed_in_ok = int(ok.BDSPPatientID.astype(int).isin(COLLAPSED_IDS).sum())
print(f"stage panel: {len(stg):,} records after excluding {n_collapsed_in_ok:,} collapsed-staging nights")
print(f"records with a complete ten-decile profile: {len(full):,} "
      f"({100 * len(full) / len(ok):.1f}% of ok)")


def quartiles(s):
    s = pd.Series(s).dropna()
    return {"n": int(len(s)), "median": float(s.median()),
            "q1": float(s.quantile(0.25)), "q3": float(s.quantile(0.75))}


def decile_rows(frame):
    out = []
    for i, t in enumerate(DEC, start=1):
        q = quartiles(frame[f"dec_mean_{t}"])
        q["decile"] = i
        out.append(q)
    return out


def stage_rows(frame):
    out = {}
    for key, lab in STAGES:
        q = quartiles(frame[f"{key}_mean"])
        q["minutes_median"] = float(frame[f"{key}_min"].dropna().median())
        q["minutes_q1"] = float(frame[f"{key}_min"].dropna().quantile(0.25))
        q["minutes_q3"] = float(frame[f"{key}_min"].dropna().quantile(0.75))
        out[lab] = q
    return out


summary = {
    "source": "oxyprofile_per_patient.csv, 2026-09-15 EC2 re-extraction (v8.2 group P, mode fleet_v8_2); bands from the v8 "
              "spo2_pct_below_90 (frozen_reference_fleet_v8_2.csv from data_frozen_v8_2026-09); stage panel excludes the "
              "collapsed-staging nights flagged by the extraction itself (rule D11)",
    "cohort_n": COHORT_N,
    "band_reference_table": V7_TABLE,
    "stage_panel_exclusion_list": COLLAPSED_LIST,
    "stage_panel_excluded_n": n_collapsed_in_ok,
    "stage_panel_n": int(len(stg)),
    "n_ok": int(len(ok)),
    "n_complete_decile_profile": int(len(full)),
    "band_labels": {k: lab for k, lab, _r in BANDS},
    "band_n_cohort": band_n,
    "band_n_complete_profile": {k: int((full.band == k).sum()) for k, _l, _r in BANDS},
    "deciles_overall": decile_rows(full),
    "deciles_by_band": {k: decile_rows(full[full.band == k]) for k, _l, _r in BANDS},
    "stages_overall": stage_rows(stg),
    "stages_by_band": {k: stage_rows(stg[stg.band == k]) for k, _l, _r in BANDS},
}

# headline quantities the legend may quote
ov = summary["deciles_overall"]
st = summary["stages_overall"]
summary["headline"] = {
    "decile_median_min": float(min(r["median"] for r in ov)),
    "decile_median_max": float(max(r["median"] for r in ov)),
    "decile_median_range_pp": float(max(r["median"] for r in ov)
                                    - min(r["median"] for r in ov)),
    "decile_of_lowest_median": int(min(ov, key=lambda r: r["median"])["decile"]),
    "stage_lowest": min(st, key=lambda k: st[k]["median"]),
    "stage_highest": max(st, key=lambda k: st[k]["median"]),
    "stage_median_spread_pp": float(max(v["median"] for v in st.values())
                                    - min(v["median"] for v in st.values())),
    "gt10_wake_minus_rem_pp": float(summary["stages_by_band"]["gt10"]["Wake"]["median"]
                                    - summary["stages_by_band"]["gt10"]["REM"]["median"]),
    "le1_wake_minus_rem_pp": float(summary["stages_by_band"]["le1"]["Wake"]["median"]
                                   - summary["stages_by_band"]["le1"]["REM"]["median"]),
}

json.dump(summary, open(os.path.join(HERE, "overnight_profile_summary.json"), "w"), indent=2)

print("\nDECILE PROFILE, whole cohort (median [IQR] mean SpO2, %)")
for r in ov:
    print(f"  decile {r['decile']:2d}  n {r['n']:,}  {r['median']:6.2f} "
          f"[{r['q1']:6.2f} to {r['q3']:6.2f}]")
print("\nDECILE PROFILE by T90 band, medians only")
for k, lab, _r in BANDS:
    rows = summary["deciles_by_band"][k]
    print(f"  {lab:<22} n {rows[0]['n']:5,}  " + " ".join(f"{r['median']:6.2f}" for r in rows))
print("\nMEAN SpO2 BY STAGE, whole cohort")
for _k, lab in STAGES:
    q = st[lab]
    print(f"  {lab:<5} n {q['n']:,}  {q['median']:6.2f} [{q['q1']:6.2f} to {q['q3']:6.2f}]   "
          f"minutes median {q['minutes_median']:6.1f}")
print("\nMEAN SpO2 BY STAGE and band, medians")
for k, lab, _r in BANDS:
    row = summary["stages_by_band"][k]
    print(f"  {lab:<22} " + "  ".join(f"{l} {row[l]['median']:6.2f}" for _kk, l in STAGES))
print("\nheadline:", json.dumps(summary["headline"], indent=1))
print("\nwrote overnight_profile_summary.json")
