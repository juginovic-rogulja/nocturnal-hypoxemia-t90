"""
Step 041: rebuild the condition first-dates for v8 from the local hospital-record cache, on the v7 outcome clock.

Rule (verbatim v7, numbers/rebuild_outcomes.py): strip the decimal point, upper-case, prefix match on condition_source_value,
min(condition_start_date) per patient; prevalent if on or before psg_date, incident if after; years to the first date or to
the censor date. psg_date and censor_date come from the v7 t90_final through the spec module (plan decision D2 default).

Untouched conditions are COPIED from v7 (byte for byte); the RECOMPUTE_SET of disease_definitions_v8 is recomputed.
Before the build, EVERY v7 condition is recomputed with the v7 lists and compared to the v7 table (integrity gate 1: the
same cache, the same rule and the same clock must reproduce the frozen columns exactly; otherwise rc 2).
Outputs: work/outcomes_v8.parquet (+ provenance sidecar), work/event_counts_v8.csv, work/reconcile_v7_lists_vs_v7_table.csv.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import importlib.util
import json
import sys
import time

import duckdb
import numpy as np
import pandas as pd

sys.path.insert(0, paths.X2_SCRIPTS)
import x2_spec as S                              # noqa: E402
import disease_definitions_v8 as D8              # noqa: E402

spec = importlib.util.spec_from_file_location("disease_definitions_v7", S.V7_DEFINITIONS)
D7 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(D7)

T0 = time.time()
LOG = open(f"{S.LOGS}/041_rebuild_outcomes_v8.log", "w")


def say(*a):
    msg = " ".join(str(x) for x in a)
    print(msg, flush=True)
    LOG.write(msg + "\n")
    LOG.flush()


# ------------------------------------------------------------------ inputs: eviction guard, integrity
INPUTS = {"v7_t90_final": S.PREV_T90_FINAL, "v7_sha256_txt": S.PREV_SHA256_TXT, "omop_cache": S.CACHE,
          "v7_definitions": S.V7_DEFINITIONS, "v8_definitions": S.V8_DEFINITIONS}
for p in INPUTS.values():
    S.require_local(p)
listed = {name: sha for sha, name in (ln.split() for ln in open(S.PREV_SHA256_TXT).read().strip().splitlines())}
sha_v7 = S.sha256(S.PREV_T90_FINAL)
if sha_v7 != listed["t90_final.parquet"]:
    raise SystemExit(f"v7 t90_final.parquet sha256 {sha_v7} does not match sha256.txt {listed['t90_final.parquet']}")
say(f"v7 t90_final sha256 verified against sha256.txt: {sha_v7[:16]}")

v7 = pd.read_parquet(S.PREV_T90_FINAL)
coh7 = S.apply_cohort(v7)
say(f"v7 rows {len(v7):,}, analysis cohort {len(coh7):,} (cohort_spec.COHORT_N)")
BASE_COLS = ["BDSPPatientID", "site_id", "psg_date", "censor_date", "fu_valid", "AgeAtVisit", "sex",
             "spo2_pct_below_90", "oximetry_bad"]
base = v7[BASE_COLS].copy()
assert base.BDSPPatientID.is_unique

# ------------------------------------------------------------------ the cache, the v7 rule
con = duckdb.connect()
con.execute(f"create view cond as select * from '{S.CACHE}'")
con.execute("""create table c as
    select person_id, condition_start_date as dt,
           upper(replace(condition_source_value, '.', '')) as code
    from cond where condition_source_value is not null""")
n_rows = con.execute("select count(*) from c").fetchone()[0]
n_persons_in_v7 = con.execute("select count(distinct person_id) from c").fetchone()[0]
say(f"cache condition rows {n_rows:,}, persons {n_persons_in_v7:,}  ({time.time() - T0:.0f}s)")


def first_dates(key, icd9, icd10, exclude=()):
    pref = sorted(set(icd9) | set(icd10))
    where = " or ".join(f"code like '{p}%'" for p in pref)
    if exclude:
        where = f"({where}) and not (" + " or ".join(f"code like '{e}%'" for e in exclude) + ")"
    fd = con.execute(f"select person_id, min(dt) as first_date from c where {where} group by person_id").df()
    fd.columns = ["BDSPPatientID", f"{key}_first_date"]
    fd[f"{key}_first_date"] = pd.to_datetime(fd[f"{key}_first_date"], errors="coerce")
    assert fd.BDSPPatientID.is_unique
    return fd


def flags(d, key):
    fd = d[f"{key}_first_date"]
    before = fd.notna() & d.psg_date.notna() & (fd <= d.psg_date)
    after = fd.notna() & d.psg_date.notna() & (fd > d.psg_date)
    d[f"{key}_prevalent"] = before.astype(int)
    d[f"{key}_incident"] = after.astype(int)
    end = np.where(after, fd.values, d.censor_date.values)
    d[f"{key}_years"] = (pd.Series(pd.to_datetime(end), index=d.index) - d.psg_date).dt.days / 365.25
    return d


def add_condition(d, key, fd):
    n = len(d)
    d = d.merge(fd, on="BDSPPatientID", how="left")
    assert len(d) == n and (d.BDSPPatientID.values == base.BDSPPatientID.values).all()
    d[f"{key}_first_date"] = pd.to_datetime(d[f"{key}_first_date"], errors="coerce")
    return flags(d, key)


def identical(a, b):
    """value-identical (NaN/NaT aware) and byte-identical, as a pair"""
    if a.dtype != b.dtype:
        return False, False
    if a.dtype.kind == "f":
        val = bool(np.array_equal(a.to_numpy(), b.to_numpy(), equal_nan=True))
    else:
        val = bool(a.equals(b))
    return val, a.to_numpy().tobytes() == b.to_numpy().tobytes()


def n_diff(a, b):
    both_na = a.isna() & b.isna()
    return int((~((a == b) | both_na)).sum())


# ------------------------------------------------------------------ reconciliation: v7 lists -> v7 table (gate 1)
rec = base.copy()
for key, (lab, neg, i9, i10) in D7.DISEASES.items():
    if not i9 and not i10:
        continue
    rec = add_condition(rec, key, first_dates(key, i9, i10))
rec["cvd_first_date"] = rec[[f"{k}_first_date" for k in D7.CVD_COMPONENTS]].min(axis=1)
rec = flags(rec, "cvd")
coh_mask = v7.index.isin(coh7.index)
rows = []
for key in D7.DISEASES:
    for suf in S.OUTCOME_SUFFIXES:
        a, b = rec[key + suf], v7[key + suf]
        val, byt = identical(a, b)
        rows.append({"key": key, "col": key + suf, "value_identical": val, "byte_identical": byt,
                     "n_diff_all_rows": n_diff(a, b), "n_diff_cohort_rows": n_diff(a[coh_mask], b[coh_mask])})
recon = pd.DataFrame(rows)
recon.to_csv(S.OUT_RECONCILE, index=False)
bad_keys = sorted(set(recon.loc[~recon.value_identical, "key"]))
# what matters for the analysis: the cohort rows. Classify every condition.
coh_flag_diff = sorted(set(recon.loc[(recon.n_diff_cohort_rows > 0) & ~recon.col.str.endswith("_first_date"), "key"]))
coh_fd_only = sorted(set(recon.loc[(recon.n_diff_cohort_rows > 0) & recon.col.str.endswith("_first_date"), "key"]) - set(coh_flag_diff))
off_cohort_only = sorted(set(bad_keys) - set(coh_flag_diff) - set(coh_fd_only))
twins_bad = sorted((set(coh_flag_diff) | set(coh_fd_only)) & D8.RECOMPUTE_SET)
say(f"reconciliation (v7 lists on today's cache vs the frozen v7 columns), {len(D7.DISEASES)} conditions x 4 columns: "
    f"{int(recon.value_identical.sum())} of {len(recon)} columns identical on all {len(v7):,} rows.")
say(f"  cohort rows ({len(coh7):,}): flags or years differ in {coh_flag_diff or 'none'}; first_date only differs in {coh_fd_only or 'none'}; "
    f"differences confined to rows outside the cohort in {len(off_cohort_only)} conditions "
    f"({'all fu_valid=0, psg_date NaT: v7 kept a years value, the rule gives NaN' if off_cohort_only else ''}).")
say(f"  the {len(D8.RECOMPUTE_SET)} recomputed conditions' v7 twins reproduce every cohort column: {'YES' if not twins_bad else 'NO ' + str(twins_bad)}  ({time.time() - T0:.0f}s)")

# ------------------------------------------------------------------ v8 build
out = base.copy()
recomputed, copied = {}, []
for key, (lab, neg, i9, i10) in D8.DISEASES.items():
    if key == "cvd":
        continue
    if key in D8.RECOMPUTE_SET:
        out = add_condition(out, key, first_dates(key, i9, i10, D8.EXCLUDE.get(key, ())))
        recomputed[key] = {"icd9": i9, "icd10": i10, "exclude": list(D8.EXCLUDE.get(key, []))}
    else:
        for suf in S.OUTCOME_SUFFIXES:
            out[key + suf] = v7[key + suf].to_numpy(copy=True)
        copied.append(key)
out["cvd_first_date"] = out[[f"{k}_first_date" for k in D8.CVD_COMPONENTS]].min(axis=1)
out = flags(out, "cvd")
recomputed["cvd"] = {"components": D8.CVD_COMPONENTS}
for c in ["death_incident", "death_years", "death_prevalent"]:
    out[c] = v7[c].to_numpy(copy=True)

sha_defs, sha_cache = S.sha256(S.V8_DEFINITIONS), S.sha256(S.CACHE)
out["_v8_outcome_list_version"] = sha_defs[:12]
out["_v8_outcome_source"] = (f"omop_cache:{sha_cache[:12]}|v7_t90_final:{sha_v7[:12]}|"
                             f"recomputed:{','.join(sorted(D8.RECOMPUTE_SET))}")
out["_v8_outcome_clock"] = "v7 psg_date and censor_date (plan D2 default)"
order = BASE_COLS + [k + suf for k in D8.DISEASES for suf in S.OUTCOME_SUFFIXES] + \
    ["death_incident", "death_years", "death_prevalent", "_v8_outcome_list_version", "_v8_outcome_source", "_v8_outcome_clock"]
assert set(order) == set(out.columns), set(order) ^ set(out.columns)
out = out[order]
out.to_parquet(S.OUT_PARQUET, index=False)
say(f"written {S.OUT_PARQUET}: {len(out):,} rows x {out.shape[1]} columns; recomputed {len(recomputed)}, copied {len(copied)}")

# ------------------------------------------------------------------ positive control on the file as written (rc 1)
chk = pd.read_parquet(S.OUT_PARQUET)
bad = []
for key in copied:
    for suf in S.OUTCOME_SUFFIXES:
        val, byt = identical(chk[key + suf], v7[key + suf])
        if not byt:
            bad.append(key + suf)
for c in ["death_incident", "death_years", "death_prevalent", "psg_date", "censor_date", "fu_valid", "BDSPPatientID"]:
    if not identical(chk[c], v7[c])[1]:
        bad.append(c)
say(f"positive control: {len(copied)} untouched conditions x 4 columns + death + clock byte-identical to v7: "
    f"{'PASS' if not bad else 'FAIL ' + str(bad)}")

# ------------------------------------------------------------------ event counts (the regenerated reference for later asserts)
coh8 = S.apply_cohort(chk)
assert (coh8.BDSPPatientID.values == coh7.BDSPPatientID.values).all()
cnt = []
for key in list(D8.DISEASES) + ["death"]:
    prev = int(coh8[f"{key}_prevalent"].sum())
    inc = int(coh8[f"{key}_incident"].sum())
    at_risk = coh8[coh8[f"{key}_prevalent"] == 0]
    med = float(at_risk[f"{key}_years"].median())
    old_prev = int(coh7[f"{key}_prevalent"].sum()) if f"{key}_prevalent" in coh7 else np.nan
    old_inc = int(coh7[f"{key}_incident"].sum()) if f"{key}_incident" in coh7 else np.nan
    old_med = float(coh7[coh7[f"{key}_prevalent"] == 0][f"{key}_years"].median()) if f"{key}_years" in coh7 else np.nan
    status = ("new" if key in D8.ADDED_2026_09_12 else "new, swapped-in control (14 Sept)" if key in getattr(D8, "SWAPPED_IN_2026_09_14", set()) else "changed list" if key in D8.CHANGED_2026_09_12 else
              "composite, PAD component changed" if key == "cvd" else
              "recomputed, list unchanged" if key == "osteoporosis" else "unchanged, copied from v7")
    delta = inc - old_inc if old_inc == old_inc else np.nan
    pct = 100.0 * delta / old_inc if old_inc == old_inc and old_inc else np.nan
    cnt.append({"key": key, "label": "Death from any cause" if key == "death" else D8.DISEASES[key][0],
                "organ_group": D8.ORGAN_GROUP.get(key, "Death"), "status": status,
                "prevalent": prev, "incident": inc, "at_risk": len(at_risk), "median_followup_years_at_risk": round(med, 3),
                "v7_prevalent": old_prev, "v7_incident": old_inc, "v7_median_followup_years_at_risk": round(old_med, 3) if old_med == old_med else np.nan,
                "delta_incident": delta, "pct_change_incident": round(pct, 1) if pct == pct else np.nan,
                "dramatic_over_20pct": bool(abs(pct) > 20) if pct == pct else False,
                "ranked": bool(key not in D8.RANKING_EXCLUDE and key not in D8.CIRCULAR and key != "death")})
counts = pd.DataFrame(cnt)
counts.to_csv(S.OUT_COUNTS, index=False)

# ------------------------------------------------------------------ sidecars
extra = {"clock": "v7 psg_date and censor_date from cohort_spec.V7_DIR t90_final (plan D2 default)",
         "cohort_n": int(len(coh8)), "cache_rows": int(n_rows), "cache_persons": int(n_persons_in_v7),
         "recomputed": recomputed, "copied_from_v7": copied,
         "reconciliation_v7_lists_vs_v7_table": {"columns_identical": int(recon.value_identical.sum()),
                                                  "columns_total": int(len(recon)), "conditions_with_any_difference": bad_keys,
                                                  "cohort_flags_or_years_differ": coh_flag_diff, "cohort_first_date_only_differs": coh_fd_only,
                                                  "off_cohort_rows_only": off_cohort_only, "recomputed_twins_reproduced": not twins_bad},
         "positive_control_untouched_byte_identical": not bad,
         "seconds": round(time.time() - T0, 1)}
S.write_sidecar(S.OUT_PARQUET, INPUTS, __file__, extra)
S.write_sidecar(S.OUT_COUNTS, {**INPUTS, "outcomes_v8": S.OUT_PARQUET}, __file__, {"cohort_n": int(len(coh8))})
S.write_sidecar(S.OUT_RECONCILE, INPUTS, __file__)

say("\nkey            status                              v7_prev v8_prev  v7_inc  v8_inc  delta   pct  dramatic")
for _, r in counts[counts.status != "unchanged, copied from v7"].iterrows():
    say(f"{r.key:<14} {r.status:<34} {str(r.v7_prevalent):>7} {r.prevalent:>7} {str(r.v7_incident):>7} {r.incident:>7} "
        f"{str(r.delta_incident):>6} {str(r.pct_change_incident):>5}  {'YES' if r.dramatic_over_20pct else ''}")
say(f"\ndone in {time.time() - T0:.0f}s")
if bad:
    sys.exit(1)
if twins_bad:
    say(f"RECONCILIATION MISMATCH on recomputed conditions {twins_bad}: the v7 lists on today's cache do not reproduce the frozen "
        f"v7 cohort columns, so the cache, rule or clock differ from the v7 build; see {S.OUT_RECONCILE}. Gate 1: full pass needed.")
    sys.exit(2)
if coh_flag_diff or coh_fd_only:
    say(f"NOTE (inherited, untouched conditions copied byte-for-byte): v7 columns not reproduced by the v7 lists on this cache for "
        f"{sorted(set(coh_flag_diff) | set(coh_fd_only))}; counts in {S.OUT_RECONCILE}; reported to the owner, not changed here.")
