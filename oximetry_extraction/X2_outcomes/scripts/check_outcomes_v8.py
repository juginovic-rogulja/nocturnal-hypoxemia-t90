"""
Checks on work/outcomes_v8.parquet: invariants, the old-vs-new table per changed outcome, the second code path (pandas prefix
match on 300 patients, psv step 041 smoke), the quantified alternatives (gate 12) and the report logs/OUTCOMES_V8.md.

rc 1 when a hard invariant fails (untouched columns not byte-identical, prevalent and incident not exclusive, flag or years
arithmetic inconsistent, second code path disagrees). Breaches of the birth and censor date invariants are counted and
reported per condition; the ones inside untouched conditions are inherited from v7 and cannot be changed in this lane, the
ones inside recomputed conditions follow the same rule by design and are listed with the alternative and its effect.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import datetime
import importlib.util
import json
import re
import sys

import duckdb
import numpy as np
import pandas as pd

sys.path.insert(0, paths.X2_SCRIPTS)
import x2_spec as S                              # noqa: E402
import disease_definitions_v8 as D8              # noqa: E402

spec = importlib.util.spec_from_file_location("disease_definitions_v7", S.V7_DEFINITIONS)
D7 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(D7)
TODAY = pd.Timestamp(datetime.date.today())
LOG = open(f"{S.LOGS}/042_check_outcomes_v8.log", "w")


def say(*a):
    msg = " ".join(str(x) for x in a)
    print(msg, flush=True)
    LOG.write(msg + "\n")
    LOG.flush()


INPUTS = {"outcomes_v8": S.OUT_PARQUET, "event_counts_v8": S.OUT_COUNTS, "v7_t90_final": S.PREV_T90_FINAL, "omop_cache": S.CACHE,
          "psg_meta_I0002": S.PSG_META["I0002"], "psg_meta_I0006": S.PSG_META["I0006"], "v8_definitions": S.V8_DEFINITIONS}
for p in INPUTS.values():
    S.require_local(p)
v8 = pd.read_parquet(S.OUT_PARQUET)
v7 = pd.read_parquet(S.PREV_T90_FINAL)
coh8, coh7 = S.apply_cohort(v8), S.apply_cohort(v7)
assert (coh8.BDSPPatientID.values == coh7.BDSPPatientID.values).all()
KEYS = list(D8.DISEASES)
UNTOUCHED = sorted((set(D7.DISEASES) - D8.RECOMPUTE_SET) & set(D8.DISEASES))   # v8.1: cataract and back_pain left the table
RECOMP = [k for k in KEYS if k in D8.RECOMPUTE_SET]
hard_fail, inv = [], {}


def bytes_eq(a, b):
    return a.dtype == b.dtype and a.to_numpy().tobytes() == b.to_numpy().tobytes()


def same(a, b):
    return bool(np.array_equal(a.to_numpy(), b.to_numpy(), equal_nan=True)) if a.dtype.kind == "f" else bool(a.equals(b))


# ------------------------------------------------------------------ I1 untouched byte-identical
bad = [k + suf for k in UNTOUCHED for suf in S.OUTCOME_SUFFIXES if not bytes_eq(v8[k + suf], v7[k + suf])]
bad += [c for c in ["death_incident", "death_years", "death_prevalent", "psg_date", "censor_date", "fu_valid", "BDSPPatientID"]
        if not bytes_eq(v8[c], v7[c])]
inv["I1_untouched_byte_identical"] = {"n_conditions": len(UNTOUCHED), "pass": not bad, "failing_columns": bad}
if bad:
    hard_fail.append("I1 " + str(bad))

# ------------------------------------------------------------------ I2 exclusivity and flag consistency (all rows, all 57)
i2 = {}
for k in KEYS:
    fd, pv, ic = v8[f"{k}_first_date"], v8[f"{k}_prevalent"], v8[f"{k}_incident"]
    both = int(((pv + ic) > 1).sum())
    exp_pv = (fd.notna() & v8.psg_date.notna() & (fd <= v8.psg_date)).astype(int)
    exp_ic = (fd.notna() & v8.psg_date.notna() & (fd > v8.psg_date)).astype(int)
    bad_rows = (exp_pv != pv) | (exp_ic != ic)
    incons = int(bad_rows.sum())
    incons_coh = int(bad_rows[v8.index.isin(coh8.index)].sum())
    if both or incons:
        i2[k] = {"prevalent_and_incident": both, "flags_inconsistent_with_first_date_all_rows": incons,
                 "flags_inconsistent_with_first_date_cohort_rows": incons_coh, "inherited": k in UNTOUCHED}
inv["I2_prevalent_incident_exclusive_and_consistent"] = {"pass_exclusive": all(v["prevalent_and_incident"] == 0 for v in i2.values()),
                                                         "conditions_with_flags_inconsistent_with_first_date": i2}
if any(v["prevalent_and_incident"] for v in i2.values()):
    hard_fail.append("I2 prevalent and incident both 1: " + str({k: v["prevalent_and_incident"] for k, v in i2.items() if v["prevalent_and_incident"]}))
if any(v["flags_inconsistent_with_first_date_all_rows"] and not v["inherited"] for v in i2.values()):
    hard_fail.append("I2 recomputed flags inconsistent with first_date")
i2_coh = {k: v["flags_inconsistent_with_first_date_cohort_rows"] for k, v in i2.items() if v["flags_inconsistent_with_first_date_cohort_rows"]}
i2_off = {k: v["flags_inconsistent_with_first_date_all_rows"] - v["flags_inconsistent_with_first_date_cohort_rows"] for k, v in i2.items()
          if v["flags_inconsistent_with_first_date_all_rows"] > v["flags_inconsistent_with_first_date_cohort_rows"]}

# ------------------------------------------------------------------ I3 years arithmetic, cohort rows, recomputed conditions
i3 = {}
for k in RECOMP:
    fd, ic, yr = coh8[f"{k}_first_date"], coh8[f"{k}_incident"], coh8[f"{k}_years"]
    end = np.where(ic == 1, fd.values, coh8.censor_date.values)
    exp = (pd.Series(pd.to_datetime(end), index=coh8.index) - coh8.psg_date).dt.days / 365.25
    i3[k] = {"years_exact": same(exp, yr), "n_years_nan": int(yr.isna().sum()), "n_years_negative": int((yr < 0).sum())}
inv["I3_years_arithmetic_recomputed_cohort"] = i3
if not all(v["years_exact"] and v["n_years_nan"] == 0 and v["n_years_negative"] == 0 for v in i3.values()):
    hard_fail.append("I3 " + str({k: v for k, v in i3.items() if not (v["years_exact"] and v["n_years_nan"] == 0 and v["n_years_negative"] == 0)}))

# ------------------------------------------------------------------ I4 birth, I5 censor, I6 future (cohort rows, every condition)
dob = pd.concat([pd.read_csv(p, usecols=["BDSPPatientID", "DateOfBirth"]) for p in S.PSG_META.values()])
dob["DateOfBirth"] = pd.to_datetime(dob.DateOfBirth, errors="coerce")
dob = dob.dropna().drop_duplicates()
conflict = dob.groupby("BDSPPatientID").DateOfBirth.nunique()
n_conflict = int((conflict > 1).sum())
dobmap = dob.groupby("BDSPPatientID").DateOfBirth.min()
birth = coh8.BDSPPatientID.map(dobmap)
covered = int(birth.notna().sum())
approx = coh8.psg_date - pd.to_timedelta(coh8.AgeAtVisit * 365.25, unit="D")
birth_used = birth.fillna(approx - pd.Timedelta(days=366))
i4, i5, i6 = {}, {}, {}
for k in KEYS:
    fd, ic = coh8[f"{k}_first_date"], coh8[f"{k}_incident"]
    nb = int((fd < birth_used).sum())
    na = int(((ic == 1) & (fd > coh8.censor_date)).sum())
    nf = int((fd > TODAY).sum())
    tag = "inherited (untouched)" if k in UNTOUCHED else "recomputed"
    if nb:
        i4[k] = {"n": nb, "status": tag}
    if na:
        i5[k] = {"n": na, "status": tag, "n_died": int(((ic == 1) & (fd > coh8.censor_date) & (coh8.death_incident == 1)).sum())}
    if nf:
        i6[k] = {"n": nf, "status": tag}
inv["I4_first_date_before_birth_cohort"] = {"dob_source": "psg_metadata DateOfBirth (I0002, I0006)", "cohort_with_dob": covered,
                                            "cohort_n": int(len(coh8)), "dob_conflicts": n_conflict,
                                            "fallback_for_uncovered": "psg_date minus AgeAtVisit minus 366 days",
                                            "pass": not i4, "breaches": i4}
inv["I5_incident_first_date_after_censor_cohort"] = {"pass": not i5, "breaches": i5,
                                                     "total_inherited": sum(v["n"] for v in i5.values() if "inherited" in v["status"]),
                                                     "total_recomputed": sum(v["n"] for v in i5.values() if v["status"] == "recomputed")}
inv["I6_first_date_after_today_cohort"] = {"pass": not i6, "breaches": i6}

# ------------------------------------------------------------------ old vs new, the recomputed conditions
def state(df, k):
    return np.where(df[f"{k}_prevalent"] == 1, "prevalent", np.where(df[f"{k}_incident"] == 1, "incident", "none"))


ov = []
for k in RECOMP:
    new_s = state(coh8, k)
    old_s = state(coh7, k) if f"{k}_prevalent" in coh7 else np.array(["none"] * len(coh8))
    trans = pd.crosstab(pd.Series(old_s, name="old"), pd.Series(new_s, name="new"))
    moves = {f"{o}->{n}": int(trans.loc[o, n]) for o in trans.index for n in trans.columns if o != n and trans.loc[o, n]}
    prev, inc = int(coh8[f"{k}_prevalent"].sum()), int(coh8[f"{k}_incident"].sum())
    med = float(coh8.loc[coh8[f"{k}_prevalent"] == 0, f"{k}_years"].median())
    if f"{k}_prevalent" in coh7:
        oprev, oinc = int(coh7[f"{k}_prevalent"].sum()), int(coh7[f"{k}_incident"].sum())
        omed = float(coh7.loc[coh7[f"{k}_prevalent"] == 0, f"{k}_years"].median())
        pct = 100.0 * (inc - oinc) / oinc
    else:
        oprev = oinc = omed = pct = np.nan
    ov.append({"key": k, "label": D8.DISEASES[k][0],
               "status": "new" if k in D8.ADDED_2026_09_12 else "changed" if k in D8.CHANGED_2026_09_12 else "recomputed, list unchanged" if k == "osteoporosis" else "composite",
               "v7_prevalent": oprev, "v8_prevalent": prev, "v7_incident": oinc, "v8_incident": inc,
               "v7_median_fu_years": round(omed, 2) if omed == omed else np.nan, "v8_median_fu_years": round(med, 2),
               "patients_changing_status": int(sum(moves.values())), "transitions": json.dumps(moves),
               "pct_change_incident": round(pct, 1) if pct == pct else np.nan, "dramatic_over_20pct": bool(abs(pct) > 20) if pct == pct else False})
oldnew = pd.DataFrame(ov)
oldnew.to_csv(S.OUT_OLD_VS_NEW, index=False)

# ------------------------------------------------------------------ second code path: pandas prefix match on 300 cohort patients
rng = np.random.default_rng(20260912)
ids = sorted(rng.choice(np.sort(coh8.BDSPPatientID.values), 300, replace=False).tolist())
con = duckdb.connect()
con.execute(f"create view cond as select * from '{S.CACHE}'")
rows = con.execute(f"select person_id, condition_start_date as dt, condition_source_value as src from cond "
                   f"where person_id in ({','.join(map(str, ids))}) and condition_source_value is not null").df()
rows["code"] = rows.src.str.replace(".", "", regex=False).str.upper()
rows["dt"] = pd.to_datetime(rows.dt)
sub = coh8[coh8.BDSPPatientID.isin(ids)].set_index("BDSPPatientID")
smoke = {}
pandas_fd = {}
for k in RECOMP:
    if k == "cvd":
        continue
    lab, neg, i9, i10 = D8.DISEASES[k]
    m = rows.code.str.startswith(tuple(i9 + i10))
    for pat in D8.EXCLUDE.get(k, []):
        m &= ~rows.code.str.match("^" + pat.replace("_", "."))
    fd = rows[m].groupby("person_id").dt.min()
    pandas_fd[k] = sub.index.to_series().map(fd)
    smoke[k] = same(pandas_fd[k].reset_index(drop=True), sub[f"{k}_first_date"].reset_index(drop=True))
comp = pd.concat([pandas_fd["pad"]] + [sub[f"{c}_first_date"] for c in D8.CVD_COMPONENTS if c != "pad"], axis=1).min(axis=1)
smoke["cvd"] = same(comp.reset_index(drop=True), sub["cvd_first_date"].reset_index(drop=True))
# flags on the pandas path
for k in RECOMP:
    fd = comp if k == "cvd" else pandas_fd[k]
    pv = (fd.notna() & sub.psg_date.notna() & (fd <= sub.psg_date)).astype(int)
    ic = (fd.notna() & sub.psg_date.notna() & (fd > sub.psg_date)).astype(int)
    end = np.where(ic == 1, fd.values, sub.censor_date.values)
    yr = (pd.Series(pd.to_datetime(end), index=sub.index) - sub.psg_date).dt.days / 365.25
    smoke[k] = smoke[k] and same(pv, sub[f"{k}_prevalent"]) and same(ic, sub[f"{k}_incident"]) and same(yr, sub[f"{k}_years"])
inv["SMOKE_second_code_path_300_patients"] = {"n_patients": len(ids), "pass": all(smoke.values()), "per_condition": smoke}
if not all(smoke.values()):
    hard_fail.append("SMOKE " + str([k for k, v in smoke.items() if not v]))

# ------------------------------------------------------------------ alternatives (gate 12), cohort incident counts
con.execute("create table c as select person_id, condition_start_date as dt, condition_source_concept_id as cid, "
            "upper(replace(condition_source_value, '.', '')) as code from cond where condition_source_value is not null")


def inc_count(where, df=coh8):
    fd = con.execute(f"select person_id as BDSPPatientID, min(dt) as fd from c where {where} group by 1").df()
    fd["fd"] = pd.to_datetime(fd.fd)
    x = df[["BDSPPatientID", "psg_date"]].merge(fd, on="BDSPPatientID", how="left")
    return int(((x.fd.notna()) & (x.fd > x.psg_date)).sum()), int(((x.fd.notna()) & (x.fd <= x.psg_date)).sum())


def where_for(k, table_code="code"):
    lab, neg, i9, i10 = D8.DISEASES[k]
    w = " or ".join(f"{table_code} like '{p}%'" for p in sorted(set(i9 + i10)))
    ex = D8.EXCLUDE.get(k, [])
    if ex:
        w = f"({w}) and not (" + " or ".join(f"{table_code} like '{e}%'" for e in ex) + ")"
    return w


alt = {}
# A1 multi-code source values split on the comma (v7 rule matches the first code only)
con.execute("create table c2 as select person_id, dt, trim(u) as code from (select person_id, dt, unnest(string_split(code, ',')) as u from c)")
a1 = {}
for k in RECOMP:
    if k == "cvd":
        continue
    fd = con.execute(f"select person_id as BDSPPatientID, min(dt) as fd from c2 where {where_for(k)} group by 1").df()
    fd["fd"] = pd.to_datetime(fd.fd)
    x = coh8[["BDSPPatientID", "psg_date"]].merge(fd, on="BDSPPatientID", how="left")
    a1[k] = {"v8_incident": int(coh8[f"{k}_incident"].sum()), "split_rule_incident": int(((x.fd.notna()) & (x.fd > x.psg_date)).sum())}
alt["A1_multicode_split_rule"] = {"note": "50,673 cache rows hold a comma-separated code list; the v7 rule sees the first code only. Not applied: "
                                          "the 45 untouched conditions must stay identical to v7.", "counts": a1}
# A2 events after the censor date not counted as incident
alt["A2_censor_capped_incident"] = {k: {"v8_incident": int(coh8[f"{k}_incident"].sum()), "after_censor": i5.get(k, {}).get("n", 0),
                                        "capped_incident": int(coh8[f"{k}_incident"].sum()) - i5.get(k, {}).get("n", 0)} for k in RECOMP}
# A3 M80 reverse placement
w_ost = " or ".join(f"code like '{p}%'" for p in ["7330", "M81"])
w_fx = " or ".join(f"code like '{p}%'" for p in D8.DISEASES["fracture"][2] + D8.DISEASES["fracture"][3] + ["M80"])
alt["A3_M80_reverse"] = {"default": "M80 in osteoporosis only", "osteoporosis_without_M80_incident_prevalent": inc_count(w_ost),
                         "fracture_with_M80_incident_prevalent": inc_count(w_fx),
                         "default_osteoporosis": [int(coh8.osteoporosis_incident.sum()), int(coh8.osteoporosis_prevalent.sum())],
                         "default_fracture": [int(coh8.fracture_incident.sum()), int(coh8.fracture_prevalent.sum())]}
# A4 vocabulary-aware matching for the two untouched conditions with cross-system collisions (falls, thyroid_dis)
con.execute("create table c3 as select *, case when cid >= 44800000 and cid < 44900000 then 'ICD9' when cid > 0 then 'ICD10' else 'unknown' end as vocab from c")
val = con.execute("select sum(code similar to '[0-9].*' and vocab='ICD9')::double/sum(code similar to '[0-9].*'), "
                  "sum(code similar to '[A-DF-UW-Z].*' and vocab='ICD10')::double/sum(code similar to '[A-DF-UW-Z].*') from c3").fetchone()
a4 = {"vocabulary_rule": "source concept id 44,800,000 to 44,899,999 = ICD-9-CM, other non-zero = ICD-10-CM",
      "validation_digit_initial_codes_ICD9_share": round(val[0], 4), "validation_letter_initial_non_EV_codes_ICD10_share": round(val[1], 4)}
if val[0] > 0.99 and val[1] > 0.99:
    for k in ["falls", "thyroid_dis"]:
        lab, neg, i9, i10 = D8.DISEASES[k]
        w = "(" + " or ".join(f"code like '{p}%'" for p in i9) + ") and vocab='ICD9' or (" + \
            " or ".join(f"code like '{p}%'" for p in i10) + ") and vocab='ICD10'"
        fd = con.execute(f"select person_id as BDSPPatientID, min(dt) as fd from c3 where {w} group by 1").df()
        fd["fd"] = pd.to_datetime(fd.fd)
        x = coh8[["BDSPPatientID", "psg_date"]].merge(fd, on="BDSPPatientID", how="left")
        a4[k] = {"v7_incident": int(coh7[f"{k}_incident"].sum()), "v7_prevalent": int(coh7[f"{k}_prevalent"].sum()),
                 "vocab_aware_incident": int(((x.fd.notna()) & (x.fd > x.psg_date)).sum()),
                 "vocab_aware_prevalent": int(((x.fd.notna()) & (x.fd <= x.psg_date)).sum())}
else:
    a4["skipped"] = "vocabulary rule not validated on the data"
alt["A4_vocabulary_aware_untouched_collisions"] = a4

# ------------------------------------------------------------------ write
inv["hard_failures"] = hard_fail
with open(S.OUT_INVARIANTS, "w") as f:
    json.dump({"invariants": inv, "alternatives": alt, "old_vs_new": ov}, f, indent=1, default=str)
S.write_sidecar(S.OUT_INVARIANTS, INPUTS, __file__)
S.write_sidecar(S.OUT_OLD_VS_NEW, INPUTS, __file__)
side = json.load(open(S.OUT_PARQUET + ".provenance.json"))
recon = side["reconciliation_v7_lists_vs_v7_table"]

L = ["# OUTCOMES_V8 (lane X2, steps 040 and 041), " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M"), "",
     "## TLDR", "",
     f"Table written: {S.OUT_PARQUET} ({len(v8):,} rows, {v8.shape[1]} columns, {len(KEYS)} conditions plus death, provenance sidecar beside it). "
     f"Clock: v7 psg_date and censor_date (plan D2 default). Cache: the April-28 local hospital-record cache ({side['cache_rows']:,} rows). "
     f"{len(UNTOUCHED)} untouched conditions copied from v7 byte for byte ({'PASS' if inv['I1_untouched_byte_identical']['pass'] else 'FAIL'}); "
     f"{len(RECOMP)} recomputed (six fixed lists, four new outcomes, osteoporosis as an identity control, the cardiovascular composite). "
     f"Hard invariants: {'all pass' if not hard_fail else 'FAIL ' + str(hard_fail)}. Second code path on 300 patients: "
     f"{'agrees' if inv['SMOKE_second_code_path_300_patients']['pass'] else 'DISAGREES'}.", "",
     "## The recomputed outcomes, old (v7) against new (v8), analysis cohort of %s" % f"{len(coh8):,}", ""]
for r in ov:
    if r["status"] == "new":
        L.append(f"- {r['label']} ({r['key']}, NEW): prevalent {r['v8_prevalent']:,}, incident {r['v8_incident']:,}, "
                 f"median follow-up at risk {r['v8_median_fu_years']} y, patients entering a status {r['patients_changing_status']:,}.")
    else:
        L.append(f"- {r['label']} ({r['key']}, {r['status']}): prevalent {r['v7_prevalent']:,} to {r['v8_prevalent']:,}, incident {r['v7_incident']:,} to "
                 f"{r['v8_incident']:,} ({r['pct_change_incident']:+.1f} percent{', DRAMATIC' if r['dramatic_over_20pct'] else ''}), median follow-up at risk "
                 f"{r['v7_median_fu_years']} to {r['v8_median_fu_years']} y, patients changing status {r['patients_changing_status']:,} ({r['transitions']}).")
L += ["", "## Invariants", "",
      f"- I1 untouched conditions byte-identical to v7 in all four columns, plus death and the clock columns: {'PASS' if inv['I1_untouched_byte_identical']['pass'] else 'FAIL ' + str(bad)} ({len(UNTOUCHED)} conditions).",
      f"- I2 prevalent and incident mutually exclusive on every row of every condition: {'PASS' if inv['I2_prevalent_incident_exclusive_and_consistent']['pass_exclusive'] else 'FAIL'}. "
      f"Flags consistent with the stored first_date on the cohort rows: " + (f"not for {i2_coh} (inherited v7 columns: obesity_hypovent's v7 first_date column is "
      "the pre-2026-08-04 one, and 2 falls plus 2 obesity_hypovent patients were censor-capped in v7; see the reconciliation)" if i2_coh else "yes for all 57")
      + (f". Outside the cohort, {len(i2_off)} untouched conditions carry v7 flags on the 3,469 rows that have no study date in the v7 table (fu_valid 0), "
         "inherited and irrelevant to the analysis." if i2_off else "") + ".",
      f"- I3 years arithmetic exact, no NaN, none negative on the cohort rows of the {len(RECOMP)} recomputed conditions: {'PASS' if all(v['years_exact'] and v['n_years_nan']==0 and v['n_years_negative']==0 for v in i3.values()) else 'FAIL'}.",
      f"- I4 no first date before birth (DateOfBirth from the two site metadata files, {covered:,} of {len(coh8):,} cohort patients covered, {n_conflict} conflicting births; "
      f"the rest use psg_date minus age minus one year): {'PASS' if not i4 else 'BREACH ' + json.dumps(i4)}.",
      f"- I5 no incident first date after the censor date: {'PASS' if not i5 else 'BREACH, ' + str(inv['I5_incident_first_date_after_censor_cohort']['total_inherited']) + ' patient-conditions inherited in untouched conditions and ' + str(inv['I5_incident_first_date_after_censor_cohort']['total_recomputed']) + ' in recomputed ones'}. "
      "The v7 rule counts an event dated after the last contact as incident (the years then run past the censor date); the same rule is kept in v8 for consistency across all 57 conditions. Per condition: "
      + (", ".join(f"{k} {v['n']} ({v['status']}, {v['n_died']} died)" for k, v in sorted(i5.items())) if i5 else "none") + ".",
      f"- I6 no first date after today: {'PASS' if not i6 else 'BREACH ' + json.dumps(i6)}.",
      "", "## Reconciliation of the v7 lists against the v7 table on this cache (integrity gate 1)", "",
      f"- {recon['columns_identical']} of {recon['columns_total']} columns identical on all rows. The v7 twins of all {len(RECOMP)} recomputed conditions reproduce every cohort column: "
      f"{'YES' if recon['recomputed_twins_reproduced'] else 'NO'}, so the cache, the rule and the clock are the ones v7 used.",
      f"- Cohort flags or years differ in {recon['cohort_flags_or_years_differ']}: obesity_hypovent's v7 first_date column still carries the pre-2026-08-04 list (278.0x) on 22,613 of 28,873 rows while its flags were rebuilt with 278.03, "
      "and both conditions have 2 cohort patients whose event after the censor date was left non-incident in v7 (the rule as written makes them incident). Untouched here, copied byte for byte, listed for the owner.",
      f"- Differences confined to rows outside the cohort in {len(recon['off_cohort_rows_only'])} conditions: 3,469 rows with fu_valid 0 and no study date keep a years value in v7 that the rule gives as NaN. No effect on the analysis cohort.",
      "", "## Alternatives (gate 12), cohort incident counts", "",
      "- A1 comma-separated code lists matched on every code, not the first: " + "; ".join(f"{k} {v['v8_incident']:,} to {v['split_rule_incident']:,}" for k, v in a1.items()) + ".",
      "- A2 events after the last contact not counted (censor-capped): " + "; ".join(f"{k} {v['v8_incident']:,} to {v['capped_incident']:,}" for k, v in alt['A2_censor_capped_incident'].items() if v['after_censor']) + ".",
      f"- A3 M80 placed in fracture instead of osteoporosis: osteoporosis incident {alt['A3_M80_reverse']['default_osteoporosis'][0]:,} to {alt['A3_M80_reverse']['osteoporosis_without_M80_incident_prevalent'][0]:,}, "
      f"fracture incident {alt['A3_M80_reverse']['default_fracture'][0]:,} to {alt['A3_M80_reverse']['fracture_with_M80_incident_prevalent'][0]:,}.",
      "- A4 inherited cross-system collisions in two untouched lists (ICD-9 E-code fall prefixes E888x/E884x also match ICD-10 E88.8x/E88.4x metabolic codes; ICD-10 E00 to E03 thyroid prefixes also match ICD-9 activity codes E001 to E030), "
      "vocabulary read from the source concept id (validated: %s of digit-initial codes are ICD-9, %s of letter-initial codes are ICD-10): " % (a4['validation_digit_initial_codes_ICD9_share'], a4['validation_letter_initial_non_EV_codes_ICD10_share'])
      + "; ".join(f"{k} incident {a4[k]['v7_incident']:,} to {a4[k]['vocab_aware_incident']:,}, prevalent {a4[k]['v7_prevalent']:,} to {a4[k]['vocab_aware_prevalent']:,}" for k in ["falls", "thyroid_dis"] if k in a4) + ". Not applied, owner's call.",
      "", "## Not done, and why", "",
      "- numbers/disease_definitions.py is not modified and numbers/audit_prefixes_v5.py does not exist on disk; the v8 module and its test live in X2_outcomes/scripts/ (lane instruction: do not touch numbers/). Phase E must import disease_definitions_v8 from there or copy it into numbers/ under the plan's step 040.",
      "- Outputs are under X2_outcomes/work/ and X2_outcomes/logs/, not V8/outcomes/ (lane instruction).",
      "- The record-frame clock at I0006 (plan B5, D2 alternative) was not applied; every v8 outcome sits on the v7 clock.",
      "", "## Files", "", f"- {S.OUT_PARQUET} (+ .provenance.json)", f"- {S.OUT_COUNTS} (+ sidecar), the regenerated reference for later asserts",
      f"- {S.OUT_OLD_VS_NEW}, {S.OUT_INVARIANTS}, {S.OUT_RECONCILE}, {S.DIFF_MD}", ""]
with open(S.REPORT_MD, "w") as f:
    f.write("\n".join(L))
say(f"check_outcomes_v8: hard failures {hard_fail or 'none'}; I4 breaches {sum(v['n'] for v in i4.values())}, I5 breaches {sum(v['n'] for v in i5.values())}, "
    f"I6 breaches {sum(v['n'] for v in i6.values())}; smoke {'PASS' if inv['SMOKE_second_code_path_300_patients']['pass'] else 'FAIL'}; report -> {S.REPORT_MD}")
sys.exit(1 if hard_fail else 0)
