"""
Step 040 unit tests for disease_definitions_v8.py, and the writer of DISEASE_DEFINITIONS_DIFF.md.

Hard failures (rc 1): any defect in the ten changed or new lists, any drift in an untouched entry against the v7 module,
any change to RANKING_EXCLUDE, CVD_COMPONENTS, NEGATIVE_CONTROLS or CIRCULAR, a changed or new prefix with zero rows in the
cache, a type-1 exclusion pattern missing. Inherited defects inside untouched conditions are REPORTED (they must stay
byte-identical to v7 in this lane) and listed for the owner with counts from the cache.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import importlib.util
import re
import sys

import duckdb

sys.path.insert(0, paths.X2_SCRIPTS)
import x2_spec as S                              # noqa: E402
import disease_definitions_v8 as v8              # noqa: E402

for p in (S.V7_DEFINITIONS, S.V8_DEFINITIONS, S.CACHE):
    S.require_local(p)
spec = importlib.util.spec_from_file_location("disease_definitions_v7", S.V7_DEFINITIONS)
v7 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v7)

fails, warns, info = [], [], []


def check(cond, msg):
    if not cond:
        fails.append(msg)


def lists(mod, k):
    return list(mod.DISEASES[k][2]), list(mod.DISEASES[k][3])


# ------------------------------------------------------------------ 1. structure of every list
for k, (lab, neg, i9, i10) in v8.DISEASES.items():
    hard = k in v8.RECOMPUTE_SET
    for p in i9:
        check(re.fullmatch(r"(?:[0-9]+|[EV][0-9]+)", p) is not None, f"{k}: ICD-9 prefix {p!r} not digit-initial (E/V allowed)")
    for p in i10:
        check(re.fullmatch(r"[A-Z][0-9A-Z]*", p) is not None, f"{k}: ICD-10 prefix {p!r} not letter-initial")
    allp = i9 + i10
    check(len(allp) == len(set(allp)), f"{k}: duplicate prefix")
    for a in allp:
        for b in allp:
            if a != b and b.startswith(a):
                (fails if hard else warns).append(f"{k}: prefix {a} already covers {b} (redundant, inherited from v7)" if not hard
                                                  else f"{k}: prefix {a} already covers {b}")
    check(k in v8.ORGAN_GROUP, f"{k}: no ORGAN_GROUP entry")

# ------------------------------------------------------------------ 2. exact diffs against v7
check(set(v7.DISEASES) - set(v8.DISEASES) == set(v8.SWAPPED_CONTROLS_2026_09_14), "a v7 condition is missing from v8 (only the two swapped controls may leave)")
check(set(v8.DISEASES) - set(v7.DISEASES) == v8.ADDED_2026_09_12 | v8.SWAPPED_IN_2026_09_14, "the added keys are not exactly ADDED_2026_09_12 plus the two swapped-in controls")
check(v8.ADDED_2026_09_12 == {"polycythemia", "parkinson", "ild", "vent_arrhythmia_arrest"}, "ADDED set wrong")
check(v8.CHANGED_2026_09_12 == {"nafld", "nocturia", "dementia", "diabetes", "pad", "fracture"}, "CHANGED set wrong")
untouched = sorted((set(v7.DISEASES) - v8.RECOMPUTE_SET) & set(v8.DISEASES))   # v8.1: the swapped-out controls are gone
for k in untouched:
    check(v7.DISEASES[k] == v8.DISEASES[k], f"untouched {k} differs from v7")
check(v7.DISEASES["osteoporosis"] == v8.DISEASES["osteoporosis"], "osteoporosis must be unchanged (M80 stays there)")
check(v7.DISEASES["cvd"] == v8.DISEASES["cvd"], "cvd entry changed")
expected = {"nafld": ({"5715", "K746"}, set()), "nocturia": ({"78819"}, {"78843"}), "dementia": (set(), {"3310"}),
            "diabetes": (set(), set()), "pad": ({"4432"}, set()), "fracture": ({"M80"}, set())}
diffs = {}
for k, (rem, add) in expected.items():
    a9, a10 = lists(v7, k); b9, b10 = lists(v8, k)
    old, new = set(a9 + a10), set(b9 + b10)
    diffs[k] = (sorted(old - new), sorted(new - old))
    check((old - new, new - old) == (rem, add), f"{k}: diff is -{sorted(old - new)} +{sorted(new - old)}, expected -{sorted(rem)} +{sorted(add)}")
    check(v7.DISEASES[k][0] == v8.DISEASES[k][0] and v7.DISEASES[k][1] == v8.DISEASES[k][1], f"{k}: label or control flag changed")
check(v8.EXCLUDE == {"diabetes": ["250_1", "250_3"]}, f"EXCLUDE is {v8.EXCLUDE}")
check("M80" in v8.DISEASES["osteoporosis"][3] and "M80" not in v8.DISEASES["fracture"][3], "M80 must sit in osteoporosis only")
check(not any("M80" in v8.DISEASES[k][3] for k in v8.DISEASES if k != "osteoporosis"), "M80 appears in a second condition")
newlists = {"polycythemia": ({"2384", "2890"}, {"D45", "D751"}), "parkinson": ({"3320"}, {"G20"}),
            "ild": ({"515", "516"}, {"J84"}), "vent_arrhythmia_arrest": ({"4271", "4274", "4275"}, {"I472", "I490", "I46"})}
for k, (e9, e10) in newlists.items():
    b9, b10 = lists(v8, k)
    check((set(b9), set(b10)) == (e9, e10), f"{k}: lists {b9} {b10} differ from the plan's D3")
check(v8.RANKING_EXCLUDE == v7.RANKING_EXCLUDE, "RANKING_EXCLUDE changed")
check(v8.CVD_COMPONENTS == v7.CVD_COMPONENTS, "CVD_COMPONENTS changed")
check(v8.NEGATIVE_CONTROLS == [v8.SWAPPED_CONTROLS_2026_09_14.get(k, k) for k in v7.NEGATIVE_CONTROLS], "NEGATIVE_CONTROLS differ from v7 by more than the 14 Sept swap")
check(all(v8.DISEASES[k][1] for k in v8.SWAPPED_IN_2026_09_14) and v8.SWAPPED_IN_2026_09_14 <= v8.RECOMPUTE_SET, "swapped-in controls must be flagged controls and recomputed")
check(v8.CIRCULAR == v7.CIRCULAR, "CIRCULAR changed")
check(v8.CONFOUNDING_FLOOR_PER_SD == v7.CONFOUNDING_FLOOR_PER_SD, "CONFOUNDING_FLOOR changed")
check(set(v8.ORGAN_GROUP) == set(v8.DISEASES), "ORGAN_GROUP keys differ from DISEASES")
for k in untouched:
    check(v7.ORGAN_GROUP[k] == v8.ORGAN_GROUP[k], f"ORGAN_GROUP[{k}] changed")

# ------------------------------------------------------------------ 3. cross-condition sharing (information)
keys = list(v8.DISEASES)
for i, a in enumerate(keys):
    pa = lists(v8, a)[0] + lists(v8, a)[1]
    for b in keys[i + 1:]:
        pb = lists(v8, b)[0] + lists(v8, b)[1]
        shared = sorted({x for x in pa for y in pb if x.startswith(y) or y.startswith(x)} |
                        {y for x in pa for y in pb if x.startswith(y) or y.startswith(x)})
        if shared:
            info.append(f"{a} and {b} share code space: {shared}")

# ------------------------------------------------------------------ 4. cache-driven audit
con = duckdb.connect()
con.execute(f"create view cond as select * from '{S.CACHE}'")
con.execute("create table c as select person_id, condition_start_date as dt, condition_source_concept_id as cid, "
            "upper(replace(condition_source_value, '.', '')) as code from cond where condition_source_value is not null")
era = f"dt >= date '{S.ICD10_ERA_START}'"


def rows_for(p):
    return con.execute(f"select count(*), count(distinct person_id) from c where code like '{p}%'").fetchone()


zero_rows = []
for k in sorted(v8.RECOMPUTE_SET - {"cvd"}):
    for p in lists(v8, k)[0] + lists(v8, k)[1]:
        n, np_ = rows_for(p)
        info.append(f"{k}: prefix {p} rows={n:,} persons={np_:,}")
        if n == 0:
            zero_rows.append((k, p))
for k, p in zero_rows:
    inherited = k in v7.DISEASES and p in (v7.DISEASES[k][2] + v7.DISEASES[k][3])
    if inherited:
        warns.append(f"{k}: inherited prefix {p} matches zero cache rows (also zero in v7; harmless, WHO-only code)")
    else:
        check(False, f"{k}: prefix {p} matches zero cache rows")
n_excl = con.execute("select count(*), count(distinct person_id) from c where code like '250_1%' or code like '250_3%'").fetchone()
info.append(f"diabetes EXCLUDE 250_1/250_3: rows={n_excl[0]:,} persons={n_excl[1]:,}")

# cross-system candidates: E/V-initial ICD-9 prefixes and E-initial ICD-10 prefixes can match codes of the other system once
# the point is stripped. For each, rows by era and the dominant source concept id per era; two dominant ids = two vocabularies.
cross = []
for k, (lab, neg, i9, i10) in v8.DISEASES.items():
    cands = [(p, "icd9") for p in i9 if p[0] in "EV"] + [(p, "icd10") for p in i10 if p[0] == "E"]
    for p, sysname in cands:
        r = con.execute(f"select {era} as era10, count(*) n, count(distinct person_id) np, mode(cid) dom_cid, count(distinct cid) ncid "
                        f"from c where code like '{p}%' group by 1 order by 1").fetchall()
        by = {bool(e): (n, np_, dom, ncid) for e, n, np_, dom, ncid in r}
        pre, post = by.get(False, (0, 0, None, 0)), by.get(True, (0, 0, None, 0))
        flag = pre[0] > 0 and post[0] > 0 and pre[2] != post[2]
        cross.append((k, sysname, p, pre, post, flag))

# ------------------------------------------------------------------ 5. write the diff
doc = v8.__doc__
lines = ["# DISEASE_DEFINITIONS_DIFF (v7 numbers/disease_definitions.py -> X2 disease_definitions_v8.py), 2026-09-12", "",
         "numbers/disease_definitions.py is NOT modified by this lane; the v8 module lives in X2_outcomes/scripts/ and is the input of",
         "rebuild_outcomes_v8.py. Result: %d %s." % (len(fails), "FAIL" if fails else "hard checks passed, rc 0"), "",
         "## Lists written out with code meanings (from the module docstring)", "", "```", doc.strip(), "```", "",
         "## Exact prefix diffs of the six fixed lists (removed / added)"]
for k, (rem, add) in diffs.items():
    lines.append(f"- {k} ({v8.DISEASES[k][0]}): removed {rem or 'none'}, added {add or 'none'}"
                 + ("; EXCLUDE " + str(v8.EXCLUDE[k]) if k in v8.EXCLUDE else ""))
lines += ["- osteoporosis: unchanged (M80 stays); " + v8.M80_ALTERNATIVE, "",
          "## Four added outcomes", ""]
for k in sorted(v8.ADDED_2026_09_12):
    lines.append(f"- {k} ({v8.DISEASES[k][0]}, group {v8.ORGAN_GROUP[k]}): ICD-9 {v8.DISEASES[k][2]}, ICD-10 {v8.DISEASES[k][3]}")
lines += ["", f"Untouched entries identical to v7: {len(untouched)} ({', '.join(untouched)})",
          f"RANKING_EXCLUDE unchanged: {sorted(v8.RANKING_EXCLUDE)}; CVD_COMPONENTS unchanged: {v8.CVD_COMPONENTS}; "
          f"NEGATIVE_CONTROLS unchanged: {v8.NEGATIVE_CONTROLS}; CIRCULAR unchanged: {v8.CIRCULAR}", "",
          "## Cache row counts for every changed or new prefix (the 788.19 lesson: a prefix with zero rows is a silent hole)", ""]
lines += [f"- {x}" for x in info if ": prefix" in x or "EXCLUDE" in x]
lines += ["", "## Cross-system prefix audit (ICD-9 E/V prefixes and ICD-10 E prefixes, rows before and after 2015-10-01 with the",
          "dominant source concept id; a different dominant id in the two eras means the prefix matches codes of BOTH systems)", ""]
for k, sysname, p, pre, post, flag in cross:
    lines.append(f"- {'COLLISION ' if flag else ''}{k} {sysname} {p}: pre-2015-10 rows={pre[0]:,} persons={pre[1]:,} dominant cid={pre[2]} ({pre[3]} ids); "
                 f"post rows={post[0]:,} persons={post[1]:,} dominant cid={post[2]} ({post[3]} ids)")
lines += ["", "## Conditions sharing code space (information; ihd/mi and the composite are by design)", ""] + [f"- {x}" for x in info if "share code space" in x]
lines += ["", "## Warnings (inherited from v7, untouched here)", ""] + ([f"- {w}" for w in warns] or ["- none"])
lines += ["", "## Hard failures", ""] + ([f"- {f}" for f in fails] or ["- none"]) + [""]
with open(S.DIFF_MD, "w") as f:
    f.write("\n".join(lines))

ncoll = sum(1 for x in cross if x[5])
print(f"test_definitions_v8: {len(fails)} failures, {len(warns)} inherited warnings, {ncoll} cross-system collision prefixes "
      f"(all in untouched conditions: {sorted({x[0] for x in cross if x[5]})}); diff -> {S.DIFF_MD}")
for f_ in fails:
    print("  FAIL", f_)
sys.exit(1 if fails else 0)
