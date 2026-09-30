"""
Materialise the two replacement negative controls into the frozen cohort.

Contact dermatitis and hemorrhoids entered the panel on 2026-08-07 but existed only inside the
25-candidate screening script, which computed them on the fly. Every downstream analysis reads
data_frozen/t90_final.parquet, so the two conditions have to be columns in that file or they can
never be fitted for the lag ladder, the stage windows, the duration cells, the subgroups or the
treatment analysis, and every floor outside the primary per-SD scale would be a maximum over
three controls instead of five.

Additive only. Six new columns per condition; no existing column is read back or rewritten. The
pre-change file is kept at data_frozen/_superseded_t90_final_2026-08-07_pre_negcontrol_v5.parquet.

THE ICD PREFIX TRAP is checked, not assumed. These archives store the source value with the
decimal point stripped, so ICD-9 E888.1 and ICD-10 E88.81 both become E8881 and a prefix cannot
tell them apart. That is what corrupted the falls definition. The guard below is the one from
negcontrols.py: no prefix may begin with E or V, no digit-initial prefix may reach a
letter-initial code, no letter-initial prefix may reach a digit-initial code, and every matched
token must be covered by a stated prefix. It is first proved to fire on E8881 itself.

Matching is on the FIRST comma-separated code of the stored value, because a small share of rows
pack several codes into one string and only the first can be reached by a prefix.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../..")))  # repository root, where paths.py lives
import paths
import sys
import warnings

warnings.filterwarnings("ignore")

import duckdb
import numpy as np
import pandas as pd

ROOT = paths.T90_ROOT
CACHE = f"{paths.OMOP_CACHE_DIR}/condition_occurrence_cohort.parquet"
sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
from disease_definitions import DISEASES, NEGATIVE_CONTROLS   # noqa: E402

NEW = ["dermatitis_contact", "haemorrhoids"]
assert all(k in NEGATIVE_CONTROLS for k in NEW)

con = duckdb.connect()
con.execute(f"""create table c as
    select person_id, condition_start_date as dt,
           upper(replace(condition_source_value, '.', '')) as code
    from '{CACHE}' where condition_source_value is not null""")
print(f"condition rows {con.execute('select count(*) from c').fetchone()[0]:,}")


def audit_prefixes(prefixes):
    """Return (first tokens reached, violations) for a prefix set."""
    viol = []
    for p in prefixes:
        if p[0] in ("E", "V"):
            viol.append(f"{p}: begins with {p[0]}, the only letters ICD-9 shares with ICD-10 "
                        "once the decimal point is stripped")
    where = " or ".join(f"code like '{p}%'" for p in prefixes)
    d = con.execute(f"select split_part(code, ',', 1) tok, count(*) n from c where {where} "
                    "group by 1 order by 1").df()
    toks = d.tok.tolist()
    for p in prefixes:
        hit = [t for t in toks if t.startswith(p)]
        if p[0].isdigit() and any(t[0].isalpha() for t in hit):
            viol.append(f"{p}: digit-initial prefix reaches letter-initial codes "
                        f"{[t for t in hit if t[0].isalpha()][:8]}")
        if p[0].isalpha() and any(t[0].isdigit() for t in hit):
            viol.append(f"{p}: letter-initial prefix reaches digit-initial codes "
                        f"{[t for t in hit if t[0].isdigit()][:8]}")
    orphan = [t for t in toks if not any(t.startswith(p) for p in prefixes)]
    if orphan:
        viol.append(f"tokens matched but not covered by any stated prefix: {orphan[:8]}")
    return toks, dict(zip(d.tok, d.n)), viol


# ---- prove the guard fires on the pair that corrupted falls ----------------------------
_, hits, trap = audit_prefixes(["E8881"])
assert trap, "the prefix guard failed on E8881, the known ICD-9 / ICD-10 collision"
print("\nprefix guard, proved on the known trap")
print(f"  E8881 reaches {hits}")
print(f"  verdict: {trap[0]}")

# ---- audit the two new prefix sets, code by code ---------------------------------------
print("\nprefix audit, every new code")
for k in NEW:
    label, is_neg, i9, i10 = DISEASES[k]
    pref = sorted(set(i9) | set(i10))
    for p in pref:
        toks, hits, viol = audit_prefixes([p])
        assert not viol, f"{k} prefix {p}: {viol}"
        print(f"  {label:<20} {p:<6} -> {len(toks):3d} distinct codes, "
              f"{sum(hits.values()):8,} rows   e.g. {toks[:5]}")
    toks, hits, viol = audit_prefixes(pref)
    assert not viol, f"{k}: {viol}"
    print(f"  {label:<20} {str(pref):<30} 0 violations, {len(toks)} codes total")

# ---- build the columns -----------------------------------------------------------------
base = pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet")
print(f"\nfrozen cohort rows {len(base):,}, columns {len(base.columns)}")
for k in NEW:
    assert f"{k}_incident" not in base.columns, f"{k} already present, refusing to overwrite"

psg = pd.to_datetime(base.psg_date, errors="coerce")
cen = pd.to_datetime(base.censor_date, errors="coerce")

for k in NEW:
    label, _, i9, i10 = DISEASES[k]
    pref = sorted(set(i9) | set(i10))
    where = " or ".join(f"code like '{p}%'" for p in pref)
    fd = con.execute(f"select person_id, min(dt) as fdt from c where {where} "
                     "group by person_id").df()
    fd.columns = ["BDSPPatientID", "fdt"]
    fd["fdt"] = pd.to_datetime(fd.fdt, errors="coerce")
    m = base[["BDSPPatientID"]].merge(fd, on="BDSPPatientID", how="left")
    f0 = m.fdt
    # Identical rule to rebuild_outcomes.py, which built the other 49 conditions in this file:
    # incident is a first record strictly after the study, with no cap at the censoring date.
    # The cap looks more correct and is what refresh_two_conditions.py used, but it is applied
    # to only two conditions in the frozen file (obesity_hypovent and falls) and NOT to back
    # pain, cataract or glaucoma. A control panel has to be internally comparable before it is
    # individually ideal, so the two replacements are built the way the three retained controls
    # were built. Applying the cap instead moves hemorrhoids from 1,445 to 1,443 events and
    # contact dermatitis not at all.
    prev = (f0.notna() & psg.notna() & (f0 <= psg)).astype(int)
    inc = (f0.notna() & psg.notna() & (f0 > psg)).astype(int)
    end = np.where(inc == 1, f0.values, cen.values)
    base[f"{k}_first_date"] = f0.values
    base[f"{k}_prevalent"] = prev.values
    base[f"{k}_incident"] = inc.values
    base[f"{k}_years"] = (pd.to_datetime(end) - psg).dt.days / 365.25
    print(f"\n{label}: prevalent {int(prev.sum()):,}  incident {int(inc.sum()):,}  "
          f"(whole frozen file, before the cohort filter)")

base.to_parquet(f"{paths.TABLES_DIR}/t90_final.parquet", index=False)
print(f"\nwritten -> data_frozen/t90_final.parquet, columns now {len(base.columns)}")

# ---- read back in the analysis cohort and check against the screening script ------------
from cohort_spec import apply_cohort   # noqa: E402
b = apply_cohort(pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet"))
screen = pd.read_csv(f"{paths.NUMBERS_DIR}/negcontrols_v4.csv").set_index("key")
print(f"\nanalysis cohort {len(b):,}")
print(f"{'condition':<20}{'events here':>13}{'events in screen':>19}{'match':>8}")
for k in NEW:
    f = b[(b[f"{k}_prevalent"] == 0) & b[f"{k}_years"].notna() & (b[f"{k}_years"] > 0)]
    ev, want = int(f[f"{k}_incident"].sum()), int(screen.loc[k, "events"])
    print(f"{DISEASES[k][0]:<20}{ev:>13,}{want:>19,}{'  yes' if ev == want else '  NO':>8}")
    assert ev == want, f"{k}: {ev} incident events here against {want} in the screen"
print("\nboth conditions reproduce the screening script's event counts exactly.")
