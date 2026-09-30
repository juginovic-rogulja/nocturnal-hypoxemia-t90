"""
The one source of truth for which family each of the ranked measurements belongs to.

Every consumer that needs a family must call family_of() and read
numbers/measure_families.csv through it. Nothing may classify a measurement by matching
substrings of its column name.

v8 (2026-09-12): the row count is cohort_spec.N_MEASURES (141: the 60 per-site standardized
twins dropped, hypoxic burden, the periodic limb movement index, sleep onset latency and REM
latency added), the file path is cohort_spec.MEASURE_FAMILIES_CSV, and "Limb movements" is a
seventh family. The family-count check is structural (the counts sum to N_MEASURES and the
family set equals FAMILY_ORDER), never a typed dictionary of counts.

WHY THIS FILE EXISTS
--------------------
Until 2026-08-07 four figure scripts each carried their own copy of a keyword cascade:

    if any(k in f for k in ("spo2","odi","sao2","desat")):        return "Oxygen"
    if any(k in f for k in ("hr_bpm","rmssd","sdnn","pnn50","lf_","hf_")): return "Heart"
    ...
    return "Sleep timing and structure"

It is first-match-wins on substrings with a catch-all default, so a name that matches nothing
becomes "Sleep timing and structure" silently and a name that collides with an unrelated word
is captured by whichever line runs first. Three of the (then) 197 were wrong:

  wholenight_n_rr   the count of R-R intervals, returned by heart_features.hrv_time_domain
                    beside RMSSD and SDNN. The cardiac keyword list carries hr_bpm, rmssd,
                    sdnn, pnn50, lf_ and hf_ but not n_rr, so it fell through to the default
                    and was drawn as sleep timing and structure, where it was the largest
                    point in the family.
  nrem_n_rr         same measurement restricted to NREM, same silent default.
  recording_dur_min epochs x 30 s, an acquisition property computed by
                    atlas_arch_pass/arch_worker.extract_arch. It matched "rdi", the
                    respiratory disturbance index, inside the letters of recoRDIng, and was
                    drawn as a breathing event, the only negative point in that family.

Four cardiac measurements are correct in the cascade only because of line order:
nrem_lf_power, wholenight_lf_power, nrem_hf_power and wholenight_hf_power all match the
brain-spectral list on "_power" and are saved only because the cardiac line runs first.
Nothing in the old code recorded that dependency. That is the reason the cascade is not
repaired in place but withdrawn.

USE
---
    from measure_family_lookup import family_of, FAMILY_ORDER
    df["family"] = df.feature.map(family_of)                       # long names, eFigure 5
    df["family"] = df.feature.map(lambda f: family_of(f, "make_v2"))   # short names, Fig 1B

Schemes are alias spellings of the same partition, so each consumer keeps its own wording
while reading one table:

    "efig05"  Oxygenation / Heart rate and variability / Brain, microstructure /
              Brain, spectral power / Breathing events / Sleep timing and structure /
              Limb movements
    "make_v2" Oxygen / Heart / Brain, microstructure / Brain, spectral / Breathing events /
              Sleep timing and structure / Limb movements
    "audit"   Oxygenation / Cardiac, HRV / EEG microstructure / EEG spectral /
              Breathing events / Sleep architecture / Limb movements
    "etable"  Oximetry / Autonomic / ... / Limb movements

An unknown feature name raises. It is never assigned a default family.
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cohort_spec import MEASURE_FAMILIES_CSV, N_MEASURES  # noqa: E402  (v8 sweep 2026-09-12)

CSV = Path(MEASURE_FAMILIES_CSV)

_SCHEME_COLUMN = {
    "efig05": "family",          # the long names eFigure 5 prints
    "family": "family",
    "make_v2": "family_make_v2",  # the short names Figure 1B prints
    "audit": "family_audit",
    "etable": "family_etable",
}

_T = pd.read_csv(CSV)
assert len(_T) == N_MEASURES, f"measure_families.csv should hold {N_MEASURES} rows, holds {len(_T)}"
assert _T.feature.is_unique, "duplicate feature in measure_families.csv"

_MAP = {s: dict(zip(_T.feature, _T[c])) for s, c in _SCHEME_COLUMN.items()}

# the panel order used by every consumer that does not sort by a statistic
FAMILY_ORDER = ["Oxygenation", "Breathing events", "Limb movements", "Heart rate and variability",
                "Brain, microstructure", "Brain, spectral power",
                "Sleep timing and structure"]

FAMILY_COUNTS = _T.family.value_counts().to_dict()
assert sum(FAMILY_COUNTS.values()) == N_MEASURES and set(FAMILY_COUNTS) == set(FAMILY_ORDER), \
    (FAMILY_COUNTS, FAMILY_ORDER)


def family_of(feature, scheme="efig05"):
    """The family of one measurement. Raises on anything not in the frozen list."""
    try:
        table = _MAP[scheme]
    except KeyError:
        raise KeyError(f"unknown family scheme {scheme!r}, "
                       f"expected one of {sorted(_SCHEME_COLUMN)}")
    try:
        return table[feature]
    except KeyError:
        raise KeyError(
            f"{feature!r} is not one of the {N_MEASURES} ranked measurements in {CSV.name}. "
            "Add it to that file by hand from the extraction code. Do NOT guess it from "
            "the column name: substring matching is what put wholenight_n_rr in sleep "
            "timing and recording_dur_min in breathing events.")


def members(family, scheme="efig05", distinct_only=False):
    """Every feature in one family, ranking order. distinct_only drops the _siteZ twins
    (none remain in v8; the flag is kept so old callers keep working)."""
    col = _SCHEME_COLUMN[scheme]
    t = _T[_T[col] == family]
    if distinct_only:
        t = t[~t.is_siteZ_duplicate.astype(bool)]
    return t.feature.tolist()


def table():
    """The whole lookup, as read from disk."""
    return _T.copy()
