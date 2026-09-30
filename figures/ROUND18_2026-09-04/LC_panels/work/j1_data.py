"""JOB 1 data: sleep duration x T90 band distribution, plus the Figure 4b control.

Cohort: data_frozen/t90_final.parquet through numbers/cohort_spec.apply_cohort (19,173).
Columns found by inspection of the frame, not guessed:
  TST_min            total sleep time on the study night, minutes  -> duration group
  spo2_pct_below_90  T90, percent of the recording below 90%       -> oxygen band
  AHI                apnea-hypopnea index                          -> apnea severity (control)

Duration cut points are the paper's own, read out of the frozen model that produced Figure 3b
(Sleep_Variability_2026-08/dur_x_oxygen_6cell_tst300/model.py):
  short  TST < 300 min          (<5 h)
  normal 360 <= TST <= 420 min  (6-7 h, inclusive both ends)
  long   TST > 420 min          (>7 h)
The 300-360 min band and the 6 recordings with no TST are outside the paper's duration groups
and are not shown, exactly as in Figure 3b.

T90 bands are Figure 4b's own: <=1, >1-5, >5-10, >10 percent of the recording.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../..")))  # repository root, where paths.py lives
import paths
import json, os, sys
from decimal import Decimal, ROUND_HALF_UP
import numpy as np
import pandas as pd

ROOT = paths.FIGURE_ROOT
sys.path.insert(0, paths.ANALYSIS_DIR)
from cohort_spec import apply_cohort, COHORT_N            # noqa: E402  # v8 sweep 2026-09-12

T90_EDGES = [(1.0, "≤1%"), (5.0, ">1–5%"), (10.0, ">5–10%"), (None, ">10%")]
T90_KEYS = [e[1] for e in T90_EDGES]
DUR_ROWS = ["<5 h", "6–7 h", ">7 h"]
AHI_ROWS = ["No apnea", "Mild", "Moderate", "Severe"]
AHI_SUB = ["(AHI <5)", "(AHI 5 to <15)", "(AHI 15 to <30)", "(AHI 30 or more)"]
PCT_BANDS = [(5.0, "<5"), (15.0, "5–15"), (35.0, "15–35"), (None, "≥35")]


def pct1(num, den):
    """Percent to one decimal, half-up, from the exact integer ratio."""
    return str((Decimal(num) * 100 / Decimal(den)).quantize(Decimal("0.1"), rounding=ROUND_HALF_UP))


def band_index(p):
    for i, (hi, _) in enumerate(PCT_BANDS):
        if hi is None or float(p) < hi:
            return i
    raise AssertionError(p)


def t90_band(v):
    for i, (hi, _) in enumerate(T90_EDGES):
        if hi is None or v <= hi:
            return i
    raise AssertionError(v)


def load():
    b = apply_cohort(pd.read_parquet(os.path.join(paths.TABLES_DIR, "t90_final.parquet")))
    assert len(b) == COHORT_N, len(b)
    assert b.spo2_pct_below_90.notna().all()
    ox = b.spo2_pct_below_90.map(t90_band)
    return b, ox


def crosstab(b, ox, group, order):
    ct = pd.crosstab(group, ox).reindex(index=order, columns=range(4)).fillna(0).astype(int)
    out = {}
    for r in order:
        n = int(ct.loc[r].sum())
        cells = []
        for c in range(4):
            k = int(ct.loc[r, c])
            p = pct1(k, n)
            cells.append(dict(band=T90_KEYS[c], n=k, pct=p,
                              label=f"{k:,} ({p})", shade=band_index(p)))
        assert sum(c["n"] for c in cells) == n
        out[r] = dict(total=n, cells=cells)
    return out


def build():
    b, ox = load()
    ahi = b.AHI
    has_ahi = ahi.notna()   # v7: the index is undefined on nights with under 60 min of true sleep
    apnea = pd.Series(np.where(ahi < 5, AHI_ROWS[0], np.where(ahi < 15, AHI_ROWS[1],
                      np.where(ahi < 30, AHI_ROWS[2], AHI_ROWS[3]))), index=b.index)
    tst = b.TST_min
    dur = pd.Series(np.where(tst < 300, DUR_ROWS[0],
                    np.where((tst >= 360) & (tst <= 420), DUR_ROWS[1],
                    np.where(tst > 420, DUR_ROWS[2], "outside"))), index=b.index)
    control = crosstab(b[has_ahi], ox[has_ahi], apnea[has_ahi], AHI_ROWS)
    new = crosstab(b, ox, dur, DUR_ROWS)
    extra = dict(
        tst_300_to_360=int(((tst >= 300) & (tst < 360)).sum()),
        tst_missing=int(tst.isna().sum()),
        shown=int(sum(new[r]["total"] for r in DUR_ROWS)),
        cohort=len(b),
        ahi_undefined=int((~has_ahi).sum()),
    )
    assert extra["shown"] + extra["tst_300_to_360"] + extra["tst_missing"] == len(b), extra
    return control, new, extra


# ---- the control: every cell of Figure 4b as it is printed on the current sheet -------------
# v7 (2026-09-07): the control is the regenerated cross-tab numbers/crosstab_v2.json (an independent computation of the
# same apnea-by-T90 grid), not the strings of the pre-correction sheet.
import json as _json
_XT = _json.load(open(os.path.join(paths.NUMBERS_DIR, "crosstab_v2.json")))["rows"]
_BN = ["0-1%", "1-5%", "5-10%", ">10%"]
CANON_4B = {r: [f"{int(x[f'n_{bn}']):,} ({x[f'pct_{bn}']:.1f})" for bn in _BN] for r, x in zip(AHI_ROWS, _XT)}
CANON_4B_SHADE = {r: [band_index(x[f"pct_{bn}"]) for bn in _BN] for r, x in zip(AHI_ROWS, _XT)}


def control_gate(control):
    bad = []
    for r in AHI_ROWS:
        got = [c["label"] for c in control[r]["cells"]]
        if got != CANON_4B[r]:
            bad.append((r, got, CANON_4B[r]))
        gs = [c["shade"] for c in control[r]["cells"]]
        if gs != CANON_4B_SHADE[r]:
            bad.append((r + " shade", gs, CANON_4B_SHADE[r]))
    return bad


if __name__ == "__main__":
    control, new, extra = build()
    bad = control_gate(control)
    print("CONTROL, Figure 4b recomputed from the cohort:")
    for r in AHI_ROWS:
        print("  %-9s n=%5d  %s" % (r, control[r]["total"],
                                    "  ".join(c["label"] for c in control[r]["cells"])))
    print("  mismatches:", bad if bad else "NONE, 16 of 16 cells and 16 of 16 fills reproduce")
    assert not bad, bad
    print("\nNEW, sleep duration x T90:")
    for r in DUR_ROWS:
        print("  %-6s n=%5d  %s" % (r, new[r]["total"],
                                    "  ".join(c["label"] for c in new[r]["cells"])))
    print("  outside the paper's duration groups:", extra)
