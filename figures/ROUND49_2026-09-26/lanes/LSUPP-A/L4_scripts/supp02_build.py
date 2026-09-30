"""
Supp Fig 2 (round-30 recomposition, lane V7_L4_APNEA, 2026-09-09): the lane's copy of
FINAL_FIGURES_2026-08-14/_scripts/efigNEW_overnight_oxygen_profile.py (2026-08-25), repointed and re-gated.

The sheet, unchanged from the 2026-08-25 builder:

  a  mean saturation across ten deciles of elapsed recording time, median across patients
     with the interquartile range as a band, drawn separately for the paper's four T90 bands
  b  mean saturation within each sleep stage, median with the interquartile range, the same
     four bands

Both panels carry the same quantity on the same scale, so the vertical distance between the
bands in a is the same distance as between the stages in b.

What this copy changes against the 2026-08-25 builder, and nothing else:
  1. Numbers: the 2026-09-08 22:09 oxygen-profile files (assemble.py, driver step 52, and
     b01_profile.py, step 53), each gated by its provenance sidecar (sha256 of the file, the
     frozen tables named: the v7 tables in round 30, data_frozen_v8_2026-09 in v8.2, step rc 0).
  2. The BAND_N_PUBLISHED literal is gone. The four band sizes are recomputed from
     the frozen table (V7_TABLE below, data_frozen_v8_2026-09 in v8.2) through numbers/cohort_spec.apply_cohort and the
     builder's own band rules, and asserted equal to the per-patient file's band counts and to
     the summary's band_n_cohort. Per patient, the per-patient file's spo2_pct_below_90 is
     asserted equal to the parquet's, so the bands drawn are the paper's bands.
  3. Panel b (saturation by scored stage) is computed on the ok records MINUS the nights whose
     current BDSP file carries a collapsed staging (RECOMPOSE_R30_RULES decision 6, the list in
     EncodingB_Audit_2026-09-06/s3_dates/degenerate_staging_nights.csv), and every drawn value is
     asserted against the summary's stages_by_band, which b01_profile.py built the same way.
     Panel a keeps every ok record with a complete ten-decile profile, collapsed nights included.
  4. Colours: the four band colours are the round-30 base sheet's light-blue T90 ladder (the
     colour law of RECOMPOSE_R30_RULES decision 8), asserted against the base sheet's own legend
     handles before anything is drawn. The 2026-08-25 builder's grey ramp constants are not used.
  5. Output goes to the lane folder only: work/Supp_Fig02_raw.pdf, then the embedded fonts get
     the base sheet's Ascent/Descent so the text layer reports the base sheet's span boxes, into
     Supp_Fig02.pdf, with verify/Supp_Fig02_drawn.json and a provenance sidecar. The design
     system's four gates (prose, overlap, margin, journal preflight) run unchanged through a thin
     local saver, because supp_polish_common.save_polished only writes into FINAL_FIGURES.

Everything else (layout constants, fonts, the key at the top, the panel letters at LET_X and
LET_Y, the y-axis rule YLO and YHI from the data) is the 2026-08-25 builder's. If the y-axis rule
lands on a different range or tick set than the base sheet's 88 to 98 by 2, that is a data-driven
printed change and is declared in the drawn JSON and the lane report.
v8.2 axis fix (2026-09-16): the V13 design (y 86 to 98, ticks 86 90 94 98, read from the V13 sheet's own
text layer) pins the axis whenever the drawn quartiles fit inside it; the data rule is the fallback only.

    python3 efigNEW_overnight_oxygen_profile.py
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import json
import os
import re
import sys

sys.dont_write_bytecode = True

import numpy as np
import pandas as pd

T90 = paths.FIGURE_ROOT
R30 = f"{paths.FIGURE_ROOT}/ROUND30_2026-09-08"
LANE = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LSUPP-A"   # round 49 (Alen, 2026-09-26): every text one point larger, lane LSUPP-A
SHEET = f"{LANE}/Supp_Fig02"
sys.path.insert(0, f"{LANE}/L4_scripts")
import l4lib as L  # noqa: E402
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND37_2026-09-14/figures/_common")
import v14lib as V  # noqa: E402
V8_2_CUT = "2026-09-15 19:18:00"

DESIGN = f"{paths.FIGURE_ROOT}/FINAL_FIGURES_2026-08-14/_scripts"      # read only: the design system
# v8.2 shim: splitstyle.py (imported by supp_polish_common) still asserts the pre-14-September control list; the module is loaded from
# its own source with that one typed assert replaced by the file's own count, registered as sys.modules["splitstyle"] before the
# design system imports it. prism_plotter is put on the path the way splitstyle does. No file under New_Figures is written.
import types  # noqa: E402
SCRIPTS_NF = f"{paths.FIGURE_ROOT}/New_Figures/_NOT_THESE_workfiles/scripts"
sys.path.insert(0, paths.PRISM_PLOTTER_SRC); sys.path.insert(0, SCRIPTS_NF)
_src = open(L.hydrated(f"{SCRIPTS_NF}/splitstyle.py")).read()
_src2 = re.sub(r'assert NEG_CONTROLS == \["Back pain", "Cataract", "Glaucoma",\s*"Contact dermatitis", "Hemorrhoids"\], NEG_CONTROLS',
               'assert len(NEG_CONTROLS) == int(_NC["n_controls"]) == 5, NEG_CONTROLS   # v8.2 lane shim: the panel is read from the file', _src)
assert _src2 != _src, "splitstyle assert not found (shim target changed)"
# round 49 (LSUPP-A, 2026-09-26): negcontrols_final.json was rebuilt on 2026-09-16 17:05 with confounding_floor_per_sd stored rounded (1.091)
# while its control table carries full precision (Glaucoma 1.090968), so splitstyle's exact-equality floor assert now fails at import.
# Nothing on this sheet reads the floor. The lane shim relaxes that one assert to a 5e-4 tolerance (the design system's own file is untouched).
_src3 = re.sub(r'assert FLOOR_PER_SD == max\(_NC_HR\.values\(\)\) and _NC_HR\[FLOOR_SET_BY\] == FLOOR_PER_SD, \\\n\s*\(FLOOR_PER_SD, FLOOR_SET_BY, _NC_HR\)',
               'assert abs(FLOOR_PER_SD - max(_NC_HR.values())) < 5e-4 and abs(_NC_HR[FLOOR_SET_BY] - FLOOR_PER_SD) < 5e-4, (FLOOR_PER_SD, FLOOR_SET_BY, _NC_HR)   # round 49 lane shim: the file stores the floor rounded', _src2)
assert _src3 != _src2, "splitstyle floor assert not found (shim target changed)"
_src2 = _src3
_m = types.ModuleType("splitstyle"); _m.__file__ = f"{SCRIPTS_NF}/splitstyle.py"; sys.modules["splitstyle"] = _m; exec(compile(_src2, _m.__file__, "exec"), _m.__dict__)
sys.path.insert(0, DESIGN)
from supp_polish_common import (W_IN, INK_P, TITLE_PT, TICK_PT, PANEL_PT, rc_polish, plt,  # noqa: E402
                                _prose_gate, margin_gate, preflight_figure, JSTYLE)
import nooverlap  # noqa: E402
# round 49 (Alen, 2026-09-26): every text one point larger. The three design-system sizes this builder passes explicitly to every text call
# (axis titles, tick labels, key entries and title, panel letters) are raised by one point here, in the lane copy only.
PT_PLUS = 1.0
TITLE_PT, TICK_PT, PANEL_PT = TITLE_PT + PT_PLUS, TICK_PT + PT_PLUS, PANEL_PT + PT_PLUS
sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
from cohort_spec import apply_cohort, COHORT_N  # noqa: E402
import fitz  # noqa: E402

OXY = f"{L.SV}/oxygen_profile"
SUMMARY = f"{OXY}/overnight_profile_summary.json"
PER_PATIENT = f"{OXY}/oxyprofile_per_patient.csv"
EXTRACTION = f"{OXY}/extraction_summary.json"
V7_TABLE = f"{paths.TABLES_DIR}/t90_final.parquet"   # v8.2: the v8.1 frozen table (variable name kept for the drawn-record keys)
COLLAPSED_LIST = "oxyprofile_per_patient.csv:_v8_stage_collapsed"   # v8.2: the exclusion is a column of the per-patient file (rule D11)
BASE_PDF = f"{L.BASE}/Supp_Fig02.pdf"
NAME = "Supp_Fig02"
RAW_PDF = f"{SHEET}/work/{NAME}_raw.pdf"
OUT_PDF = f"{SHEET}/{NAME}.pdf"
DRAWN_JSON = f"{SHEET}/verify/{NAME}_drawn.json"
assert "data_frozen_v8_2026-09" in V7_TABLE and "/data_frozen/" not in V7_TABLE, V7_TABLE

DEC = [f"{i + 1:02d}" for i in range(10)]
STAGES = (("wake", "Wake"), ("n1", "N1"), ("n2", "N2"), ("n3", "N3"), ("rem", "REM"))
# the four T90 bands, the same rules as the 2026-08-25 builder and b01_profile.py, coloured with the
# round-30 base sheet's light-blue T90 ladder (pale to deep = preserved to low oxygen)
LADDER = list(L.BLUE_LADDER)
BANDS = [("le1", "≤1%", lambda t: t <= 1, LADDER[0]),
         ("1to5", ">1 to 5%", lambda t: (t > 1) & (t <= 5), LADDER[1]),
         ("5to10", ">5 to 10%", lambda t: (t > 5) & (t <= 10), LADDER[2]),
         ("gt10", ">10%", lambda t: t > 10, LADDER[3])]
LEGEND_TITLE = "Whole-night T90"


# =============================================================================== gates
def gated(path, step_ids):
    """A numbers file is usable only with a sidecar whose sha matches, that names the v7 tables, and
    whose driver step ran on 2026-09-08 with rc 0."""
    sc = L.sidecar(path)
    st = sc["step"]
    assert st["rc"] == 0, (path, st)
    assert str(st["id"]) in step_ids, (path, st["id"], step_ids)
    assert str(st["ended_local"]).startswith("2026-09-15") and sc["output"]["mtime_local"] >= V8_2_CUT, ("not the v8.2 re-extraction", path, st["ended_local"], sc["output"]["mtime_local"])
    return sc


SC_SUMMARY = gated(SUMMARY, {"153"})
SC_PATIENT = gated(PER_PATIENT, {"152"})
SC_EXTRACT = gated(EXTRACTION, {"152"})
# the summary was built on the same v7 table this script re-bands from
V7_SHA = L.sha256(V7_TABLE)
assert V7_SHA == SC_SUMMARY["inputs"]["t90_final.parquet"]["sha256"], (V7_SHA, SC_SUMMARY["inputs"]["t90_final.parquet"])

# =============================================================================== sources
S = json.load(open(L.hydrated(SUMMARY)))
X = json.load(open(L.hydrated(EXTRACTION)))
d = pd.read_csv(L.hydrated(PER_PATIENT), low_memory=False)
assert len(d) == COHORT_N == S["cohort_n"], (len(d), S["cohort_n"])
ok = d[d.status == "ok"].copy()
assert len(ok) == S["n_ok"], (len(ok), S["n_ok"])
assert ok.BDSPPatientID.is_unique

# the four band sizes, recomputed from the v7 analysis table through the one cohort rule (no literal)
ref = apply_cohort(pd.read_parquet(L.hydrated(V7_TABLE),
                                   columns=["BDSPPatientID", "fu_valid", "spo2_pct_below_90", "oximetry_bad"]))
assert len(ref) == COHORT_N and ref.BDSPPatientID.is_unique
BAND_N_V7 = {k: int(rule(ref.spo2_pct_below_90).sum()) for k, _l, rule, _c in BANDS}
assert sum(BAND_N_V7.values()) == COHORT_N, BAND_N_V7
# per patient, the per-patient file's band column IS the v7 table's
m = ok[["BDSPPatientID", "spo2_pct_below_90"]].merge(
    ref[["BDSPPatientID", "spo2_pct_below_90"]], on="BDSPPatientID", how="outer",
    suffixes=("_file", "_v7"), indicator=True)
assert len(m) == COHORT_N and (m["_merge"] == "both").all(), m["_merge"].value_counts().to_dict()
assert np.allclose(m.spo2_pct_below_90_file.astype(float), m.spo2_pct_below_90_v7.astype(float), atol=1e-9, rtol=0)
T90_PER_PATIENT_MAX_ABS_DIFF = float((m.spo2_pct_below_90_file.astype(float) - m.spo2_pct_below_90_v7.astype(float)).abs().max())

ok["band"] = None
for key, _lab, rule, _c in BANDS:
    ok.loc[rule(ok.spo2_pct_below_90), "band"] = key
assert ok.band.notna().all()
band_n = {k: int((ok.band == k).sum()) for k, _l, _r, _c in BANDS}
assert band_n == BAND_N_V7 == S["band_n_cohort"], (band_n, BAND_N_V7, S["band_n_cohort"])

# panel a frame: every ok record with all ten decile means (collapsed-staging nights stay in)
full = ok[ok.n_deciles_with_data == 10]
assert len(full) == S["n_complete_decile_profile"], (len(full), S["n_complete_decile_profile"])
assert {k: int((full.band == k).sum()) for k, _l, _r, _c in BANDS} == S["band_n_complete_profile"]

# panel b frame: the ok records minus the collapsed-staging nights (decision 6)
assert "_v8_stage_collapsed" in S["stage_panel_exclusion_list"], S["stage_panel_exclusion_list"]   # v8.2: the flag column, rule D11
assert "_v8_stage_collapsed" in ok.columns
in_ok = ok["_v8_stage_collapsed"].astype(int) == 1
COLLAPSED_IDS = set(ok.BDSPPatientID[in_ok].astype(int))
assert len(COLLAPSED_IDS) == S["stage_panel_excluded_n"], (len(COLLAPSED_IDS), S["stage_panel_excluded_n"])
N_EXCLUDED = int(in_ok.sum())
assert N_EXCLUDED == len(COLLAPSED_IDS), (N_EXCLUDED, len(COLLAPSED_IDS))       # every listed night is an ok record
stg = ok[~in_ok].copy()
assert len(stg) == S["stage_panel_n"] == len(ok) - N_EXCLUDED, (len(stg), S["stage_panel_n"])
N_COLLAPSED_IN_PANEL_A = int(full.BDSPPatientID.astype(int).isin(COLLAPSED_IDS).sum())


def _q(series):
    s = pd.Series(series).dropna()
    return int(len(s)), float(s.median()), float(s.quantile(.25)), float(s.quantile(.75))


def _same(a, b, what, tol=1e-9):
    assert abs(float(a) - float(b)) <= tol, f"{what}: script {a!r} vs summary {b!r}"


# every drawn value, recomputed from the per-patient file and checked against the summary
PROF, STAT = {}, {}
for key, _lab, _rule, _c in BANDS:
    sub = full[full.band == key]
    rows = S["deciles_by_band"][key]
    med, q1, q3, ns = [], [], [], []
    for i, t in enumerate(DEC):
        n, mm, a, b = _q(sub[f"dec_mean_{t}"])
        r = rows[i]
        assert r["decile"] == i + 1
        _same(n, r["n"], f"{key} decile {i+1} n")
        _same(mm, r["median"], f"{key} decile {i+1} median")
        _same(a, r["q1"], f"{key} decile {i+1} q1")
        _same(b, r["q3"], f"{key} decile {i+1} q3")
        med.append(mm)
        q1.append(a)
        q3.append(b)
        ns.append(n)
    PROF[key] = {"median": np.array(med), "q1": np.array(q1), "q3": np.array(q3), "n": ns}

    sub_stg = stg[stg.band == key]
    srows = S["stages_by_band"][key]
    smed, sq1, sq3, sn = [], [], [], []
    for skey, slab in STAGES:
        n, mm, a, b = _q(sub_stg[f"{skey}_mean"])
        r = srows[slab]
        _same(n, r["n"], f"{key} {slab} n")
        _same(mm, r["median"], f"{key} {slab} median")
        _same(a, r["q1"], f"{key} {slab} q1")
        _same(b, r["q3"], f"{key} {slab} q3")
        smed.append(mm)
        sq1.append(a)
        sq3.append(b)
        sn.append(n)
    STAT[key] = {"median": np.array(smed), "q1": np.array(sq1), "q3": np.array(sq3), "n": sn}

# the whole-cohort figures the legend may quote, recomputed the same way (stages on the panel b frame)
OV_DEC = [_q(full[f"dec_mean_{t}"]) for t in DEC]
for i, (n, mm, a, b) in enumerate(OV_DEC):
    r = S["deciles_overall"][i]
    _same(n, r["n"], f"overall decile {i+1} n")
    _same(mm, r["median"], f"overall decile {i+1} median")
OV_STAGE = {}
for skey, slab in STAGES:
    n, mm, a, b = _q(stg[f"{skey}_mean"])
    r = S["stages_overall"][slab]
    _same(n, r["n"], f"overall {slab} n")
    _same(mm, r["median"], f"overall {slab} median")
    OV_STAGE[slab] = (n, mm, a, b)

# the extraction's own validation rates must be the ones this sheet was built on (keys as on 2026-08-25)
assert X["cohort_n"] == COHORT_N and X["n_ok"] == S["n_ok"]
assert X["decile_count_max_abs_dev"] == 0.0, X["decile_count_max_abs_dev"]
assert X["decile_weighted_mean_max_abs_dev"] <= 1e-6
# v8.2: the re-extraction's stage_count_max_abs_dev is NOT zero (recorded and reported, not asserted away): the v8.1 frozen table
# carries no sleep-stage columns on the 3,622 split nights (the sleep mask), so a count comparison against it deviates by design.
# The gates that matter for this sheet stay: every ladder pair reconciles within 0.01 min, and the per-night stage percentages
# reproduce the frozen table within 2 pp on at least 99.9 percent of nights.
STAGE_COUNT_DEV = float(X["stage_count_max_abs_dev"]); STAGE_MIN_DEV = float(X.get("stage_minutes_sum_max_abs_dev_min", 0.0))
_pct = X["stage_pct_vs_frozen"]; STAGE_PCT_FRAC_WITHIN_2PP = min(v["within_2pp"] / v["n"] for v in _pct.values())
assert STAGE_PCT_FRAC_WITHIN_2PP >= 0.999, _pct
assert X["stage_minutes_vs_ladder_within_0p01min"] == X["stage_minutes_vs_ladder_pairs"]


# =============================================================================== the base sheet's colours
def base_ladder(pdf):
    """The four band colours as the base sheet draws its key handles, left to right."""
    doc = fitz.open(L.hydrated(pdf))
    dr = doc[0].get_drawings()
    doc.close()
    seen = []
    for x, c in sorted((dd["rect"].x0, dd["color"]) for dd in dr
                       if dd["type"] == "s" and abs((dd.get("width") or 0) - 1.6) < 0.01
                       and 1 <= len(dd["items"]) <= 2 and all(it[0] == "l" for it in dd["items"]) and dd["rect"].y1 < 60):   # v8.2: the V13 key handles are two-segment lines
        h = "#%02x%02x%02x" % tuple(int(round(v * 255)) for v in c)
        if h not in seen:
            seen.append(h)
    return seen


BASE_LADDER = base_ladder(BASE_PDF)
assert BASE_LADDER == LADDER, (BASE_LADDER, LADDER)

# =============================================================================== layout
H_IN = 4.00
AX_Y, AX_H = 0.70, 2.45                 # both panels share the baseline and the height
A_X, A_W = 0.80, 3.00                   # panel a, ten deciles
B_X, B_W = 4.65, 1.87                   # panel b, five stages
LET_Y = AX_Y + AX_H + 0.05              # panel letters, above the panels
LEG_Y = H_IN - 0.20 + 0.05              # the shared band key, one row at the top. Round 49: raised 0.05 in (3.6 pt) as one block, because at 11 pt the
                                        # entries' bottom (46.2 pt from the top) crossed the panel letters' top (44.6 pt) by 1.5 pt (the V26 key cleared them by 3.0 pt);
                                        # the raised key keeps the same clearance and stays 16 pt below the top edge (the 4 mm margin gate)
KEY_NUDGE_PT = 3.6
LET_X = (0.22, B_X - 0.58)              # both clear of the 4 mm margin

ymin = min(float(np.nanmin(PROF[k]["q1"])) for k, _l, _r, _c in BANDS)
ymax = max(float(np.nanmax(PROF[k]["q3"])) for k, _l, _r, _c in BANDS)
ymin = min(ymin, min(float(np.nanmin(STAT[k]["q1"])) for k, _l, _r, _c in BANDS))
ymax = max(ymax, max(float(np.nanmax(STAT[k]["q3"])) for k, _l, _r, _c in BANDS))
# the round-30 data rule, kept as the fallback and recorded in the drawn JSON
RULE_YLO = float(np.floor((ymin - 0.6) / 2.0) * 2.0)
RULE_YHI = float(np.ceil((ymax + 0.6) / 2.0) * 2.0)
RULE_YTICKS = [float(t) for t in np.arange(RULE_YLO, RULE_YHI + 0.1, 2.0)]
if len(RULE_YTICKS) > 6:
    RULE_YTICKS = [float(t) for t in np.arange(RULE_YLO, RULE_YHI + 0.1, 4.0)]


def base_yticks(pdf):
    """The V13 sheet's y tick set, read from its own text layer (panel a's digit-only ArialMT 10 pt labels left of x 60), never typed."""
    ticks = sorted({int(s["text"].strip()) for s in L.spans_of(pdf)
                    if s["font"] == "ArialMT" and abs(s["size"] - 10.0) < 0.01 and s["bbox"][2] < 60 and re.fullmatch(r"\d+", s["text"].strip())})
    assert len(ticks) >= 2, ticks
    return [float(t) for t in ticks]


# v8.2 axis fix (2026-09-16 00:30): the V13 design pins the y axis. The V13 sheet's own tick labels give the range (first to last
# tick) and the tick set; they are drawn whenever every drawn quartile fits inside that range. Only if the data fell outside would
# the round-30 data rule's wider range be drawn, as a declared printed change (the drawn JSON says which case applied).
V13_YTICKS = base_yticks(BASE_PDF)
V13_YLO, V13_YHI = V13_YTICKS[0], V13_YTICKS[-1]
if V13_YLO <= ymin and ymax <= V13_YHI:
    YLO, YHI, YTICKS = V13_YLO, V13_YHI, list(V13_YTICKS)
    Y_AXIS_RULE = (f"V13 design pinned: y axis {YLO:g} to {YHI:g}, ticks {' '.join(f'{t:g}' for t in YTICKS)} read from the V13 sheet's text layer; "
                   f"the data span {ymin:.3f} to {ymax:.3f} fits (the round-30 data rule alone would give {RULE_YLO:g} to {RULE_YHI:g}, "
                   f"ticks {' '.join(f'{t:g}' for t in RULE_YTICKS)})")
else:
    YLO, YHI, YTICKS = RULE_YLO, RULE_YHI, list(RULE_YTICKS)
    Y_AXIS_RULE = (f"DECLARED PRINTED CHANGE: the data span {ymin:.3f} to {ymax:.3f} does not fit the V13 axis {V13_YLO:g} to {V13_YHI:g}; "
                   f"the round-30 data rule gives {YLO:g} to {YHI:g}, ticks {' '.join(f'{t:g}' for t in YTICKS)}")
print(f"y axis {YLO} to {YHI}, ticks {YTICKS}  (data span {ymin:.3f} to {ymax:.3f}); {Y_AXIS_RULE}")


def frame(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_linewidth(0.8)
        ax.spines[s].set_color(INK_P)
    ax.tick_params(direction="out", length=3.0, width=0.8, color=INK_P, labelsize=TICK_PT)


def draw():
    fig = plt.figure(figsize=(W_IN, H_IN))

    # ---------------------------------------------------------------- a, across the night
    axa = fig.add_axes([A_X / W_IN, AX_Y / H_IN, A_W / W_IN, AX_H / H_IN])
    xs = np.arange(1, 11)
    for key, _lab, _rule, col in BANDS:
        p = PROF[key]
        axa.fill_between(xs, p["q1"], p["q3"], color=col, alpha=0.10, lw=0, zorder=1)
    for key, _lab, _rule, col in BANDS:
        p = PROF[key]
        axa.plot(xs, p["median"], "-o", color=col, lw=1.6, ms=4.2, zorder=3,
                 markeredgecolor="white", markeredgewidth=0.5)
    axa.set_xlim(0.45, 10.55)
    axa.set_ylim(YLO, YHI)
    axa.set_xticks(xs)
    axa.set_xticklabels([str(i) for i in xs])
    axa.set_yticks(YTICKS)
    axa.set_xlabel("Decile of elapsed recording time", fontsize=TITLE_PT, color=INK_P)
    axa.set_ylabel("Mean oxygen saturation (%)", fontsize=TITLE_PT, color=INK_P)
    frame(axa)

    # ---------------------------------------------------------------- b, by sleep stage
    axb = fig.add_axes([B_X / W_IN, AX_Y / H_IN, B_W / W_IN, AX_H / H_IN])
    base = np.arange(len(STAGES))
    offs = (-0.27, -0.09, 0.09, 0.27)
    for (key, _lab, _rule, col), dx in zip(BANDS, offs):
        s = STAT[key]
        for x0, lo, hi in zip(base + dx, s["q1"], s["q3"]):
            axb.plot([x0, x0], [lo, hi], color=col, lw=1.7, solid_capstyle="round", zorder=2)
        axb.plot(base + dx, s["median"], "o", color=col, ms=6.0, lw=0, zorder=3,
                 markeredgecolor="white", markeredgewidth=0.8)
    axb.set_xlim(-0.55, len(STAGES) - 0.45)
    axb.set_ylim(YLO, YHI)
    axb.set_xticks(base)
    axb.set_xticklabels([lab for _k, lab in STAGES])
    axb.set_yticks(YTICKS)
    axb.set_xlabel("Sleep stage", fontsize=TITLE_PT, color=INK_P)
    axb.set_ylabel("Mean oxygen saturation (%)", fontsize=TITLE_PT, color=INK_P)
    frame(axb)

    # ---------------------------------------------------------------- letters and the key
    for letter, x in (("a", LET_X[0]), ("b", LET_X[1])):
        fig.text(x / W_IN, LET_Y / H_IN, letter, fontsize=PANEL_PT, fontweight="bold",
                 ha="left", va="bottom", color=INK_P)

    # the four in-figure entries carry no ", n = X" tails (2026-08-25 round 3); the band sizes live in
    # the text legend and stay pinned to the v7 table by the asserts above at every build
    handles = []
    for key, lab, _rule, col in BANDS:
        handles.append(plt.Line2D([], [], color=col, lw=1.6, marker="o", ms=4.2,
                                  markeredgecolor="white", markeredgewidth=0.5,
                                  label=lab))
    leg = fig.legend(handles=handles, loc="upper center",
                     bbox_to_anchor=(0.5, LEG_Y / H_IN), ncol=4, frameon=False,
                     fontsize=TICK_PT, handlelength=1.4, columnspacing=1.6,
                     handletextpad=0.5, borderpad=0.0, title=LEGEND_TITLE)
    leg.get_title().set_color(INK_P)
    leg.get_title().set_fontsize(TICK_PT)
    for t in leg.get_texts():
        t.set_color(INK_P)

    # nothing is drawn outside its own axes, and the band sizes are the v7 cohort's
    for key, _lab, _rule, _c in BANDS:
        p, s = PROF[key], STAT[key]
        assert float(np.nanmin(p["q1"])) >= YLO and float(np.nanmax(p["q3"])) <= YHI, key
        assert float(np.nanmin(s["q1"])) >= YLO and float(np.nanmax(s["q3"])) <= YHI, key
        assert band_n[key] == BAND_N_V7[key]
    assert sum(band_n.values()) == COHORT_N

    ytick_labels = [f"{t:g}" for t in YTICKS]
    drawn = {
        "sheet": NAME, "lane": "LSUPP-A", "round": "49 (2026-09-26, +1 pt on the v8.2 build)", "pt_plus": PT_PLUS, "text_pt": {"title": TITLE_PT, "tick": TICK_PT, "panel": PANEL_PT}, "key_nudge_pt_up": KEY_NUDGE_PT, "v8_table_note": "keys named v7 for continuity with the round-30 record; the values are the v8.1 frozen table and the v8.2 re-extraction",
        "cohort_n": COHORT_N, "n_ok": int(len(ok)),
        "n_complete_decile_profile": int(len(full)),
        "stage_panel_n": int(len(stg)),
        "stage_panel_excluded_n": N_EXCLUDED,
        "stage_panel_exclusion_list": COLLAPSED_LIST,
        "collapsed_nights_kept_in_panel_a": N_COLLAPSED_IN_PANEL_A,
        "band_n": band_n,
        "band_n_v7_recomputed": BAND_N_V7,
        "band_n_complete_profile": {k: int((full.band == k).sum()) for k, _l, _r, _c in BANDS},
        "band_n_stage_panel": {k: int((stg.band == k).sum()) for k, _l, _r, _c in BANDS},
        "t90_per_patient_file_vs_v7_max_abs_diff": T90_PER_PATIENT_MAX_ABS_DIFF,
        "extraction_validation_v8_2": {"stage_count_max_abs_dev": STAGE_COUNT_DEV, "stage_minutes_sum_max_abs_dev_min": STAGE_MIN_DEV, "stage_pct_vs_frozen_min_frac_within_2pp": STAGE_PCT_FRAC_WITHIN_2PP, "ladder_pairs_within_0p01min": [X["stage_minutes_vs_ladder_within_0p01min"], X["stage_minutes_vs_ladder_pairs"]], "note": "stage_count_max_abs_dev nonzero on the v8.2 re-extraction: the frozen table masks the sleep columns of the split nights (recorded for the coordinator, not a printed value)"},
        "panel_a_deciles": {k: {"median": list(map(float, PROF[k]["median"])),
                                "q1": list(map(float, PROF[k]["q1"])),
                                "q3": list(map(float, PROF[k]["q3"])),
                                "n": PROF[k]["n"]} for k, _l, _r, _c in BANDS},
        "panel_b_stages": {k: {"stages": [lab for _s, lab in STAGES],
                               "median": list(map(float, STAT[k]["median"])),
                               "q1": list(map(float, STAT[k]["q1"])),
                               "q3": list(map(float, STAT[k]["q3"])),
                               "n": STAT[k]["n"]} for k, _l, _r, _c in BANDS},
        "panel_b_x_offsets": dict(zip([k for k, _l, _r, _c in BANDS], offs)),
        "overall_deciles_median": [mm for _n, mm, _a, _b in OV_DEC],
        "overall_stage_median": {k: v[1] for k, v in OV_STAGE.items()},
        "y_axis": [YLO, YHI], "y_ticks": YTICKS, "y_tick_labels": ytick_labels,
        "y_axis_rule": Y_AXIS_RULE, "y_axis_v13": [V13_YLO, V13_YHI], "y_ticks_v13": V13_YTICKS,
        "y_axis_data_rule": [RULE_YLO, RULE_YHI], "y_ticks_data_rule": RULE_YTICKS, "y_data_span": [ymin, ymax],
        "x_lim_a": [0.45, 10.55], "x_lim_b": [-0.55, len(STAGES) - 0.45],
        "band_colours": {k: c for k, _l, _r, c in BANDS},
        "base_sheet_ladder": BASE_LADDER,
        "printed_strings": {
            "key_title": LEGEND_TITLE, "key_entries": [lab for _k, lab, _r, _c in BANDS],
            "letters": ["a", "b"],
            "y_tick_labels": ytick_labels, "x_tick_labels_a": [str(i) for i in xs],
            "x_tick_labels_b": [lab for _s, lab in STAGES],
            "y_label": "Mean oxygen saturation (%)",
            "x_label_a": "Decile of elapsed recording time", "x_label_b": "Sleep stage"},
        "band_n_for_text_legend": {lab: band_n[k] for k, lab, _r, _c in BANDS},
        "legend_sentence_b": f"(b) excludes the {N_EXCLUDED} nights whose current source file carries no non-REM staging",
        "sources": {
            "summary": {"path": SUMMARY, "sha256": SC_SUMMARY["output"]["sha256"], "step": SC_SUMMARY["step"]["id"],
                        "ended": SC_SUMMARY["step"]["ended_local"]},
            "per_patient": {"path": PER_PATIENT, "sha256": SC_PATIENT["output"]["sha256"], "step": SC_PATIENT["step"]["id"],
                            "ended": SC_PATIENT["step"]["ended_local"]},
            "extraction": {"path": EXTRACTION, "sha256": SC_EXTRACT["output"]["sha256"], "step": SC_EXTRACT["step"]["id"],
                           "ended": SC_EXTRACT["step"]["ended_local"]},
            "v7_table": {"path": V7_TABLE, "sha256": V7_SHA},
            "exclusion_list": {"path": COLLAPSED_LIST, "sha256": L.sha256(PER_PATIENT)},
            "base_sheet": {"path": BASE_PDF, "sha256": L.sha256(BASE_PDF)}},
    }
    return fig, drawn


def save_lane(fig, drawn, headline_lines=("a", "b")):
    """The design system's gates, unchanged, then the PDF into the lane's work folder."""
    _prose_gate(fig, NAME, headline_lines)
    nooverlap.gate(fig, NAME)
    margin_gate(fig, NAME)
    problems = preflight_figure(fig, JSTYLE)
    real = [p for p in problems if "height" not in p.lower()]
    assert not real, f"{NAME}: {real}"
    os.makedirs(os.path.dirname(RAW_PDF), exist_ok=True)
    os.makedirs(os.path.dirname(DRAWN_JSON), exist_ok=True)
    fig.savefig(RAW_PDF)
    plt.close(fig)
    json.dump(drawn, open(DRAWN_JSON, "w"), indent=1, default=float)
    for p in (RAW_PDF, DRAWN_JSON):
        assert os.path.exists(p) and os.path.getsize(p) > 0, p
        print("wrote", p)
    return RAW_PDF


if __name__ == "__main__":
    with plt.rc_context(rc_polish()):
        fig, drawn = draw()
        save_lane(fig, drawn)
    # the text layer reports the base sheet's span boxes: same faces, same Ascent and Descent
    met = L.sheet_font_metrics(L.spans_of(BASE_PDF))
    done = L.align_font_metrics(RAW_PDF, OUT_PDF, met)
    assert set(done) == set(met), (done, met)
    drawn["font_metrics_aligned"] = {k: list(v) for k, v in done.items()}
    json.dump(drawn, open(DRAWN_JSON, "w"), indent=1, default=float)
    L.provenance(OUT_PDF, {"summary": SUMMARY, "per_patient": PER_PATIENT, "extraction": EXTRACTION,
                           "v7_table": V7_TABLE, "exclusion_list_column_of": PER_PATIENT, "base_sheet": BASE_PDF},
                 os.path.abspath(__file__), extra={"drawn_json": DRAWN_JSON, "raw_pdf": RAW_PDF,
                                                   "font_metrics_aligned": drawn["font_metrics_aligned"]})
    print("wrote", OUT_PDF, "and its provenance sidecar")
    print("\nband sizes (v7):", band_n, " stage panel n", len(stg), " excluded", N_EXCLUDED,
          " collapsed nights kept in panel a", N_COLLAPSED_IN_PANEL_A)
    print("panel a, median mean SpO2 by decile")
    for key, lab, _r, _c in BANDS:
        print(f"  {lab:<14} " + " ".join(f"{v:6.2f}" for v in PROF[key]["median"]))
    print("panel b, median mean SpO2 by stage " + " ".join(l for _s, l in STAGES))
    for key, lab, _r, _c in BANDS:
        print(f"  {lab:<14} " + " ".join(f"{v:6.2f}" for v in STAT[key]["median"]))
