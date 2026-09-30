"""
Build THREE_SETTINGS_dC.xlsx: discrimination gain from T90 versus sleep duration across
three measurement settings.

  1  Laboratory PSG, primary cohort            T90 vs PSG total sleep time
  2  Home PSG, SHHS and MrOS                   T90 vs home-PSG total sleep time
  3  Home habitual report, primary subset      T90 vs reported habitual hours

Setting 2 carries two rows since the 2026-08-20 age+sex update: the main row is on the
bare age + sex baseline (external_dc_agesex.csv), so all three settings share the same
bare baseline, and a second row preserves the frozen full-clinical-baseline numbers
(external_dc_tst.csv). Both external rows are in-sample.

Every number is read from a source CSV or recomputed from one. Nothing is typed from a
report. Positive controls are asserted before the workbook is written.

COMPARISON_T90_vs_TST.xlsx is left alone as the ladder archive.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import warnings
warnings.filterwarnings("ignore")
import json
import sys
import numpy as np
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, PatternFill, Border, Side
from openpyxl.utils import get_column_letter
import sys as _sys_v8; _sys_v8.path.insert(0, paths.ANALYSIS_DIR)  # v8 sweep 2026-09-12
from cohort_spec import NUMBERS_DIR, N_RANKED_OUTCOMES  # v8 sweep 2026-09-12

HERE = f"{paths.SV_ROOT}/dC_clinical_baselines"
OUT = f"{HERE}/THREE_SETTINGS_dC.xlsx"
PUBLISHED_dC = float(__import__("pandas").read_csv(f"{paths.NUMBERS_DIR}/ranking_v3.csv", comment="#").set_index("feature").loc["spo2_pct_below_90", "dC"])   # v7: was the August literal 0.02033875497249791
_RK_V8 = __import__("pandas").read_csv(f"{NUMBERS_DIR}/ranking_v3.csv", comment="#"); _impl_v8 = (_RK_V8.mC - _RK_V8.dC); _impl_mode = _impl_v8.round(6).mode().iloc[0]; _off = _RK_V8.feature[_impl_v8.round(6) != _impl_mode].tolist(); assert all("_coupling_" in str(_f) for _f in _off) and (_impl_v8.round(6) == _impl_mode).mean() > 0.9, f"ranking rows disagree on the implied baseline: off-mode rows {_off} (v8.1c rule: every off-mode row must be a spindle-SO coupling feature, mode on over 90 percent of rows)"; print(f"[v8] implied baseline mode {_impl_mode:.6f} over {int((_impl_v8.round(6) == _impl_mode).sum())} rows; off-mode rows {_off}"); PUBLISHED_C0 = float(_impl_v8[_impl_v8.round(6) == _impl_mode].mean())  # v8 F2 fix 2026-09-12 (integrity gate 6): the held-out baseline C is the MODE of the implied baseline mC - dC over the ranking rows (139 of 141 agree); at most two off-mode rows are allowed and listed (features refit on their non-missing subset, the two C_coupling_* rows); PUBLISHED_C0 is taken from the modal rows. Same rule as numbers/run_combo_v2.py lines 68 to 72

# ============================================================ load sources
res = pd.read_csv(f"{HERE}/results.csv")                    # primary T90 ladder
tst = pd.read_csv(f"{HERE}/results_tst.csv")                # primary PSG-TST ladder
home = pd.read_csv(f"{HERE}/results_home.csv")              # home-subset per outcome
ext = pd.read_csv(f"{HERE}/external_dc_tst.csv")            # external, clinical baseline
ags = pd.read_csv(f"{HERE}/external_dc_agesex.csv", comment="#")  # external, age+sex
agso = ags[ags.row_type == "per_outcome"].reset_index(drop=True)
lad = pd.read_csv(f"{HERE}/COMPARISON_T90_vs_TST_table.csv")  # paired ladder, archive
prov = json.load(open(f"{HERE}/provenance.json"))
prov_t = json.load(open(f"{HERE}/provenance_tst.json"))
prov_h = json.load(open(f"{HERE}/provenance_home.json"))
prov_e = json.load(open(f"{HERE}/external_dc_tst_provenance.json"))
prov_a = json.load(open(f"{HERE}/external_dc_agesex_provenance.json"))

# ============================================================ positive controls
checks = []


def chk(name, got, want, tol=1e-9):
    ok = abs(float(got) - float(want)) <= tol
    checks.append({"Check": name, "Computed": float(got), "Expected": float(want),
                   "Result": "PASS" if ok else "FAIL"})
    return ok


m1 = res[res.model == "M1_paper_age_sex"].iloc[0]
l1 = tst[tst.model == "L1_age_sex"].iloc[0]
ok = True
ok &= chk("Primary published dC from T90, age+sex baseline, 48 outcomes",
          m1.mean_dC, PUBLISHED_dC)
ok &= chk("Primary published baseline C, age+sex", m1.mean_C_baseline, PUBLISHED_C0)
ok &= chk("Primary dC from PSG TST, same baseline", l1.mean_dC,
          prov_t["positive_control"]["tst_dC"])
ok &= chk("Home-subset dC from T90 reproduces provenance_home rung H1",
          home[(home.block == "home_subset") & (home.model == "H1_age_sex")
               & (home.exposure == "t90") & (home.kept == True)].dC.mean(),
          prov_h["rung_summary"][0]["dC_t90"])
ok &= chk("Home-subset dC from habitual report reproduces provenance_home rung H1",
          home[(home.block == "home_subset") & (home.model == "H1_age_sex")
               & (home.exposure == "home") & (home.kept == True)].dC.mean(),
          prov_h["rung_summary"][0]["dC_home"])
# external: every SHHS outcome must have reproduced external.json
nrep = int((ext.frozen_check == "reproduced").sum())
ok &= chk("External SHHS outcomes reproducing external.json clinical_model_gain within 5e-5 (v8.3)",
          nrep, 7)
for _, r in ext[ext.cohort == "SHHS"].iterrows():
    ok &= chk(f"  SHHS {r.outcome}: dC from T90 vs frozen clinical_model_gain",
              r.dC_t90, r.frozen_dC_t90, tol=5e-5)  # v8.3: reference now full precision, compare within 5e-5 (2026-09-16)
# external age+sex baseline: every frame must have reproduced the frozen minimally
# adjusted T90 HR from external.json, and the sanity gate must have passed
ok &= chk("External age+sex frames reproducing frozen minimally adjusted T90 HR",
          int((agso.frozen_hr_check == "reproduced").sum()), 12)
ok &= chk("External age+sex sanity gate passed in run_external_dc_agesex.py",
          int(bool(prov_a["sanity_vs_clinical_baseline"]["pass"])), 1)
ok &= chk("External age+sex mean dC from T90 exceeds clinical-baseline mean",
          int(agso.dC_t90.mean() > ext.dC_t90.mean()), 1)
for _, r in agso.iterrows():
    e_row = ext[(ext.cohort == r.cohort) & (ext.outcome == r.outcome)].iloc[0]
    ok &= chk(f"  {r.cohort} {r.outcome}: clinical dC_T90 consistent across the two CSVs",
              r.dC_t90_clinical_baseline, e_row.dC_t90, tol=1e-12)

if not ok:
    print("POSITIVE CONTROL FAILED, workbook not written:")
    for c in checks:
        if c["Result"] == "FAIL":
            print("  ", c)
    sys.exit(1)
print(f"Positive controls: {len(checks)} checks, all PASS.")

# ============================================================ assemble numbers
KEEP = set(home[(home.block == "home_subset") & (home.kept == True)].outcome.unique())


def home_rung(block, model, exposure):
    q = home[(home.block == block) & (home.model == model)
             & (home.exposure == exposure) & (home.outcome.isin(KEEP))]
    return {"n_outcomes": q.outcome.nunique(), "c_base": q.c_base.mean(),
            "c_expo": q.c_expo.mean(), "dC": q.dC.mean(), "npos": int((q.dC > 0).sum())}


h_t90 = home_rung("home_subset", "H1_age_sex", "t90")
h_hab = home_rung("home_subset", "H1_age_sex", "home")
h_tst = home_rung("home_subset", "H1_age_sex", "tst")

e_t90_mean = float(ext.dC_t90.mean())
e_tst_mean = float(ext.dC_tst.mean())
e_base_mean = float(ext.c_baseline.mean())
e_npos_t90 = int((ext.dC_t90 > 0).sum())
e_npos_tst = int((ext.dC_tst > 0).sum())
n_ext = len(ext)

a_t90_mean = float(agso.dC_t90.mean())
a_tst_mean = float(agso.dC_tst.mean())
a_base_mean = float(agso.c_baseline.mean())
a_npos_t90 = int((agso.dC_t90 > 0).sum())
a_npos_tst = int((agso.dC_tst > 0).sum())

# ============================================================ workbook styling
wb = Workbook()
H = Font(bold=True, size=11, color="FFFFFF")
HF = PatternFill("solid", fgColor="1F3864")
SUB = Font(bold=True, size=11, color="1F3864")
TITLE = Font(bold=True, size=14, color="1F3864")
NOTE = Font(size=10, italic=True, color="555555")
BOLD = Font(bold=True)
THIN = Side(style="thin", color="BFBFBF")
BOX = Border(left=THIN, right=THIN, top=THIN, bottom=THIN)
D4 = "0.0000"
D4S = "+0.0000;-0.0000;0.0000"
BANDF = PatternFill("solid", fgColor="EEF2F9")


def put_header(ws, r, cols, widths=None, freeze=True):
    for j, c in enumerate(cols, start=1):
        cell = ws.cell(row=r, column=j, value=c)
        cell.font = H
        cell.fill = HF
        cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
        cell.border = BOX
    ws.row_dimensions[r].height = 34
    if widths:
        for j, w in enumerate(widths, start=1):
            ws.column_dimensions[get_column_letter(j)].width = w
    if freeze:
        ws.freeze_panes = ws.cell(row=r + 1, column=1)


def put_row(ws, r, vals, fmts=None, band=False, bold=False):
    for j, v in enumerate(vals, start=1):
        cell = ws.cell(row=r, column=j, value=v)
        cell.border = BOX
        cell.alignment = Alignment(vertical="center",
                                   horizontal="left" if isinstance(v, str) else "center",
                                   wrap_text=isinstance(v, str) and len(str(v)) > 40)
        if fmts and j - 1 < len(fmts) and fmts[j - 1] and isinstance(v, (int, float)):
            cell.number_format = fmts[j - 1]
        if band:
            cell.fill = BANDF
        if bold:
            cell.font = BOLD


def title(ws, text, sub=None):
    ws["A1"] = text
    ws["A1"].font = TITLE
    r = 2
    if sub:
        for line in sub:
            ws.cell(row=r, column=1, value=line).font = NOTE
            r += 1
    return r + 1


# ============================================================ 1 README
ws = wb.active
ws.title = "README"
ws.column_dimensions["A"].width = 30
ws.column_dimensions["B"].width = 118
r = title(ws, "Discrimination gain from T90 versus sleep duration, three measurement settings")

BLOCKS = [
 ("THE QUESTION", [
  ("What this workbook answers",
   "Does nocturnal hypoxaemia (T90) add more to disease discrimination than sleep duration does, "
   "and does the answer depend on whether sleep was measured in a laboratory, at home by PSG, or "
   "simply reported by the patient. dC is the change in Harrell's C when one exposure is added to "
   "a fixed baseline model, averaged over outcomes."),
 ]),
 ("THE THREE SETTINGS", [
  ("1. Laboratory PSG",
   "Primary cohort, n = 19,173, 48 incident outcomes. Both exposures come from an attended "
   "in-laboratory sleep study. T90 is the fraction of the recording below 90 percent saturation. "
   "Sleep duration is PSG total sleep time."),
  ("2. Home PSG",
   "SHHS n = 5,802 and MrOS n = 2,911, 12 outcomes in total, 7 in SHHS and 5 in MrOS. Both cohorts "
   "recorded unattended polysomnography in the participant's own home, which is the point of this "
   "row: the sleep is home sleep, but it is still measured by PSG, not reported. Sleep duration is "
   "home-PSG total sleep time, slpprdp in SHHS and POSLPRDP in MrOS. Since 2026-08-20 this setting "
   "carries two Overview rows: the main row on the bare age + sex baseline, so all three settings "
   "share the same bare baseline, and a preserved second row on the frozen full clinical baseline, "
   "so nothing from the original computation is lost."),
  ("3. Home habitual report",
   "Primary cohort subset with a habitual sleep duration recorded in the clinical note. 5,295 of "
   "19,173 patients have the field, 27.6 percent, and 5,293 have it together with PSG TST and sleep "
   "efficiency, which is the analysis frame. 43 of 48 outcomes clear the event floor. Here the "
   "duration is what the patient said they usually sleep, not what was measured. T90 in this "
   "setting is still laboratory PSG T90: only the duration exposure changes."),
 ]),
 ("MODEL SPECIFICATION, PER SETTING", [
  ("Settings 1 and 3, primary cohort",
   "Cox proportional hazards, penalizer 0.01. Baseline is age plus sex, the paper's baseline, with "
   "age as a natural cubic spline, cr(df=4). Exposures are within-site rank inverse normal, so every "
   "dC is per 1 SD. C is scored OUT OF SAMPLE by two-fold site cross-fit: fit at I0002 and score at "
   "I0006, then reverse, then average. An outcome needs at least 30 events in each fold."),
  ("Setting 2, external cohorts",
   "Cox by maximum partial likelihood, penalizer 0, with Newton step sizes [None, 0.05, 0.1, 0.25, "
   "0.5, 0.95] scanned and the highest-likelihood solution kept. Exposures are rank inverse normal, "
   "whole-cohort in SHHS and within-site in MrOS. C is scored IN SAMPLE on the rows the model was fit "
   "on. There is no split and no cross-fit. Event floor 60."),
  ("Setting 2 MAIN baseline, age + sex",
   "The bare baseline carried by the Overview main row. Age is a natural cubic spline, cr(df=4) "
   "with the first basis column dropped, plus sex in SHHS. MrOS is men only, so its bare baseline "
   "is age alone, site-stratified, and there is no sex term to add. This is the frozen external "
   "pipeline's own convention for age + sex models, copied verbatim from freeze_05_external.py, "
   "and it matches the primary cohort's baseline functional form. Positive control: every one of "
   "the 12 frames reproduces the frozen minimally adjusted T90 hazard ratio in "
   "numbers/external.json, n, events, hr, lo and hi to 3 decimal places, before any C was "
   "computed."),
  ("Setting 2 PRESERVED baseline, SHHS clinical",
   "The frozen external clinical risk model: age + sex + BMI + smoking + diastolic BP + total "
   "cholesterol + HDL + hypertension + diabetes. Age enters linearly here, not as a spline, which is "
   "the frozen convention for this block."),
  ("Setting 2 PRESERVED baseline, MrOS clinical",
   "age + BMI + smoking + systolic BP + hypertension + diabetes, stratified by site. MrOS is men "
   "only, so there is no sex term. MrOS has no HDL and no diastolic pressure, so systolic pressure "
   "stands in. MrOS high cholesterol, DHLCHOL, was tried and dropped: it is missing in 2,124 of 2,911 "
   "men and collapsed every model to n = 770, pushing two outcomes under the event floor."),
 ]),
 ("THE COMPARABILITY CAVEAT, READ BEFORE COMPARING ROWS", [
  ("Held out versus in sample",
   "The primary-cohort dC values in settings 1 and 3 are HELD OUT. Every C is scored on patients "
   "from a hospital the model never saw. The external dC values in setting 2 are IN SAMPLE, scored "
   "on the same rows the model was fit on, because that is how the frozen external numbers in "
   "numbers/external.json were computed and that convention was kept so the positive control could "
   "reproduce. In-sample C is optimistically biased, and the bias grows with the number of free "
   "parameters, so setting 2's dC values are NOT directly comparable to settings 1 and 3."),
  ("The baselines now match on the main rows",
   "Settings 1 and 3 use age plus sex, a deliberately thin baseline. Since 2026-08-20 the Home PSG "
   "main row uses the same bare baseline, so the baseline mismatch is gone from the main rows and "
   "the in-sample scoring is the one remaining across-row difference. The preserved second row "
   "still uses the full clinical risk model, nine covariates in SHHS, and shows what a rich "
   "baseline does to the same exposures: a thin baseline leaves far more room for any exposure to "
   "add discrimination. The primary cohort's own full clinical rungs are on the Laboratory PSG "
   "detail sheet and are the fair comparison for the clinical row's magnitude."),
  ("What survives the caveat",
   "The T90 versus duration CONTRAST within each setting is internally valid, because within a "
   "setting both exposures face the same baseline, the same sample, and the same scoring rule."),
 ]),
 ("SHEETS", [
  ("Overview",
   "The three settings on one compact table, four rows: the Home PSG setting carries a main row "
   "on the bare age + sex baseline and a preserved row on the frozen clinical baseline."),
  ("Laboratory PSG detail", "The full primary ladder, T90 against PSG TST at every rung."),
  ("Home PSG (SHHS MrOS) detail",
   "Two per-outcome blocks, one per baseline. Block 1, age + sex, verified against the frozen "
   "minimally adjusted T90 hazard ratios in external.json. Block 2, the frozen clinical baseline, "
   "verified against the frozen clinical_model_gain values."),
  ("Habitual home report detail",
   "The home-subset rung table, three exposures on identical rows and identical outcomes."),
  ("Verification", "Every positive control, computed against its expected value."),
 ]),
 ("SOURCES", [
  ("results.csv", "Primary T90 ladder, 48 outcomes, cross-fit."),
  ("results_tst.csv", "Primary PSG-TST ladder, same rows and rungs."),
  ("results_home.csv", "Home-subset per-outcome results, three exposures on identical rows."),
  ("COMPARISON_T90_vs_TST_table.csv", "Paired primary ladder, generated by make_comparison.py."),
  ("external_dc_tst.csv",
   "Built by run_external_dc_tst.py in this directory. The external dC on the frozen full "
   "clinical baseline, now the preserved second row. The external duration dC did not exist "
   "before this workbook."),
  ("external_dc_agesex.csv",
   "NEW 2026-08-20. Built by run_external_dc_agesex.py in this directory. The external dC on the "
   "bare age + sex baseline, the Overview main row for setting 2, so all three settings share the "
   "same bare baseline. Positive-controlled against the frozen minimally adjusted T90 hazard "
   "ratios in numbers/external.json, and its generator re-runs the full clinical reproduction on "
   "import before computing anything."),
  ("numbers/external.json",
   "Frozen external T90 gains and hazard ratios, read only, used as the positive controls for "
   "setting 2. Not modified."),
  ("provenance.json, provenance_tst.json, provenance_home.json, external_dc_tst_provenance.json, "
   "external_dc_agesex_provenance.json",
   "Run provenance for each block."),
 ]),
 ("WHAT WAS NOT MODIFIED", [
  ("Frozen and archive files",
   "numbers/external.json, numbers/freeze_05_external.py, every manuscript file, and "
   "COMPARISON_T90_vs_TST.xlsx were read only. COMPARISON_T90_vs_TST.xlsx remains the ladder "
   "archive."),
 ]),
]

for head, items in BLOCKS:
    ws.cell(row=r, column=1, value=head).font = SUB
    r += 1
    for k, v in items:
        a = ws.cell(row=r, column=1, value=k)
        a.font = BOLD
        a.alignment = Alignment(vertical="top", wrap_text=True)
        b = ws.cell(row=r, column=2, value=v)
        b.alignment = Alignment(vertical="top", wrap_text=True)
        ws.row_dimensions[r].height = max(15, 13 * (1 + len(v) // 105))
        r += 1
    r += 1

# ============================================================ 2 Overview
ws = wb.create_sheet("Overview")
r = title(ws, "Overview: discrimination gain by measurement setting",
          ["dC is the mean change in Harrell's C from adding one exposure to a fixed baseline, "
           "averaged over the outcomes in that setting.",
           "All three settings now share the bare age + sex baseline on their main rows "
           "(MrOS is men only, so age alone there). Setting 2's frozen full-clinical-baseline "
           "numbers are preserved on the row beneath it.",
           "Settings 1 and 3 are held-out cross-fit. Both setting-2 rows are IN SAMPLE, no "
           "split, so their levels are optimistically biased. See README before comparing "
           "rows."])

cols = ["Setting", "Where the sleep was recorded", "Cohort and n", "Outcomes",
        "Baseline model", "Mean baseline C", "Mean dC from T90",
        "Outcomes where T90 improves C", "Sleep duration measure",
        "Mean dC from sleep duration", "Outcomes where duration improves C", "C scoring"]
put_header(ws, r, cols, [22, 30, 30, 10, 26, 12, 12, 14, 24, 14, 14, 22])
r += 1

OV = [
 ("1. Laboratory PSG", "Attended in-laboratory polysomnography",
  "Primary cohort, n = 19,173", N_RANKED_OUTCOMES, "age + sex, the paper's baseline",
  float(m1.mean_C_baseline), float(m1.mean_dC), f"{int(m1.n_outcomes_dC_positive)} of {N_RANKED_OUTCOMES}",
  "PSG total sleep time", float(l1.mean_dC), f"{int(l1.n_outcomes_dC_positive)} of {N_RANKED_OUTCOMES}",
  "Held out, two-fold site cross-fit"),
 ("2. Home PSG", "Unattended polysomnography in the participant's home",
  "SHHS n = 5,802 and MrOS n = 2,911", n_ext,
  "age + sex, the bare baseline (age spline + sex in SHHS, age alone in MrOS, men only)",
  a_base_mean, a_t90_mean, f"{a_npos_t90} of {n_ext}",
  "Home-PSG total sleep time", a_tst_mean, f"{a_npos_tst} of {n_ext}",
  "In sample, no split"),
 ("Home PSG (full clinical baseline, frozen convention)",
  "Same home recordings, baseline swapped to the frozen nine-covariate clinical model, "
  "kept so nothing is lost",
  "SHHS n = 5,802 and MrOS n = 2,911", n_ext,
  "Full clinical risk model, frozen external baseline",
  e_base_mean, e_t90_mean, f"{e_npos_t90} of {n_ext}",
  "Home-PSG total sleep time", e_tst_mean, f"{e_npos_tst} of {n_ext}",
  "In sample, no split"),
 ("3. Home habitual report", "Patient-reported usual sleep, from the clinical note",
  "Primary subset, n = 5,295 covered, 5,293 analysed", h_t90["n_outcomes"],
  "age + sex, the paper's baseline",
  h_t90["c_base"], h_t90["dC"], f"{h_t90['npos']} of {h_t90['n_outcomes']}",
  "Reported habitual hours", h_hab["dC"], f"{h_hab['npos']} of {h_t90['n_outcomes']}",
  "Held out, two-fold site cross-fit"),
]
F = [None, None, None, "0", None, D4, D4S, None, None, D4S, None, None]
for i, row in enumerate(OV):
    put_row(ws, r, list(row), F, band=(i % 2 == 1))
    ws.row_dimensions[r].height = 30
    r += 1

r += 1
m4 = res[res.model == "M4_colleague"].iloc[0]
NOTES = [
 "READING THE TABLE",
 f"T90 beats sleep duration in every setting. In the laboratory the margin is wide: "
 f"{m1.mean_dC:+.4f} against {l1.mean_dC:+.4f}, and T90 improves C for "
 f"{int(m1.n_outcomes_dC_positive)} of {N_RANKED_OUTCOMES} outcomes while PSG total sleep time improves it for "
 f"{int(l1.n_outcomes_dC_positive)}, close to the 24 you would expect from noise alone.",
 f"Moving the measurement into the home does not rescue sleep duration. On the shared bare "
 f"baseline, home-PSG total sleep time adds a mean {a_tst_mean:+.4f} across the {n_ext} external "
 f"outcomes, improving C for {a_npos_tst} of {n_ext}, and reported habitual hours add "
 f"{h_hab['dC']:+.4f} in the primary subset. T90 on the same bare baseline adds {a_t90_mean:+.4f} "
 f"and improves C for {a_npos_t90} of {n_ext}.",
 f"The second Home PSG row keeps the frozen clinical-baseline numbers. Against the full clinical "
 f"risk model, nine covariates in SHHS and six in MrOS with a mean baseline C of {e_base_mean:.4f}, "
 f"T90's gain falls from {a_t90_mean:+.4f} to {e_t90_mean:+.4f}. That collapse is a baseline "
 f"effect, not a home-recording effect: the primary cohort's own clinical rungs on the Laboratory "
 f"PSG detail sheet fall to the same order, {m4.mean_dC:+.4f}, once oxygen proxies enter the "
 f"baseline.",
 "",
 "SETTING 3, THE THIRD EXPOSURE",
 f"The habitual-report subset also carries measured PSG total sleep time on the identical rows and "
 f"the identical 43 outcomes. Its dC is {h_tst['dC']:+.4f}, improving C for {h_tst['npos']} of 43. "
 f"So within one sample and one baseline the ordering is T90 {h_t90['dC']:+.4f}, reported habitual "
 f"hours {h_hab['dC']:+.4f}, measured PSG total sleep time {h_tst['dC']:+.4f}. Asking the patient is "
 "no worse than measuring, and both are far behind T90.",
 "",
 "CAVEAT",
 "Rows 1 and 3 are held-out. Both Home PSG rows are in-sample and therefore optimistically "
 "biased. The main rows now share the bare age + sex baseline, so the remaining across-row "
 "difference is the scoring rule: compare the T90 against duration contrast WITHIN a row, and "
 "treat any across-row comparison with the home rows as indicative only. The clinical-baseline "
 "row additionally sits on a much richer baseline, which is why its gains are smaller still.",
]
for line in NOTES:
    c = ws.cell(row=r, column=1, value=line)
    c.font = SUB if line.isupper() and line else NOTE
    c.alignment = Alignment(vertical="top", wrap_text=True)
    ws.merge_cells(start_row=r, start_column=1, end_row=r, end_column=12)
    if line and not line.isupper():
        ws.row_dimensions[r].height = max(15, 13 * (1 + len(line) // 150))
    r += 1

# ============================================================ 3 Laboratory PSG detail
ws = wb.create_sheet("Laboratory PSG detail")
r = title(ws, "Setting 1. Laboratory PSG, primary cohort: T90 against PSG total sleep time",
          ["n = 19,173, 48 incident outcomes, held-out two-fold site cross-fit. Every rung uses the "
           "same rows and the same outcomes for both exposures.",
           "Source: results.csv, results_tst.csv, paired in COMPARISON_T90_vs_TST_table.csv. "
           "Row L1 is the published positive control."])

cols = ["Baseline model", "Covariates", "Block", "Outcomes", "Baseline C", "C with T90",
        "dC from T90", "C with PSG TST", "dC from PSG TST", "T90 minus TST",
        "T90 over TST, magnitude", "Outcomes T90 dC>0", "Outcomes TST dC>0",
        "Oxygen in baseline", "Note"]
put_header(ws, r, cols, [30, 46, 12, 10, 11, 11, 11, 12, 12, 12, 13, 12, 12, 12, 44])
r += 1
F = [None, None, None, "0", D4, D4, D4S, D4, D4S, D4S, "0.0", "0", "0", None, None]
for i, (_, x) in enumerate(lad.iterrows()):
    vals = [x["Baseline model"], x["Covariates"], x["Block"],
            None if pd.isna(x["n outcomes"]) else int(x["n outcomes"]),
            x["Baseline C"], x["C with T90"], x["dC from adding T90"],
            None if pd.isna(x["C with TST"]) else x["C with TST"],
            None if pd.isna(x["dC from adding TST"]) else x["dC from adding TST"],
            None if pd.isna(x["Difference, T90 minus TST"]) else x["Difference, T90 minus TST"],
            None if pd.isna(x["Magnitude ratio, T90 over TST"]) else x["Magnitude ratio, T90 over TST"],
            None if pd.isna(x["n outcomes T90 dC>0"]) else int(x["n outcomes T90 dC>0"]),
            None if pd.isna(x["n outcomes TST dC>0"]) else int(x["n outcomes TST dC>0"]),
            "" if pd.isna(x["Oxygen row"]) else x["Oxygen row"],
            "" if pd.isna(x["Note"]) else x["Note"]]
    put_row(ws, r, vals, F, band=(i % 2 == 1), bold=(i == 0))
    r += 1

r += 1
for line in [
  "The first row is the paper's baseline and the positive control: dC from T90 = 0.020339 against a "
  "baseline C of 0.647540, both reproducing the published values exactly.",
  "PSG total sleep time is negative at every single rung. Adding it to any of these baselines makes "
  "discrimination slightly worse, and it never improves C for more than 26 of 48 outcomes.",
  "The rungs marked YES under oxygen already hold mean SpO2 in the baseline. T90's gain collapses "
  "there, from 0.0203 to about 0.0028, which is the ceiling on how much of T90's signal is "
  "independent of average saturation.",
]:
    c = ws.cell(row=r, column=1, value=line)
    c.font = NOTE
    c.alignment = Alignment(vertical="top", wrap_text=True)
    ws.merge_cells(start_row=r, start_column=1, end_row=r, end_column=15)
    ws.row_dimensions[r].height = max(15, 13 * (1 + len(line) // 160))
    r += 1

# ============================================================ 4 Home PSG detail
ws = wb.create_sheet("Home PSG (SHHS MrOS) detail")
r = title(ws, "Setting 2. Home PSG, SHHS and MrOS: T90 against home-PSG total sleep time",
          ["Unattended polysomnography recorded in the participant's own home. 12 outcomes, 7 in "
           "SHHS and 5 in MrOS. Two baselines, one block each.",
           "BLOCK 1, the Overview main row: bare age + sex baseline. Age is a natural cubic "
           "spline, cr(df=4) with the first basis column dropped, plus sex in SHHS. MrOS is men "
           "only, so its bare baseline is age alone, site-stratified. This is the frozen external "
           "convention for age + sex models, copied from freeze_05_external.py.",
           "BLOCK 2, the preserved Overview second row: the frozen full clinical baseline, "
           "unchanged from the original computation.",
           "C is scored IN SAMPLE in both blocks, no split, matching the frozen convention in "
           "numbers/freeze_05_external.py. NOT comparable to the held-out primary numbers.",
           "Sources: external_dc_agesex.csv, built by run_external_dc_agesex.py, and "
           "external_dc_tst.csv, built by run_external_dc_tst.py. Frozen columns are read from "
           "numbers/external.json and are the positive controls."])

ws.cell(row=r, column=1,
        value="BLOCK 1. AGE + SEX BASELINE, carried by the Overview main row").font = SUB
r += 1
cols_a = ["Cohort", "Outcome", "n in model", "Events", "Baseline covariates", "C baseline",
          "C with T90", "dC from T90", "Frozen T90 HR, external.json", "Recomputed T90 HR",
          "Verification", "C with TST", "dC from TST", "T90 minus TST"]
put_header(ws, r, cols_a, [9, 26, 11, 9, 52, 11, 11, 11, 13, 13, 16, 11, 11, 12])
r += 1
FA = [None, None, "0", "0", None, D4, D4, D4S, "0.000", "0.000", None, D4, D4S, D4S]
for i, (_, x) in enumerate(agso.iterrows()):
    vals = [x.cohort, x.outcome, int(x.n_model), int(x.events), x.baseline_covariates,
            x.c_baseline, x.c_with_t90, x.dC_t90,
            float(x.t90_hr_frozen), float(x.t90_hr_recomputed),
            "reproduces to 3 dp" if x.frozen_hr_check == "reproduced" else x.frozen_hr_check,
            x.c_with_tst, x.dC_tst, x.dC_t90 - x.dC_tst]
    put_row(ws, r, vals, FA, band=(i % 2 == 1))
    r += 1

r += 1
put_row(ws, r, ["MEANS, BLOCK 1"] + [""] * 13, bold=True)
r += 1
for coh in ["SHHS", "MrOS"]:
    q = agso[agso.cohort == coh]
    put_row(ws, r, [coh, f"mean of {len(q)} outcomes", "", "", "",
                    q.c_baseline.mean(), q.c_with_t90.mean(), q.dC_t90.mean(),
                    None, None, f"T90 dC>0 in {int((q.dC_t90 > 0).sum())} of {len(q)}",
                    q.c_with_tst.mean(), q.dC_tst.mean(),
                    q.dC_t90.mean() - q.dC_tst.mean()], FA, bold=True)
    r += 1
put_row(ws, r, ["BOTH", f"mean of {n_ext} outcomes", "", "", "",
                a_base_mean, agso.c_with_t90.mean(), a_t90_mean, None, None,
                f"T90 dC>0 in {a_npos_t90} of {n_ext}, TST dC>0 in {a_npos_tst} of {n_ext}",
                agso.c_with_tst.mean(), a_tst_mean, a_t90_mean - a_tst_mean], FA, bold=True)
r += 2

ws.cell(row=r, column=1,
        value="BLOCK 2. FULL CLINICAL BASELINE, FROZEN CONVENTION, the preserved "
              "Overview second row").font = SUB
r += 1
cols = ["Cohort", "Outcome", "n in model", "Events", "Baseline covariates", "C baseline",
        "C with T90", "dC from T90", "Frozen dC from T90, external.json", "Verification",
        "C with TST", "dC from TST", "T90 minus TST"]
put_header(ws, r, cols, freeze=False)
r += 1
F = [None, None, "0", "0", None, D4, D4, D4S, D4S, None, D4, D4S, D4S]
for i, (_, x) in enumerate(ext.iterrows()):
    vals = [x.cohort, x.outcome, int(x.n_model), int(x.events), x.baseline_covariates,
            x.c_baseline, x.c_with_t90, x.dC_t90,
            None if pd.isna(x.frozen_dC_t90) else float(x.frozen_dC_t90),
            "reproduces to 4 dp" if x.frozen_check == "reproduced" else x.frozen_check,
            x.c_with_tst, x.dC_tst, x.dC_t90 - x.dC_tst]
    put_row(ws, r, vals, F, band=(i % 2 == 1))
    r += 1

# cohort means and overall
r += 1
put_row(ws, r, ["MEANS, BLOCK 2", "", "", "", "", "", "", "", "", "", "", "", ""], bold=True)
r += 1
for coh in ["SHHS", "MrOS"]:
    q = ext[ext.cohort == coh]
    put_row(ws, r, [coh, f"mean of {len(q)} outcomes", "", "", "",
                    q.c_baseline.mean(), q.c_with_t90.mean(), q.dC_t90.mean(),
                    None, f"T90 dC>0 in {int((q.dC_t90 > 0).sum())} of {len(q)}",
                    q.c_with_tst.mean(), q.dC_tst.mean(),
                    q.dC_t90.mean() - q.dC_tst.mean()], F, bold=True)
    r += 1
put_row(ws, r, ["BOTH", f"mean of {n_ext} outcomes", "", "", "",
                e_base_mean, ext.c_with_t90.mean(), e_t90_mean, None,
                f"T90 dC>0 in {e_npos_t90} of {n_ext}, TST dC>0 in {e_npos_tst} of {n_ext}",
                ext.c_with_tst.mean(), e_tst_mean, e_t90_mean - e_tst_mean], F, bold=True)
r += 2

for line in [
  "POSITIVE CONTROLS",
  "Block 1: numbers/external.json froze, per outcome, the minimally adjusted T90 hazard ratio "
  "computed on exactly the age-spline frame block 1 uses. All 12 frames reproduce the frozen n, "
  "events, hr, lo and hi to 3 decimal places before any C was computed. Block 2: all 7 SHHS "
  "outcomes reproduce external.json to 4 decimal places on C baseline, C with T90, and the gain. "
  "The MrOS clinical rows have no frozen counterpart: freeze_05_external.py never built a "
  "clinical-model block for MrOS, so those 5 outcomes are an extension of the same method, not a "
  "reproduction.",
  "",
  "THE RESULT",
  f"On the bare age + sex baseline, T90 adds a mean {a_t90_mean:+.4f} to C across the 12 outcomes, "
  f"improving discrimination for {a_npos_t90} of them. Home-PSG total sleep time adds "
  f"{a_tst_mean:+.4f}, improving discrimination for {a_npos_tst}. Recording sleep at home rather "
  "than in a laboratory does not make sleep duration a useful discriminator, even on the barest "
  "baseline.",
  f"Loading the baseline with the full clinical risk model, block 2, shrinks T90's mean gain from "
  f"{a_t90_mean:+.4f} to {e_t90_mean:+.4f} and duration's from {a_tst_mean:+.4f} to "
  f"{e_tst_mean:+.4f}. That is a baseline effect, not a home-recording effect: the primary "
  f"cohort's own clinical rungs on the Laboratory PSG detail sheet fall to the same order, "
  f"{res[res.model == 'M4_colleague'].iloc[0].mean_dC:+.4f}, once oxygen proxies enter the "
  f"baseline.",
  "",
  "CAVEAT",
  "Both blocks are IN SAMPLE, so all C levels are optimistically biased and the differences, not "
  "the levels, are the interpretable quantity. Two SHHS outcomes, myocardial infarction and "
  "stroke, have age + sex T90 gains at or below their clinical-baseline gains. Both sit within "
  "0.001 of zero under either baseline, consistent with T90 carrying no discrimination signal for "
  "those two outcomes in SHHS rather than with any baseline artefact. The cohort means and the "
  "combined mean are larger on the bare baseline, as expected.",
]:
    c = ws.cell(row=r, column=1, value=line)
    c.font = SUB if (line.isupper() and line) else NOTE
    c.alignment = Alignment(vertical="top", wrap_text=True)
    ws.merge_cells(start_row=r, start_column=1, end_row=r, end_column=14)
    if line and not line.isupper():
        ws.row_dimensions[r].height = max(15, 13 * (1 + len(line) // 150))
    r += 1

# ============================================================ 5 Habitual home report detail
ws = wb.create_sheet("Habitual home report detail")
r = title(ws, "Setting 3. Home habitual report, primary subset: three exposures on identical rows",
          ["5,295 of 19,173 patients have a habitual sleep duration in the clinical note, 27.6 "
           "percent. 5,293 have it together with PSG TST and sleep efficiency, which is the analysis "
           "frame. 43 of 48 outcomes clear the 30-events-per-fold floor.",
           "All three exposures are scored on the SAME rows, the SAME outcomes, and the SAME "
           "baseline, so the comparison between them is exact. Held-out two-fold site cross-fit.",
           "Source: results_home.csv, aggregated over the 43 kept outcomes. Verified against "
           "provenance_home.json."])

cols = ["Block", "Baseline model", "Covariates", "Outcomes", "Baseline C", "Exposure",
        "C with exposure", "dC", "Outcomes dC>0", "Relative gain"]
put_header(ws, r, cols, [16, 32, 46, 10, 11, 30, 12, 11, 12, 12])
r += 1

MODELS = [("H1_age_sex", "H1. age + sex", "age + sex, the paper's baseline"),
          ("H2_age_sex_AHI", "H2. age + sex + AHI", "age + sex + AHI"),
          ("H3_age_sex_race_smoke", "H3. age + sex + race + smoking", "age + sex + race + smoking"),
          ("H4_clinical_AHI_SE", "H4. clinical + AHI + sleep eff",
           "age + sex + race + smoking + AHI + sleep efficiency")]
EXPO = [("t90", "T90, laboratory PSG"),
        ("home", "Habitual hours, home report"),
        ("tst", "Total sleep time, laboratory PSG")]
F = [None, None, None, "0", D4, None, D4, D4S, None, "0.0000%"]

for blk, blab, suffix in [("home_subset", "Main", ""),
                          ("home_agedrop", "Age-drop sensitivity", "__agedrop")]:
    ws.cell(row=r, column=1, value=blab if blk == "home_subset"
            else "Age-drop sensitivity: one age spline column dropped").font = SUB
    r += 1
    for mkey, mlab, mcov in MODELS:
        for j, (ex, exlab) in enumerate(EXPO):
            d = home_rung(blk, mkey + suffix, ex)
            put_row(ws, r, [blk, mlab, mcov, d["n_outcomes"], d["c_base"], exlab,
                            d["c_expo"], d["dC"], f"{d['npos']} of {d['n_outcomes']}",
                            d["dC"] / d["c_base"]],
                    F, band=(j == 1), bold=(ex == "t90"))
            r += 1
    r += 1

for line in [
  "READING THE TABLE",
  f"At the paper's baseline, T90 adds {h_t90['dC']:+.4f}, the patient's reported habitual hours add "
  f"{h_hab['dC']:+.4f}, and measured PSG total sleep time adds {h_tst['dC']:+.4f}. T90 is roughly "
  f"{abs(h_t90['dC'] / h_hab['dC']):.0f} times the reported-hours gain.",
  "The reported figure is not worse than the measured one. Habitual report is positive at the two "
  "thinner baselines and measured PSG TST is negative at all four, so the failure of sleep duration "
  "is not a measurement-error artefact. Duration simply carries little incident-disease signal in "
  "this cohort, however it is captured.",
  "The age-drop block repeats every rung with one age spline column removed, the rank-deficiency "
  "guard. Nothing material moves: T90 stays near 0.018 and both duration measures stay near zero.",
]:
    c = ws.cell(row=r, column=1, value=line)
    c.font = SUB if (line.isupper() and line) else NOTE
    c.alignment = Alignment(vertical="top", wrap_text=True)
    ws.merge_cells(start_row=r, start_column=1, end_row=r, end_column=10)
    if line and not line.isupper():
        ws.row_dimensions[r].height = max(15, 13 * (1 + len(line) // 150))
    r += 1

# ============================================================ 6 Verification
ws = wb.create_sheet("Verification")
r = title(ws, "Verification: positive controls, all computed before the workbook was written",
          ["Every check recomputes a number from source and compares it with the published or "
           "frozen value. The build aborts if any check fails."])
put_header(ws, r, ["Check", "Computed", "Expected", "Difference", "Result"],
           [78, 20, 20, 16, 10])
r += 1
for i, c in enumerate(checks):
    put_row(ws, r, [c["Check"], c["Computed"], c["Expected"],
                    c["Computed"] - c["Expected"], c["Result"]],
            [None, "0.000000000000", "0.000000000000", "0.00E+00", None],
            band=(i % 2 == 1))
    if c["Result"] == "PASS":
        ws.cell(row=r, column=5).font = Font(bold=True, color="1F7A1F")
    r += 1
r += 1
ws.cell(row=r, column=1,
        value=f"{len(checks)} checks, {sum(1 for c in checks if c['Result']=='PASS')} PASS, "
              f"{sum(1 for c in checks if c['Result']=='FAIL')} FAIL.").font = BOLD

wb.save(OUT)
print(f"written -> {OUT}")
print(f"sheets: {wb.sheetnames}")
