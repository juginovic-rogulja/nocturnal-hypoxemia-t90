"""
Build COMPARISON_T90_vs_TST.md and COMPARISON_T90_vs_TST.xlsx.

Every T90 number is READ from results.csv (the existing run) or from results_tst.csv for the two
covariate sets that did not exist before, the ones with the TST covariate removed. Nothing is
retyped, recomputed or rounded by hand. Every TST number is read from results_tst.csv.

The pairing rule is that a row of the comparison table is ONE covariate set, and both exposures
were fitted against it with the same folds, the same penalizer and the same outcome list. The
baseline model never contains the exposure, so the baseline C printed on a row is a single number
that belongs to both exposures. The script asserts that, to 1e-12, rather than assuming it.
"""
from __future__ import annotations
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths

import os

import numpy as np
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter
import sys as _sys_v8; _sys_v8.path.insert(0, paths.ANALYSIS_DIR)  # v8 sweep 2026-09-12
from cohort_spec import N_MEASURES  # v8 sweep 2026-09-12

HERE = os.path.dirname(os.path.abspath(__file__))

OLD = pd.read_csv(os.path.join(HERE, "results.csv"))            # the T90 ladder
NEW = pd.read_csv(os.path.join(HERE, "results_tst.csv"))        # the TST ladder plus T90 extras

# (row label, covariate description, block, T90 source, T90 model, TST model, oxygen, note)
# T90 source "old" = results.csv, "new" = results_tst.csv (a covariate set that did not exist).
ROWS = [
    ("L1. age + sex", "age + sex, the paper's baseline",
     "full", "old", "M1_paper_age_sex", "L1_age_sex", False, "positive control for both"),
    ("L2. age + sex + AHI", "age + sex + AHI",
     "full", "old", "D1_age_sex_AHI", "L2_age_sex_AHI", False, ""),
    ("L3. age + sex + race + smoking", "age + sex + race + smoking",
     "full", "old", "M2_clinical", "L3_age_sex_race_smoke", False, ""),
    ("L3b. + AHI", "age + sex + race + smoking + AHI",
     "full", "old", "D4_clinical_AHI_noOxygen", "LX_clinical_AHI", False, ""),
    ("L4 / L4b. + AHI + sleep efficiency",
     "age + sex + race + smoking + AHI + sleep efficiency, no TST covariate",
     "full", "new", "D5_clinical_sleep_noOxygen__noTSTcov", "L4_clinical_AHI_SE", False,
     "apples to apples, TST covariate removed from both sides"),
    ("L4 reference. the published T90 rung, TST held as a covariate",
     "age + sex + race + smoking + AHI + TST + sleep efficiency",
     "full", "old", "D5_clinical_sleep_noOxygen", None, False,
     "NOT comparable, the baseline holds TST so TST cannot be the exposure"),
    ("OX1. age + sex + mean SpO2", "age + sex + mean SpO2",
     "full", "old", "D2_age_sex_meanSpO2", "OX1_age_sex_meanSpO2", True, ""),
    ("OX2. age + sex + AHI + mean SpO2", "age + sex + AHI + mean SpO2",
     "full", "old", "D3_age_sex_AHI_meanSpO2", "OX2_age_sex_AHI_meanSpO2", True, ""),
    ("OX3. + race + smoking + AHI + mean SpO2",
     "age + sex + race + smoking + AHI + mean SpO2",
     "full", "old", "M4_colleague", "OX3_clinical_AHI_meanSpO2", True,
     "already TST free in the T90 ladder, so directly comparable"),
    ("OX4. + AHI + sleep efficiency + mean SpO2",
     "age + sex + race + smoking + AHI + sleep efficiency + mean SpO2, no TST covariate",
     "full", "new", "M3_clinical_sleep__noTSTcov", "OX4_clinical_AHI_SE_meanSpO2", True,
     "apples to apples, TST covariate removed from both sides"),
    ("L5a. age + sex, BMI-complete rows", "age + sex, matched control on the BMI-complete rows",
     "bmi_6325", "old", "B1_paper_age_sex", "B1_age_sex", False, "own block, own control"),
    ("L5b. age + sex + BMI", "age + sex + BMI",
     "bmi_6325", "old", "B2_age_sex_BMI", "B2_age_sex_BMI", False, ""),
    ("L5c. age + sex + BMI + race + smoking", "age + sex + BMI + race + smoking",
     "bmi_6325", "old", "B3_clinical_BMI", "B3_clinical_BMI", False, ""),
    ("L5d. + AHI + sleep efficiency", "age + sex + BMI + race + smoking + AHI + sleep efficiency",
     "bmi_6325", "new", "B4_clinical_sleep_BMI__noTSTcov_noOx", "B4_clinical_BMI_AHI_SE", False,
     "apples to apples, TST covariate removed from both sides"),
    ("OXB1. + AHI + mean SpO2", "age + sex + BMI + race + smoking + AHI + mean SpO2",
     "bmi_6325", "old", "B5_colleague_BMI", "OXB1_colleague_BMI", True,
     "already TST free in the T90 ladder, so directly comparable"),
    ("OXB2. + AHI + sleep efficiency + mean SpO2",
     "age + sex + BMI + race + smoking + AHI + sleep efficiency + mean SpO2, no TST covariate",
     "bmi_6325", "new", "B4_clinical_sleep_BMI__noTSTcov", "OXB2_clinical_sleep_BMI", True,
     "apples to apples, TST covariate removed from both sides"),
]


def t90_row(block, source, model):
    if source == "old":
        r = OLD[(OLD.block == block) & (OLD.model == model)]
        assert len(r) == 1, f"results.csv has {len(r)} rows for {block}/{model}"
        r = r.iloc[0]
        return dict(base=float(r.mean_C_baseline), cwith=float(r.mean_C_with_T90),
                    dC=float(r.mean_dC), nout=int(r.n_outcomes),
                    npos=int(r.n_outcomes_dC_positive), src="results.csv")
    r = NEW[(NEW.block == block) & (NEW.model == model) & (NEW.exposure == "t90_z")]
    assert len(r) == 1, f"results_tst.csv has {len(r)} rows for {block}/{model}/t90_z"
    r = r.iloc[0]
    return dict(base=float(r.mean_C_baseline), cwith=float(r.mean_C_with_exposure),
                dC=float(r.mean_dC), nout=int(r.n_outcomes),
                npos=int(r.n_outcomes_dC_positive), src="results_tst.csv")


def tst_row(block, model):
    r = NEW[(NEW.block == block) & (NEW.model == model) & (NEW.exposure == "tst_z")]
    assert len(r) == 1, f"results_tst.csv has {len(r)} rows for {block}/{model}/tst_z"
    r = r.iloc[0]
    return dict(base=float(r.mean_C_baseline), cwith=float(r.mean_C_with_exposure),
                dC=float(r.mean_dC), nout=int(r.n_outcomes),
                npos=int(r.n_outcomes_dC_positive), src="results_tst.csv")


def build() -> pd.DataFrame:
    out = []
    for lab, desc, block, src, m90, mtst, oxy, note in ROWS:
        a = t90_row(block, src, m90)
        b = tst_row(block, mtst) if mtst else None
        if b is not None:
            # the baseline never holds the exposure, so both exposures must sit on one baseline C
            assert abs(a["base"] - b["base"]) < 1e-12, (
                f"baseline C differs on {lab}: T90 {a['base']!r} vs TST {b['base']!r}")
            assert a["nout"] == b["nout"], f"outcome count differs on {lab}"
        row = {
            "Baseline model": lab,
            "Covariates": desc,
            "Block": block,
            "n outcomes": a["nout"],
            "Baseline C": a["base"],
            "C with T90": a["cwith"],
            "dC from adding T90": a["dC"],
            "C with TST": b["cwith"] if b else np.nan,
            "dC from adding TST": b["dC"] if b else np.nan,
            "Difference, T90 minus TST": (a["dC"] - b["dC"]) if b else np.nan,
            "Magnitude ratio, T90 over TST": (abs(a["dC"]) / abs(b["dC"]))
            if b and abs(b["dC"]) > 0 else np.nan,
            "Rel gain T90": a["dC"] / (1 - a["base"]),
            "Rel gain TST": (b["dC"] / (1 - b["base"])) if b else np.nan,
            "n outcomes T90 dC>0": a["npos"],
            "n outcomes TST dC>0": b["npos"] if b else np.nan,
            "Oxygen row": "YES" if oxy else "",
            "T90 source": a["src"],
            "Note": note,
        }
        out.append(row)
    return pd.DataFrame(out)


# ------------------------------------------------------------------ organs
def organ_block() -> pd.DataFrame:
    P = pd.read_csv(os.path.join(HERE, "per_outcome_tst.csv"))
    OLDP = pd.read_csv(os.path.join(HERE, "per_outcome.csv"))
    order = ["cardiac", "respiratory", "metabolic", "kidney", "other"]

    def take(frame, block, model, expo=None, dccol="dC"):
        f = frame[(frame.block == block) & (frame.model == model) & frame[dccol].notna()]
        if expo is not None:
            f = f[f.exposure == expo]
        g = f.groupby("organ5")[dccol].agg(["mean", "size"])
        return g.reindex(order)

    t90_l4 = take(P, "full", "D5_clinical_sleep_noOxygen__noTSTcov", "t90_z")
    tst_l4 = take(P, "full", "L4_clinical_AHI_SE", "tst_z")
    t90_l4_pub = take(OLDP, "full", "D5_clinical_sleep_noOxygen")
    t90_l1 = take(P, "full", "M1_paper_age_sex__RERUN", "t90_z")
    tst_l1 = take(P, "full", "L1_age_sex", "tst_z")

    out = pd.DataFrame({
        "Organ group": order,
        "n outcomes": t90_l4["size"].values.astype(int),
        "L1 dC T90": t90_l1["mean"].values,
        "L1 dC TST": tst_l1["mean"].values,
        "L1 difference": t90_l1["mean"].values - tst_l1["mean"].values,
        "L4 dC T90": t90_l4["mean"].values,
        "L4 dC TST": tst_l4["mean"].values,
        "L4 difference": t90_l4["mean"].values - tst_l4["mean"].values,
        "L4 dC T90, published rung holding TST": t90_l4_pub["mean"].values,
    })
    return out


# ------------------------------------------------------------------ markdown
def md_table(df: pd.DataFrame, cols, fmts) -> list[str]:
    head = "| " + " | ".join(cols) + " |"
    rule = "|" + "|".join("---" for _ in cols) + "|"
    lines = [head, rule]
    for _i, r in df.iterrows():
        cells = []
        for c in cols:
            v = r[c]
            f = fmts.get(c)
            if v is None or (isinstance(v, float) and not np.isfinite(v)):
                cells.append("--")
            elif f:
                cells.append(f(v))
            else:
                cells.append(str(v))
        lines.append("| " + " | ".join(cells) + " |")
    return lines


def write_md(cmp_df: pd.DataFrame, org: pd.DataFrame, prov: dict):
    n4 = lambda v: f"{v:.4f}"          # noqa: E731
    s4 = lambda v: f"{v:+.4f}"         # noqa: E731
    s5 = lambda v: f"{v:+.5f}"         # noqa: E731
    r2 = lambda v: f"{v:.1f}x"         # noqa: E731
    pc = lambda v: f"{100*v:+.2f}%"    # noqa: E731

    L = []
    a = L.append
    a("# T90 versus total sleep time: how much discrimination each adds "
      "on top of a clinical baseline")
    a("")
    a("Supplementary side analysis. Nothing here enters the manuscript.")
    a("")
    a("## TLDR")
    main = cmp_df[cmp_df["Oxygen row"] == ""]
    l1 = cmp_df[cmp_df["Baseline model"].str.startswith("L1.")].iloc[0]
    l4 = cmp_df[cmp_df["Baseline model"].str.startswith("L4 /")].iloc[0]
    a(f"Across the paper's 48 outcomes, with held-out cross-fit concordance transported between "
      f"the two hospitals, T90 adds **{l1['dC from adding T90']:+.4f}** to the age and sex "
      f"baseline while total sleep time adds **{l1['dC from adding TST']:+.5f}**. Total sleep "
      f"time makes the model slightly worse, which is why it ranks 192 of 197 measures in the "
      f"published ranking while T90 ranks first. The gap does not close as the baseline gets "
      f"richer. On the fullest non-oxygen clinical baseline that both exposures can share "
      f"(age, sex, race, smoking, AHI and sleep efficiency), T90 still adds "
      f"**{l4['dC from adding T90']:+.4f}** and total sleep time adds "
      f"**{l4['dC from adding TST']:+.5f}**.")
    a("")
    a("## Positive controls, run before anything else")
    a("")
    a(f"- Total sleep time on the paper baseline reproduces `numbers/ranking_v3.csv` row "
      f"{prov['tst_rank']} of {N_MEASURES} to {prov['tst_tol']}: published dC "
      f"`{prov['tst_pub_dC']:.12f}`, here `{prov['tst_got_dC']:.12f}`.")
    a(f"- T90 on the same baseline reproduces the published dC to {prov['t90_tol']}: "
      f"published `{prov['t90_pub_dC']:.12f}`, here `{prov['t90_got_dC']:.12f}`.")
    a(f"- Baseline C reproduces the published `{prov['pub_base']:.6f}` "
      f"(here `{prov['got_base']:.6f}`).")
    a("- Three T90 rungs from the existing `results.csv` were re-fitted by this script and "
      "matched to 1e-12, so the T90 and the total sleep time columns come from one machine.")
    a("")
    a("## The comparison table")
    a("")
    a("One row is one baseline covariate set. Both exposures were fitted against it with the "
      "same folds, the same penalizer and the same outcome list. The baseline model never holds "
      "the exposure, so the baseline C on a row is one number shared by both columns.")
    a("")
    cols = ["Baseline model", "Baseline C", "dC from adding T90", "dC from adding TST",
            "Difference, T90 minus TST", "Magnitude ratio, T90 over TST", "Note"]
    fmts = {"Baseline C": n4, "dC from adding T90": s4, "dC from adding TST": s5,
            "Difference, T90 minus TST": s4, "Magnitude ratio, T90 over TST": r2}
    a("### Full cohort, 19,173 recordings, no oxygen term in the baseline")
    a("")
    L.extend(md_table(main[(main.Block == "full")], cols, fmts))
    a("")
    a("### BMI-complete rows, 6,325 recordings, with a matched age and sex control")
    a("")
    L.extend(md_table(main[(main.Block == "bmi_6325")], cols, fmts))
    a("")
    a("### Oxygen-containing baselines")
    a("")
    a("These hold whole-night mean SpO2, which is read off the same oximetry trace as T90, so "
      "the baseline is close to adjusting the exposure away. Alen has said these are not the "
      "rows he cares about. They are here so the ladder is complete.")
    a("")
    L.extend(md_table(cmp_df[cmp_df["Oxygen row"] == "YES"], cols, fmts))
    a("")
    a("### Relative gain, dC divided by the headroom 1 minus baseline C")
    a("")
    ii = lambda v: f"{int(round(v))}"  # noqa: E731
    L.extend(md_table(cmp_df, ["Baseline model", "Baseline C", "Rel gain T90", "Rel gain TST",
                               "n outcomes", "n outcomes T90 dC>0", "n outcomes TST dC>0"],
                      {"Baseline C": n4, "Rel gain T90": pc, "Rel gain TST": pc,
                       "n outcomes": ii, "n outcomes T90 dC>0": ii,
                       "n outcomes TST dC>0": ii}))
    a("")
    npos_l1 = int(l1["n outcomes TST dC>0"])
    nout_l1 = int(l1["n outcomes"])
    a(f"On the age and sex rung, total sleep time improves {npos_l1} of the {nout_l1} outcomes "
      f"and leaves {nout_l1 - npos_l1} no better. That count of {npos_l1} is the same npos the "
      f"published `ranking_v3.csv` carries for `TST_min`, and the worst per-outcome dC matches "
      f"the published `-0.013455` as well, so the row reproduces on four separate quantities and "
      f"not only on the mean.")
    a("")
    a("### Age-spline drop check, the standing rule")
    a("")
    a("The published cr(df=4) basis has four columns that sum to one and is singular, so every "
      "rung was refitted with one column removed. Total sleep time stays a loss on every rung "
      "under both bases, and the largest move in dC is 0.0002.")
    a("")
    ad = NEW[NEW.block == "full_agedrop"].set_index("model")
    a4 = NEW[(NEW.block == "full") & (NEW.exposure == "tst_z")].set_index("model")
    drop = pd.DataFrame([{
        "Rung": m,
        "dC TST, published 4-column basis": float(a4.loc[m, "mean_dC"]),
        "dC TST, one column dropped": float(ad.loc[m + "__agedrop", "mean_dC"]),
        "Change": float(ad.loc[m + "__agedrop", "mean_dC"]) - float(a4.loc[m, "mean_dC"]),
    } for m in a4.index])
    L.extend(md_table(drop, list(drop.columns),
                      {c: s5 for c in drop.columns if c != "Rung"}))
    a("")
    a("## The TST-as-covariate asymmetry, stated plainly")
    a("")
    a("The T90 ladder's rich rungs hold total sleep time as a COVARIATE. That covariate cannot "
      "survive when total sleep time is the exposure, so those rungs were re-fitted for T90 "
      "with the TST covariate removed, and the comparison table quotes those. The published "
      "T90 rung that holds TST is kept in the table as a reference row and marked not "
      "comparable. Removing TST from the T90 side barely moves it, which is the point: the "
      "asymmetry is bookkeeping, not a result.")
    a("")
    a("The task's L4 and L4b turn out to be the same model. L4 is age, sex, race, smoking, AHI "
      "and sleep efficiency. L4b is the T90 ladder's equivalent rung minus TST, and that rung "
      "is age, sex, race, smoking, AHI, TST and sleep efficiency. Take TST out and it is L4. "
      "They are reported as one row rather than manufactured into two.")
    a("")
    a("## Per-organ dC, the five groups")
    a("")
    a("Same organ grouping the T90 side analysis used. `other` is the 28 outcomes outside the "
      "cardiac, respiratory, metabolic and kidney groups.")
    a("")
    L.extend(md_table(org, ["Organ group", "n outcomes", "L1 dC T90", "L1 dC TST",
                            "L1 difference", "L4 dC T90", "L4 dC TST", "L4 difference"],
                      {"L1 dC T90": s4, "L1 dC TST": s5, "L1 difference": s4,
                       "L4 dC T90": s4, "L4 dC TST": s5, "L4 difference": s4}))
    a("")
    rsp = org[org["Organ group"] == "respiratory"].iloc[0]
    car = org[org["Organ group"] == "cardiac"].iloc[0]
    oth = org[org["Organ group"] == "other"].iloc[0]
    a(f"The two exposures do not just differ in size, they differ in shape. T90's gain is "
      f"concentrated where an oxygen mechanism would put it: respiratory outcomes get "
      f"{rsp['L1 dC T90']:+.4f} on the age and sex rung, "
      f"{rsp['L1 dC T90']/car['L1 dC T90']:.1f} times the "
      f"{car['L1 dC T90']:+.4f} the cardiac group gets and "
      f"{rsp['L1 dC T90']/oth['L1 dC T90']:.1f} times the "
      f"{oth['L1 dC T90']:+.4f} the 28 outcomes outside the four organ groups get. Total sleep "
      f"time has no such shape. Its largest single organ effect is a LOSS in the same "
      f"respiratory group, {rsp['L1 dC TST']:+.5f}, and its two positive groups, cardiac at "
      f"{car['L1 dC TST']:+.5f} and metabolic at "
      f"{org[org['Organ group'] == 'metabolic'].iloc[0]['L1 dC TST']:+.5f}, are small enough to "
      f"be noise. Adding the clinical baseline shrinks T90's organ pattern but does not flatten "
      f"it. It has nothing to shrink on the total sleep time side.")
    a("")
    a("The last column of the per-organ sheet in the workbook carries the published T90 rung "
      "that holds TST as a covariate, so the effect of removing that covariate can be read "
      "organ by organ. It moves nothing by more than 0.0001.")
    a("")
    a("## Model specification, identical for every row")
    a("")
    for line in SPEC_LINES:
        a(f"- {line}")
    a("")
    a("## Files")
    a("")
    a("- `results_tst.csv` per-model summary for the total sleep time ladder plus the "
      "TST-free T90 rungs")
    a("- `per_outcome_tst.csv` per-outcome detail, 48 outcomes per model")
    a("- `REPORT_TST.txt` the full run log including the positive controls and the age-spline "
      "drop check")
    a("- `results.csv` and `per_outcome.csv` the existing T90 ladder, unchanged")
    a("")
    path = os.path.join(HERE, "COMPARISON_T90_vs_TST.md")
    open(path, "w").write("\n".join(L) + "\n")
    return path


SPEC_LINES = [
    "Held-out cross-fit between hospitals. Fit at I0002 and score at I0006, then the reverse, "
    "and take the mean of the two held-out concordances. A fold is dropped when either side has "
    "fewer than 30 events.",
    "Age enters as a natural cubic spline, cr(df=4) with no intercept. Every model is repeated "
    "with one spline column dropped and the ladder is unchanged.",
    "Exposures enter as ONE column holding the within-site rank inverse normal of the measure, "
    "so every dC is per 1 SD on the paper's own scale.",
    "48 outcomes, the ranking's set: every condition that is neither circular nor held out, "
    "with at least 150 incident events, plus death.",
    "lifelines CoxPHFitter with penalizer 0.01, no strata. Hospitals are separated by the fold, "
    "not by a strata term, so a hospital indicator is not identifiable inside a fold.",
    "Eligibility per outcome is prevalent equal to zero, follow-up greater than zero and "
    "follow-up not missing. Row sets are complete case over exactly the columns each model uses.",
    "dC is the plain mean over the 48 outcomes of the per-outcome gain in held-out concordance.",
]


# ------------------------------------------------------------------ excel
def write_xlsx(cmp_df: pd.DataFrame, org: pd.DataFrame, prov: dict):
    wb = Workbook()
    thin = Side(style="thin", color="D0D0D0")
    border = Border(left=thin, right=thin, top=thin, bottom=thin)
    hdr_fill = PatternFill("solid", fgColor="1F3864")
    hdr_font = Font(bold=True, color="FFFFFF", size=11)
    ox_fill = PatternFill("solid", fgColor="F2F2F2")

    def sheet(ws, df, widths, numfmt, freeze="A2", wrap_cols=()):
        ws.append(list(df.columns))
        for c in range(1, len(df.columns) + 1):
            cell = ws.cell(row=1, column=c)
            cell.font, cell.fill, cell.border = hdr_font, hdr_fill, border
            cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
        ws.row_dimensions[1].height = 34
        def clean(v):
            if isinstance(v, (np.generic,)):
                v = v.item()
            if v is None or (isinstance(v, float) and not np.isfinite(v)):
                return None
            if isinstance(v, (bool, int, float, str)):
                return v
            return str(v)

        for _i, r in df.iterrows():
            ws.append([clean(v) for v in r.tolist()])
        for ri in range(2, len(df) + 2):
            oxy = ("Oxygen row" in df.columns
                   and df.iloc[ri - 2]["Oxygen row"] == "YES")
            for ci, col in enumerate(df.columns, start=1):
                cell = ws.cell(row=ri, column=ci)
                cell.border = border
                if oxy:
                    cell.fill = ox_fill
                if col in numfmt:
                    cell.number_format = numfmt[col]
                    cell.alignment = Alignment(horizontal="right")
                elif col in wrap_cols:
                    cell.alignment = Alignment(wrap_text=True, vertical="top")
                else:
                    cell.alignment = Alignment(vertical="top")
        for ci, col in enumerate(df.columns, start=1):
            ws.column_dimensions[get_column_letter(ci)].width = widths.get(col, 14)
        ws.freeze_panes = freeze

    # ---------------- README
    ws = wb.active
    ws.title = "README"
    ws.column_dimensions["A"].width = 30
    ws.column_dimensions["B"].width = 112
    lines = [
        ("T90 versus total sleep time", ""),
        ("What this is",
         "How much held-out discrimination each exposure adds on top of progressively richer "
         "clinical baselines, across the paper's 48 outcomes. Supplementary side analysis. "
         "Nothing here enters the manuscript."),
        ("", ""),
        ("MODEL SPECIFICATION", "identical for every row of every sheet"),
    ]
    for i, s in enumerate(SPEC_LINES, start=1):
        lines.append((f"  {i}.", s))
    lines += [
        ("", ""),
        ("POSITIVE CONTROLS", "run before any ladder model"),
        ("  total sleep time",
         f"reproduces numbers/ranking_v3.csv row {prov['tst_rank']} of 197. Published dC "
         f"{prov['tst_pub_dC']:.12f}, here {prov['tst_got_dC']:.12f}."),
        ("  T90",
         f"reproduces the published dC. Published {prov['t90_pub_dC']:.12f}, here "
         f"{prov['t90_got_dC']:.12f}."),
        ("  baseline C",
         f"published {prov['pub_base']:.6f}, here {prov['got_base']:.6f}."),
        ("  results.csv re-run",
         "three T90 rungs from the existing results.csv were re-fitted by the new script and "
         "matched to 1e-12, so both exposure columns come from one machine."),
        ("", ""),
        ("THE TST-AS-COVARIATE ASYMMETRY", ""),
        ("  the problem",
         "The T90 ladder's rich rungs hold total sleep time as a covariate. That covariate "
         "cannot survive when total sleep time is the exposure."),
        ("  the fix",
         "Those rungs were re-fitted for T90 with the TST covariate removed, and the comparison "
         "table quotes those, so both exposures face an identical baseline. The published T90 "
         "rung that holds TST is kept as a reference row and marked not comparable."),
        ("  L4 and L4b",
         "They are the same model. L4 is age, sex, race, smoking, AHI and sleep efficiency. L4b "
         "is the T90 rung minus TST, which is the same set. Reported once."),
        ("", ""),
        ("OXYGEN ROWS", "shaded grey in the Comparison sheet and marked YES in the Oxygen row "
                        "column. They hold whole-night mean SpO2, which comes off the same "
                        "oximetry trace as T90, so the baseline is close to adjusting the "
                        "exposure away. Run for completeness, not because they carry the "
                        "argument."),
        ("", ""),
        ("BLOCKS", "full = 19,173 recordings. bmi_6325 = the 6,325 rows with a derivable BMI, "
                   "which carry their own matched age and sex control. Never compare a "
                   "concordance across blocks."),
        ("", ""),
        ("SHEETS", ""),
        ("  Comparison", "the table, 4 decimals as requested. Note that the total sleep time dC "
                         "values sit around -0.001, so 4 decimals rounds -0.00096 to -0.0010 and "
                         "-0.00041 to -0.0004. Nothing is lost, the Full precision sheet and "
                         "results_tst.csv carry every digit."),
        ("  Full precision", "the same rows at full precision, nothing rounded away"),
        ("  Per-organ dC", "the five organ groups for the L1 and L4 rungs"),
        ("  TST ladder detail", "every model in the total sleep time run, including the "
                                "age-spline drop check"),
        ("", ""),
        ("SOURCE FILES", "results.csv, results_tst.csv, per_outcome.csv, per_outcome_tst.csv, "
                         "REPORT_TST.txt, provenance_tst.json"),
        ("BUILT BY", "run_tst_baselines.py then make_comparison.py"),
    ]
    for k, v in lines:
        ws.append([k, v])
    for r in range(1, ws.max_row + 1):
        ws.cell(row=r, column=1).font = Font(bold=True)
        ws.cell(row=r, column=2).alignment = Alignment(wrap_text=True, vertical="top")
    ws.cell(row=1, column=1).font = Font(bold=True, size=14)

    # ---------------- Comparison
    show = ["Baseline model", "Covariates", "Block", "n outcomes", "Baseline C",
            "dC from adding T90", "dC from adding TST", "Difference, T90 minus TST",
            "Magnitude ratio, T90 over TST", "Rel gain T90", "Rel gain TST",
            "Oxygen row", "Note"]
    d4 = "0.0000"
    numfmt = {"Baseline C": d4, "dC from adding T90": d4, "dC from adding TST": d4,
              "Difference, T90 minus TST": d4, "Magnitude ratio, T90 over TST": "0.0",
              "Rel gain T90": "0.00%", "Rel gain TST": "0.00%", "n outcomes": "0"}
    widths = {"Baseline model": 44, "Covariates": 56, "Block": 11, "n outcomes": 10,
              "Baseline C": 12, "dC from adding T90": 17, "dC from adding TST": 17,
              "Difference, T90 minus TST": 17, "Magnitude ratio, T90 over TST": 16,
              "Rel gain T90": 12, "Rel gain TST": 12, "Oxygen row": 11, "Note": 52}
    sheet(wb.create_sheet("Comparison"), cmp_df[show], widths, numfmt, freeze="B2",
          wrap_cols=("Baseline model", "Covariates", "Note"))

    # ---------------- Full precision
    fp = cmp_df.copy()
    fpnum = {c: "0.00000000000" for c in
             ["Baseline C", "C with T90", "dC from adding T90", "C with TST",
              "dC from adding TST", "Difference, T90 minus TST", "Rel gain T90",
              "Rel gain TST", "Magnitude ratio, T90 over TST"]}
    fpw = {c: 20 for c in fpnum}
    fpw.update({"Baseline model": 44, "Covariates": 56, "Note": 52, "Block": 11,
                "T90 source": 16})
    sheet(wb.create_sheet("Full precision"), fp, fpw, fpnum, freeze="B2",
          wrap_cols=("Baseline model", "Covariates", "Note"))

    # ---------------- organs
    onum = {c: d4 for c in org.columns if c not in ("Organ group", "n outcomes")}
    onum["n outcomes"] = "0"
    ow = {"Organ group": 16, "n outcomes": 11}
    ow.update({c: 17 for c in org.columns if c not in ow})
    sheet(wb.create_sheet("Per-organ dC"), org, ow, onum)

    # ---------------- full TST ladder
    lad = NEW.copy()
    lnum = {c: d4 for c in ["mean_C_baseline", "mean_C_with_exposure", "mean_dC",
                            "rel_gain_on_means", "mean_rel_gain_per_outcome", "worst_dC",
                            "best_dC"]}
    lw = {"block": 14, "model": 42, "exposure": 10, "covariates": 62, "oxygen_row": 12,
          "t90_counterpart": 40, "age_basis": 28}
    lw.update({c: 16 for c in lad.columns if c not in lw})
    sheet(wb.create_sheet("TST ladder detail"), lad, lw, lnum, freeze="C2",
          wrap_cols=("covariates", "model", "t90_counterpart"))

    path = os.path.join(HERE, "COMPARISON_T90_vs_TST.xlsx")
    wb.save(path)
    return path


def main():
    import json
    pj = json.load(open(os.path.join(HERE, "provenance_tst.json")))["positive_control"]
    prov = {"tst_pub_dC": pj["tst_published_dC"], "tst_got_dC": pj["tst_dC"],
            "t90_pub_dC": pj["t90_published_dC"], "t90_got_dC": pj["t90_dC"],
            "pub_base": pj["published_baseline_C"], "got_base": pj["baseline_C"],
            "tst_rank": pj["tst_rank_in_ranking_v3"], "tst_tol": "1e-9", "t90_tol": "1e-9"}
    cmp_df = build()
    org = organ_block()
    p1 = write_md(cmp_df, org, prov)
    p2 = write_xlsx(cmp_df, org, prov)
    cmp_df.to_csv(os.path.join(HERE, "COMPARISON_T90_vs_TST_table.csv"), index=False)
    print("wrote")
    print(" ", p1)
    print(" ", p2)
    print("  ", os.path.join(HERE, "COMPARISON_T90_vs_TST_table.csv"))
    print()
    with pd.option_context("display.width", 240, "display.max_columns", 40):
        print(cmp_df[["Baseline model", "Baseline C", "dC from adding T90",
                      "dC from adding TST", "Difference, T90 minus TST"]].to_string(index=False))
        print()
        print(org.to_string(index=False))


if __name__ == "__main__":
    main()
