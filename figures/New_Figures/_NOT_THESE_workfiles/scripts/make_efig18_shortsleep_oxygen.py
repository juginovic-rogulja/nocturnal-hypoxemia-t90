"""
eFigure 18. Short sleepers more often have low nocturnal oxygen.

Alen, 2026-08-13: this result matters for the Discussion, so it needs the paper's own styling
rather than the standalone look it was prototyped in.

The point it makes: a study of sleep duration that does not measure oxygen is not comparing
short sleepers with long sleepers on equal footing. Among otherwise healthy people, half again
as many short sleepers are hypoxemic. Any disease that low oxygen causes will therefore attach
itself to short sleep in a study blind to oxygen.

Four panels because there are two populations and two ways of measuring sleep, and the reader
should see that the home report, which was written weeks before the study by a clinician, gives
the same answer as the laboratory night that produced the oxygen value itself.

This is a PREVALENCE ratio, measured on one night. Nobody is followed forward on this sheet.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../..")))  # repository root, where paths.py lives
import paths
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import norm

_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _here)
from splitstyle import (rc, save, panel_label, headline, BLUE, GREY_PALE, INK,  # noqa: E402
                        GREY, TITLE, LABEL, TICK, ANNOT, plt)

ROOT = paths.FIGURE_ROOT
sys.path.insert(0, f"{paths.ANALYSIS_DIR}")
from cohort_spec import apply_cohort  # noqa: E402

CARD = ["copd2", "asthma", "resp_failure", "obesity_hypovent", "pulm_htn", "hf", "ihd", "mi",
        "afib", "stroke_any"]
REN = ["ckd", "cirrhosis", "cancer_any"]
MET = ["diabetes", "obesity", "nafld"]

b = apply_cohort(pd.read_parquet(f"{paths.TABLES_DIR}/t90_final.parquet"))
# 2026-08-13, second pass. The v1 extraction behind the old reported-sleep panels carried a
# 47% record-level error (read-back validated): patient handouts, PAP-compliance fields, lab
# echoes and past states scored as the patient's sleep. The v2 extraction vetoes those by
# sentence context (validated over four read-back rounds to ~13% residual error, corpus-wide
# EC2 re-scan 2026-08-13) and corrects hr+min and "between X and Y" values into hours_v2.
# Only clean records (veto == "") are used here. The earlier single-signature boilerplate
# filter is superseded by this, and the v1 downward bias (machine hours read as sleep) is
# gone: cleaned home medians run ~0.8 h longer.
_h = pd.read_parquet(f"{paths.SV_ROOT}/"
                     "habitual_sleep_v2/habitual_records_v2_full.parquet")
_u = _h[(_h.veto == "") & _h.tier.isin(["A", "B", "C"]) & _h.hours_v2.between(0.5, 14)]
_u = _u[~_u.lab_echo.fillna(False).astype(bool)]
_u = _u.rename(columns={"person_id": "BDSPPatientID"})
_per = (_u.groupby("BDSPPatientID").hours_v2.median() * 60).rename("hab").reset_index()
b = b.merge(_per, on="BDSPPatientID", how="left")
healthy = pd.Series(True, index=b.index)
for k in CARD + REN + MET:
    healthy &= (b[f"{k}_prevalent"] != 1)


def cell(frame, cut, col):
    f = frame[frame[col].notna()]
    s, r = f[f[col] < 300], f[f[col] >= 420]
    k1, n1 = int((s.spo2_pct_below_90 > cut).sum()), len(s)
    k0, n0 = int((r.spo2_pct_below_90 > cut).sum()), len(r)
    rr = (k1 / n1) / (k0 / n0)
    se = np.sqrt(1 / k1 - 1 / n1 + 1 / k0 - 1 / n0)
    p = float(2 * (1 - norm.cdf(abs(np.log(rr) / se))))
    return 100 * k1 / n1, 100 * k0 / n0, n1, n0, rr, rr * np.exp(-1.96 * se), \
        rr * np.exp(1.96 * se), p


def pfmt(p):
    return "P < .001" if p < 0.001 else "P = " + f"{p:.3f}".lstrip("0")


PANELS = [("Everyone, sleep measured in the laboratory", b, "TST_min", "A"),
          ("Free of organ and metabolic disease, measured", b[healthy], "TST_min", "B"),
          ("Everyone, sleep reported at home", b, "hab", "C"),
          ("Free of organ and metabolic disease, reported", b[healthy], "hab", "D")]

W, H = 9.6, 10.4
with plt.rc_context(rc()):
    fig = plt.figure(figsize=(W, H))
    POS = [[0.095, 0.620, 0.365, 0.255], [0.590, 0.620, 0.365, 0.255],
           [0.095, 0.145, 0.365, 0.255], [0.590, 0.145, 0.365, 0.255]]
    for (title, frame, col, letter), pos in zip(PANELS, POS):
        ax = fig.add_axes(pos)
        got = [cell(frame, c, col) for c in (5, 10)]
        xs = np.arange(2)
        ax.bar(xs - 0.19, [g[0] for g in got], 0.36, color=BLUE, linewidth=0,
               label="Sleeps under 5 h")
        ax.bar(xs + 0.19, [g[1] for g in got], 0.36, color=GREY_PALE, linewidth=0,
               label="Sleeps 7 h or more")
        for i, g in enumerate(got):
            ax.text(i - 0.19, g[0] + 0.7, f"{g[0]:.1f}", ha="center", fontsize=ANNOT,
                    fontweight="bold", color=INK)
            ax.text(i + 0.19, g[1] + 0.7, f"{g[1]:.1f}", ha="center", fontsize=ANNOT, color=INK)
            ax.text(i, -13.0, f"{g[4]:.2f} times as many\n({g[5]:.2f}-{g[6]:.2f})\n{pfmt(g[7])}",
                    ha="center", fontsize=ANNOT - 0.5, color=INK, linespacing=1.5)
        ax.set_xticks(xs)
        ax.set_xticklabels(["Below 90% for more\nthan 5% of the night",
                            "Below 90% for more\nthan 10% of the night"], fontsize=TICK)
        ax.set_ylim(0, 33)
        ax.set_ylabel("Percent of the group", fontsize=LABEL)
        for s_ in ("top", "right"):
            ax.spines[s_].set_visible(False)
        ax.set_title(f"{title}\n{got[0][2]:,} against {got[0][3]:,}", fontsize=LABEL,
                     loc="left", pad=8, linespacing=1.5)
        panel_label(ax, letter, dx=-0.155, dy=1.13)
        if letter == "A":
            ax.legend(loc="upper right", fontsize=ANNOT, frameon=False, handletextpad=0.6,
                      labelspacing=0.4)
    headline(fig, "Short sleepers more often have low nocturnal oxygen")
    save(fig, "eFigure18_short_sleep_and_oxygen", kind="Supplementary")

for title, frame, col, letter in PANELS:
    for c in (5, 10):
        g = cell(frame, c, col)
        print(f"  {letter} {title[:38]:38} T90>{c:>2}%  {g[0]:5.1f} vs {g[1]:5.1f}  "
              f"{g[4]:.2f} ({g[5]:.2f}-{g[6]:.2f})")
