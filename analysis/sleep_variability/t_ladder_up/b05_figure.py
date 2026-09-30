"""
The ladder figure: dC and per-SD FDR-significant count against the SpO2 threshold.
Working-sheet quality by instruction, NOT the publication pipeline. Every plotted
value is read from LADDER_SUMMARY.csv, nothing is typed in.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../..")))  # repository root, where paths.py lives
import paths
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import sys as _sys_v8; _sys_v8.path.insert(0, paths.ANALYSIS_DIR)  # v8 sweep 2026-09-12
from cohort_spec import N_RANKED_OUTCOMES  # v8 sweep 2026-09-12

HERE = os.path.dirname(os.path.abspath(__file__))
S = pd.read_csv(os.path.join(HERE, "LADDER_SUMMARY.csv"))

prim = S[S.threshold.isin(["T80", "T85", "T88", "T90", "T92", "T95"])].copy()
prim["x"] = prim.threshold.str[1:].astype(int)
prim = prim.sort_values("x")
sec = S[S.threshold.isin(["T88_rederiv", "T90_rederiv"])].copy()
sec["x"] = sec.threshold.str.extract(r"T(\d+)")[0].astype(int)

DARK = "#1f3b64"          # dark blue
GREY = "#8a8f98"

fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.8))

ax = axes[0]
ax.plot(prim.x, prim.dC_mean, "-o", color=DARK, ms=6, lw=1.8, zorder=3,
        label="ladder (frozen T88/T90, re-derived otherwise)")
ax.plot(sec.x, sec.dC_mean, "o", mfc="none", mec=GREY, ms=8, mew=1.6, zorder=2,
        label="same-code re-derived T88/T90")
t90 = prim[prim.threshold == "T90"]
ax.scatter(t90.x, t90.dC_mean, s=150, facecolors="none", edgecolors=DARK,
           lw=1.6, zorder=4)
ax.annotate("published T90\n+0.020339", (float(t90.x.iloc[0]), float(t90.dC_mean.iloc[0])),
            textcoords="offset points", xytext=(8, -26), fontsize=8, color=DARK)
ax.set_xlabel("SpO2 threshold (%)")
ax.set_ylabel(f"dC, mean discrimination gain over {N_RANKED_OUTCOMES} outcomes")
ax.set_title("A  Discrimination gain by threshold", loc="left", fontsize=10)
ax.legend(fontsize=7, frameon=False, loc="lower center")

ax = axes[1]
ax.plot(prim.x, prim.n_fdr_sig_of_48_persd, "-o", color=DARK, ms=6, lw=1.8, zorder=3)
ax.plot(sec.x, sec.n_fdr_sig_of_48_persd, "o", mfc="none", mec=GREY, ms=8, mew=1.6,
        zorder=2)
ax2 = ax.twinx()
ax2.plot(prim.x, prim.median_persd_hr, "--s", color=GREY, ms=4, lw=1.2, alpha=0.9)
ax2.set_ylabel("median per-SD HR (grey, dashed)", color=GREY, fontsize=8)
ax2.tick_params(axis="y", colors=GREY, labelsize=8)
ax2.spines["top"].set_visible(False)
ax.set_xlabel("SpO2 threshold (%)")
ax.set_ylabel(f"outcomes FDR-significant per SD (of {N_RANKED_OUTCOMES})")
ax.set_ylim(0, N_RANKED_OUTCOMES)
ax.set_title("B  Breadth of association by threshold", loc="left", fontsize=10)

for ax in axes:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xticks([80, 85, 88, 90, 92, 95])
fig.suptitle("Oxygen threshold ladder, T80 to T95 (side analysis, not for the manuscript)",
             fontsize=10, y=1.00)
fig.tight_layout()
out = os.path.join(HERE, "ladder_figure.pdf")
fig.savefig(out, dpi=300, bbox_inches="tight")
print(f"wrote {out}")
