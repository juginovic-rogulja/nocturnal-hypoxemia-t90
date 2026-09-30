#!$T90_PY
"""Widths of every gutter label at the round-40 size and at +1 pt against its plot edge (matplotlib, the builders' own rc)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import json, sys
import pandas as pd
L = f"{paths.FIGURE_ROOT}/ROUND50_2026-09-26/lanes/LSUPP-C"
sys.path.insert(0, f"{L}/L5b_scripts")
from splitstyle import rc, plt
NUM = paths.NUMBERS_DIR
R = dict(rc()); R.update({"savefig.bbox": "standard"})
def widths(strings, size, weight="normal"):
    with plt.rc_context(R):
        fig = plt.figure(figsize=(7.2, 2.0)); fig.canvas.draw(); r = fig.canvas.get_renderer(); out = {}
        for s in strings:
            t = fig.text(0.5, 0.5, s, fontsize=size, fontweight=weight); out[s] = t.get_window_extent(renderer=r).width / fig.dpi; t.remove()
        plt.close(fig)
    return out
def report(name, strings, gut, edge, safety, size_old, size_new):
    wo, wn = widths(strings, size_old), widths(strings, size_new)
    print(f"\n== {name}: gutter starts {gut:.3f} in, plot edge {edge:.3f} in (gate safety {safety:.2f} in), sizes {size_old} -> {size_new}")
    rows = sorted(strings, key=lambda s: -wn[s])
    for s in rows[:6]:
        print(f"   {s!r:32} old {wo[s]:.3f} in (end {gut + wo[s]:.3f}, {(edge - gut - wo[s]) * 72:+.1f} pt to the edge)  new {wn[s]:.3f} in (end {gut + wn[s]:.3f}, {(edge - gut - wn[s]) * 72:+.1f} pt to the edge)")
    n_gate = sum(1 for s in strings if gut + wn[s] > edge - safety); n_edge = sum(1 for s in strings if gut + wn[s] > edge)
    print(f"   at +1: {n_gate} of {len(strings)} fail the gate, {n_edge} cross the plot edge")
# Supp 20
D = json.load(open(f"{L}/Supp_Fig20/work/eFigure12_drawn_values_v8_2.json"))
report("Supp 20 a (12 outcomes)", [r["outcome"] for r in D["A"]], 0.18, 1.84, 0.10, 10.0, 11.0)
report("Supp 20 b (21 labels)", [r["label"] for r in D["B"]], 0.18, 1.89, 0.10, 10.0, 11.0)
report("Supp 20 c (21 diseases)", [r["disease"] for r in D["C"]], 0.18, 1.89, 0.10, 10.0, 11.0)
report("Supp 20 d (8 windows)", ["Whole recording", "Wake", "Whole sleep", "NREM", "N1", "N2", "N3", "REM"], 0.18, 1.30, 0.10, 10.0, 11.0)
# Supp 13
P = json.load(open(f"{NUM}/nonresponder_combined_v2.json"))
ORDER = [k for k, v in P["failure_rate_by_condition"].items() if v["fdr_significant"]]
report("Supp 13 a (21 rows)", ORDER, 0.18, 1.92, 0.06, 10.0, 11.0)
ROWS_C = sorted(P["independent_conditions"], key=lambda k: -P["independent_conditions"][k]["or"]) + ["Age, per year", "Male sex", "Black race"]
report("Supp 13 c (8 rows)", ROWS_C, 0.18, 1.92, 0.06, 10.0, 11.0)
report("Supp 13 d (wrapped model names)", ["Pretreatment oxygen\nalone", "and demographics", "and the\ncardiopulmonary count", "and the independent\nconditions"], 0.18, 1.74, 0.06, 10.0, 11.0)
# Supp 13 printed column against the canvas edge and the 4 mm margin
import re
from decimal import Decimal, ROUND_HALF_UP
def _r(v): return str(Decimal(repr(float(v))).quantize(Decimal("0.01"), rounding=ROUND_HALF_UP))
FULL = json.load(open(f"{L}/L5b_work/work/fullprec_nonresponder_v8_1.json")) if __import__("os").path.exists(f"{L}/L5b_work/work/fullprec_nonresponder_v8_1.json") else json.load(open(f"{paths.FIGURE_ROOT}/ROUND40_2026-09-18/lanes/LSUPP/L5b_work/fullprec_nonresponder_v8_1.json"))
colA = [f"{_r(FULL['by_condition'][k]['marginal_or'])} ({_r(FULL['by_condition'][k]['marginal_lo'])}-{_r(FULL['by_condition'][k]['marginal_hi'])}){'*' if P['failure_rate_by_condition'][k]['independent'] else ''}" for k in ORDER]
W_IN = 171 / 25.4; report("Supp 13 a printed column (x 5.38 in, canvas 6.732, 4 mm margin at 6.575)", colA + ["Odds ratio (95% CI)"], 5.38, W_IN - 4.0 / 25.4, 0.0, 9.5, 10.5)
FPC = {**FULL["joint_independent"], "Age, per year": FULL["joint_demographics"]["AgeAtVisit"], "Male sex": FULL["joint_demographics"]["male"], "Black race": FULL["joint_demographics"]["race_black"]}
colC = [f"{_r(FPC[k]['or'])} ({_r(FPC[k]['lo'])}-{_r(FPC[k]['hi'])})" for k in ROWS_C]
report("Supp 13 c printed column (x 5.38 in)", colC + ["Odds ratio (95% CI)"], 5.38, W_IN - 4.0 / 25.4, 0.0, 9.5, 10.5)
# Supp 12 rows and columns
cont = pd.read_csv(f"{paths.SV_ROOT}/pap_residual_continuous/results_continuous.csv")
RR = cont[cont.sig_fdr == True].sort_values("hr_sd", ascending=False)
report("Supp 12 rows (11)", RR.outcome.tolist(), 0.18, 2.05, 0.0, 10.0, 11.0)
DV = json.load(open(f"{L}/L5b_work/polish_out/work/eFigure21_polish_drawn_values.json"))
report("Supp 12 HR column (x 4.07 in, next column at 5.33)", [d["printed_residual"] for d in DV] + [d["printed_delta"] for d in DV] + ["HR (95% CI)"], 4.07, 5.33, 0.0, 9.5, 10.5)
report("Supp 12 q column (x 5.33 in, crop edge 425.64 pt = 5.912 in, gate at 5.842)", [d["printed_residual_q"] for d in DV] + [d["printed_delta_q"] for d in DV], 5.33, 425.64 / 72, 5.0 / 72, 9.5, 10.5)
# Supp 14 rows
DV14 = json.load(open(f"{L}/L5b_work/polish_out/work/Figure5D_polish_drawn_values.json"))
FIGW = 171 / 25.4; report("Supp 14 rows (10), gutter 4.4 mm, axes left 0.286 of 171 mm", [c["condition"] for c in DV14["conditions"]], 4.4 / 25.4, 0.286 * FIGW, 0.0, 10.0, 11.0)
