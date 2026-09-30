#!/usr/bin/env python3
"""ED_Fig03, V14 lane L2 (round 37): shared paths and the generator's rules (new module, the R30 scripts carried them inline).
Base = the V13 sheet (byte-identical to the V7 sheet, the round-27 lane C build). Numbers = numbers/results_v2.json ["graded"] (step 112,
v8.1: 58 outcomes, the five v8.1 controls carry negative_control = true) and numbers/negcontrols_v4_panel.csv (step 251, the panel of
five). graded_newcontrols_v1.json (step 158) is no longer needed: every control is in the graded block. OLD side of the delta and
the positive control: OLD = V8_RECALC_2026-09-12/compare/old_numbers/numbers/results_v2.json (the v7 snapshot, OLD side only, never a live source).
Row set: the August generator's rule (make_efig02_merged_polish.py lines 146 to 157, 230 to 244) APPLIED, which closes GAP-11 for this
sheet (the v7 lane kept the base sheet's rows because figures_alen/AlenETable_graded_all.csv, the static v7 row list, had no writer)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import csv, json, os, sys
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND37_2026-09-14/figures/_common")
from v14lib import *
from matplotlib.ticker import MaxNLocator
from statsmodels.stats.multitest import multipletests

LANE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORK, VER = os.path.join(LANE, "work"), os.path.join(LANE, "verify"); CROPS = os.path.join(VER, "crops")
for _d in (WORK, VER, CROPS): os.makedirs(_d, exist_ok=True)
BASE = f"{V13}/ED_Fig03.pdf"
DATA = f"{NUM}/results_v2.json"; PANEL_CSV = f"{NUM}/negcontrols_v4_panel.csv"
OLD_DATA = f"{OLD_V7}/numbers/results_v2.json"; OLD_PANEL = f"{OLD_V7}/numbers/negcontrols_v4_panel.csv"
BAND_KEYS = ["0-1%", "1-5%", "5-10%", ">10%"]; BANDS = ["0-1", "1-5", "5-10", ">10"]
CIRC = {"Insomnia", "Restless legs syndrome", "Nocturia", "Epilepsy"}


def graded(path):
    return json.load(open(hydrated(path)))["graded"]


def controls_of(panel_csv):
    return [r["condition"] for r in csv.DictReader(open(hydrated(panel_csv)))]


def august_rows(out, ctrl):
    """panel a = the 15 largest top-band (>10%) hazard ratios among graded conditions that are not controls and not circular,
    descending; panel b = the five controls, descending by the top-band hazard ratio."""
    cond = [k for k, v in out.items() if not v["negative_control"] and k not in CIRC and k not in ctrl]
    a = sorted(cond, key=lambda k: -out[k][">10%"]["hr"])[:15]
    b = sorted(ctrl, key=lambda k: -out[k][">10%"]["hr"])
    return a, b


def ylim_ticks(v):
    """bottom min(0.72, 0.9 x lowest CI), top max(1.05, 1.22 x highest CI); MaxNLocator(nbins=4, steps=[1,2,2.5,5,10]) on the limits,
    ticks at or above bottom + 6 percent of the span, %g labels."""
    lo = [v[b]["lo"] for b in BAND_KEYS[1:]]; hi = [v[b]["hi"] for b in BAND_KEYS[1:]]
    ymin = min(0.72, 0.9 * min(lo)); ymax = max(1.05, 1.22 * max(hi)); span = ymax - ymin
    cand = MaxNLocator(nbins=4, steps=[1, 2, 2.5, 5, 10]).tick_values(ymin, ymax)
    return ymin, ymax, [f"{float(t):g}" for t in cand if ymin + 0.06 * span <= t <= ymax]


def bh_top(out, keys=None):
    keys = list(out) if keys is None else keys
    return dict(zip(keys, (float(q) for q in multipletests([out[k][">10%"]["p"] for k in keys], method="fdr_bh")[1])))


def wrap_title(text, tw, limit):
    if tw(text) <= limit: return [text]
    words, lines, cur = text.split(" "), [], ""
    for w in words:
        trial = (cur + " " + w).strip()
        if cur and tw(trial) > limit: lines.append(cur); cur = w
        else: cur = trial
    lines.append(cur)
    return lines
