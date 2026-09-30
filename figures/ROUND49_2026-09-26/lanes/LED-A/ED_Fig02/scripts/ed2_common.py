#!/usr/bin/env python3
"""ED_Fig02, V14 lane L2 (round 37): shared paths, colours, the data derivation (the generator's row rule) and helpers.
Repointed copy of ROUND30_2026-09-08/V7_L2_T90/ED_Fig02/scripts/ed2_common.py (byte copy beside as ed2_common_PRE_V8_1.py).
V14: the base is the V13 sheet; lag_ladder.csv/json (step 130), causal_tests_all.csv (step 250) and, since round v8.2 (2026-09-15
evening), sleep_unadjusted_v1.json (step 125, re-extracted group P input, sidecar 19:29) are v8 outputs under the sidecar gate
v14lib.sidecar_v8; the v8.2 file must carry an output mtime after 2026-09-15 19:18 (V8_2_CUT) and holds every control, so panel c
draws all five controls (no NE). The controls come from numbers/disease_definitions.py (v8.1) for the live build and from the lag
ladder's own flags for the OLD snapshot (whose sleep file is the v7 snapshot's own). Edited copy backed up as ed2_common_PRE_V8_2.py."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import hashlib, json, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND37_2026-09-14/figures/_common")
import v14lib
from v14lib import hydrated, sha256, sidecar_v8, halfup
LANE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORK, VER = os.path.join(LANE, "work"), os.path.join(LANE, "verify"); CROPS = os.path.join(VER, "crops")
for _d in (WORK, VER, CROPS): os.makedirs(_d, exist_ok=True)
T90 = paths.FIGURE_ROOT
NUM = f"{paths.NUMBERS_DIR}"; OLD = f"{paths.V8_ROOT}/compare/old_numbers/numbers"     # v7 snapshot: OLD side of the delta and the positive control only
BASE = f"{paths.FIGURE_ROOT}/ROUND36_2026-09-10/L9_assembly/NEW_FINAL_SET_V13/ED_Fig02.pdf"               # the V13 design baseline
LOCK = f"{paths.FIGURE_ROOT}/ROUND30_2026-09-08/_render.lock"
ARIAL = paths.FONT_ARIAL; ARIAL_B = paths.FONT_ARIAL_BOLD
INK, BLUE, GREY, PALE, BAND = "#1a1d21", "#0288d1", "#8a9099", "#ccd1d6", "#eef0f1"
FILES = ("lag_ladder.csv", "lag_ladder.json", "causal_tests_all.csv", "sleep_unadjusted_v1.json")
V8_FILES = ("lag_ladder.csv", "lag_ladder.json", "causal_tests_all.csv", "sleep_unadjusted_v1.json")     # all four v8 outputs since v8.2
GROUP_P = ()                                                                                            # v8.2: nothing kept from v7
V8_2_CUT = "2026-09-15 19:18:00"                                                                        # the v8.2 re-extraction: sources dated after this
V8_2_FILES = ("sleep_unadjusted_v1.json", "causal_tests_all.csv")                                       # rewritten by the v8.2 chain (19:29, 23:22)


def sidecar(p):
    """Every numbers file: the v8 gate (v14lib.sidecar_v8: inputs name data_frozen_v8_2026-09 and the sidecar sha256 equals the file read now).
    The files the v8.2 chain rewrote must carry an output mtime after V8_2_CUT (refuse the older v7 group P copy)."""
    d = sidecar_v8(p); d["group_P"] = False
    if os.path.basename(p) in V8_2_FILES:
        assert d["output_mtime"] and d["output_mtime"] >= V8_2_CUT, ("NOT a v8.2 output (older than 2026-09-15 19:18): refused", p, d["output_mtime"])
        d["v8_2_gate"] = f"output mtime {d['output_mtime']} after {V8_2_CUT}"
    return d


def controls_of(root):
    """The control panel: for the live numbers the v8.1 definitions order; for the OLD snapshot the ladder's own flags in the V13 sheet's order."""
    LAG = pd.read_csv(hydrated(f"{root}/lag_ladder.csv")); lag0 = LAG[LAG.lag_years == 0]
    flagged = list(lag0[lag0.negative_control].condition)
    if root == NUM:
        sys.path.insert(0, paths.ANALYSIS_DIR)
        from disease_definitions import DISEASES, NEGATIVE_CONTROLS
        c = [DISEASES[k][0] for k in NEGATIVE_CONTROLS]; assert set(c) == set(flagged), (c, flagged); return c
    order = ["Back pain", "Cataract", "Glaucoma", "Contact dermatitis", "Hemorrhoids"]    # the V13 sheet's control rows (v7 definitions order)
    assert set(order) == set(flagged), (order, flagged); return order


def derive(root, tag):
    """The generator's row rule and lanes (eFigure07_polish.py lines 182 to 224, 437 to 509): 14 conditions ranked by the lag-0 hazard
    ratio among conditions present on all three tests, not circular, not controls, not above 2.3, then the controls. Panel c (the
    sleep-period test) comes from the same root's sleep file; a control absent from it would get c = None (drawn as NE): under v8.2
    the file holds every control, so the live build asserts c_missing == []."""
    CONTROLS = controls_of(root)
    LAG = pd.read_csv(hydrated(f"{root}/lag_ladder.csv")); CT = pd.read_csv(hydrated(f"{root}/causal_tests_all.csv"))
    WS = json.load(open(hydrated(f"{root}/sleep_unadjusted_v1.json"))); LJ = json.load(open(hydrated(f"{root}/lag_ladder.json")))   # v8.2: the sleep file of the same root (v7 snapshot for OLD, the 19:29 v8.2 file for NUM)
    lag0 = LAG[(LAG.lag_years == 0) & (~LAG.circular)]
    avail = set(lag0.condition) & set(CT.condition) & set(WS["outcomes"])
    cands = [c for c in avail if c not in CONTROLS and not bool(lag0[lag0.condition == c].negative_control.iloc[0])]
    OFFSCALE = [c for c in cands if float(lag0[lag0.condition == c].hr.iloc[0]) > 2.3]
    cands = [c for c in cands if c not in OFFSCALE]
    ranked = sorted(cands, key=lambda c: -float(lag0[lag0.condition == c].hr.iloc[0]))
    CONDS = ranked[:14]; ROWS = CONDS + CONTROLS
    excluded_by_ws = sorted(c for c in set(lag0[~lag0.negative_control].condition) - set(WS["outcomes"]))
    assert len(ROWS) == 19 and OFFSCALE == ["Obesity hypoventilation"], (len(ROWS), OFFSCALE)
    rows = []
    for c in ROWS:
        g = LAG[LAG.condition == c].set_index("lag_years"); r = CT[CT.condition == c].iloc[0]; v = WS["outcomes"].get(c)
        assert np.isfinite(r.bmi_adj_hr), c
        key = LAG[LAG.condition == c].key.iloc[0]; pj = LJ["per_condition"][key]
        for lag in (0, 1, 2, 5):
            if bool(g.loc[lag, "estimable"]):
                assert abs(pj["hr_by_lag"][str(lag)] - float(g.loc[lag, "hr"])) < 1e-9, (c, lag)
                assert abs(pj["ci_by_lag"][str(lag)][0] - float(g.loc[lag, "lo"])) < 1e-9 and abs(pj["ci_by_lag"][str(lag)][1] - float(g.loc[lag, "hi"])) < 1e-9, (c, lag)
        est5 = bool(g.loc[5, "estimable"]); inf5 = bool(g.loc[5, "lag5_informative"])
        rows.append(dict(condition=c, key=key, control=c in CONTROLS,
            a=dict(before=float(g.loc[0, "hr"]), after=float(g.loc[2, "hr"]), lo=float(g.loc[2, "lo"]), hi=float(g.loc[2, "hi"]),
                   y1=float(g.loc[1, "hr"]), y1_text=halfup(float(g.loc[1, "hr"]), 2), y5=(float(g.loc[5, "hr"]) if est5 else None),
                   y5_text=(halfup(float(g.loc[5, "hr"]), 2) if est5 else "NE"), y5_estimable=est5, y5_informative=inf5, y5_grey=(not est5) or (not inf5),
                   events={f"lag{l}": int(g.loc[l, "events"]) for l in (0, 1, 2, 5)}),
            b=dict(before=float(r.bmi_sub_hr), after=float(r.bmi_adj_hr), lo=float(r.bmi_adj_lo), hi=float(r.bmi_adj_hi), events=int(r.bmi_events)),
            c=(dict(before=float(v["sleep_alone"]["hr"]), after=float(v["sleep_adjusted_for_awake"]["hr"]), lo=float(v["sleep_adjusted_for_awake"]["lo"]),
                    hi=float(v["sleep_adjusted_for_awake"]["hi"]), events=int(v["events"])) if v is not None else None)))
    vals = [x for rw in rows for p in "abc" if rw[p] is not None for x in (rw[p]["before"], rw[p]["after"], rw[p]["lo"], rw[p]["hi"])]
    out = dict(tag=tag, root=root, rows=rows, conds=CONDS, controls=CONTROLS, offscale=OFFSCALE, ranked_candidates=ranked[:18],
               excluded_by_sleep_file=excluded_by_ws, c_missing=[rw["condition"] for rw in rows if rw["c"] is None],
               lag0_hr_ranked={c: float(lag0[lag0.condition == c].hr.iloc[0]) for c in ranked[:18]},
               generator_xlim=[min(vals) * 0.96, max(vals) * 1.04], drawn_min=min(vals), drawn_max=max(vals),
               sha256={f: sha256(f"{root}/{f}") for f in V8_FILES},
               n_ne=sum(1 for rw in rows if not rw["a"]["y5_estimable"]), n_grey5=sum(1 for rw in rows if rw["a"]["y5_estimable"] and not rw["a"]["y5_informative"]),
               ws_n=WS["n"], lag_cohort_n=int(LJ["cohort_n"]))
    if tag == "new":
        out["sidecars"] = {f: sidecar(f"{root}/{f}") for f in V8_FILES}
        assert out["c_missing"] == [], ("v8.2: the sleep file must hold every drawn control", out["c_missing"])
    return out


def load_geometry():
    return json.load(open(hydrated(os.path.join(WORK, "base_geometry.json"))))


def span_at(G, text, size=None, xmin=None, xmax=None, ymin=None, ymax=None, nth=0):
    """A visible base span by text (and optional size and region), sorted by y then x."""
    c = [s for s in G["visible_spans"] if s["text"] == text and (size is None or abs(s["size"] - size) < 0.06)
         and (xmin is None or s["origin"][0] >= xmin) and (xmax is None or s["origin"][0] <= xmax)
         and (ymin is None or s["origin"][1] >= ymin) and (ymax is None or s["origin"][1] <= ymax)]
    c = sorted(c, key=lambda s: (s["origin"][1], s["origin"][0]))
    assert len(c) > nth, (text, size, len(c))
    return c[nth]


def xmap(panel):
    return lambda v: panel["a"] + panel["perln"] * float(np.log(v))
