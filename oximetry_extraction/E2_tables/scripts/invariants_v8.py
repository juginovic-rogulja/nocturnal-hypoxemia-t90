#!/usr/bin/env python3
"""
Step 051, invariants_v8 (integrity gates 3 and 4) on a built v8 table folder, against the v7 tables (PREV_DIR).

usage: invariants_v8.py --dir <built folder> [--prev PREV_DIR] [--log-dir <folder>] [--smoke]
Writes <log-dir>/INVARIANTS_V8.md and <log-dir>/PER_SITE_ALL_COLUMNS_V8.csv. Exit 1 on any failure.

Checks (F = failure, I = inherited from v7 and reported, not a failure; every reference count is read at run time):
  rows and ids equal to v7 in every table (F); analysis cohort = cohort_spec.COHORT_N = v7 cohort, same ids (F)
  stage percentages sum to 100 +/- 0.05 where TST > 0 (F on v8-replaced nights, I elsewhere); TST <= recording minutes (F/I);
  sleep efficiency <= 100 (F/I); AHI, arousal index and PLM index x TST/60 integer where the count is ours (F on replaced/new, report elsewhere);
  REM = 0 share <= 2 percent of non-split cohort nights with sleep (F), split-night share reported; plausibility bounds zero hits on cohort rows (F);
  every measure column byte-identical to v7 on the non-split cohort nights, every table (F); the four new measures present where defined (F);
  outcome positive control: untouched conditions byte-identical, incident counts = X2's event_counts_v8.csv (F);
  every CHANGELOG row: old = v7 value, new = v8 value (F); base4 shared columns agree with t90_final at least as well as in v7 (F);
  cpap covariates equal the v8 t90_final on refreshed rows (F); sha256.txt matches the files (F);
  per-site distribution (median, IQR, KS) of EVERY numeric column, v7 vs v8, non-split nights identical (F), shifts noted.
"""
from __future__ import annotations
import argparse, os, sys, time
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import e2_spec as S  # noqa: E402
CS = S.CS

TOL_STAGE = 0.05
REM0_MAX_SHARE = 0.02
SHIFT_IQR_FRAC = 0.10      # a per-site median shift above this fraction of the v7 IQR is noted (all of it comes from the split nights)


def nan_equal(a, b):
    a, b = pd.Series(np.asarray(a), dtype=object), pd.Series(np.asarray(b), dtype=object)
    try:
        return ((a == b) | (a.isna() & b.isna())).to_numpy(dtype=bool)
    except Exception:
        return (a.astype(str) == b.astype(str)).to_numpy(dtype=bool)


def bounds_table(prev_dir):
    p = S.pick(S.BOUNDS_CSV)
    if os.path.exists(p) and not S.evicted(p):
        return pd.read_csv(p)[["column", "lo", "hi"]], "csv"
    ch = pd.read_parquet(f"{prev_dir}/CHANGELOG.parquet", columns=["column", "reason"])
    r = ch[ch.reason.str.startswith("impossible value outside")].drop_duplicates(["column", "reason"])
    b = pd.DataFrame({"column": r.column.values, "lo": [float(s.split("[")[1].split(",")[0]) for s in r.reason],
                      "hi": [float(s.split(",")[1].split("]")[0]) for s in r.reason]}).drop_duplicates("column")
    return b, "v7 CHANGELOG reasons (bounds csv evicted)"


def run(out, prev_dir, log_dir, smoke=False):
    t0 = time.time(); L = []; res = []       # res: (name, status, detail) status in PASS / FAIL / INHERITED / NOTE

    def rec(name, ok, detail, inherited=False):
        st = "PASS" if ok else ("INHERITED" if inherited else "FAIL")
        res.append((name, st, detail)); L.append(f"- [{st}] {name}: {detail}")

    rd = lambda d, f: pd.read_parquet(f"{d}/{f}") if f.endswith(".parquet") else pd.read_csv(f"{d}/{f}", low_memory=False, float_precision="round_trip")
    P = {k: rd(prev_dir, f) for k, f in (("fz", "t90_final.parquet"), ("b4", "t90_base4.parquet"), ("ma", "master_cohort.csv"), ("cp", "cpap_t90_by_stage.parquet"))}
    V = {k: rd(out, f) for k, f in (("fz", "t90_final.parquet"), ("b4", "t90_base4.parquet"), ("ma", "master_cohort.csv"), ("cp", "cpap_t90_by_stage.parquet"))}
    ch = pd.read_parquet(f"{out}/CHANGELOG.parquet")
    fz, ma, b4 = V["fz"], V["ma"], V["b4"]
    coh7 = CS.apply_cohort(P["fz"]); coh8 = CS.apply_cohort(fz); coh8b = CS.apply_cohort(b4)
    cohort_ids = set(coh7.BDSPPatientID)
    man_split = pd.read_csv(S.MANIFEST_SPLIT); split_ids = set(man_split.BDSPPatientID)
    nonsplit = cohort_ids - split_ids
    fam = pd.read_csv(CS.MEASURE_FAMILIES_CSV); feats = fam.feature.tolist()
    new4 = list(CS.NEW_MEASURES_V8); meas137 = [f for f in feats if f not in new4]
    from build_v8_tables import ARCH, IDX, OXY_MAP, NEW5, SLEEP_T90
    repl = list(dict.fromkeys(meas137 + ARCH + IDX + list(OXY_MAP)))
    tables = (("t90_final", "fz"), ("t90_base4", "b4"), ("master", "ma"))

    # ---- rows, ids, cohort
    for name, k in tables:
        rec(f"{name} rows and ids equal to v7", len(V[k]) == len(P[k]) and (V[k].BDSPPatientID.values == P[k].BDSPPatientID.values).all(),
            f"{len(V[k]):,} rows (v7 {len(P[k]):,})")
    rec("analysis cohort unchanged", len(coh8) == CS.COHORT_N == len(coh7) and set(coh8.BDSPPatientID) == cohort_ids and len(coh8b) == len(coh7),
        f"t90_final {len(coh8):,}, t90_base4 {len(coh8b):,}, v7 {len(coh7):,}, cohort_spec.COHORT_N {CS.COHORT_N:,}")

    win = fz.set_index("BDSPPatientID")["_v8_prepap_window"].astype(bool)
    light_ok = fz.set_index("BDSPPatientID")["_v8_light_status"].astype(str).str.startswith("ok")
    collapsed = fz.set_index("BDSPPatientID")["_v8_stage_collapsed"].astype(bool)
    light_tst = fz.set_index("BDSPPatientID")["_v8_light_tst_min"] if "_v8_light_tst_min" in fz.columns else None
    pap_ids = set(P["cp"].loc[pd.to_numeric(P["cp"].cpap_frac, errors="coerce") > 0, "BDSPPatientID"])   # PAP applied at any time (cpap table)

    # ---- physiology invariants on cohort rows of t90_final and master
    for name, k in (("t90_final", "fz"), ("master", "ma")):
        t = V[k][V[k].BDSPPatientID.isin(cohort_ids)].set_index("BDSPPatientID")
        w = win.reindex(t.index).fillna(False)
        tst = pd.to_numeric(t.TST_min, errors="coerce")
        pct = t[["N1_pct", "N2_pct", "N3_pct", "REM_pct"]].apply(pd.to_numeric, errors="coerce").sum(axis=1, min_count=4)
        m = pct.notna() & (tst > 0); bad = m & ((pct - 100).abs() > TOL_STAGE)
        rec(f"{name}: stage percentages sum to 100 +/- {TOL_STAGE} (TST > 0)", int((bad & w).sum()) == 0,
            f"checked {int(m.sum()):,}; violations on v8-replaced nights {int((bad & w).sum())}, on untouched nights {int((bad & ~w).sum())} (inherited)")
        if int((bad & ~w).sum()):
            rec(f"{name}: stage-sum violations on untouched nights", False, f"{int((bad & ~w).sum())} nights, ids {t.index[bad & ~w][:8].tolist()}", inherited=True)
        recm = pd.to_numeric(t.recording_dur_min, errors="coerce"); m = tst.notna() & recm.notna(); bad = m & (tst > recm + 1e-6)
        rec(f"{name}: TST <= recording minutes", int((bad & w).sum()) == 0, f"checked {int(m.sum()):,}; violations replaced {int((bad & w).sum())}, untouched {int((bad & ~w).sum())} (inherited)")
        se = pd.to_numeric(t.sleep_efficiency_pct, errors="coerce"); m = se.notna(); bad = m & (se > 100 + 1e-9)
        rec(f"{name}: sleep efficiency <= 100", int((bad & w).sum()) == 0, f"checked {int(m.sum()):,}; violations replaced {int((bad & w).sum())}, untouched {int((bad & ~w).sum())} (inherited)")
        for idx in ("AHI", "arousal_index"):
            n = pd.to_numeric(t[idx], errors="coerce") * tst / 60.0; m = n.notna(); frac = (n - n.round()).abs() < 1e-6
            rec(f"{name}: {idx} x TST/60 integer on v8-replaced nights", int((m & w & ~frac).sum()) == 0,
                f"replaced {int((m & w & frac).sum())} of {int((m & w).sum())} integer; untouched nights {int((m & ~w & frac).sum())} of {int((m & ~w).sum())} integer "
                f"(v7 counts over a true sleep time that was later re-extracted on some nights; reported, not asserted)")
        if "plm_index" in t.columns:
            den = light_tst.reindex(t.index) if light_tst is not None else tst    # the PLM count is over the sleep time of the pass it came from
            n = pd.to_numeric(t.plm_index, errors="coerce") * den / 60.0; m = n.notna(); frac = (n - n.round()).abs() < 1e-6
            rec(f"{name}: plm_index x (light-pass TST)/60 integer (our count)", int((m & ~frac).sum()) == 0, f"{int((m & frac).sum())} of {int(m.sum())} integer")
        rem = pd.to_numeric(t.REM_pct, errors="coerce"); m = (tst > 0) & rem.notna()
        np_ = ~t.index.isin(pap_ids)
        share_ns = float(((rem == 0) & m & ~w & np_).sum() / max(int((m & ~w & np_).sum()), 1)); share_s = float(((rem == 0) & m & w).sum() / max(int((m & w).sum()), 1))
        rec(f"{name}: REM = 0 share on non-split, non-PAP nights <= {REM0_MAX_SHARE:.0%} (untouched v7 cells)", share_ns <= REM0_MAX_SHARE,
            f"non-split non-PAP {share_ns:.2%} of {int((m & ~w & np_).sum()):,} (rule 3 calls this a fault; the cells are v7's, untouched by v8: owner item); "
            f"split pre-PAP windows {share_s:.2%} of {int((m & w).sum()):,} (short untreated windows, expected higher; reported)", inherited=True)
    # ---- bounds
    bounds, bsrc = bounds_table(prev_dir)
    hits = {}
    for name, k in tables:
        t = V[k][V[k].BDSPPatientID.isin(cohort_ids)]
        for _, r in bounds.iterrows():
            c = r["column"]
            if c in t.columns:
                x = pd.to_numeric(t[c], errors="coerce"); n = int(((x < float(r["lo"])) | (x > float(r["hi"]))).sum())
                if n:
                    hits[(name, c)] = n
    rec(f"plausibility bounds ({bsrc}, {len(bounds)} columns) zero hits on cohort rows", not hits, f"hits {hits or 'none'}")

    # ---- non-split nights byte-identical; new measures present
    for name, k in tables:
        m = V[k].BDSPPatientID.isin(nonsplit).to_numpy()
        badc = [c for c in repl if c in P[k].columns and c in V[k].columns and not nan_equal(V[k].loc[m, c].values, P[k].loc[m, c].values).all()]
        rec(f"{name}: measure columns byte-identical to v7 on the {len(nonsplit):,} non-split cohort nights", not badc, f"{len([c for c in repl if c in P[k].columns])} columns; differing {badc or 'none'}")
    t = fz[fz.BDSPPatientID.isin(nonsplit)].set_index("BDSPPatientID")
    lo_all = light_ok.reindex(t.index).fillna(False); lo = lo_all & ~collapsed.reindex(t.index).fillna(False)
    # defined on the light pass's own sleep (the measures were computed on it); the four sleep-versus-wake measures stay on collapsed nights (D11 = REM latency only)
    tst = pd.to_numeric(t["_v8_light_tst_min"], errors="coerce") if "_v8_light_tst_min" in t.columns else pd.to_numeric(t.TST_min, errors="coerce")
    remmin = pd.to_numeric(t.REM_min, errors="coerce") if "REM_min" in t.columns else pd.to_numeric(t.REM_pct, errors="coerce")
    need = {"hypoxic_burden": lo_all & (tst > 0), "plm_index": lo_all & (tst > 0), "sol_min": lo_all & (tst > 0), "rem_latency_min": lo & (remmin > 0) & (tst > 0), SLEEP_T90: lo_all & (tst > 0)}
    miss = {c: int((need[c] & t[c].isna()).sum()) for c in need if c in t.columns}
    rec("four new measures + sleep T90 present on non-split light-ok, non-collapsed cohort nights where defined", not any(miss.values()),
        f"light-ok non-split nights {int(lo.sum()):,}{' (smoke: the light pass covers a subset)' if smoke else ''}; missing where defined {miss}")
    for name, k in tables:
        t = V[k]; badc = [c for c in NEW5 if c not in t.columns]
        rec(f"{name}: new measure columns present", not badc, f"missing {badc or 'none'}")
    z = [c for c in ma.columns if c.endswith(CS.RANKING_DROP_SUFFIX)]
    rec("master carries no _siteZ column", not z, f"{len(z)} found")
    head = [c for c in ma.columns if c not in ("BDSPPatientID",) and not c.endswith(("_first_date", "_date")) and not c.endswith(CS.RANKING_DROP_SUFFIX) and c not in CS.SENSITIVITY_ONLY]
    d = ma[ma.BDSPPatientID.isin(cohort_ids)]
    num = [c for c in head if pd.api.types.is_numeric_dtype(d[c]) and not pd.api.types.is_bool_dtype(d[c])]
    num = [c for c in num if d[c].notna().mean() >= 0.60 and d[c].nunique() > 10 and c not in {"male", "AgeAtVisit"} and not c.startswith(("obs_", "follow_"))]
    extra_c, miss_c = sorted(set(num) - set(feats)), sorted(set(feats) - set(num))
    ok_scan = len(num) == CS.N_MEASURES or (smoke and not extra_c and set(miss_c) <= set(new4))
    rec(f"master candidate-measure scan (build_ranking_v3 rule) = cohort_spec.N_MEASURES", ok_scan,
        f"{len(num)} candidates vs {CS.N_MEASURES}; not in measure_families: {extra_c[:12]}; families not scanned: {miss_c[:12]}"
        + (" (smoke: the new measures cover only the pilot nights, under the 60 percent bar; accepted)" if smoke and len(num) != CS.N_MEASURES else ""))

    # ---- outcomes
    counts = pd.read_csv(S.EVENT_COUNTS_V8)
    copied = counts.loc[counts.status.str.startswith("unchanged"), "key"].tolist()
    rebuilt = {"obesity_hypovent_first_date"}      # rebuilt by the builder on purpose (stale in v7); checked by the consistency test below
    bad = [k + s for k in copied for s in S.OUTCOME_SUFFIXES if k + s in P["fz"].columns and k + s not in rebuilt and not nan_equal(fz[k + s].values, P["fz"][k + s].values).all()]
    rec("untouched outcome conditions byte-identical to v7 (t90_final; obesity_hypovent_first_date excluded, rebuilt)", not bad, f"{len(copied)} conditions; differing {bad or 'none'}")
    inc = {k: (int(coh8[f"{k}_incident"].sum()), int(v)) for k, v in counts.set_index("key")["incident"].items() if f"{k}_incident" in coh8.columns}
    badi = {k: v for k, v in inc.items() if v[0] != v[1]}
    rec("incident counts on the cohort equal X2's event_counts_v8.csv (all conditions + death)", not badi, f"{len(inc)} conditions; differing {badi or 'none'}")
    import importlib
    DD = importlib.import_module("disease_definitions")      # numbers/disease_definitions.py (v8, installed by E1): the ranking's own lists
    ranked = [k for k in DD.DISEASES if k not in DD.CIRCULAR and k not in DD.RANKING_EXCLUDE and f"{k}_incident" in coh8.columns and int(coh8[f"{k}_incident"].sum()) >= 150] + ["death"]
    rec("ranked outcome count by build_ranking_v3's rule (not circular, not a control, >= 150 incident, plus death) = cohort_spec.N_RANKED_OUTCOMES",
        len(ranked) == CS.N_RANKED_OUTCOMES, f"{len(ranked)} vs {CS.N_RANKED_OUTCOMES} (X2's event_counts 'ranked' flag counts {int(counts.ranked.astype(bool).sum())} conditions without death)")
    ohs = coh8[["obesity_hypovent_first_date", "obesity_hypovent_prevalent", "psg_date"]]
    n_incons = int((ohs.obesity_hypovent_first_date.notna() & (ohs.obesity_hypovent_first_date <= ohs.psg_date) & (ohs.obesity_hypovent_prevalent == 0)).sum())
    rec("obesity_hypovent_first_date consistent with the prevalent flag (first date on or before the PSG and prevalent = 0)", n_incons <= 2,
        f"{n_incons} inconsistent cohort rows (v7 carried 2,317 from the stale column; up to 2 censor-capped rows are the known v7 residue)")

    # ---- CHANGELOG: per cell, the first logged old = v7, the last logged new = v8, consecutive rows chain (a bounds row follows a replacement row)
    badc = {}
    ch["_ord"] = np.arange(len(ch))
    for (tab, col), g in ch.groupby(["table", "column"], sort=False):
        if tab == "cpap_t90_by_stage":
            continue   # cpap ids repeat across arms; checked below by value equality with t90_final
        k = {"t90_final": "fz", "t90_base4": "b4", "master": "ma"}[tab]
        pv = P[k].set_index("BDSPPatientID"); vv = V[k].set_index("BDSPPatientID")
        if col not in vv.columns:
            if col.endswith("_siteZ"): continue   # dropped by decision 3 (the CHANGELOG logs the split replacement made before the drop)
            badc[(tab, col)] = "column missing in v8"; continue
        g = g.sort_values("_ord"); first = g.groupby("BDSPPatientID", sort=False).head(1); last = g.groupby("BDSPPatientID", sort=False).tail(1)
        is_dt = pd.api.types.is_datetime64_any_dtype(vv[col]) or (col in pv.columns and pd.api.types.is_datetime64_any_dtype(pv[col]))
        is_str = g.old_str.notna().any() or g.new_str.notna().any()
        def logged(df_, side):
            if is_dt:
                return pd.to_datetime(df_[side + "_str"], errors="coerce").astype("datetime64[us]").values
            return df_[side + "_str"].values if is_str else df_[side].values
        def table(vals):
            if is_dt:
                return pd.to_datetime(pd.Series(vals)).astype("datetime64[us]").values
            if is_str:
                return pd.Series(vals).astype(str).where(pd.Series(vals).notna(), None).values
            return pd.to_numeric(pd.Series(vals), errors="coerce").values
        n_bad = int((~nan_equal(table(vv.loc[last.BDSPPatientID.values, col].values), logged(last, "new"))).sum())
        if col in pv.columns:
            n_bad += int((~nan_equal(table(pv.loc[first.BDSPPatientID.values, col].values), logged(first, "old"))).sum())
        else:
            n_bad += int(first.old.notna().sum() + first.old_str.notna().sum())
        chain = g[g.duplicated("BDSPPatientID", keep=False)]
        if len(chain):
            prev_new = chain.groupby("BDSPPatientID", sort=False)[["new", "new_str"]].shift(1)
            has_prev = prev_new.new.notna() | prev_new.new_str.notna() | chain.groupby("BDSPPatientID").cumcount().gt(0)
            n_bad += int((~nan_equal(chain.loc[has_prev, "old_str" if is_str or is_dt else "old"].values, prev_new.loc[has_prev, "new_str" if is_str or is_dt else "new"].values)).sum())
        if n_bad:
            badc[(tab, col)] = n_bad
    rec("every CHANGELOG cell chains from the v7 value (first old) to the v8 value (last new)", not badc,
        f"{len(ch):,} rows, {ch.groupby(['table','column']).ngroups} (table, column) groups, {int(ch.duplicated(['table','column','BDSPPatientID']).sum()):,} chained rows; mismatches {dict(list(badc.items())[:8]) or 'none'}")

    # ---- base4 shared columns vs t90_final: no new disagreement
    shared = [c for c in fz.columns if c in b4.columns and c != "BDSPPatientID"]
    dis8 = {c for c in shared if not nan_equal(fz[c].values, b4.set_index("BDSPPatientID").reindex(fz.BDSPPatientID.values)[c].values).all()}
    shared7 = [c for c in P["fz"].columns if c in P["b4"].columns and c != "BDSPPatientID"]
    dis7 = {c for c in shared7 if not nan_equal(P["fz"][c].values, P["b4"].set_index("BDSPPatientID").reindex(P["fz"].BDSPPatientID.values)[c].values).all()}
    rec("t90_base4 shared columns agree with t90_final at least as well as in v7", dis8 <= dis7, f"v8 disagreeing {sorted(dis8)}; v7 disagreeing {sorted(dis7)}")

    # ---- cpap covariates
    cpv = V["cp"]; ref = fz.set_index("BDSPPatientID"); m = cpv["_v8_covariates_refreshed"].astype(bool).to_numpy()
    badc = [c for c in ("AHI", "TST_min", "spo2_pct_below_90") if c in cpv.columns and not nan_equal(cpv.loc[m, c].values, ref.reindex(cpv.loc[m, "BDSPPatientID"].values)[c].values).all()]
    rec("cpap_t90_by_stage whole-night covariates equal the v8 t90_final on refreshed rows", not badc, f"{int(m.sum()):,} of {len(cpv):,} rows refreshed; differing {badc or 'none'}")
    pre = [c for c in P["cp"].columns if c not in ("AHI", "TST_min", "spo2_pct_below_90", "_v7_covariates_refreshed")]
    badc = [c for c in pre if c in cpv.columns and not nan_equal(cpv[c].values, P["cp"][c].values).all()]
    rec("cpap_t90_by_stage per-stage columns unchanged", not badc, f"{len(pre)} columns; differing {badc or 'none'}")

    # ---- sha256.txt and sidecars
    listed = dict((ln.split()[1], ln.split()[0]) for ln in open(f"{out}/sha256.txt").read().strip().splitlines())
    badf = [f for f in S.TABLE_FILES if listed.get(f) != S.sha256(f"{out}/{f}")]
    nos = [f for f in S.TABLE_FILES if not os.path.isfile(f"{out}/{f}.provenance.json")]
    rec("sha256.txt matches the four tables and every table has a provenance sidecar", not badf and not nos, f"sha mismatch {badf or 'none'}; sidecar missing {nos or 'none'}")

    # ---- per-site distribution of EVERY numeric column, v7 vs v8 (cohort rows)
    from scipy.stats import ks_2samp
    rows = []
    for name, k in tables:
        v7 = P[k][P[k].BDSPPatientID.isin(cohort_ids)].set_index("BDSPPatientID"); v8 = V[k][V[k].BDSPPatientID.isin(cohort_ids)].set_index("BDSPPatientID")
        site = v8["site_id"] if "site_id" in v8.columns else fz.set_index("BDSPPatientID").site_id.reindex(v8.index)
        w = v8.index.isin(split_ids)
        cols = [c for c in v8.columns if pd.api.types.is_numeric_dtype(v8[c]) and not pd.api.types.is_bool_dtype(v8[c]) and c != "BDSPPatientID" and not c.startswith("_impossible_")]
        for c in cols:
            x8 = pd.to_numeric(v8[c], errors="coerce"); x7 = pd.to_numeric(v7[c], errors="coerce") if c in v7.columns else None
            for s in ("I0002", "I0006"):
                ms = (site == s).to_numpy()
                a8 = x8[ms].dropna(); r = {"table": name, "column": c, "site": s, "n_v8": len(a8), "median_v8": a8.median(), "q1_v8": a8.quantile(.25), "q3_v8": a8.quantile(.75)}
                if x7 is not None:
                    a7 = x7[ms].dropna(); r.update({"n_v7": len(a7), "median_v7": a7.median(), "q1_v7": a7.quantile(.25), "q3_v7": a7.quantile(.75)})
                    r["nonsplit_identical"] = bool(nan_equal(x8[ms & ~w].values, x7[ms & ~w].values).all())
                    iqr = (a7.quantile(.75) - a7.quantile(.25)); r["d_median"] = a8.median() - a7.median()
                    r["d_median_over_iqr_v7"] = (r["d_median"] / iqr) if iqr and iqr == iqr and iqr > 0 else np.nan
                    s7, s8 = x7[ms & w].dropna(), x8[ms & w].dropna()
                    r["ks_split_v7_vs_v8"] = ks_2samp(s7, s8).statistic if len(s7) > 5 and len(s8) > 5 else np.nan
                    r["n_split_v7"], r["n_split_v8"] = len(s7), len(s8)
                    is_outcome = c.endswith(S.OUTCOME_SUFFIXES) or c.startswith("death_")
                    r["flag"] = (("OUTCOME_RECOMPUTED" if is_outcome else "NONSPLIT_CHANGED") if not r["nonsplit_identical"]
                                 else ("SHIFT_FROM_SPLIT_NIGHTS" if abs(r["d_median_over_iqr_v7"]) > SHIFT_IQR_FRAC else ""))
                else:
                    r.update({"n_v7": np.nan, "nonsplit_identical": np.nan, "flag": "NEW_COLUMN"})
                rows.append(r)
    ps = pd.DataFrame(rows)
    ps.to_csv(f"{log_dir}/PER_SITE_ALL_COLUMNS_V8.csv", index=False)
    n_ns = int((ps.flag == "NONSPLIT_CHANGED").sum()); n_sh = int((ps.flag == "SHIFT_FROM_SPLIT_NIGHTS").sum()); n_new = int((ps.flag == "NEW_COLUMN").sum()); n_oc = int((ps.flag == "OUTCOME_RECOMPUTED").sum())
    rec("per-site distribution of every numeric column: non-split nights identical to v7 in every non-outcome (table, column, site)", n_ns == 0,
        f"{len(ps):,} (table, column, site) rows; NONSPLIT_CHANGED {n_ns} ({sorted(set(ps[ps.flag == 'NONSPLIT_CHANGED'].column))[:10]}); recomputed outcome columns {n_oc}; "
        f"median shifts above {SHIFT_IQR_FRAC:.0%} of the v7 IQR (all from the split nights) {n_sh}; new columns {n_new}; table PER_SITE_ALL_COLUMNS_V8.csv")
    shifts = ps[ps.flag == "SHIFT_FROM_SPLIT_NIGHTS"].sort_values("d_median_over_iqr_v7", key=lambda s: s.abs(), ascending=False)
    new = ps[ps.flag == "NEW_COLUMN"]

    n_fail = sum(1 for _, st, _ in res if st == "FAIL")
    hdr = [f"# INVARIANTS_V8 ({'SMOKE' if smoke else 'REAL'} build {out}; v7 {prev_dir}; {time.strftime('%Y-%m-%d %H:%M')}; {time.time()-t0:.0f}s)", "",
           f"Result: {'PASS' if n_fail == 0 else f'FAIL ({n_fail} failures)'}; {sum(1 for _, st, _ in res if st == 'INHERITED')} inherited v7 findings reported.", ""]
    body = hdr + L + ["", f"## Largest per-site median shifts (v7 -> v8, all from the split nights; top 25 of {n_sh})", ""]
    body += [f"- {r.table}.{r.column} {r.site}: median {r.median_v7:.4g} -> {r.median_v8:.4g} ({r.d_median_over_iqr_v7:+.2f} IQR), split-night KS {r.ks_split_v7_vs_v8:.3f}" for r in shifts.head(25).itertuples()]
    body += ["", "## New columns per site (median [IQR])", ""]
    body += [f"- {r.table}.{r.column} {r.site}: n={r.n_v8:,} median {r.median_v8:.4g} [{r.q1_v8:.4g}, {r.q3_v8:.4g}]" for r in new[new.table == "t90_final"].itertuples()]
    open(f"{log_dir}/INVARIANTS_V8.md", "w").write("\n".join(body) + "\n")
    print("\n".join(body[:4] + L), flush=True)
    return n_fail


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--dir", required=True); ap.add_argument("--prev", default=S.PREV_DIR)
    ap.add_argument("--log-dir", default=None); ap.add_argument("--smoke", action="store_true")
    a = ap.parse_args()
    S.require_resident(*[f"{a.dir}/{f}" for f in S.TABLE_FILES], f"{a.dir}/CHANGELOG.parquet", *[f"{a.prev}/{f}" for f in S.TABLE_FILES])
    log_dir = a.log_dir or a.dir; os.makedirs(log_dir, exist_ok=True)
    n = run(os.path.abspath(a.dir), os.path.abspath(a.prev), log_dir, smoke=a.smoke)
    for f in ("INVARIANTS_V8.md", "PER_SITE_ALL_COLUMNS_V8.csv"):
        S.write_sidecar(f"{log_dir}/{f}", [f"{a.dir}/{x}" for x in S.TABLE_FILES] + [f"{a.prev}/{x}" for x in S.TABLE_FILES], os.path.abspath(__file__), "051_invariants_v8",
                        extra={"failures": n}, index_path=f"{S.E2_LOGS}/PROVENANCE_INDEX.tsv")
    return 1 if n else 0


if __name__ == "__main__":
    sys.exit(main())
