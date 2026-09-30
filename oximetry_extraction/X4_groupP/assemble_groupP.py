#!/usr/bin/env python3
"""Assemble the group P re-extraction (v8.2, 2026-09-15): fleet shards -> collapse-masked shards in the consumer folders,
provenance sidecars, invariants and the old-versus-new per-site comparison against the August extraction.
usage: assemble_groupP.py cpap | t90 | oxy      (each writes X4_groupP/logs/ASSEMBLE_<fmt>.md and .json)
  cpap : raw/{cpapA,cpapP,cpapB}/out_<s>.csv  -> cpap_stage/out{A,P,B}_<s>.csv (masked) + MASK_<arm>_<s>.json
  t90  : t90_by_stage/raw/out_<s>.csv         -> t90_by_stage/out_<s>.csv (masked) + MASK_<s>.json
  oxy  : oxygen_profile/raw/out_<s>.jsonl     -> Sleep_Variability_2026-08/oxygen_profile/_raw_fleet_v8_2/out_<s>.jsonl (masked)
The unmasked fleet output stays under raw/ (never read by a consumer)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../..")))  # repository root, where paths.py lives
import paths
import glob, json, os, sys, time
import numpy as np, pandas as pd
CCP = paths.T90_ROOT; T90 = f"{paths.T90_ROOT}"; X4 = f"{paths.V8_ROOT}/X4_groupP"
B = f"{paths.T90_ROOT}/pilots"; OXY = f"{paths.SV_ROOT}/oxygen_profile"
sys.path.insert(0, paths.X4_DIR); sys.path.insert(0, paths.ANALYSIS_DIR)
import mask_collapsed as MC
from cohort_spec import sidecar
TABLE = f"{paths.TABLES_DIR}/t90_final.parquet"
STATES = ("wake", "sleep", "nrem", "n1", "n2", "n3", "rem")


def cohort_flags():
    t = pd.read_parquet(TABLE, columns=["BDSPPatientID", "fu_valid", "spo2_pct_below_90", "oximetry_bad", "site_id", "_v8_stage_collapsed", "_v8_prepap_window"])
    t = t[(t.fu_valid == 1) & t.spo2_pct_below_90.notna() & (t.oximetry_bad == 0)]
    return set(t.BDSPPatientID.astype(int)), set(t.loc[t._v8_stage_collapsed.astype(bool), "BDSPPatientID"].astype(int)), t.set_index("BDSPPatientID").site_id.to_dict()


def per_site_compare(new, old, cols, key, site):
    rows = []
    for c in cols:
        if c not in new.columns or c not in old.columns:
            continue
        for s in sorted(set(new[site].dropna())):
            a = pd.to_numeric(new.loc[new[site] == s, c], errors="coerce").dropna(); b = pd.to_numeric(old.loc[old[site] == s, c], errors="coerce").dropna()
            if len(a) < 20 or len(b) < 20:
                continue
            rows.append({"column": c, "site": s, "n_new": len(a), "n_old": len(b), "median_new": float(a.median()), "median_old": float(b.median()),
                         "q1_new": float(a.quantile(.25)), "q1_old": float(b.quantile(.25)), "q3_new": float(a.quantile(.75)), "q3_old": float(b.quantile(.75))})
    return pd.DataFrame(rows)


def run_cpap():
    ids, coll, site_of = cohort_flags()
    out_lines, masks, all_new, all_old = [], {}, [], []
    for arm in ("A", "P", "B"):
        raw = sorted(glob.glob(f"{X4}/cpap_stage/raw/cpap{arm}/out_*.csv"))
        assert raw, f"no raw shards for cpap{arm}"
        for f in raw:
            s = os.path.basename(f).split("_")[1].split(".")[0]
            df = pd.read_csv(f, low_memory=False); df, nc = MC.mask_cpap(df)
            dst = f"{X4}/cpap_stage/out{arm}_{s}.csv"; df.to_csv(dst, index=False)
            rep = {"arm": arm, "shard": s, "records": int(len(df)), "ok": int((df.status == "ok").sum()), "collapsed": nc, "raw": f, "masked": dst}
            json.dump(rep, open(f"{X4}/cpap_stage/MASK_{arm}_{s}.json", "w"), indent=1); masks[f"{arm}_{s}"] = rep; all_new.append(df.assign(arm_=arm))
        old = pd.concat([pd.read_csv(f, low_memory=False) for f in sorted(glob.glob(f"{B}/cpap_stage/out{arm}_*.csv"))], ignore_index=True).assign(arm_=arm); all_old.append(old)
    new = pd.concat(all_new, ignore_index=True); old = pd.concat(all_old, ignore_index=True)
    wl = sum(len(pd.read_csv(f"{X4}/worklists/cpap_worklist_{a}.csv", low_memory=False)) for a in ("A", "P", "B"))
    inv = {"records": int(len(new)), "worklist_rows": int(wl), "status": new.status.value_counts().to_dict(), "collapsed": int(new._v8_stage_collapsed.sum()),
           "collapsed_in_cohort_flag": int(new.loc[new._v8_stage_collapsed == 1, "BDSPPatientID"].astype(int).isin(coll).sum()),
           "cohort_flag_collapsed_seen": int(new.BDSPPatientID.astype(int).isin(coll).sum())}
    ok = new[new.status == "ok"]
    inv["t90_in_range"] = bool(all(((ok[c].dropna() >= 0) & (ok[c].dropna() <= 100)).all() for c in ok.columns if c.endswith("_t90")))
    nc_ = ok[ok._v8_stage_collapsed == 0]
    sm = nc_[[f"all_{s}_min" for s in ("wake", "n1", "n2", "n3", "rem")]].sum(axis=1)
    inv["state_minutes_equal_staged_minutes_max_abs_diff_noncollapsed"] = float((sm - nc_.staged_min).abs().max())
    rawok = pd.concat([pd.read_csv(f, low_memory=False) for f in sorted(glob.glob(f"{X4}/cpap_stage/raw/cpap*/out_*.csv"))], ignore_index=True); rawok = rawok[rawok.status == "ok"]
    smr = rawok[[f"all_{s}_min" for s in ("wake", "n1", "n2", "n3", "rem")]].sum(axis=1)
    inv["state_minutes_equal_staged_minutes_max_abs_diff_raw_all_ok"] = float((smr - rawok.staged_min).abs().max())
    inv["collapsed_equals_flag_count"] = bool(int(new._v8_stage_collapsed.sum()) == sum(m["collapsed"] for m in masks.values()))
    inv["collapsed_stage_values_all_missing"] = bool(ok.loc[ok._v8_stage_collapsed == 1, [f"all_{s}_t90" for s in ("nrem", "n1", "n2", "n3", "rem")]].isna().all().all())
    a = ok[np.isfinite(ok.cpap_start_min)]
    inv["pre_plus_post_equals_all_max_abs_diff_min"] = float(((a.pre_any_min + a.post_any_min) - a.all_any_min).abs().max()) if len(a) else None
    inv["records_equal_worklist"] = bool(len(new) == wl)
    cols = [f"all_{s}_t90" for s in STATES] + ["all_any_t90", "cpap_start_min", "tst_min"]
    cmp = per_site_compare(ok, old[old.status == "ok"], cols, "BDSPPatientID", "SiteID"); cmp.to_csv(f"{X4}/logs/ASSEMBLE_cpap_per_site.csv", index=False)
    return inv, masks, cmp


def run_t90():
    ids, coll, site_of = cohort_flags()
    raw = sorted(glob.glob(f"{X4}/t90_by_stage/raw/out_*.csv")); assert raw, "no raw t90 shards"
    masks, frames = {}, []
    for f in raw:
        s = os.path.basename(f).split("_")[1].split(".")[0]
        df = pd.read_csv(f, low_memory=False); df, nc = MC.mask_t90(df)
        dst = f"{X4}/t90_by_stage/out_{s}.csv"; df.to_csv(dst, index=False)
        rep = {"shard": s, "records": int(len(df)), "ok": int((df.status == "ok").sum()), "collapsed": nc, "raw": f, "masked": dst}
        json.dump(rep, open(f"{X4}/t90_by_stage/MASK_{s}.json", "w"), indent=1); masks[s] = rep; frames.append(df)
    new = pd.concat(frames, ignore_index=True)
    old = pd.concat([pd.read_csv(f, low_memory=False) for f in sorted(glob.glob(f"{B}/t90_by_stage/out_*.csv"))], ignore_index=True)
    old = old[old.BDSPPatientID.isin(ids)]
    inv = {"records": int(len(new)), "unique_patients": int(new.BDSPPatientID.nunique()), "cohort": len(ids), "status": new.status.value_counts().to_dict(),
           "collapsed": int(new._v8_stage_collapsed.sum()), "collapsed_in_cohort_flag": int(new.loc[new._v8_stage_collapsed == 1, "BDSPPatientID"].astype(int).isin(coll).sum()),
           "cohort_flag_collapsed": len(coll), "missing_patients": int(len(ids - set(new.BDSPPatientID.astype(int))))}
    ok = new[new.status == "ok"]
    inv["t90_in_range"] = bool(all(((ok[c].dropna() >= 0) & (ok[c].dropna() <= 100)).all() for c in ok.columns if c.startswith("t90_")))
    nc_ = ok[ok._v8_stage_collapsed == 0]
    sm = nc_[[f"min_{s}" for s in ("wake", "n1", "n2", "n3", "rem")]].sum(axis=1)
    inv["state_minutes_vs_pap_off_min_max_abs_diff_noncollapsed"] = float((sm - nc_.pap_off_min).abs().max()) if "pap_off_min" in ok.columns else None
    inv["state_minutes_vs_pap_off_min_median_abs_diff_noncollapsed"] = float((sm - nc_.pap_off_min).abs().median()) if "pap_off_min" in ok.columns else None
    inv["sleep_minutes_vs_tst_papoff_max_abs_diff_noncollapsed"] = float((nc_[[f"min_{s}" for s in ("n1", "n2", "n3", "rem")]].sum(axis=1) - nc_.tst_papoff_min).abs().max()) if "tst_papoff_min" in ok.columns else None
    inv["collapsed_equals_flag_count"] = bool(int(new._v8_stage_collapsed.sum()) == sum(m["collapsed"] for m in masks.values()))
    inv["collapsed_stage_values_all_missing"] = bool(ok.loc[ok._v8_stage_collapsed == 1, [f"t90_{s}" for s in ("nrem", "n1", "n2", "n3", "rem")]].isna().all().all())
    cols = [f"t90_{s}" for s in ("all",) + STATES] + [f"min_{s}" for s in STATES] + ["pap_off_min", "tst_papoff_min"]
    cmp = per_site_compare(ok, old[old.status == "ok"], cols, "BDSPPatientID", "SiteID"); cmp.to_csv(f"{X4}/logs/ASSEMBLE_t90_per_site.csv", index=False)
    return inv, masks, cmp


def run_oxy():
    ids, coll, site_of = cohort_flags()
    raw = sorted(glob.glob(f"{X4}/oxygen_profile/raw/out_*.jsonl")); assert raw, "no raw oxy shards"
    dstdir = f"{OXY}/_raw_fleet_v8_2"; os.makedirs(dstdir, exist_ok=True)
    masks, rows = {}, []
    for f in raw:
        s = os.path.basename(f).split("_")[1].split(".")[0]; n = nc = 0
        dst = f"{dstdir}/out_{s}.jsonl"
        with open(f) as fi, open(dst, "w") as fo:
            for line in fi:
                if not line.strip():
                    continue
                r, fl = MC.mask_oxy_record(json.loads(line)); n += 1; nc += int(bool(fl)); rows.append(r); fo.write(json.dumps(r) + "\n")
        rep = {"shard": s, "records": n, "collapsed": nc, "raw": f, "masked": dst}
        json.dump(rep, open(f"{X4}/oxygen_profile/MASK_{s}.json", "w"), indent=1); masks[s] = rep
    new = pd.DataFrame(rows); new["BDSPPatientID"] = new.BDSPPatientID.astype(int)
    orows = []
    for f in sorted(glob.glob(f"{OXY}/_raw_fleet/*.jsonl")):
        with open(f) as fh:
            for line in fh:
                if line.strip():
                    orows.append(json.loads(line))
    old = pd.DataFrame(orows); old["BDSPPatientID"] = old.BDSPPatientID.astype(int)
    inv = {"records": int(len(new)), "unique_patients": int(new.BDSPPatientID.nunique()), "cohort": len(ids), "status": new.status.value_counts().to_dict(),
           "collapsed": int(new._v8_stage_collapsed.sum()), "collapsed_in_cohort_flag": int(new.loc[new._v8_stage_collapsed == 1, "BDSPPatientID"].isin(coll).sum()),
           "cohort_flag_collapsed": len(coll), "missing_patients": int(len(ids - set(new.BDSPPatientID)))}
    ok = new[new.status == "ok"]
    inv["t90_in_range"] = bool(all(((pd.to_numeric(ok[c], errors="coerce").dropna() >= 0) & (pd.to_numeric(ok[c], errors="coerce").dropna() <= 100)).all() for c in ok.columns if c.endswith("_t90")))
    nc_ = ok[ok._v8_stage_collapsed == 0]
    sm = nc_[[f"{s}_min" for s in ("wake", "n1", "n2", "n3", "rem")]].apply(pd.to_numeric, errors="coerce").sum(axis=1)
    inv["state_minutes_vs_staged_min_max_abs_diff_noncollapsed"] = float((sm - pd.to_numeric(nc_.staged_min, errors="coerce")).abs().max())
    inv["collapsed_equals_flag_count"] = bool(int(new._v8_stage_collapsed.sum()) == sum(m["collapsed"] for m in masks.values()))
    inv["collapsed_stage_values_all_missing"] = bool(ok.loc[ok._v8_stage_collapsed == 1, [f"{s}_mean" for s in ("n1", "n2", "n3", "rem")]].isna().all().all())
    new["site"] = new.BDSPPatientID.map(site_of); old["site"] = old.BDSPPatientID.map(site_of)
    cols = [f"{s}_t90" for s in ("rec", "wake", "sleep", "n1", "n2", "n3", "rem")] + [f"{s}_mean" for s in ("rec", "wake", "sleep", "n1", "n2", "n3", "rem")] + [f"dec_mean_{i:02d}" for i in range(1, 11)]
    cmp = per_site_compare(ok, old[old.status == "ok"], cols, "BDSPPatientID", "site"); cmp.to_csv(f"{X4}/logs/ASSEMBLE_oxy_per_site.csv", index=False)
    return inv, masks, cmp


def main():
    fmt = sys.argv[1]
    inv, masks, cmp = {"cpap": run_cpap, "t90": run_t90, "oxy": run_oxy}[fmt]()
    rep = {"format": fmt, "time": time.strftime("%Y-%m-%d %H:%M"), "invariants": inv, "shards": masks}
    json.dump(rep, open(f"{X4}/logs/ASSEMBLE_{fmt}.json", "w"), indent=1)
    L = [f"# ASSEMBLE {fmt} ({rep['time']})", "", "## Invariants"] + [f"- {k}: {v}" for k, v in inv.items()] + ["", "## Per-site medians, new vs August (q1 to q3)"]
    for r in cmp.itertuples():
        L.append(f"- {r.column} {r.site}: new {r.median_new:.3f} ({r.q1_new:.3f} to {r.q3_new:.3f}, n {r.n_new}) | August {r.median_old:.3f} ({r.q1_old:.3f} to {r.q3_old:.3f}, n {r.n_old})")
    open(f"{X4}/logs/ASSEMBLE_{fmt}.md", "w").write("\n".join(L) + "\n")
    for m in masks.values():
        sidecar(m["masked"], __file__, extra_inputs=[m["raw"]], note=f"v8.2 group P re-extraction, collapse-masked ({m['collapsed']} of {m['records']})")
    print("\n".join(L[:14]))


if __name__ == "__main__":
    main()
