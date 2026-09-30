"""
Step 020 checker. On the 300-night light pilot (sample_300.csv, split nights cut at t_end):
  * positive controls on NON-split nights: sol_min, rem_latency_min identical to reextract_FINAL; whole-recording native
    T90 identical to the v7 native value (spo2_pct_below_90_native); April-style T90 (spo2_pct_below_90) identical too;
    TST_min identical;
  * limb label vocabulary and duration distribution per site; LM against PLM-series counts per site (B3);
  * event-locked ensemble SpO2 curve per site (from the per-night curve sums) written to hb_window_calibration.json with the
    proposed W_pre / W_post per site (pre-nadir maximum to post-nadir return on the ensemble mean, see freeze_hb_definition);
  * per-site medians of every new measure;
  * on split nights: pre-PAP native T90 against cpap_t90_by_stage.pre_any_t90 and pre-PAP sleep minutes against pre_sleep_min
    (the 013 gate, run early on the 56 split nights in the sample).
Writes logs/PILOT300_LIGHT.md and the calibration json, each with a sidecar.
"""
from __future__ import annotations
import argparse, json, math, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import x1_paths as P
from provenance import write_sidecar

NEW = ["hypoxic_burden", "hb_area_pctmin", "hb_n_events_used", "plm_index", "lm_index", "plm_index_excl_resp", "lm_count", "plm_count",
       "sol_min", "rem_latency_min", "spo2_pct_below_90_sleep", "spo2_sleep_valid_min", "spo2_pct_below_90_native", "spo2_pct_below_90"]


def q(s):
    s = pd.to_numeric(s, errors="coerce")
    return f"n={int(s.notna().sum())} median={s.median():.3g} IQR=[{s.quantile(.25):.3g}, {s.quantile(.75):.3g}] NaN={int(s.isna().sum())}"


def ensemble(df, site):
    sub = df[(df.site_id == site) & df.hb_curve_count.notna()]
    if sub.empty:
        return None
    offs = np.array(json.loads(sub.iloc[0].hb_curve_offsets_sec) if isinstance(sub.iloc[0].hb_curve_offsets_sec, str) else sub.iloc[0].hb_curve_offsets_sec, float)
    S = np.zeros(offs.size); D = np.zeros(offs.size); C = np.zeros(offs.size); n_ev = 0
    for _, r in sub.iterrows():
        c = np.array(r.hb_curve_count, float); S += np.array(r.hb_curve_spo2_sum, float); D += np.array(r.hb_curve_deficit_sum, float); C += c
        n_ev += int(r.get("hb_n_events_used", 0) or 0)
    mean = np.where(C > 0, S / np.maximum(C, 1), np.nan); deficit = np.where(C > 0, D / np.maximum(C, 1), np.nan)
    # window proposal: nadir of the ensemble mean after the event end; pre = last local maximum before the nadir (searching
    # back from the nadir to -180 s); post = first point after the nadir where the mean recovers to within 0.1 point of the
    # pre-event maximum, else the post-nadir maximum.
    k0 = int(np.where(offs == 0)[0][0])
    sm = pd.Series(mean).rolling(5, center=True, min_periods=1).mean().to_numpy()
    k_nadir = k0 + int(np.nanargmin(sm[k0: k0 + 61]))          # nadir within 60 s after the end
    pre_seg = sm[: k_nadir + 1]
    k_pre = int(np.nanargmax(pre_seg))
    post_seg = sm[k_nadir:]
    k_postmax = k_nadir + int(np.nanargmax(post_seg))
    recov = np.where(post_seg >= sm[k_pre] - 0.1)[0]
    k_post = k_nadir + int(recov[0]) if recov.size else k_postmax
    return {"site": site, "n_nights": int(len(sub)), "n_events": n_ev, "offsets_sec": offs.tolist(), "mean_spo2": [None if np.isnan(v) else round(float(v), 4) for v in mean],
            "mean_deficit": [None if np.isnan(v) else round(float(v), 4) for v in deficit], "count": C.tolist(),
            "nadir_offset_sec": float(offs[k_nadir]), "pre_max_offset_sec": float(offs[k_pre]), "post_return_offset_sec": float(offs[k_post]),
            "post_max_offset_sec": float(offs[k_postmax]), "spo2_at_pre_max": float(sm[k_pre]), "spo2_at_nadir": float(sm[k_nadir]),
            "proposed_w_pre": float(-offs[k_pre]), "proposed_w_post": float(offs[k_post])}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--jsonl", default=f"{P.X1_DATA}/light_pilot300.jsonl")
    ap.add_argument("--report", default=f"{P.X1_LOGS}/PILOT300_LIGHT.md"); ap.add_argument("--calib", default=P.HB_CALIBRATION)
    a = ap.parse_args()
    P.assert_not_evicted(a.jsonl, P.REEXTRACT_FINAL, P.CPAP_STAGE_V7, P.SAMPLE_300)
    df = pd.DataFrame([json.loads(l) for l in open(a.jsonl)])
    df = df.drop_duplicates("BDSPPatientID", keep="last")
    man = pd.read_csv(P.SAMPLE_300)
    L = [f"# PILOT300_LIGHT ({len(df)} records; manifest {len(man)})", f"status: {df.status.astype(str).str[:40].value_counts().to_dict()}"]
    ok = df[df.status == "ok"].copy()
    ok["is_split"] = ok["_v8_prepap_window"].astype(bool)
    fin = pd.read_parquet(P.REEXTRACT_FINAL, columns=["BDSPPatientID", "sol_min", "rem_latency_min", "spo2_pct_below_90_native", "spo2_pct_below_90", "TST_min", "AHI", "limb_rows"])
    ns = ok[~ok.is_split].merge(fin, on="BDSPPatientID", suffixes=("", "_v7"))
    L.append(f"\n## Positive controls on {len(ns)} non-split nights (exact equality, NaN = NaN)")
    fails = {}
    for c in ("sol_min", "rem_latency_min", "spo2_pct_below_90_native", "spo2_pct_below_90", "TST_min"):
        x = pd.to_numeric(ns[c], errors="coerce"); y = pd.to_numeric(ns[f"{c}_v7"], errors="coerce")
        eq = ((x == y) | (x.isna() & y.isna()))
        fails[c] = ns.loc[~eq, ["BDSPPatientID", "site_id", c, f"{c}_v7"]]
        L.append(f"- {c}: {int(eq.sum())} of {len(ns)} identical; max |diff| {float((x - y).abs().max()):.3g}; mismatches {len(fails[c])}")
        if len(fails[c]):
            L.append(fails[c].head(10).to_string(index=False))
    L.append(f"\n## Split nights in the sample ({int(ok.is_split.sum())}): pre-PAP window against the v7 pilot table (the 013 gate, early)")
    cp = pd.read_parquet(P.CPAP_STAGE_V7, columns=["BDSPPatientID", "pre_any_t90", "pre_sleep_min", "cpap_start_min"])
    sp = ok[ok.is_split].merge(cp, on="BDSPPatientID")
    if len(sp):
        from scipy.stats import spearmanr
        x = pd.to_numeric(sp.spo2_pct_below_90_native, errors="coerce"); y = sp.pre_any_t90
        m = x.notna() & y.notna()
        rho = spearmanr(x[m], y[m]).correlation if m.sum() > 3 else np.nan
        L.append(f"- pre-PAP native T90 vs pre_any_t90: n={int(m.sum())} Spearman={rho:.4f} median|diff|={float((x - y)[m].abs().median()):.3f} "
                 f"p95|diff|={float((x - y)[m].abs().quantile(.95)):.3f} max={float((x - y)[m].abs().max()):.3f}")
        d = (pd.to_numeric(sp.TST_min, errors="coerce") - sp.pre_sleep_min).abs()
        L.append(f"- pre-PAP TST vs pre_sleep_min: n={int(d.notna().sum())} within 1 min {float((d <= 1).mean()*100):.1f}% median|diff|={float(d.median()):.3f} max={float(d.max()):.2f}")
        L.append(f"- window end vs cpap_start_min: max |diff| {float((sp._v8_window_end_min - sp.cpap_start_min).abs().max()):.4f} min")
        worst = sp.assign(dT90=(x - y).abs(), dTST=d).sort_values("dTST", ascending=False)[["BDSPPatientID", "site_id", "TST_min", "pre_sleep_min", "spo2_pct_below_90_native", "pre_any_t90", "_v8_window_end_min", "cpap_start_min"]].head(8)
        L.append(worst.to_string(index=False))
    L.append("\n## Row invariants")
    pct = ok[["N1_pct", "N2_pct", "N3_pct", "REM_pct"]].apply(pd.to_numeric, errors="coerce").sum(axis=1, min_count=4)
    L.append(f"- stage percentages sum to 100: {int(((pct - 100).abs() < 1e-6).sum())} of {int(pct.notna().sum())} with TST>0 (NaN where TST=0: {int(pct.isna().sum())}; ids {ok.loc[pct.isna(), 'BDSPPatientID'].tolist()})")
    se = pd.to_numeric(ok.sleep_efficiency_pct, errors="coerce"); L.append(f"- sleep efficiency <= 100: {int((se <= 100 + 1e-9).sum())} of {int(se.notna().sum())}")
    ahi_n = pd.to_numeric(ok.AHI, errors="coerce") * pd.to_numeric(ok.TST_min, errors="coerce") / 60
    L.append(f"- AHI x TST/60 integer: {int((ahi_n.dropna() - ahi_n.dropna().round()).abs().lt(1e-6).sum())} of {int(ahi_n.notna().sum())}")
    plm_n = pd.to_numeric(ok.plm_index, errors="coerce") * pd.to_numeric(ok.TST_min, errors="coerce") / 60
    L.append(f"- PLMI x TST/60 integer: {int((plm_n.dropna() - plm_n.dropna().round()).abs().lt(1e-6).sum())} of {int(plm_n.notna().sum())}")
    hb = pd.to_numeric(ok.hypoxic_burden, errors="coerce"); nev = pd.to_numeric(ok.hb_n_events_used, errors="coerce")
    L.append(f"- HB >= 0: {int((hb >= 0).sum())} of {int(hb.notna().sum())}; HB == 0 iff no events used: {int(((hb == 0) == (nev == 0))[hb.notna()].sum())} of {int(hb.notna().sum())}")
    L.append(f"- sol_min == sol_min_nm and rem_latency == rem_latency_nm: {int(((ok.sol_min == ok.sol_min_nm) | (ok.sol_min.isna() & ok.sol_min_nm.isna())).sum())} / {int(((ok.rem_latency_min == ok.rem_latency_min_nm) | (ok.rem_latency_min.isna() & ok.rem_latency_min_nm.isna())).sum())} of {len(ok)}")
    L.append("\n## Limb vocabulary and durations per site (B3)")
    for site in ("I0002", "I0006"):
        s = ok[ok.site_id == site]
        hist = {}
        for h in s.limb_code_hist.dropna():
            for k, v in json.loads(h).items():
                hist[k] = hist.get(k, 0) + int(v)
        em = s.limb_event_map.dropna().unique().tolist()
        L.append(f"- {site} (n={len(s)}): code counts {dict(sorted(hist.items(), key=lambda kv: -kv[1]))}")
        L.append(f"    event_map variants: {len(em)}; first: {em[0][:300] if em else None}")
        L.append(f"    limb_rows {q(s.lm_rows)}; rows with valid 0.5-10 s duration {q(s.lm_rows_dur_valid)}; file median duration {q(s.limb_dur_median_file)}; "
                 f"frac positive durations {q(s.limb_dur_frac_positive)}; frac exactly 0.5 s {q(s.limb_dur_frac_le_0p5)}")
        L.append(f"    LM index {q(s.lm_index)}; PLM index (series rule) {q(s.plm_index)}; PLM index excl. resp {q(s.plm_index_excl_resp)}; "
                 f"nights with PLMI >= 15: {int((pd.to_numeric(s.plm_index, errors='coerce') >= 15).sum())}")
    L.append("\n## Per-site distribution of every new measure (gate 4)")
    for c in NEW:
        if c in ok.columns:
            L.append(f"- {c}: I0002 {q(ok[ok.site_id=='I0002'][c])} | I0006 {q(ok[ok.site_id=='I0006'][c])}")
    L.append(f"- hb_event_table: {ok.groupby('site_id').hb_event_table.agg(lambda s: s.value_counts().to_dict()).to_dict()}")
    L.append(f"- hb_window_source: {ok.hb_window_source.value_counts().to_dict()}")
    for site in ("I0002", "I0006"):
        hist = {}
        for h in ok[ok.site_id == site].resp_code_hist.dropna():
            for k, v in json.loads(h).items():
                hist[k] = hist.get(k, 0) + int(v)
        L.append(f"- resp code histogram {site}: {dict(sorted(hist.items(), key=lambda kv: -kv[1]))}")
    # ensemble curves
    calib = {"pilot_jsonl": a.jsonl, "n_records": int(len(ok)), "per_site": {}}
    L.append("\n## Event-locked ensemble SpO2 (per site), window proposal")
    for site in ("I0002", "I0006"):
        e = ensemble(ok, site)
        if e:
            calib["per_site"][site] = e
            L.append(f"- {site}: {e['n_nights']} nights, {e['n_events']} events; pre-max at {e['pre_max_offset_sec']:+.0f} s ({e['spo2_at_pre_max']:.2f}%), "
                     f"nadir at {e['nadir_offset_sec']:+.0f} s ({e['spo2_at_nadir']:.2f}%), return at {e['post_return_offset_sec']:+.0f} s, post-max at {e['post_max_offset_sec']:+.0f} s "
                     f"-> proposed W_pre {e['proposed_w_pre']:.0f} s, W_post {e['proposed_w_post']:.0f} s")
            offs = np.array(e["offsets_sec"]); ms = e["mean_spo2"]
            L.append("    curve (s: mean SpO2): " + " ".join(f"{int(o)}:{m:.1f}" for o, m in zip(offs, ms) if m is not None and int(o) % 20 == 0))
    json.dump(calib, open(a.calib, "w"), indent=1)
    L.append("\n## Timing")
    L.append(f"- total_sec per night: median {pd.to_numeric(ok.total_sec).median():.1f}, p95 {pd.to_numeric(ok.total_sec).quantile(.95):.1f}, max {pd.to_numeric(ok.total_sec).max():.1f}")
    rep = "\n".join(L); open(a.report, "w").write(rep + "\n"); print(rep)
    idx = f"{P.X1_LOGS}/PROVENANCE_INDEX.tsv"
    write_sidecar(a.report, [a.jsonl, P.REEXTRACT_FINAL, P.CPAP_STAGE_V7, P.SAMPLE_300], __file__, "020_check_pilot_light", index_path=idx)
    write_sidecar(a.calib, [a.jsonl], __file__, "020_check_pilot_light", index_path=idx)
    n_fail = sum(len(v) for v in fails.values())
    print(f"\nPOSITIVE-CONTROL MISMATCHES: {n_fail}")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
