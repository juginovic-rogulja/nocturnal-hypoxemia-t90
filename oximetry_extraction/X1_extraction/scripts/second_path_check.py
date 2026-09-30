"""
Step 024 (gate 1, table against source by a second code path): an INDEPENDENT numpy implementation, importing neither
reex_pipeline nor new_measures nor prepap_reader, of TST, SOL, REM latency, sleep-period T90, whole-recording native T90,
hypoxic burden and the PLM index on N stratified nights read again from S3. Every value must be within 1e-6 of the pass.
The rules re-implemented from the written definitions: the stage vector at 200 Hz painted in file order on the
recording-start grid (EEG derivation length), day-wrapped starts unwrapped when start - 86400 lies in [0, D], starts beyond
that dropped, an end at or before its start replaced by start + duration (or the table's median positive duration, or the
kind's default), the pre-PAP cut = samples before round(t_end * fs) and event rows starting before t_end, SpO2 validity
50 to 100 in float32, the HB definition of hb_definition.json (or 100/100 when absent), the AASM series rule.
usage: second_path_check.py --pass-jsonl <jsonl> --manifest <csv> [--n 50] [--seed 20260912] [--hb-definition path] [--out csv]
"""
from __future__ import annotations
import argparse, json, os, sys, time
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import x1_paths as P
from provenance import write_sidecar

FS = 200
DEFAULT_DUR = {"stage": 30.0, "arousal": 3.0, "resp": 10.0, "resp_3": 10.0, "resp_4": 10.0, "limb": 0.5}
EEG = ["f3-m2", "f4-m1", "c3-m2", "c4-m1", "o1-m2", "o2-m1"]
COMPARE = ["TST_min", "sol_min", "rem_latency_min", "spo2_pct_below_90_native", "spo2_pct_below_90_sleep", "hypoxic_burden", "plm_index", "lm_index"]


def scaled(d):
    a = dict(d.attrs); x = d[:].ravel().astype(np.float64)
    if all(k in a for k in ("dig_min", "dig_max", "phys_min", "phys_max")):
        g = (float(a["phys_max"]) - float(a["phys_min"])) / (float(a["dig_max"]) - float(a["dig_min"]))
        x = x * g + (float(a["phys_min"]) - float(a["dig_min"]) * g)
    return x.astype(np.float32), float(a.get("fs", np.nan))


def table(f, kind):
    g = f[f"annotations/expert_1/{kind}"]
    return g["codes"][:].astype(np.int64), g["starts"][:].astype(np.float64), g["ends"][:].astype(np.float64), g["durations"][:].astype(np.float64)


def clean(codes, starts, ends, durs, D, kind, t_end=None):
    pos = durs[durs > 0]; med = float(np.median(pos)) if pos.size else DEFAULT_DUR[kind]
    out = []
    for c, s, e, d in zip(codes, starts, ends, durs):
        if s > D + 1.0:
            if 0.0 <= s - 86400.0 <= D:
                s2 = s - 86400.0; e2 = e - 86400.0 if (e - 86400.0) > s2 else s2 + (d if d > 0 else med); s, e = s2, e2
            else:
                continue
        elif e <= s:
            e = s + (d if d > 0 else med)
        if t_end is not None and s >= t_end:
            continue
        out.append((c, s, e))
    return out


def paint(events, n):
    out = np.zeros(n, np.int16)
    for c, s, e in events:
        a = max(int(round(s * FS)), 0); b = min(int(round(e * FS)), n)
        if b > a:
            out[a:b] = c
    return out


def one(row, hb_def, t_end):
    import s3fs, h5py
    fs = s3fs.S3FileSystem(anon=False, profile=os.environ.get("X1_AWS_PROFILE") or None)
    site, bids, ses = row["site_id"], row["bids_col"], int(row["SessionID"]); stem = f"{bids}_ses-{ses}"
    url = f"{P.S3_AP}/PSG/bids/{site}/{bids}/ses-{ses}/eeg/{stem}.h5"
    with fs.open(url, "rb", block_size=1024 * 1024) as fo, h5py.File(fo, "r") as f:
        g = f["signals"]
        eeg = [c for c in EEG if c in g and np.isfinite(float(g[c].attrs.get("fs", np.nan)))]
        if eeg:
            n_full = max(int(round(g[c].shape[0] * FS / float(g[c].attrs["fs"]))) for c in eeg)
        else:   # the reader's fallback: duration_sec attr, else the longest channel rescaled
            d = f.attrs.get("duration_sec", None)
            if d is not None and np.isfinite(float(d)):
                n_full = int(round(float(d) * FS))
            else:
                n_full = max(int(round(g[c].shape[0] * FS / float(g[c].attrs["fs"]))) for c in g if np.isfinite(float(g[c].attrs.get("fs", np.nan))))
        n = n_full if t_end is None else min(n_full, int(round(t_end * FS)))
        D_full = n_full / FS; te = None if t_end is None else n / FS
        st = paint(clean(*table(f, "stage"), D_full, "stage"), n)      # painted on the full grid rules, then the array is the window
        spo2, nfs = scaled(g["spo2"])
        if te is not None:
            spo2 = spo2[: int(round(te * nfs))]
        rk = "resp" if "annotations/expert_1/resp" in f else "resp_3"
        resp = [(c, s, min(e, n / FS)) for c, s, e in clean(*table(f, rk), D_full, rk, te) if c != 0]
        limb = [(c, s, min(e, n / FS)) for c, s, e in clean(*table(f, "limb"), D_full, "limb", te) if c != 0] if "annotations/expert_1/limb" in f else []
    # architecture
    sleep = np.isin(st, [1, 2, 3, 4]); tst_min = sleep.sum() / FS / 60.0; tst_h = tst_min / 60.0
    idx = np.where(sleep)[0]
    if idx.size:
        on = idx[0]; sol = on / FS / 60.0; rem = np.where(st[on:] == 4)[0]; reml = rem[0] / FS / 60.0 if rem.size else np.nan
    else:
        sol = reml = np.nan
    # oxygen
    if np.any(spo2 > 0) and np.nanmedian(spo2[spo2 > 0]) < 2.0:
        spo2 = spo2 * 100
    v = (spo2 >= np.float32(50)) & (spo2 <= np.float32(100))
    t90 = float(np.mean(spo2[v] < np.float32(90)) * 100) if v.any() else np.nan
    map200 = np.minimum(np.floor(np.arange(spo2.size) * (FS / nfs)).astype(np.int64), st.size - 1)
    vs = v & np.isin(st[map200], [1, 2, 3, 4])
    t90s = float(np.mean(spo2[vs] < np.float32(90)) * 100) if vs.any() else np.nan
    # hypoxic burden
    w = hb_def["per_site"][site] if hb_def and site in hb_def.get("per_site", {}) else {"w_pre": 100.0, "w_post": 100.0}
    bsec = float(hb_def.get("baseline_sec", 100.0)) if hb_def else 100.0
    ends = np.sort(np.array([e for _, _, e in resp], float)); T = spo2.size / nfs
    ends = ends[(ends > 0) & (ends <= T + 1.0)]
    x = spo2.astype(np.float64); area = 0.0
    for i, e in enumerate(ends):
        pe = ends[i - 1] if i else -np.inf; ne = ends[i + 1] if i + 1 < ends.size else np.inf
        a = max(0, int(np.floor((e - bsec) * nfs))); b = min(x.size, int(np.ceil(e * nfs)))
        seg = x[a:b][v[a:b]]
        if seg.size == 0:
            continue
        base = seg.max()
        ws = max(e - w["w_pre"], pe); we = max(e, min(e + w["w_post"], ne - w["w_pre"]))
        ia = max(0, int(round(ws * nfs))); ib = min(x.size, int(round(we * nfs)))
        if ib > ia:
            dfc = base - x[ia:ib]; dfc[~v[ia:ib]] = 0; dfc[dfc < 0] = 0; area += dfc.sum() / nfs / 60.0
    hb = area / tst_h if tst_h > 0 else np.nan
    # PLM index (AASM series rule)
    lm = sorted((s, e - s) for _, s, e in limb)
    lm = [s for s, d in lm if 0.5 <= d <= 10.0]
    kept = []
    for s in lm:
        if not kept or s - kept[-1] >= 5.0:
            kept.append(s)
    plm = 0; run = 1
    for i in range(1, len(kept)):
        gap = kept[i] - kept[i - 1]
        if 5.0 <= gap <= 90.0:
            run += 1
        else:
            plm += run if run >= 4 else 0; run = 1
    plm += run if run >= 4 else 0
    return {"BDSPPatientID": int(row["BDSPPatientID"]), "TST_min": tst_min, "sol_min": sol, "rem_latency_min": reml, "spo2_pct_below_90_native": t90,
            "spo2_pct_below_90_sleep": t90s, "hypoxic_burden": hb, "plm_index": plm / tst_h if tst_h > 0 else np.nan,
            "lm_index": len(kept) / tst_h if tst_h > 0 else np.nan, "n_resp": len(ends), "nfs": nfs}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--pass-jsonl", required=True); ap.add_argument("--manifest", required=True)
    ap.add_argument("--n", type=int, default=50); ap.add_argument("--seed", type=int, default=P.SEED)
    ap.add_argument("--hb-definition", default=P.HB_DEFINITION); ap.add_argument("--out", default=f"{P.X1_LOGS}/SECOND_PATH_50.csv")
    a = ap.parse_args()
    P.assert_not_evicted(a.pass_jsonl, a.manifest)
    hb_def = json.load(open(a.hb_definition)) if a.hb_definition and os.path.isfile(a.hb_definition) else None
    print("hb definition:", (a.hb_definition if hb_def else "none -> 100/100 default"), flush=True)
    rec = pd.DataFrame([json.loads(l) for l in open(a.pass_jsonl)]).drop_duplicates("BDSPPatientID", keep="last")
    rec = rec[rec.status == "ok"]
    man = pd.read_csv(a.manifest); man = man[man.BDSPPatientID.isin(rec.BDSPPatientID)].copy()
    man["is_split"] = man.t_end_sec.notna().astype(int) if "t_end_sec" in man else 0
    rng = np.random.default_rng(a.seed)
    parts = []
    for _, sub in man.groupby(["site_id", "is_split"]):
        k = max(1, int(round(a.n * len(sub) / len(man))))
        parts.append(sub.sample(n=min(k, len(sub)), random_state=int(rng.integers(0, 2**31 - 1))))
    pick = pd.concat(parts).head(a.n)
    rows = []; t0 = time.time()
    for i, (_, r) in enumerate(pick.iterrows(), 1):
        t_end = float(r.t_end_sec) if "t_end_sec" in r and pd.notna(r.t_end_sec) else None
        try:
            out = one(r, hb_def, t_end)
        except Exception as e:
            out = {"BDSPPatientID": int(r.BDSPPatientID), "error": f"{type(e).__name__}:{str(e)[:120]}"}
        p = rec[rec.BDSPPatientID == int(r.BDSPPatientID)].iloc[0]
        for c in COMPARE:
            out[f"{c}_pass"] = float(p[c]) if pd.notna(p[c]) else np.nan
            out[f"{c}_ok"] = (np.isnan(out.get(c, np.nan)) and np.isnan(out[f"{c}_pass"])) or (abs(out.get(c, np.nan) - out[f"{c}_pass"]) <= 1e-6)
        out["site_id"] = r.site_id; out["is_split"] = int(r.is_split); rows.append(out)
        if i % 10 == 0:
            print(f"[{i}/{len(pick)}] {(time.time()-t0)/60:.1f} min", flush=True)
    df = pd.DataFrame(rows); df.to_csv(a.out, index=False)
    okc = {c: int(df[f"{c}_ok"].sum()) for c in COMPARE}
    print(f"second path on {len(df)} nights (split {int(df.is_split.sum())}); within 1e-6:", okc)
    bad = df[~df[[f"{c}_ok" for c in COMPARE]].all(axis=1)]
    if len(bad):
        print("MISMATCHES:"); print(bad[["BDSPPatientID", "site_id", "is_split"] + sum([[c, f"{c}_pass"] for c in COMPARE if not df[f"{c}_ok"].all()], [])].to_string(index=False))
    write_sidecar(a.out, [a.pass_jsonl, a.manifest] + ([a.hb_definition] if hb_def else []), __file__, "024_second_path", index_path=f"{P.X1_LOGS}/PROVENANCE_INDEX.tsv")
    return 0 if not len(bad) else 1


if __name__ == "__main__":
    sys.exit(main())
