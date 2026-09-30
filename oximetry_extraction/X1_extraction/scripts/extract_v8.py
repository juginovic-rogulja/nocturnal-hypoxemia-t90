"""
v8 extraction of one night (or a manifest of nights) on the untreated pre-PAP window.

extract(path, t_end_sec=None, families=("full",), fileobj=None, hb_def=None)
  full  : the v7 record (extract_night's diagnostics + the four April families + lung_variants + architecture, all through
          PrePapReader) plus the new measures. Needs a local .h5.
  light : architecture + oxygen family (April-style and native) + new measures from targeted reads (spo2, the stage,
          resp/resp_3/resp_4 and limb tables); works over s3fs without downloading the file.
  ecg   : optional heart family (not run by default; decision D7).
New measures per night: hypoxic_burden (+ hb_* diagnostics and the event-locked curve sums for the calibration), plm_index,
lm_index (+ limb vocabulary), sol_min / rem_latency_min (from the architecture; a second copy from new_measures.latencies as
sol_min_nm / rem_latency_min_nm), spo2_pct_below_90_sleep (+ sleep-period oxygen family), and the window provenance
(_v8_prepap_window, _v8_window_end_min, _v8_recording_dur_min_full).

usage: extract_v8.py --manifest M --families full|light [--t-end none|column] [--source local|s3fs] [--local-dir DIR]
                     [--workers N] [--night-timeout S] [--out out.jsonl] [--resume-from out.jsonl] [--limit N] [--shard i/n]
       extract_v8.py --one <h5> [--t-end-sec S] --families full --json out.json      (the EC2 worker's per-night call)
"""
from __future__ import annotations
import argparse, json, os, signal, sys, time, traceback
from pathlib import Path
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
REEX = os.path.join(os.path.dirname(HERE), "reex_code", "scripts")
if REEX not in sys.path:
    sys.path.insert(0, REEX)
import x1_paths as P  # noqa: E402
from prepap_reader import PrePapReader  # noqa: E402
import new_measures as NM  # noqa: E402
from extract_night import diagnostics, architecture, LEGEND, SLEEP_CODES  # noqa: E402
from reex_pipeline.core.h5_reader import REGION_CHANNELS  # noqa: E402
from reex_pipeline.features.brain_features import extract_brain_features  # noqa: E402
from reex_pipeline.features.heart_features import extract_heart_features  # noqa: E402
from reex_pipeline.features.lung_features import extract_lung_features, spo2_summary, count_desats  # noqa: E402
from reex_pipeline.features.microstructure import extract_microstructure  # noqa: E402

HB_DEFAULT = {"w_pre": 100.0, "w_post": 100.0, "baseline_sec": 100.0, "source": "default_100_100 (no hb_definition.json)"}
HB_CURVE_HALF_WIDTH = 180   # seconds each side of the event end for the ensemble calibration


def load_hb_definition(path=None):
    path = path or P.HB_DEFINITION
    if path and os.path.isfile(path):
        d = json.load(open(path))
        return {"per_site": d.get("per_site", {}), "baseline_sec": float(d.get("baseline_sec", 100.0)),
                "source": f"{os.path.basename(path)} sha256:{d.get('_self_sha256_prefix', '?')}", "version": d.get("version")}
    return None


def hb_windows_for(site, hb_def):
    if hb_def and site in hb_def.get("per_site", {}):
        w = hb_def["per_site"][site]
        return float(w["w_pre"]), float(w["w_post"]), hb_def["baseline_sec"], hb_def["source"]
    return HB_DEFAULT["w_pre"], HB_DEFAULT["w_post"], HB_DEFAULT["baseline_sec"], HB_DEFAULT["source"]


def lung_variants_v8(r: PrePapReader) -> dict:
    """v7 extract_night.lung_variants with the direct dataset read replaced by reader.read_native (cut-aware)."""
    out = {}
    if r.has_channel("spo2"):
        info = r.channel_info("spo2")
        out["spo2_native_fs"] = info["fs"]; out["spo2_unit"] = info["unit"]
        x = r.read_raw("spo2")
        if np.nanmedian(x[x > 0]) < 2.0:
            x = x * 100
        s = spo2_summary(x)
        out["spo2_nadir"] = s.get("spo2_nadir", np.nan)
        xv = x[(x >= 50.0) & (x <= 100.0)]
        out["spo2_p1_fft"] = float(np.percentile(xv, 1)) if xv.size else np.nan
        out["spo2_frac_below_50"] = float(np.mean(x < 50.0)) if x.size else np.nan
        raw = r.read_native("spo2")
        if np.nanmedian(raw[raw > 0]) < 2.0:
            raw = raw * 100
        sn = spo2_summary(raw)
        out["spo2_pct_below_90_native"] = sn.get("spo2_pct_below_90", np.nan)
        out["spo2_mean_native"] = sn.get("spo2_mean", np.nan)
        out["spo2_pct_below_88_native"] = sn.get("spo2_pct_below_88", np.nan)
        out["spo2_nadir_native"] = sn.get("spo2_nadir", np.nan)
        rv = raw[(raw >= np.float32(50.0)) & (raw <= np.float32(100.0))]
        out["spo2_p1_native"] = float(np.percentile(rv, 1)) if rv.size else np.nan
        out["spo2_native_valid_min"] = float(rv.size / info["fs"] / 60.0) if np.isfinite(info["fs"]) and info["fs"] > 0 else np.nan
        nfs = int(round(info["fs"])) if np.isfinite(info["fs"]) and info["fs"] >= 1 else None
        if nfs:
            out["odi3_total_native"] = float(count_desats(raw, nfs, 3.0))
            out["odi4_total_native"] = float(count_desats(raw, nfs, 4.0))
    for kind in ("resp_3", "resp_4"):
        if r.has_table(kind):
            arr = r.annotations_array(kind)
            out[f"resp_events_any_label_{kind}"] = int(np.sum(np.diff((arr != 0).astype(np.int8)) == 1))
    return out


def new_measures(r: PrePapReader, arch: dict, site: str, hb_def) -> dict:
    out = {}
    fs = r.fs
    stage = r._stage_array()
    tst_h = float(arch.get("TST_min", np.nan)) / 60.0
    # latencies, second copy from the pure function (must equal the architecture's)
    lat = NM.latencies(stage, fs)
    out["sol_min_nm"] = lat["sol_min"]; out["rem_latency_min_nm"] = lat["rem_latency_min"]
    # respiratory events: the table the v7 AHI came from (resp at I0006, resp_3 at I0002)
    resp_key = "resp" if r.has_table("resp") else ("resp_3" if r.has_table("resp_3") else None)
    out["hb_event_table"] = resp_key
    resp_int = None
    if resp_key:
        c, s, e = r.clean_events(resp_key)
        nz = c != 0
        s, e = s[nz], e[nz]
        resp_int = np.c_[s, e] if s.size else np.zeros((0, 2))
        out["hb_n_resp_rows_nonzero"] = int(s.size)
        em = r.event_map(resp_key)
        out["resp_code_hist"] = json.dumps({str(int(k)): int(v) for k, v in zip(*np.unique(c, return_counts=True))})
        out["resp_event_map"] = json.dumps(em, sort_keys=True) if em else None
    else:
        out["hb_n_resp_rows_nonzero"] = 0
    # SpO2 native + HB + sleep-period T90
    if r.has_channel("spo2"):
        info = r.channel_info("spo2")
        raw = r.read_native("spo2")
        if np.any(raw > 0) and np.nanmedian(raw[raw > 0]) < 2.0:
            raw = raw * 100
        nfs = float(info["fs"])
        w_pre, w_post, base_sec, src = hb_windows_for(site, hb_def)
        out["hb_window_source"] = src
        if resp_key is not None and np.isfinite(nfs) and nfs > 0:
            hb = NM.hypoxic_burden(raw, nfs, e if resp_key else [], w_pre, w_post, tst_h, baseline_sec=base_sec,
                                   curve_half_width=HB_CURVE_HALF_WIDTH)
            out.update(hb)
        else:
            out["hypoxic_burden"] = np.nan; out["hb_area_pctmin"] = np.nan
        if np.isfinite(nfs) and nfs > 0:
            out.update(NM.sleep_period_t90(raw, nfs, stage, fs))
    else:
        out["hypoxic_burden"] = np.nan; out["spo2_pct_below_90_sleep"] = np.nan
    # limb movements: the whole vocabulary counts, per-night histogram kept for the per-site convention check (B3)
    if r.has_table("limb"):
        c, s, e = r.clean_events("limb")
        nz = c != 0
        c, s, e = c[nz], s[nz], e[nz]
        codes_raw, _, _, durs_raw, em = r._event_table("limb")
        out["limb_code_hist"] = json.dumps({str(int(k)): int(v) for k, v in zip(*np.unique(codes_raw, return_counts=True))})
        out["limb_event_map"] = json.dumps(em, sort_keys=True) if em else None
        pos = durs_raw[durs_raw > 0]
        out["limb_dur_median_file"] = float(np.median(pos)) if pos.size else np.nan
        out["limb_dur_frac_positive"] = float(np.mean(durs_raw > 0)) if durs_raw.size else np.nan
        out["limb_dur_frac_le_0p5"] = float(np.mean(np.abs(pos - 0.5) < 1e-6)) if pos.size else np.nan
        d = e - s
        out.update(NM.plm_index(s, d, tst_h, resp_intervals=resp_int))
    else:
        out["plm_index"] = np.nan; out["lm_index"] = np.nan; out["lm_rows"] = 0
    return out


def extract(path, t_end_sec=None, families=("full",), fileobj=None, hb_def=None, site=None):
    t0 = time.time()
    fams = set(families)
    rec = {"source": str(path), "status": "ok", "_v8_families": ",".join(sorted(fams))}
    light = "full" not in fams
    try:
        with PrePapReader(path, t_end_sec=t_end_sec, fileobj=fileobj, light=light) as r:
            site = site or ("I0006" if "I0006" in str(path) else ("I0002" if "I0002" in str(path) else "unknown"))
            if not light:
                rec.update(diagnostics(r))
                brain = extract_brain_features(r); micro = extract_microstructure(r)
                heart = extract_heart_features(r); lung = extract_lung_features(r)
                rec["n_brain"], rec["n_micro"], rec["n_heart"], rec["n_lung"] = len(brain), len(micro), len(heart), len(lung)
                rec.update(brain); rec.update(micro); rec.update(heart); rec.update(lung)
            else:
                rec["legend"] = json.dumps(r.stage_event_map(), sort_keys=True)
                rec["legend_matches"] = (r.stage_event_map() == LEGEND)
                rec["ann_groups"] = ",".join(sorted(r._f["annotations/expert_1"].keys()))
                rec["duration_sec_attr"] = float(r._f.attrs.get("duration_sec", np.nan))
                rec["n_samples"] = int(r.n_samples)
                rec.update(extract_lung_features(r))
                if "ecg" in fams:
                    try:
                        rec.update(extract_heart_features(r))
                    except KeyError:
                        pass
            rec.update(lung_variants_v8(r))
            arch = architecture(r)
            rec.update(arch)
            rec.update(new_measures(r, arch, site, hb_def))
            rec["_v8_prepap_window"] = r.t_end_sec is not None
            rec["_v8_window_end_min"] = (r.window_end_sec / 60.0) if r.t_end_sec is not None else np.nan
            rec["_v8_recording_dur_min_full"] = r.n_samples_full / r.fs / 60.0
            rec.update({f"diag_{k}": (json.dumps(v) if isinstance(v, (list, dict)) else v) for k, v in r.diag.items()})
    except Exception as e:
        rec["status"] = f"fail:{type(e).__name__}:{str(e)[:200]}"
        rec["traceback"] = traceback.format_exc()[-1500:]
    rec["runtime_sec"] = round(time.time() - t0, 2)
    return rec


# ------------------------------------------------------------------------------------------ manifest runner
class NightTimeout(Exception):
    pass


def _alarm(signum, frame):
    raise NightTimeout()


_HB_DEF = None


def _init_worker(hb_path):
    global _HB_DEF
    _HB_DEF = load_hb_definition(hb_path)


def one_night(job):
    row, fams, source, local_dir, timeout = job
    site, bids, ses = row["site_id"], row["bids_col"], int(row["SessionID"])
    stem = f"{bids}_ses-{ses}"
    t_end = row.get("t_end_sec", None)
    try:
        t_end = float(t_end) if t_end not in (None, "", "nan") and np.isfinite(float(t_end)) else None
    except Exception:
        t_end = None
    rec = {"BDSPPatientID": int(row["BDSPPatientID"]), "site_id": site, "stem": stem, "_encoding_used": row.get("_encoding_used", "")}
    t0 = time.time()
    signal.signal(signal.SIGALRM, _alarm); signal.alarm(int(timeout))
    try:
        if source == "local":
            path = os.path.join(local_dir, f"{stem}.h5")
            if not os.path.isfile(path):
                rec["status"] = "missing_local_file"; return rec
            rec.update(extract(path, t_end_sec=t_end, families=fams, hb_def=_HB_DEF, site=site))
        else:
            import s3fs
            fs = s3fs.S3FileSystem(anon=False, profile=os.environ.get("X1_AWS_PROFILE") or None)
            url = f"{P.S3_AP}/PSG/bids/{site}/{bids}/ses-{ses}/eeg/{stem}.h5"
            with fs.open(url, "rb", block_size=1024 * 1024) as fo:
                rec.update(extract(url, t_end_sec=t_end, families=fams, fileobj=fo, hb_def=_HB_DEF, site=site))
    except NightTimeout:
        rec["status"] = f"timeout_{timeout}s"
    except Exception as e:
        rec["status"] = f"fail:{type(e).__name__}:{str(e)[:200]}"
    finally:
        signal.alarm(0)
    rec["total_sec"] = round(time.time() - t0, 1)
    return rec


def run_manifest(a):
    import pandas as pd
    from multiprocessing import Pool
    df = pd.read_csv(a.manifest, low_memory=False)
    if a.shard:
        i, n = [int(x) for x in a.shard.split("/")]; df = df.iloc[i::n]
    if a.limit:
        df = df.head(a.limit)
    if a.t_end == "none" or a.t_end is None:
        df["t_end_sec"] = np.nan
    elif a.t_end != "t_end_sec":
        df["t_end_sec"] = df[a.t_end]
    done = set()
    if a.resume_from and os.path.isfile(a.resume_from):
        for line in open(a.resume_from):
            try:
                d = json.loads(line)
                if str(d.get("status", "")).startswith("ok"):
                    done.add(int(d["BDSPPatientID"]))
            except Exception:
                pass
    rows = [r for r in df.to_dict("records") if int(r["BDSPPatientID"]) not in done]
    fams = tuple(a.families.split(","))
    jobs = [(r, fams, a.source, a.local_dir, a.night_timeout) for r in rows]
    print(f"extract_v8: {len(rows)} nights ({len(done)} already done), families={fams}, source={a.source}, "
          f"workers={a.workers}, timeout={a.night_timeout}s, hb_def={a.hb_definition}", flush=True)
    t0 = time.time(); n_ok = 0
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "a") as fh, Pool(a.workers, initializer=_init_worker, initargs=(a.hb_definition,), maxtasksperchild=50) as pool:
        for i, rec in enumerate(pool.imap_unordered(one_night, jobs, chunksize=1), 1):
            fh.write(json.dumps(rec, default=str) + "\n"); fh.flush()
            n_ok += str(rec.get("status", "")).startswith("ok")
            if i % a.report_every == 0 or i == len(rows):
                el = time.time() - t0
                print(f"[{i}/{len(rows)}] {el/60:.1f} min ok={n_ok} {el/i:.2f} s/night wall", flush=True)
    print(f"DONE {n_ok} ok of {len(rows)} in {(time.time()-t0)/60:.1f} min -> {a.out}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest"); ap.add_argument("--one"); ap.add_argument("--json")
    ap.add_argument("--families", default="full"); ap.add_argument("--t-end", default="t_end_sec",
                    help="'none' or the manifest column holding the cut in seconds (default t_end_sec)")
    ap.add_argument("--t-end-sec", type=float, default=None, help="with --one: the cut in seconds")
    ap.add_argument("--source", default="local", choices=["local", "s3fs"]); ap.add_argument("--local-dir", default="")
    ap.add_argument("--workers", type=int, default=4); ap.add_argument("--night-timeout", type=int, default=1500)
    ap.add_argument("--out"); ap.add_argument("--resume-from"); ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--shard", default=""); ap.add_argument("--report-every", type=int, default=25)
    ap.add_argument("--hb-definition", default=P.HB_DEFINITION if os.path.isfile(P.HB_DEFINITION) else "")
    a = ap.parse_args()
    if a.one:
        hb = load_hb_definition(a.hb_definition) if a.hb_definition else None
        rec = extract(a.one, t_end_sec=a.t_end_sec, families=tuple(a.families.split(",")), hb_def=hb)
        if a.json:
            Path(a.json).write_text(json.dumps(rec, indent=1, default=str))
        for k in ("status", "runtime_sec", "recording_dur_min", "_v8_window_end_min", "TST_min", "sleep_efficiency_pct", "AHI",
                  "spo2_pct_below_90", "spo2_pct_below_90_native", "spo2_pct_below_90_sleep", "hypoxic_burden", "plm_index",
                  "sol_min", "rem_latency_min", "F_nrem_delta_mean", "wholenight_rmssd"):
            print(f"{k:28s} {rec.get(k)}")
        return
    if not (a.manifest and a.out):
        ap.error("--manifest and --out are required (or --one)")
    run_manifest(a)


if __name__ == "__main__":
    main()
