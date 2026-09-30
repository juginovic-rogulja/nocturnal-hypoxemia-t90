"""
Step 021: freeze the hypoxic-burden definition (decision D5) in work/hb_definition.json from the pilot calibration
(work/hb_window_calibration.json, step 020). Per site: the event table (the v7 AHI table), the baseline rule, W_pre and
W_post from the ensemble curve (rounded up to whole 10 s, floored at 30 s and capped at 180 s = the calibration half-width),
the overlap rule, the denominator. The owner sees the curve and this file before the full pass. --fixed 100 100 writes the
alternative (fixed windows) instead.
"""
from __future__ import annotations
import argparse, hashlib, json, math, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import x1_paths as P
from provenance import write_sidecar


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--calib", default=P.HB_CALIBRATION); ap.add_argument("--out", default=P.HB_DEFINITION)
    ap.add_argument("--fixed", nargs=2, type=float, default=None, help="W_pre W_post for the fixed-window alternative")
    a = ap.parse_args()
    per_site = {}
    if a.fixed:
        for site in ("I0002", "I0006"):
            per_site[site] = {"w_pre": a.fixed[0], "w_post": a.fixed[1], "source": "fixed (alternative of D5)"}
        calib_sha = None
    else:
        P.assert_not_evicted(a.calib)
        cal = json.load(open(a.calib))
        calib_sha = hashlib.sha256(open(a.calib, "rb").read()).hexdigest()
        for site, e in cal["per_site"].items():
            wp = min(180.0, max(30.0, 10.0 * math.ceil(e["proposed_w_pre"] / 10.0)))
            wq = min(180.0, max(30.0, 10.0 * math.ceil(e["proposed_w_post"] / 10.0)))
            per_site[site] = {"w_pre": wp, "w_post": wq, "ensemble_pre_max_offset_sec": e["pre_max_offset_sec"], "ensemble_nadir_offset_sec": e["nadir_offset_sec"],
                              "ensemble_post_return_offset_sec": e["post_return_offset_sec"], "ensemble_n_nights": e["n_nights"], "ensemble_n_events": e["n_events"],
                              "source": "ensemble-averaged event-locked SpO2 of the 300-night pilot (step 020)"}
    d = {"version": "v8_2026-09-12", "decision": "D5 (RERUN_PLAN_V8.md)", "reference": "Azarbarzin et al., Eur Heart J 2019;40:1149",
         "event_table": {"I0002": "resp_3 (3 percent or arousal hypopnea rule, the v7 AHI table)", "I0006": "resp (single table, RERA rows included, the v7 AHI table)"},
         "event_rows": "every row with a non-zero code in that table (apnea, hypopnea, RERA alike), after the reader's day-wrap and negative-duration rules",
         "anchor": "event END", "baseline": "maximum valid SpO2 (50 to 100 percent) in the 100 s before the event end",
         "baseline_sec": 100.0,
         "window": "[end - W_pre, end + W_post] per site; the pre-window starts no earlier than the previous event's end; the post-window ends no later than where the next event's pre-window starts and never before the event's own end (windows are disjoint, no sample counted twice)",
         "area": "sum over valid native samples of max(0, baseline - SpO2), rectangle rule at the native oximeter rate, in percent-minutes",
         "denominator": "true sleep hours (TST_min / 60, sample-level, pre-PAP window on split nights)",
         "unit": "percent-minutes per hour of sleep", "spo2_source": "signals/spo2 at its native rate (no resampling)",
         "collapsed_staging_rule": "D11: NaN where _v7_1_stage_collapsed (TST undefined)",
         "per_site": per_site, "calibration_file": a.calib if not a.fixed else None, "calibration_sha256": calib_sha,
         "frozen_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    body = json.dumps(d, indent=1, sort_keys=True)
    d["_self_sha256_prefix"] = hashlib.sha256(body.encode()).hexdigest()[:16]
    json.dump(d, open(a.out, "w"), indent=1, sort_keys=True)
    write_sidecar(a.out, [a.calib] if not a.fixed else [], __file__, "021_freeze_hb_definition", index_path=f"{P.X1_LOGS}/PROVENANCE_INDEX.tsv")
    print(json.dumps({k: v for k, v in d.items() if k in ("per_site", "window", "baseline", "denominator", "_self_sha256_prefix")}, indent=1))
    print("wrote", a.out)


if __name__ == "__main__":
    main()
