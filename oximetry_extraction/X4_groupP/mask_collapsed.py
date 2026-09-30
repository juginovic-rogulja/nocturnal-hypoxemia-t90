#!/usr/bin/env python3
"""Collapsed-staging rule for the group P re-extraction (v8.2, 2026-09-15). Rule D11 of the v8 plan, applied at
assembly from the extraction's own per-state minutes, never from a list carried over from an earlier extraction:

    a night with positive sleep whose sleep sits 100 percent in ONE sleep code (N1, N2, N3 or REM) has an impossible
    composition; its stage-specific values (nrem, n1, n2, n3, rem) are set missing and the night is flagged
    _v8_stage_collapsed = 1; sleep-versus-wake totals (wake, sleep) and whole-window values (all/any, rec, deciles) stay.

Formats
  cpap   : csv from pilots/cpap_stage/extract.py, columns {tag}_{state}_{t90,n,min,mean,nadir}, tag in all/pre/post/h1/h2.
           Collapse is judged on the whole recording (all_*_min).
  t90    : csv from pilots/t90_by_stage/extract.py, columns t90_{state}, min_{state} (+ wakepre/wakepost kept).
           Collapse is judged on min_n1/min_n2/min_n3/min_rem (the PAP-off portion, the only staging that script keeps).
  oxy    : jsonl from oxygen_profile/worker_oxyprofile.py, keys {state}_{mean,min,n,t90,valid_min,...} for
           wake/n1/n2/n3/rem/sleep/rec. Collapse is judged on n1_min/n2_min/n3_min/rem_min.
usage: mask_collapsed.py --format cpap|t90|oxy --in <file> --out <file> [--report <json>]
Writes the masked copy; the unmasked input is kept where it was (raw/). Prints the counts.
"""
import argparse, json, math
import numpy as np, pandas as pd

SLEEP = ("n1", "n2", "n3", "rem")
STAGE_ONLY = ("nrem", "n1", "n2", "n3", "rem")


def collapsed_from_minutes(mins):
    """mins: dict state -> minutes (may be NaN). Positive sleep in exactly one sleep code."""
    v = {s: (0.0 if (m is None or (isinstance(m, float) and math.isnan(m))) else float(m)) for s, m in mins.items()}
    pos = [s for s in SLEEP if v.get(s, 0.0) > 0.0]
    return len(pos) == 1 and sum(v.values()) > 0.0


def mask_cpap(df):
    tags = sorted({c.split("_")[0] for c in df.columns if c.split("_")[0] in ("all", "pre", "post", "h1", "h2")})
    flag = df.apply(lambda r: collapsed_from_minutes({s: r.get(f"all_{s}_min", np.nan) for s in SLEEP}), axis=1)
    flag &= df.get("status", pd.Series("ok", index=df.index)).eq("ok")
    for tag in tags:
        for s in STAGE_ONLY:
            for k in ("t90", "n", "min", "mean", "nadir"):
                c = f"{tag}_{s}_{k}"
                if c in df.columns:
                    df.loc[flag, c] = np.nan
    df["_v8_stage_collapsed"] = flag.astype(int)
    return df, int(flag.sum())


def mask_t90(df):
    flag = df.apply(lambda r: collapsed_from_minutes({s: r.get(f"min_{s}", np.nan) for s in SLEEP}), axis=1)
    flag &= df.get("status", pd.Series("ok", index=df.index)).eq("ok")
    for s in STAGE_ONLY:
        for c in (f"t90_{s}", f"min_{s}", f"n_{s}", f"mean_{s}", f"nadir_{s}"):
            if c in df.columns:
                df.loc[flag, c] = np.nan
    df["_v8_stage_collapsed"] = flag.astype(int)
    return df, int(flag.sum())


def mask_oxy_record(r):
    if r.get("status") != "ok":
        r["_v8_stage_collapsed"] = 0
        return r, False
    flag = collapsed_from_minutes({s: r.get(f"{s}_min") for s in SLEEP})
    if flag:
        for k in list(r.keys()):
            if k.split("_")[0] in STAGE_ONLY:
                r[k] = None
    r["_v8_stage_collapsed"] = int(flag)
    return r, flag


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--format", required=True, choices=("cpap", "t90", "oxy"))
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--report", default="")
    a = ap.parse_args()
    if a.format == "oxy":
        n = nc = 0
        with open(a.inp) as fi, open(a.out, "w") as fo:
            for line in fi:
                line = line.strip()
                if not line:
                    continue
                r = json.loads(line)
                r, f = mask_oxy_record(r)
                n += 1; nc += int(bool(f))
                fo.write(json.dumps(r) + "\n")
        rep = {"format": "oxy", "in": a.inp, "out": a.out, "records": n, "collapsed": nc}
    else:
        df = pd.read_csv(a.inp, low_memory=False)
        df, nc = (mask_cpap if a.format == "cpap" else mask_t90)(df)
        df.to_csv(a.out, index=False)
        rep = {"format": a.format, "in": a.inp, "out": a.out, "records": int(len(df)), "collapsed": nc}
    print(json.dumps(rep))
    if a.report:
        json.dump(rep, open(a.report, "w"), indent=1)


if __name__ == "__main__":
    main()
