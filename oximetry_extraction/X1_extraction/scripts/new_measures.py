"""
The four new measurements of decision 3 (plus the sleep-period T90 of decision 2), as pure numpy functions on arrays so
that unit tests, the extraction pass and the independent second path (step 024) call the same arithmetic on plain inputs.

hypoxic_burden   Azarbarzin 2019 (decision D5): per respiratory event, baseline = max valid SpO2 in the 100 s before the
                 event END; window [end - W_pre, end + W_post]; area of max(0, baseline - SpO2) over VALID native samples
                 (50 to 100 percent), rectangle rule at the native rate, percent-minutes; windows of neighbouring events
                 never overlap: the pre-window starts no earlier than the previous event's end, the post-window ends no
                 later than where the next event's pre-window starts (and never before the event's own end); summed over
                 events and divided by true sleep hours -> percent-minutes per hour.
plm_index        AASM series rule (decision D6): limb movements 0.5 to 10 s long; onsets within `merge_sec` (5 s) of the
                 previous kept onset are one movement (bilateral rule); a series = 4 or more consecutive movements with
                 onset-to-onset intervals of 5 to 90 s; PLM count = movements inside series; index = count / sleep hours.
                 Variant excluding movements within 0.5 s of a respiratory event (AASM) when resp intervals are given.
latencies        sample-level, the v7 architecture rule: SOL = first sleep sample / fs / 60; REM latency = first REM sample
                 after sleep onset, minutes from onset.
sleep_period_t90 native SpO2 samples that fall inside sleep-coded samples of the 200 Hz stage vector (decision 2).
"""
from __future__ import annotations
import numpy as np

SLEEP_CODES = (1, 2, 3, 4)
CODE_REM, CODE_WAKE = 4, 5
VALID_LO, VALID_HI = np.float32(50.0), np.float32(100.0)


def valid_mask(x32: np.ndarray) -> np.ndarray:
    """The float32 validity filter of the v7 lung code (a 1 Hz export stores 100 as 100.0000x in float64)."""
    return (x32 >= VALID_LO) & (x32 <= VALID_HI)


def hypoxic_burden(spo2: np.ndarray, fs: float, event_ends, w_pre: float, w_post: float, tst_hours: float,
                   baseline_sec: float = 100.0, curve_half_width: int = 0):
    """Returns dict with hb_area_pctmin (sum over events), hypoxic_burden (per sleep hour), event counts and, when
    curve_half_width > 0, the per-second event-locked sums (mean SpO2 and deficit) for the ensemble calibration."""
    x32 = np.asarray(spo2, dtype=np.float32)
    x = x32.astype(np.float64)
    v = valid_mask(x32)
    n = x.size
    T = n / fs
    ends = np.sort(np.asarray(list(event_ends), dtype=float))
    ends = ends[(ends > 0.0) & (ends <= T + 1.0)]
    total = 0.0; used = 0; nobase = 0; empty = 0
    K = int(curve_half_width)
    if K:
        offs = np.arange(-K, K + 1)
        c_sum = np.zeros(offs.size); c_cnt = np.zeros(offs.size); d_sum = np.zeros(offs.size)
    for i, e in enumerate(ends):
        prev_end = ends[i - 1] if i > 0 else -np.inf
        next_end = ends[i + 1] if i + 1 < ends.size else np.inf
        a = max(0, int(np.floor((e - baseline_sec) * fs))); b = min(n, int(np.ceil(e * fs)))
        seg = x[a:b][v[a:b]]
        if seg.size == 0:
            nobase += 1
            continue
        base = float(seg.max())
        ws = max(e - w_pre, prev_end)
        we = max(e, min(e + w_post, next_end - w_pre))
        ia = max(0, int(round(ws * fs))); ib = min(n, int(round(we * fs)))
        used += 1
        if ib > ia:
            seg = x[ia:ib]; vv = v[ia:ib]
            deficit = np.where(vv, base - seg, 0.0)
            deficit[deficit < 0] = 0.0
            total += float(deficit.sum()) / fs / 60.0
        else:
            empty += 1
        if K:
            idx = np.round((e + offs) * fs).astype(np.int64)
            ok = (idx >= 0) & (idx < n)
            ok[ok] &= v[idx[ok]]
            c_sum[ok] += x[idx[ok]]; c_cnt[ok] += 1
            d_sum[ok] += np.maximum(base - x[idx[ok]], 0.0)
    out = {"hb_area_pctmin": total, "hypoxic_burden": (total / tst_hours) if (tst_hours and tst_hours > 0) else np.nan,
           "hb_n_events_in_window": int(ends.size), "hb_n_events_used": int(used), "hb_n_events_no_baseline": int(nobase),
           "hb_n_events_empty_window": int(empty), "hb_w_pre": float(w_pre), "hb_w_post": float(w_post)}
    if K:
        out["hb_curve_offsets_sec"] = offs.tolist(); out["hb_curve_spo2_sum"] = c_sum.tolist()
        out["hb_curve_deficit_sum"] = d_sum.tolist(); out["hb_curve_count"] = c_cnt.tolist()
    return out


def plm_index(onsets, durations, tst_hours: float, resp_intervals=None, min_dur: float = 0.5, max_dur: float = 10.0,
              min_gap: float = 5.0, max_gap: float = 90.0, min_series: int = 4, merge_sec: float = 5.0):
    on = np.asarray(list(onsets), float); du = np.asarray(list(durations), float)
    order = np.argsort(on, kind="stable"); on = on[order]; du = du[order]
    out = {"lm_rows": int(on.size)}
    dur_ok = (du >= min_dur) & (du <= max_dur)
    out["lm_rows_dur_valid"] = int(dur_ok.sum())
    on_v = on[dur_ok]
    # bilateral / duplicate merge: an onset within merge_sec after the previous kept onset is the same movement
    kept = []
    for t in on_v:
        if kept and (t - kept[-1]) < merge_sec:
            continue
        kept.append(t)
    kept = np.asarray(kept, float)
    out["lm_count"] = int(kept.size)

    def series_count(ts):
        if ts.size < min_series:
            return 0, 0
        gaps = np.diff(ts)
        ok = (gaps >= min_gap) & (gaps <= max_gap)
        cnt = 0; nser = 0; run = 1
        for g in ok:
            if g:
                run += 1
            else:
                if run >= min_series:
                    cnt += run; nser += 1
                run = 1
        if run >= min_series:
            cnt += run; nser += 1
        return cnt, nser

    plm, nser = series_count(kept)
    h = tst_hours if (tst_hours and tst_hours > 0) else np.nan
    out.update({"plm_count": int(plm), "plm_n_series": int(nser), "plm_index": plm / h, "lm_index": kept.size / h})
    if resp_intervals is not None:
        ri = np.asarray(list(resp_intervals), float).reshape(-1, 2)
        if ri.size and kept.size:
            s = ri[:, 0] - 0.5; e = ri[:, 1] + 0.5
            assoc = np.zeros(kept.size, bool)
            for a, b in zip(s, e):
                assoc |= (kept >= a) & (kept <= b)
            kept2 = kept[~assoc]
        else:
            kept2 = kept
        plm2, nser2 = series_count(kept2)
        out.update({"lm_count_excl_resp": int(kept2.size), "plm_count_excl_resp": int(plm2), "plm_index_excl_resp": plm2 / h})
    return out


def latencies(stage: np.ndarray, fs: float):
    """stage = sample-level stage vector at fs (codes 1 N3, 2 N2, 3 N1, 4 REM, 5 Wake, 0 unscored)."""
    st = np.asarray(stage)
    sleep = np.where(np.isin(st, SLEEP_CODES))[0]
    if sleep.size == 0:
        return {"sol_min": np.nan, "rem_latency_min": np.nan}
    onset = int(sleep[0])
    rem = np.where(st[onset:] == CODE_REM)[0]
    return {"sol_min": onset / fs / 60.0, "rem_latency_min": (rem[0] / fs / 60.0) if rem.size else np.nan}


def sleep_period_t90(spo2_native: np.ndarray, fs_native: float, stage200: np.ndarray, fs200: float):
    """Native samples inside sleep-coded samples. The stage of native sample i is the 200 Hz stage at floor(i * fs200/fs_native)."""
    x32 = np.asarray(spo2_native, dtype=np.float32)
    st = np.asarray(stage200)
    n = x32.size; n200 = st.size
    if n == 0 or n200 == 0 or not (fs_native > 0):
        return {"spo2_pct_below_90_sleep": np.nan, "spo2_pct_below_88_sleep": np.nan, "spo2_mean_sleep": np.nan,
                "spo2_p1_sleep": np.nan, "spo2_sleep_valid_min": 0.0, "spo2_sleep_n_valid": 0}
    idx = np.minimum(np.floor(np.arange(n) * (fs200 / fs_native)).astype(np.int64), n200 - 1)
    sleep = np.isin(st[idx], SLEEP_CODES)
    v = valid_mask(x32) & sleep
    xv = x32[v]
    if xv.size == 0:
        return {"spo2_pct_below_90_sleep": np.nan, "spo2_pct_below_88_sleep": np.nan, "spo2_mean_sleep": np.nan,
                "spo2_p1_sleep": np.nan, "spo2_sleep_valid_min": 0.0, "spo2_sleep_n_valid": 0}
    return {"spo2_pct_below_90_sleep": float(np.mean(xv < np.float32(90.0)) * 100),
            "spo2_pct_below_88_sleep": float(np.mean(xv < np.float32(88.0)) * 100),
            "spo2_mean_sleep": float(np.mean(xv.astype(np.float64))), "spo2_p1_sleep": float(np.percentile(xv, 1)),
            "spo2_sleep_valid_min": float(xv.size / fs_native / 60.0), "spo2_sleep_n_valid": int(xv.size)}
