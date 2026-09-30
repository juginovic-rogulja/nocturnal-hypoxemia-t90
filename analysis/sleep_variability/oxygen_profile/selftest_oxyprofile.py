"""
Known-answer test of the two new blocks in worker_oxyprofile.py, plus the stage
placement and the affine they sit on.

The fleet validation proves the extraction reproduces frozen columns on real
recordings. It cannot prove the decile arithmetic or the per-stage arithmetic,
because no frozen decile or per-stage saturation column exists anywhere. This
proves both on waveforms whose answer is known by construction, and it imports
the shipped functions rather than transcribing them, so what is tested is what
runs on EC2.

    python3 selftest_oxyprofile.py

Cases
  1  ten flat deciles, distinct values, exact means and counts by construction
  2  the keep band [50,100] removes samples from a decile's DENOMINATOR
  3  a decile that is entirely dropout reports no mean, not a fabricated one
  4  a sample count not divisible by ten still partitions exactly
  5  the two identities the downstream gate checks, on every case above
  6  elapsed minutes per decile follow the sampling rate
  7  per-stage means, counts and minutes on a constructed hypnogram
  8  stage_on_grid places events on the SpO2 grid at the right samples
  9  the whole-day offset repair reproduces case 8 exactly
 10  the digital-to-physical affine, with a known linear map
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from worker_oxyprofile import (affine, decile_profile, stage_on_grid,  # noqa: E402
                               stage_profile, describe, N_DECILES)

F = []


def chk(name, got, want, tol=1e-9):
    ok = got is not None and abs(float(got) - float(want)) <= tol
    print(f"  {'PASS' if ok else 'FAIL'}  {name:<62} got {'None' if got is None else f'{float(got):12.6f}'}"
          f"  want {float(want):12.6f}")
    if not ok:
        F.append(name)


def chk_true(name, cond):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}")
    if not cond:
        F.append(name)


def run_deciles(x, fs=10.0):
    x = np.asarray(x, np.float64)
    valid = np.isfinite(x) & (x >= 50.0) & (x <= 100.0)
    out = {}
    describe(x[valid], "rec", out)
    decile_profile(x, valid, fs, out)
    return out


def identities(out, label):
    """The two identities the validation gate checks on every real record."""
    ns = [out.get(f"dec_n_{i+1:02d}", 0) for i in range(N_DECILES)]
    ms = [out.get(f"dec_mean_{i+1:02d}") for i in range(N_DECILES)]
    chk(f"{label}: decile counts sum to the whole-recording valid count",
        sum(ns), out.get("rec_n", 0))
    tot = sum(n for n in ns)
    if tot and out.get("rec_mean") is not None:
        w = sum(n * m for n, m in zip(ns, ms) if n) / tot
        chk(f"{label}: count-weighted decile mean equals the recording mean",
            w, out["rec_mean"], tol=1e-9)


print("1. ten flat deciles, 100 samples each, distinct values")
vals = [99.0, 98.0, 97.0, 96.0, 95.0, 94.0, 93.0, 92.0, 91.0, 90.0]
x1 = np.concatenate([np.full(100, v) for v in vals])
o1 = run_deciles(x1)
for i, v in enumerate(vals):
    chk(f"decile {i+1} mean is {v:.0f}", o1[f"dec_mean_{i+1:02d}"], v)
chk("decile 1 count", o1["dec_n_01"], 100)
chk("decile 10 count", o1["dec_n_10"], 100)
chk("deciles with data", o1["n_deciles_with_data"], 10)
chk_true("decile_status ok", o1["decile_status"] == "ok")

print("\n2. the keep band [50,100] drops samples from a decile's DENOMINATOR")
# decile 3 becomes 50 samples at 95 and 50 at 120 (out of band): mean 95, n 50
x2 = x1.copy()
x2[200:250] = 95.0
x2[250:300] = 120.0
o2 = run_deciles(x2)
chk("decile 3 mean over the surviving samples", o2["dec_mean_03"], 95.0)
chk("decile 3 count after the drop", o2["dec_n_03"], 50)
chk("decile 2 untouched", o2["dec_mean_02"], 98.0)
chk("decile 4 untouched", o2["dec_mean_04"], 96.0)

print("\n3. an all-dropout decile reports no mean at all")
x3 = x1.copy()
x3[500:600] = 20.0                      # decile 6 is entirely below the keep band
o3 = run_deciles(x3)
chk_true("decile 6 has no mean key", "dec_mean_06" not in o3)
chk("decile 6 count is zero", o3["dec_n_06"], 0)
chk("deciles with data drops to nine", o3["n_deciles_with_data"], 9)
chk_true("decile_status partial", o3["decile_status"] == "partial")

print("\n4. a sample count not divisible by ten still partitions exactly")
n4 = 1003
x4 = 50.0 + np.arange(n4) * (49.0 / n4)      # strictly inside the keep band
o4 = run_deciles(x4)
edges = np.rint(np.linspace(0, n4, 11)).astype(int)
for i in range(10):
    want = x4[edges[i]:edges[i + 1]].mean()
    chk(f"decile {i+1} mean matches an independent slice", o4[f"dec_mean_{i+1:02d}"], want)
chk("counts sum to every sample", sum(o4[f"dec_n_{i+1:02d}"] for i in range(10)), n4)
chk_true("spans are contiguous and non-overlapping",
         bool((np.diff(edges) > 0).all() and edges[0] == 0 and edges[-1] == n4))

print("\n5. the two identities the gate checks")
for lab, o in (("case1", o1), ("case2", o2), ("case3", o3), ("case4", o4)):
    identities(o, lab)

print("\n6. elapsed minutes per decile follow the sampling rate")
x6 = np.full(6000, 96.0)                      # 6000 samples at 10 Hz = 10 minutes
o6 = run_deciles(x6, fs=10.0)
for i in (1, 5, 10):
    chk(f"decile {i} covers one minute", o6[f"dec_min_{i:02d}"], 1.0)
chk("the ten spans cover the recording",
    sum(o6[f"dec_min_{i+1:02d}"] for i in range(10)), 10.0)

print("\n7. per-stage means, counts and minutes on a constructed hypnogram")
fs = 10.0
# 600 samples (60 s) of each code in the order N3, N2, N1, REM, Wake, at distinct saturations
codes = [(1, 97.0), (2, 96.0), (3, 95.0), (4, 88.0), (5, 98.0)]
stage7 = np.concatenate([np.full(600, c, np.int8) for c, _ in codes])
x7 = np.concatenate([np.full(600, v) for _, v in codes])
valid7 = np.isfinite(x7) & (x7 >= 50.0) & (x7 <= 100.0)
o7 = {}
stage_profile(x7, valid7, stage7, fs, o7)
for nm, (_, v) in zip(("n3", "n2", "n1", "rem", "wake"), codes):
    chk(f"{nm} mean saturation", o7[f"{nm}_mean"], v)
    chk(f"{nm} minutes", o7[f"{nm}_min"], 1.0)
    chk(f"{nm} valid count", o7[f"{nm}_n"], 600)
chk("REM percent below 90 is total", o7["rem_t90"], 100.0)
chk("N2 percent below 90 is zero", o7["n2_t90"], 0.0)
chk("total sleep time is the four sleep stages", o7["tst_min"], 4.0)
chk("staged minutes are all five", o7["staged_min"], 5.0)
chk("the five stage counts sum to the staged valid count",
    sum(o7[f"{nm}_n"] for nm in ("n3", "n2", "n1", "rem", "wake")), o7["staged_valid_n"])
chk("sleep mean is the mean over the four sleep stages",
    o7["sleep_mean"], np.mean([97.0, 96.0, 95.0, 88.0]))

print("   and the keep band applies inside a stage too")
x7b = x7.copy()
x7b[600:900] = 120.0                      # half of N2 goes out of band
valid7b = np.isfinite(x7b) & (x7b >= 50.0) & (x7b <= 100.0)
o7b = {}
stage_profile(x7b, valid7b, stage7, fs, o7b)
chk("N2 mean over the surviving half", o7b["n2_mean"], 96.0)
chk("N2 valid count halves", o7b["n2_n"], 300)
chk("N2 minutes are unchanged (scored time, not valid time)", o7b["n2_min"], 1.0)
chk("N2 valid minutes halve", o7b["n2_valid_min"], 0.5)
chk("the count identity survives",
    sum(o7b[f"{nm}_n"] for nm in ("n3", "n2", "n1", "rem", "wake")), o7b["staged_valid_n"])

print("\n8. stage_on_grid places events on the SpO2 grid at the right samples")
starts = np.array([0.0, 30.0, 60.0, 90.0], np.float64)
durs = np.array([30.0, 30.0, 30.0, 30.0], np.float64)
cc = np.array([5, 3, 2, 4], np.int32)          # Wake, N1, N2, REM
g = {"starts": starts, "durations": durs, "codes": cc}
out8 = {}
st8 = stage_on_grid(g, 120.0, fs, 1200, out8)
chk_true("an array came back", st8 is not None)
for code, lo, hi in ((5, 0, 300), (3, 300, 600), (2, 600, 900), (4, 900, 1200)):
    chk_true(f"samples {lo}-{hi} all carry code {code}", bool((st8[lo:hi] == code).all()))
chk("four rows placed", out8["n_stage_rows_placed"], 4)
chk("no row out of range", out8["stage_out_of_range_frac"], 0.0)

print("\n9. the whole-day offset repair reproduces case 8 exactly")
g9 = {"starts": starts + 86400.0, "durations": durs, "codes": cc}
out9 = {}
st9 = stage_on_grid(g9, 120.0, fs, 1200, out9)
chk("offset detected as one day", out9["stage_day_offset_days"], 1)
chk_true("the repaired stage array is identical to case 8",
         st9 is not None and bool((st9 == st8).all()))

print("\n10. the digital-to-physical affine, with a known linear map")
# digital 0..1000 maps to physical 50..100, so digital d -> 50 + d/20
dig = np.arange(0, 1000, dtype=np.float64)
xa, meta = affine(dig, {"dig_min": 0.0, "dig_max": 1000.0,
                        "phys_min": 50.0, "phys_max": 100.0})
chk("first sample maps to 50.00", xa[0], 50.0)
chk("last sample maps to 99.95", xa[-1], 99.95)
o10 = run_deciles(xa)
for i in range(10):
    want = 50.0 + (np.arange(i * 100, (i + 1) * 100) / 20.0).mean()
    chk(f"decile {i+1} mean under the affine", o10[f"dec_mean_{i+1:02d}"], want)
identities(o10, "affine")

print("\n" + "=" * 84)
print(f"VERDICT: {'ALL KNOWN-ANSWER TESTS PASS' if not F else 'FAILURES: ' + ', '.join(F)}")
print("=" * 84)
raise SystemExit(1 if F else 0)
