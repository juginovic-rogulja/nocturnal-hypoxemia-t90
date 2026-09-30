"""
Known-answer test of the T80/T85/T88/T90/T92/T95 threshold arithmetic in
worker_ladder.py, plus the A95 (percent at or above 95) complement identity.

The fleet validation proves the RULE reproduces two frozen columns on real
recordings. This proves the ARITHMETIC on waveforms whose answer is known by
construction, which the frozen columns cannot do for the new thresholds because
no frozen T80/T92/T95 exists. Extends t85_analysis/selftest_t85.py. Run with:
    python3 selftest_ladder.py
"""
import numpy as np

SPO2_LO, SPO2_HI, TARGET_FS = 50.0, 100.0, 200.0
THRESHOLDS = (80.0, 85.0, 88.0, 90.0, 92.0, 95.0)
KEYS = ("t80", "t85", "t88", "t90", "t92", "t95")


def prov(raw, fs, dig=None, phys=None):
    """The worker's prov_ block, transcribed, returning all six thresholds + a95."""
    raw = np.asarray(raw, np.float64).ravel()
    x = raw
    if dig is not None:
        dmin, dmax = dig
        pmin, pmax = phys
        r = np.clip(raw, min(dmin, dmax), max(dmin, dmax))
        x = pmin + (r - dmin) * (pmax - pmin) / (dmax - dmin)
    n_out = int(round(x.size * TARGET_FS / fs))
    if not (x.size >= 32 and 32 <= n_out <= 60_000_000):
        return None
    if abs(TARGET_FS / fs - 1.0) < 1e-12:
        y = x
    else:
        from scipy.signal import resample
        y = resample(x, n_out)
    y = np.asarray(y, np.float32).astype(np.float64)
    yv = y[np.isfinite(y) & (y >= SPO2_LO) & (y <= SPO2_HI)]
    if yv.size < 10:
        return None
    out = {k: float((yv < t).mean() * 100) for k, t in zip(KEYS, THRESHOLDS)}
    out["a95"] = float((yv >= 95.0).mean() * 100)     # independent, as the worker emits it
    out["n"] = int(yv.size)
    return out


F = []


def chk(name, got, want, tol=1e-9):
    ok = abs(got - want) <= tol
    print(f"  {'PASS' if ok else 'FAIL'}  {name:<58} got {got:12.6f}  want {want:12.6f}")
    if not ok:
        F.append(name)


print("1. flat blocks at 200 Hz, no resample, exact fractions by construction")
# 1000 samples in blocks chosen so every threshold slices a different fraction:
#   100 at 78  (below all six)
#   100 at 82  (below 85/88/90/92/95, not 80)
#   150 at 86  (below 88/90/92/95)
#   150 at 89  (below 90/92/95)
#   100 at 91  (below 92/95)
#   200 at 94  (below 95 only)
#   200 at 97  (below none)
v = np.concatenate([np.full(100, 78.0), np.full(100, 82.0), np.full(150, 86.0),
                    np.full(150, 89.0), np.full(100, 91.0), np.full(200, 94.0),
                    np.full(200, 97.0)])
r = prov(v, 200.0)
chk("t80 = 100/1000", r["t80"], 10.0)
chk("t85 = 200/1000", r["t85"], 20.0)
chk("t88 = 350/1000", r["t88"], 35.0)
chk("t90 = 500/1000", r["t90"], 50.0)
chk("t92 = 600/1000", r["t92"], 60.0)
chk("t95 = 800/1000", r["t95"], 80.0)
chk("a95 = 200/1000", r["a95"], 20.0)
chk("kept sample count", r["n"], 1000)

print("\n2. the keep band [50,100] drops out-of-range samples from the DENOMINATOR")
v2 = np.concatenate([v, np.full(500, 30.0), np.full(200, 120.0)])
r2 = prov(v2, 200.0)
chk("denominator still 1000 after dropping 700", r2["n"], 1000)
chk("t80 unchanged at 10", r2["t80"], 10.0)
chk("t92 unchanged at 60", r2["t92"], 60.0)
chk("t95 unchanged at 80", r2["t95"], 80.0)
chk("a95 unchanged at 20", r2["a95"], 20.0)

print("\n3. strictly-below, not at-or-below, at every threshold")
for thr, key in zip(THRESHOLDS, KEYS):
    ex = np.concatenate([np.full(400, thr), np.full(600, 99.0)])   # 400 sit exactly ON
    chk(f"exactly {thr:.0f} does not count toward {key}", prov(ex, 200.0)[key], 0.0)
    be = np.concatenate([np.full(400, thr - 1.0), np.full(600, 99.0)])
    chk(f"one point below {thr:.0f} does count toward {key}", prov(be, 200.0)[key], 40.0)
print("  and a95 is AT-or-above, so a sample exactly at 95 counts toward a95")
ex95 = np.concatenate([np.full(400, 95.0), np.full(600, 80.0)])
chk("exactly 95 counts toward a95", prov(ex95, 200.0)["a95"], 40.0)

print("\n4. monotone by construction on random signals, t80<=t85<=t88<=t90<=t92<=t95")
rng = np.random.RandomState(20260821)
worst = 0.0
worst_id = 0.0
for i in range(200):
    s = rng.uniform(60, 100, size=rng.randint(500, 5000))
    q = prov(s, 200.0)
    for a, b in zip(KEYS[:-1], KEYS[1:]):
        worst = max(worst, q[a] - q[b])
    worst_id = max(worst_id, abs(q["a95"] + q["t95"] - 100.0))
chk("max ordering violation over 200 random signals", worst, 0.0, tol=1e-12)
chk("max |a95 + t95 - 100| over 200 random signals", worst_id, 0.0, tol=1e-9)

print("\n5. the digital-to-physical affine, with a known linear map")
# digital 0..1000 maps to physical 50..100, so digital d -> 50 + d/20
dig = np.arange(0, 1000, dtype=np.float64)          # physical 50.00 .. 99.95
r5 = prov(dig, 200.0, dig=(0.0, 1000.0), phys=(50.0, 100.0))
# below thr means 50 + d/20 < thr, so d < 20*(thr-50)
chk("t80 with affine = 600/1000", r5["t80"], 60.0)
chk("t85 with affine = 700/1000", r5["t85"], 70.0)
chk("t88 with affine = 760/1000", r5["t88"], 76.0)
chk("t90 with affine = 800/1000", r5["t90"], 80.0)
chk("t92 with affine = 840/1000", r5["t92"], 84.0)
chk("t95 with affine = 900/1000", r5["t95"], 90.0)
chk("a95 with affine = 100/1000", r5["a95"], 10.0)

print("\n6. resampling from 1 Hz to 200 Hz preserves long flat fractions")
# 3000 s at 1 Hz: 300 s at 78, 300 s at 90.5 (below 92/95 only), 2400 s at 97.
v6 = np.concatenate([np.full(300, 78.0), np.full(300, 90.5), np.full(2400, 97.0)])
r6 = prov(v6, 1.0)
for key, want in (("t80", 10.0), ("t92", 20.0), ("t95", 20.0), ("a95", 80.0)):
    ok = abs(r6[key] - want) < 0.5
    print(f"  {'PASS' if ok else 'FAIL'}  resampled {key} within 0.5 pp of {want:.0f}       "
          f"     got {r6[key]:12.6f}  n {r6['n']:,}")
    if not ok:
        F.append(f"resample {key}")

print("\n7. the a95 complement identity on every constructed case above")
for nm, rr in (("case1", r), ("case2", r2), ("affine", r5), ("resample", r6)):
    chk(f"a95 + t95 = 100 exactly ({nm})", rr["a95"] + rr["t95"], 100.0, tol=1e-9)

print("\n" + "=" * 78)
print(f"VERDICT: {'ALL KNOWN-ANSWER TESTS PASS' if not F else 'FAILURES: ' + ', '.join(F)}")
print("=" * 78)
raise SystemExit(1 if F else 0)
