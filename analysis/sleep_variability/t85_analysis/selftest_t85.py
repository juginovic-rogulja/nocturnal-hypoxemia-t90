"""
Known-answer test of the T85/T88/T90 threshold arithmetic in worker_t85.py.

The fleet validation proves the RULE reproduces two frozen columns on real recordings.
This proves the ARITHMETIC on waveforms whose answer is known by construction, which the
frozen columns cannot do for T85 because no frozen T85 exists. Run it with:
    python3 selftest_t85.py
"""
import numpy as np

SPO2_LO, SPO2_HI, TARGET_FS = 50.0, 100.0, 200.0


def prov(raw, fs, dig=None, phys=None):
    """The worker's prov_ block, transcribed, returning t90/t88/t85."""
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
    return {"t90": float((yv < 90).mean() * 100), "t88": float((yv < 88).mean() * 100),
            "t85": float((yv < 85).mean() * 100), "n": int(yv.size)}


F = []


def chk(name, got, want, tol=1e-9):
    ok = abs(got - want) <= tol
    print(f"  {'PASS' if ok else 'FAIL'}  {name:<52} got {got:12.6f}  want {want:12.6f}")
    if not ok:
        F.append(name)


print("1. flat blocks at 200 Hz, no resample, exact fractions by construction")
# 1000 samples: 100 at 80 (below all three), 150 at 86 (below 88 and 90, not 85),
# 250 at 89 (below 90 only), 500 at 97 (below none)
v = np.concatenate([np.full(100, 80.0), np.full(150, 86.0),
                    np.full(250, 89.0), np.full(500, 97.0)])
r = prov(v, 200.0)
chk("t85 = 100/1000", r["t85"], 10.0)
chk("t88 = 250/1000", r["t88"], 25.0)
chk("t90 = 500/1000", r["t90"], 50.0)
chk("kept sample count", r["n"], 1000)

print("\n2. the keep band [50,100] drops out-of-range samples from the DENOMINATOR")
# add 500 samples at 30 (below 50, dropped) and 200 at 120 (above 100, dropped)
v2 = np.concatenate([v, np.full(500, 30.0), np.full(200, 120.0)])
r2 = prov(v2, 200.0)
chk("denominator still 1000 after dropping 700", r2["n"], 1000)
chk("t85 unchanged at 10", r2["t85"], 10.0)
chk("t90 unchanged at 50", r2["t90"], 50.0)

print("\n3. strictly-below, not at-or-below, at every threshold")
for thr, key in ((85.0, "t85"), (88.0, "t88"), (90.0, "t90")):
    ex = np.concatenate([np.full(400, thr), np.full(600, 99.0)])   # 400 sit exactly ON
    chk(f"exactly {thr:.0f} does not count toward {key}", prov(ex, 200.0)[key], 0.0)
    be = np.concatenate([np.full(400, thr - 1.0), np.full(600, 99.0)])
    chk(f"one point below {thr:.0f} does count toward {key}", prov(be, 200.0)[key], 40.0)

print("\n4. monotone by construction on random signals, t85 <= t88 <= t90")
rng = np.random.RandomState(20260821)
worst = 0.0
for i in range(200):
    s = rng.uniform(60, 100, size=rng.randint(500, 5000))
    q = prov(s, 200.0)
    worst = max(worst, q["t85"] - q["t88"], q["t88"] - q["t90"])
chk("max ordering violation over 200 random signals", worst, 0.0, tol=1e-12)

print("\n5. the digital-to-physical affine, with a known linear map")
# digital 0..1000 maps to physical 50..100, so digital d -> 50 + d/20
dig = np.arange(0, 1000, dtype=np.float64)          # physical 50.00 .. 99.95
r5 = prov(dig, 200.0, dig=(0.0, 1000.0), phys=(50.0, 100.0))
# below 85 means 50 + d/20 < 85, so d < 700 -> 700 of 1000
chk("t85 with affine = 700/1000", r5["t85"], 70.0)
chk("t88 with affine = 760/1000", r5["t88"], 76.0)
chk("t90 with affine = 800/1000", r5["t90"], 80.0)

print("\n6. resampling from 1 Hz to 200 Hz preserves a long flat fraction")
# 3000 s at 1 Hz: 600 s at 80, 2400 s at 96. Band-limited resample keeps the plateaus,
# so the fraction below 85 stays close to 20 percent.
v6 = np.concatenate([np.full(600, 80.0), np.full(2400, 96.0)])
r6 = prov(v6, 1.0)
ok = abs(r6["t85"] - 20.0) < 0.5
print(f"  {'PASS' if ok else 'FAIL'}  resampled t85 within 0.5 pp of 20            "
      f"got {r6['t85']:12.6f}  n {r6['n']:,}")
if not ok:
    F.append("resample t85")

print("\n" + "=" * 78)
print(f"VERDICT: {'ALL KNOWN-ANSWER TESTS PASS' if not F else 'FAILURES: ' + ', '.join(F)}")
print("=" * 78)
raise SystemExit(1 if F else 0)
