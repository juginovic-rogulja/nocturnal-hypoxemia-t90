"""Step 011 synthetic tests (unittest-compatible, also runnable under pytest).
(a) t_end=None reproduces extract_night.extract on a synthetic h5, key by key (exact);
    the cut reader's rasters equal the full rasters truncated; the cut changes only what it should.
(b) a 5-point 60-s square desaturation gives 5 percent-minutes per event; two events do not double count.
(c) PLMI series rule: a 4-movement train -> 4, a 3-movement train -> 0; bilateral merge; resp exclusion.
(d) SOL and REM latency on a synthetic hypnogram (pure function and the pipeline agree).
(e) sleep-period T90 on a synthetic trace.
(f) clean_event_times painted onto the grid equals the reader's rasterise (day-wrapped rows included)."""
from __future__ import annotations
import json, math, os, sys, tempfile, unittest
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)
import new_measures as NM  # noqa: E402
from prepap_reader import PrePapReader, clean_event_times  # noqa: E402
import extract_v8  # noqa: E402
from synth import make_synthetic_h5  # noqa: E402
sys.path.insert(0, extract_v8.REEX)
import extract_night  # noqa: E402


def _eq(a, b):
    if isinstance(a, float) and isinstance(b, float) and math.isnan(a) and math.isnan(b):
        return True
    if isinstance(a, (float, np.floating)) and isinstance(b, (float, np.floating)):
        return float(a) == float(b) or (math.isnan(a) and math.isnan(b))
    return a == b


class TestNoOp(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        resp = [(900.0, 920.0, 1), (1500.0, 1525.0, 9), (2400.0, 2430.0, 2), (3300.0, 3320.0, 9)]
        limb = [(700 + 20 * k, 700 + 20 * k + 2.0, 1) for k in range(6)]
        cls.p = make_synthetic_h5(os.path.join(cls.tmp, "syn.h5"), dur_sec=3600, resp=resp, limb=limb)

    def test_a_noop_reproduces_extract_night(self):
        ref = extract_night.extract(self.p)
        new = extract_v8.extract(self.p, t_end_sec=None, families=("full",))
        self.assertEqual(ref["status"], "ok"); self.assertEqual(new["status"], "ok")
        skip = {"source", "runtime_sec"}
        diffs = [k for k in ref if k not in skip and not _eq(ref[k], new.get(k, "MISSING"))]
        self.assertEqual(diffs, [], f"keys differing from extract_night: {diffs[:20]}")
        self.assertFalse(new["_v8_prepap_window"])
        self.assertAlmostEqual(new["sol_min_nm"], new["sol_min"]); self.assertAlmostEqual(new["rem_latency_min_nm"], new["rem_latency_min"])

    def test_a2_cut_rasters_equal_truncated_full(self):
        t_end = 2000.0
        with PrePapReader(self.p) as full, PrePapReader(self.p, t_end_sec=t_end) as cut:
            n_cut = cut.n_samples
            self.assertEqual(n_cut, int(round(t_end * full.fs)))
            for kind in ("stage", "arousal", "resp_3", "resp_4", "limb"):
                a = full.rasterise(kind)[:n_cut]; b = cut.rasterise(kind)
                self.assertTrue(np.array_equal(a, b), kind)
            self.assertTrue(np.array_equal(full.read_raw("spo2")[:n_cut], cut.read_raw("spo2")))
            self.assertTrue(np.array_equal(full.read_raw("f3-m2")[:n_cut], cut.read_raw("f3-m2")))
            self.assertEqual(cut.read_native("spo2").size, int(round(t_end * 25)))
            self.assertTrue(np.array_equal(full.read_native("spo2")[: int(round(t_end * 25))], cut.read_native("spo2")))
            # event table: rows starting at or after the cut are dropped, ends clipped
            c, s, e = cut.clean_events("resp_3")
            self.assertEqual(s.size, 2); self.assertTrue((s < t_end).all() and (e <= t_end).all())
        full_rec = extract_v8.extract(self.p, None, ("full",)); cut_rec = extract_v8.extract(self.p, t_end, ("full",))
        self.assertAlmostEqual(cut_rec["recording_dur_min"], t_end / 60.0)
        self.assertAlmostEqual(cut_rec["_v8_window_end_min"], t_end / 60.0)
        self.assertAlmostEqual(cut_rec["_v8_recording_dur_min_full"], full_rec["recording_dur_min"])
        self.assertTrue(cut_rec["_v8_prepap_window"])
        self.assertEqual(cut_rec["resp_3_rows"], 2); self.assertEqual(full_rec["resp_3_rows"], 4)
        self.assertLess(cut_rec["TST_min"], full_rec["TST_min"])
        self.assertEqual(cut_rec["sol_min"], full_rec["sol_min"])
        # light family on the same file agrees with the full family on the shared keys
        light = extract_v8.extract(self.p, t_end, ("light",))
        for k in ("TST_min", "spo2_pct_below_90", "spo2_pct_below_90_native", "spo2_pct_below_90_sleep", "hypoxic_burden",
                  "plm_index", "sol_min", "rem_latency_min", "AHI", "resp_3_rows"):
            self.assertTrue(_eq(light[k], cut_rec[k]), k)

    def test_f_clean_events_match_rasterise(self):
        with PrePapReader(self.p) as r:
            n = r.n_samples; fs = r.fs; D = n / fs
            for kind in ("stage", "arousal", "resp_3", "limb"):
                codes, starts, ends, durs, _ = r._event_table(kind)
                c, s, e = clean_event_times(codes, starts, ends, durs, D, kind)
                out = np.zeros(n, np.int16)
                for cc, ss, ee in zip(c, s, e):
                    a = max(int(round(ss * fs)), 0); b = min(int(round(ee * fs)), n)
                    if b > a:
                        out[a:b] = cc
                self.assertTrue(np.array_equal(out, r.rasterise(kind)), kind)
            self.assertEqual(r.diag["arousal_n_wrapped"], 3)


class TestNewMeasures(unittest.TestCase):
    def test_b_hypoxic_burden_square(self):
        fs = 25.0; T = 3000
        x = np.full(int(T * fs), 95.0)
        # event ends at 1000 s; desaturation 5 points from 990 to 1050 s (60 s)
        x[int(990 * fs):int(1050 * fs)] = 90.0
        hb = NM.hypoxic_burden(x, fs, [1000.0], 100.0, 100.0, tst_hours=1.0)
        self.assertAlmostEqual(hb["hb_area_pctmin"], 5.0, places=6); self.assertAlmostEqual(hb["hypoxic_burden"], 5.0, places=6)
        # two events 60 s apart with the same square each: 10 in total, no double counting
        x2 = np.full(int(T * fs), 95.0)
        x2[int(990 * fs):int(1050 * fs)] = 90.0; x2[int(1490 * fs):int(1550 * fs)] = 90.0
        hb2 = NM.hypoxic_burden(x2, fs, [1000.0, 1500.0], 100.0, 100.0, tst_hours=2.0)
        self.assertAlmostEqual(hb2["hb_area_pctmin"], 10.0, places=6); self.assertAlmostEqual(hb2["hypoxic_burden"], 5.0, places=6)
        # overlapping windows: events 30 s apart sharing one 60-s square -> the square is counted once
        x3 = np.full(int(T * fs), 95.0); x3[int(990 * fs):int(1050 * fs)] = 90.0
        hb3 = NM.hypoxic_burden(x3, fs, [1000.0, 1030.0], 100.0, 100.0, tst_hours=1.0)
        self.assertAlmostEqual(hb3["hb_area_pctmin"], 5.0, places=6)
        # invalid samples (0 = dropout) are not counted
        x4 = x.copy(); x4[int(1000 * fs):int(1010 * fs)] = 0.0
        hb4 = NM.hypoxic_burden(x4, fs, [1000.0], 100.0, 100.0, tst_hours=1.0)
        self.assertAlmostEqual(hb4["hb_area_pctmin"], 5.0 - 5.0 * 10 / 60, places=6)
        self.assertTrue(math.isnan(NM.hypoxic_burden(x, fs, [1000.0], 100, 100, tst_hours=0.0)["hypoxic_burden"]))
        self.assertEqual(NM.hypoxic_burden(x, fs, [], 100, 100, 1.0)["hb_area_pctmin"], 0.0)

    def test_c_plm_series_rule(self):
        on4 = [100, 120, 140, 160]; d = [2.0] * 4
        r = NM.plm_index(on4, d, tst_hours=1.0)
        self.assertEqual(r["plm_count"], 4); self.assertAlmostEqual(r["plm_index"], 4.0); self.assertEqual(r["plm_n_series"], 1)
        r3 = NM.plm_index(on4[:3], d[:3], tst_hours=1.0)
        self.assertEqual(r3["plm_count"], 0); self.assertAlmostEqual(r3["plm_index"], 0.0); self.assertAlmostEqual(r3["lm_index"], 3.0)
        # a 100-s gap breaks the series: 4 + 3 -> 4 counted
        r7 = NM.plm_index([100, 120, 140, 160, 300, 320, 340], [2.0] * 7, 1.0)
        self.assertEqual(r7["plm_count"], 4)
        # durations outside 0.5-10 s are not limb movements
        rd = NM.plm_index(on4, [2.0, 2.0, 12.0, 2.0], 1.0)
        self.assertEqual(rd["plm_count"], 0); self.assertEqual(rd["lm_count"], 3)
        # bilateral: left and right within 5 s are one movement
        rb = NM.plm_index([100, 101, 120, 121, 140, 141, 160, 161], [1.0] * 8, 1.0)
        self.assertEqual(rb["lm_count"], 4); self.assertEqual(rb["plm_count"], 4)
        # respiratory exclusion: a movement 0.3 s after an event end is excluded in the variant only
        rr = NM.plm_index(on4, d, 1.0, resp_intervals=[(90.0, 99.8)])
        self.assertEqual(rr["plm_count"], 4); self.assertEqual(rr["plm_count_excl_resp"], 0); self.assertEqual(rr["lm_count_excl_resp"], 3)
        self.assertTrue(math.isnan(NM.plm_index(on4, d, 0.0)["plm_index"]))

    def test_d_latencies(self):
        fs = 200
        ep = [5] * 6 + [3, 2, 2, 1] + [4] * 2 + [5]
        st = np.repeat(np.array(ep, np.float32), 30 * fs)
        lat = NM.latencies(st, fs)
        self.assertAlmostEqual(lat["sol_min"], 3.0); self.assertAlmostEqual(lat["rem_latency_min"], 2.0)
        self.assertTrue(math.isnan(NM.latencies(np.repeat(np.float32(5), 30 * fs), fs)["sol_min"]))
        self.assertTrue(math.isnan(NM.latencies(np.repeat(np.array([5, 2, 2], np.float32), 30 * fs), fs)["rem_latency_min"]))

    def test_e_sleep_period_t90(self):
        fs200 = 200; nfs = 25
        ep = [5] * 2 + [2] * 2 + [5] * 2     # 3 min: wake, sleep, wake
        st = np.repeat(np.array(ep, np.float32), 30 * fs200)
        x = np.full(180 * nfs, 95.0, np.float32)
        x[: 60 * nfs] = 85.0                  # low during the first wake minute only
        x[60 * nfs: 60 * nfs + 15 * nfs] = 85.0   # low during the first 15 s of sleep
        r = NM.sleep_period_t90(x, nfs, st, fs200)
        self.assertAlmostEqual(r["spo2_pct_below_90_sleep"], 25.0); self.assertAlmostEqual(r["spo2_sleep_valid_min"], 1.0)
        whole = float(np.mean(x < 90) * 100)
        self.assertGreater(whole, r["spo2_pct_below_90_sleep"])
        x[70 * nfs: 80 * nfs] = 0.0           # dropout inside sleep is excluded from the denominator
        r2 = NM.sleep_period_t90(x, nfs, st, fs200)
        self.assertAlmostEqual(r2["spo2_sleep_valid_min"], 50 / 60)


if __name__ == "__main__":
    unittest.main(verbosity=2)
