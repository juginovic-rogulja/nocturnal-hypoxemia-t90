"""Owner's look at the event-locked ensemble SpO2 per site (step 021 checkpoint): work/hb_ensemble_curve.png."""
from __future__ import annotations
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import x1_paths as P
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

cal = json.load(open(P.HB_CALIBRATION)); d = json.load(open(P.HB_DEFINITION)) if os.path.isfile(P.HB_DEFINITION) else None
fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
for ax, site in zip(axes, ("I0002", "I0006")):
    e = cal["per_site"][site]; offs = np.array(e["offsets_sec"]); m = np.array([np.nan if v is None else v for v in e["mean_spo2"]])
    ax.plot(offs, m, color="black", lw=1.5)
    ax.axvline(0, color="grey", lw=0.8, ls="--"); ax.axvline(e["nadir_offset_sec"], color="red", lw=0.8, ls=":")
    if d and site in d["per_site"]:
        w = d["per_site"][site]; ax.axvspan(-w["w_pre"], w["w_post"], color="tab:blue", alpha=0.12, label=f"frozen window [-{w['w_pre']:.0f}, +{w['w_post']:.0f}] s")
    ax.axvspan(-100, 100, color="tab:orange", alpha=0.06, label="alternative fixed 100/100 s")
    ax.set_title(f"{site}: {e['n_nights']} nights, {e['n_events']:,} events\npre-max {e['pre_max_offset_sec']:+.0f} s, nadir {e['nadir_offset_sec']:+.0f} s ({e['spo2_at_pre_max']-e['spo2_at_nadir']:.2f} pt), return {e['post_return_offset_sec']:+.0f} s", fontsize=9)
    ax.set_xlabel("seconds from respiratory event END"); ax.legend(fontsize=7, loc="lower right"); ax.grid(alpha=0.3)
axes[0].set_ylabel("ensemble mean SpO2 (%)")
fig.suptitle("Event-locked ensemble SpO2, 300-night pilot (step 020/021, decision D5)", fontsize=10)
fig.tight_layout(); out = f"{P.X1_WORK}/hb_ensemble_curve.png"; fig.savefig(out, dpi=150); print("wrote", out)
