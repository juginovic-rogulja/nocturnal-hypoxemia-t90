#!/usr/bin/env python3
"""Round 44b (2026-09-21, Alen's two further comments), on the round-40 engine (build_r40 imported):

  DIPS  every oxygen trace on every sketch has the SAME trough depth: one trace height TH_STD = 19.5 mm (Figure 6's) on the
        double-width sheets (Figure 1 f, 3 d, 4 e, 5 d, 6), and 10.3 mm on the print-width Supp 16 band (whose sheet is printed at
        1.0 where the double-width sheets print at 0.53, so the printed depth is the same 7.0 mm everywhere); inside a trace every
        dip reaches exactly the same trough (the baseline wobble is cancelled at each dip centre: trace_eq). Widths, dip counts,
        labels, arrows, organ rows, band pages and balance shifts unchanged.
  1a    Figure 1 panel a rebuilt from its round-36 builder (ROUND36 LF1G build_r36_fig1a.band, ported here) with: the leg-movement
        item folded into the sleep-structure item (one label "sleep structure and leg movements (EEG, staging, leg EMG)", no
        separate EMG trace or leader), the organ grid at scale 1.0 (= 0.8 x the round-34 1.25, Alen's 20% smaller), and the
        sheet's wording (141 measurements, 54 future diseases). Band page 968.66 x 231.114 pt as delivered.

    python3 build_r44b.py --final [--only stem,stem]
"""
from __future__ import annotations
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import argparse, fcntl, json, math, sys, time
from pathlib import Path
R40 = Path(f"{paths.FIGURE_ROOT}/ROUND40_2026-09-18/lanes/LS_SKETCH/scripts")
R44 = Path(f"{paths.FIGURE_ROOT}/ROUND44_2026-09-21/lanes/LSKETCH")
sys.path.insert(0, str(R40)); sys.path.insert(0, str(R44 / "scripts"))
import build_r40 as B                                             # noqa: E402
import pandas as pd                                               # noqa: E402
OUT = R44 / "work" / "build"; WORK = R44 / "work"
B.L.PREVIEWS = WORK / "previews"
new, finish, arrow, organ_row_f = B.new, B.finish, B.arrow, B.organ_row_f
AZ, CARD, GREEN_D, GREY, GREY_L, INK21, AMBER, MM, INK = B.AZ, B.CARD, B.GREEN_D, B.GREY, B.GREY_L, B.INK21, B.AMBER, B.MM_PER_PT, B.INK
TH_STD = 19.5                                                     # the common trace height (Figure 6's), dip trough at 0.92 TH below the trace top
TH_S16 = 10.3                                                     # Supp 16 (print-width sheet): the same printed depth
NUM = Path(paths.NUMBERS_DIR)
FAMILY_COL = {"Oxygenation": "#0288d1", "Heart rate and variability": "#8d2dfa", "Breathing events": "#d55e00", "Brain, microstructure": "#fa93a5",
              "Brain, spectral power": "#f2e279", "Sleep timing and structure": "#39c445", "Limb movements": "#8a9099"}
RED_D, RED_L = "#b3261e", "#fbe4e1"


# ================================================================ the equal-trough trace (replaces build_r40.trace for every band built here)
def trace_eq(fig, x, y, w, h, *, color, dips, ref_at=0.5, lw=0.6, label=B.LABEL90, label_pt=None, label_bold=False,
             label_gap=B.LABEL90_GAP, ref_lw=0.3, ref_dash="1.4 1.0", label_color=INK21, ref_color=INK21, label_side="left"):
    """build_r40.trace with one change: each dip's Gaussian amplitude is h*DEEP_DEPTH minus the baseline wobble at the dip centre,
    so every trough sits exactly at y + h*(0.24 + DEEP_DEPTH) = y + 0.92 h. Same return record."""
    y90 = y + h * ref_at
    n = 220
    wob = lambda t: h * 0.03 * math.sin(t * 2 * math.pi * 7)
    amp = {c: h * B.DEEP_DEPTH - wob(c) for c in dips}
    pts = []
    for i in range(n + 1):
        t = i / n
        py = y + h * 0.24 + wob(t)
        for c in dips:
            py += amp[c] * math.exp(-((t - c) / B.DEEP_SIGMA) ** 2)
        pts.append((x + w * t, py))
    fig.line(x, y90, x + w, y90, stroke=ref_color, lw=ref_lw, dash=ref_dash if B.THRESHOLD_DASH else None)
    if label_side == "left":
        lb = fig.label(label, x - label_gap, y90, pt=label_pt or fig.theme.note_pt, color=label_color, anchor="end", va="middle", bold=label_bold)
    else:
        lb = fig.label(label, x + w + label_gap, y90, pt=label_pt or fig.theme.note_pt, color=label_color, anchor="start", va="middle", bold=label_bold)
    seg, lenses = [], 0
    for px, py in pts + [(x + w, y90 - 1)]:
        if py > y90:
            seg.append((px, py))
        elif seg:
            fig.polygon([(seg[0][0], y90)] + seg + [(seg[-1][0], y90)], fill=color, opacity=0.35)
            seg = []; lenses += 1
    fig.path("M" + " L".join(f"{px:.2f} {py:.2f}" for px, py in pts), stroke=color, lw=lw)
    troughs = [round(y + h * 0.24 + wob(c) + amp[c], 4) for c in dips]
    assert max(troughs) - min(troughs) < 1e-6, troughs
    return {"y90": round(y90, 3), "n_dips": len(dips), "lenses": lenses, "color": color, "trough_y": troughs[0], "trough_depth_below_90_mm": round(troughs[0] - y90, 3),
            "trough_below_90_mm": round(max(py - y90 for _, py in pts), 3), "trace_h_mm": h, "equal_troughs": True,
            "ink_top": round(min(py for _, py in pts), 3), "ink_bottom": round(max(py for _, py in pts), 3),
            "label_box": [round(lb.x, 3), round(lb.y, 3), round(lb.x + lb.w, 3), round(lb.y + lb.h, 3)],
            "label_pt": label_pt or fig.theme.note_pt, "label_color": label_color, "ref_color": ref_color, "ref_lw": ref_lw,
            "ref_dash": ref_dash if B.THRESHOLD_DASH else None, "label_side": label_side, "trace_box": [round(x, 3), round(y, 3), round(x + w, 3), round(y + h, 3)]}


B.trace = trace_eq                                                # build_r40.trace_pair, band_4e, band_5d now draw equal troughs
trace = trace_eq


def box(p): return B.box(p)


# ================================================================ Figure 6 (as this morning, TH 19.5 already, troughs now equal)
W6 = B.W6_PT * MM; H6_MM = 100.0; PITCH_IN, BED_W = 43.0, 26.0; BED_H = BED_W * 82.0 / 122.0; ROW_DY = 24.0


def inputs_card(fig, cxi, mid, *, box_pad=8.0):
    probe = new(W6, H6_MM)
    lab_h = probe.label("sleep apnea", cxi, 0.0, anchor="middle", va="top", bold=True, color=AMBER).h
    c_top, c_bot = -PITCH_IN / 2 - BED_H / 2, PITCH_IN / 2 + 11.0 + lab_h
    y_mid_items = mid - (c_top + c_bot) / 2
    y_dur, y_ap = y_mid_items - PITCH_IN / 2, y_mid_items + PITCH_IN / 2
    k = BED_W / 20.0
    b1 = fig.place(B.B45.SLEEPER, cxi, y_dur, w=BED_W)
    l1 = fig.label("sleep duration", cxi, y_dur + 11.0, anchor="middle", va="top", bold=True, color=GREEN_D)
    B.airflow(fig, cxi, y_ap - 10.0, 16 * k, 5.5 * k, AMBER)
    af_box = [cxi - 8 * k, y_ap - 10.0 - 2.75 * k, cxi + 8 * k, y_ap - 10.0 + 2.75 * k]
    b3 = fig.place(B.B45.SLEEPER, cxi, y_ap, w=BED_W)
    l3 = fig.label("sleep apnea", cxi, y_ap + 11.0, anchor="middle", va="top", bold=True, color=AMBER)
    parts = [box(p) for p in (b1, l1, b3, l3)] + [af_box]
    cx0, cx1 = min(p[0] for p in parts), max(p[2] for p in parts); cy0, cy1 = min(p[1] for p in parts), max(p[3] for p in parts)
    half_w = max(cxi - cx0, cx1 - cxi) + box_pad
    card = [round(cxi - half_w, 3), round(cy0 - box_pad, 3), round(2 * half_w, 3), round(cy1 - cy0 + 2 * box_pad, 3)]
    assert abs((card[1] + card[3] / 2) - mid) < 0.02, (card, mid)
    return dict(card=card, y_dur=round(y_dur, 3), y_ap=round(y_ap, 3), items={"duration": box(b1), "apnea": box(b3)}, labels={"duration": box(l1), "apnea": box(l3)}, airflow=[round(v, 3) for v in af_box])


def organ_panel(fig, pcx, pcy, *, dark, pitch=25.0, scale=1.0, pad=None):
    pad = B.PANEL_PAD if pad is None else pad
    ext = B.row_extents(B.ORGANS_LIVE, pitch, scale, family=True)
    n0 = len(fig._placed)
    organ_row_f(fig, pcx - ext["ink_offset"], pcy, dark=dark, pitch=pitch, scale=scale, names=False)
    icons = fig._placed[n0:]
    assert len(icons) == 5 and not any("brain" in p.name for p in icons), [p.name for p in icons]
    x0, x1 = min(p.x for p in icons), max(p.x + p.w for p in icons)
    assert abs((x0 + x1) / 2 - pcx) < 0.01, ("ink centred", x0, x1, pcx)
    half = max(pcx - x0, x1 - pcx) + pad; ph = 2 * (B.ICON_HALF * scale + pad)
    return dict(rect=[round(pcx - half, 3), round(pcy - ph / 2, 3), round(2 * half, 3), round(ph, 3)], ink_x=[round(x0, 3), round(x1, 3)], names=[p.name for p in icons], pitch=pitch, scale=scale)


def fig6_core(fig, *, cxi, tx, organ_arrow, mid, tw=34.0, th=TH_STD, pitch=25.0, scale=1.0):
    n_parts0 = len(fig._parts)
    inp = inputs_card(fig, cxi, mid); rx1 = inp["card"][0] + inp["card"][2]
    y_top, y_bot = mid - ROW_DY, mid + ROW_DY
    tr_top = trace(fig, tx, y_top - th / 2, tw, th, color=AZ, dips=B.DIPS_PRESERVED, lw=1.0, label_pt=fig.theme.note_pt, label_bold=True, label_gap=0.7, ref_lw=0.5, ref_dash="1.4 1.2", label_side="left")
    tr_bot = trace(fig, tx, y_bot - th / 2, tw, th, color=AZ, dips=B.DIPS_LOW, lw=1.0, label_pt=fig.theme.note_pt, label_bold=True, label_gap=0.7, ref_lw=0.5, ref_dash="1.4 1.2", label_side="left")
    tips = [(tr_top["label_box"][0] - B.HEAD_CLEAR, y_top), (tr_bot["label_box"][0] - B.HEAD_CLEAR, y_bot)]
    fork = [arrow(fig, (rx1, mid), tp) for tp in tips]
    ax0 = tx + tw + B.ORGAN_ARROW_GAP; px0 = ax0 + organ_arrow + B.HEAD_CLEAR
    ext = B.row_extents(B.ORGANS_LIVE, pitch, scale, family=True); half = max(ext["left_of_centre"], ext["right_of_centre"]) + B.PANEL_PAD; pcx = px0 + half
    oa = [arrow(fig, (ax0, y_top), (px0 - B.HEAD_CLEAR, y_top)), arrow(fig, (ax0, y_bot), (px0 - B.HEAD_CLEAR, y_bot))]
    pan_top = organ_panel(fig, pcx, y_top, dark=False, pitch=pitch, scale=scale); pan_bot = organ_panel(fig, pcx, y_bot, dark=True, pitch=pitch, scale=scale)
    assert pan_top["ink_x"] == pan_bot["ink_x"]
    n1 = len(fig._parts)
    for r in (inp["card"], pan_top["rect"], pan_bot["rect"]): fig.rect(*r, rx=B.PANEL_RX, fill=CARD, stroke=GREY_L, lw=0.3)
    moved = fig._parts[n1:]; del fig._parts[n1:]; fig._parts[n_parts0:n_parts0] = moved
    for tp in tips: assert tp[0] - rx1 >= 25.0, ("fork run at least 25 mm", tp, rx1)
    right_edge = pan_top["rect"][0] + pan_top["rect"][2]
    assert right_edge < W6 - 5.0 and inp["card"][0] > 5.0, (inp["card"], right_edge)
    return dict(inputs=inp, traces={"top": tr_top, "bottom": tr_bot}, fork=fork, organ_arrows=oa, panels={"top": pan_top, "bottom": pan_bot}, rows_y=[y_top, y_bot], content_x=[inp["card"][0], round(right_edge, 3)], mid=mid)


def fig6_r44b(final, *, stem="Main_Fig6_r44b", outdir=OUT):
    fig = new(W6, H6_MM); geo = fig6_core(fig, cxi=61.0, tx=143.0, organ_arrow=20.0, mid=H6_MM / 2)
    rec = finish(fig, stem, outdir, final=final, attribution=None); rec.update(geometry=geo, page_pt=[round(W6 / MM, 2), round(H6_MM / MM, 2)], trace_h=TH_STD, moved="round 44b: equal troughs (trace_eq); layout as Main_Fig6_r44"); return rec


def lerp_hex(a, b, t):
    ca = [int(a[i:i + 2], 16) for i in (1, 3, 5)]; cb = [int(b[i:i + 2], 16) for i in (1, 3, 5)]
    return "#" + "".join(f"{round(x + (y - x) * t):02x}" for x, y in zip(ca, cb))


def triangles(fig, *, x_left, top, h, w=26.0, gap=10.0):
    rk = pd.read_csv(NUM / "ranking_v3.csv", comment="#"); fam = pd.read_csv(NUM / "measure_families.csv").set_index("feature")["family"]
    rk["family"] = rk.feature.map(fam); assert len(rk) == 141 and rk.family.notna().all()
    order = rk.sort_values("rank"); t = fig.theme; cx1 = x_left + w / 2; cx2 = cx1 + w + gap; n = len(order); step = h / n
    for i, (_, r) in enumerate(order.iterrows()):
        tt = (i + 0.5) / n; y = top + tt * h; hw = (w / 2) * (1 - tt)
        if hw > 0.15: fig.line(cx1 - hw, y, cx1 + hw, y, stroke=FAMILY_COL[r.family], lw=step * 0.92, cap="butt")
    fig.polygon([(cx1 - w / 2, top), (cx1 + w / 2, top), (cx1, top + h)], fill=None, stroke=GREY, lw=0.35)
    NS = 48
    for j in range(NS):
        t0, t1 = j / NS, (j + 1) / NS
        fig.polygon([(cx2 - (w / 2) * (1 - t0), top + t0 * h), (cx2 + (w / 2) * (1 - t0), top + t0 * h), (cx2 + (w / 2) * (1 - t1), top + t1 * h), (cx2 - (w / 2) * (1 - t1), top + t1 * h)], fill=lerp_hex(RED_D, RED_L, (t0 + t1) / 2), stroke=None)
    fig.polygon([(cx2 - w / 2, top), (cx2 + w / 2, top), (cx2, top + h)], fill=None, stroke=GREY, lw=0.35)
    t1l = fig.label("141 sleep measurements", cx1, top - 2.5, pt=t.note_pt, bold=True, color=INK21, anchor="middle", va="bottom", max_w=w + 14.0, leading=1.15)
    t2l = fig.label("disease risk", cx2, top - 2.5, pt=t.note_pt, bold=True, color=INK21, anchor="middle", va="bottom")
    r_tst = int(order[order.feature == "TST_min"]["rank"].iloc[0]); y_ox = top + (2.0 / n) * h; y_tst = top + ((r_tst - 0.5) / n) * h; hw_tst = (w / 2) * (1 - (r_tst - 0.5) / n)
    ox = fig.label("oxygen", cx1 - w / 2 - 2.0, y_ox, pt=t.note_pt, bold=True, color=AZ, anchor="end", va="middle")
    fig.line(cx1 - hw_tst - 0.6, y_tst, cx1 - w / 2 - 1.2, y_tst, stroke=GREEN_D, lw=0.35)
    sd = fig.label("sleep duration", cx1 - w / 2 - 2.0, y_tst, pt=t.note_pt, bold=True, color=GREEN_D, anchor="end", va="middle")
    pt_ = fig.label("19,173 patients", cx1, top + h + 2.0, pt=t.note_pt, color=GREY, anchor="middle", va="top")
    hi = fig.label("high", cx2 + w / 2 + 2.0, top + 2.0, pt=t.note_pt, color=GREY, anchor="start", va="middle"); lo = fig.label("low", cx2 + 2.0, top + h - 2.0, pt=t.note_pt, color=GREY, anchor="start", va="middle")
    return dict(tst_rank=r_tst, x_extent=[round(min(box(ox)[0], box(sd)[0]), 3), round(max(box(hi)[2], box(lo)[2], cx2 + w / 2), 3)], labels={k: box(v) for k, v in {"title_rank": t1l, "title_risk": t2l, "oxygen": ox, "sleep_duration": sd, "patients": pt_, "high": hi, "low": lo}.items()})


def fig6_tri_r44b(final, *, stem="Main_Fig6_tri_r44b", outdir=OUT):
    fig = new(W6, H6_MM); mid = H6_MM / 2
    tri = triangles(fig, x_left=32.0, top=mid - 30.0, h=60.0)
    geo = fig6_core(fig, cxi=130.0, tx=196.0, organ_arrow=12.0, mid=mid, tw=30.0, th=17.5, pitch=21.0, scale=0.88)
    assert tri["x_extent"][1] + 6.0 < geo["inputs"]["card"][0] and tri["x_extent"][0] > 4.0
    rec = finish(fig, stem, outdir, final=final, attribution=None); rec.update(geometry=geo, triangles=tri, page_pt=[round(W6 / MM, 2), round(H6_MM / MM, 2)], note="the variant keeps its 17.5 mm traces (the same printed depth is not required of an alternative)"); return rec


# ================================================================ Figure 1 panel f (one line), TH 19.5
def band_1f_r44b(final, *, stem="band_1f_r44b", outdir=OUT):
    W = B.B22.W1_PT * MM; H = 46.0
    fig = new(W, H); t = fig.theme; CY = H / 2
    TX, TW, TH = 21.0, 67.0, TH_STD; TY = CY - TH / 2
    r1 = fig.label(B.RANK_CLAIM_R40, TX + TW / 2, TY - 5.5, pt=t.note_pt, bold=True, color=AZ, anchor="middle", va="middle")
    tr = trace(fig, TX, TY, TW, TH, color=AZ, dips=B.DIPS_LOW, lw=0.6)
    fig.line(TX, TY + TH + 1.4, TX + TW, TY + TH + 1.4, stroke=GREY, lw=0.3)
    r3 = fig.label("one night of sleep", TX + TW / 2, TY + TH + 2.6, pt=t.note_pt, color=GREY, anchor="middle", va="top")
    OCX, PITCH = 248.0, 37.0
    ext = B.row_extents(B.ORGANS_LIVE, PITCH, 1.0, family=False); old = B.row_extents(B.ORGANS_OK_R36, PITCH, 1.0); gcx = OCX - ext["ink_offset"]
    n0 = len(fig._placed); B.organ_row_18(fig, gcx, CY, B.ORGANS_BAD, pitch=PITCH, icon_scale=1.0); icons = fig._placed[n0:]
    assert len(icons) == 5 and not any("brain" in p.name for p in icons)
    ink_x = [round(min(p.x for p in icons), 3), round(max(p.x + p.w for p in icons), 3)]
    gap_r36 = (OCX - old["left_of_centre"]) - B.ROW_X_1F_R36; ROW_X = ink_x[0] - gap_r36
    src, tip = (TX + TW + 4.0, CY), (ROW_X, CY); fa = arrow(fig, src, tip)
    ph = B.pill_h(fig); lab = B.plain_label(fig, B.LOW, (src[0] + tip[0]) / 2, CY - (ph + B.ARROW_LW / 2) - ph / 2, was_fill=AZ, was_text_color=B.WHITE)
    clear = B.box_seg_clearance(lab["box"], src, tip) - B.ARROW_LW / 2; lab.update(side="above", clearance_from_arrow_edge_mm=round(clear, 3)); assert clear >= ph - 0.01
    assert lab["box"][2] < ink_x[0] - 5.0 and lab["box"][0] > TX + TW + 4.0 and lab["box"][1] > 0.0 and box(r1)[1] > 0.0 and box(r3)[3] < H
    assert box(r1)[0] > B.B22.LETTER_ZONE_W or box(r1)[1] > B.B22.LETTER_ZONE_H
    rec = finish(fig, stem, outdir, final=final, attribution=None)
    rec.update(page_pt=[round(W / MM, 2), round(H / MM, 2)], trace=tr, arrow=fa, label=lab, organs={"row_y": CY, "ink_x": ink_x, "names": [p.name for p in icons], "pitch": PITCH}, rank_claim=B.RANK_CLAIM_R40, trace_h=TH,
               moved="round 44b: trace height 26 -> 19.5 (the common depth), equal troughs; layout as band_1f_r44"); return rec


# ================================================================ Figure 3 panel d (build_r40.band_3d with TH 19.5)
def band_3d_r44(final, *, stem="band_3d_r44", outdir=OUT, shift=None):
    d = B.SHIFT["band_3d_r40"] if shift is None else shift
    W = B.B23.W3_PT * MM; fig = new(W, B.B23.H3)
    Y_OK, Y_BAD = 23.5 - d, 50.0 - d; TW, TH = B.TW_R31, TH_STD
    Y_FLAT, Y_DIPS = Y_OK - TH / 2, Y_BAD - TH / 2; OCX, PITCH = 262.0, 28.5; TX = B.TX_3D
    bed = fig.place("sleeper-faceB", 36.0, 34.0 - d, w=54)
    dl = fig.label(B.B23.DURATIONS, 36.0, bed.y + bed.h + 0.6, bold=True, anchor="middle", va="top", color=GREEN_D)
    fa = B.fork_arrows(fig, (62.0, 34.0 - d), TX, Y_FLAT, Y_DIPS, TH)
    tp = B.trace_pair(fig, TX, Y_FLAT, Y_DIPS, TW, TH)
    brx = TX + TW + B.ORGAN_ARROW_GAP
    org = B.organ_rows(fig, OCX, (Y_OK, Y_BAD), pitch=PITCH, scale=0.95)
    tip_r34 = OCX - PITCH * 2 - 9.0; gap_r34 = org["first_icon_left_r36"] - tip_r34; tip_x = org["ink_x"][0] - gap_r34
    oa = [arrow(fig, (brx, y), (tip_x, y)) for y in (Y_OK, Y_BAD)]
    assert org["ink_x"][1] < W - 5.0 and tp["top"]["ink_top"] > 1.0 and tp["bottom"]["ink_bottom"] < B.B23.H3 - 1.0, (org["ink_x"], tp["top"]["ink_top"], tp["bottom"]["ink_bottom"])
    rec = finish(fig, stem, outdir, final=final, attribution=None)
    rec.update(shift=round(d, 3), traces=tp, fork=fa, organ_arrows=oa, organs=org, bed_box=box(bed), duration_label=box(dl), trace_h=TH, moved="round 44b: trace height 14 -> 19.5, equal troughs; nothing else (band_3d_r40 geometry)"); return rec


def band_4e_r44(final, *, stem="band_4e_r44", outdir=OUT):
    old_th = B.B45.TH; B.B45.TH = TH_STD
    try: rec = B.band_4e(final, stem=stem, outdir=outdir, shift=B.SHIFT["band_4e_r40"])
    finally: B.B45.TH = old_th
    rec.update(trace_h=TH_STD, moved="round 44b: trace height 15 -> 19.5, equal troughs; nothing else (band_4e_r40 geometry)"); return rec


def band_5d_r44(final, *, stem="band_5d_r44", outdir=OUT):
    old_th = B.B45.TH; B.B45.TH = TH_STD
    try: rec = B.band_5d(final, stem=stem, outdir=outdir, shift=B.SHIFT["band_5d_r40"])
    finally: B.B45.TH = old_th
    rec.update(trace_h=TH_STD, moved="round 44b: trace height 15 -> 19.5, equal troughs; nothing else (band_5d_r40 geometry)"); return rec


# ================================================================ Supp 16 top band (build_r40.band_s16 with TH 10.3)
def band_s16_r44(final, *, stem="band_s16_r44", outdir=OUT):
    W, H = 496.8 * MM, 79.0 * MM; fig = new(W, H); t = fig.theme; AXIS = 14.0; HX = 19.0
    B.house(fig, HX, 10.0, 20.0, 17.0); fig.place("sleeper-faceB", HX, 13.4, w=12)
    l1 = fig.label("SHHS + MrOS", HX, 19.4, pt=t.note_pt, anchor="middle", va="top", bold=True, color=INK)
    l2 = fig.label("5,802 + 2,911 people", HX, 23.0, pt=t.note_pt, anchor="middle", va="top", color=INK)
    TX, TW, TH = 61.0, 19.0, TH_S16; TY = AXIS - TH * 0.5
    tr = trace(fig, TX, TY, TW, TH, color=AZ, dips=(0.15, 0.36, 0.58, 0.80), lw=0.55)
    a1 = arrow(fig, (HX + 10.0 + 1.5, AXIS), (tr["label_box"][0] - 2.0, AXIS)); a2 = arrow(fig, (86.0, AXIS), (112.0, AXIS))
    yl = fig.label("years later", 99.0, AXIS + 2.2, pt=t.note_pt, anchor="middle", va="top", color=INK)
    OCX = 143.75
    ttl = fig.label(B.S16_TITLE.replace("predict future", "predict\nfuture"), OCX, 1.5, pt=10.0, bold=True, anchor="middle", va="top", color=INK21, leading=1.12)
    n0 = len(fig._placed); placed, bottom = B.organ_row_named(fig, OCX, AXIS, pitch=B.S16_PITCH_R40, icon_h=6.5, max_w=B.S16_NAME_MAX_W_R40, name_dy=1.5, organs=list(B.S16_ORGANS_R40)); icons = fig._placed[n0:]
    assert [p.name for p in icons] == ["heart", "human-grey"] and bottom < H - 0.5
    rec = finish(fig, stem, outdir, final=final, attribution=None)
    rec.update(trace=tr, axis=AXIS, arrows={"house_to_trace": a1, "years_later": a2}, title_box=box(ttl), years_box=box(yl), trace_h=TH, moved="round 44b: trace height 12 -> 10.3 (the same printed depth as the double-width sheets), equal troughs; nothing else (band_s16_r40)"); return rec


# ================================================================ Figure 1 panel a (build_r36_fig1a.band ported: merged leg-movement item, organs 0.8, sheet wording)
A = dict(W1_PT=968.66, CAP1="one night of sleep, 141 measurements", CAP2="Which measurement best predicts 54 future diseases and death?",
         LABEL_PT=9.0, LABEL_LEADING=1.15, TW=28.0, TH=11.0, LABEL_GAP_ABOVE=0.8, LABEL_GAP_BELOW=1.5, SLEEPER="sleeper-faceB-tight", TIGHT_W=114.0, FACEB_HEAD=(26.0, 20.5, 8.0),
         TITLE_RAISE_MM=(59.8937 - 39.49) * MM, EXT_MM=2.032, H_MM=79.5, CAL_Y=8.0, LEFT=16.68 * MM, BED_W=58.0, BED_TOP=31.0, GAP=16.0, TOP_Y=12.0, EEG_UP=7.0,
         LEADER_LW=0.4, ARROW_LW=0.9, ARROW_AIR=12.9, ARROW_MAX=60.0, ROW_GAP=2.4, COL_GAP=6.0, BLOCK_TOP_AIR=1.6, BOTTOM_AIR=2.5, ORGAN_SCALE=1.0, TITLE2_CENTRE=283.45)
A_TILES = {"spo2": "blood oxygen\n(pulse oximetry)", "air": "breathing (airflow)", "ecg": "heart rhythm (ECG)", "eeg": "brain waves (EEG)",
           "hyp": "sleep structure and\nleg movements\n(EEG, staging, leg EMG)"}
A_ROWS = [("heart and vessels", [("heart", 13.0, 0.0)]), ("lungs", [("sv-lung", 13.0, 0.0)]), ("metabolism", [("pancreas", 11.0, 0.0)]),
          ("kidney", [("kidney", 11.5, 0.0)]), ("liver", [("liver", 11.5, 0.0)]), ("brain", [("sv-brain-2", 12.0, 0.0)])]
A_ICON_W1 = {"heart": 9.75 / 13.0, "sv-lung": 13.73 / 13.0, "pancreas": 17.86 / 11.0, "kidney": 8.05 / 11.5, "liver": 11.72 / 9.8, "sv-brain-2": 11.78 / 12.0}


def band_a_r44(final, *, stem="band_a_r44", outdir=OUT, merged=True):
    """merged=False (round 44c, Alen: the merge was meant for panel b): the six round-36 items incl. the leg-movement tile, organs at 1.0."""
    L = B.L; APNEA, APNEA_D, GREY_D = L.APNEA, L.APNEA_D, L.GREY_D
    from r14_common import emg_trace
    tiles_txt = dict(A_TILES) if merged else dict(A_TILES, hyp="sleep structure\n(EEG and staging)", emg="leg movements\n(leg EMG)")
    ext, title_raise = A["EXT_MM"], A["TITLE_RAISE_MM"]; h_mm = A["H_MM"] + ext
    cal_y, bed_top, top_y = A["CAL_Y"] + ext, A["BED_TOP"] + ext, A["TOP_Y"] + ext
    title_y = cal_y + 0.1 - title_raise
    W, H = A["W1_PT"] * MM, h_mm; fig = new(W, H); t = fig.theme
    TW, TH, LABEL_PT, LEAD = A["TW"], A["TH"], A["LABEL_PT"], A["LABEL_LEADING"]
    # the merged hyp label is wider than its tile: the scene starts far enough right for the label to stay right of LEFT
    widest = max(L.text_width_mm(line, LABEL_PT, False) for line in tiles_txt["hyp"].split("\n"))
    left_x = max(A["LEFT"], A["LEFT"] + widest / 2 - TW / 2)
    k = A["BED_W"] / A["TIGHT_W"]; bed_h = 62.0 * k; bed_cy = bed_top + bed_h / 2
    bed_x = left_x + TW + A["GAP"]
    bed = fig.place(A["SLEEPER"], bed_x + A["BED_W"] / 2, bed_cy, w=A["BED_W"])
    hx, hy, hr = bed.x + A["FACEB_HEAD"][0] * k, bed.y + A["FACEB_HEAD"][1] * k, A["FACEB_HEAD"][2] * k
    icon = lambda ux, uy: (bed.x + (ux - 4.0) * k, bed.y + (uy - 14.0) * k)
    on_head = lambda deg: (hx + hr * math.cos(math.radians(deg)), hy - hr * math.sin(math.radians(deg)))
    face, head_back, head_low = on_head(52.0), on_head(180.0), on_head(215.0)
    chest, legs, hand = icon(52.0, 35.8), icon(88.0, 44.0), icon(37.0, 51.0)
    mattress_bottom = icon(0.0, 62.0)[1]; bot_y = mattress_bottom + 2.0
    eeg_ty = hy - TH / 2 - A["EEG_UP"]; top_ty = top_y + 7.4
    right_x = bed.x + A["BED_W"] + A["GAP"]; emg_ty = legs[1] - TH / 2 + 2.5
    tiles = [("air", hx - TW / 2, top_ty, APNEA_D, APNEA, "above"), ("ecg", legs[0] - TW / 2 + 3.0, top_ty, INK, GREY, "above"),
             ("eeg", left_x, eeg_ty, INK, GREY, "above"), ("spo2", bed.x + A["BED_W"] / 2 - TW / 2, bot_y, AZ, AZ, "below"), ("hyp", left_x, bot_y, INK, GREY, "above")]
    if not merged: tiles.insert(3, ("emg", right_x, emg_ty, INK, GREY, "above"))
    size_l = LABEL_PT * MM; rec_tiles = []
    for kind, x0, ty, tcol, mcol, side in tiles:
        name = tiles_txt[kind]; n_lines = name.count("\n") + 1
        block_h = size_l * 0.925 + (n_lines - 1) * LEAD * size_l; cx = x0 + TW / 2
        ly = ty - A["LABEL_GAP_ABOVE"] - block_h if side == "above" else ty + TH + A["LABEL_GAP_BELOW"]
        lb = fig.label(name, cx, ly, pt=LABEL_PT, va="top", color=tcol, leading=LEAD, anchor="middle")
        if kind == "spo2": L.spo2_trace(fig, x0, ty, TW, TH, dips=True, color=mcol, ref_at=0.55, label90=None, lw=0.5, dip_centers=(0.15, 0.36, 0.58, 0.80))
        elif kind == "air": B.airflow(fig, x0 + TW / 2, ty + TH / 2, TW, TH * 0.7, mcol, lw=0.5)
        elif kind == "ecg": fig.trace("ecg", x0, ty, TW, TH, color=mcol, lw=0.5)
        elif kind == "eeg": fig.trace("eeg", x0, ty, TW, TH, color=mcol, lw=0.45)
        elif kind == "emg": emg_trace(fig, x0, ty, TW, TH, color=mcol, lw=0.45)
        else: fig.hypnogram(x0, ty + 1.0, TW, TH - 2.0, color=mcol, rem_color=GREY_D, lw=0.5)
        rec_tiles.append({"kind": kind, "label": name, "lines": n_lines, "label_box": [round(lb.x, 2), round(lb.y, 2), round(lb.x + lb.w, 2), round(lb.y + lb.h, 2)], "trace_box": [round(x0, 2), round(ty, 2), round(x0 + TW, 2), round(ty + TH, 2)], "label_side": side})
    def leader(p, q):
        fig.line(p[0], p[1], q[0], q[1], stroke=GREY_D, lw=A["LEADER_LW"]); return [round(v, 2) for v in (*p, *q)]
    T = {r["kind"]: r for r in rec_tiles}; leads = {}
    tb = T["air"]["trace_box"]; leads["air"] = leader(((tb[0] + tb[2]) / 2 + 0.5, tb[3] + 0.8), face)
    tb = T["ecg"]["trace_box"]; leads["ecg"] = leader(((tb[0] + tb[2]) / 2, tb[3] + 0.8), chest)
    tb = T["eeg"]["trace_box"]; leads["eeg"] = leader((tb[2] + 0.8, (tb[1] + tb[3]) / 2), (head_back[0] - 0.3, hy))
    tb = T["hyp"]["trace_box"]; leads["hyp"] = leader((tb[2] + 0.8, (tb[1] + tb[3]) / 2), head_low)
    if not merged: tb = T["emg"]["trace_box"]; leads["emg"] = leader((tb[0] - 0.8, (tb[1] + tb[3]) / 2), legs)
    tb = T["spo2"]["trace_box"]; leads["spo2"] = leader((tb[0] + 6.0, tb[1] - 0.6), hand)
    scene_left = min(min(r["label_box"][0], r["trace_box"][0]) for r in rec_tiles)
    scene_right = max(max(max(r["label_box"][2], r["trace_box"][2]) for r in rec_tiles), bed.x + bed.w)
    c1 = (scene_left + scene_right) / 2; caps = {}
    lb = fig.label(A["CAP1"], c1, title_y, pt=t.label_pt, bold=True, color=INK21, anchor="middle", va="middle")
    caps[1] = {"text": A["CAP1"], "centre": round(c1, 2), "label_box": [round(lb.x, 2), round(lb.y, 2), round(lb.x + lb.w, 2), round(lb.y + lb.h, 2)]}
    rows_def = A_ROWS; cols, rows = 2, 3
    one_line_w = L.text_width_mm(A["CAP2"], t.label_pt, True); size2 = t.label_pt * MM; title_h1 = 0.925 * size2
    title_bottom = cal_y + 0.1 + title_h1 / 2; y_top = title_bottom + A["BLOCK_TOP_AIR"]; y_bot = H - A["BOTTOM_AIR"]; pitch = (y_bot - y_top) / rows
    scale_fill = (pitch - A["ROW_GAP"]) / max(ih for _n, icons in rows_def for _ic, ih, _dx in icons); scale = min(scale_fill, A["ORGAN_SCALE"])
    slot = max(A_ICON_W1[ic] * ih * scale for _n, icons in rows_def for ic, ih, _dx in icons); bw = cols * slot + (cols - 1) * A["COL_GAP"]
    c2 = A["TITLE2_CENTRE"]; block_left, block_right = c2 - bw / 2, c2 + bw / 2
    assert c2 + one_line_w / 2 <= W - 1.0 and c2 - one_line_w / 2 >= caps[1]["label_box"][2] + 3.0, ("title 2 on one line", c2, one_line_w, caps)
    rec_rows = []; icon_left, icon_right = 1e9, -1e9
    for i, (name, icons) in enumerate(rows_def):
        j, r = divmod(i, rows); cx_slot = block_left + j * (slot + A["COL_GAP"]) + slot / 2; cy = y_top + pitch * (r + 0.5); placed = []
        for ic, ih, dx in icons:
            p = fig.place(ic, cx_slot + dx * scale, cy, h=ih * scale); placed.append({"name": ic, "cx": round(p.x + p.w / 2, 2), "cy": round(p.y + p.h / 2, 2), "w": round(p.w, 2), "h": round(p.h, 2)})
            icon_left, icon_right = min(icon_left, p.x), max(icon_right, p.x + p.w)
        rec_rows.append({"name": name, "row": r, "col": j, "cy": round(cy, 2), "icons": placed})
    block_bottom = max(max(ic["cy"] + ic["h"] / 2 for ic in r["icons"]) for r in rec_rows); block_top = min(min(ic["cy"] - ic["h"] / 2 for ic in r["icons"]) for r in rec_rows)
    line_y = (y_top + y_bot) / 2; gap_mid = (scene_right + block_left) / 2; arrow_len = min(A["ARROW_MAX"], (block_left - scene_right) - 2 * A["ARROW_AIR"])
    ax0, ax1 = gap_mid - arrow_len / 2, gap_mid + arrow_len / 2
    fig.arrow((ax0, line_y), (ax1, line_y), color=INK21, lw=A["ARROW_LW"]); yl = fig.label("years", gap_mid, line_y + 2.2, pt=t.note_pt, anchor="middle", va="top", color=INK21)
    lb = fig.label(A["CAP2"], c2, title_y, pt=t.label_pt, bold=True, color=INK21, anchor="middle", va="middle", leading=1.15)
    caps[2] = {"text": A["CAP2"], "centre": round(c2, 2), "label_box": [round(lb.x, 2), round(lb.y, 2), round(lb.x + lb.w, 2), round(lb.y + lb.h, 2)]}
    assert caps[2]["label_box"][3] + A["BLOCK_TOP_AIR"] - 0.6 <= block_top and caps[1]["label_box"][1] >= 0.5 and caps[2]["label_box"][1] >= 0.5
    assert arrow_len >= 25.0, ("the years arrow keeps a run", arrow_len)
    rec = finish(fig, stem, outdir, final=final, attribution=None)
    bottom_labels = max(r["label_box"][3] for r in rec_tiles)
    assert h_mm - bottom_labels >= 2.3 and block_bottom <= h_mm - 2.0 + 1e-6 and abs((block_left + block_right) / 2 - c2) < 0.01 and caps[1]["label_box"][2] + 3.0 < caps[2]["label_box"][0] and caps[2]["label_box"][2] <= W - 0.5
    assert all(r["label_box"][0] >= A["LEFT"] - 0.01 for r in rec_tiles), ("no label left of the scene's left margin", [(r["kind"], r["label_box"][0]) for r in rec_tiles])
    rec.update({"h_mm": h_mm, "page_pt": [round(W / MM, 3), round(H / MM, 3)], "left_x": round(left_x, 3), "scene_shift_mm": round(left_x - A["LEFT"], 3), "tiles": rec_tiles, "leaders": leads, "captions": caps,
                "bed_box": box(bed), "organ_scale": round(scale, 4), "organ_scale_r36": 1.25, "organ_rows": rec_rows, "organ_block": [round(block_left, 2), round(block_top, 2), round(block_right, 2), round(block_bottom, 2)],
                "arrow_span": [round(ax0, 2), round(ax1, 2)], "years_box": box(yl), "deleted": ["the leg-movement tile (emg trace, label, leader)"] if merged else [], "merged_label": tiles_txt["hyp"], "merged": merged,
                "moved": "round 44b: leg movements folded into the sleep-structure item (one three-line label, the hypnogram kept), organs at scale 1.0 (0.8 x 1.25), wording 141 / 54; the scene starts %.2f mm right of the round-36 margin so the wider label keeps the margin" % (left_x - A["LEFT"])})
    return rec


SHEETS = {"band_a_r44": band_a_r44, "band_a_r44c": lambda final: band_a_r44(final, stem="band_a_r44c", merged=False), "band_1f_r44b": band_1f_r44b, "band_3d_r44": band_3d_r44, "band_4e_r44": band_4e_r44, "band_5d_r44": band_5d_r44,
          "band_s16_r44": band_s16_r44, "Main_Fig6_r44b": fig6_r44b, "Main_Fig6_tri_r44b": fig6_tri_r44b}

# ================================================================ round 44d: the triangle variant redrawn (Alen, 2026-09-21 13:20)
GREY_TOP, GREY_BOT = "#8f979e", "#e3e7ea"          # the other measurements: darker near the top of the ranking, lighter below
RED_D2, RED_L2 = "#a3201a", "#fbe9e7"


def triangles_v2(fig, *, x_left, top, h, w=28.0, gap=12.0):
    """Left: an inverted triangle (wide at the top) holding the 141 measurements in rank order, one line each: the oxygen measurements in
    blue at the top, every other measurement in a grey that fades downward. Right: an upright triangle (apex at the top) for disease risk,
    deepest red at the apex, palest at the base, so the top of both triangles reads 'oxygen, highest risk'."""
    rk = pd.read_csv(NUM / "ranking_v3.csv", comment="#"); fam = pd.read_csv(NUM / "measure_families.csv").set_index("feature")["family"]
    rk["family"] = rk.feature.map(fam); assert len(rk) == 141 and rk.family.notna().all()
    order = rk.sort_values("rank").reset_index(drop=True); t = fig.theme; n = len(order); step = h / n
    cx1 = x_left + w / 2; cx2 = cx1 + w + gap
    n_ox = int((order.family == "Oxygenation").sum()); last_ox = int(order.index[order.family == "Oxygenation"].max())     # rank index of the lowest oxygen measure
    for i, r in order.iterrows():
        tt = (i + 0.5) / n; y = top + tt * h; hw = (w / 2) * (1 - tt)
        col = AZ if r.family == "Oxygenation" else lerp_hex(GREY_TOP, GREY_BOT, tt)
        if hw > 0.12: fig.line(cx1 - hw, y, cx1 + hw, y, stroke=col, lw=step * 1.02, cap="butt")
    fig.polygon([(cx1 - w / 2, top), (cx1 + w / 2, top), (cx1, top + h)], fill=None, stroke=GREY, lw=0.3)
    NS = 96
    for j in range(NS):
        t0, t1 = j / NS, (j + 1) / NS
        fig.polygon([(cx2 - (w / 2) * t0, top + t0 * h), (cx2 + (w / 2) * t0, top + t0 * h), (cx2 + (w / 2) * t1, top + t1 * h), (cx2 - (w / 2) * t1, top + t1 * h)], fill=lerp_hex(RED_D2, RED_L2, (t0 + t1) / 2), stroke=None)
    fig.polygon([(cx2, top), (cx2 + w / 2, top + h), (cx2 - w / 2, top + h)], fill=None, stroke=GREY, lw=0.3)
    # titles
    t1l = fig.label("141 sleep measurements", cx1, top - 3.0, pt=t.note_pt, bold=True, color=INK21, anchor="middle", va="bottom")
    t2l = fig.label("disease risk", cx2, top - 3.0, pt=t.note_pt, bold=True, color=INK21, anchor="middle", va="bottom")
    # left-side brackets: oxygen (the top of the ranking) and the rest
    y_ox0, y_ox1 = top + 0.3, top + ((last_ox + 1) / n) * h
    xb = cx1 - w / 2 - 1.6
    fig.line(xb, y_ox0, xb, y_ox1, stroke=AZ, lw=0.6, cap="butt")
    ox = fig.label("oxygen", xb - 1.5, (y_ox0 + y_ox1) / 2, pt=t.note_pt, bold=True, color=AZ, anchor="end", va="middle")
    y_r0, y_r1 = y_ox1 + 1.0, top + h * 0.97
    hw_mid = lambda yy: (w / 2) * (1 - (yy - top) / h)
    fig.line(xb, y_r0, xb, y_r1, stroke=GREY_TOP, lw=0.6, cap="butt")
    rest = fig.label("all other\nsleep\nmeasurements", xb - 1.5, (y_r0 + y_r1) / 2, pt=t.note_pt, color=GREY_TOP, anchor="end", va="middle", leading=1.15)
    pt_ = fig.label("19,173 patients", cx1, top + h + 2.2, pt=t.note_pt, color=GREY, anchor="middle", va="top")
    # risk marks in the gap, left of the upright triangle
    hi = fig.label("high", cx2 - 4.2, top + 2.6, pt=t.note_pt, color=RED_D2, bold=True, anchor="end", va="middle")
    lo = fig.label("low", cx2 - w / 2 - 2.2, top + h - 2.0, pt=t.note_pt, color=GREY, anchor="end", va="middle")
    for lb in (hi, lo): assert lb.x > cx1 + hw_mid(lb.y + lb.h / 2) + 1.0, ("gap labels clear of the left triangle", box(lb))
    return dict(n_oxygen=n_ox, oxygen_zone_ranks=[1, last_ox + 1], x_extent=[round(min(box(ox)[0], box(rest)[0]), 3), round(cx2 + w / 2, 3)],
                labels={k: box(v) for k, v in {"title_rank": t1l, "title_risk": t2l, "oxygen": ox, "rest": rest, "patients": pt_, "high": hi, "low": lo}.items()},
                colours=dict(oxygen=AZ, rest=[GREY_TOP, GREY_BOT], risk=[RED_D2, RED_L2]), shapes="left inverted (wide top), right upright (apex top)")


def fig6_tri_r44d(final, *, stem="Main_Fig6_tri_r44d", outdir=OUT):
    fig = new(W6, H6_MM); mid = H6_MM / 2
    tri = triangles_v2(fig, x_left=34.0, top=mid - 30.0, h=60.0)
    geo = fig6_core(fig, cxi=130.0, tx=196.0, organ_arrow=12.0, mid=mid, tw=30.0, th=TH_STD, pitch=21.0, scale=0.88)
    assert tri["x_extent"][1] + 6.0 < geo["inputs"]["card"][0] and tri["x_extent"][0] > 4.0, (tri["x_extent"], geo["inputs"]["card"])
    rec = finish(fig, stem, outdir, final=final, attribution=None)
    rec.update(geometry=geo, triangles=tri, page_pt=[round(W6 / MM, 2), round(H6_MM / MM, 2)], trace_h=TH_STD, moved="round 44d: triangles redrawn (inverted ranking triangle with an oxygen cap, upright risk triangle deepest at the apex)")
    return rec


SHEETS["Main_Fig6_tri_r44d"] = fig6_tri_r44d



# ================================================================ round 44e: rounded gradient triangles (Alen, 2026-09-23: "prettier, rounded, colours inside")
def rounded_tri_path(pts, r):
    """SVG path of a polygon whose corners are true circular arcs of radius r tangent to both edges (the tangent points sit
    r / tan(angle/2) from each vertex, so an acute apex rounds visibly)."""
    import math as _m
    n = len(pts); segs = []
    for k in range(n):
        p0, p1, p2 = pts[k - 1], pts[k], pts[(k + 1) % n]
        v0 = (p0[0] - p1[0], p0[1] - p1[1]); v2 = (p2[0] - p1[0], p2[1] - p1[1])
        l0, l2 = _m.hypot(*v0), _m.hypot(*v2); u0 = (v0[0] / l0, v0[1] / l0); u2 = (v2[0] / l2, v2[1] / l2)
        ang = _m.acos(max(-1.0, min(1.0, u0[0] * u2[0] + u0[1] * u2[1]))); d = min(r / _m.tan(ang / 2), l0 / 2.2, l2 / 2.2)
        a = (p1[0] + u0[0] * d, p1[1] + u0[1] * d); b = (p1[0] + u2[0] * d, p1[1] + u2[1] * d)
        rr = d * _m.tan(ang / 2)
        cross = u0[0] * u2[1] - u0[1] * u2[0]; sweep = 1 if cross < 0 else 0
        segs.append((a, b, rr, sweep))
    out = [f"M{segs[0][0][0]:.3f} {segs[0][0][1]:.3f}"]
    for k in range(n):
        a, b, rr, sweep = segs[k]
        out.append(f"A{rr:.3f} {rr:.3f} 0 0 {sweep} {b[0]:.3f} {b[1]:.3f}")
        a_next = segs[(k + 1) % n][0]; out.append(f"L{a_next[0]:.3f} {a_next[1]:.3f}")
    return " ".join(out) + " Z"


def triangles_v3(fig, *, x_left, top, h, w=28.0, gap=12.0, r=3.0):
    """Left: rounded inverted triangle, a blue cap for the oxygen measurements over a grey fill that fades downward (all other measurements).
    Right: rounded upright triangle, deepest red at the apex fading to pale at the base (disease risk). Smooth SVG gradients, clipped fills."""
    rk = pd.read_csv(NUM / "ranking_v3.csv", comment="#"); fam = pd.read_csv(NUM / "measure_families.csv").set_index("feature")["family"]
    rk["family"] = rk.feature.map(fam); order = rk.sort_values("rank").reset_index(drop=True); n = len(order); t = fig.theme
    last_ox = int(order.index[order.family == "Oxygenation"].max()); cap_frac = (last_ox + 1) / n
    cx1 = x_left + w / 2; cx2 = cx1 + w + gap
    # the rounded apex sits e mm inside the geometric apex: extend both geometric apexes by e so the visible shapes span exactly top..top+h
    import math as _m
    e = 0.0
    for _ in range(30):
        half = _m.atan((w / 2) / (h + e)); e = r * (1 / _m.sin(half) - 1)
    L = [(cx1 - w / 2, top), (cx1 + w / 2, top), (cx1, top + h + e)]; R = [(cx2, top - e), (cx2 + w / 2, top + h), (cx2 - w / 2, top + h)]
    pL, pR = rounded_tri_path(L, r), rounded_tri_path(R, r)
    cap_h = h * cap_frac + 0.4
    fig.add(f'<defs><linearGradient id="gRest" x1="0" y1="0" x2="0" y2="1"><stop offset="0" stop-color="{GREY_TOP}"/><stop offset="1" stop-color="{GREY_BOT}"/></linearGradient>'
            f'<linearGradient id="gRisk" x1="0" y1="0" x2="0" y2="1"><stop offset="0" stop-color="{RED_D2}"/><stop offset="0.55" stop-color="#d9675c"/><stop offset="1" stop-color="{RED_L2}"/></linearGradient>'
            f'<clipPath id="clipL"><path d="{pL}"/></clipPath></defs>')
    fig.add(f'<path d="{pL}" fill="url(#gRest)" stroke="none"/>')
    fig.add(f'<rect x="{cx1 - w / 2:.3f}" y="{top:.3f}" width="{w:.3f}" height="{cap_h:.3f}" fill="{AZ}" clip-path="url(#clipL)"/>')
    fig.add(f'<path d="{pL}" fill="none" stroke="#ffffff" stroke-width="0.5"/>')
    fig.add(f'<path d="{pR}" fill="url(#gRisk)" stroke="none"/>')
    fig.add(f'<path d="{pR}" fill="none" stroke="#ffffff" stroke-width="0.5"/>')
    t1l = fig.label("141 sleep measurements", cx1, top - 3.0, pt=t.note_pt, bold=True, color=INK21, anchor="middle", va="bottom")
    t2l = fig.label("disease risk", cx2, top - 3.0, pt=t.note_pt, bold=True, color=INK21, anchor="middle", va="bottom")
    xb = cx1 - w / 2 - 1.6
    fig.line(xb, top + 0.3, xb, top + cap_h - 0.2, stroke=AZ, lw=0.7, cap="round")
    ox = fig.label("oxygen", xb - 1.5, top + cap_h / 2, pt=t.note_pt, bold=True, color=AZ, anchor="end", va="middle")
    y_r0, y_r1 = top + cap_h + 1.0, top + h - 0.8
    fig.line(xb, y_r0, xb, y_r1, stroke=INK21, lw=0.7, cap="round")
    rest = fig.label("all other\nsleep\nmeasurements", xb - 1.5, (y_r0 + y_r1) / 2, pt=t.note_pt, color=INK21, anchor="end", va="middle", leading=1.15)
    pt_ = fig.label("19,173 patients", cx1, top + h + 2.2, pt=t.note_pt, color=INK21, anchor="middle", va="top")
    hi = fig.label("high", cx2 - 4.2, top + 2.6, pt=t.note_pt, color=RED_D2, bold=True, anchor="end", va="middle")
    lo = fig.label("low", cx2 - w / 2 - 2.2, top + h - 2.0, pt=t.note_pt, color=INK21, anchor="end", va="middle")
    hw_at = lambda yy: (w / 2) * (1 - (yy - top) / h)
    for lb in (hi, lo): assert lb.x > cx1 + hw_at(lb.y + lb.h / 2) + 1.0, ("gap labels clear of the left triangle", box(lb))
    return dict(oxygen_zone_ranks=[1, last_ox + 1], cap_mm=round(cap_h, 2), corner_r=r, apex_extension_mm=round(e, 2), x_extent=[round(min(box(ox)[0], box(rest)[0]), 3), round(cx2 + w / 2, 3)],
                labels={k: box(v) for k, v in {"title_rank": t1l, "title_risk": t2l, "oxygen": ox, "rest": rest, "patients": pt_, "high": hi, "low": lo}.items()},
                colours=dict(oxygen=AZ, rest=[GREY_TOP, GREY_BOT], risk=[RED_D2, "#d9675c", RED_L2]), shapes="rounded, gradient fills (SVG linearGradient), cap clipped to the triangle")


def fig6_tri_r44e(final, *, stem="Main_Fig6_tri_r44e", outdir=OUT):
    fig = new(W6, H6_MM); mid = H6_MM / 2
    tri = triangles_v3(fig, x_left=34.0, top=mid - 23.0, h=46.0)
    geo = fig6_core(fig, cxi=130.0, tx=196.0, organ_arrow=12.0, mid=mid, tw=30.0, th=TH_STD, pitch=21.0, scale=0.88)
    assert tri["x_extent"][1] + 6.0 < geo["inputs"]["card"][0] and tri["x_extent"][0] > 4.0
    rec = finish(fig, stem, outdir, final=final, attribution=None)
    rec.update(geometry=geo, triangles=tri, page_pt=[round(W6 / MM, 2), round(H6_MM / MM, 2)], trace_h=TH_STD, moved="round 44e: rounded gradient triangles"); return rec


SHEETS["Main_Fig6_tri_r44e"] = fig6_tri_r44e


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--final", action="store_true"); ap.add_argument("--only", default=""); args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True); (WORK / "previews").mkdir(parents=True, exist_ok=True)
    lock = open(B.LOCK, "w"); fcntl.flock(lock, fcntl.LOCK_EX); t0 = time.time()
    try:
        only = [s for s in args.only.split(",") if s]
        for stem, fn in SHEETS.items():
            if only and stem not in only: continue
            r = fn(args.final); r["render"] = dict(B.RENDER_INFO)
            (WORK / f"records_r44b_{stem}.json").write_text(json.dumps(r, indent=1, default=str))
            print(json.dumps({k: r.get(k) for k in ("stem", "w_pt", "h_mm", "audit", "min_pt", "texts", "numbers")}), flush=True)
    finally:
        fcntl.flock(lock, fcntl.LOCK_UN); lock.close()
    print(f"done in {time.time() - t0:.0f} s", flush=True)


