#!/usr/bin/env python3
"""Round 44 (2026-09-21), Alen's three sketch comments, on the round-40 engine (build_r40.py imported: same build tree, icons,
colours, trace/arrow/organ-row helpers, Chrome render under the shared lock):

  Main_Fig6_r44      968.94 x H6 pt   Figure 6 as his mock: an inputs card on the LEFT with two rows (sleeper "sleep duration" in
                                      green, breathing trace + sleeper "sleep apnea" in orange; the PAP mask left out on purpose),
                                      a grey fork to two "90% SpO2" traces (two dips / six dips), a grey arrow from each trace to a
                                      card holding the four-group organ row (coloured on top, grey below). No other words.
  Main_Fig6_tri_r44  968.94 x H6 pt   the same, with two inverted triangles at the far left: the 141 measurements in rank order
                                      from the top (one line per measurement, coloured by family, "oxygen" at the top and "sleep
                                      duration" at rank 127 marked) and "disease risk" as a red gradient, darkest at the top.
  band_1f_r44        968.66 x H1 pt   Figure 1 panel f as ONE line: the rank claim, a trace with six dips over the 90% line, one
                                      grey arrow labelled "Low oxygen", the grey organ row only.

    python3 build_r44.py --final [--only stem,stem]
"""
from __future__ import annotations
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import argparse, fcntl, json, math, sys, time
from pathlib import Path
R40 = Path(f"{paths.FIGURE_ROOT}/ROUND40_2026-09-18/lanes/LS_SKETCH/scripts")
R44 = Path(f"{paths.FIGURE_ROOT}/ROUND44_2026-09-21/lanes/LSKETCH")
sys.path.insert(0, str(R40))
import build_r40 as B                                             # noqa: E402  (the engine, helpers, constants; its own outputs untouched)
import pandas as pd                                               # noqa: E402
OUT = R44 / "work" / "build"; WORK = R44 / "work"
B.L.PREVIEWS = WORK / "previews"
new, finish, trace, arrow, organ_row_f = B.new, B.finish, B.trace, B.arrow, B.organ_row_f
AZ, CARD, GREEN_D, GREY, GREY_L, INK21, AMBER, MM = B.AZ, B.CARD, B.GREEN_D, B.GREY, B.GREY_L, B.INK21, B.AMBER, B.MM_PER_PT
NUM = Path(paths.NUMBERS_DIR)
FAMILY_COL = {"Oxygenation": "#0288d1", "Heart rate and variability": "#8d2dfa", "Breathing events": "#d55e00", "Brain, microstructure": "#fa93a5",
              "Brain, spectral power": "#f2e279", "Sleep timing and structure": "#39c445", "Limb movements": "#8a9099"}   # Figure 1c's key (V14_L1_RANK/common.py)
RED_D, RED_L = "#b3261e", "#fbe4e1"
W6 = B.W6_PT * MM                                                 # 341.82 mm
H6_MM = 100.0                                                     # the new Figure 6 height (283.5 pt); the coordinator adds the 16 pt title strip
PITCH_IN, BED_W = 43.0, 26.0
BED_H = BED_W * 82.0 / 122.0
TW6, TH6 = 34.0, 19.5
ROW_DY = 24.0                                                     # the two trace rows sit MID -/+ ROW_DY
PANEL_H = 2 * (B.ICON_HALF + B.PANEL_PAD)                          # 20.5 mm


def lerp_hex(a, b, t):
    ca = [int(a[i:i + 2], 16) for i in (1, 3, 5)]; cb = [int(b[i:i + 2], 16) for i in (1, 3, 5)]
    return "#" + "".join(f"{round(x + (y - x) * t):02x}" for x, y in zip(ca, cb))


def box(p): return B.box(p)


def inputs_card(fig, cxi, mid, *, box_pad=8.0):
    """The two-row inputs card centred on (cxi, mid): sleeper + 'sleep duration' (green), breathing trace + sleeper + 'sleep apnea' (orange)."""
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
    """One card with the four-group family row, ink centred on pcx; returns the card rect and the icon extents."""
    pad = B.PANEL_PAD if pad is None else pad
    ext = B.row_extents(B.ORGANS_LIVE, pitch, scale, family=True)
    n0 = len(fig._placed)
    organ_row_f(fig, pcx - ext["ink_offset"], pcy, dark=dark, pitch=pitch, scale=scale, names=False)
    icons = fig._placed[n0:]
    assert len(icons) == 5 and not any("brain" in p.name for p in icons), [p.name for p in icons]
    x0, x1 = min(p.x for p in icons), max(p.x + p.w for p in icons)
    assert abs((x0 + x1) / 2 - pcx) < 0.01, ("ink centred", x0, x1, pcx)
    half = max(pcx - x0, x1 - pcx) + pad; ph = 2 * (B.ICON_HALF * scale + pad)
    rect = [round(pcx - half, 3), round(pcy - ph / 2, 3), round(2 * half, 3), round(ph, 3)]
    return dict(rect=rect, ink_x=[round(x0, 3), round(x1, 3)], names=[p.name for p in icons], pitch=pitch, scale=scale)


def fig6_core(fig, *, cxi, tx, organ_arrow, mid, tw=TW6, th=TH6, pitch=25.0, scale=1.0):
    """Card -> fork -> two traces -> arrows -> two organ cards. Returns the geometry; the cards are drawn first (behind)."""
    n_parts0 = len(fig._parts)
    inp = inputs_card(fig, cxi, mid)
    rx1 = inp["card"][0] + inp["card"][2]
    y_top, y_bot = mid - ROW_DY, mid + ROW_DY
    tr_top = trace(fig, tx, y_top - th / 2, tw, th, color=AZ, dips=B.DIPS_PRESERVED, lw=1.0, label_pt=fig.theme.note_pt, label_bold=True, label_gap=0.7, ref_lw=0.5, ref_dash="1.4 1.2", label_side="left")
    tr_bot = trace(fig, tx, y_bot - th / 2, tw, th, color=AZ, dips=B.DIPS_LOW, lw=1.0, label_pt=fig.theme.note_pt, label_bold=True, label_gap=0.7, ref_lw=0.5, ref_dash="1.4 1.2", label_side="left")
    tips = [(tr_top["label_box"][0] - B.HEAD_CLEAR, y_top), (tr_bot["label_box"][0] - B.HEAD_CLEAR, y_bot)]
    fork = [arrow(fig, (rx1, mid), tp) for tp in tips]
    ax0 = tx + tw + B.ORGAN_ARROW_GAP
    px0 = ax0 + organ_arrow + B.HEAD_CLEAR
    ext = B.row_extents(B.ORGANS_LIVE, pitch, scale, family=True)
    half = max(ext["left_of_centre"], ext["right_of_centre"]) + B.PANEL_PAD
    # ink-centred cards need the offset: card left edge = pcx - half; solve pcx from px0
    pcx = px0 + half
    oa = [arrow(fig, (ax0, y_top), (px0 - B.HEAD_CLEAR, y_top)), arrow(fig, (ax0, y_bot), (px0 - B.HEAD_CLEAR, y_bot))]
    pan_top = organ_panel(fig, pcx, y_top, dark=False, pitch=pitch, scale=scale)
    pan_bot = organ_panel(fig, pcx, y_bot, dark=True, pitch=pitch, scale=scale)
    assert pan_top["ink_x"] == pan_bot["ink_x"], (pan_top, pan_bot)
    n1 = len(fig._parts)
    for r in (inp["card"], pan_top["rect"], pan_bot["rect"]):
        fig.rect(*r, rx=B.PANEL_RX, fill=CARD, stroke=GREY_L, lw=0.3)
    moved = fig._parts[n1:]; del fig._parts[n1:]; fig._parts[n_parts0:n_parts0] = moved     # cards behind everything drawn since
    for tp in tips: assert tp[0] - rx1 >= 25.0, ("fork run at least 25 mm", tp, rx1)
    assert abs(tr_top["y90"] - y_top) < 0.01 and abs(tr_bot["y90"] - y_bot) < 0.01
    right_edge = pan_top["rect"][0] + pan_top["rect"][2]
    assert right_edge < W6 - 5.0 and inp["card"][0] > 5.0, ("inside the page with 5 mm margins", inp["card"], right_edge)
    return dict(inputs=inp, traces={"top": tr_top, "bottom": tr_bot}, fork=fork, organ_arrows=oa, panels={"top": pan_top, "bottom": pan_bot}, rows_y=[y_top, y_bot],
                content_x=[inp["card"][0], round(right_edge, 3)], mid=mid)


def fig6_r44(final, *, stem="Main_Fig6_r44", outdir=OUT):
    fig = new(W6, H6_MM); mid = H6_MM / 2
    geo = fig6_core(fig, cxi=61.0, tx=143.0, organ_arrow=20.0, mid=mid)
    cx_content = (geo["content_x"][0] + geo["content_x"][1]) / 2
    assert abs(cx_content - W6 / 2) < 6.0, ("content roughly centred", cx_content, W6 / 2)
    rec = finish(fig, stem, outdir, final=final, attribution=None)
    rec.update(geometry=geo, page_pt=[round(W6 / MM, 2), round(H6_MM / MM, 2)], deleted=["PAP mask and 'PAP treatment'", "'Preserved oxygen' / 'Low oxygen' labels (the mock carries none)"],
               moved="round 44: Alen's mock, inputs card on the left with two rows, fork to two traces, arrows to two organ cards")
    return rec


def triangles(fig, *, x_left, top, h, w=26.0, gap=10.0):
    """Two inverted triangles: the 141 measurements in rank order (one line each, family colours), and disease risk (red gradient)."""
    rk = pd.read_csv(NUM / "ranking_v3.csv", comment="#"); fam = pd.read_csv(NUM / "measure_families.csv").set_index("feature")["family"]
    rk["family"] = rk.feature.map(fam); assert len(rk) == 141 and rk.family.notna().all() and sorted(rk["rank"]) == list(range(1, 142))
    order = rk.sort_values("rank")
    t = fig.theme
    cx1 = x_left + w / 2; cx2 = cx1 + w + gap
    n = len(order); step = h / n
    for i, (_, r) in enumerate(order.iterrows()):
        tt = (i + 0.5) / n; y = top + tt * h; hw = (w / 2) * (1 - tt)
        if hw > 0.15: fig.line(cx1 - hw, y, cx1 + hw, y, stroke=FAMILY_COL[r.family], lw=step * 0.92, cap="butt")
    fig.polygon([(cx1 - w / 2, top), (cx1 + w / 2, top), (cx1, top + h)], fill=None, stroke=GREY, lw=0.35)
    NS = 48
    for j in range(NS):
        t0, t1 = j / NS, (j + 1) / NS
        pts = [(cx2 - (w / 2) * (1 - t0), top + t0 * h), (cx2 + (w / 2) * (1 - t0), top + t0 * h), (cx2 + (w / 2) * (1 - t1), top + t1 * h), (cx2 - (w / 2) * (1 - t1), top + t1 * h)]
        fig.polygon(pts, fill=lerp_hex(RED_D, RED_L, (t0 + t1) / 2), stroke=None)
    fig.polygon([(cx2 - w / 2, top), (cx2 + w / 2, top), (cx2, top + h)], fill=None, stroke=GREY, lw=0.35)
    t1l = fig.label("141 sleep measurements", cx1, top - 2.5, pt=t.note_pt, bold=True, color=INK21, anchor="middle", va="bottom", max_w=w + 14.0, leading=1.15)
    t2l = fig.label("disease risk", cx2, top - 2.5, pt=t.note_pt, bold=True, color=INK21, anchor="middle", va="bottom")
    # side marks on the ranking triangle: oxygen at the top, sleep duration at its rank
    r_tst = int(order[order.feature == "TST_min"]["rank"].iloc[0]); assert r_tst == 127, r_tst
    y_ox = top + (2.0 / n) * h; y_tst = top + ((r_tst - 0.5) / n) * h
    hw_tst = (w / 2) * (1 - (r_tst - 0.5) / n)
    ox = fig.label("oxygen", cx1 - w / 2 - 2.0, y_ox, pt=t.note_pt, bold=True, color=AZ, anchor="end", va="middle")
    fig.line(cx1 - hw_tst - 0.6, y_tst, cx1 - w / 2 - 1.2, y_tst, stroke=GREEN_D, lw=0.35)
    sd = fig.label("sleep duration", cx1 - w / 2 - 2.0, y_tst, pt=t.note_pt, bold=True, color=GREEN_D, anchor="end", va="middle")
    pt_ = fig.label("19,173 patients", cx1, top + h + 2.0, pt=t.note_pt, color=GREY, anchor="middle", va="top")
    hi = fig.label("high", cx2 + w / 2 + 2.0, top + 2.0, pt=t.note_pt, color=GREY, anchor="start", va="middle")
    lo = fig.label("low", cx2 + 2.0, top + h - 2.0, pt=t.note_pt, color=GREY, anchor="start", va="middle")
    return dict(ranking={"cx": cx1, "w": w, "top": top, "h": h, "n": n, "families": order.family.value_counts().to_dict(), "top5": order.feature.head(5).tolist(), "tst_rank": r_tst},
                risk={"cx": cx2, "w": w, "top": top, "h": h, "colours": [RED_D, RED_L], "slices": NS},
                labels={k: box(v) for k, v in {"title_rank": t1l, "title_risk": t2l, "oxygen": ox, "sleep_duration": sd, "patients": pt_, "high": hi, "low": lo}.items()},
                x_extent=[round(min(box(ox)[0], box(sd)[0]), 3), round(max(box(hi)[2], box(lo)[2], cx2 + w / 2), 3)])


def fig6_tri_r44(final, *, stem="Main_Fig6_tri_r44", outdir=OUT):
    fig = new(W6, H6_MM); mid = H6_MM / 2
    tri = triangles(fig, x_left=32.0, top=mid - 30.0, h=60.0)
    geo = fig6_core(fig, cxi=130.0, tx=196.0, organ_arrow=12.0, mid=mid, tw=30.0, th=17.5, pitch=21.0, scale=0.88)
    assert tri["x_extent"][1] + 6.0 < geo["inputs"]["card"][0], ("triangles clear of the card", tri["x_extent"], geo["inputs"]["card"])
    assert tri["x_extent"][0] > 4.0, tri["x_extent"]
    rec = finish(fig, stem, outdir, final=final, attribution=None)
    rec.update(geometry=geo, triangles=tri, page_pt=[round(W6 / MM, 2), round(H6_MM / MM, 2)],
               moved="round 44: the mock plus two inverted triangles (141 measurements by rank in family colours, disease risk red gradient darkest at the top)")
    return rec


def band_1f_r44(final, *, stem="band_1f_r44", outdir=OUT):
    """Figure 1 panel f as one line: rank claim, six-dip trace, one arrow 'Low oxygen', the grey organ row."""
    W = B.B22.W1_PT * MM; H = 46.0
    fig = new(W, H); t = fig.theme
    CY = H / 2
    TX, TW, TH = 21.0, 67.0, 26.0
    TY = CY - TH / 2
    r1 = fig.label(B.RANK_CLAIM_R40, TX + TW / 2, TY - 5.5, pt=t.note_pt, bold=True, color=AZ, anchor="middle", va="middle")
    tr = trace(fig, TX, TY, TW, TH, color=AZ, dips=B.DIPS_LOW, lw=0.6)
    fig.line(TX, TY + TH + 1.4, TX + TW, TY + TH + 1.4, stroke=GREY, lw=0.3)
    r3 = fig.label("one night of sleep", TX + TW / 2, TY + TH + 2.6, pt=t.note_pt, color=GREY, anchor="middle", va="top")
    OCX, PITCH = 248.0, 37.0
    ext = B.row_extents(B.ORGANS_LIVE, PITCH, 1.0, family=False)
    old = B.row_extents(B.ORGANS_OK_R36, PITCH, 1.0)
    gcx = OCX - ext["ink_offset"]
    n0 = len(fig._placed)
    B.organ_row_18(fig, gcx, CY, B.ORGANS_BAD, pitch=PITCH, icon_scale=1.0)
    icons = fig._placed[n0:]; assert len(icons) == 5 and not any("brain" in p.name for p in icons), [p.name for p in icons]
    ink_x = [round(min(p.x for p in icons), 3), round(max(p.x + p.w for p in icons), 3)]
    gap_r36 = (OCX - old["left_of_centre"]) - B.ROW_X_1F_R36              # 7.125 mm, the tip-to-first-icon gap of the two-row band
    ROW_X = ink_x[0] - gap_r36
    src, tip = (TX + TW + 4.0, CY), (ROW_X, CY)
    fa = arrow(fig, src, tip)
    ph = B.pill_h(fig); cx_ = (src[0] + tip[0]) / 2; cy_ = CY - (ph + B.ARROW_LW / 2) - ph / 2
    lab = B.plain_label(fig, B.LOW, cx_, cy_, was_fill=AZ, was_text_color=B.WHITE)
    clear = B.box_seg_clearance(lab["box"], src, tip) - B.ARROW_LW / 2
    lab.update(side="above", clearance_from_arrow_edge_mm=round(clear, 3), pill_clearance_was_mm=round(ph, 3))
    assert clear >= ph - 0.01, lab
    assert lab["box"][2] < ink_x[0] - 5.0 and lab["box"][0] > TX + TW + 4.0, ("label between the fork and the organs", lab["box"], ink_x)
    assert lab["box"][1] > 0.0 and box(r1)[1] > 0.0 and box(r3)[3] < H, ("inside the band", lab["box"], box(r1), box(r3))
    assert box(r1)[0] > B.B22.LETTER_ZONE_W or box(r1)[1] > B.B22.LETTER_ZONE_H, "the letter zone stays clear"
    assert ink_x[1] < W - 5.0, (ink_x, W)
    rec = finish(fig, stem, outdir, final=final, attribution=None)
    rec.update(page_pt=[round(W / MM, 2), round(H / MM, 2)], trace=tr, arrow=fa, label=lab, organs={"row_y": CY, "ink_x": ink_x, "names": [p.name for p in icons], "pitch": PITCH, "grid_centre": round(gcx, 3), "ink_centre": OCX},
               left_block={"rank": box(r1), "night": box(r3)}, rank_claim=B.RANK_CLAIM_R40, letter_zone=[0.0, 0.0, B.B22.LETTER_ZONE_W, B.B22.LETTER_ZONE_H],
               deleted=["the coloured organ row", "'Preserved oxygen'", "the second fork arrow"], dips=list(B.DIPS_LOW),
               moved="round 44: one line (Alen): six dips, one arrow 'Low oxygen', grey organs only; band 46 mm tall")
    return rec


SHEETS = {"Main_Fig6_r44": fig6_r44, "Main_Fig6_tri_r44": fig6_tri_r44, "band_1f_r44": band_1f_r44}

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--final", action="store_true"); ap.add_argument("--only", default=""); args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True); (WORK / "previews").mkdir(parents=True, exist_ok=True)
    lock = open(B.LOCK, "w"); fcntl.flock(lock, fcntl.LOCK_EX); t0 = time.time()
    try:
        only = [s for s in args.only.split(",") if s]
        for stem, fn in SHEETS.items():
            if only and stem not in only: continue
            r = fn(args.final); r["render"] = dict(B.RENDER_INFO)
            (WORK / f"records_r44_{stem}.json").write_text(json.dumps(r, indent=1, default=str))
            print(json.dumps({k: r[k] for k in ("stem", "w_pt", "h_mm", "audit", "min_pt", "texts", "numbers")}), flush=True)
    finally:
        fcntl.flock(lock, fcntl.LOCK_UN); lock.close()
    print(f"done in {time.time() - t0:.0f} s", flush=True)
