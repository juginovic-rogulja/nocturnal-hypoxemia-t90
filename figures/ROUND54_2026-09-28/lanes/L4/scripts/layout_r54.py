#!$T90_PY
"""ROUND 54 (2026-09-29, lane L4): the ONE layout spec of Main_Fig4 at the round-54 sizes (ticks, row labels, keys, cell values, HR and q
columns 14 pt, axis titles and column headers 15, small notes and asterisk keys 13, panel letters and the title 18). Every box, baseline
and column edge is computed HERE from the measured widths of the real strings (matplotlib, the builders' own renderer, Arial) and the
pitch rules, written to work/layout_r54.json and read by the builders, the composes and the verifiers (no script carries a page pin of
its own). The V13 pins of round 49 no longer apply: the sheet is re-laid.

Owner's decision (option A, 2026-09-29): Main_Fig4 = panels a and b side by side + the sketch band as panel c (letters a, b, c, title on
baseline 18); the counts panel (V30's panel c) is delivered as a standalone page (no letter, no title: it becomes panel c of Extended
Data Fig. 7); V30's panel d is built for the record only. Rules: panel b (42 sub-rows, the taller forest) sets the row block: sub-row
pitch 15.0 pt, pair gap 2.5 pt, the block header line one sub-row plus a pair gap; panel a's pitch is what fills the same block so that
the two axis rules sit on one line (the round-30 rule, the roles reversed); plots at least as wide as their tick labels need
(neighbouring tick labels 3 pt apart); row labels on two lines only where the sheet width forces it (the widest first, each wrap declared).
usage: layout_r54.py            (writes work/layout_r54.json and prints the geometry table)"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import json, math, os, re, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import l4lib as L

R54 = f"{paths.FIGURE_ROOT}/ROUND54_2026-09-28"; LANE = f"{R54}/lanes/L4"; WORK = f"{LANE}/Main_Fig4/work"; R49L4 = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/L4/Main_Fig4"
BAND = f"{R54}/lanes/LSKETCH/work/build/band_4e_r54.pdf"; BAND_R49 = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LSKETCH/work/build/band_4e_r49.pdf"
PAGE_W = 970.7899780273438            # V30's width, unchanged (the mediabox of V30's Main_Fig4)
RIGHT_EDGE = 956.621                  # the V13 right margin (panel b's q column ended here)
HEIGHT_GATE = 1260.0; HEIGHT_AIM = 1110.0
FS = dict(row=14.0, col=14.0, key=14.0, tick=14.0, head=15.0, block=14.0, letter=18.0, title=18.0, cell=14.0, note=13.0, clab=14.0)
SUB_PITCH_B = 15.0; PAIR_GAP_B = 2.5; PAIR_GAP_A = 3.5       # panel b's rows set the block, panel a's sub-row pitch is derived (fills the block)
KEY_PITCH = 17.5                      # key and legend line pitch at 14 pt (1.25 em)
MARK_DY = 0.35                        # marker centre above the sub-row baseline, in em (V13: 3.45 pt at 10 pt)
TICK_LEN = 3.0; TICK_GAP = 3.0; CAP = 0.716; DSC = 0.212
TITLE_BASE = 18.0; LETTER_BASE_AB = 40.0; LETTER_X_A = 14.17; LAB_X_A = 14.17
PANEL_GAP = 8.0; LAB_PLOT_GAP = 5.0; PLOT_COL_GAP = 7.0; COL_Q_GAP = 7.0; TICK_LABEL_GAP = 3.0
LEGEND_HANDLE_PT = 15.0; LEGEND_PAD_PT = 5.0   # the legend and key handles are drawn by hand (line + marker), 15 pt long, 5 pt before the text
BAND_GAP_UNDER_LETTER = 3.12          # V13: letter e baseline 1080.26, band top 1083.38


def measure():
    matplotlib, plt = L.mpl_setup()
    with plt.rc_context(L.RC):
        fig = plt.figure(figsize=(12, 12)); fig.canvas.draw(); r = fig.canvas.get_renderer()
        def tw(s, size, bold=False):
            t = fig.text(0.05, 0.5, s, fontsize=size, fontweight="bold" if bold else "normal"); w = t.get_window_extent(renderer=r).width / fig.dpi * 72; t.remove(); return float(w)
        return tw


def best_wrap(tw, label, size):
    parts = label.split(" "); best = None
    for k in range(1, len(parts)):
        l1, l2 = " ".join(parts[:k]), " ".join(parts[k:]); wmax = max(tw(l1, size), tw(l2, size))
        if best is None or wmax < best[0]: best = (wmax, l1, l2)
    return best


def tick_min_width(tw, ticks, lo, hi, size):
    need = 0.0
    for t0, t1 in zip(ticks[:-1], ticks[1:]):
        frac = (math.log(t1) - math.log(t0)) / (math.log(hi) - math.log(lo))
        need = max(need, ((tw(f"{t0:g}", size) + tw(f"{t1:g}", size)) / 2 + TICK_LABEL_GAP) / frac)
    return need


def forest_columns(tw, labels, hr_strings, q_strings, q_bold, ticks, lo, hi, wrap_pref):
    wraps = {}
    for l in wrap_pref:
        if l in labels: wraps[l] = best_wrap(tw, l, FS["row"])[1:]
    label_w = max(max(tw(ln, FS["row"]) for ln in (wraps[l] if l in wraps else (l,))) for l in labels)
    hrw = max(tw(s, FS["col"]) for s in hr_strings); qw = max(max(tw(s, FS["col"], b) for s, b in zip(q_strings, q_bold)), tw("q", FS["head"]))
    return dict(label_w=label_w, wraps=wraps, hr_w=hrw, q_w=qw, plot_min=tick_min_width(tw, ticks, lo, hi, FS["tick"]))


def rows_block(labels, n_main, y0, pitch, gap):
    """Sub-row baselines of a paired forest: wake (upper) and sleep (lower) per row, the block header line before the controls."""
    rows = []; y = y0
    for i, lab in enumerate(labels):
        if i == n_main: y += pitch + gap
        rows.append(dict(label=lab, wake=y, sleep=y + pitch, centre_base=y + pitch / 2)); y += 2 * pitch + gap
    return rows, rows[n_main]["wake"] - pitch


def gs_box(pdf):
    r = subprocess.run([paths.GS, "-q", "-dNODISPLAY", "-dNOSAFER", "-dBATCH", "-dNOPAUSE", "-dPDFINFO", pdf], capture_output=True, text=True)
    txt = r.stdout + r.stderr; m = re.search(r"MediaBox:\s*\[([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\]", txt); x0, y0, x1, y1 = (float(v) for v in m.groups())
    return len(re.findall(r"MediaBox", txt)), round(x1 - x0, 3), round(y1 - y0, 3), [x0, y0, x1, y1]


def build_spec():
    tw = measure()
    DA = json.load(open(f"{R49L4}/verify/Main_Fig4_a_drawn.json")); DB = json.load(open(f"{R49L4}/verify/Main_Fig4_b_drawn.json")); DD = json.load(open(f"{R49L4}/verify/Main_Fig4_d_drawn.json"))
    GA = DA["geometry"]; GB = DB["geometry"]
    labels_a = DA["rows_main"] + DA["rows_controls"]; labels_b = DB["rows_conditions"] + DB["rows_controls"]
    hr_a = [v["hr_ci_text"] for v in DA["values"]]; q_a = [v["q_text"] for v in DA["values"]]; qb_a = [v["sig"] for v in DA["values"]]
    hr_b = [v["hr_ci_text"] for v in DB["values"]]; q_b = [v["q_text"] for v in DB["values"]]; qb_b = [v["sig"] for v in DB["values"]]
    xlim_a = GA["xlim"]; ticks_a = GA["ticks"]; xlim_b = GB["xlim"]; ticks_b = GB["ticks"]
    # ---- widths: wrap the widest labels (the round-49 wrap of 'Ventricular arrhythmia or cardiac arrest' kept) until a | b fits the width
    fixed = LAB_X_A + LAB_PLOT_GAP + PLOT_COL_GAP + COL_Q_GAP + PANEL_GAP + LAB_PLOT_GAP + PLOT_COL_GAP + COL_Q_GAP
    order_a = sorted(labels_a, key=lambda l: -tw(l, FS["row"])); order_b = sorted([l for l in labels_b if l != "Ventricular arrhythmia or cardiac arrest"], key=lambda l: -tw(l, FS["row"]))
    wrap_a, wrap_b = [], ["Ventricular arrhythmia or cardiac arrest"]
    while True:
        ca = forest_columns(tw, labels_a, hr_a, q_a, qb_a, ticks_a, *xlim_a, wrap_a); cb = forest_columns(tw, labels_b, hr_b, q_b, qb_b, ticks_b, *xlim_b, wrap_b)
        need = fixed + ca["label_w"] + ca["plot_min"] + ca["hr_w"] + ca["q_w"] + cb["label_w"] + cb["plot_min"] + cb["hr_w"] + cb["q_w"]
        if need <= RIGHT_EDGE: break
        cand_a = [l for l in order_a if l not in wrap_a]; cand_b = [l for l in order_b if l not in wrap_b]
        if not cand_a and not cand_b: raise AssertionError(("a | b does not fit the width even with every label wrapped", need))
        if cand_a and (not cand_b or ca["label_w"] >= cb["label_w"]): wrap_a.append(cand_a[0])
        else: wrap_b.append(cand_b[0])
    spare = RIGHT_EDGE - need
    plot_a = ca["plot_min"] + spare * ca["plot_min"] / (ca["plot_min"] + cb["plot_min"]); plot_b = cb["plot_min"] + spare - (plot_a - ca["plot_min"])
    # ---- panel b first (its 42 sub-rows set the block)
    b = {}
    b["lab_x"] = LAB_X_A + ca["label_w"] + LAB_PLOT_GAP + plot_a + PLOT_COL_GAP + ca["hr_w"] + COL_Q_GAP + ca["q_w"] + PANEL_GAP
    b["plot_x0"] = b["lab_x"] + cb["label_w"] + LAB_PLOT_GAP; b["plot_x1"] = b["plot_x0"] + plot_b
    b["col_r"] = b["plot_x1"] + PLOT_COL_GAP + cb["hr_w"]; b["qcol_r"] = b["col_r"] + COL_Q_GAP + cb["q_w"]
    assert abs(b["qcol_r"] - RIGHT_EDGE) < 0.01, (b["qcol_r"], RIGHT_EDGE)
    b["legend_base"] = [LETTER_BASE_AB + 18.0 + k * KEY_PITCH for k in range(3)]; b["head_base"] = b["legend_base"][2]; b["y0"] = b["head_base"] + 17.0
    b["legend_x"] = b["lab_x"]; b["legend_end"] = b["legend_x"] + LEGEND_HANDLE_PT + LEGEND_PAD_PT + tw(DB["header_left"], FS["key"]); b["legend_fits"] = bool(b["legend_end"] <= PAGE_W - 1.5)
    assert b["legend_fits"], ("the legend's first line leaves the page", b["legend_end"])
    n_main_b = len(DB["rows_conditions"]); rows_b, b["block_base"] = rows_block(labels_b, n_main_b, b["y0"], SUB_PITCH_B, PAIR_GAP_B)
    for r in rows_b: r["lines"] = list(cb["wraps"].get(r["label"], (r["label"],)))
    b["rows"] = rows_b; b["sub_pitch"] = SUB_PITCH_B; b["pair_gap"] = PAIR_GAP_B; b["last_base"] = rows_b[-1]["sleep"]
    b["rule_y"] = b["last_base"] + DSC * FS["col"] + 5.0; b["tick_base"] = b["rule_y"] + TICK_LEN + TICK_GAP + CAP * FS["tick"]
    b["title_base"] = [b["tick_base"] + 18.0, b["tick_base"] + 36.0]; b["bottom"] = b["title_base"][1] + DSC * FS["head"]
    b["wraps"] = cb["wraps"]; b["widths"] = dict(label=cb["label_w"], plot=plot_b, plot_min=cb["plot_min"], hr=cb["hr_w"], q=cb["q_w"])
    # ---- panel a: same rule line, its sub-row pitch fills the block
    a = {}
    a["lab_x"] = LAB_X_A; a["plot_x0"] = LAB_X_A + ca["label_w"] + LAB_PLOT_GAP; a["plot_x1"] = a["plot_x0"] + plot_a
    a["col_r"] = a["plot_x1"] + PLOT_COL_GAP + ca["hr_w"]; a["qcol_r"] = a["col_r"] + COL_Q_GAP + ca["q_w"]
    assert abs(a["qcol_r"] + PANEL_GAP - b["lab_x"]) < 1e-9
    a["key_base"] = [LETTER_BASE_AB + 18.0, LETTER_BASE_AB + 18.0 + KEY_PITCH]; a["head_base"] = a["key_base"][1]; a["y0"] = a["head_base"] + 17.0
    n_main_a = len(DA["rows_main"]); n_a = len(labels_a)
    span = b["last_base"] - a["y0"]                       # 36 p + 18 g = span (17 pairs of (2p + g), one p, the block line p + g)
    pitch_a = (span - (n_a - 1 + 1) * PAIR_GAP_A) / (2 * n_a)
    assert pitch_a >= 14.5, ("panel a's derived pitch is below the 14-pt floor", pitch_a)
    rows_a, a["block_base"] = rows_block(labels_a, n_main_a, a["y0"], pitch_a, PAIR_GAP_A)
    for r in rows_a: r["lines"] = list(ca["wraps"].get(r["label"], (r["label"],)))
    a["rows"] = rows_a; a["sub_pitch"] = pitch_a; a["pair_gap"] = PAIR_GAP_A; a["last_base"] = rows_a[-1]["sleep"]
    assert abs(a["last_base"] - b["last_base"]) < 1e-6, (a["last_base"], b["last_base"])
    a["rule_y"] = b["rule_y"]; a["tick_base"] = b["tick_base"]; a["cap_base"] = list(b["title_base"]); a["bottom"] = a["cap_base"][1] + DSC * FS["head"]
    a["wraps"] = ca["wraps"]; a["widths"] = dict(label=ca["label_w"], plot=plot_a, plot_min=ca["plot_min"], hr=ca["hr_w"], q=ca["q_w"])
    AB_BOTTOM = max(a["bottom"], b["bottom"])
    a["box"] = [0.0, 0.0, a["qcol_r"] + PANEL_GAP / 2, AB_BOTTOM + 2.0]; b["box"] = [a["box"][2], 0.0, PAGE_W, AB_BOTTOM + 2.0]
    # ---- the sketch band = panel c under a | b (the coordinator's band_4e_r54.pdf, its page box read here)
    n, bw, bh, bbox = gs_box(BAND); assert n == 1 and abs(bw - PAGE_W) < 0.05, (n, bw, PAGE_W)
    e = dict(letter_top=AB_BOTTOM + 12.0); e["letter_base"] = e["letter_top"] + CAP * FS["letter"]; e["band_top"] = e["letter_base"] + BAND_GAP_UNDER_LETTER   # the letter's cap top 12 pt under the a/b ink (its ascender box 8.6 pt, the V13 gap was 8)
    e["band"] = BAND; e["band_w"] = bw; e["band_h"] = bh; e["band_mediabox"] = bbox; e["page_h"] = e["band_top"] + bh; e["letter"] = "c"
    assert e["page_h"] <= HEIGHT_AIM, ("the composed sheet is taller than the aim", e["page_h"])
    # ---- the counts panel (V30's panel c), a standalone page (no letter, no title strip): its own frame, origin (0, 0)
    c = {}
    c_row_labels = ["No or mild", "(AHI <15)", "Moderate", "(AHI 15 to <30)", "Severe", "(AHI 30 or more)"]
    c["lab_x1"] = LAB_X_A + max(tw(s, FS["clab"]) for s in c_row_labels) + 2.0; c["cell_w"] = 96.0; c["cell_gap"] = 3.6
    c["col_x"] = [c["lab_x1"] + 10.0 + i * (c["cell_w"] + c["cell_gap"]) for i in range(4)]; c["grid_x1"] = c["col_x"][3] + c["cell_w"]
    c["title_base"] = 8.0 + CAP * FS["head"]; c["header_y"] = [c["title_base"] + 6.0, c["title_base"] + 28.0]; c["row_h"] = 44.0; c["row_gap"] = 6.5
    c["row_top0"] = c["header_y"][1] + 14.0; c["grid_bottom"] = c["row_top0"] + 3 * c["row_h"] + 2 * c["row_gap"]
    c["key_base"] = [c["grid_bottom"] + 26.0, c["grid_bottom"] + 42.0]; c["bottom"] = c["key_base"][1] + DSC * FS["key"]
    c["box"] = [0.0, 0.0, c["grid_x1"] + 14.0, c["bottom"] + 8.0]
    # ---- V30's panel d, built for the record only (not delivered): its own frame
    d = {}
    d_row_labels = ["Cardiovascular", "composite", "Type 2 diabetes", "Acute kidney", "injury", "Heart failure", "Alopecia", "(negative control)", "Sleep apnea", "Oxygenation"]
    d["lab_x1"] = LAB_X_A + max(tw(s, FS["clab"]) for s in d_row_labels) + 2.0
    ci_w = max(tw(v["ci_text"], FS["cell"]) for v in DD["values"] if v["ci_text"]); hr_w = max(tw(v["hr_text"], FS["cell"]) for v in DD["values"])
    d["cell_w"] = max(ci_w, hr_w) + 2.4; d["cell_gap"] = 1.753; d["group_gap"] = 11.403
    d["grid_x0"] = d["lab_x1"] + 9.0; d["group_w"] = 2 * d["cell_w"] + d["cell_gap"]; d["group_pitch"] = d["group_w"] + d["group_gap"]; d["cell_pitch"] = d["cell_w"] + d["cell_gap"]
    d["grid_x1"] = d["grid_x0"] + 2 * d["group_pitch"] + d["group_w"]
    d["header1_y"] = [8.0, 42.0]; d["header2_y"] = [43.0, 77.0]   # two 14-pt lines (cap 10 + pitch 15 + descender 3 = 28) in 34-pt headers
    d["cell_h"] = 40.0; d["row_gap"] = 2.0; d["row0_y"] = d["header2_y"][1] + 3.0; d["row_pitch"] = d["cell_h"] + d["row_gap"]
    d["grid_bottom"] = d["row0_y"] + 4 * d["row_pitch"] + d["cell_h"]
    d["ramp"] = [d["grid_x0"] + 25.0, d["grid_bottom"] + 34.0, d["grid_x0"] + 25.0 + 124.593, d["grid_bottom"] + 44.0]
    d["ramp_title_base"] = d["grid_bottom"] + 30.0; d["ramp_tick_base"] = d["grid_bottom"] + 58.0
    d["refbox"] = [d["ramp"][2] + 30.0, d["ramp_title_base"] - 9.5, 14.5, 9.91]; d["ref_key_x"] = d["refbox"][0] + 14.5 + 5.0
    d["ref_key_base"] = d["ramp_title_base"]; d["star_key_base"] = d["ramp_title_base"] + 16.0
    d["bottom"] = d["ramp_tick_base"]; d["box"] = [0.0, 0.0, d["grid_x1"] + 14.0, d["bottom"] + 8.0]
    spec = dict(page_w=PAGE_W, right_edge=RIGHT_EDGE, sizes=FS, key_pitch=KEY_PITCH, mark_dy_em=MARK_DY, xlim_a=xlim_a, ticks_a=ticks_a, xlim_b=xlim_b, ticks_b=ticks_b,
                legend_handle_pt=LEGEND_HANDLE_PT, legend_pad_pt=LEGEND_PAD_PT, tick_label_gap=TICK_LABEL_GAP, gaps=dict(lab_plot=LAB_PLOT_GAP, plot_col=PLOT_COL_GAP, col_q=COL_Q_GAP, panel=PANEL_GAP),
                title=dict(text="Figure 4", x=LETTER_X_A, base=TITLE_BASE), letters=dict(a=dict(x=LETTER_X_A, base=LETTER_BASE_AB), b=dict(x=b["lab_x"], base=LETTER_BASE_AB), c=dict(x=LETTER_X_A, base=e["letter_base"])),
                a=a, b=b, c=c, d=d, e=e, ab_bottom=AB_BOTTOM, page_h=e["page_h"], height_gate=HEIGHT_GATE, height_aim=HEIGHT_AIM,
                arrangement="Main_Fig4 = a | b single-column paired forests side by side (b's 42 sub-rows at 15 pt set the block, a's 36 sub-rows fill it so the two axis rules share one line) + the sketch band as panel c; the counts panel standalone; V30's panel d record only (owner's option A)")
    return spec


def load(): return json.load(open(f"{WORK}/layout_r54.json"))


if __name__ == "__main__":
    os.makedirs(WORK, exist_ok=True); S = build_spec(); json.dump(S, open(f"{WORK}/layout_r54.json", "w"), indent=1, ensure_ascii=False)
    a, b, c, d, e = S["a"], S["b"], S["c"], S["d"], S["e"]
    print(f"panel a: label col {a['widths']['label']:.1f} (wraps {list(a['wraps'])}), plot {a['widths']['plot']:.1f} (min {a['widths']['plot_min']:.1f}), HR {a['widths']['hr']:.1f}, q {a['widths']['q']:.1f}; plot x {a['plot_x0']:.1f}..{a['plot_x1']:.1f}, cols {a['col_r']:.1f} / {a['qcol_r']:.1f}; sub-row pitch {a['sub_pitch']:.3f}, pair gap {a['pair_gap']}; rows {a['y0']:.1f}..{a['last_base']:.1f}, rule {a['rule_y']:.1f}, bottom {a['bottom']:.1f}")
    print(f"panel b: label col {b['widths']['label']:.1f} (wraps {list(b['wraps'])}), plot {b['widths']['plot']:.1f} (min {b['widths']['plot_min']:.1f}), HR {b['widths']['hr']:.1f}, q {b['widths']['q']:.1f}; lab x {b['lab_x']:.1f}, plot x {b['plot_x0']:.1f}..{b['plot_x1']:.1f}, cols {b['col_r']:.1f} / {b['qcol_r']:.1f}; sub-row pitch {b['sub_pitch']}, pair gap {b['pair_gap']}; rows {b['y0']:.1f}..{b['last_base']:.1f}, rule {b['rule_y']:.1f}, bottom {b['bottom']:.1f}; legend line 1 ends {b['legend_end']:.1f}")
    print(f"band (panel c): {e['band_w']} x {e['band_h']} (mediabox {e['band_mediabox']}) at top {e['band_top']:.2f}, letter c baseline {e['letter_base']:.2f} -> PAGE {S['page_w']:.3f} x {e['page_h']:.3f} (aim {HEIGHT_AIM}, gate {HEIGHT_GATE})")
    print(f"counts panel standalone: page {c['box'][2]:.1f} x {c['box'][3]:.1f}, grid x {c['col_x'][0]:.1f}..{c['grid_x1']:.1f}, y {c['header_y'][0]:.1f}..{c['grid_bottom']:.1f}; V30 panel d (record only): page {d['box'][2]:.1f} x {d['box'][3]:.1f}")
    print("wrote", f"{WORK}/layout_r54.json")
