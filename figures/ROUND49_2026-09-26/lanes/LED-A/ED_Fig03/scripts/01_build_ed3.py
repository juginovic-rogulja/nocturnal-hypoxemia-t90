#!/usr/bin/env python3
"""ROUND 49 lane LED-A copy (2026-09-26): DELTA = +1 pt inside txt() and tw(), the key sentence not drawn (round 40), no PyMuPDF render, output to
work/ED_Fig03_built.pdf (stamp_title_r49.py adds the strip). Original: ROUND37_2026-09-14/figures/lanes/V14_L2_T90/ED_Fig03/scripts/01_build_ed3.py.
ED_Fig03, V14 lane L2 (round 37), step 01 (build + compose): Extended Data Fig. 3 re-plotted from the v8.1 numbers.
Repointed copy of the round-30 builder (byte copy beside as 01_build_ed3_PRE_V8_1.py). Design and geometry: byte for byte the
round-27 lane C builder that made the V13 sheet (every drawing constant below is that builder's; text set as simple TrueType Arial
Regular and Arial Bold page fonts, the V13 sheet's own mechanism, so the text layer reads "0-1" and the fonts carry the V13 names).

What is NEW under v8.1, all data driven:
  row set     the August generator's rule APPLIED (closes GAP-11 for this sheet): panel a = the 15 largest top-band (>10%) hazard
              ratios among graded conditions that are not negative controls and not circular, descending; panel b = the five v8.1
              controls (numbers/negcontrols_v4_panel.csv, step 251) descending by the top-band hazard ratio. The positive control
              (00_probe_ed3.py) shows the same rule on the v7 snapshot reproduces the V13 sheet, so the rule IS the sheet's rule.
  values      hazard ratio and 95% CI per band from numbers/results_v2.json ["graded"]["outcomes"] (step 112) for all 20 rows
              (every v8.1 control is in the graded block; graded_newcontrols_v1.json is not read).
  y limits    bottom min(0.72, 0.9 x lowest CI), top max(1.05, 1.22 x highest CI) (the August rule)
  y ticks     MaxNLocator(nbins=4, steps=[1, 2, 2.5, 5, 10]) on the limits, ticks at or above bottom + 6% of span
  asterisks   *** where the Benjamini-Hochberg q of the top-band (>10%) contrast is below 0.001, family = every top-band contrast of
              results_v2.json graded (58 under v8.1, the round-27 rule with the file's own count); q across the non-controls, the
              15 drawn and the 5 controls are also computed and written.
  titles      wrapped by the generator's greedy rule at AXW + GAP - 5.76 pt (data driven; the V13 breaks are reproduced by it).
  key line 2  "*** top band, q below 0.001" (the V13 wording) printed only if every panel-a outcome has q < 0.001 and no control does
              (else the build stops and reports; --allow-partial draws the stars per panel without the sentence).
Usage: 01_build_ed3.py [--allow-partial]"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import fitz, gc, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ed3_common import *
DELTA = float(os.environ.get("ED3_DELTA", "1.0"))    # round 49: +1 pt on every text element, inside txt() and tw() (placement rules see the larger text)
DRAW_KEY2 = False                                     # round 40 (LNOTES): the key sentence "*** top band, q below 0.001" left for the legend
KEEP_SIZE = {"Venous thromboembolism", "Contact dermatitis"}   # round 49 exception: last-column titles that would leave the page at 10.5 pt (their lines would end at 521.5 and 522.7 pt on the 518.74 pt page): kept at 9.5 pt
KEY_MEAN_DX_R49 = 36.7 + 2.0                          # round 49 nudge: the band label "5-10" grows 2.0 pt at 10 pt, the wording moves right by that much (the V26 gap of 3.5 pt kept)

ARIAL = paths.FONT_ARIAL; ARIAL_B = paths.FONT_ARIAL_BOLD
INK = (0x1a / 255, 0x1d / 255, 0x21 / 255); WHITE = (1, 1, 1)
LADDER = [(0xb3 / 255, 0xdc / 255, 0xf2 / 255), (0x7c / 255, 0xc0 / 255, 0xe9 / 255), (0x3f / 255, 0x9f / 255, 0xd8 / 255), (0x02 / 255, 0x88 / 255, 0xd1 / 255)]
KEY_WORDING = ["1% or less", "over 1 to 5%", "over 5 to 10%", "over 10%"]
STAR, STAR_TAIL = "***", "top band, q below 0.001"
Y_LABEL = "Hazard ratio vs 0-1% (95% CI)"; X_TITLE = "Time with oxygen saturation below 90%, % of the recording"; GROUP_B = "Negative controls"
# ---- Figure 2c at authored scale (round 27, unchanged) -----------------------------------------
FIG2_AXW, FIG2_AXH, FIG2_GAP = 143.778 / 1.2, 96.768 / 1.2, 43.370 / 1.2
TITLE_DY = 4.0; XTICK_DY = 15.90 / 1.2; YTICK_DX = 7.32 / 1.2; YTICK_BASE = 0.377
XTITLE_DY = (1078.42 - 1048.65) / 1.2; KEY_DY = (1110.01 - 1078.42) / 1.2; BOTTOM = (1127.0 - 1110.01) / 1.2
KEY_ROW = 12.0; KEY_GAP = 26.0; SW_W, SW_H, SW_UP = 11.0, 8.0, 7.3; KEY_TICK_DX, KEY_MEAN_DX = 15.2, 36.7
TICK_LEN, TICK_W, SPINE_W = 2.6, 0.8, 0.8; REF_W, REF_DASH = 0.9, "[3.6 2.7] 0"; LINE_W, MARK_R, RING_W = 1.5, 2.0, 0.5; RIBBON_ALPHA = 0.16
# ---- this page (round 27, unchanged) ----------------------------------------------------------------
LEFT0, RIGHT, AXW = 45.0, 510.5, 72.3; XLAB_MIN_GAP = 2.8; GAP = (RIGHT - LEFT0 - 5 * AXW) / 4.0; PITCH = AXW + GAP
T1, ROW3_BOTTOM = 57.6, 355.0; LETTER_A, LETTER_B = (24.48, 24.8), (24.48, 388.39); YLAB_X, NC_DY, LINE_PITCH = 16.0, 20.0, 10.0
AXH = AXW * FIG2_AXH / FIG2_AXW; WRAP_LIMIT = AXW + GAP - 0.08 * 72.0
V13_BREAKS = {"Obesity hypoventilation": ["Obesity", "hypoventilation"], "Pulmonary hypertension": ["Pulmonary", "hypertension"], "Venous thromboembolism": ["Venous", "thromboembolism"]}


def load(allow_partial):
    prov = {DATA: sidecar_v8(DATA), PANEL_CSV: sidecar_v8(PANEL_CSV)}
    G = graded(DATA); out = G["outcomes"]; controls = controls_of(PANEL_CSV)
    assert len(controls) == 5 and all(c in out and out[c]["negative_control"] for c in controls), controls
    assert sorted(k for k, v in out.items() if v["negative_control"]) == sorted(controls), "the graded block's control flags are not the panel's five"
    for c in out: assert abs(out[c]["0-1%"]["hr"] - 1.0) < 1e-9, c
    base = json.load(open(hydrated(f"{WORK}/base_probe.json"))); base_rows = [p["title"] for p in base["panels"]]
    rule_a, rule_b = august_rows(out, controls); rows = rule_a + rule_b
    assert len(rows) == 20 and len(set(rows)) == 20
    qall = bh_top(out); q_nc = bh_top(out, [k for k, v in out.items() if not v["negative_control"]]); q15 = bh_top(out, rule_a); q5 = bh_top(out, rule_b)
    recs = []
    for k in rows:
        v = out[k]; hr = [1.0] + [v[b]["hr"] for b in BAND_KEYS[1:]]; lo = [None] + [v[b]["lo"] for b in BAND_KEYS[1:]]; hi = [None] + [v[b]["hi"] for b in BAND_KEYS[1:]]
        ymin, ymax, ticks = ylim_ticks(v); assert len(ticks) >= 2, (k, ticks, ymin, ymax)
        q = qall[k]
        recs.append(dict(title=k, src=DATA, keyroot=f"graded.outcomes['{k}']", hr=hr, lo=lo, hi=hi, ymin=ymin, ymax=ymax, yticks=ticks, star=STAR if q < 0.001 else "",
                         negative=bool(v["negative_control"]), events=int(v["events"]), p_by_band={b: v[b]["p"] for b in BAND_KEYS[1:]}, p_top=v[">10%"]["p"],
                         q_all=q, q_nc=q_nc.get(k), q15=q15.get(k), q5=q5.get(k)))
    all15 = all(r["q_all"] < 0.001 for r in recs[:15]); ctrl_none = all(r["q_all"] >= 0.001 for r in recs[15:])
    if not (all15 and ctrl_none) and not allow_partial:
        bad = [(r["title"], r["p_top"], r["q_all"]) for r in recs if (r["q_all"] < 0.001) != (not r["negative"])]
        raise SystemExit(f"STAR RULE CHANGES THE SHEET: {bad}. Rerun with --allow-partial after reading the report rules.")
    meta = dict(family_n=len(out), q_all=qall, q_nc=q_nc, q15=q15, q5=q5, all15_q_below_0_001=all15, controls_none_below_0_001=ctrl_none,
                rule="August rule applied (GAP-11 closed for this sheet): panel a = 15 largest top-band hazard ratios among non-control, non-circular graded conditions; panel b = the five controls, both descending",
                rows_drawn=rows, rows_v13=base_rows, entering=[r for r in rows if r not in base_rows], leaving=[r for r in base_rows if r not in rows], controls=controls)
    return G, recs, meta, prov


def build(out_pdf, drawn_path, png_path, allow_partial=False):
    G, recs, meta, prov = load(allow_partial)
    base = fitz.open(hydrated(BASE)); W, H = base[0].rect.width, base[0].rect.height; base.close()
    fR, fB = fitz.Font(fontfile=ARIAL), fitz.Font(fontfile=ARIAL_B)
    doc = fitz.open(); pg = doc.new_page(width=W, height=H)
    pg.insert_font(fontname="ArialR", fontfile=ARIAL, set_simple=True); pg.insert_font(fontname="ArialB", fontfile=ARIAL_B, set_simple=True)
    drawn = []
    def txt(x, y, s, size, bold=False, rotate=0, what=None, value=None, key=None, src=None, panel=None):
        size = size + (0.0 if (what == "title" and key in KEEP_SIZE) else DELTA)             # round 49 (KEEP_SIZE titles stay at 9.5 pt)
        pg.insert_text(fitz.Point(x, y), s, fontsize=size, fontname="ArialB" if bold else "ArialR", color=INK, rotate=rotate)
        drawn.append(dict(text=s, x=round(x, 3), baseline=round(y, 3), size=size, bold=bold, rotate=rotate, what=what, value=value, key=key, source=src, panel=panel))
        return (fB if bold else fR).text_length(s, size)
    def tw(s, size=9.5, bold=False): return (fB if bold else fR).text_length(s, size + DELTA)   # round 49: widths at the drawn size
    P = (ROW3_BOTTOM - T1 - AXH) / 2.0; row_tops = [T1 + i * P for i in range(3)]
    below = XTICK_DY + XTITLE_DY + KEY_DY + KEY_ROW + BOTTOM; b_bottom = H - below; b_top = b_bottom - AXH
    b_title_base = b_top - TITLE_DY; nc_base = b_title_base - NC_DY; xtick_b = b_bottom + XTICK_DY; xtitle_base = xtick_b + XTITLE_DY
    key1_base = xtitle_base + KEY_DY; key2_base = key1_base + KEY_ROW
    assert nc_base > LETTER_B[1] + 6.0 and row_tops[2] + AXH + XTICK_DY < LETTER_B[1] - 0.716 * (13 + DELTA) - 6.0
    for i, r in enumerate(recs):
        if i < 15: row, col = divmod(i, 5); top = row_tops[row]
        else: col = i - 15; top = b_top
        left = LEFT0 + col * PITCH; r["box"] = (left, top, left + AXW, top + AXH); r["panel"] = "a" if i < 15 else "b"; r["col"] = col
    marks = []; YTICK_NUDGES = []
    for r in recs:
        L, T, R, B = r["box"]; ymin, ymax = r["ymin"], r["ymax"]; bw = AXW / 4.0; xc = [L + bw * (k + 0.5) for k in range(4)]
        y_of = lambda v: T + (ymax - v) / (ymax - ymin) * AXH
        sh = pg.new_shape()
        for k in range(4): sh.draw_rect(fitz.Rect(L + k * bw, T, L + (k + 1) * bw, B)); sh.finish(color=None, fill=LADDER[k])
        pts = [fitz.Point(xc[1], y_of(r["hi"][1])), fitz.Point(xc[1], y_of(r["lo"][1])), fitz.Point(xc[2], y_of(r["lo"][2])), fitz.Point(xc[3], y_of(r["lo"][3])), fitz.Point(xc[3], y_of(r["hi"][3])), fitz.Point(xc[2], y_of(r["hi"][2]))]
        sh.draw_polyline(pts); sh.finish(color=None, fill=INK, fill_opacity=RIBBON_ALPHA, closePath=True); sh.commit()
        yref = y_of(1.0); sh = pg.new_shape(); sh.draw_line(fitz.Point(L, yref), fitz.Point(R, yref)); sh.finish(color=INK, width=REF_W, dashes=REF_DASH, closePath=False); sh.commit()
        sh = pg.new_shape(); sh.draw_line(fitz.Point(L, T), fitz.Point(L, B)); sh.draw_line(fitz.Point(L, B), fitz.Point(R, B)); sh.finish(color=INK, width=SPINE_W, lineCap=2, closePath=False); sh.commit()
        sh = pg.new_shape()
        for x in xc: sh.draw_line(fitz.Point(x, B), fitz.Point(x, B + TICK_LEN))
        for lab in r["yticks"]: y = y_of(float(lab)); sh.draw_line(fitz.Point(L - TICK_LEN, y), fitz.Point(L, y))
        sh.finish(color=INK, width=TICK_W, lineCap=0, closePath=False); sh.commit()
        sh = pg.new_shape(); sh.draw_polyline([fitz.Point(xc[k], y_of(r["hr"][k])) for k in range(4)]); sh.finish(color=WHITE, width=LINE_W, lineCap=1, closePath=False); sh.commit()
        sh = pg.new_shape()
        for k in range(4): sh.draw_circle(fitz.Point(xc[k], y_of(r["hr"][k])), MARK_R); sh.finish(color=WHITE, fill=INK, width=RING_W)
        sh.commit()
        for k in range(4):
            marks.append(dict(condition=r["title"], band=BAND_KEYS[k], hr=r["hr"][k], lo=r["lo"][k], hi=r["hi"][k], x_page=round(xc[k], 3), y_page=round(y_of(r["hr"][k]), 3),
                              y_lo_page=None if r["lo"][k] is None else round(y_of(r["lo"][k]), 3), y_hi_page=None if r["hi"][k] is None else round(y_of(r["hi"][k]), 3),
                              source=r["src"], key=f"{r['keyroot']}['{BAND_KEYS[k]}']"))
        tsize = 9.5 if r["title"] in KEEP_SIZE else 9.5 + DELTA; pitch = LINE_PITCH + (tsize - 9.5)     # round 49: the wrap rule and the line pitch at the drawn size
        lines = wrap_title(r["title"], lambda s: fR.text_length(s, tsize), WRAP_LIMIT)
        if r["title"] in V13_BREAKS: assert lines == V13_BREAKS[r["title"]], (r["title"], lines)
        r["title_lines"] = lines; r["title_lines_v26"] = wrap_title(r["title"], lambda s: fR.text_length(s, 9.5), WRAP_LIMIT); r["title_size"] = tsize; r["title_pitch"] = pitch
        assert L + max(fR.text_length(l, tsize) for l in lines) < W - 2.0, ("title leaves the page", r["title"], tsize)
        for j, line in enumerate(lines):
            txt(L, T - TITLE_DY - (len(lines) - 1 - j) * pitch, line, 9.5, what="title", key=r["title"], value=f"line {j + 1} of {len(lines)}", src=r["src"], panel=r["panel"])
        wl = [tw(lab, 9.5) for lab in BANDS]; lefts = [xc[k] - wl[k] / 2 for k in range(4)]
        for k in range(1, 4): lefts[k] = max(lefts[k], lefts[k - 1] + wl[k - 1] + XLAB_MIN_GAP)
        drift = (lefts[3] + wl[3] / 2) - xc[3]; lefts = [l - drift / 2 for l in lefts]
        for k, lab in enumerate(BANDS): txt(lefts[k], B + XTICK_DY, lab, 9.5, what="xtick", key=f"band label {BAND_KEYS[k]}", value=lab, src=DATA, panel=r["panel"])
        for lab in r["yticks"]:
            x_lab = L - YTICK_DX - tw(lab, 9.5); nb = (L - GAP + 0.5) if r["col"] > 0 else -1e9     # round 49 nudge: a label clears the neighbouring panel's box by 0.5 pt
            if x_lab < nb: YTICK_NUDGES.append(dict(title=r["title"], label=lab, x_rule=round(x_lab, 3), x=round(nb, 3), shift=round(nb - x_lab, 3))); x_lab = nb
            txt(x_lab, y_of(float(lab)) + YTICK_BASE * (9.5 + DELTA), lab, 9.5, what="ytick", value=float(lab), key=f"tick rule on ({r['ymin']:.4f}, {r['ymax']:.4f}) from {r['keyroot']} lo/hi", src=r["src"], panel=r["panel"])
        if r["star"]:
            pref = T + 0.17 * AXH + 0.5 * (9.5 + DELTA); y_rib = y_of(r["hi"][3]); base_y = max(min(pref, y_rib - 1.5), T + 8.3)
            txt(xc[3] - tw(STAR, 9.5, True) / 2, base_y, STAR, 9.5, bold=True, what="star", value=r["q_all"], key=f"BH q ({meta['family_n']} top-band contrasts of graded.outcomes) of {r['keyroot']}['>10%'].p = {r['p_top']:.3g}", src=DATA, panel=r["panel"])
            r["star_base"] = base_y
    txt(*LETTER_A, "a", 13, bold=True, what="letter", panel="a"); txt(*LETTER_B, "b", 13, bold=True, what="letter", panel="b")
    cx = (LEFT0 + RIGHT) / 2.0
    txt(cx - tw(GROUP_B, 10, True) / 2, nc_base, GROUP_B, 10, bold=True, what="group", panel="b")
    txt(cx - tw(X_TITLE, 11) / 2, xtitle_base, X_TITLE, 11, what="xtitle")
    block_mid = (T1 + b_bottom) / 2.0
    txt(YLAB_X, block_mid + tw(Y_LABEL, 11) / 2, Y_LABEL, 11, rotate=90, what="ylabel")
    gw = [KEY_MEAN_DX_R49 + tw(m, 9) for m in KEY_WORDING]; total = sum(gw) + 3 * KEY_GAP; x = cx - total / 2.0
    for k in range(4):
        pg.draw_rect(fitz.Rect(x, key1_base - SW_UP, x + SW_W, key1_base - SW_UP + SW_H), color=None, fill=LADDER[k])
        txt(x + KEY_TICK_DX, key1_base, BANDS[k], 9, what="key", key=f"band label {BAND_KEYS[k]}", value=BANDS[k], src=DATA)
        txt(x + KEY_MEAN_DX_R49, key1_base, KEY_WORDING[k], 9, what="key", key=f"band wording {BAND_KEYS[k]}", value=KEY_WORDING[k], src=DATA)
        x += gw[k] + KEY_GAP
    ast_w = tw(STAR, 9.5, True) + 3.0 + tw(STAR_TAIL, 9); ax0 = cx - ast_w / 2.0
    if DRAW_KEY2:                                                                            # round 40: not drawn (the rule above still gates the build)
        txt(ax0, key2_base, STAR, 9.5, bold=True, what="key")
        txt(ax0 + tw(STAR, 9.5, True) + 3.0, key2_base, STAR_TAIL, 9, what="key", key="max BH q over the 15 panel-a outcomes", value=max(r["q_all"] for r in recs[:15]), src=DATA)
    doc.subset_fonts(); n_fixed = 0
    for xref in pg.get_contents():
        st = doc.xref_stream(xref); k = st.count(b"1 J\n")
        if k: doc.update_stream(xref, st.replace(b"1 J\n", b"1 j\n1 J\n")); n_fixed += k
    assert n_fixed == 20, n_fixed
    doc.save(out_pdf, garbage=3, deflate=True); doc.close()
    # round 49: renders by Ghostscript in run_sheet.py (no PyMuPDF pixmap)
    info = dict(sheet=out_pdf, sheet_sha256=sha256(out_pdf), page=[W, H], base_sheet=BASE, base_sheet_sha256=sha256(BASE), sources=prov, band_n=G["band_n"], band_n_printed_on_sheet=False,
                rules=dict(rows=meta["rule"], ylim="bottom min(0.72, 0.9 x lowest CI), top max(1.05, 1.22 x highest CI)",
                           yticks="MaxNLocator(nbins=4, steps=[1,2,2.5,5,10]) on the limits, ticks >= bottom + 0.06 x span, %g labels",
                           star=f"*** where Benjamini-Hochberg q of the top-band contrast < 0.001, q across all {meta['family_n']} top-band contrasts of results_v2.json graded (round-27 rule)",
                           key_line_2="printed only if all 15 panel-a outcomes have q < 0.001 and no control does", title_wrap=f"greedy measured rule at {WRAP_LIMIT:.2f} pt (the V13 breaks reproduced)",
                           text="simple TrueType Arial Regular/Bold WinAnsi page fonts, the V13 sheet's mechanism; letters likewise at the V13 origins"),
                meta=meta, axes=dict(w=AXW, h=AXH, gap=GAP, pitch=PITCH), row_tops=row_tops, row_pitch=P, b_top=b_top, b_bottom=b_bottom, nc_base=nc_base, xtick_b=xtick_b,
                xtitle_base=xtitle_base, key1_base=key1_base, key2_base=key2_base, letters=dict(a=LETTER_A, b=LETTER_B),
                delta_pt=DELTA, key2_drawn=DRAW_KEY2, keep_size=sorted(KEEP_SIZE), ytick_nudges=YTICK_NUDGES, key_mean_dx=KEY_MEAN_DX_R49, key_mean_dx_v26=KEY_MEAN_DX, line_pitch_v26=LINE_PITCH, wrap_changed_at_plus1=[dict(title=r["title"], v26=r["title_lines_v26"], new=r["title_lines"]) for r in recs if r["title_lines"] != r["title_lines_v26"]],
                panels=[dict(title=r["title"], panel=r["panel"], source=r["src"], keyroot=r["keyroot"], title_lines=r["title_lines"], title_lines_v26=r["title_lines_v26"], title_size=r["title_size"], title_pitch=r["title_pitch"], box=r["box"], ymin=r["ymin"], ymax=r["ymax"], yticks=r["yticks"],
                             hr=r["hr"], lo=r["lo"], hi=r["hi"], star=r["star"], star_base=r.get("star_base"), events=r["events"], p_by_band=r["p_by_band"], p_top=r["p_top"],
                             q_all=r["q_all"], q_nc=r["q_nc"], q15=r["q15"], q5=r["q5"], negative=r["negative"]) for r in recs],
                marks=marks, text=drawn)
    json.dump(info, open(drawn_path, "w"), indent=1)
    return info


if __name__ == "__main__":
    info = build(f"{WORK}/ED_Fig03_built.pdf", f"{VER}/ED_Fig03_drawn.json", None, allow_partial="--allow-partial" in sys.argv)   # round 49: the strip is stamped afterwards into LANE/ED_Fig03.pdf
    m = info["meta"]
    print("rows drawn a:", m["rows_drawn"][:15]); print("rows drawn b:", m["rows_drawn"][15:]); print("entering:", m["entering"], "leaving:", m["leaving"])
    print("stars:", "".join("*" if p["star"] else "-" for p in info["panels"]), "family n", m["family_n"])
    print("q max over 15:", max(p["q_all"] for p in info["panels"][:15]), "q min over controls:", min(p["q_all"] for p in info["panels"][15:]))
    for p in info["panels"]:
        print(f"  {p['title']:42s} {p['title_lines']} hr {[round(v, 3) for v in p['hr']]} ticks {p['yticks']} lim {p['ymin']:.4f}-{p['ymax']:.4f} p_top {p['p_top']:.3g} q {p['q_all']:.3g} star {p['star']!r}")
    print("wrote", info["sheet"])
