#!$T90_PY
"""Round 54 sketches (2026-09-29, R54_BRIEF.md: notes 15, labels 14, captions 14 pt), from ROUND49/lanes/LSKETCH/scripts/build_r49.py
(notes 12, captions 11) on the same engine (build_r44b on build_r40 on l1_common, FigureLab/bioglyph rendered by Chrome). The same
patched sources as round 49, only the sizes and what the sizes force:
  band_a_r54        tile labels 14 (11), captions 14 (11), 'years' 15 (12); the second caption on two lines (142.9 mm on one line at
                    14 pt bold, the page has 57 mm right of its centre), the band 2.0 mm taller for the two-line label under the oxygen trace
  band_1f_r54       every label 15 (12); the trace 5.5 mm further right so the 15-pt '90% SpO2' keeps its left margin
  band_3d_r54       notes 15, caption 14; fork arrowheads 1.5 mm clear of the 15-pt '90% SpO2' (HEAD_BACK from the measured width)
  band_4e_r54       notes 15 (ramp-end labels on two lines as in round 49), caption 14; the trace pair 7.4 mm further right so the fork
                    run from the ramp keeps 25 mm; the severity caption below the taller ramp labels and the band taller for it
  band_5d_r54       notes 15, caption 14; the trace pair on the same x as 4e
  Main_Fig6_tri_r54 notes 15, card captions 14 (10), triangle labels 14 (12); triangle gap 12 -> 8, card padding 8 -> 6 mm and the card
                    4.9 mm left so the wider card, labels and the 25-mm fork run fit the unchanged width (organ arrows 12 mm as before)
usage: build_r54.py [--final] [--only stem,stem]   (Chrome renders under the ROUND30 render lock, one process)"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import argparse, fcntl, inspect, json, math, sys, time
sys.dont_write_bytecode = True                                     # no __pycache__ in the earlier rounds' script folders
from dataclasses import replace
from pathlib import Path
FIGURELAB = paths.FIGURELAB_DIR  # $T90_FIGURELAB_DIR is a broken symlink (its target moved into Downloads/Older); read from the moved folder, nothing under ~/Downloads is written
sys.path.insert(0, FIGURELAB); import bioglyph                     # noqa: E402  (cached before l1_common's own broken path insert)
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND44_2026-09-21/lanes/LSKETCH/scripts")
import build_r44b as R                                             # noqa: E402
from build_r44b import *                                           # noqa: E402,F401,F403  (new, finish, trace, arrow, box, AZ, GREY, MM, TH_STD, B, ...)
from bioglyph.textmetrics import text_width_mm, wrap_text          # noqa: E402
R54 = Path(f"{paths.FIGURE_ROOT}/ROUND54_2026-09-28/lanes/LSKETCH"); OUT = R54 / "work" / "build"; WORK = R54 / "work"
R.B.L.PREVIEWS = WORK / "previews"                                 # finish() writes its preview PNGs here (build_r44b pointed it at ROUND44's folder)
PT = 15.0     # notes (round 49: 12)
LAB = 14.0    # labels and captions (round 49: 11)
TRI_PT = 14.0 # Figure 6 triangle labels (round 49: 12, the note size)
HEAD_BACK_R54 = B.LABEL90_GAP + text_width_mm(B.LABEL90, PT, False) + B.HEAD_CLEAR   # 28.29 mm: the fork arrowheads stop 1.5 mm before the 15-pt '90% SpO2' (round 49: 18.9 + 4.5)
TX_4E5D_R54 = 138.0   # the 4e and 5d trace pairs (round 32: 130.568): the fork run from 4e's source at x 84 keeps 25 mm (84 + 25 + 28.29 = 137.3)
TX_1F_R54 = 28.5      # band 1f trace (round 49: 23): the 15-pt '90% SpO2' (25.59 mm) ends 1.2 mm before the trace and starts 1.7 mm inside the page
CAP2_MAX_W = 75.0; CAP2_LINES = ["Which parameter best predicts", "54 future diseases and death?"]
H_MM_A_R54 = 81.5     # band a H_MM (round 49: 79.5): the two-line 14-pt label under the oxygen trace keeps the 2.3 mm bottom rule
assert wrap_text("Which parameter best predicts 54 future diseases and death?", LAB, CAP2_MAX_W, True) == CAP2_LINES

# ---------------------------------------------------------------- band a (round 47/49 source patches re-applied, plus the two-line second caption)
_src = inspect.getsource(R.band_a_r44)
_old = 'widest = max(L.text_width_mm(line, LABEL_PT, False) for line in tiles_txt["hyp"].split("\\n"))'
_new = 'widest = max(L.text_width_mm(line, LABEL_PT, False) for key in ("hyp", "eeg") for line in tiles_txt[key].split("\\n"))   # round 47: the eeg label is the widest left-column label'
assert _old in _src, "band_a_r44 source changed"
_src = _src.replace(_old, _new).replace('def band_a_r44(final, *, stem="band_a_r44", outdir=OUT, merged=True):', 'def band_a_r44_r54(final, *, stem="band_a_r54", outdir=OUT, merged=False):')
_o2 = 'W, H = A["W1_PT"] * MM, h_mm; fig = new(W, H); t = fig.theme'
assert _o2 in _src, "band_a theme line"
_src = _src.replace(_o2, 'W, H = A["W1_PT"] * MM, h_mm; fig = new(W, H); fig.theme = _replace(fig.theme, note_pt=15.0, label_pt=14.0); t = fig.theme   # round 54: notes 15, captions 14')
_o3 = 'lb = fig.label(A["CAP2"], c2, title_y, pt=t.label_pt, bold=True, color=INK21, anchor="middle", va="middle", leading=1.15)'
assert _o3 in _src, "band_a caption 2 line"
_src = _src.replace(_o3, 'lb = fig.label(A["CAP2"], c2, title_y - 0.925 * t.label_pt * MM / 2, pt=t.label_pt, bold=True, color=INK21, anchor="middle", va="top", leading=1.15, max_w=A["CAP2_MAX_W"]); assert lb.lines == 2, ("caption 2 on two lines", lb)   # round 54: two lines, the first where the one-line caption sat')
_o4 = 'assert c2 + one_line_w / 2 <= W - 1.0 and c2 - one_line_w / 2 >= caps[1]["label_box"][2] + 3.0, ("title 2 on one line", c2, one_line_w, caps)'
assert _o4 in _src, "band_a one-line assert"
_src = _src.replace(_o4, 'two_line_w = max(L.text_width_mm(s, t.label_pt, True) for s in A["CAP2_LINES"]); assert c2 + two_line_w / 2 <= W - 1.0 and c2 - two_line_w / 2 >= caps[1]["label_box"][2] + 3.0, ("title 2 on two lines", c2, two_line_w, caps)')
R.__dict__["_replace"] = replace
exec(_src, R.__dict__)   # band_a_r44_r54 in build_r44b's namespace
(WORK / "band_a_r44_r54_source.py").write_text(_src)


def band_a_r54(final, *, stem="band_a_r54", outdir=OUT):
    keep = {k: R.A[k] for k in ("LABEL_PT", "LABEL_GAP_BELOW", "LABEL_LEADING", "CAP1", "CAP2", "H_MM")}; added = ("CAP2_MAX_W", "CAP2_LINES")
    R.A["CAP1"] = "one night of sleep, 141 parameters"; R.A["CAP2"] = "Which parameter best predicts 54 future diseases and death?"   # the round-49 wording
    R.A["CAP2_MAX_W"] = CAP2_MAX_W; R.A["CAP2_LINES"] = list(CAP2_LINES)
    tiles0 = R.A_TILES; R.A_TILES = dict(dict(tiles0), air="breathing\n(airflow)", ecg="heart rhythm\n(ECG)")   # round 47: two lines like the other tile labels
    R.A["LABEL_PT"] = LAB; R.A["LABEL_GAP_BELOW"] = 0.7; R.A["LABEL_LEADING"] = 1.1; R.A["H_MM"] = H_MM_A_R54
    try: rec = R.band_a_r44_r54(final, stem=stem, outdir=outdir, merged=False)
    finally:
        R.A.update(keep); [R.A.pop(k, None) for k in added]; R.A_TILES = tiles0
    bottom_labels = max(t["label_box"][3] for t in rec["tiles"]); rec["bottom_margin_mm"] = round(rec["h_mm"] - bottom_labels, 3)
    rec["caption2_lines"] = list(CAP2_LINES); rec["h_mm_r49"] = 81.532
    rec["moved"] = (f"round 54: tile labels {keep['LABEL_PT']} -> {LAB} pt, captions {keep['LABEL_PT']} -> {LAB} pt, 'years' 12 -> {PT} pt, the second caption on two lines "
                    f"(max_w {CAP2_MAX_W} mm, first line where the one-line caption sat, va top), H_MM {keep['H_MM']} -> {H_MM_A_R54} (the band {H_MM_A_R54 - keep['H_MM']:.1f} mm taller for the "
                    f"two-line label under the oxygen trace, bottom margin {rec['bottom_margin_mm']} mm, rule 2.3), the scene starts {rec['scene_shift_mm']} mm right of the round-36 margin "
                    f"(the 14-pt 'brain waves (EEG)' keeps the margin, round-47 rule). Layout of band_a_r49 otherwise"); return rec


RANK_CLAIM = "most predictive of 141 parameters"


def band_1f_r54(final, *, stem="band_1f_r54", outdir=OUT):
    W = B.B22.W1_PT * MM; H = 46.0
    fig = R.new(W, H); fig.theme = replace(fig.theme, note_pt=PT, label_pt=PT); t = fig.theme; CY = H / 2   # every 1f label at the note size, as round 49 (12) did
    TX, TW, TH = TX_1F_R54, 67.0, TH_STD; TY = CY - TH / 2
    r1 = fig.label(RANK_CLAIM, TX + TW / 2, TY - 5.5, pt=t.note_pt, bold=True, color=AZ, anchor="middle", va="middle")
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
    assert tr["label_box"][0] >= 1.0, ("the '90% SpO2' label keeps a left margin", tr["label_box"])
    rec = finish(fig, stem, outdir, final=final, attribution=None)
    rec.update(page_pt=[round(W / MM, 2), round(H / MM, 2)], trace=tr, arrow=fa, label=lab, organs={"row_y": CY, "ink_x": ink_x, "names": [p.name for p in icons], "pitch": PITCH}, rank_claim=RANK_CLAIM, trace_h=TH,
               label_pt=PT, tx=TX, moved=f"round 54: labels {PT} pt (were 12), the trace at x {TX} (was 23.0) so the wider '90% SpO2' keeps its left margin ({tr['label_box'][0]:.2f} mm). Layout of band_1f_r49 otherwise"); return rec


# ---------------------------------------------------------------- Figure 6 triangles (round 47/49 patches re-applied on triangles_v3, the labels at TRI_PT)
_tsrc = inspect.getsource(R.triangles_v3)
_told = 't1l = fig.label("141 sleep measurements", cx1, top - 3.0, pt=t.note_pt, bold=True, color=INK21, anchor="middle", va="bottom")'
_tnew = 't1l = fig.label("141 sleep measurements", cx1, top - 3.0, pt=t.note_pt, bold=True, color=INK21, anchor="middle", va="bottom", max_w=w + 6.0, leading=1.1)   # round 47: two lines'
assert _told in _tsrc, "triangles_v3 source changed"
_tsrc = _tsrc.replace(_told, _tnew).replace("def triangles_v3(fig, *, x_left, top, h, w=28.0, gap=12.0, r=3.0):", "def triangles_v5(fig, *, x_left, top, h, w=28.0, gap=12.0, r=3.0):")
assert _tsrc.count('"141 sleep measurements"') == 1 and _tsrc.count('"all other\\nsleep\\nmeasurements"') == 1, "triangle label strings"
_tsrc = _tsrc.replace('"141 sleep measurements"', '"141 PSG parameters"').replace('"all other\\nsleep\\nmeasurements"', '"all other\\nPSG\\nparameters"')   # 2026-09-26: parameters
assert _tsrc.count("top - 3.0, pt=t.note_pt, bold=True") == 2, "title lines"
_tsrc = _tsrc.replace("top - 3.0, pt=t.note_pt, bold=True", "top - 4.5, pt=t.note_pt, bold=True")   # round 47: titles 4.5 mm above the triangles
assert _tsrc.count("pt=t.note_pt") == 7, ("seven triangle labels", _tsrc.count("pt=t.note_pt"))
_tsrc = _tsrc.replace("pt=t.note_pt", "pt=TRI_PT")   # round 54: the triangle labels at 14 pt, the notes ('90% SpO2') at 15
_t1 = 't1l = fig.label("141 PSG parameters", cx1, top - 4.5, pt=TRI_PT, bold=True, color=INK21, anchor="middle", va="bottom", max_w=w + 6.0, leading=1.1)'
assert _tsrc.count(_t1) == 1, "rank title call"
_tsrc = _tsrc.replace(_t1, 't1l = fig.label("141 PSG parameters", cx1, top - 4.5 - 1.1 * TRI_PT * MM, pt=TRI_PT, bold=True, color=INK21, anchor="middle", va="bottom", max_w=w + 6.0, leading=1.1); assert t1l.lines == 2, ("rank title on two lines", t1l)')
# round 54: bioglyph's label() knows va top and middle only, so 'bottom' is the FIRST line's baseline: with two lines the second line sat one leading below
# top - 4.5, on the blue cap (a latent round-49 placement, 0.16 mm below the triangle top at 12 pt, 0.93 mm at 14). The first baseline now sits one leading
# higher so the second line's baseline is at top - 4.5, level with the one-line 'disease risk'.
R.__dict__["TRI_PT"] = TRI_PT
exec(_tsrc, R.__dict__); (WORK / "triangles_v5_source.py").write_text(_tsrc)


CARD_PAD_R54 = 6.0   # the inputs card's padding (round 49: 8.0): the 14-pt captions widen the card by 9.4 mm, 4 of them taken back here
_inputs_card_r44 = R.inputs_card
R.inputs_card = lambda fig, cxi, mid, *, box_pad=CARD_PAD_R54: _inputs_card_r44(fig, cxi, mid, box_pad=box_pad)   # fig6_core looks the helper up at call time


def fig6_tri_r54(final, *, stem="Main_Fig6_tri_r54", outdir=OUT):
    fig = R.new(R.W6, R.H6_MM); fig.theme = replace(fig.theme, note_pt=PT, label_pt=LAB); mid = R.H6_MM / 2   # R.new is the themed new (the card's probe label must measure at 14 pt too)
    tri = R.triangles_v5(fig, x_left=32.5, top=mid - 23.0, h=46.0, gap=8.0)   # round 49: x_left 34, gap 12
    lb = tri["labels"]; gap_titles = lb["title_risk"][0] - lb["title_rank"][2]
    assert gap_titles >= 2.0, ("triangle titles too close", gap_titles, lb["title_rank"], lb["title_risk"])
    assert lb["title_rank"][1] > 0.5 and lb["patients"][3] < R.H6_MM - 0.5, lb
    assert lb["title_rank"][3] <= mid - 23.0 and lb["title_risk"][3] <= mid - 23.0, ("triangle titles above the triangles' top", lb["title_rank"], lb["title_risk"], mid - 23.0)
    geo = R.fig6_core(fig, cxi=125.1, tx=202.2, organ_arrow=12.0, mid=mid, tw=30.0, th=TH_STD, pitch=21.0, scale=0.88)   # round 49: cxi 130, tx 199.5; organ arrows 12 as before
    card = geo["inputs"]["card"]; rx1 = card[0] + card[2]
    assert tri["x_extent"][1] + 6.0 < card[0] and tri["x_extent"][0] > 4.0, (tri["x_extent"], card)
    runs = [f["tip"][0] - rx1 for f in geo["fork"]]; assert min(runs) >= 25.0, runs
    for key in ("top", "bottom"):
        t_ = geo["traces"][key]; assert abs(t_["label_pt"] - PT) < 1e-6 and geo["fork"][0 if key == "top" else 1]["tip"][0] <= t_["label_box"][0] - B.HEAD_CLEAR + 1e-6, (key, t_["label_box"])
    rec = finish(fig, stem, outdir, final=final, attribution=None)
    rec.update(geometry=geo, triangles=tri, page_pt=[round(R.W6 / MM, 2), round(R.H6_MM / MM, 2)], trace_h=TH_STD, note_pt=PT, caption_pt=LAB, triangle_pt=TRI_PT, title_gap_mm=round(gap_titles, 2),
               fork_run_mm=[round(r_, 2) for r_ in runs], right_margin_mm=round(R.W6 - geo["content_x"][1], 2), left_margin_mm=round(tri["x_extent"][0], 2),
               card_pad_mm=CARD_PAD_R54,
               moved=f"round 54: notes {PT} pt (were 12), card captions {LAB} (were 10), triangle labels {TRI_PT} (were 12), the card centre at x 125.1 (was 130) with padding {CARD_PAD_R54} mm (was 8), the traces at x 202.2 (was 199.5), "
                     f"the triangle gap 8 mm (was 12), x_left 32.5 (was 34), organ arrows 12 mm as before: the wider card, the wider '90% SpO2' and the 25-mm fork run ({min(runs):.2f}) fit the unchanged width, "
                     f"the rank title's first baseline one leading higher so its second line ends level with 'disease risk' above the triangles. Layout of fig6_tri_r49 otherwise"); return rec


# ---------------------------------------------------------------- 4e: the severity caption below the taller two-line ramp labels, the band taller for it
RAMP_NONE, RAMP_SEVERE = "No apnea\n(AHI <5)", "Severe\n(AHI 30 or more)"


def _ramp_block_mm():
    pf = B.new(80.0, 40.0); pf.theme = replace(pf.theme, note_pt=PT, label_pt=LAB)
    a, b = B.severity_ladder(pf, 4.0, 10.0, 71.0, left_text=RAMP_NONE, right_text=RAMP_SEVERE, tri_h=2.0, ladder=B.AMBER_L, text_color=B.AMBER)
    return round(max(a.y + a.h, b.y + b.h) - 10.0, 3)


RAMP_BLOCK_MM = _ramp_block_mm()                                      # ladder line to the bottom of the two-line ramp-end labels at 15 pt
SEV_LABEL_Y_R54 = round(51.1 + RAMP_BLOCK_MM + 0.8, 2)               # the caption's top 0.8 mm below the ramp labels (round 49: 64.0)
DD_4E = B.SHIFT["band_4e_r40"] - B.B23.SHIFT4                         # the band's content shift (-1.13)
H4_R54 = math.ceil((SEV_LABEL_Y_R54 - DD_4E + 0.925 * LAB * MM + 2.3) * 2) / 2   # band a's 2.3 mm bottom rule under the caption, to 0.5 mm (round 49: 74.0)
_b4 = inspect.getsource(B.band_4e); assert _b4.count("SEV_LABEL_Y = 59.3") == 1 and _b4.count("W4_PT, H4 = B23.W4_PT, B23.H4") == 1, "band_4e literals"
exec(_b4.replace("SEV_LABEL_Y = 59.3", f"SEV_LABEL_Y = {SEV_LABEL_Y_R54}").replace("W4_PT, H4 = B23.W4_PT, B23.H4", f"W4_PT, H4 = B23.W4_PT, {H4_R54}").replace("def band_4e(", "def band_4e_r54src("), B.__dict__)


def themed(fn, stem, note=PT, label=LAB, two_line_ramp=False, tx=None, fork_rule=False, **kw):
    def run(final):
        o_r, o_b = R.new, B.new; o_none, o_sev, o_4e = B.B45.APNEA_NONE, B.B45.APNEA_SEVERE, B.band_4e; o_hb = B.HEAD_BACK; B.HEAD_BACK = HEAD_BACK_R54; o_tx = B.TX_4E5D
        if tx is not None: B.TX_4E5D = tx
        def new54(w, h):
            f = o_r(w, h); f.theme = replace(f.theme, note_pt=note, label_pt=label); return f
        R.new = new54; B.new = new54
        if two_line_ramp: B.B45.APNEA_NONE, B.B45.APNEA_SEVERE, B.band_4e = RAMP_NONE, RAMP_SEVERE, B.band_4e_r54src
        try: rec = fn(final, stem=stem, outdir=OUT, **kw)
        finally: R.new, B.new = o_r, o_b; B.B45.APNEA_NONE, B.B45.APNEA_SEVERE, B.band_4e = o_none, o_sev, o_4e; B.HEAD_BACK = o_hb; B.TX_4E5D = o_tx
        if fork_rule:
            fk = rec["fork"]; run_mm = round(fk["tips_x"] - fk["src"][0], 3); rec["fork_run_mm"] = run_mm; assert run_mm >= 25.0, ("fork run at least 25 mm", stem, run_mm)
            for key in ("top", "bottom"):
                lbx = rec["traces"][key]["label_box"]; assert fk["tips_x"] <= lbx[0] - B.HEAD_CLEAR + 2e-3 and abs(rec["traces"][key]["label_pt"] - note) < 1e-6, (stem, key, lbx, fk["tips_x"])
            rec["arrowhead_to_label_mm"] = round(rec["traces"]["top"]["label_box"][0] - fk["tips_x"], 3); rec["head_back_mm"] = round(HEAD_BACK_R54, 3); rec["tx"] = rec.get("tx")
        if str(rec.get("moved", "")).startswith("round 54"): return rec   # band a, 1f and Figure 6 describe their own changes (the r44 functions' round-44b notes are replaced below)
        rec["moved"] = (f"round 54: notes {note} pt, captions {label} pt (theme), fork arrowheads {HEAD_BACK_R54:.2f} mm before the trace (1.5 mm clear of the 15-pt '90% SpO2')" +
                        (f", the trace pair at x {tx} (was 130.568) for the 25-mm fork run" if tx is not None else "") +
                        (f", ramp labels on two lines, severity caption at {SEV_LABEL_Y_R54} (was 64.0), band height {H4_R54} mm (was 74.0)" if two_line_ramp else "") + ". Geometry unchanged otherwise"); return rec
    return run


SHEETS = {"band_a_r54": themed(band_a_r54, "band_a_r54"), "band_1f_r54": themed(band_1f_r54, "band_1f_r54", label=PT),
          "band_3d_r54": themed(R.band_3d_r44, "band_3d_r54", fork_rule=True), "band_4e_r54": themed(R.band_4e_r44, "band_4e_r54", two_line_ramp=True, tx=TX_4E5D_R54, fork_rule=True),
          "band_5d_r54": themed(R.band_5d_r44, "band_5d_r54", tx=TX_4E5D_R54, fork_rule=True), "Main_Fig6_tri_r54": themed(fig6_tri_r54, "Main_Fig6_tri_r54")}
BANDS = ["band_a_r54", "band_1f_r54", "band_3d_r54", "band_4e_r54", "band_5d_r54"]

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--final", action="store_true"); ap.add_argument("--only", default=""); args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True); (WORK / "previews").mkdir(parents=True, exist_ok=True)
    print(json.dumps({"head_back_r54": round(HEAD_BACK_R54, 3), "ramp_block_mm": RAMP_BLOCK_MM, "sev_label_y_r54": SEV_LABEL_Y_R54, "h4_r54": H4_R54, "tx_4e5d_r54": TX_4E5D_R54, "tx_1f_r54": TX_1F_R54, "h_mm_a": H_MM_A_R54}), flush=True)
    lock = open(B.LOCK, "w"); fcntl.flock(lock, fcntl.LOCK_EX); t0 = time.time()
    try:
        only = [s for s in args.only.split(",") if s]
        for stem, fn in SHEETS.items():
            if only and stem not in only: continue
            try: r = fn(args.final); r["render"] = dict(B.RENDER_INFO)
            except Exception as e:
                import traceback; print(json.dumps({"stem": stem, "FAILED": repr(e)[:400]}), flush=True); traceback.print_exc(); continue
            (WORK / f"records_r54_{stem}.json").write_text(json.dumps(r, indent=1, default=str))
            print(json.dumps({k: r.get(k) for k in ("stem", "w_pt", "h_mm", "audit", "min_pt", "page_pt", "fork_run_mm")}), flush=True)
    finally:
        fcntl.flock(lock, fcntl.LOCK_UN); lock.close()
    print(f"done in {time.time() - t0:.0f} s", flush=True)
