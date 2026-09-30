#!/usr/bin/env python3
"""ED_Fig07 compose, round 51 (lane LED-B, 2026-09-26), Alen's item 4: panels a and b SIDE BY SIDE (a left, b right, tops aligned) instead of stacked.
From the round-49 compose (ed7_03_compose.py beside, unedited): a new fitz page W2 x H2 with panel a (work/panel_a.pdf, unchanged) at its own size
at (0, SHIFT_A) as before, panel b (work/panel_b.pdf, unchanged) at x = X_B = the a page width (484.72) and y = 0 so that the two letters share the
baseline 21.0 and the panel tops align; the page is X_B + (panel b's rightmost ink + EDGE) wide and panel b's own height (523.98) tall (its bottom
margin EDGE), so about double the V28 width and shorter. The 'HR (95% CI)' title is stamped as in round 49 (10.5 pt, right-aligned on the estimate
column), the letters a and b with fitz.TextWriter Arial Bold 14 at (12.96, 21.0) and (X_B + 12.96, 21.0). No text size changes, no string changes."""
import gc, json, os, sys, time
import fitz
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ed7spec as S
L = S.L

LAY49 = json.load(open(S.LAYOUT)); W, H49 = LAY49["sheet_w"], LAY49["sheet_h"]; ra = LAY49["panel_a_rect"]; rb49 = LAY49["panel_b_rect"]; YB = rb49[1]; HB = rb49[3] - rb49[1]
BT = json.load(open(L.hydrated(S.BASE_TEXT))); assert abs(BT["page"][0] - W) < 0.01 and abs(BT["page"][1] - H49) < 0.01
SPV = L.dedupe(BT["spans"]); ink_span = [s for s in SPV if s["text"] == "a" and s["size"] == 13.0][0]; INK = L.hex_rgb01(ink_span["color"])
DA = json.load(open(S.DRAWN_A))
est_new = [r for r in DA["drawn_records"] if r["rule"] == "sev_hr_ci"]; head_new = [r for r in DA["drawn_records"] if r["text"] == S.BLOCK_HEAD][0]; BASE_Y = head_new["baseline"]
est_sp = [s for s in L.spans_of(L.hydrated(S.PANEL_A)) if abs(s["size"] - (9.5 + L.PT_PLUS)) < 0.01 and s["origin"][0] > 380]; assert len(est_sp) == len(est_new)
COL_X1 = max(s["bbox"][2] for s in est_sp)
# panel b's rightmost ink on its own flat page (text and line art)
db = fitz.open(L.hydrated(S.PANEL_B)); pb = db[0]; assert abs(pb.rect.width - W) < 0.01 and abs(pb.rect.height - HB) < 0.01, (pb.rect, W, HB)
b_text_right = max(s["bbox"][2] for s in L.spans_of(pb)); b_art_right = max(it["rect"].x1 for it in pb.get_drawings()); B_RIGHT = max(b_text_right, b_art_right)
X_B = W; B_CLIP_W = B_RIGHT + S.EDGE; W2 = X_B + B_CLIP_W; H2 = HB
doc = fitz.open(); page = doc.new_page(width=W2, height=H2)
da = fitz.open(L.hydrated(S.PANEL_A)); assert abs(da[0].rect.width - W) < 0.01 and abs(da[0].rect.height - (ra[3] - ra[1])) < 0.01
page.show_pdf_page(fitz.Rect(*ra), da, 0)
font = fitz.Font(fontfile=L.ARIAL); tw_ = font.text_length(S.HDR_A, fontsize=S.HDR_A_FS); x0 = COL_X1 - tw_
tw = fitz.TextWriter(page.rect); tw.append(fitz.Point(x0, BASE_Y), S.HDR_A, font=font, fontsize=S.HDR_A_FS); tw.write_text(page, color=INK)
rb = [X_B, 0.0, X_B + B_CLIP_W, H2]
page.show_pdf_page(fitz.Rect(*rb), db, 0, clip=fitz.Rect(0, 0, B_CLIP_W, HB))
letters = [["a", S.LETTER_X, S.LETTER_BASELINE], ["b", X_B + S.LETTER_X, S.LETTER_BASELINE]]
for ch, x, base in letters: L.stamp_letter(page, ch, x, base, size=LAY49["letter_size"], color=INK)
doc.save(S.OUT_PDF, garbage=3, deflate=True); doc.close(); da.close(); db.close(); gc.collect()
lay = dict(sheet_w=W2, sheet_h=H2, v28_sheet=[W, H49], panel_a_rect=ra, panel_b_rect=rb, panel_b_clip=[0, 0, B_CLIP_W, HB], panel_b_src_rect=rb49, x_b=X_B, b_right_ink=B_RIGHT, shift_a=LAY49["shift_a"], letters=letters, letter_size=LAY49["letter_size"], edge=S.EDGE)
json.dump(lay, open(f"{S.WORK}/sheet_layout_r51.json", "w"), indent=1)
log = dict(round=51, base=S.BASE_PDF, base_sha256=L.sha256(S.BASE_PDF), page=[W2, H2], v28_page=[W, H49], placements=[dict(what="panel a (round-49 panel, unchanged)", rect=ra, file=S.PANEL_A, sha256=L.sha256(S.PANEL_A)),
           dict(what=f"HR (95% CI) header, TextWriter Arial regular {S.HDR_A_FS} pt, right-aligned on the estimate column", x0=x0, x1=COL_X1, baseline=BASE_Y),
           dict(what="panel b (round-49 panel, unchanged) placed to the right, top aligned", rect=rb, clip=[0, 0, B_CLIP_W, HB], file=S.PANEL_B, sha256=L.sha256(S.PANEL_B)), dict(what="letters a, b (TextWriter Arial Bold 14 pt) at the panels' top-left corners", letters=letters)],
           output=S.OUT_PDF, output_sha256=L.sha256(S.OUT_PDF), bytes=os.path.getsize(S.OUT_PDF), written=time.strftime("%Y-%m-%d %H:%M:%S"))
json.dump(log, open(f"{S.WORK}/compose_log.json", "w"), indent=1)
print(f"wrote {S.OUT_PDF}: page {W2:.2f} x {H2:.2f} (V28 {W:.2f} x {H49:.2f}); panel b at x {X_B:.2f}, its ink to {B_RIGHT:.2f} (+ edge {S.EDGE:.2f}); letters {letters}")
