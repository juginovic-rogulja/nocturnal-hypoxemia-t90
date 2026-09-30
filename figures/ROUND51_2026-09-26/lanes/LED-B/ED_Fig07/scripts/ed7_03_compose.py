#!/usr/bin/env python3
"""Round 49 LED-B copy (2026-09-26): header 10.5 pt right-aligned on the new estimate column edge, letters 14 pt. Otherwise the round-38 copy. ED_Fig07 compose (V14 lane L4, round 37): a new fitz page of the V13 size (484.72 x 900.628 pt, LF4 layout). Panel a (work/panel_a.pdf,
the forest re-plotted from the v8.1 numbers) shown at its own size SHIFT_A pt down the page, the column title "HR (95% CI)" stamped with
fitz.TextWriter in Arial regular 9.5 pt, ink, right-aligned on the printed estimate column at the "Negative controls" header baseline
(the round-31 stamp), panel b (work/panel_b.pdf, cross9 in Fig 4d's design) below, the letters a and b stamped with fitz.TextWriter +
Arial Bold 13 pt in the V13 ink at the V13 origins. Lineage: the round-30 03_compose.py and the round-31 LF4 02_compose.py (beside as
ed7_03_compose_PRE_V8_1.py). Never insert_text with a fontfile, never fontTools, never Quartz."""
import gc, json, os, sys, time
import fitz
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ed7spec as S
L = S.L

LAY = json.load(open(S.LAYOUT)); W, H = LAY["sheet_w"], LAY["sheet_h"]; ra, rb = LAY["panel_a_rect"], LAY["panel_b_rect"]
BT = json.load(open(L.hydrated(S.BASE_TEXT))); assert abs(BT["page"][0] - W) < 0.01 and abs(BT["page"][1] - H) < 0.01, (BT["page"], W, H)
SPV = L.dedupe(BT["spans"]); ink_span = [s for s in SPV if s["text"] == "a" and s["size"] == 13.0][0]; INK = L.hex_rgb01(ink_span["color"])
DA = json.load(open(S.DRAWN_A))
# the estimate column right edge and the header baseline: from the new panel a's own drawn records (the V13 pins, checked below against V13)
est_new = [r for r in DA["drawn_records"] if r["rule"] == "sev_hr_ci"]; assert len(est_new) == DA["n_controls"] + DA["n_conditions"]
head_new = [r for r in DA["drawn_records"] if r["text"] == S.BLOCK_HEAD][0]; BASE_Y = head_new["baseline"]
est13 = [s for s in SPV if s["font"] == "ArialMT" and abs(s["size"] - 9.5) < 0.01 and s["origin"][0] > 380 and s["origin"][1] < 384]; COL_X1_13 = max(s["bbox"][2] for s in est13)
hdr13 = [s for s in SPV if s["text"].replace("\xa0", " ") == S.HDR_A and s["origin"][1] < 100]; assert len(hdr13) == 1, hdr13
assert abs(hdr13[0]["bbox"][2] - COL_X1_13) < 0.1, (hdr13[0]["bbox"], COL_X1_13)
BASE_Y_13 = hdr13[0]["origin"][1]   # round 49: the header baseline follows the new "Negative controls" baseline (matplotlib centres the taller text on the row: 0.26 pt lower), declared in the log
# round 49: the estimate column is 10.5 pt, so its right edge moves right; the header is right-aligned on the NEW column's right edge, read from the new panel's own text layer (a flat matplotlib panel)
est_sp = [s for s in L.spans_of(L.hydrated(S.PANEL_A)) if abs(s["size"] - (9.5 + L.PT_PLUS)) < 0.01 and s["origin"][0] > 380]; assert len(est_sp) == len(est_new), (len(est_sp), len(est_new))
COL_X1 = max(s["bbox"][2] for s in est_sp); assert abs(min(s["origin"][0] for s in est_sp) - min(r["x"] for r in est_new)) < 0.05
doc = fitz.open(); page = doc.new_page(width=W, height=H)
da = fitz.open(L.hydrated(S.PANEL_A)); assert abs(da[0].rect.width - W) < 0.01 and abs(da[0].rect.height - (ra[3] - ra[1])) < 0.01, da[0].rect
page.show_pdf_page(fitz.Rect(*ra), da, 0)
font = fitz.Font(fontfile=L.ARIAL); tw_ = font.text_length(S.HDR_A, fontsize=S.HDR_A_FS); x0 = COL_X1 - tw_
tw = fitz.TextWriter(page.rect); tw.append(fitz.Point(x0, BASE_Y), S.HDR_A, font=font, fontsize=S.HDR_A_FS); tw.write_text(page, color=INK)
db = fitz.open(L.hydrated(S.PANEL_B)); assert abs(db[0].rect.width - W) < 0.01 and abs(db[0].rect.height - (rb[3] - rb[1])) < 0.01, db[0].rect
page.show_pdf_page(fitz.Rect(*rb), db, 0)
for ch, x, base in LAY["letters"]: L.stamp_letter(page, ch, x, base, size=LAY["letter_size"], color=INK)
doc.save(S.OUT_PDF, garbage=3, deflate=True); doc.close(); da.close(); db.close(); gc.collect()
log = dict(base=S.BASE_PDF, base_sha256=L.sha256(S.BASE_PDF), page=[W, H], placements=[dict(what="panel a (new, v8.1)", rect=ra, file=S.PANEL_A, sha256=L.sha256(S.PANEL_A)),
           dict(what=f"HR (95% CI) header, TextWriter Arial regular {S.HDR_A_FS} pt, right-aligned on the estimate column (round 49: the column edge moved from {COL_X1_13:.2f} to {COL_X1:.2f})", x0=x0, x1=COL_X1, x1_v13=COL_X1_13, baseline=BASE_Y, baseline_v13=BASE_Y_13),
           dict(what="panel b (new, v8.1)", rect=rb, file=S.PANEL_B, sha256=L.sha256(S.PANEL_B)), dict(what="letters a, b (TextWriter Arial Bold 13 pt, V13 origins)", letters=LAY["letters"])],
           output=S.OUT_PDF, output_sha256=L.sha256(S.OUT_PDF), bytes=os.path.getsize(S.OUT_PDF), written=time.strftime("%Y-%m-%d %H:%M:%S"))
json.dump(log, open(f"{S.WORK}/compose_log.json", "w"), indent=1)
print(f"wrote {S.OUT_PDF}: page {W:.2f} x {H:.2f}; header at x {x0:.2f}..{COL_X1:.2f}, baseline {BASE_Y:.3f}; {os.path.getsize(S.OUT_PDF)/1e6:.2f} MB")
