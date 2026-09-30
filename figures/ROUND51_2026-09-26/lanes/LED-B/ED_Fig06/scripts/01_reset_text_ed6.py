#!$T90_PY
"""ED_Fig06, round 49 (lane LED-B, 2026-09-26): every text one point larger. The sheet has no from-scratch builder (the round-37 builder is an
overlay on the V13 sheet, whose text is a Chromium export; Chrome is banned in the data lanes), so the text layer of the delivered V26 sheet is
RE-SET in place: one full-page TEXT-ONLY redaction (fill=False, line art and images untouched) removes every text object, including the 55
clipped stale copies of panel a's text that Ghostscript lists (invisible on the render, LNOTES round 40 open point 1), then every one of the 116
visible strings is written back with fitz.TextWriter in the system Arial faces at its V26 size + 1 pt, in its V26 colour and weight, on its
V26 baseline, anchored the way the layout reads: centred strings on their V26 centre (panel a headers, counts, 'reference', tick labels, the
axis title, the rotated row-axis label on the grid's vertical centre), right-aligned strings on their V26 right edge (panel a row labels),
left-aligned strings at their V26 x (letters, title, outcome headers, HR column, key). Two text-only nudges so that the larger text does not
collide on panel b (declared): the 24 row labels move from x 27.36 to 21.0 (indent under the bold outcome headers 14.4 -> 8.0 pt) and the
events/group-size column with its two header lines is right-aligned at 248.0 instead of 240.4 (the plot's leftmost visible ink is at x 262.5).
Runs as a child under wd_run.sh. usage: 01_reset_text_ed6.py <in V26.pdf> <work dir> <out.pdf>"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import fitz, json, os, sys, gc, re, hashlib
IN, WORK, OUT = sys.argv[1:4]; os.makedirs(WORK, exist_ok=True)
ARIAL = paths.FONT_ARIAL; ARIALB = paths.FONT_ARIAL_BOLD
PLUS = 1.0; LAB_X_NEW = 21.0; EV_RIGHT_NEW = 248.0
assert os.stat(IN).st_blocks > 0
FR, FB = fitz.Font(fontfile=ARIAL), fitz.Font(fontfile=ARIALB)
def width(t, size, bold): return (FB if bold else FR).text_length(t, fontsize=size)
d = fitz.open(IN); p = d[0]; W, H = p.rect.width, p.rect.height; assert len(p.get_xobjects()) == 0
spans = []
for b in p.get_text("rawdict")["blocks"]:
    if b["type"] != 0: continue
    for l in b["lines"]:
        for s in l["spans"]:
            txt = "".join(ch["c"] for ch in s["chars"]).replace("\xa0", " ").replace("\xad", "-")
            if not txt.strip(): continue
            spans.append(dict(text=txt, font=s["font"], size=round(s["size"], 3), color="#%06x" % s["color"], origin=[round(v, 3) for v in s["chars"][0]["origin"]], bbox=[round(v, 3) for v in s["bbox"]], dir=[round(v, 3) for v in l["dir"]]))
assert len(spans) == 116, len(spans)
# ---- anchor rules
EV_RE = re.compile(r"^\d+/[\d,]+$")
placed = []
for s in spans:
    t, size, bold, col = s["text"], s["size"], "Bold" in s["font"], s["color"]; x, y = s["origin"]; x0, x1 = s["bbox"][0], s["bbox"][2]
    w_old = x1 - x0 if not (s["dir"] == [0.0, -1.0]) else width(t, size, bold); c_old = (x0 + x1) / 2; r_old = x1     # the V26 advance box (kerned on the Chromium-exported strings, e.g. 'Ty')
    assert abs(width(t, size, bold) - (x1 - x0)) < 1.5 or s["dir"] == [0.0, -1.0], (t, x0, x1, width(t, size, bold))
    nsize = size + PLUS; w_new = width(t, nsize, bold); rot = s["dir"] == [0.0, -1.0]
    if rot:
        kind, anchor = "centre_vertical", y - w_old / 2; nx, ny = x, anchor + w_new / 2; note = "rotated row-axis label, centred on the grid's vertical centre"
    elif t.startswith("Extended Data Fig.") or (len(t) == 1 and t in "ab" and size == 13.0) or (bold and y > 215) or abs(x - 365.44) < 0.3 or abs(x - 458.97) < 0.3:
        kind, anchor = "left", x; nx, ny = x, y; note = "left-aligned at the V26 x"
    elif y < 210 and x1 <= 108.5 and size == 10.0:
        kind, anchor = "right", r_old; nx, ny = r_old - w_new, y; note = "panel a row label, right-aligned on the V26 right edge"
    elif y < 210 or (abs(y - 626.88) < 0.5 and size == 10.0) or (size == 11.0 and y > 640):
        kind, anchor = "centre", c_old; nx, ny = c_old - w_new / 2, y; note = "centred on the V26 centre"
    elif EV_RE.match(t) or t in ("New diagnoses /", "group size"):
        kind, anchor = "right", EV_RIGHT_NEW; nx, ny = EV_RIGHT_NEW - w_new, y; note = f"events/group-size column, right-aligned at {EV_RIGHT_NEW} (V26 right edge {r_old:.2f}, nudged +{EV_RIGHT_NEW - r_old:.2f})"
    elif abs(x - 27.36) < 0.3 and size == 10.0 and y > 215:
        kind, anchor = "left", LAB_X_NEW; nx, ny = LAB_X_NEW, y; note = f"panel b row label, left-aligned at {LAB_X_NEW} (V26 x 27.36, nudged {LAB_X_NEW - 27.36:+.2f})"
    else: raise SystemExit(f"no anchor rule for {s}")
    placed.append(dict(text=t, size_old=size, size=nsize, bold=bold, color=col, x_old=x, y=y, bbox_old=s["bbox"], w_old=round(w_old, 3), w_new=round(w_new, 3), kind=kind, anchor=round(anchor, 3), x=round(nx, 3), baseline=round(ny, 3), rotated=rot, note=note))
# collision rules on panel b at the new sizes (label end vs events start on the same baseline, events end vs the plot's leftmost visible ink)
labs = [q for q in placed if q["kind"] == "left" and abs(q["x"] - LAB_X_NEW) < 0.01]; evs = [q for q in placed if q["kind"] == "right" and abs(q["anchor"] - EV_RIGHT_NEW) < 0.01 and EV_RE.match(q["text"])]
gaps = []
for q in labs:
    e = min(evs, key=lambda e: abs(e["baseline"] - q["baseline"])); assert abs(e["baseline"] - q["baseline"]) < 0.5
    gaps.append((q["text"], e["text"], round((e["x"]) - (q["x"] + q["w_new"]), 2)))
MIN_GAP = min(g[2] for g in gaps); assert MIN_GAP >= 6.0, gaps
LEFTMOST_INK = 262.5   # the CI line of 0.67 (0.35-1.27): x(0.35) on the sheet's own log axis (ticks 0.5 at 273.39 ... 4 at 336.82)
assert EV_RIGHT_NEW + 6.0 < LEFTMOST_INK
# ---- 1 full-page text-only redaction (line art NONE, images NONE): every text object leaves, the graphics stay byte for byte in effect
p.add_redact_annot(p.rect, fill=False)
assert p.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE, graphics=fitz.PDF_REDACT_LINE_ART_NONE, text=fitz.PDF_REDACT_TEXT_REMOVE)
p = d[0]; assert not p.get_text("text").strip(), "text remains after the redaction"
dropped = []
for f in p.get_fonts(full=True):
    if "+" in f[3]:            # the Chromium subset faces (PYPZVP+ArialMT, DUNYCV+Arial-BoldMT), now unused
        d.xref_set_key(p.xref, f"Resources/Font/{f[4]}", "null"); dropped.append(f[3])
d.save(f"{WORK}/ED_Fig06_NOTEXT.pdf", garbage=3, deflate=True)
# ---- 2 write every string back at +1 pt
writers = {}
for q in placed:
    key = (q["color"], q["bold"], q["rotated"]); tw = writers.setdefault(key, fitz.TextWriter(p.rect))
    tw.append(fitz.Point(q["x"], q["baseline"]), q["text"], font=FB if q["bold"] else FR, fontsize=q["size"])
for (col, bold, rot), tw in writers.items():
    rgb = tuple(int(col.lstrip("#")[i:i + 2], 16) / 255 for i in (0, 2, 4))
    if rot:
        q = [q for q in placed if q["rotated"]]; assert len(q) == 1; q = q[0]
        tw.write_text(p, color=rgb, morph=(fitz.Point(q["x"], q["baseline"]), fitz.Matrix(90)))
    else: tw.write_text(p, color=rgb)
# ToUnicode of the TextWriter faces: space glyph -> U+0020, hyphen glyph -> U+002D (the round-29b overlay idiom)
n_cmap = 0; seen = set()
for f in p.get_fonts(full=True):
    if f[3] not in ("Arial Regular", "Arial Bold") or f[0] in seen: continue
    seen.add(f[0]); tu = d.xref_get_key(f[0], "ToUnicode")
    if tu[0] != "xref": continue
    tx = int(tu[1].split()[0]); s = d.xref_stream(tx).decode("latin-1"); s2 = s.replace("<0003> <00a0>", "<0003> <0020>").replace("<0010> <00ad>", "<0010> <002d>")
    if s2 != s: d.update_stream(tx, s2.encode("latin-1")); n_cmap += 1
d.save(OUT, garbage=3, deflate=True); d.close(); gc.collect()
rec = dict(sheet="ED_Fig06", input=IN, input_sha256=hashlib.sha256(open(IN, "rb").read()).hexdigest(), output=OUT, output_sha256=hashlib.sha256(open(OUT, "rb").read()).hexdigest(), page=[W, H], plus=PLUS,
           n_visible_spans=len(spans), placed=placed, nudges=dict(row_labels_x=[27.36, LAB_X_NEW], events_right_edge=[240.4, EV_RIGHT_NEW]), panel_b_gaps=gaps, min_gap_pt=MIN_GAP, leftmost_plot_ink_x=LEFTMOST_INK,
           dropped_font_resources=dropped, cmaps_patched=n_cmap, method="full-page add_redact_annot(fill=False) + apply_redactions(images NONE, graphics NONE, text REMOVE); TextWriter Arial Regular/Bold; rotated label via morph Matrix(90)")
json.dump(rec, open(f"{WORK}/reset_record.json", "w"), indent=1, ensure_ascii=False)
print(f"placed {len(placed)} strings at +{PLUS:g} pt; panel b min gap label-to-events {MIN_GAP} pt; dropped {dropped}; cmaps patched {n_cmap}; wrote {OUT}")
