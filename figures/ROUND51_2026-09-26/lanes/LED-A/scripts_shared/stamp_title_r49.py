#!$T90_PY
"""Round 49, lane LED-A: the round-40 title strip re-applied at +1 pt (copy of ROUND40_2026-09-18/lanes/LED/scripts/stamp_title.py, SIZE 13 -> 14).
The page box grows upward by 16 pt (MediaBox y1 + 16 in PDF coordinates: the content stream is untouched, every object sits 16 pt lower in
top-down page coordinates), the title "Extended Data Fig. N" is written with fitz.TextWriter in the system Arial Bold at 14 pt, baseline
y = 12.0, x = the panel letter a's x, in the letter's colour. Placement only. Run as a watched child.
usage: stamp_title_r49.py <Sheet> <N> <in.pdf> <out.pdf> <x> <colour hex> [<strip pt, default 16>]"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import fitz, sys, os, gc, json
S, N, IN, OUT, X, COL = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4], float(sys.argv[5]), sys.argv[6]
STRIP = float(sys.argv[7]) if len(sys.argv) > 7 else 16.0
ARIALB = paths.FONT_ARIAL_BOLD; TITLE_Y = 12.0; SIZE = 14.0     # round 49: 13 + 1
assert os.stat(IN).st_blocks > 0, ("evicted", IN)
d = fitz.open(IN); assert d.page_count == 1; p = d[0]; assert p.rotation == 0
W, H = p.rect.width, p.rect.height; mb = p.mediabox
new_mb = fitz.Rect(mb.x0, mb.y0, mb.x1, mb.y1 + STRIP)
p.set_mediabox(new_mb)
for key in ("CropBox", "BleedBox", "TrimBox", "ArtBox"):
    if d.xref_get_key(p.xref, key)[0] != "null": d.xref_set_key(p.xref, key, "null")
p = d[0]
assert abs(p.rect.width - W) < 0.01 and abs(p.rect.height - (H + STRIP)) < 0.01, (tuple(p.rect), W, H)
title = f"Extended Data Fig. {N}"; font = fitz.Font(fontfile=ARIALB)
tw = fitz.TextWriter(p.rect); tw.append(fitz.Point(X, TITLE_Y), title, font=font, fontsize=SIZE)
rgb = tuple(int(COL.lstrip("#")[i:i + 2], 16) / 255 for i in (0, 2, 4))
tw.write_text(p, color=rgb)
d.save(OUT, garbage=3, deflate=True); d.close(); gc.collect()
rec = dict(sheet=S, input=IN, output=OUT, strip_pt=STRIP, page_in=[W, H], page_out=[W, H + STRIP], title=title, origin=[X, TITLE_Y], size=SIZE, font=ARIALB, colour=COL, width_pt=round(font.text_length(title, fontsize=SIZE), 3))
json.dump(rec, open(OUT[:-4] + "_stamp_record.json", "w"), indent=1); print(json.dumps(rec))
