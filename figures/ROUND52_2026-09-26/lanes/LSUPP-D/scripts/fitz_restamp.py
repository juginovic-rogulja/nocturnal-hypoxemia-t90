#!$T90_PY
"""Lane LSUPP-D child (under the watchdog, FLAT single-page files only): remove EVERY text object of the page (text-only
redaction over the whole page, or only over the given rects when whole_page is false: images and vector graphics untouched),
drop the then-unused page font resources, and re-set the planned strings with fitz.TextWriter in the sheet's Arial faces
(regular or bold) one point larger on the old string's own anchor. Anchors are computed with the SAME face metrics on both
sides: w_old = advance of the text at the old size, w = advance at the new size; align left keeps the old origin x, right keeps
ox + w_old, centre keeps ox + w_old / 2; rot -90 (a bottom-to-top axis title) keeps the baseline x and the centre of the run
(oy - w_old / 2). The TextWriter faces' ToUnicode quirk (space and hyphen mapped to U+00A0 and U+00AD) is repointed afterwards
(overlay.fix_tounicode). Records every drawn box in <out>.drawn.json.
usage: fitz_restamp.py job.json
job = {src, out, whole_page: true, redact_rects: [[x0,y0,x1,y1], ...], items: [{text, old_size, size, bold, color, ox, oy, align, rot}], metadata}"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import json, sys, fitz
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LSUPP-D/scripts")
import overlay as O
job = json.load(open(sys.argv[1]))
d = fitz.open(job["src"]); assert len(d) == 1; p = d[0]
whole = job.get("whole_page", True)
MB0 = fitz.Rect(p.mediabox); R0 = fitz.Rect(p.rect); PAD = 3000.0
if whole:
    # text that lies entirely OFF the page (ghost copies of earlier generations) is outside every on-page redaction rect: widen the
    # MediaBox temporarily so that the whole-page redaction reaches it, then restore the original boxes
    p.set_mediabox(fitz.Rect(MB0.x0 - PAD, MB0.y0 - PAD, MB0.x1 + PAD, MB0.y1 + PAD))
    p.add_redact_annot(p.rect, fill=False)
    assert p.apply_redactions(images=0, graphics=0, text=0), "apply_redactions returned False"
    p.set_mediabox(MB0)
    assert fitz.Rect(p.mediabox) == MB0 and fitz.Rect(p.cropbox) == MB0 and abs(p.rect.width - R0.width) < 1e-3 and abs(p.rect.height - R0.height) < 1e-3, (p.mediabox, p.cropbox, p.rect, R0)
if job.get("redact_rects"):
    for r in job["redact_rects"]:
        p.add_redact_annot(fitz.Rect(*r), fill=False)
    assert p.apply_redactions(images=0, graphics=0, text=0), "apply_redactions returned False"
dropped = []
if whole:
    left = [s for b in p.get_text("rawdict")["blocks"] if b["type"] == 0 for l in b["lines"] for s in l["spans"] if "".join(c["c"] for c in s["chars"]).strip()]
    assert not left, f"text left after the redaction: {[''.join(c['c'] for c in s['chars']) for s in left][:10]}"
    for f in p.get_fonts(full=True):        # no text uses the page-level fonts any more; the TextWriter adds its own
        d.xref_set_key(p.xref, f"Resources/Font/{f[4]}", "null"); dropped.append(f[3])
def rot_matrix_for(src):
    """The sense of a TextWriter morph rotation depends on the page's own content stream (a fresh page and a pdfwrite page differ):
    probe on a scratch copy of the source and return the matrix whose result reads bottom to top (line direction (0, -1))."""
    for deg in (-90, 90):
        dd = fitz.open(src); pp = dd[0]; tw = fitz.TextWriter(pp.rect); tw.append(fitz.Point(100, 100), "probe", font=O.font(False), fontsize=10)
        tw.write_text(pp, color=(0, 0, 0), morph=(fitz.Point(100, 100), fitz.Matrix(deg)))
        dirs = [l["dir"] for b in pp.get_text("rawdict")["blocks"] if b["type"] == 0 for l in b["lines"] if "probe" in "".join(c["c"] for sp in l["spans"] for c in sp["chars"])]
        dd.close()
        if dirs and dirs[0][1] < -0.5:
            return deg
    raise RuntimeError("no morph rotation reads bottom to top on this page")
ROT_DEG = rot_matrix_for(job["src"]) if any(it.get("rot", 0) == -90 for it in job.get("items", [])) else None
drawn = []
for it in job.get("items", []):
    f = O.font(it.get("bold", False)); size = float(it["size"]); osz = float(it["old_size"])
    w_old = f.text_length(it["text"], fontsize=osz); w = f.text_length(it["text"], fontsize=size)
    asc, desc = f.ascender * size, -f.descender * size
    col = O.hexrgb(it.get("color", "#1a1d21")); ox, oy = float(it["ox"]), float(it["oy"])
    if it.get("rot", 0) == -90:
        yc = oy - w_old / 2; y0 = yc + w / 2
        tw = fitz.TextWriter(p.rect); tw.append(fitz.Point(ox, y0), it["text"], font=f, fontsize=size)
        tw.write_text(p, color=col, morph=(fitz.Point(ox, y0), fitz.Matrix(ROT_DEG)))
        drawn.append(dict(it, x0=ox, y_start=round(y0, 3), yc=round(yc, 3), w_old=round(w_old, 3), width=round(w, 3), box=[round(v, 3) for v in (ox - asc, y0 - w, ox + desc, y0)]))
    else:
        al = it.get("align", "left"); ob = it.get("old_bbox")      # the old string's ACTUAL extents (PyMuPDF bbox: origin + advances as placed, kerning included)
        x1_old = ob[2] if ob else ox + w_old
        anchor = ox if al == "left" else (x1_old if al == "right" else (ox + x1_old) / 2)
        x0 = anchor if al == "left" else (anchor - w if al == "right" else anchor - w / 2)
        tw = fitz.TextWriter(p.rect); tw.append(fitz.Point(x0, oy), it["text"], font=f, fontsize=size); tw.write_text(p, color=col)
        drawn.append(dict(it, x0=round(x0, 3), anchor=round(anchor, 3), w_old=round(w_old, 3), width=round(w, 3), box=[round(v, 3) for v in (x0, oy - asc, x0 + w, oy + desc)]))
n_cmap = O.fix_tounicode(d)
if job.get("metadata"):
    d.set_metadata(job["metadata"])
d.save(job["out"], garbage=3, deflate=True); d.close()
json.dump({"dropped_fonts": dropped, "cmaps_patched": n_cmap, "rot_matrix_deg": ROT_DEG, "drawn": drawn}, open(job["out"] + ".drawn.json", "w"), indent=1)
print("restamp ok:", len(drawn), "strings,", len(dropped), "old font resources dropped,", n_cmap, "CMaps patched ->", job["out"])
