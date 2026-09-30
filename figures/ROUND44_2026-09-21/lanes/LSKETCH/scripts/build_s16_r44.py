#!$T90_PY
"""Supp_Fig16, round 40 (lane LSUPP). The V16 sheet is several generations of overlays on nested forms whose text layer still carries the
older sketch bands (and a hidden second copy of "10 of 12 / above the line"). Route, all Ghostscript and PyMuPDF placement/redaction in
children under the watchdog:
  1. Ghostscript pdfwrite re-distills the V16 sheet into ONE flat content stream (no downsampling, colours untouched, fonts embedded);
  2. PyMuPDF text-only redaction of the band area (y < 80) and of the "10 of 12 / above the line" box removes those spans (both copies);
  3. a fresh page: the flat base placed clipped to below the band (y >= 79) plus band_s16_r40.pdf (496.80 x 79.00, sha checked) at the top;
  4. panel b's key handles are checked for one marker each (blob count) and left as they are if so.
Proofs in r40_verify_s16.py: Ghostscript text layer = V16's below the band minus the two removed strings (x2), band region = the r40
band's own text, no old band string anywhere; render identical below the band; OCR read-back."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import hashlib, json, os, sys
L = f"{paths.FIGURE_ROOT}/ROUND40_2026-09-18/lanes/LSUPP"
sys.path.insert(0, f"{L}/scripts")
import lsupp_common as C
S = "Supp_Fig16"; D = f"{paths.FIGURE_ROOT}/ROUND44_2026-09-21/lanes/LSKETCH/Supp_Fig16"; W = f"{D}/work"; os.makedirs(W, exist_ok=True)
V16 = f"{C.V16}/{S}.pdf"; C.hydrated(V16)
BAND = f"{paths.FIGURE_ROOT}/ROUND44_2026-09-21/lanes/LSKETCH/work/build/band_s16_r44.pdf"; C.hydrated(BAND)
BAND_SHA = "27593c17c103a60e8421f06059d489991cd122e475c19a72940e8eae414ee669"
assert C.sha256(BAND) == BAND_SHA, ("band sha256 differs from the coordinator's", C.sha256(BAND))
bw, bh = C.page_box(BAND); assert abs(bw - 496.8) < 0.05 and abs(bh - 79.0) < 0.05, (bw, bh)
W16, H16 = C.page_box(V16); assert abs(W16 - 496.8) < 0.05 and abs(H16 - 1131.7) < 0.05, (W16, H16)
BAND_H = 79.0
FLAT = f"{W}/{S}_flat.pdf"; RED = f"{W}/{S}_flat_redacted.pdf"; NEW = f"{D}/{S}.pdf"
# 1 flatten
C.wd(f"{S}_pdfwrite", [C.GS, "-dSAFER", "-dBATCH", "-dNOPAUSE", "-dQUIET", "-sDEVICE=pdfwrite", "-dCompatibilityLevel=1.7", "-dPDFSETTINGS=/prepress",
                        "-dDownsampleColorImages=false", "-dDownsampleGrayImages=false", "-dDownsampleMonoImages=false", "-dAutoFilterColorImages=false", "-dAutoFilterGrayImages=false",
                        "-dColorImageFilter=/FlateEncode", "-dGrayImageFilter=/FlateEncode", "-dColorConversionStrategy=/LeaveColorUnchanged", "-dEmbedAllFonts=true", "-dSubsetFonts=true",
                        "-dPreserveAnnots=false", f"-sOutputFile={FLAT}", V16], max_s=600)
# 2 redact (text only): the band area and the two copies of "10 of 12" / "above the line" (V16 spans at x 147/226, y 346 and x 121/200, y 357)
TEXT_BOX = [115.0, 336.0, 270.0, 361.0]
job = {"mode": "redact", "src": FLAT, "out": RED, "rects": [[0.0, 0.0, W16, BAND_H + 1.0], TEXT_BOX]}
json.dump(job, open(f"{W}/redact_job.json", "w")); C.wd(f"{S}_redact", [C.PY, f"{L}/scripts/fitz_place.py", f"{W}/redact_job.json"])
# 3 compose: flat base below the band + the r40 band
job = {"mode": "compose", "out": NEW, "page": [W16, H16], "parts": [{"pdf": RED, "rect": [0, BAND_H, W16, H16], "clip": [0, BAND_H, W16, H16]}, {"pdf": BAND, "rect": [0, 0, bw, bh]}],   # rect = clip: show_pdf_page fits the CLIP into the rect, so a 1:1 placement needs the same box
       "letters": [], "metadata": {"title": S, "creator": "LSUPP r40 build_s16_r44.py: gs pdfwrite flatten + fitz text redaction + fitz placement (band_s16_r44)"}}
json.dump(job, open(f"{W}/compose_job.json", "w")); C.wd(f"{S}_compose", [C.PY, f"{L}/scripts/fitz_place.py", f"{W}/compose_job.json"])
json.dump({"sheet": S, "v16": V16, "v16_sha256": C.sha256(V16), "band": BAND, "band_sha256": BAND_SHA, "band_box": [bw, bh], "flat": FLAT, "flat_sha256": C.sha256(FLAT), "redacted": RED,
           "redact_rects": job and [[0.0, 0.0, W16, BAND_H + 1.0], TEXT_BOX], "new": NEW, "new_sha256": C.sha256(NEW), "page": [W16, H16], "band_h": BAND_H, "written": C.now()},
          open(f"{W}/build_record.json", "w"), indent=1)
print("wrote", NEW, "page", C.page_box(NEW))
