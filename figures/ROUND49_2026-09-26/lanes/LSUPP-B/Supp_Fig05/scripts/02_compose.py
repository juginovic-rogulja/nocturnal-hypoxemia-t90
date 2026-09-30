"""Supp Fig 5 (V14): place the four rebuilt parts on a fresh page of the V13 sheet's size and stamp the letters a to d with
fitz.TextWriter + Arial Bold 13 pt at the V13 sheet's origins. Writes Supp_Fig05.pdf and Supp_Fig05_150dpi.png."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import json
import os
import sys

LANE = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LSUPP-B/Supp_Fig05"   # R49: this lane's sheet folder (scripts/, work/, verify/ beside it)
sys.path.insert(0, f"{LANE}/scripts")
from l3b_common import compose, hydrated  # noqa: E402
import fitz  # noqa: E402

SHEET = "Supp_Fig05"
SD = LANE   # round 38: outputs in the sheet folder
WORK = f"{SD}/work"
E = json.load(open(hydrated(f"{WORK}/expected_values.json")))
# R49: the letters come from the round-38 verify record of the V13 origins (work/letters_v26.json), no PyMuPDF probe of a composed sheet
LET = json.load(open(f"{WORK}/letters_v26.json")); letters = [(l["text"], l["x"], l["y"]) for l in LET["letters"]]; pw, ph = LET["page"]
assert [l[0] for l in letters] == ["a", "b", "c", "d"], letters
parts = [(f"{WORK}/part_{p}.pdf", *E["place"][p]) for p in "abcd"]
for p, _x, _y in parts:
    assert os.path.exists(p) and os.path.getsize(p) > 0, p
    r = fitz.open(p)[0].rect
    assert abs(r.width - 484.7244) < 0.02, (p, r)
assert (round(pw, 2), round(ph, 2)) == tuple(E["page"]), (pw, ph, E["page"])
out = compose(pw, ph, parts, letters, f"{SD}/{SHEET}.pdf")   # letters at 13 + PT_PLUS = 14 pt
# R49: no render here (the verifier renders with Ghostscript)
print("wrote", out, "and", f"{SD}/{SHEET}_150dpi.png", "letters", letters)
