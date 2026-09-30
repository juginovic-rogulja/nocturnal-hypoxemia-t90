"""Supp Fig 7 (V14): place the rebuilt sheet on a fresh page of the V13 size and stamp a, b, c with TextWriter + Arial Bold 13
at the V13 sheet's letter origins (read from the V13 text layer). Writes Supp_Fig07.pdf and Supp_Fig07_150dpi.png."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import os, sys
LANE = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LSUPP-B"   # R49
sys.path.insert(0, f"{LANE}/Supp_Fig07/scripts")
from l3b_common import compose, hydrated  # noqa: E402   (R49: letters 14 pt through compose's default)
import json  # noqa: E402
import fitz  # noqa: E402
SHEET = "Supp_Fig07"; SD = f"{LANE}/{SHEET}"; WORK = f"{SD}/work"
# R40: the letters come from the round-37 verify record of the V16 sheet (a, b, c at x 12.96, y 22.99 / 258.43 / 506.84), no PyMuPDF probe of a composed sheet
LET = json.load(open(f"{WORK}/letters_v16.json")); letters = [(l["text"], l["x"], l["y"]) for l in LET["letters"]]; pw, ph = LET["page"]
assert [l[0] for l in letters] == ["a", "b", "c"], letters
part = f"{WORK}/sheet.pdf"; assert os.path.getsize(part) > 0
pd_ = fitz.open(part); pr = pd_[0].rect; pd_.close()
assert abs(pr.width - pw) < 0.05, (pr.width, pw)
out = compose(pw, pr.height, [(part, 0.0, 0.0)], letters, f"{SD}/{SHEET}.pdf")
# R40: no render here (this script runs as a child under the watchdog; the verifier renders with Ghostscript)
print("wrote", out, f"page {pw} x {pr.height} (V13 {pw} x {ph}), letters {letters}")
