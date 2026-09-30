#!$T90_PY
"""Round 54 copy (lane L5, paths only). Round 49, lane L5: Ghostscript renders under the lane watchdog. usage: 06b_render_r49.py <pdf> <png> <dpi> <name>"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import sys, os, json, time
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND54_2026-09-28/lanes/L5/scripts")
import l5_lib as L
pdf, png, dpi, name = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
t0 = time.time(); st = L.gs_render(pdf, png, dpi, name, max_s=2400)
from PIL import Image
im = Image.open(png); print(f"rendered {png}: {im.width} x {im.height} px at {dpi} dpi in {time.time() - t0:.0f} s; watchdog {st}")
json.dump(dict(pdf=pdf, pdf_sha256=L.sha256(pdf), png=png, png_sha256=L.sha256(png), dpi=dpi, px=[im.width, im.height], seconds=round(time.time() - t0, 1), watchdog=st, written=L.now()), open(png[:-4] + "_render_record.json", "w"), indent=1)
