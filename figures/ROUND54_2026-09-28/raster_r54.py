#!$T90_PY
"""Round 54: 600-dpi raster twins of the vector composes of Figures 3, 4 and 5 (the lanes L3/L4/L5 deliver <Sheet>/<Sheet>.pdf as the vector compose),
one at a time through the round-40/42 rasterizers (their own watchdogs), then copied to lanes/LRASTER/<Sheet>/ with verify/checks.txt.
usage: raster_r54.py Main_Fig3 [Main_Fig4 ...]"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../..")))  # repository root, where paths.py lives
import paths
import json, os, shutil, subprocess, sys, time, hashlib
import fitz
T90 = paths.FIGURE_ROOT; R54 = f"{paths.FIGURE_ROOT}/ROUND54_2026-09-28"; R40M = f"{paths.FIGURE_ROOT}/ROUND40_2026-09-18/lanes/LMAIN"; R42F = f"{paths.FIGURE_ROOT}/ROUND42_2026-09-19/lanes/LFIG5b"
V26 = f"{paths.FIGURE_ROOT}/ROUND47_2026-09-24/figures/NEW_FINAL_SET_V26"; PY = paths.PY
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
def pagesize(pdf):
    from pypdf import PdfReader
    b = PdfReader(pdf).pages[0].mediabox; return float(b.width), float(b.height)
for SHEET in sys.argv[1:]:
    lane = {"Main_Fig3": "L3", "Main_Fig4": "L4", "Main_Fig5": "L5"}[SHEET]; vec = f"{R54}/lanes/{lane}/{SHEET}/{SHEET}.pdf"; assert os.path.exists(vec), vec
    lchk = open(f"{R54}/lanes/{lane}/{SHEET}/verify/checks.txt").read().strip().splitlines(); assert lchk[-1] == "RESULT ALL PASS", (lane, lchk[-1])
    script, twin = (f"{R42F}/scripts/rasterize_r42b.py", f"{R42F}/{SHEET}/work") if SHEET == "Main_Fig5" else (f"{R40M}/scripts/rasterize_r40.py", f"{R40M}/{SHEET}/work")
    SD = f"{R54}/lanes/LRASTER/{SHEET}"; os.makedirs(f"{SD}/work", exist_ok=True); os.makedirs(f"{SD}/verify", exist_ok=True); checks = []
    t0 = time.time(); r = subprocess.run([PY, script, SHEET, "600", "--vector", vec, "--tag", "r54"], capture_output=True, timeout=5400); open(f"{SD}/work/raster.log", "wb").write(r.stdout + r.stderr)
    checks.append(("PASS " if r.returncode == 0 else "FAIL ") + f"600 dpi raster twin rc {r.returncode} in {time.time() - t0:.0f} s")
    if r.returncode == 0:
        for ext in (".pdf", ".tif", "_150dpi.png"): shutil.copy2(f"{twin}/{SHEET}_r54{ext}", f"{SD}/{SHEET}{ext}")
        shutil.copy2(f"{twin}/raster_record_r54.json", f"{SD}/work/raster_record.json")
        w, h = pagesize(f"{SD}/{SHEET}.pdf"); w0, h0 = pagesize(vec); checks.append(("PASS " if abs(w - w0) < 0.05 and abs(h - h0) < 0.05 else "FAIL ") + f"page {w:.2f} x {h:.2f} = the lane vector's {w0:.2f} x {h0:.2f}")
        rr = json.load(open(f"{SD}/work/raster_record.json")); checks.append(("PASS " if rr.get("size_check_within_0.5pt") in (True, "True") and str(rr.get("dpi")) == "600" else "FAIL ") + f"raster record {rr.get('width_px')} x {rr.get('height_px')} px at 600 dpi")
        d = fitz.open(f"{SD}/{SHEET}.pdf"); nimg = len(d[0].get_images()); d.close(); checks.append(("PASS " if nimg == 1 else "FAIL ") + "one image on the raster page")
        checks.append(f"PASS vector source {vec} sha256 {sha(vec)[:16]} (lane {lane} ALL PASS)")
    ok = all(c.startswith("PASS") for c in checks); open(f"{SD}/verify/checks.txt", "w").write("\n".join(checks) + ("\nRESULT ALL PASS\n" if ok else "\nRESULT FAIL\n")); print(SHEET, "ALL PASS" if ok else "FAIL", checks[-1][:80])
