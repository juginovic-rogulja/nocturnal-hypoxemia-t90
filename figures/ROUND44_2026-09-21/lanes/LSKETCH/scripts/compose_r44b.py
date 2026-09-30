#!$T90_PY
"""Round 44b compose: Main_Fig1 (band a r44 + band 1f r44b), Main_Fig6 (+ALT) from the r44b sketches, Main_Fig3 / Main_Fig4 through the
round-40 LMAIN composer and 600-dpi raster twin (--tag r44), Main_Fig5 through the round-42 LFIG5b composer (patched band sha) and its
raster twin, Supp_Fig16 through the round-40 LSUPP route (flatten + text redaction + band placement) with the r44 band. Deliverables
under lanes/LSKETCH/<Sheet>/ with verify/checks.txt. Renders: Ghostscript only (the LMAIN/LFIG5b/LSUPP helpers run their own watchdogs)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import gc, hashlib, json, os, shutil, subprocess, sys, time
import numpy as np
from PIL import Image
import fitz
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND40_2026-09-18/lanes/LMAIN/scripts")
import lmain_lib as L
T90 = paths.FIGURE_ROOT; R44 = f"{paths.FIGURE_ROOT}/ROUND44_2026-09-21"; LANE = f"{R44}/lanes/LSKETCH"; BUILD = f"{LANE}/work/build"
R40M = f"{paths.FIGURE_ROOT}/ROUND40_2026-09-18/lanes/LMAIN"; R42F = f"{paths.FIGURE_ROOT}/ROUND42_2026-09-19/lanes/LFIG5b"; V21 = f"{R44}/figures/NEW_FINAL_SET_V21"
R37W = f"{R40M}/Main_Fig1/work/r37"; GS = paths.GS; PY = paths.PY
Image.MAX_IMAGE_PIXELS = None
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
def render(pdf, png, dpi=150):
    r = subprocess.run([GS, "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=png16m", f"-r{dpi}", f"-sOutputFile={png}", pdf], capture_output=True, timeout=1800); return r.returncode == 0 and os.path.exists(png)
def pagesize(pdf):
    from pypdf import PdfReader
    b = PdfReader(pdf).pages[0].mediabox; return float(b.width), float(b.height)
def finish_checks(sheet, checks, extra=None):
    sd = f"{LANE}/{sheet}"; os.makedirs(f"{sd}/verify", exist_ok=True); os.makedirs(f"{sd}/work", exist_ok=True)
    ok = all(c.startswith("PASS") for c in checks); open(f"{sd}/verify/checks.txt", "w").write("\n".join(checks) + ("\nRESULT ALL PASS\n" if ok else "\nRESULT FAIL\n"))
    if extra: json.dump(extra, open(f"{sd}/work/compose_log.json", "w"), indent=1, default=str)
    print(sheet, "ALL PASS" if ok else "FAIL"); return ok
def run(cmd, log, timeout=3600):
    r = subprocess.run(cmd, capture_output=True, timeout=timeout); open(log, "wb").write(r.stdout + r.stderr); return r.returncode, r.stdout.decode("utf8", "replace")[-600:]
def same_size_check(sheet, out, checks):
    w, h = pagesize(out); w0, h0 = pagesize(f"{V21}/{sheet}.pdf"); checks.append(("PASS " if abs(w - w0) < 0.05 and abs(h - h0) < 0.05 else "FAIL ") + f"page {w:.2f} x {h:.2f} = V21's {w0:.2f} x {h0:.2f}")

def fig1(BAND_A=f"{BUILD}/band_a_r44.pdf", BAND_F=f"{BUILD}/band_1f_r44b.pdf", PANELS=f"{R37W}/panels_bcde_m.pdf", panels_lane=None):
    SHEET = "Main_Fig1"; SD = f"{LANE}/{SHEET}"; os.makedirs(f"{SD}/work", exist_ok=True); OUT = f"{SD}/{SHEET}.pdf"; checks = []
    REC = json.load(open(f"{R37W}/compose_record.json")); BT = json.load(open(f"{R37W}/base_text.json")); W, H_OLD = BT["page_exact"]; LET = {l["text"]: l for l in BT["letters"]}
    if panels_lane: lc = open(f"{panels_lane}/verify/checks.txt").read().strip().splitlines(); checks.append(("PASS " if lc[-1] == "RESULT ALL PASS" and not any(l.startswith("FAIL") for l in lc) else "FAIL ") + f"panels page from {panels_lane} (six families), lane ALL PASS")
    else: checks.append(("PASS " if sha(PANELS) == REC["panels"]["sha256"] else "FAIL ") + "panels page equals the round-37 record")
    SLOT_A_BOTTOM, SLOT_F_TOP = 259.736, 948.868
    da = fitz.open(BAND_A); df_ = fitz.open(BAND_F); pan = fitz.open(PANELS); Ha, Hf = da[0].rect.height, df_[0].rect.height
    checks.append(("PASS " if abs(da[0].rect.width - W) < 0.05 and abs(df_[0].rect.width - W) < 0.05 else "FAIL ") + f"band widths equal the page ({W})")
    checks.append(("PASS " if abs(Ha - 231.114) < 0.05 else "FAIL ") + f"band a height {Ha:.3f} (= the delivered band)")
    H = SLOT_F_TOP + Hf; slot_a = fitz.Rect(0, SLOT_A_BOTTOM - Ha, W, SLOT_A_BOTTOM); slot_f = fitz.Rect(0, SLOT_F_TOP, W, H)
    out = fitz.open(); pg = out.new_page(width=W, height=H); pg.show_pdf_page(slot_a, da, 0); pg.show_pdf_page(fitz.Rect(0, slot_a.y1, W, slot_f.y0), pan, 0, clip=fitz.Rect(0, slot_a.y1, W, slot_f.y0)); pg.show_pdf_page(slot_f, df_, 0)
    font = fitz.Font(fontfile=L.ARIALB); tw = fitz.TextWriter(pg.rect); letters = []
    for ch in "abcdef": tw.append((LET[ch]["bbox"][0], LET[ch]["bbox"][1] + 11.75), ch, font=font, fontsize=13.0); letters.append((ch, LET[ch]["bbox"][0], round(LET[ch]["bbox"][1] + 11.75, 2)))
    tw.write_text(pg, color=(0, 0, 0)); TITLE_X = LET["a"]["bbox"][0]; tt = fitz.TextWriter(pg.rect); tt.append(fitz.Point(TITLE_X, L.TITLE_BASELINE), "Figure 1", font=font, fontsize=L.TITLE_SIZE); tt.write_text(pg, color=L.INK_RGB)
    out.save(OUT, deflate=True, garbage=3); out.close(); da.close(); df_.close(); pan.close(); gc.collect()
    n, w, h = L.gs_pdfinfo(OUT); checks.append(("PASS " if n == 1 and abs(w - W) < 0.01 and abs(h - H) < 0.01 else "FAIL ") + f"page {w} x {h} (band f {Hf:.3f} pt from {SLOT_F_TOP})")
    same_size_check(SHEET, OUT, checks); checks.append(("PASS " if render(OUT, f"{SD}/{SHEET}_150dpi.png") else "FAIL ") + "gs 150 dpi render")
    return finish_checks(SHEET, checks, dict(band_a=dict(path=BAND_A, sha256=sha(BAND_A)), band_f=dict(path=BAND_F, sha256=sha(BAND_F)), panels=dict(path=PANELS, sha256=sha(PANELS)), letters=letters, page=[w, h], out_sha256=sha(OUT), written=L.now()))

def fig6(build, SHEET):
    SD = f"{LANE}/{SHEET}"; os.makedirs(f"{SD}/work", exist_ok=True); OUT = f"{SD}/{SHEET}.pdf"; checks = []; TOP_PAD = 16.0; TITLE_XY = (36.0, L.TITLE_BASELINE)
    n0, w0, h0 = L.gs_pdfinfo(build); checks.append(("PASS " if n0 == 1 and abs(w0 - 968.94) < 0.05 else "FAIL ") + f"build page {w0} x {h0}")
    d = fitz.open(build); p = d[0]; w, h = p.rect.width, p.rect.height; mb = p.mediabox; new = fitz.Rect(mb.x0, mb.y0, mb.x1, mb.y1 + TOP_PAD); p.set_mediabox(new)
    if abs(p.rect.height - (h + TOP_PAD)) > 0.01: p.set_cropbox(new)
    p.wrap_contents(); tw = fitz.TextWriter(p.rect); tw.append(fitz.Point(*TITLE_XY), "Figure 6", font=fitz.Font(fontfile=L.ARIALB), fontsize=L.TITLE_SIZE); tw.write_text(p, color=L.INK_RGB)
    tmp = OUT + ".part"; d.save(tmp, deflate=True, garbage=3); d.close(); os.replace(tmp, OUT); gc.collect()
    n, W, H = L.gs_pdfinfo(OUT); checks.append(("PASS " if n == 1 and abs(W - w0) < 0.01 and abs(H - (h0 + TOP_PAD)) < 0.01 else "FAIL ") + f"page {W} x {H} = build + 16 pt strip, title 'Figure 6'")
    if SHEET == "Main_Fig6": same_size_check(SHEET, OUT, checks)
    checks.append(("PASS " if render(OUT, f"{SD}/{SHEET}_150dpi.png") else "FAIL ") + "gs 150 dpi render")
    return finish_checks(SHEET, checks, dict(build=build, build_sha256=sha(build), page=[W, H], out_sha256=sha(OUT), written=L.now()))

def raster_sheet(SHEET, compose_cmd, raster_cmd, twin_dir, band):
    SD = f"{LANE}/{SHEET}"; os.makedirs(f"{SD}/work", exist_ok=True); checks = []
    rc, tail = run(compose_cmd, f"{SD}/work/compose.log"); checks.append(("PASS " if rc == 0 else "FAIL ") + f"vector compose rc {rc}" + ("" if rc == 0 else " " + tail[-300:]))
    if rc != 0: return finish_checks(SHEET, checks)
    t0 = time.time(); rc, tail = run(raster_cmd, f"{SD}/work/raster.log", timeout=5400); checks.append(("PASS " if rc == 0 else "FAIL ") + f"600 dpi raster twin rc {rc} in {time.time() - t0:.0f} s" + ("" if rc == 0 else " " + tail[-300:]))
    if rc != 0: return finish_checks(SHEET, checks)
    for ext in (".pdf", ".tif", "_150dpi.png"):
        src = f"{twin_dir}/{SHEET}_r44{ext}"; assert os.path.exists(src), src; shutil.copy2(src, f"{SD}/{SHEET}{ext}")
    shutil.copy2(f"{twin_dir}/raster_record_r44.json", f"{SD}/work/raster_record.json")
    OUT = f"{SD}/{SHEET}.pdf"; same_size_check(SHEET, OUT, checks)
    rr = json.load(open(f"{SD}/work/raster_record.json")); checks.append(("PASS " if rr.get("size_check_within_0.5pt") in (True, "True") and str(rr.get("dpi")) == "600" else "FAIL ") + f"raster record: {rr.get('width_px')} x {rr.get('height_px')} px at 600 dpi, gs {rr.get('gs_status')}")
    d = fitz.open(OUT); nimg = len(d[0].get_images()); ntext = 0; d.close(); checks.append(("PASS " if nimg == 1 else "FAIL ") + f"one image on the raster page")
    return finish_checks(SHEET, checks, dict(band=band, band_sha256=sha(band), vector=compose_cmd[-1] if "--out" not in compose_cmd else compose_cmd[compose_cmd.index("--out") + 1], out_sha256=sha(OUT), written=L.now()))

def supp16():
    SHEET = "Supp_Fig16"; SD = f"{LANE}/{SHEET}"; os.makedirs(f"{SD}/work", exist_ok=True); checks = []
    rc, tail = run([PY, f"{LANE}/scripts/build_s16_r44.py"], f"{SD}/work/build.log"); checks.append(("PASS " if rc == 0 else "FAIL ") + f"build (flatten V16 + text redaction + r44 band) rc {rc}" + ("" if rc == 0 else " " + tail[-300:]))
    if rc != 0: return finish_checks(SHEET, checks)
    OUT = f"{SD}/{SHEET}.pdf"; same_size_check(SHEET, OUT, checks)
    checks.append(("PASS " if render(OUT, f"{SD}/{SHEET}_150dpi.png") and render(f"{V21}/{SHEET}.pdf", f"{SD}/work/old_150dpi.png") and render(f"{BUILD}/band_s16_r44.pdf", f"{SD}/work/band_150dpi.png") else "FAIL ") + "gs 150 dpi renders (new, V21, band)")
    new = np.asarray(Image.open(f"{SD}/{SHEET}_150dpi.png").convert("RGB")); old = np.asarray(Image.open(f"{SD}/work/old_150dpi.png").convert("RGB")); band = np.asarray(Image.open(f"{SD}/work/band_150dpi.png").convert("RGB"))
    y0 = int(round(80 * 150 / 72)); checks.append(("PASS " if new.shape == old.shape and np.array_equal(new[y0:], old[y0:]) else "FAIL ") + f"below the band (y >= 80 pt) pixel-identical to V21 ({new.shape[1]}x{new.shape[0]})")
    yb = min(band.shape[0], int(round(79 * 150 / 72))); diff = (np.abs(new[:yb, :band.shape[1]].astype(int) - band[:yb].astype(int)).max(axis=2) > 64).mean()
    checks.append(("PASS " if diff < 0.01 else "FAIL ") + f"band region equals the r44 band's own render up to a sub-pixel render offset (strongly differing pixel share {diff:.4f}, limit 0.01)")
    crop = Image.fromarray(new[:yb]); crop = crop.resize((crop.width * 3, crop.height * 3), Image.LANCZOS); crop.save(f"{SD}/work/ocr_band.png")
    t = subprocess.run([paths.TESSERACT, f"{SD}/work/ocr_band.png", "-", "--psm", "11"], capture_output=True).stdout.decode("utf8", "replace")
    import difflib, re
    norm = lambda z: re.sub(r"[^a-z0-9%]", "", z.lower()); tn = norm(t); want = ["Does the oxygen signal predict", "future disease at home", "SHHS + MrOS", "5,802 + 2,911 people", "heart attack", "or failure", "death"]   # "years later" and "90% SpO2" sit beside arrows and tesseract skips them (seen on the crop)
    def found(w):
        k = norm(w); best = max((difflib.SequenceMatcher(None, k, tn[i:i + len(k)]).ratio() for i in range(max(1, len(tn) - len(k) + 1))), default=0); return best >= 0.8
    miss = [w for w in want if not found(w)]; checks.append(("PASS " if not miss else "FAIL ") + f"OCR reads the band's texts ({len(want) - len(miss)} of {len(want)}, fuzzy 0.8; missing {miss})")
    return finish_checks(SHEET, checks, dict(band=f"{BUILD}/band_s16_r44.pdf", band_sha256=sha(f"{BUILD}/band_s16_r44.pdf"), out_sha256=sha(OUT), written=L.now()))

if __name__ == "__main__":
    if "--supp16-only" in sys.argv: sys.exit(0 if supp16() else 1)
    if "--fig1-44c" in sys.argv: sys.exit(0 if fig1(BAND_A=f"{BUILD}/band_a_r44c.pdf", PANELS=f"{R44}/lanes/LFIG1B/Main_Fig1/work/panels_bcde_m.pdf", panels_lane=f"{R44}/lanes/LFIG1B/Main_Fig1") else 1)
    ok = fig1() & fig6(f"{BUILD}/Main_Fig6_r44b.pdf", "Main_Fig6") & fig6(f"{BUILD}/Main_Fig6_tri_r44b.pdf", "Main_Fig6_ALT")
    ok &= raster_sheet("Main_Fig3", [PY, f"{R40M}/Main_Fig3/scripts/02_compose_r40.py", "--band", f"{BUILD}/band_3d_r44.pdf", "--out", f"{LANE}/Main_Fig3/work/Main_Fig3_r44_vector.pdf"],
                       [PY, f"{R40M}/scripts/rasterize_r40.py", "Main_Fig3", "600", "--vector", f"{LANE}/Main_Fig3/work/Main_Fig3_r44_vector.pdf", "--tag", "r44"], f"{R40M}/Main_Fig3/work", f"{BUILD}/band_3d_r44.pdf")
    ok &= raster_sheet("Main_Fig4", [PY, f"{R40M}/Main_Fig4/scripts/02_compose_r40.py", "--band", f"{BUILD}/band_4e_r44.pdf", "--out", f"{LANE}/Main_Fig4/work/Main_Fig4_r44_vector.pdf"],
                       [PY, f"{R40M}/scripts/rasterize_r40.py", "Main_Fig4", "600", "--vector", f"{LANE}/Main_Fig4/work/Main_Fig4_r44_vector.pdf", "--tag", "r44"], f"{R40M}/Main_Fig4/work", f"{BUILD}/band_4e_r44.pdf")
    ok &= raster_sheet("Main_Fig5", [PY, f"{LANE}/scripts/compose_fig5_r44.py", "--band", f"{BUILD}/band_5d_r44.pdf", "--out", f"{LANE}/Main_Fig5/work/Main_Fig5_r44_vector.pdf"],
                       [PY, f"{R42F}/scripts/rasterize_r42b.py", "Main_Fig5", "600", "--vector", f"{LANE}/Main_Fig5/work/Main_Fig5_r44_vector.pdf", "--tag", "r44"], f"{R42F}/Main_Fig5/work", f"{BUILD}/band_5d_r44.pdf")
    ok &= supp16()
    sys.exit(0 if ok else 1)
