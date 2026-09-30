#!$T90_PY
"""Round 44 compose: Main_Fig1 (round-40 idiom: band a = the LFIG1A band with the organ grid 20% smaller, panels b to e 1:1 from the
round-37 page, band f = band_1f_r44 (one line, 46 mm), page height = f slot top + band height; title and letters as round 40) and
Main_Fig6 (= the r44 sketch + 16 pt title strip, "Figure 6" 13 pt at (36, 14)); the triangle variant composed the same way as
Main_Fig6_ALT. Placement only (show_pdf_page, TextWriter); renders with Ghostscript; checks.txt per sheet."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import gc, hashlib, json, os, subprocess, sys, time
import fitz
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND40_2026-09-18/lanes/LMAIN/scripts")
import lmain_lib as L
R44 = f"{paths.FIGURE_ROOT}/ROUND44_2026-09-21"; LANE = f"{R44}/lanes/LSKETCH"; BUILD = f"{LANE}/work/build"
R37W = f"{L.LANE}/Main_Fig1/work/r37"; GS = paths.GS
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
def render(pdf, png, dpi=150):
    r = subprocess.run([GS, "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=png16m", f"-r{dpi}", f"-sOutputFile={png}", pdf], capture_output=True, timeout=900); return r.returncode == 0 and os.path.exists(png)
def write_checks(sd, checks):
    ok = all(c.startswith("PASS") for c in checks); open(f"{sd}/verify/checks.txt", "w").write("\n".join(checks) + ("\nRESULT ALL PASS\n" if ok else "\nRESULT FAIL\n")); return ok

def fig1():
    SHEET = "Main_Fig1"; SD = f"{LANE}/{SHEET}"; OUT = f"{SD}/{SHEET}.pdf"; checks = []
    BAND_A = f"{R44}/lanes/LFIG1A/band_a_r44.pdf"; BAND_F = f"{BUILD}/band_1f_r44.pdf"; PANELS = f"{R37W}/panels_bcde_m.pdf"
    REC = json.load(open(f"{R37W}/compose_record.json")); BT = json.load(open(f"{R37W}/base_text.json")); W, H_OLD = BT["page_exact"]
    LET = {l["text"]: l for l in BT["letters"]}; assert set(LET) == set("abcdef")
    checks.append(("PASS " if sha(PANELS) == REC["panels"]["sha256"] else "FAIL ") + "panels page equals the round-37 record")
    lane_ok = open(f"{R44}/lanes/LFIG1A/verify/checks.txt").read().strip().endswith("RESULT ALL PASS"); checks.append(("PASS " if lane_ok else "FAIL ") + "band a lane (organ grid 0.8) ALL PASS")
    SLOT_A_BOTTOM, SLOT_F_TOP = 259.736, 948.868
    da = fitz.open(BAND_A); df_ = fitz.open(BAND_F); pan = fitz.open(PANELS)
    Ha, Hf = da[0].rect.height, df_[0].rect.height
    checks.append(("PASS " if abs(da[0].rect.width - W) < 0.01 and abs(df_[0].rect.width - W) < 0.05 and abs(pan[0].rect.width - W) < 0.01 else "FAIL ") + f"band widths equal the page ({W})")
    checks.append(("PASS " if abs(Ha - 231.114) < 0.01 else "FAIL ") + f"band a height {Ha:.3f} (= the round-37 band)")
    H = SLOT_F_TOP + Hf
    slot_a = fitz.Rect(0, SLOT_A_BOTTOM - Ha, W, SLOT_A_BOTTOM); slot_f = fitz.Rect(0, SLOT_F_TOP, W, H)
    out = fitz.open(); pg = out.new_page(width=W, height=H)
    pg.show_pdf_page(slot_a, da, 0)
    pg.show_pdf_page(fitz.Rect(0, slot_a.y1, W, slot_f.y0), pan, 0, clip=fitz.Rect(0, slot_a.y1, W, slot_f.y0))
    pg.show_pdf_page(slot_f, df_, 0)
    font = fitz.Font(fontfile=L.ARIALB); tw = fitz.TextWriter(pg.rect); letters = []
    for ch in "abcdef":
        tw.append((LET[ch]["bbox"][0], LET[ch]["bbox"][1] + 11.75), ch, font=font, fontsize=13.0); letters.append((ch, LET[ch]["bbox"][0], round(LET[ch]["bbox"][1] + 11.75, 2)))
    tw.write_text(pg, color=(0, 0, 0))
    TITLE_X = LET["a"]["bbox"][0]; tt = fitz.TextWriter(pg.rect); tt.append(fitz.Point(TITLE_X, L.TITLE_BASELINE), "Figure 1", font=font, fontsize=L.TITLE_SIZE); tt.write_text(pg, color=L.INK_RGB)
    out.save(OUT, deflate=True, garbage=3); out.close(); da.close(); df_.close(); pan.close(); gc.collect()
    n, w, h = L.gs_pdfinfo(OUT); checks.append(("PASS " if n == 1 and abs(w - W) < 0.01 and abs(h - H) < 0.01 else "FAIL ") + f"page {w} x {h} (was {W} x {H_OLD}; band f {Hf:.3f} pt in slot from {SLOT_F_TOP})")
    checks.append(("PASS " if LET["f"]["bbox"][1] + 11.75 > SLOT_F_TOP and LET["f"]["bbox"][1] + 11.75 < SLOT_F_TOP + 30 else "FAIL ") + f"letter f baseline {LET['f']['bbox'][1] + 11.75:.2f} sits in the band f slot")
    checks.append(("PASS " if render(OUT, f"{SD}/{SHEET}_150dpi.png") else "FAIL ") + "gs 150 dpi render")
    ok = write_checks(SD, checks)
    json.dump(dict(sheet=SHEET, band_a=dict(path=BAND_A, sha256=sha(BAND_A), slot=list(slot_a)), panels=dict(path=PANELS, sha256=sha(PANELS)), band_f=dict(path=BAND_F, sha256=sha(BAND_F), slot=list(slot_f)),
                   title=dict(text="Figure 1", origin=[TITLE_X, L.TITLE_BASELINE], size=L.TITLE_SIZE), letters=letters, page=[w, h], page_before=[W, H_OLD], out=OUT, out_sha256=sha(OUT), written=L.now()), open(f"{SD}/work/compose_log.json", "w"), indent=1)
    print(SHEET, "ALL PASS" if ok else "FAIL", f"page {w} x {h}"); return ok

def fig6(build, SHEET):
    SD = f"{LANE}/{SHEET}"; OUT = f"{SD}/{SHEET}.pdf"; checks = []; TOP_PAD = 16.0; TITLE_XY = (36.0, L.TITLE_BASELINE)
    n0, w0, h0 = L.gs_pdfinfo(build); checks.append(("PASS " if n0 == 1 and abs(w0 - 968.94) < 0.05 else "FAIL ") + f"build page {w0} x {h0}")
    d = fitz.open(build); p = d[0]; w, h = p.rect.width, p.rect.height; mb = p.mediabox
    new = fitz.Rect(mb.x0, mb.y0, mb.x1, mb.y1 + TOP_PAD); p.set_mediabox(new)
    if abs(p.rect.height - (h + TOP_PAD)) > 0.01: p.set_cropbox(new)
    assert abs(p.rect.height - (h + TOP_PAD)) < 0.01 and abs(p.rect.width - w) < 0.01, (p.rect, w, h)
    p.wrap_contents()
    tw = fitz.TextWriter(p.rect); tw.append(fitz.Point(*TITLE_XY), "Figure 6", font=fitz.Font(fontfile=L.ARIALB), fontsize=L.TITLE_SIZE); tw.write_text(p, color=L.INK_RGB)
    tmp = OUT + ".part"; d.save(tmp, deflate=True, garbage=3); d.close(); os.replace(tmp, OUT); gc.collect()
    n, W, H = L.gs_pdfinfo(OUT); checks.append(("PASS " if n == 1 and abs(W - w0) < 0.01 and abs(H - (h0 + TOP_PAD)) < 0.01 else "FAIL ") + f"page {W} x {H} = build + 16 pt strip, title 'Figure 6' 13 pt at {TITLE_XY}")
    checks.append(("PASS " if render(OUT, f"{SD}/{SHEET}_150dpi.png") else "FAIL ") + "gs 150 dpi render")
    ok = write_checks(SD, checks)
    json.dump(dict(sheet=SHEET, build=build, build_sha256=sha(build), page=[W, H], title=dict(text="Figure 6", origin=list(TITLE_XY), size=L.TITLE_SIZE), out=OUT, out_sha256=sha(OUT), written=L.now()), open(f"{SD}/work/compose_log.json", "w"), indent=1)
    print(SHEET, "ALL PASS" if ok else "FAIL", f"page {W} x {H}"); return ok

if __name__ == "__main__":
    ok = fig1() & fig6(f"{BUILD}/Main_Fig6_r44.pdf", "Main_Fig6") & fig6(f"{BUILD}/Main_Fig6_tri_r44.pdf", "Main_Fig6_ALT")
    sys.exit(0 if ok else 1)
