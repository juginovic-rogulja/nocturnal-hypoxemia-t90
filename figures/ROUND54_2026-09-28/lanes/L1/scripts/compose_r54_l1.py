#!$T90_PY
"""Round 54 lane L1 compose: Main_Fig1 = band a (r54) + the L1 panels page (b, c, d, e at 14 / 15 / 13 pt) + band 1f (r54), letters and the title strip
at 18 pt, through the round-44b compose logic (compose_r44b.fig1) with two round-54 differences: (1) the bands' heights are read from their page boxes:
band a keeps its bottom on the panels' top (259.736) and a taller band rises into the title strip as long as 16 pt stay above it for the 18-pt title
(its caps end at y 14); only if that fails is everything below the strip shifted down by the rest; band f sits under the panels, so the composed height
is 948.868 + shift + Hf; (2) the panel letters c, d, e follow their re-laid panels' left edges (c 483.6, d 11.17, e 473.4, V30: 496.64, 14.17, 496.69).
Usage: compose_r54_l1.py [--r49-bands]   (the round-49 bands as a provisional stand-in when BANDS_READY.txt is missing)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import gc, hashlib, json, os, subprocess, sys, time
import fitz
T90 = paths.FIGURE_ROOT; R54 = f"{paths.FIGURE_ROOT}/ROUND54_2026-09-28"; LANE = f"{R54}/lanes/L1"
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND40_2026-09-18/lanes/LMAIN/scripts"); import lmain_lib as L
R37W = f"{paths.FIGURE_ROOT}/ROUND40_2026-09-18/lanes/LMAIN/Main_Fig1/work/r37"; GS = paths.GS
SHEET = "Main_Fig1"; SD = f"{LANE}/{SHEET}"; os.makedirs(f"{SD}/work", exist_ok=True); os.makedirs(f"{SD}/verify", exist_ok=True)
PROVISIONAL = "--r49-bands" in sys.argv
BUILD = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LSKETCH/work/build" if PROVISIONAL else f"{R54}/lanes/LSKETCH/work/build"
BAND_A = f"{BUILD}/band_a_r49.pdf" if PROVISIONAL else f"{BUILD}/band_a_r54.pdf"
BAND_F = f"{BUILD}/band_1f_r49.pdf" if PROVISIONAL else f"{BUILD}/band_1f_r54.pdf"
PANELS = f"{SD}/work/panels_bcde_m.pdf"; OUT = f"{SD}/{SHEET}.pdf"
LETTER_PT = 18.0; TITLE_PT = 18.0; TITLE_BASELINE = 18.0   # round 54 (coordinator, consistency with Figures 5 and 6): the 18-pt title on baseline 18, caps from y 5.1 (LMAIN: 14.0)
SLOT_A_BOTTOM, SLOT_F_TOP, HA_V30 = 259.736, 948.868, 231.114
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
def hydrated(p):
    st = os.stat(p); assert not (st.st_size > 0 and st.st_blocks == 0), f"EVICTED (0 blocks): {p}"; return p
for p in (BAND_A, BAND_F, PANELS): hydrated(p)
if not PROVISIONAL: assert os.path.exists(f"{R54}/lanes/LSKETCH/BANDS_READY.txt"), "BANDS_READY.txt missing: compose with --r49-bands for a provisional sheet"
BT = json.load(open(f"{R37W}/base_text.json")); W, H_PANELS_PAGE = BT["page_exact"]; LET = {l["text"]: l for l in BT["letters"]}
LETTER_X = {"a": LET["a"]["bbox"][0], "b": LET["b"]["bbox"][0], "c": 483.6, "d": 11.17, "e": 473.4, "f": LET["f"]["bbox"][0]}   # round 54: c, d, e follow the re-laid panels
da = fitz.open(BAND_A); df_ = fitz.open(BAND_F); pan = fitz.open(PANELS); Ha, Hf = da[0].rect.height, df_[0].rect.height
assert abs(da[0].rect.width - W) < 0.05 and abs(df_[0].rect.width - W) < 0.05 and abs(pan[0].rect.width - W) < 0.05, (da[0].rect.width, df_[0].rect.width, pan[0].rect.width)
STRIP_MIN = TITLE_BASELINE + 0.21 * TITLE_PT + 1.0   # the title's descender (the g of Figure) plus 1 pt must stay above the band
shift = max(0.0, Ha - (SLOT_A_BOTTOM - STRIP_MIN))  # only a band taller than 243.736 pushes the panels and band f down
slot_a = fitz.Rect(0, SLOT_A_BOTTOM + shift - Ha, W, SLOT_A_BOTTOM + shift); assert slot_a.y0 >= STRIP_MIN - 1e-6, slot_a
panels_rect = fitz.Rect(0, SLOT_A_BOTTOM + shift, W, SLOT_F_TOP + shift); panels_clip = fitz.Rect(0, SLOT_A_BOTTOM, W, SLOT_F_TOP)
slot_f = fitz.Rect(0, SLOT_F_TOP + shift, W, SLOT_F_TOP + shift + Hf); H = slot_f.y1
out = fitz.open(); pg = out.new_page(width=W, height=H)
pg.show_pdf_page(slot_a, da, 0); pg.show_pdf_page(panels_rect, pan, 0, clip=panels_clip); pg.show_pdf_page(slot_f, df_, 0)
font = fitz.Font(fontfile=L.ARIALB); tw = fitz.TextWriter(pg.rect); letters = []
for ch in "abcdef":
    top = LET[ch]["bbox"][1] + (shift if ch != "a" else 0.0); base = top + 11.75   # the compose pattern's baseline offset, kept (the 18-pt caps sit 1.15 pt higher than the 14-pt ones)
    tw.append((LETTER_X[ch], base), ch, font=font, fontsize=LETTER_PT); letters.append(dict(ch=ch, x=LETTER_X[ch], baseline=round(base, 3), size=LETTER_PT))
tw.write_text(pg, color=(0, 0, 0))
tt = fitz.TextWriter(pg.rect); tt.append(fitz.Point(LETTER_X["a"], TITLE_BASELINE), "Figure 1", font=font, fontsize=TITLE_PT); tt.write_text(pg, color=L.INK_RGB)
title_w = font.text_length("Figure 1", fontsize=TITLE_PT); title_bottom = TITLE_BASELINE + 0.21 * TITLE_PT
assert title_bottom < slot_a.y0, (title_bottom, slot_a.y0)   # the title (with its descender) ends above the band a slot
tmp = OUT + ".part"; out.save(tmp, deflate=True, garbage=3); out.close(); da.close(); df_.close(); pan.close(); gc.collect(); os.replace(tmp, OUT)
n, w, h = L.gs_pdfinfo(OUT); assert n == 1 and abs(w - W) < 0.01 and abs(h - H) < 0.01, (n, w, h, W, H)
r = subprocess.run([GS, "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=png16m", "-r150", f"-sOutputFile={SD}/{SHEET}_150dpi.png", OUT], capture_output=True, timeout=1800)
assert r.returncode == 0 and os.path.getsize(f"{SD}/{SHEET}_150dpi.png") > 0, r.stderr[-300:]
log = dict(provisional_r49_bands=PROVISIONAL, band_a=dict(path=BAND_A, sha256=sha(BAND_A), height=round(Ha, 3)), band_f=dict(path=BAND_F, sha256=sha(BAND_F), height=round(Hf, 3)),
           panels=dict(path=PANELS, sha256=sha(PANELS), clip=[0, SLOT_A_BOTTOM, W, SLOT_F_TOP], placed=[0, panels_rect.y0, W, panels_rect.y1]), shift=round(shift, 3),
           slots=dict(a=[slot_a.y0, slot_a.y1], panels=[panels_rect.y0, panels_rect.y1], f=[slot_f.y0, slot_f.y1]), letters=letters, title=dict(text="Figure 1", x=LETTER_X["a"], baseline=TITLE_BASELINE, size=TITLE_PT),
           title_clearance=dict(title_x=[LETTER_X["a"], round(LETTER_X["a"] + title_w, 2)], cap_top=round(TITLE_BASELINE - 0.716 * TITLE_PT, 2), descender_bottom=round(title_bottom, 2), band_a_top=round(slot_a.y0, 3)),
           page=[round(w, 3), round(h, 3)], out_sha256=sha(OUT), written=time.strftime("%Y-%m-%d %H:%M:%S"))
json.dump(log, open(f"{SD}/work/compose_log.json", "w"), indent=1)
print(f"composed {OUT}: page {w} x {h} (band a {Ha:.3f}, shift {shift:.3f}, band f {Hf:.3f}), provisional={PROVISIONAL}")
