#!/usr/bin/env python3
"""Round 54 figure assembly (2026-09-29), V31 = V30 with the six main figures rebuilt at the round-54 text sizes (ticks 14, axis titles 15, notes 13, letters 18; Figure 1 with the family key inside panel b; Figures 3, 4, 5 as 600-dpi rasters of the lanes' vector composes). usage: --dry-run | --run"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../..")))  # repository root, where paths.py lives
import paths
import hashlib, json, os, shutil, subprocess, sys, time

T90 = paths.FIGURE_ROOT
R44 = f"{paths.FIGURE_ROOT}/ROUND44_2026-09-21"; R46 = f"{paths.FIGURE_ROOT}/ROUND46_2026-09-24"; R47 = f"{paths.FIGURE_ROOT}/ROUND47_2026-09-24"; R49 = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26"; R50 = f"{paths.FIGURE_ROOT}/ROUND50_2026-09-26"; R51 = f"{paths.FIGURE_ROOT}/ROUND51_2026-09-26"; R52 = f"{paths.FIGURE_ROOT}/ROUND52_2026-09-26"; R54 = f"{paths.FIGURE_ROOT}/ROUND54_2026-09-28"
FIG = f"{R54}/figures"
BASE = f"{R52}/figures/NEW_FINAL_SET_V30"; PNGBASE = f"{R52}/figures/word_png_v30"; BASE_SHA = f"{R52}/figures/V30_SHA256.txt"
LANES = {"L1": ["Main_Fig1"], "L2": ["Main_Fig2"], "LRASTER": ["Main_Fig3", "Main_Fig4", "Main_Fig5"], "LSKETCH": ["Main_Fig6"], "LED7C": ["ED_Fig07"]}
OUT = f"{FIG}/NEW_FINAL_SET_V31"; SEP = f"{FIG}/T90 Figures V31 - separate PDFs"; REVIEW = f"{FIG}/T90 All Figures V31.pdf"
DOCX = f"{FIG}/T90 All Figures V31.docx"; PNG = f"{FIG}/word_png_v31"
MAKE_DOCX = f"{R44}/make_docx_r44.js"; VALIDATE = paths.DOCX_VALIDATE; GS = paths.GS
DL = os.path.expanduser("~/Downloads")
for _p in (OUT, SEP, REVIEW, DOCX, PNG): assert not os.path.abspath(_p).startswith(DL), _p
ORDER = [f"Main_Fig{i}" for i in range(1, 7)] + [f"ED_Fig{i:02d}" for i in range(1, 11)] + [f"Supp_Fig{i:02d}" for i in range(1, 22)]
BOX = (36.0, 86.1, 559.3, 769.8); A4 = (595.28, 841.89); PRINT = 1.0

def label(nm):
    return (f"Figure {int(nm[8:])}" if nm.startswith("Main") else f"Extended Data Fig. {int(nm[6:])}" if nm.startswith("ED") else f"Supplementary Fig. {int(nm[8:])}")
def plain(nm):
    return (f"Figure {int(nm[8:])}.pdf" if nm.startswith("Main") else f"Extended Data Fig {int(nm[6:])}.pdf" if nm.startswith("ED") else f"Supplementary Fig {int(nm[8:])}.pdf")
def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""): h.update(c)
    return h.hexdigest()
def base_hashes():
    out = {}
    for line in open(BASE_SHA):
        parts = line.strip().split("  ", 1)
        if len(parts) == 2: out[parts[1]] = parts[0]
    return out
def rebuilt_sheets():
    got = {}
    for lane, names in LANES.items():
        for name in names:
            p = f"{R54}/lanes/{lane}/{name}/{name}.pdf"; chk = f"{R54}/lanes/{lane}/{name}/verify/checks.txt"
            assert os.path.exists(p) and os.path.exists(chk), f"{name}: lane output or verify missing"
            lines = [l.strip() for l in open(chk) if l.strip()]
            assert lines and lines[-1].startswith("RESULT ALL PASS") and not any(l.startswith("FAIL") for l in lines), f"{name}: verifier not ALL PASS"
            got[name] = p
    return got
def run_watchdog(cmd, log, minutes, rss_kb=6_000_000):
    p = subprocess.Popen(cmd, stdout=open(log, "w"), stderr=subprocess.STDOUT); t0 = time.time(); peak = 0
    while p.poll() is None:
        time.sleep(2)
        try: rss = int(subprocess.run(["ps", "-o", "rss=", "-p", str(p.pid)], capture_output=True, text=True).stdout.strip() or 0)
        except Exception: rss = 0
        peak = max(peak, rss)
        if rss > rss_kb or time.time() - t0 > minutes * 60:
            p.kill(); raise RuntimeError(f"watchdog killed {cmd[:2]} (rss {rss} KB, {time.time() - t0:.0f} s)")
    return p.returncode, peak, time.time() - t0

REVIEW_CODE = r'''
import fitz, json, sys, gc
out, order_json, box_json, a4_json, set_dir, labels_json, print_s = sys.argv[1:8]
order = json.loads(order_json); box = json.loads(box_json); a4 = json.loads(a4_json); labels = json.loads(labels_json); PRINT = float(print_s)
bw, bh = box[2] - box[0], box[3] - box[1]
rv = fitz.open(); rec = {}
for name in order:
    sh = fitz.open(f"{set_dir}/{name}.pdf"); w, h = sh[0].rect.width, sh[0].rect.height
    s = min(PRINT, bw / w, bh / h)   # round 43: print size, reduced only to fit
    dw, dh = w * s, h * s; x0 = box[0] + (bw - dw) / 2; y0 = box[1]
    pg = rv.new_page(width=a4[0], height=a4[1])
    pg.show_pdf_page(fitz.Rect(x0, y0, x0 + dw, y0 + dh), sh, 0)
    pg.insert_text((36.0, 824.6), labels[name], fontsize=8, fontname="helv", color=(0, 0, 0))
    rec[name] = dict(sheet_pt=[round(w, 2), round(h, 2)], scale=round(s, 4), at_print_size=abs(s - PRINT) < 1e-9, placed=[round(x0, 2), round(y0, 2), round(dw, 2), round(dh, 2)],
                     width_pct=round(100 * dw / bw), height_pct=round(100 * dh / bh))
    sh.close(); gc.collect()
rv.save(out, deflate=True, garbage=3); rv.close()
json.dump(dict(rule="print size, reduced to fit the box", print_scale=PRINT, box=box, sheets=rec), open(out + ".scale.json", "w"), indent=1)
print("ok")
'''

def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "--dry-run"
    hb = base_hashes(); rebuilt = rebuilt_sheets()
    plan = []
    for name in ORDER:
        src = rebuilt.get(name, f"{BASE}/{name}.pdf"); assert os.path.exists(src), src
        if name not in rebuilt: assert sha(src) == hb[f"NEW_FINAL_SET_V30/{name}.pdf"], f"{name}: V30 base sheet differs from V30_SHA256.txt"
        plan.append((name, src, "ROUND54 larger text" if name in rebuilt else "V30 unchanged"))
    print("rebuilt:", [n for n, _, w in plan if w != "V30 unchanged"], "| V30 unchanged:", sum(w == "V30 unchanged" for _, _, w in plan))
    if mode != "--run": return
    os.makedirs(FIG, exist_ok=True); os.makedirs(f"{R54}/logs", exist_ok=True)
    for d in (OUT, SEP, PNG):
        if os.path.isdir(d): shutil.rmtree(d)
        os.makedirs(d)
    manifest = {}
    for name, src, why in plan:
        dst = f"{OUT}/{name}.pdf"; shutil.copy2(src, dst); shutil.copy2(src, f"{SEP}/{plain(name)}")
        manifest[name] = {"source": src, "why": why, "sha256": sha(dst)}
    args = [paths.PY, "-c", REVIEW_CODE, REVIEW + ".tmp", json.dumps(ORDER), json.dumps(BOX), json.dumps(A4), OUT, json.dumps({n: label(n) for n in ORDER}), str(PRINT)]
    rc, peak, secs = run_watchdog(args, f"{R54}/logs/compose_v31_review.log", minutes=30)
    assert rc == 0, open(f"{R54}/logs/compose_v31_review.log").read()[-500:]
    os.replace(REVIEW + ".tmp", REVIEW); os.replace(REVIEW + ".tmp.scale.json", REVIEW + ".scale.json")
    from pypdf import PdfReader
    n_pages = len(PdfReader(REVIEW).pages); assert n_pages == 37, n_pages
    print(f"review PDF: 37 pages, {os.path.getsize(REVIEW) / 1e6:.1f} MB, child peak RSS {peak / 1024:.0f} MB, {secs:.0f} s")
    pngs = {}
    for name, src, why in plan:
        dst = f"{PNG}/{name}.png"
        if why == "V30 unchanged":
            shutil.copy2(f"{PNGBASE}/{name}.png", dst); assert sha(dst) == hb[f"word_png_v30/{name}.png"], name
            pngs[name] = {"from": "word_png_v30", "sha256": sha(dst)}
        else:
            rc, peak, secs = run_watchdog([GS, "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-dFirstPage=1", "-dLastPage=1", "-sDEVICE=png16m", "-r150", f"-sOutputFile={dst}", src], f"{R54}/logs/compose_v31_gs_{name}.log", minutes=30)
            assert rc == 0 and os.path.exists(dst), name
            pngs[name] = {"from": "gs 150 dpi", "sha256": sha(dst), "peak_rss_kb": peak, "secs": round(secs)}
    if os.path.exists(DOCX): os.remove(DOCX)
    r = subprocess.run(["node", MAKE_DOCX, PNG, DOCX], capture_output=True, text=True, cwd=os.path.dirname(MAKE_DOCX), timeout=600)
    assert r.returncode == 0 and os.path.exists(DOCX), r.stdout + r.stderr
    v = subprocess.run([paths.PY, VALIDATE, DOCX], capture_output=True, text=True, timeout=600)
    vline = (v.stdout.strip().splitlines() or ["(no output)"])[-1]
    print(f"docx: {os.path.getsize(DOCX) / 1e6:.1f} MB, validate: {vline} (rc {v.returncode})")
    lines = []
    for name in ORDER:
        lines.append(f"{manifest[name]['sha256']}  NEW_FINAL_SET_V31/{name}.pdf")
        lines.append(f"{sha(f'{SEP}/{plain(name)}')}  T90 Figures V31 - separate PDFs/{plain(name)}")
        lines.append(f"{pngs[name]['sha256']}  word_png_v31/{name}.png")
    lines.append(f"{sha(REVIEW)}  T90 All Figures V31.pdf"); lines.append(f"{sha(DOCX)}  T90 All Figures V31.docx")
    open(f"{FIG}/V31_SHA256.txt", "w").write("\n".join(lines) + "\n")
    json.dump({"sheets": manifest, "pngs": pngs, "placement_rule": "print size (scale 1.0), reduced only where the sheet would not fit the box; the Word file the same at 96/150 of the 150-dpi render",
               "time": time.strftime("%Y-%m-%d %H:%M:%S"), "script": os.path.abspath(__file__), "script_sha256": sha(__file__), "v30_set": BASE, "lane_outputs": rebuilt},
              open(f"{FIG}/V31_ASSEMBLY_MANIFEST.json", "w"), indent=1)
    print("V31 set assembled:", OUT)

if __name__ == "__main__": main()
