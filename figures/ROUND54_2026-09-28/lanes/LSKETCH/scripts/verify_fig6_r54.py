#!$T90_PY
"""Round 54 verification of the composed Main_Fig6 (lane LSKETCH), continuing Main_Fig6/verify/checks.txt after compose_r54_fig6.py:
Ghostscript txtwrite census against V30's Main_Fig6 (same strings, every string at its declared round-54 size: title 14 -> 18, notes
12 -> 15, card captions 10 -> 14, triangle labels 12 -> 14), no text below 12 pt, no two text boxes overlap, the build record's audit
clean, the placement rules (fork run at least 25 mm, arrowheads 1.5 mm clear of the '90% SpO2' labels, margins), page width = V30's,
height declared. Ghostscript only (no PyMuPDF on the composed sheet)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import json, os, sys, time
sys.dont_write_bytecode = True
T90 = paths.FIGURE_ROOT; LANE = f"{paths.FIGURE_ROOT}/ROUND54_2026-09-28/lanes/LSKETCH"; SD = f"{LANE}/Main_Fig6"; V30 = f"{paths.FIGURE_ROOT}/ROUND52_2026-09-26/figures/NEW_FINAL_SET_V30/Main_Fig6.pdf"
sys.path.insert(0, f"{LANE}/scripts"); from check_bands_r54 import census, overlaps, sha, TOL   # noqa: E402
from pypdf import PdfReader   # noqa: E402
NOTE, CAP, TRI, TITLE = 15.0, 14.0, 14.0, 18.0
EXPECT = {"Figure 6": TITLE, "90% SpO2": NOTE, "sleep duration": CAP, "sleep apnea": CAP, "141 PSG": TRI, "parameters": TRI, "disease risk": TRI, "oxygen": TRI, "all other": TRI, "PSG": TRI, "19,173 patients": TRI, "high": TRI, "low": TRI}
OLD_EXPECT = {"Figure 6": 14.0, "90% SpO2": 12.0, "sleep duration": 10.0, "sleep apnea": 10.0}   # V30's sizes for the record (the rest 12)
OUT = f"{SD}/Main_Fig6.pdf"; checks = []
P = lambda ok, msg: checks.append(("PASS " if ok else "FAIL ") + msg)
prev = []   # the compose's own checks (the first PASS/FAIL block), so that a re-run of this script neither duplicates its lines nor keeps an old header
for l in open(f"{SD}/verify/checks.txt").read().strip().splitlines():
    if l.startswith(("PASS", "FAIL")):
        if "build record:" in l: break
        prev.append(l)
rec = json.load(open(f"{LANE}/work/records_r54_Main_Fig6_tri_r54.json"))
P(rec["audit"] == [] and abs(rec["min_pt"] - 14.0) < 1e-6, f"build record: bioglyph audit clean {rec['audit']}, smallest font placed {rec['min_pt']} pt (triangle labels and card captions 14, notes 15)")
P(min(rec["fork_run_mm"]) >= 25.0, f"fork run (tip x - card right edge) {rec['fork_run_mm']} mm, rule at least 25")
g = rec["geometry"]
for key, i in (("top", 0), ("bottom", 1)):
    lb = g["traces"][key]["label_box"]; tip = g["fork"][i]["tip"][0]; P(tip <= lb[0] - 1.5 + 2e-3, f"{key} fork arrowhead at x {tip:.3f} is 1.5 mm clear of the 15-pt '90% SpO2' starting at {lb[0]:.3f}")
P(rec["left_margin_mm"] > 4.0 and rec["right_margin_mm"] > 5.0, f"content margins: triangle labels start {rec['left_margin_mm']} mm from the left edge (rule above 4), organ panels end {rec['right_margin_mm']} mm before the right edge (rule above 5)")
P(rec["title_gap_mm"] >= 2.0, f"triangle titles {rec['title_gap_mm']} mm apart (rule at least 2)")
pg = PdfReader(OUT).pages[0]; w, h = float(pg.mediabox.width), float(pg.mediabox.height); pg0 = PdfReader(V30).pages[0]; w0, h0 = float(pg0.mediabox.width), float(pg0.mediabox.height)
P(abs(w - w0) < 0.05 and abs(h - 303.464) < 0.05, f"page box {w:.3f} x {h:.3f} pt: width = V30's {w0:.3f}, height declared 303.464 (V30 {h0:.3f}: the title strip 20 pt, was 16)")
new = census(OUT, f"{SD}/verify/new_txt.xml"); old = census(V30, f"{SD}/verify/old_txt.xml")
ns, os_ = sorted(l["text"] for l in new), sorted(l["text"] for l in old)
P(ns == os_, f"census: the same {len(ns)} strings as V30's Main_Fig6 (multiset): lost {sorted(set(os_) - set(ns))} gained {sorted(set(ns) - set(os_))}")
bad = [(l["text"], round(l["size"], 3), EXPECT.get(l["text"])) for l in new if l["text"] not in EXPECT or abs(l["size"] - EXPECT[l["text"]]) > TOL]
P(not bad, f"census: every string at its declared size (title {TITLE}, notes {NOTE}, card captions {CAP}, triangle labels {TRI}, tolerance {TOL} pt): " + (f"all {len(new)} lines" if not bad else str(bad)))
badold = [(l["text"], round(l["size"], 3)) for l in old if abs(l["size"] - OLD_EXPECT.get(l["text"], 12.0)) > TOL]
P(not badold, f"census: V30's sizes as recorded (title 14, notes 12, card captions 10, triangle labels 12): {badold or 'all as recorded'}")
P(min(l["size"] for l in new) >= 12.0 - TOL, f"no text below 12 pt on the sheet: smallest {min(l['size'] for l in new):.4f} pt")
ov = overlaps(new); P(not ov, f"no two text boxes overlap on the composed sheet (txtwrite char boxes, ascent 0.718, descent 0.207, 0.5 pt): {ov or 'none'}")
ttl = [l for l in new if l["text"] == "Figure 6"][0]; P(ttl["box"][1] > 3.0 and ttl["box"][3] < 20.0 + 0.207 * 18.0 + 0.1 and ttl["box"][0] >= 35.9, f"title box {ttl['box']} (pt from the top-left): cap top {ttl['box'][1]:.2f} pt below the page edge, the descender of 'g' {ttl['box'][3] - 20.0:.2f} pt below the 20-pt strip as in round 49 (14-pt title, 16-pt strip)")
below = [l for l in new if l["text"] != "Figure 6"]; P(min(l["box"][1] for l in below) > ttl["box"][3], f"the sheet's first text below the strip starts at y {min(l['box'][1] for l in below):.2f} pt (title bottom {ttl['box'][3]:.2f})")
checks.append("INFO census lines (text, size): " + str([(l["text"], round(l["size"], 2)) for l in new]))
checks.append("INFO Figure 6 has no panel letters (one panel): only the title is stamped. The build's own texts come from the sketch engine (Chrome), the 150-dpi PNG from Ghostscript")
P(os.path.exists(f"{SD}/Main_Fig6_150dpi.png"), "Main_Fig6_150dpi.png present (Ghostscript render by the compose)")
ok = all(c.startswith("PASS") or c.startswith("INFO") for c in prev + checks)
head = [f"Round 54 lane LSKETCH, Main_Fig6, {time.strftime('%Y-%m-%d %H:%M:%S')}. NEW = {OUT} (sha256 {sha(OUT)[:16]}), OLD = V30 {V30} (sha256 {sha(V30)[:16]}). Build = work/build/Main_Fig6_tri_r54.pdf (sha256 {sha(f'{LANE}/work/build/Main_Fig6_tri_r54.pdf')[:16]}). {rec['moved']}"]
open(f"{SD}/verify/checks.txt", "w").write("\n".join(head + prev + checks) + ("\nRESULT ALL PASS\n" if ok else "\nRESULT FAIL\n"))
print("Main_Fig6 verify", "ALL PASS" if ok else "FAIL"); [print("  ", c) for c in checks if not c.startswith("PASS")]
sys.exit(0 if ok else 1)
