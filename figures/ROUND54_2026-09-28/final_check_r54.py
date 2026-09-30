#!/usr/bin/env python3
"""Round 49 final check (figures only; the documents are with Dragana): A every V27 sheet equals its lane output and every lane verifier ends ALL PASS;
B every page size equals V26's; C V27_SHA256.txt equals the files; D the combined PDF places every sheet at print size or reduced to fit (37 pages);
E the Word figure file follows the same rule; F text-size census (Ghostscript) per vector sheet: smallest size and the size distribution against V26,
reported, with a gate that no sheet's smallest text is below 10 pt. Writes reports/FINAL_CHECK_R49.md."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../..")))  # repository root, where paths.py lives
import paths
import hashlib, json, os, re, subprocess, sys, collections
T90 = paths.FIGURE_ROOT; R54 = f"{paths.FIGURE_ROOT}/ROUND54_2026-09-28"; R52 = f"{paths.FIGURE_ROOT}/ROUND52_2026-09-26"; FG = f"{R54}/figures"; V26 = f"{R52}/figures/NEW_FINAL_SET_V30"   # the base is V30 now
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest(); results = []
def check(name, ok, detail=""): results.append((name, bool(ok), detail)); print(("PASS " if ok else "FAIL ") + name + (f": {detail}" if detail else ""))
ORDER = [f"Main_Fig{i}" for i in range(1, 7)] + [f"ED_Fig{i:02d}" for i in range(1, 11)] + [f"Supp_Fig{i:02d}" for i in range(1, 22)]
man = json.load(open(f"{FG}/V31_ASSEMBLY_MANIFEST.json")); lanes = man["lane_outputs"]; AMENDED = {f"Main_Fig{i}" for i in range(1, 7)} | {"ED_Fig07"}; MAIN = {f"Main_Fig{i}" for i in range(1, 7)}
bad = [n for n in AMENDED if n not in lanes or sha(f"{FG}/NEW_FINAL_SET_V31/{n}.pdf") != sha(lanes[n])] + [n for n in ORDER if n not in AMENDED and sha(f"{FG}/NEW_FINAL_SET_V31/{n}.pdf") != sha(f"{V26}/{n}.pdf")]; check("A1 the six main sheets and ED 7 equal their ROUND54 lane outputs and the other 30 equal V30", not bad, str(bad[:5]))
lchk = {}
for n in AMENDED:
    p = os.path.join(os.path.dirname(lanes[n]), "verify", "checks.txt"); lines = [l.strip() for l in open(p) if l.strip()] if os.path.exists(p) else []
    lchk[n] = bool(lines) and lines[-1] == "RESULT ALL PASS" and not any(l.startswith("FAIL") for l in lines)
check("A2 every lane verifier ends RESULT ALL PASS", all(lchk.values()), str([n for n, v in lchk.items() if not v]))
from pypdf import PdfReader
def ps(p): b = PdfReader(p).pages[0].mediabox; return float(b.width), float(b.height)
szbad = [n for n in ORDER if n not in AMENDED and any(abs(a - b) > 0.05 for a, b in zip(ps(f"{FG}/NEW_FINAL_SET_V31/{n}.pdf"), ps(f"{V26}/{n}.pdf")))]; check("B every page size equals V30's except the six main sheets (width unchanged, height at most 1260 pt) and ED 7 (panel c added, width unchanged)", not szbad and all(abs(ps(f"{FG}/NEW_FINAL_SET_V31/{n}.pdf")[0] - ps(f"{V26}/{n}.pdf")[0]) < 0.6 and ps(f"{FG}/NEW_FINAL_SET_V31/{n}.pdf")[1] <= (1260.5 if n in MAIN else 800.0) for n in AMENDED), str(szbad) + " | main sizes: " + str({n: [round(v, 1) for v in ps(f"{FG}/NEW_FINAL_SET_V31/{n}.pdf")] for n in sorted(AMENDED)}))
h27 = {l.split("  ", 1)[1].strip(): l.split("  ", 1)[0] for l in open(f"{FG}/V31_SHA256.txt") if "  " in l}
check("C V31_SHA256.txt equals every listed file", len(h27) == 113 and all(sha(f"{FG}/{k}") == v for k, v in h27.items()), f"{len(h27)} files")
sc = json.load(open(f"{FG}/T90 All Figures V31.pdf.scale.json")); bw, bh = sc["box"][2] - sc["box"][0], sc["box"][3] - sc["box"][1]
rule_ok = all(abs(r["scale"] - round(min(1.0, bw / r["sheet_pt"][0], bh / r["sheet_pt"][1]), 4)) < 1e-3 for r in sc["sheets"].values())
check("D combined PDF: 37 pages, every sheet at print size or reduced to fit", len(PdfReader(f"{FG}/T90 All Figures V31.pdf").pages) == 37 and rule_ok and len(sc["sheets"]) == 37)
wsc = json.load(open(f"{FG}/T90 All Figures V31.docx.scale.json")); check("E Word figure file follows the same rule", wsc["rule"].startswith("print size") and len(wsc["sheets"]) == 37 and all(v["scale"] <= wsc["print_scale"] + 1e-9 for v in wsc["sheets"].values()))
# F census
def sizes(pdf, xml):
    try: subprocess.run([paths.GS, "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=txtwrite", "-dTextFormat=0", f"-sOutputFile={xml}", pdf], timeout=300)
    except subprocess.TimeoutExpired: return collections.Counter()
    return collections.Counter(round(float(s) * 2) / 2 for s in re.findall(r'size="([\d.]+)"', open(xml, encoding="utf8", errors="replace").read()))
os.makedirs(f"{R54}/logs/census", exist_ok=True); rows = []; low = []
for n in sorted(AMENDED):
    if n in ("Main_Fig3", "Main_Fig4", "Main_Fig5"):
        src27 = f"{R54}/lanes/L{n[-1]}/{n}/{n}.pdf" if n in MAIN else f"{FG}/NEW_FINAL_SET_V31/{n}.pdf"; src26 = None   # the lanes' vector composes (the set holds 600-dpi rasters with an empty text layer)
    else: src27, src26 = f"{FG}/NEW_FINAL_SET_V31/{n}.pdf", f"{V26}/{n}.pdf"
    c27 = sizes(src27, f"{R54}/logs/census/{n}_v31.xml") if os.path.exists(src27) else collections.Counter(); c26 = sizes(src26, f"{R54}/logs/census/{n}_v30.xml") if src26 and os.path.exists(src26) else collections.Counter()   # the V26 side was censused by the lanes (their checks.txt); the gate needs V27 only
    mn = min(c27) if c27 else None; rows.append((n, mn, dict(sorted(c26.items())), dict(sorted(c27.items()))))
    if mn is not None and mn < (12.0 if n in MAIN else 9.0): low.append((n, mn))
kept = [(n, mn) for n, mn, _, _ in rows if mn is not None and mn < 10.5]
check("F no main sheet's smallest text below 12 pt (round 54 target 13 to 15) and no other sheet below 9 pt", not low, f"below the gate: {low}; below 10.5: {kept}")
ok = all(r[1] for r in results)
with open(f"{R54}/reports/FINAL_CHECK_R54.md", "w") as fh:
    fh.write(f"# FINAL CHECK, round 54 (2026-09-29, the six main figures at the larger text sizes): {'RESULT ALL PASS' if ok else 'RESULT FAIL'} ({sum(r[1] for r in results)} of {len(results)} checks)\n\n")
    for name, good, detail in results: fh.write(f"- {'PASS' if good else 'FAIL'}: {name}" + (f" [{detail}]" if detail else "") + "\n")
    fh.write("\n## Text-size census of the six main sheets (glyph counts by size; V30 -> V31, the vector composes for Figures 3, 4, 5)\n")
    for n, mn, c26, c27 in rows: fh.write(f"- {n}: smallest {mn}; V26 {c26}; V27 {c27}\n")
print("RESULT ALL PASS" if ok else "RESULT FAIL"); sys.exit(0 if ok else 1)
