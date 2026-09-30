#!$T90_PY
"""Round 54 band checks (lane LSKETCH). For each of the five bands: (1) the build record's audit is empty and its smallest font is the
declared minimum, (2) the page box (pypdf mediabox) is recorded, (3) a Ghostscript txtwrite census against the round-49 band (V30's
source): the same strings (band a's second caption now on two declared lines) at the declared sizes, nothing below 12 pt, (4) no two
text boxes overlap (txtwrite char boxes rebuilt into lines), (5) the placement rules from the records (fork run at least 25 mm,
arrowheads 1.5 mm clear of the '90% SpO2' labels, ramp-end labels clear of the severity caption, band a's left-margin and bottom
rules), (6) a 150-dpi Ghostscript render to look at. Writes work/verify/<band>_checks.txt (ending RESULT ALL PASS or RESULT FAIL) and,
when all five pass, BANDS_READY.txt (written whole, then renamed into place). Ghostscript renders and censuses only."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import hashlib, html, json, os, re, subprocess, sys, time
from pypdf import PdfReader
T90 = paths.FIGURE_ROOT; LANE = f"{paths.FIGURE_ROOT}/ROUND54_2026-09-28/lanes/LSKETCH"; BUILD = f"{LANE}/work/build"; VER = f"{LANE}/work/verify"
R49B = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LSKETCH/work/build"; GS = paths.GS; TOL = 0.03
BANDS = ["band_a_r54", "band_1f_r54", "band_3d_r54", "band_4e_r54", "band_5d_r54"]
NOTE, LAB = 15.0, 14.0
EXPECT = {  # every string of the band (one census line each) with its declared round-54 size
    "band_a_r54": {**{s: LAB for s in ["breathing", "(airflow)", "heart rhythm", "(ECG)", "brain waves (EEG)", "leg movements", "(leg EMG)", "blood oxygen", "(pulse oximetry)", "sleep structure", "(EEG and staging)",
                                       "one night of sleep, 141 parameters", "Which parameter best predicts", "54 future diseases and death?"]}, "years": NOTE},
    "band_1f_r54": {"most predictive of 141 parameters": NOTE, "90% SpO2": NOTE, "one night of sleep": NOTE, "Low oxygen": NOTE},
    "band_3d_r54": {"90% SpO2": NOTE, "short, normal or long sleep": LAB},
    "band_4e_r54": {"90% SpO2": NOTE, "No apnea": NOTE, "(AHI <5)": NOTE, "Severe": NOTE, "(AHI 30 or more)": NOTE, "sleep apnea severity": LAB},
    "band_5d_r54": {"90% SpO2": NOTE, "treat sleep apnea": LAB}}
REWRAP = {"band_a_r54": {"Which parameter best predicts 54 future diseases and death?": ["Which parameter best predicts", "54 future diseases and death?"]}}   # round-49 string -> round-54 lines
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()


def census(pdf, xml):
    subprocess.run([GS, "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=txtwrite", "-dTextFormat=0", f"-sOutputFile={xml}", pdf], check=True, timeout=300)
    txt = open(xml, encoding="utf8", errors="replace").read(); chars = []
    for m in re.finditer(r'<span bbox="([^"]+)" font="([^"]*)" size="([\d.]+)">(.*?)</span>', txt, flags=re.S):   # a span holds one or more chars
        for c in re.finditer(r'<char bbox="([^"]+)" c="([^"]*)"/>', m.group(4)):
            x0, y0, x1, y1 = map(float, c.group(1).split()); chars.append((y0, x0, x1, float(m.group(3)), m.group(2), html.unescape(c.group(2))))
    chars.sort(key=lambda c: (round(c[0], 1), c[1])); lines = []
    for y, x0, x1, size, font, c in chars:
        ln = lines[-1] if lines else None
        if ln and abs(ln["y"] - y) < 0.6 and abs(ln["size"] - size) < 0.05 and x0 - ln["x1"] <= 0.9 * size:
            ln["text"] += c; ln["x1"] = max(ln["x1"], x1)
        else:
            lines.append({"text": c, "y": y, "x0": x0, "x1": x1, "size": size, "font": font})
    for ln in lines: ln["text"] = ln["text"].replace("\xa0", " ").strip(); ln["box"] = [round(ln["x0"], 2), round(ln["y"] - 0.718 * ln["size"], 2), round(ln["x1"], 2), round(ln["y"] + 0.207 * ln["size"], 2)]
    return [ln for ln in lines if ln["text"]]


def overlaps(lines, min_pt=0.5):
    out = []
    for i, a in enumerate(lines):
        for b in lines[i + 1:]:
            ox = min(a["box"][2], b["box"][2]) - max(a["box"][0], b["box"][0]); oy = min(a["box"][3], b["box"][3]) - max(a["box"][1], b["box"][1])
            if ox > min_pt and oy > min_pt: out.append((a["text"], b["text"], round(min(ox, oy), 2)))
    return out


def check_band(stem):
    pdf = f"{BUILD}/{stem}.pdf"; old = f"{R49B}/{stem.replace('_r54', '_r49')}.pdf"; rec = json.load(open(f"{LANE}/work/records_r54_{stem}.json")); checks = []
    P = lambda ok, msg: checks.append(("PASS " if ok else "FAIL ") + msg)
    P(rec["audit"] == [], f"build audit (bioglyph collisions, bounds, type floor) clean: {rec['audit']}")
    exp = EXPECT[stem]; smallest = min(exp.values()); P(abs(rec["min_pt"] - smallest) < 1e-6, f"smallest font placed {rec['min_pt']} pt = the declared minimum {smallest}")
    pg = PdfReader(pdf).pages[0]; w, h = float(pg.mediabox.width), float(pg.mediabox.height); ow, oh = (float(v) for v in (PdfReader(old).pages[0].mediabox.width, PdfReader(old).pages[0].mediabox.height))
    P(abs(w - ow) < 0.05, f"page box {w:.3f} x {h:.3f} pt (round 49: {ow:.3f} x {oh:.3f}, height {h - oh:+.3f} pt, {(h - oh) * 25.4 / 72:+.2f} mm), width unchanged")
    new = census(pdf, f"{VER}/{stem}_txt.xml"); prev = census(old, f"{VER}/{stem.replace('_r54', '_r49')}_txt.xml")
    old_strings = []
    for ln in prev: old_strings += REWRAP.get(stem, {}).get(ln["text"], [ln["text"]])
    new_strings = [ln["text"] for ln in new]
    P(sorted(old_strings) == sorted(new_strings), f"census: the round-49 band's strings (declared re-wrap applied) = the new band's strings, {len(new_strings)} lines: lost {sorted(set(old_strings) - set(new_strings))} gained {sorted(set(new_strings) - set(old_strings))}")
    bad = [(ln["text"], ln["size"], exp.get(ln["text"])) for ln in new if ln["text"] not in exp or abs(ln["size"] - exp[ln["text"]]) > TOL]
    P(not bad, f"census: every string at its declared size (notes {NOTE}, labels and captions {LAB}, tolerance {TOL} pt): " + ("all " + str(len(new)) + " lines" if not bad else str(bad)))
    P(min(ln["size"] for ln in new) >= 12.0 - TOL, f"no text below 12 pt: smallest {min(ln['size'] for ln in new):.4f}")
    ov = overlaps(new); P(not ov, f"no two text boxes overlap (txtwrite char boxes, ascent 0.718, descent 0.207, 0.5 pt): {ov or 'none'}")
    ex = {ln["text"]: (round(ln["size"], 3), ) for ln in new}; sizes_old = sorted({round(ln["size"], 2) for ln in prev}); sizes_new = sorted({round(ln["size"], 2) for ln in new})
    checks.append(f"INFO census sizes: round 49 {sizes_old} -> round 54 {sizes_new}, lines {[(ln['text'], round(ln['size'], 2)) for ln in new]}")
    if stem in ("band_3d_r54", "band_4e_r54", "band_5d_r54"):
        fk = rec["fork"]; run = fk["tips_x"] - fk["src"][0]; P(run >= 25.0, f"fork run (tips x - source x) {run:.2f} mm, rule at least 25")
        for key in ("top", "bottom"):
            lb = rec["traces"][key]["label_box"]; P(fk["tips_x"] <= lb[0] - 1.5 + 2e-3, f"{key} fork arrowhead at x {fk['tips_x']:.3f} is 1.5 mm clear of the '90% SpO2' label starting at {lb[0]:.3f} (clearance {lb[0] - fk['tips_x']:.3f})")
        oa = rec["organ_arrows"]; P(min(a["length_mm"] for a in oa) >= 20.0, f"organ arrows {[a['length_mm'] for a in oa]} mm long (head 2.72)")
    if stem == "band_4e_r54":
        ramp = [ln for ln in new if ln["text"] in ("(AHI <5)", "(AHI 30 or more)")]; cap = [ln for ln in new if ln["text"] == "sleep apnea severity"][0]
        gap = cap["box"][1] - max(ln["box"][3] for ln in ramp); P(gap > 0.5, f"severity caption box starts {gap:.2f} pt below the ramp-end labels' boxes (rule: clear)")
        P(rec["h_mm"] - rec["sev_label_box"][3] >= 2.3, f"caption bottom {rec['sev_label_box'][3]:.2f} mm, band {rec['h_mm']} mm: bottom margin {rec['h_mm'] - rec['sev_label_box'][3]:.2f} (rule 2.3)")
    if stem == "band_a_r54":
        P(all(t["label_box"][0] >= 16.68 * 0.352778 - 0.01 for t in rec["tiles"]), f"no tile label left of the scene margin 5.88 mm (leftmost {min(t['label_box'][0] for t in rec['tiles']):.2f}), scene shift {rec['scene_shift_mm']} mm")
        P(rec["bottom_margin_mm"] >= 2.3, f"lowest tile label to the band bottom {rec['bottom_margin_mm']} mm (rule 2.3)")
        c2 = rec["captions"]["2"]; P(c2["label_box"][3] < rec["organ_block"][1] - 1.0 and c2["label_box"][2] <= 341.72 - 0.5, f"second caption on two lines, box {c2['label_box']} clear of the organ block top {rec['organ_block'][1]} and of the right edge")
    if stem == "band_1f_r54":
        P(rec["trace"]["label_box"][0] >= 1.0, f"'90% SpO2' starts {rec['trace']['label_box'][0]:.2f} mm inside the page (trace at x {rec['tx']})")
        P(rec["label"]["box"][2] < rec["organs"]["ink_x"][0] - 5.0 and rec["label"]["box"][0] > rec["arrow"]["p0"][0], f"'Low oxygen' box {rec['label']['box']} over the arrow, clear of the organs")
    png = f"{VER}/{stem}_150.png"; r = subprocess.run([GS, "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=png16m", "-r150", f"-sOutputFile={png}", pdf], capture_output=True, timeout=600)
    P(r.returncode == 0 and os.path.exists(png), f"Ghostscript 150 dpi render {png}")
    ok = all(c.startswith("PASS") or c.startswith("INFO") for c in checks)
    head = [f"Round 54 lane LSKETCH, {stem}, {time.strftime('%Y-%m-%d %H:%M:%S')}. NEW = {pdf} (sha256 {sha(pdf)[:16]}), OLD = {old} (sha256 {sha(old)[:16]}, the round-49 band behind V30). {rec['moved']}"]
    open(f"{VER}/{stem}_checks.txt", "w").write("\n".join(head + checks) + ("\nRESULT ALL PASS\n" if ok else "\nRESULT FAIL\n"))
    print(stem, "ALL PASS" if ok else "FAIL"); [print("   ", c) for c in checks if not c.startswith("PASS")]
    return ok, (w, h, oh)


if __name__ == "__main__":
    os.makedirs(VER, exist_ok=True); res = {}
    for stem in BANDS: res[stem] = check_band(stem)
    if all(v[0] for v in res.values()):
        lines = [f"{stem}.pdf\t{res[stem][1][0]:.3f} x {res[stem][1][1]:.3f} pt\tsha256={sha(f'{BUILD}/{stem}.pdf')}" for stem in BANDS]
        notes = []   # exactly one line per band (name, page box, sha256), no comments: the composing lanes parse it; the notes are in REPORT.md
        tmp = f"{LANE}/BANDS_READY.txt.part"; open(tmp, "w").write("\n".join(lines + notes) + "\n"); os.replace(tmp, f"{LANE}/BANDS_READY.txt"); print("BANDS_READY.txt written")
    else:
        print("NOT READY:", [k for k, v in res.items() if not v[0]]); sys.exit(1)
