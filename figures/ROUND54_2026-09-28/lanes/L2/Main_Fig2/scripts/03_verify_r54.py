#!$T90_PY
"""Round 54 (2026-09-29), lane L2, Main_Fig2 step 3: proof of Main_Fig2.pdf (the round-54 sizes, the L-shaped layout) against the V30 sheet
and the numbers files. Adapted from the round-49 verifier (byte copy beside as 03_verify_r49_R49_ORIGINAL.py). Ghostscript only on the
composed sheets (page box, text layer, renders, the pdfwrite flatten); PyMuPDF only on the flat panel pages of this lane and of round 49
(get_drawings, rawdict text) for the vector read-back and the data-mark comparison. Writes verify/checks.txt (ends RESULT ALL PASS or
RESULT FAIL), verify/census.json, verify/values_agreement.txt, verify/Main_Fig2_drawn.json, verify/Main_Fig2_printed_values.csv,
verify/CHANGES_Main_Fig2.csv, verify/verify_record.json, verify/crops/*.png, Main_Fig2_150dpi.png (gs, 150 dpi) and work/Main_Fig2_flat.pdf.
Proves: (1) the page box is the declared one (968.66 x 1107.04, at most 1110) and gs PDFINFO agrees; (2) text census against V30: the same
multiset of strings after the five declared header joins, every string at its class size (14 ticks, labels, keys, values; 15 axis titles,
column headers, panel titles; 13 notes and the asterisk key; 18 letters and title), no text below 12 pt, the same font set; (3) letters
and title at 18 pt where the composer put them; (4) provenance; (5) the printed strings on the sheet = the verifier's own strings from the
files; (6) every marker, interval, ribbon and star read back from the flat panels' vectors agrees with the files; (7) marks preserved:
the same drawing classes and counts, the same marker sizes as the round-49 (V30) panels; (8) no glyph box overlaps another, no text
touches a data mark; (9) every panel's content inside its region; (10) every string accounted for by a drawn record; (11) the flat pdf
carries the same text layer; renders and crops."""
import csv, gc, json, math, os, re, sys, time
from collections import Counter
import numpy as np, pandas as pd
import fitz
from PIL import Image
from fig2_common import *
import v14lib
from v14lib import printed_values_csv
import lib_r54 as R

NEW = os.path.join(LANE, "Main_Fig2.pdf"); V30PDF = V30
LAY = layout(); PH = float(LAY["page_h"]); SPLIT = float(LAY["split_x"]); BC_Y = float(LAY["bc_y"])
DR = {p: json.load(open(hydrated(os.path.join(WORK, f"panel_{p}_drawn.json")))) for p in "abc"}
DR49 = {p: json.load(open(hydrated(os.path.join(WORK, f"r49_panel_{p}_drawn.json")))) for p in "abc"}
CLOG = json.load(open(os.path.join(WORK, "compose_log.json")))
T14, T15, T13, T18 = TARGET["tick"], TARGET["head"], TARGET["note"], TARGET["letter"]
ck = R.Checks(f"Round 54 lane L2, Main_Fig2, {R.now()}. INPUT = V30 {V30PDF} (sha256 {R.sha256(V30PDF)[:16]}). NEW = {NEW} (sha256 {R.sha256(NEW)[:16]}). "
              f"Change: every text element at the round-54 class size (14 pt ticks, row labels, keys, HR and q values; 15 pt axis titles, column headers, panel titles; 13 pt 'no estimate' and the asterisk key; "
              f"18 pt letters and title), the sheet re-laid as an L-shape (panel a one column on the left with the organ headers as bold rows above their blocks and the five two-line V13 headers on one line, "
              f"panel b above panel c on the right, panel c as 3 x 4 boxes), page 968.66 x {PH} (V30 1184.76). No printed number, colour, marker size, line width or data value changes.")
def check(ok, msg, detail=""): return ck.log(msg, bool(ok), detail)
def info(msg): ck.info(msg)
def w_(b): return b[2] - b[0]
def h_(b): return b[3] - b[1]
def cen(b): return ((b[0] + b[2]) / 2, (b[1] + b[3]) / 2)
rec = dict(round=54, sheet="Main_Fig2", new=NEW, new_sha256=R.sha256(NEW), v30=V30PDF, v30_sha256=R.sha256(V30PDF), written=R.now())

# ------------------------------------------------------------------ 1 page box
n1, w1, h1 = R.pypdf_mediabox(NEW); n0, w0, h0 = R.pypdf_mediabox(V30PDF); gn, gw, gh = R.gs_pdfinfo(NEW)
check(n1 == gn == 1 and abs(w1 - PAGE_W) < 0.05 and abs(h1 - PH) < 0.05 and abs(gw - w1) < 0.05 and abs(gh - h1) < 0.05 and h1 <= 1110.0 and abs(w1 - w0) < 0.05,
      f"page box is the declared one (pypdf mediabox): new {w1} x {h1} (layout.json {PAGE_W} x {PH}, at most 1110 pt, width = V30's {w0}); V30 {w0} x {h0}; gs PDFINFO {gw} x {gh}, one page")
rec["page"] = dict(new=[w1, h1], v30=[w0, h0], gs=[gw, gh])

# ------------------------------------------------------------------ 2 text census (gs txtwrite -dTextFormat=0), class rule against V30
S1 = R.spans_xml(NEW, os.path.join(WORK, "new_txtwrite.xml"), "txt_new"); S0 = R.spans_xml(V30PDF, os.path.join(WORK, "v30_txtwrite.xml"), "txt_v30")
STR1 = R.strings_r54(S1, SPLIT, BC_Y)                                        # the new sheet's strings by region a/b/c
STR0 = R.strings_r54(S0, 479.0, 780.76, region_fn=lambda x, y: "c" if y >= 780.76 else ("a" if x < 479.0 else "b"))   # V30 regions: c the full width below 780.76, a/b split at 479
V30_CAPS_C = {"Hazard ratio vs T90 0-1% (95% CI)", "Time with oxygen saturation below 90%, % of the recording"}
V30_HEADS_A = {"HR (95% CI)", "q"}
def rule(s):
    """The size class of a V30 string on the round-54 sheet."""
    t, sz, reg, bold = s["text"], s["size"], s["region"], "Bold" in s["font"]
    if t == "Figure 2" or (len(t) == 1 and t in "abc" and sz >= 13.5): return T18
    if reg == "a":
        if abs(sz - 11.0) < 0.01: return T14                                   # row labels, headers, key, ticks
        if abs(sz - 10.5) < 0.01: return T15 if t in V30_HEADS_A else T14      # column headers 15, values 14
        if abs(sz - 12.0) < 0.01: return T15                                   # caption
    if reg == "b":
        if abs(sz - 11.0) < 0.01: return T14
        if abs(sz - 12.0) < 0.01: return T15
        if abs(sz - 10.0) < 0.01 and t == "no estimate": return T13
    if reg == "c":
        if abs(sz - 14.2) < 0.01: return T15                                   # the two captions
        if abs(sz - 12.4) < 0.01:
            if t == "***": return T13 if s["y"] > 1100 else T14               # the key star 13, the panel stars 14
            return T15 if t[0].isalpha() and t not in V30_CAPS_C else T14      # panel titles 15, tick and band labels 14
        if abs(sz - 11.8) < 0.01: return T13 if t == "q < 0.001" else T14      # the asterisk key 13, the band key 14
    return None
JOINS = [(["Metabolic,", "kidney and liver"], "Metabolic, kidney and liver"), (["Infection,", "blood and other"], "Infection, blood and other"),
         (["Neurological", "and psychiatric"], "Neurological and psychiatric"), (["Musculoskeletal", "and sensory"], "Musculoskeletal and sensory"),
         (["Negative", "controls"], "Negative controls")]
CEN = R.census_r54(STR1, STR0, rule, JOINS)
json.dump(CEN, open(os.path.join(VER, "census.json"), "w"), indent=1)
check(CEN["same_strings"], f"census strings: the same multiset of strings as V30 after the five declared header joins ({CEN['n_new']} strings new, {CEN['n_old']} V30, {CEN['n_old_after_joins']} after the joins)",
      f"lost {CEN['lost'][:8]} gained {CEN['gained'][:8]}")
info("declared joins (V30 two-line headers set on one line, same words): " + "; ".join(f"{' + '.join(p)} -> '{n}'" for p, n in CEN["joins"]) + ". Declared duplicates: none (one axis per panel, every string once)")
check(CEN["sizes_ok"], f"census sizes: every V30 string reappears at its class size (14 ticks, row labels, keys, values; 15 axis titles, column headers, panel titles; 13 'no estimate' and the asterisk key; 18 letters and title); V30 sizes {CEN['sizes_old']} -> new {CEN['sizes_new']}",
      f"missing {CEN['missing'][:6]} unexpected {CEN['unexpected'][:6]} no-rule {CEN['no_rule'][:4]}")
info("census per class (region:V30 size->new size: strings): " + "; ".join(f"{k}: {v}" for k, v in CEN["per_class"].items()))
check(CEN["min_size_new"] >= 12.0 and CEN["min_size_new"] == T13, f"no text below 12 pt on the new sheet: smallest {CEN['min_size_new']} pt ('no estimate' and the asterisk key, 13 pt); no exceptions")
check(set(CEN["fonts_new"]) == set(CEN["fonts_old"]), f"census fonts: the same font set as V30 {sorted(CEN['fonts_new'])}", f"V30 {CEN['fonts_old']} new {CEN['fonts_new']}")
def letters(sp): return sorted([(s["text"], s["x0"], s["y"], s["font"], s["size"]) for s in sp if len(s["text"]) == 1 and s["text"] in "abc" and s["size"] >= 13.5 and "Bold" in s["font"]], key=lambda t: t[0])
L1 = letters(S1); LX = CLOG["letters"]
check(len(L1) == 3 and [l[0] for l in L1] == ["a", "b", "c"] and all(abs(l[4] - T18) < 0.01 and abs(l[1] - LX[l[0]][0]) < 1 and abs(l[2] - LX[l[0]][1]) < 1 for l in L1),
      f"letters a, b, c: Helvetica-Bold 18 pt at the composer's origins {[(l[0], l[1], l[2]) for l in L1]} (a top-left, b top-right, c right below b)")
T1 = [s for s in S1 if s["text"] == "Figure 2"]
check(len(T1) == 1 and abs(T1[0]["size"] - T18) < 0.01 and abs(T1[0]["x0"] - 14.2) < 1 and abs(T1[0]["y"] - 18.0) < 1 and "Bold" in T1[0]["font"],
      f"title 'Figure 2' once: {T1 and (T1[0]['font'], T1[0]['size'], T1[0]['x0'], T1[0]['y'])} (18 pt on baseline 18, cap top about 5 pt inside the page)")
rec["letters"] = L1; rec["title"] = T1[0] if T1 else None

# ------------------------------------------------------------------ 3 provenance and the verifier's own strings
FILES = ("bdsp_diseases_v3.csv", "negcontrols_v4_panel.csv", "cohorts.json", "lag_ladder.json", "results_v2.json", "fig2c_selection_v8.json")
PROV = {p: sidecar(f"{NUM}/{p}") for p in FILES}
for p, d in DR.items():
    for k, v in d.get("provenance", {}).items():
        check(v["sha256"] == PROV[k]["sha256"], f"panel {p}: {k} read now has the sha256 the builder recorded ({PROV[k]['sha256'][:12]}..., step {PROV[k]['step']} {PROV[k]['name']}, out {PROV[k]['output_mtime']})")
check(all(PROV[k]["v8_inputs"] for k in FILES), "every numbers file behind the sheet carries a v8 sidecar naming data_frozen_v8_2026-09 (v14lib.sidecar_v8)")
BD = pd.read_csv(f"{NUM}/bdsp_diseases_v3.csv"); NCP = pd.read_csv(f"{NUM}/negcontrols_v4_panel.csv")
LL = json.load(open(f"{NUM}/lag_ladder.json")); RS = json.load(open(f"{NUM}/results_v2.json"))["graded"]; SEL = json.load(open(f"{NUM}/fig2c_selection_v8.json"))
MATRIX_CK = f"{v14lib.SV}/dC_side_analyses/_work/hr_matrix_ck/spo2_pct_below_90.json"; MATRIX_CSV = f"{v14lib.SV}/dC_side_analyses/hr_matrix_141x52.csv"
CK = {r["outcome"]: r for r in json.load(open(hydrated(MATRIX_CK)))}; SC168 = sidecar(MATRIX_CSV)
NAME2KEY = {r.disease: r.key for r in BD.itertuples()}; BDI = BD.set_index("disease")
ck_drift = max(abs(float(CK[k][a]) - float(BDI.loc[c, b])) for c, k in NAME2KEY.items() if k in CK for a, b in (("HR", "t90_hr"), ("lo", "t90_lo"), ("hi", "t90_hi")))
check(len(CK) == 52 and all(r["measurement"] == "spo2_pct_below_90" for r in CK.values()) and ck_drift < 1e-5 and sum(1 for k in NAME2KEY.values() if k in CK) == 52,
      f"full-precision checkpoint (step 168, sha256 {sha256(MATRIX_CK)[:12]}...): 52 T90 rows, equal to bdsp_diseases_v3.csv at full precision on all 52 x 3 values (max |diff| {ck_drift:.2e}; content guard, round 38 rule); matrix csv step {SC168['step']} out {SC168['output_mtime']}")
def panel_a_strings(bd, ncp, ckd):
    outc = bd[~bd.negative_control].copy(); outc["q"] = bh(outc.t90_p.values); ctrl = bd[bd.negative_control].copy(); ctrl["q"] = bh(ctrl.t90_p.values)
    rows = {}; ncpi = ncp.set_index("condition")
    for r in outc.itertuples():
        v = (float(ckd[r.key]["HR"]), float(ckd[r.key]["lo"]), float(ckd[r.key]["hi"])) if r.key in ckd else (float(r.t90_hr), float(r.t90_lo), float(r.t90_hi))
        rows[r.disease] = dict(hr=hr_ci_text(*v), q=q_text(r.q), qv=float(r.q), hrv=v, ctrl=False)
    for r in ctrl.itertuples():
        n = ncpi.loc[r.disease]; v = (float(ckd[r.key]["HR"]), float(ckd[r.key]["lo"]), float(ckd[r.key]["hi"])) if r.key in ckd else (float(n.hr), float(n.lo), float(n.hi))
        rows[r.disease] = dict(hr=hr_ci_text(*v), q=q_text(r.q), qv=float(r.q), hrv=v, ctrl=True)
    return rows, halfup(ncp.hr.max(), 2), int((outc.q < 0.05).sum()), len(outc)
A_NEW, FLOOR_NEW, NSIG_NEW, NOUT_NEW = panel_a_strings(BD, NCP, CK)
check(len(A_NEW) == 58 and NOUT_NEW == 53 and len(BD[BD.negative_control]) == 5, f"primary file: {len(A_NEW)} rows = {NOUT_NEW} outcomes + 5 controls; FDR-significant {NSIG_NEW} of {NOUT_NEW}")
HRQ = re.compile(r"^\d\.\d\d \(\d\.\d\d-\d\.\d\d\)$")
sheet_hr = Counter(s["text"] for s in STR1 if s["region"] == "a" and s["size"] == T14 and HRQ.match(s["text"]))
sheet_q = Counter(s["text"] for s in STR1 if s["region"] == "a" and s["size"] == T14 and re.fullmatch(r"<0\.001|0\.\d{3,4}", s["text"]))
exp_hr = Counter(v["hr"] for v in A_NEW.values()); exp_q = Counter(v["q"] for v in A_NEW.values())
check(sheet_hr == exp_hr, f"panel a HR (95% CI) strings on the sheet = the verifier's own strings from the checkpoint and the csv (58; 52 ranked outcomes half up on the full-precision checkpoint, 6 rows half up on the csv)", f"sheet-not-expected {sorted((sheet_hr - exp_hr).elements())[:5]}, expected-not-sheet {sorted((exp_hr - sheet_hr).elements())[:5]}")
check(sheet_q == exp_q, f"panel a q strings on the sheet = BH across the 53 outcomes and across the 5 controls (58)", f"sheet-not-expected {sorted((sheet_q - exp_q).elements())[:5]}, expected-not-sheet {sorted((exp_q - sheet_q).elements())[:5]}")
fl = [s["text"] for s in STR1 if s["text"].startswith("Level the negative controls reach")]
check(len(fl) == 1 and fl[0] == f"Level the negative controls reach, {FLOOR_NEW}", f"floor key '{fl and fl[0]}' = max control hazard ratio half up ({FLOOR_NEW})")
sheet_hdr = [s["text"] for s in STR1 if s["region"] == "a" and "Bold" in s["font"] and s["size"] == T14 and s["text"][0].isalpha()]
check(Counter(sheet_hdr) == Counter(DR["a"]["r54"]["headers_one_line"].values()) and len(sheet_hdr) == 7, f"panel a: the 7 organ headers as bold rows: {sheet_hdr}")
PC = LL["per_condition"]; COH_N = int(json.load(open(f"{NUM}/cohorts.json"))["bdsp"]["n_analysis"])
panel = sorted([k for k, v in PC.items() if not v["negative_control"] and not v["circular"] and v["ci_lag0"][0] > 1], key=lambda k: -PC[k]["hr_lag0"])[:12]
ctrlk = sorted([k for k, v in PC.items() if v["negative_control"]], key=lambda k: -PC[k]["hr_lag0"])
SHORT = {"Ventricular arrhythmia or cardiac arrest": "Ventricular arrhythmia/arrest"}
exp_b = [f"{SHORT.get(PC[k]['condition'], PC[k]['condition'])} ({int(PC[k]['events_by_lag']['5']):,})" for k in panel + ctrlk]
sheet_b = [s["text"] for s in STR1 if s["region"] == "b" and s["size"] == T14 and re.search(r" \(\d[\d,]*\)$", s["text"])]
check(Counter(sheet_b) == Counter(exp_b) and DR["b"]["row_labels"] == exp_b and LL["cohort_n"] == COH_N, f"panel b 17 row labels on the sheet = the verifier's own top-12 (lag-0 lower bound above 1, by lag-0 hazard ratio) + 5 controls from lag_ladder.json; cohort n {COH_N}")
ne = sum(1 for s in STR1 if s["text"] == "no estimate")
check(ne == 3 == len(DR["b"]["noest"]) == sum(1 for k in panel + ctrlk for Lg in ("0", "1", "2", "5") if PC[k]["hr_by_lag"][Lg] is None), f"panel b 'no estimate' strings: {ne} (rows {DR['b']['noest']})")
G = RS["outcomes"]; sheet_titles = [s["text"] for s in STR1 if s["region"] == "c" and s["size"] == T15 and s["text"][0].isalpha() and not any(s["text"] in c for c in V30_CAPS_C)]   # gs splits the rotated y label into two pieces
check(Counter(sheet_titles) == Counter(SEL["selected_labels"]) and DR["c"]["order"] == SEL["selected_labels"] and sum(RS["band_n"].values()) == COH_N, f"panel c 10 titles on the sheet = fig2c_selection_v8.json selected_labels {SEL['selected_labels']}; band n sums to the cohort")
stars14 = sum(1 for s in STR1 if s["region"] == "c" and s["size"] == T14 and s["text"] == "***"); stars13 = sum(1 for s in STR1 if s["region"] == "c" and s["size"] == T13 and s["text"] == "***")
conds_nc = [k for k, v in G.items() if k not in set(NCP.condition) and not v.get("negative_control")]; q_top = dict(zip(conds_nc, bh([G[k][">10%"]["p"] for k in conds_nc])))
check(all(q_top[c] < 0.001 for c in SEL["selected_labels"]) and stars14 == 10 and stars13 == 1, f"panel c: every drawn panel's top-band BH q below 0.001 across {len(conds_nc)} conditions (largest {max(q_top[c] for c in SEL['selected_labels']):.2e}); {stars14} '***' at 14 pt on the panels + 1 at 13 pt in the key")

# ------------------------------------------------------------------ 4 flat panels: vectors and text (PyMuPDF on the flat pages only)
def spans_of(page):
    out = []
    for b in page.get_text("rawdict")["blocks"]:
        if b["type"] != 0: continue
        for l in b["lines"]:
            for s in l["spans"]:
                txt = "".join(ch["c"] for ch in s["chars"]); org = s["chars"][0]["origin"]
                out.append(dict(text=txt, font=s["font"], size=round(s["size"], 3), bbox=[round(v, 3) for v in s["bbox"]], origin=[round(org[0], 3), round(org[1], 3)], color=s.get("color"), dir=list(l.get("dir", (1, 0)))))
    return out
def drawings_of(page):
    out = []
    for d in page.get_drawings():
        fill = d.get("fill"); stroke = d.get("color")
        op = "B" if (fill is not None and stroke is not None) else ("f" if fill is not None else ("S" if stroke is not None else "n"))
        if op == "n": continue
        dash = d.get("dashes"); dash = None if (not dash or dash.strip() in ("[] 0", "[]", "")) else dash
        r = d["rect"]; npts = sum(1 + (3 if it[0] == "c" else 1) for it in d["items"] if it[0] in ("l", "c")) + sum(4 for it in d["items"] if it[0] == "re")
        out.append(dict(op=op, fill=hexcol(fill), stroke=hexcol(stroke), lw=round(d.get("width") or 0.0, 4), dash=bool(dash), bbox=[round(r.x0, 3), round(r.y0, 3), round(r.x1, 3), round(r.y1, 3)], npts=npts, alpha=round(d.get("fill_opacity") if d.get("fill_opacity") is not None else 1.0, 3)))
    return out
FLAT = {}; FLAT49 = {}; TXT = {}
for p in "abc":
    d = fitz.open(os.path.join(WORK, f"panel_{p}.pdf")); FLAT[p] = drawings_of(d[0]); TXT[p] = spans_of(d[0]); d.close()
    d = fitz.open(os.path.join(WORK, f"r49_panel_{p}.pdf")); FLAT49[p] = drawings_of(d[0]); d.close(); gc.collect()
# 4a marks preserved: the same drawing classes and counts as the round-49 (V30) panels, the same marker sizes
def klass(r): return (r["op"], r["fill"], r["stroke"], r["dash"], round(r["lw"], 2), r["npts"] if r["op"] == "B" else 0, r["alpha"])
def is_marker(r): return r["op"] == "B" and 0 < w_(r["bbox"]) < 12 and h_(r["bbox"]) < 12
class_mismatch = []; size_mismatch = []; NCLS = {}
for p in "abc":
    cn = Counter(klass(r) for r in FLAT[p]); co = Counter(klass(r) for r in FLAT49[p])
    if cn != co: class_mismatch.append((p, dict((k, v) for k, v in (cn - co).items()), dict((k, v) for k, v in (co - cn).items())))
    mn = Counter((r["fill"], round(w_(r["bbox"]), 1), round(h_(r["bbox"]), 1)) for r in FLAT[p] if is_marker(r)); mo = Counter((r["fill"], round(w_(r["bbox"]), 1), round(h_(r["bbox"]), 1)) for r in FLAT49[p] if is_marker(r))
    if mn != mo: size_mismatch.append((p, dict((k, v) for k, v in (mn - mo).items()), dict((k, v) for k, v in (mo - mn).items())))
    NCLS[p] = dict(n=len(FLAT[p]), n49=len(FLAT49[p]), n_markers=sum(mn.values()))
check(not class_mismatch, f"marks preserved: every drawing class (operation, fill, stroke, dash, line width, points, alpha) has the same count as on the round-49 (V30) panels (a {NCLS['a']['n']}, b {NCLS['b']['n']}, c {NCLS['c']['n']} drawings)", f"{class_mismatch[:3]}")
check(not size_mismatch, f"marks preserved: the same multiset of marker sizes and colours as the V30 panels (a {NCLS['a']['n_markers']}, b {NCLS['b']['n_markers']}, c {NCLS['c']['n_markers']} markers); positions follow the re-layout", f"{size_mismatch[:3]}")
rec["drawings"] = NCLS
# 4b values read back from the vectors against the numbers files
A_BOX = DR["a"]["layout"]["ax_box"]; ax = DR["a"]["axis"]
def inbox(r, box, pad=1.0): cx, cy = cen(r["bbox"]); return box[0] - pad < cx < box[2] + pad and box[1] - pad < cy < box[3] + pad
ta_n = sorted(r["bbox"][0] for r in FLAT["a"] if r["op"] in ("S", "B") and w_(r["bbox"]) < 0.05 and abs(h_(r["bbox"]) - 3.0) < 0.1 and abs(r["bbox"][1] - A_BOX[3]) < 0.5)
check(len(ta_n) == 7 and max(abs(x - y) for x, y in zip(ta_n, ax["tick_x"])) < 0.05, f"panel a tick marks: {len(ta_n)} at the builder's positions for ticks {ax['ticks']} (x range {ax['xlo']:.3f}-{ax['xhi']:.3f}, axis {ax['w_ax']} pt wide)")
def log_axis(ticks_x, vals):
    a_, b_, res = fit_log_axis(ticks_x, vals); assert res < 0.05, res
    return lambda x: math.exp((x - a_) / b_)
fa = log_axis(ta_n, ax["ticks"])
mk_a = sorted([r for r in FLAT["a"] if is_marker(r) and inbox(r, A_BOX)], key=lambda r: r["bbox"][1])
ci_a = [r for r in FLAT["a"] if r["op"] == "S" and h_(r["bbox"]) < 0.05 and abs(r["lw"] - 1.7) < 0.01 and A_BOX[1] + 1 < r["bbox"][1] < A_BOX[3] - 1 and r["bbox"][0] > A_BOX[0] - 1]
check(len(mk_a) == 58 and len(ci_a) == 58, f"panel a: {len(mk_a)} markers and {len(ci_a)} confidence-interval lines inside the axes box")
agree = [f"Main_Fig2 round 54: values read back from the flat panels' vectors against the numbers files (deviation as percent of the axis span)", ""]
worst_a = 0.0; span_a = math.log(ax["xhi"]) - math.log(ax["xlo"])
for m in mk_a:
    cx, cy = cen(m["bbox"]); row = min(DR["a"]["marks"], key=lambda mm: abs(mm["y_page"] - cy)); assert abs(row["y_page"] - cy) < 0.3, (row["condition"], cy)
    v = A_NEW[row["condition"]]["hrv"]; hr_read = fa(cx); ci = min(ci_a, key=lambda r: abs(cen(r["bbox"])[1] - cy)); lo_read, hi_read = fa(ci["bbox"][0]), fa(ci["bbox"][2])
    dev = max(abs(math.log(hr_read) - math.log(v[0])), abs(math.log(lo_read) - math.log(v[1])), abs(math.log(hi_read) - math.log(v[2]))) / span_a * 100
    worst_a = max(worst_a, dev); agree.append(f"a  {row['condition']:42s} file {v[0]:.3f} ({v[1]:.3f}-{v[2]:.3f})  read {hr_read:.3f} ({lo_read:.3f}-{hi_read:.3f})  dev {dev:.3f}%")
check(worst_a <= 0.5, f"panel a: 58 markers and CI ends read back through the tick marks agree with the checkpoint / bdsp_diseases_v3.csv / negcontrols_v4_panel.csv: max deviation {worst_a:.3f}% of the axis span")
fl_n = [r for r in FLAT["a"] if r["dash"] and h_(r["bbox"]) > 100 and r["stroke"] == GREY]; rf_n = [r for r in FLAT["a"] if r["dash"] and h_(r["bbox"]) > 100 and r["stroke"] == INK]
check(len(fl_n) == 1 and len(rf_n) == 1 and abs(fa(rf_n[0]["bbox"][0]) - 1.0) < 0.002 and abs(fa(fl_n[0]["bbox"][0]) - float(NCP.hr.max())) < 0.002, f"panel a dashed reference at 1 (read {rf_n and fa(rf_n[0]['bbox'][0]):.4f}) and dotted control level at {float(NCP.hr.max()):.4f} (read {fl_n and fa(fl_n[0]['bbox'][0]):.4f}), both the full height of the box")
B_BOX = DR["b"]["axes_box_page"]
tb_n = sorted(r["bbox"][0] for r in FLAT["b"] if r["op"] in ("S", "B") and w_(r["bbox"]) < 0.05 and abs(h_(r["bbox"]) - 3.0) < 0.1 and abs(r["bbox"][1] - B_BOX[3]) < 0.5)
check(len(tb_n) == 6 and max(abs(x - y) for x, y in zip(tb_n, DR["b"]["tick_x"])) < 0.05, f"panel b 6 tick marks (0.7 to 4) at {[round(x, 1) for x in tb_n]} (x range {DR['b']['xlim']}, axis {B_BOX[2] - B_BOX[0]:.1f} pt wide)")
fb_ = log_axis(tb_n, [0.7, 1.0, 1.5, 2.0, 3.0, 4.0])
mk_b = [r for r in FLAT["b"] if is_marker(r) and inbox(r, B_BOX)]
ci_b = [r for r in FLAT["b"] if r["op"] == "S" and h_(r["bbox"]) < 0.05 and abs(r["lw"] - 1.7) < 0.01 and B_BOX[1] < r["bbox"][1] < B_BOX[3] and r["bbox"][0] > B_BOX[0] - 1]
worst_b = 0.0; nb = 0; span_b = math.log(DR["b"]["xlim"][1]) - math.log(DR["b"]["xlim"][0])
for m in mk_b:
    cx, cy = cen(m["bbox"]); d = min([dd for dd in DR["b"]["drawn"] if dd["estimable"]], key=lambda dd: abs(dd["y_page"] - cy) + abs(dd["x_page"] - cx) / 100); assert abs(d["y_page"] - cy) < 0.3, (d, cy)
    hr = PC[d["key"]]["hr_by_lag"][str(d["lag"])]; lo, hi = PC[d["key"]]["ci_by_lag"][str(d["lag"])]
    ci = min(ci_b, key=lambda r: abs(cen(r["bbox"])[1] - cy) + abs(cen(r["bbox"])[0] - cx) / 100)
    dev = max(abs(math.log(fb_(cx)) - math.log(hr)), abs(math.log(fb_(ci["bbox"][0])) - math.log(lo)), abs(math.log(fb_(ci["bbox"][2])) - math.log(hi))) / span_b * 100
    worst_b = max(worst_b, dev); nb += 1; agree.append(f"b  {PC[d['key']]['condition']:42s} lag {d['lag']}  file {hr:.3f} ({lo:.3f}-{hi:.3f})  read {fb_(cx):.3f} ({fb_(ci['bbox'][0]):.3f}-{fb_(ci['bbox'][2]):.3f})  dev {dev:.3f}%")
check(nb == 65 and worst_b <= 0.5, f"panel b: {nb} markers and CI ends read back through the tick marks agree with lag_ladder.json: max deviation {worst_b:.3f}% of the axis span")
sub_b = sorted({round(abs(d1["y_page"] - d2["y_page"]), 2) for d1 in DR["b"]["drawn"] for d2 in DR["b"]["drawn"] if d1["key"] == d2["key"] and abs(d1["lag"] - d2["lag"]) > 0 and d1["lag"] < d2["lag"]})
check(abs(min(sub_b) - 6.44) < 0.01 and abs(DR["b"]["pitch"] - 27.0) < 0.01, f"panel b: the four lag sub-rows sit 6.44 pt apart as on V30 (row pitch {DR['b']['pitch']} pt, V30 32.2; marker sizes unchanged)")
b1_n = sorted([r["bbox"] for r in FLAT["c"] if r["op"] == "f" and r["fill"] == LADDER[0] and h_(r["bbox"]) > 50], key=lambda b: (round(b[1]), b[0])); BOXES_N = [(b[0], b[1], b[0] + 4 * w_(b), b[3]) for b in b1_n]
mk_c = [r for r in FLAT["c"] if r["op"] == "B" and r["fill"] == INK and r["stroke"] == "#ffffff" and r["npts"] > 12 and h_(r["bbox"]) < 12]
rib_n = [r for r in FLAT["c"] if r["op"] == "f" and r["fill"] == INK and abs(r["alpha"] - 0.16) < 0.01 and h_(r["bbox"]) > 2]
yt_marks = [r for r in FLAT["c"] if r["op"] in ("S", "B") and h_(r["bbox"]) < 0.05 and 2 < w_(r["bbox"]) < 8]
NS_c = [s for s in TXT["c"] if s["text"].strip()]
yt_labels = [s for s in NS_c if abs(s["size"] - T14) < 0.01 and re.fullmatch(r"[\d.]+", s["text"])]
titles_c = [s for s in NS_c if abs(s["size"] - T15) < 0.01 and s["text"][0].isalpha() and s["text"] not in V30_CAPS_C]
check(len(BOXES_N) == 10 and len(mk_c) == 40 and len(rib_n) == 10 and len(titles_c) == 10, f"panel c: {len(BOXES_N)} boxes ({DR['c']['r54']['grid'][0]} x {DR['c']['r54']['grid'][1]}, {DR['c']['r54']['box'][0]:.0f} x {DR['c']['r54']['box'][1]:.0f} pt), {len(mk_c)} ink markers, {len(rib_n)} ribbons, {len(titles_c)} titles on the flat page")
worst_c = 0.0; nc = 0; star_dev = []; YT_PAD = DR["c"]["style"]["ytick_pad"]
for bx in BOXES_N:
    t = min(titles_c, key=lambda s: abs((s["bbox"][0] + s["bbox"][2]) / 2 - (bx[0] + bx[2]) / 2) + abs(s["origin"][1] - bx[1])); cond = t["text"]
    marks = sorted([r for r in yt_marks if abs(r["bbox"][2] - bx[0]) < 0.5 and bx[1] - 1 < r["bbox"][1] < bx[3] + 1], key=lambda r: r["bbox"][1])
    pts = []
    for mrow in marks:
        yl = min(yt_labels, key=lambda s: abs((s["bbox"][1] + s["bbox"][3]) / 2 - mrow["bbox"][1]) + abs(s["bbox"][2] - (bx[0] - YT_PAD)) / 10); pts.append((mrow["bbox"][1], float(yl["text"])))
    assert len(pts) >= 2, (cond, pts)
    ys, vs = zip(*pts); bcoef, acoef = np.polyfit(ys, vs, 1); fy = lambda y: acoef + bcoef * y
    v = G[cond]; mk = sorted([r for r in mk_c if bx[0] < cen(r["bbox"])[0] < bx[2] and bx[1] - 1 < cen(r["bbox"])[1] < bx[3] + 1], key=lambda r: r["bbox"][0])
    rib = [r for r in rib_n if bx[0] < cen(r["bbox"])[0] < bx[2] and bx[1] - 1 < cen(r["bbox"])[1] < bx[3] + 1]
    span_v = (fy(bx[1]) - fy(bx[3])); hrs = [1.0] + [v[b]["hr"] for b in ("1-5%", "5-10%", ">10%")]
    dev = max(abs(fy(cen(m["bbox"])[1]) - h) for m, h in zip(mk, hrs)) / span_v * 100
    rb = rib[0]["bbox"]; lo_min = min(v[b]["lo"] for b in ("1-5%", "5-10%", ">10%")); hi_max = max(v[b]["hi"] for b in ("1-5%", "5-10%", ">10%"))
    dev_r = max(abs(fy(rb[3]) - lo_min), abs(fy(rb[1]) - hi_max)) / span_v * 100
    worst_c = max(worst_c, dev, dev_r); nc += len(mk)
    st = [s for s in NS_c if s["text"] == "***" and bx[0] < s["origin"][0] < bx[2] and bx[1] - 12 < s["origin"][1] < bx[3]]
    if st:
        sty = DR["c"]["style"]; hi4 = v[">10%"]["hi"]; dev_s = abs(fy(st[0]["origin"][1] + sty["star_dy"]) - hi4) / span_v * 100; star_dev.append((cond, round(dev_s, 3)))
    agree.append(f"c  {cond:42s} file {' '.join(f'{h:.3f}' for h in hrs)}  read {' '.join(f'{fy(cen(m['bbox'])[1]):.3f}' for m in mk)}  ribbon lo {fy(rb[3]):.3f} (file {lo_min:.3f}) hi {fy(rb[1]):.3f} (file {hi_max:.3f})  dev {max(dev, dev_r):.3f}%")
check(nc == 40 and worst_c <= 0.5, f"panel c: {nc} markers and 10 ribbon extremes read back through each panel's y ticks agree with results_v2.json graded: max deviation {worst_c:.3f}% of the panel's span")
check(len(star_dev) == 10 and max(d[1] for d in star_dev) <= 0.5, f"panel c asterisks ride on the top band's upper CI at the round-30 offsets (max deviation {max(d[1] for d in star_dev):.3f}% of span)")
check([p["condition"] for p in DR["c"]["panels"]] == SEL["selected_labels"] and all(BOXES_N[i][1] <= BOXES_N[i + 1][1] + 0.01 and (BOXES_N[i][1] < BOXES_N[i + 1][1] - 1 or BOXES_N[i][0] < BOXES_N[i + 1][0]) for i in range(9)),
      "panel c: the ten boxes keep the selection file's order, row-major in the 3 x 4 grid (the tenth box alone in the last row)")
open(os.path.join(VER, "values_agreement.txt"), "w").write("\n".join(agree) + "\n")

# ------------------------------------------------------------------ 5 sizes on the flat panels = the drawn records; glyph overlaps; text vs marks; clearances
all_recs = [t for p in "abc" for t in DR[p]["texts"]]
sz_flat = Counter(round(s["size"], 2) for p in "abc" for s in TXT[p] if s["text"].strip()); sz_rec = Counter(round(t["size"], 2) for t in all_recs)
check(set(sz_flat) == set(sz_rec) == {T13, T14, T15}, f"type sizes on the flat panels {sorted(sz_flat)} = the drawn records' sizes {sorted(sz_rec)} = the three classes 13, 14, 15 (letters and title 18 by the composer)")
S1g = R.with_glyph_boxes([dict(s) for s in STR1]); S0g = R.with_glyph_boxes([dict(s) for s in STR0])
ov_new = R.glyph_overlaps(S1g); ov_old = R.glyph_overlaps(S0g)
check(not ov_new, f"no text glyph box (cap height to descender on the gs baseline) overlaps another on the composed sheet: {len(ov_new)} pairs (V30 at its 10.44 pt pitch: {len(ov_old)})", f"{ov_new[:8]}")
# text vs data marks on the flat panels: glyph boxes against markers, CI lines, dashed rules, panel c white lines
hits = []; n_checked = 0
for p in "abc":
    marks = [r for r in FLAT[p] if is_marker(r) or (r["op"] == "S" and abs(r["lw"] - 1.7) < 0.01 and h_(r["bbox"]) < 0.05) or (r["op"] == "S" and r["dash"] and (h_(r["bbox"]) > 50 or w_(r["bbox"]) > 50)) or (r["op"] == "S" and r["stroke"] == "#ffffff")]
    for s in TXT[p]:
        if not s["text"].strip(): continue
        n_checked += 1
        if abs(s["dir"][1]) > 0.5: tb = s["bbox"]                                   # rotated label: the rawdict box
        else: tb = [s["bbox"][0], s["origin"][1] - R.CAP * s["size"], s["bbox"][2], s["origin"][1] + R.DESC * s["size"]]
        for r in marks:
            mb = r["bbox"]
            if min(tb[2], mb[2]) - max(tb[0], mb[0]) > 0.2 and min(tb[3], mb[3]) - max(tb[1], mb[1]) > 0.2: hits.append((p, s["text"], r["op"], [round(v, 1) for v in mb]))
check(not hits, f"no text glyph box touches a data mark (markers, interval lines, dashed rules, panel c lines) on the three flat panels (checked {n_checked} spans)", f"{hits[:6]}")
LABX = DR["a"]["r54"]["lab_x"]
labs_a = [s for s in TXT["a"] if abs(s["origin"][0] - LABX) < 0.5 and abs(s["size"] - T14) < 0.01]
lab_end = max(s["bbox"][2] for s in labs_a); hdr_end = max(s["bbox"][2] for s in labs_a if "Bold" in s["font"])
own = []
for m in DR["a"]["marks"]:
    lt = DR["a"]["label_short"].get(m["condition"], m["condition"]); s = [t for t in labs_a if t["text"] == lt]; assert len(s) == 1, lt
    own.append(m["x_lo"] - s[0]["bbox"][2])
hr_left = min(s["bbox"][0] for s in TXT["a"] if HRQ.match(s["text"]))
check(ax["ref_x"] - lab_end >= 3.5 and min(own) >= 2.5 and hr_left - max(m["x_hi"] for m in DR["a"]["marks"]) >= 7.0,
      f"panel a clearances: longest label or header to the dashed reference line {ax['ref_x'] - lab_end:.2f} pt (headers {ax['ref_x'] - hdr_end:.2f}), every label to its own interval at least {min(own):.2f} pt, HR column to the widest interval end {hr_left - max(m['x_hi'] for m in DR['a']['marks']):.2f} pt")
gap_a = min((ax["tick_x"][i + 1] - ax["tick_x"][i]) - (text_width_pt(None, f"{ax['ticks'][i]:g}", T14) + text_width_pt(None, f"{ax['ticks'][i + 1]:g}", T14)) / 2 for i in range(6))
check(gap_a >= 2.5, f"panel a: adjacent x tick labels clear each other by at least {gap_a:.2f} pt (0.6 and 0.8 at 14 pt on a {ax['w_ax']:.0f} pt axis)")

# ------------------------------------------------------------------ 6 regions, printed values csv
check(all(CLOG["panels"][p]["region"][0] <= CLOG["panels"][p]["content_bbox"][0] and CLOG["panels"][p]["content_bbox"][2] <= CLOG["panels"][p]["region"][2] and CLOG["panels"][p]["region"][1] <= CLOG["panels"][p]["content_bbox"][1] and CLOG["panels"][p]["content_bbox"][3] <= CLOG["panels"][p]["region"][3] for p in "abc"),
      f"every panel's content stays inside its clip region (a/b split x {SPLIT}, b/c split y {BC_Y}): a right {CLOG['panels']['a']['content_bbox'][2]:.2f}, b left {CLOG['panels']['b']['content_bbox'][0]:.2f} bottom {CLOG['panels']['b']['content_bbox'][3]:.2f}, c left {CLOG['panels']['c']['content_bbox'][0]:.2f} top {CLOG['panels']['c']['content_bbox'][1]:.2f} right {CLOG['panels']['c']['content_bbox'][2]:.2f}")
recs_out = list(all_recs)
recs_out.append(dict(text="Figure 2", x=14.2, baseline=18.0, ha="left", size=T18, source_file="static:V13", source_key="title (round 54: 18 pt on baseline 18)", source_value="Figure 2", rule="text", panel=""))
for Lt, (x, y) in CLOG["letters"].items():
    recs_out.append(dict(text=Lt, x=x, baseline=y, ha="left", size=T18, source_file="static:V13", source_key="panel letter (round 54 position and size)", source_value=Lt, rule="text", panel=Lt))
json.dump(dict(page=[PAGE_W, PH], split_x=SPLIT, bc_y=BC_Y, plus=PLUS, panels=DR, records=recs_out), open(os.path.join(VER, "Main_Fig2_drawn.json"), "w"), indent=1)
n_rows = 0; n_no = 0; pv_rows = []
for p in "abc":
    tmp = os.path.join(WORK, f"printed_values_{p}.csv")
    nr, nn = printed_values_csv(os.path.join(WORK, f"panel_{p}.pdf"), DR[p]["texts"], tmp, static_from=None)
    prow = list(csv.DictReader(open(tmp)))
    # PyMuPDF joins neighbouring spans on one baseline into one string when their gap is small against the size ('0.6 0.8' at 14 pt on the
    # panel a axis, '5-10 >10' under every panel c box): such a page string equal to the space-join of unmatched drawn records on the same
    # baseline, in x order, is accounted for by those records (match 'yes_merged_by_extractor'); nothing else is forgiven
    no_page = [r for r in prow if r["match"] == "NO" and (r["note"].startswith("string not in the drawn records") or r["note"].startswith("digit-bearing string"))]
    no_rec = [r for r in prow if r["match"] == "NO" and r["note"].startswith("drawn record not found")]
    for pg in no_page:
        y = float(pg["y"]); cands = sorted([r for r in no_rec if abs(float(r["y"]) - y) <= 1.5], key=lambda r: float(r["x"]))
        parts = pg["string"].split(" ")
        for k in range(len(cands) - len(parts) + 1):
            grp = cands[k:k + len(parts)]
            if [g["string"] for g in grp] == parts:
                pg["match"] = "yes_merged_by_extractor"; pg["note"] = "page string = the space-join of the drawn records " + " + ".join(repr(g["string"]) for g in grp) + " on one baseline (PyMuPDF span merge)"
                for g in grp: g["match"] = "yes_merged_by_extractor"; g["note"] = f"read back inside the merged page string {pg['string']!r}"; no_rec.remove(g)
                break
    nn = sum(1 for r in prow if r["match"] == "NO"); n_rows += nr; n_no += nn
    pv_rows += prow
with open(os.path.join(VER, "Main_Fig2_printed_values.csv"), "w", newline="") as f:
    wr = csv.DictWriter(f, fieldnames=list(pv_rows[0].keys())); wr.writeheader(); wr.writerows(pv_rows)
    for t in ("Figure 2", "a", "b", "c"): wr.writerow({k: "" for k in pv_rows[0].keys()} | dict(string=t, source_file="static:V13", key="stamped by the composer (title, letters)", match="static"))
check(n_no == 0, f"printed values: {n_rows} strings on the three flat panels accounted for (drawn record with source), {n_no} unaccounted (NO); title and letters stamped by the composer")

# ------------------------------------------------------------------ 7 renders and crops (Ghostscript), the flat pdf
png150 = os.path.join(LANE, "Main_Fig2_150dpi.png"); st = R.gs_render(NEW, png150, 150, "gs150_new"); info(f"gs 150 dpi render: {st}")
im = Image.open(png150); check(im.size[0] > 1000 and im.size[1] > 1000, f"150 dpi render written ({im.size[0]} x {im.size[1]} px)")
png200n = os.path.join(WORK, "new_200dpi.png"); R.gs_render(NEW, png200n, 200, "gs200_new")
CROPS_DEF = [("01_title_letter_a_key", (0, 0, 495, 130)), ("02_a_first_blocks", (0, 100, 495, 330)), ("03_a_ventricular_row_and_columns", (0, 130, 495, 200)),
             ("04_a_foot_ticks_caption", (100, 1040, 495, 1107)), ("05_b_letter_key", (485, 20, 968.66, 80)), ("06_b_rows_top", (485, 70, 968.66, 250)), ("07_b_controls_and_foot", (485, 380, 968.66, 608)),
             ("08_c_row2_titles", (485, 740, 968.66, 830)), ("09_c_last_row_key_xtitle", (485, 990, 968.66, 1107)), ("10_c_ylabel", (485, 700, 560, 1000)), ("11_a_neuro_msk_blocks", (0, 640, 495, 900))]
for name, box in CROPS_DEF: R.crop(png200n, box, 200, os.path.join(CROPS, f"{name}_NEW_200dpi.png"), f"NEW (round 54) {name}")
check(all(os.path.exists(os.path.join(CROPS, f"{n}_NEW_200dpi.png")) for n, _ in CROPS_DEF), f"{len(CROPS_DEF)} crops of the new sheet at 200 dpi written to verify/crops")
FLATPDF = os.path.join(WORK, "Main_Fig2_flat.pdf"); stf = R.gs_flatten(NEW, FLATPDF, "gs_flatten"); info(f"gs pdfwrite flatten: {stf}")
fn, fw, fh = R.gs_pdfinfo(FLATPDF); SF = R.spans_xml(FLATPDF, os.path.join(WORK, "flat_txtwrite.xml"), "txt_flat"); STRF = R.strings_r54(SF, SPLIT, BC_Y)
cf = Counter((s["text"].replace(" ", ""), s["size"]) for s in STRF); c1 = Counter((s["text"].replace(" ", ""), s["size"]) for s in STR1)   # spaces ignored: pdfwrite drops the space of the TextWriter title (Figure2)
check(fn == 1 and abs(fw - PAGE_W) < 0.05 and abs(fh - PH) < 0.05 and cf == c1, f"work/Main_Fig2_flat.pdf (gs pdfwrite 1.7 prepress, no downsampling, colours unchanged, fonts embedded): one page {fw} x {fh}, the same {sum(cf.values())} (string, size) pairs as the composed sheet (spaces ignored: the flattened title reads Figure2)", f"flat-not-sheet {sorted((cf - c1).elements())[:5]} sheet-not-flat {sorted((c1 - cf).elements())[:5]}")

# ------------------------------------------------------------------ 8 CHANGES csv
A54, B54, C54 = DR["a"]["r54"], DR["b"]["r54"], DR["c"]["r54"]
rows = [dict(sheet="Main_Fig2", panel="all", item="every text element", old="V30 sizes 10, 10.5, 11, 11.8, 12, 12.4, 14, 14.2 pt", new="14 ticks, row labels, keys, values; 15 axis titles, column headers, panel titles; 13 'no estimate' and the asterisk key; 18 letters and title", kind="size (no word or number change)", x="", y="", note="round 54: main figures print at 54 percent, +3 pt on the round-49 sizes"),
        dict(sheet="Main_Fig2", panel="all", item="page box", old="968.66 x 1184.76", new=f"968.66 x {PH}", kind="layout", x="", y="", note="width unchanged, height under the 1110 pt width-bound limit of the manuscript page"),
        dict(sheet="Main_Fig2", panel="all", item="arrangement", old="a top-left, b top-right, c full width below", new="L-shape: a one column at the left full height, b top-right, c right below b", kind="layout", x=SPLIT, y=BC_Y, note="a two-column forest plus panel b cannot fit 1260 pt (PLAN.md)"),
        dict(sheet="Main_Fig2", panel="a", item="organ headers", old="left band, two-line headers", new="bold row above each block, one line", kind="layout (declared joins)", x=A54["lab_x"], y="", note="'Metabolic,'+'kidney and liver', 'Infection,'+'blood and other', 'Neurological'+'and psychiatric', 'Musculoskeletal'+'and sensory', 'Negative'+'controls' set on one line, same words"),
        dict(sheet="Main_Fig2", panel="a", item="row pitch, block gap", old="10.44, 7.92", new=f"{DR['a']['layout']['pitch']}, {DR['a']['layout']['block_gap']}", kind="layout", x="", y="", note="14 pt glyph box 12.99 pt"),
        dict(sheet="Main_Fig2", panel="a", item="x axis width", old="144.5", new=f"{ax['w_ax']}", kind="layout", x=A_BOX[0], y="", note="same range 0.58-4.189 and ticks 0.6 to 4; the 14 pt tick labels 0.6 and 0.8 clear each other"),
        dict(sheet="Main_Fig2", panel="b", item="row pitch", old="32.2 (sub-rows 6.44)", new=f"{DR['b']['pitch']} (sub-rows {DR['b']['sub_row_pt']})", kind="layout", x="", y="", note="marker sizes unchanged, the between-row gap 12.9 -> 7.7 pt, alternate shading kept"),
        dict(sheet="Main_Fig2", panel="b", item="key", old="two columns x three rows below the letter", new="two rows beside the letter, 3 + 2 entries, same order", kind="layout", x=DR["b"]["key_rows"][0]["tx"], y=DR["b"]["key_rows"][0]["ty"], note=""),
        dict(sheet="Main_Fig2", panel="b", item="x axis width, caption anchor", old="306, centred under the axis", new=f"{B_BOX[2] - B_BOX[0]:.1f}, {B54['caption_rule']}", kind="layout", x=B54["caption_x"], y="", note="the 15 pt caption (388 pt) centred under the axis would leave the page"),
        dict(sheet="Main_Fig2", panel="c", item="grid and box", old="2 x 5, 143.8 x 96.7", new=f"{C54['grid'][0]} x {C54['grid'][1]}, {C54['box'][0]:.0f} x {C54['box'][1]:.0f}", kind="layout", x="", y="", note="row-major order unchanged; y limits, ticks, ribbons and stars follow the V30 rules on the data"),
        dict(sheet="Main_Fig2", panel="c", item="titles", old="left-aligned at the box", new="centred over the box", kind="layout", x="", y="", note="a 15 pt title left-aligned would run off the sheet in the last column"),
        dict(sheet="Main_Fig2", panel="c", item="band key and asterisk key", old="one row under the x title", new="two rows in the free slots of the last row", kind="layout", x=C54["free_x0"], y=C54["key_y"][0], note=""),
        dict(sheet="Main_Fig2", panel="title, letters", item="'Figure 2' and a, b, c", old="14 pt at (14.2, 14), (14, 40), (502, 40), (14, 803)", new=f"18 pt at (14.2, 18), {CLOG['letters']}", kind="size and position", x=14.2, y=18.0, note="title cap top about 5 pt inside the page")]
R.write_changes(os.path.join(VER, "CHANGES_Main_Fig2.csv"), rows); info(f"CHANGES_Main_Fig2.csv: {len(rows)} rows (size and layout only); declared delta of the word and numeric multisets: none (five header line joins)")
ok = ck.write(os.path.join(VER, "checks.txt"))
rec.update(census=dict(same_strings=CEN["same_strings"], sizes_ok=CEN["sizes_ok"], sizes_old=CEN["sizes_old"], sizes_new=CEN["sizes_new"], joins=CEN["joins"], min_size_new=CEN["min_size_new"]),
           values=dict(worst_a_pct=round(worst_a, 4), worst_b_pct=round(worst_b, 4), worst_c_pct=round(worst_c, 4)), overlaps=dict(new=ov_new, v30=len(ov_old)), text_vs_marks=hits, png150=png150, png150_sha256=R.sha256(png150),
           flat_pdf=FLATPDF, flat_sha256=R.sha256(FLATPDF), result="ALL PASS" if ok else "FAIL")
json.dump(rec, open(os.path.join(VER, "verify_record.json"), "w"), indent=1, default=str)
sys.exit(0 if ok else 1)
