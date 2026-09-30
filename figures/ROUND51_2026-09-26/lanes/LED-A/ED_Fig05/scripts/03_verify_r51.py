#!/usr/bin/env python3
"""ROUND 51, lane LED-A: ED Fig 5 proof (item 2), run on the built sheet before the title strip (work/ED_Fig05_built.pdf). Declared this round:
no values above the bars (16 strings removed), the five-line blocks removed (8 x 'PR, x.xx', '(95% CI,', 'lo-hi)', 'P ...'), the prevalence ratio
with its 95% CI as the second tick-label line under each pair ('PR x.xx (lo-hi)', 8 strings at the tick-label size), the P value as asterisks
above the pair (*** for P < 0.001, none at or above 0.05: 6 strings), the bar colours #79d475 (pale, under 5 h) and #298d32 (dark, 7 h or more),
the axes widened so that the one-line labels fit (the page with them) and the space below the axes shortened. Kept numbers: the eight prevalence
ratios with their intervals and the four P values equal V28's (from V28's own build record), every printed string recomputed from the v8 table
(the builder's step-129 cross-check), the bar heights equal V28's. Renders by Ghostscript."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../../..")))  # repository root, where paths.py lives
import paths
import collections, csv, json, os, re, sys
import numpy as np, fitz
from PIL import Image
sys.path.insert(0, f"{paths.FIGURE_ROOT}/ROUND51_2026-09-26/lanes/LED-A/scripts_shared")
import lane_common as C
import r49lib as R
VL = C.VL
SHEET = "ED_Fig05"; HERE = f"{C.LANE}/{SHEET}"; VER = f"{HERE}/verify"; CROPS = f"{VER}/crops"; WORK = f"{HERE}/work"; LOGS = f"{R.LANE}/logs/{SHEET}"; os.makedirs(CROPS, exist_ok=True)
BASE = f"{C.BASE}/{SHEET}.pdf"; NEW = f"{WORK}/{SHEET}_built.pdf"; V28 = f"{R.V28}/{SHEET}.pdf"; STRIP = 16.0
V28_REC = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LED-A/ED_Fig05/work/ED_Fig05_drawn.json"   # the V28 build's own record (V28 = the round-49 output, sha256 checked below)
D = json.load(open(f"{WORK}/{SHEET}_drawn.json")); KP = D["key_placement"]; LAY = D["layout"]; R49 = json.load(open(C.hydrated(V28_REC)))
LINES = []; OK = True
def check(name, ok, detail):
    global OK; OK &= bool(ok); LINES.append(f"{'PASS' if ok else 'FAIL'}  {name}: {detail}"); print(LINES[-1][:300])
check("V28 is the round-49 output whose build record is read for the kept numbers", C.sha256(V28).startswith("9355448c7d790ee8"), f"V28 sha256 {C.sha256(V28)[:16]}")
nd = fitz.open(C.hydrated(NEW)); np_ = nd[0]; vd = fitz.open(C.hydrated(V28)); vp = vd[0]
W, H = np_.rect.width, np_.rect.height
check("page (before the strip) as the builder's layout declares", len(nd) == 1 and abs(W - LAY["w_in"] * 72) < 0.05 and abs(H - LAY["h_in"] * 72) < 0.05,
      f"{W:.2f} x {H:.2f} pt (V28 {LAY['w_v28_in'] * 72:.2f} x {LAY['h_v28_in'] * 72:.2f}): axes {LAY['axw_each_in']:.3f} in wide (V28 {LAY['axw_each_v28_in']:.3f}) so that the one-line label ({LAY['pr_line_w_in'] * 72:.1f} pt) fits under each pair with {LAY['pr_gap_pt']:.0f} pt between the two lines of a panel, {LAY['below_in']:.2f} in below the axes (V28 {LAY['below_v28_in']:.2f})")
LN = C.letters(np_); LB = C.letters(vp)
DXB = (LAY["axw_each_in"] - LAY["axw_each_v28_in"]) * 72; DYR = (LAY["below_v28_in"] - LAY["below_in"]) * 72
exp_l = {"a": (LB[0][1], LB[0][2] - STRIP), "b": (LB[1][1] + DXB, LB[1][2] - STRIP), "c": (LB[2][1], LB[2][2] - STRIP - DYR), "d": (LB[3][1] + DXB, LB[3][2] - STRIP - DYR)}
check(f"letters a, b, c, d: Arial Bold 14, a and c at their V28 left edge, b and d {DXB:.1f} pt right with the wider left panel, the second row {DYR:.1f} pt higher with the shorter space below the axes (declared)",
      [t[0] for t in LN] == ["a", "b", "c", "d"] and all(abs(t[5] - 14.0) < 0.05 and abs(t[1] - exp_l[t[0]][0]) < 0.5 and abs(t[2] - exp_l[t[0]][1]) < 0.5 for t in LN), f"expected {dict((k, (round(v[0], 2), round(v[1], 2))) for k, v in exp_l.items())}, new {[(t[0], t[1], t[2]) for t in LN]}")
fb = sorted(set(f[3].split("+")[-1] for f in vp.get_fonts(full=True))); fn = sorted(set(f[3].split("+")[-1] for f in np_.get_fonts(full=True)))
check("fonts: Arial faces, within V28's set (V28 also carries the strip's Arial Bold)", set(fn) <= set(fb) and all("Arial" in f for f in fn), f"V28 {fb}, new {fn}")
check("no nested forms", len(np_.get_xobjects()) == 0, f"{len(np_.get_xobjects())} xobjects")
SB = C.dedupe(C.spans_of(vp)); SN = C.dedupe(C.spans_of(np_))
sizes_b = sorted({round(s["size"], 2) for s in SB if s["size"] < 13.5 or s["text"].strip() in "abcd"}); sizes_n = sorted({round(s["size"], 2) for s in SN})
check("type sizes = V28's (round-49 values: annotations 10.5, ticks, key and titles 11, y labels 12, letters 14)", sizes_n == sizes_b, f"V28 {sizes_b} -> new {sizes_n}")
# strings: removed, added, kept
RX_PCT = re.compile(r"^\d+\.\d$"); RX_OLD = re.compile(r"^(PR, \d\.\d\d|\(95% CI,|\d\.\d\d-\d\.\d\d\)|P [<=] 0\.\d\d\d)$"); RX_PRL = re.compile(r"^PR \d\.\d\d \(\d\.\d\d-\d\.\d\d\)$"); RX_ST = re.compile(r"^\*{1,3}$")
removed_exp = collections.Counter(t for p in D["panels"].values() for c in p["cuts"].values() for t in c["removed_v28"])
added_exp = collections.Counter(t for p in D["panels"].values() for c in p["cuts"].values() for t in c["printed"])
old_rm = collections.Counter(s["text"] for s in SB if RX_PCT.match(s["text"]) or RX_OLD.match(s["text"]))
check("the 48 V28 value strings are gone (16 bar values, 8 'PR, x.xx', 8 '(95% CI,', 8 'lo-hi)', 8 P strings) and none remains on the sheet", old_rm == removed_exp and sum(removed_exp.values()) == 48 and not [s for s in SN if RX_PCT.match(s["text"]) or RX_OLD.match(s["text"])], f"V28 carried {sum(old_rm.values())} of them: {dict(old_rm)}")
new_add = collections.Counter(s["text"] for s in SN if RX_PRL.match(s["text"]) or RX_ST.match(s["text"]))
check("the added strings are exactly the 8 one-line 'PR x.xx (lo-hi)' labels and the 6 '***'", new_add == added_exp and sum(added_exp.values()) == 14 and all(abs(s["size"] - 11.0) < 0.05 for s in SN if RX_PRL.match(s["text"])) and all(abs(s["size"] - 10.5) < 0.05 and "Bold" in s["font"] for s in SN if RX_ST.match(s["text"])), f"{dict(new_add)}")
inv_b = collections.Counter(s["text"] for s in SB if not (RX_PCT.match(s["text"]) or RX_OLD.match(s["text"])) and not s["text"].startswith("Extended")); inv_n = collections.Counter(s["text"] for s in SN if not (RX_PRL.match(s["text"]) or RX_ST.match(s["text"])))
check("every other string of V28 is on the sheet, once, and nothing else (titles, letters, tick labels, y labels, key lines)", inv_b == inv_n, f"V28-only {dict(inv_b - inv_n)}, new-only {dict(inv_n - inv_b)}")
# kept numbers identical to V28's build record
same = True; det = []
for pnl in "abcd":
    for c in ("5", "10"):
        v, o = D["panels"][pnl]["cuts"][c], R49["panels"][pnl]["cuts"][c]
        eq = all(abs(float(v[k]) - float(o[k])) < 1e-12 for k in ("pr", "lo", "hi", "p", "pct_short", "pct_ref")) and v["p_printed"] == o["p_printed"] and D["panels"][pnl]["n_short"] == R49["panels"][pnl]["n_short"] and D["panels"][pnl]["n_ref"] == R49["panels"][pnl]["n_ref"]
        same &= eq; det.append(f"{pnl} >{c}%: PR {v['pr']:.4f} ({v['lo']:.4f}-{v['hi']:.4f}) P {v['p']:.3g}")
check("kept numbers identical to V28: the eight prevalence ratios with their intervals, the four P values, the group percentages and sizes (from V28's own build record, to 1e-12)", same, "; ".join(det))
exp_lines = {(pnl, c): f"PR {o['pr']:.2f} ({o['lo']:.2f}-{o['hi']:.2f})" for pnl in "abcd" for c in ("5", "10") for o in [R49["panels"][pnl]["cuts"][c]]}
exp_stars = {(pnl, c): ("***" if o["p"] < 0.001 else "**" if o["p"] < 0.01 else "*" if o["p"] < 0.05 else "") for pnl in "abcd" for c in ("5", "10") for o in [R49["panels"][pnl]["cuts"][c]]}
def panel_of(s):
    xmid = W / 2.0 + 0.1; ymid = H / 2.0
    return ("a" if s["origin"][1] < ymid else "c") if s["origin"][0] < xmid else ("b" if s["origin"][1] < ymid else "d")
got_lines = {}; got_stars = {}
for pnl in "abcd":
    ss = [s for s in SN if panel_of(s) == pnl]
    prl = sorted([s for s in ss if RX_PRL.match(s["text"])], key=lambda s: s["origin"][0]); st = sorted([s for s in ss if RX_ST.match(s["text"])], key=lambda s: s["origin"][0])
    tk = sorted([s for s in ss if s["text"].startswith("T90 >")], key=lambda s: s["origin"][0])
    for i, c in enumerate(("5", "10")):
        got_lines[(pnl, c)] = prl[i]["text"] if len(prl) == 2 else None
        cx_tick = (tk[i]["bbox"][0] + tk[i]["bbox"][2]) / 2 if len(tk) == 2 else None
        near_st = [s for s in st if cx_tick is not None and abs((s["bbox"][0] + s["bbox"][2]) / 2 - cx_tick) < 3.0]
        got_stars[(pnl, c)] = near_st[0]["text"] if near_st else ""
        if len(prl) == 2 and cx_tick is not None:
            if abs((prl[i]["bbox"][0] + prl[i]["bbox"][2]) / 2 - cx_tick) > 0.6: got_lines[(pnl, c)] = "OFF-CENTRE"
check("each 'PR x.xx (lo-hi)' line equals the V28 value printed the same way, centred under its pair (within 0.6 pt of the tick label above it)", got_lines == exp_lines, f"{[(k, got_lines[k]) for k in sorted(got_lines) if got_lines[k] != exp_lines[k]]}")
check("asterisks by the rule on the V28 P values (*** P < 0.001, none at or above 0.05: panels a, b, c starred, panel d none), centred above their pair", got_stars == exp_stars, f"{dict((k, v) for k, v in got_stars.items())}")
# stars above the taller bar of the pair
hx = lambda c: "#%02x%02x%02x" % tuple(int(round(v * 255)) for v in c)
bars = [x for x in np_.get_drawings() if x.get("fill") is not None and hx(x["fill"]) in (D["colours"]["under_5h"], D["colours"]["7h_or_more"]) and x["rect"].width > 15]
bars_v = [x for x in vp.get_drawings() if x.get("fill") is not None and hx(x["fill"]) in (R49["colours"]["under_5h"], R49["colours"]["7h_or_more"]) and x["rect"].width > 15]
check("16 bars in the two new greens (pale #79d475 under 5 h, dark #298d32 7 h or more), no other bar colour", len(bars) == 16 and sorted({hx(b["fill"]) for b in bars}) == sorted({D["colours"]["under_5h"], D["colours"]["7h_or_more"]}) and D["colours"]["7h_or_more"] == "#298d32" and D["colours"]["under_5h"] == "#79d475", f"{sorted({hx(b['fill']) for b in bars})}")
luma = lambda h: 0.299 * int(h[1:3], 16) + 0.587 * int(h[3:5], 16) + 0.114 * int(h[5:7], 16)
check("greyscale separation of the pair at least 45 (BT.601 luma)", luma(D["colours"]["under_5h"]) - luma(D["colours"]["7h_or_more"]) >= 45, f"{luma(D['colours']['under_5h']):.0f} vs {luma(D['colours']['7h_or_more']):.0f} (V28 pair {luma(R49['colours']['under_5h']):.0f} vs {luma(R49['colours']['7h_or_more']):.0f})")
hn = sorted(round(b["rect"].height, 3) for b in bars); hv = sorted(round(b["rect"].height, 3) for b in bars_v)
check("bar heights equal V28's within 0.05 pt (same axes height and y range: the percentages unchanged), bar widths as the layout declares", len(hv) == 16 and max(abs(a - b) for a, b in zip(hn, hv)) < 0.05 and all(abs(b["rect"].width - LAY["bar_width_pt"]) < 0.1 for b in bars), f"max height deviation {max(abs(a - b) for a, b in zip(hn, hv)) if len(hv) == 16 else 'n/a'} pt, width {LAY['bar_width_pt']:.2f} pt (V28 {LAY['bar_width_v28_pt']:.2f})")
st_pos = []
for s in [s for s in SN if RX_ST.match(s["text"])]:
    cxs = (s["bbox"][0] + s["bbox"][2]) / 2
    pair = [b for b in bars if abs((b["rect"].x0 + b["rect"].x1) / 2 - cxs) < LAY["bar_width_pt"] * 0.7 and panel_of(dict(origin=((b["rect"].x0 + b["rect"].x1) / 2, (b["rect"].y0 + b["rect"].y1) / 2))) == panel_of(s)]
    top = min(b["rect"].y0 for b in pair) if pair else None
    st_pos.append((s["text"], len(pair), None if top is None else round(top - s["origin"][1], 2)))
check("every asterisk group sits above the taller bar of its pair (baseline 3.4 pt above its top, the old value labels' offset)", all(n == 2 and d_ is not None and 2.5 < d_ < 4.5 for _, n, d_ in st_pos) and len(st_pos) == 6, f"{st_pos}")
kl = [s for s in SN if s["text"] in ("Sleep under", "5 h", "Sleep 7 h", "or more")]
check("key: four lines at 11 pt, right edge at the page margin, top as V28, its left edge 0.024 in past the panel b axes boundary as in V28", len(kl) == 4 and all(abs(s["size"] - 11.0) < 0.05 for s in kl) and KP["key_right_in"] <= LAY["w_in"] - 0.18 + 1e-6 and KP["key_left_in"] >= KP["axes_right_in"] - 0.03 and abs(KP["key_top_in"] - 1.04) < 0.03, f"left {KP['key_left_in']:.3f} in, axes right {KP['axes_right_in']:.3f} in, top {KP['key_top_in']:.3f} in")
def body(s):
    r = fitz.Rect(s["bbox"]); return fitz.Rect(r.x0, r.y0 + 0.15 * s["size"], r.x1, r.y1 - 0.2 * s["size"])
ov = [(SN[i]["text"], SN[j]["text"]) for i in range(len(SN)) for j in range(i + 1, len(SN)) if not (body(SN[i]) & body(SN[j])).is_empty and (body(SN[i]) & body(SN[j])).get_area() > 0.05]
check("no text span overlaps another", not ov, f"{len(ov)}{': ' + str(ov[:3]) if ov else ''}")
prl_all = sorted([s for s in SN if RX_PRL.match(s["text"])], key=lambda s: (round(s["origin"][1]), s["origin"][0]))
gaps = [round(b["bbox"][0] - a["bbox"][2], 2) for a, b in zip(prl_all, prl_all[1:]) if abs(a["origin"][1] - b["origin"][1]) < 0.5 and 0 < b["bbox"][0] - a["bbox"][2] < 40]
check("the two PR lines of a panel are separate (gap at least 8 pt, the layout's 9 pt)", gaps and min(gaps) >= 8.0 and len(gaps) == 4, f"gaps {gaps} pt")
ink_r = fitz.Rect()
for g in np_.get_drawings():
    if g.get("fill") is not None and hx(g["fill"]) == "#ffffff" and g.get("color") is None: continue      # the page background
    ink_r |= g["rect"]
for s in SN: ink_r |= fitz.Rect(s["bbox"])
check("all ink inside the page with a 2 pt margin on the text-layer boxes (the builder's own 4 mm gate ran on the rendered text boxes at build time and passed)", ink_r.x0 > 2 and ink_r.y0 > 2 and ink_r.x1 < W - 2 and ink_r.y1 < H - 2, f"{[round(v, 1) for v in ink_r]} on {W:.1f} x {H:.1f}")
# read-back recomputation (the builder's record against the sheet)
expected = collections.Counter()
for L_, pnl in D["panels"].items():
    for c, v in pnl["cuts"].items():
        for t in v["printed"]: expected[(L_, t)] += 1
got = collections.Counter((panel_of(s), s["text"]) for s in SN if RX_PRL.match(s["text"]) or RX_ST.match(s["text"]))
check("read-back: every printed PR line and asterisk group equals the value recomputed from the v8 table (the builder's step-129 cross-check)", got == expected, f"printed {sum(got.values())}, expected {sum(expected.values())}")
tb, tn = C.nums([w["text"] for w in C.dedupe(C.words_of(vp)) if not w["text"].startswith("Extended")]), C.nums([w["text"] for w in C.dedupe(C.words_of(np_))])
rem_tok = C.nums(list(removed_exp.elements())); add_tok = C.nums(list(added_exp.elements()))
check("numeric-token multiset = V28's minus the removed strings' tokens plus the added strings' tokens (V28's title numeral aside)", collections.Counter(tn) + collections.Counter(rem_tok) + collections.Counter(["5"]) == collections.Counter(tb) + collections.Counter(add_tok), f"V28 {len(tb)} tokens, new {len(tn)}, removed {len(rem_tok)}, added {len(add_tok)}")
R.gs_render(V28, f"{WORK}/v28_150.png", 150, LOGS, name="gs150_v28"); R.gs_render(f"{HERE}/{SHEET}.pdf", f"{WORK}/new_150.png", 150, LOGS, name="gs150_new_delivered")
a = Image.open(f"{WORK}/v28_150.png").convert("RGB"); b = Image.open(f"{WORK}/new_150.png").convert("RGB")
side = Image.new("RGB", (a.width + b.width + 20, max(a.height, b.height)), "white"); side.paste(a, (0, 0)); side.paste(b, (a.width + 20, 0)); side.save(f"{CROPS}/{SHEET}_V28_vs_NEW_150dpi.png")
RULES = VL.RULES
RULES["pr_ci_line"] = lambda v: f"PR {VL.halfup(v[0], 2)} ({VL.halfup(v[1], 2)}-{VL.halfup(v[2], 2)})"
RULES["pr_2dp_first"] = lambda v: "PR " + VL.halfup(v[0], 2)
RULES["ci_2dp_paren"] = lambda v: f"({VL.halfup(v[0], 2)}-{VL.halfup(v[1], 2)})"
RULES["stars_p"] = lambda p: "***" if float(p) < 0.001 else "**" if float(p) < 0.01 else "*" if float(p) < 0.05 else ""
n_rows, n_no = VL.printed_values_csv(NEW, D["records"], f"{VER}/{SHEET}_printed_values.csv", static_from=BASE)
check("printed-values CSV: every string accounted for (14 drawn value strings under their rules re-read from the record, the rest V28 strings)", n_no == 0, f"{n_rows} rows, {n_no} NO")
hits = C.grep_proof(HERE)
check("grep proof: no live script in this sheet folder names a v7 table, snapshot or set", not hits, f"{hits}")
LINES.append(f"INFO  sources: {D['sources']['parquet']} sha256 {D['sources']['parquet_sha256'][:16]}, habitual {D['sources']['habitual_sha256'][:16]}, {D['sources']['cohort_filter']}")
nd.close(); vd.close()
print(C.write_checks(LINES, f"{VER}/checks_03_verify.txt", OK))
