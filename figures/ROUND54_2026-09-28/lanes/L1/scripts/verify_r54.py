#!$T90_PY
"""Round 54 lane L1 verifier for the composed Main_Fig1 (bands a and 1f r54 + the L1 panels at 14 / 15 / 13 pt + 18-pt letters and title) against
V30's Main_Fig1 (ROUND52 NEW_FINAL_SET_V30). Ghostscript txtwrite census only (chars assembled into words and lines, rotated spans boxed vertically).
Checks: page box, bands as declared ready, word census with the declared size mapping, nothing below 12 pt, printed values re-derived from the numbers
files (ranking_v3.csv, measure_families.csv, dC_per_disease.csv, L2_fig1d_values.json) and every value fig1_build.py recorded, no overlapping text
boxes, the key placement rule (build record and census), marks preserved against the LFIG1E record (V30's panels), letters and title at 18 pt,
banned strings, the 150-dpi render. Writes verify/checks.txt ending RESULT ALL PASS."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import collections, hashlib, html, json, os, re, subprocess, sys
import numpy as np, pandas as pd
from pypdf import PdfReader
T90 = paths.FIGURE_ROOT; R54 = f"{paths.FIGURE_ROOT}/ROUND54_2026-09-28"; LANE = f"{R54}/lanes/L1"; SD = f"{LANE}/Main_Fig1"; VER = f"{SD}/verify"
V30 = f"{paths.FIGURE_ROOT}/ROUND52_2026-09-26/figures/NEW_FINAL_SET_V30/Main_Fig1.pdf"; REF = f"{R54}/lanes/LFIG1E/Main_Fig1"   # LFIG1E = V30's panels with the key moved (its drawn record holds V30's marks)
NEW = f"{SD}/Main_Fig1.pdf"; GS = paths.GS; CAP, DESC = 0.727, 0.21
sys.path.insert(0, f"{LANE}/scripts"); import common as C
lines_out = []; ok_all = True
def check(name, good, detail=""):
    global ok_all; ok_all &= bool(good); lines_out.append(f"{'PASS' if good else 'FAIL'} {name}" + (f": {detail}" if detail else "")); print(lines_out[-1][:600])
sha = lambda p: hashlib.sha256(open(C.hydrated(p), "rb").read()).hexdigest()

# ------------------------------------------------------------------ the census
def census(pdf, xml):
    subprocess.run([GS, "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=txtwrite", "-dTextFormat=0", f"-sOutputFile={xml}", pdf], check=True, capture_output=True, timeout=600)
    raw = open(xml, encoding="utf8", errors="replace").read(); words = []; chars = []
    for m in re.finditer(r'<span bbox="([^"]+)" font="([^"]+)" size="([^"]+)">(.*?)</span>', raw, re.S):
        bx = [float(v) for v in m.group(1).split()]; size = round(float(m.group(3)) * 2) / 2; bold = "Bold" in m.group(2)
        cs = [(cb, html.unescape(c).replace("\xa0", " ")) for cb, c in re.findall(r'<char bbox="([^"]+)" c="([^"]*)"/>', m.group(4))]   # txtwrite escapes <, the minus sign, the less-or-equal sign and a no-break space
        if abs(bx[0] - bx[2]) < 0.01 and len(cs) > 1:   # rotated span: every char carries the span box, the text runs from y0 up to y1: one word entry holding the whole string
            text = " ".join("".join(c for _, c in cs).split()); ytop, ybot = min(bx[1], bx[3]), max(bx[1], bx[3])
            words.append(dict(text=text, size=size, bold=bold, rot=True, box=[bx[0] - CAP * size, ytop, bx[0] + DESC * size, ybot], y=(ytop + ybot) / 2, x=bx[0]))
            continue
        for cb, c in cs:
            b = [float(v) for v in cb.split()]; chars.append(dict(x0=b[0], x1=b[2], y=b[1], c=c, size=size, bold=bold))
    chars.sort(key=lambda c: (c["bold"], c["size"], round(c["y"]), c["x0"]))
    cur = None
    def flush():
        if cur and cur["text"]:
            words.append(dict(text=cur["text"], size=cur["size"], bold=cur["bold"], rot=False, box=[cur["x0"], cur["y"] - CAP * cur["size"], cur["x1"], cur["y"] + DESC * cur["size"]], y=cur["y"], x=cur["x0"]))
    for c in chars:
        same = cur is not None and cur["bold"] == c["bold"] and cur["size"] == c["size"] and abs(cur["y"] - c["y"]) <= 1.0 and -1.5 <= c["x0"] - cur["x1"] <= 0.2 * c["size"]   # a backwards jump is another string
        if c["c"].isspace() or not same:
            flush(); cur = None if c["c"].isspace() else dict(text=c["c"], size=c["size"], bold=c["bold"], y=c["y"], x0=c["x0"], x1=c["x1"])
            if c["c"].isspace(): continue
            continue
        cur["text"] += c["c"]; cur["x1"] = max(cur["x1"], c["x1"])
    flush()
    return words
def lines_of(words):
    """words on one baseline at one size joined left to right while the gap stays under half a size (rotated words: one line per span)"""
    out = []; groups = collections.defaultdict(list)
    for w in words: groups[(w["rot"], w["bold"], w["size"], round(w["y"]) if not w["rot"] else round(w["x"]))].append(w)
    for key, ws in groups.items():
        ws.sort(key=lambda w: (w["x"] if not w["rot"] else -w["box"][3])); cur = [ws[0]]
        for a, b in zip(ws, ws[1:]):
            gap = (b["box"][0] - a["box"][2]) if not a["rot"] else (a["box"][1] - b["box"][3])
            if gap < 0.5 * a["size"]: cur.append(b)   # a word space is 0.28 of the size, neighbouring strings on a baseline sit at least 6 pt apart
            else: out.append(cur); cur = [b]
        out.append(cur)
    return [dict(text=" ".join(w["text"] for w in ws), size=ws[0]["size"], bold=ws[0]["bold"], rot=ws[0]["rot"], y=ws[0]["y"], box=[min(w["box"][0] for w in ws), min(w["box"][1] for w in ws), max(w["box"][2] for w in ws), max(w["box"][3] for w in ws)]) for ws in out]
os.makedirs(VER, exist_ok=True)
W_NEW = census(NEW, f"{VER}/new_txt.xml"); W_OLD = census(V30, f"{VER}/v30_txt.xml")
LOG = json.load(open(f"{SD}/work/compose_log.json")); SL = LOG["slots"]
def region(w, new):
    y = w["y"]
    if new: return "a" if y < SL["panels"][0] else ("panels" if y < SL["f"][0] else "f")
    return "a" if y < 259.736 else ("panels" if y < 948.868 else "f")
for w in W_NEW: w["region"] = region(w, True)
for w in W_OLD: w["region"] = region(w, False)
L_NEW = lines_of([w for w in W_NEW if w["region"] == "panels"])

# 1 page box, bands
w_pt, h_pt = (float(v) for v in (PdfReader(NEW).pages[0].mediabox.width, PdfReader(NEW).pages[0].mediabox.height)); w30, h30 = (float(v) for v in (PdfReader(V30).pages[0].mediabox.width, PdfReader(V30).pages[0].mediabox.height))
check("page box declared: 968.66 x H, H at most 1260 (gate) and at most 1110 (aim)", abs(w_pt - 968.66) < 0.05 and h_pt <= 1260.0 and h_pt <= 1110.0, f"new {w_pt:.3f} x {h_pt:.3f}, V30 {w30:.3f} x {h30:.3f}, band a {LOG['band_a']['height']} (V30 231.114), band f {LOG['band_f']['height']}, shift {LOG['shift']}")
ready = open(f"{R54}/lanes/LSKETCH/BANDS_READY.txt").read(); declared = dict(re.findall(r"(band_\w+_r54\.pdf)\t[^\t]+\tsha256=([0-9a-f]{64})", ready))
check("the composed bands are the ones BANDS_READY.txt declares (sha256), not the round-49 stand-ins", not LOG["provisional_r49_bands"] and declared.get("band_a_r54.pdf") == LOG["band_a"]["sha256"] and declared.get("band_1f_r54.pdf") == LOG["band_f"]["sha256"] and sha(LOG["band_a"]["path"]) == LOG["band_a"]["sha256"] and sha(LOG["band_f"]["path"]) == LOG["band_f"]["sha256"], f"band a {LOG['band_a']['sha256'][:12]}, band f {LOG['band_f']['sha256'][:12]}")
check("panels page is the lane's current build (sha256)", sha(f"{SD}/work/panels_bcde_m.pdf") == LOG["panels"]["sha256"], LOG["panels"]["sha256"][:12])

# 2 census against V30 with the declared size mapping
TILE = {"Cirrhosis", "Heart", "failure", "Type", "2", "diabetes", "Hypertension", "Respiratory", "Obesity", "COPD", "Pneumonia"}
HEADER = {"Gain", "(95%", "CI)"}; LOGRANK = {"Log-rank", "P", "<", "0.001"}
def expect(w):
    s = w["size"]
    if w["bold"]: return 18.0 if s == 14.0 else s
    if s == 11.0: return 14.0
    if s == 12.0: return 14.0 if w["text"] in TILE else 15.0
    if s == 10.5: return 15.0 if w["text"] in HEADER else (13.0 if w["text"] in LOGRANK else 14.0)
    return None
def wc(ws, f=lambda w: w["size"]): return collections.Counter((t, f(dict(w, text=t))) for w in ws for t in w["text"].split())
exp = wc([w for w in W_OLD if w["region"] == "panels"], expect); got = wc([w for w in W_NEW if w["region"] == "panels"])
miss = exp - got; extra = got - exp
check("panels census: every V30 word present at its round-54 size (11 to 14, 12 to 15, tile titles 12 to 14, values 10.5 to 14, column header 10.5 to 15, log-rank 10.5 to 13, letters 14 to 18) and nothing added",
      not miss and not extra, f"{sum(exp.values())} V30 words, {sum(got.values())} new words, missing {dict(miss)}, extra {dict(extra)}")
for reg, name in (("a", "band a"), ("f", "band f")):
    eo = collections.Counter(t for w in W_OLD if w["region"] == reg for t in w["text"].split()); eg = collections.Counter(t for w in W_NEW if w["region"] == reg for t in w["text"].split())
    hist = sorted(collections.Counter(w["size"] for w in W_NEW if w["region"] == reg).items())
    check(f"{name} census: the same words as V30 (line breaks ignored) at 14, 15 or 18 pt (LSKETCH's declared sizes)", eo == eg and all(s in (14.0, 15.0, 18.0) for s, _ in hist), f"size histogram {hist}, missing {dict(eo - eg)}, extra {dict(eg - eo)}")
sizes_p = sorted({w["size"] for w in W_NEW if w["region"] == "panels" and not w["bold"]}); minsz = min(w["size"] for w in W_NEW)
check("no text below 12 pt anywhere on the sheet, panel text only at 13, 14, 15", minsz >= 12.0 and sizes_p == [13.0, 14.0, 15.0], f"smallest {minsz}, panel sizes {sizes_p}")

# 3 printed values re-derived from the numbers files
RK = pd.read_csv(C.hydrated(f"{C.NUM}/ranking_v3.csv"), comment="#"); FAM = pd.read_csv(C.hydrated(f"{C.NUM}/measure_families.csv")).set_index("feature")["family"]
RK["family"] = FAM.loc[RK.feature].values; RK.loc[RK.family == "Limb movements", "family"] = "Sleep timing and structure"   # the round-44c merge fig1_build.py applies
N = len(RK); RKI = RK.set_index("feature")
CALL = {"spo2_pct_below_90": "T90 (time below 90% saturation)", "AHI": "Apnea-hypopnea index", "N3_pct": "Deep sleep (N3)", "sleep_efficiency_pct": "Sleep efficiency", "TST_min": "Total sleep time"}
want = {f"{v}, rank {int(RKI.loc[k, 'rank'])}": (14.0, 1) for k, v in CALL.items()}
for f in sorted(set(RK.family)): want[f] = (14.0, 2)   # the key entry and the panel c row label
for s_ in (str(N), "1", "50", "100"): want[s_] = (14.0, 2)
want[f"Rank among the {N} parameters"] = (15.0, 2); want[f"Rank among the {N} parameters (rank 1 holds the largest gain)"] = (15.0, 1)
RAW = pd.read_csv(C.hydrated(f"{C.SV}/dC_per_disease/dC_per_disease.csv")); STD = RAW[(RAW.row_type == "disease") & (RAW.arm == "standard")].sort_values("dC", ascending=False)
ROWS = STD[STD.disease != "Obesity hypoventilation"].head(10)
def fmt_val(dC, lo, hi): return f"{dC:+.3f} ({lo:.3f} to {hi:.3f})".replace("-", "−")
for _, r in ROWS.iterrows(): want[str(r.disease)] = (14.0, want.get(str(r.disease), (14.0, 0))[1] + 1); want[fmt_val(float(r.dC), float(r.ci_lo), float(r.ci_hi))] = (14.0, 1)
want["Gain (95% CI)"] = (15.0, 1); want["Full follow-up"] = (14.0, 1); want["2-year landmark"] = (14.0, 1)
KM = json.load(open(C.hydrated(f"{paths.FIGURE_ROOT}/ROUND15_2026-09-04/L2_fig1d/L2_fig1d_values.json"))); tiles = KM["row1_existing"] + KM["row2_new"]
for pn in tiles:
    want[pn["label"]] = (14.0, want.get(pn["label"], (14.0, 0))[1] + 1)
    assert pn["p"] < 0.001; want["Log-rank P < 0.001"] = (13.0, 8)
    for v in pn["yticks"]: want[str(v)] = (14.0, want.get(str(v), (14.0, 0))[1] + 1)
for s_ in ("T90 ≤1%", "T90 1–10%", "T90 >10%", "Years since the sleep study", "Cumulative incidence, %"): want[s_] = (15.0 if s_[0] in "YC" else 14.0, 1)
bad = []
for s_, (sz, n_) in want.items():
    n_found = sum(len(re.findall(r"(?<!\S)" + re.escape(s_) + r"(?!\S)", l["text"])) for l in L_NEW if l["size"] == sz)
    if n_found < n_: bad.append((s_, sz, n_, n_found))
check(f"printed values re-derived from ranking_v3.csv, measure_families.csv, dC_per_disease.csv and L2_fig1d_values.json ({len(want)} strings) are on the sheet at their sizes", not bad, f"short {bad[:8]}")
DR = {k: json.load(open(f"{VER}/Main_Fig1_{k}_drawn.json")) for k in "bcde"}
recs = [(k, v) for k in "bcde" for v in DR[k]["values"] if v["kind"] in ("printed", "label", "tick")]
missing = [(k, v["text"]) for k, v in recs if not any(v["text"] == l["text"] or (" " + v["text"] + " ") in (" " + l["text"] + " ") for l in L_NEW)]
check(f"every value fig1_build.py recorded as printed, label or tick ({len(recs)}) is on the sheet", not missing, f"missing {missing[:8]}")
prov = []
for k in "bcde":
    for name, src in DR[k]["sources"].items(): prov.append((name, sha(src["path"]) == src["sha256"]))
check("numbers files unchanged since the build (sha256 of every recorded source)", all(p for _, p in prov), str([n for n, p in prov if not p]) if not all(p for _, p in prov) else ", ".join(sorted({n for n, _ in prov})))

# 4 no overlapping text boxes (every word on the sheet, integer-bbox tolerance 0.75 pt)
def depth(a, b): return min(a[2], b[2]) - max(a[0], b[0]), min(a[3], b[3]) - max(a[1], b[1])
hits = []
for i in range(len(W_NEW)):
    for j in range(i + 1, len(W_NEW)):
        dx_, dy_ = depth(W_NEW[i]["box"], W_NEW[j]["box"])
        if dx_ > 0.75 and dy_ > 0.75: hits.append((W_NEW[i]["text"], W_NEW[j]["text"], round(W_NEW[i]["y"], 1), round(dx_, 1), round(dy_, 1)))
check(f"no two text boxes overlap on the sheet ({len(W_NEW)} words, boxes cap height 0.727 and descender 0.21 of the size)", not hits, f"{hits[:8]}")

# 5 the key placement rule
g = DR["c"]["geometry"]; kc = g["key_check_r54"]; slack = g["slack_r54"]
check("key placement rule (build record): at least 6 pt from every callout box, 25 pt inside the axes' right edge, 6 pt above every bar beneath it, 4 pt under the axes top, with 2 pt of slack on every rule",
      kc["ok"] and min(slack.values()) >= 2.0, json.dumps(dict(gaps=kc["gaps_to_callouts"], key_right=kc["key_right"], key_bottom=kc["key_bottom"], bar_top=kc["highest_bar_top_under_key"], min_slack=min(slack.values()))))
key_lines = [l for l in L_NEW if l["size"] == 14.0 and abs(l["box"][0] - g["key_base0"] * 0 - 215.5) < 3.0 and l["text"] in set(RK.family)]
call_lines = [l for l in L_NEW if l["size"] == 14.0 and ", rank " in l["text"]]
kb = [min(l["box"][0] for l in key_lines), min(l["box"][1] for l in key_lines), max(l["box"][2] for l in key_lines), max(l["box"][3] for l in key_lines)]
gaps = {l["text"]: round(max(l["box"][0] - kb[2], kb[0] - l["box"][2], l["box"][1] - kb[3], kb[1] - l["box"][3]), 2) for l in call_lines}
check("key placement rule (census): the six key lines sit at least 6 pt from each of the five callout lines", len(key_lines) == 6 and len(call_lines) == 5 and min(gaps.values()) >= 6.0, f"gaps {gaps}")

# 6 marks preserved against the LFIG1E record (V30's panels)
RF = {k: json.load(open(f"{REF}/verify/Main_Fig1_{k}_drawn.json")) for k in "bcde"}
same_dots = sorted((d["rank"], d["family"]) for d in DR["b"]["dots"]) == sorted((d["rank"], d["family"]) for d in RF["b"]["dots"]) and DR["b"]["geometry"]["medians"] == RF["b"]["geometry"]["medians"]
same_bars = DR["c"]["bars"] == RF["c"]["bars"]
def marks(d, kind_key): return sorted((v["text"], json.dumps(v["value"], sort_keys=True)) for v in d["values"] if v["kind"] == "drawn" and kind_key in v["text"])
same_d = marks(DR["d"], "marker") == marks(RF["d"], "marker") and DR["d"]["geometry"]["rows"] == RF["d"]["geometry"]["rows"] and abs(DR["d"]["geometry"]["pub_dc_rule"] - RF["d"]["geometry"]["pub_dc_rule"]) < 1e-12
same_e = marks(DR["e"], "curve") == marks(RF["e"], "curve") and DR["e"]["geometry"]["tick_sets"] == RF["e"]["geometry"]["tick_sets"]
check("marks preserved: the 141 dots and medians, the 141 bars (rank, feature, family, gain, colour), the 19 forest markers with intervals, the 24 curves (8-year values, n, events) and the tick sets equal the LFIG1E (V30) record", same_dots and same_bars and same_d and same_e, f"dots {same_dots}, bars {same_bars}, forest {same_d}, tiles {same_e}")

# 7 letters and title at 18 pt, positions
bold = collections.Counter((t, w["size"]) for w in W_NEW if w["bold"] and w["size"] == 18.0 for t in w["text"].split()); bold30 = collections.Counter((t, w["size"]) for w in W_OLD if w["bold"] and w["size"] == 14.0 and w["region"] != "a" or (w["bold"] and w["text"] in ("Figure", "1", "a") and w["size"] == 14.0) for t in w["text"].split())
check("panel letters a to f and the title 'Figure 1' at 18 pt (V30 14), letters c, d, e at the re-laid panels' left edges", bold == collections.Counter({("a", 18.0): 1, ("b", 18.0): 1, ("c", 18.0): 1, ("d", 18.0): 1, ("e", 18.0): 1, ("f", 18.0): 1, ("Figure", 18.0): 1, ("1", 18.0): 1}) and [l["x"] for l in LOG["letters"]] == [14.17, 14.17, 483.6, 11.17, 473.4, 14.17],
      f"18-pt bold census {dict(bold)} (V30 14-pt bold {dict(bold30)}), letters {[(l['ch'], l['x'], l['baseline']) for l in LOG['letters']]}")
# 8 banned strings
BAD = [w["text"] for w in W_NEW if ("—" in w["text"] or ";" in w["text"])] + [t for t, b in C.banned_words([l["text"] for l in L_NEW])]
check("no em dash, semicolon or banned word in any text on the sheet", not BAD, str(BAD[:5]))
# 9 render
png = f"{SD}/Main_Fig1_150dpi.png"; check("Ghostscript 150-dpi render present", os.path.exists(png) and os.path.getsize(png) > 10000, f"{os.path.getsize(png)} bytes")
json.dump(dict(new_words=len(W_NEW), v30_words=len(W_OLD), panel_lines=[dict(text=l["text"], size=l["size"], y=round(l["y"], 2), box=[round(v, 2) for v in l["box"]]) for l in sorted(L_NEW, key=lambda l: (l["y"], l["box"][0]))],
               out_sha256=sha(NEW), v30_sha256=sha(V30)), open(f"{VER}/census_lines.json", "w"), indent=0)
lines_out.append("RESULT ALL PASS" if ok_all else "RESULT FAIL"); open(f"{VER}/checks.txt", "w").write("\n".join(lines_out) + "\n"); print(lines_out[-1])
