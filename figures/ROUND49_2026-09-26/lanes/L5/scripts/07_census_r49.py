#!$T90_PY
"""Round 49, lane L5: Ghostscript text census (txtwrite -dTextFormat=0) of the OLD sheet (V26 = the round-44 vector compose: its pieces, and
the whole vector where its census finished) and of the NEW compose, compared as multisets of strings with sizes (+1.0 rule, exceptions listed).
Strings: spans grouped by (font, size, baseline) and joined along x (a gap of at least 0.2 em = a space, a gap of at least 2.5 em = a new
string); tokens = whitespace-split strings (robust to column joins). OLD pieces: the V13 panel a text probe (work/base_text.json, the round-37
PyMuPDF read of the nested V13 sheet, x < 497 and 16 < y < 393, the letter a included), the round-42 flat page (gs), the round-44 band (gs),
the stamped 'd' 13 and 'Figure 5' 13. NEW pieces: the round-49 flow page (gs, sizes x S), the round-49 flat page (gs), the round-49 band (gs),
'a' 14, 'd' 14, 'Figure 5' 14; and the composed vector itself (gs). Writes work/census_r49.json. usage: 07_census_r49.py"""
import collections, html, json, os, re, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import l5_lib as L
SHEET = "Main_Fig5"; SD = f"{L.LANE}/{SHEET}"; WORK = f"{SD}/work"; clog = json.load(open(f"{WORK}/compose_log.json")); S_A = clog["panel_a"]["map"]["S"]


def gs_spans(pdf, xml, name, max_s=1800):
    if not (os.path.exists(xml) and os.path.getmtime(xml) > os.path.getmtime(pdf)):
        rc, st = L.wd(name, [L.GS, "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=txtwrite", "-dTextFormat=0", f"-sOutputFile={xml}", pdf], max_s=max_s)
        assert rc == 0 and os.path.exists(xml), ("gs txtwrite failed", name, st)
    return parse(xml)


def parse(xml_path):
    xml = open(xml_path, encoding="utf-8", errors="replace").read(); out = []
    for m in re.finditer(r'<span bbox="([^"]*)" font="([^"]*)" size="([^"]*)">(.*?)</span>', xml, re.S):
        bb = [float(v) for v in m.group(1).split()]; chars = re.findall(r'c="([^"]*)"', m.group(4)); txt = L.norm_text(html.unescape("".join(chars)))
        if not txt.strip(): continue
        out.append(dict(text=txt, x0=bb[0], y=bb[1], x1=bb[2], font=m.group(2).split("+")[-1].replace("-Identity-H", ""), size=float(m.group(3))))
    return out


def strings_of(spans, scale=1.0):
    """Join spans on a line (same font, size, baseline) along x; returns [(string, size, font, x0, y)]."""
    groups = collections.defaultdict(list)
    for s in spans: groups[(s["font"], round(s["size"], 2), round(s["y"], 0))].append(s)
    out = []
    for (font, size, y), ss in groups.items():
        ss.sort(key=lambda s: s["x0"]); cur = ss[0]["text"]; x0 = ss[0]["x0"]; last = ss[0]["x1"]
        for s in ss[1:]:
            gap = s["x0"] - last
            if gap >= 2.5 * size: out.append((cur, round(size * scale, 2), font, x0, y)); cur = s["text"]; x0 = s["x0"]
            elif gap >= 0.2 * size and not cur.endswith(" ") and not s["text"].startswith(" "): cur += " " + s["text"]
            else: cur += s["text"]
            last = max(last, s["x1"])
        out.append((cur, round(size * scale, 2), font, x0, y))
    return [(re.sub(r"\s+", " ", t).strip(), sz, f, x, y) for t, sz, f, x, y in out if t.strip()]


def tokens(strings): return collections.Counter((tok, sz) for t, sz, *_ in strings for tok in t.split())
def tok_sizes(strings):
    d = collections.defaultdict(list)
    for t, sz, *_ in strings:
        for tok in t.split(): d[tok].append(sz)
    return {k: sorted(v) for k, v in d.items()}


t0 = time.time(); J = json.load(open(f"{WORK}/base_text.json"))
old_a = [(re.sub(r"\s+", " ", L.norm_text(s["text"])).strip(), round(s["size"], 2), s["font"], s["x"], s["y"]) for s in J["spans"] if s["x"] < 497 and 16 < s["y"] < 393]
old_flat = strings_of(gs_spans(f"{L.LFIG5B}/{SHEET}/work/panels_bc_flat.pdf", f"{WORK}/census_old_flat_r42.xml", "census_old_flat"))
old_band = strings_of(gs_spans(f"{L.SKETCH_R44}/band_5d_r44.pdf", f"{WORK}/census_old_band_r44.xml", "census_old_band"))
old_stamp = [("d", 13.0, "Arial-BoldMT", 36.0, 744.76), ("Figure 5", 13.0, "Arial-BoldMT", 36.0, 14.0)]
OLD = old_a + old_flat + old_band + old_stamp
new_a = strings_of(gs_spans(clog["panel_a"]["page"], f"{WORK}/census_new_panel_a.xml", "census_new_a"), scale=S_A)
new_flat = strings_of(gs_spans(clog["flat"], f"{WORK}/census_new_flat.xml", "census_new_flat"))
new_band = strings_of(gs_spans(clog["band"], f"{WORK}/census_new_band.xml", "census_new_band"))
new_stamp = [("a", 14.0, "Arial-BoldMT", 36.0, 41.76), ("d", 14.0, "Arial-BoldMT", 36.0, 744.76), ("Figure 5", 14.0, "Arial-BoldMT", 36.0, 14.0)]
NEW_PIECES = new_a + new_flat + new_band + new_stamp
comp_raw = gs_spans(clog["output"], f"{WORK}/census_new_vector.xml", "census_new_vector", max_s=3600)
# the flat page is placed twice (clips b and c) and txtwrite emits the clipped-away copy too: identical spans (text, bbox, font, size) are counted once
seen = collections.Counter(); comp_spans = []
for s in comp_raw:
    k = (s["text"], s["x0"], s["y"], s["x1"], s["font"], round(s["size"], 3)); seen[k] += 1
    if seen[k] == 1: comp_spans.append(s)
DUP = collections.Counter(seen.values()); NEW_COMP = strings_of(comp_spans)
# the old whole-vector census, if the background run finished
old_vec_xml = f"{WORK}/census_r44_vector.xml"; st_file = f"{L.LANE}/logs/watchdog_census/census_r44_vector.status"
old_vec_done = os.path.exists(st_file) and "rc=0" in open(st_file).read() and os.path.exists(old_vec_xml)
OLD_COMP = strings_of(parse(old_vec_xml)) if old_vec_done else None

# ---- compare: token multisets ignoring size, then sizes per token
def words(strings): return collections.Counter(tok for t, *_ in strings for tok in t.split())
res = {}
wo, wn, wnc = words(OLD), words(NEW_PIECES), words(NEW_COMP)
res["strings_old_pieces_vs_new_pieces"] = dict(equal=wo == wn, lost=dict(wo - wn), gained=dict(wn - wo))
res["strings_new_pieces_vs_new_composed"] = dict(equal=wn == wnc, lost=dict(wn - wnc), gained=dict(wnc - wn))
# sizes: for every token, the sorted old sizes and the sorted new sizes must pair as +1.0 (tolerance 0.02), exceptions listed
to, tn, tnc = tok_sizes(OLD), tok_sizes(NEW_PIECES), tok_sizes(NEW_COMP)
def size_pairs(a, b):
    ok, bad = [], []
    for tok in sorted(set(a) | set(b)):
        sa, sb = a.get(tok, []), b.get(tok, [])
        if len(sa) != len(sb): bad.append((tok, sa, sb, "count")); continue
        for x, y in zip(sa, sb):
            (ok if abs(y - x - 1.0) <= 0.02 else bad).append((tok, x, y, round(y - x, 3)))
    return ok, bad
ok_p, bad_p = size_pairs(to, tn); ok_c, bad_c = size_pairs(to, tnc)
res["sizes_old_pieces_vs_new_pieces"] = dict(n_plus_one=len(ok_p), exceptions=bad_p)
res["sizes_old_pieces_vs_new_composed"] = dict(n_plus_one=len(ok_c), exceptions=bad_c)
res["old_vector_census"] = dict(finished=old_vec_done)
if old_vec_done:
    woc = words(OLD_COMP); res["old_vector_census"].update(strings_equal_to_old_pieces=woc == wo, lost=dict(wo - woc), gained=dict(woc - wo), n_strings=len(OLD_COMP))
    okv, badv = size_pairs(tok_sizes(OLD_COMP), tnc); res["sizes_old_composed_vs_new_composed"] = dict(n_plus_one=len(okv), exceptions=badv)
res["composed_span_multiplicity"] = dict(DUP); res["counts"] = dict(old_pieces=dict(a=len(old_a), flat=len(old_flat), band=len(old_band), stamp=len(old_stamp)), new_pieces=dict(a=len(new_a), flat=len(new_flat), band=len(new_band), stamp=len(new_stamp)), new_composed=len(NEW_COMP), old_composed=len(OLD_COMP) if OLD_COMP else None)
res["size_histogram"] = dict(old=sorted(collections.Counter(sz for _, sz, *_ in OLD).items()), new_pieces=sorted(collections.Counter(sz for _, sz, *_ in NEW_PIECES).items()), new_composed=sorted(collections.Counter(sz for _, sz, *_ in NEW_COMP).items()))
res["strings"] = dict(old=sorted(OLD, key=lambda s: (s[4], s[3])), new_composed=sorted(NEW_COMP, key=lambda s: (s[4], s[3])))
res["seconds"] = round(time.time() - t0, 1); res["written"] = L.now()
json.dump(res, open(f"{WORK}/census_r49.json", "w"), indent=1, default=str)
print(json.dumps({k: v for k, v in res.items() if k != "strings"}, indent=1, default=str)[:6000])
