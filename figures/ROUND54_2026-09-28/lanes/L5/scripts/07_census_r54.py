#!$T90_PY
"""Round 54, lane L5: Ghostscript text census (txtwrite -dTextFormat=0) of the OLD sheet and of the NEW sheet, compared as multisets of
strings and tokens, every token's size mapped old -> new by the round-54 size map (exceptions listed).
OLD = V30's Main_Fig5 (ROUND52 NEW_FINAL_SET_V30, sha256 f41fec63...) which the assembly manifests (V27 = V28 = V29 = V30) trace to the
ROUND49 LRASTER 600-dpi raster of the round-49 L5 vector compose (Main_Fig5_r49_vector.pdf, sha256 622925e0...): its text is read
  (a) from the round-49 lane's own census of that vector (census_new_vector.xml, Ghostscript, duplicate spans collapsed as in round 49) and
  (b) from a fresh Ghostscript txtwrite of a pdfwrite twin of the round-49 vector (work/old_r49_flat_twin.pdf), the two agreeing.
NEW = the round-54 vector compose, censused through the lane's flat twin work/Main_Fig5_flat.pdf (gs pdfwrite of the vector, the
coordinator's recipe), cross-checked against the pieces (the round-54 flat page, the round-54 panel a page with sizes x S, the band, the
stamped letters and title) and, where it finishes under the watchdog, against a txtwrite of the vector itself.
Strings: spans grouped by (font, size, baseline) and joined along x (a gap of at least 0.2 em = a space, of at least 2.5 em = a new string);
tokens = whitespace-split strings. Size map old -> new: 10.02 -> 14 (HR, P, q values) or 15 (the column heads, which V13 drew at the value size), 10.5 -> 14 (tick labels, key, panel c labels),
11.45 -> 14 (panel b labels) or 15 (column heads, axis titles, 'Negative controls'), 12.34 -> 14 (panel a blocks and arms), 13.48 -> 15
(panel a captions), 14 -> 18 (letters, title), band 10.99 -> 14 (caption), 11.99 -> 15 (notes). Writes work/census_r54.json.
usage: 07_census_r54.py"""
import collections, html, json, os, re, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import l5_lib as L
SHEET = "Main_Fig5"; SD = f"{L.LANE}/{SHEET}"; WORK = f"{SD}/work"; clog = json.load(open(f"{WORK}/compose_log.json")); S_A = clog["panel_a"]["map"]["S"]
R49W = f"{L.L5_R49}/Main_Fig5/work"
# old -> allowed new sizes. 10.02 was the size of the values AND of the column heads (V13 drew the heads at the value size): values -> 14, heads -> 15.
# The band is placed at its own page scale 0.99975: its 14-pt caption censuses at 13.9965 and its 15-pt notes at 14.996 (matched within 0.03).
SIZE_MAP = {10.02: {14.0, 15.0}, 10.5: {14.0}, 10.99: {14.0}, 11.45: {14.0, 15.0}, 11.99: {15.0}, 12.34: {14.0}, 13.48: {15.0}, 14.0: {18.0}}
TOL = 0.03


def gs_spans(pdf, xml, name, max_s=600):
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


def dedup(spans):
    """Identical spans (text, bbox, font, size) counted once: a flat page placed twice (two clips) is emitted twice by txtwrite."""
    seen = collections.Counter(); out = []
    for s in spans:
        k = (s["text"], s["x0"], s["y"], s["x1"], s["font"], round(s["size"], 3)); seen[k] += 1
        if seen[k] == 1: out.append(s)
    return out, collections.Counter(seen.values())


def strings_of(spans, scale=1.0):
    """Join spans on a line (same font, size, baseline) along x; returns [(string, size, font, x0, y, x1)]."""
    groups = collections.defaultdict(list)
    for s in spans: groups[(s["font"], round(s["size"], 2), round(s["y"], 0))].append(s)
    out = []
    for (font, size, y), ss in groups.items():
        ss.sort(key=lambda s: s["x0"]); cur = ss[0]["text"]; x0 = ss[0]["x0"]; last = ss[0]["x1"]
        for s in ss[1:]:
            gap = s["x0"] - last
            if gap >= 2.5 * size: out.append((cur, round(size * scale, 2), font, x0, y, last)); cur = s["text"]; x0 = s["x0"]
            elif gap >= 0.2 * size and not cur.endswith(" ") and not s["text"].startswith(" "): cur += " " + s["text"]
            else: cur += s["text"]
            last = max(last, s["x1"])
        out.append((cur, round(size * scale, 2), font, x0, y, last))
    return [(re.sub(r"\s+", " ", t).strip(), sz, f, x, y, x1) for t, sz, f, x, y, x1 in out if t.strip()]


def words(strings): return collections.Counter(tok for t, *_ in strings for tok in t.split())
def texts(strings): return collections.Counter(t for t, *_ in strings)
def tok_sizes(strings):
    d = collections.defaultdict(list)
    for t, sz, *_ in strings:
        for tok in t.split(): d[tok].append(sz)
    return {k: sorted(v) for k, v in d.items()}


def size_match(old, new):
    """Per token: the old sizes and the new sizes must pair through SIZE_MAP (a bijection; 11.45 may go to 14 or 15). Returns
    (n tokens paired, exceptions [(token, old sizes, new sizes, why)])."""
    ok = 0; bad = []
    for tok in sorted(set(old) | set(new)):
        so, sn = old.get(tok, []), new.get(tok, [])
        if len(so) != len(sn): bad.append((tok, so, sn, "count")); continue
        rem = list(sn); good = True
        single = [s for s in so if len(SIZE_MAP.get(s, set())) == 1]; multi = [s for s in so if len(SIZE_MAP.get(s, set())) > 1]; unknown = [s for s in so if s not in SIZE_MAP]
        if unknown: bad.append((tok, so, sn, f"old size not in the map {unknown}")); continue
        def take(targets):
            for t in sorted(targets):
                for i, v in enumerate(rem):
                    if abs(v - t) <= TOL: rem.pop(i); return True
            return False
        for s in single:
            if not take(SIZE_MAP[s]): good = False
        for s in multi:
            if not take(SIZE_MAP[s]): good = False
        if good and not rem: ok += len(so)
        else: bad.append((tok, so, sn, "sizes do not pair through the map"))
    return ok, bad


t0 = time.time()
# ---- OLD: the round-49 vector's own census (the round-49 lane), collapsed as in round 49, and the fresh twin census
old_raw = parse(f"{R49W}/census_new_vector.xml"); old_spans, old_dup = dedup(old_raw); OLD = strings_of(old_spans)
old_twin_pdf = f"{WORK}/old_r49_flat_twin.pdf"
if not os.path.exists(old_twin_pdf): L.gs_flatten(L.R49_VECTOR, old_twin_pdf, "flatten_r49_old")
old_twin_raw = gs_spans(old_twin_pdf, f"{WORK}/census_old_twin.xml", "txt_old_twin"); old_twin_spans, old_twin_dup = dedup(old_twin_raw); OLD_TWIN = strings_of(old_twin_spans)
# ---- NEW: the round-54 flat twin (the coordinator's census route), the pieces, the vector itself where it finishes
new_twin_pdf = f"{WORK}/Main_Fig5_flat.pdf"
if not (os.path.exists(new_twin_pdf) and os.path.getmtime(new_twin_pdf) > os.path.getmtime(clog["output"])): L.gs_flatten(clog["output"], new_twin_pdf, "flatten_r54")
new_twin_raw = gs_spans(new_twin_pdf, f"{WORK}/census_new_twin.xml", "txt_new_twin"); new_twin_spans, new_twin_dup = dedup(new_twin_raw); NEW = strings_of(new_twin_spans)
new_a = strings_of(gs_spans(clog["panel_a"]["page"], f"{WORK}/census_new_panel_a.xml", "census_new_a"), scale=S_A)
new_flat = strings_of(gs_spans(clog["flat"], f"{WORK}/census_new_flat.xml", "census_new_flat"))
new_band = strings_of(gs_spans(clog["band"], f"{WORK}/census_new_band.xml", "census_new_band"))
new_stamp = [("a", 18.0, "Arial-BoldMT", 36.0, 41.76, 48.0), ("d", 18.0, "Arial-BoldMT", 36.0, clog["letter_d"]["baseline"], 48.0), ("Figure 5", 18.0, "Arial-BoldMT", 36.0, 18.0, 106.0)]
NEW_PIECES = new_a + new_flat + new_band + new_stamp
vec_xml = f"{WORK}/census_new_vector.xml"; vec_ok = False; NEW_VEC = None
try:
    rc, st = L.wd("txt_new_vector", [L.GS, "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=txtwrite", "-dTextFormat=0", f"-sOutputFile={vec_xml}", clog["output"]], max_s=300)
    vec_ok = rc == 0 and os.path.exists(vec_xml)
    if vec_ok: NEW_VEC = strings_of(dedup(parse(vec_xml))[0])
except Exception as e: st = f"{type(e).__name__}: {e}"

# ---- compare
res = {}
wo, wot, wn, wnp = words(OLD), words(OLD_TWIN), words(NEW), words(NEW_PIECES)
res["old_sources"] = {"v30": f"{L.V30}/Main_Fig5.pdf", "v30_sha256": L.sha256(f"{L.V30}/Main_Fig5.pdf"), "chain": "V30 = V29 = V28 = V27 Main_Fig5.pdf (sha256 f41fec63...) = ROUND49 LRASTER 600-dpi raster of the round-49 L5 vector", "r49_vector": L.R49_VECTOR, "r49_vector_sha256": L.sha256(L.R49_VECTOR), "r49_census_xml": f"{R49W}/census_new_vector.xml", "old_twin": old_twin_pdf, "old_twin_sha256": L.sha256(old_twin_pdf)}
res["new_sources"] = {"vector": clog["output"], "vector_sha256": clog["output_sha256"], "twin": new_twin_pdf, "twin_sha256": L.sha256(new_twin_pdf), "band": clog["band"], "band_sha256": clog["band_sha256"]}
res["old_r49_census_vs_old_twin"] = dict(equal=wo == wot, lost=dict(wo - wot), gained=dict(wot - wo), sizes_equal=tok_sizes(OLD) == tok_sizes(OLD_TWIN))
res["strings_old_vs_new_twin"] = dict(equal=wo == wn, lost=dict(wo - wn), gained=dict(wn - wo))
res["strings_old_vs_new_twin_joined"] = dict(equal=texts(OLD) == texts(NEW), lost=dict(texts(OLD) - texts(NEW)), gained=dict(texts(NEW) - texts(OLD)))
res["strings_new_pieces_vs_new_twin"] = dict(equal=wnp == wn, lost=dict(wnp - wn), gained=dict(wn - wnp))
to, tn, tnp = tok_sizes(OLD), tok_sizes(NEW), tok_sizes(NEW_PIECES)
n_ok, bad = size_match(to, tn); res["sizes_old_vs_new_twin"] = dict(n_tokens_mapped=n_ok, exceptions=bad, size_map={str(k): sorted(v) for k, v in SIZE_MAP.items()})
n_ok_p, bad_p = size_match(to, tnp); res["sizes_old_vs_new_pieces"] = dict(n_tokens_mapped=n_ok_p, exceptions=bad_p)
res["new_vector_census"] = dict(finished=vec_ok, status=st)
if vec_ok:
    wv = words(NEW_VEC); res["new_vector_census"].update(equal_to_twin=wv == wn, lost=dict(wn - wv), gained=dict(wv - wn), sizes_equal=tok_sizes(NEW_VEC) == tn, n_strings=len(NEW_VEC))
res["span_multiplicity"] = dict(old_r49=dict(old_dup), old_twin=dict(old_twin_dup), new_twin=dict(new_twin_dup))
res["min_size_new"] = min(sz for _, sz, *_ in NEW); res["min_size_old"] = min(sz for _, sz, *_ in OLD)
res["counts"] = dict(old=len(OLD), old_twin=len(OLD_TWIN), new_twin=len(NEW), new_pieces=dict(a=len(new_a), flat=len(new_flat), band=len(new_band), stamp=len(new_stamp)), new_vector=len(NEW_VEC) if NEW_VEC else None)
res["size_histogram"] = dict(old=sorted(collections.Counter(sz for _, sz, *_ in OLD).items()), new_twin=sorted(collections.Counter(sz for _, sz, *_ in NEW).items()), new_pieces=sorted(collections.Counter(sz for _, sz, *_ in NEW_PIECES).items()))
# per string: old text, old size, new size (matched by text within the same font face where possible)
def by_text(strings):
    d = collections.defaultdict(list)
    for t, sz, f, x, y, x1 in strings: d[t].append((sz, f, x, y))
    return d
bo, bn = by_text(OLD), by_text(NEW); table = []
for t in sorted(set(bo) | set(bn), key=lambda t: (min((y for _, _, _, y in bo.get(t, bn.get(t))), default=0), t)):
    table.append(dict(text=t, old=sorted((sz, round(x), round(y)) for sz, _, x, y in bo.get(t, [])), new=sorted((sz, round(x), round(y)) for sz, _, x, y in bn.get(t, []))))
res["per_string"] = table
res["strings"] = dict(old=sorted(OLD, key=lambda s: (s[4], s[3])), new_twin=sorted(NEW, key=lambda s: (s[4], s[3])))
res["seconds"] = round(time.time() - t0, 1); res["written"] = L.now()
json.dump(res, open(f"{WORK}/census_r54.json", "w"), indent=1, default=str)
print(json.dumps({k: v for k, v in res.items() if k not in ("strings", "per_string")}, indent=1, default=str)[:7000])
