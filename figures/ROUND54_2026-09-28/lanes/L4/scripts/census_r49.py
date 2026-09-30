#!$T90_PY
"""ROUND 49 (2026-09-26, lane L4): the text census of the Main_Fig4 vector compose against the round-44 vector (the source of the V26
raster). Ghostscript txtwrite (-dTextFormat=0) on both sheets, never PyMuPDF on a composed sheet. Spans are parsed into (text, size,
bbox, font); spans on one baseline and size are joined in x order into LINES (a different span split of the same words between the two
sheets does not count). Rules: outside the band, the multiset of line strings is the same and every size is exactly +1.0, except the
declared exception regions, where the size is expected unchanged; the band (y >= band_top) is the coordinator's and is reported: its
word multiset must be the same, its sizes are listed as found.
usage: census_r49.py NEW.pdf OLD.pdf out_dir --band-top Y [--keep x0,y0,x1,y1 ...]   (a keep box = strings whose span origin lies inside stay at their size)"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import argparse, collections, html, json, os, re, subprocess, sys
GS = paths.GS
PLUS = 1.0


def txtwrite(pdf, xml):
    r = subprocess.run([GS, "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=txtwrite", "-dTextFormat=0", f"-sOutputFile={xml}", pdf], capture_output=True, text=True)
    assert r.returncode == 0 and os.path.exists(xml), ("gs txtwrite failed", pdf, r.stderr[-300:])
    return open(xml, encoding="utf-8", errors="replace").read()


def unesc(s): return html.unescape(s)   # named and numeric entities (txtwrite writes &#x2264; for the less-or-equal sign)


def spans(xml_text):
    """[(text, size, x0, baseline, x1, font, chars)] of every span with a non-blank text; chars = [(c, x0, x1)]."""
    out = []
    for m in re.finditer(r'<span bbox="([^"]*)" font="([^"]*)" size="([^"]*)">(.*?)</span>', xml_text, flags=re.S):
        bbox, font, size, body = m.groups(); x0, y0, x1, y1 = [float(v) for v in bbox.split()]
        chars = [(unesc(c), float(a), float(b)) for a, _, b, _, c in re.findall(r'<char bbox="([-\d.]+) ([-\d.]+) ([-\d.]+) ([-\d.]+)" c="([^"]*)"/>', body)]
        text = "".join(c for c, _, _ in chars)
        if not text: continue
        out.append(dict(text=text, size=round(float(size), 4), x0=x0, y=y0, x1=x1, font=font.split("+")[-1], chars=chars, blank=not text.strip()))   # blank spans kept: they carry the spaces of per-character text (the band)
    return out


def lines(sp, region=None):
    """Spans joined per (region, baseline, size) in x order: [(joined text, size, x0, x1, baseline, n_spans)]. region(span) -> key keeps
    strings of different panels on one baseline apart (Main_Fig4: the two forests share tick and sub-row baselines)."""
    rows = collections.defaultdict(list)
    for s in sp: rows[((region(s) if region else 0), round(s["y"]), s["size"])].append(s)
    out = []
    for (_, y, size), ss in rows.items():
        ss.sort(key=lambda s: s["x0"]); groups = [[ss[0]]]
        for s in ss[1:]:   # join only CONTIGUOUS spans (a kerning split has no gap; independent strings on one baseline are many points apart)
            if s["x0"] - groups[-1][-1]["x1"] <= 1.5: groups[-1].append(s)   # a kerning split is 0 to 1 pt in the integer boxes; the nearest real neighbours are 4 pt apart
            else: groups.append([s])
        for g in groups:
            t = "".join(s["text"] for s in g).strip()
            if t: out.append(dict(text=t, size=size, x0=g[0]["x0"], x1=g[-1]["x1"], y=y, n=len(g)))
    return out


def in_box(s, box): return box[0] <= s["x0"] <= box[2] and box[1] <= s["y"] <= box[3]


def words(ls): return collections.Counter(w for l in ls for w in l["text"].split())


def census(new_pdf, old_pdf, out_dir, band_top, keep_boxes=(), region=None):
    os.makedirs(out_dir, exist_ok=True)
    nsp = spans(txtwrite(new_pdf, f"{out_dir}/new_txtwrite.xml")); osp = spans(txtwrite(old_pdf, f"{out_dir}/old_txtwrite.xml"))
    nl, ol = lines(nsp, region), lines(osp, region)
    nb = [l for l in nl if l["y"] >= band_top]; ob = [l for l in ol if l["y"] >= band_top]
    nm = [l for l in nl if l["y"] < band_top]; om = [l for l in ol if l["y"] < band_top]
    # strings (line level) outside the band
    n_str = collections.Counter(l["text"] for l in nm); o_str = collections.Counter(l["text"] for l in om)
    lost = sorted((o_str - n_str).elements()); gained = sorted((n_str - o_str).elements())
    # sizes: every new line maps back to its expected old (text, size)
    kept = []; expected = collections.Counter()
    for l in nm:
        if any(in_box(l, b) for b in keep_boxes): expected[(l["text"], l["size"])] += 1; kept.append((l["text"], l["size"], l["x0"], l["y"]))
        else: expected[(l["text"], round(l["size"] - PLUS, 4))] += 1
    old_pairs = collections.Counter((l["text"], l["size"]) for l in om)
    miss = sorted((old_pairs - expected).elements()); extra = sorted((expected - old_pairs).elements())
    # the band: words the same, sizes listed
    band_words_same = words(nb) == words(ob)
    band_sizes = dict(old=sorted(collections.Counter(l["size"] for l in ob).items()), new=sorted(collections.Counter(l["size"] for l in nb).items()))
    rep = dict(new=new_pdf, old=old_pdf, band_top=band_top, keep_boxes=[list(b) for b in keep_boxes],
               n_spans=dict(new=sum(1 for x in nsp if not x["blank"]), old=sum(1 for x in osp if not x["blank"])), n_lines=dict(new=len(nl), old=len(ol)), n_lines_outside_band=dict(new=len(nm), old=len(om)),
               strings_same=(n_str == o_str), lost=lost, gained=gained, sizes_ok=(expected == old_pairs), size_missing_expected=miss, size_unexpected=extra,
               kept_at_size=kept, sizes_outside_band=dict(old=sorted(collections.Counter(l["size"] for l in om).items()), new=sorted(collections.Counter(l["size"] for l in nm).items())),
               band=dict(words_same=band_words_same, lines_old=[(l["text"], l["size"], l["x0"], l["y"]) for l in sorted(ob, key=lambda l: (l["y"], l["x0"]))],
                         lines_new=[(l["text"], l["size"], l["x0"], l["y"]) for l in sorted(nb, key=lambda l: (l["y"], l["x0"]))], sizes=band_sizes,
                         words_lost=sorted((words(ob) - words(nb)).elements()), words_gained=sorted((words(nb) - words(ob)).elements())),
               fonts=dict(old=sorted(set(s["font"] for s in osp)), new=sorted(set(s["font"] for s in nsp))))
    rep["ok"] = rep["strings_same"] and rep["sizes_ok"] and band_words_same
    json.dump(rep, open(f"{out_dir}/census.json", "w"), indent=1, ensure_ascii=False)
    return rep, nsp, osp


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("new"); ap.add_argument("old"); ap.add_argument("out"); ap.add_argument("--band-top", type=float, required=True); ap.add_argument("--keep", action="append", default=[])
    a = ap.parse_args(); boxes = [tuple(float(v) for v in k.split(",")) for k in a.keep]
    rep, _, _ = census(a.new, a.old, a.out, a.band_top, boxes)
    print(f"census: lines outside the band new {rep['n_lines_outside_band']['new']} old {rep['n_lines_outside_band']['old']}; strings same {rep['strings_same']} (lost {rep['lost'][:8]} gained {rep['gained'][:8]}); "
          f"sizes {rep['sizes_outside_band']['old']} -> {rep['sizes_outside_band']['new']}; every size +{PLUS} outside the keep boxes {rep['sizes_ok']} (kept {len(rep['kept_at_size'])} lines); band words same {rep['band']['words_same']}, band sizes {rep['band']['sizes']}")
    if not rep["sizes_ok"]: print("  missing expected", rep["size_missing_expected"][:20]); print("  unexpected", rep["size_unexpected"][:20])
    print("CENSUS PASS" if rep["ok"] else "CENSUS FAIL"); sys.exit(0 if rep["ok"] else 1)
