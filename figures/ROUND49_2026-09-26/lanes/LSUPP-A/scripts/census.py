#!/usr/bin/env python3
"""Round 49 text census (LSUPP-A): Ghostscript txtwrite (-dTextFormat=0) on the NEW sheet and on the V26 sheet, spans parsed into
(string, size), then compared: the same multiset of strings, every size exactly +1.0 (declared exceptions listed by string).
Usage: census.py NEW.pdf V26.pdf out_dir [--except "string=old_size" ...]. Prints a report and writes out_dir/census.json. Exit 0 on pass."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import os, re, sys, json, subprocess
from collections import Counter
import html
GS = paths.GS

def txtwrite(pdf, xml):
    subprocess.run([GS, "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=txtwrite", "-dTextFormat=0", f"-sOutputFile={xml}", pdf], check=True)
    return open(xml, encoding="utf-8", errors="replace").read()

def spans(xml_text):
    """Every <span size=...> as (string, size, bbox). Characters joined in order, leading and trailing blanks stripped."""
    out = []
    for m in re.finditer(r'<span bbox="([^"]*)" font="([^"]*)" size="([^"]*)">(.*?)</span>', xml_text, flags=re.S):
        bbox, font, size, body = m.groups()
        chars = re.findall(r'<char bbox="[^"]*" c="([^"]*)"/>', body)
        s = html.unescape("".join(chars)).strip()   # numeric entities too (Ghostscript writes U+2212 as &#x2212;)
        if not s: continue
        out.append((s, round(float(size), 2), font.split("+")[-1], bbox))
    return out

def lines(sp):
    """Spans joined per text line (same bbox top, in x order) so that a different span split between the two sheets does not count."""
    from collections import defaultdict
    rows = defaultdict(list)
    for s, size, font, bbox in sp:
        x0, y0, x1, y1 = [float(v) for v in bbox.split()]
        rows[(round(y0), size)].append((x0, s))
    return Counter(" ".join(t for _, t in sorted(v)) for v in rows.values())

def main():
    args = sys.argv[1:]; exc = {}
    while "--except" in args:
        i = args.index("--except"); k, v = args[i + 1].rsplit("=", 1); exc[k] = float(v); del args[i:i + 2]
    new_pdf, v26_pdf, out = args[:3]; os.makedirs(out, exist_ok=True)
    plus = 1.0
    nsp = spans(txtwrite(new_pdf, f"{out}/new_txtwrite.xml")); osp = spans(txtwrite(v26_pdf, f"{out}/v26_txtwrite.xml"))
    n_strings = Counter(s for s, _, _, _ in nsp); o_strings = Counter(s for s, _, _, _ in osp)
    same_strings = n_strings == o_strings
    lost = sorted((o_strings - n_strings).elements()); gained = sorted((n_strings - o_strings).elements())
    # per (string, size): every old (s, size) must reappear as (s, size + 1), except the declared exceptions (kept at their size)
    o_pairs = Counter((s, size) for s, size, _, _ in osp); n_pairs = Counter((s, size) for s, size, _, _ in nsp)
    expected = Counter()
    for (s, size), k in o_pairs.items():
        expected[(s, size if (s in exc and abs(exc[s] - size) < 0.01) else round(size + plus, 2))] += k
    size_ok = expected == n_pairs
    miss = sorted((expected - n_pairs).elements()); extra = sorted((n_pairs - expected).elements())
    fonts_o = Counter(f for _, _, f, _ in osp); fonts_n = Counter(f for _, _, f, _ in nsp)
    # fall back to line-level strings when span splitting differs
    line_same = lines(nsp) and (Counter(k for k in lines(nsp)) == Counter(k for k in lines(osp)))
    sizes_o = sorted({size for _, size, _, _ in osp}); sizes_n = sorted({size for _, size, _, _ in nsp})
    rep = dict(new=new_pdf, v26=v26_pdf, n_spans_new=len(nsp), n_spans_v26=len(osp), same_strings=same_strings, lost=lost, gained=gained,
               sizes_v26=sizes_o, sizes_new=sizes_n, every_size_plus_one=size_ok, missing_expected=miss, unexpected=extra, exceptions=exc,
               fonts_v26=dict(fonts_o), fonts_new=dict(fonts_n), line_level_same=bool(line_same))
    json.dump(rep, open(f"{out}/census.json", "w"), indent=1)
    print(f"census: spans new {len(nsp)} v26 {len(osp)}; strings same {same_strings} (lost {lost[:10]} gained {gained[:10]}); sizes v26 {sizes_o} -> new {sizes_n}; every size +1.0 {size_ok} exceptions {exc}")
    if not size_ok: print("  missing expected", miss[:20]); print("  unexpected", extra[:20])
    ok = same_strings and size_ok
    print("CENSUS PASS" if ok else "CENSUS FAIL"); sys.exit(0 if ok else 1)

if __name__ == "__main__": main()
