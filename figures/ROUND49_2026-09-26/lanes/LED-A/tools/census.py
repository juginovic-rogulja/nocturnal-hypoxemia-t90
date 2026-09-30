#!$T90_PY
"""Round 49, lane LED-A: Ghostscript text census of a vector sheet and the +1 pt comparison against the V26 sheet.
census(pdf, xml): gs txtwrite TextFormat=0 -> one item per <span> (text joined from its <char> elements, size, font, bbox).
compare(old_items, new_items, delta=1.0, exceptions={text: expected_new_size}): same multiset of strings, every size + delta.
usage: census.py <pdf> <xml_out> [<json_out>]     |     census.py --compare <old.json> <new.json> [delta]"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import html, json, os, re, subprocess, sys
from collections import Counter

GS = paths.GS


def census(pdf, xml_out):
    assert os.stat(pdf).st_blocks > 0, ("evicted", pdf)
    subprocess.run([GS, "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=txtwrite", "-dTextFormat=0", f"-sOutputFile={xml_out}", pdf], check=True, timeout=600)
    txt = open(xml_out, encoding="utf-8", errors="replace").read()
    items = []
    for m in re.finditer(r'<span bbox="([^"]*)" font="([^"]*)" size="([^"]*)">(.*?)</span>', txt, re.S):
        bbox = [float(v) for v in m.group(1).split()]
        chars = re.findall(r'<char bbox="([^"]*)" c="([^"]*)"/>', m.group(4))
        s = "".join(html.unescape(c) for _, c in chars)
        if not s.strip():
            continue
        items.append(dict(text=s, size=round(float(m.group(3)), 3), font=m.group(2), bbox=bbox, x=bbox[0], y=bbox[1]))
    return items


def compare(old, new, delta=1.0, exceptions=None, tol=0.011, removed=(), added=()):
    """Returns (ok, lines). Strings: same multiset (after nbsp -> space, strip) apart from the DECLARED removed/added strings (a name that
    the wrap rule splits differently at the larger size). Sizes: per string, sorted old sizes + delta == sorted new sizes (exceptions: {text: size})."""
    exceptions = exceptions or {}
    def nz(t): return t.replace("\xa0", " ").strip()
    co, cn = Counter(nz(i["text"]) for i in old), Counter(nz(i["text"]) for i in new)
    lines = []; ok = True
    missing, extra = co - cn, cn - co
    dec_rem, dec_add = Counter(nz(t) for t in removed), Counter(nz(t) for t in added)
    dec_rem, dec_add = dec_rem - dec_add, dec_add - dec_rem          # a line common to both wraps is neither removed nor added
    if missing == dec_rem and extra == dec_add:
        if dec_rem or dec_add: lines.append(f"PASS strings: the same multiset of strings apart from the declared re-wrapped lines ({sum(co.values())} on V26, {sum(cn.values())} new): removed {dict(dec_rem)}, added {dict(dec_add)}")
        else: lines.append(f"PASS strings: the same multiset of {sum(co.values())} strings ({len(co)} distinct)")
    else:
        ok = False; lines.append(f"FAIL strings: undeclared difference: missing on the new sheet {dict(missing - dec_rem)}, extra on the new sheet {dict(extra - dec_add)}, declared but not seen: removed {dict(dec_rem - missing)} added {dict(dec_add - extra)}")
    so, sn = {}, {}
    for i in old: so.setdefault(nz(i["text"]), []).append(i["size"])
    for i in new: sn.setdefault(nz(i["text"]), []).append(i["size"])
    bad = []; n_ok = 0; exc_used = []
    for t in sorted(so):
        if t not in sn: continue
        a, b = sorted(so[t]), sorted(sn[t])
        if len(a) != len(b): bad.append((t, a, b)); continue
        for x, y in zip(a, b):
            exp = exceptions.get(t, x + delta)
            if t in exceptions: exc_used.append((t, x, y))
            if abs(y - exp) > tol: bad.append((t, x, y, exp))
            else: n_ok += 1
    if bad:
        ok = False; lines.append(f"FAIL sizes: {len(bad)} string sizes off the +{delta} rule: {bad[:12]}")
    else:
        lines.append(f"PASS sizes: every one of the {n_ok} strings is exactly +{delta} pt against V26" + (f" (declared exceptions: {exc_used})" if exc_used else " (no exceptions)"))
    sizes_old = sorted(set(i["size"] for i in old)); sizes_new = sorted(set(i["size"] for i in new))
    lines.append(f"INFO sizes present: V26 {sizes_old} -> new {sizes_new}")
    return ok, lines


if __name__ == "__main__":
    if sys.argv[1] == "--compare":
        old = json.load(open(sys.argv[2])); new = json.load(open(sys.argv[3])); d = float(sys.argv[4]) if len(sys.argv) > 4 else 1.0
        ok, lines = compare(old, new, d); print("\n".join(lines)); sys.exit(0 if ok else 1)
    items = census(sys.argv[1], sys.argv[2])
    if len(sys.argv) > 3: json.dump(items, open(sys.argv[3], "w"), indent=0)
    c = Counter(i["size"] for i in items)
    print(f"{sys.argv[1]}: {len(items)} spans; sizes {dict(sorted(c.items()))}")
