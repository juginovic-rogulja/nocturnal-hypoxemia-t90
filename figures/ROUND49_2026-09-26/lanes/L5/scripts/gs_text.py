#!$T90_PY
"""Text layer of a composed (non-nested) sheet through Ghostscript txtwrite -dTextFormat=0 (spans with integer bbox x0, baseline, x1,
font, size) under the lane watchdog. Never used on the V13 base sheets or on Main_Fig5 (nested V13 panel a)."""
import html, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import l5_lib as L


def spans_xml(pdf, xml, name, max_s=600):
    L.hydrated(pdf)
    rc, st = L.wd(name, [L.GS, "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=txtwrite", "-dTextFormat=0", f"-sOutputFile={xml}", pdf], max_s=max_s)
    assert rc == 0 and os.path.exists(xml), ("gs txtwrite failed", name, st)
    return parse(xml)


def parse(xml_path):
    xml = open(xml_path, encoding="utf-8", errors="replace").read(); out = []
    for m in re.finditer(r'<span bbox="([^"]*)" font="([^"]*)" size="([^"]*)">(.*?)</span>', xml, re.S):
        bb = [float(v) for v in m.group(1).split()]; chars = re.findall(r'c="([^"]*)"', m.group(4))
        txt = L.norm_text(html.unescape("".join(chars)))
        if not txt.strip(): continue
        out.append(dict(text=txt, x0=bb[0], y=bb[1], x1=bb[2], font=m.group(2).split("+")[-1].replace("-Identity-H", ""), size=round(float(m.group(3)), 2)))
    return out


def key(s): return (s["text"].strip(), s["x0"], s["y"], s["x1"], s["font"], s["size"])


def letters(spans):
    return sorted([(s["text"].strip(), s["x0"], s["y"], s["font"], s["size"]) for s in spans if len(s["text"].strip()) == 1 and s["text"].strip() in "abcdefgh" and s["size"] >= 12.5 and "Bold" in s["font"]], key=lambda t: (t[2], t[1]))


def title(spans, n):
    return [s for s in spans if s["text"].strip() == f"Figure {n}"]


def strings(spans): return [s["text"].strip() for s in spans]
