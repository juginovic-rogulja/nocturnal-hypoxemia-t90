#!$T90_PY
import fitz, json, sys
from collections import Counter
NEW, DUMP = sys.argv[1:3]
dump = json.load(open(DUMP)); od = dump["drawings"]
d = fitz.open(NEW); p = d[0]; nd = []
for it in p.get_drawings():
    r = it["rect"]; kinds = "".join(sorted(set(k[0] for k in it["items"])))
    nd.append(dict(rect=[r.x0, r.y0, r.x1, r.y1], fill=it.get("fill"), color=it.get("color"), width=it.get("width"), dashes=it.get("dashes"), kinds=kinds, n=len(it["items"])))
def key(it): return (tuple(round(v, 2) for v in it["rect"]), tuple(round(v, 4) for v in it["fill"]) if it["fill"] else None, tuple(round(v, 4) for v in it["color"]) if it["color"] else None, round(it["width"] or 0, 3), str(it.get("dashes")), it["kinds"], it["n"])
ko, kn = Counter(key(it) for it in od), Counter(key(it) for it in nd)
lost, gained = list((ko - kn).keys()), list((kn - ko).keys())
print("n old", len(od), "n new", len(nd))
for a in lost:
    b = min(gained, key=lambda g: sum(abs(x - y) for x, y in zip(a[0], g[0])))
    print("LOST  ", a); print("GAINED", b); print()
