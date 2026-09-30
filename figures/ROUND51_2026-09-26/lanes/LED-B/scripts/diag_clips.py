#!$T90_PY
"""List the V28 ED_Fig06 drawing items that are clipped out (rect outside their effective scissor) via get_drawings(extended=True)."""
import fitz, sys, json
d = fitz.open(sys.argv[1]); p = d[0]
ext = p.get_drawings(extended=True); plain = p.get_drawings()
print("plain items", len(plain), "extended entries", len(ext), "clips", sum(1 for e in ext if e.get("type") == "clip"), "groups", sum(1 for e in ext if e.get("type") == "group"))
# effective scissor per path: walk in order, keep a stack by level
stack = {}; out = []; k = 0
for e in ext:
    lvl = e.get("level", 0)
    if e.get("type") == "clip":
        stack[lvl] = e.get("scissor"); 
        for l in [x for x in stack if x > lvl]: stack.pop(l)
        continue
    if e.get("type") == "group":
        for l in [x for x in stack if x >= lvl]: stack.pop(l)
        continue
    for l in [x for x in stack if x >= lvl]: stack.pop(l)
    sc = None
    for l in sorted(stack):
        r = stack[l]
        if r is not None: sc = r if sc is None else (sc & r)
    r = e["rect"]; vis = True if sc is None else (sc.intersects(r) and not (sc & r).is_empty)
    out.append(dict(k=k, level=lvl, rect=[round(v, 3) for v in r], scissor=None if sc is None else [round(v, 3) for v in sc], visible=vis, kinds="".join(sorted(set(x[0] for x in e["items"]))), fill=e.get("fill"), color=e.get("color"), width=e.get("width"), dashes=e.get("dashes"))); k += 1
inv = [o for o in out if not o["visible"]]
print("paths", len(out), "invisible (outside scissor)", len(inv))
for o in inv[:40]: print("  INVIS", o["k"], o["level"], o["rect"], "scissor", o["scissor"], o["kinds"], o["dashes"], o["width"])
# twin analysis: items 3.629 pt lower than another item with the same x range, kinds, fill, color, width
tw = []
for o in out:
    for q in out:
        if q is o: continue
        if abs(q["rect"][0] - o["rect"][0]) < 0.01 and abs(q["rect"][2] - o["rect"][2]) < 0.01 and abs((q["rect"][1] - o["rect"][1]) - 3.629) < 0.01 and q["kinds"] == o["kinds"] and q["fill"] == o["fill"] and q["color"] == o["color"]:
            tw.append((o["k"], q["k"], o["rect"], q["visible"]))
print("twins (upper k, lower k, upper rect, lower visible):", len(tw)); [print("  ", t) for t in tw[:40]]
json.dump(out, open(sys.argv[2], "w"), indent=0)
