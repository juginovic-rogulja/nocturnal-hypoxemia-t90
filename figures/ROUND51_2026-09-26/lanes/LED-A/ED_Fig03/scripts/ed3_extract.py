"""Read a graded-panel sheet's vectors (current ED 3 or the round-27 rebuild): per panel the axes box,
y tick marks and labels, marker centres, ribbon vertices, white/blue line, asterisk, title, and the
marker and ribbon values in data units through the panel's own y ticks."""
import fitz, gc


def spans_of(pg):
    out = []
    for b in pg.get_text("dict")["blocks"]:
        if b["type"] != 0:
            continue
        for l in b["lines"]:
            for s in l["spans"]:
                out.append(dict(bbox=[round(v, 3) for v in s["bbox"]], origin=[round(v, 3) for v in s["origin"]],
                                size=round(s["size"], 2), text=s["text"].replace("\xa0", " "), font=s["font"],
                                color=s.get("color"), dir=[round(v, 2) for v in l["dir"]]))
    return out


def dedupe_spans(spans):
    """The Illustrator-flattened sheet paints some text twice through overlapping clips, and the
    clipped copy reads back as a prefix ('**' under '***', '5-10 >1' under '5-10 >10'). Keep one span
    per (origin, size, font), the longest text."""
    best = {}
    for s in spans:
        k = (round(s["origin"][0] * 2) / 2, round(s["origin"][1] * 2) / 2, s["size"], s["font"])
        if k not in best or len(s["text"]) > len(best[k]["text"]):
            best[k] = s
    return list(best.values())


def ink(c):
    return c is not None and max(c) < 0.15


def extract(path):
    doc = fitz.open(path)
    pg = doc[0]
    W, H = pg.rect.width, pg.rect.height
    spans = dedupe_spans(spans_of(pg))
    dr = pg.get_drawings()
    # every straight ink segment on the page, whether its path holds one segment (Illustrator sheet)
    # or many (the round-27 sheet draws spines and ticks of a panel as one path)
    segs = []
    for g in dr:
        if g["color"] is not None and ink(g["color"]):
            for it in g["items"]:
                if it[0] == "l":
                    segs.append((it[1], it[2], g.get("width") or 0))
    spines = set()
    for a, b, wd in segs:
        if abs(wd - 0.8) < 0.01:
            if abs(a.y - b.y) < 0.01 and abs(a.x - b.x) > 60:
                spines.add(("h", round(min(a.x, b.x), 1), round(a.y, 1), round(max(a.x, b.x), 1)))
            if abs(a.x - b.x) < 0.01 and abs(a.y - b.y) > 40:
                spines.add(("v", round(a.x, 1), round(min(a.y, b.y), 1), round(max(a.y, b.y), 1)))
    hsp = sorted(s for s in spines if s[0] == "h")
    vsp = sorted(s for s in spines if s[0] == "v")
    panels = []
    for (_, x0, yb, x1) in hsp:
        v = [s for s in vsp if abs(s[1] - x0) < 0.6 and abs(s[3] - yb) < 0.6]
        assert len(v) == 1, (x0, yb, v)
        panels.append(dict(left=x0, right=x1, bottom=yb, top=v[0][2]))
    panels.sort(key=lambda p: (round(p["top"]), p["left"]))

    def in_box(pt, p, padx=0, pady=0):
        return p["left"] - padx <= pt[0] <= p["right"] + padx and p["top"] - pady <= pt[1] <= p["bottom"] + pady

    out = []
    for pi, p in enumerate(panels):
        L, R, T, B = p["left"], p["right"], p["top"], p["bottom"]
        tl = [s for s in spans if "Bold" not in s["font"] and s["size"] == 9.5 and abs(s["origin"][0] - L) < 0.6
              and T - 32 < s["origin"][1] < T]
        tl = sorted({(s["origin"][1], s["text"]) for s in tl})
        title = " ".join(t for _, t in tl)
        yt = set()
        for a, b, wd in segs:
            if abs(a.y - b.y) < 0.01 and abs(max(a.x, b.x) - L) < 0.6 and 2.0 < abs(a.x - b.x) < 3.0 \
                    and T - 0.1 <= a.y <= B + 0.1:
                yt.add(round(a.y, 2))
        yt = sorted(yt)
        ylab = {}
        for y in yt:
            c = [s for s in spans if "Bold" not in s["font"] and s["size"] == 9.5 and s["bbox"][2] < L
                 and s["bbox"][2] > L - 30 and abs((s["origin"][1] - 0.377 * 9.5) - y) < 1.5]
            txt = sorted({s["text"] for s in c}, key=len)
            ylab[y] = txt[-1] if txt else None
        vals = [(y, float(ylab[y])) for y in yt if ylab[y] is not None]
        assert len(vals) >= 2, (title, yt, ylab)
        n = len(vals)
        sy = sum(v[0] for v in vals); sv = sum(v[1] for v in vals)
        syy = sum(v[0] ** 2 for v in vals); syv = sum(v[0] * v[1] for v in vals)
        slope = (n * syv - sy * sv) / (n * syy - sy * sy); icpt = (sv - slope * sy) / n
        val = lambda y: icpt + slope * y
        resid = max(abs(val(y) - v) for y, v in vals)
        ymin, ymax = val(B), val(T)
        xt = sorted({round(a.x, 2) for a, b, wd in segs
                     if abs(a.x - b.x) < 0.01 and 2.0 < abs(a.y - b.y) < 3.0
                     and abs(min(a.y, b.y) - B) < 0.6 and L - 0.1 <= a.x <= R + 0.1})
        mk = {}
        for g in dr:
            if g["fill"] is not None and len(g["items"]) in (4, 8) and all(it[0] == "c" for it in g["items"]):
                r = g["rect"]; cx, cy = (r.x0 + r.x1) / 2, (r.y0 + r.y1) / 2
                if in_box((cx, cy), p, 1, 3):
                    mk[round(cx, 1)] = dict(cx=round(cx, 3), cy=round(cy, 3), r=round((r.x1 - r.x0) / 2, 3),
                                            fill=[round(c, 4) for c in g["fill"]],
                                            stroke=None if g["color"] is None else [round(c, 3) for c in g["color"]],
                                            sw=g.get("width"))
        markers = [mk[k] for k in sorted(mk)]
        rib = []
        for g in dr:
            if g["fill"] is not None and (g.get("fill_opacity") or 1) < 0.5 and in_box((g["rect"].x0, g["rect"].y0), p, 1, 1):
                pts = []
                for it in g["items"]:
                    if it[0] == "l":
                        pts += [(round(it[1].x, 3), round(it[1].y, 3)), (round(it[2].x, 3), round(it[2].y, 3))]
                dd = []
                [dd.append(q) for q in pts if q not in dd]
                key = tuple(dd)
                if key not in [tuple(r["pts"]) for r in rib]:
                    rib.append(dict(fill=[round(c, 4) for c in g["fill"]], opacity=round(g.get("fill_opacity"), 4), pts=dd))
        line = []
        for g in dr:
            if g["color"] is not None and g.get("width") and abs(g["width"] - 1.5) < 0.01 and in_box((g["rect"].x0, g["rect"].y0), p, 1, 3):
                pts = []
                for it in g["items"]:
                    if it[0] == "l":
                        pts += [(round(it[1].x, 3), round(it[1].y, 3)), (round(it[2].x, 3), round(it[2].y, 3))]
                dd = []
                [dd.append(q) for q in pts if q not in dd]
                if dd not in [l["pts"] for l in line]:
                    line.append(dict(color=[round(c, 3) for c in g["color"]], pts=dd, cap=g.get("lineCap"), join=g.get("lineJoin")))
        st = sorted({(s["origin"][0], s["origin"][1], s["text"], s.get("color")) for s in spans
                     if "*" in s["text"] and L <= s["origin"][0] <= R and T <= s["origin"][1] <= B}, key=lambda t: -len(t[2]))
        star = st[0][2] if st else ""
        star_color = st[0][3] if st else None
        ref = []
        for g in dr:
            if g["color"] is not None and g.get("dashes") and g["dashes"] not in ("[] 0", "[]", None) and in_box((g["rect"].x0, g["rect"].y0), p, 0.6, 0.6):
                ref.append(dict(y=round(g["items"][0][1].y, 3), color=[round(c, 3) for c in g["color"]], width=g.get("width"), dashes=g["dashes"]))
        band_fills = set()
        for g in dr:
            if g["fill"] is not None and len(g["items"]) == 1 and g["items"][0][0] == "re" and (g.get("fill_opacity") or 1) > 0.9:
                r = g["rect"]
                if abs(r.y0 - T) < 0.6 and abs(r.y1 - B) < 0.6 and r.x0 >= L - 0.6 and r.x1 <= R + 0.6 and (r.x1 - r.x0) < (R - L) * 0.3:
                    band_fills.add((round(r.x0, 2), round(r.x1, 2), tuple(round(c, 4) for c in g["fill"])))
        band_fills = sorted(band_fills)
        hr = [round(val(m["cy"]), 4) for m in markers]
        ribv = [[(x, round(val(y), 4)) for x, y in r["pts"]] for r in rib]
        out.append(dict(index=pi, title=title, axes=dict(left=L, right=R, top=T, bottom=B),
                        y_ticks=sorted({(y, ylab[y]) for y in yt}), y_fit_resid=round(resid, 5),
                        ymin=round(ymin, 5), ymax=round(ymax, 5), x_ticks=xt, markers=markers, marker_hr=hr,
                        ribbons=rib, ribbon_values=ribv, line=line, star=star, star_color=star_color,
                        ref=ref, ref_value=[round(val(r["y"]), 4) for r in ref], band_fills=band_fills))
    res = dict(path=path, page=[W, H], mediabox=list(pg.mediabox), panels=out, spans=spans,
               fonts=[f[3] for f in pg.get_fonts()], n_drawings=len(dr))
    doc.close(); gc.collect()
    return res
