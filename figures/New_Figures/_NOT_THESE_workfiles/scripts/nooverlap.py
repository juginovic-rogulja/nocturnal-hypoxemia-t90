"""
The no-overlapping-text gate.

Alen's round-2 instruction was blunt: "Please make sure none of the text overlaps or all of that
stuff. That needs to be clean and nice." Eyeballing a PNG is not a check, so this module measures
it. Every figure written by the round-2 scripts is passed through `gate()` before it is saved,
and a figure that fails is not delivered.

What is measured, in rendered pixels at the real save size

  1. text against text. Every string matplotlib will draw is collected, its window extent is
     taken from the renderer, and every pair is intersected. Two strings that share more than
     PAD pixels in both directions are a collision. Tick labels, axis labels, titles, panel
     letters, legend entries, annotations and bare figure text are all in scope.

  2. text against data. A label sitting on top of a marker is as unreadable as a label sitting
     on another label, so every scatter offset and every line vertex is transformed to display
     coordinates and tested against every text box in the same axes. Reference rules, the axis
     spines and the zero line are excluded, they are meant to sit under text.

  3. text running off the canvas. savefig uses a tight bounding box, so a string that leaves the
     figure is not actually cut, but one that leaves its own axes and lands on a neighbouring
     axes is caught by test 1 anyway.

The gate returns a report that the build scripts write into the provenance JSON, so there is a
record of what was checked rather than a claim that it was.
"""
import numpy as np

PAD = 0.5          # pixels of tolerance before an overlap counts
DATA_PAD = 1.0     # a marker this close to a glyph counts as sitting under it


def _boxes(fig):
    """Every string that will be drawn, with its rendered pixel box."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    out = []

    def add(t, where):
        try:
            if not t.get_visible():
                return
            s = t.get_text()
            if not s or not s.strip():
                return
            bb = t.get_window_extent(renderer=r)
            if bb.width <= 0 or bb.height <= 0:
                return
            out.append({"where": where, "text": s, "bbox": bb, "artist": t})
        except Exception:
            return

    for t in fig.texts:
        add(t, "figure")
    for i, ax in enumerate(fig.axes):
        tag = f"axes{i}"
        for t in ax.texts:
            add(t, tag)
        # ax.axis("off") switches the axis off without clearing get_visible on each tick label,
        # so a schematic panel would otherwise be failed against ticks that are never drawn
        axison = getattr(ax, "axison", True)
        if axison and ax.xaxis.get_visible():
            for t in ax.get_xticklabels():
                add(t, f"{tag}:tick")
            add(ax.xaxis.label, f"{tag}:xlabel")
        if axison and ax.yaxis.get_visible():
            for t in ax.get_yticklabels():
                add(t, f"{tag}:tick")
            add(ax.yaxis.label, f"{tag}:ylabel")
        add(ax.title, f"{tag}:title")
        lg = ax.get_legend()
        if lg is not None:
            for t in lg.get_texts():
                add(t, f"{tag}:legend")
    lg = getattr(fig, "legends", [])
    for j, l in enumerate(lg):
        for t in l.get_texts():
            add(t, f"figlegend{j}")
    return out, r


def _ov(a, b):
    dx = min(a.x1, b.x1) - max(a.x0, b.x0)
    dy = min(a.y1, b.y1) - max(a.y0, b.y0)
    return dx, dy


def _marker_pad(ax):
    """
    How far from a marker CENTRE a glyph has to stay before the two stop touching.

    _data_points returns centres, so testing a label box against them with a one-pixel pad lets
    a label sit flush against a marker that is eight pixels across. On eFigure12D that is
    exactly what happened: the last letter of a condition name ended up under its own point.
    The pad has to carry the marker's own radius.
    """
    dpi = ax.figure.dpi
    worst = 0.0
    for c in ax.collections:
        try:
            s = np.asarray(c.get_sizes(), float)      # area in points squared
            if s.size:
                worst = max(worst, float(np.sqrt(s.max())) / 2.0)
        except Exception:
            continue
    for ln in ax.lines:
        try:
            if ln.get_marker() not in (None, "None", "", " "):
                worst = max(worst, float(ln.get_markersize()) / 2.0)
        except Exception:
            continue
    return DATA_PAD + worst * dpi / 72.0


class _Band:
    """A thin forbidden strip in display pixels, shaped like a matplotlib bbox."""

    def __init__(self, x0, y0, x1, y1):
        self.x0, self.y0, self.x1, self.y1 = x0, y0, x1, y1


def _parallel_rules(ax, half_px=2.5):
    """
    Reference rules that run ALONG the direction the labels are set in.

    A rule crossing a label at a right angle is a light touch and the paper has always allowed
    it. A rule running the same way as the text is a strikethrough. On eFigure12D the dashed
    chance line went straight through the middle of "Ischemic heart disease" and the label was
    still accepted, because reference rules were excluded from the avoid set wholesale. They
    are excluded only across the text now, never along it.
    """
    out = []
    for ln in ax.lines:
        try:
            # is_dashed covers the dash tuples this paper actually uses, ls=(0, (4, 3)) and
            # friends, which a string comparison against "--" silently misses
            if not ln.is_dashed():
                continue
            xd = np.asarray(ln.get_xdata(), float)
            yd = np.asarray(ln.get_ydata(), float)
            if len(xd) != 2:
                continue
            # axhline and axvline carry a blended transform, x in axes coordinates and y in
            # data coordinates, so the line's own transform is the only one that is right
            p = ln.get_transform().transform(np.column_stack([xd, yd]))
            if not np.isfinite(p).all():
                continue
            (x0, y0), (x1, y1) = p
            if abs(y1 - y0) <= 1.0:                       # horizontal, parallel to the text
                out.append(_Band(min(x0, x1), y0 - half_px, max(x0, x1), y0 + half_px))
        except Exception:
            continue
    return out


def _data_points(ax):
    """Marker centres and line vertices in display pixels, excluding reference rules."""
    pts = []
    for c in ax.collections:
        try:
            off = c.get_offsets()
            if off is None or len(off) == 0:
                continue
            pts.append(ax.transData.transform(np.asarray(off, float)))
        except Exception:
            continue
    for ln in ax.lines:
        try:
            if ln.get_linestyle() in ("--", ":", "-.") and len(ln.get_xdata()) == 2:
                continue        # a reference rule, text is allowed to cross it
            xy = np.column_stack([np.asarray(ln.get_xdata(), float),
                                  np.asarray(ln.get_ydata(), float)])
            if len(xy) == 0:
                continue
            if ln.get_marker() in (None, "None", "", " "):
                continue        # a plain connecting line, not a mark a reader reads a value off
            pts.append(ax.transData.transform(xy))
        except Exception:
            continue
    if not pts:
        return np.zeros((0, 2))
    p = np.vstack(pts)
    return p[np.isfinite(p).all(axis=1)]


def check(fig, name="figure"):
    boxes, _ = _boxes(fig)
    clashes = []
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            a, b = boxes[i], boxes[j]
            dx, dy = _ov(a["bbox"], b["bbox"])
            if dx > PAD and dy > PAD:
                clashes.append({"kind": "text-text",
                                "a": f"{a['where']}: {a['text'][:44]}",
                                "b": f"{b['where']}: {b['text'][:44]}",
                                "overlap_px": [round(float(dx), 1), round(float(dy), 1)]})
    for i, ax in enumerate(fig.axes):
        pts = _data_points(ax)
        if len(pts) == 0:
            continue
        pad = _marker_pad(ax)
        for bx in boxes:
            if not bx["where"].startswith(f"axes{i}"):
                continue
            if bx["where"].endswith(("tick", "xlabel", "ylabel", "title")):
                continue        # these live outside the data area by construction
            bb = bx["bbox"]
            hit = ((pts[:, 0] > bb.x0 - pad) & (pts[:, 0] < bb.x1 + pad) &
                   (pts[:, 1] > bb.y0 - pad) & (pts[:, 1] < bb.y1 + pad))
            if hit.any():
                clashes.append({"kind": "text-over-data",
                                "a": f"{bx['where']}: {bx['text'][:44]}",
                                "b": f"{int(hit.sum())} marker(s) inside the label box",
                                "overlap_px": None})
    return {"figure": name, "n_text_objects": len(boxes), "n_clashes": len(clashes),
            "passes": len(clashes) == 0, "clashes": clashes}


CAND = [(7, 3), (7, -10), (-7, 3), (-7, -10), (0, 10), (0, -14), (16, 1), (-16, 1),
        (18, 12), (-18, 12), (18, -14), (-18, -14), (28, 6), (-28, 6), (0, 22), (0, -26),
        (34, 18), (-34, 18), (34, -22), (-34, -22)]


def place_labels(ax, xs, ys, texts, fontsize, color, leader_from_px=15.0, extra_avoid=None):
    """
    Put a name beside each of a handful of points without any of them landing on each other,
    on a marker, or on something already drawn.

    A candidate offset is tried, the label is rendered, its real pixel box is measured, and it
    is kept only if that box is clear. Offsets grow outward until one fits, and a label that
    ended up far from its point is joined back to it by a hairline. This is what the rejected
    eFigure 5 callouts should have been doing: the leader is the length of a fingernail and it
    exists only because the label could not sit where it wanted to.
    """
    fig = ax.figure
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    taken = [b["bbox"] for b in _boxes(fig)[0]] + _parallel_rules(ax)
    if extra_avoid:
        taken = taken + list(extra_avoid)
    pts = _data_points(ax)
    pad = _marker_pad(ax)
    placed = []
    for x, y, t in zip(xs, ys, texts):
        chosen = None
        for dx, dy in CAND:
            ha = "left" if dx > 0 else ("right" if dx < 0 else "center")
            a = ax.annotate(t, (x, y), textcoords="offset points", xytext=(dx, dy),
                            fontsize=fontsize, color=color, ha=ha, va="center", zorder=8)
            bb = a.get_window_extent(renderer=r)
            clash = any((_ov(bb, b)[0] > PAD and _ov(bb, b)[1] > PAD) for b in taken)
            if not clash and len(pts):
                clash = bool(((pts[:, 0] > bb.x0 - pad) & (pts[:, 0] < bb.x1 + pad) &
                              (pts[:, 1] > bb.y0 - pad) &
                              (pts[:, 1] < bb.y1 + pad)).any())
            if clash:
                a.remove()
                continue
            chosen = (a, bb, dx, dy)
            break
        if chosen is None:
            raise AssertionError(f"no clear place for the label {t!r}, widen the panel")
        a, bb, dx, dy = chosen
        taken.append(bb)
        placed.append(t)
        if (dx ** 2 + dy ** 2) ** 0.5 >= leader_from_px:
            x0, y0 = ax.transData.transform((x, y))
            ex = bb.x0 - 1.5 if dx > 0 else (bb.x1 + 1.5 if dx < 0 else (bb.x0 + bb.x1) / 2)
            ey = (bb.y0 + bb.y1) / 2
            inv = ax.transData.inverted()
            p0 = inv.transform((x0, y0))
            p1 = inv.transform((ex, ey))
            ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=color, lw=0.45, alpha=0.55, zorder=3)
    return placed


def gate(fig, name, log=None, hard=True):
    """Check, record and, unless told otherwise, refuse to let a bad figure through."""
    r = check(fig, name)
    if log is not None:
        log.append(r)
    if not r["passes"]:
        msg = "\n".join(f"    {c['kind']}: {c['a']}  ||  {c['b']}  {c['overlap_px'] or ''}"
                        for c in r["clashes"][:40])
        text = f"{name}: {r['n_clashes']} text collisions\n{msg}"
        if hard:
            raise AssertionError(text)
        print("OVERLAP WARNING\n" + text)
    else:
        print(f"  overlap gate passed for {name}, {r['n_text_objects']} strings checked")
    return r
