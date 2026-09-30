#!$T90_PY
"""Round 49, lane LED-A (2026-09-26): shared helpers for the +1 pt rebuilds of ED_Fig01 to ED_Fig05.
watched(cmd, log): one child process at a time, RSS polled every second, killed above RSS_KB or after max_s (the round-38 watchdog idiom).
gs_render(pdf, png, dpi): Ghostscript png16m render in a watched child (the only renderer of this lane).
page_size(pdf): pypdf mediabox. sha256(path). The brief's verification helpers: page_check, census_check (tools/census.py)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import hashlib, json, os, subprocess, sys, time

LANE = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LED-A"
V26 = f"{paths.FIGURE_ROOT}/ROUND47_2026-09-24/figures/NEW_FINAL_SET_V26"
GS = paths.GS; PY = paths.PY
DELTA = 1.0                      # the round: every text element one point larger
LOGS = f"{LANE}/logs"; os.makedirs(LOGS, exist_ok=True)
sys.path.insert(0, f"{LANE}/tools")
from census import census, compare   # noqa: E402


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""): h.update(c)
    return h.hexdigest()


def hydrated(p):
    st = os.stat(p); assert not (st.st_size > 0 and st.st_blocks == 0), f"EVICTED (0 blocks): {p}"; return p


def watched(cmd, log, max_s=900, rss_kb=4_000_000, env=None):
    """Run ONE child, poll the tree's RSS every second, SIGKILL above rss_kb or after max_s. Returns (rc, peak_kb, wall_s)."""
    with open(log, "w") as lf:
        p = subprocess.Popen(cmd, stdout=lf, stderr=subprocess.STDOUT, env=env)
    t0 = time.time(); peak = 0
    while p.poll() is None:
        time.sleep(1)
        try:
            out = subprocess.run(["ps", "-o", "rss=", "-g", str(p.pid)], capture_output=True, text=True).stdout.split()
            rss = sum(int(v) for v in out) if out else int(subprocess.run(["ps", "-o", "rss=", "-p", str(p.pid)], capture_output=True, text=True).stdout.strip() or 0)
        except Exception:
            rss = 0
        peak = max(peak, rss)
        if rss > rss_kb or time.time() - t0 > max_s:
            p.kill(); raise RuntimeError(f"watchdog killed {cmd[:3]} rss {rss} KB after {time.time() - t0:.0f} s (log {log})")
    return p.returncode, peak, time.time() - t0


def run_step(name, cmd, log_dir, max_s=900):
    """A watched child that must exit 0; prints the tail of its log; returns the log path."""
    log = f"{log_dir}/{name}.log"
    rc, peak, wall = watched(cmd, log, max_s=max_s)
    tail = open(log, errors="replace").read()
    print(f"[{name}] rc={rc} peak_rss={peak / 1024:.0f} MB wall={wall:.0f} s\n" + tail[-2500:])
    if rc != 0: raise SystemExit(f"STEP FAILED {name} (rc {rc}): see {log}")
    return log


def gs_render(pdf, png, dpi, log_dir, name=None, max_s=900):
    hydrated(pdf)
    if os.path.exists(png): os.remove(png)
    run_step(name or f"gs{dpi}_{os.path.basename(pdf)[:-4]}", [GS, "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=png16m", f"-r{dpi}", f"-sOutputFile={png}", pdf], log_dir, max_s=max_s)
    assert os.path.exists(png) and os.path.getsize(png) > 0, png
    return png


def page_size(pdf):
    from pypdf import PdfReader
    r = PdfReader(hydrated(pdf)); mb = r.pages[0].mediabox
    return len(r.pages), float(mb.width), float(mb.height)


def page_check(new_pdf, v26_pdf, tol=0.05):
    n, w, h = page_size(new_pdf); n0, w0, h0 = page_size(v26_pdf)
    ok = n == 1 == n0 and abs(w - w0) <= tol and abs(h - h0) <= tol
    return ok, f"{'PASS' if ok else 'FAIL'} page size equals V26 within {tol} pt: new {w:.3f} x {h:.3f} pt, V26 {w0:.3f} x {h0:.3f} pt, {n} page"


def census_check(new_pdf, v26_pdf, work_dir, sheet, exceptions=None, delta=DELTA, merge_gap=1.5, removed=(), added=(), v26_widened=False):
    """Ghostscript txtwrite census of both sheets, compare: same multiset of strings (apart from declared re-wrapped lines), every size + delta
    (exceptions: {text: new size}). v26_widened: the V26 sheet is the round-43 two-copy composite (ED_Fig01), deduped by the split rule."""
    rn, ro = [], []
    new = merge_touching(census(new_pdf, f"{work_dir}/{sheet}_new_census.xml"), merge_gap, rn)
    old = census(v26_pdf, f"{work_dir}/{sheet}_v26_census.xml")
    n_raw = len(old)
    if v26_widened: old = dedupe_widened(old)
    old = merge_touching(old, merge_gap, ro)
    json.dump(new, open(f"{work_dir}/{sheet}_new_census.json", "w"), indent=0); json.dump(old, open(f"{work_dir}/{sheet}_v26_census.json", "w"), indent=0)
    ok, lines = compare(old, new, delta=delta, exceptions=exceptions, removed=removed, added=added)
    lines.append(f"INFO census spans: V26 {n_raw} raw" + (f" -> {len(old)} after the widened-composite dedupe" if v26_widened else "") + f", new {len(new)}; glyph-split spans joined: V26 {ro}, new {rn}")
    return ok, lines, old, new


def dedupe_widened(items, split=160.0, g=14.4, tol=0.6):
    """The round-43 widened ED_Fig01 (two clipped show_pdf_page copies of one page) lists every string twice in gs txtwrite: at its source x and
    at x + g. The visible copy is the one at x when x < split (left clip), the one at x + g otherwise (right clip)."""
    def nz(t): return t.replace("\xa0", " ").strip()
    out = []; used = [False] * len(items)
    for i, a in enumerate(items):
        if used[i]: continue
        twin = None
        for j in range(i + 1, len(items)):
            b = items[j]
            if not used[j] and nz(b["text"]) == nz(a["text"]) and abs(b["y"] - a["y"]) < tol and abs(b["size"] - a["size"]) < 0.01 and abs(abs(b["x"] - a["x"]) - g) < tol:
                twin = j; break
        if twin is None: out.append(dict(a)); used[i] = True; continue
        lo, hi = (a, items[twin]) if a["x"] <= items[twin]["x"] else (items[twin], a)
        out.append(dict(lo if lo["x"] < split else hi)); used[i] = used[twin] = True
    return out


def merge_touching(items, gap=1.5, report=None):
    """gs txtwrite splits one string at some glyph pairs ('V' + 'entricular ...', 'T' + 'ime ...') and rounds every bbox to whole points:
    join consecutive spans of one baseline (horizontal: same y; rotated: same x, running upward), same size and font, whose gap along the
    text direction is within gap pt (whole-point rounding: a split shows as -1 to 1; the sheets' separate neighbouring strings sit 2.7 pt or more apart)."""
    def geom(i):
        x0, y0, x1, y1 = i["bbox"]
        if abs(x1 - x0) < 0.01 and abs(y1 - y0) > 0.01: return ("v", round(x0, 1), -y0, -y1)      # rotated 90: text runs upward
        return ("h", round(y0, 1), x0, x1)
    items = sorted(items, key=lambda i: (geom(i)[0], geom(i)[1], geom(i)[2]))
    out = []
    for it in items:
        if out:
            p = out[-1]; gp, gi = geom(p), geom(it)
            if gp[0] == gi[0] and abs(gp[1] - gi[1]) < 0.6 and abs(p["size"] - it["size"]) < 0.01 and p["font"] == it["font"] and gi[2] - gp[3] <= gap:      # overlapping or touching (whole-point bboxes): one string
                if report is not None: report.append((p["text"], it["text"]))
                p["text"] += it["text"]; p["bbox"] = [p["bbox"][0], p["bbox"][1], it["bbox"][2], it["bbox"][3]]; continue
        out.append(dict(it))
    return out


def positions_report(old, new, movers=(), tol=6.0):
    """Every string's origin on the new sheet against V26 (after the census merge): left-aligned strings sit at the same x, y shifts are
    declared movers. Returns (ok, lines): FAIL when an undeclared string moved more than tol pt."""
    def nz(t): return t.replace("\xa0", " ").strip()
    pool = {}
    for i in old: pool.setdefault(nz(i["text"]), []).append(i)
    worst = 0.0; bad = []; n = 0
    for i in new:
        c = pool.get(nz(i["text"]), [])
        if not c: continue
        b = min(c, key=lambda o: abs(o["x"] - i["x"]) + abs(o["y"] - i["y"]))
        dx, dy = i["x"] - b["x"], i["y"] - b["y"]; n += 1
        if nz(i["text"]) in movers: continue
        worst = max(worst, abs(dx), abs(dy))
        if abs(dx) > tol or abs(dy) > tol: bad.append((nz(i["text"])[:30], round(dx, 1), round(dy, 1)))
    ok = not bad
    return ok, [f"{'PASS' if ok else 'FAIL'} positions: {n} strings matched to V26, largest undeclared origin shift {worst:.2f} pt (limit {tol} pt, text grows from its anchor){': ' + str(bad[:8]) if bad else ''}"]


def write_checks(path, lines):
    body = [l.replace("; ", ", ").replace(";", ",").replace("—", ", ") for l in lines]      # the brief: no semicolon, no em dash in any string this lane adds
    n_fail = sum(1 for l in body if l.startswith("FAIL")); n_pass = sum(1 for l in body if l.startswith("PASS"))
    res = "RESULT ALL PASS" if n_fail == 0 and n_pass > 0 else "RESULT FAIL"
    os.makedirs(os.path.dirname(path), exist_ok=True)
    open(path, "w").write("\n".join(body + ["", f"{n_pass} PASS, {n_fail} FAIL", res]) + "\n")
    print(res); return res == "RESULT ALL PASS"
