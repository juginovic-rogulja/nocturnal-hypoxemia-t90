#!$T90_PY
"""LED-B round 51 copy (base set V28, lane ROUND51). LED-B round 49 (2026-09-26): shared helpers for the +1 pt rebuild of ED_Fig06 to ED_Fig10. Paths, eviction gate, sha256, the watchdog
runner (one PDF-opening child at a time under wd_run.sh), Ghostscript txtwrite census (parse and compare against V26), Ghostscript renders.
Never imports PyMuPDF: every PyMuPDF step is a separate child script run through wd()."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import hashlib, json, os, re, html, subprocess, time
from collections import Counter

T90 = paths.FIGURE_ROOT
LANE = f"{paths.FIGURE_ROOT}/ROUND52_2026-09-26/lanes/LED-B"   # round 52 (2026-09-26): this lane
V26 = f"{paths.FIGURE_ROOT}/ROUND51_2026-09-26/figures/NEW_FINAL_SET_V29"   # round 52: the base set is V29 (ED_Fig07 = the round-51 output); the name V26 is kept for the round-49 call sites
V28 = V26; V29 = V26
PY = paths.PY; GS = paths.GS; WD = f"{LANE}/scripts/wd_run.sh"
PT_PLUS = 1.0          # the round-49 rule: every text one point larger (Alen, 2026-09-26)
STRIP = 16.0; TITLE_SIZE = 13.0 + PT_PLUS; TITLE_Y = 12.0
V26_TITLE_X = {"ED_Fig06": 12.96, "ED_Fig07": 12.96, "ED_Fig08": 12.96, "ED_Fig09": 11.52, "ED_Fig10": 12.96}   # round-40 stamp records
INK = "#1a1d21"


def hydrated(path):
    st = os.stat(path); assert not (st.st_size > 0 and st.st_blocks == 0), f"EVICTED (0 blocks): {path}"; return path


def sha256(path):
    h = hashlib.sha256()
    with open(hydrated(path), "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""): h.update(chunk)
    return h.hexdigest()


def wd(name, max_s, cmd, rss_kb=4000000):
    """Run one command as a watchdog child (wd_run.sh: RSS and wall-time caps, one child at a time). Raises on rc != 0."""
    env = dict(os.environ, WD_LOGDIR=f"{LANE}/logs/watchdog"); os.makedirs(env["WD_LOGDIR"], exist_ok=True)
    subprocess.run([WD, name, str(max_s), str(rss_kb), "--"] + cmd, env=env, capture_output=True, text=True)
    st = open(f"{env['WD_LOGDIR']}/{name}.status").read().strip(); out = open(f"{env['WD_LOGDIR']}/{name}.out").read()
    print(f"[{name}] {st}"); print(out.rstrip()[-4000:])
    if not st.startswith("rc=0"): raise SystemExit(f"STEP FAILED {name}: {st}\n{out[-2000:]}")
    return out


def gs_render(pdf, png, dpi=150):
    r = subprocess.run([GS, "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=png16m", f"-r{dpi}", f"-sOutputFile={png}", hydrated(pdf)], capture_output=True, text=True)
    assert r.returncode == 0 and os.path.exists(png), r.stderr[-500:]; return png


def gs_txtwrite(pdf, xml):
    r = subprocess.run([GS, "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=txtwrite", "-dTextFormat=0", f"-sOutputFile={xml}", hydrated(pdf)], capture_output=True, text=True)
    assert r.returncode == 0 and os.path.exists(xml), r.stderr[-500:]; return xml


def parse_txtwrite(xml):
    """The <span> elements of a txtwrite XML: text (NBSP as space), font, size, bbox [x0, y, x1] (y = baseline, top-down points)."""
    out = []
    for m in re.finditer(r'<span bbox="([^"]*)" font="([^"]*)" size="([^"]*)">(.*?)</span>', open(xml, encoding="utf-8", errors="replace").read(), re.S):
        bbox, font, size, body = m.groups()
        chars = [html.unescape(c) for c in re.findall(r'<char bbox="[^"]*" c="([^"]*)"/>', body)]
        txt = "".join(chars).replace("\xa0", " ").replace("\xad", "-")
        if not txt.strip(): continue
        b = [float(v) for v in bbox.split()]
        cur = dict(text=txt, font=font, size=round(float(size), 3), x0=b[0], y=b[1], x1=b[2], bold="Bold" in font)
        # txtwrite splits one string at a kerning pair ("T" + "ype 2 diabetes") or a glyph-width change: merge a span that touches the
        # previous one (same font, size and baseline, x0 within 0.6 pt of the previous x1) back into it
        if out and out[-1]["font"] == font and abs(out[-1]["size"] - cur["size"]) < 0.001 and abs(out[-1]["y"] - cur["y"]) < 0.6 and -2.5 <= cur["x0"] - out[-1]["x1"] <= 0.6:   # a kerned pair overlaps by up to about 1.5 pt at 10 pt
            out[-1]["text"] += txt; out[-1]["x1"] = max(out[-1]["x1"], cur["x1"]); continue
        out.append(cur)
    return out


def census(sheet, new_pdf, old_pdf, out_dir, exceptions=None, plus=PT_PLUS):
    """Ghostscript text census of the new sheet against V26: same multiset of strings, every size exactly +plus, with declared exceptions.
    exceptions: dict(removed=[{"text","size","count","why"}], added=[{"text","size","count","why"}], sizes=[{"text","old","new","count","why"}]).
    Writes <out_dir>/census.txt and census.json; returns (ok, lines)."""
    ex = exceptions or {}; os.makedirs(out_dir, exist_ok=True)
    nx, ox = f"{out_dir}/{sheet}_NEW_txtwrite.xml", f"{out_dir}/{sheet}_V26_txtwrite.xml"
    gs_txtwrite(new_pdf, nx); gs_txtwrite(old_pdf, ox); new, old = parse_txtwrite(nx), parse_txtwrite(ox)
    lines = [f"# Ghostscript txtwrite census, {sheet}: NEW {new_pdf} ({len(new)} spans) against V26 {old_pdf} ({len(old)} spans), rule size +{plus:g}"]
    cn, co = Counter(s["text"] for s in new), Counter(s["text"] for s in old)
    exp_rem = Counter(); exp_add = Counter()
    for e in ex.get("removed", []): exp_rem[e["text"]] += e.get("count", 1)
    for e in ex.get("added", []): exp_add[e["text"]] += e.get("count", 1)
    lost, gained = co - cn, cn - co
    ok_str = (lost == exp_rem) and (gained == exp_add)
    lines.append(f"{'PASS' if ok_str else 'FAIL'}  string multiset: NEW = V26 except the declared strings (lost {sum(lost.values())}: {dict(lost)}; gained {sum(gained.values())}: {dict(gained)}; declared removed {dict(exp_rem)}, declared added {dict(exp_add)})")
    # sizes: match every V26 span (text) to a NEW span with the expected size; exceptions by text
    size_ex = {}
    for e in ex.get("sizes", []): size_ex.setdefault(e["text"], []).append(e)
    pool = Counter((s["text"], s["size"], s["bold"]) for s in new); bad = []; matched = 0; ex_used = Counter()
    rem_left = Counter(exp_rem)
    for s in old:
        want = round(s["size"] + plus, 3); alt = [e for e in size_ex.get(s["text"], []) if abs(e["old"] - s["size"]) < 0.01]
        cands = [k for k in pool if k[0] == s["text"] and k[2] == s["bold"] and pool[k] > 0 and (abs(k[1] - want) < 0.015 or any(abs(k[1] - e["new"]) < 0.015 for e in alt))]
        if not cands:
            if rem_left[s["text"]] > 0: rem_left[s["text"]] -= 1; continue
            bad.append((s["text"], s["size"], want)); continue
        k = sorted(cands, key=lambda k: abs(k[1] - want))[0]; pool[k] -= 1; matched += 1
        if abs(k[1] - want) >= 0.015: ex_used[(s["text"], s["size"], k[1])] += 1
    extra = +Counter({k: v for k, v in pool.items() if v > 0})
    extra_ok = all(exp_add[k[0]] >= v for k, v in extra.items()) and sum(extra.values()) == sum(exp_add.values())
    ok_sz = not bad and extra_ok
    lines.append(f"{'PASS' if ok_sz else 'FAIL'}  sizes: {matched} V26 spans matched by text and weight at size +{plus:g} (or a declared exception), {len(bad)} unmatched {bad[:8]}, {sum(extra.values())} spans beyond the matches {dict(extra) if extra else ''}")
    if ex_used: lines.append("INFO  size exceptions used (text, V26 size, NEW size): " + "; ".join(f"{t!r} {a:g} -> {b:g} x{n}" for (t, a, b), n in ex_used.items()))
    for e in ex.get("sizes", []): lines.append(f"INFO  declared size exception: {e['text']!r} {e['old']:g} -> {e['new']:g}: {e.get('why', '')}")
    for e in ex.get("removed", []): lines.append(f"INFO  declared removed: {e['text']!r} x{e.get('count', 1)}: {e.get('why', '')}")
    for e in ex.get("added", []): lines.append(f"INFO  declared added: {e['text']!r} x{e.get('count', 1)}: {e.get('why', '')}")
    so, sn = Counter(s["size"] for s in old), Counter(s["size"] for s in new)
    lines.append(f"INFO  size histogram V26 {dict(sorted(so.items()))} -> NEW {dict(sorted(sn.items()))}")
    fo, fn = sorted({s["font"].split("+")[-1] for s in old}), sorted({s["font"].split("+")[-1] for s in new})
    lines.append(f"INFO  fonts V26 {fo} -> NEW {fn}")
    ok = ok_str and ok_sz
    json.dump(dict(sheet=sheet, new=new_pdf, old=old_pdf, new_sha256=sha256(new_pdf), old_sha256=sha256(old_pdf), n_new=len(new), n_old=len(old), lost=dict(lost), gained=dict(gained), unmatched=bad, extra=[(list(k), v) for k, v in extra.items()], exceptions=ex, ok=ok, new_spans=new, old_spans=old), open(f"{out_dir}/census.json", "w"), indent=1, ensure_ascii=False)
    open(f"{out_dir}/census.txt", "w").write("\n".join(lines) + "\n")
    return ok, lines


def page_size_check(new_pdf, old_pdf):
    from pypdf import PdfReader
    a, b = PdfReader(hydrated(new_pdf)), PdfReader(hydrated(old_pdf)); assert len(a.pages) == 1 and len(b.pages) == 1
    ma, mb = a.pages[0].mediabox, b.pages[0].mediabox
    ok = abs(float(ma.width) - float(mb.width)) <= 0.05 and abs(float(ma.height) - float(mb.height)) <= 0.05
    return ok, (round(float(ma.width), 4), round(float(ma.height), 4)), (round(float(mb.width), 4), round(float(mb.height), 4))


def finish_checks(path, lines):
    body = [l for l in lines]; n_fail = sum(1 for l in body if l.startswith("FAIL")); n_pass = sum(1 for l in body if l.startswith("PASS"))
    res = "RESULT ALL PASS" if n_fail == 0 and n_pass > 0 else "RESULT FAIL"
    open(path, "w").write("\n".join(body + ["", f"{n_pass} PASS, {n_fail} FAIL", res]) + "\n"); print(res); return res
