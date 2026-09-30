#!/usr/bin/env python3
"""Round 49 finish step (LSUPP-A), per sheet: (1) page box equals V26's within 0.05 pt, (2) Ghostscript text census against V26 (same strings,
every size +1.0, declared exceptions), (3) the lane's own 03_verify checks all PASS (read from verify/checks.txt as the builder chain left it),
(4) the 150 dpi PNG rendered by Ghostscript. Appends the round-49 lines to verify/checks.txt and ends it with RESULT ALL PASS or RESULT FAIL.
Usage: r49_finish.py Supp_FigNN [--except "string=old_size" ...] [--eye "note"]"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import os, sys, json, subprocess, hashlib, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import census as CEN
from pypdf import PdfReader
T90 = paths.FIGURE_ROOT
LANE = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/lanes/LSUPP-A"; V26 = f"{paths.FIGURE_ROOT}/ROUND47_2026-09-24/figures/NEW_FINAL_SET_V26"; GS = paths.GS

def main():
    args = sys.argv[1:]; S = args[0]; exc = {}; eye = ""
    while "--except" in args:
        i = args.index("--except"); k, v = args[i + 1].rsplit("=", 1); exc[k] = float(v); del args[i:i + 2]
    if "--eye" in args:
        i = args.index("--eye"); eye = args[i + 1]; del args[i:i + 2]
    SD = f"{LANE}/{S}"; VER = f"{SD}/verify"; new_pdf = f"{SD}/{S}.pdf"; v26_pdf = f"{V26}/{S}.pdf"; os.makedirs(VER, exist_ok=True)
    lines = []
    def check(ok, name, detail=""):
        lines.append(f"{'PASS' if ok else 'FAIL'} {name}  {detail}"); print(lines[-1][:300]); return ok
    # the builder chain's own checks (written before this step), kept verbatim, must all PASS
    prior = []
    cp = f"{VER}/checks.txt"
    if os.path.exists(cp):
        prior = [l.rstrip("\n") for l in open(cp) if l.strip() and not l.startswith("RESULT") and not l.startswith("=== round 49")]
    # strip an earlier round-49 block (rerun)
    if "--- round 49 finish ---" in prior: prior = prior[:prior.index("--- round 49 finish ---")]
    nb = PdfReader(new_pdf).pages[0].mediabox; ob = PdfReader(v26_pdf).pages[0].mediabox
    nw, nh, ow, oh = float(nb.width), float(nb.height), float(ob.width), float(ob.height)
    check(abs(nw - ow) <= 0.05 and abs(nh - oh) <= 0.05 and len(PdfReader(new_pdf).pages) == 1, "page box equals V26 within 0.05 pt (one page)", f"new {nw:.2f} x {nh:.2f}, V26 {ow:.2f} x {oh:.2f}")
    out = f"{VER}/census"; os.makedirs(out, exist_ok=True)
    nsp = CEN.spans(CEN.txtwrite(new_pdf, f"{out}/new_txtwrite.xml")); osp = CEN.spans(CEN.txtwrite(v26_pdf, f"{out}/v26_txtwrite.xml"))
    from collections import Counter
    n_strings = Counter(s for s, _, _, _ in nsp); o_strings = Counter(s for s, _, _, _ in osp)
    lost = sorted((o_strings - n_strings).elements()); gained = sorted((n_strings - o_strings).elements())
    same = n_strings == o_strings
    if not same:   # a different span split (Ghostscript joins by spacing) is not a text change: compare per line
        same = CEN.lines(nsp) == CEN.lines(osp)
    check(same, f"census: the same multiset of strings as V26 ({len(osp)} spans on V26, {len(nsp)} new)", f"lost {lost[:8]} gained {gained[:8]}")
    o_pairs = Counter((s, size) for s, size, _, _ in osp); n_pairs = Counter((s, size) for s, size, _, _ in nsp); expected = Counter()
    for (s, size), k in o_pairs.items():
        expected[(s, size if (s in exc and abs(exc[s] - size) < 0.01) else round(size + 1.0, 2))] += k
    size_ok = expected == n_pairs
    if not size_ok and same:   # per-line sizes when the span split differs
        def line_sizes(sp):
            from collections import defaultdict
            rows = defaultdict(list)
            for s, size, f, bbox in sp:
                x0, y0, x1, y1 = [float(v) for v in bbox.split()]; rows[(round(y0), size)].append((x0, s))
            return Counter((" ".join(t for _, t in sorted(v)), k[1]) for k, v in rows.items())
        lo, ln = line_sizes(osp), line_sizes(nsp); exp2 = Counter()
        for (s, size), k in lo.items(): exp2[(s, size if (s in exc and abs(exc[s] - size) < 0.01) else round(size + 1.0, 2))] += k
        size_ok = exp2 == ln; miss = sorted((exp2 - ln).elements()); extra = sorted((ln - exp2).elements())
    else:
        miss = sorted((expected - n_pairs).elements()); extra = sorted((n_pairs - expected).elements())
    sizes_o = sorted({size for _, size, _, _ in osp}); sizes_n = sorted({size for _, size, _, _ in nsp})
    check(size_ok, f"census: every text size exactly +1.0 pt over V26 (V26 sizes {sizes_o} -> new {sizes_n}), exceptions {exc if exc else 'none'}", f"missing {miss[:8]} unexpected {extra[:8]}")
    fo = {f for _, _, f, _ in osp}; fn = {f for _, _, f, _ in nsp}
    check(fo == fn, "census: the same font faces (Arial, Arial Bold, Arial Italic subsets)", f"V26 {sorted(fo)} new {sorted(fn)}")
    prior_checks = [l for l in prior if l.startswith(("PASS", "FAIL"))]   # the L4 Checks format carries a title line and INFO lines beside the PASS/FAIL lines
    check(bool(prior_checks) and not any(l.startswith("FAIL") for l in prior_checks), f"the lane's own verify: {len(prior_checks)} checks all PASS (listed above)", "" if prior_checks else "no prior checks")
    png = f"{SD}/{S}_150dpi.png"
    r = subprocess.run([GS, "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=png16m", "-r150", "-dTextAlphaBits=4", "-dGraphicsAlphaBits=4", f"-sOutputFile={png}", new_pdf], capture_output=True, text=True)
    from PIL import Image
    im = Image.open(png) if r.returncode == 0 else None
    check(r.returncode == 0 and im is not None and abs(im.width - nw * 150 / 72) <= 1.5 and abs(im.height - nh * 150 / 72) <= 1.5, "150 dpi render by Ghostscript", f"{png} {im.size if im else r.stderr[-200:]}")
    if eye: check(True, "eye check on the 150 dpi render", eye)
    banned = [s for s, _, _, _ in nsp if ("—" in s or ";" in s)]
    check(not banned, "no em dash and no semicolon in any string", str(banned[:5]))
    sha = hashlib.sha256(open(new_pdf, "rb").read()).hexdigest()
    ok = bool(prior_checks) and not any(l.startswith("FAIL") for l in prior_checks + lines)
    body = prior + ["--- round 49 finish ---"] + lines + [f"provenance: {S}.pdf sha256 {sha}, V26 {v26_pdf}, finished {time.strftime('%Y-%m-%d %H:%M:%S')}"]
    body.append("RESULT ALL PASS" if ok else "RESULT FAIL: " + "; ".join(l.split("  ")[0][5:] for l in prior_checks + lines if l.startswith("FAIL")))
    open(cp, "w").write("\n".join(body) + "\n"); print(body[-1])
    json.dump(dict(sheet=S, new_pdf=new_pdf, sha256=sha, v26=v26_pdf, page_new=[nw, nh], page_v26=[ow, oh], exceptions=exc, sizes_v26=sizes_o, sizes_new=sizes_n, n_spans=[len(osp), len(nsp)], eye=eye), open(f"{VER}/r49_record.json", "w"), indent=1)
    sys.exit(0 if ok else 1)

if __name__ == "__main__": main()
