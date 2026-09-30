#!/usr/bin/env python3
"""Round 50 finish step (LSUPP-A), per sheet: (1) page box equals V27's within 0.05 pt, (2) Ghostscript text census against V27: the same
multiset of strings except the declared removed strings, every size unchanged, same font faces, (3) the lane's own verify checks all PASS
(read from verify/checks.txt as the builder chain left it), (4) the 150 dpi PNG rendered by Ghostscript, (5) no em dash or semicolon.
Usage: r50_finish.py Supp_FigNN [--removed "string" ...] [--eye "note"]"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../../../../..")))  # repository root, where paths.py lives
import paths
import os, sys, json, subprocess, hashlib, time
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import census as CEN
from pypdf import PdfReader
T90 = paths.FIGURE_ROOT
LANE = f"{paths.FIGURE_ROOT}/ROUND50_2026-09-26/lanes/LSUPP-A"; BASE_SET = f"{paths.FIGURE_ROOT}/ROUND49_2026-09-26/figures/NEW_FINAL_SET_V27"; GS = paths.GS

def main():
    args = sys.argv[1:]; S = args[0]; removed = []; eye = ""
    while "--removed" in args:
        i = args.index("--removed"); removed.append(args[i + 1]); del args[i:i + 2]
    if "--eye" in args:
        i = args.index("--eye"); eye = args[i + 1]; del args[i:i + 2]
    SD = f"{LANE}/{S}"; VER = f"{SD}/verify"; new_pdf = f"{SD}/{S}.pdf"; base_pdf = f"{BASE_SET}/{S}.pdf"; os.makedirs(VER, exist_ok=True)
    lines = []
    def check(ok, name, detail=""):
        lines.append(f"{'PASS' if ok else 'FAIL'} {name}  {detail}"); print(lines[-1][:300]); return ok
    cp = f"{VER}/checks.txt"; prior = []
    if os.path.exists(cp): prior = [l.rstrip("\n") for l in open(cp) if l.strip() and not l.startswith("RESULT")]
    if "--- round 50 finish ---" in prior: prior = prior[:prior.index("--- round 50 finish ---")]
    nb = PdfReader(new_pdf).pages[0].mediabox; ob = PdfReader(base_pdf).pages[0].mediabox
    nw, nh, ow, oh = float(nb.width), float(nb.height), float(ob.width), float(ob.height)
    check(abs(nw - ow) <= 0.05 and abs(nh - oh) <= 0.05 and len(PdfReader(new_pdf).pages) == 1, "page box equals V27 within 0.05 pt (one page)", f"new {nw:.2f} x {nh:.2f}, V27 {ow:.2f} x {oh:.2f}")
    out = f"{VER}/census"; os.makedirs(out, exist_ok=True)
    nsp = CEN.spans(CEN.txtwrite(new_pdf, f"{out}/new_txtwrite.xml")); osp = CEN.spans(CEN.txtwrite(base_pdf, f"{out}/v27_txtwrite.xml"))
    o_strings = Counter(s for s, _, _, _ in osp); n_strings = Counter(s for s, _, _, _ in nsp); rem = Counter(removed)
    expected = o_strings - rem
    missing_removed = sorted((rem - o_strings).elements())
    lost = sorted((expected - n_strings).elements()); gained = sorted((n_strings - expected).elements())
    check(n_strings == expected and not missing_removed, f"census: the same multiset of strings as V27 except the {sum(rem.values())} declared removed strings ({len(osp)} spans on V27, {len(nsp)} new)", f"lost {lost[:8]} gained {gained[:8]} removed-not-on-V27 {missing_removed[:8]}")
    check(all(s not in n_strings for s in rem), "census: every declared removed string is absent from the new sheet", str(sorted(s for s in rem if s in n_strings)))
    o_pairs = Counter((s, size) for s, size, _, _ in osp); n_pairs = Counter((s, size) for s, size, _, _ in nsp)
    exp_pairs = Counter(); left = Counter(rem)
    for (s, size), k in sorted(o_pairs.items()):
        drop = min(k, left[s]); left[s] -= drop
        if k - drop: exp_pairs[(s, size)] += k - drop
    sizes_o = sorted({size for _, size, _, _ in osp}); sizes_n = sorted({size for _, size, _, _ in nsp})
    check(exp_pairs == n_pairs, f"census: every remaining string at its V27 size, unchanged (V27 sizes {sizes_o}, new {sizes_n})", f"missing {sorted((exp_pairs - n_pairs).elements())[:8]} unexpected {sorted((n_pairs - exp_pairs).elements())[:8]}")
    fo = {f for _, _, f, _ in osp}; fn = {f for _, _, f, _ in nsp}
    check(fn <= fo, "census: font faces within V27's", f"V27 {sorted(fo)} new {sorted(fn)}")
    prior_checks = [l for l in prior if l.startswith(("PASS", "FAIL"))]
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
    body = prior + ["--- round 50 finish ---"] + lines + [f"provenance: {S}.pdf sha256 {sha}, V27 {base_pdf}, finished {time.strftime('%Y-%m-%d %H:%M:%S')}"]
    body.append("RESULT ALL PASS" if ok else "RESULT FAIL: " + "; ".join(l.split("  ")[0][5:] for l in prior_checks + lines if l.startswith("FAIL")))
    open(cp, "w").write("\n".join(body) + "\n"); print(body[-1])
    json.dump(dict(sheet=S, new_pdf=new_pdf, sha256=sha, v27=base_pdf, page_new=[nw, nh], page_v27=[ow, oh], removed=removed, sizes_v27=sizes_o, sizes_new=sizes_n, n_spans=[len(osp), len(nsp)], eye=eye), open(f"{VER}/r50_record.json", "w"), indent=1)
    sys.exit(0 if ok else 1)

if __name__ == "__main__": main()
