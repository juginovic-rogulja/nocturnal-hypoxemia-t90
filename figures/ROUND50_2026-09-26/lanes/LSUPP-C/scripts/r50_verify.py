"""Lane LSUPP-C (round 50): the per-sheet verifier against V27. run(spec) proves, with Ghostscript only:
  1. page box (pypdf mediabox) of NEW equals the expected box (V27's unless the spec declares a new one) within 0.05 pt;
  2. Ghostscript txtwrite census against V27 with SIZES UNCHANGED: the multiset of strings (spans of one baseline joined across
     kerning splits, split at gaps over 0.36 em) is identical except the declared removed and added strings, every (word, size) and
     every (character, size) identical except those strings, fonts Arial, no em dash or semicolon;
  3. every drawn string of the builder's record on the sheet, declared removals absent, declared wordings present;
  4. declared regions of the anti-aliased 150 dpi renders pixel-identical to V27 (threshold 32 of 255) outside their allowed boxes;
  5. data marks unchanged where the page size is unchanged (mark-colour pixels of the non-anti-aliased renders), no colour introduced;
  6. the 150 dpi deliverable render, V27 vs NEW side by side, the sheet's extra checks; checks.txt ends RESULT ALL PASS / RESULT FAIL."""
import json, os
from collections import Counter
import r49_common as C


def count_in_lines(t, lines):
    if " " not in t.strip(): return sum(1 for ln in lines for tok in ln.split() if tok == t)
    return sum(ln.count(t) for ln in lines)


def _wc(strings):
    """(word, size) and (char, size) Counters of [(string, size)]."""
    w, c = Counter(), Counter()
    for text, size in strings:
        for tok in text.split(): w[(tok, size)] += 1
        for ch in text:
            if not ch.isspace(): c[(ch, size)] += 1
    return w, c


def run(spec):
    S = spec["sheet"]; out = spec["out_dir"]; work, ver = f"{out}/work", f"{out}/verify"
    for d_ in (work, ver, f"{ver}/crops"): os.makedirs(d_, exist_ok=True)
    OLD, NEW = spec["old_pdf"], spec["new_pdf"]; C.hydrated(OLD); C.hydrated(NEW)
    ck = C.Checks(f"{S}, lane LSUPP-C round 50 (Alen's comment, sizes unchanged), {C.now()}. OLD = V27 {OLD} (sha256 {C.sha256(OLD, 16)}); NEW = {NEW} (sha256 {C.sha256(NEW, 16)}). "
                  "Text layers by Ghostscript txtwrite (TextFormat 0), renders by Ghostscript png16m, one job at a time under the watchdog; PyMuPDF for placement only (child process).")
    removed, added = list(spec.get("removed", [])), list(spec.get("added", []))
    ck.info("declared: " + spec.get("note", "") + f"; removed strings {[t for t, _s in removed]}; added strings {[t for t, _s in added]}")
    # 1 page box
    pb_new, pb_old = C.page_box(NEW), C.page_box(OLD); exp = spec.get("page_box") or list(pb_old)
    ck.log(f"page box of NEW = {exp[0]:.3f} x {exp[1]:.3f} pt ({spec.get('page_box_note') or 'as V27'})", abs(pb_new[0] - exp[0]) < 0.05 and abs(pb_new[1] - exp[1]) < 0.05, f"NEW {pb_new[0]:.3f} x {pb_new[1]:.3f}, V27 {pb_old[0]:.3f} x {pb_old[1]:.3f}")
    # 2 census
    sp_old = C.gs_text(OLD, f"{work}/{S}_V27_txtwrite.xml", f"{S}_txt_old"); sp_new = C.gs_text(NEW, f"{work}/{S}_NEW_txtwrite.xml", f"{S}_txt_new")
    json.dump({"old": [{k: v for k, v in s.items() if k != "chars"} for s in sp_old], "new": [{k: v for k, v in s.items() if k != "chars"} for s in sp_new]}, open(f"{work}/{S}_text_layers.json", "w"), indent=0)
    fonts_new = sorted({s["font"] for s in sp_new}); fonts_old = sorted({s["font"] for s in sp_old})
    ck.log("fonts of NEW within the Arial family", all("Arial" in f for f in fonts_new), f"NEW {fonts_new}; V27 {fonts_old}")
    st_old, st_new = C.string_multiset(sp_old), C.string_multiset(sp_new)
    exp_st = st_old - Counter(t for t, _s in removed) + Counter(t for t, _s in added)
    for t, _s in removed: assert st_old.get(t, 0) >= 1, ("declared removed string not on V27", t)
    ck.log("census: the multiset of strings of NEW equals V27's except the declared removed and added strings", st_new == exp_st, f"V27 {sum(st_old.values())} strings, NEW {sum(st_new.values())}; {C.fmt_delta(exp_st, st_new)}")
    rw, rc_ = _wc(removed); aw, ac = _wc(added)
    w_old, w_new = C.word_multiset(sp_old), C.word_multiset(sp_new); exp_w = w_old - rw + aw
    ck.log("census: every (word, size) of V27 reappears on NEW at the SAME size, except the declared strings (sizes unchanged this round)", w_new == exp_w, f"{sum(w_old.values())} words; {C.fmt_delta(exp_w, w_new)}")
    c_old, c_new = C.char_multiset(sp_old), C.char_multiset(sp_new); exp_c = c_old - rc_ + ac
    ck.log("census: every (character, size) of V27 reappears on NEW at the same size, except the declared strings", c_new == exp_c, f"{sum(c_old.values())} characters; {C.fmt_delta(exp_c, c_new)}")
    h_old = Counter(s["size"] for s in sp_old for _ in s["text"]); h_new = Counter(s["size"] for s in sp_new for _ in s["text"])
    ck.info(f"size histogram (characters per size): V27 {dict(sorted(h_old.items()))} -> NEW {dict(sorted(h_new.items()))}")
    alltext = " ".join(s["text"] for s in sp_new)
    ck.log("no em dash or semicolon on NEW", "\u2014" not in alltext and ";" not in alltext, "")
    # 3 drawn, absent, present
    lines_new = [t for t, _c in C.lines_of(sp_new)]; strings_new = C.strings_of(sp_new)
    drawn = spec.get("drawn", Counter()); miss = [(t, n, count_in_lines(t, lines_new)) for t, n in drawn.items() if count_in_lines(t, lines_new) < n]
    ck.log(f"every drawn string of the builder's record ({sum(drawn.values())} strings, values proved against the numbers files by the builder's own gates, rc 0) is on NEW at least as often as recorded", not miss, str(miss[:10]))
    for t, _s in removed: ck.log(f"removed string no longer drawn as a string on NEW: {t!r}", t not in strings_new, "" if t not in strings_new else "still a string of its own")
    for t in spec.get("absent", []): ck.log(f"removed wording absent from NEW: {t!r}", count_in_lines(t, lines_new) == 0 and t not in strings_new, "")
    for t in [t for t, _s in added] + list(spec.get("present", [])): ck.log(f"string present on NEW: {t!r}", count_in_lines(t, lines_new) >= 1 or t in strings_new, "")
    # 4 renders and regions
    R = {}
    R["new150"] = C.gs_render(NEW, f"{out}/{S}_150dpi.png", 150, f"{S}_new150")
    R["old150"] = C.gs_render(OLD, f"{work}/{S}_V27_150dpi.png", 150, f"{S}_old150")
    R["new_noaa"] = C.gs_render(NEW, f"{work}/{S}_NEW_150dpi_noAA.png", 150, f"{S}_new150noaa", aa=False)
    R["old_noaa"] = C.gs_render(OLD, f"{work}/{S}_V27_150dpi_noAA.png", 150, f"{S}_old150noaa", aa=False)
    ck.info("renders: " + ", ".join(os.path.basename(v) for v in R.values()) + " (Ghostscript png16m under the watchdog)")
    for rg in spec.get("regions", []):
        rd = C.region_diff(R["old150"], R["new150"], 150, rg["box"], rg.get("allowed", []), thresh=rg.get("thresh", 32))
        ck.log(f"region '{rg['name']}' {rg['box']}: NEW pixel-identical to V27 outside {len(rg.get('allowed', []))} allowed box(es) at threshold {rg.get('thresh', 32)} of 255 ({rg.get('note', '')})", rd["n_outside"] <= rg.get("max_outside", 0),
               f"{rd['n_diff']} differing px of {rd['size'][0]}x{rd['size'][1]}, {rd['n_outside']} outside, outside bbox {rd['outside_bbox_pt']}")
    if spec.get("mark_colours"):
        md = C.marks_unchanged(R["old_noaa"], R["new_noaa"], spec["mark_colours"], 150, spec.get("legend_boxes", []))
        ck.log(f"data marks unchanged: the pixels of the mark colours {spec['mark_colours']} on the non-anti-aliased renders coincide with V27's outside {len(spec.get('legend_boxes', []))} allowed box(es)", md["shape_equal"] and md["n_outside"] <= spec.get("marks_tolerance", 0), f"mark pixels V27 {md.get('n_old')}, NEW {md.get('n_new')}, differing {md.get('n_diff')}, outside {md.get('n_outside')}, bbox {md.get('outside_bbox_pt')}")
    cen_old, cen_new = C.colour_census(R["old_noaa"]), C.colour_census(R["new_noaa"]); new_only = sorted(set(cen_new) - set(cen_old))
    ck.log("no colour introduced: every colour painted on NEW (non-anti-aliased census, at least 4 px) is on V27", not new_only, f"new-only {new_only[:10]}; NEW {len(cen_new)} colours, V27 {len(cen_old)}; gone {sorted(set(cen_old) - set(cen_new))}")
    json.dump({"old": cen_old, "new": cen_new, "new_only": new_only}, open(f"{ver}/{S}_colour_census.json", "w"), indent=1)
    C.side_by_side(R["old150"], R["new150"], f"{ver}/crops/{S}_V27_vs_NEW_150dpi.png")
    ck.info(f"crops: verify/crops/{S}_V27_vs_NEW_150dpi.png (V27 left, NEW right)")
    ctx = dict(sp_old=sp_old, sp_new=sp_new, lines_new=lines_new, R=R, cen_new=cen_new, cen_old=cen_old, pb=pb_new)
    if spec.get("extra"): spec["extra"](ck, ctx)
    for n in spec.get("notes", []): ck.info(n)
    ok = ck.write(f"{ver}/checks.txt")
    json.dump({"sheet": S, "old": OLD, "old_sha256": C.sha256(OLD), "new": NEW, "new_sha256": C.sha256(NEW), "page_box_new": pb_new, "page_box_old": pb_old, "removed": removed, "added": added, "sources": spec.get("sources", {}), "written": C.now(), "result": "ALL PASS" if ok else "FAIL"},
              open(f"{ver}/{S}_verify_record.json", "w"), indent=1)
    print("RESULT ALL PASS" if ok else "RESULT FAIL"); return ok
