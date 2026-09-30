"""Lane LSUPP-C (round 49): the per-sheet verifier. run(spec) proves, with Ghostscript only:
  1. page box (pypdf mediabox) of NEW equals V26's within 0.05 pt;
  2. Ghostscript txtwrite census NEW against V26: the multiset of joined baseline strings is identical, every (word, size) and
     every (char, size) of V26 reappears at size + 1.0 except the declared kept strings (listed with their reason), fonts Arial;
  3. every drawn string of the builder's record (values proved against the numbers files by the builder's own gates) is on the
     sheet, declared removals absent, declared wordings present;
  4. data marks unchanged: the pixels of the mark colours on 150 dpi renders WITHOUT anti-aliasing coincide with V26's outside
     the declared legend boxes; no colour introduced;
  5. 150 dpi deliverable render (Ghostscript), OLD vs NEW side by side, then the sheet-specific extra checks;
  writes <out_dir>/verify/checks.txt ending RESULT ALL PASS / RESULT FAIL and a verify record json."""
import json, os
from collections import Counter
import r49_common as C


def count_in_lines(t, lines):
    if " " not in t.strip(): return sum(1 for ln in lines for tok in ln.split() if tok == t)
    return sum(ln.count(t) for ln in lines)


def box_around(spans, pred, pad_left=0.0, pad_right=4.0, pad_y=6.0):
    ss = [s for s in spans if pred(s["text"])]
    if not ss: return None
    return [min(s["x"] for s in ss) - pad_left, min(s["y"] for s in ss) - pad_y - max(s["size"] for s in ss), max(s["x1"] for s in ss) + pad_right, max(s["y"] for s in ss) + pad_y]


def run(spec):
    S = spec["sheet"]; out = spec["out_dir"]; work, ver = f"{out}/work", f"{out}/verify"
    for d_ in (work, ver, f"{ver}/crops"): os.makedirs(d_, exist_ok=True)
    OLD, NEW = spec["old_pdf"], spec["new_pdf"]; C.hydrated(OLD); C.hydrated(NEW)
    ck = C.Checks(f"{S}, lane LSUPP-C round 49 (every text +1 pt), {C.now()}. OLD = V26 {OLD} (sha256 {C.sha256(OLD, 16)}); NEW = {NEW} (sha256 {C.sha256(NEW, 16)}). "
                  "Text layers by Ghostscript txtwrite (TextFormat 0), renders by Ghostscript png16m, one job at a time under the watchdog; PyMuPDF used for placement only (child process).")
    kept = spec.get("kept", []); ck.info("declared: " + spec.get("note", "") + (f"; kept at their current size: {[(t, s, m) for t, s, m in kept]}" if kept else "; no string keeps its size"))
    # 1 page box
    pb_new, pb_old = C.page_box(NEW), C.page_box(OLD)
    ck.log(f"page box of NEW equals V26 within 0.05 pt (pypdf mediabox)", abs(pb_new[0] - pb_old[0]) < 0.05 and abs(pb_new[1] - pb_old[1]) < 0.05, f"NEW {pb_new[0]:.3f} x {pb_new[1]:.3f}, V26 {pb_old[0]:.3f} x {pb_old[1]:.3f}")
    # 2 census
    sp_old = C.gs_text(OLD, f"{work}/{S}_V26_txtwrite.xml", f"{S}_txt_old"); sp_new = C.gs_text(NEW, f"{work}/{S}_NEW_txtwrite.xml", f"{S}_txt_new")
    json.dump({"old": [{k: v for k, v in s.items() if k != "chars"} for s in sp_old], "new": [{k: v for k, v in s.items() if k != "chars"} for s in sp_new]}, open(f"{work}/{S}_text_layers.json", "w"), indent=0)
    fonts_new = sorted({s["font"] for s in sp_new}); fonts_old = sorted({s["font"] for s in sp_old})
    ck.log("fonts of NEW within the Arial family", all("Arial" in f for f in fonts_new), f"NEW {fonts_new}; V26 {fonts_old}")
    ln_old, ln_new = C.string_multiset(sp_old), C.string_multiset(sp_new)
    ck.log("census: the multiset of strings of NEW equals V26's (spans of one baseline joined across kerning splits and spaces, split at gaps over 0.36 em, which only separate text artists leave; same strings, same counts)", ln_old == ln_new, f"V26 {sum(ln_old.values())} strings, NEW {sum(ln_new.values())}; {C.fmt_delta(ln_old, ln_new)}")
    kw, kc = C.kept_counters(kept)
    w_old, w_new = C.word_multiset(sp_old), C.word_multiset(sp_new); exp_w = C.expected_plus_one(w_old, kw)
    for key, n in kw.items(): assert w_old.get(key, 0) >= n, ("kept string not on V26 at that size", key, n, w_old.get(key, 0))
    ck.log("census: every (word, size) of V26 reappears on NEW at size + 1.0 pt, except the declared kept strings", w_new == exp_w, f"{sum(w_old.values())} words; {C.fmt_delta(exp_w, w_new)}")
    c_old, c_new = C.char_multiset(sp_old), C.char_multiset(sp_new); exp_c = C.expected_plus_one(c_old, kc)
    ck.log("census: every (character, size) of V26 reappears on NEW at size + 1.0 pt, except the declared kept strings", c_new == exp_c, f"{sum(c_old.values())} characters; {C.fmt_delta(exp_c, c_new)}")
    h_old = Counter(s["size"] for s in sp_old for _ in s["text"]); h_new = Counter(s["size"] for s in sp_new for _ in s["text"])
    ck.info(f"size histogram (characters per size): V26 {dict(sorted(h_old.items()))} -> NEW {dict(sorted(h_new.items()))}")
    ck.info(f"spans V26 {len(sp_old)}, NEW {len(sp_new)}; words V26 {sum(w_old.values())}, NEW {sum(w_new.values())}")
    alltext = " ".join(s["text"] for s in sp_new)
    ck.log("no em dash or semicolon on NEW", "—" not in alltext and ";" not in alltext, "")
    # 3 drawn strings, removals, wordings
    lines_new = [t for t, _c in C.lines_of(sp_new)]
    drawn = spec.get("drawn", Counter()); miss = [(t, n, count_in_lines(t, lines_new)) for t, n in drawn.items() if count_in_lines(t, lines_new) < n]
    ck.log(f"every drawn string of the builder's record ({sum(drawn.values())} strings, values proved against the numbers files by the builder's own gates, rc 0) is on NEW at least as often as recorded", not miss, str(miss[:10]))
    for t in spec.get("absent", []): ck.log(f"removed string absent from NEW: {t!r}", count_in_lines(t, lines_new) == 0, "")
    for t in spec.get("present", []): ck.log(f"wording present on NEW: {t!r}", count_in_lines(t, lines_new) >= 1, "")
    # 4 renders and marks
    R = {}
    R["new150"] = C.gs_render(NEW, f"{out}/{S}_150dpi.png", 150, f"{S}_new150")
    R["old150"] = C.gs_render(OLD, f"{work}/{S}_V26_150dpi.png", 150, f"{S}_old150")
    R["new_noaa"] = C.gs_render(NEW, f"{work}/{S}_NEW_150dpi_noAA.png", 150, f"{S}_new150noaa", aa=False)
    R["old_noaa"] = C.gs_render(OLD, f"{work}/{S}_V26_150dpi_noAA.png", 150, f"{S}_old150noaa", aa=False)
    ck.info("renders: " + ", ".join(os.path.basename(v) for v in R.values()) + " (Ghostscript png16m under the watchdog)")
    boxes = spec.get("legend_boxes", [])
    md = C.marks_unchanged(R["old_noaa"], R["new_noaa"], spec["mark_colours"], 150, boxes)
    tol = spec.get("marks_tolerance", 0)
    ck.log(f"data marks unchanged: the pixels of the mark colours {spec['mark_colours']} on the non-anti-aliased 150 dpi renders coincide with V26's outside the {len(boxes)} declared legend box(es)" + (f", at most {tol} px tolerated" if tol else ""),
           md["shape_equal"] and md["n_outside"] <= tol, f"mark pixels V26 {md.get('n_old')}, NEW {md.get('n_new')}, differing {md.get('n_diff')}, outside the boxes {md.get('n_outside')}, outside bbox {md.get('outside_bbox_pt')}")
    cen_old, cen_new = C.colour_census(R["old_noaa"]), C.colour_census(R["new_noaa"]); new_only = sorted(set(cen_new) - set(cen_old))
    ck.log("no colour introduced: every colour painted on NEW (non-anti-aliased census, at least 4 px) is on V26", not new_only, f"new-only {new_only[:10]}; NEW {len(cen_new)} colours, V26 {len(cen_old)}")
    json.dump({"old": cen_old, "new": cen_new, "new_only": new_only}, open(f"{ver}/{S}_colour_census.json", "w"), indent=1)
    C.side_by_side(R["old150"], R["new150"], f"{ver}/crops/{S}_V26_vs_NEW_150dpi.png")
    ck.info("crops: verify/crops/" + f"{S}_V26_vs_NEW_150dpi.png (V26 left, NEW right)")
    # 5 extras
    ctx = dict(sp_old=sp_old, sp_new=sp_new, lines_new=lines_new, R=R, cen_new=cen_new, cen_old=cen_old, pb=pb_new)
    if spec.get("extra"): spec["extra"](ck, ctx)
    for n in spec.get("notes", []): ck.info(n)
    ok = ck.write(f"{ver}/checks.txt")
    json.dump({"sheet": S, "old": OLD, "old_sha256": C.sha256(OLD), "new": NEW, "new_sha256": C.sha256(NEW), "page_box_new": pb_new, "page_box_old": pb_old, "kept": kept, "sources": spec.get("sources", {}), "written": C.now(), "result": "ALL PASS" if ok else "FAIL"},
              open(f"{ver}/{S}_verify_record.json", "w"), indent=1)
    print("RESULT ALL PASS" if ok else "RESULT FAIL"); return ok
