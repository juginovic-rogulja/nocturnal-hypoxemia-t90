#!$T90_PY
"""Round 40, LMAIN: gather the deliverables' facts for LANE_REPORT.md (sha256, page boxes, dpi, checks results, OCR summary lines,
CHANGES rows) into verify/LMAIN_DELIVERABLES.json and print a markdown table. No PDF is opened here (sha256 and the records only)."""
import csv, json, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lmain_lib as L

SHEETS = ["Main_Fig1", "Main_Fig2", "Main_Fig3", "Main_Fig4", "Main_Fig5", "Main_Fig6"]
out = {}
for s in SHEETS:
    SD = f"{L.LANE}/{s}"; d = {"sheet": s, "files": {}}
    for fn in (f"{s}.pdf", f"{s}.tif", f"{s}_vector.pdf", f"{s}_150dpi.png"):
        p = f"{SD}/{fn}"
        if os.path.exists(p): d["files"][fn] = {"sha256": L.sha256(p), "bytes": os.path.getsize(p)}
    if os.path.exists(f"{SD}/{s}.pdf"): d["page_gs"] = L.gs_pdfinfo(f"{SD}/{s}.pdf")
    rr = f"{SD}/work/raster_record.json"
    if os.path.exists(rr):
        r = json.load(open(rr)); d["raster"] = {k: r[k] for k in ("dpi", "width_px", "height_px", "orig_w_pt", "orig_h_pt", "dw_pt", "dh_pt", "gs_wall_s")}
    ck = f"{SD}/verify/checks.txt"
    if os.path.exists(ck):
        lines = open(ck).read().splitlines(); d["result"] = lines[-1]; d["n_pass"] = sum(1 for l in lines if l.startswith("PASS")); d["n_fail"] = sum(1 for l in lines if l.startswith("FAIL"))
        d["ocr"] = [l for l in lines if "OCR" in l][:6]; d["fails"] = [l[:300] for l in lines if l.startswith("FAIL")]
    ch = f"{SD}/verify/CHANGES_{s}.csv"
    if os.path.exists(ch): d["changes"] = list(csv.DictReader(open(ch)))
    out[s] = d
json.dump(out, open(f"{L.LANE}/verify/LMAIN_DELIVERABLES.json", "w"), indent=1)
print("| sheet | deliverable | page (pt, gs) | dpi / px | sha256 (deliverable) | result |"); print("|---|---|---|---|---|---|")
for s in SHEETS:
    d = out[s]; f = d["files"].get(f"{s}.pdf", {}); r = d.get("raster"); pg = d.get("page_gs")
    print(f"| {s} | {'raster PDF + TIFF' if r else 'vector PDF'} | {pg[1] if pg else ''} x {pg[2] if pg else ''} | {(str(r['dpi']) + ' / ' + str(r['width_px']) + ' x ' + str(r['height_px'])) if r else 'vector'} | {f.get('sha256', '')[:16]} | {d.get('result', '')} ({d.get('n_pass', 0)} pass, {d.get('n_fail', 0)} fail) |")
for s in SHEETS:
    for l in out[s].get("ocr", []): print(s, "|", l[:400])
    for l in out[s].get("fails", []): print(s, "| FAIL", l[:300])
