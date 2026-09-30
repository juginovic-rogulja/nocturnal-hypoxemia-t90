#!$T90_PY
"""Round 40, LMAIN: rasterize a composed vector sheet with the round-38 flat pipeline (ROUND38 scripts/flat/post.py idiom, copied here):
Ghostscript png16m at DPI under the lane watchdog (one render at a time, RSS cap 4 GB, time cap --max-s), then the PNG's own zlib stream
wrapped as a lossless one-page PDF at the EXACT vector page size (FlateDecode + PNG predictors, no re-encoding) and an LZW TIFF; size
check within 0.5 pt; then a 150 dpi gs render of the raster PDF (the deliverable's PNG). Writes <SD>/<Sheet>.pdf, <Sheet>.tif,
<Sheet>_150dpi.png, work/<Sheet>_<dpi>dpi.png, work/raster_record.json.
usage: rasterize_r40.py SHEET DPI [--max-s 3600] [--vector PATH]"""
import argparse, hashlib, json, os, struct, sys, time
from PIL import Image
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lmain_lib as L
Image.MAX_IMAGE_PIXELS = None


def fmt(v): return ("%.4f" % v).rstrip("0").rstrip(".")


def png_chunks(path):
    with open(path, "rb") as f:
        assert f.read(8) == b"\x89PNG\r\n\x1a\n", "not a PNG"
        while True:
            hdr = f.read(8)
            if len(hdr) < 8: break
            ln, typ = struct.unpack(">I4s", hdr); data = f.read(ln); f.read(4)
            yield typ, data
            if typ == b"IEND": break


def build_pdf(png, pdf, w_pt, h_pt):
    idat, ihdr = [], None
    for typ, data in png_chunks(png):
        if typ == b"IHDR": ihdr = struct.unpack(">IIBBBBB", data)
        elif typ == b"IDAT": idat.append(data)
    w, h, bitdepth, colortype, comp, filt, interlace = ihdr
    assert (bitdepth, colortype, comp, filt, interlace) == (8, 2, 0, 0, 0), ("unexpected PNG layout", ihdr)
    stream = b"".join(idat)
    content = ("q %s 0 0 %s 0 0 cm /Im0 Do Q" % (fmt(w_pt), fmt(h_pt))).encode()
    objs = [b"<< /Type /Catalog /Pages 2 0 R >>", b"<< /Type /Pages /Kids [3 0 R] /Count 1 >>",
            ("<< /Type /Page /Parent 2 0 R /MediaBox [0 0 %s %s] /Resources << /XObject << /Im0 5 0 R >> >> /Contents 4 0 R >>" % (fmt(w_pt), fmt(h_pt))).encode(),
            b"<< /Length %d >>\nstream\n" % len(content) + content + b"\nendstream",
            ("<< /Type /XObject /Subtype /Image /Width %d /Height %d /ColorSpace /DeviceRGB /BitsPerComponent 8 /Filter /FlateDecode /DecodeParms << /Predictor 15 /Colors 3 /BitsPerComponent 8 /Columns %d >> /Length %d >>\nstream\n" % (w, h, w, len(stream))).encode() + stream + b"\nendstream",
            b"<< /Producer (T90 ROUND40 LMAIN raster wrapper, Ghostscript png16m render, round-38 flat pipeline) >>"]
    out = bytearray(b"%PDF-1.4\n%\xe2\xe3\xcf\xd3\n"); offsets = []
    for i, o in enumerate(objs, 1):
        offsets.append(len(out)); out += b"%d 0 obj\n" % i + o + b"\nendobj\n"
    xref = len(out); out += b"xref\n0 %d\n0000000000 65535 f \n" % (len(objs) + 1)
    for off in offsets: out += b"%010d 00000 n \n" % off
    out += b"trailer\n<< /Size %d /Root 1 0 R /Info 6 0 R >>\nstartxref\n%d\n%%%%EOF\n" % (len(objs) + 1, xref)
    with open(pdf, "wb") as f: f.write(out)
    return w, h, len(stream)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("sheet"); ap.add_argument("dpi", type=int); ap.add_argument("--max-s", type=int, default=3600); ap.add_argument("--vector", default=None); ap.add_argument("--tag", default=None, help="write work/<sheet>_<tag>.pdf/.tif/_150dpi.png instead of the deliverable (a twin for like-for-like diffs)"); a = ap.parse_args()
    sheet, dpi = a.sheet, a.dpi; SD = f"{L.LANE}/{sheet}"; WORK = f"{SD}/work"; vec = L.hydrated(a.vector or f"{SD}/{sheet}_vector.pdf")
    n, w_pt, h_pt = L.gs_pdfinfo(vec); assert n == 1, n
    png = f"{WORK}/{sheet}_{dpi}dpi{'_' + a.tag if a.tag else ''}.png"; t0 = time.time()
    st = L.gs_render(vec, png, dpi, f"{sheet}_gs{dpi}_raster{'_' + a.tag if a.tag else ''}", max_s=a.max_s); t_gs = round(time.time() - t0, 1)
    stem = f"{WORK}/{sheet}_{a.tag}" if a.tag else f"{SD}/{sheet}"; pdf = f"{stem}.pdf"; tif = f"{stem}.tif"
    w, h, nbytes = build_pdf(png, pdf, w_pt, h_pt)
    im = Image.open(png); assert im.size == (w, h) and im.mode == "RGB", (im.size, im.mode)
    im.save(tif, "TIFF", compression="tiff_lzw", dpi=(dpi, dpi)); im2 = Image.open(tif); im2.load(); im2_size = im2.size; im2 = None; im = None
    w_calc, h_calc = w / dpi * 72, h / dpi * 72
    n2, w2, h2 = L.gs_pdfinfo(pdf); assert n2 == 1 and abs(w2 - w_pt) < 0.01 and abs(h2 - h_pt) < 0.01, (n2, w2, h2, w_pt, h_pt)
    png150 = f"{stem}_150dpi.png"; st150 = L.gs_render(pdf, png150, 150, f"{sheet}_gs150_of_raster{'_' + a.tag if a.tag else ''}", max_s=900)
    res = {"sheet": sheet, "dpi": dpi, "vector": vec, "vector_sha256": L.sha256(vec), "gs_status": st, "gs_wall_s": t_gs, "png": png, "png_bytes": os.path.getsize(png), "width_px": w, "height_px": h,
           "orig_w_pt": w_pt, "orig_h_pt": h_pt, "calc_w_pt": round(w_calc, 3), "calc_h_pt": round(h_calc, 3), "dw_pt": round(w_calc - w_pt, 3), "dh_pt": round(h_calc - h_pt, 3),
           "size_check_within_0.5pt": abs(w_calc - w_pt) <= 0.5 and abs(h_calc - h_pt) <= 0.5, "pdf": pdf, "pdf_bytes": os.path.getsize(pdf), "pdf_mediabox_gs": [w2, h2], "pdf_image_stream_bytes": nbytes,
           "tif": tif, "tif_bytes": os.path.getsize(tif), "tif_reopen_size": list(im2_size), "png150": png150, "png150_status": st150,
           "sha256_png": L.sha256(png), "sha256_pdf": L.sha256(pdf), "sha256_tif": L.sha256(tif), "sha256_png150": L.sha256(png150), "written": L.now()}
    json.dump(res, open(f"{WORK}/raster_record{'_' + a.tag if a.tag else ''}.json", "w"), indent=1)
    print(json.dumps({k: res[k] for k in ("sheet", "dpi", "gs_wall_s", "width_px", "height_px", "orig_w_pt", "orig_h_pt", "dw_pt", "dh_pt", "size_check_within_0.5pt", "pdf_bytes", "tif_bytes", "sha256_pdf")}, indent=1))
    assert res["size_check_within_0.5pt"]


if __name__ == "__main__": main()
