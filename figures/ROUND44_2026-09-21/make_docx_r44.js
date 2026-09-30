// Round 43: T90 All Figures docx, one figure per A4 page, label above. Placement rule (Alen, 2026-09-20): every sheet at PRINT SIZE
// (a 150-dpi render shown at 96/150 of its pixel size), reduced only where it would not fit the page box. usage: node make_docx_r43.js PNG_DIR OUT_DOCX
const fs = require("fs"), path = require("path");
const { Document, Packer, Paragraph, TextRun, ImageRun, PageBreak, AlignmentType } = require("docx");
const [PNG_DIR, OUT] = process.argv.slice(2);
const order = [];
for (let i = 1; i <= 6; i++) order.push([`Main_Fig${i}`, `Figure ${i}`]);
for (let i = 1; i <= 10; i++) order.push([`ED_Fig${String(i).padStart(2, "0")}`, `Extended Data Fig. ${i}`]);
for (let i = 1; i <= 21; i++) order.push([`Supp_Fig${String(i).padStart(2, "0")}`, `Supplementary Fig. ${i}`]);
const PAGE_W = 11906, PAGE_H = 16838, MARGIN = 1134;
const BOX_W = Math.floor((PAGE_W - 2 * MARGIN) / 1440 * 96);
const BOX_H = Math.floor((PAGE_H - 2 * MARGIN) / 1440 * 96) - 40;
const PRINT = 96 / 150;   // a 150-dpi render at print size on a 96-dpi page
function pngSize(buf) { return { w: buf.readUInt32BE(16), h: buf.readUInt32BE(20) }; }
const bufs = {}, sizes = {};
order.forEach(([name]) => { bufs[name] = fs.readFileSync(path.join(PNG_DIR, `${name}.png`)); sizes[name] = pngSize(bufs[name]); });
const children = []; const rec = { rule: "print size, reduced to fit the box", print_scale: PRINT, box_px: [BOX_W, BOX_H], sheets: {} };
order.forEach(([name, label], idx) => {
  const { w, h } = sizes[name]; const s = Math.min(PRINT, BOX_W / w, BOX_H / h);
  rec.sheets[name] = { scale: s, at_print_size: Math.abs(s - PRINT) < 1e-9, px: [Math.round(w * s), Math.round(h * s)], width_pct: Math.round(100 * w * s / BOX_W) };
  const paras = [
    new Paragraph({ spacing: { after: 120 }, children: [new TextRun({ text: label, bold: true, size: 22 })] }),
    new Paragraph({ alignment: AlignmentType.CENTER, children: [new ImageRun({
      type: "png", data: bufs[name], transformation: { width: Math.round(w * s), height: Math.round(h * s) },
      altText: { title: label, description: `${label} of the T90 manuscript`, name: name } })] }),
  ];
  if (idx < order.length - 1) paras.push(new Paragraph({ children: [new PageBreak()] }));
  children.push(...paras);
});
const doc = new Document({
  styles: { default: { document: { run: { font: "Arial", size: 22 } } } },
  sections: [{ properties: { page: { size: { width: PAGE_W, height: PAGE_H }, margin: { top: MARGIN, right: MARGIN, bottom: MARGIN, left: MARGIN } } }, children }],
});
Packer.toBuffer(doc).then(b => { fs.writeFileSync(OUT, b); fs.writeFileSync(OUT + ".scale.json", JSON.stringify(rec, null, 1)); console.log("wrote", OUT, (b.length / 1e6).toFixed(1), "MB,", order.length, "figures, print-size rule"); });
