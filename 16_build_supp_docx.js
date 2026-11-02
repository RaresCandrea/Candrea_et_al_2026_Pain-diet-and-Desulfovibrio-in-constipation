// Render a supplementary-table JSON spec (from 15_supp_tables.R) to Word, in the style of Tables 1-2:
// A4, Times New Roman 10 pt, grey header row (D9D9D9), horizontal rules only.
// Usage: node 16_build_supp_docx.js tables/supp/S1.json tables/supp
const fs = require("fs");
const path = require("path");
const {
  Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell, WidthType, AlignmentType,
  BorderStyle, ShadingType, PageOrientation, VerticalAlign,
} = require("docx");

const [specFile, outDir] = process.argv.slice(2);
const spec = JSON.parse(fs.readFileSync(specFile, "utf8"));
const FONT = "Times New Roman";
const SIZE = 20; // half-points = 10 pt
const NOTE = 18; // 9 pt

// Markup: *italic*, _{subscript}, ^{superscript}
function runs(text, opts = {}) {
  const out = [];
  let italic = false, buf = "", i = 0;
  const flush = () => { if (buf) { out.push(new TextRun({ text: buf, font: FONT, size: opts.size || SIZE, bold: opts.bold, italics: italic })); buf = ""; } };
  while (i < text.length) {
    const c = text[i];
    if (c === "*") { flush(); italic = !italic; i++; continue; }
    if ((c === "_" || c === "^") && text[i + 1] === "{") {
      const end = text.indexOf("}", i + 2);
      flush();
      out.push(new TextRun({ text: text.slice(i + 2, end), font: FONT, size: opts.size || SIZE, bold: opts.bold, italics: italic,
                             subScript: c === "_", superScript: c === "^" }));
      i = end + 1; continue;
    }
    buf += c; i++;
  }
  flush();
  return out;
}

const rule = { style: BorderStyle.SINGLE, size: 8, color: "000000" };
const none = { style: BorderStyle.NONE, size: 0, color: "FFFFFF" };
function cell(text, width, align, { header = false, top = false, bottom = false } = {}) {
  return new TableCell({
    width: { size: width, type: WidthType.DXA },
    verticalAlign: VerticalAlign.CENTER,
    shading: header ? { fill: "D9D9D9", type: ShadingType.CLEAR, color: "auto" } : undefined,
    margins: { top: 50, bottom: 50, left: 90, right: 90 },
    borders: { top: top ? rule : none, bottom: bottom ? rule : none, left: none, right: none },
    children: [new Paragraph({ alignment: align === "left" ? AlignmentType.LEFT : AlignmentType.CENTER,
                               children: runs(String(text), { bold: header }) })],
  });
}

const children = [new Paragraph({ spacing: { after: 120 }, children: [
  new TextRun({ text: spec.label + " ", font: FONT, size: SIZE, bold: true }), ...runs(spec.title)] })];

for (const b of spec.blocks) {
  const total = b.widths.reduce((a, x) => a + x, 0);
  if (b.caption) children.push(new Paragraph({ spacing: { before: 160, after: 60 }, keepNext: true,
                                               children: runs(b.caption, { bold: true }) }));
  const head = new TableRow({ tableHeader: true, children: b.headers.map((h, j) =>
    cell(h, b.widths[j], b.align[j], { header: true, top: true, bottom: true })) });
  const body = b.rows.map((r, i) => {
    const last = i === b.rows.length - 1;
    if (r.type === "section") {
      return new TableRow({ children: [new TableCell({
        columnSpan: b.widths.length, width: { size: total, type: WidthType.DXA },
        shading: { fill: "F2F2F2", type: ShadingType.CLEAR, color: "auto" },
        margins: { top: 50, bottom: 50, left: 90, right: 90 },
        borders: { top: none, bottom: last ? rule : none, left: none, right: none },
        children: [new Paragraph({ children: runs(Array.isArray(r.cells) ? r.cells[0] : r.cells, { bold: true }) })] })] });
    }
    return new TableRow({ children: r.cells.map((t, j) => cell(t, b.widths[j], b.align[j], { bottom: last })) });
  });
  children.push(new Table({ width: { size: total, type: WidthType.DXA }, columnWidths: b.widths,
                            alignment: AlignmentType.CENTER, rows: [head, ...body] }));
}
children.push(new Paragraph({ spacing: { before: 120 }, children: runs(spec.footnote, { size: NOTE }) }));

const landscape = !!spec.landscape;
const doc = new Document({ sections: [{ properties: { page: {
  size: { width: 11909, height: 16834, orientation: landscape ? PageOrientation.LANDSCAPE : PageOrientation.PORTRAIT },
  margin: { top: 1152, right: 792, bottom: 1152, left: 792 } } }, children }] });
Packer.toBuffer(doc).then((buf) => {
  const out = path.join(outDir, spec.file + ".docx");
  fs.writeFileSync(out, buf);
  console.log("written", out);
});
