"""Fill the manuscript Table 1 / Table 2 Word templates with values from tables/T1_T2_values.csv.

Keeps the template's formatting; replaces the text of the value cells, bolds p < 0.05
and writes a log of every cell whose value changed relative to the template.
Usage: python 09_fill_tables_docx.py <template_dir> <out_dir>
"""
import csv
import sys
from pathlib import Path

import docx

tmpl_dir, out_dir = Path(sys.argv[1]), Path(sys.argv[2])
out_dir.mkdir(parents=True, exist_ok=True)
vals = {}
with open("tables/T1_T2_values.csv", encoding="utf-8") as fh:
    for r in csv.DictReader(fh):
        vals[(r["Table"], r["Label"].strip())] = r


def set_text(cell, text, bold=None):
    para = cell.paragraphs[0]
    for extra in cell.paragraphs[1:]:
        extra._element.getparent().remove(extra._element)
    runs = para.runs
    if not runs:
        para.add_run(text)
        runs = para.runs
    runs[0].text = text
    for r in runs[1:]:
        r._element.getparent().remove(r._element)
    runs[0].bold = bold


changes = []
for tab, fname in [("T1", "Table1_clinical_characteristics.docx"), ("T2", "Table2_lifestyle_diet.docx")]:
    doc = docx.Document(tmpl_dir / fname)
    table = doc.tables[0]
    used = set()
    for row in table.rows[1:]:
        cells = row.cells
        if len({c._tc for c in cells}) < 5:          # merged section header or footnote row
            continue
        label = cells[0].text.strip()
        key = (tab, label)
        if key not in vals:
            raise SystemExit(f"No value for row: {key}")
        v = vals[key]
        used.add(key)
        new = [v["All"], v["IBS_C"], v["IBS_D"], v["p"]]
        p_num = float(v["p_num"]) if v["p_num"] not in ("", "NA") else None
        for j, txt in enumerate(new, start=1):
            old = cells[j].text.strip()
            bold = True if (j == 4 and p_num is not None and p_num < 0.05) else None
            if old != txt:
                changes.append((tab, label, ["All", "IBS-C", "IBS-D", "p"][j - 1], old, txt))
            if old or txt:
                set_text(cells[j], txt, bold)
    missing = [k for k in vals if k[0] == tab and k not in used]
    if missing:
        raise SystemExit(f"Values not placed: {missing}")
    doc.save(out_dir / fname)

with open(out_dir / "changes_vs_previous_version.csv", "w", encoding="utf-8", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["Table", "Row", "Column", "Previous", "Updated"])
    w.writerows(changes)
print(f"{len(changes)} cells changed")
for c in changes:
    print(" | ".join(c))
