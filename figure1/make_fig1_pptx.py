"""Figure 1 (flow diagram) as an editable PowerPoint, reproducing the author's PDF (text, numbers, fonts, colours).
Boxes are text boxes; arrows are real connectors (editable)."""
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE, MSO_CONNECTOR
from pptx.oxml.ns import qn
from lxml import etree

FONT = "Times New Roman"
BLACK = RGBColor(0, 0, 0)
prs = Presentation()
prs.slide_width, prs.slide_height = Inches(7.2), Inches(7.35)
s = prs.slides.add_slide(prs.slide_layouts[6])


def box(x, y, w, h, lines, fill=None, align=PP_ALIGN.CENTER, lw=1.0):
    sh = s.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(x), Inches(y), Inches(w), Inches(h))
    sh.shadow.inherit = False
    if fill:
        sh.fill.solid(); sh.fill.fore_color.rgb = RGBColor.from_string(fill)
    else:
        sh.fill.solid(); sh.fill.fore_color.rgb = RGBColor(255, 255, 255)
    sh.line.color.rgb = BLACK; sh.line.width = Pt(lw)
    tf = sh.text_frame; tf.word_wrap = True; tf.vertical_anchor = MSO_ANCHOR.MIDDLE
    for m in ("margin_left", "margin_right"):
        setattr(tf, m, Inches(0.08))
    tf.margin_top = tf.margin_bottom = Inches(0.03)
    for i, (txt, bold) in enumerate(lines):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align
        r = p.add_run(); r.text = txt
        r.font.name = FONT; r.font.size = Pt(11); r.font.bold = bold; r.font.color.rgb = BLACK
    return sh


def arrow(x1, y1, x2, y2, head=True):
    c = s.shapes.add_connector(MSO_CONNECTOR.STRAIGHT, Inches(x1), Inches(y1), Inches(x2), Inches(y2))
    c.line.color.rgb = BLACK; c.line.width = Pt(1)
    if head:
        ln = c.line._get_or_add_ln()
        tail = etree.SubElement(ln, qn("a:tailEnd")); tail.set("type", "triangle"); tail.set("w", "med"); tail.set("len", "med")
    return c


X, W = 0.3, 2.9           # main column
XC = X + W / 2            # centre line
EX, EW = 3.75, 3.1        # exclusion boxes
L = PP_ALIGN.LEFT

b1 = box(X, 0.30, W, 0.95, [("Patients with gut microbiota analysis at the clinic", False), ("May 2021 – October 2026", False), ("(n=125)", False)])
e1 = box(EX, 1.38, EW, 0.70, [("Excluded (n=16)", True), ("Microbiota analysed in another laboratory or with a different method", False)], fill="F2F2F2", align=L)
b2 = box(X, 2.25, W, 0.55, [("Eligible microbiota records", False), ("(n=109)", False)])
e2 = box(EX, 2.98, EW, 0.55, [("Excluded (n=21)", True), ("No physician-confirmed diagnosis of IBS", False)], fill="F2F2F2", align=L)
b3 = box(X, 3.70, W, 0.55, [("Patients with IBS assessed for eligibility", False), ("(n=88)", False)])
e3 = box(EX, 4.43, EW, 0.95, [("Excluded (n=24)", True), ("IBS-M or IBS-U subtype (n=20)", False), ("Celiac disease (n=2)", False),
                              ("Incomplete clinical or microbiota data (n=2)", False)], fill="F2F2F2", align=L)
b4 = box(X, 5.55, W, 0.55, [("Included in the analysis", False), ("(n=64)", False)])
cw = 1.35
bc = box(X, 6.55, cw, 0.55, [("IBS-C", False), ("(n=35)", False)], fill="DCE6F1", lw=1.25)
bd = box(X + W - cw, 6.55, cw, 0.55, [("IBS-D", False), ("(n=29)", False)], fill="F8DCD9", lw=1.25)

# main flow
arrow(XC, 1.25, XC, 2.25); arrow(XC, 2.80, XC, 3.70); arrow(XC, 4.25, XC, 5.55)
# exclusions
arrow(XC, 1.73, EX, 1.73); arrow(XC, 3.255, EX, 3.255); arrow(XC, 4.905, EX, 4.905)
# split into subtypes
arrow(XC, 6.10, XC, 6.32, head=False)
arrow(X + cw / 2, 6.32, X + W - cw / 2, 6.32, head=False)
arrow(X + cw / 2, 6.32, X + cw / 2, 6.55); arrow(X + W - cw / 2, 6.32, X + W - cw / 2, 6.55)

import os; os.makedirs("figures", exist_ok=True)
prs.save("figures/Figure1_flow_diagram.pptx")
print("saved")
