"""Word helpers for the run plan, copied from SimulationPlan_v2_build/build_doc.py
(lines 1-310, 2026-11-11): page setup, styles, LaTeX -> Word equations (omath,
mathify), paragraphs, tables, boxes and figures. Import * and use doc."""


from docx import Document
from docx.enum.section import WD_ORIENT  # noqa: F401
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_BREAK
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor


import os  # noqa: E402
import re  # noqa: E402

FOOTER_DATE = os.environ.get("RUNPLAN_DATE", "2026-11-11")

from mathify import mathify  # noqa: E402
from omath import omml  # noqa: E402

WNS = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
MNS = "http://schemas.openxmlformats.org/officeDocument/2006/math"


def add_math(par, latex, size=None, color=None, bold=False):
    """Append a Word equation (OMML) built from LaTeX, styled like the text around it."""
    om = omml(latex)
    for mr in om.iter(f"{{{MNS}}}r"):
        rpr = OxmlElement("w:rPr")
        f = OxmlElement("w:rFonts")
        f.set(qn("w:ascii"), "Cambria Math")
        f.set(qn("w:hAnsi"), "Cambria Math")
        rpr.append(f)
        if bold:
            rpr.append(OxmlElement("w:b"))
        if color is not None:
            c = OxmlElement("w:color")
            c.set(qn("w:val"), str(color))
            rpr.append(c)
        sz = OxmlElement("w:sz")
        sz.set(qn("w:val"), str(int(round((size or 10) * 2))))
        rpr.append(sz)
        mrpr = mr.find(f"{{{MNS}}}rPr")
        idx = list(mr).index(mrpr) + 1 if mrpr is not None else 0
        mr.insert(idx, rpr)
    par._p.append(om)

NAVY = RGBColor(0x1F, 0x3A, 0x5F)
GRAYT = RGBColor(0x55, 0x55, 0x55)
FONT = "Arial"

doc = Document()
sec = doc.sections[0]
sec.page_width, sec.page_height = Inches(8.5), Inches(11)
for side in ("left_margin", "right_margin"):
    setattr(sec, side, Inches(0.95))
sec.top_margin, sec.bottom_margin = Inches(0.85), Inches(0.85)

st = doc.styles["Normal"]
st.font.name = FONT
st.font.size = Pt(10)
st.element.rPr.rFonts.set(qn("w:eastAsia"), FONT)
st.paragraph_format.space_after = Pt(5)
st.paragraph_format.line_spacing = 1.12
for name, size, before, after, color in [("Title", 20, 0, 4, NAVY), ("Heading 1", 14.5, 16, 6, NAVY),
                                         ("Heading 2", 12, 12, 4, NAVY), ("Heading 3", 10.5, 9, 3, NAVY)]:
    s = doc.styles[name]
    s.font.name, s.font.size, s.font.bold, s.font.color.rgb = FONT, Pt(size), True, color
    s.font.italic = False
    s.element.rPr.rFonts.set(qn("w:eastAsia"), FONT)
    s.element.rPr.rFonts.set(qn("w:ascii"), FONT)
    s.element.rPr.rFonts.set(qn("w:hAnsi"), FONT)
    s.paragraph_format.space_before, s.paragraph_format.space_after = Pt(before), Pt(after)
    s.paragraph_format.keep_with_next = True
for name in ("List Bullet", "List Number"):
    doc.styles[name].font.name = FONT
    doc.styles[name].font.size = Pt(10)
    doc.styles[name].paragraph_format.space_after = Pt(2)

# footer page numbers
fp = sec.footer.paragraphs[0]
fp.alignment = WD_ALIGN_PARAGRAPH.CENTER
r = fp.add_run(f"{os.environ.get('DOC_FOOTER', 'Run plan, paper 2 / AGU 2026')} · {FOOTER_DATE} · page ")
r.font.size, r.font.color.rgb = Pt(8), GRAYT
for t, txt in (("begin", None), (None, "PAGE"), ("end", None)):
    run = fp.add_run()
    run.font.size, run.font.color.rgb = Pt(8), GRAYT
    if t:
        el = OxmlElement("w:fldChar")
        el.set(qn("w:fldCharType"), t)
    else:
        el = OxmlElement("w:instrText")
        el.set(qn("xml:space"), "preserve")
        el.text = txt
    run._r.append(el)


# ------------------------------------------------------------------ helpers
def shade(cell, hex_fill):
    tcPr = cell._tc.get_or_add_tcPr()
    shd = OxmlElement("w:shd")
    shd.set(qn("w:val"), "clear")
    shd.set(qn("w:color"), "auto")
    shd.set(qn("w:fill"), hex_fill)
    tcPr.append(shd)


def cell_borders(cell, color="BFBFBF", left=None):
    tcPr = cell._tc.get_or_add_tcPr()
    b = OxmlElement("w:tcBorders")
    for edge in ("top", "left", "bottom", "right"):
        e = OxmlElement(f"w:{edge}")
        if edge == "left" and left:
            e.set(qn("w:val"), "single")
            e.set(qn("w:sz"), "24")
            e.set(qn("w:color"), left)
        else:
            e.set(qn("w:val"), "nil")
        b.append(e)
    tcPr.append(b)


def rich(par, text, size=None, color=None, italic=False, bold=False):
    """**bold** segments and $LaTeX$ math (Word equations) inside text; symbols are found by mathify."""
    text = mathify(text)
    on = bold
    for tok in re.split(r"(\*\*|\$[^$]+\$)", text):
        if not tok:
            continue
        if tok == "**":
            on = not on
        elif tok.startswith("$") and tok.endswith("$") and len(tok) > 1:
            add_math(par, tok[1:-1], size=size, color=color, bold=on)
        else:
            r = par.add_run(tok)
            r.bold = on
            r.italic = italic
            if size:
                r.font.size = Pt(size)
            if color:
                r.font.color.rgb = color
    return par


def P(text, size=None, after=None, italic=False, color=None, align=None):
    p = doc.add_paragraph()
    rich(p, text, size=size, italic=italic, color=color)
    if after is not None:
        p.paragraph_format.space_after = Pt(after)
    if align:
        p.alignment = align
    return p


def B(text, level=0, style="List Bullet"):
    p = doc.add_paragraph(style=style)
    rich(p, text)
    if level:
        p.paragraph_format.left_indent = Inches(0.25 + 0.25 * level)
    return p


def H(text, level):
    return doc.add_heading(text, level=level)


def set_widths(table, widths):
    """Fixed layout with explicit grid, table and cell widths (LibreOffice honors only this)."""
    table.autofit = False
    tbl = table._tbl
    tblPr = tbl.tblPr
    lay = OxmlElement("w:tblLayout")
    lay.set(qn("w:type"), "fixed")
    tblPr.append(lay)
    tw = OxmlElement("w:tblW")
    tw.set(qn("w:w"), str(int(sum(widths) * 1440)))
    tw.set(qn("w:type"), "dxa")
    for old in tblPr.findall(qn("w:tblW")):
        tblPr.remove(old)
    tblPr.append(tw)
    grid = tbl.tblGrid
    for i, gc in enumerate(grid.findall(qn("w:gridCol"))):
        gc.set(qn("w:w"), str(int(widths[i] * 1440)))
    for row in table.rows:
        for i, w in enumerate(widths):
            row.cells[i].width = Inches(w)


def T(header, rows, widths, size=8.5, head_fill="1F3A5F", zebra=True, bold_first=False, head_color=(255, 255, 255)):
    t = doc.add_table(rows=1, cols=len(header))
    t.style = "Table Grid"
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    for i, h in enumerate(header):
        c = t.rows[0].cells[i]
        c.text = ""
        rich(c.paragraphs[0], h, size=size, color=RGBColor(*head_color), bold=True)
        shade(c, head_fill)
    # repeat header row
    trPr = t.rows[0]._tr.get_or_add_trPr()
    th = OxmlElement("w:tblHeader")
    th.set(qn("w:val"), "true")
    trPr.append(th)
    for k, row in enumerate(rows):
        cells = t.add_row().cells
        for i, v in enumerate(row):
            cells[i].text = ""
            lines = str(v).split("\n")
            for j, line in enumerate(lines):
                par = cells[i].paragraphs[0] if j == 0 else cells[i].add_paragraph()
                par.paragraph_format.space_after = Pt(0)
                par.paragraph_format.line_spacing = 1.0
                rich(par, f"**{line}**" if (bold_first and i == 0 and line) else line, size=size)
            if zebra and k % 2 == 1:
                shade(cells[i], "F2F5F9")
    for row in t.rows:
        for c in row.cells:
            for par in c.paragraphs:
                par.paragraph_format.space_after = Pt(0)
    no_split(t)
    set_widths(t, widths)
    doc.add_paragraph().paragraph_format.space_after = Pt(2)
    return t


def no_split(table):
    for row in table.rows:
        trPr = row._tr.get_or_add_trPr()
        cs = OxmlElement("w:cantSplit")
        cs.set(qn("w:val"), "true")
        trPr.append(cs)


def box(title, lines, fill, edge, title_color=NAVY, width=6.6):
    t = doc.add_table(rows=1, cols=1)
    no_split(t)
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    c = t.rows[0].cells[0]
    set_widths(t, [width])
    shade(c, fill)
    cell_borders(c, left=edge)
    p = c.paragraphs[0]
    r = p.add_run(title)
    r.bold, r.font.size, r.font.color.rgb = True, Pt(9.5), title_color
    p.paragraph_format.space_after = Pt(2)
    for line in lines:
        q = c.add_paragraph()
        q.paragraph_format.space_after = Pt(2)
        rich(q, line, size=9.5)
    doc.add_paragraph().paragraph_format.space_after = Pt(1)


def question(qid, text):
    box(f"Science question {qid}", [text], "EAF2FB", "2E6DB4")


def hypothesis(qids):
    box(f"Hypothesis for {qids} (to be written by Jackson)", [" ", " "], "FFF8E1", "E0A800",
        title_color=RGBColor(0x8A, 0x63, 0x00))


def figure(path, caption, width=3.3):
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.keep_with_next = True
    p.add_run().add_picture(path, width=Inches(width))
    c = doc.add_paragraph()
    c.alignment = WD_ALIGN_PARAGRAPH.CENTER
    rich(c, caption, size=8.5, italic=True, color=GRAYT)


def figure_beside(path, caption, bullets_title, bullets, img_w=None):
    """Sketch on the left at its natural size (fonts print as drawn), setup summary on the right."""
    if img_w is None:
        img_w = float(open(path[:-4] + ".width").read())
    t = doc.add_table(rows=1, cols=2)
    no_split(t)
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    a, b = t.rows[0].cells
    set_widths(t, [img_w + 0.3, 6.6 - img_w - 0.3])
    for c in (a, b):
        cell_borders(c)
    pa = a.paragraphs[0]
    pa.alignment = WD_ALIGN_PARAGRAPH.CENTER
    pa.add_run().add_picture(path, width=Inches(img_w))
    cap = a.add_paragraph()
    cap.alignment = WD_ALIGN_PARAGRAPH.CENTER
    rich(cap, caption, size=8, italic=True, color=GRAYT)
    pb = b.paragraphs[0]
    rr = pb.add_run(bullets_title)
    rr.bold, rr.font.size, rr.font.color.rgb = True, Pt(9.5), NAVY
    for line in bullets:
        q = b.add_paragraph()
        q.paragraph_format.space_after = Pt(3)
        q.paragraph_format.left_indent = Inches(0.12)
        q.paragraph_format.first_line_indent = Inches(-0.12)
        rich(q, "• " + line, size=9)
    doc.add_paragraph().paragraph_format.space_after = Pt(1)


def page_break():
    doc.add_paragraph().add_run().add_break(WD_BREAK.PAGE)

