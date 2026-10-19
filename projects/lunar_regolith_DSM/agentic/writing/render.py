"""Render a report draft: fill {F:id} placeholders from the fact ledger, then write Markdown, a Word file (built directly, 11 pt Times) and a PDF (LaTeX, 11 pt).
Markdown subset: # title, ## heading, paragraphs, - bullets, **bold**, *italic*, `code`, pipe tables, a line `<<<PAGEBREAK>>>`.
"""
from __future__ import annotations

import re
import subprocess
import zipfile
from pathlib import Path
from xml.sax.saxutils import escape

from agentic.writing.checks import load_facts, render

FONT = "Times New Roman"; BODY_PT = 11
PAGEBREAK = "<<<PAGEBREAK>>>"


def blocks(text: str):
    """Yield ('title'|'h2'|'p'|'li'|'table'|'pagebreak', payload)."""
    lines = text.splitlines(); i = 0
    while i < len(lines):
        ln = lines[i]
        if ln.strip() == PAGEBREAK: yield "pagebreak", None
        elif ln.startswith("# "): yield "title", ln[2:]
        elif ln.startswith("## "): yield "h2", ln[3:]
        elif re.match(r"^\s*[-*]\s+", ln): yield "li", re.sub(r"^\s*[-*]\s+", "", ln)
        elif ln.startswith("|"):
            rows = []
            while i < len(lines) and lines[i].startswith("|"):
                cells = [c.strip() for c in lines[i].strip().strip("|").split("|")]
                if not all(re.fullmatch(r":?-{2,}:?", c) for c in cells): rows.append(cells)
                i += 1
            yield "table", rows; continue
        elif ln.strip(): yield "p", ln
        i += 1


# ---------------------------------------------------------------- Word
def runs(t: str, size: int, bold=False):
    out = []
    for part in re.split(r"(\*\*.+?\*\*|`[^`]+`|(?<!\*)\*(?!\*).+?(?<!\*)\*(?!\*))", t):
        if not part: continue
        b, it, mono = bold, False, False
        if part.startswith("**"): part, b = part[2:-2], True
        elif part.startswith("`"): part, mono = part[1:-1], True
        elif part.startswith("*") and part.endswith("*") and len(part) > 2: part, it = part[1:-1], True
        rpr = f'<w:rFonts w:ascii="{"Courier New" if mono else FONT}" w:hAnsi="{"Courier New" if mono else FONT}"/>' + ("<w:b/>" if b else "") + ("<w:i/>" if it else "") + f'<w:sz w:val="{size * 2}"/>'
        out.append(f'<w:r><w:rPr>{rpr}</w:rPr><w:t xml:space="preserve">{escape(part)}</w:t></w:r>')
    return "".join(out)


def para(t, size=BODY_PT, bold=False, after=80, before=0, ind=None, keep=False):
    ppr = f'<w:spacing w:before="{before}" w:after="{after}" w:line="252" w:lineRule="auto"/>' + (ind or "") + ("<w:keepNext/>" if keep else "")
    return f"<w:p><w:pPr>{ppr}</w:pPr>{runs(t, size, bold)}</w:p>"


def docx_body(text: str) -> str:
    out = []
    for kind, pl in blocks(text):
        if kind == "title": out.append(para(pl, 14, True, after=100))
        elif kind == "h2": out.append(para(pl, 12, True, after=40, before=110, keep=True))
        elif kind == "p":
            hang = '<w:ind w:left="360" w:hanging="360"/>' if re.match(r"^[A-Z][^.]+\(\d{4}\)\.", pl) else None
            out.append(para(pl, BODY_PT, ind=hang))
        elif kind == "li": out.append(para("\u2022 " + pl, BODY_PT, after=40, ind='<w:ind w:left="360" w:hanging="260"/>'))
        elif kind == "pagebreak": out.append('<w:p><w:r><w:br w:type="page"/></w:r></w:p>')
        elif kind == "table":
            n = len(pl[0]); w = int(9360 / n)
            rows = []
            for r_i, r in enumerate(pl):
                cells = "".join(f'<w:tc><w:tcPr><w:tcW w:w="{w}" w:type="dxa"/></w:tcPr>{para(c, 10, r_i == 0, after=20)}</w:tc>' for c in r)
                rows.append(f"<w:tr>{cells}</w:tr>")
            b = "".join(f'<w:{s} w:val="single" w:sz="4" w:color="808080"/>' for s in ("top", "left", "bottom", "right", "insideH", "insideV"))
            out.append(f'<w:tbl><w:tblPr><w:tblW w:w="9360" w:type="dxa"/><w:tblBorders>{b}</w:tblBorders></w:tblPr>{"".join(rows)}</w:tbl>' + para("", 4, after=0))
    return "".join(out)


def write_docx(text: str, path: Path) -> None:
    ns = 'xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"'
    doc = f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?><w:document {ns}><w:body>{docx_body(text)}<w:sectPr><w:pgSz w:w="12240" w:h="15840"/><w:pgMar w:top="1440" w:right="1440" w:bottom="1440" w:left="1440" w:header="720" w:footer="720" w:gutter="0"/></w:sectPr></w:body></w:document>'
    styles = f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?><w:styles {ns}><w:docDefaults><w:rPrDefault><w:rPr><w:rFonts w:ascii="{FONT}" w:hAnsi="{FONT}"/><w:sz w:val="{BODY_PT * 2}"/></w:rPr></w:rPrDefault></w:docDefaults><w:style w:type="paragraph" w:default="1" w:styleId="Normal"><w:name w:val="Normal"/></w:style></w:styles>'
    ct = '<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/><Override PartName="/word/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.styles+xml"/></Types>'
    rels = '<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/></Relationships>'
    drels = '<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/></Relationships>'
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml", ct); z.writestr("_rels/.rels", rels); z.writestr("word/document.xml", doc); z.writestr("word/styles.xml", styles); z.writestr("word/_rels/document.xml.rels", drels)


# ---------------------------------------------------------------- PDF
def tex_esc(s: str) -> str:
    s = s.replace("\\", "\\textbackslash{}")
    for a, b in (("&", "\\&"), ("%", "\\%"), ("#", "\\#"), ("_", "\\_"), ("$", "\\$"), ("{", "\\{"), ("}", "\\}"), ("~", "\\~{}"), ("^", "\\^{}")):
        s = s.replace(a, b)
    s = re.sub(r"\*\*(.+?)\*\*", r"\\textbf{\1}", s); s = re.sub(r"(?<!\*)\*(?!\*)(.+?)(?<!\*)\*(?!\*)", r"\\emph{\1}", s)
    s = re.sub(r"`([^`]+)`", r"\\texttt{\1}", s)
    return re.sub(r"(https?://[^\s}]+)", lambda m: "\\url{" + m.group(1).replace("\\_", "_").replace("\\&", "&").replace("\\%", "%") + "}", s)


def write_pdf(text: str, out_base: Path) -> dict:
    tex = ["\\documentclass[11pt]{article}\\usepackage[margin=1in]{geometry}\\usepackage[utf8]{inputenc}\\usepackage[T1]{fontenc}\\usepackage{mathptmx}\\usepackage{tabularx}\\usepackage[hidelinks]{hyperref}\\usepackage{url}\\urlstyle{same}\\setlength{\\parindent}{0pt}\\setlength{\\parskip}{3pt}\\pagestyle{empty}\\begin{document}"]
    for kind, pl in blocks(text):
        if kind == "title": tex.append("{\\large\\bfseries " + tex_esc(pl) + "}\\par")
        elif kind == "h2": tex.append("\\medskip{\\bfseries " + tex_esc(pl) + "}\\par\\nopagebreak")
        elif kind == "p":
            hang = re.match(r"^[A-Z][^.]+\(\d{4}\)\.", pl) is not None
            tex.append(("\\hangindent=1.5em\\hangafter=1 " if hang else "") + tex_esc(pl) + "\\par")
        elif kind == "li": tex.append("\\hangindent=1.5em\\hangafter=1 \\textbullet\\ " + tex_esc(pl) + "\\par")
        elif kind == "pagebreak": tex.append("\\newpage")
        elif kind == "table":
            n = len(pl[0]); spec = "|" + "|".join(["X"] * n) + "|"
            rows = " \\\\\\hline\n".join(" & ".join(("\\textbf{" + tex_esc(c) + "}") if i == 0 else tex_esc(c) for c in r) for i, r in enumerate(pl))
            tex.append("{\\small\\renewcommand{\\arraystretch}{1.1}\\begin{tabularx}{\\textwidth}{" + spec + "}\\hline\n" + rows + " \\\\\\hline\\end{tabularx}}\\par")
    tex.append("\\end{document}")
    t = out_base.with_suffix(".tex"); t.write_text("\n".join(tex))
    p = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "-output-directory", str(out_base.parent), str(t)], capture_output=True, text=True, stdin=subprocess.DEVNULL, timeout=120)
    pdf = out_base.with_suffix(".pdf"); res = {}
    if pdf.exists() and p.returncode == 0:
        info = subprocess.run(["mutool", "info", str(pdf)], capture_output=True, text=True).stdout
        m = re.search(r"Pages: (\d+)", info); res = {"pdf": str(pdf), "pages": int(m.group(1)) if m else None}
    else:
        res = {"pdf": "failed: " + p.stdout[-300:]}
    for ext in (".aux", ".log", ".out"):
        out_base.with_suffix(ext).unlink(missing_ok=True)
    return res


def write_all(draft: str, facts_path, out_base: Path) -> dict:
    text = render(draft, load_facts(facts_path)); out_base = Path(out_base)
    out_base.parent.mkdir(parents=True, exist_ok=True)
    md = out_base.with_suffix(".md"); md.write_text(text)
    write_docx(text, out_base.with_suffix(".docx"))
    return {"md": str(md), "docx": str(out_base.with_suffix(".docx")), **write_pdf(text, out_base)}
