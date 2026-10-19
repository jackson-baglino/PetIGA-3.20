"""Render a report draft: fill {F:id} placeholders from the fact ledger and write Markdown, a LaTeX/PDF file and a Word file (via macOS textutil)."""
from __future__ import annotations

import html
import re
import subprocess
from pathlib import Path

from agentic.writing.checks import load_facts, render


def md_inline(t: str) -> str:
    t = html.escape(t)
    t = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", t); t = re.sub(r"(?<!\*)\*(?!\*)(.+?)(?<!\*)\*(?!\*)", r"<i>\1</i>", t)
    t = re.sub(r"`([^`]+)`", r"<code>\1</code>", t)
    return re.sub(r"\[([^\]]+)\]\((https?://[^)]+)\)", r'<a href="\2">\1</a>', t)


def to_html(md: str) -> str:
    out, in_ul = ["<html><head><meta charset='utf-8'><style>body{font-family:Times New Roman;font-size:11pt;line-height:1.15} h1{font-size:15pt} h2{font-size:12pt;margin-bottom:2pt} p{margin:3pt 0}</style></head><body>"], False
    for line in md.splitlines():
        if re.match(r"^\s*[-*]\s+", line):
            if not in_ul: out.append("<ul>"); in_ul = True
            out.append("<li>" + md_inline(re.sub(r"^\s*[-*]\s+", "", line)) + "</li>"); continue
        if in_ul: out.append("</ul>"); in_ul = False
        if line.startswith("# "): out.append("<h1>" + md_inline(line[2:]) + "</h1>")
        elif line.startswith("## "): out.append("<h2>" + md_inline(line[3:]) + "</h2>")
        elif line.strip(): out.append("<p>" + md_inline(line) + "</p>")
    if in_ul: out.append("</ul>")
    return "\n".join(out) + "</body></html>"


def write_all(draft: str, facts_path, out_base: Path) -> dict:
    text = render(draft, load_facts(facts_path))
    out_base.parent.mkdir(parents=True, exist_ok=True)
    md = out_base.with_suffix(".md"); md.write_text(text)
    h = out_base.with_suffix(".html"); h.write_text(to_html(text))
    res = {"md": str(md)}
    p = subprocess.run(["textutil", "-convert", "docx", "-output", str(out_base.with_suffix(".docx")), str(h)], capture_output=True, text=True)
    res["docx"] = str(out_base.with_suffix(".docx")) if p.returncode == 0 else f"failed: {p.stderr[:80]}"
    # PDF through the browser-free route: LaTeX (pdflatex) from the same Markdown
    tex = ["\\documentclass[11pt]{article}\\usepackage[margin=1in]{geometry}\\usepackage[T1]{fontenc}\\usepackage{lmodern}\\usepackage{enumitem}\\usepackage[hidelinks]{hyperref}\\setlength{\\parindent}{0pt}\\setlength{\\parskip}{4pt}\\setlist{nosep}\\pagestyle{empty}\\begin{document}"]
    def esc(s): 
        s = s.replace("\\", "\\textbackslash{}")
        for a, b in (("&", "\\&"), ("%", "\\%"), ("#", "\\#"), ("_", "\\_"), ("$", "\\$"), ("{", "\\{"), ("}", "\\}"), ("~", "\\~{}"), ("^", "\\^{}")): s = s.replace(a, b)
        s = re.sub(r"\*\*(.+?)\*\*", r"\\textbf{\1}", s); s = re.sub(r"(?<!\*)\*(?!\*)(.+?)(?<!\*)\*(?!\*)", r"\\emph{\1}", s)
        return re.sub(r"`([^`]+)`", r"\\texttt{\1}", s)
    in_ul = False
    for line in text.splitlines():
        if re.match(r"^\s*[-*]\s+", line):
            if not in_ul: tex.append("\\begin{itemize}"); in_ul = True
            tex.append("\\item " + esc(re.sub(r"^\s*[-*]\s+", "", line))); continue
        if in_ul: tex.append("\\end{itemize}"); in_ul = False
        if line.startswith("# "): tex.append("{\\large\\bfseries " + esc(line[2:]) + "}\\par")
        elif line.startswith("## "): tex.append("\\medskip{\\bfseries " + esc(line[3:]) + "}\\par")
        elif line.strip(): tex.append(esc(line) + "\\par")
    if in_ul: tex.append("\\end{itemize}")
    tex.append("\\end{document}")
    t = out_base.with_suffix(".tex"); t.write_text("\n".join(tex))
    p = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-output-directory", str(out_base.parent), str(t)], capture_output=True, text=True)
    pdf = out_base.with_suffix(".pdf")
    if pdf.exists():
        info = subprocess.run(["mutool", "info", str(pdf)], capture_output=True, text=True).stdout
        m = re.search(r"Pages: (\d+)", info); res["pdf"] = str(pdf); res["pages"] = int(m.group(1)) if m else None
    else:
        res["pdf"] = "failed: " + p.stdout[-200:]
    return res
