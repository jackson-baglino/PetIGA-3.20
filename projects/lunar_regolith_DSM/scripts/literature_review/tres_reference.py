#!/usr/bin/env python3
"""T-res reference items from an arXiv LaTeX source (WORKFLOW.md §12).

For the pages listed in calibration/tres_pages.yaml with `source: latex`, finds the
numbered displayed equations printed on each PDF page and takes their LaTeX from the
authors' source, with the source's own macros expanded. Writes a JSON file
{page_pdf: [{label, latex}]} into an untrusted folder; prints counts only.

    tres_reference.py --tex <main.tex> --pdf <arXiv pdf> --pages 8 14 19 --out <file.json>

Numbering: equation and align environments in document order, one number per equation
and per align row, none where \\nonumber or \\notag appears; such rows belong to the
next numbered row, since a multi-line equation carries its number on its last line. A page's equations are the
"(n)" labels printed in the right quarter of the page (from mutool's structured text).
"""
import argparse
import json
import re
import subprocess
import xml.etree.ElementTree as ET


def strip_comments(tex):
    return re.sub(r"(?<!\\)%.*", "", tex)


def braced(s, i):
    """(content, index after the closing brace) for the group opening at s[i] == '{'."""
    depth = 0
    for j in range(i, len(s)):
        if s[j] == "{" and (j == 0 or s[j - 1] != "\\"):
            depth += 1
        elif s[j] == "}" and s[j - 1] != "\\":
            depth -= 1
            if depth == 0:
                return s[i + 1:j], j + 1
    raise ValueError("unbalanced braces")


def macros(tex):
    """{name: body} for argument-free \\newcommand, \\renewcommand and \\def."""
    out = {}
    for m in re.finditer(r"\\(?:re)?newcommand\*?\s*(?:\{\\([A-Za-z]+)\}|\\([A-Za-z]+))\s*(\[\d\])?\s*\{", tex):
        if m.group(3):
            continue                                        # macros with arguments are left alone
        out[m.group(1) or m.group(2)] = braced(tex, m.end() - 1)[0]
    for m in re.finditer(r"\\def\s*\\([A-Za-z]+)\s*\{", tex):
        out[m.group(1)] = braced(tex, m.end() - 1)[0]
    return out


def expand(s, defs):
    for _ in range(10):
        new = s
        for name in sorted(defs, key=len, reverse=True):
            new = re.sub(r"\\" + name + r"(?![A-Za-z])", lambda m, body=defs[name]: body, new)
        if new == s:
            return s
        s = new
    return s


def rows(body):
    """Split an align body at top-level \\\\ (not inside braces or nested environments)."""
    out, depth, env, start, i = [], 0, 0, 0, 0
    while i < len(body):
        if body.startswith("\\begin{", i):
            env += 1
        elif body.startswith("\\end{", i):
            env -= 1
        c = body[i]
        if c == "{" and body[i - 1:i] != "\\":
            depth += 1
        elif c == "}" and body[i - 1:i] != "\\":
            depth -= 1
        elif body.startswith("\\\\", i) and depth == 0 and env == 0:
            out.append(body[start:i])
            start = i + 2
            i += 2
            continue
        i += 1
    out.append(body[start:])
    return [r for r in out if r.strip()]


def numbered_equations(tex):
    """[(number, latex)] in document order."""
    doc = tex[tex.index("\\begin{document}"):]
    eqs, n = [], 0
    for m in re.finditer(r"\\begin\{(equation|align)\}(.*?)\\end\{\1\}", doc, flags=re.S):
        parts = [m.group(2)] if m.group(1) == "equation" else rows(m.group(2))
        pending = []                                        # unnumbered rows continue the next number
        for part in parts:
            part = re.sub(r"\\label\{[^}]*\}", "", part).replace("&", " ").strip()
            if re.search(r"\\(nonumber|notag)\b", part):
                pending.append(re.sub(r"\\(nonumber|notag)\b", "", part))
                continue
            n += 1
            eqs.append((n, " ".join(pending + [part])))
            pending = []
    return eqs


def printed_numbers(pdf, page):
    """Equation numbers "(n)" printed in the right quarter of a PDF page."""
    xml = subprocess.run(["mutool", "draw", "-q", "-F", "stext", "-o", "-", pdf, str(page)],
                         capture_output=True, text=True, check=True).stdout
    root = ET.fromstring(xml)
    pg = root.find("page")
    width = float(pg.get("width"))
    found = []
    for line in pg.iter("line"):
        x0 = float(line.get("bbox").split()[0])
        text = "".join(ch.get("c", "") for ch in line.iter("char")).strip()
        m = re.fullmatch(r"\((\d+)\)", text)
        if m and x0 > 0.75 * width:
            found.append(int(m.group(1)))
    return found


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--tex", required=True)
    ap.add_argument("--pdf", required=True)
    ap.add_argument("--pages", type=int, nargs="+", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    tex = strip_comments(open(a.tex).read())
    defs = macros(tex)
    eqs = dict(numbered_equations(tex))
    out, summary = {}, {"macros": len(defs), "numbered_equations": len(eqs), "pages": {}}
    for p in a.pages:
        nums = printed_numbers(a.pdf, p)
        items = [{"label": f"({k})", "latex": expand(eqs[k], defs)} for k in nums if k in eqs]
        out[str(p)] = items
        left = sum(len(re.findall(r"\\(" + "|".join(defs) + r")(?![A-Za-z])", i["latex"])) for i in items) if defs else 0
        summary["pages"][p] = {"printed_numbers": nums, "items": len(items), "unexpanded_macros": left}
    with open(a.out, "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
