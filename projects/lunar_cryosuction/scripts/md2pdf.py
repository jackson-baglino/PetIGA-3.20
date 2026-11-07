"""Render a Markdown file to PDF with headless Google Chrome.

No pandoc or Python markdown package is installed (2026-11-07), so this
converts the subset our batch reports use to HTML itself:

    # .. #### headings, paragraphs, blank-line separated
    - / * bullet lists and 1. numbered lists (one level)
    | pipe | tables |  with a |---| separator row
    ![caption](relative/path.png)  images (paths relative to the .md)
    **bold**, *italic*, `code`, [text](url), > quotes, --- rules

Anything else passes through as text. Unicode (θ, χ, µ, ≈) is fine.
The HTML is written next to the .md (so image paths resolve), printed, then
removed unless --keep-html.

    python3 scripts/md2pdf.py REPORT.md [--out REPORT.pdf] [--landscape] [--keep-html]
"""
from __future__ import annotations

import argparse
import html
import os
import re
import subprocess
import sys

CHROME = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"

CSS = """
@page { size: Letter %(orient)s; margin: 0.6in; }
body { font-family: -apple-system, "Helvetica Neue", Arial, sans-serif; font-size: 10pt;
       line-height: 1.4; color: #1c2430; }
h1 { font-size: 17pt; margin: 0 0 6pt; } h2 { font-size: 13pt; margin: 16pt 0 4pt;
  border-bottom: 1px solid #c8ced6; } h3 { font-size: 11pt; margin: 12pt 0 3pt; }
h4 { font-size: 10pt; margin: 10pt 0 2pt; }
table { border-collapse: collapse; margin: 6pt 0 10pt; font-size: 8pt; page-break-inside: auto; }
tr { page-break-inside: avoid; }
th, td { border: 1px solid #c8ced6; padding: 2pt 4pt; vertical-align: top; text-align: left; }
th { background: #eef1f5; }
code { font-family: Menlo, monospace; font-size: 8.5pt; background: #f3f5f8; padding: 0 2pt; }
img { display: block; margin: 6pt auto 2pt; }
figure { margin: 6pt 0 10pt; page-break-inside: avoid; }
figcaption { font-size: 8.5pt; color: #5b6673; text-align: center; }
blockquote { border-left: 3px solid #c8ced6; margin: 6pt 0; padding: 0 8pt; color: #3a4552; }
ul, ol { margin: 3pt 0 6pt 18pt; padding: 0; } li { margin: 1pt 0; }
"""


def inline(s: str) -> str:
    s = html.escape(s, quote=False)
    s = re.sub(r"`([^`]+)`", r"<code>\1</code>", s)
    s = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", s)
    s = re.sub(r"(?<![\w*])\*([^*\s][^*]*)\*(?![\w*])", r"<em>\1</em>", s)
    s = re.sub(r"\[([^\]]+)\]\(([^)\s]+)\)", r'<a href="\2">\1</a>', s)
    return s


def td(c: str) -> str:
    """A table cell; short ones (IDs, values with units) do not wrap."""
    style = ' style="white-space:nowrap"' if len(c) <= 14 else ""
    return f"<td{style}>{inline(c)}</td>"


def cells(line: str) -> list[str]:
    return [c.strip() for c in line.strip().strip("|").split("|")]


PAGE_IN = {"portrait": (8.5 - 1.2, 11.0 - 1.2 - 0.6), "landscape": (11.0 - 1.2, 8.5 - 1.2 - 0.6)}


def image_size(path: str, orient: str) -> tuple[float, float, float]:
    """(width_in, height_in, scale) to print a PNG at its own DPI, shrunk only
    if it does not fit the page (scale < 1 then, and fonts shrink with it)."""
    from PIL import Image
    with Image.open(path) as im:
        dpi = im.info.get("dpi", (96, 96))[0] or 96
        w, h = im.size[0] / dpi, im.size[1] / dpi
    pw, ph = PAGE_IN[orient]
    k = min(1.0, pw / w, ph / h)
    return w * k, h * k, k


def convert(md: str, base: str = ".", orient: str = "portrait") -> str:
    out: list[str] = []
    lines = md.splitlines()
    i, para = 0, []
    stack: list[tuple[int, str]] = []    # open lists: (indent, "ul"/"ol"), each with an open <li>

    def close_lists(to_indent=-1):
        while stack and stack[-1][0] > to_indent:
            out.append(f"</li></{stack.pop()[1]}>")

    def flush():
        nonlocal para
        if para:
            out.append("<p>" + inline(" ".join(para)) + "</p>")
            para = []
        close_lists()

    while i < len(lines):
        ln = lines[i]
        s = ln.strip()
        if not s:
            # A blank line ends a list only if the next text is not another item.
            j = i + 1
            while j < len(lines) and not lines[j].strip():
                j += 1
            nxt = lines[j] if j < len(lines) else ""
            if stack and re.match(r"^\s*([-*]|\d+\.)\s", nxt):
                if para:
                    out.append("<p>" + inline(" ".join(para)) + "</p>")
                    para = []
            else:
                flush()
        elif s.startswith("<!--"):
            while "-->" not in lines[i]:
                i += 1
        elif m := re.match(r"^(#{1,4})\s+(.*)", s):
            flush()
            n = len(m.group(1))
            out.append(f"<h{n}>{inline(m.group(2))}</h{n}>")
        elif re.match(r"^-{3,}$", s):
            flush()
            out.append("<hr>")
        elif m := re.match(r"^!\[([^\]]*)\]\(([^)]+)\)", s):
            flush()
            src = m.group(2)
            w, h, k = image_size(os.path.join(base, src), orient)
            if k < 0.999:
                print(f"md2pdf: {src} shrunk to {k:.2f} of its size to fit the page "
                      f"(a 7 pt font prints at {7 * k:.1f} pt)", file=sys.stderr)
            out.append(f'<figure><img src="{html.escape(src)}" style="width:{w:.3f}in;height:{h:.3f}in">'
                       f"<figcaption>{inline(m.group(1))}</figcaption></figure>")
        elif s.startswith("|") and i + 1 < len(lines) and re.match(r"^\|?\s*:?-{2,}", lines[i + 1].strip()):
            flush()
            head = cells(s)
            out.append("<table><thead><tr>" + "".join(f"<th>{inline(c)}</th>" for c in head)
                       + "</tr></thead><tbody>")
            i += 2
            while i < len(lines) and lines[i].strip().startswith("|"):
                out.append("<tr>" + "".join(td(c) for c in cells(lines[i])) + "</tr>")
                i += 1
            out.append("</tbody></table>")
            continue
        elif m := re.match(r"^([-*]|\d+\.)\s+(.*)", s):
            if para:
                out.append("<p>" + inline(" ".join(para)) + "</p>")
                para = []
            indent = len(ln) - len(ln.lstrip())
            kind = "ol" if m.group(1)[0].isdigit() else "ul"
            close_lists(indent)                    # leave deeper lists
            if stack and stack[-1][0] == indent and stack[-1][1] != kind:
                close_lists(indent - 1)            # same level, other list type
            if stack and stack[-1][0] == indent:
                out.append("</li>")                # next item of the same list
            else:
                out.append(f"<{kind}>")            # new (possibly nested) list
                stack.append((indent, kind))
            item = [m.group(2)]
            while i + 1 < len(lines) and lines[i + 1].strip() \
                    and len(lines[i + 1]) - len(lines[i + 1].lstrip()) > indent \
                    and not re.match(r"^([-*]|\d+\.)\s", lines[i + 1].strip()):
                i += 1
                item.append(lines[i].strip())
            out.append("<li>" + inline(" ".join(item)))
        elif s.startswith(">"):
            flush()
            q = []
            while i < len(lines) and lines[i].strip().startswith(">"):
                q.append(lines[i].strip()[1:].strip())
                i += 1
            out.append("<blockquote>" + inline(" ".join(q)) + "</blockquote>")
            continue
        else:
            close_lists()
            para.append(s)
        i += 1
    flush()
    return "\n".join(out)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("md")
    p.add_argument("--out")
    p.add_argument("--landscape", action="store_true")
    p.add_argument("--keep-html", action="store_true")
    a = p.parse_args(argv)
    md_path = os.path.abspath(a.md)
    out = os.path.abspath(a.out or os.path.splitext(md_path)[0] + ".pdf")
    orient = "landscape" if a.landscape else "portrait"
    body = convert(open(md_path, encoding="utf-8").read(), os.path.dirname(md_path), orient)
    css = CSS % {"orient": orient}
    page = f'<!doctype html><html><head><meta charset="utf-8"><style>{css}</style></head><body>{body}</body></html>'
    tmp = os.path.splitext(md_path)[0] + ".render.html"
    with open(tmp, "w", encoding="utf-8") as fh:
        fh.write(page)
    if not os.path.exists(CHROME):
        print(f"md2pdf: Chrome not found at {CHROME}; HTML left at {tmp}", file=sys.stderr)
        return 1
    cmd = [CHROME, "--headless=new", "--disable-gpu", "--no-pdf-header-footer",
           f"--print-to-pdf={out}", "file://" + tmp]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
    if not a.keep_html:
        os.remove(tmp)
    if r.returncode != 0 or not os.path.isfile(out):
        print(r.stderr, file=sys.stderr)
        return 1
    print(f"md2pdf: {md_path} -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
