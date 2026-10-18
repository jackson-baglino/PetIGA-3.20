#!/usr/bin/env python3
"""Turn a PDF into page-marked text plus a quality manifest, and render pages that need visual checking.

Uses mutool (same text layer the readers see). Writes <out>/<id>.txt ('=== PAGE n ===' markers) and <out>/<id>.pages.json:
  pages: n, chars, flags (garbled_glyphs, math_dense, table_like, short), render: path of a PNG when a page is flagged.
The flags do not fix the text; they route numbers and equations on those pages to an image check.

Usage: prepare_text.py --pdf FILE --id ID --out DIR [--render-dpi 90]
"""
import argparse, json, re, subprocess, sys
from pathlib import Path

PUA = re.compile(r"[-]|\(cid:\d+\)|/H\d{3,5}")
MATHCH = re.compile(r"[=∂∇∑∫√±×≈≤≥αβγδεζηθλμνξπρσφψω]")


def page_count(pdf):
    out = subprocess.run(["mutool", "info", str(pdf)], capture_output=True, text=True).stdout
    m = re.search(r"Pages:\s*(\d+)", out)
    return int(m.group(1)) if m else 0


def page_text(pdf, n):
    return subprocess.run(["mutool", "draw", "-q", "-F", "txt", "-o", "-", str(pdf), str(n)], capture_output=True, text=True).stdout


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdf", required=True); ap.add_argument("--id", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--render-dpi", type=int, default=90)
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    n = page_count(a.pdf)
    if n == 0:
        sys.exit(f"cannot read {a.pdf}")
    chunks, pages = [], []
    for k in range(1, n + 1):
        t = page_text(a.pdf, k)
        chars = len(t)
        flags = []
        if PUA.search(t): flags.append("garbled_glyphs")
        if chars and len(MATHCH.findall(t)) / chars > 0.012: flags.append("math_dense")
        if len(re.findall(r"\d", t)) / max(chars, 1) > 0.25: flags.append("table_like")
        if chars < 400: flags.append("short")
        rec = {"page": k, "chars": chars, "flags": flags}
        if {"garbled_glyphs", "math_dense", "table_like"} & set(flags):
            png = out / f"{a.id}_p{k}.png"
            subprocess.run(["mutool", "draw", "-q", "-r", str(a.render_dpi), "-o", str(png), str(a.pdf), str(k)], capture_output=True)
            rec["render"] = str(png)
        pages.append(rec)
        chunks.append(f"=== PAGE {k} ===\n{t}")
    (out / f"{a.id}.txt").write_text("\n".join(chunks))
    (out / f"{a.id}.pages.json").write_text(json.dumps({"id": a.id, "pdf": str(a.pdf), "n_pages": n, "pages": pages}, indent=1))
    flagged = sum(1 for p in pages if "render" in p)
    print(f"{a.id}: {n} pages, {sum(p['chars'] for p in pages):,} chars, {flagged} pages flagged for image check")


if __name__ == "__main__":
    main()
