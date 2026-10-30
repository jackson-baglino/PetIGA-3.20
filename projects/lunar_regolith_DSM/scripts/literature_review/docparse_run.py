#!/usr/bin/env python3
"""Run MinerU2.5-Pro on a paper's page images and write its text and layout (P1, called by prepare.py).

    venv_docparse/bin/python docparse_run.py <images dir> <outdir> <dpi>

MinerU2.5-Pro (opendatalab/MinerU2.5-Pro-2605-1.2B, Apache-2.0) runs in this process
through MLX on Apple Silicon (mineru-vl-utils, backend "mlx-engine"): a layout pass finds
the blocks of each page, and a second pass reads each block (LaTeX for equations, HTML
for tables, text otherwise). The weights come from the local Hugging Face cache; the run
is offline and nothing leaves the machine.

Writes to <outdir>:
  pages.json    {page: Markdown}; displayed equations as $$...$$ with \\tag{n} for the
                printed number, so prepare.equations reads them
  layout.json   per page: the page box in PDF points and the blocks in reading order,
                each with a type in prepare's names (Text, Table, Figure, Caption,
                Equation, SectionHeader, PageHeader, PageFooter, ...) and a box in points;
                an equation's printed number is merged into the equation's box
  tool.json     versions, model snapshot and file sizes, seconds per page
Prints nothing from the pages.
"""
import importlib.metadata
import json
import re
import sys
import time
from pathlib import Path

MODEL_REPO = "opendatalab/MinerU2.5-Pro-2605-1.2B"

TYPES = {                                   # MinerU block type -> prepare's layout type
    "text": "Text", "index": "Text", "phonetic": "Text", "unknown": "Text", "list_item": "ListItem",
    "list": "ListGroup", "ref_text": "Reference", "code": "Code", "algorithm": "Code",
    "title": "SectionHeader", "doc_title": "SectionHeader", "paragraph_title": "SectionHeader",
    "table": "Table", "image": "Figure", "chart": "Figure", "image_block": "FigureGroup",
    "equation": "Equation", "equation_block": "Equation", "formula_number": "Equation",
    "table_caption": "Caption", "image_caption": "Caption", "code_caption": "Caption", "caption": "Caption",
    "table_footnote": "Caption", "image_footnote": "Caption", "footnote": "Caption",
    "header": "PageHeader", "footer": "PageFooter", "page_number": "PageFooter",
    "page_footnote": "Footnote", "aside_text": "PageHeader",
}
TAG = re.compile(r"\\tag\*?\s*\{\s*([^}]*?)\s*\}")


def model_path():
    from huggingface_hub import snapshot_download
    return snapshot_download(MODEL_REPO, local_files_only=True)


def equation_latex(content):
    """MinerU writes \\[ ... \\tag {n} \\]; return the LaTeX inside and the number."""
    s = (content or "").strip()
    s = re.sub(r"^\\\[|\\\]$", "", s).strip()
    m = TAG.search(s)
    number = m.group(1).strip("() ") if m else None
    return TAG.sub("", s).strip(), number


def attach_numbers(blocks):
    """Merge each formula_number block into the equation level with it; returns blocks."""
    eqs = [b for b in blocks if b["type"] in ("equation", "equation_block")]
    out = []
    for b in blocks:
        if b["type"] != "formula_number":
            out.append(b)
            continue
        y = (b["bbox"][1] + b["bbox"][3]) / 2
        level = [e for e in eqs if e["bbox"][1] <= y <= e["bbox"][3]]
        num = re.sub(r"[()\s]", "", b.get("content") or "")
        if level and num:
            e = min(level, key=lambda e: abs(b["bbox"][0] - e["bbox"][2]))   # same column
            e["bbox"] = [min(e["bbox"][0], b["bbox"][0]), min(e["bbox"][1], b["bbox"][1]),
                         max(e["bbox"][2], b["bbox"][2]), max(e["bbox"][3], b["bbox"][3])]
            if not TAG.search(e.get("content") or ""):
                e["number"] = num
        else:
            out.append(b)                    # a number with no equation beside it stays a block
    return out


def markdown(blocks):
    parts = []
    for b in blocks:
        t, c = b["type"], (b.get("content") or "").strip()
        if t in ("equation", "equation_block"):
            latex, number = equation_latex(c)
            number = number or b.get("number")
            parts.append(f"$${latex}" + (f" \\tag{{{number}}}" if number else "") + "$$")
        elif TYPES.get(t) == "SectionHeader":
            parts.append(("# " if t == "doc_title" else "## ") + c)
        elif t in ("image", "chart", "image_block"):
            parts.append("[figure]")
        elif c:
            parts.append(c)
    return "\n\n".join(parts) + "\n"


def versions(path, seconds):
    v = {p: importlib.metadata.version(p) for p in ("mineru-vl-utils", "mlx", "mlx-vlm", "transformers")}
    v["model"] = {"repo": MODEL_REPO, "snapshot": Path(path).name,
                  "files": {f.name: f.stat().st_size for f in sorted(Path(path).iterdir()) if f.is_file()}}
    v["seconds_per_page"] = seconds
    return v


def main(imgdir, outdir, dpi):
    from PIL import Image
    from mineru_vl_utils import MinerUClient

    path = model_path()
    client = MinerUClient(backend="mlx-engine", model_path=path, use_tqdm=False)
    pages, layout, seconds = {}, [], {}
    for f in sorted(Path(imgdir).glob("p*.png")):
        p = int(f.stem[1:])
        t0 = time.time()
        img = Image.open(f).convert("RGB")
        blocks = attach_numbers([dict(b) for b in client.two_step_extract(img)])
        w, h = img.size[0] * 72 / dpi, img.size[1] * 72 / dpi            # page in points
        pages[p] = markdown(blocks)
        layout.append({"bbox": [0, 0, w, h],
                       "blocks": [{"type": TYPES.get(b["type"], "Text"),
                                   "bbox": [b["bbox"][0] * w, b["bbox"][1] * h, b["bbox"][2] * w, b["bbox"][3] * h]}
                                  for b in blocks]})
        seconds[p] = round(time.time() - t0, 1)
    out = Path(outdir)
    (out / "pages.json").write_text(json.dumps(pages))
    (out / "layout.json").write_text(json.dumps(layout))
    (out / "tool.json").write_text(json.dumps(versions(path, seconds), indent=1))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], float(sys.argv[3]))
