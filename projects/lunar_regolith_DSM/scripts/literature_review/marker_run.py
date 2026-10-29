#!/usr/bin/env python3
"""Run Marker once on a PDF and write its Markdown and its page layout (P1, called by prepare.py).

    marker_run.py <pdf> <outdir>

Writes <outdir>/marker.md (paginated Markdown, as marker_single writes it with the flags in
prepare.MARKER_FLAGS) and <outdir>/layout.json: for each page, the page box in PDF points
and Marker's top-level blocks in reading order, each with its type (Text, Table, Figure,
Caption, Equation, ...) and box. The document is built once and rendered twice, so the
layout is the one the text came from. Never uses an LLM service.
"""
import json
import sys
from pathlib import Path


def main(pdf, outdir):
    from marker.config.parser import ConfigParser
    from marker.converters.pdf import PdfConverter
    from marker.models import create_model_dict
    from marker.renderers.json import JSONRenderer
    from marker.renderers.markdown import MarkdownRenderer

    cp = ConfigParser({"paginate_output": True, "output_format": "markdown",
                       "disable_image_extraction": True, "disable_tqdm": True})
    converter = PdfConverter(config=cp.generate_config_dict(), artifact_dict=create_model_dict(),
                             processor_list=cp.get_processors(), renderer=cp.get_renderer())
    document = converter.build_document(pdf)
    md = converter.resolve_dependencies(MarkdownRenderer)(document)
    js = converter.resolve_dependencies(JSONRenderer)(document)
    out = Path(outdir)
    (out / "marker.md").write_text(md.markdown)
    pages = [{"bbox": p.bbox, "blocks": [{"type": b.block_type, "bbox": b.bbox} for b in (p.children or [])]}
             for p in js.children]
    (out / "layout.json").write_text(json.dumps(pages, indent=1))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
