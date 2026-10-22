#!/usr/bin/env python3
"""S1b of units-verification: a second, independent reading of every page with Marker (local OCR, no LLM).

    python scripts/units_verification/marker_extract.py --run-id <id> [--sources ID ...]

Runs `marker_single <pdf> --force_ocr --paginate_output --disable_image_extraction --output_format markdown`
(Marker 1.10.2 from venv_workflows; never --use_llm, which would send pages to an external API) on each
source that S0 recorded, splits the output into runs/units-verification/<id>/marker/<source>_p<N>.md
(N = PDF page, 1-based) and builds marker_index.json: for each literature row, the pages of its source on
which its printed symbol appears. A page Marker returns nothing for is logged and has no Marker reading.
Page Markdown is paper text and is not committed; its sha256 is, in manifest.json, with the index.
"""
import argparse
import importlib.metadata
import re
import sys
import tempfile
import time
from pathlib import Path

import unitparse
import uvlib

MARKER_FLAGS = ["--force_ocr", "--paginate_output", "--disable_image_extraction", "--output_format", "markdown",
                "--disable_tqdm"]
PAGE_SEP = re.compile(r"^\{(\d+)\}-{48}\s*$", flags=re.M)
MODEL_CACHE = Path.home() / "Library" / "Caches" / "datalab" / "models"


def split_pages(md):
    """{page (1-based): markdown} from Marker's paginated output ({0}----... separators, 0-based)."""
    parts = PAGE_SEP.split(md)
    return {int(parts[i]) + 1: parts[i + 1].strip("\n") + "\n" for i in range(1, len(parts) - 1, 2)}


def marker_versions():
    v = {p: importlib.metadata.version(p) for p in ("marker-pdf", "surya-ocr", "torch")}
    v["model_cache"] = sorted(str(p.relative_to(MODEL_CACHE)) for p in MODEL_CACHE.glob("*/*")) \
        if MODEL_CACHE.exists() else None
    return v


def run_marker(pdf, n_pages, outdir, source):
    marker = Path(sys.executable).parent / "marker_single"
    with tempfile.TemporaryDirectory() as tmp:
        t0 = time.time()
        uvlib.run([str(marker), str(pdf), "--output_dir", tmp, *MARKER_FLAGS])
        mds = list(Path(tmp).rglob("*.md"))
        if len(mds) != 1:
            raise uvlib.StageError(f"{source}: Marker wrote {len(mds)} Markdown files, expected 1")
        pages = split_pages(mds[0].read_text())
        seconds = round(time.time() - t0, 1)
    hashes, failed = {}, []
    for p in range(1, n_pages + 1):
        text = pages.get(p, "")
        if not text.strip():
            failed.append(p)
            continue
        out = outdir / f"{source}_p{p}.md"
        out.write_text(text)
        hashes[out.name] = uvlib.sha256(out)
    return hashes, failed, seconds


def build_index(rows, mdir, n_pages):
    """{symbol: {source, printed, cited_page, pages}}: pages whose Marker text contains the printed symbol."""
    index = {}
    for r in rows:
        hits = []
        if r["printed"]:
            for p in range(1, n_pages[r["source"]] + 1):
                f = mdir / f"{r['source']}_p{p}.md"
                if f.exists() and unitparse.find_symbol(unitparse.flatten_markdown(f.read_text()), r["printed"]):
                    hits.append(p)
        index[r["symbol"]] = {"source": r["source"], "printed": r["printed"], "cited_page": r["page"], "pages": hits}
    return index


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--table", type=Path, default=uvlib.UNITS_MD)
    ap.add_argument("--sources", nargs="*", help="only these sources (default: all that S0 recorded)")
    ap.add_argument("--index-only", action="store_true", help="rebuild marker_index.json from existing pages")
    ap.add_argument("--allow-dirty", action="store_true", help="development only: run with uncommitted definitions")
    a = ap.parse_args(argv)
    rdir = uvlib.run_dir(a.run_id)
    mdir = rdir / "marker"
    print(f"S1b marker_extract: run={uvlib.rel(rdir)} table={uvlib.rel(a.table)} flags={' '.join(MARKER_FLAGS)}")
    try:
        uvlib.require_committed(a.allow_dirty)
        manifest = uvlib.load_manifest(rdir)
        s0 = manifest.get("stages", {}).get("S0")
        if not s0 or s0.get("verdict") != "pass":
            raise uvlib.StageError("S0 has not passed for this run; run locate_sources.py first")
        sources = s0["sources"]
        old = manifest["stages"].get("S1b", {})
        record = {"flags": MARKER_FLAGS, "versions": marker_versions(), "pages": dict(old.get("pages", {})),
                  "failed_pages": dict(old.get("failed_pages", {})), "seconds": dict(old.get("seconds", {}))}
        mdir.mkdir(exist_ok=True)
        todo = [] if a.index_only else (a.sources or list(sources))
        for src in todo:
            v = sources[src]
            if uvlib.sha256(uvlib.abspath(v["path"])) != v["sha256"]:
                raise uvlib.StageError(f"{v['path']} changed since S0 (sha256 differs)")
            print(f"  {src}: Marker on {v['pages']} pages ...", flush=True)
            hashes, failed, seconds = run_marker(uvlib.abspath(v["path"]), v["pages"], mdir, src)
            record["pages"] = {k: h for k, h in record["pages"].items() if not k.startswith(src + "_p")} | hashes
            record["failed_pages"][src], record["seconds"][src] = failed, seconds
            print(f"  {src}: {len(hashes)} pages read, failed pages {failed or 'none'}, {seconds} s")
        rows = uvlib.literature_rows(uvlib.read_table(a.table))
        index = build_index(rows, mdir, {s: v["pages"] for s, v in sources.items()})
        uvlib.write_json(rdir / "marker_index.json", index)
    except uvlib.StageError as e:
        print(f"STOP: {e}", file=sys.stderr)
        return 2
    located = sum(1 for v in index.values() if v["pages"])
    on_cited = sum(1 for v in index.values() if v["cited_page"] in v["pages"])
    record.update({"index": "marker_index.json", "index_sha256": uvlib.sha256(rdir / "marker_index.json"),
                   "symbols_located": located, "symbols_on_cited_page": on_cited, "symbols": len(index),
                   "verdict": "pass"})
    record["pages"] = dict(sorted(record["pages"].items()))
    uvlib.update_manifest(rdir, "S1b", record, a.allow_dirty)
    print(f"index: {located}/{len(index)} printed symbols found by Marker somewhere in their source, "
          f"{on_cited}/{len(index)} on the cited page -> {uvlib.rel(rdir / 'marker_index.json')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
