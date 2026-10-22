#!/usr/bin/env python3
"""S1 of units-verification: render the page images the agents will read.

    python scripts/units_verification/render_pages.py --run-id <id> [--add SOURCE:PAGE ...] [--dpi 200]

For every literature row of the units table, renders the cited page and the pages before and after it
(`mutool draw -r 200`) to runs/units-verification/<id>/pages/<source>_p<N>.png, using the PDFs that S0
recorded. --add renders further pages (S3's re-read requests). A page rendered before is rendered again
and must give the same sha256 (a determinism check); a mismatch stops the stage.
Pages are not committed (paper material); their hashes are, in the S1 entry of manifest.json.
"""
import argparse
import sys
import tempfile
from pathlib import Path

import uvlib


def wanted_pages(rows, n_pages):
    """{(source, page)} for the cited page of each row and its neighbours, inside the document."""
    want = set()
    for r in rows:
        if r["page"] is None:
            continue
        for p in (r["page"] - 1, r["page"], r["page"] + 1):
            if 1 <= p <= n_pages[r["source"]]:
                want.add((r["source"], p))
    return want


def page_name(source, page):
    return f"{source}_p{page}.png"


def render(pdf, page, out, dpi):
    uvlib.run(["mutool", "draw", "-q", "-r", str(dpi), "-o", str(out), str(pdf), str(page)])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--table", type=Path, default=uvlib.UNITS_MD)
    ap.add_argument("--dpi", type=int, default=200)
    ap.add_argument("--add", action="append", default=[], metavar="SOURCE:PAGE")
    ap.add_argument("--allow-dirty", action="store_true", help="development only: run with uncommitted definitions")
    a = ap.parse_args(argv)
    rdir = uvlib.run_dir(a.run_id)
    print(f"S1 render_pages: run={uvlib.rel(rdir)} table={uvlib.rel(a.table)} dpi={a.dpi} add={a.add}")
    try:
        uvlib.require_committed(a.allow_dirty)
        manifest = uvlib.load_manifest(rdir)
        s0 = manifest.get("stages", {}).get("S0")
        if not s0 or s0.get("verdict") != "pass":
            raise uvlib.StageError("S0 has not passed for this run; run locate_sources.py first")
        sources = s0["sources"]
        old = manifest["stages"].get("S1", {})
        if old and old.get("dpi") != a.dpi:
            raise uvlib.StageError(f"this run's pages were rendered at {old.get('dpi')} dpi, not {a.dpi}")
        want = wanted_pages(uvlib.literature_rows(uvlib.read_table(a.table)),
                            {s: v["pages"] for s, v in sources.items()})
        for item in a.add:
            src, _, pg = item.partition(":")
            if src not in sources or not pg.isdigit() or not 1 <= int(pg) <= sources[src]["pages"]:
                raise uvlib.StageError(f"--add {item}: unknown source or page out of range")
            want.add((src, int(pg)))
        for src, v in sources.items():
            if uvlib.sha256(uvlib.abspath(v["path"])) != v["sha256"]:
                raise uvlib.StageError(f"{v['path']} changed since S0 (sha256 differs)")
        pages = dict(old.get("pages", {}))
        pdir = rdir / "pages"
        pdir.mkdir(exist_ok=True)
        new, rechecked = 0, 0
        for src, pg in sorted(want):
            name = page_name(src, pg)
            with tempfile.TemporaryDirectory() as tmp:
                out = Path(tmp) / name
                render(uvlib.abspath(sources[src]["path"]), pg, out, a.dpi)
                h = uvlib.sha256(out)
                if name in pages and pages[name] != h:
                    raise uvlib.StageError(f"{name}: re-render gives a different sha256 (not deterministic)")
                if (pdir / name).exists() and uvlib.sha256(pdir / name) != h:
                    raise uvlib.StageError(f"{name}: file on disk differs from a fresh render")
                rechecked += name in pages
                new += name not in pages
                out.replace(pdir / name)
                pages[name] = h
    except uvlib.StageError as e:
        print(f"STOP: {e}", file=sys.stderr)
        return 2
    uvlib.update_manifest(rdir, "S1", {"dpi": a.dpi, "renderer": "mutool draw", "pages": dict(sorted(pages.items())),
                                       "added": sorted(set(old.get("added", [])) | set(a.add)),
                                       "verdict": "pass"}, a.allow_dirty)
    print(f"rendered {len(want)} pages ({new} new, {rechecked} re-rendered with identical sha256) "
          f"into {uvlib.rel(pdir)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
