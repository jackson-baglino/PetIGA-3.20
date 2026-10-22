#!/usr/bin/env python3
"""S2 of units-verification: blind reads of each page by morty-unit-reader.

    python scripts/units_verification/read_pages.py --run-id <id> [--rereads] [--only SOURCE:PAGE ...]
                                                    [--repeat N --out-dir reads_repeat]

For every page cited in the units table, the reader gets the page image and the symbols cited there (id and
printed form; never the unit, value or meaning) and returns, per symbol, what the page prints.
--rereads reads instead the pages S3 asked for (rows NOT_FOUND on their cited page; compare.json).
--repeat N reads the first N cited pages again into --out-dir, for the reader-stability check.
Writes reads/<source>_p<N>.json (a later call on the same page adds its symbols) and the S2 entry of the manifest.
"""
import argparse
import json
import sys
from pathlib import Path

import agents
import uvlib

AGENT = "morty-unit-reader"
NULLABLE = {"type": ["string", "null"]}
SCHEMA = {
    "type": "object", "additionalProperties": False, "required": ["reads"],
    "properties": {"reads": {"type": "array", "items": {
        "type": "object", "additionalProperties": False,
        "required": ["symbol", "found", "printed_as", "meaning", "value_as_printed", "unit_as_printed", "location"],
        "properties": {"symbol": {"type": "string"}, "found": {"type": "boolean"}, "printed_as": NULLABLE,
                       "meaning": NULLABLE, "value_as_printed": NULLABLE, "unit_as_printed": NULLABLE,
                       "location": NULLABLE}}}}}


def task_text(symbols):
    """The reader's task: the image name and the symbols to look for (id and printed form only)."""
    lines = [{"id": s["symbol"], "look_for": s["printed"]} for s in symbols]
    return ("Page image: page.png (one page of a scientific paper).\n"
            "Report what this page prints for each of these symbols, one entry per symbol, in this order:\n"
            + json.dumps(lines, ensure_ascii=False, indent=1))


def answers_every_symbol(symbols):
    want = [s["symbol"] for s in symbols]

    def check(out):
        got = [r["symbol"] for r in out["reads"]]
        return [] if sorted(got) == sorted(want) else [f"answered {got}, asked {want}"]
    return check


def read_page(rdir, stage, source, page, symbols, out_dir="reads"):
    """One reader call for one page; merges into <out_dir>/<source>_p<page>.json. Returns the file path."""
    rdir = Path(rdir)
    img = rdir / "pages" / f"{source}_p{page}.png"
    if not img.exists():
        raise uvlib.StageError(f"{img.name} is not rendered; run render_pages.py --add {source}:{page}")
    out_file = rdir / out_dir / f"{source}_p{page}.json"
    old = json.loads(out_file.read_text()) if out_file.exists() else {"page": f"{source}_p{page}", "reads": [],
                                                                       "calls": []}
    done = {r["symbol"] for r in old["reads"]}
    symbols = [s for s in symbols if s["symbol"] not in done]
    if not symbols:
        return out_file, None
    label = f"{source}_p{page}_" + uvlib.sha256_bytes(",".join(s["symbol"] for s in symbols).encode())[:8]
    out, info = agents.call_agent(rdir, stage, AGENT, label, task_text(symbols), {"page.png": img}, SCHEMA,
                                  check=answers_every_symbol(symbols))
    old["reads"] += out["reads"]
    old["calls"].append({"symbols": [s["symbol"] for s in symbols], "page_png_sha256": uvlib.sha256(img), **info})
    uvlib.write_json(out_file, old)
    return out_file, info


def cited_groups(rows):
    groups = {}
    for r in rows:
        groups.setdefault((r["source"], r["page"]), []).append(r)
    return groups


def reread_groups(rdir, rows):
    comp = json.loads((Path(rdir) / "compare.json").read_text())
    by_sym = {r["symbol"]: r for r in rows}
    groups = {}
    for c in comp["rows"]:
        for p in c.get("reread_pages", []):
            groups.setdefault((c["source"], p), []).append(by_sym[c["symbol"]])
    return groups


def run_groups(rdir, stage, groups, out_dir="reads"):
    done, failed = [], []
    for (src, pg), syms in sorted(groups.items()):
        try:
            f, info = read_page(rdir, stage, src, pg, syms, out_dir)
            done.append({"page": f"{src}_p{pg}", "symbols": [s["symbol"] for s in syms],
                         "file": f.name, "file_sha256": uvlib.sha256(f), "called": info is not None})
            print(f"  {src} p{pg}: {len(syms)} symbols" + ("" if info else " (already read)"), flush=True)
        except agents.AgentError as e:
            failed.append({"page": f"{src}_p{pg}", "error": json.loads(str(e))})
            print(f"  {src} p{pg}: FAILED after retry, logged", flush=True)
    return done, failed


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--table", type=Path, default=uvlib.UNITS_MD)
    ap.add_argument("--rereads", action="store_true")
    ap.add_argument("--only", nargs="*", default=None, metavar="SOURCE:PAGE")
    ap.add_argument("--repeat", type=int, default=0)
    ap.add_argument("--out-dir", default="reads")
    ap.add_argument("--stage", default="S2")
    ap.add_argument("--allow-dirty", action="store_true", help="development only: run with uncommitted definitions")
    a = ap.parse_args(argv)
    rdir = uvlib.run_dir(a.run_id)
    print(f"{a.stage} read_pages: run={uvlib.rel(rdir)} table={uvlib.rel(a.table)} rereads={a.rereads} "
          f"only={a.only} repeat={a.repeat} out={a.out_dir}")
    try:
        uvlib.require_committed(a.allow_dirty)
        rows = uvlib.literature_rows(uvlib.read_table(a.table))
        groups = reread_groups(rdir, rows) if a.rereads else cited_groups(rows)
        if a.only is not None:
            keep = {(s, int(p)) for s, _, p in (x.partition(":") for x in a.only)}
            groups = {k: v for k, v in groups.items() if k in keep}
        if a.repeat:
            groups = dict(sorted(groups.items())[:a.repeat])
        done, failed = run_groups(rdir, a.stage, groups, a.out_dir)
    except uvlib.StageError as e:
        print(f"STOP: {e}", file=sys.stderr)
        return 2
    m = uvlib.load_manifest(rdir).get("stages", {}).get(a.stage, {})
    agent = agents.load_agent(AGENT)
    uvlib.update_manifest(rdir, a.stage, {
        "agent": AGENT, "model": agent["model"], "agent_sha256": agent["sha256"],
        "table": uvlib.rel(a.table), "table_sha256": uvlib.sha256(a.table),
        "pages": m.get("pages", []) + done, "failed": m.get("failed", []) + failed,
        "verdict": "pass" if not failed else "failures logged"}, a.allow_dirty)
    print(f"{len(done)} pages read, {len(failed)} failed")
    return 0 if not failed else 1


if __name__ == "__main__":
    sys.exit(main())
