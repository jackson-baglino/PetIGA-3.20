#!/usr/bin/env python3
"""S0 of units-verification: find the PDF of each source and record it in the run manifest.

    python scripts/units_verification/locate_sources.py --run-id <id> [--literature DIR] [--sources FILE]

Every PDF under the literature folder is a candidate for a source when its page-1 text contains the
source's DOI (or, if the source has no DOI on page 1, its title). Page 1 must also name the paper (title).
Several candidates for one source are one source only if every page renders to identical pixels; then
'prefer' in sources.json picks the version of record. Otherwise the stage stops and asks the user.
Writes the S0 entry of runs/units-verification/<id>/manifest.json (paths shown with ~ for the home folder).
"""
import argparse
import json
import re
import sys
import tempfile
from pathlib import Path

import uvlib


def norm(s):
    """Lower case, hyphens and line breaks dropped, single spaces: for matching titles and DOIs in page text."""
    s = s.lower().replace("‐", "-").replace("‑", "-")
    s = re.sub(r"-\s*\n\s*", "", s)          # words hyphenated across lines
    s = re.sub(r"[-\s]+", " ", s)
    return s.strip()


def page_text(pdf, page=1):
    return uvlib.run(["mutool", "draw", "-q", "-F", "txt", "-o", "-", str(pdf), str(page)])


def page_count(pdf):
    out = uvlib.run(["mutool", "info", str(pdf)])
    m = re.search(r"^Pages:\s*(\d+)", out, flags=re.M)
    if not m:
        raise uvlib.StageError(f"cannot read the page count of {pdf}")
    return int(m.group(1))


def render_hashes(pdf, dpi=50):
    """sha256 of each page rendered to raw grey pixels (PGM): equal lists mean the pages look identical."""
    with tempfile.TemporaryDirectory() as tmp:
        uvlib.run(["mutool", "draw", "-q", "-r", str(dpi), "-c", "gray", "-o", f"{tmp}/p%04d.pgm", str(pdf)])
        return [uvlib.sha256(p) for p in sorted(Path(tmp).glob("p*.pgm"))]


def matches(text, spec):
    t = norm(text)
    doi_hit = bool(spec.get("doi")) and norm(spec["doi"]) in t
    title_hit = norm(spec["title"]) in t
    return (doi_hit if spec.get("doi") else title_hit), title_hit


def locate(literature, sources):
    pdfs = sorted(p for p in Path(literature).rglob("*") if p.suffix.lower() == ".pdf" and p.is_file())
    found = {sid: [] for sid in sources}
    unreadable = []
    for pdf in pdfs:
        try:
            text = page_text(pdf)
        except uvlib.StageError as e:
            unreadable.append({"path": uvlib.rel(pdf), "error": str(e).splitlines()[0]})
            continue
        for sid, spec in sources.items():
            hit, title_hit = matches(text, spec)
            if hit:
                found[sid].append({"path": pdf, "title_on_page1": title_hit})
    result, problems = {}, []
    for sid, cands in found.items():
        spec = sources[sid]
        if not cands:
            problems.append(f"{sid}: no PDF under {uvlib.rel(literature)} has its "
                            f"{'DOI ' + spec['doi'] if spec.get('doi') else 'title'} on page 1")
            continue
        for c in cands:
            if not c["title_on_page1"]:
                problems.append(f"{sid}: {uvlib.rel(c['path'])} has the DOI but page 1 does not show the title")
        groups = {}
        for c in cands:
            c["render_sha256"] = render_hashes(c["path"])
            groups.setdefault(tuple(c["render_sha256"]), []).append(c)
        if len(groups) > 1:
            problems.append(f"{sid}: {len(cands)} candidate files render differently: "
                            + ", ".join(uvlib.rel(c["path"]) for c in cands))
            continue
        pick = cands[0]
        if len(cands) > 1:
            prefer = spec.get("prefer")
            pick = next((c for c in cands if prefer and str(c["path"]).endswith(prefer)), None)
            if pick is None:
                problems.append(f"{sid}: {len(cands)} pixel-identical copies and no 'prefer' entry naming one")
                continue
        result[sid] = {
            "path": uvlib.rel(pick["path"]),
            "sha256": uvlib.sha256(pick["path"]),
            "pages": page_count(pick["path"]),
            "doi": spec.get("doi"),
            "title": spec["title"],
            "title_on_page1": pick["title_on_page1"],
            "identical_copies": [{"path": uvlib.rel(c["path"]), "sha256": uvlib.sha256(c["path"])}
                                 for c in cands if c is not pick],
        }
    return result, problems, len(pdfs), unreadable


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--literature", type=Path, default=uvlib.LITERATURE)
    ap.add_argument("--sources", type=Path, default=uvlib.WORKFLOW_DIR / "sources.json")
    ap.add_argument("--table", type=Path, default=uvlib.UNITS_MD)
    ap.add_argument("--allow-dirty", action="store_true", help="development only: run with uncommitted definitions")
    a = ap.parse_args(argv)
    print(f"S0 locate_sources: literature={uvlib.rel(a.literature)} sources={uvlib.rel(a.sources)} "
          f"table={uvlib.rel(a.table)} run_id={a.run_id}")
    try:
        uvlib.require_committed(a.allow_dirty)
        spec = {k: v for k, v in json.loads(a.sources.read_text()).items() if not k.startswith("_")}
        needed = sorted({r["source"] for r in uvlib.literature_rows(uvlib.read_table(a.table))})
        missing_spec = [s for s in needed if s not in spec]
        if missing_spec:
            raise uvlib.StageError(f"sources.json has no entry for {missing_spec}")
        result, problems, n, unreadable = locate(a.literature, {s: spec[s] for s in needed})
    except uvlib.StageError as e:
        print(f"STOP: {e}", file=sys.stderr)
        return 2
    rdir = uvlib.run_dir(a.run_id)
    uvlib.update_manifest(rdir, "S0", {
        "inputs": {"literature": uvlib.rel(a.literature), "sources_spec": uvlib.rel(a.sources),
                   "sources_spec_sha256": uvlib.sha256(a.sources), "table": uvlib.rel(a.table),
                   "table_sha256": uvlib.sha256(a.table)},
        "pdfs_scanned": n, "unreadable_pdfs": unreadable,
        "sources": result, "problems": problems, "verdict": "pass" if not problems else "stop",
    }, a.allow_dirty)
    print(f"scanned {n} PDFs ({len(unreadable)} unreadable)")
    for sid, r in result.items():
        extra = f"; identical copies: {[c['path'] for c in r['identical_copies']]}" if r["identical_copies"] else ""
        print(f"  {sid}: {r['path']} ({r['pages']} pages, sha256 {r['sha256'][:12]}){extra}")
    if problems:
        print("STOP, ask the user:\n  " + "\n  ".join(problems), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
