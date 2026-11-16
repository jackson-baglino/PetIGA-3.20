#!/usr/bin/env python3
"""Add one paper to the library, content only (no claims), for the model line.

B1 (the literature-review line) is closed and frozen, so its text-extraction
stage cannot run outside its line. This script reuses that stage's code
unchanged (scripts/literature_review/prepare.py: 300 dpi page images, MinerU2.5-Pro
text and LaTeX, equations, sections, page map, manifest), so a paper added here
is extracted exactly like the nine B1 papers. Decided by the user 2026-11-16.

What it writes:
  library/papers/<paper_id>/source/<copy16>/<file>.pdf   the PDF as given
  library/papers/<paper_id>/content/<copy16>/            the extraction, written once
  library/papers/<paper_id>/card.v1.json                 bibliographic data (from Crossref)
  library/events.jsonl                                   one appended line: who, what, when, hashes
  runs/model-development/ingest/<paper_id>.json           counts-only summary (no paper text)

What it does not do: make claims. Claims for the model come from the model
line's own citations, each checked with library_search.py quote.

Paper text is untrusted. This script prints counts, hashes and paths only.
Run it in your own terminal or through the model line, not in a normal session
(the library is closed to normal sessions).

    venv_workflows/bin/python scripts/model_development/add_paper.py \\
        --doi 10.1103/PhysRevE.84.046707 --pdf <path to pdf> [--paper-id ID] [--dry-run]

The default paper_id is <first author>_<year>_<4 hex of sha256(doi)>, the B1 form.
After adding, rebuild the search index:
    venv_workflows/bin/python scripts/model_development/library_search.py build
"""
import argparse
import datetime
import getpass
import hashlib
import json
import re
import shutil
import sys
import tempfile
import unicodedata
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "literature_review"))


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 16), b""):
            h.update(block)
    return h.hexdigest()


def crossref(doi):
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi)
    with urllib.request.urlopen(url, timeout=20) as r:
        m = json.load(r)["message"]
    return {
        "doi": m.get("DOI"),
        "title": (m.get("title") or [""])[0],
        "authors": [{"family": a.get("family"), "given": a.get("given")} for a in m.get("author", [])],
        "container": (m.get("container-title") or [""])[0],
        "volume": m.get("volume"),
        "article": m.get("article-number") or m.get("page"),
        "year": (m.get("issued", {}).get("date-parts") or [[None]])[0][0],
        "publisher": m.get("publisher"),
        "source": "crossref",
        "fetched": datetime.datetime.now().isoformat(timespec="seconds"),
    }


INCOMING = Path.home() / ("Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/Projects/"
                          "lunar_regolith_DSM/Literature/incoming")


def pdf_from_incoming(review, doi):
    """The PDF zotero_pull.py copied for this DOI (its record: incoming/<review>/zotero_pull.json)."""
    record = INCOMING / review / "zotero_pull.json"
    rec = json.loads(record.read_text())
    items = rec.get("copies", [])  # zotero_pull.py: {"copies": [{"doi", "path" (relative to the review folder), ...}]}
    for it in items:
        if str(it.get("doi", "")).lower() == doi.lower():
            path = Path(it["path"])
            return path if path.is_absolute() else INCOMING / review / path
    sys.exit(f"no PDF recorded for {doi} in {record}")


def default_paper_id(meta):
    family = (meta["authors"][0]["family"] if meta["authors"] else "anon") or "anon"
    family = unicodedata.normalize("NFKD", family).encode("ascii", "ignore").decode().lower()
    family = re.sub(r"[^a-z]", "", family) or "anon"
    tag = hashlib.sha256(meta["doi"].lower().encode()).hexdigest()[:4]
    return f"{family}_{meta['year']}_{tag}"


def main(argv=None):
    ap = argparse.ArgumentParser(description=(__doc__ or "").split("\n\n")[0])
    ap.add_argument("--doi", required=True)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--pdf", help="path to the PDF")
    src.add_argument("--from-incoming", metavar="REVIEW",
                     help="take the PDF that zotero_pull.py recorded for this DOI under Literature/incoming/REVIEW")
    ap.add_argument("--paper-id")
    ap.add_argument("--library", default=str(ROOT / "library"))
    ap.add_argument("--dry-run", action="store_true", help="check inputs and print the plan only")
    a = ap.parse_args(argv)

    pdf = Path(a.pdf) if a.pdf else pdf_from_incoming(a.from_incoming, a.doi)
    with open(pdf, "rb") as f:
        if f.read(5) != b"%PDF-":
            sys.exit(f"{pdf.name} does not start with the PDF signature; refusing")
    copy_hash = sha256(pdf)
    meta = crossref(a.doi)
    if meta["doi"].lower() != a.doi.lower():
        sys.exit(f"Crossref returned DOI {meta['doi']} for {a.doi}; refusing")
    paper_id = a.paper_id or default_paper_id(meta)
    library = Path(a.library)
    paper_dir = library / "papers" / paper_id
    content_final = paper_dir / "content" / copy_hash[:16]
    plan = {"paper_id": paper_id, "doi": meta["doi"], "title": meta["title"], "year": meta["year"],
            "copy_sha256": copy_hash, "content": str(content_final.relative_to(ROOT)),
            "exists": content_final.exists()}
    print(json.dumps({"plan": plan}, indent=1, ensure_ascii=False))
    if a.dry_run:
        return 0
    if content_final.exists():
        sys.exit("content for this copy exists; content is written once")

    import prepare as p1  # the B1 extraction code, used unchanged

    with tempfile.TemporaryDirectory(dir=ROOT / "runs") as tmp:
        summary = p1.prepare(tmp, paper_id, str(pdf), parser="mineru")
        produced = Path(tmp) / "untrusted" / "content" / paper_id / copy_hash[:16]
        content_final.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(produced), str(content_final))

    source_dir = paper_dir / "source" / copy_hash[:16]
    source_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(pdf, source_dir / pdf.name)

    card_path = paper_dir / "card.v1.json"
    if card_path.exists():
        card = json.loads(card_path.read_text())
    else:
        card = {"paper_id": paper_id, "identifiers": {"doi": meta["doi"]}, "bibliographic": meta,
                "copies": [], "status": "content_only", "branches": [], "tags": [],
                "added_by": "model-development add_paper.py (content only, no claims)"}
    card["copies"].append({"copy_sha256": copy_hash, "file": pdf.name, "kind": "published PDF",
                           "content": str(content_final.relative_to(library)),
                           "manifest_sha256": sha256(content_final / "manifest.json")})
    card_path.write_text(json.dumps(card, indent=1, ensure_ascii=False) + "\n")

    event = {"time": datetime.datetime.now().isoformat(timespec="seconds"), "who": getpass.getuser(),
             "what": "add_paper content_only", "paper_id": paper_id, "doi": meta["doi"],
             "copy_sha256": copy_hash, "card_sha256": sha256(card_path),
             "manifest_sha256": sha256(content_final / "manifest.json"),
             "tool": "scripts/model_development/add_paper.py"}
    with open(library / "events.jsonl", "a", encoding="utf-8") as f:
        f.write(json.dumps(event) + "\n")

    out = ROOT / "runs" / "model-development" / "ingest"
    out.mkdir(parents=True, exist_ok=True)
    keep = ["pages", "equations", "numbered_equations", "headings", "failed_pages",
            "pages_ink_not_fully_covered", "page_map_offset", "seconds"]
    report = {"paper_id": paper_id, "doi": meta["doi"], "copy_sha256": copy_hash,
              **{k: summary.get(k) for k in keep}, "event": event}
    (out / f"{paper_id}.json").write_text(json.dumps(report, indent=1) + "\n")
    print(json.dumps(report, indent=1))
    print("next: venv_workflows/bin/python scripts/model_development/library_search.py build")
    return 1 if summary.get("failed_pages") else 0


if __name__ == "__main__":
    sys.exit(main())
