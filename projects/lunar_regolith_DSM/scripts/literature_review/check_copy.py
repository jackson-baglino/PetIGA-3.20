#!/usr/bin/env python3
"""A2: is this copy the version of record? (WORKFLOW.md §3, A2; B1_CONTRACT.md, "Copy kinds").

    check_copy.py --run <run folder> --paper-id <id> --pdf <file>

A script, no agent. From the PDF's text layer (mutool) it records signals on the first two pages:
- the paper's DOI from the contract, found or not;
- an arXiv stamp (arXiv:NNNN.NNNNN);
- a manuscript notice (accepted manuscript, author manuscript, HHS Public Access, preprint, submitted to);
- a publisher line (copyright sign or the publisher's name).
Kind: arxiv if an arXiv stamp; accepted_manuscript if a manuscript notice; publisher if the DOI and a publisher
line are both found and neither of the others; else unknown. Anything but publisher, or a DOI not found, opens a
blocking A2 gate for this paper: the user decides whether this copy is read (the line waits). The page text goes to
<run>/untrusted/<paper>/copy/; the summary <run>/summary_copy_<paper>.json holds only the signals and the kind.
"""
import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_b1  # noqa: E402
from open_gate import open_gate  # noqa: E402

ARXIV = re.compile(r"arXiv:\s?\d{4}\.\d{4,5}(v\d+)?", re.I)
MANUSCRIPT = re.compile(r"accepted manuscript|author manuscript|author'?s? (final|accepted)|hhs public access|"
                        r"\bpreprint\b|submitted to|manuscript received|nihms\d+", re.I)
PUBLISHER = re.compile(r"©|\(c\)\s*\d{4}|copyright|american physical society|iop publishing|john wiley|wiley "
                       r"periodicals|international press|elsevier|springer|american chemical society|"
                       r"published by", re.I)


def contract_doi(paper_id):
    text = run_b1.CONTRACT.read_text()
    sec = text[text.index("## Papers"):text.index("\n## ", text.index("## Papers") + 1)]
    for line in sec.splitlines():
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if cells and cells[0] == paper_id:
            return cells[1]
    raise SystemExit(f"{paper_id} not in the contract's Papers table")


def norm(s):
    return re.sub(r"\s+", "", s).lower()


def signals(text, doi):
    """Signals and kind from the first pages' text. Only booleans and matched stamps, never running text."""
    s = {"doi_found": norm(doi) in norm(text),
         "arxiv_stamp": (ARXIV.search(text).group(0) if ARXIV.search(text) else None),
         "manuscript_notice": bool(MANUSCRIPT.search(text)),
         "publisher_line": bool(PUBLISHER.search(text))}
    if s["arxiv_stamp"]:
        kind = "arxiv"
    elif s["manuscript_notice"]:
        kind = "accepted_manuscript"
    elif s["doi_found"] and s["publisher_line"]:
        kind = "publisher"
    else:
        kind = "unknown"
    return s, kind


def page_text(pdf, page):
    return subprocess.run(["mutool", "draw", "-q", "-F", "txt", "-o", "-", str(pdf), str(page)],
                          capture_output=True, text=True, check=True).stdout


def page_count(pdf):
    out = subprocess.run(["mutool", "show", str(pdf), "trailer/Root/Pages/Count"], capture_output=True,
                         text=True, check=True).stdout
    return int(out.strip().split()[-1])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--run", required=True)
    ap.add_argument("--paper-id", required=True)
    ap.add_argument("--pdf", required=True)
    a = ap.parse_args(argv)
    run, pdf = Path(a.run), Path(a.pdf)
    doi = contract_doi(a.paper_id)
    n = page_count(pdf)
    out = run / "untrusted" / a.paper_id / "copy"
    out.mkdir(parents=True, exist_ok=True)
    text = ""
    for p in range(1, min(n, 2) + 1):
        t = page_text(pdf, p)
        (out / f"p{p:03d}.txt").write_text(t)
        text += t + "\n"
    sig, kind = signals(text, doi)
    gate = kind != "publisher" or not sig["doi_found"]
    summary = {"paper_id": a.paper_id, "doi": doi, "pages": n, "kind": kind, **sig, "gate_opened": gate,
               "sha256": run_b1.sha256(pdf)}
    if gate:
        what = (f"Copy of {a.paper_id} is '{kind}'" + ("" if sig["doi_found"] else ", DOI not found on pages 1-2") +
                ": approve reading this copy (pages cited as PDF pages, copy named on every claim), or supply "
                "the version of record")
        open_gate(run, "A2", what, item=a.paper_id)
    (run / f"summary_copy_{a.paper_id}.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps({k: summary[k] for k in ("paper_id", "kind", "doi_found", "pages", "gate_opened")}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
