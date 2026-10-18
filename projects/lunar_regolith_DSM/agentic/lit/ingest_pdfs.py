#!/usr/bin/env python3
"""Drop-in step for PDFs you download yourself (through Caltech or Cornell access, Zotero, a library request, ...).
Scans a folder for recent PDFs, finds each paper's DOI in its first two pages, copies the file into the incoming folder as
<doi>.pdf and links it to the registry row. Nothing is fetched from the internet, no credentials are used or stored.

  python -m agentic.lit.ingest_pdfs --from ~/Downloads --hours 24 [--watch] [--out research/Literature/Wave1_incoming]
  --watch keeps polling every 15 s, so you can download in your browser and the files file themselves.
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import sqlite3
import subprocess
import time
from pathlib import Path

from agentic.lit import common

DOI = re.compile(r"\b(10\.\d{4,9}/[^\s\"<>,;()\[\]]+)", re.I)


def doi_of(pdf: Path) -> str | None:
    txt = ""
    for pg in (1, 2):
        txt += subprocess.run(["mutool", "draw", "-q", "-F", "txt", "-o", "-", str(pdf), str(pg)], capture_output=True, text=True).stdout
    # prefer a DOI that sits on the page header or footer line with 'doi'; fall back to the first match
    for m in re.finditer(r"doi[:\s/]*\s*(10\.\d{4,9}/[^\s\"<>,;()\[\]]+)", txt, re.I):
        return m.group(1).rstrip(".").lower()
    m = DOI.search(txt)
    return m.group(1).rstrip(".").lower() if m else None


def run_once(src: Path, out: Path, hours: float, con, seen: set) -> int:
    n = 0
    cutoff = time.time() - hours * 3600
    for f in sorted(src.glob("*.pdf")):
        if f in seen or f.stat().st_mtime < cutoff:
            continue
        seen.add(f)
        doi = doi_of(f)
        if not doi:
            print(f"  no DOI found in {f.name}"); continue
        row = con.execute("SELECT id, pdf_paths FROM papers WHERE doi=?", (doi,)).fetchone()
        if not row:
            print(f"  {f.name}: DOI {doi} is not in the registry (kept untouched)"); continue
        dest = out / (re.sub(r"[^A-Za-z0-9._-]+", "_", doi) + ".pdf")
        if not dest.exists():
            shutil.copy2(f, dest)
        lst = json.loads(row[1] or "[]")
        if str(dest) not in lst:
            con.execute("UPDATE papers SET pdf_paths=? WHERE id=?", (json.dumps(lst + [str(dest)]), row[0]))
            con.commit()
        print(f"  filed {f.name} -> {dest.name}  ({doi})")
        n += 1
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--from", dest="src", default=str(Path.home() / "Downloads"))
    ap.add_argument("--out", default=str(common.ROOT / "research" / "Literature" / "Wave1_incoming"))
    ap.add_argument("--hours", type=float, default=24)
    ap.add_argument("--watch", action="store_true")
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(common.STATE / "registry.sqlite")
    seen: set = set()
    total = run_once(Path(a.src), out, a.hours, con, seen)
    while a.watch:
        time.sleep(15)
        total += run_once(Path(a.src), out, a.hours, con, seen)
    print(f"filed {total} PDFs")


if __name__ == "__main__":
    main()
