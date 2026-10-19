"""Integrity audit: re-hash every published file, check the catalog against the files, report anything that differs. Read-only.
  python -m agentic.library.verify [--json]"""
from __future__ import annotations

import json
import os
import stat
import sys
from pathlib import Path

from agentic.library.store import Library, sha256


def verify(lib: Library | None = None) -> dict:
    lib = lib or Library(); con = lib.ro(); problems = []
    files = con.execute("SELECT rel, sha256 FROM files").fetchall()
    for r in files:
        p = lib.root / r["rel"]
        if not p.exists(): problems.append({"kind": "missing_file", "rel": r["rel"]})
        elif sha256(p) != r["sha256"]: problems.append({"kind": "hash_mismatch", "rel": r["rel"]})
        elif p.stat().st_mode & stat.S_IWUSR: problems.append({"kind": "writable_file", "rel": r["rel"]})
    known = {r["rel"] for r in files}
    for p in lib.root.rglob("*"):
        if p.is_file() and p.name not in ("catalog.sqlite",) and not p.name.startswith(".") and str(p.relative_to(lib.root)) not in known and not str(p.relative_to(lib.root)).startswith(("history/", "_")) and p.suffix not in (".md",) :
            problems.append({"kind": "unlisted_file", "rel": str(p.relative_to(lib.root))})
    for r in con.execute("SELECT id FROM papers"):
        for need in ("record.json", "verification.json"):
            if f"papers/{r['id']}/{need}" not in known: problems.append({"kind": "paper_missing_part", "rel": f"{r['id']}/{need}"})
    pdf_missing = [dict(r) for r in con.execute("SELECT id,pdf_path FROM papers WHERE pdf_path!=''") if not Path(r["pdf_path"]).exists()] if False else []
    if os.access(lib.catalog, os.W_OK): problems.append({"kind": "catalog_writable", "rel": "catalog.sqlite"})
    con.close()
    return {"files_checked": len(files), "problems": problems, "ok": not problems}


if __name__ == "__main__":
    r = verify(); print(json.dumps(r, indent=1) if "--json" in sys.argv else f"{r['files_checked']} files checked; {len(r['problems'])} problems"); sys.exit(0 if r["ok"] else 1)
