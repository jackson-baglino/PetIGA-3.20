#!/usr/bin/env python3
"""agentkit provenance: the PROVENANCE.md of a workflow run, written only by the launcher (doc 29 §10).

Each entry is a Markdown section that ends with a chain line,

    `chain: <sha256 of (previous chain hash + the entry's text)>`

so an edit anywhere breaks every later hash, and verify() names the first entry that no longer matches. The
final hash is also stored in manifest.json, so the whole file can be checked against the run record.

    python scripts/agentkit/provenance.py verify <run>/PROVENANCE.md
"""
import hashlib
import re
import sys
from pathlib import Path

GENESIS = "0" * 64
CHAIN = re.compile(r"^`chain: ([0-9a-f]{64})`$")


def _hash(prev, text):
    return hashlib.sha256((prev + "\n" + text).encode()).hexdigest()


class Provenance:
    def __init__(self, path, title, preamble):
        self.path = Path(path)
        self.n = 0
        self.last = GENESIS
        self.path.write_text(f"# {title}\n\n{preamble.strip()}\n")

    def add(self, heading, lines):
        """Append one entry: a heading and a list of lines (strings), then its chain line. Returns the hash."""
        self.n += 1
        text = f"## {self.n}. {heading}\n\n" + "\n".join(f"- {s}" for s in lines) + "\n"
        self.last = _hash(self.last, text)
        with self.path.open("a") as f:
            f.write("\n" + text + f"\n`chain: {self.last}`\n")
        return self.last


def verify(path):
    """(ok, message, final hash). Recomputes the chain entry by entry."""
    body = Path(path).read_text()
    parts = re.split(r"(?m)^(?=## \d+\. )", body)[1:]
    if not parts:
        return False, "no entries", None
    prev = GENESIS
    for i, part in enumerate(parts, 1):
        lines = part.rstrip("\n").split("\n")
        m = CHAIN.match(lines[-1]) if lines else None
        if not m:
            return False, f"entry {i} has no chain line", None
        text = "\n".join(lines[:-1]).rstrip("\n") + "\n"
        if _hash(prev, text) != m.group(1):
            return False, f"entry {i} was changed after it was written", None
        prev = m.group(1)
    return True, f"{len(parts)} entries, chain intact", prev


if __name__ == "__main__":
    if len(sys.argv) != 3 or sys.argv[1] != "verify":
        sys.exit("usage: provenance.py verify <PROVENANCE.md>")
    ok, msg, last = verify(sys.argv[2])
    print(f"{'OK' if ok else 'BROKEN'}: {msg}" + (f" (final {last[:12]})" if last else ""))
    sys.exit(0 if ok else 1)
