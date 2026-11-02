#!/usr/bin/env python3
"""T-det (WORKFLOW.md §12): is P1 deterministic? Compare two P1 runs of one paper, file by file.

    t_det.py --a <run folder> --b <run folder> --paper-id <id> --pdf <file>

Both runs must hold P1 content for the same copy (<run>/untrusted/content/<paper>/<sha256[:16]>/). Prints, for
each file in either folder, whether its sha256 is the same in both: names and verdicts only, never content (the
content is untrusted text). The run manifest written by P1 records times and paths, so manifest.json is compared
without its timing fields. Exit 0 if every file matches, 1 otherwise.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

VOLATILE = {"made", "seconds_per_page", "seconds", "created", "date", "time", "run", "started", "finished", "elapsed"}


def strip(obj):
    if isinstance(obj, dict):
        return {k: strip(v) for k, v in obj.items() if k not in VOLATILE}
    if isinstance(obj, list):
        return [strip(v) for v in obj]
    return obj


def digest(path):
    if path.name == "manifest.json":
        return hashlib.sha256(json.dumps(strip(json.loads(path.read_text())), sort_keys=True).encode()).hexdigest()
    return hashlib.sha256(path.read_bytes()).hexdigest()


def tree(root):
    return {str(p.relative_to(root)): digest(p) for p in sorted(root.rglob("*")) if p.is_file()}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--a", required=True)
    ap.add_argument("--b", required=True)
    ap.add_argument("--paper-id", required=True)
    ap.add_argument("--pdf", required=True)
    x = ap.parse_args(argv)
    h = hashlib.sha256(Path(x.pdf).read_bytes()).hexdigest()[:16]
    ta, tb = (tree(Path(r) / "untrusted" / "content" / x.paper_id / h) for r in (x.a, x.b))
    if not ta or not tb:
        raise SystemExit(f"no P1 content for {x.paper_id} ({h}) in {'--a' if not ta else '--b'}")
    names = sorted(set(ta) | set(tb))
    diff = [n for n in names if ta.get(n) != tb.get(n)]
    by_kind = {}
    for n in diff:
        k = n.split("/")[0] if "/" in n else Path(n).suffix or n
        by_kind[k] = by_kind.get(k, 0) + 1
    keys = []
    if "manifest.json" in diff:                  # name the differing manifest fields (P1's own metadata)
        ma, mb = (json.loads((Path(r) / "untrusted" / "content" / x.paper_id / h / "manifest.json").read_text())
                  for r in (x.a, x.b))

        def walk(u, v, path=""):
            if isinstance(u, dict) and isinstance(v, dict):
                for k in sorted(set(u) | set(v)):
                    walk(u.get(k), v.get(k), f"{path}.{k}" if path else k)
            elif u != v:
                keys.append(path)
        walk(strip(ma), strip(mb))
    print(json.dumps({"paper": x.paper_id, "files": len(names), "identical": len(names) - len(diff),
                      "manifest_fields_different": keys[:40],
                      "different": len(diff), "different_by_kind": by_kind,
                      "different_files": diff[:40]}, indent=1))
    return 0 if not diff else 1


if __name__ == "__main__":
    sys.exit(main())
