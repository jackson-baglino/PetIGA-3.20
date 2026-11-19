"""Bundle several committed files into one design file for a critic that may read only one file.

The critic (beth-critic) reads only WORKFLOW.md in its clean folder. To review a plan that lives in many files,
this script writes them, in the order of a list file, into one Markdown file: a preamble, then each file under
a header with its path, the git ref it came from, its line range and its sha256. The bundle is then passed to
run_critique.py with --design.

List file: one entry per line; blank lines and lines starting with # are skipped. An entry is
    path                       the whole file, from the working tree (must be committed and clean)
    path:START-END             lines START to END (1-based, inclusive)
    git:REF:path               the file as committed on REF (for a file on another branch)

    venv_workflows/bin/python scripts/design_critique/build_bundle.py \
        --list runs/design-critique/focus/<name>_files.txt --preamble <preamble.md> --out <bundle.md>
"""
import argparse
import hashlib
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def git(*args):
    return subprocess.run(["git", *args], cwd=ROOT, check=True, capture_output=True, text=True).stdout


def read_entry(entry):
    """(label, ref, text) for one list entry."""
    if entry.startswith("git:"):
        _, ref, path = entry.split(":", 2)
        return path, ref, git("show", f"{ref}:./{path}")
    m = re.fullmatch(r"(.+?):(\d+)-(\d+)", entry)
    path = m.group(1) if m else entry
    if git("status", "--porcelain", "--", path).strip():
        sys.exit(f"{path} has uncommitted changes; commit first")
    if not git("ls-files", "--", path).strip():
        sys.exit(f"{path} is not tracked by git")
    text = (ROOT / path).read_text()
    if m:
        lo, hi = int(m.group(2)), int(m.group(3))
        lines = text.splitlines(keepends=True)
        if hi > len(lines):
            sys.exit(f"{entry}: the file has only {len(lines)} lines")
        text = "".join(lines[lo - 1:hi])
        path = f"{path}, lines {lo}-{hi}"
    return path, "HEAD", text


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--list", required=True)
    ap.add_argument("--preamble", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    entries = [ln.strip() for ln in (ROOT / args.list).read_text().splitlines()
               if ln.strip() and not ln.strip().startswith("#")]
    head = git("rev-parse", "HEAD").strip()
    parts = [(ROOT / args.preamble).read_text().rstrip() + "\n",
             f"\nBuilt from commit {head[:12]}. Parts, in order:\n"]
    bodies = []
    for i, entry in enumerate(entries, 1):
        label, ref, text = read_entry(entry)
        sha = hashlib.sha256(text.encode()).hexdigest()
        parts.append(f"{i}. `{label}`\n")
        fence = "~~~~" if "```" in text else "```"
        bodies.append(f"\n\n# PART {i}: `{label}`\n\n(ref {ref}, sha256 {sha[:16]}, {len(text.splitlines())} lines)\n\n"
                      f"{fence}\n{text.rstrip()}\n{fence}\n")
        print(f"{i:2d}. {label}  ({len(text.splitlines())} lines, sha256 {sha[:12]})")
    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("".join(parts) + "".join(bodies))
    text = out.read_text()
    print(f"wrote {args.out}: {len(text.splitlines())} lines, {len(text)} characters")


if __name__ == "__main__":
    main()
