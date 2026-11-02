#!/usr/bin/env python3
"""Park an improvement idea in workflows/literature-review/BACKLOG.md instead of acting on it (.claude/rules/workflows.md).

    lr_backlog.py add --source <session, stage or run> --idea <text> --why <text>
    lr_backlog.py list

Ideas are read and decided at the start of the next develop phase. This only appends under "## Open".
"""
import argparse
import datetime
import sys
from pathlib import Path

BACKLOG = Path(__file__).resolve().parents[2] / "workflows" / "literature-review" / "BACKLOG.md"


def add(source, idea, why, path=BACKLOG):
    text = path.read_text()
    if "## Open\n" not in text or "\n## Decided" not in text:
        raise SystemExit(f"{path} lacks '## Open' or '## Decided'")
    entry = f"- {datetime.date.today().isoformat()} · {source} · {' '.join(idea.split())} " \
            f"Why: {' '.join(why.split())}\n"
    i = text.index("\n## Decided")
    head = text[:i].rstrip("\n") + "\n"
    path.write_text(head + entry + "\n" + text[i + 1:])
    return entry


def open_items(path=BACKLOG):
    text = path.read_text()
    sec = text[text.index("## Open\n") + 8:text.index("\n## Decided")]
    return [line for line in sec.splitlines() if line.startswith("- ")]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("add")
    p.add_argument("--source", required=True)
    p.add_argument("--idea", required=True)
    p.add_argument("--why", required=True)
    sub.add_parser("list")
    a = ap.parse_args(argv)
    if a.cmd == "add":
        print("parked:", add(a.source, a.idea, a.why).strip())
    else:
        items = open_items()
        print(f"{len(items)} open idea(s)")
        for line in items:
            print(line)
    return 0


if __name__ == "__main__":
    sys.exit(main())
