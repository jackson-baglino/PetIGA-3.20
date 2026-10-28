#!/usr/bin/env python3
"""PreToolUse hook for review mode: Read, Grep and Glob only inside an allow-list.

Permission rules cannot express an allow-list (a `!` carve-out cannot reach an
anchored rule or reopen a blocked directory; Claude Code docs, permissions,
read 2026-10-28), so this hook does it. It reads the tool call as JSON on
stdin and exits 2, which blocks the call, unless every path the call names is
inside one of ALLOWED. A call that names no path (Grep or Glob over the
working directory) is blocked too.
"""
import json
import os
import sys

PROJECT = os.environ.get("CLAUDE_PROJECT_DIR") or os.getcwd()

# Paths relative to the project root that review mode may read.
ALLOWED = [
    "runs/literature-review",          # drivers' summaries (schema-checked, no free text)
    "scripts",
    "workflows/literature-review/WORKFLOW.md",
    "workflows/literature-review/schemas",
    ".claude/permissions",
    ".claude/settings.json",
    "library/reviews",                 # briefs
    "CLAUDE.md",
    "CLAUDE.local.md",
]
# Inside runs/literature-review only these file names may be read.
RUNS_FILES = ("summary.json", "manifest.json")


def allowed(path):
    p = os.path.realpath(os.path.join(PROJECT, path))
    root = os.path.realpath(PROJECT)
    if not p.startswith(root + os.sep):
        return False
    rel = os.path.relpath(p, root)
    for a in ALLOWED:
        if rel == a or rel.startswith(a + os.sep):
            if a == "runs/literature-review":
                return os.path.basename(rel) in RUNS_FILES
            if a == "library/reviews":
                name = os.path.basename(rel)
                return name.startswith("brief.") and name.endswith(".json")
            return True
    return False


def main():
    call = json.load(sys.stdin)
    tool = call.get("tool_name", "")
    args = call.get("tool_input", {}) or {}
    if tool == "Read":
        paths = [args.get("file_path")]
    elif tool in ("Grep", "Glob"):
        paths = [args.get("path")]
    else:
        return 0
    for p in paths:
        if not p or not allowed(p):
            print(f"review mode: {tool} of {p or 'the working directory'} is outside the "
                  "allow-list (.claude/permissions/literature-review/hooks/review_read_allowlist.py)",
                  file=sys.stderr)
            return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
