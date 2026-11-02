#!/usr/bin/env python3
"""PreToolUse hook for the B1 review session: Read, Grep and Glob only on that run's evidence (WORKFLOW.md §4).

The session is started by scripts/literature_review/review_session.py, which sets LR_RUN to the run folder. The
session may read only:
- in the run: summaries, gates, the review queue, u3 decisions, and under untrusted/ each paper's read, quote
  check, checklist, verdicts, final claims and page images, the memos and the X1 comparison (not transcripts,
  logs or P1's parser output);
- in library/papers/<paper>/ for the contract's papers: cards, claims and analysis records;
- WORKFLOW.md and B1_CONTRACT.md.
Anything else, and a Grep or Glob with no path, exits 2 (blocked). The session has no Bash, no web tools and no
writes (review-session.json), so text planted in a paper can at worst mislead the discussion.
"""
import json
import os
import re
import sys

PROJECT = os.path.realpath(os.environ.get("CLAUDE_PROJECT_DIR") or os.getcwd())
RUN = os.path.realpath(os.environ.get("LR_RUN", "")) if os.environ.get("LR_RUN") else None
WF = "workflows/literature-review"
RUN_TOP = re.compile(r"(summary_[\w.]+\.json|manifest\.json|review_queue\.json|u3_decisions\.json|"
                     r"u3_missed\.json|gates/[\w.]+\.json)$")
RUN_UNTRUSTED = re.compile(r"untrusted/(memos/[\w.]+|x1/compare\.json|[a-z0-9_]+/(read|e3|checklist|verdicts|"
                           r"blind|blind_compare|final|revisions)\.json|[a-z0-9_]+/pages/p\d{3}\.png)$")
LIBRARY = re.compile(r"library/papers/([a-z0-9_]+)/(card\.v\d+\.json|claims/[\w.-]+\.json|analysis/.+\.json)$")


def contract_papers():
    try:
        text = open(os.path.join(PROJECT, WF, "B1_CONTRACT.md")).read()
    except OSError:
        return set()
    sec = text[text.index("## Papers"):]
    return set(re.findall(r"^\|\s*([a-z0-9_]+)\s*\|", sec, flags=re.M)) - {"id"}


def allowed(path):
    p = os.path.realpath(os.path.join(PROJECT, path))
    if RUN and (p == RUN or p.startswith(RUN + os.sep)):
        rel = os.path.relpath(p, RUN)
        return bool(RUN_TOP.fullmatch(rel) or RUN_UNTRUSTED.fullmatch(rel))
    if not p.startswith(PROJECT + os.sep):
        return False
    rel = os.path.relpath(p, PROJECT)
    if rel in (f"{WF}/WORKFLOW.md", f"{WF}/B1_CONTRACT.md"):
        return True
    m = LIBRARY.fullmatch(rel)
    return bool(m and m.group(1) in contract_papers())


def main():
    call = json.load(sys.stdin)
    tool = call.get("tool_name", "")
    args = call.get("tool_input", {}) or {}
    if tool not in ("Read", "Grep", "Glob"):
        return 0
    path = args.get("file_path") if tool == "Read" else args.get("path")
    if RUN is None:
        print("review session: LR_RUN is not set; start it with review_session.py", file=sys.stderr)
        return 2
    if not path or not allowed(path):
        print(f"review session: {tool} of {path or 'the working directory'} is outside this run's evidence "
              "(.claude/permissions/literature-review/hooks/review_session_allowlist.py)", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
