#!/usr/bin/env python3
"""Start the B1 review session: Claude may read one run's evidence with the user, and nothing else (WORKFLOW.md §4).

    venv_workflows/bin/python scripts/literature_review/review_session.py [--run <run folder>] [extra claude args]

Runs `claude --setting-sources "" --settings .claude/permissions/literature-review/sessions/review-session.json
--strict-mcp-config --tools Read,Grep,Glob --model claude-opus-5-5` from the project root with LR_RUN set to the run (default: the newest
B1 run). The project settings are not loaded, because their deny rules keep evidence out of normal sessions; the
review-session file instead denies Bash, web tools, edits and agents, and its hook allows Read, Grep and Glob only
on that run's evidence (review_session_allowlist.py). The user records verdicts with u3_review.py in their own
terminal, not in this session.
"""
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_b1  # noqa: E402

SETTINGS = run_b1.PROJECT / ".claude/permissions/literature-review/sessions/review-session.json"
MODEL = "claude-opus-5-5"
TOOLS = ("Read", "Grep", "Glob")   # --tools: no other tool exists in the session (checked live 2026-11-01:
#                                     without it, PushNotification, RemoteTrigger, SendMessage and others were present)


def argv_for(run, extra):
    return ["claude", "--setting-sources", "", "--settings", str(SETTINGS), "--strict-mcp-config",
            "--tools", ",".join(TOOLS), "--model", MODEL] + extra


def main(argv):
    run = None
    if "--run" in argv:
        i = argv.index("--run")
        run, argv = argv[i + 1], argv[:i] + argv[i + 2:]
    run = run or run_b1.current_run()
    if run is None:
        raise SystemExit("no B1 run: give --run")
    os.environ["LR_RUN"] = str(Path(run).resolve())
    os.environ["CLAUDE_PROJECT_DIR"] = str(run_b1.PROJECT)
    os.chdir(run_b1.PROJECT)
    os.execvp("claude", argv_for(run, argv))


if __name__ == "__main__":
    main(sys.argv[1:])
