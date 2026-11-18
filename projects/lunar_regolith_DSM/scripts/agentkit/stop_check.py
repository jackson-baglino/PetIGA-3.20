#!/usr/bin/env python3
"""agentkit stop check: the Stop hook for every script-run agent (doc 28, step 2).

Claude Code runs this when the agent tries to end its turn. If a file the work order requires is missing or
empty in the workspace, the stop is blocked (exit 2) and the reason goes back to the agent, so it finishes the
work instead of only saying it did. The launcher sets:

    AGENTKIT_WORKSPACE    the agent's working folder
    AGENTKIT_REQUIRED     JSON list of file names that must exist there, non-empty, before the agent stops

Blocks at most once: when Claude Code reports the stop hook already fired (stop_hook_active), the stop is
allowed and the launcher's own check records the missing file. On an error it allows the stop for the same
reason; the launcher, not this hook, decides the verdict.
"""
import json
import os
import sys
from pathlib import Path


def missing(workspace, required):
    return [n for n in required if not (Path(workspace) / n).is_file() or (Path(workspace) / n).stat().st_size == 0]


def main():
    try:
        event = json.load(sys.stdin)
        if event.get("stop_hook_active"):
            return 0
        gone = missing(os.environ["AGENTKIT_WORKSPACE"], json.loads(os.environ.get("AGENTKIT_REQUIRED", "[]")))
    except Exception as e:
        sys.stderr.write(f"agentkit stop check error ({e}); stop allowed, the launcher will check\n")
        return 0
    if not gone:
        return 0
    sys.stderr.write(f"Not finished: {', '.join(gone)} is missing or empty in your working folder. Write it, "
                     "then return your answer. If you cannot, return status \"blocked\" with the reason.\n")
    return 2


if __name__ == "__main__":
    sys.exit(main())
