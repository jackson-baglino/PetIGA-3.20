#!/usr/bin/env python3
"""Model-development line, Stop hook for headless agents (WORKFLOW.md §3.2).

Runs the stage's structural pre-check (MD_STOP_CHECK, a shell command set by md_agent.py for the call). If it fails,
the stop is blocked and the check's output is fed back, so the agent finishes the work instead of ending its turn.
At most two blocks per call (counted in MD_HOOK_STATE): after that the agent may stop, and the station's own
post-check reports the failure. An agent that cannot meet the criterion returns status "blocked" with the reason;
the third stop is then allowed. Does nothing in a normal session (no MD_AGENT).
"""
import json
import os
import subprocess
import sys
from pathlib import Path

MAX_BLOCKS = 2


def main():
    if not os.environ.get("MD_AGENT") or not os.environ.get("MD_STOP_CHECK"):
        return 0
    try:
        json.load(sys.stdin)
    except Exception:
        pass
    state = Path(os.environ.get("MD_HOOK_STATE", "/tmp")) / "stop_blocks"
    n = int(state.read_text()) if state.exists() else 0
    if n >= MAX_BLOCKS:
        return 0
    p = subprocess.run(os.environ["MD_STOP_CHECK"], shell=True, capture_output=True, text=True,
                       cwd=os.environ.get("MD_CWD") or None, timeout=300)
    if p.returncode == 0:
        return 0
    state.write_text(str(n + 1))
    sys.stderr.write(
        "The stage's done check fails, so the work is not finished:\n"
        + (p.stdout + p.stderr).strip()[-1500:]
        + "\nFix what it reports and check again. If you cannot meet the criterion, return status \"blocked\" "
          "with the reason and the questions for the user.\n")
    return 2


if __name__ == "__main__":
    sys.exit(main())
