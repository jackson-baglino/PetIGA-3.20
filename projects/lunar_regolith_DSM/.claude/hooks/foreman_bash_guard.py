#!/usr/bin/env python3
"""H6, beth-foreman PreToolUse Bash: the foreman may run the line and read its state, nothing else."""
import re
import sys

import lr_hook_lib as h

PY = r"(\./)?venv_workflows/bin/python scripts/literature_review/"
ALLOWED = [PY + r"run_b1\.py (init|status|run|verify)\b[^;&|<>`$]*",
           PY + r"lr_status\.py --(brief|write|check)\b[^;&|<>`$]*",
           PY + r"lr_backlog\.py (add|list)\b[^;&|<>`$]*",
           r"git (status|log)\b[^;&|<>`$]*"]

cmd = ((h.stdin_json().get("tool_input") or {}).get("command") or "").strip()
if not any(re.fullmatch(p, cmd) for p in ALLOWED):
    sys.stderr.write("beth-foreman may run only run_b1.py (init, status, run, verify), lr_status.py, "
                     "lr_backlog.py, git status and git log, one command at a time. Approvals, freezes and edits "
                     "are the user's.\n")
    sys.exit(2)
