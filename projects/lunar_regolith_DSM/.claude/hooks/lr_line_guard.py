#!/usr/bin/env python3
"""H4, PreToolUse Bash: in the execute phase, block a direct call to a station script; the line runs stations."""
import sys

import lr_hook_lib as h

data = h.stdin_json()
cmd = (data.get("tool_input") or {}).get("command") or ""
try:
    phase, version = h.run_b1().phase()
    hits = h.calls_station(cmd) if phase == "execute" else []
except Exception:
    hits = []
if hits:
    sys.stderr.write(f"Blocked: {', '.join(hits)} is a station of the B1 line ({version}, execute phase). "
                     "Run it through run_b1.py run, which checks the stages it requires.\n")
    sys.exit(2)
