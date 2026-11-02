#!/usr/bin/env python3
"""H3, PreToolUse Edit|Write: in the execute phase, ask before any edit to a frozen file."""
import json

import lr_hook_lib as h

data = h.stdin_json()
path = (data.get("tool_input") or {}).get("file_path") or ""
try:
    phase, version = h.run_b1().phase()
except Exception:
    phase, version = "develop", None
if phase == "execute" and path and h.is_frozen(path):
    print(json.dumps({"hookSpecificOutput": {
        "hookEventName": "PreToolUse", "permissionDecision": "ask",
        "permissionDecisionReason": f"{path} is frozen for {version}. Changing it ends the B1 run "
                                    "(thaw, fix, new freeze). To keep an idea for later, use /lr-backlog."}}))
