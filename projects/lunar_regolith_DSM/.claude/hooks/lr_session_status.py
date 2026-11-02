#!/usr/bin/env python3
"""H1, SessionStart: print where the B1 line stands (phase, freeze, run, next stage, waiting gates)."""
import lr_hook_lib as h

try:
    run = h.run_b1().current_run()
    print("Literature review: " + h.run_b1().brief(run) +
          ". Rules: .claude/rules/workflows.md; details: run_b1.py status.")
except Exception as e:                       # a status hook must never block a session
    print(f"Literature review status unavailable: {type(e).__name__}")
