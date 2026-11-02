#!/usr/bin/env python3
"""H2, UserPromptSubmit: one line of context while a B1 run is open; silent otherwise."""
import lr_hook_lib as h

try:
    run = h.run_b1().current_run()
    if run is not None:
        print(h.run_b1().brief(run))
except Exception:
    pass
