#!/usr/bin/env python3
"""H5, Stop: if a B1 run's state is newer than the contract's generated line status, ask for a refresh once."""
import json

import lr_hook_lib as h

data = h.stdin_json()
if not data.get("stop_hook_active"):
    try:
        import lr_status
        rb = h.run_b1()
        run = rb.current_run()
        if run is not None and lr_status.main(["--check", "--run", run]) != 0:
            print(json.dumps({"decision": "block",
                              "reason": "The B1 line status in B1_CONTRACT.md is stale. Run "
                                        "venv_workflows/bin/python scripts/literature_review/lr_status.py --write "
                                        "and commit it with the work it records."}))
    except Exception:
        pass
