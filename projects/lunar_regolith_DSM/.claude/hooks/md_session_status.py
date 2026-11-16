#!/usr/bin/env python3
"""SessionStart and UserPromptSubmit: one line on where the model line stands (B1 is closed, 2026-11-16).

Prints `run_md.py status --brief` for the most recent model-line run, or a pointer when none exists yet.
Silent on any error, so a broken run folder never blocks a session.
"""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "runs" / "model-development"
try:
    runs = [p for p in RUNS.glob("*") if (p / "state.json").exists()] if RUNS.exists() else []
    if not runs:
        print("Model line: no run yet. Design: workflows/model-development/WORKFLOW.md; "
              "start: run_md.py init --run <id> --problem regolith")
        sys.exit(0)
    p = subprocess.run([str(ROOT / "venv_workflows" / "bin" / "python"),
                        str(ROOT / "scripts" / "model_development" / "run_md.py"), "status", "--brief"],
                       capture_output=True, text=True, timeout=8)
    if p.returncode == 0 and p.stdout.strip():
        print(p.stdout.strip().splitlines()[-1])
except Exception:
    pass
