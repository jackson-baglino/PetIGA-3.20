---
name: lr-status
description: Where the B1 literature-review line stands - phase (develop or execute), freeze version, current run, next stage, gates and review items waiting for the user. Use when asked about B1 progress, what is next in the literature review, or before running or changing anything in workflows/literature-review.
allowed-tools: Bash(venv_workflows/bin/python scripts/literature_review/run_b1.py status*), Bash(venv_workflows/bin/python scripts/literature_review/lr_status.py *)
---

## Current state

!`venv_workflows/bin/python scripts/literature_review/run_b1.py status 2>&1 || venv_workflows/bin/python scripts/literature_review/lr_status.py --brief 2>&1 || true`

## How to read it

- Each row is a stage of the stage graph (`run` blocks in `workflows/literature-review/flowchart.json`):
  whether its station script exists, items done out of total (papers, or 1), its gate, and gates waiting.
- **NOT BUILT** stations are development work; the line stops there.
- **waiting** gates need the user: name each one and give the `run_b1.py approve` command; never record a
  decision yourself.
- In the **execute** phase, frozen files are not edited; ideas go to `/lr-backlog`.

Report the state in two to four sentences: phase, run, next stage, and anything waiting for the user.
