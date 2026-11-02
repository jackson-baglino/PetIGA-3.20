---
name: lr-freeze
description: Close the development phase of the B1 literature review and freeze its definitions (agents, prompts, schemas, stations, stage graph) for an execute run. The user starts this; Claude prepares and checks, the user runs the freeze command.
disable-model-invocation: true
allowed-tools: Bash(venv_workflows/bin/python -m unittest *), Bash(venv_workflows/bin/python scripts/literature_review/check_permissions.py), Bash(venv_workflows/bin/python scripts/literature_review/run_b1.py status*), Bash(git status*), Bash(git log*)
---

A freeze fixes every definition the B1 line uses (`run_b1.py`, `FROZEN`) and starts the execute phase. Only the
user freezes. Prepare it in this order and stop at the first failure:

1. `git status --short -- workflows/literature-review scripts`: must be empty (a freeze matches a commit).
2. Tests: `venv_workflows/bin/python -m unittest discover -s scripts/literature_review/tests`.
3. Permissions: `venv_workflows/bin/python scripts/literature_review/check_permissions.py`.
4. `run_b1.py status`: every stage must be built. List any that is not.
5. The calibration diagnostic on the tuning set has been rerun since the last change to an agent file, prompt,
   schema or model (WORKFLOW.md §5, "Reruns"). Say which run and its result; if you cannot tell, say so.
6. Open BACKLOG items: list them; the user decides whether any must be done before the freeze.

Then show the user what the freeze will hold and give the command for them to run:

`! venv_workflows/bin/python scripts/literature_review/run_b1.py freeze --approved-by-user`

After they run it: commit `workflows/literature-review/freeze/b1-vN.json`, tag it `lr-b1-vN`
(`git tag lr-b1-vN`), and tell them to start the run with `claude --agent beth-foreman`.
