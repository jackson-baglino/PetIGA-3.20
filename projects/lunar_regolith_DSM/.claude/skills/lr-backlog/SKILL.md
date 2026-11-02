---
name: lr-backlog
description: Park an idea for improving the literature-review workflow in BACKLOG.md instead of acting on it. Use whenever, during planned work or a run, you notice something that could be changed in an agent, prompt, schema, station or the design - and when the user says "later", "park it" or "add to the backlog".
argument-hint: "<the idea and why>"
allowed-tools: Bash(venv_workflows/bin/python scripts/literature_review/lr_backlog.py *)
---

Ideas are decided at the start of the next develop phase, not when they are noticed (`.claude/rules/workflows.md`).

1. Write the idea in one or two sentences, and why it would help, in plain words. The idea: $ARGUMENTS
2. Name the source: the session, stage or run where it came up.
3. Run `venv_workflows/bin/python scripts/literature_review/lr_backlog.py add --source "<source>" --idea "<idea>" --why "<why>"`.
4. Go back to the planned task. Mention the parked idea to the user at the end.
