---
paths:
  - "workflows/**"
  - "scripts/literature_review/**"
  - "scripts/units_verification/**"
  - "scripts/design_critique/**"
  - "runs/**"
  - ".claude/permissions/**"
---

# Agentic workflows

- **The line.** Scripts manage; agents work (`workflows/TEMPLATE.md`, "The
  line"). Stages run in the order of the workflow's stage graph, through the
  line script (literature review: `run_b1.py`; the foreman session is
  `venv_workflows/bin/python scripts/literature_review/foreman.py`). In the develop phase a station may be run by
  hand only to test it, in a dev run folder; in the execute phase never (a
  hook blocks it).
- **Phases** (`workflows/TEMPLATE.md` §7). Before editing an agent file,
  prompt, schema, station script or model map, look for
  `workflows/<name>/freeze/`. No freeze file means develop. In execute, do not
  edit; park the idea in `workflows/<name>/BACKLOG.md`.
- **Finish the planned task first.** If you see an improvement while doing
  planned work, write it in `BACKLOG.md` with the date and the reason, finish
  the task, and raise the idea with the user at the end. Do not redesign in
  the middle of a task.
- **Where agents live.** Every agent of a workflow, workers and foreman, in
  `workflows/<name>/agents/<agent>.md`; `.claude/agents/` only for general
  helpers that belong to no workflow. All use the subagent format (frontmatter + a system-prompt body with a "Your place
  in the line" section). `prompts/` holds only per-call work orders.
- **Literature review.** Models are pinned in `MODELS` in
  `scripts/literature_review/agent_call.py`; permissions are in
  `.claude/permissions/literature-review/`. After any permission change run
  `check_permissions.py` (`--write` regenerates the tables in
  `PERMISSIONS.md`; never edit them by hand).
- **Untrusted text.** A normal session never reads `library/`, `index/` or
  `runs/literature-review/**/untrusted/`. Evidence is read only in review mode
  (`claude --settings
  .claude/permissions/literature-review/sessions/review-mode.json`).
