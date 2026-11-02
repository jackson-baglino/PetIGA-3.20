---
name: new-agent
description: Create a new agent definition in this project's house format - a workflow's headless worker or foreman in workflows/<name>/agents/, or a general helper in .claude/agents/ - with the naming rule, a pinned model ID, a "Your place in the line" section, and (for literature-review workers) a permission file. Use when an agent is added or its model changes.
argument-hint: "<workflow> <purpose> <model family>"
---

Agents are designed with the user before they are built (`memory/USER.md`). Confirm the design first: role,
inputs, output schema, tools, who sends it work and who checks it.

1. **Name:** `<character>-<purpose>` from the model family (CLAUDE.md, "Agent names"): Fable → `rick`, Opus →
   `beth`, Sonnet → `morty`, Haiku → `jerry`. A model change means a new name.
2. **Where:** every agent of a workflow (headless workers and the foreman) goes in
   `workflows/<workflow>/agents/<name>.md`. A foreman is started by the workflow's launcher (for the literature
   review, `scripts/literature_review/foreman.py`, which passes the definition with `claude --agents`). Only a
   general helper that belongs to no workflow goes in `.claude/agents/<name>.md` (ask first: edits to `.claude/`
   need the user's approval).
3. **Frontmatter**, in this order: `name` (same as the file), `description` (what it does, the stage, who runs
   it), `tools` (the minimum; empty for none), `disallowedTools`, `model` (a full pinned ID such as
   `claude-opus-5-5`, never an alias), `permissionMode: dontAsk` for workers, `maxTurns`. No `memory` for
   workers (a fresh checker must be fresh) and no `effort` (workflow agents keep the default).
4. **Body** (the system prompt), following `workflows/literature-review/agents/beth-checker.md`:
   - one paragraph on the role;
   - `## Your place in the line`: who sends you work, what you return, who checks your work next, never;
   - `## When invoked`: numbered steps;
   - the rules of the job. Per-call data goes in a work order in `prompts/`, never in the body.
5. **Literature-review worker:** ask the user to approve a permission file
   `.claude/permissions/literature-review/agents/<name>.settings.json` and an `inputs.json` entry, then run
   `venv_workflows/bin/python scripts/literature_review/check_permissions.py --write`; it checks that the
   definition's tools, model and name agree with the permission file.
6. Add a test in the workflow's tests that the definition loads, then commit.
