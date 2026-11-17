---
name: new-agent
description: Create a new agent in this project's house format - a script-run worker or foreman in workflows/<workflow>/agents/<name>.md with its own folder (settings.json, env.allow, schemas, tasks, tests) next to it, or a session helper in .claude/agents/ - with the naming rule, a pinned model ID, a "Your place in the line" section, and an enforcer for every "Never". Use when an agent is added or its model changes.
argument-hint: "<workflow> <purpose> <model family>"
---

Agents are designed with the user before they are built (`memory/USER.md`). Confirm the design first: role,
inputs, output schema, tools, who sends it work and who checks it. The layout and its reasons are in
`docs/planning/28_agent_anatomy.html` and CLAUDE.md, "Where agents live".

1. **Name:** `<character>-<purpose>` from the model family (CLAUDE.md, "Agent names"): Fable → `rick`, Opus →
   `beth`, Sonnet → `morty`, Haiku → `jerry`. A model change means a new name. Check that no agent with the same
   name exists anywhere in `workflows/*/agents/` or `.claude/agents/`.
2. **Where:**
   - A script-run agent (worker or foreman) goes with the workflow that owns it:
     ```
     workflows/<workflow>/agents/<name>.md     definition: frontmatter + system prompt
     workflows/<workflow>/agents/<name>/
         settings.json    permission allow/deny, hooks, sandbox (passed with --settings)
         env.allow        environment variables passed through, one per line
         schemas/         task.schema.json (work order), result.schema.json (answer)
         tasks/           one folder per task: task.json, and grade.py when there is a known answer
         tests/           the definition loads; every "Never" is enforced; red-team probes fail
     ```
     Another workflow that needs this agent refers to it by name; never copy it.
   - Only a helper that the main session starts as a subagent goes in `.claude/agents/<name>.md`.
   - Ask the user before writing a definition or a `settings.json` (CLAUDE.md).
3. **Frontmatter**, in this order: `name` (same as the file), `description` (what it does, the stage, who runs
   it), `tools` (the minimum; empty for none), `disallowedTools`, `model` (a full pinned ID such as
   `claude-opus-5-5`, never an alias), `permissionMode: dontAsk` for workers, `maxTurns`. No `memory` for
   workers (a fresh checker must be fresh) and no `effort` (workflow agents keep the default).
4. **Body** (the system prompt), following `workflows/literature-review/agents/beth-checker.md`:
   - one paragraph on the role;
   - `## Your place in the line`: who sends you work, what you return, who checks your work next, never;
   - `## When invoked`: numbered steps;
   - the rules of the job. Per-call data goes in a work order under `tasks/`, never in the body.
5. **Enforcers.** For each "Never" in the body, name what enforces it (tool list, permission rule, hook,
   sandbox, schema, grader) in the agent's `tests/`, or mark it "prompt only".
6. Add the tests (the definition loads; name, tools and model agree with `settings.json`; each enforcer
   blocks its probe), run them, then commit.
