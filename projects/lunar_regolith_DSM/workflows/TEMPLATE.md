# WORKFLOW: <name>

Copy this file to `workflows/<name>/WORKFLOW.md` when a workflow is designed.
The flowchart is `workflows/<name>/flowchart.html`. Nothing is built until the
user approves both.

## 1. Purpose

- Question the workflow answers, and what a finished output looks like.
- Acceptance test: how we will know the output is correct.
- Baseline: what a non-agent approach would produce, for comparison.

## 2. Stages

One row per stage. Each stage ends in a file that something can check.

| # | Stage | Agent or script | Input (path, schema) | Output (path, schema) | Check | Stop / escalate |
|---|---|---|---|---|---|---|

## 3. Agents

One block per agent. The definition lives in `.claude/agents/<name>.md`.

```
name:         <kebab-case>
model:        <full model ID, pinned>
tools:        <minimal list>
skills:       <list or none>
MCP servers:  <list or none>
reads:        <exact inputs>
writes:       <exact output path and schema>
may not:      <explicit limits>
context:      fresh per task | shared (why)
```

## 4. Checks

| Output | Check | Kind (deterministic / independent agent / vote / user) | Pass criterion |
|---|---|---|---|

Measure each agent checker before trusting it: plant known errors in a sample
and record the fraction it catches.

## 5. Human decision points

Where the user approves, chooses or inspects, and what they are shown.

## 6. Reproducibility

Each run writes `runs/<name>/<run_id>/`:

- `manifest.json`: run ID, date, git commit, Claude Code version, model ID per
  agent, prompt-file hashes, input-file hashes, check verdicts.
- `transcript.jsonl`: from `claude -p ... --output-format stream-json`.
- `runs.db` (SQLite): one row per agent call (stage, agent, model, input
  hash, output path, verdict, timestamp), so a run can be resumed and audited.

Re-running a stage with the same inputs must reach the same verified
conclusions, though not identical wording.

## 7. Risks and failure modes

What can go wrong silently, and which check catches it.
