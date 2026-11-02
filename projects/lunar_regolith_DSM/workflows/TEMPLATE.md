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

One block per agent. A worker that a script calls headless is defined in
`workflows/<name>/agents/<agent>.md`; an agent the user starts or talks to in a
session (a foreman, an interactive helper) is defined in
`.claude/agents/<agent>.md`. Both use the Claude Code subagent format: YAML
frontmatter (`name`, `description`, `tools`, `disallowedTools`, `model`,
`permissionMode`, `maxTurns`) and the system prompt as the body. The body has a
short "Your place in the line" section: who sends the work order, what
it returns, who checks it next, and what it must never do.

```
name:         <rick|beth|morty|jerry>-<purpose>   (Fable|Opus|Sonnet|Haiku; must match model)
model:        <full model ID, pinned>
tools:        <minimal list>
skills:       <list or none>
MCP servers:  <list or none>
reads:        <exact inputs>
writes:       <exact output path and schema>
may not:      <explicit limits>
context:      fresh per task | shared (why)
```

### The line

Who reports to whom. Scripts manage; agents work.

```
user (owner)            approves gates, the freeze and the final output
  foreman agent         a session started with `claude --agent <foreman>`; runs the line,
                        reports, stops at blocking gates; never edits definitions
    line script         walks the stage graph in `flowchart.json`; refuses a stage whose
                        upstream stages are not done; records every call in runs.db
      station scripts   one per stage: write each work order from a fixed template,
                        validate the output against its schema, audit the transcript,
                        retry once, write report.json (counts, failures, gate items)
        workers         headless agents, one task per call, JSON out, no writes
```

List each gate as *blocking* (the line stops and the user records a decision)
or *review* (queued; the line goes on with other items).

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

Rerunning a stage with the same inputs must reach the same verified
conclusions, though not identical wording.

## 7. Phases and freeze

- **Develop.** Definitions (agents, prompts, schemas, station scripts, model
  IDs, tool versions) may change. Runs are tagged `dev` and are not results.
- **Freeze.** The user approves a freeze: a file
  `workflows/<name>/freeze/<run>-vN.json` with the hash of every definition,
  and a git tag. The phase becomes *execute*.
- **Execute.** The line script refuses to start if a frozen hash differs, and
  a hook asks before any edit to a frozen file. Ideas for improvement go to
  `workflows/<name>/BACKLOG.md` and wait for the next develop phase.
- **The one way back.** A real error found by a check stops the run. The
  phase returns to develop, the fix is made and recorded, and a new freeze
  (vN+1) is approved before the run resumes.

## 8. Risks and failure modes

What can go wrong silently, and which check catches it.
