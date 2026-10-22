# Instruction files, memory and checkpoints

## Instruction files (read by Claude Code at session start)

| File | Scope | In git | Purpose |
|---|---|---|---|
| `~/.claude/CLAUDE.md` | you, every project on this machine | no | personal habits that apply everywhere |
| `CLAUDE.md` (project root or `.claude/`) | everyone working in the project | yes | how to work in this repo |
| `CLAUDE.local.md` | you, this project | no | personal or private project rules (course rules here) |
| `CLAUDE.md` in a subfolder | that subfolder | yes | loaded only when Claude touches files there |

All found files are concatenated; nothing overrides. `@path` lines import
other files (up to four levels deep). The text is advisory context: Claude
reads it, but nothing enforces it. Rules that must hold go in
`settings.json` permissions or hooks.

**Put in CLAUDE.md:** commands Claude cannot guess (build, run, test), style
rules that differ from defaults, repo etiquette, traps, hard limits, and
pointers to where things live.

**Leave out:** anything readable from the code, tutorials, history, facts
that change weekly, and long explanations (link to a doc instead). Keep each
file under about 200 lines; every line costs context in every session.

## The project's other Markdown files

| File | Kind of memory | Written by | Loaded |
|---|---|---|---|
| `memory/MEMORY.md` | curated durable facts | agent proposes, user accepts | every session (imported) |
| `memory/USER.md` | user preferences as dated directives | same | every session (imported) |
| `memory/semantic/*.md` | semantic: units, constants | user or verified workflow | on demand |
| `memory/README.md` | the memory decision table | user | on demand |
| `PERMISSIONS.md` | spec of allowed actions | user | on demand |
| `.claude/agents/<name>.md` | procedural: one agent | designed together | when the agent runs |
| `.claude/skills/<name>/SKILL.md` | procedural: a reusable procedure | designed together | name and description always; body on use |
| `workflows/<name>/WORKFLOW.md` + `flowchart.html` | procedural: one workflow | designed together | on demand |
| `project/LOG.md` | episodic: done, remaining, bottlenecks | `wrap-up` skill | on demand |
| `README.md` | how to install and run | user | on demand |

## Two different MEMORY.md files

- **Claude Code auto-memory**:
  `~/.claude/projects/-Users-jacksonbaglino-PetIGA-3-20/memory/MEMORY.md`.
  - Claude's private notebook. Each fact is a separate file and `MEMORY.md`
    is their index. The first 200 lines (or 25 KB) of the index load every
    session.
  - Keyed by git repository, so `enceladus_DSM` sessions share it.
  - Machine-local; not part of the project record.
  - Toggle with `/memory`.
- **Project memory**: `memory/MEMORY.md`. Committed, reviewed, imported by
  `CLAUDE.md`. Survives a fresh clone and a cloud run. This is the one the
  course grades.

## Checkpoints: three mechanisms

| Mechanism | What it saves | Use |
|---|---|---|
| `/rewind` (Esc Esc) | Claude's file edits and the conversation, before each prompt | undo a wrong turn. Does not track Bash `mv`/`rm`, subagent edits or external changes; git commits are the real safety net |
| `claude --resume <id>`, `--continue` | the whole session | reopen a closed session (IDs are in `.claude/sessions.tsv`) |
| `runs.db` (SQLite) | workflow state after every step | resume, replay and audit a multi-agent run |

### runs.db

A Claude Code session remembers only its own conversation. A workflow is
several agents, each in a fresh context, so its state must live outside them,
in a file. `runs.db` is that file. It is a SQLite database, one per workflow,
with one row per agent call:

| column | example |
|---|---|
| run_id | `docs-synth-2026-10-23-a` |
| stage, agent, model | `extract`, `doc-reader`, `claude-sonnet-5-5` |
| input_hash, output_path | sha256 of inputs; `runs/.../stage2/claims.json` |
| verdict, message | `FAIL`, `eq. 4 sign differs from assembly.c:99` |
| started, finished | timestamps |

It gives three things:

1. **Resume.** Rows already marked passed are skipped; the run restarts at
   the first unfinished stage.
2. **Audit.** Which model saw which input and what the check said.
3. **Leakage control.** One `run_id` per run, so one run's state never reaches
   another (Lecture 6: one thread per revision).

LangGraph's `SqliteSaver` is a ready-made version of this. A small Python
driver that writes the rows itself does the same job without the framework.
The Lecture 6 advice (always SQLite; files that encode the memory table) is a
rule of thumb, adopted here because it is cheap and auditable.

## When to /clear or /compact

- Between unrelated tasks.
- After correcting the same mistake twice; the context now holds the wrong
  idea.
- When `/context` shows the window mostly full of old tool output.

Before either, say so; Claude runs `wrap-up`. State that is not in a file is
lost. `/compact <focus>` keeps a summary; `/clear` keeps nothing but files.
