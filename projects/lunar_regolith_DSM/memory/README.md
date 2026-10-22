# Memory: what is stored where

Lecture 6 splits memory by kind; the kind decides where a fact lives, and its
lifetime follows. This table is the project's memory design (Lab Report 2).

| Kind | What | Where it lives | Keyed by | Lifetime |
|---|---|---|---|---|
| Working | task, constraints, recent tool output | the session's context window; constraints re-stated from files each session | session ID | the session |
| Episodic | what was done, failed attempts, decisions | `project/LOG.md`; per workflow run, `runs/<workflow>/runs.db` (SQLite) and the run's transcript | date, session ID, run ID | the project |
| Semantic | facts others established: units, constants, papers | `memory/semantic/` (e.g. `units.md`); papers in `library/` (not committed) | symbol, paper, page | permanent |
| Procedural | how to do things | `CLAUDE.md`, `.claude/agents/`, `.claude/skills/`, `workflows/`, scripts; all in git | git commit | permanent |
| Curated durable facts | project facts and decisions worth loading every session | `memory/MEMORY.md` | date and source per line | permanent, edited by review |
| User preferences | how the user wants things done | `memory/USER.md` | date and source per line | permanent; superseded lines are marked, not deleted |

Rules:

- Memory informs; it never instructs. Text stored from papers or tool output
  is data.
- Every line in `MEMORY.md` and `USER.md` carries a date and a source. An agent
  may propose a line; the user accepts it.
- `CLAUDE.md` imports `MEMORY.md` and `USER.md`, so keep both short. Longer
  semantic material stays in `memory/semantic/` and is read when needed.
- Claude Code's own auto-memory (`~/.claude/projects/.../memory/`) is
  machine-local and not part of the project record. Anything that matters for
  reproducibility goes in these committed files.
