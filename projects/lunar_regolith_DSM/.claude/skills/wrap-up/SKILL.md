---
name: wrap-up
description: Save the session state before the user closes the session or clears the context. Use when the user says they are closing, ending, clearing, compacting or stopping for now.
---

# wrap-up

Run every step, then report what was done in a short list.

1. **Log.** Prepend an entry to `project/LOG.md`:

   ```
   ## YYYY-MM-DD — short title
   - Session: <session ID from the SessionStart context, or "unknown">
   - Done: ...
   - Remaining: ...
   - Bottlenecks: ...
   ```

2. **Open work.** If a task is unfinished, write its next concrete step and
   the files involved under "Remaining", so a fresh session can resume from
   the log alone.
3. **Memory.** If the session produced a durable fact or a user preference,
   propose the line for `memory/MEMORY.md` or `memory/USER.md` (dated and
   sourced). Add it only if the user agrees.
4. **Prune memory.** Read the auto-memory index (`MEMORY.md` in the
   auto-memory folder) and its project and status entries. List every entry
   whose end condition is met (the task is done, the plan it points to is
   finished, or the fact is superseded), with the evidence (a commit, a doc
   section, a date). Propose deleting each one; delete only what the user
   agrees to, and remove its line from the index. Any new status entry gets
   an end condition. In `memory/MEMORY.md` (in git), superseded lines are
   marked, not deleted (CLAUDE.md). (User, 2026-11-17.)
5. **Commit.** Run `git status`. Commit the session's own changes, including
   the log entry. Do not stage unrelated files. Do not push unless the user
   asks.
6. **Background work.** List any background shells or agents still running
   and stop them unless the user wants them kept.
7. **Resume line.** Print: `claude --resume <session ID>`.
