---
name: wrap-up
description: Save the session state before the user closes the session or clears the context. Use when the user says they are closing, ending, clearing, compacting or stopping for now.
---

# wrap-up

Run every step, then report what was done in a short list.

1. **Log.** Prepend an entry to `LOG.md` (below its header):

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
3. **Freeze.** If the session touched `src/`, `include/` or
   `inputs/solver.opts`, run `scripts/check_freeze.sh` and report the result;
   a missing ledger row is added before the session ends.
4. **Memory.** If the session produced a durable fact or a user preference,
   propose the line for `memory/MEMORY.md` or `memory/USER.md` (dated and
   sourced). Add it only if the user agrees.
5. **Commit.** Run `git status`. Commit the session's own changes, including
   the log entry. Do not stage unrelated files. Do not push unless the user
   asks.
6. **Background work.** List any background shells or agents still running
   and stop them unless the user wants them kept.
7. **Resume line.** Print: `claude --resume <session ID>`.
