# Claude Code Instructions — Permafrost Project

## Git Workflow

**Commit autonomously after every meaningful code change.** You do not need to ask permission first — committing is the default. Use your judgement on what is worth a commit: a code change, a bug fix, a refactor, a new feature, a new file all warrant a commit. A single-parameter tweak in an opts file usually does not unless it represents a deliberate design decision. When in doubt, commit.

This project's policy overrides the harness default of "never commit unless explicitly asked".

### Commit guidelines

- Stage only the files directly changed by the current task (not unrelated files).
- Write descriptive commit messages that explain *why* the change was made, not just what changed.
- Keep commits focused: one logical change per commit.
- Always check `git status` and `git diff` before committing so nothing unexpected is staged.

### When to commit

- After fixing a bug
- After refactoring a function or file
- After adding a new feature or initial condition
- After creating or updating a post-processing script
- After updating documentation (NOTES.md, TESTS.md, README, etc.)
- After any change that makes the code compile and run correctly

### Format

Commit messages should follow this style:

```
Short imperative summary (≤ 72 chars)

Optional body explaining the motivation or context.
```

### Do not skip hooks

Never use `--no-verify`. If a pre-commit hook fails, fix the underlying issue.

### Do not force-push

Never force-push to `main`. Feature branches are fine.

---

## Project context

- Language: C (PETSc / PetIGA framework), Python (post-processing)
- Build system: Make (via PetIGA)
- Main branch: `main`; active development on feature branches

## Activity log

**Before ending every session**, prepend a new entry to `.claude/ACTIVITY_LOG.md`. Format:

```
## YYYY-MM-DD — Short title

- Bullet summary of what was done and why.
- One bullet per meaningful task or change.

---
```

- Newest entries go at the top (prepend, don't append).
- Keep summaries factual and brief — focus on *what changed* and *why*.
- If nothing meaningful was done (e.g., a read-only Q&A session), still add a one-liner entry.

---

## Code style

- C99: declare loop variables inside `for (PetscInt i = ...)`.
- No C99 variable-length arrays (VLAs) — use `PetscMalloc`/`PetscFree`.
- Always `CHKERRQ(ierr)` immediately after every PETSc/MPI call.
- Prefer `PetscMin`/`PetscMax`/`PetscSqrtReal`/`PetscTanhReal` over bare C math functions.
