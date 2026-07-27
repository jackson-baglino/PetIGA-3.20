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

---

## Phase-field interface parameter (eps)

`eps` (`-eps` in opts files, `user->eps` in code) is the **decay-length scale**
of the diffuse-interface profile — it is NOT the width of the diffuse band you
see in a ParaView contour plot. Conflating the two leads to "fixing" eps by the
wrong multiplicative factor.

The equilibrium 1D profile of this model's double well is logistic:
`phi(x) = 1/(1 + exp(-x/eps))`. Its tails are long, so the band that visibly
looks diffuse on screen spans several multiples of eps, not eps itself:

- 5%–95% transition  ≈ 6·eps
- 1%–99% transition  ≈ 9.2·eps

So "N elements visibly diffuse in ParaView" corresponds to roughly
`eps ≈ N·dx / 6` to `N·dx / 9` (depending on how sharp a cutoff you're
eyeballing) — **not** `eps ≈ N·dx`.

**There is also a second, easy-to-confuse "width": the Karma-Plapp convention
`w_karma = 2√2·eps`** used inside `comp_eps.py` to size `dx` from a target
`n_per_interface`. `w_karma` is ~3.25x *narrower* than the 1%-99% band
(`w_karma · ln(99)/√2 ≈ w_1_99`) — it is a pure analytic device for relating
`eps` to physical material parameters in the sharp-interface-limit
derivation, NOT a visual or resolution-adequacy metric. **When checking
whether the mesh adequately resolves the interface, count elements across
the directly-observed phi=0.01-to-0.99 band (or measure it in ParaView the
way the user does — contour at phi=0.01 and phi=0.99, count elements between
them with Surface With Edges) — never against `n_per_interface`/Karma-width
element counts.** `comp_eps.py`'s own `--n` default (4 Karma-elements) is
calibrated to land around ~7.5-10 elements across the 1%-99% band; quoting
the Karma-element count on its own (e.g. "only 2.4 elements!") will look
alarmingly low even when the real, visible interface is adequately resolved.

Always (re)compute `eps` with `preprocess/comp_eps.py` (Kaempfer & Plapp 2009
sharp-interface bounds) for the actual domain/grain sizes/temperature in play.
Never hand-tune `eps` by visually estimating the diffuse band width, and never
reason "the interface looks 2x too wide/sharp, so halve/double eps" — that
mixes up the scale parameter with the visible band, which differ by the ~6–9x
factor above.
