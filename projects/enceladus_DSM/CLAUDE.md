# Claude Code Instructions — enceladus_DSM

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

## Running tests and experiments

These two rules apply to **every** run — quick probes and one-off diagnostics
included, not just production simulations.

### 1. Always go through the run script

Use `./scripts/Studio/run_enceladus.sh <geom> <exp> [tag] [-- extra flags]`
(or `scripts/HPC/submit_enceladus.sh` on the cluster). **Never invoke
`./enceladus_dsm` directly for anything whose result will be looked at.**

The script is not a convenience wrapper — it is what makes a run reproducible
and legible. It compiles, sizes the rank count, creates a timestamped output
folder, **stages the three `.opts` files and a copy of `src/` into it**, runs,
and then **generates the plots**. A hand-rolled `./enceladus_dsm ...` skips all
of that: no staged inputs, no provenance, no figures.

Exceptions, and they are narrow: a build smoke test that produces no result
worth reading, and unit-style gates under `studies/*/verification/` that have
their own committed driver script (those scripts *are* the reproducible record,
and they write their CSVs and figures into the repo).

### 2. Always save the results where they can be seen

Results must survive the session and must be findable. A run whose output landed
in `/tmp`, a scratchpad, or a shell variable that has since gone away did not
happen.

- Simulation output goes to `$RESULTS_BASE` via the run script — that is
  `/Users/jacksonbaglino/SimulationResults/enceladus_DSM/scratch/<geom>/<timestamp>_<exp>[_tag]/`.
- Verification and study results go **into the repo**, under
  `studies/<study>/verification/` or equivalent: the driver script, the CSV, the
  figures, and a README explaining what the numbers mean. Raw solver logs too if
  they are small.
- **Produce plots.** A CSV nobody can look at is half an answer. If the run
  script's default figures do not cover what the test is about, add a plotting
  script next to the driver and commit both.
- Say where the results are, explicitly, when reporting them.

---

## Input files: search `inputs/scratch/` before creating a new one

`inputs/` holds only `solver.opts` and `README.md` at its top level. Every
geometry file, experiment file, packing and mesh from earlier work lives in
`inputs/scratch/`, in its original sub-directory layout.

**Before writing any new input file, search the quarantine for one that already
does the job:**

```bash
find inputs/scratch -name '*<something>*'
```

- **If it exists**, pull it back out and use it — never rewrite an equivalent
  file under a new name:
  `git mv inputs/scratch/geometry/<family>/<name>.opts inputs/geometry/<family>/<name>.opts`
- **If it does not**, create a new one in the live tree.

Either way, keep the existing sub-directory scheme: geometry under
`geometry/<family>/`, experiments under `experiment/<family>/`. The run scripts
run this search automatically and print the exact `git mv` when a name misses.

The same quarantine convention applies to `preprocess/scratch/`,
`postprocess/scratch/` and `studies/scratch/`. All of them are **pending
deletion**, not archives — treat anything in them as gone unless it is pulled
back into the live tree.

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
