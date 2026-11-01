# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

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

- Language: C (PETSc 3.20 / PetIGA framework), Python (pre/post-processing)
- Build system: Make (via PetIGA)
- Main branch: `main` (single branch; older work lives under `archive/*` tags)
- Phase-field model of sublimation-driven ice metamorphism (dry snow on
  Enceladus); `README.md` has the physics summary. Parts of its layout section
  are stale (e.g. `studies/snow_thermal/` is now `studies/keff_*`).

## Commands

```bash
make                         # optimized build -> ./enceladus_dsm (needs PETSC_DIR, PETSC_ARCH, PETIGA_DIR)
make debug                   # -g3 -O0; same warnings (-Wall -Wextra -Wshadow)
./scripts/check_ic_types.sh  # every 2D/3D .opts -ic_type is in the solver dispatch; -geom_file opts carry "# DOF_GRID:"
./scripts/Studio/run_enceladus.sh <geom> <exp> [tag] [-- -flag val ...]   # local run (see rules below)
bash postprocess/run_postprocess.sh <run_folder>                          # regenerate a run's figures
source venv_enceladus/bin/activate   # Python env for preprocess/ and postprocess/ (requirements.txt)
```

There is no unit-test suite. Correctness gates are the committed driver
scripts under `studies/*/verification/` (e.g.
`studies/restart/verification/verify_restart.sh`,
`studies/keff_sharp_limit/verification/verify_keff_sharp_limit.sh`), each
with a README. A Jacobian change must be checked with `-snes_test_jacobian`
passed through the run script.

## Architecture

- **Options layering.** Every run is `inputs/solver.opts` (numerics/model
  defaults) + a geometry `.opts` (mesh, domain, `-ic_type`, `-eps`,
  `-delt_t`) + an experiment `.opts` (`-t_final`, `-temp`, `-humidity`,
  `-grad_temp0`, output cadence), later files overriding earlier ones, then
  any flags after `--`. Scripts take bare names; `scripts/lib/opts.sh`
  (`resolve_opts`) finds them in family subdirectories and refuses ambiguous
  names. `scripts/lib/alloc.sh` is the only place rank-sizing constants live.
- **Solver flow (`src/enceladus_main.c`).** Reads options, builds the IGA
  (or loads `-geom_file`), dispatches the IC by `-ic_type` (implementations in
  `initial_conditions.c`), and integrates with `TSALPHA`. Residual and
  analytic Jacobian are in `assembly.c`; constitutive laws in
  `material_properties.c`; the per-DOF SNES convergence test in
  `snes_convergence.c`. Monitors (`monitoring.c`) handle the interface-CFL
  adaptive timestep, the diagnostics table, field output, and k_eff
  sampling; a pre-step hook does phase-field bounds rollback. Restart in
  place resumes from a snapshot (`scripts/HPC/resume_batch.sh`).
- **AppCtx** (`include/enceladus_types.h`) carries all parameters and state;
  new options are parsed in `enceladus_main.c` and stored there.
- **In-line k_eff** (`-keff 1`, `keff*.c`, `keff_` options prefix): solves
  the periodic cell problem on the current phi and appends to
  `k_eff.csv`; `-keff_replay <run_dir>` recomputes it for a finished run.
  Hard-errors on non-periodic meshes, `-geom_file`, `-axisym`, `dim < 2`.
- **Sibling project.** `../lunar_regolith_DSM` was split from the same
  codebase. `assembly.c`, `material_properties.c` and `snes_convergence.c`
  started identical but have since diverged (lunar has the three-phase
  solver). When fixing a bug in shared physics, check whether the sibling
  needs the same fix and ask before changing it.
- **Run folder contents.** The run script stages the three opts files,
  `src/` and the postprocess tree into the run folder, so a downloaded run
  can be re-plotted with `run_postprocess.sh` and no source tree.
  `postprocess/pplib.py` is the shared loading library for figure scripts.
- **HPC.** SLURM scripts in `scripts/HPC/`; two or more jobs go through
  `submit_batch.sh` (not chained single submits). Results are pulled back
  with `fetch_runs.sh` / `fetch_stage.sh`.

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

## k_eff manuscript campaign

The live to-do list is `studies/keff_sintering/TODO.md` -- read it at the
start of any session on this effort. `studies/keff_sintering/RECORD.md` is
the history for writing: whenever a task or run finishes (or a plan
changes), add its HPC + local output dirs and a one-sentence result there,
in the same commit. Manuscript runs are submitted ONLY with
`scripts/HPC/submit_keff_production.sh <stage file>`, which owns every run
option; stage files list runs, never options.

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
