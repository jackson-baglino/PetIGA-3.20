# lunar_regolith_DSM — PETSc + PetIGA phase-field model

A finite-element / isogeometric (**PetIGA** on **PETSc 3.20**) phase-field model
of **sublimation-driven ice metamorphism**: coupled evolution of an ice phase
field, temperature, and water-vapor density.

**Application:** ice in **lunar regolith** — sublimation and redistribution of
ice within a regolith pore domain.

This project and its sibling `enceladus_DSM` were split on 2026-07-27 from a
single codebase (`sublimation_pf`, formerly `permafrost`) so that two papers
can evolve independently. Their solvers have since diverged; only
`src/snes_convergence.c` is still identical.

This folder is also the AM/ME/CE 295 course project (AI agents for scientific
discovery): agentic workflows that read and check the literature, and will
derive, implement and test a new coupling between the fixed regolith grains and
the evolving ice (see "Agentic workflows" below).

Agent names follow the model that runs them: `rick-*` Claude Fable, `beth-*`
Claude Opus, `morty-*` Claude Sonnet, `jerry-*` Claude Haiku. The exact model
ID is pinned in each agent's file: `workflows/<workflow>/agents/` for the
agents of a workflow, `.claude/agents/` for general helpers.

## Physics (two-phase)

Three degrees of freedom per node — ice φ, temperature T, vapor density ρ_v
(air fraction is algebraic, φ_a = 1 − φ). Allen–Cahn ice evolution with a plain
double-well, a localized sublimation source, latent-heat coupling in the
temperature equation, and a vapor-transport equation with Moure & Fu 2024
temporal scaling (ξ_v applied to diffusion *and* source together; ξ_T = 1).
Optional axisymmetric r–z mode, phase-field bounds enforcement, and an
interface-CFL adaptive timestep. (Gibbs–Thomson curvature was removed 2026-07-21.)

## Layout

```
lunar_regolith_DSM/
├─ src/              # solver: lunar_main.c (main), assembly.c (residual/Jacobian),
│                    #   initial_conditions.c, material_properties.c, monitoring.c,
│                    #   snes_convergence.c, env_helper.c
├─ include/          # headers (NASA_types.h holds AppCtx + Field)
├─ makefile          # `make` (optimized) / `make debug`; builds ./lunar_regolith_dsm
├─ inputs/
│  ├─ solver.opts            # numerical/model defaults (-dof 3, xi_v, xi_T, bounds, ...)
│  ├─ geometry/<name>.opts   # mesh, domain, IC (-ic_type), eps, delt_t
│  └─ experiment/<name>.opts # -t_final, -temp, -humidity, -grad_temp0, cadence
├─ preprocess/       # comp_eps.py -> mesh_params/ (eps and mesh engine, tests in tests/),
│                    #   build_geometry_*.py
├─ postprocess/      # plot_mass.py, make_movie.py, neck_width.py, ...
├─ scripts/
│  ├─ literature_review/   # the B1 line: run_b1.py, foreman.py, stations, tests/
│  ├─ units_verification/  # the units-verification workflow, tests/
│  ├─ agentlib.py    # shared runner for headless agents (schema, transcript audit, retry, ledger)
│  ├─ flowchart/     # render_flowchart.py: workflow flowcharts from flowchart.json
│  ├─ Studio/        # local runners (run_lunar.sh, run_batch_tests.sh)
│  ├─ HPC/           # SLURM submit/run scripts
│  ├─ lib/alloc.sh   # single source of truth for allocation constants
│  ├─ transcripts/   # short_transcript.py: prompts + final replies of a session (course assignments)
│  └─ check_ic_types.sh   # guard: validates every .opts -ic_type vs the solver
├─ studies/          # study folders (e.g. contact_angle/) with their run lists and results
├─ docs/             # design/analysis notes (see per-file banners; contact_angle/ is current)
├─ memory/           # committed project memory: MEMORY.md, USER.md, semantic/ (see memory/README.md)
├─ workflows/        # agentic workflow designs (TEMPLATE.md; one folder per workflow:
│                    #   WORKFLOW.md, flowchart.json/.html, agents/, prompts/, schemas/)
├─ runs/             # workflow run records (manifest.json, runs.db, reports); paper text stays local
├─ project/          # course project: PILLARS.md, LOG.md (running log), RETROSPECTIVE.md, guides/
├─ .claude/          # settings.json (pinned model, permissions, hooks), hooks/, rules/,
│                    #   skills/, agents/ (general helpers), permissions/
└─ PERMISSIONS.md    # what agents may and may not do, and what enforces it
```

## Build

```bash
export PETSC_DIR=/path/to/petsc PETSC_ARCH=<arch> PETIGA_DIR=/path/to/petiga
make            # optimized; produces ./lunar_regolith_dsm
make debug      # -g3 -O0
```

## Workflow tools

The agentic workflows (`workflows/`) need their own environment, kept apart
from the solver's because the PDF tools pull in PyTorch:

```bash
python3.11 -m venv venv_workflows
venv_workflows/bin/pip install -r requirements-workflows.lock
```

`requirements-workflows.txt` lists the non-Python tools too (Claude Code,
Graphviz, MuPDF, Chrome). The literature review's PDF parser, MinerU2.5-Pro,
runs offline in a third environment, `venv_docparse/` (called by
`scripts/literature_review/prepare.py`); it has no lock file yet. Marker, the
earlier parser, is still used by the units-verification workflow.

## Agentic workflows

Each workflow is designed with a flowchart before it is built
(`workflows/TEMPLATE.md`). Scripts manage and agents work: a line script walks
the stage graph in `flowchart.json`, station scripts give each headless agent
one task, validate its JSON output and record every call, and the user decides
at gates.

| Workflow | Status | Entry point |
|---|---|---|
| `workflows/units-verification/` | run 2026-10-22; `memory/semantic/units.md` verified | `scripts/units_verification/` |
| `workflows/literature-review/` (branch B1: nine fixed papers on coupling stationary grains with an evolving interface) | development phase; run contract `B1_CONTRACT.md` | `scripts/literature_review/run_b1.py` |

```bash
venv_workflows/bin/python scripts/literature_review/run_b1.py status      # where the line stands
venv_workflows/bin/python scripts/literature_review/foreman.py            # foreman session (beth-foreman)
venv_workflows/bin/python scripts/literature_review/lr_status.py --brief  # one-line status
```

- **Phases.** Develop (anything may change; runs tagged `dev`), freeze (the
  user approves; hashes of every definition and a git tag), execute (the line
  refuses changed files; ideas go to `BACKLOG.md`).
- **Reproducibility.** Every run writes `runs/<workflow>/<run_id>/manifest.json`
  (git commit, Claude Code version, model ID per agent, hashes of definitions
  and inputs) and a `runs.db` ledger with one row per agent call. Models are
  pinned by full ID.
- **Permissions and untrusted text.** Paper text is untrusted: a normal
  session cannot read `library/`, `index/` or `runs/literature-review/**/untrusted/`
  (`.claude/settings.json`); agents' permissions are in
  `.claude/permissions/literature-review/` and summarized in `PERMISSIONS.md`.
  Copyrighted papers and page images are never committed.
- **History.** `project/LOG.md` (running log), `project/RETROSPECTIVE.md`
  (pivots and lessons), `project/PILLARS.md` (project plan).

## Tests

```bash
venv_workflows/bin/python -m unittest discover -s scripts/literature_review/tests
venv_workflows/bin/python -m unittest discover -s scripts/units_verification/tests
venv_lunar/bin/python -m unittest discover -s preprocess/tests -t .
venv_lunar/bin/python -m unittest scripts/flowchart/test_render_flowchart.py
```

There is no C unit-test suite; solver checks are driver scripts under
`studies/*/verification/`, and Jacobian changes are checked with
`-snes_test_jacobian` passed through the run script.

## Run

Never invoke the binary by hand — use the run script, which assembles the three
opts files, sizes the rank count, and stages a reproducible copy of the run:

```bash
./scripts/Studio/run_lunar.sh <geometry> <experiment> [tag] [-- extra -flags]
# e.g.
./scripts/Studio/run_lunar.sh twograins_2D_L41um_eps0.49um base_T-20_h1.00_1d
```

Geometry and experiment name files in `inputs/geometry/` and
`inputs/experiment/` (without the `.opts` suffix). Extra args after the tag (or
a literal `--`) are forwarded to the executable and override the opts files.
Output lands under `~/SimulationResults/lunar_regolith_DSM/scratch/<geom>/<ts>_<exp>[_tag]/`.

On HPC, `scripts/HPC/submit_lunar.sh` computes the allocation
(`TARGET_DOFS_PER_CORE` from `scripts/lib/alloc.sh`, default 50k; `--half-cores`
halves it) and submits via `sbatch`.

## Parameters

Always (re)compute the interface width ε and mesh from the domain / grain sizes
/ temperature with `preprocess/comp_eps.py` (Kaempfer & Plapp 2009 bounds). ε is
a **decay-length scale**, not the visible diffuse-band width — see
`docs/interface_width_conventions.md`.
