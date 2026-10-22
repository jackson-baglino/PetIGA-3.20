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

Agent names follow the model that runs them: `rick-*` Claude Fable, `beth-*`
Claude Opus, `morty-*` Claude Sonnet, `jerry-*` Claude Haiku; the exact model
ID is in each agent's file under `.claude/agents/`.

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
│  ├─ Studio/        # local runners (run_lunar.sh, run_batch_tests.sh)
│  ├─ HPC/           # SLURM submit/run scripts
│  ├─ lib/alloc.sh   # single source of truth for allocation constants
│  └─ check_ic_types.sh   # guard: validates every .opts -ic_type vs the solver
├─ studies/          # study folders (e.g. contact_angle/) with their run lists and results
├─ docs/             # design/analysis notes (see per-file banners; contact_angle/ is current)
├─ memory/           # committed project memory: MEMORY.md, USER.md, semantic/ (see memory/README.md)
├─ workflows/        # agentic workflow designs (TEMPLATE.md; one folder per workflow)
├─ project/          # course project: PILLARS.md, LOG.md, guides/
├─ .claude/          # settings.json (pinned model, permissions, hooks), skills/, agents/
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
from the solver's because Marker pulls in PyTorch:

```bash
python3.11 -m venv venv_workflows
venv_workflows/bin/pip install -r requirements-workflows.lock
```

`requirements-workflows.txt` lists the non-Python tools too (Claude Code,
Graphviz, MuPDF, Chrome).

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
