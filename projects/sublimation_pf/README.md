# sublimation_pf — PETSc + PetIGA phase-field model

A finite-element / isogeometric (**PetIGA** on **PETSc 3.20**) phase-field model
of **sublimation-driven ice metamorphism**: coupled evolution of an ice phase
field, temperature, and water-vapor density. This is the **master model**,
specialized per application under `studies/`.

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
sublimation_pf/
├─ src/              # solver: permafrost2.c (main), assembly.c (residual/Jacobian),
│                    #   initial_conditions.c, material_properties.c, monitoring.c,
│                    #   snes_convergence.c, env_helper.c
├─ include/          # headers (NASA_types.h holds AppCtx + Field)
├─ makefile          # `make` (optimized) / `make debug`; builds ./permafrost
├─ inputs/
│  ├─ solver.opts            # numerical/model defaults (-dof 3, xi_v, xi_T, bounds, ...)
│  ├─ geometry/<name>.opts   # mesh, domain, IC (-ic_type), eps, delt_t
│  └─ experiment/<name>.opts # -t_final, -temp, -humidity, -grad_temp0, cadence
├─ preprocess/       # comp_eps.py (parameter engine), build_geometry_*.py, ...
├─ postprocess/      # plot_mass.py, make_movie.py, neck_width.py, ...
├─ scripts/
│  ├─ Studio/        # local runners (run_permafrost.sh, run_batch_tests.sh)
│  ├─ HPC/           # SLURM submit/run scripts
│  ├─ lib/alloc.sh   # single source of truth for allocation constants
│  └─ check_ic_types.sh   # guard: validates every .opts -ic_type vs the solver
├─ studies/          # per-paper studies (see each README)
│  ├─ icy_regolith/         # Paper 1: implicit_pore_domain + explicit_sediment_phase
│  └─ snow_thermal/         # Paper 2: DSM on packings + effective conductivity
├─ docs/             # design/analysis notes (HISTORICAL — see per-file banners)
└─ _trash/           # files staged for deletion pending review (see MANIFEST.md)
```

## Build

```bash
export PETSC_DIR=/path/to/petsc PETSC_ARCH=<arch> PETIGA_DIR=/path/to/petiga
make            # optimized; produces ./permafrost
make debug      # -g3 -O0
```

## Run

Never invoke the binary by hand — use the run script, which assembles the three
opts files, sizes the rank count, and stages a reproducible copy of the run:

```bash
./scripts/Studio/run_permafrost.sh <geometry> <experiment> [tag] [-- extra -flags]
# e.g.
./scripts/Studio/run_permafrost.sh 2D_two_ice_grains_boundary 1day_T-20_h1.00
```

Geometry and experiment name files in `inputs/geometry/` and
`inputs/experiment/` (without the `.opts` suffix). Extra args after the tag (or
a literal `--`) are forwarded to the executable and override the opts files.
Output lands under `~/SimulationResults/sublimation_pf/scratch/<geom>/<ts>_<exp>[_tag]/`.

On HPC, `scripts/HPC/submit_permafrost.sh` computes the allocation
(`TARGET_DOFS_PER_CORE` from `scripts/lib/alloc.sh`, default 50k; `--half-cores`
halves it) and submits via `sbatch`.

## Parameters

Always (re)compute the interface width ε and mesh from the domain / grain sizes
/ temperature with `preprocess/comp_eps.py` (Kaempfer & Plapp 2009 bounds). ε is
a **decay-length scale**, not the visible diffuse-band width — see
`docs/interface_width_conventions.md`.
