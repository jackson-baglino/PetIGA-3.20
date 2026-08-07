# enceladus_DSM — PETSc + PetIGA phase-field model

A finite-element / isogeometric (**PetIGA** on **PETSc 3.20**) phase-field model
of **sublimation-driven ice metamorphism**: coupled evolution of an ice phase
field, temperature, and water-vapor density.

**Application:** dry snow metamorphism of icy granular packings, applied to the
presumed state of the surface of **Enceladus** near the tiger stripes.

This project and its sibling `lunar_regolith_DSM` share the same solver. They
were split on 2026-07-27 from a single codebase (`sublimation_pf`, formerly
`permafrost`) so that two papers can evolve independently; `src/assembly.c`,
`src/material_properties.c`, and `src/snes_convergence.c` are byte-identical
between them and should stay that way unless a change is genuinely
application-specific.

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
enceladus_DSM/
├─ src/              # solver: enceladus_main.c (main), assembly.c (residual/Jacobian),
│                    #   initial_conditions.c, material_properties.c, monitoring.c,
│                    #   snes_convergence.c
│                    # k_eff: keff.c (lifecycle/guards), keff_field.c (phi projection),
│                    #   keff_cell.c (cell problem), keff_solve.c, keff_sample.c
├─ include/          # headers (enceladus_types.h holds AppCtx + Field)
├─ makefile          # `make` (optimized) / `make debug`; builds ./enceladus_dsm
├─ inputs/
│  ├─ solver.opts            # numerical/model defaults (-dof 3, xi_v, xi_T, bounds, ...)
│  ├─ geometry/<name>.opts   # mesh, domain, IC (-ic_type), eps, delt_t
│  └─ experiment/<name>.opts # -t_final, -temp, -humidity, -grad_temp0, cadence
├─ preprocess/       # comp_eps.py (parameter engine), build_geometry_*.py, ...
├─ postprocess/      # plot_mass.py, make_movie.py, neck_width.py, ...
├─ scripts/
│  ├─ Studio/        # local runners (run_enceladus.sh, run_batch_tests.sh)
│  ├─ HPC/           # SLURM submit/run scripts
│  ├─ lib/alloc.sh   # single source of truth for allocation constants
│  └─ check_ic_types.sh   # guard: validates every .opts -ic_type vs the solver
├─ studies/
│  └─ snow_thermal/  # DSM on granular packings -> effective thermal conductivity
│     └─ verification/  # k_eff analytic gate + resolution studies
└─ docs/             # design/analysis notes (HISTORICAL — see per-file banners)
```

## Effective thermal conductivity (`-keff`)

The solver computes the **effective (homogenized) thermal conductivity tensor**
in-line, by solving the periodic cell problem
`-div(k(x) ∇t_m) = div(k(x) e_m)` for each direction and averaging the resulting
flux. Merged in from `projects/effective_thermal_cond` on 2026-07-31; that
project is superseded.

It is in-line rather than a post-processing pass because the two cadences want to
be very different: a field snapshot at study resolution is hundreds of MB, so a
run writes a few hundred, while a k_eff sample is ~200 bytes. Reading k_eff off
snapshots capped k_eff(t) at the snapshot count.

```bash
# during a run, every 10 accepted steps (independent of -outp)
./scripts/Studio/run_enceladus.sh <geom> <exp> tag -- -keff 1 -keff_freq 10

# or on simulated time
... -- -keff 1 -keff_t_interv 3600

# for a run that already finished
folder=<run_dir> ./enceladus_dsm -options_file <run_dir>/solver.opts \
    -options_file <run_dir>/<geom>.opts -options_file <run_dir>/<exp>.opts \
    -keff 1 -keff_replay <run_dir>
```

Output is `k_eff.csv` in the run directory:
`step,time,k_ij…,phi_bar,k_iso,ksp_its,ksp_reason,wall_s`.

**Periodicity is required and enforced.** `-keff` hard-errors on a non-periodic
mesh, on `-geom_file`, on `-axisym`, and on `dim < 2` — each a case where the
cell problem would still return a number that means nothing.

Corrector solves default to CG + GAMG under the `keff_` option prefix; direct LU
(`-keff_ksp_type preonly -keff_pc_type lu`) is the validation reference. The
analytic gate, solver comparison and resolution studies live in
`studies/snow_thermal/verification/`.

## Build

```bash
export PETSC_DIR=/path/to/petsc PETSC_ARCH=<arch> PETIGA_DIR=/path/to/petiga
make            # optimized; produces ./enceladus_dsm
make debug      # -g3 -O0
```

## Run

Never invoke the binary by hand — use the run script, which assembles the three
opts files, sizes the rank count, and stages a reproducible copy of the run:

```bash
./scripts/Studio/run_enceladus.sh <geometry> <experiment> [tag] [-- extra -flags]
# e.g.
./scripts/Studio/run_enceladus.sh twograins_2D_L41um_eps0.49um base_T-20_h1.00_1d
```

Geometry and experiment name files in `inputs/geometry/` and
`inputs/experiment/` (without the `.opts` suffix). Extra args after the tag (or
a literal `--`) are forwarded to the executable and override the opts files.
Output lands under `~/SimulationResults/enceladus_DSM/scratch/<geom>/<ts>_<exp>[_tag]/`.

On HPC, `scripts/HPC/submit_enceladus.sh` computes the allocation
(`TARGET_DOFS_PER_CORE` from `scripts/lib/alloc.sh`, default 50k; `--half-cores`
halves it) and submits via `sbatch`.

## Parameters

Always (re)compute the interface width ε and mesh from the domain / grain sizes
/ temperature with `preprocess/comp_eps.py` (Kaempfer & Plapp 2009 bounds). ε is
a **decay-length scale**, not the visible diffuse-band width — see
`docs/interface_width_conventions.md`.
