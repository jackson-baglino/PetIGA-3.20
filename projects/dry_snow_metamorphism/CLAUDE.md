# CLAUDE.md — dry_snow_metamorphism

## Commit policy

Commit autonomously when you complete a logical unit of work. Do not wait to be asked unless the change is destructive or irreversible (branch deletes, force pushes, dropping files that exist nowhere else).

## Project overview

Phase-field model of dry snow metamorphism following Kaempfer & Plapp (2009),
implemented with PetIGA (isogeometric analysis). Three degrees of freedom per node:

| DOF | Name | Description |
|-----|------|-------------|
| 0 | φ_i | Ice phase field ∈ [0, 1] |
| 1 | T | Temperature [°C relative to T₀] |
| 2 | ρ_v | Vapor density [kg/m³] |

No sediment DOF. Keep the three-DOF structure if you add anything.

## Build

```bash
cd /Users/jacksonbaglino/PetIGA-3.20/projects/dry_snow_metamorphism
make dry_snow_metamorphism
```

Requires `PETSC_DIR`, `PETSC_ARCH`, and `PETIGA_DIR` environment variables.
IDE linters (clangd) will show false positives for missing headers — the build
itself is correct; only `make` knows the include paths.

## Running

Parameters are passed via PETSc `.opts` files, NOT environment variables:

```bash
mpiexec -np 12 ./dry_snow_metamorphism \
  -options_file inputs/solver.opts \
  -options_file inputs/geometry/grains__phi=0.24__Lxmm=2__Lymm=2__seed=7/grains.opts \
  -options_file inputs/experiment/30day_T-10_h1.00.opts \
  -output_dir /path/to/results
```

Or use the run scripts:
```bash
./scripts/studio/run_dsm.sh \
  -g inputs/grains__phi=0.24.../grains.opts \
  -e inputs/experiment/30day_T-10_h1.00.opts
```

For HPC (SLURM):
```bash
GEOMETRY_OPTS=inputs/.../grains.opts \
EXPERIMENT_OPTS=inputs/experiment/30day_T-10_h1.00.opts \
sbatch scripts/HPC/run_dsm.sh
```

## Parameter system

Three-level `.opts` hierarchy:

```
inputs/
├── solver.opts           # Numerics and solver settings — do not change per run
├── geometry/             # Mesh + grain geometry (one file per grain configuration)
│   └── *.opts
└── experiment/           # Physics + time (temperature, humidity, duration)
    └── *.opts
```

Later files override earlier ones (PETSc merges all sources). You can also
append `-flag value` to the command line to override anything.

## CRITICAL: eps is a decay-length, not the interface width

`eps` (`-eps`) is the exponential-decay-length scale of the tanh interface profile.
The *visible* diffuse band (5%–95% of φ_i) spans **≈ 6–9× eps**, not 1× eps.

Use `preprocess/comp_eps.py` to compute the correct value:
```bash
python3 preprocess/comp_eps.py --Lx 1e-3 --Ly 1e-3 --dim 2 --Rave 3e-5 --T0 -10
```

Never set `eps` to the grain radius or domain size directly.

## CRITICAL: Do NOT scale source terms with xi_T / xi_v

The flags `-xi_T` and `-xi_v` are exposed for diagnostics only. Using them to
scale the latent-heat or vapor-mass source terms **breaks conservation** and
produces unphysical results. They must remain at their default values (1.0).

## Solver

Current solver settings are in `inputs/solver.opts`:
```
-snes_rtol 1e-3 -snes_stol 1e-6 -snes_max_it 7
-ksp_gmres_restart 150 -ksp_max_it 1000
-snes_linesearch_type basic
```

Do not change these without benchmarking against a known-good case.
The time stepper is `TSALPHA` with adaptive stepping.

## Postprocessing

```bash
# Convert sol_*.dat → VTK + write dsm.pvd time series
python3 postprocess/plotDSM.py --dir /path/to/run

# Plot scalar time series (SSA_evo.dat)
python3 postprocess/plot_scalars.py --file /path/to/run/SSA_evo.dat --time-unit d

# Plot mass conservation
python3 postprocess/plot_mass.py --dir /path/to/run

# Compare multiple runs
python3 postprocess/compare_runs.py dir1 dir2 --labels "T=-20C" "T=-10C"

# Render movie (must use pvpython, not regular Python)
/Applications/ParaView-5.12.0.app/Contents/bin/pvpython \
    postprocess/make_movie.py --dir /path/to/run
```

## Preprocessing

```bash
# Compute eps and mesh sizing from grain radius
python3 preprocess/comp_eps.py --Lx 1e-3 --Ly 1e-3 --dim 2 --Rave 3e-5 --T0 -10

# Generate grains.dat + grains.opts from scratch
python3 preprocess/generateCircularGrains.py \
    --Lx 1e-3 --Ly 1e-3 --porosity 0.3 --mean-r-m 5e-5

# Generate grains.opts from an existing grains.dat
python3 preprocess/generate_opts_from_input.py grains.dat grains.opts 1e-3 1e-3
```

## File structure

```
dry_snow_metamorphism/
├── src/
│   ├── dry_snow_metamorphism.c   # Main: IGA setup, time stepping
│   ├── options_helper.c          # GetOptions() — PETSc options parsing
│   ├── monitoring.c              # Monitor(), OutputMonitor()
│   ├── grain_initialization.c    # Grain geometry setup
│   └── ...
├── include/
│   ├── NASA_types.h              # AppCtx struct definition
│   ├── NASA_main.h               # Residual/Jacobian declarations
│   └── options_helper.h          # GetOptions() declaration
├── inputs/
│   ├── solver.opts
│   ├── geometry/
│   └── experiment/
├── preprocess/
│   ├── comp_eps.py               # eps and mesh sizing calculator
│   ├── generateCircularGrains.py # Grain generator (writes grains.dat + grains.opts)
│   └── generate_opts_from_input.py
├── postprocess/
│   ├── plotDSM.py                # sol_*.dat → VTK + dsm.pvd
│   ├── plot_scalars.py           # SSA_evo.dat time series
│   ├── plot_mass.py              # Mass conservation check
│   ├── compare_runs.py           # Multi-run overlay
│   └── make_movie.py             # ParaView movie (use pvpython)
└── scripts/
    ├── studio/run_dsm.sh         # Local run script
    └── HPC/run_dsm.sh            # SLURM job script
```
