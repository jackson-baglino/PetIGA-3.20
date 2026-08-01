# effective_k_ice_homog

> # ⚠️ SUPERSEDED — 2026-07-31
>
> **This project has been merged into `projects/enceladus_DSM` and is no longer
> the way to compute k_eff.** It is kept for provenance only. Do not add
> features here; do not use it for new results.
>
> The cell problem now runs **in-line** inside the phase-field solver:
>
> ```bash
> # during a run, on a cadence independent of field output
> ./scripts/Studio/run_enceladus.sh <geom> <exp> tag -- -keff 1 -keff_freq 10
>
> # for a run that already finished — this replaces this binary entirely
> folder=<run_dir> ./enceladus_dsm \
>     -options_file <run_dir>/solver.opts \
>     -options_file <run_dir>/<geom>.opts \
>     -options_file <run_dir>/<exp>.opts \
>     -keff 1 -keff_replay <run_dir>
> ```
>
> **Why it moved.** k_eff sampling was tied to the snapshot rate, and the two
> want very different cadences: a snapshot at study resolution is hundreds of MB
> (so a run writes a few hundred) while a k_eff sample is ~200 bytes. k_eff(t)
> was needlessly sparse as a direct consequence.
>
> **What changed in the move**, all verified in
> `enceladus_DSM/studies/snow_thermal/verification/`:
>
> - **The claim below that iterative solvers are unreliable for this problem is
>   wrong.** It was a PetIGA/GAMG plumbing incompatibility, not numerics:
>   `IGACreateMat` overrides `MATOP_CREATE_VECS`/`MATOP_DUPLICATE` with
>   IGA-aware versions that hard-error when the composed `"IGA"` object is
>   absent, and AMG's coarse operators inherit the overrides without the object,
>   so `PCGAMG` died during setup on every mesh. With three lines restoring the
>   stock operations, CG+GAMG agrees with LU to all ten printed digits and is
>   17× faster at 322k unknowns. LU was never viable at the study meshes (~8M
>   unknowns per corrector: ~90 min per sample vs ~80 s).
> - The `dim`-DOF blocked system is now `dim` solves against **one shared scalar
>   operator** — the bilinear form never referenced the direction, so three
>   quarters of the matrix was structurally zero at `dim=2` (9.6 GB → 2.4 GB at
>   Nx=2829).
> - The input IGA is no longer rebuilt from re-parsed options. `field_init.c`
>   admitted its length check could not detect a `-periodic` mismatch — the one
>   that silently gives wrong answers. Reusing the producer's own IGA makes that
>   unrepresentable.
> - The cell volume is measured by quadrature rather than assumed `Lx*Ly`.
> - Fixed on the way through: `MPI_INT` → `MPIU_INT` on a `PetscInt` broadcast,
>   an 80 KB fixed stack array for the snapshot list, `sol_%05d.dat` parsing that
>   broke past 99999 steps, `MPI_Abort` in place of `SETERRQ`, and per-call
>   `fopen`/`fclose` on the output CSV.
>
> Everything below describes the retired standalone binary.

---

Computes the **effective thermal conductivity tensor** of snow/ice microstructure
via periodic homogenization. Given a phase-field snapshot of an ice grain structure
(from a separate grain growth simulation), it solves the IGA cell problem

    -div(k(x) grad t_m) = div(k(x) e_m)   in Y,   t_m periodic

for each macroscopic direction m, then integrates to obtain the Voigt-averaged
k_eff tensor:

    k_eff[i,j] = (1/|Y|) ∫_Y k(x) (∂t_i/∂x_j + δ_ij) dV

The code uses **PetIGA** (IGA finite elements) and **PETSc** for assembly and
direct sparse linear algebra.

---

## Repository structure

```
effective_thermal_cond/
├── src/
│   ├── effective_k_ice_homog.c   main() — orchestrates the computation
│   ├── assembly.c                stiffness assembly + k_eff integrand
│   ├── solver.c                  direct LU solver with DOF pinning
│   ├── field_init.c              ice field readers and analytic microstructures
│   ├── env_config.c              PETSc options parsing (ParseOptions)
│   ├── material.c                ThermalCond() interpolation
│   ├── setup.c                   IGA object creation (SetupIGA)
│   └── io_thermal.c              CSV / binary output
├── include/
│   ├── app_ctx.h                 AppCtx struct + shared macros (SQ, CU)
│   ├── assembly.h
│   ├── solver.h
│   ├── field_init.h
│   ├── env_config.h
│   ├── material.h
│   ├── setup_thermal.h
│   └── io_thermal.h
├── inputs/
│   ├── default.opts              reference options file (all parameters + defaults)
│   └── grains_template.opts     per-simulation template (copy to sim directory)
├── scripts/
│   ├── run_effective_k_ice_homog.sh    single-case run script
│   ├── batch_effective_k_ice_homog.sh  batch over subdirectories
│   ├── gen_opts.py                     convert grains.env → grains.opts
│   └── run_hpc.slurm                  Caltech Resnick SLURM script
├── postprocess/
│   ├── plot_k_eff.py
│   ├── plot_k_eff_by_temp.py
│   ├── plot_k_eff_by_phi.py
│   ├── plotSSA.py
│   ├── core/io.py
│   └── run_all_single_plots.sh
├── archive/                    old/abandoned files kept for reference
├── makefile
├── CHANGES.md
└── README.md
```

---

## Dependencies

| Dependency | Notes |
|-----------|-------|
| PETSc >= 3.20 | with SuperLU_DIST for parallel direct solves |
| PetIGA | IGA finite element library on top of PETSc |
| MPI | OpenMPI or MPICH |
| Python >= 3.9 | for post-processing and `gen_opts.py` |

Set these environment variables before building or running:

```bash
export PETSC_DIR=/path/to/petsc
export PETSC_ARCH=arch-linux-c-opt
export PETIGA_DIR=/path/to/petiga
```

---

## Build

```bash
# Optimized build
make

# Debug build (-g3 -O0, AddressSanitizer-friendly)
make debug

# Clean
make clean
```

---

## Running

### Single case

```bash
./effective_k_ice_homog \
  -options_file inputs/default.opts \
  -init_mode    layered \
  -output_dir   /tmp/keff_out
```

`-init_mode` controls the microstructure:

| Mode | Description |
|------|-------------|
| `layered` | Analytic periodic ice stripe (50 % ice, two tanh interfaces) — good for validation |
| `circle` | Single circular ice grain centred in the domain |
| `file` | Read IGA descriptor + solution vector from `-init_dir` |

### Reading from a grain growth simulation

```bash
./effective_k_ice_homog \
  -options_file /path/to/grains.opts \
  -init_mode    file \
  -init_dir     /path/to/sim_dir \
  -output_dir   /path/to/output_dir
```

`init_dir` must contain `igasol.dat` (IGA descriptor) and `sol_NNNNN.dat`
(solution vector). The first DOF of the solution is interpreted as ice fraction phi.

### Options file

All physics and mesh parameters are controlled via a PETSc `.opts` file:

```
-dim   2
-Nx    64
-Ny    64
-Lx    1.0e-3
-Ly    1.0e-3
-eps   1.0e-4
```

Copy `inputs/grains_template.opts` to your simulation directory and fill in the
values from `grains.env` (or generate it automatically — see below).

Run `./effective_k_ice_homog -help` to see all available options with defaults.

### Generating grains.opts from grains.env (legacy)

```bash
python3 scripts/gen_opts.py /path/to/sim/grains.env > /path/to/sim/grains.opts
```

### MPI parallel run

```bash
mpiexec -np 4 ./effective_k_ice_homog \
  -options_file /path/to/grains.opts \
  -init_mode    file \
  -init_dir     /path/to/sim \
  -output_dir   /path/to/out
```

Parallel direct solves use SuperLU_DIST automatically (selected when more than
one MPI rank is detected).

### Batch run over multiple simulation directories

```bash
./scripts/batch_effective_k_ice_homog.sh /path/to/parent_dir
```

The batch script iterates all subdirectories of `parent_dir`, skips those
without `sol_*.dat` files or where `k_eff.csv` already exists, and calls
`run_effective_k_ice_homog.sh` for each.

### HPC / SLURM (Caltech Resnick)

```bash
INIT_DIR=/path/to/sim OUTPUT_DIR=/path/to/out \
  sbatch scripts/run_hpc.slurm
```

Adjust `#SBATCH` directives at the top of `scripts/run_hpc.slurm` as needed.

---

## Input file format

A simulation directory should contain:

| File | Description |
|------|-------------|
| `igasol.dat` | IGA descriptor written by the grain growth code |
| `sol_NNNNN.dat` | Solution vector at time step NNNNN (PETSc binary format) |
| `grains.opts` | PETSc options file with mesh/physics parameters |
| `grains.env` | (Legacy) Shell key=value file; auto-converted to `grains.opts` |
| `SSA_evo.dat` | (Optional) SSA time series — copied to output by run script |
| `metadata.json` | (Optional) Simulation metadata — copied to output by run script |

---

## Output files

| File | Description |
|------|-------------|
| `k_eff.csv` | Effective conductivity tensor; one row per time step |
| `t_vec.dat` | Correction temperature solution vector (PETSc binary) |
| `igaice.dat` | IGA descriptor for the homogenization mesh |
| `ice_data.dat` | Ice fraction at Gauss points (diagnostic; serial only) |

`k_eff.csv` format:
```
timestep,k00,k01,k10,k11
0,1.154999,0.000000,0.000000,0.063210
```

---

## Post-processing

```bash
# Plot k_eff vs time
python3 postprocess/plot_k_eff.py /path/to/output_dir

# Plot k_eff vs temperature
python3 postprocess/plot_k_eff_by_temp.py /path/to/output_dir

# Plot k_eff vs ice fraction
python3 postprocess/plot_k_eff_by_phi.py /path/to/output_dir

# Batch all plots for a set of output directories
./postprocess/run_all_single_plots.sh
```

---

## Parallel usage notes

- **Direct solver required.** Iterative solvers (CG, GMRES, GAMG) converge to
  numerically polluted solutions for periodic homogenization problems with a
  null space. Always use `KSPPREONLY + PCLU` (direct LU).
- **Null space.** The periodic system matrix is singular (constant modes). The
  code eliminates the null space by pinning component m at node 0
  (`MatZeroRowsColumns`), making the system non-singular for direct LU.
- **SuperLU_DIST.** Selected automatically for parallel runs. MUMPS is not
  available in this PETSc build.
- **Binary I/O.** `VecView` (binary write) is a collective operation — all ranks
  must call it. Do not guard it with `if (rank == 0)`.

---

## Validation

For a layered 50% ice microstructure, the exact analytic k_eff values are:

| Direction | Formula | Value (k_ice=2.29, k_air=0.02) |
|-----------|---------|-------------------------------|
| Parallel (x) | arithmetic mean | 1.155 W/m·K |
| Perpendicular (y) | harmonic mean | 0.039 W/m·K |

Run with `-init_mode layered -eps 1e-6` and a fine mesh to approach the harmonic mean.
