# CHANGES — effective_k_ice_homog

All notable changes made on branch `feature/refactor-modular-parallel-clean`.

---

## [Unreleased]

### Critical bug fixes — incorrect k_eff values

#### 1. Remove manual dV multiplication from integrand functions (`src/assembly.c`)

**What was wrong:**  
`AssembleStiffnessMatrix` and `ComputeKeffIntegrand` were multiplying the integrand
by `dV = weight * detJac` (or variants thereof). PetIGA's `IGAPointAddMat`,
`IGAPointAddVec`, and `IGAPointAddArray` already apply `JW = detJac * weight`
internally before accumulating into the global matrix/vector. The manual
multiplication caused double-counting and produced k_eff values that were a factor
of ~(domain_volume / 2) too small.

**Fix:** Remove all manual `dV` multiplications. User integrands must return bare
(un-weighted) values; PetIGA handles the integration weight.

**Verification:**  
Layered 50% ice field, 64×64 mesh, eps=1e-4:
- `k_eff[0][0]` = 1.155000 (arithmetic mean — exact ✓)
- `k_eff[1][1]` → 0.039 as eps→0 (harmonic mean ✓)
- Off-diagonals ≈ 1e-18 ✓

#### 2. Fix IGAPointFormGrad array stride (`src/assembly.c`)

**What was wrong:**  
`ComputeKeffIntegrand` used a `PetscReal grad_t[3][3]` array. PetIGA fills gradients
with row-stride = dim (number of DOFs), which is 2 in 2D. A fixed-size `[3][3]` array
has stride 3, corrupting gradient reads for 2D problems.

**Fix:** Use a flat `PetscReal grad_t[9] = {0.0}` buffer and index as
`grad_t[i * dim + j]`.

#### 3. Fix non-periodic layered ice field (`src/field_init.c`)

**What was wrong:**  
`ComputeLayeredIceField` used a single tanh at y = Ly/2:
```
phi(0) = 1,  phi(Ly) = 0  →  NOT PERIODIC
```
The cell problem requires strictly periodic microstructure. The discontinuity at
the periodic boundary introduced a non-physical interface that corrupted k_eff[1][1].

**Fix:** Two symmetric tanh interfaces at y = Ly/4 and y = 3Ly/4:
```
phi(0) = phi(Ly) ≈ 0  (air at periodic boundaries)
phi(Ly/2) = 1          (ice in the centre)
mean phi = 0.5         (50 % ice fraction)
```

#### 4. Switch from iterative solver to direct LU (`src/solver.c`)

**What was wrong:**  
The periodic homogenization system has a constant-mode null space (one mode per DOF
component). Iterative solvers (CG, GMRES, GAMG) always converge to a numerically
polluted solution for this problem class — k_eff comes out wrong regardless of
tolerance settings.

**Fix:** Use `KSPPREONLY + PCLU` (direct LU factorization). For parallel runs with
more than one MPI process, SuperLU_DIST is selected automatically
(`PCFactorSetMatSolverType(pc, MATSOLVERSUPERLU_DIST)`).

#### 5. Fix singular system with DOF pinning (`src/solver.c`)

**What was wrong:**  
With periodic BCs, the system matrix is singular (the solution is only unique up to
an additive constant per component). Bare `PCLU` on a singular matrix fails with NaN.
The previous MatNullSpace approach does not work with direct solvers.

**Fix:** Use `MatZeroRowsColumns` to pin component m at node 0 (global DOF = m).
This sets t_m(node 0) = 0, making the system non-singular and uniquely solvable by
direct LU. After solving, a zero-mean normalization step (parallel-safe via
`MPI_Allreduce`) centres the solution.

---

### Modular refactor (branch `feature/refactor-modular-parallel-clean`)

The original `effective_k_ice_homog.c` monolith (1294 lines) was split into modules:

| Module | File(s) |
|--------|---------|
| Main entry point | `src/effective_k_ice_homog.c` |
| Material properties | `src/material.c`, `include/material.h` |
| Ice field initialization | `src/field_init.c`, `include/field_init.h` |
| PETSc options / configuration | `src/env_config.c`, `include/env_config.h` |
| Assembly (stiffness, k_eff) | `src/assembly.c`, `include/assembly.h` |
| Solver (direct LU + null space) | `src/solver.c`, `include/solver.h` |
| I/O (CSV, binary, VTK) | `src/io_thermal.c`, `include/io_thermal.h` |
| IGA setup | `src/setup.c`, `include/setup_thermal.h` |
| Shared types / macros | `include/app_ctx.h` |

---

### Replaced env-var interface with PETSc options system (`src/env_config.c`)

**What was wrong:**  
`getenv()` calls are fragile: missing variables are silently `NULL`, there is no
`-help` documentation, and SLURM/HPC environments require exporting every variable
manually.

**Fix:** `ParseOptions(AppCtx *user)` uses `PetscOptionsBegin` / `PetscOptionsEnd`
with `PetscOptionsInt`, `PetscOptionsReal`, `PetscOptionsString`, `PetscOptionsBool`.
All parameters have hardcoded defaults and appear in `-help` output automatically.
An options file (`-options_file path/to/grains.opts`) replaces `grains.env`.

---

### Efficiency: reuse KSP and matrix across file iterations (`src/solver.c`)

**What was wrong:**  
`IGACreateMat`, `IGACreateVec`, and `IGACreateKSP` were called on every iteration
of the multi-file loop — expensive for large batches.

**Fix:** Introduced `CreateSolverObjects` (call once before loop) / `Solve` (call
each iteration) / `DestroySolverObjects` (call once after loop). The direct LU
factorization is rebuilt when `KSPSetOperators` is called (because the matrix
changes each iteration), but matrix/vector memory allocation and KSP object creation
happen only once.

---

### Scripts and pipeline

- `scripts/run_effective_k_ice_homog.sh` — rewritten to accept `-options_file`
  (grains.opts preferred over grains.env); passes `-init_mode` to the executable;
  fixed `POST_SCRIPT` default (`plot_k_eff.py`, not the removed `plot_vector_field.py`)
- `scripts/batch_effective_k_ice_homog.sh` — `PARENT_DIR` is now a CLI argument
  (with a fallback default for local development)
- `scripts/gen_opts.py` — new: converts `grains.env` (shell key=value) to PETSc
  `.opts` format automatically
- `scripts/run_hpc.slurm` — new: SLURM script for Caltech Resnick HPC cluster

---

### Input files

- `inputs/default.opts` — reference options file with all defaults documented
- `inputs/grains_template.opts` — per-simulation template (copy to simulation
  directory and fill in from `grains.env` or simulation metadata)

---

### Directory cleanup

Moved to `archive/` (not deleted):
- `src/effective_k_ice.c`, `src/effective_k_ice_full.c` — old abandoned mains
- `src/effective_k_ice_homog.py` — FEniCSx prototype
- `src/setup_thermal.c`, `src/utils.c`, `src/material_properties.c` — old stubs
- `scripts/batch_keff_eval.sh` — old utility script
- `scratch/` — scratch files
- `postprocess/old/` — legacy plot scripts

---

### `.gitignore` (new)

Tracks: build artifacts (`obj/`, binary), runtime outputs (`*.dat`, `*.bin`),
Python cache, generated plots, editor artifacts.
