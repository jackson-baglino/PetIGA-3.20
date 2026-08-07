# Code Notes & To-Do List

This file summarizes the current state of the code, all cleanup and changes made,
and outstanding work items. Items marked **[HARD]** are significant engineering
efforts best deferred until the model is otherwise stable.

---

## What Was Done (Cleanup Pass, April 2026)

### Files Deleted / Moved
| Action | File | Reason |
|--------|------|--------|
| Deleted | `src/enceladus_dsm.c` | Entire file was `//`-commented dead code (old `main`) |
| Moved | `-output_path` → `scratch/` | Oddly-named leftover run summary in root |
| Removed | `.DS_Store` (×2) | macOS metadata, not source code |
| Created | `.gitignore` | Excludes `obj/`, `venv_*/`, `.DS_Store`, `outputs/`, `scratch/` |

### Dead Code Removed
| File | What was removed |
|------|-----------------|
| `src/material_properties.c` | Old commented-out `Fair()` function (replaced by Lagrange-multiplier version) |
| `src/initial_conditions.c` | ~50-line commented-out IC block from before refactoring |

### Bugs Fixed
| File | Bug | Fix |
|------|-----|-----|
| `include/NASA_types.h` | `Nx, Ny, Nz` declared as `PetscReal` | Changed to `PetscInt` |
| `src/assembly.c` | `fair_ice` declared but never assigned in `Jacobian()` | Added `Fair(user, ice, sed, NULL, &fair_ice)` call |
| `src/enceladus_main.c` | `T_BC[dim][2], LL[dim]` were VLAs — unsafe for dim=1 | Changed to fixed-size `[3][2]` and `[3]` arrays |
| `src/enceladus_main.c` | `axis1` setup ran unconditionally even for dim=1 | Wrapped in `if (dim >= 2)` |
| `src/enceladus_main.c` | `nmb` calculation wrong for dim=1 or dim=3 | Rewrote with explicit dim==1, dim==2, dim==3 branches |
| `src/enceladus_main.c` | Ny, Ly option parsing and print statements ignored dim==1 | Added `if (dim >= 2)` guards |

### Jacobian Cleanup
The large commented blocks in `Jacobian()` (temperature and vapor spatial terms)
were replaced with concise TODO comments explaining what needs to be implemented.
The analytical Jacobian is **not yet registered** — the code uses `IGAFormIJacobianFD`
(finite-difference Jacobian). See `src/assembly.c` for details.

Also added a TODO comment noting that the ice-Jacobian coefficient uses `Etai` alone
but the Residual uses `(Etai + Etaa)`. These must be made consistent before the
analytic Jacobian can be used.

### 1D Support Added
- New function `FormInitialCondition1D()` in `src/initial_conditions.c`
  - `flag_tIC == 0` → centered slab (ice in [0.35 Lx, 0.65 Lx])
  - `flag_tIC == 2` → flat interface (ice in [0, 0.5 Lx])
- Routing in `enceladus_main.c`: if `dim == 1`, calls `FormInitialCondition1D()`
- Test input files created: `inputs/tests/test_1D_IceSlab.opts` and
  `inputs/tests/test_1D_FlatInterface.opts`

---

## Outstanding Issues

### High Priority

**[1] Complete the analytic Jacobian** (`src/assembly.c:Jacobian`)  
The temperature and vapor Jacobian rows are currently stubs (time-derivative only).
The finite-difference Jacobian is used instead (`IGAFormIJacobianFD`). FD Jacobian:
- Works for correctness verification
- Roughly doubles the residual evaluations per Newton step
- Not scalable to large meshes

What needs to be done:
- Restore temperature Jacobian (`J[1][*]`) with `rho*cp` mass matrix, `thcond`
  diffusion, and latent-heat coupling terms
- Restore vapor Jacobian (`J[2][*]`) with `air_eff` mass matrix, `difvap`
  diffusion, temperature coupling, and phase-change source
- Fix coefficient inconsistency in ice Jacobian: currently `3*mob/(Etai*eps)`,
  should be `3*mob/((Etai+Etaa)*eps)` to match the Residual
- Add `Etaa` to the Jacobian's local variable list
- Re-enable with `IGASetFormIJacobian(iga, Jacobian, &user)` in `enceladus_main.c`

**[2] Validate the temperature-equation coupling in the Residual** (`src/assembly.c:Residual`)  
The old comment `// R_tem = N0[a] * tem_t;  // Does not solve temperature equation`
suggests the temperature residual went through a revision. Current active form is:
```c
R_tem = rho * cp * N0[a] * tem_t
      + xi_T * thcond * N1[a][l] * grad_tem[l]
      + xi_T * rho * lat_sub * N0[a] * air_t;
```
This looks correct for the energy equation but needs to be cross-checked against the
model derivation document.

**[3] Fix the `flag_Tdep` reset in `monitoring.c`** (`monitoring.c:70`)  
```c
user->flag_Tdep = 0;  // Disable after first call
```
This resets the flag after the first monitoring step, meaning T-dependent
Gibbs-Thomson parameters are only ever computed at `t = dt`. If the temperature
field evolves significantly, this is incorrect. Either:
- Remove the reset so parameters are updated every step (may be expensive), or
- Add a user-controllable update frequency

**[4] Grain initialization TODO** (`src/grain_initialization.c`)  
A comment marks an uncertain sediment component index:
```c
PetscInt sedComp = 0;  /* <<< TODO: set this to the correct component index */
```
Verify this is correct for the current DOF ordering (ice=0, temp=1, vapor=2).

---

### Medium Priority

**[5] IC selection is hard-coded in `enceladus_main.c`**  
The initial condition is selected by commenting/uncommenting lines. This is fragile.
Consider adding a `-ic_type` string option and dispatching via a lookup table.
Example ICs to support: `capillary`, `layered`, `enclosed`, `packed`, `flat`, `1D_slab`.

**[6] `monitoring.c` uses `getenv("folder")` for file output** — STILL OPEN, and
now in three places, since `keff.c` follows the same convention for `k_eff.csv`.
The output path is fetched via `getenv("folder")`, but `user->output_path` is
parsed from `-output_path` and read by *nothing*. These should be unified. The
k_eff work deliberately matched the existing (env-var) convention rather than
introducing a second one mid-merge; unifying is a separate change that touches
every writer at once, plus `scripts/*/run_enceladus.sh`, which currently exports
`folder` **and** passes `-output_path` to the same run.

**[7] Duplicate plotting scripts**  
`scripts/plot_fields.py` and `postprocess/plot_fields.py` appear to be the
same file (or close versions). Similarly for `plot_ssa.py` and `plot_porosity.py`.
Pick one canonical location (`scripts/`) and remove the other copies.

**[8] Virtual environments should not be in the repo**  
`venv_permafrost/` and `venv_pf311/` are committed but now excluded by `.gitignore`.
They were not removed from the working tree (they are large). To clean them from git
history: `git rm -r --cached venv_permafrost/ venv_pf311/`. This is safe since
`requirements.txt` documents the Python dependencies.

**[9] Hard-coded VLA grain arrays** (`initial_conditions.c`, `grain_initialization.c`)  
Several functions declare variable-length arrays like `PetscReal centX[3][numb_clust_sed]`.
VLAs are not standard C99/C11 for stack allocation when size is not a constant.
Consider replacing with `PetscMalloc` + `PetscFree`.

---

### Lower Priority / Future Work

**[10] 3D initial condition variants** [HARD]  
Only 2D-specific initial condition functions exist. 3D runs will use the generic
`FormInitialCondition3D` path which loads from a file. A native `LayeredPermafrost3D`
and `EnclosedPermafrost3D` would be useful for setup-free 3D testing.

**[11] Periodic BC inconsistency in 1D**  
The periodic flag correctly wraps axis0, but the index offset `k = user->p - 1` used
in the node-coordinate calculation has not been tested for 1D. Needs verification.

**[12] Automated regression test infrastructure**  
All tests in `TESTS.md` currently require manual inspection. A simple Python script
that runs a short test, parses the monitor output, and checks scalar quantities
against stored reference values would make CI possible.

**[13] `readFlag` description says "UPDATE IMPLEMENTATION!"** (`enceladus_main.c:52`)  
The flag for reading ice grains from file is set but loading logic may be incomplete.
Check `InitialIceGrains` in `grain_initialization.c` for the file-read path.

**[14] `Fwat()` is misnamed** (`material_properties.c`)  
The function is labeled "water phase" but the code uses "met" (metal/sediment) as the
inert solid. The function should be called `Fsed()` or `Fmet()` and its doc comment
updated accordingly.

**[15] Phase clamping in `Residual()` masks solver problems** (`assembly.c:54-57`)  
```c
if (ice < 0.0) ice = 0.0;
if (ice > 1.0) ice = 1.0;
```
These clamps prevent NaN propagation but break the consistency between the residual
and the Jacobian (the Jacobian does not account for these clamps). A better fix is
to improve time-step control so the phase field never goes significantly out of
[0,1], and then remove the clamps or replace them with a proper barrier potential.

---

## File Structure Reference

```
enceladus_dsm/
├── src/
│   ├── enceladus_main.c         Main program (IGA setup, time stepping, IC dispatch)
│   ├── assembly.c            Residual, Jacobian, Integration (weak forms)
│   ├── initial_conditions.c  All IC functions (2D, 3D, 1D)
│   ├── grain_initialization.c  Random grain placement (sed & ice)
│   ├── material_properties.c   Thermophysical property functions
│   ├── monitoring.c          TS monitor callbacks and output
│   ├── snes_convergence.c    Custom per-DOF SNES convergence test
│   └── env_helper.c          Environment variable parsing
├── include/
│   ├── NASA_types.h          AppCtx struct, Field typedef
│   ├── NASA_main.h           Master include
│   └── *.h                   Per-module headers
├── inputs/
│   └── tests/                Test .opts files (see TESTS.md)
├── scripts/                  Plotting and analysis scripts (canonical location)
├── postprocess/              Older copies of some scripts — consolidate into scripts/
├── preprocess/               Python grain generation tools
├── TESTS.md                  Test suite documentation
├── NOTES.md                  This file
├── README.md                 Project overview
├── makefile                  Build system (mpicc + PETSc + PetIGA)
└── requirements.txt          Python dependencies
```
