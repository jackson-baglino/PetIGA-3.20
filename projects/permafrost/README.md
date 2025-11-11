# ❄️ Permafrost Model — PETSc + PetIGA

A modular, restartable, and reproducible simulation framework built on **PETSc 3.20** and **PetIGA** for modeling multiphase diffusion and phase-field evolution in porous ice and sediment systems.

---

## 🧭 Overview

This project implements a high-performance finite-element / isogeometric analysis model of ice–air–sediment interactions using **PetIGA**.  
It’s organized for clarity, reproducibility, and easy restart of simulations.  
Each run is fully parameterized by human-readable configuration files, with automatic metadata and checkpointing.

---

## 📂 Repository Layout

```
permafrost/
├─ include/                  # Public headers
│  ├─ appctx.h        # Simulation context & options parsing
│  ├─ geometry.h       # IGA creation, BCs
│  ├─ assembly.h       # Residual & Jacobian kernels
│  ├─ io.h             # Save/load checkpoints & metadata
│  ├─ restart.h        # L² projection between spaces
│  ├─ monitors.h       # Progress/output monitors
│  └─ utils.h          # Small reusable helpers
│
├─ src/                      # Implementation files
│  ├─ main.c          # Entry point, orchestration
│  ├─ geometry.c      # Builds IGA, sets BCs
│  ├─ assembly.c      # PDE residuals & Jacobians
│  ├─ io.c            # File I/O, metadata
│  ├─ restart.c       # Projection-based restart
│  ├─ monitors.c      # Logging & CSV outputs
│  └─ utils.c         # Misc helpers
│
├─ configs/                  # Default solver/model settings
│  ├─ base.opts
│  └─ solver_strict.opts
│
├─ inputs/                   # Immutable geometry/input bundles
│  ├─ snowpack_A/
│  │  ├─ igapermafrost.dat
│  │  ├─ pf_init_t0.dat
│  │  └─ input.json
│  └─ snowpack_B/ ...
│
├─ runs/                     # One folder per experiment
│  └─ 2025-11-05_temp-240K/
│     ├─ overrides.opts
│     └─ README.md
│
├─ examples/                 # Minimal run scripts
│  ├─ 2d_steady.sh
│  ├─ 2d_restart.sh
│  └─ 2d_refine_project.sh
│
├─ test/                     # Tiny regression/unit tests
│  ├─ test_load_direct.c
│  ├─ test_project_refine.c
│  ├─ test_pf_scatter.c
│  └─ test_incompatible_vec.c
│
├─ scripts/                  # Helpers and automation
│  ├─ run_case.sh
│  ├─ sweep_temps.csv
│  └─ sweep_temps.py
│
└─ docs/                     # Human-readable documentation
   ├─ design.md
   ├─ restart.md
   └─ options.md
```

---

## ⚙️ Building

### Using Makefile
```bash
export PETSC_DIR=/path/to/petsc
export PETSC_ARCH=arch-darwin-c-opt
export PETIGA_DIR=/path/to/petiga

make clean && make
```

### Compiler flags
- **Warnings:** `-Wall -Wextra -Wpedantic -Wshadow -Wconversion`
- **Debug:** `-g -O0 -fsanitize=address,undefined`
- **Release:** `-O3 -DNDEBUG`

### Example targets
```bash
make             # build
make tests       # build and run test suite
make clean       # clean artifacts
```

---

## ▶️ Running Simulations

### Basic steady-state example
```bash
mpirun -n 4 ./permafrost   -options_file configs/base.opts   -options_file inputs/snowpack_A/geometry.opts   -options_file runs/2025-11-05_temp-240K/overrides.opts   -output_dir outputs/2025-11-05_temp-240K
```

### Restart from previous state
```bash
mpirun -n 4 ./permafrost   -options_file configs/base.opts   -initial_cond outputs/2025-11-05_temp-240K/U_t1.234000e+02.dat
```

### Project to new resolution
```bash
mpirun -n 4 ./permafrost   -options_file configs/base.opts   -project_from outputs/2025-11-05_temp-240K
```

---

## 🧩 Source Code Responsibilities

| Module | Purpose | Key Functions |
|---------|----------|---------------|
| **main.c** | Driver; parse options, orchestrate run, choose initialization path, handle restarts. | `AppCtxLoadFromOptions`, `TryLoadVec`, `ProjectOldToNew`, `SaveCheckpoint`, `SNESSolve` |
| **geometry.c/h** | Define IGA spaces and BCs. | `IGACreatePrimary`, `IGACreateSoil`, `IGAConfigureBCs` |
| **assembly.c/h** | Physics; residual & Jacobian evaluation. | `FormFunction`, `FormJacobian`, `WireSNES`, `WireTS` |
| **io.c/h** | Read/write data, checkpoints, metadata. | `SaveCheckpoint`, `TryLoadVec`, `LoadLatestCheckpoint`, `WriteRunJson` |
| **restart.c/h** | Projection restart between spaces. | `ProjectOldToNew` |
| **monitors.c/h** | Progress monitoring, CSV logging. | `RegisterMonitorsSNES`, `MonitorWriteCSV` |
| **utils.c/h** | Small generic helpers. | `REQ`, `MakeDirIfNotExist`, `JoinArgv`, etc. |
| **appctx.h** | Unified runtime configuration structure. | `AppCtxLoadFromOptions`, `AppCtxEchoResolved` |

---

## 🧠 Design Philosophy

### Separation of concerns
- **Input bundle:** fixed geometry/initial state (immutable data).
- **Configuration:** stable defaults (checked into git).
- **Run overrides:** per-experiment tunables.

### Options precedence
PETSc allows multiple options files:
```
-options_file configs/base.opts -options_file inputs/snowpack_A/geometry.opts -options_file runs/.../overrides.opts
```
Later files override earlier ones — so every parameter has a single source of truth.

### Parameter snapshot
At startup, the resolved configuration is echoed and saved as `run.json` in the output directory for full reproducibility.

---

## 🔁 Restart & Checkpoint Workflow

1. **Checkpoint files**
   - `igapermafrost.dat` — IGA definition.
   - `U_tXXXX.dat` — binary PETSc vector.
   - `run.json` — metadata.

2. **Restart**
   - Direct restart: same IGA space → `VecLoad()`.
   - Projection restart: different space → `ProjectOldToNew()` (L² projection).

3. **Validation**
   - Checks DOF, vector size, and basic metadata before accepting a restart.

---

## 🧪 Testing Suite

| Test | Purpose |
|------|----------|
| `test_load_direct` | Save + reload same-space vector; check numerical equality. |
| `test_project_refine` | Coarse → fine L² projection; verify error convergence. |
| `test_pf_scatter` | Verify PF stride scatter initialization. |
| `test_incompatible_vec` | Deliberate mismatch → expect clean error message. |

---

## 🧰 Utility Scripts

- **`scripts/run_case.sh`** – wraps a standard run; combines base + input + override options.
- **`scripts/sweep_temps.py`** – reads `sweep_temps.csv` and spawns parameter sweeps.
- **`examples/*.sh`** – ready-to-run minimal demos.

---

## 🧑‍💻 Coding Style

- Two-line function header comments:
  ```c
  // Function: FormFunction
  // Assemble residual F(U)=0 at quadrature points.
  ```
- Prefer short, modular functions (≤ 80 lines).
- PETSc naming & capitalization conventions.
- Warnings enabled (`-Wall -Wextra -Wpedantic`).
- No hidden globals; everything flows through `AppCtx`.

---

## ✅ Definition of Done (for Restart/IO)

- `SaveCheckpoint()` produces:
  - `igapermafrost.dat`, `U_t*.dat`, and `run.json`
- `-initial_cond` loads compatible space.
- `-project_from` projects between different spaces.
- All tests and examples run successfully.
- Metadata snapshot ensures complete reproducibility.

---

## 🧩 Future Extensions

- Parallel HDF5 I/O backend.
- Automated mesh refinement tests.
- YAML-based configuration parser (optional).
- Web dashboard for sweep visualization.

---

## 👩‍🔬 Citation / Credits

Developed by the **Caltech Cryo-Physics Group**  
Built with [PETSc](https://petsc.org/release/) and [PetIGA](https://github.com/dalcinl/PetIGA).
