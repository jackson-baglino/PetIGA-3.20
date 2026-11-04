# lunar_cryosuction — ice in regolith pore channels (paper 2)

Science runs on a frozen phase-field model (**PetIGA** on **PETSc 3.20**) of
**sublimation-driven ice metamorphism**: the coupled evolution of an ice phase
field, temperature and water-vapor density. The questions:

- how ice aggregates in pore channels bounded by regolith;
- how a cryosuction-like setting can be simulated, to understand ice harvesting
  on the lunar surface.

## Where this project came from

The thesis has one project per paper, as siblings under `PetIGA-3.20/projects/`:

| Paper | Project | Role |
|---|---|---|
| 1 | `enceladus_DSM` | dry snow on Enceladus |
| 2 | `lunar_cryosuction` (this one) | science on the frozen two-phase lunar model |
| 3 | `lunar_regolith_DSM` | literature review and a new model (most likely diffuse domain) |

This folder was copied from `lunar_regolith_DSM` on 2026-11-04, at git tag
`lunar_cryosuction/model-freeze-v1` (commit `147a49f0`). The copy came from
tracked files only (commit `6ebaab34`). Names and paths were changed in
`52e7ba16`. Apart from the main file's name (`lunar_main.c` →
`cryosuction_main.c`), `src/`, `include/` and `inputs/solver.opts` are
byte-identical to the tag. The model keeps developing in `lunar_regolith_DSM`;
none of that comes here.

**The model is frozen.** `FREEZE.md` gives what may change (ICs, geometry,
opt-in BC options, output, bug fixes) and what may not (physics, defaults),
and holds the ledger. `scripts/check_freeze.sh` checks both against git.

Records written before the split (`docs/`, study READMEs, `*.log`, the
`scratch/` folders, comments in `src/`) still use the old names
`lunar_regolith_DSM`, `lunar_main.c` and `run_lunar.sh`. Results made before
the split stay where those records say they are.

## Physics (two-phase)

Three degrees of freedom per node: ice φ, temperature T and vapor density ρ_v
(the air fraction is algebraic, φ_a = 1 − φ). The model has Allen–Cahn ice
evolution with a plain double well, a localized sublimation source, latent-heat
coupling in the temperature equation, and a vapor-transport equation with
Moure & Fu 2024 temporal scaling (ξ_v applied to diffusion *and* source
together; ξ_T = 1). The regolith is a wall boundary with a prescribed contact
angle (`-wall_faces`). Optional features: an axisymmetric r–z mode,
phase-field bounds enforcement, and an interface-CFL adaptive timestep.

## Layout

```
lunar_cryosuction/
├─ src/              # solver: cryosuction_main.c (main), assembly.c (residual/Jacobian),
│                    #   initial_conditions.c, material_properties.c, monitoring.c,
│                    #   snes_convergence.c, env_helper.c
├─ include/          # headers (NASA_types.h holds AppCtx + Field)
├─ makefile          # `make` (optimized) / `make debug`; builds ./lunar_cryosuction
├─ inputs/           # solver.opts, geometry/<family>/, experiment/<family>/, scratch/
├─ preprocess/       # comp_eps.py -> mesh_params/ (eps and mesh, tests in tests/), build_geometry_*.py
├─ postprocess/      # figure scripts, run_postprocess.sh, pplib.py
├─ scripts/
│  ├─ Studio/        # local runners (run_cryosuction.sh, run_batch_tests.sh)
│  ├─ HPC/           # SLURM: submit_cryosuction.sh, submit_batch.sh, run_cryosuction.sh
│  ├─ lib/           # alloc.sh (rank sizing), opts.sh (opts lookup)
│  ├─ check_freeze.sh     # the model is still the frozen model
│  └─ check_ic_types.sh   # every .opts -ic_type is in the solver dispatch
├─ studies/          # study folders with run lists, verification drivers and results
├─ docs/             # model and analysis notes (from before the split)
├─ memory/           # MEMORY.md (facts), USER.md (how the user works)
├─ FREEZE.md         # freeze rule and ledger
└─ LOG.md            # running session log
```

## Build

```bash
export PETSC_DIR=/path/to/petsc PETSC_ARCH=<arch> PETIGA_DIR=/path/to/petiga
make            # optimized; produces ./lunar_cryosuction
make debug      # -g3 -O0

python3.11 -m venv venv_cryosuction
venv_cryosuction/bin/pip install -r requirements.txt
```

On the HPC, load Python with `pyload` before creating or using the venv, and
load a gcc module before installing (igakit builds with f2py).

## Run

Never call the binary by hand. The run script puts the three opts files
together, sizes the rank count and stages a reproducible copy of the run:

```bash
scripts/check_freeze.sh
./scripts/Studio/run_cryosuction.sh <geometry> <experiment> [tag] [-- extra -flags]
```

Output goes to `~/SimulationResults/lunar_cryosuction/scratch/<geom>/<ts>_<exp>[_tag]/`.
On the HPC, `scripts/HPC/submit_cryosuction.sh` submits one run.
`scripts/HPC/submit_batch.sh --out-root /resnick/groups/rubyfu/jbaglino/simulation_outputs`
submits a batch to `<out-root>/lunar_cryosuction/`.

## Tests and checks

```bash
scripts/check_freeze.sh
scripts/check_ic_types.sh
venv_cryosuction/bin/ruff check .
scripts/lint_c.sh
venv_cryosuction/bin/python -m unittest discover -s preprocess/tests -t .
```

## Parameters

Always (re)compute the interface parameter ε and the mesh from the domain and
grain sizes, temperature and attachment coefficient with
`preprocess/comp_eps.py` (Kaempfer & Plapp 2009 bounds). ε is a decay-length
scale, not the visible width of the diffuse band; see
`docs/interface_width_conventions.md`.
