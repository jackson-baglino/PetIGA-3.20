# Sublimation phase-field simulation inputs

The opts files are split into three orthogonal categories. A run is composed
by passing one file from each category to PETSc, left to right:

```
mpiexec -np N ./permafrost \
    -options_file inputs/solver.opts \
    -options_file inputs/geometry/<geom>.opts \
    -options_file inputs/experiment/<exp>.opts
```

PETSc processes `-options_file` flags left to right; later entries override
earlier ones. The split lets you sweep one axis without touching the others
(e.g., the same geometry at different humidities, or the same experiment on
different meshes). In practice you never type this by hand — use
`./scripts/Studio/run_permafrost.sh <geom> <exp> [tag]`, which assembles the
three files and adds `-output_path`.

## What goes where

### `solver.opts` — numerical / model defaults
Things that almost never change between runs:
- DOFs (`-dof 3`: ice / temperature / vapor), polynomial order, continuity, periodicity
- Output cadence
- Model flags (BC type, temperature-dependent mobility)
- Equation time-scaling (`xi_v`, `xi_T`)
- Phase-field bounds (`phase_lo`, `phase_hi`)
- Adaptive time-stepping caps (`dtmax`, interface-CFL limiter) and SNES / KSP / PC settings

### `geometry/<name>.opts` — geometry & mesh
- IC type (`-ic_type`) — see the valid list below
- Dimension, mesh size, domain extents (`dim`, `Nx`/`Ny`/`Nz`, `Lx`/`Ly`/`Lz`)
- Grain radii / placement (`RCice`, per-grain radii)
- Interface width (`eps`) and initial time step (`delt_t`) — tied to spatial
  resolution; always (re)computed with `preprocess/comp_eps.py`
- For igakit meshes: `-geom_file <path>` plus a `# DOF_GRID: nx ny [nz]`
  comment (parsed by `run_permafrost.sh` for rank sizing — do not omit it)
- Geometry-specific solver overrides (e.g. `-ksp_type preonly -pc_type lu`
  for 1D, `-mob_sub` reductions for stiff pairs)

### `experiment/<name>.opts` — environmental conditions
The bits you change between successive runs of the same geometry:
- `t_final`
- `temp`, `humidity`, `grad_temp0`
- output cadence overrides

## Valid `-ic_type` values

The two-phase solver (`src/permafrost2.c`) dispatches exactly these in 2D/3D;
anything else aborts at startup. `scripts/check_ic_types.sh` enforces this
across every `.opts` file.

| ic_type | builder |
|---|---|
| `two_ice_grains_boundary` | `FormInitialTwoIceGrainsBoundary2D` |
| `ice_slab` | `FormInitialIceSlab2D` |
| `single_ice` | `FormInitialSingleIceGrain2D` (1D: `…1D`) |
| `multi_grains` | `FormInitialMultiGrains2D` — the multi-grain / igakit workhorse |

In 1D, any `-ic_type` other than `single_ice` falls through to
`FormInitialCondition1D` (centered slab / flat interface).

## Available geometries

- **1D:** `1D_ice_slab`, `1D_single_ice` (each with a `_hires` 2× variant)
- **2D:** `2D_single_ice`, `2D_ice_slab`, `2D_two_ice_grains_boundary`
  (+ `_refined2x` / `_refined4x`), plus the Molaro/axisym, pore-channel,
  bumpy-floor, ripening, and multi-grain families — see `inputs/geometry/`.

> **Note.** The sediment-era geometries (`*_ice_sed_pair`, `*_separated_grains`,
> `*_touching_grains`, `*_single_sed`) and the `inputs/tests/` suite were
> dropped in the 2026-06-13 two-phase fork and moved to `_trash/`. They request
> `-ic_type` values (`ice_sed_pair`, `enclosed`, `single_sed`, `ice_cap`,
> `capillary`, `contact_sed`, `slab_and_grains`) the solver no longer implements.

## Example: same geometry, two humidities

```bash
# Two-grain neck formation at saturation
./scripts/Studio/run_permafrost.sh 2D_two_ice_grains_boundary 1day_T-20_h1.00

# Same geometry, undersaturated
./scripts/Studio/run_permafrost.sh 2D_two_ice_grains_boundary 1day_T-20_h0.95
```
