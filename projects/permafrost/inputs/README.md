# Permafrost simulation inputs

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
different meshes).

## What goes where

### `solver.opts` — numerical / model defaults
Things that almost never change between runs:
- DOFs, polynomial order, continuity, periodicity
- Output cadence
- Model flags (`flag_avenue`, `n_relax`, BC type, `flag_Tdep`)
- Equation time-scaling (`xi_v`, `xi_T`)
- Phase-field bounds (`phase_lo`, `phase_hi`)
- Penalty coefficients (`difvap_pen`, `k_pen`, `k_sed_pen`, `Lambda`)
- Adaptive time-stepping caps (`dtmax`, `NRmin`, `NRmax`, `max_rej`)
- SNES / KSP / PC settings

### `geometry/<name>.opts` — geometry & mesh
- IC type (`-ic_type`)
- Dimension, mesh size, domain extents (`dim`, `Nx`/`Ny`/`Nz`, `Lx`/`Ly`/`Lz`)
- Geometric parameters (`RCice`, `RCsed`, `grain_sep`, per-grain radii)
- Interface width (`eps`) and initial time step (`delt_t`) — these are tied
  to the spatial resolution
- Geometry-specific solver overrides (e.g., `-ksp_type preonly -pc_type lu`
  for 1D tests where direct LU is fastest; `-mob_sub` reductions for the
  stiff touching/separated grain pair configurations)

### `experiment/<name>.opts` — environmental conditions
The bits you change between successive runs of the same geometry:
- `t_final`, `t_sed_freeze`
- `temp`, `humidity`, `grad_temp0`

## Available files

### Geometries
- 1D: `1D_ice_slab`, `1D_single_ice`, `1D_ice_sed_pair`, `1D_touching_grains`,
  `1D_separated_grains` (each with a `_hires` 2× variant)
- 2D: `2D_ice_slab`, `2D_single_ice`, `2D_single_sed`, `2D_ice_sed_pair`,
  `2D_touching_grains`, `2D_separated_grains` (each with a `_hires` 2× variant
  except `single_sed`)

### Experiments
- `1day_T-20_h0.95.opts` — 1 day at -20°C, undersaturated air (h = 0.95)
- `1day_T-20_h1.00.opts` — 1 day at -20°C, saturated air (h = 1.00)

## Example: same geometry, two humidities

```bash
# Run the 2D ice-sed pair at saturation
./scripts/Studio/run_batch_tests.sh \
    --geom 2D_ice_sed_pair --exp 1day_T-20_h1.00

# Same geometry, undersaturated
./scripts/Studio/run_batch_tests.sh \
    --geom 2D_ice_sed_pair --exp 1day_T-20_h0.95
```
