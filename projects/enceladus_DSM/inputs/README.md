# Sublimation phase-field simulation inputs

The opts files are split into three orthogonal categories. A run is composed
by passing one file from each category to PETSc, left to right:

```
mpiexec -np N ./enceladus_dsm \
    -options_file inputs/solver.opts \
    -options_file inputs/geometry/<geom>.opts \
    -options_file inputs/experiment/<exp>.opts
```

PETSc processes `-options_file` flags left to right; later entries override
earlier ones. The split lets you sweep one axis without touching the others
(e.g., the same geometry at different humidities, or the same experiment on
different meshes). In practice you never type this by hand — use
`./scripts/Studio/run_enceladus.sh <geom> <exp> [tag]`, which assembles the
three files and adds `-output_path`.

## Search `scratch/` before creating a new input file

Only `solver.opts` and this README live at the top of `inputs/`. Every geometry
file, experiment file, packing and mesh from earlier work sits in
`inputs/scratch/`, keeping its original sub-directory layout
(`inputs/scratch/geometry/molaro/…`, `inputs/scratch/experiment/base/…`).

**Before writing a new input file, look for one that already does the job:**

```bash
find inputs/scratch -name '*<something>*'
```

If it exists, move it back and use it — do not rewrite it under a new name:

```bash
git mv inputs/scratch/geometry/packing/<name>.opts inputs/geometry/packing/<name>.opts
```

The run scripts do this search for you: when a geometry or experiment name
misses, they check `inputs/scratch/` and print the exact `git mv` to run.

Keep the same sub-directory scheme on both sides — geometry under
`geometry/<family>/`, experiments under `experiment/<family>/`. `scratch/` is a
quarantine, not an archive: it will be cleared out wholesale, so anything worth
keeping needs to be pulled back into the live tree before then.

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
  comment (parsed by `run_enceladus.sh` for rank sizing — do not omit it)
- Geometry-specific solver overrides (e.g. `-ksp_type preonly -pc_type lu`
  for 1D, `-mob_sub` reductions for stiff pairs)

### `experiment/<name>.opts` — environmental conditions
The bits you change between successive runs of the same geometry:
- `t_final`
- `temp`, `humidity`, `grad_temp0`
- output cadence overrides

### Output cadence

Three mutually exclusive cadences; the first one set wins.

| flag | snapshots at | use when |
|---|---|---|
| `-t_out_log N` | N **log-spaced** times, `delt_t` → `t_final` | fitting a power law, or the mesh is big enough that the snapshot budget is tight |
| `-outp N` | every N-th accepted step | short runs; roughly log-spaced while the adaptive dt is still ramping |
| `-t_interv` / `-n_out` | time-uniform | watching a process that is uniform in t |

`-t_out_log` also overrides the safety valve that force-switches `-outp` to
1000 time-uniform snapshots once `t_final/(dtmax*outp)` exceeds 1000 — with an
explicit log budget there is nothing to cap. Start it later than the first
timestep with `-t_out_log_t0 <s>`.

Sizing it matters on large meshes: one `sol_*.dat` is `Nx*Ny*dof*8` bytes
(18.5 MB on the 1561x492 Molaro mesh, ~630 MB on the 9792x2671 sintering
mesh), and `plot_fields.py` roughly doubles that in `vtkOut/`.

`SSA_evo.dat` is written every step regardless of cadence, so the scalar time
series is never the thing you lose by outputting fewer fields.

## Valid `-ic_type` values

The two-phase solver (`src/<project>_main.c`) dispatches exactly these in 2D/3D;
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


## The `phi<X>_seed<N>_T-<M>` packing sweep is generated, not hand-written

`geometry/phi*_seed*_T-*.opts` and most of `experiment/` are **derived output**
of `preprocess/generate_study_opts.py`, which reads `inputs/packings/*/`
(`grains.dat` + `metadata.json` — those are the real source). Do not hand-edit
them; edit the generator and re-run it.

Only the temperatures the study drivers actually use are kept in the repo,
currently **T-20 and T-5**. The packing geometry is temperature-independent,
and because `eps` is sized from `R_feat` (the K&P Eq. 45 kinetic bound) rather
than from the temperature-dependent bound, `-eps` and the mesh come out
*identical* at every temperature — a `T-30` file differs from `T-20` only in
the `-eps_valid_temp` guard. So the extra temperatures cost 150 tracked files
and buy nothing until a driver needs them.

To bring a temperature back:

```bash
preprocess/generate_study_opts.py \
    --packings-dir inputs/packings \
    --temps -30                       # or: -5 -10 -15 -20 -25 -30 -35 -40
```

## Geometry naming scheme

Same idea as the experiment files, one family subdirectory each:

```
<family>/<family>_<dim>D[_<discriminator>]_L<size>_eps<eps>[_<special>...]
```

- **family** — `molaro`, `packing`, `rev`, `wedge`, `throat`, `regolithpore`,
  `porechannel`, `bumpyfloor`, `twograins`, `iceslab`, `singleice`,
  `multigrain`, `ripening`, `sediment`, `calib`, `tgrad`, `bump`, `sphere`,
  `sinter`. Also the subdirectory name.
- **dim** — `1D` or `2D`, from `-dim`.
- **discriminator** — only where a family needs one; `packing` carries its
  `phi<X>_seed<N>`.
- **L`<size>`** — domain extent from `-Lx`, in `um` below 1 mm and `mm` above.
- **eps`<eps>`** — the interface parameter from `-eps`, in um.
- **special** — anything set for that one geometry: `axisym`, `union`,
  `fixgeom`, `T-20pair`, `icecap`, `2band`, `90deg`, `narrowL`, `molaroscale`, …

**Vague resolution qualifiers are gone.** `_hires`, `_epsloose`, `_epsmid`,
`_epsstrict` and `_refined2x` all named a resolution without saying what it
was, and the same word meant different values in different families. They are
now the actual `eps`, so two files are directly comparable:

| was | is |
|---|---|
| `1D_ice_slab` / `1D_ice_slab_hires` | `iceslab_1D_L100um_eps0.71um` / `…_eps0.36um` |
| `2D_molaro_axisym_T-20pair_epsloose` | `molaro_2D_L386um_eps0.86um_axisym_T-20pair` |
| `2D_molaro_axisym_T-20pair_epsstrict` | `molaro_2D_L385um_eps0.35um_axisym_T-20pair` |

As with experiments you pass a **bare name** — the run scripts search the
family subdirectories. `family/name` also works.

### `geometry/meshes/` and `geometry/multigrain/`

igakit mesh `.dat` files live in `geometry/meshes/`, separate from the
hand-written `.opts` that point at them via `-geom_file`. The exception is
`geometry/multigrain/<name>/`, where each mesh keeps its own
`build_geometry.py` + `mesh.dat` + `mesh.opts` together — that layout is the
reproducible record of how the mesh was built, so it stays intact.

> Four of those `multigrain/*/mesh.dat` have never been committed
> (`compEps_304x122`, `random_bumpy_floor_v1`, `random_bumpy_floor_molaroscale`,
> `single_bump_v1`). Regenerate with the `build_geometry.py` next to each one
> before using its `mesh.opts`.

## Experiment naming scheme

Every file under `experiment/` follows one format, in one campaign subdirectory:

```
<campaign>/<campaign>_T<temp>_h<humidity>_<duration>[_G<grad>][_a<alpha_c>][_<special>...]
```

- **campaign** — what the run is for: `molaro`, `snow`, `tgrad`, `ripen`,
  `collapse`, `sinter`, `calib`, `smoke`, `dtramp`, or `base` for the generic
  ones. Also the subdirectory name.
- **T / h / duration** — the three parameters that vary in essentially every
  run. Duration is an exact integer in the largest sensible unit (`30min`,
  `6h`, `30d`, `365d`); when no unit divides evenly it falls back to compact
  seconds (`1.50e5s`).
- **G`<grad>`** — temperature-gradient magnitude in K/m, present only when
  `-grad_temp0` is nonzero.
- **a`<alpha_c>`** — the condensation coefficient, present only where a
  campaign actually varies it (`molaro`, and `base` to separate the file that
  sets kinetics from the one that does not). Where a campaign holds `alpha_c`
  fixed (`snow`, `tgrad`) it is documented in the file header instead of
  cluttering every name.
- **special** — anything set for that one run and not normally touched:
  `openBC`, `rhovHiLo`, `epsloose`, `nopenalty`, `dtmax8e2`, `cfl1e4`, …

You still pass a **bare name** to the run scripts — they search the campaign
subdirectories for you. `campaign/name` also works if you want to be explicit.

```bash
./scripts/Studio/run_enceladus.sh molaro_2D_L385um_eps0.46um_axisym molaro_T-20_h1.00_2h_a3e-2
```

### If a file has no `_a<alpha_c>` token, check before trusting its physics

`-beta_sub0` and `-d0_sub0` default to **-5 C-era values** (`9.9e5`, `1e-7`).
A file that sets neither runs those defaults whatever `-temp` says. That is
fine for smoke tests and timestep probes, and wrong for physics you intend to
believe. Those files now carry an explicit warning header naming the sibling
to use instead.

> **On the old `_arrh` suffix.** It never meant an Arrhenius law — no such
> code exists in the solver. Those files simply hard-coded the `beta_sub0` /
> `d0_sub0` that an Arrhenius *fit* produced at one temperature, so they are
> now named by the `alpha_c` they encode: `_arrh` was `alpha_c = 1.341e-2`
> at -20 C, i.e. `_a1.34e-2`.

## Example: same geometry, two humidities

```bash
# Two-grain neck formation at saturation
./scripts/Studio/run_enceladus.sh twograins_2D_L41um_eps0.49um base_T-20_h1.00_1d

# Same geometry, undersaturated
./scripts/Studio/run_enceladus.sh twograins_2D_L41um_eps0.49um base_T-20_h0.95_1d
```
