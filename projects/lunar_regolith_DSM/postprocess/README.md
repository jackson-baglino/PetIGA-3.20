# Post-processing

```bash
pip install numpy matplotlib scipy igakit
```

`igakit` is needed only by `plot_fields*.py` (VTK conversion); the rest run on
numpy + matplotlib + scipy. Everything renders through the `Agg` backend, so
these work headless on HPC.

---

## Running it

The whole sweep is one script, and it is what the run drivers call — so a local
run, an HPC run and a re-run of a downloaded folder all produce the same thing:

```bash
bash postprocess/run_postprocess.sh <run_folder>      # one run
bash run_batch_postprocess.sh <batch_folder>          # every run in a batch
```

`run_postprocess.sh` and a copy of this directory are staged into every run
folder, so a downloaded result can be re-processed with no source tree present.

## Where output goes

```
<run>/
  plots/            every figure
    porosity.png
    ssa.png
    mass.png
    mass/           per-phase mass breakdown
    timestep.png
  vtkOut/           plot_fields.py output
```

Figures go in `plots/`, and only figures. `plot_fields.py` writes to `vtkOut/`
because it is a format conversion for ParaView, not a figure.

---

## The scripts

| Script | Reads | Writes |
|---|---|---|
| `plot_porosity.py` | `SSA_evo.dat` | porosity + interface density vs time, shared axes |
| `plot_ssa.py` | `SSA_evo.dat`, `-grains_file` | ice-air surface area vs time |
| `plot_mass.py` | `igasol.dat`, `sol_*.dat` | phase mass vs time, total and per phase |
| `plot_timestep.py` | `outp.txt` | adaptive time-step history |
| `plot_fields.py` | `igasol.dat`, `sol_*.dat` | `.vts` snapshots for ParaView |
| `plot_fields_highres.py` | same | same, supersampled |
| `pplib.py` | — | shared helpers, not a plotter |

`plot_fields_highres.py` is deliberately **not** in the sweep: it is slow and
only wanted on demand.

```bash
python3 postprocess/plot_fields_highres.py --dir <run> --n-per-elem 4
```

### Two things worth knowing before reading the numbers

**Porosity is exact, not inferred from the domain size.** `Integration()` in
`src/assembly.c` accumulates `phi` and `1-phi` over the same quadrature, so
`porosity = tot_air / (tot_ice + tot_air)` needs no `-Lx`/`-Ly` and stays
correct under `-axisym`, where both integrals carry the same `2*pi*r` weight
and the ratio cancels it.

**The SSA column is not an area.** `monitoring.c` writes
`int phi^2 (1-phi)^2 dV / eps`. For the logistic equilibrium profile that
integral is `eps/6` per unit interface, so `plot_ssa.py` reports
`6 * column_0` as a length (2D) or area (3D). That factor assumes the
interface is at its equilibrium profile and resolved by the mesh — a smeared
interface reads high, so early-time values from an IC written at a different
`eps` are not trustworthy.

### `pplib.py`

Import from here rather than re-deriving; every one of these previously existed
as two to five copy-pasted variants that had drifted apart.

- `load_ssa`, and the column constants `SSA ICE TIME STEP DT AIR RHOV MASS`
- `read_opts`, `opt_float`, `grain_radii` — recover a run's settings and its
  packing from the `.opts` staged in the run folder
- `rho_vs`, `supersaturation` — mirror the solver's `RhoVS_I` exactly
- `read_vts`, `step_of`, `step_times`
- `auto_time_unit`, `in_time_unit`

---

## `scratch/`

Retired scripts, kept only so one can be pulled back out. **Pending deletion.**
See `scratch/README.md` — it names the capability that went away with the cut
(1D profiles, energy, generic scalars) and how to restore it.
