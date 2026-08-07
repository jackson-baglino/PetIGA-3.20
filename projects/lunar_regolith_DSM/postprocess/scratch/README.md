# postprocess/scratch — quarantine, not an archive

Plotting and analysis scripts that are no longer part of the standard
post-processing sweep. Kept only so one can be pulled back out if needed.

**This directory will be cleared out wholesale.**

To bring a file back: `git mv postprocess/scratch/<file> postprocess/<file>`,
then add it to the sweep in `run_postprocess.sh`.

## What is still live, and what it does

| File | Role |
|---|---|
| `plot_fields.py` | VTK conversion for 2D/3D runs → `vtkOut/`. Not a figure. |
| `plot_fields_highres.py` | Same, supersampled. Deliberately **not** in the sweep — it is slow and only wanted on demand. |
| `plot_mass.py` | Phase mass vs time → `plots/mass.png` + `plots/mass/` |
| `plot_porosity.py` | Porosity and SSA vs time, dual axis → `plots/porosity.png` |
| `plot_ssa.py` | Normalised SSA vs time → `plots/ssa.png` |
| `plot_timestep.py` | Adaptive time-step history → `plots/timestep.png` |
| `pplib.py` | Shared helpers (`load_ssa`, `auto_time_unit`). Imported by the above — infrastructure, not a plotter. |

## Capability that went away with this cut

Three scripts the runners used to call are now in here. If any of this is
wanted back, it is one `git mv` plus a line in the sweep:

- `plot_1d_profiles.py` — **all 1D post-processing.** `run_postprocess.sh` used
  to branch on `-dim 1` and call it for phase-field profiles and derived
  quantities.
- `plot_scalars.py` — generic `SSA_evo.dat` time-series. Overlaps what the
  repaired `plot_porosity.py` and `plot_ssa.py` now cover.
- `plot_energy.py` — phase-field free energy vs time, called per-run by
  `run_batch_postprocess.sh`.
