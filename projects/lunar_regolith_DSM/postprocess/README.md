# Post-processing Scripts

All scripts require Python 3 with `numpy`, `matplotlib`, `scipy`, `igakit`,
and `cmocean`.

```bash
pip install numpy matplotlib scipy igakit cmocean Pillow
```

`cmocean` and `Pillow` are optional — the scripts fall back gracefully if they
are not installed (alternative colormaps are used; GIF output is disabled).

---

## Automation

When using the Studio or HPC run scripts, 1D post-processing runs automatically
after the simulation completes.  For a run with `-dim 1` in the opts file the
following files are written to the output folder:

| Output file | Contents |
|-------------|----------|
| `phase_step_NNNNN.png` | Per-step phase field figure (one per output step) |
| `phase_animation.gif`  | Animated GIF of all phase field steps |
| `derived.png`           | Ice volume, interface position, slab width vs. time |
| `scalars.png`           | Ice volume, Δice, interface density, Δt from SSA_evo.dat |

---

## Scripts

### `plot_1d_profiles.py` — 1D field snapshots

Reads `sol_*.dat` and (optionally) `soil.dat` from a 1D run and produces:

- **Per-step phase PNGs** (default) — one figure per snapshot with three panels:
  φ_i (ice, blue), φ_s (sediment, brown), φ_a = 1−φ_i−φ_s (air, green).
  Saved as `phase_step_NNNNN.png` in `--out-dir`.

- **Thermal overlay** (`--thermal`) — all snapshots of T and ρ_v overlaid on
  two panels with a physical-time colorbar.  Uses `cmocean.thermal` for T and
  `cmocean.balance` (diverging red–blue) for ρ_v.

- **GIF animation** (`--gif`) — animated GIF of the phase fields over time.

- **Derived quantities** (`--derived`) — ice volume fraction, interface position,
  and slab width vs. physical time.

```bash
# Per-step phase PNGs (written to run directory)
python plot_1d_profiles.py --dir /path/to/run

# Also produce GIF and thermal overlay
python plot_1d_profiles.py --dir /path/to/run --gif --thermal

# Write outputs to a separate directory
python plot_1d_profiles.py --dir /path/to/run --out-dir /path/to/figs --gif

# Limit to 8 snapshots
python plot_1d_profiles.py --dir . --max-steps 8

# Derived scalar quantities
python plot_1d_profiles.py --dir . --derived --save derived.png
```

**Flags**

| Flag | Default | Description |
|------|---------|-------------|
| `--dir DIR` | `.` | Directory containing sol_*.dat files |
| `--iga FILE` | `igasol.dat` | IGA geometry file |
| `--max-steps N` | all | Limit number of snapshots loaded |
| `--out-dir DIR` | same as `--dir` | Where to write per-step PNGs |
| `--gif [PATH]` | off | Produce animated GIF (default name: `phase_animation.gif`) |
| `--thermal` | off | Produce thermal overlay figure |
| `--derived` | off | Plot derived scalar quantities instead of field profiles |
| `--save PATH` | auto | Save path for `--thermal` or `--derived` figure |

---

### `plot_scalars.py` — Scalar time series (1D / 2D / 3D)

Reads `SSA_evo.dat` and plots ice volume, change in ice volume, interface
density, and adaptive time-step size over time.

```bash
python plot_scalars.py --file /path/to/SSA_evo.dat --time-unit h --save scalars.png
```

---

### `compare_runs.py` — Multi-run comparison

Overlays scalar evolution or final 1D profiles from several output directories.
Useful for comparing different temperatures, humidities, mesh resolutions, etc.

```bash
# Compare scalars from three runs
python compare_runs.py run_T-20 run_T-10 run_T-5 \
    --labels "T=-20°C" "T=-10°C" "T=-5°C" \
    --time-unit h --save compare_temp.png

# Overlay final 1D ice profiles
python compare_runs.py run_T-20 run_T-10 --profiles --save compare_profiles.png

# Normalise by initial values
python compare_runs.py run_A run_B --normalise
```

---

### `plot_2d_snapshot.py` — 2D field visualization

Produces a 6-panel matplotlib figure (ice, sediment, air, temperature, vapor
density, supersaturation) evaluated on a regular grid from a 2D PetIGA run.
Optionally adds horizontal/vertical cross-section profiles.

```bash
# Plot step 10
python plot_2d_snapshot.py --dir /path/to/run --step 10

# Higher resolution grid, with cross-section cuts
python plot_2d_snapshot.py --dir . --step 50 --nx 300 --ny 300 --cuts --save snap.png
```

---

### `analyze_interface.py` — Quantitative interface metrics

Extracts time-series of physically meaningful scalar quantities and optionally
saves them to CSV.

**1D metrics** (reads `sol_*.dat`):
- Interface position(s) and count
- Interface width (φ_i = 0.1 → 0.9 transition distance)
- Mean ice volume fraction
- Peak and mean-air vapor supersaturation

**2D/3D metrics** (reads `SSA_evo.dat`):
- Ice volume and its change from initial
- Interface density Σ/ε and interface area Σ

```bash
# 1D analysis with CSV export
python analyze_interface.py --dir /path/to/1D/run --dim 1 \
    --save-csv metrics.csv --save-fig metrics.png

# 2D/3D analysis
python analyze_interface.py --dir /path/to/2D/run --dim 2 --eps 9.3295e-7
```

---

### Legacy scripts (already present)

| Script | Purpose |
|--------|---------|
| `plot_fields.py` | Old minimal VTK converter (superseded by `scripts/plot_fields.py`) |
| `plot_ssa.py`        | SSA time-series (hardcoded paths, use `plot_scalars.py` instead) |
| `plot_porosity.py`   | SSA + porosity dual-axis (hardcoded paths) |
| `plot_triple_well.py` | 3D surface plot of the triple-well free energy potential |
| `convert_to_stl.py`     | Convert VTK to STL |

---

---

## `pplib.py` — shared helpers

Import from `pplib` rather than re-deriving. Running any script in this
directory puts it on `sys.path`, so `from pplib import rho_vs` just works.

| Helper | Replaces |
|---|---|
| `load_ssa(path, min_cols=4)` | 5 copies with 5 different contracts |
| `rho_vs(T_C, rho_air=1.341)` | 3 copies, **all of them wrong** (see below) |
| `supersaturation(rhov, T_C)` | the `(rhov - rvs) / rvs` idiom |
| `read_vts(fn, want=None)` | 3 copies |
| `step_times(path)` | 3 copies |
| `step_of(fn)`, `auto_time_unit(t)`, `in_time_unit(t, unit)` | 2 copies each |
| `SSA, ICE, TIME, STEP, DT, AIR, RHOV, MASS` | hardcoded column indices |

`load_ssa` accepts a run directory or the file itself, pads legacy 4-column
files with a NaN `dt`, filters NaN on the first four columns only, and
deduplicates by step. It returns `None` rather than exiting, so callers decide
how to handle a missing file.

### ⚠ Supersaturation figures made before 2026-08-07 are wrong

The three previous `rho_vs` copies all computed

```python
3.25e-3 * np.exp(-6150.0 / T_K)
```

which is about **1e10 times too small** — for that Clausius-Clapeyron form the
prefactor should have been ~3e7, so the exponent looks like a typo. The error
does not cancel in `(rhov - rho_vs) / rho_vs`: a *saturated* run, which should
report 0 supersaturation, reported about **9e9** instead.

`pplib.rho_vs` now mirrors the solver's `RhoVS_I` (`src/material_properties.c`)
exactly — same K0..K5 coefficients, same `rho_air * 0.62 * Pvs / (Patm - Pvs)`.
Postprocess and solver no longer disagree. Regenerate any supersaturation
panel produced by `analyze_interface.py`, `plot_1d_profiles.py` or
`plot_2d_snapshot.py` before trusting it.

## Output File Format Reference

| File | Contents |
|------|----------|
| `igasol.dat`     | IGA geometry for the primary field (iga object) |
| `sol_NNNNN.dat`  | Solution vector at step NNNNN: (ice, T, ρ_v) per node |
| `igasoil.dat`    | IGA geometry for the sediment field |
| `soil.dat`       | Sediment phase field φ_s per node |
| `SSA_evo.dat`    | Per-step scalars: `Σ/ε   ∫φ_i   t[s]   step` |
