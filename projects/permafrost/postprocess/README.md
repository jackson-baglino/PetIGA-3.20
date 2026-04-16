# Post-processing Scripts

All scripts require Python 3 with `numpy`, `matplotlib`, `scipy`, and `igakit`.

```bash
pip install numpy matplotlib scipy igakit
```

---

## Scripts

### `plot1D_profiles.py` — 1D field snapshots

Reads `sol_*.dat` files from a 1D run and plots φ_i, T, ρ_v profiles across
multiple time snapshots overlaid in one figure.  Pass `--derived` for a second
figure showing ice volume fraction, interface position, and slab width over time.

```bash
# Field profiles (all snapshots)
python plot1D_profiles.py --dir /path/to/run

# Limit to 8 snapshots, save figure
python plot1D_profiles.py --dir . --max-steps 8 --save profiles.png

# Derived scalar quantities (interface tracking)
python plot1D_profiles.py --dir . --derived --save derived.png
```

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

### `plot2D_snapshot.py` — 2D field visualization

Produces a 6-panel matplotlib figure (ice, sediment, air, temperature, vapor
density, supersaturation) evaluated on a regular grid from a 2D PetIGA run.
Optionally adds horizontal/vertical cross-section profiles.

```bash
# Plot step 10
python plot2D_snapshot.py --dir /path/to/run --step 10

# Higher resolution grid, with cross-section cuts
python plot2D_snapshot.py --dir . --step 50 --nx 300 --ny 300 --cuts --save snap.png
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
| `plotpermafrost.py` | Old minimal VTK converter (superseded by `scripts/plotpermafrost.py`) |
| `plotSSA.py`        | SSA time-series (hardcoded paths, use `plot_scalars.py` instead) |
| `plotPorosity.py`   | SSA + porosity dual-axis (hardcoded paths) |
| `plotTripleWell.py` | 3D surface plot of the triple-well free energy potential |
| `convet2stl.py`     | Convert VTK to STL |

---

## Output File Format Reference

| File | Contents |
|------|----------|
| `igasol.dat`     | IGA geometry for the primary field (iga object) |
| `sol_NNNNN.dat`  | Solution vector at step NNNNN: (ice, T, ρ_v) per node |
| `igasoil.dat`    | IGA geometry for the sediment field |
| `soil.dat`       | Sediment phase field φ_s per node |
| `SSA_evo.dat`    | Per-step scalars: `Σ/ε   ∫φ_i   t[s]   step` |
