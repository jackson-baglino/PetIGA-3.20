# How sintering changes the effective thermal conductivity of an ice grain pack

Plan for the study: finalize configurations, run the campaign, analyse it.

Written to be executable and to say *why* at each step. Every number quoted as
"measured" comes from `studies/packing_design/`, which holds the script, the
CSV and the figure for it.

---

## Definitions used throughout

Terms introduced here are defined here, once.

| term | definition |
|---|---|
| `R_ave` | mean grain radius of the packing |
| `L/R_ave` | domain size in grain radii. The scale-invariant size — a packing at `L/R_ave = 40` is the same microstructure whether grains are 2.5 µm or 250 µm |
| `eps` (ε) | phase-field interface **decay length**, not the visible band width |
| **band** | `9.2·eps`, the width over which φ runs 0.01→0.99. The length that decides whether a gap is a channel or a weld |
| `k_eff` | effective thermal conductivity tensor from the periodic cell problem: impose a unit mean temperature gradient in each direction, solve for the periodic corrector, average the resulting flux. `src/keff_cell.c` |
| **`D_eff`** | **effective vapour diffusivity.** The *same* cell problem as `k_eff`, with the coefficient `K(φ)` replaced by `φ_air = 1 − φ` — the factor multiplying `D_v` in the solver's vapour equation (`assembly.c:215`). So `D_eff/D_0` is the fraction of bulk vapour diffusivity the pore network actually delivers across the domain: 1 for pure open air, 0 for a fully blocked pore. It replaces a yes/no percolation test, whose answer here is always "blocked", with a graded number |
| **usable pore** | pore further than `band/2` from any ice, i.e. where φ actually reaches ~0. Narrower pore conducts like solid and blocks vapour whatever the drawing shows |
| `SSA` | specific surface area, interface length per unit area. Computed as `∫|∇φ|`, which by the co-area formula equals the interface length exactly, at any `eps` |
| `L*` | kinetic crossover length `beta_HK·D_v`. Features below `L*` are attachment-limited, above it vapour-diffusion-limited. `beta_HK ∝ 1/alpha_c`, so **lower `alpha_c` means more attachment-limited** |
| **bias** | `(value at the eps in use) / (value extrapolated to eps → 0) − 1`. On a fixed packing, anything that moves with `eps` is a diffuse-band artifact; the intercept is the material |

---

## Where we are

Settled, each by measurement rather than argument:

- **Periodicity.** Fully periodic (`-periodic 1` + `--periodic xy`). Geometry
  files now carry their own `-periodic` from the packing metadata, and a
  packing the solver cannot express is refused. Removes the 90°-contact-angle
  boundary layer.
- **The metrics were wrong and are fixed.** Coordination counted at 9.5 nm
  while the band is 414 nm (so Z = 3.26, not 2.03); pore connectivity used
  4-connectivity on both phases and its cluster count diverged with raster.
- **Both phases cannot percolate in 2D.** 0 of 15 packings across porosity
  0.30–0.50. Physics of the plane, not a generator defect.
- **The biases are measured**, at `eps` = 45 nm on `L/R_ave = 20`:
  `k_eff` **+33%**, `SSA` **−11%**, `D_eff/D_0` **≤ 1.4e-3** (an upper bound —
  not raster-converged). Bias is linear in `eps`, matching the surface-excess
  theory in `calonne_to_phasefield_equivalence.tex` §6.
- **The `k_eff` bias is not a small-sample artifact.** Flat at +22.8 / +22.1 /
  +21.8% across `L/R_ave` = 10 / 20 / 40. What the bigger domain buys is
  precision: seed scatter falls 12% → 4.3%.
- **The bridges are not "a few".** 20% of all throats; contacts go 318 → 508.

Two things this implies that the campaign must respect:

1. **`k_eff(t=0)` is not usable as a baseline.** The `eps → 0` limit is tangent
   discs with *point* contacts — no more physical than the welded version. The
   band is standing in for necks that have not grown. Trust `k_eff(t)` only
   once neck radius comfortably exceeds the band.
2. **Long-range vapour redistribution is suppressed ≥100×**, so this model
   cannot produce long-range Ostwald ripening. Local sintering is unaffected
   (a neck and its feeding surfaces share a pore, and at `alpha_c = 1e-3` we
   are attachment-limited, `L*` = 144 µm). Declare it; do not claim results
   about grain-size-*distribution* evolution.

### The cost/bias trade, quantified

`bias ∝ eps` (measured) and `DOF ∝ (L/eps)²`. So **halving the bias costs 4×
the compute.** At `L/R_ave = 40` and `R_feat/R_ave = 1/50` (the convention),
`comp_eps` gives `eps` = 25 nm, `Nx` = 5657, **96M DOF, ~1200 cores per run**.
Extrapolating the measured ladder, the `k_eff` bias at `eps` = 25 nm is **~15%**.

| `R_feat/R_ave` | `eps` | `Nx` | DOF | `k_eff` bias |
|---|---|---|---|---|
| 1/25 | 50 nm | 2829 | 24M | ~33% |
| **1/50** | **25 nm** | **5657** | **96M** | **~15%** |
| 1/100 | 12.5 nm | 11314 | 384M | ~8% |

**Recommendation: stay at 1/50.** A ~15% bias that is *stable across seeds and
domain sizes* is a systematic offset that can be reported and, because it falls
as SSA falls, shrinks over the run. 384M DOF per run is not affordable across a
campaign.

---

## Part A — Finalize the initial configurations

Three questions remain. A1 needs a pilot run; A2 and A3 are cheap.

### A1. How long a run, and how big a domain? (needs the pilot, Part B stage 0)

Coarsening grows `R_ave`, so `L/R_ave` **falls during the run**. Starting at 40
and coarsening 2× ends at 20 — below the REV requirement just established. So:

- Measure the coarsening factor `g = R_ave(t_final)/R_ave(0)` in the pilot.
- Require **`L/R_ave(0) ≥ 40·g`**.
- Re-check the REV criterion on the **final** snapshot, not the first.

### A2. Porosity range — bounded above by solid percolation

Measured: the solid stops percolating between φ = 0.40 and 0.45. Above that
`k_eff` is not a conductivity of a connected medium and the homogenization
assumption fails outright.

**Use φ ∈ {0.25, 0.30, 0.35, 0.40}.** Confirm solid percolation on every
generated packing (the gate already does this).

### A3. Seeds — how many for a manuscript claim

At `L/R_ave = 40` the seed scatter in `k_eff` is 4.3% (3 seeds, so a poor
estimate of σ). For a claimed trend to survive review, the between-condition
difference must exceed the within-condition scatter. **Budget 5 seeds per
condition**, and verify against the pilot: run 5 seeds at one condition, and if
the standard error of the mean is not ≲2%, raise it.

### A4. Regenerate the production set

```bash
for phi in 0.25 0.30 0.35 0.40; do
  for s in 1 2 3 4 5; do
    venv_enceladus/bin/python preprocess/generate_packing_gravity.py \
      --Lx <40*R_ave> --porosity $phi --mean-r <R_ave> --sigma-ln 0.5 \
      --periodic xy --seed $s --band-per-mean-r <9.2*eps/R_ave> \
      --out inputs/packings/phi${phi}_seed${s}
  done
done
venv_enceladus/bin/python preprocess/generate_study_opts.py \
  --packings-dir inputs/packings/periodic_xy --out-dir inputs \
  --R-feat <R_ave/50> --temps ...
```

Accept on: solid percolates both axes; `coordination_at_band` recorded;
`max_void_radius_per_mean_r` inside its gate; porosity within tolerance.

---

## Part B — The simulation campaign

### Stage 0 — Pilot (4 runs, HPC). Do not skip; it sizes everything else.

One porosity (0.325), one temperature (−20 °C), `alpha_c = 1e-3`, 4 seeds, run
to a duration that clearly saturates. Purpose is calibration, not science:

| what it sets | measured from |
|---|---|
| `t_final` for the campaign | when `k_eff(t)` and SSA(t) flatten |
| `L/R_ave(0)` via A1 | coarsening factor `g` |
| seeds needed via A3 | scatter across the 4 |
| when `k_eff` becomes trustworthy | neck radius vs band |
| wall-clock and core-hours per run | the run itself |

**Nothing in Stage 1–2 should be submitted before the pilot is read.**

### Stage 1 — Main campaign: what changes `k_eff`?

Vary the two things expected to matter most, with statistics:

**porosity {0.25, 0.30, 0.35, 0.40} × temperature {−40, −30, −20, −10 °C} ×
5 seeds = 80 runs.**

Temperature is the strongest lever on kinetics (via `rho_vs` and `alpha_c(T)`)
and is the axis a tiger-stripe argument needs. Porosity is the strongest lever
on the initial structure.

If 80 runs at ~1200 cores is too much, drop to 3 temperatures and 4 seeds
(48 runs) and say so — do not silently reduce seeds, since the seed count is
what makes the trend a result rather than an anecdote.

### Stage 2 — Parameter sensitivity: what *else* moves it? (12–16 runs)

One-at-a-time from the Stage 1 centre point, 3 seeds each:

- **`alpha_c`** {1e-4, 1e-3, 1e-2} — the one genuine free parameter, and it
  moves `L*` across the grain scale, so it may change the *mechanism* and not
  just the rate.
- **Grain size `R_ave`** {×0.5, ×1, ×2} **at fixed `L/R_ave`** — the packings
  are then geometrically identical and only the physics changes, which isolates
  the grain-size effect from a domain-size effect. Varying `R_ave` at fixed
  physical `L` confounds all three.
- **`eps`** ×2 on one condition — the direct check that the ~15% bias behaves
  as the ladder predicts in a *running* simulation, not just at t = 0.

### Submission

Per the standing convention I will prepare the scripts and hand over the
commands; the runs are yours to launch. More than one job goes through
`scripts/HPC/submit_batch.sh`, never chained single submits. Push before
submitting so the cluster checkout matches.

---

## Part C — Analysis

### C1. The headline: does `k_eff` track SSA, and why?

The existing observation is a tight SSA–`k_eff` correlation across porosity and
temperature. The analysis has to decide whether SSA is a **cause** or a
**proxy**, because only one of those is a result.

Regress `k_eff(t)` against, in order of increasing structure:

1. porosity alone (the null model — Calonne's density parameterization)
2. porosity + SSA
3. porosity + SSA + mean neck radius `r/R`
4. porosity + `coordination_at_band`
5. porosity + chord-length statistics

**If SSA adds nothing over porosity, the correlation is a coincidence of both
tracking density.** If it does add, report the partial correlation, not the raw
one. Do this pooled across all conditions *and* within each condition — a
relation that holds within a condition is much stronger evidence than one that
only appears when porosity varies.

### C2. Anisotropy

Report `k_yy/k_xx`, not just the isotropic average. Measured at 0.88–0.92 at
`L/R_ave = 40`, consistently across seeds: a real ~10% structural anisotropy
from gravity deposition. Track whether sintering **increases or erases** it —
either is a publishable observation, and it connects to the structural
anisotropy known to control snow `k_eff` (Calonne, Löwe).

### C3. Metrics to compute per snapshot

Beyond SSA, which already exists:

| metric | why | cost |
|---|---|---|
| `k_eff` full tensor | already in the solver | free |
| neck radius `r/R` | the classic sintering coordinate; also sets when `k_eff` becomes trustworthy | existing `neck_width.py` |
| **Euler characteristic** χ of the ice | connectivity once grains merge and coordination stops being definable; standard in snow microstructure | new, cheap |
| **correlation length** ξ from `S₂(r)` | the REV criterion made quantitative — checks the domain is *still* an REV at `t_final` | new, cheap |
| `D_eff` | how much the 2D pore fragmentation throttles vapour, as a function of time | `cell_solve.py` exists |
| chord-length distributions | what Calonne/Torquato `k_eff` models are built on — the sharpest test of C1 | new |

### C4. Statistics for the manuscript

- Report mean ± standard error across seeds for every condition; never a single
  realization without its scatter.
- Show one worked example trajectory *plus* the ensemble — the example for
  mechanism, the ensemble for the claim.
- State the systematic `eps` bias (~15% on `k_eff`, −11% on SSA) as a declared
  offset with its measurement, rather than hoping nobody asks.
- State the 2D limitations explicitly: no long-range ripening; both phases
  cannot percolate; anisotropy is deposition-induced.

---

## Immediate next actions

1. **Pilot first.** Generate 4 seeds at φ = 0.325, `L/R_ave` = 40, and submit
   Stage 0. Everything downstream depends on its numbers.
2. While it runs, build the C3 metrics (χ, ξ, chord lengths) and the analysis
   harness, so the main campaign is read the moment it lands.
3. Read the pilot, fix `t_final`, `L/R_ave(0)` and the seed count, then submit
   Stage 1.
