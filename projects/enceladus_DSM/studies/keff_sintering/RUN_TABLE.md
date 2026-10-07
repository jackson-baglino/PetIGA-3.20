# k_eff campaign — what was run (for the results section)

Compiled 2026-10-06 from the stage files, `inputs/packings/keff_LR40/*/metadata.json`
and the run folders in `LOC/keff_sintering_campaign/`; completed 2026-10-07 when
3f and the redo stage (3r) landed. History and data locations: `RECORD.md`.

## Main matrix: porosity × temperature

Number of independent packings run per cell. The matrix is complete: 125 runs.

| porosity φ | −40 °C | −30 °C | −20 °C | −10 °C | −5 °C | total |
|---|---|---|---|---|---|---|
| 0.275 | 5 | 5 | 5 | 5 | 5 | 25 |
| 0.325 | 5 | 5 | 5 | 5 | 5 | 25 |
| 0.375 | 5 | 5 | 5 | 5 | 5 | 25 |
| 0.425 | 5 | 5 | 5 | 5 | 5 | 25 |
| 0.475 | 5 | 5 | 5 | 5 | 5 | 25 |
| total | 25 | 25 | 25 | 25 | 25 | 125 |

- All 125 finished with data (2026-10-07). Nine of them are reruns from the redo
  stage: φ 0.425 seed 1902 at −40 °C (hung twice before; completed cleanly), the
  five seed-1 runs at −20 °C and seed 1701 at −5 °C (resampled every step), and
  two −40 °C runs whose tables were damaged by a double submission.
- The same packing is used at every temperature in its row: a row is 5
  packings, not 25. An ordering between porosities that repeats across
  temperatures is therefore one comparison, not five.

## Supporting runs (φ 0.325, −20 °C)

| study | what varies | runs |
|---|---|---|
| domain-size convergence | L/R = 20, 30, 56, 80 (L = 1, 1.5, 2.8, 4 mm), ungated packings | 8, 6, 4, 3 |
| | L/R = 40: the five production packings | (in the matrix) |
| gate comparison | L/R 40, packings built with the homogeneity gates | 5 |
| interface width | eps = 2 µm on seeds 1701–1703 | 3 |

154 run folders with data: 125 matrix, 21 convergence, 5 gated, 3 interface-width.

## The packings (five per porosity)

| φ target | φ achieved | grains per packing | realized mean radius | contacts per grain (z at the band) | ice spans both directions |
|---|---|---|---|---|---|
| 0.275 | 0.2751 | 297–345 (mean 318) | 49.0 µm | 3.74 | 5 of 5 |
| 0.325 | 0.3250 | 271–342 (mean 296) | 49.0 µm | 3.46 | 5 of 5 |
| 0.375 | 0.3749 | 268–289 (mean 281) | 48.2 µm | 3.33 | 5 of 5 |
| 0.425 | 0.4251 | 246–264 (mean 254) | 48.7 µm | 2.84 | 5 of 5 |
| 0.475 | 0.4750 | 217–240 (mean 225) | 49.7 µm | 2.54 | 4 of 5 |

Log-normal radii, target mean 50 µm, σ_ln = 0.5; largest radius capped at
100 µm, smallest present 6–12 µm. Drop-and-roll deposition, periodic in x and
y. Accepted on target porosity, y-seam contact density (≥ 0.76 of the
interior) and solid connectivity (off at φ 0.475).

## Fixed in every matrix run

| quantity | value |
|---|---|
| domain | 2 mm × 2 mm, periodic in x and y (L/R = 40) |
| interface parameter eps | 1 µm (visible 1–99 % band ≈ 9.2 µm) |
| mesh | 2829 × 2829 elements (p = 2, C¹); 24.0 M unknowns (ice phase, temperature, vapour) |
| run length | 30 days |
| condensation coefficient α_c | 1e-3, constant |
| humidity | 1.00 (saturated), isothermal, no imposed temperature gradient |
| conductivities | k_ice = 2.29, k_air = 0.02 W/m/K; tensor conductivity law |
| k_eff sampling | every step before 11 τ_sub, then every 0.1 % drop in SSA (≤ 20 τ_sub apart) |

## What changes with temperature

| T | τ_sub | 30 days in τ_sub | time-step cap | steps per run | k_eff samples |
|---|---|---|---|---|---|
| −40 °C | 60 390 s | 43 | 0.2 d (cap) | ~213 | ~178 |
| −30 °C | 20 835 s | 124 | 0.2 d (cap) | ~213 | ~178 |
| −20 °C | 7 822 s | 331 | 2·τ_sub = 0.18 d | ~238 | ~196 |
| −10 °C | 3 164 s | 819 | 2·τ_sub | ~485 | ~271 |
| −5 °C | 2 063 s | 1256 | 2·τ_sub | ~705 | ~314 |

Steps and samples from one representative packing (seed 1702); they vary by a
few between packings.

## Caveats to carry into the text

- Every matrix run samples k_eff at every step before 11 τ_sub (the six that did
  not were rerun in 3r). The gated seed-301 supporting run still has the coarse
  early sampling; it is used only for the dtmax and gate comparisons.
- The first day or so of every curve depends on the interface width
  (`eps_sensitivity/README.md`).
- Claims are restricted to φ ≤ 0.375 (`TODO.md`, "Manuscript argument").
