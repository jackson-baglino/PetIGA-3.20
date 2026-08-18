# Molaro et al. (2019) −20 °C grain-pair replication

Two ice grains (R = 101 / 72.5 µm) sintering at −20 °C in a Peltier cryostage:
neck WIDTH 32.81 → 64.78 µm over 78 min, grains shrinking a few percent.
Source data: `inputs/validation/molaro2019_fig11_T-20.csv`
(DOI 10.1029/2018JE005773, Fig. 11 / supplemental table).

This directory holds the parameter-choice envelopes that fix the run setup
*before* any HPC time is spent, plus the batch definitions and compiled results.

## Setup, and where each number comes from

| quantity | value | set by |
|---|---|---|
| `-alpha_c0` | 0.1 | `alpha_c_estimate.csv` — saturation + data floor |
| `-eps` | 2.5852e-07 | neck fillet, ρ ≥ 6ε at their first data point |
| `safety` | 0.999 | the fillet binds first, so K&P never does |
| `-humidity` | 0.9973 centre | `vapor_bc_estimate.csv` — large-grain recession |
| domain | L/a_eff = 3 (673 × 336 µm) | domain-convergence batch |
| `-t_final` | 15 h | time for *our* neck to cross their range, not their 78 min |

## Scripts

```bash
python preprocess/estimate_alpha_c_molaro.py      # -> alpha_c_estimate.{csv,png}
python preprocess/estimate_vapor_bc_molaro.py     # -> vapor_bc_estimate.{csv,png}
```

Both import their thermodynamics from `preprocess/comp_eps.py` so there is one
source of truth for ρ_vs(T), D_v(T), d0(T) and β_HK.

## What the envelopes concluded

**α_c is exhausted as a knob near 0.1.** The attachment-only limit puts a floor of
α_c ≥ 0.139 on what the model needs to match their rate, and going from α_c = 0.1
to the physical ceiling α_c = 1 buys only 1.25× because the process is already
transport-limited (Λ = β_HK·D_v/ρ = 0.29 at α_c = 0.1). Raising α_c also *tightens*
the K&P interface-width bound, since L* ∝ 1/α_c. The residual ~50 % gap is the
surface-diffusion share — see `docs/molaro_validation_synthesis.md` §4.

**Fit the humidity to the large grain only.** Its 78-min trend is 4.4σ and its
least-squares slope (−2.93 %) reproduces the −3 % the paper's caption quotes. The
caption's −4 % for the small grain is not supported by the same table (endpoint
−0.69 %, minimum −2.40 %, fit −0.89 %, all inside its own ±2 %).

**The differential is not a ripening signal.** The imposed wall undersaturation
(1 − h = 2.7e-3) outweighs the curvature difference between the two grains
(d0·(2/R_sm − 2/R_lg) = 7.9e-6) by **343×**, so both grains sublimate at nearly
the same rate regardless of their size difference.

**The vapour IC needs no change.** `L²/D_v ≤ 1.7e-2 s`, so the quasi-steady profile
establishes instantly and the uniform `hum0·rho_vs` pore value that
`SetNodeFields` writes is correct as-is.

## Workflow

See **[RUNBOOK.md](RUNBOOK.md)** — the campaign runs **one simulation at a
time**, because each result changes what the next run should be. Grain
recession is linear in the wall undersaturation, so a measured `dR_large_pct`
gives an exact Newton step to the next humidity; a blind 5-arm grid would spend
~5x the core-hours to learn the same thing.

```bash
# Parameter envelopes (already run; regenerate if anything changes)
python preprocess/estimate_alpha_c_molaro.py
python preprocess/estimate_vapor_bc_molaro.py --L-over-a 3

# Verify the measurement tool before trusting any number it produces
python studies/molaro_2019/verification/verify_grain_shrinkage.py

# Measure a finished run (works on one run dir or a batch parent)
bash postprocess/run_batch_measure.sh <run_dir>
```

## Cost, and where it went

The first pass through this setup budgeted ~2700 core-h per arm. Two of those
numbers were wrong, and fixing them cut it to ~220:

| | first pass | corrected | why |
|---|---|---|---|
| `-dtmax` | 8.6 s (`1.1*tau_sub`) | 1.07e3 = `0.8*dt_CFL` at their last neck | `tau_sub` is the Allen-Cahn *relaxation* time, not a stability limit for an implicit solver. `-dtCFL` is the real limiter: it caps max\|dφ\| per step at 0.2 from the measured field change, i.e. `dt <= 0.8*eps/v_n` — 2 s during the fast cusp fill, 1338 s at their last neck. dtmax is now a derived backstop at 0.8x that, **5.35x looser than the tested 200 s in `inputs/solver.opts`** — rung 1 confirms dtCFL actually binds before anything else is spent. |
| `-t_final` | 15 h | 6 h | Rescaled the repo's own calibrated 28 ks (at α_c = 1.341e-2) by the **total** resistance, not by β alone: `(β'+ρ/D_v)` ratio = 0.407, so t_end ≈ 3.16 h and t\* ≈ 0.42 h. Scaling by β alone gives 1.04 h and cuts the run short. |
| DOFs/rank | 80k | 100k | Top of PETSc's healthy band, still below the 108k/rank that demonstrably ran ~7 s/step. These runs are step-limited, so a smaller allocation costs almost no wall time. dom3 goes 260 → 208 ranks. |

Steps: the interface advances `0.8*eps = 0.21 µm` per step at the CFL cap and
the neck grows 0 → 32.4 µm, so expect a few hundred steps rather than ~6300.
**Confirm on rung 1 with `plots/timestep.png` before trusting the rest.**
