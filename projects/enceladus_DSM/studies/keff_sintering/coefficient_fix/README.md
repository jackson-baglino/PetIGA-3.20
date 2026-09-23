# What the k_eff replay actually runs

## Words used precisely here

I had been using "disk" for four different things. In this directory:

| term | means |
|---|---|
| **grain** | one ice particle in the packing, `R_ave = 50 µm` |
| **cylinder benchmark** | the single-ice-cylinder verification case in `studies/keff_sharp_limit/disk/` (a 2D circle, checked against Rayleigh's square-array formula). Nothing to do with the packing. |
| **stored / on the filesystem** | files on `$SCRATCH` or in `~/SimulationResults` |
| **snapshot** | one `sol_NNNNN.dat` — the full 3-field solution at one instant, 192 MB |
| **run** | one physical trajectory: one packing seed, one temperature, 0 → 30 days |
| **leg** | one SLURM job's worth of a run. A run killed at the 24 h walltime and resumed is **one run in two legs**, in two different directories. |
| **replay** | re-reading stored snapshots and recomputing `k_eff` on them. **No time integration. Nothing is simulated.** |

## The one-sentence version

**We are not running any simulations.** The four pilot trajectories were
integrated in September and are finished. The replay re-opens their stored
snapshots and re-computes one number per snapshot, using a corrected
conductivity law. The ice never moves.

## What one job does, step by step

Each job runs `enceladus_dsm` with `-keff_replay <dir>`, which dispatches
**before** `TSSolve` (`src/enceladus_main.c:1892-1896`). So:

1. Build the same IGA mesh the original run used (2829², p=2, C=1, periodic).
2. List `sol_*.dat` in `<dir>`, sort by step, take every Nth.
3. For each one:
   a. `IGAReadVec` the stored solution; take `φ` (dof 0). **`φ` is read, never
      advanced** — no residual, no Jacobian, no Newton, no timestep.
   b. Evaluate `K` at every quadrature point from `φ` and `∇φ`, under
      `-keff_interp tensor`.
   c. Solve **two** steady-state scalar Poisson problems (one per direction
      `m`) for the periodic corrector `t_m`:
      `−∇·(K ∇t_m) = ∇·(K e_m)`, CG + GAMG, rtol 1e-10.
   d. Average the flux: `k_eff[i][j] = (1/|Y|) ∫ K (∂t_i/∂x_j + δ_ij) dV`.
   e. Append one row to `k_eff_tensor.csv`.
4. Exit.

The cost is step 3c — two Poisson solves on 8 M scalar unknowns, about 100 s
at 241 ranks. Everything else is I/O.

**The original `k_eff.csv` is not touched.** Non-default laws write
`k_eff_<law>.csv` (`src/keff.c:377-381`), so the old and new numbers sit side
by side in the same directory.

## Why 8 directories and 8 jobs, for 4 runs

There are **four runs** — pilot seeds 1–4, all at φ=0.325, `L/R_ave`=40,
`R_ave`=50 µm, `alpha_c`=1e-3, T=−20 °C, 0 → 30 days.

Every one of them **hit the 24 h walltime and was resumed**, so every one is
stored as two legs in two different places. That is a property of how the runs
were produced, not of the replay:

- leg 1 was a batch job → `$SCRATCH/enceladus_DSM/batch_2026-09-16__13.16.11_pilot_keff/<geom>__<exp>/`
- leg 2 went through `submit_enceladus.sh`, which sets no `BATCH_OUT_DIR`, so
  it landed in the single-run tree → `$SCRATCH/enceladus_DSM/<geom>/<ts>_<exp>_resume_job<id>/`

`-keff_replay` scans **one** directory for snapshots, so a leg is the unit of
work. 4 runs × 2 legs = **8 directories = 8 jobs**:

| run | leg 1 (batch dir) | leg 2 (resume dir) | join | snapshots |
|---|---|---|---|---|
| seed 1 | 0 → 14.90 d | 14.90 → 30.01 d | 14.90 d | 499 + 166 |
| seed 2 | 0 → 16.22 d | 16.22 → 30.02 d | 16.22 d | 543 + 151 |
| seed 3 | 0 → 16.13 d | 16.13 → 30.01 d | 16.13 d | 540 + 154 |
| seed 4 | 0 → 14.42 d | 14.42 → 30.01 d | 14.42 d | 483 + 174 |

`compare_laws.py` stitches each run's two CSVs back together by time, dropping
leg-1 rows at or after the join (the resume recomputes them) — the same rule
`postprocess/merge_restart_legs.py` uses on the in-line CSV.

The downloaded copies in `~/SimulationResults` look like one directory per
seed because `merge_restart_legs.py` already merged them and
`thin_snapshots.py` cut 663 snapshots to 55. **`$SCRATCH` has the unmerged,
unthinned originals.** Pricing the replay off the thinned copies was a 12×
cost error.

## Why we do not replay `arith`

Each run's existing in-line `k_eff.csv` **is** the arithmetic result — the run
wrote it while integrating. There is nothing to buy. The only reason to replay
`arith` would be to check the replay path itself reproduces the in-line
numbers, and that is one leg of one seed, not eight jobs.

## The baseline is t = 1 day, not t = 0

The initial condition is an analytic clamped sum of `tanh` profiles
(`-ic_grain_union 0`), not an equilibrated phase field. The first hours of
every run are that field relaxing onto its equilibrium profile — that is not
sintering, and it is where the error is largest. Measured on pilot seed 1, the
log-log slope of SSA(t):

| t | 0 | 1 h | 2 h | 4 h | 8 h | 12 h | 24 h | 7 d |
|---|---|---|---|---|---|---|---|---|
| `d ln SSA / d ln t` | 0.000 | −0.014 | −0.029 | −0.053 | −0.075 | −0.080 | −0.076 | −0.077 |

The slope settles by ~8 h and is flat thereafter. So **anything read at t = 0
describes the initial condition, not the model.**

This is not a small correction. On the arith curves already stored:

| baseline | seed 1 | seed 2 | seed 3 | seed 4 | mean | sd |
|---|---|---|---|---|---|---|
| t = 0 | +17.4% | +19.8% | +22.3% | +18.8% | **+19.6%** | 2.1 |
| t = 0.5 d | +15.6% | +18.1% | +20.3% | +16.8% | +17.7% | 2.0 |
| **t = 1 d** | +13.9% | +16.4% | +18.3% | +15.0% | **+15.9%** | 1.9 |
| t = 2 d | +12.4% | +15.0% | +16.7% | +13.5% | +14.4% | 1.8 |

`compare_laws.py --baseline-days` defaults to **1.0**. The +19.6% figure
quoted earlier in this campaign was from t = 0 and should not be used.

There is a second reason beyond the IC transient, and it points the same way:
the tensor law's error estimate assumes the *equilibrium logistic* profile —
that is where the antisymmetry giving `Σ_t = 0` comes from. At t = 0 the
additive IC is not that profile anywhere two grains touch, which is exactly
where the contact conduction is decided.

## Running it

```bash
# on the cluster — always dry-run first
./scripts/HPC/submit_keff_replay.sh --dry-run \
    --roots $SCRATCH/enceladus_DSM/batch_2026-09-16__13.16.11_pilot_keff \
            $SCRATCH/enceladus_DSM/packing_2D_pilot_phi0.325_Rave50um_LR40_seed*_L2mm_eps1000nm_perxy_T-20
# expect: 8 directories, 8 jobs, ~233 samples, ~$19. Drop --dry-run to submit.

# after downloading
venv_enceladus/bin/python studies/keff_sintering/coefficient_fix/compare_laws.py <batch_dir>
```

The stride is derived from `--target-samples` (default 60 per run, legs summed
so both halves share a cadence), and `--max-cost` (default $100) refuses to
submit above a ceiling.

## What the result has to show

1. **The rise from the 1-day baseline**, `arith +15.9% → tensor ?`. This is
   the correction to the campaign headline and it rescales every threshold in
   `CAMPAIGN.md`.
2. **`ksp_its`** against the cylinder benchmark's ~2.5× arith. The tensor
   operator is locally anisotropic inside the band (`K_∥/K_⊥ ≈ 29` at φ=0.5).
   A larger rise on the packing is a real campaign cost signal.
3. **`phi_bar` identical to the arith run's**, sample for sample. The
   conductivity law does not touch the ice field; a difference means the wrong
   snapshots were read.

## Why this is needed at all

`studies/keff_sharp_limit/` already verifies the tensor law against closed
forms — a planar slab (`k_⊥` exact to 4e−7) and the cylinder benchmark
(error +1.212% → +0.006%, observed order ≈ 2.5). **Neither has a contact.**
The packing's grains are seated in exact tangency, and the arithmetic law's
diffuse band bridges that contact with conducting material. The cylinder
benchmark's ~+8% bias at the pilot's `eps/R_ave = 0.02` is the
isolated-interface part of the error with the contact effect excluded by
construction — a lower bound, not an estimate. This measures the real thing.

## Files

| | |
|---|---|
| `compare_laws.py` | discovery, leg stitching, table, cross-check, figure |
| `compare_laws.csv` | one row per (seed, law), with the baseline time recorded |
| `compare_laws.png` | (a) absolute `k_eff(t)`; (b) normalised to the baseline, with the excluded IC transient dotted |

Until the replay lands, both hold the arith baseline only, and the driver says
so on stdout.
