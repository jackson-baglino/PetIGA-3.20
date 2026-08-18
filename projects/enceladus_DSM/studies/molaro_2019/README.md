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

```bash
# 0. Parameter envelopes (already run; regenerate if anything changes)
python preprocess/estimate_alpha_c_molaro.py
python preprocess/estimate_vapor_bc_molaro.py --L-over-a 3

# 1. Verify the measurement tool before trusting any number it produces
python studies/molaro_2019/verification/verify_grain_shrinkage.py

# 2. IC check -- 10 steps on HPC (the production meshes are too big locally)
./scripts/HPC/submit_enceladus.sh \
    molaro_2D_L450x225um_eps0.26um_axisym_T-20pair_tangent_dom2 \
    molaro_T-20_h1.0000_2h_a1e-1_dirichlet icheck --time=0-01:00:00 \
    -- -ts_max_steps 10 -outp 1 -t_out_log 0
# then, on the cluster: python postprocess/plot_fields.py --dir <run> --steps 0 1 10
# and rsync only those three .vts + pf.pvd back.

# 3. Domain convergence (3 arms, ~2 h each, ~1130 core-h)
./scripts/HPC/submit_batch.sh --tag molaro2019_domain \
    --tests-file studies/molaro_2019/batches/domain_check.txt

# 4. Humidity fit (5 arms, 15 h each, ~13500 core-h) -- only after step 3
./scripts/HPC/submit_batch.sh --tag molaro2019_humidity \
    --tests-file studies/molaro_2019/batches/humidity_fit.txt

# 5. Measure ON THE CLUSTER, download only summary.csv
bash postprocess/run_batch_measure.sh $SCRATCH/enceladus_DSM/batch_<...>
```

### What to check at the IC step

- the banner prints `-ic_grain_union 1` (the **default** is the additive form,
  which sums two tanh tails to a spurious bridge at exact tangency)
- φ = 0.5 is tangent at exactly one point
- the intermediate-φ band at the cusp is ~4.4ε ≈ 1.14 µm axially — that is the
  correct diffuse image of a cusp, not a defect
- by step 10 the O(ε) bridge has filled (`t_fill ≈ 2.6 s` at α_c = 0.1) with no
  spurious structure spreading along the contact plane

### What `summary.csv` answers

`t_star_s` is the clock shift — when the model neck first reaches their 32.81 µm
— and `neck_w_at_78min_um` is the model's neck one Molaro-experiment later.
Target: **64.78 µm**. `dR_large_pct` is the humidity fit target: **−2.93 %**.
