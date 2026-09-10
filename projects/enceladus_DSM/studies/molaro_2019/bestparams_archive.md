# What the archived `BestParams` runs actually did

Source: `~/SimulationResults/dry_snow_metamorphism/archived_results/BestParams/`
(NASAv2, 2024-11-26). Audited 2026-09-10, because those runs are the provenance
for the `M_0`×5 / `alph_sub`÷100 arm of the 2026-09 three-option campaign and
the claim needed checking before it was cited again.

## The parameter sets

There are **two** groups, not one, differing only in the two scale factors:

| | `mob_sub` | `alph_sub` | `humidity` |
|---|---|---|---|
| Group1 | `5 * eps/3/tau_sub` | `0.01 * lambda_sub/tau_sub` | 0.70 |
| Group2 | `5.5 * eps/3/tau_sub` | `0.0075 * lambda_sub/tau_sub` | 0.70 |

`diff` of the two `NASAv2.c` copies is exactly those two lines. Everything else
is shared, and both groups ran the same two temperatures (−20 °C and −5 °C):

```
eps = 9.096e-07     Lx x Ly = 2.424e-4 x 3.884e-4 m     Nx x Ny = 134 x 214
dim = 2 (PLANAR, not axisymmetric)   delt_t = 1e-4   t_final = 7200 s
humidity = 0.70     grad_temp0 = (0, 1e-4, 0) degC/m
dif_vap = 2.178e-5 (nominal)   d0_sub0 = 1.0e-9   beta_sub0 = 1.4e5
xi_v = 1.0e-5   xi_T = 1.0e-4   Lambd = 1.0   flag_Tdep = 1
grains: R = 101 and 73 um, centres 172.95 um apart (grainReadFile-2_Molaro.dat)
```

## Were the ×5 / ÷100 factors on the executed path? YES

**This corrects an earlier reading.** A first pass concluded they were not,
because the scaled assignments sit at `NASAv2.c:1972-73` guarded by nothing
while `flag_Tdep = 1` selects a per-quadrature-point path at `:603-604` that
carries no factors — and the log's missing `FIXED PARAMETERS` line seemed to
confirm it. That inference was wrong. Tracing the whole path:

1. `Residual`/`Jacobian` (`:292-297`, `:415-420`) branch on `flag_Tdep`:
   `1` → the arrays `user->mob[]`/`user->alph[]`, `0` → the scalars
   `user->mob_sub`/`user->alph_sub`.
2. The `flag_Tdep == 1` block lives inside **`Monitor`** (`:570`), not in the
   assembly, and its last statement is `user->flag_Tdep = 0` (`:618`).
3. PETSc calls the monitor at the initial time, before the first step. So the
   block runs **once**, fills the arrays, and flips the flag.
4. Every residual and Jacobian evaluation after that takes the `else` branch —
   **the scalars, with the ×5 and ÷100 factors applied.**

The log confirms the count directly: `b_min ... b_max` is printed once and only
once per run, and `FIXED PARAMETERS` never appears because line 1974 tests
`flag_Tdep == 0` at setup, when it is still 1. Absence of that line says
nothing about which branch the *solve* used.

**Consequence: the Libbrecht-derived per-point kinetics were computed and then
discarded.** The arrays are written, the flag flips, and nothing ever reads
them again.

## Was Libbrecht's data used to compute `beta_sub`? Computed, then thrown away

`Sigma0()` is Libbrecht's σ₀(T) table (3.0e-3 at 0 °C rising to 0.75 at
−100 °C; 3.5e-2 at −20 °C) and the code applies his law
`alpha = exp(-sigma_0 / sigma_surf)` with `sigma_surf = |rho_v - rho_vs|/rho_vs`
the *local* supersaturation, then `beta_0 = 1/(alpha * v_kin)`.

At `humidity = 0.70` the pore sits at `sigma_surf = 0.30`, which drives

| T | σ₀ | α from Libbrecht | β₀ |
|---|---|---|---|
| −20 °C | 3.5e-2 | **0.890** | 8.938e+03 s/m |
| −5 °C | 6.76e-3 | **0.978** | 2.026e+03 s/m |

The −20 °C run logs `b_min 8.94e+03 b_max 8.94e+03` — matching the computed
8.938e+03 to three figures, and confirming the reading. (`b_min == b_max`
because the IC is uniform in `rho_v`, so every quadrature point sees the same
supersaturation.)

**But that β never reached the solve.** The kinetics actually used descend from
the hardcoded `beta_sub0 = 1.4e5` at `:1968`, which corresponds to

| T | β_HK used | implied α_c |
|---|---|---|
| −20 °C | 1.293e-01 s/m | **0.057** |
| −5 °C | 5.044e-01 s/m | **0.014** |

So: Libbrecht's law is *present and evaluated*, gives α ≈ 0.89–0.98, and is
discarded; the run is governed by a fixed α_c of 0.057 (−20 °C), further
modified by ×5 on the mobility and ÷100 on the phase-change source.

## Constraint violations

**1. The interface is under-resolved.** `dx = 1.809e-6 m`, so `dx/eps = 1.99`
and only **4.6 elements** span the φ = 0.01–0.99 band. The current sizing
targets ~7.5–10 (see CLAUDE.md); the 2026-09 production mesh runs 13.2.

**2. The runs start below the neck resolution floor.** With
`sqrt(12·eps·R_ave) = 30.8 µm` as the smallest trustworthy neck *radius*, and a
reported initial neck width of 33.05 µm (radius 16.5 µm), the runs open at
**0.54× the floor** and only cross it near the end (38.2 µm radius, 1.24×).
This is the same defect that forced the 2026-08 campaign to abandon the tangent
start for a pre-necked `r = 14 µm` geometry.

**3. The initial neck is 1.75× the sharp-geometry value.** The centres are
172.95 µm apart against `R1 + R2 = 174 µm`, a 1.05 µm overlap, whose exact
lens intersection is a **18.8 µm** wide neck. The run reports 33.05 µm at
t = 0 — consistent with the additive-tanh IC bridging across a gap comparable
to its own diffuse width (9.19·eps = 8.4 µm), which is exactly why the current
geometries set `-ic_grain_union 1`.

**4. α_c sits outside the literature band on both readings.** Libbrecht's law
returns 0.89–0.98 here, ~9× above the band's 1e-1 ceiling
(`studies/alpha_c_sizing/`). The value actually used, 0.057 at −20 °C, is
inside the band — so this one is a violation only of the *intended* Libbrecht
parameterisation, not of the physics as run.

**5. Physically inconsistent settings for this problem**, each since changed:

- `Lambd = 1.0` — the triple-junction penalty is on, but this is a two-phase
  ice/air problem with no sediment and therefore no triple junction. Current
  `solver.opts` sets 0 and says why.
- `grad_temp0 = (0, 1e-4, 0)` — a transverse gradient in a nominally
  isothermal experiment.
- `humidity = 0.70` — `1 − h = 0.30`, about **100×** more undersaturated than
  the wall the 2026-09 campaign calibrated against Molaro's own measured grain
  recession (`1 − h = 2.9e-3`). Nothing in Molaro supports 0.70.
- `xi_v = 1e-5`, `xi_T = 1e-4` against the current 1e-3 and 1.0.
- **2D planar, not axisymmetric.** The grains are infinite cylinders, not
  spheres, which changes the curvature driving force at first order.

## What this means for citing these runs

The ×5 / ÷100 factors *were* the executed model, so `BestParams` is a real
tuning result and not a mis-recorded one. But it is not a like-for-like
precedent for the 2026-09 campaign: different dimensionality, an interface
under-resolved by ~3×, a start below the neck floor, an IC bridging artefact,
a wall 100× more undersaturated, and a different α_c. The 2026-09 arm-3 run
reproduced the factors in the current solver and was dropped on its own merits
(worst curve shape of four, and a neck fed by Allen–Cahn relaxation rather than
vapour — see `three_options/README.md` §6).
