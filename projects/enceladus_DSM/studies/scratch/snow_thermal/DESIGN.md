# Study design: microstructure, eps, and the temperature sweep

Decisions made while building the pipeline (2026-07-27), with the measurements
behind them. `README.md` is the overview; `PLAN.md` is the execution plan; this
file records *why* the numbers are what they are.

---

## 1. A tangent 2D packing cannot transport vapour

Jamming means the contact network spans the domain in both directions, and in
2D a spanning contact network cuts the complement into isolated cells. Measured
on a 2 mm, 433-grain, porosity-0.30 packing at exact tangency: the pore space
fell into ~1000 clusters with the largest holding **11%** of the void. Eroding
or dilating the solid did not change it, so it is topology, not rasterization.

This is the tension between "grains must touch to sinter" and "pores must
connect to diffuse", and in 2D it is an obstruction rather than a knob. It does
not arise in 3D, where both phases percolate simultaneously.

**Resolution:** `--contact-gap` in `preprocess/generate_packing.py` jams against
an inflated radius `r + gap/2` while reporting the true `r`, leaving a uniform
surface-to-surface gap at every contact. Porosity is computed from the true
radii and is unaffected.

| gap | pore percolation (largest cluster, phi 0.25-0.40) | solid |
|---|---|---|
| 0 um | never (10-50%) | percolates |
| 2 um | never (14-66%) | already broken up |
| 4 um | **always (100%, one cluster)** | broken up |
| 6 um | always (100%) | broken up |

Pore and solid percolation are mutually exclusive throughout. That is the
correct *initial* state for this study rather than a defect: an unsintered
packing has no solid conduction path, so k_eff should start near the air value
and **rise** as vapour deposition builds necks. That rise is the measurement.

**Chosen: gap = 4 um.**

## 2. Porosity floor

Jamming against `r + gap/2` costs solid fraction, so the target must satisfy

    (1 - phi) * (1 + gap/(2r))^2  <~  0.84     (2D random close packing)

| gap | penalty | phi_min at r = 45 um |
|---|---|---|
| 4 um | 1.091 | **0.230** |
| 8 um | 1.186 | 0.292 |
| 12 um | 1.284 | 0.346 |

Porosity 0.20 was confirmed unreachable: it jams at s = 0.975 and overshoots to
0.24 with overlaps remaining. Hence the **0.25 - 0.40** range. Going lower needs
larger grains (the penalty is `(1 + gap/2r)^2`), at the cost of fewer grains per
REV.

**Chosen: phi = 0.250, 0.2875, 0.325, 0.3625, 0.400** (evenly spaced), x 5
independent seeds. Every (phi, seed) cell is generated from scratch — the
previous study grew higher porosities by adding grains to a shared base, which
left configurations nearly identical and confounded the porosity trend with a
single microstructure.

## 3. eps must not be sized with a fixed v_n

K&P Eq. 45 bounds eps by `d0 / (beta_sub * v_n)`. `beta_sub` scales like
`1/rho_vs`, and rho_vs falls steeply with temperature, so a FIXED v_n makes the
bound collapse at cold temperatures. With the default `v_n = 1e-9 m/s` on a
2 mm domain:

| T [C] | eps [um] | Nx | DOF | cores @ 50k |
|---|---|---|---|---|
| -5 | 2.14 | 1321 | 5.2e6 | 105 |
| -20 | 0.86 | 3295 | 3.3e7 | 652 |
| -30 | 0.32 | 8777 | 2.3e8 | 4623 |
| -40 | 0.11 | 25438 | **1.9e9** | **38826** |

That is an artifact: colder means *slower*, so v_n should fall with the
kinetics. `comp_eps.py --vn_feature R_feat` sets
`v_n = d0/(beta_sub * R_feat)` self-consistently, which makes the Eq. 45 bound
equal `R_feat` exactly and temperature-independent. `R_feat` is the smallest
curvature feature the mesh must resolve — for a sintering packing, the neck rim
radius.

| R_feat | eps (all T) | Nx | DOF | cores |
|---|---|---|---|---|
| 2 um | 1.00 um | 2829 | 2.4e7 | 481 |
| 5 um | 2.14-2.30 um | ~1280 | 5.0e6 | ~100 |

**Chosen: R_feat = 2 um -> eps = 1.00 um at every temperature.**

## 4. Why eps = 1 um and not the binary-percolation answer

A conservative reading — require the pore to still percolate after dilating the
solid by `3*eps` (the phi_a > 0.95 region) — gives:

| gap | max dilation before the pore disconnects | implied eps |
|---|---|---|
| 4 um | ~1.0 um | 0.33 um |
| 8 um | ~2.9 um | 0.99 um |

eps = 0.33 um would mean Nx ~ 8600 and ~4400 cores/run. But that criterion is
too strict: transport is not a percolation threshold here. The vapour equation
carries `D_eff = xi_v * D_v * phi_a`, so a throttled throat still conducts.
Evaluating phi_a at the centre of the *narrowest* throat:

| gap | critical throat | eps=0.33 | eps=1.0 | eps=2.0 |
|---|---|---|---|---|
| 4 um | 1.96 um | 0.90 | **0.45** | 0.24 |
| 8 um | 5.94 um | 1.00 | 0.90 | 0.63 |

At gap = 4 um and eps = 1 um the narrowest throats run at ~45% of free-vapour
conductance, and the median throat (p50 ~ 31 um) is fully open. Throttled, not
sealed.

**This is a starting point, not a converged answer.** Phase 4 (safety-factor /
xi_v calibration) decides empirically whether eps can be relaxed toward 2 um,
which would cut cost 4x but drops the narrowest throats to 24%. The gate is
k_eff and SSA convergence, not this geometric estimate.

## 5. Temperature sweep

T = -5, -10, -15, -20, -25, -30, -35, -40 C.

Every temperature needs its **own** experiment and geometry files. The solver
rescales `beta_sub0` and `d0_sub0` by `rho_vs(T_run)/rho_ice` at runtime, but
beta_HK's `1/sqrt(T)` and d0's `1/T` are baked into the numbers in the .opts —
reusing a -20 C file at -35 C updates the density ratio and nothing else.
`-eps_valid_temp` is emitted so the solver aborts on a mismatch rather than
running silently wrong.

`alpha_c = 1.341e-2` is held fixed across the sweep (it reproduces the committed
-20 C kinetics, `beta_sub0 = 5.9216e5`). Libbrecht's temperature-dependent
parameterization is deliberately NOT used — `flag_Tdep` stays 0. A bounded
Arrhenius fit is planned separately.

## Matrix

25 packings x 8 temperatures = **200 runs**, ~481 cores each at eps = 1 um.
Grain positions are temperature-independent and reused across all temperatures,
so the temperature effect is isolated; with 5 seeds at every (phi, T) cell there
are still independent realizations to separate trend from microstructure.
