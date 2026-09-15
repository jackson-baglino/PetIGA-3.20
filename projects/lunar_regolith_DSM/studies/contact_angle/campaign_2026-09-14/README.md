# Wedge-scale campaign — batch 2026-09-14

15 runs, all numerically healthy: 61/61 distinct snapshots each, zero bounds
violations, zero divergences, all completed.

## Verdict

**The angles are fine. The growth experiment (group B) failed for a reason
worth knowing, and groups A and C stopped short of convergence.**

## 1. Equilibrium angle — ✅ both geometries

| θ_Young | channel (A) | droplet (C) |
|---|---|---|
| 30° | +0.415 | +1.200 |
| 60° | +0.098 | +0.061 |
| 90° | +0.018 | −0.001 |
| 120° | −0.010 | +0.011 |
| 150° | +0.013 | −0.106 |
| **RMS** | **0.191** | **0.540** |

Both well inside the 2° criterion, and the droplet agrees with the channel — so
the angle is set by the boundary condition, not by confinement. That was the
control the previous batch never managed to deliver.

## 2. Group B (growth vs shrinkage) — ❌ thermally throttled

The ice changed by **0.001–0.003 %** in 90 days, in every case, with no
dependence on the sign of χ. Predicted: tens of percent.

**Cause: the domain is thermally insulated** (`-flag_BC_Tfix 0`). Measured mean
temperature drift over the run, against the ΔT that exactly cancels the
Gibbs-Thomson supersaturation `d0·χ`:

| θ | χ | measured ΔT | ΔT that cancels d0·χ |
|---|---|---|---|
| 30° | −1.39e4 | **+1.493e-4 K** | +1.472e-4 |
| 60° | −8.00e3 | +8.522e-5 | +8.499e-5 |
| 90° | 0 | −1.2e-8 | 0 |
| 120° | +8.00e3 | −8.470e-5 | −8.499e-5 |
| 150° | +1.39e4 | −1.467e-4 | −1.472e-4 |

Agreement to ~1% at every angle, with the right sign throughout. The mechanism:

- a concave meniscus wants to grow, deposition releases latent heat,
- the insulated box warms, ρ_vs(T) rises,
- supersaturation falls to zero and growth stops.

Depositing **1e-5 %** of the ice releases enough heat to warm the box by
1.5e-4 K, which moves ρ_vs by exactly the 1.4e-5 supersaturation that was
driving it. The driving force annihilates itself almost immediately.

Confirmed independently: the vapour throughout the domain — including at the
"reservoir" faces — sits at `ρ_vs(1 + d0·χ)`, i.e. the whole field has
equilibrated to the meniscus rather than to the imposed boundary value.

**The sign of the drift is the physics you predicted** — concave wants to grow,
convex wants to shrink, 90° does nothing. It is throttled, not absent.

**Fix applied**: `-flag_BC_Tfix 1` added to the `grow90_*` experiment files.
With `-grad_temp0` all zero the solver pins T = temp0 on every face, giving a
thermal bath alongside the vapour reservoir. Growth needs *both* reservoirs:
mass has to come from somewhere and latent heat has to go somewhere. Pinning T
on the y-walls is harmless, unlike pinning vapour there — regolith is a thermal
bath, and it does not feed ice at the contact line.

## 3. Convergence — ⚠️ both groups short

| group | τ | t_final | t_final/τ |
|---|---|---|---|
| A/B channel | 23.5–34.0 d | 90 d | 2.6–3.8 |
| C droplet | 52.1–71.3 d | 174 d | 2.4–3.3 |

The previous 100 µm batch reached 4.2 τ and gave RMS 0.062°; this one reaches
2.6–3.8 τ and gives 0.191°. Nothing is wrong — the larger channel simply relaxes
more slowly (τ 15 → 24–34 d) and 90 days no longer buys the same convergence.
The worst case, θ=30, has a fit-window drift of 0.648° against 0.011° at θ=120.

To match the previous precision: channel **153 d** (1.32e7 s, ~7.9 h/run),
droplet **321 d** (2.77e7 s, ~16.5 h/run).

## 4. Housekeeping

The batch's own postprocessing did not produce `contact_angle.csv` for any run,
although the gate and the script were both staged correctly. Re-run locally
without trouble. Worth checking why before the next batch.
