# Ice Mass Loss Analysis: Allen-Cahn Curvature Coarsening vs. Sublimation

**Date:** 2026-05-15  
**Run analyzed:** `test_2D_IceSedPair` (86400 s at T = −20 °C, h₀ = 0.95)

---

## Observed Mass Change

Over one simulated day the domain integrals changed as follows:

| Quantity | t = 0 | t = 86400 s | Change |
|----------|-------|-------------|--------|
| TOT_ICE (∫ φᵢ dA) | 1.125 × 10⁻⁹ m² | 9.243 × 10⁻¹⁰ m² | **−17.8%** |
| TOT_AIR (∫ φₐ dA) | 5.563 × 10⁻⁹ m² | 5.776 × 10⁻⁹ m² | +3.8% |
| TOT_SED (∫ φₛ dA) | 1.125 × 10⁻⁹ m² | 1.112 × 10⁻⁹ m² | −1.2% |
| TOT_RHOV (∫ ρᵥ φₐ dA) | 4.491 × 10⁻¹² kg/m | 4.703 × 10⁻¹² kg/m | +4.7% |

---

## Root Cause: Allen-Cahn Curvature Coarsening

### The two possible mechanisms

**Sublimation** (physical): ice converts to vapor at the ice-air interface driven by
the undersaturation (ρᵥ < ρᵥₛ). Mass leaves the ice phase and enters the vapor phase.
The expected vapor mass gain equals the ice mass lost.

**Allen-Cahn curvature coarsening** (numerical): the Allen-Cahn equation minimises
surface energy by shrinking curved grains (Ostwald ripening). The ice phase field
φᵢ → 0 due to curvature, but no physical mass constraint prevents it. The "lost" ice
mass does not appear as vapor.

### Mass balance diagnostic

If the mass loss were from sublimation, the vapor mass gained would equal the ice
mass lost:

```
Expected ΔTOT_RHOV (from ice→vapor) = ρᵢ × |ΔTOT_ICE|
                                     = 919 × 2.007×10⁻¹⁰
                                     = 1.85×10⁻⁷ kg/m
Actual ΔTOT_RHOV                     = 2.12×10⁻¹³ kg/m
Ratio                                = 1.15×10⁻⁶ ≈ ρᵥₛ / ρᵢ
```

The actual vapor mass gain is **ρᵥₛ/ρᵢ ≈ 10⁻⁶ times** the expected value.
This is exactly the ratio you get if the newly-created air is filled with
saturated vapor (ρᵥₛ, maintained by the k_pen penalty), rather than receiving
the full ice mass as vapor. The ice mass simply **disappears without entering
the vapor phase**.

In parallel: ΔTOT_RHOV ≈ ρᵥₛ × ΔTOT_AIR (ratio = 1.17 ≈ 1), confirming
that the new air regions are created already at saturation — not from
sublimation transferring mass.

### Interface velocity comparison

The Allen-Cahn curvature-driven interface velocity (sharp-interface limit):

```
v_AC = mob_sub × ηᵢ / (ε × R)
     = 1.9×10⁻⁹ × 0.086 / (7.12×10⁻⁷ × 1.875×10⁻⁵)
     = 13 m/s
```

The physical diffusion-limited sublimation rate at h₀ = 0.95:

```
v_phys = Dᵥ × (1−h₀) ρᵥₛ / (ρᵢ × R)
       = 2.2×10⁻⁵ × 0.05 × 8.49×10⁻⁴ / (919 × 1.875×10⁻⁵)
       = 5.4×10⁻⁸ m/s
```

**v_AC / v_phys ≈ 250,000,000.** The Allen-Cahn equation moves the interface
250 million times faster than physical sublimation can supply mass.
The sublimation source term S_sub = αₛᵤᵦ × φᵢ² φₐ² × (ρᵥ − ρᵥₛ)/ρᵢ
cannot couple the interface motion to vapor production at anywhere near this rate.

### Humidity diagnostic (h₀ = 1.0 test)

Setting h₀ → 1.0 eliminates the sublimation driving force entirely
(ρᵥ − ρᵥₛ = 0 everywhere). Ice **still shrank**, confirming that sublimation
is not the driver — curvature alone causes the mass loss.

---

## Parameters That Control This

### `d0_sub0` — Capillary Length (primary lever)

```
mob_sub = ε / (3 τ_sub)
τ_sub ∝ ε² λ β / a₁ ∝ ε / (d₀_sub × d₀_sub)
⇒ mob_sub ∝ d₀_sub²
```

Reducing `d0_sub0` decreases `mob_sub` and hence `v_AC`. Because v_AC/v_phys ≈ 2.5×10⁸,
closing the gap entirely would require reducing `d0_sub0` by ~10⁴×, which is likely
too extreme. A 10–100× reduction (1e-11 to 1e-10 m) is a practical compromise that
slows coarsening to a rate where short simulations are meaningful.

**Result after tuning:** mass loss was no longer noticeable in the mass plot.
A new artifact appeared: slight negative φᵢ at the sediment-air triple junction.
This is a separate issue (see §below).

### `beta_sub0` — Kinetic Coefficient

Also enters `τ_sub` and inversely scales `mob_sub`. Increasing `beta_sub0` has a
similar effect to reducing `d0_sub0`. The two parameters are not independent —
both contribute to `τ_sub` through:

```
τ_sub = ε λ (β_sub/a₁ + a₂ ε/D_th + a₂ ε/Dᵥ)
```

In the current parameter regime, the `β_sub` term dominates.

### Surface Tensions (`γ_iv`, `γ_is`, `γ_sv`)

The Allen-Cahn curvature driving force scales with `ηᵢ = γ_iv + γ_is − γ_sv`.
Reducing `γ_iv` reduces ηᵢ and hence v_AC, but not enough to close the 8-order
gap. Surface tensions primarily control:

- **Contact angle** (ice wetting sediment): cos θ = (γ_sv − γ_is)/γ_iv
- **Triple-junction geometry** and energy
- **Relative driving forces** between phases

They are the correct parameters to tune for **equilibrium geometry**, not for
mass conservation.

| Parameter | Current | Effect of lowering γ_iv (e.g. 0.109 → 0.05) |
|-----------|---------|---------------------------------------------|
| ηᵢ | 0.086 | 0.027 (3× less curvature drive) |
| ηₛ | −0.020 | −0.003 (less negative) |
| ηₐ | 0.132 | 0.073 |
| Contact angle θ | 77.8° | 62.6° |
| v_AC reduction | — | 3× |

### Allen-Cahn Mobility Override

Adding a `-mob_sub_scale` factor to permafrost2.c is the most direct knob:

```c
user.mob_sub = mob_sub_scale * eps / 3.0 / tau_sub;
```

Setting `mob_sub_scale = 1e-4` reduces v_AC by 10⁴×, making it comparable to
physical sublimation rates for the grain sizes simulated.

---

## Negative φᵢ at the Sediment-Air Interface

After reducing mob_sub, slight negative φᵢ appeared at the sediment-air (triple)
junction. This is distinct from the mass-loss problem and likely has one of two causes:

1. **Triple-junction energy imbalance**: when `ηₛ < 0`, the Kim-Steinbach
   double-well for sediment is inverted, which can cause the Allen-Cahn equation
   to "push" φᵢ below zero at the sed-air junction rather than pull it toward a
   contact angle. This is consistent with the current values (ηₛ = −0.020).

2. **Curvature artifact at sharp sediment boundary**: the sediment grain is
   held fixed (frozen) while ice evolves freely. At the sed-air contact line,
   the no-flux condition on φₛ combined with the Allen-Cahn driving force on φᵢ
   can create a slight undershoot when mob_sub is very small (the interface has
   no mechanism to "escape" the corner).

Disabling the phase-change term (S_sub = 0) and observing the same negative φᵢ
suggests the AC dynamics alone drive the undershoot, ruling out the sublimation
source as the cause.

**Recommended fix**: Increase `γ_iv` relative to `γ_sv` to make ηₛ less negative
(or zero), which removes the inverted double-well for sediment at the triple junction.

---

## Summary of Recommendations

| Goal | Parameter | Direction | Notes |
|------|-----------|-----------|-------|
| Reduce AC coarsening rate | `d0_sub0` | ↓ (e.g. 1e-9 → 1e-11) | Primary lever |
| Same | `beta_sub0` | ↑ | Secondary lever |
| Correct contact angles | `γ_iv`, `γ_is`, `γ_sv` | match physical values | See material_parameters.md |
| Fix triple-junction undershoot | `γ_sv` ↑ or `γ_is` ↓ | make ηₛ ≥ 0 | Geometry fix |
| True mass conservation | Anti-trapping current | Code change | Correct long-term fix |
| Diagnostic baseline | `humidity` → 1.0 | n/a | Isolates AC vs. sublimation |

### Long-term fix

True mass conservation in sublimation phase-field models requires implementing
the **anti-trapping current** (Karma & Rappel 1998) extended to the vapor phase.
Without this, the Allen-Cahn interface velocity and the sublimation mass flux are
only decoupled in the limit where mob_sub → 0 (no AC dynamics) or mob_sub → ∞
(interface moves purely by sublimation). Neither extreme is physical.

The Karma-Rappel anti-trapping term adds a correction flux to the vapor equation
that exactly cancels the spurious mass trapping in the diffuse interface region,
enforcing the Stefan condition without requiring mob_sub to be asymptotically small.
