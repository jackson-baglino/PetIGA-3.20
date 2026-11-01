# One-sided thin-interface correction with the temperature equation: first-order results

> **Provenance (added 2026-11-01).** Copied unchanged below this note from the
> practice-round archive, `scratch/2026-10-21_practice_round/research_archive/Planning/`
> (git-ignored, local only). Other numbered docs it cites (doc 13) are still there.
> Scripts under `agentic/math/` and the `agent/math-lane` branch are in tag
> `archive/2026-10-21/agent-math-lane`.

Status: 2026-10-18. Derived by the math lane (Sonnet derivations, Opus blind replicates, symbolic checks) and checked against exact numerical solutions of the full travelling-wave problem. Not yet in the solver. Task files: `agentic/math/tasks/c2*.json` on branch `agent/math-lane`; numerical audits: `agentic/math/numerics/heat_jump.py`, `front_bvp.py`.

## 1. The question in plain words
The phase-field equations move the ice-air front at a speed that depends on a kinetic coefficient. Because the interface has a finite width (eps), the simulated front behaves as if an extra resistance were present; the thin-interface correction cancels it so that the simulated front moves as the physical (sharp-interface) front would. Published corrections (Karma and Rappel 1998) assume the diffusing field has the same diffusivity on both sides of the front. In our problem two fields diffuse, and neither is two-sided with equal properties: vapor lives only in air, and heat flows in both ice and air with very different conductivities (K_ice about 2.3 and K_air about 0.02 W/m/K, a ratio of about 115 in the code). The question was what the first-order correction is for that situation.

## 2. Model analysed
One dimension, planar front moving at constant speed V, no sediment, no curvature. Ice fraction phi obeys the Allen-Cahn equation with the source (lam/rho_i) phi^2(1-phi)^2 (rho_v - rho_s(T)); vapor obeys ((1-phi) rho_v)_t = D((1-phi) rho_v,x)_x - rho_i phi_t; temperature obeys C(phi) T_t = (K(phi) T_x)_x + rho_i L phi_t, with K and C linear in phi as in `material_properties.c`; saturation density linearised, rho_s(T) = rho_s0 + s(T - T0). Scaling: tau = O(eps^2), everything else O(1). The front is where phi = 1/2; "front temperature" means the outer temperature extrapolated to the front from the air side.

## 3. Results (to first order in eps)
1. **Vapor channel: exactly zero correction.** Integrating the one-sided vapor equation once gives rho_v' = V (rho_i - rho_v)/D wherever phi is not 1, independent of phi. So vapor has no inner structure across the front, for any eps, and contributes nothing to the correction. (Both Sonnet and Opus found this; it agrees with doc 13 and the earlier gt_deficit measurement.)
2. **Heat channel.** At leading order temperature is continuous and Ki g_minus - Ka g_plus = rL V (rL = rho_i L; g_plus, g_minus are the outer gradients on the air and ice side). At first order the outer temperature jumps across the front by Delta_T = eps ln(Ki/Ka) (g_minus - g_plus), and the flux condition gets an O(eps) correction F1 = V (g_minus - g_plus)(Ca Ki - Ka Ci) ln(Ki/Ka)/(Ki - Ka). Both vanish when Ki = Ka (and F1 also when Ci/Ki = Ca/Ka). These reproduce the exact solution to 4-5 digits at five parameter sets.
3. **Kinetic law.** rho_v(x_f) - rho_s(T_plus) = eps K1, with
   K1 = 5 rho_i tau V/(lam eps^2) - 30 s (rL V - g_plus (Ki - Ka)) Xi(Ki, Ka),
   Xi = Ki^2 (Ki^2 - 5 Ki Ka + 10 Ka^2) ln(Ki/Ka) / (30 (Ki - Ka)^5) + Ka (Ki^2 - 6 Ki Ka + Ka^2) / (20 (Ki - Ka)^4).
   For Ki = Ka = K, Xi = 47/(1800 K): the thermal term becomes (47/60) s eps rL V/K, i.e. a1 a2 with a2 = 0.15667, the value in doc 13. The first term is the bare kinetic term (the 5 is a1). The part proportional to g_plus is new: it exists only when conductivities differ and depends on the outer temperature gradient, not on the speed.
   Checked against the full coupled travelling-wave problem (front_bvp.py: phase-field, heat and vapor solved together, no asymptotics) at five parameter sets with different Ki/Ka, outer gradients, speeds and tau: agreement within 0.5 % at eps = 0.005 and converging as eps shrinks.

## 4. What this means for the solver
- `tau_vap` should be removed (exact zero), as doc 13 argued.
- The thermal term in `tau_sub` is built on a2 = 0.1581 and the mean thermal diffusivity, which is the equal-conductivity form without the Clausius-Clapeyron weight s. The correct term is the one above: proportional to s L, with a conductivity-ratio function Xi. At the code's conductivities (ratio about 115) the V-proportional part is about 5 % of what the equal-conductivity formula would give using K_air, so the present term probably over-corrects; the conversion to the solver's own variables (beta_sub, units of s) is not done yet.
- The gradient-proportional term cannot be a constant in `tau_sub`, because it depends on the local temperature gradient at the interface. Options: ignore it (size to be assessed at lunar parameters), add a source term in the weak form proportional to the local gradient of T projected on the phase-field gradient, or choose a different reference point for the front temperature (ice side or conductivity-weighted) that makes it vanish. Which is best is a design decision for you; I have not changed any solver physics.

## 5. How it was verified (and what independence that gives)
Two models derived it independently from the same task text; an algebra system recomputed every step; the adversary found no high-severity objection in stages 1b and 2; and exact numerical solutions, which share no algebra with the derivation, matched. The stage-2 numerical solution does rely on the stage-1 temperature profile only through the equations themselves, not through the derived formulas. Not independent: the choice of model equations (K and C linear in phi, linearised saturation density) and the front definition (phi = 1/2, air-side extrapolation), which change the split between terms but not the physics.

## 6. Not settled
- Terms of order eps^2 times curvature (stretching, surface diffusion; where an anti-trapping current would act) are not analysed.
- Conversion of K1 into the solver's tau_sub and beta_sub, and the size of the gradient term at lunar parameters.
- The ice-sediment interface; the effect of the numerical air floor (1e-6) and the code's d0.
- Time-dependent (non-travelling) fronts: the derivation assumed a steady speed; stage 1 kept the outer time derivatives and found they cancel.

## 7. Update (2026-10-18, later): how big is the gradient term, and does a harmonic conductivity help?

### 7.1 Where the idea came from
The effective-conductivity work in `enceladus_DSM` (`-keff_interp`, `src/keff.c`; slides in `docs/keff_tensor_slides/`) found that an arithmetic conductivity K = phi K_i + (1 - phi) K_a gives an O(eps) error for heat flowing **across** a diffuse interface, and that the "tensor law" (arithmetic along the interface, harmonic across it, i.e. series resistances) reduces the error to O(eps^2). Our front problem is the same situation: heat crosses the interface. So the natural question was whether a harmonic law removes the gradient term. The solver's heat equation currently uses the arithmetic law; the tensor law exists only in the k_eff diagnostic.

### 7.2 Harmonic result (derived by the lane, replicated by Opus, checked against the full numerical solution)
With 1/K(phi) = phi/K_i + (1-phi)/K_a (the across-the-interface part of the tensor law, which is all that matters for a planar front):
- Kinetic law: K1_harm = 5 rho_i tau V/(lam eps^2) - (s/60) [ 47 g_plus (K_a/K_i - 1) + rL V (30/K_a + 17/K_i) ].
- The O(eps) temperature jump across the front depends on V only (it no longer depends on the outer gradient), as the series-resistance picture predicts.
- **The gradient term does not vanish**, but it is bounded: the phase field reads the air-side temperature displaced into the air by (47/60)(1 - K_a/K_i) eps, which is at most 0.78 eps for any conductivity ratio. For the linear law the displacement is 4.7 eps at the code's ratio of about 114 and keeps growing with the ratio.
- Checked against the exact coupled travelling-wave solution at four parameter sets: within 1 % at eps = 0.005. Equal conductivities (K_i = K_a) give the same result for both laws.

### 7.3 Sizing at project parameters (`agentic/math/numerics/sizing_lunar.py`; every input listed in the script)
Inputs: T = -20 C, K_ice = 2.29, K_air = 0.02 W/m/K, eps from the project's geometry files (0.37 to 0.86 um), gradients 0.3, 20, 50 K/m (the project's tgrad experiments), d0 from Kaempfer-Plapp Eq. 13 (1.0e-9 m), rho_vs = 8.8e-4 kg/m3, s = 8.1e-5 kg/m3/K.
- **Curvature driving force**: 1.1 mK-equivalent at a 10 um radius, 3.7 mK at 3 um.
- **Latent-heat (speed-proportional) term**: V is about 2e-11 m/s, which gives a temperature offset of order 1e-7 K, six orders below the driving force. **It is negligible, and so is the correction the solver's tau_therm currently applies on that basis.** The vapor term is exactly zero. So for isothermal runs the correct thin-interface correction is practically nil, and the existing tau_vap and tau_therm terms, which add to the kinetic coefficient, over-correct (the 1.22 times inflation noted in doc 13).
- **Gradient term, linear law**: 0 % of the curvature force at 0.3 K/m; 3 to 7 % at 20 K/m and 8 to 18 % at 50 K/m for a 10 um grain (a third of that for a 3 um grain). Relative to the gradient-driven signal itself over a 100 um pore: 1.7 to 4 %.
- **Gradient term, harmonic law**: about six times smaller (0.5 to 1.2 % at 20 K/m; 1.3 to 3 % at 50 K/m for a 10 um grain).

### 7.4 What this changes in the options
1. Ignoring the term is acceptable for isothermal runs and for gradients up to a few K/m; at 20 to 50 K/m it is a 3 to 18 % bias in the driving force with the linear law.
2. Switching the heat equation to the harmonic/tensor law cuts the bias about six-fold and has a precedent in this project, but for a general (non-planar) interface it needs the normal n = grad(phi)/|grad(phi)| (a tensor conductivity, with the same division by |grad phi| that already makes the Jacobian delicate in the diffuse band), a new analytic Jacobian, and checks. A harmonic law applied isotropically would also change conduction **along** the interface at O(eps), so the tensor form, not a simple replacement, is the consistent choice.
3. The local-gradient source correction (Option 2 of the earlier discussion) remains the most faithful but the most intrusive.
Not yet done: the curved-interface and ice-sediment cases, and whether the tensor law interacts with the grain-scale k_eff diagnostic.
