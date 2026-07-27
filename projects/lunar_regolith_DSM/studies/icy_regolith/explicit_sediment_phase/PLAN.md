# Effort 2 — Explicit sediment phase: design & plan

Paper 1, second approach. Model the regolith **explicitly** as a fourth,
static phase field φ_s (∂φ_s/∂t = 0) with a triple-well free energy over ice /
air / sediment. Unlike Effort 1 (regolith = domain boundary), this represents
interior regolith grains of arbitrary shape as a field, so ice can metamorphose
around fully-embedded grains. Branch: `exp/regolith-explicit-sediment-phase`.

> **Status:** planning. The governing residual/Jacobian are intentionally
> **left blank below** — the user will supply the beta-eliminated
> Lagrange-multiplier form. Do not implement the weak form until they are filled in.

---

## Why the prior attempt failed (and the fix)

The earlier three-phase attempt produced a **spurious air film wetting the
ice–sediment interface** and negative φ_ice at triple junctions
(`docs/ice_mass_loss_analysis.md`, `docs/spurious_ice_sed_air_branch_log.md`).
Diagnosed cause: the surface energies were chosen ad hoc to satisfy η_i, η_s,
η_a > 0, and γ_is = 0.033 was **below** the η_s > 0 floor → inverted sediment
well → complete wetting (which *is* the air film).

**Fix — parameterize by contact angle θ, not three free energies.** The η > 0
conditions are the triangle inequalities on (γ_iv, γ_is, γ_sv). With Young's law
γ_sv − γ_is = γ_iv·cos θ, two hold automatically:

- η_i = γ_iv(1 − cos θ) > 0  (θ > 0)
- η_a = γ_iv(1 + cos θ) > 0  (θ < 180°)
- η_s = 2γ_is − γ_iv(1 − cos θ) > 0  ⇔  **γ_is > (γ_iv/2)(1 − cos θ)**  ← the one real constraint

So the free parameter is γ_is expressed as a dimensionless margin above that
floor. At θ = 77.8°, floor = 0.0435 J/m²; the old γ_is = 0.033 violated it.
See `docs/material_parameters.md` §1.2–1.3.

**Action (pre-code):** add `gammas_from_contact_angle(gamma_iv, theta, gamma_is)`
to `preprocess/comp_eps.py`, returning the η's + margin above the floor, in the
style of its existing bound diagnostics. **Re-derive the γ→η mapping against the
beta-eliminated equations when supplied** — do not assume it carries over
unchanged from Kim–Steinbach.

## Governing equations (TO BE SUPPLIED)

The beta-eliminated Lagrange-multiplier formulation. 4 DOF per node:
φ_i (ice), T, ρ_v (vapor), φ_s (sediment, static). φ_a = 1 − φ_i − φ_s.

### Free energy / triple well

<!-- TODO(user): triple-well F(φ_i, φ_s, φ_a) and the η_i/η_s/η_a mapping from
     (γ_iv, γ_is, γ_sv). DSM's Fice/Fwat/Fair + Etai/Etam/Etaa in
     dry_snow_metamorphism/src/material_properties.c are the REFERENCE ONLY. -->

### Residual R[a][k]

<!-- TODO(user): beta-eliminated weak-form residual for
       k=0 ice (Allen-Cahn + triple well + sublimation source),
       k=1 temperature,
       k=2 vapor,
       k=3 sediment  (∂φ_s/∂t = 0 → identity/no-flux; confirm the exact weak form).
     Provide sign conventions and which terms carry the xi_v / xi_T scaling. -->

### Jacobian J[a][k][b][l]

<!-- TODO(user): analytic Jacobian blocks. Note the current 2-phase code omits
     no cross-terms; keep it fully analytic. The φ_s rows/cols are trivial if
     φ_s is frozen (∂R_s/∂φ_s = mass matrix, zero coupling), but confirm. -->

## Implementation (once equations are in)

Start from the **corrected 2-phase model** and re-add the 4th DOF — NOT from
`dry_snow_metamorphism` (pre-rewrite numerics + doc rot). Use DSM only as a
cross-check for limiting cases.

Files to touch:
- `include/NASA_types.h` — add `sed` to `Field` (→ 4 fields). The `-dof` guard
  (`permafrost2.c`, added on main) auto-derives the expected dof from
  `sizeof(Field)/sizeof(PetscScalar)`, so it will start accepting `-dof 4` the
  moment `Field` gains a fourth member — no guard edit needed.
- `src/permafrost2.c` — `dof = 4`; field name for component 3; register the
  γ/θ and sediment-geometry options; set `IGASetFieldName(iga, 3, "sediment")`.
- `src/assembly.c` — the new residual/Jacobian (the blanks above); read φ_s at
  Gauss points (like DSM's `user->Phi_sed`, but re-derived).
- `src/material_properties.c` — three-phase interpolation of every property
  (thcond, cp, rho, dif_vap) over (φ_i, φ_s, φ_a).
- `src/initial_conditions.c` — a sediment-field IC (fill φ_s from a grain list).

## Method discipline (learned from the prior failure)

1. **Simplest geometry first**: one sediment slab + one ice grain. The earlier
   attempt's complex initial geometry masked the parameter problem.
2. **Verify γ's before code**: confirm γ_is clears the η_s floor for the chosen
   θ with the new `gammas_from_contact_angle` helper.
3. **Gate on mass conservation**: run the 3-phase mass-conservation check (the
   2-phase system now passes at <0.01%). `postprocess/plot_mass.py` needs a φ_s
   term added. This is the measurement that would have caught the original bug.
4. **Resist compensating machinery**: no `t_sed_freeze`, penalty terms, or
   bounds-rollback stacks (all chronicled in the historical branch log). If they
   seem necessary, that is a signal to re-check the γ's, not to add crutches.
5. **eps/temperature**: emit `-eps_valid_temp` in generated geometries; the
   solver guard will catch mismatched runs.

## First milestones

- [ ] `gammas_from_contact_angle` in comp_eps.py + a printed margin diagnostic.
- [ ] Fill the residual/Jacobian blanks from the supplied equations.
- [ ] `Field` → 4 fields; `-dof 4` accepted by the guard; builds clean.
- [ ] Slab + one-grain IC; short run; 3-phase mass conservation < ~0.01%.
- [ ] No spurious air at the ice–sediment interface at the chosen θ.
