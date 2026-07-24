# Effort 2 — Explicit sediment phase (dof=4): attempt log and why it was set aside

> ⚠️ **SET ASIDE (banner 2026-07-24).** This documents the explicit-sediment
> 3-phase route (branch `exp/regolith-explicit-sediment-phase`). The solver is
> correct and verified, but the route was paused in favor of the *implicit*
> regolith method (Effort 1) because of an intrinsic triple-junction limitation
> described in §4. The branch is left in its final experimental state (VI solver
> on, grain placed at its equilibrium contact angle); it is not deleted.

**Goal.** Model dry-snow metamorphism in icy lunar regolith by adding a 4th
degree of freedom `phi_s` (sediment) with `d phi_s/dt = 0` (frozen) and a
triple-well free energy, so ice, air/vapor, and sediment coexist explicitly.
This is the *explicit* counterpart to Effort 1, where the regolith instead
enters only as pore-space domain geometry (3 DOF, no sediment field).

---

## 1. What was built and verified

- `src/assembly_sed.c`: `Residual_A2` + analytic `Jacobian_A2` for 4 DOF
  (0=ice, 1=T, 2=rhov, 3=sed). Reduces exactly to the 2-phase kernel at
  `phi_s = 0`.
- `src/material_properties.c`: `ThermalCond3`/`HeatCap3`/`Density3` (mixtures
  over ice/sed, air = 1−ice−sed) and `TripleWell` (Adrian's `F^tri`, with the
  `Lambda` triple-junction penalty term).
- Dispatch on `dof`: `Residual`/`Jacobian` route to `_A2` when `dof==4`, else
  `_A1`. `MAX_DOF` (=4) sizes the shared buffers (no `sol[dof]` VLA).
- Surface energies parameterized so `Sigma_i`, `Sigma_a`, `Sigma_s` clear the
  `Sigma_s > 0` floor via `comp_eps.py gammas_from_contact_angle`.

**Verification (2026-07-22).** `-snes_test_jacobian` gave
`||J − Jfd|| / ||J|| ~ 1e-10` (threshold 1e-8) with sediment coupling active;
first-solve fnorm ~1e-12. Mass is conserved. The A1 (2-phase) and A2 (3-phase
at `phi_s=0`) kernels agree to ~1e-9. **The solver is correct** — the problem
below is physical/representational, not a solver bug.

---

## 2. The symptom: spurious air at the ice–sediment contact

With a real sediment slab present, a small, bounded, mass-conserving band of
`phi_a` (air) appears at the ice–sediment interface where none should exist. It
is **not** the catastrophic unbounded-air failure of the earlier scrapped
attempt (that came from an inverted sediment double-well, `Sigma_s < 0`; the
corrected `Sigma_a > 0 = air non-wetting` set removed it). It is a small
third-phase adsorption / triple-junction artifact, ~1.4e-2 in magnitude.

Two distinct locations, with different causes:

- **Deep ice–sed interface (center line):** an antisymmetric air profile
  (positive on the sediment side, zero at `phi_i = phi_s = 0.5`, negative on the
  ice side) — a *width mismatch* (the evolved ice interface is sharper than the
  frozen sediment interface).
- **Triple junctions (waterline corners, where ice/air/sediment meet):**
  both-sign violations, the largest excursions.

---

## 3. Fixes tried, in order

| # | Fix | Result |
|---|-----|--------|
| 1 | **`Lambda` triple-junction penalty sweep** {0,1,10,100} | `Lambda=1` is the sweet spot (least air while phases stay in contact). `Lambda>=10` over-stiffens and opens a spurious air *film* separating ice from sediment. Only ~12% improvement — near its useful ceiling. |
| 2 | **`-sed_width_factor`** — sharpen the frozen `phi_s` interface to match the (narrower) evolved ice interface | Helps the antisymmetric center-line air (addresses the width mismatch directly). |
| 3 | **`eps/2` hi-res run** (h/eps fixed) | Peak air ~unchanged, only the diffuse band thins ⇒ the center-line mismatch is a genuine *width-ratio* mismatch, **not** under-resolution. Confirms a finer mesh alone won't fix it. |
| 4 | **Carving construction** `phi_s = slab*(1−phi_i)` | Makes `phi_a = (1−ice)(1−slab) = 0` at the contact at t=0. Cleaned the deep interface; the violations **relocated** to the triple-junction corners. |
| 5 | **Equilibrium contact-angle placement** — put the grain centre at `cy = h_sed − RCice*cos(theta)`, `theta` from Young's law on the Sigmas. For `Sigma_i = Sigma_a` this is `theta=90°` (hemisphere, centre on the slab). Replaces the ad-hoc `ice_inset_frac=0.33`, which had imposed `theta ~ 132°` and forced a ~42° relaxation transient. | Removed the *placement*-driven relaxation, but the triple junction still overshoots. |
| 6 | **VI bound-constrained solver** (`vinewtonssls`, strict `[0,1]` on the ice DOF) | See §4 — this is where the route was stopped. |

---

## 4. Why the route was set aside (the VI run, 2026-07-24)

Run: `2D_ice_sed_contact / 1day_T-20_h1.00 / vi_theta_eq`
(`.../2026-07-24__10.58.49_1day_T-20_h1.00_vi_theta_eq`), VI on, `theta=90°`.

- VI held the **ice DOF** in `[0,1]` exactly; mass conservation was excellent
  (−0.002%) for the first 26 steps.
- At **step 26** (t ≈ 0.42 s of a 1-day run) the solve hit
  **`DIVERGED_FUNCTION_DOMAIN` at Newton iteration 0** — the residual could not
  even be evaluated at the current state. `dt` then spiraled down
  (2.1e-2 → 1.58e-2 → 1.19e-2 …) without recovering, and the run stalled.

**Root limitation.** VI (`SNESVISetVariableBounds`) can only constrain *solver
DOFs*. In this model the DOFs are ice, T, rhov, and (frozen) sediment. **Air is
a derived field, `phi_a = 1 − phi_i − phi_s`, not a DOF**, so it cannot be
VI-constrained. Even with `phi_i` clamped to `[0,1]`, the derived air overshoots
below the pointwise domain guard at the triple-junction quadrature points, and
`assembly.c`'s domain check rejects that state → `DIVERGED_FUNCTION_DOMAIN`.
Bounding the sum `phi_i + phi_s ≤ 1` is a constraint on a *combination* of DOFs,
which the simple per-stride VI box cannot express.

So the triple-junction air overshoot is **intrinsic to the frozen-sediment,
derived-air representation**, not something the numerics can be tuned around
without a fundamentally different constraint formulation (e.g. an obstacle /
constrained-sum multiphase solver). That is a much larger investment than the
implicit route, so the explicit route is paused here.

---

## 5. Final state of the branch (kept, not reverted)

- **VI solver ON** in `inputs/solver.opts` (`-vi_bounds 1`, `-vi_lo 0.0`,
  `-vi_hi 1.0`, `-snes_type vinewtonssls`). The earlier reason VI had been
  disabled — a Gibbs–Thomson `1/|∇phi|^3` Jacobian singularity — no longer
  applies because the GT curvature term was removed on 2026-07-21.
- **dof=4 VI stride-3 fix** (`src/permafrost2.c`): the VI bounds now set the
  frozen sediment stride to `(−inf, +inf)`. Previously it was left unset, and
  `IGACreateVec` zero-fills, so turning VI on with `dof=4` would have pinned
  `phi_s` to `[0,0]` and erased the sediment slab. (Latent bug, fixed here.)
- **Contact-angle grain placement** (`src/initial_conditions.c`,
  `-contact_angle_deg`): grain centre at `cy = h_sed − RCice*cos(theta)`,
  `theta` from `-contact_angle_deg` or Young's law on the Sigmas. The old
  `-ice_inset_frac` knob is deprecated (accepted but unused).

---

## 6. If this route is revived later

The blocker is representational, so the fix must change *how the phases are
constrained*, not the mesh or the Sigmas:

- A **constrained-sum / obstacle multiphase** formulation that enforces
  `phi_i, phi_s, phi_a ≥ 0` and `sum = 1` as a simplex constraint, rather than
  treating air as `1 − phi_i − phi_s` with only ice box-bounded.
- Or a genuinely **variational multiphase-field** model (N order parameters
  with a projected/normalized update) instead of the eliminated-air form.

Neither is a small change. For the icy-regolith paper, Effort 1 (implicit
pore-domain, 3 DOF) reaches the same science without a frozen-phase triple
junction, and is the active path as of 2026-07-24.
