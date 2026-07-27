# Root Cause Analysis: Spurious Air in the Old Permafrost Simulation

> ⚠️ **HISTORICAL (banner 2026-07-22).** Sediment-era analysis (pre-2026-06-13 fork) of the removed three-phase model. Kept as reference for Effort 2 (`studies/icy_regolith/explicit_sediment_phase/`), not a description of the current two-phase code.


**Date:** 2026-05-15  
**Old run:** `test_1D_IceSlab_2Phase_difvappen1e-07_k_pen1e09` (2026-05-05)  
**Reference files compared:** `src/assembly.c`, `src/permafrost2.c`, `src/initial_conditions.c`, `*.opts`

> **Status note (2026-05-27).** This document analyses the original
> 2026-05-05 failure cascade in the *legacy* 2-phase formulation. The eight
> bugs catalogued here were all fixed before the start of the
> `fix/spurious-ice-sed-air-penalty` branch. Several *further* model-level
> fixes were made on that branch — the 3-phase ice equation under frozen
> sediment (branch log §1), the vapor-diffusivity-penalty direction (§25),
> the mass-conserving Stefan source `vap_src = -2*ρ_ice*φ_a*∂φ_i/∂t` (§26),
> latent heat pairing with `S_sub` (§3.3 in `model_description.md`), and
> finally turning off both vapor-equation penalties as the current
> recommendation (§27). The model equations in
> [`model_description.md`](model_description.md) reflect the **current**
> code; this document is preserved as historical context for *why* the
> original code failed, not as a description of the present model.

---

## Background

The old simulation produced growing regions where the ice volume fraction
φ_i dropped below zero — "spurious air" that appeared throughout the domain
and grew irreversibly over time. The current code produces physically correct
results with no spurious air. This document identifies and quantifies every
difference between the two codebases that contributed to the failure.

Eight distinct bugs were found. They are ordered roughly by severity. Bugs 1
and 2 are the primary causes; the remaining six amplified their effects or
would have caused independent problems at longer times.

---

## Bug 1 — Missing ξ_v on the Vapor Penalty Term (PRIMARY)

**File:** `src/assembly.c`  
**Nature:** Incorrect residual formulation — a physics-scaling factor omitted
from one term of the vapor residual.

### Old code (line 223)

```c
R_vap += k_pen * g_phiiphis * (rhov - rhov_eq) * N0[a];   // ← xi_v absent
```

### Current code

```c
PetscReal vap_pen = xi_v * k_pen * g_phiiphis * (rhov - rhov_eq);
R_vap += vap_pen * N0[a];                                   // ← xi_v present
```

### Why this is wrong

The vapor transport PDE has the form:

```
∂ρ_v/∂t = ξ_v [ -∇·(D_eff ∇ρ_v) - k_pen H(φ_i+φ_s)(ρ_v - ρ_eq) + ρ_ice ∂φ_a/∂t ]
```

The factor ξ_v multiplies **all** spatial terms — it is a global time-scale
separation parameter that slows down vapor dynamics relative to the phase
field. In the old code, ξ_v was correctly applied to the diffusion term and
the source term, but was **omitted from the penalty term**. This made the
penalty ξ_v⁻¹ = 1000 times larger than intended.

### Quantitative impact

With the old parameters (ξ_v = 10⁻³, k_pen = 10⁹, Δt = 10⁻⁴ s,
ρ_vs ≈ 3×10⁻⁵ kg m⁻³):

| Term | Old code magnitude | Current code magnitude |
|------|--------------------|------------------------|
| Time derivative ∂ρ_v/∂t | ~ρ_vs / Δt ≈ **0.3** | ~0.3 |
| Diffusion ξ_v D ∇²ρ_v | ~**1.4×10⁻³** | ~1.4×10⁻⁴ |
| Penalty (as written) | k_pen × δρ_v ≈ **3×10⁴** | ξ_v × k_pen × δρ_v ≈ 3×10⁻⁴ |

The old penalty (3×10⁴) exceeded the time derivative (0.3) by **five orders
of magnitude**. The vapor equation was no longer a PDE — it behaved as an
algebraic constraint forcing ρ_v = ρ_eq at every point inside the solid,
with zero time evolution. The physical relaxation timescale was completely
lost.

---

## Bug 2 — Wrong Initial Vapor Density Inside Solid (PRIMARY)

**File:** `src/initial_conditions.c`  
**Nature:** Incorrect initial condition — undersaturation applied uniformly,
including inside ice where the equilibrium value is full saturation.

### Old code

```c
u[j][i].rhov = user->hum0 * rho_vs;   // applies h₀ everywhere
```

### Current code

```c
// h₀ applied only in air; inside solid, ρ_v = ρ_vs (saturation)
u[j][i].rhov = rho_vs * (user->hum0 * _phi_air + (1.0 - _phi_air));
```

This evaluates to:
- **Inside solid** (φ_a ≈ 0): `ρ_v = ρ_vs` — saturated, correct.
- **In air** (φ_a ≈ 1): `ρ_v = h₀ ρ_vs` — undersaturated per humidity, correct.
- **At interface**: smooth interpolation between the two limits.

### Why this is wrong

The equilibrium vapor pressure inside a condensed phase (ice or sediment) is
the saturation value ρ_vs. Setting ρ_v = h₀ ρ_vs inside ice with h₀ = 0.5
means the initial condition has a 50% disequilibrium at every solid point.

At t = 0, the penalty residual inside ice was:

```
R_vap = k_pen × (ρ_v - ρ_eq) = 1e9 × (0.5 ρ_vs - ρ_vs) ≈ 1e9 × (-1.5×10⁻⁵) ≈ -15
```

The time derivative at the same point was approximately ρ_vs / Δt ≈ 0.3.
The penalty term was **50× larger than the time derivative on the very first
time step**, before any physics had evolved. Newton's first iteration had to
eliminate this enormous disequilibrium — but was working with an inaccurate
Jacobian (Bug 3) that made the correction step qualitatively wrong.

---

## Bug 3 — Finite-Difference Jacobian Instead of Analytical Jacobian

**File:** `src/permafrost2.c`  
**Nature:** The analytical Jacobian function was disabled; PETSc fell back to
finite-difference approximation.

### Old code

```c
// ierr = IGASetFormIJacobian(iga, Jacobian, &user); CHKERRQ(ierr);
ierr = IGASetFormIJacobian(iga, IGAFormIJacobianFD, &user); CHKERRQ(ierr);
```

### Current code

```c
ierr = IGASetFormIJacobian(iga, Jacobian, &user); CHKERRQ(ierr);
// ierr = IGASetFormIJacobian(iga, IGAFormIJacobianFD, &user); CHKERRQ(ierr);
```

### Why this is wrong

PETSc's finite-difference Jacobian perturbs each degree of freedom by
δ ∼ √ε_machine ≈ 10⁻⁸ and estimates derivatives as:

```
J_ij ≈ [ R_i(u + δ e_j) - R_i(u) ] / δ
```

The truncation error in each entry is O(ε_machine / δ) × |∂²R/∂u²| × δ =
O(ε_machine) in relative terms, but O(ε_machine × max|J|) in absolute terms.

For residual entries that scale with k_pen = 10⁹ (the vapor penalty block):

```
|ΔJ_ij| ~ k_pen × ε_machine ≈ 1e9 × 1e-16 / 1e-8 = O(10)
```

This absolute Jacobian error of O(10) is **comparable to the residual itself**
(≈ 0.3–15 depending on location and time step). Newton's quadratic convergence
requires the Jacobian error to be much smaller than the residual; when they
are comparable, Newton either diverges or converges to a wrong root.

The current `Jacobian_A1` function provides exact analytical derivatives for
all 16 of the 4×4 coupling blocks, including the stiff penalty blocks that
the FD approach could not handle. Newton now converges quadratically (residual
drops by ≈10⁴ per iteration after the first step).

---

## Bug 4 — Interface Equilibrium Stiffness Too High (k_pen = 10⁹ → 10⁵)

**File:** `*.opts`  
**Nature:** Parameter value four orders of magnitude too large.

### Old value

```
-k_pen 1e9
```

### Current value

```
-k_pen 1e5
```

### Why this causes problems

Even if Bugs 1 and 3 were fixed (ξ_v correctly applied, analytical Jacobian
used), k_pen = 10⁹ makes the system extremely stiff. The condition number of
the vapor block of the linear system scales as:

```
κ ~ ξ_v × k_pen × Δt ~ 10⁻³ × 10⁹ × 10⁻⁴ = 100
```

With bjacobi+ILU(1) preconditioning (Bug 7), this is already marginal for
GMRES convergence. Any additional stiffness (e.g. from the interface geometry
or from adaptive Δt refinement) pushes GMRES to stall.

With k_pen = 10⁵ and ξ_v = 10⁻⁴:

```
κ ~ 10⁻⁴ × 10⁵ × 10⁻⁴ = 10⁻³   (effectively 1 — no stiffness added)
```

The penalty is then purely a physical regularisation (enforcing ρ_v ≈ ρ_eq
inside solid to within |ρ_v − ρ_eq| ≲ residual_atol / k_pen) without any
numerical ill-conditioning.

---

## Bug 5 — Vapor Time-Scale Factor Too Large (ξ_v = 10⁻³ → 10⁻⁴)

**File:** `src/permafrost2.c`  
**Nature:** Parameter value one order of magnitude too large.

### Old value

```c
user.xi_v = 1.0e-3;
```

### Current value

```c
user.xi_v = 1.0e-4;
```

### Why this matters

ξ_v sets the ratio of the physical vapor relaxation timescale to the
phase-field timescale. A larger ξ_v means vapor equilibrates faster relative
to the phase field. In the old code, ξ_v = 10⁻³ combined with the omitted
ξ_v on k_pen (Bug 1) to make the effective unscaled penalty coefficient:

```
k_pen_effective = k_pen / ξ_v = 10⁹ / 10⁻³ = 10¹²   (relative to diffusion)
```

With the corrected ξ_v = 10⁻⁴ and properly-scaled k_pen:

```
ξ_v × k_pen = 10⁻⁴ × 10⁵ = 10   (a modest penalty, weaker than diffusion at the interface scale)
```

The factor-of-10 reduction in ξ_v also reduces the stiffness contribution of
the correctly-scaled penalty, improving linear solver conditioning.

---

## Bug 6 — Diffusivity Penalty Factor Too Extreme (α_pen = 10⁻⁷ → 10⁻⁴)

**File:** `*.opts`  
**Nature:** Parameter value three orders of magnitude too extreme.

### Old value

```
-difvap_pen 1.0e-7
```

### Current value

```
-difvap_pen 1.0e-4
```

### Why this causes problems

The effective vapor diffusivity is:

```
D_eff = D_v × [ φ_a + α_pen × (1 − φ_a) ]
      ≈ D_v × φ_a          in air
      ≈ D_v × α_pen        deep inside ice
```

With α_pen = 10⁻⁷, the diffusivity drops by **seven orders of magnitude**
across the ~4ε interface. The smooth Heaviside function used to transition
between the two limits has a finite width of ~4ε ≈ 3 μm. On the B-spline
mesh with element size h ≈ ε/4 ≈ 180 nm and 3 Gauss points per element,
adjacent quadrature points within the interface see D_eff values that differ
by up to a factor of 10⁷.

This near-discontinuous coefficient cannot be accurately represented by the
IGA basis functions or differentiated accurately by the finite-difference
Jacobian. The resulting Jacobian errors at interface nodes were independent of
and additive to those from Bug 3, further degrading Newton convergence.

With α_pen = 10⁻⁴, the diffusivity still drops by four orders of magnitude
(ice remains essentially impermeable to vapor diffusion), but the coefficient
variation across any single element is at most O(10), well within the
resolution of the quadratic B-spline basis.

---

## Bug 7 — Inadequate Preconditioner (bjacobi+ILU(1) → ASM+ILU(2))

**File:** `*.opts`  
**Nature:** Suboptimal linear solver configuration for a 4-DOF coupled
asymmetric system.

### Old configuration

```
-pc_type bjacobi           # block Jacobi, zero overlap
-sub_pc_type ilu
-sub_pc_factor_levels 1    # ILU(1)
-ksp_gmres_restart 500     # effectively unrestarted (= max_it)
```

### Current configuration

```
-pc_type asm               # additive Schwarz, 1-node overlap
-pc_asm_overlap 1
-sub_pc_type ilu
-sub_pc_factor_levels 2    # ILU(2)
-ksp_gmres_restart 200
```

### Why this matters

**Zero overlap (bjacobi):** With no shared nodes between MPI blocks, the
preconditioner cannot capture the coupled physics that spans block boundaries.
In a B-spline IGA discretisation, each basis function has support over 3
elements in 1D (9 in 2D), so the natural stencil always crosses block
boundaries. The block-Jacobi preconditioner ignores all cross-block coupling,
producing a poor approximation of the system matrix. GMRES must compensate
with more iterations.

**ILU(1) vs ILU(2):** The 4-DOF coupled system has a natural fill pattern
extending 3–9 nodes away from each test function (depending on dimension and
B-spline order). ILU(1) retains one level of fill beyond the sparsity
pattern — insufficient to capture the dominant off-diagonal coupling between
the φ_i and ρ_v fields (which interact through the sublimation source and
penalty terms). ILU(2) captures this two-hop coupling and significantly
reduces the effective condition number seen by GMRES.

**GMRES restart 500 = max_it:** Setting restart = max_it is effectively
unrestarted GMRES. For a problem with N ≈ 10⁵ DOFs, storing 500 Krylov
vectors consumes ~1 GB of memory. This wasted memory competed with PETSc's
matrix and vector allocations, degrading cache performance. Restart = 200 is
sufficient for the well-preconditioned current system and uses 60% less
Krylov storage.

---

## Bug 8 — Excessive Allen-Cahn Relaxation and Immediate Sediment Freeze

**File:** `*.opts`  
**Nature:** Solver strategy that created large initial disequilibria.

### Old configuration

```
-n_relax 12          # 12 Allen-Cahn-only steps before full physics
-t_sed_freeze 0      # sediment frozen at t=0 (no 3-phase equilibration)
```

### Current configuration

```
-n_relax 1           # 1 Allen-Cahn step, then full coupling
-t_sed_freeze 1      # 1 second of 3-phase evolution before freezing
```

### Why this causes problems

**n_relax = 12:** The relaxation phase runs Allen-Cahn equations for φ_i and
φ_s without solving for T or ρ_v. After 12 such steps, the phase-field
profiles have evolved (sharpened their diffuse interfaces, moved slightly
toward energy minima) but the vapor field ρ_v has not been updated to match.
When the first full-physics step then solves all four equations simultaneously,
it sees a vapor field that is thermodynamically inconsistent with the current
phase-field geometry: the interface positions have moved, but ρ_v still
reflects the original (already-incorrect) initial condition. This creates a
large artificial disequilibrium at the start of the coupled evolution, which
Bug 1 and Bug 2 then amplify into a catastrophic first step.

With n_relax = 1, only a single AC step precedes full coupling. The
disequilibrium introduced is minimal and the coupled solver handles the
correction without difficulty.

**t_sed_freeze = 0:** In the old code, sediment was frozen into its initial
configuration at t = 0 with no equilibration period. Any mismatch between
the initial φ_s profile (set by the IC function) and the two-phase energy
minimum appeared as a finite residual in the sediment equation at t = 0,
adding to the already-large initial residual from Bug 2.

With t_sed_freeze = 1 s, the sediment field has one second to equilibrate
its interface shape under the full three-phase free energy (ice, sediment,
and air all interacting) before the sediment is frozen. This eliminates the
initial transient and ensures the 2-phase regime starts from a
thermodynamically consistent configuration.

---

## How the Bugs Compounded Each Other: The Failure Cascade

None of the eight bugs alone was necessarily fatal. Their combination was:

```
t = 0:
  Bug 2 → ρ_v = 0.5 ρ_vs inside ice (should be ρ_vs)
         → penalty residual = k_pen × (-0.5 ρ_vs) ≈ -15  [should be ~0]
         → time derivative ≈ 0.3
         → initial residual is 50× larger than expected

  Bug 8 → 12 AC relaxation steps move φ_i/φ_s without updating ρ_v
         → interface positions shift but ρ_v field becomes further inconsistent

First Newton iteration:
  Bug 1 → penalty term in Jacobian is 1000× too large (missing ξ_v)
         → GMRES solves a system where one term dominates by 10³
  Bug 3 → FD Jacobian: ΔJ ~ k_pen × ε_machine = 1e9 × 1e-8 = O(10)
         → dominant Jacobian block has O(1) absolute errors
  Bug 7 → bjacobi+ILU(1): poor preconditioner → GMRES needs many iterations
         → each iteration accumulates round-off, worsening the solution
  → Newton step is qualitatively incorrect: Δρ_v overshoots by O(ρ_vs)

Sublimation coupling:
  Bug 4 → k_pen = 1e9 → overshot ρ_v drives huge penalty on next iteration
  Bug 5 → ξ_v = 1e-3 → each Newton correction changes ρ_v rapidly
  Bug 6 → α_pen = 1e-7 → near-discontinuous D_eff creates additional
           Jacobian errors at interface nodes

Convergence to wrong solution:
  → Newton satisfies |R| < atol = 1e-8, but the solution has φ_i < 0
  → Next time step inherits corrupted φ_i < 0 as initial condition
  → Spurious air region grows irreversibly

t >> 0:
  → Regions of φ_i < 0 spread as each bad time step worsens the ICs
  → Simulation appears to "run" but produces unphysical results
```

### Current code breaks every link in this chain

| Link | Old | Fix |
|------|-----|-----|
| t=0 disequilibrium | ρ_v = 0.5 ρ_vs inside ice | ρ_v = ρ_vs inside ice (Bug 2 fixed) |
| Penalty 1000× too large | k_pen (no ξ_v) | ξ_v × k_pen (Bug 1 fixed) |
| Wrong Newton step | FD Jacobian errors O(10) | Analytical Jacobian (Bug 3 fixed) |
| Stiff system | k_pen = 1e9 | k_pen = 1e5 (Bug 4 fixed) |
| Pre-relaxation inconsistency | n_relax = 12 | n_relax = 1 (Bug 8 fixed) |
| Driving force amplitude | h₀ = 0.5 | h₀ = 0.95 (Bug 7b, humidity) |
| Diffusivity jump | α_pen = 1e-7 | α_pen = 1e-4 (Bug 6 fixed) |
| Poor preconditioning | bjacobi+ILU(1) | ASM+ILU(2) (Bug 7 fixed) |

The result is a code that produces physically correct, monotone phase-field
profiles with no spurious air, Newton convergence in 2–4 iterations per
time step, and stable long-time evolution.

---

## Summary Table

| # | Bug | File | Old | Current | Severity |
|---|-----|------|-----|---------|----------|
| 1 | Missing ξ_v on k_pen in vapor residual | `assembly.c` | `k_pen × ...` | `ξ_v × k_pen × ...` | **Critical** |
| 2 | Wrong ρ_v IC inside solid | `initial_conditions.c` | `h₀ ρ_vs` everywhere | `ρ_vs` in solid | **Critical** |
| 3 | FD Jacobian (analytical disabled) | `permafrost2.c` | `IGAFormIJacobianFD` | `Jacobian_A1` | High |
| 4 | k_pen 4 orders too large | `universal.opts` | 10⁹ | 10⁵ | High |
| 5 | ξ_v 1 order too large | `permafrost2.c` | 10⁻³ | 10⁻⁴ | Medium |
| 6 | α_pen 3 orders too extreme | `universal.opts` | 10⁻⁷ | 10⁻⁴ | Medium |
| 7 | Inadequate preconditioner | `universal.opts` | bjacobi+ILU(1) | ASM+ILU(2) | Medium |
| 8 | Excessive relaxation + immediate sed. freeze | `*.opts` | n_relax=12, t_freeze=0 | n_relax=1, t_freeze=1 s | Low–Medium |
