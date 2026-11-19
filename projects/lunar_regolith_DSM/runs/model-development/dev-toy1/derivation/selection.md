# M2 Selection: toy 1D heat conduction with a source (DRY RUN ONLY)

This is a plumbing test of the model line. It is not a result and must not be cited (spec_problem.md,
lines 3–5). It uses the approved decisions in `decisions.json`: C1 = (a), the unit case; C2 = any real q
in the code, tested with q > 0.

**Sources used here.** No paper is cited, and `citations.json` stays an empty list. The guard blocks
`library_search.py` at this stage, because the work order lists no commands. So each claim below is one of:
- a line of `spec_problem.md` (written "spec L<n>");
- a numbered step of my own derivation in this file (written "D<n>"), shown in enough detail to check by
  hand.

I could not run scratch Python at this stage either. Every number is a hand computation, and the
arithmetic is shown.

**The problem in dimensionless form** (spec L50–52):
∂θ/∂τ = ∂²θ/∂ξ² + 1 on 0 < ξ < 1, with θ(0, τ) = θ(1, τ) = 0 and θ(ξ, 0) = 0.
The steady state is θ_s(ξ) = ξ(1 − ξ)/2, whose peak is θ_s(1/2) = 1/8 (spec L46–47).

**Template items that do not apply.** The work order's template asks how grains enter and how the contact
line and the no-flux conditions are imposed. This toy has none of them. It is a homogeneous 1D rod with
no grains, no ice and no contact line (spec L108–111). Both ends are Dirichlet, so no flux condition is
imposed (spec L98–99). Each candidate says so in one line and treats the Dirichlet ends instead.

**Reference discretization** (proposed here; used only to count cost):
- N = 64 uniform intervals, so h = 1/64.
- Δτ = 10⁻³.
- τ_end = 1. This is above the steady-state time τ_99 ≈ 0.4666 (spec L42), with the margin the gate asked
  for (G:Scope, §2).

---

## Candidate: FD-BE (finite differences, backward Euler)

This is the method the brief already fixes (spec L129–130; brief lines 14–16 as cited there).

**Equations in outline.**
- Nodes ξ_j = jh for j = 0, …, N. The unknowns are θ_j^n at the interior nodes j = 1, …, N − 1.
- Second difference: (D_h θ)_j = (θ_{j−1} − 2θ_j + θ_{j+1})/h².
- Step: (θ^{n+1} − θ^n)/Δτ = D_h θ^{n+1} + 1, with θ_0 = θ_N = 0.
- Equivalently (I − Δτ D_h) θ^{n+1} = θ^n + Δτ·**1**. This is one tridiagonal solve per step.

**Grains, contact line, no-flux.** None apply (see above). The Dirichlet ends are imposed by removing
nodes 0 and N from the unknowns, with their value 0 entering the first and last rows.

**Facts derived for this candidate.**
- **D1. The second difference is exact on quadratics.** Taylor's theorem with remainder gives
  f(ξ ± h) = f ± h f′ + (h²/2) f″ ± (h³/6) f‴ + (h⁴/24) f⁗(η±).
  Adding the two expansions and subtracting 2f gives D_h f = f″ + (h²/24)[f⁗(η+) + f⁗(η−)].
  For a quadratic f, f⁗ = 0, so D_h f = f″ exactly.
- **D2. The discrete steady state equals the parabola at the nodes.**
  - By D1, θ_s restricted to the nodes satisfies D_h θ_s = θ_s″ = −1, with zero end values.
  - D_h with zero end values has eigenvectors sin(kπξ_j), k = 1, …, N − 1. Their eigenvalues are
    −λ_k, where λ_k = (4/h²) sin²(kπh/2). Each λ_k > 0, so D_h is nonsingular.
  - So the discrete steady state is unique, and it equals θ_s at every node, up to roundoff.
- **D3. Backward Euler is stable for every Δτ, and it does not oscillate.**
  - Write θ^n = θ_s + e^n. Then e^{n+1} = (I − Δτ D_h)⁻¹ e^n.
  - Mode k is multiplied by g_k = 1/(1 + Δτ λ_k) each step, and 0 < g_k < 1 for every Δτ > 0.
- **D4. Discrete maximum principle: 0 ≤ θ_j^n ≤ θ_s(ξ_j).**
  - I − Δτ D_h is an M-matrix. Its diagonal is 1 + 2Δτ/h² > 0, its off-diagonals are −Δτ/h² ≤ 0, and it is
    strictly diagonally dominant. So its inverse is entrywise nonnegative.
  - Lower bound: θ^{n+1} = (I − Δτ D_h)⁻¹(θ^n + Δτ) ≥ 0 by induction from θ^0 = 0.
  - Upper bound: θ_s − θ^n = (I − Δτ D_h)^{−n}(θ_s − θ^0). Since θ_s − θ^0 = θ_s ≥ 0, this is ≥ 0.
- **D5. Size of the spatial error in the transient.**
  - λ_k = (kπ)² sin²(z)/z² with z = kπh/2. Since sin²(z)/z² = 1 − z²/3 + O(z⁴), this gives
    λ_k = (kπ)²[1 − (kπh)²/12 + O(h⁴)].
  - For k = 1 and h = 1/64: (π/64)² = 0.0024096, and divided by 12 that is 2.008·10⁻⁴.
    So λ_1 ≈ 9.8696·(1 − 2.008·10⁻⁴) ≈ 9.8676.
- **D6. Size of the time error in the slowest mode.**
  - The exact amplitude of mode 1 is exp(−λ_1 τ). Backward Euler gives (1 + λ_1 Δτ)^{−τ/Δτ}.
  - The log of their ratio is (τ/Δτ)[λ_1Δτ − ln(1 + λ_1Δτ)] = τ λ_1² Δτ/2 + O(Δτ²).
  - The mode-1 coefficient of the transient is b_1 = 4/π³ ≈ 0.12901 (D9 below).
  - So the absolute error is about b_1 (τλ_1²Δτ/2) e^{−λ_1τ}. Its maximum over τ is at τ = 1/λ_1, where
    it equals b_1 (λ_1Δτ/2)/e.
  - For Δτ = 10⁻³: 0.12901 × 0.0049338 / 2.71828 ≈ 2.34·10⁻⁴. That is 0.19 % of the peak 1/8.

**What can go wrong.**
1. **First-order time error in the transient (D6).** A transient test at a fixed τ needs Δτ small enough.
   - Δτ = 10⁻³ gives about 2.3·10⁻⁴ in the slowest mode.
   - **Not estimated:** at early times the higher modes and the corners carry more error. At τ = 0 the
     source gives ∂θ/∂τ = 1 at the end nodes while the boundary condition holds θ at 0, so the initial and
     boundary data are incompatible at the corners. A transient test should therefore check the observed
     convergence order at τ of order 0.1 or later, not the error at the first few steps.
2. **First-order flux at the ends.** The one-sided difference (θ_{N−1} − θ_N)/h at the steady state gives
   θ_s(1 − h)/h = (1 − h)/2, an error of h/2.
   - The fix is a half-cell balance: add the source in the half cell, (h/2)·1, to the one-sided difference.
     This gives (1 − h)/2 + h/2 = 1/2, which matches the exact dimensionless flux of 1/2 (spec L100).
3. **Steady-state test cannot catch spatial error (D2).** The parabola is reproduced at the nodes
   whatever h is. So the steady test checks the plumbing, but not the order of the spatial scheme. Only the
   transient test (D5) checks that order.
4. **Variant risk: Crank–Nicolson.** It is not proposed. Its mode factor (1 − Δτλ_k/2)/(1 + Δτλ_k/2) tends
   to −1 as Δτλ_k grows. Combined with the corner incompatibility in item 1, high modes would flip sign
   from step to step. D3 shows that backward Euler has no such behavior.

**Cost for the reference case.**
- N − 1 = 63 unknowns.
- One tridiagonal solve per step, which is O(N) work.
- τ_end/Δτ = 1000 steps, so about 63,000 unknown-updates in total. This is negligible.

**Fit to this solver.** The toy's solver is a NumPy module, `toy/heat1d.py` (spec L109). The brief
specifies FD with backward Euler for it (spec L129–130), and this candidate matches that exactly. It does
not use the PetIGA C solver, which is out of scope for the toy (spec L109).

**Verification.**
- **V1. Steady state.** θ_j(τ_end) vs θ_s(ξ_j) to roundoff level, by D2. Allow for the residual transient at
  τ_end = 1: about b_1 e^{−λ_1} ≈ 0.129 × 5.2·10⁻⁵ ≈ 6.7·10⁻⁶. If roundoff level is wanted, also run a
  direct steady solve.
- **V2. Transient.** Compare against the eigenfunction series (D9) at a fixed τ, for example τ = 0.1.
  Check the observed order 2 in h with Δτ small, and order 1 in Δτ with h small.
- **V3. Bounds.** 0 ≤ θ ≤ 1/8 at every node and step, by D4.
- **V4. End flux.** The half-cell flux at the steady state equals 1/2 (item 2).

---

## Candidate: FEM-Galerkin (Galerkin finite elements: linear P1 or quadratic B-splines, backward Euler)

This candidate is included because it is the discretization family of the PetIGA solver. It would depart
from the brief's FD choice (spec L129–130).

**Equations in outline.**
- Find θ_h(τ) in V_h with θ_h = 0 at the ends, such that (∂_τθ_h, v) + (θ_h′, v′) = (1, v) for every v in
  V_h.
- With backward Euler, each step solves (M + Δτ K) θ^{n+1} = M θ^n + Δτ F.
- For P1 elements: M has diagonal 2h/3 and off-diagonals h/6, K has diagonal 2/h and off-diagonals −1/h,
  and F_i = h.

**Grains, contact line, no-flux.** None apply. The Dirichlet ends are strong, by removing the end basis
functions. A no-flux end, if one were ever wanted, would be natural: it would need no term.

**Facts derived for this candidate.**
- **D7. P1 nodal values at the steady state are exact.**
  - Galerkin orthogonality gives (θ_s − θ_h)′ ⊥ v′ for every v in V_h.
  - The Green's function G_j of −d²/dξ² with a source at ξ_j is piecewise linear, with its kink at ξ_j. So
    G_j lies in V_h.
  - Then (θ_s − θ_h)(ξ_j) = ((θ_s − θ_h)′, G_j′) = 0.
- **D7b. Quadratic B-splines reproduce the steady state exactly everywhere.** θ_s is quadratic and so lies
  in V_h. The Galerkin solution is the energy-norm best approximation, so θ_h = θ_s.
- **D8. P1 with a consistent mass matrix: the maximum-principle argument needs Δτ ≥ h²/6.**
  - The off-diagonals of M + ΔτK are h/6 − Δτ/h. They are ≤ 0 only if Δτ ≥ h²/6.
  - When that holds, the matrix is strictly diagonally dominant: 2h/3 + 2Δτ/h > 2(Δτ/h − h/6). So it is an
    M-matrix, and the argument of D4 carries over.
  - For Δτ < h²/6 the argument fails. **I have not shown that negative values actually occur.** I have
    shown only that the sufficient condition is lost.
  - For h = 1/64 the bound is h²/6 ≈ 4.07·10⁻⁵, so the reference Δτ = 10⁻³ satisfies it.
- **D8b. P1 with a lumped mass matrix is the same as FD-BE.**
  - Lumping gives M_L = hI. The step becomes (hI + ΔτK)θ = hθ^n + Δτh·**1**.
  - Dividing by h and using K/h = −D_h gives exactly the FD-BE step.

**What can go wrong.**
1. Small time steps with a consistent mass matrix lose the guaranteed bounds (D8).
2. The departure from the brief (spec L129–130) would reopen a decision the brief already made.
3. With quadratic B-splines (D7b), the steady test is exact by construction. Like FD-BE item 3, it checks
   the plumbing but not the convergence order.
4. The time error has the same first order and the same corner incompatibility as FD-BE item 1.

**Cost for the reference case.**
- P1: 63 unknowns, with a tridiagonal M + ΔτK.
- Quadratic C¹ B-splines on 64 elements: 64 + 2 = 66 basis functions, minus the 2 at the ends, gives
  64 unknowns. The matrix is pentadiagonal.
- 1000 steps, plus quadrature assembly. Still negligible.

**Fit to this solver.** It is a closer analogue of the PetIGA Galerkin path. But the toy is a NumPy module
by brief, so the extra assembly code tests nothing that the line needs from this dry run.

**Verification.** The same as V1–V4. V1 is nodally exact for P1 (D7) and exact everywhere for quadratic
splines (D7b). V3 holds only if Δτ ≥ h²/6 (D8).

---

## Candidate: Series (eigenfunction series, closed form)

**Equations in outline.**
- **D9.** Write θ = θ_s − w. Then w solves the homogeneous heat equation with w(ξ, 0) = θ_s.
- The sine coefficients of θ_s are b_k = ∫₀¹ ξ(1 − ξ) sin(kπξ) dξ. The pieces are:
  - ∫₀¹ ξ sin(kπξ) dξ = −(−1)^k/(kπ);
  - ∫₀¹ ξ² sin(kπξ) dξ = −(−1)^k/(kπ) + 2((−1)^k − 1)/(kπ)³.
  - So b_k = 2(1 − (−1)^k)/(kπ)³. This is 4/(kπ)³ for odd k and 0 for even k.
- So θ(ξ, τ) = ξ(1 − ξ)/2 − Σ_{k odd} [4/(kπ)³] sin(kπξ) e^{−k²π²τ}.
- **Check at ξ = 1/2, τ = 0:** Σ_{m≥0} (−1)^m/(2m + 1)³ = π³/32. So the sum is (4/π³)(π³/32) = 1/8 = θ_s(1/2),
  and θ(1/2, 0) = 0, as required.

**Grains, contact line, no-flux.** None apply. The Dirichlet ends hold by construction, because
sin(kπ·0) = sin(kπ) = 0.

**What can go wrong.**
1. **Truncation at τ = 0.** Keeping the odd k ≤ K leaves a tail of about Σ_{k>K, odd} 4/(kπ)³
   ≈ (4/π³)·1/(4K²) = 1/(π³K²). For K = 11 this is about 2.7·10⁻⁴. For τ > 0 the tail falls off like
   e^{−K²π²τ}.
2. **It is not a model of the general problem.** It needs constant coefficients, a uniform source and
   Dirichlet ends. It has no time stepper and no linear solve, so it does not exercise the plumbing that the
   dry run is for (spec L3–5).

**Cost for the reference case.** No unknowns to solve for. Evaluating K/2 terms at each output point,
with K = 11, is 6 terms.

**Fit to this solver.** It is not a solver. Its role is as the reference solution in V2.

**Verification.** The check at ξ = 1/2, τ = 0 above. The limit τ → ∞ gives θ_s. Truncation convergence in K.

---

## Comparison

| Criterion | FD-BE | FEM-Galerkin | Series |
|---|---|---|---|
| Fidelity to the stated problem | Full transient and steady. Steady is nodally exact (D2). Error is O(h²) + O(Δτ) (D5, D6). | Same orders. Steady is nodally exact (P1, D7) or exact (quadratic splines, D7b). | Exact up to truncation (D9). Valid only for constant coefficients. |
| Geometry flexibility | Uniform 1D grid. Enough for the toy (spec L108). | Nonuniform meshes are easy. Closer to PetIGA. | None beyond the 1D interval with Dirichlet ends. |
| Cost (reference case) | 63 unknowns, tridiagonal, 1000 steps | 63 (P1) or 64 (quadratic) unknowns, plus assembly, 1000 steps | Closed form, 6 terms |
| Fit to the solver (NumPy `toy/heat1d.py`, spec L109) | Exactly what the brief specifies (spec L129–130) | Departs from the brief | Not a solver. Used as the test oracle. |
| Robustness | Unconditionally stable. Bounds hold for every Δτ (D3, D4). | Bounds guaranteed only for Δτ ≥ h²/6 with consistent mass (D8) | Not applicable |
| Buildable and verifiable by December | Yes. About 50 lines. V1–V4. | Yes, but more code for no gain in the toy (D8b) | Yes. It is the oracle for V2. |

## Recommendation

**Choose FD-BE, and use Series as the reference solution in the transient test.**

**Reason.**
- The brief already fixes FD with backward Euler in NumPy (spec L129–130).
- FD-BE is unconditionally stable and keeps 0 ≤ θ ≤ 1/8 for every Δτ (D3, D4).
- Its steady state matches the parabola at the nodes, which is the brief's central test (D2).
- In 1D, lumped P1 FEM is the same scheme (D8b), so FEM-Galerkin adds code without adding a check.
- The one thing the steady test cannot catch, the spatial order (item 3), is covered by comparing against
  the closed-form series (V2).

**Proposed model name:** `heat1d-fd-be`, for the branch `feat/heat1d-fd-be`.

**Points for the user at G:Selection.**
- Confirm FD-BE, or reopen the brief's discretization choice in favor of FEM-Galerkin.
- Confirm the reference discretization: N = 64, Δτ = 10⁻³, τ_end = 1.
