# M4 Discretization spec: `heat1d-fd-be` (DEV TOY DRY RUN ONLY)

This is a plumbing test of the model line. It is not a result and must not be cited (problem.md, lines 3–4).

**What this spec covers.** It is the discretization of the approved strong form S29 (ledger.json; derivation.md §4).
- The scheme is the one decided at G:Selection (decisions.json, 2026-11-16T14:21:08): second-order finite differences
  in space, backward Euler (BE) in time.
- Reference discretization: N = 64, Δτ = 10⁻³, τ_end = 1.
- V1 (c): run both the τ_end comparison and a direct steady solve.
- The implementation is the NumPy module `toy/heat1d.py` in the feature worktree, with
  `solve(n, dt, t_end, params) -> (x, T)` (problem.md, lines 14–16). It is not the C solver (spec_problem.md, line 109).

**Where the C-solver items of the work order go.** The work order is written for the PetIGA solver. Each of its items is
handled as follows:
- The anchor files (`src/assembly.c`, `src/lunar_main.c`, `src/monitoring.c`, `inputs/solver.opts`) and the PETSc
  TSALPHA source were read. What I found is reported in §"Time scales and time step" → *Anchor solver: TSALPHA and
  TSAdapt*, with file and line.
- The "infra" layer (a manufactured-solution forcing option and a full-Jacobian test) is specified for the toy module,
  in the same form it would take in the C solver: the residual F(τ, U, U̇) and the Jacobian shift·∂F/∂U̇ + ∂F/∂U.
- Grains, vapor, phase field and the wall do not exist in the toy (S8, S40). Their sections say "not applicable" in one
  line each.

**Sources.** No paper is cited. Each step below is one of these:
- a ledger step (S#);
- a fact in selection.md (D#);
- a line of a named file;
- a step of this spec, numbered M4-#, with its premises and the operation named.

Numbers come from my scratch scripts: `derivation/scratch/m4_checks.py` and `derivation/scratch/m4_mms_copy.py` (the
second is a copy of `mms/mms_sources.py`; see §"Manufactured solutions"). Both were run with
`venv_workflows/bin/python -I`. They are my own checks, not evidence for anyone else.

### New symbols (discretization only)

All other symbols are those of ledger.json → `symbols`.

| Symbol | Meaning | Units |
|---|---|---|
| N | number of uniform intervals on [0, 1] (the `n` argument of `solve`) | 1 |
| h | dimensionless mesh size, 1/N | 1 |
| ξ_j | node j·h, j = 0, …, N | 1 |
| Δτ | dimensionless time step, κ·dt/L² | 1 |
| τ^n | n·Δτ | 1 |
| U_j^n | θ′ at node j and time τ^n. Unknowns: j = 1, …, N − 1. U_0 = U_N = 1 are data. | 1 |
| U̇ | time derivative of U, which the residual takes as an argument | 1 |
| g_L, g_R | Dirichlet data for θ′ at ξ = 0 and ξ = 1. Both are 1 (S29). | 1 |
| f(ξ, τ) | manufactured-solution forcing. It is ≡ 0 unless the MMS option is on. | 1 |
| F_j(τ, U, U̇) | residual function of row j, in "implicit-function" form | 1 |
| R_j(U^{n+1}) | BE residual of row j: F_j evaluated at U̇ = (U^{n+1} − U^n)/Δτ | 1 |
| shift | ∂U̇/∂U^{n+1} = 1/Δτ for BE; 0 for the steady solve | 1 |
| J | Jacobian dR/dU^{n+1} = shift·∂F/∂U̇ + ∂F/∂U, an (N−1)×(N−1) tridiagonal matrix | 1 |
| D_h | second-difference operator, (D_hU)_j = (U_{j−1} − 2U_j + U_{j+1})/h² (selection.md, line 39) | 1 |
| λ_ℓ | eigenvalues of −D_h with zero end values: (4/h²) sin²(ℓπh/2) (selection.md D2, line 54) | 1 |
| b_1 | mode-1 sine coefficient of θ_s, 4/π³. This is β_1 of S53. | 1 |
| ε_θ′ | allowed absolute time-discretization error in θ′ (the input of the dt rule) | 1 |
| θ′_m, A, B | manufactured solution and its two amplitudes | 1 |

**Why θ′ and not θ.** The unknown is θ′ = T/T_0, from S27 and S29, not θ = (T − T_0)/ΔT_q.
- θ is defined only for q ≠ 0 (S20). Decision C2 allows any real q, including q = 0.
- θ′ needs only T_0 > 0, which S7 gives.
- By S30, θ′ = 1 + Qθ. So every θ-statement carries over with a factor Q.

---

## Weak form

The steps go from S29 to the weak form that the finite-difference rows come from.

**M4-1** (from S29; operation: add a forcing term and move every term to one side). With the forcing f added, S29 reads

$$F(\xi,\tau,\theta',\partial_\tau\theta') \equiv \partial_\tau\theta' - \partial_\xi^2\theta' - Q - f(\xi,\tau) = 0,\quad 0<\xi<1,\ \tau>0;\qquad \theta'(0,\tau)=\theta'(1,\tau)=1,\quad \theta'(\xi,0)=1 \tag{M4-1}$$

In production f ≡ 0, and M4-1 is S29 exactly. Units: every term is dimensionless. The dimensional row is (k T_0/L²) × M4-1 in
W m⁻³. This factor comes from undoing S28 → S29: S28 was divided by κT_0/L², and S17 is S15 divided by ρc_p, and
ρc_p · κT_0/L² = kT_0/L² (S16).

**M4-2** (from M4-1 and S6; operation: multiply by a test function, integrate over (0, 1), integrate the ∂_ξ² term by parts).
Let w ∈ H¹₀(0, 1), so w(0) = w(1) = 0.
- Integration by parts gives −∫₀¹ ∂_ξ²θ′ w dξ = ∫₀¹ ∂_ξθ′ ∂_ξw dξ − [∂_ξθ′ w]₀¹.
- The bracket is zero because w vanishes at both ends.
- The integration by parts is valid for τ > 0, because θ′ is C² up to the ends (S6 through the affine map S27).

The result: find θ′(·, τ) with θ′ − 1 ∈ H¹₀(0, 1) such that

$$\int_0^1\Big(\partial_\tau\theta'\,w + \partial_\xi\theta'\,\partial_\xi w - (Q+f)\,w\Big)\,d\xi = 0\qquad\forall\,w\in H^1_0(0,1) \tag{M4-2}$$

The Dirichlet ends are essential conditions. No boundary integral remains, so there is no natural (flux) term. This
matches S39: the end fluxes are outputs, not imposed.

**M4-3** (from M4-2; operation: Galerkin with P1 hat functions, lumped mass, nodal quadrature of the load; then divide row
i by h). This step only shows where the FD rows come from. The FD row is then taken as the definition.
- **Trial space.** θ′_h = g_L φ_0 + g_R φ_N + Σ_{j=1}^{N−1} U_j φ_j. The test functions are w = φ_i, i = 1, …, N − 1.
- **Stiffness.** ∫φ_i′φ_i′ = 2/h and ∫φ_i′φ_{i±1}′ = −1/h.
- **Lumped mass.** h on the diagonal.
- **Load.** Nodal (trapezoid) quadrature gives ∫(Q + f)φ_i ≈ h(Q + f(ξ_i)). This is exact for the constant Q. For f the
  error is O(h²), which matches the scheme's order.
- **Divide by h.** This gives (selection.md D8b, lines 145–147):

$$\dot U_i - \frac{U_{i-1}-2U_i+U_{i+1}}{h^2} - Q - f(\xi_i,\tau) = 0,\qquad i=1,\dots,N-1,\quad U_0=g_L,\ U_N=g_R \tag{M4-3}$$

Consistency in space: D_h u = u″ + (h²/24)[u⁗(η₊) + u⁗(η₋)] (selection.md D1, lines 47–50). So the row has O(h²)
truncation error, and it is exact on quadratics.

## Discretization

**M4-4** (from M4-3; operation: backward Euler). Evaluate every term at τ^{n+1} and replace U̇ by (U^{n+1} − U^n)/Δτ:

$$\frac{U_j^{n+1}-U_j^n}{\Delta\tau} - \frac{U_{j-1}^{n+1}-2U_j^{n+1}+U_{j+1}^{n+1}}{h^2} - Q - f(\xi_j,\tau^{n+1}) = 0,\quad j=1,\dots,N-1 \tag{M4-4}$$

- **Local truncation error.** O(Δτ) from the backward difference and O(h²) from D_h.
- **Stability.** Unconditional, with no oscillation: mode ℓ is multiplied by 1/(1 + Δτλ_ℓ) ∈ (0, 1) (selection.md D3,
  lines 56–58).
- **Discrete maximum principle.** It holds for every Δτ (D4, lines 59–63).
- **Data.** The initial data are U_j^0 = 1 (S27 applied to S7). The end data are U_0^n = U_N^n = 1 for every n.

**Grid and scales used by `solve(n, dt, t_end, params)`.** `params` holds L, k, rho, cp, q, T0.

| Quantity | Formula | Source |
|---|---|---|
| N | `n` | brief signature |
| h | 1/N | |
| κ | k/(ρc_p) | S16 |
| Q | qL²/(kT_0) | S27 |
| Δτ | κ·dt/L² | S20 (τ = κt/L²) |
| number of steps | t_end/dt | must be an integer to within 10⁻¹⁰ relative, else `ValueError` |
| outputs | x_j = L·jh and T_j = T_0·U_j for j = 0, …, N, including both Dirichlet end values | S27 inverted |

I chose to require an integer number of steps rather than silently shrinking dt. The reason is that the M5 order test in
Δτ needs the step to be exactly what the caller passed. This is a proposal (question 2).

**Direct steady solve (V1 (c)).**

**M4-5** (from M4-3; operation: set U̇ = 0 and τ-independent data, with f ≡ 0). The steady rows are

$$-\frac{U_{j-1}-2U_j+U_{j+1}}{h^2} - Q = 0,\qquad U_0=U_N=1 \tag{M4-5}$$

- **Unique solution.** −D_h is nonsingular, because every λ_ℓ > 0 (selection.md D2, line 54).
- **The solution is exact at the nodes.** U_j = 1 + Qξ_j(1 − ξ_j)/2, which is S43 through S30, because D_h is exact on
  quadratics (D1, D2).
- **Scratch check.** The steady residual evaluated on these nodal values is exactly 0.0 (`m4_checks.py` item 9).

## Residual

**M4-6** (from M4-4; operation: definition). The residual has two parts.

The implicit-function residual for rows j = 1, …, N − 1, with U_0 = g_L and U_N = g_R (both 1), is

$$F_j(\tau,U,\dot U) = \underbrace{\dot U_j}_{\text{storage}} \;\underbrace{-\,\frac{U_{j-1}-2U_j+U_{j+1}}{h^2}}_{\text{conduction}} \;\underbrace{-\,Q}_{\text{source}} \;\underbrace{-\,f(\xi_j,\tau)}_{\text{MMS forcing (option)}} \tag{M4-6}$$

The BE residual that the step solves is

$$R_j(U^{n+1}) = F_j\Big(\tau^{n+1},\,U^{n+1},\,\frac{U^{n+1}-U^n}{\Delta\tau}\Big) \tag{M4-6b}$$

Rows 1 and N − 1 read the Dirichlet data inside the conduction term:
- row 1 contains −g_L/h²;
- row N − 1 contains −g_R/h².

These are listed as their own term (T3) in the table and in code_map.json, because a wrong boundary value shows up only
there.

| ID | Term in row j | Rows | Value at the reference state |
|---|---|---|---|
| T1 | storage (U_j − U_j^n)/Δτ | 1 … N−1 | 0 at the start of a step (guess U = U^n) |
| T2 | conduction, interior neighbours: −(U_{j−1} − 2U_j + U_{j+1})/h², using only unknowns | 1 … N−1 (neighbours that are unknowns) | O(1) on a smooth U |
| T3 | conduction, Dirichlet data: −g_L/h² in row 1, −g_R/h² in row N−1 | 1 and N−1 | with T2 gives O(1) |
| T4 | source −Q | 1 … N−1 | −1 (C1) |
| T5 | MMS forcing −f(ξ_j, τ^{n+1}) | 1 … N−1 | 0 when the option is off |

Every row is O(1) for the reference case, as the role brief asks:
- At the first step, with guess U = U^0 = 1, row j equals −Q − f = −1 (C1, option off).
- T2 and T3 separately reach 2/h² = 8192 for N = 64. But they cancel to O(1) on any smooth U, because D_h of a smooth
  function is O(1) (D1).

For general parameters the rows scale with Q, so a small |Q| gives small rows. This is harmless, because the solve is
direct (§Solver settings) and no absolute stopping tolerance is used.

The steady residual is M4-6 with U̇ = 0 and f = 0, which is M4-5.

## Jacobian

**M4-7** (from M4-6 and M4-6b; operation: chain rule with ∂U̇/∂U^{n+1} = 1/Δτ ≡ shift). Then J = shift·∂F/∂U̇ + ∂F/∂U.
This is the same convention as the C solver's IJacobian: `J[a][i][b][j] = dR/du + shift * dR/du_t` (src/assembly.c,
line 206), and TSALPHA passes shift_V = α_m/(α_f γ Δt) (alpha1.c, lines 65 and 386). The entries are all constant, so J
does not depend on U:

| Term | ∂R_j/∂U_{j−1} | ∂R_j/∂U_j | ∂R_j/∂U_{j+1} | Where it is a ∂F/∂U̇ entry |
|---|---|---|---|---|
| T1 storage | 0 | **shift** (= 1/Δτ; 0 for steady) | 0 | yes: ∂F_j/∂U̇_j = 1, times shift |
| T2 conduction (interior) | **−1/h²** (j ≥ 2) | **+2/h²** | **−1/h²** (j ≤ N−2) | no |
| T3 Dirichlet data | none (g_L, g_R are data, not unknowns) | 0 | 0 | no |
| T4 source | 0 | 0 | 0 | no |
| T5 MMS forcing | 0 | 0 | 0 (f depends on ξ and τ only) | no |

The resulting matrix, stored as three diagonals of length N − 1, N − 1 and N − 1 (end entries of the off-diagonals
unused), is:
- diagonal: shift + 2/h²;
- sub-diagonal and super-diagonal: −1/h².

Properties:
- **Shape.** J is symmetric and tridiagonal, as the brief asks (problem.md, line 24).
- **shift > 0.** J is strictly diagonally dominant, with diagonal minus off-diagonal sum ≥ shift > 0.
- **shift = 0 (steady).** J = −D_h is irreducibly diagonally dominant and nonsingular (D2).
- **Either case.** J is a nonsingular M-matrix (D4), so the Thomas algorithm needs no pivoting and is stable.

**Scratch check** (`m4_checks.py` items 5 and 8). J·v was compared with a central difference
[F(U + εv, U̇ + shift·εv) − F(U − εv, U̇ − shift·εv)]/(2ε), with ε = 10⁻⁶, at a random state:
- Relative error 7e-11 to 1.0e-10 for shift ∈ {0, 10, 10³}, and 8.8e-12 for shift = 10⁵.
- Split by part: ∂F/∂U gives ≤ 1.0e-10. ∂F/∂U̇ = I gives ≤ 6.3e-8, which is roundoff from rows of size about 2/h² with
  ε = 10⁻⁶.
- **Sensitivity.** Dropping the shift term gives a relative error of 9.6e-2, so the check catches that mutation.
- **Why the check is at roundoff.** F is affine in (U, U̇), so the central difference has no truncation error. Only
  roundoff is left.

## Auxiliary functions

- **Indicator ψ and its gradient.** Not applicable: there are no grains (S8). No function.
- **Initial conditions.** U_j^0 = 1 for j = 1, …, N − 1 (S7 through S27). The MMS option does not change the initial
  condition, because θ′_m(ξ, 0) = 1 (§Manufactured solutions).
- **Extension into grains.** Not applicable (S8).
- **Output masks.** None. `solve` returns x_j = L·jh and T_j = T_0·U_j at all N + 1 nodes, including the two Dirichlet end
  values.
- **End fluxes** (V4; selection.md item 2, lines 84–87). `end_fluxes(U, U_prev, dtau, h, Q, f_end)` uses the half-cell
  energy balance in θ′ form.
  - At ξ = 0: F′_0 = (U_1 − U_0)/h + (h/2)(Q + f(0, τ)) − (h/2)(U_0 − U_0^prev)/Δτ.
  - At ξ = 1: the mirror formula with U_{N−1} and U_N.
  - The last term is zero because U_0 is fixed.
  - Dimensional flux: F_0 = (kT_0/L)·F′_0 W m⁻². At steady state (kT_0/L)(Q/2) = qL/2, which is S45.
  - Scratch: the half-cell flux on the steady nodal values is 0.500000000000000 = Q/2 (`m4_checks.py` item 9).
- **Scales.** `_scales(params)` returns κ, Q, t_c = L²/κ and the factor kT_0/L².
- **Time-step helper.** `dt_for_tolerance(eps_theta_prime, n, t_end, params)` implements the dt rule below.

## Time scales and time step

### A priori scales (reference case C1, N = 64)

All values are from `m4_checks.py` items 1–2. In C1, τ and t coincide in value (t_c = 1 s).

| Scale | Formula | Value |
|---|---|---|
| Domain diffusion time t_c | L²/κ (S20) | 1 s |
| Slowest mode, continuous | L²/(π²κ) | 0.101321 s |
| Slowest mode, discrete | 1/λ_1, with λ_1 = 4N² sin²(π/(2N)) | λ_1 = 9.867623 (π² = 9.869604) |
| Mesh diffusion time | (L/N)²/κ | 2.44141e-4 s |
| Fastest discrete mode | 1/λ_max, with λ_max = λ_{N−1} | 6.10719e-5 (λ_max = 16374.13) |
| Stiffness ratio λ_max/λ_1 | | 1659 |
| Interface relaxation, kinetic, vapor diffusion | | not applicable: no interface, no vapor (S8) |

Consequences:
- An explicit method would need Δτ ≤ 2/λ_max = 1.22e-4.
- BE has no stability limit (D3). So the step is chosen for accuracy only.

### The dt rule

**M4-8** (from S53, S54 and D3; operation: compare the per-step BE amplification of mode 1 with the exact decay, then
expand to first order in Δτ).
- **Mode 1 of the transient.** Its amplitude is b_1 e^{−λ_1τ} exactly in time (S54 with the discrete λ_1), and
  b_1(1 + λ_1Δτ)^{−τ/Δτ} under BE.
- **Their ratio.** It is exp[(τ/Δτ)(λ_1Δτ − ln(1 + λ_1Δτ))] = 1 + τλ_1²Δτ/2 + O(Δτ²).
- **Error and its maximum.** So the time error in θ is about b_1(τλ_1²Δτ/2)e^{−λ_1τ}. Its maximum over τ is at
  τ = 1/λ_1, where it equals b_1λ_1Δτ/(2e) (selection.md D6, lines 69–75).
- **In θ′.** Multiply by |Q| (S30).

Setting this error equal to a tolerance ε_θ′ gives the rule.

**Inputs:**
- ε_θ′, the allowed absolute time error in θ′;
- the parameters L, k, ρ, c_p, q, T_0 (for κ and Q);
- N;
- t_end.

**Formula:**

$$\lambda_1 = 4N^2\sin^2\!\frac{\pi}{2N},\qquad \Delta\tau_{\max} = \frac{2e\,\varepsilon_{\theta'}}{|Q|\,b_1\,\lambda_1}\ \ (Q\neq0;\ \ \Delta\tau_{\max}=\tau_{end}\ \text{if}\ Q=0),\qquad b_1=\frac{4}{\pi^3}$$
$$n_{steps} = \Big\lceil \frac{\tau_{end}}{\Delta\tau_{\max}} \Big\rceil,\qquad \Delta\tau = \frac{\tau_{end}}{n_{steps}},\qquad dt = \Delta\tau\,\frac{L^2}{\kappa},\qquad \tau_{end} = \frac{\kappa\,t_{end}}{L^2} \tag{M4-8}$$

**Reference values** (`m4_checks.py` item 3):
- Δτ = 10⁻³ corresponds to ε_θ′ = 2.3415e-4 with Q = 1. That is 0.19 % of the peak rise 1/8.
- ε_θ′ = 10⁻⁴ gives Δτ_max = 4.27e-4.

**What supports the rule, and what does not.**
- **Measured, my scratch only** (`m4_checks.py` item 4). The BE error alone, defined as the fully discrete minus the
  semi-discrete solution for N = 64 with all modes, has a maximum over τ that is 0.982, 0.991, 0.995 and 0.997 times the
  mode-1 prediction, for Δτ = 4, 2, 1 and 0.5 × 10⁻³. The maxima sit at τ = 0.102 to 0.104, which is close to
  1/λ_1 = 0.1013.
- **Not shown analytically.** I have not proved that the higher modes never add to the maximum. Each mode ℓ peaks at
  τ = 1/λ_ℓ with first-order size (4/(ℓπ))Δτ/(2e). The sum of these sizes over ℓ diverges, and only the fact that they
  peak at different times keeps them from adding. The rule is a first-order estimate checked numerically, not a bound.
- **The corners.** The corner incompatibility (selection.md item 1) affects only the first few steps. The rule does not
  cover those steps, and the M5 tests should not measure there (selection.md, lines 81–83).
- **Spatial error is separate.** For mode 1 it is λ_1/π² − 1 = −2.0e-4 relative (D5).

**For the M5 tests (D-d).** At τ_end = 1, the fully discrete (N = 64, Δτ = 10⁻³) maximum of |θ − θ_s| is 7.0172e-6
(`m4_checks.py` item 9), slightly above the continuous 6.67e-6 (S57).
- A τ_end-versus-S44 tolerance at 1e-5 in θ would have a margin of only 1.4×.
- The direct steady solve (M4-5) has no such floor.

### Anchor solver: TSALPHA and TSAdapt (what I found)

This section is information for the later C-solver work. The toy does not use TSALPHA.

1. **The solver uses TSALPHA with ρ∞ = 0.5.**
   - It is set in code: `TSSetType(ts, TSALPHA)` and `TSAlphaSetRadius(ts, 0.5)` (src/lunar_main.c, lines 1357–1358).
     It is also set by option: `-ts_alpha_radius 0.5` (inputs/solver.opts, line 92).
   - With ρ∞ = 0.5: α_m = 0.5(3 − ρ∞)/(1 + ρ∞) = 5/6, α_f = 1/(1 + ρ∞) = 2/3, and γ = 1/2 + α_m − α_f = 2/3
     (alpha1.c, lines 480–482).
2. **TSALPHA's default adaptor is TSADAPTNONE.**
   - `ts->default_adapt_type = TSADAPTNONE` (alpha1.c, line 546).
   - With NONE, TSSetUp_Alpha does not allocate `vec_sol_prev` or `vec_lte_work` (alpha1.c, lines 433–437). So
     TSEvaluateWLTE_Alpha returns wlte = −1 (alpha1.c, lines 311–318): no local-truncation-error estimate is ever made.
   - TSAdaptChoose_None always accepts and keeps h (adaptnone.c, lines 6–8).
   - solver.opts sets no `-ts_adapt_type`.
   - So **error-based step control is off**. dt is set only by the patched heuristic below and by the solver's own
     monitors.
3. **The step-size control is a local patch inside TSStep_Alpha** (the lines marked `//AM`, alpha1.c, lines 188–294). It
   reads the fields `ts->adap`, `ts->NRmin`, `ts->NRmax`, `ts->factor`, `ts->dtmin`, `ts->dtmax`, `ts->maxdt_dyn` and
   `ts->snes->prev_dt_red`, which lunar_main.c sets directly (lines 1372–1379). These fields are not in stock PETSc. I could
   not read the patched header (`include/petsc/private/tsimpl.h` is outside this stage's read roots), so their
   declarations and defaults were not checked.
   - **NR, the Newton count of a step.** NR = (SNES iterations of this step), alpha1.c line 240. On a restart (the first
     step) it is NR = 1 + its/3, line 219, which averages over the three solves of TSAlpha_Restart.
   - **Shrink.** If NR ≥ NRmax, the step is **kept**:
     - `ptime` advances (line 253), and `vec_sol` already holds X1 (line 243). `th->status` was already
       `TS_STEP_COMPLETE` (line 245), so the `while` at line 198 exits.
     - dt is divided by `factor` (line 254).
     - But the code jumps to `reject_step` (line 260), which increments `ts->reject` (line 291).
     - So `TSGetStepRejections`, and the "rejections" column of solver_evo.dat (src/monitoring.c, lines 313–316), count
       these accepted-and-shrunk steps as rejections.
   - **Grow.** If NR ≤ NRmin, and the step had no rejection in this call (`rep_ts == 0`), dt is multiplied by `factor`,
     unless that would exceed dtmax (lines 278–285).
   - **Stage failure** (SNES diverged or domain error). TSAdaptCheckStage rejects the stage and multiplies dt by
     `scale_solve_failed`, which defaults to 0.25 (tsadapt.c, lines 1021–1025 and 1075). When `adap` is on, the patch then
     **overwrites** this with dt_old/factor (alpha1.c, lines 229–235 and 209–216). So the 0.25 never takes effect under
     `-adap 1`.
     - `ts->max_snes_failures = -1` (lunar_main.c, line 1379) means SNES failures never stop the run (tsadapt.c,
       line 995). Only `max_reject` does (alpha1.c, lines 295–297).
   - **Two separate dt limits.** The patch uses `ts->dtmax` and `ts->dtmin`, not the TSAdapt limits `adapt->dt_min` and
     `adapt->dt_max` (defaults 1e-20 and 1e+20, tsadapt.c, lines 1071–1072). So `-ts_adapt_dt_max` would not cap the
     patch, and `-dtmax` does not cap TSAdapt.
   - **If someone sets `-ts_adapt_type basic`.** The LTE vectors would be allocated, and the basic controller's `next_h`
     would become `next_time_step` (line 244). The patch would then multiply or divide that by `factor`. The two
     controllers would compound. I did not test this.
4. **Final time.** With `TS_EXACTFINALTIME_MATCHSTEP` (lunar_main.c, line 1349), TSAdaptChoose may shorten the last one
   or two steps to land on t_final (tsadapt.c, lines 889–905). This happens even with NONE.
5. **Monitors that also change dt after the step.**
   - InterfaceCFLMonitor applies snap-to-dtmax and the |Δφ| cap (src/monitoring.c, lines 463–503).
   - BoundsRollbackPreStep halves dt (src/monitoring.c, lines 530–534, with the halving at line 153).
   - These act after the patch's choice, so the effective dt rule of the anchor solver is the composition: patch, then
     CFL cap, then bounds rollback.
6. **What a full-Jacobian test must reproduce for TSALPHA.**
   - The SNES residual is (1/α_f)·F(t + α_fΔt, X_a, V_a) (alpha1.c, lines 64, 66 and 370–373).
   - The SNES Jacobian is IJacobian(X_a, V_a, shift_V) with shift_V = α_m/(α_f γ Δt), without the 1/α_f factor
     (alpha1.c, lines 65 and 377–386).
   - These agree, because dX_a/dX_1 = α_f and dV_a/dX_1 = α_m/(γΔt).
   - So a direct test should evaluate F and the IJacobian at a given (U, U̇) with a chosen shift. It should compare J·v with
     [F(U + εv, U̇ + shift·εv) − F(U − εv, U̇ − shift·εv)]/(2ε).
   - The existing `-test_wall_jacobian` uses shift = 0 and V = 0 (src/lunar_main.c, lines 1853–1854 and 1878), so it
     never exercises a ∂F/∂U̇ entry. Examples are the `shift * rho * cp` entry and the latent `−rho_ice * lat_sub * shift`
     entry (src/assembly.c, lines 344 and 347).
7. **Relation to the toy.** BE is generalized-α with α_m = α_f = γ = 1 (alpha1.c, lines 625–626). TSAlpha_Restart uses
   exactly this setting for its two half-steps (line 126).

## Solver settings (proposed)

| Setting | Proposal | Reason |
|---|---|---|
| Time integrator | BE, fixed step, Δτ from the caller (`dt`); helper M4-8 to choose it | Decided at G:Selection. Unconditionally stable and monotone (D3, D4). No adaptivity is needed: one transient of known decay, with no events. |
| Nonlinear solve | One Newton step per time step: U^{n+1} = U^n − J⁻¹R(U^n), then evaluate R(U^{n+1}) once as a guard | R is affine in U^{n+1} (M4-6b), so one Newton step from any guess is exact up to roundoff. This is algebraically the same as solving (I − Δτ D_h)U^{n+1} = U^n + Δτ(Q + f) + Δτ·b, with b the boundary vector. The Newton form is chosen so that the code calls the same `residual` and `jacobian` that the Jacobian test checks, as the C solver will. |
| Guard after the solve | Raise if ‖R(U^{n+1})‖_∞ > c·u·(shift + 4/h²)·max(1, ‖U^{n+1}‖_∞), where u is double-precision unit roundoff and c is a constant | Roundoff in the matrix-vector product scales with the row's largest coefficients. The constant c (I would start at 10³) is a tolerance, so it is for M5 and the verifier to set, not me. |
| Linear solver | Thomas algorithm (tridiagonal LU without pivoting), O(N) | J is a nonsingular M-matrix, and diagonally dominant in both the transient and the steady cases, so no pivoting is needed (M4-7). The brief allows NumPy only (problem.md, line 16), so no scipy banded solver. |
| Preconditioner, KSP | None | Direct solve of a 63×63 tridiagonal system. |
| Precision | float64 | |
| Steady solve | The same Thomas solve with shift = 0 and f = 0 (M4-5) | V1 (c) |

**Performance.** About 8N flops per step, so negligible. I make no speed claim, because none has been measured. M8 is
skipped for the toy (problem.md, line 26).

## Manufactured solutions

The script is `derivation/mms/mms_sources.py`. It prints f, checks the data, and prints sample values.

**M4-9** (operation: choose θ′_m, then substitute it into M4-1 to define f):

$$\theta'_m(\xi,\tau) = 1 + A\sin(\pi\xi)\,(1-e^{-\tau}) + B\,\xi(1-\xi)\sin(2\pi\tau)$$
$$f(\xi,\tau) = \partial_\tau\theta'_m - \partial_\xi^2\theta'_m - Q = A\sin(\pi\xi)\big[e^{-\tau} + \pi^2(1-e^{-\tau})\big] + 2\pi B\,\xi(1-\xi)\cos(2\pi\tau) + 2B\sin(2\pi\tau) - Q \tag{M4-9}$$

Defaults: A = 0.5, B = 0.25, and Q from the parameters (1 in C1).

**Why this θ′_m.**
- **It meets the production data exactly.** θ′_m(0, τ) = θ′_m(1, τ) = 1, and θ′_m(ξ, 0) = 1. So the MMS option changes
  only the forcing (T5). The initial condition and the end values stay as in production, and the option tests the
  production code path.
- **Its corners are compatible.** At (0, 0) and (1, 0), ∂_ξ²θ′_m + Q + f = 0, which equals the ∂_τθ′ = 0 that the end
  condition forces. So the solution is C^∞ on the closed domain and full order is expected. This differs from the
  production problem, whose corners are incompatible (selection.md item 1).
- **The quadratic part is not enough alone.** D_h is exact on the quadratic B-term (D1), so that term would not show the
  spatial order. The sin(πξ) term does.
- **Both time behaviours are present.** The B-term oscillates in τ, so the time error does not decay away. The A-term
  approaches a nonzero steady state.

**Scratch check** (`m4_mms_copy.py`, a verbatim copy of the body of `mms_sources.py`). The stage guard runs only scripts
under `derivation/scratch/`, so I did not run `mms/mms_sources.py` itself. The copy's output:
- The PDE residual is 0 symbolically.
- The end values minus 1 are 0, and the initial value minus 1 is 0. Both corner checks are 0.
- A = B = 0 gives f = −Q: the forcing cancels the source and θ′ stays 1.
- Samples with A = 0.5, B = 0.25, Q = 1:

| (ξ, τ) | θ′_m | f |
|---|---|---|
| (0.5, 0.1) | 1.084317869250300 | 0.5336200848195806 |
| (0.25, 0.5) | 1.139112419353097 | 0.2929012062596490 |

**Dimensional form**, for a solver that works in SI: q_m(x, t) = (kT_0/L²)·f(x/L, κt/L²) W m⁻³. The script prints it.

**How the option enters.**
- `params["mms"] = {"A": 0.5, "B": 0.25}` turns it on; `None`, the default, turns it off.
- When on, the term −f(ξ_j, τ^{n+1}) is subtracted in every row (T5), evaluated at the new time level as BE requires
  (M4-4).
- When off, the forcing code path is skipped entirely, not evaluated as zero. So production runs are bitwise unchanged.

**Expected orders, from my scratch run of the same scheme** (`m4_checks.py` item 7; error = max over nodes at τ_end = 0.5):

| Refinement | Errors | Observed order |
|---|---|---|
| Space, Δτ = 10⁻⁵, N = 16 → 32 | 5.2340e-4 → 1.3016e-4 | 2.008 |
| Time, N = 256, Δτ = 10⁻² → 5×10⁻³ | 7.5869e-4 → 3.7746e-4 | 1.007 |

At the reference N = 64, Δτ = 10⁻³ the error is 4.3056e-5. The test tolerances belong to M5.

## Layers

Each layer adds code behind an option that is off by default, so that production behaviour is unchanged. Each layer must
pass the listed tests, and then all earlier tests, before the next layer starts. The test IDs below are proposals; the
verifier writes the tests.

**Layer 0: infra.** Code: the T5 MMS forcing behind `params["mms"]`, and `check_jacobian`.

`check_jacobian(n, params, U=None, Udot=None, shifts=(0, 1e1, 1e3, 1e5), eps=1e-6, seed=0)` does the following:
- It evaluates F and J at a given state, or by default at a seeded random state around 1.
- It compares J·v with the central difference of F along (v, shift·v), for 5 random directions.
- It reports the relative error per block and per part. The toy has one field, so there is one block, θ′/θ′. Within it,
  the report is split three ways:
  - the boundary-adjacent rows {1, N−1} against the interior rows;
  - the full J;
  - ∂F/∂U and ∂F/∂U̇ separately.
- It has no effect on `solve`.

Tests:
- **infra-1.** With the option off, `residual` returns values bitwise identical to the residual with the T5 code removed.
  In practice: equal to T1 + T2 + T3 + T4 computed independently.
- **infra-2.** With the option on, f(ξ_j, τ) matches the samples printed by `mms_sources.py` to roundoff.
- **infra-3.** `check_jacobian` passes at shifts 0, 10, 10³ and 10⁵, in every reported group.
- **infra-4.** Sensitivity of `check_jacobian`: on a Jacobian with the shift term removed, or with an off-diagonal sign
  flipped, it reports failure.

**Layer 1: heat (the only physics layer).** Code: T1–T4, J (M4-7), the Thomas solve, BE stepping, `solve_steady`,
`end_fluxes` and `dt_for_tolerance`. Tests: the layer's tests heat-1 to heat-8, then infra-1 to infra-4 again.
- **heat-1. MMS spatial order 2** at a fixed small Δτ, and **heat-2. MMS temporal order 1** at a fixed fine N (problem.md,
  line 25).
- **heat-3.** The direct steady solve equals the parabola S44 at the nodes to roundoff (D2, M4-5).
- **heat-4.** At τ_end = 1 with the reference discretization, T is compared with S44. The tolerance must be looser than
  1e-5 in θ (D-d). The fully discrete floor is 7.0e-6 in θ.
- **heat-5.** Bounds 0 ≤ θ ≤ 1/8 at every node and step for q > 0 (D4). In θ′ this is 1 ≤ θ′ ≤ 1 + Q/8.
- **heat-6.** At the steady state, the half-cell end fluxes equal qL/2 (S45, V4).
- **heat-7.** A non-unit parameter set, for example L = 2, k = 3, ρ = 5, c_p = 7, q = 11, T_0 = 13, still reproduces S44
  with the direct steady solve. **This is the only test that can catch a scaling error.** In C1 every scale factor is 1,
  so a wrong κ, Q or kT_0/L² is invisible there.
- **heat-8.** q = 0 gives T ≡ T_0 exactly (S37), and q < 0 runs with T ≤ T_0 (C2, S38).

**Layers vapor, phase field, coupling and wall: not applicable.** There is no vapor, ice or grain (S8), and no wall
condition (S40). They are recorded as absent, with no code and no tests.

## What I could not show, or did not do

- **The dt rule is a first-order mode-1 estimate,** checked numerically within 3 % for Δτ ≤ 4e-3, not proved as a bound
  (§The dt rule).
- **The patched TSALPHA fields** `adap`, `NRmin`, `NRmax`, `factor`, `dtmin`, `dtmax` and `maxdt_dyn` are declared in a PETSc
  header I could not read, so their declarations and defaults are unchecked.
- **`mms/mms_sources.py` was not run by me as itself.** Its verbatim copy in scratch was run.
- **The new discretization symbols** (table at the top) are declared here, not in ledger.json. This stage's outputs do not
  include the ledger (question 1).
