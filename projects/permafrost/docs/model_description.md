# Permafrost Phase-Field Model: Mathematical and Numerical Description

> **Purpose:** Reference document for manuscript preparation.  
> **Code:** `PetIGA-3.20/projects/permafrost/`  
> **Last updated:** 2026-05-15

---

## 1. Overview

The model describes the coupled evolution of ice, sediment (mineral grains), and air within a permafrost pore space. The state is represented by four prognostic degrees of freedom at every point in the domain:

| DOF | Symbol | Description |
|-----|--------|-------------|
| 0 | φ_i | Ice volume fraction (phase field) |
| 1 | T | Temperature (°C) |
| 2 | ρ_v | Vapor mass density (kg m⁻³) |
| 3 | φ_s | Sediment volume fraction (phase field) |

The air volume fraction is a derived quantity:

$$\phi_a = 1 - \phi_i - \phi_s$$

The system is closed by requiring φ_i, φ_s, φ_a ∈ [0, 1] everywhere. All three phase-field variables satisfy the partition-of-unity constraint φ_i + φ_s + φ_a = 1 identically.

---

## 2. Three-Phase Free Energy

The bulk free energy follows the **Kim–Steinbach three-phase model** (Kim et al. 1999; Steinbach & Pezzolla 1999). The interfacial free energy density is:

$$\Psi(\phi_i, \phi_s, \phi_a) = \sum_{k \in \{i,s,a\}} \eta_k \phi_k^2(1-\phi_k)^2 + \Lambda\,\phi_i^2\phi_s^2\phi_a^2$$

The first sum provides individual double-well potentials for each phase; the triple-junction interaction term $\Lambda\phi_i^2\phi_s^2\phi_a^2$ penalises configurations where all three phases coexist simultaneously, confining triple lines to a narrow region near grain boundaries.

The **bulk driving force functions** entering the Allen-Cahn equations are defined as:

$$f_i(\phi_i,\phi_s) = \eta_i\phi_i(1-\phi_i)(1-2\phi_i) + 2\Lambda\,\phi_i\phi_s^2\phi_a^2$$
$$f_s(\phi_i,\phi_s) = \eta_s\phi_s(1-\phi_s)(1-2\phi_s) + 2\Lambda\,\phi_i^2\phi_s\phi_a^2$$
$$f_a(\phi_i,\phi_s) = \eta_a\phi_a(1-\phi_a)(1-2\phi_a) + 2\Lambda\,\phi_i^2\phi_s^2\phi_a$$

These are related to the partial derivatives of Ψ with respect to each phase fraction. The **phase-specific energy coefficients** η_k are derived from pairwise interface energies (J m⁻²):

$$\eta_i = \gamma_{iv} + \gamma_{is} - \gamma_{sv}$$
$$\eta_s = \gamma_{sv} + \gamma_{is} - \gamma_{iv}$$
$$\eta_a = \gamma_{iv} + \gamma_{sv} - \gamma_{is}$$

with default values γ_iv = 0.109, γ_is = 0.033, γ_sv = 0.056 J m⁻². The parameter Λ = 1.0 controls triple-junction energy.

---

## 3. Governing Equations

### 3.1 Ice Phase — Allen-Cahn Equation

The ice phase evolves by a **non-conserved Allen-Cahn** equation:

$$\frac{\partial\phi_i}{\partial t} = M\left[\varepsilon^2\nabla^2\phi_i - \frac{1}{\varepsilon}\frac{\partial\Psi^{\mathrm{eff}}}{\partial\phi_i}\right] + S_{\mathrm{sub}}$$

where M is the interface mobility (m³ s⁻¹), ε is the diffuse interface half-width (m), and S_sub is a sublimation/condensation source term (see §4).

In the **three-phase (3P) regime** (before t_sed_freeze), the full Kim-Steinbach driving force is used:

$$\frac{\partial\Psi^{\mathrm{eff}}}{\partial\phi_i} \propto \frac{(\eta_s+\eta_a)f_i - \eta_a f_s - \eta_s f_a}{E_T}$$

where $E_T = \eta_i\eta_s + \eta_i\eta_a + \eta_s\eta_a$.

In the **two-phase (2P) regime** (after t_sed_freeze, with sediment frozen), the ice evolves against the stationary sediment boundary:

$$\frac{\partial\phi_i}{\partial t} = M_{2P}\left[\varepsilon^2\nabla^2\phi_i + \eta_a\varepsilon^2\nabla^2\phi_s + \frac{f_i - f_a}{\varepsilon}\right] + S_{\mathrm{sub}}$$

with $M_{2P} = 3M/(\eta_i + \eta_a)$, and φ_s is held fixed (∂φ_s/∂t = 0 in the 2P regime).

### 3.2 Sediment Phase — Allen-Cahn Equation

In the **three-phase regime**, sediment evolves symmetrically:

$$\frac{\partial\phi_s}{\partial t} = M\left[\varepsilon^2\nabla^2\phi_s + \frac{-\eta_a f_i - \eta_i f_a + (\eta_i+\eta_a)f_s}{\varepsilon E_T}\right]$$

In the **two-phase regime** (t ≥ t_sed_freeze), sediment is frozen:

$$\frac{\partial\phi_s}{\partial t} = 0$$

This is enforced by zeroing the sediment residual to its time-derivative term only, which forces ∂φ_s/∂t = 0 without modifying the ice or vapor equations. The sediment field remains as a static background geometry for the subsequent ice evolution.

### 3.3 Temperature — Thermal Energy Balance

$$\rho c_p \frac{\partial T}{\partial t} = \nabla\cdot(\lambda\nabla T) - \rho L_{\mathrm{sub}}\frac{\partial\phi_a}{\partial t}$$

where:
- ρ = weighted mixture density (kg m⁻³)
- c_p = weighted mixture heat capacity (J kg⁻¹ K⁻¹)
- λ = weighted mixture thermal conductivity (W m⁻¹ K⁻¹)
- L_sub = 2.83 × 10⁶ J kg⁻¹ (latent heat of sublimation)

Since φ_a = 1 − φ_i − φ_s, the latent heat source is:

$$-\rho L_{\mathrm{sub}}\frac{\partial\phi_a}{\partial t} = \rho L_{\mathrm{sub}}\left(\frac{\partial\phi_i}{\partial t} + \frac{\partial\phi_s}{\partial t}\right)$$

so latent heat is released when air fraction decreases (ice or sediment grows) and absorbed when ice melts or sublimes. A scaling factor ξ_T = 10⁻² is applied to the spatial (conduction + latent heat) terms to balance the time scales of thermal and phase-field dynamics.

### 3.4 Vapor Density — Penalised Diffusion

The vapor density equation enforces two physical constraints:
1. **Diffusive transport** of vapor through the air phase.
2. **Local vapor equilibrium** at ice and sediment surfaces.

$$\frac{\partial\rho_v}{\partial t} + \nabla\cdot\left(-D_{\mathrm{eff}}\,\phi_a^{\mathrm{eff}}\,\nabla\rho_v\right) + k_{\mathrm{pen}}\,H(\phi_i+\phi_s)\,(\rho_v - \rho_v^{\mathrm{eq}}) = \xi_v\,\rho_{\mathrm{ice}}\frac{\partial\phi_a}{\partial t}$$

**Penalised diffusivity:** The effective vapor diffusivity transitions smoothly between the full air value and a reduced (penalised) value inside solid phases:

$$D_{\mathrm{eff}} = D_v H(\phi_i) + D_{\mathrm{pen}}\left[1 - H(\phi_i)\right], \quad D_{\mathrm{pen}} = \alpha_{\mathrm{pen}} D_v$$

where H(φ) = φ³(3 − 2φ) is a smooth cubic Heaviside function and α_pen = 10⁻⁴ is the penalty factor. This reduces diffusivity by four orders of magnitude inside the ice phase while maintaining a smooth, differentiable transition.

**Interface equilibrium penalty:** The term k_pen H(φ_i+φ_s)(ρ_v − ρ_v^eq) drives ρ_v toward the local thermodynamic equilibrium value inside solid phases, with k_pen = 10⁵ Pa (penalty stiffness):

$$\rho_v^{\mathrm{eq}} = (\phi_i + \phi_s)\,\rho_{vs}(T) + \phi_a\,\rho_v$$

Here ρ_{vs}(T) is the saturation vapor density over ice (see §5), and H(φ_i+φ_s) activates the penalty only in solid regions.

**Air fraction limit:** φ_a^eff = max(φ_a, φ_a^lim) with φ_a^lim = 10⁻⁶ prevents numerical singularity when φ_a → 0.

**Mass conservation source:** The right-hand side ξ_v ρ_ice ∂φ_a/∂t = −ξ_v ρ_ice(∂φ_i/∂t + ∂φ_s/∂t) accounts for the change in vapor content due to phase-field evolution. A scaling factor ξ_v = 10⁻⁴ is applied to the spatial terms.

---

## 4. Sublimation and Condensation Kinetics

The **Allen-Cahn equation for ice** includes a phase-change source term:

$$S_{\mathrm{sub}} = \frac{\alpha_{\mathrm{sub}}\,\phi_i^2\phi_a^2}{\rho_{\mathrm{ice}}}\,(\rho_v - \rho_{vs})$$

This term localises the phase change to ice-air interfaces (through the $\phi_i^2\phi_a^2$ factor, which peaks at the interface where both φ_i and φ_a are non-zero) and drives:
- **Sublimation** when ρ_v < ρ_vs: ice decreases, φ_a increases.
- **Condensation** when ρ_v > ρ_vs: ice grows, φ_a decreases.

The kinetic coefficient α_sub and mobility M are derived from the **Gibbs–Thomson relation** using capillary length d₀ and kinetic coefficient β:

$$d_0 = \frac{d_0^0}{\rho_{\mathrm{ice}}/\rho_{vs}}, \quad \beta = \frac{\beta^0}{\rho_{\mathrm{ice}}/\rho_{vs}}$$

$$\lambda_{\mathrm{sub}} = a_1 \varepsilon / d_0, \quad \tau_{\mathrm{sub}} = \varepsilon\,\lambda_{\mathrm{sub}}\left(\frac{\beta}{a_1} + a_2\frac{\varepsilon}{D_{\mathrm{th}}} + a_2\frac{\varepsilon}{D_v}\right)$$

$$M = \frac{\varepsilon}{3\tau_{\mathrm{sub}}}, \quad \alpha_{\mathrm{sub}} = \frac{10\,\lambda_{\mathrm{sub}}}{\tau_{\mathrm{sub}}}$$

with a₁ = 5.0, a₂ = 0.1581, d₀⁰ = 10⁻⁹ m, β⁰ = 1.4 × 10⁵ s m⁻¹, and D_th = ½(λ_air/ρ_air c_{p,air} + λ_ice/ρ_ice c_{p,ice}).

---

## 5. Material Properties

All effective properties use volume-fraction-weighted averages over the three phases.

### 5.1 Thermal Properties

$$\lambda(\phi_i, \phi_s) = \phi_i\lambda_i + \phi_s\lambda_s + \phi_a\lambda_a$$
$$\rho c_p(\phi_i, \phi_s) = \phi_i\rho_i c_{p,i} + \phi_s\rho_s c_{p,s} + \phi_a\rho_a c_{p,a}$$

| Property | Ice | Sediment (metal) | Air |
|----------|-----|-----------------|-----|
| λ (W m⁻¹ K⁻¹) | 2.29 | 36.0 | 0.02 |
| ρ (kg m⁻³) | 919 | 7753 | 1.341 |
| c_p (J kg⁻¹ K⁻¹) | 1960 | 486 | 1044 |

### 5.2 Vapor Diffusivity

The vapor diffusivity in air follows a power law:

$$D_v(T) = D_0 \left(\frac{T + 273.15}{273.15}\right)^{1.81}$$

with D₀ = 2.178 × 10⁻⁵ m² s⁻¹.

### 5.3 Saturation Vapor Density Over Ice

The saturation vapor pressure over ice is computed from the empirical Antoine equation (Alduchov & Eskridge):

$$P_{vs}(T) = \exp\!\left(K_0 T_K^{-1} + K_1 + K_2 T_K + K_3 T_K^2 + K_4 T_K^3 + K_5\ln T_K\right)$$

with T_K = T + 273.15 (Kelvin) and coefficients K₀ = −5865, K₁ = 22.24, K₂ = 1.375 × 10⁻², K₃ = −3.403 × 10⁻⁵, K₄ = 2.697 × 10⁻⁸, K₅ = 0.6918. The saturation vapor density is:

$$\rho_{vs}(T) = \frac{0.62\,\rho_a\,P_{vs}}{P_{\mathrm{atm}} - P_{vs}}$$

where P_atm = 101325 Pa and ρ_a = 1.341 kg m⁻³.

### 5.4 Allen-Cahn Mobility

The phase-field mobility is interpolated linearly across phases:

$$M(\phi_i, \phi_s) = M_{\mathrm{sub}}\phi_i + M_s\phi_s + M_a\phi_a$$

where M_sub, M_s, M_a are computed from the kinetic parameters. By default M_s = M_a = M_sub so the mobility is spatially uniform; it can be overridden with `-mob_sed`.

---

## 6. Phase Model Switching: Three-Phase to Two-Phase

The simulation proceeds in two stages, separated by the parameter t_sed_freeze:

**Three-phase stage** (0 ≤ t < t_sed_freeze):  
Both φ_i and φ_s evolve under the full three-phase Allen-Cahn equations. This stage allows the initial microstructure to relax toward a thermodynamically consistent configuration with proper contact angles at ice-sediment-air triple lines.

**Two-phase stage** (t ≥ t_sed_freeze):  
The sediment field is frozen (∂φ_s/∂t = 0). Ice evolves under a reduced two-phase Allen-Cahn equation that treats the frozen sediment as a fixed boundary. This separation is physically motivated: on the time scales of vapor diffusion and ice sublimation, sediment grain rearrangement is negligible.

The transition is controlled by a single parameter:
- t_sed_freeze = 0: immediately two-phase from t = 0.
- t_sed_freeze = 300 s: three seconds of three-phase relaxation before freezing.

---

## 7. Relaxation Phase

At the start of each simulation, `n_relax` Allen-Cahn relaxation steps are performed before the full physics (vapor transport, sublimation, thermal) are activated. During relaxation:
- Ice and sediment evolve under three-phase Allen-Cahn only (no sublimation, no thermal coupling).
- T and ρ_v carry only time-derivative residual terms, forcing ∂T/∂t = ∂ρ_v/∂t = 0 and keeping the fields at their initial-condition values.

This smooths the initial condition from the specified profile (which has an exact tanh interface) to a locally equilibrated diffuse interface before phase-change kinetics are engaged.

---

## 8. Initial Conditions

The diffuse phase-field initial condition for a grain of phase k centred at **x**_c with radius R is:

$$\phi_k(\mathbf{x}) = \frac{1}{2} - \frac{1}{2}\tanh\!\left(\frac{\sqrt{2}}{\varepsilon}\left(\lvert\mathbf{x} - \mathbf{x}_c\rvert - R\right)\right)$$

The factor √2/ε sets the interface width to ε (half-width at half-height of the tanh profile). Initial humidity prescribes the initial vapor density:

$$\rho_v^0 = h_0 \rho_{vs}(T^0) \cdot \phi_a + \rho_{vs}(T^0) \cdot (1 - \phi_a)$$

where h₀ ∈ (0,1] is the relative humidity parameter. In air regions φ_a ≈ 1, so ρ_v⁰ ≈ h₀ ρ_vs; inside solid phases the local vapor density is set to the saturation value ρ_vs.

Available initial-condition types:

| `-ic_type` | Description |
|-----------|-------------|
| `ice_slab` | Centered ice slab with sediment half-domain |
| `enclosed` | Pair of enclosed grains (ice shell + sediment core) |
| `single_ice` | Single pure ice grain centred in domain |
| `ice_sed_pair` | One ice grain + one sediment grain, evenly spaced |
| `slab_and_grains` | Ice slab on one side, random grains on the other |
| `capillary` | Capillary bridge geometry |
| `ice_cap` | Ice cap on sediment surface |

---

## 9. Numerical Methods

### 9.1 Isogeometric Analysis (IGA)

The spatial discretisation uses **Isogeometric Analysis (IGA)** via the PetIGA library (Dalcin et al. 2016). IGA uses B-spline basis functions that provide:
- **C¹ continuity** between elements (with polynomial order p = 2 and continuity C = 1).
- Smooth representation of field variables through the diffuse interface, important for accurate evaluation of interface gradients.
- Exact geometry representation.

The domain [0, Lx] × [0, Ly] (or 1D/3D variants) is discretised with N_x × N_y uniform elements. The basis functions N_A(ξ) are B-splines of degree p = 2 with C¹ global continuity, yielding (N + p) basis functions per direction.

The element connectivity and shape function evaluation are handled by PetIGA. The physical-space element stiffness matrices are assembled from reference-space integrals via the Gauss quadrature rule with (p+1)² = 9 quadrature points per 2D element.

### 9.2 Weak (Variational) Formulation

Multiplying each residual equation by a test function N_A and integrating by parts over the domain Ω, the element residual vector is:

**Ice (3-phase):**
$$R_A^i = \int_\Omega N_A \dot\phi_i\,d\Omega + 3M\varepsilon\int_\Omega \nabla N_A\cdot\nabla\phi_i\,d\Omega + \frac{C_3}{\varepsilon}\int_\Omega N_A\,\mu_i^{\mathrm{eff}}\,d\Omega - \int_\Omega N_A S_{\mathrm{sub}}\,d\Omega$$

where $\mu_i^{\mathrm{eff}} = [(\eta_s+\eta_a)f_i - \eta_a f_s - \eta_s f_a]\varepsilon$ and $C_3 = 3M/(E_T\varepsilon)$.

**Temperature:**
$$R_A^T = \int_\Omega \rho c_p N_A \dot T\,d\Omega + \xi_T\int_\Omega \lambda\nabla N_A\cdot\nabla T\,d\Omega + \xi_T\int_\Omega \rho L_{\mathrm{sub}} N_A \dot\phi_a\,d\Omega$$

**Vapor density:**
$$R_A^v = \int_\Omega N_A \dot\rho_v\,d\Omega + \xi_v\int_\Omega D_{\mathrm{eff}}\phi_a^{\mathrm{eff}}\nabla N_A\cdot\nabla\rho_v\,d\Omega + \xi_v\int_\Omega N_A k_{\mathrm{pen}} H(\phi_i+\phi_s)(\rho_v - \rho_v^{\mathrm{eq}})\,d\Omega - \xi_v\int_\Omega N_A\rho_{\mathrm{ice}}\dot\phi_a\,d\Omega$$

Natural boundary conditions (zero flux) are satisfied automatically by this formulation. No Dirichlet conditions are imposed on φ_i, φ_s by default; temperature and vapor density can optionally be fixed at boundaries via `-flag_BC_Tfix` and `-flag_BC_rhovfix`.

### 9.3 Analytical Jacobian

The element Jacobian J[a][i][b][j] = ∂R_a^i/∂U_b^j + shift × ∂R_a^i/∂(dU_b^j/dt) is computed analytically for Avenue 1. The full 4 × 4 block structure is assembled, including all cross-couplings:

- **[ice, ice]:** diagonal + AC stiffness + sublimation derivative d(loc)/d(φ_i)
- **[ice, T]:** sublimation via d(ρ_vs)/dT
- **[ice, ρ_v]:** sublimation linear term −α_sub loc/ρ_ice
- **[ice, sed]:** AC cross-coupling + sublimation via d(loc)/d(φ_s)
- **[T, ice]:** latent heat via shift×(−ξ_T ρ L_sub) + thermal conductivity variation
- **[T, T]:** diagonal + thermal diffusion stiffness
- **[T, sed]:** latent heat via shift×(−ξ_T ρ L_sub)
- **[ρ_v, ice]:** penalty gradient + diffusivity variation d(D_eff)/d(φ_i)
- **[ρ_v, T]:** penalty via d(ρ_vs)/dT + diffusivity via d(D_v)/dT
- **[ρ_v, ρ_v]:** diagonal + diffusion stiffness + penalty
- **[ρ_v, sed]:** penalty + (symmetry with [ρ_v, ice])
- **[sed, *]:** sediment rows match ice rows by symmetry in 3P; only mass matrix in 2P

The Jacobian is assembled once per Newton iteration and passed to the GMRES linear solver. The analytical form avoids finite-difference perturbation, which degrades in accuracy near the stiff interface penalty.

### 9.4 Time Integration: Generalized-α Method

Time integration uses the **generalized-α (TSALPHA)** scheme, which is second-order accurate and unconditionally stable for linear problems. The method introduces parameters α_m, α_f, and γ that control dissipation and accuracy. PETSc's TSALPHA with spectral radius ρ_∞ = 0.5 is used as the default.

At each time step, the system of equations to be solved is:

$$\mathbf{R}(\mathbf{U}_{n+\alpha_f},\, \dot{\mathbf{U}}_{n+\alpha_m}) = \mathbf{0}$$

where $\mathbf{U}_{n+\alpha_f} = (1-\alpha_f)\mathbf{U}_n + \alpha_f\mathbf{U}_{n+1}$ is a weighted average of the solution at the current and next time level, and similarly for the time derivative. This approach provides controllable high-frequency dissipation while retaining second-order accuracy on smooth solutions.

### 9.5 Nonlinear Solver: Newton-Krylov

At each time step the nonlinear residual R(U) = 0 is solved with **SNES newtonls** (Newton with backtracking line search). The Newton update satisfies:

$$J_n\,\delta U_n = -R(U_n), \quad U_{n+1} = U_n + \alpha_n\,\delta U_n$$

where α_n ∈ (0,1] is a backtracking step size and J_n = ∂R/∂U + shift × ∂R/∂(dU/dt) is the Jacobian including the time-derivative contribution from the generalized-α scheme.

**Solver parameters:**
- snes_max_it = 15 (maximum Newton iterations)
- snes_rtol = 10⁻⁶ (relative residual tolerance)
- snes_atol = 10⁻¹⁰ (absolute residual tolerance)
- snes_stol = 10⁻⁶ (step size tolerance)

### 9.6 Per-DOF Convergence Test

A custom convergence test (`SNESDOFConvergence`) checks each DOF independently rather than using the global residual norm. At iteration k, convergence is declared when all four DOFs satisfy at least one criterion:

$$\|\mathbf{r}_j^{(k)}\|_2 \leq r_{\mathrm{tol}}\,\|\mathbf{r}_j^{(0)}\|_2 \quad\text{(relative)}$$
$$\|\mathbf{r}_j^{(k)}\|_2 \leq a_{\mathrm{tol}} \quad\text{(absolute)}$$
$$\|\delta\mathbf{u}_j^{(k)}\|_2 \leq s_{\mathrm{tol}}\,\|\mathbf{u}_j^{(k)}\|_2 \quad\text{(step size)}$$

for j ∈ {0=ice, 1=T, 2=ρ_v, 3=sed}. A minimum of 3 Newton iterations (1 during relaxation) is required before convergence can be declared, preventing premature exit on the first Newton step when the initial residual is artificially small.

### 9.7 Linear Solver: GMRES with Additive Schwarz Preconditioning

Each Newton iteration requires solving the 4 × 4 block linear system J δU = −R. The linear solver stack is:

- **Krylov method:** GMRES (Saad & Schultz 1986)
  - ksp_rtol = 10⁻⁶, ksp_atol = 10⁻¹⁰
  - ksp_max_it = 500
  - ksp_gmres_restart = 200 (restart Krylov basis after 200 iterations)
  
  GMRES is the correct choice for this system: the Jacobian is asymmetric due to the off-diagonal coupling between DOFs (sublimation links ice to vapor; latent heat links ice to temperature). CG or MINRES would fail on a non-symmetric system.

- **Preconditioner:** Additive Schwarz Method (ASM) with ILU(2) sub-block factorisation
  - pc_type = asm, pc_asm_overlap = 1
  - sub_pc_type = ilu, sub_pc_factor_levels = 2
  
  ASM partitions the domain into MPI blocks with one layer of overlap. Within each block, ILU(2) (incomplete LU with fill level 2) is applied. This captures more of the 4-DOF coupled B-spline stencil than ILU(1) (which underfills the 9-point 2D stencil), and the overlap prevents poor conditioning at block boundaries.

- **1D override:** For 1D problems (N ≤ 760 DOFs), direct sparse LU is used:
  - ksp_type = preonly, pc_type = lu
  
  This is trivially cheap in 1D and eliminates all Krylov iteration overhead.

### 9.8 Adaptive Time Stepping

The time step is adapted based on the number of Newton iterations required in the previous solve. Three counters from the `TSADAPT` structure drive the adaptation:

| Condition | Action |
|-----------|--------|
| Newton iters ≤ NRmin = 3 | Increase dt by factor f = 10^(1/8) ≈ 1.33 |
| NRmin < Newton iters ≤ NRmax = 8 | Keep dt unchanged |
| Newton iters > NRmax = 8 | Decrease dt by factor 1/f; reject and retry |
| Rejected steps > max_rej = 10 | Abort simulation |

The default step adjustment factor is f = 10^(1/8), providing slow growth and rapid reduction. The time step is bounded by dtmin = 0.01 × Δt₀ and dtmax = 0.5 × t_out (half the output interval, to ensure snapshots land exactly at output times).

---

## 10. Computational Framework

- **PETSc** (Balay et al. 2023): nonlinear and linear solvers, time integrators, parallel data structures.
- **PetIGA** (Dalcin et al. 2016): B-spline IGA framework built on PETSc; handles element assembly, quadrature, and IGA-specific boundary conditions.
- **MPI parallelism:** domain decomposition by PetIGA; one MPI rank per process.
- **Build system:** GNU Make; configured via PETSc PETSC_DIR/PETSC_ARCH variables.

---

## 11. Diagnostics and Output

At each time step a monitor table (`SNESDOFConvergence`) prints the 2-norm of the residual and Newton step for each DOF, along with the relative and step-size convergence ratios. This allows tracking of which physical field is the convergence bottleneck.

Scalar time-series output is written to `SSA_evo.dat`, including:

| Column | Quantity |
|--------|----------|
| Time (s) | Simulation time |
| ⟨φ_i⟩ | Volume-averaged ice fraction |
| φ_i² φ_s² φ_a² integral | Triple-junction area proxy |
| ⟨φ_a⟩ | Volume-averaged air fraction |
| ⟨T⟩ | Volume-averaged temperature |
| ⟨ρ_v φ_a⟩ | Total vapor mass per unit area |
| φ_i² φ_a² integral | Ice-air interface area |
| ⟨φ_s⟩ | Volume-averaged sediment fraction |
| φ_s² φ_a² integral | Sediment-air interface area |
| φ_s² φ_i² integral | Ice-sediment contact area |

Full solution fields (φ_i, T, ρ_v, φ_s, φ_a) are written as binary PetIGA `sol_*.dat` files and can be converted to VTK format for ParaView visualisation.

---

## 12. Parameter Summary

| Parameter | Symbol | Value | Source |
|-----------|--------|-------|--------|
| B-spline degree | p | 2 | `universal.opts` |
| B-spline continuity | C | 1 (C¹) | `universal.opts` |
| Interface half-width (standard) | ε | 7.12 × 10⁻⁷ m | test opts |
| Interface half-width (hi-res) | ε | 3.56 × 10⁻⁷ m | test opts |
| Ice-vapor surface energy | γ_iv | 0.109 J m⁻² | `permafrost2.c` |
| Ice-sediment surface energy | γ_is | 0.033 J m⁻² | `permafrost2.c` |
| Sediment-vapor surface energy | γ_sv | 0.056 J m⁻² | `permafrost2.c` |
| Triple-junction penalty | Λ | 1.0 | `permafrost2.c` |
| Latent heat of sublimation | L_sub | 2.83 × 10⁶ J kg⁻¹ | `permafrost2.c` |
| Air diffusivity floor | φ_a^lim | 10⁻⁶ | `permafrost2.c` |
| Diffusivity penalty factor | α_pen | 10⁻⁴ | `universal.opts` |
| Interface equilibrium stiffness | k_pen | 10⁵ Pa | `universal.opts` |
| Vapor time scaling | ξ_v | 10⁻⁴ | `permafrost2.c` |
| Thermal time scaling | ξ_T | 10⁻² | `permafrost2.c` |
| Capillary length scale | d₀⁰ | 10⁻⁹ m | `permafrost2.c` |
| Kinetic coefficient | β⁰ | 1.4 × 10⁵ s m⁻¹ | `permafrost2.c` |
| Relaxation steps | n_relax | 1 | `universal.opts` |
| Sed. freeze time | t_sed_freeze | 1 s | `universal.opts` |

---

## 13. Why the Current Model Eliminates Spurious Air

Comparison with results from 2026-05-05 (run `test_1D_IceSlab_2Phase_difvappen1e-07_k_pen1e09`) reveals four root causes of spurious air that have since been corrected:

### 13.1 Extreme Interface Equilibrium Stiffness (k_pen)

**Old:** k_pen = 10⁹ Pa.  
**Current:** k_pen = 10⁵ Pa.

The penalty term k_pen H(φ_i+φ_s)(ρ_v − ρ_eq) enforces vapor equilibrium in solid phases. With k_pen = 10⁹, the vapor equation in the solid is dominated by a stiff algebraic constraint (residual ∼ 10⁹ × δρ_v). The Jacobian condition number scales with k_pen, making the linear system nearly singular and the Newton solver unreliable. The inaccurate linear solves caused Newton to converge to wrong roots, producing ρ_v values that overshoot or undershoot the equilibrium, which in turn drove the sublimation source term S_sub erratically. The Allen-Cahn equation for ice then received spurious driving forces, pushing φ_i below zero (creating spurious air) in regions where it should have remained solid.

With k_pen = 10⁵, the penalty is still strong enough to maintain near-equilibrium vapor in the solid (the mismatch |ρ_v − ρ_eq| < 10⁻⁵/k_pen kg m⁻³) but the linear system is well-conditioned and Newton converges cleanly.

### 13.2 Near-Discontinuous Diffusivity Transition (difvap_pen)

**Old:** α_pen = 10⁻⁷ (D_eff = D_v × 10⁻⁷ inside ice — seven orders of magnitude reduction).  
**Current:** α_pen = 10⁻⁴ (four orders of magnitude reduction).

The effective diffusivity D_eff transitions from D_v (in air) to α_pen × D_v (in ice) over the finite interface width. With α_pen = 10⁻⁷, the jump is so extreme that the smooth Heaviside H(φ_i) cannot adequately resolve it on the numerical mesh — individual Gauss points near the interface see wildly different diffusivity values depending on which side of the midpoint they fall. This creates a nearly-discontinuous coefficient in the vapor diffusion term, which the finite-difference Jacobian approximation (used in the old code) cannot handle accurately. The resulting Jacobian errors caused the linear solve to be inaccurate at the interface, allowing ρ_v gradients to develop that are inconsistent with the physical boundary condition, again driving spurious phase-field evolution.

### 13.3 Finite-Difference Jacobian → Analytical Jacobian

**Old:** `Jacobian()` returned 0, so PETSc computed the Jacobian by finite differences.  
**Current:** Full analytical Jacobian_A1, including all 4 × 4 cross-coupling blocks.

A finite-difference Jacobian perturbation δ ∼ √ε_machine ≈ 10⁻⁸ introduces O(ε_machine/δ) ≈ O(10⁻⁸) relative errors in each Jacobian entry. For entries scaled by k_pen = 10⁹, the absolute error in the [ρ_v, ρ_v] Jacobian block is O(10⁹ × 10⁻⁸) = O(10), which is comparable to the residual itself. Newton's convergence degrades from quadratic to linear or worse when the Jacobian contains O(1) errors in dominant entries.

The analytical Jacobian carries no such perturbation error. Newton converges quadratically once the linear solve is accurate, and the computed solution satisfies all four field equations simultaneously.

### 13.4 Undersaturated Initial Condition

**Old:** humidity = 0.5 (50% of saturation).  
**Current:** humidity = 0.95 (95% of saturation).

With h₀ = 0.5, the initial vapor density everywhere in the air phase is ρ_v = 0.5 ρ_vs, creating an immediate strong sublimation driving force (ρ_v − ρ_vs = −0.5 ρ_vs < 0) at every ice-air interface. With a poorly-conditioned solver, this large source term magnifies any Newton convergence error into a visible change in φ_i, and if φ_i overshoots downward (φ_i < 0) in any element, the subsequent time step inherits a corrupted initial condition. With h₀ = 0.95, the system is near equilibrium and the sublimation driving force is 10× smaller, giving the solver much more margin for error.

---

## 14. References

- Balay, S. et al. (2023). *PETSc/TAO User's Manual*. Argonne National Laboratory.
- Dalcin, L. et al. (2016). PetIGA: A framework for high-performance isogeometric analysis. *Comput. Methods Appl. Mech. Eng.* 308, 151–181.
- Kim, S.G., Kim, W.T., Suzuki, T. (1999). Phase-field model for binary alloys. *Phys. Rev. E* 60, 7186.
- Saad, Y., Schultz, M.H. (1986). GMRES: A generalized minimal residual algorithm for solving nonsymmetric linear systems. *SIAM J. Sci. Stat. Comput.* 7, 856–869.
- Steinbach, I., Pezzolla, F. (1999). A generalized field method for multiphase transformations using interface fields. *Physica D* 134, 385–393.
