# Permafrost Phase-Field Model: Mathematical and Numerical Description

> **Purpose:** Reference document for manuscript preparation.  
> **Code:** `PetIGA-3.20/projects/permafrost/`  
> **Last updated:** 2026-05-27

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

$$\rho c_p \frac{\partial T}{\partial t} \;=\; \xi_T\,\nabla\cdot(\lambda\nabla T) \;+\; \xi_T\,\rho L_{\mathrm{sub}}\,S_{\mathrm{sub}}$$

where:
- ρ = weighted mixture density (kg m⁻³)
- c_p = weighted mixture heat capacity (J kg⁻¹ K⁻¹)
- λ = weighted mixture thermal conductivity (W m⁻¹ K⁻¹)
- L_sub = 2.83 × 10⁶ J kg⁻¹ (latent heat of sublimation)
- $S_{\mathrm{sub}}$ = sublimation source defined in §4

**Pairing latent heat with $S_{\mathrm{sub}}$ rather than with $\partial\phi_a/\partial t$.** The latent-heat source is paired *directly* with the sublimation rate $S_{\mathrm{sub}}$ (the real ice ↔ vapor mass-exchange term in the ice equation), **not** with $\partial\phi_a/\partial t = -\partial\phi_i/\partial t - \partial\phi_s/\partial t$. Using $\partial\phi_a/\partial t$ conflates Allen-Cahn interface motion — mass-neutral rearrangement of the diffuse interface during initial relaxation — with real phase change. With $S_{\mathrm{sub}}$ instead, latent heat is exchanged **only** when ice and vapor actually exchange mass. This pairing is the thermal-side counterpart of the §3.4 mass-balance closure, and was a critical fix for stability during the initial AC relaxation phase. Sign convention: $S_{\mathrm{sub}} > 0$ means ice is *growing* (deposition; vapor → ice), which physically *releases* latent heat (T increases) — matching the sign of the source term above.

ξ_T = 1.0 by default (physical thermal time scaling).

### 3.4 Vapor Density — Diffusion + Stefan Source (Penalties Optional)

The vapor density equation has three terms in the residual: a diffusive transport term through air, an *optional* numerical penalty pulling $\rho_v$ toward an equilibrium anchor in deep solid, and a Stefan-condition mass-exchange source coupling ice and vapor:

$$\frac{\partial\rho_v}{\partial t} \;-\; \nabla\cdot\!\left(\xi_v\, D_{\mathrm{eff}}(\phi)\,\phi_a^{\mathrm{eff}}\,\nabla\rho_v\right) \;+\; \xi_v\, k_{\mathrm{pen}}\, G_{\mathrm{pen}}(\phi_i+\phi_s)\,(\rho_v - \rho_v^{\mathrm{eq}}) \;=\; V_{\mathrm{src}}$$

with $V_{\mathrm{src}} = -\rho_{\mathrm{ice}}\,S_{\mathrm{sub}}$.

**Stefan-condition source $V_{\mathrm{src}}$ (Moure & Fu 2024 form).** The right-hand side is paired *directly* with the sublimation kinetic source $S_{\mathrm{sub}}$ (§4):

$$V_{\mathrm{src}} \;=\; -\,\rho_{\mathrm{ice}}\,S_{\mathrm{sub}} \;=\; -\,\alpha_{\mathrm{sub}}\,\phi_i^2\,\phi_a^2\,\bigl(\rho_v - \rho_{vs}^{\mathrm{eff}}\bigr)$$

(the $\rho_{\mathrm{ice}}$ factor cancels the $1/\rho_{\mathrm{ice}}$ in $S_{\mathrm{sub}}$). This is driven by the **phase-change** term only, not by $\partial\phi_i/\partial t$. The reasons:
- **Allen-Cahn curvature motion is mass-neutral.** $\partial\phi_i/\partial t$ contains both real phase change ($S_{\mathrm{sub}}$) *and* AC curvature relaxation of the diffuse interface. The latter is a numerical artifact of the gradient flow, not physical mass exchange. Coupling it to vapor — as the earlier $V_{\mathrm{src}} = -2\rho_{\mathrm{ice}}\phi_a\,\partial\phi_i/\partial t$ form did — launders AC "mass loss" into spurious vapor (a measured ~0.6 % `TOTAL_MASS` drift during the initial AC relaxation burst).
- **Pairing with $S_{\mathrm{sub}}$ conserves the physical exchange exactly.** $R^i$ carries $-S_{\mathrm{sub}}$ and $R^v$ carries $-V_{\mathrm{src}} = +\rho_{\mathrm{ice}}S_{\mathrm{sub}}$, so every gram of ice that sublimates becomes exactly one gram of vapor, pointwise, at every quadrature point. The $\phi_i^2\phi_a^2$ localiser zeroes the source at ice-sed boundaries ($\phi_a=0$) and in bulk.
- **$\xi_v$ must scale $V_{\mathrm{src}}$ together with diffusion** (CORRECTED 2026-07-08; the earlier claim here that the source must stay unscaled was wrong). M&F 2024 eq. 26 writes the vapor equation as $(1/\xi_v)\,\partial(\phi_a\rho_v)/\partial t = \nabla\cdot(\phi_a D_v \nabla\rho_v) + \rho_i\,\partial\phi_a/\partial t$ — i.e. $\xi_v$ multiplies the *entire* right-hand side, diffusion **and** mass-exchange source. In the quasi-steady limit $\xi_v$ cancels between the two, so the vapor field and interface velocity remain physical while the fast diffusion timescale is slowed by $1/\xi_v$ (this is the whole point of the temporal scaling: large dt). Scaling only diffusion breaks that cancellation by $1/\xi_v$ (1000× at $\xi_v=10^{-3}$): the interface releases vapor at the full physical rate into a field that can only transport it at $\xi_v$ times the physical rate, producing supersaturation spikes, oscillation of the shrinking grain, and Newton failure. The conservation cost of scaling the source is an apparent ice→vapor transfer deficit of $(1-\xi_v)\times$ the vapor-mass change — but total vapor mass is $\sim\rho_v/\rho_i \sim 10^{-6}$ of ice mass, the "noticeable errors" M&F explicitly deem negligible. A raw ice+vapor `TOTAL_MASS` monitor will show this bookkeeping drift by construction; the correct conservation diagnostics under temporal scaling are a $\xi_v=1$ reference run or comparing ice-mass change against the $\xi_v$-scaled source integral.

The earlier $2\phi_a\,\partial\phi_i/\partial t$ form was introduced to *conserve* the AC mass change (branch log §26→ later iterations); it is now retired because that AC mass change is non-physical and should not be conserved by manufacturing vapor. The right levers to limit AC mass drift are reducing the mobility $M$ (`-mob_sub`) and enabling Gibbs-Thomson (§4), not coupling AC motion to vapor (branch log §28).

**Effective diffusivity (penalty-weighted linear combination).** The vapor diffusivity transitions smoothly between full physical $D_v$ in air and a reduced "penalised" value in *solid phases* (ice or sediment):

$$D_{\mathrm{eff}} \;=\; D_{\mathrm{pen}}\, G_{\mathrm{pen}}(\phi_i + \phi_s) \;+\; D_v\,\bigl[1 - G_{\mathrm{pen}}(\phi_i + \phi_s)\bigr], \qquad D_{\mathrm{pen}} = \alpha_{\mathrm{pen}} D_v$$

with $G_{\mathrm{pen}}$ the shifted smooth Heaviside defined below and $\alpha_{\mathrm{pen}}$ a dimensionless penalty factor (`-difvap_pen`). The switch argument $(\phi_i+\phi_s)$ puts the penalty on **both ice and sed**, leaving bulk air and the entire diffuse interface at full physical diffusivity. **Setting $\alpha_{\mathrm{pen}} = 1$ collapses this to $D_{\mathrm{eff}} = D_v$ everywhere** (no diffusivity penalty).

**Interface equilibrium penalty.** The term $\xi_v\,k_{\mathrm{pen}}\,G_{\mathrm{pen}}(\phi_i+\phi_s)\,(\rho_v - \rho_v^{\mathrm{eq}})$ drives $\rho_v$ toward a constant anchor in deep solid, but only where the diffuse interface has fully resolved. **Setting $k_{\mathrm{pen}} = 0$ disables this term entirely.** $G_{\mathrm{pen}}$ is a *shifted* smooth Heaviside concentrated near the bulk-solid edge:

$$G_{\mathrm{pen}}(\phi) = \begin{cases} 0 & \phi \le \phi_{\mathrm{lo}} \\ H\!\left(\dfrac{\phi - \phi_{\mathrm{lo}}}{\phi_{\mathrm{hi}} - \phi_{\mathrm{lo}}}\right) & \phi_{\mathrm{lo}} < \phi < \phi_{\mathrm{hi}} \\ 1 & \phi \ge \phi_{\mathrm{hi}} \end{cases}$$

with $\phi_{\mathrm{lo}} = 0.90$, $\phi_{\mathrm{hi}} = 1.00$ (`PenaltyWeight()` in `src/material_properties.c`) and $H(\xi) = \xi^3(3-2\xi)$ the smooth cubic Heaviside. The penalty is **off** through the entire diffuse interface (so curvature-driven phenomena emerge freely) and **on** only in the deep-solid edge of the band.

The reference anchor is a **constant** (in $\phi$):

$$\rho_v^{\mathrm{eq}}(T) \;=\; h_0\,\rho_{vs}(T)$$

where $h_0$ is the initial humidity (`-humidity`) and $\rho_{vs}(T)$ is the flat-interface saturation vapor density over ice (§5.3). An earlier formulation used $\rho_v^{\mathrm{eq}} = (\phi_i+\phi_s)\rho_{vs} + \phi_a \rho_v$, which was self-referential in air and drained vapor from air into solid along the diffuse band; the constant anchor removes that pathology (branch log §27).

**Air fraction limit.** $\phi_a^{\mathrm{eff}} = \max(\phi_a, \phi_a^{\mathrm{lim}})$ with $\phi_a^{\mathrm{lim}} = 10^{-6}$ prevents numerical singularity in pure solid where $\phi_a \to 0$.

**Current recommendation: both penalties OFF.** Once the §3.3 latent-heat pairing fix, the §3.4 mass-conserving $V_{\mathrm{src}}$, the §3.4 constant-anchor $\rho_v^{\mathrm{eq}}$, and the §3.4 diffusivity-penalty direction (penalty on **solid**, not on air) were all in place, the two penalties had nothing physical left to do. Diagnostic runs on `2D_separated_grains` at T = −20 °C, 30 days show that setting $k_{\mathrm{pen}} = 0$ and $\alpha_{\mathrm{pen}} = 1$ gives **cleaner** behaviour than any prior configuration with the penalties active (branch log §27). The penalty infrastructure is retained in the source code (zero cost when $k_{\mathrm{pen}}=0$) for clean revertability and as the natural place to plug a future physical surface-kinetics term.

---

## 4. Sublimation and Condensation Kinetics

The **Allen-Cahn equation for ice** includes a phase-change source term:

$$S_{\mathrm{sub}} \;=\; \frac{\alpha_{\mathrm{sub}}\,\phi_i^2\,\phi_a^2}{\rho_{\mathrm{ice}}}\,\bigl(\rho_v - \rho_{vs}^{\mathrm{eff}}\bigr)$$

where $\rho_{vs}^{\mathrm{eff}}$ is the Gibbs-Thomson-corrected equilibrium vapor density (below). This term localises the phase change to ice-air interfaces through the $\phi_i^2\phi_a^2$ factor, which peaks at the interface midpoint (where both $\phi_i$ and $\phi_a$ are non-zero) and **vanishes identically at ice-sed boundaries** (where $\phi_a = 0$). It drives:
- **Sublimation** when $\rho_v < \rho_{vs}^{\mathrm{eff}}$: ice decreases, $\phi_a$ increases.
- **Condensation** when $\rho_v > \rho_{vs}^{\mathrm{eff}}$: ice grows, $\phi_a$ decreases.

$S_{\mathrm{sub}}$ enters the ice residual (§3.1) as a negative contribution (subtracted from $R^i$) and the temperature residual (§3.3) as the latent-heat source.

**Mass conservation.** The mass exchange between ice and vapor is enforced by the Stefan-condition source $V_{\mathrm{src}} = -\rho_{\mathrm{ice}}\,S_{\mathrm{sub}}$ on the vapor side (§3.4) and the $-S_{\mathrm{sub}}$ term on the ice side. They are exactly equal and opposite pointwise (both carry the same $S_{\mathrm{sub}}$), so the global balance $\frac{d}{dt}\!\int\rho_{\mathrm{ice}}\phi_i + \frac{d}{dt}\!\int\rho_v = 0$ holds by construction for the physical phase change. Allen-Cahn curvature motion is deliberately *not* coupled to vapor (it is mass-neutral interface rearrangement, a numerical artifact — see §3.4).

**Gibbs-Thomson curvature dependence (active when `d0_GT` > 0).** The local equilibrium vapor density at a curved ice surface differs from the flat-interface saturation value $\rho_{vs}(T)$ by the Kelvin / Gibbs-Thomson correction:

$$\rho_{vs}^{\mathrm{eff}} \;=\; \rho_{vs}(T)\,\bigl(1 + d_0^{\mathrm{GT}} \,\kappa \bigr), \qquad \kappa \;=\; -\nabla\!\cdot\!\left(\frac{\nabla\phi_i}{|\nabla\phi_i|}\right)$$

where $d_0^{\mathrm{GT}} = \gamma_{iv} v_m / (R_g T)$ is the capillary length (~$9.6 \times 10^{-10}$ m for ice; varies weakly with $T$). This is the Lifshitz–Slyozov–Wagner (LSW) driving force for Ostwald ripening: larger grains (smaller $\kappa$) have lower $\rho_{vs}^{\mathrm{eff}}$ than smaller grains, so vapor flows small → large and grains coarsen indefinitely.

The curvature $\kappa$ is computed from the gradient and Hessian of $\phi_i$ in `Curvature()` (`material_properties.c`); `Residual_A1` and `Jacobian_A1` read the Hessian via `IGAPointFormHess`. The correction enters $S_{\mathrm{sub}}$ through $\rho_{vs}^{\mathrm{eff}}$ and is analytically differentiated in the Jacobian (`[ice,ice]` and `[vap,ice]` GT blocks via $d\kappa/d(\nabla\phi_i)$ and $d\kappa/d(\nabla^2\phi_i)$). Setting `-d0_GT 0` disables it identically ($\rho_{vs}^{\mathrm{eff}} = \rho_{vs}(T)$, flat-interface model). **Without GT the model freezes after the initial Allen-Cahn relaxation** — all grains then share one $\rho_{vs}$ and there is no driving force for ripening (branch log §28). See branch log §24–§28 for the full history.

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

**Solver parameters** (defaults in `inputs/solver.opts`):
- snes_max_it = 50 (maximum Newton iterations)
- snes_rtol = 10⁻⁶ (relative residual tolerance)
- snes_atol = 10⁻⁴ (absolute residual tolerance — aggressively relaxed; the T-equation residual at iteration 0 is ~7×10⁻⁵ in the typical configuration, so atol = 10⁻⁴ means SNES converges immediately and triggers `NRmin`-based dt-growth on every step. See §9.8.)
- snes_stol = 10⁻⁶ (step size tolerance)
- snes_linesearch_type = bt (backtracking)

### 9.6 Per-DOF Convergence Test

A custom convergence test (`SNESDOFConvergence`) checks each DOF independently rather than using the global residual norm. At iteration k, convergence is declared when all four DOFs satisfy at least one criterion:

$$\|\mathbf{r}_j^{(k)}\|_2 \leq r_{\mathrm{tol}}\,\|\mathbf{r}_j^{(0)}\|_2 \quad\text{(relative)}$$
$$\|\mathbf{r}_j^{(k)}\|_2 \leq a_{\mathrm{tol}} \quad\text{(absolute)}$$
$$\|\delta\mathbf{u}_j^{(k)}\|_2 \leq s_{\mathrm{tol}}\,\|\mathbf{u}_j^{(k)}\|_2 \quad\text{(step size)}$$

for j ∈ {0=ice, 1=T, 2=ρ_v, 3=sed}. A minimum of 3 Newton iterations (1 during relaxation) is required before convergence can be declared, preventing premature exit on the first Newton step when the initial residual is artificially small.

### 9.7 Linear Solver: GMRES with Additive Schwarz Preconditioning

Each Newton iteration requires solving the 4 × 4 block linear system J δU = −R. The linear solver stack is:

- **Krylov method:** GMRES (Saad & Schultz 1986)
  - ksp_rtol = 10⁻⁸, ksp_atol = 10⁻¹²
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
| Newton iters < NRmin = 5 | Increase dt by factor $f = 10^{1/8} \approx 1.33$ |
| NRmin ≤ Newton iters ≤ NRmax = 15 | Keep dt unchanged |
| Newton iters > NRmax = 15 | Decrease dt by factor 1/f; reject and retry |
| Rejected steps > max_rej = 25 | Abort simulation |

**Why NRmin = 5 (not the more common 3).** With the relaxed `xi_v = xi_T = 1.0` setting and `snes_atol = 10⁻⁴`, SNES converges in exactly **3** iterations on essentially every accepted step in the steady-cruise regime. With the more common `NRmin = 2`, the condition "Newton iters < NRmin" was never satisfied, so dt never grew above the initial value (~0.01 s). Setting `NRmin = 5` lets 3-iter convergence trigger dt-growth, which unlocks the dt → 10⁴ s adaptive ramp needed for multi-day integrations. This was the single largest performance unlock during the §27 stability-fix campaign.

The time step is hard-capped at `dtmax = 10⁴ s` (~ 2.8 hours), and the adaptive scheme produces naturally log-spaced snapshots when output is set to every step (`-outp 1`).

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
| B-spline degree | p | 2 | `inputs/solver.opts` |
| B-spline continuity | C | 1 (C¹) | `inputs/solver.opts` |
| Interface half-width (standard) | ε | 7.12 × 10⁻⁷ m | `inputs/geometry/*` |
| Interface half-width (hi-res) | ε | 3.56 × 10⁻⁷ m | `inputs/geometry/*_hires.opts` |
| Ice-vapor surface energy | γ_iv | 0.109 J m⁻² | `permafrost2.c` |
| Ice-sediment surface energy | γ_is | 0.033 J m⁻² | `permafrost2.c` |
| Sediment-vapor surface energy | γ_sv | 0.056 J m⁻² | `permafrost2.c` |
| Triple-junction penalty | Λ | 1.0 × 10⁴ | `inputs/solver.opts` |
| Latent heat of sublimation | L_sub | 2.83 × 10⁶ J kg⁻¹ | `permafrost2.c` |
| Air diffusivity floor | φ_a^lim | 10⁻⁶ | `permafrost2.c` |
| Diffusivity penalty factor | α_pen | 10⁻⁸ (default, with penalty) / 1.0 (penalty off) | `inputs/solver.opts` |
| Interface equilibrium stiffness | k_pen | 10³ s⁻¹ (default) / 0 (penalty off; current recommendation, §3.4) | `inputs/solver.opts` |
| Penalty ramp band | [φ_lo, φ_hi] | [0.90, 1.00] | `src/material_properties.c::PenaltyWeight` |
| Vapor time scaling | ξ_v | 1.0 (physical) | `inputs/solver.opts` |
| Thermal time scaling | ξ_T | 1.0 (physical) | `inputs/solver.opts` |
| Capillary length scale (mobility derivation) | d₀⁰ | 10⁻⁹ m | `permafrost2.c` |
| Gibbs-Thomson capillary length | d₀^GT | 0 (off) / 9.6×10⁻¹⁰ m (physical, §4) | `inputs/experiment/*.opts`, `-d0_GT` |
| Kinetic coefficient | β⁰ | 1.4 × 10⁵ s m⁻¹ | `permafrost2.c` |
| Relaxation steps | n_relax | 1 | `inputs/solver.opts` |
| Sed. freeze time | t_sed_freeze | 10 s | `inputs/experiment/*.opts` |
| Adaptive-dt growth threshold | NRmin | 5 (Newton iters strictly *below* this trigger dt growth) | `inputs/solver.opts` |
| Adaptive-dt shrink threshold | NRmax | 15 | `inputs/solver.opts` |
| dt cap | dtmax | 10⁴ s (~2.8 h) | `inputs/solver.opts` |

---

## 13. Why the Current Model Eliminates Spurious Air

> **Historical note.** This section documents the early-2026 cascade of
> bugs that produced the spurious-air-in-ice failure on the 2-phase
> formulation (run dated 2026-05-05). The fixes described here address
> that *specific* failure. Several *further* model-level fixes were
> made later on the `fix/spurious-ice-sed-air-penalty` branch — the
> 3-phase ice equation under frozen sed, the diffusivity-penalty
> direction inversion, the `vap_src ↔ sub_src` pairing, the
> Gibbs-Thomson curvature dependence — and are documented narratively in
> [`spurious_ice_sed_air_branch_log.md`](spurious_ice_sed_air_branch_log.md).
> The formal model equations in §3–§5 reflect the *current* code
> including those later fixes.

A detailed comparison of the old run (`test_1D_IceSlab_2Phase_difvappen1e-07_k_pen1e09`, 2026-05-05) with the current code reveals **eight distinct errors**, each of which would individually degrade simulation quality, and which together created a catastrophic failure cascade that drove φ_i < 0 (spurious air) throughout the domain.

### 13.0 Summary of All Changes

| Parameter / Code Element | Old | Current | Impact |
|--------------------------|-----|---------|--------|
| xi_v on k_pen in vapor residual | **absent** (bug) | `xi_v * k_pen` | Dominant: penalty 1000× too large |
| Initial ρ_v inside solid | `h₀ ρ_vs` (wrong) | `ρ_vs` (correct) | Enormous t=0 disequilibrium |
| Jacobian | FD (IGAFormIJacobianFD) | Analytical (Jacobian_A1) | O(10) errors in dominant blocks |
| k_pen | 10⁹ | 10⁵ | 4 orders of magnitude |
| ξ_v | 10⁻³ | 10⁻⁴ | Tighter time-scale separation |
| difvap_pen (α_pen) | 10⁻⁷ | 10⁻⁴ | Near-discontinuity → smooth |
| humidity (h₀) | 0.50 | 0.95 | 10× smaller initial driving force |
| n_relax | 12 | 1 | Vapor/ice coupled from step 1 |
| t_sed_freeze | 0 s | 1 s | Brief 3-phase equilibration |
| snes_atol | 10⁻⁸ | 10⁻¹⁰ | False convergence at small residual |
| ksp_rtol | 10⁻⁵ | 10⁻⁶ | Linear solve accuracy |
| Preconditioner | bjacobi + ILU(1) | ASM + ILU(2), overlap=1 | Better conditioning |

---

### 13.1 PRIMARY BUG: Missing ξ_v on the Vapor Penalty Residual

This is the single most important error. In `src/assembly.c`, the vapor residual in the old code was:

```c
/* OLD — INCORRECT */
R_vap  = N0[a] * rhov_t;                                              // time deriv (unscaled)
R_vap += xi_v * difvap * air_eff * (N1[a][l] * grad_rhov[l]);        // diffusion (xi_v present)
R_vap += k_pen * g_phiiphis * (rhov - rhov_eq) * N0[a];               // penalty (NO xi_v) ← BUG
R_vap -= xi_v * N0[a] * rho_ice * air_t;                              // source (xi_v present)
```

The current code fixes this as:

```c
/* NEW — CORRECT */
PetscReal vap_pen = xi_v * k_pen * g_phiiphis * (rhov - rhov_eq);    // xi_v present
PetscReal vap_src = xi_v * rho_ice * air_t;
R_vap  = N0[a] * rhov_t;
R_vap += xi_v * difvap * air_eff * (N1[a][l] * grad_rhov[l]);
R_vap += vap_pen * N0[a];                                              // xi_v present ← FIXED
R_vap -= vap_src * N0[a];
```

**Why this matters — quantitative analysis.**

The vapor PDE has the physical form:

$$\frac{\partial\rho_v}{\partial t} = \xi_v\left[-\nabla\cdot(D_{\mathrm{eff}}\nabla\rho_v) - k_{\mathrm{pen}} H(\phi_i+\phi_s)(\rho_v - \rho_{\mathrm{eq}}) + \rho_{\mathrm{ice}}\frac{\partial\phi_a}{\partial t}\right]$$

Every spatial term carries the factor ξ_v, which separates the vapor time scale from the phase-field time scale. The penalty should be ξ_v × k_pen. In the old code, the penalty entered without ξ_v, making it effectively ξ_v⁻¹ times larger than intended.

With the old parameters (ξ_v = 10⁻³, k_pen = 10⁹):

| Term | Old magnitude | New magnitude (ξ_v = 10⁻⁴, k_pen = 10⁵) |
|------|--------------|------------------------------------------|
| Time derivative ∂ρ_v/∂t | ~ ρ_vs / Δt ≈ 3×10⁻⁵ / 10⁻⁴ = 0.3 | ~ 0.3 |
| Diffusion ξ_v D ρ_v / ε² | ~ 10⁻³ × 2.2×10⁻⁵ × 3×10⁻⁵ / (7×10⁻⁷)² ≈ 1.4×10⁻³ | ~ 10⁻⁴ × same ≈ 1.4×10⁻⁴ |
| Penalty (as written in old code) | **k_pen × δρ_v = 10⁹ × 3×10⁻⁵ ≈ 3×10⁴** | ξ_v × k_pen × δρ_v = 10⁻⁴ × 10⁵ × 3×10⁻⁵ ≈ 3×10⁻⁴ |

The old penalty (3×10⁴) exceeded the time derivative (0.3) by **five orders of magnitude**. The vapor equation was not a PDE — it was a nearly algebraic constraint ρ_v = ρ_eq throughout the solid. The time derivative term was invisible. Any perturbation in φ_i instantly forced ρ_v to follow ρ_eq, regardless of physical plausibility.

In the current model, the penalty (3×10⁻⁴) is smaller than the time derivative (0.3) by three orders of magnitude — the vapor field evolves dynamically with a penalty that gently nudges it toward equilibrium.

---

### 13.2 Wrong Vapor Initial Condition Inside Solid

In `src/initial_conditions.c`, the old code set:

```c
/* OLD — INCORRECT for solid-phase points */
u[j][i].rhov = user->hum0 * rho_vs;   /* applies h₀ everywhere, including inside ice */
```

The current code uses:

```c
/* NEW — CORRECT */
u[j][i].rhov = rho_vs * (user->hum0 * _phi_air + (1.0 - _phi_air));
/* = rho_vs * (h₀ φ_a + φ_i + φ_s)
   → rho_vs   inside solid (φ_a ≈ 0)
   → h₀ rho_vs in air (φ_a ≈ 1) */
```

**Why this matters.** The equilibrium vapor pressure inside a solid is ρ_vs (saturation). Setting ρ_v = h₀ ρ_vs inside the solid at t=0 means the initial condition is nowhere near equilibrium inside the ice. With h₀ = 0.5, the initial disequilibrium inside ice was:

$$\rho_v - \rho_{\mathrm{eq}} = 0.5\,\rho_{vs} - \rho_{vs} = -0.5\,\rho_{vs} \approx -1.5 \times 10^{-5}\ \text{kg m}^{-3}$$

Combined with the unscaled k_pen = 10⁹, the initial penalty residual inside the solid was:

$$R_{\mathrm{vap}}^{(0)} = k_{\mathrm{pen}} \times (\rho_v - \rho_{\mathrm{eq}}) \approx 10^9 \times (-1.5 \times 10^{-5}) \approx -15\ \text{kg m}^{-3}\text{s}^{-1}$$

This residual is enormous compared to the time derivative (≈ 0.3 kg m⁻³ s⁻¹). Newton's first step would try to reduce this 50× excess by adjusting ρ_v upward by ≈ 1.5×10⁻⁵ kg m⁻³, but with an inaccurate FD Jacobian, the step overshot dramatically. The resulting corrupted ρ_v field then drove the ice Allen-Cahn equation via the sublimation source term, pushing φ_i below zero at t = Δt.

---

### 13.3 Finite-Difference Jacobian Replaced by Analytical Jacobian

**Old (`permafrost2.c`):**
```c
// ierr = IGASetFormIJacobian(iga, Jacobian, &user); CHKERRQ(ierr);
ierr = IGASetFormIJacobian(iga, IGAFormIJacobianFD, &user); CHKERRQ(ierr);
```

**Current:**
```c
ierr = IGASetFormIJacobian(iga, Jacobian, &user); CHKERRQ(ierr);
// ierr = IGASetFormIJacobian(iga, IGAFormIJacobianFD, &user); CHKERRQ(ierr);
```

PETSc's finite-difference Jacobian perturbs each DOF by δ ∼ √ε_machine ≈ 10⁻⁸ and differentiates the residual numerically. For residual entries that scale with k_pen = 10⁹, the absolute Jacobian error is:

$$\Delta J \sim \frac{\varepsilon_{\mathrm{machine}}}{\delta} \cdot k_{\mathrm{pen}} \sim \frac{10^{-16}}{10^{-8}} \cdot 10^9 = O(10)$$

This is comparable in magnitude to the actual residual (≈ 0.3–15 depending on location). Newton's convergence theory requires the Jacobian error to be small relative to the residual. When ΔJ ∼ R, the Newton step is qualitatively wrong — the correction may point in the wrong direction entirely. This explains why the old code showed Newton divergence or convergence to unphysical solutions.

The current analytical Jacobian `Jacobian_A1` covers all 16 of the 4×4 cross-coupling blocks, including the critical [ρ_v, φ_i] and [φ_i, ρ_v] off-diagonal blocks (sublimation coupling) and the [ρ_v, ρ_v] diagonal block (penalty + diffusion). With no FD perturbation, Newton converges quadratically from the second iteration onward.

---

### 13.4 Extreme Interface Equilibrium Stiffness (k_pen = 10⁹ → 10⁵)

Even with the xi_v bug fixed, k_pen = 10⁹ creates problems. The Jacobian condition number of the vapor block scales as:

$$\kappa \sim \frac{\xi_v k_{\mathrm{pen}} + \xi_v D / \varepsilon^2}{1/\Delta t} \sim \xi_v k_{\mathrm{pen}} \Delta t$$

With ξ_v = 10⁻³ and k_pen = 10⁹, κ ∼ 10⁶ × Δt. For Δt = 10⁻⁴ s, κ ∼ 100 — borderline acceptable. But with the bug present (no ξ_v on k_pen), κ ∼ k_pen × Δt = 10⁹ × 10⁻⁴ = 10⁵, which is beyond the effective range of ILU(1) preconditioning and GMRES with restart=500.

The current k_pen = 10⁵ gives ξ_v × k_pen = 10⁻⁴ × 10⁵ = 10 — three orders of magnitude smaller than diffusion at the interface scale, making the linear system well-conditioned for iterative solvers.

---

### 13.5 Vapor Time-Scale Separation (ξ_v = 10⁻³ → 10⁻⁴)

**Old:** ξ_v = 10⁻³ in `permafrost2.c`.  
**Current:** ξ_v = 10⁻⁴.

The ξ_v parameter controls the ratio of the vapor time scale to the phase-field time scale. A smaller ξ_v means vapor equilibrates more slowly relative to the phase field, giving the nonlinear solver more room to adjust both fields simultaneously without stiffness. The factor-of-10 reduction also means that the correctly-scaled penalty term ξ_v × k_pen is 10× smaller, further reducing the condition number of the vapor block.

---

### 13.6 Near-Discontinuous Diffusivity (difvap_pen = 10⁻⁷ → 10⁻⁴)

**Old:** α_pen = 10⁻⁷ — seven orders of magnitude jump in D_eff from air to ice.  
**Current:** α_pen = 10⁻⁴ — four orders of magnitude.

The smooth Heaviside H(φ_i) spans ~4ε ≈ 3×10⁻⁶ m. With α_pen = 10⁻⁷, adjacent Gauss points within this band see D_eff values differing by up to 10⁷. Even the analytical Jacobian (let alone FD) cannot accurately represent such a sharp coefficient variation on a B-spline mesh with only 3 quadrature points across the interface. The result is a poorly-resolved diffusion term that generates spurious ρ_v gradients across the interface. With α_pen = 10⁻⁴, the jump is physically motivated (ice is still essentially impermeable to vapor) but numerically resolvable on the current mesh.

---

### 13.7 Undersaturated Initial Condition for Vapor (humidity = 0.50 → 0.95)

**Old:** h₀ = 0.5, meaning ρ_v = 0.5 ρ_vs in the air phase at t=0.  
**Current:** h₀ = 0.95.

Lowering h₀ amplifies every other error. The sublimation driving force at each ice-air interface is proportional to (1 − h₀) ρ_vs. With h₀ = 0.5, the driving force is 10× larger than with h₀ = 0.95. Every Jacobian error, every convergence tolerance issue, every penalty imbalance is amplified by the same factor. Running near-equilibrium (h₀ = 0.95) is not just a physical convenience — it dramatically improves the effective tolerance of the nonlinear solver by shrinking the amplitude of all problematic source terms.

---

### 13.8 Excessive Allen-Cahn Relaxation (n_relax = 12 → 1) and Sediment Freeze Timing (t_sed_freeze = 0 → 1 s)

**Old:** 12 Allen-Cahn relaxation steps before the first full physics step; sediment frozen immediately (t_sed_freeze = 0).  
**Current:** 1 relaxation step; sediment frozen after 1 s.

With n_relax = 12, the phase fields φ_i and φ_s underwent 12 full AC updates without any coupling to the vapor equation. This allowed the phase-field profiles to evolve into positions that were geometrically reasonable (smooth interfaces) but thermodynamically inconsistent with the current ρ_v field. When the first full-physics step then tried to satisfy all four equations simultaneously, the vapor equation saw large disequilibria (created by the 12 AC steps) and the stiff penalty responded with enormous residuals.

With n_relax = 1, there is only a single AC equilibration step before full coupling begins. The initial disequilibrium is minimal, and the solver converges from a consistent starting point.

With t_sed_freeze = 0 in the old code, the three-phase free energy was immediately replaced by the simpler two-phase (sediment-frozen) form without any relaxation of the sediment interface. Any mismatch between the initial φ_s profile and the two-phase energy minimum appeared as a sudden force at t=0. With t_sed_freeze = 1 s, the sediment has one second to equilibrate its interface shape under the full three-phase energy before freezing, eliminating this initial transient.

---

### 13.9 Cascade: How the Errors Amplified Each Other

The errors did not act independently — they formed a destructive cascade:

1. **t=0:** IC set ρ_v = 0.5 ρ_vs inside ice (should be ρ_vs). Penalty residual = 10⁹ × (−0.5 ρ_vs) ≈ −15 inside ice. Time derivative ≈ 0.3. **Penalty exceeds time derivative by 50×.**

2. **FD Jacobian step:** PETSc perturbs ρ_v by 10⁻⁸ to compute ∂R/∂ρ_v. The dominant term is k_pen = 10⁹, so ΔJ ∼ 10⁹ × 10⁻⁸ = 10. But the residual is 15. **Jacobian error is comparable to residual.** The Newton correction is qualitatively wrong.

3. **Newton step:** The solver applies a correction Δρ_v that is too large (the Jacobian underestimates the curvature). ρ_v overshoots ρ_eq on the other side, creating a positive disequilibrium. The penalty now has the opposite sign and pushes ρ_v back — but again, overshoots.

4. **Sublimation coupling:** Each Newton iteration changes ρ_v by O(ρ_vs), which drives the sublimation source term S_sub = alph_sub × (ρ_v − ρ_vs) / ρ_ice. With ρ_v oscillating by O(ρ_vs), S_sub oscillates by O(alph_sub). The Allen-Cahn equation for φ_i receives this as a right-hand side forcing that pushes φ_i in alternating directions over successive Newton iterations.

5. **Divergence or wrong root:** Newton either diverges (SNES reports diverged after max_it) or converges to a φ_i < 0 solution (a spurious local minimum of the 4-DOF residual that satisfies all convergence tolerances but is physically unphysical).

6. **Time accumulation:** Once φ_i < 0 in any element at time t, the next time step starts from a corrupted initial condition. The interface reconstruction at t + Δt is inconsistent, creating new regions of spurious air that grow over time.

The current code breaks this cascade at **every link**:
- The IC formula ensures ρ_v = ρ_vs inside solid at t=0 → near-zero initial penalty residual.
- The xi_v factor reduces the penalty by 1000× → penalty comparable to time derivative.
- The analytical Jacobian eliminates FD errors → Newton converges quadratically.
- k_pen = 10⁵ instead of 10⁹ → condition number 10⁴× lower.
- h₀ = 0.95 → driving forces 10× smaller throughout.
- n_relax = 1 → no pre-relaxation that creates vapor-field inconsistency.

---

## 14. LaTeX source for the model equations

The block below is a self-contained LaTeX source for the model equations,
ready to paste into a manuscript. It assumes `amsmath` and uses standard
notation throughout. Copy from `\begin{align}` to `\end{align}` for individual
equations, or use the whole block as a model-description section.

```latex
% --------------------------------------------------------------------------
% Permafrost three-phase Allen-Cahn + Stefan-coupled sublimation model.
% Notation:
%   phi_i, phi_s, phi_a : phase fields for ice, sediment, air (phi_a = 1 - phi_i - phi_s)
%   T, rho_v            : temperature and vapor mass density
%   M                   : Allen-Cahn mobility (m^3/s)
%   epsilon             : diffuse-interface half-width (m)
%   eta_k, Lambda       : phase-specific double-well coefficients and triple-junction penalty
%   E_T                 : eta_i*eta_s + eta_i*eta_a + eta_s*eta_a
%   f_i, f_s, f_a       : bulk driving-force functions (defined inline)
%   rho, c_p, lambda    : volume-fraction-weighted mixture density, heat capacity, conductivity
%   L_sub               : latent heat of sublimation
%   D_v(T)              : vapor diffusivity in air (power law in T)
%   rho_vs(T)           : flat-interface saturation vapor density over ice
%   alpha_sub           : sublimation kinetic coefficient (derived; see §4)
%   xi_T, xi_v          : optional time-scale separation factors (both 1.0 by default)
%   k_pen, alpha_pen    : numerical penalty stiffness and diffusivity factor (both 0 / 1 by default)
%   G_pen(.)            : shifted smooth Heaviside on [phi_lo, phi_hi] = [0.90, 1.00]
%   h_0                 : initial relative humidity
%   rho_ice             : ice mass density
% --------------------------------------------------------------------------

% --- Free-energy / driving-force functions -----------------------------
\begin{align}
  \Psi(\phi_i,\phi_s,\phi_a)
    &= \sum_{k\in\{i,s,a\}} \eta_k\,\phi_k^2(1-\phi_k)^2
       \;+\; \Lambda\,\phi_i^2\phi_s^2\phi_a^2, \\[4pt]
  f_i(\phi_i,\phi_s)
    &= \eta_i\phi_i(1-\phi_i)(1-2\phi_i)
       + 2\Lambda\,\phi_i\phi_s^2\phi_a^2,\\
  f_s(\phi_i,\phi_s)
    &= \eta_s\phi_s(1-\phi_s)(1-2\phi_s)
       + 2\Lambda\,\phi_i^2\phi_s\phi_a^2,\\
  f_a(\phi_i,\phi_s)
    &= \eta_a\phi_a(1-\phi_a)(1-2\phi_a)
       + 2\Lambda\,\phi_i^2\phi_s^2\phi_a, \\[4pt]
  E_T &= \eta_i\eta_s + \eta_i\eta_a + \eta_s\eta_a.
\end{align}

% --- Ice phase: Allen-Cahn with Stefan source --------------------------
\begin{equation}
  \frac{\partial \phi_i}{\partial t}
    \;=\; M\!\left[\,
       \varepsilon^{2}\nabla^{2}\phi_i
       \;-\; \frac{1}{\varepsilon\,E_T}\bigl[(\eta_s+\eta_a)\,f_i
                                            - \eta_a\,f_s
                                            - \eta_s\,f_a\bigr]
    \right]
    \;+\; S_{\mathrm{sub}}.
  \label{eq:ice-AC}
\end{equation}

% --- Sediment phase: Allen-Cahn (frozen after t = t_sed_freeze) --------
\begin{equation}
  \frac{\partial \phi_s}{\partial t}
    \;=\; M\!\left[\,
       \varepsilon^{2}\nabla^{2}\phi_s
       \;+\; \frac{1}{\varepsilon\,E_T}\bigl[-\eta_a\,f_i
                                             - \eta_i\,f_a
                                             + (\eta_i+\eta_a)\,f_s\bigr]
    \right],
    \qquad
    \frac{\partial \phi_s}{\partial t} = 0 \;\;\text{for}\;\; t \ge t_{\mathrm{sed,freeze}}.
  \label{eq:sed-AC}
\end{equation}

% --- Temperature: energy balance, latent heat paired with S_sub --------
\begin{equation}
  \rho\,c_p\,\frac{\partial T}{\partial t}
    \;=\; \xi_T\,\nabla\!\cdot(\lambda\,\nabla T)
    \;+\; \xi_T\,\rho\,L_{\mathrm{sub}}\,S_{\mathrm{sub}}.
  \label{eq:temperature}
\end{equation}

% --- Vapor density: diffusion + (optional) penalty + Stefan source -----
\begin{equation}
  \frac{\partial \rho_v}{\partial t}
    \;-\; \nabla\!\cdot\!\bigl(\xi_v\,D_{\mathrm{eff}}\,\phi_a^{\mathrm{eff}}\,\nabla \rho_v\bigr)
    \;+\; \xi_v\,k_{\mathrm{pen}}\,G_{\mathrm{pen}}(\phi_i+\phi_s)\bigl(\rho_v - \rho_v^{\mathrm{eq}}\bigr)
    \;=\; V_{\mathrm{src}},
  \label{eq:vapor}
\end{equation}
with the Stefan-condition source (Moure \& Fu 2024; driven by the phase-change term)
\begin{equation}
  V_{\mathrm{src}} \;=\; -\,\rho_{\mathrm{ice}}\,S_{\mathrm{sub}}
                 \;=\; -\,\alpha_{\mathrm{sub}}\,\phi_i^{2}\,\phi_a^{2}\,
                       \bigl(\rho_v - \rho_{vs}^{\mathrm{eff}}\bigr),
  \label{eq:vap-src}
\end{equation}
and constant equilibrium anchor
\begin{equation}
  \rho_v^{\mathrm{eq}}(T) \;=\; h_0\,\rho_{vs}(T).
  \label{eq:rhov-eq}
\end{equation}

% --- Sublimation source: localised to ice-air diffuse band -------------
\begin{equation}
  S_{\mathrm{sub}}
    \;=\; \frac{\alpha_{\mathrm{sub}}\,\phi_i^{2}\,\phi_a^{2}}{\rho_{\mathrm{ice}}}\,
          \bigl(\rho_v - \rho_{vs}^{\mathrm{eff}}\bigr),
  \qquad
  \rho_{vs}^{\mathrm{eff}} \;=\; \rho_{vs}(T)\,\bigl(1 + d_0^{\mathrm{GT}}\,\kappa\bigr),
  \quad
  \kappa = -\nabla\!\cdot\!\frac{\nabla\phi_i}{|\nabla\phi_i|}.
  \label{eq:S-sub}
\end{equation}

% --- Effective vapor diffusivity (penalty-weighted) --------------------
\begin{equation}
  D_{\mathrm{eff}}(\phi)
    \;=\; D_{\mathrm{pen}}\,G_{\mathrm{pen}}(\phi_i + \phi_s)
       \;+\; D_v(T)\,\bigl[1 - G_{\mathrm{pen}}(\phi_i + \phi_s)\bigr],
  \qquad
  D_{\mathrm{pen}} \;=\; \alpha_{\mathrm{pen}}\,D_v(T).
  \label{eq:D-eff}
\end{equation}

% --- Shifted smooth Heaviside used in both penalty terms ---------------
\begin{equation}
  G_{\mathrm{pen}}(\phi) =
  \begin{cases}
    0 & \phi \le \phi_{\mathrm{lo}}, \\[2pt]
    H\!\left(\dfrac{\phi - \phi_{\mathrm{lo}}}{\phi_{\mathrm{hi}} - \phi_{\mathrm{lo}}}\right)
      & \phi_{\mathrm{lo}} < \phi < \phi_{\mathrm{hi}}, \\[6pt]
    1 & \phi \ge \phi_{\mathrm{hi}},
  \end{cases}
  \qquad H(\xi) = \xi^{2}(3 - 2\xi).
  \label{eq:G-pen}
\end{equation}

% --- Air-fraction floor ------------------------------------------------
\begin{equation}
  \phi_a^{\mathrm{eff}} \;=\; \max\!\bigl(\phi_a,\,\phi_a^{\mathrm{lim}}\bigr),
  \qquad \phi_a^{\mathrm{lim}} = 10^{-6}.
\end{equation}

% --- Default-configuration note (current best practice) ----------------
% The two numerical penalties in (\ref{eq:vapor}) and (\ref{eq:D-eff}) are
% retained in the code for revertability but are currently disabled by
% default:  k_pen = 0  and  alpha_pen = 1  (so D_eff = D_v(T) everywhere).
% This collapses (\ref{eq:vapor}) to a plain advection-diffusion equation
% with the Stefan source V_src as the only coupling to the phase field.
```

The `// removed` style comments at the bottom of the LaTeX block are
informational; remove them before submission.

---

## 14a. Pointwise mass balance (manuscript-grade derivation)

The Stefan source $V_{\mathrm{src}}$ in equation (\ref{eq:vap-src}) is constructed
so that the global mass balance

$$
  \frac{d}{dt}\!\int_{\Omega} \rho_{\mathrm{ice}}\,\phi_i\,d\Omega
  \;+\; \frac{d}{dt}\!\int_{\Omega} \rho_v\,d\Omega
  \;=\; 0
$$

holds pointwise to leading order in the diffuse-interface width $\varepsilon$.
Inserting (\ref{eq:vap-src}) into (\ref{eq:vapor}) and noting that the
diffusive flux and the penalty integrate to zero (the latter under the
recommendation $k_{\mathrm{pen}} = 0$), the vapor equation reads

$$
  \frac{\partial \rho_v}{\partial t}
    \;=\; -\,2\,\rho_{\mathrm{ice}}\,\phi_a\,\frac{\partial \phi_i}{\partial t}.
$$

Integrating $V_{\mathrm{src}}$ over the diffuse interface, $\phi_a$ varies
from 0 to 1 across a band of width $\sim\varepsilon$ while
$\partial\phi_i/\partial t$ is concentrated in that band, and the factor 2
normalises the integral $\int_0^1 2\phi_a\,d\phi_a = 1$, giving
$\int V_{\mathrm{src}}\,d\Omega = -\rho_{\mathrm{ice}}\,\frac{d}{dt}\!\int\phi_i\,d\Omega$
to leading order in $\varepsilon$. This is the diffuse-interface form of the
Stefan condition and closes the mass balance exactly.

---

## 15. References

- Balay, S. et al. (2023). *PETSc/TAO User's Manual*. Argonne National Laboratory.
- Dalcin, L. et al. (2016). PetIGA: A framework for high-performance isogeometric analysis. *Comput. Methods Appl. Mech. Eng.* 308, 151–181.
- Kim, S.G., Kim, W.T., Suzuki, T. (1999). Phase-field model for binary alloys. *Phys. Rev. E* 60, 7186.
- Saad, Y., Schultz, M.H. (1986). GMRES: A generalized minimal residual algorithm for solving nonsymmetric linear systems. *SIAM J. Sci. Stat. Comput.* 7, 856–869.
- Steinbach, I., Pezzolla, F. (1999). A generalized field method for multiphase transformations using interface fields. *Physica D* 134, 385–393.
