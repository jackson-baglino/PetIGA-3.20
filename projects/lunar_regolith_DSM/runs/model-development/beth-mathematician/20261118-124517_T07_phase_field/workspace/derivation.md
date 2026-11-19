# T07 phase field: evolution equation for phi

Checks: `scratch/check.py` (SymPy + pint). Every check printed the expected value.

New symbols: M(phi,T) is the Allen-Cahn mobility [m^3 J^-1 s^-1]. rho_i is the ice density [kg m^-3]. R_v = k_B/m_w is the specific gas constant of water vapor [J kg^-1 K^-1].

**S1. Variational derivative of the given F** (Euler-Lagrange operator; the gradient term is integrated by parts, and the boundary term vanishes for natural boundary conditions):

$$\frac{\delta F}{\delta\phi}=\frac{6\gamma}{\varepsilon}\phi(1-\phi)(1-2\phi)-6\gamma\varepsilon\nabla^2\phi \quad [\mathrm{J\,m^{-3}}]$$

Check: the SymPy difference is 0.

**S2. Equilibrium profile** (solve δF/δφ = 0 in 1D):

$$\phi_0=\tfrac12[1+\tanh(X/2\varepsilon)]$$

This profile solves S1. Its interface energy is $\gamma$, and φ0(ε ln 9) = 0.9, so the 10-90 % width is 2 ln 9 ε. This matches the given assumption (checked).

**S3. Driving force.** Moving vapor of mass $\rho_i$ into a unit volume of ice lowers the free energy by the volumetric driving force. This assumes ideal vapor and an incompressible ice phase.

$$\Delta g=\rho_i R_v T\ln(\rho_v/\rho_{vs}(T))\approx\rho_i R_v T\,\frac{\rho_v-\rho_{vs}(T)}{\rho_{vs}(T)}\quad[\mathrm{J\,m^{-3}}]$$

The linearization drops terms of order O(σ²), where σ = (ρ_v − ρ_vs)/ρ_vs.

**S4. Tilted free energy.** Add the bulk term −Δg h(φ), with

$$h=\phi^3(10-15\phi+6\phi^2),\qquad h'=30\phi^2(1-\phi)^2$$

The function h has h(0) = 0, h(1) = 1 and h'(0) = h'(1) = 0 (checked). So the driving force acts only where 0 < φ < 1, that is, at existing interfaces. This gives no nucleation, and bulk ice and bulk gas stay at φ = 1 and φ = 0.

$$G[\phi]=F[\phi]-\int \Delta g\,h(\phi)\,dV$$

**S5. Non-conserved relaxation (Allen-Cahn).** φ is not conserved: the mass change is supplied by the vapor field, which is a separate task. So

$$\partial_t\phi=-M(\phi,T)\,\frac{\delta G}{\delta\phi}$$

**S6. Substitute S1, S3 and S4 into S5:**

$$\partial_t\phi=-M\Big[\frac{6\gamma}{\varepsilon}\phi(1-\phi)(1-2\phi)-6\gamma\varepsilon\nabla^2\phi\Big]+M\,30\phi^2(1-\phi)^2\rho_iR_vT\frac{\rho_v-\rho_{vs}(T)}{\rho_{vs}(T)}$$

Units: every term has units of s^-1 (pint check).

**S7. Residual form** (move everything to the left side):

$$R=\partial_t\phi+M\Big[\frac{6\gamma}{\varepsilon}\phi(1-\phi)(1-2\phi)-6\gamma\varepsilon(\phi_{xx}+\phi_{yy}+\phi_{zz})\Big]-30M\rho_iR_vT\phi^2(1-\phi)^2\frac{\rho_v-\rho_{vs}(T)}{\rho_{vs}(T)}=0$$

**S8. Meaning of M (sharp-interface limit).** Use a traveling wave φ0(X − vt). Multiply S6 by φ0' and integrate over X, using ∫φ0'² dX = 1/(6ε) (checked):

$$v_n=6\varepsilon M(\Delta g-\gamma K)$$

where K is the mean curvature (sum of the principal curvatures). This is a linear kinetic law with Gibbs-Thomson. If the kinetics are Hertz-Knudsen, v_n = α_c (ρ_vs/ρ_i) √(R_vT/2π) σ, which gives

$$M=\frac{\alpha_c\,\rho_{vs}\sqrt{R_vT/2\pi}}{6\varepsilon\rho_i^2R_vT}$$

For α_c = 1, T = 200 K and ε = 1 µm, this is about 5e-10 m^3 J^-1 s^-1.
