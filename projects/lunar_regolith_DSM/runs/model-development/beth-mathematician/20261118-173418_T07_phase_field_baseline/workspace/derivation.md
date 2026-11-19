# T07_phase_field

Run `20261118-173418_T07_phase_field_baseline` · status **done** · written by `scripts/agentkit/render_derivation.py` from `answer.json`; do not edit by hand.

## Work order (verbatim)

> We are modeling the pore space of a granular material (lunar regolith) in which water ice sublimates and deposits. Treat the pore space as a three-dimensional region filled with solid ice and gas (air that contains water vapor). The local volume fraction of ice is the phase field phi(x, y, z, t), between 0 and 1: phi = 1 in ice and phi = 0 in gas. At phase coexistence the ice-gas interface is described by the free energy F[phi] given below, which you derived in an earlier task. Away from coexistence, ice deposits where the vapor density rho_v exceeds the saturation vapor density rho_vs(T) at the local temperature, and sublimates where rho_v is below it. There is no nucleation: ice grows or shrinks only where ice already meets gas. In this task rho_v and T are known fields; their equations are separate tasks. Derive the equation for how phi evolves in time, combining the relaxation of the interface toward the equilibrium described by F with the growth or sublimation driven by the departure of rho_v from rho_vs(T). State every assumption, and list every effect you leave out with the reason. Keep the coefficients general where you can, and say what each one means physically. Write the result as R = 0 in 3D Cartesian coordinates and return R.

## Symbols

| Symbol | Meaning | Units | Typical value (SI) | Source |
|---|---|---|---|---|
| $x$ | first position coordinate | $\mathrm{m}$ |  | work order |
| $y$ | second position coordinate | $\mathrm{m}$ |  | work order |
| $z$ | third position coordinate | $\mathrm{m}$ |  | work order |
| $t$ | time | $\mathrm{s}$ |  | work order |
| $\phi$ | local volume fraction of ice (phase field), between 0 and 1 | 1 |  | work order |
| $T$ | temperature | $\mathrm{K}$ |  | work order |
| $\rho_v$ | vapor density: mass of water vapor per unit volume of gas | $\mathrm{kg\,m^{-3}}$ |  | work order |
| $\rho_{vs}(T)$ | saturation vapor density over a flat ice surface at temperature T | $\mathrm{kg\,m^{-3}}$ |  | work order |
| $\gamma$ | energy per unit area of the ice-gas interface, constant | $\mathrm{J\,m^{-2}}$ |  | work order |
| $\varepsilon$ | the length that sets the width of the diffuse interface, constant | $\mathrm{m}$ |  | work order |
| $M(T)$ | Interface mobility. It sets how fast $\phi$ relaxes toward a minimum of $F$, and so how fast curvature (Gibbs–Thomson) moves the interface. A positive kinetic coefficient. | $\mathrm{m^{3}\,J^{-1}\,s^{-1}}$ | 4e-10 | introduced |
| $\alpha(T)$ | Deposition and sublimation rate coefficient: the rate of change of $\phi$ at the interface per unit vapor excess $\rho_v-\rho_{vs}$. It is related to the Hertz–Knudsen kinetic velocity $\beta$ by $\alpha=\beta/(6\varepsilon\rho_i)$. | $\mathrm{m^{3}\,kg^{-1}\,s^{-1}}$ | 2e+04 | introduced |
| $R$ | Residual of the phase-field equation | $\mathrm{s^{-1}}$ | 1 | introduced |
| $\beta$ | Kinetic velocity of the Hertz–Knudsen law, $\alpha_s\sqrt{k_B T/(2\pi m_w)}$, where $\alpha_s$ is the sticking coefficient. Used only to interpret $\alpha$. | $\mathrm{m\,s^{-1}}$ | 100 | introduced |
| $\rho_i$ | Density of ice | $\mathrm{kg\,m^{-3}}$ | 920 | introduced |
| $v$ | Normal speed of a planar front | $\mathrm{m\,s^{-1}}$ | 1e-06 | introduced |
| $\phi_0(\xi)$ | Equilibrium profile of a flat interface | 1 | 0.5 | introduced |
| $\xi$ | Coordinate normal to a flat interface | $\mathrm{m}$ | 1e-06 | introduced |

## Given

- `F = Integral(3*gamma/eps*phi**2*(1 - phi)**2 + 3*gamma*eps*(Derivative(phi, x)**2 + Derivative(phi, y)**2 + Derivative(phi, z)**2), over the volume)`

## Assumptions

- **A1.** A1: $\rho_v$ and $T$ are known, smooth fields (given).
- **A2.** A2: The flat interface of $F$ has energy $\gamma$ per unit area and a 10–90 % width of $2\ln 9\,\varepsilon$ (given; checked in S2).
- **A3.** A3: The dynamics of $\phi$ is non-conserved, linear gradient flow with a scalar, isotropic mobility $M(T)>0$.
- **A4.** A4: The deposition rate is linear in the vapor excess $\rho_v-\rho_{vs}(T)$; this is first order in the supersaturation.
- **A5.** A5: Deposition and sublimation act only where ice meets gas. This is imposed with the shape function $30\phi^2(1-\phi)^2$, which vanishes at $\phi=0$ and $\phi=1$ and integrates to 1 across the interface.

## Derivation

### S1. Variational derivative of F

The relaxation toward equilibrium is driven by $\delta F/\delta\phi$. We differentiate the double-well density and integrate the gradient term by parts. Boundary terms are dropped (no-flux or periodic boundaries).

$$
\frac{\delta F}{\delta\phi} = \frac{6\gamma}{\varepsilon}\phi(1-\phi)(1-2\phi) - 6\gamma\varepsilon\nabla^2\phi \quad [\mathrm{J\,m^{-3}}]
$$

*Uses:* given. *Operation:* Euler–Lagrange variation and integration by parts. *Checked by:* none (no tools in this run; checked by hand: d/dphi of phi^2(1-phi)^2 is 2 phi(1-phi)(1-2phi))

### S2. Check of the flat-interface profile

Setting $\delta F/\delta\phi=0$ in 1D gives the flat profile. With this profile $F$ gives energy $\gamma$ and a 10–90 % width of $2\ln 9\,\varepsilon$, so assumption A2 holds for the given $F$.

$$
\phi_0=\tfrac12\left[1+\tanh\frac{\xi}{2\varepsilon}\right],\quad \int(\phi_0')^2d\xi=\frac{1}{6\varepsilon},\quad 2\cdot 3\gamma\varepsilon\cdot\frac{1}{6\varepsilon}=\gamma,\quad \Delta\xi_{10\text{-}90}=4\varepsilon\,\mathrm{artanh}(0.8)=2\ln 9\,\varepsilon
$$

*Uses:* S1, A2. *Operation:* solve the 1D equilibrium equation and integrate. *Checked by:* none (no tools in this run; by hand: the integral of sech^4 over the real line is 4/3, and artanh(0.8)=ln 3)

### S3. Relaxation of the interface

$\phi$ is not conserved, because ice mass comes from the vapor. So the simplest dissipative dynamics is linear gradient flow (Allen–Cahn) with a positive mobility $M(T)$.

$$
\left.\frac{\partial\phi}{\partial t}\right|_{\mathrm{relax}} = -M(T)\frac{\delta F}{\delta\phi}
$$

*Uses:* S1, A3. *Operation:* linear non-conserved gradient flow. *Checked by:* none (no tools in this run; units by hand: [M][J m^-3] = s^-1)

### S4. Growth and sublimation source

Ice deposits when $\rho_v>\rho_{vs}$ and sublimates when $\rho_v<\rho_{vs}$. To first order the source is linear in $\rho_v-\rho_{vs}$. The shape function $30\phi^2(1-\phi)^2$ vanishes in pure gas and pure ice, so there is no nucleation. It integrates to 1 across the interface and does not shift the stable states $\phi=0$ and $\phi=1$.

$$
\left.\frac{\partial\phi}{\partial t}\right|_{\mathrm{dep}} = 30\,\alpha(T)\,\phi^2(1-\phi)^2\left(\rho_v-\rho_{vs}(T)\right)
$$

*Uses:* A4, A5, A1. *Operation:* constitutive choice for the source term. *Checked by:* none (no tools in this run; units by hand: m^3 kg^-1 s^-1 times kg m^-3 = s^-1)

### S5. Combined equation

Add S3 and S4 and write the Laplacian in Cartesian coordinates. Moving every term to one side gives the residual $R$.

$$
R = \partial_t\phi - M(T)\left[6\gamma\varepsilon(\phi_{xx}+\phi_{yy}+\phi_{zz}) - \frac{6\gamma}{\varepsilon}\phi(1-\phi)(1-2\phi)\right] - 30\alpha(T)\phi^2(1-\phi)^2(\rho_v-\rho_{vs}(T)) = 0
$$

*Uses:* S3, S4, S1. *Operation:* sum of S3 and S4, then substitution of S1. *Checked by:* none (no tools in this run; every term has units of s^-1)

### S6. Physical meaning of alpha

Put a planar front $\phi_0(\xi-vt)$ into S5, multiply by $\phi_0'$ and integrate over $\xi$. The $F$ terms vanish, and $\int 30\phi_0^2(1-\phi_0)^2\phi_0'\,d\xi=1$. This gives the front speed, which we match to the Hertz–Knudsen law $v=\beta(\rho_v-\rho_{vs})/\rho_i$.

$$
v = 6\varepsilon\,\alpha(T)\left(\rho_v-\rho_{vs}\right) \;\Rightarrow\; \alpha(T)=\frac{\beta}{6\varepsilon\rho_i}
$$

*Uses:* S5, S2. *Operation:* multiply by $\phi_0'$ and integrate over $\xi$. *Checked by:* none (no tools in this run; by hand: the S4 source integrates to alpha(rho_v-rho_vs), and the time derivative integrates to -v times 1/(6 eps))

## Result

$$
R = \frac{\partial \phi}{\partial t} - M(T)\left[6\gamma\varepsilon\left(\frac{\partial^2\phi}{\partial x^2}+\frac{\partial^2\phi}{\partial y^2}+\frac{\partial^2\phi}{\partial z^2}\right) - \frac{6\gamma}{\varepsilon}\phi(1-\phi)(1-2\phi)\right] - 30\,\alpha(T)\,\phi^2(1-\phi)^2\left(\rho_v-\rho_{vs}(T)\right) = 0
$$

Graded SymPy form:

```
Derivative(phi(x, y, z, t), t) - M(T(x, y, z, t))*(6*gamma*eps*(Derivative(phi(x, y, z, t), (x, 2)) + Derivative(phi(x, y, z, t), (y, 2)) + Derivative(phi(x, y, z, t), (z, 2))) - 6*gamma/eps*phi(x, y, z, t)*(1 - phi(x, y, z, t))*(1 - 2*phi(x, y, z, t))) - 30*alpha(T(x, y, z, t))*phi(x, y, z, t)**2*(1 - phi(x, y, z, t))**2*(rho_v(x, y, z, t) - rho_vs(T(x, y, z, t)))
```

## Neglected effects

| Effect | Why it is left out |
|---|---|
| Nonlinear dependence on supersaturation, for example $\ln(\rho_v/\rho_{vs})$ | Supersaturations in the pores are small, $\|\rho_v/\rho_{vs}-1\|\ll1$. The dropped terms are of order $(\rho_v/\rho_{vs}-1)^2$. |
| Anisotropy of the ice surface energy and of the attachment kinetics | $\gamma$ is constant by the work order. Faceting at the pore scale is a secondary effect for the bulk transport of ice. |
| Thermal noise (fluctuations) | The interface width $\varepsilon\sim\mu$m is much larger than the thermal length $\sqrt{k_BT/\gamma}\sim0.2$ nm. |
| Regolith grains as a third phase, including wetting by ice | The work order describes only ice and gas in the pore space. Grains would enter as boundary conditions or as a fixed field. |
| Latent heat and vapor mass exchange | They couple $\phi$ back to $T$ and $\rho_v$. Those equations are separate tasks; here the fields are given. |
| Elastic stress and mechanical motion of ice | The ice is treated as rigid and fixed in place. Stresses are small compared with the vapor-driven changes. |
| A thermodynamically consistent link between $M$ and $\alpha$ | I kept them as separate coefficients. Thermodynamic consistency would require $\alpha/M=\rho_i k_BT/(m_w\rho_{vs})$. This changes curvature effects only through the ratio $\varepsilon/$(radius of curvature). |

## Questions

- Do you want $M$ and $\alpha$ tied by $\alpha/M=\rho_i k_B T/(m_w\rho_{vs})$? This makes the curvature-driven kinetics consistent with the Gibbs–Thomson effect and the Hertz–Knudsen law. I kept them independent.
- Should the driving force be written in the normalized form $(\rho_v-\rho_{vs})/\rho_{vs}$? Then $\alpha$ would be in s$^{-1}$; the physics does not change.

