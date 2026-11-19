# T07_phase_field

Run `20261118-162816_T07_phase_field` · status **done** · written by `scripts/agentkit/render_derivation.py` from `answer.json`; do not edit by hand.

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
| $M$ | Allen-Cahn interface mobility: how fast the diffuse interface relaxes toward the minimum of F (curvature-driven motion and profile restoration), taken constant | $\mathrm{m^{3}\,J^{-1}\,s^{-1}}$ | 1e-09 | introduced |
| $\alpha$ | deposition/sublimation kinetic coefficient: rate of change of the ice fraction per unit vapor-density excess at the interface; it absorbs 1/\rho_i and the interface kinetics (sticking coefficient, molecular speed, 1/interface width), taken constant | $\mathrm{m^{3}\,kg^{-1}\,s^{-1}}$ | 1e+05 | introduced |
| $f(\phi)$ | bulk double-well density of F, f = 3\gamma/\varepsilon\,\phi^2(1-\phi)^2 | $\mathrm{J\,m^{-3}}$ | 0.1 | introduced |
| $\mu$ | variational derivative \delta F/\delta\phi (chemical potential of the phase field) | $\mathrm{J\,m^{-3}}$ | 0.1 | introduced |
| $g(\phi)$ | interface indicator weighting the source, g = \phi^2(1-\phi)^2 | 1 | 0.06 | introduced |
| $\phi_0(X)$ | equilibrium flat-interface profile along the normal coordinate X | 1 | 0.5 | introduced |
| $X$ | coordinate normal to a flat interface | $\mathrm{m}$ | 1e-07 | introduced |

## Given

- `F = Integral(3*gamma/eps*phi**2*(1 - phi)**2 + 3*gamma*eps*(Derivative(phi, x)**2 + Derivative(phi, y)**2 + Derivative(phi, z)**2), over the volume)`

## Assumptions

- **A1.** A1: $\rho_v$ and $T$ are known, smooth fields (given).
- **A2.** A2: $F$ gives a flat interface with energy $\gamma$ per area and 10-90% width $2\ln 9\,\varepsilon$ (given).
- **A3.** A3: $\phi$ is a non-conserved order parameter. Near coexistence it follows a linear gradient flow (Allen-Cahn) with a constant mobility $M$; boundaries are no-flux or periodic, so the boundary term from the integration by parts drops out.
- **A4.** A4: There is no nucleation, so the phase-change source vanishes where $\phi=0$ or $\phi=1$. It is localized at the interface by $g(\phi)=\phi^2(1-\phi)^2$; this specific shape is a modeling choice.
- **A5.** A5: The growth rate is linear in the supersaturation $\rho_v-\rho_{vs}(T)$, with a constant kinetic coefficient $\alpha$ that absorbs $1/\rho_i$ and the attachment kinetics.
- **A6.** A6: The ice is rigid and at rest; there is no advection of $\phi$.

## Derivation

### S1. Derivative of the double well

The variational derivative of $F$ needs $f'(\phi)$ for the bulk part of the integrand.

$$
f'(\phi) = \frac{6\gamma}{\varepsilon}\phi(1-\phi)(1-2\phi)
$$

*Uses:* given. *Operation:* differentiate $f=\frac{3\gamma}{\varepsilon}\phi^2(1-\phi)^2$ with respect to $\phi$. *Checked by:* `check_s1` in `scratch/check.py`: S1 df/dphi - 6g/e p(1-p)(1-2p) = 0

### S2. Variational derivative of F

Vary $\phi$ in $F$, integrate the gradient term by parts (boundary term dropped for no-flux or periodic boundaries) to get the Euler-Lagrange derivative $\mu$.

$$
\mu = \frac{\delta F}{\delta\phi} = \frac{6\gamma}{\varepsilon}\phi(1-\phi)(1-2\phi) - 6\gamma\varepsilon\left(\partial_{xx}\phi+\partial_{yy}\phi+\partial_{zz}\phi\right)
$$

*Uses:* S1, given. *Operation:* first variation and integration by parts. *Checked by:* `check_s2` in `scratch/check.py`: S2 deltaF/deltaphi - target = 0

### S3. Equilibrium profile and no-nucleation weight

The flat profile $\phi_0=\tfrac12[1+\tanh(X/2\varepsilon)]$ makes $\mu=0$, so relaxation keeps the given interface structure. The phase-change source must vanish in pure ice and pure gas (no nucleation, A4), so it is weighted by $g(\phi)=\phi^2(1-\phi)^2$, which vanishes with zero slope at $\phi=0,1$.

$$
\mu[\phi_0]=0,\quad g(0)=g(1)=g'(0)=g'(1)=0
$$

*Uses:* S2, A4. *Operation:* substitution of the tanh profile; evaluation of $g$ and $g'$ at 0 and 1. *Checked by:* `check_s3` in `scratch/check.py`: S3 deltaF/deltaphi on tanh profile = 0; g(0), g(1), g'(0), g'(1) = [0, 0, 0, 0]

### S4. Evolution law

The non-conserved phase field relaxes down the gradient of $F$ (Allen-Cahn, A3) and is driven by the supersaturation $\rho_v-\rho_{vs}(T)$ at the interface (A5). Positive supersaturation gives $\partial_t\phi>0$ (deposition); with the source off, $dF/dt=-M\int\mu^2\,dV\le 0$.

$$
\partial_t\phi = -M\mu + \alpha\,g(\phi)\left(\rho_v-\rho_{vs}(T)\right)
$$

*Uses:* S2, S3, A3, A5. *Operation:* add the gradient-flow term and the weighted source. *Checked by:* `check_s4` in `scratch/check.py`: S4 source sign at phi=1/2, rho_v-rho_vs=+1: 1

### S5. Residual in 3D Cartesian form

Substitute $\mu$ from S2 and $g$ from S3 into S4 and move everything to one side. Every term has units of s$^{-1}$.

$$
R = \partial_t\phi + M\left[\frac{6\gamma}{\varepsilon}\phi(1-\phi)(1-2\phi) - 6\gamma\varepsilon\nabla^2\phi\right] - \alpha\phi^2(1-\phi)^2\left(\rho_v-\rho_{vs}(T)\right)=0
$$

*Uses:* S4, S2, S3. *Operation:* substitution of S2 and S3 into S4. *Checked by:* `check_units` in `scratch/check.py`: units phi_t 1/s; M*g/e 1/s; M*g*e/L^2 1/s; alpha*drho 1/s

## Result

$$
R = \partial_t\phi + M\left[\frac{6\gamma}{\varepsilon}\phi(1-\phi)(1-2\phi) - 6\gamma\varepsilon\left(\partial_{xx}\phi+\partial_{yy}\phi+\partial_{zz}\phi\right)\right] - \alpha\,\phi^2(1-\phi)^2\left(\rho_v-\rho_{vs}(T)\right) = 0
$$

Graded SymPy form:

```
Derivative(phi(x, y, z, t), t) + M*(6*gamma/eps*phi(x, y, z, t)*(1 - phi(x, y, z, t))*(1 - 2*phi(x, y, z, t)) - 6*gamma*eps*(Derivative(phi(x, y, z, t), (x, 2)) + Derivative(phi(x, y, z, t), (y, 2)) + Derivative(phi(x, y, z, t), (z, 2)))) - alpha*phi(x, y, z, t)**2*(1 - phi(x, y, z, t))**2*(rho_v(x, y, z, t) - rho_vs(T(x, y, z, t)))
```

## Neglected effects

| Effect | Why it is left out |
|---|---|
| Gibbs-Thomson shift of $\rho_{vs}$ with curvature | The relaxation term $-M\mu$ already produces curvature-driven motion. The explicit shift $\rho_{vs}\,2\gamma v_m/(\kappa_B T r)$ is about $10^{-2}$ for $r\sim 0.1\,\mu$m and smaller for larger pores. |
| Dependence of $M$ and $\alpha$ on $T$ and $\phi$ | Kept constant for simplicity. They can be promoted to $M(\phi,T)$ and $\alpha(\phi,T)$ without changing the structure of the equation. |
| Anisotropy of $\gamma$ and of the attachment kinetics (facets) | Ice basal and prismatic surface energies differ by only a few percent. The model is isotropic. |
| Ice-grain (regolith) wetting and contact angle | Grains are not part of this two-phase model; they would enter as boundary conditions or a third phase. |
| Density change, mechanical stress and flow of ice | Ice is rigid at about 200 K on the time scales of interest. Stress energy is small compared with $\gamma/\varepsilon$. |
| Thermal noise in $\phi$ | $\kappa_B T$ is much smaller than $\gamma\varepsilon^2$ for $\varepsilon$ above about 1 nm, so fluctuations are negligible. |
| Coupling of the source to the vapor and heat equations (mass sink, latent heat) | Not neglected physically. It belongs to the separate $\rho_v$ and $T$ tasks, since here $\rho_v$ and $T$ are given fields. |

## Questions

- The interface weight $g(\phi)=\phi^2(1-\phi)^2$ and the constant $M$, $\alpha$ are my modeling choices (as in Moure & Fu-type ice models). If the reference model uses another weight, such as $\phi(1-\phi)$ or $30\phi^2(1-\phi)^2$, or a source divided explicitly by $\rho_i$, the coefficient $\alpha$ changes accordingly.

