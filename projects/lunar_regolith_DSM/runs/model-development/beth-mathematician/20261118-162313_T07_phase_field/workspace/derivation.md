# T07_phase_field

Run `20261118-162313_T07_phase_field` · status **done** · written by `scripts/agentkit/render_derivation.py` from `answer.json`; do not edit by hand.

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
| $M(T)$ | Allen-Cahn mobility of the phase field. It sets how fast the interface relaxes and moves. Through $v_n = 6\varepsilon M\Delta g$ it is tied to the attachment kinetics; matching the Hertz-Knudsen law gives $M = \alpha_c\,\rho_{vs}\sqrt{R_vT/2\pi}/(6\varepsilon\rho_i^2R_vT)$ | $\mathrm{m^{3}\,J^{-1}\,s^{-1}}$ | 1e-09 | introduced |
| $\rho_i$ | density of ice, taken as constant | $\mathrm{kg\,m^{-3}}$ | 920 | introduced |
| $R_v$ | specific gas constant of water vapor | $\mathrm{J\,kg^{-1}\,K^{-1}}$ | 462 | introduced |
| $h(\phi)$ | interpolation function $\phi^3(10-15\phi+6\phi^2)$, going from 0 in the gas to 1 in the ice | 1 | 0.5 | introduced |
| $\Delta g$ | bulk driving force for deposition: the drop in free energy per unit volume of ice formed from vapor at density $\rho_v$ | $\mathrm{J\,m^{-3}}$ | 1e+03 | introduced |
| $G[\phi]$ | total free energy, including the bulk drive | $\mathrm{J}$ | 1e-12 | introduced |
| $v_n$ | normal velocity of a flat interface, positive toward the gas (growth) | $\mathrm{m\,s^{-1}}$ | 1e-09 | introduced |
| $\phi_0(X)$ | 1D equilibrium profile across a flat interface, with normal coordinate $X$ | 1 | 0.5 | introduced |
| $X$ | coordinate normal to a flat interface | $\mathrm{m}$ | 1e-06 | introduced |
| $\alpha_c$ | deposition (sticking) coefficient in the Hertz-Knudsen law | 1 | 0.5 | introduced |
| $R$ | residual of the phase-field equation | $\mathrm{s^{-1}}$ | 0 | introduced |

## Given

- `F = Integral(3*gamma/eps*phi**2*(1 - phi)**2 + 3*gamma*eps*(Derivative(phi, x)**2 + Derivative(phi, y)**2 + Derivative(phi, z)**2), over the volume)`

## Assumptions

- **A1.** A1: $\rho_v(x,y,z,t)$ and $T(x,y,z,t)$ are known, smooth fields (given). In the variation of $G$ they are held fixed.
- **A2.** A2: $\phi$ is a non-conserved order parameter. Ice mass is exchanged with the vapor field, which is a separate equation, so $\phi$ follows a linear gradient (Allen-Cahn) flow with a scalar mobility $M(T)$ that does not depend on $\phi$ or on the interface orientation.
- **A3.** A3: Water vapor is an ideal gas, and ice is incompressible with constant density $\rho_i$.
- **A4.** A4: The supersaturation is small, so $\ln(\rho_v/\rho_{vs}) \approx (\rho_v-\rho_{vs})/\rho_{vs}$.
- **A5.** A5: The drive enters through $h(\phi)$ with $h'(0)=h'(1)=0$. It acts only where $0<\phi<1$ (no nucleation), and it leaves the minima of the double well at 0 and 1.
- **A6.** A6: The flat interface described by $F$ has energy $\gamma$ per unit area and a 10-90 % width of $2\ln 9\,\varepsilon$ (given).

## Derivation

### S1. Variational derivative of F

Allen-Cahn relaxation is driven by $\delta F/\delta\phi$, so we compute it first from the given functional. The Euler-Lagrange operator gives the derivative of the double well minus 6$\gamma\varepsilon$ times the Laplacian.

$$
\frac{\delta F}{\delta\phi} = \frac{6\gamma}{\varepsilon}\phi(1-\phi)(1-2\phi) - 6\gamma\varepsilon\left(\partial_{xx}\phi+\partial_{yy}\phi+\partial_{zz}\phi\right)
$$

*Uses:* given. *Operation:* Euler-Lagrange variation of F, then integration by parts with no boundary flux. *Checked by:* `check_s1` in `scratch/check.py`: S1 dF/dphi residual: 0

### S2. Interpolation function for the bulk drive

The bulk drive must change the free-energy difference between ice and gas without shifting the minima at $\phi=0,1$, and its derivative must vanish in pure gas. That is the no-nucleation requirement. The standard choice $h$ has $h(0)=0$, $h(1)=1$ and $h'(0)=h'(1)=0$.

$$
h(\phi)=\phi^3(10-15\phi+6\phi^2),\qquad h'(\phi)=30\phi^2(1-\phi)^2
$$

*Uses:* given. *Operation:* definition and differentiation. *Checked by:* `check_s2` in `scratch/check.py`: S2 h' - 30 p^2(1-p)^2: 0  h(0), h(1): 0 1  h'(0),h'(1): 0 0

### S3. Bulk driving force from the vapor supersaturation

Moving mass $\rho_i$ per unit volume from an ideal vapor at $\rho_v$ into ice in equilibrium with $\rho_{vs}(T)$ releases $\rho_i R_v T\ln(\rho_v/\rho_{vs})$ per unit volume. We linearize the logarithm, dropping terms of order $(\rho_v/\rho_{vs}-1)^2$, because the allowed operations have no log and pore supersaturations are small.

$$
\Delta g = \rho_i R_v T\ln\frac{\rho_v}{\rho_{vs}(T)} \approx \rho_i R_v T\,\frac{\rho_v-\rho_{vs}(T)}{\rho_{vs}(T)}
$$

*Uses:* A3, A4. *Operation:* ideal-gas chemical potential difference, then a first-order Taylor expansion of ln about 1. *Checked by:* `check_units` in `scratch/check.py`: units M rho_i R_v T (rho_v-rho_vs)/rho_vs 1 / second (so Delta g has units of J m^-3)

### S4. Total free energy and Allen-Cahn relaxation

Ice ($\phi=1$) is favored by $\Delta g$ when the vapor is supersaturated, so the total energy is $G = F - \int\Delta g\,h(\phi)\,dV$. $\phi$ is not conserved because mass leaves through the vapor field, so it relaxes by the non-conserved gradient flow with mobility $M$. We hold $\Delta g$ fixed in the variation because $\rho_v$ and $T$ are given fields.

$$
\partial_t\phi = -M(T)\frac{\delta G}{\delta\phi} = -M(T)\left[\frac{6\gamma}{\varepsilon}\phi(1-\phi)(1-2\phi) - 6\gamma\varepsilon\nabla^2\phi - 30\phi^2(1-\phi)^2\Delta g\right]
$$

*Uses:* S1, S2, S3, A1, A2, A5. *Operation:* non-conserved gradient flow, with S1 and S2 substituted. *Checked by:* `check_units` in `scratch/check.py`: each term has units 1 / second

### S5. Equilibrium profile reproduces gamma and the width

This checks that the interface part matches the stated assumption. With $\Delta g=0$, the flat profile $\phi_0$ solves $\delta F/\delta\phi=0$, its energy per unit area is $\gamma$, and its 10-90 % width is $4\varepsilon\ln 3 = 2\varepsilon\ln 9$.

$$
\phi_0=\tfrac12\left[1+\tanh\tfrac{X}{2\varepsilon}\right],\quad \int_{-\infty}^{\infty}(\phi_0')^2dX=\frac{1}{6\varepsilon}
$$

*Uses:* S1, A6. *Operation:* substitute the tanh profile and integrate over X. *Checked by:* `check_s5` in `scratch/check.py`: S5 profile ODE residual: 0  int phi0^2 = 1/(6*eps)  energy = gamma  10-90 width = 4*eps*log(3)

### S6. Meaning of M: the interface velocity

Multiply the 1D traveling-wave form of S4 by $\phi_0'$ and integrate across the interface. This gives the normal velocity $v_n = 6\varepsilon M\Delta g$ for a flat interface. Matching it to the Hertz-Knudsen deposition velocity $v_n=(\alpha_c/\rho_i)\sqrt{R_vT/2\pi}\,(\rho_v-\rho_{vs})$ fixes $M(T)$, as given in the symbols list.

$$
v_n\int(\phi_0')^2dX = M\Delta g\int h'(\phi_0)\phi_0'\,dX \;\Rightarrow\; v_n = 6\varepsilon M\Delta g
$$

*Uses:* S4, S5. *Operation:* multiply by phi_0' and integrate over X (solvability condition). *Checked by:* `check_s6` in `scratch/check.py`: S6 int h'(phi0) phi0' dX = 1  => v = 6*Delta_g*M*eps

### S7. Residual form in 3D Cartesian coordinates

Substitute S3 into S4, write out the Laplacian in Cartesian coordinates, and move every term to one side. The result is $R$.

$$
R = \partial_t\phi + M(T)\left[\frac{6\gamma}{\varepsilon}\phi(1-\phi)(1-2\phi) - 6\gamma\varepsilon(\partial_{xx}\phi+\partial_{yy}\phi+\partial_{zz}\phi) - 30\phi^2(1-\phi)^2\rho_iR_vT\frac{\rho_v-\rho_{vs}(T)}{\rho_{vs}(T)}\right]=0
$$

*Uses:* S3, S4. *Operation:* substitution of S3 into S4 and rearrangement. *Checked by:* `check_units` in `scratch/check.py`: dphi/dt, M gamma/eps, M gamma eps/m^2 and M rho_i R_v T (rho_v-rho_vs)/rho_vs all have units 1 / second

## Result

$$
R = \partial_t\phi + M(T)\left[\frac{6\gamma}{\varepsilon}\phi(1-\phi)(1-2\phi) - 6\gamma\varepsilon\left(\partial_{xx}\phi+\partial_{yy}\phi+\partial_{zz}\phi\right) - 30\,\phi^2(1-\phi)^2\,\rho_i R_v T\,\frac{\rho_v-\rho_{vs}(T)}{\rho_{vs}(T)}\right] = 0
$$

Graded SymPy form:

```
Derivative(phi(x, y, z, t), t) + M(T(x, y, z, t))*(6*gamma/eps*phi(x, y, z, t)*(1 - phi(x, y, z, t))*(1 - 2*phi(x, y, z, t)) - 6*gamma*eps*(Derivative(phi(x, y, z, t), (x, 2)) + Derivative(phi(x, y, z, t), (y, 2)) + Derivative(phi(x, y, z, t), (z, 2))) - 30*phi(x, y, z, t)**2*(1 - phi(x, y, z, t))**2*rho_i*R_v*T(x, y, z, t)*(rho_v(x, y, z, t) - rho_vs(T(x, y, z, t)))/rho_vs(T(x, y, z, t)))
```

## Neglected effects

| Effect | Why it is left out |
|---|---|
| Nonlinear supersaturation terms in $\ln(\rho_v/\rho_{vs})$ | They are of order $(\rho_v/\rho_{vs}-1)^2$, which is below about 1 % for supersaturations up to about 10 %. |
| Anisotropy of $\gamma$ and of the attachment kinetics (facets of ice) | The work order gives a constant, isotropic $\gamma$. Ice anisotropy at 200 K is a few percent in $\gamma$; kinetic anisotropy is larger but is not resolved at pore scale. |
| Gibbs-Thomson shift of $\rho_{vs}$ written as an explicit term | It is not left out: the curvature term comes out of $-6\gamma\varepsilon\nabla^2\phi$ in the sharp-interface limit, so adding it explicitly would count it twice. |
| Latent heat release and the vapor sink at the interface | These belong to the $T$ and $\rho_v$ equations, which are separate tasks. |
| Ice-grain (regolith) wetting, contact angle, and a third solid phase | The work order treats the pore space as two phases only. |
| Thermal noise and nucleation | No nucleation is stated. The critical nucleus energy barrier is many $k_BT$ at small supersaturation. |
| Elastic stress in the ice and flow or advection of the gas | Stresses from sublimation are small compared with $\Delta g$, and the pore gas is near rest at regolith pressures. |
| Dependence of $M$ on $\phi$ | A constant $M$ across the interface keeps the sharp-interface relation $v_n=6\varepsilon M\Delta g$ exact for a flat interface. A $\phi$-dependent $M$ changes only the $O(\varepsilon)$ corrections. |

