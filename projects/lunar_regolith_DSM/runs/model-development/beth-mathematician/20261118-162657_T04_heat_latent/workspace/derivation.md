# T04_heat_latent

Run `20261118-162657_T04_heat_latent` · status **done** · written by `scripts/agentkit/render_derivation.py` from `answer.json`; do not edit by hand.

## Work order (verbatim)

> We are modeling the pore space of a granular material (lunar regolith) in which water ice sublimates and deposits. Treat the pore space as a three-dimensional region filled with solid ice and gas (air that contains water vapor). The local volume fraction of ice is the phase field phi(x, y, z, t), between 0 and 1: phi = 1 in ice and phi = 0 in gas, with a thin, smooth transition between the two. Ice turns into vapor (sublimation) and vapor into ice (deposition), so phi changes in time. In this task phi is known; its own equation is a later task. Derive the partial differential equation for the temperature T(x, y, z, t) that expresses conservation of energy, including the heat absorbed or released by the phase change. Start from a conservation principle. State every assumption, and list every effect you leave out with the reason it is small. Keep the material properties general: do not assume a property is constant unless a step needs it, and say so when one does. Write the result as R = 0 in 3D Cartesian coordinates and return R.

## Symbols

| Symbol | Meaning | Units | Typical value (SI) | Source |
|---|---|---|---|---|
| $x$ | first position coordinate | $\mathrm{m}$ |  | work order |
| $y$ | second position coordinate | $\mathrm{m}$ |  | work order |
| $z$ | third position coordinate | $\mathrm{m}$ |  | work order |
| $t$ | time | $\mathrm{s}$ |  | work order |
| $T$ | temperature | $\mathrm{K}$ |  | work order |
| $\phi$ | local volume fraction of ice (phase field), between 0 and 1 | 1 |  | work order |
| $\rho_i$ | density of ice, taken constant (A3) | $\mathrm{kg\,m^{-3}}$ | 920 | introduced |
| $h_i(T)$ | specific enthalpy of ice | $\mathrm{J\,kg^{-1}}$ | 3e+05 | introduced |
| $h_v(T)$ | specific enthalpy of water vapor | $\mathrm{J\,kg^{-1}}$ | 3.1e+06 | introduced |
| $h_a(T)$ | specific enthalpy of the non-condensable gas (air) | $\mathrm{J\,kg^{-1}}$ | 2e+05 | introduced |
| $\rho_v$ | mass of water vapor per unit gas volume | $\mathrm{kg\,m^{-3}}$ | 1e-07 | introduced |
| $m_a$ | mass of air per unit total volume, $(1-\phi)\rho_a$ | $\mathrm{kg\,m^{-3}}$ | 0.01 | introduced |
| $\mathbf{j}_v$ | diffusive mass flux of water vapor | $\mathrm{kg\,m^{-2}\,s^{-1}}$ | 1e-09 | introduced |
| $E$ | enthalpy per unit total volume | $\mathrm{J\,m^{-3}}$ | 3e+08 | introduced |
| $\mathbf{q}$ | total energy flux (conduction plus enthalpy carried by vapor) | $\mathrm{W\,m^{-2}}$ | 1e+03 | introduced |
| $k(\phi,T)$ | effective thermal conductivity of the ice/gas mixture | $\mathrm{W\,m^{-1}\,K^{-1}}$ | 2 | introduced |
| $c_i(T)$ | specific heat of ice, $dh_i/dT$ | $\mathrm{J\,kg^{-1}\,K^{-1}}$ | 1.5e+03 | introduced |
| $c_v(T)$ | specific heat (constant pressure) of water vapor, $dh_v/dT$ | $\mathrm{J\,kg^{-1}\,K^{-1}}$ | 1.85e+03 | introduced |
| $c_a(T)$ | specific heat of air, $dh_a/dT$ | $\mathrm{J\,kg^{-1}\,K^{-1}}$ | 1e+03 | introduced |
| $C(\phi,T)$ | effective volumetric heat capacity, $\rho_i\phi c_i + m_a c_a + (1-\phi)\rho_v c_v$ | $\mathrm{J\,m^{-3}\,K^{-1}}$ | 1.4e+06 | introduced |
| $L(T)$ | specific latent heat of sublimation, $h_v(T)-h_i(T)$ | $\mathrm{J\,kg^{-1}}$ | 2.83e+06 | introduced |
| $R$ | residual of the energy equation | $\mathrm{W\,m^{-3}}$ | 0 | introduced |

## Assumptions

- **A1.** A1: $\phi$ is a known, smooth field (given).
- **A2.** A2: Pressure work, viscous heating, and kinetic energy are negligible, so the stored energy is the enthalpy.
- **A3.** A3: The ice density $\rho_i$ is constant (rigid ice, no thermal expansion). This step needs it so that $\rho_i\phi$ is the ice mass per unit volume and $\partial_t(\rho_i\phi)=\rho_i\partial_t\phi$.
- **A4.** A4: The ice is fixed and the gas has no bulk flow. Only water vapor moves, by diffusion with flux $\mathbf{j}_v$. The air mass per unit volume $m_a$ is constant in time.
- **A5.** A5: Heat is conducted by Fourier's law with an effective conductivity $k(\phi,T)$. Thermal radiation across pores is folded into $k$ or neglected.
- **A6.** A6: Ice and gas share one local temperature $T$ (local thermal equilibrium), and the vapor flux is small (see S5).

## Derivation

### S1. Enthalpy per unit volume

Energy is stored in the ice, the air and the vapor. With pressure work and kinetic energy left out, the stored energy is enthalpy, so the energy per unit total volume is the sum over the phases.

$$
E = \rho_i\,\phi\,h_i(T) + m_a\,h_a(T) + (1-\phi)\,\rho_v\,h_v(T)
$$

*Uses:* A1, A2, A3. *Operation:* sum of phase enthalpies weighted by mass per unit volume. *Checked by:* none (Definition. It cannot be checked by a tool.)

### S2. Energy conservation

For a fixed control volume, $E$ changes only by the energy flux through its surface. The flux is Fourier conduction plus the enthalpy carried by diffusing vapor. The ice is fixed in place and the air is at rest.

$$
\partial_t E + \nabla\cdot\mathbf{q} = 0,\qquad \mathbf{q} = -k(\phi,T)\nabla T + h_v(T)\,\mathbf{j}_v
$$

*Uses:* S1, A4, A5. *Operation:* conservation law in local form (divergence theorem on an arbitrary control volume). *Checked by:* none (Statement of the conservation principle. It cannot be checked by a tool.)

### S3. Water mass conservation

Total water (ice plus vapor) is conserved, and only vapor moves. This balance links the vapor flux divergence to the change in ice fraction, which brings in the phase change.

$$
\partial_t\left[\rho_i\phi + (1-\phi)\rho_v\right] + \nabla\cdot\mathbf{j}_v = 0
$$

*Uses:* A3, A4. *Operation:* conservation law in local form. *Checked by:* none (Statement of the conservation principle. It cannot be checked by a tool.)

### S4. Expand and use the mass balance

Apply the product and chain rules to $\partial_t E$ and $\nabla\cdot(h_v\mathbf{j}_v)$, with $\partial_t m_a = 0$. The terms multiplying $h_v$ are $h_v$ times the left side of S3 minus $\rho_i\partial_t\phi$, so they combine with $\rho_i h_i\partial_t\phi$ into $-\rho_i(h_v-h_i)\partial_t\phi$.

$$
C(\phi,T)\,\partial_t T - \nabla\cdot\left(k\nabla T\right) - \rho_i L(T)\,\partial_t\phi + c_v\,\mathbf{j}_v\cdot\nabla T = 0,\quad C = \rho_i\phi c_i + m_a c_a + (1-\phi)\rho_v c_v,\ L = h_v - h_i
$$

*Uses:* S2, S3, A4. *Operation:* product rule, chain rule $\partial_t h(T)=c\,\partial_t T$, substitution of S3. *Checked by:* `check_s4` in `scratch/check.py`: S4 residual (should be 0): 0

### S5. Drop vapor enthalpy advection

The term $c_v\mathbf{j}_v\cdot\nabla T$ is heat carried by the vapor flux. Compared with conduction it scales as $|\mathbf{j}_v| c_v \ell/k \approx 2\times10^{-7}\,\mathrm{m^{-1}}\cdot\ell \sim 10^{-11}$ for a pore scale $\ell\sim10^{-4}$ m, so it is dropped.

$$
C(\phi,T)\,\partial_t T = \nabla\cdot\left(k(\phi,T)\nabla T\right) + \rho_i L(T)\,\partial_t\phi
$$

*Uses:* S4, A6. *Operation:* neglect a term of relative order $|\mathbf{j}_v|c_v\ell/k$. *Checked by:* `check_size` in `scratch/check.py`: j_v c_v |gradT| = 4.2e-4 W/m^3; ratio j_v c_v/k = 2.1e-7 m^-1

### S6. Residual in Cartesian coordinates

Write the divergence in $x,y,z$ and move every term to one side. Every term has units of W m$^{-3}$. Sublimation ($\partial_t\phi<0$) cools the material, as it should.

$$
R = C\,\partial_t T - \partial_x(k\,\partial_x T) - \partial_y(k\,\partial_y T) - \partial_z(k\,\partial_z T) - \rho_i L(T)\,\partial_t\phi = 0
$$

*Uses:* S5. *Operation:* expand $\nabla\cdot$ in Cartesian coordinates and collect terms in $R$. *Checked by:* `check_units` in `scratch/check.py`: 2.8e8 W/m^3, 400 W/m^3, 2.6e9 W/m^3: all three terms in W/m^3

## Result

$$
R = C(\phi,T)\,\partial_t T - \partial_x\!\left(k(\phi,T)\,\partial_x T\right) - \partial_y\!\left(k(\phi,T)\,\partial_y T\right) - \partial_z\!\left(k(\phi,T)\,\partial_z T\right) - \rho_i\,L(T)\,\partial_t\phi = 0
$$

Graded SymPy form:

```
C(phi(x, y, z, t), T(x, y, z, t))*Derivative(T(x, y, z, t), t) - Derivative(k(phi(x, y, z, t), T(x, y, z, t))*Derivative(T(x, y, z, t), x), x) - Derivative(k(phi(x, y, z, t), T(x, y, z, t))*Derivative(T(x, y, z, t), y), y) - Derivative(k(phi(x, y, z, t), T(x, y, z, t))*Derivative(T(x, y, z, t), z), z) - rho_i*L(T(x, y, z, t))*Derivative(phi(x, y, z, t), t)
```

## Neglected effects

| Effect | Why it is left out |
|---|---|
| Enthalpy carried by vapor diffusion, $c_v\mathbf{j}_v\cdot\nabla T$ | Its ratio to conduction is about $10^{-11}$ at pore scale for $\rho_v\sim10^{-7}$ kg m$^{-3}$ near 200 K (check_size). |
| Pressure work and kinetic energy of the gas | The gas is at rest and at very low pressure. $Dp/Dt$ work is of order $\rho_v$ times the vapor heat content, which is far below $\rho_i c_i \partial_t T$. |
| Advection by bulk gas flow, and air mass change as pores open (if $m_a$ is not constant) | $m_a\ll\rho_i$ (lunar pore gas is near vacuum), so any $h_a\partial_t m_a$ term is smaller than the latent term by $m_a h_a/(\rho_i L)\lesssim10^{-6}$. |
| Thermal expansion of ice (a variable $\rho_i$) | The relative density change is about $10^{-4}$ per K times $\Delta T$, which is small over the temperature range. |
| Interface (surface) energy and curvature (Gibbs-Thomson) contribution to the energy | The surface energy per unit volume, $\gamma/\ell\sim0.1/10^{-4}=10^3$ J m$^{-3}$, is far below the latent heat $\rho_i L\sim2.6\times10^9$ J m$^{-3}$. |
| Radiative heat transfer across pores | At about 200 K, radiative conductance $4\epsilon\sigma T^3\ell\sim2\times10^{-4}$ W m$^{-1}$ K$^{-1}$ for $\ell=10^{-4}$ m. This is far below $k$ of ice ($\sim2$). It can be included in $k$ if needed. |

