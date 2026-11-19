# T04_heat_latent

Run `20261118-172700_T04_heat_latent` · status **done** · written by `scripts/agentkit/render_derivation.py` from `answer.json`; do not edit by hand.

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
| $\rho_i$ | density of ice (constant, A3) | $\mathrm{kg\,m^{-3}}$ | 920 | introduced |
| $\rho_a$ | partial density of dry air in the gas | $\mathrm{kg\,m^{-3}}$ | 0.001 | introduced |
| $\rho_v$ | partial density of water vapor in the gas | $\mathrm{kg\,m^{-3}}$ | 1e-06 | introduced |
| $h_i(T)$ | specific enthalpy of ice | $\mathrm{J\,kg^{-1}}$ | -2.8e+05 | introduced |
| $h_a(T)$ | specific enthalpy of air | $\mathrm{J\,kg^{-1}}$ | 2e+05 | introduced |
| $h_v(T)$ | specific enthalpy of water vapor | $\mathrm{J\,kg^{-1}}$ | 2.55e+06 | introduced |
| $c_i(T)$ | specific heat of ice, dh_i/dT | $\mathrm{J\,kg^{-1}\,K^{-1}}$ | 1.4e+03 | introduced |
| $c_a(T)$ | specific heat of air at constant pressure, dh_a/dT | $\mathrm{J\,kg^{-1}\,K^{-1}}$ | 1e+03 | introduced |
| $c_v(T)$ | specific heat of water vapor at constant pressure, dh_v/dT | $\mathrm{J\,kg^{-1}\,K^{-1}}$ | 1.85e+03 | introduced |
| $L(T)$ | latent heat of sublimation, h_v - h_i | $\mathrm{J\,kg^{-1}}$ | 2.83e+06 | introduced |
| $C(\phi,T)$ | effective volumetric heat capacity, \phi\rho_i c_i + (1-\phi)(\rho_a c_a + \rho_v c_v) | $\mathrm{J\,m^{-3}\,K^{-1}}$ | 1.3e+06 | introduced |
| $k(\phi,T)$ | effective thermal conductivity of the ice-gas mixture | $\mathrm{W\,m^{-1}\,K^{-1}}$ | 2 | introduced |
| $H$ | enthalpy per unit volume of the mixture | $\mathrm{J\,m^{-3}}$ | -2.6e+08 | introduced |
| $\mathbf{q}$ | conductive heat flux | $\mathrm{W\,m^{-2}}$ | 1 | introduced |
| $\mathbf{j}_v$ | mass flux of water vapor | $\mathrm{kg\,m^{-2}\,s^{-1}}$ | 1e-07 | introduced |
| $\mathbf{j}_a$ | mass flux of air | $\mathrm{kg\,m^{-2}\,s^{-1}}$ | 1e-08 | introduced |

## Assumptions

- **A1.** A1: $\phi(x,y,z,t)$ is a known, smooth field (given). A point is a mixture: ice with volume fraction $\phi$ and gas with volume fraction $1-\phi$, and both have the same local temperature $T$ (local thermal equilibrium across the diffuse interface).
- **A2.** A2: The pressure is nearly uniform and steady, so enthalpy is the conserved energy density. Pressure work $\partial_t p$, kinetic energy and viscous dissipation are neglected.
- **A3.** A3: The ice density $\rho_i$ is constant (ice is incompressible, and its thermal expansion is neglected). Because of this, the ice mass per unit volume is $\rho_i\phi$ and the latent term is $\rho_i L\,\partial_t\phi$. All other properties may depend on $\phi$ and $T$; the specific enthalpies depend only on $T$.
- **A4.** A4: Heat conduction obeys Fourier's law with an effective conductivity $k(\phi,T)$. The gas is a mixture of air and water vapor, and the energy each species carries is its enthalpy flux.
- **A5.** A5: Sensible heat that moving vapor and air carry, $(c_v\mathbf{j}_v + c_a\mathbf{j}_a)\cdot\nabla T$, is negligible compared with conduction (Péclet number about $10^{-6}$).

## Derivation

### S1. Energy balance for a control volume

Energy conservation for a fixed control volume at constant pressure (A2) says the enthalpy stored changes only by conduction and by the enthalpy that moving vapor and air carry across the boundary. The divergence theorem turns this into a local statement.

$$
\partial_t H + \nabla\cdot\left(\mathbf{q} + h_v\mathbf{j}_v + h_a\mathbf{j}_a\right) = 0,\quad H = \phi\rho_i h_i(T) + (1-\phi)\left(\rho_a h_a(T) + \rho_v h_v(T)\right),\quad \mathbf{q} = -k(\phi,T)\nabla T
$$

*Uses:* A1, A2, A4. *Operation:* conservation of enthalpy in a fixed volume, divergence theorem. *Checked by:* `check_units` in `scratch/check.py`: every term converts to W m^-3

### S2. Water mass balance

Mass conservation for the water substance: the ice lost per unit volume per unit time, $-\rho_i\partial_t\phi$, becomes vapor, which is stored in the gas or carried away by $\mathbf{j}_v$. This identity is what makes the latent heat appear.

$$
\partial_t\left[(1-\phi)\rho_v\right] + \nabla\cdot\mathbf{j}_v = -\rho_i\,\partial_t\phi
$$

*Uses:* A3, A4. *Operation:* mass conservation of water (ice + vapor). *Checked by:* none (used as a substitution inside check_s4)

### S3. Air mass balance

Air is neither made nor destroyed. When ice grows, it pushes air out of the volume, and that air carries its enthalpy with it.

$$
\partial_t\left[(1-\phi)\rho_a\right] + \nabla\cdot\mathbf{j}_a = 0
$$

*Uses:* A4. *Operation:* mass conservation of air. *Checked by:* none (used as a substitution inside check_s4)

### S4. Expand and eliminate the mass balances

Expand $\partial_t H$ and $\nabla\cdot(h\mathbf{j})$ with the chain rule, using $dh/dT = c$. Then replace $\nabla\cdot\mathbf{j}_v$ and $\nabla\cdot\mathbf{j}_a$ with S2 and S3. The enthalpy terms that multiply $\partial_t\phi$ combine to $\rho_i(h_i - h_v) = -\rho_i L$.

$$
C(\phi,T)\,\partial_t T + \left(c_v\mathbf{j}_v + c_a\mathbf{j}_a\right)\cdot\nabla T - \nabla\cdot\left(k(\phi,T)\nabla T\right) - \rho_i L(T)\,\partial_t\phi = 0,\quad C = \phi\rho_i c_i + (1-\phi)(\rho_a c_a + \rho_v c_v),\ L = h_v - h_i
$$

*Uses:* S1, S2, S3. *Operation:* product and chain rule, substitution of S2 and S3. *Checked by:* `check_s4` in `scratch/check.py`: S4 energy - expected = 0

### S5. Drop advection of sensible heat by the gas

The term $(c_v\mathbf{j}_v + c_a\mathbf{j}_a)\cdot\nabla T$ compared with conduction is a Péclet number $|\mathbf{j}|c\ell/k \sim 10^{-7}\cdot 1850\cdot 10^{-4}/0.02 \sim 10^{-6}$ for pore size $\ell \sim 100\,\mu$m, so it is dropped (A5).

$$
C(\phi,T)\,\partial_t T - \nabla\cdot\left(k(\phi,T)\nabla T\right) - \rho_i L(T)\,\partial_t\phi = 0
$$

*Uses:* S4, A5. *Operation:* neglect a term of relative order $10^{-6}$. *Checked by:* `check_units` in `scratch/check.py`: j c grad T has units of W m^-3, like the other terms; the size estimate is by hand

### S6. Write in 3D Cartesian coordinates

Write the divergence of the flux in components to get the residual R that is returned. Sublimation ($\partial_t\phi<0$) acts as a heat sink and deposition as a heat source.

$$
R = C\,\partial_t T - \partial_x(k\,\partial_x T) - \partial_y(k\,\partial_y T) - \partial_z(k\,\partial_z T) - \rho_i L(T)\,\partial_t\phi = 0
$$

*Uses:* S5. *Operation:* expand the divergence in Cartesian components. *Checked by:* none (a direct rewriting of the divergence in components; no tool check)

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
| Advection of sensible heat by vapor and air fluxes | Péclet number $\|\mathbf{j}\|c\ell/k \approx 10^{-7}\cdot1850\cdot10^{-4}/0.02 \approx 10^{-6}$. |
| Pressure work and changes in total pressure | Total pressure in the pore is nearly constant on the time scales of sublimation. $\partial_t p$ is many orders of magnitude smaller than $C\,\partial_t T$. |
| Kinetic energy and viscous dissipation of the gas | Gas speeds are tiny (mm/s or less), so $\rho u^2/2$ is negligible next to the enthalpies (about $10^6$ J/kg). |
| Thermal expansion of ice (temperature-dependent $\rho_i$) | The expansion coefficient is about $10^{-4}$ K$^{-1}$, so the extra terms are of relative order $10^{-4}\,\Delta T$. |
| Interfacial (surface) energy and curvature (Gibbs-Thomson) correction to the latent heat | $\gamma/(\rho_i L\,r) \approx 0.1/(920\cdot2.8\times10^6\cdot10^{-6}) \approx 4\times10^{-5}$ for curvature radius $r \approx 1\,\mu$m. The phase-field equation for $\phi$ is a later task. |
| Thermal radiation across pores | Effective radiative conductivity $4\sigma T^3\ell \approx 4\cdot5.7\times10^{-8}\cdot8\times10^6\cdot10^{-4} \approx 2\times10^{-4}$ W/(m K), small next to the conductivity of ice (about 2-4 W/(m K)). It matters only in nearly empty, large pores. |
| Dufour effect and Soret effect | These cross effects are second order in dilute vapor diffusion. |
| Heat exchange with the regolith grains | The task restricts the domain to the pore space. The grains enter only through boundary conditions. |

