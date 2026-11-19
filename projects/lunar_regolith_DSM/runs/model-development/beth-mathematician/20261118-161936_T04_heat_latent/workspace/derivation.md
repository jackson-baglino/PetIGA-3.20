# T04_heat_latent

Run `20261118-161936_T04_heat_latent` · status **done** · written by `scripts/agentkit/render_derivation.py` from `answer.json`; do not edit by hand.

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
| $e$ | internal energy (enthalpy) per unit volume of the ice-gas mixture | $\mathrm{J\,m^{-3}}$ | -5e+08 | introduced |
| $\mathbf{q}$ | conductive heat flux | $\mathrm{W\,m^{-2}}$ | 10 | introduced |
| $\mathbf{j}_v$ | diffusive mass flux of water vapor | $\mathrm{kg\,m^{-2}\,s^{-1}}$ | 1e-09 | introduced |
| $\dot m$ | sublimation rate: mass of ice turned into vapor per unit volume and time | $\mathrm{kg\,m^{-3}\,s^{-1}}$ | 1e-06 | introduced |
| $\rho_i$ | density of ice (constant) | $\mathrm{kg\,m^{-3}}$ | 917 | introduced |
| $\rho_v$ | mass density of water vapor in the gas | $\mathrm{kg\,m^{-3}}$ | 1e-07 | introduced |
| $h_i(T)$ | specific enthalpy of ice | $\mathrm{J\,kg^{-1}}$ | -2.8e+06 | introduced |
| $h_v(T)$ | specific enthalpy of water vapor | $\mathrm{J\,kg^{-1}}$ | 0 | introduced |
| $c_i(T)$ | specific heat of ice, dh_i/dT | $\mathrm{J\,kg^{-1}\,K^{-1}}$ | 1.5e+03 | introduced |
| $c_v(T)$ | specific heat of water vapor at constant pressure, dh_v/dT | $\mathrm{J\,kg^{-1}\,K^{-1}}$ | 1.85e+03 | introduced |
| $L_{sub}(T)$ | latent heat of sublimation, h_v(T) - h_i(T) | $\mathrm{J\,kg^{-1}}$ | 2.8e+06 | introduced |
| $C(\phi,T)$ | volumetric heat capacity of the ice-gas mixture | $\mathrm{J\,m^{-3}\,K^{-1}}$ | 1.4e+06 | introduced |
| $k(\phi,T)$ | thermal conductivity of the ice-gas mixture | $\mathrm{W\,m^{-1}\,K^{-1}}$ | 3.5 | introduced |
| $R$ | residual of the heat equation | $\mathrm{W\,m^{-3}}$ | 0 | introduced |

## Assumptions

- **A1.** A1: $\phi(x,y,z,t)$ is a known, smooth field (given), and $T$ is smooth.
- **A2.** A2: There are no external heat sources (no radiation absorbed in the pore volume, no Joule heating).
- **A3.** A3: Ice and gas share one local temperature $T$ in the diffuse interface (local thermal equilibrium).
- **A4.** A4: Pressure work and kinetic energy are negligible, so internal energy and enthalpy are interchangeable.
- **A5.** A5: The ice is rigid and at rest, with constant density $\rho_i$, so $\dot m = -\rho_i\partial_t\phi$.
- **A6.** A6: The gas is water vapor only (any inert gas is negligible in the lunar pore space); it moves only by diffusion with flux $\mathbf{j}_v$, and there is no bulk gas flow.
- **A7.** A7: Heat conduction follows Fourier's law $\mathbf{q}=-k(\phi,T)\nabla T$ with an isotropic conductivity.
- **A8.** A8: The heat carried by diffusing vapor, $c_v\mathbf{j}_v\cdot\nabla T$, is negligible next to conduction.

## Derivation

### S1. Energy conservation, integral form

Start from conservation of energy for a fixed control volume $V$ in the pore space. Energy changes only through the flux across its boundary, because there are no external heat sources.

$$
\frac{d}{dt}\int_V e\,dV = -\oint_{\partial V}\left(\mathbf{q} + h_v\,\mathbf{j}_v\right)\cdot\mathbf{n}\,dS
$$

*Uses:* given, A2. *Operation:* statement of the conservation principle. *Checked by:* none (This is a postulate (first law for an open control volume, with pressure work neglected per A4), so no tool can check it.)

### S2. Local form

Apply the divergence theorem. Because $V$ is arbitrary and the fields are smooth (A1), the integrand must vanish pointwise.

$$
\partial_t e + \nabla\cdot\left(\mathbf{q} + h_v\mathbf{j}_v\right) = 0
$$

*Uses:* S1, A1. *Operation:* divergence theorem; arbitrary V. *Checked by:* none (Standard localization argument; no symbolic content to check.)

### S3. Energy density and mass balances

Write the energy per volume as ice (fraction $\phi$) plus vapor in the gas (fraction $1-\phi$), each at the local $T$. Ice mass is lost at rate $\dot m$, which becomes vapor. Ice does not move.

$$
e = \phi\rho_i h_i(T) + (1-\phi)\rho_v h_v(T),\quad \partial_t(\phi\rho_i) = -\dot m,\quad \partial_t\left((1-\phi)\rho_v\right) + \nabla\cdot\mathbf{j}_v = \dot m
$$

*Uses:* A3, A5, A6. *Operation:* definition of e; mass balances for ice and vapor. *Checked by:* none (Definitions and modeling assumptions; their use is checked in check_s4.)

### S4. Expand and use the mass balances

Apply the product rule to $\partial_t e$ and $\nabla\cdot(h_v\mathbf{j}_v)$, then substitute the two mass balances of S3. The terms with $h_v\nabla\cdot\mathbf{j}_v$ cancel, and the enthalpy difference multiplies $\dot m$.

$$
\left[\phi\rho_i c_i + (1-\phi)\rho_v c_v\right]\partial_t T + \dot m\left(h_v - h_i\right) + c_v\,\mathbf{j}_v\cdot\nabla T + \nabla\cdot\mathbf{q} = 0
$$

*Uses:* S2, S3. *Operation:* product rule, substitution of S3. *Checked by:* `check_s4` in `scratch/check.py`: product rule ok: True; mass-balance substitution ok: True

### S5. Constitutive closures

Define the latent heat and the mixture heat capacity. Use Fourier's law for conduction with a conductivity that depends on phase and temperature. Express $\dot m$ through $\phi$, which needs constant $\rho_i$ (A5).

$$
L_{sub}(T) = h_v - h_i,\quad C(\phi,T) = \phi\rho_i c_i + (1-\phi)\rho_v c_v,\quad \mathbf{q} = -k(\phi,T)\nabla T,\quad \dot m = -\rho_i\,\partial_t\phi
$$

*Uses:* S3, A5, A7. *Operation:* definitions; Fourier's law. *Checked by:* none (Definitions; $\dot m$ follows from S3 with $\partial_t\rho_i=0$.)

### S6. Drop vapor enthalpy advection

The term $c_v\mathbf{j}_v\cdot\nabla T$ is the heat carried by diffusing vapor. It is many orders of magnitude smaller than conduction (see neglected), so drop it.

$$
C(\phi,T)\,\partial_t T + \dot m\,L_{sub}(T) + \nabla\cdot\mathbf{q} = 0
$$

*Uses:* S4, A8. *Operation:* neglect term of relative order 1e-10. *Checked by:* none (Order-of-magnitude argument, given in neglected.)

### S7. Heat equation with latent heat, residual form

Substitute the closures of S5 and write the divergence in 3D Cartesian coordinates. Deposition ($\partial_t\phi>0$) releases heat and raises $T$; sublimation absorbs it. Every term is in W m$^{-3}$.

$$
R = C(\phi,T)\partial_t T - \partial_x(k\,\partial_x T) - \partial_y(k\,\partial_y T) - \partial_z(k\,\partial_z T) - \rho_i L_{sub}(T)\,\partial_t\phi = 0
$$

*Uses:* S5, S6. *Operation:* substitution of S5; Cartesian divergence. *Checked by:* `check_s7` in `scratch/check.py`: C T_t 1400000.0 watt / meter ** 3; div k grad T 3.5 watt / meter ** 3; rho L phi_t 2567600000.0 watt / meter ** 3 (all terms in W/m^3)

## Result

$$
R = C(\phi,T)\,\partial_t T - \partial_x\!\left(k(\phi,T)\,\partial_x T\right) - \partial_y\!\left(k(\phi,T)\,\partial_y T\right) - \partial_z\!\left(k(\phi,T)\,\partial_z T\right) - \rho_i\,L_{sub}(T)\,\partial_t\phi = 0
$$

Graded SymPy form:

```
C(phi(x, y, z, t), T(x, y, z, t))*Derivative(T(x, y, z, t), t) - Derivative(k(phi(x, y, z, t), T(x, y, z, t))*Derivative(T(x, y, z, t), x), x) - Derivative(k(phi(x, y, z, t), T(x, y, z, t))*Derivative(T(x, y, z, t), y), y) - Derivative(k(phi(x, y, z, t), T(x, y, z, t))*Derivative(T(x, y, z, t), z), z) - rho_i*L_sub(T(x, y, z, t))*Derivative(phi(x, y, z, t), t)
```

## Neglected effects

| Effect | Why it is left out |
|---|---|
| Heat carried by diffusing vapor, $c_v\mathbf{j}_v\cdot\nabla T$ | $\|\mathbf{j}_v\|\sim D\rho_v/\ell \sim 10^{-5}\cdot10^{-7}/10^{-5}\approx10^{-7}$ kg m$^{-2}$ s$^{-1}$, so $c_v\|\mathbf{j}_v\|\sim2\times10^{-4}$ W m$^{-2}$ K$^{-1}$, against $k/\ell\sim3.5/10^{-5}\approx3\times10^5$ W m$^{-2}$ K$^{-1}$: a ratio near $10^{-9}$. |
| Sensible heat of the vapor inside $C$ | It is kept formally, but $(1-\phi)\rho_v c_v\sim10^{-4}$ J m$^{-3}$ K$^{-1}$ against $\rho_i c_i\sim1.4\times10^6$, so $C\approx\phi\rho_i c_i$ in practice. |
| Inert gas (air) and its enthalpy transport | Lunar pore gas is near vacuum, so its heat capacity and its enthalpy flux are negligible. A model with air would add $(1-\phi)\rho_a c_a$ to $C$ and an air-flux term like the one above. |
| Pressure (compression) work and kinetic energy | Gas pressures are tiny (vapor pressure around 1e-3 Pa at 200 K) and nothing flows fast, so $p\,dV$ work is many orders of magnitude below the latent heat. |
| Radiative heat transfer across pores | $4\epsilon\sigma T^3\ell\sim4\cdot5.7\times10^{-8}\cdot8\times10^6\cdot10^{-4}\approx2\times10^{-4}$ W m$^{-1}$ K$^{-1}$ effective conductivity at 200 K and 100 $\mu$m pores, small next to ice $k\approx3.5$. It is excluded under A2. |
| Interface (surface) energy and curvature effect on latent heat | Capillary energy $\gamma/\ell\sim0.1/10^{-5}=10^4$ J m$^{-3}$ is much less than $\rho_i L_{sub}\approx2.6\times10^9$ J m$^{-3}$. |
| Density change of ice with temperature (thermal expansion) | The volumetric expansion coefficient of ice is about $1.5\times10^{-4}$ K$^{-1}$, so $\rho_i$ changes by under 1% over tens of kelvin. Constant $\rho_i$ is needed in S5 to write $\dot m=-\rho_i\partial_t\phi$. |
| Regolith grains (solid silicate) in the energy balance | The task models only the pore space (ice + gas). The grains enter through boundary conditions, not through this equation. |

