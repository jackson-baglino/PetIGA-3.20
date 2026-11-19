# T04_heat_latent

Run `20261118-173148_T04_heat_latent_baseline` · status **done** · written by `scripts/agentkit/render_derivation.py` from `answer.json`; do not edit by hand.

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
| $\rho_i$ | density of ice (taken constant, see A4) | $\mathrm{kg\,m^{-3}}$ | 920 | introduced |
| $\rho_v$ | mass concentration of water vapor in the gas | $\mathrm{kg\,m^{-3}}$ | 2e-06 | introduced |
| $h_i(T)$ | specific enthalpy of ice | $\mathrm{J\,kg^{-1}}$ | -3e+05 | introduced |
| $h_v(T)$ | specific enthalpy of water vapor | $\mathrm{J\,kg^{-1}}$ | 2.5e+06 | introduced |
| $c_i(T)$ | specific heat of ice, dh_i/dT | $\mathrm{J\,kg^{-1}\,K^{-1}}$ | 1.4e+03 | introduced |
| $c_v(T)$ | specific heat of water vapor at constant pressure, dh_v/dT | $\mathrm{J\,kg^{-1}\,K^{-1}}$ | 1.85e+03 | introduced |
| $H$ | enthalpy per unit volume of the ice-gas mixture | $\mathrm{J\,m^{-3}}$ | 1e+06 | introduced |
| $\mathbf q$ | conductive heat flux | $\mathrm{W\,m^{-2}}$ | 10 | introduced |
| $\mathbf J_v$ | diffusive mass flux of water vapor | $\mathrm{kg\,m^{-2}\,s^{-1}}$ | 1e-08 | introduced |
| $\dot m$ | mass of vapor deposited as ice per unit volume and time (negative for sublimation) | $\mathrm{kg\,m^{-3}\,s^{-1}}$ | 1e-06 | introduced |
| $k(\phi,T)$ | effective thermal conductivity of the ice-gas mixture | $\mathrm{W\,m^{-1}\,K^{-1}}$ | 2.5 | introduced |
| $C(\phi,T)$ | volumetric heat capacity of the mixture, phi rho_i c_i + (1-phi) rho_v c_v (+ air if present) | $\mathrm{J\,m^{-3}\,K^{-1}}$ | 1.3e+06 | introduced |
| $L(T)$ | specific latent heat of sublimation, h_v - h_i | $\mathrm{J\,kg^{-1}}$ | 2.84e+06 | introduced |
| $R$ | residual of the energy equation | $\mathrm{W\,m^{-3}}$ | 1e+03 | introduced |

## Assumptions

- **A1.** A1: $\phi(x,y,z,t)$ is a known, smooth field (given), and $T$ is smooth.
- **A2.** A2: The pressure is uniform and constant, and the system is at rest, so the energy content is the enthalpy. Pressure work, viscous dissipation and kinetic energy are dropped.
- **A3.** A3: The pore space holds only ice (fraction $\phi$) and gas (fraction $1-\phi$). Ice and gas share one local temperature $T$, and the ice is immobile.
- **A4.** A4: The ice density $\rho_i$ is constant (needed in S3 to write $\dot m=\rho_i\partial_t\phi$).
- **A5.** A5: Heat conduction follows Fourier's law with an isotropic effective conductivity $k(\phi,T)$.
- **A6.** A6: Enthalpy advection by vapor diffusion, $c_v\mathbf J_v\cdot\nabla T$, is negligible compared with conduction.

## Derivation

### S1. Integral energy balance

Conservation of energy for any fixed control volume $V$ in the pore space. Pressure is taken as uniform, so the stored energy is the enthalpy $H$. Energy enters through the conductive flux $\mathbf q$ and the enthalpy carried by diffusing vapor, $h_v\mathbf J_v$.

$$
\frac{d}{dt}\int_V H\,dV = -\oint_{\partial V}\left(\mathbf q + h_v\mathbf J_v\right)\cdot\mathbf n\,dA
$$

*Uses:* A1, A2. *Operation:* first law for a fixed control volume at constant pressure. *Checked by:* none (no tools in this run)

### S2. Local energy balance

$V$ is arbitrary and the fields are smooth (A1), so applying the divergence theorem gives a pointwise equation.

$$
\partial_t H + \nabla\cdot\left(\mathbf q + h_v\mathbf J_v\right) = 0
$$

*Uses:* S1, A1. *Operation:* divergence theorem, arbitrary V. *Checked by:* none (no tools in this run)

### S3. Mass balances of ice and vapor

Ice is immobile and gains mass $\dot m$ by deposition. The vapor loses the same mass. With constant ice density (A4), the ice balance links $\dot m$ to $\partial_t\phi$.

$$
\rho_i\partial_t\phi = \dot m,\qquad \partial_t\left[(1-\phi)\rho_v\right] + \nabla\cdot\mathbf J_v = -\rho_i\partial_t\phi
$$

*Uses:* A3, A4. *Operation:* mass conservation for each water phase. *Checked by:* none (no tools in this run)

### S4. Mixture enthalpy

The phase field sets the volume fractions, so the enthalpy per volume is the fraction-weighted sum of the ice and vapor enthalpies.

$$
H = \phi\rho_i h_i(T) + (1-\phi)\rho_v h_v(T)
$$

*Uses:* A3. *Operation:* definition of H. *Checked by:* none (no tools in this run)

### S5. Time derivative of H

Product rule, using $dh_i/dT=c_i$ and $dh_v/dT=c_v$.

$$
\partial_t H = h_i\rho_i\partial_t\phi + \phi\rho_i c_i\partial_t T + h_v\partial_t\left[(1-\phi)\rho_v\right] + (1-\phi)\rho_v c_v\partial_t T
$$

*Uses:* S4, A4. *Operation:* differentiate S4 in t (product and chain rule). *Checked by:* none (no tools in this run)

### S6. Eliminate vapor storage

Substitute the vapor balance from S3 for $\partial_t[(1-\phi)\rho_v]$, and expand $\nabla\cdot(h_v\mathbf J_v)=h_v\nabla\cdot\mathbf J_v + c_v\mathbf J_v\cdot\nabla T$. The $h_v\nabla\cdot\mathbf J_v$ terms cancel.

$$
\left[\phi\rho_i c_i + (1-\phi)\rho_v c_v\right]\partial_t T - (h_v-h_i)\rho_i\partial_t\phi + c_v\mathbf J_v\cdot\nabla T + \nabla\cdot\mathbf q = 0
$$

*Uses:* S2, S3, S5. *Operation:* substitution of S3 and S5 into S2. *Checked by:* none (no tools in this run)

### S7. Fourier law, latent heat, neglect vapor enthalpy advection

Introduce $\mathbf q=-k(\phi,T)\nabla T$ (A5), $C(\phi,T)=\phi\rho_i c_i+(1-\phi)\rho_v c_v$, and $L(T)=h_v-h_i$. Drop $c_v\mathbf J_v\cdot\nabla T$: its ratio to conduction is about $10^{-8}$ or less (see neglected).

$$
C(\phi,T)\,\partial_t T = \nabla\cdot\left(k(\phi,T)\nabla T\right) + \rho_i L(T)\,\partial_t\phi
$$

*Uses:* S6, A5, A6. *Operation:* substitution of constitutive law and definitions; drop small term. *Checked by:* none (no tools in this run)

### S8. Residual in 3D Cartesian form

Write the divergence in Cartesian components and move every term to one side. Units: every term is in W m^-3.

$$
R = C\,\partial_t T - \partial_x(k\,\partial_x T) - \partial_y(k\,\partial_y T) - \partial_z(k\,\partial_z T) - \rho_i L\,\partial_t\phi = 0
$$

*Uses:* S7. *Operation:* expand divergence in x, y, z; collect as R. *Checked by:* none (no tools in this run (units by hand: J m^-3 K^-1 * K s^-1 = W m^-3; W m^-1 K^-1 * K m^-2 = W m^-3; kg m^-3 * J kg^-1 * s^-1 = W m^-3))

## Result

$$
R = C(\phi,T)\,\partial_t T - \partial_x\!\left(k(\phi,T)\,\partial_x T\right) - \partial_y\!\left(k(\phi,T)\,\partial_y T\right) - \partial_z\!\left(k(\phi,T)\,\partial_z T\right) - \rho_i L(T)\,\partial_t\phi = 0
$$

Graded SymPy form:

```
C(phi(x, y, z, t), T(x, y, z, t))*Derivative(T(x, y, z, t), t) - Derivative(k(phi(x, y, z, t), T(x, y, z, t))*Derivative(T(x, y, z, t), x), x) - Derivative(k(phi(x, y, z, t), T(x, y, z, t))*Derivative(T(x, y, z, t), y), y) - Derivative(k(phi(x, y, z, t), T(x, y, z, t))*Derivative(T(x, y, z, t), z), z) - rho_i*L(T(x, y, z, t))*Derivative(phi(x, y, z, t), t)
```

## Neglected effects

| Effect | Why it is left out |
|---|---|
| Enthalpy carried by vapor diffusion, $c_v\mathbf J_v\cdot\nabla T$ | $\mathbf J_v\sim D\nabla\rho_v$ with $\rho_v\sim10^{-6}$ kg m^-3 near 200 K. The ratio to conduction is about $c_v D\rho_v/(k\ell)\cdot\ell \sim 1850\cdot10^{-4}\cdot10^{-6}/2.5\sim10^{-7}$, and its gradient factor makes it smaller still. |
| Gas convection (bulk flow of the gas in the pores) | The pores are small (tens of micrometers) and the pressures are low, so the Peclet number is far below 1. |
| Pressure work and viscous dissipation | Pressure changes are tiny compared with $\rho c T$, and the flow is negligible. |
| Radiative heat transfer across pores | With pore size $d\sim10^{-4}$ m, $4\sigma T^3 d\sim2\times10^{-4}$ W m^-1 K^-1, much less than $k$ of ice ($\sim2.5$). |
| Interfacial (surface) energy and Gibbs-Thomson change of $L$ | $\gamma/(\rho_i L r)\sim0.1/(920\cdot2.8\times10^6\cdot10^{-6})\sim4\times10^{-5}$ relative to the latent heat. |
| Thermal expansion of ice (variable $\rho_i$) | The relative change is about $10^{-4}$ per 1 K, a negligible correction to $\rho_i\partial_t\phi$. |
| Non-equilibrium between ice and gas temperatures | The thermal diffusion time across a pore, about $d^2/\alpha\sim10^{-8}/10^{-6}=10^{-2}$ s, is far shorter than sublimation time scales. |

