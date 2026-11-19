# T05_vapor_balance

Run `20261118-172817_T05_vapor_balance` · status **done** · written by `scripts/agentkit/render_derivation.py` from `answer.json`; do not edit by hand.

## Work order (verbatim)

> We are modeling the pore space of a granular material (lunar regolith) in which water ice sublimates and deposits. Treat the pore space as a three-dimensional region filled with solid ice and gas (air that contains water vapor). The local volume fraction of ice is the phase field phi(x, y, z, t), between 0 and 1: phi = 1 in ice and phi = 0 in gas, with a thin, smooth transition between the two. Ice turns into vapor (sublimation) and vapor into ice (deposition), so phi changes in time. The vapor density rho_v(x, y, z, t) is the mass of water vapor per unit volume of gas. In this task phi and the temperature T are known; their own equations are separate tasks. Derive the partial differential equation for rho_v that expresses conservation of water mass, including the exchange of mass between ice and vapor. Start from a conservation principle. State every assumption, and list every effect you leave out with the reason it is small. Keep the material properties general: do not assume a property is constant unless a step needs it, and say so when one does. Write the result as R = 0 in 3D Cartesian coordinates and return R.

## Symbols

| Symbol | Meaning | Units | Typical value (SI) | Source |
|---|---|---|---|---|
| $x$ | first position coordinate | $\mathrm{m}$ |  | work order |
| $y$ | second position coordinate | $\mathrm{m}$ |  | work order |
| $z$ | third position coordinate | $\mathrm{m}$ |  | work order |
| $t$ | time | $\mathrm{s}$ |  | work order |
| $T$ | temperature | $\mathrm{K}$ |  | work order |
| $\phi$ | local volume fraction of ice (phase field), between 0 and 1 | 1 |  | work order |
| $\rho_v$ | vapor density: mass of water vapor per unit volume of gas | $\mathrm{kg\,m^{-3}}$ |  | work order |
| $\rho_i(T)$ | density of ice (mass of water per unit volume of ice) | $\mathrm{kg\,m^{-3}}$ | 930 | introduced |
| $D_v(T)$ | diffusivity of water vapor in the pore gas | $\mathrm{m^{2}\,s^{-1}}$ | 2e-05 | introduced |
| $m$ | total water mass (ice plus vapor) per unit total volume | $\mathrm{kg\,m^{-3}}$ | 900 | introduced |
| $\mathbf{J}$ | vapor mass flux per unit total area | $\mathrm{kg\,m^{-2}\,s^{-1}}$ | 1e-08 | introduced |
| $V$ | fixed control volume in the pore space, with boundary \partial V and outward unit normal \mathbf{n} | $\mathrm{m^{3}}$ | 1e-15 | introduced |
| $R$ | residual of the water mass balance | $\mathrm{kg\,m^{-3}\,s^{-1}}$ | 0 | introduced |

## Assumptions

- **A1.** A1: $\phi(x,y,z,t)$ and $T(x,y,z,t)$ are known, smooth fields (given).
- **A2.** A2: Water mass is conserved; it changes phase between ice and vapor but has no other sources or sinks (no chemical binding to regolith, no adsorption).
- **A3.** A3: The ice is a rigid, immobile solid with density $\rho_i(T)$ (no ice velocity, no deformation or sintering flow).
- **A4.** A4: Vapor is transported only by Fickian diffusion in the gas fraction with diffusivity $D_v(T)$; the gas has no bulk (advective) velocity.
- **A5.** A5: The effective diffusive flux in the diffuse interface is weighted by the gas fraction $1-\phi$, so no vapor flux passes through ice.

## Derivation

### S1. Integral water mass balance

Water is neither created nor destroyed; it only changes phase. So the water mass in a fixed control volume $V$ changes only by the flux $\mathbf{J}$ through its boundary.

$$
\frac{d}{dt}\int_V m\,dV = -\oint_{\partial V}\mathbf{J}\cdot\mathbf{n}\,dS
$$

*Uses:* given, A2. *Operation:* conservation of mass for a fixed control volume. *Checked by:* none (Statement of the conservation principle; no tool check.)

### S2. Water mass per unit total volume

A fraction $\phi$ of each volume is ice with density $\rho_i$, and a fraction $1-\phi$ is gas carrying vapor at density $\rho_v$ per unit gas volume. Both phases hold water, so both enter $m$; this is where the ice-vapor exchange is accounted for.

$$
m = \rho_i(T)\,\phi + (1-\phi)\,\rho_v
$$

*Uses:* A1, A3. *Operation:* volume-fraction weighting of the phase densities. *Checked by:* `check_s2` in `scratch/check.py`: S2 m at phi=1: rho_i(T(x, y, z, t))  at phi=0: rho_v(x, y, z, t)

### S3. Vapor flux

Ice is rigid and does not move, so it carries no flux. Vapor moves by Fickian diffusion in the gas only; the flux per unit total area is weighted by the gas fraction $1-\phi$, which makes it vanish inside ice.

$$
\mathbf{J} = -(1-\phi)\,D_v(T)\,\nabla\rho_v
$$

*Uses:* A3, A4, A5. *Operation:* Fick's law applied to the gas fraction. *Checked by:* `check_units` in `scratch/check.py`: units div((1-phi) D grad rho_v) kilogram / meter ** 3 / second

### S4. Local balance

Apply the divergence theorem to the surface integral of S1; since $V$ is arbitrary and the fields are smooth (A1), the integrand must vanish pointwise. Substituting S2 and S3 gives the residual.

$$
R = \partial_t\left[(1-\phi)\rho_v\right] + \partial_t\left[\rho_i(T)\phi\right] - \nabla\cdot\left[(1-\phi)D_v(T)\nabla\rho_v\right] = 0
$$

*Uses:* S1, S2, S3, A1. *Operation:* divergence theorem, arbitrary V, substitution of S2 and S3. *Checked by:* `check_s4` in `scratch/check.py`: S4 local form minus R: 0; units of all three terms: kilogram / meter ** 3 / second

### S5. Expanded form with constant ice density (optional)

Only if $\rho_i$ is taken constant, the product rule shows the familiar form: the exchange term $(\rho_i-\rho_v)\partial_t\phi$ is the source of vapor from sublimation ($\partial_t\phi<0$). This step is a rewriting for the reader; the result keeps $\rho_i(T)$ general.

$$
(1-\phi)\partial_t\rho_v + (\rho_i-\rho_v)\partial_t\phi - \nabla\cdot\left[(1-\phi)D_v\nabla\rho_v\right] = 0
$$

*Uses:* S4. *Operation:* product rule on the time derivatives, with constant $\rho_i$. *Checked by:* `check_s5` in `scratch/check.py`: S5 expanded minus R (rho_i const): 0

### S6. Cartesian form

Write the divergence and gradient of S4 in 3D Cartesian coordinates to give the requested residual.

$$
R = \partial_t[(1-\phi)\rho_v] + \partial_t[\rho_i\phi] - \sum_{\xi\in\{x,y,z\}}\partial_\xi\left[(1-\phi)D_v(T)\partial_\xi\rho_v\right]
$$

*Uses:* S4. *Operation:* expand divergence in Cartesian coordinates. *Checked by:* `check_s4` in `scratch/check.py`: R printed in expanded Cartesian form; local form minus R: 0

## Result

$$
R = \partial_t\left[(1-\phi)\rho_v\right] + \partial_t\left[\rho_i(T)\,\phi\right] - \partial_x\left[(1-\phi)D_v(T)\,\partial_x\rho_v\right] - \partial_y\left[(1-\phi)D_v(T)\,\partial_y\rho_v\right] - \partial_z\left[(1-\phi)D_v(T)\,\partial_z\rho_v\right] = 0
$$

Graded SymPy form:

```
Derivative((1 - phi(x, y, z, t))*rho_v(x, y, z, t), t) + Derivative(rho_i(T(x, y, z, t))*phi(x, y, z, t), t) - Derivative((1 - phi(x, y, z, t))*D_v(T(x, y, z, t))*Derivative(rho_v(x, y, z, t), x), x) - Derivative((1 - phi(x, y, z, t))*D_v(T(x, y, z, t))*Derivative(rho_v(x, y, z, t), y), y) - Derivative((1 - phi(x, y, z, t))*D_v(T(x, y, z, t))*Derivative(rho_v(x, y, z, t), z), z)
```

## Neglected effects

| Effect | Why it is left out |
|---|---|
| Advection of vapor by bulk gas flow, including Stefan flow caused by sublimation | The Stefan velocity is about $\rho_v/\rho_{gas}$ times the diffusive velocity. Near 200 K, $\rho_v \sim 10^{-5}$ kg m$^{-3}$, which is small, and no imposed pressure gradient drives flow in the pores. |
| Thermal diffusion (Soret effect) and pressure diffusion | The Soret coefficient for water vapor in air is of order 0.01-0.1 of the Fick flux per relative temperature change; with pore-scale $\Delta T/T \ll 1$ it is a small correction. |
| Knudsen or transition-regime transport in near-vacuum lunar pores | This is not strictly neglected: it changes the value and the dependence of $D_v$, which can depend on pore size as well as $T$. The form of the balance stays the same, but $D_v$ must be chosen for the gas-pressure regime. |
| Vapor adsorption on regolith grain surfaces and surface diffusion | The work order defines only ice and gas in the pore space, so grains are outside the domain; adsorbed water would need its own balance. |
| Ice motion (creep, sintering flow) and density change with temperature | Ice velocities are negligible on the time scales of sublimation. The thermal expansion of ice is about $10^{-4}$ K$^{-1}$; $\rho_i(T)$ is kept in the result anyway, and it can be set constant. |
| Anisotropy or tortuosity of diffusion in the diffuse interface beyond the $(1-\phi)$ weighting | These affect only the interface layer, which is thin compared with the pore size, so the error is of the order of the interface width over the pore size. |

