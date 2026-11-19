# T05_vapor_balance

Run `20261118-162055_T05_vapor_balance` · status **done** · written by `scripts/agentkit/render_derivation.py` from `answer.json`; do not edit by hand.

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
| $\rho_i(T)$ | density of ice (mass of ice per unit volume of ice) | $\mathrm{kg\,m^{-3}}$ | 930 | introduced |
| $D(\phi, T)$ | effective diffusivity of water vapor through the pore gas, per unit total cross-section (includes the gas fraction, tortuosity and any Knudsen correction; tends to 0 in the ice) | $\mathrm{m^{2}\,s^{-1}}$ | 0.0001 | introduced |
| $m_w$ | total water mass (ice plus vapor) per unit total volume | $\mathrm{kg\,m^{-3}}$ | 1 | introduced |
| $\mathbf{j}_v$ | mass flux of water vapor per unit total area | $\mathrm{kg\,m^{-2}\,s^{-1}}$ | 1e-06 | introduced |
| $\Omega$ | arbitrary fixed control volume with boundary $\partial\Omega$ and outward unit normal $\mathbf{n}$ | $\mathrm{m^{3}}$ | 1e-15 | introduced |
| $\sigma$ | mass rate of sublimation (ice to vapor) per unit total volume | $\mathrm{kg\,m^{-3}\,s^{-1}}$ | 1e-06 | introduced |

## Assumptions

- **A1.** A1: $\phi$ and $T$ are known, smooth fields (given).
- **A2.** A2: The pore space holds only ice and gas, so the volume fractions are $\phi$ and $1-\phi$; $\rho_v$ is per unit gas volume and the ice density $\rho_i(T)$ is per unit ice volume. Grains of regolith are excluded from the domain (or enter only through $D$).
- **A3.** A3: Ice is rigid and does not move; it gains or loses water only by deposition and sublimation.
- **A4.** A4: The gas has no bulk (advective or Stefan) flow; vapor moves only by diffusion.
- **A5.** A5: Vapor diffusion obeys Fick's law with an effective diffusivity $D(\phi,T)$ per unit total area, which includes the gas fraction, tortuosity and any Knudsen correction and goes to zero in ice.

## Derivation

### S1. Water mass per unit volume

Water is present as ice in the fraction $\phi$ and as vapor in the gas fraction $1-\phi$; $\rho_v$ is per unit gas volume, so it is weighted by $1-\phi$.

$$
m_w = \phi\,\rho_i(T) + (1-\phi)\,\rho_v
$$

*Uses:* given, A1, A2. *Operation:* volume-fraction weighted sum. *Checked by:* `check_s1` in `scratch/check.py`: S1 phi=1: rho_i(T)  phi=0: rho_v

### S2. Integral conservation of water

Water is neither created nor destroyed; in a fixed control volume its mass changes only by the flux through the boundary. Ice does not move (A3), so only vapor carries water across $\partial\Omega$.

$$
\frac{d}{dt}\int_\Omega m_w\,dV = -\oint_{\partial\Omega} \mathbf{j}_v\cdot\mathbf{n}\,dA
$$

*Uses:* S1, A3. *Operation:* conservation of mass over a fixed control volume. *Checked by:* none (Statement of the conservation principle; nothing to compute.)

### S3. Fickian vapor flux

With no bulk gas flow (A4), vapor moves by diffusion down its density gradient; $D$ is an effective diffusivity per unit total area that depends on $\phi$ and $T$ and vanishes in the ice.

$$
\mathbf{j}_v = -D(\phi,T)\,\nabla\rho_v
$$

*Uses:* A4, A5. *Operation:* constitutive law (Fick). *Checked by:* `check_units` in `scratch/check.py`: d_x(D d_x rho_v) has units kilogram / meter ** 3 / second, the same as d_t((1-phi) rho_v)

### S4. Local form

The divergence theorem turns the surface integral into a volume integral; the volume is fixed and arbitrary and the fields are smooth (A1), so the integrand vanishes pointwise.

$$
\frac{\partial}{\partial t}\big[(1-\phi)\rho_v\big] + \frac{\partial}{\partial t}\big[\phi\rho_i(T)\big] - \nabla\cdot\big(D(\phi,T)\nabla\rho_v\big) = 0
$$

*Uses:* S2, S3, A1. *Operation:* divergence theorem, then arbitrary $\Omega$. *Checked by:* `check_s4` in `scratch/check.py`: S4 int of div J = 0 (1D divergence-theorem test)

### S5. Exchange term made explicit

Splitting into a vapor balance and an ice balance shows the sublimation source $\sigma = -\partial_t(\phi\rho_i)$ that the ice loses and the vapor gains; with constant $\rho_i$ this is $-\rho_i\partial_t\phi$, and the vapor equation reads $(1-\phi)\partial_t\rho_v + (\rho_i-\rho_v)\partial_t\phi - \nabla\cdot(D\nabla\rho_v)=0$.

$$
\frac{\partial}{\partial t}\big[(1-\phi)\rho_v\big] - \nabla\cdot\big(D\nabla\rho_v\big) = \sigma, \qquad \sigma = -\frac{\partial}{\partial t}\big[\phi\rho_i(T)\big]
$$

*Uses:* S4. *Operation:* split into ice and vapor balances, define $\sigma$. *Checked by:* `check_s5` in `scratch/check.py`: S5 R(const rho_i) - familiar = 0

### S6. Residual in 3D Cartesian coordinates

Write the divergence and gradient of S4 component by component; R is the left side of S4.

$$
R = \partial_t\big[(1-\phi)\rho_v\big] + \partial_t\big[\phi\rho_i(T)\big] - \partial_x(D\,\partial_x\rho_v) - \partial_y(D\,\partial_y\rho_v) - \partial_z(D\,\partial_z\rho_v) = 0
$$

*Uses:* S4. *Operation:* expand $\nabla\cdot(D\nabla\rho_v)$ in Cartesian components. *Checked by:* `check_s5` in `scratch/check.py`: S5 R expanded printed in full; every term has units kg m^-3 s^-1 (check_units)

## Result

$$
R = \frac{\partial}{\partial t}\big[(1-\phi)\rho_v\big] + \frac{\partial}{\partial t}\big[\phi\,\rho_i(T)\big] - \frac{\partial}{\partial x}\Big(D(\phi,T)\frac{\partial \rho_v}{\partial x}\Big) - \frac{\partial}{\partial y}\Big(D(\phi,T)\frac{\partial \rho_v}{\partial y}\Big) - \frac{\partial}{\partial z}\Big(D(\phi,T)\frac{\partial \rho_v}{\partial z}\Big) = 0
$$

Graded SymPy form:

```
Derivative((1 - phi(x, y, z, t))*rho_v(x, y, z, t), t) + Derivative(phi(x, y, z, t)*rho_i(T(x, y, z, t)), t) - Derivative(D(phi(x, y, z, t), T(x, y, z, t))*Derivative(rho_v(x, y, z, t), x), x) - Derivative(D(phi(x, y, z, t), T(x, y, z, t))*Derivative(rho_v(x, y, z, t), y), y) - Derivative(D(phi(x, y, z, t), T(x, y, z, t))*Derivative(rho_v(x, y, z, t), z), z)
```

## Neglected effects

| Effect | Why it is left out |
|---|---|
| Advection of vapor by bulk gas flow, including Stefan flow caused by sublimation | Vapor is a trace species or the gas is rarefied; the Stefan velocity is about $\sigma/\rho_{gas}$ and gives a Peclet number $vL/D\ll 1$ at pore scale ($L\sim10^{-5}$ m). |
| Thermal diffusion (Soret effect) | The Soret coefficient for water in air is about $10^{-3}$ K$^{-1}$ or smaller, so it is small next to the Fickian flux from the much steeper dependence of saturation density on $T$. |
| Surface diffusion of adsorbed water on grains and ice | Adsorbed mass on the grains is not part of the model; its flux is small next to gas-phase transport for pore sizes much larger than the adsorbed layer. |
| Motion and creep of the ice | Ice creep velocities at 200 K are many orders of magnitude slower than sublimation fronts. |
| Temperature dependence of ice density | Kept in the result through $\rho_i(T)$; ice thermal expansion is about $1.5\times10^{-4}$ K$^{-1}$, so setting $\rho_i$ constant changes the result by under 0.1% per 10 K if a later task wants it. |
| Water bound in or adsorbed on regolith grains | Outside the ice-gas pore model by the statement. |

