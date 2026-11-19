# T05_vapor_balance

Run `20261118-173247_T05_vapor_balance_baseline` · status **done** · written by `scripts/agentkit/render_derivation.py` from `answer.json`; do not edit by hand.

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
| $\rho_i$ | density of ice (mass of ice per unit volume of ice), taken constant | $\mathrm{kg\,m^{-3}}$ | 917 | introduced |
| $D_v(T)$ | diffusivity of water vapor in the pore gas | $\mathrm{m^{2}\,s^{-1}}$ | 2e-05 | introduced |
| $m_w$ | water mass (ice plus vapor) per unit total volume | $\mathrm{kg\,m^{-3}}$ | 917 | introduced |
| $\mathbf{j}_v$ | vapor mass flux per unit total area | $\mathrm{kg\,m^{-2}\,s^{-1}}$ | 1e-07 | introduced |
| $\Omega$ | fixed control volume in the pore space, with boundary \partial\Omega and outward normal \mathbf{n} | $\mathrm{m^{3}}$ | 1e-12 | introduced |

## Assumptions

- **A1.** A1: $\phi(x,y,z,t)$ and $T(x,y,z,t)$ are known, smooth fields (given).
- **A2.** A2: Ice is a rigid, immobile solid with constant density $\rho_i$ (no ice velocity, no thermal expansion).
- **A3.** A3: Vapor moves through the gas by Fickian diffusion with diffusivity $D_v(T)$. Any dependence on gas pressure is absorbed into $D_v$.
- **A4.** A4: The gas has no net bulk flow, so there is no advection term.
- **A5.** A5: No water enters or leaves except through the vapor flux, and there is no liquid water. Water mass is conserved.

## Derivation

### S1. Water mass per unit volume

Ice fills the fraction $\phi$ of each small volume and gas fills $1-\phi$. So the water mass per unit total volume is the ice mass plus the vapor mass held in the gas.

$$
m_w = \rho_i\phi + (1-\phi)\rho_v
$$

*Uses:* given, A1, A2. *Operation:* definition of mixture mass density from volume fractions. *Checked by:* none (no tools in this run; by hand, both terms are in kg m^-3)

### S2. Integral conservation of water

Water is neither created nor destroyed (A5), so in a fixed control volume its mass changes only through flux across the boundary. Ice does not move (A2), so only vapor carries a flux.

$$
\frac{d}{dt}\int_\Omega \left[\rho_i\phi + (1-\phi)\rho_v\right]dV = -\oint_{\partial\Omega}\mathbf{j}_v\cdot\mathbf{n}\,dS
$$

*Uses:* S1, A2, A5. *Operation:* mass balance on a fixed control volume. *Checked by:* none (no tools in this run; by hand, both sides are in kg s^-1)

### S3. Local conservation law

The fields are smooth (A1), so the time derivative can move inside the integral and the divergence theorem turns the surface integral into a volume integral. Because $\Omega$ is arbitrary, the integrand must vanish.

$$
\frac{\partial}{\partial t}\left[(1-\phi)\rho_v\right] + \rho_i\frac{\partial\phi}{\partial t} + \nabla\cdot\mathbf{j}_v = 0
$$

*Uses:* S2, A1. *Operation:* divergence theorem and arbitrariness of the control volume. *Checked by:* none (no tools in this run; by hand, every term is in kg m^-3 s^-1. The constant $\rho_i$ (A2) is pulled out of $\partial_t(\rho_i\phi)$)

### S4. Fickian vapor flux in the gas fraction

Vapor diffuses only through the gas, which fills the fraction $1-\phi$ of each volume. With Fick's law in the gas (A3) and no net gas flow (A4), the flux per unit total area is the gas fraction times the Fickian flux.

$$
\mathbf{j}_v = -(1-\phi)D_v(T)\nabla\rho_v
$$

*Uses:* A3, A4. *Operation:* constitutive law (Fick) weighted by gas fraction. *Checked by:* none (no tools in this run; by hand, units are m^2 s^-1 times kg m^-4, which is kg m^-2 s^-1)

### S5. Vapor equation as R = 0

Substituting S4 into S3 and writing the divergence in Cartesian components gives the residual. The term $\rho_i\partial_t\phi$ is the exchange with ice: deposition ($\partial_t\phi>0$) removes vapor and sublimation adds it.

$$
R = \partial_t[(1-\phi)\rho_v] + \rho_i\partial_t\phi - \sum_{x_k\in\{x,y,z\}}\partial_{x_k}\left[(1-\phi)D_v(T)\partial_{x_k}\rho_v\right] = 0
$$

*Uses:* S3, S4. *Operation:* substitution of S4 into S3, Cartesian divergence. *Checked by:* none (no tools in this run; by hand, every term is in kg m^-3 s^-1)

## Result

$$
R = \frac{\partial}{\partial t}\left[(1-\phi)\rho_v\right] + \rho_i\frac{\partial \phi}{\partial t} - \frac{\partial}{\partial x}\left[(1-\phi)D_v(T)\frac{\partial \rho_v}{\partial x}\right] - \frac{\partial}{\partial y}\left[(1-\phi)D_v(T)\frac{\partial \rho_v}{\partial y}\right] - \frac{\partial}{\partial z}\left[(1-\phi)D_v(T)\frac{\partial \rho_v}{\partial z}\right] = 0
$$

Graded SymPy form:

```
Derivative((1 - phi(x, y, z, t))*rho_v(x, y, z, t), t) + rho_i*Derivative(phi(x, y, z, t), t) - Derivative((1 - phi(x, y, z, t))*D_v(T(x, y, z, t))*Derivative(rho_v(x, y, z, t), x), x) - Derivative((1 - phi(x, y, z, t))*D_v(T(x, y, z, t))*Derivative(rho_v(x, y, z, t), y), y) - Derivative((1 - phi(x, y, z, t))*D_v(T(x, y, z, t))*Derivative(rho_v(x, y, z, t), z), z)
```

## Neglected effects

| Effect | Why it is left out |
|---|---|
| Advection of vapor by bulk gas flow (including the Stefan flow at the interface) | The vapor is very dilute near 200 K (vapor pressure about 0.16 Pa, so $\rho_v$ is about 1e-6 kg m^-3), so the Stefan velocity is about $D_v\nabla\rho_v/\rho_{gas}$ and is negligible compared with diffusion. There is no imposed pressure gradient across the pores. |
| Thermal diffusion (the Soret effect) and cross-diffusion | The thermal diffusion factor for water in air is about 0.1 or less. Over pore-scale temperature differences of about 1 K at 200 K, it contributes a flux of order 1e-3 of the Fickian flux. |
| Knudsen (rarefied-gas) transport | This is valid only while the mean free path is smaller than the pore size. Under lunar near-vacuum it is not, and then $D_v$ must be replaced by a Knudsen diffusivity. The form of the equation stays the same. |
| Variation of ice density with temperature, and ice motion | The thermal expansion of ice is about 5e-5 per K, so changes in $\rho_i$ are about 1e-4 of $\rho_i$. Grains are assumed fixed. |
| Liquid water and adsorbed water on regolith grains | Near 200 K there is no liquid phase. Adsorption on the grain surfaces is outside the ice/gas model. |

## Questions

- Should the vapor flux carry the gas-fraction weighting $(1-\phi)$, as used here so that water mass is conserved exactly? The alternative is the common simplified form $\partial_t\rho_v = \nabla\cdot(D_v\nabla\rho_v) - \rho_i\partial_t\phi$.

