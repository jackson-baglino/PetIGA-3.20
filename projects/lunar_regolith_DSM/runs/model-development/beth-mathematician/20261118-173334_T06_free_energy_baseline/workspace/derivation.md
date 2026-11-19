# T06_free_energy

Run `20261118-173334_T06_free_energy_baseline` · status **done** · written by `scripts/agentkit/render_derivation.py` from `answer.json`; do not edit by hand.

## Work order (verbatim)

> We are modeling the pore space of a granular material (lunar regolith) in which water ice sublimates and deposits. Treat the pore space as a three-dimensional region filled with solid ice and gas. The local volume fraction of ice is the phase field phi(x, y, z, t), between 0 and 1: phi = 1 in ice and phi = 0 in gas, with a thin, smooth transition between the two. Construct the free energy of this system as a functional of the phase field, F[phi] = integral of f over the volume, at phase coexistence (constant temperature, no driving force for growth or sublimation). It must make bulk ice and bulk gas the stable states, and give an interface between them that has a finite width and an energy per unit area. Choose the form of f, and write its coefficients in terms of gamma and eps only (and numbers), so that a flat interface at equilibrium has energy per unit area exactly gamma and a width proportional to eps; say what multiple of eps the width is, and how you define the width. Explain why each term is there. State every assumption, and list every effect you leave out with the reason. Return the free energy density f.

## Symbols

| Symbol | Meaning | Units | Typical value (SI) | Source |
|---|---|---|---|---|
| $x$ | first position coordinate | $\mathrm{m}$ |  | work order |
| $y$ | second position coordinate | $\mathrm{m}$ |  | work order |
| $z$ | third position coordinate | $\mathrm{m}$ |  | work order |
| $t$ | time | $\mathrm{s}$ |  | work order |
| $\phi$ | local volume fraction of ice (phase field), between 0 and 1 | 1 |  | work order |
| $\gamma$ | energy per unit area of the ice-gas interface, constant | $\mathrm{J\,m^{-2}}$ |  | work order |
| $\varepsilon$ | a length that sets the width of the diffuse interface, constant | $\mathrm{m}$ |  | work order |
| $\kappa$ | gradient-energy coefficient | $\mathrm{J\,m^{-1}}$ | 2e-10 | introduced |
| $W$ | height scale of the double-well potential | $\mathrm{J\,m^{-3}}$ | 2e+05 | introduced |
| $\xi$ | coordinate normal to a flat interface, with the interface at xi = 0 | $\mathrm{m}$ | 1e-06 | introduced |
| $\phi_0(\xi)$ | equilibrium profile across a flat interface | 1 | 0.5 | introduced |
| $\delta$ | interface width, defined as 1/\|phi_0'(0)\| | $\mathrm{m}$ | 4e-06 | introduced |

## Assumptions

- **A1.** A1: the temperature is uniform and constant, at phase coexistence, so the bulk ice and bulk gas free energies are equal and the double well is symmetric.
- **A2.** A2: the interface energy $\gamma$ is isotropic, and the free energy is local, depending on $\phi$ and $\nabla\phi$ only (square-gradient approximation).
- **A3.** A3: the width scale $\varepsilon$ is defined by the tail decay $\phi_0\sim e^{-\|\xi\|/\varepsilon}$, which gives the reported width $\delta=1/\|\phi_0'(0)\|=4\varepsilon$.

## Derivation

### S1. Ansatz for the free energy density

The density is a double well, which makes $\phi=0$ and $\phi=1$ the only stable uniform states. The wells have equal depth because the temperature is at coexistence (A1), so there is no tilt term. A gradient term penalizes sharp changes in $\phi$, which gives the interface a finite width and a positive energy.

$$
f = \frac{\kappa}{2}|\nabla\phi|^2 + W\,\phi^2(1-\phi)^2
$$

*Uses:* A1, A2. *Operation:* choose the Ginzburg-Landau form. *Checked by:* none (no tools in this run; units: $\kappa|\nabla\phi|^2$ has J m^-1 m^-2 = J m^-3, and $W$ has J m^-3)

### S2. First integral for a flat interface

For a flat interface, $\phi$ depends only on $\xi$. The Euler-Lagrange equation $\kappa\phi''=W\,d[\phi^2(1-\phi)^2]/d\phi$ has a first integral. Its constant is zero because $\phi'\to0$ and the well term vanishes in both bulk phases.

$$
\frac{\kappa}{2}\phi_0'^2 = W\,\phi_0^2(1-\phi_0)^2
$$

*Uses:* S1. *Operation:* multiply the Euler-Lagrange equation by $\phi'$ and integrate once in $\xi$. *Checked by:* none (no tools in this run)

### S3. Equilibrium profile

Take the root with ice at $\xi\to-\infty$ and gas at $\xi\to+\infty$, then separate variables.

$$
\phi_0(\xi)=\frac12\left[1-\tanh\left(\frac{\xi}{2}\sqrt{\frac{2W}{\kappa}}\right)\right]
$$

*Uses:* S2. *Operation:* take the square root and integrate. *Checked by:* none (no tools in this run; by hand, $\phi_0'=-\sqrt{2W/\kappa}\,\phi_0(1-\phi_0)$ holds)

### S4. Interface energy

The excess energy per unit area is the integral of $f$ across the interface. By S2 the two terms are equal, so it equals $\int\kappa\phi_0'^2\,d\xi$. Changing the variable to $\phi$ gives a closed form.

$$
\gamma=\int_0^1\sqrt{2W\kappa}\,\phi(1-\phi)\,d\phi=\frac{\sqrt{2W\kappa}}{6}
$$

*Uses:* S2, S3. *Operation:* integrate over $\xi$ and change the variable to $\phi$. *Checked by:* none (no tools in this run; $\int_0^1\phi(1-\phi)d\phi=1/6$)

### S5. Fix the length scale

Define $\varepsilon$ as the decay length of the profile's tails, $\phi_0\sim e^{-|\xi|/\varepsilon}$. This sets $\sqrt{2W/\kappa}=1/\varepsilon$. The width is $\delta=1/|\phi_0'(0)|=4\varepsilon$; the 10-90% width is $2\ln 9\,\varepsilon$.

$$
\sqrt{\frac{2W}{\kappa}}=\frac{1}{\varepsilon},\qquad \delta=\frac{1}{|\phi_0'(0)|}=4\varepsilon
$$

*Uses:* S3. *Operation:* definition of the width. *Checked by:* none (no tools in this run; $\phi_0'(0)=-1/(4\varepsilon)$)

### S6. Solve for the coefficients

Multiplying and dividing S4 and S5 gives $2W=6\gamma/\varepsilon$ and $2/\kappa=1/(3\gamma\varepsilon)$.

$$
W=\frac{3\gamma}{\varepsilon},\qquad \kappa=6\gamma\varepsilon
$$

*Uses:* S4, S5. *Operation:* solve two equations for $W$ and $\kappa$. *Checked by:* none (no tools in this run; check: $\sqrt{2\cdot3\gamma/\varepsilon\cdot6\gamma\varepsilon}/6=\gamma$)

### S7. Free energy density

Substitute S6 into S1.

$$
f = 3\gamma\varepsilon\,|\nabla\phi|^2 + \frac{3\gamma}{\varepsilon}\,\phi^2(1-\phi)^2
$$

*Uses:* S1, S6. *Operation:* substitution of S6. *Checked by:* none (no tools in this run; units: J m^-2 m m^-2 = J m^-3, and J m^-2 m^-1 = J m^-3)

## Result

$$
f = 3\gamma\varepsilon\,|\nabla\phi|^2 + \frac{3\gamma}{\varepsilon}\,\phi^2(1-\phi)^2
$$

Graded SymPy form:

```
3*gamma*eps*(Derivative(phi(x, y, z, t), x)**2 + Derivative(phi(x, y, z, t), y)**2 + Derivative(phi(x, y, z, t), z)**2) + 3*gamma/eps*phi(x, y, z, t)**2*(1 - phi(x, y, z, t))**2
```

## Neglected effects

| Effect | Why it is left out |
|---|---|
| Driving force (supersaturation or temperature offset) | By A1 the system is at coexistence; this would add a tilt term such as $\propto\phi^2(3-2\phi)$. |
| Vapor density and heat fields | They are not part of this functional at coexistence; they couple in only away from equilibrium. |
| Anisotropy of the ice surface energy | $\gamma$ is taken to be constant, as stated. |
| Contact with the regolith grains (wetting, grain surface energies) | The pore space is treated as ice plus gas only; boundary terms would add surface integrals. |
| Elastic strain and curvature (higher-gradient) terms | These are higher order in $\varepsilon$ over the feature size, and second derivatives are excluded by the work order. |

## Questions

- The coefficients depend on how $\varepsilon$ is tied to the width. I used the tail decay length, giving $\tanh(\xi/2\varepsilon)$ and $\delta=4\varepsilon$. If the project uses $\tanh(\xi/\varepsilon)$ instead, the result becomes $f=\tfrac32\gamma\varepsilon\|\nabla\phi\|^2+\tfrac{6\gamma}{\varepsilon}\phi^2(1-\phi)^2$. Please confirm the convention.

