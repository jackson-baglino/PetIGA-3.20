# T06_free_energy

Run `20261118-172917_T06_free_energy` · status **done** · written by `scripts/agentkit/render_derivation.py` from `answer.json`; do not edit by hand.

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
| $W(\phi)$ | dimensionless double-well function, \phi^2(1-\phi)^2 | 1 | 0.0625 | introduced |
| $A$ | height coefficient of the double-well term (before fixing it) | $\mathrm{J\,m^{-3}}$ | 3e+05 | introduced |
| $K$ | gradient-energy coefficient (before fixing it) | $\mathrm{J\,m^{-1}}$ | 6e-08 | introduced |
| $n$ | coordinate normal to a flat interface | $\mathrm{m}$ | 1e-06 | introduced |
| $\phi_0(n)$ | equilibrium profile of a flat interface | 1 | 0.5 | introduced |
| $\delta$ | interface width, defined as 1/\phi_0'(0) (tangent at the midpoint) | $\mathrm{m}$ | 4e-06 | introduced |
| $F[\phi]$ | free energy functional | $\mathrm{J}$ | 1e-12 | introduced |

## Assumptions

- **A1.** A1: The temperature is uniform and constant, at phase coexistence, so the two bulk minima have equal free energy (no tilt term).
- **A2.** A2: $\gamma$ is isotropic and constant, so the gradient term is $\|\nabla\phi\|^2$ with a scalar coefficient.
- **A3.** A3: Quartic double well $W=\phi^2(1-\phi)^2$ (the simplest smooth polynomial with minima at 0 and 1); width defined as $\delta = 1/\phi_0'(0)$, the inverse slope at the midpoint.

## Derivation

### S1. Ansatz: double well plus gradient term

The bulk term $A\,W(\phi)$ with $W=\phi^2(1-\phi)^2$ has equal minima at $\phi=0$ (gas) and $\phi=1$ (ice), so both bulk states are stable and coexist (equal depth, no driving force). The barrier $W(1/2)=1/16$ penalizes mixed states, and the gradient term $\tfrac{K}{2}|\nabla\phi|^2$ penalizes sharp jumps; their balance gives a finite width and a finite energy per area.

$$
F[\phi]=\int_V \Big[A\,\phi^2(1-\phi)^2+\tfrac{K}{2}|\nabla\phi|^2\Big]\,dV
$$

*Uses:* given, A1, A2, A3. *Operation:* choice of the Ginzburg-Landau form. *Checked by:* `check_s1` in `scratch/check.py`: W' roots 1/2, 1 (and 0, dropped by a positive-symbol flag); W''(0)=W''(1)=2>0; W(1/2)=1/16

### S2. Euler-Lagrange equation for a flat interface

At equilibrium $\delta F/\delta\phi=0$. For a flat interface $\phi$ depends only on the normal coordinate $n$.

$$
K\,\phi_0''(n) = A\,W'(\phi_0),\quad \phi_0(-\infty)=0,\ \phi_0(\infty)=1
$$

*Uses:* S1. *Operation:* first variation of F, restricted to 1D. *Checked by:* `check_s3` in `scratch/check.py`: used inside check_s3 (sympy euler_equations)

### S3. First integral (equipartition)

Multiply S2 by $\phi_0'$ and integrate; the constant is zero because $\phi_0'\to0$ and $W\to0$ in the bulk.

$$
\tfrac{K}{2}\phi_0'^2 = A\,W(\phi_0)\ \Rightarrow\ \phi_0'=\sqrt{2A/K}\,\phi_0(1-\phi_0)
$$

*Uses:* S2. *Operation:* multiply by $\phi_0'$, integrate over n. *Checked by:* `check_s3` in `scratch/check.py`: dH/dx + phi' * EL = 0

### S4. Profile

Separating variables in S3 gives the logistic (tanh) profile centered at $n=0$.

$$
\phi_0(n)=\frac{1}{1+e^{-n\sqrt{2A/K}}}
$$

*Uses:* S3. *Operation:* separation of variables. *Checked by:* `check_s4` in `scratch/check.py`: EL residual of logistic profile: 0

### S5. Interface energy

The excess energy per area is $\int (A W + \tfrac K2\phi_0'^2)dn = \int 2AW\,dn$ by S3; change variable $dn=d\phi/\phi_0'$.

$$
\gamma=\int_0^1\sqrt{2AK}\,\phi(1-\phi)\,d\phi=\frac{\sqrt{2AK}}{6}
$$

*Uses:* S3. *Operation:* substitution $dn = d\phi/\phi_0'$. *Checked by:* `check_s5` in `scratch/check.py`: gamma = sqrt(2)*sqrt(A)*sqrt(K)/6

### S6. Fix A and K

Require energy $\gamma$ (S5) and set the profile's inverse decay length $\sqrt{2A/K}=1/\varepsilon$; solve the two equations.

$$
A=\frac{3\gamma}{\varepsilon},\quad K=6\gamma\varepsilon
$$

*Uses:* S4, S5. *Operation:* solve two equations for A, K. *Checked by:* `check_s6` in `scratch/check.py`: A: 3*gamma/epsilon, K: 6*epsilon*gamma

### S7. Final density, profile and width

Substitute S6 into S1. The profile is $\phi_0=\tfrac12[1+\tanh(n/2\varepsilon)]$; with width defined as $\delta=1/\phi_0'(0)$ (tangent at the midpoint) we get $\delta=4\varepsilon$ (the 10-90% width is $4\ln3\,\varepsilon\approx4.39\,\varepsilon$). Integrating f across the profile returns $\gamma$; both terms have units J m$^{-3}$.

$$
f=\frac{3\gamma}{\varepsilon}\phi^2(1-\phi)^2+3\gamma\varepsilon|\nabla\phi|^2,\quad \delta=4\varepsilon
$$

*Uses:* S1, S6. *Operation:* substitution of S6. *Checked by:* `check_s7` in `scratch/check.py`: profile == (1+tanh(x/2eps))/2: 0; width 1/phi'(0) = 4*epsilon; 10-90 width = 4*epsilon*log(3); int f dx = gamma; units of both terms = kg m^-1 s^-2 = J m^-3

## Result

$$
f(\phi,\nabla\phi) = \frac{3\gamma}{\varepsilon}\,\phi^2(1-\phi)^2 + 3\gamma\varepsilon\,|\nabla\phi|^2
$$

Graded SymPy form:

```
3*gamma/eps*phi(x, y, z, t)**2*(1 - phi(x, y, z, t))**2 + 3*gamma*eps*(Derivative(phi(x, y, z, t), x)**2 + Derivative(phi(x, y, z, t), y)**2 + Derivative(phi(x, y, z, t), z)**2)
```

## Neglected effects

| Effect | Why it is left out |
|---|---|
| Driving force (supersaturation / chemical potential tilt between phases) | Excluded by the coexistence assumption; it would add a term like $\Delta\mu\,p(\phi)$ that breaks the symmetry of the wells. |
| Anisotropy of the ice surface energy | Not asked for; gamma is given as a constant. Ice anisotropy is a few percent. |
| Regolith grains (solid wall, wetting / contact angle) | The problem treats the pore space as ice plus gas only; a third phase or a wall energy would be needed. |
| Vapor density and temperature fields and their coupling | Uniform T and coexistence; these enter the evolution model, not this equilibrium functional. |
| Elastic strain energy in ice, curvature-dependent (Gibbs-Thomson) corrections beyond the functional | Elastic energy is small for unconstrained ice at pore scale; Gibbs-Thomson follows from this functional itself for curved interfaces. |
| Higher-order gradient terms | Corrections of order $(\varepsilon/R)^2$ for interface curvature radius R much larger than $\varepsilon$. |

