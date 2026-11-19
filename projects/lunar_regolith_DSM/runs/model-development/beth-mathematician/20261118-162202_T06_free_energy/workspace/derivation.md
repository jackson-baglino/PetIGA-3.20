# T06_free_energy

Run `20261118-162202_T06_free_energy` · status **done** · written by `scripts/agentkit/render_derivation.py` from `answer.json`; do not edit by hand.

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
| $W$ | height coefficient of the double-well (bulk) term; fixed to 6 gamma/eps in S5 | $\mathrm{J\,m^{-3}}$ | 6e+05 | introduced |
| $\kappa$ | gradient-energy coefficient; fixed to 3 gamma eps in S5 | $\mathrm{J\,m^{-1}}$ | 1.5e-09 | introduced |
| $g(\phi)$ | dimensionless double-well function phi^2 (1-phi)^2 | 1 | 0.0625 | introduced |
| $X$ | coordinate normal to a flat interface, with X=0 at phi=1/2 | $\mathrm{m}$ | 1e-08 | introduced |
| $\ell$ | intrinsic decay length of the equilibrium profile, sqrt(kappa/(2W)) | $\mathrm{m}$ | 5e-09 | introduced |
| $\phi_0(X)$ | equilibrium profile across a flat interface | 1 | 0.5 | introduced |
| $\sigma$ | excess free energy per unit area of a flat equilibrium interface | $\mathrm{J\,m^{-2}}$ | 0.1 | introduced |
| $\delta$ | interface width, defined as 1 / max \|d phi_0/dX\| (the distance over which the tangent at the midpoint rises from 0 to 1) | $\mathrm{m}$ | 2e-08 | introduced |

## Assumptions

- **A1.** A1: The temperature is uniform and constant, at phase coexistence (given), so the bulk free energies of ice and gas are equal and the well has two minima of equal depth; the common bulk value is set to zero.
- **A2.** A2: The interface energy $\gamma$ is isotropic, so the gradient coefficient is a scalar constant.
- **A3.** A3: The free energy is local and depends only on $\phi$ and $\nabla\phi$ (square-gradient, Ginzburg-Landau form); the double well is the quartic $\phi^2(1-\phi)^2$.
- **A4.** A4: For calibration we use a flat, one-dimensional interface at equilibrium, far from walls, with $\phi\to0$ and $\phi\to1$ on the two sides.
- **A5.** A5: Width convention: $\delta=1/\max\|d\phi_0/dX\|=2\varepsilon$, where the profile is $\tfrac12[1+\tanh(X/\varepsilon)]$.

## Derivation

### S1. Ansatz: double well plus gradient

At coexistence the bulk energy must have two equal-depth minima at $\phi=0$ (gas) and $\phi=1$ (ice), so neither phase is favored; $g(\phi)=\phi^2(1-\phi)^2$ is the lowest-order polynomial with this property and a barrier $g(1/2)=1/16$ between them. A gradient term penalizes sharp changes, so the interface gets a finite width and a positive energy; without it the interface would be infinitely sharp.

$$
f = W\,g(\phi) + \frac{\kappa}{2}|\nabla\phi|^2,\qquad g(\phi)=\phi^2(1-\phi)^2
$$

*Uses:* A1, A2, A3. *Operation:* choose the Ginzburg-Landau form with coefficients W and kappa still free. *Checked by:* `check_s1` in `scratch/check.py`: roots of g': [0, 1/2, 1]; g(0)=g(1)=0; g''(0)=g''(1)=2; g(1/2)=1/16

### S2. Euler-Lagrange equation for a flat interface

Equilibrium means $\delta F/\delta\phi=0$. For a flat interface $\phi$ depends only on the normal coordinate $X$, with $\phi\to0$ as $X\to-\infty$ and $\phi\to1$ as $X\to+\infty$.

$$
\kappa\,\phi_0'' - W\,g'(\phi_0) = 0
$$

*Uses:* S1, A4. *Operation:* take the variational derivative of F and restrict to phi = phi_0(X). *Checked by:* `check_s3` in `scratch/check.py`: EL residual: 0 (for the profile of S3)

### S3. First integral and profile

Multiply S2 by $\phi_0'$ and integrate once; the constant is zero because $\phi_0'\to0$ and $g\to0$ in the bulk. This equipartition gives $\phi_0'=\sqrt{2W/\kappa}\,\phi_0(1-\phi_0)$, solved by a logistic profile with decay length $\ell=\sqrt{\kappa/(2W)}$.

$$
\frac{\kappa}{2}\phi_0'^2 = W g(\phi_0)\;\Rightarrow\;\phi_0(X)=\frac{1}{1+e^{-X/\ell}},\quad \ell=\sqrt{\frac{\kappa}{2W}}
$$

*Uses:* S2. *Operation:* multiply by phi_0', integrate over X, then separate variables. *Checked by:* `check_s3` in `scratch/check.py`: first integral residual: 0; EL residual: 0

### S4. Interface energy

The excess energy per area is the integral of $f$ across the interface (the bulk value is 0). With equipartition the integrand is $\kappa\phi_0'^2$, and changing variables from $X$ to $\phi$ gives a closed form.

$$
\sigma=\int_{-\infty}^{\infty}\kappa\,\phi_0'^2\,dX=\int_0^1\sqrt{2W\kappa}\,\phi(1-\phi)\,d\phi=\frac{\sqrt{2W\kappa}}{6}
$$

*Uses:* S3. *Operation:* integrate f over X, substitute dX = dphi/phi_0'. *Checked by:* `check_s4` in `scratch/check.py`: sigma = sqrt(2)*sqrt(W)*sqrt(kappa)/6; numeric sigma 0.577350269189626 equals formula 0.577350269189626 for W=3, kappa=2

### S5. Fix the coefficients

Require $\sigma=\gamma$ and choose $\ell=\varepsilon/2$, so that the profile is $\phi_0=\tfrac12[1+\tanh(X/\varepsilon)]$. Solving the two equations gives $W$ and $\kappa$ in terms of $\gamma$ and $\varepsilon$.

$$
W=\frac{6\gamma}{\varepsilon},\qquad \kappa=3\gamma\varepsilon
$$

*Uses:* S3, S4. *Operation:* solve sigma = gamma, ell = eps/2 for W and kappa. *Checked by:* `check_s5` in `scratch/check.py`: solution: [{W: 6*gamma/epsilon, kappa: 3*epsilon*gamma}]

### S6. Width and direct energy check

With these coefficients the tanh profile solves S2, and its energy integral is exactly $\gamma$. Defining the width as $\delta=1/\max|\phi_0'|$ gives $\delta=2\varepsilon$; the 10%-90% width is $2\ln 3\,\varepsilon\approx2.20\,\varepsilon$.

$$
\phi_0=\tfrac12\left[1+\tanh\frac{X}{\varepsilon}\right],\quad \delta=\frac{1}{\max|\phi_0'|}=2\varepsilon,\quad \sigma=\gamma
$$

*Uses:* S5. *Operation:* substitute S5 into S2 and S4. *Checked by:* `check_s6` in `scratch/check.py`: EL residual tanh profile: 0; 1/max slope = 2*epsilon; 10-90 width = 2*epsilon*log(3) = 2.197 eps; sigma direct = gamma

### S7. Free energy density

Substitute S5 into S1. Both terms have units J m$^{-3}$.

$$
f=\frac{6\gamma}{\varepsilon}\phi^2(1-\phi)^2+\frac{3}{2}\gamma\varepsilon|\nabla\phi|^2
$$

*Uses:* S1, S5. *Operation:* substitution of S5. *Checked by:* `check_units` in `scratch/check.py`: bulk term: kg m^-1 s^-2; gradient term: kg m^-1 s^-2; same: True

## Result

$$
f(\phi,\nabla\phi) = \frac{6\gamma}{\varepsilon}\,\phi^2(1-\phi)^2 + \frac{3}{2}\,\gamma\,\varepsilon\,|\nabla\phi|^2
$$

Graded SymPy form:

```
6*gamma/eps*phi(x, y, z, t)**2*(1 - phi(x, y, z, t))**2 + 3*gamma*eps/2*(Derivative(phi(x, y, z, t), x)**2 + Derivative(phi(x, y, z, t), y)**2 + Derivative(phi(x, y, z, t), z)**2)
```

## Neglected effects

| Effect | Why it is left out |
|---|---|
| Driving force for sublimation and deposition (tilting term, for example proportional to $(T-T_m)$ or vapor supersaturation, times $h(\phi)$) | The work order is at coexistence, so this term is zero by assumption A1. |
| Anisotropy of the ice surface energy (hexagonal ice facets) | Taken as isotropic (A2). The anisotropy of ice-vapor surface energy is a few percent, which is small compared with the geometric effects of pore confinement. |
| Regolith grain surfaces: wetting and the ice-grain and gas-grain energies | The work order treats the pore space as ice plus gas only. Contact angles would enter through boundary terms, not through $f$. |
| Coupling to temperature and vapor density fields (latent heat, vapor diffusion) | At constant uniform temperature and coexistence they do not change the equilibrium interface. They belong to the dynamics, not to this functional. |
| Elastic strain energy in ice and density change on phase change | Ice in open pores against low-pressure gas is nearly stress-free. Elastic energy in the bulk is much smaller than $\gamma/\varepsilon$ for nanometer-to-micron $\varepsilon$. |
| Curvature (Gibbs-Thomson) effects in the calibration | Calibration uses a flat interface. Curved interfaces get corrections of relative order $\varepsilon/R$, which are small when $\varepsilon$ is much smaller than the pore radius $R$. |
| Higher-order gradient terms (for example $(\nabla^2\phi)^2$) | They are excluded by the deliverable. They are of relative order $(\varepsilon/L)^2$ on scales $L\gg\varepsilon$. |

